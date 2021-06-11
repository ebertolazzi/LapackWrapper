#pragma once
#ifndef SPARSETOOL_ITERATIVE_PRECO_SSOR_HH
#define SPARSETOOL_ITERATIVE_PRECO_SSOR_HH

namespace SparseTool {

  /*
  //   #####   #####  ####### ######
  //  #     # #     # #     # #     # 
  //  #       #       #     # #     # 
  //   #####   #####  #     # ######  
  //        #       # #     # #   #   
  //  #     # #     # #     # #    #  
  //   #####   #####  ####### #     # 
  */
  //! Iterative SSOR preconditioner
  template <typename T>
  class SSORpreconditioner : public Preco<SSORpreconditioner<T> > {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    typedef SSORpreconditioner<T> SSORPRECO;
    typedef Preco<SSORPRECO>      PRECO;
    #endif

    typedef T real_type; //!< type of the elements of the preconditioner

  private:

    real_type         omega, omega1;
    integer         maxIter;

    Vector<integer> L_R;
    Vector<integer> L_J;
    Vector<real_type> L_A;

    Vector<integer> U_C;
    Vector<integer> U_I;
    Vector<real_type> U_A;

    Vector<real_type> D;

    Vector<integer> Lnnz, Unnz;

    //! build incomplete LDU decomposition with specified pattern `P` 
    template <typename MAT>
    void
    build_SOR( MAT const & A, real_type const & _omega, integer _maxIter ) {

      SPARSETOOL_ASSERT(
        A.isOrdered(),
        "SSORpreconditioner::build_SOR pattern must be ordered before use"
      )
      SPARSETOOL_ASSERT(
        A.numRows() == A.numCols(),
        "SSORpreconditioner::build_SOR only square matrix allowed"
      )
      SPARSETOOL_ASSERT(
        A.numRows() > 0,
        "SSORpreconditioner::build_SOR empty matrix"
      )

      this -> omega   = _omega;
      this -> omega1  = 1.0/_omega-1.0;
      this -> maxIter = _maxIter;

      // step 0: compute necessary memory
      PRECO::pr_size = A.numRows();
      Lnnz . resize( PRECO::pr_size );
      Unnz . resize( PRECO::pr_size );
      D    . resize( PRECO::pr_size );

      Lnnz = 0;
      Unnz = 0;
      D    = real_type(0);

      for ( A.Begin(); A.End(); A.Next() ) {
        integer i = A.row();
        integer j = A.column();
        if      ( i > j ) ++Lnnz(i);
        else if ( i < j ) ++Unnz(j);
      }

      // step 1: initialize structure
      L_R.resize( PRECO::pr_size + 1 );
      U_C.resize( PRECO::pr_size + 1 );

      L_R(0) = U_C(0) = 0;
      for ( integer i = 0; i < PRECO::pr_size; ++i ) {
        L_R(i+1) = L_R(i) + Lnnz(i);
        U_C(i+1) = U_C(i) + Unnz(i);
      }

      // step 2: allocate memory
      L_A.resize( L_R(PRECO::pr_size) );
      L_J.resize( L_R(PRECO::pr_size) );

      U_A.resize( U_C(PRECO::pr_size) );
      U_I.resize( U_C(PRECO::pr_size) );

      L_A = real_type(0);
      U_A = real_type(0);
      
      // step 3: fill structure
      for ( A.Begin(); A.End(); A.Next() ) {
        integer i = A.row();
        integer j = A.column();
        if      ( i > j ) L_J(L_R(i)+(--Lnnz(i))) = j;
        else if ( i < j ) U_I(U_C(j)+(--Unnz(j))) = i;
      }

      // step 4: sort structure
      for ( integer i = 0; i < PRECO::pr_size; ++i ) {
        std::sort( &L_J(L_R(i)), &L_J(L_R(i+1)) );
        std::sort( &U_I(U_C(i)), &U_I(U_C(i+1)) );
      }

      // insert values
      for ( A.Begin(); A.End(); A.Next() ) {
        integer i = A.row();
        integer j = A.column();
        typename MAT::real_type const val = A.value(); // (i,j);
        if ( i > j ) { // costruisco L
          integer lo  = L_R(i);
          integer hi  = L_R(i+1);
          integer len = hi - lo;
          while ( len > 0 ) {
            integer half = len / 2;
            integer mid  = lo + half;
            if ( L_J(mid) < j ) { lo = mid + 1; len -= half + 1; }
            else                  len = half;
          }
          L_A(lo) = val;
        } else if ( i < j ) { // costruisco U
          integer lo  = U_C(j);
          integer hi  = U_C(j+1);
          integer len = hi - lo;
          while ( len > 0 ) {
            integer half = len / 2;
            integer mid  = lo + half;
            if ( U_I(mid) < i ) { lo = mid + 1; len -= half + 1; }
            else                  len = half;
          }
          U_A(lo) = val;
        } else {
          D(i) = val;
          SPARSETOOL_ASSERT(
            val != real_type(0),
            "SSORpreconditioner::build_SOR D(" << i << ") = " << val
          );

        }
      }
      for ( integer i = 0; i < PRECO::pr_size; ++i )
        SPARSETOOL_ASSERT(
          D(i) != real_type(0),
          "SSORpreconditioner::D(" << i << ") = " << D(i) << " size = " << D.size()
        );
    }

  public:

    SSORpreconditioner(void) : Preco<SSORPRECO>() {}
    
    template <typename MAT>
    SSORpreconditioner( MAT const & M, real_type _omega, integer _maxIter ) : Preco<SSORPRECO>()
    { build_SOR( M, _omega, _maxIter ); }

    //! build the preconditioner from matrix `M` with pattern `P` 
    template <typename MAT>
    void
    build( MAT const & M, real_type _omega, integer _maxIter )
    { build_SOR( M, _omega, _maxIter ); }

    //!
    //! Apply preconditioner to vector `v`
    //! and store result to vector `res`.
    //!
    template <typename VECTOR>
    void
    assPreco( VECTOR & x, VECTOR const & b ) const {
      typedef typename VECTOR::real_type vType;
      integer k;
      x.setZero();
      for ( integer ii = 0; ii < maxIter; ++ii ) {
          
        // calcolo ((1/omega-1)*D-U)*x + b;;
        integer const *   pC  = U_C.data();
        integer const *   pI  = U_I.data();
        real_type const * pUA = U_A.data();

        for ( k=0; k < PRECO::pr_size; ++k ) {
          vType xk = x(k);
          x(k) = b(k) + omega1 * D(k) * xk;
          for ( integer i_cnt = pC[1] - pC[0]; i_cnt > 0; --i_cnt )
            x(*pI++) -= *pUA++ * xk;
          ++pC;
        }

        // solve (D/omega+L) x(n+1) = x(n)
        integer const *   pR  = L_R.data();
        integer const *   pJ  = L_J.data();
        real_type const * pLA = L_A.data();

        for ( k=0; k < PRECO::pr_size; ++k ) {
          vType tmp(0);
          for ( integer i_cnt = pR[1] - pR[0]; i_cnt > 0; --i_cnt )
            tmp += *pLA++ * x(*pJ++);
          x(k) = omega*(x(k)-tmp)/D(k);
          ++pR;
        };
      }

      for ( integer ii = 0; ii < maxIter; ++ii ) {

        // calcolo ((1/omega-1)*D-L)*x + b;
        integer const *   pR  = L_R.data() + PRECO::pr_size;
        integer const *   pJ  = L_J.data() + *pR;
        real_type const * pLA = L_A.data() + *pR;
        k = PRECO::pr_size;
        do {
          --k; --pR;
          vType tmp(0);
          for ( integer i_cnt = pR[1] - pR[0]; i_cnt > 0; --i_cnt )
            tmp += *--pLA * x(*--pJ);
          x(k) = omega1*D(k)*x(k)+b(k)-tmp;
        } while ( k > 0 );

        // solve (D/omega+U) x(n+1) = x(n)
        integer const *   pC  = U_C.data() + PRECO::pr_size;
        integer const *   pI  = U_I.data() + *pC;
        real_type const * pUA = U_A.data() + *pC;
        k = PRECO::pr_size;
        do {
          --k; --pC;
          vType xk = x(k)*omega/D(k);
          for ( integer i_cnt = pC[1] - pC[0]; i_cnt > 0; --i_cnt )
            x(*--pI) -= *--pUA * xk;
          x(k) = xk;
        } while ( k > 0 );
      }
    }

  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,SSORpreconditioner<TP> >
  operator / (Vector<T> const & v, SSORpreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,SSORpreconditioner<TP> >(v,P); }
  #endif

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace SparseToolLoad {
  using ::SparseTool::SSORpreconditioner;
}
#endif

#endif
