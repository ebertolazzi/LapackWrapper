#pragma once
#ifndef SPARSETOOL_ITERATIVE_PRECO_JACOBI_HH
#define SPARSETOOL_ITERATIVE_PRECO_JACOBI_HH

namespace SparseTool {

  /*
  //        #                               
  //        #   ##    ####   ####  #####  # 
  //        #  #  #  #    # #    # #    # # 
  //        # #    # #      #    # #####  # 
  //  #     # ###### #      #    # #    # # 
  //  #     # #    # #    # #    # #    # # 
  //   #####  #    #  ####   ####  #####  # 
  */
  //! Iterative JACOBI preconditioner
  template <typename T>
  class JACOBIpreconditioner : public Preco<JACOBIpreconditioner<T> > {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    typedef JACOBIpreconditioner<T> JACOBIPRECO;
    typedef Preco<JACOBIPRECO>      PRECO;
    #endif

    typedef T real_type; //!< type of the elements of the preconditioner

  private:

    real_type         omega, omega1;
    integer         maxIter;

    Vector<integer> LU_R;
    Vector<integer> LU_J;
    Vector<real_type> LU_A;

    Vector<real_type> D;

    mutable Vector<real_type> TMP;

    Vector<integer> LUnnz;

    //! build incomplete LDU decomposition with specified pattern `P` 
    template <typename MAT>
    void
    build_JACOBI( MAT const & A, real_type const & _omega, integer _maxIter ) {

      SPARSETOOL_ASSERT(
        A.isOrdered(),
        "JACOBIpreconditioner::build_JACOBI pattern must be ordered before use"
      )
      SPARSETOOL_ASSERT(
        A.numRows() == A.numCols(),
        "JACOBIpreconditioner::build_JACOBI only square matrix allowed"
      )
      SPARSETOOL_ASSERT(
        A.numRows() > 0,
        "JACOBIpreconditioner::build_JACOBI empty matrix"
      )

      this -> omega   = _omega;
      this -> omega1  = 1.0/_omega-1.0;
      this -> maxIter = _maxIter;

      // step 0: compute necessary memory
      PRECO::pr_size = A.numRows();
      LUnnz . resize( PRECO::pr_size );
      D     . resize( PRECO::pr_size );
      TMP   . resize( PRECO::pr_size );

      LUnnz . setZero();
      D     . setZero();

      for ( A.Begin(); A.End(); A.Next() ) {
        integer i = A.row();
        integer j = A.column();
        if ( i != j ) ++LUnnz(i);
      }

      // step 1: initialize structure
      LU_R.resize( PRECO::pr_size + 1 );
      LU_R(0) = 0;
      for ( integer i = 0; i < PRECO::pr_size; ++i ) LU_R(i+1) = LU_R(i) + LUnnz(i);

      // step 2: allocate memory
      LU_A.resize( LU_R(PRECO::pr_size) );
      LU_J.resize( LU_R(PRECO::pr_size) );
      LU_A.setZero();

      // step 3: fill structure
      for ( A.Begin(); A.End(); A.Next() ) {
        integer i = A.row();
        integer j = A.column();
        if ( i != j ) LU_J(LU_R(i)+(--LUnnz(i))) = j;
      }

      // step 4: sort structure
      for ( integer i = 0; i < PRECO::pr_size; ++i )
        std::sort( &LU_J(LU_R(i)), &LU_J(LU_R(i+1)) );

      // insert values
      for ( A.Begin(); A.End(); A.Next() ) {
        integer i   = A.row();
        integer j   = A.column();
        real_type val = A.value(); // (i,j);
        if ( i != j ) { // costruisco LU
          integer lo  = LU_R(i);
          integer hi  = LU_R(i+1);
          integer len = hi - lo;
          while ( len > 0 ) {
            integer half = len / 2;
            integer mid  = lo + half;
            if ( LU_J(mid) < j ) { lo = mid + 1; len -= half + 1; }
            else                   len = half;
          }
          LU_A(lo) = val;
        } else {
          D(i) = val;
        }
      }
      for ( integer i = 0; i < PRECO::pr_size; ++i )
        SPARSETOOL_ASSERT(
          D(i) != real_type(0),
          "JACOBIpreconditioner::D(" << i <<
          ") = " << D(i) << " size = " << D.size()
        );
    }

  public:

    JACOBIpreconditioner(void) : Preco<JACOBIPRECO>() {}
    
    template <typename MAT>
    JACOBIpreconditioner( MAT const & M, real_type _omega, integer _maxIter ) : Preco<JACOBIPRECO>()
    { build_JACOBI( M, _omega, _maxIter ); }

    //! build the preconditioner from matrix `M` with pattern `P` 
    template <typename MAT>
    void
    build( MAT const & M, real_type _omega, integer _maxIter )
    { build_JACOBI( M, _omega, _maxIter ); }

    //!
    //! Apply preconditioner to vector `v`
    //! and store result to vector `res`.
    //!
    template <typename VECTOR>
    void
    assPreco( VECTOR & x, VECTOR const & b ) const {
      x = omega*(b/D);
      for ( integer ii = 0; ii < maxIter; ++ii ) {

        // calcolo (D/omega)^{-1} [ ((1/omega-1)*D-L-U)*x + b ];
        integer const *   _pR = LU_R.data();
        integer const *   _pJ = LU_J.data();
        real_type const * _pA = LU_A.data();

        for ( integer k=0; k < PRECO::pr_size; ++k ) {
          real_type tmp(0);
          for ( integer i_cnt = _pR[1] - _pR[0]; i_cnt > 0; --i_cnt )
            tmp += *_pA++ * x(*_pJ++);
          TMP(k) = b(k) + omega1 * D(k) * x(k) - tmp;
          ++_pR;
        }
        x = omega * (TMP / D);
      }
    }

  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,JACOBIpreconditioner<TP> >
  operator / (Vector<T> const & v, JACOBIpreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,JACOBIpreconditioner<TP> >(v,P); }
  #endif

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace SparseToolLoad {
  using ::SparseTool::JACOBIpreconditioner;
}
#endif

#endif
