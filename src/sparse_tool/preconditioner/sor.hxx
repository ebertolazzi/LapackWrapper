#ifndef SPARSETOOL_ITERATIVE_PRECO_SOR_HH
#define SPARSETOOL_ITERATIVE_PRECO_SOR_HH

using namespace std;

namespace SparseTool {

  /*
  //   #####  ####### ######  
  //  #     # #     # #     # 
  //  #       #     # #     # 
  //   #####  #     # ######  
  //        # #     # #   #   
  //  #     # #     # #    #  
  //   #####  ####### #     # 
  */
  //! Iterative SOR preconditioner
  template <typename T>
  class SORpreconditioner : public Preco<SORpreconditioner<T> > {
  public:

    //! \cond NODOC
    typedef SORpreconditioner<T> SORPRECO;
    typedef Preco<SORPRECO>      PRECO;

    //! \endcond
    typedef T valueType; //!< type of the elements of the preconditioner

  private:

    valueType         omega, omega1;
    indexType         maxIter;

    Vector<indexType> L_R;
    Vector<indexType> L_J;
    Vector<valueType> L_A;

    Vector<indexType> U_C;
    Vector<indexType> U_I;
    Vector<valueType> U_A;

    Vector<valueType> D;

    Vector<indexType> Lnnz, Unnz;

    //! build incomplete LDU decomposition with specified pattern \c P
    template <typename MAT>
    void
    build_SOR( MAT const & A, valueType const & _omega, indexType _maxIter ) {

      SPARSETOOL_ASSERT(
        A.isOrdered(),
        "SORpreconditioner::build_SOR pattern must be ordered before use"
      )
      SPARSETOOL_ASSERT(
        A.numRows() == A.numCols(),
        "SORpreconditioner::build_SOR only square matrix allowed"
      )
      SPARSETOOL_ASSERT(
        A.numRows() > 0,
        "SORpreconditioner::build_SOR empty matrix"
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
      D    = valueType(0);

      for ( A.Begin(); A.End(); A.Next() ) {
        indexType i = A.row();
        indexType j = A.column();
        if      ( i > j ) ++Lnnz(i);
        else if ( i < j ) ++Unnz(j);
      }

      // step 1: initialize structure
      L_R.resize( PRECO::pr_size + 1 );
      U_C.resize( PRECO::pr_size + 1 );

      L_R(0) = U_C(0) = 0;
      for ( indexType i = 0; i < PRECO::pr_size; ++i ) {
        L_R(i+1) = L_R(i) + Lnnz(i);
        U_C(i+1) = U_C(i) + Unnz(i);
      }

      // step 2: allocate memory
      L_A.resize( L_R(PRECO::pr_size) );
      L_J.resize( L_R(PRECO::pr_size) );

      U_A.resize( U_C(PRECO::pr_size) );
      U_I.resize( U_C(PRECO::pr_size) );

      L_A = valueType(0);
      U_A = valueType(0);
      
      // step 3: fill structure
      for ( A.Begin(); A.End(); A.Next() ) {
        indexType i = A.row();
        indexType j = A.column();
        if      ( i > j ) L_J(L_R(i)+(--Lnnz(i))) = j;
        else if ( i < j ) U_I(U_C(j)+(--Unnz(j))) = i;
      }

      // step 4: sort structure
      for ( indexType i = 0; i < PRECO::pr_size; ++i ) {
        sort( &L_J(L_R(i)), &L_J(L_R(i+1)) );
        sort( &U_I(U_C(i)), &U_I(U_C(i+1)) );
      }

      // insert values
      for ( A.Begin(); A.End(); A.Next() ) {
        indexType i = A.row();
        indexType j = A.column();
        typename MAT::valueType const val = A.value(); // (i,j);
        if ( i > j ) { // costruisco L
          indexType lo  = L_R(i);
          indexType hi  = L_R(i+1);
          indexType len = hi - lo;
          while ( len > 0 ) {
            indexType half = len / 2;
            indexType mid  = lo + half;
            if ( L_J(mid) < j ) { lo = mid + 1; len -= half + 1; }
            else                  len = half;
          }
          L_A(lo) = val;
        } else if ( i < j ) { // costruisco U
          indexType lo  = U_C(j);
          indexType hi  = U_C(j+1);
          indexType len = hi - lo;
          while ( len > 0 ) {
            indexType half = len / 2;
            indexType mid  = lo + half;
            if ( U_I(mid) < i ) { lo = mid + 1; len -= half + 1; }
            else                  len = half;
          }
          U_A(lo) = val;
        } else {
          D(i) = val;
          SPARSETOOL_ASSERT(
            val != valueType(0),
            "SORpreconditioner::build_SOR D(" << i << ") = " << val
          );

        }
      }
      for ( indexType i = 0; i < PRECO::pr_size; ++i )
        SPARSETOOL_ASSERT(
          D(i) != valueType(0),
          "SORpreconditioner::D(" << i << ") = " << D(i) << " size = " << D.size()
        );
    }

  public:

    SORpreconditioner(void) : Preco<SORPRECO>() {}
    
    template <typename MAT>
    SORpreconditioner( MAT const & M, valueType _omega, indexType _maxIter ) : Preco<SORPRECO>() 
    { build_SOR( M, _omega, _maxIter ); }

    //! build the preconditioner from matrix \c M with pattern \c P
    template <typename MAT>
    void
    build( MAT const & M, valueType _omega, indexType _maxIter )
    { build_SOR( M, _omega, _maxIter ); }

    //! apply preconditioner to vector \c v and store result to vector \c res
    template <typename VECTOR>
    void
    assPreco( VECTOR & x, VECTOR const & b ) const {
      typedef typename VECTOR::valueType vType;
      x = vType(0);
      for ( indexType ii = 0; ii < maxIter; ++ii ) {
        // calcolo ((1/omega-1)*D-U)*x + b;
        indexType const * pC  = & U_C.front();
        indexType const * pI  = & U_I.front();
        valueType const * pUA = & U_A.front();

        for ( indexType k=0; k < PRECO::pr_size; ++k ) {
          vType xk = x(k);
          x(k) = b(k) + omega1 * D(k) * xk;
          for ( indexType i_cnt = pC[1] - pC[0]; i_cnt > 0; --i_cnt )
            x(*pI++) -= *pUA++ * xk;
          ++pC;
        }

        // solve (D/omega+L) x(n+1) = x(n)
        indexType const * pR  = & L_R.front();
        indexType const * pJ  = & L_J.front();
        valueType const * pLA = & L_A.front();

        for ( indexType k=0; k < PRECO::pr_size; ++k ) {
          vType tmp(0);
          for ( indexType i_cnt = pR[1] - pR[0]; i_cnt > 0; --i_cnt )
            tmp += *pLA++ * x(*pJ++);
          x(k) = omega*(x(k)-tmp)/D(k);
          ++pR;
        };
      }
    }

  };

  //! \cond NODOC
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,SORpreconditioner<TP> >
  operator / (Vector<T> const & v, SORpreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,SORpreconditioner<TP> >(v,P); }
  //! \endcond

}

namespace SparseToolLoad {
  using ::SparseTool::SORpreconditioner;
}

#endif
