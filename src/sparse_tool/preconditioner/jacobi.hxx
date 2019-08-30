#ifndef SPARSETOOL_ITERATIVE_PRECO_JACOBI_HH
#define SPARSETOOL_ITERATIVE_PRECO_JACOBI_HH

using namespace std;

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

    //! \cond NODOC
    typedef JACOBIpreconditioner<T> JACOBIPRECO;
    typedef Preco<JACOBIPRECO>      PRECO;

    //! \endcond
    typedef T valueType; //!< type of the elements of the preconditioner

  private:

    valueType         omega, omega1;
    indexType         maxIter;

    Vector<indexType> LU_R;
    Vector<indexType> LU_J;
    Vector<valueType> LU_A;

    Vector<valueType> D;

    mutable Vector<valueType> TMP;

    Vector<indexType> LUnnz;

    //! build incomplete LDU decomposition with specified pattern \c P
    template <typename MAT>
    void
    build_JACOBI( MAT const & A, valueType const & _omega, indexType _maxIter ) {

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
        indexType i = A.row();
        indexType j = A.column();
        if ( i != j ) ++LUnnz(i);
      }

      // step 1: initialize structure
      LU_R.resize( PRECO::pr_size + 1 );
      LU_R(0) = 0;
      for ( indexType i = 0; i < PRECO::pr_size; ++i ) LU_R(i+1) = LU_R(i) + LUnnz(i);

      // step 2: allocate memory
      LU_A.resize( LU_R(PRECO::pr_size) );
      LU_J.resize( LU_R(PRECO::pr_size) );
      LU_A.setZero();

      // step 3: fill structure
      for ( A.Begin(); A.End(); A.Next() ) {
        indexType i = A.row();
        indexType j = A.column();
        if ( i != j ) LU_J(LU_R(i)+(--LUnnz(i))) = j;
      }

      // step 4: sort structure
      for ( indexType i = 0; i < PRECO::pr_size; ++i )
        sort( &LU_J(LU_R(i)), &LU_J(LU_R(i+1)) );

      // insert values
      for ( A.Begin(); A.End(); A.Next() ) {
        indexType i   = A.row();
        indexType j   = A.column();
        valueType val = A.value(); // (i,j);
        if ( i != j ) { // costruisco LU
          indexType lo  = LU_R(i);
          indexType hi  = LU_R(i+1);
          indexType len = hi - lo;
          while ( len > 0 ) {
            indexType half = len / 2;
            indexType mid  = lo + half;
            if ( LU_J(mid) < j ) { lo = mid + 1; len -= half + 1; }
            else                   len = half;
          }
          LU_A(lo) = val;
        } else {
          D(i) = val;
        }
      }
      for ( indexType i = 0; i < PRECO::pr_size; ++i )
        SPARSETOOL_ASSERT(
          D(i) != valueType(0),
          "JACOBIpreconditioner::D(" << i <<
          ") = " << D(i) << " size = " << D.size()
        );
    }

  public:

    JACOBIpreconditioner(void) : Preco<JACOBIPRECO>() {}
    
    template <typename MAT>
    JACOBIpreconditioner( MAT const & M, valueType _omega, indexType _maxIter ) : Preco<JACOBIPRECO>()
    { build_JACOBI( M, _omega, _maxIter ); }

    //! build the preconditioner from matrix \c M with pattern \c P
    template <typename MAT>
    void
    build( MAT const & M, valueType _omega, indexType _maxIter )
    { build_JACOBI( M, _omega, _maxIter ); }

    //! apply preconditioner to vector \c v and store result to vector \c res
    template <typename VECTOR>
    void
    assPreco( VECTOR & x, VECTOR const & b ) const {
      x = omega*(b/D);
      for ( indexType ii = 0; ii < maxIter; ++ii ) {

        // calcolo (D/omega)^{-1} [ ((1/omega-1)*D-L-U)*x + b ];
        indexType const * _pR = & LU_R.front();
        indexType const * _pJ = & LU_J.front();
        valueType const * _pA = & LU_A.front();

        for ( indexType k=0; k < PRECO::pr_size; ++k ) {
          valueType tmp(0);
          for ( indexType i_cnt = _pR[1] - _pR[0]; i_cnt > 0; --i_cnt )
            tmp += *_pA++ * x(*_pJ++);
          TMP(k) = b(k) + omega1 * D(k) * x(k) - tmp;
          ++_pR;
        }
        x = omega * (TMP / D);
      }
    }

  };

  //! \cond NODOC
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,JACOBIpreconditioner<TP> >
  operator / (Vector<T> const & v, JACOBIpreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,JACOBIpreconditioner<TP> >(v,P); }
  //! \endcond

}

namespace SparseToolLoad {
  using ::SparseTool::JACOBIpreconditioner;
}

#endif
