#ifndef SPARSETOOL_ITERATIVE_PRECO_HSS_OPOLY_HH
#define SPARSETOOL_ITERATIVE_PRECO_HSS_OPOLY_HH

using namespace std;

namespace SparseTool {

  /*
  //                     _       
  //    ___  _ __   ___ | |_   _ 
  //   / _ \| '_ \ / _ \| | | | |
  //  | (_) | |_) | (_) | | |_| |
  //   \___/| .__/ \___/|_|\__, |
  //        |_|            |___/ 
  */
  //! Incomplete \c LDU preconditioner
  template <typename T>
  class HSS_OPOLY_Preconditioner : public Preco<HSS_OPOLY_Preconditioner<T> > {
  public:

    //! \cond NODOC
    typedef HSS_OPOLY_Preconditioner<T> HSS_OPOLY_PRECO;
    typedef Preco<HSS_OPOLY_PRECO>      PRECO;

    //! \endcond
    typedef T valueType; //!< type of the elements of the preconditioner
    typedef typename T::value_type rvalueType; //!< type of the elements of the preconditioner

  private:

    indexType neq, mdegree;

    Vector<indexType>  A_R;
    Vector<indexType>  A_J;
    Vector<rvalueType> A_A;
    Vector<indexType>  Annz;

    mutable Vector<valueType> s0, s1, y;

    //! build incomplete LDU decomposition with specified pattern \c P
    template <typename MAT>
    void
    build_HSS_OPOLY( MAT const & A, indexType m ) {

      SPARSETOOL_ASSERT(
        A.isOrdered(),
        "HSS_OPOLY_Preconditioner::build_LDU pattern must be ordered before use"
      )
      SPARSETOOL_ASSERT(
        A.numRows() == A.numCols(),
        "HSS_OPOLY_Preconditioner::build_LDU only square matrix allowed"
      )
      SPARSETOOL_ASSERT(
        A.numRows() > 0,
        "HSS_OPOLY_Preconditioner::build_LDU empty matrix"
      )

      mdegree = m;
      neq     = A.numRows();
      s0 . resize(neq);
      s1 . resize(neq);
      y  . resize(neq);

      // step 0: compute necessary memory
      PRECO::pr_size = A.numRows();
      Annz.resize( PRECO::pr_size );
      Annz.setZero();
        
      for ( A.Begin(); A.End(); A.Next() ) {
        indexType i = A.row();
        //indexType j = A.column();
        ++Annz(i);
      }
        
      // step 1: initialize structure
      A_R.resize( PRECO::pr_size + 1 );
      A_R(0) = 0;
      for ( indexType i = 0; i < PRECO::pr_size; ++i ) A_R(i+1) = A_R(i) + Annz(i);
        
      // step 2: allocate memory
      A_A.resize( A_R(PRECO::pr_size) );
      A_J.resize( A_R(PRECO::pr_size) );
      A_A.setZero();
        
      // step 3: fill structure
      for ( A.Begin(); A.End(); A.Next() ) {
        indexType i = A.row();
        indexType j = A.column();
        A_J(A_R(i)+(--Annz(i))) = j;
      }

      // step 4: sort structure
      for ( indexType i = 0; i < PRECO::pr_size; ++i ) sort( &A_J(A_R(i)), &A_J(A_R(i+1)) );
        
      // insert values
      for ( A.Begin(); A.End(); A.Next() ) {
        indexType i   = A.row();
        indexType j   = A.column();
        indexType lo  = A_R(i);
        indexType hi  = A_R(i+1);
        indexType len = hi - lo;
        while ( len > 0 ) {
          indexType half = len / 2;
          indexType mid  = lo + half;
          if ( A_J(mid) < j ) { lo = mid + 1; len -= half + 1; }
          else                  len = half;
        }
        A_A(lo) = A.value().real() + A.value().imag();
      }
    }

    //! apply preconditioner to vector \c v and store result in vector \c y
    void
    mulPoly( Vector<valueType> & _y, Vector<valueType> const & v ) const {
      // s0 = 1.5*v; s1 = 4*v - 10/3 * A*v
      indexType  const * pR = & A_R.front();
      indexType  const * pJ = & A_J.front();
      rvalueType const * pA = & A_A.front();
      for ( indexType k=0; k < PRECO::pr_size; ++k, ++pR ) {
        valueType Av(0,0);
        for ( indexType i_cnt = pR[1] - pR[0]; i_cnt > 0; --i_cnt )
          Av += *pA++ * v(*pJ++);
        s1(k) = 1.5 * v(k);
        _y(k) = 4.0 * v(k) - (10./3.) * Av;
      };
      for ( indexType n = 2; n <= mdegree; ++n ) {
        s0 = s1;
        s1 = _y;
        rvalueType delta = ((6*n+12)*n+4.0)/((2*n+1)*(n+2)*(n+2));
        rvalueType a = -4+(6*n+10.0)/((n+2)*(n+2));
        rvalueType b = 2-delta;
        rvalueType c = -1+delta;

        pR = & A_R.front();
        pJ = & A_J.front();
        pA = & A_A.front();
        for ( indexType k=0; k < PRECO::pr_size; ++k, ++pR ) {
          valueType As1(0,0);
          for ( indexType i_cnt = pR[1] - pR[0]; i_cnt > 0; --i_cnt )
            As1  += *pA++ * s1(*pJ++);
            _y(k) = a*(As1-v(k))+b*s1(k)+c*s0(k);
        };
      }
    }

  public:

    HSS_OPOLY_Preconditioner(void) : Preco<HSS_OPOLY_PRECO>() {}
    
    template <typename MAT>
    HSS_OPOLY_Preconditioner( MAT const & M, indexType m ) : Preco<HSS_OPOLY_PRECO>()
    { build_HSS_OPOLY( M, m ); }

    //! build the preconditioner from matrix \c M
    template <typename MAT>
    void
    build( MAT const & M, indexType m )
    { build_HSS_OPOLY(M,m); }

    //! apply preconditioner to vector \c v and store result to vector \c y
    template <typename VECTOR>
    void
    assPreco( VECTOR & _y, VECTOR const & v ) const {
      mulPoly( _y, v );
      _y *= valueType(0.5,-0.5);
    }

  };

  //! \cond NODOC
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,HSS_OPOLY_Preconditioner<TP> >
  operator / (Vector<T> const & v, HSS_OPOLY_Preconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,HSS_OPOLY_Preconditioner<TP> >(v,P); }
  //! \endcond

}

namespace SparseToolLoad {
  using ::SparseTool::HSS_OPOLY_Preconditioner;
}

#endif
