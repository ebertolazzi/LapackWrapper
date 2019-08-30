#ifndef SPARSETOOL_ITERATIVE_PRECO_HSS_CHEBYSHEV_HH
#define SPARSETOOL_ITERATIVE_PRECO_HSS_CHEBYSHEV_HH

using namespace std;

namespace SparseTool {

  /*
  //    ____ _   _ _____ ______   ______  _   _ _______     __
  //   / ___| | | | ____| __ ) \ / / ___|| | | | ____\ \   / /
  //  | |   | |_| |  _| |  _ \\ V /\___ \| |_| |  _|  \ \ / /
  //  | |___|  _  | |___| |_) || |  ___) |  _  | |___  \ V /
  //   \____|_| |_|_____|____/ |_| |____/|_| |_|_____|  \_/
  //
  */
  //! Incomplete \c LDU preconditioner
  template <typename T>
  class HSS_CHEBYSHEV_Preconditioner : public Preco<HSS_CHEBYSHEV_Preconditioner<T> > {
  public:

    //! \cond NODOC
    typedef HSS_CHEBYSHEV_Preconditioner<T> HSS_CHEBYSHEV_PRECO;
    typedef Preco<HSS_CHEBYSHEV_PRECO>      PRECO;

    //! \endcond
    typedef T valueType; //!< type of the elements of the preconditioner
    typedef typename T::value_type rvalueType; //!< type of the elements of the preconditioner

  private:

    indexType          neq, mdegree;
    rvalueType         sqrt_epsilon;

    Vector<indexType>  A_R;
    Vector<indexType>  A_J;
    Vector<rvalueType> A_A;
    Vector<indexType>  Annz;

    mutable Vector<valueType> s0, s1, y;

    //! build incomplete LDU decomposition with specified pattern \c P
    template <typename MAT>
    void
    build_HSS_CHEBYSHEV( MAT const & A, indexType m, rvalueType delta ) {
      

      SPARSETOOL_ASSERT(
        A.isOrdered(),
        "HSS_CHEBYSHEV_Preconditioner::build_LDU pattern must be ordered before use"
      )
      SPARSETOOL_ASSERT(
        A.numRows() == A.numCols(),
        "HSS_CHEBYSHEV_Preconditioner::build_LDU only square matrix allowed"
      )
      SPARSETOOL_ASSERT(
        A.numRows() > 0,
        "HSS_CHEBYSHEV_Preconditioner::build_LDU empty matrix"
      )

      mdegree = m;
      sqrt_epsilon = std::pow( (1+std::sqrt(1-delta*delta))/delta,1.0/(m+1.0) );
      sqrt_epsilon = (sqrt_epsilon-1)/(sqrt_epsilon+1);

      SPARSETOOL_ASSERT( sqrt_epsilon > 0 && sqrt_epsilon < 1,
                         "HSS_CHEBYSHEV_Preconditioner::build_LDU computed epsilon must be in the interval (0,1)" )

      neq = A.numRows();
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
      rvalueType e = sqrt_epsilon*sqrt_epsilon;
      rvalueType a = 2.0/(1.0-e);
      rvalueType b = -(1.0+e)/(1.0-e);
      rvalueType c = (sqrt_epsilon-1)/(sqrt_epsilon+1);
      // s0 = 1.5*v; s1 = 4*v - 10/3 * A*v
      indexType  const * pR = & A_R.front();
      indexType  const * pJ = & A_J.front();
      rvalueType const * pA = & A_A.front();
      for ( indexType k=0; k < PRECO::pr_size; ++k, ++pR ) {
        valueType Av(0,0);
        for ( indexType i_cnt = pR[1] - pR[0]; i_cnt > 0; --i_cnt )
          Av += *pA++ * v(*pJ++);
        s1(k) = (2.0/(1.0+e)) * v(k);
        _y(k) = (8.0*(1.0+e)/(e*(e+6.0)+1.0)) * v(k) - (8.0/(e*(e+6.0)+1.0)) * Av;
      };
      rvalueType cm1  = 1;
      rvalueType cn   = c;
      rvalueType cp1  = c*c;
      for ( indexType n = 2; n <= mdegree; ++n ) {
        cm1  = cn;
        cn   = cp1;
        cp1 *= c;
        s0   = s1;
        s1   = _y;
        rvalueType g0 = (cn+1/cn)/(cp1+1/cp1);
        rvalueType aa = 2*a*g0;
        rvalueType bb = 2*b*g0;
        rvalueType cc = -(cm1+1/cm1)/(cp1+1/cp1);
        pR = & A_R.front();
        pJ = & A_J.front();
        pA = & A_A.front();
        for ( indexType k=0; k < PRECO::pr_size; ++k, ++pR ) {
          valueType As1(0,0);
          for ( indexType i_cnt = pR[1] - pR[0]; i_cnt > 0; --i_cnt )
            As1  += *pA++ * s1(*pJ++);
            _y(k) = aa*(As1-v(k))+bb*s1(k)+cc*s0(k);
        };
      }
    }

  public:

    HSS_CHEBYSHEV_Preconditioner(void) : Preco<HSS_CHEBYSHEV_PRECO>() {}
    
    template <typename MAT>
    HSS_CHEBYSHEV_Preconditioner(
      MAT const & M,
      indexType   m,
      rvalueType  delta
    ) : Preco<HSS_CHEBYSHEV_PRECO>()
    { build_HSS_CHEBYSHEV( M, m, delta ); }

    //! build the preconditioner from matrix \c M
    template <typename MAT>
    void
    build( MAT const & M, indexType m, rvalueType delta )
    { build_HSS_CHEBYSHEV(M,m,delta); }

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
  Vector_V_div_P<Vector<T>,HSS_CHEBYSHEV_Preconditioner<TP> >
  operator / (Vector<T> const & v, HSS_CHEBYSHEV_Preconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,HSS_CHEBYSHEV_Preconditioner<TP> >(v,P); }
  //! \endcond

}

namespace SparseToolLoad {
  using ::SparseTool::HSS_CHEBYSHEV_Preconditioner;
}

#endif
