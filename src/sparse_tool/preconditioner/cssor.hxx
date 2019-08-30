#ifndef SPARSETOOL_ITERATIVE_PRECO_CSSOR_HH
#define SPARSETOOL_ITERATIVE_PRECO_CSSOR_HH

using namespace std;

namespace SparseTool {

  /*
  //   #####   #####   #####  ####### ######  
  //  #     # #     # #     # #     # #     # 
  //  #       #       #       #     # #     # 
  //  #        #####   #####  #     # ######  
  //  #             #       # #     # #   #   
  //  #     # #     # #     # #     # #    #  
  //   #####   #####   #####  ####### #     # 
  */
  //! Iterative SSOR preconditioner
  /*
  //  Precondizionatore SSOR per sistema del tipo
  //
  //  /      \ / \   / \
  //  | A -B | |x| = |b|
  //  | B  A | |y|   |c|
  //  \      / \ /   \ /
  */
  template <typename T>
  class CSSORpreconditioner : public Preco<CSSORpreconditioner<T> > {
  public:

    //! \cond NODOC
    typedef CSSORpreconditioner<T> CSSORPRECO;
    typedef Preco<CSSORPRECO>      PRECO;

    //! \endcond
    typedef typename T::value_type valueType; //!< type of the elements of the preconditioner

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

    Vector<indexType> B_R;
    Vector<indexType> B_J;
    Vector<valueType> B_A;

    Vector<indexType> Lnnz, Unnz, Bnnz;

    mutable Vector<valueType> br, bi, x, y;

    //! build incomplete LDU decomposition with specified pattern \c P
    template <typename MAT>
    void
    build_SOR( MAT const & A, valueType const & _omega, indexType _maxIter ) {

      SPARSETOOL_ASSERT(
        A.isOrdered(),
        "CSSORpreconditioner::build_SOR pattern must be ordered before use"
      )
      SPARSETOOL_ASSERT(
        A.numRows() == A.numCols(),
        "CSSORpreconditioner::build_SOR only square matrix allowed"
      )
      SPARSETOOL_ASSERT(
        A.numRows() > 0,
        "CSSORpreconditioner::build_SOR empty matrix"
      )

      this -> omega   = _omega;
      this -> omega1  = 1.0/_omega-1.0;
      this -> maxIter = _maxIter;

      // step 0: compute necessary memory
      PRECO::pr_size = A.numRows();
      Lnnz . resize( PRECO::pr_size );
      Unnz . resize( PRECO::pr_size );
      Bnnz . resize( PRECO::pr_size );
      D    . resize( PRECO::pr_size );

      Lnnz = 0;
      Unnz = 0;
      Bnnz = 0;
      D    = valueType(0);

      for ( A.Begin(); A.End(); A.Next() ) {
        indexType i = A.row();
        indexType j = A.column();
        ++Bnnz(i);
        if      ( i > j ) ++Lnnz(i);
        else if ( i < j ) ++Unnz(j);
      }

      // step 1: initialize structure
      L_R.resize( PRECO::pr_size + 1 );
      U_C.resize( PRECO::pr_size + 1 );
      B_R.resize( PRECO::pr_size + 1 );

      L_R(0) = U_C(0) = B_R(0) = 0;
      for ( indexType i = 0; i < PRECO::pr_size; ++i ) {
        L_R(i+1) = L_R(i) + Lnnz(i);
        U_C(i+1) = U_C(i) + Unnz(i);
        B_R(i+1) = B_R(i) + Bnnz(i);
      }

      // step 2: allocate memory
      L_A.resize( L_R(PRECO::pr_size) );
      L_J.resize( L_R(PRECO::pr_size) );

      U_A.resize( U_C(PRECO::pr_size) );
      U_I.resize( U_C(PRECO::pr_size) );

      B_A.resize( B_R(PRECO::pr_size) );
      B_J.resize( B_R(PRECO::pr_size) );

      L_A.setZero();
      U_A.setZero();
      B_A.setZero();
      
      // step 3: fill structure
      for ( A.Begin(); A.End(); A.Next() ) {
        indexType i = A.row();
        indexType j = A.column();
        B_J(B_R(i)+(--Bnnz(i))) = j;
        if      ( i > j ) L_J(L_R(i)+(--Lnnz(i))) = j;
        else if ( i < j ) U_I(U_C(j)+(--Unnz(j))) = i;
      }

      // step 4: sort structure
      for ( indexType i = 0; i < PRECO::pr_size; ++i ) {
        sort( &L_J(L_R(i)), &L_J(L_R(i+1)) );
        sort( &U_I(U_C(i)), &U_I(U_C(i+1)) );
        sort( &B_J(B_R(i)), &B_J(B_R(i+1)) );
      }

      // step 5: insert values
      for ( A.Begin(); A.End(); A.Next() ) {
        indexType i    = A.row();
        indexType j    = A.column();
        valueType rval = A.value().real();
        valueType ival = A.value().imag();
        { // costruisco B
          indexType lo  = B_R(i);
          indexType hi  = B_R(i+1);
          indexType len = hi - lo;
          while ( len > 0 ) {
            indexType half = len / 2;
            indexType mid  = lo + half;
            if ( B_J(mid) < j ) { lo = mid + 1; len -= half + 1; }
            else                  len = half;
          }
          B_A(lo) = ival;
        }
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
          L_A(lo) = rval;
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
          U_A(lo) = rval;
        } else {
          D(i) = rval;
        }
      }
      for ( indexType i = 0; i < PRECO::pr_size; ++i )
        SPARSETOOL_ASSERT( D(i) != valueType(0),
                            "CSSORpreconditioner::D(" << i << ") = " << D(i) << " size = " << D.size() );

      // step 6: allocate working vectors
      br . resize( PRECO::pr_size );
      bi . resize( PRECO::pr_size );
      x  . resize( PRECO::pr_size );
      y  . resize( PRECO::pr_size );
    }

  public:

    CSSORpreconditioner(void) : Preco<CSSORPRECO>() {}
    
    template <typename MAT>
    CSSORpreconditioner( MAT const & M, valueType _omega, indexType _maxIter ) : Preco<CSSORPRECO>()
    { build_SOR( M, _omega, _maxIter ); }

    //! build the preconditioner from matrix \c M with pattern \c P
    template <typename MAT>
    void
    build( MAT const & M, valueType _omega, indexType _maxIter )
    { build_SOR( M, _omega, _maxIter ); }

    //! apply preconditioner to vector \c v and store result to vector \c res
    template <typename VECTOR>
    void
    assPreco( VECTOR & xc, VECTOR const & bc ) const {
      indexType const * pC;  indexType const * pI;  valueType const * pUA;
      indexType const * pBR; indexType const * pBJ; valueType const * pBA;
      indexType const * pR;  indexType const * pJ;  valueType const * pLA;
 
      indexType k;

      // copia dato in ingresso
      for ( k=0; k < PRECO::pr_size; ++k ) {
        br(k) = bc(k).real();
        bi(k) = bc(k).imag();
      }

      x.setZero();
      y.setZero();
      for ( indexType ii = 0; ii < maxIter; ++ii ) {

        // calcolo ((1/omega-1)*D-U)*x + B*y + b --------------------
        pC  = & U_C.front();
        pI  = & U_I.front();
        pUA = & U_A.front();
        for ( k=0; k < PRECO::pr_size; ++k ) {
          valueType xk = x(k);
          x(k) = br(k) + omega1 * D(k) * xk;
          for ( indexType i_cnt = pC[1] - pC[0]; i_cnt > 0; --i_cnt )
            x(*pI++) -= *pUA++ * xk;
          ++pC;
        }

        // aggiungo B*y;
        pBR = & B_R.front();
        pBJ = & B_J.front();
        pBA = & B_A.front();
        for ( k=0; k < PRECO::pr_size; ++k ) {
          valueType tmp(0);
          for ( indexType i_cnt = pBR[1] - pBR[0]; i_cnt > 0; --i_cnt )
            tmp += *pBA++ * y(*pBJ++);
          x(k) += tmp;
          ++pBR;
        };

        // solve (D/omega+L) x(n+1) = x(n)
        pR  = & L_R.front();
        pJ  = & L_J.front();
        pLA = & L_A.front();
        for ( k=0; k < PRECO::pr_size; ++k ) {
          valueType tmp(0);
          for ( indexType i_cnt = pR[1] - pR[0]; i_cnt > 0; --i_cnt )
            tmp += *pLA++ * x(*pJ++);
          x(k) = omega*(x(k)-tmp)/D(k);
          ++pR;
        };

        // calcolo ((1/omega-1)*D-U)*y - B*x + c  ------------------------
        pC  = & U_C.front();
        pI  = & U_I.front();
        pUA = & U_A.front();
        for ( k=0; k < PRECO::pr_size; ++k ) {
          valueType yk = y(k);
          y(k) = bi(k) + omega1 * D(k) * yk;
          for ( indexType i_cnt = pC[1] - pC[0]; i_cnt > 0; --i_cnt )
            y(*pI++) -= *pUA++ * yk;
          ++pC;
        }

        // sottraggo B*x;
        pBR = & B_R.front();
        pBJ = & B_J.front();
        pBA = & B_A.front();
        for ( k=0; k < PRECO::pr_size; ++k ) {
          valueType tmp(0);
          for ( indexType i_cnt = pBR[1] - pBR[0]; i_cnt > 0; --i_cnt )
            tmp += *pBA++ * x(*pBJ++);
          y(k) -= tmp;
          ++pBR;
        };

        // solve (D/omega+L) y(n+1) = y(n)
        pR  = & L_R.front();
        pJ  = & L_J.front();
        pLA = & L_A.front();
        for ( k=0; k < PRECO::pr_size; ++k ) {
          valueType tmp(0);
          for ( indexType i_cnt = pR[1] - pR[0]; i_cnt > 0; --i_cnt )
            tmp += *pLA++ * y(*pJ++);
          y(k) = omega*(y(k)-tmp)/D(k);
          ++pR;
        };
      }

      for ( indexType ii = 0; ii < maxIter; ++ii ) {

        // calcolo ((1/omega-1)*D-L)*y - B*x + c -----------------------
        pR  = & L_R.front() + PRECO::pr_size;
        pJ  = & L_J.front() + *pR;
        pLA = & L_A.front() + *pR;
        k = PRECO::pr_size;
        do {
          --k; --pR;
          valueType tmp(0);
          for ( indexType i_cnt = pR[1] - pR[0]; i_cnt > 0; --i_cnt )
            tmp += *--pLA * y(*--pJ);
          y(k) = omega1*D(k)*y(k)+bi(k)-tmp;
        } while ( k > 0 );

        // sottraggo B*x;
        pBR = & B_R.front();
        pBJ = & B_J.front();
        pBA = & B_A.front();
        for ( k=0; k < PRECO::pr_size; ++k ) {
          valueType tmp(0);
          for ( indexType i_cnt = pBR[1] - pBR[0]; i_cnt > 0; --i_cnt )
            tmp += *pBA++ * x(*pBJ++);
          y(k) -= tmp;
          ++pBR;
        };

        // solve (D/omega+U) y(n+1) = y(n)
        pC  = & U_C.front() + PRECO::pr_size;
        pI  = & U_I.front() + *pC;
        pUA = & U_A.front() + *pC;
        k = PRECO::pr_size;
        do {
          --k; --pC;
          valueType yk = y(k)*omega/D(k);
          for ( indexType i_cnt = pC[1] - pC[0]; i_cnt > 0; --i_cnt )
            y(*--pI) -= *--pUA * yk;
          y(k) = yk;
        } while ( k > 0 );

        // calcolo ((1/omega-1)*D-L)*x + B*y + b ----------------------
        pR  = & L_R.front() + PRECO::pr_size;
        pJ  = & L_J.front() + *pR;
        pLA = & L_A.front() + *pR;
        k = PRECO::pr_size;
        do {
          --k; --pR;
          valueType tmp(0);
          for ( indexType i_cnt = pR[1] - pR[0]; i_cnt > 0; --i_cnt )
            tmp += *--pLA * x(*--pJ);
          x(k) = omega1*D(k)*x(k)+br(k)-tmp;
        } while ( k > 0 );

        // aggiungo B*y;
        pBR = & B_R.front();
        pBJ = & B_J.front();
        pBA = & B_A.front();
        for ( k=0; k < PRECO::pr_size; ++k ) {
          valueType tmp(0);
          for ( indexType i_cnt = pBR[1] - pBR[0]; i_cnt > 0; --i_cnt )
            tmp += *pBA++ * y(*pBJ++);
          x(k) += tmp;
          ++pBR;
        };

        // solve (D/omega+U) x(n+1) = x(n)
        pC  = & U_C.front() + PRECO::pr_size;
        pI  = & U_I.front() + *pC;
        pUA = & U_A.front() + *pC;
        k = PRECO::pr_size;
        do {
          --k; --pC;
          valueType xk = x(k)*omega/D(k);
          for ( indexType i_cnt = pC[1] - pC[0]; i_cnt > 0; --i_cnt )
            x(*--pI) -= *--pUA * xk;
          x(k) = xk;
        } while ( k > 0 );

      }

      // copia risulatato in uscita
      for ( k=0; k < PRECO::pr_size; ++k ) xc(k) = valueType(x(k),y(k));
    }

  };

  //! \cond NODOC
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,CSSORpreconditioner<TP> >
  operator / (Vector<T> const & v, CSSORpreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,CSSORpreconditioner<TP> >(v,P); }
  //! \endcond

}

namespace SparseToolLoad {
  using ::SparseTool::CSSORpreconditioner;
}

#endif
