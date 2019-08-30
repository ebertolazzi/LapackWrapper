#ifndef SPARSETOOL_ITERATIVE_PRECO_RILDU_HH
#define SPARSETOOL_ITERATIVE_PRECO_RILDU_HH

using namespace std;

namespace SparseTool {

  /*
  //  ######  ### #       ######  #     # 
  //  #     #  #  #       #     # #     # 
  //  #     #  #  #       #     # #     # 
  //  ######   #  #       #     # #     # 
  //  #   #    #  #       #     # #     # 
  //  #    #   #  #       #     # #     # 
  //  #     # ### ####### ######   #####  
  */
  //! Incomplete \c LDU preconditioner for the real part obly of a complex matrix
  template <typename T>
  class RILDUpreconditioner : public Preco<RILDUpreconditioner<T> > {
  public:

    //! \cond NODOC
    typedef RILDUpreconditioner<T> RILDUPRECO;
    typedef Preco<RILDUPRECO>      PRECO;

    //! \endcond
    typedef T valueType; //!< type of the elements of the preconditioner

  private:

    Vector<indexType> L_R;
    Vector<indexType> L_J;
    Vector<valueType> L_A;

    Vector<indexType> U_C;
    Vector<indexType> U_I;
    Vector<valueType> U_A;

    Vector<valueType> W;
    Vector<valueType> D;

    Vector<indexType> Lnnz, Unnz;

    //! build incomplete LDU decomposition with specified pattern \c P
    template <typename MAT, typename PAT>
    void
    build_RILDU( MAT const & A, PAT const & P ) {

      SPARSETOOL_ASSERT(
        P.isOrdered(),
        "RILDUpreconditioner::build_LDU pattern must be ordered before use"
      )
      SPARSETOOL_ASSERT(
        P.numRows() == A.numRows() && P.numCols() == A.numCols(),
        "RILDUpreconditioner::build_LDU pattern do not match matrix size"
      )
      SPARSETOOL_ASSERT(
        P.numRows() == P.numCols(),
        "RILDUpreconditioner::build_LDU only square matrix allowed"
      )
      SPARSETOOL_ASSERT(
        P.numRows() > 0,
        "RILDUpreconditioner::build_LDU empty matrix"
      )

      // step 0: compute necessary memory
      PRECO::pr_size = A.numRows();
      Lnnz.resize( PRECO::pr_size );
      Unnz.resize( PRECO::pr_size );

      Lnnz = 0;
      Unnz = 0;

      for ( P.Begin(); P.End(); P.Next() ) {
        indexType i = P.row();
        indexType j = P.column();
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

      D.resize( PRECO::pr_size );
      W.resize( PRECO::pr_size );

      D   = valueType(1);
      W   = valueType(0);
      L_A = valueType(0);
      U_A = valueType(0);
      
      // step 3: fill structure
      for ( P.Begin(); P.End(); P.Next() ) {
        indexType i = P.row();
        indexType j = P.column();
        if      ( i > j ) { indexType ii = --Lnnz(i); L_J(L_R(i)+ii) = j; }
        else if ( i < j ) { indexType jj = --Unnz(j); U_I(U_C(j)+jj) = i; }
      }

      // step 4: sort structure
      for ( indexType i = 0; i < PRECO::pr_size; ++i ) {
        sort( &L_J(L_R(i)), &L_J(L_R(i+1)) );
        sort( &U_I(U_C(i)), &U_I(U_C(i+1)) );
      }

      // insert values
      for ( A.Begin(); A.End(); A.Next() ) {
        indexType i   = A.row();
        indexType j   = A.column();
        valueType val = A.value().real(); // (i,j);
        if ( i > j ) {
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
        } else if ( i < j ) {
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
        }
      }

      /*
      //
      //  +---+---+ +---+---+ +---+---+    +-----+-------+
      //  | L | 0 | | D | 0 | | U | u |    | LDU |  LDu  |
      //  +---+---+ +---+---+ +---+---+ =  +-----+-------+
      //  |l^T| 1 | | 0 | d | | 0 | 1 |    |l^TDU|d+l^TDu|
      //  +---+---+ +---+---+ +---+---+    +-----+-------+
      //
      //      +-----+-----+
      //      | M11 | M12 |
      //  M = +-----+-----+
      //      | M21 | M22 |
      //      +-----+-----+
      //
      //  l^T D U = M21 ==> U^T D^T l = M21^T
      //
      //  l^T = D^(-1)U^(-T) M21^T
      //  u   = D^(-1)L^(-1) M12
      //  d   = M22 - lDu
      //
      */

      // build LDU decomposition
      indexType LRk1 = L_R(1);
      indexType UCk1 = U_C(1);
      for ( indexType k = 1; k < PRECO::pr_size; ++k ) {
        indexType kk;
        indexType LRk = LRk1; LRk1 = L_R(k+1);
        indexType UCk = UCk1; UCk1 = U_C(k+1);

        //  W = M21^T  ---- l^T = D^(-1)U^(-T) M21^T
        for ( kk = LRk; kk < LRk1; ++kk ) W(L_J(kk)) = L_A(kk);

        //  W = U^(-T) W ---- l^T = D^(-1)U^(-T) M21^T
        #define LDU_FAST
        #ifdef LDU_FAST
        for ( kk = LRk; kk < LRk1; ++kk ) {
          indexType j = L_J(kk);
        #else
        for ( indexType j = 1; j < k; ++j ) {
        #endif
          indexType UCj  = U_C(j);
          indexType UCj1 = U_C(j+1);
          valueType bf   = 0;
          for ( indexType jj = UCj; jj < UCj1; ++jj ) bf += W(U_I(jj))*U_A(jj);
          W(j) -= bf;
        }
        //  l^T = D^(-1) W;   W = 0---- l^T = D^(-1)U^(-T) M21^T
        #ifdef LDU_FAST
        for ( kk = LRk; kk < LRk1; ++kk ) { indexType j = L_J(kk); L_A(kk) = W(j) / D(j); W(j) = 0; }
        #else
        for ( kk = LRk; kk < LRk1; ++kk ) L_A(kk) = W(L_J(kk)) / D(L_J(kk));
        W = 0;
        #endif
        
        //  W = M12  ---- u   = D^(-1)L^(-1) M12
        for ( kk = UCk; kk < UCk1; ++kk ) W(U_I(kk)) = U_A(kk);

        //  W = L^(-1) W  ---- u   = D^(-1)L^(-1) M12
        #ifdef LDU_FAST
        for ( kk = UCk; kk < UCk1; ++kk ) {
          indexType i = U_I(kk);
        #else
        for ( indexType i = 1; i < k; ++i ) {
        #endif
          indexType LRi  = L_R(i);
          indexType LRi1 = L_R(i+1);
          valueType bf = 0;
          for ( indexType ii = LRi; ii < LRi1; ++ii ) bf += W(L_J(ii))*L_A(ii);
          W(i) -= bf;
        }

        valueType bf = 0;
        for ( kk = LRk; kk < LRk1; ++kk ) bf += L_A(kk) * W(L_J(kk));
        D(k) -= bf;

        #ifdef LDU_FAST
        for ( kk = UCk; kk < UCk1; ++kk ) { indexType i = U_I(kk); U_A(kk) = W(i) / D(i); W(i) = 0; }
        #else
        for ( kk = UCk; kk < UCk1; ++kk ) U_A(kk) = W(U_I(kk)) / D(U_I(kk));
        W = 0;
        #endif

        SPARSETOOL_ASSERT( D(k) != valueType(0), "ILDUpreconditioner found D(" << k << ") == 0!" );
      }
    }

  public:

    RILDUpreconditioner(void) : Preco<RILDUPRECO>() {}
    
    template <typename MAT>
    RILDUpreconditioner( MAT const & M ) : Preco<RILDUPRECO>() 
    { build_RILDU( M, M ); }

    template <typename MAT, typename PRE>
    RILDUpreconditioner( MAT const & M, PRE const & P ) : Preco<RILDUPRECO>() 
    { build_RILDU(M,P); }

    //! build the preconditioner from matrix \c M
    template <typename MAT>
    void
    build( MAT const & M )
    { build_RILDU(M,M); }

    //! build the preconditioner from matrix \c M with pattern \c P
    template <typename MAT, typename PRE>
    void
    build( MAT const & M, PRE const & P )
    { build_RILDU(M,P); }

    //! apply preconditioner to vector \c v and store result to vector \c res
    template <typename VECTOR>
    void
    assPreco( VECTOR & res, VECTOR const & v ) const {
      res = v;

      // solve L
      indexType const * pR  = & L_R.front();
      indexType const * pJ  = & L_J.front();
      valueType const * pLA = & L_A.front();
      indexType k;

      for ( k=1; k < PRECO::pr_size; ++k ) {
        ++pR;
        typename VECTOR::valueType tmp(0);
        for ( indexType i_cnt = pR[1] - pR[0]; i_cnt > 0; --i_cnt )
          tmp += *pLA++ * res(*pJ++);
        res(k) -= tmp;
      };

      // solve D
      for ( k = 0; k < PRECO::pr_size; ++k ) res(k) /= D(k);

      // solve U
      indexType const * pC  = & U_C.front() + PRECO::pr_size;
      indexType const * pI  = & U_I.front() + *pC;
      valueType const * pUA = & U_A.front() + *pC;

      do {
        typename VECTOR::valueType resk = res(--k);
        --pC;
        for ( indexType i_cnt = pC[1] - pC[0]; i_cnt > 0; --i_cnt )
          res(*--pI) -= *--pUA * resk;
      } while ( k > 1 );

    }

  };

  //! \cond NODOC
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,RILDUpreconditioner<TP> >
  operator / (Vector<T> const & v, RILDUpreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,RILDUpreconditioner<TP> >(v,P); }
  //! \endcond

}

namespace SparseToolLoad {
  using ::SparseTool::RILDUpreconditioner;
}

#endif
