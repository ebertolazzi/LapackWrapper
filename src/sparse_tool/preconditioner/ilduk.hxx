#ifndef SPARSETOOL_ITERATIVE_PRECO_ILDUK_HH
#define SPARSETOOL_ITERATIVE_PRECO_ILDUK_HH

using namespace std;

namespace SparseTool {

  /*
  //  ### #       ######  #     # #    # 
  //   #  #       #     # #     # #   #  
  //   #  #       #     # #     # #  #   
  //   #  #       #     # #     # ###    
  //   #  #       #     # #     # #  #   
  //   #  #       #     # #     # #   #  
  //  ### ####### ######   #####  #    #
  */
  //! Incomplete \c LDU preconditioner
  template <typename T>
  class ILDUKpreconditioner : public Preco<ILDUKpreconditioner<T> > {
  public:

    //! \cond NODOC
    typedef ILDUKpreconditioner<T> ILDUKPRECO;
    typedef Preco<ILDUKPRECO>      PRECO;

    //! \endcond
    typedef T valueType; //!< type of the element of the preconditioner

  private:

    Vector<indexType>          Lnnz, Unnz;

    Vector<valueType>          W;
    Vector<valueType>          D;

    Vector<Vector<valueType> > L_A;
    Vector<Vector<indexType> > L_J;

    Vector<Vector<valueType> > U_A;
    Vector<Vector<indexType> > U_I;

    //! build incomplete LDU decomposition with specified pattern \c P
    template <typename MAT>
    void
    build_ILDU( MAT const & A ) {

      SPARSETOOL_ASSERT(
        A.isOrdered(),
        "ILDUKpreconditioner::build_LDU pattern must be ordered before use"
      )
      SPARSETOOL_ASSERT(
        A.numRows() == A.numCols(),
        "ILDUKpreconditioner::build_LDU only square matrix allowed"
      )
      SPARSETOOL_ASSERT(
        A.numRows() > 0,
        "ILDUKpreconditioner::build_LDU empty matrix"
      )

      // step 0: count necessary memory
      PRECO::pr_size = A.numRows();
      Lnnz.resize( PRECO::pr_size );
      Unnz.resize( PRECO::pr_size );

      Lnnz = 0;
      Unnz = 0;

      for ( A.Begin(); A.End(); A.Next() ) {
        indexType i = A.row();
        indexType j = A.column();
        if      ( i > j ) ++Lnnz(i);
        else if ( i < j ) ++Unnz(j);
      }

      // step 1: allocate memory
      L_A.resize( PRECO::pr_size );
      L_J.resize( PRECO::pr_size );
      U_A.resize( PRECO::pr_size );
      U_I.resize( PRECO::pr_size );
      for ( indexType i = 0; i < PRECO::pr_size; ++i ) {
        L_A(i).reserve(Lnnz(i));
        L_J(i).reserve(Lnnz(i));
        U_A(i).reserve(Unnz(i));
        U_I(i).reserve(Unnz(i));
      }

      D.resize( PRECO::pr_size );
      W.resize( PRECO::pr_size );

      D = valueType(1);
      W = valueType(0);

      // step 2: insert values
      for ( A.Begin(); A.End(); A.Next() ) {
        indexType i = A.row();
        indexType j = A.column();
        typename MAT::valueType const val = A.value(); // (i,j);
        if      ( i > j ) { L_A(i).push_back(val); L_J(i).push_back(j); }
        else if ( i < j ) { U_A(j).push_back(val); U_I(j).push_back(i); }
        else              D(i) = val;
      }

      // SORTING (da eliminare)
      //for ( indexType i = 0; i < PRECO::pr_size; ++i ) {
      //  QuickSortI<valueType>( &L_J(i).front(), &L_A(i).front(), L_A(i).size() );
      //  QuickSortI<valueType>( &U_I(i).front(), &U_A(i).front(), U_A(i).size() );
      //}

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
      for ( indexType k = 1; k < PRECO::pr_size; ++k ) {
        indexType kk;
        Vector<valueType> & L_Ak = L_A(k);
        Vector<indexType> & L_Jk = L_J(k);
        Vector<valueType> & U_Ak = U_A(k);
        Vector<indexType> & U_Ik = U_I(k);

        // W = M21^T  ---- l^T = D^(-1)U^(-T) M21^T
        for ( kk = 0; kk < L_Ak.size(); ++kk ) W(L_Jk(kk)) = L_Ak(kk);

        // W = U^(-T) W ---- l^T = D^(-1)U^(-T) M21^T
        for ( kk = 0; kk < L_Ak.size(); ++kk ) {
          indexType              j = L_Jk(kk);
          Vector<valueType> & U_Aj = U_A(j);
          Vector<indexType> & U_Ij = U_I(j);
          valueType bf = 0;
          for ( indexType jj = 0; jj < U_Aj.size(); ++jj ) bf += W(U_Ij(jj))*U_Aj(jj);
          W(j) -= bf;
        }
        // l^T = D^(-1) W;   W = 0 ---- l^T = D^(-1)U^(-T) M21^T
        for ( kk = 0; kk < L_Ak.size(); ++kk )
          { indexType j = L_Jk(kk); L_Ak(kk) = W(j) / D(j); W(j) = 0; }

        // W = M12  ----  u = D^(-1)L^(-1) M12
        for ( kk = 0; kk < U_Ak.size(); ++kk ) W(U_Ik(kk)) = U_Ak(kk);

        // W = L^(-1) W  ----  u = D^(-1)L^(-1) M12
        for ( kk = 0; kk < U_Ak.size(); ++kk ) {
          indexType              i = U_Ik(kk);
          Vector<valueType> & L_Ai = L_A(i);
          Vector<indexType> & L_Ji = L_J(i);
          valueType bf = 0;
          for ( indexType ii = 0; ii < L_Ai.size(); ++ii ) bf += W(L_Ji(ii))*L_Ai(ii);
          W(i) -= bf;
        }

        valueType bf = 0;
        for ( kk = 0; kk < L_Ak.size(); ++kk ) bf += L_Ak(kk) * W(L_Jk(kk));
        D(k) -= bf;

        for ( kk = 0; kk < U_Ak.size(); ++kk )
          { indexType i = U_Ik(kk); U_Ak(kk) = W(i) / D(i); W(i) = 0; }

        SPARSETOOL_ASSERT( D(k) != valueType(0), "ILDUKpreconditioner found D(" << k << ") == 0!" );
      }
    }

  public:

    ILDUKpreconditioner(void) : Preco<ILDUKPRECO>() {}
    
    template <typename MAT>
    ILDUKpreconditioner( MAT const & M ) : Preco<ILDUKPRECO>() 
    { build_ILDU( M ); }

    //! build the preconditioner from matrix \c M
    template <typename MAT>
    void
    build( MAT const & M )
    { build_ILDU(M); }

    //! apply preconditioner to vector \c v and store result to vector \c res
    template <typename VECTOR>
    void
    assPreco( VECTOR & res, VECTOR const & v ) const {
      res = v;
      indexType k = 0;
      // solve L
      while ( ++k < PRECO::pr_size ) {
        indexType i_cnt = L_A(k).size();
        valueType const * pA = &L_A(k).front();
        indexType const * pJ = &L_J(k).front();
        valueType bf = 0;
        while ( i_cnt-- > 0 ) bf += *pA++ * res(*pJ++);
        res(k) -= bf;
      }

      // solve D
      for ( k = 0; k < PRECO::pr_size; ++k ) res(k) /= D(k);

      // solve U
      do {
        typename VECTOR::valueType resk = res(--k);
        indexType i_cnt = U_A(k).size();
        valueType const * pA = &U_A(k).front();
        indexType const * pI = &U_I(k).front();
        while ( i_cnt-- > 0 ) res(*pI++) -= *pA++ * resk;
      } while ( k > 1 );

    }

  };

  //! \cond NODOC
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,ILDUKpreconditioner<TP> >
  operator / (Vector<T> const & v, ILDUKpreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,ILDUKpreconditioner<TP> >(v,P); }
  //! \endcond

}

namespace SparseToolLoad {
  using ::SparseTool::ILDUKpreconditioner;
}

#endif
