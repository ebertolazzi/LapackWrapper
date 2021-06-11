#pragma once
#ifndef SPARSETOOL_ITERATIVE_PRECO_ILDUK_HH
#define SPARSETOOL_ITERATIVE_PRECO_ILDUK_HH

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

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    typedef ILDUKpreconditioner<T> ILDUKPRECO;
    typedef Preco<ILDUKPRECO>      PRECO;
    #endif

    typedef T real_type; //!< type of the element of the preconditioner

  private:

    Vector<integer>          Lnnz, Unnz;

    Vector<real_type>          W;
    Vector<real_type>          D;

    Vector<Vector<real_type> > L_A;
    Vector<Vector<integer> > L_J;

    Vector<Vector<real_type> > U_A;
    Vector<Vector<integer> > U_I;

    //! build incomplete LDU decomposition with specified pattern `P` 
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
        integer i = A.row();
        integer j = A.column();
        if      ( i > j ) ++Lnnz(i);
        else if ( i < j ) ++Unnz(j);
      }

      // step 1: allocate memory
      L_A.resize( PRECO::pr_size );
      L_J.resize( PRECO::pr_size );
      U_A.resize( PRECO::pr_size );
      U_I.resize( PRECO::pr_size );
      for ( integer i = 0; i < PRECO::pr_size; ++i ) {
        L_A(i).reserve(Lnnz(i));
        L_J(i).reserve(Lnnz(i));
        U_A(i).reserve(Unnz(i));
        U_I(i).reserve(Unnz(i));
      }

      D.resize( PRECO::pr_size );
      W.resize( PRECO::pr_size );

      D = real_type(1);
      W = real_type(0);

      // step 2: insert values
      for ( A.Begin(); A.End(); A.Next() ) {
        integer i = A.row();
        integer j = A.column();
        typename MAT::real_type const val = A.value(); // (i,j);
        if      ( i > j ) { L_A(i).push_back(val); L_J(i).push_back(j); }
        else if ( i < j ) { U_A(j).push_back(val); U_I(j).push_back(i); }
        else              D(i) = val;
      }

      // SORTING (da eliminare)
      //for ( integer i = 0; i < PRECO::pr_size; ++i ) {
      //  QuickSortI<real_type>( &L_J(i).front(), &L_A(i).front(), L_A(i).size() );
      //  QuickSortI<real_type>( &U_I(i).front(), &U_A(i).front(), U_A(i).size() );
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
      for ( integer k = 1; k < PRECO::pr_size; ++k ) {
        integer kk;
        Vector<real_type> & L_Ak = L_A(k);
        Vector<integer> & L_Jk = L_J(k);
        Vector<real_type> & U_Ak = U_A(k);
        Vector<integer> & U_Ik = U_I(k);

        // W = M21^T  ---- l^T = D^(-1)U^(-T) M21^T
        for ( kk = 0; kk < L_Ak.size(); ++kk ) W(L_Jk(kk)) = L_Ak(kk);

        // W = U^(-T) W ---- l^T = D^(-1)U^(-T) M21^T
        for ( kk = 0; kk < L_Ak.size(); ++kk ) {
          integer              j = L_Jk(kk);
          Vector<real_type> & U_Aj = U_A(j);
          Vector<integer> & U_Ij = U_I(j);
          real_type bf = 0;
          for ( integer jj = 0; jj < U_Aj.size(); ++jj ) bf += W(U_Ij(jj))*U_Aj(jj);
          W(j) -= bf;
        }
        // l^T = D^(-1) W;   W = 0 ---- l^T = D^(-1)U^(-T) M21^T
        for ( kk = 0; kk < L_Ak.size(); ++kk )
          { integer j = L_Jk(kk); L_Ak(kk) = W(j) / D(j); W(j) = 0; }

        // W = M12  ----  u = D^(-1)L^(-1) M12
        for ( kk = 0; kk < U_Ak.size(); ++kk ) W(U_Ik(kk)) = U_Ak(kk);

        // W = L^(-1) W  ----  u = D^(-1)L^(-1) M12
        for ( kk = 0; kk < U_Ak.size(); ++kk ) {
          integer              i = U_Ik(kk);
          Vector<real_type> & L_Ai = L_A(i);
          Vector<integer> & L_Ji = L_J(i);
          real_type bf = 0;
          for ( integer ii = 0; ii < L_Ai.size(); ++ii ) bf += W(L_Ji(ii))*L_Ai(ii);
          W(i) -= bf;
        }

        real_type bf = 0;
        for ( kk = 0; kk < L_Ak.size(); ++kk ) bf += L_Ak(kk) * W(L_Jk(kk));
        D(k) -= bf;

        for ( kk = 0; kk < U_Ak.size(); ++kk )
          { integer i = U_Ik(kk); U_Ak(kk) = W(i) / D(i); W(i) = 0; }

        SPARSETOOL_ASSERT( D(k) != real_type(0), "ILDUKpreconditioner found D(" << k << ") == 0!" );
      }
    }

  public:

    ILDUKpreconditioner(void) : Preco<ILDUKPRECO>() {}
    
    template <typename MAT>
    ILDUKpreconditioner( MAT const & M ) : Preco<ILDUKPRECO>() 
    { build_ILDU( M ); }

    //! build the preconditioner from matrix `M`.
    template <typename MAT>
    void
    build( MAT const & M )
    { build_ILDU(M); }

    //!
    //! Apply preconditioner to vector `v`
    //! and store result to vector `res`.
    //!
    template <typename VECTOR>
    void
    assPreco( VECTOR & res, VECTOR const & v ) const {
      res = v;
      integer k = 0;
      // solve L
      while ( ++k < PRECO::pr_size ) {
        integer i_cnt = L_A(k).size();
        real_type const * pA = &L_A(k).front();
        integer const * pJ = &L_J(k).front();
        real_type bf = 0;
        while ( i_cnt-- > 0 ) bf += *pA++ * res(*pJ++);
        res(k) -= bf;
      }

      // solve D
      for ( k = 0; k < PRECO::pr_size; ++k ) res(k) /= D(k);

      // solve U
      do {
        typename VECTOR::real_type resk = res(--k);
        integer i_cnt = U_A(k).size();
        real_type const * pA = &U_A(k).front();
        integer const * pI = &U_I(k).front();
        while ( i_cnt-- > 0 ) res(*pI++) -= *pA++ * resk;
      } while ( k > 1 );

    }

  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,ILDUKpreconditioner<TP> >
  operator / (Vector<T> const & v, ILDUKpreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,ILDUKpreconditioner<TP> >(v,P); }
  #endif

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace SparseToolLoad {
  using ::SparseTool::ILDUKpreconditioner;
}
#endif

#endif
