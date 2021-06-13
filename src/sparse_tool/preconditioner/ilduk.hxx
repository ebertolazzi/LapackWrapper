#pragma once
#ifndef SPARSETOOL_ITERATIVE_PRECO_ILDUK_HH
#define SPARSETOOL_ITERATIVE_PRECO_ILDUK_HH

namespace Sparse_tool {

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

    Vector<integer>            Lnnz, Unnz;

    Vector<real_type>          W;
    Vector<real_type>          D;

    vector<Vector<real_type> > L_A;
    vector<Vector<integer> >   L_J;
    Vector<integer>            L_size;

    vector<Vector<real_type> > U_A;
    vector<Vector<integer> >   U_I;
    Vector<integer>            U_size;

    //!
    //! Build incomplete LDU decomposition with specified pattern `P`.
    //!
    template <typename MAT>
    void
    build_ILDU( MAT const & A ) {

      UTILS_ASSERT0(
        A.isOrdered(),
        "Sparse_tool: ILDUKpreconditioner::build_LDU\n"
        "pattern must be ordered before use\n"
      );
      UTILS_ASSERT0(
        A.numRows() == A.numCols(),
        "Sparse_tool: ILDUKpreconditioner::build_LDU\n"
        "only square matrix allowed\n"
      );
      UTILS_ASSERT0(
        A.numRows() > 0,
        "Sparse_tool: ILDUKpreconditioner::build_LDU\n"
        "empty matrix\n"
      );

      // step 0: count necessary memory
      PRECO::pr_size = A.numRows();
      Lnnz.resize( PRECO::pr_size );
      Unnz.resize( PRECO::pr_size );

      Lnnz.setZero();
      Unnz.setZero();

      for ( A.Begin(); A.End(); A.Next() ) {
        integer i = A.row();
        integer j = A.column();
        if      ( i > j ) ++Lnnz(i);
        else if ( i < j ) ++Unnz(j);
      }

      // step 1: allocate memory
      L_A.resize( PRECO::pr_size );
      L_J.resize( PRECO::pr_size );
      L_size.resize( PRECO::pr_size );
      L_size.setZero();
      U_A.resize( PRECO::pr_size );
      U_I.resize( PRECO::pr_size );
      U_size.resize( PRECO::pr_size );
      U_size.setZero();
      for ( integer i = 0; i < PRECO::pr_size; ++i ) {
        L_A[i].resize(Lnnz(i));
        L_J[i].resize(Lnnz(i));
        U_A[i].resize(Unnz(i));
        U_I[i].resize(Unnz(i));
      }

      D.resize( PRECO::pr_size );
      W.resize( PRECO::pr_size );

      D.fill(1);
      W.setZero();

      // step 2: insert values
      for ( A.Begin(); A.End(); A.Next() ) {
        integer i = A.row();
        integer j = A.column();
        typename MAT::real_type const val = A.value(); // (i,j);
        if      ( i > j ) { L_A[i](L_size(i)) = val; L_J[i](L_size(i)) = j; ++L_size(i); }
        else if ( i < j ) { U_A[j](U_size(j)) = val; U_I[j](U_size(j)) = i; ++U_size(j); }
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
        Vector<real_type> & L_Ak    = L_A[k];
        Vector<integer>   & L_Jk    = L_J[k];
        integer           & L_sizek = L_size(k);
        Vector<real_type> & U_Ak    = U_A[k];
        Vector<integer>   & U_Ik    = U_I[k];
        integer           & U_sizek = U_size(k);

        // W = M21^T  ---- l^T = D^(-1)U^(-T) M21^T
        for ( kk = 0; kk < L_sizek; ++kk ) W(L_Jk(kk)) = L_Ak(kk);

        // W = U^(-T) W ---- l^T = D^(-1)U^(-T) M21^T
        for ( kk = 0; kk < L_sizek; ++kk ) {
          integer                j = L_Jk(kk);
          Vector<real_type> & U_Aj = U_A[j];
          Vector<integer>   & U_Ij = U_I[j];
          real_type bf = 0;
          for ( integer jj = 0; jj < U_Aj.size(); ++jj ) bf += W(U_Ij(jj))*U_Aj(jj);
          W(j) -= bf;
        }
        // l^T = D^(-1) W;   W = 0 ---- l^T = D^(-1)U^(-T) M21^T
        for ( kk = 0; kk < L_sizek; ++kk )
          { integer j = L_Jk(kk); L_Ak(kk) = W(j) / D(j); W(j) = 0; }

        // W = M12  ----  u = D^(-1)L^(-1) M12
        for ( kk = 0; kk < U_sizek; ++kk ) W(U_Ik(kk)) = U_Ak(kk);

        // W = L^(-1) W  ----  u = D^(-1)L^(-1) M12
        for ( kk = 0; kk < U_sizek; ++kk ) {
          integer                i = U_Ik(kk);
          Vector<real_type> & L_Ai = L_A[i];
          Vector<integer>   & L_Ji = L_J[i];
          real_type bf = 0;
          for ( integer ii = 0; ii < L_Ai.size(); ++ii ) bf += W(L_Ji(ii))*L_Ai(ii);
          W(i) -= bf;
        }

        real_type bf = 0;
        for ( kk = 0; kk < L_sizek; ++kk ) bf += L_Ak(kk) * W(L_Jk(kk));
        D(k) -= bf;

        for ( kk = 0; kk < U_sizek; ++kk )
          { integer i = U_Ik(kk); U_Ak(kk) = W(i) / D(i); W(i) = 0; }

        UTILS_ASSERT(
          D(k) != real_type(0),
          "Sparse_tool: ILDUKpreconditioner found D({}) == 0!\n", k
        );
      }
    }

  public:

    ILDUKpreconditioner(void) : Preco<ILDUKPRECO>() {}
    
    template <typename MAT>
    ILDUKpreconditioner( MAT const & M ) : Preco<ILDUKPRECO>() 
    { build_ILDU( M ); }

    //!
    //! Build the preconditioner from matrix `M`.
    //!
    template <typename MAT>
    void
    build( MAT const & M )
    { build_ILDU(M); }

    //!
    //! Apply preconditioner to vector `v`.
    //!
    template <typename VECTOR>
    void
    assPreco( VECTOR & v ) const {
      integer k = 0;
      // solve L
      while ( ++k < PRECO::pr_size ) {
        integer i_cnt = L_size(k);
        real_type const * pA = L_A[k].data();
        integer   const * pJ = L_J[k].data();
        real_type bf = 0;
        while ( i_cnt-- > 0 ) bf += *pA++ * v(*pJ++);
        v(k) -= bf;
      }

      // solve D
      v /= D;

      // solve U
      do {
        typename VECTOR::real_type vk = v(--k);
        integer i_cnt = U_size(k);
        real_type const * pA = U_A[k].data();
        integer   const * pI = U_I[k].data();
        while ( i_cnt-- > 0 ) v(*pI++) -= *pA++ * vk;
      } while ( k > 1 );

    }

    //!
    //! Apply preconditioner to vector `v`
    //! and store result to vector `res`.
    //!
    template <typename VECTOR>
    void
    assPreco( VECTOR & res, VECTOR const & v ) const {
      res = v;
      assPreco( res );
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
namespace Sparse_tool_load {
  using ::Sparse_tool::ILDUKpreconditioner;
}
#endif

#endif
