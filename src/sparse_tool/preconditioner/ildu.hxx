#pragma once
#ifndef SPARSETOOL_ITERATIVE_PRECO_ILDU_HH
#define SPARSETOOL_ITERATIVE_PRECO_ILDU_HH

namespace Sparse_tool {

  /*
  //  ###  #       ######  #     #
  //   #   #       #     # #     #
  //   #   #       #     # #     #
  //   #   #       #     # #     #
  //   #   #       #     # #     #
  //   #   #       #     # #     #
  //  ###  ####### ######   #####
  */
  //! Incomplete \c LDU preconditioner
  template <typename T>
  class ILDUpreconditioner : public Preco<ILDUpreconditioner<T> > {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    typedef ILDUpreconditioner<T> ILDUPRECO;
    typedef Preco<ILDUPRECO>      PRECO;
    #endif

    typedef T real_type; //!< type of the elements of the preconditioner

  private:

    Vector<integer>   L_R;
    Vector<integer>   L_J;
    Vector<real_type> L_A;

    Vector<integer>   U_C;
    Vector<integer>   U_I;
    Vector<real_type> U_A;

    Vector<real_type> W;
    Vector<real_type> D;

    //!
    //! Build incomplete LDU decomposition with specified pattern `P`.
    //!
    template <typename MAT, typename PAT>
    void
    build_ILDU( MAT const & A, PAT const & P ) {

      UTILS_ASSERT0(
        P.is_ordered(),
        "Sparse_tool: ILDUpreconditioner::build_ILDU\n"
        "pattern must be ordered before use\n"
      );
      UTILS_ASSERT0(
        P.nrows() == A.nrows() && P.ncols() == A.ncols(),
        "Sparse_tool: ILDUpreconditioner::build_ILDU\n"
        "pattern do not match matrix size\n"
      );
      UTILS_ASSERT0(
        P.nrows() == P.ncols(),
        "Sparse_tool: ILDUpreconditioner::build_ILDU\n"
        "only square matrix allowed\n"
      );
      UTILS_ASSERT0(
        P.nrows() > 0,
        "Sparse_tool: ILDUpreconditioner::build_ILDU\n"
        "empty matrix\n"
      );

      // step 0: compute necessary memory
      PRECO::pr_size = A.nrows();

      separate_LDU( A, L_A, L_R, L_J, D, U_A, U_I, U_C );
      W.resize(PRECO::pr_size);
      W.setZero();

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
      integer LRk1 = L_R(1);
      integer UCk1 = U_C(1);
      for ( integer k = 1; k < PRECO::pr_size; ++k ) {
        integer kk;
        integer LRk = LRk1; LRk1 = L_R(k+1);
        integer UCk = UCk1; UCk1 = U_C(k+1);

        //  W = M21^T  ---- l^T = D^(-1)U^(-T) M21^T
        for ( kk = LRk; kk < LRk1; ++kk ) W(L_J(kk)) = L_A(kk);

        //  W = U^(-T) W ---- l^T = D^(-1)U^(-T) M21^T
        //#define LDU_FAST
        #ifdef LDU_FAST
        for ( kk = LRk; kk < LRk1; ++kk ) {
          integer j = L_J(kk);
        #else
        for ( integer j = 1; j < k; ++j ) {
        #endif
          integer UCj  = U_C(j);
          integer UCj1 = U_C(j+1);
          real_type bf   = 0;
          for ( integer jj = UCj; jj < UCj1; ++jj ) bf += W(U_I(jj))*U_A(jj);
          W(j) -= bf;
        }
        //  l^T = D^(-1) W;   W = 0---- l^T = D^(-1)U^(-T) M21^T
        #ifdef LDU_FAST
        for ( kk = LRk; kk < LRk1; ++kk ) { integer j = L_J(kk); L_A(kk) = W(j) / D(j); W(j) = 0; }
        #else
        for ( kk = LRk; kk < LRk1; ++kk ) L_A(kk) = W(L_J(kk)) / D(L_J(kk));
        W = 0;
        #endif

        //  W = M12  ---- u   = D^(-1)L^(-1) M12
        for ( kk = UCk; kk < UCk1; ++kk ) W(U_I(kk)) = U_A(kk);

        //  W = L^(-1) W  ---- u   = D^(-1)L^(-1) M12
        #ifdef LDU_FAST
        for ( kk = UCk; kk < UCk1; ++kk ) {
          integer i = U_I(kk);
        #else
        for ( integer i = 1; i < k; ++i ) {
        #endif
          integer LRi  = L_R(i);
          integer LRi1 = L_R(i+1);
          real_type bf = 0;
          for ( integer ii = LRi; ii < LRi1; ++ii ) bf += W(L_J(ii))*L_A(ii);
          W(i) -= bf;
        }

        real_type bf = 0;
        for ( kk = LRk; kk < LRk1; ++kk ) bf += L_A(kk) * W(L_J(kk));
        D(k) -= bf;

        #ifdef LDU_FAST
        for ( kk = UCk; kk < UCk1; ++kk ) { integer i = U_I(kk); U_A(kk) = W(i) / D(i); W(i) = 0; }
        #else
        for ( kk = UCk; kk < UCk1; ++kk ) U_A(kk) = W(U_I(kk)) / D(U_I(kk));
        W = 0;
        #endif

        UTILS_ASSERT(
          D(k) != real_type(0),
          "Sparse_tool: ILDUpreconditioner found D({}) == 0!\n", k
        );
      }
    }

  public:

    ILDUpreconditioner(void) : Preco<ILDUPRECO>() {}

    template <typename MAT>
    ILDUpreconditioner( MAT const & M ) : Preco<ILDUPRECO>()
    { build_ILDU( M, M ); }

    template <typename MAT, typename PRE>
    ILDUpreconditioner( MAT const & M, PRE const & P ) : Preco<ILDUPRECO>()
    { build_ILDU(M,P); }

    //!
    //! Build the preconditioner from matrix `M`.
    //!
    template <typename MAT>
    void
    build( MAT const & M )
    { build_ILDU(M,M); }

    //!
    //! Build the preconditioner from matrix `M` with pattern `P`.
    //!
    template <typename MAT, typename PRE>
    void
    build( MAT const & M, PRE const & P )
    { build_ILDU(M,P); }

    //!
    //! Apply preconditioner to vector `v`
    //! and store result to vector `res`.
    //!
    template <typename VECTOR>
    void
    ass_preco( VECTOR & v ) const {
      // solve L
      integer   const * pR  = L_R.data();
      integer   const * pJ  = L_J.data();
      real_type const * pLA = L_A.data();
      integer k;

      for ( k=1; k < PRECO::pr_size; ++k ) {
        ++pR;
        typename VECTOR::real_type tt(0);
        for ( integer i_cnt = pR[1] - pR[0]; i_cnt > 0; --i_cnt )
          tt += *pLA++ * v(*pJ++);
        v(k) -= tt;
      };

      // solve D
      v /= D;

      // solve U
      integer   const * pC  = U_C.data() + PRECO::pr_size;
      integer   const * pI  = U_I.data() + *pC;
      real_type const * pUA = U_A.data() + *pC;

      do {
        typename VECTOR::real_type vk = v(--k);
        --pC;
        for ( integer i_cnt = pC[1] - pC[0]; i_cnt > 0; --i_cnt )
          v(*--pI) -= *--pUA * vk;
      } while ( k > 1 );

    }

    //!
    //! Apply preconditioner to vector `v`
    //! and store result to vector `res`.
    //!
    template <typename VECTOR>
    void
    ass_preco( VECTOR & res, VECTOR const & v ) const {
      res = v;
      this->ass_preco( res );
    }

  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,ILDUpreconditioner<TP> >
  operator / (Vector<T> const & v, ILDUpreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,ILDUpreconditioner<TP> >(v,P); }
  #endif

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::ILDUpreconditioner;
}
#endif

#endif
