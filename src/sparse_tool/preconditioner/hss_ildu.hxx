#pragma once
#ifndef SPARSETOOL_ITERATIVE_PRECO_HSS_ILDU_HH
#define SPARSETOOL_ITERATIVE_PRECO_HSS_ILDU_HH

namespace Sparse_tool {

  /*
  //  #     #  #####   #####  ### #       ######  #     #
  //  #     # #     # #     #  #  #       #     # #     #
  //  #     # #       #        #  #       #     # #     #
  //  #######  #####   #####   #  #       #     # #     #
  //  #     #       #       #  #  #       #     # #     #
  //  #     # #     # #     #  #  #       #     # #     #
  //  #     #  #####   #####  ### ####### ######   #####
  */
  //! Incomplete \c LDU preconditioner
  template <typename T>
  class HSS_ILDU_Preconditioner : public Preco<HSS_ILDU_Preconditioner<T> > {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    typedef HSS_ILDU_Preconditioner<T> HSS_ILDU_PRECO;
    typedef Preco<HSS_ILDU_PRECO>      PRECO;
    #endif

    typedef T real_type; //!< type of the elements of the preconditioner
    typedef typename T::value_type rreal_type; //!< type of the elements of the preconditioner

  private:

    integer neq;
    ILDUpreconditioner<rreal_type> preco;

    //! build incomplete LDU decomposition with specified pattern `P`
    template <typename MAT>
    void
    build_HSS_ILDU( MAT const & A ) {

      UTILS_ASSERT0(
        A.is_ordered(),
        "Sparse_tool: HSS_ILDU_Preconditioner::build_LDU\n"
        "pattern must be ordered before use\n"
      );
      UTILS_ASSERT0(
        A.nrows() == A.ncols(),
        "Sparse_tool: HSS_ILDU_Preconditioner::build_LDU\n"
        "only square matrix allowed\n"
      );
      UTILS_ASSERT0(
        A.nrows() > 0,
        "Sparse_tool: HSS_ILDU_Preconditioner::build_LDU\n"
        "empty matrix\n"
      );

      neq = A.nrows();
      CCoorMatrix<rreal_type> Amat(neq,neq,A.nnz());

      // insert values
      for ( A.Begin(); A.End(); A.Next() ) {
        integer i = A.row();
        integer j = A.column();
        Amat.insert(i,j) = A.value().real() + A.value().imag();
      }

      Amat.internal_order();
      preco.build(Amat);
    }

  public:

    HSS_ILDU_Preconditioner(void) : Preco<HSS_ILDU_PRECO>() {}

    template <typename MAT>
    HSS_ILDU_Preconditioner( MAT const & M ) : Preco<HSS_ILDU_PRECO>()
    { build_HSS_ILDU( M ); }

    //! build the preconditioner from matrix `M`.
    template <typename MAT>
    void
    build( MAT const & M )
    { build_HSS_ILDU( M ); }

    //!
    //! Apply preconditioner to vector `v`
    //! and store result to vector `res`.
    //!
    template <typename VECTOR>
    void
    ass_preco( VECTOR & _y, VECTOR const & v ) const {
      preco.ass_preco(_y,v);
      real_type cst = real_type(0.5,-0.5);
      for ( integer k=0; k < neq; ++k ) _y(k) *= cst;
    }
  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,HSS_ILDU_Preconditioner<TP> >
  operator / (Vector<T> const & v, HSS_ILDU_Preconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,HSS_ILDU_Preconditioner<TP> >(v,P); }
  #endif

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::HSS_ILDU_Preconditioner;
}
#endif

#endif
