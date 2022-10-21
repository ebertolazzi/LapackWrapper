#pragma once

#ifndef SPARSETOOL_ITERATIVE_PRECO_HSS_SSOR_HH
#define SPARSETOOL_ITERATIVE_PRECO_HSS_SSOR_HH

namespace Sparse_tool {

  /*
  //    #     #  #####   #####  ####### ######
  //    #     # #     # #     # #     # #     #
  //    #     # #       #       #     # #     #
  //    #######  #####   #####  #     # ######
  //    #     #       #       # #     # #   #
  //    #     # #     # #     # #     # #    #
  //    #     #  #####   #####  ####### #     #
  */
  //! Incomplete \c LDU preconditioner
  template <typename T>
  class HSS_SSOR_Preconditioner : public Preco<HSS_SSOR_Preconditioner<T> > {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    using HSS_SSOR_PRECO = HSS_SSOR_Preconditioner<T>;
    using PRECO          = Preco<HSS_SSOR_PRECO>;
    #endif

    using real_type  = T; //!< type of the elements of the preconditioner
    using rreal_type = typename T::value_type; //!< type of the elements of the preconditioner

  private:

    integer neq;
    SSORpreconditioner<rreal_type> preco;

    //! build incomplete LDU decomposition with specified pattern `P`
    template <typename MAT>
    void
    build_HSSOR( MAT const & A, rreal_type const & omega, integer m ) {

      UTILS_ASSERT0(
        A.is_ordered(),
        "Sparse_tool: HSS_SSOR_Preconditioner::build_LDU\n"
        "pattern must be ordered before use\n"
      );
      UTILS_ASSERT0(
        A.nrows() == A.ncols(),
        "Sparse_tool: HSS_SSOR_Preconditioner::build_LDU\n"
        "only square matrix allowed\n"
      );
      UTILS_ASSERT0(
        A.nrows() > 0,
        "Sparse_tool: HSS_SSOR_Preconditioner::build_LDU\n"
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
      preco.build(Amat,omega,m);
    }

  public:

    HSS_SSOR_Preconditioner(void) : Preco<HSS_SSOR_PRECO>() {}

    template <typename MAT>
    HSS_SSOR_Preconditioner( MAT const & M, rreal_type const & omega, integer m ) : Preco<HSS_SSOR_PRECO>()
    { build_HSSOR( M, omega, m ); }

    //! build the preconditioner from matrix `M`.
    template <typename MAT>
    void
    build( MAT const & M, rreal_type const & omega, integer m )
    { build_HSSOR( M, omega, m ); }

    //! apply preconditioner to vector `v`  and store result to vector `res`    template <typename VECTOR>
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
  Vector_V_div_P<Vector<T>,HSS_SSOR_Preconditioner<TP> >
  operator / (Vector<T> const & v, HSS_SSOR_Preconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,HSS_SSOR_Preconditioner<TP> >(v,P); }
  #endif

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::HSS_SSOR_Preconditioner;
}
#endif

#endif
