#pragma once
#ifndef SPARSETOOL_ITERATIVE_PRECO_ID_HH
#define SPARSETOOL_ITERATIVE_PRECO_ID_HH

namespace Sparse_tool {

  /*
  //  ###        ######  ######  #######  #####  ####### 
  //   #  #####  #     # #     # #       #     # #     # 
  //   #  #    # #     # #     # #       #       #     # 
  //   #  #    # ######  ######  #####   #       #     # 
  //   #  #    # #       #   #   #       #       #     # 
  //   #  #    # #       #    #  #       #     # #     # 
  //  ### #####  #       #     # #######  #####  ####### 
  */

  //!
  //! Identity matrix preconditioner.
  //!
  template <typename T>
  class IdPreconditioner : public Preco<IdPreconditioner<T> > {

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    typedef IdPreconditioner<T> IDPRECO;
    typedef Preco<IDPRECO>      PRECO;
    #endif

    typedef T real_type; //!< type of the element of the preconditioner

  public:

    IdPreconditioner(void) : Preco<IdPreconditioner<T> >() {}

    template <typename MAT>
    IdPreconditioner( Sparse<T,MAT> const & ) : IdPreconditioner<T>() {}

    //!
    //! Build the preconditioner from matrix `A`.
    //!
    template <typename MAT>
    void build( Sparse<T,MAT> const & ) {}

    //!
    //! Apply preconditioner to vector `v`
    //! and store result to vector `res`.
    //!
    template <typename VECTOR>
    void ass_preco( VECTOR & res, VECTOR const & v ) const { res = v; }

    //!
    //! Apply preconditioner to vector `v`
    //!
    template <typename VECTOR>
    void ass_preco( VECTOR & ) const { /* Nothing to do */ }
  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,IdPreconditioner<TP> >
  operator / (Vector<T> const & v, IdPreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,IdPreconditioner<TP> >(v,P); }
  #endif

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::IdPreconditioner;
}
#endif

#endif
