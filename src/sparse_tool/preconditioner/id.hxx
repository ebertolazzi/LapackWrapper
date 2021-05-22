#pragma once
#ifndef SPARSETOOL_ITERATIVE_PRECO_ID_HH
#define SPARSETOOL_ITERATIVE_PRECO_ID_HH

namespace SparseTool {

  /*
  //  ###        ######  ######  #######  #####  ####### 
  //   #  #####  #     # #     # #       #     # #     # 
  //   #  #    # #     # #     # #       #       #     # 
  //   #  #    # ######  ######  #####   #       #     # 
  //   #  #    # #       #   #   #       #       #     # 
  //   #  #    # #       #    #  #       #     # #     # 
  //  ### #####  #       #     # #######  #####  ####### 
  */
  //! Identity matrix preconditioner
  template <typename T>
  class IdPreconditioner : public Preco<IdPreconditioner<T> > {

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    typedef IdPreconditioner<T> IDPRECO;
    typedef Preco<IDPRECO>      PRECO;
    #endif

    typedef T valueType; //!< type of the element of the preconditioner

  public:

    IdPreconditioner(void) : Preco<IdPreconditioner<T> >() {}

    //! build the preconditioner from matrix \c A
    template <typename MAT>
    void build(Sparse<T,MAT> const &) {}

    //! apply preconditioner to vector \c v and store result to vector \c res
    template <typename VECTOR>
    void assPreco( VECTOR & res, VECTOR const & v ) const { res = v; }

  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,IdPreconditioner<TP> >
  operator / (Vector<T> const & v, IdPreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,IdPreconditioner<TP> >(v,P); }
  #endif

}

namespace SparseToolLoad {
  using ::SparseTool::IdPreconditioner;
}

#endif
