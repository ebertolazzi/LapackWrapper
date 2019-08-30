#ifndef SPARSETOOL_ITERATIVE_PRECO_ID_HH
#define SPARSETOOL_ITERATIVE_PRECO_ID_HH

using namespace std;

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

    // \cond NODOC
    typedef IdPreconditioner<T> IDPRECO;
    typedef Preco<IDPRECO>      PRECO;
    // \endcond
    typedef T valueType; //!< type of the element of the preconditioner

  public:

    IdPreconditioner(void) : Preco<IdPreconditioner<T> >() {}

    //! build the preconditioner from matrix \c A
    template <typename MAT>
    void build(Sparse<T,MAT> const & A) {}

    //! apply preconditioner to vector \c v and store result to vector \c res
    template <typename VECTOR>
    void assPreco( VECTOR & res, VECTOR const & v ) const { res = v; }

  };

  //! \cond NODOC
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,IdPreconditioner<TP> >
  operator / (Vector<T> const & v, IdPreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,IdPreconditioner<TP> >(v,P); }
  //! \endcond

}

namespace SparseToolLoad {
  using ::SparseTool::IdPreconditioner;
}

#endif
