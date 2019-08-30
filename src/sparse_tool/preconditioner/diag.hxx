#ifndef SPARSETOOL_ITERATIVE_PRECO_DIAG_HH
#define SPARSETOOL_ITERATIVE_PRECO_DIAG_HH

using namespace std;

namespace SparseTool {

  /*
  // ######  ######  ######  #######  #####  #######
  // #     # #     # #     # #       #     # #     #
  // #     # #     # #     # #       #       #     #
  // #     # ######  ######  #####   #       #     #
  // #     # #       #   #   #       #       #     #
  // #     # #       #    #  #       #     # #     #
  // ######  #       #     # #######  #####  #######
  */
  /*! \class Dpreconditioner
      \brief Diagonal preconditioner class
   */
  template <typename T>
  class Dpreconditioner : public Preco<Dpreconditioner<T> > {

    // \cond NODOC
    typedef Dpreconditioner<T> DPRECO;
    typedef Preco<DPRECO>      PRECO;
    // \endcond
    typedef T valueType; //!< type of the element of the preconditioner

    Vector<valueType> D;

  public:

    Dpreconditioner(void) : Preco<Dpreconditioner<T> >() {}

    template <typename MAT>
    Dpreconditioner( MAT const & M ) : Preco<DPRECO>()
    { build( M ); }

    //! build the preconditioner from matrix \c A
    template <typename MAT>
    void
    build(Sparse<T,MAT> const & A) {
      D . resize( A.numRows() );
      D = valueType(1);
      for ( A.Begin(); A.End(); A.Next() ) {
        indexType i = A.row();
        indexType j = A.column();
        if ( i == j ) D(i) = A.value();
      }
    }

    //! return the diagonal of the diagonal preconditioner
    Vector<valueType> const & GetD(void) const { return D; }

    //! apply preconditioner to vector \c v and store result to vector \c res
    template <typename VECTOR>
    void
    assPreco( VECTOR & res, VECTOR const & v ) const {
      for ( indexType i = 0; i < D.size(); ++i ) res(i) = v(i) / D(i);
    }

  };

  //! \cond NODOC
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,Dpreconditioner<TP> >
  operator / (Vector<T> const & v, Dpreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,Dpreconditioner<TP> >(v,P); }
  //! \endcond
}

namespace SparseToolLoad {
  using ::SparseTool::Dpreconditioner;
}

#endif
