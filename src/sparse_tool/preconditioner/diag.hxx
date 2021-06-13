#pragma once
#ifndef SPARSETOOL_ITERATIVE_PRECO_DIAG_HH
#define SPARSETOOL_ITERATIVE_PRECO_DIAG_HH

namespace Sparse_tool {

  /*
  // ######  ######  ######  #######  #####  #######
  // #     # #     # #     # #       #     # #     #
  // #     # #     # #     # #       #       #     #
  // #     # ######  ######  #####   #       #     #
  // #     # #       #   #   #       #       #     #
  // #     # #       #    #  #       #     # #     #
  // ######  #       #     # #######  #####  #######
  */
  //!
  //! Diagonal preconditioner class.
  //!
  template <typename T>
  class Dpreconditioner : public Preco<Dpreconditioner<T> > {

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    typedef Dpreconditioner<T> DPRECO;
    typedef Preco<DPRECO>      PRECO;
    #endif
    typedef T real_type; //!< type of the element of the preconditioner

    Vector<real_type> D;

  public:

    Dpreconditioner(void) : Preco<Dpreconditioner<T> >() {}

    template <typename MAT>
    Dpreconditioner( MAT const & M ) : Preco<DPRECO>()
    { build( M ); }

    //!
    //! Build the preconditioner from matrix `A`
    //!
    template <typename MAT>
    void
    build(Sparse<T,MAT> const & A) {
      D.resize( A.numRows() );
      D.fill(1);
      for ( A.Begin(); A.End(); A.Next() ) {
        integer i = A.row();
        integer j = A.column();
        if ( i == j ) D(i) = A.value();
      }
    }

    //! return the diagonal of the diagonal preconditioner
    Vector<real_type> const & GetD(void) const { return D; }

    //!
    //! Apply preconditioner to vector `v`
    //! and store result to vector `res`.
    //!
    template <typename VECTOR>
    void
    assPreco( VECTOR & res, VECTOR const & v ) const {
      res = v.array()/D.array();
    }

    //!
    //! Apply preconditioner to vector `v`
    //! and store result to vector `res`.
    //!
    template <typename VECTOR>
    void
    assPreco( VECTOR & inout ) const {
      inout.array() /= D.array();
    }

  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,Dpreconditioner<TP> >
  operator / (Vector<T> const & v, Dpreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,Dpreconditioner<TP> >(v,P); }
  #endif
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::Dpreconditioner;
}
#endif

#endif
