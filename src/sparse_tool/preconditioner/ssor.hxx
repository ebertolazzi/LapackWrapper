#pragma once
#ifndef SPARSETOOL_ITERATIVE_PRECO_SSOR_HH
#define SPARSETOOL_ITERATIVE_PRECO_SSOR_HH

namespace Sparse_tool {

  /*
  //   #####   #####  ####### ######
  //  #     # #     # #     # #     #
  //  #       #       #     # #     #
  //   #####   #####  #     # ######
  //        #       # #     # #   #
  //  #     # #     # #     # #    #
  //   #####   #####  ####### #     #
  */
  //! Iterative SSOR preconditioner
  template <typename T>
  class SSORpreconditioner : public Preco<SSORpreconditioner<T> > {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    typedef SSORpreconditioner<T> SSORPRECO;
    typedef Preco<SSORPRECO>      PRECO;
    #endif

    typedef T real_type; //!< type of the elements of the preconditioner

  private:

    real_type omega;

    Vector<integer>   L_R;
    Vector<integer>   L_J;
    Vector<real_type> L_A;

    Vector<integer>   U_C;
    Vector<integer>   U_I;
    Vector<real_type> U_A;

    Vector<real_type> D;

    //!
    //! Build incomplete SSOR preconditioner.
    //!
    template <typename MAT>
    void
    build_SSOR( MAT const & A, real_type const & _omega ) {
      this -> omega = _omega;
      // step 0: compute necessary memory
      PRECO::pr_size = A.nrows();
      separate_LDU( A, L_A, L_R, L_J, D, U_A, U_I, U_C );
    }

  public:

    SSORpreconditioner(void) : Preco<SSORPRECO>() {}

    template <typename MAT>
    SSORpreconditioner( MAT const & M, real_type _omega ) : Preco<SSORPRECO>()
    { build_SSOR( M, _omega ); }

    //!
    //! Build the preconditioner from matrix `M` with pattern `P`.
    //!
    template <typename MAT>
    void
    build( MAT const & M, real_type _omega )
    { build_SSOR( M, _omega ); }

    //!
    //! Change omega.
    //!
    template <typename MAT>
    void
    change_omega( real_type const & _omega ) {
      this -> omega = _omega;
    }

    //!
    //! Apply preconditioner to vector `b`
    //! and store result to vector `x`.
    //!
    template <typename VECTOR>
    void
    ass_preco( VECTOR & x, VECTOR const & b ) const {
      x = b;
      ass_preco( x );
    }

    //!
    //! Apply preconditioner to vector `x`.
    //!
    template <typename VECTOR>
    void
    ass_preco( VECTOR & x ) const {
      solve_DL( 1/omega, D, L_A, L_R, L_J, x );
      x.array() /= ((2-omega)/omega)*D.array();
      solve_DU( 1/omega, D, U_A, U_I, U_C, x );
    }

  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,SSORpreconditioner<TP> >
  operator / (Vector<T> const & v, SSORpreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,SSORpreconditioner<TP> >(v,P); }
  #endif

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::SSORpreconditioner;
}
#endif

#endif
