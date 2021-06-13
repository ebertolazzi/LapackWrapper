#pragma once
#ifndef SPARSETOOL_ITERATIVE_PRECO_SOR_HH
#define SPARSETOOL_ITERATIVE_PRECO_SOR_HH

namespace Sparse_tool {

  /*
  //   #####  ####### ######  
  //  #     # #     # #     # 
  //  #       #     # #     # 
  //   #####  #     # ######  
  //        # #     # #   #   
  //  #     # #     # #    #  
  //   #####  ####### #     # 
  */
  //!
  //! Iterative SOR preconditioner.
  //!
  template <typename T>
  class SORpreconditioner : public Preco<SORpreconditioner<T> > {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    typedef SORpreconditioner<T> SORPRECO;
    typedef Preco<SORPRECO>      PRECO;
    #endif

    typedef T real_type; //!< type of the elements of the preconditioner

  private:

    real_type omega;

    Vector<integer>   L_R, L_J;
    Vector<real_type> L_A;

    Vector<integer>   U_C, U_I;
    Vector<real_type> U_A;

    Vector<real_type> D;

    //!
    //! Build incomplete SOR preconditioner.
    //!
    template <typename MAT>
    void
    build_SOR( MAT const & A, real_type const & _omega ) {
      this -> omega = _omega;
      // step 0: compute necessary memory
      PRECO::pr_size = A.numRows();
      separate_LDU( A, L_A, L_R, L_J, D, U_A, U_I, U_C );
    }

  public:

    SORpreconditioner(void) : Preco<SORPRECO>() {}
    
    template <typename MAT>
    SORpreconditioner( MAT const & M, real_type _omega ) : Preco<SORPRECO>()
    { build_SOR( M, _omega ); }

    //!
    //! Build the preconditioner from matrix `M` with pattern `P`.
    //!
    template <typename MAT>
    void
    build( MAT const & M, real_type _omega )
    { build_SOR( M, _omega ); }

    //!
    //! Change omega.
    //!
    template <typename MAT>
    void
    change_omega( real_type const & _omega ) {
      this -> omega  = _omega;
      this -> omega1 = 1.0/_omega-1.0;
    }

    //!
    //! Apply preconditioner to vector `b`
    //! and store result to vector `x`.
    //!
    template <typename VECTOR>
    void
    assPreco( VECTOR & x, VECTOR const & b ) const {
      x = b;
      assPreco( x );
    }

    //!
    //! Apply preconditioner to vector `x`.
    //!
    template <typename VECTOR>
    void
    assPreco( VECTOR & x ) const {
      solve_DL( 1/omega, D, L_A, L_R, L_J, x );
      //solve_DU( 1/omega, D, U_A, U_I, U_C, x );
    }

  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,SORpreconditioner<TP> >
  operator / (Vector<T> const & v, SORpreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,SORpreconditioner<TP> >(v,P); }
  #endif

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::SORpreconditioner;
}
#endif

#endif
