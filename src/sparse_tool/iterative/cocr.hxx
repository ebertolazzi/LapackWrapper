#pragma once
#ifndef SPARSETOOL_ITERATIVE_COCR_dot_HH
#define SPARSETOOL_ITERATIVE_COCR_dot_HH

namespace Sparse_tool {

  /*
  //   ####   ####   ####  #####
  //  #    # #    # #    # #    #
  //  #      #    # #      #    #
  //  #      #    # #      #####
  //  #    # #    # #    # #   #
  //   ####   ####   ####  #    #
  */
  //!
  //! Preconditioned Conjugate Gradient Iterative Solver.
  //!
  //! \param A       coefficient matrix
  //! \param b       righ hand side
  //! \param x       guess and solution
  //! \param P       preconditioner
  //! \param epsi    Admitted tolerance
  //! \param maxIter maximum number of admitted iteration
  //! \param iter    total number of performed itaration
  //! \param pStream pointer to stream object for messages
  //! \return last computed residual
  //!
  //! Use preconditioned conjugate gradient (for complex symmetric system)
  //! to solve \f$ A x = b \f$.
  //!
  template <
    typename real_type,
    typename matrix_type,
    typename vector_type,
    typename preco_type
  >
  real_type
  cocr(
    matrix_type const & A,
    vector_type const & b,
    vector_type       & x,
    preco_type  const & P,
    real_type   const & epsi,
    integer             maxIter,
    integer           & iter,
    ostream_type      * pStream = nullptr
  ) {

    UTILS_ASSERT(
      A.nrows() == b.size() &&
      A.ncols() == x.size() &&
      A.nrows() == A.ncols(),
      "Sparse_tool::cocr, bad system:\n"
      "dim matrix  = {} x {}\n"
      "dim r.h.s.  = {}\n"
      "dim unknown = {}\n",
      A.nrows(), A.ncols(), b.size(), x.size()
    );

    using vType = typename vector_type::real_type;

    integer     neq = b.size();
    vector_type p(neq), q(neq), qt(neq), r(neq), rt(neq), Art(neq);
    real_type   resid;
    vType       rho, mu, alpha, beta;

    r   = b - A * x;
    rt  = r / P;
    p   = rt;
    q   = A * p;
    rho = rt.conjugate().dot(q);

    iter = 1;
    do {

      resid = rt.template lpNorm<Eigen::Infinity>();
      if ( pStream != nullptr )
        fmt::print( *pStream,
          "[cocr] iter = {:4} residual = {:.6}\n", iter, resid
        );

      if ( resid <= epsi ) break;

      qt = q / P;
      mu = qt.conjugate().dot(q);
      UTILS_ASSERT0(
        mu != vType(0),
        "Sparse_tool: COCR failed found mu == 0\n"
      );

      alpha = rho/mu;
      x   += alpha * p;
      r   -= alpha * q;
      rt  -= alpha * qt;
      Art  = A * rt;

      beta = rho;
      rho  = rt.conjugate().dot(Art);
      beta = rho/beta;
      UTILS_ASSERT0(
        beta != vType(0),
        "Sparse_tool: COCR failed found beta == 0\n"
      );

      p = rt + beta * p;
      q = Art + beta * q;

    }  while ( ++iter <= maxIter );

    return resid;
  }

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::cocr;
}
#endif

#endif

/*
// ####### ####### #######
// #       #     # #
// #       #     # #
// #####   #     # #####
// #       #     # #
// #       #     # #
// ####### ####### #
*/
