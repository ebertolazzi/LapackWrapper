#pragma once
#ifndef SPARSETOOL_ITERATIVE_CG_dot_HH
#define SPARSETOOL_ITERATIVE_CG_dot_HH

namespace Sparse_tool {

  /*
  //  #####   #####
  // #     # #     #
  // #       #
  // #       #  ####
  // #       #     #
  // #     # #     #
  //  #####   #####
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
  //! Use preconditioned conjugate gradient to solve \f$ A x = b \f$.
  //!
  template <
    typename real_type,
    typename matrix_type,
    typename vector_type,
    typename preco_type
  >
  real_type
  cg(
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
      A.numRows() == b.size() &&
      A.numCols() == x.size() &&
      A.numRows() == A.numCols(),
      "Sparse_tool::cg, bad system:\n"
      "dim matrix  = {} x {}\n"
      "dim r.h.s.  = {}\n"
      "dim unknown = {}\n",
      A.numRows(), A.numCols(), b.size(), x.size()
    );

    typedef typename vector_type::real_type vType;

    integer     neq = b.size();
    vector_type p(neq), q(neq), r(neq), Ap(neq);
    real_type   resid;
    vType       rho, rho_1;

    r   = b - A * x;
    p   = r / P;
    rho = p.dot(r);

    iter = 1;
    do {

      resid = r.template lpNorm<Eigen::Infinity>();
      if ( pStream != nullptr )
        fmt::print( *pStream,
          "[cg] iter = {:4} residual = {:.6}\n", iter, resid
        );

      if ( resid <= epsi ) break;

      Ap  = A * p;
      vType alpha = rho/Ap.dot(p);
      x += alpha * p;
      r -= alpha * Ap;
      q  = r / P;

      rho_1 = rho;
      rho   = q.dot(r);

      p = q + (rho / rho_1) * p;

    }  while ( ++iter <= maxIter );

    return resid;
  }

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::cg;
}
#endif

#endif
