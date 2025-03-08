#pragma once
#ifndef SPARSETOOL_ITERATIVE_BICGSTAB_dot_HH
#define SPARSETOOL_ITERATIVE_BICGSTAB_dot_HH

namespace Sparse_tool {

  /*
  // ######       #####   #####   #####  #######    #    #####
  // #     #  #  #     # #     # #     #    #      # #   #     #
  // #     #  #  #       #       #          #     #   #  #     #
  // ######   #  #       #  ####  #####     #    #     # ######
  // #     #  #  #       #     #       #    #    ####### #     #
  // #     #  #  #     # #     # #     #    #    #     # #     #
  // ######   #   #####   #####   #####     #    #     # ######
  */
  /*!
   *  Bi Conjugate Stabilized Conjugate Gradient Iterative Solver
   *
   *  \param A       coefficient matrix
   *  \param b       righ hand side
   *  \param x       guess and solution
   *  \param P       preconditioner
   *  \param epsi    Admitted tolerance
   *  \param maxIter maximum number of admitted iteration
   *  \param iter    total number of performed itaration
   *  \param pStream pointer to stream object for messages
   *
   *  \return last computed residual
   */
  template<
    typename real_type,
    typename integer,
    typename matrix_type,
    typename vector_type,
    typename preco_type
  >
  real_type
  bicgstab(
    matrix_type const & A,
    vector_type const & b,
    vector_type       & x,
    preco_type  const & P,
    real_type           epsi,
    integer             maxIter,
    integer           & iter,
    ostream_type      * pStream = nullptr
  ) {

    using std::abs;

    UTILS_ASSERT(
      A.nrows() == b.size() &&
      A.ncols() == x.size() &&
      A.nrows() == A.ncols(),
      "Sparse_tool::bicgstab, bad system:\n"
      "dim matrix  = {} x {}\n"
      "dim r.h.s.  = {}\n"
      "dim unknown = {}\n",
      A.nrows(), A.ncols(), b.size(), x.size()
    );

    using std::abs;
    using vType = typename vector_type::real_type;

    integer neq{ static_cast<integer>(b.size()) };
    vector_type p(neq), s(neq), t(neq), v(neq), r(neq), rtilde(neq);

    iter = 1;

    vType rhom1 = 1;
    vType alpha = 1;
    vType omega = 1;

    r = b - A * x;
    if ( pStream != nullptr )
      fmt::print( *pStream,
        "[bicgstab] initial (unpreconditioned) residual = {}\n",
        r.template lpNorm<Eigen::Infinity>()
      );

    r /= P;

    real_type resid{ r.template lpNorm<Eigen::Infinity>() };
    if ( pStream != nullptr )
      fmt::print( *pStream,
        "[bicgstab] initial (preconditioned) residual   = {}\n",
        resid
      );

    if ( resid <= epsi ) goto fine;

    rtilde = r;
    p.setZero();
    v.setZero();

    do {

      vType rho = rtilde.dot(r);
      UTILS_ASSERT0(
        abs(rho) != vType(0),
        "Sparse_tool::bicgstab, breakdown, rho == 0\n"
      );

      vType beta = (rho/rhom1) * (alpha/omega);

      p = r + beta * (p - omega * v);
      v = A * p;
      v /= P;
      alpha = rtilde.dot(v);
      UTILS_ASSERT0(
        alpha != vType(0),
        "Sparse_tool::bicgstab, breakdown, tau == 0\n"
      );
      alpha = rho / alpha;

      s = r - alpha * v;

      t = A * s;
      t /= P;

      omega = t.dot(t);
      UTILS_ASSERT0(
        abs(omega) != vType(0),
        "Sparse_tool::bicgstab, breakdown, omega == 0\n"
      );
      omega = t.dot(s) / omega;

      x = x + alpha * p + omega * s;
      r = s - omega * t;

      resid = r.template lpNorm<Eigen::Infinity>();
      if ( pStream != nullptr )
        fmt::print( *pStream, "[bicgstab] iter = {:<5} residual = {}\n", iter, resid );

      if ( resid <= epsi) goto fine;

      rhom1 = rho;

    } while ( ++iter <= maxIter && resid > epsi );

    if ( pStream != nullptr )
      fmt::print( *pStream, "[bicgstab] iter = {:<5} residual = {}\n", iter, resid );

  fine:

    return resid;
  }
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::bicgstab;
}
#endif

#endif
