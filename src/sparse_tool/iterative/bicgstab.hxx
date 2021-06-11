#pragma once
#ifndef SPARSETOOL_ITERATIVE_BICGSTAB_dot_HH
#define SPARSETOOL_ITERATIVE_BICGSTAB_dot_HH

namespace SparseTool {

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
  template <typename real_type,
            typename integer,
            typename matrix_type,
            typename vector_type,
            typename preco_type>
  real_type
  bicgstab(
    matrix_type const & A,
    vector_type const & b,
    vector_type       & x,
    preco_type  const & P,
    real_type           epsi,
    integer             maxIter,
    integer           & iter,
    ostream           * pStream = nullptr
  ) {

    using std::abs;

    SPARSETOOL_ASSERT(
      A.numRows() == b.size() &&
      A.numCols() == x.size() &&
      A.numRows() == A.numCols(),
      "Bad system in bicgstab" <<
      "dim matrix  = " << A.numRows() <<
      " x " << A.numCols() <<
      "\ndim r.h.s.  = " << b.size() <<
      "\ndim unknown = " << x.size()
    )

    using std::abs;
    typedef typename vector_type::real_type vType;

    integer neq = b.size();
    vector_type p(neq), s(neq), t(neq), v(neq), r(neq), rtilde(neq);

    iter = 1;

    vType rhom1 = 1;
    vType alpha = 1;
    vType omega = 1;

    r  = b - A * x;
    r  = r / P;

    real_type resid = r.template lpNorm<Eigen::Infinity>();
    if ( resid <= epsi ) goto fine;

    rtilde = r;
    p      = vType(0);
    v      = vType(0);

    do {

      vType rho = rtilde.dot(r);
      SPARSETOOL_ASSERT( abs(rho) != vType(0), "BiCGSTAB breakdown, rho == 0" )

      vType beta = (rho/rhom1) * (alpha/omega);
 
      p = r + beta * (p - omega * v); 
      v = A * p;
      v = v / P;
      alpha = rtilde.dot(v);
      SPARSETOOL_ASSERT( alpha != vType(0), "BiCGSTAB breakdown, tau == 0\n" )
      alpha = rho / alpha;

      s = r - alpha * v;

      t  = A * s;
      t  = t / P;

      omega = t.dot(t);
      SPARSETOOL_ASSERT( abs(omega) != vType(0), "BiCGSTAB breakdown, omega == 0\n" )
      omega = t.dot(s) / omega;

      x = x + alpha * p + omega * s;
      r = s - omega * t;

      resid = r.template lpNorm<Eigen::Infinity>();
      if ( pStream != nullptr ) (*pStream) << "iter = " << iter << " residual = " << resid << '\n';

      if ( resid <= epsi) goto fine;

      rhom1 = rho;

    } while ( ++iter <= maxIter && resid > epsi );

    if ( pStream != nullptr ) (*pStream) << "iter = " << iter << " residual = " << resid << '\n';

  fine:

    return resid;
  }
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace SparseToolLoad {
  using ::SparseTool::bicgstab;
}
#endif

#endif
