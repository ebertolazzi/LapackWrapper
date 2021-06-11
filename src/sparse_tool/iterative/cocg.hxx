#pragma once
#ifndef SPARSETOOL_ITERATIVE_COCG_dot_HH
#define SPARSETOOL_ITERATIVE_COCG_dot_HH

namespace SparseTool {

  /*
  //   ####   ####   ####   ####
  //  #    # #    # #    # #    #
  //  #      #    # #      #
  //  #      #    # #      #  ###
  //  #    # #    # #    # #    #
  //   ####   ####   ####   ####
  */
  /*!
   *  Preconditioned Conjugate Gradient Iterative Solver
   *  \param A       coefficient matrix
   *  \param b       righ hand side
   *  \param x       guess and solution
   *  \param P       preconditioner
   *  \param epsi    Admitted tolerance
   *  \param maxIter maximum number of admitted iteration
   *  \param iter    total number of performed itaration
   *  \param pStream pointer to stream object for messages
   *  \return last computed residual
   *
   *  Use preconditioned conjugate gradient (for complex symmetric system)
   *  to solve \f$ A x = b \f$.
   */
  template <typename real_type,
            typename integer,
            typename matrix_type,
            typename vector_type,
            typename preco_type>
  real_type
  cocg(
    matrix_type const & A,
    vector_type const & b,
    vector_type       & x,
    preco_type  const & P,
    real_type   const & epsi,
    integer   const   maxIter,
    integer         & iter,
    ostream           * pStream = nullptr
  ) {
    
    SPARSETOOL_ASSERT(
      A.numRows() == b.size() &&
      A.numCols() == x.size() &&
      A.numRows() == A.numCols(),
      "Bad system in cocg" <<
      "dim matrix  = " << A.numRows() <<
      " x " << A.numCols() <<
      "\ndim r.h.s.  = " << b.size() <<
      "\ndim unknown = " << x.size()
    )

    typedef typename vector_type::real_type vType;

    integer   neq = b.size();
    vector_type p(neq), q(neq), r(neq), rt(neq);
    real_type   resid;
    vType       rho, mu, alpha, beta;

    r   = b - A * x;
    p   = r / P;
    rt  = p;
    rho = rdot(r,rt);

    iter = 1;
    do {
      resid = rt.template lpNorm<Eigen::Infinity>();

      if ( pStream != nullptr )
        (*pStream) << "iter = " << iter << " residual = " << resid << '\n';

      if ( resid <= epsi ) break;

      q  = A * p;
      mu = rdot(q,p);
      SPARSETOOL_ASSERT( mu != vType(0), "COCG failed found mu == 0\n" )

      alpha = rho/mu;
      x  += alpha * p;
      r  -= alpha * q;
      rt  = r / P;

      beta = rho;
      rho  = rdot(r,rt);
      beta = rho/beta;

      p = rt + beta * p;

    }  while ( ++iter <= maxIter );

    return resid;
  }

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace SparseToolLoad {
  using ::SparseTool::cocg;
}
#endif

#endif
