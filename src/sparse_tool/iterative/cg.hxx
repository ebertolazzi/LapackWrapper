#ifndef SPARSETOOL_ITERATIVE_CG_HH
#define SPARSETOOL_ITERATIVE_CG_HH

using namespace std;

namespace SparseTool {

  /*
  //  #####   #####
  // #     # #     #
  // #       #
  // #       #  ####
  // #       #     #
  // #     # #     #
  //  #####   #####
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
   *  Use preconditioned conjugate gradient to solve \f$ A x = b \f$.
   */
  template <typename valueType,
            typename indexType,
            typename matrix_type,
            typename vector_type,
            typename preco_type>
  valueType
  cg(
    matrix_type const & A,
    vector_type const & b,
    vector_type       & x,
    preco_type  const & P,
    valueType   const & epsi,
    indexType           maxIter,
    indexType         & iter,
    ostream           * pStream = nullptr
  ) {

    SPARSETOOL_ASSERT(
      A.numRows() == b.size() &&
      A.numCols() == x.size() &&
      A.numRows() == A.numCols(),
      "Bad system in cg" <<
      "dim matrix  = " << A.numRows() <<
      " x " << A.numCols() <<
      "\ndim r.h.s.  = " << b.size() <<
      "\ndim unknown = " << x.size()
    )
    typedef typename vector_type::valueType vType;

    indexType   neq = b.size();
    vector_type p(neq), q(neq), r(neq), Ap(neq);
    valueType   resid;
    vType       rho, rho_1;

    r   = b - A * x;
    p   = r / P;
    rho = rdot(p,r);

    iter = 1;
    do {

      resid = normi(r);
      if ( pStream != nullptr )
        (*pStream) << "iter = " << iter << " residual = " << resid << '\n';

      if ( resid <= epsi ) break;

      Ap  = A * p;
      vType alpha = rho/dot(Ap,p);
      x = x + alpha * p;
      r = r - alpha * Ap;
      q = r / P;

      rho_1 = rho;
      rho   = dot(q,r);

      p = q + (rho / rho_1) * p;

    }  while ( ++iter <= maxIter );

    return resid;
  }

}

namespace SparseToolLoad {
  using ::SparseTool::cg;
}

#endif
