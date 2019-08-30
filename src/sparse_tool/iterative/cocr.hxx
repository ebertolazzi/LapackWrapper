#ifndef SPARSETOOL_ITERATIVE_COCR_HH
#define SPARSETOOL_ITERATIVE_COCR_HH

using namespace std;

namespace SparseTool {

  /*
  //                              
  //   ####   ####   ####  #####  
  //  #    # #    # #    # #    # 
  //  #      #    # #      #    # 
  //  #      #    # #      #####  
  //  #    # #    # #    # #   #  
  //   ####   ####   ####  #    #                             
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
  template <typename valueType,
            typename indexType,
            typename matrix_type,
            typename vector_type,
            typename preco_type>
  valueType
  cocr(
    matrix_type const & A,
    vector_type const & b,
    vector_type       & x,
    preco_type  const & P,
    valueType   const & epsi,
    indexType   const   maxIter,
    indexType         & iter,
    ostream           * pStream = nullptr
  ) {

    SPARSETOOL_ASSERT(
      A.numRows() == b.size() &&
      A.numCols() == x.size() &&
      A.numRows() == A.numCols(),
      "Bad system in cocr" <<
      "dim matrix  = " << A.numRows() <<
      " x " << A.numCols() <<
      "\ndim r.h.s.  = " << b.size() <<
      "\ndim unknown = " << x.size()
    )

    typedef typename vector_type::valueType vType;

    indexType   neq = b.size();
    vector_type p(neq), q(neq), qt(neq), r(neq), rt(neq), Art(neq);
    valueType   resid;
    vType       rho, mu, alpha, beta;

    r   = b - A * x;
    p   = r / P;
    rt  = p;
    Art = A * rt;
    rho = rdot(rt,Art);

    iter = 1;
    do {

      resid = normi(rt);
      if ( pStream != nullptr )
        (*pStream) << "iter = " << iter << " residual = " << resid << '\n';

      if ( resid <= epsi ) break;

      q  = A * p;
      qt = q / P;
      mu = rdot(q,qt);
      SPARSETOOL_ASSERT( mu != vType(0), "COCR failed found mu == 0\n" )

      alpha = rho/mu;
      x    = x + alpha * p;
      r    = r - alpha * q;
      rt   = r / P;
      Art  = A * rt;

      beta = rho;
      rho  = rdot(rt,Art);
      beta = rho/beta;
      SPARSETOOL_ASSERT( beta != vType(0), "COCR failed found beta == 0\n" )

      p = rt + beta * p;

    }  while ( ++iter <= maxIter );

    return resid;
  }

}

namespace SparseToolLoad {
  using ::SparseTool::cocr;
}

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
