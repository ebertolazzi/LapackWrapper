#ifndef SPARSETOOL_ITERATIVE_BICGSTAB_HH
#define SPARSETOOL_ITERATIVE_BICGSTAB_HH

using namespace std;

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
  template <typename valueType,
            typename indexType,
            typename matrix_type,
            typename vector_type,
            typename preco_type>
  valueType
  bicgstab(
    matrix_type const & A,
    vector_type const & b,
    vector_type       & x,
    preco_type  const & P,
    valueType           epsi,
    indexType           maxIter,
    indexType         & iter,
    ostream           * pStream = nullptr
  ) {

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

    using ::SparseToolFun::absval;
    typedef typename vector_type::valueType vType;

    indexType neq = b.size();
    vector_type p(neq), s(neq), t(neq), v(neq), r(neq), rtilde(neq);

    iter = 1;

    vType rhom1 = 1;
    vType alpha = 1;
    vType omega = 1;

    r  = b - A * x;
    r  = r / P;

    valueType resid = normi(r);
    if ( resid <= epsi ) goto fine;

    rtilde = r;
    p      = vType(0);
    v      = vType(0);

    do {

      vType rho = dot(rtilde,r);
      SPARSETOOL_ASSERT( absval(rho) != vType(0), "BiCGSTAB breakdown, rho == 0" )

      vType beta = (rho/rhom1) * (alpha/omega);
 
      p = r + beta * (p - omega * v); 
      v = A * p;
      v = v / P;
      alpha = dot(rtilde,v);
      SPARSETOOL_ASSERT( alpha != vType(0), "BiCGSTAB breakdown, tau == 0\n" )
      alpha = rho / alpha;

      s = r - alpha * v;

      t  = A * s;
      t  = t / P;

      omega = dot(t,t);
      SPARSETOOL_ASSERT( absval(omega) != vType(0), "BiCGSTAB breakdown, omega == 0\n" )
      omega = dot(t,s) / omega;

      x = x + alpha * p + omega * s;
      r = s - omega * t;

      resid = normi(r);
      if ( pStream != nullptr ) (*pStream) << "iter = " << iter << " residual = " << resid << '\n';

      if ( resid <= epsi) goto fine;

      rhom1 = rho;

    } while ( ++iter <= maxIter && resid > epsi );

    if ( pStream != nullptr ) (*pStream) << "iter = " << iter << " residual = " << resid << '\n';

  fine:

    return resid;
  }
}

namespace SparseToolLoad {
  using ::SparseTool::bicgstab;
}

#endif
