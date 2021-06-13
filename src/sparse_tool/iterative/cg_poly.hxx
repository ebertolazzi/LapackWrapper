#pragma once
#ifndef SPARSETOOL_ITERATIVE_CG_POLY_dot_HH
#define SPARSETOOL_ITERATIVE_CG_POLY_dot_HH

namespace Sparse_tool {

  template <typename real_type, typename integer>
  void
  mulPoly(
    integer                       mdegree,
    CRowMatrix<real_type> const & A,
    Vector<real_type>           & s0,
    Vector<real_type>           & s1,
    Vector<real_type>           & y,
    Vector<real_type>     const & v
  ) {
    real_type a0 = 3.0/2.0;
    real_type a1 = -10.0/3.0;
    real_type b1 = 4;
    // s0 = 1.5*v; s1 = 4*v - 10/3 * A*v
    s1  = a0*v;
    y   = b1*v;
    y  += a1*A*v;
    for ( integer i_n = 2; i_n <= mdegree; ++i_n ) {
      s0 = s1;
      s1 = y;
      real_type n     = i_n;
      real_type n2sq  = (n+2)*(n+2);
      real_type Delta = 2*((3*n+6)*n+2)/((2*n+1)*n2sq);
      real_type c     = Delta-1;
      real_type b     = 2-Delta;
      real_type a     = 2*(3*n+5)/n2sq-4;
      // y = a*(A*s1-v)+b*s1+c*s0;
      y  = b*s1+c*s0-a*v;
      y += a*A*s1;
    }
  }

  /*
  //  ####   ####          #####   ####  #      #   # 
  // #    # #    #         #    # #    # #       # #  
  // #      #              #    # #    # #        #   
  // #      #  ###         #####  #    # #        #   
  // #    # #    #         #      #    # #        #   
  //  ####   ####          #       ####  ######   #   
  //               #######                            
  */
  //!
  //! Preconditioned Conjugate Gradient Iterative Solver.
  //!
  //! \param A       coefficient matrix
  //! \param b       righ hand side
  //! \param x       guess and solution
  //! \param epsi    Admitted tolerance
  //! \param maxIter maximum number of admitted iteration
  //! \param mdegree degree of precondition polynomial
  //! \param iter    total number of performed itaration
  //! \param pStream pointer to stream object for messages
  //! \return last computed residual
  //!
  //! Use preconditioned conjugate gradient to solve \f$ A x = b \f$.
  //!
  template <
    typename real_type,
    typename matrix_type,
    typename vector_type
  >
  real_type
  cg_poly(
    matrix_type const & A,
    vector_type const & b,
    vector_type       & x,
    real_type   const & epsi,
    integer             maxIter,
    integer             mdegree,
    integer           & iter,
    ostream_type      * pStream = nullptr
  ) {

    UTILS_ASSERT(
      A.numRows() == b.size() &&
      A.numCols() == x.size() &&
      A.numRows() == A.numCols(),
      "Sparse_tool::cg_poly, bad system:\n"
      "dim matrix  = {} x {}\n"
      "dim r.h.s.  = {}\n"
      "dim unknown = {}\n",
      A.numRows(), A.numCols(), b.size(), x.size()
    );

    integer neq = b.size();

    // memorizzo matrice scalata
    CRowMatrix<real_type> Ascaled(A);
    Vector<real_type>     dScale(neq), s0(neq), s1(neq);

    // vettori per CG
    Vector<real_type> p(neq), q(neq), r(neq), Ap(neq);
    real_type         rho, rho_1, resid;

    // step 0: compute necessary memory
    dScale.setZero();
    for ( Ascaled.Begin(); Ascaled.End(); Ascaled.Next() ) {
      integer i = Ascaled.row();
      //integer j = A.column();
      dScale(i) += abs(Ascaled.value());
    }

    dScale = sqrt(1.01*dScale.array());

    // scale matrix values
    for ( Ascaled.Begin(); Ascaled.End(); Ascaled.Next() ) {
      integer i = Ascaled.row();
      integer j = Ascaled.column();
      Ascaled.value() /= dScale(i)*dScale(j);
    }

    // A x = b
    // D^(-1) * A * x = D^(-1) * b
    // (D^(-1) * A * D^(-1)) * D * x = D^(-1) * b

    // r  = b/dScale - A * x;
    x.array() *= dScale.array();
    r = b.array()/dScale.array();
    r -= Ascaled * x;

    // applico precondizionatore p = r / P
    Sparse_tool::mulPoly<real_type,integer>( mdegree, Ascaled, s0, s1, p, r );
    rho = p.dot(r);

    iter = 1;
    do {

      resid = r.template lpNorm<Eigen::Infinity>();
      if ( pStream != nullptr )
        fmt::print( *pStream,
          "[cg_poly] iter = {:4} residual = {:.6}\n", iter, resid
        );

      if ( resid <= epsi ) break;

      // prodotto matrice vettore
      Ap = Ascaled * p;

      real_type alpha = rho/Ap.dot(p);
      x += alpha * p;
      r -= alpha * Ap;

      // applico precondizionatore q = r / P
      Sparse_tool::mulPoly<real_type,integer>( mdegree, Ascaled, s0, s1, q, r );

      rho_1 = rho;
      rho   = q.dot(r);

      p = q + (rho / rho_1) * p;

    }  while ( ++iter <= maxIter );

    x.array() /= dScale.array();

    return resid;
  }

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::cg;
}
#endif

#endif
