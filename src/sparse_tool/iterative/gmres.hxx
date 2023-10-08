#pragma once
#ifndef SPARSETOOL_ITERATIVE_GMRES_dot_HH
#define SPARSETOOL_ITERATIVE_GMRES_dot_HH

namespace Sparse_tool {

  /*
  //   #####    #     #  #####   #####   ####
  //  #     #   ##   ##  #    #  #      #    #
  //  #         # # # #  #    #  #      #
  //  #         #  #  #  #####   ####    ####
  //  #  ####   #     #  #    #  #           #
  //  #      #  #     #  #     # #      #    #
  //   ######   #     #  #     # #####   ####
  */
  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  template <class T>
  inline
  void
  GeneratePlaneRotation(
    T const & dx,
    T const & dy,
    T       & cs,
    T       & sn
  ) {

    using std::abs;
    using std::sqrt;

    if ( dy == 0.0 ) {
      cs = 1.0;
      sn = 0.0;
    } else if ( abs(dy) > abs(dx) ) {
      T temp = dx / dy;
      sn = 1.0 / sqrt( 1.0 + temp*temp );
      cs = temp * sn;
    } else {
      T temp = dy / dx;
      cs = 1.0 / sqrt( 1.0 + temp*temp );
      sn = temp * cs;
    }
  }

  template <class T>
  inline
  void
  ApplyPlaneRotation(
    T       & dx,
    T       & dy,
    T const & cs,
    T const & sn
  ) {
    T temp = cs * dx + sn * dy;
    dy = -sn * dx + cs * dy;
    dx = temp;
  }
  #endif

  //!
  //! Generalized Minimal Residual Iterative Solver.
  //!
  //! \param A       coefficient matrix
  //! \param b       righ hand side
  //! \param x       guess and solution
  //! \param P       preconditioner
  //! \param epsi    Admitted tolerance
  //! \param m       maximum dimension of Krilov subspace
  //! \param maxIter maximum number of admitted iteration
  //! \param iter    total number of performed itaration
  //! \param pStream pointer to stream object for messages
  //! \return        last computed residual
  //!
  template<
    typename real_type,
    typename matrix_type,
    typename vector_type,
    typename preco_type
  >
  real_type
  gmres(
    matrix_type const & A,
    vector_type const & b,
    vector_type &       x,
    preco_type  const & P,
    real_type           epsi,
    integer             m, // maxSubIter
    integer             maxIter,
    integer &           iter,
    ostream_type *      pStream = nullptr
  ) {

    using std::abs;

    UTILS_ASSERT(
      A.nrows() == b.size() &&
      A.ncols() == x.size() &&
      A.nrows() == A.ncols(),
      "Sparse_tool::gmres, bad system:\n"
      "dim matrix  = {} x {}\n"
      "dim r.h.s.  = {}\n"
      "dim unknown = {}\n",
      A.nrows(), A.ncols(), b.size(), x.size()
    );

    real_type resid = 0;
    integer   m1    = m+1;
    integer   neq   = b.size();
    vector_type w(neq), r(neq), H(m1*m1), s(m1), cs(m1), sn(m1);

    Eigen::Matrix<typename vector_type::real_type,Eigen::Dynamic,Eigen::Dynamic> v(neq,m+1);

    iter = 0;
    do {

      r = b - A * x;
      if ( pStream != nullptr )
        fmt::print( *pStream,
          "[gmres] initial (unpreconditioned) residual [inf norm] = {}\n",
          r.template lpNorm<Eigen::Infinity>()
       );

      r /= P;
      if ( pStream != nullptr )
        fmt::print( *pStream,
          "[gmres] initial (preconditioned) residual [inf norm]   = {}\n",
          r.template lpNorm<Eigen::Infinity>()
        );

      real_type beta = r.norm();
      if ( beta <= epsi ) goto fine;

      typename vector_type::real_type betax = beta;
      v.col(0) = r / betax;
      s(0)     = beta;
      for ( integer k = 1; k <= m; ++k ) s(k) = 0;

      integer i{0};
      do {

        w = A * v.col(i);
        w /= P;

        integer k;
        for ( k = 0; k <= i; ++k ) {
          H(k*m1+i) = w.dot(v.col(k));
          w -= H(k*m1+i) * v.col(k);
        }

        H((i+1)*m1+i) = w.norm();
        v.col(i+1)    = w / H((i+1)*m1+i);

        for ( k = 0; k < i; ++k )
          ApplyPlaneRotation(H(k*m1+i), H((k+1)*m1+i), cs(k), sn(k));

        GeneratePlaneRotation(H(i*m1+i), H((i+1)*m1+i), cs(i), sn(i));
        ApplyPlaneRotation(H(i*m1+i), H((i+1)*m1+i), cs(i), sn(i));
        ApplyPlaneRotation(s(i), s(i+1), cs(i), sn(i));

        ++i; ++iter;
        resid = abs(s(i));
        if ( pStream != nullptr )
          fmt::print( *pStream, "[gmres] iter = {:<5} residual [2 norm] = {}\n", iter, resid );

      } while ( i < m && iter <= maxIter && resid > epsi );

      // Backsolve:
      for ( int ii = i-1; ii >= 0; --ii ) {
        s(ii) /= H(ii*m1+ii);
        for ( int jj = ii - 1; jj >= 0; --jj )
          s(jj) -= H(jj*m1+ii) * s(ii);
      }

      for ( integer jj = 0; jj < i; ++jj )
        x += v.col(jj) * s(jj);

    } while ( iter <= maxIter );

  fine:

    return resid;
  }
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::gmres;
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
