#ifndef SPARSETOOL_ITERATIVE_GMRES_HH
#define SPARSETOOL_ITERATIVE_GMRES_HH

using namespace std;

namespace SparseTool {

  /*
  //   #####    #     #  #####   #####   ####
  //  #     #   ##   ##  #    #  #      #    #
  //  #         # # # #  #    #  #      #
  //  #         #  #  #  #####   ####    ####
  //  #  ####   #     #  #    #  #           #
  //  #      #  #     #  #     # #      #    #
  //   ######   #     #  #     # #####   ####
  */
  //! \cond NODOC
 
  template <class T>
  inline
  void
  GeneratePlaneRotation(
    T const & dx,
    T const & dy,
    T       & cs,
    T       & sn
  ) {
    
    using ::SparseToolFun::absval;
    using ::SparseToolFun::sqrt;
    
    if ( dy == 0.0 ) {
      cs = 1.0;
      sn = 0.0;
    } else if ( absval(dy) > absval(dx) ) {
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
  //! \endcond
  
  /*!
   *  Generalized Minimal Residual Iterative Solver
   *  \param A       coefficient matrix
   *  \param b       righ hand side
   *  \param x       guess and solution
   *  \param P       preconditioner
   *  \param epsi    Admitted tolerance
   *  \param m       maximum dimension of Krilov subspace
   *  \param maxIter maximum number of admitted iteration
   *  \param iter    total number of performed itaration
   *  \param pStream pointer to stream object for messages
   *  \return last computed residual
   */
  template <typename valueType,
            typename indexType,
            typename matrix_type,
            typename vector_type,
            typename preco_type>
  valueType
  gmres(
    matrix_type const & A,
    vector_type const & b,
    vector_type       & x,
    preco_type  const & P,
    valueType           epsi,
    indexType           m, // maxSubIter
    indexType           maxIter,
    indexType         & iter,
    ostream           * pStream = nullptr
  ) {

    using ::SparseToolFun::absval;

    SPARSETOOL_ASSERT(
      A.numRows() == b.size() &&
      A.numCols() == x.size() &&
      A.numRows() == A.numCols(),
      "Bad system in gmres" <<
      "dim matrix  = " << A.numRows() <<
      " x " << A.numCols() <<
      "\ndim r.h.s.  = " << b.size() <<
      "\ndim unknown = " << x.size()
    )
    valueType resid = 0;
    indexType m1    = m+1;
    indexType neq   = b.size();
    vector_type w(neq), r(neq), H(m1*m1), s(m1), cs(m1), sn(m1);
    
    Vector<vector_type> v(m1);
    for ( indexType nv = 0; nv <= m; ++nv ) v(nv) . resize(neq);

    iter = 1;
    do {

      r  = b - A * x;
      r  = r / P;
      valueType beta = norm2(r);
      if ( beta <= epsi ) goto fine;
 
      typename vector_type::valueType betax = beta;
      v(0) = r / betax;
      s(0) = beta;
      for ( indexType k = 1; k <= m; ++k ) s(k) = 0;

      indexType i = 0;
      do {

        w = A * v(i);
        w = w / P;

        indexType k;
        for ( k = 0; k <= i; ++k ) {
          H(k*m1+i) = dot(w, v(k));
          w -= H(k*m1+i) * v(k);
        }

        H((i+1)*m1+i) = norm2(w);
        v(i+1)        = w / H((i+1)*m1+i);

        for ( k = 0; k < i; ++k )
          ApplyPlaneRotation(H(k*m1+i), H((k+1)*m1+i), cs(k), sn(k));

        GeneratePlaneRotation(H(i*m1+i), H((i+1)*m1+i), cs(i), sn(i));
        ApplyPlaneRotation(H(i*m1+i), H((i+1)*m1+i), cs(i), sn(i));
        ApplyPlaneRotation(s(i), s(i+1), cs(i), sn(i));
      
        ++i; ++iter;
        resid = absval(s(i));
        if ( pStream != nullptr ) (*pStream) << "iter = " << iter << " residual = " << resid << '\n';

      } while ( i < m && iter <= maxIter && resid > epsi );

      // Backsolve:  
      for ( int ii = i-1; ii >= 0; --ii ) {
        s(ii) /= H(ii*m1+ii);
        for ( int jj = ii - 1; jj >= 0; --jj )
          s(jj) -= H(jj*m1+ii) * s(ii);
      }

      for ( indexType jj = 0; jj < i; ++jj )
        x += v(jj) * s(jj);

    } while ( iter <= maxIter );

  fine:

    return resid;
  }
}

namespace SparseToolLoad {
  using ::SparseTool::gmres;
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
