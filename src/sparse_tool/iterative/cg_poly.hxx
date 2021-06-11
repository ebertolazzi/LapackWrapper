#pragma once
#ifndef SPARSETOOL_ITERATIVE_CG_POLY_dot_HH
#define SPARSETOOL_ITERATIVE_CG_POLY_dot_HH

namespace SparseTool {

  template <typename real_type, typename integer>
  void
  mulPoly(
    integer                 neq,
    integer                 mdegree,
    Vector<real_type> const & A_A,
    Vector<integer> const & A_R,
    Vector<integer> const & A_J,
    Vector<real_type>       & s0,
    Vector<real_type>       & s1,
    Vector<real_type>       & y,
    Vector<real_type> const & v
  ) {
    // s0 = 1.5*v; s1 = 4*v - 10/3 * A*v
    integer const * _pR = & A_R.front();
    integer const * _pJ = & A_J.front();
    real_type const * _pA = & A_A.front();
    for ( integer k=0; k < neq; ++k, ++_pR ) {
      real_type Av(0);
      for ( integer i_cnt = _pR[1] - _pR[0]; i_cnt > 0; --i_cnt )
        Av += *_pA++ * v(*_pJ++);
      s0(k) = 1.5 * v(k);
      s1(k) = 4.0 * v(k) - (10./3.) * Av;
    };
    real_type g = -20./9.;
    for ( integer n = 2; n <= mdegree; ++n ) {
      real_type a = 0.5+0.125/((n+0.5)*(n+1.5));
      real_type b = g*(n*(n+1.0))/(16*(n+0.5)*(n+0.5));
      g = -1.0/(a+b);
      
      // y = g1*((A-a)*s1-b*s0-1);
      _pR = & A_R.front();
      _pJ = & A_J.front();
      _pA = & A_A.front();
      for ( integer k=0; k < neq; ++k, ++_pR ) {
        real_type As1(0);
        for ( integer i_cnt = _pR[1] - _pR[0]; i_cnt > 0; --i_cnt )
          As1 += *_pA++ * s1(*_pJ++);
        y(k) = g*(As1-a*s1(k)-b*s0(k)-v(k));
      };
      s0 = s1;
      s1 = y;
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
  template <typename real_type,
            typename integer,
            typename matrix_type,
            typename vector_type>
  real_type
  cg_poly(
    matrix_type const & A,
    vector_type const & b,
    vector_type       & x,
    real_type   const & epsi,
    integer   const   maxIter,
    integer   const   mdegree,
    integer         & iter,
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
    typedef typename vector_type::real_type vType;

    integer neq = b.size();

    // memorizzo matrice scalata
    Vector<integer> A_R(neq+1), A_J, Annz(neq);
    Vector<real_type> A_A, dScale(neq), s0(neq), s1(neq);

    // vettori per CG
    Vector<real_type> p(neq), q(neq), r(neq), Ap(neq);
    real_type         rho, rho_1, resid;

    // step 0: compute necessary memory
    Annz   . setZero();
    dScale . setZero();
    for ( A.Begin(); A.End(); A.Next() ) {
      integer i = A.row();
      //integer j = A.column();
      dScale(i) += abs(A.value());
      ++Annz(i);
    }
    
    dScale = sqrt(dScale);

    // step 1: initialize structure
    A_R(0) = 0;
    for ( integer i = 0; i < neq; ++i ) A_R(i+1) = A_R(i) + Annz(i);
    
    // step 2: allocate memory
    A_A . resize( A_R(neq) );
    A_J . resize( A_R(neq) );
    A_A . setZero();
    
    // step 3: fill structure
    for ( A.Begin(); A.End(); A.Next() ) {
      integer i = A.row();
      integer j = A.column();
      A_J(A_R(i)+(--Annz(i))) = j;
    }
    
    // step 4: sort structure
    for ( integer i = 0; i < neq; ++i ) std::sort( &A_J(A_R(i)), &A_J(A_R(i+1)) );
    
    // insert values
    for ( A.Begin(); A.End(); A.Next() ) {
      integer i   = A.row();
      integer j   = A.column();
      integer lo  = A_R(i);
      integer hi  = A_R(i+1);
      integer len = hi - lo;
      while ( len > 0 ) {
        integer half = len / 2;
        integer mid  = lo + half;
        if ( A_J(mid) < j ) { lo = mid + 1; len -= half + 1; }
        else                  len = half;
      }
      A_A(lo) = A.value()/(dScale(i)*dScale(j));
    }

    // r  = b/dScale - A * x;
    integer const * _pR = & A_R.front();
    integer const * _pJ = & A_J.front();
    real_type const * _pA = & A_A.front();
    for ( integer k=0; k < neq; ++k, ++_pR ) {
      real_type tmp(0);
      for ( integer i_cnt = _pR[1] - _pR[0]; i_cnt > 0; --i_cnt )
        { tmp += *_pA++ * x(*_pJ) * dScale(*_pJ); ++_pJ; }
      r(k) = b(k)/dScale(k) - tmp;
    };

    // applico precondizionatore
    SparseTool::mulPoly<real_type,integer>( neq, mdegree, A_A, A_R, A_J, s0, s1, p, r );
    //p = r;
    //p   = r / P;
    rho = rdot(p,r);

    iter = 1;
    do {

      resid = r.template lpNorm<Eigen::Infinity>();
      if ( pStream != nullptr )
        (*pStream) << "iter = " << iter << " residual = " << resid << '\n';

      if ( resid <= epsi ) break;

      // prodotto matrice vettore
      // Ap  = A * p;
      _pR = & A_R.front();
      _pJ = & A_J.front();
      _pA = & A_A.front();
      for ( integer k=0; k < neq; ++k, ++_pR ) {
        real_type tmp(0);
        for ( integer i_cnt = _pR[1] - _pR[0]; i_cnt > 0; --i_cnt )
          tmp += *_pA++ * p(*_pJ++);
        Ap(k) = tmp;
      };

      vType alpha = rho/dot(Ap,p);
      x = x + alpha * p;
      r = r - alpha * Ap;

      // applico precondizionatore
      SparseTool::mulPoly<real_type,integer>( neq, mdegree, A_A, A_R, A_J, s0, s1, q, r );
      //q = r;
      //q = r / P;

      rho_1 = rho;
      rho   = dot(q,r);

      p = q + (rho / rho_1) * p;

    }  while ( ++iter <= maxIter );
    
    x /= dScale;

    return resid;
  }

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace SparseToolLoad {
  using ::SparseTool::cg;
}
#endif

#endif
