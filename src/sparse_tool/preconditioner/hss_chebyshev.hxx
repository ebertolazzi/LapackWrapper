#pragma once
#ifndef SPARSETOOL_ITERATIVE_PRECO_HSS_CHEBYSHEV_HH
#define SPARSETOOL_ITERATIVE_PRECO_HSS_CHEBYSHEV_HH

namespace Sparse_tool {

  /*
  //    ____ _   _ _____ ______   ______  _   _ _______     __
  //   / ___| | | | ____| __ ) \ / / ___|| | | | ____\ \   / /
  //  | |   | |_| |  _| |  _ \\ V /\___ \| |_| |  _|  \ \ / /
  //  | |___|  _  | |___| |_) || |  ___) |  _  | |___  \ V /
  //   \____|_| |_|_____|____/ |_| |____/|_| |_|_____|  \_/
  //
  */
  //! Incomplete \c LDU preconditioner
  template <typename T>
  class HSS_CHEBYSHEV_Preconditioner : public Preco<HSS_CHEBYSHEV_Preconditioner<T> > {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    using HSS_CHEBYSHEV_PRECO = HSS_CHEBYSHEV_Preconditioner<T>;
    using PRECO               = Preco<HSS_CHEBYSHEV_PRECO>;
    #endif

    using real_type  = T; //!< type of the elements of the preconditioner
    using rreal_type = typename T::value_type; //!< type of the elements of the preconditioner

  private:

    integer          neq, mdegree;
    rreal_type         sqrt_epsilon;

    Vector<integer>  A_R;
    Vector<integer>  A_J;
    Vector<rreal_type> A_A;
    Vector<integer>  Annz;

    mutable Vector<real_type> s0, s1, y;

    //! build incomplete LDU decomposition with specified pattern `P`
    template <typename MAT>
    void
    build_HSS_CHEBYSHEV( MAT const & A, integer m, rreal_type delta ) {

      UTILS_ASSERT0(
        A.is_ordered(),
        "Sparse_tool: HSS_CHEBYSHEV_Preconditioner::build_LDU\n"
        "pattern must be ordered before use\n"
      );
      UTILS_ASSERT0(
        A.nrows() == A.ncols(),
        "Sparse_tool: HSS_CHEBYSHEV_Preconditioner::build_LDU\n"
        "only square matrix allowed\n"
      );
      UTILS_ASSERT0(
        A.nrows() > 0,
        "Sparse_tool: HSS_CHEBYSHEV_Preconditioner::build_LDU\n"
        "empty matrix\n"
      );

      mdegree = m;
      sqrt_epsilon = std::pow( (1+std::sqrt(1-delta*delta))/delta,1.0/(m+1.0) );
      sqrt_epsilon = (sqrt_epsilon-1)/(sqrt_epsilon+1);

      UTILS_ASSERT0(
        sqrt_epsilon > 0 && sqrt_epsilon < 1,
        "Sparse_tool: HSS_CHEBYSHEV_Preconditioner::build_LDU\n"
        "computed epsilon must be in the interval (0,1)"
      );

      neq = A.nrows();
      s0.resize(neq);
      s1.resize(neq);
      y.resize(neq);

      // step 0: compute necessary memory
      PRECO::pr_size = A.nrows();
      Annz.resize( PRECO::pr_size );
      Annz.setZero();

      for ( A.Begin(); A.End(); A.Next() ) {
        integer i{A.row()};
        //integer j{A.column()};
        ++Annz(i);
      }

      // step 1: initialize structure
      A_R.resize( PRECO::pr_size + 1 );
      A_R(0) = 0;
      for ( integer i{0}; i < PRECO::pr_size; ++i ) A_R(i+1) = A_R(i) + Annz(i);

      // step 2: allocate memory
      A_A.resize( A_R(PRECO::pr_size) );
      A_J.resize( A_R(PRECO::pr_size) );
      A_A.setZero();

      // step 3: fill structure
      for ( A.Begin(); A.End(); A.Next() ) {
        integer i{A.row()};
        integer j{A.column()};
        A_J(A_R(i)+(--Annz(i))) = j;
      }

      // step 4: sort structure
      for ( integer i{0}; i < PRECO::pr_size; ++i ) std::sort( &A_J(A_R(i)), &A_J(A_R(i+1)) );

      // insert values
      for ( A.Begin(); A.End(); A.Next() ) {
        integer i   { A.row()    };
        integer j   { A.column() };
        integer lo  { A_R(i)     };
        integer hi  { A_R(i+1)   };
        integer len { hi - lo    };
        while ( len > 0 ) {
          integer half { len / 2 };
          integer mid  { lo + half };
          if ( A_J(mid) < j ) { lo = mid + 1; len -= half + 1; }
          else                  len = half;
        }
        A_A(lo) = A.value().real() + A.value().imag();
      }
    }

    //! apply preconditioner to vector `v`  and store result in vector \c y
    void
    multiply_by_poly( Vector<real_type> & _y, Vector<real_type> const & v ) const {
      rreal_type e { sqrt_epsilon*sqrt_epsilon };
      rreal_type a { 2.0/(1.0-e) };
      rreal_type b { -(1.0+e)/(1.0-e) };
      rreal_type c { (sqrt_epsilon-1)/(sqrt_epsilon+1) };
      // s0 = 1.5*v; s1 = 4*v - 10/3 * A*v
      integer  const *   pR { A_R.data() };
      integer  const *   pJ { A_J.data() };
      rreal_type const * pA { A_A.data() };
      for ( integer k{0}; k < PRECO::pr_size; ++k, ++pR ) {
        real_type Av{0};
        for ( integer i_cnt{pR[1] - pR[0]}; i_cnt > 0; --i_cnt )
          Av += *pA++ * v(*pJ++);
        s1(k) = (2.0/(1.0+e)) * v(k);
        _y(k) = (8.0*(1.0+e)/(e*(e+6.0)+1.0)) * v(k) - (8.0/(e*(e+6.0)+1.0)) * Av;
      };
      rreal_type cm1 { 1 };
      rreal_type cn  { c };
      rreal_type cp1 { c*c };
      for ( integer n{2}; n <= mdegree; ++n ) {
        cm1  = cn;
        cn   = cp1;
        cp1 *= c;
        s0   = s1;
        s1   = _y;
        rreal_type g0 { (cn+1/cn)/(cp1+1/cp1) };
        rreal_type aa { 2*a*g0 };
        rreal_type bb { 2*b*g0 };
        rreal_type cc { -(cm1+1/cm1)/(cp1+1/cp1) };
        pR = A_R.data();
        pJ = A_J.data();
        pA = A_A.data();
        for ( integer k{0}; k < PRECO::pr_size; ++k, ++pR ) {
          real_type As1(0,0);
          for ( integer i_cnt{pR[1] - pR[0]}; i_cnt > 0; --i_cnt )
            As1  += *pA++ * s1(*pJ++);
          _y(k) = aa*(As1-v(k))+bb*s1(k)+cc*s0(k);
        };
      }
    }

  public:

    HSS_CHEBYSHEV_Preconditioner(void) : Preco<HSS_CHEBYSHEV_PRECO>() {}

    template <typename MAT>
    HSS_CHEBYSHEV_Preconditioner(
      MAT const & M,
      integer   m,
      rreal_type  delta
    ) : Preco<HSS_CHEBYSHEV_PRECO>()
    { build_HSS_CHEBYSHEV( M, m, delta ); }

    //! build the preconditioner from matrix `M`.
    template <typename MAT>
    void
    build( MAT const & M, integer m, rreal_type delta )
    { build_HSS_CHEBYSHEV(M,m,delta); }

    //! apply preconditioner to vector `v`  and store result to vector \c y
    template <typename VECTOR>
    void
    ass_preco( VECTOR & _y, VECTOR const & v ) const {
      multiply_by_poly( _y, v );
      _y *= real_type(0.5,-0.5);
    }

  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,HSS_CHEBYSHEV_Preconditioner<TP> >
  operator / (Vector<T> const & v, HSS_CHEBYSHEV_Preconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,HSS_CHEBYSHEV_Preconditioner<TP> >(v,P); }
  #endif

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::HSS_CHEBYSHEV_Preconditioner;
}
#endif

#endif
