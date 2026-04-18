#include <lapack_wrapper/lapack_wrapper.hh>
#include <lapack_wrapper/lapack_wrapper++.hh>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#endif

using real_type = lapack_wrapper::doublereal;
using integer   = lapack_wrapper::integer;
using lapack_wrapper::Transposition;

static Utils::Console msg(&std::cout);

static constexpr real_type BIG_EPS{5e-8};

static
bool
close_enough(
  real_type a,
  real_type b,
  real_type rel = BIG_EPS,
  real_type abs = BIG_EPS
) {
  real_type scale{ std::max(std::abs(a), std::abs(b)) };
  return std::abs(a - b) <= abs + rel * scale;
}

static
void
expect_close(
  char const * what,
  real_type    a,
  real_type    b,
  real_type    tol = BIG_EPS
) {
  UTILS_ASSERT(
    close_enough(a, b, tol, tol),
    "{} mismatch: {} vs {}\n",
    what, a, b
  );
}

static
void
expect_vector_close(
  char const *      what,
  integer           n,
  real_type const * a,
  real_type const * b,
  real_type         tol = BIG_EPS
) {
  for ( integer i{0}; i < n; ++i ) expect_close( what, a[i], b[i], tol );
}

static
real_type
dense_norm1(
  integer           m,
  integer           n,
  real_type const * A,
  integer           ldA
) {
  real_type norm{0};
  for ( integer j{0}; j < n; ++j ) {
    real_type sum{0};
    for ( integer i{0}; i < m; ++i ) sum += std::abs(A[i+j*ldA]);
    if ( sum > norm ) norm = sum;
  }
  return norm;
}

static
real_type
dense_norm_inf(
  integer           m,
  integer           n,
  real_type const * A,
  integer           ldA
) {
  real_type norm{0};
  for ( integer i{0}; i < m; ++i ) {
    real_type sum{0};
    for ( integer j{0}; j < n; ++j ) sum += std::abs(A[i+j*ldA]);
    if ( sum > norm ) norm = sum;
  }
  return norm;
}

static
void
build_transposed_rhs(
  integer           m,
  integer           n,
  real_type const * A,
  integer           ldA,
  real_type const * x,
  real_type *       rhs
) {
  lapack_wrapper::gemv( Transposition::YES, m, n, 1.0, A, ldA, x, 1, 0.0, rhs, 1 );
}

static
void
build_rhs(
  integer           m,
  integer           n,
  real_type const * A,
  integer           ldA,
  real_type const * x,
  real_type *       rhs
) {
  lapack_wrapper::gemv( Transposition::NO, m, n, 1.0, A, ldA, x, 1, 0.0, rhs, 1 );
}

static
void
test_lu_and_lupq_public_api() {
  msg.green("test_lu_and_lupq_public_api\n");

  constexpr integer N{3};
  real_type A[]{
    4.0, 1.0, 2.0,
    1.0, 3.0, 0.0,
    2.0, 0.0, 5.0
  };
  real_type x0[]{ 1.0, -2.0, 0.5 };
  real_type x1[]{ -1.5, 0.25, 2.0 };
  real_type b0[N], b1[N], bt0[N], bt1[N];

  build_rhs( N, N, A, N, x0, b0 );
  build_rhs( N, N, A, N, x1, b1 );
  build_transposed_rhs( N, N, A, N, x0, bt0 );
  build_transposed_rhs( N, N, A, N, x1, bt1 );

  {
    lapack_wrapper::LU<real_type> lu;
    lapack_wrapper::MatrixWrapper<real_type> Amat(A, N, N, N);
    UTILS_ASSERT( lu.factorize( Amat ), "LU::factorize(MatrixWrapper) failed\n" );

    real_type rhs[N];
    lapack_wrapper::copy( N, b0, 1, rhs, 1 );
    UTILS_ASSERT( lu.solve( rhs ), "LU::solve failed\n" );
    expect_vector_close( "LU single solve", N, rhs, x0 );

    lapack_wrapper::copy( N, bt0, 1, rhs, 1 );
    UTILS_ASSERT( lu.t_solve( rhs ), "LU::t_solve failed\n" );
    expect_vector_close( "LU single t_solve", N, rhs, x0 );

    lapack_wrapper::Matrix<real_type> RHS(N,2);
    for ( integer i{0}; i < N; ++i ) {
      RHS(i,0) = b0[i];
      RHS(i,1) = b1[i];
    }
    UTILS_ASSERT( lu.solve( RHS ), "LU::solve(MatrixWrapper) failed\n" );
    expect_vector_close( "LU multi solve rhs0", N, RHS.data(),     x0 );
    expect_vector_close( "LU multi solve rhs1", N, RHS.data() + N, x1 );

    for ( integer i{0}; i < N; ++i ) {
      RHS(i,0) = bt0[i];
      RHS(i,1) = bt1[i];
    }
    UTILS_ASSERT( lu.t_solve( RHS ), "LU::t_solve(MatrixWrapper) failed\n" );
    expect_vector_close( "LU multi t_solve rhs0", N, RHS.data(),     x0 );
    expect_vector_close( "LU multi t_solve rhs1", N, RHS.data() + N, x1 );

    real_type cond1{ lu.cond1( dense_norm1(N, N, A, N) ) };
    real_type condI{ lu.condInf( dense_norm_inf(N, N, A, N) ) };
    UTILS_ASSERT( std::isfinite(cond1) && cond1 > 0, "LU::cond1 invalid\n" );
    UTILS_ASSERT( std::isfinite(condI) && condI > 0, "LU::condInf invalid\n" );
  }

  {
    lapack_wrapper::LU_no_alloc<real_type> lu;
    std::vector<real_type> fact( size_t(N*N) );
    std::vector<real_type> work( size_t(2*N) );
    std::vector<integer>   piv{ size_t(N) };
    std::vector<integer>   iw{ size_t(N) };

    lu.no_allocate( N, N, fact.data(), piv.data(), work.data(), iw.data() );
    UTILS_ASSERT( lu.factorize_nodim( A, N ), "LU_no_alloc::factorize_nodim failed\n" );

    real_type rhs[N];
    lapack_wrapper::copy( N, b1, 1, rhs, 1 );
    UTILS_ASSERT( lu.solve( rhs ), "LU_no_alloc::solve failed\n" );
    expect_vector_close( "LU_no_alloc solve", N, rhs, x1 );
  }

  {
    lapack_wrapper::LUPQ<real_type> lupq;
    UTILS_ASSERT( lupq.factorize( N, N, A, N ), "LUPQ::factorize failed\n" );

    real_type rhs[N];
    lapack_wrapper::copy( N, b0, 1, rhs, 1 );
    UTILS_ASSERT( lupq.solve( rhs ), "LUPQ::solve failed\n" );
    expect_vector_close( "LUPQ solve", N, rhs, x0, 1e-7 );

    lapack_wrapper::copy( N, bt0, 1, rhs, 1 );
    UTILS_ASSERT( lupq.t_solve( rhs ), "LUPQ::t_solve failed\n" );
    expect_vector_close( "LUPQ t_solve", N, rhs, x0, 1e-7 );
  }

  {
    lapack_wrapper::LUPQ_no_alloc<real_type> lupq;
    std::vector<real_type> work( size_t(N*N) );
    std::vector<integer>   iw( size_t(2*N) );

    lupq.no_allocate( N, N*N, work.data(), 2*N, iw.data() );
    UTILS_ASSERT( lupq.factorize_nodim( A, N ), "LUPQ_no_alloc::factorize_nodim failed\n" );

    real_type rhs[N];
    lapack_wrapper::copy( N, b1, 1, rhs, 1 );
    UTILS_ASSERT( lupq.solve( rhs ), "LUPQ_no_alloc::solve failed\n" );
    expect_vector_close( "LUPQ_no_alloc solve", N, rhs, x1, 1e-7 );
  }
}

static
void
test_lsc_public_api() {
  msg.green("test_lsc_public_api\n");

  constexpr integer NR{3};
  constexpr integer NC{2};

  real_type M[]{
    1.0, 1.0, 0.0,
    1.0, 0.0, 1.0
  };
  bool select[]{ true, false, false };

  lapack_wrapper::LSC<real_type> lsc;
  UTILS_ASSERT( lsc.factorize( NR, NC, M, NR, select ), "LSC::factorize failed\n" );

  real_type rhs0[]{ 1.5, 2.0, -1.0 };
  real_type sol0[]{ 2.0, -1.0 };
  real_type rhs1[]{ -3.0, -4.0, 0.5 };
  real_type sol1[]{ -4.0, 0.5 };

  {
    real_type xb[NR];
    lapack_wrapper::copy( NR, rhs0, 1, xb, 1 );
    UTILS_ASSERT( lsc.solve( xb ), "LSC::solve failed\n" );
    expect_vector_close( "LSC single solve", NC, xb, sol0 );
  }

  {
    lapack_wrapper::Matrix<real_type> XB(NR,2);
    for ( integer i{0}; i < NR; ++i ) {
      XB(i,0) = rhs0[i];
      XB(i,1) = rhs1[i];
    }
    UTILS_ASSERT( lsc.solve( XB ), "LSC::solve(MatrixWrapper) failed\n" );
    expect_vector_close( "LSC multi solve rhs0", NC, XB.data(),      sol0 );
    expect_vector_close( "LSC multi solve rhs1", NC, XB.data() + NR, sol1 );
  }

  {
    real_type I[]{
      1.0, 0.0,
      0.0, 1.0
    };
    bool constraints_only[]{ false, false };
    lapack_wrapper::LSC<real_type> only_constraints;
    UTILS_ASSERT(
      only_constraints.factorize( 2, 2, I, 2, constraints_only ),
      "LSC::factorize failed for constraints-only system\n"
    );

    real_type xb[]{ 3.0, -2.0 };
    real_type expected[]{ 3.0, -2.0 };
    UTILS_ASSERT( only_constraints.solve( xb ), "LSC::solve failed on constraints-only system\n" );
    expect_vector_close( "LSC constraints-only", 2, xb, expected );
  }

  {
    bool threw{false};
    real_type dummy[]{ 1.0, 2.0 };
    try {
      lsc.t_solve( dummy );
    } catch ( std::exception const & ) {
      threw = true;
    }
    UTILS_ASSERT( threw, "LSC::t_solve should report unsupported operation\n" );
  }
}

static
void
test_tridiagonal_public_api() {
  msg.green("test_tridiagonal_public_api\n");

  constexpr integer N{5};
  real_type L[]{ -1.0, 0.5, 1.0, -0.25 };
  real_type D[]{  4.0, 5.0, 6.0,  7.0, 8.0 };
  real_type U[]{  0.75, -1.25, 0.5, 1.5 };
  real_type A[]{
    4.0, -1.0, 0.0,  0.0,   0.0,
    0.75, 5.0, 0.5,  0.0,   0.0,
    0.0, -1.25, 6.0, 1.0,   0.0,
    0.0, 0.0,   0.5, 7.0,  -0.25,
    0.0, 0.0,   0.0, 1.5,   8.0
  };
  real_type x0[]{ 1.0, -2.0, 0.5, 3.0, -1.0 };
  real_type x1[]{ -1.5, 0.25, 2.0, -0.5, 1.0 };
  real_type b0[N], b1[N], bt0[N], bt1[N];

  build_rhs( N, N, A, N, x0, b0 );
  build_rhs( N, N, A, N, x1, b1 );
  build_transposed_rhs( N, N, A, N, x0, bt0 );
  build_transposed_rhs( N, N, A, N, x1, bt1 );

  {
    lapack_wrapper::TridiagonalLU<real_type> lu;
    UTILS_ASSERT( lu.factorize( N, L, D, U ), "TridiagonalLU::factorize failed\n" );

    real_type rhs[N];
    lapack_wrapper::copy( N, b0, 1, rhs, 1 );
    UTILS_ASSERT( lu.solve( rhs ), "TridiagonalLU::solve failed\n" );
    expect_vector_close( "TridiagonalLU solve", N, rhs, x0 );

    lapack_wrapper::copy( N, bt0, 1, rhs, 1 );
    UTILS_ASSERT( lu.t_solve( rhs ), "TridiagonalLU::t_solve failed\n" );
    expect_vector_close( "TridiagonalLU t_solve", N, rhs, x0 );

    lapack_wrapper::Matrix<real_type> RHS(N,2);
    for ( integer i{0}; i < N; ++i ) {
      RHS(i,0) = b0[i];
      RHS(i,1) = b1[i];
    }
    UTILS_ASSERT( lu.solve( RHS ), "TridiagonalLU::solve(MatrixWrapper) failed\n" );
    expect_vector_close( "TridiagonalLU multi rhs0", N, RHS.data(),     x0 );
    expect_vector_close( "TridiagonalLU multi rhs1", N, RHS.data() + N, x1 );

    for ( integer i{0}; i < N; ++i ) {
      RHS(i,0) = bt0[i];
      RHS(i,1) = bt1[i];
    }
    UTILS_ASSERT( lu.t_solve( RHS ), "TridiagonalLU::t_solve(MatrixWrapper) failed\n" );
    expect_vector_close( "TridiagonalLU multi transpose rhs0", N, RHS.data(),     x0 );
    expect_vector_close( "TridiagonalLU multi transpose rhs1", N, RHS.data() + N, x1 );

    UTILS_ASSERT(
      std::isfinite( lu.cond1( 9.5 ) ) && std::isfinite( lu.condInf( 9.75 ) ),
      "TridiagonalLU condition estimates must be finite\n"
    );
  }

  {
    real_type E[]{ -1.0, -0.5, -0.25, -0.125 };
    real_type DS[]{ 4.0, 4.5, 5.0, 5.5, 6.0 };
    real_type y0[]{ 2.0, -1.0, 0.5, 1.5, -0.25 };
    real_type y1[]{ -1.0, 0.5, 1.0, -2.0, 0.75 };
    real_type rhs0[N], rhs1[N];

    lapack_wrapper::TridiagonalSPD<real_type>::axpy( N, 1.0, E, DS, y0, 0.0, rhs0 );
    lapack_wrapper::TridiagonalSPD<real_type>::axpy( N, 1.0, E, DS, y1, 0.0, rhs1 );

    lapack_wrapper::TridiagonalSPD<real_type> spd;
    UTILS_ASSERT( spd.factorize( N, E, DS ), "TridiagonalSPD::factorize failed\n" );

    real_type rhs[N];
    lapack_wrapper::copy( N, rhs0, 1, rhs, 1 );
    UTILS_ASSERT( spd.solve( rhs ), "TridiagonalSPD::solve failed\n" );
    expect_vector_close( "TridiagonalSPD solve", N, rhs, y0 );

    lapack_wrapper::copy( N, rhs0, 1, rhs, 1 );
    UTILS_ASSERT( spd.t_solve( rhs ), "TridiagonalSPD::t_solve failed\n" );
    expect_vector_close( "TridiagonalSPD t_solve", N, rhs, y0 );

    lapack_wrapper::Matrix<real_type> RHS(N,2);
    for ( integer i{0}; i < N; ++i ) {
      RHS(i,0) = rhs0[i];
      RHS(i,1) = rhs1[i];
    }
    UTILS_ASSERT( spd.solve( RHS ), "TridiagonalSPD::solve(MatrixWrapper) failed\n" );
    expect_vector_close( "TridiagonalSPD multi rhs0", N, RHS.data(),     y0 );
    expect_vector_close( "TridiagonalSPD multi rhs1", N, RHS.data() + N, y1 );

    UTILS_ASSERT( std::isfinite( spd.cond1( 6.625 ) ), "TridiagonalSPD::cond1 invalid\n" );
  }

  {
    lapack_wrapper::TridiagonalQR<real_type> qr;
    UTILS_ASSERT( qr.factorize( N, L, D, U ), "TridiagonalQR::factorize failed\n" );

    real_type rhs[N];
    lapack_wrapper::copy( N, b0, 1, rhs, 1 );
    UTILS_ASSERT( qr.solve( rhs ), "TridiagonalQR::solve failed\n" );
    expect_vector_close( "TridiagonalQR solve", N, rhs, x0, 1e-7 );

    lapack_wrapper::copy( N, bt0, 1, rhs, 1 );
    UTILS_ASSERT( qr.t_solve( rhs ), "TridiagonalQR::t_solve failed\n" );
    expect_vector_close( "TridiagonalQR t_solve", N, rhs, x0, 1e-7 );

    lapack_wrapper::Matrix<real_type> RHS(N,2);
    for ( integer i{0}; i < N; ++i ) {
      RHS(i,0) = bt0[i];
      RHS(i,1) = bt1[i];
    }
    UTILS_ASSERT( qr.t_solve( RHS ), "TridiagonalQR::t_solve(MatrixWrapper) failed\n" );
    expect_vector_close( "TridiagonalQR multi transpose rhs0", N, RHS.data(),     x0, 1e-7 );
    expect_vector_close( "TridiagonalQR multi transpose rhs1", N, RHS.data() + N, x1, 1e-7 );

    lapack_wrapper::copy( N, b0, 1, rhs, 1 );
    qr.lsq( 1, rhs, N, 0.0 );
    expect_vector_close( "TridiagonalQR lsq", N, rhs, x0, 1e-7 );
  }

  {
    lapack_wrapper::TridiagonalQR<real_type> qr1;
    real_type dummy[]{ 0.0 };
    real_type d1[]{ 2.0 };
    real_type rhs[]{ 6.0 };
    UTILS_ASSERT( qr1.factorize( 1, dummy, d1, dummy ), "TridiagonalQR::factorize failed for N=1\n" );
    UTILS_ASSERT( qr1.solve( rhs ), "TridiagonalQR::solve failed for N=1\n" );
    expect_close( "TridiagonalQR N=1 solve", rhs[0], 3.0 );
  }
}

int
main() {
  try {
    test_lu_and_lupq_public_api();
    test_lsc_public_api();
    test_tridiagonal_public_api();
  } catch ( std::exception const & exc ) {
    msg.error( exc.what() );
    return 1;
  } catch ( ... ) {
    msg.error( "Unknown error!\n" );
    return 1;
  }

  msg.green("All extra coverage tests passed!\n");
  return 0;
}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
