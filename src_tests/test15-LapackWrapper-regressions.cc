#include <lapack_wrapper/lapack_wrapper.hh>
#include <lapack_wrapper/lapack_wrapper++.hh>

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
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
  for ( integer i{0}; i < n; ++i )
    expect_close( what, a[i], b[i], tol );
}

static
void
expect_zero_vector(
  char const *      what,
  integer           n,
  real_type const * x,
  real_type         tol = BIG_EPS
) {
  for ( integer i{0}; i < n; ++i ) {
    UTILS_ASSERT(
      std::isfinite(x[i]),
      "{} contains a non-finite value at {}\n",
      what, i
    );
    expect_close( what, x[i], real_type(0), tol );
  }
}

static
void
expect_finite_lower_3x3(
  char const *                               what,
  lapack_wrapper::QN<real_type> const &      H
) {
  for ( integer i{0}; i < 3; ++i )
    for ( integer j{0}; j <= i; ++j )
      UTILS_ASSERT(
        std::isfinite(H(i,j)),
        "{} contains a non-finite value at ({},{})\n",
        what, i, j
      );
}

static
void
test_lss_no_alloc_workspace_guard() {
  msg.green("test_lss_no_alloc_workspace_guard\n");

  constexpr integer NR{4};
  constexpr integer NC{3};
  constexpr integer NRC{NR*NC};
  constexpr integer LMIN{2*NRC + (NR < NC ? NR : NC)};

  lapack_wrapper::LSS_no_alloc<real_type> lss;
  std::vector<real_type> too_small(size_t(2*NRC));

  bool threw{false};
  try {
    lss.no_allocate( NR, NC, 2*NRC, too_small.data() );
  } catch ( std::exception const & ) {
    threw = true;
  }
  UTILS_ASSERT( threw, "LSS_no_alloc::no_allocate should reject undersized workspaces\n" );

  std::vector<real_type> storage( static_cast<size_t>(LMIN) );
  lss.no_allocate( NR, NC, LMIN, storage.data() );

  real_type A[]{
    1, 2, 0, 1,
    0, 1, 3, 2,
    2, 0, 1, 1
  };
  real_type x_true[]{ 1.5, -2.0, 0.5 };
  real_type xb[NR];
  lapack_wrapper::gemv( Transposition::NO, NR, NC, 1.0, A, NR, x_true, 1, 0.0, xb, 1 );

  UTILS_ASSERT( lss.factorize_nodim( A, NR ), "LSS_no_alloc::factorize_nodim failed\n" );
  UTILS_ASSERT( lss.solve( xb ), "LSS_no_alloc::solve failed\n" );
  expect_vector_close( "LSS_no_alloc solution", NC, xb, x_true );
}

static
void
test_lsy_workspace_and_transposed_tall_solve() {
  msg.green("test_lsy_workspace_and_transposed_tall_solve\n");

  constexpr integer NR{5};
  constexpr integer NC{3};
  constexpr integer NRC{NR*NC};

  lapack_wrapper::LSY_no_alloc<real_type> lsy_no_alloc;
  std::vector<real_type> work(size_t(2*NRC));
  std::vector<integer>   piv_small( static_cast<size_t>(NC) );

  bool threw{false};
  try {
    lsy_no_alloc.no_allocate( NR, NC, 2*NRC, work.data(), NC, piv_small.data() );
  } catch ( std::exception const & ) {
    threw = true;
  }
  UTILS_ASSERT(
    threw,
    "LSY_no_alloc::no_allocate should require a pivot workspace of size max(NR,NC)\n"
  );

  real_type A[]{
    1, 0, 2, 1, 3,
    0, 1, 1, 2, 1,
    1, 1, 0, 1, 1
  };

  real_type c1[]{ 1.0, -0.5, 2.0 };
  real_type c2[]{ -1.0, 0.75, 0.25 };
  real_type y1[NR], y2[NR];
  real_type b1[NC], b2[NC];

  lapack_wrapper::gemv( Transposition::NO,  NR, NC, 1.0, A, NR, c1, 1, 0.0, y1, 1 );
  lapack_wrapper::gemv( Transposition::NO,  NR, NC, 1.0, A, NR, c2, 1, 0.0, y2, 1 );
  lapack_wrapper::gemv( Transposition::YES, NR, NC, 1.0, A, NR, y1, 1, 0.0, b1, 1 );
  lapack_wrapper::gemv( Transposition::YES, NR, NC, 1.0, A, NR, y2, 1, 0.0, b2, 1 );

  lapack_wrapper::LSY<real_type> lsy;
  UTILS_ASSERT( lsy.factorize( NR, NC, A, NR ), "LSY::factorize failed\n" );

  real_type xb[NR]{ b1[0], b1[1], b1[2], 7.0, -3.0 };
  UTILS_ASSERT( lsy.t_solve( xb ), "LSY::t_solve failed on a tall matrix\n" );
  expect_vector_close( "LSY t_solve", NR, xb, y1 );

  lapack_wrapper::Matrix<real_type> XB(NR,2);
  XB.zero_fill();
  for ( integer i{0}; i < NC; ++i ) {
    XB(i,0) = b1[i];
    XB(i,1) = b2[i];
  }
  UTILS_ASSERT( lsy.t_solve( XB ), "LSY::t_solve(MatrixWrapper) failed\n" );
  expect_vector_close( "LSY multi t_solve rhs0", NR, XB.data(),        y1 );
  expect_vector_close( "LSY multi t_solve rhs1", NR, XB.data() + NR,   y2 );
}

static
void
test_svd_zero_and_tall_paths() {
  msg.green("test_svd_zero_and_tall_paths\n");

  {
    lapack_wrapper::SVD<real_type> svd;
    real_type A0[9]{};
    real_type xb[]{ 1.0, -2.0, 3.0 };
    real_type yb[]{ -4.0, 5.0, 6.0 };

    UTILS_ASSERT( svd.factorize( 3, 3, A0, 3 ), "SVD factorization of zero matrix failed\n" );
    UTILS_ASSERT( svd.solve( xb ),   "SVD::solve failed on zero matrix\n" );
    UTILS_ASSERT( svd.t_solve( yb ), "SVD::t_solve failed on zero matrix\n" );
    expect_zero_vector( "SVD zero solve",   3, xb );
    expect_zero_vector( "SVD zero t_solve", 3, yb );
  }

  {
    constexpr integer NR{5};
    constexpr integer NC{3};

    real_type A[]{
      1, 0, 2, 1, 3,
      0, 1, 1, 2, 1,
      1, 1, 0, 1, 1
    };

    lapack_wrapper::SVD<real_type> svd;
    UTILS_ASSERT( svd.factorize( NR, NC, A, NR ), "SVD factorization failed on tall matrix\n" );

    real_type x_true[]{ 0.5, -1.0, 2.0 };
    real_type rhs[NR];
    lapack_wrapper::gemv( Transposition::NO, NR, NC, 1.0, A, NR, x_true, 1, 0.0, rhs, 1 );

    real_type xb[NR];
    lapack_wrapper::copy( NR, rhs, 1, xb, 1 );
    UTILS_ASSERT( svd.solve( xb ), "SVD::solve failed on tall matrix\n" );
    expect_vector_close( "SVD tall solve", NC, xb, x_true );

    real_type y_true[NR];
    real_type rhs_t[NC];
    lapack_wrapper::gemv( Transposition::NO,  NR, NC, 1.0, A, NR, x_true, 1, 0.0, y_true, 1 );
    lapack_wrapper::gemv( Transposition::YES, NR, NC, 1.0, A, NR, y_true, 1, 0.0, rhs_t, 1 );

    real_type yb[NR]{ rhs_t[0], rhs_t[1], rhs_t[2], 9.0, -8.0 };
    UTILS_ASSERT( svd.t_solve( yb ), "SVD::t_solve failed on tall matrix\n" );
    expect_vector_close( "SVD tall t_solve", NR, yb, y_true );
  }
}

static
void
test_quasi_newton_degenerate_updates() {
  msg.green("test_quasi_newton_degenerate_updates\n");

  real_type x[]{ 1.0, -2.0, 0.5 };
  real_type s[]{ 1.0, 2.0, 3.0 };
  real_type y_zero[]{ 0.0, 0.0, 0.0 };
  real_type out[3];

  {
    lapack_wrapper::BFGS<real_type> bfgs;
    bfgs.allocate(3);
    bfgs.init();
    bfgs.update( y_zero, s, 1e-10 );
    expect_finite_lower_3x3( "BFGS degenerate lower", bfgs );
    bfgs.mult( x, out );
    expect_vector_close( "BFGS degenerate should preserve identity", 3, out, x );
  }

  {
    lapack_wrapper::DFP<real_type> dfp;
    dfp.allocate(3);
    dfp.init();
    dfp.update( y_zero, s, 1e-10 );
    expect_finite_lower_3x3( "DFP degenerate lower", dfp );
    dfp.mult( x, out );
    expect_vector_close( "DFP degenerate should preserve identity", 3, out, x );
  }

  {
    lapack_wrapper::SR1<real_type> sr1;
    sr1.allocate(3);
    sr1.init();
    sr1.update( s, s, 1e-10 );
    expect_finite_lower_3x3( "SR1 degenerate lower", sr1 );
    sr1.mult( x, out );
    expect_vector_close( "SR1 degenerate should preserve identity", 3, out, x );
  }
}

static
void
test_pinv_zero_rank() {
  msg.green("test_pinv_zero_rank\n");

  lapack_wrapper::PINV<real_type> pinv;
  real_type A[8]{};

  UTILS_ASSERT( pinv.factorize( 4, 2, A, 4 ), "PINV factorization of zero matrix failed\n" );

  real_type b[]{ 1.0, -2.0, 3.0, 4.0 };
  real_type x[2];
  UTILS_ASSERT( pinv.mult_inv( b, 1, x, 1 ), "PINV::mult_inv failed on zero matrix\n" );
  expect_zero_vector( "PINV mult_inv", 2, x );

  real_type bt[]{ -1.0, 5.0 };
  real_type xt[4];
  UTILS_ASSERT( pinv.t_mult_inv( bt, 1, xt, 1 ), "PINV::t_mult_inv failed on zero matrix\n" );
  expect_zero_vector( "PINV t_mult_inv", 4, xt );

  std::ostringstream stream;
  pinv.info( stream );
  UTILS_ASSERT(
    stream.str().find("rank = 0") != std::string::npos,
    "PINV info should report rank 0 for the zero matrix\n{}",
    stream.str()
  );
}

static
void
test_gsvd_sparse_matches_dense_rectangular() {
  msg.green("test_gsvd_sparse_matches_dense_rectangular\n");

  constexpr integer M{5};
  constexpr integer N{3};
  constexpr integer P{3};

  real_type A_data[]{
    1,  2,  3,  4,  5,
    6,  7,  8,  9, 10,
   11, 12, 13, 14, 15
  };
  real_type B_data[]{
    8, 3, 4,
    1, 5, 9,
    6, 7, 2
  };

  integer A_row[]{
    0,1,2,3,4,
    0,1,2,3,4,
    0,1,2,3,4
  };
  integer A_col[]{
    0,0,0,0,0,
    1,1,1,1,1,
    2,2,2,2,2
  };
  integer B_row[]{ 0,1,2,0,1,2,0,1,2 };
  integer B_col[]{ 0,0,0,1,1,1,2,2,2 };

  lapack_wrapper::MatrixWrapper<real_type> A(A_data, M, N, M);
  lapack_wrapper::MatrixWrapper<real_type> B(B_data, P, N, P);

  lapack_wrapper::GeneralizedSVD<real_type> dense(A, B);
  lapack_wrapper::GeneralizedSVD<real_type> sparse(
    M, N, P,
    15, A_data, A_row, A_col,
     9, B_data, B_row, B_col
  );

  for ( integer i{0}; i < N; ++i ) {
    expect_close( "GSVD alpha", dense.alpha(i), sparse.alpha(i), BIG_EPS );
    expect_close( "GSVD beta",  dense.beta(i),  sparse.beta(i),  BIG_EPS );
    expect_close(
      "GSVD alpha^2 + beta^2",
      dense.alpha(i) * dense.alpha(i) + dense.beta(i) * dense.beta(i),
      real_type(1),
      BIG_EPS
    );
    expect_close( "GSVD C diag", dense.getC()[i], dense.alpha(i), BIG_EPS );
    expect_close( "GSVD S diag", dense.gerS()[i], dense.beta(i), BIG_EPS );
  }

  UTILS_ASSERT( dense.getU().nrows() == M && dense.getU().ncols() == M, "bad GSVD U size\n" );
  UTILS_ASSERT( dense.getV().nrows() == P && dense.getV().ncols() == P, "bad GSVD V size\n" );
  UTILS_ASSERT( dense.getQ().nrows() == N && dense.getQ().ncols() == N, "bad GSVD Q size\n" );
  UTILS_ASSERT( !dense.info(1e-12).empty(), "GSVD info must not be empty\n" );
}

static
void
test_qr_mul_qt_matrixwrapper_dispatch() {
  msg.green("test_qr_mul_qt_matrixwrapper_dispatch\n");

  real_type A[]{
    2, 1, 0,
    1, 3, 1,
    0, 2, 4
  };
  real_type Cdata[]{
    1, 3,
    2, 4,
    5, 6
  };

  lapack_wrapper::QR<real_type> qr;
  UTILS_ASSERT( qr.factorize( 3, 3, A, 3 ), "QR factorization failed\n" );

  lapack_wrapper::Matrix<real_type> Cw(2,3);
  lapack_wrapper::Matrix<real_type> Cr(2,3);
  Cw.load( Cdata, 2 );
  Cr.load( Cdata, 2 );

  qr.mul_Qt( Cw );
  qr.mul_Qt( 2, 3, Cr.data(), Cr.ldim() );

  expect_vector_close( "QR mul_Qt(MatrixWrapper)", 6, Cw.data(), Cr.data() );
}

int
main() {
  try {
    test_lss_no_alloc_workspace_guard();
    test_lsy_workspace_and_transposed_tall_solve();
    test_svd_zero_and_tall_paths();
    test_pinv_zero_rank();
    test_gsvd_sparse_matches_dense_rectangular();
    test_qr_mul_qt_matrixwrapper_dispatch();
    test_quasi_newton_degenerate_updates();
  } catch ( std::exception const & exc ) {
    msg.error( exc.what() );
    return 1;
  } catch ( ... ) {
    msg.error( "Unknown error!\n" );
    return 1;
  }

  msg.green("All regression tests passed!\n");
  return 0;
}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
