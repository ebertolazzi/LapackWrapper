#include <lapack_wrapper/lapack_wrapper.hh>
#include <lapack_wrapper/lapack_wrapper++.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <vector>

template <typename TYPE> struct fmt::formatter<std::complex<TYPE>> : ostream_formatter {};

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

using real_type    = lapack_wrapper::doublereal;
using integer      = lapack_wrapper::integer;
using complex_type = std::complex<real_type>;
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
bool
close_enough(
  complex_type a,
  complex_type b,
  real_type    tol = BIG_EPS
) {
  return std::abs(a - b) <= tol;
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
expect_close(
  char const * what,
  complex_type a,
  complex_type b,
  real_type    tol = BIG_EPS
) {
  UTILS_ASSERT(
    close_enough(a, b, tol),
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
real_type
vector_norm(
  integer           n,
  complex_type const x[]
) {
  real_type sum{0};
  for ( integer i{0}; i < n; ++i ) sum += std::norm(x[i]);
  return std::sqrt(sum);
}

static
void
sort_complex_vector( std::vector<complex_type> & v ) {
  std::sort(
    v.begin(), v.end(),
    []( complex_type const & a, complex_type const & b ) {
      if ( a.real() < b.real() - BIG_EPS ) return true;
      if ( a.real() > b.real() + BIG_EPS ) return false;
      return a.imag() < b.imag();
    }
  );
}

static
void
test_matrix_and_diag_wrappers() {
  msg.green("test_matrix_and_diag_wrappers\n");

  lapack_wrapper::Matrix<real_type> M(3,4);
  M.fill(2.0);
  M.zero_block(1,2,1,1);
  M(0,0) = 5.0;
  M.add_to_diag(1.0);

  expect_close( "M(0,0)", M(0,0), 6.0 );
  expect_close( "M(1,1)", M(1,1), 1.0 );
  expect_close( "M(1,2)", M(1,2), 0.0 );

  lapack_wrapper::Matrix<real_type> MT(4,3);
  MT.load_transposed( M );
  for ( integer i{0}; i < M.nrows(); ++i )
    for ( integer j{0}; j < M.ncols(); ++j )
      expect_close( "load_transposed", MT(j,i), M(i,j) );

  lapack_wrapper::Matrix<real_type> Copy;
  Copy = MT;
  expect_vector_close( "Matrix copy assignment", MT.num_elems(), MT.data(), Copy.data() );

  lapack_wrapper::DiagMatrix<real_type> D(3);
  for ( integer i{0}; i < 3; ++i ) D[i] = 2.0;
  D.scale_by(0.5);
  lapack_wrapper::DiagMatrix<real_type> Dcopy(D);
  for ( integer i{0}; i < 3; ++i ) {
    expect_close( "DiagMatrix value", D[i],     1.0 );
    expect_close( "DiagMatrix copy",  Dcopy[i], 1.0 );
  }
  UTILS_ASSERT( !D.to_string().empty(), "DiagMatrix::to_string must not be empty\n" );

  real_type A_rank[]{
    1,0,0,0,
    0,2,0,0,
    0,0,0,0,
    0,0,0,0
  };
  real_type sval[3];
  integer rank{ lapack_wrapper::rankEstimate( 4, 4, A_rank, 4, real_type(1e-12), sval ) };
  UTILS_ASSERT( rank == 2, "rankEstimate should detect rank 2, found {}\n", rank );
}

static
void
test_qr_qrp_public_api() {
  msg.green("test_qr_qrp_public_api\n");

  {
    constexpr integer NR{4};
    constexpr integer NC{3};
    real_type A[]{
      1, 0, 2, 1,
      0, 1, 1, 2,
      1, 1, 0, 1
    };

    integer Lwork{ lapack_wrapper::QR_no_alloc<real_type>::get_Lwork_QR(NR,NC) + NR*NC + (NR < NC ? NR : NC) };
    std::vector<real_type> work( static_cast<size_t>(Lwork) );
    lapack_wrapper::QR_no_alloc<real_type> qr;
    qr.no_allocate( NR, NC, Lwork, work.data() );
    UTILS_ASSERT( qr.factorize_nodim( A, NR ), "QR_no_alloc::factorize_nodim failed\n" );

    lapack_wrapper::Matrix<real_type> Qr, R, Arec;
    qr.getQreduced( Qr );
    qr.getR( R );
    qr.getTau( Arec );
    UTILS_ASSERT( Arec.nrows() == 1 && Arec.ncols() == NC, "QR::getTau returned wrong size\n" );

    Arec.setup( NR, NC );
    lapack_wrapper::gemm( 1.0, Qr, R, 0.0, Arec );
    expect_vector_close( "QR reconstruction", NR*NC, Arec.data(), A, BIG_EPS );
  }

  {
    real_type A[]{
      2, 1, 0,
      1, 3, 2,
      0, 1, 4
    };

    lapack_wrapper::QR<real_type> qr;
    UTILS_ASSERT( qr.factorize( 3, 3, A, 3 ), "QR factorize failed\n" );

    lapack_wrapper::Matrix<real_type> left_w(3,2), left_r(3,2);
    lapack_wrapper::Matrix<real_type> right_w(2,3), right_r(2,3);
    left_w.fill(0.5); left_r.fill(0.5);
    right_w.fill(-1.0); right_r.fill(-1.0);

    qr.Q_mul( left_w );
    qr.Q_mul( 3, 2, left_r.data(), left_r.ldim() );
    expect_vector_close( "QR Q_mul wrapper", 6, left_w.data(), left_r.data() );

    left_w.fill(0.25); left_r.fill(0.25);
    qr.Qt_mul( left_w );
    qr.Qt_mul( 3, 2, left_r.data(), left_r.ldim() );
    expect_vector_close( "QR Qt_mul wrapper", 6, left_w.data(), left_r.data() );

    right_w.fill(1.0); right_r.fill(1.0);
    qr.mul_Q( right_w );
    qr.mul_Q( 2, 3, right_r.data(), right_r.ldim() );
    expect_vector_close( "QR mul_Q wrapper", 6, right_w.data(), right_r.data() );

    right_w.fill(-0.75); right_r.fill(-0.75);
    qr.mul_Qt( right_w );
    qr.mul_Qt( 2, 3, right_r.data(), right_r.ldim() );
    expect_vector_close( "QR mul_Qt wrapper", 6, right_w.data(), right_r.data() );

    left_w.fill(1.0); left_r.fill(1.0);
    qr.invR_mul( left_w );
    qr.invR_mul( 3, 2, left_r.data(), left_r.ldim() );
    expect_vector_close( "QR invR_mul wrapper", 6, left_w.data(), left_r.data() );

    left_w.fill(1.0); left_r.fill(1.0);
    qr.invRt_mul( left_w );
    qr.invRt_mul( 3, 2, left_r.data(), left_r.ldim() );
    expect_vector_close( "QR invRt_mul wrapper", 6, left_w.data(), left_r.data() );

    right_w.fill(1.0); right_r.fill(1.0);
    qr.mul_invR( right_w );
    qr.mul_invR( 2, 3, right_r.data(), right_r.ldim() );
    expect_vector_close( "QR mul_invR wrapper", 6, right_w.data(), right_r.data() );

    right_w.fill(1.0); right_r.fill(1.0);
    qr.mul_invRt( right_w );
    qr.mul_invRt( 2, 3, right_r.data(), right_r.ldim() );
    expect_vector_close( "QR mul_invRt wrapper", 6, right_w.data(), right_r.data() );

    real_type x1[]{ 1.0, -1.0, 0.5 };
    real_type x2[]{ -2.0, 1.5, 3.0 };
    lapack_wrapper::Matrix<real_type> RHS(3,2);
    lapack_wrapper::gemv( Transposition::NO, 3, 3, 1.0, A, 3, x1, 1, 0.0, RHS.data(),     1 );
    lapack_wrapper::gemv( Transposition::NO, 3, 3, 1.0, A, 3, x2, 1, 0.0, RHS.data() + 3, 1 );
    UTILS_ASSERT( qr.solve( RHS ), "QR::solve(MatrixWrapper) failed\n" );
    expect_vector_close( "QR solve matrix rhs0", 3, RHS.data(),     x1 );
    expect_vector_close( "QR solve matrix rhs1", 3, RHS.data() + 3, x2 );

    real_type y1[]{ 2.0, -0.5, 1.0 };
    real_type y2[]{ -1.0, 0.25, 0.5 };
    lapack_wrapper::Matrix<real_type> RT(3,2);
    lapack_wrapper::gemv( Transposition::YES, 3, 3, 1.0, A, 3, y1, 1, 0.0, RT.data(),     1 );
    lapack_wrapper::gemv( Transposition::YES, 3, 3, 1.0, A, 3, y2, 1, 0.0, RT.data() + 3, 1 );
    UTILS_ASSERT( qr.t_solve( RT ), "QR::t_solve(MatrixWrapper) failed\n" );
    expect_vector_close( "QR t_solve matrix rhs0", 3, RT.data(),     y1 );
    expect_vector_close( "QR t_solve matrix rhs1", 3, RT.data() + 3, y2 );
  }

  {
    constexpr integer NR{4};
    constexpr integer NC{4};
    real_type A[]{
      1, 0, 0, 0,
      0, 1, 0, 0,
      1, 0, 0, 0,
      0, 0, 0, 0
    };

    integer Lwork{ lapack_wrapper::QRP_no_alloc<real_type>::get_Lwork_QRP(NR,NC) + NR*NC + NR + NC };
    std::vector<real_type> work( static_cast<size_t>(Lwork) );
    std::vector<integer>   iwork( static_cast<size_t>(NC) );
    lapack_wrapper::QRP_no_alloc<real_type> qrp;
    qrp.no_allocate( NR, NC, Lwork, work.data(), NC, iwork.data() );
    UTILS_ASSERT( qrp.factorize_nodim( A, NR ), "QRP_no_alloc::factorize_nodim failed\n" );

    integer perm[NC];
    qrp.getPerm( perm );
    std::vector<real_type> v{ 1, 2, 3, 4 };
    std::vector<real_type> w(v);
    qrp.permute( w.data() );
    qrp.inv_permute( w.data() );
    expect_vector_close( "QRP permute/inv_permute", NC, w.data(), v.data() );

    lapack_wrapper::Matrix<real_type> Rows(NC,2), RowsRef(NC,2);
    lapack_wrapper::Matrix<real_type> Cols(2,NC), ColsRef(2,NC);
    for ( integer j{0}; j < Rows.ncols(); ++j )
      for ( integer i{0}; i < Rows.nrows(); ++i )
        Rows(i,j) = RowsRef(i,j) = real_type(1 + i + 10*j);
    for ( integer j{0}; j < Cols.ncols(); ++j )
      for ( integer i{0}; i < Cols.nrows(); ++i )
        Cols(i,j) = ColsRef(i,j) = real_type(1 + i + 10*j);

    qrp.permute_rows( Rows );
    qrp.permute_rows( NC, 2, RowsRef.data(), RowsRef.ldim() );
    expect_vector_close( "QRP permute_rows wrapper", Rows.num_elems(), Rows.data(), RowsRef.data() );
    qrp.inv_permute_rows( Rows );
    expect_vector_close( "QRP inv_permute_rows", Rows.num_elems(), Rows.data(), RowsRef.data() );

    qrp.permute_cols( Cols );
    qrp.permute_cols( 2, NC, ColsRef.data(), ColsRef.ldim() );
    expect_vector_close( "QRP permute_cols wrapper", Cols.num_elems(), Cols.data(), ColsRef.data() );
    qrp.inv_permute_cols( Cols );
    expect_vector_close( "QRP inv_permute_cols", Cols.num_elems(), Cols.data(), ColsRef.data() );

    integer rank{ qrp.rankEstimate(1e-12) };
    UTILS_ASSERT( rank == 2, "QRP rankEstimate should detect rank 2, found {}\n", rank );
  }
}

static
void
test_svd_public_api_and_matrix_overloads() {
  msg.green("test_svd_public_api_and_matrix_overloads\n");

  constexpr integer NR{5};
  constexpr integer NC{3};
  constexpr integer MINRC{3};
  real_type A[]{
    1, 0, 2, 1, 3,
    0, 1, 1, 2, 1,
    1, 1, 0, 1, 1
  };

  lapack_wrapper::SVD<real_type> svd;
  UTILS_ASSERT( svd.factorize( NR, NC, A, NR ), "SVD factorization failed\n" );

  for ( integer j{0}; j < MINRC; ++j ) {
    std::array<real_type,MINRC> ej{0,0,0};
    std::array<real_type,NR>    colU{0,0,0,0,0};
    std::array<real_type,NC>    colV{0,0,0};
    std::array<real_type,MINRC> recovered{0,0,0};

    ej[size_t(j)] = 1;

    svd.U_mul( 1.0, ej.data(), 1, 0.0, colU.data(), 1 );
    for ( integer i{0}; i < NR; ++i )
      expect_close( "SVD U accessor", colU[size_t(i)], svd.U(i,j) );
    svd.Ut_mul( 1.0, colU.data(), 1, 0.0, recovered.data(), 1 );
    for ( integer i{0}; i < MINRC; ++i )
      expect_close( "SVD Ut_mul", recovered[size_t(i)], i == j ? real_type(1) : real_type(0) );

    std::fill( recovered.begin(), recovered.end(), real_type(0) );
    svd.V_mul( 1.0, ej.data(), 1, 0.0, colV.data(), 1 );
    for ( integer i{0}; i < NC; ++i )
      expect_close( "SVD V accessor", colV[size_t(i)], svd.V(i,j) );
    svd.Vt_mul( 1.0, colV.data(), 1, 0.0, recovered.data(), 1 );
    for ( integer i{0}; i < MINRC; ++i )
      expect_close( "SVD Vt_mul", recovered[size_t(i)], i == j ? real_type(1) : real_type(0) );
  }

  real_type x1[]{ 1.0, -0.5, 2.0 };
  real_type x2[]{ -1.0, 1.5, 0.25 };
  lapack_wrapper::Matrix<real_type> XB(NR,2);
  lapack_wrapper::gemv( Transposition::NO, NR, NC, 1.0, A, NR, x1, 1, 0.0, XB.data(),      1 );
  lapack_wrapper::gemv( Transposition::NO, NR, NC, 1.0, A, NR, x2, 1, 0.0, XB.data() + NR, 1 );
  UTILS_ASSERT( svd.solve( XB ), "SVD::solve(MatrixWrapper) failed\n" );
  expect_vector_close( "SVD solve matrix rhs0", NC, XB.data(),      x1 );
  expect_vector_close( "SVD solve matrix rhs1", NC, XB.data() + NR, x2 );

  real_type y1[NR], y2[NR], b1[NC], b2[NC];
  lapack_wrapper::gemv( Transposition::NO,  NR, NC, 1.0, A, NR, x1, 1, 0.0, y1, 1 );
  lapack_wrapper::gemv( Transposition::NO,  NR, NC, 1.0, A, NR, x2, 1, 0.0, y2, 1 );
  lapack_wrapper::gemv( Transposition::YES, NR, NC, 1.0, A, NR, y1, 1, 0.0, b1, 1 );
  lapack_wrapper::gemv( Transposition::YES, NR, NC, 1.0, A, NR, y2, 1, 0.0, b2, 1 );

  lapack_wrapper::Matrix<real_type> YB(NR,2);
  YB.zero_fill();
  for ( integer i{0}; i < NC; ++i ) {
    YB(i,0) = b1[i];
    YB(i,1) = b2[i];
  }
  UTILS_ASSERT( svd.t_solve( YB ), "SVD::t_solve(MatrixWrapper) failed\n" );
  expect_vector_close( "SVD t_solve matrix rhs0", NR, YB.data(),      y1 );
  expect_vector_close( "SVD t_solve matrix rhs1", NR, YB.data() + NR, y2 );
}

static
void
test_eigen_dense_sparse_public_api() {
  msg.green("test_eigen_dense_sparse_public_api\n");

  real_type A[]{
    1,0,0,
    0,2,0,
    0,0,4
  };
  real_type B[]{
    1,0,0,
    0,2,0,
    0,0,3
  };

  integer rows[]{ 0,1,2 };
  integer cols[]{ 0,1,2 };
  real_type valsA[]{ 1,2,4 };
  real_type valsB[]{ 1,2,3 };

  {
    lapack_wrapper::Eigenvalues<real_type> dense(3, A, 3);
    lapack_wrapper::Eigenvalues<real_type> sparse(3, 3, valsA, rows, cols);
    std::vector<complex_type> ed, es;
    dense.getEigenvalues( ed );
    sparse.getEigenvalues( es );
    sort_complex_vector( ed );
    sort_complex_vector( es );

    std::vector<complex_type> expected{ {1,0}, {2,0}, {4,0} };
    for ( size_t i{0}; i < expected.size(); ++i ) {
      expect_close( "Eigenvalues dense",  ed[i], expected[i] );
      expect_close( "Eigenvalues sparse", es[i], expected[i] );
    }
  }

  {
    lapack_wrapper::Eigenvectors<real_type> eig(3, A, 3);
    std::vector<std::vector<complex_type>> right, left;
    eig.getRightEigenvector( right );
    eig.getLeftEigenvector( left );
    UTILS_ASSERT( right.size() == 3 && left.size() == 3, "bad eigenvector count\n" );

    for ( integer k{0}; k < 3; ++k ) {
      complex_type lambda;
      eig.getEigenvalue( k, lambda );
      complex_type rhs[3], lhs[3], tmp[3];
      for ( integer i{0}; i < 3; ++i ) {
        rhs[i] = right[size_t(k)][size_t(i)];
        lhs[i] = left[size_t(k)][size_t(i)];
        tmp[i] = complex_type(0);
      }
      tmp[0] = A[0] * rhs[0];
      tmp[1] = A[4] * rhs[1];
      tmp[2] = A[8] * rhs[2];
      for ( integer i{0}; i < 3; ++i ) tmp[i] -= lambda * rhs[i];
      UTILS_ASSERT( vector_norm(3, tmp) < BIG_EPS, "bad right eigenvector residual\n" );

      tmp[0] = A[0] * lhs[0];
      tmp[1] = A[4] * lhs[1];
      tmp[2] = A[8] * lhs[2];
      for ( integer i{0}; i < 3; ++i ) tmp[i] -= lambda * lhs[i];
      UTILS_ASSERT( vector_norm(3, tmp) < BIG_EPS, "bad left eigenvector residual\n" );
    }
  }

  {
    lapack_wrapper::GeneralizedEigenvalues<real_type> dense(3, A, 3, B, 3);
    lapack_wrapper::GeneralizedEigenvalues<real_type> sparse(
      3,
      3, valsA, rows, cols,
      3, valsB, rows, cols
    );
    std::vector<complex_type> ed, es;
    dense.getEigenvalues( ed );
    sparse.getEigenvalues( es );
    sort_complex_vector( ed );
    sort_complex_vector( es );

    std::vector<complex_type> expected{ {1,0}, {1,0}, {4.0/3.0,0} };
    sort_complex_vector( expected );
    for ( size_t i{0}; i < expected.size(); ++i ) {
      expect_close( "Generalized eigenvalues dense",  ed[i], expected[i] );
      expect_close( "Generalized eigenvalues sparse", es[i], expected[i] );
    }
  }

  {
    lapack_wrapper::GeneralizedEigenvectors<real_type> gev(3, A, 3, B, 3);
    std::vector<std::vector<complex_type>> right, left;
    gev.getRightEigenvector( right );
    gev.getLeftEigenvector( left );
    UTILS_ASSERT( right.size() == 3 && left.size() == 3, "bad generalized eigenvector count\n" );

    for ( integer k{0}; k < 3; ++k ) {
      complex_type lambda;
      gev.getEigenvalue( k, lambda );
      complex_type vr[3], vl[3], tmp[3];
      for ( integer i{0}; i < 3; ++i ) {
        vr[i] = right[size_t(k)][size_t(i)];
        vl[i] = left[size_t(k)][size_t(i)];
      }
      tmp[0] = A[0]*vr[0] - lambda*B[0]*vr[0];
      tmp[1] = A[4]*vr[1] - lambda*B[4]*vr[1];
      tmp[2] = A[8]*vr[2] - lambda*B[8]*vr[2];
      UTILS_ASSERT( vector_norm(3, tmp) < BIG_EPS, "bad generalized right eigenvector residual\n" );

      tmp[0] = A[0]*vl[0] - lambda*B[0]*vl[0];
      tmp[1] = A[4]*vl[1] - lambda*B[4]*vl[1];
      tmp[2] = A[8]*vl[2] - lambda*B[8]*vl[2];
      UTILS_ASSERT( vector_norm(3, tmp) < BIG_EPS, "bad generalized left eigenvector residual\n" );
    }
  }
}

template <typename HessianType>
static
void
check_finite_lower_3x3(
  char const *        name,
  HessianType const & H
) {
  for ( integer i{0}; i < 3; ++i )
    for ( integer j{0}; j <= i; ++j )
      UTILS_ASSERT(
        std::isfinite(H(i,j)),
        "{} contains a non-finite entry at ({},{})\n",
        name, i, j
      );
}

static
void
test_quasi_newton_public_api() {
  msg.green("test_quasi_newton_public_api\n");

  real_type x[]{ 1.0, -2.0, 0.5 };
  real_type y[]{ 11.0/3.0, -2.0, 13.0/3.0 };
  real_type s[]{ 1.0, 2.0, 3.0 };
  real_type f0[]{ 1.0, 1.0, 1.0 };
  real_type f1[]{ 2.0, 3.0, 4.0 };
  real_type x0[]{ 0.0, 0.0, 0.0 };
  real_type x1[]{ 1.0, 2.0, 3.0 };
  real_type out[3];

  lapack_wrapper::BFGS<real_type> bfgs;
  bfgs.allocate(3);
  bfgs.init();
  bfgs.mult( x, out );
  expect_vector_close( "BFGS init mult", 3, out, x );
  bfgs.minus_mult( x, out );
  for ( integer i{0}; i < 3; ++i ) expect_close( "BFGS minus_mult", out[i], -x[i] );
  bfgs.zero();
  bfgs.mult( x, out );
  for ( integer i{0}; i < 3; ++i ) expect_close( "BFGS zero mult", out[i], 0.0 );
  bfgs.init();
  bfgs.update( y, s, 1e-10 );
  UTILS_ASSERT( !bfgs.to_string().empty(), "BFGS::to_string must not be empty\n" );
  check_finite_lower_3x3( "BFGS lower", bfgs );

  lapack_wrapper::DFP<real_type> dfp;
  dfp.allocate(3);
  dfp.init();
  dfp.update( f0, f1, x0, x1, 1e-10 );
  dfp.mult( x, out );
  UTILS_ASSERT( !dfp.to_string().empty(), "DFP::to_string must not be empty\n" );
  check_finite_lower_3x3( "DFP lower", dfp );

  lapack_wrapper::SR1<real_type> sr1;
  sr1.allocate(3);
  sr1.init();
  sr1.update( f0, f1, s, 1e-10 );
  sr1.mult( x, out );
  UTILS_ASSERT( !sr1.to_string().empty(), "SR1::to_string must not be empty\n" );
  check_finite_lower_3x3( "SR1 lower", sr1 );
}

int
main() {
  try {
    test_matrix_and_diag_wrappers();
    test_qr_qrp_public_api();
    test_svd_public_api_and_matrix_overloads();
    test_eigen_dense_sparse_public_api();
    test_quasi_newton_public_api();
  } catch ( std::exception const & exc ) {
    msg.error( exc.what() );
    return 1;
  } catch ( ... ) {
    msg.error( "Unknown error!\n" );
    return 1;
  }

  msg.green("All coverage tests passed!\n");
  return 0;
}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
