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

static Utils::Console msg(&std::cout);

static constexpr real_type TEST_EPS{1e-8};

static
bool
close_enough(
  real_type a,
  real_type b,
  real_type rel = TEST_EPS,
  real_type abs = TEST_EPS
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
  real_type    tol = TEST_EPS
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
  real_type         tol = TEST_EPS
) {
  for ( integer i{0}; i < n; ++i ) expect_close( what, a[i], b[i], tol );
}

static
void
dense_matvec(
  integer           nr,
  integer           nc,
  real_type const * A,
  real_type const * x,
  real_type       * y
) {
  for ( integer i{0}; i < nr; ++i ) y[i] = 0;
  for ( integer j{0}; j < nc; ++j )
    for ( integer i{0}; i < nr; ++i )
      y[i] += A[i + j * nr] * x[j];
}

static
void
fill_rhs_matrix(
  integer                     nr,
  integer                     nc,
  std::vector<real_type> const & rhs0,
  std::vector<real_type> const & rhs1,
  lapack_wrapper::Matrix<real_type> & RHS
) {
  RHS.setup( nr, nc );
  for ( integer i{0}; i < nr; ++i ) {
    RHS(i,0) = rhs0[static_cast<size_t>(i)];
    RHS(i,1) = rhs1[static_cast<size_t>(i)];
  }
}

static
void
add_dense_block(
  integer           N,
  integer           irow,
  integer           icol,
  integer           nr,
  integer           nc,
  real_type const * B,
  integer           ldB,
  real_type       * A
) {
  for ( integer c{0}; c < nc; ++c )
    for ( integer r{0}; r < nr; ++r )
      A[(irow + r) + (icol + c) * N] = B[r + c * ldB];
}

static
void
test_banded_lu_api() {
  msg.green("test_banded_lu_api\n");

  constexpr integer N{4};
  constexpr real_type A[]{
    4, 1, 0, 0,
    1, 5, 1, 0,
    0, 1, 6, 1,
    0, 0, 1, 7
  };
  constexpr real_type B00[]{ 4, 1, 1, 5 };
  constexpr real_type B11[]{ 6, 1, 1, 7 };
  constexpr real_type x0[]{ 1.0, -1.0, 2.0, 0.5 };
  constexpr real_type x1[]{ -2.0, 0.5, 1.0, 3.0 };

  lapack_wrapper::BandedLU<real_type> blu;
  blu.setup( N, N, 1, 1 );
  blu.zero();
  blu.check( 0, 0 );
  blu.check( 3, 3 );
  blu.load_block( 2, 2, B00, 2, 0, 0 );
  blu.load_block( 2, 2, B11, 2, 2, 2 );
  blu.insert( 1, 2, 1.0, true );

  std::ostringstream stream;
  blu.dump( stream );
  UTILS_ASSERT( !stream.str().empty(), "BandedLU::dump should not be empty\n" );

  std::vector<real_type> rhs0( static_cast<size_t>(N) );
  std::vector<real_type> rhs1( static_cast<size_t>(N) );
  std::vector<real_type> axpy_rhs( static_cast<size_t>(N), 0 );

  dense_matvec( N, N, A, x0, rhs0.data() );
  dense_matvec( N, N, A, x1, rhs1.data() );
  blu.aAxpy( 1.0, x0, axpy_rhs.data() );
  expect_vector_close( "BandedLU aAxpy", N, axpy_rhs.data(), rhs0.data() );

  UTILS_ASSERT( blu.factorize(), "BandedLU::factorize failed\n" );

  std::vector<real_type> rhs_single( rhs0 );
  UTILS_ASSERT( blu.solve( rhs_single.data() ), "BandedLU::solve failed\n" );
  expect_vector_close( "BandedLU solve", N, rhs_single.data(), x0 );

  lapack_wrapper::Matrix<real_type> RHS;
  fill_rhs_matrix( N, 2, rhs0, rhs1, RHS );
  UTILS_ASSERT( blu.solve( RHS ), "BandedLU::solve(MatrixWrapper) failed\n" );
  expect_vector_close( "BandedLU multi rhs0", N, RHS.data(),     x0 );
  expect_vector_close( "BandedLU multi rhs1", N, RHS.data() + N, x1 );

  fill_rhs_matrix( N, 2, rhs0, rhs1, RHS );
  UTILS_ASSERT( blu.t_solve( RHS ), "BandedLU::t_solve(MatrixWrapper) failed\n" );
  expect_vector_close( "BandedLU transpose rhs0", N, RHS.data(),     x0 );
  expect_vector_close( "BandedLU transpose rhs1", N, RHS.data() + N, x1 );
}

static
void
exercise_banded_spd( lapack_wrapper::ULselect uplo ) {
  constexpr integer N{4};
  constexpr real_type A[]{
    4, 1, 0, 0,
    1, 5, 1, 0,
    0, 1, 6, 1,
    0, 0, 1, 7
  };
  constexpr real_type D[]{ 4, 5, 6, 7 };
  constexpr real_type OD[]{ 1, 1, 1 };
  constexpr real_type x0[]{ 1.0, -1.0, 2.0, 0.5 };
  constexpr real_type x1[]{ -2.0, 0.5, 1.0, 3.0 };

  lapack_wrapper::BandedSPD<real_type> spd;
  spd.setup( uplo, N, 1 );
  spd.zero();

  if ( uplo == lapack_wrapper::ULselect::LOWER ) {
    for ( integer j{0}; j < N; ++j ) spd(0,j) = D[j];
    for ( integer j{0}; j < N-1; ++j ) spd(1,j) = OD[j];
    UTILS_ASSERT( spd.factorize(), "BandedSPD::factorize failed for LOWER storage\n" );
  } else {
    for ( integer j{0}; j < N; ++j ) spd(1,j) = D[j];
    for ( integer j{1}; j < N; ++j ) spd(0,j) = OD[j-1];
    spd.factorize( "BandedSPD upper" );
  }

  std::vector<real_type> rhs0( static_cast<size_t>(N) );
  std::vector<real_type> rhs1( static_cast<size_t>(N) );
  dense_matvec( N, N, A, x0, rhs0.data() );
  dense_matvec( N, N, A, x1, rhs1.data() );

  std::vector<real_type> rhs_single( rhs0 );
  UTILS_ASSERT( spd.solve( rhs_single.data() ), "BandedSPD::solve failed\n" );
  expect_vector_close( "BandedSPD solve", N, rhs_single.data(), x0 );

  lapack_wrapper::Matrix<real_type> RHS;
  fill_rhs_matrix( N, 2, rhs0, rhs1, RHS );
  UTILS_ASSERT( spd.solve( RHS ), "BandedSPD::solve(MatrixWrapper) failed\n" );
  expect_vector_close( "BandedSPD multi rhs0", N, RHS.data(),     x0 );
  expect_vector_close( "BandedSPD multi rhs1", N, RHS.data() + N, x1 );

  fill_rhs_matrix( N, 2, rhs0, rhs1, RHS );
  UTILS_ASSERT( spd.t_solve( RHS ), "BandedSPD::t_solve(MatrixWrapper) failed\n" );
  expect_vector_close( "BandedSPD transpose rhs0", N, RHS.data(),     x0 );
  expect_vector_close( "BandedSPD transpose rhs1", N, RHS.data() + N, x1 );
}

static
void
test_banded_spd_api() {
  msg.green("test_banded_spd_api\n");
  exercise_banded_spd( lapack_wrapper::ULselect::LOWER );
  exercise_banded_spd( lapack_wrapper::ULselect::UPPER );
}

static
void
test_block_tridiagonal_wrapper_api() {
  msg.green("test_block_tridiagonal_wrapper_api\n");

  constexpr integer nblk{3};
  constexpr integer bs{2};
  constexpr integer N{nblk * bs};

  constexpr real_type D0[]{ 4.0, 1.0, 1.0, 3.0 };
  constexpr real_type D1[]{ 5.0, 0.0, 0.0, 4.0 };
  constexpr real_type D2[]{ 6.0, 1.0, 1.0, 5.0 };
  constexpr real_type L0[]{ 1.0, 0.0, 0.0, 1.0 };
  constexpr real_type L1[]{ 0.5, 0.0, 0.0, 0.5 };
  constexpr real_type x0[]{ 1.0, -1.0, 0.5, 2.0, -0.5, 1.5 };
  constexpr real_type x1[]{ -2.0, 0.25, 1.0, -1.5, 0.75, 2.5 };

  real_type A[N * N];
  lapack_wrapper::zero( N * N, A, 1 );
  add_dense_block( N, 0, 0, 2, 2, D0, 2, A );
  add_dense_block( N, 2, 2, 2, 2, D1, 2, A );
  add_dense_block( N, 4, 4, 2, 2, D2, 2, A );
  add_dense_block( N, 2, 0, 2, 2, L0, 2, A );
  add_dense_block( N, 0, 2, 2, 2, L0, 2, A );
  add_dense_block( N, 4, 2, 2, 2, L1, 2, A );
  add_dense_block( N, 2, 4, 2, 2, L1, 2, A );

  lapack_wrapper::BlockTridiagonalSymmetic<real_type> BT;
  BT.setup( nblk, bs );
  BT.zero();
  BT.setD( 0, D0, 2 );
  BT.setD( 1, D1, 2 );
  BT.setD( 2, D2, 2 );
  BT.setL( 0, L0, 2 );
  BT.setL( 1, L1, 2 );

  UTILS_ASSERT( BT.num_blocks() == nblk, "BlockTridiagonalSymmetic::num_blocks mismatch\n" );
  UTILS_ASSERT( BT.D_nrows(1) == bs,     "BlockTridiagonalSymmetic::D_nrows mismatch\n" );
  UTILS_ASSERT( BT.L_ncols(1) == bs,     "BlockTridiagonalSymmetic::L_ncols mismatch\n" );
  UTILS_ASSERT( BT.factorize(), "BlockTridiagonalSymmetic::factorize failed\n" );

  std::vector<real_type> rhs0( static_cast<size_t>(N) );
  std::vector<real_type> rhs1( static_cast<size_t>(N) );
  dense_matvec( N, N, A, x0, rhs0.data() );
  dense_matvec( N, N, A, x1, rhs1.data() );

  lapack_wrapper::Matrix<real_type> RHS;
  fill_rhs_matrix( N, 2, rhs0, rhs1, RHS );
  UTILS_ASSERT( BT.solve( RHS ), "BlockTridiagonalSymmetic::solve(MatrixWrapper) failed\n" );
  expect_vector_close( "BlockTridiagonal solve rhs0", N, RHS.data(),     x0 );
  expect_vector_close( "BlockTridiagonal solve rhs1", N, RHS.data() + N, x1 );

  fill_rhs_matrix( N, 2, rhs0, rhs1, RHS );
  UTILS_ASSERT( BT.t_solve( RHS ), "BlockTridiagonalSymmetic::t_solve(MatrixWrapper) failed\n" );
  expect_vector_close( "BlockTridiagonal transpose rhs0", N, RHS.data(),     x0 );
  expect_vector_close( "BlockTridiagonal transpose rhs1", N, RHS.data() + N, x1 );
}

int
main() {
  try {
    test_banded_lu_api();
    test_banded_spd_api();
    test_block_tridiagonal_wrapper_api();
    fmt::print( "All banded and block-tridiagonal coverage tests passed!\n" );
  } catch ( std::exception const & exc ) {
    msg.error( exc.what() );
    return 1;
  } catch ( ... ) {
    msg.error( "Unknown error!\n" );
    return 1;
  }
  return 0;
}
