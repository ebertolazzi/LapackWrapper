/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include <iostream>
#include <vector>
#include <lapack_wrapper/lapack_wrapper.hh>
#include <lapack_wrapper/lapack_wrapper++.hh>

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
#endif

using namespace std;
typedef double valueType;
using lapack_wrapper::integer;

static Utils::Console msg(&std::cout);

static
void
test1() {

  lapack_wrapper::QR<valueType> qr;

  integer const M   = 3;
  integer const N   = 5;
  integer const LDA = 3;
  valueType A[] = {
    0.001,      2,     3,
    0.001,  0.001,     0,
    0,      0.001,     0,
    0.001,     -1,     0,
    0.000001,   5,     3
  };

  lapack_wrapper::MatrixWrapper<valueType> Amat( A, M, N, LDA );

  msg.green(
    "\n\n\nTest1:\n\nInitial A\n" +
    lapack_wrapper::print_matrix( Amat )
  );

  msg.green( "\nDo QR factorization of A^T\n" );
  qr.t_factorize( "qr", Amat );

  valueType R[M*M];
  qr.getR( R, M );
  fmt::print( "\nR=\n{}", lapack_wrapper::print_matrix( M, M, R, M ) );

  valueType rhs[M], b[M];
  valueType x[N] = {1,2,3,4,5};
  lapack_wrapper::gemv( 1.0, Amat, x, 1, 0.0, rhs, 1 );
  lapack_wrapper::gemv( 1.0, Amat, x, 1, 0.0, b,   1 );

  fmt::print(
    "\nLS solution of A x = b\n"
    "b^T      = {}",
    lapack_wrapper::print_matrix( 1, M, b, 1 )
  );

  qr.invRt_mul( rhs, 1 );
  lapack_wrapper::copy( M, rhs, 1, x, 1 );
  lapack_wrapper::zero( N-M, x+3, 1 );
  qr.Q_mul( x );

  fmt::print(
    "x^T      = {}",
    lapack_wrapper::print_matrix( 1, M, x, 1 )
  );

  lapack_wrapper::gemv( -1.0, Amat, x, 1, 1.0, b, 1 );
  valueType res = lapack_wrapper::nrm2( M, b, 1 );
  UTILS_ASSERT( res < 1e-6, "test failed! res = {}\n", res );

  fmt::print(
    "residual = {}\n||res||_2 = {}\ndone test1\n",
    lapack_wrapper::print_matrix( 1, M, b, 1 ), res
  );
}



static
void
test2() {
  lapack_wrapper::QRP<valueType> qr;
  integer const M   = 3;
  integer const N   = 5;
  integer const LDA = 3;
  valueType A[] = {
    0.001,      2,     3,
    0.001,  0.001,     0,
    0,      0.001,     0,
    0.001,     -1,     0,
    0.000001,   5,     3
  };

  lapack_wrapper::MatrixWrapper<valueType> Amat( A, M, N, LDA );

  msg.green( "\n\n\nTest2:\n\nInitial A\n"+ lapack_wrapper::print_matrix( Amat ) );

  fmt::print("\nDo QR factorization of A^T\n");
  qr.t_factorize( "qrp", Amat );

  valueType R[M*M];
  qr.getR( R, M );
  fmt::print( "\nR=\n{}", lapack_wrapper::print_matrix( M, M, R, M ) );

  valueType rhs[M], b[M];
  valueType x[N] = {1,2,3,4,5};
  lapack_wrapper::gemv( 1.0, Amat, x, 1, 0.0, rhs, 1 );
  lapack_wrapper::gemv( 1.0, Amat, x, 1, 0.0, b,   1 );

  fmt::print(
    "\nLS solution of A x = b\n\n"
    "b^T =      {}",
    lapack_wrapper::print_matrix( 1, M, b, 1 )
  );

  qr.inv_permute( rhs );
  qr.invRt_mul( rhs, 1 );
  lapack_wrapper::copy( M,   rhs, 1, x, 1 );
  lapack_wrapper::zero( N-M, x+3, 1 );
  qr.Q_mul( x );

  fmt::print(
    "x^T =      {}",
    lapack_wrapper::print_matrix( 1, M, x, 1 )
  );

  lapack_wrapper::gemv( -1.0, Amat, x, 1, 1.0, b, 1 );
  valueType res = lapack_wrapper::nrm2( M, b, 1 );
  UTILS_ASSERT( res < 1e-6, "test failed! res = {}\n", res );

  fmt::print(
    "residual  = {}"
    "||res||_2 = {}\n"
    "done test2\n",
    lapack_wrapper::print_matrix( 1, M, b, 1 ), res
  );
}

static
void
test3() {
  lapack_wrapper::QRP<valueType> qr;
  integer const M   = 5;
  integer const N   = 5;
  integer const LDA = 5;
  valueType A[] = {
    0.001,      2,     3,     2, 3,
    0.001,  0.001,     0, 0.001, 1e-10,
    0,      0.001,     0, 0.001, 1e-12,
    0.001,     -1, 1e-12,    -1, -1e-12,
    0.000001,   5,     3,     5, 3
  };

  msg.green(
    "\n\n\nTest3:\n\nInitial A\n" + lapack_wrapper::print_matrix( M, N, A, M )
  );

  fmt::print("\nDo QR factorization of A^T\n");
  qr.t_factorize( "qr", M, N, A, LDA );

  valueType R[M*M];
  qr.getR( R, M );
  fmt::print( "\nR=\n{}", lapack_wrapper::print_matrix( M, M, R, M ) );

  valueType rhs[M], b[M];
  valueType x[N] = {1,2,3,4,5};
  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, N, 1, A, LDA, x, 1, 0, rhs, 1
  );
  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, N, 1, A, LDA, x, 1, 0, b, 1
  );

  fmt::print(
    "\nLS solution of A x = b\n\n"
    "b^T =      {}",
    lapack_wrapper::print_matrix( 1, M, b, 1 )
  );

  qr.inv_permute( rhs ); // da aggiungere!
  qr.invRt_mul( rhs, 1 );
  lapack_wrapper::copy( 3, rhs, 1, x, 1 );
  lapack_wrapper::zero( 2, x+3, 1 );
  qr.Q_mul( x );

  fmt::print(
    "x^T =      {}",
    lapack_wrapper::print_matrix( 1, M, x, 1 )
  );

  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, N, -1, A, LDA, x, 1, 1, b, 1
  );

  valueType res = lapack_wrapper::nrm2( M, b, 1 );
  UTILS_ASSERT( res < 1e-6, "test failed! res = {}\n", res );

  fmt::print(
    "residual = {}\n||res||_2 = {}\ndone test3\n",
    lapack_wrapper::print_matrix( 1, M, b, 1 ), res
  );
}

#define TEST4(NAME,F) \
  cout << "\n\nDo " << NAME << " factorization of A\n"; \
  F.factorize( NAME, M, M, A, LDA ); \
  \
  cout << NAME << " solution of A x = b\n"; \
  lapack_wrapper::copy( M, rhs, 1, x, 1 ); \
  lapack_wrapper::copy( M, rhs, 1, b, 1 ); \
  /* L.solve( x ); */ \
  F.solve( 1, x, M); \
  cout << "x^T      ="  \
       << lapack_wrapper::print_matrix( 1, M, x, 1 ); \
  \
  lapack_wrapper::gemv( lapack_wrapper::NO_TRANSPOSE, M, M, -1, A, LDA, x, 1, 1, b, 1 ); \
  cout << "residual =" \
       << lapack_wrapper::print_matrix( 1, M, b, 1 ); \
  res = lapack_wrapper::nrm2( M, b, 1 ); \
  cout << "||res||_2 = " << res << '\n'; \
  UTILS_ASSERT0( res < 1e-6, "test failed!\n" );

static
void
test4() {
  lapack_wrapper::LU<valueType>   lu;
  lapack_wrapper::LUPQ<valueType> lupq;
  lapack_wrapper::QR<valueType>   qr;
  lapack_wrapper::QRP<valueType>  qrp;
  lapack_wrapper::SVD<valueType>  svd;
  lapack_wrapper::LSS<valueType>  lss;
  lapack_wrapper::LSY<valueType>  lsy;

  integer const M   = 5;
  integer const LDA = 5;
  valueType A[] = {
    0.001,      2,     3,       2,      3,
    0.001,  0.001,     0,   0.001,  1e-10,
    0,      0.001,     0,   0.001,  1e-12,
    0.001,     -1, 1e-6+1,     -1, -1e-12,
    0.000001,   5,     3,  1e-6+5,    3+1
  };

  valueType rhs[M], b[M], res;
  valueType x[M] = {1,2,3,4,5};
  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, M, 1, A, LDA, x, 1, 0, rhs, 1
  );
  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, M, 1, A, LDA, x, 1, 0, b, 1
  );

  msg.green( fmt::format(
    "\n\n\nTest4:\n\nInitial A\n{}\n",
    lapack_wrapper::print_matrix( M, M, A, M )
  ) );

  TEST4("LU",lu);
  TEST4("LUPQ",lupq);
  TEST4("QR",qr);
  TEST4("QRP",qrp);
  TEST4("SVD",svd);
  TEST4("LSS",lss);
  TEST4("LSY",lsy);

  cout << "done test4\n";
}

#define TEST5(NAME,F) \
  cout << "\n\nDo " << NAME << " factorization of A\n"; \
  F.factorize( NAME, M, M, A, LDA ); \
  \
  cout << NAME << " solution of A x = b\n"; \
  lapack_wrapper::copy( M, rhs, 1, x, 1 ); \
  lapack_wrapper::copy( M, rhs, 1, b, 1 ); \
  /* F.t_solve( x ); */ \
  F.t_solve( 1, x, M); \
  cout << "x^T      = " \
       << lapack_wrapper::print_matrix( 1, M, x, 1 ); \
  \
  lapack_wrapper::gemv( lapack_wrapper::TRANSPOSE, M, M, -1, A, LDA, x, 1, 1, b, 1 ); \
  cout << "residual = " \
       << lapack_wrapper::print_matrix( 1, M, b, 1 ); \
  res = lapack_wrapper::nrm2( M, b, 1 ); \
  cout << "||res||_2 = " << res << '\n'; \
  UTILS_ASSERT0( res < 1e-6, "test failed!\n" );


static
void
test5() {
  lapack_wrapper::LU<valueType>   lu;
  lapack_wrapper::LUPQ<valueType> lupq;
  lapack_wrapper::QR<valueType>   qr;
  lapack_wrapper::QRP<valueType>  qrp;
  lapack_wrapper::SVD<valueType>  svd;
  lapack_wrapper::LSS<valueType>  lss;
  lapack_wrapper::LSY<valueType>  lsy;

  integer const M   = 5;
  integer const LDA = 5;
  valueType A[] = {
    0.001,      2,     3,       2,      3,
    0.001,  0.001,     0,   0.001,  1e-10,
    0,      0.001,     0,   0.001,  1e-12,
    0.001,     -1, 1e-6+1,     -1, -1e-12,
    0.000001,   5,     3,  1e-6+5,     3+1
  };

  valueType rhs[M], b[M], res;
  valueType x[M] = {1,2,3,4,5};
  lapack_wrapper::gemv(
    lapack_wrapper::TRANSPOSE, M, M, 1, A, LDA, x, 1, 0, rhs, 1
  );
  lapack_wrapper::gemv(
    lapack_wrapper::TRANSPOSE, M, M, 1, A, LDA, x, 1, 0, b, 1
  );

  msg.green( fmt::format(
    "\n\n\nTest5:\n\nInitial A\n{}\n",
    lapack_wrapper::print_matrix( M, M, A, M )
  ) );

  TEST5("LU",lu);
  TEST5("LUPQ",lupq);
  TEST5("QR",qr);
  TEST5("QRP",qrp);
  TEST5("SVD",svd);
  TEST5("LSS",lss);
  TEST5("LSY",lsy);

  cout << "\ndone test5\n";
}


static
void
test6() {
  lapack_wrapper::TridiagonalLU<valueType> lu;
  lapack_wrapper::TridiagonalQR<valueType> qr;

  integer const N = 5;
  valueType const D[] = { 1, 1, 2, -0.1, 0.1 };
  valueType const L[] = { -0.1, -1, -2, -0.1 };
  valueType const U[] = { -1, -10, -2, 0.1 };

  valueType rhs[N], b[N];
  valueType x[N] = {1,2,3,4,5};
  qr.axpy( N, 1.0, L, D, U, x, 0.0, rhs );
  qr.axpy( N, 1.0, L, D, U, x, 0.0, b   );

  msg.green( "\n\n\nTest6:\n\n" );

  cout << "\n\nDo (trid) LU  factorization of A\n";
  lu.factorize( "lu", N, L, D, U );

  cout << "(trid) LU solution of A x = b\n";
  lapack_wrapper::copy( N, rhs, 1, x, 1 );
  lapack_wrapper::copy( N, rhs, 1, b, 1 );
  lu.solve( x );
  //qr.t_solve( 1, x, M);
  cout
    << "x^T =      "
    << lapack_wrapper::print_matrix( 1, N, x, 1 );

  lu.axpy( N, -1.0, L, D, U, x, 1.0, b );
  cout
    << "residual = "
    << lapack_wrapper::print_matrix( 1, N, b, 1 );

  cout << "\n\nDo (trid) QR  factorization of A\n";
  qr.factorize( "qr", N, L, D, U );

  cout << "(trid) QR solution of A x = b\n";
  lapack_wrapper::copy( N, rhs, 1, x, 1 );
  lapack_wrapper::copy( N, rhs, 1, b, 1 );
  qr.solve( x );
  //qr.t_solve( 1, x, M);
  cout
    << "x^T      = "
    << lapack_wrapper::print_matrix( 1, N, x, 1 );

  qr.axpy( N, -1.0, L, D, U, x, 1.0, b );

  valueType res = lapack_wrapper::nrm2( N, b, 1 );
  UTILS_ASSERT( res < 1e-6, "test failed! res = {}\n", res );

  msg.green(
    fmt::format(
      "residual = {}\n||res||_2 = {}\ndone test6\n",
      lapack_wrapper::print_matrix( 1, N, b, 1 ), res
    )
  );
}

static
void
test7() {
  lapack_wrapper::QRP<valueType> qrp;

  integer const M   = 5;
  integer const LDA = 5;
  valueType A[] = {
    0.001,      2,     3,       2,      3,
    0.001,  0.001,     0,   0.001,  1e-10,
    0,      0.001,     0,   0.001,  1e-12,
    0.001,     -1, 1e-6+1,     -1, -1e-12,
    0.000001,   5,     3,  1e-6+5,    3+1
  };

  valueType rhs[M], b[M];
  valueType x[M] = {1,2,3,4,5};
  lapack_wrapper::gemv(
    lapack_wrapper::TRANSPOSE, M, M, 1, A, LDA, x, 1, 0, rhs, 1
  );
  lapack_wrapper::gemv(
    lapack_wrapper::TRANSPOSE, M, M, 1, A, LDA, x, 1, 0, b, 1
  );

  msg.green( fmt::format(
    "\n\n\nTest7:\n\nInitial A\n{}\n",
    lapack_wrapper::print_matrix( M, M, A, M )
  ) );

  cout << "\n\nDo QRP factorization of A\n";
  qrp.factorize( "qrp", M, M, A, LDA );

  cout << "QRP solution of A x = b\n";
  lapack_wrapper::copy( M, rhs, 1, x, 1 );
  lapack_wrapper::copy( M, rhs, 1, b, 1 );
  qrp.t_solve( x );
  //qrp.t_solve( 1, x, M );
  cout
    << "x^T      = "
    << lapack_wrapper::print_matrix( 1, M, x, 1 );

  lapack_wrapper::gemv(
    lapack_wrapper::TRANSPOSE, M, M, -1, A, LDA, x, 1, 1, b, 1
  );

  valueType res = lapack_wrapper::nrm2( M, b, 1 );
  UTILS_ASSERT( res < 1e-6, "test failed! res = {}\n", res );

  msg.green(
    fmt::format(
      "residual = {}\n||res||_2 = {}\ndone test7\n",
      lapack_wrapper::print_matrix( 1, M, b, 1 ), res
    )
  );
}

static
void
test8() {
  lapack_wrapper::LSC<valueType> lsc;

  integer const M   = 5;
  integer const LDA = 5;
  valueType A[] = {
    0.001,      2,     3,       2,      3,
    0.001,  0.001,     0,   0.001,  1e-10,
    0,      0.001,     0,   0.001,  1e-12,
    0.001,     -1, 1e-6+1,     -1, -1e-12,
    0.000001,   5,     3,  1e-6+5,    3+1
  };

  valueType rhs[M], b[M];
  valueType x[M] = {1,2,3,4,5};
  //bool      rselect[M] = { true, true, true, true, true };
  //bool      rselect[M] = { false, false, false, false, false };
  bool      rselect[M] = { true, false, true, true, true };
  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, M, 1, A, LDA, x, 1, 0, rhs, 1
  );
  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, M, 1, A, LDA, x, 1, 0, b, 1
  );

  msg.green( fmt::format(
    "\n\n\nTest8:\n\nInitial A\n{}\nb^T{}\n",
    lapack_wrapper::print_matrix( M, M, A, M ),
    lapack_wrapper::print_matrix( 1, M, b, 1 )
  ) );

  cout << "\n\nDo LSC factorization of A\n";
  lsc.factorize( M, M, A, LDA, rselect );

  cout << "\nLSC solution of A x = b\n";
  lapack_wrapper::copy( M, rhs, 1, x, 1 );
  lapack_wrapper::copy( M, rhs, 1, b, 1 );
  lsc.solve( x );
  //qrp.t_solve( 1, x, M );
  cout
    << "x^T      = "
    << lapack_wrapper::print_matrix( 1, M, x, 1 );

  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, M, -1, A, LDA, x, 1, 1, b, 1
  );

  valueType res = lapack_wrapper::nrm2( M, b, 1 );
  UTILS_ASSERT( res < 1e-6, "test failed! res = {}\n", res );

  msg.green(
    fmt::format(
      "residual = {}\n||res||_2 = {}\ndone test8\n",
      lapack_wrapper::print_matrix( 1, M, b, 1 ), res
    )
  );
}


static
void
test9() {
  lapack_wrapper::PINV<valueType> pinv;

  integer const M   = 5;
  integer const LDA = 5;
  valueType A[] = {
    0.001,    1,     0,      2,  1,
    1e-9, -1e+9,  1e-9,  -1e-9,  1,
    1e-9, -1e+9,     1,  -1e-9,  1,
    1e-9, -1e+9,  1e-9,      1,  1,
    //1e-9, -1e+9,  1e-9,  -1e-9,  1
    0, 0, 0, 0, 0
  };

  valueType rhs[M], b[M], x[M];
  valueType const xe[M] = {1,2,3,4,5};
  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, M, 1, A, LDA, xe, 1, 0, rhs, 1
  );
  lapack_wrapper::copy( M, rhs, 1, b, 1 );

  msg.green( fmt::format(
    "\n\n\nTest9:\n\nInitial A\n{}\nb^T\n",
    lapack_wrapper::print_matrix( M, M, A, M ),
    lapack_wrapper::print_matrix( 1, M, b, 1 )
  ) );

  cout << "\n\nDo PINV factorization of A\n";
  pinv.factorize( M, M, A, LDA );

  cout << "\nLSC solution of A x = b\n";
  lapack_wrapper::copy( M, b, 1, x, 1 );
  pinv.solve( x );
  //pinv.t_solve( 1, x, M );
  cout
    << "x^T      = "
    << lapack_wrapper::print_matrix( 1, M, x, 1 );

  valueType e[M] = { x[0]-xe[0],x[1]-xe[1],x[2]-xe[2],x[3]-xe[3],x[4]-xe[4]};

  cout
    << "e^T      = "
    << lapack_wrapper::print_matrix( 1, M, e, 1 );

  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, M, -1, A, LDA, x, 1, 1, b, 1
  );

  valueType res = lapack_wrapper::nrm2( M, b, 1 );
  UTILS_ASSERT( res < 1e-6, "test failed! res = {}\n", res );

  fmt::print(
    "residual = {}\n||res||_2 = {}\n",
    lapack_wrapper::print_matrix( 1, M, b, 1 ), res
  );

  cout << "\n\nDo PINV factorization of A\n";
  pinv.factorize( M, M-1, A, LDA );
  pinv.mult_inv( rhs, 1, x, 1 );

  lapack_wrapper::copy( M, rhs, 1, b, 1 );
  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, M-1, -1, A, LDA, x, 1, 1, b, 1
  );
  res = lapack_wrapper::nrm2( M, b, 1 );

  fmt::print(
    "x^T       = {}"
    "residual  = {}\n"
    "||res||_2 = {}\n"
    "test9\n",
    lapack_wrapper::print_matrix( 1, M-1, x, 1 ),
    lapack_wrapper::print_matrix( 1, M, b, 1 ),
    res
  );
}

static
void
test10() {
  lapack_wrapper::PINV<valueType> pinv;

  integer const M   = 6;
  integer const N   = 4;
  integer const LDA = 6;
  valueType A[] = {
    2,   1,  2,  1, 1, 1,
    2,   1,  2,  1, 1, 1,
    20, 10, 20, 10, 1, 1,
     1,  2,  3,  4, 5, 6,
  };

  valueType rhs[M], b[M], b2[2*M], x[N], x2[2*N];
  valueType const xe[N] = {1,2,3,4};
  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, N, 1, A, LDA, xe, 1, 0, rhs, 1
  );
  lapack_wrapper::copy( M, rhs, 1, b, 1 );
  lapack_wrapper::copy( M, rhs, 1, b2, 1 );
  lapack_wrapper::copy( M, rhs, 1, b2+M, 1 );
  msg.green( fmt::format(
    "\n\n\nTest10:\n\nInitial A\n{}\nb2{}\n",
    lapack_wrapper::print_matrix( M, N, A, M ),
    lapack_wrapper::print_matrix( M, 2, b2, M )
  ) );

  cout << "\n\nDo PINV factorization of A\n";
  pinv.factorize( M, N, A, LDA );

  cout << "\nPINV solution of A x = b\n";
  pinv.mult_inv( b, 1, x, 1 );
  pinv.mult_inv( 2, b2, M, x2, N );

  //pinv.t_solve( 1, x, M );
  cout
    << "x^T      = "
    << lapack_wrapper::print_matrix( 1, N, x, 1 )
    << "x2^T     = "
    << lapack_wrapper::print_matrix( 1, N, x2, 1 );

  valueType e[N] = { x[0]-xe[0], x[1]-xe[1], x[2]-xe[2] };

  cout
    << "e^T      = "
    << lapack_wrapper::print_matrix( 1, N, e, 1 );

  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, N, -1, A, LDA, x, 1, 1, b, 1
  );

  cout
    << "res^T    = "
    << lapack_wrapper::print_matrix( 1, M, b, 1 );

  valueType res = lapack_wrapper::nrm2( M, b, 1 );
  UTILS_ASSERT( res < 1e-6, "test failed! res = {}\n", res );

  msg.green(
    fmt::format(
      "residual = {}\n||res||_2 = {}\ndone test10\n",
      lapack_wrapper::print_matrix( 1, M, b, 1 ), res
    )
  );
}

static
void
test11() {
  lapack_wrapper::PINV<valueType> pinv;

  integer const M   = 6;
  integer const N   = 4;
  integer const LDA = 6;
  valueType A[] = {
    2,   1,  2,  1, 1, 1,
    2,   1,  2,  1, 1, 1,
    20, 10, 20, 10, 1, 1,
     1,  2,  3,  4, 5, 6,
  };

  valueType rhs[N], b[N], b2[2*N], x[M], x2[2*M];
  valueType const xe[M] = {1,2,3,4,5,6};
  lapack_wrapper::gemv(
    lapack_wrapper::TRANSPOSE, M, N, 1, A, LDA, xe, 1, 0, rhs, 1
  );
  lapack_wrapper::copy( N, rhs, 1, b, 1 );
  lapack_wrapper::copy( N, rhs, 1, b2, 1 );
  lapack_wrapper::copy( N, rhs, 1, b2+N, 1 );
  msg.green( fmt::format(
    "\n\n\nTest11:\n\nInitial A{}\nb2{}\n",
    lapack_wrapper::print_matrix( M, N, A, M ),
    lapack_wrapper::print_matrix( N, 2, b2, N )
  ) );

  cout << "\n\nDo PINV factorization of A\n";
  pinv.factorize( M, N, A, LDA );

  cout << "\nPINV solution of A x = b\n";
  pinv.t_mult_inv( b, 1, x, 1 );
  pinv.t_mult_inv( 2, b2, N, x2, M );

  //pinv.t_solve( 1, x, M );
  cout
    << "x^T      = "
    << lapack_wrapper::print_matrix( 1, M, x, 1 )
    << "x2^T     = "
    << lapack_wrapper::print_matrix( 1, M, x2, 1 );

  valueType e[M] = { x[0]-xe[0], x[1]-xe[1], x[2]-xe[2], x[3]-xe[3] };

  cout
    << "e^T      = "
    << lapack_wrapper::print_matrix( 1, M, e, 1 );

  lapack_wrapper::gemv(
    lapack_wrapper::TRANSPOSE, M, N, -1, A, LDA, x, 1, 1, b, 1
  );

  cout
    << "res^T    = "
    << lapack_wrapper::print_matrix( 1, N, b, 1 );

  valueType res = lapack_wrapper::nrm2( N, b, 1 );
  UTILS_ASSERT( res < 1e-6, "test failed! res = {}\n", res );

  msg.green(
    fmt::format(
      "residual = {}\n||res||_2 = {}\ndone test11\n",
      lapack_wrapper::print_matrix( 1, N, b, 1 ), res
    )
  );
}

static
void
test12() {
  lapack_wrapper::PINV<valueType> pinv;

  integer const M   = 6;
  integer const N   = 5;
  integer const LDA = M;
  valueType A[] = {
    0.001, 1e-9,  1e-9,  1e-9, 0, 0,
    1,    -1e+9, -1e+9, -1e+9, 0, 0,
    0,     1e-9,     1,  1e-9, 0, 0,
    2,    -1e-9, -1e-9,     1, 0, 0,
    1,        1,     1,     1, 0, 0
  };

  /*

format shortG;

A = [ ...
  0.001, 1e-9,  1e-9,  1e-9, 0, 0; ...
  1,    -1e+9, -1e+9, -1e+9, 0, 0; ...
  0,     1e-9,     1,  1e-9, 0, 0; ...
  2,    -1e-9, -1e-9,     1, 0, 0; ...
  1,        1,     1,     1, 0, 0; ...
].';

disp('pinv(A)');
pinv(A)

disp('pinv(A^T)^T');
pinv(A.').'

disp('pinv(A^T)^T');
pinv(A)-pinv(A.').'

p =  [2,4,3,5,1];
[Q,R]   = qr(A(:,p));
[Q1,R1,p] = qr(R(1:end-1,:).');

R
R1


*/

  valueType MAT[6*6], MAT1[6*6];

  valueType rhs[M], b[M], x[N];
  valueType const xe[M] = {1,2,3,4,5};
  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, N, 1, A, LDA, xe, 1, 0, rhs, 1
  );

  msg.green( fmt::format(
    "\n\n\nTest12:\n\nInitial A\n{}\nb^T{}\n",
    lapack_wrapper::print_matrix( M, N, A, M ),
    lapack_wrapper::print_matrix( 1, M, rhs, 1 )
  ) );

  cout << "\n\nDo PINV factorization of A\n";
  pinv.factorize( M, N, A, LDA );
  //pinv.info( std::cout );

  lapack_wrapper::geid( M, M, MAT, M );

  pinv.mult_inv( M, MAT, M, MAT1, N );
  msg.blue( fmt::format(
    "pinv(A)\n{}",
    lapack_wrapper::print_matrix( N, M, MAT1, N )
  ) );

  lapack_wrapper::geid( N, N, MAT, N );
  pinv.t_mult_inv( N, MAT, N, MAT1, M );
  msg.blue( fmt::format(
    "pinv(A)\n{}",
    lapack_wrapper::print_matrix_transpose( M, N, MAT1, M )
  ) );

  lapack_wrapper::copy( M, rhs, 1, b, 1 );
  pinv.mult_inv( b, 1, x, 1 );
  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, N, -1, A, LDA, x, 1, 1, b, 1
  );
  valueType res = lapack_wrapper::nrm2( M, b, 1 );

  fmt::print(
    "res = {}\n"
    "r = {}"
    "x = {}\n",
    res,
    lapack_wrapper::print_matrix_transpose( M, 1, b, M ),
    lapack_wrapper::print_matrix_transpose( N, 1, x, N )
  );

}


int
main() {

  try {
    /*
    test1();
    test2();
    test3();
    test4();
    test5();
    test6();
    test7();
    test8();
    test9();
    test10();
    test11();
    */
    test12();
  } catch ( exception const & exc ) {
    msg.error( exc.what() );
  } catch ( ... ) {
    msg.error("Errore Sconosciuto!\n");
  }

  cout << "\n\nAll done!\n";

  return 0;
}
