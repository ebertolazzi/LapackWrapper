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
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#endif

using namespace std;
using real_type = double;
using lapack_wrapper::integer;
using lapack_wrapper::Transposition;

static Utils::Console msg(&std::cout);

static
void
test1() {

  lapack_wrapper::QR<real_type> qr;

  integer const M   = 3;
  integer const N   = 5;
  integer const LDA = 3;
  real_type A[] = {
    0.001,      2,     3,
    0.001,  0.001,     0,
    0,      0.001,     0,
    0.001,     -1,     0,
    0.000001,   5,     3
  };

  lapack_wrapper::MatrixWrapper<real_type> Amat( A, M, N, LDA );

  msg.green(
    "\n\n\nTest1:\n\nInitial A\n" +
    lapack_wrapper::print_matrix( Amat )
  );

  msg.green( "\nDo QR factorization of A^T\n" );
  qr.t_factorize( "qr", Amat );

  real_type R[M*M];
  qr.getR( R, M );
  fmt::print( "\nR=\n{}", lapack_wrapper::print_matrix( M, M, R, M ) );

  real_type rhs[M], b[M];
  real_type x[N] = {1,2,3,4,5};
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
  real_type res = lapack_wrapper::nrm2( M, b, 1 );
  UTILS_ASSERT( res < 1e-6, "test failed! res = {}\n", res );

  fmt::print(
    "residual = {}\n||res||_2 = {}\ndone test1\n",
    lapack_wrapper::print_matrix( 1, M, b, 1 ), res
  );
  msg.green( "\n\ndone test1\n" );

}



static
void
test2() {
  lapack_wrapper::QRP<real_type> qr;
  integer const M   = 3;
  integer const N   = 5;
  integer const LDA = 3;
  real_type A[] = {
    0.001,      2,     3,
    0.001,  0.001,     0,
    0,      0.001,     0,
    0.001,     -1,     0,
    0.000001,   5,     3
  };

  lapack_wrapper::MatrixWrapper<real_type> Amat( A, M, N, LDA );

  msg.green( "\n\n\nTest2:\n\nInitial A\n"+ lapack_wrapper::print_matrix( Amat ) );

  fmt::print("\nDo QR factorization of A^T\n");
  qr.t_factorize( "qrp", Amat );

  real_type R[M*M];
  qr.getR( R, M );
  fmt::print( "\nR=\n{}", lapack_wrapper::print_matrix( M, M, R, M ) );

  real_type rhs[M], b[M];
  real_type x[N] = {1,2,3,4,5};
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
  real_type res = lapack_wrapper::nrm2( M, b, 1 );
  UTILS_ASSERT( res < 1e-6, "test failed! res = {}\n", res );

  fmt::print(
    "residual  = {}"
    "||res||_2 = {}\n"
    "done test2\n",
    lapack_wrapper::print_matrix( 1, M, b, 1 ), res
  );
  msg.green( "\n\ndone test2\n" );

}

static
void
test3() {
  lapack_wrapper::QRP<real_type> qr;
  integer const M   = 5;
  integer const N   = 5;
  integer const LDA = 5;
  real_type A[] = {
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

  real_type R[M*M];
  qr.getR( R, M );
  fmt::print( "\nR=\n{}", lapack_wrapper::print_matrix( M, M, R, M ) );

  real_type rhs[M], b[M];
  real_type x[N] = {1,2,3,4,5};
  lapack_wrapper::gemv(
    Transposition::NO, M, N, 1, A, LDA, x, 1, 0, rhs, 1
  );
  lapack_wrapper::gemv(
    Transposition::NO, M, N, 1, A, LDA, x, 1, 0, b, 1
  );

  fmt::print(
    "\nLS solution of A x = b\n\n"
    "b^T =      {}\n",
    lapack_wrapper::print_matrix( 1, M, b, 1 )
  );

  qr.inv_permute( rhs ); // da aggiungere!
  qr.invRt_mul( rhs, 1 );
  lapack_wrapper::copy( 3, rhs, 1, x, 1 );
  lapack_wrapper::zero( 2, x+3, 1 );
  qr.Q_mul( x );

  lapack_wrapper::gemv(
    Transposition::NO, M, N, -1, A, LDA, x, 1, 1, b, 1
  );

  fmt::print(
    "x^T =      {}\n"
    "residual = {}\n",
    lapack_wrapper::print_matrix( 1, M, x, 1 ),
    lapack_wrapper::print_matrix( 1, M, b, 1 )
  );

  real_type res = lapack_wrapper::nrm2( M, b, 1 );
  UTILS_ASSERT( res < 1e-6, "test failed! res = {}\n", res );

  fmt::print( "||res||_2 = {}\n", res );
  msg.green( "\n\ndone test3\n" );

}

#define TEST4(NAME,F) \
  fmt::print("\n\nDo {} factorization of A\n", NAME ); \
  F.factorize( NAME, M, M, A, LDA ); \
  \
  fmt::print("{} solution of A x = b\n", NAME); \
  lapack_wrapper::copy( M, rhs, 1, x, 1 ); \
  lapack_wrapper::copy( M, rhs, 1, b, 1 ); \
  /* L.solve( x ); */ \
  F.solve( 1, x, M); \
  fmt::print("x^T      = {}", lapack_wrapper::print_matrix( 1, M, x, 1 ) ); \
  lapack_wrapper::gemv( Transposition::NO, M, M, -1, A, LDA, x, 1, 1, b, 1 ); \
  fmt::print("residual = {}", lapack_wrapper::print_matrix( 1, M, b, 1 ) ); \
  res = lapack_wrapper::nrm2( M, b, 1 ); \
  fmt::print( "||res||_2 = {}\n", res ); \
  UTILS_ASSERT0( res < 1e-6, "test failed!\n" );

static
void
test4() {
  lapack_wrapper::LU<real_type>   lu;
  lapack_wrapper::LUPQ<real_type> lupq;
  lapack_wrapper::QR<real_type>   qr;
  lapack_wrapper::QRP<real_type>  qrp;
  lapack_wrapper::SVD<real_type>  svd;
  lapack_wrapper::LSS<real_type>  lss;
  lapack_wrapper::LSY<real_type>  lsy;

  integer const M   = 5;
  integer const LDA = 5;
  real_type A[] = {
    0.001,      2,     3,       2,      3,
    0.001,  0.001,     0,   0.001,  1e-10,
    0,      0.001,     0,   0.001,  1e-12,
    0.001,     -1, 1e-6+1,     -1, -1e-12,
    0.000001,   5,     3,  1e-6+5,    3+1
  };

  real_type rhs[M], b[M], res;
  real_type x[M] = {1,2,3,4,5};
  lapack_wrapper::gemv(
    Transposition::NO, M, M, 1, A, LDA, x, 1, 0, rhs, 1
  );
  lapack_wrapper::gemv(
    Transposition::NO, M, M, 1, A, LDA, x, 1, 0, b, 1
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

  msg.green( "\n\ndone test4\n" );

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
  lapack_wrapper::gemv( Transposition::YES, M, M, -1, A, LDA, x, 1, 1, b, 1 ); \
  cout << "residual = " \
       << lapack_wrapper::print_matrix( 1, M, b, 1 ); \
  res = lapack_wrapper::nrm2( M, b, 1 ); \
  cout << "||res||_2 = " << res << '\n'; \
  UTILS_ASSERT0( res < 1e-6, "test failed!\n" );


static
void
test5() {
  lapack_wrapper::LU<real_type>   lu;
  lapack_wrapper::LUPQ<real_type> lupq;
  lapack_wrapper::QR<real_type>   qr;
  lapack_wrapper::QRP<real_type>  qrp;
  lapack_wrapper::SVD<real_type>  svd;
  lapack_wrapper::LSS<real_type>  lss;
  lapack_wrapper::LSY<real_type>  lsy;

  integer const M   = 5;
  integer const LDA = 5;
  real_type A[] = {
    0.001,      2,     3,       2,      3,
    0.001,  0.001,     0,   0.001,  1e-10,
    0,      0.001,     0,   0.001,  1e-12,
    0.001,     -1, 1e-6+1,     -1, -1e-12,
    0.000001,   5,     3,  1e-6+5,     3+1
  };

  real_type rhs[M], b[M], res;
  real_type x[M] = {1,2,3,4,5};
  lapack_wrapper::gemv(
    Transposition::YES, M, M, 1, A, LDA, x, 1, 0, rhs, 1
  );
  lapack_wrapper::gemv(
    Transposition::YES, M, M, 1, A, LDA, x, 1, 0, b, 1
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

  msg.green( "\n\ndone test5\n" );

}


static
void
test6() {
  lapack_wrapper::TridiagonalLU<real_type> lu;
  lapack_wrapper::TridiagonalQR<real_type> qr;

  integer const N = 5;
  real_type const D[] = { 1, 1, 2, -0.1, 0.1 };
  real_type const L[] = { -0.1, -1, -2, -0.1 };
  real_type const U[] = { -1, -10, -2, 0.1 };

  real_type rhs[N], b[N];
  real_type x[N] = {1,2,3,4,5};
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

  qr.axpy( N, -1.0, L, D, U, x, 1.0, b );

  fmt::print(
    "x^T      = {}\n"
    "residual = {}\n",
    lapack_wrapper::print_matrix( 1, N, x, 1 ),
    lapack_wrapper::print_matrix( 1, N, b, 1 )
  );

  real_type res = lapack_wrapper::nrm2( N, b, 1 );
  UTILS_ASSERT( res < 1e-6, "test failed! res = {}\n", res );

  fmt::print( "||res||_2 = {}\n", res );
  msg.green( "\n\ndone test6\n" );

}

static
void
test7() {
  lapack_wrapper::QRP<real_type> qrp;

  integer const M   = 5;
  integer const LDA = 5;
  real_type A[] = {
    0.001,      2,     3,       2,      3,
    0.001,  0.001,     0,   0.001,  1e-10,
    0,      0.001,     0,   0.001,  1e-12,
    0.001,     -1, 1e-6+1,     -1, -1e-12,
    0.000001,   5,     3,  1e-6+5,    3+1
  };

  real_type rhs[M], b[M];
  real_type x[M] = {1,2,3,4,5};
  lapack_wrapper::gemv(
    Transposition::YES, M, M, 1, A, LDA, x, 1, 0, rhs, 1
  );
  lapack_wrapper::gemv(
    Transposition::YES, M, M, 1, A, LDA, x, 1, 0, b, 1
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

  lapack_wrapper::gemv(
    Transposition::YES, M, M, -1, A, LDA, x, 1, 1, b, 1
  );

  fmt::print(
    "x^T      = {}\n"
    "residual = {}\n",
    lapack_wrapper::print_matrix( 1, M, x, 1 ),
    lapack_wrapper::print_matrix( 1, M, b, 1 )
  );

  real_type res = lapack_wrapper::nrm2( M, b, 1 );
  UTILS_ASSERT( res < 1e-6, "test failed! res = {}\n", res );

  fmt::print( "||res||_2 = {}\n", res );
  msg.green( "\n\ndone test7\n" );

}

static
void
test8() {
  lapack_wrapper::LSC<real_type> lsc;

  integer const M   = 5;
  integer const LDA = 5;
  real_type A[] = {
    0.001,      2,     3,       2,      3,
    0.001,  0.001,     0,   0.001,  1e-10,
    0,      0.001,     0,   0.001,  1e-12,
    0.001,     -1, 1e-6+1,     -1, -1e-12,
    0.000001,   5,     3,  1e-6+5,    3+1
  };

  real_type rhs[M], b[M];
  real_type x[M] = {1,2,3,4,5};
  //bool      rselect[M] = { true, true, true, true, true };
  //bool      rselect[M] = { false, false, false, false, false };
  bool      rselect[M] = { true, false, true, true, true };
  lapack_wrapper::gemv(
    Transposition::NO, M, M, 1, A, LDA, x, 1, 0, rhs, 1
  );
  lapack_wrapper::gemv(
    Transposition::NO, M, M, 1, A, LDA, x, 1, 0, b, 1
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

  lapack_wrapper::gemv(
    Transposition::NO, M, M, -1, A, LDA, x, 1, 1, b, 1
  );

  fmt::print(
    "x^T      = {}\n"
    "residual = {}\n",
    lapack_wrapper::print_matrix( 1, M, x, 1 ),
    lapack_wrapper::print_matrix( 1, M, b, 1 )
  );

  real_type res = lapack_wrapper::nrm2( M, b, 1 );
  UTILS_ASSERT( res < 1e-6, "test failed! res = {}\n", res );

  fmt::print( "||res||_2 = {}\n", res );

  msg.green( "\n\ndone test8\n" );
}


static
void
test9() {
  lapack_wrapper::PINV<real_type> pinv;

  integer const M   = 5;
  integer const LDA = 5;
  real_type A[] = {
    0.001,    1,     0,      2,  1,
    1e-9, -1e+9,  1e-9,  -1e-9,  1,
    1e-9, -1e+9,     1,  -1e-9,  1,
    1e-9, -1e+9,  1e-9,      1,  1,
    //1e-9, -1e+9,  1e-9,  -1e-9,  1
    0, 0, 0, 0, 0
  };

  real_type rhs[M], b[M], x[M];
  real_type const xe[M] = {1,2,3,4,5};
  lapack_wrapper::gemv(
    Transposition::NO, M, M, 1, A, LDA, xe, 1, 0, rhs, 1
  );
  lapack_wrapper::copy( M, rhs, 1, b, 1 );

  msg.green( fmt::format(
    "\n\n\nTest9:\n\n"
    "Initial A = {}\n"
    "b^T       = {}\n",
    lapack_wrapper::print_matrix( M, M, A, M ),
    lapack_wrapper::print_matrix( 1, M, b, 1 )
  ) );

  cout << "\n\nDo PINV factorization of A\n";
  pinv.factorize( M, M, A, LDA );

  cout << "\nLSC solution of A x = b\n";
  lapack_wrapper::copy( M, b, 1, x, 1 );
  pinv.solve( x );
  //pinv.t_solve( 1, x, M );

  real_type e[M] = { x[0]-xe[0],x[1]-xe[1],x[2]-xe[2],x[3]-xe[3],x[4]-xe[4]};

  lapack_wrapper::gemv(
    Transposition::NO, M, M, -1, A, LDA, x, 1, 1, b, 1
  );

  fmt::print(
    "x^T      = {}\n"
    "e^T      = {}\n"
    "residual = {}\n",
    lapack_wrapper::print_matrix( 1, M, x, 1 ),
    lapack_wrapper::print_matrix( 1, M, e, 1 ),
    lapack_wrapper::print_matrix( 1, M, b, 1 )
  );

  real_type res = lapack_wrapper::nrm2( M, b, 1 );
  UTILS_WARNING( res < 1e-6, "\n\n\n\ntest failed! res = {}\n\n\n\n", res );

  fmt::print( "||res||_2 = {}\n", res );

  cout << "\n\nDo PINV factorization of A\n";
  pinv.factorize( M, M-1, A, LDA );
  pinv.mult_inv( rhs, 1, x, 1 );

  lapack_wrapper::copy( M, rhs, 1, b, 1 );
  lapack_wrapper::gemv(
    Transposition::NO, M, M-1, -1, A, LDA, x, 1, 1, b, 1
  );
  res = lapack_wrapper::nrm2( M, b, 1 );

  fmt::print(
    "x^T       = {}"
    "residual  = {}\n"
    "||res||_2 = {}\n",
    lapack_wrapper::print_matrix( 1, M-1, x, 1 ),
    lapack_wrapper::print_matrix( 1, M, b, 1 ),
    res
  );

  msg.green( "\n\ndone test9\n" );

}

static
void
test10() {
  lapack_wrapper::PINV<real_type> pinv;

  integer const M   = 6;
  integer const N   = 4;
  integer const LDA = 6;
  real_type A[] = {
    2,   1,  2,  1, 1, 1,
    2,   1,  2,  1, 1, 1,
    20, 10, 20, 10, 1, 1,
     1,  2,  3,  4, 5, 6,
  };

  real_type rhs[M], b[M], b2[2*M], x[N], x2[2*N];
  real_type const xe[N] = {1,2,3,4};
  lapack_wrapper::gemv(
    Transposition::NO, M, N, 1, A, LDA, xe, 1, 0, rhs, 1
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

  real_type e[N] = { x[0]-xe[0], x[1]-xe[1], x[2]-xe[2] };

  lapack_wrapper::gemv(
    Transposition::NO, M, N, -1, A, LDA, x, 1, 1, b, 1
  );

  //pinv.t_solve( 1, x, M );
  fmt::print(
    "x^T      = {}\n"
    "x2^T     = {}\n"
    "e^T      = {}\n"
    "res^T    = {}\n",
    lapack_wrapper::print_matrix( 1, N, x, 1 ),
    lapack_wrapper::print_matrix( 1, N, x2, 1 ),
    lapack_wrapper::print_matrix( 1, N, e, 1 ),
    lapack_wrapper::print_matrix( 1, M, b, 1 )
  );

  real_type res = lapack_wrapper::nrm2( M, b, 1 );
  UTILS_ASSERT( res < 1e-6, "test failed! res = {}\n", res );

  fmt::print( "||res||_2 = {}\n", res );

  msg.green( "\n\ndone test10\n" );
}

static
void
test11() {
  lapack_wrapper::PINV<real_type> pinv;

  integer const M   = 6;
  integer const N   = 4;
  integer const LDA = 6;
  real_type A[] = {
    2,   1,  2,  1, 1, 1,
    2,   1,  2,  1, 1, 1,
    20, 10, 20, 10, 1, 1,
     1,  2,  3,  4, 5, 6,
  };

  real_type rhs[N], b[N], b2[2*N], x[M], x2[2*M];
  real_type const xe[M] = {1,2,3,4,5,6};
  lapack_wrapper::gemv( Transposition::YES, M, N, 1, A, LDA, xe, 1, 0, rhs, 1 );
  lapack_wrapper::copy( N, rhs, 1, b, 1 );
  lapack_wrapper::copy( N, rhs, 1, b2, 1 );
  lapack_wrapper::copy( N, rhs, 1, b2+N, 1 );
  msg.green( fmt::format(
    "\n\n\nTest11:\n\n"
    "Initial A = {}\n"
    "b2        = {}\n",
    lapack_wrapper::print_matrix( M, N, A, M ),
    lapack_wrapper::print_matrix( N, 2, b2, N )
  ) );

  cout << "\n\nDo PINV factorization of A\n";
  pinv.factorize( M, N, A, LDA );

  cout << "\nPINV solution of A x = b\n";
  pinv.t_mult_inv( b, 1, x, 1 );
  pinv.t_mult_inv( 2, b2, N, x2, M );

  real_type e[M] = { x[0]-xe[0], x[1]-xe[1], x[2]-xe[2], x[3]-xe[3] };

  lapack_wrapper::gemv( Transposition::YES, M, N, -1, A, LDA, x, 1, 1, b, 1 );

  //pinv.t_solve( 1, x, M );
  fmt::print(
    "x^T      = {}\n"
    "x2^T     = {}\n"
    "e^T      = {}\n"
    "res^T    = {}\n",
    lapack_wrapper::print_matrix( 1, M, x, 1 ),
    lapack_wrapper::print_matrix( 1, M, x2, 1 ),
    lapack_wrapper::print_matrix( 1, M, e, 1 ),
    lapack_wrapper::print_matrix( 1, N, b, 1 )
  );

  real_type res = lapack_wrapper::nrm2( N, b, 1 );
  UTILS_WARNING(
    res < 1e-6,
    "\n\n\n\n\ntest failed! res = {}\n\n\n\n", res
  );

  fmt::print( "||res||_2 = {}\n", res );

  msg.green( "\n\ndone test11\n" );
}

static
void
test12() {
  lapack_wrapper::PINV<real_type> pinv;

  integer const M   = 8;
  integer const N   = 5;
  integer const LDA = M;
  real_type A[] = {
    0.001, 1e-9,  1e-9,  1e-9, 1, 1, 0, 0,
    1,    -1e+9, -1e+9, -1e+9, 1, 1, 0, 0,
    0,     1e-9,     1,  1e-9, 1, 1, 0, 0,
    2,    -1e-9, -1e-9,     1, 1, 1, 0, 0,
    1,        1,     1,     1, 1, 1, 0, 0
  };

  /*

format shortG;

A = [ ...
    0.001, 1e-9,  1e-9,  1e-9, 1, 1, 0, 0; ...
    1,    -1e+9, -1e+9, -1e+9, 1, 1, 0, 0; ...
    0,     1e-9,     1,  1e-9, 1, 1, 0, 0; ...
    2,    -1e-9, -1e-9,     1, 1, 1, 0, 0; ...
    1,        1,     1,     1, 1, 1, 0, 0; ...
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

  real_type MAT[M*M], MAT1[M*M], MM[M*M];

  real_type rhs[N];
  real_type const xe[M] = {1,2,3,4,5,6,7,8};
  lapack_wrapper::gemv( Transposition::YES, M, N, 1, A, LDA, xe, 1, 0, rhs, 1 );

  msg.green( fmt::format(
    "\n\n\nTest12:\n\nInitial A\n{}\n(A^T * x)^T = {}\n",
    lapack_wrapper::print_matrix( M, N, A, M ),
    lapack_wrapper::print_matrix( 1, N, rhs, 1 )
  ) );

  cout << "\n\nDo PINV factorization of A\n";
  pinv.factorize( M, N, A, LDA );
  //pinv.info( std::cout );

  lapack_wrapper::geid( M, M, MAT, M );

  pinv.mult_inv( M, MAT, M, MAT1, N );
  //for ( integer i = 0; i < M; ++i )
  //  pinv.mult_inv( MAT+M*i, 1, MAT1+N*i, 1 );

  msg.blue( fmt::format(
    "\n\npinv(A)\n{}\n\n",
    lapack_wrapper::print_matrix( N, M, MAT1, N )
  ) );

  lapack_wrapper::gemm(
    Transposition::NO,
    Transposition::NO,
    M, M, N,
    1.0, A, LDA,
    MAT1, N,
    0.0, MM, M
  );
  msg.blue( fmt::format(
    "\n\nA*pinv(A)\n{}\n\n",
    lapack_wrapper::print_matrix( M, M, MM, M )
  ) );

  lapack_wrapper::gemm(
    Transposition::NO,
    Transposition::NO,
    N, N, M,
    1.0, MAT1, N,
    A, LDA,
    0.0, MM, N
  );
  msg.blue( fmt::format(
    "\n\npinv(A)*A\n{}\n\n",
    lapack_wrapper::print_matrix( N, N, MM, N )
  ) );

  lapack_wrapper::geid( N, N, MAT, N );

  pinv.t_mult_inv( N, MAT, N, MAT1, M );
  //for ( integer i = 0; i < N; ++i )
  //  pinv.t_mult_inv( MAT+N*i, 1, MAT1+M*i, 1 );

  msg.blue( fmt::format(
    "pinv(A)^T\n{}\n\n",
    lapack_wrapper::print_matrix( M, N, MAT1, M )
  ) );

  lapack_wrapper::gemm(
    Transposition::YES,
    Transposition::NO,
    N, N, M,
    1.0, A, LDA,
    MAT1, M,
    0.0, MM, N
  );
  msg.blue( fmt::format(
    "\n\nA^T*pinv(A)^T\n{}\n\n",
    lapack_wrapper::print_matrix( N, N, MM, N )
  ) );

  lapack_wrapper::gemm(
    Transposition::NO,
    Transposition::YES,
    M, M, N,
    1.0,
    MAT1, M,
    A, LDA,
    0.0, MM, M
  );
  msg.blue( fmt::format(
    "\n\npinv(A)^T*A^T\n{}\n\n",
    lapack_wrapper::print_matrix( M, M, MM, M )
  ) );

  msg.green( "\n\ndone test12\n" );

}


int
main() {

  try {
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
    test12();
  } catch ( exception const & exc ) {
    msg.error( exc.what() );
  } catch ( ... ) {
    msg.error("Errore Sconosciuto!\n");
  }

  cout << "\n\nAll done!\n";

  return 0;
}
