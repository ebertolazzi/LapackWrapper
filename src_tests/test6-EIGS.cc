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

#include <lapack_wrapper/lapack_wrapper.hh>
#include <lapack_wrapper/lapack_wrapper++.hh>

#include <iostream>

using namespace std;

static Utils::Console msg(&std::cout);

using real_type = lapack_wrapper::doublereal;
using lapack_wrapper::Transposition;

static
void
test1() {
  // solution
  // 0, 0, 7/2+(1/2)*sqrt(73), 7/2-(1/2)*sqrt(73)])
  // 0., 0., 7.772001872, -.772001872
  real_type A[] = { 1, 2, 3, 4,
                    4, 4, 4, 4,
                    1, 2, 1, 2,
                    0, 0, 1, 1 };
  lapack_wrapper::Eigenvalues<real_type> E;
  E.setup( 4, A, 4 );
  vector<complex<real_type> > e;
  E.getEigenvalues( e );

  vector<complex<real_type> >::const_iterator it;
  for ( it = e.begin(); it != e.end(); ++it )
    cout << *it << "\n";
}

static
void
test2() {
  // solution
  // 0, 4, -2/3-I*sqrt(11)*(1/3), -2/3+I*sqrt(11)*(1/3)
  real_type A[] = { 1, 2, 3, 4,
                    4, 4, 4, 4,
                    1, 2, 1, 2,
                    0, 0, 1, 1 };
  real_type B[] = { 3, 2, 1, 1,
                    1, 1, 1, 1,
                    4, 3, 2, 1,
                    1,-1, 0, 1 };
  lapack_wrapper::GeneralizedEigenvalues<real_type> E;
  E.setup( 4, A, 4, B, 4 );
  vector<complex<real_type> > e;
  E.getEigenvalues( e );

  vector<complex<real_type> >::const_iterator it;
  for ( it = e.begin(); it != e.end(); ++it )
    cout << *it << "\n";
}

static
void
test3() {
  /*
  A = [ 1, 2, 3, 4; 4, 4, 4, 4; 1, 2, 1, 2; 0, 0, 1, 1].';
  B = [ 3, 2, 1, 1; 1, 1, 1, 1; 4, 3, 2, 1; 1,-1, 0, 1].';
  [U,V,X,C,S] = gsvd(A,B)
  %[U,V,X,C,S] = gsvd(A,B,0)
  */
  real_type A_data[] = { 1, 2, 3, 4,
                         4, 4, 4, 4,
                         1, 2, 1, 2,
                         0, 0, 1, 1 };
  real_type B_data[] = { 3, 2, 1, 1,
                         1, 1, 1, 1,
                         4, 3, 2, 1,
                         1,-1, 0, 1 };
  real_type TMP_data[16], TMP1_data[16];
  lapack_wrapper::MatrixWrapper<real_type> A(A_data,4,4,4);
  lapack_wrapper::MatrixWrapper<real_type> B(B_data,4,4,4);
  lapack_wrapper::MatrixWrapper<real_type> TMP(TMP_data,4,4,4);
  lapack_wrapper::MatrixWrapper<real_type> TMP1(TMP1_data,4,4,4);

  lapack_wrapper::GeneralizedSVD<real_type> E;
  E.setup( A, B );
  fmt::print( "{}\n", E.info( 1e-10 ) );

  lapack_wrapper::gemm(
    1.0,
    Transposition::YES, E.getU(),
    Transposition::NO,  A,
    0.0, TMP
  );
  lapack_wrapper::gemm(
    1.0,
    Transposition::NO, TMP,
    Transposition::NO, E.getQ(),
    0.0, TMP1
  );
  fmt::print( "TMP1\n{}", TMP1.to_string(1e-12) );

  lapack_wrapper::gemm(
    1.0,
    Transposition::YES, E.getV(),
    Transposition::NO,  B,
    0.0, TMP
  );
  lapack_wrapper::gemm(
    1.0,
    Transposition::NO, TMP,
    Transposition::NO, E.getQ(),
    0.0, TMP1
  );
  fmt::print( "TMP1\n{}", TMP1.to_string(1e-12) );
}

static
void
test4() {
  /*
  */
  real_type A_data[] = { 1, 2, 3, 4, 5,
                         6, 7, 8, 9, 10,
                         11,12,13,14,15 };
  real_type B_data[] = { 8, 3, 4,
                         1, 5, 9,
                         6, 7, 2 };
  real_type TMP_data[100], TMP1_data[100], TMP2_data[100], TMP3_data[100];
  lapack_wrapper::MatrixWrapper<real_type> A(A_data,5,3,5);
  lapack_wrapper::MatrixWrapper<real_type> B(B_data,3,3,3);
  lapack_wrapper::MatrixWrapper<real_type> TMP(TMP_data,5,3,5);
  lapack_wrapper::MatrixWrapper<real_type> TMP1(TMP1_data,5,3,5);
  lapack_wrapper::MatrixWrapper<real_type> TMP2(TMP2_data,3,3,3);
  lapack_wrapper::MatrixWrapper<real_type> TMP3(TMP3_data,3,3,3);

  lapack_wrapper::GeneralizedSVD<real_type> E;
  E.setup( A, B );
  fmt::print( "{}\n", E.info(1e-12) );

  lapack_wrapper::gemm(
    1.0,
    Transposition::YES, E.getU(),
    Transposition::NO,  A,
    0.0, TMP
  );
  lapack_wrapper::gemm(
    1.0,
    Transposition::NO, TMP,
    Transposition::NO, E.getQ(),
    0.0, TMP1
  );
  fmt::print( "TMP1\n{}", TMP1.to_string(1e-12) );

  lapack_wrapper::gemm(
    1.0,
    Transposition::YES, E.getV(),
    Transposition::NO,  B,
    0.0, TMP2
  );
  lapack_wrapper::gemm(
    1.0,
    Transposition::NO, TMP2,
    Transposition::NO, E.getQ(),
    0.0, TMP3
  );
  fmt::print( "TMP3\n{}", TMP3.to_string(1e-12) );
}

static
void
test5() {
  // solution
  // 0, 0, 7/2+(1/2)*sqrt(73), 7/2-(1/2)*sqrt(73)])
  real_type A_data[] = {
    2, 0, 2, 0,
    0, 2, 4, 0,
    2, 4, 4, 1,
    0, 0, 0, 4
  };
  lapack_wrapper::MatrixWrapper<real_type> A(A_data,4,4,4);
  lapack_wrapper::Eigenvectors<real_type> E;
  E.setup( 4, A_data, 4 );
  vector<complex<real_type> > e;
  E.getEigenvalues( e );

  for ( auto const & it : e )
    fmt::print( "{}\n", it );

  vector<vector< lapack_wrapper::Eigenvectors<real_type>::complex_type > > vecs;
  E.getLeftEigenvector( vecs );
  for ( size_t n = 0; n < 4; ++n ) {
    fmt::print( "vL[{}] = ", n );
    for ( size_t i = 0; i < 4; ++i )
      fmt::print( " {}", vecs[n][i] );
    fmt::print("\n");
  }
  E.getRightEigenvector( vecs );
  for ( size_t n = 0; n < 4; ++n ) {
    fmt::print( "vR[{}] = ", n );
    for ( size_t i = 0; i < 4; ++i )
      fmt::print( "{} ", vecs[n][i] );
    fmt::print("\n");
  }
}

static
void
test6() {
  // solution
  // 0, 0, 7/2+(1/2)*sqrt(73), 7/2-(1/2)*sqrt(73)])
  real_type A_data[] = { 1, 2, 3, 4,
                         4, 4, 4, 4,
                         1, 2, 1, 2,
                         0, 0, 1, 1 };
  real_type B_data[] = { 3, 2, 1, 1,
                         1, 1, 1, 1,
                         4, 3, 2, 1,
                         1,-1, 0, 1 };
  lapack_wrapper::MatrixWrapper<real_type> A(A_data,4,4,4);
  lapack_wrapper::MatrixWrapper<real_type> B(B_data,4,4,4);
  lapack_wrapper::GeneralizedEigenvectors<real_type> E;
  E.setup( 4, A_data, 4, B_data, 4 );
  vector<complex<real_type> > e;
  E.getEigenvalues( e );

  for ( auto const & it : e )
    fmt::print( "{}\n", it );

  vector<vector< lapack_wrapper::GeneralizedEigenvectors<real_type>::complex_type > > vecs;
  E.getLeftEigenvector( vecs );
  for ( size_t n = 0; n < 4; ++n ) {
    fmt::print( "vL[{}] = ", n );
    for ( size_t i = 0; i < 4; ++i )
      fmt::print( "{} ", vecs[n][i] );
    fmt::print("\n");
  }
  E.getRightEigenvector( vecs );
  for ( size_t n = 0; n < 4; ++n ) {
    fmt::print( "vR[{}] = ", n );
    for ( size_t i = 0; i < 4; ++i )
      fmt::print( "{} ", vecs[n][i] );
    fmt::print("\n");
  }
}

int
main() {
  try {
    fmt::print( "test1\n" );
    test1();
    fmt::print( "\n\ntest2\n" );
    test2();
    fmt::print( "\n\ntest3\n" );
    test3();
    fmt::print( "\n\ntest4\n" );
    test4();
    fmt::print( "\n\ntest5\n" );
    test5();
    fmt::print( "\n\ntest6\n" );
    test6();
    fmt::print( "\nAll done!\n" );
  } catch ( exception const & exc ) {
    msg.error( exc.what() );
  } catch ( ... ) {
    msg.error("Errore Sconosciuto!\n");
  }
  return 0;
}
