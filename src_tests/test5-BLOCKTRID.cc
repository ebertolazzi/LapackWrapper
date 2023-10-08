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
using lapack_wrapper::integer;
using lapack_wrapper::doublereal;

static Utils::Console msg(&std::cout);

static
void
test1() {

  lapack_wrapper::BlockTridiagonalSymmetic<doublereal> BT;
  integer rBlocks[] = { 0, 2, 5, 7};
  BT.setup( 3, rBlocks );

  doublereal D0[] = { 2, 1, 1, 1};
  doublereal D1[] = { 2, 0, 1,
                      0, 0, 0,
                      1, 0, 2};
  doublereal D2[] = { 2, 0,
                      0, 1};

  doublereal L0[] = { 1, 1, 0,
                      0, 1, 1 };
  doublereal L1[] = { 1, 1, 1,
                      1, 1, 1 };

  doublereal rhs[] = {
    5, 4, 6, 4, 6, 5, 4,
    5, 4, 6, 4, 6, 5, 4
  };

  BT.setD( 0, D0, 2 );
  BT.setD( 1, D1, 3 );
  BT.setD( 2, D2, 2 );
  BT.setL( 0, L0, 3 );
  BT.setL( 1, L1, 2 );

  BT.factorize( "BT" );

  //BT.solve( rhs );
  BT.solve( 2, rhs, 7 );

  for ( integer k = 0; k < rBlocks[3]; ++k )
    fmt::print("x[{}] = {}\n", k, rhs[7+k]);
}

static
void
test2() {

  lapack_wrapper::BlockTridiagonalSymmetic<doublereal> BT;
  integer rBlocks[] = { 0, 2, 5, 7};
  BT.setup( 3, rBlocks );

  integer ii[] = {
    1, 2, 1, 2,
    3, 5, 3, 5,
    6, 7,
    3, 4, 4, 5,
    6, 6, 6, 7, 7, 7
  };

  integer jj[] = {
    1, 1, 2, 2,
    3, 3, 5, 5,
    6, 7,
    1, 1, 2, 2,
    3, 4, 5, 3, 4, 5
  };

  doublereal vals[] = {
    2, 1, 1, 1,
    2, 1, 1, 2,
    2, 1,
    1, 1, 1, 1,
    1, 1, 1, 1, 1, 1
  };

  doublereal rhs[] = {
    5, 4, 6, 4, 6, 5, 4,
    5, 4, 6, 4, 6, 5, 4
  };

  BT.zero();
  for ( int k = 0; k < 20; ++k )
    BT(ii[k]-1,jj[k]-1) += vals[k];

  BT.factorize( "BT" );

  //BT.solve( rhs );
  BT.solve( 2, rhs, 7 );

  for ( integer k = 0; k < rBlocks[3]; ++k )
    fmt::print("x[{}] = {}\n", k, rhs[7+k]);
}

int
main() {
  try {
    fmt::print( "test1\n" );
    test1();
    fmt::print( "\n\ntest2\n" );
    test2();
    fmt::print( "All done!\n" );
  } catch ( exception const & exc ) {
    msg.error( exc.what() );
  } catch ( ... ) {
    msg.error("Errore Sconosciuto!\n");
  }
  return 0;
}
