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
#include <lapack_wrapper/TicToc.hh>

#include <iostream>

using namespace std;


int
main() {

  lapack_wrapper::doublereal s1[] = {1,2,3};
  lapack_wrapper::doublereal y1[] = {11./3.,-6./3.,13./3.};

  lapack_wrapper::doublereal s2[] = {1,0,1};
  lapack_wrapper::doublereal y2[] = {3,-2,3};

  lapack_wrapper::doublereal s3[] = {0,-1,1};
  lapack_wrapper::doublereal y3[] = {7/3.,-6/3.,8/3.};

  lapack_wrapper::BFGS<lapack_wrapper::doublereal> bfgs;

  bfgs.allocate(3);
  bfgs.init();
  cout << "\n\n";
  bfgs.print( cout );

  bfgs.update( y1, s1 );
  cout << "\n\n";
  bfgs.print( cout );

  bfgs.update( y2, s2 );
  cout << "\n\n";
  bfgs.print( cout );

  bfgs.update( y3, s3 );
  cout << "\n\n";
  bfgs.print( cout );


  for ( int i = 0; i < 90; ++i ) {
    bfgs.update( y1, s1 );
    bfgs.update( y2, s2 );
    bfgs.update( y3, s3 );
  }
  cout << "\n\n";
  bfgs.print( cout );

  cout << "All done!\n";
  return 0;
}
