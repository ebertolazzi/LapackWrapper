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

static Utils::Console msg(&std::cout);

int
main() {

  try {

    lapack_wrapper::doublereal s1[] = {1,2,3};
    lapack_wrapper::doublereal y1[] = {11./3.,-6./3.,13./3.};

    lapack_wrapper::doublereal s2[] = {1,0,1};
    lapack_wrapper::doublereal y2[] = {3,-2,3};

    lapack_wrapper::doublereal s3[] = {0,-1,1};
    lapack_wrapper::doublereal y3[] = {7/3.,-6/3.,8/3.};

    lapack_wrapper::BFGS<lapack_wrapper::doublereal> bfgs;

    lapack_wrapper::doublereal epsi = 1e-8;

    bfgs.allocate(3);
    bfgs.init();
    fmt::print( "\n\n{}\n", bfgs.to_string() );

    bfgs.update( y1, s1, epsi );
    fmt::print( "\n\n{}\n", bfgs.to_string() );

    bfgs.update( y2, s2, epsi );
    fmt::print( "\n\n{}\n", bfgs.to_string() );

    bfgs.update( y3, s3, epsi );
    fmt::print( "\n\n{}\n", bfgs.to_string() );

    for ( int i = 0; i < 90; ++i ) {
      bfgs.update( y1, s1, epsi );
      bfgs.update( y2, s2, epsi );
      bfgs.update( y3, s3, epsi );
    }
    fmt::print( "\n\n{}\n", bfgs.to_string() );

  } catch ( exception const & exc ) {
    msg.error( exc.what() );
  } catch ( ... ) {
    msg.error("Errore Sconosciuto!\n");
  }

  cout << "All done!\n";
  return 0;
}
