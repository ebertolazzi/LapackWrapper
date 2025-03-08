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

#include <iostream>

using namespace std;

static Utils::Console msg(&std::cout);

int
main() {

  try {

    lapack_wrapper::BandedLU<lapack_wrapper::doublereal> BLU;

    constexpr lapack_wrapper::integer N{10};
    lapack_wrapper::doublereal D []{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
    lapack_wrapper::doublereal L0[]{ 1,   1, 1,  1, 1,  1, 1, 1, 1 };
    lapack_wrapper::doublereal L1[]{ 1,  -1, 1, -1, 1, -1, 1, 1 };
    lapack_wrapper::doublereal L2[]{ 1,-1.4, 1, -1, 1, -1, 1 };
    lapack_wrapper::doublereal U0[]{-1,  -1,-1, -1,-1, -1,-1, -1,-1 };
    lapack_wrapper::doublereal U1[]{ 1,  -1, 1, -1, 1, -1, 1, -1 };
    lapack_wrapper::doublereal rhs[N];

    BLU.setup( N, N, 3, 2 );
    BLU.zero();

    for ( int i{0}; i < N; ++i )
      rhs[i] = BLU(i,i) = D[i];

    for ( int i{0}; i < N-1; ++i ) {
      rhs[i]   += U0[i]; BLU(i,i+1) = U0[i];
      rhs[i+1] += L0[i]; BLU(i+1,i) = L0[i];
    }

    for ( int i{0}; i < N-2; ++i ) {
      rhs[i]   += U1[i]; BLU(i,i+2) = U1[i];
      rhs[i+2] += L1[i]; BLU(i+2,i) = L1[i];
    }

    for ( int i{0}; i < N-3; ++i ) {
      rhs[i+3] += L2[i]; BLU(i+3,i) = L2[i];
    }

    BLU.factorize("BLU");
    BLU.solve(rhs);

    for ( int i{0}; i < N; ++i )
      fmt::print("x[{}] = {}\n", i, rhs[i]);

  } catch ( exception const & exc ) {
    msg.error( exc.what() );
  } catch ( ... ) {
    msg.error("Errore Sconosciuto!\n");
  }

  fmt::print("All done!\n");
  return 0;
}
