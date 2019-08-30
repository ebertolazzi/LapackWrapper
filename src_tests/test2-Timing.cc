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
#include <random>
#include <lapack_wrapper/lapack_wrapper.hh>
#include <lapack_wrapper/lapack_wrapper++.hh>
#include <lapack_wrapper/TicToc.hh>

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wmissing-noreturn"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wdeprecated"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wmissing-noreturn"
#pragma clang diagnostic ignored "-Wfloat-equal"
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wdocumentation-deprecated-sync"
#pragma clang diagnostic ignored "-Wused-but-marked-unused"
#pragma clang diagnostic ignored "-Wdeprecated"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wextra-semi"
#endif

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wc99-extensions"
#pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#pragma clang diagnostic ignored "-Wshadow"
#pragma clang diagnostic ignored "-Wunused-template"
#pragma clang diagnostic ignored "-Wcast-qual"
#endif

using namespace std;
typedef double real_type;

using lapack_wrapper::integer;

static unsigned seed1 = 2;
static std::mt19937 generator(seed1);

static
real_type
rand( real_type xmin, real_type xmax ) {
  real_type random = real_type(generator())/generator.max();
  return xmin + (xmax-xmin)*random;
}

using namespace lapack_wrapper;

#define N_TIMES 1000000

template <int N>
void
testN() {

  cout << "\nSize N = " << N << "\n";

  Malloc<real_type> baseValue("real");
  Malloc<integer>   baseIndex("integer");

  baseValue.allocate(N*N*10);
  baseIndex.allocate(N*10);

  real_type * M1 = baseValue(N*N);
  real_type * M2 = baseValue(N*N);
  real_type * M3 = baseValue(N*N);

  for ( int i = 0; i < N; ++i ) {
    for ( int j = 0; j < N; ++j ) {
      M1[i+j*N] = rand(-1,1);
      M2[i+j*N] = rand(-1,1);
      M3[i+j*N] = rand(-1,1);
    }
  }

  TicToc tm;

  // ===========================================================================

  tm.tic();
  for ( int i = 0; i < N_TIMES; ++i ) {
    gemm(
      NO_TRANSPOSE, NO_TRANSPOSE,
      N, N, N,
      -1.0, M1, N,
      M2, N,
      1.0, M3, N
    );
    copy( N*N, M3, 1, M2, 1);
  }
  tm.toc();
  cout << "MULT = " << tm.elapsed_ms() << " [ms] (lapack)\n";

  cout << "All done!\n";
}



int
main() {

  testN<2>();
  testN<3>();
  testN<4>();
  testN<5>();
  testN<6>();
  testN<7>();
  testN<8>();
  testN<16>();
  testN<32>();
  testN<64>();

  cout << "\n\nAll done!\n";

  return 0;
}
