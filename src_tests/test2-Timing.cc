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
#include <lapack_wrapper/lapack_wrapper_tmpl.hh>
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
typedef double valueType;

static unsigned seed1 = 2;
static std::mt19937 generator(seed1);

static
valueType
rand( valueType xmin, valueType xmax ) {
  valueType random = valueType(generator())/generator.max();
  return xmin + (xmax-xmin)*random;
}

using namespace lapack_wrapper;

template <int N>
void
testMM() {

  int     N_TIMES = (1000000/N);
  double  to_ps   = 1000000.0/N_TIMES;

  fmt::print("\nSize N = {}\n",N);

  Malloc<valueType> baseValue("real");
  Malloc<integer>   baseIndex("integer");

  baseValue.allocate(N*N*10);
  baseIndex.allocate(N*10);

  valueType * M1 = baseValue(N*N);
  valueType * M2 = baseValue(N*N);
  valueType * M3 = baseValue(N*N);

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
  fmt::print("(MM) MULT = {:8.4} [ps] (lapack)\n", to_ps*tm.elapsed_ms() );

  // ===========================================================================

  tm.tic();
  for ( int i = 0; i < N_TIMES; ++i ) {
    MM<valueType,N,N,N,N,N,N>::subTo(M1,M2,M3);
    memcpy( M2, M3, N*N*sizeof(valueType) );
    //Vec2<valueType,N*N,1,1>::copy(M3,M2);
  }
  tm.toc();
  fmt::print("(MM) MULT = {:8.4} [ps] (hand unrolled)\n", to_ps*tm.elapsed_ms() );

  // ===========================================================================

  fmt::print("All done!\n");
}


template <int N>
void
testMv() {

  int     N_TIMES = (1000000/N);
  double  to_ps   = 1000000.0/N_TIMES;

  fmt::print("\nSize N = {}\n",N);

  Malloc<valueType> baseValue("real");
  Malloc<integer>   baseIndex("integer");

  baseValue.allocate(N*N*10);
  baseIndex.allocate(N*10);

  valueType * M = baseValue(N*N);
  valueType * V = baseValue(N);
  valueType * R = baseValue(N);

  for ( int i = 0; i < N; ++i ) {
    V[i] = rand(-1,1);
    R[i] = rand(-1,1);
    for ( int j = 0; j < N; ++j ) {
      M[i+j*N] = rand(-1,1);
    }
  }

  TicToc tm;

  // ===========================================================================

  tm.tic();
  for ( int i = 0; i < N_TIMES; ++i ) {
    gemv(
      NO_TRANSPOSE,
      N, N,
      -1.0, M, N,
      V, 1,
      1.0, R, 1
    );
    copy( N, R, 1, V, 1);
  }
  tm.toc();
  fmt::print("(MV) MULT = {:8.4} [ps] (lapack)\n", to_ps*tm.elapsed_ms() );

  // ===========================================================================

  tm.tic();
  for ( int i = 0; i < N_TIMES; ++i ) {
    Mv<valueType,N,N,N,N,N>::subTo(M,V,R);
    memcpy( V, R, N*sizeof(valueType) );
    //Vec2<valueType,N*N,1,1>::copy(M3,M2);
  }
  tm.toc();
  fmt::print("(MV) MULT = {:8.4} [ps] (hand unrolled)\n", to_ps*tm.elapsed_ms() );

  // ===========================================================================

  fmt::print("All done!\n");
}


template <int N>
void
testCopy() {

  int     N_TIMES = (1000000/N);
  double  to_ps   = 1000000.0/N_TIMES;

  fmt::print("\nSize N = {}\n",N);

  Malloc<valueType> baseValue("real");
  Malloc<integer>   baseIndex("integer");

  baseValue.allocate(N*N*10);
  baseIndex.allocate(N*10);

  valueType * M1 = baseValue(N*N);
  valueType * M2 = baseValue(N*N);
  valueType * M3 = baseValue(N*N);

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
    gecopy( N, N, M1, N, M2, N );
    M1[0] += M2[0];
    gecopy( N, N, M1, N, M2, N );
    M1[0] += M2[0];
    gecopy( N, N, M1, N, M2, N );
    M1[0] += M2[0];
    gecopy( N, N, M1, N, M2, N );
    M1[0] += M2[0];
    gecopy( N, N, M1, N, M2, N );
    M1[0] += M2[0];
  }
  tm.toc();
  fmt::print("(MM) COPY = {:8.4} [ps] (lapack)\n", to_ps*tm.elapsed_ms() );

  fmt::print("All done!\n");
}


static
void
testMvAll() {
  testMv<2>();
  //testMv<3>();
  testMv<4>();
  //testMv<5>();
  testMv<6>();
  //testMv<7>();
  testMv<8>();
  //testMv<9>();
  //testMv<10>();
  //testMv<11>();
  testMv<12>();
  //testMv<13>();
  //testMv<14>();
  //testMv<15>();
  testMv<16>();
  //testMv<17>();
  //testMv<18>();
  //testMv<19>();
  //testMv<20>();
  testMv<100>();
}

static
void
testMMall() {
  testMM<2>();
  //testMM<3>();
  testMM<4>();
  //testMM<5>();
  testMM<6>();
  //testMM<7>();
  testMM<8>();
  //testMM<9>();
  //testMM<10>();
  //testMM<11>();
  testMM<12>();
  //testMM<13>();
  //testMM<14>();
  //testMM<15>();
  testMM<16>();
  //testMM<17>();
  //testMM<18>();
  //testMM<19>();
  //testMM<20>();
  testMM<100>();
}

static
void
testCopyAll() {
  testCopy<2>();
  //testCopy<3>();
  testCopy<4>();
  //testCopy<5>();
  testCopy<6>();
  //testCopy<7>();
  testCopy<8>();
  //testCopy<9>();
  //testCopy<10>();
  //testCopy<11>();
  testCopy<12>();
  //testCopy<13>();
  //testCopy<14>();
  //testCopy<15>();
  testCopy<16>();
  //testCopy<17>();
  //testCopy<18>();
  //testCopy<19>();
  //testCopy<20>();
  testCopy<100>();
}

int
main() {
  testMvAll();
  testMMall();
  testCopyAll();
  fmt::print("\n\nAll done!\n");
  return 0;
}
