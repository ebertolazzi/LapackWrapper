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

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "lapack_wrapper.hh"
#include "lapack_wrapper++.hh"

#include <vector>
#include <limits>
#include <string>

#ifdef LAPACK_WRAPPER_USE_OPENBLAS
  #include <omp.h>
#endif

// workaround for windows
#ifdef _MSC_VER
  #ifdef max
    #undef max
  #endif
  #ifdef min
    #undef min
  #endif
#endif

#include <algorithm>
#include <string>
#include <sstream>

namespace lapack_wrapper {

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  /*
  //              _                 __
  //    __ _  ___| |_ _ ____  __   / /   _
  //   / _` |/ _ \ __| '__\ \/ /  / / | | |
  //  | (_| |  __/ |_| |   >  <  / /| |_| |
  //   \__, |\___|\__|_|  /_/\_\/_/  \__, |
  //   |___/                         |___/
  */
  template <typename REAL>
  integer
  gtx(
    integer M,
    integer N,
    REAL    A[],
    integer LDA,
    integer IPIV[]
  ) {
    // LU DECOMPOSITION, COLUMN INTERCHANGES
    REAL * Ajj = A;
    for ( int j = 0; j < M; Ajj += LDA+1 ) {
      integer MX = iamax( N-j, Ajj, LDA );
      IPIV[j] = MX + j; // C-based
      if ( j < IPIV[j] ) swap( M, A + j*LDA, 1, A + IPIV[j]*LDA, 1 );
      if ( Utils::is_zero(Ajj[0]) ) return j;
      REAL COLM = 1/Ajj[0];
      ++j;
      scal(M-j, COLM, Ajj+1, 1);
      ger(M-j, N-j, -1.0, Ajj+1, 1, Ajj+LDA, LDA, Ajj+LDA+1, LDA);
    }
    return 0;
  }

  template <typename REAL>
  integer
  getrx(
    integer M,      // NUMBER OF ROWS IN A, UNCHANGED
    integer N,      // NUMBER OF COLUMNS IN A, UNCHANGED
    REAL    A[],    // A - ARRAY OF DIMENSION (LDA, N), OVERWRITTEN BY FACTORIZATION
                    // L - UNIT LOVER TRIANGULAR
                    // U - UPPER TRIANGULAR
    integer LDA,    // FIRST DIMENSION OF A AS DECLARED IN THE CALLING (SUB)PROGRAM, UNCHANGED
    integer IPIV[], // ARRAY OF DIMENSION (M) ON EXIT CONTAINS PIVOT INDICES
    integer MB
  ) {
    // COMPUTES LU FACTORIZATION OF M-BY-N MATRIX;
    // PARTIAL PIVOTING COLUMN INTERCHANGES
    // INFO - ERROR INDICATOR, 0 - NORMAL RETURN, POSITIVE VALUE (K) INDICATE
    // THAT U(K, K) = 0.E0 EXACTLY, ALWAYS CHECK AFTER CALL
    if ( M == 0 || N == 0 ) return 0;
    REAL * Ajj = A;
    for ( integer j = 0; j < M; j += MB, Ajj += MB*(LDA+1) ) {
      integer JB = min_index(M-j, MB);
      // FACTORIZE DIAGONAL AND SUBDIAGONAL BLOCKS AND TEST FOR SINGULARITY
      integer INFO = gtx( JB, N-j, Ajj, LDA, IPIV+j );
      if ( INFO != 0 ) return INFO + j;
      // APPLY INTERCHANGES TO PREVIOUS BLOCKS
      integer jjB = j + JB;
      REAL * Aj = A+jjB;
      for ( integer i = j; i < jjB; ++i ) {
        integer IP = (IPIV[i] += j);
        if ( i < IP ) {
          swap( j,     A  + i*LDA, 1, A  + IP*LDA, 1 );
          swap( M-jjB, Aj + i*LDA, 1, Aj + IP*LDA, 1 );
        }
      }
      // COMPUTE SUPERDIAGONAL BLOCK OF U
      trsm(
        SideMultiply::RIGHT,
        ULselect::UPPER,
        Transposition::NO,
        DiagonalType::NON_UNIT,
        M-jjB, JB, 1, Ajj, LDA, Ajj+JB, LDA
      );
      // UPDATE DIAGONAL AND SUBDIAGONAL BLOCKS
      gemm(
        Transposition::NO,
        Transposition::NO,
        M-jjB, N-jjB, JB,
        -1.0, Ajj+JB,         LDA,
              Ajj+JB*LDA,     LDA,
        1.0,  Ajj+JB*(LDA+1), LDA
      );
    }
    return 0;
  }

  template <typename REAL>
  integer
  gty(
    integer M,
    integer N,
    REAL    A[],
    integer LDA,
    integer IPIV[]
  ) {
    // LU DECOMPOSITION, ROW INTERCHANGES
    REAL * Ajj = A;
    for ( int j = 0; j < N; Ajj += LDA+1 ) {
      integer MX = iamax( M-j, Ajj, 1 );
      IPIV[j] = MX + j; // C-based
      if ( j < IPIV[j] ) swap( N, A + j, LDA, A + IPIV[j], LDA );
      if ( Utils::is_zero(Ajj[0]) ) return j;
      REAL ROWM = 1/Ajj[0];
      ++j;
      scal(M-j, ROWM, Ajj+1, 1);
      ger(M-j, N-j, -1.0, Ajj+1, 1, Ajj+LDA, LDA, Ajj+LDA+1, LDA );
    }
    return 0;
  }

  template <typename REAL>
  integer
  getry(
    integer   M,      // NUMBER OF ROWS IN A, UNCHANGED
    integer   N,      // NUMBER OF COLUMNS IN A, UNCHANGED
    REAL    * A,    // A - ARRAY OF DIMENSION (LDA, N), OVERWRITTEN BY FACTORIZATION
                    // L - UNIT LOVER TRIANGULAR
                    // U - UPPER TRIANGULAR
    integer   LDA,    // FIRST DIMENSION OF A AS DECLARED IN THE CALLING (SUB)PROGRAM, UNCHANGED
    integer * IPIV, // ARRAY OF DIMENSION (M) ON EXIT CONTAINS PIVOT INDICES
    integer   NB
  ) {
    // COMPUTES LU FACTORIZATION OF M-BY-N MATRIX;
    // PARTIAL PIVOTING ROW INTERCHANGES
    // INFO - ERROR INDICATOR, 0 - NORMAL RETURN, POSITIVE VALUE (K) INDICATE
    // THAT U(K, K) = 0.E0 EXACTLY, ALWAYS CHECK AFTER CALL
    if ( M == 0 || N == 0 ) return 0;
    REAL * Ajj = A;
    for ( integer j = 0; j < N; j += NB, Ajj += NB*(LDA+1) ) {
      integer JB = min_index(N-j, NB);
      // FACTORIZE DIAGONAL AND SUBDIAGONAL BLOCKS AND TEST FOR SINGULARITY
      integer INFO = gty( M-j, JB, Ajj, LDA, IPIV+j );
      // APPLY INTERCHANGES TO PREVIOUS BLOCKS
      integer jjB = j+JB;
      REAL * Aj = A+jjB*LDA;
      for ( integer i = j; i < jjB; ++i ) {
        integer IP = (IPIV[i] += j);
        if ( i < IP ) {
          swap( j,     A + i,  LDA, A + IP,  LDA );
          swap( N-jjB, Aj + i, LDA, Aj + IP, LDA );
        }
      }
      if ( INFO != 0 ) return INFO + j;
      // COMPUTE SUPERDIAGONAL BLOCK OF U
      trsm(
        SideMultiply::LEFT,
        ULselect::LOWER,
        Transposition::NO,
        DiagonalType::UNIT,
        JB, N-jjB, 1.0, Ajj, LDA, Ajj+JB*LDA, LDA
      );
      // UPDATE DIAGONAL AND SUBDIAGONAL BLOCKS
      gemm(
        Transposition::NO,
        Transposition::NO,
        M-jjB, N-jjB, JB,
        -1.0, Ajj+JB,         LDA,
              Ajj+JB*LDA,     LDA,
         1.0, Ajj+JB*(LDA+1), LDA
      );
    }
    return 0;
  }

  #ifndef LAPACK_WRAPPER_NO_DEBUG
  template integer gtx( integer M, integer N, float A[],  integer LDA, integer IPIV[] );
  template integer gtx( integer M, integer N, double A[], integer LDA, integer IPIV[] );
  template integer gty( integer M, integer N, float A[],  integer LDA, integer IPIV[] );
  template integer gty( integer M, integer N, double A[], integer LDA, integer IPIV[] );

  template integer getrx( integer M, integer N, float A[],  integer LDA, integer IPIV[], integer MB );
  template integer getrx( integer M, integer N, double A[], integer LDA, integer IPIV[], integer MB  );
  template integer getry( integer M, integer N, float A[],  integer LDA, integer IPIV[], integer MB  );
  template integer getry( integer M, integer N, double A[], integer LDA, integer IPIV[], integer MB  );
  #endif

  /*\
   *
   *  Solve
   *
   *  min || T x - b ||^2 + lambda ||x||^2
   *
   *
   *  / * * * * * * \
   *  |   * * * * * |
   *  |     * * * * |
   *  |       * * * |
   *  |         * * |
   *  |           * |
   *  | x - - - - - |
   *  |   x         |
   *  |     x       |
   *  |       x     |
   *  |         x   |
   *  \           x /
   *
  \*/

  template <typename T>
  void
  triTikhonov(
    integer N,
    T const Amat[],
    integer ldA,
    integer nrhs,
    T       RHS[],
    integer ldRHS,
    T       lambda
  ) {

    std::vector<T> Tmat(N*N), line(N), r(nrhs);
    T C, S;

    gecopy( N, N, Amat, ldA, &Tmat.front(), N );

    for ( integer i = 0; i < N; ++i ) {
      std::fill( line.begin()+i, line.end(), T(0) );
      std::fill( r.begin(), r.end(), T(0) );
      line[i] = lambda;
      for ( integer j = i; j < N; ++j ) {
        T * pTjj = &Tmat[j*(N+1)];
        rotg( *pTjj, line[j], C, S );
        if ( N-j-1 > 0 ) rot( N-j-1, pTjj+N, N, &line[j+1], 1, C, S );
        rot( nrhs, RHS+j, ldRHS, &r.front(), 1, C, S );
      }
    }
    // risolvo R
    trsm(
      SideMultiply::LEFT,
      ULselect::UPPER,
      Transposition::NO,
      DiagonalType::NON_UNIT,
      N, nrhs, 1.0, &Tmat.front(), N, RHS, ldRHS
    );
  }

  #ifndef LAPACK_WRAPPER_NO_DEBUG

  #if defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)    || \
      defined(LAPACK_WRAPPER_USE_BLASFEO)
  template <typename T>
  integer
  getc2_tmpl(
    integer N,
    T       A[],
    integer LDA,
    integer IPIV[],
    integer JPIV[]
  ) {
    // Set constants to control overflow
    integer INFO   = 0;
    T       EPS    = lamch<T>("P");
    T       SMLNUM = lamch<T>("S") / EPS;
    T       SMIN   = 0;
    // Factorize A using complete pivoting.
    // Set pivots less than SMIN to SMIN.
    T * Aii = A;
    for ( int II = 0; II < N-1; ++II, Aii += LDA+1 ) {
      // Find max element in matrix A
      T XMAX = 0;
      integer IPV=II, JPV=II;
      for ( int IP = II; IP < N; ++IP ) {
        for ( int JP = II; JP < N; ++JP ) {
          T absA = std::abs( A[IP+JP*LDA] );
          if ( absA > XMAX ) { XMAX = absA; IPV = IP; JPV = JP; }
        }
      }
      if ( II == 0 ) SMIN = std::max( EPS*XMAX, SMLNUM );
      // Swap rows
      IPIV[II] = IPV+1; if ( IPV != II ) swap( N, A+IPV, LDA, A+II, LDA );
      // Swap columns
      JPIV[II] = JPV+1; if ( JPV != II ) swap( N, A+JPV*LDA, 1, A+II*LDA, 1 );
      // Check for singularity
      if ( std::abs(*Aii) < SMIN ) { INFO = II+1; *Aii = SMIN; }
      for ( integer JJ = II+1; JJ < N; ++JJ ) A[JJ+II*LDA] /= *Aii;
      ger( N-II-1, N-II-1, -1, Aii+1, 1, Aii+LDA, LDA, Aii+LDA+1, LDA );
    }
    if ( std::abs(*Aii) < SMIN ) { INFO = N; *Aii = SMIN; }
    return INFO;
  }

  template integer getc2_tmpl(
    integer N,
    real    A[],
    integer LDA,
    integer IPIV[],
    integer JPIV[]
  );

  template integer getc2_tmpl(
    integer    N,
    doublereal A[],
    integer    LDA,
    integer    IPIV[],
    integer    JPIV[]
  );

  template <typename T>
  T
  gesc2_tmpl(
    integer       N,
    T       const A[],
    integer       LDA,
    T             RHS[],
    integer const IPIV[],
    integer const JPIV[]
  ) {
    // Set constants to control overflow
    T EPS    = lamch<T>("P");
    T SMLNUM = lamch<T>("S") / EPS;
    // Apply permutations IPIV to RHS
    for ( integer i = 0; i < N-1; ++i )
      if ( IPIV[i] > i+1 )
        std::swap( RHS[i], RHS[IPIV[i]-1] );
    // Solve for L part
    for ( integer i=0; i < N-1; ++i )
      for ( integer j=i+1; j < N; ++j )
        RHS[j] -= A[j+i*LDA]*RHS[i];
    // Solve for U part
    T SCALE = 1;
    // Check for scaling
    T Rmax = absmax( N, RHS, 1 );
    if ( 2*SMLNUM*Rmax > std::abs(A[(N-1)*(LDA+1)]) ) {
      T TEMP = T(0.5)/Rmax;
      scal( N, TEMP, RHS, 1 );
      SCALE *= TEMP;
    }
    integer i = N;
    while ( i-- > 0 ) {
      T TEMP = 1/A[i*(LDA+1)];
      RHS[i] *= TEMP;
      for ( integer j=i+1; j < N; ++j )
        RHS[i] -= RHS[j]*(A[i+j*LDA]*TEMP);
    }
    // Apply permutations JPIV to the solution (RHS)
    i = N-1;
    while ( i-- > 0 )
    //for ( integer i = 0; i < N-1; ++i )
      if ( JPIV[i] > i+1 )
        std::swap( RHS[i], RHS[JPIV[i]-1] );
    return SCALE;
  }

  template real gesc2_tmpl(
    integer       N,
    real    const A[],
    integer       LDA,
    real          RHS[],
    integer const IPIV[],
    integer const JPIV[]
  );

  template doublereal gesc2_tmpl(
    integer          N,
    doublereal const A[],
    integer          LDA,
    doublereal       RHS[],
    integer    const IPIV[],
    integer    const JPIV[]
  );

  #endif

  template void triTikhonov(
    integer    N,
    real const Tmat[],
    integer    LDT,
    integer    nrhs,
    real       RHS[],
    integer    ldRHS,
    real       lambda
  );

  template void triTikhonov(
    integer          N,
    doublereal const Tmat[],
    integer          LDT,
    integer          nrhs,
    doublereal       RHS[],
    integer          ldRHS,
    doublereal       lambda
  );

  #endif

} // end namespace lapack_wrapper

///
/// eof: lapack_wrapper.cc
///

