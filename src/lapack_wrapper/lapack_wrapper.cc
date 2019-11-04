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
#include "TicToc.hh"

#include <vector>
#include <limits>
#include <string>

// workaround for windows
#ifdef _MSC_VER
  #ifdef max
    #undef max
  #endif
  #ifdef min
    #undef min
  #endif
#endif

#ifdef LAPACK_WRAPPER_OS_WINDOWS
  #include <cstdlib>
  static
  char *
  basename(char const path[] ) {
    static char drive[100];
    static char dir[1024];
    static char fname[256];
    static char ext[128];
    errno_t e = _splitpath_s(
      path,
      drive, 100,
      dir,   1024,
      fname, 256,
      ext,   128
    );
    LAPACK_WRAPPER_ASSERT( e == 0, "lapack_wrapper, basename failed!" );
    return fname;
  }
/*
  #include "StackWalker.h"
  static
  inline
  void
  printStackTrace( FILE *out = stderr ) {
    fprintf( out, "stack trace:\n" );
    StackWalker sw;
    sw.ShowCallstack();
  }
*/
#else
  #include <libgen.h>
#endif

#include <algorithm>
#include <string>
#include <sstream>

#ifndef LAPACK_WRAPPER_OS_WINDOWS
#include <execinfo.h> // for backtrace
#include <dlfcn.h>    // for dladdr
#include <cxxabi.h>   // for __cxa_demangle
#endif

namespace lapack_wrapper {

  #ifndef LAPACK_WRAPPER_OS_WINDOWS
  static
  inline
  std::string
  demang( char const mangled_name[] ) {
    int status = 0 ;
    std::string retval = mangled_name;
    char * name = abi::__cxa_demangle( mangled_name, nullptr, nullptr, &status );
    if ( status == 0 ) retval = name;
    if ( name != nullptr ) std::free(name) ;
    return retval;
  }
  #endif // __linux__

  std::string
  Runtime_Error::grab_backtrace(
    std::string const & reason,
    char const          file[],
    int                 line
  ) const {

    char filename[MAXPATHLEN];
    basename_r( file, filename );
    std::ostringstream ost;
    fmt::print(
      ost, "\n{}\nOn File:{}:{}\n",
      reason, filename, line
    );

    #ifndef LAPACK_WRAPPER_OS_WINDOWS
    fmt::print(ost, "stack trace:\n");

    //  record stack trace upto 128 frames
    void *callstack[128] = {};

    // collect stack frames
    int frames = backtrace( callstack, 128);

    // get the human-readable symbols (mangled)
    char** strs = backtrace_symbols( callstack, frames);

    for ( int i = 1; i < frames; ++i) {
      char functionSymbol[1024] = {};
      char moduleName[1024]     = {};
      int  offset               = 0;
      char addr[48]             = {};

      /*

      Typically this is how the backtrace looks like:

      0   <app/lib-name>     0x0000000100000e98 _Z5tracev + 72
      1   <app/lib-name>     0x00000001000015c1 _ZNK7functorclEv + 17
      2   <app/lib-name>     0x0000000100000f71 _Z3fn0v + 17
      3   <app/lib-name>     0x0000000100000f89 _Z3fn1v + 9
      4   <app/lib-name>     0x0000000100000f99 _Z3fn2v + 9
      5   <app/lib-name>     0x0000000100000fa9 _Z3fn3v + 9
      6   <app/lib-name>     0x0000000100000fb9 _Z3fn4v + 9
      7   <app/lib-name>     0x0000000100000fc9 _Z3fn5v + 9
      8   <app/lib-name>     0x0000000100000fd9 _Z3fn6v + 9
      9   <app/lib-name>     0x0000000100001018 main + 56
      10  libdyld.dylib      0x00007fff91b647e1 start + 0

      */

      // split the string, take out chunks out of stack trace
      // we are primarily interested in module, function and address
      sscanf(
        strs[i], "%*s %s %s %s %*s %d",
        moduleName, addr, functionSymbol, &offset
      );

      //  if this is a C++ library, symbol will be demangled
      //  on success function returns 0
      //
      fmt::print( ost, "{:2} ({:30})  {} — {} + {}\n",
        i, moduleName, addr, demang( functionSymbol ), offset
      );
    }
    free(strs);
    #endif
    return ost.str();
  }

  #if 0
  std::string
  Runtime_Error::grab_backtrace(
    std::string const & reason,
    char const          file[],
    int                 line
  ) const {
    char filename[MAXPATHLEN];
    basename_r( file, filename );
    std::ostringstream ost;
    ost << reason << "\nOn File: " << filename << ":" << line << "\nStack Trace:\n";
    void *callstack[32];
    int const nMaxFrames = 32;
    int       nFrames    = ::backtrace( callstack, nMaxFrames );
    char **   symbols    = backtrace_symbols(callstack, nFrames);
    for ( int i = 0; i < nFrames; ++i ) {
      Dl_info info;
      if ( dladdr(callstack[i], &info) && info.dli_sname ) {
        unsigned offs = unsigned((char *)callstack[i] - (char *)info.dli_saddr);
        ost << symbols[i] << "\n";
        fmt::print(
          ost, "{:2} +{:04X} {}\n",
          i, offs, demang(info.dli_sname)
        );
      } else {
        ost << symbols[i] << "\n";
      }
    }
    free(symbols);
    if ( nFrames == nMaxFrames ) ost << "[truncated]\n";
    return ost.str();
  }
  #endif

  #if defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)      || \
      defined(LAPACK_WRAPPER_USE_OPENBLAS)   || \
      defined(LAPACK_WRAPPER_USE_BLASFEO)
    CBLAS_TRANSPOSE trans_cblas[3] = { CblasNoTrans, CblasTrans, CblasConjTrans };
    CBLAS_UPLO      uplo_cblas[2]  = { CblasUpper, CblasLower };
    CBLAS_DIAG      diag_cblas[2]  = { CblasUnit, CblasNonUnit };
    CBLAS_SIDE      side_cblas[2]  = { CblasLeft, CblasRight };
  #endif

  character const *trans_blas[3]   = { "N", "T", "C" };
  character const *uplo_blas[2]    = { "Upper", "Lower" };
  character const *diag_blas[2]    = { "Unit",  "NonUnit" };
  character const *side_blas[2]    = { "Left",  "Right" };
  character const *balance_blas[4] = {
    "No Balance",
    "Permute Only",
    "Scale Only",
    "Both permute and scale"
  };

  character const *job_blas[4]   = { "A", "S", "O", "N" };
  character const *sense_blas[4] = { "N", "E", "V", "B" };

  character const *direct_blas[2] = { "Forward", "Backward" };
  character const *store_blas[2]  = { "ColumnWise", "RowWise" };

  character const *mtype_blas[5]  = {
    "G", // A is a full matrix.
    "L", // A is a lower triangular matrix.
    "U", // A is an upper triangular matrix.
    "H", // A is an upper Hessenberg matrix.
    "B"  // A is a symmetric band matrix with lower bandwidth KL
  };

  character const *equilibrate_blas[4]  = { "N", "R", "C", "B" };

  //============================================================================
  /*    __                       _ _   _       _   _ 
  //   / _| ___  _   _ _ __   __| | \ | | __ _| \ | |
  //  | |_ / _ \| | | | '_ \ / _` |  \| |/ _` |  \| |
  //  |  _| (_) | |_| | | | | (_| | |\  | (_| | |\  |
  //  |_|  \___/ \__,_|_| |_|\__,_|_| \_|\__,_|_| \_|
  */
  //! check if the vector `pv` os size `DIM` contains only regular floats
  bool
  foundNaN( doublereal const pv[], integer DIM ) {
    for ( integer i = 0; i < DIM; ++i )
      if ( !isRegular(pv[i]) )
        return true;
    return false;
  }

  bool
  foundNaN( real const pv[], integer DIM ) {
    for ( integer i = 0; i < DIM; ++i )
      if ( !isRegular(pv[i]) )
        return true;
    return false;
  }

  /*       _               _    _   _       _   _ 
  //   ___| |__   ___  ___| | _| \ | | __ _| \ | |
  //  / __| '_ \ / _ \/ __| |/ /  \| |/ _` |  \| |
  // | (__| | | |  __/ (__|   <| |\  | (_| | |\  |
  //  \___|_| |_|\___|\___|_|\_\_| \_|\__,_|_| \_|
  */
  
  #define LINE_LINE_LINE_LINE "--------------------------------------------------------------------------------"

  //! check if the vector `pv` os size `DIM` contains only regular floats. If not an error is issued
  void
  checkNaN(
    doublereal const pv[],
    char       const v_name[],
    integer          DIM,
    integer          line,
    char       const file[]
  ) {
    for ( integer i = 0; i < DIM; ++i ) {
      if ( isInfinite(pv[i]) ) {
        LAPACK_WRAPPER_ERROR(
          LINE_LINE_LINE_LINE <<
          "\n(" << basename(const_cast<char*>(file)) << ':' << line <<
          ") found Infinity at " << v_name << "[" << i << "]\n" <<
          LINE_LINE_LINE_LINE
        )
      } else if ( isNaN(pv[i]) ) {
        LAPACK_WRAPPER_ERROR(
          LINE_LINE_LINE_LINE <<
          "\n(" << basename(const_cast<char*>(file)) << ':' << line <<
          ") found NaN at " << v_name << "[" << i << "]\n" <<
          LINE_LINE_LINE_LINE
        );
      }
    }
  }

  void
  checkNaN(
    real const pv[],
    char const v_name[],
    integer    DIM,
    integer    line,
    char const file[]
  ) {
    for ( integer i = 0; i < DIM; ++i ) {
      if ( isInfinite(pv[i]) ) {
        LAPACK_WRAPPER_ERROR(
          LINE_LINE_LINE_LINE <<
          "\n(" << basename(const_cast<char*>(file)) << ':' << line <<
          ") found Infinity at " << v_name << "[" << i << "]\n" <<
          LINE_LINE_LINE_LINE
        )
      } else if ( isNaN(pv[i]) ) {
        LAPACK_WRAPPER_ERROR(
          LINE_LINE_LINE_LINE <<
          "\n(" << basename(const_cast<char*>(file)) << ':' << line <<
          ") found NaN at " << v_name << "[" << i << "]\n" <<
          LINE_LINE_LINE_LINE
        );
      }
    }
  }

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
      if ( isZero(Ajj[0]) ) return j;
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
        RIGHT,
        UPPER,
        NO_TRANSPOSE,
        NON_UNIT,
        M-jjB, JB, 1, Ajj, LDA, Ajj+JB, LDA
      );
      // UPDATE DIAGONAL AND SUBDIAGONAL BLOCKS
      gemm(
        NO_TRANSPOSE,
        NO_TRANSPOSE,
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
      if ( isZero(Ajj[0]) ) return j;
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
    integer M,      // NUMBER OF ROWS IN A, UNCHANGED
    integer N,      // NUMBER OF COLUMNS IN A, UNCHANGED
    REAL    A[],    // A - ARRAY OF DIMENSION (LDA, N), OVERWRITTEN BY FACTORIZATION
                    // L - UNIT LOVER TRIANGULAR
                    // U - UPPER TRIANGULAR
    integer LDA,    // FIRST DIMENSION OF A AS DECLARED IN THE CALLING (SUB)PROGRAM, UNCHANGED
    integer IPIV[], // ARRAY OF DIMENSION (M) ON EXIT CONTAINS PIVOT INDICES
    integer NB
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
        LEFT,
        LOWER,
        NO_TRANSPOSE,
        UNIT,
        JB, N-jjB, 1.0, Ajj, LDA, Ajj+JB*LDA, LDA
      );
      // UPDATE DIAGONAL AND SUBDIAGONAL BLOCKS
      gemm(
        NO_TRANSPOSE,
        NO_TRANSPOSE,
        M-jjB, N-jjB, JB,
        -1.0, Ajj+JB,         LDA,
              Ajj+JB*LDA,     LDA,
         1.0, Ajj+JB*(LDA+1), LDA
      );
    }
    return 0;
  }
  
  template integer gtx( integer M, integer N, float A[],  integer LDA, integer IPIV[] );
  template integer gtx( integer M, integer N, double A[], integer LDA, integer IPIV[] );
  template integer gty( integer M, integer N, float A[],  integer LDA, integer IPIV[] );
  template integer gty( integer M, integer N, double A[], integer LDA, integer IPIV[] );

  template integer getrx( integer M, integer N, float A[],  integer LDA, integer IPIV[], integer MB );
  template integer getrx( integer M, integer N, double A[], integer LDA, integer IPIV[], integer MB  );
  template integer getry( integer M, integer N, float A[],  integer LDA, integer IPIV[], integer MB  );
  template integer getry( integer M, integer N, double A[], integer LDA, integer IPIV[], integer MB  );

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
      LEFT,
      UPPER,
      NO_TRANSPOSE,
      NON_UNIT,
      N, nrhs, 1.0, &Tmat.front(), N, RHS, ldRHS
    );
  }

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
    integer INFO = 0;
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

} // end namespace lapack_wrapper

///
/// eof: lapack_wrapper.cc
///

