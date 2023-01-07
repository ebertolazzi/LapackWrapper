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

///
/// file: lapack_wrapper.hh
///

#include "../lapack_wrapper_config.hh"

#ifndef LAPACK_WRAPPER_dot_HH
#define LAPACK_WRAPPER_dot_HH

#define LAPACK_WRAPPER_MAJOR_VERSION 0
#define LAPACK_WRAPPER_MINOR_VERSION 1

// possible user modification
// LAPACK_WRAPPER_USE_SYSTEM_EIGEN
// LAPACK_WRAPPER_USE_SYSTEM_SUPERLU
// LAPACK_WRAPPER_USE_SYSTEM_OPENBLAS

/*
// [sd]lamch
// [sd]copy
// [sd]swap [sd]swaps
// [sd]scal
// [sd]axpy
// [sd]zero
// [sd]rot, [sd]rotg, [sd]lartg
// [sd]lapy
// [sd]ger
// [sd]gemv
// [sd]gemm
// [sd]nrm2
// [sd]asum
// [sd]amax
// [sd]dot
// [sd]trmv, [sd]trsv, [sd]trmm, [sd]trsm
// [sd]tbmv, [sd]tbsv
// [sd]gttrf, [sd]gttrs, [sd]pttrf, [sd]pttrs
// [sd]geadd
// [sd]getrf, [sd]getrs, [sd]gesv, [sd]getc2, [sd]gesc2
// [sd]gecon, [sd]geequ, [sd]laqge
// [sd]gecopy, [sd]gezero, [sd]geid
// [sd]normInf [sd]norm1, [sd]normF, [sd]maxabs
// [sd]lascl
// [sd]gev, [sd]ggev, [sd]ggevx (autovalori)
// [sd]gesvd, [sd]gesdd
// [sd]gelsd, [sd]gelss, [sd]gelsy
// [sd]gbtrf, [sd]gbtrs
// [sd]ormqr
// [sd]larft, [sd]larfg, [sd]larfb
// [sd]geqrf, [sd]geqr2, [sd]geqp3
// [sd]tzrzf, [sd]ormrz
// [sd]laic1
// rankEstimate
*/

#include <cmath>
#include <cstring>
#include <limits>
#include <algorithm>

#ifdef LAPACK_WRAPPER_USE_LAPACK2
  #define LAPACK_WRAPPER_USE_LAPACK 1
#endif

// if not choosed the linear algebra package give an error
#if !defined(LAPACK_WRAPPER_USE_ACCELERATE) && \
    !defined(LAPACK_WRAPPER_USE_ATLAS)      && \
    !defined(LAPACK_WRAPPER_USE_OPENBLAS)   && \
    !defined(LAPACK_WRAPPER_USE_BLASFEO)    && \
    !defined(LAPACK_WRAPPER_USE_LAPACK)     && \
    !defined(LAPACK_WRAPPER_USE_MKL)
  #error "linear algebra package NOT selected"
#endif

#ifdef UTILS_OS_WINDOWS
  // se in windows includo PRIMA windows.h per evitare conflitti!
  #ifdef _MSC_VER
    #include <stdint.h>
  #else
    namespace lapack_wrapper {
      typedef          __int8  int8_t;
      typedef          __int16 int16_t;
      typedef          __int32 int32_t;
      typedef          __int64 int64_t;
      typedef unsigned __int8  uint8_t;
      typedef unsigned __int16 uint16_t;
      typedef unsigned __int32 uint32_t;
      typedef unsigned __int64 uint64_t;
    }
  #endif
#else
  #include <cstdint>
#endif

// set macro for used LAPACK name

#if defined(LAPACK_WRAPPER_USE_ACCELERATE)
  #ifdef UTILS_ARCH64
    #define LAPACK_WRAPPER_LAPACK_NAME "Accelerate OSX (x64)"
  #else
    #define LAPACK_WRAPPER_LAPACK_NAME "Accelerate OSX (x86)"
  #endif
#elif defined(LAPACK_WRAPPER_USE_LAPACK)
  #if defined(_DEBUG) || defined(DEBUG)
    #ifdef UTILS_ARCH64
      #define LAPACK_WRAPPER_LAPACK_NAME "Standard BLAS/LAPACK (x64,DEBUG)"
    #else
      #define LAPACK_WRAPPER_LAPACK_NAME "Standard BLAS/LAPACK (x86,DEBUG)"
    #endif
  #else
    #ifdef UTILS_ARCH64
      #define LAPACK_WRAPPER_LAPACK_NAME "Standard BLAS/LAPACK (x64)"
    #else
      #define LAPACK_WRAPPER_LAPACK_NAME "Standard BLAS/LAPACK (x86)"
    #endif
  #endif
#elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
  #ifdef UTILS_ARCH64
    #define LAPACK_WRAPPER_LAPACK_NAME "OPENBLAS (x64,RELEASE)"
  #else
    #define LAPACK_WRAPPER_LAPACK_NAME "OPENBLAS (x86,RELEASE)"
  #endif
#elif defined(LAPACK_WRAPPER_USE_BLASFEO)
  #ifdef UTILS_ARCH64
    #define LAPACK_WRAPPER_LAPACK_NAME "BLASFEO (x64,RELEASE)"
  #else
    #define LAPACK_WRAPPER_LAPACK_NAME "BLASFEO (x86,RELEASE)"
  #endif
#elif defined(LAPACK_WRAPPER_USE_MKL)
  #ifdef UTILS_ARCH64
    #define LAPACK_WRAPPER_LAPACK_NAME "MKL (x64)"
  #else
    #define LAPACK_WRAPPER_LAPACK_NAME "MKL (x86)"
  #endif
#else
  #define LAPACK_WRAPPER_LAPACK_NAME "UNKNOWN LAPACK LIB"
#endif

// find Headers for Lapack/Blas

#if defined(LAPACK_WRAPPER_USE_ACCELERATE)

  #include <Accelerate/Accelerate.h>
  #define CBLASNAME(A)   cblas_##A
  #define CLAPACKNAME(A) A##_

#elif defined(LAPACK_WRAPPER_USE_ATLAS)

  // atlas 3.6.0
  extern "C" {
    #include <atlas/cblas.h>
    #include <atlas/clapack.h>
  }

  #define CBLASNAME(A) cblas_##A
  #ifdef UTILS_OS_OSX
    #define CLAPACKNAME(A) A##_
  #else
    #define CLAPACKNAME(A) clapack_##A
  #endif
  #define LAPACK_F77NAME(A) A##_

  #ifndef BLASFUNC
    #define BLASFUNC(A) LAPACK_F77NAME(A)
  #endif

#elif defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
      defined(LAPACK_WRAPPER_USE_BLASFEO)
  // workaround for OPENBLAS on OSX
  #ifdef UTILS_OS_OSX
    #ifdef __clang__
      #pragma clang diagnostic ignored "-Wreserved-id-macro"
    #endif
    #define __STDC_VERSION__ __STDC__
  #endif

  #include <complex>
  #define FORCE_OPENBLAS_COMPLEX_STRUCT
  #define lapack_complex_float  std::complex<float>
  #define lapack_complex_double std::complex<double>

  #ifdef UTILS_OS_WINDOWS
    // on windows the defaul is too use local openblas
    #ifdef LAPACK_WRAPPER_USE_SYSTEM_OPENBLAS
      #include <openblas/cblas.h>
      #include <openblas/lapack.h>
    #else
      #ifdef UTILS_ARCH64
        #include "../openblas/x64/cblas.h"
        #include "../openblas/x64/lapack.h"
      #else
        #include "../openblas/x86/cblas.h"
        #include "../openblas/x86/lapack.h"
      #endif
    #endif
  #elif defined(UTILS_OS_LINUX)
    // on linux/osx the defaul is too use SYSTEM openblas
    #ifdef LAPACK_WRAPPER_DO_NOT_USE_SYSTEM_OPENBLAS
      #include "../openblas/cblas.h"
      #include "../openblas/lapack.h"
    #else
      // must use subdirectory to avoid conflict with atlas/lapack...
      #include <openblas/cblas.h>
      #include <openblas/lapack.h>
    #endif
  #else
    // on linux/osx the defaul is too use SYSTEM openblas
    #ifdef LAPACK_WRAPPER_DO_NOT_USE_SYSTEM_OPENBLAS
      #include "../openblas/cblas.h"
      #include "../openblas/lapack.h"
    #else
      // OSX   -I/usr/local/opt/openblas/include/
      #include <cblas.h>
      #include <lapacke.h>
    #endif
  #endif

  #define CBLASNAME(A)       cblas_##A
  #define LAPACK_F77NAME(A)  LAPACK_##A

  #ifdef LAPACK_WRAPPER_USE_BLASFEO
    #include <s_blas.h>
    #include <d_blas.h>
  #endif

#elif defined(LAPACK_WRAPPER_USE_LAPACK)

  #ifndef lapack_int
  #define lapack_int int32_t
  #endif

  #define LAPACK_F77NAME(A) A##_
  #define BLASFUNC(A)       A##_

#elif defined(LAPACK_WRAPPER_USE_MKL)

  #include <complex>
  #define MKL_Complex8  std::complex<float>
  #define MKL_Complex16 std::complex<double>

  #include <mkl_blas.h>
  #include <mkl_lapack.h>

#else
  #error "You must select the linear algebra packages used!"
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wpartial-availability"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#endif

namespace lapack_wrapper {

  #ifndef LAPACK_WRAPPER_NO_DEBUG

  #if defined(LAPACK_WRAPPER_USE_ACCELERATE)
    using integer          = __CLPK_integer;
    using real             = __CLPK_real;
    using doublereal       = __CLPK_doublereal;
    using character        = char;
    using return_precision = doublereal;
  #elif defined(LAPACK_WRAPPER_USE_ATLAS)
    using integer          = int;
    using real             = float;
    using doublereal       = double;
    using character        = char;
    using return_precision = doublereal;
  #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
    using integer          = lapack_int;
    using real             = float;
    using doublereal       = double;
    using character        = char;
    using return_precision = double;
  #elif defined(LAPACK_WRAPPER_USE_BLASFEO)
    using integer          = blasint;
    using real             = float;
    using doublereal       = double;
    using character        = char;
    using return_precision = double;
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)
    using integer          = lapack_int;
    using real             = float;
    using doublereal       = double;
    using character        = char;
    using return_precision = doublereal;
  #elif defined(LAPACK_WRAPPER_USE_MKL)
    using integer          = MKL_INT;
    using real             = float;
    using doublereal       = double;
    using character        = char;
    using return_precision = double;
  #else
    #error "You must select the linear algebra packages used!"
  #endif

  #endif

  extern character const *balance_blas[4];
  extern character const *equilibrate_blas[4];

  //============================================================================

  using Transposition = enum class Transposition : integer {
    NO_TRANSPOSE        = 0,
    TRANSPOSE           = 1,
    CONJUGATE_TRANSPOSE = 2
  };

  inline
  char *
  to_blas( Transposition const & T ) {
    switch( T ) {
    case Transposition::NO_TRANSPOSE:        return const_cast<char*>("NO_TRANSPOSE");
    case Transposition::TRANSPOSE:           return const_cast<char*>("TRANSPOSE");
    case Transposition::CONJUGATE_TRANSPOSE: return const_cast<char*>("CONJUGATE_TRANSPOSE");
    }
  }

  inline
  CBLAS_TRANSPOSE
  to_cblas( Transposition const & T ) {
    switch( T ) {
    case Transposition::NO_TRANSPOSE:        return CblasNoTrans;
    case Transposition::TRANSPOSE:           return CblasTrans;
    case Transposition::CONJUGATE_TRANSPOSE: return CblasConjTrans;
    }
  }

  //============================================================================

  using ULselect = enum class ULselect : integer {
    UPPER = 0,
    LOWER = 1
  };

  inline
  char *
  to_blas( ULselect const & UL ) {
    switch( UL ) {
    case ULselect::UPPER: return const_cast<char*>("UPPER");
    case ULselect::LOWER: return const_cast<char*>("LOWER");
    }
  }

  inline
  CBLAS_UPLO
  to_cblas( ULselect const & UL ) {
    switch( UL ) {
    case ULselect::UPPER: return CblasUpper;
    case ULselect::LOWER: return CblasLower;
    }
  }

  //============================================================================

  using DiagonalType = enum class DiagonalType : integer {
    UNIT     = 0,
    NON_UNIT = 1
  };

  inline
  char *
  to_blas( DiagonalType const & DT ) {
    switch( DT ) {
    case DiagonalType::UNIT:     return const_cast<char*>("UNIT");
    case DiagonalType::NON_UNIT: return const_cast<char*>("NON_UNIT");
    }
  }

  inline
  CBLAS_DIAG
  to_cblas( DiagonalType const & DT ) {
    switch( DT ) {
    case DiagonalType::UNIT:     return CblasUnit;
    case DiagonalType::NON_UNIT: return CblasNonUnit;
    }
  }

  //============================================================================

  using SideMultiply = enum class SideMultiply : integer {
    LEFT   = 0,
    RIGHT  = 1
  };

  inline
  char *
  to_blas( SideMultiply const & UL ) {
    switch( UL ) {
    case SideMultiply::LEFT:  return const_cast<char*>("LEFT");
    case SideMultiply::RIGHT: return const_cast<char*>("RIGHT");
    }
  }

  inline
  CBLAS_SIDE
  to_cblas( SideMultiply const & UL ) {
    switch( UL ) {
    case SideMultiply::LEFT:  return CblasLeft;
    case SideMultiply::RIGHT: return CblasRight;
    }
  }

  //============================================================================

  using JobType = enum class JobType : integer {
    ALL     = 0, // 'A'
    REDUCED = 1, // 'S'
    INPLACE = 2, // 'O'
    NO_JOB  = 3  // 'N'
  };

  inline
  char *
  to_blas( JobType const & JOB ) {
    switch( JOB ) {
    case JobType::ALL:     return const_cast<char*>("ALL");
    case JobType::REDUCED: return const_cast<char*>("S"); // REDUCED
    case JobType::INPLACE: return const_cast<char*>("O"); // INPLACE
    case JobType::NO_JOB:  return const_cast<char*>("NO_JOB");
    }
  }

  //============================================================================

  using SenseType = enum class SenseType : integer {
    NONE                         = 0, // 'N'
    EIGENVALUES_ONLY             = 1, // 'E'
    EIGENVECTORS_ONLY            = 2, // 'V'
    EIGENVALUES_AND_EIGENVECTORS = 3  // 'B'
  };

  inline
  char *
  to_blas( SenseType const & ST ) {
    switch( ST ) {
    case SenseType::NONE:                         return const_cast<char*>("NONE");
    case SenseType::EIGENVALUES_ONLY:             return const_cast<char*>("EIGENVALUES"); // REDUCED
    case SenseType::EIGENVECTORS_ONLY:            return const_cast<char*>("VECTORS"); // INPLACE
    case SenseType::EIGENVALUES_AND_EIGENVECTORS: return const_cast<char*>("BOTH_EIGENVALUES_AND_EIGENVECTORS");
    }
  }

  //============================================================================

  using DirectionType = enum class DirectionType : integer {
    FORWARD  = 0,
    BACKWARD = 1
  };

  inline
  char *
  to_blas( DirectionType const & DIR ) {
    switch( DIR ) {
    case DirectionType::FORWARD:  return const_cast<char*>("FORWARD");
    case DirectionType::BACKWARD: return const_cast<char*>("BACKWARD"); // REDUCED
    }
  }

  //============================================================================

  using StorageType = enum class StorageType : integer {
    COLUMNWISE = 0,
    ROWWISE    = 1
  };

  inline
  char *
  to_blas( StorageType const & ST ) {
    switch( ST ) {
    case StorageType::COLUMNWISE: return const_cast<char*>("COLUMNWISE");
    case StorageType::ROWWISE:    return const_cast<char*>("ROWWISE"); // REDUCED
    }
  }

  //============================================================================

  using MatrixType = enum class MatrixType : integer {
    FULL             = 0,
    LOWER_TRIANGULAR = 1,
    UPPER_TRIANGULAR = 2,
    HESSENBERG       = 3,
    BANDED           = 4
  };

  inline
  char *
  to_blas( MatrixType const & MT ) {
    switch( MT ) {
    case MatrixType::FULL:             return const_cast<char*>("GENERAL FULL");
    case MatrixType::LOWER_TRIANGULAR: return const_cast<char*>("LOWER_TRIANGULAR");
    case MatrixType::UPPER_TRIANGULAR: return const_cast<char*>("UPPER_TRIANGULAR");
    case MatrixType::HESSENBERG:       return const_cast<char*>("HESSENBERG");
    case MatrixType::BANDED:           return const_cast<char*>("BANDED");
    }
  }

  //============================================================================

  using BalanceType = enum {
    NO_BALANCE        = 0, // 'N'
    PERMUTE_ONLY      = 1, // 'P'
    SCALE_ONLY        = 2, // 'S'
    PERMUTE_AND_SCALE = 3  // 'B'
  };

  extern char const *BalanceType_name[];

  //============================================================================

  using EquilibrationType = enum {
    NO_EQUILIBRATE      = 0,
    EQUILIBRATE_ROWS    = 1,
    EQUILIBRATE_COLUMNS = 2,
    EQUILIBRATE_BOTH    = 3
  };

  extern char const *EquilibrationType_name[];

  //============================================================================

  static
  inline
  integer
  min_index( integer a, integer b )
  { return a < b ? a : b; }

  static
  inline
  integer
  max_index( integer a, integer b )
  { return a > b ? a : b; }

  extern "C" {
    #ifdef LAPACK_WRAPPER_USE_ACCELERATE
    int xerbla_( character const * what, integer * info, int );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    int xerbla_( character const * what, integer * info, int );
    #else
    int LAPACK_F77NAME(xerbla)( character const * what, integer * info, int );
    #endif
  };

  inline
  void
  xerbla( character const * WHAT, integer info ) {
    int len = int(strlen(WHAT));
    #ifdef LAPACK_WRAPPER_USE_ACCELERATE
    xerbla_( WHAT, &info, len );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    xerbla_( WHAT, &info, len );
    #else
    LAPACK_F77NAME(xerbla)( WHAT, &info, len );
    #endif
  }

  /*
  //   _                      _
  //  | | __ _ _ __ ___   ___| |__
  //  | |/ _` | '_ ` _ \ / __| '_ \
  //  | | (_| | | | | | | (__| | | |
  //  |_|\__,_|_| |_| |_|\___|_| |_|
  */

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {
    real
    LAPACK_F77NAME(slamch)( character * what );

    doublereal
    LAPACK_F77NAME(dlamch)( character * what );
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DLAMCH determines double precision machine parameters.
   *
   *  Arguments
   *  =========
   *
   *  CMACH   (input) CHARACTER*1
   *          Specifies the value to be returned by DLAMCH:
   *          = 'E' or 'e',   DLAMCH := eps
   *          = 'S' or 's ,   DLAMCH := sfmin
   *          = 'B' or 'b',   DLAMCH := base
   *          = 'P' or 'p',   DLAMCH := eps*base
   *          = 'N' or 'n',   DLAMCH := t
   *          = 'R' or 'r',   DLAMCH := rnd
   *          = 'M' or 'm',   DLAMCH := emin
   *          = 'U' or 'u',   DLAMCH := rmin
   *          = 'L' or 'l',   DLAMCH := emax
   *          = 'O' or 'o',   DLAMCH := rmax
   *
   *          where
   *
   *          eps   = relative machine precision
   *          sfmin = safe minimum, such that 1/sfmin does not overflow
   *          base  = base of the machine
   *          prec  = eps*base
   *          t     = number of (base) digits in the mantissa
   *          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
   *          emin  = minimum exponent before (gradual) underflow
   *          rmin  = underflow threshold - base**(emin-1)
   *          emax  = largest exponent before overflow
   *          rmax  = overflow threshold  - (base**emax)*(1-eps)
   *
   * =====================================================================
  \*/
  #ifndef LAPACK_WRAPPER_USE_MKL
  // trick for function that differs only in return type
  template <typename T> T lamch( character const WHAT[] );

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  template <>
  inline
  real
  lamch( character const WHAT[] )
  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  { return LAPACK_F77NAME(slamch)( const_cast<character*>(WHAT) ); }
  #elif defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_BLASFEO)
  { return BLASFUNC(slamch)( const_cast<character*>(WHAT) ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
  { return real(CLAPACKNAME(slamch)( const_cast<character*>(WHAT) )); }
  #else
  { return LAPACKNAME(slamch)( const_cast<character*>(WHAT) ); }
  #endif

  template <>
  inline
  doublereal
  lamch( character const WHAT[] )
  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  { return LAPACK_F77NAME(dlamch)( const_cast<character*>(WHAT) ); }
  #elif defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_BLASFEO)
  { return BLASFUNC(dlamch)( const_cast<character*>(WHAT) ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
  { return CLAPACKNAME(dlamch)( const_cast<character*>(WHAT) ); }
  #else
  { return LAPACKNAME(dlamch)( const_cast<character*>(WHAT) ); }
  #endif

  #endif

  #endif

  /*
  //   _ _
  //  (_) | __ _  ___ _ ____   __
  //  | | |/ _` |/ _ \ '_ \ \ / /
  //  | | | (_| |  __/ | | \ V /
  //  |_|_|\__,_|\___|_| |_|\_/
  */

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {
    integer
    LAPACK_F77NAME(ilaenv)(
      integer   const * ISPEC,
      character const   NAME[],
      character const   OPTS[],
      integer         * N1,
      integer         * N2,
      integer         * N3,
      integer         * N4,
      size_t          * len_NAME,
      size_t          * len_OPTS
    );
  }
  #endif

  #if defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
      defined(LAPACK_WRAPPER_USE_BLASFEO)
  // non esiste prototipo in open blas, ma supportata
  extern "C" {
    integer
    BLASFUNC(ilaenv)(
      integer   const * ISPEC,
      character const   NAME[],
      character const   OPTS[],
      integer         * N1,
      integer         * N2,
      integer         * N3,
      integer         * N4,
      size_t          * len_NAME,
      size_t          * len_OPTS
    );
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  ILAENV is called from the LAPACK routines to choose problem-dependent
   *  parameters for the local environment.  See ISPEC for a description of
   *  the parameters.
   *
   *  ILAENV returns an INTEGER
   *  if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
   *  if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
   *
   *  This version provides a set of parameters which should give good,
   *  but not optimal, performance on many of the currently available
   *  computers.  Users are encouraged to modify this subroutine to set
   *  the tuning parameters for their particular machine using the option
   *  and problem size information in the arguments.
   *
   *  This routine will not function correctly if it is converted to all
   *  lower case.  Converting it to all upper case is allowed.
   *
   *  Arguments
   *  =========
   *
   *  ISPEC   (input) INTEGER
   *          Specifies the parameter to be returned as the value of
   *          ILAENV.
   *          = 1: the optimal blocksize; if this value is 1, an unblocked
   *               algorithm will give the best performance.
   *          = 2: the minimum block size for which the block routine
   *               should be used; if the usable block size is less than
   *               this value, an unblocked routine should be used.
   *          = 3: the crossover point (in a block routine, for N less
   *               than this value, an unblocked routine should be used)
   *          = 4: the number of shifts, used in the nonsymmetric
   *               eigenvalue routines (DEPRECATED)
   *          = 5: the minimum column dimension for blocking to be used;
   *               rectangular blocks must have dimension at least k by m,
   *               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
   *          = 6: the crossover point for the SVD (when reducing an m by n
   *               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
   *               this value, a QR factorization is used first to reduce
   *               the matrix to a triangular form.)
   *          = 7: the number of processors
   *          = 8: the crossover point for the multishift QR method
   *               for nonsymmetric eigenvalue problems (DEPRECATED)
   *          = 9: maximum size of the subproblems at the bottom of the
   *               computation tree in the divide-and-conquer algorithm
   *               (used by xGELSD and xGESDD)
   *          =10: ieee NaN arithmetic can be trusted not to trap
   *          =11: infinity arithmetic can be trusted not to trap
   *          12 <= ISPEC <= 16:
   *               xHSEQR or one of its subroutines,
   *               see IPARMQ for detailed explanation
   *
   *  NAME    (input) CHARACTER*(*)
   *          The name of the calling subroutine, in either upper case or
   *          lower case.
   *
   *  OPTS    (input) CHARACTER*(*)
   *          The character options to the subroutine NAME, concatenated
   *          into a single character string.  For example, UPLO = 'U',
   *          TRANS = 'T', and DIAG = 'N' for a triangular routine would
   *          be specified as OPTS = 'UTN'.
   *
   *  N1      (input) INTEGER
   *  N2      (input) INTEGER
   *  N3      (input) INTEGER
   *  N4      (input) INTEGER
   *          Problem dimensions for the subroutine NAME; these may not all
   *          be required.
   *
   *  Further Details
   *  ===============
   *
   *  The following conventions have been used when calling ILAENV from the
   *  LAPACK routines:
   *  1)  OPTS is a concatenation of all of the character options to
   *      subroutine NAME, in the same order that they appear in the
   *      argument list for NAME, even if they are not used in determining
   *      the value of the parameter specified by ISPEC.
   *  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
   *      that they appear in the argument list for NAME.  N1 is used
   *      first, N2 second, and so on, and unused problem dimensions are
   *      passed a value of -1.
   *  3)  The parameter value returned by ILAENV is checked for validity in
   *      the calling subroutine.  For example, ILAENV is used to retrieve
   *      the optimal blocksize for STRTRI as follows:
   *
   *      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
   *      IF( NB.LE.1 ) NB = MAX( 1, N )
   *
   *  =====================================================================
  \*/

  inline
  integer
  ilaenv(
    integer         ISPEC,
    character const NAME[],
    character const OPTS[],
    integer         N1,
    integer         N2,
    integer         N3,
    integer         N4
  ) {
    #ifdef LAPACK_WRAPPER_USE_ACCELERATE
    return CLAPACKNAME(ilaenv)(
      &ISPEC,
      const_cast<character*>(NAME),
      const_cast<character*>(OPTS),
      &N1, &N2, &N3, &N4
    );
    #elif defined(LAPACK_WRAPPER_USE_LAPACK) || \
          defined(LAPACK_WRAPPER_USE_ATLAS)
    size_t len_NAME = strlen( NAME );
    size_t len_OPTS = strlen( OPTS );
    return LAPACK_F77NAME(ilaenv)(
      &ISPEC,
      const_cast<character*>(NAME),
      const_cast<character*>(OPTS),
      &N1, &N2, &N3, &N4, &len_NAME, &len_OPTS
    );
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
          defined(LAPACK_WRAPPER_USE_BLASFEO)
    size_t len_NAME = strlen( NAME );
    size_t len_OPTS = strlen( OPTS );
    return BLASFUNC(ilaenv)(
      &ISPEC,
      const_cast<character*>(NAME),
      const_cast<character*>(OPTS),
      &N1, &N2, &N3, &N4, &len_NAME, &len_OPTS
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    return ::ilaenv( &ISPEC, NAME, OPTS, &N1, &N2, &N3, &N4 );
    #else
    #error "lapack_wrapper undefined mapping!"
    return 0;
    #endif
  }

  /*\
   *  _        _    ___ ____ _
   * | |      / \  |_ _/ ___/ |
   * | |     / _ \  | | |   | |
   * | |___ / ___ \ | | |___| |
   * |_____/_/   \_\___\____|_|
   *
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {

    void
    LAPACK_F77NAME(slaic1)(
      integer const * JOB,
      integer const * J,
      real    const   X[],
      real    const * SEST,
      real    const   W[],
      real    const * GAMMA,
      real          * SESTPR,
      real          * S,
      real          * C
    );

    void
    LAPACK_F77NAME(dlaic1)(
      integer    const * JOB,
      integer    const * J,
      doublereal const   X[],
      doublereal const * SEST,
      doublereal const   W[],
      doublereal const * GAMMA,
      doublereal       * SESTPR,
      doublereal       * S,
      doublereal       * C
    );
  }
  #endif

  #if defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
      defined(LAPACK_WRAPPER_USE_BLASFEO)
  extern "C" {
    void
    BLASFUNC(slaic1)(
      integer const * JOB,
      integer const * J,
      real    const   X[],
      real    const * SEST,
      real    const   W[],
      real    const * GAMMA,
      real          * SESTPR,
      real          * S,
      real          * C
    );

    void
    BLASFUNC(dlaic1)(
      integer    const * JOB,
      integer    const * J,
      doublereal const   X[],
      doublereal const * SEST,
      doublereal const   W[],
      doublereal const * GAMMA,
      doublereal       * SESTPR,
      doublereal       * S,
      doublereal       * C
    );
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DLAIC1 applies one step of incremental condition estimation in
   *  its simplest version:
   *
   *  Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j
   *  lower triangular matrix L, such that
   *           twonorm(L*x) = sest
   *  Then DLAIC1 computes sestpr, s, c such that
   *  the vector
   *                  [ s*x ]
   *           xhat = [  c  ]
   *  is an approximate singular vector of
   *                  [ L     0  ]
   *           Lhat = [ w' gamma ]
   *  in the sense that
   *           twonorm(Lhat*xhat) = sestpr.
   *
   *  Depending on JOB, an estimate for the largest or smallest singular
   *  value is computed.
   *
   *  Note that [s c]' and sestpr**2 is an eigenpair of the system
   *
   *      diag(sest*sest, 0) + [alpha  gamma] * [ alpha ]
   *                                            [ gamma ]
   *
   *  where  alpha =  x'*w.
   *
   *  Arguments
   *  =========
   *
   *  JOB     (input) INTEGER
   *          = 1: an estimate for the largest singular value is computed.
   *          = 2: an estimate for the smallest singular value is computed.
   *
   *  J       (input) INTEGER
   *          Length of X and W
   *
   *  X       (input) DOUBLE PRECISION array, dimension (J)
   *          The j-vector x.
   *
   *  SEST    (input) DOUBLE PRECISION
   *          Estimated singular value of j by j matrix L
   *
   *  W       (input) DOUBLE PRECISION array, dimension (J)
   *          The j-vector w.
   *
   *  GAMMA   (input) DOUBLE PRECISION
   *          The diagonal element gamma.
   *
   *  SESTPR  (output) DOUBLE PRECISION
   *          Estimated singular value of (j+1) by (j+1) matrix Lhat.
   *
   *  S       (output) DOUBLE PRECISION
   *          Sine needed in forming xhat.
   *
   *  C       (output) DOUBLE PRECISION
   *          Cosine needed in forming xhat.
   *
  \*/

  inline
  void
  laic1(
    integer    JOB,
    integer    J,
    real const X[],
    real       SEST,
    real const W[],
    real       GAMMA,
    real     & SESTPR,
    real     & S,
    real     & C
  )
  #if defined(LAPACK_WRAPPER_USE_ACCELERATE)
  { CLAPACKNAME(slaic1)(
      &JOB, &J,
      const_cast<real*>(X), &SEST,
      const_cast<real*>(W), &GAMMA,
      &SESTPR, &S, &C
    );
  }
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)
  { LAPACK_F77NAME(slaic1)(
      &JOB, &J,
      const_cast<real*>(X), &SEST,
      const_cast<real*>(W), &GAMMA,
      &SESTPR, &S, &C
    );
  }
  #elif defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_BLASFEO)
  { BLASFUNC(slaic1)(
      &JOB, &J,
      const_cast<real*>(X), &SEST,
      const_cast<real*>(W), &GAMMA,
      &SESTPR, &S, &C
    );
  }
  #elif defined(LAPACK_WRAPPER_USE_ATLAS)
  { LAPACK_F77NAME(slaic1)(
      &JOB, &J,
      const_cast<real*>(X), &SEST,
      const_cast<real*>(W), &GAMMA,
      &SESTPR, &S, &C
    );
  }
  #elif defined(LAPACK_WRAPPER_USE_MKL)
  { slaic1( &JOB, &J, X, &SEST, W, &GAMMA, &SESTPR, &S, &C  ); }
  #else
  #error "lapack_wrapper undefined mapping!"
  #endif

  inline
  void
  laic1(
    integer          JOB,
    integer          J,
    doublereal const X[],
    doublereal       SEST,
    doublereal const W[],
    doublereal       GAMMA,
    doublereal     & SESTPR,
    doublereal     & S,
    doublereal     & C
  )
  #if defined(LAPACK_WRAPPER_USE_ACCELERATE)
  { CLAPACKNAME(dlaic1)(
      &JOB, &J,
      const_cast<doublereal*>(X), &SEST,
      const_cast<doublereal*>(W), &GAMMA,
      &SESTPR, &S, &C
    );
  }
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)
  { LAPACK_F77NAME(dlaic1)(
      &JOB, &J,
      const_cast<doublereal*>(X), &SEST,
      const_cast<doublereal*>(W), &GAMMA,
      &SESTPR, &S, &C
    );
  }
  #elif defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_BLASFEO)
  { BLASFUNC(dlaic1)(
      &JOB, &J,
      const_cast<doublereal*>(X), &SEST,
      const_cast<doublereal*>(W), &GAMMA,
      &SESTPR, &S, &C
    );
  }
  #elif defined(LAPACK_WRAPPER_USE_ATLAS)
  { LAPACK_F77NAME(dlaic1)(
      &JOB, &J,
      const_cast<doublereal*>(X), &SEST,
      const_cast<doublereal*>(W), &GAMMA,
      &SESTPR, &S, &C
    );
  }
  #elif defined(LAPACK_WRAPPER_USE_MKL)
  { dlaic1( &JOB, &J, X, &SEST, W, &GAMMA, &SESTPR, &S, &C ); }
  #else
  #error "lapack_wrapper undefined mapping!"
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DLARNV returns a vector of n random real numbers from a uniform or
   *  normal distribution.
   *
   *  Arguments
   *  =========
   *
   *  IDIST   (input) INTEGER
   *          Specifies the distribution of the random numbers:
   *          = 1:  uniform (0,1)
   *          = 2:  uniform (-1,1)
   *          = 3:  normal (0,1)
   *
   *  ISEED   (input/output) INTEGER array, dimension (4)
   *          On entry, the seed of the random number generator; the array
   *          elements must be between 0 and 4095, and ISEED(4) must be
   *          odd.
   *          On exit, the seed is updated.
   *
   *  N       (input) INTEGER
   *          The number of random numbers to be generated.
   *
   *  X       (output) DOUBLE PRECISION array, dimension (N)
   *          The generated random numbers.
  \*/

  // use standard Lapack routine
  #if defined(LAPACK_WRAPPER_USE_LAPACK) || defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {

    integer
    LAPACK_F77NAME(slarnv)(
      integer * IDIST,
      integer * ISEED,
      integer * N,
      real      X[]
    );

    integer
    LAPACK_F77NAME(dlarnv)(
      integer  * IDIST,
      integer  * ISEED,
      integer  * N,
      doublereal X[]
    );
  }
  #endif

  inline
  integer
  larnv(
    integer IDIST,
    integer ISEED[4],
    integer N,
    real    X[]
  ) {
    #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)  || \
        defined(LAPACK_WRAPPER_USE_BLASFEO)
    LAPACK_F77NAME(slarnv)( &IDIST, ISEED, &N, X ); // return void in openblas
    return 0;
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
    LAPACK_slarnv( &IDIST, ISEED, &N, X ); // return void in openblas
    return 0;
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    slarnv( &IDIST, ISEED, &N, X ); // return void in openblas
    return 0;
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    return CLAPACKNAME(slarnv)( &IDIST, ISEED, &N, X );
    #else
    return LAPACKNAME(slarnv)( &IDIST, ISEED, &N, X );
    #endif
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline
  integer
  larnv(
    integer    IDIST,
    integer    ISEED[4],
    integer    N,
    doublereal X[]
  ) {
    #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)  || \
        defined(LAPACK_WRAPPER_USE_BLASFEO)
    LAPACK_F77NAME(dlarnv)( &IDIST, ISEED, &N, X ); // return void in openblas
    return 0;
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
    LAPACK_dlarnv( &IDIST, ISEED, &N, X ); // return void in openblas
    return 0;
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dlarnv( &IDIST, ISEED, &N, X ); // return void in openblas
    return 0;
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    return CLAPACKNAME(dlarnv)( &IDIST, ISEED, &N, X );
    #else
    return LAPACKNAME(dlarnv)( &IDIST, ISEED, &N, X );
    #endif
  }

  /*\
   *  Purpose
   *  =======
   *
   *  DLARGE pre- and post-multiplies a real general n by n matrix A
   *  with a random orthogonal matrix: A = U*D*U'.
   *
   *  Arguments
   *  =========
   *
   *  N       (input) INTEGER
   *          The order of the matrix A.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the original n by n matrix A.
   *          On exit, A is overwritten by U*A*U' for some random
   *          orthogonal matrix U.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= N.
   *
   *  ISEED   (input/output) INTEGER array, dimension (4)
   *          On entry, the seed of the random number generator; the array
   *          elements must be between 0 and 4095, and ISEED(4) must be
   *          odd.
   *          On exit, the seed is updated.
   *
   *  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)
   *
   *  INFO    (output) INTEGER
   *          = 0: successful exit
   *          < 0: if INFO = -i, the i-th argument had an illegal value
   *
   *  =====================================================================
  \*/

  template <typename T>
  inline
  integer
  large(
    integer N,
    T       A[],
    integer LDA,
    integer ISEED[4],
    T       WORK[],
    integer LWORK
  );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}

#include "code/blas.hxx"
#include "code/general.hxx"
#include "code/general_qr.hxx"
#include "code/triangular.hxx"
#include "code/tridiagonal.hxx"
#include "code/banded.hxx"
#include "code/symmetric.hxx"
#include "code/sparse.hxx"
#include "code/wrapper.hxx"

namespace lapack_wrapper {

  inline
  void
  getranspose(
    integer    M,
    integer    N,
    real const A[],
    integer    LDA,
    real       B[],
    integer    LDB
  ) {
    for ( integer k = 0; k < N; ++k )
      copy( M, A + k*LDA, 1, B+k, LDB );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline
  void
  getranspose(
    integer          M,
    integer          N,
    doublereal const A[],
    integer          LDA,
    doublereal       B[],
    integer          LDB
  ) {
    for ( integer k = 0; k < N; ++k )
      copy( M, A + k*LDA, 1, B+k, LDB );
  }
}

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#endif

///
/// eof: lapack_wrapper.hh
///
