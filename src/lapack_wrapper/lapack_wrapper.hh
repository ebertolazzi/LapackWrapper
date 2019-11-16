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

#ifndef LAPACK_WRAPPER_HH
#define LAPACK_WRAPPER_HH

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

#ifdef LAPACK_WRAPPER_OS_WINDOWS
  // se in windows includo PRIMA windows.h per evitare conflitti!
  #ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN
  #endif
  #include <windows.h>
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
  // windows define max and min!!!
  #ifdef max
    #undef max
  #endif
  #ifdef min
    #undef min
  #endif

  // select LAPACK dll/lib for windows
  #if _MSC_VER >= 1900
    #pragma comment(lib, "legacy_stdio_definitions.lib")
  #endif

  #ifndef MINGW
    #if defined(LAPACK_WRAPPER_USE_LAPACK2)
      #if defined(_DEBUG) || defined(DEBUG)
        #ifdef LAPACK_WRAPPER_ARCH64
          #pragma comment(lib, "cbia.lib.blas.dyn.dbg.x64.12.lib")
          #pragma comment(lib, "cbia.lib.lapack.dyn.dbg.x64.12.lib")
        #else
          #pragma comment(lib, "cbia.lib.blas.dyn.dbg.x86.12.lib")
          #pragma comment(lib, "cbia.lib.lapack.dyn.dbg.x86.12.lib")
        #endif
      #else
        #ifdef LAPACK_WRAPPER_ARCH64
          #pragma comment(lib, "cbia.lib.blas.dyn.rel.x64.12.lib")
          #pragma comment(lib, "cbia.lib.lapack.dyn.rel.x64.12.lib")
        #else
          #pragma comment(lib, "cbia.lib.blas.dyn.rel.x86.12.lib")
          #pragma comment(lib, "cbia.lib.lapack.dyn.rel.x86.12.lib")
        #endif
      #endif
    #elif defined(LAPACK_WRAPPER_USE_LAPACK)
      #if defined(_DEBUG) || defined(DEBUG)
        #ifdef LAPACK_WRAPPER_ARCH64
          #pragma comment(lib, "blas_win64_MTd.lib")
          #pragma comment(lib, "lapack_win64_MTd.lib")
        #else
          #pragma comment(lib, "blas_win32_MTd.lib")
          #pragma comment(lib, "lapack_win32_MTd.lib")
        #endif
      #else
        #ifdef LAPACK_WRAPPER_ARCH64
          #pragma comment(lib, "blas_win64_MT.lib")
          #pragma comment(lib, "lapack_win64_MT.lib")
        #else
          #pragma comment(lib, "blas_win32_MT.lib")
          #pragma comment(lib, "lapack_win32_MT.lib")
        #endif
      #endif
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
          defined(LAPACK_WRAPPER_USE_BLASFEO)
      // no debug version
      #ifdef LAPACK_WRAPPER_ARCH64
        #pragma comment(lib, "libopenblas_x64.lib")
      #else
        #pragma comment(lib, "libopenblas_x86.lib")
      #endif
    #elif defined(LAPACK_WRAPPER_USE_MKL)
      #ifdef LAPACK_WRAPPER_ARCH64
        #pragma comment(lib, "mkl_intel_lp64_dll.lib")
      #else
        #pragma comment(lib, "mkl_intel_c_dll.lib")
      #endif
      #pragma comment(lib, "mkl_sequential_dll.lib")
      #pragma comment(lib, "mkl_core_dll.lib")
    #else
      #error "Only standard Lapack, Openblas, and MKL are supported for WINDOWS!"
    #endif
  #endif
#else
  #include <cstdint>
#endif

// set macro for used LAPACK name

#if defined(LAPACK_WRAPPER_USE_ACCELERATE)
  #ifdef LAPACK_WRAPPER_ARCH64
    #define LAPACK_WRAPPER_LAPACK_NAME "Accelerate OSX (x64)"
  #else
    #define LAPACK_WRAPPER_LAPACK_NAME "Accelerate OSX (x86)"
  #endif
#elif defined(LAPACK_WRAPPER_USE_LAPACK)
  #if defined(_DEBUG) || defined(DEBUG)
    #ifdef LAPACK_WRAPPER_ARCH64
      #define LAPACK_WRAPPER_LAPACK_NAME "Standard BLAS/LAPACK (x64,DEBUG)"
    #else
      #define LAPACK_WRAPPER_LAPACK_NAME "Standard BLAS/LAPACK (x86,DEBUG)"
    #endif
  #else
    #ifdef LAPACK_WRAPPER_ARCH64
      #define LAPACK_WRAPPER_LAPACK_NAME "Standard BLAS/LAPACK (x64)"
    #else
      #define LAPACK_WRAPPER_LAPACK_NAME "Standard BLAS/LAPACK (x86)"
    #endif
  #endif
#elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
  #ifdef LAPACK_WRAPPER_ARCH64
    #define LAPACK_WRAPPER_LAPACK_NAME "OPENBLAS (x64,RELEASE)"
  #else
    #define LAPACK_WRAPPER_LAPACK_NAME "OPENBLAS (x86,RELEASE)"
  #endif
#elif defined(LAPACK_WRAPPER_USE_BLASFEO)
  #ifdef LAPACK_WRAPPER_ARCH64
    #define LAPACK_WRAPPER_LAPACK_NAME "BLASFEO (x64,RELEASE)"
  #else
    #define LAPACK_WRAPPER_LAPACK_NAME "BLASFEO (x86,RELEASE)"
  #endif
#elif defined(LAPACK_WRAPPER_USE_MKL)
  #ifdef LAPACK_WRAPPER_ARCH64
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
  #ifdef LAPACK_WRAPPER_OS_OSX
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
  #ifdef LAPACK_WRAPPER_OS_OSX
    #ifdef __clang__
      #pragma clang diagnostic ignored "-Wreserved-id-macro"
    #endif
    #define __STDC_VERSION__ __STDC__
  #endif

  #include <complex>
  #define FORCE_OPENBLAS_COMPLEX_STRUCT
  #define lapack_complex_float  std::complex<float>
  #define lapack_complex_double std::complex<double>

  #ifdef LAPACK_WRAPPER_OS_WINDOWS
    // on windows the defaul is too use local openblas
    #ifdef LAPACK_WRAPPER_USE_SYSTEM_OPENBLAS
      #include <openblas/cblas.h>
      #include <openblas/lapacke.h>
    #else
      #ifdef LAPACK_WRAPPER_ARCH64
        #include "../openblas/x64/cblas.h"
        #include "../openblas/x64/lapacke.h"
      #else
        #include "../openblas/x86/cblas.h"
        #include "../openblas/x86/lapacke.h"
      #endif
    #endif
  #elif defined(LAPACK_WRAPPER_OS_LINUX)
    // on linux/osx the defaul is too use SYSTEM openblas
    #ifdef LAPACK_WRAPPER_DO_NOT_USE_SYSTEM_OPENBLAS
      #include "../openblas/cblas.h"
      #include "../openblas/lapacke.h"
    #else
      // must use subdirectory to avoid conflict with atlas/lapack...
      #include <openblas/cblas.h>
      #include <lapacke.h>
    #endif
  #else
    // on linux/osx the defaul is too use SYSTEM openblas
    #ifdef LAPACK_WRAPPER_DO_NOT_USE_SYSTEM_OPENBLAS
      #include "../openblas/cblas.h"
      #include "../openblas/lapacke.h"
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

  typedef enum {
    NO_TRANSPOSE        = 0,
    TRANSPOSE           = 1,
    CONJUGATE_TRANSPOSE = 2
  } Transposition;

  typedef enum {
    UPPER = 0,
    LOWER = 1
  } ULselect;

  typedef enum {
    UNIT     = 0,
    NON_UNIT = 1
  } DiagonalType;

  typedef enum {
    LEFT  = 0,
    RIGHT = 1
  } SideMultiply;

  typedef enum {
    NO_BALANCE        = 0, // 'N'
    PERMUTE_ONLY      = 1, // 'P'
    SCALE_ONLY        = 2, // 'S'
    PERMUTE_AND_SCALE = 3  // 'B'
  } BalanceType;

  typedef enum {
    ALL     = 0, // 'A'
    REDUCED = 1, // 'S'
    INPLACE = 2, // 'O'
    NO_JOB  = 3  // 'N'
  } JobType;

  typedef enum {
    NONE                         = 0, // 'N'
    EIGENVALUES_ONLY             = 1, // 'E'
    EIGENVECTORS_ONLY            = 2, // 'V'
    EIGENVALUES_AND_EIGENVECTORS = 3  // 'B'
  } SenseType;

  typedef enum {
    FORWARD  = 0,
    BACKWARD = 1
  } DirectionType;

  typedef enum {
    COLUMNWISE = 0,
    ROWWISE    = 1
  } StorageType;

  typedef enum {
    FULL_MATRIX             = 0,
    LOWER_TRIANGULAR_MATRIX = 1,
    UPPER_TRIANGULAR_MATRIX = 2,
    HESSENBERG_MATRIX       = 3,
    BANDED_MATRIX           = 4
  } MatrixType;

  typedef enum {
    NO_EQUILIBRATE      = 0,
    EQUILIBRATE_ROWS    = 1,
    EQUILIBRATE_COLUMNS = 2,
    EQUILIBRATE_BOTH    = 3
  } EquilibrationType;

  #if defined(LAPACK_WRAPPER_USE_ACCELERATE)
    typedef __CLPK_integer    integer;
    typedef __CLPK_real       real;
    typedef __CLPK_doublereal doublereal;
    typedef char              character;
    typedef doublereal        return_precision;
  #elif defined(LAPACK_WRAPPER_USE_ATLAS)
    typedef int        integer;
    typedef float      real;
    typedef double     doublereal;
    typedef char       character;
    typedef doublereal return_precision;
  #elif defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_BLASFEO)
    typedef blasint    integer;
    typedef float      real;
    typedef double     doublereal;
    typedef char       character;
    typedef double     return_precision;
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)
    typedef lapack_int integer;
    typedef float      real;
    typedef double     doublereal;
    typedef char       character;
    typedef doublereal return_precision;
  #elif defined(LAPACK_WRAPPER_USE_MKL)
    typedef MKL_INT integer;
    typedef float   real;
    typedef double  doublereal;
    typedef char    character;
    typedef double  return_precision;
  #else
    #error "You must select the linear algebra packages used!"
  #endif

  #if defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)      || \
      defined(LAPACK_WRAPPER_USE_OPENBLAS)   || \
      defined(LAPACK_WRAPPER_USE_BLASFEO)
    extern CBLAS_TRANSPOSE trans_cblas[3];
    extern CBLAS_UPLO      uplo_cblas[2];
    extern CBLAS_DIAG      diag_cblas[2];
    extern CBLAS_SIDE      side_cblas[2];
  #endif

  extern character const *trans_blas[3];
  extern character const *uplo_blas[2];
  extern character const *diag_blas[2];
  extern character const *side_blas[2];

  extern character const *balance_blas[4];
  extern character const *job_blas[4];
  extern character const *sense_blas[4];
  extern character const *direct_blas[2];
  extern character const *store_blas[2];
  extern character const *mtype_blas[5];
  extern character const *equilibrate_blas[4];

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
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)    || \
        defined(LAPACK_WRAPPER_USE_BLASFEO)
    LAPACK_F77NAME(slarnv)( &IDIST, ISEED, &N, X ); // return void in openblas
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
  lanrv(
    integer    IDIST,
    integer    ISEED[4],
    integer    N,
    doublereal X[]
  ) {
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)    || \
        defined(LAPACK_WRAPPER_USE_BLASFEO)
    LAPACK_F77NAME(dlarnv)( &IDIST, ISEED, &N, X ); // return void in openblas
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
  ) {
    LW_ASSERT(
      LWORK >= 2*N,
      "large, LWORK = {} must be >= {}\n", LWORK, 2*N
    );
    // Test the input arguments
    if ( N < 0 ) return -1;
    if ( LDA < max_index( 1, N ) ) return -3;
    // pre- and post-multiply A by random orthogonal matrix
    for ( integer i = N-1; i >= 0; --i ) {
      // generate random reflection
      integer Ni = N-i;
      lanrv( 3, ISEED, N, WORK );
      T WN = nrm2( Ni, WORK, 1 );
      T WA = WN;
      if ( WORK[0] < 0 ) WA = -WN;
      T TAU = 0;
      if ( WN > 0 ) {
        T WB = WORK[0] + WA;
        scal( Ni-1, 1.0/WB, WORK+1, 1 );
        WORK[0] = 1;
        TAU = WB / WA;
      }
      // multiply A(i:n,1:n) by random reflection from the left
      gemv( TRANSPOSE, Ni, N, 1.0, A+i, LDA, WORK, 1, 0.0, WORK+N, 1 );
      ger( Ni, N, -TAU, WORK, 1, WORK+N, 1, A+i, LDA );
      // multiply A(1:n,i:n) by random reflection from the right
      gemv( NO_TRANSPOSE, N, Ni, 1.0, A+i*LDA, LDA, WORK, 1, 0.0, WORK+N, 1 );
      ger( N, Ni, -TAU, WORK+N, 1, WORK, 1, A+i*LDA, LDA );
    }
    return 0;
  }

  /*\
   *  Purpose
   *  =======
   *
   *  SLAROR pre- or post-multiplies an M by N matrix A by a random
   *  orthogonal matrix U, overwriting A.  A may optionally be initialized
   *  to the identity matrix before multiplying by U.  U is generated using
   *  the method of G.W. Stewart (SIAM J. Numer. Anal. 17, 1980, 403-409).
   *
   *  Arguments
   *  =========
   *
   *  SIDE    (input) CHARACTER*1
   *          Specifies whether A is multiplied on the left or right by U.
   *          = 'L':         Multiply A on the left (premultiply) by U
   *          = 'R':         Multiply A on the right (postmultiply) by U'
   *          = 'C' or 'T':  Multiply A on the left by U and the right
   *                          by U' (Here, U' means U-transpose.)
   *
   *  INIT    (input) CHARACTER*1
   *          Specifies whether or not A should be initialized to the
   *          identity matrix.
   *          = 'I':  Initialize A to (a section of) the identity matrix
   *                   before applying U.
   *          = 'N':  No initialization.  Apply U to the input matrix A.
   *
   *          INIT = 'I' may be used to generate square or rectangular
   *          orthogonal matrices:
   *
   *          For M = N and SIDE = 'L' or 'R', the rows will be orthogonal
   *          to each other, as will the columns.
   *
   *          If M < N, SIDE = 'R' produces a dense matrix whose rows are
   *          orthogonal and whose columns are not, while SIDE = 'L'
   *          produces a matrix whose rows are orthogonal, and whose first
   *          M columns are orthogonal, and whose remaining columns are
   *          zero.
   *
   *          If M > N, SIDE = 'L' produces a dense matrix whose columns
   *          are orthogonal and whose rows are not, while SIDE = 'R'
   *          produces a matrix whose columns are orthogonal, and whose
   *          first M rows are orthogonal, and whose remaining rows are
   *          zero.
   *
   *  M       (input) INTEGER
   *          The number of rows of A.
   *
   *  N       (input) INTEGER
   *          The number of columns of A.
   *
   *  A       (input/output) REAL array, dimension (LDA, N)
   *          On entry, the array A.
   *          On exit, overwritten by U A ( if SIDE = 'L' ),
   *           or by A U ( if SIDE = 'R' ),
   *           or by U A U' ( if SIDE = 'C' or 'T').
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  ISEED   (input/output) INTEGER array, dimension (4)
   *          On entry ISEED specifies the seed of the random number
   *          generator. The array elements should be between 0 and 4095;
   *          if not they will be reduced mod 4096.  Also, ISEED(4) must
   *          be odd.  The random number generator uses a linear
   *          congruential sequence limited to small integers, and so
   *          should produce machine independent random numbers. The
   *          values of ISEED are changed on exit, and can be used in the
   *          next call to SLAROR to continue the same random number
   *          sequence.
   *
   *  X       (workspace) REAL array, dimension (3*MAX( M, N ))
   *          Workspace of length
   *              2*M + N if SIDE = 'L',
   *              2*N + M if SIDE = 'R',
   *              3*N     if SIDE = 'C' or 'T'.
   *
   *  INFO    (output) INTEGER
   *          An error flag.  It is set to:
   *          = 0:  normal return
   *          < 0:  if INFO = -k, the k-th argument had an illegal value
   *          = 1:  if the random numbers generated by SLARND are bad.
   *
   *  =====================================================================
  \*/

  template <typename T>
  inline
  integer
  laror(
    SideMultiply const & LR,
    bool                 init,
    integer              M,
    integer              N,
    T                    A[],
    integer              LDA,
    integer              ISEED[4],
    T                    WORK[],
    integer              LWORK
  ) {

    if ( LR == LEFT ) {
      LW_ASSERT(
        LWORK >= 2*M + N, "laror, LWORK = {} must be >= {}\n", LWORK, 2*M + N
      );
    } else if ( LR == RIGHT ) {
      LW_ASSERT(
        LWORK >= 2*N + M, "laror, LWORK = {} must be >= {}\n", LWORK, 2*N + M
      );
    } else {
      return -1;
    }

    if ( N == 0 || M == 0 ) return 0;

    // Check for argument errors.
    if ( M   < 0 ) return -3;
    if ( N   < 0 ) return -4;
    if ( LDA < M ) return -6;

    // Initialize A to the identity matrix if desired
    if ( init ) geid( M, N, A, LDA );

    LW_ASSERT(
      LWORK >= 2*N, "large, LWORK = {} must be >= {}\n", LWORK, 2*N
    );
    // Test the input arguments
    if ( N < 0 ) return -1;
    if ( LDA < max_index( 1, N ) ) return -3;

    if ( LR == LEFT ) {
      // pre-multiply A by random orthogonal matrix
      for ( integer i = M-1; i >= 0; --i ) {
        // generate random reflection
        integer Mi = M-i;
        lanrv( 1, ISEED, M, WORK );
        T WN = nrm2( Mi, WORK, 1 );
        T WA = WN;
        if ( WORK[0] < 0 ) WA = -WN;
        T TAU = 0;
        if ( WN > 0 ) {
          T WB = WORK[0] + WA;
          scal( Mi-1, 1.0/WB, WORK+1, 1 );
          WORK[0] = 1;
          TAU = WB / WA;
        }
        // multiply A(i:n,1:n) by random reflection from the left
        gemv( TRANSPOSE, Mi, M, 1.0, A+i, LDA, WORK, 1, 0.0, WORK+M, 1 );
        ger( Mi, N, -TAU, WORK, 1, WORK+M, 1, A+i, LDA );
      }
    } else {
      // post-multiply A by random orthogonal matrix
      for ( integer i = N-1; i >= 0; --i ) {
        // generate random reflection
        integer Ni = N-i;
        lanrv( 1, ISEED, N, WORK );
        T WN = nrm2( Ni, WORK, 1 );
        T WA = WN;
        if ( WORK[0] < 0 ) WA = -WN;
        T TAU = 0;
        if ( WN > 0 ) {
          T WB = WORK[0] + WA;
          scal( Ni-1, 1.0/WB, WORK+1, 1 );
          WORK[0] = 1;
          TAU = WB / WA;
        }
        // multiply A(1:n,i:n) by random reflection from the right
        gemv( NO_TRANSPOSE, N, Ni, 1.0, A+i*LDA, LDA, WORK, 1, 0.0, WORK+N, 1 );
        ger( N, Ni, -TAU, WORK+N, 1, WORK, 1, A+i*LDA, LDA );
      }
    }

    return 0;
  }

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
