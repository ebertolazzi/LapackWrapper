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
/// file: symmetric.hxx
///

/*
//   ____                                 _        _
//  / ___| _   _ _ __ ___  _ __ ___   ___| |_ _ __(_) ___
//  \___ \| | | | '_ ` _ \| '_ ` _ \ / _ \ __| '__| |/ __|
//   ___) | |_| | | | | | | | | | | |  __/ |_| |  | | (__
//  |____/ \__, |_| |_| |_|_| |_| |_|\___|\__|_|  |_|\___|
//         |___/
//   __  __       _        _
//  |  \/  | __ _| |_ _ __(_)_  __
//  | |\/| |/ _` | __| '__| \ \/ /
//  | |  | | (_| | |_| |  | |>  <
//  |_|  |_|\__,_|\__|_|  |_/_/\_\
*/

namespace lapack_wrapper {

  /*
  //   ___ _   _ _ __
  //  / __| | | | '__|
  //  \__ \ |_| | |
  //  |___/\__, |_|
  //       |___/
  */
  /*\
   *
   *  Purpose
   *  =======
   *
   *  DSYR   performs the symmetric rank 1 operation
   *
   *     A := alpha*x*x' + A,
   *
   *  where alpha is a real scalar, x is an n element vector and A is an
   *  n by n symmetric matrix.
   *
   *  Arguments
   *  ==========
   *
   *  UPLO   - CHARACTER*1.
   *           On entry, UPLO specifies whether the upper or lower
   *           triangular part of the array A is to be referenced as
   *           follows:
   *
   *              UPLO = 'U' or 'u'   Only the upper triangular part of A
   *                                  is to be referenced.
   *
   *              UPLO = 'L' or 'l'   Only the lower triangular part of A
   *                                  is to be referenced.
   *
   *           Unchanged on exit.
   *
   *  N      - INTEGER.
   *           On entry, N specifies the order of the matrix A.
   *           N must be at least zero.
   *           Unchanged on exit.
   *
   *  ALPHA  - DOUBLE PRECISION.
   *           On entry, ALPHA specifies the scalar alpha.
   *           Unchanged on exit.
   *
   *  X      - DOUBLE PRECISION array of dimension at least
   *           ( 1 + ( n - 1 )*abs( INCX ) ).
   *           Before entry, the incremented array X must contain the n
   *           element vector x.
   *           Unchanged on exit.
   *
   *  INCX   - INTEGER.
   *           On entry, INCX specifies the increment for the elements of
   *           X. INCX must not be zero.
   *           Unchanged on exit.
   *
   *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
   *           Before entry with  UPLO = 'U' or 'u', the leading n by n
   *           upper triangular part of the array A must contain the upper
   *           triangular part of the symmetric matrix and the strictly
   *           lower triangular part of A is not referenced. On exit, the
   *           upper triangular part of the array A is overwritten by the
   *           upper triangular part of the updated matrix.
   *           Before entry with UPLO = 'L' or 'l', the leading n by n
   *           lower triangular part of the array A must contain the lower
   *           triangular part of the symmetric matrix and the strictly
   *           upper triangular part of A is not referenced. On exit, the
   *           lower triangular part of the array A is overwritten by the
   *           lower triangular part of the updated matrix.
   *
   *  LDA    - INTEGER.
   *           On entry, LDA specifies the first dimension of A as declared
   *           in the calling (sub) program. LDA must be at least
   *           max( 1, n ).
   *           Unchanged on exit.
   *
   *
   *  Level 2 Blas routine.
   *
   *  -- Written on 22-October-1986.
   *     Jack Dongarra, Argonne National Lab.
   *     Jeremy Du Croz, Nag Central Office.
   *     Sven Hammarling, Nag Central Office.
   *     Richard Hanson, Sandia National Labs.
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  extern "C" {
    void
    BLASFUNC(ssyr)(
      char    const * UPLO,
      integer const * N,
      real    const * ALPHA,
      real    const   X[],
      integer const * INCX,
      real            A[],
      integer const * LDA
    );

    void
    BLASFUNC(dsyr)(
      char       const * UPLO,
      integer    const * N,
      doublereal const * ALPHA,
      doublereal const   X[],
      integer    const * INCX,
      doublereal         A[],
      integer    const * LDA
    );
  }
  #endif

  inline
  void
  syr(
    ULselect const & UPLO,
    integer          N,
    real             ALPHA,
    real       const X[],
    integer          INCX,
    real             A[],
    integer          LDA
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { ssyr( uplo_blas[UPLO], &N, &ALPHA, X, &INCX, A, &LDA ); }
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)
  { BLASFUNC(ssyr)(
      const_cast<character*>(uplo_blas[UPLO]),
      &N, &ALPHA,
      const_cast<real*>(X), &INCX,
      A, &LDA
    );
  }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { CBLASNAME(ssyr)( CblasColMajor, uplo_cblas[UPLO],
                     N, ALPHA, X, INCX, A, LDA ); }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  inline
  void
  syr(
    ULselect const & UPLO,
    integer          N,
    doublereal       ALPHA,
    doublereal const X[],
    integer          INCX,
    doublereal       A[],
    integer          LDA
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { dsyr( uplo_blas[UPLO], &N, &ALPHA, X, &INCX, A, &LDA ); }
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)
  { BLASFUNC(dsyr)(
      const_cast<character*>(uplo_blas[UPLO]),
      &N, &ALPHA,
      const_cast<doublereal*>(X), &INCX,
      A, &LDA
    );
  }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { CBLASNAME(dsyr)( CblasColMajor, uplo_cblas[UPLO],
                     N, ALPHA, X, INCX, A, LDA ); }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  /*
  //                 ____
  //   ___ _   _ _ _|___ \
  //  / __| | | | '__|__) |
  //  \__ \ |_| | |  / __/
  //  |___/\__, |_| |_____|
  //       |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DSYR2  performs the symmetric rank 2 operation
   *
   *     A := alpha*x*y' + alpha*y*x' + A,
   *
   *  where alpha is a scalar, x and y are n element vectors and A is an n
   *  by n symmetric matrix.
   *
   *  Arguments
   *  ==========
   *
   *  UPLO   - CHARACTER*1.
   *           On entry, UPLO specifies whether the upper or lower
   *           triangular part of the array A is to be referenced as
   *           follows:
   *
   *              UPLO = 'U' or 'u'   Only the upper triangular part of A
   *                                  is to be referenced.
   *
   *              UPLO = 'L' or 'l'   Only the lower triangular part of A
   *                                  is to be referenced.
   *
   *           Unchanged on exit.
   *
   *  N      - INTEGER.
   *           On entry, N specifies the order of the matrix A.
   *           N must be at least zero.
   *           Unchanged on exit.
   *
   *  ALPHA  - DOUBLE PRECISION.
   *           On entry, ALPHA specifies the scalar alpha.
   *           Unchanged on exit.
   *
   *  X      - DOUBLE PRECISION array of dimension at least
   *           ( 1 + ( n - 1 )*abs( INCX ) ).
   *           Before entry, the incremented array X must contain the n
   *           element vector x.
   *           Unchanged on exit.
   *
   *  INCX   - INTEGER.
   *           On entry, INCX specifies the increment for the elements of
   *           X. INCX must not be zero.
   *           Unchanged on exit.
   *
   *  Y      - DOUBLE PRECISION array of dimension at least
   *           ( 1 + ( n - 1 )*abs( INCY ) ).
   *           Before entry, the incremented array Y must contain the n
   *           element vector y.
   *           Unchanged on exit.
   *
   *  INCY   - INTEGER.
   *           On entry, INCY specifies the increment for the elements of
   *           Y. INCY must not be zero.
   *           Unchanged on exit.
   *
   *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
   *           Before entry with  UPLO = 'U' or 'u', the leading n by n
   *           upper triangular part of the array A must contain the upper
   *           triangular part of the symmetric matrix and the strictly
   *           lower triangular part of A is not referenced. On exit, the
   *           upper triangular part of the array A is overwritten by the
   *           upper triangular part of the updated matrix.
   *           Before entry with UPLO = 'L' or 'l', the leading n by n
   *           lower triangular part of the array A must contain the lower
   *           triangular part of the symmetric matrix and the strictly
   *           upper triangular part of A is not referenced. On exit, the
   *           lower triangular part of the array A is overwritten by the
   *           lower triangular part of the updated matrix.
   *
   *  LDA    - INTEGER.
   *           On entry, LDA specifies the first dimension of A as declared
   *           in the calling (sub) program. LDA must be at least
   *           max( 1, n ).
   *           Unchanged on exit.
   *
   *
   *  Level 2 Blas routine.
   *
   *  -- Written on 22-October-1986.
   *     Jack Dongarra, Argonne National Lab.
   *     Jeremy Du Croz, Nag Central Office.
   *     Sven Hammarling, Nag Central Office.
   *     Richard Hanson, Sandia National Labs.
   *
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  extern "C" {
    void
    BLASFUNC(ssyr2)(
      char    const * UPLO,
      integer const * N,
      real    const * ALPHA,
      real    const   X[],
      integer const * INCX,
      real    const   Y[],
      integer const * INCY,
      real            A[],
      integer const * LDA
    );

    void
    BLASFUNC(dsyr2)(
      char       const * UPLO,
      integer    const * N,
      doublereal const * ALPHA,
      doublereal const   X[],
      integer    const * INCX,
      doublereal const   Y[],
      integer    const * INCY,
      doublereal         A[],
      integer    const * LDA
    );
  }
  #endif

  inline
  void
  syr2(
    ULselect const & UPLO,
    integer          N,
    real             ALPHA,
    real       const X[],
    integer          INCX,
    real       const Y[],
    integer          INCY,
    real             A[],
    integer          LDA
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { ssyr2( uplo_blas[UPLO], &N, &ALPHA, X, &INCX, Y, &INCY, A, &LDA ); }
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)
  { BLASFUNC(ssyr2)(
      const_cast<character*>(uplo_blas[UPLO]),
      &N, &ALPHA,
      const_cast<real*>(X), &INCX,
      const_cast<real*>(Y), &INCY,
      A, &LDA
    );
  }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { CBLASNAME(ssyr2)( CblasColMajor, uplo_cblas[UPLO],
                      N, ALPHA, X, INCX, Y, INCY, A, LDA ); }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  inline
  void
  syr2(
    ULselect const & UPLO,
    integer          N,
    doublereal       ALPHA,
    doublereal const X[],
    integer          INCX,
    doublereal const Y[],
    integer          INCY,
    doublereal       A[],
    integer          LDA
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { dsyr2( uplo_blas[UPLO], &N, &ALPHA, X, &INCX, Y, &INCY, A, &LDA ); }
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)
  { BLASFUNC(dsyr2)(
      const_cast<character*>(uplo_blas[UPLO]),
      &N, &ALPHA,
      const_cast<doublereal*>(X), &INCX,
      const_cast<doublereal*>(Y), &INCY,
      A, &LDA
    );
  }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { CBLASNAME(dsyr2)( CblasColMajor, uplo_cblas[UPLO],
                      N, ALPHA, X, INCX, Y, INCY, A, LDA ); }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  /*
  //   ___ _   _ _ __ _____   __
  //  / __| | | | '_ ` _ \ \ / /
  //  \__ \ |_| | | | | | \ V /
  //  |___/\__, |_| |_| |_|\_/
  //       |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DSYMV  performs the matrix-vector  operation
   *
   *     y := alpha*A*x + beta*y,
   *
   *  where alpha and beta are scalars, x and y are n element vectors and
   *  A is an n by n symmetric matrix.
   *
   *  Arguments
   *  ==========
   *
   *  UPLO   - CHARACTER*1.
   *           On entry, UPLO specifies whether the upper or lower
   *           triangular part of the array A is to be referenced as
   *           follows:
   *
   *              UPLO = 'U' or 'u'   Only the upper triangular part of A
   *                                  is to be referenced.
   *
   *              UPLO = 'L' or 'l'   Only the lower triangular part of A
   *                                  is to be referenced.
   *
   *           Unchanged on exit.
   *
   *  N      - INTEGER.
   *           On entry, N specifies the order of the matrix A.
   *           N must be at least zero.
   *           Unchanged on exit.
   *
   *  ALPHA  - DOUBLE PRECISION.
   *           On entry, ALPHA specifies the scalar alpha.
   *           Unchanged on exit.
   *
   *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
   *           Before entry with  UPLO = 'U' or 'u', the leading n by n
   *           upper triangular part of the array A must contain the upper
   *           triangular part of the symmetric matrix and the strictly
   *           lower triangular part of A is not referenced.
   *           Before entry with UPLO = 'L' or 'l', the leading n by n
   *           lower triangular part of the array A must contain the lower
   *           triangular part of the symmetric matrix and the strictly
   *           upper triangular part of A is not referenced.
   *           Unchanged on exit.
   *
   *  LDA    - INTEGER.
   *           On entry, LDA specifies the first dimension of A as declared
   *           in the calling (sub) program. LDA must be at least
   *           max( 1, n ).
   *           Unchanged on exit.
   *
   *  X      - DOUBLE PRECISION array of dimension at least
   *           ( 1 + ( n - 1 )*abs( INCX ) ).
   *           Before entry, the incremented array X must contain the n
   *           element vector x.
   *           Unchanged on exit.
   *
   *  INCX   - INTEGER.
   *           On entry, INCX specifies the increment for the elements of
   *           X. INCX must not be zero.
   *           Unchanged on exit.
   *
   *  BETA   - DOUBLE PRECISION.
   *           On entry, BETA specifies the scalar beta. When BETA is
   *           supplied as zero then Y need not be set on input.
   *           Unchanged on exit.
   *
   *  Y      - DOUBLE PRECISION array of dimension at least
   *           ( 1 + ( n - 1 )*abs( INCY ) ).
   *           Before entry, the incremented array Y must contain the n
   *           element vector y. On exit, Y is overwritten by the updated
   *           vector y.
   *
   *  INCY   - INTEGER.
   *           On entry, INCY specifies the increment for the elements of
   *           Y. INCY must not be zero.
   *           Unchanged on exit.
   *
   *
   *  Level 2 Blas routine.
   *
   *  -- Written on 22-October-1986.
   *     Jack Dongarra, Argonne National Lab.
   *     Jeremy Du Croz, Nag Central Office.
   *     Sven Hammarling, Nag Central Office.
   *     Richard Hanson, Sandia National Labs.
   *
   *
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  extern "C" {
    void
    BLASFUNC(ssymv)(
      character const   UPLO[],
      integer   const * N,
      real      const * ALPHA,
      real      const   A[],
      integer   const * LDA,
      real      const   X[],
      integer   const * INCX,
      real      const * BETA,
      real              Y[],
      integer   const * INCY
    );

    void
    BLASFUNC(dsymv)(
      character  const   UPLO[],
      integer    const * N,
      doublereal const * ALPHA,
      doublereal const   A[],
      integer    const * LDA,
      doublereal const   X[],
      integer    const * INCX,
      doublereal const * BETA,
      doublereal         Y[],
      integer    const * INCY
    );
  }
  #endif

  inline
  void
  symv(
    ULselect const & UPLO,
    integer          N,
    real             ALPHA,
    real const       A[],
    integer          LDA,
    real const       X[],
    integer          INCX,
    real             BETA,
    real             Y[],
    integer          INCY
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { ssymv( uplo_blas[UPLO],
           &N, &ALPHA, A, &LDA, X, &INCX, &BETA, Y, &INCY ); }
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)
  { BLASFUNC(ssymv)(
      const_cast<character*>(uplo_blas[UPLO]),
      &N, &ALPHA,
      const_cast<real*>(A), &LDA,
      const_cast<real*>(X), &INCX,
      &BETA, Y, &INCY
    );
  }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { CBLASNAME(ssymv)( CblasColMajor, uplo_cblas[UPLO],
                      N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY ); }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  inline
  void
  symv(
    ULselect const & UPLO,
    integer          N,
    doublereal       ALPHA,
    doublereal const A[],
    integer          LDA,
    doublereal const X[],
    integer          INCX,
    doublereal       BETA,
    doublereal       Y[],
    integer          INCY
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { dsymv( uplo_blas[UPLO],
           &N, &ALPHA, A, &LDA, X, &INCX, &BETA, Y, &INCY ); }
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)
  { BLASFUNC(dsymv)(
      const_cast<character*>(uplo_blas[UPLO]),
      &N, &ALPHA,
      const_cast<doublereal*>(A), &LDA,
      const_cast<doublereal*>(X), &INCX,
      &BETA, Y, &INCY
    );
  }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { CBLASNAME(dsymv)( CblasColMajor, uplo_cblas[UPLO],
                      N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY ); }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  /*
  //   ___ _   _ _ __ ___  _ __ ___
  //  / __| | | | '_ ` _ \| '_ ` _ \
  //  \__ \ |_| | | | | | | | | | | |
  //  |___/\__, |_| |_| |_|_| |_| |_|
  //       |___/
  */
  /*\
   *
   *  Purpose
   *  =======
   *
   *  DSYMM  performs one of the matrix-matrix operations
   *
   *     C := alpha*A*B + beta*C,
   *
   *  or
   *
   *     C := alpha*B*A + beta*C,
   *
   *  where alpha and beta are scalars,  A is a symmetric matrix and  B and
   *  C are  m by n matrices.
   *
   *  Arguments
   *  ==========
   *
   *  SIDE   - CHARACTER*1.
   *           On entry,  SIDE  specifies whether  the  symmetric matrix  A
   *           appears on the  left or right  in the  operation as follows:
   *
   *              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
   *
   *              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
   *
   *           Unchanged on exit.
   *
   *  UPLO   - CHARACTER*1.
   *           On  entry,   UPLO  specifies  whether  the  upper  or  lower
   *           triangular  part  of  the  symmetric  matrix   A  is  to  be
   *           referenced as follows:
   *
   *              UPLO = 'U' or 'u'   Only the upper triangular part of the
   *                                  symmetric matrix is to be referenced.
   *
   *              UPLO = 'L' or 'l'   Only the lower triangular part of the
   *                                  symmetric matrix is to be referenced.
   *
   *           Unchanged on exit.
   *
   *  M      - INTEGER.
   *           On entry,  M  specifies the number of rows of the matrix  C.
   *           M  must be at least zero.
   *           Unchanged on exit.
   *
   *  N      - INTEGER.
   *           On entry, N specifies the number of columns of the matrix C.
   *           N  must be at least zero.
   *           Unchanged on exit.
   *
   *  ALPHA  - DOUBLE PRECISION.
   *           On entry, ALPHA specifies the scalar alpha.
   *           Unchanged on exit.
   *
   *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
   *           m  when  SIDE = 'L' or 'l'  and is  n otherwise.
   *           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
   *           the array  A  must contain the  symmetric matrix,  such that
   *           when  UPLO = 'U' or 'u', the leading m by m upper triangular
   *           part of the array  A  must contain the upper triangular part
   *           of the  symmetric matrix and the  strictly  lower triangular
   *           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
   *           the leading  m by m  lower triangular part  of the  array  A
   *           must  contain  the  lower triangular part  of the  symmetric
   *           matrix and the  strictly upper triangular part of  A  is not
   *           referenced.
   *           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
   *           the array  A  must contain the  symmetric matrix,  such that
   *           when  UPLO = 'U' or 'u', the leading n by n upper triangular
   *           part of the array  A  must contain the upper triangular part
   *           of the  symmetric matrix and the  strictly  lower triangular
   *           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
   *           the leading  n by n  lower triangular part  of the  array  A
   *           must  contain  the  lower triangular part  of the  symmetric
   *           matrix and the  strictly upper triangular part of  A  is not
   *           referenced.
   *           Unchanged on exit.
   *
   *  LDA    - INTEGER.
   *           On entry, LDA specifies the first dimension of A as declared
   *           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
   *           LDA must be at least  max( 1, m ), otherwise  LDA must be at
   *           least  max( 1, n ).
   *           Unchanged on exit.
   *
   *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
   *           Before entry, the leading  m by n part of the array  B  must
   *           contain the matrix B.
   *           Unchanged on exit.
   *
   *  LDB    - INTEGER.
   *           On entry, LDB specifies the first dimension of B as declared
   *           in  the  calling  (sub)  program.   LDB  must  be  at  least
   *           max( 1, m ).
   *           Unchanged on exit.
   *
   *  BETA   - DOUBLE PRECISION.
   *           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
   *           supplied as zero then C need not be set on input.
   *           Unchanged on exit.
   *
   *  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
   *           Before entry, the leading  m by n  part of the array  C must
   *           contain the matrix  C,  except when  beta  is zero, in which
   *           case C need not be set on entry.
   *           On exit, the array  C  is overwritten by the  m by n updated
   *           matrix.
   *
   *  LDC    - INTEGER.
   *           On entry, LDC specifies the first dimension of C as declared
   *           in  the  calling  (sub)  program.   LDC  must  be  at  least
   *           max( 1, m ).
   *           Unchanged on exit.
   *
   *
   *  Level 3 Blas routine.
   *
   *  -- Written on 8-February-1989.
   *     Jack Dongarra, Argonne National Laboratory.
   *     Iain Duff, AERE Harwell.
   *     Jeremy Du Croz, Numerical Algorithms Group Ltd.
   *     Sven Hammarling, Numerical Algorithms Group Ltd.
   *
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  extern "C" {
    void
    BLASFUNC(ssymm)(
      character const   SIDE[],
      character const   UPLO[],
      integer   const * M,
      integer   const * N,
      real      const * ALPHA,
      real      const   A[],
      integer   const * LDA,
      real      const   B[],
      integer   const * LDB,
      real      const * BETA,
      real              C[],
      integer   const * LDC
    );

    void
    BLASFUNC(dsymm)(
      character  const   SIDE[],
      character  const   UPLO[],
      integer    const * M,
      integer    const * N,
      doublereal const * ALPHA,
      doublereal const   A[],
      integer    const * LDA,
      doublereal const   B[],
      integer    const * LDB,
      doublereal const * BETA,
      doublereal         C[],
      integer    const * LDC
    );
  }
  #endif

  inline
  void
  symm(
    SideMultiply const & SIDE,
    ULselect     const & UPLO,
    integer              M,
    integer              N,
    real                 ALPHA,
    real           const A[],
    integer              LDA,
    real           const B[],
    integer              LDB,
    real                 BETA,
    real                 C[],
    integer              LDC
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { ssymm(
      side_blas[SIDE], uplo_blas[UPLO],
      &M, &N, &ALPHA, A, &LDA, B, &LDB,
      &BETA, C, &LDC
    );
  }
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)
  { BLASFUNC(ssymm)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(uplo_blas[UPLO]),
      &M, &N,
      &ALPHA, const_cast<real*>(A), &LDA,
      const_cast<real*>(B), &LDB,
      &BETA, C, &LDC
    );
  }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { CBLASNAME(ssymm)(
      CblasColMajor,
      side_cblas[SIDE],
      uplo_cblas[UPLO],
      M, N,
      ALPHA, A, LDA,
      B, LDB,
      BETA, C, LDC
    );
  }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  inline
  void
  symm(
    SideMultiply const & SIDE,
    ULselect     const & UPLO,
    integer              M,
    integer              N,
    doublereal           ALPHA,
    doublereal const     A[],
    integer              LDA,
    doublereal const     B[],
    integer              LDB,
    doublereal           BETA,
    doublereal           C[],
    integer              LDC
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { dsymm(
      side_blas[SIDE], uplo_blas[UPLO],
      &M, &N, &ALPHA, A, &LDA, B, &LDB,
      &BETA, C, &LDC
    );
  }
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)
  { BLASFUNC(dsymm)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(uplo_blas[UPLO]),
      &M, &N,
      &ALPHA, const_cast<doublereal*>(A), &LDA,
      const_cast<doublereal*>(B), &LDB,
      &BETA, C, &LDC
    );
  }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { CBLASNAME(dsymm)(
      CblasColMajor,
      side_cblas[SIDE],
      uplo_cblas[UPLO],
      M, N,
      ALPHA, A, LDA,
      B, LDB,
      BETA, C, LDC
    );
  }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

}

///
/// eof: symmetric.hxx
///
