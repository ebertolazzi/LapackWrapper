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
/// file: triangular.hxx
///

namespace lapack_wrapper {

  /*
  //   _
  //  | |_ _ __ _ __ _____   __
  //  | __| '__| '_ ` _ \ \ / /
  //  | |_| |  | | | | | \ V /
  //   \__|_|  |_| |_| |_|\_/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DTRMV  performs one of the matrix-vector operations
   *
   *     x := A*x,   or   x := A**T*x,
   *
   *  where x is an n element vector and  A is an n by n unit, or non-unit,
   *  upper or lower triangular matrix.
   *
   *  Arguments
   *  ==========
   *
   *  UPLO   - CHARACTER*1.
   *           On entry, UPLO specifies whether the matrix is an upper or
   *           lower triangular matrix as follows:
   *
   *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
   *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
   *
   *           Unchanged on exit.
   *
   *  TRANS  - CHARACTER*1.
   *           On entry, TRANS specifies the operation to be performed as
   *           follows:
   *
   *              TRANS = 'N' or 'n'   x := A*x.
   *              TRANS = 'T' or 't'   x := A**T*x.
   *              TRANS = 'C' or 'c'   x := A**T*x.
   *
   *           Unchanged on exit.
   *
   *  DIAG   - CHARACTER*1.
   *           On entry, DIAG specifies whether or not A is unit
   *           triangular as follows:
   *
   *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
   *              DIAG = 'N' or 'n'   A is not assumed to be unit triangular.
   *
   *           Unchanged on exit.
   *
   *  N      - INTEGER.
   *           On entry, N specifies the order of the matrix A.
   *           N must be at least zero.
   *           Unchanged on exit.
   *
   *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
   *           Before entry with  UPLO = 'U' or 'u', the leading n by n
   *           upper triangular part of the array A must contain the upper
   *           triangular matrix and the strictly lower triangular part of
   *           A is not referenced.
   *           Before entry with UPLO = 'L' or 'l', the leading n by n
   *           lower triangular part of the array A must contain the lower
   *           triangular matrix and the strictly upper triangular part of
   *           A is not referenced.
   *           Note that when  DIAG = 'U' or 'u', the diagonal elements of
   *           A are not referenced either, but are assumed to be unity.
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
   *           element vector x. On exit, X is overwritten with the
   *           tranformed vector x.
   *
   *  INCX   - INTEGER.
   *           On entry, INCX specifies the increment for the elements of
   *           X. INCX must not be zero.
   *           Unchanged on exit.
   *
   *  Further Details
   *  ===============
   *
   *  Level 2 Blas routine.
   *  The vector and matrix arguments are not referenced when N = 0, or M = 0
  \*/
  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  extern "C" {

    void
    BLASFUNC(strmv)(
      character const   UPLO[],
      character const   TRANS[],
      character const   DIAG[],
      integer   const * N,
      real      const   A[],
      integer   const * LDA,
      real              X[],
      integer   const * INCX
    );

    void
    BLASFUNC(dtrmv)(
      character  const   UPLO[],
      character  const   TRANS[],
      character  const   DIAG[],
      integer    const * N,
      doublereal const   A[],
      integer    const * LDA,
      doublereal         X[],
      integer    const * INCX
    );
  }
  #endif

  inline
  void
  trmv(
    ULselect      const & UPLO,
    Transposition const & TRANS,
    DiagonalType  const & DIAG,
    integer               N,
    real          const   A[],
    integer               LDA,
    real                  X[],
    integer               INCX
  ) {
  #if defined(LAPACK_WRAPPER_USE_MKL)
    strmv(
      uplo_blas[UPLO], trans_blas[TRANS], diag_blas[DIAG],
      &N, A, &LDA, X, &INCX
    );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    CBLASNAME(strmv)(
      CblasColMajor,
      uplo_cblas[UPLO],
      trans_cblas[TRANS],
      diag_cblas[DIAG],
      N, A, LDA, X, INCX
    );
  #else
    BLASFUNC(strmv)(
      const_cast<character*>(uplo_blas[UPLO]),
      const_cast<character*>(trans_blas[TRANS]),
      const_cast<character*>(diag_blas[DIAG]),
      &N, const_cast<real*>(A), &LDA, X, &INCX
    );
  #endif
  }

  inline
  void
  trmv(
    ULselect      const & UPLO,
    Transposition const & TRANS,
    DiagonalType  const & DIAG,
    integer               N,
    doublereal    const   A[],
    integer               LDA,
    doublereal            X[],
    integer               INCX
  ) {
  #if defined(LAPACK_WRAPPER_USE_MKL)
    dtrmv(
      uplo_blas[UPLO], trans_blas[TRANS], diag_blas[DIAG],
      &N, A, &LDA, X, &INCX
    );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    CBLASNAME(dtrmv)(
      CblasColMajor,
      uplo_cblas[UPLO],
      trans_cblas[TRANS],
      diag_cblas[DIAG],
      N, A, LDA, X, INCX
    );
  #else
    BLASFUNC(dtrmv)(
      const_cast<character*>(uplo_blas[UPLO]),
      const_cast<character*>(trans_blas[TRANS]),
      const_cast<character*>(diag_blas[DIAG]),
      &N, const_cast<doublereal*>(A), &LDA, X, &INCX
    );
  #endif
  }

  /*
  //   _
  //  | |_ _ __ _____   __
  //  | __| '__/ __\ \ / /
  //  | |_| |  \__ \\ V /
  //   \__|_|  |___/ \_/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DTRSV  solves one of the systems of equations
   *
   *     A*x = b,   or   A**T*x = b,
   *
   *  where b and x are n element vectors and A is an n by n unit, or
   *  non-unit, upper or lower triangular matrix.
   *
   *  No test for singularity or near-singularity is included in this
   *  routine. Such tests must be performed before calling this routine.
   *
   *  Arguments
   *  ==========
   *
   *  UPLO   - CHARACTER*1.
   *           On entry, UPLO specifies whether the matrix is an upper or
   *           lower triangular matrix as follows:
   *
   *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
   *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
   *
   *           Unchanged on exit.
   *
   *  TRANS  - CHARACTER*1.
   *           On entry, TRANS specifies the equations to be solved as
   *           follows:
   *
   *              TRANS = 'N' or 'n'   A*x = b.
   *              TRANS = 'T' or 't'   A**T*x = b.
   *              TRANS = 'C' or 'c'   A**T*x = b.
   *
   *           Unchanged on exit.
   *
   *  DIAG   - CHARACTER*1.
   *           On entry, DIAG specifies whether or not A is unit
   *           triangular as follows:
   *
   *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
   *              DIAG = 'N' or 'n'   A is not assumed to be unit triangular.
   *
   *           Unchanged on exit.
   *
   *  N      - INTEGER.
   *           On entry, N specifies the order of the matrix A.
   *           N must be at least zero.
   *           Unchanged on exit.
   *
   *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
   *           Before entry with  UPLO = 'U' or 'u', the leading n by n
   *           upper triangular part of the array A must contain the upper
   *           triangular matrix and the strictly lower triangular part of
   *           A is not referenced.
   *           Before entry with UPLO = 'L' or 'l', the leading n by n
   *           lower triangular part of the array A must contain the lower
   *           triangular matrix and the strictly upper triangular part of
   *           A is not referenced.
   *           Note that when  DIAG = 'U' or 'u', the diagonal elements of
   *           A are not referenced either, but are assumed to be unity.
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
   *           element right-hand side vector b. On exit, X is overwritten
   *           with the solution vector x.
   *
   *  INCX   - INTEGER.
   *           On entry, INCX specifies the increment for the elements of
   *           X. INCX must not be zero.
   *           Unchanged on exit.
   *
   *
   *  Level 2 Blas routine.
  \*/
  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  extern "C" {

    void
    BLASFUNC(strsv)(
      character const   UPLO[],
      character const   TRANS[],
      character const   DIAG[],
      integer   const * N,
      real      const   A[],
      integer   const * LDA,
      real              X[],
      integer   const * INCX
    );

    void
    BLASFUNC(dtrsv)(
      character  const   UPLO[],
      character  const   TRANS[],
      character  const   DIAG[],
      integer    const * N,
      doublereal const   A[],
      integer    const * LDA,
      doublereal         X[],
      integer    const * INCX
    );
  }
  #endif

  inline
  void
  trsv(
    ULselect      const & UPLO,
    Transposition const & TRANS,
    DiagonalType  const & DIAG,
    integer               N,
    real          const   A[],
    integer               LDA,
    real                  X[],
    integer               INCX
  ) {
  #if defined(LAPACK_WRAPPER_USE_MKL)
    strsv(
      uplo_blas[UPLO], trans_blas[TRANS], diag_blas[DIAG],
      &N, A, &LDA, X, &INCX
    );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    CBLASNAME(strsv)(
      CblasColMajor,
      uplo_cblas[UPLO],
      trans_cblas[TRANS],
      diag_cblas[DIAG],
      N, A, LDA, X, INCX
    );
  #else
    BLASFUNC(strsv)(
      const_cast<character*>(uplo_blas[UPLO]),
      const_cast<character*>(trans_blas[TRANS]),
      const_cast<character*>(diag_blas[DIAG]),
      &N, const_cast<real*>(A), &LDA, X, &INCX
    );
  #endif
  }

  inline
  void
  trsv(
    ULselect      const & UPLO,
    Transposition const & TRANS,
    DiagonalType  const & DIAG,
    integer               N,
    doublereal    const   A[],
    integer               LDA,
    doublereal            X[],
    integer               INCX
  ) {
  #if defined(LAPACK_WRAPPER_USE_MKL)
    dtrsv(
      uplo_blas[UPLO], trans_blas[TRANS], diag_blas[DIAG],
      &N, A, &LDA, X, &INCX
    );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    CBLASNAME(dtrsv)(
      CblasColMajor,
      uplo_cblas[UPLO],
      trans_cblas[TRANS],
      diag_cblas[DIAG],
      N, A, LDA, X, INCX
    );
  #else
    BLASFUNC(dtrsv)(
      const_cast<character*>(uplo_blas[UPLO]),
      const_cast<character*>(trans_blas[TRANS]),
      const_cast<character*>(diag_blas[DIAG]),
      &N, const_cast<doublereal*>(A), &LDA, X, &INCX
    );
  #endif
  }

  /*
  //   _
  //  | |_ _ __ _ __ ___  _ __ ___
  //  | __| '__| '_ ` _ \| '_ ` _ \
  //  | |_| |  | | | | | | | | | | |
  //   \__|_|  |_| |_| |_|_| |_| |_|
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DTRMM  performs one of the matrix-matrix operations
   *
   *     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
   *
   *  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
   *  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
   *
   *     op( A ) = A   or   op( A ) = A**T.
   *
   *  Arguments
   *  ==========
   *
   *  SIDE   - CHARACTER*1.
   *           On entry,  SIDE specifies whether  op( A ) multiplies B from
   *           the left or right as follows:
   *
   *              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
   *              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
   *
   *           Unchanged on exit.
   *
   *  UPLO   - CHARACTER*1.
   *           On entry, UPLO specifies whether the matrix A is an upper or
   *           lower triangular matrix as follows:
   *
   *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
   *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
   *
   *           Unchanged on exit.
   *
   *  TRANSA - CHARACTER*1.
   *           On entry, TRANSA specifies the form of op( A ) to be used in
   *           the matrix multiplication as follows:
   *
   *              TRANSA = 'N' or 'n'   op( A ) = A.
   *              TRANSA = 'T' or 't'   op( A ) = A**T.
   *              TRANSA = 'C' or 'c'   op( A ) = A**T.
   *
   *           Unchanged on exit.
   *
   *  DIAG   - CHARACTER*1.
   *           On entry, DIAG specifies whether or not A is unit triangular
   *           as follows:
   *
   *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
   *              DIAG = 'N' or 'n'   A is not assumed to be unit triangular.
   *
   *           Unchanged on exit.
   *
   *  M      - INTEGER.
   *           On entry, M specifies the number of rows of B. M must be at
   *           least zero.
   *           Unchanged on exit.
   *
   *  N      - INTEGER.
   *           On entry, N specifies the number of columns of B.  N must be
   *           at least zero.
   *           Unchanged on exit.
   *
   *  ALPHA  - DOUBLE PRECISION.
   *           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
   *           zero then  A is not referenced and  B need not be set before
   *           entry.
   *           Unchanged on exit.
   *
   *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
   *           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
   *           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
   *           upper triangular part of the array  A must contain the upper
   *           triangular matrix  and the strictly lower triangular part of
   *           A is not referenced.
   *           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
   *           lower triangular part of the array  A must contain the lower
   *           triangular matrix  and the strictly upper triangular part of
   *           A is not referenced.
   *           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
   *           A  are not referenced either,  but are assumed to be  unity.
   *           Unchanged on exit.
   *
   *  LDA    - INTEGER.
   *           On entry, LDA specifies the first dimension of A as declared
   *           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
   *           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
   *           then LDA must be at least max( 1, n ).
   *           Unchanged on exit.
   *
   *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
   *           Before entry,  the leading  m by n part of the array  B must
   *           contain the matrix  B,  and  on exit  is overwritten  by the
   *           transformed matrix.
   *
   *  LDB    - INTEGER.
   *           On entry, LDB specifies the first dimension of B as declared
   *           in  the  calling  (sub)  program.   LDB  must  be  at  least
   *           max( 1, m ).
   *           Unchanged on exit.
   *
   *  Further Details
   *  ===============
   *
   *  Level 3 Blas routine.
  \*/
  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  extern "C" {

    void
    BLASFUNC(strmm) (
      character const   SIDE[],     // "L" or "R"
      character const   UPLO[],     // "U" or "L"
      character const   TRANSA[],   // "N", "T", "C"
      character const   DIAG[],     // "N", "U"
      integer   const * M,
      integer   const * N,
      real      const * ALPHA,
      real      const   A[],
      integer   const * LDA,
      real              B[],
      integer   const * LDB
    );

    void
    BLASFUNC(dtrmm) (
      character  const   SIDE[],     // "L" or "R"
      character  const   UPLO[],     // "U" or "L"
      character  const   TRANSA[],   // "N", "T", "C"
      character  const   DIAG[],     // "N", "U"
      integer    const * M,
      integer    const * N,
      doublereal const * ALPHA,
      doublereal const   A[],
      integer    const * LDA,
      doublereal         B[],
      integer    const * LDB
    );
  }
  #endif

  inline
  void
  trmm(
    SideMultiply  const & SIDE,
    ULselect      const & UPLO,
    Transposition const & TRANSA,
    DiagonalType  const & DIAG,
    integer               M,
    integer               N,
    real                  alpha,
    real          const   A[],
    integer               LDA,
    real                  B[],
    integer               LDB
  ) {
  #if defined(LAPACK_WRAPPER_USE_MKL)
    strmm(
      side_blas[SIDE], uplo_blas[UPLO],
      trans_blas[TRANSA], diag_blas[DIAG],
      &M, &N, &alpha, A, &LDA, B, &LDB
    );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    CBLASNAME(strmm) (
      CblasColMajor,
      side_cblas[SIDE],
      uplo_cblas[UPLO],
      trans_cblas[TRANSA],
      diag_cblas[DIAG],
      M, N, alpha, A, LDA, B, LDB
    );
  #else
    BLASFUNC(strmm)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(uplo_blas[UPLO]),
      const_cast<character*>(trans_blas[TRANSA]),
      const_cast<character*>(diag_blas[DIAG]),
      &M, &N, &alpha,
      const_cast<real*>(A), &LDA, B, &LDB
    );
  #endif
  }

  inline
  void
  trmm(
    SideMultiply  const & SIDE,
    ULselect      const & UPLO,
    Transposition const & TRANSA,
    DiagonalType  const & DIAG,
    integer               M,
    integer               N,
    doublereal            alpha,
    doublereal    const   A[],
    integer               LDA,
    doublereal            B[],
    integer               LDB
  ) {
  #if defined(LAPACK_WRAPPER_USE_MKL)
    dtrmm(
      side_blas[SIDE], uplo_blas[UPLO],
      trans_blas[TRANSA], diag_blas[DIAG],
      &M, &N, &alpha, A, &LDA, B, &LDB
    );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    CBLASNAME(dtrmm) (
      CblasColMajor,
      side_cblas[SIDE],
      uplo_cblas[UPLO],
      trans_cblas[TRANSA],
      diag_cblas[DIAG],
      M, N, alpha, A, LDA, B, LDB
    );
  #else
    BLASFUNC(dtrmm)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(uplo_blas[UPLO]),
      const_cast<character*>(trans_blas[TRANSA]),
      const_cast<character*>(diag_blas[DIAG]),
      &M, &N, &alpha,
      const_cast<doublereal*>(A), &LDA, B, &LDB
    );
  #endif
  }

  /*
  //   _
  //  | |_ _ __ ___ _ __ ___
  //  | __| '__/ __| '_ ` _ \
  //  | |_| |  \__ \ | | | | |
  //   \__|_|  |___/_| |_| |_|
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DTRSM  solves one of the matrix equations
   *
   *     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
   *
   *  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
   *  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
   *
   *     op( A ) = A   or   op( A ) = A**T.
   *
   *  The matrix X is overwritten on B.
   *
   *  Arguments
   *  ==========
   *
   *  SIDE   - CHARACTER*1.
   *           On entry, SIDE specifies whether op( A ) appears on the left
   *           or right of X as follows:
   *
   *              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
   *              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
   *
   *           Unchanged on exit.
   *
   *  UPLO   - CHARACTER*1.
   *           On entry, UPLO specifies whether the matrix A is an upper or
   *           lower triangular matrix as follows:
   *
   *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
   *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
   *
   *           Unchanged on exit.
   *
   *  TRANSA - CHARACTER*1.
   *           On entry, TRANSA specifies the form of op( A ) to be used in
   *           the matrix multiplication as follows:
   *
   *              TRANSA = 'N' or 'n'   op( A ) = A.
   *              TRANSA = 'T' or 't'   op( A ) = A**T.
   *              TRANSA = 'C' or 'c'   op( A ) = A**T.
   *
   *           Unchanged on exit.
   *
   *  DIAG   - CHARACTER*1.
   *           On entry, DIAG specifies whether or not A is unit triangular
   *           as follows:
   *
   *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
   *              DIAG = 'N' or 'n'   A is not assumed to be unit triangular.
   *
   *           Unchanged on exit.
   *
   *  M      - INTEGER.
   *           On entry, M specifies the number of rows of B. M must be at
   *           least zero.
   *           Unchanged on exit.
   *
   *  N      - INTEGER.
   *           On entry, N specifies the number of columns of B.  N must be
   *           at least zero.
   *           Unchanged on exit.
   *
   *  ALPHA  - DOUBLE PRECISION.
   *           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
   *           zero then  A is not referenced and  B need not be set before
   *           entry.
   *           Unchanged on exit.
   *
   *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
   *           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
   *           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
   *           upper triangular part of the array  A must contain the upper
   *           triangular matrix  and the strictly lower triangular part of
   *           A is not referenced.
   *           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
   *           lower triangular part of the array  A must contain the lower
   *           triangular matrix  and the strictly upper triangular part of
   *           A is not referenced.
   *           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
   *           A  are not referenced either,  but are assumed to be  unity.
   *           Unchanged on exit.
   *
   *  LDA    - INTEGER.
   *           On entry, LDA specifies the first dimension of A as declared
   *           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
   *           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
   *           then LDA must be at least max( 1, n ).
   *           Unchanged on exit.
   *
   *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
   *           Before entry,  the leading  m by n part of the array  B must
   *           contain  the  right-hand  side  matrix  B,  and  on exit  is
   *           overwritten by the solution matrix  X.
   *
   *  LDB    - INTEGER.
   *           On entry, LDB specifies the first dimension of B as declared
   *           in  the  calling  (sub)  program.   LDB  must  be  at  least
   *           max( 1, m ).
   *           Unchanged on exit.
   *
   *  Further Details
   *  ===============
   *
   *  Level 3 Blas routine.
  \*/
  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  extern "C" {

    void
    BLASFUNC(strsm)(
      character const   SIDE[],
      character const   UPLO[],
      character const   TRANSA[],
      character const   DIAG[],
      integer   const * M,
      integer   const * N,
      real      const * alpha,
      real      const   A[],
      integer   const * LDA,
      real              B[],
      integer   const * LDB
    );

    void
    BLASFUNC(dtrsm)(
      character  const   SIDE[],
      character  const   UPLO[],
      character  const   TRANSA[],
      character  const   DIAG[],
      integer    const * M,
      integer    const * N,
      doublereal const * alpha,
      doublereal const   A[],
      integer    const * LDA,
      doublereal         B[],
      integer    const * LDB
    );
  }
  #endif

  inline
  void
  trsm(
    SideMultiply  const & SIDE,
    ULselect      const & UPLO,
    Transposition const & TRANS,
    DiagonalType  const & DIAG,
    integer               M,
    integer               N,
    real                  alpha,
    real          const   A[],
    integer               LDA,
    real                  B[],
    integer               LDB
  ) {
  #if defined(LAPACK_WRAPPER_USE_MKL)
    strsm(
      side_blas[SIDE], uplo_blas[UPLO],
      trans_blas[TRANS], diag_blas[DIAG],
      &M, &N, &alpha, A, &LDA, B, &LDB
    );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    CBLASNAME(strsm)(
      CblasColMajor,
      side_cblas[SIDE],
      uplo_cblas[UPLO],
      trans_cblas[TRANS],
      diag_cblas[DIAG],
      M, N, alpha, A, LDA, B, LDB
    );
  #else
    BLASFUNC(strsm)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(uplo_blas[UPLO]),
      const_cast<character*>(trans_blas[TRANS]),
      const_cast<character*>(diag_blas[DIAG]),
      &M, &N, &alpha,
      const_cast<real*>(A), &LDA, B, &LDB
    );
  #endif
  }

  inline
  void
  trsm(
    SideMultiply  const & SIDE,
    ULselect      const & UPLO,
    Transposition const & TRANS,
    DiagonalType  const & DIAG,
    integer               M,
    integer               N,
    doublereal            alpha,
    doublereal    const   A[],
    integer               LDA,
    doublereal            B[],
    integer               LDB
  ) {
  #if defined(LAPACK_WRAPPER_USE_MKL)
    dtrsm(
      side_blas[SIDE], uplo_blas[UPLO],
      trans_blas[TRANS], diag_blas[DIAG],
      &M, &N, &alpha, A, &LDA, B, &LDB
    );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    CBLASNAME(dtrsm)(
      CblasColMajor,
      side_cblas[SIDE],
      uplo_cblas[UPLO],
      trans_cblas[TRANS],
      diag_cblas[DIAG],
      M, N, alpha, A, LDA, B, LDB
    );
  #else
    BLASFUNC(dtrsm)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(uplo_blas[UPLO]),
      const_cast<character*>(trans_blas[TRANS]),
      const_cast<character*>(diag_blas[DIAG]),
      &M, &N, &alpha,
      const_cast<doublereal*>(A), &LDA, B, &LDB
    );
  #endif
  }
}

///
/// eof: triangular.hxx
///
