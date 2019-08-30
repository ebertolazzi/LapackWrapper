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
/// file: banded.hxx
///

namespace lapack_wrapper {

  /*
  //   ____                  _   __  __       _        _
  //  | __ )  __ _ _ __   __| | |  \/  | __ _| |_ _ __(_)_  __
  //  |  _ \ / _` | '_ \ / _` | | |\/| |/ _` | __| '__| \ \/ /
  //  | |_) | (_| | | | | (_| | | |  | | (_| | |_| |  | |>  <
  //  |____/ \__,_|_| |_|\__,_| |_|  |_|\__,_|\__|_|  |_/_/\_\
  */

  /*\
   *  Purpose
   *  =======
   *
   *  DTBMV  performs one of the matrix-vector operations
   *
   *     x := A*x,   or   x := A'*x,
   *
   *  where x is an n element vector and  A is an n by n unit, or non-unit,
   *  upper or lower triangular band matrix, with ( k + 1 ) diagonals.
   *
   *  Parameters
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
   *              TRANS = 'T' or 't'   x := A'*x.
   *              TRANS = 'C' or 'c'   x := A'*x.
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
   *  K      - INTEGER.
   *           On entry with UPLO = 'U' or 'u', K specifies the number of
   *           super-diagonals of the matrix A.
   *           On entry with UPLO = 'L' or 'l', K specifies the number of
   *           sub-diagonals of the matrix A.
   *           K must satisfy  0 .le. K.
   *           Unchanged on exit.
   *
   *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
   *           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
   *           by n part of the array A must contain the upper triangular
   *           band part of the matrix of coefficients, supplied column by
   *           column, with the leading diagonal of the matrix in row
   *           ( k + 1 ) of the array, the first super-diagonal starting at
   *           position 2 in row k, and so on. The top left k by k triangle
   *           of the array A is not referenced.
   *           The following program segment will transfer an upper
   *           triangular band matrix from conventional full matrix storage
   *           to band storage:
   *
   *                 DO 20, J = 1, N
   *                    M = K + 1 - J
   *                    DO 10, I = MAX( 1, J - K ), J
   *                       A( M + I, J ) = matrix( I, J )
   *              10    CONTINUE
   *              20 CONTINUE
   *
   *           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
   *           by n part of the array A must contain the lower triangular
   *           band part of the matrix of coefficients, supplied column by
   *           column, with the leading diagonal of the matrix in row 1 of
   *           the array, the first sub-diagonal starting at position 1 in
   *           row 2, and so on. The bottom right k by k triangle of the
   *           array A is not referenced.
   *           The following program segment will transfer a lower
   *           triangular band matrix from conventional full matrix storage
   *           to band storage:
   *
   *                 DO 20, J = 1, N
   *                    M = 1 - J
   *                    DO 10, I = J, MIN( N, J + K )
   *                       A( M + I, J ) = matrix( I, J )
   *              10    CONTINUE
   *              20 CONTINUE
   *
   *           Note that when DIAG = 'U' or 'u' the elements of the array A
   *           corresponding to the diagonal elements of the matrix are not
   *           referenced, but are assumed to be unity.
   *           Unchanged on exit.
   *
   *  LDA    - INTEGER.
   *           On entry, LDA specifies the first dimension of A as declared
   *           in the calling (sub) program. LDA must be at least
   *           ( k + 1 ).
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
    BLASFUNC(stbmv)(
      char    const * UPLO,
      char    const * TRANS,
      char    const * DIAG,
      integer const * N,
      integer const * K,
      real    const   A[],
      integer const * LDA,
      real            X[],
      integer const * INCXU
    );

    void
    BLASFUNC(dtbmv)(
      char       const * UPLO,
      char       const * TRANS,
      char       const * DIAG,
      integer    const * N,
      integer    const * K,
      doublereal const   A[],
      integer    const * LDA,
      doublereal         X[],
      integer    const * INCXU
    );
  }
  #endif

  inline
  void
  tbmv(
    ULselect      const & UPLO,
    Transposition const & TRANS,
    DiagonalType  const & DIAG,
    integer             N,
    integer             K,
    real          const A[],
    integer             ldA,
    real                xb[],
    integer             incx
  ) {
    #if defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    CBLASNAME(stbmv)(
      CblasColMajor,
      uplo_cblas[UPLO],
      trans_cblas[TRANS],
      diag_cblas[DIAG],
      N, K, A, ldA, xb, incx
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    stbmv(
      uplo_blas[UPLO], trans_blas[TRANS], diag_blas[DIAG],
      &N, &K, A, &ldA, xb, &incx
    );
    #else
    BLASFUNC(stbmv)(
      const_cast<character*>(uplo_blas[UPLO]),
      const_cast<character*>(trans_blas[TRANS]),
      const_cast<character*>(diag_blas[DIAG]),
      &N, &K,
      const_cast<real*>(A), &ldA, xb, &incx
    );
    #endif
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline
  void
  tbmv(
    ULselect      const & UPLO,
    Transposition const & TRANS,
    DiagonalType  const & DIAG,
    integer             N,
    integer             K,
    doublereal    const A[],
    integer             ldA,
    doublereal          xb[],
    integer             incx
  ) {
    #if defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    CBLASNAME(dtbmv)(
      CblasColMajor,
      uplo_cblas[UPLO],
      trans_cblas[TRANS],
      diag_cblas[DIAG],
      N, K, A, ldA, xb, incx
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dtbmv(
      uplo_blas[UPLO], trans_blas[TRANS], diag_blas[DIAG],
      &N, &K, A, &ldA, xb, &incx
    );
    #else
    BLASFUNC(dtbmv)(
      const_cast<character*>(uplo_blas[UPLO]),
      const_cast<character*>(trans_blas[TRANS]),
      const_cast<character*>(diag_blas[DIAG]),
      &N, &K,
      const_cast<doublereal*>(A), &ldA, xb, &incx
    );
    #endif
  }

  /*\
   *  Purpose
   *  =======
   *
   *  DTBSV  solves one of the systems of equations
   *
   *     A*x = b,   or   A'*x = b,
   *
   *  where b and x are n element vectors and A is an n by n unit, or
   *  non-unit, upper or lower triangular band matrix, with ( k + 1 )
   *  diagonals.
   *
   *  No test for singularity or near-singularity is included in this
   *  routine. Such tests must be performed before calling this routine.
   *
   *  Parameters
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
   *              TRANS = 'T' or 't'   A'*x = b.
   *              TRANS = 'C' or 'c'   A'*x = b.
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
   *  K      - INTEGER.
   *           On entry with UPLO = 'U' or 'u', K specifies the number of
   *           super-diagonals of the matrix A.
   *           On entry with UPLO = 'L' or 'l', K specifies the number of
   *           sub-diagonals of the matrix A.
   *           K must satisfy  0 .le. K.
   *           Unchanged on exit.
   *
   *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
   *           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
   *           by n part of the array A must contain the upper triangular
   *           band part of the matrix of coefficients, supplied column by
   *           column, with the leading diagonal of the matrix in row
   *           ( k + 1 ) of the array, the first super-diagonal starting at
   *           position 2 in row k, and so on. The top left k by k triangle
   *           of the array A is not referenced.
   *           The following program segment will transfer an upper
   *           triangular band matrix from conventional full matrix storage
   *           to band storage:
   *
   *                 DO 20, J = 1, N
   *                    M = K + 1 - J
   *                    DO 10, I = MAX( 1, J - K ), J
   *                       A( M + I, J ) = matrix( I, J )
   *              10    CONTINUE
   *              20 CONTINUE
   *
   *           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
   *           by n part of the array A must contain the lower triangular
   *           band part of the matrix of coefficients, supplied column by
   *           column, with the leading diagonal of the matrix in row 1 of
   *           the array, the first sub-diagonal starting at position 1 in
   *           row 2, and so on. The bottom right k by k triangle of the
   *           array A is not referenced.
   *           The following program segment will transfer a lower
   *           triangular band matrix from conventional full matrix storage
   *           to band storage:
   *
   *                 DO 20, J = 1, N
   *                    M = 1 - J
   *                    DO 10, I = J, MIN( N, J + K )
   *                       A( M + I, J ) = matrix( I, J )
   *              10    CONTINUE
   *              20 CONTINUE
   *
   *           Note that when DIAG = 'U' or 'u' the elements of the array A
   *           corresponding to the diagonal elements of the matrix are not
   *           referenced, but are assumed to be unity.
   *           Unchanged on exit.
   *
   *  LDA    - INTEGER.
   *           On entry, LDA specifies the first dimension of A as declared
   *           in the calling (sub) program. LDA must be at least
   *           ( k + 1 ).
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
   *
   *  -- Written on 22-October-1986.
   *     Jack Dongarra, Argonne National Lab.
   *     Jeremy Du Croz, Nag Central Office.
   *     Sven Hammarling, Nag Central Office.
   *     Richard Hanson, Sandia National Labs.
   *
  \*/
  #ifdef LAPACK_WRAPPER_USE_LAPACK
  extern "C" {

    void
    BLASFUNC(stbsv)(
      char    const * UPLO,
      char    const * TRANS,
      char    const * DIAG,
      integer const * N,
      integer const * K,
      real    const   A[],
      integer const * LDA,
      real            X[],
      integer const * INCXU
    );

    void
    BLASFUNC(dtbsv)(
      char       const * UPLO,
      char       const * TRANS,
      char       const * DIAG,
      integer    const * N,
      integer    const * K,
      doublereal const   A[],
      integer    const * LDA,
      doublereal         X[],
      integer    const * INCXU
    );
  }
  #endif

  inline
  void
  tbsv(
    ULselect      const & UPLO,
    Transposition const & TRANS,
    DiagonalType  const & DIAG,
    integer             N,
    integer             K,
    real          const A[],
    integer             ldA,
    real                xb[],
    integer             incx
  ) {
    #if defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    CBLASNAME(stbsv)(
      CblasColMajor,
      uplo_cblas[UPLO],
      trans_cblas[TRANS],
      diag_cblas[DIAG],
      N, K, A, ldA, xb, incx
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    stbsv(
      uplo_blas[UPLO], trans_blas[TRANS], diag_blas[DIAG],
      &N, &K, A, &ldA, xb, &incx
    );
    #else
    BLASFUNC(stbsv)(
      const_cast<character*>(uplo_blas[UPLO]),
      const_cast<character*>(trans_blas[TRANS]),
      const_cast<character*>(diag_blas[DIAG]),
      &N, &K,
      const_cast<real*>(A), &ldA, xb, &incx
    );
    #endif
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline
  void
  tbsv(
    ULselect      const & UPLO,
    Transposition const & TRANS,
    DiagonalType  const & DIAG,
    integer             N,
    integer             K,
    doublereal    const A[],
    integer             ldA,
    doublereal          xb[],
    integer             incx
  ) {
    #if defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    CBLASNAME(dtbsv)(
      CblasColMajor,
      uplo_cblas[UPLO],
      trans_cblas[TRANS],
      diag_cblas[DIAG],
      N, K, A, ldA, xb, incx
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dtbsv(
      uplo_blas[UPLO], trans_blas[TRANS], diag_blas[DIAG],
      &N, &K, A, &ldA, xb, &incx
    );
    #else
    BLASFUNC(dtbsv)(
      const_cast<character*>(uplo_blas[UPLO]),
      const_cast<character*>(trans_blas[TRANS]),
      const_cast<character*>(diag_blas[DIAG]),
      &N, &K,
      const_cast<doublereal*>(A), &ldA, xb, &incx
    );
    #endif
  }

  /*\
   |   ____                  _          _
   |  | __ )  __ _ _ __   __| | ___  __| |
   |  |  _ \ / _` | '_ \ / _` |/ _ \/ _` |
   |  | |_) | (_| | | | | (_| |  __/ (_| |
   |  |____/ \__,_|_| |_|\__,_|\___|\__,_|
  \*/

  /*\
   *  Purpose
   *  =======
   *
   *  DGBTRF computes an LU factorization of a real m-by-n band matrix A
   *  using partial pivoting with row interchanges.
   *
   *  This is the blocked version of the algorithm, calling Level 3 BLAS.
   *
   *  Arguments
   *  =========
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix A.  M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A.  N >= 0.
   *
   *  KL      (input) INTEGER
   *          The number of subdiagonals within the band of A.  KL >= 0.
   *
   *  KU      (input) INTEGER
   *          The number of superdiagonals within the band of A.  KU >= 0.
   *
   *  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
   *          On entry, the matrix A in band storage, in rows KL+1 to
   *          2*KL+KU+1; rows 1 to KL of the array need not be set.
   *          The j-th column of A is stored in the j-th column of the
   *          array AB as follows:
   *          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
   *
   *          On exit, details of the factorization: U is stored as an
   *          upper triangular band matrix with KL+KU superdiagonals in
   *          rows 1 to KL+KU+1, and the multipliers used during the
   *          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
   *          See below for further details.
   *
   *  LDAB    (input) INTEGER
   *          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
   *
   *  IPIV    (output) INTEGER array, dimension (min(M,N))
   *          The pivot indices; for 1 <= i <= min(M,N), row i of the
   *          matrix was interchanged with row IPIV(i).
   *
   *  INFO    (output) INTEGER
   *          = 0: successful exit
   *          < 0: if INFO = -i, the i-th argument had an illegal value
   *          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
   *               has been completed, but the factor U is exactly
   *               singular, and division by zero will occur if it is used
   *               to solve a system of equations.
   *
   *  Further Details
   *  ===============
   *
   *  The band storage scheme is illustrated by the following example, when
   *  M = N = 6, KL = 2, KU = 1:
   *
   *  On entry:                       On exit:
   *
   *      *    *    *    +    +    +       *    *    *   u14  u25  u36
   *      *    *    +    +    +    +       *    *   u13  u24  u35  u46
   *      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
   *     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
   *     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
   *     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
   *
   *  Array elements marked * are not used by the routine; elements marked
   *  + need not be set on entry, but are required by the routine to store
   *  elements of U because of fill-in resulting from the row interchanges.
   *
   *  =====================================================================
   *
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {

    void
    LAPACK_F77NAME(sgbtrf)(
      integer const * M,
      integer const * N,
      integer const * KL,
      integer const * KU,
      real            AB[],
      integer const * ldAB,
      integer         ipiv[],
      integer       * info
    );

    void
    LAPACK_F77NAME(dgbtrf)(
      integer const * M,
      integer const * N,
      integer const * KL,
      integer const * KU,
      doublereal      AB[],
      integer const * ldAB,
      integer         ipiv[],
      integer       * info
    );
  }
  #endif

  inline
  integer
  gbtrf(
    integer M,
    integer N,
    integer KL,
    integer KU,
    real    AB[],
    integer ldAB,
    integer ipiv[]
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sgbtrf)( &M, &N, &KL, &KU, AB, &ldAB, ipiv, &info );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgbtrf( &M, &N, &KL, &KU, AB, &ldAB, ipiv, &info );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgbtrf)( &M, &N, &KL, &KU, AB, &ldAB, ipiv, &info );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline
  integer
  gbtrf(
    integer    M,
    integer    N,
    integer    KL,
    integer    KU,
    doublereal AB[],
    integer    ldAB,
    integer    ipiv[]
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dgbtrf)( &M, &N, &KL, &KU, AB, &ldAB, ipiv, &info );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgbtrf( &M, &N, &KL, &KU, AB, &ldAB, ipiv, &info );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgbtrf)( &M, &N, &KL, &KU, AB, &ldAB, ipiv, &info );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  /*\
   *
   *  Purpose
   *  =======
   *
   *  DGBTRS solves a system of linear equations
   *     A * X = B  or  A**T * X = B
   *  with a general band matrix A using the LU factorization computed
   *  by DGBTRF.
   *
   *  Arguments
   *  =========
   *
   *  TRANS   (input) CHARACTER*1
   *          Specifies the form of the system of equations.
   *          = 'N':  A * X = B  (No transpose)
   *          = 'T':  A**T* X = B  (Transpose)
   *          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
   *
   *  N       (input) INTEGER
   *          The order of the matrix A.  N >= 0.
   *
   *  KL      (input) INTEGER
   *          The number of subdiagonals within the band of A.  KL >= 0.
   *
   *  KU      (input) INTEGER
   *          The number of superdiagonals within the band of A.  KU >= 0.
   *
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of columns
   *          of the matrix B.  NRHS >= 0.
   *
   *  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
   *          Details of the LU factorization of the band matrix A, as
   *          computed by DGBTRF.  U is stored as an upper triangular band
   *          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
   *          the multipliers used during the factorization are stored in
   *          rows KL+KU+2 to 2*KL+KU+1.
   *
   *  LDAB    (input) INTEGER
   *          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
   *
   *  IPIV    (input) INTEGER array, dimension (N)
   *          The pivot indices; for 1 <= i <= N, row i of the matrix was
   *          interchanged with row IPIV(i).
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
   *          On entry, the right hand side matrix B.
   *          On exit, the solution matrix X.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B.  LDB >= max(1,N).
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0: if INFO = -i, the i-th argument had an illegal value
   *
   *  =====================================================================
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {

    void
    LAPACK_F77NAME(sgbtrs)(
      character const * TRANS,
      integer   const * N,
      integer   const * KL,
      integer   const * KU,
      integer   const * nrhs,
      real      const   AB[],
      integer   const * ldAB,
      integer   const   ipiv[],
      real              B[],
      integer   const * ldB,
      integer         * info
    );

    void
    LAPACK_F77NAME(dgbtrs)(
      character  const * TRANS,
      integer    const * N,
      integer    const * KL,
      integer    const * KU,
      integer    const * nrhs,
      doublereal const   AB[],
      integer    const * ldAB,
      integer    const   ipiv[],
      doublereal         B[],
      integer    const * ldB,
      integer          * info
    );
  }
  #endif

  inline
  integer
  gbtrs(
    Transposition const & TRANS,
    integer       N,
    integer       KL,
    integer       KU,
    integer       nrhs,
    real    const AB[],
    integer       ldAB,
    integer const ipiv[],
    real          B[],
    integer       ldB
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sgbtrs)(
      const_cast<character*>(trans_blas[TRANS]),
      &N, &KL, &KU, &nrhs,
      AB, &ldAB, ipiv, B, &ldB, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgbtrs(
      trans_blas[TRANS],
      &N, &KL, &KU, &nrhs,
      AB, &ldAB, ipiv, B, &ldB, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgbtrs)(
      const_cast<char*>(trans_blas[TRANS]),
      &N, &KL, &KU, &nrhs,
      const_cast<real*>(AB), &ldAB,
      const_cast<integer*>(ipiv),
      const_cast<real*>(B), &ldB, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  inline
  integer
  gbtrs(
    Transposition const & TRANS,
    integer          N,
    integer          KL,
    integer          KU,
    integer          nrhs,
    doublereal const AB[],
    integer          ldAB,
    integer    const ipiv[],
    doublereal       B[],
    integer          ldB
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dgbtrs)(
      const_cast<character*>(trans_blas[TRANS]),
      &N, &KL, &KU, &nrhs,
      AB, &ldAB, ipiv, B, &ldB, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgbtrs(
      trans_blas[TRANS],
      &N, &KL, &KU, &nrhs, AB, &ldAB, ipiv, B, &ldB, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgbtrs)(
      const_cast<char*>(trans_blas[TRANS]),
      &N, &KL, &KU, &nrhs,
      const_cast<doublereal*>(AB), &ldAB,
      const_cast<integer*>(ipiv),
      const_cast<doublereal*>(B), &ldB, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  /*\
   *
   *  Purpose
   *  =======
   *
   *  DPBTRF computes the Cholesky factorization of a real symmetric
   *  positive definite band matrix A.
   *
   *  The factorization has the form
   *     A = U**T * U,  if UPLO = 'U', or
   *     A = L  * L**T,  if UPLO = 'L',
   *  where U is an upper triangular matrix and L is lower triangular.
   *
   *  Arguments
   *  =========
   *
   *  UPLO    (input) CHARACTER*1
   *          = 'U':  Upper triangle of A is stored;
   *          = 'L':  Lower triangle of A is stored.
   *
   *  N       (input) INTEGER
   *          The order of the matrix A.  N >= 0.
   *
   *  KD      (input) INTEGER
   *          The number of superdiagonals of the matrix A if UPLO = 'U',
   *          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
   *
   *  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
   *          On entry, the upper or lower triangle of the symmetric band
   *          matrix A, stored in the first KD+1 rows of the array.  The
   *          j-th column of A is stored in the j-th column of the array AB
   *          as follows:
   *          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
   *          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
   *
   *          On exit, if INFO = 0, the triangular factor U or L from the
   *          Cholesky factorization A = U**T*U or A = L*L**T of the band
   *          matrix A, in the same storage format as A.
   *
   *  LDAB    (input) INTEGER
   *          The leading dimension of the array AB.  LDAB >= KD+1.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *          > 0:  if INFO = i, the leading minor of order i is not
   *                positive definite, and the factorization could not be
   *                completed.
   *
   *  Further Details
   *  ===============
   *
   *  The band storage scheme is illustrated by the following example, when
   *  N = 6, KD = 2, and UPLO = 'U':
   *
   *  On entry:                       On exit:
   *
   *      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
   *      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
   *     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
   *
   *  Similarly, if UPLO = 'L' the format of A is as follows:
   *
   *  On entry:                       On exit:
   *
   *     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
   *     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
   *     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
   *
   *  Array elements marked * are not used by the routine.
   *
   *  Contributed by
   *  Peter Mayes and Giuseppe Radicati, IBM ECSEC, Rome, March 23, 1989
   *
   *  =====================================================================
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {

    void
    LAPACK_F77NAME(spbtrf)(
      character const * UPLO,
      integer   const * N,
      integer   const * KD,
      real              AB[],
      integer   const * ldAB,
      integer         * info
    );

    void
    LAPACK_F77NAME(dpbtrf)(
      character const * UPLO,
      integer   const * N,
      integer   const * KD,
      doublereal        AB[],
      integer   const * ldAB,
      integer         * info
    );
  }
  #endif

  inline
  integer
  pbtrf(
    ULselect const & UPLO,
    integer          N,
    integer          KD,
    real             AB[],
    integer          ldAB
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(spbtrf)(
      const_cast<character*>(uplo_blas[UPLO]),
      &N, &KD, AB, &ldAB, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    spbtrf( uplo_blas[UPLO], &N, &KD, AB, &ldAB, &info );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(spbtrf)(
      const_cast<character*>(uplo_blas[UPLO]),
      &N, &KD, AB, &ldAB, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline
  integer
  pbtrf(
    ULselect const & UPLO,
    integer          N,
    integer          KD,
    doublereal       AB[],
    integer          ldAB
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dpbtrf)(
      const_cast<character*>(uplo_blas[UPLO]),
      &N, &KD, AB, &ldAB, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dpbtrf( uplo_blas[UPLO], &N, &KD, AB, &ldAB, &info );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dpbtrf)(
      const_cast<character*>(uplo_blas[UPLO]),
      &N, &KD, AB, &ldAB, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  /*\
   *
   *  Purpose
   *  =======
   *
   *  DPBTRS solves a system of linear equations A*X = B with a symmetric
   *  positive definite band matrix A using the Cholesky factorization
   *  A = U**T*U or A = L*L**T computed by DPBTRF.
   *
   *  Arguments
   *  =========
   *
   *  UPLO    (input) CHARACTER*1
   *          = 'U':  Upper triangular factor stored in AB;
   *          = 'L':  Lower triangular factor stored in AB.
   *
   *  N       (input) INTEGER
   *          The order of the matrix A.  N >= 0.
   *
   *  KD      (input) INTEGER
   *          The number of superdiagonals of the matrix A if UPLO = 'U',
   *          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
   *
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of columns
   *          of the matrix B.  NRHS >= 0.
   *
   *  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
   *          The triangular factor U or L from the Cholesky factorization
   *
   *          A = U**T*U or A = L*L**T of the band matrix A, stored in the
   *
   *          first KD+1 rows of the array.  The j-th column of U or L is
   *          stored in the j-th column of the array AB as follows:
   *          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;
   *
   *          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd).
   *
   *
   *  LDAB    (input) INTEGER
   *          The leading dimension of the array AB.  LDAB >= KD+1.
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
   *          On entry, the right hand side matrix B.
   *          On exit, the solution matrix X.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B.  LDB >= max(1,N).
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *
   *  =====================================================================
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {

    void
    LAPACK_F77NAME(spbtrs)(
      character const * UPLO,
      integer   const * N,
      integer   const * KD,
      integer   const * nrhs,
      real      const   AB[],
      integer   const * ldAB,
      real              B[],
      integer   const * ldB,
      integer         * info
    );

    void
    LAPACK_F77NAME(dpbtrs)(
      character  const * UPLO,
      integer    const * N,
      integer    const * KD,
      integer    const * nrhs,
      doublereal const   AB[],
      integer    const * ldAB,
      doublereal         B[],
      integer    const * ldB,
      integer          * info
    );
  }
  #endif

  inline
  integer
  pbtrs(
    ULselect const & UPLO,
    integer          N,
    integer          KD,
    integer          nrhs,
    real     const   AB[],
    integer          ldAB,
    real             B[],
    integer          ldB
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(spbtrs)(
      const_cast<character*>(uplo_blas[UPLO]),
      &N, &KD, &nrhs, AB, &ldAB, B, &ldB, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    spbtrs( uplo_blas[UPLO], &N, &KD, &nrhs, AB, &ldAB, B, &ldB, &info );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(spbtrs)(
      const_cast<character*>(uplo_blas[UPLO]),
      &N, &KD, &nrhs,
      const_cast<real*>(AB), &ldAB,
      const_cast<real*>(B), &ldB, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline
  integer
  pbtrs(
    ULselect   const & UPLO,
    integer            N,
    integer            KD,
    integer            nrhs,
    doublereal const   AB[],
    integer            ldAB,
    doublereal         B[],
    integer            ldB
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dpbtrs)(
      const_cast<character*>(uplo_blas[UPLO]),
      &N, &KD, &nrhs, AB, &ldAB, B, &ldB, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dpbtrs( uplo_blas[UPLO], &N, &KD, &nrhs, AB, &ldAB, B, &ldB, &info );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dpbtrs)(
      const_cast<character*>(uplo_blas[UPLO]),
      &N, &KD, &nrhs,
      const_cast<doublereal*>(AB), &ldAB,
      const_cast<doublereal*>(B), &ldB, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

}
