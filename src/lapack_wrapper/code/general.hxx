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
/// file: general.hxx
///

/*
//    ____                           _   __  __       _        _
//   / ___| ___ _ __   ___ _ __ __ _| | |  \/  | __ _| |_ _ __(_)_  __
//  | |  _ / _ \ '_ \ / _ \ '__/ _` | | | |\/| |/ _` | __| '__| \ \/ /
//  | |_| |  __/ | | |  __/ | | (_| | | | |  | | (_| | |_| |  | |>  <
//   \____|\___|_| |_|\___|_|  \__,_|_| |_|  |_|\__,_|\__|_|  |_/_/\_\
*/

namespace lapack_wrapper {

  /*
  //    __ _  ___ _ __
  //   / _` |/ _ \ '__|
  //  | (_| |  __/ |
  //   \__, |\___|_|
  //   |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGER   performs the rank 1 operation
   *
   *     A := alpha*x*y' + A,
   *
   *  where alpha is a scalar, x is an m element vector, y is an n element
   *  vector and A is an m by n matrix.
   *
   *  Parameters
   *  ==========
   *
   *  M      - INTEGER.
   *           On entry, M specifies the number of rows of the matrix A.
   *           M must be at least zero.
   *           Unchanged on exit.
   *
   *  N      - INTEGER.
   *           On entry, N specifies the number of columns of the matrix A.
   *           N must be at least zero.
   *           Unchanged on exit.
   *
   *  ALPHA  - DOUBLE PRECISION.
   *           On entry, ALPHA specifies the scalar alpha.
   *           Unchanged on exit.
   *
   *  X      - DOUBLE PRECISION array of dimension at least
   *           ( 1 + ( m - 1 )*abs( INCX ) ).
   *           Before entry, the incremented array X must contain the m
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
   *           Before entry, the leading m by n part of the array A must
   *           contain the matrix of coefficients. On exit, A is
   *           overwritten by the updated matrix.
   *
   *  LDA    - INTEGER.
   *           On entry, LDA specifies the first dimension of A as declared
   *           in the calling (sub) program. LDA must be at least
   *           max( 1, m ).
   *           Unchanged on exit.
   *
   *  Level 2 Blas routine.
  \*/
  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  extern "C" {
    void
    BLASFUNC(sger)(
      integer const * M,
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
    BLASFUNC(dger)(
      integer    const * M,
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
  ger(
    integer    M,
    integer    N,
    real       ALPHA,
    real const X[],
    integer    INCX,
    real const Y[],
    integer    INCY,
    real       A[],
    integer    LDA
  ) {
  #if defined(LAPACK_WRAPPER_USE_MKL)
    sger( &M, &N, &ALPHA, X, &INCX, Y, &INCY, A, &LDA );
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)
    BLASFUNC(sger)(
      &M, &N, &ALPHA,
      const_cast<real*>(X), &INCX,
      const_cast<real*>(Y), &INCY,
      A, &LDA
    );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    CBLASNAME(sger)( CblasColMajor, M, N, ALPHA, X, INCX, Y, INCY, A, LDA );
  #else
  #error "LapackWrapper undefined mapping!"
  #endif
  }

  inline
  void
  ger(
    integer          M,
    integer          N,
    doublereal       ALPHA,
    doublereal const X[],
    integer          INCX,
    doublereal const Y[],
    integer          INCY,
    doublereal       A[],
    integer          LDA
  ) {
  #if defined(LAPACK_WRAPPER_USE_MKL)
    dger( &M, &N, &ALPHA, X, &INCX, Y, &INCY, A, &LDA );
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)
    BLASFUNC(dger)(
      &M, &N, &ALPHA,
      const_cast<doublereal*>(X), &INCX,
      const_cast<doublereal*>(Y), &INCY,
      A, &LDA
    );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    CBLASNAME(dger)( CblasColMajor, M, N, ALPHA, X, INCX, Y, INCY, A, LDA );
  #else
  #error "LapackWrapper undefined mapping!"
  #endif
  }

  /*
  //    __ _  ___ _ __ _____   __
  //   / _` |/ _ \ '_ ` _ \ \ / /
  //  | (_| |  __/ | | | | \ V /
  //   \__, |\___|_| |_| |_|\_/
  //   |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGEMV  performs one of the matrix-vector operations
   *
   *     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
   *
   *  where alpha and beta are scalars, x and y are vectors and A is an
   *  m by n matrix.
   *
   *  Parameters
   *  ==========
   *
   *  TRANS  - CHARACTER*1.
   *           On entry, TRANS specifies the operation to be performed as
   *           follows:
   *
   *              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
   *              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
   *              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
   *
   *           Unchanged on exit.
   *
   *  M      - INTEGER.
   *           On entry, M specifies the number of rows of the matrix A.
   *           M must be at least zero.
   *           Unchanged on exit.
   *
   *  N      - INTEGER.
   *           On entry, N specifies the number of columns of the matrix A.
   *           N must be at least zero.
   *           Unchanged on exit.
   *
   *  ALPHA  - DOUBLE PRECISION.
   *           On entry, ALPHA specifies the scalar alpha.
   *           Unchanged on exit.
   *
   *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
   *           Before entry, the leading m by n part of the array A must
   *           contain the matrix of coefficients.
   *           Unchanged on exit.
   *
   *  LDA    - INTEGER.
   *           On entry, LDA specifies the first dimension of A as declared
   *           in the calling (sub) program. LDA must be at least
   *           max( 1, m ).
   *           Unchanged on exit.
   *
   *  X      - DOUBLE PRECISION array of DIMENSION at least
   *           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
   *           and at least
   *           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
   *           Before entry, the incremented array X must contain the
   *           vector x.
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
   *  Y      - DOUBLE PRECISION array of DIMENSION at least
   *           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
   *           and at least
   *           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
   *           Before entry with BETA non-zero, the incremented array Y
   *           must contain the vector y. On exit, Y is overwritten by the
   *           updated vector y.
   *
   *  INCY   - INTEGER.
   *           On entry, INCY specifies the increment for the elements of
   *           Y. INCY must not be zero.
   *           Unchanged on exit.
   *
   *
   *  Level 2 Blas routine.
  \*/
  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  extern "C" {
    void
    BLASFUNC(sgemv)(
      character const   TRANS[],
      integer   const * M,
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
    BLASFUNC(dgemv)(
      character  const   TRANS[],
      integer    const * M,
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
  gemv(
    Transposition const & TRANS,
    integer               M,
    integer               N,
    real                  ALPHA,
    real const            A[],
    integer               LDA,
    real const            X[],
    integer               INCX,
    real                  BETA,
    real                  Y[],
    integer               INCY
  ) {
  #if defined(LAPACK_WRAPPER_USE_MKL)
    sgemv(
      trans_blas[TRANS],
      &M, &N, &ALPHA, A, &LDA, X, &INCX, &BETA, Y, &INCY
    );
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)
    BLASFUNC(sgemv)(
      const_cast<character*>(trans_blas[TRANS]),
      &M, &N, &ALPHA,
      const_cast<real*>(A), &LDA,
      const_cast<real*>(X), &INCX,
      &BETA, Y, &INCY
    );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    CBLASNAME(sgemv)(
      CblasColMajor, trans_cblas[TRANS],
      M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY
    );
  #else
  #error "LapackWrapper undefined mapping!"
  #endif
  }

  inline
  void
  gemv(
    Transposition const & TRANS,
    integer               M,
    integer               N,
    doublereal            ALPHA,
    doublereal const      A[],
    integer               LDA,
    doublereal const      X[],
    integer               INCX,
    doublereal            BETA,
    doublereal            Y[],
    integer               INCY
  ) {
  #if defined(LAPACK_WRAPPER_USE_MKL)
    dgemv(
      trans_blas[TRANS],
      &M, &N, &ALPHA, A, &LDA, X, &INCX, &BETA, Y, &INCY
    );
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)
    BLASFUNC(dgemv)(
      const_cast<character*>(trans_blas[TRANS]),
      &M, &N, &ALPHA,
      const_cast<doublereal*>(A), &LDA,
      const_cast<doublereal*>(X), &INCX,
      &BETA, Y, &INCY
    );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    CBLASNAME(dgemv)(
      CblasColMajor, trans_cblas[TRANS],
      M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY
    );
  #else
  #error "LapackWrapper undefined mapping!"
  #endif
  }

  /*
  //    __ _  ___ _ __ ___  _ __ ___
  //   / _` |/ _ \ '_ ` _ \| '_ ` _ \
  //  | (_| |  __/ | | | | | | | | | |
  //   \__, |\___|_| |_| |_|_| |_| |_|
  //   |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGEMM  performs one of the matrix-matrix operations
   *
   *     C := alpha*op( A )*op( B ) + beta*C,
   *
   *  where  op( X ) is one of
   *
   *     op( X ) = X   or   op( X ) = X',
   *
   *  alpha and beta are scalars, and A, B and C are matrices, with op( A )
   *  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
   *
   *  Parameters
   *  ==========
   *
   *  TRANSA - CHARACTER*1.
   *           On entry, TRANSA specifies the form of op( A ) to be used in
   *           the matrix multiplication as follows:
   *
   *              TRANSA = 'N' or 'n',  op( A ) = A.
   *              TRANSA = 'T' or 't',  op( A ) = A'.
   *              TRANSA = 'C' or 'c',  op( A ) = A'.
   *
   *           Unchanged on exit.
   *
   *  TRANSB - CHARACTER*1.
   *           On entry, TRANSB specifies the form of op( B ) to be used in
   *           the matrix multiplication as follows:
   *
   *              TRANSB = 'N' or 'n',  op( B ) = B.
   *              TRANSB = 'T' or 't',  op( B ) = B'.
   *              TRANSB = 'C' or 'c',  op( B ) = B'.
   *
   *           Unchanged on exit.
   *
   *  M      - INTEGER.
   *           On entry,  M  specifies  the number  of rows  of the  matrix
   *           op( A )  and of the  matrix  C.  M  must  be at least  zero.
   *           Unchanged on exit.
   *
   *  N      - INTEGER.
   *           On entry,  N  specifies the number  of columns of the matrix
   *           op( B ) and the number of columns of the matrix C. N must be
   *           at least zero.
   *           Unchanged on exit.
   *
   *  K      - INTEGER.
   *           On entry,  K  specifies  the number of columns of the matrix
   *           op( A ) and the number of rows of the matrix op( B ). K must
   *           be at least  zero.
   *           Unchanged on exit.
   *
   *  ALPHA  - DOUBLE PRECISION.
   *           On entry, ALPHA specifies the scalar alpha.
   *           Unchanged on exit.
   *
   *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
   *           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
   *           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
   *           part of the array  A  must contain the matrix  A,  otherwise
   *           the leading  k by m  part of the array  A  must contain  the
   *           matrix A.
   *           Unchanged on exit.
   *
   *  LDA    - INTEGER.
   *           On entry, LDA specifies the first dimension of A as declared
   *           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
   *           LDA must be at least  max( 1, m ), otherwise  LDA must be at
   *           least  max( 1, k ).
   *           Unchanged on exit.
   *
   *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
   *           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
   *           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
   *           part of the array  B  must contain the matrix  B,  otherwise
   *           the leading  n by k  part of the array  B  must contain  the
   *           matrix B.
   *           Unchanged on exit.
   *
   *  LDB    - INTEGER.
   *           On entry, LDB specifies the first dimension of B as declared
   *           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
   *           LDB must be at least  max( 1, k ), otherwise  LDB must be at
   *           least  max( 1, n ).
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
   *           On exit, the array  C  is overwritten by the  m by n  matrix
   *           ( alpha*op( A )*op( B ) + beta*C ).
   *
   *  LDC    - INTEGER.
   *           On entry, LDC specifies the first dimension of C as declared
   *           in  the  calling  (sub)  program.   LDC  must  be  at  least
   *           max( 1, m ).
   *           Unchanged on exit.
   *
   *
   *  Level 3 Blas routine.
  \*/
  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  extern "C" {
    void
    BLASFUNC(sgemm)(
      character const   TRANSA[],
      character const   TRANSB[],
      integer   const * M,
      integer   const * N,
      integer   const * K,
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
    BLASFUNC(dgemm)(
      character  const   TRANSA[],
      character  const   TRANSB[],
      integer    const * M,
      integer    const * N,
      integer    const * K,
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
  gemm(
    Transposition const & TRANSA,
    Transposition const & TRANSB,
    integer               M,
    integer               N,
    integer               K,
    real                  ALPHA,
    real            const A[],
    integer               LDA,
    real            const B[],
    integer               LDB,
    real                  BETA,
    real                  C[],
    integer               LDC
  ) {
  #if defined(LAPACK_WRAPPER_USE_MKL)
    sgemm(
      trans_blas[TRANSA], trans_blas[TRANSB],
      &M, &N, &K, &ALPHA, A, &LDA, B, &LDB,
      &BETA, C, &LDC
    );
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)
    BLASFUNC(sgemm)(
      const_cast<character*>(trans_blas[TRANSA]),
      const_cast<character*>(trans_blas[TRANSB]),
      &M, &N, &K,
      &ALPHA, const_cast<real*>(A), &LDA,
      const_cast<real*>(B), &LDB,
      &BETA, C, &LDC
    );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    CBLASNAME(sgemm)(
      CblasColMajor,
      trans_cblas[TRANSA],
      trans_cblas[TRANSB],
      M, N, K,
      ALPHA, A, LDA,
      B, LDB,
      BETA, C, LDC
    );
  #else
  #error "LapackWrapper undefined mapping!"
  #endif
  }

  inline
  void
  gemm(
    Transposition const & TRANSA,
    Transposition const & TRANSB,
    integer               M,
    integer               N,
    integer               K,
    doublereal            ALPHA,
    doublereal const      A[],
    integer               LDA,
    doublereal const      B[],
    integer               LDB,
    doublereal            BETA,
    doublereal            C[],
    integer               LDC
  ) {
  #if defined(LAPACK_WRAPPER_USE_MKL)
    dgemm(
      trans_blas[TRANSA], trans_blas[TRANSB],
      &M, &N, &K, &ALPHA, A, &LDA, B, &LDB,
      &BETA, C, &LDC
    );
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)
    BLASFUNC(dgemm)(
      const_cast<character*>(trans_blas[TRANSA]),
      const_cast<character*>(trans_blas[TRANSB]),
      &M, &N, &K,
      &ALPHA, const_cast<doublereal*>(A), &LDA,
      const_cast<doublereal*>(B), &LDB,
      &BETA, C, &LDC
    );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    CBLASNAME(dgemm)(
      CblasColMajor,
      trans_cblas[TRANSA],
      trans_cblas[TRANSB],
      M, N, K,
      ALPHA, A, LDA,
      B, LDB,
      BETA, C, LDC
    );
  #else
  #error "LapackWrapper undefined mapping!"
  #endif
  }

  /*
  //              _         __
  //    __ _  ___| |_ _ __ / _|
  //   / _` |/ _ \ __| '__| |_
  //  | (_| |  __/ |_| |  |  _|
  //   \__, |\___|\__|_|  |_|
  //   |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGETRF computes an LU factorization of a general M-by-N matrix A
   *  using partial pivoting with row interchanges.
   *
   *  The factorization has the form
   *     A = P * L * U
   *  where P is a permutation matrix, L is lower triangular with unit
   *  diagonal elements (lower trapezoidal if m > n), and U is upper
   *  triangular (upper trapezoidal if m < n).
   *
   *  This is the right-looking Level 3 BLAS version of the algorithm.
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
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the M-by-N matrix to be factored.
   *          On exit, the factors L and U from the factorization
   *          A = P*L*U; the unit diagonal elements of L are not stored.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  IPIV    (output) INTEGER array, dimension (min(M,N))
   *          The pivot indices; for 1 <= i <= min(M,N), row i of the
   *          matrix was interchanged with row IPIV(i).
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
   *                has been completed, but the factor U is exactly
   *                singular, and division by zero will occur if it is used
   *                to solve a system of equations.
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {
    void
    BLASFUNC(sgetrf)(
      integer const * N,
      integer const * M,
      real            A[],
      integer const * LDA,
      integer         IPIV[],
      integer       * INFO
    );

    void
    BLASFUNC(dgetrf)(
      integer const * N,
      integer const * M,
      doublereal      A[],
      integer const * LDA,
      integer         IPIV[],
      integer       * INFO
    );
  }
  #endif

  inline
  integer
  getrf(
    integer N,
    integer M,
    real    A[],
    integer LDA,
    integer IPIV[]
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_OPENBLAS)
    LAPACK_F77NAME(sgetrf)( &N, &M, A, &LDA, IPIV, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_LAPACK) || \
          defined(LAPACK_WRAPPER_USE_ATLAS)
    BLASFUNC(sgetrf)( &N, &M, A, &LDA, IPIV, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgetrf( &N, &M, A, &LDA, IPIV, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgetrf)( &N, &M, A, &LDA, IPIV, &INFO );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  inline
  integer
  getrf(
    integer    N,
    integer    M,
    doublereal A[],
    integer    LDA,
    integer    IPIV[]
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_OPENBLAS)
    LAPACK_F77NAME(dgetrf)( &N, &M, A, &LDA, IPIV, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_LAPACK) || \
          defined(LAPACK_WRAPPER_USE_ATLAS)
    BLASFUNC(dgetrf)( &N, &M, A, &LDA, IPIV, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgetrf( &N, &M, A, &LDA, IPIV, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgetrf)( &N, &M, A, &LDA, IPIV, &INFO );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  /*
  //              _
  //    __ _  ___| |_ _ __ ___
  //   / _` |/ _ \ __| '__/ __|
  //  | (_| |  __/ |_| |  \__ \
  //   \__, |\___|\__|_|  |___/
  //   |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGETRS solves a system of linear equations
   *     A * X = B  or  A**T * X = B
   *  with a general N-by-N matrix A using the LU factorization computed
   *  by DGETRF.
   *
   *  Arguments
   *  =========
   *
   *  TRANS   (input) CHARACTER*1
   *          Specifies the form of the system of equations:
   *          = 'N':  A * X = B  (No transpose)
   *          = 'T':  A**T* X = B  (Transpose)
   *          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
   *
   *  N       (input) INTEGER
   *          The order of the matrix A.  N >= 0.
   *
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of columns
   *          of the matrix B.  NRHS >= 0.
   *
   *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
   *          The factors L and U from the factorization A = P*L*U
   *          as computed by DGETRF.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,N).
   *
   *  IPIV    (input) INTEGER array, dimension (N)
   *          The pivot indices from DGETRF; for 1<=i<=N, row i of the
   *          matrix was interchanged with row IPIV(i).
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
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {
    void
    BLASFUNC(sgetrs)(
      character const   TRANS[],
      integer   const * N,
      integer   const * NRHS,
      real      const   A[],
      integer   const * LDA,
      integer   const   IPIV[],
      real              B[],
      integer   const * LDB,
      integer         * INFO
    );

    void
    BLASFUNC(dgetrs)(
      character  const   TRANS[],
      integer    const * N,
      integer    const * NRHS,
      doublereal const   A[],
      integer    const * LDA,
      integer    const   IPIV[],
      doublereal         B[],
      integer    const * LDB,
      integer          * INFO
    );
  }
  #endif

  inline
  integer
  getrs(
    Transposition const & TRANS,
    integer               N,
    integer               NRHS,
    real          const   A[],
    integer               LDA,
    integer       const   IPIV[],
    real                  B[],
    integer               LDB
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_OPENBLAS)
    LAPACK_F77NAME(sgetrs)(
      const_cast<character*>(trans_blas[TRANS]),
      &N, &NRHS,
      const_cast<real*>(A), &LDA,
      const_cast<integer*>(IPIV),
      B, &LDB, &INFO
    );
    #elif defined(LAPACK_WRAPPER_USE_LAPACK) || \
          defined(LAPACK_WRAPPER_USE_ATLAS)
    BLASFUNC(sgetrs)(
      const_cast<character*>(trans_blas[TRANS]),
      &N, &NRHS,
      const_cast<real*>(A), &LDA,
      const_cast<integer*>(IPIV),
      B, &LDB, &INFO
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgetrs( trans_blas[TRANS], &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgetrs)(
      const_cast<character*>(trans_blas[TRANS]),
      &N, &NRHS,
      const_cast<real*>(A), &LDA,
      const_cast<integer*>(IPIV),
      B, &LDB, &INFO
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  inline
  integer
  getrs(
    Transposition const & TRANS,
    integer               N,
    integer               NRHS,
    doublereal    const   A[],
    integer               LDA,
    integer       const   IPIV[],
    doublereal            B[],
    integer               LDB
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_OPENBLAS)
    LAPACK_F77NAME(dgetrs)(
      const_cast<character*>(trans_blas[TRANS]),
      &N, &NRHS,
      const_cast<doublereal*>(A), &LDA,
      const_cast<integer*>(IPIV),
      B, &LDB, &INFO
    );
    #elif defined(LAPACK_WRAPPER_USE_LAPACK) || \
          defined(LAPACK_WRAPPER_USE_ATLAS)
    BLASFUNC(dgetrs)(
      const_cast<character*>(trans_blas[TRANS]),
      &N, &NRHS,
      const_cast<doublereal*>(A), &LDA,
      const_cast<integer*>(IPIV),
      B, &LDB, &INFO
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgetrs( trans_blas[TRANS], &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgetrs)(
      const_cast<character*>(trans_blas[TRANS]),
      &N, &NRHS,
      const_cast<doublereal*>(A), &LDA,
      const_cast<integer*>(IPIV),
      B, &LDB, &INFO
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  /*
  //    __ _  ___  _____   __
  //   / _` |/ _ \/ __\ \ / /
  //  | (_| |  __/\__ \\ V /
  //   \__, |\___||___/ \_/
  //   |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGESV computes the solution to a real system of linear equations
   *     A * X = B,
   *  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
   *
   *  The LU decomposition with partial pivoting and row interchanges is
   *  used to factor A as
   *     A = P * L * U,
   *  where P is a permutation matrix, L is unit lower triangular, and U is
   *  upper triangular.  The factored form of A is then used to solve the
   *  system of equations A * X = B.
   *
   *  Arguments
   *  =========
   *
   *  N       (input) INTEGER
   *          The number of linear equations, i.e., the order of the
   *          matrix A.  N >= 0.
   *
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of columns
   *          of the matrix B.  NRHS >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the N-by-N coefficient matrix A.
   *          On exit, the factors L and U from the factorization
   *          A = P*L*U; the unit diagonal elements of L are not stored.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,N).
   *
   *  IPIV    (output) INTEGER array, dimension (N)
   *          The pivot indices that define the permutation matrix P;
   *          row i of the matrix was interchanged with row IPIV(i).
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
   *          On entry, the N-by-NRHS matrix of right hand side matrix B.
   *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B.  LDB >= max(1,N).
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
   *                has been completed, but the factor U is exactly
   *                singular, so the solution could not be computed.
   *
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {
    void
    BLASFUNC(sgesv)(
      integer const * N,
      integer const * NRHS,
      real            A[],
      integer const * LDA,
      integer         IPIV[],
      real            B[],
      integer const * LDB,
      integer       * INFO
    );

    void
    BLASFUNC(dgesv)(
      integer const * N,
      integer const * NRHS,
      doublereal      A[],
      integer const * LDA,
      integer         IPIV[],
      doublereal      B[],
      integer const * LDB,
      integer       * INFO
    );
  }
  #endif

  inline
  integer
  gesv(
    integer N,
    integer NRHS,
    real    A[],
    integer LDA,
    integer IPIV[],
    real    B[],
    integer LDB
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_OPENBLAS)
    LAPACK_F77NAME(sgesv)( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_LAPACK) || \
          defined(LAPACK_WRAPPER_USE_ATLAS)
    BLASFUNC(sgesv)( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgesv( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgesv)( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  inline
  integer
  gesv(
    integer    N,
    integer    NRHS,
    doublereal A[],
    integer    LDA,
    integer    IPIV[],
    doublereal B[],
    integer    LDB
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_OPENBLAS)
    LAPACK_F77NAME(dgesv)( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_LAPACK) || \
          defined(LAPACK_WRAPPER_USE_ATLAS)
    BLASFUNC(dgesv)( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgesv( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgesv)( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  /*
  //    ____ _____ _____ ____ ____
  //   / ___| ____|_   _/ ___|___ \
  //  | |  _|  _|   | || |     __) |
  //  | |_| | |___  | || |___ / __/
  //   \____|_____| |_| \____|_____|
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGETC2 computes an LU factorization with complete pivoting of the
   *  n-by-n matrix A. The factorization has the form A = P * L * U * Q,
   *  where P and Q are permutation matrices, L is lower triangular with
   *  unit diagonal elements and U is upper triangular.
   *
   *  This is the Level 2 BLAS algorithm.
   *
   *  Arguments
   *  =========
   *
   *  N       (input) INTEGER
   *          The order of the matrix A. N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
   *          On entry, the n-by-n matrix A to be factored.
   *          On exit, the factors L and U from the factorization
   *          A = P*L*U*Q; the unit diagonal elements of L are not stored.
   *          If U(k, k) appears to be less than SMIN, U(k, k) is given the
   *          value of SMIN, i.e., giving a nonsingular perturbed system.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,N).
   *
   *  IPIV    (output) INTEGER array, dimension(N).
   *          The pivot indices; for 1 <= i <= N, row i of the
   *          matrix has been interchanged with row IPIV(i).
   *
   *  JPIV    (output) INTEGER array, dimension(N).
   *          The pivot indices; for 1 <= j <= N, column j of the
   *          matrix has been interchanged with column JPIV(j).
   *
   *  INFO    (output) INTEGER
   *           = 0: successful exit
   *           > 0: if INFO = k, U(k, k) is likely to produce owerflow if
   *                we try to solve for x in Ax = b. So U is perturbed to
   *                avoid the overflow.
   *
   *  Further Details
   *  ===============
   *
   *  Based on contributions by
   *     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
   *     Umea University, S-901 87 Umea, Sweden.
   *
   *  =====================================================================
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  extern "C" {

    void
    LAPACK_F77NAME(sgetc2)(
      integer const * N,
      real          * A,
      integer const * LDA,
      integer       * IPIV,
      integer       * JPIV,
      integer       * INFO
    );

    void
    LAPACK_F77NAME(dgetc2)(
      integer    const * N,
      doublereal       * A,
      integer    const * LDA,
      integer          * IPIV,
      integer          * JPIV,
      integer          * INFO
    );
  }
  #endif

  #if defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  template <typename T>
  integer
  getc2_tmpl(
    integer N,
    T       A[],
    integer LDA,
    integer IPIV[],
    integer JPIV[]
  );
  #endif

  inline
  integer
  getc2(
    integer   N,
    real      A[],
    integer   LDA,
    integer   IPIV[],
    integer   JPIV[]
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    INFO = getc2_tmpl<real>(N,A,LDA,IPIV,JPIV);
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgetc2( &N, A, &LDA, IPIV, JPIV, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgetc2)( &N, A, &LDA, IPIV, JPIV, &INFO );
    #else
    LAPACK_F77NAME(sgetc2)( &N, A, &LDA, IPIV, JPIV, &INFO );
    #endif
    return INFO;
  }

  inline
  integer
  getc2(
    integer    N,
    doublereal A[],
    integer    LDA,
    integer    IPIV[],
    integer    JPIV[]
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    INFO = getc2_tmpl<doublereal>(N,A,LDA,IPIV,JPIV);
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgetc2( &N, A, &LDA, IPIV, JPIV, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgetc2)( &N, A, &LDA, IPIV, JPIV, &INFO );
    #else
    LAPACK_F77NAME(dgetc2)( &N, A, &LDA, IPIV, JPIV, &INFO );
    #endif
    return INFO;
  }

  /*
  //    ____ _____ ____   ____ ____
  //   / ___| ____/ ___| / ___|___ \
  //  | |  _|  _| \___ \| |     __) |
  //  | |_| | |___ ___) | |___ / __/
  //   \____|_____|____/ \____|_____|
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGESC2 solves a system of linear equations
   *
   *            A * X = scale* RHS
   *
   *  with a general N-by-N matrix A using the LU factorization with
   *  complete pivoting computed by DGETC2.
   *
   *  Arguments
   *  =========
   *
   *  N       (input) INTEGER
   *          The order of the matrix A.
   *
   *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the  LU part of the factorization of the n-by-n
   *          matrix A computed by DGETC2:  A = P * L * U * Q
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1, N).
   *
   *  RHS     (input/output) DOUBLE PRECISION array, dimension (N).
   *          On entry, the right hand side vector b.
   *          On exit, the solution vector X.
   *
   *  IPIV    (input) INTEGER array, dimension (N).
   *          The pivot indices; for 1 <= i <= N, row i of the
   *          matrix has been interchanged with row IPIV(i).
   *
   *  JPIV    (input) INTEGER array, dimension (N).
   *          The pivot indices; for 1 <= j <= N, column j of the
   *          matrix has been interchanged with column JPIV(j).
   *
   *  SCALE    (output) DOUBLE PRECISION
   *           On exit, SCALE contains the scale factor. SCALE is chosen
   *           0 <= SCALE <= 1 to prevent owerflow in the solution.
   *
   *  Further Details
   *  ===============
   *
   *  Based on contributions by
   *     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
   *     Umea University, S-901 87 Umea, Sweden.
   *
   *  =====================================================================
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  extern "C" {

    void
    LAPACK_F77NAME(sgesc2)(
      integer const * N,
      real          * A,
      integer const * LDA,
      real          * RHS,
      integer const * IPIV,
      integer const * JPIV,
      real          * SCALE
    );

    void
    LAPACK_F77NAME(dgesc2)(
      integer const * N,
      doublereal    * A,
      integer const * LDA,
      doublereal    * RHS,
      integer const * IPIV,
      integer const * JPIV,
      doublereal    * SCALE
    );
  }
  #endif

  #if defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  template <typename T>
  T
  gesc2_tmpl(
    integer       N,
    T       const A[],
    integer       LDA,
    T             RHS[],
    integer const IPIV[],
    integer const JPIV[]
  );
  #endif

  inline
  real
  gesc2(
    integer       N,
    real    const A[],
    integer       LDA,
    real          RHS[],
    integer const IPIV[],
    integer const JPIV[]
  ) {
    real SCALE = 0;
    #if defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    SCALE = gesc2_tmpl<real>( N, A, LDA, RHS, IPIV, JPIV );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgesc2( &N, A, &LDA, RHS, IPIV, JPIV, &SCALE );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgesc2)(
      &N,
      const_cast<real*>(A), &LDA, RHS,
      const_cast<integer*>(IPIV),
      const_cast<integer*>(JPIV), &SCALE
    );
    #else
    LAPACK_F77NAME(sgesc2)(
      &N,
      const_cast<real*>(A), &LDA, RHS,
      const_cast<integer*>(IPIV),
      const_cast<integer*>(JPIV), &SCALE
    );
    #endif
    return SCALE;
  }

  inline
  doublereal
  gesc2(
    integer          N,
    doublereal const A[],
    integer          LDA,
    doublereal       RHS[],
    integer    const IPIV[],
    integer    const JPIV[]
  ) {
    doublereal SCALE = 0;
    #if defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    SCALE = gesc2_tmpl<doublereal>( N, A, LDA, RHS, IPIV, JPIV );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgesc2( &N, A, &LDA, RHS, IPIV, JPIV, &SCALE );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgesc2)(
      &N,
      const_cast<doublereal*>(A), &LDA, RHS,
      const_cast<integer*>(IPIV),
      const_cast<integer*>(JPIV), &SCALE
    );
    #else
    LAPACK_F77NAME(dgesc2)(
      &N,
      const_cast<doublereal*>(A), &LDA, RHS,
      const_cast<integer*>(IPIV),
      const_cast<integer*>(JPIV), &SCALE
    );
    #endif
    return SCALE;
  }

  /*
  //    ____ _____ ____ ___  _   _
  //   / ___| ____/ ___/ _ \| \ | |
  //  | |  _|  _|| |  | | | |  \| |
  //  | |_| | |__| |__| |_| | |\  |
  //   \____|_____\____\___/|_| \_|
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGECON estimates the reciprocal of the condition number of a general
   *  real matrix A, in either the 1-norm or the infinity-norm, using
   *  the LU factorization computed by DGETRF.
   *
   *  An estimate is obtained for norm(inv(A)), and the reciprocal of the
   *  condition number is computed as
   *     RCOND = 1 / ( norm(A) * norm(inv(A)) ).
   *
   *  Arguments
   *  =========
   *
   *  NORM    (input) CHARACTER*1
   *          Specifies whether the 1-norm condition number or the
   *          infinity-norm condition number is required:
   *          = '1' or 'O':  1-norm;
   *          = 'I':         Infinity-norm.
   *
   *  N       (input) INTEGER
   *          The order of the matrix A.  N >= 0.
   *
   *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
   *          The factors L and U from the factorization A = P*L*U
   *          as computed by DGETRF.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,N).
   *
   *  ANORM   (input) DOUBLE PRECISION
   *          If NORM = '1' or 'O', the 1-norm of the original matrix A.
   *          If NORM = 'I', the infinity-norm of the original matrix A.
   *
   *  RCOND   (output) DOUBLE PRECISION
   *          The reciprocal of the condition number of the matrix A,
   *          computed as RCOND = 1/(norm(A) * norm(inv(A))).
   *
   *  WORK    (workspace) DOUBLE PRECISION array, dimension (4*N)
   *
   *  IWORK   (workspace) INTEGER array, dimension (N)
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {

    void
    LAPACK_F77NAME(sgecon)(
      character const * NORM,
      integer   const * N,
      real      const * A,
      integer   const * ldA,
      real      const * ANORM,
      real            * RCOND,
      real            * WORK,
      integer         * IWORK,
      integer         * INFO
    );

    void
    LAPACK_F77NAME(dgecon)(
      character  const * NORM,
      integer    const * N,
      doublereal const * A,
      integer    const * ldA,
      doublereal const * ANORM,
      doublereal       * RCOND,
      doublereal       * WORK,
      integer          * IWORK,
      integer          * INFO
    );
  }
  #endif

  inline
  integer
  gecon1(
    integer    N,
    real const A[],
    integer    LDA,
    real       anorm,
    real     & rcond,
    real       work[],
    integer    iwork[]
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)
    LAPACK_F77NAME(sgecon)( "1", &N, A, &LDA, &anorm, &rcond, work, iwork, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
          defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sgecon)(
      const_cast<character*>("1"),
      &N, A, &LDA, &anorm, &rcond, work, iwork, &INFO
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgecon( "1", &N, A, &LDA, &anorm, &rcond, work, iwork, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgecon)(
      const_cast<character*>("1"), &N,
      const_cast<real*>(A),
      &LDA, &anorm, &rcond, work, iwork, &INFO
    );
    #else
    LAPACK_F77NAME(sgecon)( "1", &N, A, &LDA, &anorm, &rcond, work, iwork, &INFO);
    #endif
    return INFO;
  }

  inline
  integer
  gecon1(
    integer          N,
    doublereal const A[],
    integer          LDA,
    doublereal       anorm,
    doublereal     & rcond,
    doublereal       work[],
    integer          iwork[]
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)
    LAPACK_F77NAME(dgecon)( "1", &N, A, &LDA, &anorm, &rcond, work, iwork, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
          defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dgecon)(
      const_cast<character*>("1"),
      &N, A, &LDA, &anorm, &rcond, work, iwork, &INFO
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgecon( "1", &N, A, &LDA, &anorm, &rcond, work, iwork, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgecon)(
      const_cast<character*>("1"), &N,
      const_cast<doublereal*>(A),
      &LDA, &anorm, &rcond, work, iwork, &INFO
    );
    #else
    LAPACK_F77NAME(dgecon)( "1", &N, A, &LDA, &anorm, &rcond, work, iwork, &INFO);
    #endif
    return INFO;
  }

  inline
  integer
  geconInf(
    integer    N,
    real const A[],
    integer    LDA,
    real       anorm,
    real     & rcond,
    real       work[],
    integer    iwork[]
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)
    LAPACK_F77NAME(sgecon)( "I", &N, A, &LDA, &anorm, &rcond, work, iwork, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
          defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sgecon)(
      const_cast<character*>("I"),
      &N, A, &LDA, &anorm, &rcond, work, iwork, &INFO
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgecon( "I", &N, A, &LDA, &anorm, &rcond, work, iwork, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgecon)(
      const_cast<character*>("I"), &N,
      const_cast<real*>(A),
      &LDA, &anorm, &rcond, work, iwork, &INFO
    );
    #else
    LAPACK_F77NAME(sgecon)( "I", &N, A, &LDA, &anorm, &rcond, work, iwork, &INFO);
    #endif
    return INFO;
  }

  inline
  integer
  geconInf(
    integer          N,
    doublereal const A[],
    integer          LDA,
    doublereal       anorm,
    doublereal     & rcond,
    doublereal       work[],
    integer          iwork[]
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)
    LAPACK_F77NAME(dgecon)( "I", &N, A, &LDA, &anorm, &rcond, work, iwork, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
          defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dgecon)(
      const_cast<character*>("I"),
      &N, A, &LDA, &anorm, &rcond, work, iwork, &INFO
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgecon( "I", &N, A, &LDA, &anorm, &rcond, work, iwork, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgecon)(
      const_cast<character*>("I"), &N,
      const_cast<doublereal*>(A),
      &LDA, &anorm, &rcond, work, iwork, &INFO
    );
    #else
    LAPACK_F77NAME(dgecon)( "I", &N, A, &LDA, &anorm, &rcond, work, iwork, &INFO);
    #endif
    return INFO;
  }

  /*\
   *  Purpose
   *  =======
   *
   *  DGEEQU computes row and column scalings intended to equilibrate an
   *  M-by-N matrix A and reduce its condition number.  R returns the row
   *  scale factors and C the column scale factors, chosen to try to make
   *  the largest element in each row and column of the matrix B with
   *  elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1.
   *
   *  R(i) and C(j) are restricted to be between SMLNUM = smallest safe
   *  number and BIGNUM = largest safe number.  Use of these scaling
   *  factors is not guaranteed to reduce the condition number of A but
   *  works well in practice.
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
   *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
   *          The M-by-N matrix whose equilibration factors are
   *          to be computed.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  R       (output) DOUBLE PRECISION array, dimension (M)
   *          If INFO = 0 or INFO > M, R contains the row scale factors for A.
   *
   *  C       (output) DOUBLE PRECISION array, dimension (N)
   *          If INFO = 0,  C contains the column scale factors for A.
   *
   *  ROWCND  (output) DOUBLE PRECISION
   *          If INFO = 0 or INFO > M, ROWCND contains the ratio of the
   *          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and
   *          AMAX is neither too large nor too small, it is not worth
   *          scaling by R.
   *
   *  COLCND  (output) DOUBLE PRECISION
   *          If INFO = 0, COLCND contains the ratio of the smallest
   *          C(i) to the largest C(i).  If COLCND >= 0.1, it is not
   *          worth scaling by C.
   *
   *  AMAX    (output) DOUBLE PRECISION
   *          Absolute value of largest matrix element.  If AMAX is very
   *          close to overflow or very close to underflow, the matrix
   *          should be scaled.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *          > 0:  if INFO = i,  and i is
   *                <= M:  the i-th row of A is exactly zero
   *                >  M:  the (i-M)-th column of A is exactly zero
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {

    void
    LAPACK_F77NAME(sgeequ)(
      integer const * M,
      integer const * N,
      real    const * A,
      integer const * ldA,
      real          * R,
      real          * C,
      real          * ROWCND,
      real          * COLCND,
      real          * AMAX,
      integer       * INFO
    );

    void
    LAPACK_F77NAME(dgeequ)(
      integer    const * M,
      integer    const * N,
      doublereal const * A,
      integer    const * ldA,
      doublereal       * R,
      doublereal       * C,
      doublereal       * ROWCND,
      doublereal       * COLCND,
      doublereal       * AMAX,
      integer          * INFO
    );
  }
  #endif

  inline
  integer
  geequ(
    integer    M,
    integer    N,
    real const A[],
    integer    LDA,
    real       R[],
    real       C[],
    real     & ROWCND,
    real     & COLCND,
    real     & AMAX
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)  || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    LAPACK_F77NAME(sgeequ)(
      &M, &N,
      A, &LDA,
      R, C, &ROWCND, &COLCND, &AMAX, &INFO
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgeequ( &M, &N, A, &LDA, R, C, &ROWCND, &COLCND, &AMAX, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgeequ)(
      &M, &N,
      const_cast<real*>(A), &LDA,
      R, C, &ROWCND, &COLCND, &AMAX, &INFO
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  inline
  integer
  geequ(
    integer          M,
    integer          N,
    doublereal const A[],
    integer          LDA,
    doublereal       R[],
    doublereal       C[],
    doublereal     & ROWCND,
    doublereal     & COLCND,
    doublereal     & AMAX
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)  || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    LAPACK_F77NAME(dgeequ)(
      &M, &N,
      A, &LDA,
      R, C, &ROWCND, &COLCND, &AMAX, &INFO
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgeequ( &M, &N, A, &LDA, R, C, &ROWCND, &COLCND, &AMAX, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgeequ)(
      &M, &N,
      const_cast<doublereal*>(A), &LDA,
      R, C, &ROWCND, &COLCND, &AMAX, &INFO
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  /*\
   * SUBROUTINE DLAQGE( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, EQUED )
   *
   *  Purpose
   *  =======
   *
   *  DLAQGE equilibrates a general M by N matrix A using the row and
   *  column scaling factors in the vectors R and C.
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
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the M by N matrix A.
   *          On exit, the equilibrated matrix.  See EQUED for the form of
   *          the equilibrated matrix.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(M,1).
   *
   *  R       (input) DOUBLE PRECISION array, dimension (M)
   *          The row scale factors for A.
   *
   *  C       (input) DOUBLE PRECISION array, dimension (N)
   *          The column scale factors for A.
   *
   *  ROWCND  (input) DOUBLE PRECISION
   *          Ratio of the smallest R(i) to the largest R(i).
   *
   *  COLCND  (input) DOUBLE PRECISION
   *          Ratio of the smallest C(i) to the largest C(i).
   *
   *  AMAX    (input) DOUBLE PRECISION
   *          Absolute value of largest matrix entry.
   *
   *  EQUED   (output) CHARACTER*1
   *          Specifies the form of equilibration that was done.
   *          = 'N':  No equilibration
   *          = 'R':  Row equilibration, i.e., A has been premultiplied by
   *                  diag(R).
   *          = 'C':  Column equilibration, i.e., A has been postmultiplied
   *                  by diag(C).
   *          = 'B':  Both row and column equilibration, i.e., A has been
   *                  replaced by diag(R) * A * diag(C).
   *
   *  Internal Parameters
   *  ===================
   *
   *  THRESH is a threshold value used to decide if row or column scaling
   *  should be done based on the ratio of the row or column scaling
   *  factors.  If ROWCND < THRESH, row scaling is done, and if
   *  COLCND < THRESH, column scaling is done.
   *
   *  LARGE and SMALL are threshold values used to decide if row scaling
   *  should be done based on the absolute size of the largest matrix
   *  element.  If AMAX > LARGE or AMAX < SMALL, row scaling is done.
   *
   *  =====================================================================
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {

    void
    LAPACK_F77NAME(slaqge)(
      integer   const * M,
      integer   const * N,
      real            * A,
      integer   const * ldA,
      real      const * R,
      real      const * C,
      real      const * ROWCND,
      real      const * COLCND,
      real      const * AMAX,
      character const * EQUED
    );

    void
    LAPACK_F77NAME(dlaqge)(
      integer    const * M,
      integer    const * N,
      doublereal       * A,
      integer    const * ldA,
      doublereal const * R,
      doublereal const * C,
      doublereal const * ROWCND,
      doublereal const * COLCND,
      doublereal const * AMAX,
      character  const * EQUED
    );
  }
  #endif

  #if defined(LAPACK_WRAPPER_USE_OPENBLAS)
  extern "C" {

    void
    BLASFUNC(slaqge)(
      integer   const * M,
      integer   const * N,
      real            * A,
      integer   const * ldA,
      real      const * R,
      real      const * C,
      real      const * ROWCND,
      real      const * COLCND,
      real      const * AMAX,
      character const * EQUED
    );

    void
    BLASFUNC(dlaqge)(
      integer    const * M,
      integer    const * N,
      doublereal       * A,
      integer    const * ldA,
      doublereal const * R,
      doublereal const * C,
      doublereal const * ROWCND,
      doublereal const * COLCND,
      doublereal const * AMAX,
      character  const * EQUED
    );
  }
  #endif

  inline
  void
  laqge(
    integer             M,
    integer             N,
    real                A[],
    integer             LDA,
    real const          R[],
    real const          C[],
    real                ROWCND,
    real                COLCND,
    real                AMAX,
    EquilibrationType & equ
  ) {
    #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(slaqge)(
      &M, &N, A, &LDA, R, C, &ROWCND, &COLCND, &AMAX, equilibrate_blas[equ]
    );
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
    BLASFUNC(slaqge)(
      &M, &N, A, &LDA, R, C, &ROWCND, &COLCND, &AMAX, equilibrate_blas[equ]
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    slaqge(
      &M, &N, A, &LDA, R, C, &ROWCND, &COLCND, &AMAX,
      const_cast<character*>(equilibrate_blas[equ])
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(slaqge)(
      &M, &N, A, &LDA,
      const_cast<real*>(R),
      const_cast<real*>(C),
      &ROWCND, &COLCND, &AMAX,
      const_cast<character*>(equilibrate_blas[equ])
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
  }

  inline
  void
  laqge(
    integer             M,
    integer             N,
    doublereal          A[],
    integer             LDA,
    doublereal const    R[],
    doublereal const    C[],
    doublereal          ROWCND,
    doublereal          COLCND,
    doublereal          AMAX,
    EquilibrationType & equ
  ) {
    #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dlaqge)(
      &M, &N, A, &LDA, R, C, &ROWCND, &COLCND, &AMAX, equilibrate_blas[equ]
    );
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
    BLASFUNC(dlaqge)(
      &M, &N, A, &LDA, R, C, &ROWCND, &COLCND, &AMAX, equilibrate_blas[equ]
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dlaqge(
      &M, &N, A, &LDA, R, C, &ROWCND, &COLCND, &AMAX,
      const_cast<character*>(equilibrate_blas[equ])
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dlaqge)(
      &M, &N, A, &LDA,
      const_cast<doublereal*>(R),
      const_cast<doublereal*>(C),
      &ROWCND, &COLCND, &AMAX,
      const_cast<character*>(equilibrate_blas[equ])
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
  }

  /*
  //    __ _  ___  ___ ___  _ __  _   _
  //   / _` |/ _ \/ __/ _ \| '_ \| | | |
  //  | (_| |  __/ (_| (_) | |_) | |_| |
  //   \__, |\___|\___\___/| .__/ \__, |
  //   |___/               |_|    |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DLACPY copies all or part of a two-dimensional matrix A to another
   *  matrix B.
   *
   *  Arguments
   *  =========
   *
   *  UPLO    (input) CHARACTER*1
   *          Specifies the part of the matrix A to be copied to B.
   *          = 'U':      Upper triangular part
   *          = 'L':      Lower triangular part
   *          Otherwise:  All of the matrix A
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix A.  M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A.  N >= 0.
   *
   *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
   *          The m by n matrix A.  If UPLO = 'U', only the upper triangle
   *          or trapezoid is accessed; if UPLO = 'L', only the lower
   *          triangle or trapezoid is accessed.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  B       (output) DOUBLE PRECISION array, dimension (LDB,N)
   *          On exit, B = A in the locations specified by UPLO.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B.  LDB >= max(1,M).
   *
  \*/
  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  extern "C" {
    void
    LAPACK_F77NAME(slacpy)(
      character const   UPLO[],
      integer   const * M,
      integer   const * N,
      real      const   A[],
      integer   const * LDA,
      real              B[],
      integer   const * LDB
    );

    void
    LAPACK_F77NAME(dlacpy)(
      character  const   UPLO[],
      integer    const * M,
      integer    const * N,
      doublereal const   A[],
      integer    const * LDA,
      doublereal         B[],
      integer    const * LDB
    );
  }
  #endif

  inline
  integer
  gecopy(
    integer    M,
    integer    N,
    real const A[],
    integer    LDA,
    real       B[],
    integer    LDB
  )
  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  { LAPACK_F77NAME(slacpy)( "A", &M, &N, A, &LDA, B, &LDB ); return 0; }
  #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { LAPACK_F77NAME(slacpy)(
      const_cast<character*>("A"), &M, &N, A, &LDA, B, &LDB
    );
    return 0;
  }
  #elif defined(LAPACK_WRAPPER_USE_ATLAS)
  { for ( integer j = 0; j < N; ++j ) copy( M, A+j*LDA, 1, B+j*LDB, 1 );
    return 0; }
  #elif defined(LAPACK_WRAPPER_USE_MKL)
  { slacpy( "A", &M, &N, A, &LDA, B, &LDB ); return 0; }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
  { return CLAPACKNAME(slacpy)( const_cast<character*>("A"), &M, &N,
                                const_cast<real*>(A), &LDA, B, &LDB ); }
  #else
    #error "LapackWrapper undefined mapping!"
  #endif

  inline
  integer
  gecopy(
    integer          M,
    integer          N,
    doublereal const A[],
    integer          LDA,
    doublereal       B[],
    integer          LDB
  )
  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  { LAPACK_F77NAME(dlacpy)( "A", &M, &N, A, &LDA, B, &LDB ); return 0; }
  #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { LAPACK_F77NAME(dlacpy)(
      const_cast<character*>("A"), &M, &N, A, &LDA, B, &LDB
    );
    return 0;
  }
  #elif defined(LAPACK_WRAPPER_USE_ATLAS)
  { for ( integer j = 0; j < N; ++j ) copy( M, A+j*LDA, 1, B+j*LDB, 1 );
    return 0; }
  #elif defined(LAPACK_WRAPPER_USE_MKL)
  { dlacpy( "A", &M, &N, A, &LDA, B, &LDB ); return 0; }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
  { return CLAPACKNAME(dlacpy)( const_cast<character*>("A"), &M, &N,
                                const_cast<doublereal*>(A), &LDA, B, &LDB ); }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  /*
  //    __ _  ___ _______ _ __ ___
  //   / _` |/ _ \_  / _ \ '__/ _ \
  //  | (_| |  __// /  __/ | | (_) |
  //   \__, |\___/___\___|_|  \___/
  //   |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DLASET initializes an m-by-n matrix A to BETA on the diagonal and
   *  ALPHA on the offdiagonals.
   *
   *  Arguments
   *  =========
   *
   *  UPLO    (input) CHARACTER*1
   *          Specifies the part of the matrix A to be set.
   *          = 'U':      Upper triangular part is set; the strictly lower
   *                      triangular part of A is not changed.
   *          = 'L':      Lower triangular part is set; the strictly upper
   *                      triangular part of A is not changed.
   *          Otherwise:  All of the matrix A is set.
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix A.  M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A.  N >= 0.
   *
   *  ALPHA   (input) DOUBLE PRECISION
   *          The constant to which the offdiagonal elements are to be set.
   *
   *  BETA    (input) DOUBLE PRECISION
   *          The constant to which the diagonal elements are to be set.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On exit, the leading m-by-n submatrix of A is set as follows:
   *
   *          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
   *          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
   *          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
   *
   *          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
  \*/
  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {
    void
    LAPACK_F77NAME(slaset)(
      character const   UPLO[],
      integer   const * M,
      integer   const * N,
      real      const * ALPHA,
      real      const * BETA,
      real              A[],
      integer   const * LDA
    );

    void
    LAPACK_F77NAME(dlaset)(
      character  const   UPLO[],
      integer    const * M,
      integer    const * N,
      doublereal const * ALPHA,
      doublereal const * BETA,
      doublereal         A[],
      integer    const * LDA
    );
  }
  #endif

  inline
  void
  gezero(
    integer M,
    integer N,
    real    A[],
    integer LDA
  )
  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  { real zero = 0; LAPACK_F77NAME(slaset)( "A", &M, &N, &zero, &zero, A, &LDA ); }
  #elif defined(LAPACK_WRAPPER_USE_ATLAS) || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { for ( integer j = 0; j < N; ++j ) zero( M, A+j*LDA, 1 ); }
  #elif defined(LAPACK_WRAPPER_USE_MKL)
  { real zero = 0; slaset( "A", &M, &N, &zero, &zero, A, &LDA ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
  { real zero = 0; CLAPACKNAME(slaset)( const_cast<character*>("A"), &M, &N, &zero, &zero, A, &LDA ); }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  inline
  void
  gezero(
    integer    M,
    integer    N,
    doublereal A[],
    integer    LDA
  )
  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  { doublereal zero = 0; LAPACK_F77NAME(dlaset)( "A", &M, &N, &zero, &zero, A, &LDA ); }
  #elif defined(LAPACK_WRAPPER_USE_ATLAS) || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { for ( integer j = 0; j < N; ++j ) zero( M, A+j*LDA, 1 ); }
  #elif defined(LAPACK_WRAPPER_USE_MKL)
  { doublereal zero = 0; dlaset( "A", &M, &N, &zero, &zero, A, &LDA ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
  { doublereal zero = 0; CLAPACKNAME(dlaset)( const_cast<character*>("A"), &M, &N, &zero, &zero, A, &LDA ); }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  inline
  void
  gefill(
    integer M,
    integer N,
    real    A[],
    integer LDA,
    real    value
  )
  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  { LAPACK_F77NAME(slaset)( "A", &M, &N, &value, &value, A, &LDA ); }
  #elif defined(LAPACK_WRAPPER_USE_ATLAS) || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { for ( integer j = 0; j < N; ++j ) copy( M, &value, 0, A+j*LDA, 1 ); }
  #elif defined(LAPACK_WRAPPER_USE_MKL)
  { slaset( "A", &M, &N, &value, &value, A, &LDA ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
  { CLAPACKNAME(slaset)( const_cast<character*>("A"), &M, &N, &value, &value, A, &LDA ); }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  inline
  void
  gefill(
    integer    M,
    integer    N,
    doublereal A[],
    integer    LDA,
    doublereal value
  )
  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  { LAPACK_F77NAME(dlaset)( "A", &M, &N, &value, &value, A, &LDA ); }
  #elif defined(LAPACK_WRAPPER_USE_ATLAS) || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { for ( integer j = 0; j < N; ++j ) copy( M, &value, 0, A+j*LDA, 1 ); }
  #elif defined(LAPACK_WRAPPER_USE_MKL)
  { dlaset( "A", &M, &N, &value, &value, A, &LDA ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
  { CLAPACKNAME(dlaset)( const_cast<character*>("A"), &M, &N, &value, &value, A, &LDA ); }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  /*
  //             _     _
  //   __ _  ___(_) __| |
  //  / _` |/ _ \ |/ _` |
  // | (_| |  __/ | (_| |
  //  \__, |\___|_|\__,_|
  //  |___/
  */
  inline
  void
  geid(
    integer M,
    integer N,
    real    A[],
    integer LDA,
    real    diag = 1
  )
  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  { real zero = 0; LAPACK_F77NAME(slaset)( "A", &M, &N, &zero, &diag, A, &LDA ); }
  #elif defined(LAPACK_WRAPPER_USE_ATLAS) || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { gezero(M,N,A,LDA); copy( std::min(M,N), &diag, 0, A, LDA+1); }
  #elif defined(LAPACK_WRAPPER_USE_MKL)
  { real zero = 0; slaset( "A", &M, &N, &zero, &diag, A, &LDA ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
  { real zero = 0; CLAPACKNAME(slaset)( const_cast<character*>("A"), &M, &N, &zero, &diag, A, &LDA ); }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  inline
  void
  geid(
    integer    M,
    integer    N,
    doublereal A[],
    integer    LDA,
    doublereal diag = 1
  )
  #if defined(LAPACK_WRAPPER_USE_LAPACK)
  { doublereal zero = 0; LAPACK_F77NAME(dlaset)( "A", &M, &N, &zero, &diag, A, &LDA ); }
  #elif defined(LAPACK_WRAPPER_USE_ATLAS) || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { gezero(M,N,A,LDA); copy( std::min(M,N), &diag, 0, A, LDA+1); }
  #elif defined(LAPACK_WRAPPER_USE_MKL)
  { doublereal zero = 0; dlaset( "A", &M, &N, &zero, &diag, A, &LDA ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
  { doublereal zero = 0; CLAPACKNAME(dlaset)( const_cast<character*>("A"), &M, &N, &zero, &diag, A, &LDA ); }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  /*
  //                        _     _
  //    __ _  ___  __ _  __| | __| |
  //   / _` |/ _ \/ _` |/ _` |/ _` |
  //  | (_| |  __/ (_| | (_| | (_| |
  //   \__, |\___|\__,_|\__,_|\__,_|
  //   |___/
  */
  /*
   * C <- alpha * A + beta * B
   */

  inline
  void
  geadd(
    integer    M,
    integer    N,
    real       alpha,
    real const A[], integer LDA,
    real       beta,
    real const B[], integer LDB,
    real       C[], integer LDC
  )
  #if defined(LAPACK_WRAPPER_USE_ACCELERATE) && MAC_OS_X_VERSION_MIN_REQUIRED >= MAC_OS_X_VERSION_10_10
  {
    if ( M > 0 && N > 0 )
      appleblas_sgeadd( CblasColMajor, CblasNoTrans, CblasNoTrans,
                        M, N, alpha, A, LDA, beta, B, LDB, C, LDC );
  }
  #else
  {
    integer ierr = gecopy( M, N, B, LDB, C, LDC );
    LAPACK_WRAPPER_ASSERT( ierr == 0, "geadd, ierr = " << ierr );
    real const * Aj = A;
    real       * Cj = C;
    for ( integer j = 0; j < N; ++j, Aj += LDA, Cj += LDC ) {
      lapack_wrapper::scal( M, beta,  Cj, 1 );
      lapack_wrapper::axpy( M, alpha, Aj, 1, Cj, 1 );
    }
  }
  #endif

  inline
  void
  geadd(
    integer          M,
    integer          N,
    doublereal       alpha,
    doublereal const A[], integer LDA,
    doublereal       beta,
    doublereal const B[], integer LDB,
    doublereal       C[], integer LDC
  )
  #if defined(LAPACK_WRAPPER_USE_ACCELERATE) && MAC_OS_X_VERSION_MIN_REQUIRED >= MAC_OS_X_VERSION_10_10
  {
    if ( M > 0 && N > 0 )
      appleblas_dgeadd( CblasColMajor, CblasNoTrans, CblasNoTrans,
                        M, N, alpha, A, LDA, beta, B, LDB, C, LDC );
  }
  #else
  {
    integer ierr = gecopy( M, N, B, LDB, C, LDC );
    LAPACK_WRAPPER_ASSERT( ierr == 0, "geadd, ierr = " << ierr );
    doublereal const * Aj = A;
    doublereal       * Cj = C;
    for ( integer j = 0; j < N; ++j, Aj += LDA, Cj += LDC ) {
      lapack_wrapper::scal( M, beta,  Cj, 1 );
      lapack_wrapper::axpy( M, alpha, Aj, 1, Cj, 1 );
    }
  }
  #endif

  /*
  //   __  __       _        _        _   _
  //  |  \/  | __ _| |_ _ __(_)_  __ | \ | | ___  _ __ _ __ ___
  //  | |\/| |/ _` | __| '__| \ \/ / |  \| |/ _ \| '__| '_ ` _ \
  //  | |  | | (_| | |_| |  | |>  <  | |\  | (_) | |  | | | | | |
  //  |_|  |_|\__,_|\__|_|  |_/_/\_\ |_| \_|\___/|_|  |_| |_| |_|
  */
  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  // use standard Lapack routine
  extern "C" {

    // 'M' or 'm'            maxabs
    // '1', 'O' or 'o'       norm1
    // 'I' or 'i'            normI
    // 'F', 'f', 'E' or 'e'  normF

    doublereal
    LAPACK_F77NAME(dlange)(
      character  const   NORM[],
      integer    const * M,
      integer    const * N,
      doublereal const * A,
      integer    const * LDA,
      doublereal       * WORK
    );

    real
    LAPACK_F77NAME(slange)(
      character const   NORM[],
      integer   const * M,
      integer   const * N,
      real      const * A,
      integer   const * LDA,
      real            * WORK
    );

  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DLANGE  returns the value of the one norm,  or the Frobenius norm, or
   *  the  infinity norm,  or the  element of  largest absolute value  of a
   *  real matrix A.
   *
   *  Description
   *  ===========
   *
   *  DLANGE returns the value
   *
   *     DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
   *              (
   *              ( norm1(A),         NORM = '1', 'O' or 'o'
   *              (
   *              ( normI(A),         NORM = 'I' or 'i'
   *              (
   *              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
   *
   *  where  norm1  denotes the  one norm of a matrix (maximum column sum),
   *  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
   *  normF  denotes the  Frobenius norm of a matrix (square root of sum of
   *  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
   *
   *  Arguments
   *  =========
   *
   *  NORM    (input) CHARACTER*1
   *          Specifies the value to be returned in DLANGE as described
   *          above.
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix A.  M >= 0.  When M = 0,
   *          DLANGE is set to zero.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A.  N >= 0.  When N = 0,
   *          DLANGE is set to zero.
   *
   *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
   *          The m by n matrix A.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(M,1).
   *
   *  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
   *          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
   *          referenced.
  \*/

  inline
  real
  normInf(
    integer    N,
    integer    M,
    real const A[],
    integer    LDA
  ) {
    real res = 0;
    for ( integer i=0; i < N; ++i ) {
      real normcol = lapack_wrapper::asum(M,A+i,LDA);
      if ( normcol > res ) res = normcol;
    }
    return res;
  }

  inline
  doublereal
  normInf(
    integer          N,
    integer          M,
    doublereal const A[],
    integer          LDA
  ) {
    doublereal res = 0;
    for ( integer i=0; i < N; ++i ) {
      doublereal normcol = lapack_wrapper::asum(M,A+i,LDA);
      if ( normcol > res ) res = normcol;
    }
    return res;
  }

  ///////////////////

  inline
  real
  norm1(
    integer    N,
    integer    M,
    real const A[],
    integer    LDA
  )
  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  { return LAPACK_F77NAME(slange)( "1", &N, &M, A, &LDA, nullptr ); }
  #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { return LAPACK_F77NAME(slange)(
      const_cast<character*>("1"), &N, &M, A, &LDA, nullptr
    );
  }
  #elif defined(LAPACK_WRAPPER_USE_MKL)
  { return slange( "1", &N, &M, A, &LDA, nullptr ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
  { return real(CLAPACKNAME(slange)(
      const_cast<character*>("1"), &N, &M,
      const_cast<real*>(A), &LDA, nullptr
    ));
  }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  inline
  doublereal
  norm1(
    integer          N,
    integer          M,
    doublereal const A[],
    integer          LDA
  )
  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  { return LAPACK_F77NAME(dlange)( "1", &N, &M, A, &LDA, nullptr ); }
  #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { return LAPACK_F77NAME(dlange)(
      const_cast<character*>("1"), &N, &M, A, &LDA, nullptr
    );
  }
  #elif defined(LAPACK_WRAPPER_USE_MKL)
  { return dlange( "1", &N, &M, A, &LDA, nullptr ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
  { return CLAPACKNAME(dlange)(
      const_cast<character*>("1"), &N, &M,
      const_cast<doublereal*>(A), &LDA, nullptr
    );
  }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  ///////////////////

  inline
  return_precision
  normF(
    integer    N,
    integer    M,
    real const A[],
    integer    LDA
  ) {
  #if defined(LAPACK_WRAPPER_USE_ACCELERATE)
    return CLAPACKNAME(slange)(
      const_cast<character*>("F"), &N, &M,
      const_cast<real*>(A), &LDA, nullptr
    );
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    return return_precision(LAPACK_F77NAME(slange)(
      const_cast<character*>("F"), &N, &M, A, &LDA, nullptr
    ));
  #elif defined(LAPACK_WRAPPER_USE_MKL)
    return return_precision(slange( "F", &N, &M, A, &LDA, nullptr ));
  #else
  #error "LapackWrapper undefined mapping!"
  #endif
  }

  inline
  doublereal
  normF(
    integer          N,
    integer          M,
    doublereal const A[],
    integer          LDA
  ) {
  #if defined(LAPACK_WRAPPER_USE_ACCELERATE)
    return CLAPACKNAME(dlange)(
      const_cast<character*>("F"), &N, &M,
      const_cast<doublereal*>(A), &LDA, nullptr
    );
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    return LAPACK_F77NAME(dlange)(
      const_cast<character*>("F"), &N, &M, A, &LDA, nullptr
    );
  #elif defined(LAPACK_WRAPPER_USE_MKL)
    return dlange( "F", &N, &M, A, &LDA, nullptr );
  #else
  #error "LapackWrapper undefined mapping!"
  #endif
  }

  ///////////////////

  inline
  real
  maxabs(
    integer    N,
    integer    M,
    real const A[],
    integer    LDA
  ) {
  #if defined(LAPACK_WRAPPER_USE_ACCELERATE)
    return real(CLAPACKNAME(slange)(
      const_cast<character*>("M"), &N, &M,
      const_cast<real*>(A), &LDA, nullptr
    ));
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    return LAPACK_F77NAME(slange)(
      const_cast<character*>("M"), &N, &M, A, &LDA, nullptr
    );
  #elif defined(LAPACK_WRAPPER_USE_MKL)
    return slange( "M", &N, &M, A, &LDA, nullptr );
  #else
  #error "LapackWrapper undefined mapping!"
  #endif
  }

  inline
  doublereal
  maxabs(
    integer          N,
    integer          M,
    doublereal const A[],
    integer          LDA
  ) {
  #if defined(LAPACK_WRAPPER_USE_ACCELERATE)
    return CLAPACKNAME(dlange)(
      const_cast<character*>("M"), &N, &M,
      const_cast<doublereal*>(A), &LDA, nullptr
    );

  #elif defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    return LAPACK_F77NAME(dlange)(
      const_cast<character*>("M"), &N, &M, A, &LDA, nullptr
    );
  #elif defined(LAPACK_WRAPPER_USE_MKL)
    return dlange( "M", &N, &M, A, &LDA, nullptr );
  #else
  #error "LapackWrapper undefined mapping!"
  #endif
  }

  #ifndef LAPACK_WRAPPER_OS_WINDOWS

  /*
  //   _                _
  //  | | __ _ ___  ___| |
  //  | |/ _` / __|/ __| |
  //  | | (_| \__ \ (__| |
  //  |_|\__,_|___/\___|_|
  */
  /*\
   *
   *  Purpose
   *  =======
   *
   *  DLASCL multiplies the M by N real matrix A by the real scalar
   *  CTO/CFROM.  This is done without over/underflow as long as the final
   *  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
   *  A may be full, upper triangular, lower triangular, upper Hessenberg,
   *  or banded.
   *
   *  Arguments
   *  =========
   *
   *  TYPE    (input) CHARACTER*1
   *          TYPE indices the storage type of the input matrix.
   *          = 'G':  A is a full matrix.
   *          = 'L':  A is a lower triangular matrix.
   *          = 'U':  A is an upper triangular matrix.
   *          = 'H':  A is an upper Hessenberg matrix.
   *          = 'B':  A is a symmetric band matrix with lower bandwidth KL
   *                  and upper bandwidth KU and with the only the lower
   *                  half stored.
   *          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
   *                  and upper bandwidth KU and with the only the upper
   *                  half stored.
   *          = 'Z':  A is a band matrix with lower bandwidth KL and upper
   *                  bandwidth KU.
   *
   *  KL      (input) INTEGER
   *          The lower bandwidth of A.  Referenced only if TYPE = 'B',
   *          'Q' or 'Z'.
   *
   *  KU      (input) INTEGER
   *          The upper bandwidth of A.  Referenced only if TYPE = 'B',
   *          'Q' or 'Z'.
   *
   *  CFROM   (input) DOUBLE PRECISION
   *  CTO     (input) DOUBLE PRECISION
   *          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
   *          without over/underflow if the final result CTO*A(I,J)/CFROM
   *          can be represented without over/underflow.  CFROM must be
   *          nonzero.
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix A.  M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
   *          storage type.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  INFO    (output) INTEGER
   *          0  - successful exit
   *          <0 - if INFO = -i, the i-th argument had an illegal value.
   *
   *  =====================================================================
  \*/
  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  // use standard Lapack routine
  extern "C" {
    void
    LAPACK_F77NAME(dlascl)(
      character  const   TYPE[],
      integer    const * KL,
      integer    const * KU,
      doublereal const * FROM,
      doublereal const * TO,
      integer    const * M,
      integer    const * N,
      doublereal       * A,
      integer    const * LDA,
      integer          * info
    );

    void
    LAPACK_F77NAME(slascl)(
      character const   TYPE[],
      integer   const * KL,
      integer   const * KU,
      real      const * FROM,
      real      const * TO,
      integer   const * M,
      integer   const * N,
      real            * A,
      integer   const * LDA,
      integer         * info
    );
  }
  #endif

  inline
  integer
  lascl(
    MatrixType const & TYPE,
    integer            KL,
    integer            KU,
    real               FROM,
    real               TO,
    integer            M,
    integer            N,
    real               A[],
    integer            LDA
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(slascl)(
      const_cast<character*>(mtype_blas[TYPE]),
      &KL, &KU, &FROM, &TO, &M, &N, A, &LDA, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_LAPACK) || \
          defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(slascl)(
      mtype_blas[TYPE],
      &KL, &KU, &FROM, &TO, &M, &N, A, &LDA, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
    LAPACK_F77NAME(slascl)(
      const_cast<character*>(mtype_blas[TYPE]),
      &KL, &KU, &FROM, &TO, &M, &N, A, &LDA, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    slascl( mtype_blas[TYPE], &KL, &KU, &FROM, &TO, &M, &N, A, &LDA, &info);
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  inline
  integer
  lascl(
    MatrixType const & TYPE,
    integer            KL,
    integer            KU,
    doublereal         FROM,
    doublereal         TO,
    integer            M,
    integer            N,
    doublereal         A[],
    integer            LDA
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dlascl)(
      const_cast<character*>(mtype_blas[TYPE]),
      &KL, &KU, &FROM, &TO, &M, &N, A, &LDA, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_LAPACK) || \
          defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dlascl)(
      mtype_blas[TYPE], &KL, &KU, &FROM, &TO, &M, &N, A, &LDA, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_LAPACK)  || \
          defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
          defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dlascl)(
      const_cast<character*>(mtype_blas[TYPE]),
      &KL, &KU, &FROM, &TO, &M, &N, A, &LDA, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dlascl( mtype_blas[TYPE], &KL, &KU, &FROM, &TO, &M, &N, A, &LDA, &info);
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  #endif

  /*
  //   _____ _                            _
  //  | ____(_) __ _  ___ _ ____   ____ _| |_   _  ___  ___
  //  |  _| | |/ _` |/ _ \ '_ \ \ / / _` | | | | |/ _ \/ __|
  //  | |___| | (_| |  __/ | | \ V / (_| | | |_| |  __/\__ \
  //  |_____|_|\__, |\___|_| |_|\_/ \__,_|_|\__,_|\___||___/
  //           |___/
  */
  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  // use standard Lapack routine
  extern "C" {

    // 'M' or 'm'            maxabs
    // '1', 'O' or 'o'       norm1
    // 'I' or 'i'            normI
    // 'F', 'f', 'E' or 'e'  normF

    // DGEEV - compute for an N-by-N real nonsymmetric matrix A,
    //         the eigenvalues and, optionally, the left and/or right
    //          eigenvectors

    void
    LAPACK_F77NAME(dgeev)(
      character  const   jobvl[],
      character  const   jobvr[],
      integer    const * N,
      doublereal       * A,
      integer    const * LDA,
      doublereal       * wr,
      doublereal       * wi,
      doublereal       * vl,
      integer    const * Lvl,
      doublereal       * vr,
      integer    const * Lvr,
      doublereal       * WORK,
      integer    const * Lwork,
      integer          * info
    );

    void
    LAPACK_F77NAME(sgeev)(
      character const   jobvl[],
      character const   jobvr[],
      integer   const * N,
      real            * A,
      integer   const * LDA,
      real            * wr,
      real            * wi,
      real            * vl,
      integer   const * Lvl,
      real            * vr,
      integer   const * Lvr,
      real            * WORK,
      integer   const * Lwork,
      integer         * info
    );
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DGEEV computes for an N-by-N real nonsymmetric matrix A, the
   *  eigenvalues and, optionally, the left and/or right eigenvectors.
   *
   *  The right eigenvector v(j) of A satisfies
   *                   A * v(j) = lambda(j) * v(j)
   *  where lambda(j) is its eigenvalue.
   *  The left eigenvector u(j) of A satisfies
   *                u(j)**T * A = lambda(j) * u(j)**T
   *  where u(j)**T denotes the transpose of u(j).
   *
   *  The computed eigenvectors are normalized to have Euclidean norm
   *  equal to 1 and largest component real.
   *
   *  Arguments
   *  =========
   *
   *  JOBVL   (input) CHARACTER*1
   *          = 'N': left eigenvectors of A are not computed;
   *          = 'V': left eigenvectors of A are computed.
   *
   *  JOBVR   (input) CHARACTER*1
   *          = 'N': right eigenvectors of A are not computed;
   *          = 'V': right eigenvectors of A are computed.
   *
   *  N       (input) INTEGER
   *          The order of the matrix A. N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the N-by-N matrix A.
   *          On exit, A has been overwritten.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,N).
   *
   *  WR      (output) DOUBLE PRECISION array, dimension (N)
   *  WI      (output) DOUBLE PRECISION array, dimension (N)
   *          WR and WI contain the real and imaginary parts,
   *          respectively, of the computed eigenvalues.  Complex
   *          conjugate pairs of eigenvalues appear consecutively
   *          with the eigenvalue having the positive imaginary part
   *          first.
   *
   *  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
   *          If JOBVL = 'V', the left eigenvectors u(j) are stored one
   *          after another in the columns of VL, in the same order
   *          as their eigenvalues.
   *          If JOBVL = 'N', VL is not referenced.
   *          If the j-th eigenvalue is real, then u(j) = VL(:,j),
   *          the j-th column of VL.
   *          If the j-th and (j+1)-st eigenvalues form a complex
   *          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
   *          u(j+1) = VL(:,j) - i*VL(:,j+1).
   *
   *  LDVL    (input) INTEGER
   *          The leading dimension of the array VL.  LDVL >= 1; if
   *          JOBVL = 'V', LDVL >= N.
   *
   *  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
   *          If JOBVR = 'V', the right eigenvectors v(j) are stored one
   *          after another in the columns of VR, in the same order
   *          as their eigenvalues.
   *          If JOBVR = 'N', VR is not referenced.
   *          If the j-th eigenvalue is real, then v(j) = VR(:,j),
   *          the j-th column of VR.
   *          If the j-th and (j+1)-st eigenvalues form a complex
   *          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
   *          v(j+1) = VR(:,j) - i*VR(:,j+1).
   *
   *  LDVR    (input) INTEGER
   *          The leading dimension of the array VR.  LDVR >= 1; if
   *          JOBVR = 'V', LDVR >= N.
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK.  LWORK >= max(1,3*N), and
   *          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good
   *          performance, LWORK must generally be larger.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value.
   *          > 0:  if INFO = i, the QR algorithm failed to compute all the
   *                eigenvalues, and no eigenvectors have been computed;
   *                elements i+1:N of WR and WI contain eigenvalues which
   *                have converged.
  \*/

  inline
  integer
  geev(
    bool    jobvl, // false = do not compute the left generalized eigenvectors
    bool    jobvr, // false = do not compute the right generalized eigenvectors
    integer n,     // The order of the matrices A, B, VL, and VR.  N >= 0.
    real    A[],
    integer lda,
    real    wr[],
    real    wi[],
    real    vl[],
    integer ldvl,
    real    vr[],
    integer ldvr,
    real    work[],
    integer lwork
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sgeev)(
      (jobvl?"V":"N"),
      (jobvr?"V":"N"),
      &n, A, &lda, wr, wi,
      vl, &ldvl,
      vr, &ldvr, work, &lwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
    LAPACK_F77NAME(sgeev)(
      const_cast<character*>(jobvl?"V":"N"),
      const_cast<character*>(jobvr?"V":"N"),
      &n, A, &lda, wr, wi,
      vl, &ldvl,
      vr, &ldvr, work, &lwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgeev(
      (jobvl?"V":"N"),
      (jobvr?"V":"N"),
      &n, A, &lda, wr, wi,
      vl, &ldvl, vr, &ldvr, work, &lwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgeev)(
      const_cast<character*>(jobvl?"V":"N"),
      const_cast<character*>(jobvr?"V":"N"),
      &n,
      const_cast<real*>(A), &lda, wr, wi,
      vl, &ldvl,
      vr, &ldvr,
      work, &lwork, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  inline
  integer
  geev(
    bool       jobvl, // false = do not compute the left generalized eigenvectors
    bool       jobvr, // false = do not compute the right generalized eigenvectors
    integer    n,     // The order of the matrices A, B, VL, and VR.  N >= 0.
    doublereal A[],
    integer    lda,
    doublereal wr[],
    doublereal wi[],
    doublereal vl[],
    integer    ldvl,
    doublereal vr[],
    integer    ldvr,
    doublereal work[],
    integer    lwork
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dgeev)(
      (jobvl?"V":"N"),
      (jobvr?"V":"N"),
      &n, A, &lda, wr, wi,
      vl, &ldvl,
      vr, &ldvr, work, &lwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
    LAPACK_F77NAME(dgeev)(
      const_cast<character*>(jobvl?"V":"N"),
      const_cast<character*>(jobvr?"V":"N"),
      &n, A, &lda, wr, wi,
      vl, &ldvl,
      vr, &ldvr, work, &lwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgeev(
      (jobvl?"V":"N"),
      (jobvr?"V":"N"),
      &n, A, &lda, wr, wi,
      vl, &ldvl, vr, &ldvr, work, &lwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgeev)(
      const_cast<character*>(jobvl?"V":"N"),
      const_cast<character*>(jobvr?"V":"N"),
      &n,
      const_cast<doublereal*>(A), &lda, wr, wi,
      vl, &ldvl,
      vr, &ldvr,
      work, &lwork, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  //////////////////////////////////////////////////////////////////////////////

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  // use standard Lapack routine
  extern "C" {

    // 'M' or 'm'            maxabs
    // '1', 'O' or 'o'       norm1
    // 'I' or 'i'            normI
    // 'F', 'f', 'E' or 'e'  normF

    //  DGGEV computes for a pair of N-by-N real nonsymmetric matrices (A,B)
    //  the generalized eigenvalues, and optionally, the left and/or right
    //  generalized eigenvectors.
    //
    //  A generalized eigenvalue for a pair of matrices (A,B) is a scalar
    //  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
    //  singular. It is usually represented as the pair (alpha,beta), as
    //  there is a reasonable interpretation for beta=0, and even for both
    //  being zero.

    void
    LAPACK_F77NAME(dggev)(
      character  const   jobvl[],
      character  const   jobvr[],
      integer    const * N,
      doublereal       * A,
      integer    const * LDA,
      doublereal       * B,
      integer    const * LDB,
      doublereal       * alphar,
      doublereal       * alphai,
      doublereal       * beta,
      doublereal       * vl,
      integer    const * Lvl,
      doublereal       * vr,
      integer    const * Lvr,
      doublereal       * WORK,
      integer    const * Lwork,
      integer          * info
    );

    void
    LAPACK_F77NAME(sggev)(
      character const   jobvl[],
      character const   jobvr[],
      integer   const * N,
      real            * A,
      integer   const * LDA,
      real            * B,
      integer   const * LDB,
      real            * alphar,
      real            * alphai,
      real            * beta,
      real            * vl,
      integer   const * Lvl,
      real            * vr,
      integer   const * Lvr,
      real            * WORK,
      integer   const * Lwork,
      integer         * info
    );
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DGGEV computes for a pair of N-by-N real nonsymmetric matrices (A,B)
   *  the generalized eigenvalues, and optionally, the left and/or right
   *  generalized eigenvectors.
   *
   *  A generalized eigenvalue for a pair of matrices (A,B) is a scalar
   *  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
   *  singular. It is usually represented as the pair (alpha,beta), as
   *  there is a reasonable interpretation for beta=0, and even for both
   *  being zero.
   *
   *  The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
   *  of (A,B) satisfies
   *
   *                   A * v(j) = lambda(j) * B * v(j).
   *
   *  The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
   *  of (A,B) satisfies
   *
   *                   u(j)**H * A  = lambda(j) * u(j)**H * B .
   *
   *  where u(j)**H is the conjugate-transpose of u(j).
   *
   *
   *  Arguments
   *  =========
   *
   *  JOBVL   (input) CHARACTER*1
   *          = 'N':  do not compute the left generalized eigenvectors;
   *          = 'V':  compute the left generalized eigenvectors.
   *
   *  JOBVR   (input) CHARACTER*1
   *          = 'N':  do not compute the right generalized eigenvectors;
   *          = 'V':  compute the right generalized eigenvectors.
   *
   *  N       (input) INTEGER
   *          The order of the matrices A, B, VL, and VR.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
   *          On entry, the matrix A in the pair (A,B).
   *          On exit, A has been overwritten.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of A.  LDA >= max(1,N).
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
   *          On entry, the matrix B in the pair (A,B).
   *          On exit, B has been overwritten.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of B.  LDB >= max(1,N).
   *
   *  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
   *  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
   *  BETA    (output) DOUBLE PRECISION array, dimension (N)
   *          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
   *          be the generalized eigenvalues.  If ALPHAI(j) is zero, then
   *          the j-th eigenvalue is real; if positive, then the j-th and
   *          (j+1)-st eigenvalues are a complex conjugate pair, with
   *          ALPHAI(j+1) negative.
   *
   *          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
   *          may easily over- or underflow, and BETA(j) may even be zero.
   *          Thus, the user should avoid naively computing the ratio
   *          alpha/beta.  However, ALPHAR and ALPHAI will be always less
   *          than and usually comparable with norm(A) in magnitude, and
   *          BETA always less than and usually comparable with norm(B).
   *
   *  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
   *          If JOBVL = 'V', the left eigenvectors u(j) are stored one
   *          after another in the columns of VL, in the same order as
   *          their eigenvalues. If the j-th eigenvalue is real, then
   *          u(j) = VL(:,j), the j-th column of VL. If the j-th and
   *          (j+1)-th eigenvalues form a complex conjugate pair, then
   *          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
   *          Each eigenvector is scaled so the largest component has
   *          abs(real part)+abs(imag. part)=1.
   *          Not referenced if JOBVL = 'N'.
   *
   *  LDVL    (input) INTEGER
   *          The leading dimension of the matrix VL. LDVL >= 1, and
   *          if JOBVL = 'V', LDVL >= N.
   *
   *  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
   *          If JOBVR = 'V', the right eigenvectors v(j) are stored one
   *          after another in the columns of VR, in the same order as
   *          their eigenvalues. If the j-th eigenvalue is real, then
   *          v(j) = VR(:,j), the j-th column of VR. If the j-th and
   *          (j+1)-th eigenvalues form a complex conjugate pair, then
   *          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
   *          Each eigenvector is scaled so the largest component has
   *          abs(real part)+abs(imag. part)=1.
   *          Not referenced if JOBVR = 'N'.
   *
   *  LDVR    (input) INTEGER
   *          The leading dimension of the matrix VR. LDVR >= 1, and
   *          if JOBVR = 'V', LDVR >= N.
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK.  LWORK >= max(1,8*N).
   *          For good performance, LWORK must generally be larger.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value.
   *          = 1,...,N:
   *                The QZ iteration failed.  No eigenvectors have been
   *                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
   *                should be correct for j=INFO+1,...,N.
   *          > N:  =N+1: other than QZ iteration failed in DHGEQZ.
   *                =N+2: error return from DTGEVC.
   *
  \*/

  inline
  integer
  ggev(
    bool    jobvl, // false = do not compute the left generalized eigenvectors
    bool    jobvr, // false = do not compute the right generalized eigenvectors
    integer n,     // The order of the matrices A, B, VL, and VR.  N >= 0.
    real    A[],
    integer lda,
    real    B[],
    integer ldb,
    real    alphar[],
    real    alphai[],
    real    beta[],
    real    vl[],
    integer ldvl,
    real    vr[],
    integer ldvr,
    real    work[],
    integer lwork
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sggev)(
      (jobvl?"V":"N"),
      (jobvr?"V":"N"),
      &n, A, &lda, B, &ldb, alphar, alphai, beta,
      vl, &ldvl, vr, &ldvr, work, &lwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
    LAPACK_F77NAME(sggev)(
      const_cast<character*>(jobvl?"V":"N"),
      const_cast<character*>(jobvr?"V":"N"),
      &n, A, &lda, B, &ldb, alphar, alphai, beta,
      vl, &ldvl, vr, &ldvr, work, &lwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sggev(
      (jobvl?"V":"N"),
      (jobvr?"V":"N"),
      &n, A, &lda, B, &ldb, alphar, alphai, beta,
      vl, &ldvl, vr, &ldvr, work, &lwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sggev)(
      const_cast<character*>(jobvl?"V":"N"),
      const_cast<character*>(jobvr?"V":"N"),
      &n,
      const_cast<real*>(A), &lda,
      const_cast<real*>(B), &ldb,
      alphar, alphai, beta,
      vl, &ldvl, vr, &ldvr, work, &lwork, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  inline
  integer
  ggev(
    bool       jobvl, // false = do not compute the left generalized eigenvectors
    bool       jobvr, // false = do not compute the right generalized eigenvectors
    integer    n,     // The order of the matrices A, B, VL, and VR.  N >= 0.
    doublereal A[],
    integer	   lda,
    doublereal B[],
    integer    ldb,
    doublereal alphar[],
    doublereal alphai[],
    doublereal beta[],
    doublereal vl[],
    integer    ldvl,
    doublereal vr[],
    integer    ldvr,
    doublereal work[],
    integer    lwork
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dggev)(
      (jobvl?"V":"N"),
      (jobvr?"V":"N"),
      &n, A, &lda, B, &ldb, alphar, alphai, beta,
      vl, &ldvl, vr, &ldvr, work, &lwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
    LAPACK_F77NAME(dggev)(
      const_cast<character*>(jobvl?"V":"N"),
      const_cast<character*>(jobvr?"V":"N"),
      &n, A, &lda, B, &ldb, alphar, alphai, beta,
      vl, &ldvl, vr, &ldvr, work, &lwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dggev(
      (jobvl?"V":"N"),
      (jobvr?"V":"N"),
      &n, A, &lda, B, &ldb, alphar, alphai, beta,
      vl, &ldvl, vr, &ldvr, work, &lwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dggev)(
      const_cast<character*>(jobvl?"V":"N"),
      const_cast<character*>(jobvr?"V":"N"),
      &n,
      const_cast<doublereal*>(A), &lda,
      const_cast<doublereal*>(B), &ldb,
      alphar, alphai, beta,
      vl, &ldvl, vr, &ldvr, work, &lwork, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  /*\
   *  Purpose
   *  =======
   *
   *  DGGEVX computes for a pair of N-by-N real nonsymmetric matrices (A,B)
   *  the generalized eigenvalues, and optionally, the left and/or right
   *  generalized eigenvectors.
   *
   *  Optionally also, it computes a balancing transformation to improve
   *  the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
   *  LSCALE, RSCALE, ABNRM, and BBNRM), reciprocal condition numbers for
   *  the eigenvalues (RCONDE), and reciprocal condition numbers for the
   *  right eigenvectors (RCONDV).
   *
   *  A generalized eigenvalue for a pair of matrices (A,B) is a scalar
   *  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
   *  singular. It is usually represented as the pair (alpha,beta), as
   *  there is a reasonable interpretation for beta=0, and even for both
   *  being zero.
   *
   *  The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
   *  of (A,B) satisfies
   *
   *                   A * v(j) = lambda(j) * B * v(j) .
   *
   *  The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
   *  of (A,B) satisfies
   *
   *                   u(j)**H * A  = lambda(j) * u(j)**H * B.
   *
   *  where u(j)**H is the conjugate-transpose of u(j).
   *
   *
   *  Arguments
   *  =========
   *
   *  BALANC  (input) CHARACTER*1
   *          Specifies the balance option to be performed.
   *          = 'N':  do not diagonally scale or permute;
   *          = 'P':  permute only;
   *          = 'S':  scale only;
   *          = 'B':  both permute and scale.
   *          Computed reciprocal condition numbers will be for the
   *          matrices after permuting and/or balancing. Permuting does
   *          not change condition numbers (in exact arithmetic), but
   *          balancing does.
   *
   *  JOBVL   (input) CHARACTER*1
   *          = 'N':  do not compute the left generalized eigenvectors;
   *          = 'V':  compute the left generalized eigenvectors.
   *
   *  JOBVR   (input) CHARACTER*1
   *          = 'N':  do not compute the right generalized eigenvectors;
   *          = 'V':  compute the right generalized eigenvectors.
   *
   *  SENSE   (input) CHARACTER*1
   *          Determines which reciprocal condition numbers are computed.
   *          = 'N': none are computed;
   *          = 'E': computed for eigenvalues only;
   *          = 'V': computed for eigenvectors only;
   *          = 'B': computed for eigenvalues and eigenvectors.
   *
   *  N       (input) INTEGER
   *          The order of the matrices A, B, VL, and VR.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
   *          On entry, the matrix A in the pair (A,B).
   *          On exit, A has been overwritten. If JOBVL='V' or JOBVR='V'
   *          or both, then A contains the first part of the real Schur
   *          form of the "balanced" versions of the input A and B.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of A.  LDA >= max(1,N).
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
   *          On entry, the matrix B in the pair (A,B).
   *          On exit, B has been overwritten. If JOBVL='V' or JOBVR='V'
   *          or both, then B contains the second part of the real Schur
   *          form of the "balanced" versions of the input A and B.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of B.  LDB >= max(1,N).
   *
   *  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
   *  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
   *  BETA    (output) DOUBLE PRECISION array, dimension (N)
   *          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
   *          be the generalized eigenvalues.  If ALPHAI(j) is zero, then
   *          the j-th eigenvalue is real; if positive, then the j-th and
   *          (j+1)-st eigenvalues are a complex conjugate pair, with
   *          ALPHAI(j+1) negative.
   *
   *          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
   *          may easily over- or underflow, and BETA(j) may even be zero.
   *          Thus, the user should avoid naively computing the ratio
   *          ALPHA/BETA. However, ALPHAR and ALPHAI will be always less
   *          than and usually comparable with norm(A) in magnitude, and
   *          BETA always less than and usually comparable with norm(B).
   *
   *  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
   *          If JOBVL = 'V', the left eigenvectors u(j) are stored one
   *          after another in the columns of VL, in the same order as
   *          their eigenvalues. If the j-th eigenvalue is real, then
   *          u(j) = VL(:,j), the j-th column of VL. If the j-th and
   *          (j+1)-th eigenvalues form a complex conjugate pair, then
   *          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
   *          Each eigenvector will be scaled so the largest component have
   *          abs(real part) + abs(imag. part) = 1.
   *          Not referenced if JOBVL = 'N'.
   *
   *  LDVL    (input) INTEGER
   *          The leading dimension of the matrix VL. LDVL >= 1, and
   *          if JOBVL = 'V', LDVL >= N.
   *
   *  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
   *          If JOBVR = 'V', the right eigenvectors v(j) are stored one
   *          after another in the columns of VR, in the same order as
   *          their eigenvalues. If the j-th eigenvalue is real, then
   *          v(j) = VR(:,j), the j-th column of VR. If the j-th and
   *          (j+1)-th eigenvalues form a complex conjugate pair, then
   *          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
   *          Each eigenvector will be scaled so the largest component have
   *          abs(real part) + abs(imag. part) = 1.
   *          Not referenced if JOBVR = 'N'.
   *
   *  LDVR    (input) INTEGER
   *          The leading dimension of the matrix VR. LDVR >= 1, and
   *          if JOBVR = 'V', LDVR >= N.
   *
   *  ILO     (output) INTEGER
   *  IHI     (output) INTEGER
   *          ILO and IHI are integer values such that on exit
   *          A(i,j) = 0 and B(i,j) = 0 if i > j and
   *          j = 1,...,ILO-1 or i = IHI+1,...,N.
   *          If BALANC = 'N' or 'S', ILO = 1 and IHI = N.
   *
   *  LSCALE  (output) DOUBLE PRECISION array, dimension (N)
   *          Details of the permutations and scaling factors applied
   *          to the left side of A and B.  If PL(j) is the index of the
   *          row interchanged with row j, and DL(j) is the scaling
   *          factor applied to row j, then
   *            LSCALE(j) = PL(j)  for j = 1,...,ILO-1
   *                      = DL(j)  for j = ILO,...,IHI
   *                      = PL(j)  for j = IHI+1,...,N.
   *          The order in which the interchanges are made is N to IHI+1,
   *          then 1 to ILO-1.
   *
   *  RSCALE  (output) DOUBLE PRECISION array, dimension (N)
   *          Details of the permutations and scaling factors applied
   *          to the right side of A and B.  If PR(j) is the index of the
   *          column interchanged with column j, and DR(j) is the scaling
   *          factor applied to column j, then
   *            RSCALE(j) = PR(j)  for j = 1,...,ILO-1
   *                      = DR(j)  for j = ILO,...,IHI
   *                      = PR(j)  for j = IHI+1,...,N
   *          The order in which the interchanges are made is N to IHI+1,
   *          then 1 to ILO-1.
   *
   *  ABNRM   (output) DOUBLE PRECISION
   *          The one-norm of the balanced matrix A.
   *
   *  BBNRM   (output) DOUBLE PRECISION
   *          The one-norm of the balanced matrix B.
   *
   *  RCONDE  (output) DOUBLE PRECISION array, dimension (N)
   *          If SENSE = 'E' or 'B', the reciprocal condition numbers of
   *          the eigenvalues, stored in consecutive elements of the array.
   *          For a complex conjugate pair of eigenvalues two consecutive
   *          elements of RCONDE are set to the same value. Thus RCONDE(j),
   *          RCONDV(j), and the j-th columns of VL and VR all correspond
   *          to the j-th eigenpair.
   *          If SENSE = 'N or 'V', RCONDE is not referenced.
   *
   *  RCONDV  (output) DOUBLE PRECISION array, dimension (N)
   *          If SENSE = 'V' or 'B', the estimated reciprocal condition
   *          numbers of the eigenvectors, stored in consecutive elements
   *          of the array. For a complex eigenvector two consecutive
   *          elements of RCONDV are set to the same value. If the
   *          eigenvalues cannot be reordered to compute RCONDV(j),
   *          RCONDV(j) is set to 0; this can only occur when the true
   *          value would be very small anyway.
   *          If SENSE = 'N' or 'E', RCONDV is not referenced.
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK. LWORK >= max(1,2*N).
   *          If BALANC = 'S' or 'B', or JOBVL = 'V', or JOBVR = 'V',
   *          LWORK >= max(1,6*N).
   *          If SENSE = 'E' or 'B', LWORK >= max(1,10*N).
   *          If SENSE = 'V' or 'B', LWORK >= 2*N*N+8*N+16.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  IWORK   (workspace) INTEGER array, dimension (N+6)
   *          If SENSE = 'E', IWORK is not referenced.
   *
   *  BWORK   (workspace) LOGICAL array, dimension (N)
   *          If SENSE = 'N', BWORK is not referenced.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value.
   *          = 1,...,N:
   *                The QZ iteration failed.  No eigenvectors have been
   *                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
   *                should be correct for j=INFO+1,...,N.
   *          > N:  =N+1: other than QZ iteration failed in DHGEQZ.
   *                =N+2: error return from DTGEVC.
   *
   *  Further Details
   *  ===============
   *
   *  Balancing a matrix pair (A,B) includes, first, permuting rows and
   *  columns to isolate eigenvalues, second, applying diagonal similarity
   *  transformation to the rows and columns to make the rows and columns
   *  as close in norm as possible. The computed reciprocal condition
   *  numbers correspond to the balanced matrix. Permuting rows and columns
   *  will not change the condition numbers (in exact arithmetic) but
   *  diagonal scaling will.  For further explanation of balancing, see
   *  section 4.11.1.2 of LAPACK Users' Guide.
   *
   *  An approximate error bound on the chordal distance between the i-th
   *  computed generalized eigenvalue w and the corresponding exact
   *  eigenvalue lambda is
   *
   *       chord(w, lambda) <= EPS * norm(ABNRM, BBNRM) / RCONDE(I)
   *
   *  An approximate error bound for the angle between the i-th computed
   *  eigenvector VL(i) or VR(i) is given by
   *
   *       EPS * norm(ABNRM, BBNRM) / DIF(i).
   *
   *  For further explanation of the reciprocal condition numbers RCONDE
   *  and RCONDV, see section 4.11 of LAPACK User's Guide.
   *
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  // use standard Lapack routine
  extern "C" {

    // 'M' or 'm'            maxabs
    // '1', 'O' or 'o'       norm1
    // 'I' or 'i'            normI
    // 'F', 'f', 'E' or 'e'  normF

    void
    LAPACK_F77NAME(dggevx)(
      character  const   balanc[],
      character  const   jobvl[],
      character  const   jobvr[],
      character  const   sense[],
      integer    const * N,
      doublereal       * A,
      integer    const * LDA,
      doublereal       * B,
      integer    const * LDB,
      doublereal       * alphar,
      doublereal       * alphai,
      doublereal       * beta,
      doublereal       * vl,
      integer    const * Lvl,
      doublereal       * vr,
      integer    const * Lvr,
      integer          * ILO,
      integer          * IHI,
      doublereal       * LSCALE,
      doublereal       * RSCALE,
      doublereal       * ABNRM,
      doublereal       * BBNRM,
      doublereal       * RCONDE,
      doublereal       * RCONDV,
      doublereal       * WORK,
      integer    const * Lwork,
      integer          * Iwork,
      integer          * Bwork,
      integer          * info
    );

    void
    LAPACK_F77NAME(sggevx)(
      character  const   balanc[],
      character  const   jobvl[],
      character  const   jobvr[],
      character  const   sense[],
      integer    const  * N,
      real              * A,
      integer    const  * LDA,
      real              * B,
      integer    const  * LDB,
      real              * alphar,
      real              * alphai,
      real              * beta,
      real              * vl,
      integer    const  * Lvl,
      real              * vr,
      integer    const  * Lvr,
      integer           * ILO,
      integer           * IHI,
      real              * LSCALE,
      real              * RSCALE,
      real              * ABNRM,
      real              * BBNRM,
      real              * RCONDE,
      real              * RCONDV,
      real              * WORK,
      integer     const * Lwork,
      integer           * Iwork,
      integer           * Bwork,
      integer           * info
    );
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DGGEVX computes for a pair of N-by-N real nonsymmetric matrices (A,B)
   *  the generalized eigenvalues, and optionally, the left and/or right
   *  generalized eigenvectors.
   *
   *  Optionally also, it computes a balancing transformation to improve
   *  the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
   *  LSCALE, RSCALE, ABNRM, and BBNRM), reciprocal condition numbers for
   *  the eigenvalues (RCONDE), and reciprocal condition numbers for the
   *  right eigenvectors (RCONDV).
   *
   *  A generalized eigenvalue for a pair of matrices (A,B) is a scalar
   *  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
   *  singular. It is usually represented as the pair (alpha,beta), as
   *  there is a reasonable interpretation for beta=0, and even for both
   *  being zero.
   *
   *  The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
   *  of (A,B) satisfies
   *
   *                   A * v(j) = lambda(j) * B * v(j) .
   *
   *  The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
   *  of (A,B) satisfies
   *
   *                   u(j)**H * A  = lambda(j) * u(j)**H * B.
   *
   *  where u(j)**H is the conjugate-transpose of u(j).
   *
   *
   *  Arguments
   *  =========
   *
   *  BALANC  (input) CHARACTER*1
   *          Specifies the balance option to be performed.
   *          = 'N':  do not diagonally scale or permute;
   *          = 'P':  permute only;
   *          = 'S':  scale only;
   *          = 'B':  both permute and scale.
   *          Computed reciprocal condition numbers will be for the
   *          matrices after permuting and/or balancing. Permuting does
   *          not change condition numbers (in exact arithmetic), but
   *          balancing does.
   *
   *  JOBVL   (input) CHARACTER*1
   *          = 'N':  do not compute the left generalized eigenvectors;
   *          = 'V':  compute the left generalized eigenvectors.
   *
   *  JOBVR   (input) CHARACTER*1
   *          = 'N':  do not compute the right generalized eigenvectors;
   *          = 'V':  compute the right generalized eigenvectors.
   *
   *  SENSE   (input) CHARACTER*1
   *          Determines which reciprocal condition numbers are computed.
   *          = 'N': none are computed;
   *          = 'E': computed for eigenvalues only;
   *          = 'V': computed for eigenvectors only;
   *          = 'B': computed for eigenvalues and eigenvectors.
   *
   *  N       (input) INTEGER
   *          The order of the matrices A, B, VL, and VR.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
   *          On entry, the matrix A in the pair (A,B).
   *          On exit, A has been overwritten. If JOBVL='V' or JOBVR='V'
   *          or both, then A contains the first part of the real Schur
   *          form of the "balanced" versions of the input A and B.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of A.  LDA >= max(1,N).
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
   *          On entry, the matrix B in the pair (A,B).
   *          On exit, B has been overwritten. If JOBVL='V' or JOBVR='V'
   *          or both, then B contains the second part of the real Schur
   *          form of the "balanced" versions of the input A and B.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of B.  LDB >= max(1,N).
   *
   *  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
   *  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
   *  BETA    (output) DOUBLE PRECISION array, dimension (N)
   *          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
   *          be the generalized eigenvalues.  If ALPHAI(j) is zero, then
   *          the j-th eigenvalue is real; if positive, then the j-th and
   *          (j+1)-st eigenvalues are a complex conjugate pair, with
   *          ALPHAI(j+1) negative.
   *
   *          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
   *          may easily over- or underflow, and BETA(j) may even be zero.
   *          Thus, the user should avoid naively computing the ratio
   *          ALPHA/BETA. However, ALPHAR and ALPHAI will be always less
   *          than and usually comparable with norm(A) in magnitude, and
   *          BETA always less than and usually comparable with norm(B).
   *
   *  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
   *          If JOBVL = 'V', the left eigenvectors u(j) are stored one
   *          after another in the columns of VL, in the same order as
   *          their eigenvalues. If the j-th eigenvalue is real, then
   *          u(j) = VL(:,j), the j-th column of VL. If the j-th and
   *          (j+1)-th eigenvalues form a complex conjugate pair, then
   *          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
   *          Each eigenvector will be scaled so the largest component have
   *          abs(real part) + abs(imag. part) = 1.
   *          Not referenced if JOBVL = 'N'.
   *
   *  LDVL    (input) INTEGER
   *          The leading dimension of the matrix VL. LDVL >= 1, and
   *          if JOBVL = 'V', LDVL >= N.
   *
   *  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
   *          If JOBVR = 'V', the right eigenvectors v(j) are stored one
   *          after another in the columns of VR, in the same order as
   *          their eigenvalues. If the j-th eigenvalue is real, then
   *          v(j) = VR(:,j), the j-th column of VR. If the j-th and
   *          (j+1)-th eigenvalues form a complex conjugate pair, then
   *          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
   *          Each eigenvector will be scaled so the largest component have
   *          abs(real part) + abs(imag. part) = 1.
   *          Not referenced if JOBVR = 'N'.
   *
   *  LDVR    (input) INTEGER
   *          The leading dimension of the matrix VR. LDVR >= 1, and
   *          if JOBVR = 'V', LDVR >= N.
   *
   *  ILO     (output) INTEGER
   *  IHI     (output) INTEGER
   *          ILO and IHI are integer values such that on exit
   *          A(i,j) = 0 and B(i,j) = 0 if i > j and
   *          j = 1,...,ILO-1 or i = IHI+1,...,N.
   *          If BALANC = 'N' or 'S', ILO = 1 and IHI = N.
   *
   *  LSCALE  (output) DOUBLE PRECISION array, dimension (N)
   *          Details of the permutations and scaling factors applied
   *          to the left side of A and B.  If PL(j) is the index of the
   *          row interchanged with row j, and DL(j) is the scaling
   *          factor applied to row j, then
   *            LSCALE(j) = PL(j)  for j = 1,...,ILO-1
   *                      = DL(j)  for j = ILO,...,IHI
   *                      = PL(j)  for j = IHI+1,...,N.
   *          The order in which the interchanges are made is N to IHI+1,
   *          then 1 to ILO-1.
   *
   *  RSCALE  (output) DOUBLE PRECISION array, dimension (N)
   *          Details of the permutations and scaling factors applied
   *          to the right side of A and B.  If PR(j) is the index of the
   *          column interchanged with column j, and DR(j) is the scaling
   *          factor applied to column j, then
   *            RSCALE(j) = PR(j)  for j = 1,...,ILO-1
   *                      = DR(j)  for j = ILO,...,IHI
   *                      = PR(j)  for j = IHI+1,...,N
   *          The order in which the interchanges are made is N to IHI+1,
   *          then 1 to ILO-1.
   *
   *  ABNRM   (output) DOUBLE PRECISION
   *          The one-norm of the balanced matrix A.
   *
   *  BBNRM   (output) DOUBLE PRECISION
   *          The one-norm of the balanced matrix B.
   *
   *  RCONDE  (output) DOUBLE PRECISION array, dimension (N)
   *          If SENSE = 'E' or 'B', the reciprocal condition numbers of
   *          the eigenvalues, stored in consecutive elements of the array.
   *          For a complex conjugate pair of eigenvalues two consecutive
   *          elements of RCONDE are set to the same value. Thus RCONDE(j),
   *          RCONDV(j), and the j-th columns of VL and VR all correspond
   *          to the j-th eigenpair.
   *          If SENSE = 'N or 'V', RCONDE is not referenced.
   *
   *  RCONDV  (output) DOUBLE PRECISION array, dimension (N)
   *          If SENSE = 'V' or 'B', the estimated reciprocal condition
   *          numbers of the eigenvectors, stored in consecutive elements
   *          of the array. For a complex eigenvector two consecutive
   *          elements of RCONDV are set to the same value. If the
   *          eigenvalues cannot be reordered to compute RCONDV(j),
   *          RCONDV(j) is set to 0; this can only occur when the true
   *          value would be very small anyway.
   *          If SENSE = 'N' or 'E', RCONDV is not referenced.
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK. LWORK >= max(1,2*N).
   *          If BALANC = 'S' or 'B', or JOBVL = 'V', or JOBVR = 'V',
   *          LWORK >= max(1,6*N).
   *          If SENSE = 'E' or 'B', LWORK >= max(1,10*N).
   *          If SENSE = 'V' or 'B', LWORK >= 2*N*N+8*N+16.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  IWORK   (workspace) INTEGER array, dimension (N+6)
   *          If SENSE = 'E', IWORK is not referenced.
   *
   *  BWORK   (workspace) LOGICAL array, dimension (N)
   *          If SENSE = 'N', BWORK is not referenced.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value.
   *          = 1,...,N:
   *                The QZ iteration failed.  No eigenvectors have been
   *                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
   *                should be correct for j=INFO+1,...,N.
   *          > N:  =N+1: other than QZ iteration failed in DHGEQZ.
   *                =N+2: error return from DTGEVC.
   *
   *  Further Details
   *  ===============
   *
   *  Balancing a matrix pair (A,B) includes, first, permuting rows and
   *  columns to isolate eigenvalues, second, applying diagonal similarity
   *  transformation to the rows and columns to make the rows and columns
   *  as close in norm as possible. The computed reciprocal condition
   *  numbers correspond to the balanced matrix. Permuting rows and columns
   *  will not change the condition numbers (in exact arithmetic) but
   *  diagonal scaling will.  For further explanation of balancing, see
   *  section 4.11.1.2 of LAPACK Users' Guide.
   *
   *  An approximate error bound on the chordal distance between the i-th
   *  computed generalized eigenvalue w and the corresponding exact
   *  eigenvalue lambda is
   *
   *       chord(w, lambda) <= EPS * norm(ABNRM, BBNRM) / RCONDE(I)
   *
   *  An approximate error bound for the angle between the i-th computed
   *  eigenvector VL(i) or VR(i) is given by
   *
   *       EPS * norm(ABNRM, BBNRM) / DIF(i).
   *
   *  For further explanation of the reciprocal condition numbers RCONDE
   *  and RCONDV, see section 4.11 of LAPACK User's Guide.
   *
  \*/

  inline
  integer
  ggevx(
    BalanceType const & balanc,
    bool                jobvl,
    bool                jobvr,
    SenseType   const & sense,
    integer             n,
    real                A[],
    integer             lda,
    real                B[],
    integer             ldb,
    real                alphar[],
    real                alphai[],
    real                beta[],
    real                vl[],
    integer             ldvl,
    real                vr[],
    integer             ldvr,
    integer           & ilo,
    integer           & ihi,
    real                lscale[],
    real                rscale[],
    real              & abnrm,
    real              & bbnrm,
    real                rconde[],
    real                rcondv[],
    real                work[],
    integer             lwork,
    integer             iwork[],
    integer             bwork[]
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sggevx)(
      balance_blas[balanc],
      (jobvl?"V":"N"),
      (jobvr?"V":"N"),
      sense_blas[sense],
      &n,
      A, &lda,
      B, &ldb,
      alphar, alphai, beta,
      vl, &ldvl, vr, &ldvr,
      &ilo, &ihi,
      lscale, rscale,
      &abnrm, &bbnrm,
      rconde, rcondv,
      work, &lwork, iwork, bwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
    LAPACK_F77NAME(sggevx)(
      const_cast<character*>(balance_blas[balanc]),
      const_cast<character*>(jobvl?"V":"N"),
      const_cast<character*>(jobvr?"V":"N"),
      const_cast<character*>(sense_blas[sense]),
      &n,
      A, &lda,
      B, &ldb,
      alphar, alphai, beta,
      vl, &ldvl, vr, &ldvr,
      &ilo, &ihi,
      lscale, rscale,
      &abnrm, &bbnrm,
      rconde, rcondv,
      work, &lwork, iwork, bwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sggevx(
      balance_blas[balanc],
      (jobvl?"V":"N"),
      (jobvr?"V":"N"),
      sense_blas[sense],
      &n,
      A, &lda,
      B, &ldb,
      alphar, alphai, beta,
      vl, &ldvl, vr, &ldvr,
      &ilo, &ihi,
      lscale, rscale,
      &abnrm, &bbnrm,
      rconde, rcondv,
      work, &lwork, iwork, bwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sggevx)(
      const_cast<character*>(balance_blas[balanc]),
      const_cast<character*>(jobvl?"V":"N"),
      const_cast<character*>(jobvr?"V":"N"),
      const_cast<character*>(sense_blas[sense]),
      &n,
      const_cast<real*>(A), &lda,
      const_cast<real*>(B), &ldb,
      alphar, alphai, beta,
      vl, &ldvl, vr, &ldvr,
      &ilo, &ihi,
      lscale, rscale,
      &abnrm, &bbnrm,
      rconde, rcondv,
      work, &lwork, iwork, bwork, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  inline
  integer
  ggevx(
    BalanceType const & balanc,
    bool                jobvl,
    bool                jobvr,
    SenseType   const & sense,
    integer             n,
    doublereal          A[],
    integer	            lda,
    doublereal          B[],
    integer             ldb,
    doublereal          alphar[],
    doublereal          alphai[],
    doublereal          beta[],
    doublereal          vl[],
    integer             ldvl,
    doublereal          vr[],
    integer             ldvr,
    integer           & ilo,
    integer           & ihi,
    doublereal          lscale[],
    doublereal          rscale[],
    doublereal        & abnrm,
    doublereal        & bbnrm,
    doublereal          rconde[],
    doublereal          rcondv[],
    doublereal          work[],
    integer             lwork,
    integer             iwork[],
    integer             bwork[]
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dggevx)(
      balance_blas[balanc],
      (jobvl?"V":"N"),
      (jobvr?"V":"N"),
      sense_blas[sense],
      &n,
      A, &lda,
      B, &ldb,
      alphar, alphai, beta,
      vl, &ldvl, vr, &ldvr,
      &ilo, &ihi,
      lscale, rscale,
      &abnrm, &bbnrm,
      rconde, rcondv,
      work, &lwork, iwork, bwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
    LAPACK_F77NAME(dggevx)(
      const_cast<character*>(balance_blas[balanc]),
      const_cast<character*>(jobvl?"V":"N"),
      const_cast<character*>(jobvr?"V":"N"),
      const_cast<character*>(sense_blas[sense]),
      &n,
      A, &lda,
      B, &ldb,
      alphar, alphai, beta,
      vl, &ldvl, vr, &ldvr,
      &ilo, &ihi,
      lscale, rscale,
      &abnrm, &bbnrm,
      rconde, rcondv,
      work, &lwork, iwork, bwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dggevx(
      balance_blas[balanc],
      (jobvl?"V":"N"),
      (jobvr?"V":"N"),
      sense_blas[sense],
      &n,
      A, &lda,
      B, &ldb,
      alphar, alphai, beta,
      vl, &ldvl, vr, &ldvr,
      &ilo, &ihi,
      lscale, rscale,
      &abnrm, &bbnrm,
      rconde, rcondv,
      work, &lwork, iwork, bwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dggevx)(
      const_cast<character*>(balance_blas[balanc]),
      const_cast<character*>(jobvl?"V":"N"),
      const_cast<character*>(jobvr?"V":"N"),
      const_cast<character*>(sense_blas[sense]),
      &n,
      const_cast<doublereal*>(A), &lda,
      const_cast<doublereal*>(B), &ldb,
      alphar, alphai, beta,
      vl, &ldvl, vr, &ldvr,
      &ilo, &ihi,
      lscale, rscale,
      &abnrm, &bbnrm,
      rconde, rcondv,
      work, &lwork, iwork, bwork, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  /*
  //   ____  ______     __
  //  / ___||  _ \ \   / /
  //  \___ \| | | \ \ / /
  //   ___) | |_| |\ V /
  //  |____/|____/  \_/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGESVD computes the singular value decomposition (SVD) of a real
   *  M-by-N matrix A, optionally computing the left and/or right singular
   *  vectors. The SVD is written
   *
   *       A = U * SIGMA * transpose(V)
   *
   *  where SIGMA is an M-by-N matrix which is zero except for its
   *  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
   *  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
   *  are the singular values of A; they are real and non-negative, and
   *  are returned in descending order.  The first min(m,n) columns of
   *  U and V are the left and right singular vectors of A.
   *
   *  Note that the routine returns V**T, not V.
   *
   *  Arguments
   *  =========
   *
   *  JOBU    (input) CHARACTER*1
   *          Specifies options for computing all or part of the matrix U:
   *          = 'A':  all M columns of U are returned in array U:
   *          = 'S':  the first min(m,n) columns of U (the left singular
   *                  vectors) are returned in the array U;
   *          = 'O':  the first min(m,n) columns of U (the left singular
   *                  vectors) are overwritten on the array A;
   *          = 'N':  no columns of U (no left singular vectors) are
   *                  computed.
   *
   *  JOBVT   (input) CHARACTER*1
   *          Specifies options for computing all or part of the matrix
   *          V**T:
   *          = 'A':  all N rows of V**T are returned in the array VT;
   *          = 'S':  the first min(m,n) rows of V**T (the right singular
   *                  vectors) are returned in the array VT;
   *          = 'O':  the first min(m,n) rows of V**T (the right singular
   *                  vectors) are overwritten on the array A;
   *          = 'N':  no rows of V**T (no right singular vectors) are
   *                  computed.
   *
   *          JOBVT and JOBU cannot both be 'O'.
   *
   *  M       (input) INTEGER
   *          The number of rows of the input matrix A.  M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the input matrix A.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the M-by-N matrix A.
   *          On exit,
   *          if JOBU = 'O',  A is overwritten with the first min(m,n)
   *                          columns of U (the left singular vectors,
   *                          stored columnwise);
   *          if JOBVT = 'O', A is overwritten with the first min(m,n)
   *                          rows of V**T (the right singular vectors,
   *                          stored rowwise);
   *          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
   *                          are destroyed.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
   *          The singular values of A, sorted so that S(i) >= S(i+1).
   *
   *  U       (output) DOUBLE PRECISION array, dimension (LDU,UCOL)
   *          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
   *          If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
   *          if JOBU = 'S', U contains the first min(m,n) columns of U
   *          (the left singular vectors, stored columnwise);
   *          if JOBU = 'N' or 'O', U is not referenced.
   *
   *  LDU     (input) INTEGER
   *          The leading dimension of the array U.  LDU >= 1; if
   *          JOBU = 'S' or 'A', LDU >= M.
   *
   *  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)
   *          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
   *          V**T;
   *          if JOBVT = 'S', VT contains the first min(m,n) rows of
   *          V**T (the right singular vectors, stored rowwise);
   *          if JOBVT = 'N' or 'O', VT is not referenced.
   *
   *  LDVT    (input) INTEGER
   *          The leading dimension of the array VT.  LDVT >= 1; if
   *          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
   *          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
   *          superdiagonal elements of an upper bidiagonal matrix B
   *          whose diagonal is in S (not necessarily sorted). B
   *          satisfies A = U * B * VT, so it has the same singular values
   *          as A, and singular vectors related by U and VT.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK.
   *          LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)).
   *          For good performance, LWORK should generally be larger.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit.
   *          < 0:  if INFO = -i, the i-th argument had an illegal value.
   *          > 0:  if DBDSQR did not converge, INFO specifies how many
   *                superdiagonals of an intermediate bidiagonal form B
   *                did not converge to zero. See the description of WORK
   *                above for details.
   *
   *  =====================================================================
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  // use standard Lapack routine
  extern "C" {

    void
    LAPACK_F77NAME(dgesvd)(
      character const   JOBU[],
      character const   JOBVT[],
      integer   const * M,
      integer   const * N,
      doublereal      * A,
      integer   const * LDA,
      doublereal      * S,
      doublereal      * U,
      integer   const * LDU,
      doublereal      * VT,
      integer   const * LDVT,
      doublereal      * WORK,
      integer   const * LWORK,
      integer         * info
    );

    void
    LAPACK_F77NAME(sgesvd)(
      character const   JOBU[],
      character const   JOBVT[],
      integer   const * M,
      integer   const * N,
      real            * A,
      integer   const * LDA,
      real            * S,
      real            * U,
      integer   const * LDU,
      real            * VT,
      integer   const * LDVT,
      real            * WORK,
      integer   const * LWORK,
      integer         * info
    );
  }
  #endif

  inline
  integer
  gesvd(
    JobType const & JOBU,
    JobType const & JOBVT,
    integer         M,
    integer         N,
    real            A[],
    integer         LDA,
    real            S[],
    real            U[],
    integer         LDU,
    real            VT[],
    integer         LDVT,
    real            WORK[],
    integer         LWORK
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)
    LAPACK_F77NAME(sgesvd)(
      job_blas[JOBU],
      job_blas[JOBVT],
      &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT,
      WORK, &LWORK, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
          defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sgesvd)(
      const_cast<character*>(job_blas[JOBU]),
      const_cast<character*>(job_blas[JOBVT]),
      &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT,
      WORK, &LWORK, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgesvd(
      job_blas[JOBU], job_blas[JOBVT],
      &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT,
      WORK, &LWORK, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgesvd)(
      const_cast<character*>(job_blas[JOBU]),
      const_cast<character*>(job_blas[JOBVT]),
      &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT,
      WORK, &LWORK, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  inline
  integer
  gesvd(
    JobType const & JOBU,
    JobType const & JOBVT,
    integer         M,
    integer         N,
    doublereal      A[],
    integer         LDA,
    doublereal      S[],
    doublereal      U[],
    integer         LDU,
    doublereal      VT[],
    integer         LDVT,
    doublereal      WORK[],
    integer         LWORK
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dgesvd)(
      const_cast<character*>(job_blas[JOBU]),
      const_cast<character*>(job_blas[JOBVT]),
      &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT,
      WORK, &LWORK, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgesvd(
      job_blas[JOBU], job_blas[JOBVT],
      &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT,
      WORK, &LWORK, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgesvd)(
      const_cast<character*>(job_blas[JOBU]),
      const_cast<character*>(job_blas[JOBVT]),
      &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT,
      WORK, &LWORK, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  /*
  //                       _     _
  //    __ _  ___  ___  __| | __| |
  //   / _` |/ _ \/ __|/ _` |/ _` |
  //  | (_| |  __/\__ \ (_| | (_| |
  //   \__, |\___||___/\__,_|\__,_|
  //   |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGESDD computes the singular value decomposition (SVD) of a real
   *  M-by-N matrix A, optionally computing the left and right singular
   *  vectors.  If singular vectors are desired, it uses a
   *  divide-and-conquer algorithm.
   *
   *  The SVD is written
   *
   *       A = U * SIGMA * transpose(V)
   *
   *  where SIGMA is an M-by-N matrix which is zero except for its
   *  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
   *  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
   *  are the singular values of A; they are real and non-negative, and
   *  are returned in descending order.  The first min(m,n) columns of
   *  U and V are the left and right singular vectors of A.
   *
   *  Note that the routine returns VT = V**T, not V.
   *
   *  The divide and conquer algorithm makes very mild assumptions about
   *  floating point arithmetic. It will work on machines with a guard
   *  digit in add/subtract, or on those binary machines without guard
   *  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
   *  Cray-2. It could conceivably fail on hexadecimal or decimal machines
   *  without guard digits, but we know of none.
   *
   *  Arguments
   *  =========
   *
   *  JOBZ    (input) CHARACTER*1
   *          Specifies options for computing all or part of the matrix U:
   *          = 'A':  all M columns of U and all N rows of V**T are
   *                  returned in the arrays U and VT;
   *          = 'S':  the first min(M,N) columns of U and the first
   *                  min(M,N) rows of V**T are returned in the arrays U
   *                  and VT;
   *          = 'O':  If M >= N, the first N columns of U are overwritten
   *                  on the array A and all rows of V**T are returned in
   *                  the array VT;
   *                  otherwise, all columns of U are returned in the
   *                  array U and the first M rows of V**T are overwritten
   *                  in the array A;
   *          = 'N':  no columns of U or rows of V**T are computed.
   *
   *  M       (input) INTEGER
   *          The number of rows of the input matrix A.  M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the input matrix A.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the M-by-N matrix A.
   *          On exit,
   *          if JOBZ = 'O',  A is overwritten with the first N columns
   *                          of U (the left singular vectors, stored
   *                          columnwise) if M >= N;
   *                          A is overwritten with the first M rows
   *                          of V**T (the right singular vectors, stored
   *                          rowwise) otherwise.
   *          if JOBZ .ne. 'O', the contents of A are destroyed.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
   *          The singular values of A, sorted so that S(i) >= S(i+1).
   *
   *  U       (output) DOUBLE PRECISION array, dimension (LDU,UCOL)
   *          UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M < N;
   *          UCOL = min(M,N) if JOBZ = 'S'.
   *          If JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M
   *          orthogonal matrix U;
   *          if JOBZ = 'S', U contains the first min(M,N) columns of U
   *          (the left singular vectors, stored columnwise);
   *          if JOBZ = 'O' and M >= N, or JOBZ = 'N', U is not referenced.
   *
   *  LDU     (input) INTEGER
   *          The leading dimension of the array U.  LDU >= 1; if
   *          JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M.
   *
   *  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)
   *          If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the
   *          N-by-N orthogonal matrix V**T;
   *          if JOBZ = 'S', VT contains the first min(M,N) rows of
   *          V**T (the right singular vectors, stored rowwise);
   *          if JOBZ = 'O' and M < N, or JOBZ = 'N', VT is not referenced.
   *
   *  LDVT    (input) INTEGER
   *          The leading dimension of the array VT.  LDVT >= 1; if
   *          JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N;
   *          if JOBZ = 'S', LDVT >= min(M,N).
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK. LWORK >= 1.
   *          If JOBZ = 'N',
   *            LWORK >= 3*min(M,N) + max(max(M,N),7*min(M,N)).
   *          If JOBZ = 'O',
   *            LWORK >= 3*min(M,N)*min(M,N) +
   *                     max(max(M,N),5*min(M,N)*min(M,N)+4*min(M,N)).
   *          If JOBZ = 'S' or 'A'
   *            LWORK >= 3*min(M,N)*min(M,N) +
   *                     max(max(M,N),4*min(M,N)*min(M,N)+4*min(M,N)).
   *          For good performance, LWORK should generally be larger.
   *          If LWORK = -1 but other input arguments are legal, WORK(1)
   *          returns the optimal LWORK.
   *
   *  IWORK   (workspace) INTEGER array, dimension (8*min(M,N))
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit.
   *          < 0:  if INFO = -i, the i-th argument had an illegal value.
   *          > 0:  DBDSDC did not converge, updating process failed.
   *
   *  Further Details
   *  ===============
   *
   *  Based on contributions by
   *     Ming Gu and Huan Ren, Computer Science Division, University of
   *     California at Berkeley, USA
   *
   *  =====================================================================
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  // use standard Lapack routine
  extern "C" {

    void
    LAPACK_F77NAME(dgesdd)(
      character const   JOBZ[],
      integer   const * M,
      integer   const * N,
      doublereal      * A,
      integer   const * LDA,
      doublereal      * S,
      doublereal      * U,
      integer   const * LDU,
      doublereal      * VT,
      integer   const * LDVT,
      doublereal      * WORK,
      integer   const * LWORK,
      integer   const * IWORK, // (8*min(M,N))
      integer         * info
    );

    void
    LAPACK_F77NAME(sgesdd)(
      character const   JOBZ[],
      integer   const * M,
      integer   const * N,
      real            * A,
      integer   const * LDA,
      real            * S,
      real            * U,
      integer   const * LDU,
      real            * VT,
      integer   const * LDVT,
      real            * WORK,
      integer   const * LWORK,
      integer   const * IWORK, // (8*min(M,N))
      integer         * info
    );
  }
  #endif

  inline
  integer
  gesdd(
    JobType const & JOBZ,
    integer         M,
    integer         N,
    real            A[],
    integer         LDA,
    real            S[],
    real            U[],
    integer         LDU,
    real            VT[],
    integer         LDVT,
    real            WORK[],
    integer         LWORK,
    integer         IWORK[]
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sgesdd)(
      const_cast<character*>(job_blas[JOBZ]),
      &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT,
      WORK, &LWORK, IWORK, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgesdd(
      job_blas[JOBZ],
      &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, IWORK, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgesdd)(
      const_cast<character*>(job_blas[JOBZ]),
      &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, IWORK, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  inline
  integer
  gesdd(
    JobType const & JOBZ,
    integer         M,
    integer         N,
    doublereal      A[],
    integer         LDA,
    doublereal      S[],
    doublereal      U[],
    integer         LDU,
    doublereal      VT[],
    integer         LDVT,
    doublereal      WORK[],
    integer         LWORK,
    integer         IWORK[]
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dgesdd)(
      const_cast<character*>(job_blas[JOBZ]),
      &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT,
      WORK, &LWORK, IWORK, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgesdd(
      job_blas[JOBZ],
      &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, IWORK, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgesdd)(
      const_cast<character*>(job_blas[JOBZ]),
      &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, IWORK, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }


  /*
  //    ____ _____ _     ____  ____
  //   / ___| ____| |   / ___||  _ \
  //  | |  _|  _| | |   \___ \| | | |
  //  | |_| | |___| |___ ___) | |_| |
  //   \____|_____|_____|____/|____/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGELSD computes the minimum-norm solution to a real linear least
   *  squares problem:
   *      minimize 2-norm(| b - A*x |)
   *  using the singular value decomposition (SVD) of A. A is an M-by-N
   *  matrix which may be rank-deficient.
   *
   *  Several right hand side vectors b and solution vectors x can be
   *  handled in a single call; they are stored as the columns of the
   *  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
   *  matrix X.
   *
   *  The problem is solved in three steps:
   *  (1) Reduce the coefficient matrix A to bidiagonal form with
   *      Householder transformations, reducing the original problem
   *      into a "bidiagonal least squares problem" (BLS)
   *  (2) Solve the BLS using a divide and conquer approach.
   *  (3) Apply back all the Householder tranformations to solve
   *      the original least squares problem.
   *
   *  The effective rank of A is determined by treating as zero those
   *  singular values which are less than RCOND times the largest singular
   *  value.
   *
   *  The divide and conquer algorithm makes very mild assumptions about
   *  floating point arithmetic. It will work on machines with a guard
   *  digit in add/subtract, or on those binary machines without guard
   *  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
   *  Cray-2. It could conceivably fail on hexadecimal or decimal machines
   *  without guard digits, but we know of none.
   *
   *  Arguments
   *  =========
   *
   *  M       (input) INTEGER
   *          The number of rows of A. M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of A. N >= 0.
   *
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of columns
   *          of the matrices B and X. NRHS >= 0.
   *
   *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the M-by-N matrix A.
   *          On exit, A has been destroyed.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
   *          On entry, the M-by-NRHS right hand side matrix B.
   *          On exit, B is overwritten by the N-by-NRHS solution
   *          matrix X.  If m >= n and RANK = n, the residual
   *          sum-of-squares for the solution in the i-th column is given
   *          by the sum of squares of elements n+1:m in that column.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B. LDB >= max(1,max(M,N)).
   *
   *  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
   *          The singular values of A in decreasing order.
   *          The condition number of A in the 2-norm = S(1)/S(min(m,n)).
   *
   *  RCOND   (input) DOUBLE PRECISION
   *          RCOND is used to determine the effective rank of A.
   *          Singular values S(i) <= RCOND*S(1) are treated as zero.
   *          If RCOND < 0, machine precision is used instead.
   *
   *  RANK    (output) INTEGER
   *          The effective rank of A, i.e., the number of singular values
   *          which are greater than RCOND*S(1).
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK. LWORK must be at least 1.
   *          The exact minimum amount of workspace needed depends on M,
   *          N and NRHS. As long as LWORK is at least
   *              12*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2,
   *          if M is greater than or equal to N or
   *              12*M + 2*M*SMLSIZ + 8*M*NLVL + M*NRHS + (SMLSIZ+1)**2,
   *          if M is less than N, the code will execute correctly.
   *          SMLSIZ is returned by ILAENV and is equal to the maximum
   *          size of the subproblems at the bottom of the computation
   *          tree (usually about 25), and
   *             NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )
   *          For good performance, LWORK should generally be larger.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  Further Details
   *  ===============
   *
   *  Based on contributions by
   *     Ming Gu and Ren-Cang Li, Computer Science Division, University of
   *       California at Berkeley, USA
   *     Osni Marques, LBNL/NERSC, USA
   *
   *  =====================================================================
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  // use standard Lapack routine
  extern "C" {

    void
    LAPACK_F77NAME(sgelsd)(
      integer const * m,
      integer const * n,
      integer const * nrhs,
      real    const   a[],
      integer const * lda,
      real            b[],
      integer const * ldb,
      real            s[],
      real    const * rcond,
      integer       * rank,
      real            work[],
      integer const * lwork,
      integer         iwork[],
      integer       * info
    );

    void
    LAPACK_F77NAME(dgelsd)(
      integer    const * m,
      integer    const * n,
      integer    const * nrhs,
      doublereal const   a[],
      integer    const * lda,
      doublereal         b[],
      integer    const * ldb,
      doublereal         s[],
      doublereal const * rcond,
      integer          * rank,
      doublereal         work[],
      integer    const * lwork,
      integer            iwork[],
      integer          * info
    );
  }
  #endif

  inline
  integer
  gelsd(
    integer   m,
    integer   n,
    integer   nrhs,
    real      a[],
    integer   lda,
    real      b[],
    integer   ldb,
    real      s[],
    real      rcond,
    integer & rank,
    real      work[],
    integer   lwork,
    integer   iwork[]
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sgelsd)(
      &m, &n, &nrhs, a, &lda, b, &ldb,
      s, &rcond, &rank, work, &lwork, iwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgelsd( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info);
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgelsd)( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  inline
  integer
  gelsd(
    integer    m,
    integer    n,
    integer    nrhs,
    doublereal a[],
    integer    lda,
    doublereal b[],
    integer    ldb,
    doublereal s[],
    doublereal rcond,
    integer  & rank,
    doublereal work[],
    integer    lwork,
    integer    iwork[]
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dgelsd)(
      &m, &n, &nrhs, a, &lda, b, &ldb,
      s, &rcond, &rank, work, &lwork, iwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgelsd( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info);
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgelsd)( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  /*
  //    ____ _____ _     ____ ____
  //   / ___| ____| |   / ___/ ___|
  //  | |  _|  _| | |   \___ \___ \
  //  | |_| | |___| |___ ___) |__) |
  //   \____|_____|_____|____/____/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGELSS computes the minimum norm solution to a real linear least
   *  squares problem:
   *
   *  Minimize 2-norm(| b - A*x |).
   *
   *  using the singular value decomposition (SVD) of A. A is an M-by-N
   *  matrix which may be rank-deficient.
   *
   *  Several right hand side vectors b and solution vectors x can be
   *  handled in a single call; they are stored as the columns of the
   *  M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix
   *  X.
   *
   *  The effective rank of A is determined by treating as zero those
   *  singular values which are less than RCOND times the largest singular
   *  value.
   *
   *  Arguments
   *  =========
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix A. M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A. N >= 0.
   *
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of columns
   *          of the matrices B and X. NRHS >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the M-by-N matrix A.
   *          On exit, the first min(m,n) rows of A are overwritten with
   *          its right singular vectors, stored rowwise.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
   *          On entry, the M-by-NRHS right hand side matrix B.
   *          On exit, B is overwritten by the N-by-NRHS solution
   *          matrix X.  If m >= n and RANK = n, the residual
   *          sum-of-squares for the solution in the i-th column is given
   *          by the sum of squares of elements n+1:m in that column.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B. LDB >= max(1,max(M,N)).
   *
   *  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
   *          The singular values of A in decreasing order.
   *          The condition number of A in the 2-norm = S(1)/S(min(m,n)).
   *
   *  RCOND   (input) DOUBLE PRECISION
   *          RCOND is used to determine the effective rank of A.
   *          Singular values S(i) <= RCOND*S(1) are treated as zero.
   *          If RCOND < 0, machine precision is used instead.
   *
   *  RANK    (output) INTEGER
   *          The effective rank of A, i.e., the number of singular values
   *          which are greater than RCOND*S(1).
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK. LWORK >= 1, and also:
   *          LWORK >= 3*min(M,N) + max( 2*min(M,N), max(M,N), NRHS )
   *          For good performance, LWORK should generally be larger.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value.
   *          > 0:  the algorithm for computing the SVD failed to converge;
   *                if INFO = i, i off-diagonal elements of an intermediate
   *                bidiagonal form did not converge to zero.
   *
   *  =====================================================================
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  // use standard Lapack routine
  extern "C" {

    void
    LAPACK_F77NAME(sgelss)(
      integer const * m,
      integer const * n,
      integer const * nrhs,
      real            a[],
      integer const * lda,
      real            b[],
      integer const * ldb,
      real            s[],
      real    const * rcond,
      integer       * rank,
      real          * work,
      integer const * lwork,
      integer       * info
    );

    void
    LAPACK_F77NAME(dgelss)(
      integer    const * m,
      integer    const * n,
      integer    const * nrhs,
      doublereal         a[],
      integer    const * lda,
      doublereal         b[],
      integer    const * ldb,
      doublereal         s[],
      doublereal const * rcond,
      integer          * rank,
      doublereal       * work,
      integer    const * lwork,
      integer          * info
    );
  }
  #endif

  inline
  integer
  gelss(
    integer   m,
    integer   n,
    integer   nrhs,
    real      a[],
    integer   lda,
    real      b[],
    integer   ldb,
    real      s[],
    real      rcond,
    integer & rank,
    real      work[],
    integer   lwork
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sgelss)(
      &m, &n, &nrhs, a, &lda, b, &ldb, s,
      &rcond, &rank, work, &lwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgelss(
      &m, &n, &nrhs, a, &lda, b, &ldb, s,
      &rcond, &rank, work, &lwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgelss)(
      &m, &n, &nrhs, a, &lda, b, &ldb, s,
      &rcond, &rank, work, &lwork, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  inline
  integer
  gelss(
    integer    m,
    integer    n,
    integer    nrhs,
    doublereal a[],
    integer    lda,
    doublereal b[],
    integer    ldb,
    doublereal s[],
    doublereal rcond,
    integer  & rank,
    doublereal work[],
    integer    lwork
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dgelss)(
      &m, &n, &nrhs, a, &lda, b, &ldb, s,
      &rcond, &rank, work, &lwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgelss(
      &m, &n, &nrhs, a, &lda, b, &ldb, s,
      &rcond, &rank, work, &lwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgelss)(
      &m, &n, &nrhs, a, &lda, b, &ldb, s,
      &rcond, &rank, work, &lwork, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  /*
  //              _
  //    __ _  ___| |___ _   _
  //   / _` |/ _ \ / __| | | |
  //  | (_| |  __/ \__ \ |_| |
  //   \__, |\___|_|___/\__, |
  //   |___/            |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  SGELSY computes the minimum-norm solution to a real linear least
   *  squares problem:
   *      minimize || A * X - B ||
   *  using a complete orthogonal factorization of A.  A is an M-by-N
   *  matrix which may be rank-deficient.
   *
   *  Several right hand side vectors b and solution vectors x can be
   *  handled in a single call; they are stored as the columns of the
   *  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
   *  matrix X.
   *
   *  The routine first computes a QR factorization with column pivoting:
   *      A * P = Q * [ R11 R12 ]
   *                  [  0  R22 ]
   *  with R11 defined as the largest leading submatrix whose estimated
   *  condition number is less than 1/RCOND.  The order of R11, RANK,
   *  is the effective rank of A.
   *
   *  Then, R22 is considered to be negligible, and R12 is annihilated
   *  by orthogonal transformations from the right, arriving at the
   *  complete orthogonal factorization:
   *     A * P = Q * [ T11 0 ] * Z
   *                 [  0  0 ]
   *  The minimum-norm solution is then
   *     X = P * Z' [ inv(T11)*Q1'*B ]
   *                [        0       ]
   *  where Q1 consists of the first RANK columns of Q.
   *
   *  This routine is basically identical to the original xGELSX except
   *  three differences:
   *    o The call to the subroutine xGEQPF has been substituted by the
   *      the call to the subroutine xGEQP3. This subroutine is a Blas-3
   *      version of the QR factorization with column pivoting.
   *    o Matrix B (the right hand side) is updated with Blas-3.
   *    o The permutation of matrix B (the right hand side) is faster and
   *      more simple.
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
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of
   *          columns of matrices B and X. NRHS >= 0.
   *
   *  A       (input/output) REAL array, dimension (LDA,N)
   *          On entry, the M-by-N matrix A.
   *          On exit, A has been overwritten by details of its
   *          complete orthogonal factorization.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  B       (input/output) REAL array, dimension (LDB,NRHS)
   *          On entry, the M-by-NRHS right hand side matrix B.
   *          On exit, the N-by-NRHS solution matrix X.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B. LDB >= max(1,M,N).
   *
   *  JPVT    (input/output) INTEGER array, dimension (N)
   *          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted
   *          to the front of AP, otherwise column i is a free column.
   *          On exit, if JPVT(i) = k, then the i-th column of AP
   *          was the k-th column of A.
   *
   *  RCOND   (input) REAL
   *          RCOND is used to determine the effective rank of A, which
   *          is defined as the order of the largest leading triangular
   *          submatrix R11 in the QR factorization with pivoting of A,
   *          whose estimated condition number < 1/RCOND.
   *
   *  RANK    (output) INTEGER
   *          The effective rank of A, i.e., the order of the submatrix
   *          R11.  This is the same as the order of the submatrix T11
   *          in the complete orthogonal factorization of A.
   *
   *  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK.
   *          The unblocked strategy requires that:
   *             LWORK >= MAX( MN+3*N+1, 2*MN+NRHS ),
   *          where MN = min( M, N ).
   *          The block algorithm requires that:
   *             LWORK >= MAX( MN+2*N+NB*(N+1), 2*MN+NB*NRHS ),
   *          where NB is an upper bound on the blocksize returned
   *          by ILAENV for the routines SGEQP3, STZRZF, STZRQF, SORMQR,
   *          and SORMRZ.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0: successful exit
   *          < 0: If INFO = -i, the i-th argument had an illegal value.
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  // use standard Lapack routine
  extern "C" {

    void
    LAPACK_F77NAME(sgelsy)(
      integer const * m,
      integer const * n,
      integer const * nrhs,
      real            a[],
      integer const * lda,
      real            b[],
      integer const * ldb,
      integer         jpvt[],
      real    const * rcond,
      integer       * rank,
      real          * work,
      integer const * lwork,
      integer       * info
    );

    void
    LAPACK_F77NAME(dgelsy)(
      integer    const * m,
      integer    const * n,
      integer    const * nrhs,
      doublereal         a[],
      integer    const * lda,
      doublereal         b[],
      integer    const * ldb,
      integer            jpvt[],
      doublereal const * rcond,
      integer          * rank,
      doublereal       * work,
      integer    const * lwork,
      integer          * info
    );
  }
  #endif

  inline
  integer
  gelsy(
    integer   m,
    integer   n,
    integer   nrhs,
    real      a[],
    integer   lda,
    real      b[],
    integer   ldb,
    integer   jpvt[],
    real      rcond,
    integer & rank,
    real      work[],
    integer   lwork
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)  || \
       defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
       defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sgelsy)(
      &m, &n, &nrhs, a, &lda, b, &ldb, jpvt,
      &rcond, &rank, work, &lwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgelsy(
      &m, &n, &nrhs, a, &lda, b, &ldb, jpvt,
      &rcond, &rank, work, &lwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgelsy)(
      &m, &n, &nrhs, a, &lda, b, &ldb, jpvt,
      &rcond, &rank, work, &lwork, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  inline
  integer
  gelsy(
    integer    m,
    integer    n,
    integer    nrhs,
    doublereal a[],
    integer    lda,
    doublereal b[],
    integer    ldb,
    integer    jpvt[],
    doublereal rcond,
    integer  & rank,
    doublereal work[],
    integer    lwork
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dgelsy)(
      &m, &n, &nrhs, a, &lda, b, &ldb, jpvt,
      &rcond, &rank, work, &lwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgelsy(
      &m, &n, &nrhs, a, &lda, b, &ldb, jpvt,
      &rcond, &rank, work, &lwork, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgelsy)(
      &m, &n, &nrhs, a, &lda, b, &ldb, jpvt,
      &rcond, &rank, work, &lwork, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  /*\
   *
   *  SUBROUTINE DGGSVD3( JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B,
   *                      LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK,
   *                      LWORK, IWORK, INFO )
   *
   *
   *  Purpose:
   *  =============
   *
   *  DGGSVD3 computes the generalized singular value decomposition (GSVD)
   *  of an M-by-N real matrix A and P-by-N real matrix B:
   *
   *        U**T*A*Q = D1*( 0 R ),    V**T*B*Q = D2*( 0 R )
   *
   *  where U, V and Q are orthogonal matrices.
   *  Let K+L = the effective numerical rank of the matrix (A**T,B**T)**T,
   *  then R is a K+L-by-K+L nonsingular upper triangular matrix, D1 and
   *  D2 are M-by-(K+L) and P-by-(K+L) "diagonal" matrices and of the
   *  following structures, respectively:
   *
   *  If M-K-L >= 0,
   *
   *                     K  L
   *        D1 =     K ( I  0 )
   *                 L ( 0  C )
   *             M-K-L ( 0  0 )
   *
   *                   K  L
   *        D2 =   L ( 0  S )
   *             P-L ( 0  0 )
   *
   *                 N-K-L  K    L
   *   ( 0 R ) = K (  0   R11  R12 )
   *             L (  0    0   R22 )
   *
   *  where
   *
   *     C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),
   *     S = diag( BETA(K+1),  ... , BETA(K+L) ),
   *     C**2 + S**2 = I.
   *
   *     R is stored in A(1:K+L,N-K-L+1:N) on exit.
   *
   *  If M-K-L < 0,
   *
   *                     K M-K K+L-M
   *          D1 =   K ( I  0    0   )
   *               M-K ( 0  C    0   )
   *
   *                       K M-K K+L-M
   *          D2 =   M-K ( 0  S    0  )
   *               K+L-M ( 0  0    I  )
   *                 P-L ( 0  0    0  )
   *
   *                      N-K-L  K   M-K  K+L-M
   *     ( 0 R ) =     K ( 0    R11  R12  R13  )
   *                 M-K ( 0     0   R22  R23  )
   *               K+L-M ( 0     0    0   R33  )
   *
   *  where
   *
   *    C = diag( ALPHA(K+1), ... , ALPHA(M) ),
   *    S = diag( BETA(K+1),  ... , BETA(M) ),
   *    C**2 + S**2 = I.
   *
   *    (R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N), and R33 is stored
   *    ( 0  R22 R23 )
   *    in B(M-K+1:L,N+M-K-L+1:N) on exit.
   *
   *  The routine computes C, S, R, and optionally the orthogonal
   *  transformation matrices U, V and Q.
   *
   *  In particular, if B is an N-by-N nonsingular matrix, then the GSVD of
   *  A and B implicitly gives the SVD of A*inv(B):
   *                       A*inv(B) = U*(D1*inv(D2))*V**T.
   *  If ( A**T,B**T)**T  has orthonormal columns, then the GSVD of A and B is
   *  also equal to the CS decomposition of A and B. Furthermore, the GSVD
   *  can be used to derive the solution of the eigenvalue problem:
   *                       A**T*A x = lambda* B**T*B x.
   *  In some literature, the GSVD of A and B is presented in the form
   *                   U**T*A*X = ( 0 D1 ),   V**T*B*X = ( 0 D2 )
   *  where U and V are orthogonal and X is nonsingular, D1 and D2 are
   *  ``diagonal''.  The former GSVD form can be converted to the latter
   *  form by taking the nonsingular matrix X as
   *
   *                       X = Q*( I   0    )
   *                             ( 0 inv(R) ).
   *
   *
   *  Arguments:
   *  ==========
   *
   *  [in] JOBU, CHARACTER*1
   *          = 'U':  Orthogonal matrix U is computed;
   *          = 'N':  U is not computed.
   *
   *  [in] JOBV, CHARACTER*1
   *          = 'V':  Orthogonal matrix V is computed;
   *          = 'N':  V is not computed.
   *
   *  [in] JOBQ, CHARACTER*1
   *          = 'Q':  Orthogonal matrix Q is computed;
   *          = 'N':  Q is not computed.
   *
   *  [in] M, INTEGER The number of rows of the matrix A.  M >= 0.
   *  [in] N, INTEGER The number of columns of the matrices A and B.  N >= 0.
   *
   *  [in] P, INTEGER The number of rows of the matrix B.  P >= 0.
   *  [out] K, INTEGER
   *  [out] L, INTEGER
   *        On exit, K and L specify the dimension of the subblocks
   *        described in Purpose. K + L = effective numerical rank of (A**T,B**T)**T.
   *
   *  [in,out] A, DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the M-by-N matrix A.
   *          On exit, A contains the triangular matrix R, or part of R.
   *          See Purpose for details.
   *
   *  [in] LDA, INTEGER The leading dimension of the array A. LDA >= max(1,M).
   *
   *  [in,out] B
   *          B is DOUBLE PRECISION array, dimension (LDB,N)
   *          On entry, the P-by-N matrix B.
   *          On exit, B contains the triangular matrix R if M-K-L < 0.
   *          See Purpose for details.
   *
   *  [in] LDB, INTEGER The leading dimension of the array B. LDB >= max(1,P).
   *  [out] ALPHA, DOUBLE PRECISION array, dimension (N)
   *
   *  [out] BETA, DOUBLE PRECISION array, dimension (N)
   *
   *  On exit, ALPHA and BETA contain the generalized singular value pairs of A and B;
   *    ALPHA(1:K) = 1,
   *    BETA(1:K)  = 0,
   *  and if M-K-L >= 0,
   *    ALPHA(K+1:K+L) = C,
   *    BETA(K+1:K+L)  = S,
   *  or if M-K-L < 0,
   *    ALPHA(K+1:M)=C, ALPHA(M+1:K+L)=0
   *    BETA(K+1:M) =S, BETA(M+1:K+L) =1
   *  and
   *    ALPHA(K+L+1:N) = 0
   *    BETA(K+L+1:N)  = 0
   *
   *  [out] U
   *          U is DOUBLE PRECISION array, dimension (LDU,M)
   *          If JOBU = 'U', U contains the M-by-M orthogonal matrix U.
   *          If JOBU = 'N', U is not referenced.
   *  [in] LDU, INTEGER
   *          The leading dimension of the array U. LDU >= max(1,M) if
   *          JOBU = 'U'; LDU >= 1 otherwise.
   *
   *  [out] V
   *          V is DOUBLE PRECISION array, dimension (LDV,P)
   *          If JOBV = 'V', V contains the P-by-P orthogonal matrix V.
   *          If JOBV = 'N', V is not referenced.
   *  [in] LDV
   *          LDV is INTEGER
   *          The leading dimension of the array V. LDV >= max(1,P) if
   *          JOBV = 'V'; LDV >= 1 otherwise.
   *  [out] Q
   *          Q is DOUBLE PRECISION array, dimension (LDQ,N)
   *          If JOBQ = 'Q', Q contains the N-by-N orthogonal matrix Q.
   *          If JOBQ = 'N', Q is not referenced.
   *
   *  [in] LDQ
   *          LDQ is INTEGER
   *          The leading dimension of the array Q. LDQ >= max(1,N) if
   *          JOBQ = 'Q'; LDQ >= 1 otherwise.
   *
   *  [out] WORK
   *          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *  [in] LWORK
   *          LWORK is INTEGER
   *          The dimension of the array WORK.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  [out] IWORK
   *          IWORK is INTEGER array, dimension (N)
   *          On exit, IWORK stores the sorting information. More
   *          precisely, the following loop will sort ALPHA
   *             for I = K+1, min(M,K+L)
   *                 swap ALPHA(I) and ALPHA(IWORK(I))
   *             endfor
   *          such that ALPHA(1) >= ALPHA(2) >= ... >= ALPHA(N).
   *
   *  [out] INFO
   *          INFO is INTEGER
   *          = 0:  successful exit.
   *          < 0:  if INFO = -i, the i-th argument had an illegal value.
   *          > 0:  if INFO = 1, the Jacobi-type procedure failed to
   *                converge.  For further details, see subroutine DTGSJA.
   *
   * Internal Parameters:
   * =========================
   *
   *  TOLA    DOUBLE PRECISION
   *  TOLB    DOUBLE PRECISION
   *          TOLA and TOLB are the thresholds to determine the effective
   *          rank of (A**T,B**T)**T. Generally, they are set to
   *                   TOLA = MAX(M,N)*norm(A)*MACHEPS,
   *                   TOLB = MAX(P,N)*norm(B)*MACHEPS.
   *          The size of TOLA and TOLB may affect the size of backward
   *          errors of the decomposition.
   *
   *  Authors:
   *  ========
   *
   *  \author Univ. of Tennessee
   *  \author Univ. of California Berkeley
   *  \author Univ. of Colorado Denver
   *  \author NAG Ltd.
   *
   *  \date August 2015
   *
   *  Contributors:
   *  ==================
   *
   *  Ming Gu and Huan Ren, Computer Science Division, University of
   *  California at Berkeley, USA
   *
   *
   *  Further Details:
   *  =====================
   *
   *  DGGSVD3 replaces the deprecated subroutine DGGSVD.
   *
   * =====================================================================
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  // use standard Lapack routine
  extern "C" {
    void
    LAPACK_F77NAME(sggsvd3)(
      character const JOBU[],
      character const JOBV[],
      character const JOBQ[],
      integer   const * m,
      integer   const * n,
      integer   const * p,
      integer         * k,
      integer         * l,
      real              A[],
      integer   const * ldA,
      real              B[],
      integer   const * ldB,
      real              alpha[],
      real              beta[],
      real              U[],
      integer   const * ldU,
      real              V[],
      integer   const * ldV,
      real              Q[],
      integer   const * ldQ,
      real            * work,
      integer   const * lwork,
      integer         * iwork,
      integer         * info
    );

    void
    LAPACK_F77NAME(dggsvd3)(
      character const JOBU[],
      character const JOBV[],
      character const JOBQ[],
      integer   const * m,
      integer   const * n,
      integer   const * p,
      integer         * k,
      integer         * l,
      doublereal        A[],
      integer   const * ldA,
      doublereal        B[],
      integer   const * ldB,
      doublereal        alpha[],
      doublereal        beta[],
      doublereal        U[],
      integer   const * ldU,
      doublereal        V[],
      integer   const * ldV,
      doublereal        Q[],
      integer   const * ldQ,
      doublereal      * work,
      integer   const * lwork,
      integer         * iwork,
      integer         * info
    );
  }
  #endif

  inline
  integer
  ggsvd(
    bool      JOBU,
    bool      JOBV,
    bool      JOBQ,
    integer   m,
    integer   n,
    integer   p,
    integer & k,
    integer & l,
    real      A[],
    integer   ldA,
    real      B[],
    integer   ldB,
    real      alpha[],
    real      beta[],
    real      U[],
    integer   ldU,
    real      V[],
    integer   ldV,
    real      Q[],
    integer   ldQ,
    real      work[],
    integer   lwork,
    integer   iwork[]
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sggsvd3)(
      const_cast<character*>( JOBU ? "U" : "N"),
      const_cast<character*>( JOBV ? "V" : "N"),
      const_cast<character*>( JOBQ ? "Q" : "N"),
      &m, &n, &p, &k, &l,
      A, &ldA,
      B, &ldB,
      alpha, beta,
      U, &ldU,
      V, &ldV,
      Q, &ldQ,
      work, &lwork,
      iwork,
      &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sggsvd3(
      const_cast<character*>( JOBU ? "U" : "N"),
      const_cast<character*>( JOBV ? "V" : "N"),
      const_cast<character*>( JOBQ ? "Q" : "N"),
      &m, &n, &p, &k, &l,
      A, &ldA,
      B, &ldB,
      alpha, beta,
      U, &ldU,
      V, &ldV,
      Q, &ldQ,
      work, &lwork,
      iwork,
      &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    integer lw = std::max(3*n,std::max(m,p))+n;
    if ( lwork < 0 ) {
      work[0] = real(lw);
    } else {
      LAPACK_WRAPPER_ASSERT(
        lwork >= lw, "ggsvd, lwork = " << lwork << " must be >= " << lw
      );
      CLAPACKNAME(sggsvd)(
        const_cast<character*>( JOBU ? "U" : "N"),
        const_cast<character*>( JOBV ? "V" : "N"),
        const_cast<character*>( JOBQ ? "Q" : "N"),
        &m, &n, &p, &k, &l,
        A, &ldA,
        B, &ldB,
        alpha, beta,
        U, &ldU,
        V, &ldV,
        Q, &ldQ,
        work,
        iwork,
        &info
      );
    }
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
      #ifndef LAPACK_WRAPPER_OS_OSX
      LAPACK_F77NAME(sggsvd3)(
        const_cast<character*>( JOBU ? "U" : "N"),
        const_cast<character*>( JOBV ? "V" : "N"),
        const_cast<character*>( JOBQ ? "Q" : "N"),
        &m, &n, &p, &k, &l,
        A, &ldA,
        B, &ldB,
        alpha, beta,
        U, &ldU,
        V, &ldV,
        Q, &ldQ,
        work, &lwork,
        iwork,
        &info
      );
      #else
      integer lw = std::max(3*n,std::max(m,p))+n;
      if ( lwork < 0 ) {
        work[0] = real(lw);
      } else {
        LAPACK_WRAPPER_ASSERT(
          lwork >= lw, "ggsvd, lwork = " << lwork << " must be >= " << lw
        );
        LAPACK_F77NAME(sggsvd)(
          const_cast<character*>( JOBU ? "U" : "N"),
          const_cast<character*>( JOBV ? "V" : "N"),
          const_cast<character*>( JOBQ ? "Q" : "N"),
          &m, &n, &p, &k, &l,
          A, &ldA,
          B, &ldB,
          alpha, beta,
          U, &ldU,
          V, &ldV,
          Q, &ldQ,
          work,
          iwork,
          &info
        );
      }
      #endif
    #elif defined(LAPACK_WRAPPER_USE_LAPACK)
    LAPACK_F77NAME(sggsvd3)(
      ( JOBU ? "U" : "N"),
      ( JOBV ? "V" : "N"),
      ( JOBQ ? "Q" : "N"),
      &m, &n, &p, &k, &l,
      A, &ldA,
      B, &ldB,
      alpha, beta,
      U, &ldU,
      V, &ldV,
      Q, &ldQ,
      work, &lwork,
      iwork,
      &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  inline
  integer
  ggsvd(
    bool       JOBU,
    bool       JOBV,
    bool       JOBQ,
    integer    m,
    integer    n,
    integer    p,
    integer  & k,
    integer  & l,
    doublereal A[],
    integer    ldA,
    doublereal B[],
    integer    ldB,
    doublereal alpha[],
    doublereal beta[],
    doublereal U[],
    integer    ldU,
    doublereal V[],
    integer    ldV,
    doublereal Q[],
    integer    ldQ,
    doublereal work[],
    integer    lwork,
    integer    iwork[]
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dggsvd3)(
      const_cast<character*>( JOBU ? "U" : "N"),
      const_cast<character*>( JOBV ? "V" : "N"),
      const_cast<character*>( JOBQ ? "Q" : "N"),
      &m, &n, &p, &k, &l,
      A, &ldA,
      B, &ldB,
      alpha, beta,
      U, &ldU,
      V, &ldV,
      Q, &ldQ,
      work, &lwork,
      iwork,
      &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dggsvd3(
      const_cast<character*>( JOBU ? "U" : "N"),
      const_cast<character*>( JOBV ? "V" : "N"),
      const_cast<character*>( JOBQ ? "Q" : "N"),
      &m, &n, &p, &k, &l,
      A, &ldA,
      B, &ldB,
      alpha, beta,
      U, &ldU,
      V, &ldV,
      Q, &ldQ,
      work, &lwork,
      iwork,
      &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    integer lw = std::max(3*n,std::max(m,p))+n;
    if ( lwork < 0 ) {
      work[0] = doublereal(lw);
    } else {
      LAPACK_WRAPPER_ASSERT(
        lwork >= lw, "ggsvd, lwork = " << lwork << " must be >= " << lw
      );
      CLAPACKNAME(dggsvd)(
        const_cast<character*>( JOBU ? "U" : "N"),
        const_cast<character*>( JOBV ? "V" : "N"),
        const_cast<character*>( JOBQ ? "Q" : "N"),
        &m, &n, &p, &k, &l,
        A, &ldA,
        B, &ldB,
        alpha, beta,
        U, &ldU,
        V, &ldV,
        Q, &ldQ,
        work,
        iwork,
        &info
      );
    }
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
      #ifndef LAPACK_WRAPPER_OS_OSX
      LAPACK_F77NAME(dggsvd3)(
        const_cast<character*>( JOBU ? "U" : "N"),
        const_cast<character*>( JOBV ? "V" : "N"),
        const_cast<character*>( JOBQ ? "Q" : "N"),
        &m, &n, &p, &k, &l,
        A, &ldA,
        B, &ldB,
        alpha, beta,
        U, &ldU,
        V, &ldV,
        Q, &ldQ,
        work, &lwork,
        iwork,
        &info
      );
      #else
      integer lw = std::max(3*n,std::max(m,p))+n;
      if ( lwork < 0 ) {
        work[0] = real(lw);
      } else {
        LAPACK_WRAPPER_ASSERT(
          lwork >= lw, "ggsvd, lwork = " << lwork << " must be >= " << lw
        );
        LAPACK_F77NAME(dggsvd)(
          const_cast<character*>( JOBU ? "U" : "N"),
          const_cast<character*>( JOBV ? "V" : "N"),
          const_cast<character*>( JOBQ ? "Q" : "N"),
          &m, &n, &p, &k, &l,
          A, &ldA,
          B, &ldB,
          alpha, beta,
          U, &ldU,
          V, &ldV,
          Q, &ldQ,
          work,
          iwork,
          &info
        );
      }
      #endif
    #elif defined(LAPACK_WRAPPER_USE_LAPACK)
    LAPACK_F77NAME(dggsvd3)(
      ( JOBU ? "U" : "N"),
      ( JOBV ? "V" : "N"),
      ( JOBQ ? "Q" : "N"),
      &m, &n, &p, &k, &l,
      A, &ldA,
      B, &ldB,
      alpha, beta,
      U, &ldU,
      V, &ldV,
      Q, &ldQ,
      work, &lwork,
      iwork,
      &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

}
