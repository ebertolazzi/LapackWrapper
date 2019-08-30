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
/// file: tridiagonal.hxx
///

namespace lapack_wrapper {

  /*
  //   _____     _     _ _                               _
  //  |_   _| __(_) __| (_) __ _  __ _  ___  _ __   __ _| |
  //    | || '__| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | |
  //    | || |  | | (_| | | (_| | (_| | (_) | | | | (_| | |
  //    |_||_|  |_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|
  //                             |___/
  */

  /*\
   *  Purpose
   *  =======
   *
   *  DGTTRF computes an LU factorization of a real tridiagonal matrix A
   *  using elimination with partial pivoting and row interchanges.
   *
   *  The factorization has the form
   *     A = L * U
   *  where L is a product of permutation and unit lower bidiagonal
   *  matrices and U is upper triangular with nonzeros in only the main
   *  diagonal and first two superdiagonals.
   *
   *  Arguments
   *  =========
   *
   *  N       (input) INTEGER
   *          The order of the matrix A.
   *
   *  DL      (input/output) DOUBLE PRECISION array, dimension (N-1)
   *          On entry, DL must contain the (n-1) sub-diagonal elements of A.
   *
   *          On exit, DL is overwritten by the (n-1) multipliers that
   *          define the matrix L from the LU factorization of A.
   *
   *  D       (input/output) DOUBLE PRECISION array, dimension (N)
   *          On entry, D must contain the diagonal elements of A.
   *
   *          On exit, D is overwritten by the n diagonal elements of the
   *          upper triangular matrix U from the LU factorization of A.
   *
   *  DU      (input/output) DOUBLE PRECISION array, dimension (N-1)
   *          On entry, DU must contain the (n-1) super-diagonal elements
   *          of A.
   *
   *          On exit, DU is overwritten by the (n-1) elements of the first
   *          super-diagonal of U.
   *
   *  DU2     (output) DOUBLE PRECISION array, dimension (N-2)
   *          On exit, DU2 is overwritten by the (n-2) elements of the
   *          second super-diagonal of U.
   *
   *  IPIV    (output) INTEGER array, dimension (N)
   *          The pivot indices; for 1 <= i <= n, row i of the matrix was
   *          interchanged with row IPIV(i).  IPIV(i) will always be either
   *          i or i+1; IPIV(i) = i indicates a row interchange was not
   *          required.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -k, the k-th argument had an illegal value
   *          > 0:  if INFO = k, U(k,k) is exactly zero. The factorization
   *                has been completed, but the factor U is exactly
   *                singular, and division by zero will occur if it is used
   *                to solve a system of equations.
   *
   *  =====================================================================
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {

    void
    LAPACK_F77NAME(sgttrf)(
      integer const * N,
      real          * DL,
      real          * D,
      real          * DU,
      real          * DU2,
      integer       * IPIV,
      integer       * INFO
    );

    void
    LAPACK_F77NAME(dgttrf)(
      integer const * N,
      doublereal    * DL,
      doublereal    * D,
      doublereal    * DU,
      doublereal    * DU2,
      integer       * IPIV,
      integer       * INFO
    );
  }
  #endif

  inline
  integer
  gttrf(
    integer N,
    real    DL[],
    real    D[],
    real    DU[],
    real    DU2[],
    integer IPIV[]
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgttrf)( &N, DL, D, DU, DU2, IPIV, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_LAPACK)   || \
          defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
          defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sgttrf)( &N, DL, D, DU, DU2, IPIV, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgttrf( &N, DL, D, DU, DU2, IPIV, &INFO );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  inline
  integer
  gttrf(
    integer    N,
    doublereal DL[],
    doublereal D[],
    doublereal DU[],
    doublereal DU2[],
    integer    IPIV[]
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgttrf)( &N, DL, D, DU, DU2, IPIV, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_LAPACK)   || \
          defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
          defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dgttrf)( &N, DL, D, DU, DU2, IPIV, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgttrf( &N, DL, D, DU, DU2, IPIV, &INFO );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  /*\
   *  Purpose
   *  =======
   *
   *  DGTTRS solves one of the systems of equations
   *     A*X = B  or  A**T*X = B,
   *  with a tridiagonal matrix A using the LU factorization computed
   *  by DGTTRF.
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
   *          The order of the matrix A.
   *
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of columns
   *          of the matrix B.  NRHS >= 0.
   *
   *  DL      (input) DOUBLE PRECISION array, dimension (N-1)
   *          The (n-1) multipliers that define the matrix L from the
   *          LU factorization of A.
   *
   *  D       (input) DOUBLE PRECISION array, dimension (N)
   *          The n diagonal elements of the upper triangular matrix U from
   *          the LU factorization of A.
   *
   *  DU      (input) DOUBLE PRECISION array, dimension (N-1)
   *          The (n-1) elements of the first super-diagonal of U.
   *
   *  DU2     (input) DOUBLE PRECISION array, dimension (N-2)
   *          The (n-2) elements of the second super-diagonal of U.
   *
   *  IPIV    (input) INTEGER array, dimension (N)
   *          The pivot indices; for 1 <= i <= n, row i of the matrix was
   *          interchanged with row IPIV(i).  IPIV(i) will always be either
   *          i or i+1; IPIV(i) = i indicates a row interchange was not
   *          required.
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
   *          On entry, the matrix of right hand side vectors B.
   *          On exit, B is overwritten by the solution vectors X.
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
    LAPACK_F77NAME(sgttrs)(
      character const * TRANS,
      integer   const * N,
      integer   const * NRHS,
      real      const * DL,
      real      const * D,
      real      const * DU,
      real      const * DU2,
      integer   const * IPIV,
      real            * B,
      integer   const * LDB,
      integer         * INFO
    );

    void
    LAPACK_F77NAME(dgttrs)(
      character  const * TRANS,
      integer    const * N,
      integer    const * NRHS,
      doublereal const * DL,
      doublereal const * D,
      doublereal const * DU,
      doublereal const * DU2,
      integer    const * IPIV,
      doublereal       * B,
      integer    const * LDB,
      integer          * INFO
    );
  }
  #endif

  inline
  integer
  gttrs(
    Transposition const & TRANS,
    integer               N,
    integer               NRHS,
    real          const   DL[],
    real          const   D[],
    real          const   DU[],
    real          const   DU2[],
    integer       const   IPIV[],
    real                  B[],
    integer               LDB
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sgttrs)(
      const_cast<character*>(trans_blas[TRANS]),
      &N, &NRHS, DL, D, DU, DU2, IPIV, B, &LDB, &INFO
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgttrs( trans_blas[TRANS], &N, &NRHS, DL, D, DU, DU2, IPIV, B, &LDB, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgttrs)(
      const_cast<character*>(trans_blas[TRANS]),
      &N, &NRHS,
      const_cast<real*>(DL),
      const_cast<real*>(D),
      const_cast<real*>(DU),
      const_cast<real*>(DU2),
      const_cast<integer*>(IPIV), B, &LDB, &INFO
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  inline
  integer
  gttrs(
    Transposition const & TRANS,
    integer               N,
    integer               NRHS,
    doublereal    const   DL[],
    doublereal    const   D[],
    doublereal    const   DU[],
    doublereal    const   DU2[],
    integer       const   IPIV[],
    doublereal            B[],
    integer               LDB
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dgttrs)(
      const_cast<character*>(trans_blas[TRANS]),
      &N, &NRHS, DL, D, DU, DU2, IPIV, B, &LDB, &INFO
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgttrs( trans_blas[TRANS], &N, &NRHS, DL, D, DU, DU2, IPIV, B, &LDB, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgttrs)(
      const_cast<character*>(trans_blas[TRANS]),
      &N, &NRHS,
      const_cast<doublereal*>(DL),
      const_cast<doublereal*>(D),
      const_cast<doublereal*>(DU),
      const_cast<doublereal*>(DU2),
      const_cast<integer*>(IPIV),
      B, &LDB, &INFO
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  /*\
   *
   *
   *  Purpose
   *  =======
   *
   *  DPTTRF computes the L*D*L' factorization of a real symmetric
   *  positive definite tridiagonal matrix A.  The factorization may also
   *  be regarded as having the form A = U'*D*U.
   *
   *  Arguments
   *  =========
   *
   *  N       (input) INTEGER
   *          The order of the matrix A.  N >= 0.
   *
   *  D       (input/output) DOUBLE PRECISION array, dimension (N)
   *          On entry, the n diagonal elements of the tridiagonal matrix
   *          A.  On exit, the n diagonal elements of the diagonal matrix
   *          D from the L*D*L' factorization of A.
   *
   *  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
   *          On entry, the (n-1) subdiagonal elements of the tridiagonal
   *          matrix A.  On exit, the (n-1) subdiagonal elements of the
   *          unit bidiagonal factor L from the L*D*L' factorization of A.
   *
   *          E can also be regarded as the superdiagonal of the unit
   *          bidiagonal factor U from the U'*D*U factorization of A.
   *
   *  INFO    (output) INTEGER
   *          = 0: successful exit
   *          < 0: if INFO = -k, the k-th argument had an illegal value
   *          > 0: if INFO = k, the leading minor of order k is not
   *               positive definite; if k < N, the factorization could not
   *               be completed, while if k = N, the factorization was
   *               completed, but D(N) = 0.
   *
   *  =====================================================================
   *
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {

    void
    LAPACK_F77NAME(spttrf)(
      integer const * N,
      real          * D,
      real          * E,
      integer       * INFO
    );

    void
    LAPACK_F77NAME(dpttrf)(
      integer const * N,
      doublereal    * D,
      doublereal    * E,
      integer       * INFO
    );
  }
  #endif

  inline
  integer
  pttrf(
    integer N,
    real    D[],
    real    E[]
  ) {
    integer INFO=0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(spttrf)( &N, D, E, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    spttrf( &N, D, E, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(spttrf)( &N, D, E, &INFO );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  inline
  integer
  pttrf(
    integer    N,
    doublereal D[],
    doublereal E[]
  ) {
    integer INFO=0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dpttrf)( &N, D, E, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dpttrf( &N, D, E, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dpttrf)( &N, D, E, &INFO );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  /*\
   *     ..
   *
   *  Purpose
   *  =======
   *
   *  DPTTRS solves a tridiagonal system of the form
   *     A * X = B
   *  using the L*D*L' factorization of A computed by DPTTRF.  D is a
   *  diagonal matrix specified in the vector D, L is a unit bidiagonal
   *  matrix whose subdiagonal is specified in the vector E, and X and B
   *  are N by NRHS matrices.
   *
   *  Arguments
   *  =========
   *
   *  N       (input) INTEGER
   *          The order of the tridiagonal matrix A.  N >= 0.
   *
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of columns
   *          of the matrix B.  NRHS >= 0.
   *
   *  D       (input) DOUBLE PRECISION array, dimension (N)
   *          The n diagonal elements of the diagonal matrix D from the
   *          L*D*L' factorization of A.
   *
   *  E       (input) DOUBLE PRECISION array, dimension (N-1)
   *          The (n-1) subdiagonal elements of the unit bidiagonal factor
   *
   *          L from the L*D*L' factorization of A.  E can also be regarded
   *          as the superdiagonal of the unit bidiagonal factor U from the
   *          factorization A = U'*D*U.
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
   *          On entry, the right hand side vectors B for the system of
   *          linear equations.
   *          On exit, the solution vectors, X.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B.  LDB >= max(1,N).
   *
   *  INFO    (output) INTEGER
   *          = 0: successful exit
   *          < 0: if INFO = -k, the k-th argument had an illegal value
   *
   *  =====================================================================
   *
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {

    void
    LAPACK_F77NAME(spttrs)(
      integer const * N,
      integer const * NRHS,
      real    const * D,
      real    const * E,
      real          * B,
      integer const * LDB,
      integer       * INFO
    );

    void
    LAPACK_F77NAME(dpttrs)(
      integer    const * N,
      integer    const * NRHS,
      doublereal const * D,
      doublereal const * E,
      doublereal       * B,
      integer    const * LDB,
      integer          * INFO
    );
  }
  #endif

  inline
  integer
  pttrs(
    integer    N,
    integer    NRHS,
    real const D[],
    real const E[],
    real       B[],
    integer    LDB
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(spttrs)( &N, &NRHS, D, E, B, &LDB, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    spttrs( &N, &NRHS, D, E, B, &LDB, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(spttrs)(
      &N, &NRHS,
      const_cast<real*>(D),
      const_cast<real*>(E),
      B, &LDB, &INFO
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  inline
  integer
  pttrs(
    integer          N,
    integer          NRHS,
    doublereal const D[],
    doublereal const E[],
    doublereal       B[],
    integer          LDB
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dpttrs)( &N, &NRHS, D, E, B, &LDB, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dpttrs( &N, &NRHS, D, E, B, &LDB, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dpttrs)(
      &N, &NRHS,
      const_cast<doublereal*>(D),
      const_cast<doublereal*>(E),
      B, &LDB, &INFO
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  /*\
   *  Purpose
   *  =======
   *
   *  DGTSV  solves the equation
   *
   *     A*X = B,
   *
   *  where A is an n by n tridiagonal matrix, by Gaussian elimination with
   *  partial pivoting.
   *
   *  Note that the equation  A**T*X = B  may be solved by interchanging the
   *  order of the arguments DU and DL.
   *
   *  Arguments
   *  =========
   *
   *  N       (input) INTEGER
   *          The order of the matrix A.  N >= 0.
   *
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of columns
   *          of the matrix B.  NRHS >= 0.
   *
   *  DL      (input/output) DOUBLE PRECISION array, dimension (N-1)
   *          On entry, DL must contain the (n-1) sub-diagonal elements of A.
   *
   *          On exit, DL is overwritten by the (n-2) elements of the
   *          second super-diagonal of the upper triangular matrix U from
   *          the LU factorization of A, in DL(1), ..., DL(n-2).
   *
   *  D       (input/output) DOUBLE PRECISION array, dimension (N)
   *          On entry, D must contain the diagonal elements of A.
   *
   *          On exit, D is overwritten by the n diagonal elements of U.
   *
   *  DU      (input/output) DOUBLE PRECISION array, dimension (N-1)
   *          On entry, DU must contain the (n-1) super-diagonal elements of A.
   *
   *          On exit, DU is overwritten by the (n-1) elements of the first
   *          super-diagonal of U.
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
   *          On entry, the N by NRHS matrix of right hand side matrix B.
   *          On exit, if INFO = 0, the N by NRHS solution matrix X.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B.  LDB >= max(1,N).
   *
   *  INFO    (output) INTEGER
   *          = 0: successful exit
   *          < 0: if INFO = -i, the i-th argument had an illegal value
   *          > 0: if INFO = i, U(i,i) is exactly zero, and the solution
   *               has not been computed.  The factorization has not been
   *               completed unless i = N.
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {

    void
    LAPACK_F77NAME(sgtsv)(
      integer const * N,
      integer const * NRHS,
      real          * DL,
      real          * D,
      real          * DU,
      real          * B,
      integer const * LDB,
      integer       * INFO
    );

    void
    LAPACK_F77NAME(dgtsv)(
      integer const * N,
      integer const * NRHS,
      doublereal    * DL,
      doublereal    * D,
      doublereal    * DU,
      doublereal    * B,
      integer const * LDB,
      integer       * INFO
    );
  }
  #endif

  inline
  integer
  gtsv(
    integer N,
    integer NRHS,
    real    DL[],
    real    D[],
    real    DU[],
    real    B[],
    integer LDB
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sgtsv)( &N, &NRHS, DL, D, DU, B, &LDB, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgtsv( &N, &NRHS, DL, D, DU, B, &LDB, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgtsv)( &N, &NRHS, DL, D, DU, B, &LDB, &INFO );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  inline
  integer
  gtsv(
    integer    N,
    integer    NRHS,
    doublereal DL[],
    doublereal D[],
    doublereal DU[],
    doublereal B[],
    integer    LDB
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dgtsv)( &N, &NRHS, DL, D, DU, B, &LDB, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgtsv( &N, &NRHS, DL, D, DU, B, &LDB, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgtsv)( &N, &NRHS, DL, D, DU, B, &LDB, &INFO );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  /*\
   |    ____                _ _ _   _
   |   / ___|___  _ __   __| (_) |_(_) ___  _ __
   |  | |   / _ \| '_ \ / _` | | __| |/ _ \| '_ \
   |  | |__| (_) | | | | (_| | | |_| | (_) | | | |
   |   \____\___/|_| |_|\__,_|_|\__|_|\___/|_| |_|
   |   _   _                 _
   |  | \ | |_   _ _ __ ___ | |__   ___ _ __
   |  |  \| | | | | '_ ` _ \| '_ \ / _ \ '__|
   |  | |\  | |_| | | | | | | |_) |  __/ |
   |  |_| \_|\__,_|_| |_| |_|_.__/ \___|_|
  \*/

  /*\
   *  Purpose
   *  =======
   *
   *  DGTCON estimates the reciprocal of the condition number of a real
   *  tridiagonal matrix A using the LU factorization as computed by
   *  DGTTRF.
   *
   *  An estimate is obtained for norm(inv(A)), and the reciprocal of the
   *  condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
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
   *  DL      (input) DOUBLE PRECISION array, dimension (N-1)
   *          The (n-1) multipliers that define the matrix L from the
   *          LU factorization of A as computed by DGTTRF.
   *
   *  D       (input) DOUBLE PRECISION array, dimension (N)
   *          The n diagonal elements of the upper triangular matrix U from
   *          the LU factorization of A.
   *
   *  DU      (input) DOUBLE PRECISION array, dimension (N-1)
   *          The (n-1) elements of the first superdiagonal of U.
   *
   *  DU2     (input) DOUBLE PRECISION array, dimension (N-2)
   *          The (n-2) elements of the second superdiagonal of U.
   *
   *  IPIV    (input) INTEGER array, dimension (N)
   *          The pivot indices; for 1 <= i <= n, row i of the matrix was
   *          interchanged with row IPIV(i).  IPIV(i) will always be either
   *          i or i+1; IPIV(i) = i indicates a row interchange was not
   *          required.
   *
   *  ANORM   (input) DOUBLE PRECISION
   *          If NORM = '1' or 'O', the 1-norm of the original matrix A.
   *          If NORM = 'I', the infinity-norm of the original matrix A.
   *
   *  RCOND   (output) DOUBLE PRECISION
   *          The reciprocal of the condition number of the matrix A,
   *          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
   *          estimate of the 1-norm of inv(A) computed in this routine.
   *
   *  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)
   *
   *  IWORK   (workspace) INTEGER array, dimension (N)
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
    LAPACK_F77NAME(sgtcon)(
      character const * NORM,
      integer   const * N,
      real      const * DL,
      real      const * D,
      real      const * DU,
      real      const * DU2,
      integer   const * IPIV,
      real      const * ANORM,
      real            * RCOND,
      real            * WORK,
      integer         * IWORK,
      integer         * INFO
    );

    void
    LAPACK_F77NAME(dgtcon)(
      character  const * NORM,
      integer    const * N,
      doublereal const * DL,
      doublereal const * D,
      doublereal const * DU,
      doublereal const * DU2,
      integer    const * IPIV,
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
  gtcon1(
    integer       N,
    real    const DL[],
    real    const D[],
    real    const DU[],
    real    const DU2[],
    integer const IPIV[],
    real          ANORM,
    real        & RCOND,
    real          WORK[],
    integer       IWORK[]
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sgtcon)(
      const_cast<character*>("1"),
      &N, DL, D, DU, DU2, IPIV, &ANORM, &RCOND, WORK, IWORK, &INFO
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgtcon(
      "1", &N, DL, D, DU, DU2, IPIV,
      &ANORM, &RCOND, WORK, IWORK, &INFO
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgtcon)(
      const_cast<character*>("1"), &N,
      const_cast<real*>(DL),
      const_cast<real*>(D),
      const_cast<real*>(DU),
      const_cast<real*>(DU2),
      const_cast<integer*>(IPIV),
      &ANORM, &RCOND, WORK, IWORK, &INFO
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  inline
  integer
  gtcon1(
    integer          N,
    doublereal const DL[],
    doublereal const D[],
    doublereal const DU[],
    doublereal const DU2[],
    integer    const IPIV[],
    doublereal       ANORM,
    doublereal     & RCOND,
    doublereal       WORK[],
    integer          IWORK[]
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dgtcon)(
      const_cast<character*>("1"),
      &N, DL, D, DU, DU2, IPIV, &ANORM, &RCOND, WORK, IWORK, &INFO
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgtcon(
      "1", &N, DL, D, DU, DU2, IPIV,
      &ANORM, &RCOND, WORK, IWORK, &INFO
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgtcon)(
      const_cast<character*>("1"), &N,
      const_cast<doublereal*>(DL),
      const_cast<doublereal*>(D),
      const_cast<doublereal*>(DU),
      const_cast<doublereal*>(DU2),
      const_cast<integer*>(IPIV),
      &ANORM, &RCOND, WORK, IWORK, &INFO
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  inline
  integer
  gtconInf(
    integer       N,
    real    const DL[],
    real    const D[],
    real    const DU[],
    real    const DU2[],
    integer const IPIV[],
    real          ANORM,
    real        & RCOND,
    real          WORK[],
    integer       IWORK[]
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sgtcon)(
      const_cast<character*>("I"),
      &N, DL, D, DU, DU2, IPIV, &ANORM, &RCOND, WORK, IWORK, &INFO
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgtcon(
      "I", &N, DL, D, DU, DU2, IPIV,
      &ANORM, &RCOND, WORK, IWORK, &INFO
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgtcon)(
      const_cast<character*>("I"), &N,
      const_cast<real*>(DL),
      const_cast<real*>(D),
      const_cast<real*>(DU),
      const_cast<real*>(DU2),
      const_cast<integer*>(IPIV),
      &ANORM, &RCOND, WORK, IWORK, &INFO
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  inline
  integer
  gtconInf(
    integer          N,
    doublereal const DL[],
    doublereal const D[],
    doublereal const DU[],
    doublereal const DU2[],
    integer    const IPIV[],
    doublereal       ANORM,
    doublereal     & RCOND,
    doublereal       WORK[],
    integer          IWORK[]
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dgtcon)(
      const_cast<character*>("I"),
      &N, DL, D, DU, DU2, IPIV, &ANORM, &RCOND, WORK, IWORK, &INFO
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgtcon(
      "I", &N, DL, D, DU, DU2, IPIV,
      &ANORM, &RCOND, WORK, IWORK, &INFO
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgtcon)(
      const_cast<character*>("I"), &N,
      const_cast<doublereal*>(DL),
      const_cast<doublereal*>(D),
      const_cast<doublereal*>(DU),
      const_cast<doublereal*>(DU2),
      const_cast<integer*>(IPIV),
      &ANORM, &RCOND, WORK, IWORK, &INFO
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  /*\
   *  Purpose
   *  =======
   *
   *  DPTCON computes the reciprocal of the condition number (in the
   *  1-norm) of a real symmetric positive definite tridiagonal matrix
   *  using the factorization A = L*D*L**T or A = U**T*D*U computed by
   *  DPTTRF.
   *
   *  Norm(inv(A)) is computed by a direct method, and the reciprocal of
   *  the condition number is computed as
   *               RCOND = 1 / (ANORM * norm(inv(A))).
   *
   *  Arguments
   *  =========
   *
   *  N       (input) INTEGER
   *          The order of the matrix A.  N >= 0.
   *
   *  D       (input) DOUBLE PRECISION array, dimension (N)
   *          The n diagonal elements of the diagonal matrix D from the
   *          factorization of A, as computed by DPTTRF.
   *
   *  E       (input) DOUBLE PRECISION array, dimension (N-1)
   *          The (n-1) off-diagonal elements of the unit bidiagonal factor
   *          U or L from the factorization of A,  as computed by DPTTRF.
   *
   *  ANORM   (input) DOUBLE PRECISION
   *          The 1-norm of the original matrix A.
   *
   *  RCOND   (output) DOUBLE PRECISION
   *          The reciprocal of the condition number of the matrix A,
   *          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is the
   *          1-norm of inv(A) computed in this routine.
   *
   *  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *
   *  Further Details
   *  ===============
   *
   *  The method used is described in Nicholas J. Higham, "Efficient
   *  Algorithms for Computing the Condition Number of a Tridiagonal
   *  Matrix", SIAM J. Sci. Stat. Comput., Vol. 7, No. 1, January 1986.
   *
   *  =====================================================================
  \*/


  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {

    void
    LAPACK_F77NAME(sptcon)(
      integer const * N,
      real    const * D,
      real    const * E,
      real    const * ANORM,
      real          * RCOND,
      real          * WORK,
      integer       * INFO
    );

    void
    LAPACK_F77NAME(dptcon)(
      integer    const * N,
      doublereal const * D,
      doublereal const * E,
      doublereal const * ANORM,
      doublereal       * RCOND,
      doublereal       * WORK,
      integer          * INFO
    );
  }
  #endif

  inline
  integer
  ptcon1(
    integer    N,
    real const D[],
    real const E[],
    real       ANORM,
    real     & RCOND,
    real       WORK[]
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sptcon)( &N, D, E, &ANORM, &RCOND, WORK, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sptcon( &N, D, E, &ANORM, &RCOND, WORK, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sptcon)(
      &N,
      const_cast<real*>(D),
      const_cast<real*>(E),
      &ANORM, &RCOND, WORK, &INFO
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

  inline
  integer
  ptcon1(
    integer          N,
    doublereal const D[],
    doublereal const E[],
    doublereal       ANORM,
    doublereal     & RCOND,
    doublereal       WORK[]
  ) {
    integer INFO = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dptcon)( &N, D, E, &ANORM, &RCOND, WORK, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dptcon( &N, D, E, &ANORM, &RCOND, WORK, &INFO );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dptcon)(
      &N,
      const_cast<doublereal*>(D),
      const_cast<doublereal*>(E),
      &ANORM, &RCOND, WORK, &INFO
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return INFO;
  }

}

///
/// eof: tridiagonal.hxx
///
