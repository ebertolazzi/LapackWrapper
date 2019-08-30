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
/// file: general_qr.hxx
///

namespace lapack_wrapper {

  /*
  //   _            _
  //  | | __ _ _ __| |_ __ _
  //  | |/ _` | '__| __/ _` |
  //  | | (_| | |  | || (_| |
  //  |_|\__,_|_|   \__\__, |
  //                   |___/
  */
  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {
    void
    LAPACK_F77NAME(slartg)(
      real * F,
      real * G,
      real * C,
      real * S,
      real * R
    );

    void
    LAPACK_F77NAME(dlartg)(
      doublereal * F,
      doublereal * G,
      doublereal * C,
      doublereal * S,
      doublereal * R
    );
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DLARTG generate a plane rotation so that
   *
   *     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
   *     [ -SN  CS  ]     [ G ]     [ 0 ]
   *
   *  This is a slower, more accurate version of the BLAS1 routine DROTG,
   *  with the following other differences:
   *     F and G are unchanged on return.
   *     If G=0, then CS=1 and SN=0.
   *     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
   *        floating point operations (saves work in DBDSQR when
   *        there are zeros on the diagonal).
   *
   *  If F exceeds G in magnitude, CS will be positive.
   *
   *  Arguments
   *  =========
   *
   *  F       (input) DOUBLE PRECISION
   *          The first component of vector to be rotated.
   *
   *  G       (input) DOUBLE PRECISION
   *          The second component of vector to be rotated.
   *
   *  CS      (output) DOUBLE PRECISION
   *          The cosine of the rotation.
   *
   *  SN      (output) DOUBLE PRECISION
   *          The sine of the rotation.
   *
   *  R       (output) DOUBLE PRECISION
   *          The nonzero component of the rotated vector.
   *
  \*/

  inline
  void
  lartg(
    real   F,
    real   G,
    real & C,
    real & S,
    real & R
  )
  #ifdef LAPACK_WRAPPER_USE_ACCELERATE
  { CLAPACKNAME(slartg)( &F, &G, &C, &S, &R ); }
  #elif defined(LAPACK_WRAPPER_USE_LAPACK) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
  { LAPACK_F77NAME(slartg)( &F, &G, &C, &S, &R ); }
  #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { LAPACK_F77NAME(slartgp)( &F, &G, &C, &S, &R ); }
  #elif defined(LAPACK_WRAPPER_USE_MKL)
  { slartgp( &F, &G, &C, &S, &R ); }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  inline
  void
  lartg(
    doublereal   F,
    doublereal   G,
    doublereal & C,
    doublereal & S,
    doublereal & R
  )
  #ifdef LAPACK_WRAPPER_USE_ACCELERATE
  { CLAPACKNAME(dlartg)( &F, &G, &C, &S, &R ); }
  #elif defined(LAPACK_WRAPPER_USE_LAPACK) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
  { LAPACK_F77NAME(dlartg)( &F, &G, &C, &S, &R ); }
  #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { LAPACK_F77NAME(dlartgp)( &F, &G, &C, &S, &R ); }
  #elif defined(LAPACK_WRAPPER_USE_MKL)
  { dlartgp( &F, &G, &C, &S, &R ); }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  /*
  //  _
  // | | __ _ _ __  _   _
  // | |/ _` | '_ \| | | |
  // | | (_| | |_) | |_| |
  // |_|\__,_| .__/ \__, |
  //         |_|    |___/
  */

  /*\
   *  Purpose
   *  =======
   *
   *  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
   *  overflow.
   *
   *  Arguments
   *  =========
   *
   *  X       (input) DOUBLE PRECISION
   *  Y       (input) DOUBLE PRECISION
   *          X and Y specify the values x and y.
   *  =====================================================================
   *
   *  DLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause
   *  unnecessary overflow.
   *
   *  Arguments
   *  =========
   *
   *  X       (input) DOUBLE PRECISION
   *  Y       (input) DOUBLE PRECISION
   *  Z       (input) DOUBLE PRECISION
   *          X, Y and Z specify the values x, y and z.
   *  =====================================================================
   *
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {
    real
    LAPACK_F77NAME(slapy2)(
      real const * X,
      real const * Y
    );

    doublereal
    LAPACK_F77NAME(dlapy2)(
      doublereal const * X,
      doublereal const * Y
    );

    real
    LAPACK_F77NAME(slapy3)(
      real const * X,
      real const * Y,
      real const * Z
    );

    doublereal
    LAPACK_F77NAME(dlapy3)(
      doublereal const * X,
      doublereal const * Y,
      doublereal const * Z
    );
  }
  #endif

  inline
  real
  hypot( real x, real y )
  #if defined(LAPACK_WRAPPER_USE_ACCELERATE)
  { return real(CLAPACKNAME(slapy2)( &x, &y )); }
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
  { return LAPACK_F77NAME(slapy2)( &x, &y ); }
  #elif defined(LAPACK_WRAPPER_USE_MKL)
  { return slapy2( &x, &y ); }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  inline
  real
  hypot( real x, real y, real z )
  #if defined(LAPACK_WRAPPER_USE_ACCELERATE)
  { return real(CLAPACKNAME(slapy3)( &x, &y, &z )); }
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
  { return LAPACK_F77NAME(slapy3)( &x, &y, &z ); }
  #elif defined(LAPACK_WRAPPER_USE_MKL)
  { return slapy3( &x, &y, &z ); }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  inline
  doublereal
  hypot( doublereal x, doublereal y )
  #if defined(LAPACK_WRAPPER_USE_ACCELERATE)
  { return CLAPACKNAME(dlapy2)( &x, &y ); }
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
  { return LAPACK_F77NAME(dlapy2)( &x, &y ); }
  #elif defined(LAPACK_WRAPPER_USE_MKL)
  { return dlapy2( &x, &y ); }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  inline
  doublereal
  hypot( doublereal x, doublereal y, doublereal z )
  #if defined(LAPACK_WRAPPER_USE_ACCELERATE)
  { return CLAPACKNAME(dlapy3)( &x, &y, &z ); }
  #elif defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
  { return LAPACK_F77NAME(dlapy3)( &x, &y, &z ); }
  #elif defined(LAPACK_WRAPPER_USE_MKL)
  { return dlapy3( &x, &y, &z ); }
  #else
  #error "LapackWrapper undefined mapping!"
  #endif

  /*
  //    ___  ____
  //   / _ \|  _ \
  //  | | | | |_) |
  //  | |_| |  _ <
  //   \__\_\_| \_\
  */

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  // use standard Lapack routine
  extern "C" {

    void
    LAPACK_F77NAME(sormqr)(
      character const * SIDE,
      character const * TRANS,
      integer   const * M,
      integer   const * N,
      integer   const * K,
      real      const   A[],
      integer   const * LDA,
      real      const   TAU[],
      real              C[],
      integer   const * LDC,
      real              WORK[],
      integer   const * LWORK,
      integer         * INFO
    );

    void
    LAPACK_F77NAME(dormqr)(
      character  const * SIDE,
      character  const * TRANS,
      integer    const * M,
      integer    const * N,
      integer    const * K,
      doublereal const   A[],
      integer    const * LDA,
      doublereal const   TAU[],
      doublereal         C[],
      integer    const * LDC,
      doublereal         WORK[],
      integer    const * LWORK,
      integer          * INFO );
    }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DORMQR overwrites the general real M-by-N matrix C with
   *
   *                  SIDE = 'L'     SIDE = 'R'
   *  TRANS = 'N':      Q * C          C * Q
   *  TRANS = 'T':      Q**T * C       C * Q**T
   *
   *  where Q is a real orthogonal matrix defined as the product of k
   *  elementary reflectors
   *
   *        Q = H(1) H(2) . . . H(k)
   *
   *  as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N
   *  if SIDE = 'R'.
   *
   *  Arguments
   *  =========
   *
   *  SIDE    (input) CHARACTER*1
   *          = 'L': apply Q or Q**T from the Left;
   *          = 'R': apply Q or Q**T from the Right.
   *
   *  TRANS   (input) CHARACTER*1
   *          = 'N':  No transpose, apply Q;
   *          = 'T':  Transpose, apply Q**T.
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix C. M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix C. N >= 0.
   *
   *  K       (input) INTEGER
   *          The number of elementary reflectors whose product defines
   *          the matrix Q.
   *          If SIDE = 'L', M >= K >= 0;
   *          if SIDE = 'R', N >= K >= 0.
   *
   *  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
   *          The i-th column must contain the vector which defines the
   *          elementary reflector H(i), for i = 1,2,...,k, as returned by
   *          DGEQRF in the first k columns of its array argument A.
   *          A is modified by the routine but restored on exit.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.
   *          If SIDE = 'L', LDA >= max(1,M);
   *          if SIDE = 'R', LDA >= max(1,N).
   *
   *  TAU     (input) DOUBLE PRECISION array, dimension (K)
   *          TAU(i) must contain the scalar factor of the elementary
   *          reflector H(i), as returned by DGEQRF.
   *
   *  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
   *          On entry, the M-by-N matrix C.
   *          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
   *
   *  LDC     (input) INTEGER
   *          The leading dimension of the array C. LDC >= max(1,M).
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK.
   *          If SIDE = 'L', LWORK >= max(1,N);
   *          if SIDE = 'R', LWORK >= max(1,M).
   *          For optimum performance LWORK >= N*NB if SIDE = 'L', and
   *          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
   *          blocksize.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *
  \*/

  inline
  integer
  ormqr(
    SideMultiply  const & SIDE,
    Transposition const & TRANS,
    integer               M,
    integer               N,
    integer               K,
    real const            A[],
    integer               LDA,
    real const            TAU[],
    real                  C[],
    integer               LDC,
    real                  WORK[],
    integer               LWORK
  )
  { integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sormqr)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(trans_blas[TRANS]),
      &M, &N, &K, A, &LDA,
      TAU, C, &LDC, WORK, &LWORK, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sormqr(
      side_blas[SIDE], trans_blas[TRANS],
      &M, &N, &K, A, &LDA,
      TAU, C, &LDC, WORK, &LWORK, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sormqr)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(trans_blas[TRANS]),
      &M, &N, &K, const_cast<real*>(A), &LDA,
      const_cast<real*>(TAU),
      C, &LDC, WORK, &LWORK, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  inline
  integer
  ormqr(
    SideMultiply  const & SIDE,
    Transposition const & TRANS,
    integer               M,
    integer               N,
    integer               K,
    doublereal const      A[],
    integer               LDA,
    doublereal const      TAU[],
    doublereal            C[],
    integer               LDC,
    doublereal            WORK[],
    integer               LWORK
  )
  { integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dormqr)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(trans_blas[TRANS]),
      &M, &N, &K, A, &LDA,
      TAU, C, &LDC,
      WORK, &LWORK, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dormqr(
      side_blas[SIDE], trans_blas[TRANS],
      &M, &N, &K, A, &LDA, TAU, C, &LDC,
      WORK, &LWORK, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dormqr)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(trans_blas[TRANS]),
      &M, &N, &K, const_cast<doublereal*>(A), &LDA,
      const_cast<doublereal*>(TAU), C, &LDC,
      WORK, &LWORK, &info
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
   *  DORMQR overwrites the general real M-by-N matrix C with
   *
   *                  SIDE = 'L'     SIDE = 'R'
   *  TRANS = 'N':      Q * C          C * Q
   *  TRANS = 'T':      Q**T * C       C * Q**T
   *
   *  where Q is a real orthogonal matrix defined as the product of k
   *  elementary reflectors
   *
   *        Q = H(1) H(2) . . . H(k)
   *
   *  as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N
   *  if SIDE = 'R'.
   *
   *  Arguments
   *  =========
   *
   *  SIDE    (input) CHARACTER*1
   *          = 'L': apply Q or Q**T from the Left;
   *          = 'R': apply Q or Q**T from the Right.
   *
   *  TRANS   (input) CHARACTER*1
   *          = 'N':  No transpose, apply Q;
   *          = 'T':  Transpose, apply Q**T.
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix C. M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix C. N >= 0.
   *
   *  K       (input) INTEGER
   *          The number of elementary reflectors whose product defines
   *          the matrix Q.
   *          If SIDE = 'L', M >= K >= 0;
   *          if SIDE = 'R', N >= K >= 0.
   *
   *  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
   *          The i-th column must contain the vector which defines the
   *          elementary reflector H(i), for i = 1,2,...,k, as returned by
   *          DGEQRF in the first k columns of its array argument A.
   *          A is modified by the routine but restored on exit.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.
   *          If SIDE = 'L', LDA >= max(1,M);
   *          if SIDE = 'R', LDA >= max(1,N).
   *
   *  TAU     (input) DOUBLE PRECISION array, dimension (K)
   *          TAU(i) must contain the scalar factor of the elementary
   *          reflector H(i), as returned by DGEQRF.
   *
   *  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
   *          On entry, the M-by-N matrix C.
   *          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
   *
   *  LDC     (input) INTEGER
   *          The leading dimension of the array C. LDC >= max(1,M).
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK.
   *          If SIDE = 'L', LWORK >= max(1,N);
   *          if SIDE = 'R', LWORK >= max(1,M).
   *          For optimum performance LWORK >= N*NB if SIDE = 'L', and
   *          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
   *          blocksize.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *
   *  =====================================================================
  \*/

  inline
  integer
  ormqr(
    SideMultiply  const & SIDE,
    Transposition const & TRANS,
    integer               M,
    integer               N,
    integer               K,
    real                  A[],
    integer               LDA,
    real                  TAU[],
    real                  C[],
    integer               LDC,
    real                  WORK[],
    integer               LWORK
  )
  { integer info = 0;
    #ifdef LAPACK_WRAPPER_USE_ACCELERATE
    CLAPACKNAME(sormqr)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(trans_blas[TRANS]),
      &M, &N, &K, A, &LDA, TAU, C, &LDC,
      WORK, &LWORK, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_LAPACK)   || \
          defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
          defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sormqr)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(trans_blas[TRANS]),
      &M, &N, &K, A, &LDA, TAU, C, &LDC,
      WORK, &LWORK, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sormqr(
      side_blas[SIDE], trans_blas[TRANS],
      &M, &N, &K, A, &LDA, TAU, C, &LDC,
      WORK, &LWORK, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline
  integer
  ormqr(
    SideMultiply  const & SIDE,
    Transposition const & TRANS,
    integer               M,
    integer               N,
    integer               K,
    doublereal            A[],
    integer               LDA,
    doublereal            TAU[],
    doublereal            C[],
    integer               LDC,
    doublereal            WORK[],
    integer               LWORK
  )
  { integer info = 0;
    #ifdef LAPACK_WRAPPER_USE_ACCELERATE
    CLAPACKNAME(dormqr)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(trans_blas[TRANS]),
      &M, &N, &K, A, &LDA, TAU, C, &LDC,
      WORK, &LWORK, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
          defined(LAPACK_WRAPPER_USE_LAPACK)   || \
          defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dormqr)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(trans_blas[TRANS]),
      &M, &N, &K, A, &LDA, TAU, C, &LDC,
      WORK, &LWORK, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dormqr(
      side_blas[SIDE], trans_blas[TRANS],
      &M, &N, &K, A, &LDA, TAU, C, &LDC,
      WORK, &LWORK, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  /*
  //   _             __ _
  //  | | __ _ _ __ / _| |_
  //  | |/ _` | '__| |_| __|
  //  | | (_| | |  |  _| |_
  //  |_|\__,_|_|  |_|  \__|
  */
  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  // use standard Lapack routine
  extern "C" {

    void
    LAPACK_F77NAME(slarft)(
      character const   DIRECT[],
      character const   STOREV[],
      integer   const * N,
      integer   const * K,
      real              V[],
      integer   const * LDV,
      real      const   TAU[],
      real              T[],
      integer   const * LDT
    );

    void
    LAPACK_F77NAME(dlarft)(
      character  const   DIRECT[],
      character  const   STOREV[],
      integer    const * N,
      integer    const * K,
      doublereal         V[],
      integer    const * LDV,
      doublereal const   TAU[],
      doublereal         T[],
      integer    const * LDT
    );
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  forms the triangular factor T of a real block reflector H
   *  of order n, which is defined as a product of k elementary reflectors.
   *
   *  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
   *
   *  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
   *
   *  If STOREV = 'C', the vector which defines the elementary reflector
   *  H(i) is stored in the i-th column of the array V, and
   *
   *     H  =  I - V * T * V'
   *
   *  If STOREV = 'R', the vector which defines the elementary reflector
   *  H(i) is stored in the i-th row of the array V, and
   *
   *     H  =  I - V' * T * V
   *
   *  Arguments
   *  =========
   *
   *  DIRECT = Specifies the order in which the elementary reflectors are
   *           multiplied to form the block reflector:
   *         = 'F': H = H(1) H(2) . . . H(k) (Forward)
   *         = 'B': H = H(k) . . . H(2) H(1) (Backward)
   *
   *  STOREV = Specifies how the vectors which define the elementary
   *           reflectors are stored (see also Further Details):
   *         = 'C': columnwise
   *         = 'R': rowwise
   *
   *  N = The order of the block reflector H. N >= 0.
   *  K = The order of the triangular factor T (= the number of elementary reflectors). K >= 1.
   *  V = (LDV,K) if STOREV = 'C'
   *      (LDV,N) if STOREV = 'R'
   *      The matrix V. See further details.
   *
   *  LDV = The leading dimension of the array V.
   *        If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
   *
   *  TAU = TAU(i) must contain the scalar factor of the elementary reflector H(i).
   *
   *  T = The k by k triangular factor T of the block reflector.
   *      If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
   *      lower triangular. The rest of the array is not used.
   *
   *  LDT = The leading dimension of the array T. LDT >= K.
   *
   *  Further Details
   *  ===============
   *
   *  The shape of the matrix V and the storage of the vectors which define
   *  the H(i) is best illustrated by the following example with n = 5 and
   *  k = 3. The elements equal to 1 are not stored; the corresponding
   *  array elements are modified but restored on exit. The rest of the
   *  array is not used.
   *
   *  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
   *
   *               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
   *                   ( v1  1    )                     (     1 v2 v2 v2 )
   *                   ( v1 v2  1 )                     (        1 v3 v3 )
   *                   ( v1 v2 v3 )
   *                   ( v1 v2 v3 )
   *
   *  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
   *
   *               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
   *                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
   *                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
   *                   (     1 v3 )
   *                   (        1 )
   *
  \*/

  inline
  void
  larft(
    DirectionType & DIRECT, // direct_blas[DIRECT]
    StorageType   & STOREV, // store_blas[STOREV]
    integer         N,
    integer         K,
    real            V[],
    integer         LDV,
    real            TAU[],
    real            T[],
    integer         LDT
  ) {
  #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
      defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(slarft)(
      const_cast<character*>(direct_blas[DIRECT]),
      const_cast<character*>(store_blas[STOREV]),
      &N, &K, V, &LDV, TAU, T, &LDT
    );
  #elif defined(LAPACK_WRAPPER_USE_MKL)
    slarft(
      direct_blas[DIRECT], store_blas[STOREV],
      &N, &K, V, &LDV, TAU, T, &LDT
    );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(slarft)(
      const_cast<character*>(direct_blas[DIRECT]),
      const_cast<character*>(store_blas[STOREV]),
      &N, &K, V, &LDV, TAU, T, &LDT
    );
  #else
  #error "LapackWrapper undefined mapping!"
  #endif
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline
  void
  larft(
    DirectionType & DIRECT, // direct_blas[DIRECT]
    StorageType   & STOREV, // store_blas[STOREV]
    integer         N,
    integer         K,
    doublereal      V[],
    integer         LDV,
    doublereal      TAU[],
    doublereal      T[],
    integer         LDT
  ) {
  #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
      defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dlarft)(
      const_cast<character*>(direct_blas[DIRECT]),
      const_cast<character*>(store_blas[STOREV]),
      &N, &K, V, &LDV, TAU, T, &LDT
    );
  #elif defined(LAPACK_WRAPPER_USE_MKL)
    dlarft(
      direct_blas[DIRECT], store_blas[STOREV],
      &N, &K, V, &LDV, TAU, T, &LDT
    );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dlarft)(
      const_cast<character*>(direct_blas[DIRECT]),
      const_cast<character*>(store_blas[STOREV]),
      &N, &K, V, &LDV, TAU, T, &LDT
    );
  #else
  #error "LapackWrapper undefined mapping!"
  #endif
  }

  /*
  //   _             __
  //  | | __ _ _ __ / _| __ _
  //  | |/ _` | '__| |_ / _` |
  //  | | (_| | |  |  _| (_| |
  //  |_|\__,_|_|  |_|  \__, |
  //                    |___/
  */
  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  // use standard Lapack routine
  extern "C" {

    void
    LAPACK_F77NAME(slarfg)(
      integer const * N,
      real          * ALPHA,
      real            X[],
      integer const * INCX,
      real            TAU[]
    );

    void
    LAPACK_F77NAME(dlarfg)(
      integer const * N,
      doublereal    * ALPHA,
      doublereal      X[],
      integer const * INCX,
      doublereal      TAU[]
    );
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  generates a real elementary reflector H of order n, such that
   *
   *        H * ( alpha ) = ( beta ),   H' * H = I.
   *            (   x   )   (   0  )
   *
   *  where alpha and beta are scalars, and x is an (n-1)-element real
   *  vector. H is represented in the form
   *
   *        H = I - tau * ( 1 ) * ( 1 v' ) ,
   *                      ( v )
   *
   *  where tau is a real scalar and v is a real (n-1)-element vector.
   *
   *  If the elements of x are all zero, then tau = 0 and H is taken to be
   *  the unit matrix.
   *
   *  Otherwise  1 <= tau <= 2.
   *
   *  Arguments
   *  =========
   *
   *  N     = The order of the elementary reflector.
   *  ALPHA = On entry, the value alpha.
   *          On exit, it is overwritten with the value beta.
   *  X     = On entry, the vector x.
   *          On exit, it is overwritten with the vector v.
   *  INCX = The increment between elements of X. INCX <> 0.
   *  TAU  = The value tau.
  \*/

  inline
  void
  larfg(
    integer N,
    real  & ALPHA,
    real    X[],
    integer INCX,
    real    TAU[]
  ) {
  #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
      defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(slarfg)( &N, &ALPHA, X, &INCX, TAU );
  #elif defined(LAPACK_WRAPPER_USE_MKL)
    slarfg( &N, &ALPHA, X, &INCX, TAU );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(slarfg)( &N, &ALPHA, X, &INCX, TAU );
  #else
  #error "LapackWrapper undefined mapping!"
  #endif
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline
  void
  larfg(
    integer      N,
    doublereal & ALPHA,
    doublereal   X[],
    integer      INCX,
    doublereal   TAU[]
  ) {
  #if defined(LAPACK_WRAPPER_USE_LAPACK)  || \
     defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
     defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dlarfg)( &N, &ALPHA, X, &INCX, TAU );
  #elif defined(LAPACK_WRAPPER_USE_MKL)
    dlarfg( &N, &ALPHA, X, &INCX, TAU );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dlarfg)( &N, &ALPHA, X, &INCX, TAU );
  #else
  #error "LapackWrapper undefined mapping!"
  #endif
  }

  /*
  //   _             __ _
  //  | | __ _ _ __ / _| |__
  //  | |/ _` | '__| |_| '_ \
  //  | | (_| | |  |  _| |_) |
  //  |_|\__,_|_|  |_| |_.__/
  */

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  // use standard Lapack routine
  extern "C" {

    void
    LAPACK_F77NAME(slarfb)(
      character const   SIDE[],
      character const   TRANS[],
      character const   DIRECT[],
      character const   STOREV[],
      integer   const * M,
      integer   const * N,
      integer   const * K,
      real      const   V[],
      integer   const * LDV,
      real      const   T[],
      integer   const * LDT,
      real              C[],
      integer   const * LDC,
      real              WORK[],
      integer   const * LDWORK
    );

    void
    LAPACK_F77NAME(dlarfb)(
      character  const   SIDE[],
      character  const   TRANS[],
      character  const   DIRECT[],
      character  const   STOREV[],
      integer    const * M,
      integer    const * N,
      integer    const * K,
      doublereal const   V[],
      integer    const * LDV,
      doublereal const   T[],
      integer    const * LDT,
      doublereal         C[],
      integer    const * LDC,
      doublereal         WORK[],
      integer    const * LDWORK
    );
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  applies a real block reflector H or its transpose H' to a
   *  real m by n matrix C, from either the left or the right.
   *
   *  Arguments
   *  =========
   *
   *  SIDE = 'L': apply H or H' from the Left
   *       = 'R': apply H or H' from the Right
   *
   *  TRANS = 'N': apply H (No transpose)
   *        = 'T': apply H' (Transpose)
   *
   *  DIRECT Indicates how H is formed from a product of elementary reflectors
   *       = 'F': H = H(1) H(2) . . . H(k) (Forward)
   *       = 'B': H = H(k) . . . H(2) H(1) (Backward)
   *
   *  STOREV Indicates how the vectors which define the elementary reflectors are stored:
   *       = 'C': Columnwise
   *       = 'R': Rowwise
   *
   *  M = The number of rows of the matrix C.
   *  N = The number of columns of the matrix C.
   *  K = The order of the matrix T (= the number of elementary
   *      reflectors whose product defines the block reflector).
   *  V = (LDV,K) if STOREV = 'C'
   *      (LDV,M) if STOREV = 'R' and SIDE = 'L'
   *      (LDV,N) if STOREV = 'R' and SIDE = 'R'
   *      The matrix V. See further details.
   *
   *  LDV = The leading dimension of the array V.
   *          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
   *          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
   *          if STOREV = 'R', LDV >= K.
   *
   *  T = (LDT,K) The triangular k by k matrix T in the representation of the block reflector.
   *  LDT = The leading dimension of the array T. LDT >= K.
   *  C  = (LDC,N)
   *     On entry, the m by n matrix C.
   *     On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
   *  LDC = The leading dimension of the array C. LDA >= max(1,M).
   *  WORK = (workspace) DOUBLE PRECISION array, dimension (LDWORK,K)
   *  LDWORK = The leading dimension of the array WORK.
   *     If SIDE = 'L', LDWORK >= max(1,N);
   *     if SIDE = 'R', LDWORK >= max(1,M).
   *
   *  =====================================================================
  \*/

  inline
  void
  larfb(
    SideMultiply  const & SIDE,
    Transposition const & TRANS,
    DirectionType const & DIRECT,
    StorageType   const & STOREV,
    integer               M,
    integer               N,
    integer               K,
    real          const   V[],
    integer               LDV,
    real          const   T[],
    integer               LDT,
    real                  C[],
    integer               LDC,
    real                  WORK[],
    integer               LDWORK
  ) {
  #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
      defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(slarfb)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(trans_blas[TRANS]),
      const_cast<character*>(direct_blas[DIRECT]),
      const_cast<character*>(store_blas[STOREV]),
      &M, &N, &K, V, &LDV,
      T, &LDT, C, &LDC, WORK, &LDWORK
    );
  #elif defined(LAPACK_WRAPPER_USE_MKL)
    slarfb(
      side_blas[SIDE], trans_blas[TRANS],
      direct_blas[DIRECT], store_blas[STOREV],
      &M, &N, &K, V, &LDV, T, &LDT, C, &LDC, WORK, &LDWORK
    );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(slarfb)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(trans_blas[TRANS]),
      const_cast<character*>(direct_blas[DIRECT]),
      const_cast<character*>(store_blas[STOREV]),
      &M, &N, &K,
      const_cast<real*>(V), &LDV,
      const_cast<real*>(T), &LDT,
      C, &LDC, WORK, &LDWORK
    );
  #else
  #error "LapackWrapper undefined mapping!"
  #endif
  }

  inline
  void
  larfb(
    SideMultiply  const & SIDE,
    Transposition const & TRANS,
    DirectionType const & DIRECT,
    StorageType   const & STOREV,
    integer               M,
    integer               N,
    integer               K,
    doublereal    const   V[],
    integer               LDV,
    doublereal    const   T[],
    integer               LDT,
    doublereal            C[],
    integer               LDC,
    doublereal            WORK[],
    integer               LDWORK
  ) {
  #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
      defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dlarfb)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(trans_blas[TRANS]),
      const_cast<character*>(direct_blas[DIRECT]),
      const_cast<character*>(store_blas[STOREV]),
      &M, &N, &K, V, &LDV,
      T, &LDT, C, &LDC, WORK, &LDWORK
    );
  #elif defined(LAPACK_WRAPPER_USE_MKL)
    dlarfb(
      side_blas[SIDE], trans_blas[TRANS],
      direct_blas[DIRECT], store_blas[STOREV],
      &M, &N, &K, V, &LDV, T, &LDT, C, &LDC, WORK, &LDWORK
    );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dlarfb)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(trans_blas[TRANS]),
      const_cast<character*>(direct_blas[DIRECT]),
      const_cast<character*>(store_blas[STOREV]),
      &M, &N, &K,
      const_cast<doublereal*>(V), &LDV,
      const_cast<doublereal*>(T), &LDT,
      C, &LDC, WORK, &LDWORK
    );
  #else
  #error "LapackWrapper undefined mapping!"
  #endif
  }

  /*
  //                          __
  //    __ _  ___  __ _ _ __ / _|
  //   / _` |/ _ \/ _` | '__| |_
  //  | (_| |  __/ (_| | |  |  _|
  //   \__, |\___|\__, |_|  |_|
  //   |___/         |_|
  */

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {

    void
    LAPACK_F77NAME(sgeqrf)(
      integer const * M,
      integer const * N,
      real            A[],
      integer const * LDA,
      real            TAU[],
      real            WORK[],
      integer const * LWORK,
      integer       * INFO
    );

    void
    LAPACK_F77NAME(dgeqrf)(
      integer    const * M,
      integer    const * N,
      doublereal         A[],
      integer    const * LDA,
      doublereal         TAU[],
      doublereal         WORK[],
      integer    const * LWORK,
      integer          * INFO
    );
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DGEQRF computes a QR factorization of a real M-by-N matrix A:
   *  A = Q * R.
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
   *          On entry, the M-by-N matrix A.
   *          On exit, the elements on and above the diagonal of the array
   *          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
   *          upper triangular if m >= n); the elements below the diagonal,
   *          with the array TAU, represent the orthogonal matrix Q as a
   *          product of min(m,n) elementary reflectors (see Further
   *          Details).
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
   *          The scalar factors of the elementary reflectors (see Further
   *          Details).
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK.  LWORK >= max(1,N).
   *          For optimum performance LWORK >= N*NB, where NB is
   *          the optimal blocksize.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *
   *  Further Details
   *  ===============
   *
   *  The matrix Q is represented as a product of elementary reflectors
   *
   *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
   *
   *  Each H(i) has the form
   *
   *     H(i) = I - tau * v * v'
   *
   *  where tau is a real scalar, and v is a real vector with
   *  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
   *  and tau in TAU(i).
   *
   *  =====================================================================
  \*/

  inline
  integer
  geqrf(
    integer M,
    integer N,
    real    A[],
    integer LDA,
    real    TAU[],
    real    WORK[],
    integer LWORK
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sgeqrf)( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info );
    //#elif defined(LAPACK_WRAPPER_USE_ATLAS)
    //if ( LWORK == -1 ) { info = 0; WORK[0] = 1; }
    //else info = CLAPACKNAME(sgeqrf)( CblasColMajor, M, N, A, LDA, TAU );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgeqrf( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgeqrf)( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  inline
  integer
  geqrf(
    integer    M,
    integer    N,
    doublereal A[],
    integer    LDA,
    doublereal TAU[],
    doublereal WORK[],
    integer    LWORK
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dgeqrf)( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info );
    //#elif defined(LAPACK_WRAPPER_USE_ATLAS)
    //if ( LWORK == -1 ) { info = 0; WORK[0] = 1; }
    //else info = CLAPACKNAME(dgeqrf)( CblasColMajor, M, N, A, LDA, TAU );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgeqrf( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgeqrf)( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  /*
  //                        ____
  //    __ _  ___  __ _ _ _|___ \
  //   / _` |/ _ \/ _` | '__|__) |
  //  | (_| |  __/ (_| | |  / __/
  //   \__, |\___|\__, |_| |_____|
  //   |___/         |_|
  */

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  extern "C" {

    void
    LAPACK_F77NAME(sgeqr2)(
      integer const * M,
      integer const * N,
      real            A[],
      integer const * LDA,
      real            TAU[],
      real            WORK[],
      integer       * INFO
    );

    void
    LAPACK_F77NAME(dgeqr2)(
      integer    const * M,
      integer    const * N,
      doublereal         A[],
      integer    const * LDA,
      doublereal         TAU[],
      doublereal         WORK[],
      integer          * INFO
    );
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  computes a QR factorization of a real m by n matrix A:
   *  A = Q * R.
   *
   *  Arguments
   *  =========
   *
   *  M = The number of rows of the matrix A.  M >= 0.
   *  N = The number of columns of the matrix A.  N >= 0.
   *  A = (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the m by n matrix A.
   *          On exit, the elements on and above the diagonal of the array
   *          contain the min(m,n) by n upper trapezoidal matrix R (R is
   *          upper triangular if m >= n); the elements below the diagonal,
   *          with the array TAU, represent the orthogonal matrix Q as a
   *          product of elementary reflectors (see Further Details).
   *
   *  LDA = The leading dimension of the array A.  LDA >= max(1,M).
   *  TAU = The scalar factors of the elementary reflectors (see Further Details).
   *  WORK = (workspace) DOUBLE PRECISION array, dimension (N)
   *  INFO = 0: successful exit
   *       < 0: if INFO = -i, the i-th argument had an illegal value
   *
   *  Further Details
   *  ===============
   *  The matrix Q is represented as a product of elementary reflectors
   *
   *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
   *
   *  Each H(i) has the form
   *
   *     H(i) = I - tau * v * v'
   *
   *  where tau is a real scalar, and v is a real vector with
   *  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
   *  and tau in TAU(i).
  \*/

  inline
  integer
  geqr2(
    integer M,
    integer N,
    real    A[],
    integer LDA,
    real    TAU[],
    real    WORK[]
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sgeqr2)( &M, &N, A, &LDA, TAU, WORK, &info );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgeqr2( &M, &N, A, &LDA, TAU, WORK, &info );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgeqr2)( &M, &N, A, &LDA, TAU, WORK, &info );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  inline
  integer
  geqr2(
    integer    M,
    integer    N,
    doublereal A[],
    integer    LDA,
    doublereal TAU[],
    doublereal WORK[]
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dgeqr2)( &M, &N, A, &LDA, TAU, WORK, &info );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgeqr2( &M, &N, A, &LDA, TAU, WORK, &info );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgeqr2)( &M, &N, A, &LDA, TAU, WORK, &info );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  /*
  //                         _____
  //    __ _  ___  __ _ _ __|___ /
  //   / _` |/ _ \/ _` | '_ \ |_ \
  //  | (_| |  __/ (_| | |_) |__) |
  //   \__, |\___|\__, | .__/____/
  //   |___/         |_|_|
  */

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  // use standard Lapack routine
  extern "C" {

    void
    LAPACK_F77NAME(sgeqp3)(
      integer const * M,
      integer const * N,
      real            A[],
      integer const * LDA,
      integer         JPVT[],
      real            TAU[],
      real            WORK[],
      integer const * LWORK,
      integer       * INFO
    );

    void
    LAPACK_F77NAME(dgeqp3)(
      integer   const * M,
      integer   const * N,
      doublereal        A[],
      integer   const * LDA,
      integer           JPVT[],
      doublereal        TAU[],
      doublereal        WORK[],
      integer   const * LWORK,
      integer         * INFO
    );
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DGEQP3 computes a QR factorization with column pivoting of a
   *  matrix A:  A*P = Q*R  using Level 3 BLAS.
   *
   *  Arguments
   *  =========
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix A. M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the M-by-N matrix A.
   *          On exit, the upper triangle of the array contains the
   *          min(M,N)-by-N upper trapezoidal matrix R; the elements below
   *          the diagonal, together with the array TAU, represent the
   *          orthogonal matrix Q as a product of min(M,N) elementary
   *          reflectors.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A. LDA >= max(1,M).
   *
   *  JPVT    (input/output) INTEGER array, dimension (N)
   *          On entry, if JPVT(J).ne.0, the J-th column of A is permuted
   *          to the front of A*P (a leading column); if JPVT(J)=0,
   *          the J-th column of A is a free column.
   *          On exit, if JPVT(J)=K, then the J-th column of A*P was the
   *          the K-th column of A.
   *
   *  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
   *          The scalar factors of the elementary reflectors.
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO=0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK. LWORK >= 3*N+1.
   *          For optimal performance LWORK >= 2*N+( N+1 )*NB, where NB
   *          is the optimal blocksize.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0: successful exit.
   *          < 0: if INFO = -i, the i-th argument had an illegal value.
   *
   *  Further Details
   *  ===============
   *
   *  The matrix Q is represented as a product of elementary reflectors
   *
   *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
   *
   *  Each H(i) has the form
   *
   *     H(i) = I - tau * v * v'
   *
   *  where tau is a real/complex scalar, and v is a real/complex vector
   *  with v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in
   *  A(i+1:m,i), and tau in TAU(i).
   *
   *  Based on contributions by
   *    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
   *    X. Sun, Computer Science Dept., Duke University, USA
   *
  \*/

  inline
  integer
  geqp3(
    integer M,
    integer N,
    real    A[],
    integer LDA,
    integer JPVT[],
    real    TAU[],
    real    WORK[],
    integer LWORK
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sgeqp3)( &M, &N, A, &LDA, JPVT, TAU, WORK, &LWORK, &info );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sgeqp3( &M, &N, A, &LDA, JPVT, TAU, WORK, &LWORK, &info );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sgeqp3)( &M, &N, A, &LDA, JPVT, TAU, WORK, &LWORK, &info );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline
  integer
  geqp3(
    integer    M,
    integer    N,
    doublereal A[],
    integer    LDA,
    integer    JPVT[],
    doublereal TAU[],
    doublereal WORK[],
    integer    LWORK
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dgeqp3)( &M, &N, A, &LDA, JPVT, TAU, WORK, &LWORK, &info );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dgeqp3( &M, &N, A, &LDA, JPVT, TAU, WORK, &LWORK, &info );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dgeqp3)( &M, &N, A, &LDA, JPVT, TAU, WORK, &LWORK, &info );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  /*
  //   _                 __
  //  | |_ _____ __ ____/ _|
  //  | __|_  / '__|_  / |_
  //  | |_ / /| |   / /|  _|
  //   \__/___|_|  /___|_|
  */

  /*\
   *  Purpose
   *  =======
   *
   *  DTZRZF reduces the M-by-N ( M<=N ) real upper trapezoidal matrix A
   *  to upper triangular form by means of orthogonal transformations.
   *
   *  The upper trapezoidal matrix A is factored as
   *
   *     A = ( R  0 ) * Z,
   *
   *  where Z is an N-by-N orthogonal matrix and R is an M-by-M upper
   *  triangular matrix.
   *
   *  Arguments
   *  =========
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix A.  M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A.  N >= M.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the leading M-by-N upper trapezoidal part of the
   *          array A must contain the matrix to be factorized.
   *          On exit, the leading M-by-M upper triangular part of A
   *          contains the upper triangular matrix R, and elements M+1 to
   *          N of the first M rows of A, with the array TAU, represent the
   *          orthogonal matrix Z as a product of M elementary reflectors.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  TAU     (output) DOUBLE PRECISION array, dimension (M)
   *          The scalar factors of the elementary reflectors.
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK.  LWORK >= max(1,M).
   *          For optimum performance LWORK >= M*NB, where NB is
   *          the optimal blocksize.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *
   *  Further Details
   *  ===============
   *
   *  Based on contributions by
   *    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
   *
   *  The factorization is obtained by Householder's method.  The kth
   *  transformation matrix, Z( k ), which is used to introduce zeros into
   *  the ( m - k + 1 )th row of A, is given in the form
   *
   *     Z( k ) = ( I     0   ),
   *              ( 0  T( k ) )
   *
   *  where
   *
   *     T( k ) = I - tau*u( k )*u( k )',   u( k ) = (   1    ),
   *                                                 (   0    )
   *                                                 ( z( k ) )
   *
   *  tau is a scalar and z( k ) is an ( n - m ) element vector.
   *  tau and z( k ) are chosen to annihilate the elements of the kth row
   *  of X.
   *
   *  The scalar tau is returned in the kth element of TAU and the vector
   *  u( k ) in the kth row of A, such that the elements of z( k ) are
   *  in  a( k, m + 1 ), ..., a( k, n ). The elements of R are returned in
   *  the upper triangular part of A.
   *
   *  Z is given by
   *
   *     Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ).
   *
   *  =====================================================================
  \*/

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  // use standard Lapack routine
  extern "C" {

    void
    LAPACK_F77NAME(stzrzf)(
      integer const * M,
      integer const * N,
      real          * A,
      integer const * LDA,
      real          * TAU,
      real          * WORK,
      integer const * LWORK,
      integer       * INFO
    );

    void
    LAPACK_F77NAME(dtzrzf)(
      integer    const * M,
      integer    const * N,
      doublereal       * A,
      integer    const * LDA,
      doublereal       * TAU,
      doublereal       * WORK,
      integer    const * LWORK,
      integer          * INFO
    );

  }
  #endif

  inline
  integer
  tzrzf(
    integer M,
    integer N,
    real    A[],
    integer LDA,
    real    TAU[],
    real    WORK[],
    integer LWORK
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(stzrzf)( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    stzrzf( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(stzrzf)( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline
  integer
  tzrzf(
    integer    M,
    integer    N,
    doublereal A[],
    integer    LDA,
    doublereal TAU[],
    doublereal WORK[],
    integer    LWORK
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dtzrzf)( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dtzrzf( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dtzrzf)( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  /*
  //    ___  _ __ _ __ ___  _ __ ____
  //   / _ \| '__| '_ ` _ \| '__|_  /
  //  | (_) | |  | | | | | | |   / /
  //   \___/|_|  |_| |_| |_|_|  /___|
  */

  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS)
  // use standard Lapack routine
  extern "C" {

    void
    LAPACK_F77NAME(sormrz)(
      character const * SIDE,
      character const * TRANS,
      integer   const * M,
      integer   const * N,
      integer   const * K,
      integer   const * L,
      real      const   A[],
      integer   const * LDA,
      real      const   TAU[],
      real              C[],
      integer   const * LDC,
      real              WORK[],
      integer   const * LWORK,
      integer         * INFO
    );

    void
    LAPACK_F77NAME(dormrz)(
      character  const * SIDE,
      character  const * TRANS,
      integer    const * M,
      integer    const * N,
      integer    const * K,
      integer    const * L,
      doublereal const   A[],
      integer    const * LDA,
      doublereal const   TAU[],
      doublereal         C[],
      integer    const * LDC,
      doublereal         WORK[],
      integer    const * LWORK,
      integer          * INFO
    );

  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DORMRZ overwrites the general real M-by-N matrix C with
   *
   *                  SIDE = 'L'     SIDE = 'R'
   *  TRANS = 'N':      Q * C          C * Q
   *  TRANS = 'T':      Q**T * C       C * Q**T
   *
   *  where Q is a real orthogonal matrix defined as the product of k
   *  elementary reflectors
   *
   *        Q = H(1) H(2) . . . H(k)
   *
   *  as returned by DTZRZF. Q is of order M if SIDE = 'L' and of order N
   *  if SIDE = 'R'.
   *
   *  Arguments
   *  =========
   *
   *  SIDE    (input) CHARACTER*1
   *          = 'L': apply Q or Q**T from the Left;
   *          = 'R': apply Q or Q**T from the Right.
   *
   *  TRANS   (input) CHARACTER*1
   *          = 'N':  No transpose, apply Q;
   *          = 'T':  Transpose, apply Q**T.
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix C. M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix C. N >= 0.
   *
   *  K       (input) INTEGER
   *          The number of elementary reflectors whose product defines
   *          the matrix Q.
   *          If SIDE = 'L', M >= K >= 0;
   *          if SIDE = 'R', N >= K >= 0.
   *
   *  L       (input) INTEGER
   *          The number of columns of the matrix A containing
   *          the meaningful part of the Householder reflectors.
   *          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.
   *
   *  A       (input) DOUBLE PRECISION array, dimension
   *                               (LDA,M) if SIDE = 'L',
   *                               (LDA,N) if SIDE = 'R'
   *          The i-th row must contain the vector which defines the
   *          elementary reflector H(i), for i = 1,2,...,k, as returned by
   *          DTZRZF in the last k rows of its array argument A.
   *          A is modified by the routine but restored on exit.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A. LDA >= max(1,K).
   *
   *  TAU     (input) DOUBLE PRECISION array, dimension (K)
   *          TAU(i) must contain the scalar factor of the elementary
   *          reflector H(i), as returned by DTZRZF.
   *
   *  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
   *          On entry, the M-by-N matrix C.
   *          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
   *
   *  LDC     (input) INTEGER
   *          The leading dimension of the array C. LDC >= max(1,M).
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK.
   *          If SIDE = 'L', LWORK >= max(1,N);
   *          if SIDE = 'R', LWORK >= max(1,M).
   *          For optimum performance LWORK >= N*NB if SIDE = 'L', and
   *          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
   *          blocksize.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *
   *  Further Details
   *  ===============
   *
   *  Based on contributions by
   *    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
   *
   *  =====================================================================
  \*/

  inline
  integer
  ormrz(
    SideMultiply  const & SIDE,
    Transposition const & TRANS,
    integer               M,
    integer               N,
    integer               K,
    integer               L,
    real const            A[],
    integer               LDA,
    real const            TAU[],
    real                  C[],
    integer               LDC,
    real                  WORK[],
    integer               LWORK
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(sormrz)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(trans_blas[TRANS]),
      &M, &N, &K, &L, A, &LDA, TAU,
      C, &LDC, WORK, &LWORK, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    sormrz(
      side_blas[SIDE], trans_blas[TRANS],
      &M, &N, &K, &L, A, &LDA, TAU,
      C, &LDC, WORK, &LWORK, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(sormrz)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(trans_blas[TRANS]),
      &M, &N, &K, &L,
      const_cast<real*>(A), &LDA,
      const_cast<real*>(TAU),
      C, &LDC, WORK, &LWORK, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline
  integer
  ormrz(
    SideMultiply  const & SIDE,
    Transposition const & TRANS,
    integer               M,
    integer               N,
    integer               K,
    integer               L,
    doublereal const      A[],
    integer               LDA,
    doublereal const      TAU[],
    doublereal            C[],
    integer               LDC,
    doublereal            WORK[],
    integer               LWORK
  ) {
    integer info = 0;
    #if defined(LAPACK_WRAPPER_USE_LAPACK)   || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)
    LAPACK_F77NAME(dormrz)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(trans_blas[TRANS]),
      &M, &N, &K, &L, A, &LDA, TAU,
      C, &LDC, WORK, &LWORK, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_MKL)
    dormrz(
      side_blas[SIDE], trans_blas[TRANS],
      &M, &N, &K, &L, A, &LDA, TAU,
      C, &LDC, WORK, &LWORK, &info
    );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
    CLAPACKNAME(dormrz)(
      const_cast<character*>(side_blas[SIDE]),
      const_cast<character*>(trans_blas[TRANS]),
      &M, &N, &K, &L,
      const_cast<doublereal*>(A), &LDA,
      const_cast<doublereal*>(TAU),
      C, &LDC, WORK, &LWORK, &info
    );
    #else
    #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

}
