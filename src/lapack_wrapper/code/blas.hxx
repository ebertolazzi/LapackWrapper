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
/// file: blas.hxx
///

namespace lapack_wrapper {

  /*
  //                              __   __ _ _ _
  //    ___ ___  _ __  _   _     / /  / _(_) | |
  //   / __/ _ \| '_ \| | | |   / /  | |_| | | |
  //  | (_| (_) | |_) | |_| |  / /   |  _| | | |
  //   \___\___/| .__/ \__, | /_/    |_| |_|_|_|
  //            |_|    |___/
  */
  #ifdef LAPACK_WRAPPER_USE_LAPACK
  extern "C" {
    void
    BLASFUNC(scopy)(
      integer const * N,
      real    const   X[],
      integer const * INCX,
      real            Y[],
      integer const * INCY
    );

    void
    BLASFUNC(dcopy)(
      integer    const * N,
      doublereal const   X[],
      integer    const * INCX,
      doublereal         Y[],
      integer    const * INCY
    );
  }
  #endif

  inline
  void
  copy(
    integer    N,
    real const X[],
    integer    INCX,
    real       Y[],
    integer    INCY
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { scopy( &N, X, &INCX, Y, &INCY ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { CBLASNAME(scopy)( N, X, INCX, Y, INCY ); }
  #else
  { BLASFUNC(scopy)( &N, const_cast<real*>(X), &INCX, Y, &INCY ); }
  #endif

  inline
  void
  copy(
    integer          N,
    doublereal const X[],
    integer          INCX,
    doublereal       Y[],
    integer          INCY
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { dcopy( &N, X, &INCX, Y, &INCY ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { CBLASNAME(dcopy)( N, X, INCX, Y, INCY ); }
  #else
  { BLASFUNC(dcopy)( &N, const_cast<doublereal*>(X), &INCX, Y, &INCY ); }
  #endif

  inline
  void
  fill(
    integer N,
    real    Y[],
    integer INCY,
    real    V
  )
  { copy( N, &V, 0, Y, INCY ); }

  inline
  void
  fill(
    integer    N,
    doublereal Y[],
    integer    INCY,
    doublereal V
  )
  { copy( N, &V, 0, Y, INCY ); }

  /*
  //   _____      ____ _ _ __
  //  / __\ \ /\ / / _` | '_ \
  //  \__ \\ V  V / (_| | |_) |
  //  |___/ \_/\_/ \__,_| .__/
  //                    |_|
  */
  #if defined(LAPACK_WRAPPER_USE_LAPACK) || \
      defined(LAPACK_WRAPPER_USE_ATLAS) 
  extern "C" {
    void
    BLASFUNC(sswap)(
      integer const * N,
      real            X[],
      integer const * INCX,
      real            Y[],
      integer const * INCY
    );

    void
    BLASFUNC(dswap)(
      integer const * N,
      doublereal      X[],
      integer const * INCX,
      doublereal      Y[],
      integer const * INCY
    );

    void
    BLASFUNC(slaswp)(
      integer const * NCOL,
      real            A[],
      integer const * LDA,
      integer const * K1,
      integer const * K2,
      integer const   IPIV[],
      integer const * INC
    );

    void
    BLASFUNC(dlaswp)(
      integer const * NCOL,
      doublereal      A[],
      integer const * LDA,
      integer const * K1,
      integer const * K2,
      integer const   IPIV[],
      integer const * INC
    );
  }
  #endif

  inline
  void
  swap(
    integer N,
    real    X[],
    integer INCX,
    real    Y[],
    integer INCY
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { sswap( &N, X, &INCX, Y, &INCY ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { CBLASNAME(sswap)( N, X, INCX, Y, INCY ); }
  #else
  { BLASFUNC(sswap)( &N, X, &INCX, Y, &INCY ); }
  #endif

  inline
  void
  swap(
    integer    N,
    doublereal X[],
    integer    INCX,
    doublereal Y[],
    integer    INCY
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { dswap( &N, X, &INCX, Y, &INCY ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { CBLASNAME(dswap)( N, X, INCX, Y, INCY ); }
  #else
  { BLASFUNC(dswap)( &N, X, &INCX, Y, &INCY ); }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DLASWP performs a series of row interchanges on the matrix A.
   *  One row interchange is initiated for each of rows K1 through K2 of A.
   *
   *  Arguments
   *  =========
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the matrix of column dimension N to which the row
   *          interchanges will be applied.
   *          On exit, the permuted matrix.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.
   *
   *  K1      (input) INTEGER
   *          The first element of IPIV for which a row interchange will
   *          be done.
   *
   *  K2      (input) INTEGER
   *          The last element of IPIV for which a row interchange will
   *          be done.
   *
   *  IPIV    (input) INTEGER array, dimension (K2*abs(INCX))
   *          The vector of pivot indices.  Only the elements in positions
   *          K1 through K2 of IPIV are accessed.
   *          IPIV(K) = L implies rows K and L are to be interchanged.
   *
   *  INCX    (input) INTEGER
   *          The increment between successive values of IPIV.  If IPIV
   *          is negative, the pivots are applied in reverse order.
  \*/

  inline
  integer
  swaps(
    integer       NCOL,
    real          A[],
    integer       LDA,
    integer       I1,
    integer       I2,
    integer const IPIV[],
    integer       INC
  ) {
    integer info(0);
    integer K1 = I1+1;
    integer K2 = I2+1;
    #if defined(LAPACK_WRAPPER_USE_MKL)
      slaswp( &NCOL, A, &LDA, &K1, &K2, IPIV, &INC );
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
      LAPACK_F77NAME(slaswp)(
        &NCOL, A, &LDA, &K1, &K2, const_cast<integer*>(IPIV), &INC
      );
    #elif defined(LAPACK_WRAPPER_USE_LAPACK)   || \
          defined(LAPACK_WRAPPER_USE_ATLAS)
      BLASFUNC(slaswp)(
        &NCOL, A, &LDA, &K1, &K2, const_cast<integer*>(IPIV), &INC
      );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
      info = CLAPACKNAME(slaswp)(
        &NCOL, A, &LDA, &K1, &K2, const_cast<integer*>(IPIV), &INC
      );
    #else
      #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  inline
  integer
  swaps(
    integer       NCOL,
    doublereal    A[],
    integer       LDA,
    integer       I1,
    integer       I2,
    integer const IPIV[],
    integer       INC
  ) {
    integer info(0);
    integer K1 = I1+1;
    integer K2 = I2+1;
    #if defined(LAPACK_WRAPPER_USE_MKL)
      dlaswp( &NCOL, A, &LDA, &K1, &K2, IPIV, &INC );
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS)
      LAPACK_F77NAME(dlaswp)(
        &NCOL, A, &LDA, &K1, &K2, const_cast<integer*>(IPIV), &INC
      );
    #elif defined(LAPACK_WRAPPER_USE_LAPACK) || \
          defined(LAPACK_WRAPPER_USE_ATLAS)
      BLASFUNC(dlaswp)(
        &NCOL, A, &LDA, &K1, &K2, const_cast<integer*>(IPIV), &INC
      );
    #elif defined(LAPACK_WRAPPER_USE_ACCELERATE)
      info = CLAPACKNAME(dlaswp)(
        &NCOL, A, &LDA, &K1, &K2, const_cast<integer*>(IPIV), &INC
      );
    #else
      #error "LapackWrapper undefined mapping!"
    #endif
    return info;
  }

  /*
  //                 _
  //   ___  ___ __ _| |
  //  / __|/ __/ _` | |
  //  \__ \ (_| (_| | |
  //  |___/\___\__,_|_|
  */
  #ifdef LAPACK_WRAPPER_USE_LAPACK
  extern "C" {
    void
    BLASFUNC(sscal)(
      integer const * N,
      real    const * S,
      real            X[],
      integer const * INCX
    );

    void
    BLASFUNC(dscal)(
      integer    const * N,
      doublereal const * S,
      doublereal         X[],
      integer    const * INCX
    );

    void
    BLASFUNC(srscl)(
      integer const * N,
      real    const * SA,
      real          * SX,
      integer const * INCX
    );

    void
    BLASFUNC(drscl)(
      integer    const * N,
      doublereal const * SA,
      doublereal       * SX,
      integer    const * INCX
    );
  }
  #endif

  inline
  void
  scal(
    integer N,
    real    S,
    real    X[],
    integer INCX
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { sscal( &N, &S, X, &INCX ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS) || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { CBLASNAME(sscal)( N, S, X, INCX ); }
  #else
  { BLASFUNC(sscal)( &N, &S, X, &INCX ); }
  #endif

  inline
  void
  scal(
    integer    N,
    doublereal S,
    doublereal X[],
    integer    INCX
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { dscal( &N, &S, X, &INCX ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { CBLASNAME(dscal)( N, S, X, INCX ); }
  #else
  { BLASFUNC(dscal)( &N, &S, X, &INCX ); }
  #endif

  inline
  void
  rscal(
    integer N,
    real    S,
    real    X[],
    integer INCX
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { srscl( &N, &S, X, &INCX ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { real rS = 1/S; CBLASNAME(sscal)( N, rS, X, INCX ); }
  #else
  { BLASFUNC(srscl)( &N, &S, X, &INCX ); }
  #endif

  inline
  void
  rscal(
    integer    N,
    doublereal S,
    doublereal X[],
    integer    INCX
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { drscl( &N, &S, X, &INCX ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { doublereal rS = 1/S; CBLASNAME(dscal)( N, rS, X, INCX ); }
  #else
  { BLASFUNC(drscl)( &N, &S, X, &INCX ); }
  #endif

  /*
  //    __ ___  ___ __  _   _
  //   / _` \ \/ / '_ \| | | |
  //  | (_| |>  <| |_) | |_| |
  //   \__,_/_/\_\ .__/ \__, |
  //             |_|    |___/
  */
  #ifdef LAPACK_WRAPPER_USE_LAPACK
  extern "C" {
    void
    BLASFUNC(saxpy)(
      integer const * N,
      real    const * A,
      real    const   X[],
      integer const * INCX,
      real            Y[],
      integer const * INCY
    );

    void
    BLASFUNC(daxpy)(
      integer    const * N,
      doublereal const * A,
      doublereal const   X[],
      integer    const * INCX,
      doublereal         Y[],
      integer    const * INCY
    );
  }
  #endif

  inline
  void
  axpy(
    integer    N,
    real       A,
    real const X[],
    integer    INCX,
    real       Y[],
    integer    INCY
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { saxpy( &N, &A, X, &INCX, Y, &INCY ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { CBLASNAME(saxpy)( N, A, X, INCX, Y, INCY ); }
  #else
  { BLASFUNC(saxpy)( &N, &A, const_cast<real*>(X), &INCX, Y, &INCY ); }
  #endif

  inline
  void
  axpy(
    integer          N,
    doublereal       A,
    doublereal const X[],
    integer          INCX,
    doublereal       Y[],
    integer          INCY
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { daxpy( &N, &A, X, &INCX, Y, &INCY ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { CBLASNAME(daxpy)( N, A, X, INCX, Y, INCY ); }
  #else
  { BLASFUNC(daxpy)( &N, &A, const_cast<doublereal*>(X), &INCX, Y, &INCY ); }
  #endif

  /*
  //   _______ _ __ ___
  //  |_  / _ \ '__/ _ \
  //   / /  __/ | | (_) |
  //  /___\___|_|  \___/
  */
  inline
  void
  zero(
    integer N,
    real    X[],
    integer INCX
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { real z = 0; integer iz = 0; scopy( &N, &z, &iz, X, &INCX ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { real z = 0; CBLASNAME(scopy)( N, &z, 0, X, INCX ); }
  #else
  { real z = 0; integer iz = 0; BLASFUNC(scopy)( &N, &z, &iz, X, &INCX ); }
  #endif

  inline
  void
  zero(
    integer    N,
    doublereal X[],
    integer    INCX
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { doublereal z = 0; integer iz = 0; dcopy( &N, &z, &iz, X, &INCX ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { doublereal z = 0; CBLASNAME(dcopy)( N, &z, 0, X, INCX ); }
  #else
  { doublereal z = 0; integer iz = 0; BLASFUNC(dcopy)( &N, &z, &iz, X, &INCX ); }
  #endif

  /*
  //             _
  //   _ __ ___ | |_
  //  | '__/ _ \| __|
  //  | | | (_) | |_
  //  |_|  \___/ \__|
  */
  #ifdef LAPACK_WRAPPER_USE_LAPACK
  extern "C" {
    void
    BLASFUNC(srot)(
      integer const * N,
      real            DX[],
      integer const * INCX,
      real            DY[],
      integer const * INCY,
      real    const * C,
      real    const * S
    );

    void
    BLASFUNC(drot)(
      integer    const * N,
      doublereal         DX[],
      integer    const * INCX,
      doublereal         DY[],
      integer    const * INCY,
      doublereal const * C,
      doublereal const * S
    );
  }
  #endif

  inline
  void
  rot(
    integer N,
    real    DX[],
    integer INCX,
    real    DY[],
    integer INCY,
    real    C,
    real    S
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { srot( &N, DX, &INCX, DY, &INCY, &C, &S ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { CBLASNAME(srot)( N, DX, INCX, DY, INCY, C, S ); }
  #else
  { BLASFUNC(srot)( &N, DX, &INCX, DY, &INCY, &C, &S ); }
  #endif

  inline
  void
  rot(
    integer    N,
    doublereal DX[],
    integer    INCX,
    doublereal DY[],
    integer    INCY,
    doublereal C,
    doublereal S
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { drot( &N, DX, &INCX, DY, &INCY, &C, &S ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { CBLASNAME(drot)( N, DX, INCX, DY, INCY, C, S ); }
  #else
  { BLASFUNC(drot)( &N, DX, &INCX, DY, &INCY, &C, &S ); }
  #endif

  /*
  //             _
  //   _ __ ___ | |_ __ _
  //  | '__/ _ \| __/ _` |
  //  | | | (_) | || (_| |
  //  |_|  \___/ \__\__, |
  //                |___/
  */
  /*\
   *   Construct the Givens transformation
   *
   *         ( DC  DS )
   *     G = (        ) ,    DC**2 + DS**2 = 1 ,
   *         (-DS  DC )
   *
   *     which zeros the second entry of the 2-vector  (DA,DB)**T .
   *
   *     The quantity R = (+/-)SQRT(DA**2 + DB**2) overwrites DA in
   *     storage.  The value of DB is overwritten by a value Z which
   *     allows DC and DS to be recovered by the following algorithm.
  \*/
  #ifdef LAPACK_WRAPPER_USE_LAPACK
  extern "C" {
    void
    BLASFUNC(srotg)(
      real * DX,
      real * DY,
      real * C,
      real * S
    );

    void
    BLASFUNC(drotg)(
      doublereal * DX,
      doublereal * DY,
      doublereal * C,
      doublereal * S
    );
  }
  #endif

  inline
  void
  rotg(
    real & DX,
    real & DY,
    real & C,
    real & S
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { srotg( &DX, &DY, &C, &S ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { CBLASNAME(srotg)( &DX, &DY, &C, &S ); }
  #else
  { BLASFUNC(srotg)( &DX, &DY, &C, &S ); }
  #endif

  inline
  void
  rotg(
    doublereal & DX,
    doublereal & DY,
    doublereal & C,
    doublereal & S
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { drotg( &DX, &DY, &C, &S ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { CBLASNAME(drotg)( &DX, &DY, &C, &S ); }
  #else
  { BLASFUNC(drotg)( &DX, &DY, &C, &S ); }
  #endif

  /*
  //                        ____
  //   _ __  _ __ _ __ ___ |___ \
  //  | '_ \| '__| '_ ` _ \  __) |
  //  | | | | |  | | | | | |/ __/
  //  |_| |_|_|  |_| |_| |_|_____|
  */
  /*\
   *  DNRM2 returns the euclidean norm of a vector via the function
   *  name, so that
   *
   *     DNRM2 := sqrt( x'*x )
   *
  \*/
  #ifdef LAPACK_WRAPPER_USE_LAPACK
  extern "C" {
    real
    BLASFUNC(snrm2)(
      integer const * N,
      real    const   X[],
      integer const * INCX
    );

    doublereal
    BLASFUNC(dnrm2)(
      integer    const * N,
      doublereal const   X[],
      integer    const * INCX
    );
  }
  #endif

  inline
  real
  nrm2(
    integer    N,
    real const X[],
    integer    INCX
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { return snrm2( &N, X, &INCX ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { return CBLASNAME(snrm2)( N, X, INCX ); }
  #else
  { return BLASFUNC(snrm2)( &N, const_cast<real*>(X), &INCX ); }
  #endif

  inline
  doublereal
  nrm2(
    integer          N,
    doublereal const X[],
    integer          INCX
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { return dnrm2( &N, X, &INCX ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { return CBLASNAME(dnrm2)( N, X, INCX ); }
  #else
  { return BLASFUNC(dnrm2)( &N, const_cast<doublereal*>(X), &INCX ); }
  #endif

  /*
  //    __ _ ___ _   _ _ __ ___
  //   / _` / __| | | | '_ ` _ \
  //  | (_| \__ \ |_| | | | | | |
  //   \__,_|___/\__,_|_| |_| |_|
  */
  /*\
   * Purpose
   * =======
   *
   *  DASUM takes the sum of the absolute values.
  \*/
  #ifdef LAPACK_WRAPPER_USE_LAPACK
  extern "C" {
    real
    BLASFUNC(sasum)(
      integer const * N,
      real    const   X[],
      integer const * INCX
    );

    doublereal
    BLASFUNC(dasum)(
      integer    const * N,
      doublereal const   X[],
      integer    const * INCX
    );
  }
  #endif

  inline
  real
  asum(
    integer    N,
    real const X[],
    integer    INCX
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { return sasum( &N, X, &INCX ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { return CBLASNAME(sasum)( N, X, INCX ); }
  #else
  { return BLASFUNC(sasum)( &N, const_cast<real*>(X), &INCX ); }
  #endif

  inline
  doublereal
  asum(
    integer          N,
    doublereal const X[],
    integer          INCX
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { return dasum( &N, X, &INCX ); }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { return CBLASNAME(dasum)( N, X, INCX ); }
  #else
  { return BLASFUNC(dasum)( &N, const_cast<doublereal*>(X), &INCX ); }
  #endif

  /*
  //    __ _ _ __ ___   __ ___  __
  //   / _` | '_ ` _ \ / _` \ \/ /
  //  | (_| | | | | | | (_| |>  <
  //   \__,_|_| |_| |_|\__,_/_/\_\
  */
  /*\
   *  Purpose
   *  =======
   *
   *     IDAMAX finds the index of element having max. absolute value.
  \*/
  #ifdef LAPACK_WRAPPER_USE_LAPACK
  extern "C" {
    integer
    BLASFUNC(isamax)(
      integer const * N,
      real    const   X[],
      integer const * INCX
    );

    integer
    BLASFUNC(idamax)(
      integer    const * N,
      doublereal const   X[],
      integer    const * INCX
    );
  }
  #endif

  inline
  integer
  iamax(
    integer    N,
    real const X[],
    integer    INCX
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { return isamax( &N, X, &INCX )-1; }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { return integer(CBLASNAME(isamax)( N, X, INCX )); }
  #else
  { return BLASFUNC(isamax)( &N, const_cast<real*>(X), &INCX )-1; }
  #endif

  inline
  integer
  iamax(
    integer          N,
    doublereal const X[],
    integer          INCX
  )
  #if defined(LAPACK_WRAPPER_USE_MKL)
  { return idamax( &N, X, &INCX )-1; }
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
  { return integer(CBLASNAME(idamax)( N, X, INCX )); }
  #else
  { return BLASFUNC(idamax)( &N, const_cast<doublereal*>(X), &INCX )-1; }
  #endif

  inline
  real
  absmax(
    integer    N,
    real const X[],
    integer    INCX
  )
  { real tmp = X[iamax(N,X,INCX)]; return tmp > 0 ? tmp : -tmp; }

  inline
  doublereal
  absmax(
    integer          N,
    doublereal const X[],
    integer          INCX
  )
  { doublereal tmp = X[iamax(N,X,INCX)]; return tmp > 0 ? tmp : -tmp; }

  /*
  //       _       _
  //    __| | ___ | |_
  //   / _` |/ _ \| __|
  //  | (_| | (_) | |_
  //   \__,_|\___/ \__|
  */
  /*\
   *  Purpose
   *  =======
   *
   *     DDOT forms the dot product of two vectors.
   *     uses unrolled loops for increments equal to one.
  \*/
  #ifdef LAPACK_WRAPPER_USE_LAPACK
  extern "C" {
    real
    BLASFUNC(sdot)(
      integer const * N,
      real    const   SX[],
      integer const * INCX,
      real    const   SY[],
      integer const * INCY
    );

    doublereal
    BLASFUNC(ddot)(
      integer    const * N,
      doublereal const   SX[],
      integer    const * INCX,
      doublereal const   SY[],
      integer    const * INCY
    );
  }
  #endif

  inline
  real
  dot(
    integer    N,
    real const SX[],
    integer    INCX,
    real const SY[],
    integer    INCY
  ) {
  #if defined(LAPACK_WRAPPER_USE_MKL)
    return sdot( &N, SX, &INCX, SY, &INCY );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    return CBLASNAME(sdot)( N, SX, INCX, SY, INCY );
  #else
    return BLASFUNC(sdot)(
      &N,
      const_cast<real*>(SX), &INCX,
      const_cast<real*>(SY), &INCY
    );
  #endif
  }

  inline
  doublereal
  dot(
    integer          N,
    doublereal const SX[],
    integer          INCX,
    doublereal const SY[],
    integer          INCY
  ) {
  #if defined(LAPACK_WRAPPER_USE_MKL)
    return ddot( &N, SX, &INCX, SY, &INCY );
  #elif defined(LAPACK_WRAPPER_USE_ACCELERATE) || \
        defined(LAPACK_WRAPPER_USE_ATLAS)      || \
        defined(LAPACK_WRAPPER_USE_OPENBLAS)
    return CBLASNAME(ddot)( N, SX, INCX, SY, INCY );
  #else
    return BLASFUNC(ddot)(
      &N,
      const_cast<doublereal*>(SX), &INCX,
      const_cast<doublereal*>(SY), &INCY
    );
  #endif
  }

}
