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
/// file: trid.cxx
///

namespace lapack_wrapper {

  //============================================================================

  template <typename valueType>
  inline
  void
  tridiag_axpy(
    integer         N,
    valueType       alpha,
    valueType const L[],
    valueType const D[],
    valueType const U[],
    valueType const x[],
    valueType       beta,
    valueType       y[]
  ) {
    if ( isZero(beta) ) {
      y[0] = alpha*(D[0]*x[0] + U[0] * x[1]);
      for ( integer i = 1; i < N-1; ++i )
        y[i] = alpha*(D[i]*x[i] + U[i] * x[i+1] + L[i-1] * x[i-1]);
      y[N-1] = alpha*(D[N-1]*x[N-1] + L[N-2] * x[N-2]);
    } else {
      y[0] = beta*y[0] + alpha*(D[0]*x[0] + U[0] * x[1]);
      for ( integer i = 1; i < N-1; ++i )
        y[i] = beta*y[i] + alpha*(D[i]*x[i] + U[i] * x[i+1] + L[i-1] * x[i-1]);
      y[N-1] = beta*y[N-1] + alpha*(D[N-1]*x[N-1] + L[N-2] * x[N-2]);
    }
  }

  //============================================================================
  /*\
   |   _____     _     _ _                               _ ____  ____  ____
   |  |_   _| __(_) __| (_) __ _  __ _  ___  _ __   __ _| / ___||  _ \|  _ \
   |    | || '__| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | \___ \| |_) | | | |
   |    | || |  | | (_| | | (_| | (_| | (_) | | | | (_| | |___) |  __/| |_| |
   |    |_||_|  |_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|____/|_|   |____/
   |                             |___/
  \*/
  template <typename T>
  void
  TridiagonalSPD<T>::factorize(
    char const      who[],
    integer         N,
    valueType const _L[],
    valueType const _D[]
  ) {
    if ( nRC != N ) {
      nRC = N;
      allocReals.allocate(3*N);
      L    = allocReals(N);
      D    = allocReals(N);
      WORK = allocReals(N);
    }
    copy( N, _L, 1, L, 1 );
    copy( N, _D, 1, D, 1 );
    integer info = pttrf( N, L, D );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalSPD::factorize[" << who <<
      "], return info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  T
  TridiagonalSPD<T>::cond1( valueType norm1 ) const {
    valueType rcond;
    integer info = ptcon1( nRC, D, L, norm1, rcond, WORK );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalSPD::cond1, return info = " << info
    );
    return rcond;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalSPD<T>::solve( valueType xb[] ) const {
    integer info = pttrs( nRC, 1, D, L, xb, nRC );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalSPD::solve, return info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalSPD<T>::t_solve( valueType xb[] ) const {
    integer info = pttrs( nRC, 1, D, L, xb, nRC );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalSPD::solve, return info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalSPD<T>::solve( integer nrhs, valueType xb[], integer ldXB ) const {
    integer info = pttrs( nRC, nrhs, D, L, xb, ldXB );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalSPD::solve, return info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalSPD<T>::t_solve( integer nrhs, valueType xb[], integer ldXB ) const {
    integer info = pttrs( nRC, nrhs, D, L, xb, ldXB );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalSPD::solve, return info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalSPD<T>::axpy(
    integer         N,
    valueType       alpha,
    valueType const _L[],
    valueType const _D[],
    valueType const x[],
    valueType       beta,
    valueType       y[]
  ) const {
    tridiag_axpy( N, alpha, _L, _D, _L, x, beta, y );
  }

  //============================================================================
  /*\
   |   _____     _     _ _                               _ _    _   _
   |  |_   _| __(_) __| (_) __ _  __ _  ___  _ __   __ _| | |  | | | |
   |    | || '__| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | | |  | | | |
   |    | || |  | | (_| | | (_| | (_| | (_) | | | | (_| | | |__| |_| |
   |    |_||_|  |_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|_____\___/
   |                             |___/
  \*/
  template <typename T>
  void
  TridiagonalLU<T>::factorize(
    char const      who[],
    integer         N,
    valueType const _L[],
    valueType const _D[],
    valueType const _U[]
  ) {
    if ( nRC != N ) {
      nRC = N;
      allocReals.allocate(6*N);
      allocIntegers.allocate(2*N);
      L     = allocReals(N);
      D     = allocReals(N);
      U     = allocReals(N);
      U2    = allocReals(N);
      WORK  = allocReals(2*N);
      IPIV  = allocIntegers(N);
      IWORK = allocIntegers(N);
    }
    copy( N, _L, 1, L, 1 );
    copy( N, _D, 1, D, 1 );
    copy( N, _U, 1, U, 1 );
    integer info = gttrf( N, L, D, U, U2, IPIV );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalLU::factorize[" << who <<
      "], return info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  T
  TridiagonalLU<T>::cond1( valueType norm1 ) const {
    valueType rcond;
    integer info = gtcon1( nRC, L, D, U, U2, IPIV, norm1, rcond, WORK, IWORK );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalLU::cond1, return info = " << info
    );
    return rcond;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  T
  TridiagonalLU<T>::condInf( valueType normInf ) const {
    valueType rcond;
    integer info = gtconInf( nRC, L, D, U, U2, IPIV, normInf, rcond, WORK, IWORK );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalLU::cond1, return info = " << info
    );
    return rcond;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalLU<T>::solve( valueType xb[] ) const {
    integer info = gttrs( NO_TRANSPOSE, nRC, 1, L, D, U, U2, IPIV, xb, nRC );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalLU::solve, return info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalLU<T>::t_solve( valueType xb[] ) const {
    integer info = gttrs( TRANSPOSE, nRC, 1, L, D, U, U2, IPIV, xb, nRC );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalLU::solve, return info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalLU<T>::solve( integer nrhs, valueType xb[], integer ldXB ) const {
    integer info = gttrs( NO_TRANSPOSE, nRC, nrhs, L, D, U, U2, IPIV, xb, ldXB );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalLU::solve, return info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalLU<T>::t_solve( integer nrhs, valueType xb[], integer ldXB ) const {
    integer info = gttrs( TRANSPOSE, nRC, nrhs, L, D, U, U2, IPIV, xb, ldXB );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalLU::solve, return info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalLU<T>::axpy(
    integer         N,
    valueType       alpha,
    valueType const _L[],
    valueType const _D[],
    valueType const _U[],
    valueType const x[],
    valueType       beta,
    valueType       y[]
  ) const {
    tridiag_axpy( N, alpha, _L, _D, _U, x, beta, y );
  }

  //============================================================================
  /*\
   |   _____     _     _ _                               _  ___  ____
   |  |_   _| __(_) __| (_) __ _  __ _  ___  _ __   __ _| |/ _ \|  _ \
   |    | || '__| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | | | | | |_) |
   |    | || |  | | (_| | | (_| | (_| | (_) | | | | (_| | | |_| |  _ <
   |    |_||_|  |_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|\__\_\_| \_\
   |                             |___/
  \*/

  template <typename T>
  void
  TridiagonalQR<T>::factorize(
    char const      /* who */[],
    integer         N,
    valueType const L[],
    valueType const D[],
    valueType const U[]
  ) {
    allocReals.allocate(size_t(5*(N-1)));
    nRC = N;
    C   = allocReals(size_t(N-1));
    S   = allocReals(size_t(N-1));
    BD  = allocReals(size_t(N));
    BU  = allocReals(size_t(N-1));
    BU2 = allocReals(size_t(N-2));

    /*\
      | d u       | d u @     | d u @     | d u @     | d u @     |
      | l d u     | 0 d u     | 0 d u @   | 0 d u @   | 0 d u @   |
      |   l d u   |   l d u   |   0 d u   |   0 d u @ |   0 d u @ |
      |     l d u |     l d u |     l d u |     0 d u |     0 d u |
      |       l d |       l d |       l d |       l d |       0 d |
    \*/

    lapack_wrapper::copy( N,   D, 1, BD, 1 );
    lapack_wrapper::copy( N-1, U, 1, BU, 1 );
    lapack_wrapper::zero( N-2, BU2, 1 );

    normInfA = 0;
    integer i = 0;
    for (; i < N-2; ++i ) {
      valueType Li = L[i];
      rotg( BD[i], Li, C[i], S[i] );
      rot( 1, &BU[i],  1, &BD[i+1], 1, C[i], S[i] );
      rot( 1, &BU2[i], 1, &BU[i+1], 1, C[i], S[i] );
      valueType sum = std::abs(BD[i]) + std::abs(BU[i]) + std::abs(BU2[i]);
      if ( sum > normInfA ) normInfA = sum;
    }
    valueType Li = L[i];
    rotg( BD[i], Li, C[i], S[i] );
    rot( 1, &BU[i], 1, &BD[i+1], 1, C[i], S[i] );

    valueType sum = std::abs(BD[i]) + std::abs(BU[i]);
    if ( sum > normInfA ) normInfA = sum;
    sum = std::abs(BD[i+1]);
    if ( sum > normInfA ) normInfA = sum;

    // Q A = R
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalQR<T>::Rsolve( valueType xb[] ) const {
    xb[nRC-1] /= BD[nRC-1];
    xb[nRC-2] = (xb[nRC-2]-BU[nRC-2]*xb[nRC-1])/BD[nRC-2];
    for ( integer i = nRC-3; i >= 0; --i )
      xb[i] = (xb[i]-BU[i]*xb[i+1]-BU2[i]*xb[i+2])/BD[i];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalQR<T>::RsolveTransposed( valueType xb[] ) const {
    xb[0] /= BD[0];
    xb[1] = (xb[1]-BU[0]*xb[0])/BD[1];
    for ( integer i = 2; i < nRC; ++i )
      xb[i] = (xb[i]-BU[i]*xb[i-1]-BU2[i]*xb[i-2])/BD[i];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalQR<T>::solve( valueType xb[] ) const {
    // A x = b --> Q A x = Q b --> R x = Q b
    // applico Q b
    for ( integer i = 0; i < nRC-1; ++i )
      rot( 1, &xb[i], 1, &xb[i+1], 1, C[i], S[i] );
    Rsolve( xb );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalQR<T>::t_solve( valueType xb[] ) const {
    // A^T x = b --> A^T Q^T Q x = b --> R^T Q x = b --> R^T y = b  x = Q^T y
    RsolveTransposed( xb );
    // applico Q^T b
    for ( integer i = nRC-2; i >= 0; --i )
      rot( 1, &xb[i], 1, &xb[i+1], 1, C[i], -S[i] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalQR<T>::solve( integer nrhs, valueType xb[], integer ldXB ) const {
    // A x = b --> Q A x = Q b --> R x = Q b
    // applico Q b
    for ( integer i = 0; i < nRC-1; ++i )
      rot( nrhs, &xb[i], ldXB, &xb[i+1], ldXB, C[i], S[i] );
    for ( integer i = 0; i < nrhs; ++i )
      Rsolve( xb+i*ldXB );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalQR<T>::t_solve( integer nrhs, valueType xb[], integer ldXB ) const {
    // A^T x = b --> A^T Q^T Q x = b --> R^T Q x = b --> R^T y = b  x = Q^T y
    for ( integer i = 0; i < nrhs; ++i )
      RsolveTransposed(xb+i*ldXB);
    for ( integer i = nRC-2; i >= 0; --i )
      rot( nrhs, &xb[i], ldXB, &xb[i+1], ldXB, C[i], -S[i] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalQR<T>::axpy(
    integer         N,
    valueType       alpha,
    valueType const L[],
    valueType const D[],
    valueType const U[],
    valueType const x[],
    valueType       beta,
    valueType       y[]
  ) const {
    tridiag_axpy( N, alpha, L, D, U, x, beta, y );
  }

  /*\
   *
   *  Solve
   *
   *  min || T x - b ||^2 + lambda ||x||^2
   *
   *
   *  / * * *       \
   *  |   * * *     |
   *  |     * * *   |
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
  TridiagonalQR<T>::lsq(
    integer nrhs,
    T       RHS[],
    integer ldRHS,
    T       lambda_in
  ) const {

    valueType lambda = normInfA * lambda_in;
    std::vector<T> D(nRC),
                   U(nRC-1),
                   U2(nRC-2),
                   tmp(nrhs);
    T CC, SS;

    for ( integer i = 0; i < nRC-1; ++i )
      rot( nrhs, RHS+i, ldRHS, RHS+i+1, ldRHS, C[i], S[i] );

    copy( nRC,   BD,  1, &D.front(),  1 );
    copy( nRC-1, BU,  1, &U.front(),  1 );
    copy( nRC-2, BU2, 1, &U2.front(), 1 );
    T line[3];
    integer i = 0;
    while ( i < nRC-1 ) {
      line[0] = line[2] = 0; line[1] = lambda;
      std::fill( tmp.begin(), tmp.end(), T(0) );
      integer j = i;
      while ( j < nRC-2 ) {
        line[0] = line[1];
        line[1] = line[2];
        line[2] = 0;
        rotg( D[j], line[0], CC, SS );
        rot( 1, &U[j],  1, &line[1], 1, CC, SS );
        rot( 1, &U2[j], 1, &line[2], 1, CC, SS );
        rot( nrhs, RHS+j, ldRHS, &tmp.front(), 1, CC, SS );
        ++j;
      }
      // penultima
      rotg( D[j], line[1], CC, SS );
      rot( 1, &U[j],  1, &line[2], 1, CC, SS );
      rot( nrhs, RHS+j, ldRHS, &tmp.front(), 1, CC, SS );
      ++j;

      rotg( D[j], line[2], CC, SS );
      rot( nrhs, RHS+j, ldRHS, &tmp.front(), 1, CC, SS );

      // ultima
      line[2] = lambda;
      rotg( D[j], line[2], CC, SS );
      rot( nrhs, RHS+j, ldRHS, &tmp.front(), 1, CC, SS );

      ++i; // next lambda
    }
    if ( nRC > 0 ) {
      integer j = nRC-1;
      line[0] = lambda;
      std::fill( tmp.begin(), tmp.end(), T(0) );
      rotg( D[j], line[0], CC, SS );
      rot( nrhs, RHS+j, ldRHS, &tmp.front(), 1, CC, SS );
    }

    for ( integer j = 0; j < nrhs; ++j ) {
      T * xb = RHS + j*ldRHS;
      xb[nRC-1] /= D[nRC-1];
      xb[nRC-2] = (xb[nRC-2]-U[nRC-2]*xb[nRC-1])/D[nRC-2];
      for ( integer k = nRC-3; k >= 0; --k )
        xb[k] = (xb[k]-U[k]*xb[k+1]-U2[k]*xb[k+2])/D[k];
    }
  }

}

///
/// eof: trid.cxx
///
