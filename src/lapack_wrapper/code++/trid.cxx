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
    if ( Utils::isZero(beta) ) {
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
    if ( m_nRC != N ) {
      m_allocReals.allocate(3*N);
      m_nRC  = N;
      m_L    = m_allocReals(N);
      m_D    = m_allocReals(N);
      m_WORK = m_allocReals(N);
    }
    copy( N, _L, 1, m_L, 1 );
    copy( N, _D, 1, m_D, 1 );
    integer info = pttrf( N, m_L, m_D );
    UTILS_ASSERT(
      info == 0,
      "TridiagonalSPD::factorize[{}], return info = {}\n", who, info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  TridiagonalSPD<T>::factorize(
    integer         N,
    valueType const _L[],
    valueType const _D[]
  ) {
    if ( m_nRC != N ) {
      m_allocReals.allocate(3*N);
      m_nRC  = N;
      m_L    = m_allocReals(N);
      m_D    = m_allocReals(N);
      m_WORK = m_allocReals(N);
    }
    copy( N, _L, 1, m_L, 1 );
    copy( N, _D, 1, m_D, 1 );
    integer info = pttrf( N, m_L, m_D );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  T
  TridiagonalSPD<T>::cond1( valueType norm1 ) const {
    valueType rcond;
    integer info = ptcon1( m_nRC, m_D, m_L, norm1, rcond, m_WORK );
    UTILS_ASSERT( info == 0, "TridiagonalSPD::cond1, return info = {}\n", info );
    return rcond;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  TridiagonalSPD<T>::solve( valueType xb[] ) const {
    integer info = pttrs( m_nRC, 1, m_D, m_L, xb, m_nRC );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  TridiagonalSPD<T>::t_solve( valueType xb[] ) const {
    integer info = pttrs( m_nRC, 1, m_D, m_L, xb, m_nRC );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  TridiagonalSPD<T>::solve( integer nrhs, valueType xb[], integer ldXB ) const {
    integer info = pttrs( m_nRC, nrhs, m_D, m_L, xb, ldXB );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  TridiagonalSPD<T>::t_solve( integer nrhs, valueType xb[], integer ldXB ) const {
    integer info = pttrs( m_nRC, nrhs, m_D, m_L, xb, ldXB );
    return info == 0;
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
    valueType const L[],
    valueType const D[],
    valueType const U[]
  ) {
    if ( m_nRC != N ) {
      m_nRC = N;
      m_allocReals.allocate(6*N);
      m_allocIntegers.allocate(2*N);
      m_L     = m_allocReals(N);
      m_D     = m_allocReals(N);
      m_U     = m_allocReals(N);
      m_U2    = m_allocReals(N);
      m_WORK  = m_allocReals(2*N);
      m_IPIV  = m_allocIntegers(N);
      m_IWORK = m_allocIntegers(N);
    }
    copy( N, L, 1, m_L, 1 );
    copy( N, D, 1, m_D, 1 );
    copy( N, U, 1, m_U, 1 );
    integer info = gttrf( N, m_L, m_D, m_U, m_U2, m_IPIV );
    UTILS_ASSERT(
      info == 0, "TridiagonalLU::factorize[{}], return info = {}\n", who, info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  TridiagonalLU<T>::factorize(
    integer         N,
    valueType const L[],
    valueType const D[],
    valueType const U[]
  ) {
    if ( m_nRC != N ) {
      m_nRC = N;
      m_allocReals.allocate(6*N);
      m_allocIntegers.allocate(2*N);
      m_L     = m_allocReals(N);
      m_D     = m_allocReals(N);
      m_U     = m_allocReals(N);
      m_U2    = m_allocReals(N);
      m_WORK  = m_allocReals(2*N);
      m_IPIV  = m_allocIntegers(N);
      m_IWORK = m_allocIntegers(N);
    }
    copy( N, L, 1, m_L, 1 );
    copy( N, D, 1, m_D, 1 );
    copy( N, U, 1, m_U, 1 );
    integer info = gttrf( N, m_L, m_D, m_U, m_U2, m_IPIV );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  T
  TridiagonalLU<T>::cond1( valueType norm1 ) const {
    valueType rcond;
    integer info = gtcon1(
      m_nRC, m_L, m_D, m_U, m_U2, m_IPIV, norm1, rcond, m_WORK, m_IWORK
    );
    UTILS_ASSERT( info == 0, "TridiagonalLU::cond1, return info = {}\n", info );
    return rcond;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  T
  TridiagonalLU<T>::condInf( valueType normInf ) const {
    valueType rcond;
    integer info = gtconInf(
      m_nRC, m_L, m_D, m_U, m_U2, m_IPIV, normInf, rcond, m_WORK, m_IWORK
    );
    UTILS_ASSERT( info == 0, "TridiagonalLU::cond1, return info = {}\n", info );
    return rcond;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  TridiagonalLU<T>::solve( valueType xb[] ) const {
    integer info = gttrs(
      NO_TRANSPOSE, m_nRC, 1, m_L, m_D, m_U, m_U2, m_IPIV, xb, m_nRC
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalLU<T>::solve( char const who[], valueType xb[] ) const {
    integer info = gttrs(
      NO_TRANSPOSE, m_nRC, 1, m_L, m_D, m_U, m_U2, m_IPIV, xb, m_nRC
    );
    UTILS_ASSERT(
      info == 0,
      "TridiagonalLU::solve, return info = {}\nat {}\n",
      info, who
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  TridiagonalLU<T>::t_solve( valueType xb[] ) const {
    integer info = gttrs(
      TRANSPOSE, m_nRC, 1, m_L, m_D, m_U, m_U2, m_IPIV, xb, m_nRC
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalLU<T>::t_solve( char const who[], valueType xb[] ) const {
    integer info = gttrs(
      TRANSPOSE, m_nRC, 1, m_L, m_D, m_U, m_U2, m_IPIV, xb, m_nRC
    );
    UTILS_ASSERT(
      info == 0,
      "TridiagonalLU::t_solve, return info = {}\nat {}\n",
      info, who
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  TridiagonalLU<T>::solve( integer nrhs, valueType xb[], integer ldXB ) const {
    integer info = gttrs(
      NO_TRANSPOSE, m_nRC, nrhs, m_L, m_D, m_U, m_U2, m_IPIV, xb, ldXB
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalLU<T>::solve(
    char const who[],
    integer    nrhs,
    valueType  xb[],
    integer    ldXB
  ) const {
    integer info = gttrs(
      NO_TRANSPOSE, m_nRC, nrhs, m_L, m_D, m_U, m_U2, m_IPIV, xb, ldXB
    );
    UTILS_ASSERT(
      info == 0,
      "TridiagonalLU::solve, return info = {}\nat {}\n",
      info, who
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  TridiagonalLU<T>::t_solve( integer nrhs, valueType xb[], integer ldXB ) const {
    integer info = gttrs(
      TRANSPOSE, m_nRC, nrhs, m_L, m_D, m_U, m_U2, m_IPIV, xb, ldXB
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalLU<T>::t_solve(
    char const who[],
    integer    nrhs,
    valueType  xb[],
    integer    ldXB
  ) const {
    integer info = gttrs(
      TRANSPOSE, m_nRC, nrhs, m_L, m_D, m_U, m_U2, m_IPIV, xb, ldXB
    );
    UTILS_ASSERT(
      info == 0,
      "TridiagonalLU::t_solve, return info = {}\nat {}\n",
      info, who
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalLU<T>::axpy(
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
  bool
  TridiagonalQR<T>::factorize(
    integer         N,
    valueType const L[],
    valueType const D[],
    valueType const U[]
  ) {
    m_allocReals.allocate(size_t(5*(N-1)));
    m_nRC = N;
    m_C   = m_allocReals(size_t(N-1));
    m_S   = m_allocReals(size_t(N-1));
    m_BD  = m_allocReals(size_t(N));
    m_BU  = m_allocReals(size_t(N-1));
    m_BU2 = m_allocReals(size_t(N-2));

    /*\
      | d u       | d u @     | d u @     | d u @     | d u @     |
      | l d u     | 0 d u     | 0 d u @   | 0 d u @   | 0 d u @   |
      |   l d u   |   l d u   |   0 d u   |   0 d u @ |   0 d u @ |
      |     l d u |     l d u |     l d u |     0 d u |     0 d u |
      |       l d |       l d |       l d |       l d |       0 d |
    \*/

    lapack_wrapper::copy( N,   D, 1, m_BD, 1 );
    lapack_wrapper::copy( N-1, U, 1, m_BU, 1 );
    lapack_wrapper::zero( N-2, m_BU2, 1 );

    m_normInfA = 0;
    integer i = 0;
    for (; i < N-2; ++i ) {
      valueType Li = L[i];
      rotg( m_BD[i], Li, m_C[i], m_S[i] );
      rot( 1, &m_BU[i],  1, &m_BD[i+1], 1, m_C[i], m_S[i] );
      rot( 1, &m_BU2[i], 1, &m_BU[i+1], 1, m_C[i], m_S[i] );
      valueType sum = std::abs(m_BD[i]) + std::abs(m_BU[i]) + std::abs(m_BU2[i]);
      if ( sum > m_normInfA ) m_normInfA = sum;
    }
    valueType Li = L[i];
    rotg( m_BD[i], Li, m_C[i], m_S[i] );
    rot( 1, &m_BU[i], 1, &m_BD[i+1], 1, m_C[i], m_S[i] );

    valueType sum = std::abs(m_BD[i]) + std::abs(m_BU[i]);
    if ( sum > m_normInfA ) m_normInfA = sum;
    sum = std::abs(m_BD[i+1]);
    if ( sum > m_normInfA ) m_normInfA = sum;

    // Q A = R
    return true;
  }

  template <typename T>
  void
  TridiagonalQR<T>::factorize(
    char const      who[],
    integer         N,
    valueType const L[],
    valueType const D[],
    valueType const U[]
  ) {
    bool ok = this->factorize( N, L, D, U );
    UTILS_ASSERT( ok, "TridiagonalQR<T>::factorize failed at {}\n", who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalQR<T>::Rsolve( valueType xb[] ) const {
    xb[m_nRC-1] /= m_BD[m_nRC-1];
    xb[m_nRC-2] = (xb[m_nRC-2]-m_BU[m_nRC-2]*xb[m_nRC-1])/m_BD[m_nRC-2];
    for ( integer i = m_nRC-3; i >= 0; --i )
      xb[i] = (xb[i]-m_BU[i]*xb[i+1]-m_BU2[i]*xb[i+2])/m_BD[i];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalQR<T>::RsolveTransposed( valueType xb[] ) const {
    xb[0] /= m_BD[0];
    xb[1] = (xb[1]-m_BU[0]*xb[0])/m_BD[1];
    for ( integer i = 2; i < m_nRC; ++i )
      xb[i] = (xb[i]-m_BU[i]*xb[i-1]-m_BU2[i]*xb[i-2])/m_BD[i];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  TridiagonalQR<T>::solve( valueType xb[] ) const {
    // A x = b --> Q A x = Q b --> R x = Q b
    // applico Q b
    for ( integer i = 0; i < m_nRC-1; ++i )
      rot( 1, &xb[i], 1, &xb[i+1], 1, m_C[i], m_S[i] );
    Rsolve( xb );
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  TridiagonalQR<T>::t_solve( valueType xb[] ) const {
    // A^T x = b --> A^T Q^T Q x = b --> R^T Q x = b --> R^T y = b  x = Q^T y
    RsolveTransposed( xb );
    // applico Q^T b
    for ( integer i = m_nRC-2; i >= 0; --i )
      rot( 1, &xb[i], 1, &xb[i+1], 1, m_C[i], -m_S[i] );
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  TridiagonalQR<T>::solve( integer nrhs, valueType xb[], integer ldXB ) const {
    // A x = b --> Q A x = Q b --> R x = Q b
    // applico Q b
    for ( integer i = 0; i < m_nRC-1; ++i )
      rot( nrhs, &xb[i], ldXB, &xb[i+1], ldXB, m_C[i], m_S[i] );
    for ( integer i = 0; i < nrhs; ++i )
      Rsolve( xb+i*ldXB );
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  TridiagonalQR<T>::t_solve( integer nrhs, valueType xb[], integer ldXB ) const {
    // A^T x = b --> A^T Q^T Q x = b --> R^T Q x = b --> R^T y = b  x = Q^T y
    for ( integer i = 0; i < nrhs; ++i )
      RsolveTransposed(xb+i*ldXB);
    for ( integer i = m_nRC-2; i >= 0; --i )
      rot( nrhs, &xb[i], ldXB, &xb[i+1], ldXB, m_C[i], -m_S[i] );
    return true;
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

    valueType lambda = m_normInfA * lambda_in;
    std::vector<T> D(m_nRC),
                   U(m_nRC-1),
                   U2(m_nRC-2),
                   tmp(nrhs);
    T CC, SS;

    for ( integer i = 0; i < m_nRC-1; ++i )
      rot( nrhs, RHS+i, ldRHS, RHS+i+1, ldRHS, m_C[i], m_S[i] );

    copy( m_nRC,   m_BD,  1, &D.front(),  1 );
    copy( m_nRC-1, m_BU,  1, &U.front(),  1 );
    copy( m_nRC-2, m_BU2, 1, &U2.front(), 1 );
    T line[3];
    integer i = 0;
    while ( i < m_nRC-1 ) {
      line[0] = line[2] = 0; line[1] = lambda;
      std::fill( tmp.begin(), tmp.end(), T(0) );
      integer j = i;
      while ( j < m_nRC-2 ) {
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
    if ( m_nRC > 0 ) {
      integer j = m_nRC-1;
      line[0] = lambda;
      std::fill( tmp.begin(), tmp.end(), T(0) );
      rotg( D[j], line[0], CC, SS );
      rot( nrhs, RHS+j, ldRHS, &tmp.front(), 1, CC, SS );
    }

    for ( integer j = 0; j < nrhs; ++j ) {
      T * xb = RHS + j*ldRHS;
      xb[m_nRC-1] /= D[m_nRC-1];
      xb[m_nRC-2] = (xb[m_nRC-2]-U[m_nRC-2]*xb[m_nRC-1])/D[m_nRC-2];
      for ( integer k = m_nRC-3; k >= 0; --k )
        xb[k] = (xb[k]-U[k]*xb[k+1]-U2[k]*xb[k+2])/D[k];
    }
  }

}

///
/// eof: trid.cxx
///
