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
/// file: qr.cxx
///

namespace lapack_wrapper {

  /*\
   |    ___  ____
   |   / _ \|  _ \
   |  | | | | |_) |
   |  | |_| |  _ <
   |   \__\_\_| \_\
   |
  \*/

  template <typename T>
  integer
  QR_no_alloc<T>::get_Lwork_QR( integer NR, integer NC ) const {
    real_type tmp; // get optimal allocation
    integer info = geqrf( NR, NC, nullptr, NR, nullptr, &tmp, -1 );
    UTILS_ASSERT(
      info == 0,
      "QR_no_alloc::get_Lwork_QR call lapack_wrapper::geqrf return info = {}\n", info
    );
    return std::max( integer(tmp), std::max( NR, NC ) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::no_allocate(
    integer     NR,
    integer     NC,
    integer     Lwork,
    real_type * Work
  ) {
    m_nrows      = NR;
    m_ncols      = NC;
    m_nReflector = std::min(NR,NC);
    m_LworkQR    = this->get_Lwork_QR( NR, NC );
    integer Lwmin = m_LworkQR + NR*NC + m_nReflector;
    UTILS_ASSERT(
      Lwork >= Lwmin,
      "QR_no_alloc::no_allocate( NR = {}, NC = {}, Lwork {}, .. )\n"
      "Lwork must be >= {}\n",
      NR, NC, Lwork, Lwmin
    );
    real_type * ptr = Work;
    m_Afactorized = ptr; ptr += NR*NC;
    m_Tau         = ptr; ptr += m_nReflector;
    m_WorkQR      = ptr;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::factorize_nodim(
    string_view     who,
    real_type const A[],
    integer         LDA
  ) {
    integer info = gecopy( m_nrows, m_ncols, A, LDA, m_Afactorized, m_nrows );
    UTILS_ASSERT(
      info == 0,
      "QR_no_alloc::factorize[{}] call lapack_wrapper::gecopy return info = {}\n",
      who, info
    );
    info = geqrf( m_nrows, m_ncols, m_Afactorized, m_nrows, m_Tau, m_WorkQR, m_LworkQR );
    UTILS_ASSERT(
      info == 0,
      "QR_no_alloc::factorize[{}] call lapack_wrapper::geqrf return info = {}\n",
      who, info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  QR_no_alloc<T>::factorize_nodim( real_type const A[], integer LDA ) {
    integer info = gecopy( m_nrows, m_ncols, A, LDA, m_Afactorized, m_nrows );
    bool ok = info == 0;
    if ( ok ) {
      info = geqrf( m_nrows, m_ncols, m_Afactorized, m_nrows, m_Tau, m_WorkQR, m_LworkQR );
      ok = info == 0;
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::applyQ(
    SideMultiply  SIDE,
    Transposition TRANS,
    integer       nRefl,
    integer       NR,
    integer       NC,
    real_type     C[],
    integer       ldC
  ) const {
    UTILS_ASSERT(
      (SIDE == SideMultiply::LEFT  && NR == m_nrows) ||
      (SIDE == SideMultiply::RIGHT && NC == m_nrows),
      "QR_no_alloc::applyQ( SIDE = {}, TRANS = {}, Nrefl = {}, NR = {}, NC = {}, C, ldC = {})\n"
      "original matrix is {} x {}\n",
      to_blas(SIDE),
      to_blas(TRANS),
      nRefl, NR, NC, ldC,
      m_nrows, m_ncols
    );
    integer info = ormqr(
      SIDE, TRANS,
      NR, NC,
      nRefl,  // numero riflettori usati nel prodotto Q
      m_Afactorized, m_nrows /*ldA*/,
      m_Tau,
      C, ldC,
      m_WorkQR, m_LworkQR
    );
    UTILS_ASSERT(
      info == 0,
      "QR_no_alloc::applyQ call lapack_wrapper::ormqr return info = {} Lwork = {}\n",
      info, m_LworkQR
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::Q_mul( real_type x[] ) const {
    applyQ(
      SideMultiply::LEFT,
      Transposition::NO,
      m_nReflector, m_nrows, 1, x, m_nrows
    );
  }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::Qt_mul( real_type x[] ) const {
    applyQ(
      SideMultiply::LEFT,
      Transposition::YES,
      m_nReflector, m_nrows, 1, x, m_nrows
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::Q_mul(
    integer   nr,
    integer   nc,
    real_type C[],
    integer   ldC
  ) const {
    applyQ(
      SideMultiply::LEFT,
      Transposition::NO,
      m_nReflector, nr, nc, C, ldC
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::Qt_mul(
    integer   nr,
    integer   nc,
    real_type C[],
    integer   ldC
  ) const {
    applyQ(
      SideMultiply::LEFT,
      Transposition::YES,
      m_nReflector, nr, nc, C, ldC
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::mul_Q(
    integer   nr,
    integer   nc,
    real_type C[],
    integer   ldC
  ) const {
    applyQ(
      SideMultiply::RIGHT,
      Transposition::NO,
      m_nReflector, nr, nc, C, ldC
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::mul_Qt(
    integer   nr,
    integer   nc,
    real_type C[],
    integer   ldC
  ) const {
    applyQ(
      SideMultiply::RIGHT,
      Transposition::YES,
      m_nReflector, nr, nc, C, ldC
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::invR_mul( real_type x[], integer incx ) const {
    trsv( // m_nrows = leading dimension
      ULselect::UPPER,
      Transposition::NO,
      DiagonalType::NON_UNIT,
      m_nReflector, m_Afactorized, m_nrows, x, incx
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::invR_mul( integer nr, integer nc, real_type C[], integer ldC ) const {
    trsm(
      SideMultiply::LEFT,
      ULselect::UPPER,
      Transposition::NO,
      DiagonalType::NON_UNIT,
      nr, nc, 1.0, m_Afactorized, m_nrows, C, ldC
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::invRt_mul( real_type x[], integer incx ) const {
    trsv( // m_nrows = leading dimension
      ULselect::UPPER,
      Transposition::YES,
      DiagonalType::NON_UNIT,
      m_nReflector, m_Afactorized, m_nrows, x, incx
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::invRt_mul(
    integer   nr,
    integer   nc,
    real_type C[],
    integer   ldC
  ) const {
    trsm(
      SideMultiply::LEFT,
      ULselect::UPPER,
      Transposition::YES,
      DiagonalType::NON_UNIT,
      nr, nc, 1.0, m_Afactorized, m_nrows, C, ldC
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::mul_invR(
    integer nr,
    integer nc,
    real_type C[],
    integer ldC
  ) const {
    trsm(
      SideMultiply::RIGHT,
      ULselect::UPPER,
      Transposition::NO,
      DiagonalType::NON_UNIT,
      nr, nc, 1.0, m_Afactorized, m_nrows, C, ldC
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::mul_invRt( integer nr, integer nc, real_type C[], integer ldC ) const {
    trsm(
      SideMultiply::RIGHT,
      ULselect::UPPER,
      Transposition::YES,
      DiagonalType::NON_UNIT,
      nr, nc, 1.0, m_Afactorized, m_nrows, C, ldC
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::getR( real_type R[], integer ldR ) const {
    integer minRC = std::min( m_nrows, m_ncols );
    gezero( minRC, minRC, R, ldR );
    for ( integer i{0}; i < minRC; ++i )
      for ( integer j{i}; j < minRC; ++j )
        R[i+j*ldR] = m_Afactorized[ i+j*m_nrows ];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::getRt( real_type R[], integer ldR ) const {
    integer minRC = std::min( m_nrows, m_ncols );
    gezero( minRC, minRC, R, ldR );
    for ( integer i{0}; i < minRC; ++i )
      for ( integer j{i}; j < minRC; ++j )
        R[j+i*ldR] = m_Afactorized[ i+j*m_nrows ];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::getR( Matrix<T> & R ) const {
    integer minRC = std::min( m_nrows, m_ncols );
    R.setup( minRC, minRC );
    R.zero_fill();
    for ( integer i{0}; i < minRC; ++i )
      for ( integer j{i}; j < minRC; ++j )
        R(i,j) = m_Afactorized[ i+j*m_nrows ];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::getRt( Matrix<T> & R ) const {
    integer minRC = std::min( m_nrows, m_ncols );
    R.setup( minRC, minRC );
    R.zero_fill();
    for ( integer i{0}; i < minRC; ++i )
      for ( integer j{i}; j < minRC; ++j )
        R(j,i) = m_Afactorized[ i+j*m_nrows ];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::getQ( real_type Q[], integer ldQ ) const {
    geid( m_nrows, m_nrows, Q, ldQ );
    Q_mul( m_nrows, m_nrows, Q, ldQ );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::getQ( Matrix<T> & Q ) const {
    Q.setup( m_nrows, m_nrows );
    Q.id( 1.0 );
    Q_mul( Q );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::getQreduced( real_type Q[], integer ldQ ) const {
    geid( m_nrows, m_ncols, Q, ldQ );
    Q_mul( m_nrows, m_ncols, Q, ldQ );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::getQreduced( Matrix<T> & Q ) const {
    Q.setup( m_nrows, m_ncols );
    Q.id( 1.0 );
    Q_mul( Q );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::getA( real_type A[], integer ldA ) const {
    gecopy(
      m_nrows, m_ncols, m_Afactorized, m_nrows, A, ldA
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::getA( Matrix<T> & A ) const {
    A.setup( m_nrows, m_ncols );
    gecopy(
      m_nrows, m_ncols, m_Afactorized, m_nrows,
      A.data(), A.ldim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::getTau( real_type tau[] ) const {
    copy( std::min( m_nrows, m_ncols ), m_Tau, 1, tau, 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::getTau( Matrix<T> & tau ) const {
    tau.setup( 1, std::min( m_nrows, m_ncols) );
    copy( std::min( m_nrows, m_ncols ), m_Tau, 1, tau.data(), 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  QR_no_alloc<T>::solve( real_type xb[] ) const {
    if ( m_nrows != m_ncols ) return false;
    Qt_mul(xb);
    invR_mul(xb);
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  QR_no_alloc<T>::t_solve( real_type xb[] ) const {
    if ( m_nrows != m_ncols ) return false;
    invRt_mul(xb);
    Q_mul(xb);
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  QR_no_alloc<T>::solve( integer nrhs, real_type XB[], integer ldXB ) const {
    if ( m_nrows != m_ncols ) return false;
    Qt_mul( m_nrows, nrhs, XB, ldXB );
    invR_mul( m_nrows, nrhs, XB, ldXB );
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  QR_no_alloc<T>::t_solve( integer nrhs, real_type XB[], integer ldXB ) const {
    if ( m_nrows != m_ncols ) return false;
    invRt_mul( m_nrows, nrhs, XB, ldXB );
    Q_mul( m_nrows, nrhs, XB, ldXB );
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR<T>::allocate( integer NR, integer NC ) {
    if ( m_nrows != NR || m_ncols != NC ) {
      integer Lwork = NR*NC + std::min( NR, NC ) + this->get_Lwork_QR( NR, NC );
      this->no_allocate( NR, NC, Lwork, m_allocReals.realloc( size_t(Lwork) ) );
    }
  }

  /*\
   |    ___  ____  ____
   |   / _ \|  _ \|  _ \
   |  | | | | |_) | |_) |
   |  | |_| |  _ <|  __/
   |   \__\_\_| \_\_|
   |
  \*/

  template <typename T>
  integer
  QRP_no_alloc<T>::get_Lwork_QRP( integer NR, integer NC ) const {
    real_type tmp; // get optimal allocation
    integer info = geqp3( NR, NC, nullptr, NR, nullptr, nullptr, &tmp, -1 );
    UTILS_ASSERT(
      info == 0,
      "QRP::get_Lwork_QRP call lapack_wrapper::geqp3 return info = {}\n", info
    );

    integer maxRC = std::max( NR, NC );
    integer L     = std::max( integer(tmp), maxRC );
    return L;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP_no_alloc<T>::no_allocate(
    integer     NR,
    integer     NC,
    integer     Lwork,
    real_type * Work,
    integer     Liwork,
    integer   * iWork
  ) {

    m_LworkQR = this->get_Lwork_QRP( NR, NC );

    integer Lwmin = m_LworkQR + NR*NC + NR+NC;

    UTILS_ASSERT(
      Liwork >= NC && Lwork >= Lwmin,
      "QRP_no_alloc::no_allocate( NR = {}, NC = {}, Lwork = {},..., Liwork = {}, ...)\n"
      "Lwork must be >= {} amd Liwork = NC",
      NR, NC, Lwork, Liwork, Lwmin
    );

    m_nrows      = NR;
    m_ncols      = NC;
    m_nReflector = std::min(NR,NC);

    m_JPVT = iWork;

    real_type * ptr = Work;
    m_WorkQR      = ptr; ptr += m_LworkQR;
    m_WorkPermute = ptr; ptr += std::max(NR,NC);
    m_Afactorized = ptr; ptr += NR*NC;
    m_Tau         = ptr; //ptr += m_nReflector;

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP_no_alloc<T>::factorize_nodim(
    string_view     who,
    real_type const A[],
    integer         LDA
  ) {
    // calcolo fattorizzazione QR della matrice A
    integer info = gecopy( m_nrows, m_ncols, A, LDA, m_Afactorized, m_nrows );
    UTILS_ASSERT(
      info == 0,
      "QRP_no_alloc::factorize[{}] call lapack_wrapper::gecopy return info = {}\n",
      who, info
    );
    std::fill_n( m_JPVT, m_ncols, integer(0) );
    info = geqp3(
      m_nrows, m_ncols, m_Afactorized, m_nrows,
      m_JPVT, m_Tau, m_WorkQR, m_LworkQR
    );
    UTILS_ASSERT(
      info == 0,
      "QRP_no_alloc::factorize[{}] call lapack_wrapper::geqrf return info = {}\n",
      who, info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  QRP_no_alloc<T>::factorize_nodim( real_type const A[], integer LDA ) {
    // calcolo fattorizzazione QR della matrice A
    integer info = gecopy( m_nrows, m_ncols, A, LDA, m_Afactorized, m_nrows );
    bool ok = info == 0;
    if ( ok ) {
      std::fill_n( m_JPVT, m_ncols, integer(0) );
      info = geqp3(
        m_nrows, m_ncols, m_Afactorized, m_nrows,
        m_JPVT, m_Tau, m_WorkQR, m_LworkQR
      );
      ok = info == 0;
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP_no_alloc<T>::permute( real_type x[], integer incx ) const {
    // applico permutazione
    for ( integer i{0}; i < m_ncols; ++i )
      m_WorkPermute[m_JPVT[i]-1] = x[i*incx];
    copy( m_ncols, m_WorkPermute, 1, x, incx );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP_no_alloc<T>::inv_permute( real_type x[], integer incx ) const {
    // applico permutazione
    for ( integer i{0}; i < m_ncols; ++i )
      m_WorkPermute[i] = x[incx*(m_JPVT[i]-1)];
    copy( m_ncols, m_WorkPermute, 1, x, incx );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP_no_alloc<T>::permute_rows(
    integer   nr,
    integer   nc,
    real_type C[],
    integer   ldC
  ) const {
    UTILS_ASSERT(
      nr == m_ncols,
      "QRP_no_alloc::permute_rows, bad number of row, expected {} found {}\n",
      m_ncols, nr
    );
    for ( integer j = 0; j < nc; ++j ) permute( C + ldC*j, 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP_no_alloc<T>::inv_permute_rows(
    integer   nr,
    integer   nc,
    real_type C[],
    integer   ldC
  ) const {
    UTILS_ASSERT(
      nr == m_ncols,
      "QRP_no_alloc::inv_permute_rows, bad number of row, expected {} found {}\n",
      m_ncols, nr
    );
    for ( integer j{0}; j < nc; ++j ) inv_permute( C + ldC*j, 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP_no_alloc<T>::permute_cols(
    integer   nr,
    integer   nc,
    real_type C[],
    integer   ldC
  ) const {
    UTILS_ASSERT(
      nc == m_ncols,
      "QRP_no_alloc::permute_cols, bad number of cols, expected {} found {}\n",
      m_ncols, nc
    );
    for ( integer i{0}; i < nr; ++i ) permute( C + i, ldC );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP_no_alloc<T>::inv_permute_cols(
    integer   nr,
    integer   nc,
    real_type C[],
    integer   ldC
  ) const {
    UTILS_ASSERT(
      nc == m_ncols,
      "QRP_no_alloc::inv_permute_cols, bad number of cols, expected {} found {}\n",
      m_ncols, nc
    );
    for ( integer i{0}; i < nr; ++i ) inv_permute( C + i, ldC );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  QRP_no_alloc<T>::solve( real_type xb[] ) const {
    if ( m_nrows != m_ncols ) return false;
    Qt_mul(xb);
    invR_mul(xb);
    permute(xb); // da aggiungere!
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  QRP_no_alloc<T>::t_solve( real_type xb[] ) const {
    if ( m_nrows != m_ncols ) return false;
    inv_permute(xb); // da aggiungere!
    invRt_mul(xb);
    Q_mul(xb);
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  QRP_no_alloc<T>::solve( integer nrhs, real_type XB[], integer ldXB ) const {
    if ( m_nrows != m_ncols ) return false;
    Qt_mul( m_nrows, nrhs, XB, ldXB );
    invR_mul( m_nrows, nrhs, XB, ldXB );
    permute_rows( m_nrows, nrhs, XB, ldXB ); // da aggiungere!
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  QRP_no_alloc<T>::t_solve( integer nrhs, real_type XB[], integer ldXB ) const {
    if ( m_nrows != m_ncols ) return false;
    inv_permute_rows( m_nrows, nrhs, XB, ldXB ); // da aggiungere!
    invRt_mul( m_nrows, nrhs, XB, ldXB );
    Q_mul( m_nrows, nrhs, XB, ldXB );
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP<T>::allocate( integer NR, integer NC ) {
    if ( m_nrows != NR || m_ncols != NC ) {
      integer Lwork = this->get_Lwork_QRP( NR, NC ) + NR*NC + NR+NC;
      this->no_allocate(
        NR, NC,
        Lwork, m_allocReals.realloc( size_t(Lwork) ),
        NC,    m_allocIntegers.realloc( size_t(NC) )
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP_no_alloc<T>::getPerm( integer jpvt[] ) const {
    std::copy_n( m_JPVT, m_ncols, jpvt );
  }

}

///
/// eof: qr.cxx
///
