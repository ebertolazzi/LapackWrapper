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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::factorize(
    char const      who[],
    valueType const A[],
    integer         LDA
  ) {
    integer info = gecopy(
      this->nRows, this->nCols, A, LDA, this->Afactorized, this->nRows
    );
    LW_ASSERT(
      info == 0,
      "QR_no_alloc::factorize[{}] call lapack_wrapper::gecopy return info = {}\n",
      who, info
    );
    info = geqrf(
      this->nRows, this->nCols, this->Afactorized, this->nRows,
      this->Tau, this->WorkFactorized, this->Lwork
    );
    LW_ASSERT(
      info == 0,
      "QR_no_alloc::factorize[{}] call lapack_wrapper::geqrf return info = {}\n",
      who, info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  QR_no_alloc<T>::factorize( valueType const A[], integer LDA ) {
    integer info = gecopy(
      this->nRows, this->nCols, A, LDA, this->Afactorized, this->nRows
    );
    bool ok = info == 0;
    if ( ok ) {
      info = geqrf(
        this->nRows, this->nCols, this->Afactorized, this->nRows,
        this->Tau, this->WorkFactorized, this->Lwork
      );
      ok = info == 0;
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::applyQ(
    SideMultiply  SIDE,
    Transposition TRANS,
    integer       nRefl,
    integer       NR,
    integer       NC,
    valueType     C[],
    integer       ldC
  ) const {
    LW_ASSERT(
      (SIDE == lapack_wrapper::LEFT  && NR == this->nRows) ||
      (SIDE == lapack_wrapper::RIGHT && NC == this->nRows),
      "QR_no_alloc::applyQ NR = {} NC = {} nRow = {}\n",
      NR, NC, this->nRows
    );
    integer info = ormqr(
      SIDE, TRANS,
      NR, NC,
      nRefl,  // numero riflettori usati nel prodotto Q
      this->Afactorized, this->nRows /*ldA*/,
      this->Tau,
      C, ldC,
      this->WorkFactorized, this->Lwork
    );
    LW_ASSERT(
      info == 0,
      "QR_no_alloc::applyQ call lapack_wrapper::ormqr return info = {} Lwork = {}\n",
      info, this->Lwork
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::getR( valueType R[], integer ldR ) const {
    integer minRC = std::min( this->nRows, this->nCols );
    gezero( minRC, minRC, R, ldR );
    for ( integer i = 0; i < minRC; ++i )
      for ( integer j = i; j < minRC; ++j )
        R[i+j*ldR] = this->Afactorized[ i+j*this->nRows];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::getR( Matrix<T> & R ) const {
    integer minRC = std::min( this->nRows, this->nCols );
    R.setup( minRC, minRC );
    R.zero_fill();
    for ( integer i = 0; i < minRC; ++i )
      for ( integer j = i; j < minRC; ++j )
        R(i,j) = this->Afactorized[ i+j*this->nRows];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::getQ( valueType Q[], integer ldQ ) const {
    geid( this->nRows, this->nCols, Q, ldQ );
    Q_mul( this->nRows, this->nCols , Q, ldQ );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::getQ( Matrix<T> & Q ) const {
    Q.setup( this->nRows, this->nCols );
    Q.id( 1.0 );
    Q_mul( Q );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::getA( valueType A[], integer ldA ) const {
    gecopy(
      this->nRows, this->nCols, this->Afactorized, this->nRows, A, ldA
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::getA( Matrix<T> & A ) const {
    A.setup( this->nRows, this->nCols );
    gecopy(
      this->nRows, this->nCols, this->Afactorized, this->nRows,
      A.get_data(), A.lDim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::getTau( valueType tau[] ) const {
    copy( std::min(this->nRows, this->nCols), this->Tau, 1, tau, 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::getTau( Matrix<T> & tau ) const {
    tau.setup( 1, std::min(this->nRows, this->nCols) );
    copy( std::min(this->nRows, this->nCols), this->Tau, 1, tau.get_data(), 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::solve( valueType xb[] ) const {
    LW_ASSERT0(
      this->nRows == this->nCols,
      "QR_no_alloc::solve, factored matrix must be square\n"
    );
    Qt_mul(xb);
    invR_mul(xb);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::t_solve( valueType xb[] ) const {
    LW_ASSERT0(
      this->nRows == this->nCols,
      "QR_no_alloc::solve_t, factored matrix must be square\n"
    );
    invRt_mul(xb);
    Q_mul(xb);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::solve( integer nrhs, valueType XB[], integer ldXB ) const {
    LW_ASSERT0(
      this->nRows == this->nCols,
      "QR_no_alloc::solve, factored matrix must be square\n"
    );
    Qt_mul( this->nRows, nrhs, XB, ldXB );
    invR_mul( this->nRows, nrhs, XB, ldXB );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR_no_alloc<T>::t_solve( integer nrhs, valueType XB[], integer ldXB ) const {
    LW_ASSERT0(
      this->nRows == this->nCols,
      "QR_no_alloc::solve_t, factored matrix must be square\n"
    );
    invRt_mul( this->nRows, nrhs, XB, ldXB );
    Q_mul( this->nRows, nrhs, XB, ldXB );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR<T>::allocate( integer NR, integer NC ) {
    if ( this->nRows != NR || this->nCols != NC ) {
      valueType tmp; // get optimal allocation
      integer info = geqrf( NR, NC, nullptr, NR, nullptr, &tmp, -1 );
      LW_ASSERT(
        info == 0,
        "QR::allocate call lapack_wrapper::geqrf return info = {}\n", info
      );

      integer minRC = std::min( NR, NC );
      integer maxRC = std::max( NR, NC );
      integer L     = std::max( integer(tmp), maxRC );

      this->allocReals.allocate( size_t(NR*NC+L+minRC) );

      this->nRows          = NR;
      this->nCols          = NC;
      this->nReflector     = minRC;
      this->Lwork          = L;
      this->Afactorized    = this->allocReals( size_t(NR*NC) );
      this->WorkFactorized = this->allocReals( size_t(L) );
      this->Tau            = this->allocReals( size_t(minRC) );
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
  void
  QRP_no_alloc<T>::factorize(
    char const      who[],
    valueType const A[],
    integer         LDA
  ) {
    // calcolo fattorizzazione QR della matrice A
    integer info = gecopy(
      this->nRows, this->nCols, A, LDA, this->Afactorized, this->nRows
    );
    LW_ASSERT(
      info == 0,
      "QRP_no_alloc::factorize[{}] call lapack_wrapper::gecopy return info = {}\n",
      who, info
    );
    std::fill( this->JPVT, this->JPVT+this->nCols, integer(0) );
    info = geqp3(
      this->nRows, this->nCols, this->Afactorized, this->nRows,
      this->JPVT, this->Tau, this->WorkFactorized, this->Lwork
    );
    LW_ASSERT(
      info == 0,
      "QRP_no_alloc::factorize[{}] call lapack_wrapper::geqrf return info = {}\n",
      who, info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  QRP_no_alloc<T>::factorize( valueType const A[], integer LDA ) {
    // calcolo fattorizzazione QR della matrice A
    integer info = gecopy(
      this->nRows, this->nCols, A, LDA, this->Afactorized, this->nRows
    );
    bool ok = info == 0;
    if ( ok ) {
      std::fill( this->JPVT, this->JPVT+this->nCols, integer(0) );
      info = geqp3(
        this->nRows, this->nCols, this->Afactorized, this->nRows,
        this->JPVT, this->Tau, this->WorkFactorized, this->Lwork
      );
      ok = info == 0;
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP_no_alloc<T>::permute( valueType x[], integer incx ) const {
    // applico permutazione
    for ( integer i = 0; i < this->nCols; ++i )
      this->WorkPermute[this->JPVT[i]-1] = x[i*incx];
    copy( this->nCols, this->WorkPermute, 1, x, incx );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP_no_alloc<T>::inv_permute( valueType x[], integer incx ) const {
    // applico permutazione
    for ( integer i = 0; i < this->nCols; ++i )
      this->WorkPermute[i] = x[incx*(this->JPVT[i]-1)];
    copy( this->nCols, this->WorkPermute, 1, x, incx );
  }

  template <typename T>
  void
  QRP_no_alloc<T>::permute_rows(
    integer   nr,
    integer   nc,
    valueType C[],
    integer   ldC
  ) const {
    LW_ASSERT(
      nr == this->nCols,
      "QRP_no_alloc::permute_rows, bad number of row, expected {} found {}\n",
      this->nCols, nr
    );
    for ( integer j = 0; j < nc; ++j ) permute( C + ldC*j, 1 );
  }

  template <typename T>
  void
  QRP_no_alloc<T>::inv_permute_rows(
    integer   nr,
    integer   nc,
    valueType C[],
    integer   ldC
  ) const {
    LW_ASSERT(
      nr == this->nCols,
      "QRP_no_alloc::inv_permute_rows, bad number of row, expected {} found {}\n",
      this->nCols, nr
    );
    for ( integer j = 0; j < nc; ++j ) inv_permute( C + ldC*j, 1 );
  }

  template <typename T>
  void
  QRP_no_alloc<T>::permute_cols(
    integer   nr,
    integer   nc,
    valueType C[],
    integer   ldC
  ) const {
    LW_ASSERT(
      nc == this->nCols,
      "QRP_no_alloc::permute_cols, bad number of cols, expected {} found {}\n",
      this->nCols, nc
    );
    for ( integer i = 0; i < nr; ++i ) permute( C + i, ldC );
  }

  template <typename T>
  void
  QRP_no_alloc<T>::inv_permute_cols(
    integer   nr,
    integer   nc,
    valueType C[],
    integer   ldC
  ) const {
    LW_ASSERT(
      nc == this->nCols,
      "QRP_no_alloc::inv_permute_cols, bad number of cols, expected {} found {}\n",
      this->nCols, nc
    );
    for ( integer i = 0; i < nr; ++i ) inv_permute( C + i, ldC );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP_no_alloc<T>::solve( valueType xb[] ) const {
    LW_ASSERT0(
      this->nRows == this->nCols,
      "QRP_no_alloc::solve, factored matrix must be square\n"
    );
    Qt_mul(xb);
    invR_mul(xb);
    permute(xb); // da aggiungere!
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP_no_alloc<T>::t_solve( valueType xb[] ) const {
    LW_ASSERT0(
      this->nRows == this->nCols,
      "QRP_no_alloc::solve_t, factored matrix must be square\n"
    );
    inv_permute(xb); // da aggiungere!
    invRt_mul(xb);
    Q_mul(xb);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP_no_alloc<T>::solve( integer nrhs, valueType XB[], integer ldXB ) const {
    LW_ASSERT(
      this->nRows == this->nCols,
      "QRP_no_alloc::solve, factored matrix must be square, found {} x {}\n",
      this->nRows, this->nCols
    );
    Qt_mul( this->nRows, nrhs, XB, ldXB );
    invR_mul( this->nRows, nrhs, XB, ldXB );
    permute_rows( this->nRows, nrhs, XB, ldXB ); // da aggiungere!
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP_no_alloc<T>::t_solve( integer nrhs, valueType XB[], integer ldXB ) const {
    LW_ASSERT(
      this->nRows == this->nCols,
      "QRP_no_alloc::solve_t, factored matrix must be square, found {} x {}\n",
      this->nRows, this->nCols
    );
    inv_permute_rows( this->nRows, nrhs, XB, ldXB ); // da aggiungere!
    invRt_mul( this->nRows, nrhs, XB, ldXB );
    Q_mul( this->nRows, nrhs, XB, ldXB );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP<T>::allocate( integer NR, integer NC ) {
    if ( this->nRows != NR || this->nCols != NC ) {
      valueType tmp; // get optimal allocation
      integer info = geqp3( NR, NC, nullptr, NR, nullptr, nullptr, &tmp, -1 );
      LW_ASSERT(
        info == 0,
        "QRP::allocate call lapack_wrapper::geqp3 return info = {}\n", info
      );

      integer minRC = std::min( NR, NC );
      integer maxRC = std::max( NR, NC );
      integer L     = std::max( integer(tmp), maxRC );

      this->allocReals.allocate( size_t(NR*NC+L+minRC+maxRC) );
      this->allocIntegers.allocate( size_t(NC) );

      this->nRows          = NR;
      this->nCols          = NC;
      this->nReflector     = minRC;
      this->Lwork          = L;
      this->Afactorized    = this->allocReals( size_t(NR*NC) );
      this->WorkFactorized = this->allocReals( size_t(L) );
      this->Tau            = this->allocReals( size_t(minRC) );
      this->WorkPermute    = this->allocReals( size_t(maxRC) );
      this->JPVT           = this->allocIntegers( size_t(NC) );

    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP_no_alloc<T>::getPerm( integer jpvt[] ) const {
    std::copy( this->JPVT, this->JPVT+this->nCols, jpvt );
  }

}

///
/// eof: qr.cxx
///
