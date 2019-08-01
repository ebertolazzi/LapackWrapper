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
  void
  QR<T>::setMaxNrhs( integer mnrhs ) {
    LAPACK_WRAPPER_ASSERT(
      mnrhs > 0,
      "lapack_wrapper::QR::setMaxNrhs, maxNrhs = " << mnrhs
    );
    maxNrhs = mnrhs;
  }

  template <typename T>
  void
  QR<T>::allocate( integer NR, integer NC, integer Lwrk ) {
    nRow       = NR;
    nCol       = NC;
    nReflector = std::min(nRow,nCol);
    Lwork      = Lwrk;
    allocReals.allocate( size_t(nRow*nCol+Lwork+nReflector) );
    Amat = allocReals( size_t(nRow*nCol) );
    Work = allocReals( size_t(Lwork) );
    Tau  = allocReals( size_t(nReflector) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR<T>::allocate( integer NR, integer NC ) {
    if ( nRow != NR || nCol != NC || Lwork < maxNrhs ) {
      valueType tmp; // get optimal allocation
      integer info = geqrf( NR, NC, nullptr, NR, nullptr, &tmp, -1 );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "QR::allocate call lapack_wrapper::geqrf return info = " << info
      );
      integer L = std::max( integer(tmp), maxNrhs );
      if ( L < NR ) L = NR;
      if ( L < NC ) L = NC;
      allocate( NR, NC, L );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR<T>::applyQ(
    SideMultiply  SIDE,
    Transposition TRANS,
    integer       nRefl,
    integer       NR,
    integer       NC,
    valueType     C[],
    integer       ldC
  ) const {
    LAPACK_WRAPPER_ASSERT(
      (SIDE == lapack_wrapper::LEFT  && NR == nRow) ||
      (SIDE == lapack_wrapper::RIGHT && NC == nRow),
      "QR::applyQ NR = " << NR << " NC = " << NC << " nRow = " << nRow
    );
    integer info = ormqr(
      SIDE, TRANS,
      NR, NC,
      nRefl,  // numero riflettori usati nel prodotto Q
      Amat, nRow /*ldA*/,
      Tau,
      C, ldC,
      Work, Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "QR::applyQ call lapack_wrapper::ormqr return info = " << info <<
      " Lwork = " << Lwork
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR<T>::getR( valueType R[], integer ldR ) const {
    integer minRC = std::min( nRow, nCol );
    gezero( minRC, minRC, R, ldR );
    for ( integer i = 0; i < minRC; ++i )
      for ( integer j = i; j < minRC; ++j )
        R[i+j*ldR] = Amat[ i+j*nRow];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR<T>::solve( valueType xb[] ) const {
    LAPACK_WRAPPER_ASSERT(
      nRow == nCol,
      "in QR::solve, factored matrix must be square"
    );
    Qt_mul(xb);
    invR_mul(xb);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR<T>::t_solve( valueType xb[] ) const {
    LAPACK_WRAPPER_ASSERT(
      nRow == nCol,
      "in QR::solve_t, factored matrix must be square"
    );
    invRt_mul(xb);
    Q_mul(xb);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR<T>::solve( integer nrhs, valueType XB[], integer ldXB ) const {
    LAPACK_WRAPPER_ASSERT(
      nRow == nCol,
      "in QR::solve, factored matrix must be square"
    );
    Qt_mul( nRow, nrhs, XB, ldXB );
    invR_mul( nRow, nrhs, XB, ldXB );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR<T>::t_solve( integer nrhs, valueType XB[], integer ldXB ) const {
    LAPACK_WRAPPER_ASSERT(
      nRow == nCol,
      "in QR::solve_t, factored matrix must be square"
    );
    invRt_mul( nRow, nrhs, XB, ldXB );
    Q_mul( nRow, nrhs, XB, ldXB );
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
  QRP<T>::allocate( integer NR, integer NC ) {
    if ( nRow != NR || nCol != NC || Lwork < maxNrhs ) {
      valueType tmp; // get optimal allocation
      integer info = geqp3( NR, NC, nullptr, NR, nullptr, nullptr, &tmp, -1 );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "QRP::allocate call lapack_wrapper::geqp3 return info = " << info
      );
      integer L = std::max( integer(tmp), maxNrhs );
      if ( L < NR ) L = NR;
      if ( L < NC ) L = NC;
      QR<T>::allocate( NR, NC, L );
    }
    allocIntegers.allocate(size_t(NC));
    JPVT = allocIntegers(size_t(NC));
  }

  template <typename T>
  void
  QRP<T>::permute( valueType x[] ) const {
    // applico permutazione
    for ( integer i = 0; i < nCol; ++i ) Work[JPVT[i]-1] = x[i];
    copy( nCol, Work, 1, x, 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP<T>::inv_permute( valueType x[] ) const {
    // applico permutazione
    for ( integer i = 0; i < nCol; ++i ) Work[i] = x[JPVT[i]-1];
    copy( nCol, Work, 1, x, 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP<T>::solve( valueType xb[] ) const {
    LAPACK_WRAPPER_ASSERT(
      nRow == nCol,
      "in QRP::solve, factored matrix must be square"
    );
    Qt_mul(xb);
    invR_mul(xb);
    permute(xb); // da aggiungere!
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP<T>::t_solve( valueType xb[] ) const {
    LAPACK_WRAPPER_ASSERT(
      nRow == nCol,
      "in QRP::solve_t, factored matrix must be square"
    );
    inv_permute(xb); // da aggiungere!
    invRt_mul(xb);
    Q_mul(xb);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP<T>::solve( integer nrhs, valueType XB[], integer ldXB ) const {
    LAPACK_WRAPPER_ASSERT(
      nRow == nCol,
      "in QRP::solve, factored matrix must be square"
    );
    Qt_mul( nRow, nrhs, XB, ldXB );
    invR_mul( nRow, nrhs, XB, ldXB );
    permute_rows( nRow, nrhs, XB, ldXB ); // da aggiungere!
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP<T>::t_solve( integer nrhs, valueType XB[], integer ldXB ) const {
    LAPACK_WRAPPER_ASSERT(
      nRow == nCol,
      "in QRP::solve_t, factored matrix must be square"
    );
    inv_permute_rows( nRow, nrhs, XB, ldXB ); // da aggiungere!
    invRt_mul( nRow, nrhs, XB, ldXB );
    Q_mul( nRow, nrhs, XB, ldXB );
  }

}

///
/// eof: qr.cxx
///
