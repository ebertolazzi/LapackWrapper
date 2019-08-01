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
/// file: ls.cxx
///

namespace lapack_wrapper {

  /*\
   |   _     ____ ____
   |  | |   / ___/ ___|
   |  | |   \___ \___ \
   |  | |___ ___) |__) |
   |  |_____|____/____/
   |
  \*/

  template <typename T>
  void
  LSS<T>::setMaxNrhs( integer mnrhs ) {
    LAPACK_WRAPPER_ASSERT(
      mnrhs > 0,
      "LSS::setMaxNrhs, maxNrhs = " << mnrhs
    );
    maxNrhs         = mnrhs;
    maxNrhs_changed = true;
  }

  template <typename T>
  void
  LSS<T>::allocate( integer NR, integer NC ) {

    if ( nRow != NR || nCol != NC || maxNrhs_changed ) {
      nRow = NR;
      nCol = NC;
      valueType tmp;
      integer info = gelss(
        NR, NC, maxNrhs,
        nullptr, NR, nullptr, NR, nullptr,
        rcond, rank, &tmp, -1
      );
      LAPACK_WRAPPER_ASSERT(
        info == 0, "LSS::allocate, in gelss info = " << info
      );

      Lwork = integer(tmp);
      if ( NR != NC ) {
        info = gelss(
          NC, NR, maxNrhs,
          nullptr, NC, nullptr, NC, nullptr,
          rcond, rank, &tmp, -1
        );
        LAPACK_WRAPPER_ASSERT(
          info == 0, "LSS::allocate, in gelss info = " << info
        );
        if ( Lwork < integer(tmp) ) Lwork = integer(tmp);
      }

      integer minRC = std::min(NR,NC);

      allocReals.allocate( size_t(2*NR*NC+Lwork+minRC) );
      Amat     = allocReals( size_t(2*nRow*nCol) );
      Work     = allocReals( size_t(Lwork) );
      sigma    = allocReals( size_t(minRC) );
      AmatWork = Amat+nRow*nCol;

      maxNrhs_changed = false;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSS<T>::solve( valueType xb[] ) const {
    // save matrix
    copy( nRow*nCol, Amat, 1, AmatWork, 1);
    integer info = gelss(
      nRow, nCol, 1,
      AmatWork, nRow,
      xb, nRow,
      sigma, rcond, rank,
      Work, Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0, "LSS::solve (rhs=1), in gelss info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSS<T>::t_solve( valueType xb[] ) const {
    // save matrix
    for ( integer i = 0; i < nCol; ++i )
      copy( nRow, Amat+i*nRow, 1, AmatWork+i, nCol );
    integer info = gelss(
      nCol, nRow, 1,
      AmatWork, nCol,
      xb, nCol,
      sigma, rcond, rank,
      Work, Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0, "LSS::t_solve (rhs=1), in gelss info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSS<T>::solve(
    integer   nrhs,
    valueType B[],
    integer   ldB
  ) const {
    // save matrix
    copy( nRow*nCol, Amat, 1, AmatWork, 1 );
    integer info = gelss(
      nRow, nCol, nrhs,
      AmatWork, nRow,
      B, ldB,
      sigma, rcond, rank,
      Work, Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "LSS::solve (rhs=" << nrhs << "), in gelss info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSS<T>::t_solve(
    integer   nrhs,
    valueType B[],
    integer   ldB
  ) const {
    // save matrix
    for ( integer i = 0; i < nCol; ++i )
      copy( nRow, Amat+i*nRow, 1, AmatWork+i, nCol );
    integer info = gelss(
      nCol, nRow, nrhs,
      AmatWork, nCol,
      B, ldB,
      sigma, rcond, rank,
      Work, Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "LSS::t_solve (rhs=" << nrhs << "), in gelss info = " << info
    );
  }

  /*\
   |  _     ______   __
   | | |   / ___\ \ / /
   | | |   \___ \\ V /
   | | |___ ___) || |
   | |_____|____/ |_|
   |
  \*/

  template <typename T>
  void
  LSY<T>::setMaxNrhs( integer mnrhs ) {
    LAPACK_WRAPPER_ASSERT(
      mnrhs > 0, "LSY::setMaxNrhs, maxNrhs = " << mnrhs
    );
    maxNrhs         = mnrhs;
    maxNrhs_changed = true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY<T>::allocate( integer NR, integer NC ) {

    if ( nRow != NR || nCol != NC || maxNrhs_changed ) {
      nRow = NR;
      nCol = NC;
      valueType tmp;
      integer info = gelsy(
        NR, NC, maxNrhs, nullptr, NR, nullptr, NR, nullptr,
        rcond, rank, &tmp, -1
      );
      LAPACK_WRAPPER_ASSERT(
        info == 0, "LSY::allocate, in gelss info = " << info
      );

      Lwork = integer(tmp);
      if ( NR != NC ) {
        info = gelsy(
          NC, NR, maxNrhs, nullptr, NC, nullptr, NC, nullptr,
          rcond, rank, &tmp, -1
        );
        LAPACK_WRAPPER_ASSERT(
          info == 0, "LSY::allocate, in gelss info = " << info
        );
        if ( Lwork < integer(tmp) ) Lwork = integer(tmp);
      }

      allocReals.allocate( size_t(2*NR*NC+Lwork) );
      Amat     = allocReals( size_t(2*nRow*nCol) );
      Work     = allocReals( size_t(Lwork) );
      AmatWork = Amat+nRow*nCol;
      allocInts.allocate( size_t(NC) );
      jpvt = allocInts( size_t(NC) );

      maxNrhs_changed = false;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY<T>::solve( valueType xb[] ) const {
    // save matrix
    copy( nRow*nCol, Amat, 1, AmatWork, 1);
    integer info = gelsy(
      nRow, nCol, 1,
      AmatWork, nRow,
      xb, nRow, jpvt,
      rcond, rank,
      Work, Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0, "LSS::solve (rhs=1), in gelss info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY<T>::t_solve( valueType xb[] ) const {
    // save matrix
    for ( integer i = 0; i < nCol; ++i )
      copy( nRow, Amat+i*nRow, 1, AmatWork+i, nCol );
    integer info = gelsy(
      nCol, nRow, 1,
      AmatWork, nCol,
      xb, nCol, jpvt,
      rcond, rank,
      Work, Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0, "LSS::t_solve (rhs=1), in gelss info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY<T>::solve(
    integer   nrhs,
    valueType B[],
    integer   ldB
  ) const {
    // save matrix
    copy( nRow*nCol, Amat, 1, AmatWork, 1 );
    integer info = gelsy(
      nRow, nCol, nrhs,
      AmatWork, nRow,
      B, ldB, jpvt,
      rcond, rank,
      Work, Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "LSD::solve (rhs=" << nrhs << "), in gelsd info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY<T>::t_solve(
    integer   nrhs,
    valueType B[],
    integer   ldB
  ) const {
    // save matrix
    for ( integer i = 0; i < nCol; ++i )
      copy( nRow, Amat+i*nRow, 1, AmatWork+i, nCol );
    integer info = gelsy(
      nCol, nRow, nrhs,
      AmatWork, nCol,
      B, ldB, jpvt,
      rcond, rank,
      Work, Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "LSD::t_solve (rhs=" << nrhs << "), in gelsd info = " << info
    );
  }

}

///
/// eof: ls.cxx
///
