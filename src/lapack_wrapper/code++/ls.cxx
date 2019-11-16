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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  integer
  LSS_no_alloc<T>::getL( integer NR, integer NC, integer nrhs ) const {
    valueType tmp;
    integer info = gelss(
      NR, NC, nrhs,
      nullptr, NR, nullptr, NR, nullptr,
      this->rcond, this->rank, &tmp, -1
    );
    LW_ASSERT( info == 0, "LSS_no_alloc::getL, in gelss info = {}\n", info );
    return integer(tmp);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSS_no_alloc<T>::factorize(
    char const      who[],
    valueType const A[],
    integer         LDA
  ) {
    integer info = gecopy(
      this->nRows, this->nCols, A, LDA, this->Amat, this->nRows
    );
    LW_ASSERT(
      info == 0,
      "LSS_no_alloc::factorize( who={},...) call lapack_wrapper::gecopy return info = {}\n",
      who, info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSS_no_alloc<T>::solve( valueType xb[] ) const {
    // check memory
    integer L = getL( this->nRows, this->nCols, 1 );
    if ( L > this->Lwork ) {
      allocWork.allocate( L );
      this->Lwork = L;
      this->Work  = allocWork( L );
    }
    // save matrix
    copy( this->nRows*this->nCols, this->Amat, 1, this->AmatWork, 1);
    integer info = gelss(
      this->nRows, this->nCols, 1,
      this->AmatWork, this->nRows,
      xb, this->nRows,
      this->sigma, this->rcond, this->rank,
      this->Work, this->Lwork
    );
    LW_ASSERT( info == 0, "LSS::solve (rhs=1), in gelss info = {}\n", info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSS_no_alloc<T>::t_solve( valueType xb[] ) const {
    // check memory
    integer L = getL( this->nCols, this->nRows, 1 );
    if ( L > this->Lwork ) {
      allocWork.allocate( L );
      this->Lwork = L;
      this->Work  = allocWork( L );
    }
    // save matrix
    for ( integer i = 0; i < this->nCols; ++i )
      copy( this->nRows, this->Amat+i*this->nRows, 1, this->AmatWork+i, this->nCols );
    integer info = gelss(
      this->nCols, this->nRows, 1,
      this->AmatWork, this->nCols,
      xb, this->nCols,
      this->sigma, this->rcond, this->rank,
      this->Work, this->Lwork
    );
    LW_ASSERT( info == 0, "LSS::t_solve (rhs=1), in gelss info = {}\n", info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSS_no_alloc<T>::solve(
    integer   nrhs,
    valueType B[],
    integer   ldB
  ) const {
    // check memory
    integer L = getL( this->nRows, this->nCols, nrhs );
    if ( L > this->Lwork ) {
      allocWork.allocate( L );
      this->Lwork = L;
      this->Work  = allocWork( L );
    }
    // save matrix
    copy( this->nRows*this->nCols, this->Amat, 1, this->AmatWork, 1 );
    integer info = gelss(
      this->nRows, this->nCols, nrhs,
      this->AmatWork, this->nRows,
      B, ldB,
      this->sigma, this->rcond, this->rank,
      this->Work, this->Lwork
    );
    LW_ASSERT( info == 0, "LSS::solve (rhs={}), in gelss info = {}\n", nrhs, info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSS_no_alloc<T>::t_solve(
    integer   nrhs,
    valueType B[],
    integer   ldB
  ) const {
    // check memory
    integer L = getL( this->nCols, this->nRows, nrhs );
    if ( L > this->Lwork ) {
      allocWork.allocate( L );
      this->Lwork = L;
      this->Work  = allocWork( L );
    }
    // save matrix
    for ( integer i = 0; i < this->nCols; ++i )
      copy( this->nRows, this->Amat+i*this->nRows, 1, this->AmatWork+i, this->nCols );
    integer info = gelss(
      this->nCols, this->nRows, nrhs,
      this->AmatWork, this->nCols,
      B, ldB,
      this->sigma, this->rcond, this->rank,
      this->Work, this->Lwork
    );
    LW_ASSERT(
      info == 0,"LSS::t_solve (rhs={}), in gelss info = {}\n", nrhs, info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSS<T>::allocate( integer NR, integer NC ) {

    if ( this->nRows != NR || this->nCols != NC ) {
      this->nRows   = NR;
      this->nCols   = NC;

      integer minRC = std::min(NR,NC);
      integer NRC   = NR*NC;
      allocReals.allocate( size_t(2*NRC+minRC) );

      this->no_allocate(
        NR, NC,
        allocReals( size_t(NRC) ),
        allocReals( size_t(minRC) ),
        allocReals( size_t(NRC) )
      );
    }
  }

  /*\
   |  _     ______   __
   | | |   / ___\ \ / /
   | | |   \___ \\ V /
   | | |___ ___) || |
   | |_____|____/ |_|
   |
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  integer
  LSY_no_alloc<T>::getL( integer NR, integer NC, integer nrhs ) const {
    valueType tmp;
    integer info = gelsy(
      NR, NC, nrhs,
      nullptr, NR, nullptr, NR, nullptr,
      this->rcond, this->rank, &tmp, -1
    );
    LW_ASSERT( info == 0, "LSY_no_alloc::getL, in gelss info = {}\n", info );
    return integer(tmp);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY_no_alloc<T>::factorize(
    char const      who[],
    valueType const A[],
    integer         LDA
  ) {
    integer info = gecopy(
      this->nRows, this->nCols, A, LDA, this->Amat, this->nRows
    );
    LW_ASSERT(
      info == 0,
      "LSY::factorize[{}] call lapack_wrapper::gecopy return info = {}\n",
      who, info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY_no_alloc<T>::solve( valueType xb[] ) const {
    // check memory
    integer L = getL( this->nRows, this->nCols, 1 );
    if ( L > this->Lwork ) {
      allocWork.allocate( L );
      this->Lwork = L;
      this->Work  = allocWork( L );
    }
    // save matrix
    copy( this->nRows*this->nCols, this->Amat, 1, this->AmatWork, 1 );
    integer info = gelsy(
      this->nRows, this->nCols, 1,
      this->AmatWork, this->nRows,
      xb, this->nRows, this->jpvt,
      this->rcond, this->rank,
      this->Work, this->Lwork
    );
    LW_ASSERT( info == 0, "LSS::solve (rhs=1), in gelsy info = {}\n", info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY_no_alloc<T>::t_solve( valueType xb[] ) const {
    // check memory
    integer L = getL( this->nCols, this->nRows, 1 );
    if ( L > this->Lwork ) {
      allocWork.allocate( L );
      this->Lwork = L;
      this->Work  = allocWork( L );
    }
    // save matrix
    for ( integer i = 0; i < this->nCols; ++i )
      copy( this->nRows, this->Amat+i*this->nRows, 1, this->AmatWork+i, this->nCols );
    integer info = gelsy(
      this->nCols, this->nRows, 1,
      this->AmatWork, this->nCols,
      xb, this->nCols, this->jpvt,
      this->rcond, this->rank,
      this->Work, this->Lwork
    );
    LW_ASSERT( info == 0, "LSS::t_solve (rhs=1), in gelsy info = {}\n", info  );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY_no_alloc<T>::solve(
    integer   nrhs,
    valueType B[],
    integer   ldB
  ) const {
    // check memory
    integer L = getL( this->nRows, this->nCols, nrhs );
    if ( L > this->Lwork ) {
      allocWork.allocate( L );
      this->Lwork = L;
      this->Work  = allocWork( L );
    }
    // save matrix
    copy( this->nRows*this->nCols, this->Amat, 1, this->AmatWork, 1 );
    integer info = gelsy(
      this->nRows, this->nCols, nrhs,
      this->AmatWork, this->nRows,
      B, ldB, this->jpvt,
      this->rcond, this->rank,
      this->Work, this->Lwork
    );
    LW_ASSERT(
      info == 0, "LSD::solve (rhs={}), in gelsy info = {}\n", nrhs, info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY_no_alloc<T>::t_solve(
    integer   nrhs,
    valueType B[],
    integer   ldB
  ) const {
    // check memory
    integer L = getL( this->nCols, this->nRows, nrhs );
    if ( L > this->Lwork ) {
      allocWork.allocate( L );
      this->Lwork = L;
      this->Work  = allocWork( L );
    }
    // save matrix
    for ( integer i = 0; i < this->nCols; ++i )
      copy( this->nRows, this->Amat+i*this->nRows, 1, this->AmatWork+i, this->nCols );
    integer info = gelsy(
      this->nCols, this->nRows, nrhs,
      this->AmatWork, this->nCols,
      B, ldB, this->jpvt,
      this->rcond, this->rank,
      this->Work, this->Lwork
    );
    LW_ASSERT(
      info == 0, "LSD::t_solve (rhs={}), in gelsy info = {}\n", nrhs, info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY<T>::allocate( integer NR, integer NC ) {

    if ( this->nRows != NR || this->nCols != NC ) {
      this->nRows = NR;
      this->nCols = NC;
      integer NRC = NR*NC;

      allocReals.allocate( size_t(2*NRC) );
      allocInts.allocate( size_t(NC) );

      this->no_allocate(
        NR, NC,
        allocReals( size_t(NRC) ),
        allocReals( size_t(NRC) ),
        allocInts( size_t(NC) )
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY<T>::factorize(
    char const      who[],
    integer         NR,
    integer         NC,
    valueType const A[],
    integer         LDA
  ) {
    allocate( NR, NC );
    integer info = gecopy( this->nRows, this->nCols, A, LDA, this->Amat, this->nRows );
    LW_ASSERT(
      info == 0,
      "LSY::factorize[{}] call lapack_wrapper::gecopy return info = {}\n",
      who, info
    );
  }

}

///
/// eof: ls.cxx
///
