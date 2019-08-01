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
/// file: lu.cxx
///

namespace lapack_wrapper {

  /*\
   |   _    _   _
   |  | |  | | | |
   |  | |  | | | |
   |  | |__| |_| |
   |  |_____\___/
   |
  \*/

  template <typename T>
  LU<T>::LU()
  : Factorization<T>()
  , allocReals("allocReals")
  , allocIntegers("allocIntegers")
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  LU<T>::~LU() {
    allocReals.free();
    allocIntegers.free();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU<T>::allocate( integer NR, integer NC ) {
    if ( nRow != NR || nCol != NC ) {
      nRow = NR;
      nCol = NC;
      allocReals.allocate( size_t(nRow*nCol+2*(nRow+nCol)) );
      allocIntegers.allocate( size_t(2*nRow) );
      Amat    = allocReals( size_t(nRow*nCol) );
      Work    = allocReals( size_t(2*(nRow+nCol)) );
      i_pivot = allocIntegers( size_t(nRow) );
      Iwork   = allocIntegers( size_t(nRow) );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU<T>::factorize( char const who[] ) {
    integer info = getrf( nRow, nCol, Amat, nRow, i_pivot );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "LU::factorize[" << who << "] getrf INFO = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU<T>::factorize(
    char const      who[],
    integer         NR,
    integer         NC,
    valueType const A[],
    integer LDA
  ) {
    allocate( NR, NC );
    integer info = gecopy( nRow, nCol, A, LDA, Amat, nRow );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "LU::factorize[" << who << "] gecopy INFO = " << info
    );
    factorize( who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU<T>::check_ls( char const who[] ) const {
    LAPACK_WRAPPER_ASSERT(
      nRow == nCol,
      "LU<T>::" << who << ", rectangular matrix " <<
      nRow << " x " << nCol
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU<T>::solve( valueType xb[] ) const {
    check_ls("solve");
    integer info = getrs(
      NO_TRANSPOSE, nRow, 1, Amat, nRow, i_pivot, xb, nRow
    );
    LAPACK_WRAPPER_ASSERT( info == 0, "LU::solve getrs INFO = " << info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU<T>::t_solve( valueType xb[] ) const {
    check_ls("t_solve");
    integer info = getrs(
      TRANSPOSE, nRow, 1, Amat, nRow, i_pivot, xb, nRow
    );
    LAPACK_WRAPPER_ASSERT( info == 0, "LU::t_solve getrs INFO = " << info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU<T>::solve( integer nrhs, valueType B[], integer ldB ) const {
    check_ls("solve");
    integer info = getrs(
      NO_TRANSPOSE, nRow, nrhs, Amat, nRow, i_pivot, B, ldB
    );
    LAPACK_WRAPPER_ASSERT( info == 0, "LU::solve getrs INFO = " << info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU<T>::t_solve( integer nrhs, valueType B[], integer ldB ) const {
    check_ls("t_solve");
    integer info = getrs(
      TRANSPOSE, nRow, nrhs, Amat, nRow, i_pivot, B, ldB
    );
    LAPACK_WRAPPER_ASSERT( info >= 0, "LU::t_solve getrs INFO = " << info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  typename LU<T>::valueType
  LU<T>::cond1( valueType norm1 ) const {
    valueType rcond;
    integer info = gecon1(
      nRow, Amat, nRow, norm1, rcond, Work, Iwork
    );
    LAPACK_WRAPPER_ASSERT( info == 0, "LU::cond1, gecon1 return info = " << info );
    return rcond;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  typename LU<T>::valueType
  LU<T>::condInf( valueType normInf ) const {
    valueType rcond;
    integer info = geconInf(
      nRow, Amat, nRow, normInf, rcond, Work, Iwork
    );
    LAPACK_WRAPPER_ASSERT( info == 0, "LU::condInf, geconInf return info = " << info );
    return rcond;
  }

  /*\
   |   _    _   _ ____   ___
   |  | |  | | | |  _ \ / _ \
   |  | |  | | | | |_) | | | |
   |  | |__| |_| |  __/| |_| |
   |  |_____\___/|_|    \__\_\
   |
  \*/

  template <typename T>
  LUPQ<T>::LUPQ()
  : Factorization<T>()
  , allocReals("LUPQ-allocReals")
  , allocIntegers("LUPQ-allocIntegers")
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  LUPQ<T>::~LUPQ() {
    allocReals.free();
    allocIntegers.free();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ<T>::allocate( integer NR, integer NC ) {
    LAPACK_WRAPPER_ASSERT(
      NR == NC,
      "LUPQ<T>::allocate, cannot allocate rectangular matrix " <<
      NR << " x " << NC
    );
    if ( nRow != NR || nCol != NC ) {
      nRow = NR;
      nCol = NC;
      allocReals.allocate( size_t(nRow*nCol) );
      Amat = allocReals( size_t(nRow*nCol) );
      allocIntegers.allocate( size_t(2*nRow) );
      ipiv = allocIntegers( size_t(nRow) );
      jpiv = allocIntegers( size_t(nRow) );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ<T>::check_ls( char const who[] ) const {
    LAPACK_WRAPPER_ASSERT(
      nRow == nCol,
      "LUPQ<T>::" << who << ", rectangular matrix " <<
      nRow << " x " << nCol
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ<T>::factorize( char const who[] ) {
    integer info = getc2( nRow, Amat, nRow, ipiv, jpiv );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "LUPQ::factorize[" << who << "] getrf INFO = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ<T>::factorize(
    char const      who[],
    integer         NR,
    integer         NC,
    valueType const A[],
    integer         LDA
  ) {
    LAPACK_WRAPPER_ASSERT(
      NR == NC,
      "LUPQ<T>::factorize[" << who <<
      "], cannot factorize rectangular matrix " <<
      NR << " x " << NC
    );
    allocate( NR, NC );
    integer info = gecopy( nRow, nCol, A, LDA, Amat, nRow );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "LUPQ::factorize[" << who << "] gecopy INFO = " << info
    );
    factorize( who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ<T>::solve( valueType xb[] ) const {
    // Apply permutations IPIV to RHS
    swaps( 1, xb, nRow, 0, nRow-2, ipiv, 1 );

    // Solve for L part
    trsv( LOWER, NO_TRANSPOSE, UNIT, nRow, Amat, nRow, xb, 1 );

    // Solve for U part
    trsv( UPPER, NO_TRANSPOSE, NON_UNIT, nRow, Amat, nRow, xb, 1 );

    // Apply permutations JPIV to the solution (RHS)
    swaps( 1, xb, nRow, 0, nRow-2, jpiv, -1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ<T>::solve( integer nrhs, valueType B[], integer ldB ) const {
    // Apply permutations IPIV to RHS
    swaps( nrhs, B, ldB, 0, nRow-2, ipiv, 1 );

    // Solve for L part
    trsm( LEFT, LOWER, NO_TRANSPOSE, UNIT, nRow, nrhs, 1.0, Amat, nRow, B, ldB );

    // Solve for U part
    trsm( LEFT, UPPER, NO_TRANSPOSE, NON_UNIT, nRow, nrhs, 1.0, Amat, nRow, B, ldB );

    // Apply permutations JPIV to the solution (RHS)
    swaps( nrhs, B, ldB, 0, nRow-2, jpiv, -1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ<T>::t_solve( valueType xb[] ) const {
    // Apply permutations JPIV to the solution (RHS)
    swaps( 1, xb, nRow, 0, nRow-2, jpiv, 1 );

    // Solve for U part
    trsv( UPPER, TRANSPOSE, NON_UNIT, nRow, Amat, nRow, xb, 1 );

    // Solve for L part
    trsv( LOWER, TRANSPOSE, UNIT, nRow, Amat, nRow, xb, 1 );

    // Apply permutations IPIV to RHS
    swaps( 1, xb, nRow, 0, nRow-2, ipiv, -1 );

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ<T>::t_solve( integer nrhs, valueType B[], integer ldB ) const {

    // Apply permutations JPIV to the solution (RHS)
    swaps( nrhs, B, ldB, 0, nRow-2, jpiv, 1 );

    // Solve for U part
    trsm( LEFT, UPPER, TRANSPOSE, NON_UNIT, nRow, nrhs, 1.0, Amat, nRow, B, ldB );

    // Solve for L part
    trsm( LEFT, LOWER, TRANSPOSE, UNIT, nRow, nrhs, 1.0, Amat, nRow, B, ldB );

    // Apply permutations IPIV to RHS
    swaps( nrhs, B, ldB, 0, nRow-2, ipiv, -1 );
  }

}

///
/// eof: lu.cxx
///
