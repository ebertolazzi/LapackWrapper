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
  LU_no_alloc<T>::LU_no_alloc()
  : LinearSystemSolver<T>()
  , nRows(0)
  , nCols(0)
  , Afactorized(nullptr)
  , Work(nullptr)
  , Iwork(nullptr)
  , i_pivot(nullptr)
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU_no_alloc<T>::factorize(
    char const      who[],
    valueType const A[],
    integer         LDA
  ) {
    integer info = gecopy( this->nRows, this->nCols, A, LDA, this->Afactorized, this->nRows );
    LW_ASSERT(
      info == 0,
      "LU::factorize[{}] gecopy(nRow={}, nCol={}, A, LDA={}, B, LDB={}) : INFO = {}",
      who, this->nRows, this->nCols, LDA, this->nRows, info
    );
    info = getrf(
      this->nRows, this->nCols, this->Afactorized, this->nRows, this->i_pivot
    );
    LW_ASSERT(
      info == 0, "LU::factorize[{}] getrf INFO = {}", who, info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU_no_alloc<T>::check_ls( char const who[] ) const {
    LW_ASSERT(
      this->nRows == this->nCols,
      "LU<T>::{}, rectangular matrix {} x {}", who, this->nRows, this->nCols
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU_no_alloc<T>::solve( valueType xb[] ) const {
    check_ls("solve");
    integer info = getrs(
      NO_TRANSPOSE,
      this->nRows, 1, this->Afactorized, this->nRows, this->i_pivot,
      xb, this->nRows
    );
    LW_ASSERT( info == 0, "LU::solve, getrs INFO = {}", info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU_no_alloc<T>::t_solve( valueType xb[] ) const {
    check_ls("t_solve");
    integer info = getrs(
      TRANSPOSE,
      this->nRows, 1, this->Afactorized, this->nRows, this->i_pivot,
      xb, this->nRows
    );
    LW_ASSERT( info == 0, "LU::t_solve, getrs INFO = {}", info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU_no_alloc<T>::solve( integer nrhs, valueType B[], integer ldB ) const {
    check_ls("solve");
    integer info = getrs(
      NO_TRANSPOSE,
      this->nRows, nrhs, this->Afactorized, this->nRows, this->i_pivot,
      B, ldB
    );
    LW_ASSERT( info == 0, "LU::solve getrs INFO = {}", info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU_no_alloc<T>::t_solve( integer nrhs, valueType B[], integer ldB ) const {
    check_ls("t_solve");
    integer info = getrs(
      TRANSPOSE,
      this->nRows, nrhs, this->Afactorized, this->nRows, this->i_pivot,
      B, ldB
    );
    LW_ASSERT( info >= 0, "LU::t_solve getrs INFO = {}", info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  typename LU_no_alloc<T>::valueType
  LU_no_alloc<T>::cond1( valueType norm1 ) const {
    valueType rcond;
    integer info = gecon1(
      this->nRows, this->Afactorized, this->nRows,
      norm1, rcond, this->Work, this->Iwork
    );
    LW_ASSERT( info == 0, "LU::cond1, gecon1 return info = {}", info );
    return rcond;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  typename LU_no_alloc<T>::valueType
  LU_no_alloc<T>::condInf( valueType normInf ) const {
    valueType rcond;
    integer info = geconInf(
      this->nRows, this->Afactorized, this->nRows,
      normInf, rcond, this->Work, this->Iwork
    );
    LW_ASSERT( info == 0, "LU::condInf, geconInf return info = {}", info );
    return rcond;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  LU<T>::LU()
  : LU_no_alloc<T>()
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
    if ( this->nRows != NR || this->nCols != NC ) {
      this->nRows = NR;
      this->nCols = NC;
      integer NRC = NR*NC;
      integer NW  = 2*(NR+NC);
      allocReals.allocate( size_t(NRC+NW) );
      allocIntegers.allocate( size_t(2*NR) );

      this->no_allocate(
        NR, NC,
        allocReals( size_t(NRC) ),   // Afactorized
        allocIntegers( size_t(NR) ), // i_pivot
        allocReals( size_t(NW) ),    // Work
        allocIntegers( size_t(NR) )  // Iwork
      );

    }
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
  LUPQ_no_alloc<T>::LUPQ_no_alloc()
  : LinearSystemSolver<T>()
  , nRC(0)
  , Afactorized(nullptr)
  , i_piv(nullptr)
  , j_piv(nullptr)
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ_no_alloc<T>::factorize(
    char const      who[],
    valueType const A[],
    integer         LDA
  ) {
    integer info = gecopy(
      this->nRC, this->nRC, A, LDA, this->Afactorized, this->nRC
    );
    LW_ASSERT(
      info == 0, "LUPQ_no_alloc::factorize[{}] gecopy INFO = {}", who, info
    );
    info = getc2(
      this->nRC, this->Afactorized, this->nRC, this->i_piv, this->j_piv
    );
    LW_ASSERT(
      info == 0, "LUPQ_no_alloc::factorize[{}] getc2 INFO = {}", who, info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ_no_alloc<T>::solve( valueType xb[] ) const {
    // Apply permutations IPIV to RHS
    swaps( 1, xb, this->nRC, 0, this->nRC-2, this->i_piv, 1 );

    // Solve for L part
    trsv(
      LOWER, NO_TRANSPOSE, UNIT,
      this->nRC, this->Afactorized, this->nRC, xb, 1
    );

    // Solve for U part
    trsv(
      UPPER, NO_TRANSPOSE, NON_UNIT,
      this->nRC, this->Afactorized, this->nRC, xb, 1
    );

    // Apply permutations JPIV to the solution (RHS)
    swaps( 1, xb, this->nRC, 0, this->nRC-2, this->j_piv, -1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ_no_alloc<T>::solve( integer nrhs, valueType B[], integer ldB ) const {
    // Apply permutations IPIV to RHS
    swaps( nrhs, B, ldB, 0, this->nRC-2, this->i_piv, 1 );

    // Solve for L part
    trsm(
      LEFT, LOWER, NO_TRANSPOSE, UNIT,
      this->nRC, nrhs, 1.0, this->Afactorized, this->nRC, B, ldB
    );

    // Solve for U part
    trsm(
      LEFT, UPPER, NO_TRANSPOSE, NON_UNIT,
      this->nRC, nrhs, 1.0, this->Afactorized, this->nRC, B, ldB
    );

    // Apply permutations JPIV to the solution (RHS)
    swaps( nrhs, B, ldB, 0, this->nRC-2, this->j_piv, -1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ_no_alloc<T>::t_solve( valueType xb[] ) const {
    // Apply permutations JPIV to the solution (RHS)
    swaps( 1, xb, this->nRC, 0, this->nRC-2, this->j_piv, 1 );

    // Solve for U part
    trsv(
      UPPER, TRANSPOSE, NON_UNIT,
      this->nRC, this->Afactorized, this->nRC, xb, 1
    );

    // Solve for L part
    trsv(
      LOWER, TRANSPOSE, UNIT,
      this->nRC, this->Afactorized, this->nRC, xb, 1
    );

    // Apply permutations IPIV to RHS
    swaps( 1, xb, this->nRC, 0, this->nRC-2, this->i_piv, -1 );

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ_no_alloc<T>::t_solve( integer nrhs, valueType B[], integer ldB ) const {

    // Apply permutations JPIV to the solution (RHS)
    swaps( nrhs, B, ldB, 0, this->nRC-2, this->j_piv, 1 );

    // Solve for U part
    trsm(
      LEFT, UPPER, TRANSPOSE, NON_UNIT,
      this->nRC, nrhs, 1.0, this->Afactorized, this->nRC, B, ldB
    );

    // Solve for L part
    trsm(
      LEFT, LOWER, TRANSPOSE, UNIT,
      this->nRC, nrhs, 1.0, this->Afactorized, this->nRC, B, ldB
    );

    // Apply permutations IPIV to RHS
    swaps( nrhs, B, ldB, 0, this->nRC-2, this->i_piv, -1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  LUPQ<T>::LUPQ()
  : LUPQ_no_alloc<T>()
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
  LUPQ<T>::allocate( integer NRC ) {
    if ( this->nRC != NRC ) {
      allocReals.allocate( size_t(NRC*NRC) );
      allocIntegers.allocate( size_t(2*NRC) );
      this->no_allocate(
        NRC,
        allocReals( size_t(NRC*NRC) ), //_Afactorized,
        allocIntegers( size_t(NRC) ),  // _i_piv,
        allocIntegers( size_t(NRC) )   // _j_piv
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

}

///
/// eof: lu.cxx
///
