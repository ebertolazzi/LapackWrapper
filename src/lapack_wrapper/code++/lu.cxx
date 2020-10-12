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
  , m_nRows(0)
  , m_nCols(0)
  , m_Afactorized(nullptr)
  , m_Work(nullptr)
  , m_Iwork(nullptr)
  , m_i_pivot(nullptr)
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU_no_alloc<T>::factorize(
    char const      who[],
    valueType const A[],
    integer         LDA
  ) {
    integer info = gecopy(
      m_nRows, m_nCols, A, LDA, m_Afactorized, m_nRows
    );
    LW_ASSERT(
      info == 0,
      "LU::factorize[{}] gecopy(nRow={}, nCol={}, A, LDA={}, B, LDB={}) : INFO = {}\n",
      who, m_nRows, m_nCols, LDA, m_nRows, info
    );
    info = getrf(
      m_nRows, m_nCols, m_Afactorized, m_nRows, m_i_pivot
    );
    LW_ASSERT(
      info == 0, "LU::factorize[{}] getrf INFO = {}\n", who, info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LU_no_alloc<T>::factorize( valueType const A[], integer LDA ) {
    integer info = gecopy(
      m_nRows, m_nCols, A, LDA, m_Afactorized, m_nRows
    );
    bool ok = info == 0;
    if ( ok ) {
      info = getrf(
        m_nRows, m_nCols, m_Afactorized, m_nRows, m_i_pivot
      );
      ok = info == 0;
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU_no_alloc<T>::check_ls( char const who[] ) const {
    LW_ASSERT(
      m_nRows == m_nCols,
      "LU<T>::{}, rectangular matrix {} x {}\n", who, m_nRows, m_nCols
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LU_no_alloc<T>::solve( valueType xb[] ) const {
    if ( m_nRows != m_nCols ) return false;
    integer info = getrs(
      NO_TRANSPOSE,
      m_nRows, 1, m_Afactorized, m_nRows, m_i_pivot,
      xb, m_nRows
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU_no_alloc<T>::solve( char const who[], valueType xb[] ) const {
    check_ls("solve");
    integer info = getrs(
      NO_TRANSPOSE,
      m_nRows, 1, m_Afactorized, m_nRows, m_i_pivot,
      xb, m_nRows
    );
    LW_ASSERT( info == 0, "LU::solve, getrs INFO = {}\nat {}\n", info, who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LU_no_alloc<T>::t_solve( valueType xb[] ) const {
    if ( m_nRows != m_nCols ) return false;
    integer info = getrs(
      TRANSPOSE,
      m_nRows, 1, m_Afactorized, m_nRows, m_i_pivot,
      xb, m_nRows
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU_no_alloc<T>::t_solve( char const who[], valueType xb[] ) const {
    check_ls( who );
    integer info = getrs(
      TRANSPOSE,
      m_nRows, 1, m_Afactorized, m_nRows, m_i_pivot,
      xb, m_nRows
    );
    LW_ASSERT( info == 0, "LU::t_solve, getrs INFO = {}\nat {}\n", info, who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LU_no_alloc<T>::solve( integer nrhs, valueType B[], integer ldB ) const {
    if ( m_nRows != m_nCols ) return false;
    integer info = getrs(
      NO_TRANSPOSE,
      m_nRows, nrhs, m_Afactorized, m_nRows, m_i_pivot,
      B, ldB
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU_no_alloc<T>::solve(
    char const who[],
    integer    nrhs,
    valueType  B[],
    integer    ldB
  ) const {
    check_ls(who);
    integer info = getrs(
      NO_TRANSPOSE,
      m_nRows, nrhs, m_Afactorized, m_nRows, m_i_pivot,
      B, ldB
    );
    LW_ASSERT( info == 0, "LU::solve getrs INFO = {}\nat {}\n", info, who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LU_no_alloc<T>::t_solve(
    integer    nrhs,
    valueType  B[],
    integer    ldB
  ) const {
    if ( m_nRows != m_nCols ) return false;
    integer info = getrs(
      TRANSPOSE,
      m_nRows, nrhs, m_Afactorized, m_nRows, m_i_pivot,
      B, ldB
    );
    return info >= 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU_no_alloc<T>::t_solve(
    char const who[],
    integer    nrhs,
    valueType  B[],
    integer    ldB
  ) const {
    check_ls( who );
    integer info = getrs(
      TRANSPOSE,
      m_nRows, nrhs, m_Afactorized, m_nRows, m_i_pivot,
      B, ldB
    );
    LW_ASSERT( info >= 0, "LU::t_solve getrs INFO = {}\nat {}\n", info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  typename LU_no_alloc<T>::valueType
  LU_no_alloc<T>::cond1( valueType norm1 ) const {
    valueType rcond;
    integer info = gecon1(
      m_nRows, m_Afactorized, m_nRows,
      norm1, rcond, m_Work, m_Iwork
    );
    LW_ASSERT( info == 0, "LU::cond1, gecon1 return info = {}\n", info );
    return rcond;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  typename LU_no_alloc<T>::valueType
  LU_no_alloc<T>::condInf( valueType normInf ) const {
    valueType rcond;
    integer info = geconInf(
      m_nRows, m_Afactorized, m_nRows,
      normInf, rcond, m_Work, m_Iwork
    );
    LW_ASSERT( info == 0, "LU::condInf, geconInf return info = {}\n", info );
    return rcond;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  LU<T>::LU()
  : LU_no_alloc<T>()
  , m_allocReals("allocReals")
  , m_allocIntegers("allocIntegers")
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  LU<T>::~LU() {
    m_allocReals.free();
    m_allocIntegers.free();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU<T>::allocate( integer NR, integer NC ) {
    if ( m_nRows != NR || m_nCols != NC ) {
      m_nRows = NR;
      m_nCols = NC;
      integer NRC = NR*NC;
      integer NW  = 2*(NR+NC);
      m_allocReals.allocate( size_t(NRC+NW) );
      m_allocIntegers.allocate( size_t(2*NR) );

      this->no_allocate(
        NR, NC,
        m_allocReals( size_t(NRC) ),   // Afactorized
        m_allocIntegers( size_t(NR) ), // i_pivot
        m_allocReals( size_t(NW) ),    // Work
        m_allocIntegers( size_t(NR) )  // Iwork
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
  , m_nRC(0)
  , m_Afactorized(nullptr)
  , m_i_piv(nullptr)
  , m_j_piv(nullptr)
  {}

  template <typename T>
  void
  LUPQ_no_alloc<T>::no_allocate(
    integer     NRC,
    integer     Lwork,
    valueType * Work,
    integer     Liwork,
    integer   * iWork
  ) {
    LW_ASSERT(
      Lwork >= NRC*NRC && Liwork >= 2*NRC,
      "LUPQ_no_alloc::no_allocate( NRC = {}, Lwork = {}, ..., Liwork = {}, ... )\n"
      "Lwork must be >= {} and Liwork >= {}\n",
      NRC, Lwork, Liwork, NRC*NRC, 2*NRC
    );
    m_nRC         = NRC;
    m_Afactorized = Work;
    m_i_piv       = iWork;
    m_j_piv       = iWork+NRC;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ_no_alloc<T>::factorize(
    char const      who[],
    valueType const A[],
    integer         LDA
  ) {
    integer info = gecopy(
      m_nRC, m_nRC, A, LDA, m_Afactorized, m_nRC
    );
    LW_ASSERT(
      info == 0, "LUPQ_no_alloc::factorize[{}] gecopy INFO = {}\n", who, info
    );
    info = getc2(
      m_nRC, m_Afactorized, m_nRC, m_i_piv, m_j_piv
    );
    LW_ASSERT(
      info == 0, "LUPQ_no_alloc::factorize[{}] getc2 INFO = {}\n", who, info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LUPQ_no_alloc<T>::factorize( valueType const A[], integer LDA ) {
    integer info = gecopy(
      m_nRC, m_nRC, A, LDA, m_Afactorized, m_nRC
    );
    bool ok = info == 0;
    if ( ok ) {
      info = getc2(
        m_nRC, m_Afactorized, m_nRC, m_i_piv, m_j_piv
      );
      ok = info == 0;
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LUPQ_no_alloc<T>::solve( valueType xb[] ) const {
    // Apply permutations IPIV to RHS
    swaps( 1, xb, m_nRC, 0, m_nRC-2, m_i_piv, 1 );

    // Solve for L part
    trsv(
      LOWER, NO_TRANSPOSE, UNIT,
      m_nRC, m_Afactorized, m_nRC, xb, 1
    );

    // Solve for U part
    trsv(
      UPPER, NO_TRANSPOSE, NON_UNIT,
      m_nRC, m_Afactorized, m_nRC, xb, 1
    );

    // Apply permutations JPIV to the solution (RHS)
    swaps( 1, xb, m_nRC, 0, m_nRC-2, m_j_piv, -1 );

    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LUPQ_no_alloc<T>::solve( integer nrhs, valueType B[], integer ldB ) const {
    // Apply permutations IPIV to RHS
    swaps( nrhs, B, ldB, 0, m_nRC-2, m_i_piv, 1 );

    // Solve for L part
    trsm(
      LEFT, LOWER, NO_TRANSPOSE, UNIT,
      m_nRC, nrhs, 1.0, m_Afactorized, m_nRC, B, ldB
    );

    // Solve for U part
    trsm(
      LEFT, UPPER, NO_TRANSPOSE, NON_UNIT,
      m_nRC, nrhs, 1.0, m_Afactorized, m_nRC, B, ldB
    );

    // Apply permutations JPIV to the solution (RHS)
    swaps( nrhs, B, ldB, 0, m_nRC-2, m_j_piv, -1 );

    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LUPQ_no_alloc<T>::t_solve( valueType xb[] ) const {
    // Apply permutations JPIV to the solution (RHS)
    swaps( 1, xb, m_nRC, 0, m_nRC-2, m_j_piv, 1 );

    // Solve for U part
    trsv(
      UPPER, TRANSPOSE, NON_UNIT,
      m_nRC, m_Afactorized, m_nRC, xb, 1
    );

    // Solve for L part
    trsv(
      LOWER, TRANSPOSE, UNIT,
      m_nRC, m_Afactorized, m_nRC, xb, 1
    );

    // Apply permutations IPIV to RHS
    swaps( 1, xb, m_nRC, 0, m_nRC-2, m_i_piv, -1 );

    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LUPQ_no_alloc<T>::t_solve( integer nrhs, valueType B[], integer ldB ) const {

    // Apply permutations JPIV to the solution (RHS)
    swaps( nrhs, B, ldB, 0, m_nRC-2, m_j_piv, 1 );

    // Solve for U part
    trsm(
      LEFT, UPPER, TRANSPOSE, NON_UNIT,
      m_nRC, nrhs, 1.0, m_Afactorized, m_nRC, B, ldB
    );

    // Solve for L part
    trsm(
      LEFT, LOWER, TRANSPOSE, UNIT,
      m_nRC, nrhs, 1.0, m_Afactorized, m_nRC, B, ldB
    );

    // Apply permutations IPIV to RHS
    swaps( nrhs, B, ldB, 0, m_nRC-2, m_i_piv, -1 );

    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  LUPQ<T>::LUPQ()
  : LUPQ_no_alloc<T>()
  , m_allocReals("LUPQ-allocReals")
  , m_allocIntegers("LUPQ-allocIntegers")
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  LUPQ<T>::~LUPQ() {
    m_allocReals.free();
    m_allocIntegers.free();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ<T>::allocate( integer NRC ) {
    if ( m_nRC != NRC ) {
      integer NRC2 = NRC*NRC;
      m_allocReals.allocate( size_t(NRC2) );
      m_allocIntegers.allocate( size_t(2*NRC) );
      this->no_allocate(
        NRC,
        NRC2,  m_allocReals( size_t(NRC2) ),
        2*NRC, m_allocIntegers( size_t(2*NRC) )
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

}

///
/// eof: lu.cxx
///
