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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU_no_alloc<T>::factorize_nodim(
    string_view     who,
    real_type const A[],
    integer         LDA
  ) {
    integer info{ gecopy( m_nrows, m_ncols, A, LDA, m_Afactorized, m_nrows ) };
    UTILS_ASSERT(
      info == 0,
      "LU::factorize[{}] gecopy(nRow={}, nCol={}, A, LDA={}, B, LDB={}) : INFO = {}\n",
      who, m_nrows, m_ncols, LDA, m_nrows, info
    );
    info = getrf(
      m_nrows, m_ncols, m_Afactorized, m_nrows, m_i_pivot
    );
    UTILS_ASSERT( info == 0, "LU::factorize[{}] getrf INFO = {}\n", who, info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LU_no_alloc<T>::factorize_nodim( real_type const A[], integer LDA ) {
    integer info{ gecopy( m_nrows, m_ncols, A, LDA, m_Afactorized, m_nrows ) };
    bool ok{ info == 0 };
    if ( ok ) {
      info = getrf( m_nrows, m_ncols, m_Afactorized, m_nrows, m_i_pivot );
      ok = info == 0;
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU_no_alloc<T>::check_ls( string_view who ) const {
    UTILS_ASSERT(
      m_nrows == m_ncols,
      "LU<T>::{}, rectangular matrix {} x {}\n", who, m_nrows, m_ncols
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LU_no_alloc<T>::solve( real_type xb[] ) const {
    if ( m_nrows != m_ncols ) return false;
    integer info{ getrs(
      Transposition::NO,
      m_nrows, 1, m_Afactorized, m_nrows, m_i_pivot,
      xb, m_nrows
    ) };
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU_no_alloc<T>::solve( string_view who, real_type xb[] ) const {
    check_ls("solve");
    integer info{ getrs(
      Transposition::NO,
      m_nrows, 1, m_Afactorized, m_nrows, m_i_pivot,
      xb, m_nrows
    ) };
    UTILS_ASSERT( info == 0, "LU::solve, getrs INFO = {}\nat {}\n", info, who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LU_no_alloc<T>::t_solve( real_type xb[] ) const {
    if ( m_nrows != m_ncols ) return false;
    integer info{ getrs(
      Transposition::YES,
      m_nrows, 1, m_Afactorized, m_nrows, m_i_pivot,
      xb, m_nrows
    ) };
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU_no_alloc<T>::t_solve( string_view who, real_type xb[] ) const {
    check_ls( who );
    integer info{ getrs(
      Transposition::YES,
      m_nrows, 1, m_Afactorized, m_nrows, m_i_pivot,
      xb, m_nrows
    ) };
    UTILS_ASSERT( info == 0, "LU::t_solve, getrs INFO = {}\nat {}\n", info, who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LU_no_alloc<T>::solve( integer nrhs, real_type B[], integer ldB ) const {
    if ( m_nrows != m_ncols ) return false;
    integer info{ getrs(
      Transposition::NO,
      m_nrows, nrhs, m_Afactorized, m_nrows, m_i_pivot,
      B, ldB
    ) };
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU_no_alloc<T>::solve(
    string_view who,
    integer     nrhs,
    real_type   B[],
    integer     ldB
  ) const {
    check_ls(who);
    integer info{ getrs(
      Transposition::NO,
      m_nrows, nrhs, m_Afactorized, m_nrows, m_i_pivot,
      B, ldB
    ) };
    UTILS_ASSERT( info == 0, "LU::solve getrs INFO = {}\nat {}\n", info, who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LU_no_alloc<T>::t_solve(
    integer    nrhs,
    real_type  B[],
    integer    ldB
  ) const {
    if ( m_nrows != m_ncols ) return false;
    integer info{ getrs(
      Transposition::YES,
      m_nrows, nrhs, m_Afactorized, m_nrows, m_i_pivot,
      B, ldB
    ) };
    return info >= 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU_no_alloc<T>::t_solve(
    string_view who,
    integer     nrhs,
    real_type   B[],
    integer     ldB
  ) const {
    check_ls( who );
    integer info{ getrs(
      Transposition::YES,
      m_nrows, nrhs, m_Afactorized, m_nrows, m_i_pivot,
      B, ldB
    ) };
    UTILS_ASSERT( info >= 0, "LU::t_solve getrs INFO = {}\nat {}\n", info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  typename LU_no_alloc<T>::real_type
  LU_no_alloc<T>::cond1( real_type norm1 ) const {
    real_type rcond;
    integer info{ gecon1(
      m_nrows, m_Afactorized, m_nrows,
      norm1, rcond, m_Work, m_Iwork
    ) };
    UTILS_ASSERT( info == 0, "LU::cond1, gecon1 return info = {}\n", info );
    return rcond;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  typename LU_no_alloc<T>::real_type
  LU_no_alloc<T>::condInf( real_type normInf ) const {
    real_type rcond;
    integer info { geconInf(
      m_nrows, m_Afactorized, m_nrows,
      normInf, rcond, m_Work, m_Iwork
    ) };
    UTILS_ASSERT( info == 0, "LU::condInf, geconInf return info = {}\n", info );
    return rcond;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU<T>::allocate( integer NR, integer NC ) {
    if ( m_nrows != NR || m_ncols != NC ) {
      m_nrows = NR;
      m_ncols = NC;
      integer NRC { NR*NC };
      integer NW  { 2*(NR+NC) };
      m_allocReals.reallocate( size_t(NRC+NW) );
      m_allocIntegers.reallocate( size_t(2*NR) );

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
  void
  LUPQ_no_alloc<T>::no_allocate(
    integer     NRC,
    integer     Lwork,
    real_type * Work,
    integer     Liwork,
    integer   * iWork
  ) {
    UTILS_ASSERT(
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
  LUPQ_no_alloc<T>::factorize_nodim(
    string_view     who,
    real_type const A[],
    integer         LDA
  ) {
    integer info{ gecopy( m_nRC, m_nRC, A, LDA, m_Afactorized, m_nRC ) };
    UTILS_ASSERT( info == 0, "LUPQ_no_alloc::factorize[{}] gecopy INFO = {}\n", who, info );
    info = getc2( m_nRC, m_Afactorized, m_nRC, m_i_piv, m_j_piv );
    UTILS_ASSERT( info == 0, "LUPQ_no_alloc::factorize[{}] getc2 INFO = {}\n", who, info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LUPQ_no_alloc<T>::factorize_nodim( real_type const A[], integer LDA ) {
    integer info{ gecopy( m_nRC, m_nRC, A, LDA, m_Afactorized, m_nRC ) };
    bool ok{ info == 0 };
    if ( ok ) {
      info = getc2( m_nRC, m_Afactorized, m_nRC, m_i_piv, m_j_piv );
      ok = info == 0;
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LUPQ_no_alloc<T>::solve( real_type xb[] ) const {
    // Apply permutations IPIV to RHS
    swaps( 1, xb, m_nRC, 0, m_nRC-2, m_i_piv, 1 );

    // Solve for L part
    trsv(
      ULselect::LOWER,
      Transposition::NO,
      DiagonalType::UNIT,
      m_nRC, m_Afactorized, m_nRC, xb, 1
    );

    // Solve for U part
    trsv(
      ULselect::UPPER,
      Transposition::NO,
      DiagonalType::NON_UNIT,
      m_nRC, m_Afactorized, m_nRC, xb, 1
    );

    // Apply permutations JPIV to the solution (RHS)
    swaps( 1, xb, m_nRC, 0, m_nRC-2, m_j_piv, -1 );

    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LUPQ_no_alloc<T>::solve( integer nrhs, real_type B[], integer ldB ) const {
    // Apply permutations IPIV to RHS
    swaps( nrhs, B, ldB, 0, m_nRC-2, m_i_piv, 1 );

    // Solve for L part
    trsm(
      SideMultiply::LEFT,
      ULselect::LOWER,
      Transposition::NO,
      DiagonalType::UNIT,
      m_nRC, nrhs, 1.0, m_Afactorized, m_nRC, B, ldB
    );

    // Solve for U part
    trsm(
      SideMultiply::LEFT,
      ULselect::UPPER,
      Transposition::NO,
      DiagonalType::NON_UNIT,
      m_nRC, nrhs, 1.0, m_Afactorized, m_nRC, B, ldB
    );

    // Apply permutations JPIV to the solution (RHS)
    swaps( nrhs, B, ldB, 0, m_nRC-2, m_j_piv, -1 );

    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LUPQ_no_alloc<T>::t_solve( real_type xb[] ) const {
    // Apply permutations JPIV to the solution (RHS)
    swaps( 1, xb, m_nRC, 0, m_nRC-2, m_j_piv, 1 );

    // Solve for U part
    trsv(
      ULselect::UPPER,
      Transposition::YES,
      DiagonalType::NON_UNIT,
      m_nRC, m_Afactorized, m_nRC, xb, 1
    );

    // Solve for L part
    trsv(
      ULselect::LOWER,
      Transposition::YES,
      DiagonalType::UNIT,
      m_nRC, m_Afactorized, m_nRC, xb, 1
    );

    // Apply permutations IPIV to RHS
    swaps( 1, xb, m_nRC, 0, m_nRC-2, m_i_piv, -1 );

    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LUPQ_no_alloc<T>::t_solve( integer nrhs, real_type B[], integer ldB ) const {

    // Apply permutations JPIV to the solution (RHS)
    swaps( nrhs, B, ldB, 0, m_nRC-2, m_j_piv, 1 );

    // Solve for U part
    trsm(
      SideMultiply::LEFT,
      ULselect::UPPER,
      Transposition::YES,
      DiagonalType::NON_UNIT,
      m_nRC, nrhs, 1.0, m_Afactorized, m_nRC, B, ldB
    );

    // Solve for L part
    trsm(
      SideMultiply::LEFT,
      ULselect::LOWER,
      Transposition::YES,
      DiagonalType::UNIT,
      m_nRC, nrhs, 1.0, m_Afactorized, m_nRC, B, ldB
    );

    // Apply permutations IPIV to RHS
    swaps( nrhs, B, ldB, 0, m_nRC-2, m_i_piv, -1 );

    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ<T>::allocate( integer NRC ) {
    if ( m_nRC != NRC ) {
      integer NRC2{ NRC*NRC };
      this->no_allocate(
        NRC,
        NRC2,  m_allocReals.realloc( size_t(NRC2) ),
        2*NRC, m_allocIntegers.realloc( size_t(2*NRC) )
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

}

///
/// eof: lu.cxx
///
