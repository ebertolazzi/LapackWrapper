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
/// file: band.cxx
///

namespace lapack_wrapper {

  /*\
   |   ____                  _          _ __  __       _        _
   |  | __ )  __ _ _ __   __| | ___  __| |  \/  | __ _| |_ _ __(_)_  __
   |  |  _ \ / _` | '_ \ / _` |/ _ \/ _` | |\/| |/ _` | __| '__| \ \/ /
   |  | |_) | (_| | | | | (_| |  __/ (_| | |  | | (_| | |_| |  | |>  <
   |  |____/ \__,_|_| |_|\__,_|\___|\__,_|_|  |_|\__,_|\__|_|  |_/_/\_\
  \*/

  template <typename T>
  BandedLU<T>::BandedLU()
  : m_allocReals("_BandedLU_reals")
  , m_allocIntegers("_BandedLU_integers")
  , m_m(0)
  , m_n(0)
  , m_nL(0)
  , m_nU(0)
  , m_ldAB(0)
  , m_is_factorized(false)
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  BandedLU<T>::~BandedLU()
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::setup(
    integer _m,
    integer _n,
    integer _nL,
    integer _nU
  ) {
    m_m    = _m;
    m_n    = _n;
    m_nL   = _nL;
    m_nU   = _nU;
    m_ldAB = 2*m_nL+m_nU+1;
    integer nnz = m_n*m_ldAB;
    m_allocReals.allocate( nnz );
    m_allocIntegers.allocate( m_m );
    m_AB            = m_allocReals( nnz );
    m_ipiv          = m_allocIntegers( m_m );
    m_is_factorized = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  BandedLU<T>::solve( valueType xb[] ) const {
    if ( !m_is_factorized ) return false;
    UTILS_ASSERT0( m_m == m_n, "BandedLU::solve, matrix must be square\n" );
    integer info = gbtrs(
      NO_TRANSPOSE, m_m, m_nL, m_nU, 1, m_AB, m_ldAB, m_ipiv, xb, m_m
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::solve( char const who[], valueType xb[] ) const {
    UTILS_ASSERT0( m_is_factorized, "BandedLU::solve, matrix not yet factorized\n" );
    UTILS_ASSERT0( m_m == m_n, "BandedLU::solve, matrix must be square\n" );
    integer info = gbtrs(
      NO_TRANSPOSE, m_m, m_nL, m_nU, 1, m_AB, m_ldAB, m_ipiv, xb, m_m
    );
    UTILS_ASSERT( info == 0, "BandedLU::solve, info = {}\nat {}\n", info, who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  BandedLU<T>::t_solve( valueType xb[] ) const {
    if ( !m_is_factorized ) return false;
    UTILS_ASSERT0( m_m == m_n, "BandedLU::solve, matrix must be square\n" );
    integer info = gbtrs(
      TRANSPOSE, m_m, m_nL, m_nU, 1, m_AB, m_ldAB, m_ipiv, xb, m_m
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::t_solve( char const who[], valueType xb[] ) const {
    UTILS_ASSERT0(
      m_is_factorized,
      "BandedLU::solve, matrix not yet factorized\n"
    );
    UTILS_ASSERT0(
      m_m == m_n,
      "BandedLU::solve, matrix must be square\n"
    );
    integer info = gbtrs(
      TRANSPOSE, m_m, m_nL, m_nU, 1, m_AB, m_ldAB, m_ipiv, xb, m_m
    );
    UTILS_ASSERT( info == 0, "BandedLU::t_solve, info = {}\nat {}\n", info, who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  BandedLU<T>::solve( integer nrhs, valueType B[], integer ldB ) const {
    if ( !m_is_factorized ) return false;
    UTILS_ASSERT0( m_m == m_n, "BandedLU::solve, matrix must be square\n" );
    integer info = gbtrs(
      NO_TRANSPOSE, m_m, m_nL, m_nU, nrhs, m_AB, m_ldAB, m_ipiv, B, ldB
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::solve(
    char const who[],
    integer    nrhs,
    valueType  B[],
    integer    ldB
  ) const {
    UTILS_ASSERT0( m_is_factorized, "BandedLU::solve, matrix not yet factorized\n" );
    UTILS_ASSERT0( m_m == m_n, "BandedLU::solve, matrix must be square\n" );
    integer info = gbtrs(
      NO_TRANSPOSE, m_m, m_nL, m_nU, nrhs, m_AB, m_ldAB, m_ipiv, B, ldB
    );
    UTILS_ASSERT( info == 0, "BandedLU::solve, info = {}\nat {}\n", info, who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  BandedLU<T>::t_solve( integer nrhs, valueType B[], integer ldB ) const {
    if ( !m_is_factorized ) return false;
    UTILS_ASSERT0( m_m == m_n, "BandedLU::solve, matrix must be square\n" );
    integer info = gbtrs(
      TRANSPOSE, m_m, m_nL, m_nU, nrhs, m_AB, m_ldAB, m_ipiv, B, ldB
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::t_solve(
    char const who[],
    integer    nrhs,
    valueType  B[],
    integer    ldB
  ) const {
    UTILS_ASSERT0( m_is_factorized, "BandedLU::solve, matrix not yet factorized\n" );
    UTILS_ASSERT0( m_m == m_n, "BandedLU::solve, matrix must be square\n" );
    integer info = gbtrs(
      TRANSPOSE, m_m, m_nL, m_nU, nrhs, m_AB, m_ldAB, m_ipiv, B, ldB
    );
    UTILS_ASSERT( info == 0, "BandedLU::t_solve, info = {}\nat {}\n", info, who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::factorize( char const who[] ) {
    UTILS_ASSERT(
      !m_is_factorized, "BandedLU::factorize[{}], matrix yet factorized\n", who
    );
    UTILS_ASSERT(
      m_m == m_n, "BandedLU::factorize[{}], matrix must be square\n", who
    );
    integer info = gbtrf( m_m, m_n, m_nL, m_nU, m_AB, m_ldAB, m_ipiv );
    UTILS_ASSERT(
      info == 0, "BandedLU::factorize[{}], info = {}\n", who, info
    );
    m_is_factorized = true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  BandedLU<T>::factorize() {
    if ( m_is_factorized ) return true;
    if ( m_m != m_n ) return false;
    integer info = gbtrf( m_m, m_n, m_nL, m_nU, m_AB, m_ldAB, m_ipiv );
    m_is_factorized = info == 0;
    return m_is_factorized;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::zero() {
    integer nnz = m_m*(2*m_nL+m_nU+1);
    lapack_wrapper::zero( nnz, m_AB, 1 );
    m_is_factorized = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::check( integer i, integer j ) const {
    UTILS_ASSERT(
      i >= 0 && i < m_m && j >= 0 && j < m_n,
      "BandedLU::check( {}, {} ) indices out of range: [0..{}) x [0..{})\n",
      i, j, m_m, m_n
    );
    UTILS_ASSERT(
      j >= i-m_nL && j <= i+m_nU,
      "BandedLU::check( {}, {} ) indices out of band: nLower: {} nUpper: {}\n",
      i, j, m_nL, m_nU
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::insert(
    integer   i,
    integer   j,
    valueType v,
    bool      sym
  ) {
    UTILS_ASSERT(
      i >= 0 && i < m_m && j >= 0 && j < m_n,
      "BandedLU::insert( {}, {} ) out of range: [0..{}) x [0..{})\n",
      i, j, m_m, m_n
    );
    UTILS_ASSERT(
      i-m_nL <= j && j <= i+m_nU,
      "BandedLU::insert( {}, {} ) out of band: nLower: {} nUpper: {}\n",
      i, j, m_nL, m_nU
    );
    (*this)(i,j) = v;
    if ( sym && i != j ) {
      UTILS_ASSERT(
        j-m_nL <= i && i <= j+m_nU,
        "BandedLU::insert( {}, {} ) out of band: nLower: {} nUpper: {}\n",
        i, j, m_nL, m_nU
      );
      (*this)(j,i) = v;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::load_block(
    integer         nr,
    integer         nc,
    valueType const B[],
    integer         ldB,
    integer         irow,
    integer         icol
  ) {

    UTILS_ASSERT0(
      !m_is_factorized, "BandedLU::load_block, matrix is factorized\n"
    );

    #if 1
    for ( integer r = 0; r < nr; ++r )
      for ( integer c = 0; c < nc; ++c )
        m_AB[iaddr( irow+r, icol+c )] = B[ r + c * ldB ];
    #else
    // must be checked
    check( irow,      icol      );
    check( irow+nr-1, icol+nc-1 );
    check( irow,      icol+nc-1 );
    check( irow+nr-1, icol      );

    // copy by diagonal
    for ( integer r = 0; r < nr; ++r ) {
      integer ia = iaddr( irow+r, icol );
      copy( std::min(nr-r,nc), B+r, ldB+1, AB+ia, ldAB );
    }
    for ( integer c = 1; c < nc; ++c ) {
      integer ia = iaddr( irow, icol+c );
      copy( std::min(nc-c,nr), B+c*ldB, ldB+1, AB+ia, ldAB );
    }
    #endif
  }

  // y <- beta*y + alpha*A*x
  /*
    +---------+
    | \       |
    |  \      |
    +---+-----+
  */
  template <typename T>
  void
  BandedLU<T>::aAxpy(
    valueType       alpha,
    valueType const x[],
    valueType       y[]
  ) const {
    valueType const * col = m_AB + m_nL;
    for ( integer j = 0; j < m_n; ++j, col += m_ldAB ) {
      integer imin  = j-m_nU;
      integer imax  = std::min(j+m_nL,m_m-1);
      integer imin0 = imin > 0 ? imin : 0;
      lapack_wrapper::axpy(
        imax-imin0+1,
        alpha*x[j],
        col+imin0-imin, 1,
        y+imin0,        1
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::dump( ostream_type & stream ) const {
    for ( integer i = 0; i <= m_nL+m_nU; ++i ) {
      valueType const * col = m_AB + m_nL + i;
      for ( integer j = 0; j < m_n; ++j, col += m_ldAB )
        stream << std::setw(10) << col[0] << ' ';
      stream << '\n';
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  BandedSPD<T>::BandedSPD()
  : m_allocReals("_BandedSPD_reals")
  , m_n(0)
  , m_nD(0)
  , m_ldAB(0)
  , m_is_factorized(false)
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  BandedSPD<T>::~BandedSPD()
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::setup(
    ULselect _UPLO,
    integer  _N,
    integer  _nD
  ) {
    m_UPLO = _UPLO;
    m_n    = _N;
    m_nD   = _nD;
    m_ldAB = _nD+1;
    integer nnz = m_n*m_ldAB;
    m_allocReals.allocate( nnz );
    m_AB   = m_allocReals( nnz );
    m_is_factorized = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  BandedSPD<T>::solve( valueType xb[] ) const {
    if ( !m_is_factorized ) return false;
    integer info = pbtrs( m_UPLO, m_n, m_nD, 1, m_AB, m_ldAB, xb, m_n );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::solve( char const who[], valueType xb[] ) const {
    UTILS_ASSERT0( m_is_factorized, "BandedSPD::solve, matrix not yet factorized\n" );
    integer info = pbtrs( m_UPLO, m_n, m_nD, 1, m_AB, m_ldAB, xb, m_n );
    UTILS_ASSERT( info == 0, "BandedSPD::solve, info = {}\nat {}\n", info, who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  BandedSPD<T>::t_solve( valueType xb[] ) const {
    if ( !m_is_factorized ) return false;
    integer info = pbtrs( m_UPLO, m_n, m_nD, 1, m_AB, m_ldAB, xb, m_n );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::t_solve( char const who[], valueType xb[] ) const {
    UTILS_ASSERT0( m_is_factorized, "BandedSPD::solve, matrix not yet factorized\n" );
    integer info = pbtrs( m_UPLO, m_n, m_nD, 1, m_AB, m_ldAB, xb, m_n );
    UTILS_ASSERT( info == 0, "BandedSPD::t_solve, info = {}\nat {}\n", info, who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  BandedSPD<T>::solve( integer nrhs, valueType B[], integer ldB ) const {
    if ( !m_is_factorized ) return false;
    integer info = pbtrs( m_UPLO, m_n, m_nD, nrhs, m_AB, m_ldAB, B, ldB );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::solve( char const who[], integer nrhs, valueType B[], integer ldB ) const {
    UTILS_ASSERT0( m_is_factorized, "BandedSPD::solve, matrix not yet factorized\n" );
    integer info = pbtrs( m_UPLO, m_n, m_nD, nrhs, m_AB, m_ldAB, B, ldB );
    UTILS_ASSERT( info == 0, "BandedSPD::solve, info = {}\nat {}\n", info, who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  BandedSPD<T>::t_solve( integer nrhs, valueType B[], integer ldB ) const {
    if ( !m_is_factorized ) return false;
    integer info = pbtrs( m_UPLO, m_n, m_nD, nrhs, m_AB, m_ldAB, B, ldB );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::t_solve( char const who[], integer nrhs, valueType B[], integer ldB ) const {
    UTILS_ASSERT0( m_is_factorized, "BandedSPD::solve, matrix not yet factorized\n" );
    integer info = pbtrs( m_UPLO, m_n, m_nD, nrhs, m_AB, m_ldAB, B, ldB );
    UTILS_ASSERT( info == 0, "BandedSPD::t_solve, info = {}\nat {}\n", info, who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::factorize( char const who[] ) {
    UTILS_ASSERT(
      !m_is_factorized, "BandedSPD::factorize[{}], matrix yet factorized\n", who
    );
    integer info = pbtrf( m_UPLO, m_n, m_nD, m_AB, m_ldAB );
    UTILS_ASSERT( info == 0, "BandedSPD::factorize[{}], info = {}\n", who, info );
    m_is_factorized = true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  BandedSPD<T>::factorize() {
    if ( m_is_factorized ) return true;
    integer info = pbtrf( m_UPLO, m_n, m_nD, m_AB, m_ldAB );
    m_is_factorized = info == 0;
    return m_is_factorized;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::zero() {
    lapack_wrapper::zero( m_n*m_ldAB, m_AB, 1 );
    m_is_factorized = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::insert(
    integer   i,
    integer   j,
    valueType v,
    bool      sym
  ) {
    (*this)(i,j) = v;
    if ( sym && i != j ) (*this)(j,i) = v;
  }

}

///
/// eof: band.cxx
///
