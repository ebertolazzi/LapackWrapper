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
  : allocReals("_BandedLU_reals")
  , allocIntegers("_BandedLU_integers")
  , m(0)
  , n(0)
  , nL(0)
  , nU(0)
  , ldAB(0)
  , is_factorized(false)
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
    m    = _m;
    n    = _n;
    nL   = _nL;
    nU   = _nU;
    ldAB = 2*nL+nU+1;
    integer nnz = n*ldAB;
    allocReals.allocate( nnz );
    allocIntegers.allocate(m);
    AB   = allocReals( nnz );
    ipiv = allocIntegers( m );
    is_factorized = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  BandedLU<T>::solve( valueType xb[] ) const {
    if ( !is_factorized ) return false;
    LW_ASSERT0( m == n, "BandedLU::solve, matrix must be square\n" );
    integer info = gbtrs( NO_TRANSPOSE, m, nL, nU, 1, AB, ldAB, ipiv, xb, m );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::solve( char const who[], valueType xb[] ) const {
    LW_ASSERT0( is_factorized, "BandedLU::solve, matrix not yet factorized\n" );
    LW_ASSERT0( m == n, "BandedLU::solve, matrix must be square\n" );
    integer info = gbtrs( NO_TRANSPOSE, m, nL, nU, 1, AB, ldAB, ipiv, xb, m );
    LW_ASSERT( info == 0, "BandedLU::solve, info = {}\nat {}\n", info, who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  BandedLU<T>::t_solve( valueType xb[] ) const {
    if ( !is_factorized ) return false;
    LW_ASSERT0( m == n, "BandedLU::solve, matrix must be square\n" );
    integer info = gbtrs( TRANSPOSE, m, nL, nU, 1, AB, ldAB, ipiv, xb, m );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::t_solve( char const who[], valueType xb[] ) const {
    LW_ASSERT0( is_factorized, "BandedLU::solve, matrix not yet factorized\n" );
    LW_ASSERT0( m == n, "BandedLU::solve, matrix must be square\n" );
    integer info = gbtrs( TRANSPOSE, m, nL, nU, 1, AB, ldAB, ipiv, xb, m );
    LW_ASSERT( info == 0, "BandedLU::t_solve, info = {}\nat {}\n", info, who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  BandedLU<T>::solve( integer nrhs, valueType B[], integer ldB ) const {
    if ( !is_factorized ) return false;
    LW_ASSERT0( m == n, "BandedLU::solve, matrix must be square\n" );
    integer info = gbtrs( NO_TRANSPOSE, m, nL, nU, nrhs, AB, ldAB, ipiv, B, ldB );
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
    LW_ASSERT0( is_factorized, "BandedLU::solve, matrix not yet factorized\n" );
    LW_ASSERT0( m == n, "BandedLU::solve, matrix must be square\n" );
    integer info = gbtrs( NO_TRANSPOSE, m, nL, nU, nrhs, AB, ldAB, ipiv, B, ldB );
    LW_ASSERT( info == 0, "BandedLU::solve, info = {}\nat {}\n", info, who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  BandedLU<T>::t_solve( integer nrhs, valueType B[], integer ldB ) const {
    if ( !is_factorized ) return false;
    LW_ASSERT0( m == n, "BandedLU::solve, matrix must be square\n" );
    integer info = gbtrs( TRANSPOSE, m, nL, nU, nrhs, AB, ldAB, ipiv, B, ldB );
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
    LW_ASSERT0( is_factorized, "BandedLU::solve, matrix not yet factorized\n" );
    LW_ASSERT0( m == n, "BandedLU::solve, matrix must be square\n" );
    integer info = gbtrs( TRANSPOSE, m, nL, nU, nrhs, AB, ldAB, ipiv, B, ldB );
    LW_ASSERT( info == 0, "BandedLU::t_solve, info = {}\nat {}\n", info, who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::factorize( char const who[] ) {
    LW_ASSERT(
      !is_factorized, "BandedLU::factorize[{}], matrix yet factorized\n", who
    );
    LW_ASSERT(
      m == n, "BandedLU::factorize[{}], matrix must be square\n", who
    );
    integer info = gbtrf( m, n, nL, nU, AB, ldAB, ipiv );
    LW_ASSERT(
      info == 0, "BandedLU::factorize[{}], info = {}\n", who, info
    );
    is_factorized = true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  BandedLU<T>::factorize() {
    if ( is_factorized ) return true;
    if ( m != n ) return false;
    integer info = gbtrf( m, n, nL, nU, AB, ldAB, ipiv );
    is_factorized = info == 0;
    return is_factorized;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::zero() {
    integer nnz = m*(2*nL+nU+1);
    lapack_wrapper::zero( nnz, AB, 1 );
    is_factorized = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::check( integer i, integer j ) const {
    LW_ASSERT(
      i >= 0 && i < m && j >= 0 && j < n,
      "BandedLU::check( {}, {} ) indices out of range: [0..{}) x [0..{})\n",
      i, j, m, n
    );
    LW_ASSERT(
      j >= i-nL && j <= i+nU,
      "BandedLU::check( {}, {} ) indices out of band: nLower: {} nUpper: {}\n",
      i, j, nL, nU
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
    LW_ASSERT(
      i >= 0 && i < m && j >= 0 && j < n,
      "BandedLU::insert( {}, {} ) out of range: [0..{}) x [0..{})\n",
      i, j, m, n
    );
    LW_ASSERT(
      i-nL <= j && j <= i+nU,
      "BandedLU::insert( {}, {} ) out of band: nLower: {} nUpper: {}\n",
      i, j, nL, nU
    );
    (*this)(i,j) = v;
    if ( sym && i != j ) {
      LW_ASSERT(
        j-nL <= i && i <= j+nU,
        "BandedLU::insert( {}, {} ) out of band: nLower: {} nUpper: {}\n",
        i, j, nL, nU
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

    LW_ASSERT0(
      !is_factorized, "BandedLU::load_block, matrix is factorized\n"
    );

    #if 1
    for ( integer r = 0; r < nr; ++r )
      for ( integer c = 0; c < nc; ++c )
        AB[iaddr( irow+r, icol+c )] = B[ r + c * ldB ];
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

    valueType const * col = AB + nL;
    for ( integer j = 0; j < n; ++j, col += ldAB ) {
      integer imin  = j-nU;
      integer imax  = std::min(j+nL,m-1);
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
    for ( integer i = 0; i <= nL+nU; ++i ) {
      valueType const * col = AB + nL + i;
      for ( integer j = 0; j < n; ++j, col += ldAB )
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
    LW_ASSERT0( m_is_factorized, "BandedSPD::solve, matrix not yet factorized\n" );
    integer info = pbtrs( m_UPLO, m_n, m_nD, 1, m_AB, m_ldAB, xb, m_n );
    LW_ASSERT( info == 0, "BandedSPD::solve, info = {}\nat {}\n", info, who );
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
    LW_ASSERT0( m_is_factorized, "BandedSPD::solve, matrix not yet factorized\n" );
    integer info = pbtrs( m_UPLO, m_n, m_nD, 1, m_AB, m_ldAB, xb, m_n );
    LW_ASSERT( info == 0, "BandedSPD::t_solve, info = {}\nat {}\n", info, who );
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
    LW_ASSERT0( m_is_factorized, "BandedSPD::solve, matrix not yet factorized\n" );
    integer info = pbtrs( m_UPLO, m_n, m_nD, nrhs, m_AB, m_ldAB, B, ldB );
    LW_ASSERT( info == 0, "BandedSPD::solve, info = {}\nat {}\n", info, who );
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
    LW_ASSERT0( m_is_factorized, "BandedSPD::solve, matrix not yet factorized\n" );
    integer info = pbtrs( m_UPLO, m_n, m_nD, nrhs, m_AB, m_ldAB, B, ldB );
    LW_ASSERT( info == 0, "BandedSPD::t_solve, info = {}\nat {}\n", info, who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::factorize( char const who[] ) {
    LW_ASSERT(
      !m_is_factorized, "BandedSPD::factorize[{}], matrix yet factorized\n", who
    );
    integer info = pbtrf( m_UPLO, m_n, m_nD, m_AB, m_ldAB );
    LW_ASSERT( info == 0, "BandedSPD::factorize[{}], info = {}\n", who, info );
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
