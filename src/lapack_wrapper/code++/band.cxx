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

  template <typename T>
  void
  BandedLU<T>::solve( valueType xb[] ) const {
    LAPACK_WRAPPER_ASSERT( is_factorized, "BandedLU::solve, matrix not yet factorized" );
    LAPACK_WRAPPER_ASSERT( m == n, "BandedLU::solve, matrix must be square" );
    integer info = gbtrs( NO_TRANSPOSE, m, nL, nU, 1, AB, ldAB, ipiv, xb, m );
    LAPACK_WRAPPER_ASSERT( info == 0, "BandedLU::solve, info = " << info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::t_solve( valueType xb[] ) const {
    LAPACK_WRAPPER_ASSERT( is_factorized, "BandedLU::solve, matrix not yet factorized" );
    LAPACK_WRAPPER_ASSERT( m == n, "BandedLU::solve, matrix must be square" );
    integer info = gbtrs( TRANSPOSE, m, nL, nU, 1, AB, ldAB, ipiv, xb, m );
    LAPACK_WRAPPER_ASSERT( info == 0, "BandedLU::t_solve, info = " << info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::solve( integer nrhs, valueType B[], integer ldB ) const {
    LAPACK_WRAPPER_ASSERT( is_factorized, "BandedLU::solve, matrix not yet factorized" );
    LAPACK_WRAPPER_ASSERT( m == n, "BandedLU::solve, matrix must be square" );
    integer info = gbtrs( NO_TRANSPOSE, m, nL, nU, nrhs, AB, ldAB, ipiv, B, ldB );
    LAPACK_WRAPPER_ASSERT( info == 0, "BandedLU::solve, info = " << info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::t_solve( integer nrhs, valueType B[], integer ldB ) const {
    LAPACK_WRAPPER_ASSERT( is_factorized, "BandedLU::solve, matrix not yet factorized" );
    LAPACK_WRAPPER_ASSERT( m == n, "BandedLU::solve, matrix must be square" );
    integer info = gbtrs( TRANSPOSE, m, nL, nU, nrhs, AB, ldAB, ipiv, B, ldB );
    LAPACK_WRAPPER_ASSERT( info == 0, "BandedLU::t_solve, info = " << info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::factorize( char const      who[] ) {
    LAPACK_WRAPPER_ASSERT(
      !is_factorized,
      "BandedLU::factorize[" << who << "], matrix yet factorized"
    );
    LAPACK_WRAPPER_ASSERT(
      m == n,
      "BandedLU::factorize[" << who << "], matrix must be square"
    );
    integer info = gbtrf( m, n, nL, nU, AB, ldAB, ipiv );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "BandedLU::factorize[" << who << "], info = " << info
    );
    is_factorized = true;
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
    LAPACK_WRAPPER_ASSERT(
      i >= 0 && i < m && j >= 0 && j < n,
      "BandedLU::check( " << i << " , " << j << " ) out of range"
    );
    LAPACK_WRAPPER_ASSERT(
      j >= i-nL && j <= i+nU,
      "BandedLU::check( " << i << " , " << j << " ) out of band"
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
    LAPACK_WRAPPER_ASSERT(
      i >= 0 && i < m && j >= 0 && j < n,
      "BandedLU::insert( " << i << " , " << j << " ) out of range"
    );
    LAPACK_WRAPPER_ASSERT(
      i-nL <= j && j <= i+nU,
      "BandedLU::insert( " << i << " , " << j << " ) out of band"
    );
    (*this)(i,j) = v;
    if ( sym && i != j ) {
      LAPACK_WRAPPER_ASSERT(
        j-nL <= i && i <= j+nU,
        "BandedLU::insert( " << i << " , " << j << " ) out of band"
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

    LAPACK_WRAPPER_ASSERT(
      !is_factorized,
      "BandedLU::load_block, matrix is factorized"
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
  : allocReals("_BandedSPD_reals")
  , n(0)
  , nD(0)
  , ldAB(0)
  , is_factorized(false)
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
    UPLO = _UPLO;
    n    = _N;
    nD   = _nD;
    ldAB = nD+1;
    integer nnz = n*ldAB;
    allocReals.allocate( nnz );
    AB   = allocReals( nnz );
    is_factorized = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::solve( valueType xb[] ) const {
    LAPACK_WRAPPER_ASSERT(
      is_factorized,
      "BandedSPD::solve, matrix not yet factorized"
    );
    integer info = pbtrs( UPLO, n, nD, 1, AB, ldAB, xb, n );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "BandedSPD::solve, info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::t_solve( valueType xb[] ) const {
    LAPACK_WRAPPER_ASSERT(
      is_factorized,
      "BandedSPD::solve, matrix not yet factorized"
    );
    integer info = pbtrs( UPLO, n, nD, 1, AB, ldAB, xb, n );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "BandedSPD::t_solve, info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::solve( integer nrhs, valueType B[], integer ldB ) const {
    LAPACK_WRAPPER_ASSERT(
      is_factorized,
      "BandedSPD::solve, matrix not yet factorized"
    );
    integer info = pbtrs( UPLO, n, nD, nrhs, AB, ldAB, B, ldB );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "BandedSPD::solve, info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::t_solve( integer nrhs, valueType B[], integer ldB ) const {
    LAPACK_WRAPPER_ASSERT(
      is_factorized,
      "BandedSPD::solve, matrix not yet factorized"
    );
    integer info = pbtrs( UPLO, n, nD, nrhs, AB, ldAB, B, ldB );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "BandedSPD::t_solve, info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::factorize( char const who[] ) {
    LAPACK_WRAPPER_ASSERT(
      !is_factorized,
      "BandedSPD::factorize[" << who << "], matrix yet factorized"
    );
    integer info = pbtrf( UPLO, n, nD, AB, ldAB );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "BandedSPD::factorize[" << who << "], info = " << info
    );
    is_factorized = true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::zero() {
    lapack_wrapper::zero( n*ldAB, AB, 1 );
    is_factorized = false;
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
