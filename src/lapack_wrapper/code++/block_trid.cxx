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
/// file: block_trid.cxx
///

namespace lapack_wrapper {

  //============================================================================
  /*\
   |  ___ _         _     _____    _    _ _                         _
   | | _ ) |___  __| |__ |_   _| _(_)__| (_)__ _ __ _ ___ _ _  __ _| |
   | | _ \ / _ \/ _| / /   | || '_| / _` | / _` / _` / _ \ ' \/ _` | |
   | |___/_\___/\__|_\_\   |_||_| |_\__,_|_\__,_\__, \___/_||_\__,_|_|
   |                                            |___/
   |  ___                _       _
   | / __|_  _ _ __  ___| |_ _ _(_)__
   | \__ \ || | '  \/ -_)  _| '_| / _|
   | |___/\_, |_|_|_\___|\__|_| |_\__|
   |      |__/
  \*/

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::setup(
    integer       nblks,
    integer const rBlocks[]
  ) {
    integer N = rBlocks[nblks], nrmax = 0;
    m_allocIntegers.reallocate( N + nblks+1 );
    m_allocRpointers.reallocate( 2*nblks-1 );
    m_allocIpointers.reallocate( nblks );
    m_row_blocks    = m_allocIntegers( nblks+1 );
    m_D_blocks      = m_allocRpointers( nblks );
    m_L_blocks      = m_allocRpointers( nblks-1 );
    m_B_permutation = m_allocIpointers( nblks );
    // evalute the memry usage for the L and D blocks
    integer nr0 = rBlocks[1] - rBlocks[0];
    nrmax = nr0;
    m_nnz = nr0*nr0;

    UTILS_ASSERT(
      nr0 >= 0, "BlockTridiagonalSymmetic::setup( nblks = {}, ... )\n", nr0
    );

    for ( integer i = 1; i < nblks; ++i ) {
      integer nr = rBlocks[i+1] - rBlocks[i];
      UTILS_ASSERT(
        nr >= 0, "BlockTridiagonalSymmetic::setup( nblks = {}, ... )\n", nr
      );
      m_nnz += nr*(nr+nr0);
      if ( nr > nrmax ) nrmax = nr;
      nr0 = nr;
    }
    m_allocReals.reallocate( m_nnz + nrmax * nrmax );
    nr0 = rBlocks[1] - rBlocks[0];
    m_D_blocks[0] = m_allocReals( nr0 * nr0 );
    m_B_permutation[0] = m_allocIntegers( nr0 );
    for ( integer i = 1; i < nblks; ++i ) {
      integer nr = rBlocks[i+1] - rBlocks[i];
      m_D_blocks[i]   = m_allocReals( nr * nr );
      m_L_blocks[i-1] = m_allocReals( nr * nr0 );
      nr0 = nr;
      m_B_permutation[i] = m_allocIntegers( nr );
    }

    m_Work = m_allocReals( nrmax * nrmax );

    this->zero();
    m_nBlocks = nblks;
    std::copy_n( rBlocks, nblks+1, m_row_blocks );
    m_is_factorized = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::setup(
    integer nblks, integer const block_size
  ) {
    // use a temporary vector
    std::vector<integer> rBlocks;
    rBlocks.reserve(nblks+1);
    integer n{0};
    for ( integer k = 0; k <= nblks; ++k, n += block_size )
      rBlocks.push_back(n);
    this->setup( nblks, &rBlocks.front() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::zero() { // fill to 0 all the blocks
    lapack_wrapper::zero( m_nnz, m_D_blocks[0], 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  T const &
  BlockTridiagonalSymmetic<T>::operator () ( integer ii, integer jj ) const {
    integer const * ib = m_row_blocks;
    integer const * ie = m_row_blocks+m_nBlocks+1;
    integer iBlock = integer(std::upper_bound( ib, ie, ii ) - ib) - 1;
    integer jBlock = integer(std::upper_bound( ib, ie, jj ) - ib) - 1;

    UTILS_ASSERT(
      m_row_blocks[iBlock] <= ii && ii < m_row_blocks[iBlock+1],
      "BlockTridiagonalSymmetic( i = {}, j = {} ) bad iBlock\n", ii, jj
    );
    UTILS_ASSERT(
      m_row_blocks[jBlock] <= jj && jj < m_row_blocks[jBlock+1],
      "BlockTridiagonalSymmetic( i = {}, j = {} ) bad jBlock\n", ii, jj
    );

    integer nr = m_row_blocks[iBlock+1] - m_row_blocks[iBlock];
    integer i  = ii - m_row_blocks[iBlock];
    integer j  = jj - m_row_blocks[jBlock];
    UTILS_ASSERT(
      jBlock == iBlock || jBlock+1 == iBlock,
      "BlockTridiagonalSymmetic:find( {} , {} ) --> "
      "( iBlock = {}, jBlock = {} ) --> ( i = {}, j = {} ) out of range\n",
      ii, jj, iBlock, jBlock, i, j
    );

    integer ij = i+j*nr;
    if ( iBlock == jBlock ) {
      return m_D_blocks[jBlock][ij];
    } else {
      return m_L_blocks[jBlock][ij];
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  T &
  BlockTridiagonalSymmetic<T>::operator () ( integer ii, integer jj ) {
    integer const * ib = m_row_blocks;
    integer const * ie = m_row_blocks+m_nBlocks+1;
    integer iBlock = integer(std::upper_bound( ib, ie, ii ) - ib) - 1;
    integer jBlock = integer(std::upper_bound( ib, ie, jj ) - ib) - 1;

    UTILS_ASSERT(
      m_row_blocks[iBlock] <= ii && ii < m_row_blocks[iBlock+1],
      "BlockTridiagonalSymmetic( i = {}, j = {} ) bad iBlock\n", ii, jj
    );
    UTILS_ASSERT(
      m_row_blocks[jBlock] <= jj && jj < m_row_blocks[jBlock+1],
      "BlockTridiagonalSymmetic( i = {}, j = {} ) bad jBlock\n", ii, jj
    );

    integer nr = m_row_blocks[iBlock+1] - m_row_blocks[iBlock];
    integer i  = ii - m_row_blocks[iBlock];
    integer j  = jj - m_row_blocks[jBlock];
    UTILS_ASSERT(
      jBlock == iBlock || jBlock+1 == iBlock,
      "BlockTridiagonalSymmetic:operator () ( {}, {} ) --> "
      "( iBlock = {}, jBlock = {} ) --> ( i = {}, j = {} ) out of range\n",
      ii, jj, iBlock, jBlock, i, j
    );

    integer ij = i+j*nr;
    if ( iBlock == jBlock ) {
      return m_D_blocks[jBlock][ij];
    } else {
      return m_L_blocks[jBlock][ij];
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::insert(
    integer   ii,
    integer   jj,
    real_type v,
    bool      sym
  ) {
    integer const * ib = m_row_blocks;
    integer const * ie = m_row_blocks+m_nBlocks+1;
    integer jBlock = integer(std::upper_bound( ib, ie, jj ) - ib) - 1;
    integer jmin   = m_row_blocks[jBlock];
    integer jmax   = m_row_blocks[jBlock+1];
    integer j      = jj - jmin;
    UTILS_ASSERT(
      jmin <= jj && jj < jmax && jmin <= ii,
      "BlockTridiagonalSymmetic::insert( i = {}, j = {}, v = {}, sym = {} )\n",
      ii, jj, v, sym
    );
    if ( ii < jmax ) {
      integer i  = ii - jmin;
      integer nr = jmax-jmin;
      real_type * D = m_D_blocks[jBlock];
      D[i+j*nr] = v;
      if ( sym ) D[j+i*nr] = v;
    } else {
      integer imin = jmax;
      integer imax = m_row_blocks[jBlock+2];
      integer nr   = imax-imin;
      integer i    = ii - imin;
      m_L_blocks[jBlock][i+j*nr] = v;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::setD(
    integer         n,
    real_type const data[],
    integer         ldData,
    bool            transposed
  ) {
    integer nr = this->D_nrows(n);
    if ( transposed ) {
      getranspose( nr, nr, data, ldData, m_D_blocks[n], nr );
    } else {
      integer ierr = gecopy( nr, nr, data, ldData, m_D_blocks[n], nr );
      UTILS_ASSERT(
        ierr == 0,
        "BlockTridiagonalSymmetic::setD, gecopy return ierr = {}", ierr
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::setL(
    integer         n,
    real_type const data[],
    integer         ldData,
    bool            transposed
  ) {
    integer nr = this->L_nrows(n);
    integer nc = this->L_ncols(n);
    if ( transposed ) {
      getranspose( nr, nc, data, ldData, m_L_blocks[n], nr );
    } else {
      integer ierr = gecopy( nr, nc, data, ldData, m_L_blocks[n], nr );
      UTILS_ASSERT(
        ierr == 0,
        "BlockTridiagonalSymmetic::setL, gecopy return ierr = {}\n", ierr
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::setD(
    integer           n,
    real_type const * data,
    integer           ldData,
    integer           beginRow,
    integer           beginCol,
    integer           nrow,
    integer           ncol,
    bool              transposed
  ) {
    integer ldD = this->D_nrows(n);
    real_type * D = m_D_blocks[n]+beginRow+beginCol*ldD;

    if ( transposed ) {
      getranspose( nrow, ncol, data, ldData, D, ldD );
    } else {
      integer ierr = gecopy( nrow, ncol, data, ldData, D, ldD );
      UTILS_ASSERT(
        ierr == 0,
        "BlockTridiagonalSymmetic::setD (block), gecopy return ierr = {}\n", ierr
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::setL(
    integer           n,
    real_type const * data,
    integer           ldData,
    integer           beginRow,
    integer           beginCol,
    integer           nrow,
    integer           ncol,
    bool              transposed
  ) {
    integer ldL = this->L_nrows(n);
    real_type * L = m_L_blocks[n]+beginRow+beginCol*ldL;
    if ( transposed ) {
      getranspose( nrow, ncol, data, ldData, L, ldL );
    } else {
      integer ierr = gecopy( nrow, ncol, data, ldData, L, ldL );
      UTILS_ASSERT(
        ierr == 0,
        "BlockTridiagonalSymmetic::setL (block), gecopy return ierr = {}\n", ierr
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::factorize( char const who[] ) {
    UTILS_ASSERT(
      !m_is_factorized,
      "BlockTridiagonalSymmetic::factorize[{}], already factored\n", who
    );

    integer info = lapack_wrapper::getrf(
      this->D_nrows(0), this->D_ncols(0),
      m_D_blocks[0], this->D_nrows(0),
      m_B_permutation[0]
    );

    UTILS_ASSERT(
      info == 0,
      "BlockTridiagonalSymmetic::factorize[{}] getrf INFO = {}\n", who, info
    );

    for ( integer k=1; k < m_nBlocks; ++k ) {
      integer nr0 = this->D_nrows(k-1);
      integer nr1 = this->D_nrows(k);
      real_type * L0 = m_L_blocks[k-1];
      real_type * D0 = m_D_blocks[k-1];
      real_type * D1 = m_D_blocks[k];
      // Work = L^T nr0 x nr1
      getranspose( nr1, nr0, L0, nr1, m_Work, nr0 );
      // solve
      info = getrs( Transposition::YES, nr0, nr1, D0, nr0, m_B_permutation[k-1], m_Work, nr0 );
      UTILS_ASSERT(
        info == 0,
        "BlockTridiagonalSymmetic::factorize[{}] getrs INFO = {}\n", who, info
      );

      // DD{k}   = DD{k} - (LL{k-1}*DD{k-1}) *LL{k-1}.';

      gemm(
        Transposition::YES,
        Transposition::YES,
        nr1, nr1, nr0,
        -1.0, m_Work, nr0,
        L0, nr1,
        1.0, D1, nr1
      );

      // Work --> L nr1 x nr0
      getranspose( nr0, nr1, m_Work, nr0, L0, nr1 );

      info = getrf( nr1, nr1, D1, nr1, m_B_permutation[k] );

      UTILS_ASSERT(
        info == 0,
        "BlockTridiagonalSymmetic::factorize[{}] getrf INFO = {}\n", who, info
      );

    }
    m_is_factorized = true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  BlockTridiagonalSymmetic<T>::factorize() {
    if ( m_is_factorized ) return true;

    integer info = lapack_wrapper::getrf(
      this->D_nrows(0), this->D_ncols(0),
      m_D_blocks[0], this->D_nrows(0),
      m_B_permutation[0]
    );

    if ( info != 0 ) return false;

    for ( integer k=1; k < m_nBlocks; ++k ) {
      integer nr0 = this->D_nrows(k-1);
      integer nr1 = this->D_nrows(k);
      real_type * L0 = m_L_blocks[k-1];
      real_type * D0 = m_D_blocks[k-1];
      real_type * D1 = m_D_blocks[k];
      // Work = L^T nr0 x nr1
      getranspose( nr1, nr0, L0, nr1, m_Work, nr0 );
      // solve
      info = getrs( Transposition::YES, nr0, nr1, D0, nr0, m_B_permutation[k-1], m_Work, nr0 );
      if ( info != 0 ) return false;

      // DD{k}   = DD{k} - (LL{k-1}*DD{k-1}) *LL{k-1}.';

      gemm(
        Transposition::YES,
        Transposition::YES,
        nr1, nr1, nr0,
        -1.0, m_Work, nr0,
        L0, nr1,
        1.0, D1, nr1
      );

      // Work --> L nr1 x nr0
      getranspose( nr0, nr1, m_Work, nr0, L0, nr1 );

      info = getrf( nr1, nr1, D1, nr1, m_B_permutation[k] );

      if ( info != 0 ) return false;
    }
    return m_is_factorized;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  BlockTridiagonalSymmetic<T>::solve( real_type xb[] ) const {
    if ( !m_is_factorized ) return false;

    // RR{k} = RR{k}-LL{k-1}*RR{k-1};
    integer k   = 0;
    integer nr0 = this->D_nrows(0), nr1;
    real_type * xkm1 = xb, *xk;
    while ( ++k < m_nBlocks ) {
      nr1 = this->D_nrows(k);
      real_type const * L0 = m_L_blocks[k-1];
      xk = xkm1 + nr0;
      gemv(
        Transposition::NO, nr1, nr0,
        -1.0, L0, nr1,
        xkm1, 1,
        1.0,
        xk, 1
      );
      xkm1 = xk;
      nr0  = nr1;
    }
    // RR{k} = DD{k}\RR{k};
    xk = xb;
    for ( k = 0; k < m_nBlocks; ++k ) {
      nr1 = this->D_nrows(k);
      // solve
      real_type const * D1 = m_D_blocks[k];
      integer info = getrs(
        Transposition::NO,
        nr1, 1, D1, nr1, m_B_permutation[k], xk, nr1
      );
      if ( info != 0 ) return false;
      xk += nr1;
    }
    // RR{k} = RR{k}-LL{k}.'*RR{k+1};
    nr1 = this->D_nrows(k-1);
    xk -= nr1;
    while ( --k > 0 ) {
      nr0 = this->D_nrows(k-1);
      real_type const * L0 = m_L_blocks[k-1];
      xkm1 = xk - nr0;
      gemv(
        Transposition::YES, nr1, nr0,
        -1.0, L0, nr1,
        xk, 1,
        1.0,
        xkm1, 1
      );
      xk  = xkm1;
      nr1 = nr0;
    }
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  BlockTridiagonalSymmetic<T>::solve(
    integer   nrhs,
    real_type B[],
    integer   ldB
  ) const {

    if ( !m_is_factorized ) return false;

    // RR{k} = RR{k}-LL{k-1}*RR{k-1};
    integer k   = 0;
    integer nr0 = this->D_nrows(0), nr1;
    real_type * Bkm1 = B, *Bk;
    while ( ++k < m_nBlocks ) {
      nr1 = this->D_nrows(k);
      real_type const * L0 = m_L_blocks[k-1];
      Bk = Bkm1 + nr0;
      gemm(
        Transposition::NO,
        Transposition::NO,
        nr1, nrhs, nr0,
        -1.0, L0, nr1,
        Bkm1, ldB,
        1.0,
        Bk, ldB
      );
      Bkm1 = Bk;
      nr0  = nr1;
    }
    // RR{k} = DD{k}\RR{k};
    Bk = B;
    for ( k = 0; k < m_nBlocks; ++k ) {
      nr1 = this->D_nrows(k);
      // solve
      real_type const * D1 = m_D_blocks[k];
      integer info = getrs(
        Transposition::NO,
        nr1, nrhs,
        D1, nr1, m_B_permutation[k],
        Bk, ldB
      );
      if ( info != 0 ) return false;
      Bk += nr1;
    }
    // RR{k} = RR{k}-LL{k}.'*RR{k+1};
    nr1 = this->D_nrows(k-1);
    Bk -= nr1;
    while ( --k > 0 ) {
      nr0 = this->D_nrows(k-1);
      real_type const * L0 = m_L_blocks[k-1];
      Bkm1 = Bk - nr0;
      gemm(
        Transposition::YES,
        Transposition::NO,
        nr0, nrhs, nr1,
        -1.0, L0, nr1,
        Bk, ldB,
        1.0,
        Bkm1, ldB
      );
      Bk  = Bkm1;
      nr1 = nr0;
    }
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  BlockTridiagonalSymmetic<T>::t_solve( real_type xb[] ) const
  { return BlockTridiagonalSymmetic<T>::solve( xb ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  BlockTridiagonalSymmetic<T>::t_solve(
    integer   nrhs,
    real_type xb[],
    integer   ldXB
  ) const {
    return BlockTridiagonalSymmetic<T>::solve( nrhs, xb, ldXB );
  }

}

///
/// eof: block_trid.cxx
///
