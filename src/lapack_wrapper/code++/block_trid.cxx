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
    allocIntegers.allocate( N + nblks+1 );
    allocRpointers.allocate( 2*nblks-1 );
    allocIpointers.allocate( nblks );
    this->row_blocks    = allocIntegers( nblks+1 );
    this->D_blocks      = allocRpointers( nblks );
    this->L_blocks      = allocRpointers( nblks-1 );
    this->B_permutation = allocIpointers( nblks );
    // evalute the memry usage for the L and D blocks
    integer nr0 = rBlocks[1] - rBlocks[0];
    nrmax = nr0;
    this->nnz = nr0*nr0;

    LAPACK_WRAPPER_ASSERT(
      nr0 >= 0,
      "BlockTridiagonalSymmetic::setup, nr0 = " << nr0
    );

    for ( integer i = 1; i < nblks; ++i ) {
      integer nr = rBlocks[i+1] - rBlocks[i];
      LAPACK_WRAPPER_ASSERT(
        nr >= 0,
        "BlockTridiagonalSymmetic::setup, nr = " << nr
      );
      this->nnz += nr*(nr+nr0);
      if ( nr > nrmax ) nrmax = nr;
      nr0 = nr;
    }
    allocReals.allocate( this->nnz + nrmax * nrmax );
    nr0 = rBlocks[1] - rBlocks[0];
    this->D_blocks[0] = allocReals( nr0 * nr0 );
    B_permutation[0] = allocIntegers( nr0 );
    for ( integer i = 1; i < nblks; ++i ) {
      integer nr = rBlocks[i+1] - rBlocks[i];
      this->D_blocks[i]   = allocReals( nr * nr );
      this->L_blocks[i-1] = allocReals( nr * nr0 );
      nr0 = nr;
      this->B_permutation[i] = allocIntegers( nr );
    }

    Work = allocReals( nrmax * nrmax );

    this->zero();
    this->nBlocks = nblks;
    std::copy( rBlocks, rBlocks+nblks+1, this->row_blocks );
    is_factorized = false;
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
    integer n = 0;
    for ( integer k = 0; k <= nblks; ++k, n += block_size )
      rBlocks.push_back(n);
    this->setup( nblks, &rBlocks.front() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::zero() { // fill to 0 all the blocks
    lapack_wrapper::zero( this->nnz, D_blocks[0], 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  T const &
  BlockTridiagonalSymmetic<T>::operator () ( integer ii, integer jj ) const {
    integer const * ib = this->row_blocks;
    integer const * ie = this->row_blocks+this->nBlocks+1;
    integer iBlock = integer(std::upper_bound( ib, ie, ii ) - ib) - 1;
    integer jBlock = integer(std::upper_bound( ib, ie, jj ) - ib) - 1;

    LAPACK_WRAPPER_ASSERT(
      row_blocks[iBlock] <= ii && ii < row_blocks[iBlock+1],
      "bad iBlock"
    );
    LAPACK_WRAPPER_ASSERT(
      row_blocks[jBlock] <= jj && jj < row_blocks[jBlock+1],
      "bad iBlock"
    );

    integer nr = row_blocks[iBlock+1] - row_blocks[iBlock];
    integer i  = ii - row_blocks[iBlock];
    integer j  = jj - row_blocks[jBlock];
    LAPACK_WRAPPER_ASSERT(
      jBlock == iBlock || jBlock+1 == iBlock,
      "BlockTridiagonalSymmetic:find( " << ii << " , " << jj <<
      " ) --> ( iBlock = " << iBlock << ", jBlock = " << jBlock <<
      " ) --> ( i = " << i << ", j = " << j << " ) out of range"
    );

    integer ij = i+j*nr;
    if ( iBlock == jBlock ) {
      return this->D_blocks[jBlock][ij];
    } else {
      return this->L_blocks[jBlock][ij];
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  T &
  BlockTridiagonalSymmetic<T>::operator () ( integer ii, integer jj ) {
    integer const * ib = this->row_blocks;
    integer const * ie = this->row_blocks+this->nBlocks+1;
    integer iBlock = integer(std::upper_bound( ib, ie, ii ) - ib) - 1;
    integer jBlock = integer(std::upper_bound( ib, ie, jj ) - ib) - 1;

    LAPACK_WRAPPER_ASSERT(
      row_blocks[iBlock] <= ii && ii < row_blocks[iBlock+1],
      "bad iBlock"
    );
    LAPACK_WRAPPER_ASSERT(
      row_blocks[jBlock] <= jj && jj < row_blocks[jBlock+1],
      "bad iBlock"
    );

    integer nr = row_blocks[iBlock+1] - row_blocks[iBlock];
    integer i  = ii - row_blocks[iBlock];
    integer j  = jj - row_blocks[jBlock];
    LAPACK_WRAPPER_ASSERT(
      jBlock == iBlock || jBlock+1 == iBlock,
      "BlockTridiagonalSymmetic:operator () ( " << ii << " , " << jj <<
      " ) --> ( iBlock = " << iBlock << ", jBlock = " << jBlock <<
      " ) --> ( i = " << i << ", j = " << j << " ) out of range"
    );

    integer ij = i+j*nr;
    if ( iBlock == jBlock ) {
      return this->D_blocks[jBlock][ij];
    } else {
      return this->L_blocks[jBlock][ij];
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::insert(
    integer   ii,
    integer   jj,
    valueType v,
    bool      sym
  ) {
    integer const * ib = this->row_blocks;
    integer const * ie = this->row_blocks+this->nBlocks+1;
    integer jBlock = integer(std::upper_bound( ib, ie, jj ) - ib) - 1;
    integer jmin   = row_blocks[jBlock];
    integer jmax   = row_blocks[jBlock+1];
    integer j      = jj - jmin;
    LAPACK_WRAPPER_ASSERT(
      jmin <= jj && jj < jmax && jmin <= ii,
      "bad iBlock ii = " << ii << " jj = " << jj
    );
    if ( ii < jmax ) {
      integer i  = ii - jmin;
      integer nr = jmax-jmin;
      valueType * D = this->D_blocks[jBlock];
      D[i+j*nr] = v;
      if ( sym ) D[j+i*nr] = v;
    } else {
      integer imin = jmax;
      integer imax = row_blocks[jBlock+2];
      integer nr   = imax-imin;
      integer i    = ii - imin;
      this->L_blocks[jBlock][i+j*nr] = v;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::setD(
    integer         n,
    valueType const data[],
    integer         ldData,
    bool            transposed
  ) {
    integer nr = this->DnumRows(n);
    if ( transposed ) {
      getranspose( nr, nr, data, ldData, this->D_blocks[n], nr );
    } else {
      integer ierr = gecopy( nr, nr, data, ldData, this->D_blocks[n], nr );
      LAPACK_WRAPPER_ASSERT(
        ierr == 0,
        "BlockTridiagonalSymmetic::setD, gecopy return ierr = " << ierr
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::setL(
    integer         n,
    valueType const data[],
    integer         ldData,
    bool            transposed
  ) {
    integer nr = this->LnumRows(n);
    integer nc = this->LnumCols(n);
    if ( transposed ) {
      getranspose( nr, nc, data, ldData, this->L_blocks[n], nr );
    } else {
      integer ierr = gecopy( nr, nc, data, ldData, this->L_blocks[n], nr );
      LAPACK_WRAPPER_ASSERT(
        ierr == 0,
        "BlockTridiagonalSymmetic::setL, gecopy return ierr = " << ierr
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::setD(
    integer           n,
    valueType const * data,
    integer           ldData,
    integer           beginRow,
    integer           beginCol,
    integer           nrow,
    integer           ncol,
    bool              transposed
  ) {
    integer ldD = this->DnumRows(n);
    valueType * D = D_blocks[n]+beginRow+beginCol*ldD;

    if ( transposed ) {
      getranspose( nrow, ncol, data, ldData, D, ldD );
    } else {
      integer ierr = gecopy( nrow, ncol, data, ldData, D, ldD );
      LAPACK_WRAPPER_ASSERT(
        ierr == 0,
        "BlockTridiagonalSymmetic::setD (block), gecopy return ierr = " << ierr
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::setL(
    integer           n,
    valueType const * data,
    integer           ldData,
    integer           beginRow,
    integer           beginCol,
    integer           nrow,
    integer           ncol,
    bool              transposed
  ) {
    integer ldL = this->LnumRows(n);
    valueType * L = L_blocks[n]+beginRow+beginCol*ldL;
    if ( transposed ) {
      getranspose( nrow, ncol, data, ldData, L, ldL );
    } else {
      integer ierr = gecopy( nrow, ncol, data, ldData, L, ldL );
      LAPACK_WRAPPER_ASSERT(
        ierr == 0,
        "BlockTridiagonalSymmetic::setL (block), gecopy return ierr = " << ierr
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::factorize( char const who[] ) {
    LAPACK_WRAPPER_ASSERT(
      !is_factorized,
      "BlockTridiagonalSymmetic::factorize[" << who <<
      "], already factored"
    );

    integer info = lapack_wrapper::getrf(
      this->DnumRows(0), this->DnumCols(0),
      this->D_blocks[0], this->DnumRows(0),
      B_permutation[0]
    );

    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "BlockTridiagonalSymmetic::factorize[" << who <<
      "] getrf INFO = " << info
    );

    for ( integer k=1; k < this->nBlocks; ++k ) {
      integer nr0 = this->DnumRows(k-1);
      integer nr1 = this->DnumRows(k);
      valueType * L0 = this->L_blocks[k-1];
      valueType * D0 = this->D_blocks[k-1];
      valueType * D1 = this->D_blocks[k];
      // Work = L^T nr0 x nr1
      getranspose( nr1, nr0, L0, nr1, Work, nr0 );
      // solve
      info = getrs( TRANSPOSE, nr0, nr1, D0, nr0, B_permutation[k-1], Work, nr0 );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "BlockTridiagonalSymmetic::factorize[" << who <<
        "] getrs INFO = " << info
      );

      // DD{k}   = DD{k} - (LL{k-1}*DD{k-1}) *LL{k-1}.';

      gemm(
        TRANSPOSE,
        TRANSPOSE,
        nr1, nr1, nr0,
        -1.0, Work, nr0,
        L0, nr1,
        1.0, D1, nr1
      );

      // Work --> L nr1 x nr0
      getranspose( nr0, nr1, Work, nr0, L0, nr1 );

      info = getrf( nr1, nr1, D1, nr1, B_permutation[k] );

      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "BlockTridiagonalSymmetic::factorize[" << who <<
        "] getrf INFO = " << info
      );

    }
    is_factorized = true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::solve( valueType xb[] ) const {
    LAPACK_WRAPPER_ASSERT(
      is_factorized,
      "BlockTridiagonalSymmetic::solve, matrix not factored"
    );

    // RR{k} = RR{k}-LL{k-1}*RR{k-1};
    integer k = 0;
    integer nr0 = this->DnumRows(0), nr1;
    valueType * xkm1 = xb, *xk;
    while ( ++k < this->nBlocks ) {
      nr1 = this->DnumRows(k);
      valueType const * L0 = this->L_blocks[k-1];
      xk = xkm1 + nr0;
      gemv(
        NO_TRANSPOSE, nr1, nr0,
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
    for ( k = 0; k < this->nBlocks; ++k ) {
      nr1 = this->DnumRows(k);
      // solve
      valueType const * D1 = this->D_blocks[k];
      integer info = getrs( NO_TRANSPOSE, nr1, 1, D1, nr1, B_permutation[k], xk, nr1 );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "BlockTridiagonalSymmetic::solve getrs INFO = " << info
      );
      xk += nr1;
    }
    // RR{k} = RR{k}-LL{k}.'*RR{k+1};
    nr1 = this->DnumRows(k-1);
    xk -= nr1;
    while ( --k > 0 ) {
      nr0 = this->DnumRows(k-1);
      valueType const * L0 = this->L_blocks[k-1];
      xkm1 = xk - nr0;
      gemv(
        TRANSPOSE, nr1, nr0,
        -1.0, L0, nr1,
        xk, 1,
        1.0,
        xkm1, 1
      );
      xk  = xkm1;
      nr1 = nr0;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::solve(
    integer   nrhs,
    valueType B[],
    integer   ldB
  ) const {
    LAPACK_WRAPPER_ASSERT(
      is_factorized,
      "BlockTridiagonalSymmetic::solve, matrix not factored"
    );

    // RR{k} = RR{k}-LL{k-1}*RR{k-1};
    integer k = 0;
    integer nr0 = this->DnumRows(0), nr1;
    valueType * Bkm1 = B, *Bk;
    while ( ++k < this->nBlocks ) {
      nr1 = this->DnumRows(k);
      valueType const * L0 = this->L_blocks[k-1];
      Bk = Bkm1 + nr0;
      gemm(
        NO_TRANSPOSE, NO_TRANSPOSE, nr1, nrhs, nr0,
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
    for ( k = 0; k < this->nBlocks; ++k ) {
      nr1 = this->DnumRows(k);
      // solve
      valueType const * D1 = this->D_blocks[k];
      integer info = getrs(
        NO_TRANSPOSE, nr1, nrhs,
        D1, nr1, B_permutation[k],
        Bk, ldB
      );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "BlockTridiagonalSymmetic::solve getrs INFO = " << info
      );
      Bk += nr1;
    }
    // RR{k} = RR{k}-LL{k}.'*RR{k+1};
    nr1 = this->DnumRows(k-1);
    Bk -= nr1;
    while ( --k > 0 ) {
      nr0 = this->DnumRows(k-1);
      valueType const * L0 = this->L_blocks[k-1];
      Bkm1 = Bk - nr0;
      gemm(
        TRANSPOSE, NO_TRANSPOSE, nr0, nrhs, nr1,
        -1.0, L0, nr1,
        Bk, ldB,
        1.0,
        Bkm1, ldB
      );
      Bk  = Bkm1;
      nr1 = nr0;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::t_solve( valueType xb[] ) const
  { BlockTridiagonalSymmetic<T>::solve( xb ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::t_solve(
    integer   nrhs,
    valueType xb[],
    integer   ldXB
  ) const {
    BlockTridiagonalSymmetic<T>::solve( nrhs, xb, ldXB );
  }

}

///
/// eof: block_trid.cxx
///
