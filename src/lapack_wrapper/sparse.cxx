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

namespace lapack_wrapper {

  /*\
  :|:   ____                             ____ ____ ___   ___  ____
  :|:  / ___| _ __   __ _ _ __ ___  ___ / ___/ ___/ _ \ / _ \|  _ \
  :|:  \___ \| '_ \ / _` | '__/ __|/ _ \ |  | |  | | | | | | | |_) |
  :|:   ___) | |_) | (_| | |  \__ \  __/ |__| |__| |_| | |_| |  _ <
  :|:  |____/| .__/ \__,_|_|  |___/\___|\____\____\___/ \___/|_| \_\
  :|:        |_|
  \*/

  //! Sparse Matrix Structure
  template <typename T>
  bool
  SparseCCOOR<T>::foundNaN() const
  { return lapack_wrapper::foundNaN( &this->vals.front(), this->nnz ); }

  template <typename T>
  void
  SparseCCOOR<T>::gemv(
    valueType       alpha,
    integer         DimX,
    valueType const x[],
    integer         incX,
    valueType       beta,
    integer         DimY,
    valueType       y[],
    integer         incY
  ) const {
    LAPACK_WRAPPER_ASSERT(
      DimX == this->nCols && DimY == this->nRows,
      "SparseCCOOR::gemv, bad dimensions, dimX = " << DimX <<
      ", dimY = " << DimY << " matrix is " << this->nRows <<
      " x " << this->nCols
    );
    this->y_manage( beta, DimY, y, incY );
    integer offs = this->fortran_indexing ? -1 : 0;
    for ( integer idx = 0; idx < this->nnz; ++idx ) {
      integer i = this->rows[idx] + offs;
      integer j = this->cols[idx] + offs;
      y[ i * incY ] += alpha * this->vals[idx] * x[ j * incX ];
    }
  }

  template <typename T>
  void
  SparseCCOOR<T>::gemv_Transposed(
    valueType       alpha,
    integer         DimX,
    valueType const x[],
    integer         incX,
    valueType       beta,
    integer         DimY,
    valueType       y[],
    integer         incY
  ) const {
    LAPACK_WRAPPER_ASSERT(
      DimY == this->nCols && DimX == this->nRows,
      "SparseCCOOR::gemv_Transposed, bad dimensions, dimX = " << DimX <<
      ", dimY = " << DimY << " matrix is " << this->nRows <<
      " x " << this->nCols
    );
    this->y_manage( beta, DimY, y, incY );
    integer offs = this->fortran_indexing ? -1 : 0;
    for ( integer idx = 0; idx < this->nnz; ++idx ) {
      integer j = this->rows[idx] + offs;
      integer i = this->cols[idx] + offs;
      y[ i * incY ] += alpha * this->vals[idx] * x[ j * incX ];
    }
  }

  template <typename T>
  void
  SparseCCOOR<T>::gemv_Symmetric(
    valueType       alpha,
    integer         DimX,
    valueType const x[],
    integer         incX,
    valueType       beta,
    integer         DimY,
    valueType       y[],
    integer         incY
  ) const {
    LAPACK_WRAPPER_ASSERT(
      DimY == this->nCols && DimX == this->nRows && DimX == DimY,
      "SparseCCOOR::gemv_Symmetric, bad dimensions, dimX = " << DimX <<
      ", dimY = " << DimY << " matrix is " << this->nRows <<
      " x " << this->nCols
    );
    this->y_manage( beta, DimY, y, incY );
    integer offs = this->fortran_indexing ? -1 : 0;
    for ( integer idx = 0; idx < this->nnz; ++idx ) {
      integer   i   = this->rows[idx] + offs;
      integer   j   = this->cols[idx] + offs;
      valueType tmp = alpha * this->vals[idx];
      if ( i == j ) {
        y[ i * incY ] += tmp * x[ j * incX ];
      } else {
        y[ i * incY ] += tmp * x[ j * incX ];
        y[ j * incY ] += tmp * x[ i * incX ];
      }
    }
  }

  /*\
  :|:           _    _ _ _   _               _             _   _            _
  :|:   __ _ __| |__| (_) |_(_)___ _ _  __ _| |  _ __  ___| |_| |_  ___  __| |___
  :|:  / _` / _` / _` | |  _| / _ \ ' \/ _` | | | '  \/ -_)  _| ' \/ _ \/ _` (_-<
  :|:  \__,_\__,_\__,_|_|\__|_\___/_||_\__,_|_| |_|_|_\___|\__|_||_\___/\__,_/__/
  \*/

  template <typename T>
  SparseCCOOR<T>::SparseCCOOR(
    integer N,
    integer M,
    integer reserve_nnz,
    bool    fi
  )
  : nRows(N)
  , nCols(M)
  , nnz(reserve_nnz)
  , fortran_indexing(fi)
  {
    this->vals.clear(); this->vals.reserve(reserve_nnz);
    this->rows.clear(); this->rows.reserve(reserve_nnz);
    this->cols.clear(); this->cols.reserve(reserve_nnz);
  }

  template <typename T>
  void
  SparseCCOOR<T>::clear() {
    this->nRows = 0;
    this->nCols = 0;
    this->nnz   = 0;
    this->vals.clear();
    this->rows.clear();
    this->cols.clear();
  }

  template <typename T>
  void
  SparseCCOOR<T>::setZero() {
    this->nnz = 0;
    this->vals.clear();
    this->rows.clear();
    this->cols.clear();
  }

  template <typename T>
  void
  SparseCCOOR<T>::init(
    integer N,
    integer M,
    integer reserve_nnz,
    bool    fi
  ) {
    this->fortran_indexing = fi;
    this->nRows            = N;
    this->nCols            = M;
    this->nnz              = 0;
    this->vals.clear(); this->vals.reserve(reserve_nnz);
    this->rows.clear(); this->rows.reserve(reserve_nnz);
    this->cols.clear(); this->cols.reserve(reserve_nnz);
  }

  template <typename T>
  void
  SparseCCOOR<T>::setup_as_full_column_major( integer N, integer M, bool fi ) {
    this->fortran_indexing = fi;
    this->nRows            = N;
    this->nCols            = M;
    this->nnz              = N*M;
    this->vals.resize( size_t(this->nnz) );
    this->cols.clear(); cols.reserve( size_t(this->nnz) );
    this->rows.clear(); rows.reserve( size_t(this->nnz) );
    integer offs = fi ? 1 : 0;
    for ( integer j = 0; j < M; ++j ) {
      for ( integer i = 0; i < N; ++i ) {
        this->rows.push_back( i+offs );
        this->cols.push_back( j+offs );
      }
    }
  }

  template <typename T>
  void
  SparseCCOOR<T>::setup_as_full_row_major( integer N, integer M, bool fi ) {
    this->fortran_indexing = fi;
    this->nRows            = N;
    this->nCols            = M;
    this->nnz              = N*M;
    this->vals.resize( size_t(this->nnz) );
    this->cols.clear(); cols.reserve( size_t(this->nnz) );
    this->rows.clear(); rows.reserve( size_t(this->nnz) );
    integer offs = fi ? 1 : 0;
    for ( integer i = 0; i < N; ++i ) {
      for ( integer j = 0; j < M; ++j ) {
        this->rows.push_back( i+offs );
        this->cols.push_back( j+offs );
      }
    }
  }

  template <typename T>
  void
  SparseCCOOR<T>::fill( valueType const V[], integer M ) {
    LAPACK_WRAPPER_ASSERT(
      M == this->nnz,
      "SparseCCOOR::fill(...) bad size input vector"
    );
    for ( integer i = 0; i < this->nnz; ++i )
      this->vals[i] = V[i];
  }

  template <typename T>
  void
  SparseCCOOR<T>::fill( std::vector<valueType> const & V ) {
    LAPACK_WRAPPER_ASSERT(
      V.size() == this->vals.size(),
      "SparseCCOOR::fill(...) bad size input vector"
    );
    std::copy( V.begin(), V.end(), this->vals.begin() );
  }

  template <typename T>
  void
  SparseCCOOR<T>::reserve( integer reserve_nnz ) {
    this->vals.reserve(reserve_nnz);
    this->rows.reserve(reserve_nnz);
    this->cols.reserve(reserve_nnz);
  }

  template <typename T>
  void
  SparseCCOOR<T>::to_FORTRAN_indexing() {
    if ( !this->fortran_indexing ) {
      std::vector<integer>::iterator iv;
      for ( iv = rows.begin(); iv != rows.end(); ++iv ) ++(*iv);
      for ( iv = cols.begin(); iv != cols.end(); ++iv ) ++(*iv);
      this->fortran_indexing = true;
    }
  }

  template <typename T>
  void
  SparseCCOOR<T>::to_C_indexing() {
    if ( this->fortran_indexing ) {
      std::vector<integer>::iterator iv;
      for ( iv = rows.begin(); iv != rows.end(); ++iv ) --(*iv);
      for ( iv = cols.begin(); iv != cols.end(); ++iv ) --(*iv);
      this->fortran_indexing = false;
    }
  }

  template <typename T>
  void
  SparseCCOOR<T>::push_value_C(
    integer   row,
    integer   col,
    valueType val
  ) {
    LAPACK_WRAPPER_ASSERT(
      row >= 0 && row < this->nRows && col >= 0 && col < this->nCols,
      "SparseCCOOR::push_value_C( " << row << ", " << col << ") out of bound"
    );
    if ( this->fortran_indexing ) { ++row; ++col; }
    this->vals.push_back(val);
    this->rows.push_back(row);
    this->cols.push_back(col);
    ++this->nnz;
  }

  template <typename T>
  void
  SparseCCOOR<T>::push_value_F(
    integer   row,
    integer   col,
    valueType val
  ) {
    LAPACK_WRAPPER_ASSERT(
      row > 0 && row <= this->nRows && col > 0 && col <= this->nCols,
      "SparseCCOOR::push_value_F( " << row << ", " << col << ") out of bound"
    );
    if ( !this->fortran_indexing ) { --row; --col; }
    this->vals.push_back(val);
    this->rows.push_back(row);
    this->cols.push_back(col);
    ++this->nnz;
  }

  template <typename T>
  void
  SparseCCOOR<T>::push_matrix(
    integer      row_offs,
    integer      col_offs,
    MatW const & Matrix,
    bool         transpose
  ) {
    if ( this->fortran_indexing ) { ++row_offs; ++col_offs; }
    for ( integer j = 0; j < Matrix.numCols(); ++j ) {
      for ( integer i = 0; i < Matrix.numRows(); ++i ) {
        if ( transpose ) {
          this->rows.push_back( row_offs + j );
          this->cols.push_back( col_offs + i );
        } else {
          this->rows.push_back( row_offs + i );
          this->cols.push_back( col_offs + j );
        }
        this->vals.push_back( Matrix(i,j) );
        ++this->nnz;
      }
    }
  }

  template <typename T>
  void
  SparseCCOOR<T>::push_matrix(
    integer        row_offs,
    integer        col_offs,
    Sparse const & Matrix,
    bool           transpose
  ) {
    integer const * rowsM;
    integer const * colsM;
    T       const * valsM;
    Matrix.get_data( rowsM, colsM, valsM );
    if ( transpose ) std::swap( rowsM, colsM );
    if ( Matrix.FORTRAN_indexing() ) { --row_offs; --col_offs; }
    if ( this->fortran_indexing )    { ++row_offs; ++col_offs; }
    for ( integer index = 0; index < Matrix.get_nnz(); ++index ) {
      this->rows.push_back( row_offs + rowsM[index] );
      this->cols.push_back( col_offs + colsM[index] );
      this->vals.push_back( valsM[index] );
      ++this->nnz;
    }
  }


  template <typename T>
  void
  SparseCCOOR<T>::get_matrix( MatW & M ) const {
    M.zero_fill();
    for ( integer index = 0; index < this->nnz; ++index ) {
      integer i = this->rows[index];
      integer j = this->cols[index];
      if ( fortran_indexing ) { --i; --j; }
      M(i,j) = this->vals[index];
    }
  }

  template <typename T>
  void
  SparseCCOOR<T>::get_matrix_symmetric( MatW & M ) const {
    M.zero_fill();
    for ( integer index = 0; index < this->nnz; ++index ) {
      integer i = this->rows[index];
      integer j = this->cols[index];
      if ( fortran_indexing ) { --i; --j; }
      T v = this->vals[index];
      M(i,j) = v;
      if ( i != j ) M(j,i) = v;
    }
  }

  template <typename T>
  void
  SparseCCOOR<T>::get_matrix_transposed( MatW & M ) const {
    M.zero_fill();
    for ( integer index = 0; index < this->nnz; ++index ) {
      integer i = this->rows[index];
      integer j = this->cols[index];
      if ( fortran_indexing ) { --i; --j; }
      M(j,i) = this->vals[index];
    }
  }

  template class SparseCCOOR<real>;
  template class SparseCCOOR<doublereal>;

}

///
/// eof: LapackWrapper_sparse.cxx
///
