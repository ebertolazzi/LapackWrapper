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
/// file: sparse.cxx
///

namespace lapack_wrapper {

  template <typename T>
  void
  SparseMatrixBase<T>::print( ostream_type & stream ) const {
    integer const * iRow;
    integer const * jCol;
    real    const * vals;
    this->data(iRow, jCol, vals);
    for ( size_t i = 0; i < size_t(this->get_nnz()); ++i)
      fmt::print( stream, "{}, {}, {}\n", iRow[i], jCol[i], vals[i] );
  }

  template <typename T>
  void
  SparseMatrixBase<T>::push_matrix(
    integer      row_offs,
    integer      col_offs,
    MatW const & Matrix,
    bool         transpose,
    integer      lower_upper
  ) {
    /*\
     * lower_upper 0 = full matrix
     *             1 = upper part with diagonal
     *             2 = upper part without diagonal
     *             3 = push symmetric, for each A(i,j) element an A(j,i) is inserted
     *            -3 = push anty symmetric, for each A(i,j) element an -A(j,i) is inserted
     *            -1 = lower part with diagonal
     *            -2 = lower part without diagonal
    \*/
    switch ( lower_upper ) {
    case 0:
    case 1:
    case 2:
    case -1:
    case -2:
      for ( integer j = 0; j < Matrix.numCols(); ++j ) {
        for ( integer i = 0; i < Matrix.numRows(); ++i ) {
          integer ii = i;
          integer jj = j;
          if ( transpose ) std::swap( ii, jj );
          bool do_push = true;
          switch ( lower_upper ) {
            case  1: do_push = j >= i; break;
            case  2: do_push = j  > i; break;
            case -1: do_push = j <= i; break;
            case -2: do_push = j  < i; break;
          }
          if ( do_push ) {
            ii += row_offs;
            jj += col_offs;
            this->push_value_C( ii, jj, Matrix(i,j) );
          }
        }
      }
      break;
    case 3:
      for ( integer j = 0; j < Matrix.numCols(); ++j ) {
        for ( integer i = 0; i < Matrix.numRows(); ++i ) {
          integer ii = i;
          integer jj = j;
          if ( transpose ) std::swap( ii, jj );
          T const & v = Matrix(i,j);
          this->push_value_C( ii, jj, v );
          if ( i != j ) this->push_value_C( jj, ii, v );
        }
      }
      break;
    case -3:
      for ( integer j = 0; j < Matrix.numCols(); ++j ) {
        for ( integer i = 0; i < Matrix.numRows(); ++i ) {
          integer ii = i;
          integer jj = j;
          if ( transpose ) std::swap( ii, jj );
          T const & v = Matrix(i,j);
          this->push_value_C( ii, jj, v );
          if ( i != j ) this->push_value_C( jj, ii, -v );
        }
      }
      break;
    }
  }

  template <typename T>
  void
  SparseMatrixBase<T>::push_matrix(
    integer        row_offs,
    integer        col_offs,
    Sparse const & Matrix,
    bool           transpose,
    integer        lower_upper
  ) {
    /*\
     * lower_upper 0 = full matrix
     *             1 = upper part with diagonal
     *             2 = upper part without diagonal
     *             3 = push symmetric, for each A(i,j) element an A(j,i) is inserted
     *            -3 = push anty symmetric, for each A(i,j) element an -A(j,i) is inserted
     *            -1 = lower part with diagonal
     *            -2 = lower part without diagonal
    \*/
    integer const * rowsM;
    integer const * colsM;
    T       const * valsM;
    Matrix.get_data( rowsM, colsM, valsM );
    if ( transpose ) std::swap( rowsM, colsM );
    if ( Matrix.FORTRAN_indexing() ) { --row_offs; --col_offs; }

    switch ( lower_upper ) {
    case 0:
    case 1:
    case 2:
    case -1:
    case -2:
      for ( integer index = 0; index < Matrix.get_nnz(); ++index ) {
        integer i = rowsM[index];
        integer j = colsM[index];
        T const & v = valsM[index];
        bool do_push = true;
        switch ( lower_upper ) {
          case  1: do_push = j >= i; break;
          case  2: do_push = j  > i; break;
          case -1: do_push = j <= i; break;
          case -2: do_push = j  < i; break;
        }
        if ( do_push )
          this->push_value_C( row_offs+i, col_offs+j, v );
      }
      break;
    case 3:
      for ( integer index = 0; index < Matrix.get_nnz(); ++index ) {
        integer i = rowsM[index];
        integer j = colsM[index];
        T const & v = valsM[index];
        this->push_value_C( row_offs+i, col_offs+j, v );
        if ( i != j ) this->push_value_C( row_offs+j, col_offs+i, v );
      }
      break;
    case -3:
      for ( integer index = 0; index < Matrix.get_nnz(); ++index ) {
        integer i = rowsM[index];
        integer j = colsM[index];
        T const & v = valsM[index];
        this->push_value_C( row_offs+i, col_offs+j, v );
        if ( i != j ) this->push_value_C( row_offs+j, col_offs+i, -v );
      }
      break;
    }
  }

  /*\
  :|:   ____                             ____ ____ ___   ___  ____
  :|:  / ___| _ __   __ _ _ __ ___  ___ / ___/ ___/ _ \ / _ \|  _ \
  :|:  \___ \| '_ \ / _` | '__/ __|/ _ \ |  | |  | | | | | | | |_) |
  :|:   ___) | |_) | (_| | |  \__ \  __/ |__| |__| |_| | |_| |  _ <
  :|:  |____/| .__/ \__,_|_|  |___/\___|\____\____\___/ \___/|_| \_\
  :|:        |_|
  \*/

  template <typename T>
  SparseCCOOR<T>::SparseCCOOR(
    integer N,
    integer M,
    integer reserve_nnz,
    bool    fi
  )
  : SparseMatrixBase<T>(N,M)
  , m_fortran_indexing(fi)
  , m_matrix_is_full(false)
  , m_matrix_is_row_major(false)
  {
    m_vals.clear(); m_vals.reserve(reserve_nnz);
    m_rows.clear(); m_rows.reserve(reserve_nnz);
    m_cols.clear(); m_cols.reserve(reserve_nnz);
  }

  //! Sparse Matrix Structure
  template <typename T>
  bool
  SparseCCOOR<T>::foundNaN() const
  { return Utils::foundNaN( &m_vals.front(), m_nnz ); }

  template <typename T>
  void
  SparseCCOOR<T>::gemv(
    real_type       alpha,
    integer         DimX,
    real_type const x[],
    integer         incX,
    real_type       beta,
    integer         DimY,
    real_type       y[],
    integer         incY
  ) const {
    UTILS_ASSERT(
      DimX == m_ncols && DimY == m_nrows,
      "SparseCCOOR::gemv, bad dimensions, "
      "dimX = {}, dimY = {} matrix is {} x {}\n",
      DimX, DimY, m_nrows, m_ncols
    );
    if ( m_matrix_is_full ) {
      if ( m_matrix_is_row_major ) {
        lapack_wrapper::gemv(
          TRANSPOSE, m_ncols, m_nrows,
          alpha,
          &m_vals.front(), m_ncols,
          x, 1,
          beta,
          y, 1
        );
      } else {
        lapack_wrapper::gemv(
          NO_TRANSPOSE, m_nrows, m_ncols,
          alpha,
          &m_vals.front(), m_nrows,
          x, 1,
          beta,
          y, 1
        );
      }
    } else {
      this->y_manage( beta, DimY, y, incY );
      integer offs = m_fortran_indexing ? -1 : 0;
      for ( integer idx = 0; idx < m_nnz; ++idx ) {
        integer i = m_rows[idx] + offs;
        integer j = m_cols[idx] + offs;
        y[ i * incY ] += alpha * m_vals[idx] * x[ j * incX ];
      }
    }
  }

  template <typename T>
  void
  SparseCCOOR<T>::gemv_Transposed(
    real_type       alpha,
    integer         DimX,
    real_type const x[],
    integer         incX,
    real_type       beta,
    integer         DimY,
    real_type       y[],
    integer         incY
  ) const {
    UTILS_ASSERT(
      DimY == m_ncols && DimX == m_nrows,
      "SparseCCOOR::gemv_Transposed, bad dimensions, "
      "dimX = {}, dimY = {} matrix is {} x {}\n",
      DimX, DimY, m_nrows, m_ncols
    );
    if ( m_matrix_is_full ) {
      if ( m_matrix_is_row_major ) {
        lapack_wrapper::gemv(
          NO_TRANSPOSE, m_ncols, m_nrows,
          alpha,
          &m_vals.front(), m_ncols,
          x, 1,
          beta,
          y, 1
        );
      } else {
        lapack_wrapper::gemv(
          TRANSPOSE, m_nrows, m_ncols,
          alpha,
          &m_vals.front(), m_nrows,
          x, 1,
          beta,
          y, 1
        );
      }
    } else {
      this->y_manage( beta, DimY, y, incY );
      integer offs = m_fortran_indexing ? -1 : 0;
      for ( integer idx = 0; idx < m_nnz; ++idx ) {
        integer j = m_rows[idx] + offs;
        integer i = m_cols[idx] + offs;
        y[ i * incY ] += alpha * m_vals[idx] * x[ j * incX ];
      }
    }
  }

  template <typename T>
  void
  SparseCCOOR<T>::gemv_Symmetric(
    real_type       alpha,
    integer         DimX,
    real_type const x[],
    integer         incX,
    real_type       beta,
    integer         DimY,
    real_type       y[],
    integer         incY
  ) const {
    UTILS_ASSERT(
      DimY == m_ncols && DimX == m_nrows && DimX == DimY,
      "SparseCCOOR::gemv_Symmetric, bad dimensions, "
      "dimX = {}, dimY = {} matrix is {} x {}\n",
      DimX, DimY, m_nrows, m_ncols
    );
    this->y_manage( beta, DimY, y, incY );
    integer offs = m_fortran_indexing ? -1 : 0;
    for ( integer idx = 0; idx < m_nnz; ++idx ) {
      integer   i   = m_rows[idx] + offs;
      integer   j   = m_cols[idx] + offs;
      real_type tmp = alpha * m_vals[idx];
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
  void
  SparseCCOOR<T>::clear() {
    m_nrows               = 0;
    m_ncols               = 0;
    m_nnz                 = 0;
    m_matrix_is_full      = false;
    m_matrix_is_row_major = false;
    m_vals.clear();
    m_rows.clear();
    m_cols.clear();
  }

  template <typename T>
  void
  SparseCCOOR<T>::setZero() {
    if ( m_matrix_is_full ) {
      std::fill( m_vals.begin(), m_vals.end(), T(0) );
    } else {
      m_nnz = 0;
      m_vals.clear();
      m_rows.clear();
      m_cols.clear();
    }
  }

  template <typename T>
  void
  SparseCCOOR<T>::init(
    integer N,
    integer M,
    integer reserve_nnz,
    bool    fi
  ) {
    m_fortran_indexing    = fi;
    m_nrows               = N;
    m_ncols               = M;
    m_nnz                 = 0;
    m_matrix_is_full      = false;
    m_matrix_is_row_major = false;
    m_vals.clear(); m_vals.reserve(reserve_nnz);
    m_rows.clear(); m_rows.reserve(reserve_nnz);
    m_cols.clear(); m_cols.reserve(reserve_nnz);
  }

  template <typename T>
  void
  SparseCCOOR<T>::setup_as_full_column_major( integer N, integer M, bool fi ) {
    m_fortran_indexing = fi;
    m_nrows            = N;
    m_ncols            = M;
    m_nnz              = N*M;
    m_vals.resize( size_t(m_nnz) );
    m_cols.clear(); m_cols.reserve( size_t(m_nnz) );
    m_rows.clear(); m_rows.reserve( size_t(m_nnz) );
    integer offs = fi ? 1 : 0;
    for ( integer j = 0; j < M; ++j ) {
      for ( integer i = 0; i < N; ++i ) {
        m_rows.push_back( i+offs );
        m_cols.push_back( j+offs );
      }
    }
    m_matrix_is_full      = true;
    m_matrix_is_row_major = false;
  }

  template <typename T>
  void
  SparseCCOOR<T>::setup_as_full_row_major( integer N, integer M, bool fi ) {
    m_fortran_indexing = fi;
    m_nrows            = N;
    m_ncols            = M;
    m_nnz              = N*M;
    m_vals.resize( size_t(m_nnz) );
    m_cols.clear(); m_cols.reserve( size_t(m_nnz) );
    m_rows.clear(); m_rows.reserve( size_t(m_nnz) );
    integer offs = fi ? 1 : 0;
    for ( integer i = 0; i < N; ++i ) {
      for ( integer j = 0; j < M; ++j ) {
        m_rows.push_back( i+offs );
        m_cols.push_back( j+offs );
      }
    }
    m_matrix_is_full      = true;
    m_matrix_is_row_major = true;
  }

  template <typename T>
  void
  SparseCCOOR<T>::fill( real_type const V[], integer M ) {
    UTILS_ASSERT0(
      M == m_nnz,
      "SparseCCOOR::fill(...) bad size input vector\n"
    );
    for ( integer i = 0; i < m_nnz; ++i )
      m_vals[i] = V[i];
  }

  template <typename T>
  void
  SparseCCOOR<T>::fill( std::vector<real_type> const & V ) {
    UTILS_ASSERT0(
      V.size() == m_vals.size(),
      "SparseCCOOR::fill(...) bad size input vector\n"
    );
    std::copy( V.begin(), V.end(), m_vals.begin() );
  }

  template <typename T>
  void
  SparseCCOOR<T>::reserve( integer reserve_nnz ) {
    m_vals.reserve(reserve_nnz);
    m_rows.reserve(reserve_nnz);
    m_cols.reserve(reserve_nnz);
  }

  template <typename T>
  void
  SparseCCOOR<T>::to_FORTRAN_indexing() {
    if ( !m_fortran_indexing ) {
      std::vector<integer>::iterator iv;
      for ( iv = m_rows.begin(); iv != m_rows.end(); ++iv ) ++(*iv);
      for ( iv = m_cols.begin(); iv != m_cols.end(); ++iv ) ++(*iv);
      m_fortran_indexing = true;
    }
  }

  template <typename T>
  void
  SparseCCOOR<T>::to_C_indexing() {
    if ( m_fortran_indexing ) {
      std::vector<integer>::iterator iv;
      for ( iv = m_rows.begin(); iv != m_rows.end(); ++iv ) --(*iv);
      for ( iv = m_cols.begin(); iv != m_cols.end(); ++iv ) --(*iv);
      m_fortran_indexing = false;
    }
  }

  template <typename T>
  void
  SparseCCOOR<T>::push_value_C(
    integer   row,
    integer   col,
    real_type val
  ) {
    UTILS_ASSERT(
      row >= 0 && row < m_nrows && col >= 0 && col < m_ncols,
      "SparseCCOOR::push_value_C({},{}) out of bound\n", row, col
    );
    if ( m_matrix_is_full ) {
       if ( m_matrix_is_row_major ) m_vals[ col + row * m_ncols ] =  val;
       else                         m_vals[ row + col * m_nrows ] =  val;
    } else {
      if ( m_fortran_indexing ) { ++row; ++col; }
      m_vals.push_back(val);
      m_rows.push_back(row);
      m_cols.push_back(col);
      ++m_nnz;
    }
  }

  template <typename T>
  void
  SparseCCOOR<T>::push_value_F(
    integer   row,
    integer   col,
    real_type val
  ) {
    UTILS_ASSERT(
      row > 0 && row <= m_nrows && col > 0 && col <= m_ncols,
      "SparseCCOOR::push_value_F({},{}) out of bound\n", row, col
    );
    if ( !m_fortran_indexing ) { --row; --col; }
    if ( m_matrix_is_full ) {
       if ( m_matrix_is_row_major ) m_vals[ col + row * m_ncols ] =  val;
       else                         m_vals[ row + col * m_nrows ] =  val;
    } else {
      m_vals.push_back(val);
      m_rows.push_back(row);
      m_cols.push_back(col);
      ++m_nnz;
    }
  }

  template <typename T>
  void
  SparseCCOOR<T>::get_matrix( MatW & M ) const {
    M.zero_fill();
    for ( integer index = 0; index < m_nnz; ++index ) {
      integer i = m_rows[index];
      integer j = m_cols[index];
      if ( m_fortran_indexing ) { --i; --j; }
      M(i,j) = m_vals[index];
    }
  }

  template <typename T>
  void
  SparseCCOOR<T>::get_matrix_symmetric( MatW & M ) const {
    M.zero_fill();
    for ( integer index = 0; index < m_nnz; ++index ) {
      integer i = m_rows[index];
      integer j = m_cols[index];
      if ( m_fortran_indexing ) { --i; --j; }
      T v = m_vals[index];
      M(i,j) = v;
      if ( i != j ) M(j,i) = v;
    }
  }

  template <typename T>
  void
  SparseCCOOR<T>::get_matrix_transposed( MatW & M ) const {
    M.zero_fill();
    for ( integer index = 0; index < m_nnz; ++index ) {
      integer i = m_rows[index];
      integer j = m_cols[index];
      if ( m_fortran_indexing ) { --i; --j; }
      M(j,i) = m_vals[index];
    }
  }

  template class SparseCCOOR<real>;
  template class SparseCCOOR<doublereal>;

}

///
/// eof: sparse.cxx
///
