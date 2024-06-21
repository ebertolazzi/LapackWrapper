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
/// file: wrapper.cxx
///

namespace lapack_wrapper {

  /*
  //   __  __       _        _     __        __
  //  |  \/  | __ _| |_ _ __(_)_  _\ \      / / __ __ _ _ __  _ __   ___ _ __
  //  | |\/| |/ _` | __| '__| \ \/ /\ \ /\ / / '__/ _` | '_ \| '_ \ / _ \ '__|
  //  | |  | | (_| | |_| |  | |>  <  \ V  V /| | | (_| | |_) | |_) |  __/ |
  //  |_|  |_|\__,_|\__|_|  |_/_/\_\  \_/\_/ |_|  \__,_| .__/| .__/ \___|_|
  //                                                   |_|   |_|
  */

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  template <typename T>
  MatrixWrapper<T>::MatrixWrapper(
    real_type * _data,
    integer     nr,
    integer     nc,
    integer     ld
  )
  : m_nrows(nr)
  , m_ncols(nc)
  , m_ldData(ld)
  , m_data(_data)
  {
    UTILS_ASSERT_DEBUG(
      nr >= 0 && nc >= 0 && m_ldData >= nr,
      "MatrixWrapper( data, nr={}, nc={}, ld={}) bad dimensions\n",
      nr, nc, ld
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::setup(
    real_type * _data,
    integer     nr,
    integer     nc,
    integer     ld
  ) {
    m_data   = _data;
    m_nrows  = nr;
    m_ncols  = nc;
    m_ldData = ld;
    UTILS_ASSERT_DEBUG(
      nr >= 0 && nc >= 0 && m_ldData >= nr,
      "MatrixWrapper( data, nr={}, nc={}, ld={}) bad dimensions\n",
      nr, nc, ld
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::scale_by( T sc ) {
    if ( m_ldData == m_nrows ) {
      scal( m_nrows*m_ncols, sc, m_data, 1 );
    } else {
      T * p = m_data;
      for ( integer i = 0; i < m_ncols; ++i, p += m_ldData )
        scal( m_nrows, sc, p, 1 );
    }
  }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::scale_block(
    real_type sc,
    integer   nr,
    integer   nc,
    integer   irow,
    integer   icol
  ) {
    T * p = m_data+this->iaddr(irow,icol);
    for ( integer i = 0; i < nc; ++i, p += m_ldData )
      scal( nr, sc, p, 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::check( MatW const & A ) const {
    UTILS_ASSERT(
      A.nrows() == m_nrows && A.ncols() == m_ncols,
      "MatrixWrapper::check(A) size(A) = {} x {} expected {} x {}\n",
      A.nrows(), A.ncols(), m_nrows, m_ncols
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::check( Sparse const & sp ) const {
    UTILS_ASSERT(
      sp.get_number_of_rows() <= m_nrows &&
      sp.get_number_of_cols() <= m_ncols,
      "MatrixWrapper::check(sp) size(sp) = {} x {} must be contained in {} x {}\n",
      sp.get_number_of_rows(), sp.get_number_of_cols(),
      m_nrows, m_ncols
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load( real_type const data_in[], integer ldData_in ) {
    integer info = gecopy(
      m_nrows, m_ncols, data_in, ldData_in, m_data, m_ldData
    );
    UTILS_ASSERT(
      info == 0,
      "MatrixWrapper::load call lapack_wrapper::gecopy return info = {}\n", info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load_transposed( real_type const data_in[], integer ldData_in ) {
    for ( integer i = 0; i < m_nrows; ++i )
      for ( integer j = 0; j < m_ncols; ++j )
        m_data[i+j*m_ldData] = data_in[j+i*ldData_in];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load( MatW const & A ) {
    #ifndef LAPACK_WRAPPER_NO_DEBUG
    check(A);
    #endif
    integer info = gecopy(
      A.nrows(), A.ncols(), A.data(), A.ldim(),
      m_data, m_ldData
    );
    UTILS_ASSERT(
      info == 0,
      "MatrixWrapper::load call lapack_wrapper::gecopy return info = {}\n", info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load_transposed( MatW const & A ) {
    #ifndef LAPACK_WRAPPER_NO_DEBUG
    UTILS_ASSERT(
      A.ncols() == m_nrows && A.nrows() == m_ncols,
      "MatrixWrapper::load_transposed(A) size(A) = {} x {} expected {} x {}\n",
      A.nrows(), A.ncols(), m_ncols, m_nrows
    );
    #endif
    for ( integer i=0; i < m_nrows; ++i )
      for ( integer j=0; j < m_ncols; ++j )
        m_data[iaddr(i,j)] = A(j,i);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load( Sparse const & sp ) {
    #ifndef LAPACK_WRAPPER_NO_DEBUG
    check(sp);
    #endif
    integer   const * pRows;
    integer   const * pCols;
    real_type const * pValues;
    sp.get_data( pRows, pCols, pValues );
    for ( integer idx = 0; idx < sp.get_nnz(); ++idx )
      m_data[ this->iaddr(pRows[idx],pCols[idx]) ] = pValues[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load_transposed( Sparse const & sp ) {
    #ifndef LAPACK_WRAPPER_NO_DEBUG
    check(sp);
    #endif
    integer   const * pRows;
    integer   const * pCols;
    real_type const * pValues;
    sp.get_data( pCols, pRows, pValues ); // read index transposed
    for ( integer idx = 0; idx < sp.get_nnz(); ++idx )
      m_data[ this->iaddr(pRows[idx],pCols[idx]) ] = pValues[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load( Sparse const & sp, integer i_offs, integer j_offs ) {
    #ifndef LAPACK_WRAPPER_NO_DEBUG
    check(sp);
    #endif
    integer   const * pRows;
    integer   const * pCols;
    real_type const * pValues;
    sp.get_data( pRows, pCols, pValues );
    for ( integer idx = 0; idx < sp.get_nnz(); ++idx )
      m_data[ this->iaddr(pRows[idx]+i_offs,pCols[idx]+j_offs) ] = pValues[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load_transposed( Sparse const & sp, integer i_offs, integer j_offs ) {
    #ifndef LAPACK_WRAPPER_NO_DEBUG
    check(sp);
    #endif
    integer   const * pRows;
    integer   const * pCols;
    real_type const * pValues;
    sp.get_data( pCols, pRows, pValues ); // read index transposed
    for ( integer idx = 0; idx < sp.get_nnz(); ++idx )
      m_data[ this->iaddr(pRows[idx]+i_offs,pCols[idx]+j_offs) ] = pValues[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load(
    integer   const rows[],
    integer   const cols[],
    real_type const vals[],
    integer         nnz
  ) {
    for ( integer idx = 0; idx < nnz; ++idx )
      m_data[ this->iaddr(rows[idx],cols[idx]) ] = vals[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load_symmetric(
    integer   const rows[],
    integer   const cols[],
    real_type const vals[],
    integer         nnz
  ) {
    for ( integer idx = 0; idx < nnz; ++idx ) {
      integer   ii = rows[idx];
      integer   jj = cols[idx];
      real_type rr = vals[idx];
      m_data[ this->iaddr(ii,jj) ] = rr;
      if ( ii != jj ) m_data[ this->iaddr(jj,ii) ] = rr;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load(
    integer         i_offs,
    integer         j_offs,
    integer   const rows[],
    integer   const cols[],
    real_type const vals[],
    integer         nnz
  ) {
    for ( integer idx = 0; idx < nnz; ++idx )
      m_data[ this->iaddr(i_offs+rows[idx],j_offs+cols[idx]) ] = vals[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load_symmetric(
    integer         i_offs,
    integer         j_offs,
    integer   const rows[],
    integer   const cols[],
    real_type const vals[],
    integer         nnz
  ) {
    for ( integer idx = 0; idx < nnz; ++idx ) {
      integer   ii = i_offs+rows[idx];
      integer   jj = j_offs+cols[idx];
      real_type rr = vals[idx];
      m_data[ this->iaddr(ii,jj) ] = rr;
      if ( ii != jj ) m_data[ this->iaddr(jj,ii) ] = rr;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load_sparse_column(
    integer         nnz,
    real_type const values[],
    integer   const i_row[],
    integer         j
  ) {
    for ( integer k = 0; k < nnz; ++k )
      m_data[i_row[k]+j*m_ldData] = values[k];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load_sparse_row(
    integer         nnz,
    real_type const values[],
    integer         i,
    integer   const j_col[]
  ) {
    for ( integer k = 0; k < nnz; ++k )
      m_data[i+j_col[k]*m_ldData] = values[k];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::add( real_type alpha, real_type const data_in[], integer ldData_in ) {
    geadd(
      m_nrows, m_ncols,
      alpha, data_in, ldData_in,
      1.0, m_data, m_ldData,
      m_data, m_ldData
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::add( real_type const data_in[], integer ldData_in ) {
    geadd(
      m_nrows, m_ncols,
      1.0, data_in, ldData_in,
      1.0, m_data, m_ldData,
      m_data, m_ldData
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::add( Sparse const & sp ) {
    #ifndef LAPACK_WRAPPER_NO_DEBUG
    check(sp);
    #endif
    integer   const * pRows;
    integer   const * pCols;
    real_type const * pValues;
    sp.get_data( pRows, pCols, pValues );
    for ( integer idx = 0; idx < sp.get_nnz(); ++idx )
      m_data[ this->iaddr(pRows[idx],pCols[idx]) ] += pValues[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::add( real_type alpha, Sparse const & sp ) {
    #ifndef LAPACK_WRAPPER_NO_DEBUG
    check(sp);
    #endif
    integer   const * pRows;
    integer   const * pCols;
    real_type const * pValues;
    sp.get_data( pRows, pCols, pValues );
    for ( integer idx = 0; idx < sp.get_nnz(); ++idx )
      m_data[ this->iaddr(pRows[idx],pCols[idx]) ] += alpha * pValues[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::add(
    Sparse const & sp,
    integer        i_offs,
    integer        j_offs
  ) {
    integer   const * pRows;
    integer   const * pCols;
    real_type const * pValues;
    sp.get_data( pRows, pCols, pValues );
    for ( integer idx = 0; idx < sp.get_nnz(); ++idx )
      m_data[ this->iaddr(pRows[idx]+i_offs,pCols[idx]+j_offs) ] += pValues[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::add(
    real_type      alpha,
    Sparse const & sp,
    integer        i_offs,
    integer        j_offs
  ) {
    integer   const * pRows;
    integer   const * pCols;
    real_type const * pValues;
    sp.get_data( pRows, pCols, pValues );
    for ( integer idx = 0; idx < sp.get_nnz(); ++idx )
      m_data[ this->iaddr(pRows[idx]+i_offs,pCols[idx]+j_offs) ] += alpha * pValues[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*\
  :|:   __  __       _        _
  :|:  |  \/  | __ _| |_ _ __(_)_  __
  :|:  | |\/| |/ _` | __| '__| \ \/ /
  :|:  | |  | | (_| | |_| |  | |>  <
  :|:  |_|  |_|\__,_|\__|_|  |_/_/\_\
  \*/

  template <typename T>
  Matrix<T>::Matrix()
  : MatrixWrapper<T>()
  , m_mem("Matrix")
  {}

  template <typename T>
  Matrix<T>::Matrix( Matrix<T> const & rhs )
  : MatrixWrapper<T>( nullptr, rhs.m_nrows, rhs.m_ncols, rhs.m_nrows )
  , m_mem("Matrix")
  {
    size_t sz = rhs.m_nrows*rhs.m_ncols;
    m_data = m_mem.realloc( sz );
    gecopy(
      rhs.m_nrows, rhs.m_ncols, rhs.m_data, rhs.m_ldData,
      m_data, m_ldData
    );
  }

  template <typename T>
  Matrix<T>::Matrix( integer nr, integer nc )
  : MatrixWrapper<T>( nullptr, nr, nc, nr )
  , m_mem("Matrix")
  {
    m_data = m_mem.realloc( size_t(nr*nc) );
  }

  template <typename T>
  void
  Matrix<T>::setup( integer nr, integer nc ) {
    this->MatrixWrapper<T>::setup(
      m_mem.realloc( size_t(nr*nc) ), nr, nc, nr
    );
  }

  template <typename T>
  Matrix<T> const &
  Matrix<T>::operator = ( Matrix<T> const & rhs ) {
    size_t sz = rhs.m_nrows*rhs.m_ncols;
    m_data   = m_mem.realloc( sz );
    m_nrows  = rhs.m_nrows;
    m_ncols  = rhs.m_ncols;
    m_ldData = rhs.m_nrows;
    gecopy(
      rhs.m_nrows, rhs.m_ncols, rhs.m_data, rhs.m_ldData,
      m_data, m_ldData
    );
    return *this;
  }

  /*\
  :|:   ____  _             __  __       _        _
  :|:  |  _ \(_) __ _  __ _|  \/  | __ _| |_ _ __(_)_  __
  :|:  | | | | |/ _` |/ _` | |\/| |/ _` | __| '__| \ \/ /
  :|:  | |_| | | (_| | (_| | |  | | (_| | |_| |  | |>  <
  :|:  |____/|_|\__,_|\__, |_|  |_|\__,_|\__|_|  |_/_/\_\
  :|:                 |___/
  \*/

  template <typename T>
  DiagMatrix<T>::DiagMatrix( DiagMatrix<T> const & D )
  : MatrixWrapper<T>( nullptr, D.dim )
  {
    this->data = m_mem.realloc( size_t(D.dim) );
    std::copy_n( D.data, D.dim, this->data );
  }

  template <typename T>
  DiagMatrix<T>::DiagMatrix( integer dim )
  : DiagMatrixWrapper<T>( nullptr, dim )
  {
    m_data = m_mem.realloc( size_t(m_dim) );
  }

  template <typename T>
  typename DiagMatrixWrapper<T>::DMatW &
  DiagMatrixWrapper<T>::operator = ( DMatW const & COPY ) {
    UTILS_ASSERT0(
      m_dim == COPY.m_dim,
      "DiagMatrixWrapper operator = bad matrix dimensions\n"
    );
    std::copy_n( COPY.m_data, COPY.m_dim, m_data );
    return *this;
  }

  template <typename T>
  string
  DiagMatrixWrapper<T>::to_string() const {
    string res{""};
    for ( integer i{0}; i < m_dim; ++i ) {
      for ( integer j{0}; j < m_dim; ++j ) {
        if ( i != j ) res += ".              ";
        else          res += fmt::format( "{:14} ", m_data[i] );
      }
      res += '\n';
    }
    return res;
  }

  template <typename T>
  void
  DiagMatrix<T>::setup( integer dim ) {
    m_dim  = dim;
    m_data = m_mem.realloc( size_t(m_dim) );
  }

  template <typename T>
  DiagMatrix<T> const &
  DiagMatrix<T>::operator = ( DiagMatrix<T> const & rhs ) {
    m_dim  = rhs.m_dim;
    m_data = m_mem.realloc( size_t(m_dim) );
    std::copy_n( rhs.m_data, rhs.m_dim, m_data );
  }

  template class MatrixWrapper<real>;
  template class MatrixWrapper<doublereal>;

  template class Matrix<real>;
  template class Matrix<doublereal>;

}

///
/// eof: wrapper.cxx
///
