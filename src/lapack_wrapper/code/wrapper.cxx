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
    valueType * _data,
    integer     nr,
    integer     nc,
    integer     ld
  )
  : nRows(nr)
  , nCols(nc)
  , ldData(ld)
  , data(_data)
  {
    LW_ASSERT_DEBUG(
      nr >= 0 && nc >= 0 && this->ldData >= nr,
      "MatrixWrapper( data, nr={}, nc={}, ld={}) bad dimensions\n",
      nr, nc, ld
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::setup(
    valueType * _data,
    integer     nr,
    integer     nc,
    integer     ld
  ) {
    this->data   = _data;
    this->nRows  = nr;
    this->nCols  = nc;
    this->ldData = ld;
    LW_ASSERT_DEBUG(
      nr >= 0 && nc >= 0 && this->ldData >= nr,
      "MatrixWrapper( data, nr={}, nc={}, ld={}) bad dimensions\n",
      nr, nc, ld
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::scale_by( T sc ) {
    if ( this->ldData == this->nRows ) {
      scal( this->nRows*this->nCols, sc, this->data, 1 );
    } else {
      T * p = this->data;
      for ( integer i = 0; i < this->nCols; ++i, p += this->ldData )
        scal( this->nRows, sc, p, 1 );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::check( MatW const & A ) const {
    LW_ASSERT(
      A.numRows() == this->nRows && A.numCols() == this->nCols,
      "MatrixWrapper::check(A) size(A) = {} x {} expected {} x {}\n",
      A.numRows(), A.numCols(), this->nRows, this->nCols
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::check( Sparse const & sp ) const {
    LW_ASSERT(
      sp.get_number_of_rows() <= this->nRows &&
      sp.get_number_of_cols() <= this->nCols,
      "MatrixWrapper::check(sp) size(sp) = {} x {} must be contained in {} x {}\n",
      sp.get_number_of_rows(), sp.get_number_of_cols(),
      this->nRows, this->nCols
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load( valueType const data_in[], integer ldData_in ) {
    integer info = gecopy(
      this->nRows, this->nCols, data_in, ldData_in, this->data, this->ldData
    );
    LW_ASSERT(
      info == 0,
      "MatrixWrapper::load call lapack_wrapper::gecopy return info = {}\n", info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load_transposed( valueType const data_in[], integer ldData_in ) {
    for ( integer i = 0; i < this->nRows; ++i )
      for ( integer j = 0; j < this->nCols; ++j )
        this->data[i+j*this->ldData] = data_in[j+i*ldData_in];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load( MatW const & A ) {
    #ifndef LAPACK_WRAPPER_NO_DEBUG
    check(A);
    #endif
    integer info = gecopy(
      A.numRows(), A.numCols(), A.get_data(), A.lDim(),
      this->data, this->ldData
    );
    LW_ASSERT(
      info == 0,
      "MatrixWrapper::load call lapack_wrapper::gecopy return info = {}\n", info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load_transposed( MatW const & A ) {
    #ifndef LAPACK_WRAPPER_NO_DEBUG
    LW_ASSERT(
      A.numCols() == this->nRows && A.numRows() == this->nCols,
      "MatrixWrapper::load_transposed(A) size(A) = {} x {} expected {} x {}\n",
      A.numRows(), A.numCols(), this->nCols, this->nRows
    );
    #endif
    for ( integer i=0; i < this->nRows; ++i )
      for ( integer j=0; j < this->nCols; ++j )
        this->data[iaddr(i,j)] = A(j,i);
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
    valueType const * pValues;
    sp.get_data( pRows, pCols, pValues );
    for ( integer idx = 0; idx < sp.get_nnz(); ++idx )
      this->data[ this->iaddr(pRows[idx],pCols[idx]) ] = pValues[idx];
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
    valueType const * pValues;
    sp.get_data( pCols, pRows, pValues ); // read index transposed
    for ( integer idx = 0; idx < sp.get_nnz(); ++idx )
      this->data[ this->iaddr(pRows[idx],pCols[idx]) ] = pValues[idx];
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
    valueType const * pValues;
    sp.get_data( pRows, pCols, pValues );
    for ( integer idx = 0; idx < sp.get_nnz(); ++idx )
      this->data[ this->iaddr(pRows[idx]+i_offs,pCols[idx]+j_offs) ] = pValues[idx];
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
    valueType const * pValues;
    sp.get_data( pCols, pRows, pValues ); // read index transposed
    for ( integer idx = 0; idx < sp.get_nnz(); ++idx )
      this->data[ this->iaddr(pRows[idx]+i_offs,pCols[idx]+j_offs) ] = pValues[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load(
    integer   const rows[],
    integer   const cols[],
    valueType const vals[],
    integer         nnz
  ) {
    for ( integer idx = 0; idx < nnz; ++idx )
      this->data[ this->iaddr(rows[idx],cols[idx]) ] = vals[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load_symmetric(
    integer   const rows[],
    integer   const cols[],
    valueType const vals[],
    integer         nnz
  ) {
    for ( integer idx = 0; idx < nnz; ++idx ) {
      integer   ii = rows[idx];
      integer   jj = cols[idx];
      valueType rr = vals[idx];
      this->data[ this->iaddr(ii,jj) ] = rr;
      if ( ii != jj ) this->data[ this->iaddr(jj,ii) ] = rr;
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
    valueType const vals[],
    integer         nnz
  ) {
    for ( integer idx = 0; idx < nnz; ++idx )
      this->data[ this->iaddr(i_offs+rows[idx],j_offs+cols[idx]) ] = vals[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load_symmetric(
    integer         i_offs,
    integer         j_offs,
    integer   const rows[],
    integer   const cols[],
    valueType const vals[],
    integer         nnz
  ) {
    for ( integer idx = 0; idx < nnz; ++idx ) {
      integer   ii = i_offs+rows[idx];
      integer   jj = j_offs+cols[idx];
      valueType rr = vals[idx];
      this->data[ this->iaddr(ii,jj) ] = rr;
      if ( ii != jj ) this->data[ this->iaddr(jj,ii) ] = rr;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load_sparse_column(
    integer         nnz,
    valueType const values[],
    integer   const i_row[],
    integer         j
  ) {
    for ( integer k = 0; k < nnz; ++k )
      data[i_row[k]+j*ldData] = values[k];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load_sparse_row(
    integer         nnz,
    valueType const values[],
    integer         i,
    integer   const j_col[]
  ) {
    for ( integer k = 0; k < nnz; ++k )
      data[i+j_col[k]*ldData] = values[k];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::add( valueType alpha, valueType const data_in[], integer ldData_in ) {
    geadd(
      this->nRows, this->nCols,
      alpha, data_in, ldData_in,
      1.0, this->data, this->ldData,
      this->data, this->ldData
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::add( valueType const data_in[], integer ldData_in ) {
    geadd(
      this->nRows, this->nCols,
      1.0, data_in, ldData_in,
      1.0, this->data, this->ldData,
      this->data, this->ldData
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
    valueType const * pValues;
    sp.get_data( pRows, pCols, pValues );
    for ( integer idx = 0; idx < sp.get_nnz(); ++idx )
      this->data[ this->iaddr(pRows[idx],pCols[idx]) ]
        += pValues[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::add( valueType alpha, Sparse const & sp ) {
    #ifndef LAPACK_WRAPPER_NO_DEBUG
    check(sp);
    #endif
    integer   const * pRows;
    integer   const * pCols;
    valueType const * pValues;
    sp.get_data( pRows, pCols, pValues );
    for ( integer idx = 0; idx < sp.get_nnz(); ++idx )
      this->data[ this->iaddr(pRows[idx],pCols[idx]) ]
        += alpha * pValues[idx];
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
    valueType const * pValues;
    sp.get_data( pRows, pCols, pValues );
    for ( integer idx = 0; idx < sp.get_nnz(); ++idx )
      this->data[ this->iaddr(pRows[idx]+i_offs,pCols[idx]+j_offs) ]
        += pValues[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::add(
    valueType      alpha,
    Sparse const & sp,
    integer        i_offs,
    integer        j_offs
  ) {
    integer   const * pRows;
    integer   const * pCols;
    valueType const * pValues;
    sp.get_data( pRows, pCols, pValues );
    for ( integer idx = 0; idx < sp.get_nnz(); ++idx )
      this->data[ this->iaddr(pRows[idx]+i_offs,pCols[idx]+j_offs) ]
        += alpha * pValues[idx];
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
  , mem("Matrix")
  {}

  template <typename T>
  Matrix<T>::Matrix( Matrix<T> const & rhs )
  : MatrixWrapper<T>( nullptr, rhs.nRows, rhs.nCols, rhs.nRows )
  , mem("Matrix")
  {
    size_t sz = rhs.nRows*rhs.nCols;
    mem.allocate( sz );
    this->data = mem( sz );
    gecopy( rhs.nRows,  rhs.nCols, rhs.data, rhs.ldData,
            this->data, this->ldData );
  }

  template <typename T>
  Matrix<T>::Matrix( integer nr, integer nc )
  : MatrixWrapper<T>( nullptr, nr, nc, nr )
  , mem("Matrix")
  {
    mem.allocate( size_t(nr*nc) );
    this->data = mem( size_t(nr*nc) );
  }

  template <typename T>
  void
  Matrix<T>::setup( integer nr, integer nc ) {
    mem.allocate( size_t(nr*nc) );
    this->MatrixWrapper<T>::setup( mem( size_t(nr*nc) ), nr, nc, nr );
  }

  template <typename T>
  Matrix<T> const &
  Matrix<T>::operator = ( Matrix<T> const & rhs ) {
    size_t sz = rhs.nRows*rhs.nCols;
    mem.allocate( sz );
    this->data   = mem( sz );
    this->nRows  = rhs.nRows;
    this->nCols  = rhs.nCols;
    this->ldData = rhs.nRows;
    gecopy( rhs.nRows,  rhs.nCols, rhs.data, rhs.ldData,
            this->data, this->ldData );
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
  DiagMatrix<T>::DiagMatrix()
  : DiagMatrixWrapper<T>()
  , mem("DiagMatrix")
  {}

  template <typename T>
  DiagMatrix<T>::DiagMatrix( DiagMatrix<T> const & D )
  : MatrixWrapper<T>( nullptr, D.dim )
  , mem("DiagMatrix")
  {
    mem.allocate( size_t(D.dim) );
    this->data = mem( size_t(D.dim) );
    std::copy_n( D.data, D.dim, this->data );
  }

  template <typename T>
  DiagMatrix<T>::DiagMatrix( integer _dim )
  : DiagMatrixWrapper<T>( nullptr, _dim )
  , mem("DiagMatrix")
  {
    mem.allocate( size_t(this->dim) );
    this->data = mem( size_t(this->dim) );
  }

  template <typename T>
  void
  DiagMatrix<T>::setup( integer _dim ) {
    this->dim = _dim;
    mem.allocate( size_t(this->dim) );
    this->dim  = _dim;
    this->data = mem( size_t(this->dim) );
  }

  template <typename T>
  DiagMatrix<T> const &
  DiagMatrix<T>::operator = ( DiagMatrix<T> const & rhs ) {
    mem.allocate( size_t(rhs.dim) );
    this->dim  = rhs.dim;
    this->data = mem( size_t(rhs.dim) );
    std::copy_n( rhs.data, rhs.dim, this->data );
  }

  template class MatrixWrapper<real>;
  template class MatrixWrapper<doublereal>;

  template class Matrix<real>;
  template class Matrix<doublereal>;

}

///
/// eof: wrapper.cxx
///
