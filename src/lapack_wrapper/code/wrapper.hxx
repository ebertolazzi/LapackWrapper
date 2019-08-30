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
/// file: wrapper.hxx
///

#ifdef __clang__
#pragma clang diagnostic ignored "-Wpadded"
#endif

namespace lapack_wrapper {

  /*\
  :|:   ___  _           __  __      _       _    __      __
  :|:  |   \(_)__ _ __ _|  \/  |__ _| |_ _ _(_)_ _\ \    / / _ __ _ _ __ _ __  ___ _ _
  :|:  | |) | / _` / _` | |\/| / _` |  _| '_| \ \ /\ \/\/ / '_/ _` | '_ \ '_ \/ -_) '_|
  :|:  |___/|_\__,_\__, |_|  |_\__,_|\__|_| |_/_\_\ \_/\_/|_| \__,_| .__/ .__/\___|_|
  :|:              |___/                                           |_|  |_|
  \*/
  template <typename T>
  class DiagMatrixWrapper {

  public:
    typedef T                    valueType;
    typedef DiagMatrixWrapper<T> DMatW;

  private:
    integer     dim;
    valueType * data;
  public:

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: build an empy wrapper
    \*/
    explicit
    DiagMatrixWrapper( )
    : dim(0)
    , data(nullptr)
    {}

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: Map a piece of memory into a matrix object
    :|:
    :|: \param _data  pointer of the memory to be mapped as a matrix
    :|: \param _dim   number of elements on the diagonal of the mapped matrix
    \*/
    explicit
    DiagMatrixWrapper( valueType * _data, integer _dim ) {
      this->data = _data;
      this->dim  = _dim;
    }

    integer getDim()   const { return this->dim;}  //!< Number of elements
    integer numElems() const { return this->dim;}  //!< Number of elements

    valueType const * get_data() const { return this->data; }
    valueType       * get_data()       { return this->data; }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: Map a piece of memory into a matrix object
    :|:
    :|: \param _data pointer of the memory to be mapped as a matrix
    :|: \param _dim  dimension the mapped diagonal matrix
    \*/
    void
    setup( valueType * _data, integer _dim ) {
      this->data = _data;
      this->dim  = _dim;
    }

    valueType const &
    operator [] ( integer i ) const
    { return this->data[i]; }

    valueType &
    operator [] ( integer i )
    { return this->data[i]; }

    DMatW const & operator = (valueType v) {
      fill( this->dim, this->data, v );
      return *this;
    }

    DMatW & operator = ( DMatW const & COPY ) {
      LAPACK_WRAPPER_ASSERT(
        this->dim == COPY.dim,
        "DiagMatrixWrapper operator = bad matrix dimensions"
      );
      std::copy( COPY.data, COPY.data+COPY.dim, this->data );
      return *this;
    }

    void
    print( ostream_type & stream ) const {
      for ( integer i = 0; i < this->dim; ++i ) {
        for ( integer j = 0; j < this->dim; ++j ) {
          stream << std::setw(14);
          if ( i != j ) stream << '.';
          else          stream << this->data[i];
          stream << ' ';
        }
        stream << '\n';
      }
    }

  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  inline
  ostream_type &
  operator << ( ostream_type & stream, DiagMatrixWrapper<T> const & A ) {
    A.print( stream );
    return stream;
  }

  /*
  //   __  __       _        _     __        __
  //  |  \/  | __ _| |_ _ __(_)_  _\ \      / / __ __ _ _ __  _ __   ___ _ __
  //  | |\/| |/ _` | __| '__| \ \/ /\ \ /\ / / '__/ _` | '_ \| '_ \ / _ \ '__|
  //  | |  | | (_| | |_| |  | |>  <  \ V  V /| | | (_| | |_) | |_) |  __/ |
  //  |_|  |_|\__,_|\__|_|  |_/_/\_\  \_/\_/ |_|  \__,_| .__/| .__/ \___|_|
  //                                                   |_|   |_|
  */

  template <typename T>
  class MatrixWrapper {

    typedef T                    valueType;
    typedef MatrixWrapper<T>     MatW;
    typedef DiagMatrixWrapper<T> DiagW;
    typedef SparseMatrixBase<T>  Sparse;

  protected:

    integer     nRows; //!< Number of rows
    integer     nCols; //!< Number of columns
    integer     ldData;  //!< Leadind dimension
    valueType * data;

    #if defined(DEBUG) || defined(_DEBUG)
    integer
    iaddr( integer i,  integer j ) const {
      LAPACK_WRAPPER_ASSERT(
        i >= 0 && i < this->nRows && j >= 0 && j < this->nCols,
        "iaddr(" << i << ", " << j << ") out of range [0," <<
        this->nRows << ") x [0," << this->nCols << ")"
      );
      return i + j*this->ldData;
    }
    #else
    integer
    iaddr( integer i,  integer j ) const
    { return i + j*this->ldData; }
    #endif

    void
    check( Sparse const & sp ) const;

    void
    check( MatW const & M ) const;

  public:

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: build an empy wrapper
    \*/
    explicit
    MatrixWrapper( )
    : nRows(0)
    , nCols(0)
    , ldData(0)
    , data(nullptr)
    {
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: Map a piece of memory into a matrix object
    :|:
    :|: \param _data  pointer of the memory to be mapped as a matrix
    :|: \param nr     number of rows of the mapped matrix
    :|: \param nc     number of columns of the mapped matrix
    :|: \param ld     leading dimension of the matrix (Fortran addressing 0 based)
    \*/
    explicit
    MatrixWrapper(
      valueType * _data,
      integer     nr,
      integer     nc,
      integer     ld
    );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    explicit
    MatrixWrapper( MatrixWrapper<T> const & M )
    : nRows(M.nRows)
    , nCols(M.nCols)
    , ldData(M.ldData)
    , data(M.data)
    {
    }

    MatrixWrapper<T> const &
    operator = ( MatrixWrapper<T> const & rhs ) {
      this->nRows  = rhs.nRows;
      this->nCols  = rhs.nCols;
      this->ldData = rhs.ldData;
      this->data   = rhs.data;
      return *this;
    }

    integer numRows()  const { return this->nRows; }  //!< Number of rows
    integer numCols()  const { return this->nCols; }  //!< Number of columns
    integer lDim()     const { return this->ldData; } //!< Leading dimension
    integer numElems() const { return this->nRows*this->nCols; }  //!< Number of elements

    valueType const * get_data() const { return this->data; }
    valueType       * get_data()       { return this->data; }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: Map a piece of memory into a matrix object
    :|:
    :|: \param _data  pointer of the memory to be mapped as a matrix
    :|: \param nr     number of rows of the mapped matrix
    :|: \param nc     number of columns of the mapped matrix
    :|: \param ld     leading dimension of the matrix (Fortran addressing 0 based)
    \*/
    void
    setup(
      valueType * _data,
      integer     nr,
      integer     nc,
      integer     ld
    );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: Access element (i,j) of the matrix
    :|:
    :|: \param i row of the element
    :|: \param j column of the element
    \*/
    valueType const &
    operator () ( integer i,  integer j ) const
    { return this->data[this->iaddr(i,j)]; }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: Access element (i,j) of the matrix
    :|:
    :|: \param i row of the element
    :|: \param j column of the element
    \*/
    valueType &
    operator () ( integer i,  integer j )
    { return this->data[this->iaddr(i,j)]; }

    /*!
    :|: Fill the matrix with zeros
    \*/
    void
    zero_fill() {
      gezero( this->nRows, this->nCols, this->data, this->ldData );
    }

    /*!
    :|: Fill the matrix with zeros
    \*/
    void
    setZero() { zero_fill(); }

    /*!
    :|: Fill the matrix with value `val`
    \*/
    void
    fill( valueType val ) {
      gefill( this->nRows, this->nCols, this->data, this->ldData, val );
    }

    void
    scale_by( valueType sc );

    /*!
    :|: Zeroes a rectangular block of the stored matrix
    :|: staring at `(irow,icol)` position
    :|:
    :|: \param[in] nr    number of rows of the block to be zeroed
    :|: \param[in] nc    number of columns of the block to be zeroed
    :|: \param[in] irow  starting row
    :|: \param[in] icol  stating column
    \*/
    void
    zero_block(
      integer nr,
      integer nc,
      integer irow,
      integer icol
    ) {
      gezero( nr, nc, this->data + this->iaddr(irow,icol), this->ldData );
    }

    /*!
    :|: Initialize the matrix as an identity matrix
    \*/
    void
    id( valueType dg ) {
      geid( this->nRows, this->nCols, this->data, this->ldData, dg );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: Initialize matrix
    :|:
    :|: \param data   pointer of memory with data to be copied
    :|: \param ldData leading dimension of the memory to be copied
    :|:
    \*/
    void load( valueType const data[], integer ldData );

    /*!
    :|: Initialize matrix
    :|:
    :|: \param A initialize matrix with the matrix of the object A
    :|:
    \*/
    void load( MatW const & A );

    /*!
    :|:  Copy a matrix to a rectangular block of the stored matrix
    :|:  staring at `(irow,icol)` position
    :|:
    :|:  \param[in] B     reference to the matrix of the object B
    :|:  \param[in] irow  starting row
    :|:  \param[in] icol  stating column
    \*/
    void
    load_block(
      MatW const & B,
      integer      irow = 0,
      integer      icol = 0
    ) {
      this->load_block(
        B.nRows, B.nCols, B.data, B.ldData, irow, icol
      );
    }

    /*!
    :|:  Copy a matrix to a rectangular block of the stored matrix
    :|:  staring at `(irow,icol)` position
    :|:
    :|:  \param[in] B     reference to the matrix of the object B
    :|:  \param[in] irow  starting row
    :|:  \param[in] icol  stating column
    \*/
    void
    load_block_transposed(
      MatW const & B,
      integer      irow = 0,
      integer      icol = 0
    ) {
      this->load_block_transposed(
        B.nRows, B.nCols, B.data, B.ldData, irow, icol
      );
    }

    /*!
    :|:  Copy a matrix to a rectangular block of the stored matrix
    :|:  staring at `(irow,icol)` position
    :|:
    :|:  \param[in] nr    number of rows of the block to be zeroed
    :|:  \param[in] nc    number of columns of the block to be zeroed
    :|:  \param[in] B     pointer to memory storing the input matrix `B`
    :|:  \param[in] ldB   leading dimension of the matrix `B`
    :|:  \param[in] irow  starting row
    :|:  \param[in] icol  stating column
    \*/
    void
    load_block(
      integer         nr,
      integer         nc,
      valueType const B[],
      integer         ldB,
      integer         irow = 0,
      integer         icol = 0
    ) {
      LAPACK_WRAPPER_ASSERT(
        irow + nr <= this->nRows &&
        icol + nc <= this->nCols &&
        irow >= 0 && icol >= 0,
        "load_block( nr = " << nr << " nc = " << nc <<
        ",..., irow = " << irow << ", icol = " << icol <<
        ") bad parameters"
      );
      integer info = gecopy(
        nr, nc,
        B, ldB,
        this->data + this->iaddr(irow,icol), this->ldData
      );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "load_block call lapack_wrapper::gecopy return info = " << info
      );
    }

    /*!
    :|:  Copy a matrix to a rectangular block of the stored matrix
    :|:  staring at `(irow,icol)` position
    :|:
    :|:  \param[in] nr    number of rows of the block to be zeroed
    :|:  \param[in] nc    number of columns of the block to be zeroed
    :|:  \param[in] B     pointer to memory storing the input matrix `B`
    :|:  \param[in] ldB   leading dimension of the matrix `B`
    :|:  \param[in] irow  starting row
    :|:  \param[in] icol  stating column
    \*/
    void
    load_block_transposed(
      integer         nr,
      integer         nc,
      valueType const B[],
      integer         ldB,
      integer         irow = 0,
      integer         icol = 0
    ) {
      LAPACK_WRAPPER_ASSERT(
        irow + nc <= this->nRows &&
        icol + nr <= this->nCols &&
        irow >= 0 && icol >= 0,
        "load_block_transpose( nr = " << nr << " nc = " << nc <<
        ",..., irow = " << irow << ", icol = " << icol <<
        ") bad parameters"
      );
      valueType const * pd = B;
      valueType       * pp = this->data + this->iaddr(irow,icol);
      for ( integer i = 0; i < nc; ++i, pd += ldB, ++pp )
        lapack_wrapper::copy( nr, pd, 1, pp, this->ldData );
    }

    /*!
    :|:  Copy a matrix to a rectangular block of the stored matrix
    :|:  staring at `(irow,icol)` position
    :|:
    :|:  \param[in] n     number of element of the diagonal
    :|:  \param[in] D     pointer to memory storing the diagonal
    :|:  \param[in] irow  starting row
    :|:  \param[in] icol  stating column
    \*/
    void
    load_diagonal_block(
      integer         n,
      valueType const D[],
      integer         irow = 0,
      integer         icol = 0
    ) {
      LAPACK_WRAPPER_ASSERT(
        irow + n <= this->nRows &&
        icol + n <= this->nCols &&
        irow >= 0 && icol >= 0,
        "load_diagonal_block( n = " << n <<
        ",..., irow = " << irow << ", icol = " << icol <<
        ") bad parameters"
      );
      for ( integer i = 0; i < n; ++i )
        this->data[this->iaddr(irow+i,icol+i)] = D[i];
    }

    /*!
    :|:  Copy a matrix to a rectangular block of the stored matrix
    :|:  staring at `(irow,icol)` position
    :|:
    :|:  \param[in] D     reference to the matrix of the object D
    :|:  \param[in] irow  starting row
    :|:  \param[in] icol  stating column
    \*/
    void
    load_block(
      DiagW const & D,
      integer       irow = 0,
      integer       icol = 0
    ) {
      this->load_diagonal_block( D.getDim(), D.get_data(), irow, icol );
    }

    /*!
    :|:  Copy vector `column` to the `icol`th column of the internal stored matrix
    :|:  \param[in] column the column vector
    :|:  \param[in] icol   the column to be changed
    \*/
    void
    load_column( valueType const column[], integer icol )
    { copy( this->nRows, column, 1, this->data + icol * this->ldData, 1 ); }

    /*!
    :|:  Copy vector `row` to the `irow`th row of the internal stored matrix
    :|:  \param[in] row  the row vector
    :|:  \param[in] irow the row to be changed
    \*/
    void
    load_row( valueType const row[], integer irow )
    { copy( this->nCols, row, 1, this->data + irow, this->ldData ); }

    /*!
    :|: Initialize matrix to 0 and next copy sparse matrix
    :|:
    :|: \param sp sparse matrix to be copied
    :|:
    \*/
    void load0( Sparse const & sp );

    void
    load0(
      integer   const rows[],
      integer   const cols[],
      valueType const vals[],
      integer         nnz
    );

    /*!
    :|: Copy sparse matrix into the object
    :|:
    :|: \param sp sparse matrix to be copied
    :|:
    \*/
    void load( Sparse const & sp );
    void load_transposed( Sparse const & sp );

    void
    load(
      integer   const rows[],
      integer   const cols[],
      valueType const vals[],
      integer         nnz
    );

    /*!
    :|: Copy sparse matrix into the object
    :|:
    :|: \param sp     sparse matrix to be copied
    :|: \param i_offs offset added to to the row indices
    :|: \param j_offs offset added to to the column indices
    :|:
    \*/
    void load( Sparse const & sp, integer i_offs, integer j_offs );
    void load_transposed( Sparse const & sp, integer i_offs, integer j_offs );

    void
    load(
      integer         i_offs,
      integer         j_offs,
      integer   const rows[],
      integer   const cols[],
      valueType const vals[],
      integer         nnz
    );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: Add a block of memory to the matrix scaled by alpha
    :|:
    :|: \param alpha  scaling factor
    :|: \param data   pointer of memory with data to be copied
    :|: \param ldData leading dimension of the memory to be copied
    :|:
    \*/
    void add( valueType alpha, valueType const data[], integer ldData );

    /*!
    :|: Add a block of memory to the matrix
    :|:
    :|: \param data   pointer of memory with data to be copied
    :|: \param ldData leading dimension of the memory to be copied
    :|:
    \*/
    void add( valueType const data[], integer ldData );

    // R <- R + alpha * A
    /*!
    :|:
    :|: \param A first mapped matrix
    :|:
    \*/
    void
    add( valueType alpha, MatW const & A ) {
      #if defined(DEBUG) || defined(_DEBUG)
      LAPACK_WRAPPER_ASSERT(
        nRows == A.nRows && nCols == A.nCols,
        "add(..) incompatible dimensions"
      );
      #endif
      this->add( alpha, A.data, A.ldData );
    }

    void
    add( MatW const & A ) {
      #if defined(DEBUG) || defined(_DEBUG)
      LAPACK_WRAPPER_ASSERT(
        nRows == A.nRows && nCols == A.nCols,
        "add(..) incompatible dimensions"
      );
      #endif
      this->add( A.data, A.ldData );
    }

    /*!
    :|: Add sparse matrix to the object
    :|:
    :|: \param sp sparse matrix to be copied
    :|:
    \*/
    void add( Sparse const & sp );

    /*!
    :|: Add sparse matrix to the object
    :|:
    :|: \param sp     sparse matrix to be copied
    :|: \param i_offs offset added to to the row indices
    :|: \param j_offs offset added to to the column indices
    :|:
    \*/
    void add( Sparse const & sp, integer i_offs, integer j_offs );

    /*!
    :|: Add sparse multiplied by `alpha` matrix to the object
    :|:
    :|: \param alpha  scalar used to multiply the sparse matrix
    :|: \param sp     sparse matrix to be copied
    :|:
    \*/
    void add( valueType alpha, Sparse const & sp );

    /*!
    :|: Add sparse multiplied by `alpha` matrix to the object
    :|:
    :|: \param alpha  scalar used to multiply the sparse matrix
    :|: \param sp     sparse matrix to be copied
    :|: \param i_offs offset added to to the row indices
    :|: \param j_offs offset added to to the column indices
    :|:
    \*/
    void
    add(
      valueType      alpha,
      Sparse const & sp,
      integer        i_offs,
      integer        j_offs
    );

    // alpha*A + beta*B -> C
    /*!
    :|: Add a linear combination of two matrix and assign to a third one
    :|:
    :|: \param alpha scalar used to multiply the matrix `A`
    :|: \param A     first mapped matrix
    :|: \param beta  scalar used to multiply the matrix `B`
    :|: \param B     second mapped matrix
    :|: \param C     mapped matrix which store the result
    :|:
    \*/
    friend
    void
    add(
      valueType    alpha,
      MatW const & A,
      valueType    beta,
      MatW const & B,
      MatW       & C
    ) {
      geadd(
        C.nRows, C.nCols,
        alpha, A.data, A.ldData,
        beta,  B.data, B.ldData,
        C.data, C.ldData
      );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: make a view of a block
    :|:
    :|: \param i_offs initial row of the block to be extracted
    :|: \param j_offs initial column of the block to be extracted
    :|: \param nrow   number of rows of the block to be extracted
    :|: \param ncol   number of columns of the block to be extracted
    :|: \param to     mapped matrix which store the result
    :|:
    \*/
    void
    view_block(
      integer i_offs,
      integer j_offs,
      integer nrow,
      integer ncol,
      MatW &  to
    ) {
      #if defined(DEBUG) || defined(_DEBUG)
      LAPACK_WRAPPER_ASSERT(
        i_offs >= 0 && i_offs+nrow <= this->nRows &&
        j_offs >= 0 && j_offs+ncol <= this->nCols,
        "view_block(i_offs=" << i_offs << ", j_offs=" << j_offs <<
        ", nrow=" << nrow << ", ncol=" << ncol <<
        ") will be out of range [0," <<
        this->nRows << ") x [0," << this->nCols << ")"
      );
      #endif
      to.setup(
        this->data+this->iaddr(i_offs,j_offs),
        nrow, ncol, this->ldData
      );
    }

    /*!
    :|: Extract a matrix in transposed form
    :|:
    :|: \param out mapped matrix which store the result
    :|:
    \*/
    void
    get_transposed( MatW & out ) {
      #if defined(DEBUG) || defined(_DEBUG)
      LAPACK_WRAPPER_ASSERT(
        this->nRows == out.nCols && this->nCols == out.nRows,
        "get_transposed(...) incompatible matrices"
      );
      #endif
      valueType const * pc = this->data;
      for ( integer i = 0; i < this->nCols; ++i, pc += this->ldData )
        lapack_wrapper::copy( this->nRows, pc, 1, out.data + i, out.ldData );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: Extract a block
    :|:
    :|: \param i_offs initial row of the block to be extracted
    :|: \param j_offs initial column of the block to be extracted
    :|: \param to     mapped matrix which store the result
    :|:
    \*/
    void
    get_block(
      MatW &  to,
      integer i_offs,
      integer j_offs
    ) {
      #if defined(DEBUG) || defined(_DEBUG)
      LAPACK_WRAPPER_ASSERT(
        i_offs >= 0 && i_offs+to.nRows <= this->nRows &&
        j_offs >= 0 && j_offs+to.nCols <= this->nCols,
        "get_block(i_offs=" << i_offs << ", j_offs=" << j_offs <<
        ", ..) will be out of range [0," <<
        this->nRows << ") x [0," << this->nCols << ")"
      );
      #endif
      lapack_wrapper::gecopy(
        to.nRows, to.nCols,
        this->data+this->iaddr(i_offs,j_offs), this->ldData,
        to.data, to.ldData
      );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: Extract a block
    :|:
    :|: \param i_offs initial row of the block to be extracted
    :|: \param j_offs initial column of the block to be extracted
    :|: \param to     mapped matrix which store the result
    :|:
    \*/
    void
    get_block_transposed(
      MatW &  to,
      integer i_offs,
      integer j_offs
    ) {
      #if defined(DEBUG) || defined(_DEBUG)
      LAPACK_WRAPPER_ASSERT(
        i_offs >= 0 && i_offs+to.nCols <= this->nRows &&
        j_offs >= 0 && j_offs+to.nRows <= this->nCols,
        "get_block_transposed(i_offs=" << i_offs << ", j_offs=" << j_offs <<
        ", ..) will be out of range [0," <<
        this->nRows << ") x [0," << this->nCols << ")"
      );
      #endif
      valueType const * pc = this->data+this->iaddr(i_offs,j_offs);
      for ( integer i = 0; i < to.nCols; ++i, pc += this->ldData )
        lapack_wrapper::copy( to.nRows, pc, 1, to.data + i, to.ldData );
    }

    void
    print( ostream_type & stream ) const {
      for ( integer i = 0; i < this->nRows; ++i ) {
        for ( integer j = 0; j < this->nCols; ++j )
          stream << std::setw(14) << (*this)(i,j) << ' ';
        stream << '\n';
      }
    }

    void
    print0( ostream_type & stream, valueType eps ) const {
      for ( integer i = 0; i < this->nRows; ++i ) {
        for ( integer j = 0; j < this->nCols; ++j ) {
          valueType aij = (*this)(i,j);
          stream << std::setw(14);
          if ( std::abs( aij ) < eps ) stream << '.';
          else                         stream << aij;
          stream << ' ';
        }
        stream << '\n';
      }
    }
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  template <typename T>
  inline
  ostream_type &
  operator << ( ostream_type & stream, MatrixWrapper<T> const & A ) {
    A.print(stream);
    return stream;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // c = beta*c + alpha*A*v
  /*!
  :|: Perform matrix vector multiplication `c = beta*c + alpha*A*v`
  :|:
  :|: \param alpha  matrix `A` is multiplied by `alpha`
  :|: \param TRANSA choose if multiply with the transpose of `A`
  :|: \param A      matrix used in the multiplication
  :|: \param v      vector to be multiplied
  :|: \param incv   stride of the vector `v`
  :|: \param beta   scalar used to multiply `c`
  :|: \param c      result vector
  :|: \param incc   stride of the vector `c`
  :|:
  \*/
  template <typename T>
  inline
  void
  gemv(
    T                        alpha,
    Transposition    const & TRANSA,
    MatrixWrapper<T> const & A,
    T const                  v[],
    integer                  incv,
    T                        beta,
    T                        c[],
    integer                  incc
  ) {
    lapack_wrapper::gemv(
      TRANSA,
      A.numRows(), A.numCols(),
      alpha,
      A.get_data(), A.lDim(),
      v, incv,
      beta,
      c, incc
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // c = beta*c + alpha*A*v
  /*!
  :|: Perform matrix vector multiplication `c = beta*c + alpha*A*v`
  :|:
  :|: \param alpha  matrix `A` is multiplied by `alpha`
  :|: \param A      matrix used in the multiplication
  :|: \param v      vector to be multiplied
  :|: \param incv   stride of the vector `v`
  :|: \param beta   scalar used to multiply `c`
  :|: \param c      result vector
  :|: \param incc   stride of the vector `c`
  :|:
  \*/
  template <typename T>
  inline
  void
  gemv(
    T                        alpha,
    MatrixWrapper<T> const & A,
    T const                  v[],
    integer                  incv,
    T                        beta,
    T                        c[],
    integer                  incc
  ) {
    lapack_wrapper::gemv(
      NO_TRANSPOSE,
      A.numRows(), A.numCols(),
      alpha,
      A.get_data(), A.lDim(),
      v, incv,
      beta,
      c, incc
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*!
  :|: Perform matrix matrix multiplication `C = beta*C + alpha*A*B`
  :|:
  :|: \param where  message added in case of error
  :|: \param alpha  matrix `A` is multiplied by `alpha`
  :|: \param A      matrix used in the multiplication
  :|: \param B      matrix used in the multiplication
  :|: \param beta   scalar used to multiply `C`
  :|: \param C      result matrix
  :|:
  \*/
  // C = beta*C + alpha*A*B
  template <typename T>
  inline
  void
  gemm(
    char const               where[],
    T                        alpha,
    MatrixWrapper<T> const & A,
    MatrixWrapper<T> const & B,
    T                        beta,
    MatrixWrapper<T>       & C
  ) {
    LAPACK_WRAPPER_ASSERT(
      A.numCols() == B.numRows() &&
      A.numRows() == C.numRows() &&
      B.numCols() == C.numCols(),
      "gemm, at `" << where << "' inconsistent dimensions: " <<
      "\nA = " << A.numRows() << " x " << A.numCols() <<
      "\nB = " << B.numRows() << " x " << B.numCols() <<
      "\nC = " << C.numRows() << " x " << C.numCols()
    );
    lapack_wrapper::gemm(
      NO_TRANSPOSE,
      NO_TRANSPOSE,
      A.numRows(), B.numCols(), A.numCols(),
      alpha,
      A.get_data(), A.lDim(),
      B.get_data(), B.lDim(),
      beta,
      C.get_data(), C.lDim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*!
  :|: Perform matrix matrix multiplication `C = beta*C + alpha*A*B`
  :|:
  :|: \param alpha  matrix `A` is multiplied by `alpha`
  :|: \param A      matrix used in the multiplication
  :|: \param B      matrix used in the multiplication
  :|: \param beta   scalar used to multiply `C`
  :|: \param C      result matrix
  :|:
  \*/
 // C = beta*C + alpha*A*B (NO DEBUG)
  template <typename T>
  inline
  void
  gemm(
    T                        alpha,
    MatrixWrapper<T> const & A,
    MatrixWrapper<T> const & B,
    T                        beta,
    MatrixWrapper<T>       & C
  ) {
    LAPACK_WRAPPER_ASSERT(
      A.numCols() == B.numRows() &&
      A.numRows() == C.numRows() &&
      B.numCols() == C.numCols(),
      "gemm, inconsistent dimensions: " <<
      "\nA = " << A.numRows() << " x " << A.numCols() <<
      "\nB = " << B.numRows() << " x " << B.numCols() <<
      "\nC = " << C.numRows() << " x " << C.numCols()
    );
    lapack_wrapper::gemm(
      NO_TRANSPOSE,
      NO_TRANSPOSE,
      A.numRows(), B.numCols(), A.numCols(),
      alpha,
      A.get_data(), A.lDim(),
      B.get_data(), B.lDim(),
      beta,
      C.get_data(), C.lDim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*!
  :|: Perform matrix matrix multiplication `C = beta*C + alpha*A*B`
  :|:
  :|: \param where  message added in case of error
  :|: \param alpha  matrix `A` is multiplied by `alpha`
  :|: \param TRANSA choose `A` transposed or not
  :|: \param A      matrix used in the multiplication
  :|: \param TRANSB choose `B` transposed or not
  :|: \param B      matrix used in the multiplication
  :|: \param beta   scalar used to multiply `C`
  :|: \param C      result matrix
  :|:
  \*/
  // C = beta*C + alpha*A*B
  template <typename T>
  inline
  void
  gemm(
    char const               where[],
    T                        alpha,
    Transposition    const & TRANSA,
    MatrixWrapper<T> const & A,
    Transposition    const & TRANSB,
    MatrixWrapper<T> const & B,
    T                        beta,
    MatrixWrapper<T>       & C
  ) {
    integer Ar = TRANSA == NO_TRANSPOSE ? A.numRows() : A.numCols();
    integer Ac = TRANSA == NO_TRANSPOSE ? A.numCols() : A.numRows();
    integer Br = TRANSB == NO_TRANSPOSE ? B.numRows() : B.numCols();
    integer Bc = TRANSB == NO_TRANSPOSE ? B.numCols() : B.numRows();
    LAPACK_WRAPPER_ASSERT(
      C.numRows() == Ar && C.numCols() == Bc && Ac == Br,
      "gemm, at `" << where << "' inconsistent dimensions: " <<
      "\nA = " << A.numRows() << " x " << A.numCols() <<
      "\nB = " << B.numRows() << " x " << B.numCols() <<
      "\nC = " << C.numRows() << " x " << C.numCols() <<
      "\nA " << (NO_TRANSPOSE?"NO":"") << " transposed" <<
      "\nB " << (NO_TRANSPOSE?"NO":"") << " transposed"
    );
    lapack_wrapper::gemm(
      TRANSA,
      TRANSB,
      C.numRows(), C.numCols(), Ac,
      alpha,
      A.get_data(), A.lDim(),
      B.get_data(), B.lDim(),
      beta,
      C.get_data(), C.lDim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // C = beta*C + alpha*A*B (NO DEBUG)
  /*!
  :|: Perform matrix matrix multiplication `C = beta*C + alpha*A*B`
  :|:
  :|: \param alpha  matrix `A` is multiplied by `alpha`
  :|: \param TRANSA choose `A` transposed or not
  :|: \param A      matrix used in the multiplication
  :|: \param TRANSB choose `B` transposed or not
  :|: \param B      matrix used in the multiplication
  :|: \param beta   scalar used to multiply `C`
  :|: \param C      result matrix
  :|:
  \*/
  template <typename T>
  inline
  void
  gemm(
    T                        alpha,
    Transposition    const & TRANSA,
    MatrixWrapper<T> const & A,
    Transposition    const & TRANSB,
    MatrixWrapper<T> const & B,
    T                        beta,
    MatrixWrapper<T>       & C
  ) {
    integer Ar = TRANSA == NO_TRANSPOSE ? A.numRows() : A.numCols();
    integer Ac = TRANSA == NO_TRANSPOSE ? A.numCols() : A.numRows();
    integer Br = TRANSB == NO_TRANSPOSE ? B.numRows() : B.numCols();
    integer Bc = TRANSB == NO_TRANSPOSE ? B.numCols() : B.numRows();
    LAPACK_WRAPPER_ASSERT(
      C.numRows() == Ar && C.numCols() == Bc && Ac == Br,
      "gemm, inconsistent dimensions: " <<
      "\nA = " << A.numRows() << " x " << A.numCols() <<
      "\nB = " << B.numRows() << " x " << B.numCols() <<
      "\nC = " << C.numRows() << " x " << C.numCols() <<
      "\nA " << (NO_TRANSPOSE?"NO":"") << " transposed" <<
      "\nB " << (NO_TRANSPOSE?"NO":"") << " transposed"
    );
    lapack_wrapper::gemm(
      TRANSA,
      TRANSB,
      C.numRows(), C.numCols(), Ac,
      alpha,
      A.get_data(), A.lDim(),
      B.get_data(), B.lDim(),
      beta,
      C.get_data(), C.lDim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // x = A * x or x = A^T * x
  /*!
  :|: Perform matrix vector multiplication `x = A * x` or `x = A^T * x`
  :|:
  :|: \param UPLO   choose upper or lower part
  :|: \param TRANS  choose `A` transposed or not
  :|: \param DIAG   use or not diagonal elements to 1
  :|: \param x      vector used in the multiplcation and result
  :|: \param incx   stride of vector `x`
  :|:
  \*/
  template <typename T>
  inline
  void
  trmv(
    ULselect         const & UPLO,
    Transposition    const & TRANS,
    DiagonalType     const & DIAG,
    MatrixWrapper<T> const & A,
    T                        x[],
    integer                  incx
  ) {
    LAPACK_WRAPPER_ASSERT(
      A.numRows() == A.numCols(),
      "trmv, matrix is " << A.numRows() << " x " <<
      A.numRows() << " expected square"
    );
    lapack_wrapper::trmv(
      UPLO, TRANS, DIAG, A.numRows(), A.get_data(), A.lDim(), x, incx
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // x = A * x or x = A^T * x
  /*!
  :|: Perform matrix vector multiplication `x = A * x` or `x = A^T * x`
  :|:
  :|: \param where  message added in case of error
  :|: \param UPLO   choose upper or lower part
  :|: \param TRANS  choose `A` transposed or not
  :|: \param DIAG   use or not diagonal elements to 1
  :|: \param A      matrix to multiplied
  :|: \param x      vector used in the multiplcation and result
  :|: \param incx   stride of vector `x`
  :|:
  \*/
  template <typename T>
  inline
  void
  trmv(
    char const               where[],
    ULselect         const & UPLO,
    Transposition    const & TRANS,
    DiagonalType     const & DIAG,
    MatrixWrapper<T> const & A,
    T                        x[],
    integer                  incx
  ) {
    LAPACK_WRAPPER_ASSERT(
      A.numRows() == A.numCols(),
      "trmv at `" << where << "` matrix is " <<
      A.numRows() << " x " << A.numRows() << " expected square"
    );
    lapack_wrapper::trmv(
      UPLO, TRANS, DIAG, A.numRows(), A.get_data(), A.lDim(), x, incx
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // x = A^(-1) * x or x = A^(-T) * x
  /*!
  :|: Perform matrix vector multiplication `x = A * x` or `x = A^T * x`
  :|:
  :|: \param UPLO   choose upper or lower part
  :|: \param TRANS  choose `A` transposed or not
  :|: \param DIAG   use or not diagonal elements to 1
  :|: \param A      matrix to multiplied
  :|: \param x      vector used in the multiplcation and result
  :|: \param incx   stride of vector `x`
  :|:
  \*/
  template <typename T>
  inline
  void
  trsv(
    ULselect         const & UPLO,
    Transposition    const & TRANS,
    DiagonalType     const & DIAG,
    MatrixWrapper<T> const & A,
    T                        x[],
    integer                  incx
  ) {
    LAPACK_WRAPPER_ASSERT(
      A.numRows() == A.numCols(),
      "trsv, matrix is " << A.numRows() << " x " <<
      A.numRows() << " expected square"
    );
    lapack_wrapper::trsv(
      UPLO, TRANS, DIAG, A.numRows(), A.get_data(), A.lDim(), x, incx
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // x = A^(-1) * x or x = A^(-T) * x
  /*!
  :|: Perform matrix vector multiplication `x = A^(-1) * x` or `x = A^(-T) * x`
  :|:
  :|: \param where  message added in case of error
  :|: \param UPLO   choose upper or lower part
  :|: \param TRANS  choose `A` transposed or not
  :|: \param DIAG   use or not diagonal elements to 1
  :|: \param A      matrix to multiplied
  :|: \param x      vector used in the multiplcation and result
  :|: \param incx   stride of vector `x`
  :|:
  \*/
  template <typename T>
  inline
  void
  trsv(
    char const               where[],
    ULselect         const & UPLO,
    Transposition    const & TRANS,
    DiagonalType     const & DIAG,
    MatrixWrapper<T> const & A,
    T                        x[],
    integer                  incx
  ) {
    LAPACK_WRAPPER_ASSERT(
      A.numRows() == A.numCols(),
      "trsv at `" << where << "` matrix is " <<
      A.numRows() << " x " << A.numRows() << " expected square"
    );
    lapack_wrapper::trsv(
      UPLO, TRANS, DIAG, A.numRows(), A.get_data(), A.lDim(), x, incx
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*!
  :|: Perform matrix matrix multiplication
  :|:
  :|:    B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
  :|:    op( A ) = A   or   op( A ) = A**T.
  :|:
  :|: \param SIDE   multiply on left or right
  :|: \param UPLO   choose upper or lower part
  :|: \param TRANS  choose `A` transposed or not
  :|: \param DIAG   use or not diagonal elements to 1
  :|: \param alpha  scalar used in the multiplication
  :|: \param A      matrix to be multiplied
  :|: \param B      matrix to be multiplied
  :|:
  \*/
  template <typename T>
  inline
  void
  trmm(
    SideMultiply     const & SIDE,
    ULselect         const & UPLO,
    Transposition    const & TRANS,
    DiagonalType     const & DIAG,
    T                        alpha,
    MatrixWrapper<T> const & A,
    MatrixWrapper<T>       & B
  ) {
    LAPACK_WRAPPER_ASSERT(
      A.numRows() == A.numCols(),
      "trmm, matrix is " << A.numRows() << " x " <<
      A.numRows() << " expected square"
    );
    lapack_wrapper::trmm(
      SIDE, UPLO, TRANS, DIAG,
      B.numRows(), B.numCols(), alpha,
      A.get_data(), A.lDim(),
      B.get_data(), B.lDim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*!
  :|: Perform matrix matrix multiplication
  :|:
  :|:    B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
  :|:    op( A ) = A   or   op( A ) = A**T.
  :|:
  :|: \param where  message added in case of error
  :|: \param SIDE   multiply on left or right
  :|: \param UPLO   choose upper or lower part
  :|: \param TRANS  choose `A` transposed or not
  :|: \param DIAG   use or not diagonal elements to 1
  :|: \param alpha  scalar used in the multiplication
  :|: \param A      matrix to be multiplied
  :|: \param B      matrix to be multiplied
  :|:
  \*/
  template <typename T>
  inline
  void
  trmm(
    char const               where[],
    SideMultiply     const & SIDE,
    ULselect         const & UPLO,
    Transposition    const & TRANS,
    DiagonalType     const & DIAG,
    T                        alpha,
    MatrixWrapper<T> const & A,
    MatrixWrapper<T>       & B
  ) {
    LAPACK_WRAPPER_ASSERT(
      A.numRows() == A.numCols(),
      "trmm, at `" << where <<
      "` matrix is " << A.numRows() << " x " <<
      A.numRows() << " expected square"
    );
    lapack_wrapper::trmm(
      SIDE, UPLO, TRANS, DIAG,
      B.numRows(), B.numCols(), alpha,
      A.get_data(), A.lDim(),
      B.get_data(), B.lDim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*!
  :|: Perform matrix matrix multiplication
  :|:
  :|:    B := alpha*op( A^(-1) )*B,   or   B := alpha*B*op( A^(-1) ),
  :|:    op( A ) = A   or   op( A ) = A**T.
  :|:
  :|: \param SIDE   multiply on left or right
  :|: \param UPLO   choose upper or lower part
  :|: \param TRANS  choose `A` transposed or not
  :|: \param DIAG   use or not diagonal elements to 1
  :|: \param alpha  scalar used in the multiplication
  :|: \param A      matrix to be multiplied
  :|: \param B      matrix to be multiplied
  :|:
  \*/
  template <typename T>
  inline
  void
  trsm(
    SideMultiply     const & SIDE,
    ULselect         const & UPLO,
    Transposition    const & TRANS,
    DiagonalType     const & DIAG,
    T                        alpha,
    MatrixWrapper<T> const & A,
    MatrixWrapper<T>       & B
  ) {
    LAPACK_WRAPPER_ASSERT(
      A.numRows() == A.numCols(),
      "trsm, matrix is " << A.numRows() << " x " <<
      A.numRows() << " expected square"
    );
    lapack_wrapper::trsm(
      SIDE, UPLO, TRANS, DIAG,
      B.numRows(), B.numCols(), alpha,
      A.get_data(), A.lDim(),
      B.get_data(), B.lDim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*!
  :|: Perform matrix matrix multiplication
  :|:
  :|:    B := alpha*op( A^(-1) )*B,   or   B := alpha*B*op( A^(-1) ),
  :|:    op( A ) = A   or   op( A ) = A**T.
  :|:
  :|: \param where  message added in case of error
  :|: \param SIDE   multiply on left or right
  :|: \param UPLO   choose upper or lower part
  :|: \param TRANS  choose `A` transposed or not
  :|: \param DIAG   use or not diagonal elements to 1
  :|: \param alpha  scalar used in the multiplication
  :|: \param A      matrix to be multiplied
  :|: \param B      matrix to be multiplied
  :|:
  \*/
  template <typename T>
  inline
  void
  trsm(
    char const               where[],
    SideMultiply     const & SIDE,
    ULselect         const & UPLO,
    Transposition    const & TRANS,
    DiagonalType     const & DIAG,
    T                        alpha,
    MatrixWrapper<T> const & A,
    MatrixWrapper<T>       & B
  ) {
    LAPACK_WRAPPER_ASSERT(
      A.numRows() == A.numCols(),
      "trmm, at `" << where <<
      "` matrix is " << A.numRows() << " x " <<
      A.numRows() << " expected square"
    );
    lapack_wrapper::trsm(
      SIDE, UPLO, TRANS, DIAG,
      B.numRows(), B.numCols(), alpha,
      A.get_data(), A.lDim(),
      B.get_data(), B.lDim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  inline
  T
  normInf( MatrixWrapper<T> const & A ) {
    return lapack_wrapper::normInf( A.numRows(), A.numCols(), A.get_data(), A.lDim() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  inline
  T
  norm1( MatrixWrapper<T> const & A ) {
    return lapack_wrapper::norm1( A.numRows(), A.numCols(), A.get_data(), A.lDim() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  inline
  T
  normF( MatrixWrapper<T> const & A ) {
    return lapack_wrapper::normF( A.numRows(), A.numCols(), A.get_data(), A.lDim() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  inline
  T
  maxabs( MatrixWrapper<T> const & A ) {
    return lapack_wrapper::maxabs( A.numRows(), A.numCols(), A.get_data(), A.lDim() );
  }

  // explicit instantiation declaration to suppress warnings

  #ifdef LAPACK_WRAPPER_USE_CXX11
  extern template class MatrixWrapper<real>;
  extern template class MatrixWrapper<doublereal>;
  #endif

}

///
/// eof: wrapper.hxx
///
