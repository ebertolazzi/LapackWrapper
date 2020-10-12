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
    integer     m_dim;
    valueType * m_data;
  public:

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: build an empy wrapper
    \*/
    explicit
    DiagMatrixWrapper( )
    : m_dim(0)
    , m_data(nullptr)
    {}

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
     * Map a piece of memory into a matrix object
     * \param _data  pointer of the memory to be mapped as a matrix
     * \param _dim   number of elements on the diagonal of the mapped matrix
     */
    explicit
    DiagMatrixWrapper( valueType * _data, integer _dim ) {
      this->data = _data;
      this->dim  = _dim;
    }

    integer getDim()   const { return m_dim;}  //!< Number of elements
    integer numElems() const { return m_dim;}  //!< Number of elements

    valueType const * data() const { return m_data; }
    valueType       * data()       { return m_data; }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
     * Map a piece of memory into a matrix object
     *
     * \param _data pointer of the memory to be mapped as a matrix
     * \param _dim  dimension the mapped diagonal matrix
     */
    void
    setup( valueType * _data, integer _dim ) {
      m_data = _data;
      m_dim  = _dim;
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
      LW_ASSERT0(
        this->dim == COPY.dim,
        "DiagMatrixWrapper operator = bad matrix dimensions\n"
      );
      std::copy( COPY.data, COPY.data+COPY.dim, this->data );
      return *this;
    }

    void
    print( ostream_type & stream ) const {
      for ( integer i = 0; i < m_dim; ++i ) {
        for ( integer j = 0; j < m_dim; ++j ) {
          stream << std::setw(14);
          if ( i != j ) stream << '.';
          else          stream << m_data[i];
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

  //! base class to view a piece of memory as a matrix
  template <typename T>
  class MatrixWrapper {

    typedef T                    valueType;
    typedef MatrixWrapper<T>     MatW;
    typedef DiagMatrixWrapper<T> DiagW;
    typedef SparseMatrixBase<T>  Sparse;

  protected:

    integer     m_nRows;   //!< Number of rows
    integer     m_nCols;   //!< Number of columns
    integer     m_ldData;  //!< Leadind dimension
    valueType * m_data;    //!< pointer to matrix data

    #if defined(DEBUG) || defined(_DEBUG)
    integer
    iaddr( integer i,  integer j ) const {
      LW_ASSERT_TRACE(
        i >= 0 && i < m_nRows && j >= 0 && j < m_nCols,
        "iaddr({},{}) out of range [0,{}) x [0,{})\n",
        i, j, m_nRows, m_nCols
      );
      return i + j*m_ldData;
    }
    #else
    integer
    iaddr( integer i,  integer j ) const
    { return i + j*m_ldData; }
    #endif

    void
    check( Sparse const & sp ) const;

    void
    check( MatW const & M ) const;

  public:

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
     * build an empy wrapper
     */
    explicit
    MatrixWrapper( )
    : m_nRows(0)
    , m_nCols(0)
    , m_ldData(0)
    , m_data(nullptr)
    {
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
     *  Map a piece of memory into a matrix object
     *
     *  \param _data  pointer of the memory to be mapped as a matrix
     *  \param nr     number of rows of the mapped matrix
     *  \param nc     number of columns of the mapped matrix
     *  \param ld     leading dimension of the matrix (Fortran addressing 0 based)
     */
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
    : m_nRows(M.m_nRows)
    , m_nCols(M.m_nCols)
    , m_ldData(M.m_ldData)
    , m_data(M.m_data)
    { }

    MatrixWrapper<T> const &
    operator = ( MatrixWrapper<T> const & rhs ) {
      m_nRows  = rhs.m_nRows;
      m_nCols  = rhs.m_nCols;
      m_ldData = rhs.m_ldData;
      m_data   = rhs.m_data;
      return *this;
    }

    integer numRows()  const { return m_nRows; }  //!< Number of rows
    integer numCols()  const { return m_nCols; }  //!< Number of columns
    integer lDim()     const { return m_ldData; } //!< Leading dimension
    integer numElems() const { return m_nRows*m_nCols; } //!< Number of elements

    valueType const * data() const { return m_data; }
    valueType       * data()       { return m_data; }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
     * Map a piece of memory into a matrix object
     *
     * \param _data  pointer of the memory to be mapped as a matrix
     * \param nr     number of rows of the mapped matrix
     * \param nc     number of columns of the mapped matrix
     * \param ld     leading dimension of the matrix (Fortran addressing 0 based)
     */
    void
    setup(
      valueType * _data,
      integer     nr,
      integer     nc,
      integer     ld
    );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
     * Access element (i,j) of the matrix
     *
     * \param[in] i row of the element
     * \param[in] j column of the element
     */
    valueType const &
    operator () ( integer i,  integer j ) const
    { return m_data[this->iaddr(i,j)]; }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
     * Access element (i,j) of the matrix
     *
     * \param[in] i row of the element
     * \param[in] j column of the element
     */
    valueType &
    operator () ( integer i,  integer j )
    { return m_data[this->iaddr(i,j)]; }

    /*!
     * Fill the matrix with zeros
     */
    void
    zero_fill() {
      gezero( m_nRows, m_nCols, m_data, m_ldData );
    }

    /*!
     * Fill the matrix with zeros
     */
    void
    setZero() { zero_fill(); }

    /*!
     * Fill the matrix with value `val`
     * \param[in] val value used to fill matrix
     */
    void
    fill( valueType val ) {
      gefill( m_nRows, m_nCols, m_data, m_ldData, val );
    }

    /*!
     * Scale the matrix with value `sc` (multiply all elements by `sc`)
     * \param[in] sc value used to scale matrix
     */
    void
    scale_by( valueType sc );

    /*!
     * Zeroes a rectangular block of the stored matrix
     * staring at `(irow,icol)` position
     *
     * \param[in] nr    number of rows of the block to be zeroed
     * \param[in] nc    number of columns of the block to be zeroed
     * \param[in] irow  starting row
     * \param[in] icol  stating column
     */
    void
    zero_block(
      integer nr,
      integer nc,
      integer irow,
      integer icol
    ) {
      gezero( nr, nc, m_data + this->iaddr(irow,icol), m_ldData );
    }

    /*!
     * Initialize the matrix as an identity matrix
     */
    void
    id( valueType dg ) {
      geid( m_nRows, m_nCols, m_data, m_ldData, dg );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void
    add_to_diag( valueType mu ) {
      axpy(
        std::min(m_nRows,m_nCols),
        1.0, &mu, 0,
        m_data, m_ldData+1
      );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
     * Initialize matrix
     *
     * \param[in] data   pointer of memory with data to be copied
     * \param[in] ldData leading dimension of the memory to be copied
     *
     */
    void load( valueType const data[], integer ldData );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
     * Initialize matrix
     *
     * \param[in] data   pointer of memory with data to be copied
     * \param[in] ldData leading dimension of the memory to be copied
     *
     */
    void load_transposed( valueType const data[], integer ldData );

    /*!
     * Initialize matrix
     *
     * \param[in] A initialize matrix with the matrix of the object A
     *
     */
    void load( MatW const & A );

    /*!
     * Initialize matrix
     *
     * \param[in] A initialize matrix with the matrix of the object A
     *
     */
    void load_transposed( MatW const & A );

    /*!
     *  Copy a matrix to a rectangular block of the stored matrix
     *  staring at `(irow,icol)` position
     *
     *  \param[in] B     reference to the matrix of the object B
     *  \param[in] irow  starting row
     *  \param[in] icol  stating column
     */
    void
    load_block(
      MatW const & B,
      integer      irow = 0,
      integer      icol = 0
    ) {
      this->load_block(
        B.m_nRows, B.m_nCols, B.m_data, B.m_ldData, irow, icol
      );
    }

    /*!
     *  Copy a matrix to a rectangular block of the stored matrix
     *  staring at `(irow,icol)` position
     *
     *  \param[in] B     reference to the matrix of the object B
     *  \param[in] irow  starting row
     *  \param[in] icol  stating column
     */
    void
    load_block_transposed(
      MatW const & B,
      integer      irow = 0,
      integer      icol = 0
    ) {
      this->load_block_transposed(
        B.m_nRows, B.m_nCols, B.m_data, B.m_ldData, irow, icol
      );
    }

    /*!
     *  Copy a matrix to a rectangular block of the stored matrix
     *  staring at `(irow,icol)` position
     *
     *  \param[in] nr    number of rows of the block to be zeroed
     *  \param[in] nc    number of columns of the block to be zeroed
     *  \param[in] B     pointer to memory storing the input matrix `B`
     *  \param[in] ldB   leading dimension of the matrix `B`
     *  \param[in] irow  starting row
     *  \param[in] icol  stating column
     */
    void
    load_block(
      integer         nr,
      integer         nc,
      valueType const B[],
      integer         ldB,
      integer         irow = 0,
      integer         icol = 0
    ) {
      LW_ASSERT(
        irow + nr <= m_nRows &&
        icol + nc <= m_nCols &&
        irow >= 0 && icol >= 0,
        "load_block( nr = {} nc = {},..., irow = {}, icol = {}) bad parameters\n",
        nr, nc, irow, icol
      );
      integer info = gecopy(
        nr, nc,
        B, ldB,
        m_data + this->iaddr(irow,icol), m_ldData
      );
      LW_ASSERT(
        info == 0,
        "load_block call lapack_wrapper::gecopy return info = {}\n", info
      );
    }

    /*!
     *  Copy a matrix to a rectangular block of the stored matrix
     *  staring at `(irow,icol)` position
     *
     *  \param[in] nr    number of rows of the block to be zeroed
     *  \param[in] nc    number of columns of the block to be zeroed
     *  \param[in] B     pointer to memory storing the input matrix `B`
     *  \param[in] ldB   leading dimension of the matrix `B`
     *  \param[in] irow  starting row
     *  \param[in] icol  stating column
     */
    void
    load_block_transposed(
      integer         nr,
      integer         nc,
      valueType const B[],
      integer         ldB,
      integer         irow = 0,
      integer         icol = 0
    ) {
      LW_ASSERT(
        irow + nc <= m_nRows &&
        icol + nr <= m_nCols &&
        irow >= 0 && icol >= 0,
        "load_block_transpose( nr = {}, nc = {},..., irow = {}, icol = {} ) "
        "bad parameters\n",
        nr, nc, irow, icol
      );
      valueType const * pd = B;
      valueType       * pp = m_data + this->iaddr(irow,icol);
      for ( integer i = 0; i < nc; ++i, pd += ldB, ++pp )
        lapack_wrapper::copy( nr, pd, 1, pp, m_ldData );
    }

    /*!
     *  Copy a matrix to a rectangular block of the stored matrix
     *  staring at `(irow,icol)` position
     *
     *  \param[in] n     number of element of the diagonal
     *  \param[in] D     pointer to memory storing the diagonal
     *  \param[in] irow  starting row
     *  \param[in] icol  stating column
     */
    void
    load_diagonal_block(
      integer         n,
      valueType const D[],
      integer         irow = 0,
      integer         icol = 0
    ) {
      LW_ASSERT(
        irow + n <= m_nRows &&
        icol + n <= m_nCols &&
        irow >= 0 && icol >= 0,
        "load_diagonal_block( n = {},..., irow = {}, icol = {}) "
        "bad parameters\n",
        n, irow, icol
      );
      for ( integer i = 0; i < n; ++i )
        m_data[this->iaddr(irow+i,icol+i)] = D[i];
    }

    /*!
     *  Copy a matrix to a rectangular block of the stored matrix
     *  staring at `(irow,icol)` position
     *
     *  \param[in] D     reference to the matrix of the object D
     *  \param[in] irow  starting row
     *  \param[in] icol  stating column
     */
    void
    load_block(
      DiagW const & D,
      integer       irow = 0,
      integer       icol = 0
    ) {
      this->load_diagonal_block( D.getDim(), D.data(), irow, icol );
    }

    /*!
     *  Copy vector `column` to the `icol`th column of the internal stored matrix
     *  \param[in] column the column vector
     *  \param[in] icol   the column to be changed
     */
    void
    load_column( valueType const column[], integer icol )
    { copy( m_nRows, column, 1, m_data + icol * m_ldData, 1 ); }

    /*!
     *  Copy vector `row` to the `irow`th row of the internal stored matrix
     *  \param[in] row  the row vector
     *  \param[in] irow the row to be changed
     */
    void
    load_row( valueType const row[], integer irow )
    { copy( m_nCols, row, 1, m_data + irow, m_ldData ); }

    void
    load_sparse_column(
      integer         nnz,
      valueType const values[],
      integer   const i_row[],
      integer         j
    );

    void
    load_sparse_row(
      integer         nnz,
      valueType const values[],
      integer         i,
      integer   const j_col[]
    );

    /*!
     * Copy sparse matrix into the object
     *
     * \param[in] sp sparse matrix to be copied
     *
     */
    void load( Sparse const & sp );
    /*!
     * Copy sparse matrix into the object transposing it
     *
     * \param[in] sp sparse matrix to be copied
     *
     */
    void load_transposed( Sparse const & sp );

    /*!
     * Copy sparse matrix into the object
     *
     * \param[in] rows row indices of the sparse matrix
     * \param[in] cols column indices of the sparse matrix
     * \param[in] vals values of the sparse matrix
     * \param[in] nnz  number of nonzeros elements
     *
     */
    void
    load(
      integer   const rows[],
      integer   const cols[],
      valueType const vals[],
      integer         nnz
    );

    void
    load0(
      integer   const rows[],
      integer   const cols[],
      valueType const vals[],
      integer         nnz
    ) {
      this->zero_fill();
      this->load( rows, cols, vals, nnz );
    }

    /*!
     * Copy sparse matrix into the object
     *
     * \param[in] rows row indices of the sparse matrix
     * \param[in] cols column indices of the sparse matrix
     * \param[in] vals values of the sparse matrix
     * \param[in] nnz  number of nonzeros elements
     *
     */
    void
    load_symmetric(
      integer   const rows[],
      integer   const cols[],
      valueType const vals[],
      integer         nnz
    );

    void
    load0_symmetric(
      integer   const rows[],
      integer   const cols[],
      valueType const vals[],
      integer         nnz
    ) {
      this->zero_fill();
      this->load_symmetric( rows, cols, vals, nnz );
    }

    /*!
     * Copy sparse matrix into the object
     *
     * \param sp     sparse matrix to be copied
     * \param i_offs offset added to to the row indices
     * \param j_offs offset added to to the column indices
     *
     */
    void
    load( Sparse const & sp, integer i_offs, integer j_offs );

    void
    load0( Sparse const & sp, integer i_offs, integer j_offs ) {
      this->zero_fill();
      this->load( sp, i_offs, j_offs );
    }

    /*!
     * Copy sparse matrix into the object transposing it
     *
     * \param sp     sparse matrix to be copied
     * \param i_offs offset added to to the row indices
     * \param j_offs offset added to to the column indices
     *
     */
    void
    load_transposed( Sparse const & sp, integer i_offs, integer j_offs );

    void
    load0_transposed( Sparse const & sp, integer i_offs, integer j_offs ) {
      this->zero_fill();
      this->load_transposed( sp, i_offs, j_offs );
    }

    /*!
     * Copy sparse matrix into the object
     *
     * \param[in] i_offs row index offset
     * \param[in] j_offs column index offset
     * \param[in] rows   row indices of the sparse matrix
     * \param[in] cols   column indices of the sparse matrix
     * \param[in] vals   values of the sparse matrix
     * \param[in] nnz    number of nonzeros elements
     *
     */
    void
    load(
      integer         i_offs,
      integer         j_offs,
      integer   const rows[],
      integer   const cols[],
      valueType const vals[],
      integer         nnz
    );

    void
    load0(
      integer         i_offs,
      integer         j_offs,
      integer   const rows[],
      integer   const cols[],
      valueType const vals[],
      integer         nnz
    ) {
      this->zero_fill();
      this->load( i_offs, j_offs, rows, cols, vals, nnz );
    }

    /*!
     * Copy sparse matrix into the object
     *
     * \param[in] i_offs row index offset
     * \param[in] j_offs column index offset
     * \param[in] rows   row indices of the sparse matrix
     * \param[in] cols   column indices of the sparse matrix
     * \param[in] vals   values of the sparse matrix
     * \param[in] nnz    number of nonzeros elements
     *
     */
    void
    load_symmetric(
      integer         i_offs,
      integer         j_offs,
      integer   const rows[],
      integer   const cols[],
      valueType const vals[],
      integer         nnz
    );

    void
    load0_symmetric(
      integer         i_offs,
      integer         j_offs,
      integer   const rows[],
      integer   const cols[],
      valueType const vals[],
      integer         nnz
    ) {
      this->zero_fill();
      this->load_symmetric( i_offs, j_offs, rows, cols, vals, nnz );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
     * Add a block of memory to the matrix scaled by alpha
     *
     * \param alpha  scaling factor
     * \param data   pointer of memory with data to be copied
     * \param ldData leading dimension of the memory to be copied
     *
     */
    void add( valueType alpha, valueType const data[], integer ldData );

    /*!
     * Add a block of memory to the matrix
     *
     * \param data   pointer of memory with data to be copied
     * \param ldData leading dimension of the memory to be copied
     *
     */
    void add( valueType const data[], integer ldData );

    // R <- R + alpha * A
    /*!
     *
     * \param A first mapped matrix
     *
     */
    void
    add( valueType alpha, MatW const & A ) {
      #if defined(DEBUG) || defined(_DEBUG)
      LW_ASSERT0(
        m_nRows == A.m_nRows && m_nCols == A.m_nCols,
        "add(..) incompatible dimensions\n"
      );
      #endif
      this->add( alpha, A.m_data, A.m_ldData );
    }

    void
    add( MatW const & A ) {
      #if defined(DEBUG) || defined(_DEBUG)
      LW_ASSERT0(
        m_nRows == A.m_nRows && m_nCols == A.m_nCols,
        "add(..) incompatible dimensions\n"
      );
      #endif
      this->add( A.m_data, A.m_ldData );
    }

    /*!
     * Add sparse matrix to the object
     *
     * \param sp sparse matrix to be copied
     *
     */
    void add( Sparse const & sp );

    /*!
     * Add sparse matrix to the object
     *
     * \param sp     sparse matrix to be copied
     * \param i_offs offset added to to the row indices
     * \param j_offs offset added to to the column indices
     *
     */
    void add( Sparse const & sp, integer i_offs, integer j_offs );

    /*!
     * Add sparse multiplied by `alpha` matrix to the object
     *
     * \param alpha  scalar used to multiply the sparse matrix
     * \param sp     sparse matrix to be copied
     *
     */
    void add( valueType alpha, Sparse const & sp );

    /*!
     * Add sparse multiplied by `alpha` matrix to the object
     *
     * \param alpha  scalar used to multiply the sparse matrix
     * \param sp     sparse matrix to be copied
     * \param i_offs offset added to to the row indices
     * \param j_offs offset added to to the column indices
     *
     */
    void
    add(
      valueType      alpha,
      Sparse const & sp,
      integer        i_offs,
      integer        j_offs
    );

    // alpha*A + beta*B -> C
    /*!
     * Add a linear combination of two matrix and assign to a third one
     *
     * \param alpha scalar used to multiply the matrix `A`
     * \param A     first mapped matrix
     * \param beta  scalar used to multiply the matrix `B`
     * \param B     second mapped matrix
     * \param C     mapped matrix which store the result
     *
     */
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
     * make a view of a block
     *
     * \param i_offs initial row of the block to be extracted
     * \param j_offs initial column of the block to be extracted
     * \param nrow   number of rows of the block to be extracted
     * \param ncol   number of columns of the block to be extracted
     * \param to     mapped matrix which store the result
     *
     */
    void
    view_block(
      integer i_offs,
      integer j_offs,
      integer nrow,
      integer ncol,
      MatW &  to
    ) {
      if ( nrow > 0 && ncol > 0 ) {
        #if defined(DEBUG) || defined(_DEBUG)
        LW_ASSERT_TRACE(
          i_offs >= 0 && i_offs+nrow <= m_nRows &&
          j_offs >= 0 && j_offs+ncol <= m_nCols,
          "view_block(i_offs={}, j_offs={}, nrow={}, ncol={}) "
          "will be out of range [0,{}) x [0,{})\n",
          i_offs, j_offs, nrow, ncol, m_nRows, m_nCols
        );
        #endif
        to.setup(
          m_data+this->iaddr(i_offs,j_offs),
          nrow, ncol, m_ldData
        );
      }
    }

    /*!
     * Extract a matrix in transposed form
     *
     * \param out mapped matrix which store the result
     *
     */
    void
    get_transposed( MatW & out ) {
      #if defined(DEBUG) || defined(_DEBUG)
      LW_ASSERT0(
        m_nRows == out.m_nCols && m_nCols == out.m_nRows,
        "get_transposed(...) incompatible matrices\n"
      );
      #endif
      valueType const * pc = m_data;
      for ( integer i = 0; i < m_nCols; ++i, pc += m_ldData )
        lapack_wrapper::copy( m_nRows, pc, 1, out.m_data + i, out.m_ldData );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
     * Extract a block
     *
     * \param i_offs initial row of the block to be extracted
     * \param j_offs initial column of the block to be extracted
     * \param to     mapped matrix which store the result
     *
     */
    void
    get_block(
      MatW &  to,
      integer i_offs,
      integer j_offs
    ) {
      #if defined(DEBUG) || defined(_DEBUG)
      LW_ASSERT_TRACE(
        i_offs >= 0 && i_offs+to.m_nRows <= m_nRows &&
        j_offs >= 0 && j_offs+to.m_nCols <= m_nCols,
        "get_block(i_offs={}, j_offs={}, ...) "
        "will be out of range [0,{}) x [0,{})\n",
        i_offs, j_offs, m_nRows, m_nCols
      );
      #endif
      lapack_wrapper::gecopy(
        to.m_nRows, to.m_nCols,
        m_data+this->iaddr(i_offs,j_offs), m_ldData,
        to.m_data, to.m_ldData
      );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
     * Extract a block
     *
     * \param i_offs initial row of the block to be extracted
     * \param j_offs initial column of the block to be extracted
     * \param to     mapped matrix which store the result
     *
     */
    void
    get_block_transposed(
      MatW &  to,
      integer i_offs,
      integer j_offs
    ) {
      #if defined(DEBUG) || defined(_DEBUG)
      LW_ASSERT_TRACE(
        i_offs >= 0 && i_offs+to.m_nCols <= m_nRows &&
        j_offs >= 0 && j_offs+to.m_nRows <= m_nCols,
        "get_block_transposed(i_offs={}, j_offs={}, ...) "
        "will be out of range [0,{}) x [0,{})\n",
        i_offs, j_offs, m_nRows, m_nCols
      );
      #endif
      valueType const * pc = m_data+this->iaddr(i_offs,j_offs);
      for ( integer i = 0; i < to.m_nCols; ++i, pc += m_ldData )
        lapack_wrapper::copy( to.m_nRows, pc, 1, to.m_data + i, to.m_ldData );
    }

    void
    print( ostream_type & stream ) const {
      for ( integer i = 0; i < m_nRows; ++i ) {
        for ( integer j = 0; j < m_nCols; ++j )
          stream << std::setw(14) << (*this)(i,j) << ' ';
        stream << '\n';
      }
    }

    void
    print0( ostream_type & stream, valueType eps ) const {
      for ( integer i = 0; i < m_nRows; ++i ) {
        for ( integer j = 0; j < m_nCols; ++j ) {
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
   * Perform matrix vector multiplication `c = beta*c + alpha*A*v`
   *
   * \param alpha  matrix `A` is multiplied by `alpha`
   * \param TRANSA choose if multiply with the transpose of `A`
   * \param A      matrix used in the multiplication
   * \param v      vector to be multiplied
   * \param incv   stride of the vector `v`
   * \param beta   scalar used to multiply `c`
   * \param c      result vector
   * \param incc   stride of the vector `c`
   *
   */
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
      A.data(), A.lDim(),
      v, incv,
      beta,
      c, incc
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // c = beta*c + alpha*A*v
  /*!
   * Perform matrix vector multiplication `c = beta*c + alpha*A*v`
   *
   * \param alpha  matrix `A` is multiplied by `alpha`
   * \param A      matrix used in the multiplication
   * \param v      vector to be multiplied
   * \param incv   stride of the vector `v`
   * \param beta   scalar used to multiply `c`
   * \param c      result vector
   * \param incc   stride of the vector `c`
   *
   */
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
      A.data(), A.lDim(),
      v, incv,
      beta,
      c, incc
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*!
   * Perform matrix matrix multiplication `C = beta*C + alpha*A*B`
   *
   * \param where  message added in case of error
   * \param alpha  matrix `A` is multiplied by `alpha`
   * \param A      matrix used in the multiplication
   * \param B      matrix used in the multiplication
   * \param beta   scalar used to multiply `C`
   * \param C      result matrix
   *
   */
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
    LW_ASSERT_DEBUG(
      A.numCols() == B.numRows() &&
      A.numRows() == C.numRows() &&
      B.numCols() == C.numCols(),
      "gemm, at `{}' inconsistent dimensions:\n"
      "A = {} x {}\nB = {} x {}\nC = {} x {}\n",
      where,
      A.numRows(), A.numCols(),
      B.numRows(), B.numCols(),
      C.numRows(), C.numCols()
    );
    lapack_wrapper::gemm(
      NO_TRANSPOSE,
      NO_TRANSPOSE,
      A.numRows(), B.numCols(), A.numCols(),
      alpha,
      A.data(), A.lDim(),
      B.data(), B.lDim(),
      beta,
      C.data(), C.lDim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*!
   * Perform matrix matrix multiplication `C = beta*C + alpha*A*B`
   *
   * \param alpha  matrix `A` is multiplied by `alpha`
   * \param A      matrix used in the multiplication
   * \param B      matrix used in the multiplication
   * \param beta   scalar used to multiply `C`
   * \param C      result matrix
   *
   */
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
    LW_ASSERT_DEBUG(
      A.numCols() == B.numRows() &&
      A.numRows() == C.numRows() &&
      B.numCols() == C.numCols(),
      "gemm, inconsistent dimensions:\n"
      "A = {} x {}\nB = {} x {}\nC = {} x {}\n",
      A.numRows(), A.numCols(), B.numRows(), B.numCols(), C.numRows(), C.numCols()
    );
    lapack_wrapper::gemm(
      NO_TRANSPOSE,
      NO_TRANSPOSE,
      A.numRows(), B.numCols(), A.numCols(),
      alpha,
      A.data(), A.lDim(),
      B.data(), B.lDim(),
      beta,
      C.data(), C.lDim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*!
   * Perform matrix matrix multiplication `C = beta*C + alpha*A*B`
   *
   * \param where  message added in case of error
   * \param alpha  matrix `A` is multiplied by `alpha`
   * \param TRANSA choose `A` transposed or not
   * \param A      matrix used in the multiplication
   * \param TRANSB choose `B` transposed or not
   * \param B      matrix used in the multiplication
   * \param beta   scalar used to multiply `C`
   * \param C      result matrix
   *
   */
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
    LW_ASSERT_DEBUG(
      C.numRows() == Ar && C.numCols() == Bc && Ac == Br,
      "gemm, at `{}' inconsistent dimensions:"
      "\nA = {} x {}\nB = {} x {}\nC = {} x {}"
      "\nA {} transposed\nB {} transposed\n",
      where,
      A.numRows(), A.numCols(),
      B.numRows(), B.numCols(),
      C.numRows(), C.numCols(),
      (NO_TRANSPOSE?"NO":""),
      (NO_TRANSPOSE?"NO":"")
    );
    lapack_wrapper::gemm(
      TRANSA,
      TRANSB,
      C.numRows(), C.numCols(), Ac,
      alpha,
      A.data(), A.lDim(),
      B.data(), B.lDim(),
      beta,
      C.data(), C.lDim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // C = beta*C + alpha*A*B (NO DEBUG)
  /*!
   * Perform matrix matrix multiplication `C = beta*C + alpha*A*B`
   *
   * \param alpha  matrix `A` is multiplied by `alpha`
   * \param TRANSA choose `A` transposed or not
   * \param A      matrix used in the multiplication
   * \param TRANSB choose `B` transposed or not
   * \param B      matrix used in the multiplication
   * \param beta   scalar used to multiply `C`
   * \param C      result matrix
   *
   */
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
    LW_ASSERT_DEBUG(
      C.numRows() == Ar && C.numCols() == Bc && Ac == Br,
      "gemm, inconsistent dimensions:"
      "\nA = {} x {}\nB = {} x {}\nC = {} x {}"
      "\nA {} transposed\nB {} transposed\n",
      A.numRows(), A.numCols(),
      B.numRows(), B.numCols(),
      C.numRows(), C.numCols(),
      (NO_TRANSPOSE?"NO":""),
      (NO_TRANSPOSE?"NO":"")
    );
    lapack_wrapper::gemm(
      TRANSA,
      TRANSB,
      C.numRows(), C.numCols(), Ac,
      alpha,
      A.data(), A.lDim(),
      B.data(), B.lDim(),
      beta,
      C.data(), C.lDim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // x = A * x or x = A^T * x
  /*!
   * Perform matrix vector multiplication `x = A * x` or `x = A^T * x`
   *
   * \param UPLO   choose upper or lower part
   * \param TRANS  choose `A` transposed or not
   * \param DIAG   use or not diagonal elements to 1
   * \param x      vector used in the multiplcation and result
   * \param incx   stride of vector `x`
   *
   */
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
    LW_ASSERT_DEBUG(
      A.numRows() == A.numCols(),
      "trmv, matrix is {} x {} expected square\n",
      A.numRows(), A.numCols()
    );
    lapack_wrapper::trmv(
      UPLO, TRANS, DIAG, A.numRows(), A.data(), A.lDim(), x, incx
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // x = A * x or x = A^T * x
  /*!
   * Perform matrix vector multiplication `x = A * x` or `x = A^T * x`
   *
   * \param where  message added in case of error
   * \param UPLO   choose upper or lower part
   * \param TRANS  choose `A` transposed or not
   * \param DIAG   use or not diagonal elements to 1
   * \param A      matrix to multiplied
   * \param x      vector used in the multiplcation and result
   * \param incx   stride of vector `x`
   *
   */
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
    LW_ASSERT_DEBUG(
      A.numRows() == A.numCols(),
      "trmv at `{}` matrix is {} x {} expected square\n",
      where, A.numRows(), A.numRows()
    );
    lapack_wrapper::trmv(
      UPLO, TRANS, DIAG, A.numRows(), A.data(), A.lDim(), x, incx
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // x = A^(-1) * x or x = A^(-T) * x
  /*!
   * Perform matrix vector multiplication `x = A * x` or `x = A^T * x`
   *
   * \param UPLO   choose upper or lower part
   * \param TRANS  choose `A` transposed or not
   * \param DIAG   use or not diagonal elements to 1
   * \param A      matrix to multiplied
   * \param x      vector used in the multiplcation and result
   * \param incx   stride of vector `x`
   *
   */
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
    LW_ASSERT_DEBUG(
      A.numRows() == A.numCols(),
      "trsv, matrix is {} x {} expected square\n",
      A.numRows(), A.numCols()
    );
    lapack_wrapper::trsv(
      UPLO, TRANS, DIAG, A.numRows(), A.data(), A.lDim(), x, incx
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // x = A^(-1) * x or x = A^(-T) * x
  /*!
   * Perform matrix vector multiplication `x = A^(-1) * x` or `x = A^(-T) * x`
   *
   * \param where  message added in case of error
   * \param UPLO   choose upper or lower part
   * \param TRANS  choose `A` transposed or not
   * \param DIAG   use or not diagonal elements to 1
   * \param A      matrix to multiplied
   * \param x      vector used in the multiplcation and result
   * \param incx   stride of vector `x`
   *
   */
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
    LW_ASSERT_DEBUG(
      A.numRows() == A.numCols(),
      "trsv at `{}` matrix is {} x {} expected square\n",
      where, A.numRows(), A.numRows()
    );
    lapack_wrapper::trsv(
      UPLO, TRANS, DIAG, A.numRows(), A.data(), A.lDim(), x, incx
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*!
   * Perform matrix matrix multiplication
   *
   *    B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
   *    op( A ) = A   or   op( A ) = A**T.
   *
   * \param SIDE   multiply on left or right
   * \param UPLO   choose upper or lower part
   * \param TRANS  choose `A` transposed or not
   * \param DIAG   use or not diagonal elements to 1
   * \param alpha  scalar used in the multiplication
   * \param A      matrix to be multiplied
   * \param B      matrix to be multiplied
   *
   */
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
    LW_ASSERT_DEBUG(
      A.numRows() == A.numCols(),
      "trmm, matrix is {} x {} expected square\n",
      A.numRows(), A.numCols()
    );
    lapack_wrapper::trmm(
      SIDE, UPLO, TRANS, DIAG,
      B.numRows(), B.numCols(), alpha,
      A.data(), A.lDim(),
      B.data(), B.lDim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*!
   * Perform matrix matrix multiplication
   *
   *    B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
   *    op( A ) = A   or   op( A ) = A**T.
   *
   * \param where  message added in case of error
   * \param SIDE   multiply on left or right
   * \param UPLO   choose upper or lower part
   * \param TRANS  choose `A` transposed or not
   * \param DIAG   use or not diagonal elements to 1
   * \param alpha  scalar used in the multiplication
   * \param A      matrix to be multiplied
   * \param B      matrix to be multiplied
   *
   */
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
    LW_ASSERT_DEBUG(
      A.numRows() == A.numCols(),
      "trmm, at `{}` matrix is {} x {} expected square\n",
      where, A.numRows(), A.numRows()
    );
    lapack_wrapper::trmm(
      SIDE, UPLO, TRANS, DIAG,
      B.numRows(), B.numCols(), alpha,
      A.data(), A.lDim(),
      B.data(), B.lDim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*!
   * Perform matrix matrix multiplication
   *
   *    B := alpha*op( A^(-1) )*B,   or   B := alpha*B*op( A^(-1) ),
   *    op( A ) = A   or   op( A ) = A**T.
   *
   * \param SIDE   multiply on left or right
   * \param UPLO   choose upper or lower part
   * \param TRANS  choose `A` transposed or not
   * \param DIAG   use or not diagonal elements to 1
   * \param alpha  scalar used in the multiplication
   * \param A      matrix to be multiplied
   * \param B      matrix to be multiplied
   *
   */
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
    LW_ASSERT_DEBUG(
      A.numRows() == A.numCols(),
      "trsm, matrix is {} x {} expected square\n",
      A.numRows(), A.numCols()
    );
    lapack_wrapper::trsm(
      SIDE, UPLO, TRANS, DIAG,
      B.numRows(), B.numCols(), alpha,
      A.data(), A.lDim(),
      B.data(), B.lDim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*!
   * Perform matrix matrix multiplication
   *
   *    B := alpha*op( A^(-1) )*B,   or   B := alpha*B*op( A^(-1) ),
   *    op( A ) = A   or   op( A ) = A**T.
   *
   * \param where  message added in case of error
   * \param SIDE   multiply on left or right
   * \param UPLO   choose upper or lower part
   * \param TRANS  choose `A` transposed or not
   * \param DIAG   use or not diagonal elements to 1
   * \param alpha  scalar used in the multiplication
   * \param A      matrix to be multiplied
   * \param B      matrix to be multiplied
   *
   */
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
    LW_ASSERT_DEBUG(
      A.numRows() == A.numCols(),
      "trmm, at `{}` matrix is {} x {} expected square\n",
      where, A.numRows(), A.numRows()
    );
    lapack_wrapper::trsm(
      SIDE, UPLO, TRANS, DIAG,
      B.numRows(), B.numCols(), alpha,
      A.data(), A.lDim(),
      B.data(), B.lDim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  inline
  T
  normInf( MatrixWrapper<T> const & A ) {
    return lapack_wrapper::normInf( A.numRows(), A.numCols(), A.data(), A.lDim() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  inline
  T
  norm1( MatrixWrapper<T> const & A ) {
    return lapack_wrapper::norm1( A.numRows(), A.numCols(), A.data(), A.lDim() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  inline
  T
  normF( MatrixWrapper<T> const & A ) {
    return lapack_wrapper::normF( A.numRows(), A.numCols(), A.data(), A.lDim() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  inline
  T
  maxabs( MatrixWrapper<T> const & A ) {
    return lapack_wrapper::maxabs( A.numRows(), A.numCols(), A.data(), A.lDim() );
  }

  // explicit instantiation declaration to suppress warnings
  extern template class MatrixWrapper<real>;
  extern template class MatrixWrapper<doublereal>;

}

///
/// eof: wrapper.hxx
///
