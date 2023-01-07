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
    using real_type = T;
    using DMatW     = DiagMatrixWrapper<T>;

  protected:
    integer     m_dim{0};
    real_type * m_data{nullptr};
  public:

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //!
    //! Build an empty wrapper.
    //!
    explicit DiagMatrixWrapper( ) = default;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //!
    //! Map a piece of memory into a matrix object.
    //!
    //! \param data pointer of the memory to be mapped as a matrix
    //! \param dim  number of elements on the diagonal of the mapped matrix
    //!
    explicit
    DiagMatrixWrapper( real_type * data, integer dim )
    : m_dim(dim)
    , m_data(data)
    { }

    integer getDim()   const { return m_dim;}  //!< Number of elements
    integer numElems() const { return m_dim;}  //!< Number of elements

    real_type const * data() const { return m_data; }
    real_type       * data()       { return m_data; }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //!
    //! Map a piece of memory into a matrix object.
    //!
    //! \param data pointer of the memory to be mapped as a matrix
    //! \param dim  dimension the mapped diagonal matrix
    //!
    void
    setup( real_type * data, integer dim ) {
      m_data = data;
      m_dim  = dim;
    }

    real_type const &
    operator [] ( integer i ) const
    { return m_data[i]; }

    real_type &
    operator [] ( integer i )
    { return m_data[i]; }

    DMatW const & operator = (real_type v) {
      fill( m_dim, m_data, v );
      return *this;
    }

    DMatW & operator = ( DMatW const & COPY ) {
      UTILS_ASSERT0(
        m_dim == COPY.m_dim,
        "DiagMatrixWrapper operator = bad matrix dimensions\n"
      );
      std::copy_n( COPY.m_data, COPY.m_dim, m_data );
      return *this;
    }

    string
    to_string() const {
      string res = "";
      for ( integer i = 0; i < m_dim; ++i ) {
        for ( integer j = 0; j < m_dim; ++j ) {
          if ( i != j ) res += ".              ";
          else          res += fmt::format( "{:14} ", m_data[i] );
        }
        res += '\n';
      }
      return res;
    }

  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  inline
  ostream_type &
  operator << ( ostream_type & stream, DiagMatrixWrapper<T> const & A ) {
    stream << A.to_string();
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

  //!
  //! Base class to view a piece of memory as a matrix.
  //!
  template <typename T>
  class MatrixWrapper {

    using real_type = T;
    using MatW      = MatrixWrapper<T>;
    using DiagW     = DiagMatrixWrapper<T>;
    using Sparse    = SparseMatrixBase<T>;

  protected:

    integer     m_nrows;   //!< Number of rows
    integer     m_ncols;   //!< Number of columns
    integer     m_ldData;  //!< Leadind dimension
    real_type * m_data;    //!< pointer to matrix data

    #if defined(DEBUG) || defined(_DEBUG)
    integer
    iaddr( integer i,  integer j ) const {
      UTILS_ASSERT_TRACE(
        i >= 0 && i < m_nrows && j >= 0 && j < m_ncols,
        "iaddr({},{}) out of range [0,{}) x [0,{})\n",
        i, j, m_nrows, m_ncols
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
    //!
    //! Build an empty wrapper.
    //!
    explicit
    MatrixWrapper( )
    : m_nrows(0)
    , m_ncols(0)
    , m_ldData(0)
    , m_data(nullptr)
    {
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //!
    //! Map a piece of memory into a matrix object.
    //!
    //! \param data  pointer of the memory to be mapped as a matrix
    //! \param nr    number of rows of the mapped matrix
    //! \param nc    number of columns of the mapped matrix
    //! \param ld    leading dimension of the matrix (Fortran addressing 0 based)
    //!
    explicit
    MatrixWrapper(
      real_type * data,
      integer     nr,
      integer     nc,
      integer     ld
    );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    explicit
    MatrixWrapper( MatrixWrapper<T> const & M )
    : m_nrows(M.m_nrows)
    , m_ncols(M.m_ncols)
    , m_ldData(M.m_ldData)
    , m_data(M.m_data)
    { }

    MatrixWrapper<T> const &
    operator = ( MatrixWrapper<T> const & rhs ) {
      m_nrows  = rhs.m_nrows;
      m_ncols  = rhs.m_ncols;
      m_ldData = rhs.m_ldData;
      m_data   = rhs.m_data;
      return *this;
    }

    integer nrows()     const { return m_nrows; }  //!< Number of rows
    integer ncols()     const { return m_ncols; }  //!< Number of columns
    integer ldim()      const { return m_ldData; } //!< Leading dimension
    integer num_elems() const { return m_nrows*m_ncols; } //!< Number of elements

    // ALIAS
    #ifdef LAPACK_WRAPPER_USE_ALIAS
    integer numRows()  const { return nrows(); }  //!< Number of rows
    integer numCols()  const { return ncols(); }  //!< Number of columns
    integer ldim()     const { return lDim(); } //!< Leading dimension
    integer numElems() const { return num_elems(); } //!< Number of elements
    #endif

    real_type const * data() const { return m_data; }
    real_type       * data()       { return m_data; }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //!
    //! Map a piece of memory into a matrix object.
    //!
    //! \param data  pointer of the memory to be mapped as a matrix
    //! \param nr    number of rows of the mapped matrix
    //! \param nc    number of columns of the mapped matrix
    //! \param ld    leading dimension of the matrix (Fortran addressing 0 based)
    //!
    void
    setup(
      real_type * data,
      integer     nr,
      integer     nc,
      integer     ld
    );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //!
    //! Access element (i,j) of the matrix.
    //!
    //! \param[in] i row of the element
    //! \param[in] j column of the element
    //!
    real_type const &
    operator () ( integer i, integer j ) const
    { return m_data[this->iaddr(i,j)]; }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //!
    //! Access element (i,j) of the matrix.
    //!
    //! \param[in] i row of the element
    //! \param[in] j column of the element
    //!
    real_type &
    operator () ( integer i, integer j )
    { return m_data[this->iaddr(i,j)]; }

    //!
    //! Fill the matrix with zeros.
    //!
    void
    zero_fill() {
      gezero( m_nrows, m_ncols, m_data, m_ldData );
    }

    //!
    //! Fill the matrix with zeros.
    //!
    void
    setZero() { zero_fill(); }

    //!
    //! Fill the matrix with value `val`.
    //!
    //! \param[in] val value used to fill matrix
    //!
    void
    fill( real_type val ) {
      gefill( m_nrows, m_ncols, m_data, m_ldData, val );
    }

    //!
    //! Scale the matrix with value `sc` (multiply all elements by `sc`).
    //!
    //! \param[in] sc value used to scale matrix
    //!
    void
    scale_by( real_type sc );

    //!
    //! Scale the matrix with value `sc`
    //! (multiply all elements of the block by `sc`).
    //!
    //! \param[in] sc    value used to scale matrix
    //! \param[in] nr    number of rows of the block to be zeroed
    //! \param[in] nc    number of columns of the block to be zeroed
    //! \param[in] irow  starting row
    //! \param[in] icol  stating column
    //!
    void
    scale_block(
      real_type sc,
      integer   nr,
      integer   nc,
      integer   irow,
      integer   icol
    );

    //!
    //! Zeroes a rectangular block of the stored matrix
    //! staring at `(irow,icol)` position.
    //!
    //! \param[in] nr    number of rows of the block to be zeroed
    //! \param[in] nc    number of columns of the block to be zeroed
    //! \param[in] irow  starting row
    //! \param[in] icol  stating column
    //!
    void
    zero_block(
      integer nr,
      integer nc,
      integer irow,
      integer icol
    ) {
      gezero( nr, nc, m_data + this->iaddr(irow,icol), m_ldData );
    }

    //!
    //! Initialize the matrix as an identity matrix.
    //!
    void
    id( real_type dg ) {
      geid( m_nrows, m_ncols, m_data, m_ldData, dg );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void
    add_to_diag( real_type mu ) {
      axpy(
        std::min(m_nrows,m_ncols),
        1.0, &mu, 0,
        m_data, m_ldData+1
      );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //!
    //! Initialize matrix.
    //!
    //! \param[in] data   pointer of memory with data to be copied
    //! \param[in] ldData leading dimension of the memory to be copied
    //!
    //!
    void load( real_type const data[], integer ldData );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //!
    //! Initialize matrix.
    //!
    //! \param[in] data   pointer of memory with data to be copied
    //! \param[in] ldData leading dimension of the memory to be copied
    //!
    //!
    void load_transposed( real_type const data[], integer ldData );

    //!
    //! Initialize matrix using another matrix.
    //!
    //! \param[in] A initialize matrix with the matrix of the object A
    //!
    void load( MatW const & A );

    //!
    //! Initialize matrix using another matrix (transposed).
    //!
    //! \param[in] A initialize matrix with the matrix of the object A
    //!
    void load_transposed( MatW const & A );

    //!
    //! Copy a matrix to a rectangular block of the stored matrix
    //! staring at `(irow,icol)` position.
    //!
    //! \param[in] B     reference to the matrix of the object B
    //! \param[in] irow  starting row
    //! \param[in] icol  stating column
    //!
    void
    load_block(
      MatW const & B,
      integer      irow = 0,
      integer      icol = 0
    ) {
      this->load_block(
        B.m_nrows, B.m_ncols, B.m_data, B.m_ldData, irow, icol
      );
    }

    //!
    //! Copy a matrix to a rectangular block of the stored matrix
    //! staring at `(irow,icol)` position.
    //!
    //! \param[in] B     reference to the matrix of the object B
    //! \param[in] irow  starting row
    //! \param[in] icol  stating column
    //!
    void
    load_block_transposed(
      MatW const & B,
      integer      irow = 0,
      integer      icol = 0
    ) {
      this->load_block_transposed(
        B.m_nrows, B.m_ncols, B.m_data, B.m_ldData, irow, icol
      );
    }

    //!
    //! Copy a matrix to a rectangular block of the stored matrix
    //! staring at `(irow,icol)` position.
    //!
    //! \param[in] nr    number of rows of the block to be zeroed
    //! \param[in] nc    number of columns of the block to be zeroed
    //! \param[in] B     pointer to memory storing the input matrix `B`
    //! \param[in] ldB   leading dimension of the matrix `B`
    //! \param[in] irow  starting row
    //! \param[in] icol  stating column
    //!
    void
    load_block(
      integer         nr,
      integer         nc,
      real_type const B[],
      integer         ldB,
      integer         irow = 0,
      integer         icol = 0
    ) {
      integer inr = irow + nr;
      integer inc = icol + nc;
      UTILS_ASSERT(
        inr <= m_nrows &&
        inc <= m_ncols &&
        irow >= 0 && icol >= 0,
        "load_block( nr = {} nc = {}, B, irow = {}, icol = {}) bad parameters\n"
        "must be irow + nr ({}) <= nrows ({})\n"
        "must be icol + nc ({}) <= ncols ({})\n",
        nr, nc, irow, icol, inr, m_nrows, inc, m_ncols
      );
      integer info = gecopy(
        nr, nc,
        B, ldB,
        m_data + this->iaddr(irow,icol), m_ldData
      );
      UTILS_ASSERT(
        info == 0,
        "load_block call lapack_wrapper::gecopy return info = {}\n", info
      );
    }

    //!
    //! Copy a matrix to a rectangular block of the stored matrix
    //! staring at `(irow,icol)` position.
    //!
    //! \param[in] nr    number of rows of the block to be zeroed
    //! \param[in] nc    number of columns of the block to be zeroed
    //! \param[in] B     pointer to memory storing the input matrix `B`
    //! \param[in] ldB   leading dimension of the matrix `B`
    //! \param[in] irow  starting row
    //! \param[in] icol  stating column
    //!
    void
    load_block_transposed(
      integer         nr,
      integer         nc,
      real_type const B[],
      integer         ldB,
      integer         irow = 0,
      integer         icol = 0
    ) {
      UTILS_ASSERT(
        irow + nc <= m_nrows &&
        icol + nr <= m_ncols &&
        irow >= 0 && icol >= 0,
        "load_block_transpose( nr = {}, nc = {},..., irow = {}, icol = {} ) "
        "bad parameters\n",
        nr, nc, irow, icol
      );
      real_type const * pd = B;
      real_type       * pp = m_data + this->iaddr(irow,icol);
      for ( integer i = 0; i < nc; ++i, pd += ldB, ++pp )
        lapack_wrapper::copy( nr, pd, 1, pp, m_ldData );
    }

    //!
    //! Copy a matrix to a rectangular block of the stored matrix
    //! staring at `(irow,icol)` position
    //!
    //! \param[in] n     number of element of the diagonal
    //! \param[in] D     pointer to memory storing the diagonal
    //! \param[in] irow  starting row
    //! \param[in] icol  stating column
    //!
    void
    load_diagonal_block(
      integer         n,
      real_type const D[],
      integer         irow = 0,
      integer         icol = 0
    ) {
      UTILS_ASSERT(
        irow + n <= m_nrows &&
        icol + n <= m_ncols &&
        irow >= 0 && icol >= 0,
        "load_diagonal_block( n = {},..., irow = {}, icol = {}) "
        "bad parameters\n",
        n, irow, icol
      );
      for ( integer i = 0; i < n; ++i )
        m_data[this->iaddr(irow+i,icol+i)] = D[i];
    }

    //!
    //! Copy a matrix to a rectangular block of the stored matrix
    //! staring at `(irow,icol)` position.
    //!
    //! \param[in] D     reference to the matrix of the object D
    //! \param[in] irow  starting row
    //! \param[in] icol  stating column
    //!
    void
    load_block(
      DiagW const & D,
      integer       irow = 0,
      integer       icol = 0
    ) {
      this->load_diagonal_block( D.getDim(), D.data(), irow, icol );
    }

    //!
    //! Copy vector `column` to the `icol-th` column of
    //! the internal stored matrix.
    //!
    //! \param[in] column the column vector
    //! \param[in] icol   the column to be changed
    //!
    void
    load_column( real_type const column[], integer icol )
    { copy( m_nrows, column, 1, m_data + icol * m_ldData, 1 ); }

    //!
    //! Copy vector `row` to the `irow-th`
    //! row of the internal stored matrix.
    //!
    //! \param[in] row  the row vector
    //! \param[in] irow the row to be changed
    //!
    void
    load_row( real_type const row[], integer irow )
    { copy( m_ncols, row, 1, m_data + irow, m_ldData ); }

    void
    load_sparse_column(
      integer         nnz,
      real_type const values[],
      integer   const i_row[],
      integer         j
    );

    void
    load_sparse_row(
      integer         nnz,
      real_type const values[],
      integer         i,
      integer   const j_col[]
    );

    //!
    //! Copy sparse matrix into the object.
    //!
    //! \param[in] sp sparse matrix to be copied
    //!
    void load( Sparse const & sp );

    //!
    //! Copy sparse matrix into the object transposing it.
    //!
    //! \param[in] sp sparse matrix to be copied
    //!
    void load_transposed( Sparse const & sp );

    //!
    //! Copy sparse matrix into the object
    //!
    //! \param[in] rows row indices of the sparse matrix
    //! \param[in] cols column indices of the sparse matrix
    //! \param[in] vals values of the sparse matrix
    //! \param[in] nnz  number of nonzeros elements
    //!
    void
    load(
      integer   const rows[],
      integer   const cols[],
      real_type const vals[],
      integer         nnz
    );

    void
    load0(
      integer   const rows[],
      integer   const cols[],
      real_type const vals[],
      integer         nnz
    ) {
      this->zero_fill();
      this->load( rows, cols, vals, nnz );
    }

    //!
    //! Copy sparse matrix into the object.
    //!
    //! \param[in] rows row indices of the sparse matrix
    //! \param[in] cols column indices of the sparse matrix
    //! \param[in] vals values of the sparse matrix
    //! \param[in] nnz  number of nonzeros elements
    //!
    void
    load_symmetric(
      integer   const rows[],
      integer   const cols[],
      real_type const vals[],
      integer         nnz
    );

    void
    load0_symmetric(
      integer   const rows[],
      integer   const cols[],
      real_type const vals[],
      integer         nnz
    ) {
      this->zero_fill();
      this->load_symmetric( rows, cols, vals, nnz );
    }

    //!
    //! Copy sparse matrix into the object.
    //!
    //! \param sp     sparse matrix to be copied
    //! \param i_offs offset added to to the row indices
    //! \param j_offs offset added to to the column indices
    //!
    void
    load( Sparse const & sp, integer i_offs, integer j_offs );

    void
    load0( Sparse const & sp, integer i_offs, integer j_offs ) {
      this->zero_fill();
      this->load( sp, i_offs, j_offs );
    }

    //!
    //! Copy sparse matrix into the object transposing it
    //!
    //! \param sp     sparse matrix to be copied
    //! \param i_offs offset added to to the row indices
    //! \param j_offs offset added to to the column indices
    //!
    void
    load_transposed( Sparse const & sp, integer i_offs, integer j_offs );

    void
    load0_transposed( Sparse const & sp, integer i_offs, integer j_offs ) {
      this->zero_fill();
      this->load_transposed( sp, i_offs, j_offs );
    }

    //!
    //! Copy sparse matrix into the object
    //!
    //! \param[in] i_offs row index offset
    //! \param[in] j_offs column index offset
    //! \param[in] rows   row indices of the sparse matrix
    //! \param[in] cols   column indices of the sparse matrix
    //! \param[in] vals   values of the sparse matrix
    //! \param[in] nnz    number of nonzeros elements
    //!
    void
    load(
      integer         i_offs,
      integer         j_offs,
      integer   const rows[],
      integer   const cols[],
      real_type const vals[],
      integer         nnz
    );

    void
    load0(
      integer         i_offs,
      integer         j_offs,
      integer   const rows[],
      integer   const cols[],
      real_type const vals[],
      integer         nnz
    ) {
      this->zero_fill();
      this->load( i_offs, j_offs, rows, cols, vals, nnz );
    }

    //!
    //! Copy sparse matrix into the object
    //!
    //! \param[in] i_offs row index offset
    //! \param[in] j_offs column index offset
    //! \param[in] rows   row indices of the sparse matrix
    //! \param[in] cols   column indices of the sparse matrix
    //! \param[in] vals   values of the sparse matrix
    //! \param[in] nnz    number of nonzeros elements
    //!
    void
    load_symmetric(
      integer         i_offs,
      integer         j_offs,
      integer   const rows[],
      integer   const cols[],
      real_type const vals[],
      integer         nnz
    );

    void
    load0_symmetric(
      integer         i_offs,
      integer         j_offs,
      integer   const rows[],
      integer   const cols[],
      real_type const vals[],
      integer         nnz
    ) {
      this->zero_fill();
      this->load_symmetric( i_offs, j_offs, rows, cols, vals, nnz );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //!
    //! Add a block of memory to the matrix scaled by alpha.
    //!
    //! \param alpha  scaling factor
    //! \param data   pointer of memory with data to be copied
    //! \param ldData leading dimension of the memory to be copied
    //!
    void add( real_type alpha, real_type const data[], integer ldData );

    //!
    //! Add a block of memory to the matrix.
    //!
    //! \param data   pointer of memory with data to be copied
    //! \param ldData leading dimension of the memory to be copied
    //!
    void add( real_type const data[], integer ldData );

    //
    //! Compute `R <- R + alpha * A`
    //!
    //! \param alpha scale used in multiplication
    //! \param A     first mapped matrix
    //!
    void
    add( real_type alpha, MatW const & A ) {
      #if defined(DEBUG) || defined(_DEBUG)
      UTILS_ASSERT0(
        m_nrows == A.m_nrows && m_ncols == A.m_ncols,
        "add(..) incompatible dimensions\n"
      );
      #endif
      this->add( alpha, A.m_data, A.m_ldData );
    }

    void
    add( MatW const & A ) {
      #if defined(DEBUG) || defined(_DEBUG)
      UTILS_ASSERT0(
        m_nrows == A.m_nrows && m_ncols == A.m_ncols,
        "add(..) incompatible dimensions\n"
      );
      #endif
      this->add( A.m_data, A.m_ldData );
    }

    //!
    //! Add sparse matrix to the object
    //!
    //! \param sp sparse matrix to be copied
    //!
    void add( Sparse const & sp );

    //!
    //! Add sparse matrix to the object
    //!
    //! \param sp     sparse matrix to be copied
    //! \param i_offs offset added to to the row indices
    //! \param j_offs offset added to to the column indices
    //!
    void add( Sparse const & sp, integer i_offs, integer j_offs );

    //!
    //! Add sparse multiplied by `alpha` matrix to the object
    //!
    //! \param alpha  scalar used to multiply the sparse matrix
    //! \param sp     sparse matrix to be copied
    //!
    void add( real_type alpha, Sparse const & sp );

    //!
    //! Add sparse multiplied by `alpha` matrix to the object
    //!
    //! \param alpha  scalar used to multiply the sparse matrix
    //! \param sp     sparse matrix to be copied
    //! \param i_offs offset added to to the row indices
    //! \param j_offs offset added to to the column indices
    //!
    void
    add(
      real_type      alpha,
      Sparse const & sp,
      integer        i_offs,
      integer        j_offs
    );

    // alpha*A + beta*B -> C
    //!
    //! Add a linear combination of two matrix and assign to a third one.
    //!
    //! \param alpha scalar used to multiply the matrix `A`
    //! \param A     first mapped matrix
    //! \param beta  scalar used to multiply the matrix `B`
    //! \param B     second mapped matrix
    //! \param C     mapped matrix which store the result
    //!
    friend
    void
    add(
      real_type    alpha,
      MatW const & A,
      real_type    beta,
      MatW const & B,
      MatW       & C
    ) {
      geadd(
        C.m_nrows, C.m_ncols,
        alpha, A.m_data, A.m_ldData,
        beta,  B.m_data, B.m_ldData,
        C.m_data, C.m_ldData
      );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //!
    //! Make a view of a block.
    //!
    //! \param i_offs initial row of the block to be extracted
    //! \param j_offs initial column of the block to be extracted
    //! \param nrow   number of rows of the block to be extracted
    //! \param ncol   number of columns of the block to be extracted
    //! \param to     mapped matrix which store the result
    //!
    void
    view_block(
      integer i_offs,
      integer j_offs,
      integer nrow,
      integer ncol,
      MatW &  to
    ) {
      #if defined(DEBUG) || defined(_DEBUG)
      UTILS_ASSERT_TRACE(
        i_offs >= 0 && i_offs+nrow <= m_nrows &&
        j_offs >= 0 && j_offs+ncol <= m_ncols,
        "view_block(i_offs={}, j_offs={}, nrow={}, ncol={}) "
        "will be out of range [0,{}) x [0,{})\n",
        i_offs, j_offs, nrow, ncol, m_nrows, m_ncols
      );
      #endif
      integer offs = 0;
      if ( m_nrows > 0 && m_ncols > 0 &&
           nrow    > 0 && ncol    > 0 ) offs = this->iaddr(i_offs,j_offs);
      to.setup( m_data+offs, nrow, ncol, m_ldData );
    }

    //!
    //! Extract a matrix in transposed form
    //!
    //! \param out mapped matrix which store the result
    //!
    void
    get_transposed( MatW & out ) {
      #if defined(DEBUG) || defined(_DEBUG)
      UTILS_ASSERT0(
        m_nrows == out.m_ncols && m_ncols == out.m_nrows,
        "get_transposed(...) incompatible matrices\n"
      );
      #endif
      real_type const * pc = m_data;
      for ( integer i = 0; i < m_ncols; ++i, pc += m_ldData )
        lapack_wrapper::copy( m_nrows, pc, 1, out.m_data + i, out.m_ldData );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //!
    //! Extract a block
    //!
    //! \param i_offs initial row of the block to be extracted
    //! \param j_offs initial column of the block to be extracted
    //! \param to     mapped matrix which store the result
    //!
    void
    get_block(
      MatW &  to,
      integer i_offs,
      integer j_offs
    ) {
      #if defined(DEBUG) || defined(_DEBUG)
      UTILS_ASSERT_TRACE(
        i_offs >= 0 && i_offs+to.m_nrows <= m_nrows &&
        j_offs >= 0 && j_offs+to.m_ncols <= m_ncols,
        "get_block(i_offs={}, j_offs={}, ...) "
        "will be out of range [0,{}) x [0,{})\n",
        i_offs, j_offs, m_nrows, m_ncols
      );
      #endif
      lapack_wrapper::gecopy(
        to.m_nrows, to.m_ncols,
        m_data+this->iaddr(i_offs,j_offs), m_ldData,
        to.m_data, to.m_ldData
      );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //!
    //! Extract a block.
    //!
    //! \param i_offs initial row of the block to be extracted
    //! \param j_offs initial column of the block to be extracted
    //! \param to     mapped matrix which store the result
    //!
    void
    get_block_transposed(
      MatW &  to,
      integer i_offs,
      integer j_offs
    ) {
      #if defined(DEBUG) || defined(_DEBUG)
      UTILS_ASSERT_TRACE(
        i_offs >= 0 && i_offs+to.m_ncols <= m_nrows &&
        j_offs >= 0 && j_offs+to.m_nrows <= m_ncols,
        "get_block_transposed(i_offs={}, j_offs={}, ...) "
        "will be out of range [0,{}) x [0,{})\n",
        i_offs, j_offs, m_nrows, m_ncols
      );
      #endif
      real_type const * pc = m_data+this->iaddr(i_offs,j_offs);
      for ( integer i = 0; i < to.m_ncols; ++i, pc += m_ldData )
        lapack_wrapper::copy( to.m_nrows, pc, 1, to.m_data + i, to.m_ldData );
    }

    string
    to_string( real_type eps = 0 ) const {
      string res = "";
      for ( integer i = 0; i < m_nrows; ++i ) {
        for ( integer j = 0; j < m_ncols; ++j ) {
          real_type aij = (*this)(i,j);
          if ( std::abs( aij ) < eps ) res += ".              ";
          else                         res += fmt::format( "{:14} ", aij );
        }
        res += '\n';
      }
      return res;
    }

  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  template <typename T>
  inline
  ostream_type &
  operator << ( ostream_type & stream, MatrixWrapper<T> const & A ) {
    stream << A.to_string();
    return stream;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // c = beta*c + alpha*A*v
  //!
  //! Perform matrix vector multiplication `c = beta*c + alpha*A*v`.
  //!
  //! \param alpha  matrix `A` is multiplied by `alpha`
  //! \param TRANSA choose if multiply with the transpose of `A`
  //! \param A      matrix used in the multiplication
  //! \param v      vector to be multiplied
  //! \param incv   stride of the vector `v`
  //! \param beta   scalar used to multiply `c`
  //! \param c      result vector
  //! \param incc   stride of the vector `c`
  //!
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
      A.nrows(), A.ncols(),
      alpha,
      A.data(), A.ldim(),
      v, incv,
      beta,
      c, incc
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // c = beta*c + alpha*A*v
  //!
  //! Perform matrix vector multiplication `c = beta*c + alpha*A*v`.
  //!
  //! \param alpha  matrix `A` is multiplied by `alpha`
  //! \param A      matrix used in the multiplication
  //! \param v      vector to be multiplied
  //! \param incv   stride of the vector `v`
  //! \param beta   scalar used to multiply `c`
  //! \param c      result vector
  //! \param incc   stride of the vector `c`
  //!
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
      Transposition::NO_TRANSPOSE,
      A.nrows(), A.ncols(),
      alpha,
      A.data(), A.ldim(),
      v, incv,
      beta,
      c, incc
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //! Perform matrix matrix multiplication `C = beta*C + alpha*A*B`.
  //!
  //! \param where  message added in case of error
  //! \param alpha  matrix `A` is multiplied by `alpha`
  //! \param A      matrix used in the multiplication
  //! \param B      matrix used in the multiplication
  //! \param beta   scalar used to multiply `C`
  //! \param C      result matrix
  //!
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
    UTILS_ASSERT_DEBUG(
      A.ncols() == B.nrows() &&
      A.nrows() == C.nrows() &&
      B.ncols() == C.ncols(),
      "gemm, at `{}' inconsistent dimensions:\n"
      "A = {} x {}\nB = {} x {}\nC = {} x {}\n",
      where,
      A.nrows(), A.ncols(),
      B.nrows(), B.ncols(),
      C.nrows(), C.ncols()
    );
    lapack_wrapper::gemm(
      Transposition::NO_TRANSPOSE,
      Transposition::NO_TRANSPOSE,
      A.nrows(), B.ncols(), A.ncols(),
      alpha,
      A.data(), A.ldim(),
      B.data(), B.ldim(),
      beta,
      C.data(), C.ldim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //! Perform matrix matrix multiplication `C = beta*C + alpha*A*B`.
  //!
  //! \param alpha  matrix `A` is multiplied by `alpha`
  //! \param A      matrix used in the multiplication
  //! \param B      matrix used in the multiplication
  //! \param beta   scalar used to multiply `C`
  //! \param C      result matrix
  //!
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
    UTILS_ASSERT_DEBUG(
      A.ncols() == B.nrows() &&
      A.nrows() == C.nrows() &&
      B.ncols() == C.ncols(),
      "gemm, inconsistent dimensions:\n"
      "A = {} x {}\nB = {} x {}\nC = {} x {}\n",
      A.nrows(), A.ncols(), B.nrows(), B.ncols(), C.nrows(), C.ncols()
    );
    lapack_wrapper::gemm(
      Transposition::NO_TRANSPOSE,
      Transposition::NO_TRANSPOSE,
      A.nrows(), B.ncols(), A.ncols(),
      alpha,
      A.data(), A.ldim(),
      B.data(), B.ldim(),
      beta,
      C.data(), C.ldim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //! Perform matrix matrix multiplication `C = beta*C + alpha*A*B`.
  //!
  //! \param where  message added in case of error
  //! \param alpha  matrix `A` is multiplied by `alpha`
  //! \param TRANSA choose `A` transposed or not
  //! \param A      matrix used in the multiplication
  //! \param TRANSB choose `B` transposed or not
  //! \param B      matrix used in the multiplication
  //! \param beta   scalar used to multiply `C`
  //! \param C      result matrix
  //!
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
    integer Ar = TRANSA == Transposition::NO_TRANSPOSE ? A.nrows() : A.ncols();
    integer Ac = TRANSA == Transposition::NO_TRANSPOSE ? A.ncols() : A.nrows();
    integer Br = TRANSB == Transposition::NO_TRANSPOSE ? B.nrows() : B.ncols();
    integer Bc = TRANSB == Transposition::NO_TRANSPOSE ? B.ncols() : B.nrows();
    UTILS_ASSERT_DEBUG(
      C.nrows() == Ar && C.ncols() == Bc && Ac == Br,
      "gemm, at `{}' inconsistent dimensions:"
      "\nA = {} x {}\nB = {} x {}\nC = {} x {}"
      "\nA {}\nB {}\n",
      where,
      A.nrows(), A.ncols(),
      B.nrows(), B.ncols(),
      C.nrows(), C.ncols(),
      to_blas(TRANSA),
      to_blas(TRANSB)
    );
    lapack_wrapper::gemm(
      TRANSA,
      TRANSB,
      C.nrows(), C.ncols(), Ac,
      alpha,
      A.data(), A.ldim(),
      B.data(), B.ldim(),
      beta,
      C.data(), C.ldim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // C = beta*C + alpha*A*B (NO DEBUG)
  //!
  //! Perform matrix matrix multiplication `C = beta*C + alpha*A*B`.
  //!
  //! \param alpha  matrix `A` is multiplied by `alpha`
  //! \param TRANSA choose `A` transposed or not
  //! \param A      matrix used in the multiplication
  //! \param TRANSB choose `B` transposed or not
  //! \param B      matrix used in the multiplication
  //! \param beta   scalar used to multiply `C`
  //! \param C      result matrix
  //!
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
    integer Ar = TRANSA == Transposition::NO_TRANSPOSE ? A.nrows() : A.ncols();
    integer Ac = TRANSA == Transposition::NO_TRANSPOSE ? A.ncols() : A.nrows();
    integer Br = TRANSB == Transposition::NO_TRANSPOSE ? B.nrows() : B.ncols();
    integer Bc = TRANSB == Transposition::NO_TRANSPOSE ? B.ncols() : B.nrows();
    UTILS_ASSERT_DEBUG(
      C.nrows() == Ar && C.ncols() == Bc && Ac == Br,
      "gemm, inconsistent dimensions:"
      "\nA = {} x {}\nB = {} x {}\nC = {} x {}"
      "\nA {}\nB {}\n",
      A.nrows(), A.ncols(),
      B.nrows(), B.ncols(),
      C.nrows(), C.ncols(),
      to_blas(TRANSA),
      to_blas(TRANSB)
    );
    lapack_wrapper::gemm(
      TRANSA,
      TRANSB,
      C.nrows(), C.ncols(), Ac,
      alpha,
      A.data(), A.ldim(),
      B.data(), B.ldim(),
      beta,
      C.data(), C.ldim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // x = A * x or x = A^T * x
  //!
  //! Perform matrix vector multiplication `x = A * x` or `x = A^T * x`.
  //!
  //! \param[in]     UPLO   choose upper or lower part
  //! \param[in]     TRANS  choose `A` transposed or not
  //! \param[in]     DIAG   use or not diagonal elements to 1
  //! \param[in]     A      input matrix
  //! \param[in,out] x      vector used in the multiplcation and result
  //! \param[in]     incx   stride of vector `x`
  //!
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
    UTILS_ASSERT_DEBUG(
      A.nrows() == A.ncols(),
      "trmv, matrix is {} x {} expected square\n",
      A.nrows(), A.ncols()
    );
    lapack_wrapper::trmv(
      UPLO, TRANS, DIAG, A.nrows(), A.data(), A.ldim(), x, incx
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // x = A * x or x = A^T * x
  //!
  //! Perform matrix vector multiplication `x = A * x` or `x = A^T * x`.
  //!
  //! \param where  message added in case of error
  //! \param UPLO   choose upper or lower part
  //! \param TRANS  choose `A` transposed or not
  //! \param DIAG   use or not diagonal elements to 1
  //! \param A      matrix to multiplied
  //! \param x      vector used in the multiplcation and result
  //! \param incx   stride of vector `x`
  //!
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
    UTILS_ASSERT_DEBUG(
      A.nrows() == A.ncols(),
      "trmv at `{}` matrix is {} x {} expected square\n",
      where, A.nrows(), A.ncols()
    );
    lapack_wrapper::trmv(
      UPLO, TRANS, DIAG, A.nrows(), A.data(), A.ldim(), x, incx
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // x = A^(-1) * x or x = A^(-T) * x
  //!
  //! Perform matrix vector multiplication `x = A * x` or `x = A^T * x`.
  //!
  //! \param UPLO   choose upper or lower part
  //! \param TRANS  choose `A` transposed or not
  //! \param DIAG   use or not diagonal elements to 1
  //! \param A      matrix to multiplied
  //! \param x      vector used in the multiplcation and result
  //! \param incx   stride of vector `x`
  //!
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
    UTILS_ASSERT_DEBUG(
      A.nrows() == A.ncols(),
      "trsv, matrix is {} x {} expected square\n",
      A.nrows(), A.ncols()
    );
    lapack_wrapper::trsv(
      UPLO, TRANS, DIAG, A.nrows(), A.data(), A.ldim(), x, incx
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // x = A^(-1) * x or x = A^(-T) * x
  //!
  //! Perform matrix vector multiplication
  //! `x = A^(-1) * x` or `x = A^(-T) * x`.
  //!
  //! \param where  message added in case of error
  //! \param UPLO   choose upper or lower part
  //! \param TRANS  choose `A` transposed or not
  //! \param DIAG   use or not diagonal elements to 1
  //! \param A      matrix to multiplied
  //! \param x      vector used in the multiplcation and result
  //! \param incx   stride of vector `x`
  //!
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
    UTILS_ASSERT_DEBUG(
      A.nrows() == A.ncols(),
      "trsv at `{}` matrix is {} x {} expected square\n",
      where, A.nrows(), A.ncols()
    );
    lapack_wrapper::trsv(
      UPLO, TRANS, DIAG, A.nrows(), A.data(), A.ldim(), x, incx
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //! Perform matrix matrix multiplication
  //!
  //! \code{.unparsed}
  //!    B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
  //!    op( A ) = A   or   op( A ) = A**T.
  //! \endcode
  //!
  //! \param SIDE   multiply on left or right
  //! \param UPLO   choose upper or lower part
  //! \param TRANS  choose `A` transposed or not
  //! \param DIAG   use or not diagonal elements to 1
  //! \param alpha  scalar used in the multiplication
  //! \param A      matrix to be multiplied
  //! \param B      matrix to be multiplied
  //!
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
    UTILS_ASSERT_DEBUG(
      A.nrows() == A.ncols(),
      "trmm, matrix is {} x {} expected square\n",
      A.nrows(), A.ncols()
    );
    lapack_wrapper::trmm(
      SIDE, UPLO, TRANS, DIAG,
      B.nrows(), B.ncols(), alpha,
      A.data(), A.ldim(),
      B.data(), B.ldim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //! Perform matrix matrix multiplication
  //!
  //! \code{.unparsed}
  //!    B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
  //!    op( A ) = A   or   op( A ) = A**T.
  //! \endcode
  //!
  //! \param where  message added in case of error
  //! \param SIDE   multiply on left or right
  //! \param UPLO   choose upper or lower part
  //! \param TRANS  choose `A` transposed or not
  //! \param DIAG   use or not diagonal elements to 1
  //! \param alpha  scalar used in the multiplication
  //! \param A      matrix to be multiplied
  //! \param B      matrix to be multiplied
  //!
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
    UTILS_ASSERT_DEBUG(
      A.nrows() == A.ncols(),
      "trmm, at `{}` matrix is {} x {} expected square\n",
      where, A.nrows(), A.ncols()
    );
    lapack_wrapper::trmm(
      SIDE, UPLO, TRANS, DIAG,
      B.nrows(), B.ncols(), alpha,
      A.data(), A.ldim(),
      B.data(), B.ldim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //! Perform matrix matrix multiplication
  //!
  //! \code{.unparsed}
  //!    B := alpha*op( A^(-1) )*B,   or   B := alpha*B*op( A^(-1) ),
  //!    op( A ) = A   or   op( A ) = A**T.
  //! \endcode
  //!
  //! \param SIDE   multiply on left or right
  //! \param UPLO   choose upper or lower part
  //! \param TRANS  choose `A` transposed or not
  //! \param DIAG   use or not diagonal elements to 1
  //! \param alpha  scalar used in the multiplication
  //! \param A      matrix to be multiplied
  //! \param B      matrix to be multiplied
  //!
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
    UTILS_ASSERT_DEBUG(
      A.nrows() == A.ncols(),
      "trsm, matrix is {} x {} expected square\n",
      A.nrows(), A.ncols()
    );
    lapack_wrapper::trsm(
      SIDE, UPLO, TRANS, DIAG,
      B.nrows(), B.ncols(), alpha,
      A.data(), A.ldim(),
      B.data(), B.ldim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //! Perform matrix matrix multiplication.
  //!
  //! \code{.unparsed}
  //!    B := alpha*op( A^(-1) )*B,   or   B := alpha*B*op( A^(-1) ),
  //!    op( A ) = A   or   op( A ) = A**T.
  //! \endcode
  //!
  //! \param where  message added in case of error
  //! \param SIDE   multiply on left or right
  //! \param UPLO   choose upper or lower part
  //! \param TRANS  choose `A` transposed or not
  //! \param DIAG   use or not diagonal elements to 1
  //! \param alpha  scalar used in the multiplication
  //! \param A      matrix to be multiplied
  //! \param B      matrix to be multiplied
  //!
  //!
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
    UTILS_ASSERT_DEBUG(
      A.nrows() == A.ncols(),
      "trmm, at `{}` matrix is {} x {} expected square\n",
      where, A.nrows(), A.ncols()
    );
    lapack_wrapper::trsm(
      SIDE, UPLO, TRANS, DIAG,
      B.nrows(), B.ncols(), alpha,
      A.data(), A.ldim(),
      B.data(), B.ldim()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  inline
  T
  normInf( MatrixWrapper<T> const & A ) {
    return lapack_wrapper::normInf( A.nrows(), A.ncols(), A.data(), A.ldim() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  inline
  T
  norm1( MatrixWrapper<T> const & A ) {
    return lapack_wrapper::norm1( A.nrows(), A.ncols(), A.data(), A.ldim() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  inline
  T
  normF( MatrixWrapper<T> const & A ) {
    return lapack_wrapper::normF( A.nrows(), A.ncols(), A.data(), A.ldim() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  inline
  T
  maxabs( MatrixWrapper<T> const & A ) {
    return lapack_wrapper::maxabs( A.nrows(), A.ncols(), A.data(), A.ldim() );
  }

  // explicit instantiation declaration to suppress warnings
  extern template class MatrixWrapper<real>;
  extern template class MatrixWrapper<doublereal>;

}

///
/// eof: wrapper.hxx
///
