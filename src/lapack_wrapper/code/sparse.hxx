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
/// file: sparse.hxx
///

#include <vector>
#include <string>

namespace lapack_wrapper {

  using std::vector;
  using std::string;

  template <typename T> class MatrixWrapper;

  template <typename T>
  class SparseMatrixBase {
  public:
    using real_type = T;
    using Sparse    = SparseMatrixBase<T>;
    using MatW      = MatrixWrapper<T>;

  protected:
    integer m_nrows; //!< Number of rows
    integer m_ncols; //!< Number of columns
    integer m_nnz;   //!< Total number of nonzeros

    //!
    //! Protected Constructor of the class SparseMatrixBase.
    //!
    SparseMatrixBase()
    : m_nrows(0)
    , m_ncols(0)
    , m_nnz(0)
    {}

    SparseMatrixBase( integer N, integer M )
    : m_nrows(N)
    , m_ncols(M)
    , m_nnz(0)
    {}

    void
    y_manage(
      real_type beta,
      integer   DimY,
      real_type y[],
      integer   incY
    ) const {
      if ( Utils::is_zero(beta) ) {
        for ( integer i = 0; i < DimY; ++i ) y[i*incY] = 0;
      } else if ( !Utils::is_zero(beta-1) ) {
        for ( integer i = 0; i < DimY; ++i ) y[i*incY] *= beta;
      }
    }

  public:

    //!
    //! Virtual Destructor of the class SparseMatrixBase.
    //!
    virtual
    ~SparseMatrixBase() {}

    //!
    //! Initialize sparse matrix.
    //!
    //! \param[in] N           number of rows
    //! \param[in] M           number of columns
    //! \param[in] reserve_nnz estimated number of nonzeros
    //! \param[in] fi          if true use FORTRAN 1-based indices
    //!
    virtual
    void
    init(
      integer N,
      integer M,
      integer reserve_nnz,
      bool    fi = false
    ) = 0;

    //!
    //! Destroy the sparse matrix.
    //!
    virtual
    void
    clear() = 0;

    //!
    //! Clean to zero the sparse matrix.
    //!
    virtual
    void
    setZero() = 0;

    //!
    //! Return true if index are 1-based.
    //!
    virtual
    bool
    FORTRAN_indexing() const = 0;

    integer get_number_of_rows() const { return m_nrows; }
    integer get_number_of_cols() const { return m_ncols; }
    integer get_nnz()            const { return m_nnz; }

    //!
    //! Returns the number of nonzeroes and dimension of the sparse matrix.
    //!
    //! \param[out] nrows Row dimension of the matrix
    //! \param[out] ncols Column dimension of the matrix
    //! \param[out] nnz     the number of nonzeroes
    //!
    virtual
    void
    get_info(
      integer & nrows,
      integer & ncols,
      integer & nnz
    ) const {
      nrows = m_nrows;
      ncols = m_ncols;
      nnz   = m_nnz;
    }

    //!
    //! Returns the sparse matrix.
    //!
    //! **Example of usage**
    //!
    //! \code
    //! int_type nnz;
    //! int_type const * rows;
    //! int_type const * cols;
    //! real     const * values;
    //! ptr->get_data( nnz, &rows, &cols, &values );
    //! \endcode
    //!
    //! \param[out] pRows   vector of pointers of row indices where to data will be copied.
    //! \param[out] pCols   vector of pointers where to store column indices.
    //! \param[out] pValues vector of pointers where to store values.
    //!
    virtual
    void
    get_data(
      integer   const * & pRows,
      integer   const * & pCols,
      real_type const * & pValues
    ) const = 0;

    //!
    //! Transpose the matrix.
    //!
    virtual
    void
    transpose() = 0;

    //!
    //! Convert the matrix to 1-based index.
    //!
    virtual
    void
    to_FORTRAN_indexing() = 0;

    //!
    //! Convert the matrix to 0-based index.
    //!
    virtual
    void
    to_C_indexing() = 0;

    //!
    //! Insert a value into matrix using 0-based index.
    //!
    //! \param[in] row the row index
    //! \param[in] col the column index
    //! \param[in] val the inserted value
    //!
    virtual
    void
    push_value_C(
      integer   row,
      integer   col,
      real_type val
    ) = 0;

    //!
    //! Insert a value into matrix using 1-based index.
    //!
    //! \param[in] row the row index
    //! \param[in] col the column index
    //! \param[in] val the inserted value
    //!
    virtual
    void
    push_value_F(
      integer   row,
      integer   col,
      real_type val
    ) = 0;

    //!
    //! Insert a matrix into the sparse matrix.
    //!
    //! \param[in] row_offs    offset to the row indices
    //! \param[in] col_offs    offset to the column indices
    //! \param[in] Matrix      the (full) matrix to be inserted
    //! \param[in] transpose   if true matrix is inserted transposed
    //! \param[in] lower_upper 0 = full matrix
    //!                        1 = upper part with diagonal
    //!                        2 = upper part without diagonal
    //!                        -1 = lower part with diagonal
    //!                        -2 = lower part without diagonal
    //!
    virtual
    void
    push_matrix(
      integer      row_offs,
      integer      col_offs,
      MatW const & Matrix,
      bool         transpose   = false,
      integer      lower_upper = 0
    );

    //!
    //! Insert a matrix into the sparse matrix.
    //!
    //! \param[in] row_offs  offset to the row indices
    //! \param[in] col_offs  offset to the column indices
    //! \param[in] Matrix    the (sparsse) matrix to be inserted
    //! \param[in] transpose if true matrix is inserted transposed
    //! \param[in] lower_upper 0 = full matrix
    //!                        1 = upper part with diagonal
    //!                        2 = upper part without diagonal
    //!                        3 = push symmetric, for each A(i,j) element an A(j,i) is inserted
    //!                       -3 = push anty symmetric, for each A(i,j) element an -A(j,i) is inserted
    //!                       -1 = lower part with diagonal
    //!                       -2 = lower part without diagonal
    //!
    virtual
    void
    push_matrix(
      integer        row_offs,
      integer        col_offs,
      Sparse const & Matrix,
      bool           transpose   = false,
      integer        lower_upper = 0
    );

    //!
    //! Extract the sparse matrix to a full matrix.
    //!
    //! \param[out] M the output matrix
    //!
    virtual
    void
    get_matrix( MatW & M ) const = 0;

    //!
    //! Extract the sparse matrix to a full matrix
    //! filling symmetrically.
    //!
    //! \param[out] M the output matrix
    //!
    virtual
    void
    get_matrix_symmetric( MatW & M ) const = 0;

    //!
    //! Extract the sparse matrix to a full matrix
    //! trasposing the result.
    //!
    //! \param[out] M the output matrix
    //!
    virtual
    void
    get_matrix_transposed( MatW & M ) const = 0;

    //!
    //! \return `true` if Inf or NaN are found in the sparse data
    //!
    virtual
    bool
    foundNaN() const = 0;

    //!
    //! Calls the blas-Routine (\f$y = \beta y + \alpha A x + y\f$).
    //!
    //! \param[in]     alpha Scalar \f$\alpha \f$.
    //! \param[in]     DimX  Dimension of the vector x.
    //! \param[in]     x     Vector x.
    //! \param[in]     incX  Stride of vector x.
    //! \param[in]     DimY  Dimension of the vector y.
    //! \param[in]     beta  Scalar \f$\beta \f$.
    //! \param[in,out] y     In-/Output vector y.
    //! \param[in]     incY  Stride of vector y.
    //!
    //!
    virtual
    void
    gemv(
      real_type       alpha,
      integer         DimX,
      real_type const x[],
      integer         incX,
      real_type       beta,
      integer         DimY,
      real_type       y[],
      integer         incY
    ) const = 0;

    //!
    //! Calls the blas-Routine (\f$y = \beta y + \alpha A^T x + y\f$).
    //!
    //! \param[in]     alpha Scalar \f$\alpha \f$.
    //! \param[in]     DimX  Dimension of the vector x.
    //! \param[in]     x     Vector x.
    //! \param[in]     incX  Stride of vector x.
    //! \param[in]     DimY  Dimension of the vector y.
    //! \param[in]     beta  Scalar \f$\beta \f$.
    //! \param[in,out] y     In-/Output vector y.
    //! \param[in]     incY  Stride of vector y.
    //!
    virtual
    void
    gemv_Transposed(
      real_type       alpha,
      integer         DimX,
      real_type const x[],
      integer         incX,
      real_type       beta,
      integer         DimY,
      real_type       y[],
      integer         incY
    ) const = 0;

    //!
    //! Calls the blas-Routine (\f$y = \beta y + \alpha (A+A^T-\textrm{diag}(A)) x + y\f$).
    //!
    //! \param[in]     alpha Scalar \f$\alpha \f$.
    //! \param[in]     DimX  Dimension of the vector x.
    //! \param[in]     x     Vector x.
    //! \param[in]     incX  Stride of vector x.
    //! \param[in]     DimY  Dimension of the vector y.
    //! \param[in]     beta  Scalar \f$\beta \f$.
    //! \param[in,out] y     In-/Output vector y.
    //! \param[in]     incY  Stride of vector y.
    //!
    virtual
    void
    gemv_Symmetric(
      real_type       alpha,
      integer         DimX,
      real_type const x[],
      integer         incX,
      real_type       beta,
      integer         DimY,
      real_type       y[],
      integer         incY
    ) const = 0;

    virtual
    void
    get_full_view( MatrixWrapper<real_type> & ) {
      UTILS_ERROR0( "get_full_view not defined\n");
    }

    //!
    //! Print the sparse matrix to a stream element by element:
    //!
    //! row, col, value
    //!
    //! \param[in] stream  Stream where to write elements.
    //!
    string to_string() const;

  };

  /*
  //   ____                             ____ ____ ___   ___  ____
  //  / ___| _ __   __ _ _ __ ___  ___ / ___/ ___/ _ \ / _ \|  _ \
  //  \___ \| '_ \ / _` | '__/ __|/ _ \ |  | |  | | | | | | | |_) |
  //   ___) | |_) | (_| | |  \__ \  __/ |__| |__| |_| | |_| |  _ <
  //  |____/| .__/ \__,_|_|  |___/\___|\____\____\___/ \___/|_| \_\
  //        |_|
  */

  //! Sparse Matrix Structure
  template <typename real>
  class SparseCCOOR : public SparseMatrixBase<real> {
  public:
    using Sparse    = SparseMatrixBase<real>;
    using MatW      = MatrixWrapper<real>;
    using real_type = typename Sparse::real_type;

  protected:
    std::vector<real_type> m_vals;  //!< the values of the sparse matrix
    std::vector<integer>   m_rows;  //!< the rows index
    std::vector<integer>   m_cols;  //!< the columns index

    bool m_fortran_indexing;
    bool m_matrix_is_full;
    bool m_matrix_is_row_major;

  public:

    using SparseMatrixBase<real>::m_nrows;
    using SparseMatrixBase<real>::m_ncols;
    using SparseMatrixBase<real>::m_nnz;

    SparseCCOOR()
    : SparseMatrixBase<real>()
    , m_fortran_indexing(false)
    , m_matrix_is_full(false)
    , m_matrix_is_row_major(false)
    {}

    SparseCCOOR(
      integer N,
      integer M,
      integer reserve_nnz,
      bool    fi
    );

    ~SparseCCOOR() override
    {}

    void clear() override;
    void setZero() override;

    void
    init(
      integer N,
      integer M,
      integer reserve_nnz,
      bool    fi = false
    ) override;

    bool
    FORTRAN_indexing() const override
    { return m_fortran_indexing; }

    void
    transpose() override {
      m_rows.swap(m_cols);
      std::swap( m_nrows, m_ncols );
    }

    void to_FORTRAN_indexing() override;
    void to_C_indexing() override;

    void
    push_value_C(
      integer   row,
      integer   col,
      real_type val
    ) override;

    void
    push_value_F(
      integer   row,
      integer   col,
      real_type val
    ) override;

    void get_matrix( MatW & M ) const override;
    void get_matrix_symmetric( MatW & M ) const override;
    void get_matrix_transposed( MatW & M ) const override;
    bool foundNaN() const override;

    void
    get_data(
      integer   const * & pRows,
      integer   const * & pCols,
      real_type const * & pValues
    ) const override {
      pRows   = &m_rows.front();
      pCols   = &m_cols.front();
      pValues = &m_vals.front();
    }

    void
    gemv(
      real_type       alpha,
      integer         DimX,
      real_type const x[],
      integer         incX,
      real_type       beta,
      integer         DimY,
      real_type       y[],
      integer         incY
    ) const override;

    void
    gemv_Transposed(
      real_type       alpha,
      integer         DimX,
      real_type const x[],
      integer         incX,
      real_type       beta,
      integer         DimY,
      real_type       y[],
      integer         incY
    ) const override;

    void
    gemv_Symmetric(
      real_type       alpha,
      integer         DimX,
      real_type const x[],
      integer         incX,
      real_type       beta,
      integer         DimY,
      real_type       y[],
      integer         incY
    ) const override;

    /*\
    :|:           _    _ _ _   _               _             _   _            _
    :|:   __ _ __| |__| (_) |_(_)___ _ _  __ _| |  _ __  ___| |_| |_  ___  __| |___
    :|:  / _` / _` / _` | |  _| / _ \ ' \/ _` | | | '  \/ -_)  _| ' \/ _ \/ _` (_-<
    :|:  \__,_\__,_\__,_|_|\__|_\___/_||_\__,_|_| |_|_|_\___|\__|_||_\___/\__,_/__/
    \*/

    template <typename i_type, typename r_type>
    void
    export_data(
      integer NNZ,
      i_type  i[],
      i_type  j[],
      r_type  val[],
      bool    fi
    ) {
      UTILS_ASSERT(
        NNZ == m_rows.size() &&
        NNZ == m_cols.size() &&
        NNZ == m_vals.size(),
        "export_data, bad dimension\n"
        "NNZ = {}, |rows| = {}, |cols| = {}, |vals| = {}\n",
        NNZ, m_rows.size(), m_cols.size(), m_vals.size()
      );
      integer offs{0};
      if ( fi ) ++offs;
      if ( m_fortran_indexing ) --offs;
      for ( integer index = 0; index < NNZ; ++index ) {
        i[index]   = i_type(m_rows[index]+offs);
        j[index]   = i_type(m_cols[index]+offs);
        val[index] = r_type(m_vals[index]);
      }
    }

    template <typename i_type, typename r_type>
    void
    export_data(
      std::vector<i_type> & i,
      std::vector<i_type> & j,
      std::vector<r_type> & v,
      bool                  fi
    ) {
      integer offs{0};
      if ( fi ) ++offs;
      if ( m_fortran_indexing ) --offs;
      i.clear(); i.reserve( m_nnz );
      j.clear(); j.reserve( m_nnz );
      v.clear(); v.reserve( m_nnz );
      for ( integer index = 0; index < m_nnz; ++index ) {
        i.push_back( i_type(m_rows[index]+offs) );
        j.push_back( i_type(m_cols[index]+offs) );
        v.push_back( r_type(m_vals[index]) );
      }
    }

    void
    setup_as_full_row_major(
      integer N,
      integer M,
      bool    fi = false
    );

    void
    setup_as_full_column_major(
      integer N,
      integer M,
      bool    fi = false
    );

    void
    fill( real_type const V[], integer M );

    void
    fill( std::vector<real_type> const & V );

    void
    reserve( integer reserve_nnz );

    void
    get_full_view( MatrixWrapper<real_type> & MW ) override {
      UTILS_ASSERT0(
        m_matrix_is_full,
        "get_full_view, matrix is sparse\n"
      );
      if ( m_matrix_is_row_major ) {
        MW.setup( &m_vals.front(), m_ncols, m_nrows, m_ncols );
      } else {
        MW.setup( &m_vals.front(), m_nrows, m_ncols, m_nrows );
      }
    }
  };

  // explicit instantiation declaration to suppress warnings
  extern template class SparseCCOOR<real>;
  extern template class SparseCCOOR<doublereal>;

}

///
/// eof: sparse.hxx
///
