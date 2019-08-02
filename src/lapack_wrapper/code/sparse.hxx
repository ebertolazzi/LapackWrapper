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

namespace lapack_wrapper {

  template <typename T> class MatrixWrapper;

  template <typename T>
  class SparseMatrixBase {
  public:
    typedef T                   valueType;
    typedef SparseMatrixBase<T> Sparse;
    typedef MatrixWrapper<T>    MatW;

  protected:
    integer nRows; //!< Number of rows
    integer nCols; //!< Number of columns
    integer nnz;   //!< Total number of nonzeros

    /*!
     * \brief SparseMatrixBase:
     *        Protected Constructor of the class SparseMatrixBase.
    \*/
    SparseMatrixBase()
    : nRows(0)
    , nCols(0)
    , nnz(0)
    {}

    SparseMatrixBase( integer N, integer M )
    : nRows(N)
    , nCols(M)
    , nnz(0)
    {}

    void
    y_manage(
      valueType beta,
      integer   DimY,
      valueType y[],
      integer   incY
    ) const {
      if ( isZero(beta) ) {
        for ( integer i = 0; i < DimY; ++i ) y[i*incY] = 0;
      } else if ( !isZero(beta-1) ) {
        for ( integer i = 0; i < DimY; ++i ) y[i*incY] *= beta;
      }
    }

  public:

    /*!
     * \brief ~SparseMatrixBase:
     *        Virtual Destructor of the class SparseMatrixBase.
    \*/
    virtual
    ~SparseMatrixBase() {}

    /*!
     *  \brief Initialize sparse matrix
     *
     *  \param[in] N           number of rows
     *  \param[in] M           number of columns
     *  \param[in] reserve_nnz estimated number of nonzeros
     *  \param[in] fi          if true use FORTRAN 1-based indices
    \*/
    virtual
    void
    init(
      integer N,
      integer M,
      integer reserve_nnz,
      bool    fi = false
    ) LAPACK_WRAPPER_PURE_VIRTUAL;

    /*!
     *  \brief Destroy the sparse matrix
     *
    \*/
    virtual
    void
    clear() LAPACK_WRAPPER_PURE_VIRTUAL;

    /*!
     *  \brief Clean to zero the sparse matrix
     *
    \*/
    virtual
    void
    setZero() LAPACK_WRAPPER_PURE_VIRTUAL;

    /*!
     *  \brief Return true if index are 1-based
    \*/
    virtual
    bool
    FORTRAN_indexing() const LAPACK_WRAPPER_PURE_VIRTUAL;

    integer
    get_number_of_rows() const
    { return this->nRows; }

    integer
    get_number_of_cols() const
    { return this->nCols; }

    integer
    get_nnz() const
    { return this->nnz; }

    /*!
     * \brief get_info:
     *        Returns the number of nonzeroes and dimension of the sparse matrix.
     * \param[out] _numRows Row dimension of the matrix
     * \param[out] _numCols Column dimension of the matrix
     * \param[out] _nnz     the number of nonzeroes
     *
    \*/
    virtual
    void
    get_info(
      integer & _numRows,
      integer & _numCols,
      integer & _nnz
    ) const {
      _numRows = this->nRows;
      _numCols = this->nCols;
      _nnz     = this->nnz;
    }

    /*!
     * \brief get_data:
     *        Returns the sparse matrix.
     *
     * Example of usage
     *
     * \code
     * int_type nnz;
     * int_type const * rows;
     * int_type const * cols;
     * real     const * values;
     * ptr->get_data( nnz, &rows, &cols, &values );
     * \endcode
     *
     * \param[out] pRows   vector of pointers of row indices where to data will be copied.
     * \param[out] pCols   vector of pointers where to store column indices.
     * \param[out] pValues vector of pointers where to store values.
     *
    \*/

    virtual
    void
    get_data(
      integer   const * & pRows,
      integer   const * & pCols,
      valueType const * & pValues
    ) const LAPACK_WRAPPER_PURE_VIRTUAL;

    /*!
     *  \brief transpose the matrix
    \*/
    virtual
    void
    transpose() LAPACK_WRAPPER_PURE_VIRTUAL;

    /*!
     *  \brief convert the matrix to 1-based index
    \*/
    virtual
    void
    to_FORTRAN_indexing() LAPACK_WRAPPER_PURE_VIRTUAL;

    /*!
     *  \brief convert the matrix to 0-based index
    \*/
    virtual
    void
    to_C_indexing() LAPACK_WRAPPER_PURE_VIRTUAL;

    /*!
     *  \brief insert a value into matrix using 0-based index
     *
     *  \param[in] row the row index
     *  \param[in] col the column index
     *  \param[in] val the inserted value
     *
    \*/
    virtual
    void
    push_value_C(
      integer   row,
      integer   col,
      valueType val
    ) LAPACK_WRAPPER_PURE_VIRTUAL;

    /*!
     *  \brief insert a value into matrix using 1-based index
     *
     *  \param[in] row the row index
     *  \param[in] col the column index
     *  \param[in] val the inserted value
     *
    \*/
    virtual
    void
    push_value_F(
      integer   row,
      integer   col,
      valueType val
    ) LAPACK_WRAPPER_PURE_VIRTUAL;

    /*!
     *  \brief insert a matrix into the sparse matrix
     *
     *  \param[in] row_offs    offset to the row indices
     *  \param[in] col_offs    offset to the column indices
     *  \param[in] Matrix      the (full) matrix to be inserted
     *  \param[in] transpose   if true matrix is inserted transposed
     *  \param[in] lower_upper 0 = full matrix
     *                         1 = upper part with diagonal
     *                         2 = upper part without diagonal
     *                         -1 = lower part with diagonal
     *                         -2 = lower part without diagonal
     *
    \*/
    virtual
    void
    push_matrix(
      integer      row_offs,
      integer      col_offs,
      MatW const & Matrix,
      bool         transpose   = false,
      integer      lower_upper = 0
    );

    /*!
     *  \brief insert a matrix into the sparse matrix
     *
     *  \param[in] row_offs  offset to the row indices
     *  \param[in] col_offs  offset to the column indices
     *  \param[in] Matrix    the (sparsse) matrix to be inserted
     *  \param[in] transpose if true matrix is inserted transposed
     *  \param[in] lower_upper 0 = full matrix
     *                         1 = upper part with diagonal
     *                         2 = upper part without diagonal
     *                         3 = push symmetric, for each A(i,j) element an A(j,i) is inserted
     *                        -3 = push anty symmetric, for each A(i,j) element an -A(j,i) is inserted
     *                        -1 = lower part with diagonal
     *                        -2 = lower part without diagonal
     *
    \*/
    virtual
    void
    push_matrix(
      integer        row_offs,
      integer        col_offs,
      Sparse const & Matrix,
      bool           transpose   = false,
      integer        lower_upper = 0
    );

    /*!
     *  \brief extract the sparse matrix to a full matrix
     *
     *  \param[out] M the output matrix
     *
    \*/
    virtual
    void
    get_matrix( MatW & M ) const LAPACK_WRAPPER_PURE_VIRTUAL;

    /*!
     *  \brief extract the sparse matrix to a full matrix
     *         filling symmetrically
     *
     *  \param[out] M the output matrix
     *
    \*/
    virtual
    void
    get_matrix_symmetric( MatW & M ) const LAPACK_WRAPPER_PURE_VIRTUAL;

    /*!
     *  \brief extract the sparse matrix to a full matrix
     *         trasposing the result
     *
     *  \param[out] M the output matrix
     *
    \*/
    virtual
    void
    get_matrix_transposed( MatW & M ) const LAPACK_WRAPPER_PURE_VIRTUAL;

    /*!
     *  \return true if Inf or NaN are found in the sparse data
    \*/
    virtual
    bool
    foundNaN() const LAPACK_WRAPPER_PURE_VIRTUAL;

    /*!
     * \brief gemv:
     *        Calls the blas-Routine (\f$y = \beta y + \alpha A x + y\f$).
     *
     * \param[in]     alpha         Scalar \f$\alpha \f$.
     * \param[in]     DimX          Dimension of the vector x.
     * \param[in]     x             Vector x.
     * \param[in]     incX          Stride of vector x.
     * \param[in]     DimY          Dimension of the vector y.
     * \param[in]     beta          Scalar \f$\beta \f$.
     * \param[in,out] y             In-/Output vector y.
     * \param[in]     incY          Stride of vector y.
     *
    \*/
    virtual
    void
    gemv(
      valueType       alpha,
      integer         DimX,
      valueType const x[],
      integer         incX,
      valueType       beta,
      integer         DimY,
      valueType       y[],
      integer         incY
    ) const LAPACK_WRAPPER_PURE_VIRTUAL;

    /*!
     * \brief gemv:
     *        Calls the blas-Routine (\f$y = \beta y + \alpha A^T x + y\f$).
     *
     * \param[in]     alpha Scalar \f$\alpha \f$.
     * \param[in]     DimX  Dimension of the vector x.
     * \param[in]     x     Vector x.
     * \param[in]     incX  Stride of vector x.
     * \param[in]     DimY  Dimension of the vector y.
     * \param[in]     beta  Scalar \f$\beta \f$.
     * \param[in,out] y     In-/Output vector y.
     * \param[in]     incY  Stride of vector y.
     *
    \*/

    virtual
    void
    gemv_Transposed(
      valueType       alpha,
      integer         DimX,
      valueType const x[],
      integer         incX,
      valueType       beta,
      integer         DimY,
      valueType       y[],
      integer         incY
    ) const LAPACK_WRAPPER_PURE_VIRTUAL;

    /*!
     * \brief gemv:
     *        Calls the blas-Routine (\f$y = \beta y + \alpha (A+A^T-\textrm{diag}(A)) x + y\f$).
     *
     * \param[in]     alpha Scalar \f$\alpha \f$.
     * \param[in]     DimX  Dimension of the vector x.
     * \param[in]     x     Vector x.
     * \param[in]     incX  Stride of vector x.
     * \param[in]     DimY  Dimension of the vector y.
     * \param[in]     beta  Scalar \f$\beta \f$.
     * \param[in,out] y     In-/Output vector y.
     * \param[in]     incY  Stride of vector y.
     *
    \*/
    virtual
    void
    gemv_Symmetric(
      valueType       alpha,
      integer         DimX,
      valueType const x[],
      integer         incX,
      valueType       beta,
      integer         DimY,
      valueType       y[],
      integer         incY
    ) const LAPACK_WRAPPER_PURE_VIRTUAL;

    virtual
    void
    get_full_view( MatrixWrapper<valueType> & ) {
      LAPACK_WRAPPER_ERROR( "get_full_view not defined");
    }

    /*!
     * \brief print:
     *        Print the sparse matrix to a stream element by element
     *
     *        row, col, value
     *
     *
     * \param[in] stream  Stream where to write elements.
     *
     */
    void print( ostream_type & stream ) const;

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
    typedef SparseMatrixBase<real>     Sparse;
    typedef MatrixWrapper<real>        MatW;
    typedef typename Sparse::valueType valueType;

  protected:
    std::vector<valueType> vals;  //!< the values of the sparse matrix
    std::vector<integer>   rows;  //!< the rows index
    std::vector<integer>   cols;  //!< the columns index

    bool fortran_indexing;
    bool matrix_is_full;
    bool matrix_is_row_major;

  public:

    using SparseMatrixBase<real>::nRows;
    using SparseMatrixBase<real>::nCols;
    using SparseMatrixBase<real>::nnz;

    SparseCCOOR()
    : SparseMatrixBase<real>()
    , fortran_indexing(false)
    , matrix_is_full(false)
    , matrix_is_row_major(false)
    {}

    SparseCCOOR(
      integer N,
      integer M,
      integer reserve_nnz,
      bool    fi
    );

    virtual
    ~SparseCCOOR() LAPACK_WRAPPER_OVERRIDE
    {}

    virtual
    void
    clear() LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    setZero() LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    init(
      integer N,
      integer M,
      integer reserve_nnz,
      bool    fi = false
    ) LAPACK_WRAPPER_OVERRIDE;

    virtual
    bool
    FORTRAN_indexing() const LAPACK_WRAPPER_OVERRIDE
    { return this->fortran_indexing; }

    virtual
    void
    transpose() LAPACK_WRAPPER_OVERRIDE {
      this->rows.swap(this->cols);
      std::swap( this->nRows, this->nCols );
    }

    virtual
    void
    to_FORTRAN_indexing() LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    to_C_indexing() LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    push_value_C(
      integer   row,
      integer   col,
      valueType val
    ) LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    push_value_F(
      integer   row,
      integer   col,
      valueType val
    ) LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    get_matrix( MatW & M ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    get_matrix_symmetric( MatW & M ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    get_matrix_transposed( MatW & M ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    bool
    foundNaN() const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    get_data(
      integer   const * & pRows,
      integer   const * & pCols,
      valueType const * & pValues
    ) const LAPACK_WRAPPER_OVERRIDE {
      pRows   = &this->rows.front();
      pCols   = &this->cols.front();
      pValues = &this->vals.front();
    }

    virtual
    void
    gemv(
      valueType       alpha,
      integer         DimX,
      valueType const x[],
      integer         incX,
      valueType       beta,
      integer         DimY,
      valueType       y[],
      integer         incY
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    gemv_Transposed(
      valueType       alpha,
      integer         DimX,
      valueType const x[],
      integer         incX,
      valueType       beta,
      integer         DimY,
      valueType       y[],
      integer         incY
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    gemv_Symmetric(
      valueType       alpha,
      integer         DimX,
      valueType const x[],
      integer         incX,
      valueType       beta,
      integer         DimY,
      valueType       y[],
      integer         incY
    ) const LAPACK_WRAPPER_OVERRIDE;

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
      LAPACK_WRAPPER_ASSERT(
        NNZ == this->rows.size() &&
        NNZ == this->cols.size() &&
        NNZ == this->vals.size(),
        "export_data, bad dimension"
      );
      integer offs = 0;
      if ( fi ) ++offs;
      if ( this->fortran_indexing ) --offs;
      for ( integer index = 0; index < NNZ; ++index ) {
        i[index]   = i_type(this->rows[index]+offs);
        j[index]   = i_type(this->cols[index]+offs);
        val[index] = r_type(this->vals[index]);
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
      integer offs = 0;
      if ( fi ) ++offs;
      if ( this->fortran_indexing ) --offs;
      i.clear(); i.reserve( this->nnz );
      j.clear(); j.reserve( this->nnz );
      v.clear(); v.reserve( this->nnz );
      for ( integer index = 0; index < this->nnz; ++index ) {
        i.push_back( i_type(this->rows[index]+offs) );
        j.push_back( i_type(this->cols[index]+offs) );
        v.push_back( r_type(this->vals[index]) );
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
    fill( valueType const V[], integer M );

    void
    fill( std::vector<valueType> const & V );

    void
    reserve( integer reserve_nnz );

    void
    get_full_view( MatrixWrapper<valueType> & MW ) LAPACK_WRAPPER_OVERRIDE {
      LAPACK_WRAPPER_ASSERT(
        this->matrix_is_full,
        "get_full_view, matrix is sparse"
      );
      if ( this->matrix_is_row_major ) {
        MW.setup( &this->vals.front(), nCols, nRows, nCols );
      } else {
        MW.setup( &this->vals.front(), nRows, nCols, nRows );
      }
    }
  };

  // explicit instantiation declaration to suppress warnings

  #ifdef LAPACK_WRAPPER_USE_CXX11
  extern template class SparseCCOOR<real>;
  extern template class SparseCCOOR<doublereal>;
  #endif

}

///
/// eof: sparse.hxx
///
