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
/// file: lapack_wrapper++.hh
///

#ifndef LAPACK_WRAPPERPP_HH
#define LAPACK_WRAPPERPP_HH

#include "lapack_wrapper.hh"
#include "lapack_wrapper++.hh"
#include "TicToc.hh"

#include <complex>

#ifdef __GNUC__ 
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpadded"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#endif

/*\
:|:    ____            _       _             __
:|:   / ___| _     _  (_)_ __ | |_ ___ _ __ / _| __ _  ___ ___
:|:  | |   _| |_ _| |_| | '_ \| __/ _ \ '__| |_ / _` |/ __/ _ \
:|:  | |__|_   _|_   _| | | | | ||  __/ |  |  _| (_| | (_|  __/
:|:   \____||_|   |_| |_|_| |_|\__\___|_|  |_|  \__,_|\___\___|
\*/

namespace lapack_wrapper {

  /*\
  :|:   __  __       _ _
  :|:  |  \/  | __ _| | | ___   ___
  :|:  | |\/| |/ _` | | |/ _ \ / __|
  :|:  | |  | | (_| | | | (_) | (__
  :|:  |_|  |_|\__,_|_|_|\___/ \___|
  \*/

  //! Allocate memory
  template <typename T>
  class Malloc {
  public:
    typedef T valueType;

  private:

    std::string _name;
    size_t      numTotValues;
    size_t      numTotReserved;
    size_t      numAllocated;
    valueType * pMalloc;

    Malloc(Malloc<T> const &); // blocco costruttore di copia
    Malloc<T> const & operator = (Malloc<T> &) const; // blocco copia

  public:

    //! malloc object constructor
    explicit
    Malloc( std::string const & __name )
    : _name(__name)
    , numTotValues(0)
    , numTotReserved(0)
    , numAllocated(0)
    , pMalloc(nullptr)
    { }

    //! malloc object destructor
    ~Malloc() { free(); }

    //! allocate memory for `n` objects
    void
    allocate( size_t n ) {
      using std::cerr;
      using std::exit;
      try {
        if ( n > numTotReserved ) {
          delete [] pMalloc;
          numTotValues   = n;
          numTotReserved = n + (n>>3); // 12% more values
          pMalloc = new T[numTotReserved];
        }
      }
      catch ( std::exception const & exc ) {
        cerr
          << "Memory allocation failed: " << exc.what()
          << "\nTry to allocate " << n << " bytes for " << _name
          << '\n';
        exit(0);
      }
      catch (...) {
        cerr
          << "Malloc allocation failed for " << _name
          << ": memory exausted\n";
        exit(0);
      }
      numTotValues = n;
      numAllocated = 0;
    }

    //! free memory
    void
    free(void) {
      if ( pMalloc != nullptr ) {
        delete [] pMalloc; pMalloc = nullptr;
        numTotValues   = 0;
        numTotReserved = 0;
        numAllocated   = 0;
      }
    }

    //! number of objects allocated
    size_t size(void) const { return numTotValues; }

    //! get pointer of allocated memory for `sz` objets
    T * operator () ( size_t sz ) {
      size_t offs = numAllocated;
      numAllocated += sz;
      if ( numAllocated > numTotValues ) {
        using std::cerr;
        using std::exit;
        cerr
          << "\nMalloc<" << _name
          << ">::operator () (" << sz << ") -- Malloc EXAUSTED\n";
        exit(0);
      }
      return pMalloc + offs;
    }
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  inline
  void
  print_matrix(
    ostream_type & stream,
    integer        nr,
    integer        nc,
    t_Value const  A[],
    integer        ldA
  ) {
    for ( integer i = 0; i < nr; ++i ) {
      for ( integer j = 0; j < nc; ++j )
        stream << std::setw(14) << A[i+j*ldA] << " ";
      stream << '\n';
    }
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
  class Matrix : public MatrixWrapper<T> {
    lapack_wrapper::Malloc<T> mem;
  public:
    Matrix();
    Matrix( Matrix<T> const & M );
    Matrix( integer nr, integer nc );
    void setup( integer nr, integer nc );
    Matrix<T> const & operator = ( Matrix<T> const & rhs );
  };

  /*\
  :|:   ____  _             __  __       _        _
  :|:  |  _ \(_) __ _  __ _|  \/  | __ _| |_ _ __(_)_  __
  :|:  | | | | |/ _` |/ _` | |\/| |/ _` | __| '__| \ \/ /
  :|:  | |_| | | (_| | (_| | |  | | (_| | |_| |  | |>  <
  :|:  |____/|_|\__,_|\__, |_|  |_|\__,_|\__|_|  |_/_/\_\
  :|:                 |___/
  \*/
  template <typename T>
  class DiagMatrix : public DiagMatrixWrapper<T> {
    Malloc<T> mem;
  public:
    DiagMatrix();
    DiagMatrix( DiagMatrix const & D );

    DiagMatrix( integer _dim );
    void setup( integer _dim );

    void
    scale_by( T sc )
    { scal( this->dim, sc, this->data, 1 ); }

    DiagMatrix<T> const & operator = ( DiagMatrix<T> const & rhs );
  };

  // explicit instantiation declaration to suppress warnings

  #ifdef LAPACK_WRAPPER_USE_CXX11
  extern template class Matrix<real>;
  extern template class Matrix<doublereal>;
  extern template class DiagMatrix<real>;
  extern template class DiagMatrix<doublereal>;
  #endif

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  template <typename T>
  integer
  rankEstimate(
    integer M,
    integer N,
    T       A[],
    integer LDA,
    T       RCOND,
    T       SVAL[3]
  );

  //! base class for linear system solver
  template <typename T>
  class LinearSystemSolver {
  public:
    typedef T valueType;

  public:

    LinearSystemSolver()
    {}

    virtual ~LinearSystemSolver()
    {}

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    void
    solve( valueType xb[] ) const LAPACK_WRAPPER_PURE_VIRTUAL;

    virtual
    void
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_PURE_VIRTUAL;

    virtual
    void
    solve( integer nrhs, valueType B[], integer ldB ) const {
      for ( integer i = 0; i < nrhs; ++i )
        solve( B + i*ldB );
    }

    virtual
    void
    t_solve( integer nrhs, valueType B[], integer ldB ) const {
      for ( integer i = 0; i < nrhs; ++i )
        t_solve( B + i*ldB );
    }

    void
    solve( MatrixWrapper<valueType>  & M )
    { solve( M.numCols(), M.get_data(),  M.lDim() ); }

    void
    t_solve( MatrixWrapper<valueType>  & M )
    { t_solve( M.numCols(), M.get_data(),  M.lDim() ); }

  };

  /*\
  :|:   _____          _             _          _   _
  :|:  |  ___|_ _  ___| |_ ___  _ __(_)______ _| |_(_) ___  _ __
  :|:  | |_ / _` |/ __| __/ _ \| '__| |_  / _` | __| |/ _ \| '_ \
  :|:  |  _| (_| | (__| || (_) | |  | |/ / (_| | |_| | (_) | | | |
  :|:  |_|  \__,_|\___|\__\___/|_|  |_/___\__,_|\__|_|\___/|_| |_|
  \*/

  template <typename T>
  class Factorization : public LinearSystemSolver<T> {
  public:
    typedef T valueType;

  protected:

    valueType * Amat;
    integer     nRow;
    integer     nCol;

  public:

    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;

    Factorization()
    : nRow(0)
    , nCol(0)
    {}

    virtual
    ~Factorization()
    {}

    integer           numRow()   const { return nRow; }
    integer           numCol()   const { return nCol; }
    valueType *       Apointer()       { return Amat; }
    valueType const * Apointer() const { return Amat; }

    void
    zero()
    { gezero( this->nRow, this->nCol, this->Amat, this->nRow ); }

    /*!
    :|:  Zeroes a rectangular block of the stored matrix
    :|:  staring at `(irow,icol)` position
    :|:
    :|:  \param[in] nr    number of rows of the block to be zeroed
    :|:  \param[in] nc    number of columns of the block to be zeroed
    :|:  \param[in] irow  starting row
    :|:  \param[in] icol  stating column
    \*/
    void
    zero_block(
      integer nr,
      integer nc,
      integer irow,
      integer icol
    ) {
      gezero( nr, nc, Amat + irow + icol * nRow, nRow );
    }

    /*!
    :|:  Copy a matrix to a rectangular block of the stored matrix
    :|:  staring at `(irow,icol)` position
    :|:
    :|:  \param[in] B     matrix wrapper of the input matrix `B`
    \*/
    void
    load( MatrixWrapper<T> const & B ) {
      allocate( B.numRows(), B.numCols() );
      integer info = gecopy(
        B.numRows(), B.numCols(),
        B.get_data(), B.lDim(),
        Amat, nRow
      );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "block_load call lapack_wrapper::gecopy return info = " << info
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
      integer info = gecopy(
        nr, nc,
        B, ldB,
        Amat + irow + icol * nRow, nRow
      );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "load_block call lapack_wrapper::gecopy return info = " << info
      );
    }

    /*!
    :|:  Copy vector `column` to the `icol`th column of the internal stored matrix
    :|:  \param[in] column the column vector
    :|:  \param[in] icol   the column to be changed
    \*/
    void
    load_column( valueType const column[], integer icol )
    { copy( this->nRow, column, 1, this->Amat + icol * this->nRow, 1 ); }

    /*!
    :|:  Copy vector `row` to the `irow`th row of the internal stored matrix
    :|:  \param[in] row  the row vector
    :|:  \param[in] irow the row to be changed
    \*/
    void
    load_row( valueType const row[], integer irow )
    { copy( this->nCol, row, 1, this->Amat + irow, this->nRow ); }

    /*!
    :|:  Copy vector element of a sparse vector to a column of
    :|:  the internal stored matrix
    :|:  \param[in] nnz    number of nonzeros of the columns
    :|:  \param[in] values the values of the sparse vector
    :|:  \param[in] row    index position of the values of the sparse vector
    :|:  \param[in] icol   the column to be changed
    :|:  \param[in] offs   offset for the index, 0 for C based vector -1 for FORTRAN based vector
    \*/
    void
    load_sparse_column(
      integer         nnz,
      valueType const values[],
      integer   const row[],
      integer         icol,
      integer         offs = 0
    ) {
      valueType * Acol = Amat + icol * nRow;
      lapack_wrapper::zero( nRow, Acol, 1 );
      for ( integer i = 0; i < nnz; ++i ) Acol[row[i]+offs] = values[i];
    }

    /*!
    :|:  Copy vector element of a sparse vector to a row of the
    :|:  internal stored matrix
    :|:  \param[in] nnz    number of nonzeros of the columns
    :|:  \param[in] values the values of the sparse vector
    :|:  \param[in] col    index position of the values of the sparse vector
    :|:  \param[in] irow   the column to be changed
    :|:  \param[in] offs   offset for the index, 0 for C based vector -1 for FORTRAN based vector
    \*/
    void
    load_sparse_row(
      integer         nnz,
      valueType const values[],
      integer   const col[],
      integer         irow,
      integer         offs = 0
    ) {
      valueType * Arow = Amat + irow;
      lapack_wrapper::zero( nRow, Arow, nRow );
      for ( integer i = 0; i < nnz; ++i ) Arow[col[i]+offs] = values[i];
    }

    /*!
    :|:  Copy a sparse matrix into the internal stored matrix
    :|:  \param[in] nnz    number of nonzeros of the columns
    :|:  \param[in] values the values of the sparse vector
    :|:  \param[in] row    index row position of the values of the sparse vector
    :|:  \param[in] col    index column position of the values of the sparse vector
    \*/
    void
    load_sparse(
      integer         nnz,
      valueType const values[],
      integer   const row[],
      integer   const col[]
    ) {
      lapack_wrapper::zero( nRow*nCol, Amat, 1 );
      for ( integer i = 0; i < nnz; ++i )
        Amat[row[i] + col[i] * nRow] = values[i];
    }

    /*!
    :|:  Copy a sparse matrix into the internal stored matrix
    :|:  \param[in] nnz    number of nonzeros of the columns
    :|:  \param[in] values the values of the sparse vector
    :|:  \param[in] row    index row position of the values of the sparse vector
    :|:  \param[in] r_offs offset for the index, 0 for C based vector -1 for FORTRAN based vector
    :|:  \param[in] col    index column position of the values of the sparse vector
    :|:  \param[in] c_offs offset for the index, 0 for C based vector -1 for FORTRAN based vector
    \*/
    void
    load_sparse(
      integer         nnz,
      valueType const values[],
      integer   const row[], integer r_offs,
      integer   const col[], integer c_offs
    ) {
      lapack_wrapper::zero( nRow*nCol, Amat, 1 );
      for ( integer i = 0; i < nnz; ++i )
        Amat[row[i]+r_offs + (col[i]+c_offs) * nRow] = values[i];
    }

    /*!
    :|:  Copy a sparse matrix into the internal stored matrix.
    :|:  The matrix is assumed symmetric and only the lower or upper
    :|:  part is passed.
    :|:
    :|:  \param[in] nnz    number of nonzeros of the columns
    :|:  \param[in] values the values of the sparse vector
    :|:  \param[in] row    index row position of the values of the sparse vector
    :|:  \param[in] col    index column position of the values of the sparse vector
    \*/
    void
    load_sparse_sym(
      integer         nnz,
      valueType const values[],
      integer   const row[],
      integer   const col[]
    ) {
      lapack_wrapper::zero( nRow*nCol, Amat, 1 );
      for ( integer i = 0; i < nnz; ++i ) {
        integer ii = row[i];
        integer jj = col[i];
        Amat[ii + jj * nRow] = values[i];
        if ( ii != jj ) Amat[ jj + ii * nRow] = values[i];
      }
    }

    /*!
    :|:  Copy a sparse matrix into the internal stored matrix.
    :|:  The matrix is assumed symmetric and only the lower or upper
    :|:  part is passed.
    :|:
    :|:  \param[in] nnz    number of nonzeros of the columns
    :|:  \param[in] values the values of the sparse vector
    :|:  \param[in] row    index row position of the values of the sparse vector
    :|:  \param[in] r_offs offset for the index, 0 for C based vector -1 for FORTRAN based vector
    :|:  \param[in] col    index column position of the values of the sparse vector
    :|:  \param[in] c_offs offset for the index, 0 for C based vector -1 for FORTRAN based vector
    \*/
    void
    load_sparse_sym(
      integer         nnz,
      valueType const values[],
      integer   const row[], integer r_offs,
      integer   const col[], integer c_offs
    ) {
      lapack_wrapper::zero( nRow*nCol, Amat, 1 );
      for ( integer i = 0; i < nnz; ++i ) {
        integer ii = row[i]+r_offs;
        integer jj = col[i]+c_offs;
        Amat[ii + jj * nRow] = values[i];
        if ( ii != jj ) Amat[ jj + ii * nRow] = values[i];
      }
    }

    void
    add_to_diag( valueType mu ) {
      lapack_wrapper::axpy( std::min(nRow,nCol), 1.0, &mu, 0, Amat, nRow+1 );
    }

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    void
    allocate( integer NR, integer NC ) LAPACK_WRAPPER_PURE_VIRTUAL;

    virtual
    void
    factorize( char const who[] ) LAPACK_WRAPPER_PURE_VIRTUAL;

    virtual
    void
    factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) LAPACK_WRAPPER_PURE_VIRTUAL;

    void
    factorize( char const who[], MatrixWrapper<valueType> const & M )
    { factorize( who, M.numRows(), M.numCols(), M.get_data(), M.lDim() ); }

  };

}

#include "code++/lu.hxx"
#include "code++/qr.hxx"
#include "code++/svd.hxx"
#include "code++/ls.hxx"
#include "code++/trid.hxx"
#include "code++/band.hxx"
#include "code++/qn.hxx"
#include "code++/eig.hxx"
#include "code++/block_trid.hxx"

namespace lapack_wrapper {


  //============================================================================
  /*\
  :|:   _     ____   ___   ____
  :|:  | |   / ___| / _ \ / ___|
  :|:  | |   \___ \| | | | |
  :|:  | |___ ___) | |_| | |___
  :|:  |_____|____/ \__\_\\____|
  \*/
  /*!  least square constained
  :|:
  :|:  minimize    (1/2) * || A * x - b ||^2
  :|:
  :|:  subject to  B * x = c
  :|:
  :|:  L = (1/2) * || A * x - b ||^2 + lambda * ( B * x - c )
  :|:
  :|:  is equivalent to solve linear system
  :|:
  :|:  / A^T*A  B^T \  /   x    \   / A^T b \
  :|:  |            |  |        | = |       |
  :|:  \ B       0  /  \ lambda /   \   c   /
  :|:
  :|:  is equivalent to solve linear system
  :|:
  :|:  / 0   A^T  B^T \  /   x    \   / A^T b \
  :|:  |              |  |        |   |       |
  :|:  | A   -I    0  |  | omega  | = |   0   |
  :|:  |              |  |        |   |       |
  :|:  \ B    0    0  /  \ lambda /   \   c   /
  :|:
  \*/

  #if 0

  template <typename T>
  class LSQC {
  public:
    typedef T valueType;

    Malloc<valueType> allocReals;
    Malloc<integer>   allocIntegers;

    integer     NR1, NR2;
    integer   * to_row1;
    integer   * to_row2;
    valueType * rhs;
    LU<T> lu;

  public:

    LSQC()
    : allocReals("LSQC-allocReals")
    , allocIntegers("LSQC-allocIntegers")
    , to_row1(nullptr)
    , to_row2(nullptr)
    , rhs(nullptr)
    {}

    ~LSQC() {}

    void
    factorize(
      integer         NR,
      integer         NC,
      valueType const M[],
      integer         ldM,
      bool const      row_select[] // row selected for minimization
    );

    void
    solve( valueType xb[] ) const;

    void
    solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const;

  };

  #endif

  // explicit instantiation declaration to suppress warnings

  #ifdef LAPACK_WRAPPER_USE_CXX11

  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  #endif

  extern template class LU<real>;
  extern template class LU<doublereal>;

  extern template class LUPQ<real>;
  extern template class LUPQ<doublereal>;

  extern template class QR<real>;
  extern template class QR<doublereal>;

  extern template class QRP<real>;
  extern template class QRP<doublereal>;

  extern template class SVD<real>;
  extern template class SVD<doublereal>;

  extern template class LSS<real>;
  extern template class LSS<doublereal>;

  extern template class LSY<real>;
  extern template class LSY<doublereal>;

  extern template class TridiagonalSPD<real>;
  extern template class TridiagonalSPD<doublereal>;

  extern template class TridiagonalLU<real>;
  extern template class TridiagonalLU<doublereal>;

  extern template class TridiagonalQR<real>;
  extern template class TridiagonalQR<doublereal>;

  extern template class BlockTridiagonalSymmetic<real>;
  extern template class BlockTridiagonalSymmetic<doublereal>;

  extern template class BandedLU<real>;
  extern template class BandedLU<doublereal>;

  extern template class BandedSPD<real>;
  extern template class BandedSPD<doublereal>;

  extern template class QN<real>;
  extern template class QN<doublereal>;

  extern template class BFGS<real>;
  extern template class BFGS<doublereal>;

  extern template class DFP<real>;
  extern template class DFP<doublereal>;

  extern template class Eigenvalues<real>;
  extern template class Eigenvalues<doublereal>;

  extern template class Eigenvectors<real>;
  extern template class Eigenvectors<doublereal>;

  extern template class GeneralizedEigenvalues<real>;
  extern template class GeneralizedEigenvalues<doublereal>;

  extern template class GeneralizedEigenvectors<real>;
  extern template class GeneralizedEigenvectors<doublereal>;

  extern template class GeneralizedSVD<real>;
  extern template class GeneralizedSVD<doublereal>;

  #endif

} // end namespace lapack_wrapper

#ifdef __GNUC__ 
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif

#endif

///
/// eof: lapack_wrapper++.hh
///
