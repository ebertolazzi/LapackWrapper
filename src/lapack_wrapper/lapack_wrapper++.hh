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
    void scale_by( T );
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

  //============================================================================
  /*\
  :|:   _    _   _
  :|:  | |  | | | |
  :|:  | |  | | | |
  :|:  | |__| |_| |
  :|:  |_____\___/
  \*/

  template <typename T>
  class LU : public Factorization<T> {
  public:
    typedef typename Factorization<T>::valueType valueType;

  private:

    valueType * Work;
    integer   * Iwork;
    integer   * i_pivot;

    Malloc<valueType> allocReals;
    Malloc<integer>   allocIntegers;

    void check_ls( char const who[] ) const;

  public:

    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;
    using Factorization<T>::factorize;
    using Factorization<T>::solve;
    using Factorization<T>::t_solve;

    using Factorization<T>::nRow;
    using Factorization<T>::nCol;
    using Factorization<T>::Amat;

    LU();
    virtual ~LU() LAPACK_WRAPPER_OVERRIDE;

    valueType cond1( valueType norm1 ) const;
    valueType condInf( valueType normInf ) const;

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    void
    allocate( integer NR, integer NC ) LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    factorize( char const who[] ) LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const LAPACK_WRAPPER_OVERRIDE;

  };

  //============================================================================
  /*\
  :|:   _    _   _ ____   ___
  :|:  | |  | | | |  _ \ / _ \
  :|:  | |  | | | | |_) | | | |
  :|:  | |__| |_| |  __/| |_| |
  :|:  |_____\___/|_|    \__\_\
  \*/
  template <typename T>
  class LUPQ : public Factorization<T> {
  public:
    typedef typename Factorization<T>::valueType valueType;

  private:

    integer * ipiv;
    integer * jpiv;

    Malloc<valueType> allocReals;
    Malloc<integer>   allocIntegers;

    void check_ls( char const who[] ) const;

  public:

    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;
    using Factorization<T>::factorize;
    using Factorization<T>::solve;
    using Factorization<T>::t_solve;

    using Factorization<T>::nRow;
    using Factorization<T>::nCol;
    using Factorization<T>::Amat;

    LUPQ();
    virtual ~LUPQ() LAPACK_WRAPPER_OVERRIDE;

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    void
    allocate( integer NR, integer NC ) LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    factorize( char const who[] ) LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const LAPACK_WRAPPER_OVERRIDE;

  };

  //============================================================================
  /*\
  :|:    ___  ____
  :|:   / _ \|  _ \
  :|:  | | | | |_) |
  :|:  | |_| |  _ <
  :|:   \__\_\_| \_\
  \*/

  template <typename T>
  class QR : public Factorization<T> {
  public:
    typedef typename Factorization<T>::valueType valueType;

  private:

    Malloc<valueType> allocReals;

  protected:

    valueType * Work;
    valueType * Tau;

    integer nReflector;
    integer Lwork;
    integer maxNrhs;

  public:

    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;
    using Factorization<T>::factorize;
    using Factorization<T>::solve;
    using Factorization<T>::t_solve;

    using Factorization<T>::nRow;
    using Factorization<T>::nCol;
    using Factorization<T>::Amat;

    QR()
    : Factorization<T>()
    , allocReals("QR-allocReals")
    , nReflector(0)
    , Lwork(0)
    , maxNrhs(1)
    {}

    QR( integer nr, integer nc )
    : Factorization<T>()
    , allocReals("QR-allocReals")
    , nReflector(0)
    , Lwork(0)
    , maxNrhs(1)
    { this->allocate(nr,nc); }

    virtual
    ~QR() LAPACK_WRAPPER_OVERRIDE
    { allocReals.free(); }

    void
    setMaxNrhs( integer mnrhs );

    void
    allocate( integer nr, integer nc, integer Lwrk );

    /*\
    :|:  overwrites the general real M-by-N matrix C with
    :|:
    :|:  SIDE = 'L'     SIDE = 'R'
    :|:  TRANS = 'N':      Q * C          C * Q
    :|:  TRANS = 'T':      Q**T * C       C * Q**T
    \*/
    void
    applyQ(
      SideMultiply  SIDE,
      Transposition TRANS,
      integer       nRefl,
      integer       nr,
      integer       nc,
      valueType     C[],
      integer       ldC
    ) const;

    //! x <- Q*x
    void
    Q_mul( valueType x[] ) const
    { applyQ( LEFT, NO_TRANSPOSE, nReflector, nRow, 1, x, nRow ); }

    //! x <- Q'*x
    void
    Qt_mul( valueType x[] ) const
    { applyQ( LEFT, TRANSPOSE, nReflector, nRow, 1, x, nRow ); }

    //! C <- Q*C
    void
    Q_mul( integer nr, integer nc, valueType C[], integer ldC ) const
    { applyQ( LEFT, NO_TRANSPOSE, nReflector, nr, nc, C, ldC ); }

    //! C <- Q'*C
    void
    Qt_mul( integer nr, integer nc, valueType C[], integer ldC ) const
    { applyQ( LEFT, TRANSPOSE, nReflector, nr, nc, C, ldC ); }

    //! C <- C*Q
    void
    mul_Q( integer nr, integer nc, valueType C[], integer ldC ) const
    { applyQ( RIGHT, NO_TRANSPOSE, nReflector, nr, nc, C, ldC ); }

    //! C <- C*Q'
    void
    mul_Qt( integer nr, integer nc, valueType C[], integer ldC ) const
    { applyQ( RIGHT, TRANSPOSE, nReflector, nr, nc, C, ldC ); }

    // -------------------------------------------------------------------------

    void
    Rsolve(
      Transposition TRANS,
      integer       rk,
      valueType     x[],
      integer       incx
    ) const {
      trsv( UPPER, TRANS, NON_UNIT, rk, Amat, nRow, x, incx );
    }

    void
    Rsolve(
      SideMultiply  SIDE,
      Transposition TRANS,
      integer       nr,
      integer       nc,
      valueType     alpha,
      valueType     Bmat[],
      integer       ldB
    ) const {
      trsm( SIDE, UPPER, TRANS, NON_UNIT,
            nr, nc, alpha, Amat, nRow, Bmat, ldB );
    }

    //! x <- R^(-1) * x
    void
    invR_mul( valueType x[], integer incx = 1 ) const
    { trsv( UPPER, NO_TRANSPOSE, NON_UNIT,
            nReflector, Amat, nRow, x, incx ); }

    //! x <- R^(-T) * x
    void
    invRt_mul( valueType x[], integer incx = 1 ) const
    { trsv( UPPER, TRANSPOSE, NON_UNIT,
            nReflector, Amat, nRow, x, incx ); }

    //! C <- R^(-1) * C
    void
    invR_mul( integer nr, integer nc, valueType C[], integer ldC ) const
    { Rsolve( LEFT, NO_TRANSPOSE, nr, nc, 1.0, C, ldC ); }

    //! C <- R^(-T) * C
    void
    invRt_mul( integer nr, integer nc, valueType C[], integer ldC ) const
    { Rsolve( LEFT, TRANSPOSE, nr, nc, 1.0, C, ldC ); }

    //! C <- C * R^(-1)
    void
    mul_invR( integer nr, integer nc, valueType C[], integer ldC ) const
    { Rsolve( RIGHT, NO_TRANSPOSE, nr, nc, 1.0, C, ldC ); }

    void
    mul_invRt( integer nr, integer nc, valueType C[], integer ldC ) const
    { Rsolve( RIGHT, TRANSPOSE,  nr, nc, 1.0, C, ldC ); }

    void
    getR( valueType R[], integer ldR ) const;

    // -------------------------------------------------------------------------
    // dummy routines
    void permute( valueType [] ) const {}
    void inv_permute( valueType [] ) const {}
    void permute_rows( integer, integer, valueType [], integer ) const {}
    void inv_permute_rows( integer, integer, valueType [], integer ) const {}

    /*!
      Do QR factorization of the transpose of a rectangular matrix
      \param NR  number of rows of the matrix
      \param NC  number of columns of the matrix
      \param A   pointer to the matrix
      \param LDA Leading dimension of the matrix
    \*/
    void
    t_factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) {
      allocate( NC, NR );
      for ( integer i = 0; i < NR; ++i )
        copy( NC, A+i, LDA, Amat + i*nRow, 1 );
      factorize( who );
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
    allocate( integer nr, integer nc ) LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    factorize( char const who[] ) LAPACK_WRAPPER_OVERRIDE {
      integer info = geqrf(
        nRow, nCol, Amat, nRow, Tau, Work, Lwork
      );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "QR::factorize[" << who << 
        "] call lapack_wrapper::geqrf return info = " << info
      );
    }

    /*!
    :|:  Do QR factorization of a rectangular matrix
    :|:  \param NR  number of rows of the matrix
    :|:  \param NC  number of columns of the matrix
    :|:  \param A   pointer to the matrix
    :|:  \param LDA Leading dimension of the matrix
    \*/
    virtual
    void
    factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) LAPACK_WRAPPER_OVERRIDE {
      allocate( NR, NC );
      integer info = gecopy( NR, NC, A, LDA, Amat, nRow );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "QR::factorize[" << who << 
        "] call lapack_wrapper::gecopy return info = " << info
      );
      factorize( who );
    }

    /*!
    :|:  In case of QR factorization of a square matrix solve the
    :|:  linear system \f$ QR x = b \f$
    :|:  \param xb on input the rhs of linear system on output the solution
    \*/
    virtual
    void
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    /*!
     |  In case of QR factorization of a square matrix solve the
     |  linear system \f$ (QR)^T x = b \f$
     |  \param xb on input the rhs of linear system on output the solution
    \*/
    virtual
    void
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const LAPACK_WRAPPER_OVERRIDE;

  };

  //============================================================================
  /*\
   |    ___  ____  ____
   |   / _ \|  _ \|  _ \
   |  | | | | |_) | |_) |
   |  | |_| |  _ <|  __/
   |   \__\_\_| \_\_|
  \*/
  template <typename T>
  class QRP : public QR<T> {
  public:
    typedef typename QR<T>::valueType valueType;

  private:
    Malloc<integer> allocIntegers;
    integer * JPVT;

  public:

    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;
    using Factorization<T>::factorize;
    using Factorization<T>::solve;
    using Factorization<T>::t_solve;

    using Factorization<T>::nRow;
    using Factorization<T>::nCol;
    using Factorization<T>::Amat;

    using QR<T>::setMaxNrhs;

    using QR<T>::Lwork;
    using QR<T>::Work;
    using QR<T>::maxNrhs;
    using QR<T>::Tau;
    using QR<T>::Q_mul;
    using QR<T>::Qt_mul;
    using QR<T>::invR_mul;
    using QR<T>::invRt_mul;

    QRP()
    : QR<T>(), allocIntegers("QRP-allocIntegers")
    {}

    QRP( integer nr, integer nc )
    : QR<T>(), allocIntegers("QRP-allocIntegers")
    { this->allocate(nr,nc); }

    virtual
    ~QRP() LAPACK_WRAPPER_OVERRIDE
    { allocIntegers.free(); }

    // -------------------------------------------------------------------------
    void
    permute( valueType x[] ) const;

    void
    inv_permute( valueType x[] ) const;

    void
    permute_rows(
      integer   nr,
      integer   nc,
      valueType C[],
      integer   ldC
    ) const {
      LAPACK_WRAPPER_ASSERT(
        nr == nRow,
        "QRP::permute_rows, bad number of row, expected " <<
        nRow << " find " << nr
      );
      for ( integer j = 0; j < nc; ++j ) permute( C + ldC*j );
    }

    void
    inv_permute_rows(
      integer   nr,
      integer   nc,
      valueType C[],
      integer   ldC
    ) const {
      LAPACK_WRAPPER_ASSERT(
        nr == nRow,
        "QRP::permute_rows, bad number of row, expected " <<
        nRow << " find " << nr
      );
      for ( integer j = 0; j < nc; ++j ) inv_permute( C + ldC*j );
    }

    integer
    rankEstimate( valueType rcond ) const {
      valueType SVAL[3];
      return lapack_wrapper::rankEstimate(
        nRow, nCol, Amat, nRow, rcond, SVAL
      );
    }

    /*!
    :|:  Do QR factorization with column pivoting of the transpose of a rectangular matrix
    :|:  \param NR  number of rows of the matrix
    :|:  \param NC  number of columns of the matrix
    :|:  \param A   pointer to the matrix
    :|:  \param LDA Leading dimension of the matrix
    \*/
    void
    t_factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) {
      // calcolo fattorizzazione QR della matrice A
      this->allocate( NC, NR );
      for ( integer i = 0; i < NR; ++i )
        copy( NC, A+i, LDA, Amat + i*nRow, 1 );
      factorize( who );
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
    allocate( integer NR, integer NC ) LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    factorize( char const who[] ) LAPACK_WRAPPER_OVERRIDE {
      std::fill( JPVT, JPVT + nCol, 0 );
      integer info = geqp3(
        nRow, nCol, Amat, nRow, JPVT, Tau, Work, Lwork
      );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "QRP::factorize[" << who << 
        "] call lapack_wrapper::geqrf return info = " << info
      );
    }

    /*!
    :|:  Do QR factorization with column pivoting of a rectangular matrix
    :|:  \param NR  number of rows of the matrix
    :|:  \param NC  number of columns of the matrix
    :|:  \param A   pointer to the matrix
    :|:  \param LDA Leading dimension of the matrix
    \*/
    virtual
    void
    factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) LAPACK_WRAPPER_OVERRIDE {
      // calcolo fattorizzazione QR della matrice A
      allocate( NC, NR );
      integer info = gecopy( NR, NC, A, LDA, Amat, nRow );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "QR::factorize[" << who << 
        "] call lapack_wrapper::gecopy return info = " << info
      );
      factorize( who );
    }

    /*!
    :|:  In case of QR factorization of a square matrix solve the
    :|:  linear system \f$ QR x = b \f$
    :|:  \param xb on input the rhs of linear system on output the solution
    \*/
    virtual
    void
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    /*!
    :|:  In case of QR factorization of a square matrix solve the
    :|:  linear system \f$ (QR)^T x = b \f$
    :|:  \param xb on input the rhs of linear system on output the solution
    \*/
    virtual
    void
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const LAPACK_WRAPPER_OVERRIDE;
  };

  //============================================================================
  /*\
  :|:   ______     ______
  :|:  / ___\ \   / /  _ \
  :|:  \___ \\ \ / /| | | |
  :|:   ___) |\ V / | |_| |
  :|:  |____/  \_/  |____/
  \*/
  template <typename T>
  class SVD : public Factorization<T> {
  public:
    typedef typename Factorization<T>::valueType valueType;

  protected:

    Malloc<valueType> allocReals;
    Malloc<integer>   allocIntegers;

    valueType * Work;
    valueType * Umat;
    valueType * VTmat;
    valueType * Svec;
    integer   * IWork;

    valueType   rcond;

    integer     minRC, Lwork;

    typedef enum { USE_GESVD = 0, USE_GESDD = 1 } SVD_USED;
    SVD_USED    svd_used;

  public:

    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;
    using Factorization<T>::factorize;
    using Factorization<T>::solve;
    using Factorization<T>::t_solve;

    using Factorization<T>::nRow;
    using Factorization<T>::nCol;
    using Factorization<T>::Amat;

    SVD( SVD_USED _svd_used = USE_GESVD )
    : Factorization<T>()
    , allocReals("SVD-allocReals")
    , allocIntegers("SVD-allocIntegers")
    , rcond(machineEps<valueType>())
    , Lwork(0)
    , svd_used(_svd_used)
    {}

    virtual
    ~SVD() LAPACK_WRAPPER_OVERRIDE
    { allocReals.free(); allocIntegers.free(); }

    void
    setRcond( valueType r )
    { rcond = r; }

    valueType U( integer i, integer j ) const { return Umat[i+j*nRow]; }
    valueType V( integer i, integer j ) const { return VTmat[j+i*nCol]; }
    valueType sigma( integer i ) const { return Svec[i]; }

    //! y <- alpha * U * x + beta * y
    void
    U_mul(
      valueType       alpha,
      valueType const x[],
      integer         incx,
      valueType       beta,
      valueType       y[],
      integer incy
    ) const {
      gemv(
        NO_TRANSPOSE,
        nRow, minRC,
        alpha, Umat, nRow,
        x, incx,
        beta, y, incy
      );
    }

    //! y <- alpha * U' * x + beta * y
    void
    Ut_mul(
      valueType       alpha,
      valueType const x[],
      integer         incx,
      valueType       beta,
      valueType       y[],
      integer         incy
    ) const {
      gemv(
        TRANSPOSE,
        nRow, minRC,
        alpha, Umat, nRow,
        x, incx,
        beta, y, incy
      );
    }

    //! y <- alpha * V * x + beta * y
    void
    V_mul(
      valueType       alpha,
      valueType const x[],
      integer         incx,
      valueType       beta,
      valueType       y[],
      integer         incy
    ) const {
      gemv(
        TRANSPOSE,
        minRC, nCol,
        alpha, VTmat, nRow,
        x, incx,
        beta, y, incy
      );
    }

    //! y <- alpha * V' * x + beta * y
    void
    Vt_mul(
      valueType       alpha,
      valueType const x[],
      integer         incx,
      valueType       beta,
      valueType       y[],
      integer         incy
    ) const {
      gemv(
        NO_TRANSPOSE,
        minRC, nCol,
        alpha, VTmat, nRow,
        x, incx,
        beta, y, incy
      );
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
    allocate( integer NR, integer NC ) LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    factorize( char const who[] ) LAPACK_WRAPPER_OVERRIDE;

    /*!
    :|:  Do SVD factorization of a rectangular matrix
    :|:  \param NR  number of rows of the matrix
    :|:  \param NC  number of columns of the matrix
    :|:  \param A   pointer to the matrix
    :|:  \param LDA Leading dimension of the matrix
    \*/
    virtual
    void
    factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) LAPACK_WRAPPER_OVERRIDE {
      allocate( NR, NC );
      integer info = gecopy( nRow, nCol, A, LDA, Amat, nRow );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "SVD::factorize[" << who << 
        "] call lapack_wrapper::gecopy return info = " << info
      );
      factorize( who );
    }

    virtual
    void
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

  };

  //============================================================================
  /*\
  :|:   _     ____ ____
  :|:  | |   / ___/ ___|
  :|:  | |   \___ \___ \
  :|:  | |___ ___) |__) |
  :|:  |_____|____/____/
  \*/

  template <typename T>
  class LSS : public Factorization<T> {
  public:
    typedef typename Factorization<T>::valueType valueType;

  protected:

    Malloc<valueType> allocReals;

    mutable valueType * Work;
    mutable valueType * sigma;
    mutable valueType * AmatWork;
    mutable integer     rank;

    valueType rcond;
    integer   Lwork;
    integer   maxNrhs;
    bool      maxNrhs_changed;

  public:

    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;
    using Factorization<T>::factorize;
    using Factorization<T>::solve;
    using Factorization<T>::t_solve;

    using Factorization<T>::nRow;
    using Factorization<T>::nCol;
    using Factorization<T>::Amat;

    explicit
    LSS()
    : Factorization<T>()
    , allocReals("LSS-allocReals")
    , Work(nullptr)
    , sigma(nullptr)
    , AmatWork(nullptr)
    , rank(0)
    , rcond(-1)
    , Lwork(0)
    , maxNrhs(1)
    , maxNrhs_changed(true)
    {}

    void
    setMaxNrhs( integer mnrhs );

    virtual
    ~LSS() LAPACK_WRAPPER_OVERRIDE
    { allocReals.free(); }

    void
    setRcond( valueType r )
    { rcond = r; }

    integer
    getRank() const
    { return rank; }

    valueType
    getSigma( integer i ) const
    { return sigma[i]; }

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    void
    allocate( integer NR, integer NC ) LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    factorize( char const [] ) LAPACK_WRAPPER_OVERRIDE {
      // nothing to do
    }

    /*!
    :|:  Do SVD factorization of a rectangular matrix
    :|:  \param NR  number of rows of the matrix
    :|:  \param NC  number of columns of the matrix
    :|:  \param A   pointer to the matrix
    :|:  \param LDA Leading dimension of the matrix
    \*/
    virtual
    void
    factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) LAPACK_WRAPPER_OVERRIDE {
      allocate( NR, NC );
      integer info = gecopy( nRow, nCol, A, LDA, Amat, nRow );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "LSS::factorize[" << who << 
        "] call lapack_wrapper::gecopy return info = " << info
      );
      factorize( who );
    }

    virtual
    void
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    solve( integer nrhs, valueType B[], integer ldB ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve( integer nrhs, valueType B[], integer ldB ) const LAPACK_WRAPPER_OVERRIDE;

  };

  //============================================================================
  /*\
  :|:  _     ______   __
  :|: | |   / ___\ \ / /
  :|: | |   \___ \\ V /
  :|: | |___ ___) || |
  :|: |_____|____/ |_|
  \*/
  template <typename T>
  class LSY : public Factorization<T> {
  public:
    typedef typename Factorization<T>::valueType valueType;

  protected:

    Malloc<valueType> allocReals;
    Malloc<integer>   allocInts;

    mutable valueType * Work;
    mutable valueType * AmatWork;
    mutable integer   * jpvt;
    mutable integer     rank;

    valueType rcond;
    integer   Lwork;
    integer   maxNrhs;
    bool      maxNrhs_changed;

  public:

    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;
    using Factorization<T>::factorize;
    using Factorization<T>::solve;
    using Factorization<T>::t_solve;

    using Factorization<T>::nRow;
    using Factorization<T>::nCol;
    using Factorization<T>::Amat;

    explicit
    LSY()
    : Factorization<T>()
    , allocReals("LSY-allocReals")
    , allocInts("LSY-allocInts")
    , Work(nullptr)
    , AmatWork(nullptr)
    , rank(0)
    , rcond(-1)
    , Lwork(0)
    , maxNrhs(1)
    , maxNrhs_changed(true)
    {}

    virtual
    ~LSY() LAPACK_WRAPPER_OVERRIDE
    { allocReals.free(); allocInts.free(); }

    void
    setMaxNrhs( integer mnrhs );

    void
    setRcond( valueType r )
    { rcond = r; }

    integer
    getRank() const
    { return rank; }

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    void
    allocate( integer NR, integer NC ) LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    factorize( char const [] ) LAPACK_WRAPPER_OVERRIDE {
      // nothing to do
    }

    /*!
    :|:  Do SVD factorization of a rectangular matrix
    :|:  \param NR  number of rows of the matrix
    :|:  \param NC  number of columns of the matrix
    :|:  \param A   pointer to the matrix
    :|:  \param LDA Leading dimension of the matrix
    \*/
    virtual
    void
    factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) LAPACK_WRAPPER_OVERRIDE {
      allocate( NR, NC );
      integer info = gecopy( nRow, nCol, A, LDA, Amat, nRow );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "LSY::factorize[" << who << 
        "] call lapack_wrapper::gecopy return info = " << info
      );
      factorize( who );
    }

    virtual
    void
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    solve( integer nrhs, valueType B[], integer ldB ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve( integer nrhs, valueType B[], integer ldB ) const LAPACK_WRAPPER_OVERRIDE;

  };

  //============================================================================
  /*\
  :|:   _____     _     _ _                               _ ____  ____  ____
  :|:  |_   _| __(_) __| (_) __ _  __ _  ___  _ __   __ _| / ___||  _ \|  _ \
  :|:    | || '__| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | \___ \| |_) | | | |
  :|:    | || |  | | (_| | | (_| | (_| | (_) | | | | (_| | |___) |  __/| |_| |
  :|:    |_||_|  |_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|____/|_|   |____/
  :|:                             |___/
  \*/

  template <typename T>
  class TridiagonalSPD : public LinearSystemSolver<T> {
  public:
    typedef T valueType;

  private:

    Malloc<valueType> allocReals;

    valueType * L;
    valueType * D;
    valueType * WORK;
    integer     nRC;

  public:

    TridiagonalSPD()
    : allocReals("allocReals")
    , nRC(0)
    {}

    virtual
    ~TridiagonalSPD() LAPACK_WRAPPER_OVERRIDE
    { allocReals.free(); }

    valueType cond1( valueType norm1 ) const;

    void
    factorize(
      char const      who[],
      integer         N,
      valueType const _L[],
      valueType const _D[]
    );

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    void
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    solve(
      integer   nrhs,
      valueType xb[],
      integer   ldXB
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve(
      integer   nrhs,
      valueType xb[],
      integer   ldXB
    ) const LAPACK_WRAPPER_OVERRIDE;

    /*\
    :|:     _
    :|:    / \  _   ___  __
    :|:   / _ \| | | \ \/ /
    :|:  / ___ \ |_| |>  <
    :|: /_/   \_\__,_/_/\_\
    :|:
    \*/

    void
    axpy(
      integer         N,
      valueType       alpha,
      valueType const L[],
      valueType const D[],
      valueType const x[],
      valueType       beta,
      valueType       y[]
    ) const;

  };

  //============================================================================
  /*\
  :|:   _____     _     _ _                               _ _    _   _
  :|:  |_   _| __(_) __| (_) __ _  __ _  ___  _ __   __ _| | |  | | | |
  :|:    | || '__| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | | |  | | | |
  :|:    | || |  | | (_| | | (_| | (_| | (_) | | | | (_| | | |__| |_| |
  :|:    |_||_|  |_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|_____\___/
  :|:                             |___/
  \*/

  template <typename T>
  class TridiagonalLU : public LinearSystemSolver<T> {
  public:
    typedef T valueType;

  private:

    Malloc<valueType> allocReals;
    Malloc<integer>   allocIntegers;

    valueType * L;
    valueType * D;
    valueType * U;
    valueType * U2;
    valueType * WORK;
    integer   * IPIV;
    integer   * IWORK;

    integer     nRC;

  public:

    TridiagonalLU()
    : allocReals("TridiagonalLU-allocReals")
    , allocIntegers("TridiagonalLU-allocIntegers")
    , nRC(0)
    {}

    virtual
    ~TridiagonalLU() LAPACK_WRAPPER_OVERRIDE {
      allocReals.free();
      allocIntegers.free();
    }

    valueType cond1( valueType norm1 ) const;
    valueType condInf( valueType normInf ) const;

    void
    factorize(
      char const      who[],
      integer         N,
      valueType const _L[],
      valueType const _D[],
      valueType const _U[]
    );

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    void
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    solve(
      integer   nrhs,
      valueType xb[],
      integer   ldXB
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve(
      integer   nrhs,
      valueType xb[],
      integer   ldXB
    ) const LAPACK_WRAPPER_OVERRIDE;

    /*\
    :|:     _
    :|:    / \  _   ___  __
    :|:   / _ \| | | \ \/ /
    :|:  / ___ \ |_| |>  <
    :|: /_/   \_\__,_/_/\_\
    :|:
    \*/

    void
    axpy(
      integer         N,
      valueType       alpha,
      valueType const L[],
      valueType const D[],
      valueType const U[],
      valueType const x[],
      valueType       beta,
      valueType       y[]
    ) const;

  };

  //============================================================================
  /*\
  :|:   _____     _     _ _                               _  ___  ____
  :|:  |_   _| __(_) __| (_) __ _  __ _  ___  _ __   __ _| |/ _ \|  _ \
  :|:    | || '__| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | | | | | |_) |
  :|:    | || |  | | (_| | | (_| | (_| | (_) | | | | (_| | | |_| |  _ <
  :|:    |_||_|  |_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|\__\_\_| \_\
  :|:                             |___/
  \*/
  template <typename T>
  class TridiagonalQR : public LinearSystemSolver<T> {
  public:
    typedef T valueType;

  private:

    Malloc<valueType> allocReals;

    valueType * C;   // rotazioni givens
    valueType * S;
    valueType * BD;  // band triangular matrix
    valueType * BU;  // band triangular matrix
    valueType * BU2; // band triangular matrix

    valueType   normInfA;
    integer     nRC;

    void Rsolve( valueType xb[] ) const;
    void RsolveTransposed( valueType xb[] ) const;

  public:

    TridiagonalQR()
    : allocReals("allocReals")
    , nRC(0)
    {}

    virtual
    ~TridiagonalQR() LAPACK_WRAPPER_OVERRIDE {
      allocReals.free();
    }

    void
    factorize(
      char const      who[],
      integer         N,
      valueType const L[],
      valueType const D[],
      valueType const U[]
    );

    void
    lsq(
      integer nrhs,
      T       RHS[],
      integer ldRHS,
      T       lambda
    ) const;

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    void
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    solve(
      integer   nrhs,
      valueType xb[],
      integer   ldXB
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve(
      integer   nrhs,
      valueType xb[],
      integer   ldXB
    ) const LAPACK_WRAPPER_OVERRIDE;

    /*\
    :|:     _
    :|:    / \  _   ___  __
    :|:   / _ \| | | \ \/ /
    :|:  / ___ \ |_| |>  <
    :|: /_/   \_\__,_/_/\_\
    :|:
    \*/

    void
    axpy(
      integer         N,
      valueType       alpha,
      valueType const L[],
      valueType const D[],
      valueType const U[],
      valueType const x[],
      valueType       beta,
      valueType       y[]
    ) const;

  };

  //============================================================================
  /*\
  :|:  ___ _         _     _____    _    _ _                         _
  :|: | _ ) |___  __| |__ |_   _| _(_)__| (_)__ _ __ _ ___ _ _  __ _| |
  :|: | _ \ / _ \/ _| / /   | || '_| / _` | / _` / _` / _ \ ' \/ _` | |
  :|: |___/_\___/\__|_\_\   |_||_| |_\__,_|_\__,_\__, \___/_||_\__,_|_|
  :|:                                            |___/
  :|:  ___                _       _
  :|: / __|_  _ _ __  ___| |_ _ _(_)__
  :|: \__ \ || | '  \/ -_)  _| '_| / _|
  :|: |___/\_, |_|_|_\___|\__|_| |_\__|
  :|:      |__/
  \*/

  template <typename T>
  class BlockTridiagonalSymmetic : public LinearSystemSolver<T> {
  public:
    typedef T valueType;

  private:

    Malloc<valueType>  allocReals;
    Malloc<integer>    allocIntegers;
    Malloc<valueType*> allocRpointers;
    Malloc<integer*>   allocIpointers;

    integer      nBlocks, nnz;
    valueType ** D_blocks;
    valueType ** L_blocks;
    valueType *  Work;
    integer   ** B_permutation;
    integer   *  row_blocks;
    bool         is_factorized;

    void
    find(
      integer   ii,
      integer   jj,
      integer & iBlock,
      integer & jBlock,
      integer & ij,
      integer & ji
    ) const ;

  public:

    BlockTridiagonalSymmetic()
    : allocReals("BlockTridiagonalSymmetic-allocReals")
    , allocIntegers("BlockTridiagonalSymmetic-allocIntegers")
    , allocRpointers("BlockTridiagonalSymmetic-allocRpointers")
    , allocIpointers("BlockTridiagonalSymmetic-allocIpointers")
    , nBlocks(0)
    , nnz(0)
    , is_factorized(false)
    {}

    virtual
    ~BlockTridiagonalSymmetic() LAPACK_WRAPPER_OVERRIDE {
      allocReals.free();
      allocIntegers.free();
      allocRpointers.free();
      allocIpointers.free();
    }

    void
    setup( integer nblks, integer const rBlocks[] );

    void
    setup( integer nblks, integer const block_size );

    void
    zero();

    integer numBlocks() const { return nBlocks; }

    integer DnumRows( integer n ) const { return row_blocks[n+1] - row_blocks[n]; }
    integer DnumCols( integer n ) const { return row_blocks[n+1] - row_blocks[n]; }

    integer LnumRows( integer n ) const { return DnumRows(n+1); }
    integer LnumCols( integer n ) const { return DnumCols(n); }

    void
    setD(
      integer         n,
      valueType const data[],
      integer         ldData,
      bool            transposed=false
    );

    void
    setL(
      integer         n,
      valueType const data[],
      integer         ldData,
      bool            transposed=false
    );

    void
    setD(
      integer         n,
      valueType const data[],
      integer         ldData,
      integer         beginRow,
      integer         beginCol,
      integer         nrow,
      integer         ncol,
      bool            transposed=false
    );

    void
    setL(
      integer         n,
      valueType const data[],
      integer         ldData,
      integer         beginRow,
      integer         beginCol,
      integer         nrow,
      integer         ncol,
      bool            transposed=false
    );

    void
    insert( integer i, integer j, valueType v, bool sym );

    valueType const & operator () ( integer ii, integer jj ) const;
    valueType       & operator () ( integer ii, integer jj );

    void
    factorize( char const who[] );

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    void
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    solve(
      integer   nrhs,
      valueType xb[],
      integer   ldXB
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve(
      integer   nrhs,
      valueType xb[],
      integer   ldXB
    ) const LAPACK_WRAPPER_OVERRIDE;

  };

  //============================================================================
  /*\
  :|:   ____                  _          _ _    _   _
  :|:  | __ )  __ _ _ __   __| | ___  __| | |  | | | |
  :|:  |  _ \ / _` | '_ \ / _` |/ _ \/ _` | |  | | | |
  :|:  | |_) | (_| | | | | (_| |  __/ (_| | |__| |_| |
  :|:  |____/ \__,_|_| |_|\__,_|\___|\__,_|_____\___/
  \*/
  //! base class for linear systema solver
  template <typename T>
  class BandedLU : public LinearSystemSolver<T> {
  public:
    typedef T valueType;

    Malloc<valueType> allocReals;
    Malloc<integer>   allocIntegers;

    integer     m, n, nL, nU, ldAB;
    integer   * ipiv;
    valueType * AB;

    bool is_factorized;

  public:

    BandedLU();
    virtual ~BandedLU() LAPACK_WRAPPER_OVERRIDE;

    void
    setup(
      integer M,  // number of rows
      integer N,  // number of columns
      integer nL, // number of lower diagonal
      integer nU  // number of upper diagonal
    );

    integer
    iaddr( integer i, integer j ) const {
      integer d = (i-j+nL+nU);
      return d+j*ldAB;
    }

    void
    check( integer i, integer j ) const;

    valueType const &
    operator () ( integer i, integer j ) const
    { return AB[iaddr(i,j)]; }

    valueType &
    operator () ( integer i, integer j )
    { return AB[iaddr(i,j)]; }

    void
    insert( integer i, integer j, valueType v, bool sym );

    void zero();

    void
    load_block(
      integer         nr,
      integer         nc,
      valueType const B[],
      integer         ldB,
      integer         irow,
      integer         icol
    );

    // do internal factorization, to be executed (only once) before to call solve or t_solve
    void factorize( char const who[] );

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    void
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const LAPACK_WRAPPER_OVERRIDE;

    /*\
    :|:     _
    :|:    / \  _   ___  __
    :|:   / _ \| | | \ \/ /
    :|:  / ___ \ |_| |>  <
    :|: /_/   \_\__,_/_/\_\
    :|:
    \*/

    // y <- beta*y + alpha*A*x
    void
    aAxpy(
      valueType       alpha,
      valueType const x[],
      valueType       y[]
    ) const;

    void
    dump( ostream_type & stream ) const;

  };

  //============================================================================
  /*\
  :|:   ____                  _          _ ____  ____  ____
  :|:  | __ )  __ _ _ __   __| | ___  __| / ___||  _ \|  _ \
  :|:  |  _ \ / _` | '_ \ / _` |/ _ \/ _` \___ \| |_) | | | |
  :|:  | |_) | (_| | | | | (_| |  __/ (_| |___) |  __/| |_| |
  :|:  |____/ \__,_|_| |_|\__,_|\___|\__,_|____/|_|   |____/
  \*/
  //! base class for linear systema solver
  template <typename T>
  class BandedSPD : public LinearSystemSolver<T> {
  public:
    typedef T valueType;

    Malloc<valueType> allocReals;

    integer     n, nD, ldAB;
    valueType * AB;
    ULselect    UPLO;
    bool is_factorized;

  public:

    BandedSPD();

    virtual
    ~BandedSPD() LAPACK_WRAPPER_OVERRIDE;

    void
    setup(
      ULselect UPLO,
      integer  N,    // number of rows and columns
      integer  nD    // number of upper diagonal
    );

    valueType const &
    operator () ( integer i, integer j ) const
    { return AB[i+j*ldAB]; }

    valueType &
    operator () ( integer i, integer j )
    { return AB[i+j*ldAB]; }

    void
    insert( integer i, integer j, valueType v, bool sym );

    void zero();

    // do internal fatcorization, to be executed (only once) before to call solve or t_solve
    void factorize( char const who[] );

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    void
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const LAPACK_WRAPPER_OVERRIDE;

    /*\
    :|:     _
    :|:    / \  _   ___  __
    :|:   / _ \| | | \ \/ /
    :|:  / ___ \ |_| |>  <
    :|: /_/   \_\__,_/_/\_\
    :|:
    \*/

    /* not yet available
    // y <- beta*y + alpha*A*x
    void
    aAxpy( valueType       alpha,
           valueType const x[],
           valueType       y[] ) const;

    void
    dump( ostream_type & stream ) const;
    */

  };

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

  /*\
  :|:    ___                  _ _   _               _
  :|:   / _ \ _   _  __ _ ___(_) \ | | _____      _| |_ ___  _ __
  :|:  | | | | | | |/ _` / __| |  \| |/ _ \ \ /\ / / __/ _ \| '_ \
  :|:  | |_| | |_| | (_| \__ \ | |\  |  __/\ V  V /| || (_) | | | |
  :|:   \__\_\\__,_|\__,_|___/_|_| \_|\___| \_/\_/  \__\___/|_| |_|
  \*/

  template <typename T>
  class QN {
  public:
    typedef T valueType;

  protected:
    Malloc<valueType> allocReals;

    integer     n;
    valueType * H;
    valueType * s;
    valueType * y;
    valueType * z;

  public:

    QN()
    : allocReals("QN reals")
    , n(0)
    , H(nullptr)
    , s(nullptr)
    , y(nullptr)
    , z(nullptr)
    {}

    virtual
    ~QN() {}

    void
    allocate( integer N );

    valueType const &
    operator () ( integer i, integer j ) const
    { return H[i+j*n]; }

    valueType &
    operator () ( integer i, integer j )
    { return H[i+j*n]; }

    void
    zero()
    { lapack_wrapper::gezero( n, n, H, n ); }

    void
    init()
    { lapack_wrapper::geid( n, n, H, n ); }

    // r <- beta * r + alpha * H * x
    void
    mult(
      valueType       alpha,
      valueType const x[],
      integer         inc_x,
      valueType       beta,
      valueType       r[],
      integer         inc_r
    ) const {
      symv(
        LOWER,
        n,
        alpha, H, n,
        x, inc_x,
        beta,
        r, inc_r
      );
    }

    // r <- H * x
    void
    mult( valueType const x[], valueType r[] ) const {
      mult( valueType(1), x, 1, valueType(0), r, 1 );
    }

    void
    update(
      valueType const f0[],
      valueType const f1[],
      valueType const ss[]
    ) {
      copy( n,     f1, 1, y, 1 );
      axpy( n, -1, f0, 1, y, 1 );
      update( y, ss );
    }

    void
    update(
      valueType const f0[],
      valueType const f1[],
      valueType const x0[],
      valueType const x1[]
    ) {
      copy( n,     f1, 1, y, 1 );
      axpy( n, -1, f0, 1, y, 1 );
      copy( n,     x1, 1, s, 1 );
      axpy( n, -1, x0, 1, s, 1 );
      update( y, s );
    }

    void
    print( ostream_type & stream ) const;

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    void
    update(
      valueType const y[],
      valueType const s[]
    ) LAPACK_WRAPPER_PURE_VIRTUAL;

  };

  /*\
  :|:  ____  _____ ____ ____
  :|: | __ )|  ___/ ___/ ___|
  :|: |  _ \| |_ | |  _\___ \
  :|: | |_) |  _|| |_| |___) |
  :|: |____/|_|   \____|____/
  \*/

  template <typename T>
  class BFGS : public QN<T> {
  public:
    typedef T valueType;

    using QN<T>::allocate;
    using QN<T>::zero;
    using QN<T>::init;
    using QN<T>::mult;
    using QN<T>::update;

    using QN<T>::n;
    using QN<T>::H;
    using QN<T>::z;

    BFGS() {}

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    void
    update(
      valueType const y[],
      valueType const s[]
    ) LAPACK_WRAPPER_OVERRIDE;

  };

  /*\
  :|:   ____  _____ ____
  :|:  |  _ \|  ___|  _ \
  :|:  | | | | |_  | |_) |
  :|:  | |_| |  _| |  __/
  :|:  |____/|_|   |_|
  \*/

  template <typename T>
  class DFP : public QN<T> {
  public:
    typedef T valueType;

    using QN<T>::allocate;
    using QN<T>::zero;
    using QN<T>::init;
    using QN<T>::mult;
    using QN<T>::update;

    using QN<T>::n;
    using QN<T>::H;
    using QN<T>::z;

    DFP() {}

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    void
    update(
      valueType const y[],
      valueType const s[]
    ) LAPACK_WRAPPER_OVERRIDE;
  };

  /*\
  :|:   _____ _                            _
  :|:  | ____(_) __ _  ___ _ ____   ____ _| |_   _  ___  ___
  :|:  |  _| | |/ _` |/ _ \ '_ \ \ / / _` | | | | |/ _ \/ __|
  :|:  | |___| | (_| |  __/ | | \ V / (_| | | |_| |  __/\__ \
  :|:  |_____|_|\__, |\___|_| |_|\_/ \__,_|_|\__,_|\___||___/
  :|:           |___/
  \*/

  //! Sparse Matrix Structure
  template <typename T>
  class Eigenvalues {
  public:
    typedef T                valueType;
    typedef std::complex<T>  complexType;
    typedef MatrixWrapper<T> MatW;
    typedef SparseCCOOR<T>   Sparse;

  private:
    Malloc<valueType> mem_real;

    integer     N, Lwork;
    valueType * Re;
    valueType * Im;
    valueType * Work;
    valueType * A_saved;

    void allocate( integer N );
    void compute( );

  public:

    Eigenvalues();
    Eigenvalues( integer NRC, valueType const data[], integer ldData );
    Eigenvalues( MatW const & M );
    Eigenvalues(
      integer         NRC,
      integer         nnz,
      valueType const values[],
      integer   const row[],
      integer   const col[]
    );

    void setup( integer NRC, valueType const data[], integer ldData );
    void setup( MatW const & M );

    void
    setup(
      integer         NRC,
      integer         nnz,
      valueType const values[],
      integer   const row[],
      integer   const col[]
    );

    void getEigenvalue( integer n, valueType & re, valueType & im ) const;
    void getEigenvalue( integer n, complexType & eig ) const;

    void getEigenvalues( std::vector<valueType> & re, std::vector<valueType> & im ) const;
    void getEigenvalues( std::vector<complexType> & eigs ) const;

  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  class Eigenvectors {
  public:
    typedef T                valueType;
    typedef std::complex<T>  complexType;
    typedef MatrixWrapper<T> MatW;
    typedef SparseCCOOR<T>   Sparse;

  private:
    Malloc<valueType> mem_real;

    integer     N, Lwork;
    valueType * Re;
    valueType * Im;
    valueType * A_saved;
    valueType * VL;
    valueType * VR;
    valueType * Work;

    void allocate( integer N );
    void compute( );

  public:

    Eigenvectors();

    Eigenvectors( integer NRC, valueType const A[], integer ldA );

    Eigenvectors( MatW const & A );

    Eigenvectors(
      integer         NRC,
      integer         A_nnz,
      valueType const A_values[],
      integer   const A_row[],
      integer   const A_col[]
    );

    void
    setup( integer NRC, valueType const A[], integer ldA );

    void setup( MatW const & A );

    void
    setup(
      integer         NRC,
      integer         A_nnz,
      valueType const A_values[],
      integer   const A_row[],
      integer   const A_col[]
    );

    void getEigenvalue( integer n, valueType & re, valueType & im ) const;
    void getEigenvalue( integer n, complexType & eig ) const;

    void getEigenvalues( std::vector<valueType> & re, std::vector<valueType> & im ) const;
    void getEigenvalues( std::vector<complexType> & eigs ) const;

    void getLeftEigenvector( std::vector<std::vector<complexType> > & vecs ) const;
    void getRightEigenvector( std::vector<std::vector<complexType> > & vecs ) const;
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  class GeneralizedEigenvalues {
  public:
    typedef T                valueType;
    typedef std::complex<T>  complexType;
    typedef MatrixWrapper<T> MatW;
    typedef SparseCCOOR<T>   Sparse;

  private:
    Malloc<valueType> mem_real;

    integer     N, Lwork;
    valueType * alphaRe;
    valueType * alphaIm;
    valueType * beta;
    valueType * Work;
    valueType * A_saved;
    valueType * B_saved;

    void allocate( integer N );
    void compute( );

  public:

    GeneralizedEigenvalues();

    GeneralizedEigenvalues(
      integer NRC,
      valueType const A[], integer ldA,
      valueType const B[], integer ldB
    );

    GeneralizedEigenvalues( MatW const & A, MatW const & B );

    GeneralizedEigenvalues(
      integer         NRC,
      integer         A_nnz,
      valueType const A_values[],
      integer   const A_row[],
      integer   const A_col[],
      integer         B_nnz,
      valueType const B_values[],
      integer   const B_row[],
      integer   const B_col[]
    );

    void
    setup(
      integer NRC,
      valueType const A[], integer ldA,
      valueType const B[], integer ldB
    );

    void setup( MatW const & A, MatW const & B );

    void
    setup(
      integer         NRC,
      integer         A_nnz,
      valueType const A_values[],
      integer   const A_row[],
      integer   const A_col[],
      integer         B_nnz,
      valueType const B_values[],
      integer   const B_row[],
      integer   const B_col[]
    );

    void getEigenvalue( integer n, valueType & re, valueType & im ) const;
    void getEigenvalue( integer n, complexType & eig ) const;

    void getEigenvalues( std::vector<valueType> & re, std::vector<valueType> & im ) const;
    void getEigenvalues( std::vector<complexType> & eigs ) const;

  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  class GeneralizedEigenvectors {
  public:
    typedef T                valueType;
    typedef std::complex<T>  complexType;
    typedef MatrixWrapper<T> MatW;
    typedef SparseCCOOR<T>   Sparse;

  private:
    Malloc<valueType> mem_real;
    Malloc<integer>   mem_int;

    integer     N, Lwork, ilo, ihi;
    valueType   abnorm, bbnorm;
    valueType * alphaRe;
    valueType * alphaIm;
    valueType * beta;
    valueType * A_saved;
    valueType * B_saved;
    valueType * VL;
    valueType * VR;
    valueType * lscale;
    valueType * rscale;
    valueType * rconde;
    valueType * rcondv;
    valueType * Work;
    integer   * iWork;
    integer   * bWork;

    void allocate( integer N );
    void compute( );

  public:

    GeneralizedEigenvectors();

    GeneralizedEigenvectors(
      integer NRC,
      valueType const A[], integer ldA,
      valueType const B[], integer ldB
    );

    GeneralizedEigenvectors( MatW const & A, MatW const & B );

    GeneralizedEigenvectors(
      integer         NRC,
      integer         A_nnz,
      valueType const A_values[],
      integer   const A_row[],
      integer   const A_col[],
      integer         B_nnz,
      valueType const B_values[],
      integer   const B_row[],
      integer   const B_col[]
    );

    void
    setup(
      integer NRC,
      valueType const A[], integer ldA,
      valueType const B[], integer ldB
    );

    void setup( MatW const & A, MatW const & B );

    void
    setup(
      integer         NRC,
      integer         A_nnz,
      valueType const A_values[],
      integer   const A_row[],
      integer   const A_col[],
      integer         B_nnz,
      valueType const B_values[],
      integer   const B_row[],
      integer   const B_col[]
    );

    void getEigenvalue( integer n, valueType & re, valueType & im ) const;
    void getEigenvalue( integer n, complexType & eig ) const;

    void getEigenvalues( std::vector<valueType> & re, std::vector<valueType> & im ) const;
    void getEigenvalues( std::vector<complexType> & eigs ) const;

    void getLeftEigenvector( std::vector<std::vector<complexType> > & vecs ) const;
    void getRightEigenvector( std::vector<std::vector<complexType> > & vecs ) const;

    valueType balancedAnorm1() const { return this->abnorm; }
    valueType balancedBnorm1() const { return this->bbnorm; }
    valueType const * getLscale() const { return this->lscale; }
    valueType const * getRscale() const { return this->rscale; }
    valueType const * RcondEigenvalues()  const { return this->rconde; }
    valueType const * RcondEigenvectors() const { return this->rcondv; }
  };

  /*\
  :|:
  :|:    ____                           _ _             _ ______     ______
  :|:   / ___| ___ _ __   ___ _ __ __ _| (_)_______  __| / ___\ \   / /  _ \
  :|:  | |  _ / _ \ '_ \ / _ \ '__/ _` | | |_  / _ \/ _` \___ \\ \ / /| | | |
  :|:  | |_| |  __/ | | |  __/ | | (_| | | |/ /  __/ (_| |___) |\ V / | |_| |
  :|:   \____|\___|_| |_|\___|_|  \__,_|_|_/___\___|\__,_|____/  \_/  |____/
  :|:
  \*/

  template <typename T>
  class GeneralizedSVD {

    typedef T                valueType;
    typedef MatrixWrapper<T> MatW;
    typedef SparseCCOOR<T>   Sparse;

    Malloc<valueType> mem_real;
    Malloc<integer>   mem_int;

    integer     M, N, P, K, L, Lwork;
    valueType * Work;
    integer   * IWork;
    valueType * alpha_saved;
    valueType * beta_saved;
    valueType * A_saved;
    valueType * B_saved;
    valueType * U_saved;
    valueType * V_saved;
    valueType * Q_saved;

    MatrixWrapper<T>     U, V, Q, R;
    DiagMatrixWrapper<T> Dalpha, Dbeta;

    void allocate( integer N, integer M, integer P );
    void compute( );

    // R is stored in A(1:K+L,N-K-L+1:N) on exit.
    valueType &
    A( integer i, integer j )
    { return this->A_saved[i+j*this->M]; }

    valueType &
    B( integer i, integer j )
    { return this->B_saved[i+j*this->P]; }

  public:

    GeneralizedSVD();

    GeneralizedSVD(
      integer         m,
      integer         n,
      integer         p,
      valueType const A[], integer ldA,
      valueType const B[], integer ldB
    );

    GeneralizedSVD( MatW const & A, MatW const & B );

    GeneralizedSVD(
      integer         m,
      integer         n,
      integer         p,
      integer         A_nnz,
      valueType const A_values[],
      integer   const A_row[],
      integer   const A_col[],
      integer         B_nnz,
      valueType const B_values[],
      integer   const B_row[],
      integer   const B_col[]
    );

    void
    setup(
      integer         m,
      integer         n,
      integer         p,
      valueType const A[], integer ldA,
      valueType const B[], integer ldB
    );

    void setup( MatW const & A, MatW const & B );

    void
    setup(
      integer         m,
      integer         n,
      integer         p,
      integer         A_nnz,
      valueType const A_values[],
      integer   const A_row[],
      integer   const A_col[],
      integer         B_nnz,
      valueType const B_values[],
      integer   const B_row[],
      integer   const B_col[]
    );

    MatrixWrapper<T> const & getU() const { return this->U; }
    MatrixWrapper<T> const & getV() const { return this->V; }
    MatrixWrapper<T> const & getQ() const { return this->Q; }
    MatrixWrapper<T> const & getR() const { return this->R; }
    DiagMatrixWrapper<T> const & getC() const { return this->Dalpha; }
    DiagMatrixWrapper<T> const & gerS() const { return this->Dbeta; }

    valueType const & alpha( integer i ) const { return alpha_saved[i]; }
    valueType const & beta( integer i ) const { return beta_saved[i]; }

    valueType const * getAlpha() const { return alpha_saved; }
    valueType const * getBeta()  const { return beta_saved; }

    void info( ostream_type & stream, valueType eps = 0 ) const;

  };

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
