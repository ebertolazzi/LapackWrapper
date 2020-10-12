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

    std::string m_name;
    size_t      m_numTotValues;
    size_t      m_numTotReserved;
    size_t      m_numAllocated;
    valueType * m_pMalloc;

    Malloc(Malloc<T> const &); // blocco costruttore di copia
    Malloc<T> const & operator = (Malloc<T> &) const; // blocco copia

  public:

    //! malloc object constructor
    explicit
    Malloc( std::string const & __name )
    : m_name(__name)
    , m_numTotValues(0)
    , m_numTotReserved(0)
    , m_numAllocated(0)
    , m_pMalloc(nullptr)
    { }

    //! malloc object destructor
    ~Malloc() { free(); }

    //! allocate memory for `n` objects
    void
    allocate( size_t n ) {
      try {
        if ( n > m_numTotReserved ) {
          delete [] m_pMalloc;
          m_numTotValues   = n;
          m_numTotReserved = n + (n>>3); // 12% more values
          m_pMalloc        = new T[m_numTotReserved];
        }
      }
      catch ( std::exception const & exc ) {
        std::string reason = fmt::format(
          "Memory allocation failed: {}\nTry to allocate {} bytes for {}\n",
          exc.what(), n, m_name
        );
        printTrace( __LINE__, __FILE__, reason, std::cerr );
        std::exit(0);
      }
      catch (...) {
        std::string reason = fmt::format(
          "Memory allocation failed for {}: memory exausted\n", m_name
        );
        printTrace( __LINE__, __FILE__, reason, std::cerr );
        std::exit(0);
      }
      m_numTotValues = n;
      m_numAllocated = 0;
    }

    //! free memory
    void
    free(void) {
      if ( m_pMalloc != nullptr ) {
        delete [] m_pMalloc; m_pMalloc = nullptr;
        m_numTotValues   = 0;
        m_numTotReserved = 0;
        m_numAllocated   = 0;
      }
    }

    //! number of objects allocated
    size_t size(void) const { return m_numTotValues; }

    //! get pointer of allocated memory for `sz` objets
    T * operator () ( size_t sz ) {
      size_t offs = m_numAllocated;
      m_numAllocated += sz;
      if ( m_numAllocated > m_numTotValues ) {
        std::string reason = fmt::format(
          "nMalloc<{}>::operator () ({}) -- Memory EXAUSTED\n", m_name, sz
        );
        printTrace( __LINE__, __FILE__, reason, std::cerr );
        std::exit(0);
      }
      return m_pMalloc + offs;
    }

    bool is_empty() const { return m_numAllocated <= m_numTotValues; }

    void
    must_be_empty( char const where[] ) const {
      if ( m_numAllocated < m_numTotValues ) {
        std::string tmp = fmt::format(
          "in {} {}: not fully used!\nUnused: {} values\n",
          m_name, where, m_numTotValues - m_numAllocated
        );
        printTrace( __LINE__,__FILE__, tmp, std::cerr );
      }
      if ( m_numAllocated > m_numTotValues ) {
        std::string tmp = fmt::format(
          "in {} {}: too much used!\nMore used: {} values\n",
          m_name, where, m_numAllocated - m_numTotValues
        );
        printTrace( __LINE__,__FILE__, tmp, std::cerr );
      }
    }

  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  inline
  std::string
  print_matrix(
    integer       nr,
    integer       nc,
    t_Value const A[],
    integer       ldA
  ) {
    std::string out;
    for ( integer i = 0; i < nr; ++i ) {
      for ( integer j = 0; j < nc; ++j )
        out += fmt::format("{:16.8} ",A[i+j*ldA]);
      out += '\n';
    }
    return out;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  inline
  std::string
  print_matrix( MatrixWrapper<t_Value> const & Amat ) {
    return print_matrix(
      Amat.numRows(), Amat.numCols(), Amat.get_data(), Amat.lDim()
    );
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
    lapack_wrapper::Malloc<T> m_mem;
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
    Malloc<T> m_mem;
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

  extern template class Matrix<real>;
  extern template class Matrix<doublereal>;
  extern template class DiagMatrix<real>;
  extern template class DiagMatrix<doublereal>;

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
    bool
    solve( valueType xb[] ) const LAPACK_WRAPPER_PURE_VIRTUAL;

    virtual
    bool
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_PURE_VIRTUAL;

    virtual
    bool
    solve( integer nrhs, valueType B[], integer ldB ) const {
      for ( integer i = 0; i < nrhs; ++i )
        if ( !solve( B + i*ldB ) ) return false;
      return true;
    }

    virtual
    bool
    t_solve( integer nrhs, valueType B[], integer ldB ) const {
      for ( integer i = 0; i < nrhs; ++i )
        if ( !t_solve( B + i*ldB ) ) return false;
      return true;
    }

    virtual
    void
    solve( char const who[], valueType xb[] ) const {
      LW_ASSERT(
        this->solve( xb ),
        "LinearSystemSolver::solve failed at {}\n", who
      );
    }

    virtual
    void
    t_solve( char const who[], valueType xb[] ) const {
      LW_ASSERT(
        this->t_solve( xb ),
        "LinearSystemSolver::t_solve failed at {}\n", who
      );
    }

    virtual
    void
    solve( char const who[], integer nrhs, valueType B[], integer ldB ) const {
      LW_ASSERT(
        this->solve( nrhs, B, ldB ),
        "LinearSystemSolver::solve failed at {}\n", who
      );
    }

    virtual
    void
    t_solve( char const who[], integer nrhs, valueType B[], integer ldB ) const {
      LW_ASSERT(
        this->t_solve( nrhs, B, ldB ),
        "LinearSystemSolver::solve failed at {}\n", who
      );
    }

    virtual
    void
    factorize(
      char const      /* who */ [],
      integer         /* NR  */   ,
      integer         /* NC  */   ,
      valueType const /* A   */ [],
      integer         /* LDA */
    ) {
      LW_ERROR0( "LinearSystemSolver::factorize, not defined in derived class\n" );
    }

    virtual
    bool
    factorize(
      integer         /* NR  */   ,
      integer         /* NC  */   ,
      valueType const /* A   */ [],
      integer         /* LDA */
    ) {
      return false;
    }

    // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    bool
    solve( MatrixWrapper<valueType> & M )
    { return solve( M.numCols(), M.get_data(),  M.lDim() ); }

    void
    solve( char const who[], MatrixWrapper<valueType> & M )
    { this->solve( who, M.numCols(), M.get_data(),  M.lDim() ); }

    bool
    t_solve( MatrixWrapper<valueType> & M )
    { return t_solve( M.numCols(), M.get_data(),  M.lDim() ); }

    void
    t_solve( char const who[], MatrixWrapper<valueType> & M )
    { this->t_solve( who, M.numCols(), M.get_data(),  M.lDim() ); }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    factorize( char const who[], MatrixWrapper<valueType> const & M ) {
      this->factorize( who, M.numRows(), M.numCols(), M.get_data(), M.lDim() );
    }

    bool
    factorize( MatrixWrapper<valueType> const & M ) {
      return this->factorize( M.numRows(), M.numCols(), M.get_data(), M.lDim() );
    }

    void
    t_factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) {
      Matrix<valueType> M(NC,NR);
      for ( integer i = 0; i < NR; ++i )
        for ( integer j = 0; j < NC; ++j )
          M(j,i) = A[i+j*LDA];
      this->factorize( who, M );
    }

    bool
    t_factorize(
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) {
      Matrix<valueType> M(NC,NR);
      for ( integer i = 0; i < NR; ++i )
        for ( integer j = 0; j < NC; ++j )
          M(j,i) = A[i+j*LDA];
      return this->factorize( M );
    }

    void
    t_factorize( char const who[], MatrixWrapper<valueType> const & M ) {
      this->t_factorize( who, M.numRows(), M.numCols(), M.get_data(), M.lDim() );
    }

    bool
    t_factorize( MatrixWrapper<valueType> const & M ) {
      return this->t_factorize( M.numRows(), M.numCols(), M.get_data(), M.lDim() );
    }

  };

}

#include "code++/lu.hxx"
#include "code++/qr.hxx"
#include "code++/svd.hxx"
#include "code++/ls.hxx"
#include "code++/lsc.hxx"
#include "code++/trid.hxx"
#include "code++/band.hxx"
#include "code++/qn.hxx"
#include "code++/eig.hxx"
#include "code++/block_trid.hxx"

namespace lapack_wrapper {

  // explicit instantiation declaration to suppress warnings

  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  #endif

  extern template class LU_no_alloc<real>;
  extern template class LU_no_alloc<doublereal>;
  extern template class LU<real>;
  extern template class LU<doublereal>;

  extern template class LUPQ_no_alloc<real>;
  extern template class LUPQ_no_alloc<doublereal>;
  extern template class LUPQ<real>;
  extern template class LUPQ<doublereal>;

  extern template class QR_no_alloc<real>;
  extern template class QR_no_alloc<doublereal>;
  extern template class QR<real>;
  extern template class QR<doublereal>;

  extern template class QRP_no_alloc<real>;
  extern template class QRP_no_alloc<doublereal>;
  extern template class QRP<real>;
  extern template class QRP<doublereal>;

  extern template class SVD_no_alloc<real>;
  extern template class SVD_no_alloc<doublereal>;
  extern template class SVD<real>;
  extern template class SVD<doublereal>;

  extern template class LSS_no_alloc<real>;
  extern template class LSS_no_alloc<doublereal>;
  extern template class LSS<real>;
  extern template class LSS<doublereal>;

  extern template class LSY_no_alloc<real>;
  extern template class LSY_no_alloc<doublereal>;
  extern template class LSY<real>;
  extern template class LSY<doublereal>;

  extern template class LSC<real>;
  extern template class LSC<doublereal>;

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
