/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2019                                                      |
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
/// file: svd.hxx
///

namespace lapack_wrapper {

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

}

///
/// eof: svd.hxx
///
