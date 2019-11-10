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
  class SVD_no_alloc : public LinearSystemSolver<T> {
  public:
    typedef typename LinearSystemSolver<T>::valueType valueType;
    typedef enum { USE_GESVD = 0, USE_GESDD = 1 } SVD_USED;

  protected:

    integer nRows;
    integer nCols;

    valueType * Afactorized;
    valueType * Work;
    valueType * Umat;
    valueType * VTmat;
    valueType * Svec;
    integer   * IWork;

    valueType rcond;

    integer minRC;
    integer Lwork;

    SVD_USED svd_used;

  public:

    using LinearSystemSolver<T>::factorize;
    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;

    SVD_no_alloc( SVD_USED _svd_used = USE_GESVD )
    : LinearSystemSolver<T>()
    , nRows(0)
    , nCols(0)
    , Afactorized(nullptr)
    , Work(nullptr)
    , Umat(nullptr)
    , VTmat(nullptr)
    , Svec(nullptr)
    , IWork(nullptr)
    , rcond(machineEps<valueType>())
    , minRC(0)
    , Lwork(0)
    , svd_used(_svd_used)
    {}

    virtual
    ~SVD_no_alloc() LAPACK_WRAPPER_OVERRIDE
    {}

    void
    no_allocate(
      integer     NR,
      integer     NC,
      integer     L,
      valueType * _Afactorized,
      valueType * _Work,
      valueType * _Umat,
      valueType * _VTmat,
      valueType * _Svec,
      integer   * _IWork
    ) {
      this->nRows       = NR;
      this->nCols       = NC;
      this->Lwork       = L;
      this->minRC       = std::min( NR, NC );
      this->Afactorized = _Afactorized;
      this->Work        = _Work;
      this->Umat        = _Umat;
      this->VTmat       = _VTmat;
      this->Svec        = _Svec;
      this->IWork       = _IWork;
    }

    /*!
     *  Do SVD factorization of a rectangular matrix
     *  \param A   pointer to the matrix
     *  \param LDA Leading dimension of the matrix
     */
    virtual
    void
    factorize(
      char const      who[],
      valueType const A[],
      integer         LDA
    );

    void
    setRcond( valueType r )
    { rcond = r; }

    valueType U    ( integer i, integer j ) const { return Umat[i+j*this->nRows]; }
    valueType V    ( integer i, integer j ) const { return VTmat[j+i*this->nCols]; }
    valueType sigma( integer i )            const { return Svec[i]; }

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
        this->nRows, this->minRC,
        alpha, this->Umat, this->nRows,
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
        this->nRows, this->minRC,
        alpha, this->Umat, this->nRows,
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
        this->minRC, this->nCols,
        alpha, this->VTmat, this->nRows,
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
        this->minRC, this->nCols,
        alpha, this->VTmat, this->nRows,
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
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

  };

  //============================================================================

  template <typename T>
  class SVD : public SVD_no_alloc<T> {
  public:
    typedef typename SVD_no_alloc<T>::valueType valueType;
    typedef typename SVD_no_alloc<T>::SVD_USED  SVD_USED;

  protected:

    Malloc<valueType> allocReals;
    Malloc<integer>   allocIntegers;

  public:

    using SVD_no_alloc<T>::solve;
    using SVD_no_alloc<T>::t_solve;
    using SVD_no_alloc<T>::factorize;
    using SVD_no_alloc<T>::setRcond;
    using SVD_no_alloc<T>::U;
    using SVD_no_alloc<T>::V;
    using SVD_no_alloc<T>::sigma;
    using SVD_no_alloc<T>::U_mul;
    using SVD_no_alloc<T>::Ut_mul;
    using SVD_no_alloc<T>::V_mul;
    using SVD_no_alloc<T>::Vt_mul;
    using SVD_no_alloc<T>::no_allocate;

    SVD( SVD_USED _svd_used = SVD_no_alloc<T>::USE_GESVD )
    : SVD_no_alloc<T>(_svd_used)
    , allocReals("SVD-allocReals")
    , allocIntegers("SVD-allocIntegers")
    {}

    virtual
    ~SVD() LAPACK_WRAPPER_OVERRIDE
    { allocReals.free(); allocIntegers.free(); }

    void
    allocate( integer NR, integer NC );

    /*!
     *  Do SVD factorization of a rectangular matrix
     *  \param NR  number of rows of the matrix
     *  \param NC  number of columns of the matrix
     *  \param A   pointer to the matrix
     *  \param LDA Leading dimension of the matrix
     */
    virtual
    void
    factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) LAPACK_WRAPPER_OVERRIDE {
      this->allocate( NR, NC );
      this->factorize( who, A, LDA );
    }

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
