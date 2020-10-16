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

    integer     m_nRows;
    integer     m_nCols;

    valueType * m_Afactorized;
    valueType * m_WorkSVD;
    valueType * m_Umat;
    valueType * m_VTmat;
    valueType * m_Svec;
    integer   * m_IWorkSVD;

    valueType   m_rcond;

    integer     m_minRC;
    integer     m_LworkSVD;

    SVD_USED    m_svd_used;

  public:

    using LinearSystemSolver<T>::factorize;
    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;

    SVD_no_alloc( SVD_USED _svd_used = USE_GESVD )
    : LinearSystemSolver<T>()
    , m_nRows(0)
    , m_nCols(0)
    , m_Afactorized(nullptr)
    , m_WorkSVD(nullptr)
    , m_Umat(nullptr)
    , m_VTmat(nullptr)
    , m_Svec(nullptr)
    , m_IWorkSVD(nullptr)
    , m_rcond(machineEps<valueType>())
    , m_minRC(0)
    , m_LworkSVD(0)
    , m_svd_used(_svd_used)
    {}

    virtual
    ~SVD_no_alloc() LAPACK_WRAPPER_OVERRIDE
    {}

    integer
    get_Lwork( integer NR, integer NC ) const;

    void
    no_allocate(
      integer     NR,
      integer     NC,
      integer     Lwork,
      valueType * Work,
      integer     Liwork,
      integer   * iWork
    );

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

    virtual
    bool
    factorize( valueType const A[], integer LDA );

    void
    setRcond( valueType r )
    { m_rcond = r; }

    valueType U    ( integer i, integer j ) const { return m_Umat[i+j*m_nRows]; }
    valueType V    ( integer i, integer j ) const { return m_VTmat[j+i*m_nCols]; }
    valueType sigma( integer i )            const { return m_Svec[i]; }

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
        m_nRows, m_minRC,
        alpha, m_Umat, m_nRows,
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
        m_nRows, m_minRC,
        alpha, m_Umat, m_nRows,
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
        m_minRC, m_nCols,
        alpha, m_VTmat, m_nRows,
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
        m_minRC, m_nCols,
        alpha, m_VTmat, m_nRows,
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
    bool
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    bool
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

  };

  //============================================================================

  template <typename T>
  class SVD : public SVD_no_alloc<T> {
  public:
    typedef typename SVD_no_alloc<T>::valueType valueType;
    typedef typename SVD_no_alloc<T>::SVD_USED  SVD_USED;

  protected:

    Malloc<valueType> m_allocReals;
    Malloc<integer>   m_allocIntegers;

  public:

    using SVD_no_alloc<T>::m_nRows;
    using SVD_no_alloc<T>::m_nCols;
    using SVD_no_alloc<T>::m_Afactorized;
    using SVD_no_alloc<T>::m_WorkSVD;
    using SVD_no_alloc<T>::m_Umat;
    using SVD_no_alloc<T>::m_VTmat;
    using SVD_no_alloc<T>::m_Svec;
    using SVD_no_alloc<T>::m_IWorkSVD;
    using SVD_no_alloc<T>::m_rcond;
    using SVD_no_alloc<T>::m_minRC;
    using SVD_no_alloc<T>::m_LworkSVD;
    using SVD_no_alloc<T>::m_svd_used;

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

    SVD( SVD_USED _svd_used = SVD_no_alloc<T>::USE_GESVD );

    virtual
    ~SVD() LAPACK_WRAPPER_OVERRIDE;

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

    virtual
    bool
    factorize(
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) LAPACK_WRAPPER_OVERRIDE {
      this->allocate( NR, NC );
      return this->factorize( A, LDA );
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

    Malloc<valueType> m_mem_real;
    Malloc<integer>   m_mem_int;

    integer     m_M;
    integer     m_N;
    integer     m_P;
    integer     m_K;
    integer     m_L;
    integer     m_Lwork;
    valueType * m_Work;
    integer   * m_IWork;
    valueType * m_alpha_saved;
    valueType * m_beta_saved;
    valueType * m_A_saved;
    valueType * m_B_saved;
    valueType * m_U_saved;
    valueType * m_V_saved;
    valueType * m_Q_saved;

    MatrixWrapper<T>     m_U, m_V, m_Q, m_R;
    DiagMatrixWrapper<T> m_Dalpha, m_Dbeta;

    void allocate( integer N, integer M, integer P );
    void compute( );

    // R is stored in A(1:K+L,N-K-L+1:N) on exit.
    valueType &
    A( integer i, integer j )
    { return m_A_saved[i+j*m_M]; }

    valueType &
    B( integer i, integer j )
    { return m_B_saved[i+j*m_P]; }

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

    MatrixWrapper<T>     const & getU() const { return m_U; }
    MatrixWrapper<T>     const & getV() const { return m_V; }
    MatrixWrapper<T>     const & getQ() const { return m_Q; }
    MatrixWrapper<T>     const & getR() const { return m_R; }
    DiagMatrixWrapper<T> const & getC() const { return m_Dalpha; }
    DiagMatrixWrapper<T> const & gerS() const { return m_Dbeta; }

    valueType const & alpha( integer i ) const { return m_alpha_saved[i]; }
    valueType const & beta( integer i ) const { return m_beta_saved[i]; }

    valueType const * getAlpha() const { return m_alpha_saved; }
    valueType const * getBeta()  const { return m_beta_saved; }

    void info( ostream_type & stream, valueType eps = 0 ) const;

  };

}

///
/// eof: svd.hxx
///
