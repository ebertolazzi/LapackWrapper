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
 |      Università degli Studi di Trento                                    |
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
    using real_type = typename LinearSystemSolver<T>::real_type;
    using SVD_USED  = enum { USE_GESVD = 0, USE_GESDD = 1 };

  protected:

    integer     m_nrows{0};
    integer     m_ncols{0};

    real_type * m_Afactorized{nullptr};
    real_type * m_WorkSVD{nullptr};
    real_type * m_Umat{nullptr};
    real_type * m_VTmat{nullptr};
    real_type * m_Svec{nullptr};
    integer   * m_IWorkSVD{nullptr};

    real_type   m_rcond{Utils::machine_eps<real_type>()};

    integer     m_minRC{0};
    integer     m_LworkSVD{0};

    SVD_USED    m_svd_used;

  public:

    using LinearSystemSolver<T>::factorize;
    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;

    explicit
    SVD_no_alloc( SVD_USED _svd_used = USE_GESVD )
    : LinearSystemSolver<T>()
    , m_svd_used(_svd_used)
    {}

    static integer get_Lwork( integer NR, integer NC );

    void
    no_allocate(
      integer     NR,
      integer     NC,
      integer     Lwork,
      real_type * Work,
      integer     Liwork,
      integer   * iWork
    );

    //!
    //! Do SVD factorization of a rectangular matrix.
    //!
    //! \param who string used in error message
    //! \param A   pointer to the matrix
    //! \param LDA Leading dimension of the matrix
    //!
    virtual
    void
    factorize_nodim(
      string_view     who,
      real_type const A[],
      integer         LDA
    );

    virtual
    bool
    factorize_nodim( real_type const A[], integer LDA );

    void
    setRcond( real_type r )
    { m_rcond = r; }

    real_type U    ( integer i, integer j ) const { return m_Umat[i+j*m_nrows]; }
    real_type V    ( integer i, integer j ) const { return m_VTmat[j+i*m_ncols]; }
    real_type sigma( integer i )            const { return m_Svec[i]; }

    //! `y <- alpha * U * x + beta * y`
    void
    U_mul(
      real_type       alpha,
      real_type const x[],
      integer         incx,
      real_type       beta,
      real_type       y[],
      integer incy
    ) const {
      gemv(
        Transposition::NO,
        m_nrows, m_minRC,
        alpha, m_Umat, m_nrows,
        x, incx,
        beta, y, incy
      );
    }

    //! `y <- alpha * U' * x + beta * y`
    void
    Ut_mul(
      real_type       alpha,
      real_type const x[],
      integer         incx,
      real_type       beta,
      real_type       y[],
      integer         incy
    ) const {
      gemv(
        Transposition::YES,
        m_nrows, m_minRC,
        alpha, m_Umat, m_nrows,
        x, incx,
        beta, y, incy
      );
    }

    //! `y <- alpha * V * x + beta * y`
    void
    V_mul(
      real_type       alpha,
      real_type const x[],
      integer         incx,
      real_type       beta,
      real_type       y[],
      integer         incy
    ) const {
      gemv(
        Transposition::YES,
        m_minRC, m_ncols,
        alpha, m_VTmat, m_nrows,
        x, incx,
        beta, y, incy
      );
    }

    //! `y <- alpha * V' * x + beta * y`
    void
    Vt_mul(
      real_type       alpha,
      real_type const x[],
      integer         incx,
      real_type       beta,
      real_type       y[],
      integer         incy
    ) const {
      gemv(
        Transposition::NO,
        m_minRC, m_ncols,
        alpha, m_VTmat, m_nrows,
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

    bool solve( real_type xb[] ) const override;
    bool t_solve( real_type xb[] ) const override;

  };

  //============================================================================

  template <typename T>
  class SVD : public SVD_no_alloc<T> {
  public:
    using real_type = typename SVD_no_alloc<T>::real_type;
    using SVD_USED  = typename SVD_no_alloc<T>::SVD_USED;

  protected:

    Malloc<real_type> m_allocReals{"GeneralizedSVD(real)"};
    Malloc<integer>   m_allocIntegers{"GeneralizedSVD(int)"};

  public:

    using SVD_no_alloc<T>::m_nrows;
    using SVD_no_alloc<T>::m_ncols;
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

    explicit
    SVD( SVD_USED _svd_used = SVD_no_alloc<T>::USE_GESVD )
    : SVD_no_alloc<T>(_svd_used) {}

    void
    allocate( integer NR, integer NC );

    //!
    //! Do SVD factorization of a rectangular matrix.
    //!
    //! \param who string used in error message
    //! \param NR  number of rows of the matrix
    //! \param NC  number of columns of the matrix
    //! \param A   pointer to the matrix
    //! \param LDA Leading dimension of the matrix
    //!
    void
    factorize(
      string_view     who,
      integer         NR,
      integer         NC,
      real_type const A[],
      integer         LDA
    ) override {
      this->allocate( NR, NC );
      this->factorize_nodim( who, A, LDA );
    }

    bool
    factorize(
      integer         NR,
      integer         NC,
      real_type const A[],
      integer         LDA
    ) override {
      this->allocate( NR, NC );
      return this->factorize_nodim( A, LDA );
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

    using real_type = T;
    using MatW      = MatrixWrapper<T>;
    using Sparse    = SparseCCOOR<T>;

    Malloc<real_type> m_mem_real{"GeneralizedSVD(real)"};
    Malloc<integer>   m_mem_int{"GeneralizedSVD(int)"};

    integer     m_M{0};
    integer     m_N{0};
    integer     m_P{0};
    integer     m_K{0};
    integer     m_L{0};
    integer     m_Lwork{0};
    real_type * m_Work{nullptr};
    integer   * m_IWork{nullptr};
    real_type * m_alpha_saved{nullptr};
    real_type * m_beta_saved{nullptr};
    real_type * m_A_saved{nullptr};
    real_type * m_B_saved{nullptr};
    real_type * m_U_saved{nullptr};
    real_type * m_V_saved{nullptr};
    real_type * m_Q_saved{nullptr};

    MatrixWrapper<T>     m_U, m_V, m_Q, m_R;
    DiagMatrixWrapper<T> m_Dalpha, m_Dbeta;

    void allocate( integer N, integer M, integer P );
    void compute( );

    // R is stored in A(1:K+L,N-K-L+1:N) on exit.
    real_type &
    A( integer i, integer j )
    { return m_A_saved[i+j*m_M]; }

    real_type &
    B( integer i, integer j )
    { return m_B_saved[i+j*m_P]; }

  public:

    GeneralizedSVD() = default;

    GeneralizedSVD(
      integer         m,
      integer         n,
      integer         p,
      real_type const A[], integer ldA,
      real_type const B[], integer ldB
    ) {
      this->setup( m, n, p, A, ldA, B, ldB );
    }

    GeneralizedSVD( MatW const & A, MatW const & B ) {
      this->setup( A, B );
    }

    GeneralizedSVD(
      integer         m,
      integer         n,
      integer         p,
      integer         A_nnz,
      real_type const A_values[],
      integer   const A_row[],
      integer   const A_col[],
      integer         B_nnz,
      real_type const B_values[],
      integer   const B_row[],
      integer   const B_col[]
    ) {
      this->setup(
        m, n, p,
        A_nnz, A_values, A_row, A_col,
        B_nnz, B_values, B_row, B_col
      );
    }

    void
    setup(
      integer         m,
      integer         n,
      integer         p,
      real_type const A[], integer ldA,
      real_type const B[], integer ldB
    );

    void setup( MatW const & A, MatW const & B );

    void
    setup(
      integer         m,
      integer         n,
      integer         p,
      integer         A_nnz,
      real_type const A_values[],
      integer   const A_row[],
      integer   const A_col[],
      integer         B_nnz,
      real_type const B_values[],
      integer   const B_row[],
      integer   const B_col[]
    );

    MatrixWrapper<T>     const & getU() const { return m_U; }
    MatrixWrapper<T>     const & getV() const { return m_V; }
    MatrixWrapper<T>     const & getQ() const { return m_Q; }
    MatrixWrapper<T>     const & getR() const { return m_R; }
    DiagMatrixWrapper<T> const & getC() const { return m_Dalpha; }
    DiagMatrixWrapper<T> const & gerS() const { return m_Dbeta; }

    real_type const & alpha( integer i ) const { return m_alpha_saved[i]; }
    real_type const & beta( integer i ) const { return m_beta_saved[i]; }

    real_type const * getAlpha() const { return m_alpha_saved; }
    real_type const * getBeta()  const { return m_beta_saved; }

    string info( real_type eps = 0 ) const;

  };

}

///
/// eof: svd.hxx
///
