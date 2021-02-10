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
/// file: qr.hxx
///

namespace lapack_wrapper {

  //============================================================================
  /*\
  :|:    ___  ____
  :|:   / _ \|  _ \
  :|:  | | | | |_) |
  :|:  | |_| |  _ <
  :|:   \__\_\_| \_\
  \*/

  template <typename T>
  class QR_no_alloc : public LinearSystemSolver<T> {
  public:
    typedef typename LinearSystemSolver<T>::valueType valueType;

  protected:

    integer m_nrows;
    integer m_ncols;
    integer m_nReflector;
    integer m_LworkQR;

    valueType * m_Afactorized;
    valueType * m_WorkQR;
    valueType * m_Tau;

  public:

    using LinearSystemSolver<T>::factorize;
    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;

    QR_no_alloc()
    : LinearSystemSolver<T>()
    , m_nrows(0)
    , m_ncols(0)
    , m_nReflector(0)
    , m_LworkQR(0)
    , m_Afactorized(nullptr)
    , m_WorkQR(nullptr)
    , m_Tau(nullptr)
    {}

    virtual
    ~QR_no_alloc() UTILS_OVERRIDE {}

    integer
    get_Lwork_QR( integer NR, integer NC ) const;

    void
    no_allocate(
      integer     NR,
      integer     NC,
      integer     Lwork,
      valueType * Work
    );

    /*!
     * Do QR factorization of a rectangular matrix
     * \param A   pointer to the matrix
     * \param LDA Leading dimension of the matrix
     */
    void
    factorize_nodim(
      char const      who[],
      valueType const A[],
      integer         LDA
    );

    bool
    factorize_nodim( valueType const A[], integer LDA );

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
    void Q_mul( valueType x[] ) const;

    //! x <- Q'*x
    void Qt_mul( valueType x[] ) const;

    //! C <- Q*C
    void Q_mul( integer nr, integer nc, valueType C[], integer ldC ) const;

    void
    Q_mul( MatrixWrapper<T> & C ) const {
      this->Q_mul( C.numRows(), C.numCols(), C.data(), C.lDim() );
    }

    //! C <- Q'*C
    void Qt_mul( integer nr, integer nc, valueType C[], integer ldC ) const;

    void
    Qt_mul( MatrixWrapper<T> & C ) const {
      this->Qt_mul( C.numRows(), C.numCols(), C.data(), C.lDim() );
    }

    //! C <- C*Q
    void mul_Q( integer nr, integer nc, valueType C[], integer ldC ) const;

    void
    mul_Q( MatrixWrapper<T> & C ) const {
      this->mul_Q( C.numRows(), C.numCols(), C.data(), C.lDim() );
    }

    //! C <- C*Q'
    void mul_Qt( integer nr, integer nc, valueType C[], integer ldC ) const;

    void
    mul_Qt( MatrixWrapper<T> & C ) const {
      this->mul_Q( C.numRows(), C.numCols(), C.data(), C.lDim() );
    }

    // -------------------------------------------------------------------------

    //! x <- R^(-1) * x
    void invR_mul( valueType x[], integer incx = 1 ) const;

    //! C <- R^(-1) * C
    void invR_mul( integer nr, integer nc, valueType C[], integer ldC ) const;

    void
    invR_mul( MatrixWrapper<valueType> & C ) const {
      this->invR_mul( C.numRows(), C.numCols(), C.data(), C.lDim() );
    }

    // -------------------------------------------------------------------------

    //! x <- R^(-T) * x
    void invRt_mul( valueType x[], integer incx = 1 ) const;

    //! C <- R^(-T) * C
    void invRt_mul( integer nr, integer nc, valueType C[], integer ldC ) const;

    void
    invRt_mul( MatrixWrapper<valueType> & C ) const {
      this->invRt_mul( C.numRows(), C.numCols(), C.data(), C.lDim() );
    }

    // -------------------------------------------------------------------------

    //! C <- C * R^(-1)
    void mul_invR( integer nr, integer nc, valueType C[], integer ldC ) const;

    void
    mul_invR( MatrixWrapper<valueType> & C ) const {
      this->mul_invR( C.numRows(), C.numCols(), C.data(), C.lDim() );
    }

    // -------------------------------------------------------------------------

    void mul_invRt( integer nr, integer nc, valueType C[], integer ldC ) const;

    void
    mul_invRt( MatrixWrapper<valueType> & C ) const {
      this->mul_invRt( C.numRows(), C.numCols(), C.data(), C.lDim() );
    }

    // -------------------------------------------------------------------------

    void getR( valueType R[], integer ldR ) const;
    void getR( Matrix<T> & R ) const;
    void getRt( valueType R[], integer ldR ) const;
    void getRt( Matrix<T> & R ) const;

    // -------------------------------------------------------------------------

    void getQ( valueType Q[], integer ldQ ) const;
    void getQ( Matrix<T> & Q ) const;

    // -------------------------------------------------------------------------

    void getQreduced( valueType Q[], integer ldQ ) const;
    void getQreduced( Matrix<T> & Q ) const;

    // -------------------------------------------------------------------------

    void getA( valueType A[], integer ldA ) const;
    void getA( Matrix<T> & A ) const;

    // -------------------------------------------------------------------------

    void getTau( valueType tau[] ) const;
    void getTau( Matrix<T> & tau ) const;

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    /*!
     *   In case of QR factorization of a square matrix solve the
     *   linear system \f$ QR x = b \f$
     * param xb on input the rhs of linear system on output the solution
     */
    virtual
    bool
    solve( valueType xb[] ) const UTILS_OVERRIDE;

    /*!
     *  In case of QR factorization of a square matrix solve the
     *  linear system \f$ (QR)^T x = b \f$
     *  \param xb on input the rhs of linear system on output the solution
     */
    virtual
    bool
    t_solve( valueType xb[] ) const UTILS_OVERRIDE;

    virtual
    bool
    solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const UTILS_OVERRIDE;

    virtual
    bool
    t_solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const UTILS_OVERRIDE;

  };

  //============================================================================
  //============================================================================
  //============================================================================

  template <typename T>
  class QR : public QR_no_alloc<T> {
  public:
    typedef typename QR_no_alloc<T>::valueType valueType;

  protected:

    Malloc<valueType> m_allocReals;

  public:

    using QR_no_alloc<T>::m_nrows;
    using QR_no_alloc<T>::m_ncols;
    using QR_no_alloc<T>::m_nReflector;
    using QR_no_alloc<T>::m_LworkQR;

    using QR_no_alloc<T>::m_Afactorized;
    using QR_no_alloc<T>::m_WorkQR;
    using QR_no_alloc<T>::m_Tau;

    using QR_no_alloc<T>::solve;
    using QR_no_alloc<T>::t_solve;
    using QR_no_alloc<T>::Q_mul;
    using QR_no_alloc<T>::Qt_mul;
    using QR_no_alloc<T>::mul_Q;
    using QR_no_alloc<T>::mul_Qt;
    using QR_no_alloc<T>::invR_mul;
    using QR_no_alloc<T>::invRt_mul;
    using QR_no_alloc<T>::getR;
    using QR_no_alloc<T>::getQ;
    using QR_no_alloc<T>::getA;
    using QR_no_alloc<T>::getTau;
    using QR_no_alloc<T>::get_Lwork_QR;
    using QR_no_alloc<T>::no_allocate;
    using QR_no_alloc<T>::factorize;

    QR()
    : QR_no_alloc<T>()
    , m_allocReals("QR-allocReals")
    {}

    virtual
    ~QR() UTILS_OVERRIDE
    { m_allocReals.free(); }

    void
    allocate( integer nr, integer nc );

    /*!
     * Do QR factorization of a rectangular matrix
     * \param NR  number of rows of the matrix
     * \param NC  number of columns of the matrix
     * \param A   pointer to the matrix
     * \param LDA Leading dimension of the matrix
     */
    void
    factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) UTILS_OVERRIDE {
      this->allocate( NR, NC );
      this->factorize_nodim( who, A, LDA );
    }

    bool
    factorize(
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) UTILS_OVERRIDE {
      this->allocate( NR, NC );
      return this->factorize_nodim( A, LDA );
    }

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
  class QRP_no_alloc : public QR_no_alloc<T> {
  public:
    typedef typename QR_no_alloc<T>::valueType valueType;

  protected:
    integer   * m_JPVT;
    valueType * m_WorkPermute;

  public:

    using QR_no_alloc<T>::m_nrows;
    using QR_no_alloc<T>::m_ncols;
    using QR_no_alloc<T>::m_nReflector;
    using QR_no_alloc<T>::m_LworkQR;

    using QR_no_alloc<T>::m_Afactorized;
    using QR_no_alloc<T>::m_WorkQR;
    using QR_no_alloc<T>::m_Tau;

    using QR_no_alloc<T>::solve;
    using QR_no_alloc<T>::t_solve;
    using QR_no_alloc<T>::Q_mul;
    using QR_no_alloc<T>::Qt_mul;
    using QR_no_alloc<T>::mul_Q;
    using QR_no_alloc<T>::mul_Qt;
    using QR_no_alloc<T>::invR_mul;
    using QR_no_alloc<T>::invRt_mul;
    using QR_no_alloc<T>::getR;
    using QR_no_alloc<T>::getQ;
    using QR_no_alloc<T>::getA;
    using QR_no_alloc<T>::getTau;
    using QR_no_alloc<T>::no_allocate;
    using QR_no_alloc<T>::factorize;

    QRP_no_alloc()
    : QR_no_alloc<T>()
    , m_JPVT(nullptr)
    {}

    virtual
    ~QRP_no_alloc() UTILS_OVERRIDE
    {}

    integer
    get_Lwork_QRP( integer NR, integer NC ) const;

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
     *  Do QR factorization with column pivoting of a rectangular matrix
     *  \param A   pointer to the matrix
     *  \param LDA Leading dimension of the matrix
     */
    void
    factorize_nodim(
      char const      who[],
      valueType const A[],
      integer         LDA
    );

    bool
    factorize_nodim( valueType const A[], integer LDA );

    // -------------------------------------------------------------------------

    void
    getPerm( integer jpvt[] ) const;

    // -------------------------------------------------------------------------
    void
    permute( valueType x[], integer incx = 1 ) const;

    void
    inv_permute( valueType x[], integer incx = 1 ) const;

    void
    permute_rows(
      integer   nr,
      integer   nc,
      valueType C[],
      integer   ldC
    ) const;

    void
    permute_rows( MatrixWrapper<T> & M ) const
    { permute_rows( M.numRows(), M.numCols(), M.data(), M.lDim() ); }

    void
    inv_permute_rows(
      integer   nr,
      integer   nc,
      valueType C[],
      integer   ldC
    ) const;

    void
    inv_permute_rows( MatrixWrapper<T> & M ) const
    { inv_permute_rows( M.numRows(), M.numCols(), M.data(), M.lDim() ); }

    void
    permute_cols(
      integer   nr,
      integer   nc,
      valueType C[],
      integer   ldC
    ) const;

    void
    permute_cols( MatrixWrapper<T> & M ) const
    { permute_cols( M.numRows(), M.numCols(), M.data(), M.lDim() ); }

    void
    inv_permute_cols(
      integer   nr,
      integer   nc,
      valueType C[],
      integer   ldC
    ) const;

    void
    inv_permute_cols( MatrixWrapper<T> & M ) const
    { inv_permute_cols( M.numRows(), M.numCols(), M.data(), M.lDim() ); }

    integer
    rankEstimate( valueType rcond ) const {
      valueType SVAL[3];
      return lapack_wrapper::rankEstimate(
        m_nrows, m_ncols, m_Afactorized, m_nrows, rcond, SVAL
      );
    }

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    /*!
     *  In case of QR factorization of a square matrix solve the
     *  linear system \f$ QR x = b \f$
     *  \param xb on input the rhs of linear system on output the solution
     */
    virtual
    bool
    solve( valueType xb[] ) const UTILS_OVERRIDE;

    /*!
     *  In case of QR factorization of a square matrix solve the
     *  linear system \f$ (QR)^T x = b \f$
     *  \param xb on input the rhs of linear system on output the solution
     */
    virtual
    bool
    t_solve( valueType xb[] ) const UTILS_OVERRIDE;

    virtual
    bool
    solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const UTILS_OVERRIDE;

    virtual
    bool
    t_solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const UTILS_OVERRIDE;
  };

  //============================================================================
  //============================================================================
  //============================================================================

  template <typename T>
  class QRP : public QRP_no_alloc<T> {
  public:
    typedef typename QRP_no_alloc<T>::valueType valueType;

  private:
    Malloc<valueType> m_allocReals;
    Malloc<integer>   m_allocIntegers;

  public:

    using QRP_no_alloc<T>::m_nrows;
    using QRP_no_alloc<T>::m_ncols;
    using QRP_no_alloc<T>::m_nReflector;
    using QRP_no_alloc<T>::m_LworkQR;

    using QRP_no_alloc<T>::m_Afactorized;
    using QRP_no_alloc<T>::m_WorkQR;
    using QRP_no_alloc<T>::m_Tau;

    using QRP_no_alloc<T>::solve;
    using QRP_no_alloc<T>::t_solve;
    using QRP_no_alloc<T>::Q_mul;
    using QRP_no_alloc<T>::Qt_mul;
    using QRP_no_alloc<T>::mul_Q;
    using QRP_no_alloc<T>::mul_Qt;
    using QRP_no_alloc<T>::invR_mul;
    using QRP_no_alloc<T>::invRt_mul;
    using QRP_no_alloc<T>::getR;
    using QRP_no_alloc<T>::getQ;
    using QRP_no_alloc<T>::getA;
    using QRP_no_alloc<T>::getTau;
    using QRP_no_alloc<T>::no_allocate;
    using QRP_no_alloc<T>::factorize;

    QRP()
    : QRP_no_alloc<T>()
    , m_allocReals("QRP-allocReals")
    , m_allocIntegers("QRP-allocIntegers")
    {}

    virtual
    ~QRP() UTILS_OVERRIDE
    { m_allocIntegers.free(); }

    void
    allocate( integer nr, integer nc );

    /*!
     *  Do QR factorization with column pivoting of a rectangular matrix
     *  \param NR  number of rows of the matrix
     *  \param NC  number of columns of the matrix
     *  \param A   pointer to the matrix
     *  \param LDA Leading dimension of the matrix
     */
    void
    factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) UTILS_OVERRIDE {
      this->allocate( NR, NC );
      this->QRP_no_alloc<T>::factorize_nodim( who, A, LDA );
    }

    /*!
     *  Do QR factorization with column pivoting of a rectangular matrix
     *  \param NR  number of rows of the matrix
     *  \param NC  number of columns of the matrix
     *  \param A   pointer to the matrix
     *  \param LDA Leading dimension of the matrix
     */
    bool
    factorize(
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) UTILS_OVERRIDE {
      this->allocate( NR, NC );
      return this->QRP_no_alloc<T>::factorize_nodim( A, LDA );
    }

  };

}

///
/// eof: qr.hxx
///
