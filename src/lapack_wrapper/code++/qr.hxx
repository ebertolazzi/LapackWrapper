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

    integer nRows;
    integer nCols;
    integer nReflector;
    integer Lwork;

    valueType * Afactorized;
    valueType * WorkFactorized;
    valueType * Tau;

  public:

    using LinearSystemSolver<T>::factorize;
    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;

    QR_no_alloc()
    : LinearSystemSolver<T>()
    , nRows(0)
    , nCols(0)
    , nReflector(0)
    , Lwork(0)
    , Afactorized(nullptr)
    , WorkFactorized(nullptr)
    , Tau(nullptr)
    {}

    virtual
    ~QR_no_alloc() LAPACK_WRAPPER_OVERRIDE {}

    void
    no_allocate(
      integer     NR,
      integer     NC,
      integer     _Lwork,
      valueType * _Afactorized,
      valueType * _WorkFactorized,
      valueType * _Tau
    ) {
      this->nRows          = NR;
      this->nCols          = NC;
      this->nReflector     = std::min(NR,NC);
      this->Lwork          = _Lwork;
      this->Afactorized    = _Afactorized;
      this->WorkFactorized = _WorkFactorized;
      this->Tau            = _Tau;
    }

    /*!
     * Do QR factorization of a rectangular matrix
     * \param A   pointer to the matrix
     * \param LDA Leading dimension of the matrix
     */
    void
    factorize(
      char const      who[],
      valueType const A[],
      integer         LDA
    );

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
    Q_mul( valueType x[] ) const {
      applyQ(
        LEFT, NO_TRANSPOSE, this->nReflector, this->nRows, 1, x, this->nRows
      );
    }

    //! x <- Q'*x
    void
    Qt_mul( valueType x[] ) const
    { applyQ( LEFT, TRANSPOSE, nReflector, this->nRows, 1, x, this->nRows ); }

    //! C <- Q*C
    void
    Q_mul( integer nr, integer nc, valueType C[], integer ldC ) const
    { applyQ( LEFT, NO_TRANSPOSE, this->nReflector, nr, nc, C, ldC ); }

    void
    Q_mul( MatrixWrapper<T> & C ) const {
      applyQ(
        LEFT, NO_TRANSPOSE, this->nReflector,
        C.numRows(), C.numCols(), C.get_data(), C.lDim()
      );
    }

    //! C <- Q'*C
    void
    Qt_mul( integer nr, integer nc, valueType C[], integer ldC ) const
    { applyQ( LEFT, TRANSPOSE, this->nReflector, nr, nc, C, ldC ); }

    void
    Qt_mul( MatrixWrapper<T> & C ) const {
      applyQ(
        LEFT, TRANSPOSE, this->nReflector,
        C.numRows(), C.numCols(), C.get_data(), C.lDim()
      );
    }

    //! C <- C*Q
    void
    mul_Q( integer nr, integer nc, valueType C[], integer ldC ) const
    { applyQ( RIGHT, NO_TRANSPOSE, this->nReflector, nr, nc, C, ldC ); }

    void
    mul_Q( MatrixWrapper<T> & C ) const {
      applyQ(
        RIGHT, NO_TRANSPOSE, this->nReflector,
        C.numRows(), C.numCols(), C.get_data(), C.lDim()
      );
    }

    //! C <- C*Q'
    void
    mul_Qt( integer nr, integer nc, valueType C[], integer ldC ) const
    { applyQ( RIGHT, TRANSPOSE, this->nReflector, nr, nc, C, ldC ); }

    void
    mul_Qt( MatrixWrapper<T> & C ) const {
      applyQ(
        RIGHT, TRANSPOSE, this->nReflector,
        C.numRows(), C.numCols(), C.get_data(), C.lDim()
      );
    }

    // -------------------------------------------------------------------------

    void
    Rsolve(
      Transposition TRANS,
      integer       rk,
      valueType     x[],
      integer       incx
    ) const {
      trsv( UPPER, TRANS, NON_UNIT, rk, this->Afactorized, this->nRows, x, incx );
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
      trsm(
        SIDE, UPPER, TRANS, NON_UNIT,
        nr, nc, alpha, this->Afactorized, this->nRows, Bmat, ldB
      );
    }

    void
    Rsolve(
      SideMultiply               SIDE,
      Transposition              TRANS,
      valueType                  alpha,
      MatrixWrapper<valueType> & Bmat
    ) const {
      trsm(
        SIDE, UPPER, TRANS, NON_UNIT,
        Bmat.numRows(), Bmat.numCols(),
        alpha, this->Afactorized, this->nRows,
        Bmat.get_data(), Bmat.lDim()
      );
    }

    // -------------------------------------------------------------------------

    //! x <- R^(-1) * x
    void
    invR_mul( valueType x[], integer incx = 1 ) const {
      trsv(
        UPPER, NO_TRANSPOSE, NON_UNIT,
        nReflector, this->Afactorized, this->nRows, x, incx
      );
    }

    //! C <- R^(-1) * C
    void
    invR_mul( integer nr, integer nc, valueType C[], integer ldC ) const
    { Rsolve( LEFT, NO_TRANSPOSE, nr, nc, 1.0, C, ldC ); }

    void
    invR_mul( MatrixWrapper<valueType> & C ) const {
      Rsolve(
        LEFT, NO_TRANSPOSE,
        C.numRows(), C.numCols(), 1.0, C.get_data(), C.lDim()
      );
    }

    // -------------------------------------------------------------------------

    //! x <- R^(-T) * x
    void
    invRt_mul( valueType x[], integer incx = 1 ) const {
      trsv(
        UPPER, TRANSPOSE, NON_UNIT,
        nReflector, this->Afactorized, this->nRows, x, incx
      );
    }

    //! C <- R^(-T) * C
    void
    invRt_mul( integer nr, integer nc, valueType C[], integer ldC ) const
    { Rsolve( LEFT, TRANSPOSE, nr, nc, 1.0, C, ldC ); }

    void
    invRt_mul( MatrixWrapper<valueType> & C ) const {
      Rsolve(
        LEFT, TRANSPOSE, C.numRows(), C.numCols(), 1.0, C.get_data(), C.lDim()
      );
    }

    // -------------------------------------------------------------------------

    //! C <- C * R^(-1)
    void
    mul_invR( integer nr, integer nc, valueType C[], integer ldC ) const
    { Rsolve( RIGHT, NO_TRANSPOSE, nr, nc, 1.0, C, ldC ); }

    void
    mul_invR( MatrixWrapper<valueType> & C ) const {
      Rsolve(
        RIGHT, NO_TRANSPOSE,
        C.numRows(), C.numCols(), 1.0, C.get_data(), C.lDim()
      );
    }

    // -------------------------------------------------------------------------

    void
    mul_invRt( integer nr, integer nc, valueType C[], integer ldC ) const
    { Rsolve( RIGHT, TRANSPOSE, nr, nc, 1.0, C, ldC ); }

    void
    mul_invRt( MatrixWrapper<valueType> & C ) const
    { Rsolve( RIGHT, TRANSPOSE, C.numRows(), C.numCols(), 1.0, C.get_data(), C.lDim( ) ); }

    // -------------------------------------------------------------------------

    void
    getR( valueType R[], integer ldR ) const;

    void
    getR( Matrix<T> & R ) const;

    // -------------------------------------------------------------------------

    void
    getQ( valueType Q[], integer ldQ ) const;

    void
    getQ( Matrix<T> & Q ) const;

    // -------------------------------------------------------------------------

    void
    getA( valueType A[], integer ldA ) const;

    void
    getA( Matrix<T> & A ) const;

    // -------------------------------------------------------------------------

    void
    getTau( valueType tau[] ) const;

    void
    getTau( Matrix<T> & tau ) const;

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
    void
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    /*!
     *  In case of QR factorization of a square matrix solve the
     *  linear system \f$ (QR)^T x = b \f$
     *  \param xb on input the rhs of linear system on output the solution
     */
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
  //============================================================================
  //============================================================================

  template <typename T>
  class QR : public QR_no_alloc<T> {
  public:
    typedef typename QR_no_alloc<T>::valueType valueType;

  protected:

    Malloc<valueType> allocReals;

  public:

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

    QR()
    : QR_no_alloc<T>()
    , allocReals("QR-allocReals")
    {}

    virtual
    ~QR() LAPACK_WRAPPER_OVERRIDE
    { allocReals.free(); }

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
    ) LAPACK_WRAPPER_OVERRIDE {
      this->allocate( NR, NC );
      this->factorize( who, A, LDA );
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
    integer   * JPVT;
    valueType * WorkPermute;

  public:

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
    , JPVT(nullptr)
    {}

    virtual
    ~QRP_no_alloc() LAPACK_WRAPPER_OVERRIDE
    {}

    void
    no_allocate(
      integer     NR,
      integer     NC,
      integer     _Lwork,
      valueType * _Afactorized,
      valueType * _WorkFactorized,
      valueType * _Tau,
      valueType * _WorkPermute,
      integer   * _JPVT
    ) {
      this->no_allocate( NR, NC, _Lwork, _Afactorized, _WorkFactorized, _Tau );
      this->WorkPermute = _WorkPermute;
      this->JPVT        = _JPVT;
    }

    /*!
     *  Do QR factorization with column pivoting of a rectangular matrix
     *  \param A   pointer to the matrix
     *  \param LDA Leading dimension of the matrix
     */
    void
    factorize(
      char const      who[],
      valueType const A[],
      integer         LDA
    );

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
    { permute_rows( M.numRows(), M.numCols(), M.get_data(), M.lDim() ); }

    void
    inv_permute_rows(
      integer   nr,
      integer   nc,
      valueType C[],
      integer   ldC
    ) const;

    void
    inv_permute_rows( MatrixWrapper<T> & M ) const
    { inv_permute_rows( M.numRows(), M.numCols(), M.get_data(), M.lDim() ); }

    void
    permute_cols(
      integer   nr,
      integer   nc,
      valueType C[],
      integer   ldC
    ) const;

    void
    permute_cols( MatrixWrapper<T> & M ) const
    { permute_cols( M.numRows(), M.numCols(), M.get_data(), M.lDim() ); }

    void
    inv_permute_cols(
      integer   nr,
      integer   nc,
      valueType C[],
      integer   ldC
    ) const;

    void
    inv_permute_cols( MatrixWrapper<T> & M ) const
    { inv_permute_cols( M.numRows(), M.numCols(), M.get_data(), M.lDim() ); }

    integer
    rankEstimate( valueType rcond ) const {
      valueType SVAL[3];
      return lapack_wrapper::rankEstimate(
        this->nRows, this->nCols, this->Afactorized, this->nRows, rcond, SVAL
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
    void
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    /*!
     *  In case of QR factorization of a square matrix solve the
     *  linear system \f$ (QR)^T x = b \f$
     *  \param xb on input the rhs of linear system on output the solution
     */
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
  //============================================================================
  //============================================================================

  template <typename T>
  class QRP : public QRP_no_alloc<T> {
  public:
    typedef typename QRP_no_alloc<T>::valueType valueType;

  private:
    Malloc<valueType> allocReals;
    Malloc<integer>   allocIntegers;

  public:

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
    , allocReals("QRP-allocReals")
    , allocIntegers("QRP-allocIntegers")
    {}

    virtual
    ~QRP() LAPACK_WRAPPER_OVERRIDE
    { allocIntegers.free(); }

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
    ) LAPACK_WRAPPER_OVERRIDE {
      this->allocate( NR, NC );
      this->QRP_no_alloc<T>::factorize( who, A, LDA );
    }

  };

}

///
/// eof: qr.hxx
///
