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
    using real_type = typename LinearSystemSolver<T>::real_type;

  protected:

    integer m_nrows{0};
    integer m_ncols{0};
    integer m_nReflector{0};
    integer m_LworkQR{0};

    real_type * m_Afactorized{nullptr};
    real_type * m_WorkQR{nullptr};
    real_type * m_Tau{nullptr};

  public:

    using LinearSystemSolver<T>::factorize;
    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;

    QR_no_alloc() : LinearSystemSolver<T>() {}

    static integer get_Lwork_QR( integer NR, integer NC );

    void
    no_allocate(
      integer     NR,
      integer     NC,
      integer     Lwork,
      real_type * Work
    );

    //!
    //! Do QR factorization of a rectangular matrix.
    //!
    //! \param who string used in error message
    //! \param A   pointer to the matrix
    //! \param LDA Leading dimension of the matrix
    //!
    void
    factorize_nodim(
      string_view     who,
      real_type const A[],
      integer         LDA
    );

    bool
    factorize_nodim( real_type const A[], integer LDA );

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
      real_type     C[],
      integer       ldC
    ) const;

    //! `x <- Q*x`
    void Q_mul( real_type x[] ) const;

    //! `x <- Q'*x`
    void Qt_mul( real_type x[] ) const;

    //! `C <- Q*C`
    void Q_mul( integer nr, integer nc, real_type C[], integer ldC ) const;

    void
    Q_mul( MatrixWrapper<T> & C ) const {
      this->Q_mul( C.nrows(), C.ncols(), C.data(), C.ldim() );
    }

    //! `C <- Q'*C`
    void Qt_mul( integer nr, integer nc, real_type C[], integer ldC ) const;

    void
    Qt_mul( MatrixWrapper<T> & C ) const {
      this->Qt_mul( C.nrows(), C.ncols(), C.data(), C.ldim() );
    }

    //! `C <- C*Q`
    void mul_Q( integer nr, integer nc, real_type C[], integer ldC ) const;

    void
    mul_Q( MatrixWrapper<T> & C ) const {
      this->mul_Q( C.nrows(), C.ncols(), C.data(), C.ldim() );
    }

    //! `C <- C*Q'`
    void mul_Qt( integer nr, integer nc, real_type C[], integer ldC ) const;

    void
    mul_Qt( MatrixWrapper<T> & C ) const {
      this->mul_Q( C.nrows(), C.ncols(), C.data(), C.ldim() );
    }

    // -------------------------------------------------------------------------

    //! `x <- R^(-1) * x`
    void invR_mul( real_type x[], integer incx = 1 ) const;

    //! `C <- R^(-1) * C`
    void invR_mul( integer nr, integer nc, real_type C[], integer ldC ) const;

    void
    invR_mul( MatrixWrapper<real_type> & C ) const {
      this->invR_mul( C.nrows(), C.ncols(), C.data(), C.ldim() );
    }

    // -------------------------------------------------------------------------

    //! `x <- R^(-T) * x`
    void invRt_mul( real_type x[], integer incx = 1 ) const;

    //! `C <- R^(-T) * C`
    void invRt_mul( integer nr, integer nc, real_type C[], integer ldC ) const;

    void
    invRt_mul( MatrixWrapper<real_type> & C ) const {
      this->invRt_mul( C.nrows(), C.ncols(), C.data(), C.ldim() );
    }

    // -------------------------------------------------------------------------

    //! `C <- C * R^(-1)`
    void mul_invR( integer nr, integer nc, real_type C[], integer ldC ) const;

    void
    mul_invR( MatrixWrapper<real_type> & C ) const {
      this->mul_invR( C.nrows(), C.ncols(), C.data(), C.ldim() );
    }

    // -------------------------------------------------------------------------

    void mul_invRt( integer nr, integer nc, real_type C[], integer ldC ) const;

    void
    mul_invRt( MatrixWrapper<real_type> & C ) const {
      this->mul_invRt( C.nrows(), C.ncols(), C.data(), C.ldim() );
    }

    // -------------------------------------------------------------------------

    void getR( real_type R[], integer ldR ) const;
    void getR( Matrix<T> & R ) const;
    void getRt( real_type R[], integer ldR ) const;
    void getRt( Matrix<T> & R ) const;

    // -------------------------------------------------------------------------

    void getQ( real_type Q[], integer ldQ ) const;
    void getQ( Matrix<T> & Q ) const;

    // -------------------------------------------------------------------------

    void getQreduced( real_type Q[], integer ldQ ) const;
    void getQreduced( Matrix<T> & Q ) const;

    // -------------------------------------------------------------------------

    void getA( real_type A[], integer ldA ) const;
    void getA( Matrix<T> & A ) const;

    // -------------------------------------------------------------------------

    void getTau( real_type tau[] ) const;
    void getTau( Matrix<T> & tau ) const;

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    //!
    //! In case of QR factorization of a square matrix solve the
    //! linear system \f$ QR x = b \f$.
    //!
    //! \param xb on input the rhs of linear system on output the solution
    //!
    bool
    solve( real_type xb[] ) const override;

    //!
    //! In case of QR factorization of a square matrix solve the
    //! linear system \f$ (QR)^T x = b \f$.
    //!
    //! \param xb on input the rhs of linear system on output the solution
    //!
    bool
    t_solve( real_type xb[] ) const override;

    bool
    solve(
      integer   nrhs,
      real_type B[],
      integer   ldB
    ) const override;

    bool
    t_solve(
      integer   nrhs,
      real_type B[],
      integer   ldB
    ) const override;

  };

  //============================================================================
  //============================================================================
  //============================================================================

  template <typename T>
  class QR : public QR_no_alloc<T> {
  public:
    using real_type = typename QR_no_alloc<T>::real_type;

  protected:

    Malloc<real_type> m_allocReals{"QR-allocReals"};

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

    QR() : QR_no_alloc<T>() {}

    ~QR() override { m_allocReals.free(); }

    void
    allocate( integer nr, integer nc );

    //!
    //! Do QR factorization of a rectangular matrix.
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
    using real_type = typename QR_no_alloc<T>::real_type;

  protected:
    integer   * m_JPVT{nullptr};
    real_type * m_WorkPermute{nullptr};

  private:
    using QR_no_alloc<T>::get_Lwork_QR;

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
    using QR_no_alloc<T>::factorize;

    QRP_no_alloc() : QR_no_alloc<T>() {}

    static integer get_Lwork_QRP( integer NR, integer NC );

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
    //! Do QR factorization with column pivoting of a rectangular matrix.
    //!
    //! \param who string used in error message
    //! \param A   pointer to the matrix
    //! \param LDA Leading dimension of the matrix
    //!
    void
    factorize_nodim(
      string_view     who,
      real_type const A[],
      integer         LDA
    );

    bool
    factorize_nodim( real_type const A[], integer LDA );

    // -------------------------------------------------------------------------

    void
    getPerm( integer jpvt[] ) const;

    // -------------------------------------------------------------------------
    void
    permute( real_type x[], integer incx = 1 ) const;

    void
    inv_permute( real_type x[], integer incx = 1 ) const;

    void
    permute_rows(
      integer   nr,
      integer   nc,
      real_type C[],
      integer   ldC
    ) const;

    void
    permute_rows( MatrixWrapper<T> & M ) const
    { permute_rows( M.nrows(), M.ncols(), M.data(), M.ldim() ); }

    void
    inv_permute_rows(
      integer   nr,
      integer   nc,
      real_type C[],
      integer   ldC
    ) const;

    void
    inv_permute_rows( MatrixWrapper<T> & M ) const
    { inv_permute_rows( M.nrows(), M.ncols(), M.data(), M.ldim() ); }

    void
    permute_cols(
      integer   nr,
      integer   nc,
      real_type C[],
      integer   ldC
    ) const;

    void
    permute_cols( MatrixWrapper<T> & M ) const
    { permute_cols( M.nrows(), M.ncols(), M.data(), M.ldim() ); }

    void
    inv_permute_cols(
      integer   nr,
      integer   nc,
      real_type C[],
      integer   ldC
    ) const;

    void
    inv_permute_cols( MatrixWrapper<T> & M ) const
    { inv_permute_cols( M.nrows(), M.ncols(), M.data(), M.ldim() ); }

    integer
    rankEstimate( real_type rcond ) const {
      real_type SVAL[3];
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

    //!
    //! In case of QR factorization of a square matrix solve the
    //! linear system \f$ QR x = b \f$.
    //!
    //! \param xb on input the rhs of linear system on output the solution
    //!
    bool
    solve( real_type xb[] ) const override;

    //!
    //! In case of QR factorization of a square matrix solve the
    //! linear system \f$ (QR)^T x = b \f$.
    //!
    //! \param xb on input the rhs of linear system on output the solution
    //!
    bool
    t_solve( real_type xb[] ) const override;

    bool
    solve(
      integer   nrhs,
      real_type B[],
      integer   ldB
    ) const override;

    bool
    t_solve(
      integer   nrhs,
      real_type B[],
      integer   ldB
    ) const override;
  };

  //============================================================================
  //============================================================================
  //============================================================================

  template <typename T>
  class QRP : public QRP_no_alloc<T> {
  public:
    using real_type = typename QRP_no_alloc<T>::real_type;

  private:
    Malloc<real_type> m_allocReals{"QRP-allocReals"};
    Malloc<integer>   m_allocIntegers{"QRP-allocIntegers"};

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

    QRP() : QRP_no_alloc<T>() {}

    ~QRP() override { m_allocIntegers.free(); }

    void
    allocate( integer nr, integer nc );

    //!
    //! Do QR factorization with column pivoting of a rectangular matrix.
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
      this->QRP_no_alloc<T>::factorize_nodim( who, A, LDA );
    }

    //!
    //! Do QR factorization with column pivoting of a rectangular matrix.
    //!
    //! \param NR  number of rows of the matrix
    //! \param NC  number of columns of the matrix
    //! \param A   pointer to the matrix
    //! \param LDA Leading dimension of the matrix
    //!
    bool
    factorize(
      integer         NR,
      integer         NC,
      real_type const A[],
      integer         LDA
    ) override {
      this->allocate( NR, NC );
      return this->QRP_no_alloc<T>::factorize_nodim( A, LDA );
    }

  };

}

///
/// eof: qr.hxx
///
