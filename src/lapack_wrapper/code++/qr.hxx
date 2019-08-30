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
      std::fill( JPVT, JPVT+nCol, integer(0) );
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

}

///
/// eof: qr.hxx
///
