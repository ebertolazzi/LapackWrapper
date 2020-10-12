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
  :|:   _ __ (_)_ ____   __
  :|:  | '_ \| | '_ \ \ / /
  :|:  | |_) | | | | \ V /
  :|:  | .__/|_|_| |_|\_/
  :|:  |_|
  \*/

  template <typename T>
  class PINV_no_alloc : public LinearSystemSolver<T> {
  public:
    typedef typename LinearSystemSolver<T>::valueType valueType;

  protected:

    integer     m_nRows;
    integer     m_nCols;
    integer     m_rank;
    integer     m_minRC;
    valueType   m_epsi;
    valueType * m_Rscale;
    valueType * m_Cscale;
    valueType * m_Ascaled;
    valueType * m_Rt;

    integer     m_LWorkQR2;
    valueType * m_WorkQR2;
    integer   * m_iWorkQR2;

    QRP_no_alloc<valueType> m_QR1;
    //QRP_no_alloc<valueType> m_QR2;
    QR_no_alloc<valueType> m_QR2;

    EquilibrationType m_equ;

  public:

    using LinearSystemSolver<T>::factorize;
    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;

    PINV_no_alloc();

    virtual
    ~PINV_no_alloc() LAPACK_WRAPPER_OVERRIDE {}

    void
    no_allocate(
      integer     NR,
      integer     NC,
      integer     Lwork,
      valueType * Work,
      integer     Liwork,
      integer   * iWork
    );

    void
    get_required_allocation(
      integer   NR,
      integer   NC,
      integer & Lwork,
      integer & Liwork
    ) const;

    void
    factorize( char const who[], valueType const A[], integer LDA );

    bool
    factorize( valueType const A[], integer LDA );

    bool
    mult_inv( valueType const b[], valueType x[] ) const;

    bool
    t_mult_inv( valueType const b[], valueType x[] ) const;

    bool
    mult_inv(
      integer         nrhs,
      valueType const B[],
      integer         ldB,
      valueType       X[],
      integer         ldX
    ) const;

    bool
    t_mult_inv(
      integer         nrhs,
      valueType const B[],
      integer         ldB,
      valueType       X[],
      integer         ldX
    ) const;

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
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE {
      if ( m_nRows != m_nCols ) return false;
      return mult_inv( xb, xb );
    }

    /*!
     *  In case of QR factorization of a square matrix solve the
     *  linear system \f$ (QR)^T x = b \f$
     *  \param xb on input the rhs of linear system on output the solution
     */
    virtual
    bool
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE {
      if ( m_nRows != m_nCols ) return false;
      return t_mult_inv( xb, xb );
    }

    virtual
    bool
    solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const LAPACK_WRAPPER_OVERRIDE {
      if ( m_nRows != m_nCols ) return false;
      return mult_inv( nrhs, B, ldB, B, ldB );
    }

    virtual
    bool
    t_solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const LAPACK_WRAPPER_OVERRIDE {
      if ( m_nRows != m_nCols ) return false;
      return t_mult_inv( nrhs, B, ldB, B, ldB );
    }

  };

  //============================================================================
  //============================================================================
  //============================================================================

  template <typename T>
  class PINV : public PINV_no_alloc<T> {
  public:
    typedef typename PINV_no_alloc<T>::valueType valueType;

  protected:

    Malloc<valueType> m_allocReals;
    Malloc<integer>   m_allocIntegers;

  public:

    using PINV_no_alloc<T>::mult_inv;
    using PINV_no_alloc<T>::t_mult_inv;
    using PINV_no_alloc<T>::solve;
    using PINV_no_alloc<T>::t_solve;
    using PINV_no_alloc<T>::get_required_allocation;
    using PINV_no_alloc<T>::no_allocate;
    using PINV_no_alloc<T>::factorize;

    PINV()
    : PINV_no_alloc<T>()
    , m_allocReals("PINV-allocReals")
    , m_allocIntegers("PINV-allocIntegers")
    {}

    virtual
    ~PINV() LAPACK_WRAPPER_OVERRIDE
    { m_allocReals.free(); }

    void
    allocate( integer nr, integer nc );

    /*!
     *  \param NR  number of rows of the matrix
     *  \param NC  number of columns of the matrix
     *  \param A   pointer to the matrix
     *  \param LDA Leading dimension of the matrix
     */
    void
    factorize(
      char      const who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) LAPACK_WRAPPER_OVERRIDE {
      this->allocate( NR, NC );
      this->PINV_no_alloc<T>::factorize( who, A, LDA );
    }

    /*!
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
    ) LAPACK_WRAPPER_OVERRIDE {
      this->allocate( NR, NC );
      return this->PINV_no_alloc<T>::factorize( A, LDA );
    }

  };

}

///
/// eof: pinv.hxx
///
