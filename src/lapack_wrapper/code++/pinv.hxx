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
    using real_type = typename LinearSystemSolver<T>::real_type;

  protected:

    mutable Malloc<real_type> m_alloc_work;
    mutable integer           L_mm_work;
    mutable real_type       * mm_work;

    integer     m_nrows;
    integer     m_ncols;
    integer     m_rank;
    integer     m_rcmax;
    real_type   m_epsi;
    real_type * m_A_factored;
    real_type * m_Rt;

    integer     m_LWorkQR2;
    real_type * m_WorkQR2;

    QRP_no_alloc<real_type> m_QR1;
    QR_no_alloc<real_type>  m_QR2;

    void factorize( char const who[] );
    bool factorize();

  public:

    using LinearSystemSolver<T>::factorize;
    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;

    PINV_no_alloc();
    ~PINV_no_alloc() override {}

    void
    no_allocate(
      integer     NR,
      integer     NC,
      integer     Lwork,
      real_type * Work,
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

    void factorize_nodim( char const who[], real_type const A[], integer LDA );
    bool factorize_nodim( real_type const A[], integer LDA );
    void t_factorize_nodim( char const who[], real_type const A[], integer LDA );
    bool t_factorize_nodim( real_type const A[], integer LDA );

    bool
    mult_inv(
      real_type const b[],
      integer         incb,
      real_type       x[],
      integer         incx
    ) const;

    bool
    t_mult_inv(
      real_type const b[],
      integer         incb,
      real_type       x[],
      integer         incx
    ) const;

    bool
    mult_inv(
      integer         nrhs,
      real_type const B[],
      integer         ldB,
      real_type       X[],
      integer         ldX
    ) const;

    bool
    t_mult_inv(
      integer         nrhs,
      real_type const B[],
      integer         ldB,
      real_type       X[],
      integer         ldX
    ) const;

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
    solve( real_type xb[] ) const override {
      if ( m_nrows != m_ncols ) return false;
      return mult_inv( xb, 1, xb, 1 );
    }

    //!
    //! In case of QR factorization of a square matrix solve the
    //! linear system \f$ (QR)^T x = b \f$.
    //!
    //! \param xb on input the rhs of linear system on output the solution
    //!
    bool
    t_solve( real_type xb[] ) const override {
      if ( m_nrows != m_ncols ) return false;
      return t_mult_inv( xb, 1, xb, 1 );
    }

    bool
    solve(
      integer   nrhs,
      real_type B[],
      integer   ldB
    ) const override {
      if ( m_nrows != m_ncols ) return false;
      return mult_inv( nrhs, B, ldB, B, ldB );
    }

    bool
    t_solve(
      integer   nrhs,
      real_type B[],
      integer   ldB
    ) const override {
      if ( m_nrows != m_ncols ) return false;
      return t_mult_inv( nrhs, B, ldB, B, ldB );
    }

    void
    info( ostream_type & stream ) const;

  };

  //============================================================================
  //============================================================================
  //============================================================================

  template <typename T>
  class PINV : public PINV_no_alloc<T> {
  public:
    using real_type = typename PINV_no_alloc<T>::real_type;

  protected:

    Malloc<real_type> m_allocReals;
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

    ~PINV() override
    { m_allocReals.free(); }

    void
    allocate( integer nr, integer nc );

    //!
    //! \param who string used in error message
    //! \param NR  number of rows of the matrix
    //! \param NC  number of columns of the matrix
    //! \param A   pointer to the matrix
    //! \param LDA Leading dimension of the matrix
    //!
    void
    factorize(
      char      const who[],
      integer         NR,
      integer         NC,
      real_type const A[],
      integer         LDA
    ) override {
      this->allocate( NR, NC );
      this->PINV_no_alloc<T>::factorize_nodim( who, A, LDA );
    }

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
      return this->PINV_no_alloc<T>::factorize_nodim( A, LDA );
    }

    //!
    //! \param who string used in error message
    //! \param NR  number of rows of the matrix
    //! \param NC  number of columns of the matrix
    //! \param A   pointer to the matrix
    //! \param LDA Leading dimension of the matrix
    //!
    void
    t_factorize(
      char      const who[],
      integer         NR,
      integer         NC,
      real_type const A[],
      integer         LDA
    ) override {
      this->allocate( NC, NR );
      this->PINV_no_alloc<T>::t_factorize_nodim( who, A, LDA );
    }

    //!
    //! \param NR  number of rows of the matrix
    //! \param NC  number of columns of the matrix
    //! \param A   pointer to the matrix
    //! \param LDA Leading dimension of the matrix
    //!
    bool
    t_factorize(
      integer         NR,
      integer         NC,
      real_type const A[],
      integer         LDA
    ) override {
      this->allocate( NC, NR );
      return this->PINV_no_alloc<T>::t_factorize_nodim( A, LDA );
    }

  };

}

///
/// eof: pinv.hxx
///
