/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
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
/// file: qr.cxx
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
  PINV_no_alloc<T>::PINV_no_alloc()
  : LinearSystemSolver<T>()
  , m_allocReals("PINV_no_alloc")
  , L_mm_work(0)
  , mm_work(nullptr)
  , m_nRows(0)
  , m_nCols(0)
  , m_rank(0)
  , m_Rscale(nullptr)
  , m_Cscale(nullptr)
  , m_Ascaled(nullptr)
  , m_Rt(nullptr)
  {
    m_epsi = std::pow( std::numeric_limits<T>::epsilon(), T(0.65) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  PINV_no_alloc<T>::get_required_allocation(
    integer   NR,
    integer   NC,
    integer & Lwork,
    integer & Liwork
  ) const {
    integer minRC = std::min( NR, NC );
    integer L1    = m_QR1.get_Lwork( NR, NC ) + (NR+1)*(NC+1);
    integer L2    = m_QR2.get_Lwork( minRC, minRC ) + (minRC+1)*(minRC+1) + 2*minRC;
    Lwork  = L1+L2 + (NR+1)*(NC+1) + minRC*(minRC+2);
    Liwork = NC+minRC;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  PINV_no_alloc<T>::no_allocate(
    integer     NR,
    integer     NC,
    integer     Lwork,
    valueType * Work,
    integer     Liwork,
    integer   * iWork
  ) {

    LW_ASSERT(
      NR >= NC,
      "PINV_no_alloc::no_allocate( NR = {}, NC = {},...)\n"
      "must be NR >= NC\n",
      NR, NC
    );

    integer L1   = m_QR1.get_Lwork( NR, NC ) + (NR+1)*(NC+1);
    integer L2   = m_QR2.get_Lwork( NC, NC ) + (NC+1)*(NC+1);
    integer Ltot = L1 + L2 + (NR+1)*(NC+1) + NC*NC;
    LW_ASSERT(
      Lwork >= Ltot && Liwork >= NC,
      "PINV_no_alloc::no_allocate( NR = {}, NC = {}, Lwork = {},..., Liwork = {},...)\n"
      "Lwork must be >= {} and Liwork >= {}\n",
      NR, NC, Lwork, Liwork, Ltot, NC
    );
    valueType * ptr = Work;
    m_QR1.no_allocate( NR, NC, L1, ptr, NC, iWork ); ptr += L1;
    m_Rscale  = ptr; ptr += NR;
    m_Cscale  = ptr; ptr += NC;
    m_Ascaled = ptr; ptr += NR*NC;
    m_Rt      = ptr; ptr += NC*NC;

    m_LWorkQR2 = L2;
    m_WorkQR2  = ptr;

    m_nRows    = NR;
    m_nCols    = NC;
    m_rank     = 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  PINV_no_alloc<T>::factorize(
    char      const who[],
    valueType const A[],
    integer         LDA
  ) {
    // copy matrix and scale
    integer info = gecopy( m_nRows, m_nCols, A, LDA, m_Ascaled, m_nRows );
    LW_ASSERT(
      info == 0,
      "{}, call to gecopy( N = {}, M = {}, A, LDA = {}, B, LDB = {}\n"
      "return info = {}\n",
      who, m_nRows, m_nCols, LDA, m_nRows, info
    );

    // balance matrix
    valueType ROWCND, COLCND, AMAX;
    ROWCND = COLCND = 1; AMAX = 0;
    fill( m_nRows, m_Rscale, 1, 1 );
    fill( m_nCols, m_Cscale, 1, 1 );
    info = geequ(
      m_nRows, m_nCols, m_Ascaled, m_nRows, m_Rscale, m_Cscale,
      ROWCND, COLCND, AMAX
    );

    LW_ASSERT(
      info >= 0,
      "{}, call to geequ( N = {}, M = {}, A, LDA = {}, ...)\n"
      "return info = {}\n",
      who, m_nRows, m_nCols, LDA, info
    );

    // equilibrate
    laqge(
      m_nRows, m_nCols, m_Ascaled, m_nRows, m_Rscale, m_Cscale,
      ROWCND, COLCND, AMAX, m_equ
    );

    // perform QR factorization
    m_QR1.factorize( who, m_Ascaled, m_nRows );

    // evaluate rank
    m_QR1.getRt( m_Rt, m_nCols );

    m_rank = m_nCols;
    valueType threshold = absmax( m_nCols, m_Rt, m_nCols+1 ) * m_epsi;
    for ( integer i = 1; i < m_rank; ++i ) {
      if ( std::abs(m_Rt[i*(m_nCols+1)]) < threshold )
        { m_rank = i; break; }
    }

    if ( m_rank < m_nCols ) {
      // allocate for QR factorization
      m_QR2.no_allocate( m_nCols, m_rank, m_LWorkQR2, m_WorkQR2 );
      // perform QR factorization
      m_QR2.factorize( who, m_Rt, m_nCols );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  PINV_no_alloc<T>::factorize(
    valueType const A[],
    integer         LDA
  ) {
    // copy matrix and scale
    integer info = gecopy( m_nRows, m_nCols, A, LDA, m_Ascaled, m_nRows );
    if ( info != 0 ) return false;

    // balance matrix
    valueType ROWCND, COLCND, AMAX;
    ROWCND = COLCND = 1; AMAX = 0;
    fill( m_nRows, m_Rscale, 1, 1 );
    fill( m_nCols, m_Cscale, 1, 1 );
    info = geequ(
      m_nRows, m_nCols, m_Ascaled, m_nRows, m_Rscale, m_Cscale,
      ROWCND, COLCND, AMAX
    );

    if ( info < 0 ) return false;

    // equilibrate
    laqge(
      m_nRows, m_nCols, m_Ascaled, m_nRows, m_Rscale, m_Cscale,
      ROWCND, COLCND, AMAX, m_equ
    );

    // perform QR factorization
    bool ok = m_QR1.factorize( m_Ascaled, m_nRows );
    if ( !ok ) return false;

    // evaluate rank
    m_QR1.getRt( m_Rt, m_nCols );

    m_rank = m_nCols;
    valueType threshold = absmax( m_nCols, m_Rt, m_nCols+1 ) * m_epsi;
    for ( integer i = 1; i < m_rank; ++i ) {
      if ( std::abs(m_Rt[i*(m_nCols+1)]) < threshold )
        { m_rank = i; break; }
    }

    if ( m_rank < m_nCols ) {
      // allocate for QR factorization
      m_QR2.no_allocate( m_nCols, m_rank, m_LWorkQR2, m_WorkQR2 );
      // perform QR factorization
      ok = m_QR2.factorize( m_Rt, m_nCols );
      if ( !ok ) return false;
    }

    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  PINV_no_alloc<T>::mult_inv(
    valueType const b[],
    integer         incb,
    valueType       x[],
    integer         incx
  ) const {

    if ( L_mm_work < m_nRows ) {
      L_mm_work = m_nRows;
      m_allocReals.free();
      m_allocReals.allocate( size_t(L_mm_work) );
      mm_work = m_allocReals( size_t(L_mm_work) );
    }

    // (S*A*C)*P = Q*R --> A = S^(-1) * Q * R * P^(-1) * C^(-1)
    // A^+ = C P R^(-1) Q^T S
    copy( m_nRows, b, incb, mm_work, 1 );
    if ( m_equ == EQUILIBRATE_BOTH || m_equ == EQUILIBRATE_ROWS )
      lascl2( m_nRows, 1, m_Rscale, mm_work, 1 );

    m_QR1.Qt_mul( mm_work );

    if ( m_rank < m_nCols ) {
      /*
      // computation of R^+
      //        / L1 | 0 \
      //  R^T = |    |   | = ( Q1 * R1 | Q2 * 0 )
      //        \ M1 | 0 /
      //
      //  / L1 \+
      //  |    |  =  ( R1^T Q1^T Q1 R1 )^+ ( L1^T M1^T )
      //  \ M1 /
      //          =  ( R1^T R1 )^+ ( L1^T M1^T )
      //          =  R1^(-1) R1^(-T) ( L1^T M1^T )
      //          =  R1^(-1) R1^(-T) R1^T Q1^T
      //          =  R1^(-1) Q1^T
      //
      //            / L1 | 0 \+  / R1^(-1) Q1^T \
      //  (R^T)^+ = |    |   | = |              |
      //            \ M1 | 0 /   \   0 * Q2^2   /
      //                                            / R1^(-T)   0 \
      //  R^+ = ( Q1 * R1^-T,  Q2 * 0 ) = ( Q1 Q2 ) |             |
      //                                            \   0       0 /
      */
      m_QR2.invRt_mul( mm_work );
      std::fill_n( mm_work+m_rank, m_nRows-m_rank, 0 );
      m_QR2.Q_mul( mm_work );
    } else {
      m_QR1.invR_mul( mm_work );
    }

    m_QR1.permute( mm_work );

    if ( m_equ == EQUILIBRATE_BOTH || m_equ == EQUILIBRATE_COLUMNS )
      lascl2( m_nCols, 1, m_Cscale, mm_work, 1 );

    copy( m_nCols, mm_work, 1, x, incx );
    return true;
  }

  template <typename T>
  bool
  PINV_no_alloc<T>::t_mult_inv(
    valueType const b[],
    integer         incb,
    valueType       x[],
    integer         incx
  ) const {

    if ( L_mm_work < m_nRows ) {
      L_mm_work = m_nRows;
      m_allocReals.free();
      m_allocReals.allocate( size_t(L_mm_work) );
      mm_work = m_allocReals( size_t(L_mm_work) );
    }

    copy( m_nCols, b, incb, mm_work, 1 );

    if ( m_equ == EQUILIBRATE_BOTH || m_equ == EQUILIBRATE_COLUMNS )
      lascl2( m_nCols, 1, m_Cscale, mm_work, 1 );

    m_QR1.inv_permute( mm_work );

    if ( m_rank < m_nCols ) {
      m_QR2.Qt_mul( mm_work );
      m_QR2.invR_mul( mm_work );
      std::fill_n( mm_work+m_rank, m_nRows-m_rank, 0 );
    } else {
      m_QR1.invRt_mul( mm_work );
    }

    m_QR1.Q_mul( mm_work );

    if ( m_equ == EQUILIBRATE_BOTH || m_equ == EQUILIBRATE_ROWS )
      lascl2( m_nRows, 1, m_Rscale, mm_work, 1 );

    copy( m_nRows, mm_work, 1, x, incx );

    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  PINV_no_alloc<T>::mult_inv(
    integer         nrhs,
    valueType const B[],
    integer         ldB,
    valueType       X[],
    integer         ldX
  ) const {

    integer NR = m_nRows*nrhs;
    if ( L_mm_work < NR ) {
      L_mm_work = NR;
      m_allocReals.free();
      m_allocReals.allocate( size_t(L_mm_work) );
      mm_work = m_allocReals( size_t(L_mm_work) );
    }

    // (S*A*C)*P = Q*R --> A = S^(-1) * Q * R * P^(-1) * C^(-1)
    // A^+ = C P R^(-1) Q^T S
    gecopy( m_nRows, nrhs, B, ldB, mm_work, m_nRows );
    if ( m_equ == EQUILIBRATE_BOTH || m_equ == EQUILIBRATE_ROWS )
      lascl2( m_nRows, nrhs, m_Rscale, mm_work, m_nRows );

    m_QR1.Qt_mul( m_nRows, nrhs, mm_work, m_nRows );

    if ( m_rank < m_nCols ) {
      m_QR2.invRt_mul( m_rank, nrhs, mm_work, m_nRows );
      gezero( m_nRows-m_rank, nrhs, mm_work+m_rank, m_nRows );
      m_QR2.Q_mul( m_nCols, nrhs, mm_work, m_nRows );
    } else {
      m_QR1.invR_mul( m_nRows, nrhs, mm_work, m_nRows );
    }

    m_QR1.permute_rows( m_nCols, nrhs, mm_work, m_nRows );

    if ( m_equ == EQUILIBRATE_BOTH || m_equ == EQUILIBRATE_COLUMNS )
      lascl2( m_nCols, nrhs, m_Cscale, mm_work, m_nRows );

    gecopy( m_nCols, nrhs, mm_work, m_nRows, X, ldX );

    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  PINV_no_alloc<T>::t_mult_inv(
    integer         nrhs,
    valueType const B[],
    integer         ldB,
    valueType       X[],
    integer         ldX
  ) const {

    integer NR = m_nRows*nrhs;
    if ( L_mm_work < NR ) {
      L_mm_work = NR;
      m_allocReals.free();
      m_allocReals.allocate( size_t(L_mm_work) );
      mm_work = m_allocReals( size_t(L_mm_work) );
    }

    gecopy( m_nCols, nrhs, B, ldB, mm_work, m_nRows );

    if ( m_equ == EQUILIBRATE_BOTH || m_equ == EQUILIBRATE_COLUMNS )
      lascl2( m_nCols, nrhs, m_Cscale, mm_work, m_nRows );

    m_QR1.inv_permute_rows( m_nCols, nrhs, mm_work, m_nRows );

    if ( m_rank < m_nCols ) {
      m_QR2.Qt_mul( m_nCols, nrhs, mm_work, m_nRows );
      m_QR2.invR_mul( m_rank, nrhs, mm_work, m_nRows );
      gezero( m_nRows-m_rank, nrhs, mm_work+m_rank, m_nRows );
    } else {
      m_QR1.invRt_mul( m_nRows, nrhs, mm_work, m_nRows );
    }

    m_QR1.Q_mul( m_nRows, nrhs, mm_work, m_nRows );

    if ( m_equ == EQUILIBRATE_BOTH || m_equ == EQUILIBRATE_ROWS )
      lascl2( m_nRows, nrhs, m_Rscale, mm_work, m_nRows );

    gecopy( m_nRows, nrhs, mm_work, m_nRows, X, ldX );

    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  PINV<T>::allocate( integer NR, integer NC ) {
    integer Ltot, Litot;
    this->get_required_allocation( NR, NC, Ltot, Litot );
    m_allocReals.allocate( size_t(Ltot) );
    m_allocIntegers.allocate( size_t(Litot) );
    this->no_allocate(
      NR, NC,
      Ltot,  m_allocReals( size_t(Ltot) ),
      Litot, m_allocIntegers( size_t(Litot) )
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

}

///
/// eof: pinv.cxx
///
