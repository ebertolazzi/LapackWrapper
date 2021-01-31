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
  , m_nrows(0)
  , m_ncols(0)
  , m_rank(0)
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

    UTILS_ASSERT(
      NR >= NC,
      "PINV_no_alloc::no_allocate( NR = {}, NC = {},...)\n"
      "must be NR >= NC\n",
      NR, NC
    );

    integer L1   = m_QR1.get_Lwork( NR, NC ) + (NR+1)*(NC+1);
    integer L2   = m_QR2.get_Lwork( NC, NC ) + (NC+1)*(NC+1);
    integer Ltot = L1 + L2 + 2*NC*NC;
    UTILS_ASSERT(
      Lwork >= Ltot && Liwork >= NC,
      "PINV_no_alloc::no_allocate( NR = {}, NC = {}, Lwork = {},..., Liwork = {},...)\n"
      "Lwork must be >= {} and Liwork >= {}\n",
      NR, NC, Lwork, Liwork, Ltot, NC
    );
    valueType * ptr = Work;
    m_QR1.no_allocate( NR, NC, L1, ptr, NC, iWork ); ptr += L1;
    m_Ascaled = ptr; ptr += NR*NC;
    m_Rt      = ptr; ptr += NC*NC;

    m_LWorkQR2 = L2;
    m_WorkQR2  = ptr;

    m_nrows    = NR;
    m_ncols    = NC;
    m_rank     = 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  PINV_no_alloc<T>::factorize_nodim(
    char      const who[],
    valueType const A[],
    integer         LDA
  ) {
    // copy matrix and scale
    integer info = gecopy( m_nrows, m_ncols, A, LDA, m_Ascaled, m_nrows );
    UTILS_ASSERT(
      info == 0,
      "{}, call to gecopy( N = {}, M = {}, A, LDA = {}, B, LDB = {}\n"
      "return info = {}\n",
      who, m_nrows, m_ncols, LDA, m_nrows, info
    );
    this->factorize( who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  PINV_no_alloc<T>::t_factorize_nodim(
    char      const who[],
    valueType const A[],
    integer         LDA
  ) {
    for ( integer i = 0; i < m_nrows; ++i )
      for ( integer j = 0; j < m_ncols; ++j )
        m_Ascaled[i+j*m_nrows] = A[j+i*LDA];
    this->factorize( who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  PINV_no_alloc<T>::factorize( char const who[] ) {

    // perform QR factorization
    m_QR1.factorize_nodim( who, m_Ascaled, m_nrows );

    // evaluate rank
    m_QR1.getRt( m_Rt, m_ncols );

    m_rank = m_ncols;
    valueType threshold = absmax( m_ncols, m_Rt, m_ncols+1 ) * m_epsi;
    for ( integer i = 1; i < m_rank; ++i ) {
      if ( std::abs(m_Rt[i*(m_ncols+1)]) < threshold )
        { m_rank = i; break; }
    }

    if ( m_rank < m_ncols ) {
      // allocate for QR factorization
      m_QR2.no_allocate( m_ncols, m_rank, m_LWorkQR2, m_WorkQR2 );
      // perform QR factorization
      m_QR2.factorize_nodim( who, m_Rt, m_ncols );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  PINV_no_alloc<T>::factorize_nodim(
    valueType const A[],
    integer         LDA
  ) {
    // copy matrix and scale
    integer info = gecopy( m_nrows, m_ncols, A, LDA, m_Ascaled, m_nrows );
    if ( info != 0 ) return false;
    return this->factorize();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  PINV_no_alloc<T>::t_factorize_nodim(
    valueType const A[],
    integer         LDA
  ) {
    for ( integer i = 0; i < m_nrows; ++i )
      for ( integer j = 0; j < m_ncols; ++j )
        m_Ascaled[i+j*m_nrows] = A[j+i*LDA];
    return this->factorize();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  PINV_no_alloc<T>::factorize() {

    // perform QR factorization
    bool ok = m_QR1.factorize_nodim( m_Ascaled, m_nrows );
    if ( !ok ) return false;

    // evaluate rank
    m_QR1.getRt( m_Rt, m_ncols );

    m_rank = m_ncols;
    valueType threshold = absmax( m_ncols, m_Rt, m_ncols+1 ) * m_epsi;
    for ( integer i = 1; i < m_rank; ++i ) {
      if ( std::abs(m_Rt[i*(m_ncols+1)]) < threshold )
        { m_rank = i; break; }
    }

    if ( m_rank < m_ncols ) {
      // allocate for QR factorization
      m_QR2.no_allocate( m_ncols, m_rank, m_LWorkQR2, m_WorkQR2 );
      // perform QR factorization
      ok = m_QR2.factorize_nodim( m_Rt, m_ncols );
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

    if ( L_mm_work < m_nrows ) {
      L_mm_work = m_nrows;
      m_allocReals.free();
      m_allocReals.allocate( size_t(L_mm_work) );
      mm_work = m_allocReals( size_t(L_mm_work) );
    }

    // A*P = Q*R --> A = Q * R * P^(-1)
    // A^+ = P R^(-1) Q^T
    copy( m_nrows, b, incb, mm_work, 1 );

    m_QR1.Qt_mul( mm_work );

    if ( m_rank < m_ncols ) {
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
      std::fill_n( mm_work+m_rank, m_nrows-m_rank, 0 );
      m_QR2.Q_mul( mm_work );
    } else {
      m_QR1.invR_mul( mm_work );
    }

    m_QR1.permute( mm_work );

    copy( m_ncols, mm_work, 1, x, incx );
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

    integer NR = m_nrows*nrhs;
    if ( L_mm_work < NR ) {
      L_mm_work = NR;
      m_allocReals.free();
      m_allocReals.allocate( size_t(L_mm_work) );
      mm_work = m_allocReals( size_t(L_mm_work) );
    }

    gecopy( m_nrows, nrhs, B, ldB, mm_work, m_nrows );

    // A*P = Q*R --> A =  Q * R * P^(-1)
    // A^+ = P R^(-1) Q^T

    m_QR1.Qt_mul( m_nrows, nrhs, mm_work, m_nrows );

    if ( m_rank < m_ncols ) {
      m_QR2.invRt_mul( m_rank, nrhs, mm_work, m_nrows );
      gezero( m_nrows-m_rank, nrhs, mm_work+m_rank, m_nrows );
      m_QR2.Q_mul( m_ncols, nrhs, mm_work, m_nrows );
    } else {
      m_QR1.invR_mul( m_nrows, nrhs, mm_work, m_nrows );
    }

    m_QR1.permute_rows( m_ncols, nrhs, mm_work, m_nrows );

    gecopy( m_ncols, nrhs, mm_work, m_nrows, X, ldX );

    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  PINV_no_alloc<T>::t_mult_inv(
    valueType const b[],
    integer         incb,
    valueType       x[],
    integer         incx
  ) const {

    if ( L_mm_work < m_nrows ) {
      L_mm_work = m_nrows;
      m_allocReals.free();
      m_allocReals.allocate( size_t(L_mm_work) );
      mm_work = m_allocReals( size_t(L_mm_work) );
    }

    copy( m_ncols, b, incb, mm_work, 1 );

    m_QR1.inv_permute( mm_work );

    if ( m_rank < m_ncols ) {
      m_QR2.Qt_mul( mm_work );
      m_QR2.invR_mul( mm_work );
      std::fill_n( mm_work+m_rank, m_nrows-m_rank, 0 );
    } else {
      m_QR1.invRt_mul( mm_work );
    }

    m_QR1.Q_mul( mm_work );

    copy( m_nrows, mm_work, 1, x, incx );

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

    integer NR = m_nrows*nrhs;
    if ( L_mm_work < NR ) {
      L_mm_work = NR;
      m_allocReals.free();
      m_allocReals.allocate( size_t(L_mm_work) );
      mm_work = m_allocReals( size_t(L_mm_work) );
    }

    gecopy( m_ncols, nrhs, B, ldB, mm_work, m_nrows );

    m_QR1.inv_permute_rows( m_ncols, nrhs, mm_work, m_nrows );

    if ( m_rank < m_ncols ) {
      m_QR2.Qt_mul( m_ncols, nrhs, mm_work, m_nrows );
      gezero( m_nrows-m_rank, nrhs, mm_work+m_rank, m_nrows );
      m_QR2.invR_mul( m_rank, nrhs, mm_work, m_nrows );
    } else {
      m_QR1.invRt_mul( m_ncols, nrhs, mm_work, m_nrows );
    }

    m_QR1.Q_mul( m_nrows, nrhs, mm_work, m_nrows );

    gecopy( m_nrows, nrhs, mm_work, m_nrows, X, ldX );

    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  PINV_no_alloc<T>::info( ostream_type & stream ) const {
    Malloc<integer> memi("PINV_no_alloc<T>::info");
    memi.allocate( size_t(m_ncols) );
    integer * jpvt = memi( size_t(m_ncols) );
    m_QR1.getPerm( jpvt );
    fmt::print( stream,
      "PINV INFO\n"
      "size = {} x {}\n"
      "rank = {}\n"
      "perm = {}\n"
      "R\n{}\n",
      m_nrows, m_ncols, m_rank,
      print_matrix_transpose2( m_ncols, 1, jpvt, m_ncols ),
      print_matrix_transpose( m_ncols, m_ncols, m_Rt, m_ncols )
    );

    if ( m_rank < m_ncols ) {
      Malloc<valueType> mem("PINV_no_alloc<T>::info");
      size_t dim = size_t( m_rank * m_rank );
      mem.allocate( dim );
      valueType * R1 = mem( dim );
      m_QR2.getR( R1, m_rank );
      fmt::print( stream,
        "R1\n{}\n",
        print_matrix( m_rank, m_rank, R1, m_rank )
      );
    }
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
