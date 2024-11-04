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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  PINV_no_alloc<T>::get_required_allocation(
    integer   NR,
    integer   NC,
    integer & Lwork,
    integer & Liwork
  ) const {
    integer L1 = m_QRP1.get_Lwork_QRP( NR, NC ) + (NR+1)*(NC+1);
    integer L2 = m_QR2.get_Lwork_QR( NC, NC )   + (NC+1)*(NC+1);
    Lwork  = L1 + L2 + (NR+NC)*NC;
    Liwork = NC;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  PINV_no_alloc<T>::no_allocate(
    integer     NR,
    integer     NC,
    integer     Lwork,
    real_type * Work,
    integer     Liwork,
    integer   * iWork
  ) {

    UTILS_ASSERT(
      NR >= NC,
      "PINV_no_alloc::no_allocate( NR = {}, NC = {},...)\n"
      "must be NR >= NC\n",
      NR, NC
    );

    m_LWorkQRP1 = m_QRP1.get_Lwork_QRP( NR, NC ) + (NR+1)*(NC+1);
    m_LWorkQR2  = m_QR2.get_Lwork_QR( NC, NC )   + (NC+1)*(NC+1);
    integer Ltot = m_LWorkQRP1 + m_LWorkQR2 + (NR+NC)*NC;
    UTILS_ASSERT(
      Lwork >= Ltot && Liwork >= NC,
      "PINV_no_alloc::no_allocate( NR = {}, NC = {}, Lwork = {},..., Liwork = {},...)\n"
      "Lwork must be >= {} and Liwork >= {}\n",
      NR, NC, Lwork, Liwork, Ltot, NC
    );
    real_type * ptr = Work;
    m_A_factored = ptr; ptr += NR*NC;
    m_Rt         = ptr; ptr += NC*NC;
    m_QRP1.no_allocate( NR, NC, m_LWorkQRP1, ptr, Liwork, iWork ); ptr += m_LWorkQRP1;
    m_WorkQR2    = ptr;

    m_nrows = NR;
    m_ncols = NC;
    m_rank  = 0;
    m_rcmax = NR > NC ? NR : NC;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  PINV_no_alloc<T>::factorize_nodim(
    char      const who[],
    real_type const A[],
    integer         LDA
  ) {
    // copy matrix and scale
    integer info = gecopy( m_nrows, m_ncols, A, LDA, m_A_factored, m_nrows );
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
    real_type const A[],
    integer         LDA
  ) {
    for ( integer i = 0; i < m_nrows; ++i )
      for ( integer j = 0; j < m_ncols; ++j )
        m_A_factored[i+j*m_nrows] = A[j+i*LDA];
    this->factorize( who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  PINV_no_alloc<T>::factorize( char const who[] ) {

    // perform QR factorization
    m_QRP1.factorize_nodim( who, m_A_factored, m_nrows );

    // evaluate rank
    m_QRP1.getRt( m_Rt, m_ncols );

    m_rank = m_ncols;
    real_type threshold = absmax( m_ncols, m_Rt, m_ncols+1 ) * m_epsi;
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
    real_type const A[],
    integer         LDA
  ) {
    // copy matrix and scale
    integer info = gecopy( m_nrows, m_ncols, A, LDA, m_A_factored, m_nrows );
    if ( info != 0 ) return false;
    return this->factorize();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  PINV_no_alloc<T>::t_factorize_nodim(
    real_type const A[],
    integer         LDA
  ) {
    for ( integer i = 0; i < m_nrows; ++i )
      for ( integer j = 0; j < m_ncols; ++j )
        m_A_factored[i+j*m_nrows] = A[j+i*LDA];
    return this->factorize();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  PINV_no_alloc<T>::factorize() {

    // perform QR factorization
    bool ok = m_QRP1.factorize_nodim( m_A_factored, m_nrows );
    if ( !ok ) return false;

    // evaluate rank
    m_QRP1.getRt( m_Rt, m_ncols );

    m_rank = m_ncols;
    real_type threshold = absmax( m_ncols, m_Rt, m_ncols+1 ) * m_epsi;
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
    real_type const b[],
    integer         incb,
    real_type       x[],
    integer         incx
  ) const {

    if ( L_mm_work < m_rcmax ) {
      L_mm_work = m_rcmax;
      mm_work   = m_alloc_work.realloc( size_t(L_mm_work) );
    }

    // A*P = Q*R --> A = Q * R * P^(-1)
    // A^+ = P R^(-1) Q^T
    copy( m_nrows, b, incb, mm_work, 1 );

    m_QRP1.Qt_mul( mm_work );

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
      std::fill_n( mm_work+m_rank, m_nrows-m_rank, real_type(0) );
      m_QR2.Q_mul( mm_work );
    } else {
      m_QRP1.invR_mul( mm_work );
    }

    m_QRP1.permute( mm_work );

    copy( m_ncols, mm_work, 1, x, incx );
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  PINV_no_alloc<T>::mult_inv(
    integer         nrhs,
    real_type const B[],
    integer         ldB,
    real_type       X[],
    integer         ldX
  ) const {

    UTILS_ASSERT(
      ldB >= m_nrows && ldX >= m_ncols,
      "PINV_no_alloc::mult_inv( nrhs = {}, B, ldB = {}, X, ldX = {}\n"
      "ldB must be >= {} and ldX must be >= {}\n",
      nrhs, ldB, ldX, m_nrows, m_ncols
    );

    integer msize = m_rcmax*nrhs;
    if ( L_mm_work < msize ) {
      L_mm_work = msize;
      mm_work   = m_alloc_work.realloc( size_t(L_mm_work) );
    }

    gecopy( m_nrows, nrhs, B, ldB, mm_work, m_nrows );

    // A*P = Q*R --> A =  Q * R * P^(-1)
    // A^+ = P R^(-1) Q^T

    m_QRP1.Qt_mul( m_nrows, nrhs, mm_work, m_nrows );

    #if 0
    fmt::print(
      "Q^T rhs\n{}",
      lapack_wrapper::print_matrix( m_nrows, nrhs, mm_work, m_nrows )
    );
    #endif

    if ( m_rank < m_ncols ) {
      m_QR2.invRt_mul( m_rank, nrhs, mm_work, m_nrows );
      gezero( m_nrows-m_rank, nrhs, mm_work+m_rank, m_nrows );
      m_QR2.Q_mul( m_ncols, nrhs, mm_work, m_nrows );
    } else {
      m_QRP1.invR_mul( m_ncols, nrhs, mm_work, m_nrows );
    }

    #if 0
    fmt::print(
    "R^(-1) Q^T rhs\n{}",
      lapack_wrapper::print_matrix( m_ncols, nrhs, mm_work, m_nrows )
    );
    #endif

    m_QRP1.permute_rows( m_ncols, nrhs, mm_work, m_nrows );

    gecopy( m_ncols, nrhs, mm_work, m_nrows, X, ldX );

    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  PINV_no_alloc<T>::t_mult_inv(
    real_type const b[],
    integer         incb,
    real_type       x[],
    integer         incx
  ) const {

    if ( L_mm_work < m_rcmax ) {
      L_mm_work = m_rcmax;
      mm_work   = m_alloc_work.realloc( size_t(L_mm_work) );
    }

    copy( m_ncols, b, incb, mm_work, 1 );
    if ( m_nrows > m_ncols ) std::fill_n( mm_work+m_ncols, m_nrows-m_ncols, real_type(0) );

    m_QRP1.inv_permute( mm_work );

    if ( m_rank < m_ncols ) {
      m_QR2.Qt_mul( mm_work );
      m_QR2.invR_mul( mm_work );
      std::fill_n( mm_work+m_rank, m_nrows-m_rank, real_type(0) );
    } else {
      m_QRP1.invRt_mul( mm_work );
    }

    m_QRP1.Q_mul( mm_work );

    copy( m_nrows, mm_work, 1, x, incx );

    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  PINV_no_alloc<T>::t_mult_inv(
    integer         nrhs,
    real_type const B[],
    integer         ldB,
    real_type       X[],
    integer         ldX
  ) const {

    UTILS_ASSERT(
      ldB >= m_ncols && ldX >= m_nrows,
      "PINV_no_alloc::t_mult_inv( nrhs = {}, B, ldB = {}, X, ldX = {}\n"
      "ldB must be >= {} and ldX must be >= {}\n",
      nrhs, ldB, ldX, m_ncols, m_nrows
    );

    integer msize = m_rcmax*nrhs;
    if ( L_mm_work < msize ) {
      L_mm_work = msize;
      mm_work   = m_alloc_work.realloc( size_t(L_mm_work) );
    }

    gecopy( m_ncols, nrhs, B, ldB, mm_work, m_nrows );
    if ( m_nrows > m_ncols )
      gezero( m_nrows-m_ncols, nrhs, mm_work+m_ncols, m_nrows );

    m_QRP1.inv_permute_rows( m_ncols, nrhs, mm_work, m_nrows );

    if ( m_rank < m_ncols ) {
      m_QR2.Qt_mul( m_ncols, nrhs, mm_work, m_nrows );
      gezero( m_nrows-m_rank, nrhs, mm_work+m_rank, m_nrows );
      m_QR2.invR_mul( m_rank, nrhs, mm_work, m_nrows );
    } else {
      m_QRP1.invRt_mul( m_ncols, nrhs, mm_work, m_nrows );
    }

    m_QRP1.Q_mul( m_nrows, nrhs, mm_work, m_nrows );

    gecopy( m_nrows, nrhs, mm_work, m_nrows, X, ldX );

    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  PINV_no_alloc<T>::info( ostream_type & stream ) const {
    Malloc<integer> memi("PINV_no_alloc<T>::info");
    integer * jpvt = memi.malloc( size_t(m_ncols) );
    m_QRP1.getPerm( jpvt );
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
      Malloc<real_type> mem("PINV_no_alloc<T>::info");
      size_t dim = size_t( m_rank * m_rank );
      real_type * R1 = mem.malloc( dim );
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
    this->no_allocate(
      NR, NC,
      Ltot,  m_allocReals.realloc( size_t(Ltot) ),
      Litot, m_allocIntegers.realloc( size_t(Litot) )
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

}

///
/// eof: pinv.cxx
///
