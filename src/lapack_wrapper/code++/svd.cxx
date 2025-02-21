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
/// file: svd.cxx
///

namespace lapack_wrapper {

  /*\
   |   ______     ______
   |  / ___\ \   / /  _ \
   |  \___ \\ \ / /| | | |
   |   ___) |\ V / | |_| |
   |  |____/  \_/  |____/
   |
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  integer
  SVD_no_alloc<T>::get_Lwork( integer NR, integer NC ) {
    integer minRC = std::min(NR,NC);
    real_type tmp;
    integer info = gesvd(
      JobType::REDUCED,
      JobType::REDUCED,
      NR, NC,
      nullptr, NR,
      nullptr,
      nullptr, NR,
      nullptr, minRC,
      &tmp, -1
    );
    UTILS_ASSERT( info == 0, "SVD::allocate, in gesvd info = {}\n", info );
    integer L = integer(tmp);
    info = gesdd(
      JobType::REDUCED,
      NR, NC,
      nullptr, NR,
      nullptr,
      nullptr, NR,
      nullptr, minRC,
      &tmp, -1, nullptr
    );
    UTILS_ASSERT( info == 0, "SVD::allocate, in gesdd info = {}\n", info );
    if ( integer(tmp) > L ) L = integer(tmp);
    return L;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  SVD_no_alloc<T>::no_allocate(
    integer     NR,
    integer     NC,
    integer     Lwork,
    real_type * Work,
    integer     Liwork,
    integer   * iWork
  ) {
    m_nrows = NR;
    m_ncols = NC;
    m_minRC = std::min( NR, NC );
    integer ibf   = (m_minRC*(NR+NC+1)+NR*NC);
    integer Lmin  = ibf+this->get_Lwork( NR, NC );
    integer Limin = 8*m_minRC;
    m_LworkSVD = Lwork - ibf;
    UTILS_ASSERT(
      Lwork >= Lmin && Liwork >= Limin,
      "SVD_no_alloc::no_allocate( NR = {}, NC = {}, Lwork = {},..., Liwork = {},...)\n"
      "Lwork must be >= {} and Liwork >= {}\n",
      NR, NC, Lwork, Liwork, Lmin, Limin
    );
    real_type * ptr = Work;
    m_Afactorized = ptr; ptr += NR*NC;
    m_Umat        = ptr; ptr += m_minRC*NR;
    m_VTmat       = ptr; ptr += m_minRC*NC;
    m_Svec        = ptr; ptr += m_minRC;
    m_WorkSVD     = ptr;
    m_IWorkSVD    = iWork;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  SVD_no_alloc<T>::factorize_nodim(
    string_view     who,
    real_type const A[],
    integer         LDA
  ) {
    integer info = gecopy( m_nrows, m_ncols, A, LDA, m_Afactorized, m_nrows );
    UTILS_ASSERT(
      info == 0,
      "SVD_no_alloc::factorize[{}] call lapack_wrapper::gecopy return info = {}\n",
      who, info
    );
    switch ( m_svd_used ) {
    case USE_GESVD:
      info = gesvd(
        JobType::REDUCED,
        JobType::REDUCED,
        m_nrows, m_ncols, m_Afactorized, m_nrows,
        m_Svec,
        m_Umat,    m_nrows,
        m_VTmat,   m_minRC,
        m_WorkSVD, m_LworkSVD
      );
      UTILS_ASSERT(
        info == 0,
        "SVD_no_alloc::factorize[{}] call lapack_wrapper::gesvd return info = {}\n",
        who, info
      );
      break;
    case USE_GESDD:
      info = gesdd(
        JobType::REDUCED,
        m_nrows, m_ncols, m_Afactorized, m_nrows,
        m_Svec,
        m_Umat,    m_nrows,
        m_VTmat,   m_minRC,
        m_WorkSVD, m_LworkSVD, m_IWorkSVD
      );
      UTILS_ASSERT(
        info == 0,
        "SVD_no_alloc::factorize[{}] call lapack_wrapper::gesdd return info = {}\n",
        who, info
      );
      break;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  SVD_no_alloc<T>::factorize_nodim( real_type const A[], integer LDA ) {
    integer info = gecopy(
      m_nrows, m_ncols, A, LDA, m_Afactorized, m_nrows
    );
    if ( info != 0 ) return false;
    switch ( m_svd_used ) {
    case USE_GESVD:
      info = gesvd(
        JobType::REDUCED,
        JobType::REDUCED,
        m_nrows, m_ncols, m_Afactorized, m_nrows,
        m_Svec,
        m_Umat,    m_nrows,
        m_VTmat,   m_minRC,
        m_WorkSVD, m_LworkSVD
      );
      break;
    case USE_GESDD:
      info = gesdd(
        JobType::REDUCED,
        m_nrows, m_ncols, m_Afactorized, m_nrows,
        m_Svec,
        m_Umat,    m_nrows,
        m_VTmat,   m_minRC,
        m_WorkSVD, m_LworkSVD, m_IWorkSVD
      );
      break;
    }
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  SVD_no_alloc<T>::solve( real_type xb[] ) const {
    // A = U*S*VT
    // U*S*VT*x=b --> VT^T S^+ U^T b
    // U  nRow x minRC
    // VT minRC x nCol
    real_type smin = m_rcond*m_Svec[0];
    Ut_mul( 1.0, xb, 1, 0.0, m_WorkSVD, 1 );
    for ( integer i = 0; i < m_minRC; ++i ) m_WorkSVD[i] /= std::max(m_Svec[i],smin);
    V_mul( 1.0, m_WorkSVD, 1, 0.0, xb, 1 );
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  SVD_no_alloc<T>::t_solve( real_type xb[] ) const {
    // A = U*S*VT
    // U*S*VT*x=b --> VT^T S^+ U^T b
    // U  nRow x minRC
    // VT minRC x nCol
    real_type smin = m_rcond*m_Svec[0];
    Vt_mul( 1.0, xb, 1, 0.0, m_WorkSVD, 1 );
    for ( integer i = 0; i < m_minRC; ++i ) m_WorkSVD[i] /= std::max(m_Svec[i],smin);
    U_mul( 1.0, m_WorkSVD, 1, 0.0, xb, 1 );
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  SVD<T>::allocate( integer NR, integer NC ) {
    if ( m_nrows != NR || m_ncols != NC ) {
      integer minRC = std::min(NR,NC);
      integer L     = NR*NC+minRC*(NR+NC+1)+this->get_Lwork( NR, NC );
      integer L1    = 8*minRC;
      this->no_allocate(
        NR, NC,
        L,  m_allocReals.realloc( size_t(L) ),
        L1, m_allocIntegers.realloc( size_t(L1) )
      );
    }
  }

  /*\
  :|:
  :|:    ____                           _ _             _ ______     ______
  :|:   / ___| ___ _ __   ___ _ __ __ _| (_)_______  __| / ___\ \   / /  _ \
  :|:  | |  _ / _ \ '_ \ / _ \ '__/ _` | | |_  / _ \/ _` \___ \\ \ / /| | | |
  :|:  | |_| |  __/ | | |  __/ | | (_| | | |/ /  __/ (_| |___) |\ V / | |_| |
  :|:   \____|\___|_| |_|\___|_|  \__,_|_|_/___\___|\__,_|____/  \_/  |____/
  :|:
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedSVD<T>::allocate( integer m, integer n, integer p ) {
    integer k, l;
    real    wL;
    integer info = ggsvd(
      true, true, true, m, n, p, k, l,
      nullptr, m,
      nullptr, p,
      nullptr,
      nullptr,
      nullptr, m,
      nullptr, p,
      nullptr, n,
      &wL,
      -1,
      nullptr
    );
    UTILS_ASSERT(
      info == 0,
      "GeneralizedSVD<T>::allocate(m={},n={},p={}) failed, info = {}\n",
      m, n, p, info
    );
    m_M     = m;
    m_N     = n;
    m_P     = p;
    m_Lwork = integer(wL);
    m_IWork = m_mem_int.realloc( size_t(n) );

    m_mem_real.reallocate( m_Lwork + (m+p+2)*n + m*m + p*p + n*n );
    m_Work        = m_mem_real( m_Lwork );
    m_alpha_saved = m_mem_real( n );
    m_beta_saved  = m_mem_real( n );
    m_A_saved     = m_mem_real( m*n );
    m_B_saved     = m_mem_real( p*n );
    m_U_saved     = m_mem_real( m*m );
    m_V_saved     = m_mem_real( p*p );
    m_Q_saved     = m_mem_real( n*n );

    m_U.setup( m_U_saved, m, m, m );
    m_V.setup( m_V_saved, p, p, p );
    m_Q.setup( m_Q_saved, n, n, n );
    m_Dbeta.setup( m_beta_saved, n );
    m_Dalpha.setup( m_alpha_saved, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedSVD<T>::compute() {
    integer info = ggsvd(
      true, true, true,
      m_M, m_N, m_P, m_K, m_L,
      m_A_saved, m_M,
      m_B_saved, m_P,
      m_alpha_saved,
      m_beta_saved,
      m_U_saved, m_M,
      m_V_saved, m_P,
      m_Q_saved, m_N,
      m_Work,
      m_Lwork,
      m_IWork
    );
    UTILS_ASSERT(
      info == 0, "GeneralizedSVD<T>::compute() failed, info = {}\n", info
    );
    // R is stored in A(1:K+L,N-K-L+1:N) on exit.
    m_R.setup(
      m_A_saved + m_M * ( m_N-m_K-m_L ),
      m_N,
      m_K+m_L,
      m_M
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedSVD<T>::setup(
    integer         m,
    integer         n,
    integer         p,
    real_type const A[], integer ldA,
    real_type const B[], integer ldB
  ) {
    this->allocate( m, n, p );
    integer info = gecopy( m, n, A, ldA, m_A_saved, m );
    UTILS_ASSERT(
      info == 0,
      "GeneralizedSVD<T>::setup(...) failed to copy A, info = {}\n", info
    );
    info = gecopy( p, n, B, ldB, m_B_saved, p );
    UTILS_ASSERT(
      info == 0,
      "GeneralizedSVD<T>::setup(...) failed to copy B, info = {}\n", info
    );
    compute();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedSVD<T>::setup( MatW const & A, MatW const & B ) {
    integer m = A.nrows();
    integer n = A.ncols();
    integer p = B.nrows();
    UTILS_ASSERT(
      n == B.ncols(),
      "GeneralizedSVD<T>::setup( A, B ) incompatible matrices\n"
      "A is {} x {} and B is {} x {}\n",
      A.nrows(), A.ncols(), B.nrows(), B.ncols()
    );
    this->setup( m, n, p, A.data(), A.ldim(), B.data(), B.ldim() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedSVD<T>::setup(
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
    this->allocate( m, n, p );
    lapack_wrapper::zero( m_N*m_M, m_A_saved, 1 );
    lapack_wrapper::zero( m_P*m_M, m_B_saved, 1 );
    for ( integer k = 0; k < A_nnz; ++k ) A(A_row[k],A_col[k]) = A_values[k];
    for ( integer k = 0; k < B_nnz; ++k ) B(B_row[k],B_col[k]) = B_values[k];
    compute();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  string
  GeneralizedSVD<T>::info( real_type eps ) const {
    string res = fmt::format(
      "A = {} x {}\n"
      "B = {} x {}\n",
      m_M, m_N, m_P, m_N
    );
    for ( integer i{0}; i < m_N; ++i ) {
      T a = m_alpha_saved[i];
      T b = m_beta_saved[i];
      res += fmt::format(
        "alpha[{}]={}, beta[{}]={}, alpha^2+beta^2 = {}\n",
        i, a, i, b, a*a+b*b
      );
    }
    res += fmt::format(
      "U\n{}"
      "V\n{}"
      "Q\n{}"
      "R\n{}",
      m_U.to_string( eps ),
      m_V.to_string( eps ),
      m_Q.to_string( eps ),
      m_R.to_string( eps )
    );
    return res;
  }

}

///
/// eof: svd.cxx
///
