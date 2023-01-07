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
/// file: eig.cxx
///

namespace lapack_wrapper {

  /*\
  :|:   _____ _                            _
  :|:  | ____(_) __ _  ___ _ ____   ____ _| |_   _  ___  ___
  :|:  |  _| | |/ _` |/ _ \ '_ \ \ / / _` | | | | |/ _ \/ __|
  :|:  | |___| | (_| |  __/ | | \ V / (_| | | |_| |  __/\__ \
  :|:  |_____|_|\__, |\___|_| |_|\_/ \__,_|_|\__,_|\___||___/
  :|:           |___/
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Eigenvalues<T>::allocate( integer Nin ) {
    m_N = Nin;
    // calcolo memoria ottimale
    real_type Lworkdummy;
    integer info = geev(
      false, false, Nin, nullptr, Nin,
      nullptr, nullptr, nullptr, Nin, nullptr, Nin, &Lworkdummy, -1
    );
    UTILS_ASSERT(
      info == 0, "Eigenvalues<T>::allocate, call geev return info = {}\n", info
    );
    m_Lwork = integer( Lworkdummy );
    m_mem.reallocate( size_t( m_Lwork + (2+m_N) * m_N) );
    m_Re      = m_mem( size_t(m_N) );
    m_Im      = m_mem( size_t(m_N) );
    m_Work    = m_mem( size_t(m_Lwork) );
    m_A_saved = m_mem( size_t(m_N*m_N) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Eigenvalues<T>::compute( ) {
    integer info = geev(
      false, false, m_N, m_A_saved, m_N,
      m_Re, m_Im, nullptr, m_N, nullptr, m_N,
      m_Work, m_Lwork
    );
    UTILS_ASSERT(
      info == 0, "Eigenvalues<T>::compute, call geev return info = {}\n", info
    );
  }


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Eigenvalues<T>::setup(
    integer         NRC,
    real_type const data[],
    integer         ldData
  ) {
    this->allocate( NRC );
    integer info = gecopy( NRC, NRC, data, ldData, m_A_saved, NRC );
    UTILS_ASSERT(
      info == 0, "Eigenvalues<T>::setup, call gecopy return info = {}\n", info
    );
    this->compute();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Eigenvalues<T>::setup( MatW const & M ) {
    this->allocate( M.nrows() );
    integer info = gecopy( m_N, m_N, M.data(), M.ldim(), m_A_saved, m_N );
    UTILS_ASSERT(
      info == 0, "Eigenvalues<T>::setup, call gecopy return info = {}\n", info
    );
    this->compute();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Eigenvalues<T>::setup(
    integer         NRC,
    integer         nnz,
    real_type const values[],
    integer   const row[],
    integer   const col[]
  ) {
    this->allocate( NRC );
    lapack_wrapper::zero( NRC*NRC, m_A_saved, 1 );
    for ( integer i = 0; i < nnz; ++i )
      m_A_saved[row[i]+col[i]*NRC] += values[i];
    this->compute();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Eigenvalues<T>::getEigenvalue(
    integer n, real_type & re, real_type & im
  ) const {
    re = m_Re[n];
    im = m_Im[n];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Eigenvalues<T>::getEigenvalue(
    integer n, std::complex<real_type> & eig
  ) const {
    eig = std::complex<real_type>( m_Re[n], m_Im[n] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Eigenvalues<T>::getEigenvalues(
    std::vector<real_type> & re,
    std::vector<real_type> & im
  ) const {
    re.clear(); re.reserve( m_N );
    im.clear(); im.reserve( m_N );
    for ( int i = 0; i < m_N; ++i ) {
      re.push_back( m_Re[i] );
      im.push_back( m_Im[i] );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Eigenvalues<T>::getEigenvalues(
    std::vector<std::complex<real_type> > & eigs
  ) const {
    eigs.clear(); eigs.reserve( m_N );
    for ( int i = 0; i < m_N; ++i )
      eigs.push_back( std::complex<real_type>( m_Re[i], m_Im[i]) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Eigenvectors<T>::allocate( integer Nin ) {
    m_N = Nin;
    // calcolo memoria ottimale
    integer doLwork    = -1;
    T       Lworkdummy = 1;
    integer info = geev(
      true, true,
      m_N,
      nullptr, m_N,
      nullptr, nullptr,
      m_VL, m_N,
      m_VR, m_N,
      &Lworkdummy, doLwork
    );
    UTILS_ASSERT(
      info == 0, "Eigenvectors::allocate, call geev return info = {}\n", info
    );
    m_Lwork = integer(Lworkdummy);
    m_mem.reallocate( size_t( m_Lwork + (2+3*m_N) * m_N) );
    m_Re      = m_mem( size_t(m_N) );
    m_Im      = m_mem( size_t(m_N) );
    m_A_saved = m_mem( size_t(m_N*m_N) );
    m_VL      = m_mem( size_t(m_N*m_N) );
    m_VR      = m_mem( size_t(m_N*m_N) );
    m_Work    = m_mem( size_t(m_Lwork) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Eigenvectors<T>::compute( ) {
    integer info = geev(
      m_VL != nullptr,
      m_VR != nullptr,
      m_N,
      m_A_saved, m_N,
      m_Re,      m_Im,
      m_VL,      m_N,
      m_VR,      m_N,
      m_Work,    m_Lwork
    );
    UTILS_ASSERT(
      info == 0,
      "GeneralizedEigenvectors::compute, call geev return info = {}\n", info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Eigenvectors<T>::setup(
    integer         NRC,
    real_type const A_data[],
    integer         ldA
  ) {
    this->allocate( NRC );
    integer info = gecopy( NRC, NRC, A_data, ldA, m_A_saved, NRC );
    UTILS_ASSERT(
      info == 0,
      "Eigenvectors::setup, call gecopy return info = {}\n", info
    );
    this->compute();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Eigenvectors<T>::setup( MatW const & A ) {
    this->allocate( A.nrows() );
    integer info = gecopy( m_N, m_N, A.data(), A.ldim(),  m_A_saved, m_N );
    UTILS_ASSERT(
      info == 0,
      "Eigenvectors::setup, call gecopy return info = {}\n", info
    );
    this->compute();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Eigenvectors<T>::setup(
    integer         NRC,
    integer         A_nnz,
    real_type const A_values[],
    integer   const A_row[],
    integer   const A_col[]
  ) {
    this->allocate( NRC );
    lapack_wrapper::zero( NRC*NRC, m_A_saved, 1 );
    for ( integer i = 0; i < A_nnz; ++i )
      m_A_saved[A_row[i]+A_col[i]*NRC] += A_values[i];
    this->compute();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Eigenvectors<T>::getEigenvalue(
    integer n, real_type & re, real_type & im
  ) const {
    re = m_Re[n];
    im = m_Im[n];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Eigenvectors<T>::getEigenvalue(
    integer n, std::complex<real_type> & eig
  ) const {
    eig = std::complex<real_type>( m_Re[n], m_Im[n] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Eigenvectors<T>::getEigenvalues(
    std::vector<real_type> & re,
    std::vector<real_type> & im
  ) const {
    re.clear(); re.reserve( m_N );
    im.clear(); im.reserve( m_N );
    for ( int i = 0; i < m_N; ++i ) {
      re.push_back( m_Re[i] );
      im.push_back( m_Im[i] );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Eigenvectors<T>::getEigenvalues(
    std::vector<std::complex<real_type> > & eigs
  ) const {
    eigs.clear(); eigs.reserve( m_N );
    for ( int i = 0;i < m_N; ++i )
      eigs.push_back( std::complex<real_type>( m_Re[i], m_Im[i]) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Eigenvectors<T>::getLeftEigenvector(
    std::vector<std::vector<complex_type> > & vecs
  ) const {
    vecs.resize( size_t(m_N) );
    for ( integer n = 0; n < m_N; ++n ) {
      std::vector<complex_type> & v = vecs[n];
      v.clear(); v.reserve( m_N );
      T const * vr = m_VL + n * m_N;
      if ( m_Im[n] > 0 ) {
        std::vector<complex_type> & v1 = vecs[++n];
        v1.clear(); v1.reserve( m_N );
        T const * vi = vr + m_N;
        for ( integer j = 0; j < m_N; ++j ) {
          // salvo vettore "gia" coniugato
          v.push_back( complex_type( vr[j], -vi[j] ) );
          v1.push_back( complex_type( vr[j], vi[j] ) );
        }
      } else {
        for ( integer j = 0; j < m_N; ++j )
          v.push_back( complex_type( vr[j], 0 ) );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Eigenvectors<T>::getRightEigenvector(
    std::vector<std::vector<complex_type> > & vecs
  ) const {
    vecs.resize( size_t(m_N) );
    for ( integer n = 0; n < m_N; ++n ) {
      std::vector<complex_type> & v = vecs[n];
      v.clear(); v.reserve( m_N );
      T const * vr = m_VR + n * m_N;
      if ( m_Im[n] > 0 ) {
        std::vector<complex_type> & v1 = vecs[++n];
        v1.clear(); v1.reserve( m_N );
        T const * vi = vr + m_N;
        for ( integer j = 0; j < m_N; ++j ) {
          v.push_back( complex_type( vr[j], vi[j] ) );
          v1.push_back( complex_type( vr[j], -vi[j] ) );
        }
      } else {
        for ( integer j = 0; j < m_N; ++j )
          v.push_back( complex_type( vr[j], 0 ) );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedEigenvalues<T>::allocate( integer Nin ) {
    m_N = Nin;
    // calcolo memoria ottimale
    real_type Lworkdummy;
    integer info = ggev(
      false, false, Nin, nullptr, Nin, nullptr, Nin,
      nullptr, nullptr, nullptr,
      nullptr, Nin, nullptr, Nin, &Lworkdummy, -1
    );
    UTILS_ASSERT(
      info == 0,
      "GeneralizedEigenvalues::allocate, call geev return info = {}\n", info
    );
    m_Lwork = integer(Lworkdummy);
    m_mem.reallocate( size_t( m_Lwork + (3+2*m_N) * m_N) );
    m_alphaRe = m_mem( size_t(m_N) );
    m_alphaIm = m_mem( size_t(m_N) );
    m_beta    = m_mem( size_t(m_N) );
    m_Work    = m_mem( size_t(m_Lwork) );
    m_A_saved = m_mem( size_t(m_N*m_N) );
    m_B_saved = m_mem( size_t(m_N*m_N) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedEigenvalues<T>::compute( ) {
    integer info = ggev(
      false, false,
      m_N,
      m_A_saved, m_N,
      m_B_saved, m_N,
      m_alphaRe, m_alphaIm, m_beta,
      nullptr, m_N, nullptr, m_N,
      m_Work, m_Lwork
    );
    UTILS_ASSERT(
      info == 0,
      "GeneralizedEigenvalues::compute, call geev return info = {}\n", info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedEigenvalues<T>::setup(
    integer         NRC,
    real_type const A_data[],
    integer         ldA,
    real_type const B_data[],
    integer         ldB
  ) {
    this->allocate( NRC );
    integer info1 = gecopy( NRC, NRC, A_data, ldA, m_A_saved, NRC );
    integer info2 = gecopy( NRC, NRC, B_data, ldB, m_B_saved, NRC );
    UTILS_ASSERT(
      info1 == 0 && info2 == 0,
      "GeneralizedEigenvalues::setup, call gecopy return info1 = {}, info2 = {}\n",
      info1, info2
    );
    this->compute();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedEigenvalues<T>::setup( MatW const & A, MatW const & B ) {
    this->allocate( A.nrows() );
    integer info1 = gecopy( m_N, m_N, A.data(), A.ldim(), m_A_saved, m_N );
    integer info2 = gecopy( m_N, m_N, B.data(), B.ldim(), m_B_saved, m_N );
    UTILS_ASSERT(
      info1 == 0 && info2 == 0,
      "GeneralizedEigenvalues::setup, call gecopy return info1 = {}, info2 = {}\n",
      info1, info2
    );
    this->compute();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedEigenvalues<T>::setup(
    integer         NRC,
    integer         A_nnz,
    real_type const A_values[],
    integer   const A_row[],
    integer   const A_col[],
    integer         B_nnz,
    real_type const B_values[],
    integer   const B_row[],
    integer   const B_col[]
  ) {
    this->allocate( NRC );
    lapack_wrapper::zero( NRC*NRC, m_A_saved, 1 );
    lapack_wrapper::zero( NRC*NRC, m_B_saved, 1 );
    for ( integer i = 0; i < A_nnz; ++i )
      m_A_saved[A_row[i]+A_col[i]*NRC] += A_values[i];
    for ( integer i = 0; i < B_nnz; ++i )
      m_B_saved[B_row[i]+B_col[i]*NRC] += B_values[i];
    this->compute();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedEigenvalues<T>::getEigenvalue(
    integer n, real_type & re, real_type & im
  ) const {
    re = m_alphaRe[n]/m_beta[n];
    im = m_alphaIm[n]/m_beta[n];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedEigenvalues<T>::getEigenvalue(
    integer n, std::complex<real_type> & eig
  ) const {
    eig = std::complex<real_type>( m_alphaRe[n], m_alphaIm[n] ) / m_beta[n];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedEigenvalues<T>::getEigenvalues(
    std::vector<real_type> & re,
    std::vector<real_type> & im
  ) const {
    re.clear(); re.reserve( m_N );
    im.clear(); im.reserve( m_N );
    for ( int i = 0;i < m_N; ++i ) {
      re.push_back( m_alphaRe[i]/m_beta[i] );
      im.push_back( m_alphaIm[i]/m_beta[i] );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedEigenvalues<T>::getEigenvalues(
    std::vector<std::complex<real_type> > & eigs
  ) const {
    eigs.clear(); eigs.reserve( m_N );
    for ( int i = 0;i < m_N; ++i )
      eigs.push_back(
        std::complex<real_type>( m_alphaRe[i], m_alphaIm[i] ) / m_beta[i]
      );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedEigenvectors<T>::allocate( integer Nin ) {
    m_N = Nin;
    // calcolo memoria ottimale
    integer doLwork    = -1;
    T       Lworkdummy = 1;
    integer info = ggevx(
      BalanceType::PERMUTE_AND_SCALE,
      false, false,
      SenseType::EIGENVALUES_AND_EIGENVECTORS,
      m_N, nullptr, m_N, nullptr, m_N,
      nullptr, nullptr, nullptr, m_VL, m_N, m_VR, m_N,
      m_ilo, m_ihi,
      nullptr, nullptr,
      m_abnorm, m_bbnorm,
      nullptr, nullptr,
      &Lworkdummy, doLwork,
      nullptr, nullptr
    );
    UTILS_ASSERT(
      info == 0,
      "GeneralizedEigenvectors::allocate, call geev return info = {}\n", info
    );
    m_Lwork = integer(Lworkdummy);
    m_mem_real.reallocate( size_t( m_Lwork + (7+4*m_N) * m_N) );
    m_mem_int.reallocate( size_t( 2*m_N + 6 ) );
    m_alphaRe = m_mem_real( size_t(m_N) );
    m_alphaIm = m_mem_real( size_t(m_N) );
    m_beta    = m_mem_real( size_t(m_N) );
    m_A_saved = m_mem_real( size_t(m_N*m_N) );
    m_B_saved = m_mem_real( size_t(m_N*m_N) );
    m_VL      = m_mem_real( size_t(m_N*m_N) );
    m_VR      = m_mem_real( size_t(m_N*m_N) );
    m_lscale  = m_mem_real( size_t(m_N) );
    m_rscale  = m_mem_real( size_t(m_N) );
    m_rconde  = m_mem_real( size_t(m_N) );
    m_rcondv  = m_mem_real( size_t(m_N) );
    m_Work    = m_mem_real( size_t(m_Lwork) );
    m_iWork   = m_mem_int( size_t(m_N+6) );
    m_bWork   = m_mem_int( size_t(m_N) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedEigenvectors<T>::compute( ) {
    integer info = ggevx(
      BalanceType::PERMUTE_ONLY, // ATTENZIONE LA SCALATURA NON FUNZIONA
      m_VL != nullptr,
      m_VR != nullptr,
      SenseType::EIGENVALUES_AND_EIGENVECTORS,
      m_N,
      m_A_saved, m_N,
      m_B_saved, m_N,
      m_alphaRe, m_alphaIm, m_beta,
      m_VL,      m_N,
      m_VR,      m_N,
      m_ilo,     m_ihi,
      m_lscale,  m_rscale,
      m_abnorm,  m_bbnorm,
      m_rconde,  m_rcondv,
      m_Work,    m_Lwork,
      m_iWork,   m_bWork
    );
    UTILS_ASSERT(
      info == 0,
      "GeneralizedEigenvectors::compute, call ggevx return info = {}\n", info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedEigenvectors<T>::setup(
    integer         NRC,
    real_type const A_data[],
    integer         ldA,
    real_type const B_data[],
    integer         ldB
  ) {
    this->allocate( NRC );
    integer info1 = gecopy( NRC, NRC, A_data, ldA, m_A_saved, NRC );
    integer info2 = gecopy( NRC, NRC, B_data, ldB, m_B_saved, NRC );
    UTILS_ASSERT(
      info1 == 0 && info2 == 0,
      "GeneralizedEigenvectors::setup, call gecopy return info1 = {}, info2 = {}\n",
      info1, info2
    );
    this->compute();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedEigenvectors<T>::setup( MatW const & A, MatW const & B ) {
    this->allocate( A.nrows() );
    integer info1 = gecopy( m_N, m_N, A.data(), A.ldim(), m_A_saved, m_N );
    integer info2 = gecopy( m_N, m_N, B.data(), B.ldim(), m_B_saved, m_N );
    UTILS_ASSERT(
      info1 == 0 && info2 == 0,
      "GeneralizedEigenvectors::setup, call gecopy return info1 = {}, info2 = {}\n",
      info1, info2
    );
    this->compute();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedEigenvectors<T>::setup(
    integer         NRC,
    integer         A_nnz,
    real_type const A_values[],
    integer   const A_row[],
    integer   const A_col[],
    integer         B_nnz,
    real_type const B_values[],
    integer   const B_row[],
    integer   const B_col[]
  ) {
    this->allocate( NRC );
    lapack_wrapper::zero( NRC*NRC, m_A_saved, 1 );
    lapack_wrapper::zero( NRC*NRC, m_B_saved, 1 );
    for ( integer i = 0; i < A_nnz; ++i )
      m_A_saved[A_row[i]+A_col[i]*NRC] += A_values[i];
    for ( integer i = 0; i < B_nnz; ++i )
      m_B_saved[B_row[i]+B_col[i]*NRC] += B_values[i];
    this->compute();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedEigenvectors<T>::getEigenvalue(
    integer n, real_type & re, real_type & im
  ) const {
    re = m_alphaRe[n] / m_beta[n];
    im = m_alphaIm[n] / m_beta[n];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedEigenvectors<T>::getEigenvalue(
    integer n, std::complex<real_type> & eig
  ) const {
    eig = std::complex<real_type>( m_alphaRe[n], m_alphaIm[n] ) / m_beta[n];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedEigenvectors<T>::getEigenvalues(
    std::vector<real_type> & re,
    std::vector<real_type> & im
  ) const {
    re.clear(); re.reserve( m_N );
    im.clear(); im.reserve( m_N );
    for ( int i = 0; i < m_N; ++i ) {
      re.push_back( m_alphaRe[i] / m_beta[i] );
      im.push_back( m_alphaIm[i] / m_beta[i] );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedEigenvectors<T>::getEigenvalues(
    std::vector<std::complex<real_type> > & eigs
  ) const {
    eigs.clear(); eigs.reserve( m_N );
    for ( int i = 0; i < m_N; ++i )
      eigs.push_back(
        std::complex<real_type>( m_alphaRe[i], m_alphaIm[i] ) / m_beta[i]
      );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedEigenvectors<T>::getLeftEigenvector(
    std::vector<std::vector<complex_type> > & vecs
  ) const {
    vecs.resize( size_t(m_N) );
    for ( integer n = 0; n < m_N; ++n ) {
      std::vector<complex_type> & v = vecs[n];
      v.clear(); v.reserve( size_t(m_N) );
      T const * vr = m_VL + n * m_N;
      if ( m_alphaIm[n] > 0 ) {
        std::vector<complex_type> & v1 = vecs[++n];
        v1.clear(); v1.reserve( m_N );
        T const * vi = vr + m_N;
        for ( integer j = 0; j < m_N; ++j ) {
          // salvo vettore "gia" coniugato
          v.push_back( complex_type( vr[j], -vi[j] ) );
          v1.push_back( complex_type( vr[j], vi[j] ) );
        }
      } else {
        for ( integer j = 0; j < m_N; ++j )
          v.push_back( complex_type( vr[j], 0 ) );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  GeneralizedEigenvectors<T>::getRightEigenvector(
    std::vector<std::vector<complex_type> > & vecs
  ) const {
    vecs.resize( size_t(m_N) );
    for ( integer n = 0; n < m_N; ++n ) {
      std::vector<complex_type> & v = vecs[n];
      v.clear(); v.reserve( size_t(m_N) );
      T const * vr = m_VR + n * m_N;
      if ( m_alphaIm[n] > 0 ) {
        std::vector<complex_type> & v1 = vecs[++n];
        v1.clear(); v1.reserve( m_N );
        T const * vi = vr + m_N;
        for ( integer j = 0; j < m_N; ++j ) {
          v.push_back( complex_type( vr[j], vi[j] ) );
          v1.push_back( complex_type( vr[j], -vi[j] ) );
        }
      } else {
        for ( integer j = 0; j < m_N; ++j )
          v.push_back( complex_type( vr[j], 0 ) );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

}

///
/// eof: eig.cxx
///
