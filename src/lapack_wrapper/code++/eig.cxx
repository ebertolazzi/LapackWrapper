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

  template <typename T>
  Eigenvalues<T>::Eigenvalues()
  : mem_real("Eigenvalues::mem_real")
  , N(0)
  , Re(nullptr)
  , Im(nullptr)
  , Work(nullptr)
  , A_saved(nullptr)
  {}

  template <typename T>
  void
  Eigenvalues<T>::allocate( integer Nin ) {
    this->N = Nin;
    // calcolo memoria ottimale
    valueType Lworkdummy;
    integer info = geev(
      false, false, Nin, nullptr, Nin,
      nullptr, nullptr, nullptr, Nin, nullptr, Nin, &Lworkdummy, -1
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "Eigenvalues<T>::allocate, call geev return info = " << info
    );
    this->Lwork = integer(Lworkdummy);
    this->mem_real.allocate( size_t( this->Lwork + (2+this->N) * this->N) );
    this->Re      = this->mem_real( size_t(this->N) );
    this->Im      = this->mem_real( size_t(this->N) );
    this->Work    = this->mem_real( size_t(this->Lwork) );
    this->A_saved = this->mem_real( size_t(this->N*this->N) );
  }

  template <typename T>
  void
  Eigenvalues<T>::compute( ) {
    integer info = geev(
      false, false, this->N, this->A_saved, this->N,
      this->Re, this->Im, nullptr, this->N, nullptr, this->N,
      this->Work, this->Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "Eigenvalues<T>::compute, call geev return info = " << info
    );
  }

  template <typename T>
  Eigenvalues<T>::Eigenvalues(
    integer         NRC,
    valueType const data[],
    integer         ldData
  )
  : mem_real("Eigenvalues::mem_real")
  , N(0)
  , Re(nullptr)
  , Im(nullptr)
  , Work(nullptr)
  , A_saved(nullptr)
  {
    this->setup( NRC, data, ldData );
  }

  template <typename T>
  Eigenvalues<T>::Eigenvalues( MatW const & M )
  : mem_real("Eigenvalues::mem_real")
  , N(0)
  , Re(nullptr)
  , Im(nullptr)
  , Work(nullptr)
  , A_saved(nullptr)
  {
    this->setup( M );
  }

  template <typename T>
  Eigenvalues<T>::Eigenvalues(
    integer         NRC,
    integer         nnz,
    valueType const values[],
    integer   const row[],
    integer   const col[]
  )
  : mem_real("Eigenvalues::mem_real")
  , N(0)
  , Re(nullptr)
  , Im(nullptr)
  , Work(nullptr)
  , A_saved(nullptr)
  {
    this->setup( NRC, nnz, values, row, col );
  }

  template <typename T>
  void
  Eigenvalues<T>::setup(
    integer         NRC,
    valueType const data[],
    integer         ldData
  ) {
    this->allocate( NRC );
    integer info = gecopy( NRC, NRC, data, ldData, A_saved, NRC );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "Eigenvalues<T>::setup, call gecopy return info = " << info
    );
    this->compute();
  }

  template <typename T>
  void
  Eigenvalues<T>::setup( MatW const & M ) {
    this->allocate( M.numRows() );
    integer info = gecopy( this->N, this->N, M.get_data(), M.lDim(),  this->A_saved, this->N );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "Eigenvalues<T>::setup, call gecopy return info = " << info
    );
    this->compute();
  }

  template <typename T>
  void
  Eigenvalues<T>::setup(
    integer         NRC,
    integer         nnz,
    valueType const values[],
    integer   const row[],
    integer   const col[]
  ) {
    this->allocate( NRC );
    lapack_wrapper::zero( NRC*NRC, this->A_saved, 1 );
    for ( integer i = 0; i < nnz; ++i )  this->A_saved[row[i]+col[i]*NRC] += values[i];
    this->compute();
  }

  template <typename T>
  void
  Eigenvalues<T>::getEigenvalue(
    integer n, valueType & re, valueType & im
  ) const {
    re = Re[n];
    im = Im[n];
  }

  template <typename T>
  void
  Eigenvalues<T>::getEigenvalue(
    integer n, std::complex<valueType> & eig
  ) const {
    eig = std::complex<valueType>(Re[n],Im[n]);
  }

  template <typename T>
  void
  Eigenvalues<T>::getEigenvalues(
    std::vector<valueType> & re,
    std::vector<valueType> & im
  ) const {
    re.clear(); re.reserve( this->N );
    im.clear(); im.reserve( this->N );
    for ( int i = 0;i < this->N; ++i ) {
      re.push_back( Re[i] );
      im.push_back( Im[i] );
    }
  }

  template <typename T>
  void
  Eigenvalues<T>::getEigenvalues(
    std::vector<std::complex<valueType> > & eigs
  ) const {
    eigs.clear(); eigs.reserve( this->N );
    for ( int i = 0;i < this->N; ++i )
      eigs.push_back( std::complex<valueType>(Re[i],Im[i]) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  Eigenvectors<T>::Eigenvectors()
  : mem_real("Eigenvectors::mem_real")
  , N(0)
  , Re(nullptr)
  , Im(nullptr)
  , A_saved(nullptr)
  , VL(nullptr)
  , VR(nullptr)
  , Work(nullptr)
  {}

  template <typename T>
  void
  Eigenvectors<T>::allocate( integer Nin ) {
    this->N = Nin;
    // calcolo memoria ottimale
    integer doLwork    = -1;
    T       Lworkdummy = 1;
    integer info = geev(
      true, true,
      this->N,
      nullptr, this->N,
      nullptr, nullptr,
      this->VL, this->N,
      this->VR, this->N,
      &Lworkdummy, doLwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "Eigenvectors::allocate, call geev return info = " << info
    );
    this->Lwork = integer(Lworkdummy);
    this->mem_real.allocate( size_t( this->Lwork + (2+3*this->N) * this->N) );
    this->Re      = this->mem_real( size_t(this->N) );
    this->Im      = this->mem_real( size_t(this->N) );
    this->A_saved = this->mem_real( size_t(this->N*this->N) );
    this->VL      = this->mem_real( size_t(this->N*this->N) );
    this->VR      = this->mem_real( size_t(this->N*this->N) );
    this->Work    = this->mem_real( size_t(this->Lwork) );
  }

  template <typename T>
  void
  Eigenvectors<T>::compute( ) {
    integer info = geev(
      this->VL != nullptr,
      this->VR != nullptr,
      this->N,
      this->A_saved, this->N,
      this->Re,      this->Im,
      this->VL,      this->N,
      this->VR,      this->N,
      this->Work,    this->Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "GeneralizedEigenvectors::compute, call ggevx return info = " << info
    );
  }

  template <typename T>
  Eigenvectors<T>::Eigenvectors(
    integer         NRC,
    valueType const A_data[],
    integer         ldA
  )
  : mem_real("Eigenvectors::mem_real")
  , N(0)
  , Re(nullptr)
  , Im(nullptr)
  , A_saved(nullptr)
  , VL(nullptr)
  , VR(nullptr)
  , Work(nullptr)
  {
    this->setup( NRC, A_data, ldA );
  }

  template <typename T>
  Eigenvectors<T>::Eigenvectors( MatW const & A )
  : mem_real("Eigenvectors::mem_real")
  , N(0)
  , Re(nullptr)
  , Im(nullptr)
  , A_saved(nullptr)
  , VL(nullptr)
  , VR(nullptr)
  , Work(nullptr)
  {
    this->setup( A );
  }

  template <typename T>
  Eigenvectors<T>::Eigenvectors(
    integer         NRC,
    integer         A_nnz,
    valueType const A_values[],
    integer   const A_row[],
    integer   const A_col[]
  )
  : mem_real("Eigenvectors::mem_real")
  , N(0)
  , Re(nullptr)
  , Im(nullptr)
  , A_saved(nullptr)
  , VL(nullptr)
  , VR(nullptr)
  , Work(nullptr)
  {
    this->setup( NRC, A_nnz, A_values, A_row, A_col );
  }

  template <typename T>
  void
  Eigenvectors<T>::setup(
    integer         NRC,
    valueType const A_data[],
    integer         ldA
  ) {
    this->allocate( NRC );
    integer info = gecopy( NRC, NRC, A_data, ldA, A_saved, NRC );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "Eigenvectors::setup, call gecopy return info = " << info
    );
    this->compute();
  }

  template <typename T>
  void
  Eigenvectors<T>::setup( MatW const & A ) {
    this->allocate( A.numRows() );
    integer info = gecopy( this->N, this->N, A.get_data(), A.lDim(),  this->A_saved, this->N );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "Eigenvectors::setup, call gecopy return info = " << info
    );
    this->compute();
  }

  template <typename T>
  void
  Eigenvectors<T>::setup(
    integer         NRC,
    integer         A_nnz,
    valueType const A_values[],
    integer   const A_row[],
    integer   const A_col[]
  ) {
    this->allocate( NRC );
    lapack_wrapper::zero( NRC*NRC, this->A_saved, 1 );
    for ( integer i = 0; i < A_nnz; ++i )
      this->A_saved[A_row[i]+A_col[i]*NRC] += A_values[i];
    this->compute();
  }

  template <typename T>
  void
  Eigenvectors<T>::getEigenvalue(
    integer n, valueType & re, valueType & im
  ) const {
    re = Re[n];
    im = Im[n];
  }

  template <typename T>
  void
  Eigenvectors<T>::getEigenvalue(
    integer n, std::complex<valueType> & eig
  ) const {
    eig = std::complex<valueType>(Re[n],Im[n]);
  }

  template <typename T>
  void
  Eigenvectors<T>::getEigenvalues(
    std::vector<valueType> & re,
    std::vector<valueType> & im
  ) const {
    re.clear(); re.reserve( this->N );
    im.clear(); im.reserve( this->N );
    for ( int i = 0;i < this->N; ++i ) {
      re.push_back( Re[i] );
      im.push_back( Im[i] );
    }
  }

  template <typename T>
  void
  Eigenvectors<T>::getEigenvalues(
    std::vector<std::complex<valueType> > & eigs
  ) const {
    eigs.clear(); eigs.reserve( this->N );
    for ( int i = 0;i < this->N; ++i )
      eigs.push_back( std::complex<valueType>(Re[i],Im[i]) );
  }

  template <typename T>
  void
  Eigenvectors<T>::getLeftEigenvector(
    std::vector<std::vector<complexType> > & vecs
  ) const {
    vecs.resize( size_t(this->N) );
    for ( integer n = 0; n < this->N; ++n ) {
      std::vector<complexType> & v = vecs[n];
      v.clear(); v.reserve( this->N );
      T const * vr = this->VL + n * this->N;
      if ( Im[n] > 0 ) {
        std::vector<complexType> & v1 = vecs[++n];
        v1.clear(); v1.reserve( this->N );
        T const * vi = vr + this->N;
        for ( integer j = 0; j < this->N; ++j ) {
          // salvo vettore "gia" coniugato
          v.push_back( complexType( vr[j], -vi[j] ) );
          v1.push_back( complexType( vr[j], vi[j] ) );
        }
      } else {
        for ( integer j = 0; j < this->N; ++j )
          v.push_back( complexType( vr[j], 0 ) );
      }
    }
  }

  template <typename T>
  void
  Eigenvectors<T>::getRightEigenvector(
    std::vector<std::vector<complexType> > & vecs
  ) const {
    vecs.resize( size_t(this->N) );
    for ( integer n = 0; n < this->N; ++n ) {
      std::vector<complexType> & v = vecs[n];
      v.clear(); v.reserve( this->N );
      T const * vr = this->VR + n * this->N;
      if ( Im[n] > 0 ) {
        std::vector<complexType> & v1 = vecs[++n];
        v1.clear(); v1.reserve( this->N );
        T const * vi = vr + this->N;
        for ( integer j = 0; j < this->N; ++j ) {
          v.push_back( complexType( vr[j], vi[j] ) );
          v1.push_back( complexType( vr[j], -vi[j] ) );
        }
      } else {
        for ( integer j = 0; j < this->N; ++j )
          v.push_back( complexType( vr[j], 0 ) );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  GeneralizedEigenvalues<T>::GeneralizedEigenvalues()
  : mem_real("GeneralizedEigenvalues::mem_real")
  , N(0)
  , alphaRe(nullptr)
  , alphaIm(nullptr)
  , beta(nullptr)
  , Work(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  {}

  template <typename T>
  void
  GeneralizedEigenvalues<T>::allocate( integer Nin ) {
    this->N = Nin;
    // calcolo memoria ottimale
    valueType Lworkdummy;
    integer info = ggev(
      false, false, Nin, nullptr, Nin, nullptr, Nin,
      nullptr, nullptr, nullptr,
      nullptr, Nin, nullptr, Nin, &Lworkdummy, -1
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "GeneralizedEigenvalues::allocate, call geev return info = " << info
    );
    this->Lwork = integer(Lworkdummy);
    this->mem_real.allocate( size_t( this->Lwork + (3+2*this->N) * this->N) );
    this->alphaRe = this->mem_real( size_t(this->N) );
    this->alphaIm = this->mem_real( size_t(this->N) );
    this->beta    = this->mem_real( size_t(this->N) );
    this->Work    = this->mem_real( size_t(this->Lwork) );
    this->A_saved = this->mem_real( size_t(this->N*this->N) );
    this->B_saved = this->mem_real( size_t(this->N*this->N) );
  }

  template <typename T>
  void
  GeneralizedEigenvalues<T>::compute( ) {
    integer info = ggev(
      false, false,
      this->N,
      this->A_saved, this->N,
      this->B_saved, this->N,
      this->alphaRe, this->alphaIm, this->beta,
      nullptr, this->N, nullptr, this->N,
      this->Work, this->Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "GeneralizedEigenvalues::compute, call geev return info = " << info
    );
  }

  template <typename T>
  GeneralizedEigenvalues<T>::GeneralizedEigenvalues(
    integer         NRC,
    valueType const A_data[],
    integer         ldA,
    valueType const B_data[],
    integer         ldB
  )
  : mem_real("GeneralizedEigenvalues::mem_real")
  , N(0)
  , alphaRe(nullptr)
  , alphaIm(nullptr)
  , beta(nullptr)
  , Work(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  {
    this->setup( NRC, A_data, ldA, B_data, ldB );
  }

  template <typename T>
  GeneralizedEigenvalues<T>::GeneralizedEigenvalues(
    MatW const & A, MatW const & B
  )
  : mem_real("GeneralizedEigenvalues::mem_real")
  , N(0)
  , alphaRe(nullptr)
  , alphaIm(nullptr)
  , beta(nullptr)
  , Work(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  {
    this->setup( A, B );
  }

  template <typename T>
  GeneralizedEigenvalues<T>::GeneralizedEigenvalues(
    integer         NRC,
    integer         A_nnz,
    valueType const A_values[],
    integer   const A_row[],
    integer   const A_col[],
    integer         B_nnz,
    valueType const B_values[],
    integer   const B_row[],
    integer   const B_col[]
  )
  : mem_real("GeneralizedEigenvalues::mem_real")
  , N(0)
  , alphaRe(nullptr)
  , alphaIm(nullptr)
  , beta(nullptr)
  , Work(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  {
    this->setup(
      NRC,
      A_nnz, A_values, A_row, A_col,
      B_nnz, B_values, B_row, B_col
    );
  }

  template <typename T>
  void
  GeneralizedEigenvalues<T>::setup(
    integer         NRC,
    valueType const A_data[],
    integer         ldA,
    valueType const B_data[],
    integer         ldB
  ) {
    this->allocate( NRC );
    integer info1 = gecopy( NRC, NRC, A_data, ldA, A_saved, NRC );
    integer info2 = gecopy( NRC, NRC, B_data, ldB, B_saved, NRC );
    LAPACK_WRAPPER_ASSERT(
      info1 == 0 && info2 == 0,
      "GeneralizedEigenvalues::setup, call gecopy return info1 = " << info1 <<
      ", info2 = " << info2
    );
    this->compute();
  }

  template <typename T>
  void
  GeneralizedEigenvalues<T>::setup( MatW const & A, MatW const & B ) {
    this->allocate( A.numRows() );
    integer info1 = gecopy( this->N, this->N, A.get_data(), A.lDim(),  this->A_saved, this->N );
    integer info2 = gecopy( this->N, this->N, B.get_data(), B.lDim(),  this->B_saved, this->N );
    LAPACK_WRAPPER_ASSERT(
      info1 == 0 && info2 == 0,
      "GeneralizedEigenvalues::setup, call gecopy return info1 = " << info1 <<
      ", info2 = " << info2
    );
    this->compute();
  }

  template <typename T>
  void
  GeneralizedEigenvalues<T>::setup(
    integer         NRC,
    integer         A_nnz,
    valueType const A_values[],
    integer   const A_row[],
    integer   const A_col[],
    integer         B_nnz,
    valueType const B_values[],
    integer   const B_row[],
    integer   const B_col[]
  ) {
    this->allocate( NRC );
    lapack_wrapper::zero( NRC*NRC, this->A_saved, 1 );
    lapack_wrapper::zero( NRC*NRC, this->B_saved, 1 );
    for ( integer i = 0; i < A_nnz; ++i )
      this->A_saved[A_row[i]+A_col[i]*NRC] += A_values[i];
    for ( integer i = 0; i < B_nnz; ++i )
      this->B_saved[B_row[i]+B_col[i]*NRC] += B_values[i];
    this->compute();
  }

  template <typename T>
  void
  GeneralizedEigenvalues<T>::getEigenvalue(
    integer n, valueType & re, valueType & im
  ) const {
    re = alphaRe[n]/beta[n];
    im = alphaIm[n]/beta[n];
  }

  template <typename T>
  void
  GeneralizedEigenvalues<T>::getEigenvalue(
    integer n, std::complex<valueType> & eig
  ) const {
    eig = std::complex<valueType>(alphaRe[n],alphaIm[n])/beta[n];
  }

  template <typename T>
  void
  GeneralizedEigenvalues<T>::getEigenvalues(
    std::vector<valueType> & re,
    std::vector<valueType> & im
  ) const {
    re.clear(); re.reserve( this->N );
    im.clear(); im.reserve( this->N );
    for ( int i = 0;i < this->N; ++i ) {
      re.push_back( alphaRe[i]/beta[i] );
      im.push_back( alphaIm[i]/beta[i] );
    }
  }

  template <typename T>
  void
  GeneralizedEigenvalues<T>::getEigenvalues(
    std::vector<std::complex<valueType> > & eigs
  ) const {
    eigs.clear(); eigs.reserve( this->N );
    for ( int i = 0;i < this->N; ++i )
      eigs.push_back( std::complex<valueType>(alphaRe[i],alphaIm[i])/beta[i] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  GeneralizedEigenvectors<T>::GeneralizedEigenvectors()
  : mem_real("GeneralizedEigenvectors::mem_real")
  , mem_int("GeneralizedEigenvectors::mem_int")
  , N(0)
  , alphaRe(nullptr)
  , alphaIm(nullptr)
  , beta(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  , VL(nullptr)
  , VR(nullptr)
  , lscale(nullptr)
  , rscale(nullptr)
  , rconde(nullptr)
  , rcondv(nullptr)
  , Work(nullptr)
  , iWork(nullptr)
  , bWork(nullptr)
  {}

  template <typename T>
  void
  GeneralizedEigenvectors<T>::allocate( integer Nin ) {
    this->N = Nin;
    // calcolo memoria ottimale
    integer doLwork    = -1;
    T       Lworkdummy = 1;
    integer info = ggevx(
      lapack_wrapper::PERMUTE_AND_SCALE,
      false, false,
      lapack_wrapper::EIGENVALUES_AND_EIGENVECTORS,
      this->N, nullptr, this->N, nullptr, this->N,
      nullptr, nullptr, nullptr, this->VL, this->N, this->VR, this->N,
      ilo, ihi,
      nullptr, nullptr,
      abnorm,  bbnorm,
      nullptr, nullptr,
      &Lworkdummy, doLwork,
      nullptr, nullptr
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "GeneralizedEigenvectors::allocate, call geev return info = " << info
    );
    this->Lwork = integer(Lworkdummy);
    this->mem_real.allocate( size_t( this->Lwork + (7+4*this->N) * this->N) );
    this->mem_int.allocate( size_t( 2*this->N + 6 ) );
    this->alphaRe = this->mem_real( size_t(this->N) );
    this->alphaIm = this->mem_real( size_t(this->N) );
    this->beta    = this->mem_real( size_t(this->N) );
    this->A_saved = this->mem_real( size_t(this->N*this->N) );
    this->B_saved = this->mem_real( size_t(this->N*this->N) );
    this->VL      = this->mem_real( size_t(this->N*this->N) );
    this->VR      = this->mem_real( size_t(this->N*this->N) );
    this->lscale  = this->mem_real( size_t(this->N) );
    this->rscale  = this->mem_real( size_t(this->N) );
    this->rconde  = this->mem_real( size_t(this->N) );
    this->rcondv  = this->mem_real( size_t(this->N) );
    this->Work    = this->mem_real( size_t(this->Lwork) );
    this->iWork   = this->mem_int( size_t(this->N+6) );
    this->bWork   = this->mem_int( size_t(this->N) );
  }

  template <typename T>
  void
  GeneralizedEigenvectors<T>::compute( ) {
    integer info = ggevx(
      lapack_wrapper::PERMUTE_ONLY, // "B", // ATTENZIONE LA SCALATURA NON FUNZIONA
      this->VL != nullptr,
      this->VR != nullptr,
      lapack_wrapper::EIGENVALUES_AND_EIGENVECTORS,
      this->N,
      this->A_saved, this->N,
      this->B_saved, this->N,
      this->alphaRe, this->alphaIm, this->beta,
      this->VL, this->N,
      this->VR, this->N,
      ilo, ihi,
      this->lscale, this->rscale,
      abnorm,       bbnorm,
      this->rconde, this->rcondv,
      this->Work,   this->Lwork,
      this->iWork,  this->bWork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "GeneralizedEigenvectors::compute, call ggevx return info = " << info
    );
  }

  template <typename T>
  GeneralizedEigenvectors<T>::GeneralizedEigenvectors(
    integer         NRC,
    valueType const A_data[],
    integer         ldA,
    valueType const B_data[],
    integer         ldB
  )
  : mem_real("GeneralizedEigenvectors::mem_real")
  , mem_int("GeneralizedEigenvectors::mem_int")
  , N(0)
  , alphaRe(nullptr)
  , alphaIm(nullptr)
  , beta(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  , VL(nullptr)
  , VR(nullptr)
  , lscale(nullptr)
  , rscale(nullptr)
  , rconde(nullptr)
  , rcondv(nullptr)
  , Work(nullptr)
  , iWork(nullptr)
  , bWork(nullptr)
  {
    this->setup( NRC, A_data, ldA, B_data, ldB );
  }

  template <typename T>
  GeneralizedEigenvectors<T>::GeneralizedEigenvectors(
    MatW const & A, MatW const & B
  )
  : mem_real("GeneralizedEigenvectors::mem_real")
  , mem_int("GeneralizedEigenvectors::mem_int")
  , N(0)
  , alphaRe(nullptr)
  , alphaIm(nullptr)
  , beta(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  , VL(nullptr)
  , VR(nullptr)
  , lscale(nullptr)
  , rscale(nullptr)
  , rconde(nullptr)
  , rcondv(nullptr)
  , Work(nullptr)
  , iWork(nullptr)
  , bWork(nullptr)
  {
    this->setup( A, B );
  }

  template <typename T>
  GeneralizedEigenvectors<T>::GeneralizedEigenvectors(
    integer         NRC,
    integer         A_nnz,
    valueType const A_values[],
    integer   const A_row[],
    integer   const A_col[],
    integer         B_nnz,
    valueType const B_values[],
    integer   const B_row[],
    integer   const B_col[]
  )
  : mem_real("GeneralizedEigenvectors::mem_real")
  , mem_int("GeneralizedEigenvectors::mem_int")
  , N(0)
  , alphaRe(nullptr)
  , alphaIm(nullptr)
  , beta(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  , VL(nullptr)
  , VR(nullptr)
  , lscale(nullptr)
  , rscale(nullptr)
  , rconde(nullptr)
  , rcondv(nullptr)
  , Work(nullptr)
  , iWork(nullptr)
  , bWork(nullptr)
  {
    this->setup(
      NRC,
      A_nnz, A_values, A_row, A_col,
      B_nnz, B_values, B_row, B_col
    );
  }

  template <typename T>
  void
  GeneralizedEigenvectors<T>::setup(
    integer         NRC,
    valueType const A_data[],
    integer         ldA,
    valueType const B_data[],
    integer         ldB
  ) {
    this->allocate( NRC );
    integer info1 = gecopy( NRC, NRC, A_data, ldA, A_saved, NRC );
    integer info2 = gecopy( NRC, NRC, B_data, ldB, B_saved, NRC );
    LAPACK_WRAPPER_ASSERT(
      info1 == 0 && info2 == 0,
      "GeneralizedEigenvectors::setup, call gecopy return info1 = " << info1 <<
      ", info2 = " << info2
    );
    this->compute();
  }

  template <typename T>
  void
  GeneralizedEigenvectors<T>::setup( MatW const & A, MatW const & B ) {
    this->allocate( A.numRows() );
    integer info1 = gecopy( this->N, this->N, A.get_data(), A.lDim(),  this->A_saved, this->N );
    integer info2 = gecopy( this->N, this->N, B.get_data(), B.lDim(),  this->B_saved, this->N );
    LAPACK_WRAPPER_ASSERT(
      info1 == 0 && info2 == 0,
      "GeneralizedEigenvectors::setup, call gecopy return info1 = " << info1 <<
      ", info2 = " << info2
    );
    this->compute();
  }

  template <typename T>
  void
  GeneralizedEigenvectors<T>::setup(
    integer         NRC,
    integer         A_nnz,
    valueType const A_values[],
    integer   const A_row[],
    integer   const A_col[],
    integer         B_nnz,
    valueType const B_values[],
    integer   const B_row[],
    integer   const B_col[]
  ) {
    this->allocate( NRC );
    lapack_wrapper::zero( NRC*NRC, this->A_saved, 1 );
    lapack_wrapper::zero( NRC*NRC, this->B_saved, 1 );
    for ( integer i = 0; i < A_nnz; ++i )
      this->A_saved[A_row[i]+A_col[i]*NRC] += A_values[i];
    for ( integer i = 0; i < B_nnz; ++i )
      this->B_saved[B_row[i]+B_col[i]*NRC] += B_values[i];
    this->compute();
  }

  template <typename T>
  void
  GeneralizedEigenvectors<T>::getEigenvalue(
    integer n, valueType & re, valueType & im
  ) const {
    re = alphaRe[n]/beta[n];
    im = alphaIm[n]/beta[n];
  }

  template <typename T>
  void
  GeneralizedEigenvectors<T>::getEigenvalue(
    integer n, std::complex<valueType> & eig
  ) const {
    eig = std::complex<valueType>(alphaRe[n],alphaIm[n])/beta[n];
  }

  template <typename T>
  void
  GeneralizedEigenvectors<T>::getEigenvalues(
    std::vector<valueType> & re,
    std::vector<valueType> & im
  ) const {
    re.clear(); re.reserve( this->N );
    im.clear(); im.reserve( this->N );
    for ( int i = 0;i < this->N; ++i ) {
      re.push_back( alphaRe[i]/beta[i] );
      im.push_back( alphaIm[i]/beta[i] );
    }
  }

  template <typename T>
  void
  GeneralizedEigenvectors<T>::getEigenvalues(
    std::vector<std::complex<valueType> > & eigs
  ) const {
    eigs.clear(); eigs.reserve( this->N );
    for ( int i = 0;i < this->N; ++i )
      eigs.push_back( std::complex<valueType>(alphaRe[i],alphaIm[i])/beta[i] );
  }

  template <typename T>
  void
  GeneralizedEigenvectors<T>::getLeftEigenvector(
    std::vector<std::vector<complexType> > & vecs
  ) const {
    vecs.resize( size_t(this->N) );
    for ( integer n = 0; n < this->N; ++n ) {
      std::vector<complexType> & v = vecs[n];
      v.clear(); v.reserve( this->N );
      T const * vr = this->VL + n * this->N;
      if ( alphaIm[n] > 0 ) {
        std::vector<complexType> & v1 = vecs[++n];
        v1.clear(); v1.reserve( this->N );
        T const * vi = vr + this->N;
        for ( integer j = 0; j < this->N; ++j ) {
          // salvo vettore "gia" coniugato
          v.push_back( complexType( vr[j], -vi[j] ) );
          v1.push_back( complexType( vr[j], vi[j] ) );
        }
      } else {
        for ( integer j = 0; j < this->N; ++j )
          v.push_back( complexType( vr[j], 0 ) );
      }
    }
  }

  template <typename T>
  void
  GeneralizedEigenvectors<T>::getRightEigenvector(
    std::vector<std::vector<complexType> > & vecs
  ) const {
    vecs.resize( size_t(this->N) );
    for ( integer n = 0; n < this->N; ++n ) {
      std::vector<complexType> & v = vecs[n];
      v.clear(); v.reserve( this->N );
      T const * vr = this->VR + n * this->N;
      if ( alphaIm[n] > 0 ) {
        std::vector<complexType> & v1 = vecs[++n];
        v1.clear(); v1.reserve( this->N );
        T const * vi = vr + this->N;
        for ( integer j = 0; j < this->N; ++j ) {
          v.push_back( complexType( vr[j], vi[j] ) );
          v1.push_back( complexType( vr[j], -vi[j] ) );
        }
      } else {
        for ( integer j = 0; j < this->N; ++j )
          v.push_back( complexType( vr[j], 0 ) );
      }
    }
  }

}

///
/// eof: eig.cxx
///
