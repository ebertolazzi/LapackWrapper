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

  template <typename T>
  void
  SVD<T>::allocate( integer NR, integer NC ) {

    if ( nRow != NR || nCol != NC ) {
      nRow  = NR;
      nCol  = NC;
      minRC = std::min(NR,NC);
      valueType tmp;
      integer info = gesvd(
        REDUCED, REDUCED,
        NR, NC,
        nullptr, NR,
        nullptr,
        nullptr, NR,
        nullptr, minRC,
        &tmp, -1
      );
      LAPACK_WRAPPER_ASSERT(
        info == 0, "SVD::allocate, in gesvd info = " << info
      );
      Lwork = integer(tmp);
      info = gesdd(
        REDUCED,
        NR, NC,
        nullptr, NR,
        nullptr,
        nullptr, NR,
        nullptr, minRC,
        &tmp, -1, nullptr
      );
      if ( integer(tmp) > Lwork ) Lwork = integer(tmp);
      allocReals.allocate( size_t(nRow*nCol+minRC*(nRow+nCol+1)+Lwork) );
      Amat = allocReals( size_t(nRow*nCol));
      Svec  = allocReals( size_t(minRC) );
      Umat  = allocReals( size_t(minRC*nRow) );
      VTmat = allocReals( size_t(minRC*nCol) );
      Work  = allocReals( size_t(Lwork) );
      allocIntegers.allocate( size_t(8*minRC) );
      IWork = allocIntegers( size_t(8*minRC) );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  SVD<T>::factorize( char const who[] ) {
    integer info;
    switch ( svd_used ) {
    case USE_GESVD:
      info = gesvd(
        REDUCED,
        REDUCED,
        nRow, nCol, Amat, nRow,
        Svec,
        Umat, nRow,
        VTmat, minRC,
        Work, Lwork
      );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "SVD::factorize[" << who <<
        "] call lapack_wrapper::gesvd return info = " << info
      );
      break;
    case USE_GESDD:
      info = gesdd(
        REDUCED,
        nRow, nCol, Amat, nRow,
        Svec,
        Umat, nRow,
        VTmat, minRC,
        Work, Lwork, IWork
      );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "SVD::factorize[" << who <<
        "] call lapack_wrapper::gesdd return info = " << info
      );
      break;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  SVD<T>::solve( valueType xb[] ) const {
    // A = U*S*VT
    // U*S*VT*x=b --> VT^T S^+ U^T b
    // U  nRow x minRC
    // VT minRC x nCol
    valueType smin = rcond*Svec[0];
    Ut_mul( 1.0, xb, 1, 0.0, Work, 1 );
    for ( integer i = 0; i < minRC; ++i ) Work[i] /= std::max(Svec[i],smin);
    V_mul( 1.0, Work, 1, 0.0, xb, 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  SVD<T>::t_solve( valueType xb[] ) const {
    // A = U*S*VT
    // U*S*VT*x=b --> VT^T S^+ U^T b
    // U  nRow x minRC
    // VT minRC x nCol
    valueType smin = rcond*Svec[0];
    Vt_mul( 1.0, xb, 1, 0.0, Work, 1 );
    for ( integer i = 0; i < minRC; ++i ) Work[i] /= std::max(Svec[i],smin);
    U_mul( 1.0, Work, 1, 0.0, xb, 1 );
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

  template <typename T>
  GeneralizedSVD<T>::GeneralizedSVD()
  : mem_real("GeneralizedSVD(real)")
  , mem_int("GeneralizedSVD(int)")
  , M(0)
  , N(0)
  , P(0)
  , K(0)
  , L(0)
  , Lwork(0)
  , Work(nullptr)
  , IWork(nullptr)
  , alpha_saved(nullptr)
  , beta_saved(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  , U_saved(nullptr)
  , V_saved(nullptr)
  , Q_saved(nullptr)
  {}

  template <typename T>
  GeneralizedSVD<T>::GeneralizedSVD(
    integer         m,
    integer         n,
    integer         p,
    valueType const A[], integer ldA_in,
    valueType const B[], integer ldB_in
  )
  : mem_real("GeneralizedSVD(real)")
  , mem_int("GeneralizedSVD(int)")
  , M(0)
  , N(0)
  , P(0)
  , K(0)
  , L(0)
  , Lwork(0)
  , Work(nullptr)
  , IWork(nullptr)
  , alpha_saved(nullptr)
  , beta_saved(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  , U_saved(nullptr)
  , V_saved(nullptr)
  , Q_saved(nullptr)
  { this->setup( m, n, p, A, ldA_in, B, ldB_in ); }

  template <typename T>
  GeneralizedSVD<T>::GeneralizedSVD( MatW const & A, MatW const & B )
  : mem_real("GeneralizedSVD(real)")
  , mem_int("GeneralizedSVD(int)")
  , M(0)
  , N(0)
  , P(0)
  , K(0)
  , L(0)
  , Lwork(0)
  , Work(nullptr)
  , IWork(nullptr)
  , alpha_saved(nullptr)
  , beta_saved(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  , U_saved(nullptr)
  , V_saved(nullptr)
  , Q_saved(nullptr)
  { this->setup( A, B ); }

  template <typename T>
  GeneralizedSVD<T>::GeneralizedSVD(
    integer         m,
    integer         n,
    integer         p,
    integer         A_nnz,
    valueType const A_values[],
    integer   const A_row[],
    integer   const A_col[],
    integer         B_nnz,
    valueType const B_values[],
    integer   const B_row[],
    integer   const B_col[]
  )
  : mem_real("GeneralizedSVD(real)")
  , mem_int("GeneralizedSVD(int)")
  , M(0)
  , N(0)
  , P(0)
  , K(0)
  , L(0)
  , Lwork(0)
  , Work(nullptr)
  , IWork(nullptr)
  , alpha_saved(nullptr)
  , beta_saved(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  , U_saved(nullptr)
  , V_saved(nullptr)
  , Q_saved(nullptr)
  { this->setup(
      m, n, p,
      A_nnz, A_values, A_row, A_col,
      B_nnz, B_values, B_row, B_col
    );
  }

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
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "GeneralizedSVD<T>::allocate(m=" << m << ",n" << n << ",p=" << p <<
      ") failed, info = " << info
    );
    this->M     = m;
    this->N     = n;
    this->P     = p;
    this->Lwork = integer(wL);

    this->mem_int.allocate( n );
    this->IWork = mem_int( n );

    this->mem_real.allocate( Lwork + (m+p+2)*n + m*m + p*p + n*n );
    this->Work        = mem_real( Lwork );
    this->alpha_saved = mem_real( n );
    this->beta_saved  = mem_real( n );
    this->A_saved     = mem_real( m*n );
    this->B_saved     = mem_real( p*n );
    this->U_saved     = mem_real( m*m );
    this->V_saved     = mem_real( p*p );
    this->Q_saved     = mem_real( n*n );

    this->U.setup( this->U_saved, m, m, m );
    this->V.setup( this->V_saved, p, p, p );
    this->Q.setup( this->Q_saved, n, n, n );
    this->Dbeta.setup( this->beta_saved, n );
    this->Dalpha.setup( this->alpha_saved, n );
  }

  template <typename T>
  void
  GeneralizedSVD<T>::compute() {
    integer info = ggsvd(
      true, true, true,
      this->M, this->N, this->P, this->K, this->L,
      this->A_saved, this->M,
      this->B_saved, this->P,
      this->alpha_saved,
      this->beta_saved,
      this->U_saved, this->M,
      this->V_saved, this->P,
      this->Q_saved, this->N,
      Work,
      Lwork,
      IWork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "GeneralizedSVD<T>::compute() failed, info = " << info
    );
    // R is stored in A(1:K+L,N-K-L+1:N) on exit.
    this->R.setup( this->A_saved + this->M * ( this->N-this->K-this->L ),
                   this->N, this->K+this->L, this->M );
  }

  template <typename T>
  void
  GeneralizedSVD<T>::setup(
    integer         m,
    integer         n,
    integer         p,
    valueType const A[], integer ldA,
    valueType const B[], integer ldB
  ) {
    this->allocate( m, n, p );
    integer info = gecopy( m, n, A, ldA, A_saved, m );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "GeneralizedSVD<T>::setup(...) failed to copy A, info = " << info
    );
    info = gecopy( p, n, B, ldB, B_saved, p );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "GeneralizedSVD<T>::setup(...) failed to copy B, info = " << info
    );
    compute();
  }

  template <typename T>
  void
  GeneralizedSVD<T>::setup( MatW const & A, MatW const & B ) {
    integer m = A.numRows();
    integer n = A.numCols();
    integer p = B.numRows();
    LAPACK_WRAPPER_ASSERT(
      n == B.numCols(),
      "GeneralizedSVD<T>::setup( A, B ) incompatible matrices\n" <<
      "A is " << A.numRows() << " x " << A.numCols() << "\n"
      "B is " << B.numRows() << " x " << B.numCols()
    );
    this->setup( m, n, p, A.get_data(), A.lDim(), B.get_data(), B.lDim() );
  }

  template <typename T>
  void
  GeneralizedSVD<T>::setup(
    integer         m,
    integer         n,
    integer         p,
    integer         A_nnz,
    valueType const A_values[],
    integer   const A_row[],
    integer   const A_col[],
    integer         B_nnz,
    valueType const B_values[],
    integer   const B_row[],
    integer   const B_col[]
  ) {
    this->allocate( m, n, p );
    lapack_wrapper::zero( this->N*this->M, this->A_saved, 1 );
    lapack_wrapper::zero( this->P*this->M, this->B_saved, 1 );
    for ( integer k = 0; k < A_nnz; ++k ) A(A_row[k],A_col[k]) = A_values[k];
    for ( integer k = 0; k < B_nnz; ++k ) B(B_row[k],B_col[k]) = B_values[k];
    compute();
  }

  template <typename T>
  void
  GeneralizedSVD<T>::info( ostream_type & stream, valueType eps ) const {
    stream
      << "A = " << this->M << " x " << this->N << '\n'
      << "B = " << this->P << " x " << this->N << '\n';
    for ( integer i = 0; i < this->N; ++i ) {
      T a =  this->alpha_saved[i];
      T b =  this->beta_saved[i];
      stream
        << "alpha[" << i << "]=" << std::setw(14) << a
        << ", beta[" << i << "]=" << std::setw(14) << b
        << ", alpha^2+beta^2 = " << a*a+b*b << '\n';
    }
    stream << "U\n"; this->U.print0( stream, eps );
    stream << "V\n"; this->V.print0( stream, eps );
    stream << "Q\n"; this->Q.print0( stream, eps );
    stream << "R\n"; this->R.print0( stream, eps );
    stream << "C\n"; this->Dalpha.print( stream );
    stream << "S\n"; this->Dbeta.print( stream );
    stream << '\n';
  }

}

///
/// eof: svd.cxx
///
