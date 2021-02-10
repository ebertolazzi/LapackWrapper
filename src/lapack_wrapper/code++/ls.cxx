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
/// file: ls.cxx
///

namespace lapack_wrapper {

  /*\
   |   _     ____ ____
   |  | |   / ___/ ___|
   |  | |   \___ \___ \
   |  | |___ ___) |__) |
   |  |_____|____/____/
   |
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSS_no_alloc<T>::no_allocate(
    integer     NR,
    integer     NC,
    integer     Lwork,
    valueType * Work
  ) {
    integer minRC  = std::min(NR,NC);
    integer NRC    = NR*NC;
    integer Lwmin  = 2*NRC+minRC;
    UTILS_ASSERT(
      Lwork >= 2*NRC,
      "LSS_no_alloc( NR = {}, NC = {}, LWork = {}, ...)\n"
      "LWork must be >= {}\n",
      NR, NC, Lwork, Lwmin
    );
    m_nrows    = NR;
    m_ncols    = NC;
    m_Amat     = Work;
    m_AmatWork = Work+NRC;
    m_sigma    = Work+2*NRC;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  integer
  LSS_no_alloc<T>::getL( integer NR, integer NC, integer nrhs ) const {
    valueType tmp;
    integer NRC = NR > NC ? NR : NC;
    integer info = gelss(
      NR, NC, nrhs,
      nullptr, NR, nullptr, NRC, nullptr,
      m_rcond, m_rank, &tmp, -1
    );
    UTILS_ASSERT(
      info == 0,
      "LSS_no_alloc::getL( NR={}, NC={}, nrhs={} )\n"
      "call of gelss return info = {}\n",
      NR, NC, nrhs, info
    );
    return integer(tmp);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSS_no_alloc<T>::factorize_nodim(
    char const      who[],
    valueType const A[],
    integer         LDA
  ) {
    integer info = gecopy(
      m_nrows, m_ncols, A, LDA, m_Amat, m_nrows
    );
    UTILS_ASSERT(
      info == 0,
      "LSS_no_alloc::factorize( who={},...)"
      " call lapack_wrapper::gecopy return info = {}\n",
      who, info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LSS_no_alloc<T>::factorize_nodim( valueType const A[], integer LDA ) {
    integer info = gecopy(
      m_nrows, m_ncols, A, LDA, m_Amat, m_nrows
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LSS_no_alloc<T>::solve( valueType xb[] ) const {
    // check memory
    integer L = getL( m_nrows, m_ncols, 1 );
    if ( L > m_Lwork ) {
      m_Lwork = L;
      m_Work  = m_allocWork.malloc( size_t(L) );
    }
    // save matrix
    copy( m_nrows*m_ncols, m_Amat, 1, m_AmatWork, 1);
    integer nrc = m_nrows > m_ncols ? m_nrows : m_ncols;
    integer info = gelss(
      m_nrows, m_ncols, 1,
      m_AmatWork, m_nrows,
      xb, nrc,
      m_sigma, m_rcond, m_rank,
      m_Work, m_Lwork
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LSS_no_alloc<T>::t_solve( valueType xb[] ) const {
    // check memory
    integer L = getL( m_ncols, m_nrows, 1 );
    if ( L > m_Lwork ) {
      m_Lwork = L;
      m_Work  = m_allocWork.malloc( size_t(L) );
    }
    // save matrix
    for ( integer i = 0; i < m_ncols; ++i )
      copy( m_nrows, m_Amat+i*m_nrows, 1, m_AmatWork+i, m_ncols );
    integer nrc = m_nrows > m_ncols ? m_nrows : m_ncols;
    integer info = gelss(
      m_ncols, m_nrows, 1,
      m_AmatWork, m_ncols,
      xb, nrc,
      m_sigma, m_rcond, m_rank,
      m_Work, m_Lwork
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LSS_no_alloc<T>::solve(
    integer   nrhs,
    valueType B[],
    integer   ldB
  ) const {
    // check memory
    integer L = getL( m_nrows, m_ncols, nrhs );
    if ( L > m_Lwork ) {
      m_Lwork = L;
      m_Work  = m_allocWork.malloc( size_t(L) );
    }
    // save matrix
    copy( m_nrows*m_ncols, m_Amat, 1, m_AmatWork, 1 );
    integer info = gelss(
      m_nrows, m_ncols, nrhs,
      m_AmatWork, m_nrows,
      B, ldB,
      m_sigma, m_rcond, m_rank,
      m_Work, m_Lwork
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LSS_no_alloc<T>::t_solve(
    integer   nrhs,
    valueType B[],
    integer   ldB
  ) const {
    // check memory
    integer L = getL( m_ncols, m_nrows, nrhs );
    if ( L > m_Lwork ) {
      m_Lwork = L;
      m_Work  = m_allocWork.malloc( size_t(L) );
    }
    // save matrix
    for ( integer i = 0; i < m_ncols; ++i )
      copy( m_nrows, m_Amat+i*m_nrows, 1, m_AmatWork+i, m_ncols );

    integer info = gelss(
      m_ncols, m_nrows, nrhs,
      m_AmatWork, m_ncols,
      B, ldB,
      m_sigma, m_rcond, m_rank,
      m_Work, m_Lwork
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSS<T>::allocate( integer NR, integer NC ) {
    if ( m_nrows != NR || m_ncols != NC ) {
      integer Lwork = 2*NR*NC+std::min(NR,NC);
      this->no_allocate(
        NR, NC, Lwork, m_allocReals.malloc( size_t(Lwork) )
      );
    }
  }

  /*\
   |  _     ______   __
   | | |   / ___\ \ / /
   | | |   \___ \\ V /
   | | |___ ___) || |
   | |_____|____/ |_|
   |
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY_no_alloc<T>::no_allocate(
    integer     NR,
    integer     NC,
    integer     Lwork,
    valueType * Work,
    integer     Liwork,
    integer   * iWork
  ) {
    integer NRC   = NR*NC;
    integer Lwmin = 2*NRC;
    UTILS_ASSERT(
      Lwork >= 2*NRC && Liwork >= NC,
      "LSS_no_alloc( NR = {}, NC = {}, LWork = {}, ..., Liwork = {})\n"
      "LWork must be >= {} and Liwork must be >= {}\n",
      NR, NC, Lwork, Liwork, Lwmin, NC
    );
    m_nrows    = NR;
    m_ncols    = NC;
    m_Amat     = Work;
    m_AmatWork = Work+NRC;
    m_jpvt     = iWork;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  integer
  LSY_no_alloc<T>::getL( integer NR, integer NC, integer nrhs ) const {
    valueType tmp;
    integer NRC = NR > NC ? NR : NC;
    integer info = gelsy(
      NR, NC, nrhs,
      nullptr, NR, nullptr, NRC, nullptr,
      m_rcond, m_rank, &tmp, -1
    );
    UTILS_ASSERT(
      info == 0,
      "LSY_no_alloc::getL( NR={}, NC={}, nrhs={} )\n"
      "call of gelss return info = {}\n",
      NR, NC, nrhs, info
    );
    return integer(tmp);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY_no_alloc<T>::factorize_nodim(
    char const      who[],
    valueType const A[],
    integer         LDA
  ) {
    integer info = gecopy(
      m_nrows, m_ncols, A, LDA, m_Amat, m_nrows
    );
    UTILS_ASSERT(
      info == 0,
      "LSY::factorize[{}] call lapack_wrapper::gecopy return info = {}\n",
      who, info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LSY_no_alloc<T>::factorize_nodim( valueType const A[], integer LDA ) {
    integer info = gecopy(
      m_nrows, m_ncols, A, LDA, m_Amat, m_nrows
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LSY_no_alloc<T>::solve( valueType xb[] ) const {
    // check memory
    integer L = getL( m_nrows, m_ncols, 1 );
    if ( L > m_Lwork ) {
      m_Lwork = L;
      m_Work  = m_allocWork.malloc( size_t(L) );
    }
    // save matrix
    copy( m_nrows*m_ncols, m_Amat, 1, m_AmatWork, 1 );
    integer nrc = m_nrows > m_ncols ? m_nrows : m_ncols;
    integer info = gelsy(
      m_nrows, m_ncols, 1,
      m_AmatWork, m_nrows,
      xb, nrc, m_jpvt,
      m_rcond, m_rank,
      m_Work, m_Lwork
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LSY_no_alloc<T>::t_solve( valueType xb[] ) const {
    // check memory
    integer L = getL( m_ncols, m_nrows, 1 );
    if ( L > m_Lwork ) {
      m_Lwork = L;
      m_Work  = m_allocWork.malloc( size_t( L ) );
    }
    // save matrix
    for ( integer i = 0; i < m_ncols; ++i )
      copy( m_nrows, m_Amat+i*m_nrows, 1, m_AmatWork+i, m_ncols );
    integer nrc = m_nrows > m_ncols ? m_nrows : m_ncols;
    integer info = gelsy(
      m_ncols, m_nrows, 1,
      m_AmatWork, m_ncols,
      xb, nrc, m_jpvt,
      m_rcond, m_rank,
      m_Work, m_Lwork
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LSY_no_alloc<T>::solve(
    integer   nrhs,
    valueType B[],
    integer   ldB
  ) const {
    // check memory
    integer L = getL( m_nrows, m_ncols, nrhs );
    if ( L > m_Lwork ) {
      m_Lwork = L;
      m_Work  = m_allocWork.malloc( size_t(L) );
    }
    // save matrix
    copy( m_nrows*m_ncols, m_Amat, 1, m_AmatWork, 1 );
    integer info = gelsy(
      m_nrows, m_ncols, nrhs,
      m_AmatWork, m_nrows,
      B, ldB, m_jpvt,
      m_rcond, m_rank,
      m_Work, m_Lwork
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LSY_no_alloc<T>::t_solve(
    integer   nrhs,
    valueType B[],
    integer   ldB
  ) const {
    // check memory
    integer L = getL( m_ncols, m_nrows, nrhs );
    if ( L > m_Lwork ) {
      m_Lwork = L;
      m_Work  = m_allocWork.malloc( size_t(L) );
    }
    // save matrix
    for ( integer i = 0; i < m_ncols; ++i )
      copy( m_nrows, m_Amat+i*m_nrows, 1, m_AmatWork+i, m_ncols );
    integer info = gelsy(
      m_ncols, m_nrows, nrhs,
      m_AmatWork, m_ncols,
      B, ldB, m_jpvt,
      m_rcond, m_rank,
      m_Work, m_Lwork
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY<T>::allocate( integer NR, integer NC ) {
    if ( m_nrows != NR || m_ncols != NC ) {
      integer Lwork = 2*NR*NC;
      this->no_allocate(
        NR, NC,
        Lwork, m_allocReals.malloc( size_t(Lwork) ),
        NC,    m_allocInts.malloc( size_t(NC) )
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY<T>::factorize(
    char const      who[],
    integer         NR,
    integer         NC,
    valueType const A[],
    integer         LDA
  ) {
    allocate( NR, NC );
    integer info = gecopy(
      m_nrows, m_ncols, A, LDA, m_Amat, m_nrows
    );
    UTILS_ASSERT(
      info == 0,
      "LSY::factorize[{}] call lapack_wrapper::gecopy return info = {}\n",
      who, info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LSY<T>::factorize(
    integer         NR,
    integer         NC,
    valueType const A[],
    integer         LDA
  ) {
    allocate( NR, NC );
    integer info = gecopy(
      m_nrows, m_ncols, A, LDA, m_Amat, m_nrows
    );
    return info == 0;
  }

}

///
/// eof: ls.cxx
///
