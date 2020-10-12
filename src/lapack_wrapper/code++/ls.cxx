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
    LW_ASSERT(
      Lwork >= 2*NRC,
      "LSS_no_alloc( NR = {}, NC = {}, LWork = {}, ...)\n"
      "LWork must be >= {}\n",
      NR, NC, Lwork, Lwmin
    );
    m_nRows    = NR;
    m_nCols    = NC;
    m_Amat     = Work;
    m_AmatWork = Work+NRC;
    m_sigma    = Work+2*NRC;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  integer
  LSS_no_alloc<T>::getL( integer NR, integer NC, integer nrhs ) const {
    valueType tmp;
    integer info = gelss(
      NR, NC, nrhs,
      nullptr, NR, nullptr, NR, nullptr,
      m_rcond, m_rank, &tmp, -1
    );
    LW_ASSERT( info == 0, "LSS_no_alloc::getL, in gelss info = {}\n", info );
    return integer(tmp);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSS_no_alloc<T>::factorize(
    char const      who[],
    valueType const A[],
    integer         LDA
  ) {
    integer info = gecopy(
      m_nRows, m_nCols, A, LDA, m_Amat, m_nRows
    );
    LW_ASSERT(
      info == 0,
      "LSS_no_alloc::factorize( who={},...) call lapack_wrapper::gecopy return info = {}\n",
      who, info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LSS_no_alloc<T>::factorize( valueType const A[], integer LDA ) {
    integer info = gecopy(
      m_nRows, m_nCols, A, LDA, m_Amat, m_nRows
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LSS_no_alloc<T>::solve( valueType xb[] ) const {
    // check memory
    integer L = getL( m_nRows, m_nCols, 1 );
    if ( L > m_Lwork ) {
      m_allocWork.allocate( L );
      m_Lwork = L;
      m_Work  = m_allocWork( L );
    }
    // save matrix
    copy( m_nRows*m_nCols, m_Amat, 1, m_AmatWork, 1);
    integer info = gelss(
      m_nRows, m_nCols, 1,
      m_AmatWork, m_nRows,
      xb, m_nRows,
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
    integer L = getL( m_nCols, m_nRows, 1 );
    if ( L > m_Lwork ) {
      m_allocWork.allocate( L );
      m_Lwork = L;
      m_Work  = m_allocWork( L );
    }
    // save matrix
    for ( integer i = 0; i < m_nCols; ++i )
      copy( m_nRows, m_Amat+i*m_nRows, 1, m_AmatWork+i, m_nCols );
    integer info = gelss(
      m_nCols, m_nRows, 1,
      m_AmatWork, m_nCols,
      xb, m_nCols,
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
    integer L = getL( m_nRows, m_nCols, nrhs );
    if ( L > m_Lwork ) {
      m_allocWork.allocate( L );
      m_Lwork = L;
      m_Work  = m_allocWork( L );
    }
    // save matrix
    copy( m_nRows*m_nCols, m_Amat, 1, m_AmatWork, 1 );
    integer info = gelss(
      m_nRows, m_nCols, nrhs,
      m_AmatWork, m_nRows,
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
    integer L = getL( m_nCols, m_nRows, nrhs );
    if ( L > m_Lwork ) {
      m_allocWork.allocate( L );
      m_Lwork = L;
      m_Work  = m_allocWork( L );
    }
    // save matrix
    for ( integer i = 0; i < m_nCols; ++i )
      copy( m_nRows, m_Amat+i*m_nRows, 1, m_AmatWork+i, m_nCols );

    integer info = gelss(
      m_nCols, m_nRows, nrhs,
      m_AmatWork, m_nCols,
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
    if ( m_nRows != NR || m_nCols != NC ) {
      integer Lwork = 2*NR*NC+std::min(NR,NC);
      m_allocReals.allocate( size_t(Lwork) );
      this->no_allocate( NR, NC, Lwork, m_allocReals( size_t(Lwork) ) );
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
    LW_ASSERT(
      Lwork >= 2*NRC && Liwork >= NC,
      "LSS_no_alloc( NR = {}, NC = {}, LWork = {}, ..., Liwork = {})\n"
      "LWork must be >= {} and Liwork must be >= {}\n",
      NR, NC, Lwork, Liwork, Lwmin, NC
    );
    m_nRows    = NR;
    m_nCols    = NC;
    m_Amat     = Work;
    m_AmatWork = Work+NRC;
    m_jpvt     = iWork;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  integer
  LSY_no_alloc<T>::getL( integer NR, integer NC, integer nrhs ) const {
    valueType tmp;
    integer info = gelsy(
      NR, NC, nrhs,
      nullptr, NR, nullptr, NR, nullptr,
      m_rcond, m_rank, &tmp, -1
    );
    LW_ASSERT( info == 0, "LSY_no_alloc::getL, in gelss info = {}\n", info );
    return integer(tmp);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY_no_alloc<T>::factorize(
    char const      who[],
    valueType const A[],
    integer         LDA
  ) {
    integer info = gecopy(
      m_nRows, m_nCols, A, LDA, m_Amat, m_nRows
    );
    LW_ASSERT(
      info == 0,
      "LSY::factorize[{}] call lapack_wrapper::gecopy return info = {}\n",
      who, info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LSY_no_alloc<T>::factorize( valueType const A[], integer LDA ) {
    integer info = gecopy(
      m_nRows, m_nCols, A, LDA, m_Amat, m_nRows
    );
    return info == 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LSY_no_alloc<T>::solve( valueType xb[] ) const {
    // check memory
    integer L = getL( m_nRows, m_nCols, 1 );
    if ( L > m_Lwork ) {
      m_allocWork.allocate( L );
      m_Lwork = L;
      m_Work  = m_allocWork( L );
    }
    // save matrix
    copy( m_nRows*m_nCols, m_Amat, 1, m_AmatWork, 1 );
    integer info = gelsy(
      m_nRows, m_nCols, 1,
      m_AmatWork, m_nRows,
      xb, m_nRows, m_jpvt,
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
    integer L = getL( m_nCols, m_nRows, 1 );
    if ( L > m_Lwork ) {
      m_allocWork.allocate( L );
      m_Lwork = L;
      m_Work  = m_allocWork( L );
    }
    // save matrix
    for ( integer i = 0; i < m_nCols; ++i )
      copy( m_nRows, m_Amat+i*m_nRows, 1, m_AmatWork+i, m_nCols );
    integer info = gelsy(
      m_nCols, m_nRows, 1,
      m_AmatWork, m_nCols,
      xb, m_nCols, m_jpvt,
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
    integer L = getL( m_nRows, m_nCols, nrhs );
    if ( L > m_Lwork ) {
      m_allocWork.allocate( L );
      m_Lwork = L;
      m_Work  = m_allocWork( L );
    }
    // save matrix
    copy( m_nRows*m_nCols, m_Amat, 1, m_AmatWork, 1 );
    integer info = gelsy(
      m_nRows, m_nCols, nrhs,
      m_AmatWork, m_nRows,
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
    integer L = getL( m_nCols, m_nRows, nrhs );
    if ( L > m_Lwork ) {
      m_allocWork.allocate( L );
      m_Lwork = L;
      m_Work  = m_allocWork( L );
    }
    // save matrix
    for ( integer i = 0; i < m_nCols; ++i )
      copy( m_nRows, m_Amat+i*m_nRows, 1, m_AmatWork+i, m_nCols );
    integer info = gelsy(
      m_nCols, m_nRows, nrhs,
      m_AmatWork, m_nCols,
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
    if ( m_nRows != NR || m_nCols != NC ) {
      integer Lwork = 2*NR*NC;
      m_allocReals.allocate( size_t(Lwork) );
      m_allocInts.allocate( size_t(NC) );
      this->no_allocate(
        NR, NC,
        Lwork, m_allocReals( size_t(Lwork) ),
        NC, m_allocInts( size_t(NC) )
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
      m_nRows, m_nCols, A, LDA, m_Amat, m_nRows
    );
    LW_ASSERT(
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
      m_nRows, m_nCols, A, LDA, m_Amat, m_nRows
    );
    return info == 0;
  }

}

///
/// eof: ls.cxx
///
