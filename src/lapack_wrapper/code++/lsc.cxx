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

  //============================================================================
  /*\
  :|:   _     ____   ____
  :|:  | |   / ___| / ___|
  :|:  | |   \___ \| |
  :|:  | |___ ___) | |___
  :|:  |_____|____/ \____|
  \*/
  /*!  least square constained
  :|:
  :|:  minimize    (1/2) * || A * x - b ||^2
  :|:
  :|:  subject to  B * x = c
  :|:
  :|:  L = (1/2) * || A * x - b ||^2 + lambda * ( B * x - c )
  :|:
  :|:  is equivalent to solve linear system
  :|:
  :|:  / A^T*A  B^T \  /   x    \   / A^T b \
  :|:  |            |  |        | = |       |
  :|:  \ B       0  /  \ lambda /   \   c   /
  :|:
  :|:  is equivalent to solve linear system
  :|:
  :|:  / 0   A^T  B^T \  /   x    \   / A^T b \
  :|:  |              |  |        |   |       |
  :|:  | A   -I    0  |  | omega  | = |   0   |
  :|:  |              |  |        |   |       |
  :|:  \ B    0    0  /  \ lambda /   \   c   /
  :|:
  \*/

  template <typename T>
  LSC<T>::LSC()
  : LinearSystemSolver<T>()
  , m_allocReals( "LSC-allocReals" )
  , m_allocIntegers( "LSC-allocIntegers" )
  , m_NR(0)
  , m_NC(0)
  , m_NRA(0)
  , m_NRB(0)
  , m_to_rowA(nullptr)
  , m_to_rowB(nullptr)
  , m_Amat(nullptr)
  , m_lu()
  , m_allocWorks( "LSC-allocWorks" )
  , m_work(nullptr)
  , m_rhs(nullptr)
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LSC<T>::factorize(
    integer         NR,
    integer         NC,
    valueType const M[],
    integer         ldM,
    bool const      row_select[] // row selected for minimization
  ) {
    m_NR = NR;
    m_NC = NC;
    // dimension of the big system
    integer NN  = NR+NC;
    integer NN2 = NN*NN;

    // count row
    m_NRA = m_NRB = 0;
    for ( integer i = 0; i < NR; ++i ) {
      if ( row_select[i] ) ++m_NRA; else ++m_NRB;
    }

    m_allocReals.allocate( size_t(m_NRA*NC) );
    m_Amat = m_allocReals( size_t(m_NRA*NC) );

    m_allocIntegers.allocate( size_t(NR) );
    m_to_rowA = m_allocIntegers( size_t(m_NRA) );
    m_to_rowB = m_allocIntegers( size_t(m_NRB) );

    m_allocWorks.allocate( size_t(NN2) );
    m_work = m_allocWorks( size_t(NN2) );

    integer LDA = m_NRA;
    m_NRA = m_NRB = 0;
    for ( integer i = 0; i < NR; ++i ) {
      if ( row_select[i] ) {
        lapack_wrapper::copy( NC, M+i, ldM, m_Amat+m_NRA, LDA );
        m_to_rowA[m_NRA++] = i;
      } else {
        m_to_rowB[m_NRB++] = i;
      }
    }

    /*\
    :|:  / 0   A^T  B^T \
    :|:  |              |
    :|:  | A   -I    0  |
    :|:  |              |
    :|:  \ B    0    0  /
    \*/

    lapack_wrapper::zero( NN2, m_work, 1 );

    // A
    for ( integer i = 0; i < m_NRA; ++i ) {
      /*\
      :|:  / 0   A^T \
      :|:  |         |
      :|:  \ A   -I  /
      \*/
      integer           NCi = NC + i;
      valueType const * Mi  = M + m_to_rowA[i];;
      lapack_wrapper::copy( NC, Mi, ldM, m_work+NCi,    NN );
      lapack_wrapper::copy( NC, Mi, ldM, m_work+NCi*NN, 1  );
      m_work[ NCi*(NN+1) ] = -1;
    }
    // B
    for ( integer i = 0; i < m_NRB; ++i ) {
      integer           NCi = NC + m_NRA + i;
      valueType const * Mi  = M + m_to_rowB[i];
      lapack_wrapper::copy( NC, Mi, ldM, m_work+NCi,    NN );
      lapack_wrapper::copy( NC, Mi, ldM, m_work+NCi*NN, 1  );
    }
    return m_lu.factorize( NN, NN, m_work, NN );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LSC<T>::solve( valueType xb[] ) const {

    m_allocWorks.allocate( size_t(m_NRA+m_NR+m_NC) );
    m_work = m_allocWorks( size_t(m_NRA) );
    m_rhs  = m_allocWorks( size_t(m_NR+m_NC) );

    /*\
    :|:  / 0   A^T  B^T \  /   x    \   / A^T b \
    :|:  |              |  |        |   |       |
    :|:  | A   -I    0  |  | omega  | = |   0   |
    :|:  |              |  |        |   |       |
    :|:  \ B    0    0  /  \ lambda /   \   c   /
    \*/
    // A^T b
    if ( m_NRA > 0 ) {
      for ( integer i = 0; i < m_NRA; ++i ) m_work[i] = xb[m_to_rowA[i]];
      lapack_wrapper::gemv(
        lapack_wrapper::TRANSPOSE, m_NRA, m_NC,
        1.0, m_Amat, m_NRA,
        m_work, 1,
        0.0, m_rhs, 1
      );
      lapack_wrapper::zero( m_NRA, m_rhs+m_NC, 1 );
    }

    for ( integer i = 0; i < m_NRB; ++i )
      m_rhs[m_NC+m_NRA+i] = xb[m_to_rowB[i]];

    bool ok = m_lu.solve( m_rhs );
    if ( ok ) lapack_wrapper::copy( m_NC, m_rhs, 1, xb, 1 );
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  bool
  LSC<T>::solve(
    integer   nrhs,
    valueType B[],
    integer   ldB
  ) const {

    m_allocWorks.allocate( size_t( (m_NRA+m_NR+m_NC)*nrhs ) );
    m_work = m_allocWorks( size_t( m_NRA*nrhs ) );
    m_rhs  = m_allocWorks( size_t( (m_NR+m_NC)*nrhs ) );

    /*\
    :|:  / 0   A^T  B^T \  /   x    \   / A^T b \
    :|:  |              |  |        |   |       |
    :|:  | A   -I    0  |  | omega  | = |   0   |
    :|:  |              |  |        |   |       |
    :|:  \ B    0    0  /  \ lambda /   \   c   /
    \*/
    // A^T b
    for ( integer i = 0; i < m_NRA; ++i )
      lapack_wrapper::copy( nrhs, B+m_to_rowA[i], ldB, m_work+i, m_NRA );

/*
    lapack_wrapper::gemm(
      lapack_wrapper::TRANSPOSE,
      lapack_wrapper::NO_TRANSPOSE,
      m_NRA, m_NC,
      integer               K,
      1.0, m_Amat, m_NRA,
      B, ldB,
      0.0,
      real                  C[],
      integer               LDC
    );



    lapack_wrapper::gemv(
      lapack_wrapper::TRANSPOSE, m_NRA, m_NC,
      1.0, m_Amat, m_NRA,
      m_work, 1,
      0.0, m_rhs, 1
    );

    lapack_wrapper::zero( m_NRA, m_rhs+m_NC, 1 );

    for ( integer i = 0; i < m_NRB; ++i )
      m_rhs[m_NC+m_NRA+i] = xb[m_to_rowB[i]];
*/

    return true;
  }
}

///
/// eof: lsc.cxx
///
