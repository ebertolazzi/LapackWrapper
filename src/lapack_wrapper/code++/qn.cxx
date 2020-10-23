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
/// file: qn.cxx
///

namespace lapack_wrapper {

  /*\
   |    ___                  _ _   _               _
   |   / _ \ _   _  __ _ ___(_) \ | | _____      _| |_ ___  _ __
   |  | | | | | | |/ _` / __| |  \| |/ _ \ \ /\ / / __/ _ \| '_ \
   |  | |_| | |_| | (_| \__ \ | |\  |  __/\ V  V /| || (_) | | | |
   |   \__\_\\__,_|\__,_|___/_|_| \_|\___| \_/\_/  \__\___/|_| |_|
  \*/

  template <typename T>
  void
  QN<T>::allocate( integer N ) {
    UTILS_ASSERT(
      N > 0 && N <= 1000,
      "QN<T>::allocate, N = {} must be > 0 and <= 1000\n", N
    );
    m_n = N;
    m_allocReals.allocate( size_t(m_n*(m_n+3)) );
    m_H = m_allocReals( size_t(m_n*m_n) );
    m_s = m_allocReals( size_t(m_n) );
    m_y = m_allocReals( size_t(m_n) );
    m_z = m_allocReals( size_t(m_n) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QN<T>::print( ostream_type & stream ) const {
    for ( integer i = 0; i < m_n; ++i ) {
      for ( integer j = 0; j < m_n; ++j )
        fmt::print( stream, "{:14} ", m_H[i+j*m_n] );
      stream << '\n';
    }
  }

  /*\
   |  ____  _____ ____ ____
   | | __ )|  ___/ ___/ ___|
   | |  _ \| |_ | |  _\___ \
   | | |_) |  _|| |_| |___) |
   | |____/|_|   \____|____/
  \*/

  template <typename T>
  void
  BFGS<T>::update(
    valueType const y[],
    valueType const s[]
  ) {
    valueType sy = dot( m_n, s, 1, y, 1 );
    if ( sy > 0 ) {
      mult( y, m_z );
      valueType yHy = dot( m_n, m_z, 1, y, 1 );
      syr( LOWER, m_n, (1+yHy/sy)/sy, s, 1, m_H, m_n );
      syr2( LOWER, m_n, -1/sy, s, 1, m_z, 1, m_H, m_n );
    }
  }

  /*\
   |   ____  _____ ____
   |  |  _ \|  ___|  _ \
   |  | | | | |_  | |_) |
   |  | |_| |  _| |  __/
   |  |____/|_|   |_|
  \*/

  template <typename T>
  void
  DFP<T>::update(
    valueType const y[],
    valueType const s[]
  ) {
    valueType sy = dot( m_n, s, 1, y, 1 );
    if ( sy > 0 ) {
      mult( y, m_z );
      valueType yHy = dot( m_n, m_z, 1, y, 1 );
      syr( LOWER, m_n, -1/yHy, m_z, 1, m_H, m_n );
      syr( LOWER, m_n,   1/sy, s, 1, m_H, m_n );
    }
  }

}

///
/// eof: qn.cxx
///
