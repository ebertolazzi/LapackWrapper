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
    m_dim = N;
    m_allocReals.allocate( size_t(m_dim*(m_dim+3)) );
    m_H = m_allocReals( size_t(m_dim*m_dim) );
    m_s = m_allocReals( size_t(m_dim) );
    m_y = m_allocReals( size_t(m_dim) );
    m_z = m_allocReals( size_t(m_dim) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QN<T>::print( ostream_type & stream ) const {
    for ( integer i = 0; i < m_dim; ++i ) {
      for ( integer j = 0; j < m_dim; ++j )
        fmt::print( stream, "{:14} ", m_H[i+j*m_dim] );
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
    valueType const s[],
    valueType       epsi
  ) {
    valueType sy   = dot( m_dim, s, 1, y, 1 );
    valueType slen = nrm2( m_dim, s, 1 );
    valueType ylen = nrm2( m_dim, y, 1 );
    if ( sy >= epsi*slen*ylen ) {
      mult( y, m_z );
      valueType yHy = dot( m_dim, m_z, 1, y, 1 );
      syr( LOWER, m_dim, (1+yHy/sy)/sy, s, 1, m_H, m_dim );
      syr2( LOWER, m_dim, -1/sy, s, 1, m_z, 1, m_H, m_dim );
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
    valueType const s[],
    valueType       epsi
  ) {
    valueType sy   = dot( m_dim, s, 1, y, 1 );
    valueType slen = nrm2( m_dim, s, 1 );
    valueType ylen = nrm2( m_dim, y, 1 );
    if ( sy >= epsi*slen*ylen ) {
      mult( y, m_z );
      valueType yHy = dot( m_dim, m_z, 1, y, 1 );
      syr( LOWER, m_dim, -1/yHy, m_z, 1, m_H, m_dim );
      syr( LOWER, m_dim,   1/sy, s,   1, m_H, m_dim );
    }
  }

  /*\
  :|:   ____  ____  _
  :|:  / ___||  _ \/ |
  :|:  \___ \| |_) | |
  :|:   ___) |  _ <| |
  :|:  |____/|_| \_\_|
  \*/

  template <typename T>
  void
  SR1<T>::update(
    valueType const y[],
    valueType const s[],
    valueType       epsi
  ) {
    // z <- y - H * s
    copy( m_dim, y, 1, m_z, 1 );
    mult( valueType(-1), s, 1, valueType(1), m_z, 1 );
    valueType sz   = dot( m_dim, s, 1, m_z, 1 );
    valueType slen = nrm2( m_dim, s, 1 );
    valueType zlen = nrm2( m_dim, m_z, 1 );
    if ( std::abs(sz) >= epsi*slen*zlen )
      syr( LOWER, m_dim, 1/sz, m_z, 1, m_H, m_dim );
  }
}

///
/// eof: qn.cxx
///
