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
    LAPACK_WRAPPER_ASSERT(
      N > 0 && N <= 1000,
      "QN<T>::allocate, N = " << N << " must be > 0 and <= 1000"
    );
    n = N;
    allocReals.allocate( size_t(n*(n+3)) );
    H = allocReals( size_t(n*n) );
    s = allocReals( size_t(n) );
    y = allocReals( size_t(n) );
    z = allocReals( size_t(n) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QN<T>::print( ostream_type & stream ) const {
    for ( integer i = 0; i < this->n; ++i ) {
      for ( integer j = 0; j < this->n; ++j )
        stream
          << std::setw(14)
          << this->H[i+j*this->n]
          << ' ';
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
    valueType sy = dot( n, s, 1, y, 1 );
    if ( sy > 0 ) {
      mult( y, z );
      valueType yHy = dot( n, z, 1, y, 1 );
      syr( LOWER, n, (1+yHy/sy)/sy, s, 1, H, n );
      syr2( LOWER, n, -1/sy, s, 1, z, 1, H, n );
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
    valueType sy = dot( n, s, 1, y, 1 );
    if ( sy > 0 ) {
      mult( y, z );
      valueType yHy = dot( n, z, 1, y, 1 );
      syr( LOWER, n, -1/yHy, z, 1, H, n );
      syr( LOWER, n,   1/sy, s, 1, H, n );
    }
  }

}

///
/// eof: qn.cxx
///
