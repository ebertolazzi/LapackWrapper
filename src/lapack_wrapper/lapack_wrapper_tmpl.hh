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
/// file: Alglin_tmpl.hh
///

#pragma once

#ifndef LAPACK_WRAPPER_TMPL_HH
#define LAPACK_WRAPPER_TMPL_HH

#include "lapack_wrapper.hh"

namespace lapack_wrapper {

  /*\
   |  __     __        _
   |  \ \   / /__  ___/ |
   |   \ \ / / _ \/ __| |
   |    \ V /  __/ (__| |
   |     \_/ \___|\___|_|
  \*/

  template <typename t_Value, int N, int STEP>
  class Vec1 {
  public:
    static inline void zero( t_Value * A ) {
      A[0] = 0;
      Vec1<t_Value,N-1,STEP>::zero( A+STEP );
    }
    static inline void scal( t_Value s, t_Value * x ) {
      x[0] *= s;
      Vec1<t_Value,N-1,STEP>::scal(s,x+STEP);
    }
    static inline t_Value sum( t_Value const * x ) {
      return x[0] + Vec1<t_Value,N-1,STEP>::sum(x+STEP);
    }
    static inline t_Value asum( t_Value const * x ) {
      return std::abs(x[0]) + Vec1<t_Value,N-1,STEP>::asum(x+STEP);
    }
    static inline t_Value max( t_Value const * x ) {
      return std::max( x[0], Vec1<t_Value,N-1,STEP>::max(x+STEP) );
    }
    static inline t_Value min( t_Value const * x ) {
      return std::min( x[0], Vec1<t_Value,N-1,STEP>::min(x+STEP) );
    }
    static inline t_Value amax( t_Value const * x ) {
      return std::max( std::abs(x[0]), Vec1<t_Value,N-1,STEP>::amax(x+STEP) );
    }
    static inline void amax( t_Value const * x, t_Value & am, integer & ipos ) {
      integer ipos1;
      t_Value am1;
      am = std::abs(x[0]);
      Vec1<t_Value,N-1,STEP>::amax(x+STEP,am1,ipos1);
      if ( am < am1 ) { am = am1; ipos = ipos1+1; }
      else ipos = 0;
    }
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value, int STEP>
  class Vec1<t_Value,1,STEP> {
  public:
    static inline void zero( t_Value * A ) { A[0] = 0; }
    static inline void scal( t_Value s, t_Value * A ) { A[0] *= s; }
    static inline t_Value sum( t_Value const * x ) { return x[0]; }
    static inline t_Value asum( t_Value const * x ) { return std::abs(x[0]); }
    static inline t_Value max( t_Value const * x ) { return x[0]; }
    static inline t_Value min( t_Value const * x ) { return x[0]; }
    static inline t_Value amax( t_Value const * x ) { return std::abs(x[0]); }
    static inline void amax( t_Value const * x, t_Value & am, integer & ipos )
    { am = std::abs(x[0]); ipos = 0; }
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value, int STEP>
  class Vec1<t_Value,0,STEP> {
  public:
    static inline void zero( t_Value * ) {}
    static inline void scal( t_Value, t_Value * ) {}
    static inline t_Value sum ( t_Value const * ) { return 0; }
    static inline t_Value asum( t_Value const * ) { return 0; }
    static inline t_Value max ( t_Value const * ) { return 0; }
    static inline t_Value min ( t_Value const * ) { return 0; }
    static inline t_Value amax( t_Value const * ) { return 0; }
    static inline void amax( t_Value const *, t_Value & am, integer & ipos )
    { am = 0; ipos = -1; }
  };

  /*\
   |  __     __        ____
   |  \ \   / /__  ___|___ \
   |   \ \ / / _ \/ __| __) |
   |    \ V /  __/ (__ / __/
   |     \_/ \___|\___|_____|
  \*/

  template <typename t_Value, int N, int STEPA = 1, int STEPB = 1>
  class Vec2 {
  public:
    static inline void copy( t_Value const * A, t_Value * B ) {
      B[0] = A[0];
      Vec2<t_Value,N-1,STEPA,STEPB>::copy( A+STEPA, B+STEPB );
    }
    static inline void swap( t_Value * A, t_Value * B ) {
      std::swap(B[0],A[0]);
      Vec2<t_Value,N-1,STEPA,STEPB>::swap( A+STEPA, B+STEPB );
    }
    static inline t_Value dot( t_Value const * A, t_Value const * B ) {
      return B[0]*A[0]+Vec2<t_Value,N-1,STEPA,STEPB>::dot( A+STEPA, B+STEPB );
    }
    static inline void addto( t_Value const * x, t_Value * y ) {
      y[0] += x[0];
      Vec2<t_Value,N-1,STEPA,STEPB>::addto(x+STEPA,y+STEPB);
    }
    static inline void subto( t_Value const * x, t_Value * y ) {
      y[0] -= x[0];
      Vec2<t_Value,N-1,STEPA,STEPB>::subto(x+STEPA,y+STEPB);
    }
    static inline void assax( t_Value a, t_Value const * x, t_Value * y ) {
      y[0] = a*x[0];
      Vec2<t_Value,N-1,STEPA,STEPB>::axy(a,x+STEPA,y+STEPB);
    }
    static inline void axpy( t_Value a, t_Value const * x, t_Value * y ) {
      y[0] += a*x[0];
      Vec2<t_Value,N-1,STEPA,STEPB>::axpy(a,x+STEPA,y+STEPB);
    }
    static inline void axpby( t_Value a, t_Value const * x, t_Value b, t_Value * y ) {
      y[0] = b*y[0] + a*x[0];
      Vec2<t_Value,N-1,STEPA,STEPB>::axpby(a,x+STEPA,b,y+STEPB);
    }
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value, int STEPA, int STEPB>
  class Vec2<t_Value,0,STEPA,STEPB> {
  public:
    static inline void copy( t_Value const *, t_Value * ) {}
    static inline void swap( t_Value *, t_Value * ) {}
    static inline t_Value dot( t_Value const *, t_Value const * ) { return 0; }
    static inline void addto( t_Value const *, t_Value * ) {}
    static inline void subto( t_Value const *, t_Value * ) {}
    static inline void assax( t_Value, t_Value const *, t_Value * ) {}
    static inline void axpy( t_Value, t_Value const *, t_Value * ) {}
    static inline void axpby( t_Value, t_Value const *, t_Value, t_Value *) {}
  };

  /*
  //   _               _
  //  | |    ___  ___ | |_   _____
  //  | |   / __|/ _ \| \ \ / / _ \
  //  | |___\__ \ (_) | |\ V /  __/
  //  |_____|___/\___/|_| \_/ \___|
  */
  template <typename t_Value, int N, int LDL>
  class LsolveUnit {
  public:
    static inline void eval( t_Value const * L, t_Value * x ) {
      Vec2<t_Value,N-1,1,1>::axpy(-x[0],L+1,x+1);
      LsolveUnit<t_Value,N-1,LDL>::eval(L+LDL+1,x+1);
    }
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value, int LDL>
  class LsolveUnit<t_Value,1,LDL> {
  public:
    static inline void eval( t_Value const *, t_Value * ) {}
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value, int N, int LDL>
  class Lsolve {
  public:
    static inline void eval( t_Value const * L, t_Value * x ) {
      x[0] /= L[0];
      Vec2<t_Value,N-1,1,1>::axpy(-x[0],L+1,x+1);
      Lsolve<t_Value,N-1,LDL>::eval(L+LDL+1,x+1);
    }
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value, int LDL>
  class Lsolve<t_Value,1,LDL> {
  public:
    static inline void eval( t_Value const * L, t_Value * x ) { x[0] /= L[0]; }
  };

  /*
  //   _   _           _
  //  | | | |___  ___ | |_   _____
  //  | | | / __|/ _ \| \ \ / / _ \
  //  | |_| \__ \ (_) | |\ V /  __/
  //   \___/|___/\___/|_| \_/ \___|
  */
  template <typename t_Value, int N, int LDU>
  class UsolveUnit {
  public:
    static inline void eval( t_Value const * U, t_Value * x ) {
      Vec2<t_Value,N-1,1,1>::axpy(-x[N-1],U+(N-1)*LDU,x);
      UsolveUnit<t_Value,N-1,LDU>::eval(U,x);
    }
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value, int LDU>
  class UsolveUnit<t_Value,1,LDU> {
  public:
    static void eval( t_Value const *, t_Value * ) {}
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value, int N, int LDU>
  class Usolve {
  public:
    static inline void eval( t_Value const * U, t_Value * x ) {
      x[N-1] /= U[(N-1)*(LDU+1)];
      Vec2<t_Value,N-1,1,1>::axpy(-x[N-1],U+(N-1)*LDU,x);
      Usolve<t_Value,N-1,LDU>::eval(U,x);
    }
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value, int LDU>
  class Usolve<t_Value,1,LDU> {
  public:
    static inline void eval( t_Value const * U, t_Value * x ) {
      x[0] /= U[0];
    }
  };
  
  /*
  //   __  __
  //  |  \/  |_   __
  //  | |\/| \ \ / /
  //  | |  | |\ V /
  //  |_|  |_| \_/
  */

  template <typename t_Value, int M, int N, int LDM, int INCX, int INCR>
  class Mv {
  public:
    static inline void ass( t_Value const * Mat, t_Value const * x, t_Value * res ) {
      Vec2<t_Value,M,1,INCR>::assax( x[0], Mat, res );
      Mv<t_Value,M,N-1,LDM,INCX,INCR>::addTo( Mat+LDM, x+INCX, res );
    }
    static inline void addTo( t_Value const * Mat, t_Value const * x, t_Value * res ) {
      Vec2<t_Value,M,1,INCR>::axpy( x[0], Mat, res );
      Mv<t_Value,M,N-1,LDM,INCX,INCR>::addTo( Mat+LDM, x+INCX, res );
    }
    static inline void subTo( t_Value const * Mat, t_Value const * x, t_Value * res ) {
      Vec2<t_Value,M,1,INCR>::axpy( -x[0], Mat, res );
      Mv<t_Value,M,N-1,LDM,INCX,INCR>::subTo( Mat+LDM, x+INCX, res );
    }
    static inline void aMxpby( t_Value a, t_Value const * Mat, t_Value const * x,
                               t_Value b, t_Value * res ) {
      Vec2<t_Value,M,1,INCR>::axby( a*x[0], Mat, b, res );
      Mv<t_Value,M,N-1,LDM,INCX,INCR>::aMxpby( a, Mat+LDM, x+INCX, b, res );
    }
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value, int N, int LDM, int INCX, int INCR>
  class Mv<t_Value,0,N,LDM,INCX,INCR> {
  public:
    static inline void ass( t_Value const *, t_Value const *, t_Value * ) { }
    static inline void addTo( t_Value const *, t_Value const *, t_Value * ) { }
    static inline void subTo( t_Value const *, t_Value const *, t_Value * ) { }
    static inline void aMxpby( t_Value, t_Value const *, t_Value const *,
                               t_Value, t_Value * ) { }
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value, int M, int LDM, int INCX, int INCR>
  class Mv<t_Value,M,0,LDM,INCX,INCR> {
  public:
    static inline void ass( t_Value const *, t_Value const *, t_Value * ) { }
    static inline void addTo( t_Value const *, t_Value const *, t_Value * ) { }
    static inline void subTo( t_Value const *, t_Value const *, t_Value * ) { }
    static inline void aMxpby( t_Value, t_Value const *, t_Value const *,
                               t_Value, t_Value * ) { }
  };
  
  /*
  //   __  __ __  __
  //  |  \/  |  \/  |
  //  | |\/| | |\/| |
  //  | |  | | |  | |
  //  |_|  |_|_|  |_|
  //
  //  C  = A*B
  //  C += A*B
  //  C -= A*B
  //  C  = a*C + b*A*B
  */

  template <typename t_Value, int M, int N, int K, int LDA, int LDB, int LDC>
  class MM {
  public:
    static inline void ass( t_Value const * A, t_Value const * B, t_Value * C ) {
      for ( integer i = 0; i < N; ++i, B += LDB, C += LDC )
        Mv<t_Value,M,K,LDA,1,1>::ass( A, B, C );
    }
    static inline void addTo( t_Value const * A, t_Value const * B, t_Value * C ) {
      for ( integer i = 0; i < N; ++i, B += LDB, C += LDC )
        Mv<t_Value,M,K,LDA,1,1>::addTo( A, B, C );
    }
    static inline void subTo( t_Value const * A, t_Value const * B, t_Value * C ) {
      for ( integer i = 0; i < N; ++i, B += LDB, C += LDC )
        Mv<t_Value,M,K,LDA,1,1>::subTo( A, B, C );
    }
    static inline void aMxpby( t_Value a, t_Value const * A, t_Value const * B,
                               t_Value b, t_Value * C ) {
      for ( integer i = 0; i < N; ++i, B += LDB, C += LDC )
        Mv<t_Value,M,K,LDA,1,1>::aMxpby( a, A, B, b, C );
    }
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value, int M, int K, int LDA, int LDB, int LDC>
  class MM<t_Value,M,0,K,LDA,LDB,LDC> {
  public:
    static inline void ass( t_Value const *, t_Value const *, t_Value * ) {}
    static inline void addTo( t_Value const *, t_Value const *, t_Value * ) {}
    static inline void subTo( t_Value const *, t_Value const *, t_Value * ) {}
    static inline void aMxpby( t_Value, t_Value const *, t_Value const *, t_Value, t_Value *) {}
  };

  /*
  //   ____             _    _
  //  |  _ \ __ _ _ __ | | _/ |
  //  | |_) / _` | '_ \| |/ / |
  //  |  _ < (_| | | | |   <| |
  //  |_| \_\__,_|_| |_|_|\_\_|
  //
  //  M = u*v^T
  //  M += u*v^T
  //  M -= u*v^T
  //  M  = b*M + a*u*v^T
  //
  */

  template <typename t_Value, int M, int N, int LDM, int INCX, int INCY>
  class Rank1 {
  public:
    static inline void ass( t_Value const * x, t_Value const * y, t_Value * Mat ) {
      Vec2<t_Value,M,INCX,1>::assax( y[0], x, Mat );
      Rank1<t_Value,M,N-1,LDM,INCX,INCY>::ass(x,y+INCY,Mat+LDM);
    }
    static inline void addTo( t_Value const * x, t_Value const * y, t_Value * Mat ) {
      Vec2<t_Value,M,INCX,1>::axpy( y[0], x, Mat );
      Rank1<t_Value,M,N-1,LDM,INCX,INCY>::addTo(x,y+INCY,Mat+LDM);
    }
    static inline void subTo( t_Value const * x, t_Value const * y, t_Value * Mat ) {
      Vec2<t_Value,M,INCX,1>::axpy( -y[0], x, Mat );
      Rank1<t_Value,M,N-1,LDM,INCX,INCY>::subTo(x,y+INCY,Mat+LDM);
    }
    static inline void auvpbM( t_Value a, t_Value const * x, t_Value const * y,
                               t_Value b, t_Value * Mat ) {
      Vec2<t_Value,M,INCX,1>::axpby( a*y[0], x, b, Mat );
      Rank1<t_Value,M,N-1,LDM,INCX,INCY>::auvpbM(x,y+INCY,Mat+LDM);
    }
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value, int M, int LDM, int INCX, int INCY>
  class Rank1<t_Value,M,0,LDM,INCX,INCY> {
  public:
    static inline void ass( t_Value const *, t_Value const *, t_Value * ) { }
    static inline void addTo( t_Value const *, t_Value const *, t_Value * ) { }
    static inline void subTo( t_Value const *, t_Value const *, t_Value * ) { }
    static inline void auvpbM( t_Value, t_Value const *, t_Value const *,
                               t_Value, t_Value * ) { }
  };


} // end namespace alglin

#endif

///
/// eof: Alglin_tmpl.hh
///

