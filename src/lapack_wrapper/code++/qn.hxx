/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2019                                                      |
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
/// file: qn.hxx
///

namespace lapack_wrapper {

  /*\
  :|:    ___                  _ _   _               _
  :|:   / _ \ _   _  __ _ ___(_) \ | | _____      _| |_ ___  _ __
  :|:  | | | | | | |/ _` / __| |  \| |/ _ \ \ /\ / / __/ _ \| '_ \
  :|:  | |_| | |_| | (_| \__ \ | |\  |  __/\ V  V /| || (_) | | | |
  :|:   \__\_\\__,_|\__,_|___/_|_| \_|\___| \_/\_/  \__\___/|_| |_|
  \*/

  template <typename T>
  class QN {
  public:
    typedef T valueType;

  protected:
    Malloc<valueType> m_allocReals;

    integer     m_dim;
    valueType * m_H;
    valueType * m_s;
    valueType * m_y;
    valueType * m_z;

  public:

    QN()
    : m_allocReals("QN reals")
    , m_dim(0)
    , m_H(nullptr)
    , m_s(nullptr)
    , m_y(nullptr)
    , m_z(nullptr)
    {}

    virtual
    ~QN() {}

    void
    allocate( integer N );

    valueType const &
    operator () ( integer i, integer j ) const
    { return m_H[i+j*m_dim]; }

    valueType &
    operator () ( integer i, integer j )
    { return m_H[i+j*m_dim]; }

    void
    zero()
    { lapack_wrapper::gezero( m_dim, m_dim, m_H, m_dim ); }

    void
    init()
    { lapack_wrapper::geid( m_dim, m_dim, m_H, m_dim ); }

    // r <- beta * r + alpha * H * x
    void
    mult(
      valueType       alpha,
      valueType const x[],
      integer         inc_x,
      valueType       beta,
      valueType       r[],
      integer         inc_r
    ) const {
      symv(
        LOWER,
        m_dim,
        alpha, m_H, m_dim,
        x, inc_x,
        beta,
        r, inc_r
      );
    }

    // r <- H * x
    void
    mult( valueType const x[], valueType r[] ) const {
      mult( valueType(1), x, 1, valueType(0), r, 1 );
    }

    // r <- -H * x
    void
    minus_mult( valueType const x[], valueType r[] ) const {
      mult( valueType(-1), x, 1, valueType(0), r, 1 );
    }

    void
    update(
      valueType const f0[],
      valueType const f1[],
      valueType const ss[],
      valueType       epsi
    ) {
      copy( m_dim,     f1, 1, m_y, 1 );
      axpy( m_dim, -1, f0, 1, m_y, 1 );
      update( m_y, ss, epsi );
    }

    void
    update(
      valueType const f0[],
      valueType const f1[],
      valueType const x0[],
      valueType const x1[],
      valueType       epsi
    ) {
      copy( m_dim,     f1, 1, m_y, 1 );
      axpy( m_dim, -1, f0, 1, m_y, 1 );
      copy( m_dim,     x1, 1, m_s, 1 );
      axpy( m_dim, -1, x0, 1, m_s, 1 );
      update( m_y, m_s, epsi );
    }

    void
    print( ostream_type & stream ) const;

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    void
    update(
      valueType const y[],
      valueType const s[],
      valueType       epsi
    ) UTILS_PURE_VIRTUAL;

  };

  /*\
  :|:  ____  _____ ____ ____
  :|: | __ )|  ___/ ___/ ___|
  :|: |  _ \| |_ | |  _\___ \
  :|: | |_) |  _|| |_| |___) |
  :|: |____/|_|   \____|____/
  \*/

  template <typename T>
  class BFGS : public QN<T> {
  public:
    typedef T valueType;

    using QN<T>::allocate;
    using QN<T>::zero;
    using QN<T>::init;
    using QN<T>::mult;
    using QN<T>::update;

    using QN<T>::m_dim;
    using QN<T>::m_H;
    using QN<T>::m_z;

    BFGS() {}

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    void
    update(
      valueType const y[],
      valueType const s[],
      valueType       epsi
    ) UTILS_OVERRIDE;

  };

  /*\
  :|:   ____  _____ ____
  :|:  |  _ \|  ___|  _ \
  :|:  | | | | |_  | |_) |
  :|:  | |_| |  _| |  __/
  :|:  |____/|_|   |_|
  \*/

  template <typename T>
  class DFP : public QN<T> {
  public:
    typedef T valueType;

    using QN<T>::allocate;
    using QN<T>::zero;
    using QN<T>::init;
    using QN<T>::mult;
    using QN<T>::update;

    using QN<T>::m_dim;
    using QN<T>::m_H;
    using QN<T>::m_z;

    DFP() {}

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    void
    update(
      valueType const y[],
      valueType const s[],
      valueType       epsi
    ) UTILS_OVERRIDE;
  };

  /*\
  :|:   ____  ____  _
  :|:  / ___||  _ \/ |
  :|:  \___ \| |_) | |
  :|:   ___) |  _ <| |
  :|:  |____/|_| \_\_|
  \*/

  template <typename T>
  class SR1 : public QN<T> {
  public:
    typedef T valueType;

    using QN<T>::allocate;
    using QN<T>::zero;
    using QN<T>::init;
    using QN<T>::mult;
    using QN<T>::update;

    using QN<T>::m_dim;
    using QN<T>::m_H;
    using QN<T>::m_z;

    SR1() {}

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    void
    update(
      valueType const y[],
      valueType const s[],
      valueType       epsi
    ) UTILS_OVERRIDE;
  };

}

///
/// eof: qn.hxx
///
