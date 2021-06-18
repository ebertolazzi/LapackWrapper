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
/// file: trid.hxx
///

namespace lapack_wrapper {

  //============================================================================
  /*\
  :|:   _____     _     _ _                               _ ____  ____  ____
  :|:  |_   _| __(_) __| (_) __ _  __ _  ___  _ __   __ _| / ___||  _ \|  _ \
  :|:    | || '__| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | \___ \| |_) | | | |
  :|:    | || |  | | (_| | | (_| | (_| | (_) | | | | (_| | |___) |  __/| |_| |
  :|:    |_||_|  |_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|____/|_|   |____/
  :|:                             |___/
  \*/

  template <typename T>
  class TridiagonalSPD : public LinearSystemSolver<T> {
  public:
    typedef T real_type;

  private:

    Malloc<real_type> m_allocReals;

    real_type * m_L;
    real_type * m_D;
    real_type * m_WORK;
    integer     m_nRC;

  public:

    using LinearSystemSolver<T>::factorize;

    TridiagonalSPD()
    : m_allocReals("allocReals")
    , m_nRC(0)
    {}

    ~TridiagonalSPD() override
    { m_allocReals.free(); }

    real_type cond1( real_type norm1 ) const;

    void
    factorize(
      char const      who[],
      integer         N,
      real_type const L[],
      real_type const D[]
    );

    bool
    factorize(
      integer         N,
      real_type const L[],
      real_type const D[]
    );

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    bool solve( real_type xb[] ) const override;
    bool t_solve( real_type xb[] ) const override;

    bool
    solve(
      integer   nrhs,
      real_type xb[],
      integer   ldXB
    ) const override;

    bool
    t_solve(
      integer   nrhs,
      real_type xb[],
      integer   ldXB
    ) const override;

    /*\
    :|:     _
    :|:    / \  _   ___  __
    :|:   / _ \| | | \ \/ /
    :|:  / ___ \ |_| |>  <
    :|: /_/   \_\__,_/_/\_\
    :|:
    \*/

    void
    axpy(
      integer         N,
      real_type       alpha,
      real_type const L[],
      real_type const D[],
      real_type const x[],
      real_type       beta,
      real_type       y[]
    ) const;

  };

  //============================================================================
  /*\
  :|:   _____     _     _ _                               _ _    _   _
  :|:  |_   _| __(_) __| (_) __ _  __ _  ___  _ __   __ _| | |  | | | |
  :|:    | || '__| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | | |  | | | |
  :|:    | || |  | | (_| | | (_| | (_| | (_) | | | | (_| | | |__| |_| |
  :|:    |_||_|  |_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|_____\___/
  :|:                             |___/
  \*/

  template <typename T>
  class TridiagonalLU : public LinearSystemSolver<T> {
  public:
    typedef T real_type;

  private:

    Malloc<real_type> m_allocReals;
    Malloc<integer>   m_allocIntegers;

    real_type * m_L;
    real_type * m_D;
    real_type * m_U;
    real_type * m_U2;
    real_type * m_WORK;
    integer   * m_IPIV;
    integer   * m_IWORK;

    integer     m_nRC;

  public:

    using LinearSystemSolver<T>::factorize;

    TridiagonalLU()
    : m_allocReals("TridiagonalLU-allocReals")
    , m_allocIntegers("TridiagonalLU-allocIntegers")
    , m_nRC(0)
    {}

    ~TridiagonalLU() override {
      m_allocReals.free();
      m_allocIntegers.free();
    }

    real_type cond1( real_type norm1 ) const;
    real_type condInf( real_type normInf ) const;

    void
    factorize(
      char const      who[],
      integer         N,
      real_type const L[],
      real_type const D[],
      real_type const U[]
    );

    bool
    factorize(
      integer         N,
      real_type const L[],
      real_type const D[],
      real_type const U[]
    );

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    bool solve( real_type xb[] ) const override;
    void solve( char const who[], real_type xb[] ) const override;
    bool t_solve( real_type xb[] ) const override;
    void t_solve( char const who[], real_type xb[] ) const override;

    bool
    solve(
      integer   nrhs,
      real_type xb[],
      integer   ldXB
    ) const override;

    void
    solve(
      char const who[],
      integer    nrhs,
      real_type  xb[],
      integer    ldXB
    ) const override;

    bool
    t_solve(
      integer   nrhs,
      real_type xb[],
      integer   ldXB
    ) const override;

    void
    t_solve(
      char const who[],
      integer    nrhs,
      real_type  xb[],
      integer    ldXB
    ) const override;

    /*\
    :|:     _
    :|:    / \  _   ___  __
    :|:   / _ \| | | \ \/ /
    :|:  / ___ \ |_| |>  <
    :|: /_/   \_\__,_/_/\_\
    :|:
    \*/

    void
    axpy(
      integer         N,
      real_type       alpha,
      real_type const L[],
      real_type const D[],
      real_type const U[],
      real_type const x[],
      real_type       beta,
      real_type       y[]
    ) const;

  };

  //============================================================================
  /*\
  :|:   _____     _     _ _                               _  ___  ____
  :|:  |_   _| __(_) __| (_) __ _  __ _  ___  _ __   __ _| |/ _ \|  _ \
  :|:    | || '__| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | | | | | |_) |
  :|:    | || |  | | (_| | | (_| | (_| | (_) | | | | (_| | | |_| |  _ <
  :|:    |_||_|  |_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|\__\_\_| \_\
  :|:                             |___/
  \*/
  template <typename T>
  class TridiagonalQR : public LinearSystemSolver<T> {
  public:
    typedef T real_type;

  private:

    Malloc<real_type> m_allocReals;

    real_type * m_C;   // rotazioni givens
    real_type * m_S;
    real_type * m_BD;  // band triangular matrix
    real_type * m_BU;  // band triangular matrix
    real_type * m_BU2; // band triangular matrix

    real_type   m_normInfA;
    integer     m_nRC;

    void Rsolve( real_type xb[] ) const;
    void RsolveTransposed( real_type xb[] ) const;

  public:

    using LinearSystemSolver<T>::factorize;

    TridiagonalQR()
    : m_allocReals("allocReals")
    , m_nRC(0)
    {}

    ~TridiagonalQR() override {
      m_allocReals.free();
    }

    void
    factorize(
      char const      who[],
      integer         N,
      real_type const L[],
      real_type const D[],
      real_type const U[]
    );

    bool
    factorize(
      integer         N,
      real_type const L[],
      real_type const D[],
      real_type const U[]
    );

    void
    lsq(
      integer nrhs,
      T       RHS[],
      integer ldRHS,
      T       lambda
    ) const;

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    bool solve( real_type xb[] ) const override;
    bool t_solve( real_type xb[] ) const override;

    bool
    solve(
      integer   nrhs,
      real_type xb[],
      integer   ldXB
    ) const override;

    bool
    t_solve(
      integer   nrhs,
      real_type xb[],
      integer   ldXB
    ) const override;

    /*\
    :|:     _
    :|:    / \  _   ___  __
    :|:   / _ \| | | \ \/ /
    :|:  / ___ \ |_| |>  <
    :|: /_/   \_\__,_/_/\_\
    :|:
    \*/

    void
    axpy(
      integer         N,
      real_type       alpha,
      real_type const L[],
      real_type const D[],
      real_type const U[],
      real_type const x[],
      real_type       beta,
      real_type       y[]
    ) const;

  };

}

///
/// eof: trid.hxx
///
