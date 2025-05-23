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
 |      Università degli Studi di Trento                                    |
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
    using real_type = T;

  private:

    Malloc<real_type> m_allocReals{"allocReals"};

    real_type * m_L{nullptr};
    real_type * m_D{nullptr};
    real_type * m_WORK{nullptr};
    integer     m_nRC{0};

  public:

    using LinearSystemSolver<T>::factorize;

    TridiagonalSPD() = default;

    ~TridiagonalSPD() override
    { m_allocReals.free(); }

    real_type cond1( real_type norm1 ) const;

    void
    factorize(
      string_view     who,
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

    static
    void
    axpy(
      integer         N,
      real_type       alpha,
      real_type const L[],
      real_type const D[],
      real_type const x[],
      real_type       beta,
      real_type       y[]
    );

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
    using real_type = T;

  private:

    Malloc<real_type> m_allocReals{"TridiagonalLU-allocReals"};
    Malloc<integer>   m_allocIntegers{"TridiagonalLU-allocIntegers"};

    real_type * m_L{nullptr};
    real_type * m_D{nullptr};
    real_type * m_U{nullptr};
    real_type * m_U2{nullptr};
    real_type * m_WORK{nullptr};
    integer   * m_IPIV{nullptr};
    integer   * m_IWORK{nullptr};

    integer m_nRC{0};

  public:

    using LinearSystemSolver<T>::factorize;

    TridiagonalLU() = default;

    ~TridiagonalLU() override {
      m_allocReals.free();
      m_allocIntegers.free();
    }

    real_type cond1( real_type norm1 ) const;
    real_type condInf( real_type normInf ) const;

    void
    factorize(
      string_view     who,
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
    void solve( string_view who, real_type xb[] ) const override;
    bool t_solve( real_type xb[] ) const override;
    void t_solve( string_view who, real_type xb[] ) const override;

    bool
    solve(
      integer   nrhs,
      real_type xb[],
      integer   ldXB
    ) const override;

    void
    solve(
      string_view who,
      integer     nrhs,
      real_type   xb[],
      integer     ldXB
    ) const override;

    bool
    t_solve(
      integer   nrhs,
      real_type xb[],
      integer   ldXB
    ) const override;

    void
    t_solve(
      string_view who,
      integer     nrhs,
      real_type   xb[],
      integer     ldXB
    ) const override;

    /*\
    :|:     _
    :|:    / \  _   ___  __
    :|:   / _ \| | | \ \/ /
    :|:  / ___ \ |_| |>  <
    :|: /_/   \_\__,_/_/\_\
    :|:
    \*/

    static
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
    );

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
    using real_type = T;

  private:

    Malloc<real_type> m_allocReals{"allocReals"};

    real_type * m_C{nullptr};   // rotazioni givens
    real_type * m_S{nullptr};
    real_type * m_BD{nullptr};  // band triangular matrix
    real_type * m_BU{nullptr};  // band triangular matrix
    real_type * m_BU2{nullptr}; // band triangular matrix

    real_type   m_normInfA{0};
    integer     m_nRC{0};

    void Rsolve( real_type xb[] ) const;
    void RsolveTransposed( real_type xb[] ) const;

  public:

    using LinearSystemSolver<T>::factorize;

    TridiagonalQR() = default;

    ~TridiagonalQR() override {
      m_allocReals.free();
    }

    void
    factorize(
      string_view     who,
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

    static
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
    );

  };

}

///
/// eof: trid.hxx
///
