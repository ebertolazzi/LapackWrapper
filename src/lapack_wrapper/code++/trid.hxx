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
    typedef T valueType;

  private:

    Malloc<valueType> m_allocReals;

    valueType * m_L;
    valueType * m_D;
    valueType * m_WORK;
    integer     m_nRC;

  public:

    using LinearSystemSolver<T>::factorize;

    TridiagonalSPD()
    : m_allocReals("allocReals")
    , m_nRC(0)
    {}

    virtual
    ~TridiagonalSPD() LAPACK_WRAPPER_OVERRIDE
    { m_allocReals.free(); }

    valueType cond1( valueType norm1 ) const;

    void
    factorize(
      char const      who[],
      integer         N,
      valueType const _L[],
      valueType const _D[]
    );

    bool
    factorize(
      integer         N,
      valueType const _L[],
      valueType const _D[]
    );

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    bool
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    bool
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    bool
    solve(
      integer   nrhs,
      valueType xb[],
      integer   ldXB
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    bool
    t_solve(
      integer   nrhs,
      valueType xb[],
      integer   ldXB
    ) const LAPACK_WRAPPER_OVERRIDE;

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
      valueType       alpha,
      valueType const L[],
      valueType const D[],
      valueType const x[],
      valueType       beta,
      valueType       y[]
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
    typedef T valueType;

  private:

    Malloc<valueType> m_allocReals;
    Malloc<integer>   m_allocIntegers;

    valueType * m_L;
    valueType * m_D;
    valueType * m_U;
    valueType * m_U2;
    valueType * m_WORK;
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

    virtual
    ~TridiagonalLU() LAPACK_WRAPPER_OVERRIDE {
      m_allocReals.free();
      m_allocIntegers.free();
    }

    valueType cond1( valueType norm1 ) const;
    valueType condInf( valueType normInf ) const;

    void
    factorize(
      char const      who[],
      integer         N,
      valueType const _L[],
      valueType const _D[],
      valueType const _U[]
    );

    bool
    factorize(
      integer         N,
      valueType const _L[],
      valueType const _D[],
      valueType const _U[]
    );

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    bool
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    solve( char const who[], valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    bool
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve( char const who[], valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    bool
    solve(
      integer   nrhs,
      valueType xb[],
      integer   ldXB
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    solve(
      char const who[],
      integer    nrhs,
      valueType  xb[],
      integer    ldXB
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    bool
    t_solve(
      integer   nrhs,
      valueType xb[],
      integer   ldXB
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve(
      char const who[],
      integer    nrhs,
      valueType  xb[],
      integer    ldXB
    ) const LAPACK_WRAPPER_OVERRIDE;

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
      valueType       alpha,
      valueType const L[],
      valueType const D[],
      valueType const U[],
      valueType const x[],
      valueType       beta,
      valueType       y[]
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
    typedef T valueType;

  private:

    Malloc<valueType> m_allocReals;

    valueType * m_C;   // rotazioni givens
    valueType * m_S;
    valueType * m_BD;  // band triangular matrix
    valueType * m_BU;  // band triangular matrix
    valueType * m_BU2; // band triangular matrix

    valueType   m_normInfA;
    integer     m_nRC;

    void Rsolve( valueType xb[] ) const;
    void RsolveTransposed( valueType xb[] ) const;

  public:

    using LinearSystemSolver<T>::factorize;

    TridiagonalQR()
    : m_allocReals("allocReals")
    , m_nRC(0)
    {}

    virtual
    ~TridiagonalQR() LAPACK_WRAPPER_OVERRIDE {
      m_allocReals.free();
    }

    void
    factorize(
      char const      who[],
      integer         N,
      valueType const L[],
      valueType const D[],
      valueType const U[]
    );

    bool
    factorize(
      integer         N,
      valueType const L[],
      valueType const D[],
      valueType const U[]
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

    virtual
    bool
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    bool
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    bool
    solve(
      integer   nrhs,
      valueType xb[],
      integer   ldXB
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    bool
    t_solve(
      integer   nrhs,
      valueType xb[],
      integer   ldXB
    ) const LAPACK_WRAPPER_OVERRIDE;

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
      valueType       alpha,
      valueType const L[],
      valueType const D[],
      valueType const U[],
      valueType const x[],
      valueType       beta,
      valueType       y[]
    ) const;

  };

}

///
/// eof: trid.hxx
///
