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

    Malloc<valueType> allocReals;

    valueType * L;
    valueType * D;
    valueType * WORK;
    integer     nRC;

  public:

    TridiagonalSPD()
    : allocReals("allocReals")
    , nRC(0)
    {}

    virtual
    ~TridiagonalSPD() LAPACK_WRAPPER_OVERRIDE
    { allocReals.free(); }

    valueType cond1( valueType norm1 ) const;

    void
    factorize(
      char const      who[],
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
    void
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    solve(
      integer   nrhs,
      valueType xb[],
      integer   ldXB
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
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

    Malloc<valueType> allocReals;
    Malloc<integer>   allocIntegers;

    valueType * L;
    valueType * D;
    valueType * U;
    valueType * U2;
    valueType * WORK;
    integer   * IPIV;
    integer   * IWORK;

    integer     nRC;

  public:

    TridiagonalLU()
    : allocReals("TridiagonalLU-allocReals")
    , allocIntegers("TridiagonalLU-allocIntegers")
    , nRC(0)
    {}

    virtual
    ~TridiagonalLU() LAPACK_WRAPPER_OVERRIDE {
      allocReals.free();
      allocIntegers.free();
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

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    void
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    solve(
      integer   nrhs,
      valueType xb[],
      integer   ldXB
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
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

    Malloc<valueType> allocReals;

    valueType * C;   // rotazioni givens
    valueType * S;
    valueType * BD;  // band triangular matrix
    valueType * BU;  // band triangular matrix
    valueType * BU2; // band triangular matrix

    valueType   normInfA;
    integer     nRC;

    void Rsolve( valueType xb[] ) const;
    void RsolveTransposed( valueType xb[] ) const;

  public:

    TridiagonalQR()
    : allocReals("allocReals")
    , nRC(0)
    {}

    virtual
    ~TridiagonalQR() LAPACK_WRAPPER_OVERRIDE {
      allocReals.free();
    }

    void
    factorize(
      char const      who[],
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
    void
    solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve( valueType xb[] ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    solve(
      integer   nrhs,
      valueType xb[],
      integer   ldXB
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
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
