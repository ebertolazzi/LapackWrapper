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
/// file: lu.hxx
///

namespace lapack_wrapper {

  //============================================================================
  /*\
  :|:   _    _   _
  :|:  | |  | | | |
  :|:  | |  | | | |
  :|:  | |__| |_| |
  :|:  |_____\___/
  \*/

  template <typename T>
  class LU_no_alloc : public LinearSystemSolver<T> {
  public:
    typedef typename LinearSystemSolver<T>::valueType valueType;

  protected:

    integer nRows;
    integer nCols;

    valueType * Afactorized;
    valueType * Work;
    integer   * Iwork;
    integer   * i_pivot;

    void check_ls( char const who[] ) const;

  public:

    using LinearSystemSolver<T>::factorize;
    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;

    LU_no_alloc();
    virtual ~LU_no_alloc() LAPACK_WRAPPER_OVERRIDE {}

    void
    no_allocate(
      integer     NR,
      integer     NC,
      valueType * _Afactorized,
      integer   * _i_pivot,
      valueType * _Work  = nullptr,
      integer   * _Iwork = nullptr
    ) {
      this->nRows       = NR;
      this->nCols       = NC;
      this->Afactorized = _Afactorized;
      this->i_pivot     = _i_pivot;
      // only for condition number
      this->Work        = _Work;
      this->Iwork       = _Iwork;
    }

    void
    factorize(
      char const      who[],
      valueType const A[],
      integer         LDA
    );

    valueType cond1( valueType norm1 ) const;
    valueType condInf( valueType normInf ) const;

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
      valueType B[],
      integer   ldB
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const LAPACK_WRAPPER_OVERRIDE;

  };

  //============================================================================

  template <typename T>
  class LU : public LU_no_alloc<T> {
  public:
    typedef typename LU_no_alloc<T>::valueType valueType;

  private:

    Malloc<valueType> allocReals;
    Malloc<integer>   allocIntegers;

  public:

    using LU_no_alloc<T>::factorize;
    using LU_no_alloc<T>::solve;
    using LU_no_alloc<T>::t_solve;
    using LU_no_alloc<T>::cond1;
    using LU_no_alloc<T>::condInf;

    LU();
    virtual ~LU() LAPACK_WRAPPER_OVERRIDE;

    void allocate( integer NR, integer NC );

    void
    factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) LAPACK_WRAPPER_OVERRIDE {
      this->allocate( NR, NC );
      this->factorize( who, A, LDA );
    }

  };

  //============================================================================

  /*\
  :|:   _    _   _ ____   ___
  :|:  | |  | | | |  _ \ / _ \
  :|:  | |  | | | | |_) | | | |
  :|:  | |__| |_| |  __/| |_| |
  :|:  |_____\___/|_|    \__\_\
  \*/

  template <typename T>
  class LUPQ_no_alloc : public LinearSystemSolver<T> {
  public:
    typedef typename LinearSystemSolver<T>::valueType valueType;

  protected:

    integer     nRC;
    valueType * Afactorized;
    integer   * i_piv;
    integer   * j_piv;

  public:

    using LinearSystemSolver<T>::factorize;
    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;

    LUPQ_no_alloc();
    virtual ~LUPQ_no_alloc() LAPACK_WRAPPER_OVERRIDE {}

    void
    no_allocate(
      integer     NRC,
      valueType * _Afactorized,
      integer   * _i_piv,
      integer   * _j_piv
    ) {
      this->nRC         = NRC;
      this->Afactorized = _Afactorized;
      this->i_piv       = _i_piv;
      this->j_piv       = _j_piv;
    }

    void
    factorize(
      char const      who[],
      valueType const A[],
      integer         LDA
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
      valueType B[],
      integer   ldB
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const LAPACK_WRAPPER_OVERRIDE;

  };

  //============================================================================

  template <typename T>
  class LUPQ : public LUPQ_no_alloc<T> {
  public:
    typedef typename LUPQ_no_alloc<T>::valueType valueType;

  private:

    Malloc<valueType> allocReals;
    Malloc<integer>   allocIntegers;

  public:

    using LUPQ_no_alloc<T>::solve;
    using LUPQ_no_alloc<T>::t_solve;
    using LUPQ_no_alloc<T>::factorize;

    LUPQ();
    virtual ~LUPQ() LAPACK_WRAPPER_OVERRIDE;

    void allocate( integer NRC );

    void
    factorize(
      char const      who[],
      integer         NRC,
      valueType const A[],
      integer         LDA
    ) {
      this->allocate( NRC );
      this->factorize( who, A, LDA );
    }

    void
    factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) LAPACK_WRAPPER_OVERRIDE {
      LW_ASSERT(
        NR == NC, "LUPQ::factorize, non square matrix: {} x {}", NR, NC
      );
      factorize( who, NR, A, LDA );
    }

  };

}

///
/// eof: lu.hxx
///
