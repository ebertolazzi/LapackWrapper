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

    integer     m_nrows;
    integer     m_ncols;

    valueType * m_Afactorized;
    valueType * m_Work;
    integer   * m_Iwork;
    integer   * m_i_pivot;

    void check_ls( char const who[] ) const;

  public:

    using LinearSystemSolver<T>::factorize;
    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;

    LU_no_alloc();
    virtual ~LU_no_alloc() UTILS_OVERRIDE {}

    void
    no_allocate(
      integer     NR,
      integer     NC,
      valueType * _Afactorized,
      integer   * _i_pivot,
      valueType * _Work  = nullptr,
      integer   * _Iwork = nullptr
    ) {
      m_nrows       = NR;
      m_ncols       = NC;
      m_Afactorized = _Afactorized;
      m_i_pivot     = _i_pivot;
      // only for condition number
      m_Work        = _Work;
      m_Iwork       = _Iwork;
    }

    void
    factorize( char const who[], valueType const A[], integer LDA );

    bool
    factorize( valueType const A[], integer LDA );

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
    bool
    solve( valueType xb[] ) const UTILS_OVERRIDE;

    virtual
    void
    solve( char const who[], valueType xb[] ) const UTILS_OVERRIDE;

    virtual
    bool
    t_solve( valueType xb[] ) const UTILS_OVERRIDE;

    virtual
    void
    t_solve( char const who[], valueType xb[] ) const UTILS_OVERRIDE;

    virtual
    bool
    solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const UTILS_OVERRIDE;

    virtual
    void
    solve(
      char const who[],
      integer    nrhs,
      valueType  B[],
      integer    ldB
    ) const UTILS_OVERRIDE;

    virtual
    bool
    t_solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const UTILS_OVERRIDE;

    virtual
    void
    t_solve(
      char const who[],
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const UTILS_OVERRIDE;

  };

  //============================================================================

  template <typename T>
  class LU : public LU_no_alloc<T> {
  public:
    typedef typename LU_no_alloc<T>::valueType valueType;

  private:

    Malloc<valueType> m_allocReals;
    Malloc<integer>   m_allocIntegers;

  public:

    using LU_no_alloc<T>::m_nrows;
    using LU_no_alloc<T>::m_ncols;
    using LU_no_alloc<T>::m_Afactorized;
    using LU_no_alloc<T>::m_Work;
    using LU_no_alloc<T>::m_Iwork;
    using LU_no_alloc<T>::m_i_pivot;

    using LU_no_alloc<T>::factorize;
    using LU_no_alloc<T>::solve;
    using LU_no_alloc<T>::t_solve;
    using LU_no_alloc<T>::cond1;
    using LU_no_alloc<T>::condInf;

    LU();
    virtual ~LU() UTILS_OVERRIDE;

    void allocate( integer NR, integer NC );

    void
    factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) UTILS_OVERRIDE {
      this->allocate( NR, NC );
      this->factorize( who, A, LDA );
    }

    bool
    factorize(
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) UTILS_OVERRIDE {
      this->allocate( NR, NC );
      return this->factorize( A, LDA );
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

    integer     m_nRC;
    valueType * m_Afactorized;
    integer   * m_i_piv;
    integer   * m_j_piv;

  public:

    using LinearSystemSolver<T>::factorize;
    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;

    LUPQ_no_alloc();
    virtual ~LUPQ_no_alloc() UTILS_OVERRIDE {}

    void
    no_allocate(
      integer     NRC,
      integer     Lwork,
      valueType * Work,
      integer     Liwork,
      integer   * iWork
    );

    void
    factorize(
      char const      who[],
      valueType const A[],
      integer         LDA
    );

    bool
    factorize( valueType const A[], integer LDA );

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    bool
    solve( valueType xb[] ) const UTILS_OVERRIDE;

    virtual
    bool
    t_solve( valueType xb[] ) const UTILS_OVERRIDE;

    virtual
    bool
    solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const UTILS_OVERRIDE;

    virtual
    bool
    t_solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const UTILS_OVERRIDE;

  };

  //============================================================================

  template <typename T>
  class LUPQ : public LUPQ_no_alloc<T> {
  public:
    typedef typename LUPQ_no_alloc<T>::valueType valueType;

  private:

    Malloc<valueType> m_allocReals;
    Malloc<integer>   m_allocIntegers;

  public:

    using LUPQ_no_alloc<T>::m_nRC;
    using LUPQ_no_alloc<T>::m_Afactorized;
    using LUPQ_no_alloc<T>::m_i_piv;
    using LUPQ_no_alloc<T>::m_j_piv;

    using LUPQ_no_alloc<T>::solve;
    using LUPQ_no_alloc<T>::t_solve;
    using LUPQ_no_alloc<T>::factorize;

    LUPQ();
    virtual ~LUPQ() UTILS_OVERRIDE;

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
    ) UTILS_OVERRIDE {
      UTILS_ASSERT(
        NR == NC, "LUPQ::factorize, non square matrix: {} x {}\n", NR, NC
      );
      factorize( who, NR, A, LDA );
    }

    bool
    factorize(
      integer         NRC,
      valueType const A[],
      integer         LDA
    ) {
      this->allocate( NRC );
      return this->factorize( A, LDA );
    }

    bool
    factorize(
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) UTILS_OVERRIDE {
      if ( NR != NC ) return false;
      return factorize( NR, A, LDA );
    }

  };

}

///
/// eof: lu.hxx
///
