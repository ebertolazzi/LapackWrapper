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
/// file: ls.hxx
///

namespace lapack_wrapper {

  //============================================================================
  /*\
  :|:   _     ____ ____
  :|:  | |   / ___/ ___|
  :|:  | |   \___ \___ \
  :|:  | |___ ___) |__) |
  :|:  |_____|____/____/
  \*/

  template <typename T>
  class LSS_no_alloc : public LinearSystemSolver<T> {
  public:
    typedef typename LinearSystemSolver<T>::valueType valueType;

  protected:

    mutable Malloc<valueType> allocWork;
    mutable valueType       * Work;
    mutable integer           Lwork;

    integer nRows;
    integer nCols;

            valueType * Amat;
    mutable valueType * sigma;
    mutable valueType * AmatWork;
    mutable integer     rank;

    valueType rcond;

  public:

    using LinearSystemSolver<T>::factorize;
    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;

    explicit
    LSS_no_alloc()
    : LinearSystemSolver<T>()
    , allocWork("LSS_no_alloc-allocReals")
    , Work(nullptr)
    , Lwork(0)
    , nRows(0)
    , nCols(0)
    , Amat(nullptr)
    , sigma(nullptr)
    , AmatWork(nullptr)
    , rank(0)
    , rcond(-1)
    {}

    virtual
    ~LSS_no_alloc() LAPACK_WRAPPER_OVERRIDE
    { allocWork.free(); }

    integer getL( integer NR, integer NC, integer nrhs ) const;

    void
    no_allocate(
      integer     NR,
      integer     NC,
      valueType * _Amat,
      valueType * _sigma,
      valueType * _AmatWork
    ) {
      this->nRows    = NR;
      this->nCols    = NC;
      this->Amat     = _Amat;
      this->sigma    = _sigma;
      this->AmatWork = _AmatWork;
    }

    /*!
     *  Do SVD factorization of a rectangular matrix
     *  \param who  string with the name of the calling routine
     *  \param A     pointer to the matrix
     *  \param LDA Leading dimension of the matrix
     */
    void
    factorize(
      char const      who[],
      valueType const A[],
      integer         LDA
    );

    void      setRcond( valueType r )     { this->rcond = r; }
    integer   getRank()             const { return this->rank; }
    valueType getSigma( integer i ) const { return this->sigma[i]; }

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
    solve( integer nrhs, valueType B[], integer ldB ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve( integer nrhs, valueType B[], integer ldB ) const LAPACK_WRAPPER_OVERRIDE;

  };

  //============================================================================
  //============================================================================
  //============================================================================

  template <typename T>
  class LSS : public LSS_no_alloc<T> {
  public:
    typedef typename LSS_no_alloc<T>::valueType valueType;

  protected:

    Malloc<valueType> allocReals;

  public:

    using LSS_no_alloc<T>::solve;
    using LSS_no_alloc<T>::t_solve;
    using LSS_no_alloc<T>::setRcond;
    using LSS_no_alloc<T>::getRank;
    using LSS_no_alloc<T>::getSigma;
    using LSS_no_alloc<T>::factorize;

    explicit
    LSS()
    : LSS_no_alloc<T>()
    , allocReals("LSS-allocReals")
    {}

    virtual
    ~LSS() LAPACK_WRAPPER_OVERRIDE
    { allocReals.free(); }

    void
    allocate( integer NR, integer NC );

    /*!
     *  Do SVD factorization of a rectangular matrix
     *  \param who  string with the name of the calling routine
     *  \param NR   number of rows of the matrix
     *  \param NC   number of columns of the matrix
     *  \param A     pointer to the matrix
     *  \param LDA Leading dimension of the matrix
     */
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
  :|:  _     ______   __
  :|: | |   / ___\ \ / /
  :|: | |   \___ \\ V /
  :|: | |___ ___) || |
  :|: |_____|____/ |_|
  \*/
  template <typename T>
  class LSY_no_alloc : public LinearSystemSolver<T> {
  public:
    typedef typename LinearSystemSolver<T>::valueType valueType;

  protected:

    mutable Malloc<valueType> allocWork;
    mutable valueType       * Work;
    mutable integer           Lwork;

    integer nRows;
    integer nCols;

    valueType * Amat;

    mutable valueType * AmatWork;
    mutable integer   * jpvt;
    mutable integer     rank;

    valueType rcond;

  public:

    using LinearSystemSolver<T>::factorize;
    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;

    explicit
    LSY_no_alloc()
    : LinearSystemSolver<T>()
    , allocWork("LSY_no_alloc")
    , Work(nullptr)
    , Lwork(0)
    , nRows(0)
    , nCols(0)
    , Amat(nullptr)
    , AmatWork(nullptr)
    , jpvt(nullptr)
    , rank(0)
    , rcond(-1)
    {}

    virtual
    ~LSY_no_alloc() LAPACK_WRAPPER_OVERRIDE
    {}

    integer getL( integer NR, integer NC, integer nrhs ) const;

    void
    no_allocate(
      integer     NR,
      integer     NC,
      valueType * _Amat,
      valueType * _AmatWork,
      integer   * _jpvt
    ) {
      this->nRows    = NR;
      this->nCols    = NC;
      this->Amat     = _Amat;
      this->AmatWork = _AmatWork;
      this->jpvt     = _jpvt;
    }

    /*!
     *  Do SVD factorization of a rectangular matrix
     *  \param A   pointer to the matrix
     *  \param LDA Leading dimension of the matrix
     */
    void
    factorize(
      char const      who[],
      valueType const A[],
      integer         LDA
    );

    void
    setRcond( valueType r )
    { rcond = r; }

    integer
    getRank() const
    { return rank; }

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
    solve( integer nrhs, valueType B[], integer ldB ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    t_solve( integer nrhs, valueType B[], integer ldB ) const LAPACK_WRAPPER_OVERRIDE;

  };

  //============================================================================
  //============================================================================
  //============================================================================

  template <typename T>
  class LSY : public LSY_no_alloc<T> {
  public:
    typedef typename LSY_no_alloc<T>::valueType valueType;

  protected:

    Malloc<valueType> allocReals;
    Malloc<integer>   allocInts;

  public:

    using LSY_no_alloc<T>::factorize;
    using LSY_no_alloc<T>::solve;
    using LSY_no_alloc<T>::t_solve;
    using LSY_no_alloc<T>::setRcond;
    using LSY_no_alloc<T>::getRank;

    explicit
    LSY()
    : LSY_no_alloc<T>()
    , allocReals("LSY-allocReals")
    , allocInts("LSY-allocInts")
    {}

    virtual
    ~LSY() LAPACK_WRAPPER_OVERRIDE
    { allocReals.free(); allocInts.free(); }

    void
    allocate( integer NR, integer NC );

    /*!
     *  Do SVD factorization of a rectangular matrix
     *  \param NR  number of rows of the matrix
     *  \param NC  number of columns of the matrix
     *  \param A   pointer to the matrix
     *  \param LDA Leading dimension of the matrix
     */
    void
    factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) LAPACK_WRAPPER_OVERRIDE;

  };

}

///
/// eof: ls.hxx
///
