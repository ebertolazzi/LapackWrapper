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

    mutable Malloc<valueType> m_allocWork;
    mutable valueType       * m_Work;
    mutable integer           m_Lwork;

    integer m_nrows;
    integer m_ncols;

            valueType * m_Amat;
    mutable valueType * m_sigma;
    mutable valueType * m_AmatWork;
    mutable integer     m_rank;

    valueType m_rcond;

  public:

    using LinearSystemSolver<T>::factorize;
    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;

    explicit
    LSS_no_alloc()
    : LinearSystemSolver<T>()
    , m_allocWork("LSS_no_alloc-allocReals")
    , m_Work(nullptr)
    , m_Lwork(0)
    , m_nrows(0)
    , m_ncols(0)
    , m_Amat(nullptr)
    , m_sigma(nullptr)
    , m_AmatWork(nullptr)
    , m_rank(0)
    , m_rcond(-1)
    {}

    ~LSS_no_alloc() override
    { m_allocWork.free(); }

    integer getL( integer NR, integer NC, integer nrhs ) const;

    void
    no_allocate(
      integer     NR,
      integer     NC,
      integer     Lwork,
      valueType * Work
    );

    //!
    //! Do SVD factorization of a rectangular matrix.
    //!
    //! \param who  string with the name of the calling routine
    //! \param A     pointer to the matrix
    //! \param LDA Leading dimension of the matrix
    //!
    void
    factorize_nodim(
      char const      who[],
      valueType const A[],
      integer         LDA
    );

    bool
    factorize_nodim( valueType const A[], integer LDA );

    void      setRcond( valueType r )     { m_rcond = r; }
    integer   getRank()             const { return m_rank; }
    valueType getSigma( integer i ) const { return m_sigma[i]; }

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    bool solve( valueType xb[] ) const override;
    bool t_solve( valueType xb[] ) const override;
    bool solve( integer nrhs, valueType B[], integer ldB ) const override;
    bool t_solve( integer nrhs, valueType B[], integer ldB ) const override;

  };

  //============================================================================
  //============================================================================
  //============================================================================

  template <typename T>
  class LSS : public LSS_no_alloc<T> {
  public:
    typedef typename LSS_no_alloc<T>::valueType valueType;

  protected:

    Malloc<valueType> m_allocReals;

  public:

    using LSS_no_alloc<T>::m_allocWork;
    using LSS_no_alloc<T>::m_Work;
    using LSS_no_alloc<T>::m_Lwork;
    using LSS_no_alloc<T>::m_nrows;
    using LSS_no_alloc<T>::m_ncols;
    using LSS_no_alloc<T>::m_Amat;
    using LSS_no_alloc<T>::m_sigma;
    using LSS_no_alloc<T>::m_AmatWork;
    using LSS_no_alloc<T>::m_rank;
    using LSS_no_alloc<T>::m_rcond;

    using LSS_no_alloc<T>::solve;
    using LSS_no_alloc<T>::t_solve;
    using LSS_no_alloc<T>::setRcond;
    using LSS_no_alloc<T>::getRank;
    using LSS_no_alloc<T>::getSigma;
    using LSS_no_alloc<T>::factorize;

    explicit
    LSS()
    : LSS_no_alloc<T>()
    , m_allocReals("LSS-allocReals")
    {}

    ~LSS() override
    { m_allocReals.free(); }

    void
    allocate( integer NR, integer NC );

    //!
    //! Do SVD factorization of a rectangular matrix.
    //!
    //! \param who  string with the name of the calling routine
    //! \param NR   number of rows of the matrix
    //! \param NC   number of columns of the matrix
    //! \param A     pointer to the matrix
    //! \param LDA Leading dimension of the matrix
    //!
    void
    factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) override {
      this->allocate( NR, NC );
      this->factorize_nodim( who, A, LDA );
    }

    bool
    factorize(
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) override {
      this->allocate( NR, NC );
      return this->factorize_nodim( A, LDA );
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

    mutable Malloc<valueType> m_allocWork;
    mutable valueType       * m_Work;
    mutable integer           m_Lwork;

    integer m_nrows;
    integer m_ncols;

    valueType * m_Amat;

    mutable valueType * m_AmatWork;
    mutable integer   * m_jpvt;
    mutable integer     m_rank;

    valueType m_rcond;

  public:

    using LinearSystemSolver<T>::factorize;
    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;

    explicit
    LSY_no_alloc()
    : LinearSystemSolver<T>()
    , m_allocWork("LSY_no_alloc")
    , m_Work(nullptr)
    , m_Lwork(0)
    , m_nrows(0)
    , m_ncols(0)
    , m_Amat(nullptr)
    , m_AmatWork(nullptr)
    , m_jpvt(nullptr)
    , m_rank(0)
    , m_rcond(-1)
    {}

    ~LSY_no_alloc() override
    {}

    integer getL( integer NR, integer NC, integer nrhs ) const;

    void
    no_allocate(
      integer     NR,
      integer     NC,
      integer     Lwork,
      valueType * Work,
      integer     Liwork,
      integer   * iWork
    );

    //!
    //! Do SVD factorization of a rectangular matrix.
    //!
    //! \param who string used in error message
    //! \param A   pointer to the matrix
    //! \param LDA Leading dimension of the matrix
    //!
    void
    factorize_nodim(
      char const      who[],
      valueType const A[],
      integer         LDA
    );

    bool
    factorize_nodim( valueType const A[], integer LDA );

    void
    setRcond( valueType r )
    { m_rcond = r; }

    integer
    getRank() const
    { return m_rank; }

    /*\
    :|:         _      _               _
    :|:  __   _(_)_ __| |_ _   _  __ _| |___
    :|:  \ \ / / | '__| __| | | |/ _` | / __|
    :|:   \ V /| | |  | |_| |_| | (_| | \__ \
    :|:    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    bool solve( valueType xb[] ) const override;
    bool t_solve( valueType xb[] ) const override;
    bool solve( integer nrhs, valueType B[], integer ldB ) const override;
    bool t_solve( integer nrhs, valueType B[], integer ldB ) const override;

  };

  //============================================================================
  //============================================================================
  //============================================================================

  template <typename T>
  class LSY : public LSY_no_alloc<T> {
  public:
    typedef typename LSY_no_alloc<T>::valueType valueType;

  protected:

    Malloc<valueType> m_allocReals;
    Malloc<integer>   m_allocInts;

  public:

    using LSY_no_alloc<T>::m_allocWork;
    using LSY_no_alloc<T>::m_Work;
    using LSY_no_alloc<T>::m_Lwork;
    using LSY_no_alloc<T>::m_nrows;
    using LSY_no_alloc<T>::m_ncols;
    using LSY_no_alloc<T>::m_Amat;
    using LSY_no_alloc<T>::m_AmatWork;
    using LSY_no_alloc<T>::m_rank;
    using LSY_no_alloc<T>::m_rcond;

    using LSY_no_alloc<T>::factorize;
    using LSY_no_alloc<T>::solve;
    using LSY_no_alloc<T>::t_solve;
    using LSY_no_alloc<T>::setRcond;
    using LSY_no_alloc<T>::getRank;

    explicit
    LSY()
    : LSY_no_alloc<T>()
    , m_allocReals("LSY-allocReals")
    , m_allocInts("LSY-allocInts")
    {}

    ~LSY() override
    { m_allocReals.free(); m_allocInts.free(); }

    void
    allocate( integer NR, integer NC );

    //!
    //! Do SVD factorization of a rectangular matrix.
    //!
    //! \param who string used in error message
    //! \param NR  number of rows of the matrix
    //! \param NC  number of columns of the matrix
    //! \param A   pointer to the matrix
    //! \param LDA Leading dimension of the matrix
    //!
    void
    factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) override;

    bool
    factorize(
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) override;

  };

}

///
/// eof: ls.hxx
///
