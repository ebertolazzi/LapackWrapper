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
  :|:  _     ______   __
  :|: | |   / ___\ \ / /
  :|: | |   \___ \\ V /
  :|: | |___ ___) || |
  :|: |_____|____/ |_|
  \*/
  template <typename T>
  class LSY : public Factorization<T> {
  public:
    typedef typename Factorization<T>::valueType valueType;

  protected:

    Malloc<valueType> allocReals;
    Malloc<integer>   allocInts;

    mutable valueType * Work;
    mutable valueType * AmatWork;
    mutable integer   * jpvt;
    mutable integer     rank;

    valueType rcond;
    integer   Lwork;
    integer   maxNrhs;
    bool      maxNrhs_changed;

  public:

    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;
    using Factorization<T>::factorize;
    using Factorization<T>::solve;
    using Factorization<T>::t_solve;

    using Factorization<T>::nRow;
    using Factorization<T>::nCol;
    using Factorization<T>::Amat;

    explicit
    LSY()
    : Factorization<T>()
    , allocReals("LSY-allocReals")
    , allocInts("LSY-allocInts")
    , Work(nullptr)
    , AmatWork(nullptr)
    , rank(0)
    , rcond(-1)
    , Lwork(0)
    , maxNrhs(1)
    , maxNrhs_changed(true)
    {}

    virtual
    ~LSY() LAPACK_WRAPPER_OVERRIDE
    { allocReals.free(); allocInts.free(); }

    void
    setMaxNrhs( integer mnrhs );

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
    allocate( integer NR, integer NC ) LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    factorize( char const [] ) LAPACK_WRAPPER_OVERRIDE {
      // nothing to do
    }

    /*!
     *  Do SVD factorization of a rectangular matrix
     *  \param NR  number of rows of the matrix
     *  \param NC  number of columns of the matrix
     *  \param A   pointer to the matrix
     *  \param LDA Leading dimension of the matrix
     */
    virtual
    void
    factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) LAPACK_WRAPPER_OVERRIDE {
      allocate( NR, NC );
      integer info = gecopy( nRow, nCol, A, LDA, Amat, nRow );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "LSY::factorize[" << who <<
        "] call lapack_wrapper::gecopy return info = " << info
      );
      factorize( who );
    }

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

}

///
/// eof: ls.hxx
///
