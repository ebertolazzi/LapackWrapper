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
/// file: band.hxx
///

namespace lapack_wrapper {

  //============================================================================
  /*\
  :|:   ____                  _          _ _    _   _
  :|:  | __ )  __ _ _ __   __| | ___  __| | |  | | | |
  :|:  |  _ \ / _` | '_ \ / _` |/ _ \/ _` | |  | | | |
  :|:  | |_) | (_| | | | | (_| |  __/ (_| | |__| |_| |
  :|:  |____/ \__,_|_| |_|\__,_|\___|\__,_|_____\___/
  \*/
  //! base class for linear systema solver
  template <typename T>
  class BandedLU : public LinearSystemSolver<T> {
  public:
    typedef T valueType;

    Malloc<valueType> allocReals;
    Malloc<integer>   allocIntegers;

    integer     m, n, nL, nU, ldAB;
    integer   * ipiv;
    valueType * AB;

    bool is_factorized;

  public:

    BandedLU();
    virtual ~BandedLU() LAPACK_WRAPPER_OVERRIDE;

    void
    setup(
      integer M,  // number of rows
      integer N,  // number of columns
      integer nL, // number of lower diagonal
      integer nU  // number of upper diagonal
    );

    integer
    iaddr( integer i, integer j ) const {
      integer d = (i-j+nL+nU);
      return d+j*ldAB;
    }

    void
    check( integer i, integer j ) const;

    valueType const &
    operator () ( integer i, integer j ) const
    { return AB[iaddr(i,j)]; }

    valueType &
    operator () ( integer i, integer j )
    { return AB[iaddr(i,j)]; }

    void
    insert( integer i, integer j, valueType v, bool sym );

    void zero();

    void
    load_block(
      integer         nr,
      integer         nc,
      valueType const B[],
      integer         ldB,
      integer         irow,
      integer         icol
    );

    // do internal factorization, to be executed (only once) before to call solve or t_solve
    void factorize( char const who[] );

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

    /*\
    :|:     _
    :|:    / \  _   ___  __
    :|:   / _ \| | | \ \/ /
    :|:  / ___ \ |_| |>  <
    :|: /_/   \_\__,_/_/\_\
    :|:
    \*/

    // y <- beta*y + alpha*A*x
    void
    aAxpy(
      valueType       alpha,
      valueType const x[],
      valueType       y[]
    ) const;

    void
    dump( ostream_type & stream ) const;

  };

  //============================================================================
  /*\
  :|:   ____                  _          _ ____  ____  ____
  :|:  | __ )  __ _ _ __   __| | ___  __| / ___||  _ \|  _ \
  :|:  |  _ \ / _` | '_ \ / _` |/ _ \/ _` \___ \| |_) | | | |
  :|:  | |_) | (_| | | | | (_| |  __/ (_| |___) |  __/| |_| |
  :|:  |____/ \__,_|_| |_|\__,_|\___|\__,_|____/|_|   |____/
  \*/
  //! base class for linear systema solver
  template <typename T>
  class BandedSPD : public LinearSystemSolver<T> {
  public:
    typedef T valueType;

    Malloc<valueType> allocReals;

    integer     n, nD, ldAB;
    valueType * AB;
    ULselect    UPLO;
    bool is_factorized;

  public:

    BandedSPD();

    virtual
    ~BandedSPD() LAPACK_WRAPPER_OVERRIDE;

    void
    setup(
      ULselect UPLO,
      integer  N,    // number of rows and columns
      integer  nD    // number of upper diagonal
    );

    valueType const &
    operator () ( integer i, integer j ) const
    { return AB[i+j*ldAB]; }

    valueType &
    operator () ( integer i, integer j )
    { return AB[i+j*ldAB]; }

    void
    insert( integer i, integer j, valueType v, bool sym );

    void zero();

    // do internal fatcorization, to be executed (only once) before to call solve or t_solve
    void factorize( char const who[] );

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

    /*\
    :|:     _
    :|:    / \  _   ___  __
    :|:   / _ \| | | \ \/ /
    :|:  / ___ \ |_| |>  <
    :|: /_/   \_\__,_/_/\_\
    :|:
    \*/

    /* not yet available
    // y <- beta*y + alpha*A*x
    void
    aAxpy( valueType       alpha,
           valueType const x[],
           valueType       y[] ) const;

    void
    dump( ostream_type & stream ) const;
    */
  };

}

///
/// eof: band.hxx
///
