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
    using real_type = T;

    Malloc<real_type> m_allocReals{"BandedLU_reals"};
    Malloc<integer>   m_allocIntegers{"BandedLU_integers"};

    integer     m_m{0};
    integer     m_n{0};
    integer     m_nL{0};
    integer     m_nU{0};
    integer     m_ldAB{0};
    integer   * m_ipiv{nullptr};
    real_type * m_AB{nullptr};
    bool        m_is_factorized{false};

  public:

    using LinearSystemSolver<T>::factorize;

    BandedLU() = default;
    ~BandedLU() override = default;

    void
    setup(
      integer M,  // number of rows
      integer N,  // number of columns
      integer nL, // number of lower diagonal
      integer nU  // number of upper diagonal
    );

    integer
    iaddr( integer i, integer j ) const {
      integer d{ i - j + m_nL + m_nU };
      return d+j*m_ldAB;
    }

    void
    check( integer i, integer j ) const;

    real_type const &
    operator () ( integer i, integer j ) const
    { return m_AB[iaddr(i,j)]; }

    real_type &
    operator () ( integer i, integer j )
    { return m_AB[iaddr(i,j)]; }

    void
    insert( integer i, integer j, real_type v, bool sym );

    void zero();

    void
    load_block(
      integer         nr,
      integer         nc,
      real_type const B[],
      integer         ldB,
      integer         irow,
      integer         icol
    );

    // do internal factorization, to be executed (only once) before to call solve or t_solve
    void factorize( string_view who );

    bool factorize();

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
      real_type B[],
      integer   ldB
    ) const override;

    void
    solve(
      string_view who,
      integer     nrhs,
      real_type   B[],
      integer     ldB
    ) const override;

    bool
    t_solve(
      integer   nrhs,
      real_type B[],
      integer   ldB
    ) const override;

    void
    t_solve(
      string_view who,
      integer     nrhs,
      real_type   B[],
      integer     ldB
    ) const override;

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
      real_type       alpha,
      real_type const x[],
      real_type       y[]
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
    using real_type = T;

    Malloc<real_type> m_allocReals;

    integer     m_n{0}, m_nD{0}, m_ldAB{0};
    real_type * m_AB{nullptr};
    ULselect    m_UPLO;
    bool        m_is_factorized{false};

  public:

    using LinearSystemSolver<T>::factorize;

    BandedSPD();
    ~BandedSPD() override;

    void
    setup(
      ULselect UPLO,
      integer  N,    // number of rows and columns
      integer  nD    // number of upper diagonal
    );

    real_type const &
    operator () ( integer i, integer j ) const
    { return m_AB[i+j*m_ldAB]; }

    real_type &
    operator () ( integer i, integer j )
    { return m_AB[i+j*m_ldAB]; }

    void
    insert( integer i, integer j, real_type v, bool sym );

    void zero();

    // do internal fatcorization, to be executed (only once) before to call solve or t_solve
    void factorize( string_view who );

    bool factorize( );

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
      real_type B[],
      integer   ldB
    ) const override;

    void
    solve(
      string_view who,
      integer     nrhs,
      real_type   B[],
      integer     ldB
    ) const override;

    bool
    t_solve(
      integer   nrhs,
      real_type B[],
      integer   ldB
    ) const override;

    void
    t_solve(
      string_view who,
      integer     nrhs,
      real_type   B[],
      integer     ldB
    ) const override;

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
    aAxpy( real_type       alpha,
           real_type const x[],
           real_type       y[] ) const;

    void
    dump( ostream_type & stream ) const;
    */
  };

}

///
/// eof: band.hxx
///
