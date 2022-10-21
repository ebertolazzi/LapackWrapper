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
/// file: block_trid.hxx
///

namespace lapack_wrapper {

  //============================================================================
  /*\
  :|:  ___ _         _     _____    _    _ _                         _
  :|: | _ ) |___  __| |__ |_   _| _(_)__| (_)__ _ __ _ ___ _ _  __ _| |
  :|: | _ \ / _ \/ _| / /   | || '_| / _` | / _` / _` / _ \ ' \/ _` | |
  :|: |___/_\___/\__|_\_\   |_||_| |_\__,_|_\__,_\__, \___/_||_\__,_|_|
  :|:                                            |___/
  :|:  ___                _       _
  :|: / __|_  _ _ __  ___| |_ _ _(_)__
  :|: \__ \ || | '  \/ -_)  _| '_| / _|
  :|: |___/\_, |_|_|_\___|\__|_| |_\__|
  :|:      |__/
  \*/

  template <typename T>
  class BlockTridiagonalSymmetic : public LinearSystemSolver<T> {
  public:
    using real_type = T;

  private:

    Malloc<real_type>  m_allocReals;
    Malloc<integer>    m_allocIntegers;
    Malloc<real_type*> m_allocRpointers;
    Malloc<integer*>   m_allocIpointers;

    integer      m_nBlocks{0};
    integer      m_nnz{0};
    real_type ** m_D_blocks{nullptr};
    real_type ** m_L_blocks{nullptr};
    real_type *  m_Work{nullptr};
    integer   ** m_B_permutation{nullptr};
    integer   *  m_row_blocks{nullptr};
    bool         m_is_factorized{false};

  public:

    using LinearSystemSolver<T>::factorize;

    BlockTridiagonalSymmetic()
    : m_allocReals("BlockTridiagonalSymmetic-allocReals")
    , m_allocIntegers("BlockTridiagonalSymmetic-allocIntegers")
    , m_allocRpointers("BlockTridiagonalSymmetic-allocRpointers")
    , m_allocIpointers("BlockTridiagonalSymmetic-allocIpointers")
    {}

    ~BlockTridiagonalSymmetic() override {
      m_allocReals.free();
      m_allocIntegers.free();
      m_allocRpointers.free();
      m_allocIpointers.free();
    }

    void
    setup( integer nblks, integer const rBlocks[] );

    void
    setup( integer nblks, integer const block_size );

    void
    zero();

    integer num_blocks() const { return m_nBlocks; }

    integer D_nrows( integer n ) const { return m_row_blocks[n+1] - m_row_blocks[n]; }
    integer D_ncols( integer n ) const { return m_row_blocks[n+1] - m_row_blocks[n]; }

    integer L_nrows( integer n ) const { return D_nrows(n+1); }
    integer L_ncols( integer n ) const { return D_ncols(n); }

    // ALIAS
    #ifdef LAPACK_WRAPPER_USE_ALIAS
    integer numBlocks() const { return m_nBlocks; }

    integer DnumRows( integer n ) const { return D_nrows(n); }
    integer DnumCols( integer n ) const { return D_ncols(n); }

    integer LnumRows( integer n ) const { return L_nrows(n); }
    integer LnumCols( integer n ) const { return L_ncols(n); }
    #endif

    void
    setD(
      integer         n,
      real_type const data[],
      integer         ldData,
      bool            transposed=false
    );

    void
    setL(
      integer         n,
      real_type const data[],
      integer         ldData,
      bool            transposed=false
    );

    void
    setD(
      integer         n,
      real_type const data[],
      integer         ldData,
      integer         beginRow,
      integer         beginCol,
      integer         nrow,
      integer         ncol,
      bool            transposed=false
    );

    void
    setL(
      integer         n,
      real_type const data[],
      integer         ldData,
      integer         beginRow,
      integer         beginCol,
      integer         nrow,
      integer         ncol,
      bool            transposed=false
    );

    void
    insert( integer i, integer j, real_type v, bool sym );

    real_type const & operator () ( integer ii, integer jj ) const;
    real_type       & operator () ( integer ii, integer jj );

    void factorize( char const who[] );
    bool factorize( );

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

  };

}

///
/// eof: block_trid.hxx
///
