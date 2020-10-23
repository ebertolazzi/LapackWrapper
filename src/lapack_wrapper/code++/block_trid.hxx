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
    typedef T valueType;

  private:

    Malloc<valueType>  m_allocReals;
    Malloc<integer>    m_allocIntegers;
    Malloc<valueType*> m_allocRpointers;
    Malloc<integer*>   m_allocIpointers;

    integer      m_nBlocks;
    integer      m_nnz;
    valueType ** m_D_blocks;
    valueType ** m_L_blocks;
    valueType *  m_Work;
    integer   ** m_B_permutation;
    integer   *  m_row_blocks;
    bool         m_is_factorized;

  public:

    using LinearSystemSolver<T>::factorize;

    BlockTridiagonalSymmetic()
    : m_allocReals("BlockTridiagonalSymmetic-allocReals")
    , m_allocIntegers("BlockTridiagonalSymmetic-allocIntegers")
    , m_allocRpointers("BlockTridiagonalSymmetic-allocRpointers")
    , m_allocIpointers("BlockTridiagonalSymmetic-allocIpointers")
    , m_nBlocks(0)
    , m_nnz(0)
    , m_is_factorized(false)
    {}

    virtual
    ~BlockTridiagonalSymmetic() UTILS_OVERRIDE {
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

    integer numBlocks() const { return m_nBlocks; }

    integer DnumRows( integer n ) const { return m_row_blocks[n+1] - m_row_blocks[n]; }
    integer DnumCols( integer n ) const { return m_row_blocks[n+1] - m_row_blocks[n]; }

    integer LnumRows( integer n ) const { return DnumRows(n+1); }
    integer LnumCols( integer n ) const { return DnumCols(n); }

    void
    setD(
      integer         n,
      valueType const data[],
      integer         ldData,
      bool            transposed=false
    );

    void
    setL(
      integer         n,
      valueType const data[],
      integer         ldData,
      bool            transposed=false
    );

    void
    setD(
      integer         n,
      valueType const data[],
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
      valueType const data[],
      integer         ldData,
      integer         beginRow,
      integer         beginCol,
      integer         nrow,
      integer         ncol,
      bool            transposed=false
    );

    void
    insert( integer i, integer j, valueType v, bool sym );

    valueType const & operator () ( integer ii, integer jj ) const;
    valueType       & operator () ( integer ii, integer jj );

    void factorize( char const who[] );
    bool factorize( );

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
      valueType xb[],
      integer   ldXB
    ) const UTILS_OVERRIDE;

    virtual
    bool
    t_solve(
      integer   nrhs,
      valueType xb[],
      integer   ldXB
    ) const UTILS_OVERRIDE;

  };

}

///
/// eof: block_trid.hxx
///
