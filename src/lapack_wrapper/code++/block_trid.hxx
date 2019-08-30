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

    Malloc<valueType>  allocReals;
    Malloc<integer>    allocIntegers;
    Malloc<valueType*> allocRpointers;
    Malloc<integer*>   allocIpointers;

    integer      nBlocks, nnz;
    valueType ** D_blocks;
    valueType ** L_blocks;
    valueType *  Work;
    integer   ** B_permutation;
    integer   *  row_blocks;
    bool         is_factorized;

  public:

    BlockTridiagonalSymmetic()
    : allocReals("BlockTridiagonalSymmetic-allocReals")
    , allocIntegers("BlockTridiagonalSymmetic-allocIntegers")
    , allocRpointers("BlockTridiagonalSymmetic-allocRpointers")
    , allocIpointers("BlockTridiagonalSymmetic-allocIpointers")
    , nBlocks(0)
    , nnz(0)
    , is_factorized(false)
    {}

    virtual
    ~BlockTridiagonalSymmetic() LAPACK_WRAPPER_OVERRIDE {
      allocReals.free();
      allocIntegers.free();
      allocRpointers.free();
      allocIpointers.free();
    }

    void
    setup( integer nblks, integer const rBlocks[] );

    void
    setup( integer nblks, integer const block_size );

    void
    zero();

    integer numBlocks() const { return nBlocks; }

    integer DnumRows( integer n ) const { return row_blocks[n+1] - row_blocks[n]; }
    integer DnumCols( integer n ) const { return row_blocks[n+1] - row_blocks[n]; }

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

    void
    factorize( char const who[] );

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

  };

}

///
/// eof: block_trid.hxx
///
