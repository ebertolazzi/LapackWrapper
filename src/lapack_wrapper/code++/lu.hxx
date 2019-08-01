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
  class LU : public Factorization<T> {
  public:
    typedef typename Factorization<T>::valueType valueType;

  private:

    valueType * Work;
    integer   * Iwork;
    integer   * i_pivot;

    Malloc<valueType> allocReals;
    Malloc<integer>   allocIntegers;

    void check_ls( char const who[] ) const;

  public:

    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;
    using Factorization<T>::factorize;
    using Factorization<T>::solve;
    using Factorization<T>::t_solve;

    using Factorization<T>::nRow;
    using Factorization<T>::nCol;
    using Factorization<T>::Amat;

    LU();
    virtual ~LU() LAPACK_WRAPPER_OVERRIDE;

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
    allocate( integer NR, integer NC ) LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    factorize( char const who[] ) LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) LAPACK_WRAPPER_OVERRIDE;

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
  /*\
  :|:   _    _   _ ____   ___
  :|:  | |  | | | |  _ \ / _ \
  :|:  | |  | | | | |_) | | | |
  :|:  | |__| |_| |  __/| |_| |
  :|:  |_____\___/|_|    \__\_\
  \*/
  template <typename T>
  class LUPQ : public Factorization<T> {
  public:
    typedef typename Factorization<T>::valueType valueType;

  private:

    integer * ipiv;
    integer * jpiv;

    Malloc<valueType> allocReals;
    Malloc<integer>   allocIntegers;

    void check_ls( char const who[] ) const;

  public:

    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;
    using Factorization<T>::factorize;
    using Factorization<T>::solve;
    using Factorization<T>::t_solve;

    using Factorization<T>::nRow;
    using Factorization<T>::nCol;
    using Factorization<T>::Amat;

    LUPQ();
    virtual ~LUPQ() LAPACK_WRAPPER_OVERRIDE;

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
    factorize( char const who[] ) LAPACK_WRAPPER_OVERRIDE;

    virtual
    void
    factorize(
      char const      who[],
      integer         NR,
      integer         NC,
      valueType const A[],
      integer         LDA
    ) LAPACK_WRAPPER_OVERRIDE;

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

}

///
/// eof: lu.hxx
///
