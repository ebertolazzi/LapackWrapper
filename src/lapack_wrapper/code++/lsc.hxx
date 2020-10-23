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
  :|:   _     ____   ____
  :|:  | |   / ___| / ___|
  :|:  | |   \___ \| |
  :|:  | |___ ___) | |___
  :|:  |_____|____/ \____|
  \*/
  /*!  least square constained
  :|:
  :|:  minimize    (1/2) * || A * x - b ||^2
  :|:
  :|:  subject to  B * x = c
  :|:
  :|:  L = (1/2) * || A * x - b ||^2 + lambda * ( B * x - c )
  :|:
  :|:  is equivalent to solve linear system
  :|:
  :|:  / A^T*A  B^T \  /   x    \   / A^T b \
  :|:  |            |  |        | = |       |
  :|:  \ B       0  /  \ lambda /   \   c   /
  :|:
  :|:  is equivalent to solve linear system
  :|:
  :|:  / 0   A^T  B^T \  /   x    \   / A^T b \
  :|:  |              |  |        |   |       |
  :|:  | A   -I    0  |  | omega  | = |   0   |
  :|:  |              |  |        |   |       |
  :|:  \ B    0    0  /  \ lambda /   \   c   /
  :|:
  \*/

  template <typename T>
  class LSC : public LinearSystemSolver<T> {
  public:
    typedef T valueType;

    Malloc<valueType> m_allocReals;
    Malloc<integer>   m_allocIntegers;

    integer     m_NR,  m_NC;
    integer     m_NRA, m_NRB;
    integer   * m_to_rowA;
    integer   * m_to_rowB;
    valueType * m_Amat;
    LU<T>       m_lu;

    mutable Malloc<valueType> m_allocWorks;
    mutable valueType *       m_work;
    mutable valueType *       m_rhs;

  public:

    LSC();
    ~LSC() UTILS_OVERRIDE {}

    void
    factorize(
      char const []      /* who */,
      integer            /* NR  */,
      integer            /* NC  */,
      valueType const [] /* M[] */,
      integer            /* ldM */
    ) UTILS_OVERRIDE {
      UTILS_ERROR0("LSC::factorize, virtual one not defined!");
    }

    bool
    factorize(
      integer            /* NR  */,
      integer            /* NC  */,
      valueType const [] /* M[] */,
      integer            /* ldM */
    ) UTILS_OVERRIDE {
      UTILS_ERROR0("LSC::factorize, virtual one not defined!");
    }

    bool
    factorize(
      integer         NR,
      integer         NC,
      valueType const M[],
      integer         ldM,
      bool const      row_select[] // row selected for minimization
    );

    bool
    solve( valueType xb[] ) const UTILS_OVERRIDE;

    bool
    solve(
      integer   nrhs,
      valueType B[],
      integer   ldB
    ) const UTILS_OVERRIDE;

    bool
    t_solve( valueType [] ) const UTILS_OVERRIDE {
      UTILS_ERROR0("LSC::t_solve not defined!");
    }

  };

}

///
/// eof: lsc.hxx
///
