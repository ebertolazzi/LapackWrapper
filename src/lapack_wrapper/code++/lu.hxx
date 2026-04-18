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
/// file: lu.hxx
///

/*!
 * \file lu.hxx
 * \brief Dense LU-based direct solvers.
 *
 * The classes declared here provide LU factorizations for general rectangular
 * and square dense matrices, both with and without complete pivoting, with
 * interfaces that support preallocated workspaces or internally managed memory.
 */

namespace lapack_wrapper {

  //============================================================================
  /*\
  :|:   _    _   _
  :|:  | |  | | | |
  :|:  | |  | | | |
  :|:  | |__| |_| |
  :|:  |_____\___/
  \*/

  /*!
   * \brief LU factorization wrapper using user-provided work buffers.
   *
   * The class stores the factorized matrix together with the row pivot vector
   * in external memory supplied through `no_allocate()`. It is useful when the
   * caller wants full control over allocations.
   */
  template <typename T>
  class LU_no_alloc : public LinearSystemSolver<T> {
  public:
    using real_type = typename LinearSystemSolver<T>::real_type;

  protected:

    integer     m_nrows{0};
    integer     m_ncols{0};

    real_type * m_Afactorized{nullptr};
    real_type * m_Work{nullptr};
    integer   * m_Iwork{nullptr};
    integer   * m_i_pivot{nullptr};

    void check_ls( string_view who ) const;

  public:

    using LinearSystemSolver<T>::factorize;
    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;

    LU_no_alloc() : LinearSystemSolver<T>() {}

    virtual ~LU_no_alloc() override = default;

    void
    no_allocate(
      integer     NR,
      integer     NC,
      real_type * _Afactorized,
      integer   * _i_pivot,
      real_type * _Work  = nullptr,
      integer   * _Iwork = nullptr
    ) {
      m_nrows       = NR;
      m_ncols       = NC;
      m_Afactorized = _Afactorized;
      m_i_pivot     = _i_pivot;
      // only for condition number
      m_Work        = _Work;
      m_Iwork       = _Iwork;
    }

    void
    factorize_nodim( string_view who, real_type const A[], integer LDA );

    bool
    factorize_nodim( real_type const A[], integer LDA );

    real_type cond1( real_type norm1 ) const;
    real_type condInf( real_type normInf ) const;

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

  };

  //============================================================================

  /*!
   * \brief Owning LU solver with internally managed storage.
   *
   * `LU` allocates the memory required by `LU_no_alloc` and exposes the same
   * dense factorization and solve interface.
   */
  template <typename T>
  class LU : public LU_no_alloc<T> {
  public:
    using real_type = typename LU_no_alloc<T>::real_type;

  private:

    Malloc<real_type> m_allocReals{"LUPQ-allocReals"};
    Malloc<integer>   m_allocIntegers{"LUPQ-allocIntegers"};

  public:

    using LU_no_alloc<T>::m_nrows;
    using LU_no_alloc<T>::m_ncols;
    using LU_no_alloc<T>::m_Afactorized;
    using LU_no_alloc<T>::m_Work;
    using LU_no_alloc<T>::m_Iwork;
    using LU_no_alloc<T>::m_i_pivot;

    using LU_no_alloc<T>::factorize;
    using LU_no_alloc<T>::solve;
    using LU_no_alloc<T>::t_solve;
    using LU_no_alloc<T>::cond1;
    using LU_no_alloc<T>::condInf;

    LU() : LU_no_alloc<T>() {}

    void allocate( integer NR, integer NC );

    void
    factorize(
      string_view     who,
      integer         NR,
      integer         NC,
      real_type const A[],
      integer         LDA
    ) override {
      this->allocate( NR, NC );
      this->factorize_nodim( who, A, LDA );
    }

    bool
    factorize(
      integer         NR,
      integer         NC,
      real_type const A[],
      integer         LDA
    ) override {
      this->allocate( NR, NC );
      return this->factorize_nodim( A, LDA );
    }

  };

  //============================================================================

  /*\
  :|:   _    _   _ ____   ___
  :|:  | |  | | | |  _ \ / _ \
  :|:  | |  | | | | |_) | | | |
  :|:  | |__| |_| |  __/| |_| |
  :|:  |_____\___/|_|    \__\_\
  \*/

  /*!
   * \brief LU factorization with complete pivoting on caller-owned storage.
   *
   * The factorization is intended for square matrices and tracks both row and
   * column permutations, which makes it more robust on ill-conditioned or
   * nearly singular systems than plain partial-pivot LU.
   */
  template <typename T>
  class LUPQ_no_alloc : public LinearSystemSolver<T> {
  public:
    using real_type = typename LinearSystemSolver<T>::real_type;

  protected:

    integer     m_nRC{0};
    real_type * m_Afactorized{nullptr};
    integer   * m_i_piv{nullptr};
    integer   * m_j_piv{nullptr};

  public:

    using LinearSystemSolver<T>::factorize;
    using LinearSystemSolver<T>::solve;
    using LinearSystemSolver<T>::t_solve;

    LUPQ_no_alloc() = default;

    void
    no_allocate(
      integer     NRC,
      integer     Lwork,
      real_type * Work,
      integer     Liwork,
      integer   * iWork
    );

    void
    factorize_nodim(
      string_view     who,
      real_type const A[],
      integer         LDA
    );

    bool
    factorize_nodim( real_type const A[], integer LDA );

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
      real_type B[],
      integer   ldB
    ) const override;

    bool
    t_solve(
      integer   nrhs,
      real_type B[],
      integer   ldB
    ) const override;

  };

  //============================================================================

  /*!
   * \brief Owning complete-pivot LU solver.
   *
   * This class manages the workspace required by `LUPQ_no_alloc` and offers
   * the same factorize/solve API with automatic allocation.
   */
  template <typename T>
  class LUPQ : public LUPQ_no_alloc<T> {
  public:
    using real_type = typename LUPQ_no_alloc<T>::real_type;

  private:

    Malloc<real_type> m_allocReals{""};
    Malloc<integer>   m_allocIntegers{""};

  public:

    using LUPQ_no_alloc<T>::m_nRC;
    using LUPQ_no_alloc<T>::m_Afactorized;
    using LUPQ_no_alloc<T>::m_i_piv;
    using LUPQ_no_alloc<T>::m_j_piv;

    using LUPQ_no_alloc<T>::solve;
    using LUPQ_no_alloc<T>::t_solve;
    using LUPQ_no_alloc<T>::factorize;

    LUPQ() = default;

    ~LUPQ() override {
      m_allocReals.free();
      m_allocIntegers.free();
    }

    void allocate( integer NRC );

    void
    factorize(
      string_view     who,
      integer         NRC,
      real_type const A[],
      integer         LDA
    ) {
      this->allocate( NRC );
      this->factorize_nodim( who, A, LDA );
    }

    void
    factorize(
      string_view     who,
      integer         NR,
      integer         NC,
      real_type const A[],
      integer         LDA
    ) override {
      UTILS_ASSERT( NR == NC, "LUPQ::factorize, non square matrix: {} x {}\n", NR, NC );
      factorize( who, NR, A, LDA );
    }

    bool
    factorize(
      integer         NRC,
      real_type const A[],
      integer         LDA
    ) {
      this->allocate( NRC );
      return this->factorize_nodim( A, LDA );
    }

    bool
    factorize(
      integer         NR,
      integer         NC,
      real_type const A[],
      integer         LDA
    ) override {
      if ( NR != NC ) return false;
      return factorize( NR, A, LDA );
    }

  };

}

///
/// eof: lu.hxx
///
