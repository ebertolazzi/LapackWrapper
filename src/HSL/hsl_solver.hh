/*
 * \file hslsolver.h
 * Header definitions of the class HSLsolver.
 *
 * \author A. Huber and E.Bertolazzi
 * \since  13.11.2018
 */

#ifndef HSL_SOLVER_dot_HH
#define HSL_SOLVER_dot_HH

#include "hsl.hh"
#include <string>
#include <algorithm>

namespace lapack_wrapper {

  template <typename real>
  class HSLsolver {
  protected:
    integer m_nrows{0};
    integer m_ncols{0};
    integer m_nnz{0};

    /// True if the HSL solver is initialized.
    bool m_isInitialized{false};
    /// True if the HSL solver is factorized.
    bool m_isFactorized{false};

    mutable std::string m_last_error{""};

    /**
     * \brief HSLsolver:
     *        Constructor of the class HSLsolver.
     */
    HSLsolver() = default;

  public:

    /**
     * \brief ~HSLsolver:
     *        Virtual destructor of the class HSLsolver.
     */
    virtual ~HSLsolver() = default;

    /*\
    :|:   ___      _    _ _      __  __           _
    :|:  | _ \_  _| |__| (_)__  |  \/  |___ _ __ | |__  ___ _ _ ___
    :|:  |  _/ || | '_ \ | / _| | |\/| / -_) '  \| '_ \/ -_) '_(_-<
    :|:  |_|  \_,_|_.__/_|_\__| |_|  |_\___|_|_|_|_.__/\___|_| /__/
    \*/

    /**
     * \brief init:
     *        Initializes the HSL solver and does a symbolic solve.
     *
     * \param[in] nnz               Number of non zeros.
     * \param[in] N_Row             Row dimension of the matrix.
     * \param[in] N_Col             Col dimension of the matrix.
     * \param[in] i_Row             Rows of the matrix.
     * \param[in] j_Col             Cols of the matrix.
     * \param[in] isFortranIndexing True for Fortran indexing and false for C indexing.
     *
     * \return True if HSL could be initialized.
     *
     */
    virtual
    bool
    init(
      integer       nnz,
      integer       N_Row,
      integer       N_Col,
      integer const i_Row[],
      integer const j_Col[],
      bool          isFortranIndexing,
      bool          isStoredSymmetric
    ) = 0;

    /**
     * \brief factorize:
     *        Factorizes the linear system \f$ Ax = b\f$ for \f$ A\f$ for a sparse matrix.
     *
     * \param[in] ArrayA Values of the sparse matrix.
     *
     * \return True if the linear system \f$ Ax = b\f$ for \f$ A\f$ could be solved.
     *
     */
    virtual
    bool
    factorize( real const ArrayA[] ) = 0;

    /**
     * \brief solve:
     *        Solves the linear system \f$ Ax = b\f$ for \f$ A\f$ for a sparse matrix.
     *
     * \param[in]  nrhs  number of right-hand-side of the linear system.
     * \param[in]  RHS   right-hand-side(s) of the linear system.
     * \param[in]  ldRHS leading dimension of RHS
     * \param[out] X     solution of the linear system.
     * \param[in]  ldX   leading dimension of X
     *
     * \return True if the linear system \f$ Ax = b\f$ for \f$ A\f$ could be solved.
     *
     */
    virtual
    bool
    solve(
      integer    nrhs,
      real const RHS[],
      integer    ldRHS,
      real       X[],
      integer    ldX
    ) const = 0;

    /**
     * \brief solve:
     *        Solves the linear system \f$ Ax = b\f$ for \f$ A\f$ for a sparse matrix.
     *
     * \param[in]  RHS   right-hand-side(s) of the linear system.
     * \param[out] X     solution of the linear system.
     *
     * \return True if the linear system \f$ Ax = b\f$ for \f$ A\f$ could be solved.
     *
     */
    bool
    solve( real const RHS[], real X[] ) const {
      return solve( 1, RHS, m_nrows, X, m_nrows );
    }

    /**
     * \brief solve_transposed:
     *        Solves the linear system \f$ A^\top x = b\f$ for \f$ A\f$ for a sparse matrix.
     *
     * \param[in]  nrhs  number of right-hand-side of the linear system.
     * \param[in]  RHS   right-hand-side(s) of the linear system.
     * \param[in]  ldRHS leading dimension of RHS
     * \param[out] X     solution of the linear system.
     * \param[in]  ldX   leading dimension of X
     *
     * \return True if the linear system \f$ A^\top x = b\f$ for \f$ A\f$ could be solved.
     *
     */
    virtual
    bool
    solve_transposed(
      integer    nrhs,
      real const RHS[],
      integer    ldRHS,
      real       X[],
      integer    ldX
    ) const = 0;

    /**
     * \brief solve_transposed:
     *        Solves the linear system \f$ A^\top x = b\f$ for \f$ A\f$
     *        for a sparse matrix.
     *
     * \param[in]  RHS   right-hand-side(s) of the linear system.
     * \param[out] X     solution of the linear system.
     *
     * \return True if the linear system \f$ A^\top x = b\f$ for \f$ A\f$ could be solved.
     *
     */
    bool
    solve_transposed( real const RHS[], real X[] ) const {
      return solve_transposed( 1, RHS, m_ncols, X, m_ncols );
    }

    /**
     * \brief Checks if matrix stored is initialized
     *
     * \return Return true if matrix stored is initialized
     *
     */
    bool
    is_initialized()
    { return m_isInitialized; }

    /**
     * \brief Checks if matrix stored is facorized
     *
     * \return Return true if matrix stored is factorized
     *
     */
    bool
    is_factorized()
    { return m_isFactorized; }
  };

  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  #pragma clang diagnostic ignored "-Wweak-template-vtables"
  #endif

  extern template class HSLsolver<float>;
  extern template class HSLsolver<double>;

  #ifdef __clang__
  #pragma clang diagnostic pop
  #endif

} // namespace lapack_wrapper

#endif
