/**
 * \file hslsolver.h
 * Header definitions of the class HSLsolver.
 *
 * \author A. Huber and E.Bertolazzi
 * \since  13.11.2018
 */

#ifndef HSLSOLVER_H
#define HSLSOLVER_H

#include "hsl.h"
#include "../lapack_wrapper/lapack_wrapper_config.hh"
#include <string>
#include <algorithm>

namespace lapack_wrapper {

  template <typename real>
  class HSLsolver {
  protected:
    int numRows;
    int numCols;
    int nnz;

    /// True if the HSL solver is initialized.
    bool isInitialized;
    /// True if the HSL solver is factorized.
    bool isFactorized;

    mutable std::string last_error;

    /**
     * \brief HSLsolver:
     *        Constructor of the class HSLsolver.
     */
    HSLsolver();

  public:

    /**
     * \brief ~HSLsolver:
     *        Virtual destructor of the class HSLsolver.
     */
    virtual ~HSLsolver();

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
      int       nnz,
      int       N_Row,
      int       N_Col,
      int const i_Row[],
      int const j_Col[],
      bool      isFortranIndexing
    ) LAPACK_WRAPPER_PURE_VIRTUAL;

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
    factorize( real const ArrayA[] ) LAPACK_WRAPPER_PURE_VIRTUAL;

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
      int        nrhs,
      real const RHS[],
      int        ldRHS,
      real       X[],
      int        ldX
    ) const LAPACK_WRAPPER_PURE_VIRTUAL;

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
    solve(
      real const RHS[],
      real       X[]
    ) const {
      return solve( 1, RHS, this->numRows, X, this->numRows );
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
      int        nrhs,
      real const RHS[],
      int        ldRHS,
      real       X[],
      int        ldX
    ) const LAPACK_WRAPPER_PURE_VIRTUAL;

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
    solve_transposed(
      real const RHS[],
      real       X[]
    ) const {
      return solve_transposed( 1, RHS, this->numCols, X, this->numCols );
    }

    /**
     * \brief checkInitialized:
     *        Checks the flag isInitialized.
     *
     * \return Returns the flag isInitialized.
     *
     */
    bool
    checkInitialized()
    { return this->isInitialized; }

    /**
     * \brief checkFactorized:
     *        Checks the flag isFactorized.
     *
     * \return Returns the flag isFactorized.
     *
     */
    bool
    checkFactorized()
    { return this->isFactorized; }
  };

  #ifdef LAPACK_WRAPPER_USE_CXX11

  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  #endif

  extern template class HSLsolver<float>;
  extern template class HSLsolver<double>;

  #ifdef __clang__
  #pragma clang diagnostic pop
  #endif

  #endif
} // namespace lapack_wrapper

#endif // HSLSOLVER_H
