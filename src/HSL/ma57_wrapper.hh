/*
 * \file ma57_wrapper.hh
 * Header definitions of the class MA57.
 *
 * \author A. Huber and E.Bertolazzi
 * \since  13.11.2018
 */

#pragma once

#ifndef MA57_WRAPPER_dot_HH
#define MA57_WRAPPER_dot_HH

#include "hsl_solver.hh"
#include <cmath>
#include <vector>

namespace lapack_wrapper {

  template <typename real>
  class MA57 : public HSLsolver<real> {
  private:
    // Memory vectors for int.
    std::vector<integer> m_i_Row;
    std::vector<integer> m_j_Col;
    std::vector<integer> m_iKeep;
    std::vector<integer> m_ifact;
    mutable std::vector<integer> m_iPivotSeq;

    // Memory vectors for real.
    std::vector<real> m_a_stored;
    std::vector<real> m_fact;
    mutable std::vector<real> m_RHSRefinement;
    mutable std::vector<real> m_Work;
    mutable std::vector<real> m_residualVec;

    // MA57 Storage:
    /// Array of length 5 with real control parameters.
    real m_cntl[5];
    /// Array of length 20 with int control parameters.
    integer m_icntl[20];
    /// Integer Info variables.
    mutable integer m_iinfo[40];
    /// Real Info variables.
    mutable real m_rinfo[20];

    // Factors in MA57:
    // Workspace MA57:

    /// Integer specifying task (1 for solve AX = B).
    integer m_job;

    /// Specifies whether to refine iteratively.
    bool m_doRefinement;

    /// Residuum tolerance.
    real m_tolRes;
    /// Maximum number of refinements.
    integer m_MaxRefinements;

    /**
     * \brief getResidual:
     *        Get the maximum from the residual vector.
     *
     * \param[in] _residualVec The residual vector.
     *
     * \return The maximum from the residual vector.
     */
    real
    getResidual( std::vector<real> const & _residualVec ) const {
      typename std::vector<real>::const_iterator it = _residualVec.begin();
      real ret = real(0);
      while ( it != _residualVec.end() ) {
        real av = std::abs(*it++);
        if ( av > ret ) ret = av;
      }
      return ret;
    }

  public:

    using HSLsolver<real>::m_nrows;
    using HSLsolver<real>::m_ncols;
    using HSLsolver<real>::m_nnz;
    using HSLsolver<real>::m_isInitialized;
    using HSLsolver<real>::m_isFactorized;
    using HSLsolver<real>::m_last_error;

    /**
     * \brief MA57:
     *        Constructor of the class MA57.
     */
    MA57() : HSLsolver<real>() {
      // Default values:
      m_doRefinement   = true;
      m_tolRes         = real(1e-11);
      m_MaxRefinements = 100;
    }

    /**
     * \brief ~MA57:
     *        Virtual destructor of the class MA57.
     */
    virtual ~MA57() override { }

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
     * \param[in] Nnz               Number of non zeros.
     * \param[in] N_Row             Row dimension of the matrix.
     * \param[in] N_Col             Col dimension of the matrix.
     * \param[in] i_Row             Rows of the matrix.
     * \param[in] j_Col             Cols of the matrix.
     * \param[in] isFortranIndexing True for Fortran indexing and false for C indexing.
     *
     * \return True if HSL could be initialized.
     *
     */
    bool
    init(
      integer       Nnz,
      integer       N_Row,
      integer       N_Col,
      integer const i_Row[],
      integer const j_Col[],
      bool          isFortranIndexing,
      bool          isStoredSymmetric
    ) override;

    /**
     * \brief factorize:
     *        Factorizes the linear system \f$ Ax = b\f$ for \f$ A\f$ for a sparse matrix.
     *
     * \param[in] ArrayA Values of the sparse matrix.
     *
     * \return True if the linear system \f$ Ax = b\f$ for \f$ A\f$ could be solved.
     *
     */
    bool
    factorize( real const ArrayA[] ) override;

    bool
    solve(
      integer    nrhs,
      real const RHS[],
      integer    ldRHS,
      real       X[],
      integer    ldX
    ) const override;

    bool
    solve_transposed(
      integer    nrhs,
      real const RHS[],
      integer    ldRHS,
      real       X[],
      integer    ldX
    ) const override {
      // Symmetric matrix:
      return this->solve(nrhs, RHS, ldRHS, X, ldX);
    }

    /**
     * \brief setSolverMode:
     *        Switches to a reinfined solver mode.
     *
     * \param[in] DoRefinement    Tells wether so activate reinfined solver mode.
     * \param[in] TolRes          Residuum tolerance.
     * \param[in] maxRefinements  Maximum number of refinements.
     *
     */
    void
    setSolverMode(
      bool    DoRefinement,
      real    TolRes         = real(1e-11),
      integer maxRefinements = 100
    ) {
      m_doRefinement   = DoRefinement;
      m_tolRes         = TolRes;
      m_MaxRefinements = maxRefinements;
    }
  };

  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  #pragma clang diagnostic ignored "-Wweak-template-vtables"
  #endif

  extern template class MA57<float>;
  extern template class MA57<double>;

  #ifdef __clang__
  #pragma clang diagnostic pop
  #endif

} // namespace lapack_wrapper

#endif // MA57_WRAPPER_HH
