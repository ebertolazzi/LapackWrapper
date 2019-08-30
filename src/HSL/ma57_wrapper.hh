/**
 * \file ma57_wrapper.hh
 * Header definitions of the class MA57.
 *
 * \author A. Huber and E.Bertolazzi
 * \since  13.11.2018
 */

#ifndef MA57_WRAPPER_HH
#define MA57_WRAPPER_HH

#include "hsl_solver.hh"
#include <cmath>
#include <vector>

namespace lapack_wrapper {

  template <typename real>
  class MA57 : public HSLsolver<real> {
  private:
    // Memory vectors for int.
    std::vector<int> i_Row;
    std::vector<int> j_Col;
    std::vector<int> iKeep;
    std::vector<int> ifact;
    mutable std::vector<int> iPivotSeq;

    // Memory vectors for real.
    std::vector<real> a_stored;
    std::vector<real> fact;
    mutable std::vector<real> RHSReinfinement;
    mutable std::vector<real> Work;
    mutable std::vector<real> residualVec;

    // MA57 Storage:
    /// Array of length 5 with real control parameters.
    real cntl[5];
    /// Array of length 20 with int control parameters.
    int icntl[20];
    /// Integer Info variables.
    mutable int iinfo[40];
    /// Real Info variables.
    mutable real rinfo[20];

    // Factors in MA57:
    // Workspace MA57:

    /// Integer specifying task (1 for solve AX = B).
    int job;

    /// Specifies whether to refine iteratively.
    bool doRefinement;

    /// Residuum tolerance.
    real tolRes;
    /// Maximum number of refinements.
    int MaxRefinements;

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

    /**
     * \brief MA57:
     *        Constructor of the class MA57.
     */
    MA57() : HSLsolver<real>() {
      // Default values:
      this->doRefinement   = true;
      this->tolRes         = real(1e-11);
      this->MaxRefinements = 100;
    }

    /**
     * \brief ~MA57:
     *        Virtual destructor of the class MA57.
     */
    virtual ~MA57() LAPACK_WRAPPER_OVERRIDE { }

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
    ) LAPACK_WRAPPER_OVERRIDE;

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
    factorize( real const ArrayA[] ) LAPACK_WRAPPER_OVERRIDE;

    virtual
    bool
    solve(
      int        nrhs,
      real const RHS[],
      int        ldRHS,
      real       X[],
      int        ldX
    ) const LAPACK_WRAPPER_OVERRIDE;

    virtual
    bool
    solve_transposed(
      int        nrhs,
      real const RHS[],
      int        ldRHS,
      real       X[],
      int        ldX
    ) const LAPACK_WRAPPER_OVERRIDE {
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
      bool DoRefinement,
      real TolRes         = real(1e-11),
      int  maxRefinements = 100
    ) {
      this->doRefinement   = DoRefinement;
      this->tolRes         = TolRes;
      this->MaxRefinements = maxRefinements;
    }
  };

  #ifdef LAPACK_WRAPPER_USE_CXX11

  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  #endif

  extern template class MA57<float>;
  extern template class MA57<double>;

  #ifdef __clang__
  #pragma clang diagnostic pop
  #endif

  #endif

} // namespace lapack_wrapper

#endif // MA57_WRAPPER_HH
