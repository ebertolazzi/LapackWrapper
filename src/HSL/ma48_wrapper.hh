/**
 * \file ma48_wrapper.hh
 * Header definitions of the class CPPMA48.
 *
 * \author A. Huber and E.Bertolazzi
 * \since 15.11.2018
 */

#ifndef MA48_WRAPPER_HH
#define MA48_WRAPPER_HH

#include "hsl_solver.hh"
#include <vector>

namespace lapack_wrapper {

  template <typename real>
  class MA48 : public HSLsolver<real> {
  private:
    // Memory vectors for int.
    std::vector<int> i_Row_stored;
    std::vector<int> j_Col_stored;
    std::vector<int> irn;
    std::vector<int> jcn;
    std::vector<int> iw;
    std::vector<int> keep;

    // Memory vectors for real.
    std::vector<real> a;
    std::vector<real> w;

    /// Integer specifying task (1 for solve AX = B).
    int job;

    /// Length of the arrays a, irn and jcn (mostly 3 * DimSparse).
    int la;

    /// Array of length 5 with real control parameters.
    real cntl[10];

    /// Array of length 20 with int control parameters.
    int icntl[20];

    /// Error Array von MA48.
    mutable real ma48errors[3];

    /// Integer Info Variables.
    mutable int info[20];
    /// Real Info Variables.
    mutable real rinfo[10];

    void load_error_string( int info ) const;

  public:

    /**
     * \brief MA48:
     *        Constructor of the class MA48.
     */
    MA48() {}

    /**
      * \brief ~MA48:
      *        Virtual destructor of the class MA48.
      */
    virtual ~MA48() LAPACK_WRAPPER_OVERRIDE { }

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

    bool
    solve(
      real const RHS[],
      real       X[],
      bool       itrans
    ) const;

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
    ) const LAPACK_WRAPPER_OVERRIDE;
  };

  #ifdef LAPACK_WRAPPER_USE_CXX11

  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  #endif

  extern template class MA48<float>;
  extern template class MA48<double>;

  #ifdef __clang__
  #pragma clang diagnostic pop
  #endif

  #endif

} // namespace lapack_wrapper

#endif // MA48_WRAPPER_HH
