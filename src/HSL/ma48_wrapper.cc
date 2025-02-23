/*
 * \file ma48_wrapper.cc
 * Implementation of the class MA48.
 *
 * \author A. Huber and E.Bertolazzi
 * \since 15.11.2018
 */

#include "ma48_wrapper.hh"
#include "hsl.hh"
#include <iostream>

// windows workaround
#ifdef max
  #undef max
#endif

#ifdef min
  #undef min
#endif

namespace lapack_wrapper {

  /*\
  :|:   ___      _    _ _      __  __           _
  :|:  | _ \_  _| |__| (_)__  |  \/  |___ _ __ | |__  ___ _ _ ___
  :|:  |  _/ || | '_ \ | / _| | |\/| / -_) '  \| '_ \/ -_) '_(_-<
  :|:  |_|  \_,_|_.__/_|_\__| |_|  |_\___|_|_|_|_.__/\___|_| /__/
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename real>
  void
  MA48<real>::load_error_string( int const info ) const {
    switch (info) {
    case -1:
      m_last_error = "MA48::factorize, Number of row or columns < 1";
      break;
    case -2:
      m_last_error = "MA48::factorize, Number of nonzeros < 1";
      break;
    case -3:
      m_last_error = "MA48::factorize, Correction of LA failed!";
      break;
    case -4:
      m_last_error = "MA48::factorize, matrix is structurally rank deficient";
      break;
    case -5:
      m_last_error = "MA48::factorize, faulty permutation inputin KEEP when JOB=2.";
      break;
    case -6:
      m_last_error = "MA48::factorize, JOB has a value less than 1 or greater than 3";
      break;
    case -7:
      m_last_error = "MA48::factorize, On a call with JOB = 2, the matrix entries are unsuitable for the pivot sequence chosen on the previous call";
      break;
    case -8:
      m_last_error = "MA48::factorize, Iterative refinement has not converged";
      break;
    case -9:
      m_last_error = "MA48::factorize, A problem has occurred in the calculation of matrix norms using MC71A/AD";
      break;
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
      m_last_error = "CPPMA48::factorize";
      if ((m_info[0] & 0x01) != 0)
        m_last_error += "\nOne or more row or columns indices are out of range or one or more entries are for the same position in the matrix, or both are true";
      if ((m_info[0] & 0x02) != 0)
        m_last_error += "\nThe matrix is rank deficient";
      if ((m_info[0] & 0x04) != 0)
        m_last_error += "\nNot possible to choose all pivots from diagonal";
      break;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename real>
  bool
  MA48<real>::init(
    int  const Nnz,
    int  const N_Row,
    int  const N_Col,
    int  const i_Row[],
    int  const j_Col[],
    bool const isFortranIndexing,
    bool const isStoredSymmetric
  ) {
    m_isFactorized = false;
    m_nnz          = Nnz;
    m_nrows        = N_Row;
    m_ncols        = N_Col;
    m_la           = 3 * m_nnz;
    m_job          = 1;
    m_a.resize(static_cast<size_t>(m_la));
    m_irn.resize(static_cast<size_t>(m_la));
    m_jcn.resize(static_cast<size_t>(m_la));
    m_i_Row_stored.resize(static_cast<size_t>(m_nnz));
    m_j_Col_stored.resize(static_cast<size_t>(m_nnz));

    // Copy data:
    copy_n( i_Row, m_nnz, m_i_Row_stored.begin() );
    copy_n( j_Col, m_nnz, m_j_Col_stored.begin() );
    if ( !isFortranIndexing) {
      // Correct fortran indexing:
      for ( integer i{0}; i < m_nnz; ++i ) {
        ++m_i_Row_stored[i];
        ++m_j_Col_stored[i];
      }
    }

    std::copy( m_i_Row_stored.begin(), m_i_Row_stored.end(), m_irn.begin() );
    std::copy( m_j_Col_stored.begin(), m_j_Col_stored.end(), m_jcn.begin() );

    // Initialize MA48:
    lapack_wrapper::ma48i<real>( m_cntl, m_icntl );

    int ihlp = m_icntl[5];
    if ( ihlp <= 0 ) {
      //ihlp            = 1;
      m_last_error    = "CPPMA48::init, wrong value of icntl[5] in MA48. Set to 1!";
      m_isInitialized = false;
      return false;
    }
    int const N_RC_Max = std::max( m_nrows, m_ncols );
    int const N_keep   = m_nrows + 5 * m_ncols + 4 * (m_ncols / ihlp) + 7;
    int const N_iw     = m_nrows * 6 + m_ncols * 3;
    int const N_w      = 4 * N_RC_Max;
    m_keep.resize(static_cast<size_t>(N_keep));
    m_iw.resize(static_cast<size_t>(N_iw));
    m_w.resize(static_cast<size_t>(N_w));
    m_isInitialized = true;
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename real>
  bool
  MA48<real>::factorize( real const ArrayA[] ) {
    // Check valid call:
    if (!m_isInitialized) {
      m_last_error   = "MA48::factorize, The function factorize can only be called after calling init!";
      m_isFactorized = false;
      return false;
    }

    // Copy Memory:
    std::copy_n( ArrayA, m_nnz, m_a.begin() );
    // Pivoting:
    lapack_wrapper::ma48a<real>(
      m_nrows,
      m_ncols,
      m_nnz,
      m_job,
      m_la,
      &m_a.front(),
      &m_irn.front(),
      &m_jcn.front(),
      &m_keep.front(),
      m_cntl,
      m_icntl,
      &m_iw.front(),
      m_info,
      m_rinfo
    );

    if ( m_info[3] > m_la ) {
      // Increase internal workspace:
      m_la = m_info[3];
      m_a.resize(static_cast<size_t>(m_la));
      m_irn.resize(static_cast<size_t>(m_la));
      m_jcn.resize(static_cast<size_t>(m_la));
      std::copy_n( ArrayA, m_nnz, m_a.begin());
      std::copy(
        m_i_Row_stored.begin(),
        m_i_Row_stored.end(),
        m_irn.begin()
      );
      std::copy(
        m_j_Col_stored.begin(),
        m_j_Col_stored.end(),
        m_jcn.begin()
      );

      lapack_wrapper::ma48a<real>(
        m_nrows,
        m_ncols,
        m_nnz,
        m_job,
        m_la,
        &m_a.front(),
        &m_irn.front(),
        &m_jcn.front(),
        &m_keep.front(),
        m_cntl,
        m_icntl,
        &m_iw.front(),
        m_info,
        m_rinfo
      );
      // Errors::set_Warning("Workspace of MA48 to small! Correcting!");
    }

    load_error_string( m_info[0] );

    m_isFactorized = false;
    if ( m_info[0] != 0 ) return false;

    // Factorize:
    lapack_wrapper::ma48b<real>(
      m_nrows,
      m_ncols,
      m_nnz,
      m_job,
      m_la,
      &m_a.front(),
      &m_irn.front(),
      &m_jcn.front(),
      &m_keep.front(),
      m_cntl,
      m_icntl,
      &m_w.front(),
      &m_iw.front(),
      m_info,
      m_rinfo
    );

    load_error_string( m_info[0] );

    m_isFactorized = m_info[0] == 0;
    return m_isFactorized;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename real>
  bool
  MA48<real>::solve( real const RHS[], real X[], bool const itrans ) const {
    // Check valid call:
    if (!m_isInitialized) {
      m_last_error = "MA48::solve, Can not solve uninitialized system!";
      return false;
    }
    if (!m_isFactorized) {
      m_last_error = "MA48::solve, Can not solve unfactorized system!";
      return false;
    }

    // Solve with right hand side:
    lapack_wrapper::ma48c<real>(
      m_nrows,
      m_ncols,
      itrans ? 1 : 0,
      m_job,
      m_la,
      &m_a.front(),
      &m_irn.front(),
      &m_keep.front(),
      m_cntl,
      m_icntl,
      RHS,
      X,
      ma48errors,
      &m_w.front(),
      &m_iw.front(),
      m_info
    );

    load_error_string( m_info[0] );
    return m_info[0] == 0;
  }

  template <typename real>
  bool
  MA48<real>::solve(
    int  const nrhs,
    real const RHS[],
    int  const ldRHS,
    real       X[],
    int  const ldX
  ) const {
    bool ok{ true };
    if ( RHS == X ) {
      Malloc<real> mem("MA48<real>::solve");
      real * tmpRHS = mem.realloc( m_nrows );
      for ( int j{0}; ok && j < nrhs; ++j ) {
        real const * RHSj{ RHS + j * ldRHS };
        std::copy_n( RHSj, m_nrows, tmpRHS );
        ok = this->solve( tmpRHS, X + j * ldX, false );
      }
    } else {
      for ( int j{0}; ok && j < nrhs; ++j )
        ok = this->solve(RHS + j * ldRHS, X + j * ldX, false);
    }
    return ok;
  }

  template <typename real>
  bool
  MA48<real>::solve_transposed(
    int  const nrhs,
    real const RHS[],
    int  const ldRHS,
    real       X[],
    int  const ldX
  ) const {
    bool ok{true};
    if ( RHS == X ) {
      vector<real> tmpRHS;
      tmpRHS.resize(static_cast<size_t>(m_ncols));
      for ( int j{0}; ok && j < nrhs; ++j ) {
        real const * RHSj{ RHS + j * ldRHS };
        copy_n( RHSj, m_ncols, tmpRHS.begin() );
        ok = this->solve(&tmpRHS.front(), X + j * ldX, true);
      }
    } else {
      for ( int j{0}; ok && j < nrhs; ++j )
        ok = this->solve(RHS + j * ldRHS, X + j * ldX, true);
    }
    return ok;
  }

  template class MA48<double>;
  template class MA48<float>;

} // namespace lapack_wrapper
