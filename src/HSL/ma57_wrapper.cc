/*
 * \file ma57_wrapper.cc
 * Implementation of the class CPPMA57.
 *
 * \author A. Huber and E.Bertolazzi
 * \since 13.11.2018
 */

#include "ma57_wrapper.hh"
#include <iostream>

namespace lapack_wrapper {

  using std::copy_n;

  /*\
  :|:   ___      _    _ _      __  __           _
  :|:  | _ \_  _| |__| (_)__  |  \/  |___ _ __ | |__  ___ _ _ ___
  :|:  |  _/ || | '_ \ | / _| | |\/| / -_) '  \| '_ \/ -_) '_(_-<
  :|:  |_|  \_,_|_.__/_|_\__| |_|  |_\___|_|_|_|_.__/\___|_| /__/
  \*/
  template <typename real>
  bool
  MA57<real>::init(
    int  const Nnz,
    int  const N_Row,
    int  const N_Col,
    int  const i_Row[],
    int  const j_Col[],
    bool const isFortranIndexing,
    bool const isStoredSymmetric
  ) {
    // Set dimension:
    m_nnz = Nnz;
    if (N_Row != N_Col) {
      m_isInitialized = false;
      m_isFactorized  = false;
      m_last_error    = "CPPMA57<real>::init N_Row != N_Col!";
      return false;
    }
    m_nrows = m_ncols = N_Row;
    // allocate memory:
    m_i_Row.resize(static_cast<size_t>(m_nnz));
    m_j_Col.resize(static_cast<size_t>(m_nnz));
    m_iKeep.resize(static_cast<size_t>(5 * m_nrows + 2 * m_nnz + 42));
    m_iPivotSeq.resize(static_cast<size_t>(5 * m_nrows));

    // Copy row and column arrays:
    std::copy_n(i_Row, m_nnz, m_i_Row.begin());
    std::copy_n(j_Col, m_nnz, m_j_Col.begin());

    if (!isFortranIndexing) {
      // Correct fortran indexing:
      for ( integer i{0}; i < m_nnz; ++i ) {
        ++m_i_Row[i];
        ++m_j_Col[i];
      }
    }

    m_job = 1;
    // Initialize MA57:
    lapack_wrapper::ma57i<real>( m_cntl, m_icntl );
    // Analyse the structure of the linear system:
    int const NiKeep{static_cast<int>(m_iKeep.size())};
    lapack_wrapper::ma57a<real>(
      m_nrows,
      m_nnz,
      &m_i_Row.front(),
      &m_j_Col.front(),
      NiKeep,
      &m_iKeep.front(),
      &m_iPivotSeq.front(),
      m_icntl,
      m_iinfo,
      m_rinfo
    );
    // Restize memory for factorization:
    m_fact.resize(static_cast<size_t>(2 * m_iinfo[8]));
    m_ifact.resize(static_cast<size_t>(2 * m_iinfo[9]));
    // Restize memory for workspace:
    m_Work.resize(static_cast<size_t>(3 * m_nnz));
    // allocate reintinmet memory:
    m_RHSRefinement.resize(static_cast<size_t>(m_nrows));
    m_residualVec.resize(static_cast<size_t>(m_nrows));
    m_isInitialized = true;
    m_isFactorized  = false;
    return true;
  }

  template <typename real>
  bool
  MA57<real>::factorize( real const ArrayA[] ) {
    if (!m_isInitialized) {
      m_last_error   = "The function factorize can only be called after calling CPPMA57<real>::init!";
      m_isFactorized = false;
      return false;
    }

    m_a_stored.resize(static_cast<size_t>(m_nnz));
    std::copy_n(ArrayA, m_nnz, m_a_stored.begin());

    // Now factorize the system:
    int const NiKeep = static_cast<int>(m_iKeep.size());
    int const lifact = static_cast<int>(m_ifact.size());
    int const lfact  = static_cast<int>(m_fact.size());
    lapack_wrapper::ma57b<real>(
      m_nrows,
      m_nnz,
      &m_a_stored.front(),
      &m_fact.front(),
      lfact,
      &m_ifact.front(),
      lifact,
      NiKeep,
      &m_iKeep.front(),
      &m_iPivotSeq.front(),
      m_icntl,
      m_cntl,
      m_iinfo,
      m_rinfo
    );

    m_isFactorized = m_iinfo[0] == 0;
    return m_isFactorized;
  }

  template <typename real>
  bool
  MA57<real>::solve(
    int  const nrhs,
    real const RHS[],
    int  const ldRHS,
    real       X[],
    int  const ldX
  ) const {
    // Check valid call:
    if (!m_isInitialized) {
      m_last_error = "Can not solve uninitialized system with CPPMA57::solve!";
      return false;
    }
    if (!m_isFactorized) {
      m_last_error = "Can not solve unfactorized system with CPPMA57::solve!";
      return false;
    }

    // Solve system:
    if ( m_doRefinement ) {
      std::copy_n( RHS, m_nrows, m_RHSRefinement.begin() );
      real residual = 10.0;
      int  localjob = 0;
      for ( int counter = 0;
            counter < m_MaxRefinements && std::abs(residual) > m_tolRes;
            ++counter ) {
        // Löse mit rechter Seite: (iterative reinfinement)
        int lifact = static_cast<int>(m_ifact.size());
        int lfact  = static_cast<int>(m_fact.size());
        lapack_wrapper::ma57d<real>(
          localjob,
          m_nrows,
          static_cast<int>(m_a_stored.size()),
          &m_a_stored.front(),
          &m_i_Row.front(),
          &m_j_Col.front(),
          &m_fact.front(),
          lfact,
          &m_ifact.front(),
          lifact,
          &m_RHSRefinement.front(),
          X,
          &m_residualVec.front(),
          &m_Work.front(),
          &m_iPivotSeq.front(),
          m_icntl,
          m_cntl,
          m_iinfo,
          m_rinfo
        );
        residual = this->getResidual(m_residualVec);
        localjob = 2;
        if(m_iinfo[0] < 0) return false;
        if (counter >= m_MaxRefinements - 1) {
          m_last_error = "Refinement of HSL MA57 failed (Maxiter)";
          if ( residual > real(1e-5) ) return false;
          fmt::print("Residual: {:g}\n", residual);
        }
      }
    } else {
      // Solve with right side:
      if (RHS != X)
        for (int j{0}; j < nrhs; ++j)
           copy_n(RHS + j * ldRHS, m_nrows, X + j * ldX);
      int const lfact  { static_cast<int>(m_fact.size())  };
      int const lifact { static_cast<int>(m_ifact.size()) };
      int const NWork  { static_cast<int>(m_Work.size())  };
      lapack_wrapper::ma57c<real>(
        m_job,
        m_nrows,
        &m_fact.front(),
        lfact,
        &m_ifact.front(),
        lifact,
        nrhs,
        X,
        ldX,
        &m_Work.front(),
        NWork,
        &m_iPivotSeq.front(),
        m_icntl,
        m_iinfo
      );
    }
    return m_iinfo[0] == 0;
  }

  template class MA57<double>;
  template class MA57<float>;

} // namespace lapack_wrapper
