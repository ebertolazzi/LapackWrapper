/**
 * \file ma57_wrapper.cc
 * Implementation of the class CPPMA57.
 *
 * \author A. Huber and E.Bertolazzi
 * \since 13.11.2018
 */

#include "ma57_wrapper.hh"
#include <iostream>

namespace lapack_wrapper {

  /*\
  :|:   ___      _    _ _      __  __           _
  :|:  | _ \_  _| |__| (_)__  |  \/  |___ _ __ | |__  ___ _ _ ___
  :|:  |  _/ || | '_ \ | / _| | |\/| / -_) '  \| '_ \/ -_) '_(_-<
  :|:  |_|  \_,_|_.__/_|_\__| |_|  |_\___|_|_|_|_.__/\___|_| /__/
  \*/  
  template <typename real>
  bool
  MA57<real>::init(
    int       Nnz,
    int       N_Row,
    int       N_Col,
    int const iRow[],
    int const jCol[],
    bool      isFortranIndexing
  ) {
    // Set dimension:
    this->nnz = Nnz;
    if (N_Row != N_Col) {
      this->isInitialized = false;
      this->isFactorized  = false;
      this->last_error    = "CPPMA57<real>::init N_Row != N_Col!";
      return false;
    }
    this->numRows = this->numCols = N_Row;
    // allocate memory:
    this->i_Row.resize(size_t(this->nnz));
    this->j_Col.resize(size_t(this->nnz));
    this->iKeep.resize(size_t(5 * this->numRows + 2 * this->nnz + 42));
    this->iPivotSeq.resize(size_t(5 * this->numRows));

    // Copy row and column arrays:
    std::copy_n(iRow, this->nnz, this->i_Row.begin());
    std::copy_n(jCol, this->nnz, this->j_Col.begin());

    if (!isFortranIndexing) {
      // Correct fortran indexing:
      for ( size_t i = 0; i < size_t(this->nnz); ++i ) {
        ++this->i_Row[i];
        ++this->j_Col[i];
      }
    }

    this->job = 1;
    // Initialize MA57:
    HSL::ma57i<real>(this->cntl, this->icntl);
    // Analyse the structure of the linear system:
    int NiKeep = int(this->iKeep.size());
    HSL::ma57a<real>(this->numRows, this->nnz, &this->i_Row.front(), &this->j_Col.front(), NiKeep, &this->iKeep.front(), &this->iPivotSeq.front(), this->icntl, this->iinfo, this->rinfo);
    // Restize memory for factorization:
    this->fact.resize(size_t(2 * this->iinfo[8]));
    this->ifact.resize(size_t(2 * this->iinfo[9]));
    // Restize memory for workspace:
    this->Work.resize(size_t(3 * this->nnz));
    // allocate reintinmet memory:
    this->RHSReinfinement.resize(size_t(this->numRows));
    this->residualVec.resize(size_t(this->numRows));
    this->isInitialized = true;
    this->isFactorized  = false;
    return true;
  }

  template <typename real>
  bool
  MA57<real>::factorize( real const ArrayA[] ) {
    if (!this->isInitialized) {
      this->last_error   = "The function factorize can only be called after calling CPPMA57<real>::init!";
      this->isFactorized = false;
      return false;
    }

    a_stored.resize(size_t(this->nnz));
    std::copy_n(ArrayA, this->nnz, a_stored.begin());

    // Now factorize the system:
    int NiKeep = int(this->iKeep.size());
    int lifact = int(this->ifact.size());
    int lfact  = int(this->fact.size());
    HSL::ma57b<real>(
      this->numRows,
      this->nnz,
      &a_stored.front(),
      &this->fact.front(),
      lfact,
      &this->ifact.front(),
      lifact,
      NiKeep,
      &this->iKeep.front(),
      &this->iPivotSeq.front(),
      this->icntl,
      this->cntl,
      this->iinfo,
      this->rinfo
    );

    this->isFactorized = this->iinfo[0] == 0;
    return this->isFactorized;
  }

  template <typename real>
  bool
  MA57<real>::solve(
    int        nrhs,
    real const RHS[],
    int        ldRHS,
    real       X[],
    int        ldX
  ) const {
    // Check valid call:
    if (!this->isInitialized) {
      this->last_error = "Can not solve uninitialized system with CPPMA57::solve!";
      return false;
    }
    if (!this->isFactorized) {
      this->last_error = "Can not solve unfactorized system with CPPMA57::solve!";
      return false;
    }

    // Solve system:
    if (this->doRefinement) {
      std::copy( RHS, RHS+this->numRows, this->RHSReinfinement.begin() );
      real residual = 10.0;
      int  localjob = 0;
      for ( int counter = 0;
            counter < this->MaxRefinements && std::abs(residual) > this->tolRes;
            ++counter ) {
        // Löse mit rechter Seite: (iterative reinfinement)
        int lifact = int(this->ifact.size());
        int lfact  = int(this->fact.size());
        HSL::ma57d<real>(
          localjob,
          this->numRows,
          int(a_stored.size()),
          &a_stored.front(),
          &this->i_Row.front(),
          &this->j_Col.front(),
          &this->fact.front(),
          lfact,
          &this->ifact.front(),
          lifact,
          &this->RHSReinfinement.front(),
          X,
          &this->residualVec.front(),
          &this->Work.front(),
          &this->iPivotSeq.front(),
          this->icntl,
          this->cntl,
          this->iinfo,
          this->rinfo
        );
        residual = this->getResidual(this->residualVec);
        localjob = 2;
        if(this->iinfo[0] < 0) return false;
        if (counter >= this->MaxRefinements - 1) {
          this->last_error = "Refinement of HSL MA57 failed (Maxiter)";
          if ( residual > real(1e-5) ) return false;
          std::cout
            << std::scientific << "Residual: " << residual
            << std::fixed << std::endl;
        }
      }
    } else {
      // Solve with right side:
      if (RHS != X)
        for (int j = 0; j < nrhs; ++j)
           std::copy_n(RHS + j * ldRHS, this->numRows, X + j * ldX);
      int lfact  = int(this->fact.size());
      int lifact = int(this->ifact.size());
      int NWork  = int(this->Work.size());
      HSL::ma57c<real>(
        this->job,
        this->numRows,
        &this->fact.front(),
        lfact,
        &this->ifact.front(),
        lifact,
        nrhs,
        X,
        ldX,
        &this->Work.front(),
        NWork,
        &this->iPivotSeq.front(),
        this->icntl,
        this->iinfo
      );
    }
    return this->iinfo[0] == 0;
  }

  template class MA57<double>;
  template class MA57<float>;

} // namespace lapack_wrapper
