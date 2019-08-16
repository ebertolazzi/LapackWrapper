/**
 * \file ma48_wrapper.cc
 * Implementation of the class MA48.
 *
 * \author A. Huber and E.Bertolazzi
 * \since 15.11.2018
 */

#include "ma48_wrapper.hh"
#include "hsl.h"
#include <iostream>

namespace lapack_wrapper {

  /*\
  :|:   ___      _    _ _      __  __           _
  :|:  | _ \_  _| |__| (_)__  |  \/  |___ _ __ | |__  ___ _ _ ___
  :|:  |  _/ || | '_ \ | / _| | |\/| / -_) '  \| '_ \/ -_) '_(_-<
  :|:  |_|  \_,_|_.__/_|_\__| |_|  |_\___|_|_|_|_.__/\___|_| /__/
  \*/

  template <>
  void
  MA48<float>::ma48i(float _cntl[], int _icntl[]) const {
    HSL_F77NAME(ma48i)(_cntl, _icntl);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <>
  void
  MA48<double>::ma48i(double _cntl[], int _icntl[]) const {
    HSL_F77NAME(ma48id)(_cntl, _icntl);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <>
  void
  MA48<float>::ma48a(
    int         _nrow,
    int         _ncol,
    int         _nz,
    int         _job,
    int         _la,
    float       _a[],
    int         _irn[],
    int         _jcn[],
    int         _keep[],
    float const _cntl[10],
    int const   _icntl[20],
    int         _iw[],
    int         _info[20],
    float       _rinfo[10]
  ) const {
    HSL_F77NAME(ma48a)(
      &_nrow, &_ncol, &_nz, &_job, &_la, _a, _irn, _jcn, _keep,
      _cntl, _icntl, _iw, _info, _rinfo
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <>
  void
  MA48<double>::ma48a(
    int          _nrow,
    int          _ncol,
    int          _nz,
    int          _job,
    int          _la,
    double       _a[],
    int          _irn[],
    int          _jcn[],
    int          _keep[],
    double const _cntl[10],
    int const    _icntl[20],
    int          _iw[],
    int          _info[20],
    double       _rinfo[10]
  ) const {
    HSL_F77NAME(ma48ad)(
      &_nrow, &_ncol, &_nz, &_job, &_la, _a, _irn, _jcn, _keep,
      _cntl, _icntl, _iw, _info, _rinfo
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <>
  void
  MA48<float>::ma48b(
    int const   _nrow,
    int const   _ncol,
    int const   _nz,
    int const   _job,
    int const   _la,
    float       _a[],
    int         _irn[],
    int         _jcn[],
    int const   _keep[],
    float const _cntl[10],
    int const   _icntl[20],
    float       _w[],
    int         _iw[],
    int         _info[20],
    float       _rinfo[10]
  ) const {
    HSL_F77NAME(ma48b)(
      &_nrow, &_ncol, &_nz, &_job, &_la, _a, _irn, _jcn, _keep,
      _cntl, _icntl, _w, _iw, _info, _rinfo
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <>
  void
  MA48<double>::ma48b(
    int const    _nrow,
    int const    _ncol,
    int const    _nz,
    int const    _job,
    int const    _la,
    double       _a[],
    int          _irn[],
    int          _jcn[],
    int const    _keep[],
    double const _cntl[10],
    int const    _icntl[20],
    double       _w[],
    int          _iw[],
    int          _info[20],
    double       _rinfo[10]
  ) const {
    HSL_F77NAME(ma48bd)(
      &_nrow, &_ncol, &_nz, &_job, &_la, _a, _irn, _jcn, _keep,
      _cntl, _icntl, _w, _iw, _info, _rinfo
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <>
  void
  MA48<float>::ma48c(
    int         _nrow,
    int         _ncol,
    int         _itrans,
    int         _job,
    int         _la,
    float const _a[],
    int const   _irn[],
    int const   _keep[],
    float const _cntl[10],
    int const   _icntl[20],
    float const _rhs[],
    float       _x[],
    float       _errors[3],
    float const _w[],
    int const   _iw[],
    int         _info[20]
  ) const {
    HSL_F77NAME(ma48c)(
      &_nrow, &_ncol, &_itrans, &_job, &_la, _a, _irn, _keep,
      _cntl, _icntl, _rhs, _x, _errors, _w, _iw, _info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <>
  void
  MA48<double>::ma48c(
    int          _nrow,
    int          _ncol,
    int          _itrans,
    int          _job,
    int          _la,
    double const _a[],
    int const    _irn[],
    int const    _keep[],
    double const _cntl[10],
    int const    _icntl[20],
    double const _rhs[],
    double       _x[],
    double       _errors[3],
    double const _w[],
    int const    _iw[],
    int          _info[20]
  ) const {
    HSL_F77NAME(ma48cd)(
      &_nrow, &_ncol, &_itrans, &_job, &_la, _a, _irn, _keep,
      _cntl, _icntl, _rhs, _x, _errors, _w, _iw, _info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename real>
  void
  MA48<real>::load_error_string( int _info ) const {
    switch (_info) {
    case -1:
      this->last_error = "MA48::factorize, Number of row or columns < 1";
      break;
    case -2:
      this->last_error = "MA48::factorize, Number of nonzeros < 1";
      break;
    case -3:
      this->last_error = "MA48::factorize, Correction of LA failed!";
      break;
    case -4:
      this->last_error = "MA48::factorize, matrix is structurally rank deficient";
      break;
    case -5:
      this->last_error = "MA48::factorize, faulty permutation inputin KEEP when JOB=2.";
      break;
    case -6:
      this->last_error = "MA48::factorize, JOB has a value less than 1 or greater than 3";
      break;
    case -7:
      this->last_error = "MA48::factorize, On a call with JOB = 2, the matrix entries are unsuitable for the pivot sequence chosen on the previous call";
      break;
    case -8:
      this->last_error = "MA48::factorize, Iterative refinement has not converged";
      break;
    case -9:
      this->last_error = "MA48::factorize, A problem has occurred in the calculation of matrix norms using MC71A/AD";
      break;
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
      this->last_error = "CPPMA48::factorize";
      if ((this->info[0] & 0x01) != 0)
        this->last_error += "\nOne or more row or columns indices are out of range or one or more entries are for the same position in the matrix, or both are true";
      if ((this->info[0] & 0x02) != 0)
        this->last_error += "\nThe matrix is rank deficient";
      if ((this->info[0] & 0x04) != 0)
        this->last_error += "\nNot possible to choose all pivots from diagonal";
      break;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename real>
  bool
  MA48<real>::init(
    int       Nnz,
    int       N_Row,
    int       N_Col,
    int const i_Row[],
    int const j_Col[],
    bool      isFortranIndexing
  ) {
    this->isFactorized = false;
    this->nnz          = Nnz;
    this->numRows      = N_Row;
    this->numCols      = N_Col;
    this->la           = 3 * this->nnz;
    this->job          = 1;
    this->a.resize(size_t(this->la));
    this->irn.resize(size_t(this->la));
    this->jcn.resize(size_t(this->la));
    this->i_Row_stored.resize(size_t(this->nnz));
    this->j_Col_stored.resize(size_t(this->nnz));

    // Copy data:
    std::copy(i_Row, i_Row+this->nnz, this->i_Row_stored.begin());
    std::copy(j_Col, j_Col+this->nnz, this->j_Col_stored.begin());
    if ( !isFortranIndexing) {
      // Correct fortran indexing:
      for ( size_t i = 0; i < size_t(this->nnz); ++i ) {
        ++this->i_Row_stored[i];
        ++this->j_Col_stored[i];
      }
    }

    std::copy(i_Row_stored.begin(), i_Row_stored.end(), this->irn.begin());
    std::copy(j_Col_stored.begin(), j_Col_stored.end(), this->jcn.begin());

    // Initialize MA48:
    this->ma48i(this->cntl, this->icntl);

    int ihlp = this->icntl[5];
    if ( ihlp <= 0 ) {
      ihlp                = 1;
      this->last_error    = "CPPMA48::init, wrong value of icntl[5] in MA48. Set to 1!";
      this->isInitialized = false;
      return false;
    }
    int N_RC_Max = std::max(this->numRows, this->numCols);
    int N_keep   = this->numRows + 5 * this->numCols + 4 * (this->numCols / ihlp) + 7;
    int N_iw     = this->numRows * 6 + this->numCols * 3;
    int N_w      = 3 * N_RC_Max + N_RC_Max;
    this->keep.resize(size_t(N_keep));
    this->iw.resize(size_t(N_iw));
    this->w.resize(size_t(N_w));
    this->isInitialized = true;
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename real>
  bool
  MA48<real>::factorize( real const ArrayA[] ) {
    // Check valid call:
    if (!this->isInitialized) {
      this->last_error   = "MA48::factorize, The function factorize can only be called after calling init!";
      this->isFactorized = false;
      return false;
    }

    // Copy Memory:
    std::copy( ArrayA, ArrayA+this->nnz, this->a.begin() );
    // Pivoting:
    this->ma48a(
      this->numRows,
      this->numCols,
      this->nnz,
      this->job,
      this->la,
      &this->a.front(),
      &this->irn.front(),
      &this->jcn.front(),
      &this->keep.front(),
      this->cntl,
      this->icntl,
      &this->iw.front(),
      this->info,
      this->rinfo
    );

    if ( this->info[3] > this->la ) {
      // Increase internal workspace:
      this->la = this->info[3];
      this->a.resize(size_t(this->la));
      this->irn.resize(size_t(this->la));
      this->jcn.resize(size_t(this->la));
      std::copy(ArrayA, ArrayA+this->nnz, this->a.begin());
      std::copy(this->i_Row_stored.begin(),
                this->i_Row_stored.end(),
                this->irn.begin());
      std::copy(this->j_Col_stored.begin(),
                this->j_Col_stored.end(),
                this->jcn.begin());

      this->ma48a(
        this->numRows,
        this->numCols,
        this->nnz,
        this->job,
        this->la,
        &this->a.front(),
        &this->irn.front(),
        &this->jcn.front(),
        &this->keep.front(),
        this->cntl,
        this->icntl,
        &this->iw.front(),
        this->info,
        this->rinfo
      );
      // Errors::set_Warning("Workspace of MA48 to small! Correcting!");
    }

    load_error_string(this->info[0]);

    this->isFactorized = false;
    if ( this->info[0] != 0 ) return false;

    // Factorize:
    this->ma48b(
      this->numRows,
      this->numCols,
      this->nnz,
      this->job,
      this->la,
      &this->a.front(),
      &this->irn.front(),
      &this->jcn.front(),
      &this->keep.front(),
      this->cntl,
      this->icntl,
      &this->w.front(),
      &this->iw.front(),
      this->info,
      this->rinfo
    );

    load_error_string(this->info[0]);

    this->isFactorized = this->info[0] == 0;
    return this->isFactorized;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename real>
  bool
  MA48<real>::solve( real const RHS[], real X[], bool transposed ) const {
    // Check valid call:
    if (!this->isInitialized) {
      this->last_error = "MA48::solve, Can not solve uninitialized system!";
      return false;
    }
    if (!this->isFactorized) {
      this->last_error = "MA48::solve, Can not solve unfactorized system!";
      return false;
    }

    // Solve with right hand side:
    this->ma48c(
      this->numRows,
      this->numCols,
      transposed ? 1 : 0,
      this->job,
      this->la,
      &this->a.front(),
      &this->irn.front(),
      &this->keep.front(),
      this->cntl,
      this->icntl,
      RHS,
      X,
      this->ma48errors,
      &this->w.front(),
      &this->iw.front(),
      this->info
    );

    load_error_string(this->info[0]);
    return this->info[0] == 0;
  }

  template <typename real>
  bool
  MA48<real>::solve(
    int        nrhs,
    real const RHS[],
    int        ldRHS,
    real       X[],
    int        ldX
  ) const {
    bool ok = true;
    if ( RHS == X ) {
      std::vector<real> tmpRHS;
      tmpRHS.resize(size_t(this->numRows));
      for ( int j = 0; ok && j < nrhs; ++j ) {
        real const * RHSj = RHS + j * ldRHS;
        std::copy( RHSj, RHSj+this->numRows, tmpRHS.begin() );
        ok = this->solve( &tmpRHS.front(), X + j * ldX, false );
      }
    } else {
      for ( int j = 0; ok && j < nrhs; ++j )
        ok = this->solve(RHS + j * ldRHS, X + j * ldX, false);
    }
    return ok;
  }

  template <typename real>
  bool
  MA48<real>::solve_transposed(
    int        nrhs,
    real const RHS[],
    int        ldRHS,
    real       X[],
    int        ldX
  ) const {
    bool ok = true;
    if ( RHS == X ) {
      std::vector<real> tmpRHS;
      tmpRHS.resize(size_t(this->numCols));
      for ( int j = 0; ok && j < nrhs; ++j ) {
        real const * RHSj = RHS + j * ldRHS;
        std::copy( RHSj, RHSj+this->numCols, tmpRHS.begin() );
        ok = this->solve(&tmpRHS.front(), X + j * ldX, true);
      }
    } else {
      for ( int j = 0; ok && j < nrhs; ++j )
        ok = this->solve(RHS + j * ldRHS, X + j * ldX, true);
    }
    return ok;
  }

  template class MA48<double>;
  template class MA48<float>;

} // namespace lapack_wrapper
