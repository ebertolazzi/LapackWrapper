/**
 * \file hslsolver.cc
 * Header definitions of the class HSLsolver.
 *
 * \author A. Huber and E.Bertolazzi
 * \since  13.11.2018
 */

#include "hsl.h"
#include <iostream>
#include <iomanip>
#include <stdexcept>

#ifndef HSL_ERROR
  #define HSL_ERROR(MSG) {                             \
    std::cerr << "in file: " << __FILE__ << "\nline: " \
              << __LINE__ << '\n' << MSG << '\n';      \
    std::exit(1);                                      \
  }
#endif

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wmissing-noreturn"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wmissing-noreturn"
#endif

namespace HSL {

  HSL_extern
  void
  HSL_F77NAME(mc47a)(
    integer const * N,
    integer const * NE,
    integer         PE[],
    integer         IW[],
    integer const * IWLEN,
    integer const * MP,
    integer const   INFO[]
  ) {
    HSL_ERROR("mc47a not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma41i)(
    single  CNTL[10],
    integer ICNTL[20],
    integer KEEP[50]
  ) {
    HSL_ERROR("ma41i not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma41a)(
    integer const * JOB,
    integer const * N,
    integer const * NZ,
    integer         IRN[],
    integer         JCN[],
    single          ASPK[],
    single          RHS[],
    single          COLSCA[],
    single          ROWSCA[],
    integer         KEEP[50],
    integer         IS[],
    integer const * MAXIS,
    single          S[],
    integer const * MAXS,
    single          CNTL[10],
    integer         ICNTL[20],
    integer         INFO[20],
    single          RINFO[20]
  ) {
    HSL_ERROR("ma41a not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma41id)(
    real    CNTL[10],
    integer ICNTL[20],
    integer KEEP[50]
  ) {
    HSL_ERROR("ma41id not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma41ad)(
    integer const * JOB,
    integer const * N,
    integer const * NZ,
    integer         IRN[],
    integer         JCN[],
    real            ASPK[],
    real            RHS[],
    real            COLSCA[],
    real            ROWSCA[],
    integer         KEEP[50],
    integer         IS[],
    integer const * MAXIS,
    real            S[],
    integer const * MAXS,
    real            CNTL[10],
    integer         ICNTL[20],
    integer         INFO[20],
    real            RINFO[20]
  ) {
    HSL_ERROR("ma41ad not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma48id)( real cntl[5], integer icntl[20] ) {
    HSL_ERROR("ma48id not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma48ad)(
    integer const * nrow,
    integer const * ncol,
    integer const * nz,
    integer const * job,
    integer const * la,
    real    const   a[],
    integer         irn[],
    integer         jcn[],
    integer         keep[],
    real    const   cntl[10],
    integer const   icntl[20],
    integer         iw[],
    integer         info[20],
    real            rinfo[10]
  ) {
    HSL_ERROR("ma48ad not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma48bd)(
    integer const * nrow,
    integer const * ncol,
    integer const * nz,
    integer const * job,
    integer const * la,
    real            a[],
    integer         irn[],
    integer         jcn[],
    integer const   keep[],
    real    const   cntl[10],
    integer const   icntl[20],
    real            w[],
    integer         iw[],
    integer         info[20],
    real            rinfo[10]
  ) {
    HSL_ERROR("ma48bd not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma48cd)(
    integer const * nrow,
    integer const * ncol,
    integer const * itrans,
    integer const * job,
    integer const * la,
    real    const   a[],
    integer const   irn[],
    integer const   keep[],
    real    const   cntl[10],
    integer const   icntl[20],
    real    const   rhs[],
    real            x[],
    real            errors[3],
    real    const   w[],
    integer const   iw[],
    integer         info[20]
  ) {
    HSL_ERROR("ma48cd not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma48i)( single cntl[5], integer icntl[20] ) {
    HSL_ERROR("ma48i not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma48a)(
    integer const * nrow,
    integer const * ncol,
    integer const * nz,
    integer const * job,
    integer const * la,
    single          a[],
    integer         irn[],
    integer         jcn[],
    integer         keep[],
    single  const   cntl[10],
    integer const   icntl[20],
    integer         iw[],
    integer         info[20],
    single          rinfo[10]
  ) {
    HSL_ERROR("ma48a not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma48b)(
    integer const * nrow,
    integer const * ncol,
    integer const * nz,
    integer const * job,
    integer const * la,
    single          a[],
    integer         irn[],
    integer         jcn[],
    integer const   keep[],
    single  const   cntl[10],
    integer const   icntl[20],
    single          w[],
    integer         iw[],
    integer         info[20],
    single          rinfo[10]
  ) {
    HSL_ERROR("ma48b not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma48c)(
    integer const * nrow,
    integer const * ncol,
    integer const * itrans,
    integer const * job,
    integer const * la,
    single  const   a[],
    integer const   irn[],
    integer const   keep[],
    single  const   cntl[10],
    integer const   icntl[20],
    single  const   rhs[],
    single          x[],
    single          errors[3],
    single  const   w[],
    integer const   iw[],
    integer         info[20]
  ) {
    HSL_ERROR("ma48c not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma57id)( real cntl[5], integer icntl[20] ) {
    HSL_ERROR("ma57id not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma57ad)(
    integer const * n,
    integer const * nz,
    integer const   irn[],
    integer const   jcn[],
    integer const * lkeep,
    integer         keep[],
    integer         iw[],
    integer const   icntl[20],
    integer         info[40],
    real            rinfo[20]
  ) {
    HSL_ERROR("ma57ad not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma57bd)(
    integer const * n,
    integer const * nz,
    real    const   a[],
    real            fact[],
    integer const * lfact,
    integer         ifact[],
    integer const * lifact,
    integer const * lkeep,
    integer const   keep[],
    integer         iw[],
    integer const   icntl[20],
    real    const   cntl[5],
    integer         info[40],
    real            rinfo[20]
  ) {
    HSL_ERROR("ma57bd not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma57cd)(
    integer const * job,
    integer const * n,
    real    const   fact[],
    integer const * lfact,
    integer const   ifact[],
    integer const * lifact,
    integer const * nrhs,
    real            rhs[],
    integer const * lrhs,
    real            w[],
    integer const * lw,
    integer         iw[],
    integer const   icntl[20],
    integer         info[40]
  ) {
    HSL_ERROR("ma57cd not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma57dd)(
    integer const * job,
    integer const * n,
    integer const * ne,
    real    const   a[],
    integer const   irn[],
    integer const   jcn[],
    real    const   fact[],
    integer const * lfact,
    integer const   ifact[],
    integer const * lifact,
    real    const   rhs[],
    real            x[],
    real            resid[],
    real            w[],
    integer         iw[],
    integer const   icntl[20],
    real    const   cntl[5],
    integer         info[40],
    real            rinfo[20]
  ) {
    HSL_ERROR("ma57dd not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma57i)( single cntl[5], integer icntl[20] ) {
    HSL_ERROR("ma57i not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma57a)(
    integer const * n,
    integer const * nz,
    integer const   irn[],
    integer const   jcn[],
    integer const * lkeep,
    integer         keep[],
    integer         iw[],
    integer const   icntl[20],
    integer         info[40],
    single          rinfo[20]
  ) {
    HSL_ERROR("ma57a not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma57b)(
    integer const * n,
    integer const * nz,
    single  const   a[],
    single          fact[],
    integer const * lfact,
    integer const   ifact[],
    integer const * lifact,
    integer const * lkeep,
    integer const   keep[],
    integer         iw[],
    integer const   icntl[20],
    single  const   cntl[5],
    integer         info[40],
    single          rinfo[20]
  ) {
    HSL_ERROR("ma57b not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma57c)(
    integer const * job,
    integer const * n,
    single  const   fact[],
    integer const * lfact,
    integer const   ifact[],
    integer const * lifact,
    integer const * nrhs,
    single          rhs[],
    integer const * ldRhs,
    single          w[],
    integer const * lw,
    integer         iw[],
    integer const   icntl[20],
    integer         info[40]
  ) {
    HSL_ERROR("ma57c not defined");
  }

  HSL_extern
  void
  HSL_F77NAME(ma57d)(
    integer const * job,
    integer const * n,
    integer const * ne,
    single  const   a[],
    integer const   irn[],
    integer const   jcn[],
    single  const   fact[],
    integer const * lfact,
    integer const   ifact[],
    integer const * lifact,
    single  const   rhs[],
    single          x[],
    single        * resid,
    single          w[],
    integer         iw[],
    integer const   icntl[20],
    single  const   cntl[5],
    integer         info[40],
    single          rinfo[20]
  ) {
    HSL_ERROR("ma57d not defined");
  }

} // namespace HSL

//
// EOF: hsl_solver.cc
//
