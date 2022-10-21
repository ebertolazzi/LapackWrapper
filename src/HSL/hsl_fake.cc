/*
 * \file hslsolver.cc
 * Header definitions of the class HSLsolver.
 *
 * \author A. Huber and E.Bertolazzi
 * \since  13.11.2018
 */

#include "hsl.hh"
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

extern "C" {

  using integer = lapack_wrapper::integer;

  void
  HSL_F77NAME(mc47a)(
    integer const * /*N*/,
    integer const * /*NE*/,
    integer         /*PE*/[],
    integer         /*IW*/[],
    integer const * /*IWLEN*/,
    integer const * /*MP*/,
    integer const   /*INFO*/[]
  ) {
    HSL_ERROR("mc47a not defined");
  }

  void
  HSL_F77NAME(ma41i)(
    float   /*CNTL*/[10],
    integer /*ICNTL*/[20],
    integer /*KEEP*/[50]
  ) {
    HSL_ERROR("ma41i not defined");
  }

  void
  HSL_F77NAME(ma41a)(
    integer const * /*JOB*/,
    integer const * /*N*/,
    integer const * /*NZ*/,
    integer         /*IRN*/[],
    integer         /*JCN*/[],
    float           /*ASPK*/[],
    float           /*RHS*/[],
    float           /*COLSCA*/[],
    float           /*ROWSCA*/[],
    integer         /*KEEP*/[50],
    integer         /*IS*/[],
    integer const * /*MAXIS*/,
    float           /*S*/[],
    integer const * /*MAXS*/,
    float           /*CNTL*/[10],
    integer         /*ICNTL*/[20],
    integer         /*INFO*/[20],
    float           /*RINFO*/[20]
  ) {
    HSL_ERROR("ma41a not defined");
  }

  void
  HSL_F77NAME(ma41id)(
    double  /*CNTL*/[10],
    integer /*ICNTL*/[20],
    integer /*KEEP*/[50]
  ) {
    HSL_ERROR("ma41id not defined");
  }

  void
  HSL_F77NAME(ma41ad)(
    integer const * /*JOB*/,
    integer const * /*N*/,
    integer const * /*NZ*/,
    integer         /*IRN*/[],
    integer         /*JCN*/[],
    double          /*ASPK*/[],
    double          /*RHS*/[],
    double          /*COLSCA*/[],
    double          /*ROWSCA*/[],
    integer         /*KEEP*/[50],
    integer         /*IS*/[],
    integer const * /*MAXIS*/,
    double          /*S*/[],
    integer const * /*MAXS*/,
    double          /*CNTL*/[10],
    integer         /*ICNTL*/[20],
    integer         /*INFO*/[20],
    double          /*RINFO*/[20]
  ) {
    HSL_ERROR("ma41ad not defined");
  }

  void
  HSL_F77NAME(ma48id)(
    double  /*cntl*/[5],
    integer /*icntl*/[20]
  ) {
    HSL_ERROR("ma48id not defined");
  }

  void
  HSL_F77NAME(ma48ad)(
    integer const * /*nrow*/,
    integer const * /*ncol*/,
    integer const * /*nz*/,
    integer const * /*job*/,
    integer const * /*la*/,
    double  const   /*a*/[],
    integer         /*irn*/[],
    integer         /*jcn*/[],
    integer         /*keep*/[],
    double  const   /*cntl*/[10],
    integer const   /*icntl*/[20],
    integer         /*iw*/[],
    integer         /*info*/[20],
    double          /*rinfo*/[10]
  ) {
    HSL_ERROR("ma48ad not defined");
  }

  void
  HSL_F77NAME(ma48bd)(
    integer const * /*nrow*/,
    integer const * /*ncol*/,
    integer const * /*nz*/,
    integer const * /*job*/,
    integer const * /*la*/,
    double          /*a*/[],
    integer         /*irn*/[],
    integer         /*jcn*/[],
    integer const   /*keep*/[],
    double  const   /*cntl*/[10],
    integer const   /*icntl*/[20],
    double          /*w*/[],
    integer         /*iw*/[],
    integer         /*info*/[20],
    double          /*rinfo*/[10]
  ) {
    HSL_ERROR("ma48bd not defined");
  }

  void
  HSL_F77NAME(ma48cd)(
    integer const * /*nrow*/,
    integer const * /*ncol*/,
    integer const * /*itrans*/,
    integer const * /*job*/,
    integer const * /*la*/,
    double  const   /*a*/[],
    integer const   /*irn*/[],
    integer const   /*keep*/[],
    double  const   /*cntl*/[10],
    integer const   /*icntl*/[20],
    double  const   /*rhs*/[],
    double          /*x*/[],
    double          /*errors*/[3],
    double  const   /*w*/[],
    integer const   /*iw*/[],
    integer         /*info*/[20]
  ) {
    HSL_ERROR("ma48cd not defined");
  }

  void
  HSL_F77NAME(ma48i)(
    float   /*cntl*/[5],
    integer /*icntl*/[20]
  ) {
    HSL_ERROR("ma48i not defined");
  }

  void
  HSL_F77NAME(ma48a)(
    integer const * /*nrow*/,
    integer const * /*ncol*/,
    integer const * /*nz*/,
    integer const * /*job*/,
    integer const * /*la*/,
    float           /*a*/[],
    integer         /*irn*/[],
    integer         /*jcn*/[],
    integer         /*keep*/[],
    float   const   /*cntl*/[10],
    integer const   /*icntl*/[20],
    integer         /*iw*/[],
    integer         /*info*/[20],
    float           /*rinfo*/[10]
  ) {
    HSL_ERROR("ma48a not defined");
  }

  void
  HSL_F77NAME(ma48b)(
    integer const * /*nrow*/,
    integer const * /*ncol*/,
    integer const * /*nz*/,
    integer const * /*job*/,
    integer const * /*la*/,
    float           /*a*/[],
    integer         /*irn*/[],
    integer         /*jcn*/[],
    integer const   /*keep*/[],
    float   const   /*cntl*/[10],
    integer const   /*icntl*/[20],
    float           /*w*/[],
    integer         /*iw*/[],
    integer         /*info*/[20],
    float           /*rinfo*/[10]
  ) {
    HSL_ERROR("ma48b not defined");
  }

  void
  HSL_F77NAME(ma48c)(
    integer const * /*nrow*/,
    integer const * /*ncol*/,
    integer const * /*itrans*/,
    integer const * /*job*/,
    integer const * /*la*/,
    float   const   /*a*/[],
    integer const   /*irn*/[],
    integer const   /*keep*/[],
    float   const   /*cntl*/[10],
    integer const   /*icntl*/[20],
    float   const   /*rhs*/[],
    float           /*x*/[],
    float           /*errors*/[3],
    float   const   /*w*/[],
    integer const   /*iw*/[],
    integer         /*info*/[20]
  ) {
    HSL_ERROR("ma48c not defined");
  }

  void
  HSL_F77NAME(ma57id)(
    double  /*cntl*/[5],
    integer /*icntl*/[20]
  ) {
    HSL_ERROR("ma57id not defined");
  }

  void
  HSL_F77NAME(ma57ad)(
    integer const * /*n*/,
    integer const * /*nz*/,
    integer const   /*irn*/[],
    integer const   /*jcn*/[],
    integer const * /*lkeep*/,
    integer         /*keep*/[],
    integer         /*iw*/[],
    integer const   /*icntl*/[20],
    integer         /*info*/[40],
    double          /*rinfo*/[20]
  ) {
    HSL_ERROR("ma57ad not defined");
  }

  void
  HSL_F77NAME(ma57bd)(
    integer const * /*n*/,
    integer const * /*nz*/,
    double  const   /*a*/[],
    double          /*fact*/[],
    integer const * /*lfact*/,
    integer         /*ifact*/[],
    integer const * /*lifact*/,
    integer const * /*lkeep*/,
    integer const   /*keep*/[],
    integer         /*iw*/[],
    integer const   /*icntl*/[20],
    double  const   /*cntl*/[5],
    integer         /*info*/[40],
    double          /*rinfo*/[20]
  ) {
    HSL_ERROR("ma57bd not defined");
  }

  void
  HSL_F77NAME(ma57cd)(
    integer const * /*job*/,
    integer const * /*n*/,
    double  const   /*fact*/[],
    integer const * /*lfact*/,
    integer const   /*ifact*/[],
    integer const * /*lifact*/,
    integer const * /*nrhs*/,
    double          /*rhs*/[],
    integer const * /*lrhs*/,
    double          /*w*/[],
    integer const * /*lw*/,
    integer         /*iw*/[],
    integer const   /*icntl*/[20],
    integer         /*info*/[40]
  ) {
    HSL_ERROR("ma57cd not defined");
  }

  void
  HSL_F77NAME(ma57dd)(
    integer const * /*job*/,
    integer const * /*n*/,
    integer const * /*ne*/,
    double  const   /*a*/[],
    integer const   /*irn*/[],
    integer const   /*jcn*/[],
    double  const   /*fact*/[],
    integer const * /*lfact*/,
    integer const   /*ifact*/[],
    integer const * /*lifact*/,
    double  const   /*rhs*/[],
    double          /*x*/[],
    double          /*resid*/[],
    double          /*w*/[],
    integer         /*iw*/[],
    integer const   /*icntl*/[20],
    double  const   /*cntl*/[5],
    integer         /*info*/[40],
    double          /*rinfo*/[20]
  ) {
    HSL_ERROR("ma57dd not defined");
  }

  void
  HSL_F77NAME(ma57i)(
    float   /*cntl*/[5],
    integer /*icntl*/[20]
  ) {
    HSL_ERROR("ma57i not defined");
  }

  void
  HSL_F77NAME(ma57a)(
    integer const * /*n*/,
    integer const * /*nz*/,
    integer const   /*irn*/[],
    integer const   /*jcn*/[],
    integer const * /*lkeep*/,
    integer         /*keep*/[],
    integer         /*iw*/[],
    integer const   /*icntl*/[20],
    integer         /*info*/[40],
    float           /*rinfo*/[20]
  ) {
    HSL_ERROR("ma57a not defined");
  }

  void
  HSL_F77NAME(ma57b)(
    integer const * /*n*/,
    integer const * /*nz*/,
    float   const   /*a*/[],
    float           /*fact*/[],
    integer const * /*lfact*/,
    integer const   /*ifact*/[],
    integer const * /*lifact*/,
    integer const * /*lkeep*/,
    integer const   /*keep*/[],
    integer         /*iw*/[],
    integer const   /*icntl*/[20],
    float   const   /*cntl*/[5],
    integer         /*info*/[40],
    float           /*rinfo*/[20]
  ) {
    HSL_ERROR("ma57b not defined");
  }

  void
  HSL_F77NAME(ma57c)(
    integer const * /*job*/,
    integer const * /*n*/,
    float   const   /*fact*/[],
    integer const * /*lfact*/,
    integer const   /*ifact*/[],
    integer const * /*lifact*/,
    integer const * /*nrhs*/,
    float           /*rhs*/[],
    integer const * /*ldRhs*/,
    float           /*w*/[],
    integer const * /*lw*/,
    integer         /*iw*/[],
    integer const   /*icntl*/[20],
    integer         /*info*/[40]
  ) {
    HSL_ERROR("ma57c not defined");
  }

  void
  HSL_F77NAME(ma57d)(
    integer const * /*job*/,
    integer const * /*n*/,
    integer const * /*ne*/,
    float   const   /*a*/[],
    integer const   /*irn*/[],
    integer const   /*jcn*/[],
    float   const   /*fact*/[],
    integer const * /*lfact*/,
    integer const   /*ifact*/[],
    integer const * /*lifact*/,
    float   const   /*rhs*/[],
    float           /*x*/[],
    float         * /*resid*/,
    float           /*w*/[],
    integer         /*iw*/[],
    integer const   /*icntl*/[20],
    float   const   /*cntl*/[5],
    integer         /*info*/[40],
    float           /*rinfo*/[20]
  ) {
    HSL_ERROR("ma57d not defined");
  }

}

//
// EOF: hsl_fake.cc
//
