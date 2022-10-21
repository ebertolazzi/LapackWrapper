/*
 * \file hsl.h
 * Header Definition for Fortran-Routines of MA57 and MA48.
 *
 * \author A. Huber and E.Bertolazzi
 * \since 13.11.2018
 */

#ifndef HSL_dot_HH
#define HSL_dot_HH

#include "../lapack_wrapper/lapack_wrapper.hh"

#ifndef HSL_F77NAME
#define HSL_F77NAME(A) A##_
#endif

#define HSL_extern extern "C"

namespace lapack_wrapper {

  using logical = bool;

  /*\
   |   __  __    _   _  _   _
   |  |  \/  |  / \ | || | / |
   |  | |\/| | / _ \| || |_| |
   |  | |  | |/ ___ \__   _| |
   |  |_|  |_/_/   \_\ |_| |_|
  \*/

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
  );

  HSL_extern
  void
  HSL_F77NAME(ma41i)(
    float   CNTL[10],
    integer ICNTL[20],
    integer KEEP[50]
  );

  HSL_extern
  void
  HSL_F77NAME(ma41a)(
    integer const * JOB,
    integer const * N,
    integer const * NZ,
    integer         IRN[],
    integer         JCN[],
    float           ASPK[],
    float           RHS[],
    float           COLSCA[],
    float           ROWSCA[],
    integer         KEEP[50],
    integer         IS[],
    integer const * MAXIS,
    float           S[],
    integer const * MAXS,
    float           CNTL[10],
    integer         ICNTL[20],
    integer         INFO[20],
    float           RINFO[20]
  );

  HSL_extern
  void
  HSL_F77NAME(ma41id)(
    double  CNTL[10],
    integer ICNTL[20],
    integer KEEP[50]
  );

  HSL_extern
  void
  HSL_F77NAME(ma41ad)(
    integer const * JOB,
    integer const * N,
    integer const * NZ,
    integer         IRN[],
    integer         JCN[],
    double          ASPK[],
    double          RHS[],
    double          COLSCA[],
    double          ROWSCA[],
    integer         KEEP[50],
    integer         IS[],
    integer const * MAXIS,
    double          S[],
    integer const * MAXS,
    double          CNTL[10],
    integer         ICNTL[20],
    integer         INFO[20],
    double          RINFO[20]
  );

  template <typename T>
  void
  ma41i( T CNTL[10], integer ICNTL[20], integer KEEP[50] );

  template <>
  inline
  void
  ma41i<float>(
    float   CNTL[10],
    integer ICNTL[20],
    integer KEEP[50]
  ) {
    HSL_F77NAME(ma41i)(CNTL, ICNTL, KEEP);
    ICNTL[7]  = 0; // do scaling
    ICNTL[9]  = 5; // iterative refinement
    ICNTL[10] = 1; // do analisys
  }

  template <>
  inline
  void
  ma41i<double>(
    double  CNTL[10],
    integer ICNTL[20],
    integer KEEP[50]
  ) {
    HSL_F77NAME(ma41id)(CNTL, ICNTL, KEEP);
    ICNTL[7]  = 0; // do scaling
    ICNTL[9]  = 5; // iterative refinement
    ICNTL[10] = 1; // do analisys
  }

  template <typename T>
  void
  ma41a(
    integer JOB,
    integer N,
    integer NZ,
    integer IRN[],
    integer JCN[],
    T       ASPK[],
    T       RHS[],
    T       COLSCA[],
    T       ROWSCA[],
    integer KEEP[50],
    integer IS[],
    integer MAXIS,
    T       S[],
    integer MAXS,
    T       CNTL[10],
    integer ICNTL[20],
    integer INFO[20],
    T       RINFO[20]
  );

  template <>
  inline
  void
  ma41a<float>(
    integer JOB,
    integer N,
    integer NZ,
    integer IRN[],
    integer JCN[],
    float   ASPK[],
    float   RHS[],
    float   COLSCA[],
    float   ROWSCA[],
    integer KEEP[50],
    integer IS[],
    integer MAXIS,
    float   S[],
    integer MAXS,
    float   CNTL[10],
    integer ICNTL[20],
    integer INFO[20],
    float   RINFO[20]
  ) {
    HSL_F77NAME(ma41a)(
      &JOB, &N, &NZ, IRN, JCN, ASPK,
      RHS, COLSCA, ROWSCA, KEEP, IS, &MAXIS, S, &MAXS,
      CNTL, ICNTL, INFO, RINFO
    );
  }

  template <>
  inline
  void
  ma41a<double>(
    integer JOB,
    integer N,
    integer NZ,
    integer IRN[],
    integer JCN[],
    double  ASPK[],
    double  RHS[],
    double  COLSCA[],
    double  ROWSCA[],
    integer KEEP[50],
    integer IS[],
    integer MAXIS,
    double  S[],
    integer MAXS,
    double  CNTL[10],
    integer ICNTL[20],
    integer INFO[20],
    double  RINFO[20]
  ) {
    HSL_F77NAME(ma41ad)(
      &JOB, &N, &NZ, IRN, JCN, ASPK,
      RHS, COLSCA, ROWSCA, KEEP, IS, &MAXIS, S, &MAXS,
      CNTL, ICNTL, INFO, RINFO
    );
  }

  /*\
   |   __  __    _   _  _    ___
   |  |  \/  |  / \ | || |  ( _ )
   |  | |\/| | / _ \| || |_ / _ \
   |  | |  | |/ ___ \__   _| (_) |
   |  |_|  |_/_/   \_\ |_|  \___/
  \*/
  /**
   * \brief ma48id_:
   * Sets the control parameters of MA48.
   *
   * \param[out] cntl  Control parameter for double.
   * \param[out] icntl Control parameter for int.
   *
   */
  HSL_extern
  void
  HSL_F77NAME(ma48id)( double cntl[5], integer icntl[20] );

  /**
   * \brief ma48ad_:
   * Symbolic decomposition and choice of the pivot elements.
   *
   * \param[in]     nrow  Number of lines.
   * \param[in]     ncol  Number of columns.
   * \param[in]     nz    Number of elements of the spar matrix.
   * \param[in]     job   Type of linear equation system solution (here \f$ job = 1 \f$).
   * \param[in]     la    Auxiliary dimension \f$ la \f$.
   * \param[in,out] a     Value array of the spar matrix.
   * \param[in,out] irn   row array of the sparse matrix.
   * \param[in,out] jcn   column array of the sparse matrix.
   * \param[out]    keep  Array of integers for pivot strategy.
   * \param[in]     cntl  Control parameter for double.
   * \param[in]     icntl Control parameter for int.
   * \param[out]    iw    Workspace for int.
   * \param[out]    info  Info-Array for int.
   * \param[out]    rinfo Info-Array for double.
   *
   */
  HSL_extern
  void
  HSL_F77NAME(ma48ad)(
    integer const * nrow,
    integer const * ncol,
    integer const * nz,
    integer const * job,
    integer const * la,
    double  const   a[],
    integer         irn[],
    integer         jcn[],
    integer         keep[],
    double  const   cntl[10],
    integer const   icntl[20],
    integer         iw[],
    integer         info[20],
    double          rinfo[10]
  );

  /**
   * \brief ma48bd_:
   * Factorizes the savings matrix.
   *
   * \param[in]     nrow  Number of lines.
   * \param[in]     ncol  Number of columns.
   * \param[in]     nz    Number of elements of the spar matrix.
   * \param[in]     job   Type of linear equation system solution (here \f$ job = 1 \f$).
   * \param[in]     la    Auxyliary dimension \f$ la \f$.
   * \param[in,out] a     Value array of the spar matrix.
   * \param[in,out] irn   Row array of the sparse matrix.
   * \param[in,out] jcn   Column array of the sparse matrix.
   * \param[in]     keep  Array of integers for pivot strategy.
   * \param[in]     cntl  Control parameter for double.
   * \param[in]     icntl Control parameter for int.
   * \param[out]    w     Workspace for double.
   * \param[out]    iw    Workspace for int.
   * \param[out]    info  Info-Array for int.
   * \param[out]    rinfo Info-Array for double.
   *
   */
  HSL_extern
  void
  HSL_F77NAME(ma48bd)(
    integer const * nrow,
    integer const * ncol,
    integer const * nz,
    integer const * job,
    integer const * la,
    double          a[],
    integer         irn[],
    integer         jcn[],
    integer const   keep[],
    double  const   cntl[10],
    integer const   icntl[20],
    double          w[],
    integer         iw[],
    integer         info[20],
    double          rinfo[10]
  );

  /**
   * \brief ma48cd_:
   * Solves the linear equation system with respect to a right-hand side.
   *
   * \param[in]  nrow    Number of lines.
   * \param[in]  ncol    Number of columns.
   * \param[in]  itrans  Indicates whether the spar matrix is ​​transposed.
   * \param[in]  job     Type of linear equation system solution (here \ f $ job = 1 \ f $).
   * \param[in]  la      Auxyliary dimension \f$ la \f$.
   * \param[in]  a       Value array of the spar matrix.
   * \param[in]  irn     Row array of the spar matrix.
   * \param[in]  keep    Array of integers for pivot strategy.
   * \param[in]  cntl    Control parameter for double.
   * \param[in]  icntl   Control parameter for int.
   * \param[in]  rhs     Right side of the linear system of equations.
   * \param[out] x       Solution of the linear system of equations.
   * \param[out] errors  Errorarray of the solver.
   * \param[in]  w       Workspace for double.
   * \param[in]  iw      Workspace for int.
   * \param[out] info    Info-Array for int.
   *
   */
  HSL_extern
  void
  HSL_F77NAME(ma48cd)(
    integer const * nrow,
    integer const * ncol,
    integer const * itrans,
    integer const * job,
    integer const * la,
    double  const   a[],
    integer const   irn[],
    integer const   keep[],
    double  const   cntl[10],
    integer const   icntl[20],
    double  const   rhs[],
    double          x[],
    double          errors[3],
    double  const   w[],
    integer const   iw[],
    integer         info[20]
  );

  /**
   * \brief ma48i_:
   * Sets the control parameters of MA48.
   *
   * \param[out] cntl  Control parameter for double.
   * \param[out] icntl Control parameter for int.
   *
   */
  HSL_extern
  void
  HSL_F77NAME(ma48i)( float cntl[5], integer icntl[20] );

  /**
   * \brief ma48a_:
   * Symbolic decomposition and choice of the pivot elements.
   *
   * \param[in]     nrow  Number of lines.
   * \param[in]     ncol  Number of columns.
   * \param[in]     nz    Number of elements of the spar matrix.
   * \param[in]     job   Type of linear equation system solution (here \f$ job = 1 \f$).
   * \param[in]     la    Auxyliary dimension \f$ la \f$.
   * \param[in,out] a     Value array of the spar matrix.
   * \param[in,out] irn   Row array of the sparse matrix.
   * \param[in,out] jcn   Column array of the sparse matrix.
   * \param[out]    keep  Array of integers for pivot strategy.
   * \param[in]     cntl  Control parameter for double.
   * \param[in]     icntl Control parameter for int.
   * \param[out]    iw    Workspace for int.
   * \param[out]    info  Info-Array for int.
   * \param[out]    rinfo Info-Array for double.
   *
   */
  HSL_extern
  void
  HSL_F77NAME(ma48a)(
    integer const * nrow,
    integer const * ncol,
    integer const * nz,
    integer const * job,
    integer const * la,
    float           a[],
    integer         irn[],
    integer         jcn[],
    integer         keep[],
    float   const   cntl[10],
    integer const   icntl[20],
    integer         iw[],
    integer         info[20],
    float           rinfo[10]
  );

  /**
   * \brief ma48b_:
   * Factorizes the savings matrix.
   *
   * \param[in]     nrow  Number of lines.
   * \param[in]     ncol  Number of columns.
   * \param[in]     nz    Number of elements of the spar matrix.
   * \param[in]     job   Type of linear equation system solution (here \f$ job = 1 \f$).
   * \param[in]     la    Auxyliary dimension \f$ la \f$.
   * \param[in,out] a     Value array of the spar matrix.
   * \param[in,out] irn   Row array of the sparse matrix.
   * \param[in,out] jcn   Column array of the sparse matrix.
   * \param[in]     keep  Array of integers for pivot strategy.
   * \param[in]     cntl  Control parameter for double.
   * \param[in]     icntl Control parameter for int.
   * \param[out]    w     Workspace for double.
   * \param[out]    iw    Workspace for int.
   * \param[out]    info  Info-Array for int.
   * \param[out]    rinfo Info-Array for double.
   *
   */
  HSL_extern
  void
  HSL_F77NAME(ma48b)(
    integer const * nrow,
    integer const * ncol,
    integer const * nz,
    integer const * job,
    integer const * la,
    float           a[],
    integer         irn[],
    integer         jcn[],
    integer const   keep[],
    float   const   cntl[10],
    integer const   icntl[20],
    float           w[],
    integer         iw[],
    integer         info[20],
    float           rinfo[10]
  );

  /**
   * \brief ma48c_:
   * Solves the linear equation system with respect to a right-hand side.
   *
   * \param[in]  nrow   Number of lines.
   * \param[in]  ncol   Number of columns.
   * \param[in]  itrans Indicates whether the spar matrix is ​​transposed.
   * \param[in]  job    Type of linear equation system solution (here \f$ job = 1 \f$).
   * \param[in]  la     Auxyliary dimension \f$ la \f$.
   * \param[in]  a      Value array of the spar matrix.
   * \param[in]  line   Array of the spar matrix.
   * \param[in]  keep   Array of integers for pivot strategy.
   * \param[in]  cntl   Control parameter for double.
   * \param[in]  icntl  Control parameter for int.
   * \param[in]  rhs    Right side of the linear system of equations.
   * \param[out] x      Solution of the linear system of equations.
   * \param[out] errors Error array of the solver.
   * \param[in]  w      Workspace for double.
   * \param[in]  iw     Workspace for int.
   * \param[out] info   Info-Array for int.
   *
   */
  HSL_extern
  void
  HSL_F77NAME(ma48c)(
    integer const * nrow,
    integer const * ncol,
    integer const * itrans,
    integer const * job,
    integer const * la,
    float   const   a[],
    integer const   irn[],
    integer const   keep[],
    float   const   cntl[10],
    integer const   icntl[20],
    float   const   rhs[],
    float           x[],
    float           errors[3],
    float   const   w[],
    integer const   iw[],
    integer         info[20]
  );

  template <typename T>
  void
  ma48i( T cntl[5], integer icntl[20] );

  template <>
  inline
  void
  ma48i<float>( float cntl[5], integer icntl[20] ) {
    HSL_F77NAME(ma48i)( cntl, icntl );
  }

  template <>
  inline
  void
  ma48i<double>( double cntl[5], integer icntl[20] ) {
    HSL_F77NAME(ma48id)( cntl, icntl );
  }

  template <typename T>
  void
  ma48a(
    integer       nrow,
    integer       ncol,
    integer       nz,
    integer       job,
    integer       la,
    T             a[],
    integer       irn[],
    integer       jcn[],
    integer       keep[],
    T       const cntl[10],
    integer const icntl[20],
    integer       iw[],
    integer       info[20],
    T             rinfo[10]
  );

  template <>
  inline
  void
  ma48a<float>(
    integer       nrow,
    integer       ncol,
    integer       nz,
    integer       job,
    integer       la,
    float         a[],
    integer       irn[],
    integer       jcn[],
    integer       keep[],
    float   const cntl[10],
    integer const icntl[20],
    integer       iw[],
    integer       info[20],
    float         rinfo[10]
  ) {
    HSL_F77NAME(ma48a)(
      &nrow, &ncol, &nz, &job, &la,
      a, irn, jcn, keep, cntl, icntl, iw, info, rinfo
    );
  }

  template <>
  inline
  void
  ma48a<double>(
    integer       nrow,
    integer       ncol,
    integer       nz,
    integer       job,
    integer       la,
    double        a[],
    integer       irn[],
    integer       jcn[],
    integer       keep[],
    double  const cntl[10],
    integer const icntl[20],
    integer       iw[],
    integer       info[20],
    double        rinfo[10]
  ) {
    HSL_F77NAME(ma48ad)(
      &nrow, &ncol, &nz, &job, &la,
      a, irn, jcn, keep, cntl, icntl, iw, info, rinfo
    );
  }

  template <typename T>
  void
  ma48b(
    integer       nrow,
    integer       ncol,
    integer       nz,
    integer       job,
    integer       la,
    T              a[],
    integer       irn[],
    integer       jcn[],
    integer const keep[],
    T       const cntl[10],
    integer const icntl[20],
    T             w[],
    integer       iw[],
    integer       info[20],
    T             rinfo[10]
  );

  template <>
  inline
  void
  ma48b<float>(
    integer       nrow,
    integer       ncol,
    integer       nz,
    integer       job,
    integer       la,
    float         a[],
    integer       irn[],
    integer       jcn[],
    integer const keep[],
    float   const cntl[10],
    integer const icntl[20],
    float         w[],
    integer       iw[],
    integer       info[20],
    float         rinfo[10]
  ) {
    HSL_F77NAME(ma48b)(
      &nrow, &ncol, &nz, &job, &la,
      a, irn, jcn, keep, cntl, icntl, w, iw, info, rinfo
    );
  }

  template <>
  inline
  void
  ma48b<double>(
    integer       nrow,
    integer       ncol,
    integer       nz,
    integer       job,
    integer       la,
    double        a[],
    integer       irn[],
    integer       jcn[],
    integer const keep[],
    double  const cntl[10],
    integer const icntl[20],
    double        w[],
    integer       iw[],
    integer       info[20],
    double        rinfo[10]
  ) {
    HSL_F77NAME(ma48bd)(
      &nrow, &ncol, &nz, &job, &la,
      a, irn, jcn, keep, cntl, icntl, w, iw, info, rinfo
    );
  }

  template <typename T>
  void
  ma48c(
    integer       nrow,
    integer       ncol,
    integer       itrans,
    integer       job,
    integer       la,
    T       const a[],
    integer const irn[],
    integer const keep[],
    T       const cntl[10],
    integer const icntl[20],
    T       const rhs[],
    T             x[],
    T             errors[3],
    T       const w[],
    integer const iw[],
    integer       info[20]
  );

  template <>
  inline
  void
  ma48c<float>(
    integer       nrow,
    integer       ncol,
    integer       itrans,
    integer       job,
    integer       la,
    float   const a[],
    integer const irn[],
    integer const keep[],
    float   const cntl[10],
    integer const icntl[20],
    float   const rhs[],
    float         x[],
    float         errors[3],
    float   const w[],
    integer const iw[],
    integer       info[20]
  ) {
    HSL_F77NAME(ma48c)(
      &nrow, &ncol, &itrans, &job, &la,
      a, irn, keep, cntl, icntl, rhs, x, errors,
      w, iw, info
    );
  }

  template <>
  inline
  void
  ma48c<double>(
    integer       nrow,
    integer       ncol,
    integer       itrans,
    integer       job,
    integer       la,
    double  const a[],
    integer const irn[],
    integer const keep[],
    double  const cntl[10],
    integer const icntl[20],
    double  const rhs[],
    double        x[],
    double        errors[3],
    double  const w[],
    integer const iw[],
    integer       info[20]
  ) {
    HSL_F77NAME(ma48cd)(
      &nrow, &ncol, &itrans, &job, &la,
      a, irn, keep, cntl, icntl, rhs, x, errors,
      w, iw, info
    );
  }

  /*\
   |   __  __    _    ____ _____
   |  |  \/  |  / \  | ___|___  |
   |  | |\/| | / _ \ |___ \  / /
   |  | |  | |/ ___ \ ___) |/ /
   |  |_|  |_/_/   \_\____//_/
  \*/

  /**
   * \brief ma57id:
   *        Get the control parameters of MA57.
   *
   * \param[out] cntl   Control parameter for double.
   * \param[out] icntl  Control parameter for int.
   *
   */
  HSL_extern
  void
  HSL_F77NAME(ma57id)( double cntl[5], integer icntl[20] );

  /**
   * \brief ma57ad_:
   *        Symbolic decomposition and choice of the pivot elements
   *
   * \param[in]  n     order of the matrix.
   * \param[in]  nz    Number of sparse elements.
   * \param[in]  irn   line array of the sparse matrix.
   * \param[in]  jcn   column array of the sparse matrix.
   * \param[in]  lkeep Dimension of the keep array (\f$ lkeep \geq 5 n + nz + \max (n, nz) +42 \f$).
   * \param[out] keep  array of integers for pivot strategy.
   * \param[out] iw    Workspace of integers (\f$ 5n \f$).
   * \param[in]  icntl Control parameter for int.
   * \param[out] info  Info-Array for int.
   * \param[out] rinfo info array for double.
   *
   */
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
    double          rinfo[20]
  );

  /**
   * \brief ma57bd:
   * Factorizes the savings matrix.
   *
   * \param[in]  n      Order of the matrix.
   * \param[in]  nz     Number of sparse elements.
   * \param[in]  a      Value array of the sparse matrix.
   * \param[out] fact   Factors of the saving matrix for double.
   * \param[in]  lfact  Dimension of the factor array for double.
   * \param[out] ifact  Factors of the spar matrix for int.
   * \param[in]  lifact Dimension of the factor array for int.
   * \param[in]  lkeep  Dimension of the keep array (\f$ lkeep \geq 5 n + nz + max (n, nz) +42 \f$).
   * \param[in]  keep   array of integers for pivot strategy.
   * \param[out] iw     Workspace of integers (\f$ 5n \f$).
   * \param[in]  icntl  Control parameter for int.
   * \param[in]  cntl   Control parameter for double.
   * \param[out] info   Info-Array for int.
   * \param[out] rinfo  Info-Array for double.
   *
   */
  HSL_extern
  void
  HSL_F77NAME(ma57bd)(
    integer const * n,
    integer const * nz,
    double  const   a[],
    double          fact[],
    integer const * lfact,
    integer         ifact[],
    integer const * lifact,
    integer const * lkeep,
    integer const   keep[],
    integer         iw[],
    integer const   icntl[20],
    double  const   cntl[5],
    integer         info[40],
    double          rinfo[20]
  );

  /**
   * \brief ma57cd_:
   * Solves the linear equation system with respect to a right-hand side.
   *
   * \param[in]     job    Type of linear equation system solution (here \f$ job = 1 \f$).
   * \param[in]     n      Order of the matrix.
   * \param[in]     fact   Factors of the sparing matrix for double.
   * \param[in]     lfact  Dimension of the factor array for double.
   * \param[in]     ifact  Factors of the spar matrix for int.
   * \param[in]     lifact Dimension of the factor array for int.
   * \param[in]     nrhs   Number of right sides.
   * \param[in,out] rhs    Right side of the linear system of equations.
   * \param[in]     lrhs   Line dimension of the right side.
   * \param[out]    w      Workspace for double.
   * \param[in]     lw     Dimension of the workspace for double.
   * \param[out]    iw     Workspace for int.
   * \param[in]     icntl  Control parameter for int.
   * \param[out]    info   Info-Array for int.
   *
   */
  HSL_extern
  void
  HSL_F77NAME(ma57cd)(
    integer const * job,
    integer const * n,
    double  const   fact[],
    integer const * lfact,
    integer const   ifact[],
    integer const * lifact,
    integer const * nrhs,
    double          rhs[],
    integer const * lrhs,
    double          w[],
    integer const * lw,
    integer         iw[],
    integer const   icntl[20],
    integer         info[40]
  );

  /**
   * \brief ma57dd_:
   * Solves the linear equation system with respect to a right-hand side
   * with iterative refinement.
   *
   * \param[in]  job    Type of linear equation system solution.
   * \param[in]  n      Order of the matrix.
   * \param[in]  ne     Number of non-zero elements in the source matrix.
   * \param[in]  a      Value array of the source matrix.
   * \param[in]  irn    Index array of the lines of the source matrix.
   * \param[in]  jcn    Index array of the columns of the source matrix.
   * \param[in]  fact   factors of the sparing matrix for double.
   * \param[in]  lfact  Dimension of the factor array for double.
   * \param[in]  ifact  Factors of the spar matrix for int.
   * \param[in]  lifact Dimension of the factor array for int.
   * \param[in]  rhs    Right side of the linear system of equations.
   * \param[out] x      Solution of linear equation system.
   * \param[out] resid  Residuum of the solution of the linear system of equations.
   * \param[out] w      Workspace for double.
   * \param[out] iw     Workspace for int.
   * \param[in]  icntl  Control parameter for int.
   * \param[in]  cntl   Control parameter for double.
   * \param[out] info   Info-Array for int.
   * \param[out] rinfo  Info-Array for double.
   *
   */
  HSL_extern
  void
  HSL_F77NAME(ma57dd)(
    integer const * job,
    integer const * n,
    integer const * ne,
    double  const   a[],
    integer const   irn[],
    integer const   jcn[],
    double  const   fact[],
    integer const * lfact,
    integer const   ifact[],
    integer const * lifact,
    double  const   rhs[],
    double          x[],
    double          resid[],
    double          w[],
    integer         iw[],
    integer const   icntl[20],
    double  const   cntl[5],
    integer         info[40],
    double          rinfo[20]
  );

  /**
   * \brief ma57i_:
   * Sets the control parameters of MA57.
   *
   * \param[out] cntl  Control parameter for double.
   * \param[out] icntl Control parameter for int.
   *
   */
  HSL_extern
  void
  HSL_F77NAME(ma57i)( float cntl[5], integer icntl[20] );

  /**
   * \brief ma57a_:
   * Symbolic decomposition and choice of the pivot elements.
   *
   * \param[in]  n     Order of the matrix.
   * \param[in]  nz    Number of sparse elements.
   * \param[in]  line  Array of the spar matrix.
   * \param[in]  jcn   column array of the spar matrix.
   * \param[in]  lkeep Dimension of the keep array (\f$ lkeep \geq 5 n + nz + max (n, nz) +42 \f$).
   * \param[out] keep  Array of integers for pivot strategy.
   * \param[out] iw    Workspace of integers (\f$ 5n \f$).
   * \param[in]  icntl Control parameter for int.
   * \param[out] info  Info-Array for int.
   * \param[out] rinfo Info-Array for double.
   *
   */
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
    float           rinfo[20]
  );

  /**
   * \brief ma57b_:
   * Factorizes the savings matrix.
   *
   * \param[in]  n      Order of the matrix.
   * \param[in]  nz     Number of sparse elements.
   * \param[in]  a      Value array of the spar matrix.
   * \param[out] fact   Factors of the saving matrix for double.
   * \param[in]  lfact  Dimension of the factor array for double.
   * \param[out] ifact  Factors of the spar matrix for int.
   * \param[in]  lifact Dimension of the factor array for int.
   * \param[in]  lkeep  Dimension of the keep array (\f$ lkeep \geq 5 n + nz + max (n, nz) +42 \f$).
   * \param[in]  keep   Array of integers for pivot strategy.
   * \param[out] iw     Workspace of integers (\f$ 5n \f$).
   * \param[in]  icntl  Control parameter for int.
   * \param[in]  cntl   Control parameter for double.
   * \param[out] info   Info-Array for int.
   * \param[out] rinfo  Info-Array for double.
   *
   */
  HSL_extern
  void
  HSL_F77NAME(ma57b)(
    integer const * n,
    integer const * nz,
    float   const   a[],
    float           fact[],
    integer const * lfact,
    integer const   ifact[],
    integer const * lifact,
    integer const * lkeep,
    integer const   keep[],
    integer         iw[],
    integer const   icntl[20],
    float   const   cntl[5],
    integer         info[40],
    float           rinfo[20]
  );

  /**
   * \brief ma57c_:
   * Solves the linear equation system with respect to a right-hand side.
   *
   * \param[in]     job    Type of linear equation system solution (here \f$ job = 1 \f$).
   * \param[in]     n      Order of the matrix.
   * \param[in]     fact   Factors of the sparing matrix for double.
   * \param[in]     lfact  Dimension of the factor array for double.
   * \param[in]     ifact  Factors of the spar matrix for int.
   * \param[in]     lifact Dimension of the factor array for int.
   * \param[in]     nrhs   Number of right sides.
   * \param[in,out] rhs    Right side of the linear system of equations.
   * \param[in]     ldRhs  Right-sided row dimension.
   * \param[out]    w      Workspace for double.
   * \param[in]     lw     Dimension of the workspace for double.
   * \param[out]    iw     Workspace for int.
   * \param[in]     icntl  Control parameter for int.
   * \param[out]    info   Info-Array for int.
   *
   */
  HSL_extern
  void
  HSL_F77NAME(ma57c)(
    integer const * job,
    integer const * n,
    float   const   fact[],
    integer const * lfact,
    integer const   ifact[],
    integer const * lifact,
    integer const * nrhs,
    float           rhs[],
    integer const * ldRhs,
    float           w[],
    integer const * lw,
    integer         iw[],
    integer const   icntl[20],
    integer         info[40]
  );

  /**
   * \brief ma57d_:
   * Solves the linear equation system with respect to a right-hand side
   * with iterative refinement.
   *
   * \param[in]  job    Type of linear equation system solution.
   * \param[in]  n      Order of the matrix.
   * \param[in]  ne     Number of non-zero elements in the source matrix.
   * \param[in]  a      Value array of the source matrix.
   * \param[in]  irn    Index array of the lines of the source matrix.
   * \param[in]  jcn    Index array of the columns of the source matrix.
   * \param[in]  fact   Factors of the sparing matrix for double.
   * \param[in]  lfact  Dimension of the factor array for double.
   * \param[in]  ifact  Factors of the spar matrix for int.
   * \param[in]  lifact Dimension of the factor array for int.
   * \param[in]  rhs    Right side of the linear system of equations.
   * \param[out] x      Solution of linear equation system.
   * \param[out] resid  Residuum of the solution of the linear system of equations.
   * \param[out] w      Workspace for double.
   * \param[out] iw     Workspace for int.
   * \param[in]  icntl  Control parameter for int.
   * \param[in]  cntl   Control parameter for double.
   * \param[out] info   Info-Array for int.
   * \param[out] rinfo  Info-Array for double.
   *
   */
  HSL_extern
  void
  HSL_F77NAME(ma57d)(
    integer const * job,
    integer const * n,
    integer const * ne,
    float   const   a[],
    integer const   irn[],
    integer const   jcn[],
    float   const   fact[],
    integer const * lfact,
    integer const   ifact[],
    integer const * lifact,
    float   const   rhs[],
    float           x[],
    float         * resid,
    float           w[],
    integer         iw[],
    integer const   icntl[20],
    float   const   cntl[5],
    integer         info[40],
    float           rinfo[20]
  );

  template <typename T>
  void
  ma57i( T _cntl[], integer _icntl[] );

  template <>
  inline
  void
  ma57i<float>( float _cntl[], integer _icntl[] ) {
    HSL_F77NAME(ma57i)(_cntl, _icntl);
  }

  template <>
  inline
  void
  ma57i<double>( double _cntl[], integer _icntl[] ) {
    HSL_F77NAME(ma57id)(_cntl, _icntl);
  }

  template <typename T>
  void
  ma57a(
    integer       _n,
    integer       _nz,
    integer const _irn[],
    integer const _jcn[],
    integer       _lkeep,
    integer       _keep[],
    integer       _iw[],
    integer const _icntl[20],
    integer       _info[40],
    T            _rinfo[20]
  );

  template <>
  inline
  void
  ma57a<float>(
    integer       _n,
    integer       _nz,
    integer const _irn[],
    integer const _jcn[],
    integer       _lkeep,
    integer       _keep[],
    integer       _iw[],
    integer const _icntl[20],
    integer       _info[40],
    float        _rinfo[20]
  ) {
    HSL_F77NAME(ma57a)(
      &_n, &_nz, _irn, _jcn, &_lkeep, _keep, _iw, _icntl, _info, _rinfo
    );
  }

  template <>
  inline
  void
  ma57a<double>(
    integer       _n,
    integer       _nz,
    integer const _irn[],
    integer const _jcn[],
    integer       _lkeep,
    integer       _keep[],
    integer       _iw[],
    integer const _icntl[20],
    integer       _info[40],
    double        _rinfo[20]
  ) {
    HSL_F77NAME(ma57ad)(
      &_n, &_nz, _irn, _jcn, &_lkeep, _keep, _iw, _icntl, _info, _rinfo
    );
  }

  template <typename T>
  void
  ma57b(
    integer       _n,
    integer       _nz,
    T       const _a[],
    T             _fact[],
    integer       _lfact,
    integer       _ifact[],
    integer       _lifact,
    integer       _lkeep,
    integer const _keep[],
    integer       _iw[],
    integer const _icntl[20],
    T       const _cntl[5],
    integer       _info[40],
    T            _rinfo[20]
  );

  template <>
  inline
  void
  ma57b<float>(
    integer       _n,
    integer       _nz,
    float   const _a[],
    float         _fact[],
    integer       _lfact,
    integer       _ifact[],
    integer       _lifact,
    integer       _lkeep,
    integer const _keep[],
    integer       _iw[],
    integer const _icntl[20],
    float   const _cntl[5],
    integer       _info[40],
    float         _rinfo[20]
  ) {
    HSL_F77NAME(ma57b)(
      &_n, &_nz, _a, _fact, &_lfact, _ifact, &_lifact, &_lkeep, _keep,
      _iw, _icntl, _cntl, _info, _rinfo
    );
  }

  template <>
  inline
  void
  ma57b<double>(
    integer       _n,
    integer       _nz,
    double  const _a[],
    double        _fact[],
    integer       _lfact,
    integer       _ifact[],
    integer       _lifact,
    integer       _lkeep,
    integer const _keep[],
    integer       _iw[],
    integer const _icntl[20],
    double  const _cntl[5],
    integer       _info[40],
    double        _rinfo[20]
  ) {
    HSL_F77NAME(ma57bd)(
      &_n, &_nz, _a, _fact, &_lfact, _ifact, &_lifact, &_lkeep, _keep,
      _iw, _icntl, _cntl, _info, _rinfo
    );
  }

  template <typename T>
  void
  ma57c(
    integer       _job,
    integer       _n,
    T       const _fact[],
    integer       _lfact,
    integer const _ifact[],
    integer       _lifact,
    integer       _nrhs,
    T             _rhs[],
    integer       _lrhs,
    T             _w[],
    integer       _lw,
    integer       _iw[],
    integer const _icntl[20],
    integer       _info[40]
  );

  template <>
  inline
  void
  ma57c<float>(
    integer       _job,
    integer       _n,
    float   const _fact[],
    integer       _lfact,
    integer const _ifact[],
    integer       _lifact,
    integer       _nrhs,
    float         _rhs[],
    integer       _lrhs,
    float         _w[],
    integer       _lw,
    integer       _iw[],
    integer const _icntl[20],
    integer       _info[40]
  ) {
    HSL_F77NAME(ma57c)(
      &_job, &_n, _fact, &_lfact, _ifact, &_lifact,
      &_nrhs, _rhs, &_lrhs, _w, &_lw, _iw, _icntl, _info
    );
  }

  template <>
  inline
  void
  ma57c<double>(
    integer       _job,
    integer       _n,
    double  const _fact[],
    integer       _lfact,
    integer const _ifact[],
    integer       _lifact,
    integer       _nrhs,
    double        _rhs[],
    integer       _lrhs,
    double        _w[],
    integer       _lw,
    integer       _iw[],
    integer const _icntl[20],
    integer       _info[40]
  ) {
    HSL_F77NAME(ma57cd)(
      &_job, &_n, _fact, &_lfact, _ifact, &_lifact, &_nrhs, _rhs, &_lrhs,
      _w, &_lw, _iw, _icntl, _info
    );
  }

  template <typename T>
  void
  ma57d(
    integer       _job,
    integer       _n,
    integer       _ne,
    T       const _a[],
    integer const _irn[],
    integer const _jcn[],
    T       const _fact[],
    integer       _lfact,
    integer const _ifact[],
    integer       _lifact,
    T       const _rhs[],
    T             _x[],
    T             _resid[],
    T             _w[],
    integer       _iw[],
    integer const _icntl[20],
    T       const _cntl[5],
    integer       _info[40],
    T             _rinfo[20]
  );

  template <>
  inline
  void
  ma57d<float>(
    integer       _job,
    integer       _n,
    integer       _ne,
    float   const _a[],
    integer const _irn[],
    integer const _jcn[],
    float   const _fact[],
    integer       _lfact,
    integer const _ifact[],
    integer       _lifact,
    float   const _rhs[],
    float         _x[],
    float         _resid[],
    float         _w[],
    integer       _iw[],
    integer const _icntl[20],
    float   const _cntl[5],
    integer       _info[40],
    float         _rinfo[20]
  ) {
    HSL_F77NAME(ma57d)(
      &_job, &_n, &_ne, _a, _irn, _jcn, _fact, &_lfact, _ifact, &_lifact,
      _rhs, _x, _resid, _w, _iw, _icntl, _cntl, _info, _rinfo
    );
  }

  template <>
  inline
  void
  ma57d<double>(
    integer       _job,
    integer       _n,
    integer       _ne,
    double  const _a[],
    integer const _irn[],
    integer const _jcn[],
    double  const _fact[],
    integer       _lfact,
    integer const _ifact[],
    integer       _lifact,
    double  const _rhs[],
    double        _x[],
    double        _resid[],
    double        _w[],
    integer       _iw[],
    integer const _icntl[20],
    double  const _cntl[5],
    integer       _info[40],
    double        _rinfo[20]
  ) {
    HSL_F77NAME(ma57dd)(
      &_job, &_n, &_ne, _a, _irn, _jcn, _fact, &_lfact, _ifact, &_lifact,
      _rhs, _x, _resid, _w, _iw, _icntl, _cntl, _info, _rinfo
    );
  }

} // end namesapace


#endif // HSL_dot_HH
