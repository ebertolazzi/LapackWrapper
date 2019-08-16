/* (C) 2018, Andreas Huber, UniBw München LRT 1, Germany. All rights reserved. */

/**
 * @file hsl.h
 * Header Definition for Fortran-Routines of MA57 and MA48.
 *
 * @author A. Huber and E.Bertolazzi
 * @since 13.11.2018
 */

#ifndef HSL_H
#define HSL_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef HSL_F77NAME
#define HSL_F77NAME(A) A##_
#endif

/**
 * \brief ma57id:
 *        Get the control parameters of MA57.
 *
 * \param[out] cntl   Control parameter for double.
 * \param[out] icntl  Control parameter for int.
 *
 */
extern
void
HSL_F77NAME(ma57id)( double cntl[5], int icntl[20] );

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
extern
void
HSL_F77NAME(ma57ad)(
  int const * n,
  int const * nz,
  int const   irn[],
  int const   jcn[],
  int const * lkeep,
  int         keep[],
  int         iw[],
  int const   icntl[20],
  int         info[40],
  double      rinfo[20]
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
extern
void
HSL_F77NAME(ma57bd)(
  int const *  n,
  int const *  nz,
  double const a[],
  double       fact[],
  int const *  lfact,
  int          ifact[],
  int const *  lifact,
  int const *  lkeep,
  int const    keep[],
  int          iw[],
  int const    icntl[20],
  double const cntl[5],
  int          info[40],
  double       rinfo[20]
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
extern
void
HSL_F77NAME(ma57cd)(
  int const *  job,
  int const *  n,
  double const fact[],
  int const *  lfact,
  int const    ifact[],
  int const *  lifact,
  int const *  nrhs,
  double       rhs[],
  int const *  lrhs,
  double       w[],
  int const *  lw,
  int          iw[],
  int const    icntl[20],
  int          info[40]
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
extern
void
HSL_F77NAME(ma57dd)(
  int const *  job,
  int const *  n,
  int const *  ne,
  double const a[],
  int const    irn[],
  int const    jcn[],
  double const fact[],
  int const *  lfact,
  int const    ifact[],
  int const *  lifact,
  double const rhs[],
  double       x[],
  double       resid[],
  double       w[],
  int          iw[],
  int const    icntl[20],
  double const cntl[5],
  int          info[40],
  double       rinfo[20]
);

/**
 * \brief ma57i_:
 * Sets the control parameters of MA57.
 *
 * \param[out] cntl  Control parameter for double.
 * \param[out] icntl Control parameter for int.
 *
 */
extern
void
HSL_F77NAME(ma57i)( float cntl[5], int icntl[20] );

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
extern
void
HSL_F77NAME(ma57a)(
  int const * n,
  int const * nz,
  int const   irn[],
  int const   jcn[],
  int const * lkeep,
  int         keep[],
  int         iw[],
  int const   icntl[20],
  int         info[40],
  float       rinfo[20]
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
extern
void
HSL_F77NAME(ma57b)(
  int const * n,
  int const * nz,
  float const a[],
  float       fact[],
  int const * lfact,
  int const   ifact[],
  int const * lifact,
  int const * lkeep,
  int const   keep[],
  int         iw[],
  int const   icntl[20],
  float const cntl[5],
  int         info[40],
  float       rinfo[20]
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
extern
void
HSL_F77NAME(ma57c)(
  int const * job,
  int const * n,
  float const fact[],
  int const * lfact,
  int const   ifact[],
  int const * lifact,
  int const * nrhs,
  float       rhs[],
  int const * ldRhs,
  float       w[],
  int const * lw,
  int         iw[],
  int const   icntl[20],
  int         info[40]
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
extern
void
HSL_F77NAME(ma57d)(
  int const * job,
  int const * n,
  int const * ne,
  float const a[],
  int const   irn[],
  int const   jcn[],
  float const fact[],
  int const * lfact,
  int const   ifact[],
  int const * lifact,
  float const rhs[],
  float       x[],
  float *     resid,
  float       w[],
  int         iw[],
  int const   icntl[20],
  float const cntl[5],
  int         info[40],
  float       rinfo[20]
);

/**
 * \brief ma48id_:
 * Sets the control parameters of MA48.
 *
 * \param[out] cntl  Control parameter for double.
 * \param[out] icntl Control parameter for int.
 *
 */
extern
void
HSL_F77NAME(ma48id)( double cntl[5], int icntl[20] );

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
extern
void
HSL_F77NAME(ma48ad)(
  int const *  nrow,
  int const *  ncol,
  int const *  nz,
  int const *  job,
  int const *  la,
  double const a[],
  int          irn[],
  int          jcn[],
  int          keep[],
  double const cntl[10],
  int const    icntl[20],
  int          iw[],
  int          info[20],
  double       rinfo[10]
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
extern
void
HSL_F77NAME(ma48bd)(
  int const *  nrow,
  int const *  ncol,
  int const *  nz,
  int const *  job,
  int const *  la,
  double       a[],
  int          irn[],
  int          jcn[],
  int const    keep[],
  double const cntl[10],
  int const    icntl[20],
  double       w[],
  int          iw[],
  int          info[20],
  double       rinfo[10]
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
extern
void
HSL_F77NAME(ma48cd)(
  int const *  nrow,
  int const *  ncol,
  int const *  itrans,
  int const *  job,
  int const *  la,
  double const a[],
  int const    irn[],
  int const    keep[],
  double const cntl[10],
  int const    icntl[20],
  double const rhs[],
  double       x[],
  double       errors[3],
  double const w[],
  int const    iw[],
  int          info[20]
);

/**
 * \brief ma48i_:
 * Sets the control parameters of MA48.
 *
 * \param[out] cntl  Control parameter for double.
 * \param[out] icntl Control parameter for int.
 *
 */
extern
void
HSL_F77NAME(ma48i)( float cntl[5], int icntl[20] );

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
extern
void
HSL_F77NAME(ma48a)(
  int const * nrow,
  int const * ncol,
  int const * nz,
  int const * job,
  int const * la,
  float       a[],
  int         irn[],
  int         jcn[],
  int         keep[],
  float const cntl[10],
  int const   icntl[20],
  int         iw[],
  int         info[20],
  float       rinfo[10]
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
extern
void
HSL_F77NAME(ma48b)(
  int const * nrow,
  int const * ncol,
  int const * nz,
  int const * job,
  int const * la,
  float       a[],
  int         irn[],
  int         jcn[],
  int const   keep[],
  float const cntl[10],
  int const   icntl[20],
  float       w[],
  int         iw[],
  int         info[20],
  float       rinfo[10]
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
extern
void
HSL_F77NAME(ma48c)(
  int const * nrow,
  int const * ncol,
  int const * itrans,
  int const * job,
  int const * la,
  float const a[],
  int const   irn[],
  int const   keep[],
  float const cntl[10],
  int const   icntl[20],
  float const rhs[],
  float       x[],
  float       errors[3],
  float const w[],
  int const   iw[],
  int         info[20]
);

#ifdef __cplusplus
}
#endif

#endif // HSL_H
