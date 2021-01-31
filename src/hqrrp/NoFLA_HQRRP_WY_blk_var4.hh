
#include <iostream>
#include <vector>
#include <lapack_wrapper/lapack_wrapper.hh>
#include <lapack_wrapper/lapack_wrapper++.hh>

namespace lapack_wrapper {

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  dgeqp4(
    integer & m,
    integer & n,
    double    A[],
    integer & lda,
    integer   jpvt[],
    double    tau[],
    double    work[],
    integer & lwork,
    integer & info
  );

  int
  NoFLA_HQRRP_WY_blk_var4(
    int      m_A,
    int      n_A,
    double * buff_A,
    int      ldim_A,
    int    * buff_jpvt,
    double * buff_tau,
    int      nb_alg,
    int      pp,
    int      panel_pivoting
  );
}
