/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Sparse_tool  : DRIVER FOR TESTING THE TOOLKIT INTERFACING WITH MA41     |
 |                                                                          |
 |  date         : 2011, 18 July                                            |
 |  version      : 1.0.                                                     |
 |  file         : MA41.cc                                                  |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Dipartimento di Ingegneria Industriale                   |
 |                 Università degli Studi di Trento                         |
 |                 email: enrico.bertolazzi@unitn.it                        |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "ma41.hh"

namespace Sparse_tool {

  template <typename T>
  void
  MA41<T>::msg_info() const {
    fmt::print("Real    space needs for factors    = {}\n", INFO[2] );
    fmt::print("Integer space needs for factors    = {}\n", INFO[3] );
    fmt::print("Maximum frontal size               = {}\n", INFO[4] );
    fmt::print("Nodes in the tree                  = {}\n", INFO[5] );
    fmt::print("Minimum MAXIS                      = {}\n", INFO[6] );
    fmt::print("Minimum MAXS                       = {}\n", INFO[7] );
    fmt::print("Real    Space for LU factorizarion = {}\n", INFO[8] );
    fmt::print("Integer Space for LU factorizarion = {}\n", INFO[9] );
    fmt::print("Integer Space for LU factorizarion = {}\n", INFO[10] );
    fmt::print("Largest frontal matrix             = {}\n", INFO[11] );
  }

  template <typename T>
  void
  MA41<T>::msg_infor() const {
    fmt::print("NFLOPs for LU factorization    = {}\n", RINFO[0]);
    fmt::print("NFLOPs for assembly process    = {}\n", RINFO[1]);
    fmt::print("NFLOPs for elimination process = {}\n", RINFO[2]);
    if ( ICNTL[10] > 0 ) {
      fmt::print("Infinity norm of the matrix    = {}\n", RINFO[3]);
      fmt::print("Infinity norm of the solution  = {}\n", RINFO[4]);
      fmt::print("Norm of scaled residual        = {}\n", RINFO[5]);
    }
  }

  template <typename T>
  void
  MA41<T>::msg_error() const {
    if ( INFO[0] >= 0 ) return;
    fmt::print(stderr,"\n***** ERROR *****\n");
    fmt::print(stderr,"dim = {}\n", N);
    fmt::print(stderr,"nnz = {}\n", NE);
    switch ( INFO[0] ) {
    case -1:
      fmt::print(stderr,"Value of N is out of range N = {}\n", INFO[1]);
      break;
    case -2:
      fmt::print(stderr,"Value of NE is out of range NE = {}\n", INFO[1]);
      break;
    case -3:
      fmt::print(stderr,"JOB has wrong value or analisys was not performed prior to factorization. JOB = {}\n", INFO[1]);
      break;
    case -4:
      fmt::print(stderr,"Error in permutation error\n");
      break;
    case -5:
      fmt::print(stderr,"Not enought space to preprocess the input matrix\n");
      fmt::print(stderr,"MAXS should be increased to at least {}\n", INFO[1]);
      break;
    case -6:
      fmt::print(stderr,"The matrix is structurally singular\n");
      fmt::print(stderr,"Estimated rank = {}\n", INFO[1]);
      break;
    case -7:
      fmt::print(stderr,"Error from analysis\n");
      fmt::print(stderr,"MAXIS should be increased to at least {}\n", INFO[1]);
      break;
    case -8:
      fmt::print(stderr,"Error from numerical factorization\n");
      fmt::print(stderr,"MAXIS should be increased to at least {}\n", INFO[1]);
      break;
    case -9:
      fmt::print(stderr,"Error from numerical factorization\n");
      fmt::print(stderr,"MAXS should be increased to at least {}\n", INFO[1]);
      break;
    case -10:
      fmt::print(stderr,"The matrix is numerically singular\n");
      fmt::print(stderr,"Estimated rank = {}\n", INFO[1]);
      break;
    case -11:
      fmt::print(stderr,"Error from the solution phase\n");
      fmt::print(stderr,"MAXS should be increased to at least {}\n", INFO[1]);
      break;
    case -12:
      fmt::print(stderr,"Not enought space to postprocess the solution\n");
      fmt::print(stderr,"MAXS should be increased to at least {}\n", INFO[1]);
      break;
    }
    fmt::print(stderr,"***** ERROR *****\n");
  }

  template <typename T>
  void
  MA41<T>::load_matrix( integer const nr ) {
    if ( this->verbose ) fmt::print("ma41::load_matrix(...)\n");
    this->N  = nr;
    this->NE = integer(this->A.size());

    this->COLSCA . resize(this->N);
    this->ROWSCA . resize(this->N);
    this->IS     . resize(2 * this->NE + 11 * this->N + 1);
    this->S      . resize(2 * this->NE + 11 * this->N + 1);

    fill(this->COLSCA.begin(), this->COLSCA.end(), 0 );
    fill(this->ROWSCA.begin(), this->ROWSCA.end(), 0 );
    fill(this->IS.begin(),     this->IS.end(),     0 );
    fill(this->S.begin(),      this->S.end(),      0 );

    // set pars
    HSL::ma41i<real_type>(this->CNTL, this->ICNTL, this->KEEP);

    if ( this->verbose ) fmt::print("done\n");
  }

  template <typename T>
  void
  MA41<T>::symbfac() {
    if ( this->verbose ) fmt::print("ma41::symbfac() perform analisys\n");

    // analysis phase
    vector<integer> tmpIS(2 * NE + 12 * N + 1);
    HSL::ma41a<T>(
      ANALISYS, this->N, this->NE,
      &this->IRN.front(), &this->JCN.front(), &this->A.front(),
      nullptr, &this->COLSCA.front(), &this->ROWSCA.front(), this->KEEP,
      &tmpIS.front(), integer(tmpIS.size()),
      &this->S.front(), integer(this->S.size()),
      this->CNTL, this->ICNTL, this->INFO, this->RINFO
    );

    if ( this->verbose ) msg_error();

    IS.resize(10*INFO[6]);
    S.resize(10*INFO[7]);

    copy( tmpIS.begin(), tmpIS.end(), IS.begin() );

    if ( this->verbose )
      fmt::print("ma41_wrapper::symbfac() do factorization\n");

    // factorization phase
    // -------------------
    HSL::ma41a<real_type>(
      FACTORIZATION, this->N, this->NE,
      &this->IRN.front(), &this->JCN.front(), &this->A.front(),
      nullptr, &this->COLSCA.front(), &this->ROWSCA.front(), this->KEEP,
      &this->IS.front(), integer(this->IS.size()),
      &this->S.front(), integer(this->S.size()),
      this->CNTL, this->ICNTL, this->INFO, this->RINFO
    );

    if ( this->verbose ) {
      msg_info();
      msg_error();
      fmt::print("done\n");
    }

  }

  template <typename T>
  void
  MA41<T>::init(void) {
    this->IRN.clear();
    this->JCN.clear();
    this->A.clear();
  }

  template <typename T>
  void
  MA41<T>::insert(
    integer   i,
    integer   j,
    real_type a
  ) {
    UTILS_ASSERT(
      i >= 0 && i < this->N && j >= 0 && j < this->N,
      "Sparse_tool: MA41::insert( {}, {}, a ) out of matrix\n", i, j
    )
    this->IRN.push_back(i+1);
    this->JCN.push_back(j+1);
    this->A.push_back(a);
  }

  template <typename T>
  void
  MA41<T>::setup( integer nr, bool v ) {
    verbose = v;
    load_matrix(nr);
    symbfac();
  }

  template <typename T>
  int
  MA41<T>::solve( real_type RHS[], bool transpose ) {
    this->ICNTL[8] = transpose ? 0 : 1;
    HSL::ma41a<real_type>(
      SOLVE, this->N, this->NE,
      &this->IRN.front(), &this->JCN.front(), &this->A.front(),
      RHS, &this->COLSCA.front(), &this->ROWSCA.front(), this->KEEP,
      &this->IS.front(), integer(this->IS.size()),
      &this->S.front(), integer(this->S.size()),
      this->CNTL, this->ICNTL, this->INFO, this->RINFO
    );

    if ( this->verbose ) {
      msg_infor();
      msg_error();
    }
    return this->INFO[0];
  }

  template <typename T>
  int
  MA41<T>::solve( T const b[], T x[], bool transpose ) {
    std::copy_n( b, this->N, x );
    return this->solve( x, transpose );
  }

  template class MA41<float>;
  template class MA41<double>;

}
