/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Sparse_tool  : DRIVER FOR TESTING THE TOOLKIT INTERFACING WITH MA41     |
 |                                                                          |
 |  date         : 2011, 18 July                                            |
 |  version      : 1.0.                                                     |
 |  file         : MA48.cc                                                  |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Dipartimento di Ingegneria Industriale                   |
 |                 Universita` degli Studi di Trento                        |
 |                 email: enrico.bertolazzi@unitn.it                        |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "ma48.hh"

namespace Sparse_tool {

  template <typename T>
  void
  MA48<T>::msg_info() const {
    fmt::print( "Estimated minimum space needed for factors = {}\n", INFO[2]);
    fmt::print( "Estimated space needed for factor          = {}\n", INFO[3]);
    fmt::print( "Estimate of the rank of the matrix         = {}\n", INFO[4]);
    fmt::print( "Entries dropped from the data structure    = {}\n", INFO[5]);
    fmt::print( "Order of the largest non-triangular block  = {}\n", INFO[6]);
    fmt::print( "Sum of orders of non-triangular blocks     = {}\n", INFO[7]);
    fmt::print( "# of entries in the non-triangular blocks  = {}\n", INFO[8]);
    fmt::print( "structural rank of the matrix.             = {}\n", INFO[9]);
    fmt::print( "# of multiple entries in the input matrix  = {}\n", INFO[10]);
    fmt::print( "# of entries with out-of-range indices     = {}\n", INFO[11]);
  }

  template <typename T>
  void
  MA48<T>::msg_infor() const {
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
  MA48<T>::msg_error() const {
    if ( INFO[0] >= 0 ) return;
    fmt::print(stderr,"***** ERROR *****\n");
    fmt::print(stderr,"nRow = {}\n", this->nRow);
    fmt::print(stderr,"nCol = {}\n", this->nCol);
    fmt::print(stderr,"nnz  = {}\n", this->nnz);
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
      fmt::print(stderr,"MAXS should be increased to at least {}\n",INFO[1]);
      break;
    case -6:
      fmt::print(stderr,"The matrix is structurally singular\n");
      fmt::print(stderr,"Estimated rank = {}\n",INFO[1]);
      break;
    case -7:
      fmt::print(stderr,"Error from analysis\n");
      fmt::print(stderr,"MAXIS should be increased to at least {}\n",INFO[1]);
      break;
    case -8:
      fmt::print(stderr,"Error from numerical factorization\n");
      fmt::print(stderr,"MAXIS should be increased to at least {}\n",INFO[1]);
      break;
    case -9:
      fmt::print(stderr,"Error from numerical factorization\n");
      fmt::print(stderr,"MAXS should be increased to at least {}\n", INFO[1]);
      break;
    case -10:
      fmt::print(stderr,"The matrix is numerically singular\n");
      fmt::print(stderr,"Estimated rank = {}\n",INFO[1]);
      break;
    case -11:
      fmt::print(stderr,"Error from the solution phase\n");
      fmt::print(stderr,"MAXS should be increased to at least {}\n",INFO[1]);
      break;
    case -12:
      fmt::print(stderr,"Not enought space to postprocess the solution\n");
      fmt::print(stderr,"MAXS should be increased to at least {}\n",INFO[1]);
      break;
    }
    fmt::print(stderr,"***** ERROR *****\n");
  }

  template <typename T>
  int
  MA48<T>::factorize() {
    if ( this->verbose ) fmt::print("ma48::factorize() perform analisys\n");

    // analysis phase
    A_work.resize( 3*A.size() );
    std::fill( A_work.begin(), A_work.end(), 0 );
    std::copy( A.begin(), A.end(), A_work.begin() );

    IRN_work.resize( 3*A.size() );
    std::fill( IRN_work.begin(), IRN_work.end(), 0 );
    std::copy( IRN.begin(), IRN.end(), IRN_work.begin() );

    JCN_work.resize( 3*A.size() );
    std::fill( JCN_work.begin(), JCN_work.end(), 0 );
    std::copy( JCN.begin(), JCN.end(), JCN_work.begin() );

    HSL::ma48a<T>(
      this->nRow, this->nCol, this->nnz, 1,
      integer(A_work.size()), &A_work.front(),
      &IRN_work.front(), &JCN_work.front(), &KEEP.front(),
      this->CNTL, this->ICNTL, &IW.front(),
      this->INFO, this->RINFO
    );

    for (
      unsigned kkk{0};
      this->INFO[3] > A_work.size() && kkk < 10;
      ++kkk
    ) {
      integer newsz{ this->INFO[3] };
      // Increase internal workspace:
      A_work.resize( newsz );
      std::fill( A_work.begin(), A_work.end(), 0 );
      std::copy( A.begin(), A.end(), A_work.begin() );

      IRN_work.resize( newsz );
      std::fill( IRN_work.begin(), IRN_work.end(), 0 );
      std::copy( IRN.begin(), IRN.end(), IRN_work.begin() );

      JCN_work.resize( newsz );
      std::fill( JCN_work.begin(), JCN_work.end(), 0 );
      std::copy( JCN.begin(), JCN.end(), JCN_work.begin() );

      HSL::ma48a<T>(
        this->nRow, this->nCol, this->nnz, 1,
        integer(A_work.size()), &A_work.front(),
        &IRN_work.front(), &JCN_work.front(), &KEEP.front(),
        this->CNTL, this->ICNTL, &IW.front(),
        this->INFO, this->RINFO
      );
    }

    if ( this->verbose ) msg_error();

    if ( this->INFO[0] < 0 ) return this->INFO[0];

    // factorization phase
    // -------------------
    HSL::ma48b<real_type>(
      this->nRow, this->nCol, this->nnz, 1,
      integer(A_work.size()), &A_work.front(),
      &IRN_work.front(), &JCN_work.front(), &KEEP.front(),
      this->CNTL, this->ICNTL,
      &W.front(), &IW.front(),
      this->INFO, this->RINFO
    );

    if ( this->verbose ) {
      msg_info();
      msg_error();
      fmt::print("done\n");
    }
    return this->INFO[0];
  }

  template <typename T>
  void
  MA48<T>::init(void) {
    IRN.clear();
    JCN.clear();
    A.clear();
  }

  template <typename T>
  void
  MA48<T>::insert(
    integer   i,
    integer   j,
    real_type a
  ) {
    UTILS_ASSERT(
      i >= 0 && i < this->nRow && j >= 0 && j < this->nCol,
      "Sparsetool: MA48::insert( {}, {}, a ) out of matrix\n", i, j
    );
    IRN.push_back(i+1);
    JCN.push_back(j+1);
    A.push_back(a);
  }

  template <typename T>
  int
  MA48<T>::setup( bool v ) {
    this->verbose = v;
    if ( this->verbose ) fmt::print("ma48::setup()\n");
    // set pars
    HSL::ma48i<real_type>(this->CNTL, this->ICNTL);
    //this->ICNTL[6] = 0; // handle structurally rank deficient matrices.
    this->KEEP.resize(this->nRow+9*this->nCol+7);
    this->IW.resize(6*this->nRow+3*this->nCol);
    this->W.resize(this->nRow);
    if ( this->verbose ) fmt::print("ma48::setup() done\n");
    return factorize();
  }

  template <typename T>
  void
  MA48<T>::solve( real_type X[], real_type const RHS[], bool transposed ) {
    real_type ERROR[3];
    integer   NM{std::max(this->nRow,this->nCol)};
    W.resize(2*NM);
    HSL::ma48c<real_type>(
      this->nRow, this->nCol, (transposed ? 1 : 0), 1,
      integer(A_work.size()), &A_work.front(),
      &IRN_work.front(), &KEEP.front(),
      this->CNTL, this->ICNTL,
      RHS, X, ERROR,
      &W.front(), &IW.front(), this->INFO
    );

    if ( this->verbose ) {
      msg_infor();
      msg_error();
    }
  }

  template class MA48<float>;
  template class MA48<double>;

}
