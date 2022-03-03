/*
 * \file ma48_wrapper.cc
 * Implementation of the class MA48.
 *
 * \author A. Huber and E.Bertolazzi
 * \since 15.11.2018
 */

#include "ma41_wrapper.hh"

#ifndef F77NAME
  #define F77NAME(A) A##_
#endif

namespace lapack_wrapper {

  typedef bool logical;

  template <typename real>
  void
  MA41<real>::msg_info( ostream_type & s ) const {
    fmt::print( s,
      "Real    space needs for factors    = {}\n"
      "Integer space needs for factors    = {}\n"
      "Maximum frontal size               = {}\n"
      "Nodes in the tree                  = {}\n"
      "Minimum MAXIS                      = {}\n"
      "Minimum MAXS                       = {}\n"
      "Real    Space for LU factorizarion = {}\n"
      "Integer Space for LU factorizarion = {}\n"
      "Integer Space for LU factorizarion = {}\n"
      "Largest frontal matrix             = {}\n",
      INFO[2], INFO[3], INFO[4], INFO[5], INFO[6],
      INFO[7], INFO[8], INFO[9], INFO[10], INFO[11]
    );
  }

  template <typename real>
  void
  MA41<real>::msg_infor( ostream_type & s ) const {
    fmt::print( s,
      "NFLOPs for LU factorization    = {}\n"
      "NFLOPs for assembly process    = {}\n"
      "NFLOPs for elimination process = {}\n",
      RINFO[0], RINFO[1], RINFO[2]
    );
    if ( ICNTL[10] > 0 ) {
      fmt::print( s,
        "Infinity norm of the matrix    = {}\n"
        "Infinity norm of the solution  = {}\n"
        "Norm of scaled residual        = {}\n",
        RINFO[3], RINFO[4], RINFO[5]
      );
    }
  }

  template <typename real>
  void
  MA41<real>::msg_error( ostream_type & s ) const {
    if ( INFO[0] >= 0 ) return;
    fmt::print( s,
      "***** ERROR *****\n"
      "dim = {}\n"
      "nnz = {}\n",
      N, NE
    );
    switch ( INFO[0] ) {
    case -1:
      fmt::print( s, "Value of N is out of range N = {}\n", INFO[1] );
      break ;
    case -2:
      fmt::print( s, "Value of NE is out of range NE = {}\n", INFO[1] );
      break ;
    case -3:
      fmt::print( s, "JOB has wrong value or analisys was not "
                     "performed prior to factorization.\n"
                     "JOB = {}\n", INFO[1]);
      break ;
    case -4:
      fmt::print( s, "Error in permutation error\n" );
      break ;
    case -5:
      fmt::print( s,
        "Not enought space to preprocess the input matrix\n"
        "MAXS should be increased to at least {}\n",
        INFO[1]
      );
      break ;
    case -6:
      fmt::print( s,
        "The matrix is structurally singular\n"
        "Estimated rank = {}\n",
        INFO[1]
      );
      break ;
    case -7:
      fmt::print( s,
        "Error from analysis\n"
        "MAXIS should be increased to at least {}\n",
        INFO[1]
      );
      break ;
    case -8:
      fmt::print( s,
        "Error from numerical factorization\n"
        "MAXIS should be increased to at least {}\n",
        INFO[1]
      );
      break ;
    case -9:
      fmt::print( s,
        "Error from numerical factorization\n"
        "MAXS should be increased to at least {}\n",
        INFO[1]
      );
      break ;
    case -10:
      fmt::print( s,
        "The matrix is numerically singular\n"
        "Estimated rank = {}\n",
        INFO[1]
      );
      break ;
    case -11:
      fmt::print( s,
        "Error from the solution phase\n"
        "MAXS should be increased to at least {}\n",
        INFO[1]
      );
      break ;
    case -12:
      fmt::print( s,
        "Not enought space to postprocess the solution\n"
        "MAXS should be increased to at least {}\n",
        INFO[1]
      );
      break ;
    }
    fmt::print( s, "***** ERROR *****\n");
  }

  // private functionalities
  template <typename real>
  void
  MA41<real>::load_matrix( integer nr ) {
    if ( verbose ) fmt::print( "MA41::load_matrix(...)\n");
    N  = nr ;
    NE = A . size() ;

    COLSCA . resize(N) ;
    ROWSCA . resize(N) ;
    IS     . resize(2 * NE + 11 * N + 1) ;
    S      . resize(2 * NE + 11 * N + 1) ;

    fill(COLSCA . begin(), COLSCA . end(), 0 ) ;
    fill(ROWSCA . begin(), ROWSCA . end(), 0 ) ;
    fill(IS     . begin(), IS     . end(), 0 ) ;
    fill(S      . begin(), S      . end(), 0 ) ;
     // set pars
    lapack_wrapper::ma41i<real>(CNTL, ICNTL, KEEP) ;
    if ( verbose ) fmt::print( "done\n" );
  }

  template <typename real>
  void
  MA41<real>::symbfac() {
    if ( verbose ) fmt::print( "MA41::symbfac() perform analisys\n" );

    // analysis phase
    vector_integer tmpIS(2 * NE + 12 * N + 1);
    integer JOB = ANALISYS;
    lapack_wrapper::ma41a<real>(
      JOB, N, NE,
      IRN.data(),
      JCN.data(),
      A.data(),
      NULL,
      COLSCA.data(),
      ROWSCA.data(),
      KEEP,
      tmpIS.data(), tmpIS.size(),
      S.data(),     S.size(),
      CNTL, ICNTL, INFO, RINFO
    );

    if ( verbose ) msg_error(std::cerr);

    IS . resize(10*INFO[6]);
    S  . resize(10*INFO[7]);

    copy( tmpIS.begin(), tmpIS.end(), IS.begin() ) ;

    if ( verbose ) fmt::print( "MA41::symbfac() do factorization\n" );

    // factorization phase
    // -------------------
    JOB = FACTORIZATION;
    lapack_wrapper::ma41a<real>(
      JOB, N, NE,
      IRN.data(),
      JCN.data(),
      A.data(),
      NULL,
      COLSCA.data(),
      ROWSCA.data(),
      KEEP,
      IS.data(), IS.size(),
      S.data(),  S.size(),
      CNTL, ICNTL, INFO, RINFO
    );

    if ( verbose ) {
      msg_info( std::cout );
      msg_error(  std::cout );
      std::cout << "done\n";
    }
  }

  template <typename real>
  void
  MA41<real>::solve( real RHS[] ) {
    integer JOB = SOLVE;
    lapack_wrapper::ma41a<real>(
      JOB, N, NE,
      IRN.data(),
      JCN.data(),
      A.data(),
      RHS,
      COLSCA.data(),
      ROWSCA.data(),
      KEEP,
      IS.data(), IS.size(),
      S.data(),  S.size(),
      CNTL, ICNTL, INFO, RINFO
    );

    if ( verbose ) {
      msg_infor( std::cout ) ;
      msg_error( std::cout ) ;
    }
  };

  template class MA41<float>;
  template class MA41<double>;
}
