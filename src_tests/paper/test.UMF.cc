/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  SparseTool   : DRIVER FOR TESTING THE TOOLKIT INTERFACING WITH UMFPACK  |
 |                                                                          |
 |  date         : 2008, 14 April                                           |
 |  version      : 1.0                                                      |
 |  file         : test.UMF.cc                                              |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Department of Mechanics and Structures Engineering       |
 |                 University of Trento                                     |
 |                 via Mesiano 77, I -- 38050 Trento, Italy                 |
 |                 email : enrico.bertolazzi@ing.unitn.it                   |
 |                                                                          |
 |  purpose:                                                                |
 |                                                                          |
 |    Test the interface with UMFPACK reading a matrix from a MatrixMarket  |
 |    file and solving a linear system.                                     |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#define SPARSETOOL_DEBUG
#include "SparseTool.hh"
#include "SparseToolUmf.hh"
#include "SparseToolExtra.hh"
#include <fstream>

using namespace SparseToolLoad ;
using namespace std ;

void
testUMF( string const & mm_file ) {
  Timing              tm ;
  MatrixMarket        mm ;
  UMF                 umf ;
  CCoorMatrix<double> A ;
  Vector<double>      x, rhs, exact, resid ;

  mm . read( mm_file ) ;
  cout << mm ;
  
  cout << "load matrix..." << flush ;
  mm . load( A ) ;
  cout << "done\n" ;

  //Spy( mm_file + ".eps" , A, 15.0 ) ;

  exact . resize( A . numRows() ) ;
  x     . resize( A . numRows() ) ;
  rhs   . resize( A . numRows() ) ;
  resid . resize( A . numRows() ) ;

  cout << "factorize (UMF) ..." << flush ;
  tm . start() ;
  cout << umf . load(A) ;
  tm . stop() ;
  cout << " " << tm . milliseconds() << "[ms] done\n" ;

  exact = 1 ;
  rhs   = A * exact ;

  cout << "solve (UMF) ... " << flush ;
  tm . start() ;
  cout << umf . solve(rhs,x) ;
  tm . stop() ;
  cout << " " << tm . milliseconds()  << "[ms] done\n" ;
  
  resid = rhs - A*x ;

  cout << "\nerror    (UFM) = " << dist2( x, exact )
       << "\nresidual (UFM) = " << normi( resid )
       << "\n" ;
  
}


int
main() {
  char *rMatrix[] = { "hor__131.mtx",
                      "af23560.mtx",
                      "plat1919.mtx",
                      NULL } ;
  for ( char **p = rMatrix ; *p != NULL ; ++p ) testUMF( string("mm/")+*p) ;
  return 0 ;
}
