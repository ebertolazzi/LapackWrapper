/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  SparseTool   : DRIVER FOR TESTING THE TOOLKIT INTERFACING WITH SuperLU  |
 |                                                                          |
 |  date         : 2008, 14 April                                           |
 |  version      : 1.0                                                      |
 |  file         : test.SLU.cc                                              |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Department of Mechanics and Structures Engineering       |
 |                 University of Trento                                     |
 |                 via Mesiano 77, I -- 38050 Trento, Italy                 |
 |                 email : enrico.bertolazzi@ing.unitn.it                   |
 |                                                                          |
 |  purpose:                                                                |
 |                                                                          |
 |    Test the interface with SuperLU reading a matrix from a MatrixMarket  |
 |    file and solving a linear system.                                     |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#define SPARSETOOL_DEBUG
#include "SparseTool.hh"
#include "SparseToolSuperLU.hh"
#include "SparseToolExtra.hh"
#include <fstream>

using namespace SparseToolLoad ;
using namespace std ;

void
testSLU( string const & mm_file ) {
  Timing              tm ;
  MatrixMarket        mm ;
  superLU             slu ;
  CCoorMatrix<double> A ;
  Vector<double>      x, rhs, exact, resid ;

  mm . read( mm_file ) ;
  cout << mm ;
  
  cout << "load matrix..." << flush ;
  mm . load( A ) ;
  cout << "done\n" ;

  Spy( mm_file + ".eps" , A, 15.0 ) ;

  exact . resize( A . numRows() ) ;
  x     . resize( A . numRows() ) ;
  rhs   . resize( A . numRows() ) ;
  resid . resize( A . numRows() ) ;

  cout << "factorize (SLU) ..." << flush ;
  tm . start() ;
  cout << slu . load(A) ;
  tm . stop() ;
  cout << " " << tm . milliseconds() << "[ms] done\n" ;

  exact = 1 ;
  rhs = A * exact ;

  cout << "solve (SLU) ... " << flush ;
  tm . start() ;
  cout << slu . solve(rhs,x) << " " ;
  tm . stop() ;
  cout << tm . milliseconds()  << "[ms] done\n" ;
  
  resid = rhs - A*x ;

  cout << "\nerror    (SLU) = " << dist2( x, exact )
       << "\nresidual (SLU) = " << normi( resid )
       << "\n" ;
  
}

int
main() {
  char *rMatrix[] = { "hor__131.mtx",
                      "af23560.mtx",
                      "plat1919.mtx",
                      NULL } ;
  for ( char **p = rMatrix ; *p != NULL ; ++p ) testSLU( string("mm/")+*p) ;
  return 0 ;
}
