/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  SparseTool   : DRIVER FOR TESTING THE TOOLKIT INTERFACING WITH MA41     |
 |                                                                          |
 |  date         : 2008, 22 April                                           |
 |  version      : 1.0                                                      |
 |  file         : test.MA41.cc                                             |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Department of Mechanics and Structures Engineering       |
 |                 University of Trento                                     |
 |                 via Mesiano 77, I -- 38050 Trento, Italy                 |
 |                 email : enrico.bertolazzi@ing.unitn.it                   |
 |                                                                          |
 |  purpose:                                                                |
 |                                                                          |
 |    Test the interface with MA41 reading a matrix from a MatrixMarket     |
 |    file and solving a linear system.                                     |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#define SPARSETOOL_DEBUG
#include "SparseTool.hh"
#include "SparseToolMA41.hh"
#include "SparseToolExtra.hh"
#include <fstream>

using namespace SparseToolLoad;
using namespace std;

void
testMA41( string const & mm_file ) {
  Timing              tm;
  MatrixMarket        mm;
  MA41                ma41;
  CCoorMatrix<double> A;
  Vector<double>      x, rhs, exact, resid;

  mm.read( mm_file );
  cout << mm;

  cout << "load matrix..." << flush;
  mm . load( A );
  cout << "done\n";

  //Spy( mm_file + ".eps" , A, 15.0 );

  exact . resize( A.nrows() );
  x     . resize( A.nrows() );
  rhs   . resize( A.nrows() );
  resid . resize( A.nrows() );

  cout << "factorize (MA41) ..." << flush;
  tm . start();
  cout << ma41 . load(A);
  tm . stop();
  cout << " " << tm . milliseconds() << "[ms] done\n";

  exact = 1;
  rhs   = A * exact;

  cout << "solve (MA41) ... " << flush;
  tm . start();
  cout << ma41 . solve(rhs,x);
  tm . stop();
  cout << " " << tm . milliseconds()  << "[ms] done\n";

  resid = rhs - A*x;

  cout << "\nerror    (MA41) = " << dist2( x, exact )
       << "\nresidual (MA41) = " << normi( resid )
       << "\n";

}


int
main() {
  char *rMatrix[] = { "hor__131.mtx",
                      "af23560.mtx",
                      "plat1919.mtx",
                      NULL };
  for ( char **p = rMatrix; *p != NULL; ++p ) testMA41( string("mm/")+*p);
  return 0;
}
