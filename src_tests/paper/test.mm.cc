/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  SparseTool   : DRIVER FOR TESTING THE TOOLKIT INTERFACING WITH          |
 |                 Matrix Market                                            |
 |                                                                          |
 |  date         : 2008, 14 April                                           |
 |  version      : 1.0                                                      |
 |  file         : test.mm.cc                                               |
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

#include "SparseTool.hh"

using namespace SparseToolLoad ;
using namespace std ;

void
testMM( string const & mm_file ) {
  MatrixMarket                  mm ;
  SparsePattern                 sp ;
  CCoorMatrix<int>              Ai ;
  CCoorMatrix<double>           Ar ;
  CCoorMatrix<complex<double> > Ac ;

  cout << "Reading: " << mm_file << "..." << flush ;
  mm . read( mm_file ) ;
  cout << "done\n" << mm << flush ;
  
  cout << "\nLoading Matrix..." << flush ;
  switch ( mm . value_type() ) {
    case MM_PATTERN: mm . load(sp) ; break ;
    case MM_INTEGER: mm . load(Ai) ; break ;
    case MM_REAL:    mm . load(Ar) ; break ;
    case MM_COMPLEX: mm . load(Ac) ; break ;
  }

  cout << "done\nWriting Matrix..." << flush ;
  switch ( mm . value_type() ) {
    case MM_PATTERN:
      MatrixMarketSaveToFile( mm_file + ".txt", sp, mm . matrix_type() ) ;
    break ;
    case MM_INTEGER:
      MatrixMarketSaveToFile( mm_file + ".txt",
                              Ai, mm . value_type(), mm . matrix_type() ) ;
    break ;
    case MM_REAL:
      MatrixMarketSaveToFile( mm_file + ".txt",
                              Ar, mm . value_type(), mm . matrix_type() ) ;
    break ;
    case MM_COMPLEX:
      MatrixMarketSaveToFile( mm_file + ".txt",
                              Ac, mm . value_type(), mm . matrix_type() ) ;
    break ;
  }
  cout << "done\n" ;

}

int
main() {

  char *Matrix[] = { "hor__131.mtx",             // real matrices
                     "af23560.mtx",
                     "plat1919.mtx",
                     "dwg961b.mtx",              // complex matrices
                     "bcsstk30.mtx",             // pattern
                     NULL } ;

  for ( char **p = Matrix ; *p != NULL ; ++p ) testMM( string("mm/")+*p) ;  
  return 0 ;
}