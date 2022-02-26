/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  SparseTool   : DRIVER FOR CHECKING THE SPY COMMAND                      |
 |                                                                          |
 |  date         : 2008, 12 April                                           |
 |  version      : 1.0                                                      |
 |  file         : test.paper.spy.cc                                        |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Department of Mechanics and Structures Engineering       |
 |                 University of Trento                                     |
 |                 via Mesiano 77, I -- 38050 Trento, Italy                 |
 |                 email : enrico.bertolazzi@ing.unitn.it                   |
 |                                                                          |
 |  purpose:                                                                |
 |                                                                          |
 |    Generate a sparse pattern and save a silhouette of the pattern        |
 |    to the file : S.eps.                                                  |
 |    The position for an horizontal and a vertical row are contained in    |
 |    the vectors rr and rc.                                                |
 |                                                                          |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "SparseTool.hh"

using namespace SparseToolLoad ;
using namespace std ;

int main() {
  SparsePattern sp(10,20) ;
  sp.insert(0,0).insert(1,1).insert(2,2).insert(3,3).insert(3,15);
  sp.insert(4,3).insert(5,10).insert(5,11).insert(5,12).insert(6,14);
  sp.insert(7,13).insert(8,6).insert(8,14).insert(9,14).insert(9,19);
  sp.internal_order() ;
  vector<integer> rc, rr ; rc.push_back(3); rr.push_back(5);
  Spy( "S.eps", sp, 12.0, &rr, &rc );
  return 0 ;
}
