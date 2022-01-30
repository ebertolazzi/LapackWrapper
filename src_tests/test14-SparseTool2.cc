/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  SparseTool: DRIVER FOR CHECKING THE PACKAGE WITH A LEAST SQUARE PROBLEM |
 |                                                                          |
 |  date         : 2008, 10 April                                           |
 |  version      : 1.0                                                      |
 |  file         : test.paper.2.cc                                          |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Department of Mechanics and Structures Engineering       |
 |                 University of Trento                                     |
 |                 via Mesiano 77, I -- 38050 Trento, Italy                 |
 |                 email : enrico.bertolazzi@ing.unitn.it                   |
 |                                                                          |
 |  purpose:                                                                |
 |    Solve the overdetermined linear system Ax = b with the method of      |
 |    least squares.                                                        |
 |    The problem is solved using the augmented block linear system:        |
 |                                                                          |
 |       /         \  /   \     /   \                                       |
 |       | I     A |  | x |     | b |                                       |
 |       |         |  |   |   = |   |                                       |
 |       | A^T   0 |  | r |     | 0 |                                       |
 |       \         /  \   /     \   /                                       |
 |                                                                          |
 |   where the coefficient matrix A and r.h.s. vector b are:                |
 |                                                                          |
 |       /         \         /   \                                          |
 |       | 1       |         | 1 |          /    \                          |
 |       |   1     |         | 1 |          | x0 |                          |
 |       |     1   |         | 1 |          | x1 |                          |
 |   A = |       1 |     b = | 1 |      x = | x2 |                          |
 |       | 1 1 1 1 |         | 1 |          | x3 |                          |
 |       | 1 2     |         | 1 |          \    /                          |
 |       |     1 2 |         | 1 |                                          |
 |       \         /         \   /                                          |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#define SPARSETOOL_DEBUG
#include <sparse_tool/sparse_tool.hh>
#include <sparse_tool/sparse_tool_extra.hh>
#include <sparse_tool/sparse_tool_iterative.hh>
#include <sparse_tool/sparse_tool_matrix_market.hh>

static Utils::Console msg(&std::cout);

using namespace ::Sparse_tool_load;
using           ::Sparse_tool::integer;
using namespace ::std;

typedef typename Vector<double>::V_base VBASE;

int
main() {
  integer iter, nr = 7, nc = 4, nnz = 12;
  integer I[] = { 0, 4, 5, 1, 4, 5, 2, 4, 6, 3, 4, 6 };
  integer J[] = { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3 };
  integer V[] = { 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2 };

  try {
    CCoorMatrix<double> A(nr,nc),   M(nr+nc,nr+nc);
    Vector<double>      rhs(nr+nc), X(nr+nc), bf1(nr), bf2(nc);
  
    // build matrix A
    for ( integer i = 0; i < nnz; ++i ) A.insert(I[i],J[i]) = V[i];
    A.internalOrder();
  
    // build block matrix
    //   /           \
    //   |  I     A  |
    //   |           |
    //   |  A^T   0  |
    //   \           /
    for ( A.Begin(); A.End(); A.Next() ) {
      integer i = A.row(), j = A.column();
      M.insert(i,j+nr) = A.value();
      M.insert(j+nr,i) = A.value();
    }
    for ( integer i = 0; i < nr; ++i ) M.insert(i,i) = 1;
    M.internalOrder();
  
    Eigen::Map<VBASE> b1(NULL,0), b2(NULL,0), r(NULL,0), x(NULL,0);
    new (&b1) Eigen::Map<VBASE>(rhs.data(), nr);    b1.array() = 1;
    new (&b2) Eigen::Map<VBASE>(rhs.data()+nr, nc); b2.array() = 0;
    new (&r)  Eigen::Map<VBASE>(X.data(), nr);      r.array()  = 0;
    new (&x)  Eigen::Map<VBASE>(X.data()+nr, nc);   x.array()  = 0;

    //b1.slice( rhs, 0, nr ) = 1; b2.slice( rhs, nr, nc+nc ) = 0;
    //r.slice(  X,   0, nr ) = 0; x.slice(  X,   nr, nr+nc ) = 0;

    ILDUpreconditioner<double> P(M);
    //IdPreconditioner<double>   P(M);
    //Dpreconditioner<double>    P(M);
    //SORpreconditioner<double>  P(M,1.2);
    //SSORpreconditioner<double> P(M,1.2);

    double res = bicgstab( M, rhs, X, P, 1E-15, 100, iter, &cout ); // solve extended linear system
  
    bf1 = b1 - (A*x); bf2 = A^bf1; // Chech the projected residual A^T(b-A*x)
    fmt::print(
      "res      = {}\n"
      "solution = {}\n"
      "residual = {}\n",
      res, x, bf2.template lpNorm<Eigen::Infinity>()
    );
  } catch ( exception const & exc ) {
    msg.error( exc.what() );
  } catch ( ... ) {
    msg.error("Errore Sconosciuto!\n");
  }
  return 0;
}
