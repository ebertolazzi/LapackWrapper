/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  SparseTool : DRIVER FOR CHECKING THE PACKAGE WITH AN ELLIPTIC PROBLEM   |
 |                                                                          |
 |  date         : 2008, 10 April                                           |
 |  version      : 1.0                                                      |
 |  file         : test.paper.1.cc                                          |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Department of Mechanics and Structures Engineering       |
 |                 University of Trento                                     |
 |                 via Mesiano 77, I -- 38050 Trento, Italy                 |
 |                 email : enrico.bertolazzi@ing.unitn.it                   |
 |                                                                          |
 |  purpose:                                                                |
 |    Build and solve the linear system associated to the finite diference  |
 |    approximation of the Laplacian in [0,1] x [0,1].                      |
 |                                                                          |
 |    u  (x,y)   + u  (x,y) = 1  in [0,1] x [0,1]                           |
 |     xx           yy                                                      |
 |                                                                          |
 |    u(x,y) = 0  on the border of [0,1] x [0,1]                            |
 |                                                                          |
 |    the finite difference will be:                                        |
 |                                                                          |
 |    -4*u    + u      + u     + u      + u       = h^2   i,j=1,2,...,n-1   |
 |        i,j    i-1,j    i+1,j   i,j-1    i,j+1                            |
 |                                                                          |
 |     u    = u     = u    = u    = 0                                       |
 |      0,j     n,j    i,0    i,n                                           |
 |                                                                          |
 |    where: h = 1/n                                                        |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#define SPARSETOOL_DEBUG
#include <sparse_tool/sparse_tool.hh>
#include <sparse_tool/sparse_tool_extra.hh>
#include <sparse_tool/sparse_tool_iterative.hh>
#include <sparse_tool/sparse_tool_matrix_market.hh>

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wshorten-64-to-32"
#endif

static Utils::Console msg(&std::cout);

using namespace ::Sparse_tool_load;
using           ::Sparse_tool::integer;
using namespace ::std;

inline
integer
ipos( integer const i, integer const j, integer const nx, integer /* ny */ ) {
  return (i-1)+(j-1)*(nx-1);
}

int
main() {
  try {
    integer           iter;
    constexpr integer maxiter{ 1000};
    constexpr integer n{100};
    constexpr integer neq{(n-1)*(n-1)};
    CCoorMatrix<double>        A(neq,neq,5*neq); // initialize matrix A with reseved rom for 5*neq nonzeros
    Vector<double>             b(neq), x(neq), r(neq);   // initialize r.h.s. and solution vector
    //ILDUpreconditioner<double> P; // P will be the incomplete LDU preconditioner
    //IdPreconditioner<double>   P; // P will be the incomplete LDU preconditioner
    //Dpreconditioner<double>   P; // P will be the incomplete LDU preconditioner
    //SORpreconditioner<double>   P; // P will be the incomplete LDU preconditioner
    SSORpreconditioner<double>   P; // P will be the incomplete LDU preconditioner
    b = -1.0/(n*n); // h^2
    for ( integer i{1}; i<n; ++i ) {
      for ( integer j{1}; j<n; ++j ) {
        integer const ii{ ipos( i, j, n, n ) };
        A.insert(ii,ii) = 4;
        if ( i > 1   ) A.insert(ii,ipos( i-1, j, n, n )) = -1;
        if ( i < n-1 ) A.insert(ii,ipos( i+1, j, n, n )) = -1;
        if ( j > 1   ) A.insert(ii,ipos( i, j-1, n, n )) = -1;
        if ( j < n-1 ) A.insert(ii,ipos( i, j+1, n, n )) = -1;
      }
    }
    A.internal_order(); // order internal structure of the sparse matrice
    //P.build(A);       // build incomplete LDU preconditioner
    P.build(A,1.6);     // build incomplete LDU preconditioner
    //double res = bicgstab( A, b, x, P, 1E-14, maxiter, iter, &cout );
    //double res = gmres( A, b, x, P, 1E-14, 50, maxiter, iter, &cout );
    //double res = cg( A, b, x, P, 1E-14, maxiter, iter, &cout );
    //double res = cg_poly( A, b, x, 1E-14, maxiter, 10, iter, &cout );
    double res = cocg( A, b, x, P, 1E-14, maxiter, iter, &cout );
    //double res = cocr( A, b, x, P, 1E-14, maxiter, iter, &cout );
    fmt::print( "res = {}\n", res );

    r = b - A*x;
    fmt::print(
      "Verify Residual = {}\n",
      r.lpNorm<Eigen::Infinity>()
    );

  } catch ( exception const & exc ) {
    msg.error( exc.what() );
  } catch ( ... ) {
    msg.error("Errore Sconosciuto!\n");
  }
  return 0;
}
