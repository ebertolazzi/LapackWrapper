/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Sparse_tool  : DRIVER FOR TESTING THE TOOLKIT INTERFACING WITH UMFPACK  |
 |                                                                          |
 |  date         : 2008, 14 April                                           |
 |  version      : 1.0                                                      |
 |  file         : test.UMF.cc                                              |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Dipartimento di Ingegneria Industriale                   |
 |                 Università degli Studi di Trento                         |
 |                 email : enrico.bertolazzi@unitn.it                       |
 |                                                                          |
 |  purpose:                                                                |
 |                                                                          |
 |    Test the interface with UMFPACK reading a matrix from a MatrixMarket  |
 |    file and solving a linear system.                                     |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#define SPARSETOOL_DEBUG
#include <sparse_tool/sparse_tool.hh>
#include <sparse_tool/sparse_tool_extra.hh>
#include <sparse_tool/sparse_tool_iterative.hh>
#include <sparse_tool/sparse_tool_matrix_market.hh>

//#include <sparse_tool/interfaces/MA41.hh>
//#include <sparse_tool/interfaces/MA48.hh>
//#include <sparse_tool/interfaces/SuperLU.hh>
//#include <sparse_tool/interfaces/mkl_pardiso.hh>
//#include <sparse_tool/interfaces/UMF.hh>

#include <fstream>
#include <iostream>

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

#ifdef __clang__
#pragma clang diagnostic pop
#endif

using namespace Sparse_tool_load;
using namespace std;

using Utils::istream_type;
using Utils::ostream_type;

static
void
testSparseTool( istream_type & mm_file ) {
  Utils::TicToc              tm;
  MatrixMarket               mm;
  IdPreconditioner<double>   preco;
  //Dpreconditioner<double>    preco;
  //SORpreconditioner<double>  preco;
  //SSORpreconditioner<double>  preco;
  //ILDUpreconditioner<double> preco;
  //ILDUKpreconditioner<double> preco;
  //ILDUiterPreconditioner<double> preco;
  //RILDUpreconditioner<double> preco;
  //SuperLUpreconditioner<double> preco;
  //CSSORpreconditioner<double> preco;
  //SORpreconditioner<double> preco;
  //JACOBIpreconditioner<double> preco;
  //UMFpreconditioner
  //HSS_CGSSOR_Preconditioner
  //HSS_CHEBYSHEV_Preconditioner
  //HSS_ILDU_Preconditioner
  //HSS_OPOLY_Preconditioner
  //HSS_OPOLY_SSOR_Preconditioner
  //HSS_SSOR_Preconditioner
  CCoorMatrix<double>        A;
  CRowMatrix<double>         B;
  SparsePattern              SP;
  Vector<double>             x, rhs, exact, resid;

  fmt::print("read matrix...");
  mm.read( mm_file, A );
  fmt::print("done\n{}",mm);

  SP.resize( A );
  B.resize( SP );
  //Spy( mm_file + ".eps" , A, 15.0 );

  exact.resize( A.nrows() );
  x.resize( A.nrows() );
  rhs.resize( A.nrows() );
  resid.resize( A.nrows() );

  exact.fill(1);
  rhs = A * exact;

#if 1

  x.setZero();
  fmt::print("factorize (preco) ...");
  tm.tic();
  //preco.build(A,1.5,1);
  preco.build(A);
  //preco.build(A,1.2);
  tm.toc();
  fmt::print(" {} [ms] done\n", tm.elapsed_ms());
  fmt::print("solve  (bicgstab) ...\n");
  tm.tic();

  double               epsi       { 1e-15 };
  Sparse_tool::integer maxIter    { 1000 };
  Sparse_tool::integer maxSubIter { 50 };
  Sparse_tool::integer iter;
  //double  res = bicgstab( A, rhs, x, preco, epsi, maxIter, iter, &cout );
  double  res = gmres( A, rhs, x, preco, epsi, maxSubIter, maxIter, iter, &cout );

  tm.toc();
  fmt::print(" {} [ms] done\n",tm.elapsed_ms());
  fmt::print("res = {}\n",res);

  resid = rhs - A*x;

  fmt::print("error    (bicgstab) = {}\n", (x-exact).norm() );
  fmt::print("residual (bicgstab) = {}\n", resid.lpNorm<Eigen::Infinity>() );
#endif

#if 0
  x.setZero();
  MA41<double> ma41;
  ma41.load( A );
  fmt::print("\nsolve    (MA41) ... ");
  tm.tic();
  ma41.solve( rhs, x );
  tm.toc();
  fmt::print(" {} [ms] done\n", tm.elapsed_ms());

  resid = rhs - A*x;
  fmt::print("error    (MA41) = {}\n", (x-exact).norm() );
  fmt::print("residual (MA41) = {}\n", resid.lpNorm<Eigen::Infinity>() );
#endif

#if 0
  x.setZero();
  MA48<double> ma48;
  ma48.load( A );
  fmt::print("\nsolve    (MA48) ... ");
  tm.tic();
  ma48.solve( rhs, x );
  tm.toc();
  fmt::print(" {} [ms] done\n", tm.elapsed_ms() );

  resid = rhs - A*x;
  fmt::print("error    (MA48) = {}\n", (x-exact).norm() );
  fmt::print("residual (MA48) = {}\n", resid.lpNorm<Eigen::Infinity>() );
#endif

#if 0
  x.setZero();
  mkl_PardisoRealU pardiso;
  pardiso.load( A );
  //pardiso.check_matrix();
  pardiso.factorize();
  fmt::print("\nsolve    (pardiso) ... ");
  tm.tic();
  pardiso.solve( rhs, x );
  tm.toc();
  fmt::print(" {} [ms] done\n", tm.elapsed_ms());

  resid = rhs - A*x;
  fmt::print("error    (pardiso) = {}\n", (x-exact).norm() );
  fmt::print("residual (pardiso) = {}\n", resid.lpNorm<Eigen::Infinity>() );
#endif

#if 0
  x.setZero();
  SuperLU<double> superlu;
  superlu.load( A );
  fmt::print("\nsolve    (superlu) ... ");
  tm.tic();
  superlu.solve( rhs, x );
  tm.toc();
  fmt::print(" {} [ms] done\n", tm.elapsed_ms());

  resid = rhs - A*x;
  fmt::print("error    (superlu) = {}\n", (x-exact).norm() );
  fmt::print("residual (superlu) = {}\n", resid.lpNorm<Eigen::Infinity>() );
#endif

#if 0
  x.setZero();
  UMF<double> umf;
  umf.load( A );
  fmt::print("\nsolve    (UMF) ... ");
  tm.tic();
  umf.solve( rhs, x );
  tm.toc();
  fmt::print(" {} [ms] done\n", tm.elapsed_ms());

  resid = rhs - A*x;
  fmt::print("error    (UMF) = {}\n", (x-exact).norm() );
  fmt::print("residual (UMF) = {}\n", resid.lpNorm<Eigen::Infinity> );
#endif

  resid = rhs - A*exact;
  fmt::print( "\n\nresidual (exact) = {}\n\n", resid.lpNorm<Eigen::Infinity>() );
}

int
main() {
  char const * rMatrix[]{
    //"af23560.mtx", // MA48 fails
    //"memplus.mtx", //
    //"fidap005.mtx",
    "ASIC_100k.mtx",           // 99340 (ok) 200 iter
    //"ASIC_320ks.mtx",          // 321671 (ok) 200 iter
    //"ASIC_680k.mtx",           // 682862 (ok) 200 iter
    //"af23560.mtx",             // 23560 (no)
    //"audikw_1.mtx",            // 943695 (ok) 11 iter
    //"barrier2-2.mtx",          // 113076 (no)
    //"bmw3_2.mtx",              // 227362 (ok) 8 iter
    //"cage15.mtxz",              // In reading Matrix Market File, bad pattern index on line 38025901
    //"CO.mtx",                  // 221119 (ok) 10 iter
    //"dwg961b.mtx",             // 961 (ok) 14
    //"ecology2.mtx",            // 999999 (ok) 53 iter
    //"fidapm05.mtx",            // 42 (no)
    //"memchip.mtx",             // 2707524 (ok) 200 iter
    //"hor__131.mtx",            // 434 (ok) 200 iter
    //"ldoor.mtx",               // 952203 (ok) 11 iter
    //"para-9.mtx",              // 155924 (no)
    //"parabolic_fem.mtx",       // 525825 (ok) 54 iter
    //"plat1919.mtx",            // 1919 (ok) 200 iter
    //"s3dkq4m2.mtx",            // 90449 (ok) 10 iter
    nullptr
  };

  for ( char const **p = rMatrix; *p != nullptr; ++p ) {
    string fname = string("mm/")+*p;
    ifstream file( fname.data() );
    if ( !file.good() ) {
      cerr << "Cannot open file: " << fname << "\n";
      exit(0);
    }
    testSparseTool( file );
    file.close();
  }
  fmt::print("\nAll Done Folks\n\n");
  return 0;
}
