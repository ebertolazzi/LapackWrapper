/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  SparseTool   : DRIVER FOR TESTING ALL FEATURE OF THE TOOLKIT            |
 |                                                                          |
 |  date         : 2008, 1 April                                            |
 |  version      : 1.0                                                      |
 |  file         : SparseTool_timing.cc                                     |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Department of Mechanics and Structures Engineering       |
 |                 University of Trento                                     |
 |                 via Mesiano 77, I -- 38050 Trento, Italy                 |
 |                 email : enrico.bertolazzi@ing.unitn.it                   |
 |                                                                          |
 |  purpose:                                                                |
 |                                                                          |
 |    Test the toolkit for near all the operators and functions             |
 |                                                                          |
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

using namespace ::Sparse_tool_load;
using           ::Sparse_tool::integer;
using namespace ::std;

typedef double Real;

#include <Utils.hh>

inline
Real power2( Real const &  a)
{ return a*a; }

#include <fstream>

#define LOOP(N)   for ( integer i = 0; i < N; ++i )
#define REPEAT(N) for ( integer iii = 0; iii < N; ++iii )

typedef double Real;

static
void
out(
  integer N,
  integer nnz,
  Real    ta,
  Real    tb
) {
  fmt::print(
    "N={:<3} TIME ({:.5}:{:.5}) FP({:.5}:{:.5}) (sparselib:hand made)\n",
    N, ta, tb, 1e-3*N*nnz/ta, 1e-3*N*nnz/tb
  );
}

// test CRow  ********************************************************

static
void 
test_CRow(
  integer N,
  CRowMatrix<Real> const & crow,
  Vector<Real>     const & v
) {

  Utils::TicToc tm;
  Vector<Real> res( v.size() );

  tm.tic();
  REPEAT(N) res = crow * v;
  tm.toc();
  Real timea = tm.elapsed_ms();

  tm.tic();
  REPEAT(N) {
    integer const * R = crow.getR().data();
    integer const * J = crow.getJ().data();
    Real    const * A = crow.getA().data();

    for ( integer i = 0; i < crow.numRows(); ++i, ++R ) {
      integer n = R[1]-R[0];
      Real tmp = 0;
      while ( n-- > 0 ) tmp += *A++ * v(*J++);
      res(i) = tmp;
    }
  }
  tm.toc();
  Real timeb = tm.elapsed_ms();

  cout << "CRowMatrix:  ";
  out(N, crow.nnz(), timea, timeb);

}

// test CCol ********************************************************

static
void
test_CCol(
  integer                  N,
  CColMatrix<Real> const & ccol,
  Vector<Real>     const & v
) {

  Utils::TicToc tm;
  Vector<Real> res( v.size() );

  tm.tic();
  REPEAT(N) res = ccol * v;
  tm.toc();
  Real timea = tm.elapsed_ms();

  tm.tic();
  REPEAT(N) {
    integer const * I = ccol.getI().data();
    integer const * C = ccol.getC().data();
    Real    const * A = ccol.getA().data();

    res = 0;
    for ( integer i = 0; i < ccol.numCols(); ++i, ++C ) {
      Real tmp = v(i);
      integer n = C[1]-C[0];
      while ( n-- > 0 ) res(*I++) += *A++ * tmp;
    }
  }
  tm.toc();
  Real timeb = tm.elapsed_ms();

  cout << "CColMatrix:  ";
  out( N, ccol.nnz(), timea, timeb);
}

// test CCoor ********************************************************
static
void
test_CCoor(
  integer                   N,
  CCoorMatrix<Real> const & ccoor,
  Vector<Real>      const & v
) {

  Utils::TicToc tm;
  Vector<Real> res( v.size() );
  
  tm.tic();
  REPEAT(N) { res = ccoor * v; }
  tm.toc();
  Real timea = tm.elapsed_ms();
  
  tm.tic();
  REPEAT(N) {
    integer const * I = ccoor.getI().data();
    integer const * J = ccoor.getJ().data();
    Real    const * A = ccoor.getA().data();
    res = 0;
    for ( integer i = 0; i < ccoor.nnz(); ++i )
      res(*I++) += *A++ * v(*J++);
  }
  tm.toc();
  Real timeb = tm.elapsed_ms();

  cout << "CCoorMatrix: ";
  out( N, ccoor.nnz(), timea, timeb );
}

static
void
test_timing() {
  Utils::TicToc tm;

  integer const n = 10000;

  Real scal = 1.23;
  Vector<Real> a(n), b(n), res(n), res1(n);

  LOOP(n) {
    a[i] = ((i+1)%7)  * 2.134 + 1;
    b[i] = ((i+2)%31) * 3.234 + 1;
  }

  res = a + b;
  LOOP(n) res1[i] = a[i] + b[i];
  cout << "test vector + vector " << (res-res1).norm() << '\n';

  res = a - b;
  LOOP(n) res1[i] = a[i] - b[i];
  cout << "test vector - vector " << (res-res1).norm() << '\n';

  res = a.array() * b.array();
  LOOP(n) res1[i] = a[i] * b[i];
  cout << "test vector * vector " << (res-res1).norm() << '\n';

  res = a.array() / b.array();
  LOOP(n) res1[i] = a[i] / b[i];
  cout << "test vector / vector " << (res-res1).norm() << '\n';

  res = scal + a.array();
  LOOP(n) res1[i] = scal + a[i];
  cout << "test scal + vector " << (res-res1).norm() << '\n';

  res = scal - a.array();
  LOOP(n) res1[i] = scal - a[i];
  cout << "test scal - vector " << (res-res1).norm() << '\n';
  
  res = scal * a;
  LOOP(n) res1[i] = scal * a[i];
  cout << "test scal * vector " << (res-res1).norm() << '\n';

  res = a + scal * b;
  LOOP(n) res1[i] = a[i] + scal * b[i];
  cout << "test vector + scal * vector " << (res-res1).norm() << '\n';

  res = a - scal * b;
  LOOP(n) res1[i] = a[i] - scal * b[i];
  cout << "test vector + scal * vector " << (res-res1).norm() << '\n';
 
  cout << "\nMATRIX TESTS\n";
       
  char const *rMatrix[] = { "hor__131.mtx.gz", "af23560.mtx.gz", "plat1919.mtx.gz", NULL };
  char const **p        = rMatrix;
  
  Vector<Real> v;
  
  for (; *p != NULL; ++p ) {
  
    cout << "TEST matrix = " << *p << '\n';
    CCoorMatrix<Real> ccoor;
    CRowMatrix<Real>  crow;
    CColMatrix<Real>  ccol;

    MatrixMarket mm;
    cout << "load matrix..." << flush;
    mm.read( string("mm/") + *p, ccoor);
    cout << "done\n";

    mm.info(cout);

    Spy( string("mm/") + *p + ".eps", ccoor, 15.0 );

    tm.tic();
    crow = ccoor;
    tm.toc();
    Real timea = tm.elapsed_ms();

    tm.tic();
    ccol = ccoor;
    tm.toc();
    Real timeb = tm.elapsed_ms();

    fmt::print(
      "TO CROW {:.5} [ms] \n"
      "TO CCOL {:.5} [ms] \n",
      timea, timeb
    );

    integer const nr = ccoor.numRows();

    v.resize(nr);
    res.resize(nr);

    v = 1;
    v[nr/2] = 234;
    v[nr/4] = -3;
    v[(nr*3)/4 ] = -12;

    res = ccoor * v;
    res = res - crow * v;
    fmt::print( "CCoor - CRow = {:.5}\n", res.template lpNorm<Eigen::Infinity>() );

    res = ccoor * v;
    res = res - ccol * v;
    fmt::print( "CCoor - CCol = {:.5}\n", res.template lpNorm<Eigen::Infinity>() );
  
    integer cicle_repeat = 1 + 400000 / nr;
    test_CRow (cicle_repeat, crow,  v);
    test_CCol (cicle_repeat, ccol,  v);
    test_CCoor(cicle_repeat, ccoor, v);
  }

  cout << "\nTEST_TIME ALL DONE\n\n";
}


int
main() {
  test_timing();
  cout << "\nALL DONE FOLKS!\n\n";
  return 0;
}
