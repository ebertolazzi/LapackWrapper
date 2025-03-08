/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  SparseTool   : DRIVER FOR TESTING PERFORMANCE of SparseTool             |
 |                                                                          |
 |  date         : 2008, 14 April                                           |
 |  version      : 1.0                                                      |
 |  file         : test.performance.cc                                      |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Department of Mechanics and Structures Engineering       |
 |                 University of Trento                                     |
 |                 via Mesiano 77, I -- 38050 Trento, Italy                 |
 |                 email : enrico.bertolazzi@ing.unitn.it                   |
 |                                                                          |
 |  purpose:                                                                |
 |                                                                          |
 |    Compare the performace of SparseTool with SPARSKIT and Sparse BLAS    |
 |                                                                          |
\*--------------------------------------------------------------------------*/


//#define SPARSETOOL_DEBUG
#include "SparseTool.hh"
#include "SparseToolExtra.hh"
#include "blas_enum.h"
#include "blas_sparse_proto.h"

#define F77NAME(A) A##_

extern "C" {
  void F77NAME(amux)( int          * n,
                      double const * x,
                      double       * y,
                      double const * a,
                      int    const * ja,
                      int    const * ia );
}

using namespace ::SparseToolLoad;
using namespace ::std;

#define LOOP(N)   for ( integer i{0}; i < N; ++i )
#define REPEAT(N) for ( integer iii{0}; iii < N; ++iii )

using Real = double;

blas_sparse_matrix mHandle;

Vector<Real> A;
Vector<int>  R;
Vector<int>  J;

static
void
out(
  integer      N,
  integer      nnz,
  Real const & ta,
  Real const & tb,
  Real const & tc
) {
  cout << "N=" << setw(3) << N
       << " TIME ("
       << setw(10) << ta << ":"
       << setw(10) << tb << ":"
       << setw(10) << tc << " ) FP("
       << setw(10) << 1e-3*nnz/ta
       << ":"
       << setw(10) << 1e-3*nnz/tb
       << ":"
       << setw(10) << 1e-3*nnz/tc
       << ")  (sparselib:BLAS:SPARSKIT)\n";
}

// test CRow  ********************************************************

static
void
test_CRow(
  integer                  N,
  blas_sparse_matrix       mHandle,
  CRowMatrix<Real> const & crow,
  Vector<Real>     const & v
) {

  Timing tm;
  Vector<Real> res( v . size() );

  sleep(1);
  tm . start();
  REPEAT(N) res = crow * v;
  tm . stop();
  Real timea = tm . milliseconds();

  int n = v.size();
  double const * pv { &v(0) };
  double       * pr { &res(0) };
  double const * pA { &A(0) };
  int    const * pR { &R(0) };
  int    const * pJ { &J(0) };

  sleep(1);
  tm . start();
  REPEAT(N) BLAS_dusmv( blas_no_trans, 1.0, mHandle, pv, 1, pr, 1 );
  tm . stop();
  Real timeb = tm . milliseconds();

  sleep(1);
  tm . start();
  REPEAT(N) F77NAME(amux)( &n, pv, pr, pA, pJ, pR );
  tm . stop();
  Real timec = tm . milliseconds();

  cout << "CRowMatrix:  ";
  out(N, crow . nnz(), timea/N, timeb/N, timec/N );

}

// test CCol ********************************************************

static
void
test_CCol(
  integer                  N,
  blas_sparse_matrix       mHandle,
  CColMatrix<Real> const & ccol,
  Vector<Real>     const & v
) {

  Timing tm;
  Vector<Real> res( v.size() );

  sleep(1);
  tm.start();
  REPEAT(N) res = ccol * v;
  tm.stop();
  Real timea = tm.milliseconds();

  int n = v.size();
  double const * pv { &v(0) };
  double       * pr { &res(0) };
  double const * pA { &A(0) };
  int    const * pR { &R(0) };
  int    const * pJ { &J(0) };

  sleep(1);
  tm . start();
  REPEAT(N) BLAS_dusmv( blas_no_trans, 1.0, mHandle, pv, 1, pr, 1 );
  tm . stop();
  Real timeb = tm . milliseconds();

  sleep(1);
  tm . start();
  REPEAT(N) F77NAME(amux)( &n, pv, pr, pA, pJ, pR );
  tm . stop();
  Real timec = tm . milliseconds();

  cout << "CColMatrix:  ";
  out(N, ccol . nnz(), timea/N, timeb/N, timec/N );
}

// test CCoor ********************************************************
static
void
test_CCoor(
  integer N,
  blas_sparse_matrix        mHandle,
  CCoorMatrix<Real> const & ccoor,
  Vector<Real>      const & v
) {

  Timing tm;
  Vector<Real> res( v.size() );

  sleep(1);
  tm.start();
  REPEAT(N) { res = ccoor * v; }
  tm.stop();
  Real timea = tm.milliseconds();

  int n = v.size();
  double const * pv { &v(0)   };
  double       * pr { &res(0) };
  double const * pA { &A(0)   };
  int    const * pR { &R(0)   };
  int    const * pJ { &J(0)   };

  sleep(1);
  tm.start();
  REPEAT(N) BLAS_dusmv( blas_no_trans, 1.0, mHandle, &v(0), 1, &res(0), 1 );
  tm.stop();
  Real timeb = tm.milliseconds();

  sleep(1);
  tm.start();
  REPEAT(N) F77NAME(amux)( &n, pv, pr, pA, pJ, pR );
  tm.stop();
  Real timec = tm . milliseconds();

  cout << "CCoorMatrix: ";
  out( N, ccoor . nnz(), timea/N, timeb/N, timec/N );
}

int
main() {
  Timing tm;

  cout << "\nMATRIX TESTS\n";

  char const *rMatrix[] = { "hor__131.mtx", "af23560.mtx", "plat1919.mtx", NULL };
  char const **p = rMatrix;

  Vector<Real> v, res;

  for (; *p != NULL; ++p ) {

    cout << "TEST matrix = " << *p << '\n';
    CCoorMatrix<Real> ccoor;
    CRowMatrix<Real>  crow;
    CColMatrix<Real>  ccol;

    MatrixMarket mm;
    mm . read( string("mm/") + *p );
    cout << mm;

    cout << "load matrix..." << flush;
    mm . load( ccoor );
    cout << "done\n";

    Spy( string("mm/") + *p + ".eps", ccoor, 15.0 );

    tm.start();
    crow = ccoor;
    tm.stop();
    Real timea = tm.milliseconds();

    tm.start();
    ccol = ccoor;
    tm.stop();
    Real timeb = tm.milliseconds();

    tm.start();
    blas_sparse_matrix mHandle = BLAS_duscr_begin( ccoor.nrows(), ccoor.ncols() );
    BLAS_duscr_insert_entries(
      mHandle, ccoor.nnz(),
      ccoor.getA().data(),
      (int*)ccoor.getI().data(),
      (int*)ccoor.getJ().data()
    );
    BLAS_suscr_end( mHandle );
    tm.stop();
    Real timec = tm.milliseconds();

    cout << "TO CROW  " << timea << "[ms]" << '\n';
    cout << "TO CCOL  " << timeb << "[ms]" << '\n';
    cout << "TO SBLAS " << timec << "[ms]" << '\n';

    A.resize(ccoor.nnz());
    R.resize(ccoor.nrows()+1);
    J.resize(ccoor.nnz());
    std::copy( crow.getA().begin(), crow.getA().end(), A.begin() );
    std::copy( crow.getR().begin(), crow.getR().end(), R.begin() );
    std::copy( crow.getJ().begin(), crow.getJ().end(), J.begin() );
    for ( integer i{0}; i < R.size(); ++i ) ++R[i];
    for ( integer i{0}; i < J.size(); ++i ) ++J[i];

    integer const nr{ ccoor.nrows() };

    v.resize(n);
    res.resize(n);

    v = 1;
    v[nr/2] = 234;
    v[nr/4] = -3;
    v[(nr*3)/4 ] = -12;

    res = ccoor * v;
    res -= crow * v;
    cout << "CCoor - CRow = " << normi(res) << '\n';

    res = ccoor * v;
    res -= ccol * v;
    cout << "CCoor - CCol = " << normi(res) << '\n';

    res = ccoor * v;
    BLAS_dusmv( blas_no_trans, -1.0, mHandle, &v(0), 1, &res(0), 1 );
    cout << "CCoor - SBLAS = " << normi(res) << '\n';

    F77NAME(amux)( (int*)&n, &v(0), &res(0), &A(0), &J(0), &R(0) );
    res -= ccoor * v;
    cout << "CCoor - SPARSKIT = " << normi(res) << '\n';

    integer cicle_repeat{ 1 + 1000000 / nr };
    test_CRow (cicle_repeat, mHandle, crow,  v);
    test_CCol (cicle_repeat, mHandle, ccol,  v);
    test_CCoor(cicle_repeat, mHandle, ccoor, v);
  }

  cout << "\nTEST_TIME ALL DONE\n\n";

}
