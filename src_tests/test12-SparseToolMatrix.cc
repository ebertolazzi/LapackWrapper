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

using Real = double;

class all_matrix_test {
  Vector<Real>      a, b, c, res, rres;
  Real              scalar{1};
  CRowMatrix<Real>  crow;
  CColMatrix<Real>  ccol;
  CCoorMatrix<Real> ccoor;

  void
  test_diff(Vector<Real> const & x, Vector<Real> const & y) {
    Real err{0};
    for ( integer i{0}; i < x.size(); ++i ) {
      Real bf = x[i] - y[i];
      err += bf > 0 ? bf : -bf;
    }
    test_diff(err);
  }

  void
  test_diff( Real const & err ) {
    if ( err < 1e-10 ) {
      cout << " OK! " << err << '\n';
    } else {
      cout << " ************** ERR = " << err << '\n';
      for ( integer i{0}; i < res.size(); ++i )
        cout << i << " err = "
             << setw(15) << setprecision(10) << res[i] - rres[i] << " "
             << setw(15) << setprecision(10) << res[i] << " "
             << setw(15) << setprecision(10) << rres[i] << '\n';
      exit(0);
    }
  }

  void
  add_S_mul_M_mul_V( Vector<Real> & Res, Real const s, Vector<Real> const & av ) {
    for ( integer k{0}; k < ccoor.nnz(); ++k ) {
      integer i{ ccoor.getI()(k) };
      integer j{ ccoor.getJ()(k) };
      Real    v{ ccoor.getA()(k) };
      Res(i) += s * v * av(j);
    }
  }

  void
  add_S_mul_Mt_mul_V( Vector<Real> & Res, Real const s, Vector<Real> const & av ) {
    for ( integer k{0}; k < ccoor.nnz(); ++k ) {
      integer i{ ccoor.getI()(k) };
      integer j{ ccoor.getJ()(k) };
      Real    v{ ccoor.getA()(k) };
      Res(j) += s * v * av(i);
    }
  }

public:

  all_matrix_test() : a(10), b(10), c(10), res(10), rres(10) {
    ccoor.resize( 10, 10 );
    for ( integer i{0}; i < 10; ++i ) {
      ccoor.insert(i,i)        = 10;
      ccoor.insert(i,(i+4)%10) = -1;
      ccoor.insert((i+4)%10,i) = -1.3;
      ccoor.insert(i,(i+7)%10) = -2;
      ccoor.insert((i+7)%10,i) = -2.3;
    }
    ccoor.internal_order();
    crow = ccoor;
    ccol = ccoor;
  }

  ~all_matrix_test() {}

  void preco(void) {
    integer i;
    for ( i = 0; i < a.size(); ++i ) a[i] = 1.2+i*sqrt(2.0);
    for ( i = 0; i < b.size(); ++i ) b[i] = 1+i*sqrt(3.0);
    c = a;
    scalar = sqrt(2.0);
  }

  void do_all_tests(void) {
    cout << "*******************************\n"
         << "*******************************\n"
         << "MATRIX TESTS\n";

    test001(); test002(); test003(); test004(); test005();
    test006(); test007(); test008(); test009(); test010();

    cout << "MATRIX TESTS: ALL DONE\n\n\n";

  }

  void test001();
  void test002();
  void test003();
  void test004();
  void test005();
  void test006();
  void test007();
  void test008();
  void test009();
  void test010();

  void test011();
  void test012();
  void test013();
  void test014();
  void test015();
  void test016();
  void test017();
  void test018();
  void test019();
  void test020();

  void test021();
  void test022();
  void test023();
  void test024();
  void test025();
  void test026();
  void test027();
  void test028();
  void test029();
  void test030();

  void test031();
  void test032();
  void test033();
  void test034();
  void test035();
  void test036();
  void test037();
  void test038();
  void test039();
  void test040();

  void test041();
  void test042();
  void test043();
  void test044();
  void test045();
  void test046();
  void test047();
  void test048();
  void test049();
  void test050();

  void test051();
  void test052();
  void test053();
  void test054();
  void test055();
  void test056();
  void test057();
  void test058();
  void test059();
  void test060();

  void test061();
  void test062();

};

void all_matrix_test::test001(void) {
  preco();

  cout << "test(1)\n";

  cout << "vector =  ccoor * vector";
  res  = ccoor * a;
  rres = 0;
  add_S_mul_M_mul_V( rres, 1.0, a );
  test_diff(res, rres);

  cout << "vector =  ccol  * vector";
  res  = ccol * a;
  rres = 0;
  add_S_mul_M_mul_V( rres, 1.0, a );
  test_diff(res, rres);

  cout << "vector =  crow  * vector";
  res  = crow * a;
  rres = 0;
  add_S_mul_M_mul_V( rres, 1.0, a );
  test_diff(res, rres);

  //////////////////////////////////

  cout << "vector += ccoor * vector";
  res   = 1;
  res  += ccoor * a;
  rres  = 1;
  add_S_mul_M_mul_V( rres, 1.0, a );
  test_diff(res, rres);

  cout << "vector += ccol  * vector";
  res  = 1;
  res += ccol * a;
  rres = 1;
  add_S_mul_M_mul_V( rres, 1.0, a );
  test_diff(res, rres);

  cout << "vector += crow  * vector";
  res  = 1;
  res += crow * a;
  rres = 1;
  add_S_mul_M_mul_V( rres, 1.0, a );
  test_diff(res, rres);

  //////////////////////////////////

  cout << "vector -= ccoor * vector";
  res   = 1;
  res  -= ccoor * a;
  rres  = 1;
  add_S_mul_M_mul_V( rres, -1.0, a );
  test_diff(res, rres);

  cout << "vector -= ccol  * vector";
  res  = 1;
  res -= ccol * a;
  rres = 1;
  add_S_mul_M_mul_V( rres, -1.0, a );
  test_diff(res, rres);

  cout << "vector -= crow  * vector";
  res  = 1;
  res -= crow * a;
  rres = 1;
  add_S_mul_M_mul_V( rres, -1.0, a );
  test_diff(res, rres);

}

void all_matrix_test::test002(void) {
  preco();

  cout << "test(2)\n";

  cout << "vector =  ccoor ^ vector";
  res  = ccoor ^ a;
  rres = 0;
  add_S_mul_Mt_mul_V( rres, 1.0, a );
  test_diff(res, rres);

  cout << "vector =  ccol  ^ vector";
  res  = ccol ^ a;
  rres = 0;
  add_S_mul_Mt_mul_V( rres, 1.0, a );
  test_diff(res, rres);

  cout << "vector =  crow  ^ vector";
  res  = crow ^ a;
  rres = 0;
  add_S_mul_Mt_mul_V( rres, 1.0, a );
  test_diff(res, rres);

  //////////////////////////////////

  cout << "vector += ccoor ^ vector";
  res   = 1;
  res  += ccoor ^ a;
  rres  = 1;
  add_S_mul_Mt_mul_V( rres, 1.0, a );
  test_diff(res, rres);

  cout << "vector += ccol  ^ vector";
  res  = 1;
  res += ccol ^ a;
  rres = 1;
  add_S_mul_Mt_mul_V( rres, 1.0, a );
  test_diff(res, rres);

  cout << "vector += crow  ^ vector";
  res  = 1;
  res += crow ^ a;
  rres = 1;
  add_S_mul_Mt_mul_V( rres, 1.0, a );
  test_diff(res, rres);

  //////////////////////////////////

  cout << "vector -= ccoor ^ vector";
  res   = 1;
  res  -= ccoor ^ a;
  rres  = 1;
  add_S_mul_Mt_mul_V( rres, -1.0, a );
  test_diff(res, rres);

  cout << "vector -= ccol  ^ vector";
  res  = 1;
  res -= ccol ^ a;
  rres = 1;
  add_S_mul_Mt_mul_V( rres, -1.0, a );
  test_diff(res, rres);

  cout << "vector -= crow  ^ vector";
  res  = 1;
  res -= crow ^ a;
  rres = 1;
  add_S_mul_Mt_mul_V( rres, -1.0, a );
  test_diff(res, rres);

}


void all_matrix_test::test003(void) {
  preco();

  cout << "test(3)\n";

  cout << "vector =  scalar * ccoor * vector";
  res  = scalar * ccoor * a;
  rres = 0;
  add_S_mul_M_mul_V( rres, scalar, a );
  test_diff(res, rres);

  cout << "vector =  scalar * ccol  * vector";
  res  = scalar * ccol * a;
  rres = 0;
  add_S_mul_M_mul_V( rres, scalar, a );
  test_diff(res, rres);

  cout << "vector =  scalar * crow  * vector";
  res  = scalar * crow * a;
  rres = 0;
  add_S_mul_M_mul_V( rres, scalar, a );
  test_diff(res, rres);

  //////////////////////////////////

  cout << "vector += scalar * ccoor * vector";
  res   = 1;
  res  += scalar * ccoor * a;
  rres  = 1;
  add_S_mul_M_mul_V( rres, scalar, a );
  test_diff(res, rres);

  cout << "vector += scalar * ccol  * vector";
  res  = 1;
  res += scalar * ccol * a;
  rres = 1;
  add_S_mul_M_mul_V( rres, scalar, a );
  test_diff(res, rres);

  cout << "vector += scalar * crow  * vector";
  res  = 1;
  res += scalar * crow * a;
  rres = 1;
  add_S_mul_M_mul_V( rres, scalar, a );
  test_diff(res, rres);

  //////////////////////////////////

  cout << "vector -= scalar * ccoor * vector";
  res   = 1;
  res  -= scalar * ccoor * a;
  rres  = 1;
  add_S_mul_M_mul_V( rres, -scalar, a );
  test_diff(res, rres);

  cout << "vector -= scalar * ccol  * vector";
  res  = 1;
  res -= scalar * ccol * a;
  rres = 1;
  add_S_mul_M_mul_V( rres, -scalar, a );
  test_diff(res, rres);

  cout << "vector -= scalar * crow  * vector";
  res  = 1;
  res -= scalar * crow * a;
  rres = 1;
  add_S_mul_M_mul_V( rres, -scalar, a );
  test_diff(res, rres);

}


void all_matrix_test::test004(void) {
  preco();

  cout << "test(4)\n";

  cout << "vector =  scalar * ccoor ^ vector";
  res  = scalar * ccoor ^ a;
  rres = 0;
  add_S_mul_Mt_mul_V( rres, scalar, a );
  test_diff(res, rres);

  cout << "vector =  scalar * ccol  ^ vector";
  res  = scalar * ccol ^ a;
  rres = 0;
  add_S_mul_Mt_mul_V( rres, scalar, a );
  test_diff(res, rres);

  cout << "vector =  scalar * crow  ^ vector";
  res  = scalar * crow ^ a;
  rres = 0;
  add_S_mul_Mt_mul_V( rres, scalar, a );
  test_diff(res, rres);

  //////////////////////////////////

  cout << "vector += scalar * ccoor ^ vector";
  res   = 1;
  res  += scalar * ccoor ^ a;
  rres  = 1;
  add_S_mul_Mt_mul_V( rres, scalar, a );
  test_diff(res, rres);

  cout << "vector += scalar * ccol  ^ vector";
  res  = 1;
  res += scalar * ccol ^ a;
  rres = 1;
  add_S_mul_Mt_mul_V( rres, scalar, a );
  test_diff(res, rres);

  cout << "vector += scalar * crow  ^ vector";
  res  = 1;
  res += scalar * crow ^ a;
  rres = 1;
  add_S_mul_Mt_mul_V( rres, scalar, a );
  test_diff(res, rres);

  //////////////////////////////////

  cout << "vector -= scalar * ccoor ^ vector";
  res   = 1;
  res  -= scalar * ccoor ^ a;
  rres  = 1;
  add_S_mul_Mt_mul_V( rres, -scalar, a );
  test_diff(res, rres);

  cout << "vector -= scalar * ccol  ^ vector";
  res  = 1;
  res -= scalar * ccol ^ a;
  rres = 1;
  add_S_mul_Mt_mul_V( rres, -scalar, a );
  test_diff(res, rres);

  cout << "vector -= scalar * crow  ^ vector";
  res  = 1;
  res -= scalar * crow ^ a;
  rres = 1;
  add_S_mul_Mt_mul_V( rres, -scalar, a );
  test_diff(res, rres);

}



void all_matrix_test::test005(void) {
  preco();

  cout << "test(5)\n";

  cout << "vector =  vector + ccoor * vector";
  res  = b + ccoor * a;
  rres = b;
  add_S_mul_M_mul_V( rres, +1, a );
  test_diff(res, rres);

  cout << "vector =  vector + ccol  * vector";
  res  = b +  ccol * a;
  rres = b;
  add_S_mul_M_mul_V( rres, +1, a );
  test_diff(res, rres);

  cout << "vector =  vector + crow  * vector";
  res  = b + crow * a;
  rres = b;
  add_S_mul_M_mul_V( rres, +1, a );
  test_diff(res, rres);

  //////////////////////////////////

  cout << "vector += vector + ccoor * vector";
  res   = 1;
  res  += b + ccoor * a;
  rres  = 1.0 + b.array();
  add_S_mul_M_mul_V( rres, +1, a );
  test_diff(res, rres);

  cout << "vector += vector + ccol  * vector";
  res  = 1;
  res += b + ccol * a;
  rres = 1.0 + b.array();
  add_S_mul_M_mul_V( rres, +1, a );
  test_diff(res, rres);

  cout << "vector += vector + crow  * vector";
  res  = 1;
  res += b + crow * a;
  rres = 1.0 + b.array();
  add_S_mul_M_mul_V( rres, +1, a );
  test_diff(res, rres);

  //////////////////////////////////

  cout << "vector -= vector + ccoor * vector";
  res   = 1;
  res  -= b + ccoor * a;
  rres  = 1.0 - b.array();
  add_S_mul_M_mul_V( rres, -1, a );
  test_diff(res, rres);

  cout << "vector -= vector + ccol  * vector";
  res  = 1;
  res -= b + ccol * a;
  rres = 1.0 - b.array();
  add_S_mul_M_mul_V( rres, -1, a );
  test_diff(res, rres);

  cout << "vector -= vector + crow  * vector";
  res  = 1;
  res -= b + crow * a;
  rres = 1.0 - b.array();
  add_S_mul_M_mul_V( rres, -1, a );
  test_diff(res, rres);

}


void all_matrix_test::test006(void) {
  preco();

  cout << "test(6)\n";

  cout << "vector =  vector + (ccoor ^ vector)";
  res  = b + (ccoor ^ a);
  rres = b;
  add_S_mul_Mt_mul_V( rres, +1, a );
  test_diff(res, rres);

  cout << "vector =  vector + (ccol  ^ vector)";
  res  = b +  (ccol ^ a);
  rres = b;
  add_S_mul_Mt_mul_V( rres, +1, a );
  test_diff(res, rres);

  cout << "vector =  vector + (crow  ^ vector)";
  res  = b + (crow ^ a);
  rres = b;
  add_S_mul_Mt_mul_V( rres, +1, a );
  test_diff(res, rres);

  //////////////////////////////////

  cout << "vector += vector + (ccoor ^ vector)";
  res   = 1;
  res  += b + (ccoor ^ a);
  rres  = 1.0 + b.array();
  add_S_mul_Mt_mul_V( rres, +1, a );
  test_diff(res, rres);

  cout << "vector += vector + (ccol  ^ vector)";
  res  = 1;
  res += b + (ccol ^ a);
  rres = 1.0 + b.array();
  add_S_mul_Mt_mul_V( rres, +1, a );
  test_diff(res, rres);

  cout << "vector += vector + (crow  ^ vector)";
  res  = 1;
  res += b + (crow ^ a);
  rres = 1.0 + b.array();
  add_S_mul_Mt_mul_V( rres, +1, a );
  test_diff(res, rres);

  //////////////////////////////////

  cout << "vector -= vector + (ccoor ^ vector)";
  res   = 1;
  res  -= b + (ccoor ^ a);
  rres  = 1.0 - b.array();
  add_S_mul_Mt_mul_V( rres, -1, a );
  test_diff(res, rres);

  cout << "vector -= vector + (ccol  ^ vector)";
  res  = 1;
  res -= b + (ccol ^ a);
  rres = 1.0 - b.array();
  add_S_mul_Mt_mul_V( rres, -1, a );
  test_diff(res, rres);

  cout << "vector -= vector + (crow  ^ vector)";
  res  = 1;
  res -= b + (crow ^ a);
  rres = 1.0 - b.array();
  add_S_mul_Mt_mul_V( rres, -1, a );
  test_diff(res, rres);

}



void all_matrix_test::test007(void) {
  preco();

  cout << "test(7)\n";

  cout << "vector =  vector - ccoor * vector";
  res  = b - ccoor * a;
  rres = b;
  add_S_mul_M_mul_V( rres, -1, a );
  test_diff(res, rres);

  cout << "vector =  vector - ccol  * vector";
  res  = b - ccol * a;
  rres = b;
  add_S_mul_M_mul_V( rres, -1, a );
  test_diff(res, rres);

  cout << "vector =  vector - crow  * vector";
  res  = b - crow * a;
  rres = b;
  add_S_mul_M_mul_V( rres, -1, a );
  test_diff(res, rres);

  //////////////////////////////////

  cout << "vector += vector - ccoor * vector";
  res   = 1;
  res  += b - ccoor * a;
  rres  = 1.0 + b.array();
  add_S_mul_M_mul_V( rres, -1, a );
  test_diff(res, rres);

  cout << "vector += vector - ccol  * vector";
  res  = 1;
  res += b - ccol * a;
  rres = 1.0 + b.array();
  add_S_mul_M_mul_V( rres, -1, a );
  test_diff(res, rres);

  cout << "vector += vector - crow  * vector";
  res  = 1;
  res += b - crow * a;
  rres = 1.0 + b.array();
  add_S_mul_M_mul_V( rres, -1, a );
  test_diff(res, rres);

  //////////////////////////////////

  cout << "vector -= vector - ccoor * vector";
  res   = 1;
  res  -= b - ccoor * a;
  rres  = 1.0 - b.array();
  add_S_mul_M_mul_V( rres, +1, a );
  test_diff(res, rres);

  cout << "vector -= vector - ccol  * vector";
  res  = 1;
  res -= b - ccol * a;
  rres = 1.0 - b.array();
  add_S_mul_M_mul_V( rres, +1, a );
  test_diff(res, rres);

  cout << "vector -= vector - crow  * vector";
  res  = 1;
  res -= b - crow * a;
  rres = 1.0 - b.array();
  add_S_mul_M_mul_V( rres, +1, a );
  test_diff(res, rres);

}


void all_matrix_test::test008(void) {
  preco();

  cout << "test(8)\n";

  cout << "vector =  vector - (ccoor ^ vector)";
  res  = b - (ccoor ^ a);
  rres = b;
  add_S_mul_Mt_mul_V( rres, -1, a );
  test_diff(res, rres);

  cout << "vector =  vector - (ccol  ^ vector)";
  res  = b - (ccol ^ a);
  rres = b;
  add_S_mul_Mt_mul_V( rres, -1, a );
  test_diff(res, rres);

  cout << "vector =  vector - (crow  ^ vector)";
  res  = b - (crow ^ a);
  rres = b;
  add_S_mul_Mt_mul_V( rres, -1, a );
  test_diff(res, rres);

  //////////////////////////////////

  cout << "vector += vector - (ccoor ^ vector)";
  res   = 1;
  res  += b - (ccoor ^ a);
  rres  = 1.0 + b.array();
  add_S_mul_Mt_mul_V( rres, -1, a );
  test_diff(res, rres);

  cout << "vector += vector - (ccol  ^ vector)";
  res  = 1;
  res += b - (ccol ^ a);
  rres = 1.0 + b.array();
  add_S_mul_Mt_mul_V( rres, -1, a );
  test_diff(res, rres);

  cout << "vector += vector - (crow  ^ vector)";
  res  = 1;
  res += b - (crow ^ a);
  rres = 1.0 + b.array();
  add_S_mul_Mt_mul_V( rres, -1, a );
  test_diff(res, rres);

  //////////////////////////////////

  cout << "vector -= vector - (ccoor ^ vector)";
  res   = 1;
  res  -= b - (ccoor ^ a);
  rres  = 1.0 - b.array();
  add_S_mul_Mt_mul_V( rres, +1, a );
  test_diff(res, rres);

  cout << "vector -= vector - (ccol  ^ vector)";
  res  = 1;
  res -= b - (ccol ^ a);
  rres = 1.0 - b.array();
  add_S_mul_Mt_mul_V( rres, +1, a );
  test_diff(res, rres);

  cout << "vector -= vector - (crow  ^ vector)";
  res  = 1;
  res -= b - (crow ^ a);
  rres = 1.0 - b.array();
  add_S_mul_Mt_mul_V( rres, +1, a );
  test_diff(res, rres);

}


void all_matrix_test::test009(void) {
  preco();

  cout << "test(9)\n";

  cout << "vector =  vector + scalar * ccoor * vector";
  res  = b + scalar * ccoor * a;
  rres = b;
  add_S_mul_M_mul_V( rres, scalar, a );
  test_diff(res, rres);

  cout << "vector =  vector + scalar * ccol  * vector";
  res  = b + scalar * ccol * a;
  rres = b;
  add_S_mul_M_mul_V( rres, scalar, a );
  test_diff(res, rres);

  cout << "vector =  vector + scalar * crow  * vector";
  res  = b + scalar * crow * a;
  rres = b;
  add_S_mul_M_mul_V( rres, scalar, a );
  test_diff(res, rres);

  //////////////////////////////////

  cout << "vector += vector + scalar * ccoor * vector";
  res   = 1;
  res  += b + scalar * ccoor * a;
  rres  = 1.0 + b.array();
  add_S_mul_M_mul_V( rres, scalar, a );
  test_diff(res, rres);

  cout << "vector += vector + scalar * ccol  * vector";
  res  = 1;
  res += b + scalar * ccol * a;
  rres = 1.0 + b.array();
  add_S_mul_M_mul_V( rres, scalar, a );
  test_diff(res, rres);

  cout << "vector += vector + scalar * crow  * vector";
  res  = 1;
  res += b + scalar * crow * a;
  rres = 1.0 + b.array();
  add_S_mul_M_mul_V( rres, scalar, a );
  test_diff(res, rres);

  //////////////////////////////////

  cout << "vector -= vector + scalar * ccoor * vector";
  res   = 1;
  res  -= b + scalar * ccoor * a;
  rres  = 1.0 - b.array();
  add_S_mul_M_mul_V( rres, -scalar, a );
  test_diff(res, rres);

  cout << "vector -= vector + scalar * ccol  * vector";
  res  = 1;
  res -= b + scalar * ccol * a;
  rres = 1.0 - b.array();
  add_S_mul_M_mul_V( rres, -scalar, a );
  test_diff(res, rres);

  cout << "vector -= vector + scalar * crow  * vector";
  res  = 1;
  res -= b + scalar * crow * a;
  rres = 1.0 - b.array();
  add_S_mul_M_mul_V( rres, -scalar, a );
  test_diff(res, rres);

}



void all_matrix_test::test010(void) {
  preco();

  cout << "test(10)\n";

  cout << "vector =  vector + scalar * (ccoor ^ vector)";
  res  = b + scalar * (ccoor ^ a);
  rres = b;
  add_S_mul_Mt_mul_V( rres, scalar, a );
  test_diff(res, rres);

  cout << "vector =  vector + scalar * (ccol  ^ vector)";
  res  = b + scalar * (ccol ^ a);
  rres = b;
  add_S_mul_Mt_mul_V( rres, scalar, a );
  test_diff(res, rres);

  cout << "vector =  vector + scalar * (crow  ^ vector)";
  res  = b + scalar * (crow ^ a);
  rres = b;
  add_S_mul_Mt_mul_V( rres, scalar, a );
  test_diff(res, rres);

  //////////////////////////////////

  cout << "vector += vector + scalar * (ccoor ^ vector)";
  res   = 1;
  res  += b + scalar * (ccoor ^ a);
  rres  = 1.0 + b.array();
  add_S_mul_Mt_mul_V( rres, scalar, a );
  test_diff(res, rres);

  cout << "vector += vector + scalar * (ccol  ^ vector)";
  res  = 1;
  res += b + scalar * (ccol ^ a);
  rres = 1.0 + b.array();
  add_S_mul_Mt_mul_V( rres, scalar, a );
  test_diff(res, rres);

  cout << "vector += vector + scalar * (crow  ^ vector)";
  res  = 1;
  res += b + scalar * (crow ^ a);
  rres = 1.0 + b.array();
  add_S_mul_Mt_mul_V( rres, scalar, a );
  test_diff(res, rres);

  //////////////////////////////////

  cout << "vector -= vector + scalar * (ccoor ^ vector)";
  res   = 1;
  res  -= b + scalar * (ccoor ^ a);
  rres  = 1.0 - b.array();
  add_S_mul_Mt_mul_V( rres, -scalar, a );
  test_diff(res, rres);

  cout << "vector -= vector + scalar * (ccol  ^ vector)";
  res  = 1;
  res -= b + scalar * (ccol ^ a);
  rres = 1.0 - b.array();
  add_S_mul_Mt_mul_V( rres, -scalar, a );
  test_diff(res, rres);

  cout << "vector -= vector + scalar * (crow  ^ vector)";
  res  = 1;
  res -= b + scalar * (crow ^ a);
  rres = 1.0 - b.array();
  add_S_mul_Mt_mul_V( rres, -scalar, a );
  test_diff(res, rres);

}

//********************************************************************
//********************************************************************
//********************************************************************
//********************************************************************

#define LOOP(N)   for ( integer i{0}; i < N; ++i )
#define REPEAT(N) for ( integer iii{0}; iii < N; ++iii )

using Real = double;

static
void
out(
  integer const   N,
  integer const   nnz,
  Real    const & ta,
  Real    const & tb
) {
  fmt::print(
    "N={:3} TIME ({:.5}:{:.5}) FP({:.5}:{:.5}) (sparselib:hand made)\n",
    N, ta, tb, 1e-3*N*nnz/ta, 1e-3*N*nnz/tb
  );
}

// test CRow  ********************************************************

static
void
test_CRow(
  integer                  N,
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
    integer const * R { crow.getR().data() };
    integer const * J { crow.getJ().data() };
    Real    const * A { crow.getA().data() };

    for ( integer i{0}; i < crow.nrows(); ++i, ++R ) {
      integer n{ R[1]-R[0] };
      Real tmp{0};
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
  Real timea{ tm.elapsed_ms() };

  tm.tic();
  REPEAT(N) {
    integer const * I { ccol.getI().data() };
    integer const * C { ccol.getC().data() };
    Real    const * A { ccol.getA().data() };

    res = 0;
    for ( integer i{0}; i < ccol.ncols(); ++i, ++C ) {
      Real    tmp { v(i) };
      integer n   { C[1]-C[0] };
      while ( n-- > 0 ) res(*I++) += *A++ * tmp;
    }
  }
  tm.toc();
  Real timeb = tm.elapsed_ms();

  cout << "CColMatrix:  ";
  out( N, ccol.nnz(), timea, timeb );
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
  Real timea{ tm.elapsed_ms() };

  tm.tic();
  REPEAT(N) {
    integer const * I { ccoor.getI().data() };
    integer const * J { ccoor.getJ().data() };
    Real    const * A { ccoor.getA().data() };
    res = 0;
    for ( integer i{0}; i < ccoor.nnz(); ++i )
      res(*I++) += *A++ * v(*J++);
  }
  tm.toc();
  Real timeb{ tm.elapsed_ms() };

  cout << "CCoorMatrix: ";
  out( N, ccoor.nnz(), timea, timeb );
}

static
void
test_timing() {
  Utils::TicToc tm;

  constexpr integer n{10000};

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

  char const *rMatrix[] = {
    "hor__131.mtx",
    "af23560.mtx",
    "plat1919.mtx",
    nullptr
  };
  char const **p = rMatrix;

  Vector<Real> v;

  for (; *p != nullptr; ++p ) {

    fmt::print( "TEST matrix = {}\n", *p );
    CCoorMatrix<Real> ccoor;
    CRowMatrix<Real>  crow;
    CColMatrix<Real>  ccol;

    MatrixMarket mm;
    cout << "load matrix..." << flush;
    mm.read( fmt::format( "mm/{}", *p ), ccoor );
    cout << "done\n";
    mm.info(cout);

    Spy( fmt::format("mm/{}.eps", *p), ccoor, 15.0 );

    tm.tic();
    crow = ccoor;
    tm.toc();
    Real timea = tm.elapsed_ms();

    tm.tic();
    ccol = ccoor;
    tm.toc();
    Real timeb = tm.elapsed_ms();

    fmt::print( "TO CROW {}[ms]\n", timea );
    fmt::print( "TO CCOL {}[ms]\n", timeb );

    integer const nr{ ccoor.nrows() };

    v.resize(nr);
    res.resize(nr);

    v = 1;
    v[nr/2] = 234;
    v[nr/4] = -3;
    v[(nr*3)/4 ] = -12;

    res = ccoor * v;
    res = res - crow * v;
    fmt::print( "CCoor - CRow = {}\n", res.lpNorm<Eigen::Infinity>() );

    res = ccoor * v;
    res = res - ccol * v;
    fmt::print( "CCoor - CCol = {}\n", res.lpNorm<Eigen::Infinity>() );

    integer cicle_repeat{1 + 400000 / nr};
    test_CRow (cicle_repeat, crow,  v);
    test_CCol (cicle_repeat, ccol,  v);
    test_CCoor(cicle_repeat, ccoor, v);
  }

  cout << "\nTEST_TIME ALL DONE\n\n";
}

//********************************************************************
//********************************************************************
//********************************************************************
//********************************************************************

int
main() {
  try {
    all_matrix_test allm;
    allm.do_all_tests();
    test_timing();
  } catch ( exception const & exc ) {
    msg.error( exc.what() );
  } catch ( ... ) {
    msg.error("Errore Sconosciuto!\n");
  }
  cout << "\nALL DONE FOLKS!\n\n";
  return 0;
}
