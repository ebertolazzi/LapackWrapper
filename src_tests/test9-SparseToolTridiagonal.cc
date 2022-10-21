#define SPARSETOOL_DEBUG
#include <sparse_tool/sparse_tool.hh>
#include <sparse_tool/sparse_tool_extra.hh>
#include <sparse_tool/sparse_tool_iterative.hh>
#include <sparse_tool/sparse_tool_matrix_market.hh>

//#include <sparse_tool/interfaces/MA41.hh>
//#include <sparse_tool/interfaces/MA48.hh>
#include <sparse_tool/interfaces/SuperLU.hh>
//#include <sparse_tool/interfaces/mkl_pardiso.hh>
//#include <sparse_tool/interfaces/UMF.hh>

#include <fstream>
#include <iostream>
#include <cmath>

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshorten-64-to-32"
#endif

#include "Utils.hh"

#ifdef __clang__
#pragma clang diagnostic pop
#endif

using namespace Sparse_tool_load;
using namespace std;

using real_type = double;
using integer   = int;

static Utils::Console msg(&std::cout);

class all_trid_test {

  Vector<real_type>     a, b, RES, RRES;
  TridMatrix<real_type> A;
  real_type             scalar;

  void
  test_diff( Vector<real_type> const & x, Vector<real_type> const & y ) {
    real_type err = (x-y).template lpNorm<Eigen::Infinity>();
    test_diff(err);
  }

  void
  test_diff( real_type const & err ) {
    if ( err < 1e-6 ) {
      cout << " OK!\n\n";
    } else {
      cout << " ************** ERR =\n\n";
      for ( integer i = 0; i < RES.size(); ++i )
        cout << i << " err = "
             << setw(15) << setprecision(10) << RES[i] - RRES[i] << " "
             << setw(15) << setprecision(10) << RES[i] << ' '
             << setw(15) << setprecision(10) << RRES[i] << '\n';
      exit(0);
    }
  }

public:

  all_trid_test( integer n ) {
    cout << "all_trid_test_allocating..." << flush;
    A.resize( n );
    a.resize( n );
    b.resize( n );
    RES.resize( n );
    RRES.resize( n );
    cout << "done\n";
  }

  ~all_trid_test() { }

  void
  preco(void) {
    Sparse_tool::integer i;
    for ( i = 0; i < a.size(); ++i ) a[i] = 1.2+i*sqrt(2.0);
    for ( i = 0; i < b.size(); ++i ) b[i] = 1+i*sqrt(3.0);

    A = a;
    for ( i = 1; i+1 < A.nrows(); ++i )
      A(i,i-1) = A(i,i+1) = 1;
    scalar = sqrt(2.0);
  }

  void do_all_tests(void) {
    cout << "*******************************\n"
         << "*******************************\n";

    test001();
    test002();
    test003();
    test004();

    cout << "ALL DONE\n";

  }

  void test001();
  void test002();
  void test003();
  void test004();

};

void
all_trid_test::test001(void) {
  preco();
  cout << "\ntest(1)\n";
  Sparse_tool::integer N(A.nrows());

  RES = A * a;
  for ( Sparse_tool::integer i = 0; i < N; ++i ) {
    RRES[i] = A(i,i) * a[i];
    if ( i > 0   ) RRES[i] += A(i,i-1) * a[i-1];
    if ( i < N-1 ) RRES[i] += A(i,i+1) * a[i+1];
  }
  test_diff(RES, RRES);
}

void
all_trid_test::test002(void) {
  preco();
  cout << "\ntest(2)\n";
  Sparse_tool::integer N(A.nrows());
  RES = b + A * a;
  for ( Sparse_tool::integer i = 0; i < N; ++i ) {
    RRES[i] = b[i] + A(i,i) * a[i];
    if ( i > 0   ) RRES[i] += A(i,i-1) * a[i-1];
    if ( i < N-1 ) RRES[i] += A(i,i+1) * a[i+1];
  }
  test_diff(RES, RRES);
}

void
all_trid_test::test003(void) {
  preco();
  cout << "\ntest(3)\n";
  Sparse_tool::integer N(A.nrows());
  RES = b - A * a;
  for ( Sparse_tool::integer i = 0; i < N; ++i ) {
    RRES[i] = b[i] - A(i,i) * a[i];
    if ( i > 0   ) RRES[i] -= A(i,i-1) * a[i-1];
    if ( i < N-1 ) RRES[i] -= A(i,i+1) * a[i+1];
  }
  test_diff(RES, RRES);
}

void
all_trid_test::test004(void) {
  preco();
  cout << "\ntest(4)\n";
  b   = A * a;
  RES = b / A;
  test_diff(RES, a);
}

//********************************************************************
//********************************************************************
//********************************************************************
//********************************************************************

int
main() {
  try {
    all_trid_test allm(4);
    allm.do_all_tests();
  } catch ( exception const & exc ) {
    msg.error( exc.what() );
  } catch ( ... ) {
    msg.error("Errore Sconosciuto!\n");
  }
  cout << "\nTEST_TRID ALL DONE\n\n";
  return 0;
}
