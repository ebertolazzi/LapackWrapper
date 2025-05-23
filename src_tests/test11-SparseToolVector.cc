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

//#include <sparse_tool/interfaces/MA41.hh>
//#include <sparse_tool/interfaces/MA48.hh>
//#include <sparse_tool/interfaces/SuperLU.hh>
//#include <sparse_tool/interfaces/mkl_pardiso.hh>
//#include <sparse_tool/interfaces/UMF.hh>

using namespace ::Sparse_tool_load;
using           ::Sparse_tool::integer;
using namespace ::std;

using Real = double;

static Utils::Console msg(&std::cout);

inline
Real
power2( Real const a )
{ return a*a; }

class all_vector_test {

  Vector<Real> a, b, c, res, rres;
  Real         scalar{1};

  void
  test_diff(Vector<Real> const & x, Vector<Real> const & y) {
    Real err{ 0 };
    for ( integer i{0}; i < x.size(); ++i ) {
      Real const bf{ x[i] - y[i] };
      err += bf > 0 ? bf : -bf;
    }
    test_diff(err);
  }

  void
  test_diff(Real const & err) {
    if ( err < 1e-10 ) {
      cout << " OK! " << err << '\n';
    } else {
      cout << " ************** ERR = " << err << '\n';
      for ( integer i{0}; i < res.size(); ++i )
        cout << i << " err = "
             << setw(15) << setprecision(10) << res[i] - rres[i] << ' '
             << setw(15) << setprecision(10) << res[i]  << ' '
             << setw(15) << setprecision(10) << rres[i] << '\n';
      exit(0);
    }
  }

public:

  explicit
  all_vector_test( integer const n ) : a(n), b(n), c(n), res(n), rres(n) { }

  ~all_vector_test() = default;

  void
  preco() {
    integer adim{ static_cast<integer>(a.size()) };
    integer bdim{ static_cast<integer>(b.size()) };
    for ( integer i{0}; i < adim; ++i ) a[i] = 1.2+i*sqrt(2.0);
    for ( integer i{0}; i < bdim; ++i ) b[i] = 1+i*sqrt(3.0);
    c = a;
    scalar = sqrt(2.0);
  }

  void
  do_all_tests() {
    cout << "*******************************\n"
         << "*******************************\n"
         << "VECTOR TESTS\n";

    test001(); test002(); test003(); test004(); test005();
    test006(); test007(); test008(); test009(); test010();

    test011(); test012(); test013(); test014(); test015();
    test016(); test017(); test018(); test019(); test020();

    test021(); test022(); test023(); test024(); test025();
    test026(); test027(); test028(); test029(); test030();

    test031(); test032(); test033(); test034(); test035();
    test036(); test037(); test038(); test039(); test030();

    test041(); test042(); test043(); test044(); test045();
    test046(); test047(); test048(); test049(); test050();

    test051(); test052(); test053(); test054(); test055();
    test056(); test057(); test057(); test058(); test059();
    test060(); test061(); test062(); test063(); test064();

    cout << "VECTOR TESTS: ALL DONE\n\n\n";

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
  void test063();
  void test064();

};


void
all_vector_test::test001() {
  cout << "test(1)\n";

  scalar = 1.2345;

  cout << "vector = scalar";
  res = scalar;
  for ( integer i{0}; i < a.size(); ++i) rres[i] = scalar;
  test_diff(res, rres);

  preco();
  cout << "vector += scalar";
  res += scalar;
  for ( integer i{0}; i < a.size(); ++i) rres[i] += scalar;
  test_diff(res, rres);

  cout << "vector -= scalar";
  res -= scalar;
  for ( integer i{0}; i < a.size(); ++i) rres[i] -= scalar;
  test_diff(res, rres);

  cout << "vector *= scalar";
  res *= scalar;
  for ( integer i{0}; i < a.size(); ++i) rres[i] *= scalar;
  test_diff(res, rres);

  cout << "vector /= scalar";
  res /= scalar;
  for ( integer i{0}; i < a.size(); ++i) rres[i] /= scalar;
  test_diff(res, rres);

}

#define DO_TESTS(TESTNAME, OP, OPC)                      \
  cout << "vector  = " << TESTNAME;                      \
  res = OP;                                              \
  for ( integer i{0}; i < sz; ++i) rres[i] = OPC;        \
  test_diff(res,rres);                                   \
                                                         \
  cout << "vector += " << TESTNAME;                      \
  res += OP;                                             \
  for ( integer i{0}; i < sz; ++i) rres[i] += OPC;       \
  test_diff(res,rres);                                   \
                                                         \
  cout << "vector -= " << TESTNAME;                      \
  res -= OP;                                             \
  for ( integer i{0}; i < sz; ++i) rres[i] -= OPC;       \
  test_diff(res,rres);                                   \
                                                         \
  cout << "vector *= " << TESTNAME;                      \
  res.array() *= (OP).array();                           \
  for ( integer i{0}; i < sz; ++i) rres[i] *= OPC;       \
  test_diff(res,rres);                                   \
                                                         \
  cout << "vector /= " << TESTNAME;                      \
  res.array() /= (OP).array();                           \
  for ( integer i{0}; i < sz; ++i) rres[i] /= OPC;       \
  test_diff(res,rres)

# define MIN2(A,B)   static_cast<integer>( std::min(A.size(),B.size()) )
# define MIN3(A,B,C) static_cast<integer>( std::min( std::min( A.size(), B.size() ), C.size() ) )

void
all_vector_test::test002() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(2)\n";
  DO_TESTS("vector", a, a[i]);
}


void
all_vector_test::test003() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(3)\n";
  DO_TESTS("- vector", - a, - a[i]);
}


void
all_vector_test::test004() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(4)\n";
  DO_TESTS("vector + vector", a + b, a[i] + b[i]);
}


void
all_vector_test::test005() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(5)\n";
  DO_TESTS("vector - vector", a - b, a[i] - b[i]);
}


void
all_vector_test::test006() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(6)\n";
  DO_TESTS("vector * vector", a.array() * b.array(), a[i] * b[i]);
}


void
all_vector_test::test007() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(7)\n";
  DO_TESTS("vector / vector", a.array() / b.array(), a[i] / b[i]);
}


void
all_vector_test::test008() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(8)\n";
  DO_TESTS("scalar + vector", scalar + a.array(), scalar + a[i]);
}


void
all_vector_test::test009() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(9)\n";
  DO_TESTS("scalar - vector", scalar - a.array(), scalar - a[i]);
}


void
all_vector_test::test010() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(10)\n";
  DO_TESTS("scalar * vector", scalar * a, scalar * a[i]);
}


void
all_vector_test::test011() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(11)\n";
  DO_TESTS("scalar / vector", scalar / a.array(), scalar / a[i]);
}


void
all_vector_test::test012() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(12)\n";
  DO_TESTS("vector + scalar", a.array() + scalar, a[i] + scalar);
}


void
all_vector_test::test013() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(13)\n";
  DO_TESTS("vector - scalar", a.array() - scalar, a[i] - scalar);
}


void
all_vector_test::test014() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(14)\n";
  DO_TESTS("vector * scalar", a * scalar, a[i] * scalar);
}


void
all_vector_test::test015() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(15)\n";
  DO_TESTS("vector / scalar", a / scalar, a[i] / scalar);
}

  // combined test

void
all_vector_test::test016() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(16)\n";
  DO_TESTS("vector + scalar + vector", a.array() + scalar + b.array(), a[i] + scalar + b[i]);
}


void
all_vector_test::test017() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(17)\n";
  DO_TESTS("vector - scalar + vector", a.array() - scalar + b.array(), a[i] - scalar + b[i]);
}


void
all_vector_test::test018() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(18)\n";
  DO_TESTS("vector * scalar + vector", a.array() * scalar + b.array(), a[i] * scalar + b[i]);
}


void
all_vector_test::test019() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(19)\n";
  DO_TESTS("vector / scalar + vector", a / scalar + b, a[i] / scalar + b[i]);
}
  //-----------------

void
all_vector_test::test020() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(20)\n";
  DO_TESTS("vector + scalar - vector", a.array() + scalar - b.array(), a[i] + scalar - b[i]);
}


void
all_vector_test::test021() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(21)\n";
  DO_TESTS("vector - scalar - vector", a.array() - scalar - b.array(), a[i] - scalar - b[i]);
}


void
all_vector_test::test022() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(22)\n";
  DO_TESTS("vector * scalar - vector", a * scalar - b, a[i] * scalar - b[i]);
}


void
all_vector_test::test023() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(23)\n";
  DO_TESTS("vector / scalar - vector", a / scalar - b, a[i] / scalar - b[i]);
}
  //-----------------

void
all_vector_test::test024() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(24)\n";
  DO_TESTS("vector + scalar * vector", a + scalar * b, a[i] + scalar * b[i] );
}


void
all_vector_test::test025() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(25)\n";
  DO_TESTS("vector - scalar * vector", a - scalar * b, a[i] - scalar * b[i]);
}


void
all_vector_test::test026() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(26)\n";
  DO_TESTS("vector * scalar * vector", a.array() * scalar * b.array(), a[i] * scalar * b[i]);
}


void
all_vector_test::test027() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(27)\n";
  DO_TESTS("vector / scalar * vector", a.array() / scalar * b.array(), a[i] / scalar * b[i]);
}
  //-----------------

void
all_vector_test::test028() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(28)\n";
  DO_TESTS("vector + scalar / vector", a.array() + scalar / b.array(), a[i] + scalar / b[i]);
}


void
all_vector_test::test029() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(29)\n";
  DO_TESTS("vector - scalar / vector", a.array() - scalar / b.array(), a[i] - scalar / b[i]);
}


void
all_vector_test::test030() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(30)\n";
  DO_TESTS("vector * scalar / vector", a.array() * scalar / b.array(), a[i] * scalar / b[i]);
}


void
all_vector_test::test031() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(31)\n";
  DO_TESTS("vector / scalar / vector", a.array() / scalar / b.array(), a[i] / scalar / b[i]);
}

// FUNCTIONS

void
all_vector_test::test032() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(32)\n";
  DO_TESTS("abs( vector )", a.array().abs(), abs(a[i]) );
}

void
all_vector_test::test033() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(33)\n";
  DO_TESTS("sin( vector )", a.array().sin(), sin(a[i]) );
}

void
all_vector_test::test034() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(34)\n";
  DO_TESTS("cos( vector )", a.array().cos(), cos(a[i]) );
}

void
all_vector_test::test035() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(35)\n";
  DO_TESTS("tan( vector )", a.array().tan(), tan(a[i]) );
}

void
all_vector_test::test036() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(36)\n";
  DO_TESTS("asin( vector / (1+abs(vector)) )",
	     asin(a.array()/(1.0+abs(a.array()))), asin(a[i]/(1.0+abs(a[i]))) );
}

void
all_vector_test::test037() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(37)\n";
  DO_TESTS("acos( vector / (1+abs(vector)) )",
	     acos(a.array()/(1.0+abs(a.array()))), acos(a[i]/(1.0+abs(a[i]))) );
}

void
all_vector_test::test038() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(38)\n";
  DO_TESTS("atan( vector )", atan(a.array()), atan(a[i]) );
}

void
all_vector_test::test039() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(39)\n";
  DO_TESTS("sinh( vector )", sinh(a.array()), sinh(a[i]) );
}

void
all_vector_test::test040() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(40)\n";
  DO_TESTS("cosh( vector )", cosh(a.array()), cosh(a[i]) );
}

void
all_vector_test::test041() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(41)\n";
  DO_TESTS("tanh( vector )", tanh(a.array()), tanh(a[i]) );
}

void
all_vector_test::test042() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(42)\n";
  DO_TESTS("sqrt( abs(vector) )", sqrt(abs(a.array())), sqrt(abs(a[i])) );
}

void
all_vector_test::test043() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(43)\n";
  DO_TESTS("ceil( vector )", ceil(a.array()), ceil(a[i]) );
}

void
all_vector_test::test044() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(44)\n";
  DO_TESTS("floor( vector )", floor(a.array()), floor(a[i]) );
}

void
all_vector_test::test045() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(45)\n";
  DO_TESTS("log( 1+abs(vector) )", log(1.0+abs(a.array())), log(1.0+abs(a[i])) );
}

void
all_vector_test::test046() {
  preco();
  integer const sz{ MIN2(a,res) };
  cout << "\ntest(46)\n";
  DO_TESTS("log10( 1+abs(vector) )", log10(1.0+abs(a.array())), log10(1.0+abs(a[i])) );
}

void
all_vector_test::test047() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(47)\n";
  DO_TESTS("max(vector,vector)", a.array().max(b.array()), a[i] > b[i] ? a[i] : b[i] );
}

void
all_vector_test::test048() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(48)\n";
  DO_TESTS("min(vector,vector)", a.array().min(b.array()), a[i] < b[i] ? a[i] : b[i] );
}

void
all_vector_test::test049() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(49)\n";
  DO_TESTS("pow(1.0+abs(vector),1.0+abs(vector))",
           pow(1.0+abs(a.array()),1.0+abs(b.array())), pow(1.0+abs(a[i]),1.0+abs(b[i])) );
}

void
all_vector_test::test050() {
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(50)\n";
  DO_TESTS(
    "atan2(vector,vector)",
    a.binaryExpr(b, [] (double aa, double bb) { return std::atan2(aa,bb);} ),
    atan2(a[i],b[i])
  );
}

void
all_vector_test::test051() {
  preco();
  cout << "\ntest(51) norm1(vector)" << flush;
  Real const bfa{a.array().abs().sum()};
  Real bfb{0};
  for ( integer i{0}; i < a.size(); ++i) bfb += abs(a[i]);
  test_diff(abs(bfa-bfb));
}

void
all_vector_test::test052() {
  preco();
  cout << "\ntest(52) norm2(vector)" << flush;
  Real const bfa{ a.norm() };
  Real       bfb{ 0 };
  for ( integer i{0}; i < a.size(); ++i) bfb += power2(a[i]);
  test_diff(abs(bfa-sqrt(bfb)));
}

void
all_vector_test::test053() {
  preco();
  cout << "\ntest(53) normi(vector)" << flush;
  Real const bfa{ a.lpNorm<Eigen::Infinity>()};
  Real       bfb{ 0 };
  for ( integer i{0}; i < a.size(); ++i)
    if ( bfb < abs(a[i]) ) bfb = abs(a[i]);
  test_diff(abs(bfa-bfb));
}

void
all_vector_test::test054() {
  preco();
  cout << "\ntest(54) maxval(vector)" << flush;
  Real const bfa{ a.maxCoeff()};
  Real       bfb{ a[0] };
  for ( integer i{0}; i < a.size(); ++i)
    if ( bfb < a[i] ) bfb = a[i];
  test_diff(abs(bfa-bfb));
}

void
all_vector_test::test055() {
  preco();
  cout << "\ntest(55) minval(vector)" << flush;
  Real const bfa{ a.minCoeff() };
  Real       bfb{ a[0] };
  for ( integer i{0}; i < a.size(); ++i)
    if ( bfb > a[i] ) bfb = a[i];
  test_diff(abs(bfa-bfb));
}

void
all_vector_test::test056() {
  preco();
  cout << "\ntest(56) dot(vector,vector)" << flush;
  Real const bfa{ a.dot(b) };
  Real       bfb{ 0 };
  for ( integer i{0}; i < MIN2(a,b); ++i) bfb += a[i]*b[i];
  test_diff(abs(bfa-bfb));
}

void
all_vector_test::test057() { // specializations
  preco();
  cout << "\ntest(57) norm1(a) - sum(abs(a)) = " << flush;
  test_diff( a.lpNorm<1>() - abs(a.array()).sum() );
}

void
all_vector_test::test058() {
  preco();
  cout << "\ntest(58) dist(vector,vector)" << flush;
  Real const bfa{ (a-b).norm() };
  Real       bfb{ 0 };
  for ( integer i{0}; i < MIN2(a,b); ++i) bfb += power2(a[i]-b[i]);
  test_diff(abs(bfa-sqrt(bfb)));
}

void
all_vector_test::test059() {
  preco();
  cout << "\ntest(59) dist2(vector,vector)" << flush;
  Real const bfa{ (a-b).norm() };
  Real       bfb{ 0 };
  for ( integer i{0}; i < MIN2(a,b); ++i) bfb += power2(a[i]-b[i]);
  test_diff(abs(bfa-sqrt(bfb)));
}

void
all_vector_test::test060() {
  preco();
  cout << "\ntest(60) normp(vector,scalar)" << flush;
  Real const bfa{ pow(a.array().pow(scalar).sum(),1/scalar) };
  Real       bfb{ 0 };
  for ( integer i{0}; i < a.size(); ++i) {
    Real const bf{abs(a[i])};
    bfb += pow(bf,scalar);
  }
  test_diff( abs(bfa-pow(bfb,1/scalar)) );
}

void
all_vector_test::test061() { // specializations
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(61)\n";
  DO_TESTS("vector * scalar + vector",
	    b * scalar + c, b[i] * scalar + c[i]);
}

void
all_vector_test::test062() { // specializations
  preco();
  integer const sz{ MIN3(a,b,res) };
  cout << "\ntest(62)\n";
  DO_TESTS("vector * scalar + vector",
	    b * 2.0 + c, b[i] * 2.0 + c[i]);
}

void
all_vector_test::test063() {
  preco();
  cout << "\ntest(63) sum(1.0+abs(vector))" << flush;
  Real const bfa{ (1.0+abs(b.array())).sum() };
  Real       bfb{ 0 };
  for ( integer i{0}; i < MIN2(a,b); ++i) bfb += 1.0+abs(b[i]);
  test_diff(abs(bfa-bfb));
}

void
all_vector_test::test064() {
  preco();
  cout << "\ntest(64) prod(1.0+abs(vector))" << flush;
  Real const bfa{ (1.0+abs(b.array())).prod() };
  Real       bfb{ 1 };
  for ( integer i{0}; i < MIN2(a,b); ++i) bfb *= 1.0+abs(b[i]);
  test_diff(abs(bfa-bfb));
}

//********************************************************************
//********************************************************************
//********************************************************************
//********************************************************************

int
main() {
  try {
    all_vector_test allv(3);
    allv.do_all_tests();
  } catch ( exception const & exc ) {
    msg.error( exc.what() );
  } catch ( ... ) {
    msg.error("Errore Sconosciuto!\n");
  }
  return 0;
}
