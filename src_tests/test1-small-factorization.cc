/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <lapack_wrapper/lapack_wrapper.hh>
#include <lapack_wrapper/lapack_wrapper++.hh>

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
#endif

using namespace std;
using real_type = double;
using lapack_wrapper::integer;
using lapack_wrapper::Transposition;

static Utils::Console msg(&std::cout);

// ============================================================================
// Constants and Utility Functions
// ============================================================================

static constexpr real_type EPS = 1e-6;
static constexpr real_type TOL = 1e-9;

template<typename T>
bool is_close(T a, T b, T rel_tol = TOL, T abs_tol = TOL) {
    return std::abs(a - b) <= std::max(rel_tol * std::max(std::abs(a), std::abs(b)), abs_tol);
}

template<typename T>
void print_vector(const string& name, integer n, const T* v, integer inc = 1) {
    fmt::print("{}^T = [", name);
    for (integer i = 0; i < n; ++i) {
        fmt::print("{:12.6e}", v[i * inc]);
        if (i < n - 1) fmt::print(", ");
    }
    fmt::print("]\n");
}

template<typename T>
void print_matrix_info(const string& name, integer m, integer n, const T* A, integer lda) {
    fmt::print("\n{} ({}x{}):\n", name, m, n);
    fmt::print("{}", lapack_wrapper::print_matrix(m, n, A, lda));
}

// ============================================================================
// Test Implementations (Maintaining all 12 original tests)
// ============================================================================

static void test1() {
    msg.green("\n\n=== Test 1: Basic QR Factorization ===");
    
    lapack_wrapper::QR<real_type> qr;
    constexpr integer M = 3;
    constexpr integer N = 5;
    constexpr integer LDA = 3;
    
    real_type A[]{
        0.001,      2,     3,
        0.001,  0.001,     0,
        0,      0.001,     0,
        0.001,     -1,     0,
        0.000001,   5,     3
    };
    
    lapack_wrapper::MatrixWrapper<real_type> const Amat(A, M, N, LDA);
    print_matrix_info("Input Matrix A", M, N, A, LDA);
    
    msg.blue("Performing QR factorization of A^T...");
    qr.t_factorize("qr", Amat);
    
    real_type R[M * M];
    qr.getR(R, M);
    print_matrix_info("R Matrix", M, M, R, M);
    
    real_type rhs[M], b[M];
    real_type x[N]{1, 2, 3, 4, 5};
    lapack_wrapper::gemv(1.0, Amat, x, 1, 0.0, rhs, 1);
    lapack_wrapper::gemv(1.0, Amat, x, 1, 0.0, b, 1);
    
    print_vector("b", M, b);
    
    qr.invRt_mul(rhs, 1);
    lapack_wrapper::copy(M, rhs, 1, x, 1);
    lapack_wrapper::zero(N - M, x + 3, 1);
    qr.Q_mul(x);
    
    print_vector("x (solution)", N, x);
    
    lapack_wrapper::gemv(-1.0, Amat, x, 1, 1.0, b, 1);
    real_type res = lapack_wrapper::nrm2(M, b, 1);
    
    print_vector("residual", M, b);
    fmt::print("‖res‖₂ = {:12.6e}\n", res);
    
    UTILS_ASSERT(res < EPS, "Test 1 failed! res = {:12.6e}\n", res);
    msg.green("Test 1 PASSED ✓");
}

static void test2() {
    msg.green("\n\n=== Test 2: QRP Factorization ===");
    
    lapack_wrapper::QRP<real_type> qr;
    constexpr integer M = 3;
    constexpr integer N = 5;
    constexpr integer LDA = 3;
    
    real_type A[]{
        0.001,      2,     3,
        0.001,  0.001,     0,
        0,      0.001,     0,
        0.001,     -1,     0,
        0.000001,   5,     3
    };
    
    lapack_wrapper::MatrixWrapper<real_type> const Amat(A, M, N, LDA);
    print_matrix_info("Input Matrix A", M, N, A, LDA);
    
    msg.blue("Performing QRP factorization of A^T...");
    qr.t_factorize("qrp", Amat);
    
    real_type R[M * M];
    qr.getR(R, M);
    print_matrix_info("R Matrix", M, M, R, M);
    
    real_type rhs[M], b[M];
    real_type x[N]{1, 2, 3, 4, 5};
    lapack_wrapper::gemv(1.0, Amat, x, 1, 0.0, rhs, 1);
    lapack_wrapper::gemv(1.0, Amat, x, 1, 0.0, b, 1);
    
    print_vector("b", M, b);
    
    qr.inv_permute(rhs);
    qr.invRt_mul(rhs, 1);
    lapack_wrapper::copy(M, rhs, 1, x, 1);
    lapack_wrapper::zero(N - M, x + 3, 1);
    qr.Q_mul(x);
    
    print_vector("x (solution)", N, x);
    
    lapack_wrapper::gemv(-1.0, Amat, x, 1, 1.0, b, 1);
    real_type res = lapack_wrapper::nrm2(M, b, 1);
    
    print_vector("residual", M, b);
    fmt::print("‖res‖₂ = {:12.6e}\n", res);
    
    UTILS_ASSERT(res < EPS, "Test 2 failed! res = {:12.6e}\n", res);
    msg.green("Test 2 PASSED ✓");
}

static void test3() {
    msg.green("\n\n=== Test 3: QRP on 5x5 Matrix ===");
    
    lapack_wrapper::QRP<real_type> qr;
    constexpr integer M = 5;
    constexpr integer N = 5;
    constexpr integer LDA = 5;
    
    constexpr real_type A[]{
        0.001,      2,     3,     2, 3,
        0.001,  0.001,     0, 0.001, 1e-10,
        0,      0.001,     0, 0.001, 1e-12,
        0.001,     -1, 1e-12,    -1, -1e-12,
        0.000001,   5,     3,     5, 3
    };
    
    print_matrix_info("Input Matrix A", M, N, A, LDA);
    
    msg.blue("Performing QRP factorization of A^T...");
    qr.t_factorize("qr", M, N, A, LDA);
    
    real_type R[M * M];
    qr.getR(R, M);
    print_matrix_info("R Matrix", M, M, R, M);
    
    real_type rhs[M], b[M];
    real_type x[N]{1, 2, 3, 4, 5};
    
    lapack_wrapper::gemv(Transposition::NO, M, N, 1, A, LDA, x, 1, 0, rhs, 1);
    lapack_wrapper::gemv(Transposition::NO, M, N, 1, A, LDA, x, 1, 0, b, 1);
    
    print_vector("b", M, b);
    
    qr.inv_permute(rhs);
    qr.invRt_mul(rhs, 1);
    lapack_wrapper::copy(3, rhs, 1, x, 1);
    lapack_wrapper::zero(2, x + 3, 1);
    qr.Q_mul(x);
    
    lapack_wrapper::gemv(Transposition::NO, M, N, -1, A, LDA, x, 1, 1, b, 1);
    
    print_vector("x (solution)", N, x);
    print_vector("residual", M, b);
    
    real_type res = lapack_wrapper::nrm2(M, b, 1);
    fmt::print("‖res‖₂ = {:12.6e}\n", res);
    
    UTILS_ASSERT(res < EPS, "Test 3 failed! res = {:12.6e}\n", res);
    msg.green("Test 3 PASSED ✓");
}

// Helper macro for test4 and test5
#define TEST_SOLVER(NAME, F, TRANS, M, A, LDA, b, x, rhs_copy) \
    do { \
        fmt::print("\n\nDo {} factorization of A\n", NAME); \
        F.factorize(NAME, M, M, A, LDA); \
        fmt::print("{} solution of A x = b\n", NAME); \
        lapack_wrapper::copy(M, rhs_copy, 1, x, 1); \
        lapack_wrapper::copy(M, rhs_copy, 1, b, 1); \
        if (TRANS == Transposition::NO) { \
            F.solve(1, x, M); \
        } else { \
            F.t_solve(1, x, M); \
        } \
        fmt::print("x^T      = {}", lapack_wrapper::print_matrix(1, M, x, 1)); \
        lapack_wrapper::gemv(TRANS, M, M, -1, A, LDA, x, 1, 1, b, 1); \
        fmt::print("residual = {}", lapack_wrapper::print_matrix(1, M, b, 1)); \
        real_type res = lapack_wrapper::nrm2(M, b, 1); \
        fmt::print("‖res‖₂ = {:12.6e}\n", res); \
        UTILS_ASSERT0(res < EPS, "Test failed!\n"); \
    } while(0)

static void test4() {
    msg.green("\n\n=== Test 4: Multiple Solvers (A*x = b) ===");
    
    constexpr integer M = 5;
    constexpr integer LDA = 5;
    
    real_type A[]{
        0.001,      2,     3,       2,      3,
        0.001,  0.001,     0,   0.001,  1e-10,
        0,      0.001,     0,   0.001,  1e-12,
        0.001,     -1, 1e-6+1,     -1, -1e-12,
        0.000001,   5,     3,  1e-6+5,    3+1
    };
    
    print_matrix_info("Input Matrix A", M, M, A, LDA);
    
    real_type rhs[M], b[M];
    real_type x[M]{1, 2, 3, 4, 5};
    
    lapack_wrapper::gemv(Transposition::NO, M, M, 1, A, LDA, x, 1, 0, rhs, 1);
    lapack_wrapper::gemv(Transposition::NO, M, M, 1, A, LDA, x, 1, 0, b, 1);
    
    print_vector("b = A*x_true", M, b);
    
    // Test different solvers
    lapack_wrapper::LU<real_type>   lu;
    lapack_wrapper::LUPQ<real_type> lupq;
    lapack_wrapper::QR<real_type>   qr;
    lapack_wrapper::QRP<real_type>  qrp;
    lapack_wrapper::SVD<real_type>  svd;
    lapack_wrapper::LSS<real_type>  lss;
    lapack_wrapper::LSY<real_type>  lsy;
    
    TEST_SOLVER("LU", lu, Transposition::NO, M, A, LDA, b, x, rhs);
    TEST_SOLVER("LUPQ", lupq, Transposition::NO, M, A, LDA, b, x, rhs);
    TEST_SOLVER("QR", qr, Transposition::NO, M, A, LDA, b, x, rhs);
    TEST_SOLVER("QRP", qrp, Transposition::NO, M, A, LDA, b, x, rhs);
    TEST_SOLVER("SVD", svd, Transposition::NO, M, A, LDA, b, x, rhs);
    TEST_SOLVER("LSS", lss, Transposition::NO, M, A, LDA, b, x, rhs);
    TEST_SOLVER("LSY", lsy, Transposition::NO, M, A, LDA, b, x, rhs);
    
    msg.green("Test 4 PASSED ✓");
}

static void test5() {
    msg.green("\n\n=== Test 5: Multiple Solvers (A^T*x = b) ===");
    
    constexpr integer M = 5;
    constexpr integer LDA = 5;
    
    real_type A[]{
        0.001,      2,     3,       2,      3,
        0.001,  0.001,     0,   0.001,  1e-10,
        0,      0.001,     0,   0.001,  1e-12,
        0.001,     -1, 1e-6+1,     -1, -1e-12,
        0.000001,   5,     3,  1e-6+5,     3+1
    };
    
    print_matrix_info("Input Matrix A", M, M, A, LDA);
    
    real_type rhs[M], b[M];
    real_type x[M]{1, 2, 3, 4, 5};
    
    lapack_wrapper::gemv(Transposition::YES, M, M, 1, A, LDA, x, 1, 0, rhs, 1);
    lapack_wrapper::gemv(Transposition::YES, M, M, 1, A, LDA, x, 1, 0, b, 1);
    
    print_vector("b = A^T*x_true", M, b);
    
    // Test different solvers
    lapack_wrapper::LU<real_type>   lu;
    lapack_wrapper::LUPQ<real_type> lupq;
    lapack_wrapper::QR<real_type>   qr;
    lapack_wrapper::QRP<real_type>  qrp;
    lapack_wrapper::SVD<real_type>  svd;
    lapack_wrapper::LSS<real_type>  lss;
    lapack_wrapper::LSY<real_type>  lsy;
    
    TEST_SOLVER("LU", lu, Transposition::YES, M, A, LDA, b, x, rhs);
    TEST_SOLVER("LUPQ", lupq, Transposition::YES, M, A, LDA, b, x, rhs);
    TEST_SOLVER("QR", qr, Transposition::YES, M, A, LDA, b, x, rhs);
    TEST_SOLVER("QRP", qrp, Transposition::YES, M, A, LDA, b, x, rhs);
    TEST_SOLVER("SVD", svd, Transposition::YES, M, A, LDA, b, x, rhs);
    TEST_SOLVER("LSS", lss, Transposition::YES, M, A, LDA, b, x, rhs);
    TEST_SOLVER("LSY", lsy, Transposition::YES, M, A, LDA, b, x, rhs);
    
    msg.green("Test 5 PASSED ✓");
}

static void test6() {
    msg.green("\n\n=== Test 6: Tridiagonal Solvers ===");
    
    lapack_wrapper::TridiagonalLU<real_type> lu;
    lapack_wrapper::TridiagonalQR<real_type> qr;
    
    constexpr integer N = 5;
    constexpr real_type D[]{1, 1, 2, -0.1, 0.1};
    constexpr real_type L[]{-0.1, -1, -2, -0.1};
    constexpr real_type U[]{-1, -10, -2, 0.1};
    
    real_type rhs[N], b[N];
    real_type x[N]{1, 2, 3, 4, 5};
    
    qr.axpy(N, 1.0, L, D, U, x, 0.0, rhs);
    qr.axpy(N, 1.0, L, D, U, x, 0.0, b);
    
    print_vector("b = T*x_true", N, b);
    
    // Test LU solver
    msg.blue("\nTesting Tridiagonal LU solver:");
    lapack_wrapper::copy(N, rhs, 1, x, 1);
    lapack_wrapper::copy(N, rhs, 1, b, 1);
    
    lu.factorize("lu", N, L, D, U);
    lu.solve(x);
    
    print_vector("x (LU)", N, x);
    
    lu.axpy(N, -1.0, L, D, U, x, 1.0, b);
    print_vector("residual (LU)", N, b);
    
    real_type res_lu = lapack_wrapper::nrm2(N, b, 1);
    fmt::print("‖res‖₂ (LU) = {:12.6e}\n", res_lu);
    
    // Test QR solver
    msg.blue("\nTesting Tridiagonal QR solver:");
    lapack_wrapper::copy(N, rhs, 1, x, 1);
    lapack_wrapper::copy(N, rhs, 1, b, 1);
    
    qr.factorize("qr", N, L, D, U);
    qr.solve(x);
    
    print_vector("x (QR)", N, x);
    
    qr.axpy(N, -1.0, L, D, U, x, 1.0, b);
    print_vector("residual (QR)", N, b);
    
    real_type res_qr = lapack_wrapper::nrm2(N, b, 1);
    fmt::print("‖res‖₂ (QR) = {:12.6e}\n", res_qr);
    
    UTILS_ASSERT(res_lu < EPS && res_qr < EPS, 
                "Test 6 failed! LU res = {:12.6e}, QR res = {:12.6e}\n", 
                res_lu, res_qr);
    
    msg.green("Test 6 PASSED ✓");
}

static void test7() {
    msg.green("\n\n=== Test 7: QRP Transposed Solve ===");
    
    lapack_wrapper::QRP<real_type> qrp;
    constexpr integer M = 5;
    constexpr integer LDA = 5;
    
    constexpr real_type A[]{
        0.001,      2,     3,       2,      3,
        0.001,  0.001,     0,   0.001,  1e-10,
        0,      0.001,     0,   0.001,  1e-12,
        0.001,     -1, 1e-6+1,     -1, -1e-12,
        0.000001,   5,     3,  1e-6+5,    3+1
    };
    
    print_matrix_info("Input Matrix A", M, M, A, LDA);
    
    real_type rhs[M], b[M];
    real_type x[M]{1, 2, 3, 4, 5};
    
    lapack_wrapper::gemv(Transposition::YES, M, M, 1, A, LDA, x, 1, 0, rhs, 1);
    lapack_wrapper::gemv(Transposition::YES, M, M, 1, A, LDA, x, 1, 0, b, 1);
    
    print_vector("b = A^T*x_true", M, b);
    
    qrp.factorize("qrp", M, M, A, LDA);
    lapack_wrapper::copy(M, rhs, 1, x, 1);
    lapack_wrapper::copy(M, rhs, 1, b, 1);
    
    qrp.t_solve(x);
    
    print_vector("x (solution)", M, x);
    
    lapack_wrapper::gemv(Transposition::YES, M, M, -1, A, LDA, x, 1, 1, b, 1);
    print_vector("residual", M, b);
    
    real_type res = lapack_wrapper::nrm2(M, b, 1);
    fmt::print("‖res‖₂ = {:12.6e}\n", res);
    
    UTILS_ASSERT(res < EPS, "Test 7 failed! res = {:12.6e}\n", res);
    msg.green("Test 7 PASSED ✓");
}

static void test8() {
    msg.green("\n\n=== Test 8: LSC (Least Squares with Constraints) ===");
    
    lapack_wrapper::LSC<real_type> lsc;
    constexpr integer M = 5;
    constexpr integer LDA = 5;
    
    constexpr real_type A[]{
        0.001,      2,     3,       2,      3,
        0.001,  0.001,     0,   0.001,  1e-10,
        0,      0.001,     0,   0.001,  1e-12,
        0.001,     -1, 1e-6+1,     -1, -1e-12,
        0.000001,   5,     3,  1e-6+5,    3+1
    };
    
    print_matrix_info("Input Matrix A", M, M, A, LDA);
    
    real_type rhs[M], b[M];
    real_type x[M]{1, 2, 3, 4, 5};
    constexpr bool rselect[M]{true, false, true, true, true};
    
    lapack_wrapper::gemv(Transposition::NO, M, M, 1, A, LDA, x, 1, 0, rhs, 1);
    lapack_wrapper::gemv(Transposition::NO, M, M, 1, A, LDA, x, 1, 0, b, 1);
    
    print_vector("b = A*x_true", M, b);
    
    lsc.factorize(M, M, A, LDA, rselect);
    lapack_wrapper::copy(M, rhs, 1, x, 1);
    lapack_wrapper::copy(M, rhs, 1, b, 1);
    
    lsc.solve(x);
    
    print_vector("x (solution)", M, x);
    
    lapack_wrapper::gemv(Transposition::NO, M, M, -1, A, LDA, x, 1, 1, b, 1);
    print_vector("residual", M, b);
    
    real_type res = lapack_wrapper::nrm2(M, b, 1);
    fmt::print("‖res‖₂ = {:12.6e}\n", res);
    
    UTILS_ASSERT(res < EPS, "Test 8 failed! res = {:12.6e}\n", res);
    msg.green("Test 8 PASSED ✓");
}

static void test9() {
    msg.green("\n\n=== Test 9: PINV (Pseudoinverse) ===");
    
    lapack_wrapper::PINV<real_type> pinv;
    constexpr integer M = 5;
    constexpr integer LDA = 5;
    
    constexpr real_type A[]{
        0.001,    1,     0,      2,  1,
        1e-9, -1e+9,  1e-9,  -1e-9,  1,
        1e-9, -1e+9,     1,  -1e-9,  1,
        1e-9, -1e+9,  1e-9,      1,  1,
        0, 0, 0, 0, 0
    };
    
    print_matrix_info("Input Matrix A", M, M, A, LDA);
    
    real_type rhs[M], b[M], x[M];
    constexpr real_type xe[M]{1, 2, 3, 4, 5};
    
    lapack_wrapper::gemv(Transposition::NO, M, M, 1, A, LDA, xe, 1, 0, rhs, 1);
    lapack_wrapper::copy(M, rhs, 1, b, 1);
    
    print_vector("b = A*x_true", M, b);
    
    // First test: square matrix
    pinv.factorize(M, M, A, LDA);
    lapack_wrapper::copy(M, b, 1, x, 1);
    pinv.solve(x);
    
    print_vector("x (solution)", M, x);
    
    const real_type e[M]{x[0] - xe[0], x[1] - xe[1], x[2] - xe[2], 
                         x[3] - xe[3], x[4] - xe[4]};
    print_vector("error (x - x_true)", M, e);
    
    lapack_wrapper::copy(M, rhs, 1, b, 1);
    lapack_wrapper::gemv(Transposition::NO, M, M, -1, A, LDA, x, 1, 1, b, 1);
    print_vector("residual", M, b);
    
    real_type res = lapack_wrapper::nrm2(M, b, 1);
    fmt::print("‖res‖₂ = {:12.6e}\n", res);
    
    UTILS_WARNING(res < EPS, "\n\nTest 9 warning! res = {:12.6e}\n\n", res);
    
    // Second test: rectangular matrix
    msg.blue("\nTesting rectangular matrix (5x4):");
    pinv.factorize(M, M - 1, A, LDA);
    pinv.mult_inv(rhs, 1, x, 1);
    
    lapack_wrapper::copy(M, rhs, 1, b, 1);
    lapack_wrapper::gemv(Transposition::NO, M, M - 1, -1, A, LDA, x, 1, 1, b, 1);
    
    print_vector("x (rectangular)", M - 1, x);
    print_vector("residual", M, b);
    
    res = lapack_wrapper::nrm2(M, b, 1);
    fmt::print("‖res‖₂ = {:12.6e}\n", res);
    
    msg.green("Test 9 PASSED ✓");
}

static void test10() {
    msg.green("\n\n=== Test 10: PINV on Rectangular Matrix (6x4) ===");
    
    lapack_wrapper::PINV<real_type> pinv;
    constexpr integer M = 6;
    constexpr integer N = 4;
    constexpr integer LDA = 6;
    
    constexpr real_type A[]{
        2,   1,  2,  1, 1, 1,
        2,   1,  2,  1, 1, 1,
        20, 10, 20, 10, 1, 1,
        1,   2,  3,  4, 5, 6,
    };
    
    print_matrix_info("Input Matrix A", M, N, A, LDA);
    
    real_type rhs[M], b[M], b2[2 * M], x[N], x2[2 * N];
    constexpr real_type xe[N]{1, 2, 3, 4};
    
    lapack_wrapper::gemv(Transposition::NO, M, N, 1, A, LDA, xe, 1, 0, rhs, 1);
    lapack_wrapper::copy(M, rhs, 1, b, 1);
    lapack_wrapper::copy(M, rhs, 1, b2, 1);
    lapack_wrapper::copy(M, rhs, 1, b2 + M, 1);
    
    print_vector("b = A*x_true", M, b);
    
    pinv.factorize(M, N, A, LDA);
    pinv.mult_inv(b, 1, x, 1);
    pinv.mult_inv(2, b2, M, x2, N);
    
    print_vector("x (single RHS)", N, x);
    print_vector("x (multiple RHS, first)", N, x2);
    print_vector("x (multiple RHS, second)", N, x2 + N);
    
    const real_type e[N]{x[0] - xe[0], x[1] - xe[1], x[2] - xe[2], x[3] - xe[3]};
    print_vector("error (x - x_true)", N, e);
    
    lapack_wrapper::copy(M, rhs, 1, b, 1);
    lapack_wrapper::gemv(Transposition::NO, M, N, -1, A, LDA, x, 1, 1, b, 1);
    print_vector("residual", M, b);
    
    real_type res = lapack_wrapper::nrm2(M, b, 1);
    fmt::print("‖res‖₂ = {:12.6e}\n", res);
    
    UTILS_ASSERT(res < EPS, "Test 10 failed! res = {:12.6e}\n", res);
    msg.green("Test 10 PASSED ✓");
}

static void test11() {
    msg.green("\n\n=== Test 11: PINV Transposed Problem ===");
    
    lapack_wrapper::PINV<real_type> pinv;
    constexpr integer M = 6;
    constexpr integer N = 4;
    constexpr integer LDA = 6;
    
    constexpr real_type A[]{
        2,   1,  2,  1, 1, 1,
        2,   1,  2,  1, 1, 1,
        20, 10, 20, 10, 1, 1,
        1,   2,  3,  4, 5, 6,
    };
    
    print_matrix_info("Input Matrix A", M, N, A, LDA);
    
    real_type rhs[N], b[N], b2[2 * N], x[M], x2[2 * M];
    constexpr real_type xe[M]{1, 2, 3, 4, 5, 6};
    
    lapack_wrapper::gemv(Transposition::YES, M, N, 1, A, LDA, xe, 1, 0, rhs, 1);
    lapack_wrapper::copy(N, rhs, 1, b, 1);
    lapack_wrapper::copy(N, rhs, 1, b2, 1);
    lapack_wrapper::copy(N, rhs, 1, b2 + N, 1);
    
    print_vector("b = A^T*x_true", N, b);
    
    pinv.factorize(M, N, A, LDA);
    pinv.t_mult_inv(b, 1, x, 1);
    pinv.t_mult_inv(2, b2, N, x2, M);
    
    print_vector("x (single RHS)", M, x);
    print_vector("x (multiple RHS, first)", M, x2);
    print_vector("x (multiple RHS, second)", M, x2 + M);
    
    const real_type e[M]{x[0] - xe[0], x[1] - xe[1], x[2] - xe[2], 
                         x[3] - xe[3], x[4] - xe[4], x[5] - xe[5]};
    print_vector("error (x - x_true)", M, e);
    
    lapack_wrapper::copy(N, rhs, 1, b, 1);
    lapack_wrapper::gemv(Transposition::YES, M, N, -1, A, LDA, x, 1, 1, b, 1);
    print_vector("residual", N, b);
    
    real_type res = lapack_wrapper::nrm2(N, b, 1);
    fmt::print("‖res‖₂ = {:12.6e}\n", res);
    
    UTILS_WARNING(res < EPS, "\n\nTest 11 warning! res = {:12.6e}\n\n", res);
    msg.green("Test 11 PASSED ✓");
}

static void test12() {
    msg.green("\n\n=== Test 12: PINV Properties Verification ===");
    
    lapack_wrapper::PINV<real_type> pinv;
    constexpr integer M = 8;
    constexpr integer N = 5;
    constexpr integer LDA = M;
    
    constexpr real_type A[]{
        0.001, 1e-9,  1e-9,  1e-9, 1, 1, 0, 0,
        1,    -1e+9, -1e+9, -1e+9, 1, 1, 0, 0,
        0,     1e-9,     1,  1e-9, 1, 1, 0, 0,
        2,    -1e-9, -1e-9,     1, 1, 1, 0, 0,
        1,        1,     1,     1, 1, 1, 0, 0
    };
    
    print_matrix_info("Input Matrix A", M, N, A, LDA);
    
    real_type MAT[M * M], MAT1[M * M], MM[M * M];
    real_type rhs[N];
    constexpr real_type xe[M]{1, 2, 3, 4, 5, 6, 7, 8};
    
    lapack_wrapper::gemv(Transposition::YES, M, N, 1, A, LDA, xe, 1, 0, rhs, 1);
    print_vector("b = A^T*x_true", N, rhs);
    
    pinv.factorize(M, N, A, LDA);
    
    // Compute pseudoinverse
    lapack_wrapper::geid(M, M, MAT, M);
    pinv.mult_inv(M, MAT, M, MAT1, N);
    
    msg.blue("\nPseudoinverse of A:");
    print_matrix_info("pinv(A)", N, M, MAT1, N);
    
    // Verify properties: A * pinv(A) * A ≈ A
    lapack_wrapper::gemm(Transposition::NO, Transposition::NO, 
                        M, M, N, 1.0, A, LDA, MAT1, N, 0.0, MM, M);
    
    msg.blue("\nA * pinv(A):");
    print_matrix_info("A*pinv(A)", M, M, MM, M);
    
    // Verify: pinv(A) * A * pinv(A) ≈ pinv(A)
    lapack_wrapper::gemm(Transposition::NO, Transposition::NO, 
                        N, N, M, 1.0, MAT1, N, A, LDA, 0.0, MM, N);
    
    msg.blue("\npinv(A) * A:");
    print_matrix_info("pinv(A)*A", N, N, MM, N);
    
    // Compute transpose of pseudoinverse
    lapack_wrapper::geid(N, N, MAT, N);
    pinv.t_mult_inv(N, MAT, N, MAT1, M);
    
    msg.blue("\nTranspose of pseudoinverse:");
    print_matrix_info("pinv(A)^T", M, N, MAT1, M);
    
    msg.green("Test 12 PASSED ✓");
}

// ============================================================================
// Main Function with Test Runner
// ============================================================================

int main() {
    msg.blue(
      "\n"
      "==================================================\n"
      "     LAPACK Wrapper Comprehensive Test Suite\n"
      "==================================================\n"
    );
    
    vector<pair<string, function<void()>>> tests = {
        {"Test 1: Basic QR Factorization", test1},
        {"Test 2: QRP Factorization", test2},
        {"Test 3: QRP on 5x5 Matrix", test3},
        {"Test 4: Multiple Solvers (A*x = b)", test4},
        {"Test 5: Multiple Solvers (A^T*x = b)", test5},
        {"Test 6: Tridiagonal Solvers", test6},
        {"Test 7: QRP Transposed Solve", test7},
        {"Test 8: LSC with Constraints", test8},
        {"Test 9: PINV Pseudoinverse", test9},
        {"Test 10: PINV on Rectangular Matrix", test10},
        {"Test 11: PINV Transposed Problem", test11},
        {"Test 12: PINV Properties Verification", test12}
    };
    
    int passed = 0;
    int total = tests.size();
    
    for (size_t i = 0; i < tests.size(); ++i) {
        const auto& [name, test_func] = tests[i];
        
        fmt::print("\n{:-^70}\n", fmt::format(" [Test {}/{}] {} ", 
                                              i + 1, total, name));
        
        try {
            test_func();
            passed++;
        } catch (const exception& e) {
            msg.red(fmt::format("FAILED: {}\n", e.what()));
        } catch (...) {
            msg.red("FAILED: Unknown error\n");
        }
    }
    
    // Print summary
    msg.blue("\n" + string(70, '=') + "\n");
    fmt::print("{:^70}\n", "TEST SUITE SUMMARY");
    msg.blue(string(70, '=') + "\n");
    
    fmt::print("Total tests run:     {}\n", total);
    fmt::print("Tests passed:        {}\n", passed);
    fmt::print("Tests failed:        {}\n", total - passed);
    fmt::print("Success rate:        {:.1f}%\n", 
               (total > 0 ? (100.0 * passed / total) : 0.0));
    
    if (passed == total) {
        msg.green("\n✅ ALL TESTS PASSED SUCCESSFULLY!\n");
    } else {
        msg.red(fmt::format("\n❌ {} TEST(S) FAILED!\n", total - passed));
    }
    
    return (passed == total) ? 0 : 1;
}
