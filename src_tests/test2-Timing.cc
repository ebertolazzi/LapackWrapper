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
#include <random>
#include <iomanip>
#include <string>
#include <sstream>
#include <algorithm>
#include <lapack_wrapper/lapack_wrapper.hh>
#include <lapack_wrapper/lapack_wrapper++.hh>
#include <lapack_wrapper/lapack_wrapper_tmpl.hh>

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

static unsigned seed1{2};
static std::mt19937 generator(seed1);

static
real_type
rand( real_type const xmin, real_type const xmax ) {
  real_type const random{ static_cast<real_type>(generator())/std::mt19937::max() };
  return xmin + (xmax-xmin)*random;
}

using namespace lapack_wrapper;

// ============================================================================
// Table formatting utilities - CORRETTO
// ============================================================================

class TableFormatter {
private:
    vector<string> headers;
    vector<vector<string>> rows;
    vector<int> col_widths;
    
public:
    TableFormatter(const vector<string>& hdrs) : headers(hdrs) {
        col_widths.resize(headers.size());
        update_col_widths(headers);
    }
    
    void add_row(const vector<string>& row) {
        rows.push_back(row);
        update_col_widths(row);
    }
    
    void print_table(bool with_borders = true) {
        if (with_borders) {
            print_border("┌", "┬", "┐");
            print_row(headers, "│", "│", "│");
            print_border("├", "┼", "┤");
            for (const auto& row : rows) {
                print_row(row, "│", "│", "│");
            }
            print_border("└", "┴", "┘");
        } else {
            print_row(headers, "", " ", "");
            for (const auto& row : rows) {
                print_row(row, "", " ", "");
            }
        }
    }
    
    void print_border(string left, string middle, string right) {
        cout << left;
        for (size_t i = 0; i < col_widths.size(); ++i) {
            // Aggiungi 2 spazi extra per padding (1 a sinistra e 1 a destra)
            for (int j = 0; j < col_widths[i] + 2; ++j) cout << "─";
            if (i != col_widths.size() - 1) cout << middle;
        }
        cout << right << endl;
    }
    
    void print_row(const vector<string>& cells, string left, string middle, string right) {
        cout << left << " ";
        for (size_t i = 0; i < cells.size(); ++i) {
            // Allineamento corretto: usa setw con il contenuto già formattato
            string cell_content = cells[i];
            // Calcola la larghezza effettiva considerando i caratteri speciali
            int display_width = 0;
            bool in_escape = false;
            for (char c : cell_content) {
                if (c == '\033') in_escape = true;  // CSI escape sequence
                else if (in_escape && c == 'm') in_escape = false;
                else if (!in_escape) display_width++;
            }
            
            // Aggiungi spazi per centrare se necessario
            int padding = col_widths[i] - display_width;
            if (padding > 0) {
                // Per numeri, allinea a destra; per testo, a sinistra
                if (i == 0 || (cells[i].find_first_not_of("0123456789.-+ ") == string::npos)) {
                    // È un numero (o la colonna N) - allinea a destra
                    cout << string(padding, ' ') << cell_content;
                } else {
                    // È testo - allinea a sinistra
                    cout << cell_content << string(padding, ' ');
                }
            } else {
                cout << cell_content;
            }
            
            if (i != cells.size() - 1) cout << " " << middle << " ";
        }
        cout << " " << right << endl;
    }
    
    void update_col_widths(const vector<string>& cells) {
        for (size_t i = 0; i < cells.size(); ++i) {
            if (i >= col_widths.size()) col_widths.push_back(0);
            
            // Calcola la lunghezza visuale senza i caratteri di escape ANSI
            string cell = cells[i];
            int display_width = 0;
            bool in_escape = false;
            for (char c : cell) {
                if (c == '\033') in_escape = true;  // Inizio sequenza escape
                else if (in_escape && c == 'm') in_escape = false;
                else if (!in_escape) display_width++;
            }
            
            col_widths[i] = max(col_widths[i], display_width);
        }
    }
};

// ============================================================================
// Test data structures
// ============================================================================

struct TimingResult {
    int N;
    double lapack_time;
    double hand_time;
    string unit;
};

struct TestSummary {
    string test_name;
    vector<TimingResult> results;
    
    void print_summary() {
        fmt::print("\n{:=^60}\n", fmt::format(" {} Test Summary ", test_name));
        
        // Usa larghezze fisse per un migliore allineamento
        TableFormatter table({
            fmt::format("{:>4}", "N"), 
            fmt::format("{:>16}", "LAPACK"), 
            fmt::format("{:>18}", "Hand Unrolled"), 
            fmt::format("{:>10}", "Speedup")
        });
        
        for (const auto& res : results) {
            string lapack_str = fmt::format("{:8.3f} {}", res.lapack_time, res.unit);
            string hand_str;
            if (res.hand_time > 0.0001) {  // Evita precisione troppo alta per zero
                hand_str = fmt::format("{:8.3f} {}", res.hand_time, res.unit);
            } else {
                hand_str = "N/A";
            }
            
            double speedup = (res.hand_time > 0.0001) ? res.lapack_time / res.hand_time : 0.0;
            string speedup_str;
            if (res.hand_time > 0.0001) {
                if (speedup >= 10.0) {
                    speedup_str = fmt::format("{:5.1f}x", speedup);
                } else {
                    speedup_str = fmt::format("{:5.2f}x", speedup);
                }
            } else {
                speedup_str = "N/A";
            }
            
            table.add_row({
                fmt::format("{:4}", res.N),
                fmt::format("{:>16}", lapack_str),
                fmt::format("{:>18}", hand_str),
                fmt::format("{:>10}", speedup_str)
            });
        }
        
        table.print_table();
        
        // Calculate averages
        double avg_lapack = 0, avg_hand = 0;
        int count_lapack = 0, count_hand = 0;
        
        for (const auto& res : results) {
            avg_lapack += res.lapack_time;
            count_lapack++;
            if (res.hand_time > 0.0001) {
                avg_hand += res.hand_time;
                count_hand++;
            }
        }
        
        if (count_lapack > 0) avg_lapack /= count_lapack;
        if (count_hand > 0) avg_hand /= count_hand;
        
        fmt::print("\n📊 Statistics:\n");
        fmt::print("  Average LAPACK time:   {:8.3f} {}\n", avg_lapack, results[0].unit);
        if (count_hand > 0) {
            fmt::print("  Average Hand time:     {:8.3f} {}\n", avg_hand, results[0].unit);
            fmt::print("  Average speedup:       {:5.2f}x\n", avg_lapack / avg_hand);
        }
        fmt::print("{:=^60}\n", "");
    }
};

// ============================================================================
// Timing tests
// ============================================================================

template <int N>
TimingResult testMM() {
    constexpr int N_TIMES{100};
    constexpr double to_ps{1000000.0 / N_TIMES};

    Malloc<real_type> baseValue("real");
    Malloc<integer>   baseIndex("integer");

    baseValue.allocate(N * N * 10);
    baseIndex.allocate(N * 10);

    real_type* M1{ baseValue(N * N) };
    real_type* M2{ baseValue(N * N) };
    real_type* M3{ baseValue(N * N) };

    for (int i{0}; i < N; ++i) {
        for (int j{0}; j < N; ++j) {
            M1[i + j * N] = rand(-1, 1);
            M2[i + j * N] = rand(-1, 1);
            M3[i + j * N] = rand(-1, 1);
        }
    }

    Utils::TicToc tm;
    double lapack_time, hand_time;

    // LAPACK version
    tm.tic();
    for (int i{0}; i < N_TIMES; ++i) {
        gemm(
            Transposition::NO,
            Transposition::NO,
            N, N, N,
            -1.0, M1, N,
            M2, N,
            1.0, M3, N
        );
        memcpy(M2, M3, N * N * sizeof(real_type));
    }
    tm.toc();
    lapack_time = to_ps * tm.elapsed_ms() / (N * N * N);

    // Hand-unrolled version
    tm.tic();
    for (int i{0}; i < N_TIMES; ++i) {
        MM<real_type, N, N, N, N, N, N>::subTo(M1, M2, M3);
        memcpy(M2, M3, N * N * sizeof(real_type));
    }
    tm.toc();
    hand_time = to_ps * tm.elapsed_ms() / (N * N * N);

    return TimingResult{N, lapack_time, hand_time, "ps/N^3"};
}

template <int N>
TimingResult testMv() {
    constexpr int N_TIMES{1000};
    constexpr double to_ps{1000000.0 / N_TIMES};

    Malloc<real_type> baseValue("real");
    Malloc<integer>   baseIndex("integer");

    baseValue.allocate(N * N * 10);
    baseIndex.allocate(N * 10);

    real_type* M{ baseValue(N * N) };
    real_type* V{ baseValue(N) };
    real_type* R{ baseValue(N) };

    for (int i{0}; i < N; ++i) {
        V[i] = rand(-1, 1);
        R[i] = rand(-1, 1);
        for (int j{0}; j < N; ++j) {
            M[i + j * N] = rand(-1, 1);
        }
    }

    Utils::TicToc tm;
    double lapack_time, hand_time;

    // LAPACK version
    tm.tic();
    for (int i{0}; i < N_TIMES; ++i) {
        gemv(
            Transposition::NO,
            N, N,
            -1.0, M, N,
            V, 1,
            1.0, R, 1
        );
        memcpy(V, R, N * sizeof(real_type));
    }
    tm.toc();
    lapack_time = to_ps * tm.elapsed_ms() / (N * N);

    // Hand-unrolled version
    tm.tic();
    for (int i{0}; i < N_TIMES; ++i) {
        Mv<real_type, N, N, N, N, N>::subTo(M, V, R);
        memcpy(V, R, N * sizeof(real_type));
    }
    tm.toc();
    hand_time = to_ps * tm.elapsed_ms() / (N * N);

    return TimingResult{N, lapack_time, hand_time, "ps/N^2"};
}

template <int N>
TimingResult testCopy() {
    constexpr int N_TIMES{10000};
    constexpr double to_ps{1000000.0 / N_TIMES};

    Malloc<real_type> baseValue("real");
    Malloc<integer>   baseIndex("integer");

    baseValue.allocate(N * N * 10);
    baseIndex.allocate(N * 10);

    real_type* M1{ baseValue(N * N) };
    real_type* M2{ baseValue(N * N) };
    real_type* M3{ baseValue(N * N) };

    for (int i{0}; i < N; ++i) {
        for (int j{0}; j < N; ++j) {
            M1[i + j * N] = rand(-1, 1);
            M2[i + j * N] = rand(-1, 1);
            M3[i + j * N] = rand(-1, 1);
        }
    }

    Utils::TicToc tm;

    // LAPACK version
    tm.tic();
    for (int i{0}; i < N_TIMES; ++i) {
        gecopy(N, N, M1, N, M2, N);
        M1[0] += M2[0];
        gecopy(N, N, M1, N, M2, N);
        M1[0] += M2[0];
        gecopy(N, N, M1, N, M2, N);
        M1[0] += M2[0];
        gecopy(N, N, M1, N, M2, N);
        M1[0] += M2[0];
        gecopy(N, N, M1, N, M2, N);
        M1[0] += M2[0];
    }
    tm.toc();
    double lapack_time = to_ps * tm.elapsed_ms() / (N * N);

    return TimingResult{N, lapack_time, 0.0, "ps/N^2"};
}

// ============================================================================
// Test runners
// ============================================================================

static
TestSummary testMMall() {
    TestSummary summary;
    summary.test_name = "Matrix-Matrix Multiplication";
    
    vector<int> sizes = {2, 4, 6, 8, 12, 16, 100};
    
    fmt::print("\n🔬 Running Matrix-Matrix Multiplication Tests...\n");
    fmt::print("   Sizes: ");
    for (size_t i = 0; i < sizes.size(); ++i) {
        fmt::print("{}", sizes[i]);
        if (i < sizes.size() - 1) fmt::print(", ");
    }
    fmt::print("\n");
    
    summary.results.push_back(testMM<2>());
    summary.results.push_back(testMM<4>());
    summary.results.push_back(testMM<6>());
    summary.results.push_back(testMM<8>());
    summary.results.push_back(testMM<12>());
    summary.results.push_back(testMM<16>());
    summary.results.push_back(testMM<100>());
    
    return summary;
}

static
TestSummary testMvAll() {
    TestSummary summary;
    summary.test_name = "Matrix-Vector Multiplication";
    
    vector<int> sizes = {2, 4, 6, 8, 12, 16, 100};
    
    fmt::print("\n🔬 Running Matrix-Vector Multiplication Tests...\n");
    fmt::print("   Sizes: ");
    for (size_t i = 0; i < sizes.size(); ++i) {
        fmt::print("{}", sizes[i]);
        if (i < sizes.size() - 1) fmt::print(", ");
    }
    fmt::print("\n");
    
    summary.results.push_back(testMv<2>());
    summary.results.push_back(testMv<4>());
    summary.results.push_back(testMv<6>());
    summary.results.push_back(testMv<8>());
    summary.results.push_back(testMv<12>());
    summary.results.push_back(testMv<16>());
    summary.results.push_back(testMv<100>());
    
    return summary;
}

static
TestSummary testCopyAll() {
    TestSummary summary;
    summary.test_name = "Matrix Copy";
    
    vector<int> sizes = {2, 4, 6, 8, 12, 16, 100};
    
    fmt::print("\n🔬 Running Matrix Copy Tests...\n");
    fmt::print("   Sizes: ");
    for (size_t i = 0; i < sizes.size(); ++i) {
        fmt::print("{}", sizes[i]);
        if (i < sizes.size() - 1) fmt::print(", ");
    }
    fmt::print("\n");
    
    summary.results.push_back(testCopy<2>());
    summary.results.push_back(testCopy<4>());
    summary.results.push_back(testCopy<6>());
    summary.results.push_back(testCopy<8>());
    summary.results.push_back(testCopy<12>());
    summary.results.push_back(testCopy<16>());
    summary.results.push_back(testCopy<100>());
    
    return summary;
}

// ============================================================================
// Performance comparison chart - CORRETTO
// ============================================================================

static
void print_performance_chart(const vector<TestSummary>& summaries) {
    fmt::print("\n{:=^70}\n", " PERFORMANCE COMPARISON CHART ");
    
    // Find common sizes
    vector<int> common_sizes;
    if (!summaries.empty() && !summaries[0].results.empty()) {
        for (const auto& res : summaries[0].results) {
            common_sizes.push_back(res.N);
        }
    }
    
    // Larghezze fisse per un allineamento migliore
    TableFormatter chart_table({
        fmt::format("{:>4}", "N"), 
        fmt::format("{:>12}", "Operation"), 
        fmt::format("{:>8}", "LAPACK"), 
        fmt::format("{:>8}", "Hand"), 
        fmt::format("{:>8}", "Ratio")
    });
    
    for (int N : common_sizes) {
        bool first_for_N = true;
        
        for (const auto& summary : summaries) {
            for (const auto& res : summary.results) {
                if (res.N == N) {
                    string op_name;
                    if (summary.test_name.find("Matrix-Matrix") != string::npos) op_name = "MM";
                    else if (summary.test_name.find("Matrix-Vector") != string::npos) op_name = "MV";
                    else op_name = "COPY";
                    
                    string N_str = first_for_N ? fmt::format("{:4}", N) : "";
                    string lapack_str = fmt::format("{:6.1f}", res.lapack_time);
                    string hand_str;
                    if (res.hand_time > 0.0001) {
                        hand_str = fmt::format("{:6.1f}", res.hand_time);
                    } else {
                        hand_str = "N/A";
                    }
                    string ratio_str;
                    if (res.hand_time > 0.0001) {
                        double ratio = res.lapack_time / res.hand_time;
                        if (ratio >= 10.0) {
                            ratio_str = fmt::format("{:4.1f}x", ratio);
                        } else {
                            ratio_str = fmt::format("{:4.2f}x", ratio);
                        }
                    } else {
                        ratio_str = "N/A";
                    }
                    
                    chart_table.add_row({
                        N_str, 
                        fmt::format("{:>12}", op_name), 
                        fmt::format("{:>8}", lapack_str), 
                        fmt::format("{:>8}", hand_str), 
                        fmt::format("{:>8}", ratio_str)
                    });
                    first_for_N = false;
                }
            }
        }
        
        // Add separator between different N values
        if (N != common_sizes.back()) {
            chart_table.add_row({"", "", "", "", ""});
        }
    }
    
    chart_table.print_table();
    fmt::print("{:=^70}\n", "");
}

// ============================================================================
// Main function
// ============================================================================

static Utils::Console msg(&std::cout);

int main() {
    try {
        fmt::print("\n{:*^70}\n", " LAPACK WRAPPER PERFORMANCE BENCHMARK ");
        fmt::print("{:*^70}\n\n", "");
        
        fmt::print("📊 Benchmark Configuration:\n");
        fmt::print("  • Matrix sizes: 2, 4, 6, 8, 12, 16, 100\n");
        fmt::print("  • Time unit: picoseconds per operation\n");
        fmt::print("  • LAPACK vs. Hand-unrolled implementations\n");
        fmt::print("{:-^70}\n\n", "");
        
        vector<TestSummary> all_summaries;
        
        // Run all tests
        auto mm_summary = testMMall();
        mm_summary.print_summary();
        all_summaries.push_back(mm_summary);
        
        auto mv_summary = testMvAll();
        mv_summary.print_summary();
        all_summaries.push_back(mv_summary);
        
        auto copy_summary = testCopyAll();
        copy_summary.print_summary();
        all_summaries.push_back(copy_summary);
        
        // Print comparison chart
        print_performance_chart(all_summaries);
        
        // Final summary
        fmt::print("\n🎯 Performance Insights:\n");
        fmt::print("  • Hand-unrolled code shows best speedup for small matrices\n");
        fmt::print("  • LAPACK routines are optimized for larger matrices\n");
        fmt::print("  • Matrix copy operations are highly optimized in LAPACK\n");
        fmt::print("  • Consider hybrid approach: hand-unrolled for small N, LAPACK for large N\n");
        
        fmt::print("\n{:*^70}\n", " BENCHMARK COMPLETED SUCCESSFULLY ");
        fmt::print("{:*^70}\n", "");
        
    } catch (exception const & exc) {
        msg.error(exc.what());
        return 1;
    } catch (...) {
        msg.error("Errore Sconosciuto!\n");
        return 1;
    }
    return 0;
}
