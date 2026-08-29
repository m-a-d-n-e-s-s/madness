/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory
  Copyright (C) 2026 MADNESS developers

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

#include <madness/madness_config.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <complex>

#include <madness/world/safempi.h>
#include <madness/world/posixmem.h>
#include <madness/tensor/cblas.h>
#include <madness/tensor/tensor.h>
#include <madness/tensor/mxm.h>
#include <madness/tensor/mTxmq.h>

using namespace madness;

static double ran() {
    static unsigned long seed = 76521;
    seed = seed * 1812433253 + 12345;
    return ((double)(seed & 0x7fffffff)) * 4.6566128752458e-10;
}

static void ran_fill(size_t n, double* a) {
    for (size_t i = 0; i < n; ++i) a[i] = ran();
}

static void test_shape(long ni, long nj, long nk, long ldb,
                       const double* a, const double* b, double* c_test, double* c_ref,
                       long& num_tested, long& num_failed) {
    if (ldb == -1) ldb = nj;
    long n_c = (ni > 0 && nj > 0) ? ni * nj : 1;
    for (long i = 0; i < n_c; ++i) {
        c_test[i] = -999.0;
        c_ref[i] = -999.0;
    }

    mTxmq_reference(ni, nj, nk, c_ref, a, b, ldb);
    mTxmq(ni, nj, nk, c_test, a, b, ldb);

    num_tested++;
    if (ni <= 0 || nj <= 0) return;

    for (long i = 0; i < ni * nj; ++i) {
        double diff = std::abs(c_test[i] - c_ref[i]);
        double ref_val = std::abs(c_ref[i]);
        double tol = 1e-12 + 1e-12 * ref_val;
        if (diff > tol) {
            std::cerr << "FAILED test_mtxmq shape (" << ni << ", " << nj << ", " << nk
                      << ", ldb=" << ldb << ") at elem " << i
                      << ": test=" << c_test[i] << " ref=" << c_ref[i]
                      << " diff=" << diff << " tol=" << tol << "\n";
            num_failed++;
            return;
        }
    }
}

static void run_benchmarks(long nimax, long njmax, long nkmax,
                           const double* a, const double* b, double* c) {
    std::cout << "\n============================== mTxmq Benchmarks ==============================\n";
    std::cout << "               Shape (ixjxk)      mTxmq (GF/s)      BLAS (GF/s)      Speedup\n";
    std::cout << "------------------------------------------------------------------------------\n";

    struct BenchShape { const char* label; long ni, nj, nk; };
    std::vector<BenchShape> shapes = {
        {"1D Tensor (k=6)",    6, 6, 6},
        {"1D Tensor (k=8)",    8, 8, 8},
        {"1D Tensor (k=10)",  10, 10, 10},
        {"2D Tensor (k=6)",   36, 6, 6},
        {"2D Tensor (k=8)",   64, 8, 8},
        {"2D Tensor (k=10)", 100, 10, 10},
        {"3D Tensor (k=6)",  216, 6, 6},
        {"3D Tensor (k=8)",  512, 8, 8},
        {"Separated Op (r=4)", 64, 4, 8},
        {"Separated Op (r=8)", 64, 8, 8},
        {"Separated Op (r=12)", 64, 12, 8},
    };

    for (const auto& s : shapes) {
        if (s.ni * s.nk > nimax * nkmax || s.nk * s.nj > nkmax * njmax || s.ni * s.nj > nimax * njmax)
            continue;

        double nflop = 2.0 * s.ni * s.nj * s.nk;
        int nloops = std::max(50, int(2e6 / (nflop + 1)));

        // Warmup
        for (int l = 0; l < 10; ++l) mTxmq(s.ni, s.nj, s.nk, c, a, b);

        // Benchmark mTxmq
        double t_fast = 1e9;
        for (int rep = 0; rep < 5; ++rep) {
            double t0 = SafeMPI::Wtime();
            for (int l = 0; l < nloops; ++l) {
                mTxmq(s.ni, s.nj, s.nk, c, a, b);
            }
            double elapsed = SafeMPI::Wtime() - t0;
            if (elapsed < t_fast) t_fast = elapsed;
        }
        double gflops_fast = 1e-9 * nflop * nloops / t_fast;

        // Benchmark BLAS
        double t_blas = 1e9;
        for (int rep = 0; rep < 5; ++rep) {
            double t0 = SafeMPI::Wtime();
            for (int l = 0; l < nloops; ++l) {
                cblas::gemm(cblas::NoTrans, cblas::Trans, s.nj, s.ni, s.nk, 1.0, b, s.nj, a, s.ni, 0.0, c, s.nj);
            }
            double elapsed = SafeMPI::Wtime() - t0;
            if (elapsed < t_blas) t_blas = elapsed;
        }
        double gflops_blas = 1e-9 * nflop * nloops / t_blas;

        char shape_str[32];
        snprintf(shape_str, sizeof(shape_str), "(%ld x %ld x %ld)", s.ni, s.nj, s.nk);
        printf("%-20s %-16s %10.2f        %10.2f        %6.2fx\n",
               s.label, shape_str, gflops_fast, gflops_blas, gflops_fast / (gflops_blas > 0 ? gflops_blas : 1.0));
    }
    std::cout << "==============================================================================\n\n";
}

int main(int argc, char* argv[]) {
    SafeMPI::Init_thread(argc, argv, MPI_THREAD_SINGLE);

    bool run_bench = false;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--benchmark") == 0 || strcmp(argv[i], "-b") == 0) {
            run_bench = true;
        }
    }

    const long nimax = 600;
    const long njmax = 64;
    const long nkmax = 64;
    const long ldbmax = 128;

    double *a = nullptr, *b = nullptr, *c_test = nullptr, *c_ref = nullptr;
    if (posix_memalign((void**)&a, 64, nimax * nkmax * sizeof(double)) != 0) return 1;
    if (posix_memalign((void**)&b, 64, nkmax * ldbmax * sizeof(double)) != 0) return 1;
    if (posix_memalign((void**)&c_test, 64, nimax * njmax * sizeof(double)) != 0) return 1;
    if (posix_memalign((void**)&c_ref, 64, nimax * njmax * sizeof(double)) != 0) return 1;

    ran_fill(nimax * nkmax, a);
    ran_fill(nkmax * ldbmax, b);

    long num_tested = 0;
    long num_failed = 0;

    std::cout << "Testing mTxmq (double precision)...\n";

    // 1. Full dense sweep of small matrices: 1 <= ni <= 24, 1 <= nj <= 12, 1 <= nk <= 12
    for (long ni = 1; ni <= 24; ++ni) {
        for (long nj = 1; nj <= 12; ++nj) {
            for (long nk = 1; nk <= 12; ++nk) {
                test_shape(ni, nj, nk, nj, a, b, c_test, c_ref, num_tested, num_failed);
            }
        }
    }

    // 2. Strided sweep with ldb == nk: 1 <= ni <= 24, 1 <= nj <= 12, 1 <= nk <= 12
    for (long ni = 1; ni <= 24; ++ni) {
        for (long nj = 1; nj <= 12; ++nj) {
            for (long nk = 1; nk <= 12; ++nk) {
                if (nk >= nj) {
                    test_shape(ni, nj, nk, nk, a, b, c_test, c_ref, num_tested, num_failed);
                }
            }
        }
    }

    // 3. Representative MADNESS tensor shapes
    std::vector<long> nis = {1, 2, 4, 8, 16, 32, 64, 100, 128, 216, 256, 400, 512};
    std::vector<long> njs = {1, 2, 4, 6, 8, 10, 12, 16, 20, 24, 30};
    std::vector<long> nks = {1, 2, 4, 6, 8, 10, 12, 16, 20, 24, 30};

    for (long ni : nis) {
        for (long nj : njs) {
            for (long nk : nks) {
                // Dense
                test_shape(ni, nj, nk, nj, a, b, c_test, c_ref, num_tested, num_failed);
                // Strided (ldb = nk or arbitrary)
                if (nk >= nj) test_shape(ni, nj, nk, nk, a, b, c_test, c_ref, num_tested, num_failed);
                test_shape(ni, nj, nk, nj + 4, a, b, c_test, c_ref, num_tested, num_failed);
            }
        }
    }

    // 4. Edge cases
    test_shape(0, 8, 8, 8, a, b, c_test, c_ref, num_tested, num_failed);
    test_shape(8, 0, 8, 8, a, b, c_test, c_ref, num_tested, num_failed);
    test_shape(8, 8, 0, 8, a, b, c_test, c_ref, num_tested, num_failed);

    std::cout << "Tested " << num_tested << " shapes: "
              << (num_failed == 0 ? "ALL PASSED" : "FAILURES DETECTED") << "!\n";

    if (run_bench) {
        run_benchmarks(nimax, njmax, nkmax, a, b, c_test);
    }

    free(a);
    free(b);
    free(c_test);
    free(c_ref);

    SafeMPI::Finalize();
    return num_failed == 0 ? 0 : 1;
}
