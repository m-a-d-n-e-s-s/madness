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

static void ran_fill_d(size_t n, double* a) {
    for (size_t i = 0; i < n; ++i) a[i] = ran();
}

static void ran_fill_f(size_t n, float* a) {
    for (size_t i = 0; i < n; ++i) a[i] = static_cast<float>(ran());
}

static void ran_fill_cd(size_t n, std::complex<double>* a) {
    for (size_t i = 0; i < n; ++i) a[i] = std::complex<double>(ran(), ran());
}

static void ran_fill_cf(size_t n, std::complex<float>* a) {
    for (size_t i = 0; i < n; ++i) a[i] = std::complex<float>(static_cast<float>(ran()), static_cast<float>(ran()));
}

template <typename T>
static void test_shape_real(long ni, long nj, long nk, long ldb,
                            const T* a, const T* b, T* c_test, T* c_ref,
                            long& num_tested, long& num_failed, const char* tname, double tol) {
    if (ldb == -1) ldb = nj;
    long n_c = (ni > 0 && nj > 0) ? ni * nj : 1;
    for (long i = 0; i < n_c; ++i) {
        c_test[i] = static_cast<T>(-999.0);
        c_ref[i] = static_cast<T>(-999.0);
    }

    mTxmq_reference(ni, nj, nk, c_ref, a, b, ldb);
    mTxmq(ni, nj, nk, c_test, a, b, ldb);

    num_tested++;
    if (ni <= 0 || nj <= 0) return;

    for (long i = 0; i < ni * nj; ++i) {
        double diff = std::abs(c_test[i] - c_ref[i]);
        double ref_val = std::abs(c_ref[i]);
        double cur_tol = tol + tol * ref_val;
        if (diff > cur_tol) {
            std::cerr << "FAILED test_mtxmq<" << tname << "> shape (" << ni << ", " << nj << ", " << nk
                      << ", ldb=" << ldb << ") at elem " << i
                      << ": test=" << c_test[i] << " ref=" << c_ref[i]
                      << " diff=" << diff << " tol=" << cur_tol << "\n";
            num_failed++;
            return;
        }
    }
}

template <typename cT, typename rT>
static void test_shape_mixed(long ni, long nj, long nk, long ldb,
                             const cT* a, const rT* b, cT* c_test, cT* c_ref,
                             long& num_tested, long& num_failed, const char* tname, double tol) {
    if (ldb == -1) ldb = nj;
    long n_c = (ni > 0 && nj > 0) ? ni * nj : 1;
    for (long i = 0; i < n_c; ++i) {
        c_test[i] = cT(-999.0, -999.0);
        c_ref[i] = cT(-999.0, -999.0);
    }

    mTxmq_reference(ni, nj, nk, c_ref, a, b, ldb);
    mTxmq(ni, nj, nk, c_test, a, b, ldb);

    num_tested++;
    if (ni <= 0 || nj <= 0) return;

    for (long i = 0; i < ni * nj; ++i) {
        double diff = std::abs(c_test[i] - c_ref[i]);
        double ref_val = std::abs(c_ref[i]);
        double cur_tol = tol + tol * ref_val;
        if (diff > cur_tol) {
            std::cerr << "FAILED test_mtxmq<" << tname << "> shape (" << ni << ", " << nj << ", " << nk
                      << ", ldb=" << ldb << ") at elem " << i
                      << ": test=" << c_test[i] << " ref=" << c_ref[i]
                      << " diff=" << diff << " tol=" << cur_tol << "\n";
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

    double *a_d = nullptr, *b_d = nullptr, *c_test_d = nullptr, *c_ref_d = nullptr;
    float *a_f = nullptr, *b_f = nullptr, *c_test_f = nullptr, *c_ref_f = nullptr;
    std::complex<double> *a_cd = nullptr, *c_test_cd = nullptr, *c_ref_cd = nullptr;
    std::complex<float> *a_cf = nullptr, *c_test_cf = nullptr, *c_ref_cf = nullptr;

    if (posix_memalign((void**)&a_d, 64, nimax * nkmax * sizeof(double)) != 0) return 1;
    if (posix_memalign((void**)&b_d, 64, nkmax * ldbmax * sizeof(double)) != 0) return 1;
    if (posix_memalign((void**)&c_test_d, 64, nimax * njmax * sizeof(double)) != 0) return 1;
    if (posix_memalign((void**)&c_ref_d, 64, nimax * njmax * sizeof(double)) != 0) return 1;

    if (posix_memalign((void**)&a_f, 64, nimax * nkmax * sizeof(float)) != 0) return 1;
    if (posix_memalign((void**)&b_f, 64, nkmax * ldbmax * sizeof(float)) != 0) return 1;
    if (posix_memalign((void**)&c_test_f, 64, nimax * njmax * sizeof(float)) != 0) return 1;
    if (posix_memalign((void**)&c_ref_f, 64, nimax * njmax * sizeof(float)) != 0) return 1;

    if (posix_memalign((void**)&a_cd, 64, nimax * nkmax * sizeof(std::complex<double>)) != 0) return 1;
    if (posix_memalign((void**)&c_test_cd, 64, nimax * njmax * sizeof(std::complex<double>)) != 0) return 1;
    if (posix_memalign((void**)&c_ref_cd, 64, nimax * njmax * sizeof(std::complex<double>)) != 0) return 1;

    if (posix_memalign((void**)&a_cf, 64, nimax * nkmax * sizeof(std::complex<float>)) != 0) return 1;
    if (posix_memalign((void**)&c_test_cf, 64, nimax * njmax * sizeof(std::complex<float>)) != 0) return 1;
    if (posix_memalign((void**)&c_ref_cf, 64, nimax * njmax * sizeof(std::complex<float>)) != 0) return 1;

    ran_fill_d(nimax * nkmax, a_d);
    ran_fill_d(nkmax * ldbmax, b_d);
    ran_fill_f(nimax * nkmax, a_f);
    ran_fill_f(nkmax * ldbmax, b_f);
    ran_fill_cd(nimax * nkmax, a_cd);
    ran_fill_cf(nimax * nkmax, a_cf);

    long num_tested = 0;
    long num_failed = 0;

    std::cout << "Testing mTxmq<double>...\n";
    for (long ni = 1; ni <= 24; ++ni) {
        for (long nj = 1; nj <= 12; ++nj) {
            for (long nk = 1; nk <= 12; ++nk) {
                test_shape_real(ni, nj, nk, nj, a_d, b_d, c_test_d, c_ref_d, num_tested, num_failed, "double", 1e-12);
            }
        }
    }

    std::vector<long> nis = {1, 2, 4, 8, 16, 32, 64, 100, 128, 216, 256, 400, 512};
    std::vector<long> njs = {1, 2, 4, 6, 8, 10, 12, 16, 20, 24, 32};
    std::vector<long> nks = {1, 2, 4, 6, 8, 10, 12, 16, 20, 24, 32};

    for (long ni : nis) {
        for (long nj : njs) {
            for (long nk : nks) {
                test_shape_real(ni, nj, nk, nj, a_d, b_d, c_test_d, c_ref_d, num_tested, num_failed, "double", 1e-12);
                if (nk >= nj) test_shape_real(ni, nj, nk, nk, a_d, b_d, c_test_d, c_ref_d, num_tested, num_failed, "double", 1e-12);
            }
        }
    }

    test_shape_real(0, 8, 8, 8, a_d, b_d, c_test_d, c_ref_d, num_tested, num_failed, "double", 1e-12);
    test_shape_real(8, 0, 8, 8, a_d, b_d, c_test_d, c_ref_d, num_tested, num_failed, "double", 1e-12);
    test_shape_real(8, 8, 0, 8, a_d, b_d, c_test_d, c_ref_d, num_tested, num_failed, "double", 1e-12);

    std::cout << "Testing mTxmq<float>...\n";
    for (long ni : {1, 4, 8, 16, 32, 64, 128, 216}) {
        for (long nj : {1, 4, 6, 8, 12, 16}) {
            for (long nk : {1, 4, 6, 8, 12, 16}) {
                test_shape_real(ni, nj, nk, nj, a_f, b_f, c_test_f, c_ref_f, num_tested, num_failed, "float", 1e-5);
                if (nk >= nj) test_shape_real(ni, nj, nk, nk, a_f, b_f, c_test_f, c_ref_f, num_tested, num_failed, "float", 1e-5);
            }
        }
    }

    std::cout << "Testing mTxmq(complex<double>*, complex<double>*, double*)...\n";
    for (long ni : {1, 4, 8, 16, 32, 64, 128}) {
        for (long nj : {1, 4, 6, 8, 12}) {
            for (long nk : {1, 4, 6, 8, 12}) {
                test_shape_mixed(ni, nj, nk, nj, a_cd, b_d, c_test_cd, c_ref_cd, num_tested, num_failed, "complex<double>*double", 1e-12);
                if (nk >= nj) test_shape_mixed(ni, nj, nk, nk, a_cd, b_d, c_test_cd, c_ref_cd, num_tested, num_failed, "complex<double>*double", 1e-12);
            }
        }
    }

    std::cout << "Testing mTxmq(complex<float>*, complex<float>*, float*)...\n";
    for (long ni : {1, 4, 8, 16, 32, 64}) {
        for (long nj : {1, 4, 6, 8}) {
            for (long nk : {1, 4, 6, 8}) {
                test_shape_mixed(ni, nj, nk, nj, a_cf, b_f, c_test_cf, c_ref_cf, num_tested, num_failed, "complex<float>*float", 1e-5);
            }
        }
    }

    std::cout << "Tested " << num_tested << " shapes across double, float, and mixed types: "
              << (num_failed == 0 ? "ALL PASSED" : "FAILURES DETECTED") << "!\n";

    if (run_bench) {
        run_benchmarks(nimax, njmax, nkmax, a_d, b_d, c_test_d);
    }

    free(a_d); free(b_d); free(c_test_d); free(c_ref_d);
    free(a_f); free(b_f); free(c_test_f); free(c_ref_f);
    free(a_cd); free(c_test_cd); free(c_ref_cd);
    free(a_cf); free(c_test_cf); free(c_ref_cf);

    SafeMPI::Finalize();
    return num_failed == 0 ? 0 : 1;
}
