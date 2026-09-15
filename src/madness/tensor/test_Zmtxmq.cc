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

typedef std::complex<double> double_complex;
typedef std::complex<float> float_complex;

static double ran() {
    static unsigned long seed = 76521;
    seed = seed * 1812433253 + 12345;
    return ((double)(seed & 0x7fffffff)) * 4.6566128752458e-10;
}

static double_complex ran_cd() {
    return double_complex(ran(), ran());
}

static float_complex ran_cf() {
    return float_complex(static_cast<float>(ran()), static_cast<float>(ran()));
}

template <typename T>
static void test_shape_c(long ni, long nj, long nk, long ldb,
                         const T* a, const T* b,
                         T* c_test, T* c_ref,
                         long& num_tested, long& num_failed, const char* tname, double tol) {
    if (ldb == -1) ldb = nj;
    long n_c = (ni > 0 && nj > 0) ? ni * nj : 1;
    for (long i = 0; i < n_c; ++i) {
        c_test[i] = T(-999.0, -999.0);
        c_ref[i] = T(-999.0, -999.0);
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
            std::cerr << "FAILED test_Zmtxmq<" << tname << "> shape (" << ni << ", " << nj << ", " << nk
                      << ", ldb=" << ldb << ") at elem " << i
                      << ": test=" << c_test[i] << " ref=" << c_ref[i]
                      << " diff=" << diff << " tol=" << cur_tol << "\n";
            num_failed++;
            return;
        }
    }
}

int main(int argc, char* argv[]) {
    SafeMPI::Init_thread(argc, argv, MPI_THREAD_SINGLE);

    const long nimax = 600;
    const long njmax = 64;
    const long nkmax = 64;
    const long ldbmax = 128;

    double_complex *a_cd = nullptr, *b_cd = nullptr, *c_test_cd = nullptr, *c_ref_cd = nullptr;
    float_complex *a_cf = nullptr, *b_cf = nullptr, *c_test_cf = nullptr, *c_ref_cf = nullptr;

    if (posix_memalign((void**)&a_cd, 64, nimax * nkmax * sizeof(double_complex)) != 0) return 1;
    if (posix_memalign((void**)&b_cd, 64, nkmax * ldbmax * sizeof(double_complex)) != 0) return 1;
    if (posix_memalign((void**)&c_test_cd, 64, nimax * njmax * sizeof(double_complex)) != 0) return 1;
    if (posix_memalign((void**)&c_ref_cd, 64, nimax * njmax * sizeof(double_complex)) != 0) return 1;

    if (posix_memalign((void**)&a_cf, 64, nimax * nkmax * sizeof(float_complex)) != 0) return 1;
    if (posix_memalign((void**)&b_cf, 64, nkmax * ldbmax * sizeof(float_complex)) != 0) return 1;
    if (posix_memalign((void**)&c_test_cf, 64, nimax * njmax * sizeof(float_complex)) != 0) return 1;
    if (posix_memalign((void**)&c_ref_cf, 64, nimax * njmax * sizeof(float_complex)) != 0) return 1;

    for (size_t i = 0; i < nimax * nkmax; ++i) a_cd[i] = ran_cd();
    for (size_t i = 0; i < nkmax * ldbmax; ++i) b_cd[i] = ran_cd();
    for (size_t i = 0; i < nimax * nkmax; ++i) a_cf[i] = ran_cf();
    for (size_t i = 0; i < nkmax * ldbmax; ++i) b_cf[i] = ran_cf();

    long num_tested = 0;
    long num_failed = 0;

    std::cout << "Testing mTxmq<complex<double>>...\n";
    for (long ni = 1; ni <= 24; ++ni) {
        for (long nj = 1; nj <= 12; ++nj) {
            for (long nk = 1; nk <= 12; ++nk) {
                test_shape_c(ni, nj, nk, nj, a_cd, b_cd, c_test_cd, c_ref_cd, num_tested, num_failed, "complex<double>", 1e-12);
            }
        }
    }

    std::vector<long> nis = {1, 2, 4, 8, 16, 32, 64, 100, 128, 216, 256, 400};
    std::vector<long> njs = {1, 2, 4, 6, 8, 10, 12, 16, 20, 24};
    std::vector<long> nks = {1, 2, 4, 6, 8, 10, 12, 16, 20, 24};

    for (long ni : nis) {
        for (long nj : njs) {
            for (long nk : nks) {
                test_shape_c(ni, nj, nk, nj, a_cd, b_cd, c_test_cd, c_ref_cd, num_tested, num_failed, "complex<double>", 1e-12);
                if (nk >= nj) test_shape_c(ni, nj, nk, nk, a_cd, b_cd, c_test_cd, c_ref_cd, num_tested, num_failed, "complex<double>", 1e-12);
            }
        }
    }

    test_shape_c(0, 8, 8, 8, a_cd, b_cd, c_test_cd, c_ref_cd, num_tested, num_failed, "complex<double>", 1e-12);
    test_shape_c(8, 0, 8, 8, a_cd, b_cd, c_test_cd, c_ref_cd, num_tested, num_failed, "complex<double>", 1e-12);
    test_shape_c(8, 8, 0, 8, a_cd, b_cd, c_test_cd, c_ref_cd, num_tested, num_failed, "complex<double>", 1e-12);

    std::cout << "Testing mTxmq<complex<float>>...\n";
    for (long ni : {1, 4, 8, 16, 32, 64, 128}) {
        for (long nj : {1, 4, 6, 8, 12}) {
            for (long nk : {1, 4, 6, 8, 12}) {
                test_shape_c(ni, nj, nk, nj, a_cf, b_cf, c_test_cf, c_ref_cf, num_tested, num_failed, "complex<float>", 1e-5);
                if (nk >= nj) test_shape_c(ni, nj, nk, nk, a_cf, b_cf, c_test_cf, c_ref_cf, num_tested, num_failed, "complex<float>", 1e-5);
            }
        }
    }

    std::cout << "Tested " << num_tested << " complex shapes across double and float: "
              << (num_failed == 0 ? "ALL PASSED" : "FAILURES DETECTED") << "!\n";

    free(a_cd); free(b_cd); free(c_test_cd); free(c_ref_cd);
    free(a_cf); free(b_cf); free(c_test_cf); free(c_ref_cf);

    SafeMPI::Finalize();
    return num_failed == 0 ? 0 : 1;
}
