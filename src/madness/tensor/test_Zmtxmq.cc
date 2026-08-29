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

static double_complex ran_c() {
    static unsigned long seed = 76521;
    seed = seed * 1812433253 + 12345;
    double d1 = double(seed & 0x7fffffff) * 4.6566128752458e-10;
    seed = seed * 1812433253 + 12345;
    double d2 = double(seed & 0x7fffffff) * 4.6566128752458e-10;
    return double_complex(d1, d2);
}

static void ran_fill_c(size_t n, double_complex* a) {
    for (size_t i = 0; i < n; ++i) a[i] = ran_c();
}

static void test_shape_c(long ni, long nj, long nk, long ldb,
                         const double_complex* a, const double_complex* b,
                         double_complex* c_test, double_complex* c_ref,
                         long& num_tested, long& num_failed) {
    if (ldb == -1) ldb = nj;
    long n_c = (ni > 0 && nj > 0) ? ni * nj : 1;
    for (long i = 0; i < n_c; ++i) {
        c_test[i] = double_complex(-999.0, -999.0);
        c_ref[i] = double_complex(-999.0, -999.0);
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
            std::cerr << "FAILED test_Zmtxmq shape (" << ni << ", " << nj << ", " << nk
                      << ", ldb=" << ldb << ") at elem " << i
                      << ": test=" << c_test[i] << " ref=" << c_ref[i]
                      << " diff=" << diff << " tol=" << tol << "\n";
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

    double_complex *a = nullptr, *b = nullptr, *c_test = nullptr, *c_ref = nullptr;
    if (posix_memalign((void**)&a, 64, nimax * nkmax * sizeof(double_complex)) != 0) return 1;
    if (posix_memalign((void**)&b, 64, nkmax * ldbmax * sizeof(double_complex)) != 0) return 1;
    if (posix_memalign((void**)&c_test, 64, nimax * njmax * sizeof(double_complex)) != 0) return 1;
    if (posix_memalign((void**)&c_ref, 64, nimax * njmax * sizeof(double_complex)) != 0) return 1;

    ran_fill_c(nimax * nkmax, a);
    ran_fill_c(nkmax * ldbmax, b);

    long num_tested = 0;
    long num_failed = 0;

    std::cout << "Testing mTxmq (complex double precision)...\n";

    // 1. Full sweep of small matrices: 1 <= ni <= 24, 1 <= nj <= 12, 1 <= nk <= 12
    for (long ni = 1; ni <= 24; ++ni) {
        for (long nj = 1; nj <= 12; ++nj) {
            for (long nk = 1; nk <= 12; ++nk) {
                test_shape_c(ni, nj, nk, nj, a, b, c_test, c_ref, num_tested, num_failed);
            }
        }
    }

    // 2. Representative MADNESS tensor shapes
    std::vector<long> nis = {1, 2, 4, 8, 16, 32, 64, 100, 128, 216, 256, 400};
    std::vector<long> njs = {1, 2, 4, 6, 8, 10, 12, 16, 20, 24};
    std::vector<long> nks = {1, 2, 4, 6, 8, 10, 12, 16, 20, 24};

    for (long ni : nis) {
        for (long nj : njs) {
            for (long nk : nks) {
                test_shape_c(ni, nj, nk, nj, a, b, c_test, c_ref, num_tested, num_failed);
                if (nk >= nj) test_shape_c(ni, nj, nk, nk, a, b, c_test, c_ref, num_tested, num_failed);
            }
        }
    }

    // 3. Edge cases
    test_shape_c(0, 8, 8, 8, a, b, c_test, c_ref, num_tested, num_failed);
    test_shape_c(8, 0, 8, 8, a, b, c_test, c_ref, num_tested, num_failed);
    test_shape_c(8, 8, 0, 8, a, b, c_test, c_ref, num_tested, num_failed);

    std::cout << "Tested " << num_tested << " complex shapes: "
              << (num_failed == 0 ? "ALL PASSED" : "FAILURES DETECTED") << "!\n";

    free(a);
    free(b);
    free(c_test);
    free(c_ref);

    SafeMPI::Finalize();
    return num_failed == 0 ? 0 : 1;
}
