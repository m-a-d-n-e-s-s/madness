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
#include <madness/world/madness_exception.h>
#include <madness/tensor/mTxmq.h>
#include <madness/tensor/cblas.h>
#include <complex>
#include <cstring>
#include <atomic>
#include <mutex>
#include <cassert>
#include <cstdint>

#ifdef HAVE_MTXMQ
#include <libxsmm.h>
#endif

#if defined(MTXMQ_PROFILE)
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#endif

namespace madness {

#ifdef HAVE_MTXMQ

namespace {

using namespace mTxmq_detail;

// Dual 3D Lock-Free Atomic Cache Tables for double:
// 1. Dense Table: ldb == dimj
// 2. Strided Table: ldb == dimk
static std::atomic<libxsmm_gemmfunction> s_dense_table[MAX_DIMI + 1][MAX_DIMJ + 1][MAX_DIMK + 1] = {{{}}};
static std::atomic<libxsmm_gemmfunction> s_strided_table[MAX_DIMI + 1][MAX_DIMJ + 1][MAX_DIMK + 1] = {{{}}};

// Dual 3D Lock-Free Atomic Cache Tables for float:
static std::atomic<libxsmm_gemmfunction> s_float_dense_table[MAX_DIMI + 1][MAX_DIMJ + 1][MAX_DIMK + 1] = {{{}}};
static std::atomic<libxsmm_gemmfunction> s_float_strided_table[MAX_DIMI + 1][MAX_DIMJ + 1][MAX_DIMK + 1] = {{{}}};

static std::atomic<bool> s_initialized{false};
static std::mutex s_init_mutex;

#if defined(MTXMQ_PROFILE)
static std::atomic<uint64_t> s_dense_counts[MAX_DIMI + 1][MAX_DIMJ + 1][MAX_DIMK + 1] = {{{}}};
static std::atomic<uint64_t> s_strided_counts[MAX_DIMI + 1][MAX_DIMJ + 1][MAX_DIMK + 1] = {{{}}};

struct GeneralShapeKey {
    long dimi, dimj, dimk, ldb;
    int tier; // 1 = Table, 2 = Dynamic JIT, 3 = Fallback BLAS
    bool operator<(const GeneralShapeKey& o) const {
        if (dimi != o.dimi) return dimi < o.dimi;
        if (dimj != o.dimj) return dimj < o.dimj;
        if (dimk != o.dimk) return dimk < o.dimk;
        if (ldb != o.ldb) return ldb < o.ldb;
        return tier < o.tier;
    }
};

static std::map<GeneralShapeKey, uint64_t> s_general_counts;
static std::mutex s_profile_mutex;
static std::atomic<bool> s_has_calls{false};
static std::atomic<bool> s_profile_printed{false};
#endif

void init_kernels_internal() {
    std::lock_guard<std::mutex> lock(s_init_mutex);
    if (s_initialized.load(std::memory_order_relaxed)) {
        return;
    }

    libxsmm_init();
    s_initialized.store(true, std::memory_order_release);
}

__attribute__((constructor))
void mTxmq_auto_init() {
    init_kernels_internal();
}

inline libxsmm_gemmfunction dispatch_and_cache(long dimi, long dimj, long dimk, long ldb,
                                               libxsmm_datatype dtype,
                                               std::atomic<libxsmm_gemmfunction>& slot) {
    const libxsmm_gemm_shape shape = libxsmm_create_gemm_shape(
        static_cast<libxsmm_blasint>(dimj),
        static_cast<libxsmm_blasint>(dimi),
        static_cast<libxsmm_blasint>(dimk),
        static_cast<libxsmm_blasint>(ldb),
        static_cast<libxsmm_blasint>(dimi),
        static_cast<libxsmm_blasint>(dimj),
        dtype, dtype, dtype, dtype
    );
    const libxsmm_bitfield flags = LIBXSMM_GEMM_FLAG_TRANS_B | LIBXSMM_GEMM_FLAG_BETA_0;
    libxsmm_gemmfunction kernel = libxsmm_dispatch_gemm(shape, flags, LIBXSMM_GEMM_PREFETCH_NONE);
    slot.store(kernel, std::memory_order_release);
    return kernel;
}

#if defined(MTXMQ_PROFILE)
std::string format_with_commas(uint64_t val) {
    std::string s = std::to_string(val);
    int n = static_cast<int>(s.length()) - 3;
    while (n > 0) {
        s.insert(n, ",");
        n -= 3;
    }
    return s;
}

struct ShapeRecord {
    long dimi, dimj, dimk, ldb;
    int tier;
    uint64_t count;
    double flops;
};

void render_profile_stream(std::ostream& os) {
    std::vector<ShapeRecord> records;
    uint64_t total_calls = 0;
    double total_flops = 0;
    uint64_t t1 = 0, t2 = 0, t3 = 0;

    for (long i = 1; i <= MAX_DIMI; ++i) {
        for (long j = 1; j <= MAX_DIMJ; ++j) {
            for (long k = 1; k <= MAX_DIMK; ++k) {
                uint64_t cnt = s_dense_counts[i][j][k].load(std::memory_order_relaxed);
                if (cnt > 0) {
                    double flops = 2.0 * i * j * k * cnt;
                    records.push_back({i, j, k, j, 1, cnt, flops});
                    total_calls += cnt;
                    total_flops += flops;
                    t1 += cnt;
                }
            }
        }
    }

    for (long i = 1; i <= MAX_DIMI; ++i) {
        for (long j = 1; j <= MAX_DIMJ; ++j) {
            for (long k = 1; k <= MAX_DIMK; ++k) {
                if (k == j) continue;
                uint64_t cnt = s_strided_counts[i][j][k].load(std::memory_order_relaxed);
                if (cnt > 0) {
                    double flops = 2.0 * i * j * k * cnt;
                    records.push_back({i, j, k, k, 1, cnt, flops});
                    total_calls += cnt;
                    total_flops += flops;
                    t1 += cnt;
                }
            }
        }
    }

    {
        std::lock_guard<std::mutex> lock(s_profile_mutex);
        for (const auto& [key, cnt] : s_general_counts) {
            if (cnt > 0) {
                double flops = 2.0 * key.dimi * key.dimj * key.dimk * cnt;
                records.push_back({key.dimi, key.dimj, key.dimk, key.ldb, key.tier, cnt, flops});
                total_calls += cnt;
                total_flops += flops;
                if (key.tier == 2) t2 += cnt;
                else if (key.tier == 3) t3 += cnt;
                else t1 += cnt;
            }
        }
    }

    if (total_calls == 0) return;

    std::sort(records.begin(), records.end(), [](const ShapeRecord& a, const ShapeRecord& b) {
        if (a.count != b.count) return a.count > b.count;
        return a.flops > b.flops;
    });

    os << "\n========================================================================================================================\n";
    os << "                                      madness::mTxmq Execution & Shape Call Profile                                     \n";
    os << "========================================================================================================================\n";
    os << std::left  << std::setw(6)  << " Rank"
       << std::setw(21) << "Matrix Shape (ixjxk)"
       << std::setw(6)  << "ldb"
       << std::setw(9)  << "Tier"
       << std::right
       << std::setw(14) << "Calls"
       << std::setw(10) << "% Calls"
       << std::setw(10) << "Cum %"
       << std::setw(12) << "FLOP/call"
       << std::setw(15) << "GFLOP"
       << std::setw(10) << "% FLOP" << "\n";
    os << std::string(120, '-') << "\n";

    const size_t displayed = std::min<size_t>(records.size(), 40);
    double cum_calls = 0.0;
    for (size_t idx = 0; idx < displayed; ++idx) {
        const ShapeRecord& r = records[idx];
        const double pct_calls = (100.0 * r.count) / total_calls;
        cum_calls += pct_calls;
        const uint64_t flops_per_call =
            static_cast<uint64_t>(2) * r.dimi * r.dimj * r.dimk;
        const double gflops = r.flops * 1e-9;
        const double pct_flops = total_flops > 0 ? (100.0 * r.flops) / total_flops : 0.0;
        const char* tier_str = (r.tier == 1) ? "Table" : ((r.tier == 2) ? "DynJIT" : "BLAS");
        char shape_buf[32];
        snprintf(shape_buf, sizeof(shape_buf), "(%ld x %ld x %ld)", r.dimi, r.dimj, r.dimk);

        os << std::left  << std::setw(6)  << (idx + 1)
           << std::setw(21) << shape_buf
           << std::setw(6)  << r.ldb
           << std::setw(9)  << tier_str
           << std::right
           << std::setw(14) << format_with_commas(r.count)
           << std::fixed << std::setprecision(2)
           << std::setw(9)  << pct_calls << "%"
           << std::setw(9)  << cum_calls << "%"
           << std::setw(12) << format_with_commas(flops_per_call)
           << std::setw(15) << std::setprecision(3) << gflops
           << std::setw(9)  << std::setprecision(2) << pct_flops << "%\n";
    }

    if (records.size() > displayed) {
        uint64_t other_count = 0;
        double other_flops = 0;
        for (size_t idx = displayed; idx < records.size(); ++idx) {
            other_count += records[idx].count;
            other_flops += records[idx].flops;
        }
        double pct_calls = (100.0 * other_count) / total_calls;
        double pct_flops = total_flops > 0 ? (100.0 * other_flops) / total_flops : 0.0;
        std::string other_label = "Other (" + std::to_string(records.size() - displayed) + " shapes)";

        os << std::left << "   -  "
           << std::setw(21) << other_label
           << std::setw(6)  << "-"
           << std::setw(9)  << "Various"
           << std::right
           << std::setw(14) << format_with_commas(other_count)
           << std::fixed << std::setprecision(2)
           << std::setw(9)  << pct_calls << "%"
           << std::setw(9)  << "100.00%"
           << std::setw(12) << "-"
           << std::setw(15) << std::setprecision(3) << (other_flops * 1e-9)
           << std::setw(9)  << std::setprecision(2) << pct_flops << "%\n";
    }

    os << std::string(120, '=') << "\n";
    os << " Total Calls: " << format_with_commas(total_calls)
       << " | Total Compute: " << std::fixed << std::setprecision(3) << (total_flops * 1e-9) << " GFLOP"
       << " | Unique Shapes: " << records.size() << "\n";
    os << " Dispatch Breakdown: Tier 1 (Fast Table): " << format_with_commas(t1) << " (" << std::setprecision(1) << (100.0 * t1 / total_calls) << "%)"
       << " | Tier 2 (Dynamic JIT): " << format_with_commas(t2) << " (" << std::setprecision(1) << (100.0 * t2 / total_calls) << "%)"
       << " | Tier 3 (BLAS Fallback): " << format_with_commas(t3) << " (" << std::setprecision(1) << (100.0 * t3 / total_calls) << "%)\n";
    os << std::string(120, '=') << "\n\n";
    os.flush();
}

struct ProfileAutoReporter {
    ~ProfileAutoReporter() {
        if (!s_profile_printed.load(std::memory_order_relaxed)) {
            print_mtxmq_profile();
        }
    }
};

static ProfileAutoReporter s_auto_reporter;
#endif // MTXMQ_PROFILE

} // anonymous namespace

void print_mtxmq_profile() {
#if defined(MTXMQ_PROFILE)
    if (!s_has_calls.load(std::memory_order_relaxed)) return;

    const char* env_prof = std::getenv("MTXMQ_PROFILE");
    if (env_prof) {
        std::string val(env_prof);
        std::transform(val.begin(), val.end(), val.begin(), ::tolower);
        if (val == "0" || val == "off" || val == "false" || val == "no") {
            return;
        }
    }

    s_profile_printed.store(true, std::memory_order_relaxed);

    const char* env_file = std::getenv("MTXMQ_PROFILE_FILE");
    if (env_file && env_file[0] != '\0') {
        std::ofstream ofs(env_file);
        if (ofs.is_open()) {
            render_profile_stream(ofs);
            return;
        }
    }

    render_profile_stream(std::cout);
#endif
}

void reset_mtxmq_profile() {
#if defined(MTXMQ_PROFILE)
    for (long i = 1; i <= MAX_DIMI; ++i) {
        for (long j = 1; j <= MAX_DIMJ; ++j) {
            for (long k = 1; k <= MAX_DIMK; ++k) {
                s_dense_counts[i][j][k].store(0, std::memory_order_relaxed);
                s_strided_counts[i][j][k].store(0, std::memory_order_relaxed);
            }
        }
    }
    {
        std::lock_guard<std::mutex> lock(s_profile_mutex);
        s_general_counts.clear();
    }
    s_has_calls.store(false, std::memory_order_relaxed);
    s_profile_printed.store(false, std::memory_order_relaxed);
#endif
}

void mTxmq_init() {
    if (__builtin_expect(!s_initialized.load(std::memory_order_acquire), 0)) {
        init_kernels_internal();
    }
}

// Explicit specialization for double precision with LIBXSMM
template <>
void mTxmq(long dimi, long dimj, long dimk,
           double* MADNESS_RESTRICT c,
           const double* a,
           const double* b,
           long ldb) {
    if (ldb == -1) ldb = dimj;
    MADNESS_CHECK(ldb >= dimj);
    if (__builtin_expect(dimi <= 0 || dimj <= 0, 0)) return;
    if (__builtin_expect(dimk <= 0, 0)) {
        for (long i = 0; i < dimi * dimj; ++i) c[i] = 0.0;
        return;
    }

    if (__builtin_expect(!s_initialized.load(std::memory_order_acquire), 0)) {
        init_kernels_internal();
    }

    if (dimi <= MAX_DIMI && dimj <= MAX_DIMJ && dimk <= MAX_DIMK) {
        if (ldb == dimj) {
            libxsmm_gemmfunction kernel = s_dense_table[dimi][dimj][dimk].load(std::memory_order_relaxed);
            if (__builtin_expect(kernel == nullptr, 0)) {
                kernel = dispatch_and_cache(dimi, dimj, dimk, ldb, LIBXSMM_DATATYPE_F64, s_dense_table[dimi][dimj][dimk]);
            }
            if (__builtin_expect(kernel != nullptr, 1)) {
#if defined(MTXMQ_PROFILE)
                s_dense_counts[dimi][dimj][dimk].fetch_add(1, std::memory_order_relaxed);
                if (__builtin_expect(!s_has_calls.load(std::memory_order_relaxed), 0)) {
                    s_has_calls.store(true, std::memory_order_relaxed);
                }
#endif
                libxsmm_gemm_param param;
                param.a.primary = const_cast<double*>(b);
                param.b.primary = const_cast<double*>(a);
                param.c.primary = c;
                kernel(&param);
                return;
            }
        } else if (ldb == dimk) {
            libxsmm_gemmfunction kernel = s_strided_table[dimi][dimj][dimk].load(std::memory_order_relaxed);
            if (__builtin_expect(kernel == nullptr, 0)) {
                kernel = dispatch_and_cache(dimi, dimj, dimk, ldb, LIBXSMM_DATATYPE_F64, s_strided_table[dimi][dimj][dimk]);
            }
            if (__builtin_expect(kernel != nullptr, 1)) {
#if defined(MTXMQ_PROFILE)
                s_strided_counts[dimi][dimj][dimk].fetch_add(1, std::memory_order_relaxed);
                if (__builtin_expect(!s_has_calls.load(std::memory_order_relaxed), 0)) {
                    s_has_calls.store(true, std::memory_order_relaxed);
                }
#endif
                libxsmm_gemm_param param;
                param.a.primary = const_cast<double*>(b);
                param.b.primary = const_cast<double*>(a);
                param.c.primary = c;
                kernel(&param);
                return;
            }
        }
    }

    const libxsmm_gemm_shape shape = libxsmm_create_gemm_shape(
        static_cast<libxsmm_blasint>(dimj),
        static_cast<libxsmm_blasint>(dimi),
        static_cast<libxsmm_blasint>(dimk),
        static_cast<libxsmm_blasint>(ldb),
        static_cast<libxsmm_blasint>(dimi),
        static_cast<libxsmm_blasint>(dimj),
        LIBXSMM_DATATYPE_F64, LIBXSMM_DATATYPE_F64,
        LIBXSMM_DATATYPE_F64, LIBXSMM_DATATYPE_F64
    );
    const libxsmm_bitfield flags = LIBXSMM_GEMM_FLAG_TRANS_B | LIBXSMM_GEMM_FLAG_BETA_0;
    libxsmm_gemmfunction kernel = libxsmm_dispatch_gemm(shape, flags, LIBXSMM_GEMM_PREFETCH_NONE);

    if (kernel != nullptr) {
#if defined(MTXMQ_PROFILE)
        {
            std::lock_guard<std::mutex> lock(s_profile_mutex);
            s_general_counts[{dimi, dimj, dimk, ldb, 2}]++;
        }
        if (__builtin_expect(!s_has_calls.load(std::memory_order_relaxed), 0)) {
            s_has_calls.store(true, std::memory_order_relaxed);
        }
#endif
        libxsmm_gemm_param param;
        param.a.primary = const_cast<double*>(b);
        param.b.primary = const_cast<double*>(a);
        param.c.primary = c;
        kernel(&param);
        return;
    }

#if defined(MTXMQ_PROFILE)
    {
        std::lock_guard<std::mutex> lock(s_profile_mutex);
        s_general_counts[{dimi, dimj, dimk, ldb, 3}]++;
    }
    if (__builtin_expect(!s_has_calls.load(std::memory_order_relaxed), 0)) {
        s_has_calls.store(true, std::memory_order_relaxed);
    }
#endif

    cblas::gemm(cblas::NoTrans, cblas::Trans, dimj, dimi, dimk, 1.0, b, ldb, a, dimi, 0.0, c, dimj);
}

// Explicit specialization for float precision with LIBXSMM
template <>
void mTxmq(long dimi, long dimj, long dimk,
           float* MADNESS_RESTRICT c,
           const float* a,
           const float* b,
           long ldb) {
    if (ldb == -1) ldb = dimj;
    MADNESS_CHECK(ldb >= dimj);
    if (__builtin_expect(dimi <= 0 || dimj <= 0, 0)) return;
    if (__builtin_expect(dimk <= 0, 0)) {
        for (long i = 0; i < dimi * dimj; ++i) c[i] = 0.0f;
        return;
    }

    if (__builtin_expect(!s_initialized.load(std::memory_order_acquire), 0)) {
        init_kernels_internal();
    }

    if (dimi <= MAX_DIMI && dimj <= MAX_DIMJ && dimk <= MAX_DIMK) {
        if (ldb == dimj) {
            libxsmm_gemmfunction kernel = s_float_dense_table[dimi][dimj][dimk].load(std::memory_order_relaxed);
            if (__builtin_expect(kernel == nullptr, 0)) {
                kernel = dispatch_and_cache(dimi, dimj, dimk, ldb, LIBXSMM_DATATYPE_F32, s_float_dense_table[dimi][dimj][dimk]);
            }
            if (__builtin_expect(kernel != nullptr, 1)) {
                libxsmm_gemm_param param;
                param.a.primary = const_cast<float*>(b);
                param.b.primary = const_cast<float*>(a);
                param.c.primary = c;
                kernel(&param);
                return;
            }
        } else if (ldb == dimk) {
            libxsmm_gemmfunction kernel = s_float_strided_table[dimi][dimj][dimk].load(std::memory_order_relaxed);
            if (__builtin_expect(kernel == nullptr, 0)) {
                kernel = dispatch_and_cache(dimi, dimj, dimk, ldb, LIBXSMM_DATATYPE_F32, s_float_strided_table[dimi][dimj][dimk]);
            }
            if (__builtin_expect(kernel != nullptr, 1)) {
                libxsmm_gemm_param param;
                param.a.primary = const_cast<float*>(b);
                param.b.primary = const_cast<float*>(a);
                param.c.primary = c;
                kernel(&param);
                return;
            }
        }
    }

    const libxsmm_gemm_shape shape = libxsmm_create_gemm_shape(
        static_cast<libxsmm_blasint>(dimj),
        static_cast<libxsmm_blasint>(dimi),
        static_cast<libxsmm_blasint>(dimk),
        static_cast<libxsmm_blasint>(ldb),
        static_cast<libxsmm_blasint>(dimi),
        static_cast<libxsmm_blasint>(dimj),
        LIBXSMM_DATATYPE_F32, LIBXSMM_DATATYPE_F32,
        LIBXSMM_DATATYPE_F32, LIBXSMM_DATATYPE_F32
    );
    const libxsmm_bitfield flags = LIBXSMM_GEMM_FLAG_TRANS_B | LIBXSMM_GEMM_FLAG_BETA_0;
    libxsmm_gemmfunction kernel = libxsmm_dispatch_gemm(shape, flags, LIBXSMM_GEMM_PREFETCH_NONE);

    if (kernel != nullptr) {
        libxsmm_gemm_param param;
        param.a.primary = const_cast<float*>(b);
        param.b.primary = const_cast<float*>(a);
        param.c.primary = c;
        kernel(&param);
        return;
    }

    cblas::gemm(cblas::NoTrans, cblas::Trans, dimj, dimi, dimk, 1.0f, b, ldb, a, dimi, 0.0f, c, dimj);
}

#else // !HAVE_MTXMQ (Fallback implementations when LIBXSMM is disabled)

void print_mtxmq_profile() {}
void reset_mtxmq_profile() {}
void mTxmq_init() {}

template <>
void mTxmq(long dimi, long dimj, long dimk,
           double* MADNESS_RESTRICT c,
           const double* a,
           const double* b,
           long ldb) {
    if (ldb == -1) ldb = dimj;
    MADNESS_CHECK(ldb >= dimj);
    if (__builtin_expect(dimi <= 0 || dimj <= 0, 0)) return;
    if (__builtin_expect(dimk <= 0, 0)) {
        for (long i = 0; i < dimi * dimj; ++i) c[i] = 0.0;
        return;
    }
    cblas::gemm(cblas::NoTrans, cblas::Trans, dimj, dimi, dimk, 1.0, b, ldb, a, dimi, 0.0, c, dimj);
}

template <>
void mTxmq(long dimi, long dimj, long dimk,
           float* MADNESS_RESTRICT c,
           const float* a,
           const float* b,
           long ldb) {
    if (ldb == -1) ldb = dimj;
    MADNESS_CHECK(ldb >= dimj);
    if (__builtin_expect(dimi <= 0 || dimj <= 0, 0)) return;
    if (__builtin_expect(dimk <= 0, 0)) {
        for (long i = 0; i < dimi * dimj; ++i) c[i] = 0.0f;
        return;
    }
    cblas::gemm(cblas::NoTrans, cblas::Trans, dimj, dimi, dimk, 1.0f, b, ldb, a, dimi, 0.0f, c, dimj);
}

#endif // HAVE_MTXMQ

// Explicit specialization for std::complex<double> (Complex * Complex)
template <>
void mTxmq(long dimi, long dimj, long dimk,
           std::complex<double>* MADNESS_RESTRICT c,
           const std::complex<double>* a,
           const std::complex<double>* b,
           long ldb) {
    if (ldb == -1) ldb = dimj;
    MADNESS_CHECK(ldb >= dimj);
    if (__builtin_expect(dimi <= 0 || dimj <= 0, 0)) return;
    if (__builtin_expect(dimk <= 0, 0)) {
        for (long i = 0; i < dimi * dimj; ++i) c[i] = std::complex<double>(0.0, 0.0);
        return;
    }
    const std::complex<double> one(1.0, 0.0);
    const std::complex<double> zero(0.0, 0.0);
    cblas::gemm(cblas::NoTrans, cblas::Trans, dimj, dimi, dimk, one, b, ldb, a, dimi, zero, c, dimj);
}

// Explicit specialization for std::complex<float> (Complex * Complex)
template <>
void mTxmq(long dimi, long dimj, long dimk,
           std::complex<float>* MADNESS_RESTRICT c,
           const std::complex<float>* a,
           const std::complex<float>* b,
           long ldb) {
    if (ldb == -1) ldb = dimj;
    MADNESS_CHECK(ldb >= dimj);
    if (__builtin_expect(dimi <= 0 || dimj <= 0, 0)) return;
    if (__builtin_expect(dimk <= 0, 0)) {
        for (long i = 0; i < dimi * dimj; ++i) c[i] = std::complex<float>(0.0f, 0.0f);
        return;
    }
    const std::complex<float> one(1.0f, 0.0f);
    const std::complex<float> zero(0.0f, 0.0f);
    cblas::gemm(cblas::NoTrans, cblas::Trans, dimj, dimi, dimk, one, b, ldb, a, dimi, zero, c, dimj);
}

} // namespace madness
