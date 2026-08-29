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
#include <madness/tensor/mTxmq.h>
#include <madness/tensor/mxm.h>
#include <madness/tensor/cblas.h>

#ifdef HAVE_MTXMQ

#include <libxsmm.h>
#include <atomic>
#include <mutex>
#include <cassert>
#include <cstdint>

#if defined(MTXMQ_PROFILE)
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <fstream>
#include <cstdlib>
#endif

namespace madness {

namespace {

using namespace mTxmq_detail;

// Dual 3D Lock-Free Atomic Cache Tables:
// 1. Dense Table: ldb == dimj (1.91 MB)
// 2. Strided Table: ldb == dimk (1.91 MB)
static std::atomic<libxsmm_gemmfunction> s_dense_table[MAX_DIMI + 1][MAX_DIMJ + 1][MAX_DIMK + 1] = {{{}}};
static std::atomic<libxsmm_gemmfunction> s_strided_table[MAX_DIMI + 1][MAX_DIMJ + 1][MAX_DIMK + 1] = {{{}}};

static std::atomic<bool> s_initialized{false};
static std::mutex s_init_mutex;

#if defined(MTXMQ_PROFILE)
// Call statistics counters for dense and strided domains
static std::atomic<uint64_t> s_dense_counts[MAX_DIMI + 1][MAX_DIMJ + 1][MAX_DIMK + 1] = {{{}}};
static std::atomic<uint64_t> s_strided_counts[MAX_DIMI + 1][MAX_DIMJ + 1][MAX_DIMK + 1] = {{{}}};

// Call statistics for general / out-of-domain shapes
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

// Automatically initialize runtime on library load
__attribute__((constructor))
void mTxmq_auto_init() {
    init_kernels_internal();
}

inline libxsmm_gemmfunction dispatch_and_cache(long dimi, long dimj, long dimk, long ldb,
                                               std::atomic<libxsmm_gemmfunction>& slot) {
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
    uint64_t t1 = 0;
    uint64_t t2 = 0;
    uint64_t t3 = 0;

    // Collect dense table counts (Tier 1, ldb == dimj)
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

    // Collect strided table counts (Tier 1, ldb == dimk != dimj)
    for (long i = 1; i <= MAX_DIMI; ++i) {
        for (long j = 1; j <= MAX_DIMJ; ++j) {
            for (long k = 1; k <= MAX_DIMK; ++k) {
                if (k == j) continue; // already counted in dense table
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

    // Collect general / fallback counts (Tier 2 & 3)
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

    // Sort descending by call count
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
       << std::setw(12) << "FLOPs/call"
       << std::setw(15) << "Total GFLOP"
       << std::setw(10) << "% GFLOP\n";
    os << std::string(120, '-') << "\n";

    size_t max_rows = 50;
    const char* env_max = std::getenv("MTXMQ_PROFILE_MAX_ROWS");
    if (env_max) {
        std::string s_max(env_max);
        if (s_max == "all" || s_max == "0") {
            max_rows = records.size();
        } else {
            max_rows = std::stoul(s_max);
        }
    }

    double cum_calls = 0.0;
    size_t displayed = std::min(records.size(), max_rows);
    int rank = 1;

    for (size_t idx = 0; idx < displayed; ++idx) {
        const auto& r = records[idx];
        double pct_calls = (100.0 * r.count) / total_calls;
        cum_calls += pct_calls;
        double gflops = r.flops * 1e-9;
        double pct_flops = total_flops > 0 ? (100.0 * r.flops) / total_flops : 0.0;
        uint64_t flops_per_call = 2ULL * r.dimi * r.dimj * r.dimk;

        std::string shape_str = "(" + std::to_string(r.dimi) + " x " + std::to_string(r.dimj) + " x " + std::to_string(r.dimk) + ")";
        std::string tier_str = "Tier " + std::to_string(r.tier);

        os << std::left << " " << std::setw(5) << rank++
           << std::setw(21) << shape_str
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

// Automatic exit reporter singleton
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

// Explicit specialization for double precision in the madness namespace
template <>
void mTxmq(long dimi, long dimj, long dimk,
           double* MADNESS_RESTRICT c,
           const double* a,
           const double* b,
           long ldb) {
    if (ldb == -1) ldb = dimj;
    if (__builtin_expect(dimi <= 0 || dimj <= 0, 0)) return;
    if (__builtin_expect(dimk <= 0, 0)) {
        for (long i = 0; i < dimi * dimj; ++i) c[i] = 0.0;
        return;
    }

    // Ensure LIBXSMM runtime is initialized with acquire barrier
    if (__builtin_expect(!s_initialized.load(std::memory_order_acquire), 0)) {
        init_kernels_internal();
    }

    // Fast path: Primary domain (dimi <= 400, dimj <= 24, dimk <= 24)
    if (dimi <= MAX_DIMI && dimj <= MAX_DIMJ && dimk <= MAX_DIMK) {
        if (ldb == dimj) {
            // Dense Table Lookup (ldb == dimj)
            libxsmm_gemmfunction kernel = s_dense_table[dimi][dimj][dimk].load(std::memory_order_relaxed);
            if (__builtin_expect(kernel == nullptr, 0)) {
                kernel = dispatch_and_cache(dimi, dimj, dimk, ldb, s_dense_table[dimi][dimj][dimk]);
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
            // Strided Table Lookup (ldb == dimk)
            libxsmm_gemmfunction kernel = s_strided_table[dimi][dimj][dimk].load(std::memory_order_relaxed);
            if (__builtin_expect(kernel == nullptr, 0)) {
                kernel = dispatch_and_cache(dimi, dimj, dimk, ldb, s_strided_table[dimi][dimj][dimk]);
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

    // Dynamic dispatch path: for general dimensions or arbitrary strides
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

    // Fallback path: BLAS fallback
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

} // namespace madness

#endif // HAVE_MTXMQ
