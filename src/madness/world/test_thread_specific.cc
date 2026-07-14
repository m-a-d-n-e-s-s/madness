/// \file test_thread_specific.cc
/// \brief Unit tests for madness::detail::thread_specific — the reclaimable
///        per-thread storage backing the eval scratch buffers.
///
/// Covers: single-thread identity and seeding, clear()-reseeds (including that
/// clear() invalidates a thread's cached fast-path slot via the generation
/// stamp), and one distinct instance per concurrent thread.

#include <madness/world/MADworld.h>
#include <madness/world/thread_specific.h>

#include <atomic>
#include <cstdio>
#include <mutex>
#include <set>
#include <thread>
#include <utility>
#include <vector>

using madness::detail::thread_specific;

namespace {

int check(bool ok, const char* what) {
    if (!ok) std::printf("  FAIL: %s\n", what);
    return ok ? 0 : 1;
}

int test_single_thread() {
    int e = 0;
    thread_specific<int> ts(7);

    // Seeded from the init value; local() returns the same object each call.
    int& a = ts.local();
    e += check(a == 7, "seed value");
    a = 42;
    e += check(ts.local() == 42, "local() is stable (cached fast path)");
    e += check(ts.size() == 1, "one item after first touch");

    // clear() frees the item; the next local() reseeds.  The calling thread
    // still holds a cached slot from before clear(), so this also exercises the
    // generation stamp that must force it back through the slow path.
    ts.local() = 99;
    ts.clear();
    e += check(ts.size() == 0, "no items after clear()");
    e += check(ts.local() == 7, "reseed after clear() despite stale cached slot");
    e += check(ts.local() == 7, "reseeded item is stable");
    return e;
}

int test_multi_thread() {
    constexpr int N = 8;
    thread_specific<std::pair<std::thread::id, int>> ts;

    std::mutex m;
    std::set<void*> addrs;
    std::vector<std::thread> pool;
    for (int i = 0; i < N; ++i) {
        pool.emplace_back([&, i] {
            auto& item = ts.local();
            item = {std::this_thread::get_id(), i};
            void* p = &ts.local();  // fast path must return the same object
            std::lock_guard<std::mutex> lk(m);
            addrs.insert(p);
        });
    }
    for (auto& th : pool) th.join();

    int e = 0;
    e += check(addrs.size() == std::size_t(N), "distinct item per concurrent thread");
    e += check(ts.size() == std::size_t(N), "one map entry per thread (survives thread exit)");
    return e;
}

}  // namespace

int main(int argc, char** argv) {
    madness::initialize(argc, argv);
    int errors = 0;
    errors += test_single_thread();
    errors += test_multi_thread();
    if (errors == 0) std::printf("test_thread_specific: all passed\n");
    else             std::printf("test_thread_specific: %d failures\n", errors);
    madness::finalize();
    return errors;
}
