/// \file test_exchange_batch_cache.cc
/// \brief tests the exchange operator's batch record keys and its bounded batch cache

#include <madness/chem/exchangeoperator.h>
#include <madness/world/test_utilities.h>

#include <set>
#include <string>
#include <vector>

using namespace madness;

using cacheT = ExchangeBatchLRU<long, std::string>;

int test_record_keys() {
    test_output t1("testing exchange batch record keys");

    const long salt = 0x1234;
    const Batch_1D r0(0, 5), r1(5, 10);

    t1.checkpoint(exchange_batch_record_key(salt, EXCHANGE_BATCH_VF, r0)
                          == exchange_batch_record_key(salt, EXCHANGE_BATCH_VF, r0),
                  "the same batch derives the same key, so store and fetch agree");

    // the three operand vectors and every batch range must map to distinct records --
    // a collision would silently serve one operand's coefficients as another's
    std::set<long> keys;
    long ncase = 0;
    for (int dim : {EXCHANGE_BATCH_VF, EXCHANGE_BATCH_BRA, EXCHANGE_BATCH_KET}) {
        for (long begin = 0; begin < 40; ++begin) {
            for (long len = 1; len <= 4; ++len) {
                keys.insert(exchange_batch_record_key(salt, dim, Batch_1D(begin, begin + len)));
                ++ncase;
            }
        }
    }
    t1.checkpoint(long(keys.size()) == ncase, "distinct dimension/range pairs give distinct keys");

    t1.checkpoint(exchange_batch_record_key(salt, EXCHANGE_BATCH_VF, r0)
                          != exchange_batch_record_key(salt + 1, EXCHANGE_BATCH_VF, r0),
                  "a different salt gives a different key, so records of one application "
                  "cannot be served to another");
    t1.checkpoint(exchange_batch_record_key(salt, EXCHANGE_BATCH_VF, r0)
                          != exchange_batch_record_key(salt, EXCHANGE_BATCH_VF, r1),
                  "adjacent ranges of the same operand give different keys");

    return t1.end();
}

int test_batch_cache() {
    test_output t1("testing the exchange batch cache");

    cacheT cache;
    cache.set_transient_capacity(2);

    t1.checkpoint(cache.find(1) == nullptr, "a miss on an empty cache returns nothing");

    cache.insert(1, "one", false);
    const std::string* hit = cache.find(1);
    t1.checkpoint(hit != nullptr and *hit == "one", "an inserted batch is found again");

    // transient entries are bounded: inserting a third evicts the least recently used
    cache.insert(2, "two", false);
    cache.insert(3, "three", false);
    t1.checkpoint(cache.n_transient() == 2 and not cache.contains(1) and
                  cache.contains(2) and cache.contains(3),
                  "transient entries are held to the capacity, evicting least-recently-used");

    // ... and a hit renews an entry, so it is no longer the eviction candidate
    cache.find(2);
    cache.insert(4, "four", false);
    t1.checkpoint(cache.contains(2) and cache.contains(4) and not cache.contains(3),
                  "a hit renews an entry, so the untouched one is evicted instead");

    // owned batches are pinned: they survive any amount of transient churn, which is the
    // point -- a rank reuses its own batch across all of its tasks
    cacheT pinned;
    pinned.set_transient_capacity(1);
    pinned.insert(100, "owned-a", true);
    pinned.insert(101, "owned-b", true);
    for (long k = 0; k < 20; ++k) pinned.insert(k, "transient", false);
    t1.checkpoint(pinned.contains(100) and pinned.contains(101),
                  "pinned batches are never evicted, however many transients churn through");
    t1.checkpoint(pinned.n_transient() == 1 and pinned.size() == 3,
                  "pinning does not enlarge the transient budget");
    const std::string* owned = pinned.find(100);
    t1.checkpoint(owned != nullptr and *owned == "owned-a", "a pinned batch is still readable");

    // capacity 0 would mean no transient could ever be held, so it is raised to one
    cacheT tiny;
    tiny.set_transient_capacity(0);
    tiny.insert(7, "seven", false);
    t1.checkpoint(tiny.transient_capacity() == 1 and tiny.contains(7),
                  "a zero capacity is raised to one rather than evicting on insert");

    // references must survive promotion and insertion of other entries, since a caller
    // holds the fetched batch while the next fetch happens
    cacheT stable;
    stable.set_transient_capacity(4);
    const std::string& ref = stable.insert(1, "first", true);
    stable.insert(2, "second", false);
    stable.find(2);
    stable.insert(3, "third", false);
    t1.checkpoint(ref == "first", "a returned reference stays valid as other entries come and go");

    cache.clear();
    t1.checkpoint(cache.size() == 0 and cache.transient_capacity() == 2,
                  "clear drops the entries but keeps the configured capacity");

    return t1.end();
}

int main(int argc, char** argv) {
    madness::initialize(argc, argv);
    int result = 0;
    result += test_record_keys();
    result += test_batch_cache();
    if (result == 0) print("\n --> all tests \033[32m passed \033[0m \n");
    else print("\n --> all tests \033[31m failed \033[0m \n");
    madness::finalize();
    return result;
}
