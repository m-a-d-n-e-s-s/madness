/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

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

/// \file test_dc_archive_chunk.cc
/// \brief Regression tests for storing >INT_MAX data through the parallel-container archives.
///
/// Two independent paths are exercised, both relevant to the cloud (Cloud uses
/// ContainerRecord{Output,Input}Archive, and its StoreFunction policy serializes a
/// function's coefficient WorldContainer through ParallelOutputArchive):
///
///   test_record_chunk:   a single record whose serialized size exceeds INT_MAX,
///                        stored/loaded via ContainerRecord{Output,Input}Archive.
///                        Exercises the chunking added in commit 1953b906.
///                        EXPECTED: passes (chunking is correct).
///
///   test_container_store: a WorldContainer whose *total* serialized size exceeds
///                        INT_MAX, stored via ParallelOutputArchive<VectorOutputArchive>
///                        (the path Cloud's StoreFunction uses for Function coeffs).
///                        This hits the int-typed size/offset arithmetic in
///                        worlddc.h ArchiveStoreImpl<ParallelOutputArchive,WorldContainer>:
///                          const int size = local_size;                 // truncates
///                          size_t total_size = offsets.back()+sizes.back(); // int overflow
///                          new unsigned char[total_size];               // -> std::bad_alloc
///                        EXPECTED (before fix): std::bad_alloc. (after fix): byte-exact.
///
///   test_adjacent_keys, test_lookalike, test_unframed, test_erase: the chunk framing
///                        must not touch neighbouring records, must not misread ordinary
///                        records that look like headers, must leave unframed (batch) records
///                        as the bare payload, and must free chunks on erase.
///
/// Run -np 1 for the simplest reproduction (single rank => local_size == total > INT_MAX);
/// -np 2 also exercises the MPI_Gather/Gatherv displacement path.

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <madness/world/MADworld.h>
#include <madness/world/worlddc.h>
#include <madness/world/vector_archive.h>
#include <madness/world/parallel_archive.h>
#include <madness/world/parallel_dc_archive.h>

using namespace madness;

namespace {

using keyT = madness::archive::ContainerRecordOutputArchive::keyT;          // long
using containerT = WorldContainer<keyT, std::vector<unsigned char>>;

/// pmap that puts every key on the last rank, so rank 0 always communicates
class LastRankPmap : public WorldDCPmapInterface<keyT> {
    const ProcessID dest;
public:
    explicit LastRankPmap(World& world) : dest(world.size() - 1) {}
    ProcessID owner(const keyT& /*key*/) const { return dest; }
};

/// deterministic, position-dependent byte pattern (so we need not keep the source around)
inline unsigned char pat(std::size_t i) {
    std::uint64_t x = i * 1099511628211ull + 1469598103934665603ull;
    x ^= x >> 33;
    return static_cast<unsigned char>(x >> 24);
}

constexpr std::size_t INT_MAX_SZ = static_cast<std::size_t>(std::numeric_limits<int>::max());

// ---------------------------------------------------------------------------
// Test 1: single oversized record through the chunking archive (commit 1953b906)
// ---------------------------------------------------------------------------
bool test_record_chunk(World& world, std::size_t nbytes) {
    const bool io = (world.rank() == 0);
    if (io) {
        print("\n[test_record_chunk] ranks", world.size(), " payload", nbytes,
              " >INT_MAX?", nbytes > INT_MAX_SZ);
    }

    std::shared_ptr<WorldDCPmapInterface<keyT>> pmap(new LastRankPmap(world));
    containerT container(world, pmap);
    const keyT record = 1;

    {
        std::vector<unsigned char> data;
        if (io) { data.resize(nbytes); for (std::size_t i = 0; i < nbytes; ++i) data[i] = pat(i); }
        madness::archive::ContainerRecordOutputArchive ar(world, container, record);
        if (io) ar & data;
    }
    world.gop.fence();

    std::vector<unsigned char> readback;
    {
        madness::archive::ContainerRecordInputArchive ar(world, container, record);
        if (io) ar & readback;
    }
    world.gop.fence();

    int ok = 1;
    if (io) {
        if (readback.size() != nbytes) { print("  FAIL: size", readback.size(), "!=", nbytes); ok = 0; }
        else for (std::size_t i = 0; i < nbytes; ++i)
            if (readback[i] != pat(i)) { print("  FAIL: byte", i); ok = 0; break; }
        if (ok) print("  test_record_chunk: byte-exact OK");
    }
    world.gop.broadcast(ok, 0);
    return ok == 1;
}

/// pmap that puts every key on rank 0, so its stores and finds stay local
class FirstRankPmap : public WorldDCPmapInterface<keyT> {
public:
    ProcessID owner(const keyT& /*key*/) const { return 0; }
};

/// store `data` (significant on rank 0 only) as one record
void store_record(World& world, containerT& c, keyT record,
                  const std::vector<unsigned char>& data, bool framed = true) {
    {
        madness::archive::ContainerRecordOutputArchive ar(world, c, record, framed);
        if (world.rank() == 0) ar & data;
    }
    world.gop.fence();
}

/// load one record; meaningful on rank 0 only
std::vector<unsigned char> load_record(World& world, containerT& c, keyT record) {
    std::vector<unsigned char> readback;
    {
        madness::archive::ContainerRecordInputArchive ar(world, c, record);
        if (world.rank() == 0) ar & readback;
    }
    world.gop.fence();
    return readback;
}

std::vector<unsigned char> pattern(std::size_t n, std::size_t offset = 0) {
    std::vector<unsigned char> d(n);
    for (std::size_t i = 0; i < n; ++i) d[i] = pat(i + offset);
    return d;
}

bool collective_ok(World& world, int ok) {
    world.gop.broadcast(ok, 0);
    return ok == 1;
}

/// global number of entries in `c`
std::size_t global_size(World& world, containerT& c) {
    std::size_t n = c.size();
    world.gop.sum(n);
    return n;
}

// ---------------------------------------------------------------------------
// Test: oversized record at k, ordinary record at k+1
// ---------------------------------------------------------------------------
bool test_adjacent_keys(World& world, std::size_t nbytes) {
    const bool io = (world.rank() == 0);
    if (io) print("\n[test_adjacent_keys] ranks", world.size(), " payload", nbytes);

    std::shared_ptr<WorldDCPmapInterface<keyT>> pmap(new LastRankPmap(world));
    containerT container(world, pmap);
    const keyT big = 4711, small = 4712;                 // as in test_cloud's batch records
    const std::vector<unsigned char> small_data = pattern(1000, 12345);

    // both orders: the neighbour written after the big record, and before it
    store_record(world, container, big, io ? pattern(nbytes) : std::vector<unsigned char>{});
    store_record(world, container, small, io ? small_data : std::vector<unsigned char>{});
    store_record(world, container, big - 1, io ? small_data : std::vector<unsigned char>{});
    store_record(world, container, big, io ? pattern(nbytes) : std::vector<unsigned char>{});

    std::vector<unsigned char> rb_big = load_record(world, container, big);
    std::vector<unsigned char> rb_small = load_record(world, container, small);
    std::vector<unsigned char> rb_prev = load_record(world, container, big - 1);

    int ok = 1;
    if (io) {
        if (rb_small != small_data) { print("  FAIL: record", small, "was overwritten"); ok = 0; }
        if (rb_prev != small_data) { print("  FAIL: record", big - 1, "was overwritten"); ok = 0; }
        if (rb_big.size() != nbytes) { print("  FAIL: big record size", rb_big.size()); ok = 0; }
        else for (std::size_t i = 0; i < nbytes; ++i)
            if (rb_big[i] != pat(i)) { print("  FAIL: big record byte", i); ok = 0; break; }
        if (ok) print("  test_adjacent_keys: OK");
    }
    return collective_ok(world, ok);
}

// ---------------------------------------------------------------------------
// Test: ordinary records that look like chunk headers
// ---------------------------------------------------------------------------
bool test_lookalike(World& world) {
    const bool io = (world.rank() == 0);
    if (io) print("\n[test_lookalike] ranks", world.size());

    std::shared_ptr<WorldDCPmapInterface<keyT>> pmap(new LastRankPmap(world));
    containerT container(world, pmap);

    // the old in-band header [0x4D41444E43484B53][8], with an 8-byte record at key+1,
    // stored as raw 16 bytes (store(ptr,n) writes no length prefix)
    std::vector<unsigned char> old_header(16);
    const std::uint64_t magic = 0x4D41444E43484B53ull, total = 8;
    std::memcpy(old_header.data(), &magic, 8);
    std::memcpy(old_header.data() + 8, &total, 8);
    // the new header layout as a payload: [RECORD_CHUNKED][total][nchunk]
    std::vector<unsigned char> new_header(1 + 16, 0);
    new_header[0] = madness::archive::RECORD_CHUNKED;
    new_header[1] = 1;

    const keyT k = 100;
    for (const auto* payload : {&old_header, &new_header}) {
        {
            madness::archive::ContainerRecordOutputArchive ar(world, container, k);
            if (io) ar.store(payload->data(), long(payload->size()));
        }
        world.gop.fence();
        store_record(world, container, k + 1, io ? pattern(8) : std::vector<unsigned char>{});
        std::vector<unsigned char> rb(payload->size());
        {
            madness::archive::ContainerRecordInputArchive ar(world, container, k);
            if (io) ar.load(rb.data(), long(rb.size()));
        }
        world.gop.fence();
        if (io and rb != *payload) {
            print("  FAIL: look-alike record of", payload->size(), "bytes not read back verbatim");
            return collective_ok(world, 0);
        }
    }
    if (io) print("  test_lookalike: OK");
    return collective_ok(world, 1);
}

// ---------------------------------------------------------------------------
// Test: unframed records (the cloud's batch container) are the bare payload
// ---------------------------------------------------------------------------
bool test_unframed(World& world, std::size_t nbytes) {
    const bool io = (world.rank() == 0);
    if (io) print("\n[test_unframed] ranks", world.size(), " payload", nbytes);

    // owner-local, as store_batch does on the owner: no RMI, so no INT_MAX limit
    std::shared_ptr<WorldDCPmapInterface<keyT>> pmap(new FirstRankPmap());
    containerT container(world, pmap);
    const keyT record = 7;
    {
        madness::archive::ContainerRecordOutputArchive ar(world, container, record, /*framed=*/false);
        if (io) { std::vector<unsigned char> d = pattern(nbytes); ar.store(d.data(), long(d.size())); }
    }
    world.gop.fence();

    int ok = 1;
    if (io) {
        // read the way BatchTransport does: the raw value at `record`
        typename containerT::const_accessor acc;
        if (not container.find(acc, record)) { print("  FAIL: record missing"); ok = 0; }
        else if (acc->second.size() != nbytes) {
            print("  FAIL: stored", acc->second.size(), "bytes, expected the bare", nbytes); ok = 0;
        } else {
            for (std::size_t i = 0; i < nbytes; ++i)
                if (acc->second[i] != pat(i)) { print("  FAIL: byte", i); ok = 0; break; }
        }
        if (ok and container.size() != 1) { print("  FAIL: extra records", container.size()); ok = 0; }
        if (ok) print("  test_unframed: OK");
    }
    return collective_ok(world, ok);
}

// ---------------------------------------------------------------------------
// Test: erasing a chunked record frees its chunks
// ---------------------------------------------------------------------------
bool test_erase(World& world, std::size_t nbytes) {
    const bool io = (world.rank() == 0);
    if (io) print("\n[test_erase] ranks", world.size(), " payload", nbytes);

    std::shared_ptr<WorldDCPmapInterface<keyT>> pmap(new LastRankPmap(world));
    containerT container(world, pmap);
    store_record(world, container, 1, io ? pattern(nbytes) : std::vector<unsigned char>{});
    store_record(world, container, 2, io ? pattern(10) : std::vector<unsigned char>{});
    const std::size_t before = global_size(world, container);
    if (io) madness::archive::erase_container_record(container, 1);
    world.gop.fence();
    const std::size_t after = global_size(world, container);

    int ok = 1;
    if (io) {
        if (before < 3) { print("  FAIL: expected a chunked record, container holds", before); ok = 0; }
        if (after != 1) { print("  FAIL: after erase the container holds", after, "records"); ok = 0; }
        if (ok) print("  test_erase: OK (", before, "->", after, "records )");
    }
    return collective_ok(world, ok);
}

// ---------------------------------------------------------------------------
// Test 2: WorldContainer whose total serialized size exceeds INT_MAX, stored
//         through ParallelOutputArchive (the Cloud StoreFunction path).
// ---------------------------------------------------------------------------
bool test_container_store(World& world, std::size_t total_bytes) {
    const bool io = (world.rank() == 0);
    if (io) {
        print("\n[test_container_store] ranks", world.size(), " total", total_bytes,
              " >INT_MAX?", total_bytes > INT_MAX_SZ);
    }

    // a handful of entries, each well under INT_MAX, summing to > INT_MAX.
    const std::size_t nentry = 5;
    const std::size_t per = total_bytes / nentry;

    containerT src(world);
    if (io) {
        for (std::size_t e = 0; e < nentry; ++e) {
            std::vector<unsigned char> d(per);
            for (std::size_t i = 0; i < per; ++i) d[i] = pat(e * per + i);
            src.replace(static_cast<keyT>(e), d);
        }
    }
    world.gop.fence();

    // store: dispatches to ArchiveStoreImpl<ParallelOutputArchive<VectorOutputArchive>,WorldContainer>
    std::vector<unsigned char> buf;
    {
        madness::archive::VectorOutputArchive var(buf);
        madness::archive::ParallelOutputArchive<madness::archive::VectorOutputArchive> par(world, var);
        par & src;
    }
    world.gop.fence();
    if (io) print("  serialized buffer bytes:", buf.size());

    // load back and verify
    containerT dst(world);
    {
        madness::archive::VectorInputArchive var(buf);
        madness::archive::ParallelInputArchive<madness::archive::VectorInputArchive> par(world, var);
        par & dst;
    }
    world.gop.fence();

    int ok = 1;
    if (io) {
        for (std::size_t e = 0; e < nentry && ok; ++e) {
            auto it = dst.find(static_cast<keyT>(e)).get();
            if (it == dst.end()) { print("  FAIL: entry", e, "missing"); ok = 0; break; }
            const auto& d = it->second;
            if (d.size() != per) { print("  FAIL: entry", e, "size", d.size(), "!=", per); ok = 0; break; }
            for (std::size_t i = 0; i < per; ++i)
                if (d[i] != pat(e * per + i)) { print("  FAIL: entry", e, "byte", i); ok = 0; break; }
        }
        if (ok) print("  test_container_store: byte-exact OK");
    }
    world.gop.broadcast(ok, 0);
    return ok == 1;
}

} // namespace

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    bool passed = true;
    try {
        std::size_t record_bytes = 2200000000ull;   // > INT_MAX (single record / chunking)
        std::size_t total_bytes  = 2400000000ull;   // > INT_MAX (container total / gather)
        if (const char* e = std::getenv("RECORD_SIZE")) record_bytes = std::strtoull(e, nullptr, 10);
        if (const char* e = std::getenv("TOTAL_SIZE"))  total_bytes  = std::strtoull(e, nullptr, 10);

        const char* which = std::getenv("WHICH");          // one test name, or all
        if (!which || std::string(which) == "record")
            passed = test_record_chunk(world, record_bytes) && passed;
        if (!which || std::string(which) == "adjacent")
            passed = test_adjacent_keys(world, record_bytes) && passed;
        if (!which || std::string(which) == "lookalike")
            passed = test_lookalike(world) && passed;
        if (!which || std::string(which) == "unframed")
            passed = test_unframed(world, record_bytes) && passed;
        if (!which || std::string(which) == "erase")
            passed = test_erase(world, record_bytes) && passed;
        if (!which || std::string(which) == "container")
            passed = test_container_store(world, total_bytes) && passed;
    }
    catch (const SafeMPI::Exception& e) { error("caught an MPI exception"); }
    catch (const madness::MadnessException& e) { print(e); error("caught a MADNESS exception"); }
    catch (const std::exception& e) {
        if (world.rank() == 0) print("\n>>> caught std::exception:", e.what(),
                                     "(this is the failure under investigation)");
        passed = false;
    }
    catch (...) { if (world.rank()==0) print("caught unhandled exception"); passed = false; }

    if (world.rank() == 0)
        print("\ntest_dc_archive_chunk", passed ? "passed" : "FAILED");

    finalize();
    return passed ? 0 : 1;
}
