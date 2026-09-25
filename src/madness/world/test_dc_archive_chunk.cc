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
/// Run -np 1 for the simplest reproduction (single rank => local_size == total > INT_MAX);
/// -np 2 also exercises the MPI_Gather/Gatherv displacement path.

#include <cstdint>
#include <cstdlib>
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

        const char* which = std::getenv("WHICH");          // "record", "container", or both
        if (!which || std::string(which) == "record")
            passed = test_record_chunk(world, record_bytes) && passed;
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
