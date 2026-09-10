#include <madness/world/MADworld.h>
#include <madness/world/atomicint.h>
#include <madness/world/deferred_cleanup.h>
#include <madness/world/worlddc.h>

#include <memory>
#include <type_traits>

using namespace madness;

namespace {

AtomicInt output_serializations;
AtomicInt packed_payload_starts;

class RoundRobinPmap final : public WorldDCPmapInterface<int> {
    ProcessID nproc_;

public:
    explicit RoundRobinPmap(World& world) : nproc_(world.size()) {}

    [[nodiscard]] ProcessID owner(const int& key) const override {
        return key % nproc_;
    }
};

// Returns the archive position of the first value in a packed
// vector<pair<int, CountingValue>>. A count-only archive measures the length
// and the first key, so the test does not assume their sizes or padding.
std::size_t first_packed_value_offset() {
    archive::BufferOutputArchive probe;
    const std::size_t length = 0;
    const int key = 0;
    probe & length & key;
    return probe.size();
}

// The `= -1` initializer is necessary. Without it, the type is trivially
// serializable, the archive writes raw bytes, and the counters stay at 0.
struct CountingValue {
    int value = -1;

    CountingValue() = default;
    explicit CountingValue(int value) : value(value) {}

    template <typename Archive>
    void serialize(const Archive& ar) {
        if constexpr (is_output_archive_v<Archive>) {
            ++output_serializations;

            // A value at this position starts a packed vector. The count
            // shows one vector per root, not one broadcast per entry.
            if constexpr (std::is_same_v<Archive, archive::BufferOutputArchive>) {
                static const std::size_t first_value_offset = first_packed_value_offset();
                if (ar.size() == first_value_offset)
                    ++packed_payload_starts;
            }
        }
        ar & value;
    }
};

// Replicates the container and makes sure that each root sent its shard one
// time. Each key k must hold value_scale * k + value_offset.
void replicate_and_check_sends(World& world, WorldContainer<int, CountingValue>& container,
                               int nkeys, int value_scale, int value_offset) {
    world.gop.fence();
    output_serializations = 0;
    packed_payload_starts = 0;

    container.replicate(true);

    MADNESS_CHECK(container.is_replicated());
    MADNESS_CHECK(static_cast<int>(container.size()) == nkeys);
    for (int key = 0; key < nkeys; ++key)
        MADNESS_CHECK(container.find(key).get()->second.value == value_scale * key + value_offset);

    int total_output_serializations = output_serializations;
    int total_packed_payload_starts = packed_payload_starts;
    world.gop.sum(total_output_serializations);
    world.gop.sum(total_packed_payload_starts);
    const int expected = world.size() == 1 ? 0 : 2 * nkeys;
    const int expected_payload_starts = world.size() == 1 ? 0 : 2 * world.size();
    MADNESS_CHECK(total_output_serializations == expected);
    MADNESS_CHECK(total_packed_payload_starts == expected_payload_starts);
    world.gop.fence();
}

void test_snapshot_replication(World& world) {
    const int nkeys = 8 * world.size();
    auto pmap = std::make_shared<RoundRobinPmap>(world);
    WorldContainer<int, CountingValue> container(world, pmap);

    for (int key = 0; key < nkeys; ++key) {
        if (container.owner(key) == world.rank())
            container.replace(key, CountingValue(3 * key + 1));
    }
    world.gop.fence();
    output_serializations = 0;
    packed_payload_starts = 0;

    container.replicate(true);

    MADNESS_CHECK(container.is_replicated());
    MADNESS_CHECK(static_cast<int>(container.size()) == nkeys);
    for (int key = 0; key < nkeys; ++key)
        MADNESS_CHECK(container.find(key).get()->second.value == 3 * key + 1);

    int total_output_serializations = output_serializations;
    int total_packed_payload_starts = packed_payload_starts;
    world.gop.sum(total_output_serializations);
    world.gop.sum(total_packed_payload_starts);
    const int expected = world.size() == 1 ? 0 : 2 * nkeys;
    const int expected_payload_starts = world.size() == 1 ? 0 : 2 * world.size();
    MADNESS_CHECK(total_output_serializations == expected);
    MADNESS_CHECK(total_packed_payload_starts == expected_payload_starts);
    world.gop.fence();
}

void test_idempotent_rank_replication(World& world) {
    const int nkeys = 8 * world.size();
    auto pmap = std::make_shared<RoundRobinPmap>(world);
    WorldContainer<int, CountingValue> container(world, pmap);

    for (int key = 0; key < nkeys; ++key) {
        if (container.owner(key) == world.rank())
            container.replace(key, CountingValue(5 * key + 2));
    }
    world.gop.fence();
    container.replicate(true);
    output_serializations = 0;
    packed_payload_starts = 0;

    container.replicate(true);

    int total_output_serializations = output_serializations;
    int total_packed_payload_starts = packed_payload_starts;
    world.gop.sum(total_output_serializations);
    world.gop.sum(total_packed_payload_starts);
    MADNESS_CHECK(total_output_serializations == 0);
    MADNESS_CHECK(total_packed_payload_starts == 0);
    MADNESS_CHECK(static_cast<int>(container.size()) == nkeys);
    for (int key = 0; key < nkeys; ++key)
        MADNESS_CHECK(container.find(key).get()->second.value == 5 * key + 2);
    world.gop.fence();
}

AtomicInt cleanup_count;

struct CleanupMarker {
    ~CleanupMarker() { ++cleanup_count; }
};

void test_map_only_reset_does_not_claim_replication(World& world) {
    const int nkeys = 8 * world.size();
    auto pmap = std::make_shared<RoundRobinPmap>(world);
    WorldContainer<int, CountingValue> container(world, pmap);

    for (int key = 0; key < nkeys; ++key) {
        if (container.owner(key) == world.rank())
            container.replace(key, CountingValue(7 * key + 3));
    }
    world.gop.fence();

    // A local pmap is not a completed replication. Each rank holds only its
    // shard, so replicate() must send all shards.
    container.reset_pmap_to_local();
    replicate_and_check_sends(world, container, nkeys, 7, 3);
}

void test_clear_resets_rank_replication(World& world) {
    const int nkeys = 8 * world.size();
    auto pmap = std::make_shared<RoundRobinPmap>(world);
    WorldContainer<int, CountingValue> container(world, pmap);

    for (int key = 0; key < nkeys; ++key) {
        if (container.owner(key) == world.rank())
            container.replace(key, CountingValue(9 * key + 1));
    }
    world.gop.fence();
    container.replicate(true);

    // All ranks clear the container, and only rank 0 fills it again. The next
    // replicate() must send the rank 0 entries, not take the fast path.
    container.clear();
    if (world.rank() == 0) {
        for (int key = 0; key < nkeys; ++key)
            container.replace(key, CountingValue(11 * key + 4));
    }
    world.gop.fence();

    container.replicate(true);

    MADNESS_CHECK(container.is_replicated());
    MADNESS_CHECK(static_cast<int>(container.size()) == nkeys);
    for (int key = 0; key < nkeys; ++key)
        MADNESS_CHECK(container.find(key).get()->second.value == 11 * key + 4);
    world.gop.fence();
}

void test_partial_clear_replicates_on_every_rank(World& world) {
    // With one rank, no other rank can send the cleared entries.
    if (world.size() == 1) return;

    const int nkeys = 8 * world.size();
    auto pmap = std::make_shared<RoundRobinPmap>(world);
    WorldContainer<int, CountingValue> container(world, pmap);

    for (int key = 0; key < nkeys; ++key) {
        if (container.owner(key) == world.rank())
            container.replace(key, CountingValue(13 * key + 6));
    }
    world.gop.fence();
    container.replicate(true);

    // clear() is local, so only rank 0 resets its flag. If each rank chooses
    // the fast path alone, rank 0 waits in a broadcast that other ranks skip.
    if (world.rank() == 0) container.clear();

    container.replicate(true);

    MADNESS_CHECK(static_cast<int>(container.size()) == nkeys);
    for (int key = 0; key < nkeys; ++key)
        MADNESS_CHECK(container.find(key).get()->second.value == 13 * key + 6);
    world.gop.fence();
}

void test_redistribute_resets_rank_replication(World& world) {
    const int nkeys = 8 * world.size();
    auto pmap = std::make_shared<RoundRobinPmap>(world);
    WorldContainer<int, CountingValue> container(world, pmap);

    for (int key = 0; key < nkeys; ++key) {
        if (container.owner(key) == world.rank())
            container.replace(key, CountingValue(15 * key + 2));
    }
    world.gop.fence();
    container.replicate(true);

    // redistribute() replaces the container pmap during the call. This copy
    // keeps the old pmap alive until the call returns.
    auto replicated_pmap = container.get_pmap();
    replicated_pmap->redistribute(world, pmap);
    MADNESS_CHECK(container.is_distributed());

    replicate_and_check_sends(world, container, nkeys, 15, 2);
}

void test_coalesced_redistribute_resets_rank_replication(World& world) {
    const int nkeys = 8 * world.size();
    auto pmap = std::make_shared<RoundRobinPmap>(world);
    WorldContainer<int, CountingValue> container(world, pmap);

    for (int key = 0; key < nkeys; ++key) {
        if (container.owner(key) == world.rank())
            container.replace(key, CountingValue(17 * key + 5));
    }
    world.gop.fence();
    container.replicate(true);

    // The fence order is the same as in test_dc.cc.
    world.gop.fence();
    container.redistribute_coalesced_phase1(pmap);
    world.gop.fence();
    container.redistribute_coalesced_phase2(3);
    world.gop.fence();
    MADNESS_CHECK(container.is_distributed());

    replicate_and_check_sends(world, container, nkeys, 17, 5);
}

void test_idempotent_fence_argument(World& world) {
    auto pmap = std::make_shared<RoundRobinPmap>(world);
    WorldContainer<int, CountingValue> container(world, pmap);
    container.replicate(true);

    cleanup_count = 0;
    auto marker = std::make_shared<CleanupMarker>();
    detail::deferred_cleanup(world, marker, true);
    marker.reset();

    container.replicate(false);
    MADNESS_CHECK(static_cast<int>(cleanup_count) == 0);

    container.replicate(true);
    MADNESS_CHECK(static_cast<int>(cleanup_count) == 1);
}

}  // namespace

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    test_snapshot_replication(world);
    test_idempotent_rank_replication(world);
    test_map_only_reset_does_not_claim_replication(world);
    test_clear_resets_rank_replication(world);
    test_partial_clear_replicates_on_every_rank(world);
    test_redistribute_resets_rank_replication(world);
    test_coalesced_redistribute_resets_rank_replication(world);
    test_idempotent_fence_argument(world);
    finalize();
    return 0;
}
