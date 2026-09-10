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

}  // namespace

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    test_snapshot_replication(world);
    finalize();
    return 0;
}
