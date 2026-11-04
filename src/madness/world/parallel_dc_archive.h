#ifndef MAD_PARALLEL_DC_ARCHIVE_H_INCLUDED
#define MAD_PARALLEL_DC_ARCHIVE_H_INCLUDED

#include <madness/world/MADworld.h>
#include <madness/world/worlddc.h>
#include <madness/world/vector_archive.h>
#include <madness/world/worldhash.h>

#include <algorithm>   // std::min
#include <cstdint>     // std::uint64_t
#include <cstring>     // std::memcpy
#include <limits>      // std::numeric_limits

namespace madness {



    namespace archive {
        
        /// Tag byte that starts every framed record (see ContainerRecordOutputArchive::close).
        ///
        /// A record larger than the RMI INT_MAX message limit cannot be sent through
        /// `WorldContainer::replace` in one piece, so it is split: the record at `key`
        /// becomes a header `[RECORD_CHUNKED][total][nchunk]` and the payload is stored as
        /// `nchunk` records `[RECORD_CHUNK][key][i][slice]` at `container_record_chunk_key(key, i)`.
        /// Any other record is `[RECORD_INLINE][payload]`, so the tag alone tells the reader
        /// which it is.
        ///
        /// Chunk keys share the container's key space, so a chunk can still collide with an
        /// unrelated record. That is unlikely (64-bit keys) but not prevented: ruling it out
        /// would need an existence check on the remote owner, and WorldContainer's only remote
        /// lookup, find(), copies the value, which may itself exceed 2 GiB. A collision is
        /// detected instead: each chunk names its record and index, so the reader throws
        /// rather than returning another record's bytes.
        ///
        /// Records are write-once: a key holding a chunked record must be freed with
        /// erase_container_record() before it is written again, or the old chunks stay
        /// allocated (not checked, for the same reason). The cloud never rewrites a key.
        ///
        /// An unframed record (`framed=false`) is the bare payload, never chunked: for
        /// containers whose bytes are read directly rather than through
        /// ContainerRecordInputArchive (the cloud's batch container, which has its own
        /// chunked transport). Framed and unframed records must not share a container.
        enum ContainerRecordTag : unsigned char {
            RECORD_INLINE  = 0x01,
            RECORD_CHUNKED = 0x02,
            RECORD_CHUNK   = 0x03
        };

        /// Largest payload sent through a single WorldContainer::replace.
        ///
        /// The worldrmi assert is on the *total* message size (serialized key + value +
        /// AM header), so we leave 1 MiB of headroom below INT_MAX for that overhead.
        static constexpr std::size_t CONTAINER_RECORD_MAX_CHUNK =
            static_cast<std::size_t>(std::numeric_limits<int>::max()) - (std::size_t(1) << 20);

        /// Key of chunk `ichunk` of the chunked record at `key`: hashed, not `key+1, key+2, …`,
        /// so a chunk cannot land on a neighbouring record whose key a caller chose sequentially.
        ///
        /// A full 64-bit mix (splitmix64's finalizer, a bijection), not worldhash's hash_value,
        /// which folds a long to 32 bits: two records whose keys agreed in those 32 bits would
        /// get identical chunk keys.
        inline long container_record_chunk_key(const long key, const std::uint64_t ichunk) {
            auto mix = [](std::uint64_t z) {
                z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
                z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
                return z ^ (z >> 31);
            };
            return long(mix(mix(std::uint64_t(key) ^ 0x4D41444E43484B53ull) + ichunk)); // "MADNCHKS"
        }

        /// size of a chunk's prefix `[RECORD_CHUNK][key][ichunk]`
        static constexpr std::size_t CONTAINER_RECORD_CHUNK_PREFIX = 1 + sizeof(long) + sizeof(std::uint64_t);

        /// parse a chunk header; false if `raw` is not one
        ///
        /// Not in a nested `detail` namespace: code that says `using namespace madness;
        /// using namespace archive;` would then find `detail::` ambiguous.
        inline bool parse_container_record_chunk_header(const std::vector<unsigned char>& raw,
                                                        std::uint64_t& total, std::uint64_t& nchunk) {
            if (raw.size() != 1 + 2 * sizeof(std::uint64_t) or raw[0] != RECORD_CHUNKED) return false;
            std::memcpy(&total, raw.data() + 1, sizeof(total));
            std::memcpy(&nchunk, raw.data() + 1 + sizeof(total), sizeof(nchunk));
            return true;
        }

        /// Erase a framed record and, if it is chunked, its chunks; erasing only `key` would
        /// leave the chunks allocated.
        inline void erase_container_record(WorldContainer<long, std::vector<unsigned char>>& dc,
                                           const long key) {
            auto it = dc.find(key).get();
            if (it == dc.end()) return;
            std::uint64_t total = 0, nchunk = 0;
            const bool chunked = parse_container_record_chunk_header(it->second, total, nchunk);
            dc.erase(key);
            if (chunked)
                for (std::uint64_t c = 0; c < nchunk; ++c) dc.erase(container_record_chunk_key(key, c));
        }

        class ContainerRecordOutputArchive : public BaseOutputArchive {
        public:
            using keyT = long;
        private:
            using containerT = WorldContainer<keyT,std::vector<unsigned char>>;
            World& subworld;
            keyT key;
            containerT& dc; // lifetime???
            const bool framed;
            bool closed = false;     // the destructor closes too; write the record only once
            std::vector<unsigned char> v;
            VectorOutputArchive ar;

        public:

            /// @param[in] framed  false writes the bare payload, never chunked (see ContainerRecordTag)
            ContainerRecordOutputArchive(World& subworld, containerT& dc, const keyT& key,
                                         const bool framed = true)
                : subworld(subworld)
                , key(key)
                , dc(dc)
                , framed(framed)
                , v()
                , ar(v)
            {
                if (framed and subworld.rank() == 0) {
                    const unsigned char tag = RECORD_INLINE;
                    ar.store(&tag, 1);
                }
            }

            ~ContainerRecordOutputArchive()
            {
                close();
            }

            VectorOutputArchive& get_archive() {
                return ar;
            }
            template <class T>
            inline
            typename std::enable_if< madness::is_trivially_serializable<T>::value, void >::type
            store(const T* t, long n) const {
                MADNESS_CHECK(subworld.rank() == 0);
                ar.store(t, n);
            }

            void open() {}

            void flush() {}

            void close() {
                if (closed) return;
                closed = true;
                if (subworld.rank() != 0) return;

                // common case: the serialized record fits in a single RMI message
                if (not framed or v.size() <= CONTAINER_RECORD_MAX_CHUNK) {
                    dc.replace(key, v);
                    return;
                }

                // oversized record: write a header at `key` and split the payload (which
                // follows the inline tag at v[0]) into tagged chunks (see ContainerRecordTag).
                const std::uint64_t total = v.size() - 1;
                const std::size_t slice = CONTAINER_RECORD_MAX_CHUNK - CONTAINER_RECORD_CHUNK_PREFIX;
                const std::uint64_t nchunk = (total + slice - 1) / slice;
                std::vector<unsigned char> header(1 + 2 * sizeof(std::uint64_t));
                header[0] = RECORD_CHUNKED;
                std::memcpy(header.data() + 1, &total, sizeof(total));
                std::memcpy(header.data() + 1 + sizeof(total), &nchunk, sizeof(nchunk));
                dc.replace(key, header);

                for (std::uint64_t c = 0; c < nchunk; ++c) {
                    const std::size_t start = 1 + c * slice;
                    const std::size_t end = std::min<std::size_t>(start + slice, v.size());
                    std::vector<unsigned char> chunk(CONTAINER_RECORD_CHUNK_PREFIX);
                    chunk[0] = RECORD_CHUNK;
                    std::memcpy(chunk.data() + 1, &key, sizeof(key));
                    std::memcpy(chunk.data() + 1 + sizeof(key), &c, sizeof(c));
                    chunk.insert(chunk.end(), v.begin() + start, v.begin() + end);
                    dc.replace(container_record_chunk_key(key, c), chunk);
                }
            }
        };
        
        class ContainerRecordInputArchive : public BaseInputArchive {
            using keyT = long;
            using containerT = WorldContainer<keyT,std::vector<unsigned char>>;
            ProcessID rank;
            std::vector<unsigned char> v;
            VectorInputArchive ar;

            // reassemble a chunked record: concatenate the payloads of its `nchunk` chunks
            // (see ContainerRecordTag).
            static std::vector<unsigned char>
            load_chunks(const containerT& dc, const keyT& key, std::uint64_t total, std::uint64_t nchunk) {
                std::vector<unsigned char> buf;
                buf.reserve(total);
                for (std::uint64_t c = 0; c < nchunk; ++c) {
                    containerT::const_iterator it = dc.find(container_record_chunk_key(key, c)).get();
                    MADNESS_CHECK_THROW(it != dc.end(), "ContainerRecordInputArchive: chunk record not found");
                    const std::vector<unsigned char>& chunk = it->second;
                    keyT parent = 0;
                    std::uint64_t index = 0;
                    const bool is_chunk = chunk.size() >= CONTAINER_RECORD_CHUNK_PREFIX and chunk[0] == RECORD_CHUNK;
                    if (is_chunk) {
                        std::memcpy(&parent, chunk.data() + 1, sizeof(parent));
                        std::memcpy(&index, chunk.data() + 1 + sizeof(parent), sizeof(index));
                    }
                    MADNESS_CHECK_THROW(is_chunk and parent == key and index == c,
                        "ContainerRecordInputArchive: chunk overwritten by another record (key collision)");
                    buf.insert(buf.end(), chunk.begin() + CONTAINER_RECORD_CHUNK_PREFIX, chunk.end());
                }
                MADNESS_CHECK_THROW(buf.size() == total,
                    "ContainerRecordInputArchive: chunked record size mismatch");
                return buf;
            }

            // fetch the record at `key`, transparently reassembling chunked records
            static std::vector<unsigned char>
            load_raw(const World& subworld, const containerT& dc, const keyT& key) {
                containerT::const_iterator it = dc.find(key).get();
                if (it == dc.end()) {
                    std::cout << "key " << key << " in world " << subworld.id()
                              << "dc.world " << dc.get_world().id() << std::endl;
                    MADNESS_EXCEPTION("record not found", key);
                }
                const std::vector<unsigned char>& raw = it->second;

                // the tag byte says whether this is the record itself or a chunk header
                if (not raw.empty() and raw[0] == RECORD_INLINE)
                    return std::vector<unsigned char>(raw.begin() + 1, raw.end());
                std::uint64_t total = 0, nchunk = 0;
                MADNESS_CHECK_THROW(parse_container_record_chunk_header(raw, total, nchunk),
                    "ContainerRecordInputArchive: untagged record (written with framed=false?)");
                return load_chunks(dc, key, total, nchunk);
            }

        public:
            ContainerRecordInputArchive(World& subworld, const containerT& dc, const keyT& key)
                : rank(subworld.rank())
                , v(rank == 0 ? load_raw(subworld, dc, key) : std::vector<unsigned char>{})
                , ar(v)
            {}

            ~ContainerRecordInputArchive()
            {}
            
            template <class T>
            inline
            typename std::enable_if< madness::is_trivially_serializable<T>::value, void >::type
            load(T* t, long n) const {
                MADNESS_CHECK(rank == 0);
                ar.load(t,n);
            }
            
            void open() {}
            
            void flush() {}
            
            void close() {}
        };


        /// Implementation of functions for storing the pre/postamble in ContainerRecord archives.

        /// \attention No type checking over Vector buffers, for efficiency.
        /// \tparam T The data type.
        template <class T>
        struct ArchivePrePostImpl<ContainerRecordOutputArchive,T> {
            /// Store the preamble.

            /// \param[in] ar The archive.
            static void preamble_store(const ContainerRecordOutputArchive& ar) {};

            /// Store the postamble.

            /// \param[in] ar The archive.
            static inline void postamble_store(const ContainerRecordOutputArchive& ar) {};
        };

        /// Implementation of functions for loading the pre/postamble in ContainerRecord archives.

        /// \attention No type checking over ContainerRecord buffers, for efficiency.
        /// \tparam T The data type.
        template <class T>
        struct ArchivePrePostImpl<ContainerRecordInputArchive,T> {
            /// Load the preamble.

            /// \param[in] ar The archive.
            static inline void preamble_load(const ContainerRecordInputArchive& ar) {};

            /// Load the postamble.

            /// \param[in] ar The archive.
            static inline void postamble_load(const ContainerRecordInputArchive& ar) {};
        };

        // Forward storing to VectorOutputArchive
        template <class keyT, class valueT>
        struct ArchiveStoreImpl< ParallelOutputArchive<ContainerRecordOutputArchive>, WorldContainer<keyT,valueT> > {
            static void store(const ParallelOutputArchive<ContainerRecordOutputArchive>& ar, const WorldContainer<keyT,valueT>& t) {
                std::vector<unsigned char> v;
                VectorOutputArchive dummyar(v,0);
                const int me = ar.get_world()->rank();

                // Need to pass local archive by reference
                ParallelOutputArchive<VectorOutputArchive> par(*(ar.get_world()), (me==0) ? ar.local_archive().get_archive() : dummyar);
                par & t;

            }
        };

        

    }


}

#endif
