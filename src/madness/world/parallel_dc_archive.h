#ifndef MAD_PARALLEL_DC_ARCHIVE_H_INCLUDED
#define MAD_PARALLEL_DC_ARCHIVE_H_INCLUDED

#include <madness/world/MADworld.h>
#include <madness/world/worlddc.h>
#include <madness/world/vector_archive.h>

#include <algorithm>   // std::min
#include <cstdint>     // std::uint64_t
#include <cstring>     // std::memcpy
#include <limits>      // std::numeric_limits

namespace madness {



    namespace archive {
        
        /// Marker for a chunked record header (see ContainerRecordOutputArchive::close).
        ///
        /// A record larger than the RMI INT_MAX message limit cannot be sent through
        /// `WorldContainer::replace` in one piece, so it is split: the record at `key`
        /// becomes a fixed 16-byte header `[CHUNK_MAGIC][total]` and the payload is
        /// stored as `ceil(total/CONTAINER_RECORD_MAX_CHUNK)` records at `key+1, key+2, …`.
        /// The magic distinguishes the header from a legitimate record; load additionally
        /// verifies that chunk `key+1` exists before trusting it.
        static constexpr std::uint64_t CONTAINER_RECORD_CHUNK_MAGIC = 0x4D41444E43484B53ull; // "MADNCHKS"

        /// Largest payload sent through a single WorldContainer::replace.
        ///
        /// The worldrmi assert is on the *total* message size (serialized key + value +
        /// AM header), so we leave 1 MiB of headroom below INT_MAX for that overhead.
        static constexpr std::size_t CONTAINER_RECORD_MAX_CHUNK =
            static_cast<std::size_t>(std::numeric_limits<int>::max()) - (std::size_t(1) << 20);

        class ContainerRecordOutputArchive : public BaseOutputArchive {
        public:
            using keyT = long;
        private:
            using containerT = WorldContainer<keyT,std::vector<unsigned char>>;
            World& subworld;
            keyT key;
            containerT& dc; // lifetime???
            std::vector<unsigned char> v;
            VectorOutputArchive ar;

        public:

            ContainerRecordOutputArchive(World& subworld, containerT& dc, const keyT& key)
                : subworld(subworld)
                , key(key)
                , dc(dc)
                , v()
                , ar(v)
            {}

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
                if (subworld.rank() != 0) return;

                // common case: the serialized record fits in a single RMI message
                if (v.size() <= CONTAINER_RECORD_MAX_CHUNK) {
                    dc.replace(key, v);
                    return;
                }

                // oversized record: write a header at `key` and split the payload into
                // chunks at `key+1, key+2, …` (see CONTAINER_RECORD_CHUNK_MAGIC).
                const std::uint64_t total = v.size();
                const std::uint64_t magic = CONTAINER_RECORD_CHUNK_MAGIC;
                std::vector<unsigned char> header(2 * sizeof(std::uint64_t));
                std::memcpy(header.data(), &magic, sizeof(magic));
                std::memcpy(header.data() + sizeof(magic), &total, sizeof(total));
                dc.replace(key, header);

                keyT chunk_key = key + 1;
                for (std::size_t start = 0; start < v.size();
                     start += CONTAINER_RECORD_MAX_CHUNK, ++chunk_key) {
                    const std::size_t end = std::min(start + CONTAINER_RECORD_MAX_CHUNK, v.size());
                    dc.replace(chunk_key, std::vector<unsigned char>(v.begin() + start, v.begin() + end));
                }
            }
        };
        
        class ContainerRecordInputArchive : public BaseInputArchive {
            using keyT = long;
            using containerT = WorldContainer<keyT,std::vector<unsigned char>>;
            ProcessID rank;
            std::vector<unsigned char> v;
            VectorInputArchive ar;

            // reassemble a chunked record: concatenate the records at key+1, key+2, …
            // until `total` bytes have been collected (see CONTAINER_RECORD_CHUNK_MAGIC).
            static std::vector<unsigned char>
            load_chunks(const containerT& dc, const keyT& key, std::uint64_t total) {
                std::vector<unsigned char> buf;
                buf.reserve(total);
                for (keyT chunk_key = key + 1; buf.size() < total; ++chunk_key) {
                    containerT::const_iterator it = dc.find(chunk_key).get();
                    MADNESS_CHECK_THROW(it != dc.end(),
                        "ContainerRecordInputArchive: chunk record not found");
                    buf.insert(buf.end(), it->second.begin(), it->second.end());
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

                // a chunked header is exactly [magic][total] (16 bytes); confirm the
                // magic and that the first chunk really exists before trusting it, so a
                // legitimate 16-byte record can never be misread as a header.
                if (raw.size() == 2 * sizeof(std::uint64_t)) {
                    std::uint64_t magic = 0, total = 0;
                    std::memcpy(&magic, raw.data(), sizeof(magic));
                    if (magic == CONTAINER_RECORD_CHUNK_MAGIC) {
                        std::memcpy(&total, raw.data() + sizeof(magic), sizeof(total));
                        if (dc.find(key + 1).get() != dc.end()) {
                            return load_chunks(dc, key, total);
                        }
                    }
                }
                return raw;
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
