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

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/

/**
 \file mpi_archive.h
 \brief Implements archives to serialize data for MPI.
 \ingroup serialization
*/

#ifndef MADNESS_WORLD_MPI_ARCHIVE_H__INCLUDED
#define MADNESS_WORLD_MPI_ARCHIVE_H__INCLUDED

#include <type_traits>
#include <madness/world/archive.h>
#include <madness/world/world.h>
#include <madness/world/vector_archive.h>

namespace madness {
    namespace archive {

        /// \addtogroup serialization
        /// @{

        /// Archive allowing serialization and point-to-point communication between processes with MPI.
        class MPIRawOutputArchive : public BaseOutputArchive {
            mutable World* world; ///< The world.
            ProcessID dest; ///< The destination process.
            int tag; ///< MPI communication tag.

        public:
            /// Construct an archive for sending data via MPI.

            /// \param[in] world The world.
            /// \param[in] dest The destination process.
            /// \param[in] tag MPI communication tag.
            MPIRawOutputArchive(World& world, const ProcessID& dest, int tag=SafeMPI::MPIAR_TAG)
                    : world(&world), dest(dest), tag(tag) {};

            /// Serialize data and send it to the destination process.

            /// The function only appears (due to \c enable_if) if \c T is
            /// fundamental.
            /// \tparam T The data type to be sent.
            /// \param[in] t Pointer to the data to be sent.
            /// \param[in] n The number of data items to be sent.
            template <class T>
            inline
            typename std::enable_if< is_trivially_serializable<T>::value, void >::type
            store(const T* t, long n) const {
                if (n > 0) {
                    world->mpi.Send(t, n, dest, tag);
                }
            }
        };

        /// Archive allowing deserialization and point-to-point communication between processes with MPI.
        class MPIRawInputArchive : public BaseInputArchive {
            mutable World* world; ///< The world.
            ProcessID src; ///< The source process.
            int tag; ///< MPI communication tag.

        public:
            /// Construct an archive for receiving data via MPI.

            /// \todo Descriptions needed.
            /// \param[in] world The world.
            /// \param[in] src The source process.
            /// \param[in] tag MPI communication tag.
            MPIRawInputArchive(World& world, const ProcessID& src, int tag=SafeMPI::MPIAR_TAG)
                    : world(&world), src(src), tag(tag) {};

            /// Receive data from the source process and deserialize it.

            /// The function only appears (due to \c enable_if) if \c T is
            /// fundamental.
            /// \tparam T The data type to receive.
            /// \param[out] t Pointer to where the data should be stored.
            /// \param[in] n The number of data items to receive.
            template <class T>
            inline
            typename std::enable_if< is_trivially_serializable<T>::value, void >::type
            load(T* t, long n) const {
                if (n > 0) {
                    world->mpi.Recv(t, n, src, tag);
                }
            }
        };

        /// Archive allowing buffering, serialization of data, and point-to-point communication between processes with MPI.
        class MPIOutputArchive : public BaseOutputArchive {
            mutable World* world; ///< The world.
            ProcessID dest; ///< The destination process.
            int tag; ///< MPI communication tag.
            const std::size_t bufsize; ///< Size of the buffer.
            mutable std::vector<unsigned char> v; ///< The buffer.
            madness::archive::VectorOutputArchive var; ///< Archive for storing the buffer.

        public:
            /// Construct an archive for sending data via MPI.

            /// \param[in] world The world.
            /// \param[in] dest The destination process.
            /// \param[in] tag MPI communication tag.
            MPIOutputArchive(World& world, const ProcessID& dest, int tag=SafeMPI::MPIAR_TAG)
                    : world(&world), dest(dest), tag(tag), bufsize(1024*1024), v(), var(v) {
                v.reserve(2*bufsize);
            };

            /// Serialize data and store it in the buffer.

            /// The function only appears (due to \c enable_if) if \c T is
            /// fundamental.
            /// \tparam T The data type to be serialized.
            /// \param[in] t Pointer to the data.
            /// \param[in] n Number of data items to serialize.
            template <class T>
            inline
            typename std::enable_if< is_trivially_serializable<T>::value, void >::type
            store(const T* t, long n) const {
                if (v.size() > bufsize) flush();
                if (n > 0) {
                    var.store(t, n);
                    if (v.size() > bufsize) flush();
                }
            }

            /// Send all data in the buffer to the destination process.

            /// \todo Check out the "?? why ??" comment.
            void flush() const {
                if (v.size()) {
                    world->mpi.Send(v.size(), dest, tag);
                    world->mpi.Send(&v[0], v.size(), dest, tag);
                    v.clear();
                    if (v.capacity() < 2*bufsize)
                        v.reserve(2*bufsize); // ?? why ??
                }
            };

            /// Close the archive (i.e., send any data in the buffer).
            void close() {
                flush();
            };

            /// Destructor. Close the archive first, which may entail sending data.
            ~MPIOutputArchive() {
                close();
            };
        };

        /// Archive allowing buffering, deserialization of data, and point-to-point communication between processes with MPI.
        class MPIInputArchive : public BaseInputArchive {
            mutable World* world; ///< The world.
            ProcessID src; ///< The source process.
            int tag; ///< MPI communication tag.
            mutable std::vector<unsigned char> v; ///< The buffer.
            madness::archive::VectorInputArchive var; ///< Archive for loading the buffer.

        public:
            /// Construct an archive for receiving data via MPI.

            /// \param[in] world The world.
            /// \param[in] src The source process.
            /// \param[in] tag MPI communication tag.
            MPIInputArchive(World& world, const ProcessID& src, int tag=SafeMPI::MPIAR_TAG)
                    : world(&world), src(src), tag(tag), v(), var(v) {};

            /// Deserialize data and store it in the buffer.

            /// The function only appears (due to \c enable_if) if \c T is
            /// fundamental.
            /// \tparam T The data type to be deserialized.
            /// \param[out] t Pointer to the data.
            /// \param[in] n Number of data items to serialize.
            template <class T>
            inline
            typename std::enable_if< is_trivially_serializable<T>::value, void >::type
            load(T* t, long n) const {
                if (n > 0) {
                    if (!var.nbyte_avail()) {
                        var.rewind();
                        std::size_t m;
                        world->mpi.Recv(m, src, tag);
                        v.resize(m);
                        world->mpi.Recv(v.data(), m, src, tag);
                    }
                    var.load(t, n);
                }
            }
        };

        /// Implementation of functions for storing the pre/postamble in MPI archives.

        /// \attention No type checking over MPI streams, for efficiency.
        /// \tparam T The data type.
        template <class T>
        struct ArchivePrePostImpl<MPIRawOutputArchive,T> {
            /// Store the preamble.

            /// \param[in] ar The archive.
            static void preamble_store(const MPIRawOutputArchive& ar) {};

            /// Store the postamble.

            /// \param[in] ar The archive.
            static inline void postamble_store(const MPIRawOutputArchive& ar) {};
        };

        /// Implementation of functions for loading the pre/postamble in MPI archives.

        /// \attention No type checking over MPI streams, for efficiency.
        /// \tparam T The data type.
        template <class T>
        struct ArchivePrePostImpl<MPIRawInputArchive,T> {
            /// Load the preamble.

            /// \param[in] ar The archive.
            static inline void preamble_load(const MPIRawInputArchive& ar) {};

            /// Load the postamble.

            /// \param[in] ar The archive.
            static inline void postamble_load(const MPIRawInputArchive& ar) {};
        };

        /// Implementation of functions for storing the pre/postamble in MPI archives.

        /// \attention No type checking over MPI streams, for efficiency.
        /// \tparam T The data type.
        template <class T>
        struct ArchivePrePostImpl<MPIOutputArchive,T> {
            /// Store the preamble.

            /// \param[in] ar The archive.
            static void preamble_store(const MPIOutputArchive& ar) {};

            /// Store the postamble.

            /// \param[in] ar The archive.
            static inline void postamble_store(const MPIOutputArchive& ar) {};
        };

        /// Implementation of functions for loading the pre/postamble in MPI archives.

        /// \attention No type checking over MPI streams, for efficiency.
        /// \tparam T The data type.
        template <class T>
        struct ArchivePrePostImpl<MPIInputArchive,T> {
            /// Load the preamble.

            /// \param[in] ar The archive.
            static inline void preamble_load(const MPIInputArchive& ar) {};

            /// Load the postamble.

            /// \param[in] ar The archive.
            static inline void postamble_load(const MPIInputArchive& ar) {};
        };

        /// @}
    }
}
#endif // MADNESS_WORLD_MPI_ARCHIVE_H__INCLUDED
