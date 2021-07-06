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

#ifndef MADNESS_WORLD_PARALLEL_ARCHIVE_H__INCLUDED
#define MADNESS_WORLD_PARALLEL_ARCHIVE_H__INCLUDED

/**
 \file parallel_archive.h
 \brief Implements \c ParallelInputArchive and \c ParallelOutputArchive for parallel serialization of data.
 \ingroup serialization
*/

#include <type_traits>
#include <madness/world/archive.h>
#include <madness/world/binary_fstream_archive.h>
#include <madness/world/world.h>
#include <madness/world/worldgop.h>

#include <unistd.h>
#include <cstring>
#include <cstdio>

namespace madness {
    namespace archive {

        /// \addtogroup serialization
        /// @{

        /// Objects that implement their own parallel archive interface should derive from this class.
        class ParallelSerializableObject {};


        /// Base class for input and output parallel archives.

        /// \tparam Archive The local archive. Only tested for \c BinaryFstreamInputArchive and \c BinaryFstreamOutputArchive.
        /// \todo Should this class derive from \c BaseArchive?
        template <typename Archive>
        class BaseParallelArchive {
            World* world; ///< The world.
            mutable Archive ar; ///< The local archive.
            int nio; ///< Number of I/O nodes (always includes node zero).
            bool do_fence; ///< If true (default), a read/write of parallel objects fences before and after I/O.
            char fname[256]; ///< Name of the archive.
            int nclient; ///< Number of clients of this node, including self. Zero if not I/O node.

        public:
            static const bool is_parallel_archive = true; ///< Mark this class as a parallel archive.

            /// Default constructor.
            BaseParallelArchive()
                : world(nullptr), ar(), nio(0), do_fence(true) {}

            /// Returns the process doing I/O for given node.

            /// Currently assigned in a round-robin fashion to the first
            /// \c nio processes, except on IBM BG/P where we use every 64th.
            /// \param[in] rank The node to check.
            /// \return The process doing I/O for process \c rank.
            ProcessID io_node(ProcessID rank) const {
                return rank%nio;
            }

            /// Returns the process doing I/O for this node.

            /// \return The process doing I/O for this node.
            ProcessID my_io_node() const {
                MADNESS_ASSERT(world);
                return io_node(world->rank());
            }

            /// Returns the number of I/O clients for this node, including self (zero if not an I/O node).

            /// \return The number of I/O clients for this node, including self (zero if not an I/O node).
            int num_io_clients() const {
                MADNESS_ASSERT(world);
                return nclient;
            }

            /// Returns true if this node is doing physical I/O.

            /// \return True if this node is doing physical I/O.
            bool is_io_node() const {
                MADNESS_ASSERT(world);
                return world->rank() == my_io_node();
            }

            /// Returns a pointer to the world.

            /// \return A pointer to the world.
            World* get_world() const {
                MADNESS_ASSERT(world);
                return world;
            }

            /// Opens the parallel archive.

            /// \attention When writing to a new archive, the number of writers
            /// specified is used. When reading from an existing archive,
            /// the number of `ionode`s is adjusted to to be the same as
            /// the number that wrote the original archive. Presently,
            /// we don't have logic to handle reading an archive using
            /// fewer processes originally used to write it. If you
            /// want to fix this have a look in worlddc.h for the only
            /// spot that currently needs changing to make that work.
            ///
            /// \note The default number of I/O nodes is one and there is an
            /// arbitrary maximum of 50 set. On IBM BG/P the maximum
            /// is nproc/64.
            /// \param[in] world The world.
            /// \param[in] filename Name of the file.
            /// \param[in] nwriter The number of writers.
            void open(World& world, const char* filename, int nwriter=1) {
                this->world = &world;
                nio = nwriter;
#if defined(HAVE_IBMBGP) || defined(HAVE_IBMBGQ)
                /* Jeff believes that BG is designed to handle up to *
                 * one file per node and I assume no more than 8 ppn */
                int maxio = world.size()/8;
#else
                int maxio = 50;
#endif
                if (nio > maxio) nio = maxio; // Sanity?
                if (nio > world.size()) nio = world.size();

                MADNESS_ASSERT(filename);
                MADNESS_ASSERT(strlen(filename)-1<sizeof(fname));
                strcpy(fname,filename); // Save the filename for later
                char buf[256];
                MADNESS_ASSERT(strlen(filename)+7 <= sizeof(buf));
                sprintf(buf, "%s.%5.5d", filename, world.rank());

                if (world.rank() == 0) {
                    ar.open(buf);
                    ar & nio; // read/write nio from/to the archive
                    MADNESS_ASSERT(nio <= world.size());
                }

                // Ensure all agree on value of nio that may also have changed if reading
                world.gop.broadcast(nio, 0);

                // Other reader/writers can now open the local archive
                if (is_io_node() && world.rank()) {
                    ar.open(buf);
                }

                // Count #client
                ProcessID me = world.rank();
                nclient=0;
                for (ProcessID p=0; p<world.size(); ++p) if (io_node(p) == me) ++nclient;

//                 if (is_io_node()) {
//                     madness::print("I am an IO node with",nclient,"clients and file",buf);
//                 }
//                 else {
//                     madness::print("I am a client served by",my_io_node(),fname);
//                 }
            }

            /// Returns true if the named, unopened archive exists on disk with read access.

            /// This is a collective operation.
            /// \param[in] world The world.
            /// \param[in] filename Name of the file.
            /// \return True if the named, unopened archive exists and is readable.
            static bool exists(World& world, const char* filename) {
                char buf[256];
                MADNESS_ASSERT(strlen(filename)+7 <= sizeof(buf));
                sprintf(buf, "%s.%5.5d", filename, world.rank());
                bool status;
                if (world.rank() == 0)
                    status = (access(buf, F_OK|R_OK) == 0);

                world.gop.broadcast(status);

                return status;
            }

            /// Closes the parallel archive.
            void close() {
                MADNESS_ASSERT(world);
                if (is_io_node()) ar.close();
            }

            /// Returns a reference to the local archive.

            /// \throw MadnessException If not an I/O node.
            /// \return A reference to the local archive.
            Archive& local_archive() const {
                MADNESS_ASSERT(world);
                MADNESS_ASSERT(is_io_node());
                return ar;
            }

            /// Same as `world.gop.broadcast_serializable(obj, root)`.

            /// \tparam objT Type of object to broadcast.
            /// \param[in] obj The object to broadcast.
            /// \param[in] root The root process for broadcasting.
            template <typename objT>
            void broadcast(objT& obj, ProcessID root) const {
                get_world()->gop.broadcast_serializable(obj, root);
            }

            /// Deletes the files associated with the archive of the given name.

            /// Presently assumes a shared file system since process zero does the
            /// deleting.
            /// \param[in] world The world.
            /// \param[in] filename Base name of the file.
            static void remove(World& world, const char* filename) {
                if (world.rank() == 0) {
                    char buf[268];
                    MADNESS_ASSERT(strlen(filename)+7 <= sizeof(buf));
                    for (ProcessID p=0; p<world.size(); ++p) {
                        sprintf(buf, "%s.%5.5d", filename, p);
                        if (::remove(buf)) break;
                    }
                }
            }

            /// Removes the files associated with the current archive.
            void remove() {
                MADNESS_ASSERT(world);
                remove(*world, fname);
            }

            /// Check if we should fence around a read/write operation.

            /// \return True if we should fence; false otherwise.
            bool dofence() const {
                return this->do_fence;
            }

            /// Set the flag for fencing around a read/write operation.

            /// \param[in] dofence True if we should fence; false otherwise.
            void set_dofence(bool dofence) {
                do_fence = dofence;
            }
        };


        /// An archive for storing local or parallel data wrapping a \c BinaryFstreamOutputArchive.

        /// \note Writes of process-local objects only store the data from process zero.
        ///
        /// \note Writes of parallel containers (presently only \c WorldContainer) store all data.
        ///
        /// Each of the server or I/O nodes creates a
        /// \c BinaryFstreamOutputArchive with the name `filename.rank`. Client
        /// processes send their data to servers in a round-robin fashion.
        ///
        /// Process zero records the number of writers so that, when the archive is opened
        /// for reading, the number of readers is forced to match.
        class ParallelOutputArchive : public BaseParallelArchive<BinaryFstreamOutputArchive>, public BaseOutputArchive {
        public:
            /// Default constructor.
            ParallelOutputArchive() {}

            /// Creates a parallel archive for output with given base filename and number of I/O nodes.

            /// \param[in] world The world.
            /// \param[in] filename Base name of the file.
            /// \param[in] nio The number of I/O nodes.
            ParallelOutputArchive(World& world, const char* filename, int nio=1)  {
                open(world, filename, nio);
            }

            /// Flush any data in the archive.
            void flush() {
                if (is_io_node()) local_archive().flush();
            }
        };

        /// An archive for storing local or parallel data, wrapping a \c BinaryFstreamInputArchive.

        /// \note Reads of process-local objects load the values originally stored by process zero,
        /// which is then broadcast to all processes.
        ///
        /// \note Reads of parallel containers (presently only \c WorldContainer) load all data.
        ///
        /// The number of I/O nodes or readers is presently ignored. It is
        /// forced to be the same as the original number of writers and,
        /// therefore, you cannot presently read an archive from a parallel job
        /// with fewer total processes than the number of writers.
        class ParallelInputArchive : public BaseParallelArchive<BinaryFstreamInputArchive>, public  BaseInputArchive {
        public:
            /// Default constructor.
            ParallelInputArchive() {}

            /// Creates a parallel archive for input.

            /// \param[in] world The world.
            /// \param[in] filename Base name of the file.
            /// \param[in] nio The number of writers. Ignored, see above.
            ParallelInputArchive(World& world, const char* filename, int nio=1) {
                open(world, filename, nio);
            }
        };

        /// Disable type info for parallel output archives.

        /// \tparam T The data type.
        template <class T>
        struct ArchivePrePostImpl<ParallelOutputArchive,T> {
            /// Store the preamble for this data type in the parallel archive.

            /// \param[in] ar The archive.
            static void preamble_store(const ParallelOutputArchive& ar) {}

            /// Store the postamble for this data type in the parallel archive.

            /// \param[in] ar The archive.
            static inline void postamble_store(const ParallelOutputArchive& ar) {}
        };

        /// Disable type info for parallel input archives.

        /// \tparam T The data type.
        template <class T>
        struct ArchivePrePostImpl<ParallelInputArchive,T> {
            /// Load the preamble for this data type in the parallel archive.

            /// \param[in] ar The archive.
            static inline void preamble_load(const ParallelInputArchive& ar) {}

            /// Load the postamble for this data type in the parallel archive.

            /// \param[in] ar The archive.
            static inline void postamble_load(const ParallelInputArchive& ar) {}
        };

        /// Specialization of \c ArchiveImpl for parallel output archives.

        /// \attention No type-checking is performed.
        /// \tparam T The data type.
        template <class T>
        struct ArchiveImpl<ParallelOutputArchive, T, std::enable_if_t<!std::is_base_of_v<ParallelSerializableObject, T>>> {
            /// Store the data in the archive.

            /// Serial objects write only from process 0.
            ///
            /// \param[in] ar The parallel archive.
            /// \param[in] t The serial data.
            /// \return The parallel archive.
            static inline
            const ParallelOutputArchive&
            wrap_store(const ParallelOutputArchive& ar, const T& t) {
                if (ar.get_world()->rank()==0) {
                    ar.local_archive() & t;
                }
                return ar;
            }
        };

        /// Specialization of \c ArchiveImpl for parallel input archives.

        /// \attention No type-checking is performed.
        /// \tparam T The data type.
        template <class T>
        struct ArchiveImpl<ParallelInputArchive, T, std::enable_if_t<!std::is_base_of_v<ParallelSerializableObject, T>>> {
            /// Load the data from the archive.

            /// Serial objects are read only from process 0 and then broadcasted.
            ///
            /// \param[in] ar The parallel archive.
            /// \param[out] t Where to put the loaded data.
            /// \return The parallel archive.
            static inline
            const ParallelInputArchive&
            wrap_load(const ParallelInputArchive& ar, const T& t) {
                if (ar.get_world()->rank()==0) {
                    ar.local_archive() & t;
                }
                ar.broadcast(const_cast<T&>(t), 0);
                return ar;
            }
        };


        /// Write the archive array only from process zero.

        /// \tparam T The array data type.
        template <class T>
        struct ArchiveImpl< ParallelOutputArchive, archive_array<T> > {
            /// Store the \c archive_array in the parallel archive.

            /// \param[in] ar The parallel archive.
            /// \param[in] t The array to store.
            /// \return The parallel archive.
            static inline const ParallelOutputArchive& wrap_store(const ParallelOutputArchive& ar, const archive_array<T>& t) {
                if (ar.get_world()->rank() == 0) ar.local_archive() & t;
                return ar;
            }
        };

        /// Read the archive array and broadcast.

        /// \tparam T The array data type.
        template <class T>
        struct ArchiveImpl< ParallelInputArchive, archive_array<T> > {
            /// Load the \c archive_array from the parallel archive and broadcast it.

            /// \param[in] ar The parallel archive.
            /// \param[out] t Where to put the loaded array.
            /// \return The parallel archive.
            static inline const ParallelInputArchive& wrap_load(const ParallelInputArchive& ar, const archive_array<T>& t) {
                if (ar.get_world()->rank() == 0) ar.local_archive() & t;
                ar.broadcast(t, 0);
                return ar;
            }
        };

        /// @}
    }
}

#endif // MADNESS_WORLD_PARALLEL_ARCHIVE_H__INCLUDED
