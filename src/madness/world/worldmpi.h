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

#ifndef MADNESS_WORLD_WORLDMPI_H__INCLUDED
#define MADNESS_WORLD_WORLDMPI_H__INCLUDED

/**
 \file worldmpi.h
 \brief Implements \c WorldMpiInterface.
 \ingroup mpi
*/

/*
// If include mpi.h BEFORE stdio/iostream should not need undefs
#ifdef SEEK_CUR
#undef SEEK_CUR
#endif
#ifdef SEEK_SET
#undef SEEK_SET
#endif
#ifdef SEEK_END
#undef SEEK_END
#endif
*/

#include <type_traits>
#include <madness/world/safempi.h>
#include <madness/world/worldtypes.h>
#include <cstdlib>

/// \addtogroup mpi
/// @{

#ifdef MADNESS_USE_BSEND_ACKS
/// \todo Verify: Size of the acknowledgment buffer.
#define MADNESS_ACK_BUFF_SIZE 1000
#endif // MADNESS_USE_BSEND_ACKS

/// String description of the MPI thread level.

/// \param[in] level The MPI thread level.
/// \return A string description of the thread level.
#define MPI_THREAD_STRING(level)  \
        ( level==MPI_THREAD_SERIALIZED ? "THREAD_SERIALIZED" : \
            ( level==MPI_THREAD_MULTIPLE ? "THREAD_MULTIPLE" : \
                ( level==MPI_THREAD_FUNNELED ? "THREAD_FUNNELED" : \
                    ( level==MPI_THREAD_SINGLE ? "THREAD_SINGLE" : "THREAD_UNKNOWN" ) ) ) )

namespace madness {

    // Forward declarations
    class World;
    World& initialize(int&, char**&, const SafeMPI::Intracomm&);
    void finalize();

    /// \todo Brief description needed.
    static const Tag DYNAMIC_TAG_BASE = 1024;

    namespace detail {

        class WorldMpiRuntime;

        /// MPI singleton that manages MPI setup and teardown for MADNESS.

        /// MADNESS will call \c WorldMpi::initialize and \c WorldMpi::finalize
        /// to setup and teardown the MPI runtime.
        class WorldMpi {
        private:
            // Friends of MpiWorld
            friend class WorldMpiRuntime;

            /// Pointer to help MADNESS manage MPI.

            /// This shared pointer is used to manage the lifetime of the MPI
            /// within MADNESS. It ensures that MPI is destroyed only after the
            /// last world object is destroyed.
            static std::shared_ptr<WorldMpi> world_mpi;
            static bool own_mpi; ///< \todo Brief description needed.

#ifdef MADNESS_USE_BSEND_ACKS
            /// Acknowledgment buffer.
            static char* mpi_ack_buffer[MADNESS_ACK_BUFF_SIZE];
#endif // MADNESS_USE_BSEND_ACKS

            /// \c WorldMpi constructor.

            /// Initialize the MPI runtime for MADNESS.
            /// \note This constructor is private to prevent incorrect
            ///    initialization. The user should call \c initialize.
            /// \param[in,out] argc The number of command-line arguments to process.
            /// \param[in,out] argv The command-line arguments.
            /// \param requested The requested thread support for MPI runtime
            WorldMpi(int& argc, char**& argv, int requested) {
                if(own_mpi) {
                    // Assume that MADNESS is managing MPI.
                    SafeMPI::Init_thread(argc, argv, requested);
                } else {
                    // MPI has already been initialized, so it is the user's
                    // responsibility to manage MPI and MADNESS world objects.
                    SafeMPI::detail::init_comm_world();
                }

#ifdef MADNESS_USE_BSEND_ACKS
                // Register the acknowlegement buffer for RMI
                SafeMPI::Attach_buffer(mpi_ack_buffer, MADNESS_ACK_BUFF_SIZE);
#endif // MADNESS_USE_BSEND_ACKS
            }

            // Not allowed
            WorldMpi(const WorldMpi&) = delete;
            WorldMpi& operator=(const WorldMpi&) = delete;

        public:

            /// \c WorldMpi destructor.

            /// This will teardown the MPI, SafeMPI.
            ~WorldMpi() {
#ifdef MADNESS_USE_BSEND_ACKS
                // Unregister the acknowlegement buffer for RMI
                void* buff = nullptr;
                SafeMPI::Detach_buffer(buff);
#endif // MADNESS_USE_BSEND_ACKS

                // Teardown MPI/SafeMPI
                if(own_mpi) {
                    const int result = SafeMPI::Finalize();

                    // Check that MPI exited cleanly.
                    if(result != MPI_SUCCESS) {
                        // Print the error message returned by MPI_Finalize().
                        char mpi_error_string[MPI_MAX_ERROR_STRING];
                        int len = 0;
                        if(MPI_Error_string(result, mpi_error_string, &len) != MPI_SUCCESS) {
                                std::strncpy(mpi_error_string, "UNKNOWN MPI ERROR!", MPI_MAX_ERROR_STRING);
                        }
                        std::cout << "!! MPI Error: " << mpi_error_string << "\n";
                    }
                }
            }

            /// Initialize the MPI runtime.

            /// This function starts the MPI runtime. If MPI is already running,
            /// then MADNESS delegate responsibility for MPI to the user. In
            /// either case, MPI thread support is checked to make sure MPI will
            /// play nice with MADNESS.
            ///
            /// \throw madness::Exception When MADNESS has already been initialized.
            /// \throw madness::Exception When MPI has already been finalized.
            /// \throw SafeMPI::Exception When an MPI error occurs.
            ///
            /// \param[in,out] argc The number of command line arguments.
            /// \param[in,out] argv The values of command line arguments.
            /// \param[in] requested The requested thread support for MPI runtime.
            static void initialize(int& argc, char**& argv, int requested) {
                // Check that world_mpi has not been initialized yet and that
                // MPI has not been finalized
                MADNESS_ASSERT(! world_mpi);
                MADNESS_ASSERT(! SafeMPI::Is_finalized());

                // Check for ownership of MPI (user or MADNESS runtime).
                own_mpi = ! SafeMPI::Is_initialized();

                // Initialize the MPI runtime
                world_mpi.reset(new WorldMpi(argc, argv, requested));

                // Check that the thread support provided by MPI matches the
                // requested and required thread support.
                const int provided = SafeMPI::Query_thread();
                const int rank = SafeMPI::COMM_WORLD.Get_rank();
                if((provided < requested) && (rank == 0)) {
                    std::cout << "!! Error: MPI_Init_thread did not provide requested functionality: "
                              << MPI_THREAD_STRING(requested) << " (" << MPI_THREAD_STRING(provided) << "). \n"
                              << "!! Error: The MPI standard makes no guarantee about the correctness of a program in such circumstances. \n"
                              << "!! Error: Please reconfigure your MPI to provide the proper thread support. \n"
                              << std::endl;
                    MPI_Abort(MPI_COMM_WORLD, 1);
                } else if((provided > requested) && (rank == 0)) {
                    std::cout << "!! Warning: MPI_Init_thread provided more than the requested functionality: "
                              << MPI_THREAD_STRING(requested) << " (" << MPI_THREAD_STRING(provided) << "). \n"
                              << "!! Warning: You are likely using an MPI implementation with mediocre thread support. \n"
                              << std::endl;
                }

#if defined(MVAPICH2_VERSION)
                // Check that MVAPICH2 has has the correct thread affinity
                char * mv2_string = nullptr;
                int mv2_affinity = 1; /* this is the default behavior of MVAPICH2 */

                if ((mv2_string = getenv("MV2_ENABLE_AFFINITY")) != nullptr) {
                    mv2_affinity = atoi(mv2_string);
                }

                if (mv2_affinity!=0) {
                    std::cout << "!! Error: You are using MVAPICH2 with affinity enabled, probably by default. \n"
                              << "!! Error: This will cause catastrophic performance issues in MADNESS. \n"
                              << "!! Error: Rerun your job with MV2_ENABLE_AFFINITY=0 \n"
                              << std::endl;
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
#endif // defined(MVAPICH2_VERSION)
            }

            /// Finalize the MPI runtime.

            /// This function starts the teardown process of the MPI runtime.
            /// The actual \c MPI_Finalize will only be called when all the
            /// objects using MPI have been destroyed.
            static void finalize() {
                world_mpi.reset();
            }
        }; // class WorldMpi

        /// MPI runtime reference counter.

        /// This object is used to manage the lifetime of the MPI runtime by
        /// holding a reference to the `WorldMpi::world_mpi` pointer.
        class WorldMpiRuntime {
        private:
            /// A pointer to `WorldMpi::world_mpi`. Used to help manage the lifetime of MPI.
            std::shared_ptr<WorldMpi> world_mpi;

        public:
            /// Constructor.
            WorldMpiRuntime() : world_mpi(WorldMpi::world_mpi) { }

            /// Destructor.
            ~WorldMpiRuntime() { world_mpi.reset(); }
        }; // class WorldMpiInstance

    } // namespace detail


    /// This class wraps/extends the MPI interface for \c World.
    class WorldMpiInterface
        : private detail::WorldMpiRuntime, public SafeMPI::Intracomm
    {

        // Not allowed
        WorldMpiInterface(const WorldMpiInterface&) = delete;
        WorldMpiInterface& operator=(const WorldMpiInterface&) = delete;

    public:
        /// Constructs an interface in the specified \c SafeMPI communicator.

        /// \todo Verify this documentation.
        /// \param[in] comm The communicator.
        WorldMpiInterface(const SafeMPI::Intracomm& comm) :
            detail::WorldMpiRuntime(), SafeMPI::Intracomm(comm)
        { }

        ~WorldMpiInterface() = default;

        /// Returns the associated \c SafeMPI communicator.

        /// \return The associated \c SafeMPI communicator.
        SafeMPI::Intracomm& comm() {
            return *static_cast<SafeMPI::Intracomm*>(this);
        }

        using SafeMPI::Intracomm::Isend;
        using SafeMPI::Intracomm::Irecv;
        using SafeMPI::Intracomm::Send;
        using SafeMPI::Intracomm::Recv;
        using SafeMPI::Intracomm::Bcast;

        // !! All of the routines below call the protected interfaces provided above.
        // !! Please ensure any additional routines follow this convention.

        /// Isend one element.

        /// \note Disabled for pointers to reduce accidental misuse.
        /// \tparam T The type of data to send.
        /// \param[in] datum The element to send.
        /// \param[in] dest The destination process.
        /// \param[in] tag The MPI tag.
        template <typename T>
        typename std::enable_if<!std::is_pointer<T>::value, SafeMPI::Request>::type
        Isend(const T& datum, int dest, int tag=SafeMPI::DEFAULT_SEND_RECV_TAG) const {
            return SafeMPI::Intracomm::Isend(&datum, sizeof(T), MPI_BYTE, dest, tag);
        }

        /// Async receive data of up to \c count elements from process \c source.

        /// \tparam T The type of data to receive.
        /// \param[out] buf Where to put the received data.
        /// \param[in] count The number of data elements to receive.
        /// \param[in] source The source process.
        /// \param[in] tag The MPI tag.
        template <typename T>
        SafeMPI::Request
        Irecv(T* buf, int count, int source, int tag=SafeMPI::DEFAULT_SEND_RECV_TAG) const {
            return SafeMPI::Intracomm::Irecv(buf, count*sizeof(T), MPI_BYTE, source, tag);
        }


        /// Async receive datum from process \c source with default `tag=1`.

        /// \tparam T The type of data to receive.
        /// \param[out] buf Where to put the received datum.
        /// \param[in] source The source process.
        /// \param[in] tag The MPI tag.
        template <typename T>
        typename std::enable_if<!std::is_pointer<T>::value, SafeMPI::Request>::type
        Irecv(T& buf, int source, int tag=SafeMPI::DEFAULT_SEND_RECV_TAG) const {
            return SafeMPI::Intracomm::Irecv(&buf, sizeof(T), MPI_BYTE, source, tag);
        }


        /// Send array of \c lenbuf elements to process \c dest.

        /// \tparam T The type of data to send.
        /// \param[in] buf Pointer to the data.
        /// \param[in] lenbuf The number of data elements to send.
        /// \param[in] dest The destination process.
        /// \param[in] tag The MPI tag.
        template <class T>
        void Send(const T* buf, long lenbuf, int dest, int tag=SafeMPI::DEFAULT_SEND_RECV_TAG) const {
            SafeMPI::Intracomm::Send((void*)buf, lenbuf*sizeof(T), MPI_BYTE, dest, tag);
        }


        /// Send element to process \c dest with default `tag=1001`.

        /// \note Disabled for pointers to reduce accidental misuse.
        /// \tparam T The type of data to send.
        /// \param[in] datum The data element to send.
        /// \param[in] dest The destination process.
        /// \param[in] tag The MPI tag.
        template <typename T>
        typename std::enable_if<!std::is_pointer<T>::value, void>::type
        Send(const T& datum, int dest, int tag=SafeMPI::DEFAULT_SEND_RECV_TAG) const {
            SafeMPI::Intracomm::Send((void*)&datum, sizeof(T), MPI_BYTE, dest, tag);
        }


        /// Receive data of up to \c lenbuf elements from process \c src.

        /// \tparam T The type of data to receive.
        /// \param[out] buf Where to put the received data.
        /// \param[in] lenbuf The maximum number of data elements to receive.
        /// \param[in] src The source process.
        /// \param[in] tag The MPI tag.
        template <typename T>
        void Recv(T* buf, long lenbuf, int src, int tag) const {
            SafeMPI::Intracomm::Recv(buf, lenbuf*sizeof(T), MPI_BYTE, src, tag);
        }

        /// Receive data of up to \c lenbuf elements from process \c dest with \c status.

        /// \tparam T The type of data to receive.
        /// \param[out] buf Where to put the received data.
        /// \param[in] lenbuf The maximum number of data elements to receive.
        /// \param[in] src The source process.
        /// \param[in] tag The MPI tag.
        /// \param[in] status The status.
        template <typename T>
        void Recv(T* buf, long lenbuf, int src, int tag, SafeMPI::Status& status) const {
            SafeMPI::Intracomm::Recv(buf, lenbuf*sizeof(T), MPI_BYTE, src, tag, status);
        }


        /// Receive datum from process \c src.

        /// \tparam T The type of data to receive.
        /// \param[out] buf The received datum.
        /// \param[in] src The source process.
        /// \param[in] tag The MPI tag.
        template <typename T>
        typename std::enable_if<!std::is_pointer<T>::value, void>::type
        Recv(T& buf, int src, int tag=SafeMPI::DEFAULT_SEND_RECV_TAG) const {
            SafeMPI::Intracomm::Recv(&buf, sizeof(T), MPI_BYTE, src, tag);
        }


        /// MPI broadcast an array of \c count elements.

        /// \note Read documentation about interaction of MPI collectives and
        ///    AM/task handling.
        /// \tparam T The type of data to broadcast.
        /// \param[in,out] buffer The data to send (if this is the \c root
        ///    process); otherwise, a buffer to receive the data.
        /// \param[in] count The number of data elements being broadcast.
        /// \param[in] root The process that is sending the data to other
        ///    processes.
        template <typename T>
        void Bcast(T* buffer, int count, int root) const {
            SafeMPI::Intracomm::Bcast(buffer,count*sizeof(T),MPI_BYTE,root);
        }


        /// MPI broadcast a datum.

        /// \note Read documentation about interaction of MPI collectives and
        ///    AM/task handling.
        /// \tparam T The type of data to broadcast.
        /// \param[in,out] buffer The datum to send (if this is the \c root
        ///    process); otherwise, the received datum.
        /// \param[in] root The process that is sending the data to other
        ///    processes.
        template <typename T>
        typename std::enable_if<!std::is_pointer<T>::value, void>::type
        Bcast(T& buffer, int root) const {
            SafeMPI::Intracomm::Bcast(&buffer, sizeof(T), MPI_BYTE,root);
        }

        /// Access the rank of this process.

        /// \return The rank of this process.
        int rank() const { return SafeMPI::Intracomm::Get_rank(); }

        /// Access the total number of processes.

        /// \return The number of processes.
        int nproc() const { return SafeMPI::Intracomm::Get_size(); }

        /// Access the total number of processes.

        /// \return The number of processes.
        int size() const { return SafeMPI::Intracomm::Get_size(); }
    }; // class WorldMpiInterface

}

/// @}

#endif // MADNESS_WORLD_WORLDMPI_H__INCLUDED
