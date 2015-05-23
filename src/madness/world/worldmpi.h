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

/// \file worldmpi.h
/// \brief Implements WorldMpiInterface
/// \addtogroup mpi
///@{

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


#ifdef MADNESS_USE_BSEND_ACKS
#define MADNESS_ACK_BUFF_SIZE 1000
#endif // MADNESS_USE_BSEND_ACKS

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

    static const Tag DYNAMIC_TAG_BASE = 1024;

    namespace detail {

        class WorldMpiRuntime;

        /// MPI singleton that manages MPI setup and teardown for MADNESS

        /// MADNESS will call \c WorldMpi::initialize and \c WorldMpi::finalize
        /// to setup and teardown the MPI runtime.
        class WorldMpi {
        private:
            // Friends of MpiWorld
            friend class WorldMpiRuntime;

            // This shared pointer is used to manage the lifetime of the MPI
            // within MADNESS. It ensures that MPI is destroyed only after the
            // last world object is destroyed.
            static std::shared_ptr<WorldMpi> world_mpi;
            static bool own_mpi;

#ifdef MADNESS_USE_BSEND_ACKS
            static char* mpi_ack_buffer[MADNESS_ACK_BUFF_SIZE];
#endif // MADNESS_USE_BSEND_ACKS

            /// WorldMpi constructor

            /// Initialize the MPI runtime for MADNESS.
            /// \note This constructor is private to prevent incorrect
            /// initialization.
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
            WorldMpi(const WorldMpi&);
            WorldMpi& operator=(const WorldMpi&);

        public:

            /// WorldMpi destructor

            /// This will teardown the MPI, SafeMPI
            ~WorldMpi() {
#ifdef MADNESS_USE_BSEND_ACKS
                // Unregister the acknowlegement buffer for RMI
                void* buff = NULL;
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

            /// Initialize the MPI runtime

            /// This function starts the MPI runtime. If MPI is already running,
            /// then MADNESS delegate responsibility for MPI to the user. In
            /// either case, MPI thread support is checked to make sure MPI will
            /// play nice with MADNESS.
            /// \param argc The number of command line arguments
            /// \param argv The values of command line arguments
            /// \param requested The requested thread support for MPI runtime
            /// \throw madness::Exception When MADNESS has already been initialized.
            /// \throw madness::Exception When MPI has already been finalized.
            /// \throw SafeMPI::Exception When an MPI error occurs.
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
                char * mv2_string = NULL;
                int mv2_affinity = 1; /* this is the default behavior of MVAPICH2 */

                if ((mv2_string = getenv("MV2_ENABLE_AFFINITY")) != NULL) {
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

            /// Finalize the MPI runtime

            /// This function starts the teardown process of the MPI runtime.
            /// The actual MPI_Finalize will only be called when all the objects
            /// using MPI have been destroyed.
            static void finalize() {
                world_mpi.reset();
            }
        }; // class WorldMpi

        /// MPI runtime reference counter

        /// This object is used to manage the lifetime of the MPI runtime by
        /// holding a reference to the WorldMpi::world_mpi pointer.
        class WorldMpiRuntime {
        private:
            std::shared_ptr<WorldMpi> world_mpi;

        public:
            WorldMpiRuntime() : world_mpi(WorldMpi::world_mpi) { }
            ~WorldMpiRuntime() { world_mpi.reset(); }
        }; // class WorldMpiInstance

    } // namespace detail


    /// This class wraps/extends the MPI interface for World
    class WorldMpiInterface : private detail::WorldMpiRuntime, public SafeMPI::Intracomm {

        // Not allowed
        WorldMpiInterface(const WorldMpiInterface&);
        WorldMpiInterface& operator=(const WorldMpiInterface&);

    public:
        WorldMpiInterface(const SafeMPI::Intracomm& comm) :
            detail::WorldMpiRuntime(), SafeMPI::Intracomm(comm)
        { }

        ~WorldMpiInterface() { }

        /// Returns the associated SafeMPI communicator
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
        /// Isend one element ... disabled for pointers to reduce accidental misuse.
        template <typename T>
        typename std::enable_if<!std::is_pointer<T>::value, SafeMPI::Request>::type
        Isend(const T& datum, int dest, int tag=SafeMPI::DEFAULT_SEND_RECV_TAG) const {
            return SafeMPI::Intracomm::Isend(&datum, sizeof(T), MPI_BYTE, dest, tag);
        }

        /// Async receive data of up to lenbuf elements from process dest
        template <typename T>
        SafeMPI::Request
        Irecv(T* buf, int count, int source, int tag=SafeMPI::DEFAULT_SEND_RECV_TAG) const {
            return SafeMPI::Intracomm::Irecv(buf, count*sizeof(T), MPI_BYTE, source, tag);
        }


        /// Async receive datum from process dest with default tag=1
        template <typename T>
        typename std::enable_if<!std::is_pointer<T>::value, SafeMPI::Request>::type
        Irecv(T& buf, int source, int tag=SafeMPI::DEFAULT_SEND_RECV_TAG) const {
            return SafeMPI::Intracomm::Irecv(&buf, sizeof(T), MPI_BYTE, source, tag);
        }


        /// Send array of lenbuf elements to process dest
        template <class T>
        void Send(const T* buf, long lenbuf, int dest, int tag=SafeMPI::DEFAULT_SEND_RECV_TAG) const {
            SafeMPI::Intracomm::Send((void*)buf, lenbuf*sizeof(T), MPI_BYTE, dest, tag);
        }


        /// Send element to process dest with default tag=1001

        /// Disabled for pointers to reduce accidental misuse.
        template <typename T>
        typename std::enable_if<!std::is_pointer<T>::value, void>::type
        Send(const T& datum, int dest, int tag=SafeMPI::DEFAULT_SEND_RECV_TAG) const {
            SafeMPI::Intracomm::Send((void*)&datum, sizeof(T), MPI_BYTE, dest, tag);
        }


        /// Receive data of up to lenbuf elements from process dest
        template <typename T>
        void Recv(T* buf, long lenbuf, int src, int tag) const {
            SafeMPI::Intracomm::Recv(buf, lenbuf*sizeof(T), MPI_BYTE, src, tag);
        }

        /// Receive data of up to lenbuf elements from process dest with status
        template <typename T>
        void Recv(T* buf, long lenbuf, int src, int tag, SafeMPI::Status& status) const {
            SafeMPI::Intracomm::Recv(buf, lenbuf*sizeof(T), MPI_BYTE, src, tag, status);
        }


        /// Receive datum from process src
        template <typename T>
        typename std::enable_if<!std::is_pointer<T>::value, void>::type
        Recv(T& buf, int src, int tag=SafeMPI::DEFAULT_SEND_RECV_TAG) const {
            SafeMPI::Intracomm::Recv(&buf, sizeof(T), MPI_BYTE, src, tag);
        }


        /// MPI broadcast an array of count elements

        /// NB.  Read documentation about interaction of MPI collectives and AM/task handling.
        template <typename T>
        void Bcast(T* buffer, int count, int root) const {
            SafeMPI::Intracomm::Bcast(buffer,count*sizeof(T),MPI_BYTE,root);
        }


        /// MPI broadcast a datum

        /// NB.  Read documentation about interaction of MPI collectives and AM/task handling.
        template <typename T>
        typename std::enable_if<!std::is_pointer<T>::value, void>::type
        Bcast(T& buffer, int root) const {
            SafeMPI::Intracomm::Bcast(&buffer, sizeof(T), MPI_BYTE,root);
        }

        int rank() const { return SafeMPI::Intracomm::Get_rank(); }

        int nproc() const { return SafeMPI::Intracomm::Get_size(); }

        int size() const { return SafeMPI::Intracomm::Get_size(); }
    }; // class WorldMpiInterface

}

///@}

#endif // MADNESS_WORLD_WORLDMPI_H__INCLUDED
