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

#ifndef MADNESS_WORLD_WORLDAM_H__INCLUDED
#define MADNESS_WORLD_WORLDAM_H__INCLUDED

/// \file worldam.h
/// \brief Implements active message layer for World on top of RMI layer

#include <madness/world/buffer_archive.h>
#include <madness/world/worldrmi.h>
#include <madness/world/world.h>
#include <vector>
#include <cstddef>
#include <memory>
#include <pthread.h>

namespace madness {

    /*
      The RMI layer just does transport and does not know about World
      or even necessarily about MPI.  It also has no buffering or
      virtualization of resources.  In particular, we must be careful
      about having too many outstanding sends and active message
      handlers must be careful about what they do in order to avoid
      deadlock --- especially problematic is a handler trying to send
      a message.

      The WorldAM class provides a World-aware RMI capability that
      limits the number of outstanding sends and can optionally manage
      buffers.  The issue of what handlers can do safely is handled by
      the integration of WorldAM and the task queues ... if you have
      an operation that might

      - send messages

      - take a long time

      - consume a lot of stack/heap (e.g., recursive algorithm)

      then the right thing to do is to send a task rather than
      an active message.
     */

    template <class Derived> class WorldObject;

    class AmArg;
    /// Type of AM handler functions
    typedef void (*am_handlerT)(const AmArg&);

    /// World active message that extends an RMI message
    class AmArg {
    private:
        friend class WorldAmInterface;
        template <class Derived> friend class WorldObject;

        friend AmArg* alloc_am_arg(std::size_t nbyte);

        unsigned char header[RMI::HEADER_LEN]; // !!!!!!!!!  MUST BE FIRST !!!!!!!!!!
        std::size_t nbyte;      // Size of user payload
        unsigned long worldid;  // Id of associated world
        std::ptrdiff_t func;    // User function to call, as a relative fn ptr (see archive::to_rel_fn_ptr)
        ProcessID src;          // Rank of process sending the message
        unsigned int flags;     // Misc. bit flags

        // On 32 bit machine AmArg is HEADER_LEN+4+4+4+4+4=84 bytes
        // On 64 bit machine AmArg is HEADER_LEN+8+8+8+4+4=96 bytes

        // No copy constructor or assignment
        AmArg(const AmArg&);
        AmArg& operator=(const AmArg&);

        void set_src(ProcessID source) { src = source; }

        void set_worldid(unsigned long id) { worldid = id; }

        void set_func(am_handlerT handler) {
            MADNESS_ASSERT(handler);
            func = archive::to_rel_fn_ptr(handler);
        }

        void set_size(std::size_t numbyte) { nbyte = numbyte; }

        void set_pending() { flags |= 0x1ul; }

        bool is_pending() const { return flags & 0x1ul; }

        void clear_flags() { flags = 0; }

        am_handlerT get_func() const { return archive::to_abs_fn_ptr<am_handlerT>(func); }

        archive::BufferInputArchive make_input_arch() const {
            return archive::BufferInputArchive(buf(),size());
        }

        archive::BufferOutputArchive make_output_arch() const {
            return archive::BufferOutputArchive(buf(),size());
        }

    public:
        AmArg() {}

        /// Returns a pointer to the user's payload (aligned in same way as AmArg)
        unsigned char* buf() const { return (unsigned char*)(this) + sizeof(AmArg); }

        /// Returns the size of the user's payload
        std::size_t size() const { return nbyte; }

        /// Used to deserialize arguments from incoming message
        template <typename T>
        archive::BufferInputArchive operator&(T& t) const {
            return make_input_arch() & t;
        }

        /// Used to serialize arguments into outgoing message
        template <typename T>
        archive::BufferOutputArchive operator&(const T& t) const {
            return make_output_arch() & t;
        }

        /// For incoming AM gives the source process
        ProcessID get_src() const { return src; }

        // This is not inline in order to keep World opaque.
        /// For incoming AM gives the associated world
        World* get_world() const { return World::world_from_id(worldid); }

        /// Return the world id
        unsigned long get_worldid() const { return worldid; }
    };


    /// Allocates a new AmArg with nbytes of user data ... delete with free_am_arg
    inline AmArg* alloc_am_arg(std::size_t nbyte) {
        std::size_t narg = 1 + (nbyte+sizeof(AmArg)-1)/sizeof(AmArg);
        AmArg *arg = new AmArg[narg];
        arg->set_size(nbyte);
        return arg;
    }


    inline AmArg* copy_am_arg(const AmArg& arg) {
        AmArg* r = alloc_am_arg(arg.size());
        memcpy(reinterpret_cast<void*>(r), &arg, arg.size()+sizeof(AmArg));
        return r;
    }

    /// Frees an AmArg allocated with alloc_am_arg
    inline void free_am_arg(AmArg* arg) {
        //std::cout << " freeing amarg " << (void*)(arg) << " " << pthread_self() << std::endl;
        delete [] arg;
    }

    /// Terminate argument serialization
    template <typename Archive>
    inline void serialize_am_args(Archive&&) { }

    /// Argument serialization
    template <typename Archive, typename T, typename... argT>
    inline void serialize_am_args(Archive&& archive, T&& t, argT&&... args) {
        serialize_am_args(archive & t, std::forward<argT>(args)...);
    }

    /// Convenience template for serializing arguments into a new AmArg
    template <typename... argT>
    inline AmArg* new_am_arg(const argT&... args) {
        // compute size
        archive::BufferOutputArchive count;
        serialize_am_args(count, args...);

        // Serialize arguments
        AmArg* am_args = alloc_am_arg(count.size());
        serialize_am_args(*am_args, args...);
        return am_args;
    }


    /// Implements AM interface
    class WorldAmInterface : private SCALABLE_MUTEX_TYPE {
        friend class WorldGopInterface;
        friend class World;
    private:

#ifdef HAVE_CRAYXT
        static const int DEFAULT_NSEND = 512;
#else
        static const int DEFAULT_NSEND = 128;
#endif

        class SendReq : public SPINLOCK_TYPE, public RMISendReq {
            AmArg* buf;
            RMI::Request req;
            void free() {if (buf) {free_am_arg(buf); buf=0;}}
        public:
            SendReq() : buf(0) {}
            SendReq(AmArg* b, const RMI::Request& r) : buf(b), req(r) {} // lock is NOT set
            void set(AmArg* b, const RMI::Request& r) {buf=b; req=r;} // assumes lock held
            bool TestAndFree() { // assumes lock held if necessary
                if (buf) {
                    bool ok = req.Test(); 
                    if (ok) free(); 
                    return ok;
                }
                else {
                    return true;
                }
            }
            ~SendReq() {free();}
        };

        // Multiple threads are making their way thru here ... must be careful
        // to ensure updates are atomic and consistent

        int nsend;                          ///< Max no. of pending sends
        std::unique_ptr<SendReq []> send_req; ///< Send requests and managed buffers 
        unsigned long worldid;              ///< The world which contains this instance of WorldAmInterface
        const ProcessID rank;
        const int nproc;
        volatile int cur_msg;               ///< Index of next buffer to attempt to use
        volatile unsigned long nsent;       ///< Counts no. of AM sent for purpose of termination detection
        volatile unsigned long nrecv;       ///< Counts no. of AM received for purpose of termination detection

        std::vector<int> map_to_comm_world; ///< Maps rank in current MPI communicator to SafeMPI::COMM_WORLD

        /// This handles all incoming RMI messages for all instances
        static void handler(void *buf, std::size_t nbyte) {
            // It will be singled threaded since only the RMI receiver
            // thread will invoke it ... however note that nrecv will
            // be read by the main thread during fence operations.
            AmArg* arg = static_cast<AmArg*>(buf);
            am_handlerT func = arg->get_func();
            World* w = arg->get_world();
            MADNESS_ASSERT(arg->size() + sizeof(AmArg) == nbyte);
            MADNESS_ASSERT(w);
            MADNESS_ASSERT(func);
            func(*arg);
            w->am.nrecv++;  // Must be AFTER execution of the function
        }

    public:
        WorldAmInterface(World& world);

        virtual ~WorldAmInterface();

        /// Currently a noop
        void fence() {}

        /// Sends a managed non-blocking active message
        void send(ProcessID dest, am_handlerT op, const AmArg* arg,
                  const int attr=RMI::ATTR_ORDERED)
        {
            // Setup the header
            {
                AmArg* argx = const_cast<AmArg*>(arg);

                argx->set_worldid(worldid);
                argx->set_src(rank);
                argx->set_func(op);
                argx->clear_flags(); // Is this the right place for this?
            }

            // Sanity check
            MADNESS_ASSERT(arg->get_world());
            MADNESS_ASSERT(arg->get_func());

            // Map dest from world's communicator to comm_world
            dest = map_to_comm_world[dest];

            // Remaining code refactored to avoid blocking with lock
            // and to enable finer grained calls into MPI send

            if (RMI::get_this_thread_is_server()) {

                // The server thread may be executing a handler that
                // is trying to send a message (e.g., assigning a
                // remote future).  However, the server must not block
                // in here.  It also needs to send messages in order,
                // thus it puts all send_reqs onto a queue which it
                // processes in its main loop using the RMI::send
                // interface.

                lock(); nsent++; unlock(); // This world must still keep track of messages

                RMI::send_req.emplace_back(std::make_unique<SendReq>((AmArg*)(arg), RMI::isend(arg, arg->size()+sizeof(AmArg), dest, handler, attr)));

                //std::cout << "sending message from server " << (void*)(arg) << " " << pthread_self() << " " << p <<  std::endl;

                return;
            }


            // Find a free buffer oldest first (in order to assist
            // with flow control).  Exit loop with a lock on buffer.

            // This design will need nsend >= nthreads
            int i=-1;
            while (i == -1) {
                lock();   // << Protect cur_msg and nsent;
                if (send_req[cur_msg].try_lock()) { // << matching unlock at end of routine
                    i = cur_msg;
                    cur_msg = (cur_msg + 1) % nsend;
                    nsent++;
                }
                unlock(); // << Protect cur_msg and nsent;
            }                


            // If the buffer is still in-use wait for it to complete
            while (!send_req[i].TestAndFree()) {
                // If the oldest message has still not completed then
                // there is likely severe network or end-point
                // congestion, so pause for 100us in a rather
                // arbitrary attempt to decrease the injection rate.
                // Both the server thread and this call to Test()
                // should ensure progress.
                myusleep(100);
            }

            // Buffer is now free but still locked by me
            send_req[i].set((AmArg*)(arg), RMI::isend(arg, arg->size()+sizeof(AmArg), dest, handler, attr));
            send_req[i].unlock(); // << matches try_lock above
        }

        /// Frees as many send buffers as possible, returning the number that are free
        int free_managed_buffers() {
            int nfree = 0;
            for (int i=0; i<nsend; i++) {
                if (send_req[i].try_lock()) { // Someone may be trying to put a message into this buffer
                    if (send_req[i].TestAndFree()) nfree++;
                    send_req[i].unlock();     // matching unlock
                }
            }
            return nfree;
        }

    };
}

#endif // MADNESS_WORLD_WORLDAM_H__INCLUDED
