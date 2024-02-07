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

#ifndef MADNESS_WORLD_WORLDRMI_H__INCLUDED
#define MADNESS_WORLD_WORLDRMI_H__INCLUDED

#include <madness/world/safempi.h>
#include <madness/world/thread.h>
#include <madness/world/worldtypes.h>
#include <madness/world/archive.h>
#include <sstream>
#include <utility>
#include <list>
#include <memory>
#include <tuple>
#include <pthread.h>
#include <madness/world/print.h>

/*
  There is just one server thread and it is the only one
  messing with the recv buffers, so there is no need for
  mutex on recv related data.

  Multiple threads (including the server) may send hence
  we need to be careful about send-related data.

  When MPI is initialized we need to use init_thread with
  multiple required.

  This RMI service operates only in (a clone of) COMM_WORLD.  It easy enough
  to extend to other communicators but the point is to have
  only one server thread for all possible uses.  You just
  have to translate rank_in_comm into rank_in_world by
  getting the groups from both communicators using
  MPI_Comm_group and then creating a map from ranks in
  comm to ranks in world using MPI_Group_translate_ranks.

  The class is a singleton ... i.e., there is only one instance of it
  that is made the first time that you call RMI::instance().

  Handler routines should have this type

  typedef void (*rmi_handlerT)(void* buf, size_t nbyte);

  There are few user accessible routines.

  RMI::Request RMI::isend(const void* buf, size_t nbyte, int dest,
                          rmi_handlerT func, unsigned int attr=0)
  - to send an asynchronous message
  - RMI::Request has the same interface as SafeMPI::Request
  (right now it is a SafeMPI::Request but this is not guaranteed)

  void RMI::begin()
  - to start the server thread

  void RMI::end()
  - to terminate the server thread

  bool RMI::get_debug()
  - to get the debug flag

  void RMI::set_debug(bool)
  - to set the debug flag

*/

/**
 \file worldrmi.h
 \brief Lowest level API for sending active messages --- you should probably be looking at worldam.h instead.
 \addtogroup parallel_runtime
 */

namespace madness {

    /// This is the generic low-level interface for a message handler
    typedef void (*rmi_handlerT)(void* buf, size_t nbyte);
    typedef std::ptrdiff_t rel_fn_ptr_t;

    struct qmsg {
        typedef uint16_t counterT;  //!< counter for ordered messages
        typedef uint32_t attrT;  //!< attributes of the message; high 16 bits are the counter
        size_t len;
        rmi_handlerT func;
        int i;               // buffer index
        ProcessID src;
        attrT attr;
        counterT count;

        qmsg(size_t len, rmi_handlerT func, int i, int src, attrT attr, counterT count)
            : len(len), func(func), i(i), src(src), attr(attr), count(count) {}

        // N.B. since msg counters in same batch might wrap around 0, need to sort buckets defined by the 2 highest
        // bits ... basically we want 11xxx < 00xxx < 01xxx < 10xxx < 11xxx ... assume we only have messages from
        // at most 2 adjacent buckets, thus only when comparing counters from buckets 00 and 11 reverse the order
        // P.S. we thus assume we won't have to deal with msg sequences > 2^14 (per rank)
        friend inline bool operator<(const qmsg& a, const qmsg& b) {
            const auto a_src = a.src;
            const auto b_src = b.src;
            if (a_src == b_src) {
              const auto a_count_bk = a.count >> 14;
              const auto b_count_bk = b.count >> 14;
              if (a_count_bk == 0b00 && b_count_bk == 0b11) {
                return false;
              } else if (a_count_bk == 0b11 && b_count_bk == 0b00) {
                return true;
              } else {
                return a.count < b.count;
              }
            } else {
              return a_src < b_src;
            }
          }

        qmsg() {}
    }; // struct qmsg


    // Holds message passing statistics
    struct RMIStats {
        uint64_t nmsg_sent;
        uint64_t nbyte_sent;
        uint64_t nmsg_recv;
        uint64_t nbyte_recv;
        uint64_t max_serv_send_q;

        RMIStats()
            : nmsg_sent(0), nbyte_sent(0), nmsg_recv(0), nbyte_recv(0), max_serv_send_q(0) {}
    };

    /// This for RMI server thread to manage lifetime of WorldAM messages that it is sending
    struct RMISendReq {
        virtual bool TestAndFree() = 0;
        virtual ~RMISendReq() {} // ESSENTIAL!!
    };

    /// This class implements the communications server thread and provides the only send interface
    class RMI  {
        typedef qmsg::counterT counterT;
        typedef qmsg::attrT attrT;

        /// @return reference to the boolean variable indicating whether this thread is the server thread
        static bool& is_server_thread_accessor();

    public:

        typedef SafeMPI::Request Request;

        // Choose header length to hold at least sizeof(header) and
        // also to ensure good alignment of the user payload.
        static const size_t ALIGNMENT = 64;
        static const size_t HEADER_LEN = ALIGNMENT;
        static const attrT ATTR_UNORDERED=0x0;
        static const attrT ATTR_ORDERED=0x1;

        static int testsome_backoff_us;

        static void set_this_thread_is_server(bool flag = true) { is_server_thread_accessor() = flag;}
        static bool get_this_thread_is_server() {return is_server_thread_accessor();}

        static std::list< std::unique_ptr<RMISendReq> > send_req; // List of outstanding world active messages sent by the server

    private:

        static void clear_send_req() {
            //std::cout << "clearing server messages " << pthread_self() << std::endl;
            stats.max_serv_send_q = std::max(stats.max_serv_send_q,uint64_t(send_req.size()));
            auto it=send_req.begin();
            while (it != send_req.end()) {
                if ((*it)->TestAndFree()) 
                    it = send_req.erase(it);
                else 
                    ++it;
            }
        }

        class RmiTask
#if HAVE_INTEL_TBB
                : private madness::Mutex
#else
                : public madness::ThreadBase, private madness::Mutex
#endif // HAVE_INTEL_TBB
        {

        public:

            struct header {
                rel_fn_ptr_t func;
                attrT attr;
            }; // struct header

            /// q of huge messages, each msg = {source,nbytes,tag}
            std::list< std::tuple<int,size_t,int> > hugeq;

            SafeMPI::Intracomm comm;
            const int nproc;            // No. of processes in comm world
            const ProcessID rank;       // Rank of this process
            std::atomic<bool> finished;     // True if finished ... atomic seems preferable to volatile
            std::unique_ptr<counterT[]> send_counters; // used to be volatile but no need
            std::unique_ptr<counterT[]> recv_counters;
            std::size_t max_msg_len_;
            std::size_t nrecv_;
            long nssend_;
            std::size_t maxq_;
            std::unique_ptr<void*[]> recv_buf; // Will be at least ALIGNMENT aligned ... +1 for huge messages
            std::unique_ptr<SafeMPI::Request[]> recv_req;

            std::unique_ptr<SafeMPI::Status[]> status;
            std::unique_ptr<int[]> ind;
            std::unique_ptr<qmsg[]> q;
            int n_in_q;

            static inline bool is_ordered(attrT attr) { return attr & ATTR_ORDERED; }

            void process_some();

            RmiTask(const SafeMPI::Intracomm& comm = SafeMPI::COMM_WORLD);
            virtual ~RmiTask();

            static void set_rmi_task_is_running(bool flag = true);

#if HAVE_INTEL_TBB
            void run() {
  	        ::madness::binder.bind();
                set_rmi_task_is_running(true);
                RMI::set_this_thread_is_server(true);

                while (! finished) process_some();

                RMI::set_this_thread_is_server(false);
                set_rmi_task_is_running(false);

                finished = false;  // to ensure that RmiTask::exit() that
                                   // triggered the exit proceeds to completion
            }
#else
            void run() {
  	        ::madness::binder.bind();
                RMI::set_this_thread_is_server(true);
                try {
                    while (! finished) process_some();
                    finished = false;
                } catch(...) {
                    delete this;
                    RMI::set_this_thread_is_server(false);
                    throw;
                }
                RMI::set_this_thread_is_server(false);
            }
#endif // HAVE_INTEL_TBB

            void exit() {
                if (debugging)
                  print_error(rank, ":RMI: sending exit request to server thread\n");

                // Set finished flag
                finished = true;
                while(finished)
                    myusleep(1000);
            }

            static void huge_msg_handler(void *buf, size_t nbytein);

            Request isend(const void* buf, size_t nbyte, ProcessID dest, rmi_handlerT func, attrT attr);

            void post_pending_huge_msg();

            void post_recv_buf(int i);

        private:

            /// thread-safely round-robins through tags in [first_tag, first_tag+period) range
            /// @returns new tag to be used in messaging
            int unique_tag() const;
            /// the period of tags returned by unique_tag()
            /// @warning this bounds how many huge messages each RmiTask will be able to process
            static constexpr int unique_tag_period() { return 2048; }

        }; // class RmiTask


        static std::unique_ptr<RmiTask> task_ptr;    // Pointer to the singleton instance
        static RMIStats stats;
        static bool debugging;    // True if debugging ... used to be volatile but no need

        static const size_t DEFAULT_MAX_MSG_LEN = 3*512*1024;  //!< the default size of recv buffers, in bytes; the actual size can be configured by the user via envvar MAD_BUFFER_SIZE
        static const int DEFAULT_NRECV = 128;  //!< the default # of recv buffers; the actual number can be configured by the user via envvar MAD_RECV_BUFFERS

        // Not allowed
        RMI(const RMI&);
        RMI& operator=(const RMI&);

    public:

        /// Returns the size of recv buffers, in bytes

        /// @return The size of recv buffers, in bytes
        /// @note The default value is given by RMI::DEFAULT_MAX_MSG_LEN, can be overridden at runtime by the user via environment variable MAD_BUFFER_SIZE.
        /// @warning Cannot be smaller than 1024 bytes.
        static std::size_t max_msg_len() {
            MADNESS_ASSERT(task_ptr);
            return task_ptr->max_msg_len_;
        }
        static std::size_t maxq() {
            MADNESS_ASSERT(task_ptr);
            return task_ptr->maxq_;
        }

        /// Returns the number of recv buffers

        /// @return The number of recv buffers
        /// @note The default value is given by RMI::DEFAULT_NRECV, can be overridden at runtime by the user via environment variable MAD_RECV_BUFFERS
        /// @warning Cannot be smaller than 32.
        static std::size_t nrecv() {
            MADNESS_ASSERT(task_ptr);
            return task_ptr->nrecv_;
        }

        /// Send a remote method invocation (again you should probably be looking at worldam.h instead)

        /// @param[in] buf Pointer to the data buffer (do not modify until send is completed)
        /// @param[in] nbyte Size of the data in bytes
        /// @param[in] dest Process to receive the message
        /// @param[in] func The function to handle the message on the remote end
        /// @param[in] attr Attributes of the message (ATTR_UNORDERED or ATTR_ORDERED)
        /// @return The status as an RMI::Request that presently is a SafeMPI::Request
        static Request
        isend(const void* buf, size_t nbyte, ProcessID dest, rmi_handlerT func, unsigned int attr=ATTR_UNORDERED) {
            if(!task_ptr) {
              print_error(
                  "!! MADNESS RMI error: Attempting to send a message when the RMI thread is not running\n"
                  "!! MADNESS RMI error: This typically occurs when an active message is sent or a remote task is spawned after calling madness::finalize()\n");
              MADNESS_EXCEPTION("!! MADNESS error: The RMI thread is not running", (task_ptr != nullptr));
            }
            return task_ptr->isend(buf, nbyte, dest, func, attr);
        }

        /// will complain to std::cerr and throw if ASLR is on by making
        /// sure that address of this function matches across @p comm
        /// @param[in] comm the communicator
        static void assert_aslr_off(const SafeMPI::Intracomm& comm = SafeMPI::COMM_WORLD);

        static void begin(const SafeMPI::Intracomm& comm = SafeMPI::COMM_WORLD);

        static void end() {
            if(task_ptr) {
                task_ptr->exit();
                //exit insures that RMI task is completed, therefore it is OK to delete it
                task_ptr = nullptr;
            }
        }

        static void set_debug(bool status) { debugging = status; }

        static bool get_debug() { return debugging; }

        static const RMIStats& get_stats() { return stats; }
    }; // class RMI

} // namespace madness

#endif // MADNESS_WORLD_WORLDRMI_H__INCLUDED
