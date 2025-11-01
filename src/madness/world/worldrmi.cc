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

  $Id$
*/

#include <madness/world/worldrmi.h>
#include <madness/world/posixmem.h>
#include <madness/world/timers.h>
#include <madness/world/units.h>
#include <iostream>
#include <algorithm>
#include <utility>
#include <sstream>
#include <list>
#include <memory>
#include <atomic>
#include <madness/world/safempi.h>
#include <madness/world/archive.h>

namespace madness {

    std::unique_ptr<RMI::RmiTask> RMI::task_ptr = nullptr;
    RMIStats RMI::stats;
    bool RMI::debugging = false;
    std::list< std::unique_ptr<RMISendReq> > RMI::send_req;

    bool& RMI::is_server_thread_accessor() {
      static thread_local bool is_server_thread = false;
      return is_server_thread;
    }


    void RMI::RmiTask::process_some() {

        const bool print_debug_info = RMI::debugging;

        if (print_debug_info && n_in_q)
            print_error(rank, ":RMI: about to call Testsome with ", n_in_q, " messages in the queue\n");

        // If MPI is not safe for simultaneous entry by multiple threads we
        // cannot call Waitsome ... have to poll via Testsome

        // Now that the server thread doing other stuff (including being
        // responsible for its own outbound messages) we have to poll.
        int narrived = 0, iterations = 0;

        MutexWaiter waiter;
        while((narrived == 0) && (iterations < 1000)) {
          narrived = SafeMPI::Request::Testsome(maxq_, recv_req.get(), ind.get(), status.get());
          if (narrived) break;
          ++iterations;
          clear_send_req();
          myusleep(RMI::testsome_backoff_us);
        }

#ifndef HAVE_CRAYXT
        waiter.reset();
#endif

        if (print_debug_info && narrived > 0)
            print_error(rank, ":RMI: ", narrived, " messages just arrived\n");

        if (narrived) {
            for (int m=0; m<narrived; ++m) {
                const int src = status[m].Get_source();
                const size_t len = status[m].Get_count(MPI_BYTE);
                const int i = ind[m];

                ++(RMI::stats.nmsg_recv);
                RMI::stats.nbyte_recv += len;

                const header* h = (const header*)(recv_buf[i]);
                rmi_handlerT func = archive::to_abs_fn_ptr<rmi_handlerT>(h->func);
                const attrT attr = h->attr;
                const counterT count = (attr>>16); //&&0xffff;

                if (!is_ordered(attr) || count==recv_counters[src]) {
                    // Unordered and in order messages should be digested as soon as possible.
                    if (print_debug_info)
                      print_error(rank, ":RMI: invoking from=", src,
                                  " nbyte=", len, " func=", func,
                                  " ordered=", is_ordered(attr),
                                  " count=", count, "\n");

                    if (is_ordered(attr)) ++(recv_counters[src]);
                    func(recv_buf[i], len);
                    post_recv_buf(i);
                }
                else {
                  if (print_debug_info)
                    print_error(rank, ":RMI: enqueing from=", src,
                                " nbyte=", len, " func=", func,
                                " ordered=", is_ordered(attr),
                                " fromcount=", count,
                                " herecount=", int(recv_counters[src]), "\n");
                  // Shove it in the queue
                  const int n = n_in_q++;
                  if (n >= (int)maxq_)
                    MADNESS_EXCEPTION(
                        "RMI:server: overflowed out-of-order message q\n", n);
                  q[n] = qmsg(len, func, i, src, attr, count);
                }
            }

            // Only ordered messages can end up in the queue due to
            // out-of-order receipt or order of recv buffer processing.

            // Sort queued messages by source and ascending (modulo overflow) recv count
            std::sort(q.get(),q.get()+n_in_q);

            // Loop thru messages ... since we have sorted only one pass
            // is necessary and if we cannot process a message we
            // save it at the beginning of the queue
            int nleftover = 0;
            for (int m=0; m<n_in_q; ++m) {
                const int src = q[m].src;
                if (q[m].count == recv_counters[src]) {
                  if (print_debug_info)
                    print_error(rank, ":RMI: queue invoking from=", src,
                                " nbyte=", q[m].len, " func=", q[m].func,
                                " ordered=", is_ordered(q[m].attr),
                                " count=", q[m].count, "\n");

                  ++(recv_counters[src]);
                  q[m].func(recv_buf[q[m].i], q[m].len);
                  post_recv_buf(q[m].i);
                }
                else {
                    q[nleftover++] = q[m];
                    if (print_debug_info)
                      print_error(rank,
                                  ":RMI: queue pending out of order from=", src,
                                  " nbyte=", q[m].len, " func=", q[m].func,
                                  " ordered=", is_ordered(q[m].attr),
                                  " count=", q[m].count, "\n");
                }
            }
            n_in_q = nleftover;

            post_pending_huge_msg();

#if !(defined(HAVE_INTEL_TBB) || defined(HAVE_PARSEC))
            // Since this thread never waits or tries to run a task we
            // must manually flush the thread-local prebuffer that
            // aggregates task submission.
            ThreadPool::instance()->flush_prebuf();
#endif
            clear_send_req();
        }
    }

    void RMI::RmiTask::post_pending_huge_msg() {
        if (recv_buf[nrecv_]) return;      // Message already pending
        if (!hugeq.empty()) {
            const auto& hugemsg = hugeq.front();
            const int src = std::get<0>(hugemsg);
            const size_t nbyte = std::get<1>(hugemsg);
            const int tag = std::get<2>(hugemsg);
            hugeq.pop_front();
            if (posix_memalign(&recv_buf[nrecv_], ALIGNMENT, nbyte))
                MADNESS_EXCEPTION("RMI: failed allocating huge message", 1);
            recv_req[nrecv_] = comm.Irecv(recv_buf[nrecv_], nbyte, MPI_BYTE, src, tag);
            int nada=0;
            // make unique tags to ensure that ack msgs do not collide with normal recv msgs
#ifdef MADNESS_USE_BSEND_ACKS
            comm.Bsend(&nada, sizeof(nada), MPI_BYTE, src, tag + unique_tag_period());
#else
            comm.Send(&nada, sizeof(nada), MPI_BYTE, src, tag + unique_tag_period());
#endif // MADNESS_USE_BSEND_ACKS
        }
    }

    void RMI::RmiTask::post_recv_buf(int i) {
        if (i < (int)nrecv_) {
            recv_req[i] = comm.Irecv(recv_buf[i], max_msg_len_, MPI_BYTE, MPI_ANY_SOURCE, SafeMPI::RMI_TAG);
        }
        else if (i == (int)nrecv_) {
            free(recv_buf[i]);
            recv_buf[i] = 0;
            post_pending_huge_msg();
        }
        else {
            MADNESS_EXCEPTION("RMI::post_recv_buf: confusion", i);
        }
    }

    RMI::RmiTask::~RmiTask() {
        //         if (!SafeMPI::Is_finalized()) {
        //             for (int i=0; i<nrecv_; ++i) {
        //                 if (!recv_req[i].Test())
        //                     recv_req[i].Cancel();
        //             }
        //         }
        //for (int i=0; i<nrecv_; ++i) free(recv_buf[i]);
    }

    static std::atomic<bool> rmi_task_is_running = false;

    RMI::RmiTask::RmiTask(const SafeMPI::Intracomm& _comm)
            : comm(_comm.Clone())
            , nproc(comm.Get_size())
            , rank(comm.Get_rank())
            , finished(false)
            , send_counters(new counterT[nproc])
            , recv_counters(new counterT[nproc])
            , max_msg_len_(DEFAULT_MAX_MSG_LEN)
            , nrecv_(DEFAULT_NRECV)
            , maxq_(DEFAULT_NRECV + 1)
            , recv_buf()
            , recv_req()
            , status()
            , ind()
            , q()
            , n_in_q(0)
    {
        // Get the maximum buffer size from the MAD_BUFFER_SIZE environment
        // variable.
        const char* mad_buffer_size = getenv("MAD_BUFFER_SIZE");
        if(mad_buffer_size) {
            // Convert the string into bytes
            max_msg_len_ = cstr_to_memory_size(mad_buffer_size);
            // Check that the size of the receive buffers is reasonable.
            if(max_msg_len_ < 1024) {
                max_msg_len_ = DEFAULT_MAX_MSG_LEN; // = 3*512*1024
                print_error(
                    "!!! WARNING: MAD_BUFFER_SIZE must be at least 1024 "
                    "bytes.\n",
                    "!!! WARNING: Increasing MAD_BUFFER_SIZE to the default "
                    "size, ",
                    max_msg_len_, " bytes.\n");
            }
            // Check that the buffer has the correct alignment
            const std::size_t unaligned = max_msg_len_ % ALIGNMENT;
            if(unaligned != 0)
                max_msg_len_ += ALIGNMENT - unaligned;
        }

        // Get the number of receive buffers from the MAD_RECV_BUFFERS
        // environment variable.
        const char* mad_recv_buffs = getenv("MAD_RECV_BUFFERS");
        if(mad_recv_buffs) {
            std::stringstream ss(mad_recv_buffs);
            ss >> nrecv_;
            // Check that the number of receive buffers is reasonable.
            if(nrecv_ < 32) {
                nrecv_ = DEFAULT_NRECV;
                print_error(
                    "!!! WARNING: MAD_RECV_BUFFERS must be at least 32.\n",
                    "!!! WARNING: Increasing MAD_RECV_BUFFERS to ", nrecv_,
                    ".\n");
            }
            maxq_ = nrecv_ + 1;
        }

        // Get environment variable controlling use of synchronous send (MAD_NSSEND)
        // negative=sends synchronous message every MAD_RECV_BUFFER sends (default)
        //        0=never send synchronous message
        //      n>0=sends synchronous message every n sends (n=1 always uses ssend)
        nssend_ = nrecv_;
        const char* mad_nssend = getenv("MAD_NSSEND");
        if (mad_nssend) {
            std::stringstream ss(mad_nssend);
            ss >> nssend_;
            if (nssend_ < 0) {
                nssend_ = nrecv_;
            }
        }

        // Allocate memory for receive buffer and requests
        recv_buf.reset(new void*[maxq_]);
        recv_req.reset(new Request[maxq_]);

        // Initialize the send/recv counts
        std::fill_n(send_counters.get(), nproc, 0);
        std::fill_n(recv_counters.get(), nproc, 0);

        // Allocate buffers for message tracking
        status.reset(new SafeMPI::Status[maxq_]);
        ind.reset(new int[maxq_]);
        q.reset(new qmsg[maxq_]);
        MADNESS_ASSERT(maxq_ <= 1<<14);  // 16 bit task counter is sufficient to ensure that up to 2^14 tasks can be
                                         // pending per rank .. although maxq_ controls the TOTAL queue size, do this
                                         // purely as a reminder

        // Allocate receive buffers
        if(nproc > 1) {
            for(int i = 0; i < (int)nrecv_; ++i) {
                if(posix_memalign(&recv_buf[i], ALIGNMENT, max_msg_len_))
                    MADNESS_EXCEPTION("RMI:initialize:failed allocating aligned recv buffer", 1);
                post_recv_buf(i);
            }
            recv_buf[nrecv_] = 0;
        }
    }


    void RMI::RmiTask::huge_msg_handler(void *buf, size_t /*nbytein*/) {
        const size_t* info = (size_t *)(buf);
        int nword = HEADER_LEN/sizeof(size_t);
        const int src = info[nword];
        const size_t nbyte = info[nword+1];
        const int tag = info[nword+2];

        // extra dose of paranoia: assert that we never process so many huge messages
        // that the tag wraparound somewhere becomes possible ...
        // the worst case is where only one node sends huge messages to every node in the communicator
        // AND it has enough threads to use up all tags
        // NB list::size() is O(1) in c++11, but O(N) in older libstdc++
        bool OK = (ThreadPool::size() < size_t(RMI::RmiTask::unique_tag_period()) ||
                   RMI::task_ptr->hugeq.size() <
                   std::size_t(RMI::RmiTask::unique_tag_period() / RMI::task_ptr->comm.Get_size()));
        if (!OK) MADNESS_EXCEPTION("huge_msg_handler paranoid test failing", RMI::RmiTask::unique_tag_period());
        RMI::task_ptr->hugeq.push_back(std::make_tuple(src, nbyte, tag));
        RMI::task_ptr->post_pending_huge_msg();
    }

    namespace detail {
    void compare_fn_addresses(void* addresses_in, void* addresses_inout,
                              int* len, MPI_Datatype* type) {
      MADNESS_ASSERT(*type == MPI_LONG);
      long* in = static_cast<long*>(addresses_in);
      long* inout = static_cast<long*>(addresses_inout);
      int n = *len;
      // produce zero if addresses do not match; zero address trumps everything else
      for(int i=0; i!=n; ++i) {
        if (in[i] == 0 || inout[i] == 0 || in[i] != inout[i]) inout[i] = 0;
      }
    }
    }  // namespace detail

    void RMI::assert_aslr_off(const SafeMPI::Intracomm& comm) {
      static_assert(sizeof(long) >= sizeof(std::ptrdiff_t), "std::ptrdiff_t must not exceed the width of long");
      long my_address = reinterpret_cast<long>(&assert_aslr_off);
      MPI_Op compare_fn_addresses_op = MPI_OP_NULL;
      MPI_Op_create(reinterpret_cast<MPI_User_function*>(&detail::compare_fn_addresses), 1, &compare_fn_addresses_op);
      long zero_if_addresses_differ = 0;
      comm.Reduce(&my_address, &zero_if_addresses_differ, 1, MPI_LONG, compare_fn_addresses_op, 0);
      if (comm.Get_rank() == 0) {
        if (zero_if_addresses_differ == 0) {
          MADNESS_EXCEPTION("Address Space Layout Randomization (ASLR) detected and MADNESS could not work around it, possibly due to the use of Microsoft C++ ABI by this platform (please contact the developers to request such support); please turn off ASLR, disable by providing appropriate linker flags (see MADNESS_DISABLEPIE_LINKER_FLAG), or compile your application statically",0);
        }
        MADNESS_ASSERT(zero_if_addresses_differ == my_address);
      }
      SafeMPI::Op_free(compare_fn_addresses_op);
    }

    void RMI::begin(const SafeMPI::Intracomm& comm) {

      // if MADNESS was told to assume ASLR is disabled *completely*, check ... RMI can work with ALSR on, but that requires loading all RMI-using code as one giant shared library
#ifdef MADNESS_ASSUMES_ASLR_DISABLED
            // complain loudly and throw if ASLR is on ... RMI requires ASLR to be off
            assert_aslr_off(comm);
#endif

            testsome_backoff_us = 5;
            const char* buf = getenv("MAD_BACKOFF_US");
            if (buf) {
                std::stringstream ss(buf);
                ss >> testsome_backoff_us;
                if (testsome_backoff_us < 0) testsome_backoff_us = 0;
                if (testsome_backoff_us > 100) testsome_backoff_us = 100;
            }

            MADNESS_ASSERT(task_ptr == nullptr);
            task_ptr.reset(new RmiTask(comm));

#if HAVE_INTEL_TBB
            ThreadPool::tbb_arena->enqueue([]{
                task_ptr->run();
            });

            //TODO: is it needed ?
            task_ptr->comm.Barrier();

            while (!rmi_task_is_running) {
              myusleep(100000);
            }
#else
            task_ptr->start();
#endif // HAVE_INTEL_TBB
        }


    void RMI::RmiTask::set_rmi_task_is_running(bool flag) {
        rmi_task_is_running = flag; // Yipeeeeeeeeeeeeeeeeeeeeee ... fighting TBB laziness
    }

    RMI::Request
    RMI::RmiTask::RmiTask::isend(const void* buf, size_t nbyte, ProcessID dest, rmi_handlerT func, attrT attr) {

        MADNESS_ASSERT(nbyte <= std::numeric_limits<int>::max());

        int tag = SafeMPI::RMI_TAG;
        static std::size_t numsent = 0; // for tracking synchronous sends

        if (nbyte > max_msg_len_) {
            // Huge message protocol ... send message to dest indicating size and origin of huge message.
            // Remote end posts a buffer then acks the request.  This end can then send.
            const int nword = HEADER_LEN/sizeof(size_t);
            size_t info[nword+3];
            info[nword  ] = rank;
            info[nword+1] = nbyte;
            tag = unique_tag();
            info[nword+2] = tag;

            int ack;
            // make unique tags to ensure that ack msgs do not collide with normal recv msgs
            Request req_ack = comm.Irecv(&ack, sizeof(ack), MPI_BYTE, dest, tag + unique_tag_period());
            Request req_send = isend(info, sizeof(info), dest, RMI::RmiTask::huge_msg_handler, ATTR_UNORDERED);

            MutexWaiter waiter;
            while (!req_send.Test()) waiter.wait();
            waiter.reset();
            while (!req_ack.Test()) waiter.wait();
        }
        else if (nbyte < HEADER_LEN) {
            MADNESS_EXCEPTION("RMI::isend --- your buffer is too small to hold the header", static_cast<int>(nbyte));
        }

        if (RMI::debugging)
          print_error(rank, ":RMI: sending buf=", buf, " nbyte=", nbyte,
                      " dest=", dest, " func=", func,
                      " ordered=", is_ordered(attr),
                      " count=", int(send_counters[dest]), "\n");

        // Since most uses are ordered and we need the mutex to accumulate stats
        // we presently always get the lock
        lock();

        // If ordering need the mutex to enclose sending the message
        // otherwise there is a livelock scenario due to a starved thread
        // holding an early counter.
        if (is_ordered(attr)) {
            //lock();
            attr |= ((send_counters[dest]++)<<16);
        }

        header* h = (header*)(buf);
        h->func = archive::to_rel_fn_ptr(func);
        h->attr = attr;

        ++(RMI::stats.nmsg_sent);
        RMI::stats.nbyte_sent += nbyte;


        numsent++;
        Request result;
        if (nssend_ && numsent==std::size_t(nssend_)) {
            result = comm.Issend(buf, nbyte, MPI_BYTE, dest, tag);
            numsent %= nssend_;
        }
        else {
            result = comm.Isend(buf, nbyte, MPI_BYTE, dest, tag);
        }

        unlock();

        return result;
    }

    int RMI::RmiTask::unique_tag() const {
        constexpr int first_tag = 4096;
        static int tag = first_tag;
        /// should be able to do this with atomics
        lock();
        tag = (tag == first_tag+unique_tag_period()-1) ? first_tag : tag + 1;
        const int result = tag;
        unlock();
        return result;
    }

  int RMI::testsome_backoff_us = 2;

} // namespace madness
