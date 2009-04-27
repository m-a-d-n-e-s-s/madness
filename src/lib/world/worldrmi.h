#ifndef WORLD_RMI_H
#define WORLD_RMI_H

#include <world/safempi.h>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <world/worldthread.h>
#include <world/worldtypes.h>

/*
  There is just one server thread and it is the only one
  messing with the recv buffers, so there is no need for
  mutex on recv related data.

  Multiple threads (including the server) may send hence
  we need to be careful about send-related data.

  When MPI is initialized we need to use init_thread with
  multiple required.

  This RMI service operates only in COMM_WORLD.  It easy enough
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
  - RMI::Request has the same interface as MPI::Request
  (right now it is an MPI::Request but this is not guaranteed)

  void RMI::begin()
  - to start the server thread

  void RMI::end()
  - to terminate the server thread

  bool RMI::get_debug()
  - to get the debug flag

  void RMI::set_debug(bool)
  - to set the debug flag

*/

#if ON_A_MAC
#include <sys/errno.h>
static inline int posix_memalign(void **memptr, std::size_t alignment, std::size_t size) {
    *memptr=malloc(size);
    if (*memptr) return 0;
    else return ENOMEM;
}
#elif MISSING_POSIX_MEMALIGN_PROTO
extern "C"  int posix_memalign(void **memptr, std::size_t alignment, std::size_t size);
#endif

namespace madness {

    // This is the generic low-level interface for a message handler
    typedef void (*rmi_handlerT)(void* buf, size_t nbyte);

    struct qmsg {
        typedef uint16_t counterT;
        typedef uint32_t attrT;
        size_t len;
        rmi_handlerT func; 
        int i;               // buffer index
        ProcessID src;
        attrT attr;
        counterT count;
        
        qmsg(size_t len, rmi_handlerT func, int i, int src, attrT attr, counterT count)
            : len(len), func(func), i(i), src(src), attr(attr), count(count) {}
        
        bool operator<(const qmsg& other) const {
            return count < other.count;
        }
        
        qmsg() {}
    };


    // Holds message passing statistics
    struct RMIStats {
        uint64_t nmsg_sent;
        uint64_t nbyte_sent;
        uint64_t nmsg_recv;
        uint64_t nbyte_recv;

        RMIStats()
                : nmsg_sent(0), nbyte_sent(0), nmsg_recv(0), nbyte_recv(0) {}
    };

    class RMI : private madness::ThreadBase , madness::Mutex {
        typedef uint16_t counterT;
        typedef uint32_t attrT;
    public:
        // Choose header length to hold at least sizeof(header) and
        // also to ensure good alignment of the user payload.
        static const size_t ALIGNMENT = 64;
        static const size_t HEADER_LEN = ALIGNMENT;
        static const size_t MAX_MSG_LEN = 256*1024;
        //static const size_t MAX_MSG_LEN = 1024*1024;

        static const attrT ATTR_UNORDERED=0x0;
        static const attrT ATTR_ORDERED=0x1;

        typedef SafeMPI::Request Request;

    private:
#ifdef HAVE_CRAYXT
        static const int NRECV=256;
#else
        static const int NRECV=32;
#endif
        static const int MAXQ=NRECV+1;

        std::list< std::pair<int,size_t> > hugeq; // q for incoming huge messages

        RMIStats stats;
        SafeMPI::Intracomm comm;
        const int nproc;            // No. of processes in comm world
        const ProcessID rank;       // Rank of this process
        volatile bool debugging;    // True if debugging
        volatile bool finished;     // True if finished

        volatile counterT* send_counters;
        counterT* recv_counters;
        unsigned char* recv_buf[NRECV+1]; // Will be at least ALIGNMENT aligned ... +1 for huge messages
        SafeMPI::Request recv_req[NRECV+1];

        static RMI* instance_ptr;    // Pointer to the singleton instance

        static bool is_ordered(attrT attr) {
            return attr & ATTR_ORDERED;
        }

        struct header {
            rmi_handlerT func;
            attrT attr;
        };

        void run() {
            ThreadBase::set_affinity(1); // The RMI thread is logical thread 1
            if (debugging)
                std::cerr << rank << ":RMI: server thread is running" << std::endl;
            // The RMI server thread spends its life in here
            MPI::Status status[NRECV+1];
            int ind[NRECV+1];
            qmsg q[MAXQ];
            int n_in_q = 0;

            while (1) {

                if (debugging && n_in_q)
                    std::cerr << rank << ":RMI: about to call Waitsome with "
                              << n_in_q << " messages in the queue" << std::endl;

                // If MPI is not safe for simultaneous entry by multiple threads we
                // cannot call Waitsome ... have to poll via Testsome
                int narrived;

                MutexWaiter waiter;
                while (!(narrived = SafeMPI::Request::Testsome(NRECV+1, recv_req, ind, status))) {
                    if (finished) return;
#ifdef HAVE_CRAYXT
                    myusleep(5);
#else
                    waiter.wait();
#endif
                }

#ifndef HAVE_CRAYXT
                waiter.reset();
#endif

                if (debugging)
                    std::cerr << rank << ":RMI: " << narrived
                              << " messages just arrived" << std::endl;

                if (narrived) {
                    for (int m=0; m<narrived; m++) {
                        int src = status[m].Get_source();
                        size_t len = status[m].Get_count(MPI::BYTE);
                        int i = ind[m]; 

                        stats.nmsg_recv++;
                        stats.nbyte_recv += len;

                        const header* h = (const header*)(recv_buf[i]);
                        rmi_handlerT func = h->func;
                        attrT attr = h->attr;
                        counterT count = (attr>>16); //&&0xffff;

                        if (!is_ordered(attr) || count==recv_counters[src]) {
                            // Unordered and in order messages should be digested as soon as possible.
                            if (debugging)
                                std::cerr << rank
                                          << ":RMI: invoking from=" << src
                                          << " nbyte=" << len
                                          << " func=" << (void*)(func)
                                          << " ordered=" << is_ordered(attr)
                                          << " count=" << count
                                          << std::endl;

                            if (is_ordered(attr)) recv_counters[src]++;
                            func(recv_buf[i], len);
                            post_recv_buf(i);
                        }
                        else {
                            if (debugging)
                                std::cerr << rank
                                          << ":RMI: enqueing from=" << src
                                          << " nbyte=" << len
                                          << " func=" << (void*)(func)
                                          << " ordered=" << is_ordered(attr)
                                          << " fromcount=" << count
                                          << " herecount=" << int(recv_counters[src])
                                          << std::endl;
                            // Shove it in the queue
                            int n = n_in_q++;
                            if (n >= MAXQ) throw "RMI:server: overflowed out-of-order message q\n";
                            q[n] = qmsg(len, func, i, src, attr, count);
                        }
                    }

                    // Only ordered messages can end up in the queue due to
                    // out-of-order receipt or order of recv buffer processing.
                    
                    // Sort queued messages by ascending recv count
                    std::sort(&q[0],&q[0]+n_in_q);
                    
                    // Loop thru messages ... since we have sorted only one pass
                    // is necessary and if we cannot process a message we
                    // save it at the beginning of the queue
                    int nleftover = 0;
                    for (int m=0; m<n_in_q; m++) {
                        int src = q[m].src;
                        if (q[m].count == recv_counters[src]) {
                            if (debugging)
                                std::cerr << rank
                                          << ":RMI: queue invoking from=" << src
                                          << " nbyte=" << q[m].len
                                          << " func=" << (void*)(q[m].func)
                                          << " ordered=" << is_ordered(q[m].attr)
                                          << " count=" << q[m].count
                                          << std::endl;
                            
                            recv_counters[src]++;
                            q[m].func(recv_buf[q[m].i], q[m].len);
                            post_recv_buf(q[m].i);
                        }
                        else {
                            q[nleftover++] = q[m];
                            if (debugging)
                                std::cerr << rank
                                          << ":RMI: queue pending out of order from=" << src
                                          << " nbyte=" << q[m].len
                                          << " func=" << (void*)(q[m].func)
                                          << " ordered=" << is_ordered(q[m].attr)
                                          << " count=" << q[m].count
                                          << std::endl;
                        }
                    }
                    n_in_q = nleftover;

                    post_pending_huge_msg();
                }
            }
        }

        void post_pending_huge_msg() {
            if (recv_buf[NRECV]) return;      // Message already pending
            if (!hugeq.empty()) {
                int src = hugeq.front().first;
                size_t nbyte = hugeq.front().second;
                hugeq.pop_front();
                if (posix_memalign((void **)(recv_buf+NRECV), ALIGNMENT, nbyte))
                    throw "RMI: failed allocating huge message";
                recv_req[NRECV] = comm.Irecv(recv_buf[NRECV], nbyte, MPI::BYTE, src, SafeMPI::RMI_HUGE_DAT_TAG);
                int nada=0;
                comm.Send(&nada, sizeof(nada), MPI::BYTE, src, SafeMPI::RMI_HUGE_ACK_TAG);
            }
        }

        void post_recv_buf(int i) {
            if (i < NRECV) {
                recv_req[i] = comm.Irecv(recv_buf[i], MAX_MSG_LEN, MPI::BYTE, MPI::ANY_SOURCE, SafeMPI::RMI_TAG);
            }
            else if (i == NRECV) {
                free(recv_buf[i]);
                recv_buf[i] = 0;
                post_pending_huge_msg();
            }
            else {
                throw "RMI::post_recv_buf: confusion";
            }
        }

        ~RMI() {
            //delete send_counters;
            //delete recv_counters;
            //         if (!MPI::Is_finalized()) {
            //             for (int i=0; i<NRECV; i++) {
            //                 if (!recv_req[i].Test())
            //                     recv_req[i].Cancel();
            //             }
            //         }
            //for (int i=0; i<NRECV; i++) free(recv_buf[i]);
        }

        RMI()
                : comm(MPI::COMM_WORLD)
                , nproc(comm.Get_size())
                , rank(comm.Get_rank())
                , debugging(false)
                , finished(false)
                , send_counters(new unsigned short[nproc])
                , recv_counters(new unsigned short[nproc]) {
            for (int i=0; i<nproc; i++) send_counters[i] = 0;
            for (int i=0; i<nproc; i++) recv_counters[i] = 0;
            if (nproc > 1) {
                for (int i=0; i<NRECV; i++) {
                    if (posix_memalign((void**)(recv_buf+i), ALIGNMENT, MAX_MSG_LEN))
                        throw "RMI:initialize:failed allocating aligned recv buffer";
                    post_recv_buf(i);
                }
                recv_buf[NRECV] = 0;
                start();
            }
        }

        static RMI* instance() {
            if (!instance_ptr) {
                instance_ptr = new RMI();
            }
            return instance_ptr;
        }

        static void huge_msg_handler(void *buf, size_t nbytein) {
            const size_t* info = (size_t *)(buf);
            int nword = HEADER_LEN/sizeof(size_t);
            int src = info[nword];
            size_t nbyte = info[nword+1];

            instance()->hugeq.push_back(std::make_pair(src,nbyte));
            instance()->post_pending_huge_msg();
        }

        Request private_isend(const void* buf, size_t nbyte, ProcessID dest, rmi_handlerT func, attrT attr) {
            int tag = SafeMPI::RMI_TAG;

            if (nbyte > MAX_MSG_LEN) {
                // Huge message protocol ... send message to dest indicating size and origin of huge message.
                // Remote end posts a buffer then acks the request.  This end can then send.
                int nword = HEADER_LEN/sizeof(size_t);
                size_t info[nword+2];
                info[nword  ] = rank;
                info[nword+1] = nbyte;

                int ack;
                Request req_ack = comm.Irecv(&ack, sizeof(ack), MPI::BYTE, dest, SafeMPI::RMI_HUGE_ACK_TAG);
                Request req_send = private_isend(info, sizeof(info), dest, RMI::huge_msg_handler, ATTR_UNORDERED);

                MutexWaiter waiter;
                while (!req_send.Test()) waiter.wait();
                waiter.reset();
                while (!req_ack.Test()) waiter.wait();

                tag = SafeMPI::RMI_HUGE_DAT_TAG;
            }
            else if (nbyte < HEADER_LEN) {
                throw "RMI::isend --- your buffer is too small to hold the header";
            }

            if (debugging)
                std::cerr << instance_ptr->rank
                          << ":RMI: sending buf=" << buf
                          << " nbyte=" << nbyte
                          << " dest=" << dest
                          << " func=" << (void*)(func)
                          << " ordered=" << is_ordered(attr)
                          << " count=" << int(send_counters[dest])
                          << std::endl;

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
            h->func = func;
            h->attr = attr;

            stats.nmsg_sent++;
            stats.nbyte_sent += nbyte;

            Request result = comm.Isend(buf, nbyte, MPI::BYTE, dest, tag);

            //if (is_ordered(attr)) unlock();
            unlock();

            return result;
        }

        void private_exit() {
            if (debugging)
                std::cerr << instance_ptr->rank << ":RMI: sending exit request to server thread" << std::endl;

            finished = true;
            usleep(10000);

            //delete this;
        }

    public:
        static Request isend(const void* buf, size_t nbyte, ProcessID dest, rmi_handlerT func, unsigned int attr=ATTR_UNORDERED) {
            return instance()->private_isend(buf, nbyte, dest, func, attr);
        }

        static void end() {
            if (instance_ptr) instance_ptr->private_exit();
        }

        static void begin() {
            instance();
        }

        static void set_debug(bool status) {
            instance()->debugging = status;
        }

        static bool get_debug() {
            return instance()->debugging;
        }

        static const RMIStats& get_stats() {
            return instance()->stats;
        }
    };
}

#endif
