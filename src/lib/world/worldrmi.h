#ifndef WORLD_RMI_H
#define WORLD_RMI_H

#include <world/safempi.h>
#include <cstdlib>
#include <iostream>
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

namespace madness {

    // This is the generic low-level interface for a message handler
    typedef void (*rmi_handlerT)(void* buf, size_t nbyte);


    class RMI : private madness::ThreadBase , madness::Mutex {
    public:
        // Choose header length to hold at least sizeof(header) and
        // also to ensure good alignment of the user payload.
        static const size_t ALIGNMENT = 64;
        static const size_t HEADER_LEN = ALIGNMENT;
        static const size_t MAX_MSG_LEN = 256*1024;

        static const unsigned int ATTR_UNORDERED=0x0;
        static const unsigned int ATTR_ORDERED=0x1;
        
        typedef SafeMPI::Request Request;

    private:
        static const int NRECV=32;
        static const int MAXQ=4*NRECV;

        SafeMPI::Intracomm comm;
        const Tag tag;              // The MPI message type 
        const int nproc;            // No. of processes in comm world
        const ProcessID rank;       // Rank of this process
        volatile bool debugging;    // True if debugging
        volatile bool finished;     // True if finished

        volatile unsigned char* send_counters;
        unsigned char* recv_counters;
        unsigned char* recv_buf[NRECV]; // Will be at least ALIGNMENT aligned
        SafeMPI::Request recv_req[NRECV];
    
        static RMI* instance_ptr;    // Pointer to the singleton instance

        static bool is_ordered(unsigned int attr) {
            return attr & ATTR_ORDERED;
        }

        struct header {
            rmi_handlerT func;
            unsigned int attr;
        };

        void run() {
            ThreadBase::set_affinity(1); // The RMI thread is logical thread 1
            if (debugging)
                std::cerr << rank << ":RMI: server thread is running" << std::endl;
            // The RMI server thread spends its life in here
            MPI::Status status[NRECV];
            int ind[NRECV];
            struct qmsg {
                size_t len;
                rmi_handlerT func;
                int i;
                ProcessID src;
                unsigned int attr;
                int count;

                qmsg(size_t len, rmi_handlerT func, int i, int src, unsigned int attr, int count) 
                    : len(len), func(func), i(i), src(src), attr(attr), count(count) {}

                qmsg() {}
            };
            qmsg q[MAXQ];
            int n_in_q = 0;
            
            while (1) {

                if (debugging && n_in_q)
                    std::cerr << rank << ":RMI: about to call Waitsome with " 
                              << n_in_q << " messages in the queue" << std::endl;

                // If MPI is not safe for simultaneous entry by multiple threads we
                // cannot call Waitsome ... have to poll via Testsome
                int narrived;
                while (!(narrived = SafeMPI::Request::Testsome(NRECV, recv_req, ind, status))) {
                    //usleep(100);
                    if (finished) return;
                }

                if (debugging)
                    std::cerr << rank << ":RMI: " << narrived 
                              << " messages just arrived" << std::endl;

                if (narrived) {
                    for (int m=0; m<narrived; m++) {
                        int src = status[m].Get_source();
                        size_t len = status[m].Get_count(MPI::BYTE);
                        int i = ind[m];
                    
                        const header* h = (const header*)(recv_buf[i]);
                        rmi_handlerT func = h->func;
                        unsigned int attr = h->attr;
                        int count = (attr>>16); //&&0xff;

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
                }

                // Only ordered messages can end up in the queue due to
                // out-of-order receipt or recv buffer processing.
                int ndone;
                do {
                    ndone = 0;
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
                            ndone++;
                            q[m].func(recv_buf[q[m].i], q[m].len);
                            post_recv_buf(q[m].i);
                        
                            // Replace msg just processed with one at end (if there)
                            n_in_q--;
                            if (m != n_in_q) {
                                q[m] = q[n_in_q];
                                m--; // Since for loop will increment it
                            }
                        }
                    }
                } while (ndone);
            }
        }

        void post_recv_buf(int i) {
            recv_req[i] = comm.Irecv(recv_buf[i], MAX_MSG_LEN, MPI::BYTE, MPI::ANY_SOURCE, tag);
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
            , tag(1023)
            , nproc(comm.Get_size())
            , rank(comm.Get_rank())
            , debugging(false)
            , finished(false)
            , send_counters(new unsigned char[nproc])
            , recv_counters(new unsigned char[nproc])
        {
            for (int i=0; i<nproc; i++) send_counters[i] = 0;
            for (int i=0; i<nproc; i++) recv_counters[i] = 0;
            for (int i=0; i<NRECV; i++) {
                if (posix_memalign((void**)(recv_buf+i), ALIGNMENT, MAX_MSG_LEN)) 
                    throw "RMI:initialize:failed allocating aligned recv buffer";
                post_recv_buf(i);
            }
            start();
        }

        static RMI* instance() {
            if (!instance_ptr) {
                instance_ptr = new RMI();
            }
            return instance_ptr;
        }
    
        Request private_isend(const void* buf, size_t nbyte, ProcessID dest, rmi_handlerT func, unsigned int attr) {
            if (debugging) 
                std::cerr << instance_ptr->rank 
                          << ":RMI: sending buf=" << buf 
                          << " nbyte=" << nbyte 
                          << " dest=" << dest 
                          << " func=" << (void*)(func) 
                          << " ordered=" << is_ordered(attr)
                          << " count=" << int(send_counters[dest])
                          << std::endl;

            if (nbyte < HEADER_LEN) 
                throw "RMI::isend --- your buffer is too small to hold the header";
        
            // If ordering need the mutex to enclose sending the message
            // otherwise there is a livelock scenario due to a starved thread
            // holding an early counter.
            if (is_ordered(attr)) {
                lock();
                attr |= ((send_counters[dest]++)<<16);
            }

            header* h = (header*)(buf);
            h->func = func;
            h->attr = attr;

            Request result = comm.Isend(buf, nbyte, MPI::BYTE, dest, tag);

            if (is_ordered(attr)) unlock();

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
    };
}
        
#endif
