/*
  This file is part of MADNESS.
  
  Copyright (C) <2007> <Oak Ridge National Laboratory>
  
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

  
#ifndef WORLDAM_H
#define WORLDAM_H

/// \file worldam.h
/// \brief Implements active message layer for World

namespace madness {
    
    typedef unsigned long ulong;

    class AmArg;
    
    /// Type of AM handler functions for short messages
    typedef void (*am_handlerT)(World&, ProcessID, const AmArg&);
    
    /// Type of AM handler functions for long messages
    typedef void (*am_long_handlerT)(World&, ProcessID, void *buf, size_t nbyte);
    

    /// Holds arguments sent via short active messages

    /// Total size is 128 bytes on a 64-bit machine and accords with
    /// gasnet minimum size for 64-bit architectures.
    class AmArg {
        friend class WorldAmInterface;
    private:
        ulong flags;   ///< System managed state 
        mutable am_handlerT function; ///< System managed handler
    public:
        ulong buf[14]; ///< System dependent length

        ulong& operator[](int i) {
            return buf[i];
        };

        ulong operator[](int i) const {
            return buf[i];
        };

        
        AmArg() : flags(0), function(0) {};
        AmArg(ulong arg0) 
            : flags(0), function(0) {
            buf[0]=arg0;
        };
        AmArg(ulong arg0, ulong arg1) 
            : flags(0), function(0) {
            buf[0]=arg0; buf[1]=arg1;
        };
        AmArg(ulong arg0, ulong arg1, ulong arg2) 
            : flags(0), function(0) {
            buf[0]=arg0; buf[1]=arg1; buf[2]=arg2;
        };
        AmArg(ulong arg0, ulong arg1, ulong arg2, ulong arg3) 
            : flags(0), function(0) {
            buf[0]=arg0; buf[1]=arg1; buf[2]=arg2; buf[3]=arg3;
        };
        AmArg(ulong arg0, ulong arg1, ulong arg2, ulong arg3, ulong arg4) 
            : flags(0), function(0) {
            buf[0]=arg0; buf[1]=arg1; buf[2]=arg2; buf[3]=arg3; buf[4]=arg4;
        };
        AmArg(ulong arg0, ulong arg1, ulong arg2, ulong arg3, ulong arg4, ulong arg5) 
            : flags(0), function(0) {
            buf[0]=arg0; buf[1]=arg1; buf[2]=arg2; buf[3]=arg3; buf[4]=arg4; buf[5]=arg5;
        };
        AmArg(ulong arg0, ulong arg1, ulong arg2, ulong arg3, ulong arg4, ulong arg5, ulong arg6) 
            : flags(0), function(0) {
            buf[0]=arg0; buf[1]=arg1; buf[2]=arg2; buf[3]=arg3; buf[4]=arg4; buf[5]=arg5; buf[6]=arg6; 
        };
        AmArg(ulong arg0, ulong arg1, ulong arg2, ulong arg3, ulong arg4, ulong arg5, ulong arg6, ulong arg7)  
            : flags(0), function(0) {
            buf[0]=arg0; buf[1]=arg1; buf[2]=arg2; buf[3]=arg3; buf[4]=arg4; buf[5]=arg5; buf[6]=arg6; buf[7]=arg7; 
        };

        template <typename A>
        AmArg(const A& a) : flags(0), function(0) {
            this->stuff(a);
        }

        template <typename A, typename B>
        AmArg(const A& a, const B& b) : flags(0), function(0) {
            this->stuff(a,b);
        }

        template <typename A, typename B, typename C>
        AmArg(const A& a, const B& b, const C& c) : flags(0), function(0) {
            this->stuff(a,b,c);
        }

        template <typename A, typename B, typename C, typename D>
        AmArg(const A& a, const B& b, const C& c, const D& d) : flags(0), function(0) {
            this->stuff(a,b,c,d);
        }

        template <typename A, typename B, typename C, typename D, typename E>
        AmArg(const A& a, const B& b, const C& c, const D& d, const E& e) : flags(0), function(0) {
            this->stuff(a,b,c,d,e);
        }

        template <typename T>
        inline operator T () const {
            T t;
            this->unstuff(t);
            return t;
        }

        template <typename A>
        inline void stuff(const A& a) {
            BufferOutputArchive ar(buf,sizeof(buf));
            ar & a;
        }


        /// Convenience template for serializing arguments into AmArg
        template <typename A, typename B>
        inline void stuff(const A& a, const B& b) {
            BufferOutputArchive ar(buf,sizeof(buf));
            ar & a & b;
        }


        /// Convenience template for serializing arguments into AmArg
        template <typename A, typename B, typename C>
        inline void stuff(const A& a, const B& b, const C& c) {
            BufferOutputArchive ar(buf,sizeof(buf));
            ar & a & b & c;
        }
        
        
        /// Convenience template for serializing arguments into AmArg
        template <typename A, typename B, typename C, typename D>
        inline void stuff(const A& a, const B& b, const C& c, const D& d) {
            BufferOutputArchive ar(buf,sizeof(buf));
            ar & a & b & c & d;
        }
        

        /// Convenience template for serializing arguments into AmArg
        template <typename A, typename B, typename C, typename D, typename E>
        inline void stuff(const A& a, const B& b, const C& c, const D& d, const E& e) {
            BufferOutputArchive ar(buf,sizeof(buf));
            ar & a & b & c & d & e;
        }
        
        
        /// Convenience template for deserializing arguments from AmArg
        template <typename A>
        inline void unstuff(A& a) const {
            BufferInputArchive ar(buf,sizeof(buf));
            ar & a;
        }
        
        
        /// Convenience template for deserializing arguments from AmArg
        template <typename A, typename B>
        inline void unstuff(A& a, B& b) const {
            BufferInputArchive ar(buf,sizeof(buf));
            ar & a & b;
        }
        
        
        /// Convenience template for deserializing arguments from AmArg
        template <typename A, typename B, typename C>
        inline void unstuff(A& a, B& b, C& c) const {
            BufferInputArchive ar(buf,sizeof(buf));
            ar & a & b & c;
        }
        
        
        /// Convenience template for deserializing arguments from AmArg
        template <typename A, typename B, typename C, typename D>
        inline void unstuff(A& a, B& b, C& c, D& d) const {
            BufferInputArchive ar(buf,sizeof(buf));
            ar & a & b & c & d;
        }


        /// Convenience template for deserializing arguments from AmArg
        template <typename A, typename B, typename C, typename D, typename E>
        inline void unstuff(A& a, B& b, C& c, D& d, E& e) const {
            BufferInputArchive ar(buf,sizeof(buf));
            ar & a & b & c & d & e;
        }
    };

    
    class LongAmArg;


    /// Central location for freeing managed long buffers
    inline void free_long_am_arg(unsigned long *buf) {
        delete [] buf;
    };

    /// Convenience routine to allocate a LongAmArg but using "new (unsigned long)[]"
    LongAmArg* alloc_long_am_arg(std::size_t nbyte = 0);

    class WorldGopInterface;

    /// Implements AM interface
    class WorldAmInterface {
        friend class WorldGopInterface;
    public:
        // LONG_MSG_LEN and LONG_MSG_HEADER_LEN are assumed to be an integral no. of unsigned longs
        static const int SHORT_MSG_LEN = sizeof(AmArg);                             ///< Length of short message
        static const int SHORT_MSG_HEADER_LEN = 2*sizeof(unsigned long);            ///< Length of header in short message
        static const int SHORT_MSG_USER_LEN = SHORT_MSG_LEN-SHORT_MSG_HEADER_LEN;   ///< Length of user data in short message
        static const int LONG_MSG_HEADER_LEN = 4*sizeof(unsigned long);             ///< No. of bytes reserved for long message header
        static const int LONG_MSG_LEN = 9*128*1024;                                 ///< Max length of long messages
        static const int LONG_MSG_USER_LEN = LONG_MSG_LEN-LONG_MSG_HEADER_LEN;      ///< Length of user data in long messages
        
    private:
        
        // Structure to record received but not yet processed messages
        struct qmsg {
            AmArg arg;
            ProcessID src;
            unsigned long *buf;
            size_t nbyte;
            int i;
            int count;
            bool is_managed;
            bool is_short;
            
            qmsg(ProcessID src, const AmArg& arg) 
                : arg(arg)
                , src(src)
                , buf(0)
                , nbyte(0)
                , i(0)
                , count(arg.flags&COUNT_MASK)
                , is_managed(false)
                , is_short(true)
            {};
            
            qmsg(ProcessID src, unsigned long* buf, size_t nbyte, int i, bool is_managed)
                : arg()
                , src(src)
                , buf(buf)
                , nbyte(nbyte)
                , i(i)
                , count(buf[0]&COUNT_MASK)
                , is_managed(is_managed)
                , is_short(false)
            {};

            void printout(std::ostream& s) const {
                if (is_short) {
                    unsigned long flags = arg.flags;
                    bool isbcast = flags & BCAST_MASK;
                    int count = flags & COUNT_MASK;
                    s << "qmsg(short, src=" << src 
                      << ", isbcast=" << isbcast 
                      << ", function=" << (void *) arg.function 
                      << ", count=" << count 
                      << ")\n";
                }
                else {
                    unsigned long flags = buf[0];
                    int count = flags & COUNT_MASK;
                    am_long_handlerT function = (am_long_handlerT) (buf[1]);
                    bool isbcast = flags & BCAST_MASK;
                    s << "qmsg(long, src=" << src 
                      << ", isbcast=" << isbcast 
                      << ", nbyte=" << nbyte 
                      << ", buf=" << (void *) buf 
                      << ", i=" << i
                      << ", managed=" << is_managed 
                      << ", function=" << (void *) function 
                      << ", count=" << count 
                      << ")\n";
                }
            };

        };

        // Masks used to set flags field of AmArgs
        // bits 0-7 = 1 byte counter used to ensure message ordering
        // bit-8 = unused
        // bit-9 = true if broadcast, false if point-to-point
        // bits 10 and 11 currently free
        // top 20 bits hold root of broadcast tree, if any
        static const unsigned long COUNT_MASK = 0xff;
        static const unsigned long BCAST_MASK = 0x1ul<<9;
        
        static const int NSHORT_RECV = 64;         ///< No. of posted short recv buffers
        static const int NLONG_RECV = 512;          ///< No. of posted long recv buffers
        static const int NRECV =  NSHORT_RECV + NLONG_RECV;
        static const int LOG2_NSEND = 8;
        static const int NSEND = 1<<LOG2_NSEND;    ///< Max no. of outstanding short+long Isends
        
        MPI::Request recv_handle[NRECV];  ///< Handles for AM Irecv
        mutable MPI::Request send_handle[NSEND];  ///< Handles for AM Isend
        AmArg recv_buf[NSHORT_RECV];      ///< Buffers for short AM Irecv
        AmArg send_buf[NSEND];            ///< Buffers for short AM Isend
        unsigned char* managed_send_buf[NSEND];     ///< Managed buffers for long AM Isend
        unsigned long* long_recv_buf[NLONG_RECV];   ///< Buffers for long AM Irecv
        
        // !! Careful if you reorder these to not break the constructor member initialization
        World& world; ///< The world which contains this instance of WorldAmInterface
        WorldMpiInterface& mpi;  ///< The corresponding mpi instance for brevity
        const ProcessID rank;
        const int nproc;
        const Tag short_tag; ///< Reserved tag used for short active messages over MPI
        const Tag long_tag;  ///< Reserved tag used for long active messages over MPI
        unsigned long nsent; ///< Counts no. of AM sent
        unsigned long nrecv; ///< Counts no. of AM received
        long suspended;      ///< If non-zero, AM processing is suspended
        bool debug;                   ///< If true, print debug info
        int nfree_long_recv_buf; ///< Tracks number of free long recv buffers

        std::list<qmsg*> msgq;      ///< Queue of messages received but not yet processed

        std::size_t msgq_count;     ///< For stats only: #msg currently in q
        std::size_t msgq_count_max; ///< For stats only: lifetime max of msgq_count
        std::size_t msgq_mem;       ///< For stats only: amount of memory used by q
        std::size_t msgq_mem_max;   ///< For stats only: lifetime max of msgq_mem
        std::size_t ncall_free_recv_buf; ///< For stats only: #times free_long_recv_bufs called
        std::size_t nbytes_sent, nbytes_recv; ///< For stats only: data volume
        std::size_t nsend_req_wait; ///< Number of times ran out of send requests

        unsigned char* recv_counters; ///< Used to impose sequential consistency
        unsigned char* send_counters; ///< Used to impose sequential consistency
        
        /// Private: (re)Posts the active message Irecv
        inline void post_recv(int i) {
            // Requests are consecutive in order short then long.
            if (i < NSHORT_RECV) {
                recv_handle[i] = mpi.Irecv(recv_buf+i, sizeof(AmArg), 
                                           MPI::BYTE, MPI::ANY_SOURCE, short_tag);
            }
            else {
                int j = i-NSHORT_RECV;
                recv_handle[i] = mpi.Irecv(long_recv_buf[j], LONG_MSG_LEN, 
                                           MPI::BYTE, MPI::ANY_SOURCE, long_tag);
                nfree_long_recv_buf++;
                MADNESS_ASSERT(nfree_long_recv_buf<=NLONG_RECV);
            }
        };    
        
        /// Private:
        inline void poll_short_msg_action(ProcessID src, const AmArg& arg) {
            PROFILE_MEMBER_FUNC(WorldAM);
            bool isbcast = arg.flags & BCAST_MASK;
            if (debug) print("World:",rank,"got short AM from",src,
                             "isbcast",isbcast,
                             "function",(void *) arg.function, "count", arg.flags&COUNT_MASK);
            
            if (isbcast) 
                _broadcast(arg.function, arg, arg.flags>>12);
            
            arg.function(world, src, arg);
        };
        
        /// Private:
        inline void poll_long_msg_action(ProcessID src, size_t nbyte, unsigned long* buf) {
            PROFILE_MEMBER_FUNC(WorldAM);
            unsigned long flags = buf[0];
            am_long_handlerT function = (am_long_handlerT) (buf[1]);
            bool isbcast = flags & BCAST_MASK;
            if (debug) print("World:",rank,"got  long AM from",src,
                             "nbyte",nbyte,"isbcast",isbcast,
                             "function",(void *) function, "count", flags&COUNT_MASK);
            if (isbcast) 
                _broadcast_long(function, buf, nbyte, flags>>12);
            
            function(world, src, buf, nbyte);
        };

        inline void free_managed_send_buf(int i) {
            PROFILE_MEMBER_FUNC(WorldAM);
            if (managed_send_buf[i]) {
                free_long_am_arg((unsigned long*) managed_send_buf[i]);
                managed_send_buf[i] = 0;
            }
        };


        /// Private: Finds/waits for a free send request
        inline int get_free_send_request() {
            PROFILE_MEMBER_FUNC(WorldAM);
            static int cur_msg = 0;
            static const int MASK = (1<<LOG2_NSEND)-1; // this is why we need a power of 2 messages

            // Must call poll here to keep pulling messages off the
            // network to avoid dead/livelock but don't don't need to
            // poll in other worlds or run tasks/AM.
            if (!send_handle[cur_msg].Test()) {
                nsend_req_wait++;
                World::poll_all(true);
                while (!send_handle[cur_msg].Test())  World::poll_all(true);
            }

            free_managed_send_buf(cur_msg);

            int i = cur_msg;
            cur_msg = (cur_msg + 1) & MASK;
            return i;
        };

//         /// Private: Finds/waits for a free send request
//         inline int get_free_send_request() {
//             PROFILE_MEMBER_FUNC(WorldAM);
//             // Must call poll here to keep pulling messages off the
//             // network to avoid dead/livelock but don't don't need to
//             // poll in other worlds or run tasks/AM.
//             int i;

//             // Was suspending in here but there is no need AS LONG as we
//             // suspend recursive processing of AM which is now a formal 
//             // part of the model.
//             while (!MPI::Request::Testany(NSEND, send_handle, i)) World::poll_all(true);
//             if (i == MPI::UNDEFINED) i = 0;
//             free_managed_send_buf(i);
//             return i;
//         };


        /// Private: push an active message onto the queue

        /// Separate routine to make capturing statistics easier
        inline void msgq_push_back(qmsg* qm) {
            PROFILE_MEMBER_FUNC(WorldAM);
            msgq.push_back(qm);
            msgq_count++;
            msgq_count_max = std::max(msgq_count,msgq_count_max);
            if (qm->is_managed) {
                msgq_mem += qm->nbyte;
                msgq_mem_max = std::max(msgq_mem,msgq_mem_max);
            }
        }


        /// Private: Sends a short active message setting internal flags
        inline void _send(ProcessID dest, am_handlerT op, const AmArg& arg, ProcessID root) {
            PROFILE_MEMBER_FUNC(WorldAM);
            unsigned long flags = send_counters[dest]++;
            if (root >= 0) flags |= (root << 12) | BCAST_MASK;
            if (debug) print("World:",rank,"sending short AM to",dest,
                             "op",(void *) op,"count",flags&COUNT_MASK);
            if (dest == rank) {  
                // Calls to self are appended to the queue
                AmArg* p = const_cast<AmArg*>(&arg);
                p->flags = flags;
                p->function = op;
                msgq_push_back(new qmsg(rank, arg));
            }
            else {
                int i = get_free_send_request(); 
                send_buf[i] = arg; // Must copy since using a non-blocking send
                send_buf[i].flags = flags;
                send_buf[i].function = op;
                send_handle[i] = mpi.Isend(send_buf+i, sizeof(AmArg), MPI::BYTE, dest, short_tag);
            }
            nsent++;
            nbytes_sent += sizeof(AmArg);
            //suspend();
            World::poll_all();
            //resume();
        };
        

        /// Private: Sends a long active message setting internal flags
        inline int _send_long(ProcessID dest, am_long_handlerT op, void *buf, size_t nbyte, 
                              ProcessID root, bool managed) {
            PROFILE_MEMBER_FUNC(WorldAM);
            unsigned long flags = send_counters[dest]++;
            int i;
            if (root >= 0) flags |= (root << 12) | BCAST_MASK;
            unsigned long* u = (unsigned long *)buf;
            u[0] = flags;  // Want flags at front to get at counter easily
            u[1] = (unsigned long) op;
            if (debug) 
                print("World:",rank,"sending long AM to",dest,"op",(void *) op,
                      "size",nbyte,"count",flags&COUNT_MASK);
            if (dest == rank) {  
                i = -1;
                // To enable wait_for_complete of a local long AM
                // without executing everything in the call stack
                // (which might violate ordering) we must copy the
                // darn buffer.  Various optimizations can sometimes
                // circumvent this but for now go with slow but
                // correct.  Actually it turns out this is not such
                // a big issue since WorldObj, WorldDC and WorldTask
                // all separately treat local actions so the only
                // stuff coming in here is a direct user call which
                // is not really expected.
                if (!managed) {
                    unsigned long *tmp = (unsigned long*) alloc_long_am_arg(nbyte);
                    memcpy(tmp, u, nbyte);
                    u = tmp;
                }
                msgq_push_back(new qmsg(rank, u, nbyte, i, true));
            }
            else {
                i = get_free_send_request();
                send_handle[i] = mpi.Isend(buf, nbyte, MPI::BYTE, dest, long_tag);
                if (managed) managed_send_buf[i] = (unsigned char*) buf;
            }
            nsent++;
            nbytes_sent += nbyte;
            
            //suspend();
            World::poll_all();
            //resume();
            return i;
        };

        /// Private: broadcasts a short active message to all other nodes
        inline void _broadcast(am_handlerT op, const AmArg& arg, ProcessID root) {
            ProcessID parent, child0, child1;
            mpi.binary_tree_info(root, parent, child0, child1);

            if (child0 != -1) _send(child0, op, arg, root);
            if (child1 != -1) _send(child1, op, arg, root);
        };


        /// Private: broadcasts a long active message to all other nodes
        inline void _broadcast_long(am_long_handlerT op, void *buf, size_t nbyte, ProcessID root) {
            ProcessID parent, child0, child1;
            mpi.binary_tree_info(root, parent, child0, child1);

            // Unlike short messages, long messages are not copied before
            // sending.  In this circumstance we know that the entire buffer
            // (header and user data) will be the same for both operations
            // so it is OK to reuse the buffer and am.fence after both.

            // !!!!!!!!!!!!!!! Inspection of the _send_long code suggests
            // that the above comment is not true.  Perhaps that code has
            // changed?  Anyway we seem to need to wait for completion here.
            if (child0 != -1) _send_long(child0, op, buf, nbyte, root, false);
            if (child1 != -1) _send_long(child1, op, buf, nbyte, root, false);
            fence();
        };

        long pull_msgs_into_q() {
            PROFILE_MEMBER_FUNC(WorldAM);
            MPI::Status status[NRECV];
            int ind[NRECV];
            long narrived = MPI::Request::Testsome(NRECV, recv_handle, ind, status);
            MADNESS_ASSERT(narrived != MPI::UNDEFINED);
            for (long m=0; m<narrived; m++) {
                ProcessID src = status[m].Get_source();
                int i = ind[m];
                if (i < NSHORT_RECV) {
                    msgq_push_back(new qmsg(src,recv_buf[i]));
                    post_recv(i);  // Repost short msg buffers ASAP
                }
                else {
                    size_t nbyte = status[m].Get_count(MPI::BYTE);
                    msgq_push_back(new qmsg(src, long_recv_buf[i-NSHORT_RECV], nbyte, i, false));
                    nfree_long_recv_buf--;
                    MADNESS_ASSERT(nfree_long_recv_buf>=0);
                }
            }

            return narrived;
        };


        /// Private: Copy data of long pending messages to free up recv buffers

        /// Need to do this in order to avoid deadlock when being flooded.
        void free_long_recv_bufs() {
            PROFILE_MEMBER_FUNC(WorldAM);
            ncall_free_recv_buf++; // Count times called
            // Iterate backwards thru msg list since in use buffers
            // must be at the end.
            for (std::list<qmsg*>::iterator it= msgq.end(); 
                 it != msgq.begin();) {
                --it;
                qmsg* msg = *it;
                if (!msg->is_short && !msg->is_managed) {
                    int i=msg->i;
                    //print("freeing buffer",i); std::cout.flush();
                    unsigned long* buf = (unsigned long*) alloc_long_am_arg(msg->nbyte);
                    memcpy((void*) buf, (void*) long_recv_buf[i-NSHORT_RECV], msg->nbyte);
                    msg->buf = buf;
                    msg->is_managed = true;
                    msg->i = -1;
                    post_recv(i);
                    msgq_mem += msg->nbyte;
                    msgq_mem_max = std::max(msgq_mem, msgq_mem_max);

                    if (nfree_long_recv_buf == NLONG_RECV) break;
                }
            }
        };

        void process_messages_in_q() {
            PROFILE_MEMBER_FUNC(WorldAM);
            // All incoming messages regardless of origin are put into this queue.
            //
            // Sequential consistency must be maintained using the message counters.
            //
            // Recursive/reentrant calls must be expected and therefore every time
            // we process a message we must assume that all iterators, etc., are
            // invalid and must iterate again from the beginning.

            qmsg *msg;
            do {
                if (msgq.empty()) return;
                // Attempt to find an in-order message ... as soon as
                // one is found remove it from the queue and exit the
                // search.  A null pointer for msg indicates no message.
                msg = 0;
                for (std::list<qmsg*>::iterator it = msgq.begin(); it != msgq.end(); ++it) {
                    qmsg *tmp = *it;
                    if (tmp->count == recv_counters[tmp->src]) {
                        msg = tmp;
                        msgq.erase(it);
                        break;
                    }
                }
                if (msg) {
                    suspend();
                    if (msg->is_short) {
                        poll_short_msg_action(msg->src, msg->arg);
                        nbytes_recv += sizeof(AmArg);

                        // Short message receive already reposted.
                    }
                    else {
                        poll_long_msg_action(msg->src, msg->nbyte, msg->buf);
                        if (msg->i >= 0) post_recv(msg->i);
                        if (msg->is_managed) free_long_am_arg(msg->buf);
                        nbytes_recv += msg->nbyte;
                    }
                    recv_counters[msg->src]++;
                    nrecv++;
                    msgq_count--;
                    if (msg->is_managed) msgq_mem -= msg->nbyte;
                    delete msg;
                    resume(false);
                }
            } while (msg);
        }


    public:
        // !! Careful.  This constructor is called with an only partially
        // constructed world object and the ONLY thing you can rely on
        // actually existing at this point is the MPI interface.
        WorldAmInterface(World& world) 
            : world(world)
            , mpi(world.mpi)
            , rank(mpi.rank())
            , nproc(mpi.nproc())
            , short_tag(mpi.unique_reserved_tag())
            , long_tag(mpi.unique_reserved_tag())
            , nsent(0)
            , nrecv(0)
            , suspended(0)
            , debug(false)
            , nfree_long_recv_buf(0)
            , msgq()
            , msgq_count(0)
            , msgq_count_max(0)
            , msgq_mem(0)
            , msgq_mem_max(0)
            , nbytes_sent(0)
            , nbytes_recv(0)
            , nsend_req_wait(0)
        {
            // Allocate buffers for long message receives
            for (int i=0; i<NLONG_RECV; i++) 
                long_recv_buf[i] = new unsigned long [LONG_MSG_LEN/sizeof(unsigned long)];

            // Allocate counters used to ensure messages are processed
            // in the order sent.  Sigh.  Having multiple irecv
            // buffers for efficiency means that messages may be
            // encountered out-of-order.  Thus, we have to include a 
            // counter in each message to enable ordering.
            recv_counters = new unsigned char[world.mpi.nproc()];
            send_counters = new unsigned char[world.mpi.nproc()];
            for (int i=0; i<world.mpi.nproc(); i++) 
                recv_counters[i] = send_counters[i] = 0;
            

            // Post Irecv for all short and long AM buffers
            for (int i=0; i<NRECV; i++)
                post_recv(i);
            for (int i=0; i<NSEND; i++)
                managed_send_buf[i] = 0;
        };


        virtual ~WorldAmInterface() {
            if (!MPI::Is_finalized()) fence();
            suspend();
            if (!MPI::Is_finalized()) {
                for (int i=0; i<NRECV; i++) {
                    recv_handle[i].Cancel();
                    recv_handle[i].Wait(); // Required?
                }
            }
            delete [] recv_counters;
            delete [] send_counters;
            for (int i=0; i<NLONG_RECV; i++) {
                delete [] long_recv_buf[i];
            }
        };
        
        
        /// Set debug flag to new value and return old value
        bool set_debug(bool value) {
            bool status = debug;
            debug = value;
            return status;
        };
        

        /// Suspends processing of AM
        inline void suspend() {
            suspended++;
        };
        

        /// Resumes processing of AM

        /// By default will poll for and process AM immediately after resumption
        inline void resume(bool poll=true) {
            suspended--;
            MADNESS_ASSERT(suspended >= 0);
            if (poll && !suspended) World::poll_all();
        };

        /// Mostly for debugging AM itself.  Print pending message q
        void printq() const {
            print(world.rank(),"World: number of messages in local q", msgq.size());
            for (std::list<qmsg*>::const_iterator it = msgq.begin(); it != msgq.end(); ++it) (*it)->printout(std::cout);
            std::cout.flush();
        };            

        /// Mostly for debugging AM itself.  Print bunch of stuff.
        void print_stuff() const {
            print("nsent",nsent,"nrecv",nrecv);
            std::cout << "send_counters = [";
            for (int i=0; i<world.size(); i++) {
                std::cout << (int) send_counters[i]; 
                if (i != (world.size()-1)) std::cout << ", ";
                else std::cout << "]\n";
            }
            std::cout << "recv_counters = [";
            for (int i=0; i<world.size(); i++) {
                std::cout << (int) recv_counters[i]; 
                if (i != (world.size()-1)) std::cout << ", ";
                else std::cout << "]\n";
            }
            for (int i=0; i<NSEND; i++) {
                if (send_handle[i] == MPI::Request()) {
                    print("send handle", i, "is null");
                }
                else {
                    print("send handle", i, "is active");
                }
            }
            printq();
        };

        /// dummy routine for jaguar
        static inline void usleep(int) {};


        struct FenceProbe {
            mutable MPI::Request* send_handle;
            mutable std::list<qmsg*>* plist;
            FenceProbe(MPI::Request* send_handle, std::list<qmsg*>* plist)
                : send_handle(send_handle), plist(plist) {};
            bool operator()() const {
                return MPI::Request::Testall(NSEND, send_handle) && plist->empty();
            };
        };

        /// Ensures all local AM operations have completed

        /// After this call completes, all outgoing messages will
        /// have been sent so that buffers may be reused, and all
        /// incoming messages will have executed.
        ///
        /// Tasks and other activities are executed while waiting.
        ///
        /// Calling this from an AM handler will cause a hang.
        inline void fence() {
            MADNESS_ASSERT(!suspended);
            World::await(FenceProbe(send_handle, &msgq));
            for (int i=0; i<NSEND; i++) free_managed_send_buf(i);
        };

        /// Wait for completion of a send request so you can reuse your buffer
        void wait_to_complete(int handle) {
            if (handle >= 0) {
                World::await(send_handle[handle]);
            }
        };

        void print_stats() const {
            std::cout << "\n    Active Message Statistics\n";
            std::cout << "    -------------------------\n";
            std::cout << "    Sent " << nsent << " messages (" << (nbytes_sent>>20) << " Mbytes)\n";
            std::cout << "    Recv " << nrecv << " messages (" << (nbytes_recv>>20) << " Mbytes)\n";
            std::cout << "    Current pending message q " << msgq_count << " messages (" << (msgq_mem>>20) << " Mbytes)\n";
            std::cout << "    Maximum pending message q " << msgq_count_max << " messages (" << (msgq_mem_max>>20) << " Mbytes)\n";
            std::cout << "    Number of times ran out of recv bufs " << ncall_free_recv_buf << "\n";
            std::cout << "    Number of times ran out of send reqs " << nsend_req_wait << "\n";
            std::cout.flush();
        }

        /// AM polling function (if any) for this world instance

        /// The AM processing guarantees sequential consistency so the
        /// messages are processed in the order sent.
        void poll() {
            PROFILE_MEMBER_FUNC(WorldAM);
            // We suck and process messages until nothing more arrives.
            // We don't wait for all messages to be processed.

            // Suspending now only stops processing of active messages.
            // Sucking of messages from the network happens no matter
            // what in order to avoid lots of problems.  Of course
            // we don't have an unbounded memory so don't get too
            // carried away.

            // Call here in hope of processing pending local messages
            // to free memory
            if (!suspended) process_messages_in_q();

            long narrived;
            do {
                narrived = pull_msgs_into_q();
                if (!suspended) process_messages_in_q();
                if (nfree_long_recv_buf == 0) free_long_recv_bufs();
            } while (narrived);
        };


        /// Sends a short active message
        inline void send(ProcessID dest, am_handlerT op, const AmArg& arg) {
            _send(dest, op, arg, -1);
        };
        
        
        /// Sends a long active message (non-blocking ... see details)

        /// A long active message is a contiguous buffer up to
        /// long_msg_len() bytes long.  The first
        /// long_msg_header_len() bytes of this message are
        /// reserved for internal system use, leaving the following
        /// bytes for use by the application.
        ///
        /// NB.  The send is non-blocking and the buffer cannot be
        /// reused or freed without first calling either
        /// am.wait_to_complete(handle) (which waits for just this
        /// message to be sent) or am.fence() (which waits for all
        /// messages and local actions to complete).  
        ///
        /// Returns an integer that can be used to wait for completion
        /// of the send so that the buffer can be reused.
        inline int send_long(ProcessID dest, am_long_handlerT op, void *buf, size_t nbyte) {
            return _send_long(dest, op, buf, nbyte, -1, false);
        }

        /// Sends a long active message passing ownership of buffer to World

        /// The buffer, assumed to be allocated with 
        /// \code
        ///     new unsigned char[length]
        /// \endcode
        /// will be automatically deleted when the AM is actually
        /// sent, hence unlike send_long, there is no need to call
        /// am.fence.  The user should NOT do anything else with
        /// the buffer after calling this routine.  See send_long for
        /// more info about the long AM interface.
        inline void send_long_managed(ProcessID dest, am_long_handlerT op, void *buf, size_t nbyte) {
            _send_long(dest, op, buf, nbyte, -1, true);
        }

        inline void send_long_managed(ProcessID dest, am_long_handlerT op, const std::pair<void*,size_t>& arg) {
            _send_long(dest, op, arg.first, arg.second, -1, true);
        }


        /// Broadcasts a short active message to all other nodes
        inline void broadcast(am_handlerT op, const AmArg& arg) {
            _broadcast(op, arg, rank);
        }


        /// Broadcasts a long active message to all other nodes

        /// This routine internally calls am.fence to ensure
        /// that the broadcast is complete (therefore, unlike
        /// send_long(), the buffer is available for immediate
        /// reuse).
        inline void broadcast(am_long_handlerT op, void *buf, size_t nbyte) {
            _broadcast_long(op, buf, nbyte, rank);
        }

           
        /// Send AM to p and recveive in reply an MPI message from q with given tag
        
        /// Request+reply is a common idiom so a convenience routine is provided
        inline void send_recv(ProcessID p, am_handlerT handler, const AmArg& arg, 
                                 void* buf, long count, ProcessID q, Tag tag) {
            MPI::Request req = mpi.Irecv(buf, count, MPI::BYTE, q, tag);
            send(p,handler,arg);
            World::await(req); 
        };
    };


    /// Convenience class for using arguments for long active messages
    struct LongAmArg {
        unsigned char header[WorldAmInterface::LONG_MSG_HEADER_LEN];
        unsigned char buf[8]; // Payload

        /// Convenience template for deserializing arguments from LongAmArg
        template <typename A>
        inline void unstuff(size_t nbyte, A& a) const {
            BufferInputArchive ar(buf,nbyte-sizeof(header));
            ar & a;
        }
        
        
        /// Convenience template for deserializing arguments from LongAmArg
        template <typename A, typename B>
        inline void unstuff(size_t nbyte, A& a, B& b) const {
            BufferInputArchive ar(buf,nbyte-sizeof(header));
            ar & a & b;
        }
        
        
        /// Convenience template for deserializing arguments from LongAmArg
        template <typename A, typename B, typename C>
        inline void unstuff(size_t nbyte, A& a, B& b, C& c) const {
            BufferInputArchive ar(buf,nbyte-sizeof(header));
            ar & a & b & c;
        }
        
        
        /// Convenience template for deserializing arguments from LongAmArg
        template <typename A, typename B, typename C, typename D>
        inline void unstuff(size_t nbyte, A& a, B& b, C& c, D& d) const {
            BufferInputArchive ar(buf,nbyte-sizeof(header));
            ar & a & b & c & d;
        }
        
        
        /// Convenience template for deserializing arguments from LongAmArg
        template <typename A, typename B, typename C, typename D, typename E>
        inline void unstuff(size_t nbyte, A& a, B& b, C& c, D& d, E& e) const {
            BufferInputArchive ar(buf,nbyte-sizeof(header));
            ar & a & b & c & d & e;
        }
        
        
        /// Convenience template for deserializing arguments from LongAmArg
        template <typename A, typename B, typename C, typename D, typename E, typename F>
        inline void unstuff(size_t nbyte, A& a, B& b, C& c, D& d, E& e, F& f) const {
            BufferInputArchive ar(buf,nbyte-sizeof(header));
            ar & a & b & c & d & e & f;
        }
        
        
        /// Convenience template for deserializing arguments from LongAmArg
        template <typename A, typename B, typename C, typename D, typename E, typename F, typename G>
        inline void unstuff(size_t nbyte, A& a, B& b, C& c, D& d, E& e, F& f, G& g) const {
            BufferInputArchive ar(buf,nbyte-sizeof(header));
            ar & a & b & c & d & e & f & g;
        }
        
        
        /// Convenience template for deserializing arguments from LongAmArg
        template <typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H>
        inline void unstuff(size_t nbyte, A& a, B& b, C& c, D& d, E& e, F& f, G& g, H& h) const {
            BufferInputArchive ar(buf,nbyte-sizeof(header));
            ar & a & b & c & d & e & f & g & h;
        }


        /// Convenience template for deserializing arguments from LongAmArg
        template <typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H, typename I>
        inline void unstuff(size_t nbyte, A& a, B& b, C& c, D& d, E& e, F& f, G& g, H& h, I& i) const {
            BufferInputArchive ar(buf,nbyte-sizeof(header));
            ar & a & b & c & d & e & f & g & h & i;
        }


        /// Convenience template for deserializing arguments from LongAmArg
        template <typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H, typename I, typename J>
        inline void unstuff(size_t nbyte, A& a, B& b, C& c, D& d, E& e, F& f, G& g, H& h, I& i, J& j) const {
            BufferInputArchive ar(buf,nbyte-sizeof(header));
            ar & a & b & c & d & e & f & g & h & i & j;
        }

    };


    inline LongAmArg* alloc_long_am_arg(std::size_t nbyte) {
        if (nbyte == 0) nbyte = WorldAmInterface::LONG_MSG_LEN;
        LongAmArg *buf = (LongAmArg*) new unsigned long[(nbyte-1)/sizeof(long)+1];
        return buf;
    };

    /// Convenience template for serializing arguments into LongAmArg
    template <typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H, typename I, typename J>
    inline std::pair<void*, std::size_t> new_long_am_arg(const A& a, const B& b, const C& c, const D& d, 
                                                              const E& e, const F& f, const G& g, const H& h, 
                                                              const I& i, const J& j) {
        BufferOutputArchive count;
        count & a & b & c & d & e & f & g & h & i & j;
        std::size_t size = count.size();
        LongAmArg* arg = alloc_long_am_arg(size+WorldAmInterface::LONG_MSG_HEADER_LEN);
        BufferOutputArchive ar(arg->buf,size);
        ar & a & b & c & d & e & f & g & h & i & j;
        size += WorldAmInterface::LONG_MSG_HEADER_LEN;
        return std::pair<void*, std::size_t>((void *) arg, size);
    }

    /// Convenience template for serializing arguments into LongAmArg
    template <typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H, typename I>
    inline std::pair<void*, std::size_t> new_long_am_arg(const A& a, const B& b, const C& c, const D& d, 
                                                              const E& e, const F& f, const G& g, const H& h, 
                                                              const I& i) {
        BufferOutputArchive count;
        count & a & b & c & d & e & f & g & h & i;
        std::size_t size = count.size();
        LongAmArg* arg = alloc_long_am_arg(size+WorldAmInterface::LONG_MSG_HEADER_LEN);
        BufferOutputArchive ar(arg->buf,size);
        ar & a & b & c & d & e & f & g & h & i;
        size += WorldAmInterface::LONG_MSG_HEADER_LEN;
        return std::pair<void*, std::size_t>((void *) arg, size);
    }

    /// Convenience template for serializing arguments into LongAmArg
    template <typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H>
    inline std::pair<void*, std::size_t> new_long_am_arg(const A& a, const B& b, const C& c, const D& d, 
                                                              const E& e, const F& f, const G& g, const H& h) {
        BufferOutputArchive count;
        count & a & b & c & d & e & f & g & h;
        std::size_t size = count.size();
        LongAmArg* arg = alloc_long_am_arg(size+WorldAmInterface::LONG_MSG_HEADER_LEN);
        BufferOutputArchive ar(arg->buf,size);
        ar & a & b & c & d & e & f & g & h;
        size += WorldAmInterface::LONG_MSG_HEADER_LEN;
        return std::pair<void*, std::size_t>((void *) arg, size);
    }

    /// Convenience template for serializing arguments into LongAmArg
    template <typename A, typename B, typename C, typename D, typename E, typename F, typename G>
    inline std::pair<void*, std::size_t> new_long_am_arg(const A& a, const B& b, const C& c, const D& d, 
                                                              const E& e, const F& f, const G& g) {
        BufferOutputArchive count;
        count & a & b & c & d & e & f & g;
        std::size_t size = count.size();
        LongAmArg* arg = alloc_long_am_arg(size+WorldAmInterface::LONG_MSG_HEADER_LEN);
        BufferOutputArchive ar(arg->buf,size);
        ar & a & b & c & d & e & f & g;
        size += WorldAmInterface::LONG_MSG_HEADER_LEN;
        return std::pair<void*, std::size_t>((void *) arg, size);
    }

    /// Convenience template for serializing arguments into LongAmArg
    template <typename A, typename B, typename C, typename D, typename E, typename F>
    inline std::pair<void*, std::size_t> new_long_am_arg(const A& a, const B& b, const C& c, const D& d, 
                                                              const E& e, const F& f) {
        BufferOutputArchive count;
        count & a & b & c & d & e & f;
        std::size_t size = count.size();
        LongAmArg* arg = alloc_long_am_arg(size+WorldAmInterface::LONG_MSG_HEADER_LEN);
        BufferOutputArchive ar(arg->buf,size);
        ar & a & b & c & d & e & f;
        size += WorldAmInterface::LONG_MSG_HEADER_LEN;
        return std::pair<void*, std::size_t>((void *) arg, size);
    }

    /// Convenience template for serializing arguments into LongAmArg
    template <typename A, typename B, typename C, typename D, typename E>
    inline std::pair<void*, std::size_t> new_long_am_arg(const A& a, const B& b, const C& c, const D& d, 
                                                              const E& e) {
        BufferOutputArchive count;
        count & a & b & c & d & e;
        std::size_t size = count.size();
        LongAmArg* arg = alloc_long_am_arg(size+WorldAmInterface::LONG_MSG_HEADER_LEN);
        BufferOutputArchive ar(arg->buf,size);
        ar & a & b & c & d & e;
        size += WorldAmInterface::LONG_MSG_HEADER_LEN;
        return std::pair<void*, std::size_t>((void *) arg, size);
    }

    /// Convenience template for serializing arguments into LongAmArg
    template <typename A, typename B, typename C, typename D>
    inline std::pair<void*, std::size_t> new_long_am_arg(const A& a, const B& b, const C& c, const D& d) {
        BufferOutputArchive count;
        count & a & b & c & d;
        std::size_t size = count.size();
        LongAmArg* arg = alloc_long_am_arg(size+WorldAmInterface::LONG_MSG_HEADER_LEN);
        BufferOutputArchive ar(arg->buf,size);
        ar & a & b & c & d;
        size += WorldAmInterface::LONG_MSG_HEADER_LEN;
        return std::pair<void*, std::size_t>((void *) arg, size);
    }

    /// Convenience template for serializing arguments into LongAmArg
    template <typename A, typename B, typename C>
    inline std::pair<void*, std::size_t> new_long_am_arg(const A& a, const B& b, const C& c) {
        BufferOutputArchive count;
        count & a & b & c;
        std::size_t size = count.size();
        LongAmArg* arg = alloc_long_am_arg(size+WorldAmInterface::LONG_MSG_HEADER_LEN);
        BufferOutputArchive ar(arg->buf,size);
        ar & a & b & c;
        size += WorldAmInterface::LONG_MSG_HEADER_LEN;
        return std::pair<void*, std::size_t>((void *) arg, size);
    }

    /// Convenience template for serializing arguments into LongAmArg
    template <typename A, typename B>
    inline std::pair<void*, std::size_t> new_long_am_arg(const A& a, const B& b) {
        BufferOutputArchive count;
        count & a & b;
        std::size_t size = count.size();
        LongAmArg* arg = alloc_long_am_arg(size+WorldAmInterface::LONG_MSG_HEADER_LEN);
        BufferOutputArchive ar(arg->buf,size);
        ar & a & b;
        size += WorldAmInterface::LONG_MSG_HEADER_LEN;
        return std::pair<void*, std::size_t>((void *) arg, size);
    }

    /// Convenience template for serializing arguments into LongAmArg
    template <typename A>
    inline std::pair<void*, std::size_t> new_long_am_arg(const A& a) {
        BufferOutputArchive count;
        count & a;
        std::size_t size = count.size();
        LongAmArg* arg = alloc_long_am_arg(size+WorldAmInterface::LONG_MSG_HEADER_LEN);
        BufferOutputArchive ar(arg->buf,size);
        ar & a;
        size += WorldAmInterface::LONG_MSG_HEADER_LEN;
        return std::pair<void*, std::size_t>((void *) arg, size);
    }

}



#endif
