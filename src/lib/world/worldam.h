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
        ulong flags;   //< System managed state 
        mutable am_handlerT function; //< System managed handler
    public:
        ulong buf[14]; //< System dependent length
        
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


    class WorldGopInterface;

    /// Implements AM interface
    class WorldAmInterface {
        friend class WorldGopInterface;
    public:
        // LONG_MSG_LEN and LONG_HEADER_LEN are assumed to be an integral no. of unsigned longs
        static const int SHORT_MSG_LEN = sizeof(AmArg);                             //< Length of short message
        static const int SHORT_MSG_HEADER_LEN = 2*sizeof(unsigned long);            //< Length of header in short message
        static const int SHORT_MSG_USER_LEN = SHORT_MSG_LEN-SHORT_MSG_HEADER_LEN;   //< Length of user data in short message
        static const int LONG_MSG_HEADER_LEN = 4*sizeof(unsigned long);             //< No. of bytes reserved for long message header
        static const int LONG_MSG_LEN = 128*1024;                                   //< Max length of long messages
        static const int LONG_MSG_USER_LEN = LONG_MSG_LEN-LONG_MSG_HEADER_LEN;      //< Length of user data in long messages
        
    private:
        // Masks used to set flags field of AmArgs
        // bits 0-7 = 1 byte counter used to ensure message ordering
        // bit-8 = unused
        // bit-9 = true if broadcast, false if point-to-point
        // bits 10 and 11 currently free
        // top 20 bits hold root of broadcast tree, if any
        static const unsigned long COUNT_MASK = 0xff;
        static const unsigned long BCAST_MASK = 0x1ul<<9;
        
        static const int NSHORT_RECV = 8;         //< No. of posted short recv buffers
        static const int NLONG_RECV = 8;          //< No. of posted long recv buffers
        static const int NRECV =  NSHORT_RECV + NLONG_RECV;
        static const int NSEND = 10;              //< Max no. of outstanding short+long Isends
        
        MPI::Request recv_handle[NRECV];  //< Handles for AM Irecv
        mutable MPI::Request send_handle[NSEND];  //< Handles for AM Isend
        AmArg recv_buf[NSHORT_RECV];      //< Buffers for short AM Irecv
        AmArg send_buf[NSEND];            //< Buffers for short AM Isend
        unsigned char* long_send_buf[NSEND];        //< Managed buffers for long AM Isend
        unsigned long* long_recv_buf[NLONG_RECV];   //< Buffers for long AM Irecv
	unsigned long* flag_ptrs[NRECV];  //< Points to flags in incoming message headers
        
        // !! Careful if you reorder these to not break the constructor member initialization
        World& world; //< The world which contains this instance of WorldAmInterface
        WorldMpiInterface& mpi;  //< The corresponding mpi instance for brevity
        const ProcessID rank;
        const int nproc;
        const Tag short_tag; //< Reserved tag used for short active messages over MPI
        const Tag long_tag;  //< Reserved tag used for long active messages over MPI
        volatile unsigned long nsent; //< Counts no. of AM sent
        volatile unsigned long nrecv; //< Counts no. of AM received
        volatile long suspended;      //< If non-zero, AM processing is suspended
        bool debug;                   //< If true, print debug info
        unsigned char* recv_counters; //< Used to impose sequential consistency
        unsigned char* send_counters; //< Used to impose sequential consistency
        
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
            }
        };    
        
        /// Private:
        
        /// Note that this might be called directly for local messages so
        /// don't increment AM recv'd count in here.
        inline void poll_short_msg_action(ProcessID src, const AmArg& arg) {
            bool isbcast = arg.flags & BCAST_MASK;
            if (debug) print("World:",rank,"poll got short msg from",src,
                             "isbcast",isbcast,
                             "function",(void *) arg.function);
            
            if (isbcast) 
                _broadcast(arg.function, arg, arg.flags>>12);
            
            arg.function(world, src, arg);
        };
        
        /// Private:
        
        /// Note that this might be called directly for local messages so
        /// don't increment AM recv'd count in here.
        inline void poll_long_msg_action(ProcessID src, size_t nbyte, unsigned long* buf) {
            unsigned long flags = buf[0];
            am_long_handlerT function = (am_long_handlerT) (buf[1]);
            bool isbcast = flags & BCAST_MASK;
            if (debug) print("World:",rank,"poll got long msg from",src,
                             "nbyte",nbyte,"isbcast",isbcast,
                             "function",(void *) function);
            if (isbcast) 
                _broadcast_long(function, buf, nbyte, flags>>12);
            
            function(world, src, buf, nbyte);
        };

        inline void free_long_send_buf(int i) {
            if (long_send_buf[i]) {
                delete [] long_send_buf[i];
                long_send_buf[i] = 0;
            }
        };
        
        // Used to handle out-of-order active messages
        struct qmsg {
            ProcessID src;
            int msgcount;
            int i;
            size_t nbyte;
            AmArg tmp;
        };

        /// Checks that AM msg is in order and then processes it

        /// Returns true if the message was handled
        inline bool process_am(const qmsg& msg) {
            ProcessID src = msg.src;
            int target_count = recv_counters[src];
            if (msg.msgcount != target_count) return false;
            recv_counters[src]++;

            int i = msg.i;
            if (i < NSHORT_RECV) {
                poll_short_msg_action(src, msg.tmp);
            }
            else {
	        poll_long_msg_action(src,msg.nbyte, long_recv_buf[i-NSHORT_RECV]);
                post_recv(i);
            }
            nrecv++;                   // Must be after handler invocation
            return true;
        };
        

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
        {
            // Allocate buffers for long message receives
            for (int i=0; i<NLONG_RECV; i++) 
                long_recv_buf[i] = new unsigned long [LONG_MSG_LEN/sizeof(unsigned long)];

	    // Make pointers to flags in short and long incoming messages
            for (int i=0; i<NSHORT_RECV; i++) 
                flag_ptrs[i] = &recv_buf[i].flags;
            for (int i=0; i<NLONG_RECV; i++) 
	        flag_ptrs[i+NSHORT_RECV] = long_recv_buf[i];

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
                long_send_buf[i] = 0;
            
        };
        
        virtual ~WorldAmInterface() {
//             if (!MPI::Is_finalized()) {
//                 fence();
//                 suspend();
//                 for (int i=0; i<NRECV; i++) {
//                     recv_handle[i].Cancel();
//                     recv_handle[i].Wait(); // Required?
//                 }
//                 delete [] recv_counters;
//                 delete [] send_counters;
//                 for (int i=0; i<NLONG_RECV; i++) {
//                     delete [] long_recv_buf[i];
//                 }
//             }
            suspend();
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
        inline void resume() {
            suspended--;
            if (suspended < 0) {
                print("AM: suspended was negative --- mismatched suspend/resume pairing");
                suspended = 0;
            };
        };
        

        /// dummy routine for jaguar
        static inline void usleep(int) {};

	int get_msg_counter(int i) const {
  	    return *(flag_ptrs[i])&COUNT_MASK;
	};
	

        /// Ensures all AM send operations (short and long) have completed
        inline void fence() {
            // Must call poll here to avoid livelock with two
            // processes trying to talk to each other when both
            // are out of recv buffers.  However, don't need to
            // poll in other worlds or run tasks.
            while (!MPI::Request::Testall(NSEND, send_handle)) 
                poll();
            for (int i=0; i<NSEND; i++) free_long_send_buf(i);
        };
        
        /// AM polling function (if any) for this world instance

        /// The AM processing guarantees sequential consistency so the
        /// messages are processed in the order sent.  Not sure if
        /// this guarantee can or should be continued once we go to
        /// multi-core dispatch.
        void poll() {
            //if (debug) print("World:",rank,"in poll with suspended =",suspended);
            if (suspended) return;
            suspend();

            // Well, we have had to switch back to using testany/all since
            // at least in MPICH 2-1.0.3 Get_status does not ensure
            // progress.  

            // Testall also seems to be broken so now we just call
            // test manually.

            // If the message is out of order, it must already be in a
            // buffer due to the non-overtaking property of MPI
            // messages of the same type from the same sender.  Grab
            // all ready messages via testall, extract the necessary
            // info, and then process them in order of msgcount to
            // ensure ordering.

            // Well either I misunderstood the non-overtaking guarantee
            // or I should be asking for my money back (MPICH is free).
            // Hence, the logic below does no longer assumes this.

            // Catch is that a handler may invoke poll() which can
            // then grab other messages which are then out of order
            // w.r.t. the ones here.  Thus, if we want handlers to be
            // able to poll in this world we would have to make this
            // internal queue static.  Doable (bounded buffer
            // problem?), but for now we simply disable polling in
            // this world while in the AM handlers which is perfectly
            // reasonable.

            const int MAXINQ = 10*NRECV;
            qmsg q[NRECV];

            MPI::Status status;
            int narrived;// Counts no. messages arriving each pass thru while loop
            int nq = 0;  // Counts no. of pending out-of-order messages
            do {
                narrived = 0;
                // Get messages from MPI and process those that are in-order
                for (int i=0; i<NRECV; i++) {
                    // request could be null for out-of-order long messages
                    if (recv_handle[i] != MPI::REQUEST_NULL && 
                        recv_handle[i].Test(status)) {
                        narrived++;
                        if (nq >= MAXINQ) error("AM: exceeded max OOO q");
                        q[nq].src = status.Get_source();
                        MADNESS_ASSERT(q[nq].src>=0 && q[nq].src<nproc);
                        q[nq].msgcount = get_msg_counter(i);
                        q[nq].i = i;
                        q[nq].nbyte = status.Get_count(MPI::BYTE);
                        if (i < NSHORT_RECV) {  // Repost short msg buffers ASAP even if OOO
                            q[nq].tmp = recv_buf[i];
                            post_recv(i);
                        }
                        if (!process_am(q[nq])) {
                            if (debug) print("AM OOO recvd",nq,"src",q[nq].src,"msgcount",
                                             q[nq].msgcount,int(recv_counters[q[nq].src]),
                                             "i",q[nq].i,"nb",q[nq].nbyte);
                            nq++;
                        }
                    }
                }
                
                // Now take a look at OOO messages.
                if (nq) {
                    if (debug) {
                        print("AM: OOOq",nq);
                        for (int msg=0; msg<nq; msg++) {
                            print("    ",msg,"src",q[msg].src,"msgcount",
                                  q[msg].msgcount,int(recv_counters[q[msg].src]),
                                  "i",q[msg].i,"nb",q[msg].nbyte);
                        }
                    }

                    for (int msg=nq-1; msg>=0; msg--) {
                        if(process_am(q[msg])) {
                            nq--;
                            if (msg != nq) q[msg] = q[nq];
                        }
                    }

                    if (nq) {
                        if (debug) print("AM: OOO second pass",nq);
                        for (int msg=nq-1; msg>=0; msg--) {
                            if(process_am(q[msg])) {
                                nq--;
                                if (msg != nq) q[msg] = q[nq];
                            }
                        }
                        if (nq) {
                            if (debug) print("AM: OOO third pass - yielding",nq);
                        }
                    }
                }
                

            } while (narrived + nq);
            resume();
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
        /// reused or freed without first calling fence().
        inline void send_long(ProcessID dest, am_long_handlerT op, void *buf, size_t nbyte) {
            _send_long(dest, op, buf, nbyte, -1, false);
        }


        /// Sends a long active message passing ownership of buffer to World

        /// The buffer, assumed to be allocated with 
        /// \code
        ///     new unsigned char[length]
        /// \endcode
        /// will be automatically deleted when the AM is actually
        /// sent, hence unlike send_long, there is no need to call
        /// am_fence.  The user should NOT doing anything else with
        /// the buffer after calling this routine.  See send_long for
        /// more info about the long AM interface.
        inline void send_long_managed(ProcessID dest, am_long_handlerT op, void *buf, size_t nbyte) {
            _send_long(dest, op, buf, nbyte, -1, true);
        }


        /// Broadcasts a short active message to all other nodes
        inline void broadcast(am_handlerT op, const AmArg& arg) {
            _broadcast(op, arg, rank);
        }


        /// Broadcasts a long active message to all other nodes

        /// This routine internally calls fence to ensure
        /// that the broadcast is complete (therefore, unlike
        /// send_long(), the buffer is available for immediate
        /// reuse).
        inline void broadcast(am_long_handlerT op, void *buf, size_t nbyte) {
            _broadcast_long(op, buf, nbyte, rank);
        }

        
    private:
        /// Private: Finds/waits for a free send request
        inline int get_free_send_request() {
            // Testany() might be really slow.  It may be better to
            // keep track of free handles ourselves and only call
            // testany when full.  On the other hand, calling test
            // could be necessary to make progress.
            int i;
            while (!MPI::Request::Testany(NSEND, send_handle, i)) {
                // Must call poll here to avoid livelock with two
                // processes trying to talk to each other when both
                // are out of recv buffers.  However, don't need to
                // poll in other worlds or run tasks.
                poll();
            }
            if (i == MPI::UNDEFINED) i = 0;
            else free_long_send_buf(i);
            return i;
        };


        /// Private: Sends a short active message setting internal flags
        inline void _send(ProcessID dest, am_handlerT op, const AmArg& arg, ProcessID root) {
            unsigned long flags = send_counters[dest]++;
            if (root >= 0) flags |= (root << 12) | BCAST_MASK;
            int i = get_free_send_request();
            send_buf[i] = arg; // Must copy since using a non-blocking send
            send_buf[i].flags = flags;
            send_buf[i].function = op;
            if (debug) print(rank,"sending short AM to",dest,"op",(void *) op,"flags",root);
            if (dest == rank) {  // Don't send local messages ... just invoke action
                
                // !!!! NEED TO BREAK LONG TREES OF CALLS THAT CAN EXPLODE THE STACK

                poll_short_msg_action(dest, send_buf[i]);
            }
            else {
                send_handle[i] = mpi.Isend(send_buf+i, sizeof(AmArg), MPI::BYTE, dest, short_tag);
                nsent++;
            }
            World::poll_all();
        };
        

        /// Private: Sends a long active message setting internal flags
        inline void _send_long(ProcessID dest, am_long_handlerT op, void *buf, size_t nbyte, 
                               ProcessID root, bool managed) {
            unsigned long flags = send_counters[dest]++;
            if (root >= 0) flags |= (root << 12) | BCAST_MASK;
            unsigned long* u = (unsigned long *)buf;
            u[0] = flags;  // Want flags at front to get at counter easily
            u[1] = (unsigned long) op;
            if (debug) print(rank,"sending long AM to",dest,"op",(void *) op,"size",nbyte,"flags",root);
            if (dest == rank) {  // Don't send local messages ... just invoke action
                
                // !!!! NEED TO BREAK LONG TREES OF CALLS THAT CAN EXPLODE THE STACK

                poll_long_msg_action(dest, nbyte, (unsigned long*) buf);
            }
            else {
                int i = get_free_send_request();
                send_handle[i] = mpi.Isend(buf, nbyte, MPI::BYTE, dest, long_tag);
                if (managed) long_send_buf[i] = (unsigned char*) buf;
                nsent++;
            }
            World::poll_all();
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
            // so it is OK to reuse the buffer and fence after both.
            if (child0 != -1) _send_long(child0, op, buf, nbyte, root, false);
            if (child1 != -1) _send_long(child1, op, buf, nbyte, root, false);
            fence();
        };


    public:
        
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
        unsigned char buf[WorldAmInterface::LONG_MSG_USER_LEN];

        /// Convenience template for serializing arguments into LongAmArg
        
        /// Returns total size of active message include header
        template <typename A>
        inline size_t stuff(const A& a) {
            BufferOutputArchive ar(buf,sizeof(buf));
            ar & a;
            return ar.size()+sizeof(header);
        }


        /// Convenience template for serializing arguments into LongAmArg
        
        /// Returns size of active message include header
        template <typename A, typename B>
        inline size_t stuff(const A& a, const B& b) {
            BufferOutputArchive ar(buf,sizeof(buf));
            ar & a & b;
            return ar.size()+sizeof(header);
        }


        /// Convenience template for serializing arguments into LongAmArg
        
        /// Returns size of active message include header
        template <typename A, typename B, typename C>
        inline size_t stuff(const A& a, const B& b, const C& c) {
            BufferOutputArchive ar(buf,sizeof(buf));
            ar & a & b & c;
            return ar.size()+sizeof(header);
        }
        
        
        /// Convenience template for serializing arguments into LongAmArg
        
        /// Returns size of active message include header
        template <typename A, typename B, typename C, typename D>
        inline size_t stuff(const A& a, const B& b, const C& c, const D& d) {
            BufferOutputArchive ar(buf,sizeof(buf));
            ar & a & b & c & d;
            return ar.size()+sizeof(header);
        }
        
        
        /// Convenience template for serializing arguments into LongAmArg
        
        /// Returns size of active message include header
        template <typename A, typename B, typename C, typename D, typename E>
        inline size_t stuff(const A& a, const B& b, const C& c, const D& d, const E& e) {
            BufferOutputArchive ar(buf,sizeof(buf));
            ar & a & b & c & d & e;
            return ar.size()+sizeof(header);
        }
        
        
        /// Convenience template for serializing arguments into LongAmArg
        
        /// Returns size of active message include header
        template <typename A, typename B, typename C, typename D, typename E, typename F>
        inline size_t stuff(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f) {
            BufferOutputArchive ar(buf,sizeof(buf));
            ar & a & b & c & d & e & f;
            return ar.size()+sizeof(header);
        }
        
        
        /// Convenience template for serializing arguments into LongAmArg
        
        /// Returns size of active message include header
        template <typename A, typename B, typename C, typename D, typename E, typename F, typename G>
        inline size_t stuff(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g) {
            BufferOutputArchive ar(buf,sizeof(buf));
            ar & a & b & c & d & e & f & g;
            return ar.size()+sizeof(header);
        }
        
        
        /// Convenience template for serializing arguments into LongAmArg
        
        /// Returns size of active message include header
        template <typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H>
        inline size_t stuff(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g, const H& h) {
            BufferOutputArchive ar(buf,sizeof(buf));
            ar & a & b & c & d & e & f & g & h;
            return ar.size()+sizeof(header);
        }
        
        
        /// Convenience template for serializing arguments into LongAmArg
        
        /// Returns size of active message include header
        template <typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H, typename I>
        inline size_t stuff(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g, const H& h, const I& i) {
            BufferOutputArchive ar(buf,sizeof(buf));
            ar & a & b & c & d & e & f & g & h & i;
            return ar.size()+sizeof(header);
        }
        
        
        /// Convenience template for serializing arguments into LongAmArg
        
        /// Returns size of active message include header
        template <typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H, typename I, typename J>
        inline size_t stuff(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g, const H& h, const I& i, const J& j) {
            BufferOutputArchive ar(buf,sizeof(buf));
            ar & a & b & c & d & e & f & g & h & i & j;
            return ar.size()+sizeof(header);
        }
        
        
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

}



#endif
