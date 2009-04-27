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


#ifndef WORLDGOP_H
#define WORLDGOP_H

/// \file worldgop.h
/// \brief Implements global operations

/// If you can recall the Intel hypercubes, their comm lib used GOP as
/// the abbreviation.



namespace madness {
    template <typename T>
    struct WorldSumOp {
        inline T operator()(const T& a, const T& b) const {
            return a+b;
        };
    };

    template <typename T>
    struct WorldMultOp {
        inline T operator()(const T& a, const T& b) const {
            return a*b;
        };
    };

    template <typename T>
    struct WorldMaxOp {
        inline T operator()(const T& a, const T& b) const {
            return a>b? a : b;
        };
    };

    template <typename T>
    struct WorldAbsMaxOp {
        inline T operator()(const T& a, const T& b) const {
            return abs(a)>abs(b)? abs(a) : abs(b);
        };
    };

    template <typename T>
    struct WorldMinOp {
        inline T operator()(const T& a, const T& b) const {
            return a<b? a : b;
        };
    };

    template <typename T>
    struct WorldAbsMinOp {
        inline T operator()(const T& a, const T& b) const {
            return abs(a)<abs(b)? abs(a) : abs(b);
        };
    };


    /// Provides collectives that interoperate with the AM and task interfaces

    /// If native AM interoperates with MPI we probably should map these to MPI.
    class WorldGopInterface {
    private:
        World& world;
        WorldMpiInterface& mpi;
        WorldAmInterface& am;
        WorldTaskQueue& taskq;
        ProcessID rank;
        const int nproc;
        bool debug;
    public:

        // In the World constructor can ONLY rely on MPI and MPI being initialized
        WorldGopInterface(World& world)
                : world(world)
                , mpi(world.mpi)
                , am(world.am)
                , taskq(world.taskq)
                , rank(world.mpi.rank())
                , nproc(world.mpi.size())
                , debug(false) {};


        /// Set debug flag to new value and return old value
        bool set_debug(bool value) {
            bool status = debug;
            debug = value;
            return status;
        };

        /// Synchronizes all processes in communicator ... does NOT fence pending AM or tasks
        void barrier() {
            long i = rank;
            sum(i);
            if (i != nproc*(nproc-1)/2) error("bad value after sum in barrier");
        };


        /// Synchronizes all processes in communicator AND globally ensures no pending AM or tasks

        /// Runs Dykstra-like termination algorithm on binary tree by
        /// locally ensuring ntask=0 and all am sent and processed,
        /// and then participating in a global sum of nsent and nrecv.
        /// Then globally checks that nsent=nrecv and that both are
        /// constant over two traversals.  We are then we are sure
        /// that all tasks and AM are processed and there no AM in
        /// flight.
        void fence() {
            PROFILE_MEMBER_FUNC(WorldGopInterface);
            unsigned long nsent_prev=0, nrecv_prev=1; // invalid initial condition
            SafeMPI::Request req0, req1;
            ProcessID parent, child0, child1;
            mpi.binary_tree_info(0, parent, child0, child1);
            Tag gfence_tag = mpi.unique_tag();
            int npass = 0;

            while (1) {
                uint64_t sum0[2]={0,0}, sum1[2]={0,0}, sum[2];
                if (child0 != -1) req0 = mpi.Irecv((void*) &sum0, sizeof(sum0), MPI::BYTE, child0, gfence_tag);
                if (child1 != -1) req1 = mpi.Irecv((void*) &sum1, sizeof(sum1), MPI::BYTE, child1, gfence_tag);
                world.taskq.fence();
                if (child0 != -1) World::await(req0);
                if (child1 != -1) World::await(req1);

                bool finished;
                uint64_t ntask1, nsent1, nrecv1, ntask2, nsent2, nrecv2;
                do {
                    taskq.fence();

                    // Since the number of outstanding tasks and number of AM sent/recv
                    // don't share a critical section read each twice and ensure they
                    // are unchanged to ensure that are consistent ... they don't have
                    // to be current.
                    ntask1 = taskq.size();
                    nsent1 = am.nsent;
                    nrecv1 = am.nrecv;

                    ntask2 = taskq.size();
                    nsent2 = am.nsent;
                    nrecv2 = am.nrecv;

                    finished = (ntask2==0) && (ntask1==0) && (nsent1==nsent2) && (nrecv1==nrecv2);
                }
                while (!finished);

                sum[0] = sum0[0] + sum1[0] + nsent2; // Must use values read above
                sum[1] = sum0[1] + sum1[1] + nrecv2;

                if (parent != -1) {
                    req0 = mpi.Isend(&sum, sizeof(sum), MPI::BYTE, parent, gfence_tag);
                    World::await(req0);
                }

                // While we are probably idle free unused communication buffers
                world.am.free_managed_buffers();

                bool dowork = (npass==0) || (ThreadPool::size()==0);
                broadcast(&sum, sizeof(sum), 0, dowork);
                npass++;

                //madness::print("GOPFENCE", npass, sum[0], nsent_prev, sum[1], nrecv_prev);

                if (sum[0]==sum[1] && sum[0]==nsent_prev && sum[1]==nrecv_prev) break;

                nsent_prev = sum[0];
                nrecv_prev = sum[1];

            };
            world.do_deferred_cleanup();
        };


        /// Broadcasts bytes from process root while still processing AM & tasks

        /// Optimizations can be added for long messages
        void broadcast(void* buf, size_t nbyte, ProcessID root, bool dowork = true) {
            SafeMPI::Request req0, req1;
            ProcessID parent, child0, child1;
            mpi.binary_tree_info(root, parent, child0, child1);
            Tag bcast_tag = mpi.unique_tag();

            //print("BCAST TAG", bcast_tag);

            if (parent != -1) {
                req0 = mpi.Irecv(buf, nbyte, MPI::BYTE, parent, bcast_tag);
                World::await(req0, dowork);
            }

            if (child0 != -1) req0 = mpi.Isend(buf, nbyte, MPI::BYTE, child0, bcast_tag);
            if (child1 != -1) req1 = mpi.Isend(buf, nbyte, MPI::BYTE, child1, bcast_tag);

            if (child0 != -1) World::await(req0, dowork);
            if (child1 != -1) World::await(req1, dowork);
        };


        /// Broadcasts typed contiguous data from process root while still processing AM & tasks

        /// Optimizations can be added for long messages
        template <typename T>
        inline void broadcast(T* buf, size_t nelem, ProcessID root) {
            broadcast((void *) buf, nelem*sizeof(T), root);
        }

        /// Broadcast of a scalar from node 0 to all other nodes
        template <typename T>
        void broadcast(T& t) {
            broadcast(&t, 1, 0);
        }

        /// Broadcast of a scalar from node root to all other nodes
        template <typename T>
        void broadcast(T& t, ProcessID root) {
            broadcast(&t, 1, root);
        }

        /// Broadcast a serializable object

        /// Current dumb version assumes object fits in 1MB
        /// ... you are free to add intelligence.
        template <typename objT>
        void broadcast_serializable(objT& obj, ProcessID root) {
            size_t BUFLEN;
            if (world.rank() == root) {
                BufferOutputArchive count;
                count & obj;
                BUFLEN = count.size();
            }
            broadcast(BUFLEN, root);

            unsigned char* buf = new unsigned char[BUFLEN];
            if (world.rank() == root) {
                BufferOutputArchive ar(buf,BUFLEN);
                ar & obj;
            }
            broadcast(buf, BUFLEN, root);
            if (world.rank() != root) {
                BufferInputArchive ar(buf,BUFLEN);
                ar & obj;
            }
            delete [] buf;
        }

        /// Inplace global reduction (like MPI all_reduce) while still processing AM & tasks

        /// Optimizations can be added for long messages and to reduce the memory footprint
        template <typename T, class opT>
        void reduce(T* buf, size_t nelem, opT op) {
            SafeMPI::Request req0, req1;
            ProcessID parent, child0, child1;
            mpi.binary_tree_info(0, parent, child0, child1);
            Tag gsum_tag = mpi.unique_tag();

            T* buf0 = new T[nelem];
            T* buf1 = new T[nelem];

            if (child0 != -1) req0 = mpi.Irecv(buf0, nelem*sizeof(T), MPI::BYTE, child0, gsum_tag);
            if (child1 != -1) req1 = mpi.Irecv(buf1, nelem*sizeof(T), MPI::BYTE, child1, gsum_tag);

            if (child0 != -1) {
                World::await(req0);
                for (long i=0; i<(long)nelem; i++) buf[i] = op(buf[i],buf0[i]);
            }
            if (child1 != -1) {
                World::await(req1);
                for (long i=0; i<(long)nelem; i++) buf[i] = op(buf[i],buf1[i]);
            }

            delete [] buf0;
            delete [] buf1;

            if (parent != -1) {
                req0 = mpi.Isend(buf, nelem*sizeof(T), MPI::BYTE, parent, gsum_tag);
                World::await(req0);
            }

            broadcast(buf, nelem, 0);
        }

        /// Inplace global sum while still processing AM & tasks
        template <typename T>
        inline void sum(T* buf, size_t nelem) {
            reduce< T, WorldSumOp<T> >(buf, nelem, WorldSumOp<T>());
        }

        /// Inplace global min while still processing AM & tasks
        template <typename T>
        inline void min(T* buf, size_t nelem) {
            reduce< T, WorldMinOp<T> >(buf, nelem, WorldMinOp<T>());
        }

        /// Inplace global max while still processing AM & tasks
        template <typename T>
        inline void max(T* buf, size_t nelem) {
            reduce< T, WorldMaxOp<T> >(buf, nelem, WorldMaxOp<T>());
        }

        /// Inplace global absmin while still processing AM & tasks
        template <typename T>
        inline void absmin(T* buf, size_t nelem) {
            reduce< T, WorldAbsMinOp<T> >(buf, nelem, WorldAbsMinOp<T>());
        }

        /// Inplace global absmax while still processing AM & tasks
        template <typename T>
        inline void absmax(T* buf, size_t nelem) {
            reduce< T, WorldAbsMaxOp<T> >(buf, nelem, WorldAbsMaxOp<T>());
        }

        /// Inplace global product while still processing AM & tasks
        template <typename T>
        inline void product(T* buf, size_t nelem) {
            reduce< T, WorldMultOp<T> >(buf, nelem, WorldMultOp<T>());
        }

        /// Global sum of a scalar while still processing AM & tasks
        template <typename T>
        void sum(T& a) {
            sum(&a, 1);
        }

        /// Global max of a scalar while still processing AM & tasks
        template <typename T>
        void max(T& a) {
            max(&a, 1);
        }

        /// Global min of a scalar while still processing AM & tasks
        template <typename T>
        void min(T& a) {
            min(&a, 1);
        }

        /// Concatenate an STL vector of serializable stuff onto node 0
        template <typename T>
        std::vector<T> concat0(const std::vector<T>& v, size_t bufsz=1024*1024) {
            SafeMPI::Request req0, req1;
            ProcessID parent, child0, child1;
            mpi.binary_tree_info(0, parent, child0, child1);
            Tag gsum_tag = mpi.unique_tag();

            unsigned char* buf0 = new unsigned char[bufsz];
            unsigned char* buf1 = new unsigned char[bufsz];

            if (child0 != -1) req0 = mpi.Irecv(buf0, bufsz, MPI::BYTE, child0, gsum_tag);
            if (child1 != -1) req1 = mpi.Irecv(buf1, bufsz, MPI::BYTE, child1, gsum_tag);

            std::vector<T> left, right;
            if (child0 != -1) {
                World::await(req0);
                BufferInputArchive ar(buf0, bufsz);
                ar & left;
            }
            if (child1 != -1) {
                World::await(req1);
                BufferInputArchive ar(buf1, bufsz);
                ar & right;
                for (unsigned int i=0; i<right.size(); i++) left.push_back(right[i]);
            }

            for (unsigned int i=0; i<v.size(); i++) left.push_back(v[i]);

            if (parent != -1) {
                BufferOutputArchive ar(buf0, bufsz);
                ar & left;
                req0 = mpi.Isend(buf0, ar.size(), MPI::BYTE, parent, gsum_tag);
                World::await(req0);
            }

            delete [] buf0;
            delete [] buf1;

            if (parent == -1) return left;
            else return std::vector<T>();
        }
    };
}


#endif
