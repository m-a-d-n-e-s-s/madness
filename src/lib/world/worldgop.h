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

  $LastChangedDate$
  $Rev$
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
        inline T operator()(const T& a, const T& b) const {return a+b;};
    };
    
    template <typename T>
    struct WorldMultOp {
        inline T operator()(const T& a, const T& b) const {return a*b;};
    };
    
    template <typename T>
    struct WorldMaxOp {
        inline T operator()(const T& a, const T& b) const {return a>b? a : b;};
    };
    
    template <typename T>
    struct WorldAbsMaxOp {
        inline T operator()(const T& a, const T& b) const {return abs(a)>abs(b)? abs(a) : abs(b);};
    };
    
    template <typename T>
    struct WorldMinOp {
        inline T operator()(const T& a, const T& b) const {return a<b? a : b;};
    };
    
    template <typename T>
    struct WorldAbsMinOp {
        inline T operator()(const T& a, const T& b) const {return abs(a)<abs(b)? abs(a) : abs(b);};
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
        int nproc;
        const Tag bcast_tag;    //< Reserved tag used for broadcasting
        const Tag gsum_tag;     //< Reserved tag used for up-tree part of global sum
        const Tag gfence_tag;   //< Reserved tag used for up-tree part of global fence
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
            , bcast_tag(mpi.unique_reserved_tag())
            , gsum_tag(mpi.unique_reserved_tag())
            , gfence_tag(mpi.unique_reserved_tag())
            , debug(false)
            {};

        
        /// Set debug flag to new value and return old value
        bool set_debug(bool value) {
            bool status = debug;
            debug = value;
            return status;
        };
        
        /// Synchronizes all processes in communicator ... does NOT fence pending AM or tasks
        void barrier() {
            long i = rank;
            global_sum(i);
            if (i != nproc*(nproc-1)/2) error("bad value after sum in barrier");
        };


        /// Synchronizes all processes in communicator AND globally ensures no pending AM or tasks
        
        /// Runs Dykstra-like termination algorithm on binary tree by
        /// summing (nsent-nrecv) over all proceses, applying
        /// local fence before computing local sum.  If the global sum
        /// is zero, then we are sure that all tasks and AM are processed
        /// and there no AM in flight.
        ///
        /// By default fences both tasks and AM.  If the optional argument amonly is
        /// set to true, if fences only the AM.
        void fence(bool amonly = false) {
            long sum, sum0=0, sum1=0;
            MPI::Request req0, req1;
            ProcessID parent, child0, child1;
            mpi.binary_tree_info(0, parent, child0, child1);
            do {
                if (child0 != -1) req0 = mpi.Irecv(&sum0, sizeof(sum0), MPI::BYTE, child0, gfence_tag);
                if (child1 != -1) req1 = mpi.Irecv(&sum1, sizeof(sum1), MPI::BYTE, child1, gfence_tag);
                if (child0 != -1) World::await(req0);
                if (child1 != -1) World::await(req1);
                
                if (amonly)
                    am.fence();
                else
                    taskq.fence();

                sum = sum0 + sum1 + am.nsent - am.nrecv;
                
                if (parent != -1) {
                    req0 = mpi.Isend(&sum, sizeof(sum), MPI::BYTE, parent, gfence_tag);
                    World::await(req0);
                }
                
                broadcast(sum);
                if (debug) {
                    print(rank,"fence:",sum0,sum1,sum);
                    //usleep(1000000);
                }
            } while (sum);

            // If deferred cleanup occured we need another fence, but
            // it will be much cheaper the second time since everyone
            // is already synchronized.
            //
            // Uh?  Why do we need another fence?  Commented this out
            // until I can reconvince myself it really is necessary.
            //if (world.do_deferred_cleanup()) fence();
            world.do_deferred_cleanup();
        };
        
        
        /// Broadcasts bytes from process root while still processing AM & tasks
        
        /// Optimizations can be added for long messages
        void broadcast(void* buf, size_t nbyte, ProcessID root) {
            MPI::Request req0, req1;
            ProcessID parent, child0, child1;
            mpi.binary_tree_info(root, parent, child0, child1);
            
            if (parent != -1) {
                req0 = mpi.Irecv(buf, nbyte, MPI::BYTE, parent, bcast_tag);
                World::await(req0);
            }
            
            if (child0 != -1) req0 = mpi.Isend(buf, nbyte, MPI::BYTE, child0, bcast_tag);
            if (child1 != -1) req1 = mpi.Isend(buf, nbyte, MPI::BYTE, child1, bcast_tag);
            
            if (child0 != -1) World::await(req0);
            if (child1 != -1) World::await(req1);
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


        /// Inplace global reduction (like MPI all_reduce) while still processing AM & tasks
        
        /// Optimizations can be added for long messages
        template <typename T, class opT>
        void global_reduce(T* buf, size_t nelem, opT op) {
            MPI::Request req0, req1;
            ProcessID parent, child0, child1;
            mpi.binary_tree_info(0, parent, child0, child1);
            
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
        inline void global_sum(T* buf, size_t nelem) {
            global_reduce< T, WorldSumOp<T> >(buf, nelem, WorldSumOp<T>());
        }
        
        /// Inplace global min while still processing AM & tasks
        template <typename T>
        inline void global_min(T* buf, size_t nelem) {
            global_reduce< T, WorldMinOp<T> >(buf, nelem, WorldMinOp<T>());
        }
        
        /// Inplace global max while still processing AM & tasks
        template <typename T>
        inline void global_max(T* buf, size_t nelem) {
            global_reduce< T, WorldMaxOp<T> >(buf, nelem, WorldMaxOp<T>());
        }
        
        /// Inplace global absmin while still processing AM & tasks
        template <typename T>
            inline void global_absmin(T* buf, size_t nelem) {
            global_reduce< T, WorldAbsMinOp<T> >(buf, nelem, WorldAbsMinOp<T>());
        }
        
        /// Inplace global absmax while still processing AM & tasks
        template <typename T>
            inline void global_absmax(T* buf, size_t nelem) {
            global_reduce< T, WorldAbsMaxOp<T> >(buf, nelem, WorldAbsMaxOp<T>());
        }
        
        /// Inplace global product while still processing AM & tasks
        template <typename T>
            inline void global_product(T* buf, size_t nelem) {
            global_reduce< T, WorldMultOp<T> >(buf, nelem, WorldMultOp<T>());
        }
        
        /// Global sum of a scalar while still processing AM & tasks
        template <typename T>
            void global_sum(T& a) {global_sum(&a, 1);}
        
    };        
}


#endif
