#ifndef MADNESS_SYSTOLIC_H
#define MADNESS_SYSTOLIC_H

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

#include <madness/world/MADworld.h>
#include <utility>
#include <madness/tensor/tensor.h>
#include <madness/tensor/distributed_matrix.h>

#ifdef HAVE_INTEL_TBB
# include <tbb/parallel_for.h>
#endif

namespace madness {

    /// Base class for parallel algorithms that employ a systolic loop to generate all row pairs in parallel
    template <typename T>
    class SystolicMatrixAlgorithm : public TaskInterface {
    private:
        DistributedMatrix<T>& A;
        const int64_t nproc;            ///< No. of processes with rows of the matrix (not size of world)
        const int64_t coldim;           ///< A(coldim,rowdim)
        const int64_t rowdim;           ///< A(coldim,rowdim)
        const int64_t nlocal;           ///< No. of local pairs
        const ProcessID rank;           ///< Rank of current process
        const int tag;                  ///< MPI tag to be used for messages
        std::vector<T*> iptr, jptr;     ///< Indirection for implementing cyclic buffer !! SHOULD BE VOLATILE ?????
        std::vector<int64_t> map;       ///< Used to keep track of actual row indices

#ifdef HAVE_INTEL_TBB
        void iteration(const int nthread) {

            // Parallel initialization hook
            tbb::parallel_for(0, nthread, [=](const int id) {
                this->start_iteration_hook(TaskThreadEnv(nthread, id));
            });

            if (nlocal > 0) {
//                int64_t ilo, ihi;
//                A.local_colrange(ilo, ihi);

                const int neven = coldim + (coldim&0x1);
                const int pairlo = rank*A.coltile()/2;

                for (int loop=0; loop<(neven-1); ++loop) {

                    // This loop is parallelized over threads
                    tbb::parallel_for(0, nthread,
                        [this,neven,pairlo,nthread,loop](const int id) {
                            for (int pair=id; pair<nlocal; pair+=nthread) {

                                int rp = neven/2-1-(pair+pairlo);
                                int iii = (rp+loop)%(neven-1);
                                int jjj = (2*neven-2-rp+loop)%(neven-1);
                                if (rp == 0) jjj = neven-1;

                                iii = map[iii];
                                jjj = map[jjj];

                                if (jptr[pair]) {
                                    this->kernel(iii, jjj, iptr[pair], jptr[pair]);
                                }
                            }
                        });

                    cycle();

                }
            }

            // Parallel finalization hook
            tbb::parallel_for(0, nthread, [=](const int id) {
                this->end_iteration_hook(TaskThreadEnv(nthread, id));
            });

        }
#else
        void iteration(const TaskThreadEnv& env) {

            env.barrier();
            start_iteration_hook(env);
            env.barrier();

            if (nlocal > 0) {
                int64_t ilo, ihi;
                A.local_colrange(ilo, ihi);

                int neven = coldim + (coldim&0x1);

                int pairlo = rank*A.coltile()/2;

                int threadid = env.id();
                int nthread = env.nthread();

                for (int loop=0; loop<(neven-1); ++loop) {

                    // This loop is parallelized over threads
                    for (int pair=env.id(); pair<nlocal; pair+=nthread) {

                        int rp = neven/2-1-(pair+pairlo);
                        int iii = (rp+loop)%(neven-1);
                        int jjj = (2*neven-2-rp+loop)%(neven-1);
                        if (rp == 0) jjj = neven-1;

                        iii = map[iii];
                        jjj = map[jjj];

                        if (jptr[pair]) {
                            kernel(iii, jjj, iptr[pair], jptr[pair]);
                        }
                    }
                    env.barrier();

                    if (threadid == 0) cycle();

                    env.barrier();
                }
            }

            end_iteration_hook(env);

            env.barrier();
        }
#endif // HAVE_INTEL_TBB

        /// Call this after iterating to restore correct order of rows in original matrix

        /// At the end of each iteration the matrix rows are logically back in
        /// their correct order.  However, due to indirection to reduce data motion,
        /// if the local column dimension is not a factor of the number of cycles
        /// the underlying data may be in a different order.  This restores sanity.
        ///
        /// Only one thread should invoke this routine
        void unshuffle() {
            if (nlocal <= 0) return;
            Tensor<T>& t = A.data();
            Tensor<T> tmp(2L, t.dims(), false);
            T* tp = tmp.ptr();
            for (int64_t i=0; i<nlocal; ++i) {
                memcpy(tp+i*rowdim, iptr[i], rowdim*sizeof(T));
                if (jptr[i]) {
                    memcpy(tp+(i+nlocal)*rowdim, jptr[i], rowdim*sizeof(T));
                }
                iptr[i] = &t(i,0);
                jptr[i] = &t(i+nlocal,0);
            }
            memcpy(t.ptr(), tmp.ptr(), t.size()*sizeof(T));

            if (rank==(nproc-1) && (coldim&0x1)) jptr[nlocal-1] = 0;
        }

        /// Cycles data around the loop ... only one thread should invoke this
        void cycle() {
            if (coldim <= 2) return; // No cycling necessary
            if (nlocal <= 0) {       // Nothing local
                MADNESS_ASSERT(rank >= nproc);
                return;
            }

            // Check assumption that tiling put incomplete tile at the end
            MADNESS_ASSERT(A.local_coldim() == A.coltile()  ||  rank == (nproc-1));

            const ProcessID left = rank-1; //Invalid values are not used
            const ProcessID right = rank+1;

            /*
              Consider matrix (10,*) distributed with coltile=4 over
              three processors.

              .   0 1 2 3      4 5 6 7      8 9

              This is divided up as follows into this initial
              configuration for the loop

              .            P=0          P=1         P=2
              .                  msg          msg
              .   i    -->0-->1  -->   4-->5  -->    8  -->
              .       ^                                   |  msg
              .       |                         <---------
              .   j    <--2<--3  <--   6<--7  <--|   9
              .                  msg          msg

              The first and last processes in the loop have to wrap ... others
              just pass left and right.  Note that 9 stays put.

              Note that the algorithm is assuming distribution puts equal
              amount of data on all nodes except the last.

              The i data is considered as flowing to the right.
              The j data is considered as flowing to the left.


              Hence, we should explore the pairs in this order
              (n-1 sets of n/2 pairs)

              .          P=0         P=1        P=2
              .          0  1        4  5       8
              .          2  3        6  7       9

              .          2  0        1  4       5
              .          3  6        7  8       9

              .          3  2        0  1       4
              .          6  7        8  5       9

              .          6  3        2  0       1
              .          7  8        5  4       9

              .          7  6        3  2       0
              .          8  5        4  1       9

              .          8  7        6  3       2
              .          5  4        1  0       9

              .          5  8        7  6       3
              .          4  1        0  2       9

              .          4  5        8  7       6
              .          1  0        2  3       9

              .          1  4        5  8       7
              .          0  2        3  6       9
            */

            // Copy end elements before they are overwritten
            T* ilast  = iptr[nlocal-1];
            T* jfirst = jptr[0];

            // Cycle local pointers
            for (int64_t i=0; i<nlocal-1; ++i) {
                iptr[nlocal-i-1] = iptr[nlocal-i-2];
                jptr[i] = jptr[i+1];
            }

            World& world = A.get_world();

            if (nproc == 1) {
                iptr[0] = jfirst;
                jptr[nlocal-2] = ilast;
            }
            else if (rank == 0) {
                iptr[0] = jfirst;
                world.mpi.Send(ilast, rowdim, right, tag);
                jptr[nlocal-1] = ilast;
                world.mpi.Recv(ilast, rowdim, right, tag);
            }
            else if (rank == (nproc-1)) {
                if (nlocal > 1) {
                    iptr[0] = jfirst;
                    jptr[nlocal-2] = ilast;
                }
                std::vector<T> buf(rowdim);
                SafeMPI::Request req = world.mpi.Irecv(&buf[0], rowdim, left, tag);
                world.mpi.Send(iptr[0], rowdim, left, tag);
                world.await(req,false);
                std::memcpy(iptr[0], &buf[0], rowdim*sizeof(T));
            }
            else {
                std::vector<T> buf1(rowdim);
                std::vector<T> buf2(rowdim);
                SafeMPI::Request req1 = world.mpi.Irecv(&buf1[0], rowdim, left, tag);
                SafeMPI::Request req2 = world.mpi.Irecv(&buf2[0], rowdim, right, tag);
                world.mpi.Send( ilast, rowdim, right, tag);
                world.mpi.Send(jfirst, rowdim,  left, tag);
                world.await(req1,false);
                world.await(req2,false);
                std::memcpy(ilast, &buf2[0], rowdim*sizeof(T)); //world.mpi.Recv( ilast, rowdim, right, tag);
                std::memcpy(jfirst, &buf1[0], rowdim*sizeof(T)); //world.mpi.Recv(jfirst, rowdim,  left, tag);

                iptr[0] = jfirst;
                jptr[nlocal-1] = ilast;
            }
        }

        /// Get the task id

        /// \param id The id to set for this task
        virtual void get_id(std::pair<void*,unsigned short>& id) const {
            PoolTaskInterface::make_id(id, *this);
        }

    public:
        /// A must be a column distributed matrix with an even column tile >= 2

        /// It is assumed that it is the main thread invoking this.
        /// @param[in,out] A The matrix on which the algorithm is performed and modified in-place
        /// @param[in] tag The MPI tag used for communication (obtain from \c world.mpi.comm().unique_tag() )
        /// @param[in] nthread The number of local threads to use (default is main thread all threads in the pool)
        SystolicMatrixAlgorithm(DistributedMatrix<T>& A, int tag, int nthread=ThreadPool::size()+1)
            : A(A)
            , nproc(A.process_coldim()*A.process_rowdim())
            , coldim(A.coldim())
            , rowdim(A.rowdim())
            , nlocal((A.local_coldim()+1)/2)
            , rank(A.get_world().rank())
            , tag(tag)
            , iptr(nlocal)
            , jptr(nlocal)
            , map(coldim+(coldim&0x1))
        {
            TaskInterface::set_nthread(nthread);

            MADNESS_ASSERT(A.is_column_distributed() && (nproc==1 || (A.coltile()&0x1)==0));

            // Initialize vectors of pointers to matrix rows)
            Tensor<T>& t = A.data();

            //madness::print(nproc, coldim, rowdim, nlocal, rank, tag);

            for (int64_t i=0; i<nlocal; ++i) {
                iptr[i] = &t(i,0);
                jptr[i] = &t(i+nlocal,0);
            }

            // If no. of rows is odd, last process should have an empty last row
            if (rank==(nproc-1) && (coldim&0x1)) jptr[nlocal-1] = 0;

            // Initialize map from logical index order to actual index order

            int neven = (coldim+1)/2;
            int ii=0;
            for (ProcessID p=0; p<nproc; ++p) {
                int64_t lo, hi;
                A.get_colrange(p, lo, hi);
                int p_nlocal = (hi - lo + 2)/2;
                //print("I think process",p,"has",lo,hi,p_nlocal);
                for (int i=0; i<p_nlocal; ++i) {
                    map[ii+i] = lo+i;
                    //map[coldim-ii-nlocal+i] = lo+i+nlocal;
                    map[ii+i+neven] = lo+i+p_nlocal;
                }
                ii += p_nlocal;
            }

            std::reverse(map.begin(),map.begin()+neven);

            //print("MAP", map);
        }

        virtual ~SystolicMatrixAlgorithm() {}

        /// Threadsafe routine to apply the operation to rows i and j of the matrix

        /// @param[in] i First row index in the matrix
        /// @param[in] j Second row index in the matrix
        /// @param[in] rowi Pointer to row \c i of the matrix (to be modified by kernel in-place)
        /// @param[in] rowj Pointer to row \c j of the matrix (to be modified by kernel in-place)
        virtual void kernel(int i, int j, T* rowi, T* rowj) = 0;


        /// Invoked simultaneously by all threads after each sweep to test for convergence

        /// There is a thread barrier before and after the invocation of this routine
        /// @param[in] env The madness thread environment in case synchronization between threads is needed during computation of the convergence condition.
        virtual bool converged(const TaskThreadEnv& env) const = 0;


        /// Invoked by all threads at the start of each iteration

        /// There is a thread barrier before and after the invocation of this routine
        /// @param[in] env The madness thread environment in case synchronization between threads is needed during startup.
        virtual void start_iteration_hook(const TaskThreadEnv& env) {}


        /// Invoked by all threads at the end of each iteration before convergence test

        /// There is a thread barrier before and after the invocation of this routine.
        /// Note that the \c converged() method is \c const whereas this can modify the class.
        /// @param[in] env The madness thread environment in case synchronization between threads is needed during startup.
        virtual void end_iteration_hook(const TaskThreadEnv& env) {}


#ifdef HAVE_INTEL_TBB
        /// Invoked by the task queue to run the algorithm with multiple threads

        /// This is a collective call ... all processes in world should submit this task
        void run(World& world, const TaskThreadEnv& env) {
            const int nthread = env.nthread();
            bool done = false;
            do {
                iteration(nthread);
                done = tbb::parallel_reduce(tbb::blocked_range<int>(0,nthread), true,
                    [=] (const tbb::blocked_range<int>& range, bool init) -> bool {
                        for(int id = range.begin(); id < range.end(); ++id)
                            init = init &&
                                this->converged(TaskThreadEnv(nthread, id));
                        return init;
                    },
                    [] (const bool l, const bool r) { return l && r; });

            } while (!done);

            unshuffle();
        }
#else
        /// Invoked by the task queue to run the algorithm with multiple threads

        /// This is a collective call ... all processes in world should submit this task
        void run(World& world, const TaskThreadEnv& env) {
            do {
                iteration(env);
            } while (!converged(env));

            if (env.id() == 0) unshuffle();

            env.barrier();
        }
#endif // HAVE_INTEL_TBB

        /// Invoked by the user to run the algorithm with one thread mostly for debugging

        /// This is a collective call ... all processes in world should call this routine.
        void solve_sequential() {
            run(A.get_world(), TaskThreadEnv(1,0,0));
        }

        /// Returns length of row
        int64_t get_rowdim() const {return rowdim;}


        /// Returns length of column
        int64_t get_coldim() const {return coldim;}

        /// Returns a reference to the world
        World& get_world() const {
            return A.get_world();
        }

        /// Returns rank of this process in the world
        ProcessID get_rank() const {
            return rank;
        }
    };
}

#endif
