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

#include <limits>
#include <madness/world/worldgop.h>
#include <madness/world/MADworld.h>
#ifdef MADNESS_HAS_GOOGLE_PERF_MINIMAL
#include <gperftools/malloc_extension.h>
#endif
namespace madness {


    /// Synchronizes all processes in communicator AND globally ensures no pending AM or tasks

    /// Runs Dykstra-like termination algorithm on binary tree by
    /// locally ensuring ntask=0 and all am sent and processed,
    /// and then participating in a global sum of nsent and nrecv.
    /// Then globally checks that nsent=nrecv and that both are
    /// constant over two traversals.  We are then we are sure
    /// that all tasks and AM are processed and there no AM in
    /// flight.
    void WorldGopInterface::fence(bool debug) {
        PROFILE_MEMBER_FUNC(WorldGopInterface);
        unsigned long nsent_prev=0, nrecv_prev=1; // invalid initial condition
        SafeMPI::Request req0, req1;
        ProcessID parent, child0, child1;
        world_.mpi.binary_tree_info(0, parent, child0, child1);
        Tag gfence_tag = world_.mpi.unique_tag();
        Tag bcast_tag = world_.mpi.unique_tag();
        int npass = 0;

        //double start = wall_time();

      if (debug)
        madness::print(world_.rank(), ": WORLD.GOP.FENCE: entering fence loop, gfence_tag=", gfence_tag, " bcast_tag=", bcast_tag);

      while (1) {
            uint64_t sum0[2]={0,0}, sum1[2]={0,0}, sum[2];
            if (child0 != -1) req0 = world_.mpi.Irecv((void*) &sum0, sizeof(sum0), MPI_BYTE, child0, gfence_tag);
            if (child1 != -1) req1 = world_.mpi.Irecv((void*) &sum1, sizeof(sum1), MPI_BYTE, child1, gfence_tag);
            world_.taskq.fence();
            if (child0 != -1) World::await(req0);
            if (child1 != -1) World::await(req1);

            if (debug && (child0 != -1 || child1 != -1))
              madness::print(world_.rank(), ": WORLD.GOP.FENCE: npass=", npass, " received messages from children={", child0, ",", child1, "} gfence_tag=", gfence_tag);

            bool finished;
            uint64_t ntask1, nsent1, nrecv1, ntask2, nsent2, nrecv2;
            do {
                world_.taskq.fence();

                // Since the number of outstanding tasks and number of AM sent/recv
                // don't share a critical section read each twice and ensure they
                // are unchanged to ensure that are consistent ... they don't have
                // to be current.

                ntask1 = world_.taskq.size();
                nsent1 = world_.am.nsent;
                nrecv1 = world_.am.nrecv;

                __asm__ __volatile__ (" " : : : "memory");

                ntask2 = world_.taskq.size();
                nsent2 = world_.am.nsent;
                nrecv2 = world_.am.nrecv;

                __asm__ __volatile__ (" " : : : "memory");

                finished = (ntask2==0) && (ntask1==0) && (nsent1==nsent2) && (nrecv1==nrecv2);
            }
            while (!finished);

            sum[0] = sum0[0] + sum1[0] + nsent2; // Must use values read above
            sum[1] = sum0[1] + sum1[1] + nrecv2;

            if (parent != -1) {
                req0 = world_.mpi.Isend(&sum, sizeof(sum), MPI_BYTE, parent, gfence_tag);
                if (debug)
                  madness::print(world_.rank(), ": WORLD.GOP.FENCE: npass=", npass, " sent message to parent=", parent, " gfence_tag=", gfence_tag);
                World::await(req0);
                if (debug)
                  madness::print(world_.rank(), ": WORLD.GOP.FENCE: npass=", npass, " parent=", parent, ", confirmed receipt");
            }

            // While we are probably idle free unused communication buffers
            //world_.am.free_managed_buffers();

            //bool dowork = (npass==0) || (ThreadPool::size()==0);
            bool dowork = true;
            broadcast(&sum, sizeof(sum), 0, dowork, bcast_tag);
            ++npass;

            if (debug)
              madness::print(world_.rank(), ": WORLD.GOP.FENCE: npass=", npass, " sum0=", sum[0], " nsent_prev=", nsent_prev, " sum1=", sum[1], " nrecv_prev=", nrecv_prev);

            if (sum[0]==sum[1] && sum[0]==nsent_prev && sum[1]==nrecv_prev) {
              if (debug)
                madness::print(world_.rank(), ": WORLD.GOP.FENCE: npass=", npass, " exiting fence loop");
              break;
            }

//                 if (wall_time() - start > 1200.0) {
//                     std::cout << rank() << " FENCE " << nsent2 << " "
//                         << nsent_prev << " " << nrecv2 << " " << nrecv_prev
//                         << " " << sum[0] << " " << sum[1] << " " << npass
//                         << " " << taskq.size() << std::endl;
//                     std::cout.flush();
//                     //myusleep(1000);
//                     MADNESS_ASSERT(0);
//                 }

            nsent_prev = sum[0];
            nrecv_prev = sum[1];

        };
        world_.am.free_managed_buffers(); // free up communication buffers
        deferred_->do_cleanup();
#ifdef MADNESS_HAS_GOOGLE_PERF_MINIMAL
        MallocExtension::instance()->ReleaseFreeMemory();
//        print("clearing memory");
#endif
      if (debug)
        madness::print(world_.rank(), ": WORLD.GOP.FENCE: done with fence in ", npass, (npass > 1 ? " loops" : " loop"));
    }


    /// Broadcasts bytes from process root while still processing AM & tasks

    /// Optimizations can be added for long messages
    void WorldGopInterface::broadcast(void* buf, size_t nbyte, ProcessID root, bool dowork, Tag bcast_tag) {
        MADNESS_ASSERT(nbyte <= size_t(std::numeric_limits<int>::max()));
        SafeMPI::Request req0, req1;
        ProcessID parent, child0, child1;
        world_.mpi.binary_tree_info(root, parent, child0, child1);
        if(bcast_tag < 0)
            bcast_tag = world_.mpi.unique_tag();

        //print("BCAST TAG", bcast_tag);

        if (parent != -1) {
            req0 = world_.mpi.Irecv(buf, nbyte, MPI_BYTE, parent, bcast_tag);
            World::await(req0, dowork);
        }

        if (child0 != -1) req0 = world_.mpi.Isend(buf, nbyte, MPI_BYTE, child0, bcast_tag);
        if (child1 != -1) req1 = world_.mpi.Isend(buf, nbyte, MPI_BYTE, child1, bcast_tag);

        if (child0 != -1) World::await(req0, dowork);
        if (child1 != -1) World::await(req1, dowork);
    }

} // namespace madness
