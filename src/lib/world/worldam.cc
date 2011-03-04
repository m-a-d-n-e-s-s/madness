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


  $Id: $
*/

#include <world/worldam.h>
#include <world/world.h>
#include <world/worldmpi.h>

namespace madness {

    void AmArg::set_world(World* world) const {
        worldid = world->id();
    }

    archive::BufferInputArchive AmArg::make_input_arch() const {
        return archive::BufferInputArchive(buf(),size());
    }

    archive::BufferOutputArchive AmArg::make_output_arch() const {
        return archive::BufferOutputArchive(buf(),size());
    }

    World* AmArg::get_world() const {
        return World::world_from_id(worldid);
    }


    void WorldAmInterface::free_managed_send_buf(int i) {
        // WE ASSUME WE ARE INSIDE A CRITICAL SECTION WHEN IN HERE
        if (managed_send_buf[i]) {
            free_am_arg(managed_send_buf[i]);
            managed_send_buf[i] = 0;
        }
    }

    int WorldAmInterface::get_free_send_request() {
        // WE ASSUME WE ARE INSIDE A CRITICAL SECTION WHEN IN HERE
//             // Sequentially loop looking for next free request.
//             while (!send_req[cur_msg].Test()) {
//                 cur_msg++;
//                 if (cur_msg >= NSEND) cur_msg = 0;
//                 myusleep(5);
//             }

        // Wait for oldest request to complete
        while (!send_req[cur_msg].Test()) {
            // If the oldest message has still not completed then there is likely
            // severe network or end-point congestion, so pause for 100us in a rather
            // abitrary attempt to decreate the injection rate.  The server thread
            // is still polling every 1us (which is required to suck data off the net
            // and by some engines to ensure progress on sends).
            myusleep(100);
        }

        free_managed_send_buf(cur_msg);
        int result = cur_msg;
        cur_msg++;
        if (cur_msg >= NSEND) cur_msg = 0;

        return result;
    }

    void WorldAmInterface::handler(void *buf, std::size_t nbyte) {
        // It will be singled threaded since only the RMI receiver
        // thread will invoke it ... however note that nrecv will
        // be read by the main thread during fence operations.
        AmArg* arg = (AmArg*)(buf);
        am_handlerT func = arg->get_func();
        World* world = arg->get_world();
        MADNESS_ASSERT(arg->size() + sizeof(AmArg) == nbyte);
        MADNESS_ASSERT(world);
        MADNESS_ASSERT(func);
        func(*arg);
        world->am.nrecv++;  // Must be AFTER execution of the function
    }

    RMI::Request WorldAmInterface::isend(ProcessID dest, am_handlerT op, const AmArg* arg, int attr, bool managed) {
        arg->set_world(&world);
        arg->set_src(rank);
        arg->set_func(op);
        arg->clear_flags(); // Is this the right place for this?

        MADNESS_ASSERT(arg->get_world());
        MADNESS_ASSERT(arg->get_func());

        // Map dest from world's communicator to comm_world
        dest = map_to_comm_world[dest];

        lock();    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        nsent++;
        int i = get_free_send_request();
        send_req[i] = RMI::isend(arg, arg->size()+sizeof(AmArg), dest, handler, attr);
        if (managed) managed_send_buf[i] = (AmArg*)(arg);
        unlock();  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        return send_req[i];
    }

    WorldAmInterface::WorldAmInterface(World& world)
            : world(world)
            , rank(world.mpi.Get_rank())
            , nproc(world.mpi.Get_size())
            , cur_msg(0)
            , nsent(0)
            , nrecv(0)
            , map_to_comm_world(nproc)
    {
        lock();
        for (int i=0; i<NSEND; ++i) managed_send_buf[i] = 0;

        std::vector<int> fred(nproc);
        for (int i=0; i<nproc; ++i) fred[i] = i;
        MPI::Group::Translate_ranks(world.mpi.comm().Get_group(), nproc, &fred[0],
                                    MPI::COMM_WORLD.Get_group(), &map_to_comm_world[0]);

        // for (int i=0; i<nproc; ++i) {
        //     std::cout << "map " << i << " " << map_to_comm_world[i] << std::endl;
        // }

        unlock();
    }

    WorldAmInterface::~WorldAmInterface() {
        for (int i=0; i<NSEND; ++i) {
            while (!send_req[i].Test()) {
                myusleep(100);
            }
            free_managed_send_buf(i);
        }
    }
//    RMI::Request WorldAmInterface::isend(ProcessID dest, am_handlerT op, const AmArg* arg, int attr) {
//        std::cerr << "ISEND_ING AM\n";
//        return isend(dest, op, arg, attr, false);
//    }

    void WorldAmInterface::send(ProcessID dest, am_handlerT op, const AmArg* arg, int attr) {
        isend(dest, op, arg, attr, true);
    }

    void WorldAmInterface::free_managed_buffers() {
        int ind[NSEND];
        lock(); // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        int n = SafeMPI::Request::Testsome(NSEND, send_req, ind);
        if (n != MPI_UNDEFINED) {
            for (int i=0; i<n; ++i) {
                free_managed_send_buf(ind[i]);
            }
        }
        unlock(); // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    }

} // namespace madness
