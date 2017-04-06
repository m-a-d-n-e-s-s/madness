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

#include <madness/world/worldam.h>
#include <madness/world/MADworld.h>
#include <madness/world/worldmpi.h>
#include <sstream>

namespace madness {



    WorldAmInterface::WorldAmInterface(World& world)
            : nsend(DEFAULT_NSEND)
            , send_req(nullptr)
            , worldid(0) // worldid is initialized in the World constructor
            , rank(world.mpi.Get_rank())
            , nproc(world.mpi.Get_size())
            , cur_msg(0)
            , nsent(0)
            , nrecv(0)
            , map_to_comm_world(nproc)
    {
        lock();

        // Initialize the number of send buffers
        const char* mad_send_buffs = getenv("MAD_SEND_BUFFERS");
        if(mad_send_buffs) {
            std::stringstream ss(mad_send_buffs);
            ss >> nsend;
            // Check that the number of send buffers is reasonable.
            if(nsend < 32) {
                nsend = DEFAULT_NSEND;
                std::cerr << "!!! WARNING: MAD_SEND_BUFFERS must be at least 32.\n"
                          << "!!! WARNING: Increasing MAD_SEND_BUFFERS to " << nsend << ".\n";
            }
        }

        // Allocate send buffers and requests
        send_req.reset(new SendReq[nsend]);

        for (int i=0; i<nsend; ++i) send_req[i].set((AmArg*) 0,RMI::Request());

        std::vector<int> fred(nproc);
        for (int i=0; i<nproc; ++i) fred[i] = i;
        world.mpi.comm().Get_group().Translate_ranks(nproc,
                                                     &fred[0], SafeMPI::COMM_WORLD.Get_group(),
                                                     &map_to_comm_world[0]);

        // for (int i=0; i<nproc; ++i) {
        //     std::cout << "map " << i << " " << map_to_comm_world[i] << std::endl;
        // }

        unlock();
    }

    WorldAmInterface::~WorldAmInterface() {
        if(!SafeMPI::Is_finalized()) {
            while (free_managed_buffers() != nsend) myusleep(100);
        }
        // otherwise the send buffers are freed when the WorldAMInterface::send_req is freed
    }

} // namespace madness
