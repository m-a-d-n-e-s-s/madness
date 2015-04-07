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

#include <madness/world/safempi.h>
#include <madness/world/madness_exception.h>

namespace SafeMPI {

    madness::SCALABLE_MUTEX_TYPE charon;

    void Intracomm::binary_tree_info(int root, int& parent, int& child0, int& child1) {
        const int np = Get_size();
        const int me = (Get_rank() + np - root) % np; // Renumber processes so root has me=0
        parent = (me == 0 ? -1 : (((me - 1) >> 1) + root) % np); // Parent in binary tree
        child0 = (me << 1) + 1 + root; // Left child
        child1 = child0 + 1; // Right child
        const int np_plus_root = np + root;
        if(child0 < np_plus_root)
            child0 %= np;
        else
            child0 = -1;
        if(child1 < np_plus_root)
            child1 %= np;
        else
            child1 = -1;
    }

    // The logic here is to intended to cause an error if someone tries to use
    // this constructor from somewhere else.

    struct Intracomm::WorldInitObject {
        MPI_Comm comm_world() const { return MPI_COMM_WORLD; }
    };

    Intracomm::Intracomm(const WorldInitObject& init) :
        pimpl(new Impl(init.comm_world(), -1, -1, false))
    { }

    Intracomm COMM_WORLD = Intracomm::WorldInitObject();

} //namespace SafeMPI
