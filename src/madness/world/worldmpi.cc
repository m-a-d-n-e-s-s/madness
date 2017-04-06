/*
  This file is part of MADNESS.

  Copyright (C) 2013  Virginia Tech

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

/**
 \file worldmpi.cc
 \brief Several variables needed for \c WorldMPI.
 \ingroup mpi
*/

#include <madness/world/worldmpi.h>

namespace madness {
    namespace detail {

        /// \addtogroup mpi
        /// @{

        std::shared_ptr<WorldMpi> WorldMpi::world_mpi;

        /// Flag storing if MADNESS is responsible for MPI.

        /// \todo Verify the above brief description.
        bool WorldMpi::own_mpi = false;

#ifdef MADNESS_USE_BSEND_ACKS
        /// MPI buffer.
        char* WorldMpi::mpi_ack_buffer[MADNESS_ACK_BUFF_SIZE];
#endif // MADNESS_USE_BSEND_ACKS

        /// @}

    } // namespace detail
} // namespace madness
