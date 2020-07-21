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

/**
 \file worldinit.h
 \brief Declares the functions that initialize the parallel runtime.
 \ingroup world
*/

#ifndef MADNESS_WORLD_WORLDINIT_H__INCLUDED
#define MADNESS_WORLD_WORLDINIT_H__INCLUDED

// Madness world header files needed by world
#include <madness/world/worldmpi.h>

/// \addtogroup world
/// @{

namespace madness {

    class World;
    class WorldTaskQueue;
    class WorldAmInterface;
    class WorldGopInterface;

    /// redirects standard output and error to rank-specific files

    /// @param[in] world the World object that determines rank of this process
    /// @param[in] split if true, write standard output to log.<rank> and standard error to err.<rank>,
    ///            otherwise write both standard output and error to log.<rank>. The default is false.
    void redirectio(const World& world, bool split = false);

    /// Initialize the MADNESS runtime.

    /// Call this once at the very top of your main program to initialize the
    /// MADNESS runtime. This function should be called instead of \c MPI_Init()
    /// or \c MPI_Init_thread().
    /// \param[in,out] argc Application argument count.
    /// \param[in,out] argv Application argument values.
    /// \param[in] quiet If false, will announce to \c std::cout on rank 0 when
    ///            the runtime has been initialized.
    /// \return A reference to the default \c World, which is constructed with
    ///     \c MPI_COMM_WORLD.
    World& initialize(int& argc, char**& argv, bool quiet = false);

    /// Initialize the MADNESS runtime.

    /// Call this once at the very top of your main program to initialize the
    /// MADNESS runtime. This function should be called instead of \c MPI_Init()
    /// or \c MPI_Init_thread().
    /// \param[in,out] argc Application argument count.
    /// \param[in,out] argv Application argument values.
    /// \param comm The communicator that should be used to construct the
    ///     \c World object.
    /// \param[in] quiet If false, will announce to \c std::cout on rank 0 when
    ///            the runtime has been initialized.
    /// \return A reference to the \c World constructed with \c comm.
    World& initialize(int& argc, char**& argv, const SafeMPI::Intracomm& comm,
        bool quiet = false);

    /// Initialize the MADNESS runtime.

    /// Call this once at the very top of your main program to initialize the
    /// MADNESS runtime. This function should be called instead of \c MPI_Init()
    /// or \c MPI_Init_thread().
    /// \param[in,out] argc Application argument count.
    /// \param[in,out] argv Application argument values.
    /// \param comm The MPI communicator that should be used to construct the
    ///     \c World object.
    /// \param[in] quiet If false, will announce to \c std::cout on rank 0 when
    ///            the runtime has been initialized.
    /// \return A reference to the World constructed with \c comm.
    World& initialize(int& argc, char**& argv, const MPI_Comm& comm,
        bool quiet = false);

    /// Call this once at the very end of your main program instead of MPI_Finalize().
    void finalize();

    /// Check if the MADNESS runtime has been initialized (and not subsequently finalized).

    /// @return true if \c madness::initialize had been called more recently than \c madness::finalize, false otherwise.
    bool initialized();

    /// Check if the MADNESS runtime was initialized for quiet operation.

    /// @return true if \c madness::initialize was called with \c quiet=true .
    bool quiet();

} // namespace madness

/// @}

#endif // MADNESS_WORLD_WORLDINIT_H__INCLUDED
