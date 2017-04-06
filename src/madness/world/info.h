/*
  This file is part of MADNESS.

  Copyright (C) 2015 Stony Brook University

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

#ifndef MADNESS_WORLD_INFO_H__INCLUDED
#define MADNESS_WORLD_INFO_H__INCLUDED

/**
 \file info.h
 \brief Defines functions that give information on this version of MADNESS.
 \ingroup configuration
*/

namespace madness {
    namespace info {

        /// Get the git commit number for this version.

        /// \return The git commit number.
        const char* git_commit();

        /// Get the MADNESS version number.

        /// \return The MADNESS version number.
        const char* version();

    } // namespace info
} // namespace madness

#endif // MADNESS_WORLD_INFO_H__INCLUDED
