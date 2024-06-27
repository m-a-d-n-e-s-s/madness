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
#include <iostream>
using std::cout;
using std::endl;
#include <fstream>
#include <cstdio>

static std::ofstream fout;
namespace madness {
    void redirectio(const World& world, bool split) {
        constexpr std::size_t bufsize=256;
        char filename[bufsize];
        std::snprintf(filename,bufsize,"log.%5.5d",world.mpi.rank());
        char errfilename[bufsize];
        std::snprintf(errfilename,bufsize,"%s.%5.5d", (split ? "err" : "log"), world.mpi.rank());
        if (!freopen(   filename, "w", stdout)) MADNESS_EXCEPTION("reopening stdout failed", 0);
        if (!freopen(errfilename, "w", stderr)) MADNESS_EXCEPTION("reopening stderr failed", 0);
	      std::cout.sync_with_stdio(true);
	      std::cerr.sync_with_stdio(true);
    }
}




