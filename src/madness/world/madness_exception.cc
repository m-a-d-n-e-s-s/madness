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
 \file madness_exception.cc
 \brief Implementation of \c madness::MadnessException.
 \ingroup libraries
*/

#include <madness/world/madness_exception.h>
#include <cstdlib>
#include <sstream>
#include <iostream>

namespace madness {

    std::ostream& operator<<(std::ostream& out, const MadnessException& e) {
        out << "\nMadnessException : \n\n";
        if (e.msg) out << "msg=" << e.msg << " :\n ";
        if (e.assertion) out << "assertion=" << e.assertion << " :\n ";
        out << "value=" << e.value << " : ";
        if (e.line) out << "line=" << e.line << " : ";
        if (e.function) out << "function=" << e.function << " : ";
        if (e.filename) out << "filename='" << e.filename << "'";
        out << std::endl;
        return out;
    }

    void exception_break(bool message) {
        if(message)
            std::cerr << "A madness exception occurred. Place a break point at madness::exception_break to debug.\n";
    }

} // namespace madness
