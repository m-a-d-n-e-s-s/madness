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
 \file print.cc
 \brief Implementation of functions defined in print.h.
 \ingroup libraries
*/

#include <madness/world/print.h>
#include <cstring>

namespace madness {

    /// big section heading
    void print_header1(const std::string& s) {
        auto line=std::string(80,'=');
        print("");
        print(line);
        print("  ",s);
        print(line,"\n");
    };


    /// medium section heading
    void print_header2(const std::string& s) {
        auto line=std::string(80,'-');
        print("");
        print(line);
        print("  ",s);
        print(line,"\n");
    };

    /// small section heading
    void print_header3(const std::string& s) {
        auto line=std::string(80,'-');
        print("");
        print("  ",s);
        print(line,"\n");
    };

    void printf_msg_energy_time(const std::string msg, const double energy, const double time) {
        printf("%50s: %12.8f at time %8.1fs\n",msg.c_str(),energy,time);
    }


void print_justified(const char* s, int column, bool underline) {
        for (int i=0; i<column; ++i) std::cout << " ";
        std::cout << s << ENDL;
        if (underline) {
            for (int i=0; i<column; ++i) std::cout << " ";
            for (unsigned int i=0; i<std::strlen(s); ++i) std::cout << "-";
            std::cout << ENDL;
        }
    }

void print_centered(const char* s, int column, bool underline) {
    print_justified(s, column-std::strlen(s)/2, underline);
}

// Used to restore cout after redirection (cf io_redirect and io_redirect_cout)
std::streambuf* madness::io_redirect::stream_buffer_cout_default = nullptr;

} // namespace madness
