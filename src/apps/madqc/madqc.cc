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

#include "madchem.h"

using namespace madness;



int main(int argc, char** argv) {


    World& world=initialize(argc, argv,false);
    if (world.rank() == 0) {
        print_header1("MADQC -- numerical quantum chemistry in MADNESS");
    }

    startup(world,argc,argv,true);
    std::cout.precision(6);
    if (world.rank()==0) print(info::print_revision_information());


    commandlineparser parser(argc,argv);

    // create parameter classes
    // 1. read in all input blocks independently
    // 2. set up parameter logic
    // 2a from the model downstream
    // 2b from the task downstream


    // read input file
    // read into parameter handler




    // create class corresponding to qc model





    // check for the existence of the input file


    finalize();
    return 0;
}
