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

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES

/*!
  \file examples/oep.cc
  \brief optimized effective potentials for DFT
*/


#include <chem/oep.h>
#include <madness/misc/gitinfo.h>


/// Create a specific test input for the Be atom
void write_test_input() {

	std::string filename = "test_input";

	std::ofstream of(filename);
	of << "\ngeometry\n";
	of << "  Be 0.0 0.0 0.0\n";
	of << "end\n\n\n";

    of << "dft\n";
    of << "  local canon\n";
    of << "  xc hf\n";
    of << "  maxiter 100\n";
    of << "  econv 1.0e-6\n";
    of << "  dconv 3.0e-4\n";
    of << "  L 20\n";
    of << "  k 9\n";
    of << "  no_orient 1\n";
    of << "  ncf slater 2.0\n";
    of << "end\n\n\n";

    of << "oep\n";
    of << "  model [mrks]\n";
    of << "  maxiter 2\n";
    of << "  density_threshold_high 1.0e-4\n";
    of << "  density_threshold_low 1.0e-7\n";
    of << "end\n\n\n";

    of << "plot\n";
    of << "  line x3\n";
    of << "  plane x1 x3\n";
    of << "  zoom 1\n";
    of << "  npoints 100\n";
    of << "end\n\n";

    of.close();

}


int main(int argc, char** argv) {

    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    if (world.rank() == 0) {
    	print("\n  OEP -- optimized effective potentials for DFT  \n");
    	printf("starting at time %.1f\n", wall_time());
    }
    startup(world, argc, argv,true);
    std::cout.precision(6);

    if (world.rank()==0) print(info::print_revision_information());

    commandlineparser parser(argc,argv);

    // to allow to test the program
    bool test = parser.key_exists("test");
    bool analyze = parser.key_exists("analyze");

    // create test input file if program is tested
    if (test) {
    	write_test_input();
    	parser.set_keyval("input","test_input");
    }

    // do approximate OEP calculation or test the program
    std::shared_ptr<OEP> oep(new OEP(world, parser));
    oep->print_parameters({"reference","oep","oep_calc"});
    int ierr=0;
    if (test) ierr=oep->test_oep();
    else if (analyze)  oep->analyze();
    else oep->value();

    finalize();
    return ierr;
}
