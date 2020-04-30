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
    of << "  model mrks\n";
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

    // to allow to test the program
    bool test = false;

    // parse command line arguments
    std::vector<std::string> allArgs(argv, argv + argc);
    for (auto& a : allArgs) {
        std::replace_copy(a.begin(), a.end(), a.begin(), '=', ' ');
        std::replace_copy(a.begin(), a.end(), a.begin(), '-', ' ');
        std::string key;
        std::stringstream sa(a);
        sa >> key;
        if (key == "test") test = true;
    }

    // create test input file if program is tested
    std::string input = "input";
    if (test) {
    	write_test_input();
    	input = "test_input";
    }

    std::shared_ptr<SCF> calc(new SCF(world, input.c_str())); /// see constructor in SCF.h

    if (world.rank() == 0) {
        calc->molecule.print();
//        print("\n");
//        calc->param.print("oep");
    }

	// set reference orbitals to canonical by default
    std::string arg="canon";
	calc->param.set_derived_value("localize",arg);

    std::shared_ptr<OEP> oep(new OEP(world, calc, input));

    vecfuncT HF_nemos;
    tensorT HF_orbens;

    /// TODO: find a way to save eigenvalues and implement restart options
//    const std::string saved_nemos = "HF_nemos";
//    const std::string saved_orbens = "HF_orbens";
//    std::ifstream f1(saved_nemos.c_str());
//    std::ifstream f2(saved_orbens.c_str());
//    if (f1.good() and f2.good()) { // if file HF_nemos and HF_orbens exist
//    	load_function(world, HF_nemos, saved_nemos);
//    	// load_tensor(... HF_orbens, saved_orbens ...);
//    }
//    else {
//    	const double energy = oep->value();
//    	HF_nemos = copy(world, oep->get_calc()->amo);
//    	HF_orbens = copy(oep->get_calc()->aeps);
//    	save_function(HF_nemos, saved_nemos);
//    	// save_tensor(... HF_orbens, saved_orbens ...);
//    }

    const double energy = oep->value();

    if (world.rank() == 0) {
        printf("final energy   %12.8f\n", energy);
        printf("finished at time %.1f\n", wall_time());
    }

    // save converged HF MOs and orbital energies
    HF_nemos = copy(world, oep->get_calc()->amo);
    HF_orbens = copy(oep->get_calc()->aeps);

    if (test) printf("\n   +++ starting test of the OEP program +++\n\n");
    else printf("\n   +++ starting approximate OEP iterative calculation +++\n\n");

    // read additional OEP parameters from same input file used for SCF calculation (see above)
//    std::ifstream in(input.c_str());
//    oep->read_oep_param(in);

    // do approximate OEP calculation or test the program
    if (test) oep->test_oep(HF_nemos, HF_orbens);
    else oep->solve_oep(HF_nemos, HF_orbens);

    finalize();
    return 0;
}
