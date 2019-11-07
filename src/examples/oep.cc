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


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    if (world.rank() == 0) {
    	print("\n  OEP -- optimized effective potentials for DFT  \n");
    	printf("starting at time %.1f\n", wall_time());
    }
    startup(world, argc, argv);
    std::cout.precision(6);

    const std::string input = "input";
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

    // OEP model final energy
    printf("\n   +++ starting approximate OEP iterative calculation +++\n\n");

    // read additional OEP parameters from same input file used for SCF calculation (see above)
//    std::ifstream in(input.c_str());
//    oep->read_oep_param(in);

    oep->solve_oep(HF_nemos, HF_orbens);

    finalize();
    return 0;
}
