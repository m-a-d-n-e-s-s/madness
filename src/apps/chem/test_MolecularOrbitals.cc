/*
 * MolecularOrbitals_test.cpp
 *
 *  Created on: 11 Jul 2019
 *      Author: fbischoff
 */

#include<chem/MolecularOrbitals.h>
#include <chem/SCF.h>
#include <chem/write_test_input.h>

using namespace madness;

int compare_calc_and_mos(World& world, const SCF& calc, const MolecularOrbitals<double,3>& amo) {
	int success=0;
	double eps_error=(amo.get_eps()-calc.aeps).normf();
	double mo_error=norm2(world,amo.get_mos()-calc.amo);
	double occ_error=(amo.get_occ()-calc.aocc).normf();
	double set_error=(amo.get_localize_sets() != calc.aset);
	print("errors", eps_error, mo_error, occ_error, set_error);

	if (eps_error>1.e-12) success++;
	if (mo_error>1.e-12) success++;
	if (occ_error>1.e-12) success++;
	if (set_error>1.e-12) success++;
	return success;
}

int test_read_restartdata(World& world) {
        //int success=0;
	CalculationParameters param1;
	param1.set_user_defined_value("maxiter",2);
	param1.set_user_defined_value("protocol",std::vector<double>({1.e-4}));

	// write restart file
	write_test_input test_input(param1);
	commandlineparser parser;
	parser.set_keyval("input",test_input.filename());
	SCF calc(world,parser);
    calc.set_protocol<3>(world, 1e-4);
	MolecularEnergy ME(world, calc);
	//double energy=ME.value(calc.molecule.get_all_coords().flat()); // ugh!
	ME.value(calc.molecule.get_all_coords().flat()); // ugh!

	// this has the derived parameters as well
	CalculationParameters param=calc.param;

	// read restart file (shorthand notation)
	{
		auto [amo, bmo]=MolecularOrbitals<double,3>::read_restartdata(world,
				param.prefix()+".restartdata",calc.molecule, param.nmo_alpha(), param.nmo_beta());
	}
	// read restart file (less shorthand notation, but with autocomplete enabled)
	{
		auto mos=MolecularOrbitals<double,3>::read_restartdata(world,
				"restartdata",calc.molecule, param.nmo_alpha(), param.nmo_beta());
		MolecularOrbitals<double,3> amo=mos.first,bmo=mos.second;
	}
	// read restart file (long notation)
	std::pair<MolecularOrbitals<double,3>, MolecularOrbitals<double,3> > mos;
	mos=MolecularOrbitals<double,3>::read_restartdata(world,
			param.prefix()+".restartdata",calc.molecule, param.nmo_alpha(), param.nmo_beta());
	MolecularOrbitals<double,3> amo,bmo;
	amo=mos.first;
	bmo=mos.second;

	return compare_calc_and_mos(world,calc,amo);
}

int test_read_restartaodata(World& world) {
	CalculationParameters param1;
	param1.set_user_defined_value("maxiter",2);
	param1.set_user_defined_value("protocol",std::vector<double>({1.e-3}));

	// write restart file
    write_test_input test_input(param1);
    commandlineparser parser;
    parser.set_keyval("input",test_input.filename());
    SCF calc(world,parser);
    calc.set_protocol<3>(world, 1e-4);
	MolecularEnergy ME(world, calc);
	//double energy=ME.value(calc.molecule.get_all_coords().flat()); // ugh!
	ME.value(calc.molecule.get_all_coords().flat()); // ugh!

	// this has the derived parameters as well
	CalculationParameters param=calc.param;

	// read MOs from file and save them as AO projection
	std::vector<Function<double,3> > aos=SCF::project_ao_basis_only(world, calc.aobasis, calc.molecule);
	auto mos=MolecularOrbitals<double,3>::read_restartdata(world,
			param.prefix()+".restartdata",calc.molecule, param.nmo_alpha(), param.nmo_beta());
	MolecularOrbitals<double,3> amo=mos.first,bmo=mos.second;
	MolecularOrbitals<double,3>::save_restartaodata(world,calc.molecule,amo,bmo,calc.aobasis);

	// read AO projections
	auto mos1=MolecularOrbitals<double,3>::read_restartaodata(world, calc.molecule, param.have_beta());
	MolecularOrbitals<double,3> amo1=mos1.first,bmo1=mos.second;

	// accept error in the MOs (by construction)
	int success=similar(amo,amo1) and similar(bmo,bmo1);
	return success;
}



int main(int argc, char** argv) {
	World& world=madness::initialize(argc, argv);
	int result=0;
	world.gop.fence();
	startup(world,argc,argv);

	result+=test_read_restartdata(world);
	result+=test_read_restartaodata(world);
	print("result",result);
	madness::finalize();
	return result;
}
