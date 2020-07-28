/*
 * MolecularOrbitals_test.cpp
 *
 *  Created on: 11 Jul 2019
 *      Author: fbischoff
 */

#include<chem/MolecularOrbitals.h>
#include <chem/SCF.h>

using namespace madness;


/// will write a test input and remove it from disk upon destruction
struct write_test_input {

    double eprec=1.e-3; // was 1e-4 ... trying to make test faster

    std::string filename_;
    write_test_input(const CalculationParameters& param, const std::string& mol="lih") : filename_("test_MO_input") {
    	std::ofstream of(filename_);
        of << "dft\n";
        of << param.print_to_string(true);
        of << "end\n";

        if (mol=="lih") {
            of << "geometry\n";
            of << "eprec " << eprec << std::endl;
            of << "Li 0.0    0.0 0.0\n";
            of << "H  1.4375 0.0 0.0\n";
            of << "end\n";
        } else if (mol=="hf") {
            //double eprec=1.e-5; // trying to make test faster
            of << "geometry\n";
            of << "eprec " << eprec << std::endl;
            of << "F  0.1    0.0 0.2\n";
            of << "H  1.4375 0.0 0.0\n";
            of << "end\n";
        }
        of.close();
    }

    ~write_test_input() {
        std::remove(filename_.c_str());
    }

    std::string filename() const {return filename_;}
};


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
	SCF calc(world,test_input.filename().c_str());
        calc.set_protocol<3>(world, 1e-4);
	MolecularEnergy ME(world, calc);
	//double energy=ME.value(calc.molecule.get_all_coords().flat()); // ugh!
	ME.value(calc.molecule.get_all_coords().flat()); // ugh!

	// this has the derived parameters as well
	CalculationParameters param=calc.param;

	// read restart file (shorthand notation)
	{
		auto [amo, bmo]=MolecularOrbitals<double,3>::read_restartdata(world,calc.molecule, param.nmo_alpha(), param.nmo_beta());
	}
	// read restart file (less shorthand notation, but with autocomplete enabled)
	{
		auto mos=MolecularOrbitals<double,3>::read_restartdata(world,calc.molecule, param.nmo_alpha(), param.nmo_beta());
		MolecularOrbitals<double,3> amo=mos.first,bmo=mos.second;
	}
	// read restart file (long notation)
	std::pair<MolecularOrbitals<double,3>, MolecularOrbitals<double,3> > mos;
	mos=MolecularOrbitals<double,3>::read_restartdata(world,calc.molecule, param.nmo_alpha(), param.nmo_beta());
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
	SCF calc(world,test_input.filename().c_str());
        calc.set_protocol<3>(world, 1e-4);
	MolecularEnergy ME(world, calc);
	//double energy=ME.value(calc.molecule.get_all_coords().flat()); // ugh!
	ME.value(calc.molecule.get_all_coords().flat()); // ugh!

	// this has the derived parameters as well
	CalculationParameters param=calc.param;

	// read MOs from file and save them as AO projection
	std::vector<Function<double,3> > aos=SCF::project_ao_basis_only(world, calc.aobasis, calc.molecule);
	auto mos=MolecularOrbitals<double,3>::read_restartdata(world,calc.molecule, param.nmo_alpha(), param.nmo_beta());
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
