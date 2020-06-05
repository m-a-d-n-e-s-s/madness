/*
 * tequila.cpp
 *
 * Madness tequila interface and
 * an empirical proove that I suck at naming things
 *
 *  Created on: Mar. 5, 2020
 *      Author: jsk
 */

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <iomanip>
#include <chem/SCF.h>
#include <chem/nemo.h>
#include <chem/projector.h>

#include "NumCpp.hpp"

using namespace madness;

// DEFINE PARAMETER TAGS FOR THE INPUT FILE
const std::string TAG_TQ = "tequila";

std::vector<Tensor<double> > solve_two_electron_cid(const vecfuncT& occ, const vecfuncT& virt, madness::Nemo& nemo){
	/// Primitive CI doubles solver for a closed shell two electron system
	/// returns the one and two electron density matrices in determinantal basis
	const double thresh=1.e-5;
	//madness::MADNESS_ASSERT(occ.size() ==1);
	//madness::MADNESS_ASSERT(virt.size() ==1);

	const auto mo = occ[0];
	const auto vi = virt[0];
	World& world=mo.world();

	madness::Nuclear V(world, &nemo);
	madness::Kinetic<double, 3> T(world);
	auto gop = std::shared_ptr < real_convolution_3d > (madness::CoulombOperatorPtr(world, 1.e-6,1.e-6));

	Tensor<double> detmat(2,2);
	detmat(0,0) = 2.0*(V(mo,mo) + T(mo,mo)) + (mo*mo).inner((*gop)(mo*mo));
	detmat(1,1) = 2.0*(V(vi,vi) + T(vi,vi)) + (vi*vi).inner((*gop)(vi*vi));
	detmat(0,1) = (mo*vi).inner((*gop)(vi*mo));
	detmat(1,0) =  detmat(0,1);

	Tensor<double> U, evals;
	syev(detmat, U, evals);
	if(world.rank()==0) std::cout << "CI Matrix:\n" << detmat << "\n";
	if(world.rank()==0) std::cout << "CI Results:\n" << evals << "\n" << U << "\n";

	// compute the reduced one particle density matrix of the ground state
	// no single exictations in this example so it is diagonal
	Tensor<double> Dij(2,2);
	Dij(0,0) = U(0,0)*U(0,0);
	Dij(1,1) = U(1,0)*U(1,0);

	if(world.rank()==0) std::cout << "One RDM:\n" << Dij.flat() << "\n";

	// two particle density matrix
	Tensor<double> Dijkl(2,2,2,2);
	Dijkl(0,0,0,0) = U(0,0)*U(0,0) ;
	Dijkl(1,1,1,1) = U(1,0)*U(1,0) ;
	Dijkl(1,1,0,0) = U(0,0)*U(1,0) ;
	Dijkl(0,0,1,1) = U(1,0)*U(0,0) ;

	if(world.rank()==0) std::cout << "Two RDM:\n" << Dijkl.flat() << "\n";


	return {Dij, Dijkl};

}

int main(int argc, char** argv) {
	initialize(argc, argv);
	World world(SafeMPI::COMM_WORLD);
	if (world.rank() == 0) printf("starting at time %.1f\n", wall_time());
	const double time_start = wall_time();
	std::cout.precision(6);

	startup(world,argc,argv,true);
	//print_meminfo(world.rank(), "startup");

	// Get the name of the input file (if given)
	const std::string input = (argc > 1) ? argv[1] : "input";

	// Compute the SCF Reference
	const double time_scf_start = wall_time();
	std::shared_ptr<SCF> calc(new SCF(world, input));
	Nemo nemo(world, calc, input);
	nemo.get_calc()->param.print();
	const double scf_energy = nemo.value();
	if (world.rank() == 0) print("nemo energy: ", scf_energy);
	if (world.rank() == 0) printf(" at time %.1f\n", wall_time());
	const double time_scf_end = wall_time();

	// extract what we need from nemo
	const auto amo = nemo.get_calc()->amo;
	const auto bmo = nemo.get_calc()->bmo;
	const auto atomic_guess = nemo.get_calc() -> ao;
	const double nuc_rep = nemo.get_calc() -> molecule.nuclear_repulsion_energy();
	madness::Nuclear V(world, &nemo);
	madness::Kinetic<double, 3> T(world);
	auto gop = std::shared_ptr < real_convolution_3d > (madness::CoulombOperatorPtr(world, 1.e-6,1.e-6));
	QProjector<double,3> Q(world, amo);

	auto atomic_mo = atomic_guess[0] + atomic_guess[1];
	auto atomic_vi = atomic_guess[0] - atomic_guess[1];
	atomic_mo.scale(1.0/atomic_mo.norm2());
	atomic_vi.scale(1.0/atomic_vi.norm2());
	auto mo = amo[0];
	const double hf_energy_atomic = 2.0*(V(atomic_mo,atomic_mo) + T(atomic_mo,atomic_mo)) + (atomic_mo*atomic_mo).inner((*gop)(atomic_mo*atomic_mo)) + nuc_rep;
	std::cout << "Atomic HF Energy = " << hf_energy_atomic << "\n";
	std::cout << "Atomic HF Energy = " << hf_energy_atomic  - nuc_rep << " (no nuclear repulsion added)\n";

	const double hf_energy = 2.0*(V(mo,mo) + T(mo,mo)) + (mo*mo).inner((*gop)(mo*mo)) + nuc_rep;
	std::cout << "HF Energy = " << hf_energy << "\n";

	auto vi = Q(atomic_vi);
	atomic_vi.scale(1.0/atomic_vi.norm2());

	Tensor<double> hij,gijkl,dij,dijkl;
	for (auto iter=0; iter<20; iter++){
		std::cout << "Iteration " << iter << "\n";

//		std::cout << "DOING ATOMIC TEST\n";
//		vecfuncT all_basis_functions = {atomic_mo, atomic_vi};
		vecfuncT all_basis_functions = {mo, vi};
		//vecfuncT all_basis_functions = {atomic_mo, atomic_vi};
		const auto S = madness::matrix_inner(world, all_basis_functions, all_basis_functions, true);

		if(world.rank() ==0 ){
			std::cout << "S:\n";
			std::cout << S << "\n";
		}

		const auto rdms = solve_two_electron_cid(vecfuncT(1,mo), vecfuncT(1,vi), nemo);
		dij = rdms[0];
		dijkl = rdms[1];

		const Tensor<double> h = V(all_basis_functions, all_basis_functions) + T(all_basis_functions, all_basis_functions);
		std::cout << "hcore\n" << h << "\n";
		Tensor<double> g(all_basis_functions.size(), all_basis_functions.size(),all_basis_functions.size(),all_basis_functions.size());
		double energy1 = 0.0;
		double energy2 = 0.0;
		const auto& f = all_basis_functions;
		for (auto i=0; i<all_basis_functions.size(); i++){
			for (auto j=0; j<all_basis_functions.size(); j++){
				energy1 += 2.0*h(i,j)*dij(i,j);
				for (auto k=0; k<all_basis_functions.size(); k++){
					for (auto l=0; l<all_basis_functions.size(); l++){
						g(i,j,k,l) = (f[i]*f[k]).inner((*gop)(f[j]*f[l]));
						energy2 += g(i,j,k,l)*dijkl(i,j,k,l);
					}
				}
			}
		}
		const double energy = energy1 + energy2;

		if(world.rank()==0){
			std::cout << "energy1 = " << energy1 << "\nenergy2 = " << energy2 << "\n";
			std::cout << "energy over rdms = " << energy << "\n";
			std::cout << "Total energy with nuclear repulsion = " << energy + nuc_rep << "\n";
		}
		auto g11 = (*gop)(vi*vi);
		auto g01 = (*gop)(mo*vi);
		// CI solver was wrong before, re-derive equations
		auto rhs = V(vi)*dij(1,1) + Q(dijkl(1,1,1,1)*g11*vi + dijkl(1,1,0,0)*g01*mo) - dij(1,1)*h(0,1)*mo;
		const double lambda = dij(1,1)*T(vi,vi) + vi.inner(rhs) ;
		std::cout << "lambda = " << lambda << "\n";
		auto BSH = BSHOperator3D(world, sqrt(-2.0 * lambda/dij(1,1)), 1.e-6, 1.e-6);
		auto vi2 = Q(BSH(-2.0/dij(1,1)*rhs));
		vi2.scale(1.0/vi2.norm2());
		auto res = vi2 - vi;
		auto err = res.norm2();
		std::cout << "err=" << err << "\n";
		vi = vi2;

		gijkl = g;
		hij = h;

		if(std::fabs(err) < 1.e-5) break;
	}

	{
		auto h = hij.flat();
		nc::NdArray<double> hh(h.ptr(), h.size(), 1);
		hh.tofile("hij.bin", "");
	}{
		auto h = dij.flat();
		nc::NdArray<double> hh(h.ptr(), h.size(), 1);
		hh.tofile("dij.bin", "");
	}{
		auto g = gijkl.flat();
		nc::NdArray<double> gg(g.ptr(), g.size(), 1);
		gg.tofile("gijkl.bin", "");
	}{
		auto g = dijkl.flat();
		nc::NdArray<double> gg(g.ptr(), g.size(), 1);
		gg.tofile("dijkl.bin", "");
	}


	print_stats(world);
	finalize();
	return 0;
}







