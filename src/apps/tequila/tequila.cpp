/*
 * tequila.cpp
 *
 *  Created on: Mar. 9, 2020
 *      Author: jsk
 */

//#include <cnpy.h>
#include "/home/jsk/install/cnpy/include/cnpy.h" // integrate into cmake!
#include "NumCpp.hpp" // integrate into cmake (inluce only just pass cxx flag -I/home/jsk/install/numcpp/include/
#include <iomanip>
#include <chem/SCF.h>
#include <chem/nemo.h>
#include <chem/projector.h>
#include <nonlinsol.h>

using namespace madness;

Tensor<double> load_tensor(const std::string& name){
	cnpy::NpyArray data = cnpy::npy_load(name+".npy");
	Tensor<double> result(data.num_vals);
	auto data_ptr = data.data<double>();
	UNARYITERATOR(double, result, *_p0 = data_ptr[_i]);
	std::vector<long> shape;
	for(const auto& i:data.shape){
		shape.push_back(long(i));
	}

	result = result.reshape(shape);
	return result;
}

template<typename opT>
Tensor<double> compute_one_electron_tensor(const vecfuncT& basis, const opT& op){
	return op(basis,basis);
}

template<typename opT>
Tensor<double> compute_two_electron_tensor(World& world, const vecfuncT& basis, const opT& op){
	// <ij|kl> = <ji|lk> = <kl|ij> = <lk|ji>
	Tensor<double> g(basis.size(), basis.size(),basis.size(),basis.size());
	double energy1 = 0.0;
	double energy2 = 0.0;
	const auto& f = basis;
	for (auto i=0; i<basis.size(); i++){
		const auto im_ik = f[i]*basis;
		for (auto j=0; j<basis.size(); j++){
			const auto im_jl = apply(world, op, f[j]*basis);
			for (auto k=0; k<basis.size(); k++){
				for (auto l=0; l<basis.size(); l++){
					const auto tmp=(f[i]*f[k]).inner(op(f[j]*f[l]));
					g(i,j,k,l) = tmp;
					g(j,i,l,k) = tmp;
					g(k,l,i,j) = tmp;
					g(l,k,j,i) = tmp;
				}
			}
		}
	}
	return g;
}

double compute_energy(const Tensor<double>& hij, const Tensor<double>& gijkl, const Tensor<double>& dij, const Tensor<double>& dijkl){
	std::cout << "dij.ndim()=" << dij.flat().ndim() << "\ndijkl.ndim()" << dijkl.flat().ndim() << "\n";
	std::cout << "hij.ndim()=" << hij.flat().ndim() << "\ngijkl.ndim()" << gijkl.flat().ndim() << "\n";

	const auto one_particle_energy = hij.flat().trace(dij.flat());
	std::cout << "one particle energy = " << one_particle_energy << "\n";
	const auto two_particle_energy = gijkl.flat().trace(dijkl.flat());
	std::cout << "two particle energy = " << two_particle_energy << "\n";

	return one_particle_energy + 0.5*two_particle_energy;
}

Nemo compute_scf(World& world, const std::string& input="input"){
	const double time_scf_start = wall_time();
	std::shared_ptr<SCF> calc(new SCF(world, input));
	Nemo nemo(world, calc, input);
	nemo.get_calc()->param.print();
	const double scf_energy = nemo.value();
	if (world.rank() == 0) print("nemo energy: ", scf_energy);
	if (world.rank() == 0) printf(" at time %.1f\n", wall_time());
	const double time_scf_end = wall_time();
	return nemo;
}

vecfuncT compute_initial_guess(const Nemo& nemo){
	// replace this with guess factory of pno or cis

	// extract what we need from nemo
	const auto amo = nemo.get_calc()->amo;
	MADNESS_ASSERT(amo.size()==1);
	const auto atomic_guess = nemo.get_calc() -> ao;


	QProjector<double,3> Q(nemo.world, amo);

	auto atomic_mo = atomic_guess[0] + atomic_guess[1];
	auto atomic_vi = atomic_guess[0] - atomic_guess[1];
	atomic_mo.scale(1.0/atomic_mo.norm2());
	atomic_vi.scale(1.0/atomic_vi.norm2());

	auto vi = Q(atomic_vi);
	atomic_vi.scale(1.0/atomic_vi.norm2());

	vecfuncT basis = {amo[0], atomic_vi};
	return basis;
}

int main(int argc, char** argv) {
	initialize(argc, argv);
	World world(SafeMPI::COMM_WORLD);
	if (world.rank() == 0) printf("starting at time %.1f\n", wall_time());
	const double time_start = wall_time();
	std::cout.precision(6);

	startup(world,argc,argv,true);
	//print_meminfo(world.rank(), "startup");

	// read in the spin integrated density matrices
	auto dij = load_tensor("dij");
	auto dijkl = load_tensor("dijkl");

	std::cout << "dij\n" << dij << '\n';
	std::cout << "dijkl " << dijkl.dims()  << " " << dijkl.ndim()<< "\n" << dijkl << '\n';

	auto nemo = compute_scf(world);
	const double nuc_rep = nemo.get_calc() -> molecule.nuclear_repulsion_energy();
	auto basis = compute_initial_guess(nemo);

	auto T = Kinetic<double, 3>(world);
	auto Vnuc =  Nuclear(world, &nemo);
	auto gop = CoulombOperator(world, 1.e-6,1.e-6);
	QProjector<double,3> Q(nemo.world, nemo.get_calc()->amo);
	// will only optimize the virtuals here
	NonlinearSolver solver;
	for(auto i=0; i<100; i++){
		auto vi = basis.back();
		auto mo = basis.front();
		auto h = Vnuc(basis,basis) + T(basis,basis);
		auto g11 = gop(vi*vi);
		auto g01 = gop(mo*vi);
		auto V = Vnuc(vi)*dij(1,1) + Q(dijkl(1,1,1,1)*g11*vi + dijkl(1,1,0,0)*g01*mo) - dij(1,1)*h(0,1)*mo;
		const double lambda = dij(1,1)*T(vi,vi) + vi.inner(V) ;
		auto BSH = BSHOperator3D(world, sqrt(-2.0 * lambda/dij(1,1)), 1.e-6, 1.e-6);
		auto vi2 = BSH(-2.0/dij(1,1)*V);
		auto res = vi - vi2;
		auto err = res.norm2();
		std::cout << "err=" << err << " lambda=" << lambda << "\n";
		//vi = vi2;
		vi = solver.update(vi, res);

		vi = Q(vi);
		const auto norm = vi.norm2();
		std::cout << "norm of updated function " << norm << "\n";
		vi.scale(1.0/norm);
		basis = {mo, vi};
		if(std::fabs(err) < 1.e-3) break;
	}

	auto hij = compute_one_electron_tensor(basis, T);
	hij += compute_one_electron_tensor(basis, Vnuc);
	const auto gijkl = compute_two_electron_tensor(world, basis, gop);
	const auto energy = compute_energy(hij, gijkl, dij, dijkl);
	std::cout << "final energy = " << 2.0* energy + nuc_rep << "\n"; // factor 2 currently missing in density matrices (or because its only AB)

//	{
//		auto h = hij.flat();
//		nc::NdArray<double> hh(h.ptr(), h.size(), 1);
//		hh.tofile("hij.bin", "");
//	}{
//		auto h = dij.flat();
//		nc::NdArray<double> hh(h.ptr(), h.size(), 1);
//		hh.tofile("dij.bin", "");
//	}{
//		auto g = gijkl.flat();
//		nc::NdArray<double> gg(g.ptr(), g.size(), 1);
//		gg.tofile("gijkl.bin", "");
//	}{
//		auto g = dijkl.flat();
//		nc::NdArray<double> gg(g.ptr(), g.size(), 1);
//		gg.tofile("dijkl.bin", "");
//	}



	print_stats(world);
	finalize();
}
