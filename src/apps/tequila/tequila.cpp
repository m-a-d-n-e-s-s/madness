/*
 * tequila.cpp
 *
 *  Created on: Mar. 9, 2020
 *      Author: jsk
 */

//#include <cnpy.h>
#include "/home/jsk/install/cnpy/include/cnpy.h" // integrate into cmake!
#include "NumCpp.hpp" // integrate into cmake (include only just pass cxx flag -I/home/jsk/install/numcpp/include/
#include <iomanip>
#include <chem/SCF.h>
#include <chem/nemo.h>
#include <chem/projector.h>
#include <nonlinsol.h>
#include <chem/MolecularOrbitals.h>

using namespace madness;


//template<typename T, std::size_t NDIM>
//struct allocator {
//	World& world;
//	const int n;
//
//	/// @param[in]	world	the world
//	/// @param[in]	nn		the number of functions in a given vector
//	allocator(World& world, const int nn) :
//			world(world), n(nn) {
//	}
//
//	/// allocate a vector of n empty functions
//	std::vector<Function<T, NDIM> > operator()() {
//		return zero_functions<T, NDIM>(world, n);
//	}
//};

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
	std::cout << atomic_guess.size() << " atomic guess orbitals\n";

	QProjector<double,3> Q(nemo.world, amo);

	auto atomic_mo = atomic_guess[0] + atomic_guess[1];
	auto atomic_vi = atomic_guess[0] - atomic_guess[1];
	atomic_mo.scale(1.0/atomic_mo.norm2());
	atomic_vi.scale(1.0/atomic_vi.norm2());

	auto vi = Q(atomic_vi);
	vi.scale(1.0/vi.norm2());

	auto tmp = amo[0].inner(vi);
	std::cout << "overlap of guess = " << tmp << "\n";
	vecfuncT basis = {amo[0], vi};
	return basis;
}

void print_orbital_information(World& world, const MolecularOrbitals<double, 3>& orbitals){
	if (world.rank()==0){
		std::cout << orbitals.get_mos().size() << " Molecular orbitals\n";
		std::cout << "irreps: " << orbitals.get_irreps() << "\n";
		std::cout << "occ   : " << orbitals.get_occ() << "\n";
		std::cout << "set   : " << orbitals.get_localize_sets() << "\n";
	}
}

int main(int argc, char** argv) {
	initialize(argc, argv);
	World world(SafeMPI::COMM_WORLD);
	if (world.rank() == 0) printf("starting at time %.1f\n", wall_time());
	const double time_start = wall_time();
	std::cout.precision(6);

	startup(world,argc,argv,true);
	//print_meminfo(world.rank(), "startup");

	const std::string input = "input";

	// read in the spin integrated density matrices
	auto dij = load_tensor("dij");
	auto dijkl = load_tensor("dijkl");

	std::cout << "dij\n" << dij << '\n';
	std::cout << "dijkl " << dijkl.dims()  << " " << dijkl.ndim()<< "\n" << dijkl << '\n';

	std::shared_ptr<SCF> scf_ptr(new SCF(world, input));
	Nemo nemo(world, scf_ptr, input);
	try{
		nemo.get_calc() -> load_mos(world);
		nemo.get_calc() -> ao = nemo.get_calc() -> project_ao_basis(world, nemo.get_calc() -> aobasis);
		nemo.get_calc() -> make_nuclear_potential(world);
	}catch(MadnessException& e){
		nemo.value();
	}

	std::vector<std::string> irreps;
	auto symmetry_projector = nemo.get_symmetry_projector();
	symmetry_projector(nemo.get_calc() -> amo, irreps);
	auto alphas = MolecularOrbitals<double, 3>(nemo.get_calc()->amo, nemo.get_calc()->aeps, irreps, nemo.get_calc()->aocc, nemo.get_calc()->aset);
	auto betas = MolecularOrbitals<double, 3>(nemo.get_calc()->amo, nemo.get_calc()->aeps, irreps, nemo.get_calc()->aocc, nemo.get_calc()->aset);
	const double nuc_rep = nemo.get_calc() -> molecule.nuclear_repulsion_energy();
	auto T = Kinetic<double, 3>(world);
	auto Vnuc =  Nuclear(world, nemo.get_calc().get());
	auto gop = CoulombOperator(world, 1.e-6,1.e-6);
	QProjector<double,3> Q(nemo.world, nemo.get_calc()->amo);

	// currently only closed shell
	MADNESS_ASSERT(betas.get_mos().size()==0);
	print_orbital_information(world, alphas);

	std::cout << "compute initial guess\n";
	auto basis = compute_initial_guess(nemo);
	std::cout << "done\n";

//	std::vector<std::string> virreps;
//	vecfuncT guessvirt = {basis.back()};
//	std::cout << "guessvirt.size() = " << guessvirt.size() << "\n";
//	guessvirt = symmetry_projector(guessvirt, virreps);
//	guessvirt = vecfuncT({guessvirt.back()});
//	virreps = std::vector<std::string>({"b1u"});
//	guessvirt = symmetry_projector(guessvirt, virreps);
//	auto virtuals = MolecularOrbitals<double,3>(guessvirt, Tensor<double>(1), virreps, Tensor<double>(1), std::vector<int>(1,0));
//	print_orbital_information(world, virtuals);

	std::vector<std::string> airreps;
	//basis = symmetry_projector(basis, airreps);
	auto bas = MolecularOrbitals<double, 3>(basis, Tensor<double>(basis.size()), airreps, Tensor<double>(basis.size()), std::vector<int>(basis.size(),0));
	print_orbital_information(world, bas);

	// will only optimize the virtuals here
	const double thresh = 1.e-4;
	auto size = bas.get_mos().size();
	typedef allocator<double, 3> allocT;
	typedef XNonlinearSolver<std::vector<Function<double, 3> >, double, allocT> solverT;
	allocT alloc(world, size);
	solverT solver(allocT(world, size));
	bool canonicalize = true;
	for(auto i=0; i<10; i++){
		// hardcoded for 2e example
//		auto mo = alphas.get_mos().front();
//		auto vi = virtuals.get_mos().front();
//		auto h = Vnuc(basis,basis) + T(basis,basis);
//		auto g11 = gop(vi*vi);
//		auto g01 = gop(mo*vi);
//		auto Vx = Vnuc(vi)*dij(1,1) + Q(dijkl(1,1,1,1)*g11*vi + dijkl(1,1,0,0)*g01*mo) - dij(1,1)*h(0,1)*mo;
//		const double lambdax = dij(1,1)*T(vi,vi) + vi.inner(Vx) ;
//		auto BSH = BSHOperator3D(world, sqrt(-2.0 * lambdax/dij(1,1)), 1.e-6, 1.e-6);
//		auto vi2 = BSH(-2.0/dij(1,1)*Vx);
//		auto resx = vi - vi2;
//		auto err = resx.norm2();
//		std::cout << "err=" << err << " lambda=" << lambdax << "\n";
//		//vi = vi2;
//		//vi = solver.update(vi, res);
//		vi = Q(vi);
//		const auto norm = vi.norm2();
//		std::cout << "norm of updated function " << norm << "\n";
//		vi.scale(1.0/norm);
//		Tensor<double> tmp(1);
//		tmp[0] = lambdax;
//		vecfuncT tmpvi(1,vi);
//		virtuals.update_mos_and_eps(tmpvi, tmp);
//		if(std::fabs(err) < 1.e-3) break;
//		basis = {mo, vi};

		//vecfuncT hop = T(bas.get_mos()) + Vnuc(bas.get_mos()); // not supported currently; would need for non-diagonal dij

		auto S = matrix_inner(world, bas.get_mos(), bas.get_mos());
		std::cout << "Overlap\n" << S << "\n";

		std::vector<vecfuncT> glm;

		for (auto l=0;l<size;++l){
			const auto tmp = truncate(apply(world,gop, truncate(bas.get_mos()[l]*bas.get_mos() , thresh)),thresh);
			glm.push_back(tmp);
		}

		vecfuncT V;
		Tensor<double> kinetic_part(size,size); // diagonal matrix
		for (auto a=0; a<size; ++a){
			real_function_3d Va = Vnuc(bas.get_mos()[a])*dij[a,a];
			vecfuncT moa;
			for(auto k=0;k<size; ++k){
				if(k!=a) moa.push_back(bas.get_mos()[k]);
			}
			auto Qa = QProjector<double,3>(world, moa);
			kinetic_part(a,a) = dij(a,a)*T(bas.get_mos()[a], bas.get_mos()[a]);
			for(auto k=0;k<size; ++k){
				//Va += (1.0 - a==k)*dij(a,k)*hop[k]; // assuming diagonal dij for now
				for(auto l=0;l<size; ++l){
					for(auto m=0;m<size; ++m){
						Va += 0.5*Qa((dijkl(a,k,l,m) + dijkl(k,a,l,m))*glm[l][m]*bas.get_mos()[k]);
					}
				}
			}
			V.push_back(Va);
		}

		auto Lambda = kinetic_part + matrix_inner(world, bas.get_mos(), V, true);
		std::vector<double> eps(size);
		if (canonicalize){
			Tensor<double> U, lambda;
			syev(Lambda, U, lambda);
			std::cout << "lambdas : " << lambda << "\n";
			bas.update_mos(transform(world, bas.get_mos(), U));
			V = transform(world, V, U);
			for (int i = 0;i < size;++i) {
				eps[i] = lambda(i)/dij(i,i);
			}
		}else{
			for (auto a=0; a<size; ++a){
				eps[i] = Lambda(i,i)/dij(i,i);
				for(auto k=0;k<size; ++k){
					V[a] -= Lambda(a,k)*bas.get_mos()[k];
				}
				V[a] = -2.0/dij(a,a)*V[a];
			}
			std::cout << "Lambda  :\n" << Lambda << "\n";
		}




		std::vector<poperatorT> ops(size);
		for (int i = 0;i < size;++i) {
			if (eps[i] > 0) {
				MADNESS_EXCEPTION("Positive Eigenvalue for BSH Operator?????", 1);
			}
			ops[i] = poperatorT(BSHOperatorPtr3D(world, sqrt(-2.0 * eps[i]), 1.e-6, 1.e-6));
		}

//		auto res = bas.get_mos() - apply(world, ops, V);
//		auto updated = solver.update(bas.get_mos(), res);
		auto updated = apply(world,ops,V);
		for (auto& u: updated){
			u.scale(1.0/u.norm2());
		}
		auto res = bas.get_mos() - updated;

		auto errv = norm2s(world, res);
		auto maxerr = std::max_element(std::begin(errv), std::end(errv));
		std::cout << "errv : " << errv << " | max " << *maxerr <<  "\n";
		if (*maxerr < 1.e-3) break;

		bas.update_mos(updated);


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
