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

/// \file SCF.cc
/// \brief Molecular HF and DFT code
/// \defgroup moldft The molecular density funcitonal and Hartree-Fock code


//#define WORLD_INSTANTIATE_STATIC_TEMPLATES


#include "SCF.h"
#include <cmath>
#include <madness/mra/qmprop.h>
#include <chem/nemo.h>
#include <chem/SCFOperators.h>
#include <madness/world/worldmem.h>
#include <chem/projector.h>

namespace madness {

//    // moved to vmra.h
//    template <typename T, std::size_t NDIM>
//    DistributedMatrix<T> matrix_inner(const DistributedMatrixDistribution& d,
//                                      const std::vector< Function<T,NDIM> >& f,
//                                      const std::vector< Function<T,NDIM> >& g,
//                                      bool sym=false)


template <typename T, std::size_t NDIM>
static void verify_tree(World& world, const std::vector< Function<T,NDIM> >& v) {
	for (unsigned int i=0; i<v.size(); i++) {
		v[i].verify_tree();
	}
}


template<int NDIM>
struct unaryexp {
	void operator()(const Key<NDIM>& key, Tensor<double_complex>& t) const {
		//vzExp(t.size, t.ptr(), t.ptr());
		UNARY_OPTIMIZED_ITERATOR(double_complex, t, *_p0 = exp(*_p0););
	}
	template <typename Archive>
	void serialize(Archive& ar) {}
};




static double rsquared(const coordT& r) {
	return r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
}

// Returns exp(-I*t*V)
static Function<double_complex,3> make_exp(double t, const Function<double,3>& v) {
	v.reconstruct();
	Function<double_complex,3> expV = double_complex(0.0,-t)*v;
	expV.unaryop(unaryexp<3>());
	//expV.truncate(); expV.reconstruct();
	return expV;
}

// Timer modified to correctly nest
static bool print_timings=false;
static std::vector<double> ttt, sss;
static void START_TIMER(World& world) {
	world.gop.fence(); ttt.push_back(wall_time()); sss.push_back(cpu_time());
}

static double pop(std::vector<double>& v) {
	MADNESS_ASSERT(v.size());
	double x=v.back();
	v.pop_back();
	return x;
}
static void END_TIMER(World& world, const char* msg) {
	double wall=wall_time()-pop(ttt), cpu=cpu_time()-pop(sss);
	if (world.rank()==0 and print_timings) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, cpu, wall);
}


/// Given overlap matrix, return rotation with 3rd order error to orthonormalize the vectors
tensorT Q3(const tensorT& s) {
	tensorT Q = inner(s,s);
	Q.gaxpy(0.2,s,-2.0/3.0);
	for (int i=0; i<s.dim(0); ++i) Q(i,i) += 1.0;
	return Q.scale(15.0/8.0);
}

/// Given overlap matrix, return rotation with 2nd order error to orthonormalize the vectors
tensorT Q2(const tensorT& s) {
	tensorT Q = -0.5*s;
	for (int i=0; i<s.dim(0); ++i) Q(i,i) += 1.5;
	return Q;
}


//    SCF::SCF(World & world, const char *filename) : SCF(world, (world.rank() == 0 ? std::make_shared<std::ifstream>(filename) : nullptr)){
//    }

/// collective constructor, reads \c input on rank 0, broadcasts to all
SCF::SCF(World& world, const std::string& inputfile) : param(CalculationParameters()) {
	FunctionDefaults<3>::set_truncate_mode(1);
	PROFILE_MEMBER_FUNC(SCF);

	if (world.rank() == 0) {

		// read input parameters from the input file
		param.read(world,inputfile,"dft");

		std::ifstream ifile(inputfile);
		molecule.read(ifile);

		// set derived parameters for the molecule

		//if psp_calc is true, set all atoms to PS atoms
		//if not, check whether some atoms are PS atoms or if this a pure AE calculation
		if (param.get<bool>("psp_calc")) {
			for (size_t iatom = 0; iatom < molecule.natom(); iatom++) {
				molecule.set_pseudo_atom(iatom,true);
			}
		}

		//modify atomic charge for complete PSP calc or individual PS atoms
		for (size_t iatom = 0; iatom < molecule.natom(); iatom++) {
			if (molecule.get_pseudo_atom(iatom)){
				unsigned int an=molecule.get_atom_number(iatom);
				double zeff=get_charge_from_file("gth.xml",an);
				molecule.set_atom_charge(iatom,zeff);
			}
		}

		if (param.core_type() != "none") {
			molecule.read_core_file(param.core_type());
		}

		if(not param.no_orient()) molecule.orient();

		reset_aobasis(param.aobasis());
		param.set_derived_values(molecule,aobasis);

	}
	world.gop.broadcast_serializable(molecule, 0);
	world.gop.broadcast_serializable(param, 0);
	world.gop.broadcast_serializable(aobasis, 0);

	if (param.print_level()>2) print_timings=true;

	xc.initialize(param.xc(), !param.spin_restricted(), world, param.print_level()>1);
	//xc.plot();

	FunctionDefaults < 3 > ::set_cubic_cell(-param.L(), param.L());
	//set_protocol < 3 > (world, param.econv());
	FunctionDefaults<3>::set_truncate_mode(1);

}


void SCF::save_mos(World& world) {
	PROFILE_MEMBER_FUNC(SCF);
	archive::ParallelOutputArchive ar(world, "restartdata", param.get<int>("nio"));
	ar & current_energy & param.spin_restricted();
	ar & (unsigned int) (amo.size());
	ar & aeps & aocc & aset;
	for (unsigned int i = 0; i < amo.size(); ++i)
		ar & amo[i];
	if (!param.spin_restricted()) {
		ar & (unsigned int) (bmo.size());
		ar & beps & bocc & bset;
		for (unsigned int i = 0; i < bmo.size(); ++i)
			ar & bmo[i];
	}

	tensorT Saoamo = matrix_inner(world, ao, amo);
	tensorT Saobmo = (!param.spin_restricted()) ? matrix_inner(world, ao, bmo) : tensorT();
	if (world.rank() == 0) {
		archive::BinaryFstreamOutputArchive arao("restartaodata");
		arao << Saoamo << aeps << aocc << aset;
		if (!param.spin_restricted()) arao << Saobmo << beps << bocc << bset;
	}
}

void SCF::load_mos(World& world) {
	PROFILE_MEMBER_FUNC(SCF);
	//        const double trantol = vtol / std::min(30.0, double(param.nalpha));
	const double thresh = FunctionDefaults < 3 > ::get_thresh();
	const int k = FunctionDefaults < 3 > ::get_k();
	unsigned int nmo = 0;
	bool spinrest = false;
	amo.clear();
	bmo.clear();

	archive::ParallelInputArchive ar(world, "restartdata");

	/*
          File format:

          bool spinrestricted --> if true only alpha orbitals are present

          unsigned int nmo_alpha;
          Tensor<double> aeps;
          Tensor<double> aocc;
          vector<int> aset;
          for i from 0 to nalpha-1:
          .   Function<double,3> amo[i]

          repeat for beta if !spinrestricted

	 */

	// LOTS OF LOGIC MISSING HERE TO CHANGE OCCUPATION NO., SET,
	// EPS, SWAP, ... sigh
	ar & current_energy & spinrest;

	ar & nmo;
	MADNESS_ASSERT(nmo >= unsigned(param.nmo_alpha()));
	ar & aeps & aocc & aset;
	amo.resize(nmo);
	for (unsigned int i = 0; i < amo.size(); ++i)
		ar & amo[i];
	unsigned int n_core = molecule.n_core_orb_all();
	if (nmo > unsigned(param.nmo_alpha())) {
		aset = vector<int>(aset.begin() + n_core,
				aset.begin() + n_core + param.nmo_alpha());
		amo = vecfuncT(amo.begin() + n_core,
				amo.begin() + n_core + param.nmo_alpha());
		aeps = copy(aeps(Slice(n_core, n_core + param.nmo_alpha() - 1)));
		aocc = copy(aocc(Slice(n_core, n_core + param.nmo_alpha() - 1)));
	}

	if (amo[0].k() != k) {
		reconstruct(world, amo);
		for (unsigned int i = 0; i < amo.size(); ++i)
			amo[i] = madness::project(amo[i], k, thresh, false);
		world.gop.fence();
	}
	set_thresh(world,amo,thresh);

	//        normalize(world, amo);
	//        amo = transform(world, amo, Q3(matrix_inner(world, amo, amo)), trantol, true);
	//        truncate(world, amo);
	//        normalize(world, amo);

	if (!param.spin_restricted()) {

		if (spinrest) { // Only alpha spin orbitals were on disk
			MADNESS_ASSERT(param.nmo_alpha() >= param.nmo_beta());
			bmo.resize(param.nmo_beta());
			bset.resize(param.nmo_beta());
			beps = copy(aeps(Slice(0, param.nmo_beta() - 1)));
			bocc = copy(aocc(Slice(0, param.nmo_beta() - 1)));
			for (int i = 0; i < param.nmo_beta(); ++i)
				bmo[i] = copy(amo[i]);
		} else {
			ar & nmo;
			ar & beps & bocc & bset;

			bmo.resize(nmo);
			for (unsigned int i = 0; i < bmo.size(); ++i)
				ar & bmo[i];

			if (nmo > unsigned(param.nmo_beta())) {
				bset = vector<int>(bset.begin() + n_core,
						bset.begin() + n_core + param.nmo_beta());
				bmo = vecfuncT(bmo.begin() + n_core,
						bmo.begin() + n_core + param.nmo_beta());
				beps = copy(beps(Slice(n_core, n_core + param.nmo_beta() - 1)));
				bocc = copy(bocc(Slice(n_core, n_core + param.nmo_beta() - 1)));
			}

			if (bmo[0].k() != k) {
				reconstruct(world, bmo);
				for (unsigned int i = 0; i < bmo.size(); ++i)
					bmo[i] = madness::project(bmo[i], k, thresh, false);
				world.gop.fence();
			}
			set_thresh(world,amo,thresh);

			//                normalize(world, bmo);
			//                bmo = transform(world, bmo, Q3(matrix_inner(world, bmo, bmo)), trantol, true);
			//                truncate(world, bmo);
			//                normalize(world, bmo);

		}
	}
}

void SCF::do_plots(World& world) {
	PROFILE_MEMBER_FUNC(SCF);
	START_TIMER(world);

	std::vector<long> npt(3, static_cast<long>(param.get<int>("npt_plot")));

	if (param.plot_cell().size() == 0)
		param.plot_cell() = copy(FunctionDefaults < 3 > ::get_cell());

	if (param.get<bool>("plotdens") || param.get<bool>("plotcoul")) {
		functionT rho;
		rho = make_density(world, aocc, amo);

		if (param.spin_restricted()) {
			rho.scale(2.0);
		} else {
			functionT rhob = make_density(world, bocc, bmo);
			functionT rho_spin = rho - rhob;
			rho += rhob;
			plotdx(rho_spin, "spin_density.dx", param.plot_cell(), npt, true);

		}
		plotdx(rho, "total_density.dx", param.plot_cell(), npt, true);
		if (param.get<bool>("plotcoul")) {
			real_function_3d vnuc = potentialmanager->vnuclear();
			functionT vlocl = vnuc + apply(*coulop, rho);
			vlocl.truncate();
			vlocl.reconstruct();
			plotdx(vlocl, "coulomb.dx", param.plot_cell(), npt, true);
		}
	}

	for (int i = param.get<int>("plotlo"); i <= param.get<int>("plothi"); ++i) {
		char fname[256];
		if (i < param.nalpha()) {
			sprintf(fname, "amo-%5.5d.dx", i);
			plotdx(amo[i], fname, param.plot_cell(), npt, true);
		}
		if (!param.spin_restricted() && i < param.nbeta()) {
			sprintf(fname, "bmo-%5.5d.dx", i);
			plotdx(bmo[i], fname, param.plot_cell(), npt, true);
		}
	}
	END_TIMER(world, "plotting");
}

void SCF::project(World & world) {
	PROFILE_MEMBER_FUNC(SCF);
	reconstruct(world, amo);
	for (unsigned int i = 0; i < amo.size(); ++i) {
		amo[i] = madness::project(amo[i], FunctionDefaults < 3 > ::get_k(),
				FunctionDefaults < 3 > ::get_thresh(), false);
	}
	world.gop.fence();
	truncate(world, amo);
	normalize(world, amo);
	if (param.nbeta() && !param.spin_restricted()) {
		reconstruct(world, bmo);
		for (unsigned int i = 0; i < bmo.size(); ++i) {
			bmo[i] = madness::project(bmo[i], FunctionDefaults < 3 > ::get_k(),
					FunctionDefaults < 3 > ::get_thresh(), false);
		}
		world.gop.fence();
		truncate(world, bmo);
		normalize(world, bmo);
	}
}

void SCF::make_nuclear_potential(World & world) {
	PROFILE_MEMBER_FUNC(SCF);
	START_TIMER(world);
	potentialmanager = std::shared_ptr < PotentialManager
			> (new PotentialManager(molecule, param.core_type()));
	gthpseudopotential = std::shared_ptr<GTHPseudopotential<double>
	>(new GTHPseudopotential<double>(world, molecule));

	if (!param.pure_ae()){
		gthpseudopotential->make_pseudo_potential(world);}
	if (!param.psp_calc()) {
		potentialmanager->make_nuclear_potential(world);}
	END_TIMER(world, "Project vnuclear");
}

vecfuncT SCF::project_ao_basis(World & world, const AtomicBasisSet& aobasis) {
	PROFILE_MEMBER_FUNC(SCF);
	// Make at_to_bf, at_nbf ... map from atom to first bf on atom, and nbf/atom
	aobasis.atoms_to_bfn(molecule, at_to_bf, at_nbf);

	return SCF::project_ao_basis_only(world,aobasis, molecule);
}

vecfuncT SCF::project_ao_basis_only(World & world, const AtomicBasisSet& aobasis,
		const Molecule& molecule) {
	vecfuncT ao = vecfuncT(aobasis.nbf(molecule));
	for (int i = 0; i < aobasis.nbf(molecule); ++i) {
		functorT aofunc(new AtomicBasisFunctor(
						aobasis.get_atomic_basis_function(molecule, i)));
		ao[i] = factoryT(world).functor(aofunc).truncate_on_project().nofence().truncate_mode(1);
	}
	world.gop.fence();
	truncate(world, ao);
	normalize(world, ao);
	return ao;
}




distmatT SCF::localize_PM(World & world, const vecfuncT & mo,
		const std::vector<int> & set, const double thresh,
		const double thetamax, const bool randomize,
		const bool doprint) const {
	PROFILE_MEMBER_FUNC(SCF);
	START_TIMER(world);
	distmatT dUT = distributed_localize_PM(world, mo, ao, set, at_to_bf, at_nbf,
			thresh, thetamax, randomize, doprint);
	END_TIMER(world, "Pipek-Mezy distributed ");
	//print(UT);

	return dUT;
}

distmatT SCF::localize_new(World & world, const vecfuncT & mo,
		const std::vector<int> & set, double thresh,
		const double thetamax, const bool randomize,
		const bool doprint) const {
	// PROFILE_MEMBER_FUNC(SCF);
	START_TIMER(world);
	int nmo = mo.size();
	int nao = ao.size();

	tensorT C = matrix_inner(world, mo, ao);
	std::vector<int> at_to_bf, at_nbf; // OVERRIDE DATA IN CLASS OBJ TO USE ATOMS OR SHELLS FOR TESTING

	bool use_atomic_evecs = true;
	if (use_atomic_evecs) {
		// Transform from AOs to orthonormal atomic eigenfunctions
		int ilo = 0;
		for (size_t iat=0; iat<molecule.natom(); ++iat) {
			const tensorT& avec = aobasis.get_avec(molecule, iat);
			int ihi = ilo+avec.dim(1);
			Slice s(ilo,ihi-1);
			C(_,s) = inner(C(_,s),avec);

			// generate shell dimensions for atomic eigenfunctions
			// ... this relies upon spherical symmetry being enforced
			// when making atomic states
			const tensorT& aeps = aobasis.get_aeps(molecule, iat);
			//print(aeps);
			double prev = aeps(0L);
			int start = 0;
			int i; // used after loop
			for (i=0; i<aeps.dim(0); ++i) {
				//print(" ... ", i, prev, aeps(i), (std::abs(aeps(i)-prev) > 1e-2*std::abs(prev)));
				if (std::abs(aeps(i)-prev) > 1e-2*std::abs(prev)) {
					at_to_bf.push_back(ilo+start);
					at_nbf.push_back(i-start);
					//print("    ", start, i-start);
					start = i;
				}
				prev = aeps(i);
			}
			at_to_bf.push_back(ilo+start);
			at_nbf.push_back(i-start);
			//print("    ", start, i-start);
			ilo = ihi;
		}
		MADNESS_ASSERT(ilo==nao);
		MADNESS_ASSERT(std::accumulate(at_nbf.begin(),at_nbf.end(),0)==nao);
		MADNESS_ASSERT(at_to_bf.back()+at_nbf.back()==nao);
		//print(at_to_bf, at_nbf);
	} 
	else {
		aobasis.shells_to_bfn(molecule, at_to_bf, at_nbf);
		//aobasis.atoms_to_bfn(molecule, at_to_bf, at_nbf);
	}

	// Below here atoms may be shells or atoms --- by default shells

	int natom = at_to_bf.size();

	tensorT U(nmo, nmo);
	for (int i = 0; i < nmo; ++i) U(i, i) = 1.0;


	default_random_generator.setstate(182041+world.rank()*10101); // To help with reproducibility for debugging, etc.

	if (world.rank() == 0) {
		//MKL_Set_Num_Threads_Local(16);

		tensorT Q(nmo,natom);
		double breaksym = 1e-3;
		auto QQ = [&at_to_bf, &at_nbf,&breaksym](const tensorT& C, int i, int j, int a) -> double {
			int lo = at_to_bf[a], nbf = at_nbf[a];
			const double* Ci = &C(i,lo);
			const double* Cj = &C(j,lo);
			double qij = 0.0;
			for(int mu=0; mu<nbf; ++mu) qij += Ci[mu] * Cj[mu];
			return qij*(1.0+breaksym*a); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! break symmetry
		};

		auto makeGW = [&Q,&nmo,&natom,&QQ](const tensorT& C, double& W, tensorT& g) -> void {
			W = 0.0;
			for (int i=0; i<nmo; ++i) {
				for (int a=0; a<natom; ++a) {
					Q(i,a) = QQ(C,i,i,a);
					W += Q(i,a)*Q(i,a);
				}
			}

			for (int i = 0; i < nmo; ++i) {
				for (int j = 0; j < i; ++j) {
					double Qiiij = 0.0, Qijjj = 0.0;
					for (int a=0; a<natom; ++a) {
						double Qija = QQ(C,i,j,a);
						Qijjj += Qija*Q(j,a);
						Qiiij += Qija*Q(i,a);
					}
					g(j,i) = Qiiij - Qijjj;
					g(i,j) = - g(j,i);
				}
			}
		};

		tensorT xprev; // previous search direction
		tensorT gprev; // previous gradient
		bool rprev=true; // if true previous iteration restricted step or did incomplete search (so don't do conjugate)
		const int N = (nmo*(nmo-1))/2; // number of independent variables
		for (int iter = 0; iter < 1200; ++iter) {
			tensorT g(nmo,nmo);
			double W;

			makeGW(C,W,g);

			if (randomize && iter == 0) {
				for (int i=0; i<nmo; ++i) {
					for (int j=0; j<i; ++j) {
						g(i,j) += 0.1*(RandomValue<double>() - 0.5);
						g(j,i) = - g(i,j);
					}
				}
			}

			double maxg = g.absmax();
			if (doprint) printf("iteration %d W=%.8f maxg=%.2e\n", iter, W, maxg);
			if (maxg < thresh) break;

			// construct search direction using conjugate gradient approach
			tensorT x = copy(g);
			if (!rprev) { // Only apply conjugacy if did LS with real gradient
				double gamma = g.trace(g-gprev)/gprev.trace(gprev);
				if (doprint) print("gamma", gamma);
				x.gaxpy(1.0,xprev,gamma);
			}

			// Perform the line search.
			rprev = false;
			double dxgrad = x.trace(g)*2.0;  // 2*2 = 4 which should be prefactor on integrals in gradient
			if (dxgrad < 0 || ((iter+1)%N)==0) {
				if (doprint) print("resetting since dxgrad -ve or due to dimension", dxgrad, iter, N);
				x = copy(g);
				dxgrad = x.trace(g)*2.0;
			}
			xprev = x; // Save for next iteration
			gprev = copy(g);

			double mu = 0.01/std::max(0.1,maxg); // Restrict intial step mu by size of max gradient
			tensorT dU = matrix_exponential(x*mu);
			tensorT newC = inner(dU,C,0,0);
			double newW;
			makeGW(newC,newW,g);
			double dxgnew = x.trace(g)*2.0;

			if (randomize && iter==0) {
				rprev = true; // since did not use real gradient
			}
			else { // perform quadratic fit using f(0), df(0)/dx=dxgrad, f(mu) --- actually now use f(0), df(0)/dx, df(mu)/dx for better accuracy
				double f0 = W;
				double f1 = newW;
				//double hess = 2.0*(f1-f0-mu*dxgrad)/(mu*mu);
				double hess = (dxgnew-dxgrad)/mu; // Near convergence this is more accurate
				if (hess >= 0) {
					if (doprint) print("+ve hessian", hess);
					hess = -2.0*dxgrad; // force a bigish step to get out of bad region
					rprev = true; // since did not do line search
				}
				double mu2 = -dxgrad/hess;
				if (mu2*maxg > 0.25) {
					mu2 = 0.25/maxg; // pi/6 = 0.524, pi/4=0.785
					rprev = true; // since did not do line search
				}
				double f2p = f0 + dxgrad*mu2 + 0.5*hess*mu2*mu2;
				if (doprint) print(f0,f1,f0-f1,f2p,"dxg", dxgrad,"hess", hess, "mu", mu, "mu2", mu2);
				mu = mu2;
			}

			if (maxg < 10*thresh) {
			  breaksym = 1e-5;
			  rprev = true; // since just messed up the gradient
			}

			dU = matrix_exponential(x*mu);
			U = inner(U,dU,1,0);
			C = inner(dU,C,0,0);
		}
		bool switched = true;
		while (switched) {
			switched = false;
			for (int i = 0; i < nmo; i++) {
				for (int j = i + 1; j < nmo; j++) {
					if (set[i] == set[j]) {
						double sold = U(i, i) * U(i, i) + U(j, j) * U(j, j);
						double snew = U(i, j) * U(i, j) + U(j, i) * U(j, i);
						if (snew > sold) {
							tensorT tmp = copy(U(_, i));
							U(_, i) = U(_, j);
							U(_, j) = tmp;
							switched = true;
						}
					}
				}
			}
		}

		// Fix phases.
		for (int i = 0; i < nmo; ++i) {
			if (U(i, i) < 0.0)
				U(_, i).scale(-1.0);
		}
		//MKL_Set_Num_Threads_Local(1);
	}
	//done:
	world.gop.broadcast(U.ptr(), U.size(), 0);

	DistributedMatrix<double> dUT = column_distributed_matrix<double>(world, nmo, nmo);
	dUT.copy_from_replicated(transpose(U));

	// distmatT dUT = distributed_localize_PM(world, mo, ao, set, at_to_bf, at_nbf,
	//                                        thresh, thetamax, randomize, doprint);
	//print(UT);
	END_TIMER(world, "Pipek-Mezy new ");
	return dUT;
}

void SCF::analyze_vectors(World& world, const vecfuncT & mo, const tensorT& occ,
		const tensorT& energy, const std::vector<int>& set) {
	START_TIMER(world);
	PROFILE_MEMBER_FUNC(SCF);
	tensorT Saomo = matrix_inner(world, ao, mo);
	tensorT Saoao = matrix_inner(world, ao, ao, true);
	int nmo1 = mo.size();
	tensorT rsq, dip(3, nmo1);
	{
		functionT frsq = factoryT(world).f(rsquared).initial_level(4);
		rsq = inner(world, mo, mul_sparse(world, frsq, mo, vtol));
		for (int axis = 0; axis < 3; ++axis) {
			functionT fdip = factoryT(world).functor(
					functorT(new DipoleFunctor(axis))).initial_level(4);
			dip(axis, _) = inner(world, mo, mul_sparse(world, fdip, mo, vtol));
			for (int i = 0; i < nmo1; ++i)
				rsq(i) -= dip(axis, i) * dip(axis, i);

		}
	}
	tensorT C;
	END_TIMER(world, "Analyze vectors");

	START_TIMER(world);
	gesvp(world, Saoao, Saomo, C);
	END_TIMER(world, "Compute eigen gesv analyze vectors");
	C = transpose(C);
	long nmo = mo.size();
	size_t ncoeff = 0;
	for (long i = 0; i < nmo; ++i) {
		size_t ncoeffi = mo[i].size();
		ncoeff += ncoeffi;
		if (world.rank() == 0 and (param.print_level()>1)) {
			printf("  MO%4ld : ", i);
			if (set.size())
				printf("set=%d : ", set[i]);

			if (occ.size())
				printf("occ=%.2f : ", occ(i));

			if (energy.size())
				printf("energy=%13.8f : ", energy(i));

			printf("ncoeff=%.2e:",(double) ncoeffi);

			printf("center=(%.2f,%.2f,%.2f) : radius=%.2f\n", dip(0, i),
					dip(1, i), dip(2, i), sqrt(rsq(i)));
			aobasis.print_anal(molecule, C(i, _));
			printf("total number of coefficients = %.8e\n\n", double(ncoeff));
		}
	}
}

distmatT SCF::localize_boys(World & world, const vecfuncT & mo,
		const std::vector<int> & set, double thresh,
		const double thetamax, const bool randomize, const bool doprint) const {
	START_TIMER(world);
	long nmo = mo.size();
	tensorT dip(nmo, nmo, 3);
	for (int axis = 0; axis < 3; ++axis) {
		functionT fdip = factoryT(world).functor(functorT(new DipoleFunctor(axis))).initial_level(4);
		dip(_, _, axis) = matrix_inner(world, mo, mul_sparse(world, fdip, mo, vtol), true);
	}
	//print("dip\n", dip);
	//print("tolloc", thresh, "thetamax", thetamax);
	if (thresh < 1e-6) thresh = 1e-6; //<<<<<<<<<<<<<<<<<<<<< need to implement new line search like in pm routine
	tensorT U(nmo, nmo);
	default_random_generator.setstate(182041+world.rank()*10101); // To help with reproducibility for debugging, etc.
	if (world.rank() == 0) {
		for (long i = 0; i < nmo; ++i)
			U(i, i) = 1.0;

		tensorT xprev; // previous search direction
		tensorT gprev; // previous gradient
		bool rprev=true; // if true previous iteration restricted step or did incomplete search (so don't do conjugate)
		const int N = (nmo*(nmo-1))/2;
		for (long iter = 0; iter < 1200; ++iter) {
			tensorT g(nmo,nmo);
			double W = 0.0;
			// cannot restrict size of individual gradients if want to do line search --- should instead modify line search direction
			for (long i = 0; i < nmo; ++i) {
				W += DIP(dip, i, i, i, i);
				for (long j = 0; j < i; ++j) {
					g(j,i) = (DIP(dip, i, i, i, j) - DIP(dip, j, j, j, i));
					if (randomize && iter == 0) g(j,i) += 0.1*(RandomValue<double>() - 0.5);
					g(i,j) = - g(j,i);
				}
			}
			double maxg = g.absmax();
			if (doprint)
				printf("iteration %ld W=%.8f maxg=%.2e\n", iter, W, maxg);
			if (maxg < thresh) break;

			// construct search direction using conjugate gradient approach
			tensorT x = copy(g);
			if (!rprev) { // Only apply conjugacy if did LS with real gradient
				double gamma = g.trace(g-gprev)/gprev.trace(gprev);
				if (doprint) print("gamma", gamma);
				x.gaxpy(1.0,xprev,gamma);
			}

			// Perform the line search.
			rprev = false;
			double dxgrad = x.trace(g)*2.0;
			if (dxgrad < 0 || ((iter+1)%N)==0) {
				if (doprint) print("resetting since dxgrad -ve or due to dimension", dxgrad, iter, N);
				x = copy(g);
				dxgrad = x.trace(g)*2.0; // 2*2 = 4 which should be prefactor on integrals in gradient
			}
			xprev = x; // Save for next iteration, noting shallow copy
			gprev = g;

			double mu = 0.01/std::max(0.1,maxg); // Restrict intial step mu by size of max gradient
			tensorT dU = matrix_exponential(x*mu);
			tensorT newdip = inner(dU,dip,0,1); // can optimize this since only want (ii|ii)
			newdip = inner(dU,newdip,0,1);
			double newW = 0.0;
			for (long i = 0; i < nmo; ++i) {
				newW += DIP(newdip, i, i, i, i);
			}

			if (randomize && iter==0) {
				rprev = true; // since did not use real gradient
			}
			else { // perform quadratic fit using f(0), df(0)/dx=dxgrad, f(mu)
				double f0 = W;
				double f1 = newW;
				double hess = 2.0*(f1-f0-mu*dxgrad)/(mu*mu);
				if (hess >= 0) {
					if (doprint) print("+ve hessian", hess);
					hess = -2.0*dxgrad; // force a bigish step to get out of bad region
					rprev = true; // since did not do line search
				}
				double mu2 = -dxgrad/hess;
				if (mu2*maxg > 0.5) {
					mu2 = 0.5/maxg; // pi/6 = 0.524, pi/4=0.785
					rprev = true; // since did not do line search
				}
				double f2p = f0 + dxgrad*mu2 + 0.5*hess*mu2*mu2;
				if (doprint) print(f0,f1,f2p,"dxg", dxgrad,"hess", hess, "mu", mu, "mu2", mu2);
				mu = mu2;
			}

			dU = matrix_exponential(x*mu);
			U = inner(U,dU,1,0);
			dip = inner(dU,dip,0,1);
			dip = inner(dU,dip,0,1);
		}

		bool switched = true;
		while (switched) {
			switched = false;
			for (int i = 0; i < nmo; i++) {
				for (int j = i + 1; j < nmo; j++) {
					if (set[i] == set[j]) {
						double sold = U(i, i) * U(i, i) + U(j, j) * U(j, j);
						double snew = U(i, j) * U(i, j) + U(j, i) * U(j, i);
						if (snew > sold) {
							tensorT tmp = copy(U(_, i));
							U(_, i) = U(_, j);
							U(_, j) = tmp;
							switched = true;
						}
					}
				}
			}
		}

		// Fix phases.
		for (long i = 0; i < nmo; ++i) {
			if (U(i, i) < 0.0)
				U(_, i).scale(-1.0);
		}

	}

	world.gop.broadcast(U.ptr(), U.size(), 0);

	DistributedMatrix<double> dUT = column_distributed_matrix<double>(world, nmo, nmo);
	dUT.copy_from_replicated(transpose(U));

	END_TIMER(world, "Boys localize");
	return dUT;
}

// this version is faster than the previous version on BG/Q
distmatT SCF::kinetic_energy_matrix(World & world, const vecfuncT & v) const {
	PROFILE_MEMBER_FUNC(SCF);
	int n = v.size();
	distmatT r = column_distributed_matrix<double>(world, n, n);
	START_TIMER(world);
	reconstruct(world, v);
	END_TIMER(world, "KEmat reconstruct");
	START_TIMER(world);
	vecfuncT dvx = apply(world, *(gradop[0]), v, false);
	vecfuncT dvy = apply(world, *(gradop[1]), v, false);
	vecfuncT dvz = apply(world, *(gradop[2]), v, false);
	world.gop.fence();
	END_TIMER(world, "KEmat differentiate");
	START_TIMER(world);
	compress(world,dvx,false);
	compress(world,dvy,false);
	compress(world,dvz,false);
	world.gop.fence();
	END_TIMER(world, "KEmat compress");
	START_TIMER(world);
	r += matrix_inner(r.distribution(), dvx, dvx, true);
	r += matrix_inner(r.distribution(), dvy, dvy, true);
	r += matrix_inner(r.distribution(), dvz, dvz, true);
	END_TIMER(world, "KEmat inner products");
	r *= 0.5;
	//tensorT p(v.size(),v.size());
	//r.copy_to_replicated(p);
	return r;
}

distmatT SCF::kinetic_energy_matrix(World & world, const vecfuncT & vbra, const vecfuncT & vket) const {
	PROFILE_MEMBER_FUNC(SCF);
	MADNESS_ASSERT(vbra.size() == vket.size());
	int n = vbra.size();
	distmatT r = column_distributed_matrix<double>(world, n, n);
	reconstruct(world, vbra);
	reconstruct(world, vket);
	vecfuncT dvx_bra = apply(world, *(gradop[0]), vbra, false);
	vecfuncT dvy_bra = apply(world, *(gradop[1]), vbra, false);
	vecfuncT dvz_bra = apply(world, *(gradop[2]), vbra, false);
	vecfuncT dvx_ket = apply(world, *(gradop[0]), vket, false);
	vecfuncT dvy_ket = apply(world, *(gradop[1]), vket, false);
	vecfuncT dvz_ket = apply(world, *(gradop[2]), vket, false);
	world.gop.fence();
	compress(world,dvx_bra,false);
	compress(world,dvy_bra,false);
	compress(world,dvz_bra,false);
	compress(world,dvx_ket,false);
	compress(world,dvy_ket,false);
	compress(world,dvz_ket,false);
	world.gop.fence();
	r += matrix_inner(r.distribution(), dvx_bra, dvx_ket, true);
	r += matrix_inner(r.distribution(), dvy_bra, dvy_ket, true);
	r += matrix_inner(r.distribution(), dvz_bra, dvz_ket, true);
	r *= 0.5;
	return r;
}

vecfuncT SCF::core_projection(World & world, const vecfuncT & psi,
		const bool include_Bc) {
	PROFILE_MEMBER_FUNC(SCF);
	int npsi = psi.size();
	if (npsi == 0)
		return psi;
	size_t natom = molecule.natom();
	vecfuncT proj = zero_functions_compressed<double, 3>(world, npsi);
	tensorT overlap_sum(static_cast<long>(npsi));

	for (size_t i = 0; i < natom; ++i) {
		Atom at = molecule.get_atom(i);
		unsigned int atn = at.atomic_number;
		unsigned int nshell = molecule.n_core_orb(atn);
		if (nshell == 0)
			continue;
		for (unsigned int c = 0; c < nshell; ++c) {
			unsigned int l = molecule.get_core_l(atn, c);
			int max_m = (l + 1) * (l + 2) / 2;
			nshell -= max_m - 1;
			for (int m = 0; m < max_m; ++m) {
				functionT core = factoryT(world).functor(
						functorT(new CoreOrbitalFunctor(molecule, i, c, m)));
				tensorT overlap = inner(world, core, psi);
				overlap_sum += overlap;
				for (int j = 0; j < npsi; ++j) {
					if (include_Bc)
						overlap[j] *= molecule.get_core_bc(atn, c);
					proj[j] += core.scale(overlap[j]);
				}
			}
		}
		world.gop.fence();
	}
	if (world.rank() == 0 and param.print_level()>3)
		print("sum_k <core_k|psi_i>:", overlap_sum);
	return proj;
}

double SCF::core_projector_derivative(World & world, const vecfuncT & mo,
		const tensorT & occ, int atom, int axis) {
	PROFILE_MEMBER_FUNC(SCF);
	vecfuncT cores, dcores;
	std::vector<double> bc;
	unsigned int atn = molecule.get_atom(atom).atomic_number;
	unsigned int ncore = molecule.n_core_orb(atn);

	// projecting core & d/dx core
	for (unsigned int c = 0; c < ncore; ++c) {
		unsigned int l = molecule.get_core_l(atn, c);
		int max_m = (l + 1) * (l + 2) / 2;
		for (int m = 0; m < max_m; ++m) {
			functorT func = functorT(
					new CoreOrbitalFunctor(molecule, atom, c, m));
			cores.push_back(
					functionT(
							factoryT(world).functor(func).truncate_on_project()));
			func = functorT(
					new CoreOrbitalDerivativeFunctor(molecule, atom, axis, c,
							m));
			dcores.push_back(
					functionT(
							factoryT(world).functor(func).truncate_on_project()));
			bc.push_back(molecule.get_core_bc(atn, c));
		}
	}

	// calc \sum_i occ_i <psi_i|(\sum_c Bc d/dx |core><core|)|psi_i>
	double r = 0.0;
	for (unsigned int c = 0; c < cores.size(); ++c) {
		double rcore = 0.0;
		tensorT rcores = inner(world, cores[c], mo);
		tensorT rdcores = inner(world, dcores[c], mo);
		for (unsigned int i = 0; i < mo.size(); ++i) {
			rcore += rdcores[i] * rcores[i] * occ[i];
		}
		r += 2.0 * bc[c] * rcore;
	}

	return r;
}

bool SCF::restart_aos(World& world) {
	tensorT Saoamo, Saobmo;
	bool OK = true;
	if (world.rank() == 0) {
		try {
			archive::BinaryFstreamInputArchive arao("restartaodata");
			arao >> Saoamo >> aeps >> aocc >> aset;
			if (Saoamo.dim(0) != int(ao.size()) || Saoamo.dim(1) != param.nmo_alpha()) {
				print(" AO alpha restart data size mismatch --- starting from atomic guess instead", Saoamo.dim(0), ao.size(), Saoamo.dim(1), param.nmo_alpha());
				OK = false;
			}
			if (!param.spin_restricted()) {
				arao >> Saobmo >> beps >> bocc >> bset;
				if (Saobmo.dim(0) != int(ao.size()) || Saobmo.dim(1) != param.nmo_beta()) {
					print(" AO beta restart data size mismatch --- starting from atomic guess instead", Saobmo.dim(0), ao.size(), Saobmo.dim(1), param.nmo_beta());
					OK = false;
				}
			}
			print("\nRestarting from AO projections on disk\n");
		}
		catch (...) {
			print("\nAO restart file open/reading failed --- starting from atomic guess instead\n");
			OK=false;
		}
	}
	int fred = OK;
	world.gop.broadcast(fred, 0);
	OK = fred;
	if (!OK) return false;

	world.gop.broadcast_serializable(Saoamo, 0);
	if (!param.spin_restricted()) world.gop.broadcast_serializable(Saobmo, 0);

	tensorT S = matrix_inner(world, ao, ao), c;

	gesvp(world, S, Saoamo, c);
	amo = transform(world, ao, c, vtol, true);
	truncate(world, amo);
	orthonormalize(world, amo, param.nalpha());

	if (!param.spin_restricted()) {
		gesvp(world, S, Saobmo, c);
		bmo = transform(world, ao, c, vtol, true);
		truncate(world, bmo);
		orthonormalize(world, bmo, param.nbeta());
	}

	return true;
}

void SCF::initial_guess(World & world) {
	PROFILE_MEMBER_FUNC(SCF);
	START_TIMER(world);
	if (param.restart()) {
		load_mos(world);
	} else {


		// recalculate initial guess density matrix without core orbitals
		if (!param.pure_ae()){
			for (size_t iatom = 0; iatom < molecule.natom(); iatom++) {
				if (molecule.get_pseudo_atom(iatom)){
					double zeff=molecule.get_atom_charge(iatom);
					int atn=molecule.get_atom_number(iatom);
					aobasis.modify_dmat_psp(atn,zeff);
				}
			}
		}

		// Use the initial density and potential to generate a better process map
		functionT rho =
				factoryT(world).functor(
						functorT(
								new MolecularGuessDensityFunctor(molecule,
										aobasis))).truncate_on_project();
		double nel = rho.trace();
		if (world.rank() == 0 and param.print_level()>3)
			print("guess dens trace", nel);
		END_TIMER(world, "guess density");
		rho.scale(std::round(nel)/nel);

		if (world.size() > 1) {
			START_TIMER(world);
			LoadBalanceDeux < 3 > lb(world);
			real_function_3d vnuc;
			if (param.psp_calc()){
				vnuc = gthpseudopotential->vlocalpot();}
			else if (param.pure_ae()){
				vnuc = potentialmanager->vnuclear();}
			else {
				vnuc = potentialmanager->vnuclear();
				vnuc = vnuc + gthpseudopotential->vlocalpot();}

			lb.add_tree(vnuc,
					lbcost<double, 3>(param.vnucextra() * 1.0, param.vnucextra() * 8.0), false);
			lb.add_tree(rho, lbcost<double, 3>(1.0, 8.0), true);

			FunctionDefaults < 3 > ::redistribute(world, lb.load_balance(param.get<int>("loadbalparts")));
			END_TIMER(world, "guess loadbal");
		}

		// Diag approximate fock matrix to get initial mos
		functionT vlocal;
		if (param.nalpha() + param.nbeta() > 1) {
			START_TIMER(world);
			real_function_3d vnuc;
			if (param.psp_calc()){
				vnuc = gthpseudopotential->vlocalpot();}
			else if (param.pure_ae()){
				vnuc = potentialmanager->vnuclear();}
			else {
				vnuc = potentialmanager->vnuclear();
				vnuc = vnuc + gthpseudopotential->vlocalpot();}
			vlocal = vnuc + apply(*coulop, rho);
			END_TIMER(world, "guess Coulomb potn");
			START_TIMER(world);
			vlocal = vlocal + make_lda_potential(world, rho);
			vlocal.truncate();
			END_TIMER(world, "guess lda potn");
		} else {
			real_function_3d vnuc;
			if (param.psp_calc()){
				vnuc = gthpseudopotential->vlocalpot();}
			else if (param.pure_ae()){
				vnuc = potentialmanager->vnuclear();}
			else {
				vnuc = potentialmanager->vnuclear();
				vnuc = vnuc + gthpseudopotential->vlocalpot();}
			vlocal = vnuc;
		}
		rho.clear();
		vlocal.reconstruct();
		if (world.size() > 1) {
			START_TIMER(world);
			LoadBalanceDeux < 3 > lb(world);
			real_function_3d vnuc;
			if (param.psp_calc()){
				vnuc = gthpseudopotential->vlocalpot();}
			else if (param.pure_ae()){
				vnuc = potentialmanager->vnuclear();}
			else {
				vnuc = potentialmanager->vnuclear();
				vnuc = vnuc + gthpseudopotential->vlocalpot();}
			lb.add_tree(vnuc,
					lbcost<double, 3>(param.vnucextra() * 1.0, param.vnucextra() * 8.0), false);
			for (unsigned int i = 0; i < ao.size(); ++i) {
				lb.add_tree(ao[i], lbcost<double, 3>(1.0, 8.0), false);
			}
			FunctionDefaults < 3 > ::redistribute(world, lb.load_balance(param.get<int>("loadbalparts")));
			END_TIMER(world, "guess loadbal");
		}
		START_TIMER(world);
		tensorT overlap = matrix_inner(world, ao, ao, true);
		END_TIMER(world, "guess overlap");
		START_TIMER(world);

		tensorT kinetic(ao.size(),ao.size());
		{
			distmatT dkinetic = kinetic_energy_matrix(world, ao);
			dkinetic.copy_to_replicated(kinetic);
		}
		END_TIMER(world, "guess Kinet potn");

		START_TIMER(world);
		reconstruct(world, ao);
		vlocal.reconstruct();
		vecfuncT vpsi;

		//debug plots:
		/*{
                int npt=1001;
                functionT rhotmp =
                    factoryT(world).functor(
                                            functorT(
                                                     new MolecularGuessDensityFunctor(molecule,
                                                                                      aobasis))).truncate_on_project();
                functionT vlda=make_lda_potential(world, rhotmp);
                functionT coul=apply(*coulop, rhotmp);
                plot_line("vlocal.dat",npt, {0.0,0.0,-50.0}, {0.0,0.0,50.0}, vlocal);
                plot_line("vcoul.dat",npt, {0.0,0.0,-50.0}, {0.0,0.0,50.0}, vcoul);    
                plot_line("vlda.dat",npt, {0.0,0.0,-50.0}, {0.0,0.0,50.0}, vlda);
                plot_line("dens.dat",npt, {0.0,0.0,-50.0}, {0.0,0.0,50.0}, rhotmp);

                if (!param.pure_ae && !param.psp_calc){
                    real_function_3d vloc_ae;
                    vloc_ae = potentialmanager->vnuclear();
                    vloc_ae.reconstruct();
                    plot_line("vlocal_ae.dat",npt, {0.0,0.0,-50.0}, {0.0,0.0,50.0}, vloc_ae);
                    real_function_3d vloc_psp;
                    vloc_psp = gthpseudopotential->vlocalpot();
                    vloc_psp.reconstruct();
                    plot_line("vlocal_psp.dat",npt, {0.0,0.0,-50.0}, {0.0,0.0,50.0}, vloc_psp);
                }
            }*/

		//vlocal treated in psp includes psp and ae contribution so don't need separate clause for mixed psp/AE
		if (!param.pure_ae()) {
			double enl;
			tensorT occ = tensorT(ao.size());
			for (int i = 0;i < param.nalpha();++i) {
				occ[i] = 1.0;
			}
			for (int i = param.nalpha();size_t(i) < ao.size();++i) {
				occ[i] = 0.0;
			}
			vpsi = gthpseudopotential->apply_potential(world, vlocal, ao, occ, enl);
		}
		else {
			vpsi = mul_sparse(world, vlocal, ao, vtol);
		}

		compress(world, vpsi);
		truncate(world, vpsi);
		compress(world, ao);
		tensorT potential = matrix_inner(world, vpsi, ao, true);
		vpsi.clear();
		tensorT fock = kinetic + potential;
		fock = 0.5 * (fock + transpose(fock));
		tensorT c, e;

		//debug printing
		/*double ep = 0.0;
            double ek = 0.0;
            for(int i = 0;i < ao.size();++i){
                ep += potential(i, i);
                ek += kinetic(i, i);
                std::cout << "pot/kin " << i << "  " << potential(i,i) << "  "<< kinetic(i,i) << std::endl;
            }

            if(world.rank() == 0){
                printf("\n              epot, ekin, efock %16.8f  %16.8f  %16.8f\n", ek, ep, ek+ep);
		 */

		END_TIMER(world, "guess fock");

		START_TIMER(world);
		sygvp(world, fock, overlap, 1, c, e);
		END_TIMER(world, "guess eigen sol");
		print_meminfo(world.rank(), "guess eigen sol");

		// NAR 7/5/2013
		// commented out because it generated a lot of output
		// if(world.rank() == 0 && 0){
		//   print("initial eigenvalues");
		//   print(e);
		//   print("\n\nWSTHORNTON: initial eigenvectors");
		//   print(c);
		// }

		START_TIMER(world);
		compress(world, ao);

		unsigned int ncore = 0;
		if (param.core_type() != "none") {
			ncore = molecule.n_core_orb_all();
		}

		amo = transform(world, ao, c(_, Slice(ncore, ncore + param.nmo_alpha() - 1)), vtol, true);
		truncate(world, amo);
		normalize(world, amo);
		aeps = e(Slice(ncore, ncore + param.nmo_alpha() - 1));

		aocc = tensorT(param.nmo_alpha());
		for (int i = 0; i < param.nalpha(); ++i)
			aocc[i] = 1.0;

		if (world.rank()==0 and param.print_level()>3) print("grouping alpha orbitals into sets");
		aset=group_orbital_sets(world,aeps,aocc,param.nmo_alpha());

		if (param.nbeta() && !param.spin_restricted()) {
			bmo = transform(world, ao, c(_, Slice(ncore, ncore + param.nmo_beta() - 1)), vtol, true);
			truncate(world, bmo);
			normalize(world, bmo);
			beps = e(Slice(ncore, ncore + param.nmo_beta() - 1));
			bocc = tensorT(param.nmo_beta());
			for (int i = 0; i < param.nbeta(); ++i)
				bocc[i] = 1.0;

			if (world.rank()==0 and param.print_level()>3) print("grouping beta orbitals into sets");
			bset=group_orbital_sets(world,beps,bocc,param.nmo_beta());

		}
		END_TIMER(world, "guess orbital grouping");
	}
}

/// group orbitals into sets of similar orbital energies for localization

/// @param[in]	eps	orbital energies
/// @param[in]	occ	occupation numbers
/// @param[in]	nmo number of MOs for the given spin
/// @return		vector of length nmo with the set index for each MO
std::vector<int> SCF::group_orbital_sets(World& world, const tensorT& eps,
		const tensorT& occ, const int nmo) const {
	PROFILE_MEMBER_FUNC(SCF);

	std::vector<int> set = std::vector<int>(static_cast<size_t>(nmo), 0);
	for (int i = 1; i < nmo; ++i) {
		set[i] = set[i - 1];
		// Only the new/boys localizers can tolerate not separating out the core orbitals
		if (param.localize_pm() && (eps[i] - eps[i - 1] > 1.5 || occ[i] != 1.0)) ++(set[i]);
	}

	// pretty print out
	int lo=0;
	int iset=0;
	for (size_t i=0; i<set.size(); ++i) {
		if (iset!=set[i]) {
			if (world.rank()==0 and (param.print_level()>3)) print("set ",iset++,"  ",lo," - ", i-1);
			lo=i;
		}
	}
	if (world.rank()==0 and (param.print_level()>3)) print("set ",iset,"  ",lo," - ", nmo-1);
	return set;
}


void SCF::initial_load_bal(World & world) {
	PROFILE_MEMBER_FUNC(SCF);
	LoadBalanceDeux < 3 > lb(world);
	real_function_3d vnuc;
	if (param.psp_calc()){
		vnuc = gthpseudopotential->vlocalpot();}
	else if (param.pure_ae()){
		vnuc = potentialmanager->vnuclear();}
	else {
		vnuc = potentialmanager->vnuclear();
		vnuc = vnuc + gthpseudopotential->vlocalpot();}
	lb.add_tree(vnuc, lbcost<double, 3>(param.vnucextra() * 1.0, param.vnucextra() * 8.0));

	FunctionDefaults < 3 > ::redistribute(world, lb.load_balance(param.loadbalparts()));
}

functionT SCF::make_density(World & world, const tensorT & occ,
		const vecfuncT & v) const {
	PROFILE_MEMBER_FUNC(SCF);
	vecfuncT vsq = square(world, v);
	compress(world, vsq);
	functionT rho = factoryT(world);
	rho.compress();
	for (unsigned int i = 0; i < vsq.size(); ++i) {
		if (occ[i]) rho.gaxpy(1.0, vsq[i], occ[i], false);
	}
	world.gop.fence();
	vsq.clear();
	return rho;
}

functionT SCF::make_density(World & world, const tensorT & occ,
		const cvecfuncT & v) {
	PROFILE_MEMBER_FUNC(SCF);
	reconstruct(world, v); // For max parallelism
	std::vector < functionT > vsq(v.size());
	for (unsigned int i = 0; i < v.size(); i++) {
		vsq[i] = abssq(v[i], false);
	}
	world.gop.fence();

	compress(world, vsq); // since will be using gaxpy for accumulation
	functionT rho = factoryT(world);
	rho.compress();

	for (unsigned int i = 0; i < vsq.size(); ++i) {
		if (occ[i])
			rho.gaxpy(1.0, vsq[i], occ[i], false);

	}
	world.gop.fence();
	vsq.clear();
	rho.truncate();

	return rho;
}

std::vector<poperatorT> SCF::make_bsh_operators(World& world, const tensorT& evals) const {
	PROFILE_MEMBER_FUNC(SCF);
	int nmo = evals.dim(0);
	std::vector < poperatorT > ops(nmo);
	double tol = FunctionDefaults < 3 > ::get_thresh();
	for (int i = 0; i < nmo; ++i) {
		double eps = evals(i);
		if (eps > 0) {
			if (world.rank() == 0 and (param.print_level()>3)) {
				print("bsh: warning: positive eigenvalue", i, eps);
			}
			eps = -0.1;
		}

		ops[i] = poperatorT(
				BSHOperatorPtr3D(world, sqrt(-2.0 * eps), param.lo(), tol));
	}

	return ops;
}

std::vector<poperatorT> SCF::make_gradbsh_operators(World& world,
		const tensorT& evals, const int axis) const {
	PROFILE_MEMBER_FUNC(SCF);
	int nmo = evals.dim(0);
	std::vector < poperatorT > ops(nmo);
	double tol = FunctionDefaults < 3 > ::get_thresh();
	for (int i = 0; i < nmo; ++i) {
		double eps = evals(i);
		if (eps > 0) {
			if (world.rank() == 0 and (param.print_level()>3)) {
				print("bsh: warning: positive eigenvalue", i, eps);
			}
			eps = -0.1;
		}

		ops[i] = GradBSHOperator(world, sqrt(-2.0 * eps),
				param.lo(), tol)[axis];
	}

	return ops;
}


/// apply the HF exchange on a set of orbitals

/// @param[in]  world   the world
/// @param[in]  occ     occupation numbers
/// @param[in]  psi     the orbitals in the exchange operator
/// @param[in]  f       the orbitals |i> that the operator is applied on
/// @return     a vector of orbitals  K| i>
//    vecfuncT SCF::apply_hf_exchange(World & world, const tensorT & occ,
//                                    const vecfuncT & psi, const vecfuncT & f) const {
//        PROFILE_MEMBER_FUNC(SCF);
//        const bool same = (&psi == &f);
//        int nocc = psi.size();
//        int nf = f.size();
//        double tol = FunctionDefaults < 3 > ::get_thresh(); /// Important this is consistent with Coulomb
//        vecfuncT Kf = zero_functions_compressed<double, 3>(world, nf);
//        reconstruct(world, psi);
//        norm_tree(world, psi);
//        if (!same) {
//            reconstruct(world, f);
//            norm_tree(world, f);
//        }
//
//        //         // Smaller memory algorithm ... possible 2x saving using i-j sym
//        //         for(int i=0; i<nocc; ++i){
//        //             if(occ[i] > 0.0){
//        //                 vecfuncT psif = mul_sparse(world, psi[i], f, tol); /// was vtol
//        //                 truncate(world, psif);
//        //                 psif = apply(world, *coulop, psif);
//        //                 truncate(world, psif);
//        //                 psif = mul_sparse(world, psi[i], psif, tol); /// was vtol
//        //                 gaxpy(world, 1.0, Kf, occ[i], psif);
//        //             }
//        //         }
//
//        // Larger memory algorithm ... use i-j sym if psi==f
//        vecfuncT psif;
//        for (int i = 0; i < nocc; ++i) {
//            int jtop = nf;
//            if (same)
//                jtop = i + 1;
//            for (int j = 0; j < jtop; ++j) {
//                psif.push_back(mul_sparse(psi[i], f[j], tol, false));
//            }
//        }
//
//        world.gop.fence();
//        truncate(world, psif);
//        psif = apply(world, *coulop, psif);
//        truncate(world, psif, tol);
//        reconstruct(world, psif);
//        norm_tree(world, psif);
//        vecfuncT psipsif = zero_functions<double, 3>(world, nf * nocc);
//        int ij = 0;
//        for (int i = 0; i < nocc; ++i) {
//            int jtop = nf;
//            if (same)
//                jtop = i + 1;
//            for (int j = 0; j < jtop; ++j, ++ij) {
//                psipsif[i * nf + j] = mul_sparse(psif[ij], psi[i], false);
//                if (same && i != j) {
//                    psipsif[j * nf + i] = mul_sparse(psif[ij], psi[j], false);
//                }
//            }
//        }
//        world.gop.fence();
//        psif.clear();
//        world.gop.fence();
//        compress(world, psipsif);
//        for (int i = 0; i < nocc; ++i) {
//            for (int j = 0; j < nf; ++j) {
//                Kf[j].gaxpy(1.0, psipsif[i * nf + j], occ[i], false);
//            }
//        }
//        world.gop.fence();
//        psipsif.clear();
//        world.gop.fence();
//
//        truncate(world, Kf, tol);
//        return Kf;
//    }
//
// Used only for initial guess that is always spin-restricted LDA
functionT SCF::make_lda_potential(World & world, const functionT & arho) {
	PROFILE_MEMBER_FUNC(SCF);
	functionT vlda = copy(arho);
	vlda.reconstruct();
	vlda.unaryop(xc_lda_potential());
	return vlda;
}

vecfuncT SCF::apply_potential(World & world, const tensorT & occ,
		const vecfuncT & amo,
		const functionT & vlocal, double & exc, double & enl, int ispin) {
	PROFILE_MEMBER_FUNC(SCF);
	functionT vloc = copy(vlocal);
	exc = 0.0;
	enl = 0.0;

	// compute the local DFT potential for the MOs
	if (xc.is_dft() && !(xc.hf_exchange_coefficient() == 1.0)) {
		START_TIMER(world);

		XCOperator xcoperator(world,this,ispin,param.dft_deriv());
		if (ispin==0) exc=xcoperator.compute_xc_energy();
		vloc+=xcoperator.make_xc_potential();

		END_TIMER(world, "DFT potential");
	}

	vloc.truncate();

	START_TIMER(world);
	vecfuncT Vpsi;
	if (!param.pure_ae()){
		Vpsi = gthpseudopotential->apply_potential(world, vloc, amo, occ, enl);}
	else {
		Vpsi = mul_sparse(world, vloc, amo, vtol);}

	END_TIMER(world, "V*psi");
	print_meminfo(world.rank(), "V*psi");
	if (xc.hf_exchange_coefficient()) {
		START_TIMER(world);
		//            vecfuncT Kamo = apply_hf_exchange(world, occ, amo, amo);
		Exchange<double,3> K=Exchange<double,3>(world,this,ispin).small_memory(false).same(true);
		vecfuncT Kamo=K(amo);
		tensorT excv = inner(world, Kamo, amo);
		double exchf = 0.0;
		for (unsigned long i = 0; i < amo.size(); ++i) {
			exchf -= 0.5 * excv[i] * occ[i];
		}
		if (!xc.is_spin_polarized())
			exchf *= 2.0;
		gaxpy(world, 1.0, Vpsi, -xc.hf_exchange_coefficient(), Kamo);
		Kamo.clear();
		END_TIMER(world, "HF exchange");
		exc = exchf * xc.hf_exchange_coefficient() + exc;
	}
	// need to come back to this for psp - when is this used?
	if (param.pure_ae()){
		potentialmanager->apply_nonlocal_potential(world, amo, Vpsi);}

	START_TIMER(world);
	truncate(world, Vpsi);
	END_TIMER(world, "Truncate Vpsi");
	print_meminfo(world.rank(), "Truncate Vpsi");
	world.gop.fence();
	return Vpsi;
}

tensorT SCF::derivatives(World & world, const functionT& rho) const {
	PROFILE_MEMBER_FUNC(SCF);
	START_TIMER(world);

	vecfuncT dv(molecule.natom() * 3);
	vecfuncT du = zero_functions<double, 3>(world, molecule.natom() * 3);
	tensorT rc(molecule.natom() * 3);
	for (size_t atom = 0; atom < molecule.natom(); ++atom) {
		for (int axis = 0; axis < 3; ++axis) {
			functorT func(new MolecularDerivativeFunctor(molecule, atom, axis));
			dv[atom * 3 + axis] =
					functionT(
							factoryT(world).functor(func).nofence().truncate_on_project().truncate_mode(0));
			if (param.core_type() != "none"
					&& molecule.is_potential_defined_atom(atom)) {
				// core potential contribution
				func = functorT(
						new CorePotentialDerivativeFunctor(molecule, atom,
								axis));
				du[atom * 3 + axis] = functionT(
						factoryT(world).functor(func).truncate_on_project());

				// core projector contribution
				rc[atom * 3 + axis] =
						potentialmanager->core_projector_derivative(world, amo,
								aocc, atom, axis);
				if (!param.spin_restricted()) {
					if (param.nbeta())
						rc[atom * 3 + axis] +=
								potentialmanager->core_projector_derivative(
										world, bmo, bocc, atom, axis);
				} else {
					rc[atom * 3 + axis] *= 2 * 2;
					// because of 2 electrons in each valence orbital bra+ket
				}
			}
		}
	}

	world.gop.fence();
	tensorT r = inner(world, rho, dv);
	world.gop.fence();
	tensorT ru = inner(world, rho, du);
	dv.clear();
	du.clear();
	world.gop.fence();
	tensorT ra(r.size());
	for (size_t atom = 0; atom < molecule.natom(); ++atom) {
		for (int axis = 0; axis < 3; ++axis) {
			ra[atom * 3 + axis] = molecule.nuclear_repulsion_derivative(atom,
					axis);
		}
	}
	//if (world.rank() == 0) print("derivatives:\n", r, ru, rc, ra);
	r += ra + ru + rc;
	END_TIMER(world, "derivatives");

	if (world.rank() == 0 and (param.print_level()>1)) {
		print("\n Derivatives (a.u.)\n -----------\n");
		print(
				"  atom        x            y            z          dE/dx        dE/dy        dE/dz");
		print(
				" ------ ------------ ------------ ------------ ------------ ------------ ------------");
		for (size_t i = 0; i < molecule.natom(); ++i) {
			const Atom& atom = molecule.get_atom(i);
			printf(" %5d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n", int(i),
					atom.x, atom.y, atom.z, r[i * 3 + 0], r[i * 3 + 1],
					r[i * 3 + 2]);
		}
	}
	return r;
}

void SCF::dipole_matrix_elements(World& world, const vecfuncT & mo, const tensorT& occ,
		const tensorT& energy, int spin) {
	START_TIMER(world);
	int nmo = mo.size();
	tensorT mat_el(3, nmo, nmo);
	for (int axis = 0; axis < 3; ++axis) {
		functionT fdip = factoryT(world).functor(
				functorT(new DipoleFunctor(axis)));
		mat_el(axis, _, _) = matrix_inner(world, mo, mul_sparse(world, fdip, mo, vtol), true);
	}

	double ha2ev=27.211396132;
	FILE *f=0;
	if (spin==0){
		f = fopen("mat_els_alpha.dat", "w");}
	else{
		f = fopen("mat_els_beta.dat", "w");}
	fprintf(f, "#initial | Energy (eV) | final  | Energy (eV) | Matrix el.  | Trans. E (eV)\n");
	fprintf(f, "%4i  %4i\n", nmo, nmo);
	fprintf(f, "%2i\n", 1);
	fprintf(f, "%13.8f\n", 0.0);
	for (int axis = 0; axis < 3; ++axis) {
		fprintf(f, "# Cartesian component %2i\n", axis+1);
		for (int i = 0; i < nmo; ++i) {
			for (int j = 0; j < nmo; ++j) {
				fprintf(f, "%4i\t %13.8f\t %4i\t %13.8f\t %13.8f\t %13.8f\n", i+1, energy(i)*ha2ev, j+1, energy(j)*ha2ev, mat_el(axis,i,j), (energy(j)-energy(i))*ha2ev);
			}
		}
	}
	fclose(f);
	END_TIMER(world, "Matrix elements");
}

tensorT SCF::dipole(World & world, const functionT& rho) const {
	PROFILE_MEMBER_FUNC(SCF);
	START_TIMER(world);
	tensorT mu(3);

	for (unsigned int axis = 0; axis < 3; ++axis) {
		std::vector<int> x(3ul, 0);
		x[axis] = true;
		functionT dipolefunc = factoryT(world)
                    		.functor(functorT(new MomentFunctor(x)));
		mu[axis] = -dipolefunc.inner(rho);
		mu[axis] += molecule.nuclear_dipole(axis);
	}

	if (world.rank() == 0 and (param.print_level()>1)) {
		print("\n Dipole Moment (a.u.)\n -----------\n");
		print("     x: ", mu[0]);
		print("     y: ", mu[1]);
		print("     z: ", mu[2]);
		print(" Total Dipole Moment: ", mu.normf(),"\n");
	}
	END_TIMER(world, "dipole");

	return mu;
}

void SCF::vector_stats(const std::vector<double> & v, double & rms,
		double & maxabsval) const {
	PROFILE_MEMBER_FUNC(SCF);
	rms = 0.0;
	maxabsval = v[0];
	for (unsigned int i = 0; i < v.size(); ++i) {
		rms += v[i] * v[i];
		maxabsval = std::max<double>(maxabsval, std::abs(v[i]));
	}
	rms = sqrt(rms / v.size());
}

vecfuncT SCF::compute_residual(World & world, tensorT & occ, tensorT & fock,
		const vecfuncT & psi, vecfuncT & Vpsi, double & err) {

	START_TIMER(world);
	PROFILE_MEMBER_FUNC(SCF);
	double trantol = vtol / std::min(30.0, double(psi.size()));
	int nmo = psi.size();

	tensorT eps(nmo);
	for (int i = 0; i < nmo; ++i) {
		eps(i) = std::min(-0.05, fock(i, i));
		fock(i, i) -= eps(i);
	}
	vecfuncT fpsi = transform(world, psi, fock, trantol, true);

	for (int i = 0; i < nmo; ++i) { // Undo the damage
		fock(i, i) += eps(i);
	}

	gaxpy(world, 1.0, Vpsi, -1.0, fpsi);
	fpsi.clear();
	std::vector<double> fac(nmo, -2.0);
	scale(world, Vpsi, fac);
	std::vector < poperatorT > ops = make_bsh_operators(world, eps);
	set_thresh(world, Vpsi, FunctionDefaults < 3 > ::get_thresh());
	END_TIMER(world, "Compute residual stuff");

	START_TIMER(world);
	vecfuncT new_psi = apply(world, ops, Vpsi);
	END_TIMER(world, "Apply BSH");
	ops.clear();
	Vpsi.clear();
	world.gop.fence();

	// Thought it was a bad idea to truncate *before* computing the residual
	// but simple tests suggest otherwise ... no more iterations and
	// reduced iteration time from truncating.
	START_TIMER(world);
	truncate(world, new_psi);
	END_TIMER(world, "Truncate new psi");

	START_TIMER(world);
	vecfuncT r = sub(world, psi, new_psi);
	std::vector<double> rnorm = norm2s(world, r);
	if (world.rank() == 0 and (param.print_level()>1))
		print("residuals", rnorm);
	double rms, maxval;
	vector_stats(rnorm, rms, maxval);
	err = maxval;
	if (world.rank() == 0 and (param.print_level()>1))
		print("BSH residual: rms", rms, "   max", maxval);
	END_TIMER(world, "BSH residual");
	return r;
}

tensorT SCF::make_fock_matrix(World & world, const vecfuncT & psi,
		const vecfuncT & Vpsi, const tensorT & occ, double & ekinetic) const {
	PROFILE_MEMBER_FUNC(SCF);
	START_TIMER(world);
	tensorT pe = matrix_inner(world, Vpsi, psi, true);
	END_TIMER(world, "PE matrix");

	std::shared_ptr< WorldDCPmapInterface< Key<3> > > oldpmap = FunctionDefaults<3>::get_pmap();
	vecfuncT psicopy=psi; // Functions are shallow copy so this is lightweight
	if (world.size() > 1) {
		START_TIMER(world);
		LoadBalanceDeux < 3 > lb(world);
		for (unsigned int i = 0; i < psi.size(); ++i) {
			lb.add_tree(psi[i], lbcost<double, 3>(1.0, 8.0), false);
		}
		world.gop.fence();
		END_TIMER(world, "KE compute loadbal");

		START_TIMER(world);
		std::shared_ptr< WorldDCPmapInterface< Key<3> > > newpmap = lb.load_balance(param.loadbalparts());
		FunctionDefaults<3>::set_pmap(newpmap);

		world.gop.fence();
		for (unsigned int i=0; i<psi.size(); ++i) psicopy[i] = copy(psi[i],newpmap,false);
		world.gop.fence();
		END_TIMER(world, "KE redist");
	}

	START_TIMER(world);
	tensorT ke(psi.size(),psi.size());
	{
		distmatT k = kinetic_energy_matrix(world, psicopy);
		k.copy_to_replicated(ke); // !!!!!!!! ugh
	}
	END_TIMER(world, "KE matrix");

	psicopy.clear();
	if (world.size() > 1) {
		FunctionDefaults<3>::set_pmap(oldpmap); // ! DON'T FORGET !
	}

	START_TIMER(world);
	int nocc = occ.size();
	ekinetic = 0.0;
	for (int i = 0; i < nocc; ++i) {
		ekinetic += occ[i] * ke(i, i);
	}
	ke += pe;
	pe = tensorT();
	ke.gaxpy(0.5, transpose(ke), 0.5);
	END_TIMER(world, "Make fock matrix rest");
	return ke;
}

/// Compute the two-electron integrals over the provided set of orbitals

/// Returned is a *replicated* tensor of \f$(ij|kl)\f$ with \f$i>=j\f$
/// and \f$k>=l\f$.  The symmetry \f$(ij|kl)=(kl|ij)\f$ is enforced.
Tensor<double> SCF::twoint(World& world, const vecfuncT& psi) const {
	PROFILE_MEMBER_FUNC(SCF);
	double tol = FunctionDefaults < 3 > ::get_thresh(); /// Important this is consistent with Coulomb
	reconstruct(world, psi);
	norm_tree(world, psi);

	// Efficient version would use mul_sparse vector interface
	vecfuncT pairs;
	for (unsigned int i = 0; i < psi.size(); ++i) {
		for (unsigned int j = 0; j <= i; ++j) {
			pairs.push_back(mul_sparse(psi[i], psi[j], tol, false));
		}
	}

	world.gop.fence();
	truncate(world, pairs);
	vecfuncT Vpairs = apply(world, *coulop, pairs);

	return matrix_inner(world, pairs, Vpairs, true);
}

tensorT SCF::matrix_exponential(const tensorT& A) const {
	PROFILE_MEMBER_FUNC(SCF);
	MADNESS_ASSERT(A.dim(0) == A.dim(1));

	// Power iteration to estimate the 2-norm of the matrix. Used
	// to use Frobenius or 1-norms but neither were very tight.
	double anorm;
	{
		tensorT x(A.dim(0));
		x.fillrandom(); x.scale(1.0/x.normf());
		double prev = 0.0;
		for (int i=0; i<100; i++) {
			tensorT xnew = inner(A,inner(A,x,1,0),0,0);
			anorm = std::sqrt(std::abs((x.trace(xnew))));
			double err = std::abs(prev-anorm)/anorm;
			//print(i,anorm,err,A.normf());
			if (err < 0.01) break; // just need 1-2 digits
			x = xnew.scale(1.0/xnew.normf());
			prev = anorm;
		}
	}

	// Scale A by a power of 2 until it is "small"
	int n = 0;
	double scale = 1.0;
	while (anorm * scale > 0.089) { // so that 9th order expansion is accurate to 1e-15
		++n;
		scale *= 0.5;
	}
	tensorT B = scale * A;    // B = A*2^-n

	// Make identity
	tensorT I = tensorT(2, B.dims());
	for (int i = 0; i < I.dim(0); ++i) I(i, i) = 1.0;

	// Compute exp(B) using Taylor series optimized to reduce cost --- Chebyshev is only a minor improvement
	tensorT expB;
	if (anorm > 0.24e-1) {
		tensorT B2 = inner(B,B);
		tensorT B4 = inner(B2,B2);
		tensorT B6 = inner(B4,B2);
		expB = I + inner(B,B6+42.*B4+840.*B2+5040.*I).scale(1./5040.) + inner(B2,B6+56.*B4+1680.*B2+20160.*I).scale(1./40320.);
	}
	else if (anorm > 0.26e-2) {
		tensorT B2 = inner(B,B);
		tensorT B4 = inner(B2,B2);
		expB = I + inner(B,42.*B4+840.*B2+5040.*I).scale(1./5040.) + inner(B2,56.*B4+1680.*B2+20160.*I).scale(1./40320.);
	}
	else if (anorm > 0.18e-4) {
		tensorT B2 = inner(B,B);
		expB = I + inner(B,840.*B2+5040.*I).scale(1./5040.) + inner(B2,1680.*B2+20160.*I).scale(1./40320.);
	}
	else if (anorm > 4.5e-8) {
		expB = I + B + inner(B,B).scale(0.5);
	} 
	else {
		expB = I + B;
	} 

	// // Old algorithm
	// tensorT oldexpB = copy(I);
	// const double tol = 1e-13;
	// int k = 1;
	// tensorT term = B;
	// while (term.normf() > tol) {
	//     oldexpB += term;
	//     term = inner(term, B);
	//     ++k;
	//     term.scale(1.0 / k);
	// }
	// Error check for validation
	// double err = (expB-oldexpB).normf();
	// print("matxerr", anorm, err);

	// Repeatedly square to recover exp(A)
	while (n--) expB = inner(expB, expB);

	return expB;
}

/// compute the unitary transformation that diagonalizes the fock matrix

/// @param[in]  world   the world
/// @param[in]  overlap the overlap matrix of the orbitals
/// @param[inout]       fock    the fock matrix; diagonal upon exit
/// @param[out] evals   the orbital energies
/// @param[in]  occ     the occupation numbers
/// @param[in]  thresh_degenerate       threshold for orbitals being degenerate
/// @return             the unitary matrix U: U^T F U = evals
tensorT SCF::get_fock_transformation(World& world, const tensorT& overlap,
		tensorT& fock, tensorT& evals, const tensorT& occ,
		const double thresh_degenerate) const {
	PROFILE_MEMBER_FUNC(SCF);

	START_TIMER(world);
	tensorT U;
	sygvp(world, fock, overlap, 1, U, evals);
	END_TIMER(world, "Diagonalization Fock-mat w sygv");

	long nmo = fock.dim(0);

	START_TIMER(world);
	// Within blocks with the same occupation number attempt to
	// keep orbitals in the same order (to avoid confusing the
	// non-linear solver).
	// !!!!!!!!!!!!!!!!! NEED TO RESTRICT TO OCCUPIED STATES?
	bool switched = true;
	while (switched) {
		switched = false;
		for (int i = 0; i < nmo; i++) {
			for (int j = i + 1; j < nmo; j++) {
				if (occ(i) == occ(j)) {
					double sold = U(i, i) * U(i, i) + U(j, j) * U(j, j);
					double snew = U(i, j) * U(i, j) + U(j, i) * U(j, i);
					if (snew > sold) {
						tensorT tmp = copy(U(_, i));
						U(_, i) = U(_, j);
						U(_, j) = tmp;
						std::swap(evals[i], evals[j]);
						switched = true;
					}
				}
			}
		}
	}

	// Fix phases.
	for (long i = 0; i < nmo; ++i)
		if (U(i, i) < 0.0)
			U(_, i).scale(-1.0);

	// Rotations between effectively degenerate states confound
	// the non-linear equation solver ... undo these rotations
	long ilo = 0; // first element of cluster
	while (ilo < nmo - 1) {
		long ihi = ilo;
		while (fabs(evals[ilo] - evals[ihi + 1])
				< thresh_degenerate * 10.0 * std::max(fabs(evals[ilo]), 1.0)) {
			++ihi;
			if (ihi == nmo - 1)
				break;
		}
		long nclus = ihi - ilo + 1;
		if (nclus > 1) {
			//print("   found cluster", ilo, ihi);
			tensorT q = copy(U(Slice(ilo, ihi), Slice(ilo, ihi)));
			//print(q);
			// Special code just for nclus=2
			// double c = 0.5*(q(0,0) + q(1,1));
			// double s = 0.5*(q(0,1) - q(1,0));
			// double r = sqrt(c*c + s*s);
			// c /= r;
			// s /= r;
			// q(0,0) = q(1,1) = c;
			// q(0,1) = -s;
			// q(1,0) = s;

			// Polar Decomposition
			tensorT VH(nclus, nclus);
			tensorT W(nclus, nclus);
			Tensor<double> sigma(nclus);

			svd(q, W, sigma, VH);
			q = transpose(inner(W,VH)).conj();
			U(_, Slice(ilo, ihi)) = inner(U(_, Slice(ilo, ihi)), q);
		}
		ilo = ihi + 1;
	}

	world.gop.broadcast(U.ptr(), U.size(), 0);
	world.gop.broadcast(evals.ptr(), evals.size(), 0);

	fock = 0;
	for (unsigned int i = 0; i < nmo; ++i)
		fock(i, i) = evals(i);
	return U;
}

/// diagonalize the fock matrix, taking care of degenerate states

/// Vpsi is passed in to make sure orbitals and Vpsi are in phase
/// @param[in]  world   the world
/// @param[inout]       fock    the fock matrix (diagonal upon exit)
/// @param[inout]       psi             the orbitals
/// @param[inout]       Vpsi    the orbital times the potential
/// @param[out] evals   the orbital energies
/// @param[in]  occ             occupation numbers
/// @param[in]  thresh  threshold for rotation and truncation
/// @return             the unitary matrix U: U^T F U = evals
tensorT SCF::diag_fock_matrix(World& world, tensorT& fock, vecfuncT& psi,
		vecfuncT& Vpsi, tensorT& evals, const tensorT& occ,
		const double thresh) const {
	PROFILE_MEMBER_FUNC(SCF);

	// compute the unitary transformation matrix U that diagonalizes
	// the fock matrix
	tensorT overlap = matrix_inner(world, psi, psi, true);
	tensorT U = get_fock_transformation(world, overlap, fock, evals, occ,
			thresh);

	//eliminate mixing between occ and unocc
	int nmo=U.dim(0);
	for (int i=0; i<param.nalpha(); ++i){
		//make virt orthog to occ without changing occ states
		for (int j=param.nalpha(); j<nmo; ++j){
			U(j,i)=0.0;
		}
	}

	// transform the orbitals and the orbitals times the potential
	Vpsi = transform(world, Vpsi, U, vtol / std::min(30.0, double(psi.size())),
			false);
	psi = transform(world, psi, U,
			FunctionDefaults < 3
			> ::get_thresh() / std::min(30.0, double(psi.size())),
			true);
	truncate(world, Vpsi, vtol, false);
	truncate(world, psi);
	normalize(world, psi);

	END_TIMER(world, "Diagonalization rest");
	return U;
}

void SCF::loadbal(World & world, functionT & arho, functionT & brho,
		functionT & arho_old, functionT & brho_old, subspaceT & subspace) {
	if (world.size() == 1)
		return;

	LoadBalanceDeux < 3 > lb(world);
	real_function_3d vnuc;
	if (param.psp_calc()){
		vnuc = gthpseudopotential->vlocalpot();}
	else if (param.pure_ae()){
		vnuc = potentialmanager->vnuclear();}
	else {
		vnuc = potentialmanager->vnuclear();
		vnuc = vnuc + gthpseudopotential->vlocalpot();}
	lb.add_tree(vnuc, lbcost<double, 3>(param.vnucextra() * 1.0, param.vnucextra() * 8.0),
			false);
	lb.add_tree(arho, lbcost<double, 3>(1.0, 8.0), false);
	for (unsigned int i = 0; i < amo.size(); ++i) {
		lb.add_tree(amo[i], lbcost<double, 3>(1.0, 8.0), false);
	}
	if (param.nbeta() && !param.spin_restricted()) {
		lb.add_tree(brho, lbcost<double, 3>(1.0, 8.0), false);
		for (unsigned int i = 0; i < bmo.size(); ++i) {
			lb.add_tree(bmo[i], lbcost<double, 3>(1.0, 8.0), false);
		}
	}
	world.gop.fence();

	FunctionDefaults < 3 > ::redistribute(world, lb.load_balance(param.loadbalparts())); // 6.0 needs retuning after param.vnucextra

	world.gop.fence();
}

void SCF::rotate_subspace(World& world, const tensorT& U, subspaceT& subspace,
		int lo, int nfunc, double trantol) const {
	PROFILE_MEMBER_FUNC(SCF);
	for (unsigned int iter = 0; iter < subspace.size(); ++iter) {
		vecfuncT& v = subspace[iter].first;
		vecfuncT& r = subspace[iter].second;
		vecfuncT vnew = transform(world, vecfuncT(&v[lo], &v[lo + nfunc]), U, trantol, false);
		vecfuncT rnew = transform(world, vecfuncT(&r[lo], &r[lo + nfunc]), U, trantol, false);
		world.gop.fence();
		for (int i=0; i<nfunc; i++) {
			v[i] = vnew[i];
			r[i] = rnew[i];
		}
	}
	world.gop.fence();
}

void SCF::rotate_subspace(World& world, const distmatT& dUT, subspaceT& subspace,
		int lo, int nfunc, double trantol) const {
	PROFILE_MEMBER_FUNC(SCF);
	for (unsigned int iter = 0; iter < subspace.size(); ++iter) {
		vecfuncT& v = subspace[iter].first;
		vecfuncT& r = subspace[iter].second;
		vecfuncT vnew = transform(world, vecfuncT(&v[lo], &v[lo + nfunc]), dUT, false);
		vecfuncT rnew = transform(world, vecfuncT(&r[lo], &r[lo + nfunc]), dUT, false);
		world.gop.fence();
		for (int i=0; i<nfunc; i++) {
			v[i] = vnew[i];
			r[i] = rnew[i];
		}
	}
	world.gop.fence();
}

void SCF::update_subspace(World & world, vecfuncT & Vpsia, vecfuncT & Vpsib,
		tensorT & focka, tensorT & fockb, subspaceT & subspace, tensorT & Q,
		double & bsh_residual, double & update_residual) {
	PROFILE_MEMBER_FUNC(SCF);
	double aerr = 0.0, berr = 0.0;
	vecfuncT vm = amo;

	// Orbitals with occ!=1.0 exactly must be solved for as eigenfunctions
			// so zero out off diagonal lagrange multipliers
			for (int i = 0; i < param.nmo_alpha(); i++) {
				if (aocc[i] != 1.0) {
					double tmp = focka(i, i);
					focka(i, _) = 0.0;
					focka(_, i) = 0.0;
					focka(i, i) = tmp;
				}
			}

			vecfuncT rm = compute_residual(world, aocc, focka, amo, Vpsia, aerr);
			if (param.nbeta() != 0 && !param.spin_restricted()) {
				for (int i = 0; i < param.nmo_beta(); i++) {
					if (bocc[i] != 1.0) {
						double tmp = fockb(i, i);
						fockb(i, _) = 0.0;
						fockb(_, i) = 0.0;
						fockb(i, i) = tmp;
					}
				}

				vecfuncT br = compute_residual(world, bocc, fockb, bmo, Vpsib, berr);
				vm.insert(vm.end(), bmo.begin(), bmo.end());
				rm.insert(rm.end(), br.begin(), br.end());
			}

			START_TIMER(world);
			bsh_residual = std::max(aerr, berr);
			world.gop.broadcast(bsh_residual, 0);
			compress(world, vm, false);
			compress(world, rm, false);
			world.gop.fence();

			restart:
			subspace.push_back(pairvecfuncT(vm, rm));
			int m = subspace.size();
			tensorT ms(m);
			tensorT sm(m);
			for (int s = 0; s < m; ++s) {
				const vecfuncT & vs = subspace[s].first;
				const vecfuncT & rs = subspace[s].second;
				for (unsigned int i = 0; i < vm.size(); ++i) {
					ms[s] += vm[i].inner_local(rs[i]);
					sm[s] += vs[i].inner_local(rm[i]);
				}
			}

			world.gop.sum(ms.ptr(), m);
			world.gop.sum(sm.ptr(), m);
			tensorT newQ(m, m);
			if (m > 1)
				newQ(Slice(0, -2), Slice(0, -2)) = Q;

			newQ(m - 1, _) = ms;
			newQ(_, m - 1) = sm;
			Q = newQ;
			//if (world.rank() == 0) { print("kain Q"); print(Q); }
			tensorT c;
			//if (world.rank() == 0) {
				double rcond = 1e-12;
				while (1) {
					c = KAIN(Q, rcond);
					if (world.rank() == 0 and (param.print_level()>3)) print("kain c:", c);
					//if (std::abs(c[m - 1]) < 5.0) { // was 3
						if (c.absmax() < 3.0) { // was 3
							break;
						} else if (rcond < 0.01) {
							if (world.rank() == 0  and (param.print_level()>3)) print("Increasing subspace singular value threshold ", c[m - 1], rcond);
							rcond *= 100;
						} else {
							//print("Forcing full step due to subspace malfunction");
							// c = 0.0;
							// c[m - 1] = 1.0;
							// break;
							if (world.rank() == 0 and (param.print_level()>3)) print("Restarting KAIN due to subspace malfunction");
							Q = tensorT();
							subspace.clear();
							goto restart; // fortran hat on ...
						}
				}
				//}
				END_TIMER(world, "Update subspace stuff");

				world.gop.broadcast_serializable(c, 0); // make sure everyone has same data
				if (world.rank() == 0 and (param.print_level()>3)) {
					print("Subspace solution", c);
				}
				START_TIMER(world);
				vecfuncT amo_new = zero_functions_compressed<double, 3>(world, amo.size(), false);
				vecfuncT bmo_new = zero_functions_compressed<double, 3>(world, bmo.size(), false);
				world.gop.fence();
				for (unsigned int m = 0; m < subspace.size(); ++m) {
					const vecfuncT & vm = subspace[m].first;
					const vecfuncT & rm = subspace[m].second;
					const vecfuncT vma(vm.begin(), vm.begin() + amo.size());
					const vecfuncT rma(rm.begin(), rm.begin() + amo.size());
					const vecfuncT vmb(vm.end() - bmo.size(), vm.end());
					const vecfuncT rmb(rm.end() - bmo.size(), rm.end());
					gaxpy(world, 1.0, amo_new, c(m), vma, false);
					gaxpy(world, 1.0, amo_new, -c(m), rma, false);
					gaxpy(world, 1.0, bmo_new, c(m), vmb, false);
					gaxpy(world, 1.0, bmo_new, -c(m), rmb, false);
				}
				world.gop.fence();
				END_TIMER(world, "Subspace transform");
				if (param.maxsub() <= 1) {
					subspace.clear();
				} else if (subspace.size() == size_t(param.maxsub())) {
					subspace.erase(subspace.begin());
					Q = Q(Slice(1, -1), Slice(1, -1));
				}

				do_step_restriction(world, amo, amo_new, "alpha");
				orthonormalize(world, amo_new, param.nalpha());
				amo = amo_new;

				if (!param.spin_restricted() && param.nbeta() != 0) {
					do_step_restriction(world, bmo, bmo_new, "beta");
					orthonormalize(world, bmo_new, param.nbeta());
					bmo = bmo_new;
				} else {
					bmo = amo;
				}
}

/// perform step restriction following the KAIN solver

/// Limit maximum step size to make convergence more robust
/// @param[in]          world   the world
/// @param[in]          mo              vector of orbitals from previous iteration
/// @param[inout]       new_mo  vector of orbitals from the KAIN solver
/// @param[in]          spin    "alpha" or "beta" for user information
/// @return                     max residual
double SCF::do_step_restriction(World& world, const vecfuncT& mo, vecfuncT& mo_new,
		std::string spin) const {
	PROFILE_MEMBER_FUNC(SCF);
	std::vector<double> anorm = norm2s(world, sub(world, mo, mo_new));
	int nres = 0;
	for (unsigned int i = 0; i < mo.size(); ++i) {
		if (anorm[i] > param.maxrotn()) {
			double s = param.maxrotn() / anorm[i];
			++nres;
			if (world.rank() == 0) {
				if (nres == 1 and (param.print_level()>1))
					printf("  restricting step for %s orbitals:", spin.c_str());
				printf(" %d", i);
			}
			mo_new[i].gaxpy(s, mo[i], 1.0 - s, false);
		}
	}
	if (nres > 0 && world.rank() == 0  and (param.print_level()>1))
		printf("\n");

	world.gop.fence();
	double rms, maxval;
	vector_stats(anorm, rms, maxval);
	if (world.rank() == 0  and (param.print_level()>1))
		print("Norm of vector changes", spin, ": rms", rms, "   max", maxval);
	return maxval;
}

/// orthonormalize the vectors (symmetric in occupied spaced, gramm-schmidt for virt to occ)

/// @param[in]          world   the world
/// @param[inout]       amo_new the vectors to be orthonormalized
void SCF::orthonormalize(World& world, vecfuncT& amo_new, int nocc) const {
	PROFILE_MEMBER_FUNC(SCF);
	START_TIMER(world);
	double trantol = vtol / std::min(30.0, double(amo_new.size()));
	normalize(world, amo_new);
	double maxq;
	do {
		tensorT Q = Q2(matrix_inner(world, amo_new, amo_new)); // Q3(matrix_inner(world, amo_new, amo_new))
		maxq = 0.0;
		for (int j=1; j<Q.dim(0); j++)
			for (int i=0; i<j; i++)
				maxq = std::max(std::abs(Q(j,i)),maxq);

		Q.screen(trantol); // Is this really needed? Just for speed.

		//make virt orthog to occ without changing occ states --- ASSUMES symmetric form for Q2
		for (int j=nocc; j<Q.dim(0); ++j) {
			for (int i=0; i<nocc; ++i) {
				Q(j,i)=0.0;
				Q(i,j)*=2.0;
			}
		}

		amo_new = transform(world, amo_new,
				Q, trantol, true);
		truncate(world, amo_new);
		if (world.rank() == 0 and (param.print_level()>3)) print("ORTHOG2a: maxq trantol", maxq, trantol);
		//print(Q);

	} while (maxq>0.01);
	normalize(world, amo_new);

	END_TIMER(world, "Orthonormalize");

}

/// orthonormalize the vectors ignoring occupied/virtual distinctions

/// @param[in]          world   the world
/// @param[inout]       amo_new the vectors to be orthonormalized
void SCF::orthonormalize(World& world, vecfuncT& amo_new) const {
	PROFILE_MEMBER_FUNC(SCF);
	START_TIMER(world);
	double trantol = vtol / std::min(30.0, double(amo.size()));
	normalize(world, amo_new);
	double maxq;
	do {
		tensorT Q = Q2(matrix_inner(world, amo_new, amo_new)); // Q3(matrix_inner(world, amo_new, amo_new))
		maxq = 0.0;
		for (int j=1; j<Q.dim(0); j++)
			for (int i=0; i<j; i++)
				maxq = std::max(std::abs(Q(j,i)),maxq);

		//Q.screen(trantol); // ???? Is this really needed?
		amo_new = transform(world, amo_new,
				Q, trantol, true);
		truncate(world, amo_new);
		if (world.rank() == 0  and (param.print_level()>3)) print("ORTHOG2b: maxq trantol", maxq, trantol);
		//print(Q);

	} while (maxq>0.01);
	normalize(world, amo_new);
	END_TIMER(world, "Orthonormalize");
}


void SCF::propagate(World& world, double omega, int step0) {
	PROFILE_MEMBER_FUNC(SCF);
	// Load molecular orbitals
	set_protocol < 3 > (world, 1e-4);
	make_nuclear_potential(world);
	initial_load_bal(world);
	load_mos(world);

	int nstep = 1000;
	double time_step = 0.05;

	double strength = 0.1;

	// temporary way of doing this for now
	//      VextCosFunctor<double> Vext(world,new DipoleFunctor(2),omega);
	functionT fdipx =
			factoryT(world).functor(functorT(new DipoleFunctor(0))).initial_level(
					4);
	functionT fdipy =
			factoryT(world).functor(functorT(new DipoleFunctor(1))).initial_level(
					4);
	functionT fdipz =
			factoryT(world).functor(functorT(new DipoleFunctor(2))).initial_level(
					4);

	world.gop.broadcast(time_step);
	world.gop.broadcast(nstep);

	// Need complex orbitals :(
			double thresh = 1e-4;
	cvecfuncT camo = zero_functions<double_complex, 3>(world, param.nalpha());
	cvecfuncT cbmo = zero_functions<double_complex, 3>(world, param.nbeta());
	for (int iorb = 0; iorb < param.nalpha(); iorb++) {
		camo[iorb] = std::exp(double_complex(0.0, 2 * constants::pi * strength))
		* amo[iorb];
		camo[iorb].truncate(thresh);
	}
	if (!param.spin_restricted() && param.nbeta()) {
		for (int iorb = 0; iorb < param.nbeta(); iorb++) {
			cbmo[iorb] = std::exp(
					double_complex(0.0, 2 * constants::pi * strength))
			* bmo[iorb];
			cbmo[iorb].truncate(thresh);
		}
	}

	// Create free particle propagator
	// Have no idea what to set "c" to
	double c = 20.0;
	printf("Creating G\n");
	Convolution1D < double_complex > *G = qm_1d_free_particle_propagator(
			FunctionDefaults < 3 > ::get_k(), c, 0.5 * time_step,
			2.0 * param.L());
	printf("Done creating G\n");

	// Start iteration over time
	for (int step = 0; step < nstep; step++) {
		//        if (world.rank() == 0) printf("Iterating step %d:\n\n", step);
		double t = time_step * step;
		//        iterate_trotter(world, G, Vext, camo, cbmo, t, time_step);
		iterate_trotter(world, G, camo, cbmo, t, time_step, thresh);
		functionT arho = make_density(world, aocc, camo);
		functionT brho =
				(!param.spin_restricted() && param.nbeta()) ?
						make_density(world, aocc, camo) : copy(arho);
		functionT rho = arho + brho;
		double xval = inner(fdipx, rho);
		double yval = inner(fdipy, rho);
		double zval = inner(fdipz, rho);
		if (world.rank() == 0)
			printf("%15.7f%15.7f%15.7f%15.7f\n", t, xval, yval, zval);
	}
}

complex_functionT APPLY(const complex_operatorT* q1d,
		const complex_functionT& psi) {
	complex_functionT r = psi; // Shallow copy violates constness !!!!!!!!!!!!!!!!!
	coordT lo, hi;
	lo[2] = -10;
	hi[2] = +10;

	r.reconstruct();
	r.broaden();
	r.broaden();
	r.broaden();
	r.broaden();
	r = apply_1d_realspace_push(*q1d, r, 2);
	r.sum_down();
	r = apply_1d_realspace_push(*q1d, r, 1);
	r.sum_down();
	r = apply_1d_realspace_push(*q1d, r, 0);
	r.sum_down();

	return r;
}

void SCF::iterate_trotter(World& world, Convolution1D<double_complex>* G,
		cvecfuncT& camo, cvecfuncT& cbmo, double t, double time_step,
		double thresh) {
	PROFILE_MEMBER_FUNC(SCF);

	// first kinetic energy apply
	cvecfuncT camo2 = zero_functions<double_complex, 3>(world, param.nalpha());
	cvecfuncT cbmo2 = zero_functions<double_complex, 3>(world, param.nbeta());
	for (int iorb = 0; iorb < param.nalpha(); iorb++) {
		//        if (world.rank()) printf("Apply free-particle Green's function to alpha orbital %d\n", iorb);
		camo2[iorb] = APPLY(G, camo[iorb]);
		camo2[iorb].truncate(thresh);
	}
	if (!param.spin_restricted() && param.nbeta()) {
		for (int iorb = 0; iorb < param.nbeta(); iorb++) {
			cbmo2[iorb] = APPLY(G, cbmo[iorb]);
			cbmo2[iorb].truncate(thresh);
		}
	}
	// Construct new density
	//      START_TIMER(world);
	functionT arho = make_density(world, aocc, amo), brho;

	if (param.nbeta()) {
		if (param.spin_restricted()) {
			brho = arho;
		} else {
			brho = make_density(world, bocc, bmo);
		}
	} else {
		brho = functionT(world); // zero
	}
	functionT rho = arho + brho;
	//      END_TIMER(world, "Make densities");

	// Do RPA only for now
	real_function_3d vnuc = potentialmanager->vnuclear();
	functionT vlocal = vnuc;
	//      START_TIMER(world);
	functionT vcoul = apply(*coulop, rho);
	//      END_TIMER(world, "Coulomb");
	//      vlocal += vcoul + Vext(t+0.5*time_step);
	//      vlocal += vcoul + std::cos(0.1*(t+0.5*time_step))*fdip;

	// exponentiate potential
	//      if (world.rank()) printf("Apply Kohn-Sham potential to orbitals\n");
	complex_functionT expV = make_exp(time_step, vlocal);
	cvecfuncT camo3 = mul_sparse(world, expV, camo2, vtol, false);
	world.gop.fence();

	// second kinetic energy apply
	for (int iorb = 0; iorb < param.nalpha(); iorb++) {
		//        if (world.rank() == 0) printf("Apply free-particle Green's function to alpha orbital %d\n", iorb);
		camo3[iorb].truncate(thresh);
		camo[iorb] = APPLY(G, camo3[iorb]);
		camo[iorb].truncate();
	}
	if (!param.spin_restricted() && param.nbeta()) {
     	        cvecfuncT cbmo3 = mul_sparse(world, expV, cbmo2, vtol); // Removed nofence --- must fence here

		// second kinetic energy apply
		for (int iorb = 0; iorb < param.nbeta(); iorb++) {
			cbmo[iorb] = APPLY(G, cbmo3[iorb]);
			cbmo[iorb].truncate();
		}
	}
}

// For given protocol, solve the DFT/HF/response equations
void SCF::solve(World & world) {
	PROFILE_MEMBER_FUNC(SCF);
	functionT arho_old, brho_old;
	const double dconv = std::max(FunctionDefaults < 3 > ::get_thresh(),
			param.dconv());
	const double trantol = vtol / std::min(30.0, double(amo.size()));
	const double tolloc = 1e-6; // was std::min(1e-6,0.01*dconv) but now trying to avoid unnecessary change
	double update_residual = 0.0, bsh_residual = 0.0;
	subspaceT subspace;
	tensorT Q;
	bool do_this_iter = true;
	bool converged = false;

	// Shrink subspace until stop localizing/canonicalizing--- probably not a good idea
	// int maxsub_save = param.maxsub;
	// param.maxsub = 2;

	for (int iter = 0; iter < param.maxiter(); ++iter) {
		if (world.rank() == 0 and (param.print_level()>1))
			printf("\nIteration %d at time %.1fs\n\n", iter, wall_time());

		// if (iter > 0 && update_residual < 0.1) {
		//     //do_this_iter = false;
		//     param.maxsub = maxsub_save;
		// }

		if (param.do_localize() && do_this_iter) {
			distmatT dUT;
			if (param.localize_pm()) {
				dUT = localize_PM(world, amo, aset, tolloc, 0.1, iter == 0, false);
			}
			else if (param.localize_method()=="new") {
				dUT = localize_new(world, amo, aset, tolloc, 0.1, iter == 0, false);
			}
			else if (param.localize_method()=="boys") {
				dUT = localize_boys(world, amo, aset, tolloc, 0.1, iter == 0, false);
			}
			else
				throw "localization confusion";

			dUT.data().screen(trantol);
			START_TIMER(world);
			amo = transform(world, amo, dUT);
			truncate(world, amo);
			normalize(world, amo);
			if (!param.spin_restricted() && param.nbeta() != 0) {
				if (param.localize_pm()) {
					dUT = localize_PM(world, bmo, bset, tolloc, 0.1, iter == 0, false);
				}
				else if (param.localize_method()=="new") {
					dUT = localize_new(world, bmo, bset, tolloc, 0.1, iter == 0, false);
				}
				else {
					dUT = localize_boys(world, bmo, bset, tolloc, 0.1, iter == 0, false);
				}

				START_TIMER(world);
				dUT.data().screen(trantol);
				bmo = transform(world, bmo, dUT);
				truncate(world, bmo);
				normalize(world, bmo);
				END_TIMER(world, "Rotate subspace");
			}
		}

		START_TIMER(world);
		functionT arho = make_density(world, aocc, amo), brho;

		if (param.nbeta()) {
			if (param.spin_restricted()) {
				brho = arho;
			} else {
				brho = make_density(world, bocc, bmo);
			}
		} else {
			brho = functionT(world); // zero
		}
		END_TIMER(world, "Make densities");
		print_meminfo(world.rank(), "Make densities");

		if (iter < 2 || (iter % 10) == 0) {
			START_TIMER(world);
			loadbal(world, arho, brho, arho_old, brho_old, subspace);
			END_TIMER(world, "Load balancing");
			print_meminfo(world.rank(), "Load balancing");
		}
		double da = 0.0, db = 0.0;
		if (iter > 0) {
			da = (arho - arho_old).norm2();
			db = (brho - brho_old).norm2();
			if (world.rank() == 0 and (param.print_level()>2))
				print("delta rho", da, db, "residuals", bsh_residual,
						update_residual);

		}

		START_TIMER(world);
		arho_old = arho;
		brho_old = brho;
		functionT rho = arho + brho;
		rho.truncate();

		real_function_3d vnuc;
		if (param.psp_calc()){
			vnuc = gthpseudopotential->vlocalpot();}
		else if (param.pure_ae()){
			vnuc = potentialmanager->vnuclear();}
		else {
			vnuc = potentialmanager->vnuclear();
			vnuc = vnuc + gthpseudopotential->vlocalpot();}
		double enuclear = inner(rho, vnuc);
		END_TIMER(world, "Nuclear energy");

		START_TIMER(world);
		functionT vcoul = apply(*coulop, rho);
		functionT vlocal;
		END_TIMER(world, "Coulomb");
		print_meminfo(world.rank(), "Coulomb");

		double ecoulomb = 0.5 * inner(rho, vcoul);
		rho.clear(false);
		vlocal = vcoul + vnuc;

		// compute the contribution of the solvent to the local potential
		double epcm=0.0;
		if (param.pcm_data() != "none") {
			START_TIMER(world);
			functionT vpcm=pcm.compute_pcm_potential(vcoul);
			vlocal+=vpcm;
			epcm=pcm.compute_pcm_energy();
			END_TIMER(world, "PCM");
			print_meminfo(world.rank(), "PCM");
		}

		vcoul.clear(false);
		vlocal.truncate();
		double exca = 0.0, excb = 0.0;

		double enla = 0.0, enlb = 0.0;
		vecfuncT Vpsia = apply_potential(world, aocc, amo, vlocal, exca, enla, 0);
		vecfuncT Vpsib;
		if (!param.spin_restricted() && param.nbeta()) {
			Vpsib = apply_potential(world, bocc, bmo, vlocal, excb, enlb, 1);
		}
		else if (param.nbeta() != 0) {
			enlb = enla;
		}

		double ekina = 0.0, ekinb = 0.0;
		tensorT focka = make_fock_matrix(world, amo, Vpsia, aocc, ekina);
		tensorT fockb = focka;

		if (!param.spin_restricted() && param.nbeta() != 0)
			fockb = make_fock_matrix(world, bmo, Vpsib, bocc, ekinb);
		else if (param.nbeta() != 0) {
			ekinb = ekina;
		}

		if (!param.do_localize() && do_this_iter) {
			tensorT U = diag_fock_matrix(world, focka, amo, Vpsia, aeps, aocc,
					FunctionDefaults < 3 > ::get_thresh());
			//rotate_subspace(world, U, subspace, 0, amo.size(), trantol); ??
					if (!param.spin_restricted() && param.nbeta() != 0) {
						U = diag_fock_matrix(world, fockb, bmo, Vpsib, beps, bocc,
								FunctionDefaults < 3 > ::get_thresh());
						//rotate_subspace(world, U, subspace, amo.size(), bmo.size(),trantol);
					}
		}

		double enrep = molecule.nuclear_repulsion_energy();
		double ekinetic = ekina + ekinb;
		double enonlocal = enla + enlb;
		double exc = exca + excb;
		double etot = ekinetic + enuclear + ecoulomb + exc + enrep + enonlocal + epcm;
		current_energy = etot;
		//esol = etot;

		if (world.rank() == 0  and (param.print_level()>1)) {
			//lots of dps for testing Exc stuff
			/*printf("\n              kinetic %32.24f\n", ekinetic);
                printf("         nonlocal psp %32.24f\n", enonlocal);
                printf("   nuclear attraction %32.24f\n", enuclear);
                printf("              coulomb %32.24f\n", ecoulomb);
                printf(" exchange-correlation %32.24f\n", exc);
                printf("    nuclear-repulsion %32.24f\n", enrep);
                printf("                total %32.24f\n\n", etot);*/

			printf("\n              kinetic %16.8f\n", ekinetic);
			printf("         nonlocal psp %16.8f\n", enonlocal);
			printf("   nuclear attraction %16.8f\n", enuclear);
			printf("              coulomb %16.8f\n", ecoulomb);
			printf("                  PCM %16.8f\n", epcm);
			printf(" exchange-correlation %16.8f\n", exc);
			printf("    nuclear-repulsion %16.8f\n", enrep);
			printf("                total %16.8f\n\n", etot);
		}

		if (iter > 0) {
			//print("##convergence criteria: density delta=", da < dconv * molecule.natom() && db < dconv * molecule.natom(), ", bsh_residual=", (param.conv_only_dens || bsh_residual < 5.0*dconv));
			if (da < dconv * std::max(size_t(5),molecule.natom()) && db < dconv * std::max(size_t(5),molecule.natom())
			&& (param.get<bool>("conv_only_dens") || bsh_residual < 5.0 * dconv)) converged=true;
			// previous conv was too tight for small systems
			// if (da < dconv * molecule.natom() && db < dconv * molecule.natom()
			//     && (param.conv_only_dens || bsh_residual < 5.0 * dconv)) converged=true;

			// do diagonalization etc if this is the last iteration, even if the calculation didn't converge
			if (converged || iter==param.maxiter()-1) {
				if (world.rank() == 0 && converged and (param.print_level()>1)) {
					print("\nConverged!\n");
				}

				// Diagonalize to get the eigenvalues and if desired the final eigenvectors
				tensorT U;
				START_TIMER(world);
				tensorT overlap = matrix_inner(world, amo, amo, true);
				END_TIMER(world, "Overlap");

				START_TIMER(world);
				sygvp(world, focka, overlap, 1, U, aeps);
				END_TIMER(world, "focka eigen sol");

				if (!param.do_localize()) {
					START_TIMER(world);
					amo = transform(world, amo, U, trantol, true);
					truncate(world, amo);
					normalize(world, amo);
					END_TIMER(world, "Transform MOs");
				}
				if (param.nbeta() != 0 && !param.spin_restricted()) {

					START_TIMER(world);
					overlap = matrix_inner(world, bmo, bmo, true);
					END_TIMER(world, "Overlap");

					START_TIMER(world);
					sygvp(world, fockb, overlap, 1, U, beps);
					END_TIMER(world, "fockb eigen sol");

					if (!param.do_localize()) {
						START_TIMER(world);
						bmo = transform(world, bmo, U, trantol, true);
						truncate(world, bmo);
						normalize(world, bmo);
						END_TIMER(world, "Transform MOs");
					}
				}

				if (world.rank() == 0 and (param.print_level()>1)) {
					print(" ");
					print("alpha eigenvalues");
					print (aeps);
					if (param.nbeta() != 0 && !param.spin_restricted()) {
						print("beta eigenvalues");
						print (beps);
					}


					// write eigenvalues etc to a file at the same time for plotting DOS etc.
					FILE *f=0;
					if (param.nbeta() != 0 && !param.spin_restricted()) {
						f = fopen("energies_alpha.dat", "w");}
					else{
						f = fopen("energies.dat", "w");}

					long nmo = amo.size();
					fprintf(f, "# %8li\n", nmo);
					for (long i = 0; i < nmo; ++i) {
						fprintf(f, "%13.8f\n", aeps(i));
					}
					fclose(f);

					if (param.nbeta() != 0 && !param.spin_restricted()) {
						long nmo = bmo.size();
						FILE *f=0;
						f = fopen("energies_beta.dat", "w");

						fprintf(f, "# %8li\n", nmo);
						for (long i = 0; i < nmo; ++i) {
							fprintf(f, "%13.8f\t", beps(i));
						}
						fclose(f);
					}

				}

				if (param.do_localize()) {
					// Restore the diagonal elements for the analysis
					for (unsigned int i = 0; i < amo.size(); ++i)
						aeps[i] = focka(i, i);
					if (param.nbeta() != 0 && !param.spin_restricted())
						for (unsigned int i = 0; i < bmo.size(); ++i)
							beps[i] = fockb(i, i);
				}

				break;
			}

		}

		update_subspace(world, Vpsia, Vpsib, focka, fockb, subspace, Q,
				bsh_residual, update_residual);

	}

	// compute the dipole moment
	functionT rho = make_density(world, aocc, amo);
	if (!param.spin_restricted()) {
		if (param.nbeta())
			rho += make_density(world, bocc, bmo);
	} else {
		rho.scale(2.0);
	}
	dipole(world,rho);

	if (world.rank() == 0 and (param.print_level()>1)) {
		if (param.do_localize())
			print(
					"Orbitals are localized - energies are diagonal Fock matrix elements\n");
		else
			print("Orbitals are eigenvectors - energies are eigenvalues\n");
		print("Analysis of alpha MO vectors");
	}

	analyze_vectors(world, amo, aocc, aeps);
	if (param.nbeta() != 0 && !param.spin_restricted()) {
		if (world.rank() == 0 and (param.print_level()>1))
			print("Analysis of beta MO vectors");

		analyze_vectors(world, bmo, bocc, beps);
	}

	if (param.get<bool>("print_dipole_matels")) {
		dipole_matrix_elements(world, amo, aocc, aeps, 0);
		if (param.nbeta() != 0 && !param.spin_restricted()) {
			dipole_matrix_elements(world, bmo, bocc, beps, 1);
		}
	}
}        // end solve function


//vama polarizability
void SCF::update_response_subspace(World & world,
		vecfuncT & ax, vecfuncT & ay,
		vecfuncT & bx, vecfuncT & by,
		vecfuncT & rax, vecfuncT & ray,
		vecfuncT & rbx, vecfuncT & rby,
		subspaceT & subspace, tensorT & Q, double & update_residual)
{
	vecfuncT vm = ax;
	vm.insert(vm.end(), ay.begin(), ay.end());

	vecfuncT rm = rax;
	rm.insert(rm.end(), ray.begin(), ray.end());

	if(param.nbeta() != 0 && !param.spin_restricted()){
		vm.insert(vm.end(), bx.begin(), bx.end());
		vm.insert(vm.end(), by.begin(), by.end());
		rm.insert(rm.end(), rbx.begin(), rbx.end());
		rm.insert(rm.end(), rby.begin(), rby.end());
	}

	compress(world, vm, false);
	compress(world, rm, false);
	world.gop.fence();
	subspace.push_back(pairvecfuncT(vm, rm));
	int m = subspace.size();
	tensorT ms(m);
	tensorT sm(m);
	for(int s = 0;s < m;++s){
		const vecfuncT & vs = subspace[s].first;
		const vecfuncT & rs = subspace[s].second;
		for(unsigned int i = 0;i < vm.size();++i){
			ms[s] += vm[i].inner_local(rs[i]);
			sm[s] += vs[i].inner_local(rm[i]);
		}
	}

	world.gop.sum(ms.ptr(), m);
	world.gop.sum(sm.ptr(), m);
	tensorT newQ(m, m);
	if(m > 1)
		newQ(Slice(0, -2), Slice(0, -2)) = Q;

	newQ(m - 1, _) = ms;
	newQ(_, m - 1) = sm;
	Q = newQ;
	//if (world.rank() == 0) { print("kain Q"); print(Q); }
	tensorT c;
	if(world.rank() == 0){
		double rcond = 1e-12;
		while(1){
			c = KAIN(Q, rcond);
			//if (world.rank() == 0) print("kain c:", c);
			if(std::abs(c[m - 1]) < 3.0){
				break;
			} else if(rcond < 0.01){
				if (param.print_level()>3) print("Increasing subspace singular value threshold ", c[m - 1], rcond);
				rcond *= 100;
			} else {
				if(param.print_level()>3) print("Forcing full step due to subspace malfunction");
				c = 0.0;
				c[m - 1] = 1.0;
				break;
			}
		}
	}

	world.gop.broadcast_serializable(c, 0);
	if(world.rank() == 0  and (param.print_level()>3)){
		print("Response Subspace solution", c);
	}
	START_TIMER(world);
	vecfuncT ax_new = zero_functions_compressed<double,3>(world, ax.size());
	vecfuncT ay_new = zero_functions_compressed<double,3>(world, ay.size());
	vecfuncT bx_new = zero_functions_compressed<double,3>(world, bx.size());
	vecfuncT by_new = zero_functions_compressed<double,3>(world, by.size());

	for(unsigned int m = 0;m < subspace.size();++m){
		const vecfuncT & vm = subspace[m].first;
		const vecfuncT & rm = subspace[m].second;
		const vecfuncT vmax(vm.begin(), vm.begin() + ax.size());
		const vecfuncT rmax(rm.begin(), rm.begin() + rax.size());
		const vecfuncT vmay(vm.begin() + ax.size(), vm.begin() + ax.size() + ay.size());
		const vecfuncT rmay(rm.begin() + rax.size(), rm.begin() + rax.size() + ray.size());
		gaxpy(world, 1.0, ax_new, c(m), vmax, false);
		gaxpy(world, 1.0, ax_new, -c(m), rmax, false);
		gaxpy(world, 1.0, ay_new, c(m), vmay, false);
		gaxpy(world, 1.0, ay_new, -c(m), rmay, false);
		//if(param.nbeta != 0 && !param.spin_restricted){
			const vecfuncT vmbx(vm.end() - by.size() - bx.size(), vm.end() - by.size());
			const vecfuncT rmbx(rm.end() - rby.size() - rbx.size(), rm.end() - rby.size());
			const vecfuncT vmby(vm.end() - by.size(), vm.end());
			const vecfuncT rmby(rm.end() - rby.size(), rm.end());
			gaxpy(world, 1.0, bx_new, c(m), vmbx, false);
			gaxpy(world, 1.0, bx_new, -c(m), rmbx, false);
			gaxpy(world, 1.0, by_new, c(m), vmby, false);
			gaxpy(world, 1.0, by_new, -c(m), rmby, false);
			//}
}
	world.gop.fence();
	END_TIMER(world, "Subspace transform");
	if(param.maxsub() <= 1){
		subspace.clear();
	} else if(subspace.size() == size_t(param.maxsub())){
		subspace.erase(subspace.begin());
		Q = Q(Slice(1, -1), Slice(1, -1));
	}

	std::vector<double> axnorm = norm2s(world, sub(world, ax, ax_new));
	std::vector<double> aynorm = norm2s(world, sub(world, ay, ay_new));
	std::vector<double> bxnorm = norm2s(world, sub(world, bx, bx_new));
	std::vector<double> bynorm = norm2s(world, sub(world, by, by_new));
	int nres = 0;
	for(unsigned int i = 0;i < ax.size();++i){
		if(axnorm[i] > param.maxrotn()){
			double s = param.maxrotn() / axnorm[i];
			++nres;
			if(world.rank() == 0  and (param.print_level()>3)){
				if(nres == 1)
					printf("  restricting step for alpha orbitals:");

				printf(" %d", i);
			}
			ax_new[i].gaxpy(s, ax[i], 1.0 - s, false);
		}

	}
	if(nres > 0 && world.rank() == 0  and (param.print_level()>3))
		printf("\n");

	nres = 0;
	for(unsigned int i = 0;i < ay.size();++i){
		if(aynorm[i] > param.maxrotn()){
			double s = param.maxrotn() / aynorm[i];
			++nres;
			if(world.rank() == 0  and (param.print_level()>3)){
				if(nres == 1)
					printf("  restricting step for alpha orbitals:");

				printf(" %d", i);
			}
			ay_new[i].gaxpy(s, ay[i], 1.0 - s, false);
		}

	}
	if(nres > 0 && world.rank() == 0  and (param.print_level()>3))
		printf("\n");

	//if(param.nbeta != 0 && !param.spin_restricted){
	nres = 0;
	for(unsigned int i = 0;i < bx.size();++i){
		if(bxnorm[i] > param.maxrotn()){
			double s = param.maxrotn() / bxnorm[i];
			++nres;
			if(world.rank() == 0  and (param.print_level()>3)){
				if(nres == 1)
					printf("  restricting step for  beta orbitals:");

				printf(" %d", i);
			}
			bx_new[i].gaxpy(s, bx[i], 1.0 - s, false);
		}

	}
	if(nres > 0 && world.rank() == 0  and (param.print_level()>3))
		printf("\n");

	nres = 0;
	for(unsigned int i = 0;i < by.size();++i){
		if(bynorm[i] > param.maxrotn()){
			double s = param.maxrotn() / bynorm[i];
			++nres;
			if(world.rank() == 0  and (param.print_level()>3)){
				if(nres == 1)
					printf("  restricting step for  beta orbitals:");

				printf(" %d", i);
			}
			by_new[i].gaxpy(s, by[i], 1.0 - s, false);
		}

	}
	if(nres > 0 && world.rank() == 0  and (param.print_level()>3))
		printf("\n");
	//}

	world.gop.fence();
	double rms, maxval_x, maxval_y, maxval_b;
	vector_stats(axnorm, rms, maxval_x);
	vector_stats(aynorm, rms, maxval_y);

	update_residual = std::max(maxval_x, maxval_y);

	if(bxnorm.size()){
		vector_stats(bxnorm, rms, maxval_x);
		vector_stats(bynorm, rms, maxval_y);

		maxval_b = std::max(maxval_x, maxval_y);
		update_residual = std::max(update_residual, maxval_b);
	}
	//START_TIMER(world);
	//double trantol = vtol / std::min(30.0, double(ax.size()));
	//normalize(world, ax_new);
	//normalize(world, ay_new);
	//ax_new = transform(world, ax_new, Q3(matrix_inner(world, ax_new, ax_new)), trantol, true);
	//ay_new = transform(world, ay_new, Q3(matrix_inner(world, ay_new, ay_new)), trantol, true);
	truncate(world, ax_new);
	truncate(world, ay_new);
	//normalize(world, ax_new);
	//normalize(world, ay_new);
	if(param.nbeta() != 0  && !param.spin_restricted()){
		//normalize(world, bx_new);
		//normalize(world, by_new);
		//bx_new = transform(world, bx_new, Q3(matrix_inner(world, bx_new, bx_new)), trantol, true);
		//by_new = transform(world, by_new, Q3(matrix_inner(world, by_new, by_new)), trantol, true);
		truncate(world, bx_new);
		truncate(world, by_new);
		//normalize(world, bx_new);
		//normalize(world, by_new);
	}
	//END_TIMER(world, "Orthonormalize");
	ax = ax_new;
	ay = ay_new;
	bx = bx_new;
	by = by_new;
}

vecfuncT SCF::apply_potential_response(World & world, const vecfuncT & dmo,
		const XCOperator& xcop,  const functionT & vlocal, int ispin)
{
	functionT vloc = copy(vlocal);

	if (xc.is_dft() && !(xc.hf_exchange_coefficient()==1.0)) {
		START_TIMER(world);

		//            XCOperator xcoperator(world,this,ispin);
		//            if (ispin==0) exc=xcoperator.compute_xc_energy();
		xcop.set_ispin(ispin);
		vloc += xcop.make_xc_potential();

		// TODO: fbischoff thinks this is double-counting the gga potential part
		//
		//#ifdef MADNESS_HAS_LIBXC
		//            if (xc.is_gga() ) {
		//
		//                functionT vsigaa = xcoperator.make_xc_potential();
		//                functionT vsigab;
		//                if (xc.is_spin_polarized() && param.nbeta != 0)// V_ab
		//                    vsigab = xcoperator.make_xc_potential();
		//
		//                for (int axis=0; axis<3; axis++) {
		//                    functionT gradn = delrho[axis + 3*ispin];
		//                    functionT ddel = vsigaa*gradn;
		//                    if (xc.is_spin_polarized() && param.nbeta != 0) {
		//                        functionT vsab = vsigab*delrho[axis + 3*(1-ispin)];
		//                        ddel = ddel + vsab;
		//                    }
		//                    ddel.scale(xc.is_spin_polarized() ? 2.0 : 4.0);
		//                    Derivative<double,3> D = free_space_derivative<double,3>(world, axis);
		//                    functionT vxc2=D(ddel);
		//                    vloc = vloc - vxc2;//.truncate();
		//                }
		//            }
		//#endif
		END_TIMER(world, "DFT potential");
	}



	START_TIMER(world);
	vecfuncT Vdmo = mul_sparse(world, vloc, dmo, vtol);
	END_TIMER(world, "V*dmo");
	print_meminfo(world.rank(), "V*dmo");
	if(xc.hf_exchange_coefficient()){
		START_TIMER(world);
		vecfuncT Kdmo;
		Exchange<double,3> K=Exchange<double,3>(world,this,ispin).small_memory(false).same(false);
		if(ispin == 0)
			Kdmo=K(amo);
		if(ispin == 1)
			Kdmo=K(bmo);
		//tensorT excv = inner(world, Kdmo, dmo);
		//double exchf = 0.0;
		//for(unsigned long i = 0;i < dmo.size();++i){
		//    exchf -= 0.5 * excv[i] * occ[i];
		//}
		//if (!xc.is_spin_polarized()) exchf *= 2.0;
		gaxpy(world, 1.0, Vdmo, -xc.hf_exchange_coefficient(), Kdmo);
		Kdmo.clear();
		END_TIMER(world, "HF exchange");
		//exc = exchf* xc.hf_exchange_coefficient() + exc;
	}
	if (param.pure_ae())
		potentialmanager->apply_nonlocal_potential(world, amo, Vdmo);

	truncate(world, Vdmo);

	print_meminfo(world.rank(), "Truncate Vdmo");
	world.gop.fence();
	return Vdmo;
}

void SCF::this_axis(World & world, int axis)
{
	print("\n");
	if (world.rank() == 0) {
		if(axis == 0)
			print(" AXIS of frequency = x");

		else if(axis == 1)
			print(" AXIS of frequency = y");

		else if(axis == 2)
			print(" AXIS of frequency = z");
	}
}

vecfuncT SCF::calc_dipole_mo(World & world,  vecfuncT & mo, const int axis)
{
	//START_TIMER(world);

	vecfuncT dipolemo = zero_functions<double,3>(world, mo.size());

	std::vector<int> f(3, 0);
	f[axis] = true;
	//print("f = ", f[0]," ",  f[1], " ", f[2]);
	functionT dipolefunc = factoryT(world).functor(functorT(new MomentFunctor(f)));
	reconstruct(world, mo);

	// dipolefunc * mo[iter]
	for(size_t p=0; p<mo.size(); ++p)
		dipolemo[p] =  mul_sparse(dipolefunc, mo[p],false);
	world.gop.fence(); // Must fence here

	//END_TIMER(world, "Make perturbation");
	print_meminfo(world.rank(), "Make perturbation");

	truncate(world, dipolemo);
	return dipolemo;
}

void SCF::calc_freq(World & world, double & omega, tensorT & ak, tensorT & bk, int sign)
{

	for(int i=0; i<param.nalpha(); ++i){
		ak[i] = sqrt(-2.0 * (aeps[i] + sign * omega));
		if (world.rank() == 0)
			print(" kxy(alpha) [", i, "] : sqrt(-2 * (eps +/- omega)) = ", ak[i]);
	}
	if(!param.spin_restricted() && param.nbeta() != 0) {
		for(int i=0; i<param.nbeta(); ++i){
			bk[i] = sqrt(-2.0 * (beps[i] + sign * omega));
			if (world.rank() == 0)
				if (world.rank() == 0)
					print(" kxy(beta) [", i, "]: sqrt(-2 * (eps +/- omega)) = ", bk[i]);
		}
	}
}

void SCF::make_BSHOperatorPtr(World & world, tensorT & ak, tensorT & bk,
		std::vector<poperatorT> & aop, std::vector<poperatorT> & bop)
{
	//START_TIMER(world);
	double tol = FunctionDefaults < 3 > ::get_thresh();

	for(int i=0; i<param.nalpha(); ++i) {
		// thresh tol : 1e-6
		aop[i] = poperatorT(BSHOperatorPtr3D(world, ak[i], param.lo() , tol));
	}
	if(!param.spin_restricted() && param.nbeta() != 0) {
		for(int i=0; i<param.nbeta(); ++i) {
			bop[i] = poperatorT(BSHOperatorPtr3D(world, bk[i], param.lo(), tol));
		}
	}
	//END_TIMER(world, "Make BSHOp");
	print_meminfo(world.rank(), "Make BSHOp");
}

functionT SCF::make_density_ground(World & world, functionT & arho, functionT & brho)
{

	//START_TIMER(world);

	functionT rho = factoryT(world);

	arho = make_density(world, aocc, amo);
	if (!param.spin_restricted()) {
		brho = make_density(world, bocc, bmo);
	}
	else {
		brho = arho;
	}

	rho = arho + brho;
	rho.truncate();

	//END_TIMER(world, "Make densities");
	print_meminfo(world.rank(), "Make densities");

	return rho;
}

functionT SCF::make_derivative_density(World & world, const vecfuncT & mo,
		const tensorT & occ ,
		const vecfuncT & x, const vecfuncT & y)
{
	functionT drho = factoryT(world);
	drho.compress();
	for(size_t i=0; i<mo.size(); ++i) {
		functionT rhoi = mo[i] * x[i] + mo[i] * y[i];
		rhoi.compress();
		if(occ[i])
			drho.gaxpy(1.0, rhoi, occ[i], false);
		// drho += (mo[i] * x[i]) + (mo[i] * y[i]);
	}
	world.gop.fence();
	//drho.truncate();
	return drho;
}


functionT SCF::calc_exchange_function(World & world,  const int & p,
		const vecfuncT & dmo1,  const vecfuncT & dmo2,
		const vecfuncT & mo, int & spin)
{

	functionT dKmo = factoryT(world);
	reconstruct(world, mo);
	reconstruct(world, dmo1);
	reconstruct(world, dmo2);

	functionT k1 = factoryT(world);
	functionT k2 = factoryT(world);
	for(size_t i=0; i<mo.size(); ++i) {
		k1 = apply(*coulop, ( mo[i] * mo[p] )) * dmo1[i];
		k2 = apply(*coulop, ( mo[p] * dmo2[i] )) * mo[i];
		dKmo = dKmo - (k1 + k2);
		k1.clear(false);
		k2.clear(false);
	}
	dKmo.truncate();
	return dKmo;
}


/// param[in]   drho    the perturbed density
vecfuncT SCF::calc_xc_function(World & world, XCOperator& xc_alda,
		const vecfuncT & mo,  const functionT & drho)
{
	START_TIMER(world);
	reconstruct(world, mo);

	functionT dJ = apply(*coulop, drho);
	dJ.truncate();

	functionT vloc = dJ;

	// TODO openshell ?
			if (xc.is_dft() && xc.hf_exchange_coefficient() !=1.0) {
				vloc =  dJ +  xc_alda.apply_xc_kernel(drho);
			}

			vecfuncT Vxcmo = mul_sparse(world, vloc, mo, vtol);
			truncate(world, Vxcmo);

			END_TIMER(world, "Calc calc_xc_function ");
			print_meminfo(world.rank(), "Calc calc_xc_function");
			return Vxcmo;
}

/// @param[in]  drho    the perturbed density
vecfuncT SCF::calc_djkmo(World & world, XCOperator& xc_alda, const vecfuncT & dmo1,
		const vecfuncT & dmo2,  const functionT & drho, const vecfuncT & mo,
		const functionT & drhos,
		int  spin)
{

	vecfuncT djkmo = zero_functions<double,3>(world, mo.size());
	// TODO becareful with drhos , drho
	// open shell drhoa !=drhob

	vecfuncT dkxcmo = calc_xc_function(world, xc_alda, mo, drho);
	//TODO hybrdid functs: not sure if should i have to apply
	if(xc.hf_exchange_coefficient() == 1.0){
		//if(xc.hf_exchange_coefficient()){
		START_TIMER(world);
		for(size_t p=0; p<mo.size(); ++p) {
			djkmo[p] = calc_exchange_function(world, p, dmo1, dmo2, mo,spin);
			//add a fraction only
			djkmo[p].scale(xc.hf_exchange_coefficient());
		}
		END_TIMER(world, "Calc calc_exchange_function ");
		print_meminfo(world.rank(), "Calc calc_exchange_function");
	}
	gaxpy(world, 1.0, djkmo, 1.0, dkxcmo);
	truncate(world, djkmo);

	return djkmo;
}

vecfuncT SCF::calc_rhs(World & world, const vecfuncT & mo ,
		const vecfuncT & Vdmo,
		const vecfuncT & dipolemo, const vecfuncT & djkmo )
{
	//START_TIMER(world);
	vecfuncT rhs;

	// the projector on the unperturbed density
	Projector<double,3> rho0(mo);

	vecfuncT gp = add(world, dipolemo, djkmo);
	for (size_t i=0; i<Vdmo.size(); ++i) {
		functionT gp1 =  gp[i];
		gp1 = gp1 - rho0(gp1);
		gp1 = Vdmo[i] + gp1 ;
		rhs.push_back(gp1);
	}

	//END_TIMER(world, "Sum rhs response");
	print_meminfo(world.rank(), "Sum rhs response");
	truncate(world, rhs);

	return rhs;
}


void SCF::calc_response_function(World & world, vecfuncT & dmo,
		std::vector<poperatorT> & op, vecfuncT & rhs)
{
	// new response function
	// BSHOperatorPrt3D : op
	dmo = apply(world, op, rhs);
	scale(world, dmo, -2.0);
	truncate(world, dmo);
}

// orthogonalization
void SCF::orthogonalize_response(World & world, vecfuncT & dmo, vecfuncT & mo )
{
	reconstruct(world, dmo);
	for(size_t i=0; i<mo.size(); ++i){
		for (size_t j=0; j<mo.size(); ++j){
			// new_x = new_x - < psi | new_x > * psi
			dmo[i] = dmo[i] - dmo[i].inner(mo[j])*mo[j];
		}
	}
}


//vama ugly ! alpha_ij(w) = - sum(m occ) [<psi_m(0)|r_i|psi_mj(1)(w)> + <psi_mj(1)(-w)|r_i|psi_m(0)>]

void SCF::dpolar(World & world, tensorT & polar, functionT & drho, const int axis)
{
	for(int i=0; i<3; ++i) {
		std::vector<int> f(3, 0);
		f[i] = true;
		functionT dipolefunc = factoryT(world).functor(functorT(new MomentFunctor(f)));
		polar(axis, i) = -2.0 * dipolefunc.inner(drho);
	}
}

void SCF::calc_dpolar(World & world,
		const vecfuncT & ax, const vecfuncT & ay,
		const vecfuncT & bx, const vecfuncT & by,
		const int axis,
		tensorT & Dpolar_total, tensorT & Dpolar_alpha, tensorT & Dpolar_beta)
{
	double Dpolar_average = 0.0;
	double Dpolar_iso = 0.0;

	//START_TIMER(world);
	// derivative density matrix
	functionT drhoa = make_derivative_density(world, amo, aocc, ax, ay);
	functionT drhob;
	if(!param.spin_restricted())
		drhob = make_derivative_density(world, bmo, bocc, bx, by );
	else
		drhob = drhoa;

	functionT drho = drhoa + drhob;

	dpolar(world, Dpolar_alpha, drhoa, axis);
	dpolar(world, Dpolar_beta,  drhob, axis);
	dpolar(world, Dpolar_total, drho,  axis);

	for(int i=0; i<3; ++i)
		Dpolar_total(axis, i) = 0.5 * Dpolar_total(axis, i);

	drhoa.clear(false);
	drhob.clear(false);
	drho.clear(false);


	if (world.rank() == 0 ) {
		printf("Dynamic Polarizability alpha ( Frequency = %.6f, axis %d )\n", param.response_freq(), axis);
		for(unsigned int i=0; i<3; ++i)
			printf(" \t %.6f ", Dpolar_alpha(axis,i));
		printf("\n");


		if(param.nbeta() != 0) {
			printf("Dynamic Polarizability beta ( Frequency = %.6f, axis %d )\n", param.response_freq(), axis);
			for(unsigned int i=0; i<3; ++i)
				printf(" \t %.6f ", Dpolar_beta(axis,i));
			print("\n");
		}

	}

	// last round
	if(axis == 2) {
		//diagonalize
		tensorT V, epolar, eapolar, ebpolar;
		syev(Dpolar_alpha, V, eapolar);
		syev(Dpolar_total, V, epolar);
		if(param.nbeta() != 0)
			syev(Dpolar_beta, V, ebpolar);
		for(unsigned int i=0; i<3; ++i)
			Dpolar_average = Dpolar_average + epolar[i];
		Dpolar_average = Dpolar_average /3.0;
		Dpolar_iso= sqrt(.5)*sqrt( std::pow(Dpolar_alpha(0,0) -  Dpolar_alpha(1,1),2) +
				std::pow(Dpolar_alpha(1,1) -  Dpolar_alpha(2,2),2) +
				std::pow(Dpolar_alpha(2,2) -  Dpolar_alpha(0,0),2));

		if (world.rank() == 0) {
			print("Total Dynamic Polarizability Tensor ( Frequency = ", param.response_freq(), ")\n");
			print(Dpolar_total);
			printf("\tEigenvalues = ");
			printf("\t %.6f \t %.6f \t %.6f \n", epolar[0], epolar[1], epolar[2]);
			printf("\tIsotropic   = \t %.6f \n", Dpolar_average);
			printf("\tAnisotropic = \t %.6f \n", Dpolar_iso);
			printf("\n");
			printf("\n");
		}
	}
	//END_TIMER(world, "Calc D polar");
	print_meminfo(world.rank(), "Calc D polar");
	// end of solving polarizability
}
double SCF::residual_response(World & world, const vecfuncT & x,const  vecfuncT & y,
		const vecfuncT & x_old, const vecfuncT & y_old,
		vecfuncT & rx, vecfuncT & ry)
{
	double residual = 0.0;

	//START_TIMER(world);
	rx = sub(world, x_old, x);
	ry = sub(world, y_old, y);
	std::vector<double> rnormx = norm2s(world, rx);
	std::vector<double> rnormy = norm2s(world, ry);

	double rms, maxval_x, maxval_y;
	vector_stats(rnormx, rms, maxval_x);
	vector_stats(rnormy, rms, maxval_y);
	residual = std::max(maxval_x, maxval_y);

	//END_TIMER(world, "Residual X,Y");
	print_meminfo(world.rank(), "Residual X,Y");

	return residual;
}

/// Calculates the dynamic polarizability of the current system
///
/// This is the only external function for a polarizability calcualtion.
/// The input parameters (such as frequency of perturbing radiation) are
/// all input via the input file (which is parsed on creation of a
/// calculation object, see SCF.h), but are listed here for completeness:
///
/// bool response;                    // response function calculation
/// double response_freq;             // Frequency for calculation response function
/// std::vector<bool> response_axis;  // Calculation protocol
/// bool nonrotate;                   // If true do not molcule orient
/// double rconv;                     // Response convergence
/// double efield;                    // eps for finite field
/// double efield_axis;               // eps for finite field axis

void SCF::polarizability(World & world)
{
	if(world.rank() == 0) {
		print("\n\n\n");
		print(" ------------------------------------------------------------------------------");
		print(" |                MADNESS RESPONSE                                            |");
		print(" ------------------------------------------------------------------------------");
		print(" \n\n");
	}


	// TODO move this  X:axis=0, Y:axis=1, Z:axis=2
	double omega = param.response_freq();

	if (world.rank() == 0) {
		print(" eps_alpha");
		print(aeps);
		if(!param.spin_restricted() && param.nbeta() != 0) print(" eps_beta  = ", beps);

		print(" Frequency for response function (omega)= ", omega);
		print(" Number of alpha orbitals = ", param.nalpha());
		print(" Number of beta orbitals = ", param.nbeta());
	}

	// START_TIMER(world);
	// Green's function
	tensorT akx(param.nalpha());
	tensorT aky(param.nalpha());
	tensorT bkx(param.nbeta());
	tensorT bky(param.nbeta());

	// combine frequency term and eigenvalues
	calc_freq(world, omega, akx, bkx, 1);
	if(omega != 0.0)
		calc_freq(world, omega, aky, bky, -1);
	print_meminfo(world.rank(), "Make frequency term");

	// make density matrix
	functionT arho ;
	functionT brho ;
	functionT rho = make_density_ground(world, arho, brho);

	// vlocal = vnuc + 2*J
	functionT vlocal;
	{
		functionT vnuc;
		// TODO vnuc = potentialmanager->vnuclear();
		if (param.psp_calc()){
			vnuc = gthpseudopotential->vlocalpot();
		}
		else if (param.pure_ae()){
			vnuc = potentialmanager->vnuclear();
		}
		else {
			vnuc = potentialmanager->vnuclear();
			vnuc = vnuc + gthpseudopotential->vlocalpot();
		}
		START_TIMER(world);
		functionT vcoul = apply(*coulop, rho);
		vlocal = vcoul + vnuc;
		END_TIMER(world, "Calc vlocal");
	}
	vlocal.reconstruct();
	vlocal.truncate();

	// BSHOperatorPtr
	std::vector<poperatorT> aopx(param.nalpha());
	std::vector<poperatorT> bopx(param.nbeta());
	std::vector<poperatorT> aopy(param.nalpha());
	std::vector<poperatorT> bopy(param.nbeta());
	make_BSHOperatorPtr(world, akx, bkx, aopx, bopx);

	if(omega != 0.0)
		make_BSHOperatorPtr(world, aky, bky, aopy, bopy);

	tensorT Dpolar_total(3, 3), Dpolar_alpha(3, 3), Dpolar_beta(3, 3);

	double update_residual = 0.0;
	const double rconv = std::max(FunctionDefaults<3>::get_thresh(), param.get<double>("rconv"));
	//        int maxsub_save = param.maxsub();

	for (size_t axis=0; axis<param.response_axis().size(); axis++) {
		if(!param.response_axis()[axis]) continue;

		subspaceT subspace;
		tensorT Q;
		if (world.rank() == 0) {
			this_axis(world, axis);
		}

		// perturbation
		vecfuncT dipoleamo = zero_functions<double,3>(world, param.nalpha());
		vecfuncT dipolebmo = zero_functions<double,3>(world, param.nbeta());

		// make response function x, y
		vecfuncT ax = zero_functions<double,3>(world, param.nalpha());
		vecfuncT ay = zero_functions<double,3>(world, param.nalpha());
		vecfuncT bx = zero_functions<double,3>(world, param.nbeta());
		vecfuncT by = zero_functions<double,3>(world, param.nbeta());

		// old response function
		vecfuncT ax_old = zero_functions<double,3>(world, param.nalpha());
		vecfuncT ay_old = zero_functions<double,3>(world, param.nalpha());
		vecfuncT bx_old = zero_functions<double,3>(world, param.nbeta());
		vecfuncT by_old = zero_functions<double,3>(world, param.nbeta());

		vecfuncT axrhs = zero_functions<double,3>(world, param.nalpha());
		vecfuncT ayrhs = zero_functions<double,3>(world, param.nalpha());
		vecfuncT bxrhs = zero_functions<double,3>(world, param.nbeta());
		vecfuncT byrhs = zero_functions<double,3>(world, param.nbeta());

		vecfuncT aVx;
		vecfuncT bVx;

		vecfuncT aVy;
		vecfuncT bVy;

		// make (dJ-dK)*2*mo
				vecfuncT djkamox;
		vecfuncT djkamoy;

		vecfuncT djkbmox;
		vecfuncT djkbmoy;

		// ri * psi_0
		dipoleamo = calc_dipole_mo(world, amo, axis);
		if(!param.spin_restricted() && param.nbeta() != 0) {
			dipolebmo = calc_dipole_mo(world, bmo, axis);
		}
		else {
			dipolebmo = dipoleamo;
		}

		//guess : drho=rho_0=sum[rho_i]=sum[psi_i^2]
											functionT drhoa = make_derivative_density( world, amo , aocc, dipoleamo, dipoleamo);
											drhoa.reconstruct();
											functionT drhob;
											if(!param.spin_restricted() && param.nbeta() != 0) {
												drhob = make_derivative_density( world, bmo, aocc, dipolebmo, dipolebmo );
												drhob.reconstruct();
											} else {
												drhob = drhoa;
											}
											functionT drho = drhoa + drhob;

											// construct xc operator only once since the ground state density
											// will not change during the iterations.
											XCOperator xcop(world,this,arho,brho);

											// construct xc operator for acting on the perturbed density --
											// use only the LDA approximation
											XCOperator xc_alda(world, "LDA", not param.spin_restricted(), arho, brho);

											for(int iter = 0; iter < param.maxiter(); ++iter) {
												if(world.rank() == 0)
													printf("\nIteration %d at time %.1fs\n\n", iter, wall_time());

												double residual = 0.0;

												//                if (iter > 0 && update_residual < 0.1) {
												//                    //do_this_iter = false;
												//                    param.maxsub = maxsub_save;
												//                }


												if(iter == 0) {
													// iter = 0 initial_guess
													aVx = apply_potential_response(world, dipoleamo, xcop, vlocal,  0);
													djkamox = calc_djkmo(world, xc_alda, dipoleamo, dipoleamo, drho, amo, drhoa,  0);
													axrhs = calc_rhs(world, amo,  aVx, dipoleamo, djkamox);

													if(!param.spin_restricted() && param.nbeta() != 0) {
														bVx = apply_potential_response(world, dipolebmo, xcop, vlocal, 0);
														djkbmox = calc_djkmo(world, xc_alda, dipolebmo, dipolebmo, drho, bmo, drhob,  0);
														bxrhs = calc_rhs(world, bmo,  bVx, dipolebmo, djkbmox);
													}

													if(omega != 0.0) {
														aVy = apply_potential_response(world, dipoleamo, xcop, vlocal, 0);
														djkamoy = calc_djkmo(world, xc_alda, dipoleamo, dipoleamo, drho, amo, drhoa,  0);
														ayrhs = calc_rhs(world, amo,  aVy, dipoleamo, djkamoy);
														if(!param.spin_restricted() && param.nbeta() != 0) {
															bVy = apply_potential_response(world, dipolebmo, xcop, vlocal,  0);
															djkbmoy = calc_djkmo(world, xc_alda, dipolebmo, dipolebmo, drho, bmo, drhob,  0);
															byrhs = calc_rhs(world, bmo,  bVy, dipolebmo, djkbmoy);
														}
													}
												}

												else{
													drhoa = make_derivative_density( world, amo, aocc, ax_old, ay_old );
													drhoa.reconstruct();
													if(!param.spin_restricted() && param.nbeta() != 0) {
														drhob = make_derivative_density( world, bmo, bocc, bx_old, by_old );
														drhob.reconstruct();
													} else {
														drhob = drhoa;
													}
													drho = drhoa + drhob;

													// calculate (dJ-dK)*2*mo
															aVx = apply_potential_response(world, ax_old, xcop, vlocal, 0);
															// make potential * wave function
															djkamox = calc_djkmo(world, xc_alda, ax_old, ay_old, drho, amo, drhoa,  0);
															// axrhs = -2.0 * (aVx + dipoleamo + duamo)
															axrhs = calc_rhs(world, amo,  aVx, dipoleamo, djkamox);

															if(!param.spin_restricted() && param.nbeta() != 0) {
																bVx = apply_potential_response(world, bx_old, xcop, vlocal,  1);
																djkbmox = calc_djkmo(world, xc_alda, bx_old, by_old, drho, bmo, drhob , 1);
																bxrhs = calc_rhs(world, bmo, bVx, dipolebmo, djkbmox);
															}

															if(omega != 0.0) {
																aVy = apply_potential_response(world, ay_old, xcop,  vlocal, 0);
																djkamoy = calc_djkmo(world, xc_alda, ay_old, ax_old,  drho, amo, drhoa, 0);
																// bxrhs = -2.0 * (bVx + dipolebmo + dubmo)
																ayrhs = calc_rhs(world, amo, aVy, dipoleamo, djkamoy);

																if(!param.spin_restricted() && param.nbeta() != 0) {
																	bVy = apply_potential_response(world, by_old, xcop,  vlocal, 1);
																	djkbmoy = calc_djkmo(world, xc_alda, by_old, bx_old, drho, amo, drhob, 1);
																	byrhs = calc_rhs(world, bmo, bVy, dipolebmo, djkbmoy);
																}
															}
															aVx.clear();
															bVx.clear();
															aVy.clear();
															bVy.clear();
															djkamox.clear();
															djkamoy.clear();
															djkbmox.clear();
															djkbmoy.clear();

												}

												//START_TIMER(world);
												// ax_new = G * axrhs;
												calc_response_function(world, ax, aopx, axrhs);
												orthogonalize_response(world, ax, amo);
												truncate(world, ax);
												axrhs.clear();
												if(!param.spin_restricted() && param.nbeta() != 0) {
													// bx_new = G * bxrhs;
													calc_response_function(world, bx, bopx, bxrhs);
													orthogonalize_response(world, bx, bmo);
													truncate(world, bx);
													bxrhs.clear();
												}
												else {
													bx = ax;
												}

												if(omega != 0.0){
													calc_response_function(world, ay, aopy, ayrhs);
													orthogonalize_response(world, ay, amo);
													truncate(world, ay);
													ayrhs.clear();

													if(!param.spin_restricted() && param.nbeta() != 0) {
														calc_response_function(world, by, bopy, byrhs);
														orthogonalize_response(world, by, bmo);
														truncate(world, by);
														byrhs.clear();
													}
													else {
														by = ay;
													}
												}
												else {
													ay = ax;
													by = bx;
												}
												//END_TIMER(world, "Make response func");
												print_meminfo(world.rank(), "Make response func");

												if(iter > 0) {
													// START_TIMER(world);
													residual = 0.0;

													vecfuncT rax = zero_functions<double,3>(world, param.nalpha()); //residual alpha x
													vecfuncT ray = zero_functions<double,3>(world, param.nalpha()); //residual alpha y
													vecfuncT rbx = zero_functions<double,3>(world, param.nbeta());  //residual beta x
													vecfuncT rby = zero_functions<double,3>(world, param.nbeta());  //residual beta y

													double aresidual =  residual_response(world, ax, ay, ax_old, ay_old, rax, ray);
													double bresidual = 0.0;
													world.gop.fence();
													if(!param.spin_restricted() && param.nbeta() != 0) {
														bresidual = aresidual + residual_response(world, bx, by, bx_old, by_old, rbx, rby);
														residual = std::max(aresidual, bresidual);
														world.gop.fence();
													}
													else {
														residual = aresidual;
													}

													if (world.rank() == 0)
														print("\nresiduals_response (first) = ", residual);
													residual = 0.0;

													double nx,ny;
													////////UPDATE
													nx=norm2(world, ax);
													if (world.rank() == 0)
														print("CURRENT_X_norm2() = ", nx);
													update_response_subspace(world, ax, ay, bx, by, rax, ray, rbx, rby, subspace, Q, update_residual);

													nx = norm2(world, ax);
													ny = norm2(world, ay);
													if (world.rank() == 0) {
														print("new X (alpha) norm2() = ", nx);
														print("new Y (alpha) norm2() = ", ny);
													}

													aresidual = residual_response(world, ax, ay, ax_old, ay_old, rax, ray);
													bresidual = 0.0;

													if(!param.spin_restricted() && param.nbeta() != 0) {
														bresidual = residual_response(world, bx, by, bx_old, by_old, rbx, rby);
														residual = std::max(aresidual, bresidual);
													}
													else residual = aresidual;

													double thresh = rconv *(param.nalpha() + param.nbeta())*2;
													if (world.rank() == 0) {
														print("\nresiduals_response (final) = ", residual);
														print("rconv *(param.nalpha + param.nbeta)*2", thresh);
													}

													//  END_TIMER(world, "Update response func");
													print_meminfo(world.rank(), "Update response func");

													if( residual < (rconv *(param.nalpha() + param.nbeta())*2))
													{
														if (world.rank() == 0) {
															print("\nConverged response function!!\n");
															print("\n\n\n");
															print(" ------------------------------------------------------------------------------");
															print(" |                  MADNESS CALCULATION POLARIZABILITY                        |");
															print(" ------------------------------------------------------------------------------");
															print(" \n\n");
														}
														break;
													}
												}

												ax_old = ax;
												ay_old = ay;
												bx_old = bx;
												by_old = by;
												ax.clear();
												ay.clear();
												bx.clear();
												by.clear();

											} //end iteration
											//END_TIMER(world, "Make response func");
											print_meminfo(world.rank(), "Make response func");

											calc_dpolar(world, ax_old, ay_old, bx_old, by_old, axis, Dpolar_total, Dpolar_alpha, Dpolar_beta);

#if  0
//hyper polarizability
											for (int p=0; p < ax_old.size(); p++){
												axx.push_back(ax_old[p]);
												ayx.push_back(ay_old[p]);
												if(!param.spin_restricted && param.nbeta != 0) {
													bxx.push_back(bx_old[p]);
													byx.push_back(by_old[p]);
												}
											}
#endif 
ax_old.clear();
ay_old.clear();
bx_old.clear();
by_old.clear();

dipoleamo.clear();
dipolebmo.clear();
	} //end axis

}
//vama polarizability


}


