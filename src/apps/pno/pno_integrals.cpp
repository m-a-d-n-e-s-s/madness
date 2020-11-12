/*
 * pno_integrals.hpp
 *
 *  Created on: Feb. 28, 2020
 *      Author: jsk


 CURRENTLY TO PLAY AROUND
 A LOT IS HARDCODED FOR H2
 SO BE CAREFUL WITH DIFFERENT INPUTS


 */


//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/world/print.h>
#include <iomanip>
#include <madness/mra/vmra.h>
#include <chem/SCF.h>
#include <chem/nemo.h>
#include <chem/PNO.h>
#include <vector>
#include <fstream>
#include <tuple>

using namespace madness;

// DEFINE PARAMETER TAGS FOR THE INPUT FILE
const std::string TAG_PNO = "pno";
const std::string TAG_F12 = "f12";
const std::string TAG_PNOInt = "pnoint";
const std::string TAG_CP = "computeprotocol";

// this needs to be added to include
#include "NumCpp.hpp"

template<typename T1, typename T2>
std::ostream& operator << (std::ostream& os, const std::pair<T1,T2>& v){
    os << "(" << v.first << "," << v.second << ")";
    return os;
}


template<typename T>
std::ostream& operator << (std::ostream& os, const std::vector<T>& v){
    os << "[";
    for(auto i=0; i<v.size();++i){
        os << v[i] << " ";
    }
    os << "]";
    return os;
}

// Function for orthogonalization of a basis
std::vector<real_function_3d> orthonormalize_basis(const std::vector<real_function_3d> &basis, const std::string orthogonalization,
		const double thresh, World& world, const bool print_out) {

	// compute overlap, to be passed in orthonormalization routines and potentially printed
	auto S = madness::matrix_inner(world, basis, basis, true);

	auto out_basis = basis;
	if (orthogonalization == "cholesky"){
		out_basis = madness::orthonormalize_cd(basis, S);
	}else if(orthogonalization == "symmetric"){
		out_basis = madness::orthonormalize_symmetric(basis, S);
	}else if(orthogonalization == "canonical") {
		out_basis = madness::orthonormalize_canonical(basis, S, thresh);
	}else if(orthogonalization == "rr_cholesky") {
		if(world.rank()==0) std::cout << std::endl <<  "Be cautious:" << std::endl;
		if(world.rank()==0) std::cout << "Rank reduced cholesky changes ordering in basis via pivoting. Make sure, that this does not interfere with your application (e.g. active-space)." << std::endl;
		out_basis = madness::orthonormalize_rrcd(basis, S, thresh);
	}else{
		MADNESS_EXCEPTION("unknown orthonormalization method", 1);
	}

	if (print_out) {
		const auto new_S = madness::matrix_inner(world, out_basis, out_basis);
		if(world.rank()==0) {
			std::cout << "Overlap Matrix before " << orthogonalization << "\n";
			std::cout << S;
			std::cout << "Overlap Matrix after " << orthogonalization << "\n";
			std::cout << new_S;
		}
	}

	return out_basis;
}

// --------------------------------------------------------------------------------------------- //
int main(int argc, char** argv) {
	madness::initialize(argc, argv);
	{
		World world(SafeMPI::COMM_WORLD);
		if (world.rank() == 0) printf("starting at time %.1f\n", wall_time());
		const double time_start = wall_time();
		std::cout.precision(6);

		startup(world,argc,argv,true);
		print_meminfo(world.rank(), "startup");

		// Get the name of the input file (if given)
		const std::string input = (argc > 1) ? argv[1] : "input";

		if(world.rank()==0){
			std::cout << "\n\n";
			std::cout << "-------------------------------------------------------------------------------------\n";
			std::cout << "SOLVING MRA-PNO-F12 as described in \n";
			std::cout << "J.S. Kottmann, F.A. Bischoff, E.F. Valeev\n";
			std::cout << "Direct determination of optimal pair-natural orbitals in a real-space representation:\n";
			std::cout << "the second-order Møller-Plesset energy\n";
			std::cout << "Journal of Chemical Physics ... 2020\n";
			std::cout << "-------------------------------------------------------------------------------------\n";
			std::cout << "\n\n";

			std::cout << "This script will run PNO-MP2 and print out tensors in binary\n";
			//std::cout << "Call as: pno_integrals inputfile orthogonalization basis_size";
			std::cout << "Call as: pno_integrals inputfile";
			std::cout << "input is " << input << "\n";
			//std::cout << "orthogonalization is " << orthogonalization << "\n";
			//std::cout << "basis size is " << basis_size << "\n";
			//std::cout << "using CABS-option: " << cabs_option << "\n";
			//std::cout << "with pno-OBS-size: " << pno_obs_size << "\n";

			//std::cout << "only diag is " << only_diag << "\n";
			//std::cout << "cherry_pick is " << cherry_pick << "\n";
		}

		// Compute the SCF Reference
		const double time_scf_start = wall_time();
		std::shared_ptr<SCF> calc(new SCF(world, input));
		Nemo nemo(world, calc, input);
		nemo.get_calc()->param.print();
		const double scf_energy = nemo.value();
		if (world.rank() == 0) print("nemo energy: ", scf_energy);
		if (world.rank() == 0) printf(" at time %.1f\n", wall_time());
		const double time_scf_end = wall_time();

		// Compute MRA-PNO-MP2-F12
		const double time_pno_start = wall_time();
		PNOParameters parameters(world,input,nemo.get_calc()->molecule,TAG_PNO);
		F12Parameters paramf12(world, input, parameters, TAG_F12);
		PNOIntParameters paramsint(world, input, parameters, TAG_PNOInt);
		paramsint.print("PNO Integrals evaluated as:\npnoint","end");
		PNO pno(world, nemo, parameters, paramf12);
		pno.solve();
		const double time_pno_end = wall_time();


		if(world.rank()==0){
			std::cout << std::setfill(' ');
			std::cout << "\n\n\n";
			std::cout << "--------------------------------------------------\n";
			std::cout << "MRA-PNO-MP2-F12 ended \n";
			std::cout << "--------------------------------------------------\n";
			std::cout << std::setw(25) << "time scf" << " = " << time_scf_end - time_scf_start << "\n";
			std::cout << std::setw(25) << "energy scf" << " = " << scf_energy << "\n";
			std::cout << "--------------------------------------------------\n";
		}

		if(world.rank()==0){
			std::cout << "restarting PNO to reload pairs that converged before and were frozen\n";
		}
		pno.param.set_user_defined_value<std::string>("restart", "all");
		pno.param.set_user_defined_value<std::string>("no_opt", "all");
		pno.param.set_user_defined_value<std::string>("no_guess", "all");
		pno.param.set_user_defined_value<std::string>("adaptive_solver", "none");
		std::vector<PNOPairs> all_pairs;
		pno.solve(all_pairs);

		double mp2_energy = 0.0;
		if(world.rank()==0) std::cout<< std::setw(25) << "time pno" << " = " << time_pno_end - time_pno_start << "\n";
		for(const auto& pairs: all_pairs){
			if(pairs.type == MP2_PAIRTYPE){
				mp2_energy = pairs.energies.total_energy();
			}
			std::pair<size_t, size_t> ranks= pno.get_average_rank(pairs.pno_ij);
			if(world.rank()==0){
				std::string name;
				std::stringstream ss;
				ss << pairs.type;
				ss >> name;
				std::cout<< std::setw(25) << "energy "+name << " = " << pairs.energies.total_energy() << "\n";
				std::cout<< std::setw(25) << "average pno rank " + name << " = " << ranks.first << "\n";
				std::cout<< std::setw(25) << "max pno rank " + name << " = " << ranks.second << "\n";
			}
		}
		if(world.rank()==0 and mp2_energy != 0.0){
			std::cout << "--------------------------------------------------\n";
			std::cout<< std::setw(25) << "energy(total)" << " = " << scf_energy + mp2_energy << "\n";
			std::cout << "--------------------------------------------------\n";
			std::cout << "\n\n\n";
		}

		// compute orthogonalized mp2 basis and print out hamiltonian tensors
		std::cout << std::setprecision(8);
		std::cout << std::fixed;
		std::cout << std::showpos;

		const std::string orthogonalization = paramsint.orthogonalization();
		if (world.rank()==0) std::cout << "Orthonormalization technique used: " << orthogonalization << std::endl;
		const bool orthogonalize = orthogonalization != "none";
		const double h_thresh = 1.e-4;
		const double thresh = parameters.thresh();

		const std::string cabs_option = paramsint.cabs_option();
		bool cabs_switch = false;
		if (cabs_option=="pno" || cabs_option=="gbs" || cabs_option=="mixed")
			cabs_switch = true;
		if (world.rank()==0) {
			if (cabs_switch) std::cout << "CABS option " << cabs_option << std::endl;
			else if (!cabs_switch) std::cout << "No CABS used." << std::endl;
		}
		const bool only_diag = paramsint.only_diag();
		const int basis_size = paramsint.n_pno();
		int pno_cabs_size = paramsint.pno_cabs_size();

		if(world.rank()==0) std::cout << "Tightening thresholds to 1.e-6 for post-processing\n";
		FunctionDefaults<3>::set_thresh(1.e-6);

		vecfuncT reference = nemo.get_calc()->amo;
		vecfuncT obs_pnos;
		std::vector<real_function_3d> rest_pnos;
		std::vector<double> occ;
		std::vector<double> rest_occ;
		std::vector<std::pair<size_t,size_t>> pno_ids;
		std::vector<std::pair<size_t,size_t>> rest_ids;
		std::string name;
		for(auto& pairs: all_pairs){
			const auto& pno_ij = pairs.pno_ij;
			const auto& rdm_evals = pairs.rdm_evals_ij;
			const bool is_gs = pairs.type == MP2_PAIRTYPE;

			if (not is_gs){
				const auto& x = pairs.cis.x;
				reference.insert(reference.end(), x.begin(), x.end());
				name = "ex" + std::to_string(pairs.cis.number);
			}else{
				name = "gs";
			}

			std::vector<real_function_3d> all_current_pnos;
			// collect PNOs from all pairs and sort by occupation number, keeping pair information via name
			for(ElectronPairIterator it=pno.pit();it;++it){
				if (only_diag and not it.diagonal()){
					if(world.rank()==0) std::cout << "skipping pair " << it.name() << "\n";
					continue;
				}else{
					const auto& pair = pno_ij[it.ij()];
					all_current_pnos.insert(all_current_pnos.end(), pair.begin(), pair.end());
					for (auto ii=0; ii<rdm_evals[it.ij()].size();++ii){
						occ.push_back(rdm_evals[it.ij()][ii]);
						pno_ids.push_back(std::make_pair(it.i(),it.j()));  // for each eigenvalue ~ PNO, store pair affiliation
					}
				}
			}

			std::vector<std::tuple<double, real_function_3d, std::pair<size_t,size_t> > > zipped;
			for (auto i=0; i< all_current_pnos.size(); ++i){
				zipped.push_back(std::make_tuple(occ[i], all_current_pnos[i], pno_ids[i]));
			}

			std::sort(zipped.begin(), zipped.end(), [](const auto& i, const auto& j) { return std::get<0>(i) > std::get<0>(j); });

			std::vector<double> unzipped_first;
			std::vector<real_function_3d> unzipped_second;
			std::vector<std::pair<size_t,size_t> > unzipped_third;
			for (auto i=0; i<basis_size;++i){
				unzipped_first.push_back(std::get<0>(zipped[i]));
				unzipped_second.push_back(std::get<1>(zipped[i]));
				unzipped_third.push_back(std::get<2>(zipped[i]));
			}
			occ = unzipped_first;
			all_current_pnos = unzipped_second;
			pno_ids = unzipped_third;

			// if CABS shall be used, keep the remaining PNOs
			unzipped_first.clear();
			unzipped_second.clear();
			unzipped_third.clear();
			if (cabs_option=="pno" || cabs_option=="mixed") {
				if (pno_cabs_size==-1) pno_cabs_size = int(zipped.size());
				for (int i=basis_size; i<pno_cabs_size; ++i) {
					unzipped_first.push_back(std::get<0>(zipped[i]));
					unzipped_second.push_back(std::get<1>(zipped[i]));
					unzipped_third.push_back(std::get<2>(zipped[i]));
				}
			}
			rest_occ = unzipped_first;
			rest_pnos = unzipped_second;
			rest_ids = unzipped_third;


			obs_pnos.insert(obs_pnos.end(), all_current_pnos.begin(), all_current_pnos.end());
		}

		// reference projector (not automatically fullfilled for CIS)
		// projection is to keep the reference and CIS orbitals untouched in the orthogonalization
		madness::QProjector<double, 3> Q(world, reference);
		obs_pnos = Q(obs_pnos);

		auto basis = obs_pnos;

		// compute overlap of all PNOs before orthogonalization
		if (paramsint.print_pno_overlap()) {
			const auto S = madness::matrix_inner(world, obs_pnos, obs_pnos, true);
			if(world.rank()==0) std::cout << "Overlap Matrix of all PNOs:\n";
			for (int i=0;i<obs_pnos.size();++i){
				for (int j=0;j<obs_pnos.size();++j){
					if(world.rank()==0) std::cout << S(i,j) << " ";
				}
				if(world.rank()==0) std::cout << "\n";
			}
		}

		if(world.rank()==0){
			std::cout << "Before cherry pick" << std::endl;
			std::cout << "rank is " << world.rank() << "\n";
			std::cout << "all used occupation numbers:\n" << occ << std::endl;
			std::cout << "corresponding to active pairs:\n" << pno_ids << std::endl;
			if (cabs_option=="pno" || cabs_option=="mixed") {
				std::cout << "add remaining occupation numbers for cabs:\n" << rest_occ << std::endl;
				std::cout << "corresponding to pairs:\n" << rest_ids << std::endl;
			}
		}

		// cherry pick PNOs
		std::vector<int> cherry_pick = paramsint.cherry_pick();
		if(not cherry_pick.empty()){
			vecfuncT cp;
			std::vector<double> cp_occ;
			std::vector<std::pair<size_t,size_t> > cp_pno_ids;
			if(world.rank()==0){
				std::cout << "Cherry picking orbitals: " << cherry_pick << " from pno basis." << std::endl;
			}
			for(auto i: cherry_pick){ // TODO: whatever does not get cherry-picked, add to CABS, if CABS!
				cp.push_back(basis[i]);
				cp_occ.push_back(occ[i]);
				cp_pno_ids.push_back(pno_ids[i]);
			}
			basis = cp;
			occ = cp_occ;
			pno_ids = cp_pno_ids;
		}
		if(world.rank()==0){
			std::cout << "After cherry pick" << std::endl;
			std::cout << "all used occupation numbers:\n" << occ << std::endl;
			std::cout << "corresponding to active pairs:\n" << pno_ids << std::endl;
			if (cabs_option=="pno" || cabs_option=="mixed") {
				std::cout << "add remaining occupation numbers for cabs:\n" << rest_occ << std::endl;
				std::cout << "corresponding to pairs:\n" << rest_ids << std::endl;
			}
		}

		// orthogonalize orbital basis
		if (orthogonalize){
			basis = orthonormalize_basis(basis, orthogonalization, thresh, world, paramsint.print_pno_overlap());
		}

		// will include CIS orbitals for excited states
		if(world.rank()==0) std::cout << "Adding " << reference.size() << " Reference orbitals\n";
		basis.insert(basis.begin(), reference.begin(), reference.end());

		// include virtual orbitals if demanded
		// not the most elegant solution ... but lets see if we need this first
		vecfuncT virtuals;
		std::vector<double> veps;
		if (paramsint.n_virt() > 0){
			const auto refsize = reference.size();
			SCF calcx(world, input);
			calcx.param.set_user_defined_value("nvalpha", paramsint.n_virt());
			calcx.param.set_user_defined_value("nvbeta", paramsint.n_virt());
			calcx.param.set_user_defined_value("restart", false);
			calcx.param.set_user_defined_value("no_compute", false);
	        calcx.param.set_user_defined_value("nmo_alpha",calc->param.nalpha() + paramsint.n_virt());
	        calcx.param.set_user_defined_value("nmo_beta",calc->param.nbeta() + paramsint.n_virt());
	        calcx.initial_guess(world);
			calcx.param.print();
			calcx.amo = calc->amo;
			auto F = Fock(world, calc.get());
			auto Fmat = F(calc->ao, calc->ao);
			Tensor<double> U, evals;
			syev(Fmat, U, evals);
			auto vguess = madness::transform(world, calc->ao, U);
			for (auto i=0;i<paramsint.n_virt();++i){
				calcx.amo.push_back(vguess[i]);
			}
			calcx.aocc = tensorT(calcx.param.nmo_alpha());
			for (int i = 0; i < calcx.param.nalpha(); ++i)
				calcx.aocc[i] = 1.0;

			MolecularEnergy E(world, calcx);
			double energy=E.value(calcx.molecule.get_all_coords().flat());
			for (auto i=refsize; i<refsize+paramsint.n_virt(); ++i){
				virtuals.push_back(calcx.amo[i]);
				veps.push_back(calcx.aeps(i));
			}
			auto QB = QProjector<double,3>(world, basis);
			virtuals = QB(virtuals);
			madness::normalize(world, virtuals);
			basis.insert(basis.end(), virtuals.begin(), virtuals.end());
			if(world.rank() ==0) std::cout << "added " << paramsint.n_virt() << " virtuals\n";
		}

		if (parameters.save_pnos()) {
			for(auto i=0; i<basis.size(); ++i){
				madness::save(basis[i], "mra_orbital_"+std::to_string(i));
			}
		}

		// Save pair information to file after having picked the cherries
		std::ofstream pairwriter ("pnoinfo.txt", std::ofstream::out|std::ofstream::trunc);
		pairwriter << "MADNESS MRA-PNO INFORMATION" << std::endl;
		pairwriter << "pairinfo=";
		for(auto i=0; i<reference.size(); ++i) {
			pairwriter << i << ",";
		}
		// print total indices of pairs
		const auto nfreeze = parameters.freeze();
		for(auto i=0; i<pno_ids.size(); ++i) {
			pairwriter << pno_ids[i].first+nfreeze << "." << pno_ids[i].second+nfreeze;
			if (!(i==pno_ids.size()-1))pairwriter << ",";
		}
		// print info of virtuals if there
		for (auto i=0; i<virtuals.size(); ++i){
			pairwriter << ",";
			pairwriter << -(i+1);
		}

		pairwriter << std::endl;
		pairwriter << "nuclear_repulsion=" << nemo.get_calc()->molecule.nuclear_repulsion_energy() << std::endl;
		pairwriter << "occinfo=";
		for (auto i=0; i<reference.size();++i){
			pairwriter << nemo.get_calc()-> aeps(i);
			pairwriter << ",";
		}
		for (auto i=0; i<occ.size();++i){
			pairwriter << occ[i];
			if (!(i==occ.size()-1))pairwriter << ",";
		}
		for (auto i=0; i<veps.size();++i){
			pairwriter << ",";
			pairwriter << veps[i];
		}
		pairwriter.close();


		auto size_obs = basis.size();

		// if desired, build a CABS using either a Gaussian basis set or the remaining PNOs
		if (cabs_switch) {
			if(world.rank()==0) std::cout << "Trying to use external CABS basis..." << std::endl;
			std::vector<real_function_3d> cabs;
			// Get CABS from rest-pnos or external basis
			if (cabs_option == "gbs" || cabs_option == "mixed") {
				cabs = pno.f12.read_cabs_from_file(paramf12.auxbas_file());
			}
			if (cabs_option == "pno" || cabs_option == "mixed") {
				cabs.insert(cabs.begin(), rest_pnos.begin(), rest_pnos.end());
			}
			// set up CABS
			if (!cabs.empty()) {
				if(world.rank()==0) std::cout << "Found CABS..." << std::endl;
				MyTimer time2 = MyTimer(world).start();
				// Project out reference
				//cabs = Q(cabs);
				// Project out {pno} + ref
				auto pno_plus_ref = basis;
				if(world.rank()==0) std::cout << "\tProject out PNO + ref" << std::endl;
				madness::QProjector<double, 3> Qpno(world, pno_plus_ref);
				cabs = Qpno(cabs);
				// Orthonormalize {cabs}
				if(world.rank()==0) std::cout << "\tOrthonormalize CABS" << std::endl;
				cabs = orthonormalize_basis(cabs, paramsint.cabs_orthogonalization(), paramsint.cabs_thresh(),
						world, paramsint.print_pno_overlap());
				time2.stop().print("Made pair specific ABS from PNOS and " + std::to_string(cabs.size()) + " functions");

				// Truncate cabs
				madness::truncate(world, cabs, thresh);

				// Merge {cabs} + {pno}
				if(world.rank()==0) std::cout << "Adding {cabs} to {pno+ref}.\n";

				basis.clear();
				basis = cabs;
				basis.insert(basis.begin(), pno_plus_ref.begin(), pno_plus_ref.end());
			}
			else if (cabs.empty()) {
				if(world.rank()==0) std::cout << "CABS basis set not found..." << std::endl; 
			}
		}

		auto size_full = basis.size();
		if(world.rank()==0) std::cout << "Size before adding CABS: " << size_obs << ".\n";
		if(world.rank()==0) std::cout << "Size after adding CABS: " << size_full << ".\n";

		// Build tensors
		// define Coulomb and F12-operator (as SlaterOperator exp[-mu*r] )
		auto gop = std::shared_ptr < real_convolution_3d > (madness::CoulombOperatorPtr(world, 1.e-6, parameters.op_thresh()));
		auto fop = SlaterOperator(world, paramsint.gamma(), 1.e-6, parameters.op_thresh());

		madness::Tensor<double> g(size_full, size_obs, size_full, size_obs); // using mulliken notation since thats more efficient to compute here: Tensor is (pq|g|rs) = <pr|g|qs>
		madness::Tensor<double> f(size_full, size_obs, size_full, size_obs);

		std::vector<vecfuncT> PQ;
		for (const auto& x : basis){
			PQ.push_back(madness::truncate(x*basis,thresh));
		}
		std::vector<vecfuncT> GPQ;
		std::vector<vecfuncT> FPQ;
		for (const auto& x : basis){
			GPQ.push_back(madness::truncate(madness::apply(world, *gop, madness::truncate(x*basis,thresh)), thresh));
			if (cabs_switch){
				FPQ.push_back(madness::truncate(madness::apply(world,  fop, madness::truncate(x*basis,thresh)), thresh));
			}
		}

		auto J = madness::Coulomb(world, &nemo);
		auto K = madness::Exchange<double, 3>(world, &nemo, 0);
		auto Jmat = J(basis, basis);
		auto Kmat = K(basis, basis);

		// Build tensors using symmetries (both Coulomb  and SlaterOperator are symmetric)  TODO maybe in separate function...
		// pqrs = qprs = pqsr = qpsr && pqrs = rspq for full-size tensors
		int non_zero=0, non_zero_f=0; // TODO check counting! probably far from correct as of now
		if(world.rank() ==0 ) std::cout << "Compute G and F Tensor:\n";
		for (int p=0; p<size_full; p++){
			for (int q=0; q<size_obs; q++){
				if (PQ[p][q].norm2() < h_thresh) continue;
				// (NN|NN)-part, full symmetries (N<=size_obs)
				if ((p < size_obs) && (q < size_obs)) {
					if (p>=q) {
						for (int r=0; r<size_obs; r++){
							for (int s=0; s<=r; s++){

								if(paramsint.hardcore_boson()){
									const auto c1 = p==q and r==s;
									const auto c2 = p==r and q==s;
									const auto c3 = p==s and q==r;
									if (not(c1 or c2 or c3)){
										continue;
									}
								}

								if ((p*size_obs+q >= r*size_obs+s)) {
									// g-tensor
									if (GPQ[r][s].norm2() >= h_thresh){
										g(p,q,r,s) = PQ[p][q].inner(GPQ[r][s]);
										// symm
										g(r,s,p,q) = g(p,q,r,s);
										if(std::fabs(g(p,q,r,s)) > h_thresh ){
											if(world.rank()==0 and basis.size() < 3) std::cout << " g " << p  << " "<<  q  << " "<<  r  << " "<<  s << " = " << g(p,q,r,s) << "\n";
											non_zero += (p!=r && q!=s)? 2 : 1;
										}
									}
									// f12-tensor, only if cabs are used
									if (cabs_switch and FPQ[r][s].norm2() >= h_thresh) {
										f(p,q,r,s) = PQ[p][q].inner(FPQ[r][s]);
										// symm
										f(r,s,p,q) = f(p,q,r,s);
										if(std::fabs(f(p,q,r,s)) > h_thresh ){
											non_zero_f += (p!=r && q!=s)? 2 : 1;
										}
									}
								}
							}
						}
					}
					// (NN|MN), M>size_obs
					for (int r=size_obs; r<size_full; r++){
						for (int s=0; s<size_obs; s++){
							// g-tensor
							if (GPQ[r][s].norm2() >= h_thresh)  {
								g(p,q,r,s) = PQ[p][q].inner(GPQ[r][s]);
								if(std::fabs(g(p,q,r,s)) > h_thresh ){
									if(world.rank()==0 and basis.size() < 3) std::cout << " g " << p  << " "<<  q  << " "<<  r  << " "<<  s << " = " << g(p,q,r,s) << "\n";
									++non_zero;
								}
							}
							// f12-tensor
							if (cabs_switch and FPQ[r][s].norm2() >= h_thresh) {
								//f(p,q,r,s) = PQ[p][q].inner(FPQ[r][s]);
								f(p,q,r,s) = PQ[p][q].inner(FPQ[r][s]);
								if(std::fabs(f(p,q,r,s)) > h_thresh ){
									++non_zero_f;
								}
							}
						}
					}
				}
				else if (p >= size_obs){
					// (MN|MN), M>N
					for (int r=size_obs; r<size_full; r++){
						for (int s=0; s<=size_obs; s++){
							if ((p*size_full+q >= r*size_full+s)) {
								// g-tensor
								if (GPQ[r][s].norm2() >= h_thresh)  {
									g(p,q,r,s) = PQ[p][q].inner(GPQ[r][s]);
									g(r,s,p,q) = g(p,q,r,s);
									if(std::fabs(g(p,q,r,s)) > h_thresh ){
										if(world.rank()==0 and basis.size() < 3) std::cout << " g " << p  << " "<<  q  << " "<<  r  << " "<<  s << " = " << g(p,q,r,s) << "\n";
										++non_zero;
									}
								}
								// f12-tensor
								if (cabs_switch and FPQ[r][s].norm2() >= h_thresh) {
									f(p,q,r,s) = PQ[p][q].inner(FPQ[r][s]);
									f(r,s,p,q) = f(p,q,r,s);
									if(std::fabs(f(p,q,r,s)) > h_thresh ){
										++non_zero_f;
									}
								}
							}
						}
					}
					// (MN|NN)
					for (int r=0; r<size_obs; r++){
						for (int s=0; s<size_obs; s++){
							// g-tensor
							if (GPQ[r][s].norm2() >= h_thresh)  {
								g(p,q,r,s) = PQ[p][q].inner(GPQ[r][s]);
								if(std::fabs(g(p,q,r,s)) > h_thresh ){
									if(world.rank()==0 and basis.size() < 3) std::cout << " g " << p  << " "<<  q  << " "<<  r  << " "<<  s << " = " << g(p,q,r,s) << "\n";
									++non_zero;
								}
							}
							// f12-tensor
							if (cabs_switch and FPQ[r][s].norm2() >= h_thresh) {
								f(p,q,r,s) = PQ[p][q].inner(FPQ[r][s]);
								if(std::fabs(f(p,q,r,s)) > h_thresh ){
									++non_zero_f;
								}
							}
						}
					}
				}
			}
		}
		// Assemble full g,f from symmetries (misses pqrs->qpsr)
		if (world.rank()==0) std::cout << "Assemble g,f using symmetries..." << std::endl;
		for (int p=0; p<size_obs; p++) {
			for (int q=0; q<=p; q++) {
				if (PQ[p][q].norm2() < h_thresh) continue; for (int r=0; r<size_obs; r++){
					for (int s=0; s<=r; s++){
						// g-tensor
						g(p,q,s,r) = g(p,q,r,s);
						g(q,p,r,s) = g(p,q,r,s);
						g(q,p,s,r) = g(p,q,r,s);
						if      (p==q && r==s) non_zero += 0;
						else if (p==q && r!=s) non_zero += 1;
						else if (p!=q && r==s) non_zero += 1;
						else if (p!=q && r!=s) non_zero += 3;
						non_zero += (p==q && r==s) ? 0 : 3;
						// f12-tensor
						if (cabs_switch) {
							f(p,q,s,r) = f(p,q,r,s);
							f(q,p,r,s) = f(p,q,r,s);
							f(q,p,s,r) = f(p,q,r,s);
							if      (p==q && r==s) non_zero_f += 0;
							else if (p==q && r!=s) non_zero_f += 1;
							else if (p!=q && r==s) non_zero_f += 1;
							else if (p!=q && r!=s) non_zero_f += 3;
						}
					}
				}
			}
		}


		int non_zero_h = 0;

		Tensor<double> h;
		{
			auto T = madness::Kinetic<double, 3>(world);
			auto V = madness::Nuclear(world, &nemo);
			h = T(basis,basis) + V(basis,basis);
		}
		for (int i=0;i<basis.size();++i){
			for (int j=0;j<basis.size();++j){
				if(std::fabs(h(i,j)) > h_thresh){
					if(world.rank()==0 and basis.size() < 3) std::cout << " h " << i << " "<< j << " "<< h(i,j) << "\n";
					++non_zero_h;
				}
			}
		}

		// print out number of non-zero elements
		if(world.rank()==0){
			std::cout << "non zero elements:\n g : " << non_zero << "\n h : " << non_zero_h;
			if (cabs_switch) std::cout << "non zero elements:\n f : " << non_zero_f;
			std::cout << std::endl;
		}

		// save integrals to binary files
		h = h.flat();
		nc::NdArray<double> hh(h.ptr(), h.size(), 1);
		hh.tofile(name+"_htensor.bin", "");

		g = g.flat();
		nc::NdArray<double> gg(g.ptr(), g.size(), 1);
		gg.tofile(name+"_gtensor.bin", "");


		if (cabs_switch) {
			f = f.flat();
			nc::NdArray<double> ff(f.ptr(), f.size(), 1);
			ff.tofile(name+"_f12tensor.bin", "");
		}

		if (not orthogonalize || paramsint.print_pno_overlap()){
			auto S = madness::matrix_inner(world, basis, basis, true);
			if(paramsint.print_pno_overlap()) {
				if(world.rank()==0) std::cout << "Overlap over whole basis\n" << S << "\n";
			}
			if (not orthogonalize) {
				S = S.flat();
				nc::NdArray<double> gg(S.ptr(), S.size(), 1);
				gg.tofile(name+"_overlap.bin", "");
			}
		}

		auto Fop =  madness::Fock(world, &nemo);
		auto F = Fop(basis, basis);
		if(world.rank()==0) std::cout << "F\n" << F << "\n";

		world.gop.fence();
		if (world.rank() == 0) printf("finished at time %.1f\n", wall_time());

		print_stats(world);
	}
	finalize();

	return 0;
}
