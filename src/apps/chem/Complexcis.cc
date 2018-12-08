/*
 * Complexcis.cpp
 *
 *  Created on: 21 Nov 2018
 *      Author: fbischoff
 */

#include <chem/Complexcis.h>
#include <chem/GuessFactory.h>


namespace madness {

double Complex_cis::value() {
	std::vector<std::vector<complex_function_3d> > xa=make_guess("alpha");
	std::vector<std::vector<complex_function_3d> > xb=make_guess("beta");
	return 0.0;
}

std::vector<std::vector<complex_function_3d> > Complex_cis::read_guess(const std::string spin) const {

	int nmo= (spin=="alpha") ? nemo.cparam.nalpha : nemo.cparam.nbeta;
	std::vector<real_function_3d> real_mo=zero_functions<double,3>(world,nmo);

//	// load the converged orbitals
//    for (std::size_t imo = 0; imo < nmo; ++imo) {
//    	print("loading mos ",spin,imo);
//    	load(real_mo[imo], "nemo" + stringify(imo));
//    }
//    return convert<double,double_complex,3>(world,real_mo);

}

std::vector<std::vector<complex_function_3d> > Complex_cis::make_guess(const std::string spin) const {

	std::vector<complex_function_3d> virtuals;
	std::vector<complex_function_3d> active_mo=nemo.amo;
	std::vector<coord_3d> centers = guessfactory::compute_centroids(active_mo);


	// prepare the list of excitation operators and copied seeds
	std::vector<std::pair<std::vector<complex_function_3d>, std::string> > exlist;
	{
//		std::vector<std::string> exop_strings=cis_param.exops();
		std::vector<std::string> exop_strings=(guessfactory::make_predefined_exop_strings(cis_param.guess_excitation_operators()));
		for(const auto ex: exop_strings){
			std::vector<complex_function_3d> cseed=copy(world,active_mo,false);
			exlist.push_back(std::make_pair(cseed,ex));
		}
	}
	world.gop.fence();
	std::cout << "will create " << exlist.size()*centers.size() << " virtuals, from " << centers.size()
			<< " seeds and " << exlist.size() << " excitation operators"   << std::endl;

	// create the virtuals by unary operations: multiply excitation operators with seeds
	for(auto it:exlist){
		virtuals=append(virtuals,guessfactory::apply_trigonometric_exop(it.first,it.second,centers,false));
	}
	world.gop.fence();

	// remove linear dependencies
	const size_t spre=virtuals.size();
	virtuals=orthonormalize_canonical(virtuals,1.e-6);
	if(virtuals.size()!=spre) std::cout << "removed " << spre-virtuals.size() << " virtuals due to linear dependencies" << std::endl;

	std::vector<std::vector<complex_function_3d> > guess;
	return guess;

}

Tensor<double_complex> Complex_cis::make_CIS_matrix(std::vector<complex_function_3d> virtuals) const {

}



} /* namespace madness */
