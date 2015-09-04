/*
 * tdhf_DFT.cc
 *
 *  Created on: Jun 30, 2014
 *      Author: kottmanj
 */


// DFT functions for tdhf_CIS class
#include<chem/TDA_XC.h>
//#include<madness/mra/mra.h>
using namespace madness;
// Dirac/Slater exchange

// For testing purposes, this is far from beeing numerically stable
// only for GGA
vecfuncT TDA_DFT::evaluate_derivatives()const{
	if(not xcfunctional_.is_gga()) MADNESS_EXCEPTION("get_derivatives function works only for GGA",1);
	vecfuncT derivatives;
	vecfuncT densities;
	densities.push_back(rho_);
	densities.push_back(sigma_);

	for(int what=0;what<3;what++){
		real_function_3d tmp = multiop_values<double,get_derivatives,3>(get_derivatives(xcfunctional_,0,what),densities);
		derivatives.push_back(tmp);
	}
	// GGA does not work
	MADNESS_EXCEPTION("GGA does not work yet, sorry",1);
}

real_function_3d TDA_DFT::convolution_with_kernel(real_function_3d &perturbed_density)const{

	// The LIBXC interface uses the same interface as moldft
	// it is necessary to give the unperturbed density als alpha density (so no factor of 2)
	// the perturbed density has the factor 2 already
	// for GGA we therefore need the factor 2 here when calculating the gradient
	// the sigma_ variable stored here does not need the factor two (will be added by the libxc interface)

	if(xcfunctional_.is_lda()){
	std::vector<real_function_3d> densities;
	densities.push_back(rho_);
	densities.push_back(perturbed_density);
	reconstruct(world, densities);
	refine_to_common_level(world,densities);
	real_function_3d vxc = multiop_values<double,perturbed_vxc,3>(perturbed_vxc(xcfunctional_,0,0),densities);
	return vxc;
	}
	if(xcfunctional_.is_gga()){

		std::vector<real_function_3d> densities;
		densities.push_back(rho_);
		densities.push_back(sigma_);

		// Make the gradient of the perturbed and unperturbed density
		vecfuncT grad_rhoprime;
		std::vector< std::shared_ptr<real_derivative_3d> > gradop;
		gradop =gradient_operator<double,3>(world);
		perturbed_density.reconstruct();
		for (int axis = 0; axis < 3; ++axis){
			grad_rhoprime.push_back((*gradop[axis])(perturbed_density, false)); // delrho
		}
		world.gop.fence(); // NECESSARY

		vecfuncT grad_rho;
		std::vector< std::shared_ptr<real_derivative_3d> > gradop2;
		gradop2 =gradient_operator<double,3>(world);
		rho_.reconstruct();
		for (int axis = 0; axis < 3; ++axis){
			grad_rho.push_back((*gradop2[axis])(rho_, false)); // delrho
		}
		world.gop.fence(); // NECESSARY
		// The unperturbed density is not scaled correctly because it is the spin alhpa density
		// so the gradient has to be scaled with a factor of 2
		scale(world,grad_rho,2.0);

		// Make the perturbed sigma : grad_rho*grad_perturbed_density
		real_function_3d sigma_prime = grad_rho[0]*grad_rhoprime[0] + grad_rho[1]*grad_rhoprime[1] + grad_rho[2]*grad_rhoprime[2];

		densities.push_back(sigma_prime);
		densities.push_back(perturbed_density);
		reconstruct(world,densities);
		refine_to_common_level(world,densities);

		// Make the first part of yanais formula (13)
		// d2fdrho2*rhoprime + 2*d2f/drhodsigma*sigmaprime
		real_function_3d part1 = multiop_values<double,perturbed_vxc,3>(perturbed_vxc(xcfunctional_,0,0),densities);

		// Make the second part of yanais formula (13)
		// -div(grad(rho)*(2*d2fdrhodsigma*rhoprime+4*d2fdsigma2*sigmaprime))

		// first evaluate 2*d2fdrhodsigma*rhoprime+4*d2fdsigma2*sigmaprime
		real_function_3d tmp = multiop_values<double,perturbed_vxc,3>(perturbed_vxc(xcfunctional_,0,1),densities);

		// Then multiply with grad_rho and take the divergence
		vecfuncT vec_tmp = mul(world,tmp,grad_rho);
		real_function_3d part2 = real_factory_3d(world);
		for(int axis=0;axis<3;axis++){
			Derivative<double,3> D = free_space_derivative<double,3>(world, axis);
			part2 -= D(vec_tmp[axis]);
		}

		// Make the last part of yanais formula (13) for closed shell
		// -div(2*dfdsigma*grad_rhoprime)
		// Same game as before
		real_function_3d tmp2 = multiop_values<double,perturbed_vxc,3>(perturbed_vxc(xcfunctional_,0,2),densities);
		tmp2.scale(2.0);
		vecfuncT vec_tmp2 = mul(world,tmp2,grad_rhoprime);
		real_function_3d part3 = real_factory_3d(world);
		for(int axis=0;axis<3;axis++){
			Derivative<double,3> D = free_space_derivative<double,3>(world, axis);
			part3 -= D(vec_tmp2[axis]);
		}

		real_function_3d result = part1+part2+part3;
		result.truncate();

		/// DEBUG
		plot_plane(world,result,"perturbed_potential");
		plot_plane(world,part1,"perturbed_potential_p1");
		plot_plane(world,part2,"perturbed_potential_p2");
		plot_plane(world,perturbed_density,"perturbed_density");
		plot_plane(world,sigma_prime,"sigmaprime");
		plot_plane(world,grad_rhoprime[0],"grad_rhoprime_x");
		plot_plane(world,part3,"perturbed_potential_p3");
		plot_plane(world,grad_rho[0],"grad_rho_x");
		/// DEBUG END

		return result;


	}
	MADNESS_EXCEPTION("Reached end of convolution with kernel function",1);
}

// evaluates fxc*active_mo instead of fxc*perturbed density
vecfuncT TDA_DFT::multiply_with_kernel(vecfuncT &active_mo)const{

	if(xcfunctional_.is_lda()){
		vecfuncT result;
		for(size_t i=0;i<active_mo.size();i++){
			vecfuncT carrier;
			carrier.push_back(rho_);
			carrier.push_back(active_mo[i]);
			reconstruct(world,carrier);
//            active_mo[i].refine_to_common_level(carrier);
			refine_to_common_level(world,carrier);
			real_function_3d tmp = multiop_values<double,perturbed_vxc,3>(perturbed_vxc(xcfunctional_,0,0),carrier);
			result.push_back(tmp);
		}
		plot_plane(world,result.back(),"fxcxlastorbital");
		return result;
	}


	MADNESS_EXCEPTION("Reached end of convolution with kernel 2 function",1);
}

vecfuncT TDA_DFT::apply_kernel(const vecfuncT & x) const{

	if(xcfunctional_.is_gga()) MADNESS_EXCEPTION("apply kernel function not useful with GGA, was for testing purposes only",1);

	// make perturbed density (factor 2 is for closed shell)
	real_function_3d perturbed_density=real_factory_3d(world);
	for(size_t i=0;i<calc.amo.size();i++){
		perturbed_density += calc.amo[i]*x[i];
	}
	perturbed_density.scale(2.0);
	// make a vector which containes as first entry the unperturbed and as last two entries one active orbital and the perturbed density
	// everything in between can be used later for GGA
	vecfuncT applied_kernel;
	for(size_t i=0;i<calc.amo.size();i++){
	std::vector<real_function_3d> densities;
	densities.push_back(rho_);
	if(xcfunctional_.is_gga()) densities.push_back(sigma_);
	densities.push_back(perturbed_density);
	densities.push_back(calc.amo[i]);
	reconstruct(world, densities);
//    perturbed_density.refine_to_common_level(densities);
	refine_to_common_level(world,densities);
	real_function_3d tmp = multiop_values<double,apply_kernel_functor,3>(apply_kernel_functor(xcfunctional_,0,0),densities);
	applied_kernel.push_back(tmp);
	}
	return applied_kernel;
}

real_function_3d TDA_DFT::make_unperturbed_vxc(const real_function_3d &rho)const{

	std::vector<real_function_3d> density;
	density.push_back(rho);

	if(xcfunctional_.is_gga()) density.push_back(sigma_);

	reconstruct(world,density);
//    density[0].refine_to_common_level(density);
	refine_to_common_level(world,density);
	real_function_3d vxc = multiop_values<double,unperturbed_vxc,3>(unperturbed_vxc(xcfunctional_,0,0),density);

	return vxc;

}

real_function_3d TDA_DFT::make_lda_kernel(const real_function_3d &rho)const{
	std::vector<real_function_3d> density;
	density.push_back(rho);
	reconstruct(world,density);
//    density[0].refine_to_common_level(world,density);
	refine_to_common_level(world,density);
	real_function_3d fxc = multiop_values<double,make_fxc,3>(make_fxc(xcfunctional_,0,0),density);

	return fxc;
}



