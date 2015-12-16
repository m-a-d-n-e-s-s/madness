/*
 * tdhf_DFT.h
 *
 *  Created on: Jun 30, 2014
 *      Author: kottmanj
 */

#ifndef TDHF_DFT_H_
#define TDHF_DFT_H_

#include <chem/projector.h>
//#include <examples/mp2.h>

// LIBXC
//#ifdef MADNESS_HAS_LIBXC
//#include <xc.h>
//#endif

#include<examples/nonlinsol.h>
#include<chem/SCF.h>
#include <madness/mra/operator.h>
#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/mra/lbdeux.h>
#include <madness/misc/ran.h>

#include <cmath>       /* pow */



namespace madness{

// Stolen from xc_functional helper structure of lib/chem/xcfunctional.h
struct unperturbed_vxc{
	unperturbed_vxc(const XCfunctional& xc, int ispin,int what)
    : xc(&xc), what(what), ispin(ispin){}

    const XCfunctional* xc;
    const int what;
    const int ispin;

    madness::Tensor<double> operator()(const madness::Key<3> & key, const std::vector< madness::Tensor<double> >& t) const
    {
        MADNESS_ASSERT(xc);
        madness::Tensor<double> r; // debug = xc->vxc(t, ispin, what);
        return r;
    }
};

struct make_fxc{
	make_fxc(const XCfunctional& xc, int ispin,int what)
    : xc(&xc), what(what), ispin(ispin){}

    const XCfunctional* xc;
    const int what;
    const int ispin;

    madness::Tensor<double> operator()(const madness::Key<3> & key, const std::vector< madness::Tensor<double> >& t) const
    {
        MADNESS_ASSERT(xc);
        madness::Tensor<double> r; // debug = xc->fxc(t, ispin, what);
        return r;
    }
};

struct get_derivatives{
	get_derivatives(const XCfunctional& xc, int ispin, int what)
	: xc(&xc), what(what), ispin(ispin){}

	    const XCfunctional* xc;
	    const int what;
	    const int ispin;

	    madness::Tensor<double> operator()(const madness::Key<3> & key, const std::vector< madness::Tensor<double> >& t) const
	        {
	    	MADNESS_ASSERT(xc);
	    	madness::Tensor<double> result;

// debug
//	    	if(what == 0){
//	    		// gives back d2fdrho2
//	    		result=xc->fxc(t,ispin,0);
//	    	}
//	    	else if(what == 1){
//	    		// gives back d2fdrhodsigma
//	    		result=xc->fxc(t,ispin,1);
//	    	}
//	    	else if(what == 2){
//	    		// gives back d2fdsigma2
//	    		result=xc->fxc(t,ispin,2);
//	    	}
//	    	else MADNESS_EXCEPTION("what can only be in the range of 0-2",1);

	    	return result;
	        }
};

struct perturbed_vxc{
	perturbed_vxc(const XCfunctional& xc, int ispin,int what)
    : xc(&xc), what(what), ispin(ispin){}

    const XCfunctional* xc;
    const int what;
    const int ispin;

    madness::Tensor<double> operator()(const madness::Key<3> & key, const std::vector< madness::Tensor<double> >& t) const
    {
        MADNESS_ASSERT(xc);
        // Make fxc kernel
        madness::Tensor<double> result;
        if(xc->is_lda()){
        	std::vector< madness::Tensor<double> > rho;
        	rho.push_back(t[0]);
        	madness::Tensor<double> fxc; // debug = xc->fxc(rho, ispin, 0);



        	// multiply kernel with perturbed density (last element of t-vector)
        	fxc.emul(t.back());

        	result = fxc;

        }
        // For GGA the uncomming vector of tensors should contain:
        // t[0] = unperturbed density
        // t[1] = unperturbed sigma (gradient of rho * gradient of rho)
        // t[2] = perturbed sigma (gradient of rho times gradient of rhoprime)
        // t[3] = perturbed density (rhoprime)

         if(xc->is_gga()){
        	std::vector< madness::Tensor<double> > rho_and_sigma;

        	madness::Tensor<double> d2fdrho2;
        	madness::Tensor<double> d2fdrhodsigma;
        	madness::Tensor<double> d2fdsigma2;
        	rho_and_sigma.push_back(t[0]);
        	rho_and_sigma.push_back(t[1]);

// debug
//        	if(what==0){
//        		// Create the first term of yanais formula (13) for closed shells
//        		// d2f/drho2*rhoprime + 2*d2f/drhodsigma*(grad_rho*grad_rhoprime)
//        		d2fdrho2 = xc->fxc(rho_and_sigma,ispin,0);
//        		d2fdrhodsigma = xc->fxc(rho_and_sigma,ispin,1);
//
//        		result = d2fdrho2;
//        		result.emul(t[3]);
//        		d2fdrhodsigma.emul(t[2]);
//        		result.gaxpy(1.0,d2fdrhodsigma,2.0);
//
//        	}
//        	else if(what==1){
//        		// Create part of the second term of yanais formula (13) for closed shell
//        		// 2*d2fdrhodsigma * rhoprime + 4*d2fdisgma2(grad_rho*grad_rhoprime);
//        		// This has to be multiplied with grad_rho and the divergence should be taken afterwards
//        		d2fdrhodsigma = xc->fxc(rho_and_sigma,ispin,1);
//        		d2fdsigma2 = xc->fxc(rho_and_sigma,ispin,2);
//
//        		result = d2fdrhodsigma;
//        		result.emul(t[3]);
//        		d2fdsigma2.emul(t[2]);
//        		result.gaxpy(2.0,d2fdsigma2,4.0);
//
//        	}
//        	else if(what==2){
//        		// Create the last part of yanais formula (13) for closed shell
//        		// df/dsigma
//        		// This has to be multiplied with 2, contracted with grad_rhoprime and the divergence taken afterwards
//        		result = xc->vxc(rho_and_sigma,ispin,1);
//
//        	} else {
//        	    MADNESS_EXCEPTION("What parameter of convolute_with_kernel was not from 0-2",1);
//        	}



        }
        return result;
    }
};

struct apply_kernel_functor{
	apply_kernel_functor(const XCfunctional & xc, int ispin, int what): xc(&xc),what(what),ispin(ispin){}
	const XCfunctional * xc;
	const int what;
	const int ispin;

	madness::Tensor<double> operator()(const madness::Key<3> & key, const std::vector<madness::Tensor<double> >& t)const
	{
		MADNESS_ASSERT(xc);
		// Make fxc kernel with rho_ (first entry of t)
		std::vector<madness::Tensor<double> > rho;
		rho.push_back(t[0]);
		if(xc->is_gga()) rho.push_back(t[1]);
		madness::Tensor<double> fxc; // debug = xc->fxc(rho, ispin, what);
		// multiply the kernel with the density and the active mo (product is the last entry of t)
		fxc.emul(t[1]);
		fxc.emul(t[2]);
		return fxc;
	}
};



class TDA_DFT {

public:
	static double munge_density(double rho){
        if (fabs(rho) <= 1.e-8) rho=1.e-8;
        return rho;
	}

	TDA_DFT(World &world,const SCF &calc): world(world),calc(calc), xcfunctional_(calc.xc){

		// Potential check
		std::cout << "Used Potential:\n";
		std::cout << "lda: " << xcfunctional_.is_lda();
		std::cout << "gga: " << xcfunctional_.is_gga();
		std::cout << " spin polarized: " << xcfunctional_.is_spin_polarized() << std::endl;

		// Make unperturbed Density and unperturbed xc-potential
		// The density for spin unrestricted and closed shell is arho (without factor 2, because this is added in the XCfunctional class automatically)
		real_function_3d rho = real_factory_3d(world);
		for(size_t i=0;i<calc.amo.size();i++){
			rho +=calc.amo[i]*calc.amo[i];
		}

		plot_plane(world,rho,"rho");
		rho_=rho;
		rho.unaryop(munge_density);


//		// Make the contracted gradient of the unperturbed density sum_i (d/dxi rho)^2
//		if(xcfunctional_.is_gga()){
//			vecfuncT delrho;
//			std::vector< std::shared_ptr<real_derivative_3d> > gradop;
//			gradop =gradient_operator<double,3>(world);
//			rho.reconstruct();
//			for (int axis = 0; axis < 3; ++axis){
//				delrho.push_back((*gradop[axis])(rho, false)); // delrho
//			}
//			world.gop.fence(); // NECESSARY
//
//			sigma_=	delrho[0] * delrho[0] + delrho[1] * delrho[1]+ delrho[2] * delrho[2];     // sigma_aa
//		}

	}
	///Disable copy constructor
	TDA_DFT(const TDA_DFT &other);


	/// The world
	World &world;

	void print_information();
	vecfuncT evaluate_derivatives()const;
	real_function_3d get_perturbed_potential(const real_function_3d &perturbed_density){return V(perturbed_density);}
	real_function_3d convolution_with_kernel(real_function_3d &perturbed_density)const;
	vecfuncT apply_kernel(const vecfuncT &x)const;
	real_function_3d get_unperturbed_vxc()const{
		if(vxc_.is_initialized())
		return vxc_;
		else{
			real_function_3d vxc=make_unperturbed_vxc(rho_);
			vxc_ = vxc;
			return vxc_;
		}
	}
	XCfunctional get_xcfunctional()const{return xcfunctional_;}
	real_function_3d calculate_unperturbed_vxc(const real_function_3d &rho)const{return make_unperturbed_vxc(rho);}
	vecfuncT get_lda_intermediate(vecfuncT & mos)const{return multiply_with_kernel(mos);}
	real_function_3d get_fxc()const{
		if(fxc_.is_initialized())
		return fxc_;
		else{
				real_function_3d fxc = make_lda_kernel(rho_);
				fxc.truncate();
				fxc_=fxc;
				return fxc_;
			}
		}
	// Test which type of functional is used
	bool is_gga()const{
		return xcfunctional_.is_gga();
	}

private:
	/// The SCF calculation as reference (make density and assign xc potential)
	/// Try to use the XCfunctional already implemented and used by moldft (now SCF)
	const SCF &calc;

	/// The XCfunctional class
	const XCfunctional &xcfunctional_;
	/// The unperturbed density
	real_function_3d rho_;
	/// The contracted gradient of the unperturbed density
	real_function_3d sigma_;
	/// unperturbed exchange correlation potential
	mutable real_function_3d vxc_;
	/// The LDA kernel
	mutable real_function_3d fxc_;
	/// Type of exchange and correlation
	real_function_3d vx_dirac(const real_function_3d &perturbed_density);
	real_function_3d V(const real_function_3d &perturbed_density);
	real_function_3d make_unperturbed_vxc(const real_function_3d &rho)const;
	real_function_3d make_lda_kernel(const real_function_3d &rho)const;
	vecfuncT multiply_with_kernel(vecfuncT &active_mo)const;



};

} // madness namesapce
#endif /* TDHF_DFT_H_ */
