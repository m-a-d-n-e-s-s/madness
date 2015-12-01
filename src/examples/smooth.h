/*
 * smooth.h
 *
 *  Created on: Nov 25, 2015
 *      Author: kottmanj
 */

#ifndef SRC_EXAMPLES_SMOOTH_H_
#define SRC_EXAMPLES_SMOOTH_H_

#include <iostream>
#include <string>
#include <chem/nemo.h>
#include <madness/mra/funcplot.h>
#include <cmath>
#include <chem/SCFOperators.h>

double estimate_area(double x) {
	double thresh =FunctionDefaults<3>::get_thresh();
	if(fabs(x)<thresh) return 0.0;
	else return 1.0;
    }

static double test_1d_functor(const coord_1d &x){
	return exp(-fabs(x[0]))+0.1*sin(30.0*x[0]);
}

static double xf(const coord_3d &x){
	return x[1];
}
static double r2(const coord_3d &x){
	return x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
}

static double slater_functor(const coord_3d &x){
	double rsq = r2(x);
	double r = sqrt(rsq);
	return exp(-r);
}

static double mask_factor = 0.75;
static double mask_functor(const coord_3d &x){
	double r = sqrt(r2(x));
	return 1-erf(mask_factor*r);
}
static double inv_mask_functor(const coord_3d &x){
	double r = sqrt(r2(x));
	return erf(mask_factor*r);
}

struct merging_operator {
	merging_operator() : thresh_(FunctionDefaults<3>::get_thresh()) {}
	merging_operator(const double thresh) : thresh_(thresh) {}
    void operator()(const Key<3>& key,
                    real_tensor U,
                    const real_tensor& function,
                    const real_tensor& smoothed_function) const {
        ITERATOR(U,
                 double f = function(IND);
                 double sf = smoothed_function(IND);
                 if (fabs(f)<thresh_)
                     U(IND) = sf;
                 else
                     U(IND) = f;
                 );
    }

    template <typename Archive>
    void serialize(Archive& ar) {}
    double thresh_;
};

struct munging_operator {
	munging_operator() : thresh_(FunctionDefaults<3>::get_thresh()) {}
	munging_operator(const double thresh) : thresh_(thresh) {}
    void operator()(const Key<3>& key,
                    real_tensor U,
                    const real_tensor& function,
                    const real_tensor& smoothed_function) const {
        ITERATOR(U,
                 double f = function(IND);
                 double sf = smoothed_function(IND);
                 if (fabs(f)<thresh_)
                     U(IND) = 0.0;
                 else
                     U(IND) = f;
                 );
    }

    template <typename Archive>
    void serialize(Archive& ar) {}
    double thresh_;
};

struct asymptotic_density{
public:
	asymptotic_density(const double &ie):ionization_energy_(ie), center_of_charge_(0.0){}
	asymptotic_density(const double &ie, const coord_3d &coc):ionization_energy_(ie), center_of_charge_(coc){}
		double operator ()(const coord_3d &x1)const{
		coord_3d x;
		x[0] = x1[0]-center_of_charge_[0];
		x[1] = x1[1]-center_of_charge_[1];
		x[2] = x1[2]-center_of_charge_[2];
		double r = sqrt((x[0]*x[0]+x[1]*x[1]+x[2]*x[2]));
		double result= exp(-2.0*sqrt(2.0*fabs(ionization_energy_))*r);
		return result;
	}

	const double ionization_energy_;
	const coord_3d center_of_charge_;
};

struct asymptotic_slater_kernel{
public:
	asymptotic_slater_kernel(const asymptotic_density &tmp): asymptotic_rho_(tmp){}

	double operator()(const coord_3d &x)const{
		double prefactor = 1.0/M_PI;
		double r = sqrt((x[0]*x[0]+x[1]*x[1]+x[2]*x[2]));
		double kernel = prefactor*exp(2.0/3.0*2.0*sqrt(2.0*fabs(asymptotic_rho_.ionization_energy_))*r);//pow(rho,2.0/3.0);
		double rho = asymptotic_rho_(x);
 		return kernel*rho;
	}

private:
	const asymptotic_density asymptotic_rho_;
};

struct slater_kernel {
	slater_kernel() : thresh_(FunctionDefaults<3>::get_thresh()) {}
	slater_kernel(const double thresh) : thresh_(thresh) {}
    void operator()(const Key<3>& key,
                    real_tensor U,
                    const real_tensor& rho,
                    const real_tensor& pt_rho) const {
    	double safety = thresh_*0.01;
    	double safety_kernel = pow(safety,-2.0/3.0);
    	double safety_kernel_safety_pd = pre*pow(safety,1.0/3.0);
        ITERATOR(U,
                 double d = rho(IND);
                 double pd = pt_rho(IND);
                 if(d<safety and fabs(pd)<safety) U(IND)=safety_kernel_safety_pd;
                 else if(d<safety)  U(IND)= pre*safety_kernel*pd;
                 else U(IND) = pre*pow(d,-2.0/3.0)*pd;
                 );
    }

    template <typename Archive>
    void serialize(Archive& ar) {}
    double thresh_;
    double pre = 1.0/M_PI;
};

static double analytical_slater_functor(double rho){
	double tmp;
	if(rho < 0.0) tmp = 0.0;
	else tmp = pow(rho,1.0/3.0);
	double prefactor = 1.0/M_PI;
	return prefactor*tmp;
}



class smooth{
public:
	smooth(World &world,const Nemo &nemo): world(world), asymptotic_density_functor(make_asymtotic_density_functor(nemo)) {
		const vecfuncT nemos = nemo.get_calc()->amo;
		ionization_energy_ = nemo.get_calc()->aeps(nemos.size()-1);
		real_function_3d nemo_rho = real_factory_3d(world);
		R_ = nemo.nuclear_correlation -> function();
		R2_= nemo.nuclear_correlation -> square();
		for(auto nemo:nemos){
			nemo_rho = nemo*nemo;
		}
		nemo_rho.truncate();
		nemo_rho_ = copy(nemo_rho);
		rho_ = R2_*nemo_rho;
		plot_plane(world,rho_,"density");
		plot_plane(world,nemo_rho_,"nemo_density");
		mask_ = real_factory_3d(world).f(mask_functor);
		inv_mask_ = real_factory_3d(world).f(inv_mask_functor);
		make_plots(mask_,"mask");
		make_plots(inv_mask_,"inv_mask");
		//test
		real_function_3d dm = rho_*mask_;
		real_function_3d dm2= rho_*inv_mask_;
		real_function_3d d = dm+dm2;
		double diff = (d-rho_).norm2();
		std::cout << "Mask test: diff =" << diff << std::endl;
	}

	void estimate_diffuseness(const real_function_3d &f, const std::string & msg="function")const{
		real_function_3d range_of_f = copy(f);
		range_of_f.unaryop(estimate_area);
		dirac_smoothing(range_of_f,0.04,msg);
		make_plots(range_of_f,msg);
		double box_width = FunctionDefaults<3>::get_cell_min_width();
		double box_volume = box_width*3.0;
		if(world.rank()==0) std::cout << "range_of_" << msg << " = " << range_of_f.norm2() << std::endl;
		if(world.rank()==0) std::cout << "range_of_" << msg << " compared to full box = " << range_of_f.norm2()/box_volume << std::endl;
	}

	real_function_3d dirac_smoothing(const real_function_3d &f, const double &eps=0.04,const std::string &name="function")const{
		output("\nSmoothing " + name + " with eps " + stringify(eps));
		double exponent = 1.0/(2.0*eps);
		Tensor<double> coeffs(1), exponents(1);
	    exponents(0L) =  exponent;
	    coeffs(0L) = gauss_norm(exponent);
	    SeparatedConvolution<double,3> op(world, coeffs, exponents);
	    real_function_3d smoothed_f = apply(op,f);
	    double diff = (f-smoothed_f).norm2();
	    output("||" + name + " - smoothed_" + name + "||="+stringify(diff));
	    if(diff>10.0*FunctionDefaults<3>::get_thresh()) output("Smoothing Difference far above the accuracy threshold\n");
	    else if(diff>FunctionDefaults<3>::get_thresh()) output("Smoothing Difference above the accuracy threshold\n");
	    else output("Smoothing did not affect accuracy\n");
	    make_plots(f,name);
	    smoothed_f = binary_op(smoothed_f, smoothed_f, munging_operator(FunctionDefaults<3>::get_thresh()));
	    make_plots(smoothed_f,"smoothed_"+name);
	    real_function_3d merged = merge_functions(f,smoothed_f,name);
	    make_plots(merged,"merged_"+name);
	    return merged;

	}

	real_function_3d merge_functions(const real_function_3d &f, const real_function_3d &sf,const std::string &name="f")const{
//		real_function_3d result =  binary_op(f, sf, merging_operator(FunctionDefaults<3>::get_thresh()));
//		output("||f-munged_f||="+stringify(diff));

		real_function_3d part1 = mask_*f;
		real_function_3d part2 = inv_mask_*sf;
		real_function_3d result = part1 + part2;

		make_plots(part1,"masked_"+name);
		make_plots(part2,"masked_s"+name);

		double diff = (f-result).norm2();
		if(diff>FunctionDefaults<3>::get_thresh()) output("Munging Error is above the accuracy threshold!\n");

		return result;
	}

	void make_plots(const real_function_3d &f,const std::string &name="function")const{
		double width = FunctionDefaults<3>::get_cell_min_width()/2.0 - 1.e-3;
		plot_plane(world,f,name);
		coord_3d start(0.0); start[0]=-width;
		coord_3d end(0.0); end[0]=width;
		plot_line(("line_"+name).c_str(),1000,start,end,f);
	}

	void make_smooth_gradient()const{
		double eps =0.1;
		std::vector < std::shared_ptr<real_derivative_3d> > gradop;
		gradop = gradient_operator<double, 3>(world);
		vecfuncT result;
		// make gradient with density
		real_function_3d smooth_density = dirac_smoothing(rho_,eps,"density");
		real_function_3d dx_sd = (*gradop[0])(smooth_density);
		make_plots(dx_sd,"dx_smoothed_density");
		real_function_3d dx_d = (*gradop[0])(rho_);
		make_plots(dx_d,"dx_density");
		real_function_3d s_dx_sd = dirac_smoothing(dx_sd,eps,"s_dx_density");
		make_plots(s_dx_sd,"smoothed_dx_smoothed_density");
		real_function_3d s_dx_d = dirac_smoothing(dx_d,eps,"smoothed_dx_density");

		// make gradient with nemo_density:  dx(density) = dx(R2*nemo_density) = dx(R2)*nemo_density + R2*dx(nemo_density)
	}

	void make_smooth_slater_kernel()const{
		real_function_3d r2f = real_factory_3d(world).f(r2);
		real_function_3d pt_rho = r2f*rho_;
		asymptotic_slater_kernel asymptotic_slater_functor(asymptotic_density_functor);
		real_function_3d asymptotic_slater_kernel = real_factory_3d(world).functor(asymptotic_slater_functor);
		asymptotic_slater_kernel = asymptotic_slater_kernel*r2f;
		make_plots(asymptotic_slater_kernel,"asymptotic_slater_kernel");

		real_function_3d analytical_slater_kernel = copy(rho_);
		analytical_slater_kernel.unaryop(analytical_slater_functor);
		analytical_slater_kernel = analytical_slater_kernel*r2f;
		make_plots(analytical_slater_kernel,"analytical_slater_kernel");

		real_function_3d numerical_slater_kernel = binary_op(rho_, pt_rho, slater_kernel(FunctionDefaults<3>::get_thresh()));
		make_plots(numerical_slater_kernel,"numerical_slater_kernel");
	}

	void make_smooth_XC_kernel(const std::string &xc_data)const{
		double homo = ionization_energy_;
		asymptotic_density tmp(homo);
		real_function_3d asymptotic_rho = real_factory_3d(world).functor(tmp);
		make_plots(rho_,"rho");
		make_plots(asymptotic_rho,"asymptotic_rho");

		real_function_3d ana_rho = real_factory_3d(world).f(slater_functor);
		double norm = ana_rho.norm2();
		ana_rho.scale(1.0/norm);
		ana_rho = (ana_rho*ana_rho);

		estimate_diffuseness(ana_rho,"range_of_ana_rho");
		estimate_diffuseness(rho_,"range_of_rho");

		real_function_3d x = real_factory_3d(world).f(r2);
		real_function_3d pt_rho = x*ana_rho;
		// make XCOperator with std density and smoothed_density
		double eps = 0.04;
		real_function_3d smooth_rho = ana_rho;// dirac_smoothing(ana_rho,eps,"rho");
		pt_rho = dirac_smoothing(pt_rho,eps,"pt_rho");
		XCOperator xcop(world,xc_data,false,smooth_rho,smooth_rho);
		real_function_3d xc_ptrho = xcop.apply_xc_kernel(pt_rho);
		real_function_3d s_xc_ptrho = dirac_smoothing(xc_ptrho,eps,"xc_ptrho");
	}

	void test()const{
		output("Smoothing Nemo_Density");
		double eps = 0.1;
		for(size_t i=0;i<10;i++){
		dirac_smoothing(nemo_rho_,eps,"nemo_density_"+stringify(eps));
		eps -= 0.01;
		}
		output("Smoothing Density");
		eps = 0.1;
		for(size_t i=0;i<10;i++){
		dirac_smoothing(rho_,eps,"density_"+stringify(eps));
		eps -= 0.01;
		}
	}

	void test_1d()const{
		FunctionDefaults<1>::set_thresh(1.e-5);
		FunctionDefaults<1>::set_k(8);
		FunctionDefaults<1>::set_cubic_cell(-25.0,25.0);
		real_function_1d f = real_factory_1d(world).f(test_1d_functor);
		Tensor<double> coeffs(1), exponents(1);
	    exponents(0L) =  1.0/(2.0*0.04);
	    coeffs(0L) = pow((1.0/(2.0*0.04))/M_PI,0.5);
	    SeparatedConvolution<double,1> op(world, coeffs, exponents);
	    real_function_1d sf = apply(op,f);
	    double diff = (f-sf).norm2();
	    output("||f - smoothed_f||_1D="+stringify(diff));
	    coord_1d start; start[0]=-10.0;
	    coord_1d end; end[0]=10.0;
	    plot_line("smoothed_1d",1000,start,end,sf);
	    plot_line("noise_1d",1000,start,end,f);
	}

	double gauss_norm(const double &exponent)const{
		return pow(exponent/M_PI,0.5*3.0);
	}

	void output(const std::string &msg)const{
		if(world.rank()==0) std::cout << msg << std::endl;
	}

private:
	World &world;
	real_function_3d rho_;
	real_function_3d nemo_rho_;
	real_function_3d R_;
	real_function_3d R2_;
	real_function_3d mask_;
	real_function_3d inv_mask_;
	double ionization_energy_;
	const asymptotic_density asymptotic_density_functor;
	asymptotic_density make_asymtotic_density_functor(const Nemo &nemo)const{
		size_t nocc = nemo.get_calc()->amo.size();
		double ehomo = nemo.get_calc()->aeps(nocc-1);
		const std::vector<Atom>& atoms = nemo.get_calc()->molecule.get_atoms();
		double total_nuclear_charge = 0.0;
		coord_3d mx(0.0);
		for(auto x:atoms){
			total_nuclear_charge += x.q;
			mx[0] += x.q*x.x;
			mx[1] += x.q*x.y;
			mx[2] += x.q*x.z;
		}
		mx *= 1.0/total_nuclear_charge;
		std::cout << "Center of Charge is\n" << mx[0] << ", " << mx[1] << ", " << mx[2] << std::endl;
		asymptotic_density tmp(ehomo,mx);
		return tmp;
	}
};



#endif /* SRC_EXAMPLES_SMOOTH_H_ */
