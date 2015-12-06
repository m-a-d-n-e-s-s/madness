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
#include <madness/mra/function_common_data.h>

static double RHOMIN = 1.e-15;

double make_log(double x){
	if(x<0.0) return log(RHOMIN);
	return log(x);
}

double make_exp(double x){
	return exp(x);
}

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

static double mask_factor = 5.0;
static double cutoff_radius = 5.0;
//static double mask_functor(const coord_3d &x){
//	double r = sqrt(r2(x));
//	//return 1-erf(mask_factor*r);
//	return 0.5*(1.0 - tanh(mask_factor*(r-cutoff_radius)));
//}
static double mask_functor_box(const coord_3d &x){
	double r = sqrt(r2(x));
	//return 1-erf(mask_factor*r);
	return 0.5*(1.0 - tanh(1.0*(r-15)));
}
//static double inv_mask_functor(const coord_3d &x){
//	double r = sqrt(r2(x));
//	return 0.5*(1.0 + tanh(mask_factor*(r-cutoff_radius)));
//}

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

struct asymptotic_slater : public FunctionFunctorInterface<double,3>{
public:
	asymptotic_slater(const double &ipot_): ipot(ipot_){}
	const double ipot;
	const double prefactor = 1.0/M_PI;

	double operator()(const coord_3d &x)const{
		const double r = sqrt((x[0]*x[0]+x[1]*x[1]+x[2]*x[2]));
		const double kernel = prefactor*exp(2.0/3.0*2.0*sqrt(2.0*fabs(ipot))*r);
		return kernel;
	}
};

struct apply_kernel_helper{

public:
	apply_kernel_helper(const double &ipot_): ipot(ipot_), asl(ipot_){}

	struct slater_kernel{
		double operator()(const double & rho)const{
			double tmp = rho;
			if(tmp<0.0) tmp=1.e-10;
			return prefactor*pow(tmp,-2.0/3.0);
		}
		const double prefactor = 1.0/M_PI;
	};



	madness::Tensor<double> slater_apply(const std::vector< madness::Tensor<double> >& t,const madness::Key<3> & key,const FunctionCommonData<double,3> cdata) const{
	madness::Tensor<double> rho = t.front();
	madness::Tensor<double> pt_rho = t.back();
    const long & k = FunctionDefaults<3>::get_k();
	const Tensor<double>& quad_x = cdata.quad_x;

	// if all points of the perturbed density are below the threshold: Set the result to zero
	bool pt_zero = true;
	if(pt_rho.absmax()>safety_thresh_ptrho) pt_zero = false;

	// fast return if pt_rho is zero
	if(pt_zero){
		std::vector<long> kdims(3,k);
		Tensor<double> zero_result(kdims,true);
		return zero_result;
	}

	// if all points of the unperturbed density are below the threshold: Use asymptotical density
	bool rho_ana = true;
	//if(rho.absmax()>safety_thresh_rho) rho_ana = false;

	// analytical kernel based on asymptotical density
	Tensor<double> kernel_val(k,k,k);
	if(rho_ana){
		madness::fcube<double,3>(key, asl, quad_x,kernel_val);
	}else{
		kernel_val = copy(rho);
		slater_kernel slater_tmp;
		kernel_val.unaryop<slater_kernel>(slater_tmp);
	}

	// multiply kernel values with pt_rho values
	Tensor<double> result = kernel_val.emul(pt_rho);
	return result;

	}

	const double prefactor = 1.0/M_PI;
	const double ipot;
	asymptotic_slater asl;
	const double safety_thresh_rho = FunctionDefaults<3>::get_thresh()*0.01;
	const double safety_thresh_ptrho = FunctionDefaults<3>::get_thresh()*0.01;
};

struct slater_kernel_apply {

    const FunctionCommonData<double,3>& cdata;

    const apply_kernel_helper *xc;
    slater_kernel_apply(const apply_kernel_helper &xc_) : cdata(FunctionCommonData<double,3>::get(FunctionDefaults<3>::get_k())), xc(&xc_) {}

    madness::Tensor<double> operator()(const madness::Key<3> & key,
            const std::vector< madness::Tensor<double> >& t) const {
        madness::Tensor<double> r = xc->slater_apply(t,key,cdata);
        return r;
    }
};


template<std::size_t NDIM>
static double make_radius(const Vector<double,NDIM> &x){
	double r2=0.0;
	for(double xi:x) r2+=xi*xi;
	return sqrt(r2);
}
template<std::size_t NDIM>
static double mask_munging(const double x){
	if(fabs(x)<0.1*FunctionDefaults<NDIM>::get_thresh()) return 0.0;
	else return x;
}
template<std::size_t NDIM>
static double unitfunctor(const Vector<double,NDIM> &x){
	return 1.0;
}

template<typename T, std::size_t NDIM>
class smooth{
public:
	smooth(World&world):
		world(world),
		box_mask_(make_mask(1.0,0.75*FunctionDefaults<NDIM>::get_cell_min_width())),
		mask_(make_mask(1.0,0.5*FunctionDefaults<NDIM>::get_cell_min_width())),
		inv_mask_(make_inv_mask(1.0,0.5*FunctionDefaults<NDIM>::get_cell_min_width())){}
	smooth(World&world, const double &box_mask_factor, const double box_mask_cutoff):
		world(world),
		box_mask_(make_mask(box_mask_factor,box_mask_cutoff)),
		mask_(make_mask(1.0,0.5*FunctionDefaults<NDIM>::get_cell_min_width())),
		inv_mask_(make_inv_mask(1.0,0.5*FunctionDefaults<NDIM>::get_cell_min_width())){}

	struct mask_functor{
	public:
		mask_functor(const double &factor, const double &cutoff): factor_(factor), cutoff_radius_(cutoff), shift_(0.0){}
		mask_functor(const double &factor, const double &cutoff, const Vector<double,NDIM> shift): factor_(factor), cutoff_radius_(cutoff), shift_(shift){std::cout << "shift is " << shift_ << std::endl;}
		double operator()(const Vector<double,NDIM> &x)const{
			const Vector<double,3> shifted_x = x-shift_;
			const double r = make_radius<NDIM>(shifted_x);
			return (1.0-tanh(r));
		}
	private:
		const double factor_;
		const double cutoff_radius_;
		const Vector<double,NDIM> shift_;
	};
	struct inv_mask_functor{
	public:
		inv_mask_functor(const double &factor, const double &cutoff): factor_(factor), cutoff_radius_(cutoff), shift_(0.0){}
		inv_mask_functor(const double &factor, const double &cutoff, const Vector<double,NDIM> shift): factor_(factor), cutoff_radius_(cutoff), shift_(shift){}
		double operator()(const Vector<double,NDIM> &x)const{
			const Vector<double,3> shifted_x = x-shift_;
			const double r = make_radius<NDIM>(shifted_x);
			return 0.5*(1.0 + tanh(factor_*(r-cutoff_radius_)));
		}
	private:
		const double factor_;
		const double cutoff_radius_;
		const Vector<double,NDIM> shift_;
	};

	Function<T,NDIM> make_mask(const double &factor, const double &cutoff)const{
		Vector<double,NDIM> coord(0.0);
		return make_mask(factor,cutoff,coord);
	}
	Function<T,NDIM> make_inv_mask(const double &factor, const double &cutoff)const{
		Vector<double,NDIM> coord(0.0);
		return make_inv_mask(factor,cutoff,coord);
	}
	Function<T,NDIM> make_mask(const double &factor, const double &cutoff,const Vector<double,3> &shift)const{
		mask_functor mask_functor_object(factor,cutoff,shift);
		Function<T,NDIM> mask = FunctionFactory<T,NDIM>(world).functor(mask_functor_object);
		// set values below the threshold to zero
		mask.unaryop(mask_munging<NDIM>);
		return mask;
	}

	Function<T,NDIM> make_inv_mask(const double &factor, const double &cutoff, const Vector<double,3> &shift)const{
		//inv_mask_functor inv_mask_functor_object(factor,cutoff,shift);
		//Function<T,NDIM> inv_mask = FunctionFactory<T,NDIM>(world).functor(inv_mask_functor_object);
		//return inv_mask;
		Function<T,NDIM> unit = FunctionFactory<T,NDIM>(world).f(unitfunctor);
		return unit-mask_;
	}

	void set_simple_mask(const double &factor, const double &cutoff){
		mask_ = make_mask(factor,cutoff);
		inv_mask_ = make_inv_mask(factor,cutoff);
	}

	// place a simple mask on every atom center
	void set_molecule_mask(const double &factor, const double &cutoff, const std::vector<Atom> &atoms){
		if(world.rank()==0) std::cout << "Making Mask for " << atoms.size() << "-atomic molecule with exponent=" << factor << " and cutoff radius =" << cutoff << "\n";
		if(NDIM != 3) MADNESS_EXCEPTION("set_molecule_mask of smooth class only works for 3D functions",1);
		Function<T,NDIM> mask = FunctionFactory<T,NDIM>(world);
		Function<T,NDIM> inv_mask = FunctionFactory<T,NDIM>(world);
		for(Atom a:atoms){
			const Vector<double,3> shift = a.get_coords();
			std::cout << "ashift is " << shift << std::endl;
			Function<T,NDIM> simple_mask_on_a = make_mask(factor,cutoff,shift);
			make_plots(simple_mask_on_a,"mask_on_atom");
			mask += simple_mask_on_a;
//			Function<T,NDIM> simple_inv_mask_on_a = make_inv_mask(factor,cutoff);
//			inv_mask += simple_inv_mask_on_a;
		}
		make_plots(mask,"molecule_mask");
		make_plots(inv_mask,"molecule_inv_mask");
		mask_ = mask;
		inv_mask_ = make_inv_mask(0.0,0.0);
	}

	Function<T,NDIM> make_explog(const Function<T,NDIM> &f, const std::string & msg="function")const{
		// make logf
		Function<T,NDIM> logf = copy(f);
		logf.unaryop(make_log);

		// smooth logf
		Function<T,NDIM> smoothed_logf = gaussian_smoothing(logf,0.5);

		// reconstruct f via f=exp^(log(f));
		Function<T,NDIM> explogf = copy(smoothed_logf);
		explogf.unaryop(make_exp);

		// apply the box_mask to eliminate weird boundary effects
		explogf = box_mask_*explogf;

		// make some plots
		make_plots(f,msg);
		make_plots(smoothed_logf,"smoothed_log"+msg);
		make_plots(logf,"log"+msg);
		make_plots(explogf,"explog_"+msg);
		return explogf;
	}

	void estimate_diffuseness(const Function<T,NDIM> &f, const std::string & msg="function")const{
		Function<T,NDIM> range_of_f = copy(f);
		range_of_f.unaryop(estimate_area);
		gaussian_smoothing(range_of_f,0.04,msg);
		make_plots(range_of_f,msg);
		double box_width = FunctionDefaults<3>::get_cell_min_width();
		double box_volume = box_width*3.0;
		if(world.rank()==0) std::cout << "range_of_" << msg << " = " << range_of_f.norm2() << std::endl;
		if(world.rank()==0) std::cout << "range_of_" << msg << " compared to full box = " << range_of_f.norm2()/box_volume << std::endl;
	}

	Function<T,NDIM> gaussian_smoothing(const Function<T,NDIM> &f, const double &eps=0.04,const std::string &name="nooutput")const{
		if(name!="noputput") output("\nSmoothing " + name + " with eps " + stringify(eps));
		double exponent = 1.0/(2.0*eps);
		Tensor<double> coeffs(1), exponents(1);
	    exponents(0L) =  exponent;
	    coeffs(0L) = gauss_norm(exponent);
	    SeparatedConvolution<T,NDIM> op(world, coeffs, exponents);
	    Function<T,NDIM> smoothed_f = apply(op,f);
	    make_plots(smoothed_f,"smoothed_"+name);
	    Function<T,NDIM> merged = merge_functions(f,smoothed_f,name);
	    make_plots(merged,"smoothed_and_merged_"+name);

	    double diff_smoothed = (f-smoothed_f).norm2();
	    double diff_merged = (f-merged).norm2();
	    if(world.rank()==0 and name!="noputput"){
	    	std::cout << "||f-smoothed_f||=" << diff_smoothed << std::endl;
	    	std::cout << "||f-merged_f||  =" << diff_merged << std::endl;
	    }
	    return merged;

	}

	Function<T,NDIM> merge_functions(const Function<T,NDIM> &f, const Function<T,NDIM> &sf,const std::string &name="f")const{
		Function<T,NDIM> part1 = mask_*f;
		Function<T,NDIM> part2 = inv_mask_*sf;
		Function<T,NDIM> result = part1 + part2;
		make_plots(part1,"masked_"+name);
		make_plots(part2,"masked_s"+name);
		return result;
	}

	void make_plots(const Function<T,NDIM> &f,const std::string &name="function")const{
		double width = FunctionDefaults<3>::get_cell_min_width()/2.0 - 1.e-3;
		plot_plane(world,f,name);
		coord_3d start(0.0); start[0]=-width;
		coord_3d end(0.0); end[0]=width;
		plot_line(("line_"+name).c_str(),1000,start,end,f);
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
	const Function<T,NDIM> box_mask_;
	Function<T,NDIM> mask_;
	Function<T,NDIM> inv_mask_;
};



#endif /* SRC_EXAMPLES_SMOOTH_H_ */
