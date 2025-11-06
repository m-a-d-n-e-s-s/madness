/*
 * test_BSHApply.cc
 *
 *  Created on: May 1, 2020
 *      Author: fbischoff
 */

#include<madness/mra/mra.h>
#include<madness/chem/BSHApply.h>
#include<madness/world/test_utilities.h>

using namespace madness;

int factorial(int n) {
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

double hermite(int n, double x) {
	if (n==0) return 1;
	if (n==1) return 2.0*x;
	if (n==2) return 4.0*x*x - 2;
	else throw std::runtime_error("hermite number too large..");
}

template<std::size_t NDIM>
struct HO_function {
	double alpha;
	double normalization;
	std::vector<int> n;
	HO_function(const double m, const double k, const std::vector<int> n)
		: alpha(sqrt(m*k)),n(n) {

		normalization=std::pow(alpha/constants::pi,0.25*double(NDIM));
		for (size_t idim=0; idim<NDIM; ++idim) {
			double term2=1.0/sqrt(std::pow(2.0,double(n[idim])*factorial(n[idim])));
			normalization*=term2;
		}
	}

	double operator()(const Vector<double,NDIM>& r) const {
		double result=normalization;
		for (size_t i=0; i<NDIM; ++i) result*=hermite(n[i],r[i]*sqrt(alpha));
		result*=exp(-0.5*alpha*inner(r,r));
		return result;
	}
};

template<typename T, std::size_t NDIM>
T compute_energy(const Function<double,NDIM>& potential, const Function<T,NDIM>& f) {
	std::vector<Function<T,NDIM> > g=grad(f);
	return inner(f,potential*f) + 0.5*inner(g,g);
}

double compute_energy(const std::vector<int>& n) {
	double ntotal=0.0;
	for (size_t i=0; i<n.size(); ++i) ntotal+=0.5+n[i];
	return ntotal;
}

/// compute the Fock (aka coupling) matrix
template<typename T, std::size_t NDIM>
Tensor<T> compute_fock_matrix(World& world, const Function<double,NDIM>& potential,
		const std::vector<Function<T,NDIM> > f) {

	Tensor<T> fock=matrix_inner(world,f,potential*f);

	for (size_t i=0; i<f.size(); ++i) {
		std::vector<Function<T,NDIM> > gi=grad(f[i]);
		for (size_t j=0; j<f.size(); ++j) {
			std::vector<Function<T,NDIM> > gj=grad(f[j]);
			fock(i,j)+=0.5*inner(gi,gj);
		}
	}
	return fock;
}

Tensor<double_complex> su_complex(const double theta) {
	double_complex z1=std::polar(sqrt(0.5),theta);
	double_complex z2=std::polar(sqrt(0.5),2.0*theta);
	Tensor<double_complex> su(2,2);
	su(0,0)=z1;
	su(0,1)=z2;
	su(1,0)=-std::conj(z2);
	su(1,1)=std::conj(z1);
	Tensor<double_complex> one=inner(su,conj(su),0,0);
	if (std::abs(one.sum()-double_complex(2.0,0.0))>1.e-5) {
		print("error",one.sum());
		print(su);
		print(one);
	}
	return su;
}

Tensor<double> su_real(const double theta) {
	Tensor<double> su(2,2);
	su(0,0)=cos(theta);
	su(0,1)=-sin(theta);
	su(1,0)=sin(theta);
	su(1,1)=cos(theta);
	return su;
}

template<std::size_t NDIM>
std::vector<Function<double,NDIM> > mix(const std::vector<Function<double,NDIM> >& vf, const double theta) {
	MADNESS_ASSERT(vf.size()==2);
	return transform(vf.front().world(),vf,su_real(theta));
}

template<std::size_t NDIM>
std::vector<Function<double_complex,NDIM> > mix(const std::vector<Function<double_complex,NDIM> >& vf, const double theta) {
	MADNESS_ASSERT(vf.size()==2);
	return transform(vf.front().world(),vf,su_complex(theta));
}

/// test the n-dimensional harmonic oscillator
template<typename T, std::size_t NDIM>
int test_converged_function(World& world, double shift, bool coupling) {
	int result=0;
	test_output test("testing NDIM= "+std::to_string(NDIM)+" coupling= "+std::to_string(coupling)+" shift= "+std::to_string(shift));

	// the potential is V(r) = 0.5 k r^2 = 0.5 omega^2 r^2 (m=1)
	const double k=1.0;
	const double m=1.0;

	Function<double,NDIM> potential=FunctionFactory<double,NDIM>(world)
			.functor([&k](const Vector<double,NDIM>& coord_nd){return 0.5*k*inner(coord_nd,coord_nd);});

	// set up and test HO function 1
	std::vector<int> n0(NDIM,0);
	Function<T,NDIM> f0=FunctionFactory<T,NDIM>(world)
			.functor(HO_function<NDIM>(m,k,n0));

	double norm0=f0.norm2();
	T e0=compute_energy(potential,f0);
	if(fabs(norm0-1.0)>1.e-4) MADNESS_EXCEPTION("failure 0",1);
	if(std::abs(T(e0)-T(compute_energy(n0)))>1.e-4) {
		print(e0);
		print(n0,compute_energy(n0));
		MADNESS_EXCEPTION("failure 1",1);
	}

	// set up and test HO function 2
	std::vector<int> n1(NDIM,1);
	Function<T,NDIM> f1=FunctionFactory<T,NDIM>(world)
			.functor(HO_function<NDIM>(m,k,n1));

	double norm1=f1.norm2();
	T e1=compute_energy(potential,f1);
	if(fabs(norm1-1.0)>1.e-4) MADNESS_EXCEPTION("failure 2",1);
	if(std::abs(T(e1)-T(compute_energy(n1)))>1.e-4) {
		print(e1);
		print(n1,compute_energy(n1));
		MADNESS_EXCEPTION("failure 3",1);
	}


	std::vector<Function<T,NDIM> > vf={f0,f1};


	if (coupling) vf=mix(vf,0.3);
	Tensor<T> fock=compute_fock_matrix<T,NDIM>(world,potential,vf);

	// apply the BSH operator
	BSHApply<T,NDIM> bsh_apply(world);
	bsh_apply.levelshift=shift;
	auto [residual,eps_update]=bsh_apply(vf,fock,potential*vf);

	// test residual norm
	double rnorm=norm2(world,residual);
	test.logger<< "rnorm " << rnorm;
	const double thresh=FunctionDefaults<NDIM>::get_thresh();
	bool success=rnorm<sqrt(thresh)*0.1;

	// test orbital energy update
	test.logger << "eps update " <<	 eps_update;
	double loose=(shift!=0.0) ? 10.0 : 1.0; 	 // loosen the threshold if shifts are used
	//success=success and eps_update.normf()*loose;

	test.end(success);
	if (not success) print(fock);
    if (!success) result++;

    return result;
}

template<typename T, std::size_t NDIM>
int test_convergence(World& world, double shift, bool coupling) {
	int result=0;
	test_output test("testing convergence= "+std::to_string(NDIM)+" coupling= "+std::to_string(coupling)+" shift= "+std::to_string(shift));

	// the potential is V(r) = 0.5 k r^2 = 0.5 omega^2 r^2 (m=1)
	const double k_pot=0.05;
	const double m=1.0;
	const double global_shift=-12.0;

	Function<double,NDIM> potential=FunctionFactory<double,NDIM>(world)
			.functor([&k_pot](const Vector<double,NDIM>& coord_nd){return 0.5*k_pot*inner(coord_nd,coord_nd);});

	potential=potential+global_shift;
	// set up and test HO function 1
	std::vector<int> n0(NDIM,0);
	const double k_wf=0.15;
	Function<T,NDIM> f0=FunctionFactory<T,NDIM>(world).functor(HO_function<NDIM>(m,k_wf,n0));
	std::vector<Function<T,NDIM> > vf(1,f0);
	normalize(world,vf);
	Tensor<double> fock=compute_fock_matrix(world,potential,vf);

	double delta=0.0;
	for (int i=0; i<100; ++i) {

		BSHApply<T,NDIM> bsh_apply(world);
		bsh_apply.levelshift=shift;
		auto [residual,eps_update]=bsh_apply(vf,fock,potential*vf);
		vf-=residual;
		double rnorm=norm2(world,residual);
		test.logger <<  "rnorm " << rnorm << " eps " << fock(0,0)-global_shift << " eps-update " <<eps_update[0] << std::endl;
		normalize(world,vf);
		fock(0,0)-=eps_update[0];
		delta=eps_update[0];
	}

//	test.print_and_clear_log();

	bool success=fabs(delta)<1.e-5;
	test.end(success);
    if (!success) result++;
	return result;
}

int main(int argc, char** argv) {
    madness::initialize(argc, argv);

    madness::World world(SafeMPI::COMM_WORLD);
    world.gop.fence();
    startup(world,argc,argv);
    const int k=7;
    const double thresh=1.e-5;
    const double L=24.0;
    FunctionDefaults<1>::set_cubic_cell(-L,L);
    FunctionDefaults<2>::set_cubic_cell(-L,L);
    FunctionDefaults<3>::set_cubic_cell(-L,L);
    FunctionDefaults<1>::set_thresh(thresh);
    FunctionDefaults<2>::set_thresh(thresh);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<1>::set_k(k);
    FunctionDefaults<2>::set_k(k);
    FunctionDefaults<3>::set_k(k);


    int result=0;
    for (bool coupling : {true,false}) {
    	for (double shift : {0.0, -5.0}) {
    	    result+=test_converged_function<double,1>(world,shift,coupling);
    	    result+=test_converged_function<double,2>(world,shift,coupling);
    	    result+=test_converged_function<double,3>(world,shift,coupling);
#ifndef USE_GENTENSOR
    	    result+=test_converged_function<double_complex,1>(world,shift,coupling);
    	    result+=test_converged_function<double_complex,2>(world,shift,coupling);
//    	    result+=test_converged_function<double_complex,3>(world,shift,coupling);
#endif
    	    result+=test_convergence<double,2>(world,shift,coupling);
    	}
    }

    print("result",result);
    madness::finalize();
    return result;

}

