/*
 * BSHApply.h
 *
 *  Created on: Feb 26, 2020
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_BSHAPPLY_H_
#define SRC_APPS_CHEM_BSHAPPLY_H_

#include<chem/MolecularOrbitals.h>

/// apply the BSH operator on a vector of functions with corresponding potentials
template<typename T, std::size_t NDIM>
class BSHApply {

public:
	World& world;
	double levelshift=0.0;
	double lo=1.e-6;
	double bshtol=1.e-5;
	bool printme=false;
	bool destroy_Vpsi=false;
	Function<double,NDIM> metric;

public:
	BSHApply(World& world) : world(world),
			metric() {
		bshtol = FunctionDefaults <NDIM> ::get_thresh();
		printme = printme and world.rank()==0;
	}


	/// apply the BSH operator on the vector of functions and respective potentials

	/// @param[in]	arg		the MO structure holding MOs and orbital energies
	/// @param[in]	Vpsi	vector of functions V*MOs
	/// @return		an MO structure holding the residuals and the orbital energy updates
	std::tuple<std::vector<Function<T,NDIM> >, Tensor<double> > operator()(
			const std::vector<Function<T,NDIM> > psi,
			const Tensor<T> eps,
			const std::vector<Function<T,NDIM> >& Vpsi) const {

		double cpu0=cpu_time();

	    std::vector < std::shared_ptr<SeparatedConvolution<double,NDIM> > > ops(psi.size());
	    for (int i=0; i<eps.size(); ++i) {
	    	ops[i]=std::shared_ptr<SeparatedConvolution<double,NDIM> >(
	    			BSHOperatorPtr<NDIM>(world,
	    					sqrt(-2.*std::min(-0.05,eps(i)+levelshift)),
	    					lo, bshtol));
	    	ops[i]->destructive()=destroy_Vpsi;
	    }

	    const std::vector<Function<T,NDIM> > tmp = apply(world,ops,-2.0*Vpsi-2.0*levelshift*psi);
	    const std::vector<Function<T,NDIM> > res=truncate(psi-tmp,FunctionDefaults<NDIM>::get_thresh()*0.1);
	    const std::vector<Function<T,NDIM> > bra_res=make_bra(res);

	    // update eps
	    Tensor<double> norms=real(inner(world,make_bra(tmp),tmp));
	    Tensor<T> rnorms=inner(world,bra_res,res);
	    for (int i=0; i<norms.size(); ++i) {
	    	norms(i)=sqrt(norms(i));
	    	rnorms(i)=sqrt(rnorms(i));
	    }
	    if (printme) {
	    	print("norm2(tmp)",norms);
	    	print("norm2(res)",rnorms);
	    }

	    Tensor<T> in=inner(world,Vpsi,bra_res);	// no shift here!
	    Tensor<double> delta_eps(eps.size());
	    for (int i=0; i<eps.size(); ++i) delta_eps(i)=std::real(in(i))/(norms[i]*norms[i]);

	    if (printme) print("orbital energy update",delta_eps);
	    double cpu1=cpu_time();
	    if (printme) printf("time in BSHApply()  %8.4fs\n",cpu1-cpu0);

	    return std::make_tuple(res,delta_eps);
	}


	/// apply the BSH operator on the vector of functions and respective potentials

	/// @param[in]	arg		the MO structure holding MOs and orbital energies
	/// @param[in]	Vpsi	vector of functions V*MOs
	/// @return		an MO structure holding the residuals and the orbital energy updates
	MolecularOrbitals<T,NDIM> operator()(const MolecularOrbitals<T,NDIM>& arg,
			const std::vector<Function<T,NDIM> >& Vpsi) const {
		std::tuple<std::vector<Function<T,NDIM>, Tensor<T> > > result;
		result=this->operator()(arg.get_mos(),arg.get_eps(),Vpsi);
		return result;
	}

	std::vector<Function<T,NDIM> > make_bra(const std::vector<Function<T,NDIM> >& rhs) const {
		if (metric.is_initialized()) return metric*rhs;
		return rhs;
	}

};


#endif /* SRC_APPS_CHEM_BSHAPPLY_H_ */
