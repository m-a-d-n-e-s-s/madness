/*
 * BSHApply.h
 *
 *  Created on: Feb 26, 2020
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_BSHAPPLY_H_
#define SRC_APPS_CHEM_BSHAPPLY_H_

#include<madness/chem/MolecularOrbitals.h>
#include<madness/world/MADworld.h>

namespace madness {

/// apply the BSH operator on a vector of functions with corresponding potentials

/// this class
///  - constructs the bsh operator with the appropriate exponents
///  - performs a level shift if necessary
///  - adds coupling terms:  ( T - fock(i,i) ) psi_i  = -V psi_i + \sum_{j\neq i} psi_j fock(j,i)
/// TODO: adding a level shift seems to make the operation less precise, why??
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

	/// @param[in]	eps		orbital energies or the square fock matrix
	/// @param[in]	Vpsi	vector of functions V*MOs
	/// @return		an MO structure holding the residuals and the orbital energy updates
	std::tuple<std::vector<Function<T,NDIM> >, Tensor<double> > operator()(
			const std::vector<Function<T,NDIM> > psi,
			const Tensor<T> eps,
			const std::vector<Function<T,NDIM> >& Vpsi1) const {

		double cpu0=cpu_time();
		std::vector<Function<T,NDIM> > Vpsi=copy(world,Vpsi1);

		// add coupling between the functions (i.e. in case of localized orbitals)
		// ( T + V ) \psi_i	= \sum_j \psi_j F_{ji}
		// ( T - F_{ii}) \psi_i = -V \psi_i + \sum_{j\neq i} \psi_jF_{ji}
		// no coupling means F_{ij} =\eps_i \delta_{ij}, and the coupling term vanishes
		Vpsi-=add_coupling_and_levelshift(psi,eps);

	    std::vector < std::shared_ptr<SeparatedConvolution<double,NDIM> > > ops(psi.size());
	    for (int i=0; i<eps.dim(0); ++i) {
	    	T e_i= (eps.ndim()==2) ? eps(i,i) : eps(i);
	    	ops[i]=std::shared_ptr<SeparatedConvolution<double,NDIM> >(
	    			BSHOperatorPtr<NDIM>(world, sqrt(-2.0*eps_in_green(e_i)), lo, bshtol));
	    	ops[i]->destructive()=true;
	    }

	    const std::vector<Function<T,NDIM> > tmp = apply(world,ops,-2.0*Vpsi);
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
	    Tensor<double> delta_eps(psi.size());
	    for (size_t i=0; i<psi.size(); ++i) delta_eps(i)=std::real(in(i))/(norms[i]*norms[i]);

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

	/// limit the orbital energy (diagonal fock matrix element) entering the Green's function parameter mu
	double eps_in_green(const T eps) const {

    	if (std::imag(eps)>1.e-12)
    		MADNESS_EXCEPTION("complex orbital energies in BSHApply",1);
		return std::min(-0.05,std::real(eps)+levelshift);
	}

	std::vector<Function<T,NDIM> > add_coupling_and_levelshift(const std::vector<Function<T,NDIM> > psi,
			const Tensor<T> fock1) const {

		// check dimensions
   	        bool consistent=(psi.size()==size_t(fock1.dim(0)));
		if ((fock1.ndim()==2) and not (psi.size()==size_t(fock1.dim(1)))) consistent=false;

		if (not consistent) {
			print("Fock matrix dimensions",fock1.ndim(), fock1.dim(0), fock1.dim(1));
			print("number of orbitals",psi.size());
			MADNESS_EXCEPTION("confused Fock matrix/orbital energies in BSHApply - 1",1);
		}

        bool do_coupling=(fock1.ndim()==2);

		// subtract the BSH energy (aka Fock diagonal elements) from the rhs
		// ( T - fock(i,i) ) psi_i  = -V psi_i + \sum_{j\neq i} psi_j fock(j,i)
		// if there is no level shift and the orbital energies are large enough the
		// diagonal should simply vanish.
		if (do_coupling) {
			Tensor<T> fock=copy(fock1);
			for (int i=0; i<fock.dim(0); ++i) {
				fock(i,i)-=eps_in_green(fock(i,i));
			}
			return transform(world, psi, fock);

		} else  {
			std::vector<T> eps(psi.size());
			if (fock1.ndim()==1)
				for (int i=0; i<fock1.dim(0); ++i) eps[i]=fock1(i)-eps_in_green(fock1(i));
			if (fock1.ndim()==2)
				for (int i=0; i<fock1.dim(0); ++i) eps[i]=fock1(i,i)-eps_in_green(fock1(i,i));

			std::vector<Function<T,NDIM> > result=copy(world,psi);
			scale(world,result,eps);
			return result;
		}

	}

};
} // namespace madness


#endif /* SRC_APPS_CHEM_BSHAPPLY_H_ */
