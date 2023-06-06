/*
 * F12Potentials.h
 *
 *  Created on: Aug 4, 2017
 *      Author: kottmanj
 */

#ifndef F12POTENTIALS_H_
#define F12POTENTIALS_H_

#include<madness/chem/PNOGuessFunctions.h>
#include<madness/chem/PNOParameters.h>
#include<madness/chem/PNOStructures.h>

namespace madness {

/// Class that provides all necessary F12 Potentials and Integrals
class F12Potentials {

public:
	World& world;															///< MPI Communicator
	const Nemo& nemo;														///< all SCF data
	const BasisFunctions& basis;											///< all methods to create basis functions
	F12Parameters param;													///< parameters
	std::shared_ptr<real_convolution_3d> fop;								///< operator for convolutions with the correlation factor: f12 = 1/(2*\gamma)*(1-exp(-\gamma*r_{12}))
	std::shared_ptr<real_convolution_3d> slaterop;							///< operator for convolutions with the slater factor exp(-\gamma*r_{12})
	std::shared_ptr<real_convolution_3d> slaterop_sq;						///< operator for convolutions with the slater factor exp(-2.0*\gamma*r_{12})
	std::shared_ptr<real_convolution_3d> coulombop;							///< operator for convolutions with the coulomb operator
	std::shared_ptr<real_convolution_3d> bshop;   							///< operator for convolutions with the bsh operator (parametrized with the exponent of the corrfac: gamma)
	std::shared_ptr<real_convolution_3d> bshop_sq;   						///< operator for convolutions with the bsh operator (parametrized with the exponent of the corrfac: 2*gamma)
	std::vector< std::shared_ptr< real_derivative_3d > > gradop;			///< Gradient Operator
	std::vector< std::shared_ptr< real_convolution_3d > > slatergradop; 	/// Convolution with the gradient of the slater operator
	std::vector< std::shared_ptr< real_convolution_3d > > slatergradop_sq;  /// Convolution with the gradient of the squared (\gamma*2) slater operator
	vector_real_function_3d mos;											///< all molecular orbitals (needed for Q projection parts)
	vector_real_function_3d acmos;											///< Active Molecular Orbitals
	vector_real_function_3d acKmos; 										///< Intermediate for K(active_mos) (takes approx. as much memory than the mos itself. Since #mos<<<#PNOs this should not be a problem
	ParametrizedExchange K;													///< Exchange Operator
	QProjector<double, 3> Q;												///< Projector on virtual space: Q = 1-|k><k|


	F12Potentials(World& world, const Nemo& nemo,const BasisFunctions& basis,const F12Parameters& pp);

	vector_real_function_3d initialize_active_mos(const madness::Nemo& nemo)const{
		vector_real_function_3d ac;
		const vector_real_function_3d& mos=nemo.get_calc()->amo;
		for(size_t i=param.freeze();i<mos.size();++i) ac.push_back(mos[i]);
		MADNESS_ASSERT(ac.size()==mos.size()-param.freeze());
		return ac;
	}

	/// Convenience: Get number of pairs
	size_t npairs()const{
		const size_t n=acmos.size();//mos.size()-param.freeze();
		return n*(n+1)/2;
	}


	/// Convenience: Get electron pair iterator
	ElectronPairIterator pit()const{ return ElectronPairIterator(acmos.size()+param.freeze(),param.freeze());}

	/// apply the regularized potential
	/// \f$ result = <bra|Ue-[f,K]|ket1,ket2> \f$
	/// \f$ <bra|Ue|ket1,ket2> \f$
	/// \f$ <bra|[f,K]|ket1,ket2> = <bra[f,K1]|ket1,ket2> + <bra[f,K2]|ket1,ket2> \f$
	/// \f$ <bra|[f,K1]|ket1,ket2> = <bra|f|K(ket1)>*|ket2> - <K(bra)|f|ket1>*|ket2>
	/// \f$ <bra|[f,K2]|ket1,ket2> = <bra|f|ket1>*|K(ket2)> - |K(<bra|f|ket1>*ket2)>
	vector_real_function_3d apply_regularized_potential(const real_function_3d& moi, const real_function_3d& moj,const real_function_3d Ki, const real_function_3d& Kj,
			const vector_real_function_3d& bra, const vector_real_function_3d& Kpno) const;

	/// compute <ab|g|ij> with the regularized coulomb operator g
	Tensor<double> compute_regularized_fluctuation_matrix(const ElectronPairIterator& it, const vector_real_function_3d& pnos, const vector_real_function_3d& Kpnos)const;
	Tensor<double> compute_regularized_fluctuation_matrix(const std::pair<vector_real_function_3d,vector_real_function_3d> & KPNOA,const std::pair<vector_real_function_3d,vector_real_function_3d> & KPNOB, const std::pair<real_function_3d,real_function_3d>& MKI, const std::pair<real_function_3d,real_function_3d> MKJ)const;

	/// Compute the <ij|op|ab> integrals
	std::valarray<Tensor<double> > compute_ijab_integrals(
			const std::valarray<vector_real_function_3d>& functions,
			const std::shared_ptr<real_convolution_3d>& op) const;

	/// Compute the <xy|op|ab> integrals, op==NULL pointer means overlap: <xy|ab>
	Tensor<double> compute_xyab_integrals(const real_function_3d& x, const real_function_3d& y, const vector_real_function_3d& a, const vector_real_function_3d& b, const std::shared_ptr<real_convolution_3d>& op,
			const bool diagonal = false) const;

	/// Compute the <ij|g|ab> integrals
	std::valarray<Tensor<double> > compute_gijab_integrals(const std::valarray<vector_real_function_3d>& functions)const{
		return compute_ijab_integrals(functions,coulombop);
	}

	/// Compute the <ij|f|ab> integrals
	std::valarray<Tensor<double> > compute_fijab_integrals(const std::valarray<vector_real_function_3d>& functions)const{
		return compute_ijab_integrals(functions,fop);
	}
	Tensor<double> compute_fijab_integrals(const size_t &i, const size_t& j,const vector_real_function_3d& functions)const{
		return compute_xyab_integrals(acmos[i], acmos[j], functions, functions, fop, i==j);
	}

	/// Compute the <ij|fQg|ij> and |ij|fQg|ji> integrals
	std::valarray<double> compute_fQg_integrals() const;

	/// Compute the <ij|fQU|ij> and |ij|fQU|ji> integrals
	std::valarray<double> compute_fQU_integrals() const;

	/// Compute the fQc integrals with ABS used as: Q12 = Q*|a><a|
	/// K1f: <i|fja|x>, with fja = <a|f|j> and x=QK(<a|f|j>*i)
	/// K2f: <i|fja|y>, with y = Q(<Ka|f|j>*i
	/// fK2: <i|fja|z>, with z = Q(<a|f|Kj>*i) (note: actually no abs needed for fK1)
	/// fK1: <i|fja|w>, with w = Q(<a|f|j>*Ki)
	std::valarray<double> compute_fQc_integrals(
			const vector_real_function_3d& Kmos,
			const std::valarray<vector_real_function_3d>& functions) const;
	std::pair<double, double> compute_fQc_integrals_ij(
			const vector_real_function_3d& Kmos,
			const vector_real_function_3d& functions,
			const ElectronPairIterator& it,
			const bool& use_no_intermediates = false) const;


	/// Compute the <ab|[f,K]|ij> integrals
	std::valarray<Tensor<double> > compute_cijab_integrals(const std::valarray<vector_real_function_3d>& functions) const;
	Tensor<double> compute_cijab_integrals(const size_t &i, const size_t& j, const real_function_3d& Ki, const real_function_3d& Kj , const vector_real_function_3d& functions)const;

	/// Compute the <ab|Ue|ij> Integrals
	Tensor<double> compute_uijab_integrals(const size_t&i, const size_t&j, const vector_real_function_3d& functions)const{
		const vector_real_function_3d vUij = convolve_with_U(functions,acmos[i],acmos[j]);
		const Tensor<double> result = matrix_inner(world,functions,vUij);
		return result;
	}

	/// convenience for formated printing
	/// @param[in] es: singlet pair energies
	/// @param[in] et: triplet pair energies
	void print_pair_energies(const PNOPairs& pairs,
			const std::string& msg = "") const;
	void print_pair_energies(const std::valarray<double>& es, const std::valarray<double>& et, const std::string& msg = "",const PairType& type=UNKNOWN_PAIRTYPE) const;

	/// print_f12_energies with PairEnergies POD
	void print_f12_energies(const PairEnergies& pe, const std::string& msg="")const{
		print_pair_energies(pe.eijs_f12,pe.eijt_f12,msg);
	}

	// compute the f12 energy of the MP2 equivalent of CIS(D) correction which is the S2b and S2c part
	// we do not use hylleras type right now, but in analogy to 6D MRA just the f12Q12g12 part (this might not be completely correct)
	// result: (2<xij| + <jxi)( f12Q12g12|x_i,j> + f12Q12g12|ix_j> - fQOxg|ij> - fOxQg|ij>
	PairEnergies compute_cispd_f12_energies(const vector_real_function_3d xcis)const{
		PairEnergies result(npairs());
		double f12_energy=0.0;
		for(ElectronPairIterator it=pit();it;++it){
			const real_function_3d& moi=acmos[it.i()];
			const real_function_3d& moj=acmos[it.j()];
			const real_function_3d& xi=xcis[it.i()];
			const real_function_3d& xj=xcis[it.j()];


			// part 1. S2b part
			double part1=0.0;
			if(it.diagonal()){
				part1 = 2.0*(compute_fQg_integral(xi,moj,xi,moj)+ compute_fQg_integral(xi,moj,moi,xj)); // this is the |xi,j> + |i,xj> part
				// now the Ox part: -(O1x+O2x)|ij>
				vector_real_function_3d p1,p2;
				p1 = apply(world,*fop,acmos*moj)*moi;
				p1=append(p1,xcis);
				p2 = xcis;
				p2=append(p2,apply(world,*fop,acmos*moi)*moj);
				part1 -= 2.0*(madness::inner(world,xi*p1,apply(world,*coulombop,p2*moj)).sum());
			}
			else{
				MADNESS_EXCEPTION("Still wrong ... also part2 ",1);
				part1 = 2.0*(compute_fQg_integral(moi,moj,xi,moj) + compute_fQg_integral(moi,moj,moi,xj) - compute_fQg_integral(moj,moi,xi,moj) + compute_fQg_integral(moj,moi,moi,xj)) ;
			}

			// part2: S2c part
			double part2=0.0;
			if(it.diagonal()){

				const vector_real_function_3d p1=xcis;
				const vector_real_function_3d p2=Q(apply(world,*coulombop,moi*acmos)*moj);
				part2 = 2.0*( madness::inner(world,p1*xi,apply(world,*fop,p2*moj)).sum() + madness::inner(world,p1*moi,apply(world,*fop,p2*xj)).sum()  ); // |xi,j> + |i,xj> part

				// now projector part
				// transform mos: |ti> = sum_k |k> <x_k|x_i> = Ox^t|x>
				Projector<double,3> Oxt(xcis,acmos);
				// OQ part
				vector_real_function_3d pp1=Oxt(p1);
				vector_real_function_3d pp2=p2;
				// QO part
				pp1=append(pp1,p1);
				pp2=append(pp2,Oxt(p2));
				part2 -= 2.0*( madness::inner(world,p1*moi,apply(world,*fop,p2*moj))).sum();
			}else{
				MADNESS_EXCEPTION("Still wrong ... also part1 ",1);
				// use this intermediates and avoid calculating them twice
				const vector_real_function_3d fi=apply(world,*fop,moi*acmos);
				const vector_real_function_3d fj=apply(world,*fop,moj*acmos);
				{
				        //const vector_real_function_3d& p1=acmos;
					const vector_real_function_3d& p2=Q(apply(world,*coulombop,moi*acmos)*moj);
					part2 += 2.0*madness::inner(world,fi,moj*p2).sum() - madness::inner(world,fj,moi*p2).sum();
				}
				{
					const vector_real_function_3d& p1=Q(apply(world,*coulombop,moj*acmos)*moi);
					// const vector_real_function_3d& p2=acmos;
					part2 += 2.0*madness::inner(world,moi*p1,fj).sum() - madness::inner(world,moj*p1,fi).sum();
				}
			}
			double factor=1.0;
			if(it.diagonal()) factor=0.5; // cancels the factor 2 above ... probably stupid to implement that way
			const double result_ij = factor*(part1 - part2);

			std::cout << "\npart1=" << part1 << "\npart2=" << part2 << "\n";

			f12_energy += result_ij;
		}
		result.energy_f12=f12_energy;
		return result;
	}


	// efficient computation using Q12 = 1 - O1 - O2 + O12 = 1 - O1q2 - q1O2 and q1 = 1 - 0.5*O1
	double compute_fQg_integral(const real_function_3d bra1,
			const real_function_3d& bra2, const real_function_3d& ket1,
			const real_function_3d& ket2, const bool& diagonal = false) const;

	// Call the correct functions depending which parameters are set
	PairEnergies compute_f12_energies(
			const std::valarray<vector_real_function_3d>& pnos) const;
	PairEnergies compute_f12_energies() const;

	/// return 2.0 <ij|fQg|ij> - <ji|fQg|ij> as well as singlet and triplets
	PairEnergies compute_projected_f12_energies() const;

	std::vector<real_function_3d> read_cabs_from_file(const std::string& filename)const;

	PairEnergies compute_hylleraas_f12_energies(
			const std::valarray<vector_real_function_3d>& pnos) const;


	PairEnergies compute_f12_pair_energy(const std::valarray<vector_real_function_3d>& pnos, const std::valarray<vector_real_function_3d>& abs) const;

	/// Compute the f12 pair energy if Ansatz 2 is used
	/// This means: Regularized MRA-PNO-MP2-F12 calculation
	/// We use RI for the integrals with the exchange commutator
	/// in two steps: 1. PNOs as approx for Q, 2. Correction with ABS set (either auto generated or from file)
	/// EF12 =<ij|(fQg + fQU - fQ[f,K])(2|ij> - |ji>)
	PairEnergies compute_f12_pair_energies(const std::valarray<vector_real_function_3d>& abs) const;

	/// Compute the f12 pair energy if Ansatz 1 is used
	/// This means: Normal MRA-PNO-MP2 calculation + F12 Correction
	/// Old Code
	/// Currently: Neglect Commutator if no abs is given
	/// if abs_u is true then fQU part is also neglected (RI with PNOs as virtuals: fQU ~= fVU => fQU - fVU ~= 0)
	/// EF12 = <ij|(fQg-fVg + fUg - fVg)|2ij-ji>
	PairEnergies compute_f12_correction(
			const std::valarray<vector_real_function_3d>& pnos, const std::valarray<vector_real_function_3d>& abs) const;

	/// Compute the applied Ue potential
	/// \f$ <bra|Ue|ket1,ket2>_1 =  <bra|Ulocal|ket1,ket2> + <bra|Unlocal|ket1,ket2> \f$
	real_function_3d convolve_with_U(const real_function_3d& bra,
			const real_function_3d& ket1, const real_function_3d& ket2,
			const bool symmetric = false) const;
	vector_real_function_3d convolve_with_U(const vector_real_function_3d& bra,
			const real_function_3d& ket1, const real_function_3d& ket2,
			const bool symmetric = false) const;

	/// @return \f$ \int dr2 (1-e^{-gamma*r12})/r12 + gamma/2* e^{-\gamma*r12} *function(r2)  \f$
	/// we use: 1-e^{-gamma*r12})/r12 = 2*gamma*f12*g12
	real_function_3d convolve_with_local_U(
			const real_function_3d& function) const{
		vector_real_function_3d dummy(1,function);
		return convolve_with_local_U(dummy).front();
	}

	vector_real_function_3d convolve_with_local_U(const vector_real_function_3d& functions)const;

	/// @return convolution with the slater potential which is
	/// \f$ V_{slater} = \exp{(-\gamma*r_{12})}  \f$
	/// @param[in] function -> operator is applied onto this funciton
	/// if nemo-style orbitals are used the bra element has to be multiplied with R^2 beforehand
	vector_real_function_3d convolve_with_slater_potential(const vector_real_function_3d& functions)const{
		return apply(world,*slaterop,functions);
	}


	/// @return \f$ \int dx2 f12*g12*f(2)  \f$
	/// @param[in] function: function on which f*g operator is applied
	/// right now we use: fg = 1/(2*gamma)*(g - 4*pi*G(gamma))
	/// todo: implement fg operator directly (may be cheaper)
	real_function_3d convolve_with_fg(const real_function_3d& function)const;
	vector_real_function_3d convolve_with_fg(const vector_real_function_3d& functions)const;

	/// @return \f$ -1.0/2.0*<bra(1)|e^{-gamma*r_{12}}/r_{12}* \vec{r_{12}}\cdot (\nabla_1 - \nabla_2)[axis]|ket1,ket2>_1  \f$
	/// @param[in] bra: no derivatives (for nemo-style orbitals this should be multiplied with R^2 beforehand)
	/// @param[in] ket1: this is the part (for electron 1) where the derivative operators will partially act on and over which the integration variable will run
	/// @param[in] ket2: this is the part (for electron 2) where the derivative operators will partially act on, no integration over this variable
	/// @param[in] symmetric: if ket1 and ket2 are the same functions (little speedup because the gradient has to be calculated only once)
	/// @param[in] squared: factor for gamma: a*\gamma. Default is just 1.0. if true then a=2
	/// We use the following relations:
	///   1.0/2.0*e^{-\gamma*r_{12}}*\vec{r_{12}}/r_{12}
	///=  1.0/(2.0*\gamma) \nabla_2(e^{-\gamma*r_{12}})
	///= -1.0/(2.0*\gamma) \nabla_1(e^{-\gamma*r_{12}})
	///= -1.0/(2.0*\gamma) \nabla_{12}(e^{-\gamma*r_{12}} )
	/// so we can use the gradient of the SlaterOperator Grad(e^{-\gamma*r}
	/// with this we have
	/// \f$ -1.0/2.0*<bra|e^{-gamma*r12}/r12* \vec{r12}\cdot (\nabla_1 - \nabla_2)[axis]|ket>
	/// =    1.0/(2.0*\gamma)*<bra|Grad(e^{-\gamma*r_{12}})\cdot (\nabla_1 - \nabla_2)[axis]|ket1,ket2>_{particle}   \f$
	/// todo: Partial integration may lead to more efficient scheme without the Gradient of the SlaterOperator
	real_function_3d convolve_with_nonlocal_U(const real_function_3d& bra, const real_function_3d& ket1,const real_function_3d& ket2, const bool symmetric=false, const bool& squared=false)const;
	vector_real_function_3d convolve_with_nonlocal_U(const vector_real_function_3d& bra, const real_function_3d& ket1,const real_function_3d& ket2, const bool symmetric=false, const bool& squared=false)const;

	/// Convolution with f*U
	/// This breaks apart into the convolution with:
	/// \f$ f*U = 1/(2\gamma)*(U - e^{-\gamma*r_{12}} \f$
	/// local part: (could be further simplified by combining the BSH Operators of Ulocal the expression below)
	/// \f$ f*U^{\mbox{loc}} =  1/(2\gamma)*(U^{\mbox{loc}} - 4\pi*G(\gamma) + 4\pi*G(2\gamma) - \gamma/2 (e^{-2\gamma*r_{12}} )
	/// nonlocal part:
	/// \f$ f*U^{\mbox{nloc}} = 1/(2\gamma)*(U^{\mbox{nloc}}(\gamma) - U^{\mbox{nloc}}(2\gamma))
	vector_real_function_3d convolve_with_fU(const vector_real_function_3d& bra, const real_function_3d& ket1,
			const real_function_3d& ket2, const bool symmetric = false) const;
	real_function_3d convolve_with_fU(const real_function_3d& bra,const real_function_3d& ket1,
			const real_function_3d& ket2, const bool symmetric = false) const{
		vector_real_function_3d dummy(1,bra);
		return convolve_with_fU(dummy,ket1,ket2,symmetric).front();
	}

	/// convolve with the gradient of the slater potential or the squared slater potential
	vector_real_function_3d convolve_with_gradslater(const vector_real_function_3d functions, const size_t axis, const bool& squared=false)const{
		if(squared){
			return apply(world,*slatergradop_sq[axis],functions);
		}else{
			return apply(world,*slatergradop[axis],functions);
		}

	}

	/// Inner product for X_ab Tensors. Convenience
	double inner(const Tensor<double>&A, const Tensor<double>&B)const{
		return A.trace(B);
	}


};

} /* namespace madness */

#endif /* F12POTENTIALS_H_ */
