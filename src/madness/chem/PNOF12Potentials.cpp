/*
 * F12Potentials.cc
 *
 *  Created on: Aug 4, 2017
 *      Author: kottmanj
 */

#include<madness/chem/PNOF12Potentials.h>

namespace madness {

/// Factory function generating operator for convolution with grad(e^{-\gamma*r}) in 3D

/// Returns a 3-vector containing the convolution operator for the
/// x, y, and z components of grad(e^{-\gamma*r})
/// maybe move to madness library (operator.h) at some point
static
inline
std::vector< std::shared_ptr< SeparatedConvolution<double,3> > >
GradSlaterOperator(World& world,
		double gamma,
		double lo,
		double eps,
		const BoundaryConditions<3>& bc=FunctionDefaults<3>::get_bc(),
		int k=FunctionDefaults<3>::get_k())
		{
	typedef SeparatedConvolution<double,3> real_convolution_3d;
	typedef std::shared_ptr<real_convolution_3d> real_convolution_3d_ptr;
	const double pi = constants::pi;
	const Tensor<double> width = FunctionDefaults<3>::get_cell_width();
	double hi = width.normf(); // Diagonal width of cell
	const bool isperiodicsum = (bc(0,0)==BC_PERIODIC);
	if (isperiodicsum) hi *= 100; // Extend range for periodic summation

	GFit<double,3> fit=GFit<double,3>::SlaterFit(gamma,lo,hi,eps,false);
	Tensor<double> coeff=fit.coeffs();
	Tensor<double> expnt=fit.exponents();

	if (bc(0,0) == BC_PERIODIC) {
		fit.truncate_periodic_expansion(coeff, expnt, width.max(), true);
	}

	int rank = coeff.dim(0);

	std::vector<real_convolution_3d_ptr> gradG(3);
  const auto lattice_ranges = bc_lattice_ranges<3>(bc);

	for (int dir=0; dir<3; dir++) {
		std::vector< ConvolutionND<double,3> > ops(rank);
		for (int mu=0; mu<rank; mu++) {
			// We cache the normalized operator so the factor is the value we must multiply
			// by to recover the coeff we want.
			double c = std::pow(sqrt(expnt(mu)/pi),3); // Normalization coeff
			ops[mu].setfac(coeff(mu)/c/width[dir]);

			for (int d=0; d<3; d++) {
				if (d != dir)
					ops[mu].setop(d,GaussianConvolution1DCache<double>::get(k, expnt(mu)*width[d]*width[d], 0, lattice_ranges[dir]));
			}
			ops[mu].setop(dir,GaussianConvolution1DCache<double>::get(k, expnt(mu)*width[dir]*width[dir], 1, lattice_ranges[dir]));
		}
		gradG[dir] = real_convolution_3d_ptr(new SeparatedConvolution<double,3>(world, ops));
	}

	return gradG;
		}


F12Potentials::F12Potentials(World& world,const Nemo& nemo, const BasisFunctions& basis, const F12Parameters& pp) :
														world(world),
														nemo(nemo),basis(basis),
														param(pp),
														mos(nemo.get_calc()->amo),
														acmos(initialize_active_mos(nemo)),
														K(ParametrizedExchange(world, nemo, pp.exchange())),
														Q(nemo.get_calc()->amo) {
	const double lo = 1.e-6;
	const double eps = param.op_thresh();
	coulombop = std::shared_ptr < real_convolution_3d > (CoulombOperatorPtr(world, lo, eps));
	bshop = std::shared_ptr < real_convolution_3d > (BSHOperatorPtr3D(world, param.gamma(), lo, eps));
	bshop_sq = std::shared_ptr < real_convolution_3d > (BSHOperatorPtr3D(world, param.gamma() * 2.0, lo, eps));
	fop = std::shared_ptr < real_convolution_3d > (SlaterF12OperatorPtr(world, param.gamma(), lo, eps));
	slaterop = std::shared_ptr < real_convolution_3d > (SlaterF12OperatorPtr(world, param.gamma(), lo, eps));
	slaterop_sq = std::shared_ptr < real_convolution_3d > (SlaterF12OperatorPtr(world, param.gamma() * 2.0, lo, eps));
	fop = std::shared_ptr < real_convolution_3d > (SlaterF12OperatorPtr(world, param.gamma(), lo, eps));

	//test
	{
		MyTimer time = MyTimer(world).start();
		real_convolution_3d slater = SlaterOperator(world, param.gamma(), lo, eps);
		real_convolution_3d slater_sq = SlaterOperator(world, param.gamma() * 2.0, lo, eps);
		const real_function_3d mo = nemo.get_calc()->amo.front();
		const real_function_3d momo = (mo * mo).truncate();
		const double test1 = slater(momo).norm2();
		const double test2 = (*slaterop)(momo).norm2();
		const double test3 = (*fop)(momo).norm2();
		const double test4 = (momo).trace() * pow(FunctionDefaults < 3 > ::get_cell_min_width(), 3.0);
		const double factor = 1.0 / (param.gamma() * 2.0);
		MADNESS_ASSERT((test1 - test2) < 0.1 * FunctionDefaults < 3 > ::get_thresh());
		MADNESS_ASSERT(test3 - factor * (test4 - test2) < 0.1 * FunctionDefaults < 3 > ::get_thresh());
		if (world.rank() == 0)
			std::cout << "F12 Operators initialized right\n";
		time.stop().print("Testing F12 Operators");

		MyTimer timeK = MyTimer(world).start();
		MADNESS_ASSERT(acmos.size()==nemo.get_calc()->amo.size()-param.freeze());
		acKmos = K(acmos);
		timeK.stop().print("Initializing (active) Exchange Orbitals");

		// give information on the used memory
		const double size_m=get_size(world,acmos);
		const double size_k=get_size(world,acKmos);
		if(world.rank()==0) std::cout << "Size of active mos and their exchange potentials:\n"
				<< "K(acmos)=" << size_k << " Gbyte\n"
				<< "acmos   =" << size_m << " Gbyte\n";

	}

	gradop = gradient_operator<double, 3>(world);
	slatergradop = GradSlaterOperator(world, param.gamma(), lo, eps);
	slatergradop_sq = GradSlaterOperator(world, param.gamma() * 2.0, lo, eps);

}

Tensor<double> F12Potentials::compute_xyab_integrals(const real_function_3d& x, const real_function_3d& y, const vector_real_function_3d& a, const vector_real_function_3d& b, const std::shared_ptr<real_convolution_3d>& op,
		const bool diagonal) const {
	const double muleps=0.0;
	vector_real_function_3d iv = mul_sparse(world, x, a,  muleps);
	truncate(world, iv);
	vector_real_function_3d jv;
	if (diagonal)
		jv = iv;
	else {
		jv = mul_sparse(world, y, b,  muleps);
		truncate(world, jv);
	}
	vector_real_function_3d opiv;
	if (op != NULL)
		opiv = apply(world, *op, iv);
	else
		opiv = iv;

	Tensor<double> result_ij = matrix_inner(world, opiv, jv, diagonal);
	return result_ij;
}

PairEnergies F12Potentials::compute_f12_pair_energy(const std::valarray<vector_real_function_3d>& pnos, const std::valarray<vector_real_function_3d>& cabs) const {
	// F12 is used to regularize
	// so here we compute the f12 part which is missing from the regularized mp2 energies
	if(world.rank()==0) std::cout << "Computing F12 Energy (Ansatz 2)\n";
	std::valarray<vector_real_function_3d> abs(pnos.size());
	if(cabs.size()!=pnos.size()){
		if(world.rank()==0) std::cout << "No CABS given using only the PNOs as ABS\n";
		abs = pnos;
	}
	else{
		for(ElectronPairIterator it=pit();it;++it){
			vector_real_function_3d tmp = pnos[it.ij()];
			for(const auto& x:cabs[it.ij()]) tmp.push_back(x);
			abs[it.ij()]=tmp;
			if(world.rank()==0) std::cout << "Constructing ABS for pair "
					<< it.name() << " from " << pnos[it.ij()].size() << " PNOs and "
					<< cabs[it.ij()].size() << " CABS functions\n";
		}
	}
	return compute_f12_pair_energies(abs);
}

PairEnergies F12Potentials::compute_f12_correction(const std::valarray<vector_real_function_3d>& pnos,const std::valarray<vector_real_function_3d>& abs) const {

	PairEnergies result(npairs());

	const std::valarray<double> fQg = compute_fQg_integrals();
	std::valarray<double> fQU;
	if (!param.abs_u())
		fQU = compute_fQU_integrals();

	std::valarray < Tensor<double> > fijab = compute_fijab_integrals(pnos);
	std::valarray < Tensor<double> > gijab = compute_gijab_integrals(pnos);

	std::valarray < Tensor<double> > cijab_abs = compute_cijab_integrals(abs);
	std::valarray < Tensor<double> > cijab = compute_cijab_integrals(pnos);
	std::valarray < Tensor<double> > fijab_abs = compute_fijab_integrals(abs);

	std::valarray<double> es(npairs()); // singlet f12 energies
	std::valarray<double> et(npairs()); // triplet f12 energies
	std::valarray<double> eij(npairs());// total f12 energies
	std::valarray<double> esg(npairs()); // singlet f12 energies: fQg part
	std::valarray<double> esu(npairs()); // singlet f12 energies: fQU part
	std::valarray<double> esc(npairs()); // singlet f12 energies: Commutator part
	std::valarray<double> etg(npairs()); // triplet f12 energies: fQg part
	std::valarray<double> etu(npairs()); // triplet f12 energies: fQU part
	std::valarray<double> etc(npairs()); // triplet f12 energies: Commutator part
	// same for the V12 part of the projector
	std::valarray<double> esvg(npairs()); // singlet f12 energies: fQg part
	std::valarray<double> esvu(npairs()); // singlet f12 energies: fQU part
	std::valarray<double> esvc(npairs()); // singlet f12 energies: Commutator part
	std::valarray<double> etvg(npairs()); // triplet f12 energies: fQg part
	std::valarray<double> etvu(npairs()); // triplet f12 energies: fQU part
	std::valarray<double> etvc(npairs()); // triplet f12 energies: Commutator part
	double est=0.0;
	double ett=0.0;
	for (ElectronPairIterator it=pit(); it; ++it) {
		const double fvg = inner(fijab[it.ij()], gijab[it.ij()]);
		const double gij = fQg[it.ij()];
		double gji = 0.0;
		double Uij = 0.0;
		double Uji = 0.0;
		double fvu = 0.0;
		double fvux = 0.0;
		double fvgx = 0.0;
		double cij =0.0;
		double fvc =0.0;
		double cji = 0.0;
		double fvcx =0.0;
		Tensor<double> uijab = compute_uijab_integrals(it.i(), it.j(), pnos[it.ij()]);
		fvu = inner(fijab[it.ij()], uijab);
		Tensor<double> uijab_abs;
		if(param.abs_u()){
			uijab_abs = compute_uijab_integrals(it.i(),it.j(),abs[it.ij()]);
			Uij = inner(fijab_abs[it.ij()],uijab_abs);
		}
		else{
			Uij = fQU[it.ij()];
		}
		if(param.abs_c()){
			// neglect if there is no abs
			cij =0.0;
			fvc =0.0;
		}
		else{
			cij = inner(fijab_abs[it.ij()],cijab_abs[it.ij()]);
			fvc = inner(fijab[it.ij()],cijab[it.ij()]);
		}
		if (it.diagonal()) {

			// additional information
			esg[it.ij()]=gij;
			esu[it.ij()]=Uij;
			esc[it.ij()]=cij;
			etg[it.ij()]=0.0;
			etu[it.ij()]=0.0;
			etc[it.ij()]=0.0;
			esvg[it.ij()]= -1.0*(fvg);
			esvu[it.ij()]= -1.0*(fvu);
			esvc[it.ij()]= -1.0*(fvc);
			etvg[it.ij()]=0.0;
			etvu[it.ij()]=0.0;
			etvc[it.ij()]=0.0;

			es[it.ij()] = gij-fvg + Uij-fvu - cij + fvc;
			et[it.ij()] = 0.0;
		} else {

			if(param.abs_c()){
				// neglect if there is no abs
				cij =0.0;
				fvc =0.0;
			}
			else{
				cji  = inner(fijab_abs[it.ij()],transpose(cijab_abs[it.ij()]));
				fvcx = inner(fijab[it.ij()],transpose(cijab[it.ij()]));
			}
			Tensor<double> gjiab = transpose(gijab[it.ij()].reshape(pnos[it.ij()].size(), pnos[it.ij()].size()));
			fvgx = inner(fijab[it.ij()], gjiab);
			gji = fQg[it.ij() + it.npairs()];

			Tensor<double> ujiab = transpose(uijab.reshape(pnos[it.ij()].size(), pnos[it.ij()].size()));
			fvux = inner(fijab[it.ij()], ujiab);
			if (param.abs_u()) {
				Uij = inner(fijab_abs[it.ij()],transpose(uijab_abs));
			}else{
				Uij = fQU[it.ij() + it.npairs()];
			}

			// additional information factor 1/2 missing because ij <-> ji pairs
			esg[it.ij()]=gij+gji;
			esu[it.ij()]=Uij+Uji;
			esc[it.ij()]=cij+cji;
			etg[it.ij()]=3.0*(gij-gji);
			etu[it.ij()]=3.0*(Uij-Uji);
			etc[it.ij()]=3.0*(cij-cji);
			esvg[it.ij()]= -1.0*(fvg + fvgx);
			esvu[it.ij()]= -1.0*(fvu + fvux);
			esvc[it.ij()]= -1.0*(fvc + fvcx);
			etvg[it.ij()]= -1.0*(fvg - fvgx);
			etvu[it.ij()]= -1.0*(fvu - fvux);
			etvc[it.ij()]= -1.0*(fvc - fvcx);
			// real singlet and triplet energies
			es[it.ij()] = ((gij-fvg + gji-fvgx) + (Uij-fvu + Uji-fvux) - (cij - fvc + cji - fvcx));
			et[it.ij()] = 3.0* ((gij-fvg - gji+fvgx) + (Uij-fvu - Uji+fvux) - (cij - fvc - cji + fvcx));
		}
		est += es[it.ij()];
		ett += et[it.ij()];
		eij[it.ij()] = es[it.ij()] + et[it.ij()];

	}

	// information about the different parts
	print_pair_energies(esg,etg  ,"fQg energy");
	print_pair_energies(esu,etu  ,"fQU energy");
	print_pair_energies(esc,etc  ,"fQc energy");
	print_pair_energies(esvg,etvg,"fVg energy");
	print_pair_energies(esvu,etvu,"fVU energy");
	print_pair_energies(esvc,etvc,"fVc energy");


	// the real deal
	print_pair_energies(es,et,"F12 energy");

	result.eijs_f12 = es;
	result.eijt_f12 = et;
	result.energy_f12=est+ett;

	return result;

}

PairEnergies F12Potentials::compute_f12_pair_energies(const std::valarray<vector_real_function_3d>& abs) const {

	PairEnergies result(npairs());

	const std::valarray<double> fQg = compute_fQg_integrals();
	const std::valarray<double> fQU= compute_fQU_integrals();\
	//const vector_real_function_3d Kmos = K(mos);
	const std::valarray<double> fQc= compute_fQc_integrals(acKmos,abs);
	double energy=0.0;
	//std::valarray < Tensor<double> > fijab = compute_fijab_integrals(abs);
	//std::valarray < Tensor<double> > cijab = compute_cijab_integrals(abs);
	std::valarray<double> es(npairs()); // singlet f12 energies
	std::valarray<double> et(npairs()); // triplet f12 energies

	std::valarray<double> esg(npairs()); // singlet f12 energies: fQg part
	std::valarray<double> esu(npairs()); // singlet f12 energies: fQU part
	std::valarray<double> esc(npairs()); // singlet f12 energies: Commutator part
	std::valarray<double> etg(npairs()); // triplet f12 energies: fQg part
	std::valarray<double> etu(npairs()); // triplet f12 energies: fQU part
	std::valarray<double> etc(npairs()); // triplet f12 energies: Commutator part

	for (ElectronPairIterator it=pit(); it; ++it) {
		const double cij = fQc[it.ij()];
		const double gij = fQg[it.ij()];
		const double Uij = fQU[it.ij()];
		//double fvu = 0.0;
		//double fvux= 0.0;
		double cji = fQc[it.ij()+it.npairs()];
		double gji = fQg[it.ij()+it.npairs()];
		double Uji = fQU[it.ij()+it.npairs()];
		if (it.diagonal()) {
			// additional information
			esg[it.ij()]=gij;
			esu[it.ij()]=Uij;
			esc[it.ij()]=cij;
			etg[it.ij()]=0.0;
			etu[it.ij()]=0.0;
			etc[it.ij()]=0.0;
			// singlet and triplet energies
			es[it.ij()] = gij+ Uij - cij;
			et[it.ij()] = 0.0;
		} else {
			// additional information
			esg[it.ij()]=gij+gji;
			esu[it.ij()]=Uij+Uji;
			esc[it.ij()]=cij+cji;
			etg[it.ij()]=3.0*(gij-gji);
			etu[it.ij()]=3.0*(Uij-Uji);
			etc[it.ij()]=3.0*(cij-cji);
			// singlet and triplet energies
			es[it.ij()] = ((gij + gji) + (Uij + Uji) - (cij + cji));
			et[it.ij()] = 3.0*((gij - gji) + (Uij - Uji) - (cij - cji));
		}
		energy+= es[it.ij()] + et[it.ij()];
	}


	// information about the different parts
	print_pair_energies(esg,etg  ,"fQg energy");
	print_pair_energies(esu,etu  ,"fQU energy");
	print_pair_energies(esc,etc  ,"fQc energy");

	// the real deal
	print_pair_energies(es,et,"F12 energy");
	result.eijs_f12=es;
	result.eijt_f12=et;
	result.energy_f12=energy;

	return result;
}

void F12Potentials::print_pair_energies(const std::valarray<double>& es, const std::valarray<double>& et, const std::string& msg, const PairType& type) const{
  //double est = 0.0;
  //double ett = 0.0;
  //double est_f12 = 0.0;
  //double ett_f12 = 0.0;
	if (world.rank() == 0) {
		std::cout << std::setfill(' ');
		std::cout << std::scientific << std::setprecision(5);
		double est = 0.0, ett = 0.0; // total singlet and triplet energies
		print("=======================================================");
		print("\n" + msg + "\n");
		print("=======================================================");
		if(type==CISPD_PAIRTYPE)
			std::cout << " "  << std::setw(4) << "i" << " " << std::setw(4) << "j" << " " << std::setw(12) << "CORR GS" << " " << std::setw(12) << "CORR ES" << " " << std::setw(12) << "Total" << "\n";
		else
			std::cout << " "  << std::setw(4) << "i" << " " << std::setw(4) << "j" << " " << std::setw(12) << "singlet" << " " << std::setw(12) << "triplet" << " " << std::setw(12) << "Total" << "\n";
		print("=======================================================");
		for (ElectronPairIterator it=pit(); it; ++it) {
			std::cout << std::setfill(' ') << std::fixed   << " "
					<< std::setw(4) << it.i()+param.freeze() <<" "
					<< std::setw(4) << it.j()+param.freeze() <<" "
					<< std::scientific << std::setprecision(5) <<  " "
					<< std::setw(12)<<  es[it.ij()]  <<   " "
					<< std::setw(12)<<  et[it.ij()]  <<  " "
					<< std::setw(12)<<  (es[it.ij()] + et[it.ij()]) << "\n";
			est += es[it.ij()];
			ett += et[it.ij()];
		}
		print("=======================================================");
		if(type==CISPD_PAIRTYPE){
			std::cout <<msg + " (Sum GS) : " << std::setw(12)<< est << "\n";
			std::cout <<msg + " (Sum ES) : " << std::setw(12)<< ett << "\n";
			std::cout <<msg + " (∆CIS(D)): " << std::setw(12)<< est + ett  << "\n";
		}else{
			std::cout <<msg + " (singlet): " << std::setw(12)<< est << "\n";
			std::cout <<msg + " (triplet): " << std::setw(12)<< ett << "\n";
			std::cout <<msg + " (total)  : " << std::setw(12)<< est + ett  << "\n";
		}
		print("=======================================================");
	}
}


vector_real_function_3d F12Potentials::apply_regularized_potential(const real_function_3d& ket1, const real_function_3d& ket2,const real_function_3d Ki, const real_function_3d& Kj, const vector_real_function_3d& bra,const vector_real_function_3d& Kpno) const {
  //const double muleps = 0.0;//FunctionDefaults < 3 > ::get_thresh();

	MyTimer timeU = MyTimer(world).start();
	const vector_real_function_3d Uepart = convolve_with_U(bra, ket1, ket2);
	timeU.stop().print("Vreg:Ue");

	MyTimer timeK2 = MyTimer(world).start();
	vector_real_function_3d comm = zero_functions<double,3>(world,Uepart.size());
	if(param.exchange()=="full"){
		// K1 part
		vector_real_function_3d Kv = Kpno;
		if(Kv.empty()){
			MyTimer tk=MyTimer(world).start();
			Kv=K(bra);
			tk.stop().print("recompute Kv");
		}
		const vector_real_function_3d Kvi = Kv*ket1;
		const vector_real_function_3d vKi = bra*Ki;
		vector_real_function_3d difftmp=Kvi-vKi;
		//truncate(world,difftmp);
		const vector_real_function_3d KffK1 = apply(world,*fop,difftmp)*ket2; // saves application of 1 f12 convolution

		// K2 part
		const vector_real_function_3d vi = bra*ket1;
		const vector_real_function_3d fvi = apply(world,*fop,vi);
		const vector_real_function_3d KffK2 = K(fvi*ket2)-fvi*Kj; // not faster

		comm=KffK1+KffK2;
	}else if (world.rank()==0){
		std::cout << "exchange is " << param.exchange() << " neglecting [K,f] part\n";
	}
	timeK2.stop().print("Vreg:K2");
	//	const double err = norm2(world,((K1f - fK1 + K2f - fK2) - test));
	//	std::cout << "error=" << err << "\n";
	vector_real_function_3d result = Uepart - comm;//(K1f - fK1 + K2f - fK2);
	truncate(world,result);


	return result;

	MADNESS_EXCEPTION("should not end up here", 1);
}

/// Ue part is not really faster than with the fluctuation potential but saves memory
/// [K,f] part is faster as over the fluctuation potential
Tensor<double> F12Potentials::compute_regularized_fluctuation_matrix(const ElectronPairIterator& it, const vector_real_function_3d& pnos, const vector_real_function_3d& Kpnos)const{
	const real_function_3d& Ki=acKmos[it.i()];
	const real_function_3d& Kj=acKmos[it.j()];
	const real_function_3d& moi=acmos[it.i()];
	const real_function_3d& moj=acmos[it.j()];


	MyTimer timeIm = MyTimer(world).start();
	const vector_real_function_3d vj = mul(world,moj,pnos,false);
	vector_real_function_3d vi;
	if(it.diagonal()) vi=vj;
	else vi = mul(world,moi,pnos,false);
	const vector_real_function_3d di=grad(moi); // move to intermediates
	const vector_real_function_3d dj=grad(moj); // move to intermediates
	world.gop.fence();
	timeIm.stop().print("Wij: Intermediates");

	Tensor<double> Uepart(pnos.size(),pnos.size());
	MyTimer timeU = MyTimer(world).start();
	if(false){
		// use this untill bug is fixed
		Uepart=matrix_inner(world,pnos,convolve_with_U(pnos,moi,moj),it.diagonal());
		//Uepart = 0.5*(matrix_inner(world,pnos,convolve_with_U(pnos,moi,moj),it.diagonal()) + matrix_inner(world,pnos,convolve_with_U(pnos,moj,moi),it.diagonal()));
	}else{
		// for the local part there is currently no cheaper option
		const vector_real_function_3d ul = convolve_with_local_U(vi);
		Tensor<double> Uel=matrix_inner(world,vj,ul);

		if(it.diagonal()){
			const vector_real_function_3d di=grad(moi); // move to intermediates

			Tensor<double> Uenl(pnos.size(),pnos.size());
			for(size_t axis=0;axis<3;++axis){
				const vector_real_function_3d ket=convolve_with_gradslater(vj,axis);
				const vector_real_function_3d bra=mul(world,di[axis],pnos,false);
				Uenl += matrix_inner(world,bra,ket);
			}

			const double prefactor = 1.0 / (2.0 * (param.gamma()));
			Uenl = prefactor*(Uenl+transpose(Uenl));

			Uepart = Uel + Uenl;

		}else{

			MyTimer timer_nl=MyTimer(world).start();
			Tensor<double> Uenl1(pnos.size(),pnos.size());
			Tensor<double> Uenl2(pnos.size(),pnos.size());
			for(size_t axis=0;axis<3;++axis){
				const vector_real_function_3d ket=convolve_with_gradslater(vj,axis);
				const vector_real_function_3d bra=mul(world,di[axis],pnos);
				Uenl1 += matrix_inner(world,bra,ket);
				const vector_real_function_3d ket2=convolve_with_gradslater(vi,axis);
				const vector_real_function_3d bra2=mul(world,dj[axis],pnos);
				Uenl2 += matrix_inner(world,bra2,ket2);
			}
			timer_nl.stop().print("Wij: Ueij nonloc");

			const double prefactor = 1.0 / (2.0 * (param.gamma()));
			Tensor<double> Uenl = prefactor*(transpose(Uenl1)+Uenl2);

			Uepart = Uel + Uenl;
		}
	}
	timeU.stop().print("Wij: Ueij");


	Tensor<double> Kpart(pnos.size(),pnos.size());
	MyTimer timeK = MyTimer(world).start();
	{
		if(it.diagonal()){
			const vector_real_function_3d vKi = mul(world,Ki,pnos,false);
			const vector_real_function_3d iKv = mul(world,moi,Kpnos,false);
			world.gop.fence();
			// K1 part
			vector_real_function_3d ket1=apply(world,*fop,(iKv-vKi));
			Tensor<double> K1 = matrix_inner(world,vj,ket1);
			Kpart = K1+transpose(K1);

		}else{
			const vector_real_function_3d vKi = mul(world,Ki,pnos,false);
			const vector_real_function_3d iKv = mul(world,moi,Kpnos,false);
			const vector_real_function_3d vKj = mul(world,Kj,pnos,false);
			const vector_real_function_3d jKv = mul(world,moj,Kpnos,false);
			world.gop.fence();
			// K1 part
			vector_real_function_3d ket1=apply(world,*fop,(iKv-vKi));
			Tensor<double> K1 = matrix_inner(world,vj,ket1);
			ket1.clear();
			// K2 part
			vector_real_function_3d ket2=apply(world,*fop,(jKv-vKj));
			Tensor<double> K2 = matrix_inner(world,vi,ket2);
			ket2.clear();
			Kpart = K1+transpose(K2);
		}
	}
	timeK.stop().print("Wij: Kij");
	return transpose(Uepart-Kpart); // we computed <ba|ij> instead of <ab|ij> (historical reasons)
}
Tensor<double> F12Potentials::compute_regularized_fluctuation_matrix(const std::pair<vector_real_function_3d,vector_real_function_3d> & KPNOA, const std::pair<vector_real_function_3d,vector_real_function_3d> & KPNOB, const std::pair<real_function_3d,real_function_3d>& MKI, const std::pair<real_function_3d,real_function_3d> MKJ)const{

	const vector_real_function_3d& a=KPNOA.first;
	const vector_real_function_3d& Ka=KPNOA.second;
	bool same_pnos=false;
	vector_real_function_3d b=KPNOB.first;
	vector_real_function_3d Kb=KPNOB.second;
	if(KPNOB.first.empty()){
		same_pnos=true;
		b=a;
		Kb=Ka;
	}
	const real_function_3d& moi=MKI.first;
	const real_function_3d& Ki=MKI.second;
	real_function_3d moj=MKJ.first;
	real_function_3d Kj=MKJ.second;
	bool diagonal=false;
	if(MKJ.first.is_initialized()==false){
		diagonal=true;
		moj=moi;
		Kj=Ki;
	}

	MyTimer timeIm = MyTimer(world).start();
	const vector_real_function_3d bj = mul(world,moj,b,false);
	vector_real_function_3d ai;
	if(diagonal and same_pnos) ai=bj;
	else ai = mul(world,moi,a,false);
	const vector_real_function_3d di=grad(moi);
	vector_real_function_3d dj;
	if(diagonal and same_pnos) dj=di;
	else dj=grad(moj);
	world.gop.fence();
	timeIm.stop().print("Wij: Intermediates");

	Tensor<double> Uepart(a.size(),b.size());
	MyTimer timeU = MyTimer(world).start();
	{
		// for the local part there is currently no cheaper option
		const vector_real_function_3d ul = convolve_with_local_U(ai);
		Tensor<double> Uel=matrix_inner(world,bj,ul);

		if(diagonal and same_pnos){
			const vector_real_function_3d& pnos=a;
			//const vector_real_function_3d& vi=ai;
			const vector_real_function_3d& vj=bj;
			Tensor<double> Uenl(pnos.size(),pnos.size());
			for(size_t axis=0;axis<3;++axis){
				const vector_real_function_3d ket=convolve_with_gradslater(vj,axis);
				const vector_real_function_3d bra=mul(world,di[axis],pnos,false);
				Uenl += matrix_inner(world,bra,ket);
			}

			const double prefactor = 1.0 / (2.0 * (param.gamma()));
			Uenl = prefactor*(Uenl+transpose(Uenl));

			Uepart = Uel + Uenl;

		}else{

			MyTimer timer_nl=MyTimer(world).start();
			Tensor<double> Uenl1(a.size(),b.size());
			Tensor<double> Uenl2(b.size(),a.size());
			for(size_t axis=0;axis<3;++axis){
				const vector_real_function_3d ket=convolve_with_gradslater(bj,axis);
				const vector_real_function_3d bra=mul(world,di[axis],a);
				Uenl1 += matrix_inner(world,bra,ket);
				const vector_real_function_3d ket2=convolve_with_gradslater(ai,axis);
				const vector_real_function_3d bra2=mul(world,dj[axis],b);
				Uenl2 += matrix_inner(world,bra2,ket2);
			}
			timer_nl.stop().print("Wij: Ueij nonloc");

			const double prefactor = 1.0 / (2.0 * (param.gamma()));
			Tensor<double> Uenl = prefactor*(transpose(Uenl1)+Uenl2);

			Uepart = Uel + Uenl;
		}
	}
	timeU.stop().print("Wij: Ueij");

	Tensor<double> Kpart(a.size(),b.size());
	MyTimer timeK = MyTimer(world).start();
	{
		if(diagonal and same_pnos){
			const vector_real_function_3d aKi = mul(world,Ki,a,false);
			const vector_real_function_3d iKa = mul(world,moi,Ka,false);
			world.gop.fence();
			// K1 part
			vector_real_function_3d ket1=apply(world,*fop,(iKa-aKi))*(-1.0/(2.0*param.gamma()));
			Tensor<double> K1 = matrix_inner(world,bj,ket1);
			Kpart = K1+transpose(K1);

		}else{
			const vector_real_function_3d aKi = mul(world,Ki,a,false);
			const vector_real_function_3d iKa = mul(world,moi,Ka,false);
			const vector_real_function_3d bKj = mul(world,Kj,b,false);
			const vector_real_function_3d jKb = mul(world,moj,Kb,false);
			world.gop.fence();
			// K1 part
			vector_real_function_3d ket1=apply(world,*fop,(iKa-aKi));
			Tensor<double> K1 = matrix_inner(world,bj,ket1);
			ket1.clear();
			// K2 part
			vector_real_function_3d ket2=apply(world,*fop,(jKb-bKj));
			Tensor<double> K2 = matrix_inner(world,ai,ket2);
			ket2.clear();
			Kpart = K1+transpose(K2);
		}
	}
	timeK.stop().print("Wij: Kij");
	return transpose(Uepart-Kpart); // this whole routine computed the <ba|op|ij> matrix
}

std::valarray<Tensor<double> > F12Potentials::compute_ijab_integrals(const std::valarray<vector_real_function_3d>& functions, const std::shared_ptr<real_convolution_3d>& op) const {
	if(functions.size()!=npairs()) return std::valarray<Tensor<double> >();
	std::valarray < Tensor<double> > result(npairs());
	for (ElectronPairIterator it=pit(); it; ++it) {
		result[it.ij()] = compute_xyab_integrals(acmos[it.i()], acmos[it.j()], functions[it.ij()], functions[it.ij()], op, it.diagonal());
	}
	return result;
}

std::valarray<double> F12Potentials::compute_fQg_integrals() const {
	// first npair are the <ij|fQg|ij> second the <ij|fQg|ji>
	std::valarray<double> result(2 * npairs());
	const double muleps = 0.0;
	for (ElectronPairIterator it=pit(); it; ++it) {
		const size_t i = it.i();
		const size_t j = it.j();
		const size_t ij = it.ij();
		const size_t ji = it.ij() + it.npairs();
		vector_real_function_3d fim = apply(world, *fop, mul_sparse(world, acmos[i], mos, muleps)); // all mos !
		truncate(world, fim);
		vector_real_function_3d gim = apply(world, *coulombop, mul_sparse(world, acmos[i], mos, muleps)); // all mos !
		truncate(world, gim);
		vector_real_function_3d jfim = mul_sparse(world, acmos[j], fim, muleps);
		vector_real_function_3d jgim = mul_sparse(world, acmos[j], gim, muleps);
		truncate(world, jfim);
		truncate(world, jgim);
		const double ijfO1gij = madness::inner(jfim, jgim);
		Tensor<double> ijfmn = matrix_inner(world, jfim, mos); // all mos !
		Tensor<double> ijgmn = matrix_inner(world, jgim, mos); // all mos !
		const double ijfO12gij = inner(ijfmn, ijgmn);
		if (it.diagonal()) {
			const real_function_3d ii = (acmos[i] * acmos[i]).truncate();
			const double ijfgij = madness::inner(ii, convolve_with_fg(ii));
			const double result_ij = ijfgij - 2.0 * ijfO1gij + ijfO12gij;
			//std::cout << "\npart1=" << ijfgij << "\npart2=" <<(ijfO1gij-0.5*ijfO12gij) << "\npart3=" << (ijfO1gij-0.5*ijfO12gij) << "\ndbpart=" <<ijfO1gij  << "\ndbpart2=" << 0.5*ijfO12gij <<  "\n";
			result[ij] = result_ij;
			result[ji] = result_ij;
		} else {
			const real_function_3d ii = (acmos[i] * acmos[i]).truncate();
			const real_function_3d jj = (acmos[j] * acmos[j]).truncate();
			const double ijfgij = madness::inner(ii, convolve_with_fg(jj));
			vector_real_function_3d fjm = apply(world, *fop, mul_sparse(world, acmos[j], mos,  muleps));// all mos !
			vector_real_function_3d gjm = apply(world, *coulombop, mul_sparse(world, acmos[j], mos,  muleps));// all mos !
			vector_real_function_3d ifjm = mul_sparse(world, acmos[i], fjm,  muleps);
			vector_real_function_3d igjm = mul_sparse(world, acmos[i], gjm,  muleps);
			const double ijfO2gij = madness::inner(ifjm, igjm);
			const double result_ij = ijfgij - ijfO1gij - ijfO2gij + ijfO12gij;
			result[ij] = result_ij;
			const real_function_3d mos_ij = (acmos[i] * acmos[j]).truncate();
			const double ijfgji = madness::inner(mos_ij, convolve_with_fg(mos_ij));
			const double ijfO1gji = madness::inner(jfim, igjm);
			double ijfO2gji = madness::inner(ifjm, jgim);
			const double ijfO12gji = inner(ijfmn, transpose(ijgmn));
			double result_ji = ijfgji - ijfO1gji - ijfO2gji + ijfO12gji;
			result[ji] = result_ji;
		}
	}
	return result;
}

std::valarray<double> F12Potentials::compute_fQU_integrals() const {
	// first npair are the <ij|fQg|ij> second the <ij|fQg|ji>
	std::valarray<double> result(2.0 * npairs());
	const double muleps=0.0;
	for (ElectronPairIterator it=pit(); it; ++it) {
		const size_t i = it.i();
		const size_t j = it.j();
		const size_t ij = it.ij();
		const size_t ji = it.ij() + it.npairs();
		vector_real_function_3d fim = apply(world, *fop, mul_sparse(world, acmos[i], mos, muleps));
		truncate(world, fim);
		vector_real_function_3d jfim = mul_sparse(world, acmos[j], fim, muleps);
		vector_real_function_3d jUim = convolve_with_U(mos, acmos[i], acmos[j]);
		truncate(world, jfim);
		truncate(world, jUim);
		const double ijfO1Uij = madness::inner(jfim, jUim);
		Tensor<double> ijfmn = matrix_inner(world, jfim, mos);
		Tensor<double> ijUmn = matrix_inner(world, jUim, mos);
		const double ijfO12Uij = inner(ijfmn, ijUmn);
		const double ijfUij = madness::inner(acmos[j],convolve_with_fU(acmos[i], acmos[i], acmos[j]));
		if (it.diagonal()) {
			const double result_ij = ijfUij - 2.0 * ijfO1Uij + ijfO12Uij;
			result[ij] = result_ij;
			result[ji] = result_ij;
			if(param.debug() and world.rank()==0){
				std::cout << "fQU debug for pair " << it.ij() << "=(" << it.i() << "," << it.j() << ")\n";
				std::cout << "ijfUij   =" << ijfUij  <<"\n";
				std::cout << "ijfO1Uij =" << ijfO1Uij <<"\n";
				std::cout << "ijfO12Uij=" << ijfO12Uij<<"\n";
			}
		} else {
			vector_real_function_3d fjm = apply(world, *fop, mul_sparse(world, acmos[j], mos, muleps));
			vector_real_function_3d ifjm = mul_sparse(world, acmos[i], fjm, muleps);
			vector_real_function_3d iUjm = convolve_with_U(mos, acmos[j], acmos[i]);
			const double ijfO2Uij = madness::inner(ifjm, iUjm);
			const double result_ij = ijfUij - ijfO1Uij - ijfO2Uij + ijfO12Uij;
			result[ij] = result_ij;
			const double ijfUji = madness::inner(acmos[j], convolve_with_fU(acmos[i], acmos[j], acmos[i]));
			const double ijfO1Uji = madness::inner(jfim, iUjm);
			const double ijfO2Uji = madness::inner(ifjm, jUim);
			const double ijfO12Uji = inner(ijfmn, transpose(ijUmn));
			double result_ji = ijfUji - ijfO1Uji - ijfO2Uji + ijfO12Uji;
			result[ji] = result_ji;
			if(param.debug() and world.rank()==0){
				std::cout << "fQU debug for pair " << it.ij() << "=(" << it.i() << "," << it.j() << ")\n";
				std::cout << "ijfUij   =" << ijfUij   <<"\n";
				std::cout << "ijfO1Uij =" << ijfO1Uij <<"\n";
				std::cout << "ijfO2Uij =" << ijfO2Uij <<"\n";
				std::cout << "ijfO12Uij=" << ijfO12Uij<<"\n";
			}
		}
	}
	return result;
}

/// Compute the <ab|[f,K]|ij> integrals
std::valarray<Tensor<double> > F12Potentials::compute_cijab_integrals(const std::valarray<vector_real_function_3d>& functions) const {

	if(functions.size()!=npairs()) return std::valarray<Tensor<double> >();

	//vector_real_function_3d Kmos = K(mos);

	std::valarray < Tensor<double> > result(2.0*npairs());

	for (ElectronPairIterator it=pit(); it; ++it) {
		result[it.ij()]=compute_cijab_integrals(it.i(),it.j(), acKmos[it.i()], acKmos[it.j()],functions[it.ij()]);
		result[it.ij()+npairs()]=compute_cijab_integrals(it.j(),it.i(), acKmos[it.j()], acKmos[it.i()],functions[it.ij()]);
	}
	return result;
}
Tensor<double> F12Potentials::compute_cijab_integrals(const size_t &i, const size_t& j, const real_function_3d& Ki, const real_function_3d& Kj , const vector_real_function_3d& functions)const{
	const Tensor<double> fK1 = compute_xyab_integrals(Ki, acmos[j], functions, functions, fop, false); // never symmetric since Ki!=i
	const Tensor<double> fK2 = compute_xyab_integrals(acmos[i], Kj, functions, functions, fop, false);

	const vector_real_function_3d Kv = K(functions);

	const Tensor<double> K1f = compute_xyab_integrals(acmos[i], acmos[j], Kv, functions, fop);
	const Tensor<double> K2f = compute_xyab_integrals(acmos[i], acmos[j], functions, Kv, fop);

	const Tensor<double> result_ij =  K1f + K2f - fK1 - fK2;

	return result_ij;
}

real_function_3d F12Potentials::convolve_with_U(const real_function_3d& bra, const real_function_3d& ket1, const real_function_3d& ket2, const bool symmetric) const {
	const real_function_3d local_part = convolve_with_local_U((bra * ket1).truncate()) * ket2;
	const real_function_3d nonlocal_part = convolve_with_nonlocal_U(bra, ket1, ket2, symmetric);
	const real_function_3d result = (local_part + nonlocal_part).truncate();
	return result;
}
vector_real_function_3d F12Potentials::convolve_with_U(const vector_real_function_3d& bra, const real_function_3d& ket1, const real_function_3d& ket2, const bool symmetric) const {
	const double muleps=0.0;
	MyTimer time_local=MyTimer(world).start();
	const vector_real_function_3d braket1 = mul_sparse(world, ket1, bra, muleps);
	const vector_real_function_3d ul = convolve_with_local_U(braket1);
	const vector_real_function_3d local_part = mul_sparse(world, ket2, ul, muleps);
	time_local.stop().print("Ue: local");

	MyTimer time_nonlocal=MyTimer(world).start();
	const vector_real_function_3d nonlocal_part = convolve_with_nonlocal_U(bra, ket1, ket2, symmetric);
	time_nonlocal.stop().print("Ue: non-local");
	vector_real_function_3d result = (local_part + nonlocal_part);
	truncate(world, result);
	return result;
}

vector_real_function_3d F12Potentials::convolve_with_local_U(const vector_real_function_3d& functions) const {
	const double factor1 = 2.0 * param.gamma();
	const vector_real_function_3d part1 = factor1 * convolve_with_fg(functions);
	const double factor2 = param.gamma() / 2.0;
	const vector_real_function_3d part2 = factor2 * convolve_with_slater_potential(functions);
	vector_real_function_3d result = (part1 + part2);
	truncate(world, result);
	return result;
}

/// @return \f$ \int dx2 f12*g12*f(2)  \f$
/// @param[in] function: function on which f*g operator is applied
/// right now we use: fg = 1/(2*gamma)*(g - 4*pi*G(gamma))
/// todo: implement fg operator (may be cheaper)
real_function_3d F12Potentials::convolve_with_fg(const real_function_3d& function) const {
	const double prefactor = 1.0 / (2.0 * param.gamma());
	const double bshfactor = 4.0 * constants::pi;
	const real_function_3d result = prefactor * ((*coulombop)(function) - bshfactor * (*bshop)(function)).truncate();
	return result;
}
vector_real_function_3d F12Potentials::convolve_with_fg(const vector_real_function_3d& functions) const {
	const double prefactor = 1.0 / (2.0 * param.gamma());
	const double bshfactor = 4.0 * constants::pi;
	const vector_real_function_3d gpart = apply(world, *coulombop, functions);
	const vector_real_function_3d bshpart = apply(world, *bshop, functions);
	vector_real_function_3d result = prefactor *( gpart - bshfactor * bshpart);
	truncate(world, result);
	return result;
}

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
real_function_3d F12Potentials::convolve_with_nonlocal_U(const real_function_3d& bra, const real_function_3d& ket1, const real_function_3d& ket2, const bool symmetric, const bool& squared) const {
	vector_real_function_3d dummy(1, bra);
	return convolve_with_nonlocal_U(dummy, ket1, ket2, symmetric, squared).front();
}

vector_real_function_3d F12Potentials::convolve_with_nonlocal_U(const vector_real_function_3d& bra, const real_function_3d& ket1, const real_function_3d& ket2, const bool symmetric, const bool& squared) const {
	const double muleps=0.0;
	double factor = 1.0;
	if (squared)
		factor = 2.0;
	const double prefactor = 1.0 / (2.0 * (param.gamma() * factor));

	const std::vector<std::shared_ptr<Derivative<double, 3> > >& D = gradop;

	// make gradient of ket1 and ket2
	std::vector < real_function_3d > Dket1 = zero_functions<double, 3>(world, 3);
	std::vector < real_function_3d > Dket2 = zero_functions<double, 3>(world, 3);
	for (size_t axis = 0; axis < 3; ++axis) {
		Dket1[axis] = (*D[axis])(ket1);
		if (symmetric)
			Dket2[axis] = Dket1[axis];
		else
			Dket2[axis] = (*D[axis])(ket2);
	}

	// part1: <bra|Dsop|Dket1>*|ket2>
	vector_real_function_3d part1 = zero_functions<double, 3>(world, bra.size());
	for (size_t axis = 0; axis < 3; ++axis) {
		vector_real_function_3d braDket1 = mul_sparse(world, Dket1[axis], bra,  muleps);
		vector_real_function_3d tmp1 = convolve_with_gradslater(braDket1, axis, squared);
		vector_real_function_3d tmp2 = mul_sparse(world, ket2, tmp1, muleps);
		part1 += tmp2;
	}

	// part2: <bra|Dsop|ket1>*|Dket2>
	vector_real_function_3d part2 = zero_functions<double, 3>(world, bra.size());
	vector_real_function_3d braket1 = mul_sparse(world, ket1, bra,  muleps);
	for (size_t axis = 0; axis < 3; ++axis) {
		vector_real_function_3d tmp1 = convolve_with_gradslater(braket1, axis, squared);
		vector_real_function_3d tmp2 = mul_sparse(world, Dket2[axis], tmp1,  muleps);
		part2 += tmp2;
	}

	vector_real_function_3d result = prefactor * (part1 - part2);
	truncate(world, result);
	return (-1.0) * result; // check sign -> correct since it should be -part1+part2
}

vector_real_function_3d F12Potentials::convolve_with_fU(const vector_real_function_3d& bra, const real_function_3d& ket1, const real_function_3d& ket2, const bool symmetric) const {
	const double prefactor = 1.0 / (2.0 * param.gamma());
	const double muleps=0.0;
	vector_real_function_3d result = zero_functions<double, 3>(world, bra.size());
	// local part
	{
		const vector_real_function_3d ket1bras = mul_sparse(world, ket1, bra,  muleps);
		const double bshfactor = 4.0 * constants::pi;
		const vector_real_function_3d part1 = convolve_with_local_U(ket1bras);
		const vector_real_function_3d part2 = apply(world, *bshop, ket1bras);
		const vector_real_function_3d part3 = apply(world, *bshop_sq, ket1bras);
		const vector_real_function_3d part4 = apply(world, *slaterop_sq, ket1bras);
		const vector_real_function_3d local_part = part1 - bshfactor * part2 + bshfactor * part3 - param.gamma() / 2.0 * part4;
		result += mul_sparse(world, ket2, local_part,  muleps);
		if(param.debug()){
			const double db1 = norm2(world,part1);
			const double db2 = norm2(world,part2);
			const double db3 = norm2(world,part3);
			const double db4 = norm2(world,part4);
			if(world.rank()==0) std::cout << "convolve_with_fU norm of local parts:\n" << db1 << ", " << db2 << ", " << db3 << ", " << db4 << "\n";
		}
	}
	// nonlocal part
	{
		const vector_real_function_3d part1 = convolve_with_nonlocal_U(bra, ket1, ket2, symmetric, false);
		const vector_real_function_3d part2 = convolve_with_nonlocal_U(bra, ket1, ket2, symmetric, true);
		const vector_real_function_3d nonlocal_part = (part1 - part2);
		result += nonlocal_part;
		if(param.debug()){
			const double db1 = norm2(world,part1);
			const double db2 = norm2(world,part2);
			if(world.rank()==0) std::cout << "convolve_with_fU norm of nonlocal parts:\n" << db1 << ", " << db2 << "\n";
		}
	}
	scale(world, result, prefactor);
	truncate(world, result);
	return result;
}

std::valarray<double> F12Potentials::compute_fQc_integrals(
		const vector_real_function_3d& Kmos,
		const std::valarray<vector_real_function_3d>& functions) const {
	//if(functions.size()!=npairs()) return std::valarray<double>();
	std::valarray<double> result(2.0 * npairs());
	for (ElectronPairIterator it = pit(); it; ++it) {
		std::pair<double, double> rij = compute_fQc_integrals_ij(Kmos,
				functions[it.ij()], it);
		result[it.ij()] = rij.first;
		result[it.ij() + it.npairs()] = rij.second;
	}
	return result;
}

std::pair<double, double> F12Potentials::compute_fQc_integrals_ij(
		const vector_real_function_3d& Kmos,
		const vector_real_function_3d& functions,
		const ElectronPairIterator& it,
		const bool& use_no_intermediates) const {
	const vector_real_function_3d& v = functions;
	const real_function_3d& moi = acmos[it.i()];
	const real_function_3d& moj = acmos[it.j()];
	const real_function_3d& Ki = acKmos[it.i()];
	const real_function_3d& Kj = acKmos[it.j()];
	if (use_no_intermediates) {
		// alternative with flexible: Q12=OxQ or QOx -> more accurate this way
		// additionally no large intermediates are calculated
		// in the long run it is probably better to use intermediates but restricted to small blocks
		// so just feed the virtuals in chunks
		// keeping this code to debug and to compare performance
		if (world.rank() == 0)
			std::cout
			<< "!!!fQc is calculated without intermediates! Is this by choice?\n";

		// Part 1 K1, Q12=OxQ, ij part
		// <j|(ifx)*Q*((Kxfi)-(xfKi))|j>
		// ji part
		// <i|(jfx)*Q*((Kxfi)-(xfKi))|j> or  <j|(ifx)*Q*((Kxfj)-(xfKj))|i>
		// Part 2 K2, Q12=QOx, ij part
		// <i|(jfx)*Q*((Kxfj)-(xfKj))|i>
		// ji part
		// <j|(ifx)*Q*((Kxfj)-(xfKj))|i> or <i|(jfx)*Q*((Kxfi)-(xfKi))|j>
		double result1_ij = 0.0;
		double result1_ji = 0.0;
		double result2_ij = 0.0;
		double result2_ji = 0.0;
		for (const auto& x : v) {
			const real_function_3d braijx = Q((*fop)(moi * x) * moj);
			const real_function_3d ket1 = ((*fop)((K(x) * moi - x * Ki)) * moj);
			const double ijpart = madness::inner(braijx, ket1);
			result1_ij += ijpart;
			if (it.diagonal()) {
				// in that case all parts are the same
				result2_ij += ijpart;
				result1_ji += ijpart;
				result2_ji += ijpart;
			} else {
				const real_function_3d brajix = Q((*fop)(moj * x) * moi);
				result1_ji += madness::inner(brajix, ket1);
				const real_function_3d ket2 = ((*fop)((K(x) * moj - x * Kj))
						* moi);
				result2_ij += madness::inner(brajix, ket2); //braji is correct
				result2_ji += madness::inner(braijx, ket2); //braij is correct
			}
		}
		if (param.debug()) {
			if (world.rank() == 0)
				std::cout << "fQc_1_ij=" << result1_ij << "\n" << "fQc_1_ji="
				<< result1_ji << "\n" << "fQc_2_ij=" << result2_ij
				<< "\n" << "fQc_2_ji=" << result2_ji << "\n";
		}
		if (it.diagonal())
			MADNESS_ASSERT(
					result1_ij == result2_ij && result1_ij == result1_ji
					&& result1_ji == result2_ji);

		const double resultij = result1_ij + result2_ij;
		const double resultji = result1_ji + result2_ji;
		return std::make_pair(resultij, resultji);
	} else {
		// like before but with intermediates, less fences, should be faster
		// calulate intermediates
		const vector_real_function_3d Kv = K(v);
		{
			const double tmpsize = get_size(world, Kv);
			if (world.rank() == 0)
				std::cout << "Kv intermediate: " << tmpsize << " Gbyte\n";
		}
		const vector_real_function_3d braij = Q(
				apply(world, *fop, moi * v) * moj);
		{
			const double tmpsize = get_size(world, braij);
			if (world.rank() == 0)
				std::cout << "Q(braij) intermediate: " << tmpsize << " Gbyte\n";
		}            // Part 1 K1, Q12=OxQ, ij part
		// <j|(ifx)*Q*((Kxfi)-(xfKi))|j>
		// ji part
		// <i|(jfx)*Q*((Kxfi)-(xfKi))|j> or  <j|(ifx)*Q*((Kxfj)-(xfKj))|i>
		// Part 2 K2, Q12=QOx, ij part
		// <i|(jfx)*Q*((Kxfj)-(xfKj))|i>
		// ji part
		// <j|(ifx)*Q*((Kxfj)-(xfKj))|i> or <i|(jfx)*Q*((Kxfi)-(xfKi))|j>
		double result1_ij = 0.0;
		double result1_ji = 0.0;
		double result2_ij = 0.0;
		double result2_ji = 0.0;
		if (it.diagonal()) {
			const vector_real_function_3d Kvi = mul(world, moi, Kv, false);
			const vector_real_function_3d vKi = mul(world, Ki, v, false);
			world.gop.fence();
			const double tmpsize1 = get_size(world, braij);
			if (world.rank() == 0)
				std::cout << "Kvi intermediate: " << tmpsize1 << " Gbyte\n";

			const double tmpsize2 = get_size(world, braij);
			if (world.rank() == 0)
				std::cout << "vKi intermediate: " << tmpsize2 << " Gbyte\n";

			const vector_real_function_3d ket1 = (apply(world, *fop, Kvi - vKi))
											* moj;
			result1_ij = madness::inner(braij, ket1);
			result2_ij = result1_ij;
			result1_ji = result1_ij;
			result2_ji = result1_ij;
		} else {
			const vector_real_function_3d braji = Q(
					apply(world, *fop, moj * v) * moi);
			{
				const double tmpsize = get_size(world, braij);
				if (world.rank() == 0)
					std::cout << "Q(braij) intermediate: " << tmpsize
					<< " Gbyte\n";
			}
			const vector_real_function_3d Kvi = mul(world, moi, Kv, false);
			const vector_real_function_3d vKi = mul(world, Ki, v, false);
			const vector_real_function_3d Kvj = mul(world, moj, Kv, false);
			const vector_real_function_3d vKj = mul(world, Kj, v, false);
			world.gop.fence(); //important fence
			const double tmpsize1 = get_size(world, Kvi);
			if (world.rank() == 0)
				std::cout << "Kvi intermediate: " << tmpsize1 << " Gbyte\n";

			const double tmpsize2 = get_size(world, vKi);
			if (world.rank() == 0)
				std::cout << "vKi intermediate: " << tmpsize2 << " Gbyte\n";

			const double tmpsize3 = get_size(world, Kvj);
			if (world.rank() == 0)
				std::cout << "Kvj intermediate: " << tmpsize3 << " Gbyte\n";

			const double tmpsize4 = get_size(world, vKj);
			if (world.rank() == 0)
				std::cout << "vKj intermediate: " << tmpsize4 << " Gbyte\n";

			const vector_real_function_3d ket1 = (apply(world, *fop, Kvi - vKi))
											* moj;
			const vector_real_function_3d ket2 = (apply(world, *fop, Kvj - vKj))
											* moi;
			result1_ij = madness::inner(braij, ket1);
			result2_ij = madness::inner(braji, ket2); // braji is correct
			result1_ji = madness::inner(braji, ket1);
			result2_ji = madness::inner(braij, ket2); // braij is correct
		}
		if (param.debug()) {
			if (world.rank() == 0)
				std::cout << "fQc_1_ij=" << result1_ij << "\n" << "fQc_1_ji="
				<< result1_ji << "\n" << "fQc_2_ij=" << result2_ij
				<< "\n" << "fQc_2_ji=" << result2_ji << "\n";
		}
		if (it.diagonal())
			MADNESS_ASSERT(
					result1_ij == result2_ij && result1_ij == result1_ji
					&& result1_ji == result2_ji);

		const double resultij = result1_ij + result2_ij;
		const double resultji = result1_ji + result2_ji;
		return std::make_pair(resultij, resultji);
	}
}

void F12Potentials::print_pair_energies(const PNOPairs& pairs,
		const std::string& msg) const {
	std::cout << std::scientific << std::setprecision(5) << std::setfill(' ');
	print_pair_energies(pairs.energies.eijs, pairs.energies.eijt, msg,
			pairs.type);
	if (world.rank() == 0 && pairs.type == CISPD_PAIRTYPE) {
		std::cout << "-------------------------------------------\n";
		std::cout << "GS correction reg=" << std::setw(12)
		<< pairs.energies.eijs.sum() << "\n";
		std::cout << "GS correction f12=" << std::setw(12)
		<< pairs.energies.eijs_f12.sum() << "\n";
		std::cout << "GS correction tot=" << std::setw(12)
		<< pairs.energies.eijs.sum() + pairs.energies.eijs_f12.sum()
		<< "\n";
		std::cout << "ES correction reg=" << std::setw(12)
		<< pairs.energies.eijt.sum() << "\n";
		std::cout << "ES correction f12=" << std::setw(12)
		<< pairs.energies.eijt_f12.sum() << "\n";
		std::cout << "ES correction tot=" << std::setw(12)
		<< pairs.energies.eijt.sum() + pairs.energies.eijt_f12.sum()
		<< "\n";
		std::cout << "-------------------------------------------\n";
		std::cout << "∆CIS(D) reg      =" << std::setw(12)
		<< pairs.energies.energy << "\n";
		std::cout << "∆CIS(D) f12      =" << std::setw(12)
		<< pairs.energies.energy_f12 << "\n";
		;
		std::cout << "∆CIS(D) tot      =" << std::setw(12)
		<< pairs.energies.total_energy() << "\n";
		;
		std::cout << "-------------------------------------------\n";
	} else if (world.rank() == 0) {
		std::cout << "-------------------------------------------\n";
		std::cout << std::scientific << std::setprecision(5);
		std::cout << "Total Reg Energy  =" << std::setw(12)
		<< pairs.energies.energy << "\n";
		std::cout << "Total F12 Energy  =" << std::setw(12)
		<< pairs.energies.energy_f12 << "\n";
		;
		std::cout << "Sum Energies Total=" << std::setw(12)
		<< pairs.energies.total_energy() << "\n";
		;
		std::cout << "-------------------------------------------\n";
	}
}

double F12Potentials::compute_fQg_integral(const real_function_3d bra1,
		const real_function_3d& bra2, const real_function_3d& ket1,
		const real_function_3d& ket2, const bool& diagonal) const {
	Projector<double, 3> O(mos);
	const double part1 = madness::inner(bra1 * ket1,
			convolve_with_fg(bra2 * ket2));
	double part2 = 0.0;
	double part3 = 0.0;
	{
		vector_real_function_3d O1gij;
		for (const auto& k : mos)
			O1gij.push_back((*coulombop)(k * ket1) * ket2);
		vector_real_function_3d O1q2gij = O1gij - 0.5 * O(O1gij); // this is particle 2, particle 1 is just acmos
		const vector_real_function_3d& p2 = O1q2gij;
		const vector_real_function_3d& p1 = mos;
		MADNESS_ASSERT(p1.size() == p2.size());
		// convolute
		part2 =
				madness::inner(world, bra2 * p2, apply(world, *fop, bra1 * p1)).sum();
	}
	if (!diagonal) {
		vector_real_function_3d O2gij;
		for (const auto& k : mos)
			O2gij.push_back((*coulombop)(k * ket2) * ket1);
		vector_real_function_3d q1O2gij = O2gij - 0.5 * O(O2gij); // this is particle 1, particle 2 is just acmos
		const vector_real_function_3d& p1 = q1O2gij;
		const vector_real_function_3d& p2 = mos;
		MADNESS_ASSERT(p1.size() == p2.size());
		// convolute
		part3 =
				madness::inner(world, bra1 * p1, apply(world, *fop, bra2 * p2)).sum();
	} else
		part3 = part2;

	//std::cout << "\npart1=" << part1 << "\npart2=" << part2 << "\npart3=" << part3 << "\n";
	return part1 - part2 - part3;
}

PairEnergies F12Potentials::compute_f12_energies() const {
	if (world.rank() == 0)
		std::cout << "computing " << (param.energytype()) << " f12 energies\n";

	switch (param.energytype()) {
	case PROJECTED_ENERGYTYPE:
		return compute_projected_f12_energies();
	case HYLLERAAS_ENERGYTYPE: {
		if (world.rank() == 0)
			std::cout << "!!!WARNING no PNOS were given!!!\n";
		std::valarray<vector_real_function_3d> pnos(
				zero_functions<double, 3>(world, 1), npairs());
		return compute_hylleraas_f12_energies(pnos);
	}
	default: {
		MADNESS_EXCEPTION("should not end up here", 1);
		return PairEnergies(0);
	}
	}
}

PairEnergies F12Potentials::compute_f12_energies(
		const std::valarray<vector_real_function_3d>& pnos) const {
	if (world.rank() == 0)
		std::cout << "computing " << (param.energytype()) << " f12 energies\n";

	switch (param.energytype()) {
	case PROJECTED_ENERGYTYPE:
		return compute_projected_f12_energies();
	case HYLLERAAS_ENERGYTYPE:
		return compute_hylleraas_f12_energies(pnos);
	default: {
		MADNESS_EXCEPTION("should not end up here", 1);
		return PairEnergies(0); // silence compiler
	}
	}
}

PairEnergies F12Potentials::compute_projected_f12_energies() const {
	MyTimer time = MyTimer(world).start();
	std::valarray<double> fQg = compute_fQg_integrals();
	std::valarray<double> es(npairs());
	std::valarray<double> et(npairs());
	double energy = 0.0;
	double check = 0.0;
	for (ElectronPairIterator it = pit(); it; ++it) {
		double sfactor = 1.0;
		if (it.diagonal())
			sfactor = 0.5;

		const double tfactor = 3.0;
		es[it.ij()] = sfactor * (fQg[it.ij()] + fQg[it.ij() + npairs()]);
		et[it.ij()] = tfactor * (fQg[it.ij()] - fQg[it.ij() + npairs()]);
		check += sfactor * 2.0
				* (2.0 * (fQg[it.ij()]) - fQg[it.ij() + npairs()]);
	}
	energy = es.sum() + et.sum();
	MADNESS_ASSERT(fabs(check - energy) < 1.e-5);
	PairEnergies result(npairs());
	result.eijs_f12 = es;
	result.eijt_f12 = et;
	result.energy_f12 = energy;
	result.update();
	time.stop().print("compute_projected_f12_energies");
	return result;
}

std::vector<real_function_3d> F12Potentials::read_cabs_from_file(const std::string& filename)const{
	MyTimer time_1 = MyTimer(world).start();
	if(world.rank()==0) std::cout <<"Reading CABS from external file:" << filename << "\n";
	// Read CABS Exponents from file
	std::map<std::string, std::vector<std::vector<double> > > exponents =
			basis.read_basis_from_file(param.auxbas_file(),
					nemo.get_calc()->molecule.get_atoms());
	if(world.rank()==0) std::cout << "Exponents from file:" << param.auxbas_file() << "\n";
	for (const auto& x : exponents) {
		if (world.rank() == 0) {
			std::cout << x.first << " Exponents\n";
			for (size_t l = 0; l < x.second.size(); ++l) {
				std::cout << "l=" << l << "\n";
                                using madness::operators::operator<<;
				std::cout << x.second[l] << "\n";
			}
		}
	}
	time_1.stop().print("Read Exponents From File");
	// make CABS virtuals
	vector_real_function_3d cabsf;
	MyTimer time_2 = MyTimer(world).start();
	for (const madness::Atom& atom : nemo.get_calc()->molecule.get_atoms()) {
		std::vector<std::vector<double> > exp_atom = exponents.at(
				atomic_number_to_symbol(atom.atomic_number));
		for (size_t l = 0; l < exp_atom.size(); ++l) {
			for (const double e : exp_atom[l]) {
				cabsf = append(cabsf,
						basis.guess_virtual_gaussian_shell(atom, l, e));
			}
		}
	}
	time_2.stop().print("Creating CABS Basis");
	MyTimer time_3 = MyTimer(world).start();
	cabsf = Q(cabsf);
	time_3.stop().print("Apply Q");
	print_size(world, cabsf, "cabs functions");
	return cabsf;
}

PairEnergies F12Potentials::compute_hylleraas_f12_energies(
		const std::valarray<vector_real_function_3d>& pnos) const {
	MyTimer time = MyTimer(world).start();

	std::valarray<vector_real_function_3d> abs_ij(pnos.size());
	vector_real_function_3d cabs;
	if (param.auxbas_file() != "none"){
		cabs = read_cabs_from_file(param.auxbas_file());
	}else if (param.get<std::string>("auxbas") != "none"){
		cabs = basis.guess_virtuals_internal(param.auxbas());
	}
	if (not cabs.empty()){
		MyTimer time2 = MyTimer(world).start();
		// project out pnos from cabs
		cabs = Q(cabs);
		cabs = orthonormalize_rrcd(cabs,0.1 * param.thresh());
		for (ElectronPairIterator it = pit(); it; ++it) {
			// right now this will make the same guess for all pairs
			//const vector_real_function_3d tmp=guess_virtuals(param.abs);
			QProjector<double, 3> Qpno( pnos[it.ij()]);
			const vector_real_function_3d tmp = Qpno(cabs);
			abs_ij[it.ij()] = tmp;
		}
		time2.stop().print("Make pair specific ABS from PNOS and " + std::to_string(cabs.size()) + " functions");
	}
	PairEnergies result = compute_f12_pair_energy(pnos, abs_ij);
	time.stop().print("compute hylleraas f12 energy (total) ");
	return result;
}


} /* namespace madness */
