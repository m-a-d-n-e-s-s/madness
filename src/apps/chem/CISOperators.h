/*
 * CISOperators.h
 *
 *  Created on: Aug 28, 2015
 *      Author: kottmanj
 */

#ifndef CISOPERATORS_H_
#define CISOPERATORS_H_
#include <chem/nemo.h>
#include <chem/SCFOperators.h>
namespace madness {

typedef std::vector<Function<double, 3> > vecfuncT;

class CIS_Operators {
public:
	//	CC_3D_Operator(World&world, const Nemo &nemo): world(world), mo_ket_(nemo.get_calc() -> amo),R2(init_R2(nemo.nuclear_correlation -> square())){
	//		poisson = std::shared_ptr<real_convolution_3d>(CoulombOperatorPtr(world, nemo.get_calc() -> param.lo, nemo.get_calc() ->param.econv));
	//		mo_bra_ = mul(world,R2,mo_ket_);
	//		exchange_intermediate_ = make_exchange_intermediate();
	//		if(mo_bra_.empty()) std::cout << "\n\n!!!!!WARNING: mo_bra_ vector is empty!!!!!\n\n";
	//		if(mo_ket_.empty()) std::cout << "\n\n!!!!!WARNING: mo_ket_ vector is empty!!!!!\n\n";
	//		//madness::Nuclear U(world,nemo);
	//		//nuclear_potential_ =U;
	//	}
	// If other mos than the one in the nemo struct are needed (e.g. if lower thresh is demanded -> guess calculations)
	CIS_Operators(World&world, const Nemo &nemo, const vecfuncT &mos);
	/// Make shure that R2 gets the right thresh and is constant
	real_function_3d init_R2(const Nemo &nemo) const;

	// Make the derivative of R2
	vecfuncT get_gradient(const real_function_3d f) const;

	void sanitycheck() const;

	void memory_information(const vecfuncT &v, const std::string &msg =
			"vectorfunction size is: ") const;

	std::vector<vecfuncT> make_exchange_intermediate() const;

public:
	// The nuclear potential is missing (or the U potential for the regularized approach)

	// Closed Shell Triplet CIS potential without the nuclear potential
	// returns (2J - K)x + S3C_X
	// which is in components VCIS_j =  2*\sum_i <i|r12|i> |x_j> - \sum_i <i|r12|x_j> |i> - Q\sum_i <i|r12|j> |x_i>
	vecfuncT get_CIS_potential_triplet(const vecfuncT &x) const {
		return add(world, fock_residue_closed_shell(x), S3C_X(x));
	}
	// Closed Shell Singlet CIS potential without the nuclear potential
	// returns (2J - K)x + 2*S3C_C + S3C_X
	// which is in components VCIS_j =  2*\sum_i <i|r12|i> |x_j> - \sum_i <i|r12|x_j> |i> + 2*Q\sum_i <i|r12|x_i> |j> - Q\sum_i <i|r12|j> |x_i>
	vecfuncT get_CIS_potential_singlet(const vecfuncT &x) const {
		vecfuncT S3CC = S3C_C(x);
		scale(world, S3CC, 2.0);
		vecfuncT S3CX = S3C_X(x);
		return add(world, fock_residue_closed_shell(x), add(world, S3CX, S3CC));
	}

	// Closed Shell Singlet TDA potential without the nuclear potential
	vecfuncT get_TDA_potential_singlet(const vecfuncT &x)const{
		vecfuncT S3CC = S3C_C(x); // Hartree Potential
		scale(world,S3CC,2.0);
		START_TIMER();
		real_function_3d pert_rho = make_density(mo_ket_,x);
		END_TIMER("computing untruncated perturbed density");
		START_TIMER();
		real_function_3d fxc =  xcoperator.apply_xc_kernel(pert_rho);// the exchange correlation kernel
		vecfuncT applied_fxc = mul(world,fxc,mo_ket_);
		plot_plane(world,pert_rho,"pert_rho");
		plot_plane(world,fxc,"fxc");
		plot_plane(world,applied_fxc.back(),"fxc_x_homo");
		Q(applied_fxc);
		END_TIMER("applied fxc");
		vecfuncT KSx= KS_residue_closed_shell(x);
		vecfuncT tda_pot = add(world,S3CC,applied_fxc);
		vecfuncT result = add(world,KSx,tda_pot);
		plot_plane(world,KSx.back(),"KS_residue_xhomo");
		plot_plane(world,result.back(),"smooth_potential_inops_homo");
		return result;
	}

	// Test TDA potential
	void test_tda(const bool dft,const Nemo &nemo)const{
		// F|phi_i> = epsilon_i|phi_i> ?
		vecfuncT mo = copy(world,mo_ket_);
		const double old_thresh = FunctionDefaults<3>::get_thresh();
		FunctionDefaults<3>::set_thresh(old_thresh*0.1);
		set_thresh(world,mo,FunctionDefaults<3>::get_thresh());
		vecfuncT V = get_SCF_potential(mo,dft);
		// make <i|T|i>
		std::vector < std::shared_ptr<real_derivative_3d> > gradop =gradient_operator<double, 3>(world);
		std::vector<double> T(mo.size());
		for(size_t i=0;i<mo.size();i++){
			double Ti=0.0;
			for(size_t x=0;x<3;x++){
				real_function_3d grad_ket  = (*gradop[x])(mo[i]);
				real_function_3d grad_bra  = (*gradop[x])(mo[i]*nemo.nuclear_correlation -> square());
				Ti += 0.5*grad_ket.inner(grad_bra);
			}
			T[i]=Ti;
		}
		// nuclear potential
		vecfuncT Ux;
		{
			Nuclear U(world,&nemo);
			Ux = U(mo);
			truncate(world,Ux);
		}
		V = add(world,V,Ux);

		std::vector<vecfuncT> tmp(1,mo);
		std::cout << "Kinetik1\n" << get_matrix_kinetic(tmp);
		std::cout << "\nKinetik2\n" << T << std::endl;
		std::vector<double> eps(mo_ket_.size());
		for(size_t i=0;i<eps.size();i++){
			eps[i] = T[i] + (mo[i]*nemo.nuclear_correlation->square()).inner(V[i]);
		}

		std::cout << "Test SCF results:\n";
		std::cout << eps << std::endl;

		std::cout << "\n\n comparison with moldft\n\n";
		for(size_t i=0;i<eps.size();i++){
			std::cout << eps[i] << " | " << orbital_energies_(i) << " | difference = " << eps[i]-orbital_energies_(i)  <<  std::endl;
		}
		FunctionDefaults<3>::set_thresh(old_thresh);
	}


	// Closed Shell potential for Virtuals (or SCF MOs)
	vecfuncT get_SCF_potential(const vecfuncT &x, const bool dft) const {
		if(dft){
			std::cout << " Making KS residue\n ";
			return KS_residue_closed_shell(x);
		}
		else{
			std::cout << " Making HF residue\n ";
			return fock_residue_closed_shell(x);
		}
	}

	// get the ground state density
	real_function_3d make_density() const {
		return make_density(mo_ket_, mo_ket_);
	}

	// Make a density out of two vectorfunctions f and g
	// density = \sum_i |f_i><g_i|
	real_function_3d make_density(const vecfuncT &bra_, const vecfuncT &ket_)const{
		MADNESS_ASSERT(bra_.size()==ket_.size());
		vecfuncT bra = copy(world,bra_);
		vecfuncT ket = copy(world,ket_);
		real_function_3d density = real_factory_3d(world);
		for(size_t i=0;i<bra.size();i++){
			bra[i].refine();
			ket[i].refine();
			density += bra[i]*ket[i];
		}
		if(use_nuclear_correlation_factor_){
			std::cout << "Making Density*R2" << std::endl;
			density = density*R2;
		}
		return density;
	}

	// The Fock operator is partitioned into F = T + Vn + R
	// the fock residue R= 2J-K for closed shell is computed here
	// J_j = \sum_i <i|r12|i> |tau>
	// K_j = \sum_i <i|r12|tau_j> |i>
	vecfuncT fock_residue_closed_shell(const vecfuncT &tau) const;

	// The same residue for the case that the Fock operator is the Kohn-Sham Operator
	vecfuncT KS_residue_closed_shell(
			const vecfuncT &tau) const;

	// Kinetik energy
	// -1/2 <x|R2Nabla2|x> = +1/2 <Nabla R2 x | Nabla x> = grad(R2x)*grad(x)
	// grad(R2x) = GradR2*x + R2*gradx
	// grad(R2x)*grad(y) = GradR2*x*Grady + R2*Gradx*Grady
	double get_matrix_element_kinetic_energy(const vecfuncT &ket,
			const vecfuncT &bra) const;
	double get_matrix_element_kinetic_2(const vecfuncT &bra,
			const vecfuncT &ket) const;

	// Kinetic part of the CIS perturbed fock matrix
	Tensor<double> get_matrix_kinetic(const std::vector<vecfuncT> &x) const;

	// Diagrammatic Potentials:

	// The coulomb Term of the S3C diagram: Positive sign
	// \     /
	//  \---/  = Q\sum_j(<j|g12|tau_j>)|i>
	//  _\_/_
	vecfuncT S3C_C(const vecfuncT &tau) const;

	// The Exchange Term of the S3C diagram: Negative sign
	// \  /
	//  \/...   = Q\sum_j(<j|g12|i>|tau_j>)
	//     / \
	//    _\_/_
	vecfuncT S3C_X(const vecfuncT &tau) const;

	// Project out the occupied space
	void Q(vecfuncT &f) const ;
	void Q(real_function_3d &f) const;

	// Make an inner product between vecfunctions
	double make_inner_product(const vecfuncT &bra, const vecfuncT &ket) const;
	// inner product between functions
	double make_inner_product(const real_function_3d &bra,
			const real_function_3d &ket) const;
	// inner product between function and vecfunction
	double make_inner_product(const real_function_3d &bra,
			const vecfuncT &ket) const;

private:
	bool use_timer_ = true;
	World &world;
	XCOperator xcoperator;
	bool use_nuclear_correlation_factor_;
	vecfuncT mo_bra_, mo_ket_;
	Tensor<double> orbital_energies_;
	/// The squared nuclear correlation factor and its derivative;
	const real_function_3d R2;
	vecfuncT dR2;
	std::vector<vecfuncT> exchange_intermediate_;
	std::shared_ptr<real_convolution_3d> poisson;
	//	Nuclear nuclear_potential_;
	void error(const std::string &msg) const {
		std::cout
				<< "\n\n\n !!!! ERROR IN CC_3D_OPERATOR CLASS:\n ERROR MESSAGE IS: "
				<< msg << "\n";
		MADNESS_EXCEPTION("!!!!ERROR IN CC_3D_OPERATOR CLASS!!!!", 1);
	}
	// Timer
	mutable double ttt, sss;
	void START_TIMER() const {
		if (use_timer_)
			world.gop.fence();
		ttt = wall_time();
		sss = cpu_time();
	}

	void END_TIMER(const std::string msg) const {
		if (use_timer_)
			END_TIMER(msg.c_str());
	}

	void END_TIMER(const char* msg) const {
		if (use_timer_) {
			ttt = wall_time() - ttt;
			sss = cpu_time() - sss;
			if (world.rank() == 0)
				printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
		}
	}

};


} /* namespace madness */

#endif /* CISOPERATORS_H_ */
