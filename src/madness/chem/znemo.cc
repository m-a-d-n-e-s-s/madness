/*
 * Nemocomplex.cc
 *
 *  Created on: 14 Nov 2018
 *      Author: fbischoff
 */

#include <madness/mra/mra.h>
#include "znemo.h"
#include <madness/world/test_utilities.h>
#include<madness/chem/masks_and_boxes.h>
#include<madness/chem/MolecularOrbitals.h>
#include <madness/mra/magnetic_kernels.h>
#include <madness/chem/RestartPlan.h>



namespace madness {

namespace {
/// true if the molecule is invariant under inversion through the origin (same charge at r and -r)
bool centrosymmetric(const Molecule& mol) {
	for (const Atom& a : mol.get_atoms()) {
		bool found=false;
		for (const Atom& b : mol.get_atoms()) {
			if (a.q==b.q and (a.get_coords()+b.get_coords()).normf()<1.e-6) {found=true; break;}
		}
		if (not found) return false;
	}
	return true;
}

/// <f|P f>/<f|f> with P f(r) = f(-r): +1 gerade, -1 ungerade; the field Hamiltonian commutes with P in the
/// symmetric and the Landau gauge (A(-r) = -A(r)), and the Landau phase e^{i(B/2)xy} is even
Tensor<double> inversion_parity(const std::vector<complex_function_3d>& f) {
	Tensor<double> p(long(f.size()));
	for (std::size_t i=0; i<f.size(); ++i) {
		const complex_function_3d pf=mirror(f[i],std::vector<long>{-1,-1,-1});
		p(long(i))=real(inner(f[i],pf))/real(inner(f[i],f[i]));
	}
	return p;
}
}

Znemo::Znemo(World& world, const commandlineparser& parser) : NemoBase(world), mol(world,parser), param(world,parser),
                    nemo_param(world,parser), calc_param(world,parser) {
	calc_param.set_derived_value("k",9);
	calc_param.read_input_and_commandline_options(world,parser,"dft");

    FunctionDefaults<3>::set_k(calc_param.k());
    // truncation threshold: the last value of the dft block's protocol (znemo runs a single protocol step),
    // independent of the energy convergence econv
    FunctionDefaults<3>::set_thresh(calc_param.protocol().back());
    FunctionDefaults<3>::set_refine(true);
    // initial level of the original magnetic-kernel codes; truncation tightened by 2^-n on fine levels
    // (mode 1): the variational energy measures the orbital noise through |pi phi|^2, i.e. derivatives,
    // so absolute truncation (mode 0) leaves it scattering at strong field (mode 2 makes the kernel
    // application far more expensive)
    FunctionDefaults<3>::set_initial_level(2);
    FunctionDefaults<3>::set_truncate_mode(1);

    aobasis.read_file(calc_param.aobasis());
//    cparam.set_molecular_info(mol, aobasis, 0);
    calc_param.set_derived_values(mol);
    calc_param.set_derived_value("econv",0.1*FunctionDefaults<3>::get_thresh());

    FunctionDefaults<3>::set_cubic_cell(-calc_param.L(), calc_param.L());
    calc_param.set_derived_value("spin_restricted",false);

	print_info=printleveler(calc_param.print_level());

	if (world.rank()==0 and print_info.print_setup()) {
		param.print("complex","end");
		calc_param.print("dft","end");
		mol.print();
	}

	update_smoothing_radii();
	MagneticKernelParameters kp(world,parser);
	set_kernel_derived_values(kp,FunctionDefaults<3>::get_thresh());
	if (world.rank()==0) kp.print("magnetic_kernel","end");
	mkernel=make_magnetic_kernel(world,param.kernel(),param.physical_B(),kp);

	B={0,0,param.physical_B()};


    // set the linear moments
    for (int i=0; i<3; ++i) rvec.push_back(real_factory_3d(world).functor([&i] (const coord_3d& r) {return r[i];}));

	coulop=std::shared_ptr<real_convolution_3d>(CoulombOperatorPtr(world,calc_param.lo(),0.1*FunctionDefaults<3>::get_thresh()));
};

bool Znemo::need_recompute_factors_and_potentials(const double thresh) const {
	bool need=NemoBase::need_recompute_factors_and_potentials(thresh);
	if ((not (potentialmanager.get() and potentialmanager->vnuclear().is_initialized()))
	        		or (potentialmanager->vnuclear().thresh()>thresh)) need=true;
    return need;
}

void Znemo::invalidate_factors_and_potentials() {
    NemoBase::invalidate_factors_and_potentials();
	potentialmanager.reset();
	// mkernel holds the input parameters: rebuilt in recompute_factors_and_potentials when they change
}

void Znemo::recompute_factors_and_potentials(const double thresh) {

	print("recomputing factors and potentials in znemo");

	// the guess is read from a previous nemo calculation
	// make sure the molecule was not reoriented there
	if (not mol.parameters.no_orient()) {
		MADNESS_EXCEPTION("the molecule of the reference calculation was reoriented\n\n",1);
	}

	update_smoothing_radii();
	potentialmanager=std::shared_ptr<PotentialManager>(new PotentialManager(mol, molecule().parameters.core_type()));
	potentialmanager->make_nuclear_potential(world);
	// the magnetic kernels act on the orbitals, or with a nuclear correlation factor (nemo block, ncf) on the nemos F = psi/R: then
	// R^{-1} (pi^2/2) R F = pi^2/2 F + U1.grad F + i U1.A F - (grad^2 R/2R) F, U1 = -grad R/R, so the kernel
	// G_B(E) is unchanged and its right-hand side is U2 F + U1.grad F + i U1.A F + (J - K) F (Nuclear(ncf),
	// zeeman_R_comm), with the metric R^2 in the Fock matrix and the energy
    construct_nuclear_correlation_factor(mol, potentialmanager, nemo_param.ncf());

	MagneticKernelParameters kp=mkernel->parameters();
	set_kernel_derived_values(kp,thresh);
	if (kp!=mkernel->parameters()) {
		if (world.rank()==0) kp.print("magnetic_kernel","end");
		mkernel=make_magnetic_kernel(world,param.kernel(),param.physical_B(),kp);
	}
	if (world.rank()==0) print("magnetic kernel",mkernel->name(),"for B =",param.physical_B());

    // (re) construct the Poisson solver
	coulop=std::shared_ptr<real_convolution_3d>(CoulombOperatorPtr(world,get_calc_param().lo(),0.1*thresh));

}

void Znemo::update_smoothing_radii() {
	mol.update_rcut_with_eprec(mol.get_eprec(),param.physical_B());
	if (world.rank()==0) {
		std::vector<double> c;
		for (double rc : mol.get_rcut()) c.push_back(1.0/rc);
		print("nuclear smoothing radii c",c,"eprec",mol.get_eprec(),"B",param.physical_B());
	}
}

void Znemo::set_kernel_derived_values(MagneticKernelParameters& kp, const double thresh) const {
	std::vector<MagneticKernelNucleus> nuclei;
	const std::vector<double> rcut=mol.get_rcut();
	for (std::size_t i=0; i<mol.natom(); ++i) {
		const Atom& a=mol.get_atom(i);
		if (a.q<=0.0) continue;
		nuclei.push_back({a.q,a.x,1.0/rcut[i],std::pow(a.q,3)*nuclear_density_field_factor(a.q,param.physical_B())/constants::pi});
	}
	kp.set_derived_values(param.physical_B(),thresh,nuclei,nemo_param.ncf().first!="none");
	kp.verbose=(get_calc_param().print_level()>=3);
}



/// compute the molecular energy
double Znemo::value(const Tensor<double>& x) {

	// Znemo's energy expression (compute_energy) has no dispersion term, and its
	// gradient runs through NemoBase::compute_gradient, which deliberately has
	// none either -- D3 is not parameterized for a magnetic-field Hamiltonian.
	DispersionCorrection::reject(get_calc_param().dispersion(), "znemo");

	if ((x-mol.get_all_coords()).normf()>1.e-12) invalidate_factors_and_potentials();

	// compute the molecular potential
	mol.set_all_coords(x.reshape(mol.natom(),3));
	if (need_recompute_factors_and_potentials(FunctionDefaults<3>::get_thresh()))
		recompute_factors_and_potentials(FunctionDefaults<3>::get_thresh());

	if (get_initial_orbitals()) {
		iterate();
		save_orbitals("final");
	}
	save_orbitals();

	double energy=analyze();
    output_calc_info_schema(energy);

	return energy;
}

void Znemo::output_calc_info_schema(const double& energy) const {
    nlohmann::json j;
    j["scf_eigenvalues_a"]=tensor_to_json(aeps);
    j["scf_eigenvalues_b"]=tensor_to_json(beps);
    j["model"]="UHF";
    j["B"]=this->B;
    j["driver"]="energy";
    j["return_energy"]=energy;
    j["kernel"]=mkernel->name();
    j["kernel_energy"]=kernel_energy;
    update_schema(get_calc_param().prefix()+".calc_info", j);
}

Tensor<double> Znemo::gradient(const Tensor<double>& x) {

	const real_function_3d rhonemo=compute_nemo_density(amo,bmo);
	const Tensor<double> grad=NemoBase::compute_gradient(rhonemo,mol);

    if (world.rank() == 0) {
        print("\n Derivatives (a.u.)\n -----------\n");
        print("  atom        x            y            z          dE/dx        dE/dy        dE/dz");
        print(" ------ ------------ ------------ ------------ ------------ ------------ ------------");
        for (size_t i = 0; i < mol.natom(); ++i) {
            const Atom& atom = mol.get_atom(i);
            printf(" %5d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n", int(i),
                   atom.x, atom.y, atom.z, grad[i * 3 + 0], grad[i * 3 + 1],
                   grad[i * 3 + 2]);
        }
    }
    return grad;
}


void Znemo::iterate() {

	double thresh=FunctionDefaults<3>::get_thresh();
	double energy=1.e10;
	double oldenergy=0.0;

	// the diamagnetic box

//	XNonlinearSolver<std::vector<complex_function_3d> ,double_complex, allocator<double_complex,3> > solvera(allocator<double_complex,3> (world,amo.size()));
//	XNonlinearSolver<std::vector<complex_function_3d> ,double_complex, allocator<double_complex,3> > solverb(allocator<double_complex,3> (world,bmo.size()));
    auto solvera= nonlinear_vector_solver<double_complex,3>(world,amo.size());
    auto solverb= nonlinear_vector_solver<double_complex,3>(world,bmo.size());
	solvera.set_maxsub(get_calc_param().maxsub()); // @suppress("Method cannot be resolved")
	solvera.do_print=(get_calc_param().print_level()>2);
	solverb.set_maxsub(get_calc_param().maxsub());
	solverb.do_print=(get_calc_param().print_level()>2);

	// set end of iteration cycles for intermediate calculations
	int maxiter = get_calc_param().maxiter();
	double dconv=get_calc_param().dconv();
	double na=1.0,nb=1.0;	// residual norms

	solvera.clear_subspace();
	solverb.clear_subspace();
	bool converged=false;

	// magnetic kernels: kernel energies E_i (without the spin Zeeman term), started from the Fock matrix
	// and then updated by E - Re<V psi|psi - psi~>/|psi~|^2 as in the original drivers: the Fock energies of
	// poor first iterates may lie above the Landau threshold, where G(E) has an unbounded range along B
	Tensor<double> akeps, bkeps;
	double kernel_energy_old=1.e10;

	// iterate the SCF cycles
	for (int iter=0; iter<maxiter; ++iter) {
//		if (iter%10==0) save_orbitals("iteration"+stringify(iter));
		save_orbitals();

		// compute the density
		real_function_3d density=compute_density(amo,bmo);

		// compute the dual of the orbitals
		std::vector<complex_function_3d> R2amo=make_bra(amo);
		std::vector<complex_function_3d> R2bmo=make_bra(bmo);

		// compute the fock matrix
		std::vector<complex_function_3d> Vnemoa, Vnemob;
		Tensor<double_complex> focka, fockb(0l,0l);
		potentials apot(world,amo.size()), bpot(world,bmo.size());

		apot=compute_potentials(amo, density, amo);		// 4.53
		Vnemoa=assemble_rhs(apot,amo);
		truncate(world,Vnemoa,thresh*0.1);
		focka=kinetic_pi_matrix(amo) + compute_vmat(amo,apot);

		if (have_beta()) {
			bpot=compute_potentials(bmo, density, bmo);
			Vnemob=assemble_rhs(bpot,bmo);
			truncate(world,Vnemob,thresh*0.1);
			fockb=kinetic_pi_matrix(bmo) + compute_vmat(bmo,bpot);
		}

		if (world.rank()==0 and (print_info.print_convergence())) {
			print("Fock matrix");
			print(focka);
			print(fockb);
		}

		oldenergy=energy;
		energy=compute_energy(amo,apot,bmo,bpot,print_info.print_convergence());
		current_energy=energy;

		Tensor<double_complex> ovlp=matrix_inner(world,R2amo,amo);
		canonicalize(amo,Vnemoa,apot,solvera,focka,ovlp);
		truncate(world,Vnemoa);
		if (have_beta()) {
			Tensor<double_complex> ovlp=matrix_inner(world,R2bmo,bmo);
			canonicalize(bmo,Vnemob,bpot,solverb,fockb,ovlp);
		}

		if (world.rank()==0 and (print_info.print_convergence())) print("using fock matrix for the orbital energies");
		for (int i=0; i<focka.dim(0); ++i) aeps(i)=real(focka(i,i));
		for (int i=0; i<fockb.dim(0); ++i) beps(i)=real(fockb(i,i));


		if (world.rank()==0 and (print_info.print_convergence())) {
			print("orbital energies alpha",aeps);
			print("orbital energies beta ",beps);
		}
		if ((world.rank() == 0) and print_info.print_energies()) {
			printf("finished iteration %2d at time %8.1fs with energy, norms %12.8f %12.8f %12.8f\n",
					iter, wall_time(), energy, na, nb);
		}

		// test the kernel-consistent energy, as the original drivers did -- the variational energy (derivatives of
		// orbitals with strong phases) scatters at strong field
		const double denergy=kernel_energy-kernel_energy_old;
		if (std::abs(denergy)<get_calc_param().econv() and (sqrt(na*na+nb*nb)<dconv)) {
			print("energy converged");
			converged_for_thresh=FunctionDefaults<3>::get_thresh();
			converged_for_dconv=dconv;
			save_orbitals("converged");
			converged=true;
		}
		if (converged) break;

		// compute the residual of the Greens' function
		// the spin Zeeman term +-B/2 is a constant per spin: it is not part of the orbital energies (and drops out
		// of the iteration), only of the total energy (compute_energy)
		if (ncf->type()!=NuclearCorrelationFactor::None) {
			std::vector<complex_function_3d> g=Vnemoa;
			for (auto& f : Vnemob) g.push_back(f);
			Tensor<double> E(long(g.size()));
			for (long i=0; i<aeps.size(); ++i) E(i)=aeps(i);
			for (long i=0; i<beps.size(); ++i) E(aeps.size()+i)=beps(i);
			adapt_kernel_tmin(g,E,iter==0);
		}
		if (iter==0) {
			akeps=copy(aeps);	// Tensor assignment is shallow: the kernel energies must not alias the Fock diagonal
			bkeps=copy(beps);
		}
		std::vector<complex_function_3d> resa=compute_residuals(Vnemoa,amo,akeps);

		Tensor<double> normsa=real(inner(world,make_bra(resa),resa));
		na=sqrt(normsa.sumsq());

		// truncate the update at 0.1*thresh and not again after orthonormalization -- each
		// truncation at thresh perturbs the converged orbitals, which the variational energy picks up through
		// |pi phi|^2, and the kernel energies E_i stall at a noise floor above econv
		std::vector<complex_function_3d> amo_new=solvera.update(amo,resa,0.01,3);
		amo_new=truncate(amo_new,0.1*thresh);

		do_step_restriction(amo,amo_new);

		amo=orthonormalize(amo_new);

		if (have_beta()) {
			std::vector<complex_function_3d> resb=compute_residuals(Vnemob,bmo,bkeps);
			Tensor<double> normsb=real(inner(world,make_bra(resb),resb));
			nb=sqrt(normsb.sumsq());

			std::vector<complex_function_3d> bmo_new=solverb.update(bmo,resb,0.01,3);
			bmo_new=truncate(bmo_new,0.1*thresh);
			do_step_restriction(bmo,bmo_new);
			bmo=orthonormalize(bmo_new);
		} else {
			nb=0.0;
		}
		kernel_energy_old=kernel_energy;
		// E_var + sum_i (E_i - F_ii) = sum_i E_i - (J - K) + E_nuc + E_spin-Zeeman
		kernel_energy=energy;
		for (long i=0; i<akeps.size(); ++i) kernel_energy+=akeps(i)-aeps(i);
		for (long i=0; i<bkeps.size(); ++i) kernel_energy+=bkeps(i)-beps(i);
		if ((world.rank()==0) and print_info.print_energies())
			printf("kernel-consistent energy %12.8f\n",kernel_energy);
	}
}


std::vector<complex_function_3d> Znemo::orthonormalize(const std::vector<complex_function_3d>& mo) const {
	std::vector<complex_function_3d> result=copy(world,mo);
	NemoBase::orthonormalize(result,R);
	return result;
}

std::vector<complex_function_3d> Znemo::normalize(const std::vector<complex_function_3d>& mo) const {
	std::vector<complex_function_3d> result=copy(world,mo);
	NemoBase::normalize(result,R);
	return result;
}


real_function_3d Znemo::compute_density(const std::vector<complex_function_3d>& amo,
		const std::vector<complex_function_3d>& bmo) const {
	return (compute_nemo_density(amo,bmo)*R_square).truncate();
}

real_function_3d Znemo::compute_spin_density(const std::vector<complex_function_3d>& amo,
		const std::vector<complex_function_3d>& bmo) const {
	return (compute_nemo_spin_density(amo,bmo)*R_square).truncate();
}

std::vector<complex_function_3d> Znemo::make_bra(const std::vector<complex_function_3d>& rhs) const {
	return truncate(R_square*rhs);
}

std::vector<complex_function_3d> Znemo::assemble_rhs(const potentials& pot,
		const std::vector<complex_function_3d>& mo) const {
	// the kinetic energy and the diamagnetic term are in the Green's function, the
	// spin Zeeman term in its energy argument (compute_residuals); the paramagnetic term only if the
	// kernel does not contain it (mehler, symmetric gauge)
	std::vector<complex_function_3d> v=pot.vnuc_mo-pot.K_mo+pot.J_mo+pot.zeeman_R_comm;
	if (not mkernel->absorbs_paramagnetic()) {
		MADNESS_CHECK_THROW(mkernel->gauge()==MagneticGauge::symmetric,"paramagnetic term needs the symmetric gauge");
		// B-spline derivative: with the Alpert derivative the iteration diverged at B = 100
		v+=0.5*B[2]*Lz<double_complex,3>(world,true)(mo);
	}
	return v;
}

std::vector<complex_function_3d> Znemo::pi_apply(const std::vector<complex_function_3d>& f, const int axis) const {
	const MagneticKernel& K=*mkernel;
	complex_derivative_3d D=free_space_derivative<double_complex,3>(world,axis);
	std::vector<complex_function_3d> pf=apply(world,D,f);
	scale(world,pf,double_complex(0.0,-1.0));
	// symmetric gauge: A_x, A_y; Landau gauge: A_x only
	const bool have_A=(param.physical_B()!=0.0) and (axis==0 or (axis==1 and K.gauge()==MagneticGauge::symmetric));
	if (have_A) {
		real_function_3d A=real_factory_3d(world).functor([&K,axis](const coord_3d& r) {return K.vector_potential(r)[axis];});
		pf=pf+A*f;
	}
	return pf;
}

/// <psi_i|pi^2/2|psi_j> = (1/2) <pi psi_i|pi psi_j>, Hermitian; with a correlation factor (psi = R F) without
/// differentiating the cusp of R: pi (R F) = R (pi + i U1) F, U1 = -grad R/R, so the matrix is
/// (1/2) sum_d <R^2 w_id|w_jd>, w_d = (pi_d + i U1_d) F (contains L_z and A^2 through pi)
Tensor<double_complex> Znemo::kinetic_pi_matrix(const std::vector<complex_function_3d>& mo) const {
	Tensor<double_complex> t(mo.size(),mo.size());
	if (mo.size()==0) return t;
	const bool nemo=kernel_nemo();
	for (int axis=0; axis<3; ++axis) {
		std::vector<complex_function_3d> w=pi_apply(mo,axis);
		if (nemo) w=w+double_complex(0.0,1.0)*(ncf->U1(axis)*mo);
		t+=0.5*matrix_inner(world,nemo ? make_bra(w) : w,w);
	}
	return t;
}

/// nemo right-hand sides g = U F are cusp-free: tmin = MagneticKernelParameters::tmin_cusp_free(max_i C_i),
/// C = <g|T - E|g> = (1/2)|pi g|^2 - E |g|^2, not below the bare value lo^2/4. C grows as the orbitals contract,
/// so it is re-evaluated in every iteration: tmin may grow in the first one only, and shrinks later when needed.
void Znemo::adapt_kernel_tmin(const std::vector<complex_function_3d>& g, const Tensor<double>& E, const bool first) {
	if (ncf->type()==NuclearCorrelationFactor::None or g.empty()) return;
	if (mkernel->parameters().is_user_defined("tmin")) return;
	double cmax=0.0;
	std::vector<double> pn(g.size(),0.0);
	for (int axis=0; axis<3; ++axis) {
		const std::vector<double> n=norm2s(world,pi_apply(g,axis));
		for (std::size_t i=0; i<g.size(); ++i) pn[i]+=n[i]*n[i];
	}
	const std::vector<double> gn=norm2s(world,g);
	for (std::size_t i=0; i<g.size(); ++i) cmax=std::max(cmax,0.5*pn[i]-std::min(E(long(i)),mkernel->bottom())*gn[i]*gn[i]);
	const MagneticKernelParameters& p=mkernel->parameters();
	const double tmin_bare=0.25*p.lo()*p.lo();
	const double tmin_now=p.tmin();
	const double tmin=std::max(tmin_bare,p.tmin_cusp_free(cmax));
	if (world.rank()==0 and print_info.print_convergence())
		print("nemo: C = <g|T-E|g> =",cmax,"-> tmin",tmin,"(current",tmin_now,", bare",tmin_bare,")");
	if ((first and tmin>1.01*tmin_now) or tmin<0.7*tmin_now) {
		if (world.rank()==0) print("nemo: kernel rebuilt with tmin",tmin);
		MagneticKernelParameters q=p;
		q.set_derived_value("tmin",tmin);
		mkernel=make_magnetic_kernel(world,param.kernel(),param.physical_B(),q);
	}
}

double Znemo::analyze() {

	if (need_recompute_factors_and_potentials(FunctionDefaults<3>::get_thresh()))
		recompute_factors_and_potentials(FunctionDefaults<3>::get_thresh());

	real_function_3d density=compute_density(amo,bmo);
	real_function_3d spindensity=compute_spin_density(amo,bmo);

	potentials apot(world,amo.size()), bpot(world,bmo.size());
	apot=compute_potentials(amo,density,amo);
	if (have_beta()) {
		bpot=compute_potentials(bmo,density,bmo);
	}

	double energy=compute_energy(amo,apot,bmo,bpot,true);

	if (world.rank()==0) {
		print("orbital energies alpha",aeps);
		print("orbital energies beta ",beps);
	}

	// the analysis below assumes the symmetric gauge: transform Landau-gauge orbitals, restore at the end
	const std::vector<complex_function_3d> amo_kernel=amo, bmo_kernel=bmo;
	if (mkernel->gauge()!=MagneticGauge::symmetric) {
		for (auto& mo : amo) mo=mkernel->to_symmetric_gauge(mo).truncate();
		for (auto& mo : bmo) mo=mkernel->to_symmetric_gauge(mo).truncate();
		density=compute_density(amo,bmo);
		spindensity=compute_spin_density(amo,bmo);
	}

	// compute the current density
	std::vector<real_function_3d> j=compute_current_density(amo,bmo);
	save(j[0],"j0");
	save(j[1],"j1");
	save(j[2],"j2");

	Lz<double_complex,3> lz(world);
	// compute the expectation values of the Lz operator
	std::vector<complex_function_3d> lzamo=0.5*B[2]*lz(amo);
	std::vector<complex_function_3d> lzbmo=0.5*B[2]*lz(bmo);
	Tensor<double_complex> lza_exp=inner(world,amo,lzamo);
	Tensor<double_complex> lzb_exp=inner(world,bmo,lzbmo);
	print("< amo | lz | amo >",lza_exp);
	print("< bmo | lz | bmo >",lzb_exp);
	if (centrosymmetric(mol)) {
		print("inversion parity <P> amo (g +1, u -1)",inversion_parity(amo));
		print("inversion parity <P> bmo (g +1, u -1)",inversion_parity(bmo));
	}

    std::vector<real_function_3d> A=compute_magnetic_vector_potential(world,B);

    Tensor<double> p_exp=compute_kinetic_momentum();
    Tensor<double> A_exp=compute_magnetic_potential_expectation(A);

    print("<p>       ",p_exp);
    print("<A>       ",A_exp);
    print("(p-eA)    ",p_exp+A_exp);

    // compute the standard kinetic gauge origin, defined as the gauge origin, where the
    // expectation value of the kinetic momentum p vanishes
    const long nmo=amo.size()+bmo.size();
    Tensor<double> v=-1.0/nmo*p_exp;
    Tensor<double> S=compute_standard_gauge_shift(p_exp);
    print("standard gauge shift S",S);

    // compute the standardized components of the canonical momentum square
    const double v2=v.trace(v);	// term B^2(Sx^2+Sy^2)
    const double vp=v.trace(p_exp);
    const double vA=v.trace(A_exp);
    print("v2, vp, vA", v2, vp, vA);

    print("expectation values in standard gauge");
    print("Delta 1/2 <p^2>  ",0.5*(nmo*v2 + 2.0*vp));

    print("saving orbitals for plotting");
	save_orbitals("plot");
	save(density,"density");
	save(spindensity,"spindensity");

	auto components=std::vector<std::vector<real_function_3d>>({real(amo),imag(amo),real(bmo),imag(bmo)});
	auto component_names=std::vector<std::string>({"real_amo","imag_amo","real_bmo","imag_bmo"});
	std::vector<real_function_3d> real_aos=SCF::project_ao_basis_only(world, aobasis, mol);
	for (size_t i=0; i<components.size(); ++i) {
		if (world.rank()==0) print("analysis of MO component ",component_names[i]);
		if (components[i].size()>0) SCF::analyze_vectors(world, components[i], real_aos,
			FunctionDefaults<3>::get_thresh()*0.1, molecule(), get_calc_param().print_level(), aobasis);

	}
	amo=amo_kernel;
	bmo=bmo_kernel;

    return energy;
}


double Znemo::compute_energy(const std::vector<complex_function_3d>& amo, const Znemo::potentials& apot,
		const std::vector<complex_function_3d>& bmo, const Znemo::potentials& bpot, const bool do_print) const {

	timer tenergy(world,print_info.print_timings());
    double fac= get_calc_param().spin_restricted() ? 2.0 : 1.0;
    std::vector<complex_function_3d> dia2amo=make_bra(amo);
    std::vector<complex_function_3d> dia2bmo=make_bra(bmo);


    // 1/2 |pi phi|^2 contains the paramagnetic and the diamagnetic terms
    double_complex kinetic=0.0;
    Tensor<double_complex> ta=kinetic_pi_matrix(amo), tb=kinetic_pi_matrix(bmo);
    for (long i=0; i<ta.dim(0); ++i) kinetic+=fac*ta(i,i);
    for (long i=0; i<tb.dim(0); ++i) kinetic+=fac*tb(i,i);

    const bool sym=kernel_nemo();	// symmetric kinetic term: bare V, no R commutator (see compute_vmat)
    double_complex nuclear_potential=sym
    		? fac*(inner(world,dia2amo,apot.vbare_mo).sum()+inner(world,dia2bmo,bpot.vbare_mo).sum())
    		: fac*(inner(world,dia2amo,apot.vnuc_mo).sum()+inner(world,dia2bmo,bpot.vnuc_mo).sum());
    // spin Zeeman term +B/2 per alpha and -B/2 per beta electron: in the total energy only
    const double_complex spin_zeeman=0.5*B[2]*fac*(double(amo.size())-double(bmo.size()));
    double_complex coulomb=fac*0.5*(inner(world,dia2amo,apot.J_mo).sum()+inner(world,dia2bmo,bpot.J_mo).sum());
    double_complex exchange=fac*0.5*(inner(world,dia2amo,apot.K_mo).sum()+inner(world,dia2bmo,bpot.K_mo).sum());
    double_complex zeeman_R_comm=sym ? 0.0 : fac*(inner(world,dia2amo,apot.zeeman_R_comm).sum()+inner(world,dia2bmo,bpot.zeeman_R_comm).sum());

    double_complex energy=kinetic + nuclear_potential + mol.nuclear_repulsion_energy() +
    		spin_zeeman + coulomb - exchange + zeeman_R_comm;

	if (world.rank()==0 and do_print) {
		printf("  kinetic energy        %12.8f \n", real(kinetic));
		printf("  nuclear potential     %12.8f \n", real(nuclear_potential));
		printf("  nuclear repulsion     %12.8f \n", mol.nuclear_repulsion_energy());
		printf("  orbital zeeman R comm %12.8f \n", real(zeeman_R_comm));
		printf("  spin zeeman term      %12.8f \n", real(spin_zeeman));
		printf("  Coulomb               %12.8f \n", real(coulomb));
		printf("  exchange              %12.8f \n", real(exchange));
		printf("  total                 %12.8f \n", real(energy));
	}

	if(fabs(imag(energy))>1.e-7) {


		print("real part");
		printf("  kinetic energy      %12.8f \n", real(kinetic));
		printf("  nuclear potential   %12.8f \n", real(nuclear_potential));
		printf("  nuclear repulsion   %12.8f \n", mol.nuclear_repulsion_energy());
		printf("  orbital zeeman R comm %12.8f \n", real(zeeman_R_comm));
		printf("  spin zeeman term    %12.8f \n", real(spin_zeeman));
		printf("  Coulomb             %12.8f \n", real(coulomb));
		printf("  exchange            %12.8f \n", real(exchange));
		printf("  total               %12.8f \n", real(energy));

		print("imaginary part");
		printf("  kinetic energy      %12.8f \n", imag(kinetic));
		printf("  nuclear potential   %12.8f \n", imag(nuclear_potential));
		printf("  nuclear repulsion   %12.8f \n", 0.0);
		printf("  orbital zeeman R comm %12.8f \n", imag(zeeman_R_comm));
		printf("  spin zeeman term    %12.8f \n", imag(spin_zeeman));
		printf("  Coulomb             %12.8f \n", imag(coulomb));
		printf("  exchange            %12.8f \n", imag(exchange));
		printf("  total               %12.8f \n", imag(energy));

		print("imaginary part of the energy",energy.imag());

		MADNESS_EXCEPTION("complex energy computation.. ",1);
	}

	tenergy.end("compute_energy");
	return real(energy);
}

/// following Lazzeretti, J. Mol. Struct, 313 (1994)

/// see also the first three terms of Eq. (47) of Hodge, Am. J. Phys., 82 (2014)
std::vector<real_function_3d> Znemo::compute_current_density(
		const std::vector<complex_function_3d>& alpha_mo,
		const std::vector<complex_function_3d>& beta_mo) const {

	// compute density and spin density
	real_function_3d adens=real_factory_3d(world);
	real_function_3d bdens=real_factory_3d(world);
	for (auto& mo : alpha_mo) adens+=abs_square(mo);
	for (auto& mo : beta_mo) bdens+=abs_square(mo);
	real_function_3d density=adens+bdens;
	real_function_3d spin_density=adens-bdens;

//	double ndiff=(spin_density*R_square).trace();
//	print("nalpha - nbeta",ndiff);
//	double nel=(density*R_square).trace();
//	print("nelectron",nel);

	std::vector<real_function_3d> vspin_density=zero_functions_compressed<double,3>(world,3);;
	vspin_density[2]=spin_density;

	// compute first contribution to current density from orbitals: Re (\Psi \hat p \Psi)
	// p = -i del
	// Re(-i z) + Re( -i a - i^2 b) = b  =Im (z)
	std::vector<real_function_3d> j=zero_functions_compressed<double,3>(world,3);
	for (auto& mo : alpha_mo) j+=imag(conj(mo)*grad(mo));
	for (auto& mo : beta_mo) j+=imag(conj(mo)*grad(mo));
//
//	real_function_3d null1=div(j*R_square);
//	double n31=null1.norm2();
//	print("div(j)",n31);

    // compute the magnetic potential and density contribution: A \rho
	std::vector<real_function_3d> A=compute_magnetic_vector_potential(world,B);
//	auto vAdens=A*density;
	j+=A*density;
//	real_function_3d null2=div(vAdens*R_square);
//	double n32=null2.norm2();
//	print("div(A*dens)",n32);

	// spin density contribution: rot(\Psi \hat S \Psi)
	j+=0.5*rot(vspin_density);
	j[0]+=0.5*2.0*(-1.0)*ncf->U1(1)*spin_density;
	j[1]-=0.5*2.0*(-1.0)*ncf->U1(0)*spin_density;

//	real_function_3d null3=div(j*R_square);
//	double n33=null3.norm2();
//	print("div(j)",n33);

//	save(j[0],"j0reg");
//	save(j[1],"j1reg");
//	save(j[2],"j2reg");

//	auto divj=div(j);
//	save(divj,"div_reg_j");
//	auto term1=2.0*dot(world,ncf->U1vec(),j);
//	save(term1,"U1term");

	// sanity check
	real_function_3d null=R_square*(-2.0*dot(world,ncf->U1vec(),j) + div(j,true));
	save(null,"null");
	double n3=null.norm2();
	print("final div(j)",n3);

	j=truncate(j*R_square);
//	auto jnorm=norm2s(world,j);
//	print("norm(j)",jnorm);

	return j;
}


void Znemo::test_compute_current_density() const {

	Lz<double_complex,3> lz(world);
	complex_function_3d pp=complex_factory_3d(world).functor(p_orbital(1,1.0,{0,0,0}));
	double norm=pp.norm2();
	pp.scale(1/norm);
	double_complex l=inner(amo[0],0.5*B[2]*lz(amo[0]));
	print("<p+|Lz|p+>",l);

//	std::vector<complex_function_3d> vpp(1,amo);
	std::vector<real_function_3d> j=compute_current_density(amo,std::vector<complex_function_3d>());

	save(j[0],"j0");
	save(j[1],"j1");
	save(j[2],"j2");

}




void Znemo::do_step_restriction(const std::vector<complex_function_3d>& mo,
		std::vector<complex_function_3d>& mo_new) const {
    PROFILE_MEMBER_FUNC(SCF);
    std::vector<double> anorm = norm2s(world, sub(world, mo, mo_new));
    int nres = 0;
    for (unsigned int i = 0; i < mo.size(); ++i) {
        if (anorm[i] > get_calc_param().maxrotn()) {
            double s = get_calc_param().maxrotn() / anorm[i];
            ++nres;
            if (world.rank() == 0 and print_info.print_convergence()) {
                if (nres == 1)
                    printf("  restricting step ");
                printf(" %d", i);
            }
            mo_new[i].gaxpy(s, mo[i], 1.0 - s, false);
        }
    }
    if (nres > 0 && (world.rank() == 0) and print_info.print_convergence())
        printf("\n");

    world.gop.fence();
}


/// the restart sequence of moldft and nemo (RestartPlan): the restartdata of a previous znemo calculation if
/// the dft restart and the archive allow it, otherwise the orbitals of the previous geometry or the hcore guess
/// @return false if the orbitals read are the answer (restart read_only, or auto with a converged archive)
bool Znemo::get_initial_orbitals() {

	// znemo reads its own restartdata only: no AO projections, no NWChem movecs (as nemo)
	RestartPlan plan=make_restart_plan(world,restart_mode_from_string(get_calc_param().restart()),get_calc_param(),
			mol,Representation::znemo,RestartCapabilities::restartdata_only(),restart_ncf(),param.physical_B());

	if (plan.needs_load()) {
		MADNESS_CHECK(plan.source==RestartSource::restartdata);
		if (world.rank()==0) print("reading the orbitals from",get_calc_param().prefix()+".restartdata");
		load_restartdata();
	} else if (amo.empty() and bmo.empty()) {
		auto zmos=hcore_guess();
		amo=zmos.first.get_mos();
		bmo=zmos.second.get_mos();
		converged_for_thresh=converged_for_dconv=1.e10;
	} else {
		// a later geometry: the orbitals of the previous one are the better guess
		if (world.rank()==0) print("starting from the orbitals of the previous geometry");
		converged_for_thresh=converged_for_dconv=1.e10;
	}

	// reading can invalidate the plan: orbitals reprojected on reading are no longer converged
	const double target_thresh=get_calc_param().protocol().back();
	const double target_dconv=std::max(target_thresh,get_calc_param().dconv());
	if (not plan.iterate and plan.mode==RestartMode::automatic
			and not (converged_for_thresh<=target_thresh and converged_for_dconv<=target_dconv)) {
		plan.iterate=true;
		if (world.rank()==0) print("the orbitals had to be reprojected on reading and are no longer converged; iterating");
	}
	if (not plan.iterate and world.rank()==0) print("not solving the SCF equations:",plan.why);

	MADNESS_CHECK_THROW(long(amo.size())==get_calc_param().nalpha() and long(bmo.size())==get_calc_param().nbeta(),
			"initial orbitals: wrong number of orbitals");
	if (aeps.size()!=long(amo.size())) aeps=Tensor<double>(amo.size());
	if (beps.size()!=long(bmo.size())) beps=Tensor<double>(bmo.size());

	if (plan.iterate) {
		amo=orthonormalize(amo);
		bmo=orthonormalize(bmo);
	}
	set_thresh(world,amo,FunctionDefaults<3>::get_thresh());
	set_thresh(world,bmo,FunctionDefaults<3>::get_thresh());
	return plan.iterate;
}


/// write <prefix>.restartdata: the RestartMetadata header (representation znemo, with B and eprec) and one
/// block per spin as in SCF::save_mos
void Znemo::save_orbitals() const {
	if (not get_calc_param().save()) return;
	const std::string name=get_calc_param().prefix()+".restartdata";
	archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(world,name.c_str(),get_calc_param().get<int>("nio"));

	RestartMetadata meta;
	meta.current_energy=current_energy;
	meta.spin_restricted=false;
	meta.L=get_calc_param().L();
	meta.k=FunctionDefaults<3>::get_k();
	meta.molecule=mol;
	meta.xc=get_calc_param().xc();
	meta.localize=get_calc_param().localize_method();
	meta.converged_for_thresh=converged_for_thresh;
	meta.converged_for_dconv=converged_for_dconv;
	meta.representation=Representation::znemo;
	meta.ncf=restart_ncf();
	meta.eprec=mol.parameters.eprec();
	meta.madness_version=MADNESS_PACKAGE_VERSION;
	meta.magnetic_field=param.physical_B();
	meta.write(ar);

	auto block=[&](const std::vector<complex_function_3d>& mo, const Tensor<double>& eps) {
		Tensor<double> occ(long(mo.size()));
		occ=1.0;
		ar & (unsigned int) (mo.size());
		ar & eps & occ & std::vector<int>(mo.size(),0);
		for (const auto& f : mo) ar & f;
	};
	block(amo,aeps);
	block(bmo,beps);
}


/// read <prefix>.restartdata; orbitals at another k or a looser thresh are reprojected and no longer converged
void Znemo::load_restartdata() {
	const std::string name=get_calc_param().prefix()+".restartdata";
	archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(world,name.c_str());
	RestartMetadata meta;
	meta.read(ar);
	MADNESS_CHECK_THROW(meta.representation_matches(Representation::znemo),"restartdata does not hold znemo orbitals");
	MADNESS_CHECK_THROW(meta.L==get_calc_param().L(),"inconsistent box size in restartdata");
	current_energy=meta.current_energy;
	converged_for_thresh=meta.converged_for_thresh;
	converged_for_dconv=meta.converged_for_dconv;

	MolecularOrbitals<double_complex,3> a, b;
	a.load_mos(ar,mol,get_calc_param().nalpha());
	b.load_mos(ar,mol,get_calc_param().nbeta());
	amo=a.get_mos();
	bmo=b.get_mos();
	aeps=copy(a.get_eps());
	beps=copy(b.get_eps());

	const int k=FunctionDefaults<3>::get_k();
	const double thresh=FunctionDefaults<3>::get_thresh();
	bool redo=false;
	for (auto* mo : {&amo,&bmo}) {
		for (auto& f : *mo) {
			if (f.k()!=k) {
				f=project(f,k,thresh,false);
				redo=true;
			} else if (f.thresh()*0.999>thresh) {
				redo=true;
			}
		}
	}
	world.gop.fence();
	if (redo) {
		if (world.rank()==0) print("restartdata: reprojected to k",k,"thresh",thresh);
		converged_for_thresh=converged_for_dconv=1.e10;
	}
}





namespace {
/// phi(r) on center A, magnetically translated (London phase e^{-i(B/2)(A_x y - A_y x)}) and Landau-squeezed
/// (e^{-(|B|/8) rho_A^2}); in the Landau gauge times e^{+i(B/2) x y}
class magnetic_ao_functor : public FunctionFunctorInterface<double_complex,3> {
	std::function<double_complex(const coord_3d&)> phi;
	coord_3d A;
	double B;
	bool landau;
public:
	magnetic_ao_functor(std::function<double_complex(const coord_3d&)> phi, const coord_3d& A, double B, bool landau)
		: phi(phi), A(A), B(B), landau(landau) {}
	double_complex operator()(const coord_3d& r) const {
		const double dx=r[0]-A[0], dy=r[1]-A[1];
		double_complex out=phi(r)*std::exp(-0.125*std::abs(B)*(dx*dx+dy*dy))
				*std::exp(double_complex(0.0,-0.5*B*(A[0]*r[1]-A[1]*r[0])));
		if (landau) out*=std::exp(double_complex(0.0,0.5*B*r[0]*r[1]));
		return out;
	}
	std::vector<coord_3d> special_points() const {return {A};}
};
}

std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> >
Znemo::hcore_guess() const {

	// basis: the AOs of aobasis and, in a field, the s and p_{-sgn B} Slater functions e^{-kappa r},
	// (x - i sgn(B) y) e^{-kappa r} of each nucleus, kappa = sqrt(2) Z times each guess_slater_scale (the field
	// compresses the atom along B: with sqrt(2) Z alone the guess is too diffuse at |B| ~ 100, its orbital
	// energies lie above |B|/2 and pi_u/pi_g of H2 are not resolved); all magnetically translated and
	// Landau-squeezed.
	// The Hamiltonian is the one-electron part of the Fock matrix of the solver (kernel gauge).
	const double Bz=param.physical_B();
	const double sgn=(Bz<0.0) ? -1.0 : 1.0;
	const bool landau=mkernel->gauge()==MagneticGauge::landau;
	std::vector<std::function<double_complex(const coord_3d&)> > phi;
	std::vector<coord_3d> center;
	for (int i=0; i<aobasis.nbf(mol); ++i) {
		const AtomicBasisFunction ao=aobasis.get_atomic_basis_function(mol,i);
		phi.push_back([ao](const coord_3d& r) {return double_complex(ao(r[0],r[1],r[2]),0.0);});
		center.push_back(ao.get_coords_vec());
	}
	if (Bz!=0.0) {
		for (const Atom& a : mol.get_atoms()) {
			const coord_3d A=a.get_coords();
			for (const double scale : param.guess_slater_scale()) {
				const double kap=std::sqrt(2.0)*a.q*scale;
				auto s=[A,kap](const coord_3d& r) {return double_complex(std::exp(-kap*(r-A).normf()),0.0);};
				phi.push_back(s);
				phi.push_back([A,sgn,s](const coord_3d& r) {return double_complex(r[0]-A[0],-sgn*(r[1]-A[1]))*s(r);});
				center.push_back(A);
				center.push_back(A);
			}
		}
	}
	auto project=[&](const bool landau_phase) {
		std::vector<complex_function_3d> f(phi.size());
		for (std::size_t i=0; i<phi.size(); ++i) {
			f[i]=complex_factory_3d(world).functor(std::shared_ptr<FunctionFunctorInterface<double_complex,3> >(
					new magnetic_ao_functor(phi[i],center[i],Bz,landau_phase)));
			f[i].scale(1.0/f[i].norm2());
		}
		return truncate(f);
	};
	std::vector<complex_function_3d> aos=project(landau);
	const long nbf=aos.size();
	std::vector<complex_function_3d> nemo_aos=(ncf->type()==NuclearCorrelationFactor::None)
			? aos : truncate(aos*ncf->inverse());

	// one-electron Hamiltonian and overlap as in the Fock matrix of the iteration
	std::vector<complex_function_3d> bra=make_bra(nemo_aos);
	Tensor<double_complex> h=kinetic_pi_matrix(nemo_aos)+compute_vmat(nemo_aos,one_electron_potentials(nemo_aos));
	Tensor<double_complex> S=matrix_inner(world,bra,nemo_aos);
	h=0.5*(h+conj_transpose(h));
	S=0.5*(S+conj_transpose(S));

	// canonical orthogonalization (the Slater functions overlap with the AOs)
	auto [s,U]=syev(S);
	std::vector<long> keep;
	for (long i=0; i<nbf; ++i) if (s(i)>1.e-6*s(nbf-1)) keep.push_back(i);
	Tensor<double_complex> X(nbf,long(keep.size()));
	for (std::size_t j=0; j<keep.size(); ++j) X(_,long(j))=U(_,keep[j])*(1.0/std::sqrt(s(keep[j])));
	auto [e,Cp]=syev(Tensor<double_complex>(inner(conj_transpose(X),inner(h,X))));
	const Tensor<double_complex> C=inner(X,Cp);
	const long nmo=C.dim(1);

	// <L_z> of the eigenvectors (symmetric gauge, on the orbitals)
	Tensor<double> lzexp(nmo);
	{
		std::vector<complex_function_3d> sym_aos=landau ? project(false) : aos;
		Lz<double_complex,3> lz(world,false);
		const Tensor<double_complex> Lmat=matrix_inner(world,sym_aos,lz(sym_aos));
		const Tensor<double_complex> Smat=matrix_inner(world,sym_aos,sym_aos);
		for (long k=0; k<nmo; ++k) {
			const Tensor<double_complex> c=copy(C(_,k)), cc=conj(c);
			lzexp(k)=real(cc.trace(inner(Lmat,c)))/real(cc.trace(inner(Smat,c)));
		}
	}

	// inversion parity <P> of the eigenvectors (centrosymmetric molecules only)
	const bool have_parity=centrosymmetric(mol);
	Tensor<double> parity(nmo);
	if (have_parity) {
		std::vector<complex_function_3d> paos(aos.size());
		for (std::size_t i=0; i<aos.size(); ++i) paos[i]=mirror(aos[i],std::vector<long>{-1,-1,-1});
		const Tensor<double_complex> Pmat=matrix_inner(world,aos,paos);
		const Tensor<double_complex> Smat=matrix_inner(world,aos,aos);
		for (long k=0; k<nmo; ++k) {
			const Tensor<double_complex> c=copy(C(_,k)), cc=conj(c);
			parity(k)=real(cc.trace(inner(Pmat,c)))/real(cc.trace(inner(Smat,c)));
		}
	}

	// occupied orbitals: aufbau, or the lowest with the requested L_z
	const long nocc=std::max(get_calc_param().nalpha(),get_calc_param().nbeta());
	MADNESS_CHECK_THROW(nocc<=nmo,"hcore guess: more occupied orbitals than basis functions");
	const std::vector<int> want=param.guess_lz();
	std::vector<long> pick;
	if (want.empty()) {
		for (long k=0; k<nocc; ++k) pick.push_back(k);
	} else {
		MADNESS_CHECK_THROW(long(want.size())==nocc,"hcore guess: guess_lz needs one L_z per occupied orbital");
		for (int m : want) {
			long k=0;
			for (; k<nmo; ++k) if (std::abs(lzexp(k)-m)<0.25 and std::find(pick.begin(),pick.end(),k)==pick.end()) break;
			MADNESS_CHECK_THROW(k<nmo,"hcore guess: no orbital with the requested L_z (is B along the molecular axis?)");
			pick.push_back(k);
		}
	}
	if (world.rank()==0) {
		print("hcore guess:",nbf,"functions,",nmo,"after canonical orthogonalization (aobasis",get_calc_param().aobasis()+")");
		print(have_parity ? "     k      energy        <Lz>       <P>" : "     k      energy        <Lz>");
		for (long k=0; k<std::min(nmo,nocc+6); ++k) {
			const bool occ=std::find(pick.begin(),pick.end(),k)!=pick.end();
			if (have_parity) printf("  %4ld  %12.6f  %8.4f  %8.4f  %s\n",k,e(k),lzexp(k),parity(k),occ ? "occupied" : "");
			else printf("  %4ld  %12.6f  %8.4f  %s\n",k,e(k),lzexp(k),occ ? "occupied" : "");
		}
	}

	std::vector<complex_function_3d> mo;
	for (long k : pick) mo.push_back(transform(world,nemo_aos,copy(C(_,Slice(k,k))))[0]);
	mo=truncate(mo);
	std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> > zmos;
	zmos.first.set_mos(std::vector<complex_function_3d>(mo.begin(),mo.begin()+get_calc_param().nalpha()));
	zmos.second.set_mos(std::vector<complex_function_3d>(mo.begin(),mo.begin()+get_calc_param().nbeta()));
	return zmos;
}

Znemo::potentials Znemo::one_electron_potentials(const std::vector<complex_function_3d>& rhs) const {
	Nuclear<double_complex,3> nuc(world,ncf);

	potentials pot(world,rhs.size());
	pot.vnuc_mo=nuc(rhs);

	// the paramagnetic and diamagnetic terms enter through kinetic_pi_matrix (and, for a kernel without the
	// paramagnetic term, through assemble_rhs)
	if (kernel_nemo()) pot.vbare_mo=truncate(potentialmanager->vnuclear()*rhs);

	// i U1.A F: symmetric gauge about r = 0 i (B/2)(x U1_y - y U1_x); Landau gauge of the kernel (A = (-B y, 0, 0))
	// -i B y U1_x
	if (mkernel->gauge()==MagneticGauge::landau)
		pot.zeeman_R_comm = double_complex(0.0,-1.0)*B[2]*(rvec[1]*ncf->U1(0))*rhs;
	else
		pot.zeeman_R_comm = double_complex(0.0,0.5)*B[2]*(rvec[0]* ncf->U1(1)-rvec[1]* ncf->U1(0))*rhs;

	truncate(world,pot.vnuc_mo);
	return pot;
}

/// compute the potential operators applied on the orbitals
Znemo::potentials Znemo::compute_potentials(const std::vector<complex_function_3d>& mo,
		const real_function_3d& density,
		const std::vector<complex_function_3d>& rhs) const {

	timer potential_timer(world,print_info.print_timings());
	std::vector<complex_function_3d> dia2mo=make_bra(mo);

	// prepare exchange operator
	Exchange<double_complex,3> K(world,get_calc_param().lo());
	Tensor<double> occ(mo.size());
	occ=1.0;
    K.set_bra_and_ket(conj(world, dia2mo), mo);

	potentials pot=one_electron_potentials(rhs);
	pot.J_mo=(*coulop)(density)*rhs;
	pot.K_mo=K(rhs);

	truncate(world,pot.J_mo);
	truncate(world,pot.K_mo);
	truncate(world,pot.vnuc_mo);
	potential_timer.end("compute_potentials");
	return pot;
};


Tensor<double_complex> Znemo::compute_vmat(const std::vector<complex_function_3d>& mo,
		const potentials& pot) const {

	std::vector<complex_function_3d> dia2mo=make_bra(mo);
	truncate(world,dia2mo);
	// kernels with a correlation factor: the symmetric kinetic term contains U1.grad + i U1.A + U2 - V
	const bool sym=kernel_nemo();
	Tensor<double_complex> Vnucmat=matrix_inner(world,dia2mo,sym ? pot.vbare_mo : pot.vnuc_mo);
	Tensor<double_complex> Kmat=matrix_inner(world,dia2mo,pot.K_mo);
	Tensor<double_complex> Jmat=matrix_inner(world,dia2mo,pot.J_mo);
	Tensor<double_complex> ZRcom_mat=sym ? Tensor<double_complex>(mo.size(),mo.size()) : matrix_inner(world,dia2mo,pot.zeeman_R_comm);

	Tensor<double_complex> vmat=Vnucmat-Kmat+Jmat+ZRcom_mat;
	return vmat;
};


std::vector<complex_function_3d>
Znemo::compute_residuals(
		const std::vector<complex_function_3d>& Vpsi,
		const std::vector<complex_function_3d>& psi,
		Tensor<double>& eps) const {
	timer residual_timer(world,print_info.print_timings());

	double tol = FunctionDefaults < 3 > ::get_thresh();

    // tmp = -(T_k - E)^{-1} Vpsi
    // G(E) needs E < |B|/2, and its range along B grows as 1/sqrt(2(|B|/2 - E)) (and with it the memory: at
    // B = 10, E = B/2 - 1e-3 exceeded 45 GB). An energy bound by at least 1e-3 is used as it is; one closer to
    // the threshold is no bound-state estimate (a guess) and is replaced by |B|/2 - 0.05 (decay length 3 bohr)
    const double bottom=mkernel->bottom();
    Tensor<double> E(eps.size());
    for (int i=0; i<eps.size(); ++i) E(i)=(eps(i)<bottom-1.e-3) ? eps(i) : bottom-0.05;
    if ((world.rank()==0) and (eps.size()>0) and (eps.max()>=bottom-1.e-3))
    	print("kernel energies at the Landau threshold |B|/2 replaced by |B|/2 - 0.05:",eps,"->",E);
    for (int i=0; i<eps.size(); ++i) eps(i)=E(i);	// the update below starts from the energy used
    if ((world.rank()==0) and (print_info.print_convergence())) print("kernel energies",E);
    const std::vector<complex_function_3d> tmp=mkernel->apply(-1.0*Vpsi,E);

    // truncating the residual at thresh zeroes it (and the energy update below) before the fixed point is
    // reached, at strong field by up to ~1e-3 Eh in a gauge-dependent way
    const std::vector<complex_function_3d> res=truncate(psi-tmp,0.1*tol);
    const std::vector<complex_function_3d> bra_res=make_bra(res);

    // update eps
    Tensor<double> norms=real(inner(world,make_bra(tmp),tmp));
    Tensor<double_complex> rnorms=inner(world,bra_res,res);
    for (int i=0; i<norms.size(); ++i) {
    	norms(i)=sqrt(norms(i));
    	rnorms(i)=sqrt(rnorms(i));
    }
    if ((world.rank()==0) and (print_info.print_convergence())) {
    	print("norm2(tmp)",norms);
    	print("norm2(res)",rnorms);
    }
    Tensor<double_complex> in=inner(world,Vpsi,bra_res);	// no shift here!
    Tensor<double> delta_eps(eps.size());
    for (int i=0; i<eps.size(); ++i) delta_eps(i)=real(in(i))/(norms[i]*norms[i]);
    eps-=delta_eps;

    if ((world.rank()==0) and (print_info.print_convergence())) {
    	print("orbital energy update",delta_eps);
    }
    residual_timer.end("compute_residuals");
    return res;

}


void
Znemo::canonicalize(std::vector<complex_function_3d>& amo,
		std::vector<complex_function_3d>& vnemo,
		potentials& pot,
		XNonlinearSolver<std::vector<complex_function_3d> ,double_complex, vector_function_allocator<double_complex,3> >& solver,
		Tensor<double_complex> fock, Tensor<double_complex> ovlp) const {

    Tensor<double_complex> U;
    Tensor<double> evals;
    sygv(fock, ovlp, 1, U, evals);
    // Fix phases.
    for (size_t i = 0; i < amo.size(); ++i) if (real(U(i, i)) < 0.0) U(_, i).scale(-1.0);

    fock = 0.0;
    for (unsigned int i = 0; i < amo.size(); ++i) fock(i, i) = evals(i);
    amo = transform(world, amo, U);
    vnemo = transform(world, vnemo, U);
    rotate_subspace(U,solver);
    pot.transform(U);

}

void Znemo::save_orbitals(std::string suffix) const
{
	suffix = "_" + suffix;
	for (size_t i = 0; i < amo.size(); ++i)
		save(amo[i], "amo" + stringify(i) + suffix);
	for (size_t i = 0; i < bmo.size(); ++i)
		save(bmo[i], "bmo" + stringify(i) + suffix);
	for (size_t i = 0; i < amo.size(); ++i)
		save(real(amo[i]), "amo_real" + stringify(i) + suffix);
	for (size_t i = 0; i < amo.size(); ++i)
		save(imag(amo[i]), "amo_imag" + stringify(i) + suffix);
	for (size_t i = 0; i < bmo.size(); ++i)
		save(real(bmo[i]), "bmo_real" + stringify(i) + suffix);
	for (size_t i = 0; i < bmo.size(); ++i)
		save(imag(bmo[i]), "bmo_imag" + stringify(i) + suffix);
	for (size_t i = 0; i < amo.size(); ++i)
		save(madness::abs(amo[i]), "absamo" + stringify(i) + suffix);
	for (size_t i = 0; i < bmo.size(); ++i)
		save(madness::abs(bmo[i]), "absbmo" + stringify(i) + suffix);
}

} // namespace madness
