/*
 This file is part of MADNESS.

 Copyright (C) 2007,2010 Oak Ridge National Laboratory

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

 For more information please contact:

 Robert J. Harrison
 Oak Ridge National Laboratory
 One Bethel Valley Road
 P.O. Box 2008, MS-6367

 email: harrisonrj@ornl.gov
 tel:   865-241-3937
 fax:   865-572-0680
 */

/*!
 \file examples/nemo.h
 \brief solve the HF equations using numerical exponential MOs

 The source is
 <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
 /trunk/src/apps/examples/nemo.h>here</a>.

 */

#ifndef NEMO_H_
#define NEMO_H_

#include <madness/mra/mra.h>
#include <madness/mra/funcplot.h>
#include <madness/mra/operator.h>
#include <madness/mra/lbdeux.h>
#include<madness/chem/SCF.h>
#include<madness/chem/CalculationParameters.h>
#include<madness/chem/SCFProtocol.h>
#include<madness/chem/correlationfactor.h>
#include<madness/chem/molecular_optimizer.h>
#include <madness/mra/nonlinsol.h>
#include <madness/mra/vmra.h>
#include<madness/chem/pcm.h>
#include<madness/chem/AC.h>
#include<madness/chem/pointgroupsymmetry.h>
#include"madness/mra/commandlineparser.h"
#include<madness/chem/QCPropertyInterface.h>
#include <madness/world/timing_utilities.h>

namespace madness {

class PNO;
class OEP;


class NemoBase : public MolecularOptimizationTargetInterface {

public:

	NemoBase(World& w) : world(w) {}

    virtual ~NemoBase() {}

    virtual std::shared_ptr<Fock<double,3>> make_fock_operator() const {
	    MADNESS_EXCEPTION("implement make_fock operator for your derived NemoBase class",1);
	    return std::shared_ptr<Fock<double,3>>();
	}

        /// create an instance of the derived object based on the input parameters
	std::shared_ptr<NuclearCorrelationFactor> get_ncf_ptr() const {
		return ncf;
	}

	/// normalize the nemos
	template<typename T, std::size_t NDIM>
	void static normalize(std::vector<Function<T,NDIM> >& nemo,
			const Function<double,NDIM> metric=Function<double,NDIM>()) {

        if (nemo.size()==0) return;
        World& world=nemo[0].world();
		// compute the norm of the reconstructed orbitals, includes the factor
		std::vector<Function<T,NDIM> > mos = (metric.is_initialized()) ? metric*nemo : nemo;
		std::vector<double> norms = norm2s(world, mos);

		// scale the nemos, excludes the nuclear correlation factor
		std::vector<double> invnorm(norms.size());
		for (std::size_t i = 0; i < norms.size(); ++i)
			invnorm[i] = 1.0 / norms[i];
		scale(world, nemo, invnorm);
	}

	template<typename T>
    static Tensor<T> Q2(const Tensor<T>& s) {
		Tensor<T> Q = -0.5*s;
        for (int i=0; i<s.dim(0); ++i) Q(i,i) += 1.5;
        return Q;
    }

	/// orthonormalize the vectors
	template<typename T, std::size_t NDIM>
	void orthonormalize(std::vector<Function<T,NDIM> >& nemo,
			const Function<double,NDIM> metric=Function<double,NDIM>(),
			const double trantol=FunctionDefaults<NDIM>::get_thresh()*0.01) const {

		if (nemo.size()==0) return;
	    normalize(nemo,metric);
	    double maxq;
	    do {
			std::vector<Function<T,NDIM> > Rnemo = (metric.is_initialized()) ? metric*nemo : nemo;
	        Tensor<T> Q = Q2(matrix_inner(world, Rnemo, Rnemo));
	        maxq=0.0;
	        for (int i=0; i<Q.dim(0); ++i)
	            for (int j=0; j<i; ++j)
	                maxq = std::max(maxq,std::abs(Q(i,j)));

	        Q.screen(trantol); // ???? Is this really needed?
	        nemo = transform(world, nemo, Q, trantol, true);
	        truncate(world, nemo);
//	        if (world.rank() == 0) print("ORTHOG2: maxq trantol", maxq, trantol);

	    } while (maxq>0.01);
	    normalize(nemo,metric);
	}

	template<typename T, std::size_t NDIM>
	Function<typename Tensor<T>::scalar_type,NDIM> compute_density(const std::vector<Function<T,NDIM> > nemo) const {
		return sum(world,abssq(world,nemo)).truncate();
	}

    virtual bool need_recompute_factors_and_potentials(const double thresh) const {
        bool need=false;
        if ((not R.is_initialized()) or (R.thresh()>thresh)) need=true;
        if (not ncf) need=true;
        if ((not R_square.is_initialized()) or (R_square.thresh()>thresh)) need=true;
        return need;
    };

    virtual void invalidate_factors_and_potentials() {
        R.clear();
        R_square.clear();
        ncf.reset();
    };

    void construct_nuclear_correlation_factor(const Molecule& molecule,
			const std::shared_ptr<PotentialManager> pm,
			const std::pair<std::string,double> ncf_parameter) {

	    // construct the nuclear correlation factor:
	    if (not ncf) {
	    	ncf=create_nuclear_correlation_factor(world, molecule, pm, ncf_parameter);
	    }

	    // re-project the ncf
	    ncf->initialize(FunctionDefaults<3>::get_thresh());
	    R = ncf->function();
	    R.set_thresh(FunctionDefaults<3>::get_thresh());
	    R_square = ncf->square();
	    R_square.set_thresh(FunctionDefaults<3>::get_thresh());
	}

	/// compute the nuclear gradients
	Tensor<double> compute_gradient(const real_function_3d& rhonemo,
			const Molecule& molecule) const;

    /// compute kinetic energy as square of the "analytical" expectation value

	/// @param[in]	the nemo orbitals F
	/// @return		T = 1/2 \sum_i \int R^2 U1.U1 F^2 + 2 R^2 U1.grad(F) + R^2 grad(F)^2
	template<typename T, std::size_t NDIM>
	double compute_kinetic_energy(const std::vector<Function<T,NDIM> >& nemo) const {

		// T = 0.5\sum_i \int R^2 U1.U1 F^2 - 2 R^2 U1.grad(F) F + R^2 grad(F)^2
		//   = 0.5 (<U1.U1 | rho > + <R^2|grad(F)^2> - 2<R^2 | U1.grad(F) >)
		// note: U1=-grad(R)/R
		//auto id=nemo.front().world().id();
		//auto id1=R_square.world().id();
		//auto worldid=world.id();
		world.gop.fence();
		real_function_3d dens=dot(world,nemo,nemo)*R_square;
	    real_function_3d U1dotU1=real_factory_3d(world)
	    		.functor(NuclearCorrelationFactor::U1_dot_U1_functor(ncf.get()));
	    double ke1=inner(dens,U1dotU1);

	    double ke2=0.0;
	    double ke3=0.0;
	    //double ke3_real=0.0;
	    //double ke3_imag=0.0;

	    for (size_t axis = 0; axis < NDIM; axis++) {
	        real_derivative_3d D = free_space_derivative<double, NDIM>(world, axis);
	        const std::vector<Function<T,NDIM> > dnemo = apply(world, D, nemo);

	        real_function_3d term2=dot(world,dnemo,nemo)*ncf->U1(axis);
	        double tmp=-2.0*inner(R_square,term2);
//	        ke3_real -=2.0*std::real(tmp);
//	        ke3_imag -=2.0*std::imag(tmp);
	        ke3 +=tmp;

	        const real_function_3d term1=dot(world,dnemo,dnemo);
	        world.gop.fence();
	        ke2 += inner(term1,R_square);

	    }
//	    if (ke3_imag>1.e-8) {
//	    	print("kinetic energy, imaginary part: ",ke3_imag);
//	    	MADNESS_EXCEPTION("imaginary kinetic energy",1);
//	    }
//	    double ke=2.0*(ke1+ke2+ke3_real); // closed shell
	    double ke=2.0*(ke1+ke2+ke3); // closed shell
	    return 0.5*ke;
	}

    /// compute kinetic energy as square of the "analytical" derivative of the orbitals

	/// @param[in]	the nemo orbitals F
	/// @return		T = 1/2 \sum_i || grad(R)*F_i + R*grad(F_i)||^2
	template<typename T, std::size_t NDIM>
	double compute_kinetic_energy1(const std::vector<Function<T,NDIM> >& nemo) const {
        timer timer1(world);
		double ke=0.0;
		for (int i=0; i<nemo.size(); ++i) {
			double fnorm2=norm2(world,-1.0*R*ncf->U1vec()*nemo[i] + R*grad(nemo[i]));
			ke+=2.0*fnorm2*fnorm2;
		}
		timer1.end("compute_kinetic_energy1");
		return 0.5*ke;
	}


    /// compute kinetic energy as square of the "analytical" derivative of the orbitals

	/// @param[in]	the nemo orbitals F
	/// @return		T = 1/2 \sum_i || grad(R)*F_i + R*grad(F_i)||^2
	template<typename T, std::size_t NDIM>
	double compute_kinetic_energy1a(const std::vector<Function<T,NDIM> >& nemo) const {
        timer timer1(world);
		double ke=0.0;
		for (int i=0; i<NDIM; ++i) {
	        std::vector< std::shared_ptr< Derivative<T,NDIM> > > grad=
	                gradient_operator<T,NDIM>(world);
			double fnorm2=norm2(world,R*(-1.0*ncf->U1(i)*nemo + apply(world,*(grad[i]),nemo)));
			ke+=2.0*fnorm2*fnorm2;
		}
		timer1.end("compute_kinetic_energy1a");
		return 0.5*ke;
	}

    /// compute kinetic energy as direct derivative of the orbitals (probably imprecise)

	/// @param[in]	the nemo orbitals F
	/// @return		T = 1/2 \sum_i || grad(R*F_i)||^2
	template<typename T, std::size_t NDIM>
    double compute_kinetic_energy2(const std::vector<Function<T,NDIM> >& nemo) const {

    	// it's ok to use phi here, no regularization necessary for this eigenvalue
    	double E_kin = 0.0;
    	for (int axis = 0; axis < 3; axis++) {
    		real_derivative_3d D = free_space_derivative<double, 3>(world, axis);
    		const vecfuncT dphi = apply(world, D,R*nemo);
    		E_kin += 0.5 * (inner(world, dphi, dphi)).sum();
    		// -1/2 sum <Psi|Nabla^2|Psi> = 1/2 sum <NablaPsi|NablaPsi>   (integration by parts)
    	}
    	E_kin *= 2.0; // 2 because closed shell
    	return E_kin;
    }


	bool check_convergence(const std::vector<double> energies,
			const std::vector<double> oldenergies, const double bsh_norm,
			const double delta_density, const CalculationParameters& param,
			const double econv, const double dconv) const {

        double maxenergychange=fabs(energies.size()-oldenergies.size());	// >0 if oldenergyvec not initialized
        for (auto iter1=energies.begin(), iter2=oldenergies.begin();
        		(iter1!=energies.end() and iter2!=oldenergies.end()); iter1++, iter2++) {
        	maxenergychange=std::max(maxenergychange,fabs(*iter1 - *iter2));
        }
        double delta_energy=fabs(energies[0]-oldenergies[0]);

		bool bsh_conv=param.converge_bsh_residual() ? bsh_norm<dconv : true;
		bool total_energy_conv=param.converge_total_energy() ? delta_energy<econv : true;
		bool each_energy_conv=param.converge_each_energy() ? maxenergychange<econv*3.0 : true;
		bool density_conv=param.converge_density() ? delta_density<dconv : true;

		if (world.rank()==0 and param.print_level()>2) {
			std::stringstream line;
			line << "convergence: bshresidual, energy change, max energy change, density change "
					<< std::scientific << std::setprecision(1)
					<< bsh_norm << " " << delta_energy << " "
					<< maxenergychange << " " << delta_density;
			print(line.str());
		}

		return (bsh_conv and density_conv and each_energy_conv and total_energy_conv);
	}

	World& world;

	/// the nuclear correlation factor
	std::shared_ptr<NuclearCorrelationFactor> ncf;

	/// the nuclear correlation factor
	real_function_3d R;

    /// the square of the nuclear correlation factor
    real_function_3d R_square;


};


/// The Nemo class
class Nemo: public NemoBase, public QCPropertyInterface {
	typedef std::shared_ptr<real_convolution_3d> poperatorT;
	friend class PNO;
	friend class TDHF;

public:
	/// class holding parameters for a nemo calculation beyond the standard dft parameters from moldft
	struct NemoCalculationParameters : public QCCalculationParametersBase {
		static constexpr char const* tag = "nemo";

        NemoCalculationParameters(World& world, const commandlineparser& parser) {
            initialize_nemo_parameters();
        	read_input_and_commandline_options(world,parser,tag);
        }

		NemoCalculationParameters() {
			initialize_nemo_parameters();
		}

		std::string get_tag() const override {
        	return std::string("dft");
        }

		void initialize_nemo_parameters() {
        	// check if parameters are initialized for a nemo calculation already
        	if (parameter_exists("ncf")) return;
			initialize<std::pair<std::string,double> > ("ncf",{"slater",2.0},"nuclear correlation factor");
			initialize<bool> ("hessian",false,"compute the hessian matrix");
			initialize<bool> ("read_cphf",false,"read the converged orbital response for nuclear displacements from file");
			initialize<bool> ("restart_cphf",false,"read the guess orbital response for nuclear displacements from file");
			initialize<bool> ("purify_hessian",false,"symmetrize the hessian matrix based on atomic charges");
		}

		std::pair<std::string,double> ncf() const {return get<std::pair<std::string,double> >("ncf");}
		bool hessian() const {return get<bool>("hessian");}

	};


public:
	std::filesystem::path work_dir;

	/// ctor

	/// @param[in]	world1	the world
	/// @param[in]	calc	the SCF
//	Nemo(World& world1, std::shared_ptr<SCF> calc, const std::string inputfile);

    Nemo(World& world, const commandlineparser& parser);

	Nemo(World& world, const CalculationParameters& param, const NemoCalculationParameters& nemo_param,
		const Molecule& molecule);

    std::string name() const {return "nemo";}
    bool selftest() {return false;}

    static void help() {
        print_header2("help page for NEMO");
        print("The nemo code computes Hartree-Fock and DFT energies, gradients and hessians using a nuclear correlation factor");
        print("that regularizes the singular nuclear potential. SCF orbitals for the basis for post-SCF calculations like");
        print("excitation energies (cis), correlation energies (cc2), local potentials (oep), etc\n");
        print("A nemo calculation input is mostly identical to a moldft calculation input, but it uses the additional input");
        print("parameter ncf (nuclear correlation factor)\n");
        print("You can print all available calculation parameters by running\n");
        print("nemo --print_parameters\n");
        print("You can perform a simple calculation by running\n");
        print("nemo --geometry=h2o.xyz\n");
        print("provided you have an xyz file in your directory.");

    }

    static void print_parameters() {
        NemoCalculationParameters param;
        print("default parameters for the nemo program are");
        param.print("nemo","end");
        print("\n\nthe molecular geometry must be specified in a separate block:");
        Molecule::print_parameters();
    }


	bool check_converged(const Tensor<double>& x) const {
    	double xsq = x.sumsq();
    	return (xsq == coords_sum);
	}

	virtual double value() {return value(calc->molecule.get_all_coords());}

	virtual double value(const Tensor<double>& x);

	void load_mos(World& w) {
		calc->load_mos(w);
	}

	/// compute dipole moment and gradient at the current geometry
	virtual nlohmann::json analyze() const;

	/// compute the nuclear gradients
	Tensor<double> gradient(const Tensor<double>& x);

	bool provides_gradient() const {return true;}

	/// returns the molecular hessian matrix at structure x
	Tensor<double> hessian(const Tensor<double>& x);

	/// construct the fock operator based on the calculation parameters (K or XC?)
	virtual std::shared_ptr<Fock<double,3>> make_fock_operator() const;

	/// purify and symmetrize the hessian

	/// The hessian should be symmetric, but it is not, because
	/// \f[
	///  \langle i^{Y_B}|H^{X_A}|i\rangle \neq \langle i|H^{X_A}|i^{Y_B}\rangle
	/// \f]
	/// does holds analytically, but not numerically. If the two numbers
	/// differ, pick the more trustworthy, which is the one with a heavy
	/// atom causing the perturbed density and the light atom being the
	/// nuclear singularity.
	/// @param[in]  hessian the raw hessian
	/// @return     a symmetrized hessian
	Tensor<double> purify_hessian(const Tensor<double>& hessian) const;


	/// solve the CPHF equations for the nuclear displacements

	/// this function computes that part of the orbital response that is
	/// orthogonal to the occupied space. If no NCF's are used this
	/// corresponds to the normal response. If NCF's are used the part
	/// parallel to the occupied space must be added!
	/// \f[
	///     F^X = F^\perp + F^\parallel
	/// \f]
	/// cf parallel_CPHF()
	/// @param[in]  iatom   the atom A to be moved
	/// @param[in]  iaxis   the coordinate X of iatom to be moved
	/// @return     \ket{i^X} or \ket{F^\perp}
	vecfuncT solve_cphf(const size_t iatom, const int iaxis, const Tensor<double> fock,
	        const vecfuncT& guess, const vecfuncT& rhsconst,
	        const Tensor<double> incomplete_hessian, const vecfuncT& parallel,
	        const SCFProtocol& p, const std::string& xc_data) const;

	/// solve the CPHF equation for all displacements

	/// this function computes the nemo response F^X
    /// \f[
    ///     F^X = F^\perp + F^\parallel
    /// \f]
	/// To reconstruct the unregularized orbital response (not recommended):
	/// \f[
	///   i^X   = R^X F + R F^X
	/// \f]
	/// The orbital response i^X constructed in this way is automatically
	/// orthogonal to the occupied space because of the parallel term F^\parallel
	/// @return a vector of the nemo response F^X for all displacements
	std::vector<vecfuncT> compute_all_cphf();

    /// this function computes that part of the orbital response that is
    /// parallel to the occupied space.
    /// \f[
    ///     F^X = F^\perp + F^\parallel
    /// \f]
    /// If no NCF's are used F^\parallel vanishes.
    /// If NCF's are used this term does not vanish because the derivatives of
    /// the NCF does not vanish, and it is given by
    /// \f[
    ///  F_i^\parallel = -\frac{1}{2}\sum_k|F_k ><F_k | (R^2)^X | F_i>
    /// \f]
    vecfuncT compute_cphf_parallel_term(const size_t iatom, const int iaxis) const;

    /// compute the IR intensities in the double harmonic approximation

    /// use the projected normal modes; units are km/mol
    /// @param[in]  normalmodes the normal modes
    /// @param[in]  dens_pt the perturbed densities for each nuclear displacement
    Tensor<double> compute_IR_intensities(const Tensor<double>& normalmodes,
            const vecfuncT& dens_pt) const;

	std::shared_ptr<SCF> get_calc() const {return calc;}

    const NemoCalculationParameters& get_nemo_param() const {return nemo_param;}
	const CalculationParameters& get_calc_param() const {return calc->param;}

	PCM get_pcm()const{return pcm;}

	/// compute the Fock matrix from scratch
	tensorT compute_fock_matrix(const vecfuncT& nemo, const tensorT& occ) const;

	/// return a reference to the molecule
	Molecule& molecule() {return calc->molecule;}

    /// return a reference to the molecule
    Molecule& molecule() const {
        return calc->molecule;
    }

    /// make the density (alpha or beta)
    real_function_3d make_density(const Tensor<double>& occ,
            const vecfuncT& nemo) const;

    /// make the density using different bra and ket vectors

    /// e.g. for computing the perturbed density \sum_i \phi_i \phi_i^X
    /// or when using nemos: \sum_i R2nemo_i nemo_i
    real_function_3d make_density(const tensorT & occ,
            const vecfuncT& bra, const vecfuncT& ket, const bool refine=false) const;

    /// make the derivative of the density

    /// \f$ \nabla\rho = 2R^X R \rho_R + R^2\nabla \rho_R \f$
    /// @param[in]  rhonemo    the regularized density
    /// @param[in]  axis       the component of the nabla operator
    /// @return     the gradient of the *reconstructed* density
    real_function_3d make_ddensity(const real_function_3d& rhonemo,
            const int axis) const;

    /// compute the reduced densities sigma (gamma) for GGA functionals
    real_function_3d make_sigma(const real_function_3d& rho1,
            const real_function_3d& rho2) const;


    /// the Laplacian of the density

    /// The Laplacian should currently only be used for subsequent convolution
    /// with a Green's function (which is reasonably stable), but not on its own!
    ///
    /// The Laplacian of the cuspy density is numerically fairly unstable:
    ///  - a singular term may be rewritten using the nuclear potential (see below)
    ///  - the Laplacian of the regularized density is still very noisy
    ///
    /// It may be computed as
    /// \f[
    ///   \Delta \rho = \Delta (R^2 \rho_R)
    ///          = \Delta (R^2) \rho_R + 2\nabla R \nabla \rho_R + R^2 \Delta \rho_R
    ///          = 2 R^2 U1^2 \rho_R -4 R^2 ( U-V ) \rho_R + R^2 \Delta\rho_R
    /// \f]
    /// where we can use the identity
    /// \f[
    ///   U=V + R^{-1}[T,R]
    ///   -2 R (U-V) = \Delta R + 2\nabla R\dot \nabla
    /// \f]
    /// first term comes from the definition of the U potential as the commutator
    /// over the kinetic energy (aka the Laplacian)
    /// @param[in]  rhonemo    the regularized density \rho_R
    /// @return     the laplacian of the reconstructed density \Delta (R^2\rho_R)
    real_function_3d make_laplacian_density(const real_function_3d& rhonemo) const;

    /// compute the kinetic energy potential using Eq. (16) of
    /// R. A. King and N. C. Handy, “Kinetic energy functionals from the Kohn–Sham potential,”
    /// Phys. Chem. Chem. Phys., vol. 2, no. 22, pp. 5049–5056, 2000.
    real_function_3d kinetic_energy_potential(const vecfuncT& nemo) const;


    /// smooth a function by projecting it onto k-1 and then average with k

    /// kept it here for further testing
    static void smoothen(real_function_3d& f) {
        int k=f.get_impl()->get_k();
        real_function_3d fproj=project(f,k-1);
        real_function_3d freproj=project(fproj,k);
        f=0.5*(f+freproj);
    }

protected:

	std::shared_ptr<SCF> calc;

public:
    const NemoCalculationParameters nemo_param;

protected:
    projector_irrep symmetry_projector;

public:

    /// return the symmetry_projector
    projector_irrep get_symmetry_projector() const {
    	return symmetry_projector;
    }

private:

	/// sum of square of coords at last solved geometry
	mutable double coords_sum;

protected:
	/// a poisson solver
	std::shared_ptr<real_convolution_3d> poisson;

	/// asymptotic correction for DFT
	AC<3> ac;

//    /// apply the AC scheme of Tozer/Handy with the multipole approximation
//    Function<double,3> apply_ac(const Function<double,3>& vxc)const{
//    	return ac.apply(vxc);
//    }
//
//    /// apply the AC scheme of Tozer/Handy using the hartree potential
//    Function<double,3> apply_ac(const Function<double,3>& vxc, const Function<double,3>& vhartree)const{
//    	return ac.apply(vxc,vhartree);
//    }
//
//    /// apply the AC scheme of Tozer/Handy
//    Function<double,3> apply_ac_excited(Function<double,3>& vxc, const Function<double,3>& vhartree)const{
//    	return ac.apply_potential(vxc,vhartree);
//    }

private:
	/// polarizable continuum model
	PCM pcm;
//	AC<3> ac;

protected:
	/// adapt the thresholds consistently to a common value
    void set_protocol(const double thresh) {

        calc->set_protocol<3>(world,thresh);

        if (need_recompute_factors_and_potentials(thresh)) {
            timer timer1(world,get_calc_param().print_level()>2);
            get_calc()->make_nuclear_potential(world);
            construct_nuclear_correlation_factor(calc->molecule, calc->potentialmanager, get_nemo_param().ncf());
            timer1.end("reproject ncf");
        }

        // (re) construct the Poisson solver
        poisson = std::shared_ptr<real_convolution_3d>(
                CoulombOperatorPtr(world, get_calc_param().lo(), FunctionDefaults<3>::get_thresh()));

        // set thresholds for the MOs
        set_thresh(world,calc->amo,thresh);
        set_thresh(world,calc->bmo,thresh);

    }

	/// solve the HF equations
	double solve(const SCFProtocol& proto);

    /// given nemos, compute the HF energy using the regularized expressions for T and V
    std::vector<double> compute_energy_regularized(const vecfuncT& nemo, const vecfuncT& Jnemo,
            const vecfuncT& Knemo, const vecfuncT& Unemo) const;

	/// compute the reconstructed orbitals, and all potentials applied on nemo

	/// to use these potentials in the fock matrix computation they must
	/// be multiplied by the nuclear correlation factor
	/// @param[in]	nemo	the nemo orbitals
	/// @param[out]	Jnemo	Coulomb operator applied on the nemos
	/// @param[out]	Knemo	exchange operator applied on the nemos
	/// @param[out]	pcmnemo	PCM (solvent) potential applied on the nemos
	/// @param[out]	Unemo	regularized nuclear potential applied on the nemos
	void compute_nemo_potentials(const vecfuncT& nemo,
			vecfuncT& Jnemo, vecfuncT& Knemo, vecfuncT& xcnemo, vecfuncT& pcmnemo,
			vecfuncT& Unemo) const;

	/// return the Coulomb potential
	real_function_3d get_coulomb_potential(const vecfuncT& psi) const;

	/// compute the incomplete hessian

	/// incomplete hessian is the nuclear-nuclear contribution, and the
	/// contribution from the second derivative of the nuclear potential,
	/// and also the derivative of the nuclear correlation factor.
	/// i.e. all contributions that *do not* contain the regularized perturbed
	/// density, but it will contain parts of the perturbed density
	Tensor<double> make_incomplete_hessian() const;

    /// compute the complementary incomplete hessian

	/// @param[in]  xi the response functions including the parallel part
    Tensor<double> make_incomplete_hessian_response_part(
            const std::vector<vecfuncT>& xi) const;

	/// compute the constant term for the CPHF equations

	/// mainly all terms with the nuclear correlation factor's derivatives
	vecfuncT make_cphf_constant_term(const size_t iatom, const int iaxis,
	        const vecfuncT& R2nemo, const real_function_3d& rhonemo) const;

public:

	bool is_dft() const {return calc->xc.is_dft();}

	bool do_pcm() const {return get_calc_param().pcm_data() != "none";}
	
	bool do_ac() const {return get_calc_param().ac_data() != "none";}

	AC<3> get_ac() const {return ac;}

	bool do_symmetry() const {return (symmetry_projector.get_pointgroup()!="C1");}

protected:

	/// localize the nemo orbitals
    vecfuncT localize(const vecfuncT& nemo, const double dconv, const bool randomize) const;
protected:
	/// return the threshold for vanishing elements in orbital rotations
    double trantol() const {
        return calc->vtol / std::min(30.0, double(get_calc()->amo.size()));
    }

    void make_plots(const real_function_3d &f,const std::string &name="function")const{
        double width = FunctionDefaults<3>::get_cell_min_width()/2.0 - 1.e-3;
        plot_plane(world,f,name);
        coord_3d start(0.0); start[0]=-width;
        coord_3d end(0.0); end[0]=width;
        plot_line(("line_"+name).c_str(),1000,start,end,f);
    }

    /// save a function
    template<typename T, size_t NDIM>
    void save_function(const std::vector<Function<T,NDIM> >& f, const std::string name) const;

    /// load a function
    template<typename T, size_t NDIM>
    void load_function(std::vector<Function<T,NDIM> >& f, const std::string name) const;

};

/// save a function
template<typename T, size_t NDIM>
void Nemo::save_function(const std::vector<Function<T,NDIM> >& f, const std::string name) const {
    if (world.rank()==0) print("saving vector of functions",name);
    archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(world, name.c_str(), 1);
    ar & f.size();
    for (const Function<T,NDIM>& ff:f)  ar & ff;
}

/// load a function
template<typename T, size_t NDIM>
void Nemo::load_function(std::vector<Function<T,NDIM> >& f, const std::string name) const {
    if (world.rank()==0) print("loading vector of functions",name);
    archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(world, name.c_str(), 1);
    std::size_t fsize=0;
    ar & fsize;
    f.resize(fsize);
    for (std::size_t i=0; i<fsize; ++i) ar & f[i];
}

}

#endif /* NEMO_H_ */

