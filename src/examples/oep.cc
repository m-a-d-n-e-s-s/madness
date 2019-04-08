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

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES

/*!
  \file examples/oep.cc
  \brief optimized effective potentials for DFT
*/

#include <chem/nemo.h>
#include <chem/cheminfo.h>
#include <chem/SCFOperators.h>
#include <chem/projector.h>

using namespace madness;

struct dens_inv{

	double threshold;

	dens_inv(const double thresh = 1.0e-8) { // same default value as for dens_thresh
		threshold = thresh;
	}

    /// @param[out] U   result
    /// @param[in]  t   numerator
    /// @param[in]  inv density to be inverted >0.0
    void operator()(const Key<3>& key, Tensor<double>& U, const Tensor<double>& t,
            const Tensor<double>& inv) const {
        ITERATOR(
            U, double d = t(IND);
            double p = std::max(inv(IND), threshold);
            U(IND) = d/p;
        );
   }
    template <typename Archive>
    void serialize(Archive& ar) {}

};


/// TODO: change dens_inv such that Slater potential is munged unsing the following form:

/*
/// Class to compute terms of the potential
struct xc_potential {
    const XCfunctional* xc;
    const int ispin;

    xc_potential(const XCfunctional& xc, int ispin) : xc(&xc), ispin(ispin)
    {}



    std::vector<madness::Tensor<double> > operator()(const madness::Key<3> & key,
            const std::vector< madness::Tensor<double> >& t) const {
        MADNESS_ASSERT(xc);
        std::vector<madness::Tensor<double> > r = xc->vxc(t, ispin);
        return r;
    }
};
*/


struct binary_munge{

	double longrangevalue;
	double threshold;

	binary_munge(const double thresh = 1.0e-8, const double lrv = 0.0) { // same default value as for munge_thresh
		longrangevalue = lrv;
		threshold = thresh;
	}

    /// @param[out] U   result
    /// @param[in]  f   function to be munged
    /// @param[in]  r   refdensity
    void operator()(const Key<3>& key, Tensor<double>& U, const Tensor<double>& f,
            const Tensor<double>& refdens) const {
        ITERATOR(
            U, double r = refdens(IND);
            double ff = f(IND);
            U(IND) = (r > threshold) ? ff : longrangevalue;
        );
    }

    template <typename Archive>
    void serialize(Archive& ar) {}

};


/// simple structure to take the pointwise logarithm of a function, shifted by +14
struct logme{
    typedef double resultT;
    struct logme1 {
        double operator()(const double& val) {return log(std::max(1.e-14,val))+14.0;}
    };
    Tensor<double> operator()(const Key<3>& key, const Tensor<double>& val) const {
        Tensor<double> result=copy(val);
        logme1 op;
        return result.unaryop(op);
    }

    template <typename Archive>
    void serialize(Archive& ar) {}
};


class OEP : public Nemo {
	typedef std::shared_ptr<real_convolution_3d> poperatorT;

private:

	double dens_thresh = 1.0e-8;   // default 1.0e-8
	double munge_thresh = 1.0e-8;  // default 1.0e-8
	std::vector<bool> oep_model = {false, false, false};

	void set_model_oaep() {oep_model[0] = true;}
	void unset_model_oaep() {oep_model[0] = false;}
	bool is_oaep() const {return oep_model[0];}

	void set_model_ocep() {oep_model[1] = true;}
	void unset_model_ocep() {oep_model[1] = false;}
	bool is_ocep() const {return oep_model[1];}

	void set_model_dcep() {oep_model[2] = true;}
	void unset_model_cdep() {oep_model[2] = false;}
	bool is_dcep() const {return oep_model[2];}

	/// returns true if all members of a vector of booleans are true, otherwise false
	bool IsAlltrue(std::vector<bool> vec) {
		for (int i = 0; i < vec.size(); i++) {
			if (!vec[i]) return false;
		}
		return true;
	}

public:

	OEP(World& world, const std::shared_ptr<SCF> calc) : Nemo(world, calc) {}

	void read_oep_param(std::istream& in) {
        position_stream(in, "oep");
        std::string str;

        while (in >> str) {
            if (str == "end") {
                break;
            }
            else if (str == "oaep") {
                set_model_oaep();
            }
            else if (str == "ocep") {
            	set_model_ocep();
            }
            else if (str == "dcep") {
            	set_model_dcep();
            }
            else if (str == "density_threshold") {
            	in >> dens_thresh;
            }
            else if (str == "munge_threshold") {
            	in >> munge_thresh;
            }
            else {
                print("oep: unrecognized input keyword:", str);
                MADNESS_EXCEPTION("input error",0);
            }
        }
	}


    void solve_oep(const vecfuncT& HF_nemo, const tensorT& HF_eigvals) {

        /// Iterative energy calculation for approximate OEP with EXACT EXCHANGE functional
    	/// for other functionals, slater potential must be modified
    	/// HF orbitals and eigenvalues are used as the guess here
    	/// note that KS_nemo is a reference and changes oep->get_calc()->amo orbitals
    	/// same for orbital energies (eigenvalues) KS_eigvals which is oep->get_calc()->aeps

    	double energy = 0.0;
    	bool update_converged = false; // checks micro iterations and if true, Voep gets updated
    	bool potential_converged = false; // checks macro iterations for convergence of orbital energies
    	unsigned int update_counter = 0; // counts amount of updates in OCEP cycle

    	// exit if no OEP model was set
    	if (!is_oaep() and !is_ocep() and !is_dcep()) {
    		printf("No approximate OEP model selected!\n");
    		printf("please choose oaep/ocep/dcep\n");
    		return;
    	}

    	// for only OAEP calculation, potential (macro) is already converged by construction
    	if (is_oaep() and !is_ocep() and !is_dcep()) potential_converged = true;

		// compute Slater potential Vs and average IHF from HF orbitals and eigenvalues
    	const real_function_3d Vs = compute_slater_potential(HF_nemo, homo_ind(HF_eigvals));
    	const real_function_3d IHF = compute_average_I(HF_nemo, HF_eigvals);
    	save(IHF, "IHF");

    	// set KS_nemo as reference to MOs
    	vecfuncT& KS_nemo = calc->amo;
    	tensorT& KS_eigvals = calc->aeps; // 1d tensor of same length as KS_nemo

    	// all necessary operators applied on nemos (Knemo is used later):
    	vecfuncT Jnemo, Unemo, Vnemo, Knemo;
    	real_function_3d Voep = Vs;

    	// define the solver
    	typedef allocator<double, 3> allocT;
    	typedef XNonlinearSolver<vecfunc<double, 3>, double, allocT> solverT;
    	allocT alloc(world, KS_nemo.size());
    	solverT solver(allocT(world, KS_nemo.size()));


    	// iterations until self-consistency
    	for (int iter = 0; iter < calc->param.maxiter; ++iter) {
//    		save_function(KS_nemo, "KS_nemo_it_"+stringify(iter));
    		for (int i = 0; i < KS_nemo.size(); i++) {
//    		save(KS_nemo[i], "KS_it_"+stringify(iter)+"_nemo_"+stringify(i));
    		}
    		vecfuncT R2KS_nemo = R_square*KS_nemo;
    		truncate(world, R2KS_nemo);

//    		real_function_3d rho = compute_density(KS_nemo);
//    		save(rho, "density_it_"+stringify(iter));

    		if (is_ocep() or is_dcep()) { // only update if OCEP and/or DCEP is enabled

        		// is computed (= updated) only if calculation was converged for smooth convergence
        		if (update_converged) {
        			update_counter++;
        			if (!is_oaep() or update_counter > 1)
        				printf("\n\n     *** updating OCEP potential ***\n\n");

            		// compute OCEP potential from current nemos and eigenvalues
        			real_function_3d corr = compute_OCEP_correction(HF_nemo, HF_eigvals, KS_nemo, KS_eigvals);
        			Voep = Vs + corr;
        			save(corr, "OCEP_correction_update_"+stringify(update_counter));
        			save(Voep, "OCEP_potential_update_"+stringify(update_counter));

        			if (is_oaep() and update_counter == 1) {
        				printf("\n\n     *** V_OAEP converged ***\n");
        				printf("\n  saving V_OCEP with converged OAEP orbitals and eigenvalues\n");
        		    	save(compute_density(KS_nemo), "density_update_1");
        		    	save(compute_average_I(KS_nemo, KS_eigvals), "IKS_update_1");
        		    	printf("     done\n\n");

        	    		printf("\nfinal shifted OAEP orbital energies:\n");
        	    		for (long i = KS_eigvals.size() - 1; i >= 0; i--) {
        	    			printf(" e%2.2lu = %12.8f\n", i, KS_eigvals(i) + homo_diff(HF_eigvals, KS_eigvals));
        	    		}
        	    		printf("HF/KS HOMO energy difference is %12.8f\n", homo_diff(HF_eigvals, KS_eigvals));

        		    	printf("\nFINAL OAEP ENERGY Evir:");
        		    	double Evir = compute_energy(R*KS_nemo, R*Jnemo, Vs, Knemo, true);

        		    	printf("\nFINAL OAEP ENERGY Econv:");
        		    	compute_exchange_potential(KS_nemo, Knemo);
        		    	double Econv = compute_energy(R*KS_nemo, R*Jnemo, Vs, R*Knemo, false);

        		    	printf("OAEP  Evir = %15.8f  Eh", Evir);
        		    	printf("\nOAEP Econv = %15.8f  Eh", Econv);
        		    	printf("\nOAEP DEvir = %15.8f mEh\n", (Evir - Econv)*1000);

        		    	printf("\n\n\n     *** updating OCEP potential ***\n\n");
        			}

        		}

    		}

    		// compute parts of the Fock matrix J, Unuc and Voep
    		compute_nemo_potentials(KS_nemo, Jnemo, Unemo, Voep, Vnemo);

    		// compute Fock matrix F = J + Voep + Vnuc and kinetic energy
    		vecfuncT Fnemo = Jnemo + Vnemo + Unemo;
    		truncate(world, Fnemo);
    		tensorT F = matrix_inner(world, R2KS_nemo, Fnemo, false); // matrix_inner gives 2d tensor
    		Kinetic<double,3> T(world);
    		F += T(R2KS_nemo, KS_nemo); // 2d tensor = Fock-matrix  // R_square in bra, no R in ket

    		// report the off-diagonal Fock matrix elements because canonical orbitals are used
            tensorT F_offdiag = copy(F);
            for (int i = 0; i < F.dim(0); ++i) F_offdiag(i, i) = 0.0;
            double max_F_offidag = F_offdiag.absmax();
            if (world.rank() == 0) print("F max off-diagonal ", max_F_offidag);

    		// compute new (current) energy
            double old_energy = energy;
            printf("energy contributions of iteration %3u:", iter);
    		energy = compute_energy(R*KS_nemo, R*Jnemo, Voep, Knemo, true); // Knemo is not used here
    		// compute_exchange_potential(KS_nemo, Knemo);
    		// energy = compute_energy(R*KS_nemo, R*Jnemo, Voep, R*Knemo, false);
    		// there should be no difference between these two methods, because energy is only needed
    		// for checking convergence threshold; but: Evir should be much faster because K is expensive

    		tensorT old_eigvals = copy(KS_eigvals);

            // diagonalize the Fock matrix to get the eigenvalues and eigenvectors (canonical)
    		// FC = epsilonSC and X^dSX with transform matrix X, see Szabo/Ostlund (3.159) and (3.165)
            tensorT X; // must be formed from R*nemos but can then be used for nemos also
            tensorT overlap = matrix_inner(world, R*KS_nemo, R*KS_nemo, true);
            X = calc->get_fock_transformation(world, overlap, F, KS_eigvals, calc->aocc,
            		FunctionDefaults<3>::get_thresh());
            KS_nemo = transform(world, KS_nemo, X, trantol(), true);
            rotate_subspace(world, X, solver, 0, KS_nemo.size());

            truncate(world, KS_nemo);
            normalize(KS_nemo);

    		// calculate new orbital energies (current eigenvalues from Fock-matrix)
    		for (int i = 0; i < KS_nemo.size(); ++i) {
    			KS_eigvals(i) = std::min(-0.05, F(i, i)); // orbital energy is set to -0.05 if it was above
    		}

    		/// TODO: Question: is this necessary in our programme or even bad?
    		// if requested: subtract orbital shift from orbital energies
    		if (calc->param.orbitalshift > 0.0) {
    			if (world.rank() == 0) print("shifting orbitals by ",
    					calc->param.orbitalshift, " to lower energies");
    			KS_eigvals -= calc->param.orbitalshift;
    		}

    		if (is_ocep() or is_dcep()) { // evaluation only necessary if OCEP and/or DCEP is enabled

        		potential_converged = false; // to ensure convergence is checked again in any case
        		// evaluate convergence of orbital energies if Voep was updated in this iteration (see above)
        		if (update_converged) {
        			// build vector of convergence information of every orbital energy
        			std::vector<bool> conv(KS_eigvals.size());
        			for (int i = 0; i < KS_eigvals.size(); i++) {
        				if (fabs(KS_eigvals(i) - old_eigvals(i)) < calc->param.dconv) conv[i] = true;
        				else conv[i] = false;
        			}

        			if (IsAlltrue(conv)) potential_converged = true; // converged if all are converged
        			update_converged = false; // reset micro iterations so it is checked again at the end
        		}

    		}

    		// print orbital energies:
    		printf("\norbital energies of iteration %3u:\n", iter);
    		for (long i = KS_eigvals.size() - 1; i >= 0; i--) {
    			printf(" e%2.2lu = %12.8f\n", i, KS_eigvals(i));
    		}
    		printf("HF/KS HOMO energy difference is %12.8f\n", homo_diff(HF_eigvals, KS_eigvals));

    		// construct the BSH operators ops
    		std::vector<poperatorT> G = calc->make_bsh_operators(world, KS_eigvals);

    		// remember Fock matrix * nemos from above; make sure it's in phase with nemo (transform)
    		Fnemo = transform(world, Fnemo, X, trantol(), true);
    		truncate(world, Fnemo);

    		// apply the BSH operators G (here ops) on the wave function
    		scale(world, Fnemo, -2.0);
    		vecfuncT GFnemo = apply(world, G, Fnemo);
    		truncate(world, GFnemo);

    		double n1 = norm2(world, KS_nemo);
    		double n2 = norm2(world, GFnemo);
    		print("\nnorm of nemo and GFnemo, ratio ", n1, n2, n1/n2);

    		// compute the residuals
    		vecfuncT residual = KS_nemo - GFnemo;
    		const double norm = norm2(world, residual) / sqrt(KS_nemo.size());

    		// KAIN (helps to converge)
    		vecfuncT nemo_new;
    		if (norm < 5.0e-1) {
    			nemo_new = (solver.update(KS_nemo, residual)).x;
    		} else {
    			nemo_new = GFnemo;
    		}
    		truncate(world, nemo_new);
    		normalize(nemo_new);

    		// What is step restriction?
    		calc->do_step_restriction(world, KS_nemo, nemo_new, "ab spin case");
    		orthonormalize(nemo_new);
    		KS_nemo = nemo_new;

    		// evaluate convergence via norm error and energy difference
    		if ((norm < calc->param.dconv) and (fabs(energy - old_energy) < calc->param.econv))
    			update_converged = true;

    		if (calc->param.save) calc->save_mos(world);

    		if (world.rank() == 0) {
    			printf("\nfinished iteration %2d at time %8.1fs with energy %12.8f\n", iter, wall_time(), energy);
    			print("current residual norm ", norm, "\n");
    		}

    		if (update_converged and potential_converged) break;

    	}

    	if (update_converged and potential_converged) {
    		if (world.rank() == 0) print("\nIterations converged\n");
    	}
    	else {
    		if (world.rank() == 0) print("\nIterations failed\n");
    		energy = 0.0;
    	}

    	// calculate and print all final numbers
    	printf("\n  computing final IKS and density\n");
    	real_function_3d IKS = compute_average_I(KS_nemo, KS_eigvals);
    	real_function_3d rho = compute_density(KS_nemo);
    	save(rho, "density_final");
    	save(IKS, "IKS_final");
    	printf("     done\n");

    	if (is_oaep() and !is_ocep() and !is_dcep()){
    		printf("\n  computing V_OCEP with converged OAEP orbitals and eigenvalues\n");
        	real_function_3d correction = compute_OCEP_correction(HF_nemo, HF_eigvals, KS_nemo, KS_eigvals);
        	real_function_3d ocep_oaep_pot = Vs + correction;
        	save(correction, "OCEP_correction");
        	save(ocep_oaep_pot, "OCEP_potential_with_OAEP_orbs");
    	}
    	if (is_ocep() and !is_dcep()) {
    		printf("\n  computing final V_OCEP with converged OCEP orbitals and eigenvalues\n");
        	real_function_3d correction_final = compute_OCEP_correction(HF_nemo, HF_eigvals, KS_nemo, KS_eigvals);
        	Voep = Vs + correction_final;
        	save(correction_final, "OCEP_correction_final");
        	save(Voep, "OCEP_potential_final");
    	}
    	if (is_dcep()) {
    		printf("\n  computing final V_DCEP with converged DCEP orbitals and eigenvalues\n");
    		// ...
    	}
    	printf("     done\n\n");

    	// print final orbital energies
    	if (is_oaep() and !is_ocep() and !is_dcep()){
    		printf("\nfinal shifted OAEP orbital energies:\n");
    		for (long i = KS_eigvals.size() - 1; i >= 0; i--) {
    			printf(" e%2.2lu = %12.8f\n", i, KS_eigvals(i) + homo_diff(HF_eigvals, KS_eigvals));
    		}
    		printf("HF/KS HOMO energy difference is %12.8f\n", homo_diff(HF_eigvals, KS_eigvals));
    	}
    	if (is_ocep() and !is_dcep()) {
    		printf("\nfinal shifted OCEP orbital energies:\n");
    		for (long i = KS_eigvals.size() - 1; i >= 0; i--) {
    			printf(" e%2.2lu = %12.8f\n", i, KS_eigvals(i));
    		}
    		printf("HF/KS HOMO energy difference is %12.8f\n", homo_diff(HF_eigvals, KS_eigvals));
    	}
    	if (is_dcep()) {
    		printf("\nfinal shifted DCEP orbital energies:\n");
    		for (long i = KS_eigvals.size() - 1; i >= 0; i--) {
    			printf(" e%2.2lu = %12.8f\n", i, KS_eigvals(i));
    		}
    		printf("HF/KS HOMO energy difference is %12.8f\n", homo_diff(HF_eigvals, KS_eigvals));
    	}

    	if (is_oaep() and !is_ocep() and !is_dcep()) printf("\nFINAL OAEP ENERGY Evir:");
    	if (is_ocep() and !is_dcep()) {
    		if (update_converged and potential_converged)
    			printf("\nV_OCEP converged after %3u updates\n", update_counter);
    		printf("\nFINAL OCEP ENERGY Evir:");
    	}
    	if (is_dcep()) printf("\nFINAL DCEP ENERGY Evir:");
    	double Evir = compute_energy(R*KS_nemo, R*Jnemo, Voep, Knemo, true); // Knemo is not used here

    	if (is_oaep() and !is_ocep() and !is_dcep()) printf("\nFINAL OAEP ENERGY Econv:");
    	if (is_ocep() and !is_dcep()) printf("\nFINAL OCEP ENERGY Econv:");
    	if (is_dcep()) printf("\nFINAL DCEP ENERGY Econv:");
    	compute_exchange_potential(KS_nemo, Knemo);
    	double Econv = compute_energy(R*KS_nemo, R*Jnemo, Voep, R*Knemo, false); // Voep is not used here

    	printf("      Evir = %15.8f  Eh", Evir);
    	printf("\n     Econv = %15.8f  Eh", Econv);
    	printf("\n     DEvir = %15.8f mEh\n", (Evir - Econv)*1000);

    }

    /// get index of HOMO from a given set of orbital energies
    long homo_ind(const tensorT orbens) const {
    	long index;
    	double en_homo = orbens.max(&index);
    	return index;
    }

    /// get difference of HF and KS HOMO energies as HOMO_KS - HOMO_HF
    double homo_diff(const tensorT ev1, const tensorT ev2) const {
    	return ev1(homo_ind(ev1)) - ev2(homo_ind(ev2));
    }

    /// TODO: write a function that prints orbital energies

    /// compute density from orbitals with ragularization (Bischoff, 2014_1, equation (19))
    real_function_3d compute_density(const vecfuncT& nemo) const {
    	real_function_3d density = 2.0*R_square*dot(world, nemo, nemo); // 2 because closed shell
    	return density;
    }

    /// compute Slater potential (Kohut, 2014, equation (15))
    real_function_3d compute_slater_potential(const vecfuncT& nemo, const long homo_ind) const {

        Exchange K(world, this, 0); // no - in K here, so factor -1 must be included at the end
        vecfuncT Knemo = K(nemo);
        real_function_3d numerator = 2.0*R_square*dot(world, nemo, Knemo); // 2 because closed shell
        real_function_3d rho = compute_density(nemo);

        // dividing by rho: the minimum value for rho is dens_thresh
        real_function_3d Vs = -1.0*binary_op(numerator, rho, dens_inv(dens_thresh));
        save(Vs, "Slaterpotential_nolra");

        // long-range asymptotic behavior for Slater potential is \int 1/|r-r'| * |phi_HOMO|^2 dr'
        // in order to compute this lra, use Coulomb potential with only HOMO density (= |phi_HOMO|^2)
        Coulomb J(world, this);
        real_function_3d lra = -1.0*J.compute_potential(R_square*square(nemo[homo_ind]));
//        real_function_3d lra = (-0.5/nemo.size())*J.compute_potential(this);

//        for (int i = 0; i < nemo.size(); i++) {
//        	real_function_3d product = R_square*square(nemo[i]);
//        	real_function_3d potential = -1.0*J.compute_potential(R_square*square(nemo[i]));
//        	save(R*nemo[i], "phi"+stringify(i));
//        	save(product, "phi"+stringify(i)+"phi"+stringify(i));
//        	save(potential, "int_phi"+stringify(i)+"phi"+stringify(i));
//        }

        // use apply function from adiabatic correction (see AC.h, nemo.h and nemo.cc) with own potentials
        Vs = ac.apply(Vs, lra);

        save(lra, "lra_slater");
        save(Vs, "Slaterpotential");
        return Vs;

    }

    /// compute average ionization energy I (eigenvalues as 1d tensor)
    real_function_3d compute_average_I(const vecfuncT& nemo, const tensorT eigvals) const {

    	// transform tensor eigvals to vector epsilon
		std::vector<double> epsilon(eigvals.size());
		for (int i = 0; i < eigvals.size(); i++) epsilon[i] = eigvals(i);

		vecfuncT nemo_square = square(world, nemo); // |nemo|^2
		scale(world, nemo_square, epsilon); // epsilon*|nemo|^2
		real_function_3d numerator = 2.0*R_square*sum(world, nemo_square); // 2 because closed shell
        real_function_3d rho = compute_density(nemo);

        // like Kohut, 2014, equations (21) and (25)
        real_function_3d I = -1.0*binary_op(numerator, rho, dens_inv(dens_thresh));

    	// munge I for long-range asymptotic behavior which is -epsilon_HOMO
       	real_function_3d homo_func = real_factory_3d(world).functor([] (const coord_3d& r) {return 1.0;});
       	homo_func.scale(-1.0*eigvals(homo_ind(eigvals)));
       	I = ac.apply(I, homo_func);

//    	I = binary_op(I, rho, binary_munge(munge_thresh, -1.0*eigvals(homo_ind(eigvals))));

        return I;

    }

    /// compute OCEP potential from OAEP Slater potential Vs
    real_function_3d compute_OCEP_correction(const vecfuncT& nemoHF, const tensorT eigvalsHF,
    		const vecfuncT& nemoKS, const tensorT eigvalsKS) const {

    	// Hartree-Fock and Kohn-Sham average ionization energie
    	real_function_3d IHF = compute_average_I(nemoHF, eigvalsHF);
    	real_function_3d IKS = compute_average_I(nemoKS, eigvalsKS);

    	// calculate correction IHF - IKS like Kohut, 2014, equation (26)
    	// and shift potential so that HOMO_HF = HOMO_KS, so potential += (HOMO_HF - HOMO_KS)
    	real_function_3d correction = IHF - IKS + homo_diff(eigvalsHF, eigvalsKS);
    	return correction;

    }

    /// compute all potentials from given nemos except kinetic energy
    void compute_nemo_potentials(const vecfuncT& nemo, vecfuncT& Jnemo, vecfuncT& Unemo,
    		const real_function_3d V, vecfuncT& Vnemo) const {

    	// compute Coulomb part
    	Coulomb J = Coulomb(world, this);
    	Jnemo = J(nemo);
    	truncate(world, Jnemo);

    	// compute nuclear potential part
    	Nuclear Unuc(world, this->nuclear_correlation);
    	Unemo = Unuc(nemo);

    	// compute approximate OEP exchange potential part
    	Vnemo = V*nemo;

    }

    /// compute exchange potential (needed for Econv)
    void compute_exchange_potential(const vecfuncT& nemo, vecfuncT& Knemo) const {

    	Exchange K = Exchange(world, this, 0);
    	Knemo = K(nemo);
    	truncate(world, Knemo);

    }

    /// compute energy from given nemos and given OEP model for exchange
    /// for example Slater potential for OAEP
    double compute_energy(const vecfuncT phi, const vecfuncT Jphi, const real_function_3d Vx,
    		const vecfuncT Kphi, const bool vir) const {

    	// compute kinetic energy
    	// it's ok to use phi here, no regularization necessary for this eigenvalue
    	double E_kin = 0.0;
    	for (int axis = 0; axis < 3; ++axis) {
    		real_derivative_3d D = free_space_derivative<double, 3>(world, axis);
    		const vecfuncT dphi = apply(world, D, phi);
    		E_kin += 0.5 * (inner(world, dphi, dphi)).sum();
    		// -1/2 sum <Psi|Nabla^2|Psi> = 1/2 sum <NablaPsi|NablaPsi>   (integration by parts)
    	}
    	E_kin *= 2.0; // because a = b closed shell

    	// check vir: if true -> compute Evir using Levy-Perdew virial relation (Kohut_2014, (43))
    	// if false -> compute exchange energy using the expectation value of the exchange operator
    	double E_X;
    	if (vir) {
        	/// make vector of functions r = (x, y, z)
        	auto monomial_x = [] (const coord_3d& r) {return r[0];};
        	auto monomial_y = [] (const coord_3d& r) {return r[1];};
        	auto monomial_z = [] (const coord_3d& r) {return r[2];};
        	vecfuncT r(3);
        	r[0]=real_factory_3d(world).functor(monomial_x);
        	r[1]=real_factory_3d(world).functor(monomial_y);
        	r[2]=real_factory_3d(world).functor(monomial_z);
        	/// for density note that phi should be R*nemo, so no R_square is needed
        	real_function_3d rho = 2.0*dot(world, phi, phi); // 2 because closed shell
        	real_function_3d rhoterm = 3*rho + dot(world, r, grad(rho));
        	E_X = inner(Vx, rhoterm);
    	}
    	else {
    		// this uses the oep->get_calc()->amo orbitals in the current form
    		E_X = -1.0*inner(world, phi, Kphi).sum();
    	}

    	// compute external potential (nuclear attraction)
    	real_function_3d Vext = calc->potentialmanager->vnuclear();
    	const vecfuncT Vextphi = Vext*phi;

    	// compute remaining energies: nuclear attraction, Coulomb, nuclear repulsion
    	// computed as expectation values (see Szabo, Ostlund (3.81))
    	const double E_ext = 2.0 * inner(world, phi, Vextphi).sum(); // because a = b closed shell
    	const double E_J = inner(world, phi, Jphi).sum();
    	const double E_nuc = calc->molecule.nuclear_repulsion_energy();
    	double energy = E_kin + E_ext + E_J + E_X + E_nuc;

    	if (world.rank() == 0) {
    		printf("\n                       kinetic energy %15.8f\n", E_kin);
    		printf("   electron-nuclear attraction energy %15.8f\n", E_ext);
    		printf("                       Coulomb energy %15.8f\n", E_J);
    		if (vir) printf(" exchange energy (exchange potential) %15.8f\n", E_X);
    		else printf("  exchange energy (exchange operator) %15.8f\n", E_X);
    		printf("     nuclear-nuclear repulsion energy %15.8f\n", E_nuc);
    		printf("                         total energy %15.8f\n\n", energy);
//            printf("    works for exact exchange functional only...\n");
    	}

    	return energy;

    }


    //
    /// compute the kinetic energy potential using the Kohut trick Eq. (30)
    real_function_3d kinetic_energy_potential2(const vecfuncT& nemo) const {

        // the density
        real_function_3d rhonemo=dot(world,nemo,nemo);
        real_function_3d rho=R_square*rhonemo;

        // compute tauL: Eq. (20) of Kohut
        vecfuncT Rnemo=R*nemo;
        Laplacian<double,3> Laplace(world,0.0);
        vecfuncT laplace_Rnemo=Laplace(Rnemo);
        real_function_3d tauL=-0.5*dot(world,laplace_Rnemo,Rnemo);
        save(tauL,"tauL");
        real_function_3d tauL_div_rho=binary_op(tauL,rho,dens_inv(dens_thresh));
        save(tauL_div_rho,"tauL_div_rho");

        // compute tau = | grad(mo) |^2
        //  = dR.dR * |nemo|^2 + 2*dR.grad(nemo) + R^2*|grad(nemo)|^2
        real_function_3d tau=real_factory_3d(world).compressed();
        vecfuncT dR=this->nuclear_correlation->U1vec();

        NuclearCorrelationFactor::U1_dot_U1_functor u1_dot_u1(nuclear_correlation.get());
        const real_function_3d U1dot=real_factory_3d(world).functor(u1_dot_u1).truncate_on_project();

        for (const real_function_3d& n : nemo) {
            vecfuncT gradnemo=grad(n);
            tau+=U1dot*n*n;
            tau-=2.0*n*dot(world,dR,gradnemo);  // grad(R) = -U1
            tau+=dot(world,gradnemo,gradnemo);
        }
        tau=0.5*tau*R_square;
        real_function_3d tau_div_rho=binary_op(tau,rho,dens_inv(dens_thresh));
        save(tau_div_rho,"tau_div_rho");

        // compute the laplacian of the density with the log trick
        real_function_3d logdensa=unary_op(rho,logme());
        save(logdensa,"logdensa");
        vecfuncT gradzeta=grad(logdensa);
        madness::save_function(gradzeta,"gradzeta");
        real_function_3d d2rho_div_rho=div(gradzeta) + dot(world,gradzeta,gradzeta);
        save(d2rho_div_rho,"d2rho_div_rho");
        real_function_3d result=tau_div_rho-0.25*d2rho_div_rho;
        save(result,"kinetic2");
        return result;
    }

    //
    real_function_3d kinetic_energy_potential(const vecfuncT& nemo) const {

        const Nuclear U_op(world,this->nuclear_correlation);
        const Nuclear V_op(world,this->get_calc().get());

        const vecfuncT Vnemo=V_op(nemo);  // eprec is important here!
        const vecfuncT Unemo=U_op(nemo);

        // nabla^2 nemo
        Laplacian<double,3> Laplace(world,0.0);
        vecfuncT laplace_nemo=Laplace(nemo);

        vecfuncT tmp=Unemo-this->nuclear_correlation->U2()*nemo;
        laplace_nemo-=2.0*tmp;
        vecfuncT D2Rnemo=R*laplace_nemo;

//        // double check result: recompute the density from its laplacian
//        vecfuncT nemo_rec=apply(world,*poisson,D2Rnemo);
//        scale(world,nemo_rec,-1./(4.*constants::pi));
        vecfuncT Rnemo=mul(world,R,nemo);
//        vecfuncT diff=sub(world,Rnemo,nemo_rec);
//        double dnorm=norm2(world,diff);
//        print("dnorm of laplacian phi ",dnorm);

        // compute \sum_i \phi_i \Delta \phi_i
        real_function_3d phiD2phi=dot(world,Rnemo,D2Rnemo);
        save(phiD2phi,"phiD2phi");

        // compute \sum_i \phi_i \epsilon_i \phi_i
        vecfuncT R2nemo=R_square*nemo;
        const real_function_3d rho=2.0*dot(world,nemo,R2nemo);
        const real_function_3d rhonemo=2.0*dot(world,nemo,nemo);

        // compute T1
        double T1 = 0.0;
        for (int axis = 0; axis < 3; axis++) {
            real_derivative_3d D = free_space_derivative<double, 3>(world, axis);
            const vecfuncT dnemo = apply(world, D, Rnemo);
            T1 += 2.0* 0.5 * (inner(world, dnemo, dnemo)).sum();
        }
        printf("T1 %16.8f \n",T1);

        std::vector<double> eps(nemo.size());
        for (std::size_t i=0; i<eps.size(); ++i) eps[i]=calc->aeps(i);
        scale(world,R2nemo,eps);
        real_function_3d phiepsilonphi=dot(world,R2nemo,nemo);

        // divide by the density
        real_function_3d numerator=2.0*(-0.5*phiD2phi);   // fac 2 for sum over spin orbitals
//        real_function_3d numerator=2.0*(-0.5*phiD2phi-phiepsilonphi);   // fac 2 for sum over spin orbitals
        real_function_3d kinetic1=binary_op(numerator,rho,dens_inv(dens_thresh));
        save(kinetic1,"kinetic1");
        real_function_3d nu_bar=kinetic1 + 2.0*(-0.5)*binary_op(phiepsilonphi,rho,dens_inv(dens_thresh));

        // reintroduce the nuclear potential *after* smoothing
        real_function_3d uvnuc=calc->potentialmanager->vnuclear()-nuclear_correlation->U2();
        nu_bar=nu_bar-uvnuc;

        // compute T2
        vecfuncT dipole(3), drho(3);
        dipole[0]=real_factory_3d(world).functor(MomentFunctor(1,0,0));
        dipole[1]=real_factory_3d(world).functor(MomentFunctor(0,1,0));
        dipole[2]=real_factory_3d(world).functor(MomentFunctor(0,0,1));
        drho[0]=make_ddensity(rhonemo,0);
        drho[1]=make_ddensity(rhonemo,1);
        drho[2]=make_ddensity(rhonemo,2);
        real_function_3d one=real_factory_3d(world).functor(MomentFunctor(0,0,0));
        real_function_3d arg=3.0*rho+dot(world,dipole,drho);
        double T2=0.5*inner(nu_bar,arg);
        printf("T2 %16.8f \n",T2);

        Coulomb J(world,this);
        XCOperator xc(world,this);

        real_function_3d vne=calc->potentialmanager->vnuclear();
        real_function_3d vj=J.potential();
        real_function_3d vxc=xc.make_xc_potential();

        real_function_3d euler=nu_bar + vne + vj + vxc;
        double eulernorm=euler.norm2();
        print("eulernorm ",eulernorm);
        save(euler,"euler");
        save(vne,"vne");
        save(vj,"vj");
        save(vxc,"vxc");

        return nu_bar;
    }

    //
    /// compute the laplacian of the density using the log trick
    real_function_3d make_laplacian_density_oep(const real_function_3d& rhonemo) const {


        // U1^2 operator
        NuclearCorrelationFactor::U1_dot_U1_functor u1_dot_u1(nuclear_correlation.get());
        const real_function_3d U1dot=real_factory_3d(world).functor(u1_dot_u1).truncate_on_project();
        real_function_3d result=(2.0*U1dot).truncate();

        // U2 operator
        const real_function_3d V=calc->potentialmanager->vnuclear();
        const real_function_3d U2=nuclear_correlation->U2();

        real_function_3d term2=2.0*(U2-V).truncate();
        result-=term2;


        // compute the laplacian of the density with the log trick
        real_function_3d logdensa=unary_op(rhonemo,logme());
        vecfuncT gradzeta=grad(logdensa);

        // zeta contribution
        real_function_3d d2rho_div_rho=div(gradzeta) + dot(world,gradzeta,gradzeta);
        result+=d2rho_div_rho;

        real_function_3d term3=4.0*dot(world,nuclear_correlation->U1vec(),gradzeta);
        result-=term3;

        result=(R_square*rhonemo*result).truncate();
        save(result,"d2rho");

        // double check result: recompute the density from its laplacian
        real_function_3d rho_rec=-1./(4.*constants::pi)*(*poisson)(result);
        save(rho_rec,"rho_reconstructed");

        real_function_3d rho=rhonemo*R_square;
        save(rho,"rho");
        real_function_3d laplace_rho=div(grad(rho));
        save(laplace_rho,"d2rho_direct");


        return result;
    }

};



int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    if (world.rank() == 0) {
    	print("\n  OEP -- optimized effective potentials for DFT  \n");
    	printf("starting at time %.1f\n", wall_time());
    }
    startup(world, argc, argv);
    std::cout.precision(6);

    const std::string input = "input";
    std::shared_ptr<SCF> calc(new SCF(world, input.c_str())); /// see constructor in SCF.h

    if (world.rank() == 0) {
        calc->molecule.print();
        print("\n");
        calc->param.print(world);
    }

    std::shared_ptr<OEP> oep(new OEP(world, calc));

    vecfuncT HF_nemos;
    tensorT HF_orbens;

    /// TODO: find a way to save eigenvalues and implement restart options
//    const std::string saved_nemos = "HF_nemos";
//    const std::string saved_orbens = "HF_orbens";
//    std::ifstream f1(saved_nemos.c_str());
//    std::ifstream f2(saved_orbens.c_str());
//    if (f1.good() and f2.good()) { /// if file HF_nemos and HF_orbens exist
//    	load_function(world, HF_nemos, saved_nemos);
//    	// load_tensor(... HF_orbens, saved_orbens ...);
//    }
//    else {
//    	const double energy = oep->value();
//    	HF_nemos = copy(world, oep->get_calc()->amo);
//    	HF_orbens = copy(oep->get_calc()->aeps);
//    	save_function(HF_nemos, saved_nemos);
//    	// save_tensor(... HF_orbens, saved_orbens ...);
//    }

    const double energy = oep->value();

    if (world.rank() == 0) {
        printf("final energy   %12.8f\n", energy);
        printf("finished at time %.1f\n", wall_time());
    }

    // save converged HF MOs and orbital energies
    HF_nemos = copy(world, oep->get_calc()->amo);
    HF_orbens = copy(oep->get_calc()->aeps);


//    oep->kinetic_energy_potential(oep->get_calc()->amo);
//    oep->kinetic_energy_potential2(oep->get_calc()->amo);
//    real_function_3d rhonemo = dot(world,oep->get_calc()->amo,oep->get_calc()->amo);
//    oep->make_laplacian_density_oep(rhonemo);


    // OEP model final energy
    printf("\n   +++ starting approximate OEP iterative calculation +++\n\n");

    // read additional OEP parameters from same input file used for SCF calculation (see above)
    std::ifstream in(input.c_str());
    oep->read_oep_param(in);

    oep->solve_oep(HF_nemos, HF_orbens);


    finalize();
    return 0;
}
