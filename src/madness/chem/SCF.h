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


  $Id$
 */

/// \file moldft.cc
/// \brief Molecular HF and DFT code
/// \defgroup moldft The molecular density functional and Hartree-Fock code


#ifndef MADNESS_CHEM_SCF_H__INCLUDED
#define MADNESS_CHEM_SCF_H__INCLUDED

#include <memory>

#include<madness/chem/molecular_functors.h>
#include <madness/mra/mra.h>

#include<madness/chem/CalculationParameters.h>
#include"madness/mra/commandlineparser.h"
#include<madness/chem/molecule.h>
#include<madness/chem/molecularbasis.h>
#include<madness/chem/corepotential.h>
#include<madness/chem/xcfunctional.h>
#include<madness/chem/potentialmanager.h>
#include<madness/chem/gth_pseudopotential.h>
#include<madness/chem/SCFOperators.h>
#include <madness/tensor/solvers.h>
#include <madness/tensor/distributed_matrix.h>
#include<madness/chem/pcm.h>
#include<madness/chem/QCPropertyInterface.h>

#include <madness/tensor/tensor_json.hpp>
#include <memory>

namespace madness {

typedef std::shared_ptr<WorldDCPmapInterface<Key<3> > > pmapT;
typedef Vector<double, 3> coordT;
typedef std::shared_ptr<FunctionFunctorInterface<double, 3> > functorT;
typedef Function<double, 3> functionT;
typedef std::vector<functionT> vecfuncT;
typedef std::pair<vecfuncT, vecfuncT> pairvecfuncT;
typedef std::vector<pairvecfuncT> subspaceT;
typedef Tensor<double> tensorT;
typedef DistributedMatrix<double> distmatT;
typedef FunctionFactory<double, 3> factoryT;
typedef SeparatedConvolution<double, 3> operatorT;
typedef std::shared_ptr<operatorT> poperatorT;
typedef Function<std::complex<double>, 3> complex_functionT;
typedef std::vector<complex_functionT> cvecfuncT;
typedef Convolution1D<double_complex> complex_operatorT;


template<typename T, int NDIM>
struct lbcost {
    double leaf_value;
    double parent_value;

    lbcost(double leaf_value = 1.0, double parent_value = 0.0) : leaf_value(leaf_value), parent_value(parent_value) {}

    double operator()(const Key<NDIM>& key, const FunctionNode<T, NDIM>& node) const {
        if (key.level() < 1) {
            return 100.0 * (leaf_value + parent_value);
        } else if (node.is_leaf()) {
            return leaf_value;
        } else {
            return parent_value;
        }
    }
};


inline double mask1(double x) {
    /* Iterated first beta function to switch smoothly
           from 0->1 in [0,1].  n iterations produce 2*n-1
           zero derivatives at the end points. Order of polyn
           is 3^n.

           Currently use one iteration so that first deriv.
           is zero at interior boundary and is exactly representable
           by low order multiwavelet without refinement */

    x = (x * x * (3. - 2. * x));
    return x;
}

static double mask3(const coordT& ruser) {
    coordT rsim;
    user_to_sim(ruser, rsim);
    double x = rsim[0], y = rsim[1], z = rsim[2];
    double lo = 0.0625, hi = 1.0 - lo, result = 1.0;
    double rlo = 1.0 / lo;

    if (x < lo)
        result *= mask1(x * rlo);
    else if (x > hi)
        result *= mask1((1.0 - x) * rlo);
    if (y < lo)
        result *= mask1(y * rlo);
    else if (y > hi)
        result *= mask1((1.0 - y) * rlo);
    if (z < lo)
        result *= mask1(z * rlo);
    else if (z > hi)
        result *= mask1((1.0 - z) * rlo);

    return result;
}

/// A MADNESS functor to compute either x, y, or z
class DipoleFunctor : public FunctionFunctorInterface<double, 3> {
private:
    const int axis;
public:
    DipoleFunctor(int axis) : axis(axis) {}

    double operator()(const coordT& x) const {
        return x[axis];
    }
};


/// A MADNESS functor to compute the cartesian moment x^i * y^j * z^k (i, j, k integer and >= 0)
class MomentFunctor : public FunctionFunctorInterface<double, 3> {
private:
    const int i, j, k;
public:
    MomentFunctor(int i, int j, int k) : i(i), j(j), k(k) {}

    MomentFunctor(const std::vector<int>& x) : i(x[0]), j(x[1]), k(x[2]) {}

    double operator()(const coordT& r) const {
        double xi = 1.0, yj = 1.0, zk = 1.0;
        for (int p = 0; p < i; ++p) xi *= r[0];
        for (int p = 0; p < j; ++p) yj *= r[1];
        for (int p = 0; p < k; ++p) zk *= r[2];
        return xi * yj * zk;
    }
};

    class scf_data {

        std::map<std::string, std::vector<double>> e_data;
        int iter;
    public:

        scf_data();

        void to_json(json &j) const;

        void print_data();

        void add_data(std::map<std::string, double> values);
    };


class SCF {
public:
    std::shared_ptr<PotentialManager> potentialmanager;
    std::shared_ptr<GTHPseudopotential<double> > gthpseudopotential;
    Molecule molecule;
    CalculationParameters param;
    XCfunctional xc;
    PCM pcm;
    AtomicBasisSet aobasis;
    functionT mask;

    scf_data e_data;

    /// alpha and beta molecular orbitals
    vecfuncT amo, bmo;

    /// sets of orbitals grouped by their orbital energies (for localization?)
    /// only orbitals within the same set will be mixed to localize
    std::vector<int> aset, bset;

    /// MRA projection of the minimal basis set
    vecfuncT ao;

    std::vector<int> at_to_bf, at_nbf;

    /// occupation numbers for alpha and beta orbitals
    tensorT aocc, bocc;

    /// orbital energies for alpha and beta orbitals
    tensorT aeps, beps;
    poperatorT coulop;
    std::vector<std::shared_ptr<real_derivative_3d> > gradop;
    double vtol;
    double current_energy;
    double converged_for_thresh=1.e10;    ///< mos are converged for this threshold
    //double esol;//etot;
    //double vacuo_energy;

    /// collective constructor for SCF uses contents of file \c filename and broadcasts to all nodes
//	SCF(World & world, const char *filename);
    /// collective constructor for SCF uses contents of stream \c input and broadcasts to all nodes
//	SCF(World & world, std::shared_ptr<std::istream> input);
//	SCF(World& world, const std::string& inputfile);
    SCF(World& world, const commandlineparser& parser);

    void copy_data(World& world, const SCF& other);

    static void help() {
        print_header2("help page for MOLDFT ");
        print("The moldft code computes Hartree-Fock and DFT energies and gradients, It is the fastest code in MADNESS");
        print("and considered the reference implementation. No nuclear correlation factor can be used");
        print("SCF orbitals are the basis for post-SCF calculations like");
        print("excitation energies (cis), correlation energies (cc2), local potentials (oep), etc\n\n");
        print("You can print all available calculation parameters by running\n");
        print("moldft --print_parameters\n");
        print("You can perform a simple calculation by running\n");
        print("moldft --geometry=h2o.xyz\n");
        print("provided you have an xyz file in your directory.");

    }

    static void print_parameters() {
        CalculationParameters param;
        print("default parameters for the moldft program are");
        param.print("dft", "end");
        print("\n\nthe molecular geometry must be specified in a separate block:");
        Molecule::print_parameters();
    }

    void set_print_timings(const bool value);

    template<std::size_t NDIM>
    void set_protocol(World& world, double thresh) {
        int k;
        // Allow for imprecise conversion of threshold
        if (thresh >= 0.9e-2)
            k = 4;
        else if (thresh >= 0.9e-4)
            k = 6;
        else if (thresh >= 0.9e-6)
            k = 8;
        else if (thresh >= 0.9e-8)
            k = 10;
        else
            k = 12;

        // k defaults to make sense with thresh, override by providing k in input file
        if (param.k() == -1) {
            FunctionDefaults<NDIM>::set_k(k);
            //        	param.k=k;
        } else {
            FunctionDefaults<NDIM>::set_k(param.k());
        }
        // don't forget to adapt the molecular smoothing parameter!! NO ... it is independent
        //        molecule.set_eprec(std::min(thresh,molecule.get_eprec()));
        FunctionDefaults<NDIM>::set_thresh(thresh);
        FunctionDefaults<NDIM>::set_refine(true);
        FunctionDefaults<NDIM>::set_initial_level(2);
        //        FunctionDefaults<NDIM>::set_truncate_mode(1);
        FunctionDefaults<NDIM>::set_autorefine(false);
        FunctionDefaults<NDIM>::set_apply_randomize(false);
        FunctionDefaults<NDIM>::set_project_randomize(false);
        FunctionDefaults<NDIM>::set_cubic_cell(-param.L(), param.L());
        GaussianConvolution1DCache<double>::map.clear();
        double safety = 0.1;
        vtol = FunctionDefaults<NDIM>::get_thresh() * safety;
        coulop = poperatorT(CoulombOperatorPtr(world, param.lo(), thresh));
        gradop = gradient_operator<double, 3>(world);

        // Update coefficients if using a different derivative
        if (param.deriv() == "bspline") {
            for (int i = 0; i < 3; ++i) (*gradop[i]).set_bspline1();
        } else if (param.deriv() == "ble") {
            for (int i = 0; i < 3; ++i) (*gradop[i]).set_ble1();
        }

        mask = functionT(factoryT(world).f(mask3).initial_level(4).norefine());
        if (world.rank() == 0 and param.print_level() > 1) {
            print("\nSolving NDIM=", NDIM, " with thresh", thresh, "    k",
                  FunctionDefaults<NDIM>::get_k(), "   conv", std::max(thresh, param.dconv()), "\n");
        }
    }

    /// getter for the molecular orbitals, alpha spin
    const vecfuncT& get_amo() const { return amo; }

    /// getter for the molecular orbitals, beta spin
    const vecfuncT& get_bmo() const { return bmo; }

    /// getter for the occupation numbers, alpha spin
    const tensorT& get_aocc() const { return aocc; }

    /// getter for the occupation numbers, alpha spin
    const tensorT& get_bocc() const { return bocc; }

    bool is_spin_restricted() const { return param.get<bool>("spin_restricted"); }

    void save_mos(World& world);

    void load_mos(World& world);

    bool restart_aos(World& world);

    void do_plots(World& world);

    void project(World& world);

    void make_nuclear_potential(World& world);

    vecfuncT project_ao_basis(World& world, const AtomicBasisSet& aobasis);

    static vecfuncT project_ao_basis_only(World& world, const AtomicBasisSet& aobasis,
                                          const Molecule& molecule);

    void reset_aobasis(const std::string& aobasisname) {
        aobasis = AtomicBasisSet(); // reset
        aobasis.read_file(aobasisname);
    }

    /// group orbitals into sets of similar orbital energies for localization

    /// @param[in]	eps	orbital energies
    /// @param[in]	occ	occupation numbers
    /// @param[in]	nmo number of MOs for the given spin
    /// @return		vector of length nmo with the set index for each MO
    std::vector<int> group_orbital_sets(World& world, const tensorT& eps,
                                        const tensorT& occ, const int nmo) const;

    static void analyze_vectors(World& world, const vecfuncT& mo,
            const vecfuncT& ao, double vtol,
            const Molecule& molecule, const int print_level,
            const AtomicBasisSet& aobasis, const tensorT& occ = tensorT(),
            const tensorT& energy = tensorT(), const std::vector<int>& set = std::vector<int>());

    distmatT kinetic_energy_matrix(World& world, const vecfuncT& v) const;

    void get_initial_orbitals(World& world);

    void initial_guess(World& world);

    void initial_guess_from_nwchem(World& world);

    void initial_load_bal(World& world);

    functionT make_density(World& world, const tensorT& occ, const vecfuncT& v) const;

    functionT make_density(World& world, const tensorT& occ, const cvecfuncT& v);

    std::vector<poperatorT> make_bsh_operators(World& world, const tensorT& evals) const;

    // Used only for initial guess that is always spin-restricted LDA
    static functionT make_lda_potential(World& world, const functionT& arho);


    //    functionT make_dft_potential(World & world, const vecfuncT& vf, int ispin, int what)
    //    {
    //        return multiop_values<double, xc_potential, 3>(xc_potential(xc, ispin, what), vf);
    //    }

    double make_dft_energy(World& world, const vecfuncT& vf, int ispin) {
        functionT vlda = multiop_values<double, xc_functional, 3>(xc_functional(xc), vf);
        return vlda.trace();
    }

    vecfuncT apply_potential(World& world, const tensorT& occ,
                             const vecfuncT& amo,
                             const functionT& vlocal, double& exc, double& enl, int ispin);

    tensorT derivatives(World& world, const functionT& rho) const;

    /// compute the total dipole moment of the molecule

    /// @param[in]  rho the total (alpha + beta) density
    /// @return     the x,y,z components of the el. + nucl. dipole moment
    tensorT dipole(World& world, const functionT& rho) const;

    void vector_stats(const std::vector<double>& v, double& rms,
                      double& maxabsval) const;

    vecfuncT compute_residual(World& world, tensorT& occ, tensorT& fock,
                              const vecfuncT& psi, vecfuncT& Vpsi, double& err);

    tensorT make_fock_matrix(World& world, const vecfuncT& psi,
                             const vecfuncT& Vpsi, const tensorT& occ,
                             double& ekinetic) const;

    /// make the Coulomb potential given the total density
    functionT make_coulomb_potential(const functionT& rho) const {
        return apply(*coulop, rho);
    }

    /// Compute the two-electron integrals over the provided set of orbitals

    /// Returned is a *replicated* tensor of \f$(ij|kl)\f$ with \f$i>=j\f$
    /// and \f$k>=l\f$.  The symmetry \f$(ij|kl)=(kl|ij)\f$ is enforced.
    Tensor<double> twoint(World& world, const vecfuncT& psi) const;

    /// compute the unitary transformation that diagonalizes the fock matrix

    /// @param[in]	world	the world
    /// @param[in]	overlap	the overlap matrix of the orbitals
    /// @param[in,out]	fock	the fock matrix; diagonal upon exit
    /// @param[out]	evals	the orbital energies
    /// @param[in]	occ	the occupation numbers
    /// @param[in]	thresh_degenerate	threshold for orbitals being degenerate
    /// @return		the unitary matrix U: U^T F U = evals
    tensorT get_fock_transformation(World& world, const tensorT& overlap,
                                    tensorT& fock, tensorT& evals, const tensorT& occ,
                                    const double thresh_degenerate) const;


    /// diagonalize the fock matrix, taking care of degenerate states

    /// Vpsi is passed in to make sure orbitals and Vpsi are in phase
    /// @param[in]	world	the world
    /// @param[in,out]	fock	the fock matrix (diagonal upon exit)
    /// @param[in,out]	psi		the orbitals
    /// @param[in,out]	Vpsi	the orbital times the potential
    /// @param[out]	evals	the orbital energies
    /// @param[in]	occ		occupation numbers
    /// @param[in]	thresh	threshold for rotation and truncation
    /// @return		the unitary matrix U: U^T F U = evals
    tensorT diag_fock_matrix(World& world, tensorT& fock,
                             vecfuncT& psi, vecfuncT& Vpsi, tensorT& evals,
                             const tensorT& occ, const double thresh) const;


    void loadbal(World& world, functionT& arho, functionT& brho, functionT& arho_old,
                 functionT& brho_old, subspaceT& subspace);


    void rotate_subspace(World& world, const tensorT& U, subspaceT& subspace,
                         int lo, int nfunc, double trantol) const;

    void rotate_subspace(World& world, const distmatT& U, subspaceT& subspace,
                         int lo, int nfunc, double trantol) const;

    void update_subspace(World& world,
                         vecfuncT& Vpsia, vecfuncT& Vpsib,
                         tensorT& focka, tensorT& fockb,
                         subspaceT& subspace, tensorT& Q,
                         double& bsh_residual, double& update_residual);

    /// perform step restriction following the KAIN solver

    /// undo the rotation from the KAIN solver if the rotation exceeds the
    /// maxrotn parameter
    /// @param[in]		world	the world
    /// @param[in]		mo		vector of orbitals from previous iteration
    /// @param[in,out]	mo_new	vector of orbitals from the KAIN solver
    /// @param[in]		spin	"alpha" or "beta" for user information
    /// @return			max residual
    double do_step_restriction(World& world, const vecfuncT& mo,
                               vecfuncT& mo_new, std::string spin) const;

    /// orthonormalize the vectors

    /// @param[in]		world	the world
    /// @param[in,out]	amo_new	the vectors to be orthonormalized
    void orthonormalize(World& world, vecfuncT& amo_new) const;

    void orthonormalize(World& world, vecfuncT& amo_new, int nocc) const;

    complex_functionT APPLY(const complex_operatorT *q1d, const complex_functionT& psi) {
        complex_functionT r = psi;  // Shallow copy violates constness !!!!!!!!!!!!!!!!!
        coordT lo, hi;
        lo[2] = -10;
        hi[2] = +10;

        r.reconstruct();
        r.broaden();
        r.broaden();
        r.broaden();
        r.broaden();
        r = apply_1d_realspace_push(*q1d, r, 2);
        r.sum_down();
        r = apply_1d_realspace_push(*q1d, r, 1);
        r.sum_down();
        r = apply_1d_realspace_push(*q1d, r, 0);
        r.sum_down();

        return r;
    }

    // For given protocol, solve the DFT/HF/response equations
    void solve(World& world);

    void output_calc_info_schema() const;

//    void output_scf_info_schema(const std::map<std::string, double> &vals,
//                                const tensorT &dipole_T) const;

};

// Computes molecular energy as a function of the geometry
// This is cludgy ... need better factorization of functionality
// between calculation, main program and this ... or just merge it all.
class MolecularEnergy : public OptimizationTargetInterface, public QCPropertyInterface {
    World& world;
    SCF& calc;
    mutable double coords_sum;     // sum of square of coords at last solved geometry

public:
    MolecularEnergy(World& world, SCF& calc)
            : world(world), calc(calc), coords_sum(-1.0) {}

    std::string name() const { return "Molecularenerg"; }

    bool selftest() { return true; }

    bool provides_gradient() const { return true; }

    double value(const Tensor<double>& x) {
        double xsq = x.sumsq();
        if (xsq == coords_sum) {
            return calc.current_energy;
        }
        calc.molecule.set_all_coords(x.reshape(calc.molecule.natom(), 3));
        coords_sum = xsq;

        // read converged wave function from disk if there is one
        if (calc.param.no_compute()) {
            calc.load_mos(world);
            calc.make_nuclear_potential(world);
            calc.ao = calc.project_ao_basis(world, calc.aobasis);
            return calc.current_energy;
        }

        // initialize the PCM solver for this geometry
        if (calc.param.pcm_data() != "none") {
            calc.pcm = PCM(world, calc.molecule, calc.param.pcm_data(), true);
        }

        calc.get_initial_orbitals(world);

        // AOs are needed for final analysis, and for localization
        calc.reset_aobasis("sto-3g");
        calc.ao.clear(); world.gop.fence(); 
        calc.ao = calc.project_ao_basis(world, calc.aobasis);

        // The below is missing convergence test logic, etc.

        // Make the nuclear potential, initial orbitals, etc.
        for (unsigned int proto = 0; proto < calc.param.protocol().size(); proto++) {

            int nvalpha = calc.param.nmo_alpha() - calc.param.nalpha();
            int nvbeta = calc.param.nmo_beta() - calc.param.nbeta();
            int nvalpha_start, nv_old;

            //repeat with gradually decreasing nvirt, only for first protocol
            if (proto == 0 && nvalpha > 0) {
                nvalpha_start = nvalpha * calc.param.nv_factor();
            } else {
                nvalpha_start = nvalpha;
            }

            nv_old = nvalpha_start;

            for (int nv = nvalpha_start; nv >= nvalpha; nv -= nvalpha) {

                if (nv > 0 && world.rank() == 0) std::cout << "Running with " << nv << " virtual states" << std::endl;

                calc.param.set_user_defined_value("nmo_alpha", calc.param.nalpha() + nv);
                // check whether this is sensible for spin restricted case
                if (calc.param.nbeta() && !calc.param.spin_restricted()) {
                    if (nvbeta == nvalpha) {
                        calc.param.set_user_defined_value("nmo_beta", calc.param.nbeta() + nv);
                    } else {
                        calc.param.set_user_defined_value("nmo_beta", calc.param.nbeta() + nv + nvbeta - nvalpha);
                    }
                }

                calc.set_protocol<3>(world, calc.param.protocol()[proto]);
                calc.make_nuclear_potential(world);

                if (nv != nv_old) {
                    calc.amo.resize(calc.param.nmo_alpha());
                    calc.bmo.resize(calc.param.nmo_beta());

                    calc.aocc = tensorT(calc.param.nmo_alpha());
                    for (int i = 0; i < calc.param.nalpha(); ++i)
                        calc.aocc[i] = 1.0;

                    calc.bocc = tensorT(calc.param.nmo_beta());
                    for (int i = 0; i < calc.param.nbeta(); ++i)
                        calc.bocc[i] = 1.0;

                    // might need to resize aset, bset, but for the moment this doesn't seem to be necessary

                }

                // project orbitals into higher k
                if (proto > 0) calc.project(world);

                // If the basis for the inital guess was not sto-3g
                // switch to sto-3g since this is needed for analysis
                // of the MOs and orbital localization
                // Only do this if not starting from NWChem.
                // analysis will be done on NWChem orbitals.

                if (calc.param.aobasis() != "sto-3g") { // was also  && calc.param.nwfile() == "none"
                    calc.reset_aobasis("sto-3g");
                }
                calc.ao.clear(); world.gop.fence(); 
                calc.ao = calc.project_ao_basis(world, calc.aobasis);
                calc.solve(world);

                if (calc.param.save())
                    calc.save_mos(world);

                nv_old = nv;
                // exit loop over decreasing nvirt if nvirt=0
                if (nv == 0) break;

            }

        }
        return calc.current_energy;
    }

    madness::Tensor<double> gradient(const Tensor<double>& x) {
        value(x); // Ensures DFT equations are solved at this geometry

        functionT rho = calc.make_density(world, calc.aocc, calc.amo);
        functionT brho = rho;
        if (!calc.param.spin_restricted())
            brho = calc.make_density(world, calc.bocc, calc.bmo);
        rho.gaxpy(1.0, brho, 1.0);

        return calc.derivatives(world, rho);
    }


    void energy_and_gradient(const Molecule& molecule, double& energy, Tensor<double>& gradient) {
        value(molecule.get_all_coords().flat()); // Ensures DFT equations are solved at this geometry

        functionT rho = calc.make_density(world, calc.aocc, calc.amo);
        functionT brho = rho;
        if (!calc.param.spin_restricted())
            brho = calc.make_density(world, calc.bocc, calc.bmo);
        rho.gaxpy(1.0, brho, 1.0);

        energy = calc.current_energy;
        gradient = calc.derivatives(world, rho);
    }

    void output_calc_info_schema() {
        nlohmann::json j = {};
        vec_pair_ints int_vals;
        vec_pair_T<double> double_vals;
        vec_pair_tensor_T<double> double_tensor_vals;

        CalculationParameters param = calc.param;

        nlohmann::json calc_precision={ };
        calc_precision["eprec"]=calc.molecule.parameters.eprec();
        calc_precision["dconv"]=calc.param.dconv();
        calc_precision["econv"]=calc.param.econv();
        calc_precision["thresh"]=FunctionDefaults<3>::get_thresh();
        calc_precision["k"]=FunctionDefaults<3>::get_k();

        auto mol_json=this->calc.molecule.to_json();

        int_vals.push_back({"calcinfo_nmo", param.nmo_alpha() + param.nmo_beta()});
        int_vals.push_back({"calcinfo_nalpha", param.nalpha()});
        int_vals.push_back({"calcinfo_nbeta", param.nbeta()});
        int_vals.push_back({"calcinfo_natom", calc.molecule.natom()});


        to_json(j, int_vals);
        double_vals.push_back({"return_energy", value(calc.molecule.get_all_coords().flat())});
        to_json(j, double_vals);
        double_tensor_vals.push_back({"scf_eigenvalues_a", calc.aeps});
        if (param.nbeta() != 0 && !param.spin_restricted()) {
            double_tensor_vals.push_back({"scf_eigenvalues_b", calc.beps});
        }

        to_json(j, double_tensor_vals);
        param.to_json(j);
        calc.e_data.to_json(j);

        j["precision"]=calc_precision;
        j["molecule"]=mol_json;

        output_schema(param.prefix()+".calc_info", j);
    }


};
}

#endif /* SCF_H_ */

