/*
 * electronic_correlation_factor.h
 *
 *  Created on: Jul 9, 2015
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_ELECTRONIC_CORRELATION_FACTOR_H_
#define SRC_APPS_CHEM_ELECTRONIC_CORRELATION_FACTOR_H_


#include <madness/mra/mra.h>
#include <madness/mra/lbdeux.h>
#include<madness/chem/molecule.h>
#include<madness/world/timing_utilities.h>
#include <madness/chem/operator_diagnostics.h>
#include <iomanip>

#include "CCStructures.h"

namespace madness {
    /// a class holding the electronic correlation factor for R12 theory
    class CorrelationFactor {
        World& world;
        double _gamma; ///< the correlation factor exp(-gamma r12)
        double dcut; ///< the cutoff for the 1/r potential
        double lo; ///< smallest length scale to be resolved
        int truncate_mode=1;

    public:
        /// ctor, use negative gamma for linear correlation factor r12
        CorrelationFactor(World& world) : world(world), _gamma(-1.0), dcut(1.e-10),
                                          lo(1.e-10) {}

        /// ctor, use negative gamma for linear correlation factor r12
        CorrelationFactor(World& world, const double& gamma, const double dcut,
                          const Molecule& molecule) : world(world), _gamma(gamma), dcut(dcut) {
            lo = 1.e-6; //lo = molecule.smallest_length_scale();
            //        if (world.rank()==0) {
            //            if (gamma>0.0) print("constructed correlation factor with gamma=",gamma);
            //            else if (gamma==0.0) print("constructed linear correlation factor");
            //        }
        }

        /// ctor, use negative gamma for linear correlation factor r12
        CorrelationFactor(World& world, const double& gamma, const double dcut,
                          const double lo) : world(world), _gamma(gamma), dcut(dcut), lo(lo) {
            //        if (world.rank()==0) {
            //            if (gamma>0.0) print("constructed correlation factor with gamma=",gamma);
            //            else if (gamma==0.0) print("constructed linear correlation factor");
            //        }
        }

        /// copy ctor
        CorrelationFactor(const CorrelationFactor& other) : world(other.world) {
            _gamma = other._gamma;
            dcut = other.dcut;
            lo = other.lo;
        }

        /// assignment; assume other's world is this world
        CorrelationFactor& operator=(const CorrelationFactor& other) {
            _gamma = other._gamma;
            dcut = other.dcut;
            lo = other.lo;
            return *this;
        }

        /// return the exponent of this correlation factor
        double gamma() const { return _gamma; }

        void set_truncate_mode(int mode) {
            truncate_mode=mode;
        }

        /// return the value of the correlation factor
        double operator()(const coord_6d& r) const {
            const double rr = r12(r);
            if (_gamma > 0.0) return (1.0 - exp(-_gamma * rr)) / (2.0 * _gamma);
            return 0.5 * rr;
        }

        /// compute the local part of the Ue potential
        real_function_6d apply_U_local(const real_function_3d& phi_i, const real_function_3d& phi_j,
            const real_convolution_6d& op_mod, double thresh) const {
            if (world.rank()==0) print("applying local part of Ue term");
            double fac=1.0;
            if (truncate_mode<0) fac = 0.1;
            print("truncate mode",truncate_mode);
            print("factor for thresh",fac);
            fg_ func(_gamma, dcut);
            real_function_6d fg3 = real_factory_6d(world).functor(func).is_on_demand();
            real_function_6d mul = CompositeFactory<double, 6, 3>(world)
                                   .g12(fg3).particle1(copy(phi_i)).particle2(copy(phi_j)).thresh(thresh*fac)
                                    .truncate_mode(truncate_mode);
            mul.fill_cuspy_tree(op_mod).truncate(FunctionDefaults<6>::get_thresh()*0.3);
            mul.print_size("local Ue|ij>");
            return mul;
        }

        /// compute the semi-local part of the Ue potential containing derivatives
        real_function_6d apply_U_semilocal(const real_function_3d& phi_i, const real_function_3d& phi_j,
                                        const real_convolution_6d& op_mod, double thresh, const bool symmetric = false) const {
            if (world.rank()==0) print("applying semi-local part of Ue term");
            double fac=1.0;
            if (truncate_mode<0) fac = 0.1;
            print("truncate mode",truncate_mode);
            print("factor for thresh",fac);

            auto dphi_i=grad(phi_i);
            auto dphi_j=grad(phi_j);

            real_function_6d result = real_factory_6d(world).truncate_mode(truncate_mode).thresh(thresh*fac);
            for (int axis = 0; axis < 3; ++axis) {
                real_function_6d u = U1(axis);

                std::vector<real_function_3d> p1({copy(dphi_i[axis]),copy(phi_i)});
                std::vector<real_function_3d> p2({copy(phi_j),-1.0*copy(dphi_j[axis])});

                real_function_6d tmp = CompositeFactory<double, 6, 3>(world)
                       .g12(u).particle1(p1).particle2(p2).thresh(thresh*fac)
                                    .truncate_mode(truncate_mode);
                tmp.fill_cuspy_tree(op_mod).truncate();
                result = result + tmp;
                world.gop.fence();
            }
            result.truncate(FunctionDefaults<6>::get_thresh()*0.3).reduce_rank();
            result.print_size("semi-local Ue|ij>");
            return result;
        }

        /// compute the double commutator electronic/nuclear correlation factor:  R^{-1}[[T,f],R]

        /// @param[in]  U1nuc the U1 potential of the nuclear correlation factor
        real_function_6d apply_U_mixed_commutator(const real_function_3d& phi_i, const real_function_3d& phi_j,
                                        const real_convolution_6d& op_mod, double thresh,
                                        const std::vector<real_function_3d> U1nuc)const {
            if (world.rank()==0) print("applying mixed double commutator of the Ue term");
            double fac=1.0;
            if (truncate_mode<0) fac = 0.1;
            print("truncate mode",truncate_mode);
            print("factor for thresh",fac);

            auto dphi_i=U1nuc*phi_i;
            auto dphi_j=U1nuc*phi_j;

            real_function_6d result = real_factory_6d(world).truncate_mode(truncate_mode).thresh(thresh*fac);
            for (int axis = 0; axis < 3; ++axis) {
                real_function_6d u = U1(axis);

                std::vector<real_function_3d> p1({copy(dphi_i[axis]),copy(phi_i)});
                std::vector<real_function_3d> p2({copy(phi_j),-1.0*copy(dphi_j[axis])});

                real_function_6d tmp = CompositeFactory<double, 6, 3>(world)
                       .g12(u).particle1(p1).particle2(p2).thresh(thresh*fac)
                                    .truncate_mode(truncate_mode);
                tmp.fill_cuspy_tree(op_mod).truncate();
                result = result + tmp;
                world.gop.fence();
            }
            result.truncate(FunctionDefaults<6>::get_thresh()*0.3).reduce_rank();
            result.print_size("mixed commutator of Ue: R^{-1} [[T,f12],R] |ij>");
            return result;
        }


        /// apply Kutzelnigg's regularized potential to an orbital product
        real_function_6d apply_U(const real_function_3d& phi_i, const real_function_3d& phi_j,
                                 const real_convolution_6d& op_mod, const std::vector<real_function_3d>& U1nuc,
                                 const std::vector<real_function_3d> ao=std::vector<real_function_3d>()) const {
            MADNESS_CHECK_THROW(op_mod.modified(), "ElectronicCorrelationFactor::apply_U, op_mod must be in modified_NS form");
            const double thresh = FunctionDefaults<6>::get_thresh();

            real_function_6d local=apply_U_local(phi_i, phi_j,op_mod,thresh);
            real_function_6d semilocal=apply_U_semilocal(phi_i, phi_j,op_mod,thresh);
            real_function_6d mixed=apply_U_mixed_commutator(phi_i, phi_j,op_mod,thresh,U1nuc);
            real_function_6d result=local+semilocal-mixed;
            result.truncate(FunctionDefaults<6>::get_thresh()*0.3).reduce_rank();

            // do some diagnostics
            if (ao.size()>0) {
                auto ao_ortho=orthonormalize_canonical(ao);
                auto dm = diagnose_Ue(local,semilocal,phi_i,phi_j,ao,mixed,U1nuc);
                std::cout << std::scientific << std::setprecision(8);
                print("diagnosis for Ue term:");
                dm.print_report("Ue");
            }

            result.print_size("Ue|ij>");
            return result;
        }

        /// return the U1 term of the correlation function
        real_function_6d U1(const int axis) const {
            U func(_gamma, axis, dcut);
            const real_function_6d u1 = real_factory_6d(world)
                                        .functor(func).is_on_demand();
            return u1;
        }

        /// return the U1 term of the correlation function
        real_function_6d U2() const {
            if (world.rank() == 0) print("U2 for the electronic correlation factor");
            if (world.rank() == 0) print("is expensive -- do you really need it??");
            MADNESS_EXCEPTION("U2() not implemented, since it might be expensive", 1);
            return real_factory_6d(world);
        }

        /// return the correlation factor as on-demand function
        real_function_6d f() const {
            //        real_function_6d tmp=real_factory_6d(world).functor2(*this).is_on_demand();
            double thresh = FunctionDefaults<3>::get_thresh();
            real_function_6d tmp = TwoElectronFactory(world)
                                   .dcut(dcut).gamma(_gamma).f12().thresh(thresh);
            return tmp;
        }

        /// return f^2 as on-demand function
        real_function_6d f2() const {
            f2_ func(_gamma);
            real_function_6d tmp = real_factory_6d(world).functor(func).is_on_demand();
            return tmp;
        }

        /// return fg+sth as on-demand function
        real_function_6d fg() const {
            fg_ func(_gamma, dcut);
            real_function_6d tmp = real_factory_6d(world).functor(func).is_on_demand();
            return tmp;
        }

        /// return f/r as on-demand function
        real_function_6d f_over_r() const {
            f_over_r_ func(_gamma, dcut);
            real_function_6d tmp = real_factory_6d(world).functor(func).is_on_demand();
            return tmp;
        }

        /// return (\nabla f)^2 as on-demand functions
        real_function_6d nablaf2() const {
            nablaf2_ func(_gamma);
            real_function_6d tmp = real_factory_6d(world).functor(func).is_on_demand();
            return tmp;
        }

        /// compute the error in Uphi by comparing to reference values: project onto a (small) aobasis
        /// @param[in]  Uphi_mixed   U_mixed|ij> (optional; pass default-constructed to skip)
        /// @param[in]  U1nuc        nuclear U1 potentials (required when Uphi_mixed is provided)
        /// @return DiagnosticMatrix with entries "local", "semilocal" (and "mixed" if applicable).
        ///         ref=3D analytic formula, result=6D projection <ab|U|ij>, error=||ref-result||.
        DiagnosticMatrix diagnose_Ue(const real_function_6d& Uphi_local,
                             const real_function_6d& Uphi_semilocal,
                             const real_function_3d& phi_i, const real_function_3d& phi_j,
                             const std::vector<real_function_3d>& aobasis,
                             const real_function_6d& Uphi_mixed = real_function_6d(),
                             const std::vector<real_function_3d>& U1nuc = {}) const {
            const bool do_mixed = Uphi_mixed.is_initialized() && !U1nuc.empty();
            double wall0=wall_time();
            timer t(world);

            DiagnosticMatrix dm(world, aobasis);
            dm.init("local");
            dm.init("semilocal");
            if (do_mixed) dm.init("mixed");

            // Fill result: 6D projections <ab|U|ij>
            dm.entries["local"].result     = dm.project_ab(Uphi_local);
            dm.entries["semilocal"].result = dm.project_ab(Uphi_semilocal);
            if (do_mixed) dm.entries["mixed"].result = dm.project_ab(Uphi_mixed);
            t.tag("compute 6D projections");


            // reference for the local part: <ab|U_local|ij>

            // Fill ref: 3D analytic formulas for each piece.
            //
            // correlation factor f12 = 1/(2 gamma) (1 - exp(-gamma r12))
            // local part (see fg_ functor): U_local = (1-exp(-gamma r12))/r12 + (gamma/2) exp(-gamma r12)
            // MADNESS operator conventions (see operator.h):
            //   OT_FG12(gamma) kernel = (1 - exp(-gamma r))/(2 gamma * r)
            //   OT_SLATER(gamma) kernel = exp(-gamma r)
            //   OT_BSH(gamma) kernel = exp(-gamma r)/(4 pi r)
            // therefore:
            //   <ab | U_local | ij> = 2 gamma * inner(b*j, fg(a*i)) + (gamma/2) * inner(b*j, slater(a*i))
            // auto fg = SeparatedConvolution<double, 3>(
                // world, OperatorInfo(_gamma, dcut, FunctionDefaults<3>::get_thresh(), OT_FG12));
            auto fg=std::shared_ptr<SeparatedConvolution<double,3>>(FGOperatorPtr(world, 1.0, 1.e-6, FunctionDefaults<3>::get_thresh()));
            auto slater = SeparatedConvolution<double, 3>(
                world, OperatorInfo(_gamma, dcut, FunctionDefaults<3>::get_thresh(), OT_SLATER));
            auto bsh = SeparatedConvolution<double, 3>(
                world, OperatorInfo(_gamma, dcut, FunctionDefaults<3>::get_thresh(), OT_BSH));

            // set up local part as CCPairFunction
            std::vector<CCPairFunction<double,6>> U_local;
            U_local.push_back(CCPairFunction<double,6>(fg,phi_i,phi_j));

//            auto& ref_local = dm.entries["local"].ref;
//            for (int a = 0; a < (int)aobasis.size(); ++a) {
//                real_function_3d ket = 2.0 * _gamma * fg(aobasis[a] * phi_i)
//                                     + 0.5 * _gamma * slater(aobasis[a] * phi_i);
//                for (int b = 0; b < (int)aobasis.size(); ++b)
//                    ref_local(a, b) = inner(phi_j * aobasis[b], ket);
//            }
            t.tag("3D ref for local Ue");

            // semi-local part: <ab|U_sl|ij> = -1/2 <ab|bsh|r12.nabla_12|ij>
            //   r12.nabla_12 |ij> = +A - B - C + D  (see comments in old diagnose_Ue)
            std::vector<real_function_3d> r(3);
            for (int axis = 0; axis < 3; ++axis)
                r[axis] = real_factory_3d(world).functor([&axis](const coord_3d& xyz) {
                    return xyz[axis];
                });

            auto dphi_i = grad(phi_i);
            auto dphi_j = grad(phi_j);
            auto rphi_i = phi_i * r;
            auto rphi_j = phi_j * r;
            real_function_3d itilde = dot(world, r, grad(phi_i));
            real_function_3d jtilde = dot(world, r, grad(phi_j));

            auto& ref_sl = dm.entries["semilocal"].ref;
            for (int a = 0; a < (int)aobasis.size(); ++a) {
                auto term1 = bsh(aobasis[a] * itilde);   // bsh(a*~i)  -> A
                auto term2 = bsh(aobasis[a] * rphi_i);   // bsh(a*r*i) -> B
                auto term3 = bsh(aobasis[a] * dphi_i);   // bsh(a*Di)  -> C
                auto term4 = bsh(aobasis[a] * phi_i);    // bsh(a*i)   -> D
                for (int b = 0; b < (int)aobasis.size(); ++b) {
                    double tmp = 0.0;
                    tmp += inner(aobasis[b] * phi_j,  term1);  // +A
                    tmp -= inner(aobasis[b] * dphi_j, term2);  // -B
                    tmp -= inner(aobasis[b] * rphi_j, term3);  // -C
                    tmp += inner(aobasis[b] * jtilde, term4);  // +D
                    ref_sl(a, b) = -2.0 * constants::pi * tmp;
                }
            }
            t.tag("3D ref for semilocal Ue");

            // mixed commutator: same as semilocal with ∇φ → U1nuc*φ
            if (do_mixed) {
                auto Wphii = U1nuc * phi_i;
                auto Wphij = U1nuc * phi_j;
                real_function_3d itilde_nuc = dot(world, r, Wphii);
                real_function_3d jtilde_nuc = dot(world, r, Wphij);
                auto& ref_mx = dm.entries["mixed"].ref;
                for (int a = 0; a < dm.nbasis(); ++a) {
                    auto term1 = bsh(aobasis[a] * itilde_nuc);
                    auto term2 = bsh(aobasis[a] * rphi_i);
                    auto term3 = bsh(aobasis[a] * Wphii);
                    auto term4 = bsh(aobasis[a] * phi_i);
                    for (int b = 0; b < dm.nbasis(); ++b) {
                        double tmp = 0.0;
                        tmp += inner(aobasis[b] * phi_j,      term1);
                        tmp -= inner(aobasis[b] * Wphij,      term2);
                        tmp -= inner(aobasis[b] * rphi_j,     term3);
                        tmp += inner(aobasis[b] * jtilde_nuc, term4);
                        ref_mx(a, b) = -2.0 * constants::pi * tmp;
                    }
                }
                t.tag("3D ref for mixed Ue");
            }

            dm.compute_errors();
            dm.time = wall_time() - wall0;
            return dm;
        }

        /// compute <ab|(T-E)^{-1} Ue|ij> via 3D Gaussian quadrature of the Schwinger representation
        /// and compare against the 6D reference obtained by projecting G U_local/U_semilocal |ij>.
        ///
        /// The 6D Green's function G = (T-E)^{-1} is written via the Schwinger proper-time expansion
        ///   G(r1,r2;r1',r2') = int_0^inf e^{tE} h_t(r1-r1') h_t(r2-r2') dt
        /// where h_t is the 3D heat kernel.  Discretizing via the BSH Gaussian fit with
        /// mu = sqrt(-2E) gives the separable approximation
        ///   G ≈ sum_n w_n^{6d}  g_{alpha_n}(r1-r1') g_{alpha_n}(r2-r2')
        /// with 6D weights  w_n^{6d} = c_n^{bsh} * (alpha_n/pi)^{3/2}.
        ///
        /// The matrix elements then reduce to
        ///   M(a,b) ≈ sum_n w_n^{6d} <ã_n b̃_n | Ue | ij>
        /// where ã_n = g_{alpha_n} * a (3D Gaussian convolution).
        ///
        /// @param[in]  GUphi_local      G U_local    |ij> as a 6D MRA function (reference)
        /// @param[in]  GUphi_semilocal  G U_semilocal|ij> as a 6D MRA function (reference)
        /// @param[in]  phi_i    orbital i (particle 1)
        /// @param[in]  phi_j    orbital j (particle 2)
        /// @param[in]  aobasis  projector basis {a}
        /// @param[in]  energy   two-particle energy E = eps_i + eps_j  (must be < 0)
        /// @param[in]  GUphi_mixed  G U_mixed|ij> (optional; pass default-constructed to skip)
        /// @param[in]  U1nuc        nuclear U1 potentials (required when GUphi_mixed is provided)
        /// @return DiagnosticMatrix with entries "Glocal", "Gsemilocal" (and "Gmixed" if applicable).
        ///         ref=3D Schwinger quadrature, result=6D projection <ab|G Ue|ij>, error=||ref-result||.
        DiagnosticMatrix diagnose_GUe(const real_function_6d& GUphi_local,
                                    const real_function_6d& GUphi_semilocal,
                                    const real_function_3d& phi_i,
                                    const real_function_3d& phi_j,
                                    const std::vector<real_function_3d>& aobasis,
                                    const double energy,
                                    const real_function_6d& GUphi_mixed = real_function_6d(),
                                    const std::vector<real_function_3d>& U1nuc = {}) const {
            MADNESS_CHECK_THROW(energy < 0.0, "diagnose_GUe: energy must be negative");
            const bool do_mixed = GUphi_mixed.is_initialized() && !U1nuc.empty();
            const double thresh = FunctionDefaults<3>::get_thresh();
            const int nbasis = aobasis.size();
            const double wall0 = wall_time();

            DiagnosticMatrix dm(world, aobasis);
            dm.init("Glocal");
            dm.init("Gsemilocal");
            if (do_mixed) dm.init("Gmixed");

            // Fill result: project 6D functions G Ue |ij> onto <ab|
            dm.entries["Glocal"].result     = dm.project_ab(GUphi_local);
            dm.entries["Gsemilocal"].result = dm.project_ab(GUphi_semilocal);
            if (do_mixed) dm.entries["Gmixed"].result = dm.project_ab(GUphi_mixed);

            // Fill ref: 3D Schwinger quadrature  <ab|G Ue|ij> ≈ Σ_n w_n <ã_n b̃_n|Ue|ij>
            const double mu = sqrt(-2.0 * energy);
            const double hi = FunctionDefaults<3>::get_cell_width().normf();
            auto fit   = GFit<double,3>::BSHFit(mu, lo, hi, thresh);
            auto c3d   = fit.coeffs();
            auto alpha = fit.exponents();
            const int nfit = c3d.dim(0);

            std::vector<real_function_3d> rvec(3);
            for (int ax = 0; ax < 3; ++ax)
                rvec[ax] = real_factory_3d(world).functor(
                    [ax](const coord_3d& xyz){ return xyz[ax]; });

            auto dphi_i = grad(phi_i);
            auto dphi_j = grad(phi_j);
            auto rphi_i = phi_i * rvec;
            auto rphi_j = phi_j * rvec;
            real_function_3d itilde = dot(world, rvec, grad(phi_i));
            real_function_3d jtilde = dot(world, rvec, grad(phi_j));

            std::vector<real_function_3d> Wphii, Wphij;
            real_function_3d itilde_nuc, jtilde_nuc;
            if (do_mixed) {
                Wphii = U1nuc * phi_i;
                Wphij = U1nuc * phi_j;
                itilde_nuc = dot(world, rvec, Wphii);
                jtilde_nuc = dot(world, rvec, Wphij);
            }

            auto fg_op     = SeparatedConvolution<double,3>(world, OperatorInfo(_gamma, dcut, thresh, OT_FG12));
            auto slater_op = SeparatedConvolution<double,3>(world, OperatorInfo(_gamma, dcut, thresh, OT_SLATER));
            auto bsh_op    = SeparatedConvolution<double,3>(world, OperatorInfo(_gamma, dcut, thresh, OT_BSH));

            auto& ref_local = dm.entries["Glocal"].ref;
            auto& ref_sl    = dm.entries["Gsemilocal"].ref;

            for (int n = 0; n < nfit; ++n) {
                const double an  = alpha[n];
                const double w6d = c3d[n] * std::pow(an / constants::pi, 1.5);

                auto gauss = SeparatedConvolution<double,3>(world, OperatorInfo(an, lo, thresh, OT_GAUSS));
                std::vector<real_function_3d> conv(nbasis);
                for (int a = 0; a < nbasis; ++a)
                    conv[a] = gauss(aobasis[a]);

                // Local: <ã b̃|U_local|ij> = inner(b̃*j, 2γ*fg(ã*i) + γ/2*slater(ã*i))
                for (int a = 0; a < nbasis; ++a) {
                    real_function_3d ket = 2.0 * _gamma * fg_op(conv[a] * phi_i)
                                         + 0.5 * _gamma * slater_op(conv[a] * phi_i);
                    for (int b = 0; b < nbasis; ++b)
                        ref_local(a, b) += w6d * inner(phi_j * conv[b], ket);
                }

                // Semi-local: <ã b̃|U_sl|ij> = -2π × (+A -B -C +D)
                for (int a = 0; a < nbasis; ++a) {
                    auto term1 = bsh_op(conv[a] * itilde);
                    auto term2 = bsh_op(conv[a] * rphi_i);
                    auto term3 = bsh_op(conv[a] * dphi_i);
                    auto term4 = bsh_op(conv[a] * phi_i);
                    for (int b = 0; b < nbasis; ++b) {
                        double tmp = 0.0;
                        tmp += inner(conv[b] * phi_j,  term1);
                        tmp -= inner(conv[b] * dphi_j, term2);
                        tmp -= inner(conv[b] * rphi_j, term3);
                        tmp += inner(conv[b] * jtilde, term4);
                        ref_sl(a, b) += w6d * (-2.0 * constants::pi * tmp);
                    }
                }

                // Mixed commutator
                if (do_mixed) {
                    auto& ref_mx = dm.entries["Gmixed"].ref;
                    for (int a = 0; a < nbasis; ++a) {
                        auto term1 = bsh_op(conv[a] * itilde_nuc);
                        auto term2 = bsh_op(conv[a] * rphi_i);
                        auto term3 = bsh_op(conv[a] * Wphii);
                        auto term4 = bsh_op(conv[a] * phi_i);
                        for (int b = 0; b < nbasis; ++b) {
                            double tmp = 0.0;
                            tmp += inner(conv[b] * phi_j,      term1);
                            tmp -= inner(conv[b] * Wphij,      term2);
                            tmp -= inner(conv[b] * rphi_j,     term3);
                            tmp += inner(conv[b] * jtilde_nuc, term4);
                            ref_mx(a, b) += w6d * (-2.0 * constants::pi * tmp);
                        }
                    }
                }
            }

            dm.compute_errors();
            dm.time = wall_time() - wall0;
            return dm;
        }

    private:
        /// functor for the local potential (1-f12)/r12 + sth (doubly connected term of the commutator)

        /// TODO: turn this into coeffs directly
        struct fg_ : FunctionFunctorInterface<double, 6> {
            double gamma;
            double dcut;

            fg_(double gamma, double dcut) : gamma(gamma), dcut(dcut) {
                MADNESS_ASSERT(gamma>0.0);
            }

            double operator()(const coord_6d& r) const {
                const double rr = r12(r);
                const double e = exp(-gamma * rr);
                if (rr < 5.e-2) {
                    double value = 1.5 * gamma - 1. * std::pow(gamma, 2) * rr + 0.41666666666666663 * std::pow(gamma, 3)
                        * std::pow(rr, 2) -
                        0.125 * std::pow(gamma, 4) * std::pow(rr, 3) + 0.029166666666666667 * std::pow(gamma, 5) *
                        std::pow(rr, 4)
                        + 0.005555555555555556 * std::pow(gamma, 6) * std::pow(rr, 5) +
                        0.0008928571428571429 * std::pow(gamma, 7) * std::pow(rr, 6);
                    return value;
                }
                else {
                    // return (1.0-e)*u(rr,dcut) + 0.5*gamma*e;
                    return (1.0 - e) / rr + 0.5 * gamma * e;
                }
            }
        };

        /// functor for the local potential (1-f12)/r12
        struct f_over_r_ : FunctionFunctorInterface<double, 6> {
            double gamma;
            double dcut;

            f_over_r_(double gamma, double dcut) : gamma(gamma), dcut(dcut) {
                MADNESS_ASSERT(gamma>0.0);
            }

            double operator()(const coord_6d& r) const {
                const double rr = r12(r);
                const double e = exp(-gamma * rr);
                return (1.0 - e) * u(rr, dcut) / (2.0 * gamma);
            }
        };

        /// functor for the local part of the regularized potential: f12/r12*(r1-r2)(D1-D2)
        struct U : FunctionFunctorInterface<double, 6> {
            double gamma;
            int axis;
            double dcut;

            U(double gamma, int axis, double dcut) : gamma(gamma), axis(axis),
                                                     dcut(dcut) {
                MADNESS_ASSERT(axis>=0 and axis<3);
            }

            double operator()(const coord_6d& r) const {
                const double rr = r12(r);
                const coord_3d vr12{r[0] - r[3], r[1] - r[4], r[2] - r[5]};
                const coord_3d N = unitvec(vr12);
                if (gamma > 0.0) return -0.5 * exp(-gamma * rr) * N[axis];
                MADNESS_EXCEPTION("no gamma in electronic corrfac::U1", 1);
                //          const double rr=r12(r);
                //            const double g12=u(rr,dcut);
                //            double a=0.5;
                //            if (gamma>0.0) a=0.5*exp(-gamma*rr);
                //            return -a*x12(r,axis) * g12;
            }
        };

        /// functor for the local potential (1-f12)^2
        struct f2_ : FunctionFunctorInterface<double, 6> {
            double gamma;
            f2_(double gamma) : gamma(gamma) { MADNESS_ASSERT(gamma>0.0); }

            double operator()(const coord_6d& r) const {
                const double rr = r12(r);
                const double e = exp(-gamma * rr);
                const double f = (1.0 - e) / (2.0 * gamma);
                return f * f;
            }
        };

        /// functor for the local potential (\nabla f)^2
        struct nablaf2_ : FunctionFunctorInterface<double, 6> {
            double gamma;

            nablaf2_(double gamma) : gamma(gamma) {
                MADNESS_ASSERT(gamma>0.0);
                MADNESS_ASSERT(gamma==1.0);
            }

            double operator()(const coord_6d& r) const {
                const double rr = r12(r);
                const double f = exp(-2.0 * gamma * rr) / (4.0 * gamma * gamma);
                return f;
            }
        };

        /// Smoothed 1/r potential (c is the smoothing distance)
        static double u(double r, double c) {
            r = r / c;
            double r2 = r * r, pot;
            if (r > 6.5) {
                pot = 1.0 / r;
            }
            else if (r > 1e-2) {
                pot = erf(r) / r + exp(-r2) * 0.56418958354775630;
            }
            else {
                pot = 1.6925687506432689 - r2 * (0.94031597257959381 - r2 * (0.39493270848342941 - 0.12089776790309064 *
                    r2));
            }
            return pot / c;
        }

        static double r12(const coord_6d& r) {
            const double x12 = r[0] - r[3];
            const double y12 = r[1] - r[4];
            const double z12 = r[2] - r[5];
            const double r12 = sqrt(x12 * x12 + y12 * y12 + z12 * z12);
            return r12;
        }

        static double x12(const coord_6d& r, const int axis) {
            return r[axis] - r[axis + 3];
        }

        static coord_3d smoothed_unitvec(const coord_3d& xyz, double smoothing) {
            //        if (smoothing==0.0) smoothing=molecule.get_eprec();
            // TODO:need to test this
            // reduce the smoothing for the unitvector
            //if (not (this->type()==None or this->type()==Two)) smoothing=sqrt(smoothing);
            const double r = xyz.normf();
            const double cutoff = smoothing;
            if (r > cutoff) {
                return 1.0 / r * xyz;
            }
            else {
                const double xi = r / cutoff;
                const double xi2 = xi * xi;
                const double xi3 = xi * xi * xi;
                //            const double nu21=0.5+1./32.*(45.*xi - 50.*xi3 + 21.*xi*xi*xi*xi*xi);
                const double nu22 = 0.5 + 1. / 64. * (105 * xi - 175 * xi3 + 147 * xi2 * xi3 - 45 * xi3 * xi3 * xi);
                //            const double nu40=0.5 + 1./128.*(225 *xi - 350 *xi3 + 189*xi2*xi3);
                const double kk = 2. * nu22 - 1.0;
                return kk / (r + 1.e-15) * xyz;
            }
        }
    };

    /// a class holding the electronic correlation factor for R12 theory
    /// CorrelationFactor2 = (1-0.5*exp(-gamma*r12), gamma=0.5
    /// (CorrelationFactor + 1)*2.0 = CorrelationFactor2 (currently CorrelationFactor2 is only implemented for gamma=0.5 so use this gamma also on CorrelationFactor
    class CorrelationFactor2 {
        World& world;
        double _gamma; ///< the correlation factor exp(-gamma r12)
        typedef std::shared_ptr<FunctionFunctorInterface<double, 6>> functorT;

    public:
        double dcut; ///< the cutoff for the 1/r potential
        double lo; ///< smallest length scale to be resolved
        double vtol; ///< initial projection threshold


        /// ctor, use negative gamma for linear correlation factor r12
        CorrelationFactor2(World& world) : world(world), _gamma(0.5), dcut(1.e-10),
                                           lo(1.e-10), vtol(FunctionDefaults<3>::get_thresh() * 0.1) {
            MADNESS_ASSERT(_gamma==0.5);
        }

        /// return the exponent of this correlation factor
        double gamma() const { return _gamma; }

        real_function_6d function() const {
            functorT R = functorT(new R_functor(_gamma, 1));
            return real_factory_6d(world).functor(R).is_on_demand();
        }

        real_function_6d square() const {
            functorT R2 = functorT(new R_functor(_gamma, 2));
            return real_factory_6d(world).functor(R2).is_on_demand();
        }

        real_function_6d inverse() const {
            functorT R = functorT(new R_functor(_gamma, -1));
            return real_factory_6d(world).functor(R).is_on_demand();
        }

        /// return the U1 term of the correlation function
        real_function_6d U1(const int axis) const {
            functorT U1f = functorT(new U1_functor(_gamma, axis));
            return real_factory_6d(world).functor(U1f).is_on_demand();
        }

        /// return the U2 term of the correlation function
        real_function_6d U2() const {
            functorT U2f = functorT(new U2_functor(_gamma));
            return real_factory_6d(world).functor(U2f).is_on_demand();
        }

        /// apply Kutzelnigg's regularized potential to an orbital product
        real_function_6d apply_U(const real_function_6d& psi, const double eps) const {
            const double bsh_thresh = 1.e-7;

            real_function_6d result = real_factory_6d(world);

            real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2 * eps), lo, bsh_thresh);
            op_mod.modified() = true;

            for (int axis = 0; axis < 3; ++axis) {
                //if (world.rank()==0) print("working on axis",axis);
                real_derivative_6d D1 = free_space_derivative<double, 6>(world, axis);
                real_derivative_6d D2 = free_space_derivative<double, 6>(world, axis + 3);
                const real_function_6d Drhs1 = D1(psi).truncate();
                const real_function_6d Drhs2 = D2(psi).truncate();

                const real_function_6d u1 = U1(axis);

                real_function_6d tmp1 = CompositeFactory<double, 6, 3>(world)
                                        .g12(u1).ket(copy(Drhs1));
                tmp1.fill_cuspy_tree(op_mod).truncate();

                real_function_6d tmp2 = CompositeFactory<double, 6, 3>(world)
                                        .g12(u1).ket(copy(Drhs2));
                tmp2.fill_cuspy_tree(op_mod).truncate();
                // if (world.rank()==0) print("done with fill_tree");

                result = result + (tmp1 - tmp2).truncate();
                tmp1.clear();
                tmp2.clear();
                world.gop.fence();
                result.truncate().reduce_rank();

                // if (world.rank()==0) printf("done with multiplication with U at ime %.1f\n",wall_time());
                // result.print_size("result");
            }

            real_function_6d u2 = U2();
            real_function_6d r2 = CompositeFactory<double, 6, 3>(world).ket(copy(psi))
                                                                       .g12(u2);
            r2.fill_tree(op_mod);
            result = (result + r2).truncate();
            return result;
        }

    private:
        /// functor for the correlation factor R
        class R_functor : public FunctionFunctorInterface<double, 6> {
            double gamma;
            int exponent;

        public:
            R_functor(double gamma, int e = 1) : gamma(gamma), exponent(e) {
                MADNESS_ASSERT(gamma==0.5);
            }

            // only valid for gamma=1
            double operator()(const coord_6d& r) const {
                const double rr = r12(r);
                double val = (1.0 - 0.5 * exp(-gamma * rr));
                if (exponent == 1) return val;
                else if (exponent == 2) return val * val;
                else if (exponent == -1) return 1.0 / val;
                else {
                    MADNESS_EXCEPTION("fancy exponent in correlationfactor2", 1);
                }
            }
        };

        /// functor for the U2 local potential
        class U2_functor : public FunctionFunctorInterface<double, 6> {
            double gamma;

        public:
            U2_functor(double gamma) : gamma(gamma) {
                MADNESS_ASSERT(gamma==0.5);
            }

            // only valid for gamma=1
            double operator()(const coord_6d& r) const {
                const double rr = r12(r);
                // Taylor expansion for small r
                if (rr < 1.e-4) { // valid for gamma==0.5, otherwise singular
                    return (5. / 4.0 - rr + (35.0 * rr * rr) / 48.0 - (101.0 * rr * rr * rr) / 192.0);
                }
                const double egr = exp(-gamma * rr);
                return -(-8. * egr + 8.0 + rr * egr) / (4.0 * rr * egr - 8 * rr);
            }
        };

        /// functor for the U1 = -\frac{\vec\nabla_1 f_{12}}{f_{12}}  potential

        /// the potential is given by
        /// U1 = -\frac{\vec\nabla_1 f_{12}}{f_{12}}
        ///    =  \frac{e^{-r12/2}{4-2e^{-r12/2}} \vec unitvec
        /// the derivative operators are not included
        class U1_functor : public FunctionFunctorInterface<double, 6> {
            double gamma;
            int axis;

        public:
            U1_functor(double gamma, int axis) : gamma(gamma), axis(axis) {
                MADNESS_ASSERT(gamma==0.5);
                MADNESS_ASSERT(axis<3);
            }

            double operator()(const coord_6d& r) const {
                const double rr = r12(r);
                const coord_3d vr12{r[0] - r[3], r[1] - r[4], r[2] - r[5]};
                const coord_3d N = unitvec(vr12);
                // Taylor expansion for small r
                double val;
                if (rr < 1.e-4) { // valid for gamma==0.5, otherwise singular
                    val = 0.5 - 0.5 * rr + 0.125 * (3. * rr * rr) - (13. * rr * rr * rr) / 48.0;
                }
                else {
                    const double egr = exp(-gamma * rr);
                    val = egr / (4.0 - 2.0 * egr);
                }
                // NOTE the sign
                return -val * N[axis];
            }
        };

        /// helper function
        static double r12(const coord_6d& r) {
            const double x12 = r[0] - r[3];
            const double y12 = r[1] - r[4];
            const double z12 = r[2] - r[5];
            const double r12 = sqrt(x12 * x12 + y12 * y12 + z12 * z12);
            return r12;
        }
    };
}


#endif /* SRC_APPS_CHEM_ELECTRONIC_CORRELATION_FACTOR_H_ */
