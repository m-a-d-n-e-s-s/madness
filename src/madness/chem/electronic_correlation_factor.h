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
#include <iomanip>

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
            real_function_6d result=local+semilocal;
            result-=apply_U_mixed_commutator(phi_i, phi_j,op_mod,thresh,U1nuc);
            result.truncate(FunctionDefaults<6>::get_thresh()*0.3).reduce_rank();

            // do some diagnostics
            if (ao.size()>0) {
                auto ao_ortho=orthonormalize_canonical(ao);
                auto diag = diagnose_Ue(local,semilocal,phi_i,phi_j,ao);
                std::cout << std::scientific << std::setprecision(8);
                print("diagnosis for Ue term:");
                print("local part error:      ", diag.error_local);
                print("semi-local part error: ", diag.error_semilocal);
                print("mixed-comm part error: ", " -- missing");
                print("time in diagnostics:   ", diag.time_ref);
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

        struct Diagnostics {
            Tensor<double> ref_local; // matrix elements for local part
            Tensor<double> ref_semilocal; // matrix elements for semi-local part containing r12 . nabla12
            Tensor<double> result_local; // matrix elements for local part
            Tensor<double> result_semilocal; // matrix elements for semi-local part containing r12 . nabla12
            double error_local=0.0;
            double error_semilocal=0.0;
            double time_ref; // computation time for this
        };

        /// matrix elements <ab|(T-E)^{-1} Ue|ij>: 6D reference vs 3D Schwinger quadrature
        struct DiagnosticsGUe {
            Tensor<double> ref_local;         ///< <ab| G U_local    |ij> from 6D projection
            Tensor<double> ref_semilocal;     ///< <ab| G U_semilocal|ij> from 6D projection
            Tensor<double> result_local;      ///< <ab| G U_local    |ij> from 3D quadrature
            Tensor<double> result_semilocal;  ///< <ab| G U_semilocal|ij> from 3D quadrature
            double error_local = 0.0;
            double error_semilocal = 0.0;
            double time = 0.0;
        };

        /// compute the error in Uphi by comparing to reference values: project onto a (small) aobasis
        Diagnostics diagnose_Ue(const real_function_6d& Uphi_local,
                             const real_function_6d& Uphi_semilocal,
                             const real_function_3d& phi_i, const real_function_3d& phi_j,
                             const std::vector<real_function_3d>& aobasis) const {
            Diagnostics diag;
            double wall0=wall_time();
            timer t(world);
            diag.ref_local=Tensor<double>(aobasis.size(), aobasis.size());
            diag.ref_semilocal=Tensor<double>(aobasis.size(), aobasis.size());
            diag.result_local=Tensor<double>(aobasis.size(), aobasis.size());
            diag.result_semilocal=Tensor<double>(aobasis.size(), aobasis.size());

            // matrix elements for U
            for (int a = 0; a < aobasis.size(); ++a) {
                for (int b = 0; b < aobasis.size(); ++b) {
                    real_function_6d ab=CompositeFactory<double,6,3>(world).particle1(aobasis[a]).particle2(aobasis[b]).is_on_demand();
                    diag.ref_local(a, b) = inner(ab,Uphi_local);
                    diag.ref_semilocal(a, b) = inner(ab,Uphi_semilocal);
                }
            }
            t.tag("compute results");


            // reference for the local part: <ab|U_local|ij>

            // correlation factor is f12 = 1/(2 gamma) (1 - exp(-gamma r12))
            // local part as projected by apply_U_local (see fg_ functor):
            //   U_local = (1 - exp(-gamma r12))/r12 + (gamma/2) exp(-gamma r12)
            // MADNESS operator conventions (see operator.h):
            //   OT_FG12(gamma) kernel = (1 - exp(-gamma r))/(2 gamma * r)   (the 1/(2 gamma) is included)
            //   OT_SLATER(gamma) kernel = exp(-gamma r)
            //   OT_BSH(gamma) kernel = exp(-gamma r)/(4 pi r)               (the 1/(4 pi) is included)
            // therefore:
            //   <ab | (1 - exp(-gamma r12))/r12 | ij> = 2 gamma * inner(b*j, OT_FG12(a*i))
            //   <ab | exp(-gamma r12) | ij>           =          inner(b*j, OT_SLATER(a*i))
            //   <ab | exp(-gamma r12)/r12 | f g>      = 4 pi   * inner(b*g, OT_BSH(a*f))
            auto fg = SeparatedConvolution<double, 3>(
                world, OperatorInfo(_gamma, dcut, FunctionDefaults<3>::get_thresh(), OT_FG12));
            auto slater = SeparatedConvolution<double, 3>(
                world, OperatorInfo(_gamma, dcut, FunctionDefaults<3>::get_thresh(), OT_SLATER));
            auto bsh = SeparatedConvolution<double, 3>(
                world, OperatorInfo(_gamma, dcut, FunctionDefaults<3>::get_thresh(), OT_BSH));

            // matrix elements for local part:
            //   <ab | U_local | ij> = 2 gamma * inner(b*j, fg(a*i)) + (gamma/2) * inner(b*j, slater(a*i))
            for (int a = 0; a < aobasis.size(); ++a) {
                real_function_3d ket = 2.0 * _gamma * fg(aobasis[a] * phi_i)
                                     + 0.5 * _gamma * slater(aobasis[a] * phi_i);
                for (int b = 0; b < aobasis.size(); ++b) {
                    diag.result_local(a, b) = inner(phi_j * aobasis[b], ket);
                }
            }
            t.tag("compute reference for local Ue");

            // matrix elements for semi-local part:
            //   <ab | U_sl | ij> = -1/2 <ab | bsh | r12.nabla_12 |ij>>
            // expand r12.nabla_12 = (r1-r2).(D1-D2) acting on i(r1) j(r2):
            //   r12.nabla_12 |ij> = +A - B - C + D
            // with the four separable contributions
            //   A = (r.D phi_i)(r1) * phi_j(r2)                        =: ~i * j
            //   B = sum_alpha (r_alpha phi_i)(r1) * (D_alpha phi_j)(r2)
            //   C = sum_alpha (D_alpha phi_i)(r1) * (r_alpha phi_j)(r2)
            //   D = phi_i(r1) * (r.D phi_j)(r2)                        =: i * ~j
            // for any separable f(r1) g(r2): <ab|bsh|f g> = inner(b*g, bsh(a*f))
            // (and bsh is symmetric, so swapping a<->b labels in the inner is fine)
            std::vector<real_function_3d> r(3);
            for (int axis = 0; axis < 3; ++axis) r[axis] = real_factory_3d(world).functor([&axis](const coord_3d& xyz) {
                return xyz[axis];
            });

            auto dphi_i = grad(phi_i);
            auto dphi_j = grad(phi_j);
            auto rphi_i = phi_i * r;
            auto rphi_j = phi_j * r;
            real_function_3d itilde=dot(world,r,grad(phi_i));
            real_function_3d jtilde=dot(world,r,grad(phi_j));
            // particle 1 carries (a, i, ...), particle 2 carries (b, j, ...).
            // bsh acts on the particle-1 half; inner contracts with the particle-2 half.
            // For i==j this is symmetric in i<->j so swapping the assignment is invisible;
            // for i!=j swapping gives the transposed matrix elements -- keep the assignment
            // that puts ~i / Di / r*i / i on the bsh side.
            for (int a = 0; a < aobasis.size(); ++a) {
                auto term1=bsh(aobasis[a]*itilde);           // bsh(a*~i)       -> contracts to A
                auto term2=bsh(aobasis[a]*rphi_i);           // bsh(a*r*i)      -> contracts to B
                auto term3=bsh(aobasis[a]*dphi_i);           // bsh(a*Di)       -> contracts to C
                auto term4=bsh(aobasis[a]*phi_i);            // bsh(a*i)        -> contracts to D
                for (int b = 0; b < aobasis.size(); ++b) {
                    double tmp=0.0;
                    tmp += inner(aobasis[b]*phi_j,  term1);  // +A
                    tmp -= inner(aobasis[b]*dphi_j, term2);  // -B (3 axes summed)
                    tmp -= inner(aobasis[b]*rphi_j, term3);  // -C (3 axes summed)
                    tmp += inner(aobasis[b]*jtilde, term4);  // +D
                    // -1/2 prefactor on r12.nabla_12, times 4*pi to undo OT_BSH's 1/(4 pi)
                    diag.result_semilocal(a, b) = -2.0 * constants::pi * tmp;
                }
            }
            t.tag("compute reference for semi-local Ue");

            diag.error_local=(diag.ref_local-diag.result_local).normf();
            diag.error_semilocal=(diag.ref_semilocal-diag.result_semilocal).normf();
            diag.time_ref=wall_time()-wall0;
            return diag;
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
        DiagnosticsGUe diagnose_GUe(const real_function_6d& GUphi_local,
                                    const real_function_6d& GUphi_semilocal,
                                    const real_function_3d& phi_i,
                                    const real_function_3d& phi_j,
                                    const std::vector<real_function_3d>& aobasis,
                                    const double energy) const {
            MADNESS_CHECK_THROW(energy < 0.0, "diagnose_GUe: energy must be negative");
            const double thresh = FunctionDefaults<3>::get_thresh();
            const int nbasis = aobasis.size();

            DiagnosticsGUe diag;
            diag.ref_local        = Tensor<double>(nbasis, nbasis);
            diag.ref_semilocal    = Tensor<double>(nbasis, nbasis);
            diag.result_local     = Tensor<double>(nbasis, nbasis);
            diag.result_semilocal = Tensor<double>(nbasis, nbasis);
            const double wall0 = wall_time();

            // Reference: project G Ue |ij> onto the AO basis pair |ab>
            for (int a = 0; a < nbasis; ++a) {
                for (int b = 0; b < nbasis; ++b) {
                    real_function_6d ab = CompositeFactory<double,6,3>(world)
                        .particle1(copy(aobasis[a])).particle2(copy(aobasis[b])).is_on_demand();
                    diag.ref_local(a, b)     = inner(ab, GUphi_local);
                    diag.ref_semilocal(a, b) = inner(ab, GUphi_semilocal);
                }
            }

            // mu for the 6D Green's function G = (T - E)^{-1}
            const double mu = sqrt(-2.0 * energy);

            // hi: diagonal cell width (same convention as SeparatedConvolution::make_coeff_for_operator)
            const double hi = FunctionDefaults<3>::get_cell_width().normf();

            // BSH 3D Gaussian fit: exp(-mu r)/(4pi r) = sum_n c_n exp(-alpha_n r^2)
            auto fit   = GFit<double,3>::BSHFit(mu, lo, hi, thresh);
            auto c3d   = fit.coeffs();      // 3D expansion coefficients
            auto alpha = fit.exponents();   // Gaussian exponents
            const int nfit = c3d.dim(0);

            // Precompute r-vectors, orbital derivatives, and Euler-weighted orbitals
            std::vector<real_function_3d> rvec(3);
            for (int ax = 0; ax < 3; ++ax)
                rvec[ax] = real_factory_3d(world).functor(
                    [ax](const coord_3d& xyz){ return xyz[ax]; });

            auto dphi_i = grad(phi_i);
            auto dphi_j = grad(phi_j);
            auto rphi_i = phi_i * rvec;
            auto rphi_j = phi_j * rvec;
            real_function_3d itilde = dot(world, rvec, grad(phi_i)); // r.nabla phi_i
            real_function_3d jtilde = dot(world, rvec, grad(phi_j)); // r.nabla phi_j

            // Ue operators — same normalization conventions as in diagnose():
            //   OT_FG12(gamma):   (1-exp(-gamma r))/(2 gamma r)
            //   OT_SLATER(gamma): exp(-gamma r)
            //   OT_BSH(gamma):    exp(-gamma r)/(4 pi r)
            auto fg_op     = SeparatedConvolution<double,3>(world, OperatorInfo(_gamma, dcut, thresh, OT_FG12));
            auto slater_op = SeparatedConvolution<double,3>(world, OperatorInfo(_gamma, dcut, thresh, OT_SLATER));
            auto bsh_op    = SeparatedConvolution<double,3>(world, OperatorInfo(_gamma, dcut, thresh, OT_BSH));

            // Accumulate over Gaussian quadrature nodes
            for (int n = 0; n < nfit; ++n) {
                const double an  = alpha[n];
                const double w6d = c3d[n] * std::pow(an / constants::pi, 1.5);

                // Convolve each AO with exp(-alpha_n r^2)
                auto gauss = SeparatedConvolution<double,3>(world, OperatorInfo(an, lo, thresh, OT_GAUSS));
                std::vector<real_function_3d> conv(nbasis);
                for (int a = 0; a < nbasis; ++a)
                    conv[a] = gauss(aobasis[a]);

                // Local part: <ã_a b̃_b | U_local | ij>
                //   = inner(b̃_b * j,  2*gamma*FG12(ã_a*i) + gamma/2*Slater(ã_a*i))
                for (int a = 0; a < nbasis; ++a) {
                    real_function_3d ket = 2.0 * _gamma * fg_op(conv[a] * phi_i)
                                         + 0.5 * _gamma * slater_op(conv[a] * phi_i);
                    for (int b = 0; b < nbasis; ++b)
                        diag.result_local(a, b) += w6d * inner(phi_j * conv[b], ket);
                }

                // Semi-local part: <ã_a b̃_b | U_sl | ij>
                //   U_sl = -1/2 * BSH_kernel * r12.nabla12
                //   r12.nabla12 |ij> = +A - B - C + D  (see diagnose() for notation)
                //   prefactor -1/2 * 4pi (undo OT_BSH's 1/(4pi)) = -2pi
                for (int a = 0; a < nbasis; ++a) {
                    auto term1 = bsh_op(conv[a] * itilde);      // bsh(ã * ~i)
                    auto term2 = bsh_op(conv[a] * rphi_i);      // bsh(ã * r*i)  [vector]
                    auto term3 = bsh_op(conv[a] * dphi_i);      // bsh(ã * Di)   [vector]
                    auto term4 = bsh_op(conv[a] * phi_i);       // bsh(ã * i)
                    for (int b = 0; b < nbasis; ++b) {
                        double tmp = 0.0;
                        tmp += inner(conv[b] * phi_j,  term1);  // +A
                        tmp -= inner(conv[b] * dphi_j, term2);  // -B  (3 components summed)
                        tmp -= inner(conv[b] * rphi_j, term3);  // -C  (3 components summed)
                        tmp += inner(conv[b] * jtilde, term4);  // +D
                        diag.result_semilocal(a, b) += w6d * (-2.0 * constants::pi * tmp);
                    }
                }
            }

            diag.error_local     = (diag.ref_local     - diag.result_local).normf();
            diag.error_semilocal = (diag.ref_semilocal - diag.result_semilocal).normf();
            diag.time = wall_time() - wall0;
            return diag;
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
