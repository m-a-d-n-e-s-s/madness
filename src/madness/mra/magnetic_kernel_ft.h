/// \file magnetic_kernel_ft.h
/// \brief Fourier approach: (T_B - E)^{-1} by a Fourier separation in y (A. Armen, branch magnetic-fourier)
///
/// T_B = T_perp(x, y) - 1/2 d^2/dz^2, so e^{-tau T_B} = U_perp(tau) K_z(tau), and U_perp separates exactly as
///
///   U_perp = e^{iBxy/2} (1/2pi) int dq e^{iq(y-y')} g(x-X_q) N(x-x') g(x'-X_q) e^{-iBx'y'/2},
///
///   X_q = -q/B,  g(u) = exp(-m u^2),  N(d) = sqrt(c/pi) exp(-c d^2),  m = (B/2) tanh(B tau/2),  c = B/(2 sinh(B tau)),
///
/// i.e. the Mehler factorization with frequency B in x for each y wavenumber q, with the mask written as
/// m (x - X_q)^2 = m x^2 + b q x + d q^2, b = tanh(B tau/2), d = tanh(B tau/2)/(2B).  Applied to f:
///
///   1. ft = e^{-iBxy/2} f (on the cylinder rho < rcut, erfc taper)
///   2. F_q(x, z) = int dy e^{-iqy} ft(x, y, z)                 (project_out, per q)
///   3. Psi_q = sum_tau w g_q(x) * [(N_x K_z) * (g_q(x) F_q)]   (x mask, 2D separable convolution, x mask)
///   4. result = e^{iBxy/2} (1/P) sum_q e^{iqy} Psi_q(x, z)     (hartree_product)
///
/// with q_j = 2 pi j/P (y treated as periodic with period P) and a trapezoid rule in ln tau from tmin to
/// tmax_factor/(|B|/2 - E), [0, tmin] added as tmin f.  Steps 2 and 4 work in the axis order (x, z, y).
/// Gauge origin r = 0; any orientation of the molecule; the default rcut assumes the nuclei near the z axis.
/// Port of MagneticGreen (magnetic_green.h on magnetic-fourier) onto magnetic_primitives.h.

#ifndef MADNESS_MRA_MAGNETIC_KERNEL_FT_H__INCLUDED
#define MADNESS_MRA_MAGNETIC_KERNEL_FT_H__INCLUDED

#include <madness/mra/magnetic_kernel.h>

namespace madness {

class MagneticKernelFT : public MagneticKernel {
public:
    using operator2dT = SeparatedConvolution<double, 2>;

    MagneticKernelFT(World& world, const MagneticKernelParameters& p) : MagneticKernel(p), world_(world) {
        MADNESS_CHECK_THROW(p_.ox == 0.0 && p_.oy == 0.0, "MagneticKernelFT: the gauge origin must be r = 0");
        magnetic::set_lowdim_defaults_from_3d();
        const Tensor<double> cell = FunctionDefaults<3>::get_cell();
        lo_ = cell(0, 0);
        hi_ = cell(0, 1);
        const double B = p_.B;
        const double thr = FunctionDefaults<3>::get_thresh();
        // the phase is stored only on the cylinder rho < rcut (+ erfc taper): f and G f decay as
        // e^{-|B| rho^2/4} across B, and a full-cell 3D phase refines everywhere (memory)
        rcut_ = p_.ft_rcut > 0 ? p_.ft_rcut : (B == 0.0 ? 1e10 : 1.5 * std::sqrt(8.0 * std::log(1.0 / thr) / std::abs(B)));
        rtaper_ = p_.ft_rtaper > 0 ? p_.ft_rtaper : std::min(0.5, rcut_ / 4);
        // ft lives on |x|, |y| < rext: y-wavenumbers up to |B| rext/2 + those of f, and the period
        // P = 4 rext keeps the y-images of ft and of G ft apart (P <= 2 x cell width)
        const double rext = rcut_ + 3 * rtaper_;
        const double L = std::min(std::max(std::abs(lo_), std::abs(hi_)), rext);
        period_ = p_.ft_period > 0 ? p_.ft_period : std::min(2.0 * (hi_ - lo_), 4.0 * rext);
        const double qmax = std::abs(B) * L / 2 + p_.ft_qextra;
        const long J = long(std::ceil(qmax * period_ / constants::pi / 2));
        for (long j = -J; j <= J; ++j) q_.push_back(2 * constants::pi * j / period_);

        phase_.refinement = magnetic::PhaseRefinement::box_phase;
        phase_.dphi = p_.ft_dphi;
        phase_.rext = rext;
        phase_.tol = 0.1 * thr;
        phase_.rcut = rcut_;
        phase_.rtaper = rtaper_;

        // step 4 uses e^{iqy} windowed to |y| < rcut + 9 rtaper: every term of the q-sum is then local in y
        const double yw = rcut_ + 9 * rtaper_, rt = rtaper_;
        for (double q : q_) {
            // e^{-iqy} of step 2, kept redundant (the state project_out reads) for the whole run
            complex_function_1d pw = complex_factory_1d(world_).functor(
                [q](const coord_1d& r) { return std::exp(double_complex(0.0, q * r[0])); });
            planewave_conj_.push_back(conj(pw));
            planewave_conj_.back().change_tree_state(redundant);
            planewave_win_.push_back(complex_factory_1d(world_).functor([q, yw, rt](const coord_1d& r) {
                return 0.5 * std::erfc((std::abs(r[0]) - yw) / rt) * std::exp(double_complex(0.0, q * r[0]));
            }));
        }
        if (p_.verbose && world_.rank() == 0)
            print("MagneticKernelFT: B =", B, " cell [", lo_, ",", hi_, "]  period =", period_, " modes =", q_.size(),
                  " qmax =", q_.back(), " rcut =", rcut_, " rtaper =", rtaper_);
    }

    bool absorbs_paramagnetic() const override { return true; }
    std::string name() const override { return "ft"; }

    complex_function_3d apply(const complex_function_3d& f, const double E) const override {
        MADNESS_CHECK_THROW(E < bottom(), "MagneticKernelFT: require E < |B|/2");
        const double gap = bottom() - E;
        const magnetic::TauQuadrature quad = magnetic::trapezoid_log(p_.ft_h, p_.ft_tmin, p_.ft_tmax_factor / gap);
        std::vector<Node> nodes;
        for (const auto& qn : quad.nodes) {
            Node nd = node(qn.tau);
            // weight of the normalized kernel sqrt(c/pi) sqrt(1/(2 pi tau)), in log form
            nd.log_w = qn.log_dw + E * qn.tau + 0.5 * (nd.log_c - std::log(constants::pi))
                     + 0.5 * (-std::log(2 * qn.tau) - std::log(constants::pi));
            if (nd.log_w > -60.0) nodes.push_back(nd);
        }
        if (p_.verbose && world_.rank() == 0)
            print("MagneticKernelFT: E =", E, " nodes =", nodes.size(), "of", quad.nodes.size());
        return apply_nodes(f, nodes, quad.tmin);
    }

    const std::vector<double>& qgrid() const { return q_; }

private:
    struct Node {
        double tau, m, b, d, log_c, log_w;
    };

    static double tanh_over_x(double x) { return x == 0.0 ? 1.0 : std::tanh(x) / x; }

    /// mask exponents m, b, d and log c of node tau (finite at B = 0)
    Node node(const double tau) const {
        const double B = p_.B, x = 0.5 * B * tau;
        Node nd;
        nd.tau = tau;
        nd.m = magnetic::mehler_m(B, tau);
        nd.b = std::tanh(x);
        nd.d = 0.25 * tau * tanh_over_x(x);
        nd.log_c = magnetic::mehler_log_c(B, tau);
        nd.log_w = 0.0;
        return nd;
    }

    /// step 2 for the modes j0..j1-1: F_q(x, z) = int dy e^{-iqy} fxzy(x, z, y), as Function::project_out
    /// (fxzy reconstructed, e^{-iqy} redundant), with one fence for all modes
    std::vector<complex_function_2d> project_out_y(const complex_function_3d& fxzy, std::size_t j0, std::size_t j1) const {
        std::vector<complex_function_2d> F;
        for (std::size_t j = j0; j < j1; ++j)
            F.push_back(FunctionFactory<double_complex, 2>(world_).k(planewave_conj_[j].k())
                        .thresh(planewave_conj_[j].thresh()).nofence());
        world_.gop.fence();
        for (std::size_t j = j0; j < j1; ++j)
            fxzy.get_impl()->project_out(F[j - j0].get_impl().get(), planewave_conj_[j].get_impl().get(), 1, false);
        world_.gop.fence();
        for (auto& r : F) {
            r.get_impl()->trickle_down(false);
            r.get_impl()->set_tree_state(reconstructed);
        }
        world_.gop.fence();
        return F;
    }

    complex_function_3d apply_nodes(const complex_function_3d& f, const std::vector<Node>& nodes, const double s0) const {
        const double thresh = FunctionDefaults<3>::get_thresh();
        const double fq_tol = (p_.ft_fq_tol > 0 ? p_.ft_fq_tol : 0.01 * thresh) * f.norm2();
        std::vector<std::shared_ptr<operator2dT>> ops;
        // unnormalized exp(-c dx^2 - dz^2/(2 tau)) on the (x, z) plane
        for (const Node& nd : nodes)
            ops.push_back(magnetic::gauss_product_operator<2>(world_, {nd.log_c, -std::log(2 * nd.tau)}));

        // step 1, and (x, y, z) -> (x, z, y)
        complex_function_3d ft = copy(f);
        magnetic::mul_bilinear_phase_inplace(ft, phase_.with_kappa(-0.5 * p_.B));
        complex_function_3d fxzy = mapdim(ft, std::vector<long>{0, 2, 1});
        fxzy.reconstruct();

        // the modes are processed in batches: each step below is one vector operation over the batch
        std::vector<complex_function_2d> left;
        std::vector<complex_function_1d> right;
        long napplied = 0, nskip = 0;
        const double wall0 = wall_time();
        double wall_proj = 0.0;
        const std::size_t nq = q_.size();
        const std::size_t batch = p_.ft_qbatch > 0 ? std::size_t(p_.ft_qbatch) : nq;
        for (std::size_t j0 = 0; j0 < nq; j0 += batch) {
            const std::size_t j1 = std::min(nq, j0 + batch);
            // step 2: F_q(x, z) = int dy e^{-iqy} ft   (project_out does not conjugate)
            const double wp = wall_time();
            std::vector<complex_function_2d> Fall = project_out_y(fxzy, j0, j1);
            const std::vector<double> nall = norm2s(world_, Fall);
            wall_proj += wall_time() - wp;
            std::vector<std::size_t> jq;
            std::vector<complex_function_2d> F;
            std::vector<double> nF;
            for (std::size_t i = 0; i < Fall.size(); ++i) {
                if (nall[i] < fq_tol) continue;
                jq.push_back(j0 + i);
                F.push_back(Fall[i]);
                nF.push_back(nall[i]);
            }
            Fall.clear();
            if (F.empty()) continue;
            // step 3
            std::vector<complex_function_2d> Psi(F.size());
            for (std::size_t i = 0; i < F.size(); ++i) Psi[i] = mul(s0, F[i], false);
            world_.gop.fence();
            compress(world_, Psi);
            for (std::size_t k = 0; k < nodes.size(); ++k) {
                const double m = nodes[k].m, b = nodes[k].b, d = nodes[k].d;
                std::vector<real_function_2d> g;
                std::vector<complex_function_2d> Fk, Pk;
                for (std::size_t i = 0; i < F.size(); ++i) {
                    const double q = q_[jq[i]];
                    auto expo = [m, b, d, q](double x) { return m * x * x + b * q * x + d * q * q; };
                    double emin = std::min(expo(lo_), expo(hi_));
                    if (m > 0) {
                        const double xs = -b * q / (2 * m);
                        if (xs > lo_ && xs < hi_) emin = std::min(emin, expo(xs));
                    }
                    // ||g K g F|| <= max|g|^2 ||K||_1 ||F||, ||K||_1 = sqrt(pi/c) sqrt(2 pi tau)
                    const double log_bound = nodes[k].log_w + 0.5 * (std::log(constants::pi) - nodes[k].log_c)
                                           + 0.5 * std::log(2 * constants::pi * nodes[k].tau) - 2.0 * emin + std::log(nF[i]);
                    if (log_bound < std::log(fq_tol / nodes.size())) { ++nskip; continue; }
                    ++napplied;
                    g.push_back(real_factory_2d(world_).functor(
                        [expo](const coord_2d& r) { return std::exp(-expo(r[0])); }).nofence());
                    Fk.push_back(F[i]);
                    Pk.push_back(Psi[i]);     // shares the implementation: gaxpy below updates Psi[i]
                }
                if (g.empty()) continue;
                world_.gop.fence();
                std::vector<complex_function_2d> t = mul(world_, g, Fk);
                t = madness::apply(world_, *ops[k], t);
                t = mul(world_, g, t);
                compress(world_, t);
                gaxpy(world_, 1.0, Pk, std::exp(nodes[k].log_w), t);
            }
            truncate(world_, Psi);
            for (std::size_t i = 0; i < Psi.size(); ++i) {
                left.push_back(Psi[i]);
                right.push_back(planewave_win_[jq[i]]);
            }
        }
        if (p_.verbose && world_.rank() == 0)
            print("MagneticKernelFT: kept", left.size(), "of", q_.size(), "modes, applied", napplied,
                  "(node, q) terms, skipped", nskip, " wall: project", wall_proj, "s, modes", wall_time() - wall0, "s");
        // step 4, and (x, z, y) -> (x, y, z)
        complex_function_3d rxzy = hartree_product(left, right);
        complex_function_3d result = mapdim(rxzy, std::vector<long>{0, 2, 1});
        result.scale(1.0 / period_);
        magnetic::mul_bilinear_phase_inplace(result, phase_.with_kappa(0.5 * p_.B));
        return result;
    }

    World& world_;
    double lo_, hi_, period_, rcut_, rtaper_;
    std::vector<double> q_;
    std::vector<complex_function_1d> planewave_conj_, planewave_win_;
    magnetic::BilinearPhase phase_;
};

}  // namespace madness

#endif  // MADNESS_MRA_MAGNETIC_KERNEL_FT_H__INCLUDED
