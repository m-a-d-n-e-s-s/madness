/// \file magnetic_kernel_chain.h
/// \brief Real-space magnetic chain: (T_B - E)^{-1} from the exact splitting of the heat kernel
///
///   e^{-tau T_perp} = e^{-a pi_y^2} e^{-s pi_x^2} e^{-a pi_y^2},   a = tanh(B tau/2)/(2B), s = sinh(B tau)/(2B),
///
/// with each e^{-c pi^2} a Gaussian convolution conjugated by a phase, and the free z motion merged into
/// the x convolution:
///
///   e^{-tau T} f = e^{-iBxy/2} G_y(2a) [e^{+iBxy} G_xz(2s, tau) [e^{-iBxy} G_y(2a) [e^{+iBxy/2} f]]].
///
/// The bracket without its outer pair e^{-+iBxy/2} is the Landau-gauge propagator (chi = (B/2) xy), so
/// the Landau kernel is this kernel without the outer phases.  G_E = tmin + int_tmin^tmax dtau e^{E tau}
/// e^{-tau T} on trapezoid nodes in ln tau (MagneticKernel::tau_quadrature), the chain propagated incrementally
/// from node to node, sub-steps with |B| dtau <= max_btau, the Landau zero point e^{|B| tau/2} carried by the propagated
/// function, and only the free z motion beyond tau = -ln(eps)/|B|.
///
/// MADNESS layout: internal axis order (y, x, z), so the y convolution is a 1D low-dimensional apply on
/// axis 0 and the x+z convolution a 2D one on axes (1, 2); input and output are in (x, y, z).
/// C++ port of MagneticGreensFunction3D / MagneticHeatKernel3D (hydrogen_magnetic_pymadness.py) and of the
/// Landau-gauge Gop of h2p_landau.py; see magnetic_kernels.md.

#ifndef MADNESS_MRA_MAGNETIC_KERNEL_CHAIN_H__INCLUDED
#define MADNESS_MRA_MAGNETIC_KERNEL_CHAIN_H__INCLUDED

#include <madness/mra/magnetic_kernel.h>

namespace madness {

class MagneticKernelChain : public MagneticKernel {
public:
    MagneticKernelChain(World& world, const double B, const MagneticKernelParameters& p, const MagneticGauge gauge)
        : MagneticKernel(B, p), world_(world), gauge_(gauge) {
        magnetic::set_lowdim_defaults_from_3d();
        // lowest-Landau-level switch: relative to the lowest Landau level the excited ones decay as e^{-|B| n tau},
        // so beyond t_lll = -ln(eps)/|B| they are below eps and the transverse propagator is the LLL projector
        // (times the zero point e^{-|B| tau/2}, carried by S): only the free z motion is propagated
        t_lll_ = p_.t_lll();
    }

    bool absorbs_paramagnetic() const override { return true; }
    MagneticGauge gauge() const override { return gauge_; }
    std::string name() const override { return gauge_ == MagneticGauge::landau ? "chain_landau" : "chain_symmetric"; }

    complex_function_3d apply(const complex_function_3d& f, const double E) const override {
        const double B = B_;
        const double gap = bottom() - E;
        MADNESS_CHECK_THROW(gap > 0.0, "MagneticKernelChain: require E < |B|/2");
        const magnetic::TauQuadrature quad = tau_quadrature(gap);
        const std::vector<long> swap_xy{1, 0, 2};
        const complex_function_3d a = mapdim(f, swap_xy);

        const bool sym = gauge_ == MagneticGauge::symmetric;
        complex_function_3d out0 = quad.tmin * a;
        complex_function_3d S = sym ? mulphase(a, 0.5 * B) : a;
        complex_function_3d out;
        double tprev = 0.0;
        for (const auto& nd : quad.nodes) {
            const double t = nd.tau;
            // S carries the Landau zero point e^{|B| tau/2}
            const double w = std::exp(nd.log_dw + (E - bottom()) * t);
            // full chain up to t_lll, then the free z motion only (the z heat kernel over t - tm)
            const double tm = std::min(std::max(tprev, t_lll_), t);
            if (tm > tprev) S = inner_chain(S, tm - tprev, true);
            if (t > tm) S = conv_z(S, t - tm);
            accumulate(out, S, w);
            tprev = t;
        }
        out.truncate(thresh());
        complex_function_3d result = sym ? out0 + mulphase(out, -0.5 * B) : out0 + out;
        result.truncate(thresh());
        return mapdim(result, swap_xy);
    }

    /// e^{+iBxy/2} e^{-tau T} e^{-iBxy/2} a in the internal order (the chain without its outer phases);
    /// rescale: times e^{|B| tau/2}, per sub-step
    complex_function_3d inner_chain(const complex_function_3d& a, const double tau, const bool rescale) const {
        const double B = B_;
        const long nstep = std::max(1L, long(std::ceil(std::abs(B) * tau / p_.max_btau() - 1e-12)));
        const auto [A, S] = magnetic::splitting_coefficients(B, tau / nstep);
        complex_function_3d f = conv_y(a, 2 * A);
        for (long istep = 0; istep < nstep; ++istep) {
            f = mulphase(f, -B);
            f = conv_xz(f, 2 * S, tau / nstep);
            f = mulphase(f, B);
            f = conv_y(f, istep == nstep - 1 ? 2 * A : 4 * A);
            if (rescale) f.scale(std::exp(0.5 * std::abs(B) * tau / nstep));
        }
        return f;
    }

private:
    double thresh() const { return FunctionDefaults<3>::get_thresh(); }

    /// a e^{i kappa x y}, truncated (internal order)
    complex_function_3d mulphase(const complex_function_3d& a, const double kappa) const {
        if (kappa == 0.0) return a;
        complex_function_3d r = magnetic::mul_bilinear_phase(a, phase_.with_kappa(kappa));
        r.truncate(thresh());
        return r;
    }
    complex_function_3d conv_y(const complex_function_3d& a, const double var) const {
        complex_function_3d r = magnetic::gauss_conv_1d(a, 0, var);
        r.truncate(thresh());
        return r;
    }
    complex_function_3d conv_z(const complex_function_3d& a, const double var) const {
        complex_function_3d r = magnetic::gauss_conv_1d(a, 2, var);
        r.truncate(thresh());
        return r;
    }
    complex_function_3d conv_xz(const complex_function_3d& a, const double var_x, const double var_z) const {
        complex_function_3d r = magnetic::gauss_conv_2d(a, 1, var_x, var_z);
        r.truncate(thresh());
        return r;
    }

    /// out <- out + w S: in place if both are compressed, else out of place (S unchanged)
    static void accumulate(complex_function_3d& out, const complex_function_3d& S, const double w) {
        if (!out.is_initialized()) {
            out = w * S;
            return;
        }
        if (out.is_compressed() && S.is_compressed()) {
            out.gaxpy(1.0, S, w);
            return;
        }
        out = out + w * S;
    }

    World& world_;
    MagneticGauge gauge_;
    double t_lll_;      ///< lowest-Landau-level switch -ln(eps)/|B|
    magnetic::BilinearPhase phase_;
};

}  // namespace madness

#endif  // MADNESS_MRA_MAGNETIC_KERNEL_CHAIN_H__INCLUDED
