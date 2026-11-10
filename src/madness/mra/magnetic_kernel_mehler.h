/// \file magnetic_kernel_mehler.h
/// \brief Mehler kernel: (T_HO - E)^{-1} with T_HO = -1/2 nabla^2 + (B^2/8) rho~^2, the diamagnetic part of T_B
///
/// T_HO is a 2D harmonic oscillator of frequency w = |B|/2 in (x~, y~) plus the free z motion; its heat
/// kernel factorizes exactly at rank 1 (Mehler),
///
///   <r|e^{-tau T_HO}|r'> = (c/pi) (2 pi tau)^{-1/2} g(rho~) exp(-c |rho - rho'|^2 - (z-z')^2/(2 tau)) g(rho~'),
///   g = exp(-m rho~^2),  m = (w/2) tanh(w tau/2),  c = w/(2 sinh(w tau)),
///
/// so G_E = tmin + sum_n w_n g_n K_n g_n with one 3D Gaussian convolution K_n per node, Gauss-Legendre in
/// ln tau on [tmin, -ln(eps)/(w - E)] (the nodes of the chain kernel).  inf spec T_HO = w = |B|/2 = inf spec T_B.
/// The paramagnetic (B/2) Lz is not in T_HO: it stays on the right-hand side of the iteration.  On functions with
/// Lz = 0 about the gauge origin the Mehler kernel equals G_B.  3D version of MehlerGreensFunction2D
/// (landau_mehler_pymadness.py); see mehler_kernel_factorization.md and magnetic_kernels.md.

#ifndef MADNESS_MRA_MAGNETIC_KERNEL_MEHLER_H__INCLUDED
#define MADNESS_MRA_MAGNETIC_KERNEL_MEHLER_H__INCLUDED

#include <madness/mra/magnetic_kernel.h>

namespace madness {

class MagneticKernelMehler : public MagneticKernel {
public:
    MagneticKernelMehler(World& world, const MagneticKernelParameters& p) : MagneticKernel(p), world_(world) {}

    bool absorbs_paramagnetic() const override { return false; }
    std::string name() const override { return "mehler"; }

    complex_function_3d apply(const complex_function_3d& f, const double E) const override {
        const double w = bottom();
        const double gap = w - E;
        MADNESS_CHECK_THROW(gap > 0.0, "MagneticKernelMehler: require E < |B|/2");
        const magnetic::TauQuadrature quad = magnetic::gauss_legendre_log(p_.nodes, p_.tmin, -std::log(p_.eps) / gap);
        const double thresh = FunctionDefaults<3>::get_thresh();
        const double ox = p_.ox, oy = p_.oy;

        complex_function_3d out = quad.tmin * f;
        out.compress();
        for (const auto& nd : quad.nodes) {
            const double tau = nd.tau;
            const double m = magnetic::mehler_m(w, tau), log_c = magnetic::mehler_log_c(w, tau);
            // normalized kernel (c/pi) (2 pi tau)^{-1/2}, in log form
            const double log_w = nd.log_dw + E * tau + log_c - std::log(constants::pi) - 0.5 * std::log(2 * constants::pi * tau);
            if (log_w < -60.0) continue;
            const auto op = magnetic::gauss_product_operator<3>(world_, {log_c, log_c, -std::log(2 * tau)});
            real_function_3d g = real_factory_3d(world_).functor([m, ox, oy](const coord_3d& r) {
                const double xt = r[0] - ox, yt = r[1] - oy;
                return std::exp(-m * (xt * xt + yt * yt));
            });
            complex_function_3d t = g * f;
            t.truncate(thresh);
            magnetic::broaden_until_stable(t, p_.nbroaden);
            t = madness::apply(*op, t);
            t = g * t;
            t.truncate(thresh);
            t.compress();
            out.gaxpy(1.0, t, std::exp(log_w));
        }
        out.truncate(thresh);
        return out;
    }

private:
    World& world_;
};

}  // namespace madness

#endif  // MADNESS_MRA_MAGNETIC_KERNEL_MEHLER_H__INCLUDED
