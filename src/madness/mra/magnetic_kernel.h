/// \file magnetic_kernel.h
/// \brief Common interface of the Green's-function kernels for an electron in a uniform magnetic field
///
/// Conventions (all kernels): atomic units, B along +z, pi = p + A, gauge origin O = (ox, oy),
///
///   T_B = 1/2 (p + A)^2 = -1/2 nabla^2 + (B/2) Lz + (B^2/8) rho^2      (symmetric gauge, A = (B/2)(-y~, x~, 0))
///
/// with x~ = x - ox, y~ = y - oy.  Each kernel applies (T_k - E)^{-1} for E < bottom() = |B|/2, where
///   - ft, chain_symmetric:  T_k = T_B in the symmetric gauge,
///   - chain_landau:         T_k = T_B in the Landau gauge A_L = (-B y~, 0, 0) = A - grad chi, chi = (B/2) x~ y~,
///                           acting on phi_L = e^{+i chi} phi_sym,
///   - mehler:               T_k = -1/2 nabla^2 + (B^2/8) rho^2  (the paramagnetic (B/2) Lz stays on the right-hand side).
/// Theory, provenance and parameters: magnetic_kernels.md.

#ifndef MADNESS_MRA_MAGNETIC_KERNEL_H__INCLUDED
#define MADNESS_MRA_MAGNETIC_KERNEL_H__INCLUDED

#include <madness/mra/magnetic_primitives.h>

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace madness {

enum class MagneticGauge { symmetric, landau };

/// parameters of all magnetic kernels; the defaults are those of the original implementations
struct MagneticKernelParameters {
    std::string kind = "chain_symmetric";   ///< ft | chain_symmetric | chain_landau | mehler
    double B = 0.0;                         ///< field along z
    double ox = 0.0, oy = 0.0;              ///< gauge origin (the FT kernel requires 0, 0)
    bool verbose = false;

    // chain and mehler: Gauss-Legendre in ln tau on [tmin, -ln(eps)/(|B|/2 - E)]
    int nodes = 40;
    double tmin = 1e-4;
    double eps = 1e-10;
    double nlll = 28.0;         ///< chain: beyond nlll/|B| only the free z motion is propagated (0: off)
    double max_btau = 2.0;      ///< chain: sub-steps with |B| dtau <= max_btau
    double phase_tolfac = 0.1;  ///< chain: phase refinement tolerance factor
    bool phase_recursive = true;
    int nbroaden = 0;           ///< chain, mehler: broaden passes before each convolution (narrow kernels couple
                                ///  neighbours through edge terms that cancel only between boxes of equal level)

    // ft: trapezoid in ln tau with step ft_h on [ft_tmin, ft_tmax_factor/(|B|/2 - E)]
    double ft_h = 0.25;
    double ft_tmin = 1e-8;
    double ft_tmax_factor = 45.0;
    double ft_period = 0.0;     ///< y-period P; 0 -> min(2 x cell width, 4 rext)
    double ft_qextra = 10.0;    ///< |q| <= |B| L/2 + qextra
    double ft_rcut = 0.0;       ///< phase cut off (erfc taper) at rho > rcut; 0 -> 1.5 sqrt(8 ln(1/thresh)/|B|)
    double ft_rtaper = 0.0;     ///< taper width; 0 -> min(0.5, rcut/4)
    double ft_dphi = 3.0;       ///< refine until the phase varies by <= dphi per box
    double ft_fq_tol = 0.0;     ///< drop modes with ||F_q|| < fq_tol ||f||; 0 -> 0.01 thresh
    long ft_qbatch = 64;        ///< q-modes processed together; 0 -> all
};

class MagneticKernel {
public:
    explicit MagneticKernel(const MagneticKernelParameters& p) : p_(p) {}
    virtual ~MagneticKernel() = default;

    /// (T_k - E)^{-1} f, E < bottom()
    virtual complex_function_3d apply(const complex_function_3d& f, double E) const = 0;

    /// (T_k - E_i)^{-1} f_i
    std::vector<complex_function_3d> apply(const std::vector<complex_function_3d>& f, const Tensor<double>& E) const {
        MADNESS_CHECK_THROW(long(f.size()) == E.size(), "MagneticKernel::apply: one energy per function");
        std::vector<complex_function_3d> result;
        for (std::size_t i = 0; i < f.size(); ++i) result.push_back(apply(f[i], E(long(i))));
        return result;
    }

    /// inf spec T_k = |B|/2
    double bottom() const { return 0.5 * std::abs(p_.B); }
    /// is (B/2) Lz part of T_k (else it belongs on the right-hand side)?
    virtual bool absorbs_paramagnetic() const = 0;
    /// is (B^2/8) rho^2 part of T_k?
    virtual bool absorbs_diamagnetic() const { return true; }
    virtual MagneticGauge gauge() const { return MagneticGauge::symmetric; }
    virtual std::string name() const = 0;

    const MagneticKernelParameters& parameters() const { return p_; }

    /// vector potential of the kernel's gauge at r
    std::array<double, 3> vector_potential(const coord_3d& r) const {
        const double xt = r[0] - p_.ox, yt = r[1] - p_.oy;
        if (gauge() == MagneticGauge::landau) return {-p_.B * yt, 0.0, 0.0};
        return {-0.5 * p_.B * yt, 0.5 * p_.B * xt, 0.0};
    }

    /// e^{+i chi} f (to the kernel's gauge from the symmetric one), chi = (B/2) x~ y~; identity for the
    /// symmetric gauge.  sign = -1: back to the symmetric gauge.  Untruncated.
    complex_function_3d from_symmetric_gauge(const complex_function_3d& f, const double sign = 1.0) const {
        if (gauge() == MagneticGauge::symmetric || p_.B == 0.0) return copy(f);
        magnetic::BilinearPhase ph(sign * 0.5 * p_.B, p_.ox, p_.oy);
        ph.tolfac = p_.phase_tolfac;
        ph.recursive = p_.phase_recursive;
        return magnetic::mul_bilinear_phase(f, ph);
    }
    complex_function_3d to_symmetric_gauge(const complex_function_3d& f) const { return from_symmetric_gauge(f, -1.0); }

protected:
    MagneticKernelParameters p_;
};

}  // namespace madness

#endif  // MADNESS_MRA_MAGNETIC_KERNEL_H__INCLUDED
