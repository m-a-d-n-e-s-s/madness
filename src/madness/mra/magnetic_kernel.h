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
#include <madness/mra/QCCalculationParametersBase.h>

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace madness {

enum class MagneticGauge { symmetric, landau };

/// a nucleus as seen by the derived kernel parameters
struct MagneticKernelNucleus {
    double Z = 0.0;     ///< nuclear charge
    double x = 0.0;     ///< x coordinate (the ft y-wavenumbers carry its phase |B x|/2)
    double c = 1.0;     ///< smoothing radius of the nuclear potential
    double rho0 = 0.0;  ///< 1s density at the nucleus in the field, Z^3 g(|B|/Z^2)/pi
};

/// parameters of all magnetic kernels, input block "magnetic_kernel"

/// The kind and the field are not parameters of the kernel but of the physical problem (znemo: kernel and
/// physical_B of the "complex" block): make_magnetic_kernel takes them, set_derived_values the field.
/// set_derived_values sets all derived parameters (precedence: user-defined > derived > default).
struct MagneticKernelParameters : public QCCalculationParametersBase {
    static constexpr char const* tag = "magnetic_kernel";

    bool verbose = false;

    MagneticKernelParameters() {
        // proper-time quadrature of all kernels: derived from lo and eps, as the Coulomb and BSH fits (gfit.h)
        initialize<std::vector<double>>("origin", {0.0, 0.0, 0.0}, "gauge origin; only x, y matter (ft: 0)");
        initialize<double>("lo", 1.e-4, "smallest length scale: the smallest nuclear smoothing radius (derived)");
        initialize<double>("eps", 1.e-7, "precision of (T_k - E)^{-1}: 1e-3 thresh (derived)");
        initialize<double>("tau_step", 0.5, "step h of the trapezoid rule in ln tau (derived)");
        initialize<double>("tmin", 2.5e-9, "proper-time cut-off: int_0^tmin as tmin f (derived; nemo: adapted)");
        initialize<double>("t_lll", 1.e300, "chain: only the lowest Landau level beyond t_lll (derived)");
        initialize<double>("max_btau", 2.0, "chain: sub-steps with |B| dtau <= max_btau");
        initialize<double>("phase_tolfac", 0.1, "chain, ft: the phase product errs by <= phase_tolfac x the truncation "
                "tolerance per box (10x tighter costs ~1/3 refinement level at k = 10: uncritical)");
        initialize<bool>("phase_recursive", true, "chain, ft: recursive phase refinement");
        initialize<int>("nbroaden", 0, "chain, mehler: broaden passes before each convolution");
        // ft
        initialize<double>("ft_qextra", -1.0, "ft: |q| <= |B| L/2 + qextra; <0: from the nuclear cusps (ft_qmax)");
        initialize<double>("ft_qmax", 0.0, "ft: y-wavenumber window |q| <= qmax (derived)");
        initialize<double>("ft_rcut", 1.e10, "ft: radius of the cylinder holding f and G f (derived)");
        initialize<double>("ft_rtaper", 0.5, "ft: edge width of the y-window (derived)");
        initialize<double>("ft_period", 0.0, "ft: y-period (derived)");
        initialize<double>("ft_fq_tol", 1.e-6, "ft: drop modes with ||F_q|| < fq_tol ||f|| (derived)");
        initialize<long>("ft_qbatch", 64, "ft: q-modes processed together; 0: all");
    }

    MagneticKernelParameters(World& world, const commandlineparser& parser) : MagneticKernelParameters() {
        read_input_and_commandline_options(world, parser, tag);
    }

    MagneticKernelParameters(const MagneticKernelParameters& other) = default;
    MagneticKernelParameters& operator=(const MagneticKernelParameters& other) = default;

    std::string get_tag() const override { return tag; }

    /// set all derived parameters; requires the cubic cell of FunctionDefaults<3>

    /// \param B           field along z
    /// \param thresh      truncation threshold of the orbitals
    /// \param nuclei      the nuclei; none: lo as set, the ft window from qextra = 10 (smooth f)
    /// \param cusp_free   the right-hand side carries no nuclear cusp (nemo, F = psi/R)
    void set_derived_values(const double B, const double thresh,
                            const std::vector<MagneticKernelNucleus>& nuclei = {}, const bool cusp_free = false) {
        const double absB = std::abs(B);

        // lo = min_A c_A, capped at 1
        if (not nuclei.empty()) {
            double l = 1.0;
            for (const auto& n : nuclei) l = std::min(l, n.c);
            set_derived_value("lo", l);
        }
        // eps = 1e-3 thresh, as the BSH operators
        set_derived_value("eps", 1.e-3 * thresh);
        // h: per spectral component kappa the integrand in u = ln tau is e^{u - kappa e^u}; by Poisson summation the
        // relative error is 2 sum_k |Gamma(1 + 2 pi i k/h)| ~ (4 pi/sqrt h) e^{-pi^2/h} = eps, for all kappa
        double h = 0.5;
        for (int i = 0; i < 5; ++i) h = constants::pi * constants::pi / std::log(4.0 * constants::pi / (eps() * std::sqrt(h)));
        set_derived_value("tau_step", h);
        // tmin = lo^2/4: tmin f shifts the energy by -1.79 (tmin/lo^2)^2 (~ -10%) of the nuclear smoothing error of
        // radius lo (both scale as Z^2 rho(0) lo^3); a cusp-free right-hand side allows more, tmin_cusp_free
        set_derived_value("tmin", 0.25 * lo() * lo());
        // chain: beyond t_lll = -ln(eps)/|B| the excited Landau levels, decaying as e^{-|B| tau} relative to the
        // lowest, are below eps; only the free z motion is propagated
        set_derived_value("t_lll", absB > 0.0 ? -std::log(eps()) / absB : 1.e300);

        // ft: f and G f decay as e^{-|B| rho^2/4}: cylinder rcut = 1.5 sqrt(8 ln(1/thresh)/|B|), edge width
        // min(0.5, rcut/4)
        set_derived_value("ft_rcut", absB > 0.0 ? 1.5 * std::sqrt(8.0 * std::log(1.0 / thresh) / absB) : 1.e10);
        set_derived_value("ft_rtaper", std::min(0.5, ft_rcut() / 4.0));
        // ft lives on |x|, |y| < rext = rcut + 3 rtaper: the period P = min(2 x cell width, 4 rext) keeps the
        // y-images of ft and G ft apart; y-wavenumbers up to |B| L/2, L = min(cell half width, rext)
        const Tensor<double> cell = FunctionDefaults<3>::get_cell();
        const double rext = ft_rcut() + 3.0 * ft_rtaper();
        const double L = std::min(std::max(std::abs(cell(0, 0)), std::abs(cell(0, 1))), rext);
        set_derived_value("ft_period", std::min(2.0 * (cell(0, 1) - cell(0, 0)), 4.0 * rext));
        // plus those of f: qextra >= 0 as given; else the cusp momenta of the right-hand side at each nucleus,
        // shifted by the phase |B x_A|/2, with S = sum_A n_A Z_A^2 rho_A(0), n_A = min(2, Z_A):
        //  - bare V psi: the cut |q| > Q raises the energy by (4/3) S/Q^3; set to a third of the ft precision floor
        //    79 uEh S^0.75 (thresh/1e-4)^0.69 (fit, H at B = 1-100, thresh 1e-4-1e-5; errors add in quadrature)
        //  - cusp-free U F: the leading non-smooth term U1.grad F ~ Z r^.d costs (32/45) Z^2 |d|^2/Q^5, |d|^2 <=
        //    2 n_A rho_A(0) (p-like, worst case): Q = [64 S/(45 dE)]^(1/5), dE = 0.1 thresh
        double qmax = 0.5 * absB * L + (ft_qextra() >= 0.0 ? ft_qextra() : (nuclei.empty() ? 10.0 : 0.0));
        if (ft_qextra() < 0.0 and not nuclei.empty()) {
            double S = 0.0, xmax = 0.0;
            for (const auto& n : nuclei) {
                S += std::min(2.0, n.Z) * n.Z * n.Z * n.rho0;
                xmax = std::max(xmax, std::abs(n.x));
            }
            const double floor = 79.e-6 * std::pow(S, 0.75) * std::pow(thresh / 1.e-4, 0.69);
            const double Q = cusp_free ? std::pow(64.0 * S / (45.0 * 0.1 * thresh), 0.2) : std::cbrt(4.0 * S / floor);
            qmax = std::max(qmax, Q + 0.5 * absB * xmax);
        }
        set_derived_value("ft_qmax", qmax);
        // modes with ||F_q|| < 0.01 thresh ||f|| are dropped
        set_derived_value("ft_fq_tol", 0.01 * thresh);
    }

    /// tmin for a cusp-free right-hand side g: int_0^tmin e^{-tau (T-E)} g ~ tmin g shifts the energy by
    /// -(tmin^2/2) C, C = <g|T - E|g> finite (it diverges with the cusp): tmin = sqrt(2 eps/C), eps = 1e-3 thresh
    /// (H, B = 100, chain, thresh 1e-4: 0.1 / 0.01 / 0.001 thresh -> E_k -257 / +19.7 / +4.6 uEh, bare +2.2)
    double tmin_cusp_free(const double C) const { return std::sqrt(2.0 * eps() / std::max(C, 1.e-10)); }

    coord_3d origin() const {
        const auto o = get<std::vector<double>>("origin");
        MADNESS_CHECK_THROW(o.size() == 3, "MagneticKernelParameters: origin needs 3 coordinates");
        return coord_3d(o);
    }
    double lo() const { return get<double>("lo"); }
    double eps() const { return get<double>("eps"); }
    double tau_step() const { return get<double>("tau_step"); }
    double tmin() const { return get<double>("tmin"); }
    double t_lll() const { return get<double>("t_lll"); }
    double max_btau() const { return get<double>("max_btau"); }
    double phase_tolfac() const { return get<double>("phase_tolfac"); }
    bool phase_recursive() const { return get<bool>("phase_recursive"); }
    int nbroaden() const { return get<int>("nbroaden"); }
    double ft_qextra() const { return get<double>("ft_qextra"); }
    double ft_qmax() const { return get<double>("ft_qmax"); }
    double ft_rcut() const { return get<double>("ft_rcut"); }
    double ft_rtaper() const { return get<double>("ft_rtaper"); }
    double ft_period() const { return get<double>("ft_period"); }
    double ft_fq_tol() const { return get<double>("ft_fq_tol"); }
    long ft_qbatch() const { return get<long>("ft_qbatch"); }

    /// set a parameter from its string representation (user-defined precedence; Python bindings, tests)
    void set_user_value(const std::string& key, const std::string& value) {
        MADNESS_CHECK_THROW(parameter_exists(key), ("MagneticKernelParameters: unknown parameter " + key).c_str());
        MADNESS_CHECK_THROW(try_setting_any<all_parameter_types>(key, value),
                            ("MagneticKernelParameters: cannot set " + key + " to " + value).c_str());
    }
};

class MagneticKernel {
public:
    MagneticKernel(const double B, const MagneticKernelParameters& p) : B_(B), p_(p) {}
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
    double bottom() const { return 0.5 * std::abs(B_); }
    /// field along z
    double B() const { return B_; }
    /// is (B/2) Lz part of T_k (else it belongs on the right-hand side)?
    virtual bool absorbs_paramagnetic() const = 0;
    /// is (B^2/8) rho^2 part of T_k?
    virtual bool absorbs_diamagnetic() const { return true; }
    virtual MagneticGauge gauge() const { return MagneticGauge::symmetric; }
    virtual std::string name() const = 0;

    const MagneticKernelParameters& parameters() const { return p_; }

    /// vector potential of the kernel's gauge at r
    std::array<double, 3> vector_potential(const coord_3d& r) const {
        const coord_3d o = p_.origin();
        const double xt = r[0] - o[0], yt = r[1] - o[1];
        if (gauge() == MagneticGauge::landau) return {-B_ * yt, 0.0, 0.0};
        return {-0.5 * B_ * yt, 0.5 * B_ * xt, 0.0};
    }

    /// e^{+i chi} f (to the kernel's gauge from the symmetric one), chi = (B/2) x~ y~; identity for the
    /// symmetric gauge.  sign = -1: back to the symmetric gauge.  Untruncated.
    complex_function_3d from_symmetric_gauge(const complex_function_3d& f, const double sign = 1.0) const {
        if (gauge() == MagneticGauge::symmetric || B_ == 0.0) return copy(f);
        magnetic::BilinearPhase ph(sign * 0.5 * B_, p_.origin()[0], p_.origin()[1]);
        ph.tolfac = p_.phase_tolfac();
        ph.recursive = p_.phase_recursive();
        return magnetic::mul_bilinear_phase(f, ph);
    }
    complex_function_3d to_symmetric_gauge(const complex_function_3d& f) const { return from_symmetric_gauge(f, -1.0); }

    /// (T_k - E)^{-1} = int_0^inf e^{-tau (T_k - E)} dtau as a trapezoid rule in u = ln tau on [tmin, tmax], step
    /// tau_step (both from MagneticKernelParameters::set_derived_values), [0, tmin] as tmin f:
    ///  - the nodes cut below ln tmin sum to tmin (h/2) coth(h/2): the weight of tmin f
    ///  - tmax = -ln(eps)/(|B|/2 - E): the slowest component has decayed to eps
    magnetic::TauQuadrature tau_quadrature_to(const double tmax) const {
        const double h = p_.tau_step();
        magnetic::TauQuadrature q = magnetic::trapezoid_log(h, p_.tmin(), tmax);
        q.tmin *= 0.5 * h / std::tanh(0.5 * h);
        return q;
    }
    magnetic::TauQuadrature tau_quadrature(const double gap) const { return tau_quadrature_to(-std::log(p_.eps()) / gap); }

protected:
    double B_;
    MagneticKernelParameters p_;
};

}  // namespace madness

#endif  // MADNESS_MRA_MAGNETIC_KERNEL_H__INCLUDED
