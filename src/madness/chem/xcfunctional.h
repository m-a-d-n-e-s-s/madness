#ifndef MADNESS_CHEM_XCFUNCTIONAL_H__INCLUDED
#define MADNESS_CHEM_XCFUNCTIONAL_H__INCLUDED

#include <madness/madness_config.h>

MADNESS_PRAGMA_GCC(diagnostic push)
MADNESS_PRAGMA_GCC(diagnostic ignored "-Wcomment")
MADNESS_PRAGMA_CLANG(diagnostic push)
MADNESS_PRAGMA_CLANG(diagnostic ignored "-Wcomment")

/// \file moldft/xcfunctional.h
/// \brief Defines interface for DFT XC functionals
/// \ingroup chemistry

#include <madness/tensor/tensor.h>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <string>
#include <utility>
#include <madness/mra/key.h>
#include <madness/world/MADworld.h>
#include <madness/mra/function_common_data.h>
#include <madness/mra/function_interface.h>

#ifdef MADNESS_HAS_LIBXC
#include <xc.h>
#endif

namespace madness {

/// TEMPORARY investigation switches for the XC path (see XC-NEMO-MEMORY-INVESTIGATION.md)

/// All MRA functions carry noise at O(thresh); a derivative operator amplifies
/// that noise by 2^n at level n, so every numerical derivative on the XC path is
/// a noise amplifier. These switches select reformulations that keep the
/// derivative off the noisy/cusped object, in the same spirit as the existing
/// zeta = grad log(rho) representation. They exist so both forms can be measured
/// with one binary and are meant to disappear once the better form is settled.

/// MAD_XC_FACTORED_GGA: form the semilocal potential as
/// -div(rho B) = -rho (div B + 2 de/dsigma |zeta|^2), B = 2 de/dsigma zeta,
/// so the differentiated object carries no factor of rho and the noise of the
/// divergence is damped by rho instead of spread over the whole box.
inline bool xc_factored_gga_potential() {
    static const bool flag = []() {
        const char* e = std::getenv("MAD_XC_FACTORED_GGA");
        return (e != nullptr) and (std::string(e) != "0");
    }();
    return flag;
}

/// MAD_XC_NEMO_ZETA: with a nuclear correlation factor, build
/// zeta = grad log(rho) = grad log(rho_reg) - 2 U1 from the *regularized*
/// density and the analytic U1, so the nuclear cusp of rho = R^2 rho_reg never
/// goes under a numerical derivative.
inline bool xc_nemo_zeta() {
    static const bool flag = []() {
        const char* e = std::getenv("MAD_XC_NEMO_ZETA");
        return (e != nullptr) and (std::string(e) != "0");
    }();
    return flag;
}

/// MAD_XC_WEAK_GGA: abandon the multiplicative xc potential.

/// The semilocal potential is -div(X), X = 2 de/dsigma grad(rho), and X has a jump
/// discontinuity at every nucleus (zeta = grad log rho -> -2Z r_hat, whose
/// Cartesian components flip sign across the origin). Differentiating that jump is
/// what produces the +-8e4 excursions; no algebraic rearrangement of the
/// multiplicative form avoids it, because div(X) *is* a derivative of X.
///
/// The weak form never differentiates X:
///     <phi|v|psi> = int (df/drho) phi psi + int X . grad(phi psi)
/// and in a Green's-function code the same holds for the orbital update, because a
/// radial convolution commutes with the gradient:
///     G * (psi div X) = div(G * (X psi)) - G * (X . grad psi)
/// so the divergence acts on G*(X psi), which is C^1, and the jump is only ever
/// convolved. Same for the meta-gga term -1/2 div(v_tau grad psi).
inline bool xc_weak_gga() {
    static const bool flag = []() {
        const char* e = std::getenv("MAD_XC_WEAK_GGA");
        return (e != nullptr) and (std::string(e) != "0");
    }();
    return flag;
}

/// MAD_XC_ONELEVEL_DIV: algorithm A'' of XC-POTENTIAL-ALGORITHMS.md 1G.3.1.

/// Expand the semilocal divergence one level only,
///     div(X) = grad(u).grad(rho) + u laplacian(rho),   u = 2 de/dsigma
/// with laplacian(rho) = div(rho zeta) -- a single divergence of the
/// semi-analytic gradient, which is the route measured to conserve. grad(u) is
/// taken numerically from the projected u; the Hessian and libxc's second
/// derivatives are deliberately NOT used.
///
/// This does not remove the 1/r pole (u is finite at a nucleus while
/// laplacian(rho) ~ -4 Z rho(0)/r), but it isolates the pole in one term with a
/// known coefficient, which is what the analytic subtraction of 1G.4 needs. What
/// it costs is the structural int div(X) = 0 that the plain divergence has for
/// free -- and measuring that loss is the point of the switch. Spin-restricted
/// only.
inline bool xc_onelevel_div() {
    static const bool flag = []() {
        const char* e = std::getenv("MAD_XC_ONELEVEL_DIV");
        return (e != nullptr) and (std::string(e) != "0");
    }();
    return flag;
}

/// MAD_XC_A2_MUNGED_GRAD: force algorithm A'' to use the *munged* grad(rho) that
/// the reference flux is built from, instead of make_ddensity's unmunged one.

/// With this set, u*grad(rho) == X holds exactly, so any remaining difference
/// between A'' and algorithm A is A'''s own error rather than the two computing
/// slightly different operators in the tail.
/// MAD_XC_PLOT: if set, plot_plane() the assembled xc potential under this tag.
/// Written every time make_xc_potential runs, so the file left behind is the
/// converged one. The plotting window comes from a `plot ... end` block in a file
/// named `input` in the run directory (plot_plane's default). Diagnostic only.
inline std::string xc_plot_tag() {
    static std::string tag = []{
        const char* e = std::getenv("MAD_XC_PLOT");
        return std::string(e ? e : "");
    }();
    return tag;
}

inline bool xc_a2_munged_grad() {
    static const bool flag = []() {
        const char* e = std::getenv("MAD_XC_A2_MUNGED_GRAD");
        return (e != nullptr) and (std::string(e) != "0");
    }();
    return flag;
}

/// MAD_XC_FLUX_NUMDIV: in weak form, evaluate 2 div(G*Y) with a numerical ABGV
/// divergence of the Green's-function output instead of the analytic gradient of
/// the kernel.

/// The default is the analytic route, sum_a (d_a G) * Y_a, via
/// GradBSHOperator: differentiation of a convolution may be moved onto either
/// factor, and moving it onto the kernel means no numerical derivative appears
/// in the flux term at all. The numerical route differentiates a function that
/// madness computed to bshtol, and a derivative amplifies node-level error by
/// 2^n -- at thresh 1e-8 and depth 22 that is an absolute error of order 0.04 in
/// the orbital update, which is bounded (unlike the pole it replaced) but is the
/// leading suspect for the weak form needing more SCF iterations. Kept as a
/// switch so the two can be compared on one binary.
inline bool xc_flux_numerical_div() {
    static const bool flag = []() {
        const char* e = std::getenv("MAD_XC_FLUX_NUMDIV");
        return (e != nullptr) and (std::string(e) != "0");
    }();
    return flag;
}

/// MAD_XC_SMOOTH_VTAU: replace the hard binary_munge screen on de/dtau by the C^2
/// ramp of XCfunctional::gga_ramp. The hard screen puts a jump of up to O(10^2) on
/// the rho = ggatol iso-surface, and apply_tau_term differentiates exactly that
/// product -- a jump is unrepresentable in a multiwavelet basis and its derivative
/// refines without bound.
inline bool xc_smooth_vtau() {
    static const bool flag = []() {
        const char* e = std::getenv("MAD_XC_SMOOTH_VTAU");
        return (e != nullptr) and (std::string(e) != "0");
    }();
    return flag;
}

/// MAD_XC_FLUX_TRUNC: truncate the semilocal flux at this multiple of thresh before
/// taking its divergence (0 = off, the current behaviour).

/// refine_to_common_level lifts every xc_arg to tau's refinement level, so on the
/// nemo path the flux is represented at depth 20-22 where moldft has 10-13. Its
/// genuine content needs neither: the flux norm is identical at every depth. What
/// the extra levels carry is O(thresh) noise, and d/dx multiplies a level-n
/// coefficient by 2^n -- measured, div(flux)'s pointwise max grows like 2^depth
/// while its norm does not move at all.
inline double xc_flux_truncation() {
    static const double fac = []() {
        const char* e = std::getenv("MAD_XC_FLUX_TRUNC");
        return (e == nullptr) ? 0.0 : std::atof(e);
    }();
    return fac;
}

/// MAD_XC_TAU_NO_U1: DIAGNOSTIC ONLY, gives a wrong tau. Drops the two U1 terms of
/// the product rule in set_tau, whose per-axis products are what force tau to depth
/// 20-22 (U1 ~ x/r is non-smooth at the origin componentwise). Use it to test
/// whether the common refinement level is what drives the runaway.
inline bool xc_tau_no_u1() {
    static const bool flag = []() {
        const char* e = std::getenv("MAD_XC_TAU_NO_U1");
        return (e != nullptr) and (std::string(e) != "0");
    }();
    return flag;
}

/// MAD_XC_TAU_U1_POINTWISE: evaluate the two U1 terms of tau's product rule from
/// the analytic functor at the quadrature points of the (shallow) orbital tree,
/// instead of projecting U1 and |U1|^2 into MRA Functions and multiplying.
/// U1 ~ x/r is non-smooth at the origin componentwise, so U1_x as a Function needs
/// depth ~18 and every product with it inherits that depth; refine_to_common_level
/// then lifts every xc intermediate to it. Never forming the product is phase 1 of
/// XC-NEMO-MEMORY-INVESTIGATION.md 11.0. Unlike xc_tau_no_u1() this is not a
/// diagnostic -- it computes the same tau, only on a tree set by the orbitals.
inline bool xc_tau_u1_pointwise() {
    static const bool flag = []() {
        const char* e = std::getenv("MAD_XC_TAU_U1_POINTWISE");
        return (e != nullptr) and (std::string(e) != "0");
    }();
    return flag;
}

/// MAD_XC_NEMO_SIGMA: build grad(rho) and the sigma matrix from the split form

/// With rho = R^2 n and U1 = -grad(R)/R,
///   grad(rho_s) = R^2 (grad n_s - 2 U1 n_s) = 2 R^2 (G_s - U1 n_s),   G_s = 1/2 grad n_s
///   sigma_st    = 4 R^4 [ G_s.G_t - n_t U1.G_s - n_s U1.G_t + n_s n_t |U1|^2 ]
/// with the smooth pieces from MRA and U1 from its analytic functor, so nothing
/// involving U1 is differentiated numerically or projected. The default path
/// instead takes zeta = grad(log rho) -- a numerical derivative of a function with
/// a kink at every nucleus -- and forms sigma = rho^2 (zeta.zeta).
///
/// Measured on LiH/TPSS: the two agree to 4-5 digits outside 1e-3 bohr, but inside
/// that the default form errs by 30-70 %, with the Kato ratio |grad rho|/(2 Z rho)
/// swinging 0.73 to 1.31 where the split form holds 1.00-1.08. At Li that drives
/// z = tau_W/tau to 1.72 and fires the von Weizsaecker floor at the nucleus.
///
/// |U1|^2 must come from U1_dot_U1_functor, not from summing the squares of the
/// components: the vector's direction is smoothed over eprec and vanishes at r = 0,
/// which destroys the diagonal, while the scalar functor keeps (S'/S)^2 exact.
///
/// On by default *when the regularized pieces are present*, which is the pointwise
/// tau route; set to 0 to fall back and isolate the two changes.
inline bool xc_nemo_sigma() {
    static const bool flag = []() {
        const char* e = std::getenv("MAD_XC_NEMO_SIGMA");
        return (e == nullptr) or (std::string(e) != "0");
    }();
    return flag;
}

/// MAD_XC_EXPORT_DFDS: hand out 2 de/dsigma (same spin) as an extra pointwise
/// output of vxc, so a diagnostic can look at it directly. It is the coefficient
/// whose product with grad(rho) is the semilocal flux, hence the origin of the 1/r
/// in div(X); nothing in the potential depends on it being exported. Diagnostic
/// only -- it changes the size of vxc's result vector, so leave it off in
/// production and never index the result vector by a hardcoded tail offset.
inline bool xc_export_dfds() {
    static const bool flag = []() {
        const char* e = std::getenv("MAD_XC_EXPORT_DFDS");
        return (e != nullptr) and (std::string(e) != "0");
    }();
    return flag;
}

/// index of the 2 de/dsigma slot in vxc's result vector, or -1 when not exported

/// Shared by the producer (XCfunctional::vxc) and the consumer
/// (XCOperator::make_xc_potential) so the two cannot drift apart. The layout is
/// [0] local, [1..3] same-spin flux, [4..6] cross-spin flux if polarized,
/// [next] de/dtau if meta, [next] 2 de/dsigma if exported, then the A'' tail.
inline int xc_dfds_slot(const bool spin_polarized, const bool needs_tau,
                        const bool needs_sigma) {
    if (not (xc_export_dfds() and needs_sigma)) return -1;
    return (spin_polarized ? 7 : 4) + (needs_tau ? 1 : 0);
}

/// MAD_XC_PROBE: 1 = per-term tree/depth/norm census, 2 = also pointwise min/max
inline int xc_probe_level() {
    static const int level = []() {
        const char* e = std::getenv("MAD_XC_PROBE");
        return (e == nullptr) ? 0 : std::atoi(e);
    }();
    return level;
}

/// Compute the spin-restricted LDA potential using unaryop (only for the initial guess)
struct xc_lda_potential {
    xc_lda_potential() {}

    void operator()(const Key<3> & key, Tensor<double>& t) const
    {
        int x_rks_s__(const double *r__, double *f, double * dfdra);
        int c_rks_vwn5__(const double *r__, double *f, double * dfdra);
        double* rho = t.ptr();
        for (int i=0; i<t.size(); i++) {
            double r = std::max(rho[i],1e-12);
            double q, dq1, dq2;
            x_rks_s__(&r, &q, &dq1);
            c_rks_vwn5__(&r, &q, &dq2);
            rho[i] = dq1 + dq2;
        }
    }
};

/// Simplified interface to XC functionals
class XCfunctional {
public:

    /// The ordering of the intermediates is fixed, but the code can handle
    /// non-initialized functions, so if e.g. no GGA is requested, all the
    /// corresponding vector components may be left empty.
    ///
    /// Note the additional quantities \f$ \zeta \f$ and \f$ \chi \f$, which are defined as
    /// \f[
    /// \rho = \exp(\zeta)
    /// \f]
    /// and thus the derivative of rho is given by
    /// \f[
    /// \nabla_x\rho = \exp(\zeta)\nabla_x\zeta = \rho \nabla_x\zeta
    /// \f]
    /// The reduced gradients \sigma may then be expressed as
    /// \f[
    ///   \sigma = |\nabla\rho|^2 = |\rho|^2 |\nabla\zeta|^2 = |\rho|^2 \chi
    /// \f]
    enum xc_arg {
        enum_rhoa=0,            ///< alpha density \f$ \rho_\alpha \f$
        enum_rhob=1,            ///< beta density \f$ \rho_\beta \f$
        enum_rho_pt=2,          ///< perturbed density (CPHF, TDKS) \f$ \rho_{pt} \f$
        enum_taua=3,            ///< alpha kinetic energy density \f$ \tau_\alpha = \frac{1}{2}\sum_i|\nabla\psi_{i\alpha}|^2 \f$
        enum_taub=4,            ///< beta kinetic energy density \f$ \tau_\beta \f$
        enum_saa=10,            ///< \f$ \sigma_{aa} = \nabla \rho_{\alpha}.\nabla \rho_{\alpha} \f$
        enum_sab=11,            ///< \f$ \sigma_{ab} = \nabla \rho_{\alpha}.\nabla \rho_{\beta} \f$
        enum_sbb=12,            ///< \f$ \sigma_{bb} = \nabla \rho_{\beta}.\nabla \rho_{\beta} \f$
        enum_sigtot=13,         ///< \f$ \sigma = \nabla \rho.\nabla \rho \f$
        enum_sigma_pta_div_rho=14,      ///< \f$ \zeta_{\alpha}.\nabla\rho_{pt} \f$
        enum_sigma_ptb_div_rho=15,      ///< \f$ \zeta_{\beta}.\nabla\rho_{pt} \f$
        enum_zetaa_x=16,        ///< \f$ \zeta_{a,x}=\partial/{\partial x} \ln(\rho_a)  \f$
        enum_zetaa_y=17,        ///< \f$ \zeta_{a,y}=\partial/{\partial y} \ln(\rho_a)  \f$
        enum_zetaa_z=18,        ///< \f$ \zeta_{a,z}=\partial/{\partial z} \ln(\rho_a)  \f$
        enum_zetab_x=19,        ///< \f$ \zeta_{b,x} = \partial/{\partial x} \ln(\rho_b)  \f$
        enum_zetab_y=20,        ///< \f$ \zeta_{b,y} = \partial/{\partial y} \ln(\rho_b)  \f$
        enum_zetab_z=21,        ///< \f$ \zeta_{b,z} = \partial/{\partial z} \ln(\rho_b)  \f$
        // Slots 22-24 held chi_st = zeta_s.zeta_t as three separately represented
        // functions. They are gone: a projected product is not pointwise consistent
        // with the zeta components it is built from, so the sigma matrix handed to
        // libxc was not the Gram matrix of the density gradients -- chi_aa, a sum of
        // squares, could come out negative near the nuclear cusp, and the total sigma
        // followed it. make_libxc_args contracts zeta pointwise instead. Left as a
        // hole rather than reused, so the surviving indices keep their meaning.
        enum_ddens_ptx=25,      ///< \f$ \nabla\rho_{pt}\f$
        enum_ddens_pty=26,      ///< \f$ \nabla\rho_{pt}\f$
        enum_ddens_ptz=27,      ///< \f$ \nabla\rho_{pt}\f$

        // ---- the nemo meta-gga decomposition, kept apart on purpose ----------
        //
        // With psi = R F the product rule gives
        //   |grad psi|^2 = R^2 ( |grad F|^2 - 2 F U1.grad F + |U1|^2 F^2 ).
        // The first group below is everything in that expression which is SMOOTH:
        // F is cusp-free by construction, so |grad F|^2, n and G are shallow (depth
        // 8-9 on LiH) and belong in MRA. R^2 is smooth too.
        //
        // U1 = -grad(R)/R is not. U1_x ~ x/r is non-smooth at every nucleus
        // componentwise, its direction smoothed only over eprec, and as a Function
        // it costs depth ~20. Carrying it in MRA is bad twice over: the depth taxes
        // every intermediate through refine_to_common_level, and any *product* with
        // it has to be projected onto a fixed tree, which is where the oscillations
        // come from. So it is never a Function. The second group is filled by the
        // op from the analytic functor at the quadrature points of the box it is
        // already working on -- exactly as make_libxc_args contracts zeta pointwise
        // rather than carrying chi. See nemo_u1_functors.
        enum_nemo_R2=28,        ///< \f$ R^2 \f$, the ncf squared          [MRA, smooth]
        enum_gradfa=29,         ///< \f$ \sum_i w_i|\nabla F_{i\alpha}|^2 \f$  [MRA, smooth]
        enum_gradfb=30,         ///< beta counterpart                        [MRA, smooth]
        enum_na=31,             ///< \f$ n_\alpha=\sum_i w_iF_{i\alpha}^2 \f$   [MRA, smooth]
        enum_nb=32,             ///< beta counterpart                        [MRA, smooth]
        enum_Ga_x=33,           ///< \f$ G_{\alpha,x}=\sum_i w_iF_i\partial_xF_i \f$ [MRA, smooth]
        enum_Ga_y=34,           ///< \f$ G_{\alpha,y} \f$                    [MRA, smooth]
        enum_Ga_z=35,           ///< \f$ G_{\alpha,z} \f$                    [MRA, smooth]
        enum_Gb_x=36,           ///< beta counterpart                        [MRA, smooth]
        enum_Gb_y=37,           ///< beta counterpart                        [MRA, smooth]
        enum_Gb_z=38,           ///< beta counterpart                        [MRA, smooth]

        enum_u1_x=39,           ///< \f$ U_{1,x} \f$        [functor, never MRA]
        enum_u1_y=40,           ///< \f$ U_{1,y} \f$        [functor, never MRA]
        enum_u1_z=41,           ///< \f$ U_{1,z} \f$        [functor, never MRA]
        enum_u1sq=42            ///< \f$ |\mathbf U_1|^2 \f$ [functor, never MRA]
    };
    const static int number_xc_args=43;     ///< max number of intermediates

    /// return the munging threshold for the density
    double get_rhotol() const {return rhotol;}

    /// return the binary munging threshold for the final result in the GGA potential/kernel

    /// the GGA potential will be munged based on the smallness of the original
    /// density, which we call binary munging
    double get_ggatol() const {return ggatol;}

    /// return the floor for the kinetic energy density

    /// meta-gga functionals build the iso-orbital indicators alpha and z with tau
    /// in the denominator, so tau needs a floor well above libxc's own 1e-20
    /// default for a real-space code, where tau is a numerical derivative
    double get_tautol() const {return tautol;}

    /// C^2 ramp from 0 (rho <= rhotol) to 1 (rho >= ggatol), smooth in log(rho)

    /// The semilocal flux 2 de/dsigma grad(rho) is switched off by munge(), which
    /// sets the density hard to zero below rhotol -- a jump discontinuity in the
    /// flux at the rho = rhotol iso-surface. That jump is negligible while the
    /// flux carries a factor rho, but the factored form divides it out and
    /// de/dsigma grows like rho^{-1/3} in the tail, so the jump becomes O(10^2)
    /// and its divergence is unrepresentable. This ramp replaces it with a C^2
    /// switch over the three decades rhotol..ggatol, which is where the tail
    /// screening of the meta-gga term already lives.
    double gga_ramp(const double rho) const {
        if (rho >= ggatol) return 1.0;
        if (rho <= rhotol) return 0.0;
        const double t = std::log(rho/rhotol)/std::log(ggatol/rhotol);
        return t*t*t*(t*(6.0*t - 15.0) + 10.0);
    }

protected:

    bool spin_polarized;        ///< True if the functional is spin polarized
    double hf_coeff;            ///< Factor multiplying HF exchange (+1.0 gives HF)
    double rhomin, rhotol;      ///< See initialize and munge*
    double ggatol;              ///< See initialize and munge*
    double tautol;              ///< floor for the kinetic energy density, see initialize

#ifdef MADNESS_HAS_LIBXC
    std::vector< std::pair<xc_func_type*,double> > funcs;
#endif

    /// convert the raw density (gradient) data to be used by the xc operators

    /// Involves mainly munging of the densities and multiplying with 2
    /// if the calculation is spin-restricted.
    /// Response densities and density gradients are munged based on the
    /// value of the ground state density, since they may become negative
    /// and may also be much more diffuse.
    /// dimensions of the output tensors are for spin-restricted and unrestricted
    /// (with np the number of grid points in the box):
    /// rho(np) or rho(2*np)
    /// sigma(np) sigma(3*np)
    /// rho_pt(np)
    /// sigma_pt(2*np)
    /// @param[in]  t       input density (gradients)
    /// @param[out] rho     ground state (spin) density, properly munged
    /// @param[out] sigma   ground state (spin) density gradients, properly munged
    /// @param[out] tau     ground state (spin) kinetic energy density, properly munged
    /// @param[out] rho_pt  response density, properly munged (no spin)
    /// @param[out] sigma_pt  response (spin) density gradients, properly munged
    /// @param[out] drho    density derivative, constructed from rho and zeta
    /// @param[out] drho_pt response density derivative directly from xc_args
    /// @param[in]  need_response   flag if rho_pt and sigma_pt need to be calculated
    void make_libxc_args(const std::vector< madness::Tensor<double> >& t,
                         madness::Tensor<double>& rho,
                         madness::Tensor<double>& sigma,
                         madness::Tensor<double>& tau,
                         madness::Tensor<double>& rho_pt,
                         madness::Tensor<double>& sigma_pt,
                         std::vector<madness::Tensor<double> >& drho,
                         std::vector<madness::Tensor<double> >& drho_pt,
                         const bool need_response) const;

    /// the number of xc kernel derivatives (lda: 0, gga: 1, etc)
    int nderiv;


    /// Smoothly switches between constant (x<xmin) and linear function (x>xmax)

    /// \f[
    /// f(x,x_{\mathrm{min}},x_{\mathrm{max}}) = \left\{
    ///   \begin{array}{ll}
    ///     x_{\mathrm{min}}                       & x < x_{\mathrm{min}} \\
    ///     p(x,x_{\mathrm{min}},x_{\mathrm{max}}) & x_{\mathrm{min}} \leq x_{\mathrm{max}} \\
    ///     x                                      & x_{\mathrm{max}} < x
    ///   \end{array}
    /// \right.
    /// \f]
    /// where \f$p(x)\f$ is the unique quintic polynomial that
    /// satisfies \f$p(x_{min})=x_{min}\f$, \f$p(x_{max})=x_{max}\f$,
    /// \f$dp(x_{max})/dx=1\f$, and
    /// \f$dp(x_{min})/dx=d^2p(x_{min})/dx^2=d^2p(x_{max})/dx^2=0\f$.
    static void polyn(const double x, double& p, double& dpdx) {
        // All of the static const stuff is evaluated at compile time

        static const double xmin = 1.e-6; // <<<< MINIMUM VALUE OF DENSITY
        static const double xmax = 5.e-5;  // <<<< DENSITY SMOOTHLY MODIFIED BELOW THIS VALUE

        static const double xmax2 = xmax*xmax;
        static const double xmax3 = xmax2*xmax;
        static const double xmin2 = xmin*xmin;
        static const double xmin3 = xmin2*xmin;
        static const double r = 1.0/((xmax-xmin)*(-xmin3+(3.0*xmin2+(-3.0*xmin+xmax)*xmax)*xmax));
        static const double a0 = xmax3*xmin*(xmax-4.0*xmin)*r;
        static const double a = xmin2*(xmin2+(-4.0*xmin+18.0*xmax)*xmax)*r;
        static const double b = -6.0*xmin*xmax*(3.0*xmax+2.0*xmin)*r;
        static const double c = (4.0*xmin2+(20.0*xmin+6.0*xmax)*xmax)*r;
        static const double d = -(8.0*xmax+7.0*xmin)*r;
        static const double e = 3.0*r;

        if (x > xmax) {
            p = x;
            dpdx = 1.0;
        }
        else if (x < xmin) {
            p = xmin;
            dpdx = 0.0;
        }
        else {
            p = a0+(a+(b+(c+(d+e*x)*x)*x)*x)*x;
            dpdx = a+(2.0*b+(3.0*c+(4.0*d+5.0*e*x)*x)*x)*x;
        }
    }
public:
    static double munge_old(double rho) {
        double p, dpdx;
        polyn(rho, p, dpdx);
        return p;
    }

private:

    /// simple munging for the density only (LDA)
    double munge(double rho) const {
    	if (rho <= rhotol) rho=rhomin;
        return rho;
    }

    /// munge rho if refrho is small

    /// special case for perturbed densities, which might be negative and diffuse.
    /// Munge rho (e.g. the perturbed density) if the reference density refrho
    /// e.g. the ground state density is small. Only where the reference density
    /// is large enough DFT is numerically well-defined.
    /// @param[in]  rho     number to be munged
    /// @param[in]  refrho  reference value for munging
    /// @param[in]  thresh  threshold for munging
    double binary_munge(double rho, double refrho, const double thresh) const {
        if (refrho<thresh) rho=rhomin;
        return rho;
    }

public:
    /// Default constructor is required
    XCfunctional();

    /// Initialize the object from the user input data

    /// @param[in] input_line User input line (without beginning XC keyword)
    /// @param[in] polarized Boolean flag indicating if the calculation is spin-polarized
    void initialize(const std::string& input_line, bool polarized, World& world,
            const bool verbose=false);

    /// Destructor
    ~XCfunctional();

    /// Returns true if the potential is lda
    bool is_lda() const;

    /// Returns true if the potential is gga (needs first derivatives)
    bool is_gga() const;

    /// Returns true if the potential is meta gga (needs the kinetic energy density)
    bool is_meta() const;

    /// Returns true if the functional needs the reduced density gradients sigma

    /// True for gga AND meta-gga -- a meta-gga needs the density gradients as
    /// well as the kinetic energy density. Use this, not is_gga(), to decide
    /// whether the gradient intermediates have to be computed.
    bool needs_sigma() const;

    /// Returns true if the functional needs the kinetic energy density tau
    bool needs_tau() const;

    /// Returns true if there is a DFT functional (false probably means Hatree-Fock exchange only)
    bool is_dft() const;

    /// Returns true if the functional is spin_polarized
    bool is_spin_polarized() const
    {
        return spin_polarized;
    }

    /// Returns true if the second derivative of the functional is available (not yet supported)
    bool has_fxc() const;

    /// Returns true if the third derivative of the functional is available (not yet supported)
    bool has_kxc() const;

    /// Returns the value of the hf exact exchange coefficient
    double hf_exchange_coefficient() const
    {
        return hf_coeff;
    }

    /// Computes the energy functional at given points

    /// This uses the convention that the total energy is
    /// \f$ E[\rho] = \int \epsilon[\rho(x)] dx\f$
    /// Any HF exchange contribution must be separately computed. Items in the
    /// vector argument \c t are interpreted similarly to the xc_arg enum.
    /// @param[in] t The input densities and derivatives as required by the functional
    /// @return The exchange-correlation energy functional
    madness::Tensor<double> exc(const std::vector< madness::Tensor<double> >& t) const;

    /// Computes components of the potential (derivative of the energy functional) at np points

    /// Any HF exchange contribution must be separately computed. Items in the
    /// vector argument \c t are interpreted similarly to the xc_arg enum.
    ///
    /// We define \f$ \sigma_{\mu \nu} = \nabla \rho_{\mu} . \nabla \rho_{\nu} \f$
    /// with \f$ \mu, \nu = \alpha\f$ or \f$ \beta \f$.
    ///
    /// For unpolarized GGA, matrix elements of the potential are
    /// \f[
    ///   < \phi | \hat V | \psi > = \int \left( \frac{\partial \epsilon}{\partial \rho} \phi \psi
    ///                  +  \left( 2 \frac{\partial \epsilon}{\partial \sigma} \right)
    ///                  \nabla \rho \cdot \nabla \left( \phi \psi \right) \right) dx
    /// \f]
    ///
    /// For polarized GGA, matrix elements of the potential are
    /// \f[
    ///   < \phi_{\alpha} | \hat V | \psi_{\alpha} > = \int \left( \frac{\partial \epsilon}{\partial \rho_{\alpha}} \phi \psi
    ///            +  \left( 2 \frac{\partial \epsilon}{\partial \sigma_{\alpha \alpha}} \nabla \rho_{\alpha}
    ///            + \frac{\partial \epsilon}{\partial \sigma_{\alpha \beta}} \nabla \rho_{\beta}  \right) . \nabla \left( \phi \psi \right) \right) dx
    /// \f]
    ///
    /// Integrating the above by parts and assuming free-space or periodic boundary conditions
    /// we obtain that the local multiplicative form of the GGA potential is
    /// \f[
    ///    V_{\alpha} =  \frac{\partial \epsilon}{\partial \rho_{\alpha}}
    ///                  - \left(\nabla . \left(2 \frac{\partial \epsilon}{\partial \sigma_{\alpha \alpha}} \nabla \rho_{\alpha}
    ///                  + \frac{\partial \epsilon}{\partial \sigma_{\alpha \beta}} \nabla \rho_{\beta}  \right)  \right)
    /// \f]
    ///
    /// Return the following quantities for RHF: (see Yanai2005, Eq. (12))
    /// \f{eqnarray*}{
    ///     \mbox{result[0]}    &:& \qquad \frac{\partial \epsilon}{\partial \rho} \\
    ///     \mbox{result[1-3]}  &:& \qquad 2 \rho \frac{\partial \epsilon}{\partial \sigma} \nabla\rho
    /// \f}
    /// and for UHF same-spin and other-spin quantities
    /// \f{eqnarray*}{
    ///     \mbox{result[0]}    &:& \qquad \frac{\partial \epsilon}{\partial \rho_{\alpha}} \\
    ///     \mbox{result[1-3]}  &:& \qquad \rho_\alpha \frac{\partial \epsilon}{\partial \sigma_{\alpha \alpha}} \nabla\rho_\alpha\\
    ///     \mbox{result[4-6]}  &:& \qquad \rho_\alpha \frac{\partial \epsilon}{\partial \sigma_{\alpha \beta}} \nabla\rho_\beta
    /// \f}
    /// @param[in] t The input densities and derivatives as required by the functional
    /// @param[in] ispin Specifies which component of the potential is to be computed as described above
    /// @return the requested quantity, based on ispin (0: same spin, 1: other spin)
    std::vector<madness::Tensor<double> > vxc(const std::vector< madness::Tensor<double> >& t,
            const int ispin) const;


    /// compute the second derivative of the XC energy wrt the density and apply

    /// Return the following quantities (RHF only) (see Yanai2005, Eq. (13))
    /// \f{eqnarray*}{
    ///     \mbox{result[0]}    &:& \qquad \frac{\partial^2 \epsilon}{\partial \rho^2} \rho_\mathrm{pt}
    ///                             + 2.0 * \frac{\partial^2 \epsilon}{\partial \rho\partial\sigma}\sigma_\mathrm{pt}\\
    ///     \mbox{result[1-3]}  &:& \qquad 2.0 * \frac{\partial\epsilon}{\partial\sigma}\nabla\rho_\mathrm{pt}
    ///                             + 2.0 * \frac{\partial^2\epsilon}{\partial\rho\partial\sigma} \rho_\mathrm{pt}\nabla\rho
    ///                             + 4.0 * \frac{\partial^2\epsilon}{\partial^2\sigma} \sigma_\mathrm{pt}\nabla\rho
    /// \f}
    /// @param[in]  t   The input densities and derivatives as required by the functional,
    ///                 as in the xc_arg enum
    /// @param[in]  ispin not referenced since only RHF is implemented, always 0
    /// @return a vector of Functions containing the contributions to the kernel apply
    std::vector<madness::Tensor<double> > fxc_apply(
            const std::vector< madness::Tensor<double> >& t, const int ispin) const;


    /// Crude function to plot the energy and potential functionals
    void plot() const {
        long npt = 1001;
        double lo=1e-6, hi=1e+1, s=std::pow(hi/lo, 1.0/(npt-1));

        madness::Tensor<double> rho(npt);
        for (int i=0; i<npt; i++) {
            rho[i] = lo;
            lo *= s;
        }
        std::vector< madness::Tensor<double> > t(13);
        t[enum_rhoa]=(rho);
        if (is_spin_polarized()) t[enum_rhob]=(rho);
//        if (is_gga()) t[enum_saa]=madness::Tensor<double>(npt); // sigma_aa=0
        if (needs_sigma()) t[enum_saa]=0.5*rho; // sigma_aa=0
        if (needs_tau()) t[enum_taua]=0.5*rho;
        madness::Tensor<double> f  = exc(t); //pending UGHHHHH
        std::vector<madness::Tensor<double> > va = vxc(t,0);
        for (long i=0; i<npt; i++) {
            printf("%.3e %.3e %.3e\n", rho[i], f[i], va[0][i]);
        }
    }
};

/// Class to compute the energy functional
/// the cuspy half of the nemo tau decomposition, supplied pointwise

/// Holds the four analytic ncf quantities that must never become MRA functions --
/// \f$ U_{1,x}, U_{1,y}, U_{1,z}, |\mathbf U_1|^2 \f$ -- and writes their values
/// into the argument vector at the quadrature points of whatever box the caller is
/// operating on. Empty unless a nuclear correlation factor is in play, in which case
/// active() is true and the four enum_u1* slots are filled.
///
/// The point is that nothing is projected. A product of U1 with anything, formed as
/// a Function, has to be represented on some tree; on a tree too coarse for U1's
/// eprec-scale structure the polynomial fit rings across the whole box. Evaluating
/// U1 here instead means its values go straight into the functional's pointwise
/// arithmetic and only the *potential* is ever projected -- which the existing
/// machinery already does, and already has to.
///
/// \f$ |\mathbf U_1|^2 \f$ comes from its own functor rather than from summing the
/// squares of the three components: that functor treats its diagonal specially,
/// because smoothed_unitvec has norm < 1 inside eprec while the exact diagonal is
/// \f$ (S'/S)^2 \f$.
struct nemo_u1_functors {
    typedef madness::FunctionFunctorInterface<double,3> functorT;

    nemo_u1_functors() : cdata(madness::FunctionCommonData<double,3>::get(
                                       madness::FunctionDefaults<3>::get_k())) {}

    /// @param[in] u1  x, y, z components of U1 followed by |U1|^2 -- four functors
    explicit nemo_u1_functors(const std::vector<std::shared_ptr<functorT> >& u1)
            : f(u1), cdata(madness::FunctionCommonData<double,3>::get(
                                   madness::FunctionDefaults<3>::get_k())) {
        MADNESS_CHECK_THROW(f.empty() or f.size()==4,
                            "nemo_u1_functors wants U1_{x,y,z} and |U1|^2, in that order");
    }

    bool active() const {return f.size()==4;}

    /// write U1 and |U1|^2 at this box's quadrature points into t[enum_u1*]
    void append(const madness::Key<3>& key,
                std::vector<madness::Tensor<double> >& t) const {
        if (not active()) return;
        if (long(t.size()) < XCfunctional::number_xc_args)
            t.resize(XCfunctional::number_xc_args);

        const madness::Tensor<double>& qx = cdata.quad_x;
        const long npt = qx.dim(0);
        // cdata was captured at construction from FunctionDefaults; if the functions
        // actually carry a different k the quadrature points below are the wrong
        // ones and every U1 value lands at the wrong place. Silent, and it would
        // look like a physics error, so check rather than trust.
        if (t[XCfunctional::enum_rhoa].size())
            MADNESS_CHECK_THROW(t[XCfunctional::enum_rhoa].dim(0) == npt,
                                "nemo_u1_functors: quadrature order does not match "
                                "the xc arguments -- k changed after construction");
        const double h = std::pow(0.5, double(key.level()));
        const madness::Tensor<double>& cell = madness::FunctionDefaults<3>::get_cell();
        const madness::Tensor<double>& cw = madness::FunctionDefaults<3>::get_cell_width();

        const long dims[3] = {npt, npt, npt};
        madness::Tensor<double> v[4];
        double* p[4];
        for (int q = 0; q < 4; ++q) {
            v[q] = madness::Tensor<double>(3L, dims);
            p[q] = v[q].ptr();
        }

        // the same box-to-user-coordinate construction fcube() uses, written out so
        // this header needs no mraimpl.h
        long idx = 0;
        madness::Vector<double,3> c;
        for (long i = 0; i < npt; ++i) {
            c[0] = cell(0,0) + h*cw[0]*(key.translation()[0] + qx(i));
            for (long j = 0; j < npt; ++j) {
                c[1] = cell(1,0) + h*cw[1]*(key.translation()[1] + qx(j));
                for (long k = 0; k < npt; ++k, ++idx) {
                    c[2] = cell(2,0) + h*cw[2]*(key.translation()[2] + qx(k));
                    for (int q = 0; q < 4; ++q) p[q][idx] = (*f[q])(c);
                }
            }
        }
        t[XCfunctional::enum_u1_x]  = v[0];
        t[XCfunctional::enum_u1_y]  = v[1];
        t[XCfunctional::enum_u1_z]  = v[2];
        t[XCfunctional::enum_u1sq]  = v[3];
    }

    std::vector<std::shared_ptr<functorT> > f;
    madness::FunctionCommonData<double,3> cdata;
};


struct xc_functional {
    const XCfunctional* xc;
    nemo_u1_functors u1;      ///< empty without a nuclear correlation factor

    xc_functional(const XCfunctional& xc) : xc(&xc) {}
    xc_functional(const XCfunctional& xc, const nemo_u1_functors& u1) : xc(&xc), u1(u1) {}

    madness::Tensor<double> operator()(const madness::Key<3> & key,
            const std::vector< madness::Tensor<double> >& t) const {
        MADNESS_ASSERT(xc);
        if (not u1.active()) return xc->exc(t);
        std::vector<madness::Tensor<double> > tt(t);   // Tensor copy is shallow
        u1.append(key, tt);
        return xc->exc(tt);
    }
};


/// Class to compute terms of the potential
struct xc_potential {
    const XCfunctional* xc;
    const int ispin;
    nemo_u1_functors u1;      ///< empty without a nuclear correlation factor

    xc_potential(const XCfunctional& xc, int ispin) : xc(&xc), ispin(ispin)
    {}
    xc_potential(const XCfunctional& xc, int ispin, const nemo_u1_functors& u1)
            : xc(&xc), ispin(ispin), u1(u1)
    {}

    std::size_t get_result_size() const {
        // local terms, same spin
        if (xc->is_lda()) return 1;
        // local term + 3x semilocal terms (x,y,z), same spin
        std::size_t result_size=4;
        // 3x semilocal terms (x,y,z) for the opposite spin
        if (xc->is_spin_polarized()) result_size+=3;
        // de/dtau, same spin -- the non-multiplicative meta-gga term
        if (xc->needs_tau()) result_size+=1;
        // MAD_XC_EXPORT_DFDS: 2 de/dsigma, diagnostic. Sits here, before the A''
        // block, so that block keeps its nr-4 .. nr-1 tail indexing.
        if (xc_export_dfds() and xc->needs_sigma()) result_size+=1;
        // algorithm A'': u = 2 de/dsigma and the munged grad(rho), appended last.
        // grad(rho) has to come from here rather than be rebuilt as rho*zeta at
        // the MRA level: the flux uses the *munged* density, which vanishes below
        // rhotol, and an unmunged rho*zeta carries zeta's tail garbage all the way
        // to the box wall, so its divergence picks up a surface term.
        if (xc_onelevel_div() and xc->needs_sigma()) result_size+=4;
        return result_size;
    }

    std::vector<madness::Tensor<double> > operator()(const madness::Key<3> & key,
            const std::vector< madness::Tensor<double> >& t) const {
        MADNESS_ASSERT(xc);
        if (not u1.active()) return xc->vxc(t, ispin);
        // U1 is cuspy, so it arrives here as values rather than as a Function --
        // see nemo_u1_functors. Nothing involving it is ever projected.
        std::vector<madness::Tensor<double> > tt(t);   // Tensor copy is shallow
        u1.append(key, tt);
        return xc->vxc(tt, ispin);
    }
};


/// Class to compute terms of the kernel
struct xc_kernel_apply {
    const XCfunctional* xc;
    const int ispin;
    const FunctionCommonData<double,3>& cdata;

    xc_kernel_apply(const XCfunctional& xc, int ispin) : xc(&xc), ispin(ispin),
          cdata(FunctionCommonData<double,3>::get(FunctionDefaults<3>::get_k())) {
        MADNESS_ASSERT(ispin==0);   // closed shell only!
    }

    std::size_t get_result_size() const {
        // all spin-restricted
        if (xc->is_gga()) return 4; // local terms,  3x semilocal terms (x,y,z)
        return 1;   // local terms only
    }

    std::vector<madness::Tensor<double> > operator()(const madness::Key<3> & key,
            const std::vector< madness::Tensor<double> >& t) const {
        MADNESS_ASSERT(xc);
        std::vector<madness::Tensor<double> > r = xc->fxc_apply(t, ispin);
        return r;
    }
};

}

MADNESS_PRAGMA_CLANG(diagnostic pop)
MADNESS_PRAGMA_GCC(diagnostic pop)

#endif
