#ifndef MADNESS_CHEM_XCFUNCTIONAL_H__INCLUDED
#define MADNESS_CHEM_XCFUNCTIONAL_H__INCLUDED

#include <madness/madness_config.h>

/// \file moldft/xcfunctional.h
/// \brief Defines interface for DFT XC functionals
/// \ingroup chemistry

#include <madness/tensor/tensor.h>
#include <vector>
#include <algorithm>
#include <utility>
#include <madness/mra/key.h>
#include <madness/world/MADworld.h>
#include <madness/mra/function_common_data.h>

#ifdef MADNESS_HAS_LIBXC
#include <xc.h>
#endif

namespace madness {
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
    /// \f$
    /// \rho = exp(\zeta)
    /// \f$
    /// and thus the derivative of rho is given by
    /// \f$
    /// \nabla_x\rho = exp(\zeta)\nabla_x\zeta = \rho \nabla_x\zeta
    /// \f$
    /// The reduced gradients \sigma may then be expressed as
    /// \f$
    ///   \sigma = |\nabla\rho|^2 = |\rho|^2 |\nabla\zeta|^2 = |\rho|^2 \chi
    /// \f$
    enum xc_arg {
        enum_rhoa=0,            ///< alpha density \f$ \rho_\alpha \f$
        enum_rhob=1,            ///< beta density \f$ \rho_\beta \f$
        enum_rho_pt=2,          ///< perturbed density (CPHF, TDKS) \f$ \rho_{pt} \f$
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
        enum_chi_aa=22,         ///< \f$ \chi_{aa} = \nabla \zeta_{\alpha}.\nabla \zeta_{\alpha} \f$
        enum_chi_ab=23,         ///< \f$ \chi_{ab} = \nabla \zeta_{\alpha}.\nabla \zeta_{\beta} \f$
        enum_chi_bb=24          ///< \f$ \chi_{bb} = \nabla \zeta_{\beta}.\nabla \zeta_{\beta} \f$
    };
    const static int number_xc_args=25;     ///< max number of intermediates


    /// which contribution is requested from the XCfunctional

    /// Three types of contributions are required:
    /// those that are second derivatives of the xc kernel and that are local
    /// those that are first/second derivatives of the kernel and of which the
    /// derivative will be taken later on (-> semi-local).
    /// See vxc() and fxc_apply() for more details.
    enum xc_contribution {
        potential_rho,              ///< potential df/drho
        potential_same_spin,        ///< potential df/dsigma_aa
        potential_mixed_spin,       ///< potential df/dsigma_ab
        kernel_second_local,        ///< kernel (2nd derivative)
        kernel_second_semilocal,    ///< kernel (2nd derivative)
        kernel_first_semilocal      ///< kernel (2nd derivative)
    };

    /// return the munging threshold for the density
    double get_rhotol() const {return rhotol;}

    /// return the binary munging threshold for the final result in the GGA potential/kernel

    /// the GGA potential will be munged based on the smallness of the original
    /// density, which we call binary munging
    double get_ggatol() const {return ggatol;}

protected:

    bool spin_polarized;        ///< True if the functional is spin polarized
    double hf_coeff;            ///< Factor multiplying HF exchange (+1.0 gives HF)
    double rhomin, rhotol;      ///< See initialize and munge*
    double ggatol;              ///< See initialize and munge*

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
    /// @param[out] rho_pt  response density, properly munged (no spin)
    /// @param[out] sigma_pt  response (spin) density gradients, properly munged
    /// @param[in]  need_response   flag if rho_pt and sigma_pt need to be calculated
    void make_libxc_args(const std::vector< madness::Tensor<double> >& t,
                         madness::Tensor<double>& rho,
                         madness::Tensor<double>& sigma,
                         madness::Tensor<double>& rho_pt,
                         madness::Tensor<double>& sigma_pt,
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
    double binary_munge(double rho, double refrho) const {
        if (refrho<rhotol) rho=rhomin;
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

    /// Returns true if the potential is meta gga (needs second derivatives ... not yet supported)
    bool is_meta() const;

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
    /// \f$
    ///   < \phi | \hat V | \psi > = \int \left( \frac{\partial \epsilon}{\partial \rho} \phi \psi
    ///                  +  \left( 2 \frac{\partial \epsilon}{\partial \sigma} \right)
    ///                  \nabla \rho \cdot \nabla \left( \phi \psi \right) \right) dx
    /// \f$
    ///
    /// For polarized GGA, matrix elements of the potential are
    /// \f$
    ///   < \phi_{\alpha} | \hat V | \psi_{\alpha} > = \int \left( \frac{\partial \epsilon}{\partial \rho_{\alpha}} \phi \psi
    ///            +  \left( 2 \frac{\partial \epsilon}{\partial \sigma_{\alpha \alpha}} \nabla \rho_{\alpha}
    ///            + \frac{\partial \epsilon}{\partial \sigma_{\alpha \beta}} \nabla \rho_{\beta}  \right) . \nabla \left( \phi \psi \right) \right) dx
    /// \f$
    ///
    /// Integrating the above by parts and assuming free-space or periodic boundary conditions
    /// we obtain that the local multiplicative form of the GGA potential is
    /// \f$
    ///    V_{\alpha} =  \frac{\partial \epsilon}{\partial \rho_{\alpha}}
    ///                  - \left(\nabla . \left(2 \frac{\partial \epsilon}{\partial \sigma_{\alpha \alpha}} \nabla \rho_{\alpha}
    ///                  + \frac{\partial \epsilon}{\partial \sigma_{\alpha \beta}} \nabla \rho_{\beta}  \right)  \right)
    /// \f$
    ///
    /// Until we get a madness::Function operation that can produce
    /// multiple results we need to compute components of the
    /// functional and potential separately.
    ///
    /// Return the following quantities depending on the xc_contribution request.
    /// note the additional \f$\rho\f$ in some terms.
    /// \f{align}{
    ///     \text{potentialrho, RHF}      : &\qquad \frac{\partial \epsilon}{\partial \rho} \\
    ///     \text{potentialrho, UHF}      : &\qquad \frac{\partial \epsilon}{\partial \rho_{\alpha}} \\
    ///     \text{potentialsamespin, RHF} : &\qquad \rho \frac{\partial \epsilon}{\partial \sigma} \\
    ///     \text{potentialsamespin, UHF} : &\qquad \rho_\alpha \frac{\partial \epsilon}{\partial \sigma_{\alpha \alpha}}\\
    ///     \text{potentialotherspin, UHF}: &\qquad \rho_\alpha \frac{\partial \epsilon}{\partial \sigma_{\alpha \beta}}
    /// \f}
    /// @param[in] t The input densities and derivatives as required by the functional
    /// @param[in] ispin Specifies which component of the potential is to be computed as described above
    /// @return the requested quantity, based on ispin and xc_contrib
    madness::Tensor<double> vxc(const std::vector< madness::Tensor<double> >& t,
            const int ispin, const xc_contribution xc_contrib) const;


    /// compute the second derivative of the XC energy wrt the density and apply

    /// Return the following quantities depending on the xc_contribution request.
    /// \f{align}{
    ///     \text{kernelsecondlocal}      : &\qquad \frac{\partial^2 \epsilon}{\partial \rho^2} \rho_\mathrm{pt}
    ///                                           + 2.0 * \frac{\partial^2 \epsilon}{\partial \rho\partial\sigma}\sigma_\mathrm{pt}\\
    ///     \text{kernelfirstsemilocal}   : &\qquad 2.0 * \frac{\partial\epsilon}{\partial\sigma}\\
    ///     \text{kernelsecondsemilocal}  : &\qquad 2.0 * \frac{\partial^2\epsilon}{\partial\rho\partial\sigma} \rho_\mathrm{pt}
    ///                                           + 4.0 * \frac{\partial^2\epsilon}{\partial^2\sigma} \sigma_\mathrm{pt}
    /// \f}
    /// semilocal means that the functional derivative has to be multiplied with
    /// the gradient of the (perturbed) density and then the div is taken.
    /// @param[in]  t   The input densities and derivatives as required by the functional,
    ///                 as in the xc_arg enum
    /// @param[in]  ispin not referenced since only RHF is implemented, always 0
    /// @param[in]  xc_contrib  the contribution of the functional derivatives which is requested
    /// @return The component specified by the \c xc_contrib parameter
    madness::Tensor<double> fxc_apply(const std::vector< madness::Tensor<double> >& t,
            const int ispin, const xc_contribution xc_contrib) const;

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
        if (is_gga()) t[enum_saa]=0.5*rho; // sigma_aa=0
        madness::Tensor<double> f  = exc(t); //pending UGHHHHH
        madness::Tensor<double> va = vxc(t,0,XCfunctional::xc_contribution::potential_rho);
        madness::Tensor<double> vb = vxc(t,0,XCfunctional::xc_contribution::potential_same_spin);
        for (long i=0; i<npt; i++) {
            printf("%.3e %.3e %.3e %.3e\n", rho[i], f[i], va[i], vb[i]);
        }
    }
};

/// Class to compute the energy functional
struct xc_functional {
    const XCfunctional* xc;

    xc_functional(const XCfunctional& xc) : xc(&xc) {}

    madness::Tensor<double> operator()(const madness::Key<3> & key,
            const std::vector< madness::Tensor<double> >& t) const {
        MADNESS_ASSERT(xc);
        return xc->exc(t);
    }
};

/// Class to compute terms of the potential
struct xc_potential {
    const XCfunctional* xc;
    const XCfunctional::xc_contribution what;
    const int ispin;

    xc_potential(const XCfunctional& xc, int ispin, XCfunctional::xc_contribution what)
        : xc(&xc), what(what), ispin(ispin)
    {}

    madness::Tensor<double> operator()(const madness::Key<3> & key,
            const std::vector< madness::Tensor<double> >& t) const {
        MADNESS_ASSERT(xc);
        madness::Tensor<double> r = xc->vxc(t, ispin, what);
        return r;
    }
};


/// Class to compute terms of the kernel
struct xc_kernel_apply {
    const XCfunctional* xc;
    const int ispin;
    const XCfunctional::xc_contribution xc_contrib;
    const FunctionCommonData<double,3>& cdata;

    xc_kernel_apply(const XCfunctional& xc, int ispin,
            const XCfunctional::xc_contribution xc_contrib)
        : xc(&xc), ispin(ispin), xc_contrib(xc_contrib), cdata(FunctionCommonData<double,3>::get(FunctionDefaults<3>::get_k())) {}

    madness::Tensor<double> operator()(const madness::Key<3> & key,
            const std::vector< madness::Tensor<double> >& t) const {
        MADNESS_ASSERT(xc);
        madness::Tensor<double> r = xc->fxc_apply(t, ispin, xc_contrib);
        return r;
    }
};

}
#endif
