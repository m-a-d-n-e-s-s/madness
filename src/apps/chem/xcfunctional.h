#ifndef MADNESS_CHEM_XCFUNCTIONAL_H__INCLUDED
#define MADNESS_CHEM_XCFUNCTIONAL_H__INCLUDED

#include <madness/madness_config.h>

/// \file moldft/xcfunctional.h
/// \brief Defines interface for DFT XC functionals
/// \ingroup moldft

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

    /// the ordering of the intermediates is fixed, but the code can handle
    /// non-initialized functions, so if e.g. no GGA is requested, all the
    /// corresponding vector components may be left empty.
    enum xc_arg {
        enum_rhoa=0,            ///< alpha density \f$ \rho_\alpha \f$
        enum_rhob=1,            ///< \f$ \rho_\beta \f$
        enum_rho_pt=2,          ///< perturbed density (CPHF, TDKS) \f$ \rho_{pt} \f$
        enum_drhoa_x=4,         ///< \f$ \partial/{\partial x} \rho_{\alpha} \f$
        enum_drhoa_y=5,         ///< \f$ \partial/{\partial y} \rho_{\alpha} \f$
        enum_drhoa_z=6,         ///< \f$ \partial/{\partial z} \rho_{\alpha} \f$
        enum_drhob_x=7,         ///< \f$ \partial/{\partial x} \rho_{\beta} \f$
        enum_drhob_y=8,         ///< \f$ \partial/{\partial y} \rho_{\beta} \f$
        enum_drhob_z=9,         ///< \f$ \partial/{\partial z} \rho_{\beta} \f$
        enum_saa=10,            ///< \f$ \sigma_{aa} = \nabla \rho_{\alpha}.\nabla \rho_{\alpha} \f$
        enum_sab=11,            ///< \f$ \sigma_{ab} = \nabla \rho_{\alpha}.\nabla \rho_{\beta} \f$
        enum_sbb=12,            ///< \f$ \sigma_{bb} = \nabla \rho_{\beta}.\nabla \rho_{\beta} \f$
        enum_sigtot=13,         ///< \f$ \sigma = \nabla \rho.\nabla \rho \f$
        enum_sigma_pta=14,      ///< \f$ \nabla\rho_{\alpha}.\nabla\rho_{pt} \f$
        enum_sigma_ptb=15       ///< \f$ \nabla\rho_{\beta}.\nabla\rho_{pt} \f$
    };
    const static int number_xc_args=16;     ///< max number of intermediates


    /// which contribution is requested from the XCfunctional

    /// three types of contributions are required:
    /// those that are second derivatives of the xc kernel and that are local
    /// those that are first/second derivatives of the kernel and of which the
    /// derivative will be taken later on (-> semi-local)
    enum xc_contribution {
        potential_rho,              // potential df/drho
        potential_same_spin,        // potential df/dsigma_aa
        potential_mixed_spin,       // potential df/dsigma_ab
        kernel_second_local,        // kernel (2nd derivative)
        kernel_second_semilocal,    // kernel (2nd derivative)
        kernel_first_semilocal      // kernel (2nd derivative)
    };

    /// different munging for potential and for kernel
    enum munging_type {xc_potential, xc_kernel=1};

    double get_rhotol() const {return rhotol;}
    double get_ggatol() const {return ggatol;}

protected:

    bool spin_polarized;        ///< True if the functional is spin polarized
    double hf_coeff;            ///< Factor multiplying HF exchange (+1.0 gives HF)
    double rhomin, rhotol, sigmin, sigtol; // See initialize and munge*
    double ggatol; // See initialize and munge*
    double munge_ratio; // See initialize and munge*

#ifdef MADNESS_HAS_LIBXC
    std::vector< std::pair<xc_func_type*,double> > funcs;

    void make_libxc_args(const std::vector< madness::Tensor<double> >& t,
                         madness::Tensor<double>& rho,
                         madness::Tensor<double>& sigma,
                         const munging_type& munging) const;
    void make_libxc_args_old(const std::vector< madness::Tensor<double> >& t,
                         madness::Tensor<double>& rho,
                         madness::Tensor<double>& sigma,
                         const munging_type& munging) const;

    /// the number of xc kernel derivatives (lda: 0, gga: 1, etc)
    int nderiv;
#endif


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

    struct munger {
        munger(const double rhotol1, const double rhomin1)
            : rhotol(rhotol1), rhomin(rhomin1) {}
        double operator()(double rho) const {
//            if(rho<1.e-7) rho = 1.e-7;
            if (fabs(rho) <= rhotol) rho=rhomin;
            return rho;
        }
        double rhotol,rhomin;
    };

    /// simple munging for the density only (LDA)
    double munge(double rho) const {
//        if(rho<1.e-7) rho = 1.e-7;
    	if (fabs(rho) <= rhotol) rho=rhomin;
        return rho;
    }

    /// similar to the Laura's ratio thresholding, but might be more robust
    void munge_xc_kernel(double& rho, double& sigma) const {
        if (sigma<0.0) sigma=sigmin;
        if (rho < rhotol) rho=rhomin;                   // 1.e-8 or so
        if (rho < 1.e-2) sigma=munge_ratio*rho*rho;
    }

    /// new 'ratio' threshold'
    /// still need to ensure rho and sigma don't go negative
    void munge_laura(double& rho, double& sigma) const {
        if (rho < 0.0 || sigma < 0.0 ||
            (rho<1e-2 && (sigma/(rho*rho)>10000.0))) {
            rho = rhomin;
            sigma = sigmin;
        }
//        if ((rho<rhotol) or (sigma<(rhotol*rhotol*100))) {
//            rho=rhomin;
//            sigma=sigmin;
//        }
//        if (sigma<sigtol) sigma=sigmin;
//        if (sigma/(rho*rho)) {
//            rho=rhomin;
//            sigma=sigmin;
//        }
    }

    void munge2(double& rho, double& sigma, const munging_type& munging) const {
        if (munging==xc_potential) {
            munge_laura(rho,sigma);
        } else if (munging==xc_kernel) {
            munge_xc_kernel(rho,sigma);
        } else {
            MADNESS_EXCEPTION("unknown munging type in xcfunctional.h",1);
        }
        return;

        /*if ( (0.5 * log10(sigma) - 2) > log10(rho) || rho < 0.0 || sigma < 0.0){
           //std::cout << "rho,sig " << rho << " " << sigma << " " << rhomin << " " << sigmin << std::endl;
           rho=rhomin;
           sigma=sigmin;
        }*/
    }

    void munge5(double& rhoa, double& rhob, double& saa, double& sab,
            double& sbb, const munging_type munging) const {
        munge2(rhoa, saa, munging);
        munge2(rhob, sbb, munging);
        if (rhoa==rhomin || rhob==rhomin) sab=sigmin;
    }

    /// munge rho and sigma to physical values

    /// since we use ratio thresholding we need the ratio of the density
    /// to the reducued density gradient sigma. There is no such thing
    /// for the mixed-spin sigma, therefore we use the identity
    /// \f[
    ///   \sigma_{total} = \nabla rho . \nabla \rho
    ///                  = \sigma_{aa} + 2\sigma_{ab} + \sigma_{bb}
    /// \f]
    /// so we can reconstruct sigma_ab from the spin densities (and density
    /// gradients) and the total densities (and density gradients)
    /// @param[in,out]  rhoa alpha spin density
    /// @param[in,out]  rhob beta spin density
    /// @param[in,out]  rho  total density
    /// @param[in,out]  saa  alpha spin reduced density gradient
    /// @param[in,out]  sab  mixed spin reduced density gradient
    /// @param[in,out]  sbb  beta spin reduced density gradient
    /// @param[in,out]  stot total reduced density gradient
    /// @param[in]  munging munging type
    void munge7(double& rhoa, double& rhob, double& rho,
            double& saa, double& sab, double& sbb, double stot,
            const munging_type munging) const {
        munge2(rhoa, saa, munging);
        munge2(rhob, sbb, munging);
        munge2(rho, stot, munging);
        sab=std::max(sigmin,0.5*(stot-saa-sbb));
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
    ///   < \phi | \hat V | \psi > = \int \left( \frac{\partial \epsilon}{\partial \rho_{\alpha}} \phi \psi
    ///                  +  \left( 2 \frac{\partial \epsilon}{\partial \sigma_{\alpha \alpha}} + \frac{\partial \epsilon}{\partial \sigma_{\alpha \beta}} \right) \nabla \rho_{\alpha} . \nabla \left( \phi \psi \right) \right) dx
    /// \f$
    ///
    /// For polarized GGA, matrix elements of the potential are
    /// \f$
    ///   < \phi_{\alpha} | \hat V | \psi_{\alpha} > = \int \left( \frac{\partial \epsilon}{\partial \rho_{\alpha}} \phi \psi
    ///                  +  \left( 2 \frac{\partial \epsilon}{\partial \sigma_{\alpha \alpha}} \nabla \rho_{\alpha}  + \frac{\partial \epsilon}{\partial \sigma_{\alpha \beta}} \nabla \rho_{\beta}  \right) . \nabla \left( \phi \psi \right) \right) dx
    /// \f$
    ///
    /// Integrating the above by parts and assuming free-space or periodic boundary conditions
    /// we obtain that the local multiplicative form of the GGA potential is
    /// \f$
    ///    V_{\alpha} =  \frac{\partial \epsilon}{\partial \rho_{\alpha}} - \left(\nabla . \left(2 \frac{\partial \epsilon}{\partial \sigma_{\alpha \alpha}} \nabla \rho_{\alpha}  + \frac{\partial \epsilon}{\partial \sigma_{\alpha \beta}} \nabla \rho_{\beta}  \right)  \right)
    /// \f$
    ///
    /// Until we get a madness::Function operation that can produce
    /// multiple results we need to compute components of the
    /// functional and potential separately:
    /// @param[in] t The input densities and derivatives as required by the functional
    /// @param[in] ispin Specifies which component of the potential is to be computed as described above
    /// @return The component specified by the \c ispin parameter
    /// \todo Matt changed the parameter "what" (no longer listed below) to "ispin". Is this correct?
    madness::Tensor<double> vxc(const std::vector< madness::Tensor<double> >& t,
            const int ispin, const xc_contribution xc_contrib) const;

    /// compute the second derivative of the XC energy wrt the density and apply

    /// apply the kernel on the fly on the provided (response) density
    /// @param[in] t The input densities and derivatives as required by the functional
    /// @return The component specified by the \c what parameter
    madness::Tensor<double> fxc_apply(const std::vector< madness::Tensor<double> >& t,
            const int ispin, const xc_contribution xc_contrib) const;

    madness::Tensor<double> fxc_old(const std::vector< madness::Tensor<double> >& t,
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
