#ifndef MOLDFT_XCMOLDFT_H
#define MOLDFT_XCMOLDFT_H

#include <madness_config.h>

/// \file moldft/xcfunctional.h
/// \brief Defines interface for DFT XC functionals
/// \ingroup moldft

#include <tensor/tensor.h>
#include <vector>
#include <algorithm>
#include <utility>
#include <mra/key.h>

#ifdef MADNESS_HAS_LIBXC
#include <xc.h>
#endif

/// Compute the spin-restricted LDA potential using unaryop (only for the initial guess) 
struct xc_lda_potential {
    xc_lda_potential() {}
    
    void operator()(const madness::Key<3> & key, madness::Tensor<double>& t) const 
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
protected:
    bool spin_polarized;        ///< True if the functional is spin polarized
    double hf_coeff;            ///< Factor multiplying HF exchange (+1.0 gives HF)

#ifdef MADNESS_HAS_LIBXC
    std::vector< std::pair<xc_func_type*,double> > funcs;
    int nderiv;
#endif

    // For x<xmax, smoothly restricts x to be greater than or equal to xmin>0
    static double munge(double x, double xmin, double xmax) {
        // A quintic polyn that smoothly interpolates between
        //
        // x=0    where it has value xmin and zero slope, and
        //
        // x=xmax where it has value xmax, unit slope, and zero seond and third derivatives
        
        if (x > xmax) {
            return x; // most probable case
        }
        else if (x <= 0.0) {
            return xmin;
        }
        else {
            double xmax2 = xmax*xmax;
            double xmax3 = xmax2*xmax;
            double xmax4 = xmax2*xmax2;
            double xmax5 = xmax3*xmax2;

            return xmin+((-10.0*xmin+4.0*xmax)/xmax2+(-(-20.0*xmin+6.0*xmax)/xmax3+((-15.0*xmin+4.0*xmax)/xmax4-(-4.0*xmin+xmax)*x/xmax5)*x)*x)*x*x;
        }
    }

    // Uses same smoothing as above but smooths f (=x), and |del f|^2
    static void munge2(double& x, double& delxsq, double xmin, double xmax) {
        
        if (x > xmax) {
            return;
        }
        else if (x <= 0.0) {
            x = xmin;
            delxsq = 0.0;
        }
        else {
            double xmax2 = xmax*xmax;
            double xmax3 = xmax2*xmax;
            double xmax4 = xmax2*xmax2;
            double xmax5 = xmax3*xmax2;

            x = xmin+((-10.0*xmin+4.0*xmax)/xmax2+(-(-20.0*xmin+6.0*xmax)/xmax3+((-15.0*xmin+4.0*xmax)/xmax4-(-4.0*xmin+xmax)*x/xmax5)*x)*x)*x*x;

            // This is just the derivative of the above 
            double d = ((-20.0*xmin+8.0*xmax)/xmax2+((60.0*xmin-18.0*xmax)/xmax3+((-60.0*xmin+16.0*xmax)/xmax4+(20.0*xmin-5*xmax)*x/xmax5)*x)*x)*x;

            delxsq *= d*d;
        }
    }

    // Smoothly maps the density to a value slightly greater than zero
    static double munge_rho(double r)  {
        return munge(r, 1e-10, 1e-8);
        return r;
    }

    // Smooths density and also returns consistent scaling of sigma = |grad rho|^2
    static void munge_rho_sig(double& r, double& sig)  {
        return munge2(r, sig, 1e-10, 1e-8);
    }

public:
    /// Default constructor is required
    XCfunctional();

    /// Initialize the object from the user input data

    /// @param[in] input_line User input line (without beginning XC keyword)
    /// @param[in] polarized Boolean flag indicating if the calculation is spin-polarized
    void initialize(const std::string& input_line, bool polarized);
    
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

    /// This uses the convention that the total energy is \f$ E[\rho] = \int \epsilon[\rho(x)] dx\f$
    ///
    /// Any HF exchange contribution must be separately computed.
    ///
    /// Items in the vector argument \c t are interpreted as follows
    ///  - Spin un-polarized
    ///    - \c t[0] = \f$ \rho_{\alpha} \f$
    ///    - \c t[1/2/3] = \f$ \nabla \rho_{\alpha} \f$ (GGA only)
    ///  - Spin polarized
    ///    - \c t[0] = \f$ \rho_{\alpha} \f$
    ///    - \c t[1] = \f$ \rho_{\beta} \f$
    ///    - \c t[2/3/4] = \f$ \nabla \rho_{\alpha} \f$ (GGA only)
    ///    - \c t[5/6/7] = \f$ \nabla \rho_{\beta} \f$  (GGA only)
    ///
    /// @param t The input densities and derivatives as required by the functional
    /// @return The exchange-correlation energy functional
    madness::Tensor<double> exc(const std::vector< madness::Tensor<double> >& t) const;
    
    /// Computes components of the potential (derivative of the energy functional) at np points

    /// Any HF exchange contribution must be separately computed.
    ///
    /// See the documenation of the \c exc() method for contents of the input \c t[] argument
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
    /// 
    /// Until we get a madness::Function operation that can produce
    /// multiple results we need to compute components of the
    /// functional and potential separately:
    ///
    /// - Spin un-polarized
    ///   - \c what=0 \f$ \frac{\partial \epsilon}{\partial \rho_{\alpha}}\f$ 
    ///   - \c what=1 \f$ 2 \frac{\partial \epsilon}{\partial \sigma_{\alpha \alpha}} + \frac{\partial \epsilon}{\partial \sigma_{\alpha \beta}}\f$ (GGA only)
    /// - Spin polarized
    ///   - \c what=0 \f$ \frac{\partial \epsilon}{\partial \rho_{\alpha}}\f$ 
    ///   - \c what=1 \f$ \frac{\partial \epsilon}{\partial \rho_{\beta}}\f$ 
    ///   - \c what=2 \f$ \frac{\partial \epsilon}{\partial \sigma_{\alpha \alpha}} \f$
    ///   - \c what=3 \f$ \frac{\partial \epsilon}{\partial \sigma_{\alpha \beta}} \f$
    ///   - \c what=4 \f$ \frac{\partial \epsilon}{\partial \sigma_{\beta \beta}} \f$
    ///
    /// @param[in] t The input densities and derivatives as required by the functional
    /// @param[in] what Specifies which component of the potential is to be computed as described above
    /// @return The component specified by the \c what parameter
    madness::Tensor<double> vxc(const std::vector< madness::Tensor<double> >& t, const int what=0) const;

    /// Crude function to plot the energy and potential functionals
    void plot() const {
        long npt = 1001;
        double lo=1e-6, hi=1e+1, s=std::pow(hi/lo, 1.0/(npt-1));

        madness::Tensor<double> rho(npt);
        for (int i=0; i<npt; i++) {
            rho[i] = lo;
            lo *= s;
        }
        std::vector< madness::Tensor<double> > t;
        t.push_back(rho);
        if (is_spin_polarized()) t.push_back(rho);
        madness::Tensor<double> f  = exc(t);
        madness::Tensor<double> va = vxc(t,0);
        madness::Tensor<double> vb = vxc(t,1);
        for (long i=0; i<npt; i++) {
            printf("%.3e %.3e %.3e %.3e\n", rho[i], f[i], va[i], vb[i]);
        }
    }
};

/// Class to compute the energy functional
struct xc_functional {
    const XCfunctional* xc;

    xc_functional(const XCfunctional& xc) 
        : xc(&xc)
    {}
    
    madness::Tensor<double> operator()(const madness::Key<3> & key, const std::vector< madness::Tensor<double> >& t) const 
    {
        MADNESS_ASSERT(xc);
        return xc->exc(t);
    }
};

/// Class to compute terms of the potential
struct xc_potential {
    const XCfunctional* xc;
    const int what;
    
    xc_potential(const XCfunctional& xc, int what) 
        : xc(&xc), what(what)
    {}
    
    madness::Tensor<double> operator()(const madness::Key<3> & key, const std::vector< madness::Tensor<double> >& t) const 
    {
        MADNESS_ASSERT(xc);
        madness::Tensor<double> r = xc->vxc(t, what);
        return r;
    }
};

#endif
