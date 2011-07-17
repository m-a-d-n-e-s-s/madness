#ifndef MOLDFT_XCMOLDFT_H
#define MOLDFT_XCMOLDFT_H

// libxc is currently only working for spin-restricted LDA
//#define HAVE_LIBXC 1

/// \file moldft/xcfunctional.h
/// \brief Defines interface for DFT XC functionals
/// \ingroup moldft

#include <tensor/tensor.h>
#include <vector>
#include <algorithm>
#include <utility>
#include <mra/key.h>

#ifdef HAVE_LIBXC
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
    bool spin_polarized;
    double hf_coeff;

#ifdef HAVE_LIBXC
    std::vector< std::pair<xc_func_type*,double> > funcs;
    int nderiv;
#endif

    static double munge_rho(double r)  {
        if(r < 1e-12)  r = 1e-12;
        return r;
    }

public:
    /// Default constructor is required
    XCfunctional();

    /// Initialize the object from the user input data
    void initialize(const std::string& input_line, bool polarized);
    
    ~XCfunctional();

    /// returns true if the potential is lda
    bool is_lda() const;
        
    /// returns true if the potential is gga (needs first derivatives)
    bool is_gga() const;
    
    /// returns true if the potential is meta gga (needs second derivatives)
    bool is_meta() const;

    /// returns true if there is a functional (false probably means Hatree-Fock exchange only)
    bool is_dft() const;
    
    /// returns true if the functional is spin_polarized
    bool is_spin_polarized() const 
    {
        return spin_polarized;
    }

    /// returns true if the second derivative of the functional is available
    bool has_fxc() const;

    /// returns true if the third derivative of the functional is available
    bool has_kxc() const;

    /// returns the value of the hf exact exchange coefficient
    double hf_exchange_coefficient() const 
    {
        return hf_coeff;
    }

    /// computes the energy functional at given points

    /// lda     ---  t[0] = alpha density
    /// lda-pol ---  t[0] = alpha density, t[1] = beta density
    ///
    /// any hf exchange contribution must be separately computed.
    madness::Tensor<double> exc(const std::vector< madness::Tensor<double> >& t) const;
    
    /// computes the potential (derivative of the energy functional) at np points

    /// lda     ---  t[0] = alpha density
    /// lda-pol ---  t[0] = alpha density, t[1] = beta density
    ///
    /// any hf exchange contribution must be separately computed.
    madness::Tensor<double> vxc(const std::vector< madness::Tensor<double> >& t, const int ispin) const;

    void plot() const {
        long npt = 1001;
        double lo=1e-6, hi=1e+1, s=std::pow(hi/lo, 1.0/(npt-1));
        madness::Tensor<double> r(npt);
        for (int i=0; i<npt; i++) {
            r[i] = lo;
            lo *= s;
        }
        std::vector< madness::Tensor<double> > t;
        t.push_back(r);
        madness::Tensor<double> f = exc(t);
        for (long i=0; i<npt; i++) {
            printf("%.3e %.3e\n", r[i], f[i]);
        }
    }
};

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

struct xc_potential {
    const XCfunctional* xc;
    const int ispin;
    
    xc_potential(const XCfunctional& xc, int ispin) 
        : xc(&xc), ispin(ispin)
    {}
    
    madness::Tensor<double> operator()(const madness::Key<3> & key, const std::vector< madness::Tensor<double> >& t) const 
    {
        MADNESS_ASSERT(xc);
        madness::Tensor<double> r = xc->vxc(t, ispin);
        return r;
    }
};

#endif
