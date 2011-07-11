#ifndef MOLDFT_XCMOLDFT_H
#define MOLDFT_XCMOLDFT_H

/// \file moldft/xcfunctional.h
/// \brief Defines interface for DFT XC functionals
/// \ingroup moldft

#include <iostream>
#include <sstream>

extern int x_rks_s__(const double *rho, double *f, double *dfdra);

extern int c_rks_vwn5__(const double *rho, double *f, double *dfdra);

extern int x_uks_s__(double *ra, double *rb, double *f,
                     double *dfdra, double *dfdrb);

extern int c_uks_vwn5__(double *ra, double *rb, double *
                        f, double *dfdra, double *dfdrb);


/// Simplified interface to XC functionals
class XCfunctional {
protected:
    bool initialized;
    bool spin_polarized;
    int nderiv;
    double hf_coeff;

    static double munge_rho(double r)  {
        if(r < 1e-12)  r = 1e-12;
        return r;
    }

public:
    /// Default constructor is required
    XCfunctional()
        : initialized(false)
    {}

    /// Initialize the object from the user input data
    void initialize(const std::string& input_line, bool polarized) 
    {
        spin_polarized = polarized;
        nderiv = 0;
        
        std::stringstream s(input_line);
        std::string token;
        while (s >> token) {
            if (token == "lda") {
                hf_coeff = 0.0;
            }
            else if (token == "hf") {
                hf_coeff = 1.0;
            }
        }

        initialized = true;
    }
    
    ~XCfunctional() 
    {}

    /// returns true if the potential is lda
    bool is_lda() const {
        return (hf_coeff == 0.0);
    }
        
    /// returns true if the potential is gga (needs first derivatives)
    bool is_gga() const {
        return false;
    }
    
    /// returns true if the potential is meta gga (needs second derivatives)
    bool is_meta_gga() const {
        return false;
    }
    
    /// returns true if the functional is spin_polarized
    bool is_spin_polarized() const 
    {
        return spin_polarized;
    }

    /// returns true if the second derivative of the functional is available
    bool has_fxc() const 
    {
        return false;
    }

    /// returns true if the third derivative of the functional is available
    bool has_kxc() const
    {
        return false;
    }

    /// returns the value of the hf exact exchange coefficient
    double hf_exchange_coefficient() const 
    {
        return hf_coeff;
    }

    /// computes the energy functional at np points

    /// array arguments are only referenced if they are used, so when
    /// evaluating an lda you need not provide derivatives, and for a
    /// gga need not provide second derivatives.
    ///
    /// for non-polarized functionals pass the values for the alpha
    /// spin component ... the beta values are not referenced.
    /// 
    /// any hf exchange contribution must be separately computed.
    void exc(double* f,
             int np, 
             const double* arho, const double* brho, 
             const double* darho=0, const double* dbrho=0, 
             const double* d2arho=0, const double* d2brho=0) const 
    {
        if (hf_coeff == 1.0) {
            for (int i=0; i<np; i++) {
                f[i] = 0.0;
            }
        }
        else if (spin_polarized) {
            throw "not yet";
        }
        else {
            double q1, q2, dq;
            for (int i=0; i<np; i++) {
                double r = munge_rho(2.0 * arho[i]); 
                x_rks_s__(&r, &q1, &dq);
                c_rks_vwn5__(&r, &q2, &dq);
                f[i] = q1 + q2;
            }
        }
    }
    

    /// computes the potential (derivative of the energy functional) at np points

    /// array arguments are only referenced if they are used, so when
    /// evaluating an lda you need not provide derivatives, and for a
    /// gga need not provide second derivatives.
    ///
    /// for non-polarized functionals pass the values for the alpha
    /// spin component ... the beta values are not referenced.
    /// 
    /// any hf exchange contribution must be separately computed.
    void vxc(double* f, 
             int np, 
             const double* arho, const double* brho, 
             const double* darho=0, const double* dbrho=0, 
             const double* d2arho=0, const double* d2brho=0) const 
    {
        if (hf_coeff == 1.0) {
            for (int i=0; i<np; i++) {
                f[i] = 0.0;
            }
        }
        else if (spin_polarized) {
            throw "not yet";
        }
        else {
            for (int i=0; i<np; i++) {
                double r = munge_rho(2.0 * arho[i]); 
                double q, dq1, dq2;
                x_rks_s__(&r, &q, &dq1);
                c_rks_vwn5__(&r, &q, &dq2); 
                f[i] = dq1 + dq2;
            }
        }
    }
};


struct xc_lda_functional {
    const XCfunctional* xc;
    
    xc_lda_functional(const XCfunctional& xc) 
        : xc(&xc)
    {}

    void operator()(const Key<3> & key, Tensor<double>& t) const 
    {
        xc->exc(t.ptr(), t.size(), t.ptr(), t.ptr());
        t.scale(0.5);
    }
};

struct xc_lda_potential {
    const XCfunctional* xc;

    xc_lda_potential(const XCfunctional& xc) 
        : xc(&xc)
    {}

    void operator()(const Key<3> & key, Tensor<double>& t) const 
    {
        xc->vxc(t.ptr(), t.size(), t.ptr(), t.ptr());
    }
};



#endif
