#ifndef MOLDFT_XCMOLDFT_H
#define MOLDFT_XCMOLDFT_H

/// \file moldft/xcfunctional.h
/// \brief Defines interface for DFT XC functionals
/// \ingroup moldft

#include <iostream>
#include <sstream>
#include <tensor/tensor.h>

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
    bool is_meta() const {
        return false;
    }

    /// returns true if there is a functional (false probably means Hatree-Fock exchange only)
    bool is_dft() const {
        return (is_lda() || is_gga() || is_meta());
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

    /// computes the energy functional at given points

    /// lda     ---  t[0] = alpha density
    /// lda-pol ---  t[0] = alpha density, t[1] = beta density
    ///
    /// any hf exchange contribution must be separately computed.
    madness::Tensor<double> exc(const std::vector< madness::Tensor<double> >& t) const 
    {
        if (spin_polarized) {
            throw "not yet";
        }
        else {
            madness::Tensor<double> result(3L, t[0].dims(), false);
            double* f = result.ptr();
            const double* arho = t[0].ptr();

            double q1, q2, dq;
            for (unsigned int i=0; i<result.size(); i++) {
                double r = munge_rho(2.0 * arho[i]); 
                x_rks_s__(&r, &q1, &dq);
                c_rks_vwn5__(&r, &q2, &dq);
                f[i] = q1 + q2;
            }
            return result;
        }
    }
    
    
    /// computes the potential (derivative of the energy functional) at np points

    /// lda     ---  t[0] = alpha density
    /// lda-pol ---  t[0] = alpha density, t[1] = beta density
    ///
    /// any hf exchange contribution must be separately computed.
    madness::Tensor<double> vxc(const std::vector< madness::Tensor<double> >& t, int ispin) const 
    {
        if (spin_polarized) {
            throw "not yet";
        }
        else {
            madness::Tensor<double> result(3L, t[0].dims(), false);
            double* f = result.ptr();
            const double* arho = t[0].ptr();
            for (unsigned int i=0; i<result.size(); i++) {
                double r = munge_rho(2.0 * arho[i]); 
                double q, dq1, dq2;
                x_rks_s__(&r, &q, &dq1);
                c_rks_vwn5__(&r, &q, &dq2); 
                f[i] = dq1 + dq2;
            }
            return result;
        }
    }
};

struct xc_functional {
    const XCfunctional* xc;

    xc_functional(const XCfunctional& xc) 
        : xc(&xc)
    {}
    
    madness::Tensor<double> operator()(const Key<3> & key, const std::vector< madness::Tensor<double> >& t) const 
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

    madness::Tensor<double> operator()(const Key<3> & key, const std::vector< madness::Tensor<double> >& t) const 
    {
        MADNESS_ASSERT(xc);
        madness::Tensor<double> r = xc->vxc(t, ispin);
        //std::cout << key << " " << t[0].sumsq() << " " << r.sumsq() << std::endl;
        return r;
    }
};




#endif
