
#include <moldft/xcfunctional.h>

#ifndef HAVE_LIBXC

#include <tensor/tensor.h>
#include <sstream>

int x_rks_s__(const double *r__, double *f, double * dfdra);
int c_rks_vwn5__(const double *r__, double *f, double * dfdra);
int x_uks_s__(double *ra, double *rb, double *f, double *dfdra, double *dfdrb);
int c_uks_vwn5__(double *ra, double *rb, double *f, double *dfdra, double *dfdrb);

XCfunctional::XCfunctional() {}

void XCfunctional::initialize(const std::string& input_line, bool polarized) 
{
    spin_polarized = polarized;
    
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
}

XCfunctional::~XCfunctional() {}

bool XCfunctional::is_lda() const {
    return (hf_coeff == 0.0);
}
        
bool XCfunctional::is_gga() const {
    return false;
}
    
bool XCfunctional::is_meta() const {
    return false;
}

bool XCfunctional::is_dft() const {
    return (is_lda() || is_gga() || is_meta());
}
    
bool XCfunctional::has_fxc() const 
{
    return false;
}

/// returns true if the third derivative of the functional is available
bool XCfunctional::has_kxc() const
{
    return false;
}

madness::Tensor<double> XCfunctional::exc(const std::vector< madness::Tensor<double> >& t) const 
{
    const double* arho = t[0].ptr();
    madness::Tensor<double> result(3L, t[0].dims(), false);
    double* f = result.ptr();
    if (spin_polarized) {
        MADNESS_ASSERT(t.size() == 2);
        const double* brho = t[1].ptr();
        for (unsigned int i=0; i<result.size(); i++) {
            double ra = munge_rho(arho[i]); 
            double rb = munge_rho(brho[i]); 
            double xf, cf, xdfdr[2], cdfdr[2];
            
            x_uks_s__(&ra, &rb, &xf, xdfdr, xdfdr+1);
            c_uks_vwn5__(&ra, &rb, &cf, cdfdr, cdfdr+1);
            
            f[i] = xf + cf;
        }
    }
    else {
        MADNESS_ASSERT(t.size() == 1);
        double q1, q2, dq;
        for (unsigned int i=0; i<result.size(); i++) {
            double r = munge_rho(2.0 * arho[i]); 
            x_rks_s__(&r, &q1, &dq);
            c_rks_vwn5__(&r, &q2, &dq);
            f[i] = q1 + q2; 
        }
    }
    return result;
}

madness::Tensor<double> XCfunctional::vxc(const std::vector< madness::Tensor<double> >& t, const int ispin) const 
{
    const double* arho = t[0].ptr();
    madness::Tensor<double> result(3L, t[0].dims(), false);
    double* f = result.ptr();
    
    if (spin_polarized) {
        MADNESS_ASSERT(t.size() == 2);
        const double* brho = t[1].ptr();
        for (unsigned int i=0; i<result.size(); i++) {
            double ra = munge_rho(arho[i]); 
            double rb = munge_rho(brho[i]); 
            double xf, cf, xdfdr[2], cdfdr[2];
            
            x_uks_s__(&ra, &rb, &xf, xdfdr, xdfdr+1);
            c_uks_vwn5__(&ra, &rb, &cf, cdfdr, cdfdr+1);
            
            f[i] = xdfdr[ispin] + cdfdr[ispin];
        }
    }
    else {
        MADNESS_ASSERT(t.size() == 1);
        const double* arho = t[0].ptr();
        for (unsigned int i=0; i<result.size(); i++) {
            double r = munge_rho(2.0 * arho[i]); 
            double q, dq1, dq2;
            x_rks_s__(&r, &q, &dq1);
            c_rks_vwn5__(&r, &q, &dq2); 
            f[i] = dq1 + dq2;
        }
    }
    return result;
}


#endif
