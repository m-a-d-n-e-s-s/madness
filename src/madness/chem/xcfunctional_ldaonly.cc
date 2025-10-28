#include <madness/madness_config.h>

#ifndef MADNESS_HAS_LIBXC

#include<madness/chem/xcfunctional.h>
#include <madness/tensor/tensor.h>
#include <sstream>
#include <cmath>
#include <madness/world/MADworld.h>

namespace madness {


int x_rks_s__(const double *r__, double *f, double * dfdra);
int c_rks_vwn5__(const double *r__, double *f, double * dfdra);
int x_uks_s__(double *ra, double *rb, double *f, double *dfdra, double *dfdrb);
int c_uks_vwn5__(double *ra, double *rb, double *f, double *dfdra, double *dfdrb);

XCfunctional::XCfunctional() : hf_coeff(0.0) {
    rhotol=1e-7; rhomin=1e-12; // default values
}

void XCfunctional::initialize(const std::string& input_line, bool polarized,
        World& world, bool verbose) {
    rhotol=1e-7; rhomin=1e-12; // default values

    spin_polarized = polarized;

    std::stringstream s(input_line);
    std::string token;
    bool found_valid_token = false;
    while (s >> token) {
        std::transform(token.begin(), token.end(), token.begin(), ::toupper);
        if (token == "LDA") {
            hf_coeff = 0.0;
            found_valid_token = true;
        }
        else if (token == "RHOMIN") {
            s >> rhomin;
        }
        else if (token == "RHOTOL") {
            s >> rhotol;
        }
        else if (token == "HF") {
            hf_coeff = 1.0;
            found_valid_token = true;
        }
    }
    if (not found_valid_token)
        throw "XCfunctional(ldaonly)::initialize() -- did not find a valid XC functional";
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

bool XCfunctional::has_kxc() const
{
    return false;
}

madness::Tensor<double> XCfunctional::exc(const std::vector< madness::Tensor<double> >& t) const
{
    const double* arho = t[0].ptr();
    madness::Tensor<double> result(3L, t[0].dims(), false);
    double* f = result.ptr();
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wtautological-constant-compare"
#endif
    auto isnan = [](double v) { return std::isnan(v); };
#ifdef __clang__
#pragma clang diagnostic pop
#endif
    if (spin_polarized) {
        const double* brho = t[1].ptr();
        for (unsigned int i=0; i<result.size(); i++) {
            double ra = munge(arho[i]);
            double rb = munge(brho[i]);
            double xf, cf, xdfdr[2], cdfdr[2];

            x_uks_s__(&ra, &rb, &xf, xdfdr, xdfdr+1);
            c_uks_vwn5__(&ra, &rb, &cf, cdfdr, cdfdr+1);

            f[i] = xf + cf;
            if (isnan(f[i])) {
                print("bad 1?", ra, rb);
                throw "numerical error in lda functional";
            }
        }
    }
    else {
        double q1, q2, dq;
        for (unsigned int i=0; i<result.size(); i++) {
            double r = munge(2.0 * arho[i]);
            x_rks_s__(&r, &q1, &dq);
            c_rks_vwn5__(&r, &q2, &dq);
            f[i] = q1 + q2;
            if (isnan(f[i])) {
                print("bad? 2", r);
                throw "numerical error in lda functional";
            }
        }
    }
    return result;
}

std::vector<madness::Tensor<double> > XCfunctional::vxc(const std::vector< madness::Tensor<double> >& t,
        const int ispin) const
{
    //MADNESS_ASSERT(what == 0);
    const double* arho = t[0].ptr();
    std::vector<madness::Tensor<double> > result(1);
    result[0]=madness::Tensor<double>(3L, t[0].dims(), false);
    double* f = result[0].ptr();
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wtautological-constant-compare"
#endif
    auto isnan = [](double v) { return std::isnan(v); };
#ifdef __clang__
#pragma clang diagnostic pop
#endif

    if (spin_polarized) {
        const double* brho = t[1].ptr();
        for (unsigned int i=0; i<result[0].size(); i++) {
            double ra = munge(arho[i]);
            double rb = munge(brho[i]);
            double xf, cf, xdfdr[2], cdfdr[2];

            x_uks_s__(&ra, &rb, &xf, xdfdr, xdfdr+1);
            c_uks_vwn5__(&ra, &rb, &cf, cdfdr, cdfdr+1);

//            f[i] = xdfdr[what] + cdfdr[what];
            f[i] = xdfdr[ispin] + cdfdr[ispin];
            if (isnan(f[i])) {
                print("bad? 3", ra, rb);
                throw "numerical error in lda functional";
            }
        }
    }
    else {
        const double* arho = t[0].ptr();
        for (unsigned int i=0; i<result[0].size(); i++) {
            double r = munge(2.0 * arho[i]);
            double q, dq1, dq2;
            x_rks_s__(&r, &q, &dq1);
            c_rks_vwn5__(&r, &q, &dq2);
            f[i] = dq1 + dq2;
            if (isnan(f[i])) {
                print("bad? 4", r);
                throw "numerical error in lda functional";
            }
        }
    }
    return result;
}

std::vector<madness::Tensor<double> > XCfunctional::fxc_apply(const std::vector<Tensor<double> >& t,
        const int ispin) const{
	MADNESS_EXCEPTION("fxc_apply not implemented in xcfunctional_ldaonly.cc... use libxc",1);
}

void XCfunctional::make_libxc_args(const std::vector< madness::Tensor<double> >& t,
                        madness::Tensor<double>& rho,
                        madness::Tensor<double>& sigma,
                        madness::Tensor<double>& rho_pt,
                        madness::Tensor<double>& sigma_pt,
                        std::vector<madness::Tensor<double> >& drho,
                        std::vector<madness::Tensor<double> >& drho_pt,
                        const bool need_response) const {
    MADNESS_EXCEPTION("no make_libxc_args without libxc",1);
}

}
#endif
