#include <madness_config.h>

#ifndef MADNESS_HAS_LIBXC

//#ifndef HAS_LIBMADXC

#include <DFcode/xcfunctional.h>
#include <tensor/tensor.h>
#include <sstream>
#include <libMADxc.h>


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

bool XCfunctional::has_kxc() const
{
    return false;
}

madness::Tensor<double> XCfunctional::exc(const std::vector< madness::Tensor<double> >& t, const int ispin) const
{

    madness::Tensor<double> result(3L, t[0].dims(), false);
    const int  np = result.size();

    int deriv = 0;

    if (spin_polarized) {
        MADNESS_ASSERT(t.size() == 2);

        madness::Tensor<double> sigmaaa1(np);
        madness::Tensor<double> sigmabb1(np);
        madness::Tensor<double> sigmaab1(np);

        madness::Tensor<double> vrhoa(np);
        madness::Tensor<double> vrhob(np);

        madness::Tensor<double> vsigmaaa(np);
        madness::Tensor<double> vsigmabb(np);
        madness::Tensor<double> vsigmaab(np);

        madness::Tensor<double> v2rhoa2(np);
        madness::Tensor<double> v2rhob2(np);
        madness::Tensor<double> v2rhoab(np);

        madness::Tensor<double> v2rhoasigmaaa(np);
        madness::Tensor<double> v2rhoasigmaab(np);
        madness::Tensor<double> v2rhoasigmabb(np);

        madness::Tensor<double> v2rhobsigmabb(np);
        madness::Tensor<double> v2rhobsigmaab(np);
        madness::Tensor<double> v2rhobsigmaaa(np);

        madness::Tensor<double> v2sigmaaa2(np);
        madness::Tensor<double> v2sigmaaaab(np);
        madness::Tensor<double> v2sigmaaabb(np);

        madness::Tensor<double> v2sigmaab2(np);
        madness::Tensor<double> v2sigmaabbb(np);
        madness::Tensor<double> v2sigmabb2(np);

        const double * restrict rhoa = t[0].ptr();
        const double * restrict rhob = t[1].ptr();

        madness::Tensor<double> rho1(np);
        madness::Tensor<double> rho2(np);
        double * restrict densa = rho1.ptr();
        double * restrict densb = rho2.ptr();

        for (long j=0; j<np; j++) {
              double ra=rhoa[j];
              double rb=rhob[j];
              densa[j]=munge(ra);
              densb[j]=munge(rb);
        }

        madness::Tensor<double> qx(3L, t[0].dims(), false);
        madness::Tensor<double> qc(3L, t[0].dims(), false);

        uks_x_lda_(&deriv, &np, densa, densb,
                    sigmaaa1.ptr(),sigmabb1.ptr(),sigmaab1.ptr(),
                    qx.ptr(),vrhoa.ptr(),vrhob.ptr(),vsigmaaa.ptr(),vsigmabb.ptr(),vsigmaab.ptr(),
                    v2rhoa2.ptr(),v2rhob2.ptr(),v2rhoab.ptr(),
                    v2rhoasigmaaa.ptr(),v2rhoasigmaab.ptr(),v2rhoasigmabb.ptr(),
                    v2rhobsigmabb.ptr(),v2rhobsigmaab.ptr(),v2rhobsigmaaa.ptr(),
                    v2sigmaaa2.ptr(),v2sigmaaaab.ptr(),v2sigmaaabb.ptr(),
                    v2sigmaab2.ptr(),v2sigmaabbb.ptr(),v2sigmabb2.ptr());
        uks_c_vwn5_(&deriv, &np, densa, densb,
                    sigmaaa1.ptr(),sigmabb1.ptr(),sigmaab1.ptr(),
                    qc.ptr(),vrhoa.ptr(),vrhob.ptr(),vsigmaaa.ptr(),vsigmabb.ptr(),vsigmaab.ptr(),
                    v2rhoa2.ptr(),v2rhob2.ptr(),v2rhoab.ptr(),
                    v2rhoasigmaaa.ptr(),v2rhoasigmaab.ptr(),v2rhoasigmabb.ptr(),
                    v2rhobsigmabb.ptr(),v2rhobsigmaab.ptr(),v2rhobsigmaaa.ptr(),
                    v2sigmaaa2.ptr(),v2sigmaaaab.ptr(),v2sigmaaabb.ptr(),
                    v2sigmaab2.ptr(),v2sigmaabbb.ptr(),v2sigmabb2.ptr());
        result = qx + qc;
    }
    else {
        MADNESS_ASSERT(t.size() == 1);

        madness::Tensor<double> sigmaaa(np);
        madness::Tensor<double> vrhoa(np);
        madness::Tensor<double> vsigmaaa(np);
        madness::Tensor<double> v2rhoa2(np);
        madness::Tensor<double> v2rhoasigmaaa(np);
        madness::Tensor<double> v2sigmaaa2(np);

        const double * restrict rhoa = t[0].ptr();

        madness::Tensor<double> rho(np);
        double * restrict dens = rho.ptr();

        for (long j=0; j<np; j++) {
              double ra=rhoa[j];
              dens[j]=munge(2.0*ra);
        }

        madness::Tensor<double> qx(3L, t[0].dims(), false);
        madness::Tensor<double> qc(3L, t[0].dims(), false);

        rks_x_lda_( &deriv, &np , dens, sigmaaa.ptr(), qx.ptr(), 
                    vrhoa.ptr(), vsigmaaa.ptr(), v2rhoa2.ptr(), v2rhoasigmaaa.ptr(), v2sigmaaa2.ptr());

        rks_c_vwn5_( &deriv, &np , dens, sigmaaa.ptr(), qc.ptr(), 
                    vrhoa.ptr(), vsigmaaa.ptr(), v2rhoa2.ptr(), v2rhoasigmaaa.ptr(), v2sigmaaa2.ptr());

        result = qx + qc;
    }
    return result;
}

madness::Tensor<double> XCfunctional::vxc(const std::vector< madness::Tensor<double> >& t, const int ispin, const int what) const
{
    //MADNESS_ASSERT(what == 0);

    madness::Tensor<double> result(3L, t[0].dims(), false);

    const int  np = result.size();

    int deriv = 1;
    
    if (spin_polarized) {
        MADNESS_ASSERT(t.size() == 2);
        madness::Tensor<double> sigmaaa1(np);
        madness::Tensor<double> sigmabb1(np);
        madness::Tensor<double> sigmaab1(np);

        madness::Tensor<double> vrhoax(3L, t[0].dims(), false);
        madness::Tensor<double> vrhobx(3L, t[0].dims(), false);
        madness::Tensor<double> vrhoac(3L, t[0].dims(), false);
        madness::Tensor<double> vrhobc(3L, t[0].dims(), false);

        madness::Tensor<double> vsigmaaa(np);
        madness::Tensor<double> vsigmabb(np);
        madness::Tensor<double> vsigmaab(np);

        madness::Tensor<double> v2rhoa2(np);
        madness::Tensor<double> v2rhob2(np);
        madness::Tensor<double> v2rhoab(np);

        madness::Tensor<double> v2rhoasigmaaa(np);
        madness::Tensor<double> v2rhoasigmaab(np);
        madness::Tensor<double> v2rhoasigmabb(np);

        madness::Tensor<double> v2rhobsigmabb(np);
        madness::Tensor<double> v2rhobsigmaab(np);
        madness::Tensor<double> v2rhobsigmaaa(np);

        madness::Tensor<double> v2sigmaaa2(np);
        madness::Tensor<double> v2sigmaaaab(np);
        madness::Tensor<double> v2sigmaaabb(np);

        madness::Tensor<double> v2sigmaab2(np);
        madness::Tensor<double> v2sigmaabbb(np);
        madness::Tensor<double> v2sigmabb2(np);

        const double * restrict rhoa = t[0].ptr();
        const double * restrict rhob = t[1].ptr();

        madness::Tensor<double> rho1(np);
        madness::Tensor<double> rho2(np);
        double * restrict densa = rho1.ptr();
        double * restrict densb = rho2.ptr();

        for (long j=0; j<np; j++) {
              double ra=rhoa[j];
              double rb=rhob[j];
              densa[j]=munge(ra);
              densb[j]=munge(rb);
        }

        madness::Tensor<double> qx(3L, t[0].dims(), false);
        madness::Tensor<double> qc(3L, t[0].dims(), false);

        uks_x_lda_(&deriv, &np, densa, densb,
                    sigmaaa1.ptr(),sigmabb1.ptr(),sigmaab1.ptr(),
                    qx.ptr(),vrhoax.ptr(),vrhobx.ptr(),vsigmaaa.ptr(),vsigmabb.ptr(),vsigmaab.ptr(),
                    v2rhoa2.ptr(),v2rhob2.ptr(),v2rhoab.ptr(),
                    v2rhoasigmaaa.ptr(),v2rhoasigmaab.ptr(),v2rhoasigmabb.ptr(),
                    v2rhobsigmabb.ptr(),v2rhobsigmaab.ptr(),v2rhobsigmaaa.ptr(),
                    v2sigmaaa2.ptr(),v2sigmaaaab.ptr(),v2sigmaaabb.ptr(),
                    v2sigmaab2.ptr(),v2sigmaabbb.ptr(),v2sigmabb2.ptr());
        uks_c_vwn5_(&deriv, &np, densa, densb,
                    sigmaaa1.ptr(),sigmabb1.ptr(),sigmaab1.ptr(),
                    qc.ptr(),vrhoac.ptr(),vrhobc.ptr(),vsigmaaa.ptr(),vsigmabb.ptr(),vsigmaab.ptr(),
                    v2rhoa2.ptr(),v2rhob2.ptr(),v2rhoab.ptr(),
                    v2rhoasigmaaa.ptr(),v2rhoasigmaab.ptr(),v2rhoasigmabb.ptr(),
                    v2rhobsigmabb.ptr(),v2rhobsigmaab.ptr(),v2rhobsigmaaa.ptr(),
                    v2sigmaaa2.ptr(),v2sigmaaaab.ptr(),v2sigmaaabb.ptr(),
                    v2sigmaab2.ptr(),v2sigmaabbb.ptr(),v2sigmabb2.ptr());

        result = vrhoax + vrhoac + vrhobx + vrhobc ;
    }
    else {
        MADNESS_ASSERT(t.size() == 1);

        madness::Tensor<double> sigmaaa(np);
        madness::Tensor<double> vsigmaaa(np);
        madness::Tensor<double> v2rhoa2(np);
        madness::Tensor<double> v2rhoasigmaaa(np);
        madness::Tensor<double> v2sigmaaa2(np);

        madness::Tensor<double> qx(3L, t[0].dims(), false);
        madness::Tensor<double> qc(3L, t[0].dims(), false);

        madness::Tensor<double> vrhoax(3L, t[0].dims(), false);
        madness::Tensor<double> vrhoac(3L, t[0].dims(), false);

        const double * restrict rhoa = t[0].ptr();

        madness::Tensor<double> rho(np);
        double * restrict dens = rho.ptr();

        for (long j=0; j<np; j++) {
              double ra=rhoa[j];
              dens[j]=munge(2.0*ra);
        }

        rks_x_lda_( &deriv, &np , dens, sigmaaa.ptr(), qx.ptr(), 
                    vrhoax.ptr(), vsigmaaa.ptr(), v2rhoa2.ptr(), v2rhoasigmaaa.ptr(), v2sigmaaa2.ptr());

        rks_c_vwn5_( &deriv, &np , dens, sigmaaa.ptr(), qc.ptr(), 
                    vrhoac.ptr(), vsigmaaa.ptr(), v2rhoa2.ptr(), v2rhoasigmaaa.ptr(), v2sigmaaa2.ptr());

        result = vrhoax + vrhoac;
    }
    return result;
}


#endif
