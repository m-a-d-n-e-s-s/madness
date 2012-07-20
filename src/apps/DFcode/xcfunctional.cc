#include <madness_config.h>

#ifndef MADNESS_HAS_LIBXC

//#ifndef HAS_LIBMADXC

#include <DFcode/xcfunctional.h>
#include <tensor/tensor.h>
#include <sstream>

#include <libMADxc.h>


XCfunctional::XCfunctional() {}
// This function is here because compiling with -ffast-math breaks the C
// function isnan. Also, there is a bug in some compilers where isnan is
// undefined when <cmath> is included.
namespace {
    inline int isnan_x(double x) {
        volatile double y = x;
        return x != y;
    }
}

void XCfunctional::initialize(const std::string& input_line, bool polarized) 
{
    spin_polarized = polarized;
    
 
    std::stringstream s(input_line);
    std::string token;
    while (s >> token) {
        if (token == "lda") {
            lda = true;
        }
        if (token == "pbe") {
            gga = true;
        }
        else if (token == "hf") {
            hf_coeff = 1.0;
        }
    }
    hf_coeff=0.0;
    lda=false;
    gga=false;

    gga=true;
}

XCfunctional::~XCfunctional() {}

bool XCfunctional::is_lda() const {
    return (lda);
}
        
bool XCfunctional::is_gga() const {
    return (gga);
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

/// Allocates rho (and if GGA also sigma) and copies data from t[] into rho and sigma.
void XCfunctional::make_xc_args(const std::vector< madness::Tensor<double> >& t,
                                std::vector<madness::Tensor<double> >& rho, 
                                std::vector<madness::Tensor<double> >& sigma, 
                                std::vector<madness::Tensor<double> >& delrho) const
{
    const int np = t[0].size();
    if (spin_polarized) {
        if (is_lda()) {
            MADNESS_ASSERT(t.size() == 2);
            const double * restrict rhoa = t[0].ptr();
            const double * restrict rhob = t[1].ptr();

            madness::Tensor<double> rhoa_t (np);
            madness::Tensor<double> rhob_t (np);
            double * restrict densa = rhoa_t.ptr();
            double * restrict densb = rhob_t.ptr();

            for (long i=0; i<np; i++) {
                densa[i] = munge(rhoa[i]);
                densb[i] = munge(rhob[i]);
            }
            
           rho.push_back(rhoa_t);
           rho.push_back(rhob_t);
        }
        else if (is_gga()) {
           MADNESS_ASSERT(t.size() == 5);
//vama3            const double * restrict rhoa  = t[0].ptr();
//vama3            const double * restrict rhob  = t[1].ptr();
//vama3
//vama3            const double * restrict sigaa = t[2].ptr();
//vama3            const double * restrict sigab = t[3].ptr();
//vama3            const double * restrict sigbb = t[4].ptr();
//vama3
//vama3            const double * restrict drhoax = t[5].ptr();
//vama3            const double * restrict drhoay = t[6].ptr();
//vama3            const double * restrict drhoaz = t[7].ptr();
//vama3
//vama3            const double * restrict drhobx = t[8].ptr();
//vama3            const double * restrict drhoby = t[9].ptr();
//vama3            const double * restrict drhobz = t[10].ptr();
//vama3
//vama3            rho   = madness::Tensor<double>(np*2L);
//vama3            sigma = madness::Tensor<double>(np*3L);
//vama3            delrho  = madness::Tensor<double>(np*6L);
//vama3
//vama3            double * restrict ddens = delrho.ptr();
//vama3            double * restrict dens = rho.ptr();
//vama3            double * restrict sig  = sigma.ptr();
//vama3            for (long i=0; i<np; i++) {
//vama3                double ra=rhoa[i], rb=rhob[i], saa=sigaa[i], sab=sigab[i], sbb=sigbb[i];
//vama3
//vama3                double dax=drhoax[i], day=drhoay[i], daz=drhoaz[i];
//vama3                double dbx=drhobx[i], dby=drhoby[i], dbz=drhobz[i];
//vama3            //    munge5(ra, rb, saa, sab, sbb);
//vama3                munge5_der(ra, rb, saa, sab, sbb, dax, day, daz, dbx, dby, dbz);
//vama3                dens[2*i  ] = ra;
//vama3                dens[2*i+1] = rb;
//vama3
//vama3                sig[3*i  ] = saa;
//vama3                sig[3*i+1] = sab;
//vama3                sig[3*i+2] = sbb;
//vama3
//vama3                ddens[6*i  ] = dax;
//vama3                ddens[6*i+1] = day;
//vama3                ddens[6*i+2] = daz;
//vama3
//vama3                ddens[6*i+3] = dbx;
//vama3                ddens[6*i+4] = dby;
//vama3                ddens[6*i+5] = dbz;
//vama3            }
        }
        else {
            throw "not yet";
       }
    }
    else {
        if (is_lda()) {
            MADNESS_ASSERT(t.size() == 1);

            const double * restrict rhoa = t[0].ptr();

            madness::Tensor<double> rhoa_t (np);
            double * restrict densa = rhoa_t.ptr();

            for (long i=0; i<np; i++) {
                densa[i] = munge(2*rhoa[i]);
            }
            
           rho.push_back(rhoa_t);
        }
        else if (is_gga()) {
//            MADNESS_ASSERT(t.size() == 2);
            const double * restrict rhoa = t[0].ptr();
            const double * restrict sigaa = t[1].ptr();
            const double * restrict drhoax = t[2].ptr();
            const double * restrict drhoay = t[3].ptr();
            const double * restrict drhoaz = t[4].ptr();

            madness::Tensor<double> rhoa_t (np);
            madness::Tensor<double> sigmaa_t (np);
            madness::Tensor<double> delrhox_t (np);
            madness::Tensor<double> delrhoy_t (np);
            madness::Tensor<double> delrhoz_t (np);

            double * restrict densa = rhoa_t.ptr();
            double * restrict sig_aa = sigmaa_t.ptr();
            double * restrict ddensax = delrhox_t.ptr();
            double * restrict ddensay = delrhoy_t.ptr();
            double * restrict ddensaz = delrhoz_t.ptr();

            for (long i=0; i<np; i++) {
                double ra=2.0*rhoa[i], saa=4.0*sigaa[i];
                double dax=2.0*drhoax[i];
                double day=2.0*drhoay[i];
                double daz=2.0*drhoaz[i];
                munge_der(ra, saa, dax, day, daz);
                densa[i] = ra;
                sig_aa[i] = saa;
                ddensax[i] = dax;
                ddensay[i] = day;
                ddensaz[i] = daz;
            }
           rho.push_back(rhoa_t);
           sigma.push_back(sigmaa_t);
           delrho.push_back(delrhox_t);
           delrho.push_back(delrhoy_t);
           delrho.push_back(delrhoz_t);
        }
        else {
            throw "not yet";
        }
    }
}
madness::Tensor<double> XCfunctional::exc(const std::vector< madness::Tensor<double> >& t, const int ispin) const
{

    std::vector<madness::Tensor<double> > rho;
    std::vector<madness::Tensor<double> > sigma;
    std::vector<madness::Tensor<double> > delrho;
//    madness::Tensor<double> rho;
//    madness::Tensor<double> sigma;
//    madness::Tensor<double> delrho;

    make_xc_args(t, rho, sigma, delrho);

    const int  np = t[0].size();

    madness::Tensor<double> result(3L, t[0].dims(), false);
    double * restrict res = result.ptr();

//    madness::Tensor<double> tmp(3L, t[0].dims(), false);
//    double * restrict work = tmp.ptr();

    int deriv = 0;

//vama1    madness::Tensor<double> sigmaaa1(np);
//vama1    madness::Tensor<double> sigmabb1(np);
//vama1    madness::Tensor<double> sigmaab1(np);
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

    madness::Tensor<double> vrhoa(np);
    madness::Tensor<double> vrhob(np);

    make_xc_args(t, rho, sigma, delrho);

    if (spin_polarized) {
        MADNESS_ASSERT(t.size() == 2);

        const double * restrict rhoa = rho[0].ptr();
        const double * restrict rhob = rho[1].ptr();

        if (lda) {

            madness::Tensor<double> sigmaaa1(np);
            madness::Tensor<double> sigmabb1(np);
            madness::Tensor<double> sigmaab1(np);

    
            madness::Tensor<double> qx(3L, t[0].dims(), false);
            madness::Tensor<double> qc(3L, t[0].dims(), false);
    
            uks_x_lda_(&deriv, &np, rhoa, rhob,
                        sigmaaa1.ptr(),sigmabb1.ptr(),sigmaab1.ptr(),
                        qx.ptr(),vrhoa.ptr(),vrhob.ptr(),vsigmaaa.ptr(),vsigmabb.ptr(),vsigmaab.ptr(),
                        v2rhoa2.ptr(),v2rhob2.ptr(),v2rhoab.ptr(),
                        v2rhoasigmaaa.ptr(),v2rhoasigmaab.ptr(),v2rhoasigmabb.ptr(),
                        v2rhobsigmabb.ptr(),v2rhobsigmaab.ptr(),v2rhobsigmaaa.ptr(),
                        v2sigmaaa2.ptr(),v2sigmaaaab.ptr(),v2sigmaaabb.ptr(),
                        v2sigmaab2.ptr(),v2sigmaabbb.ptr(),v2sigmabb2.ptr());
            uks_c_vwn5_(&deriv, &np, rhoa, rhob,
                        sigmaaa1.ptr(),sigmabb1.ptr(),sigmaab1.ptr(),
                        qc.ptr(),vrhoa.ptr(),vrhob.ptr(),vsigmaaa.ptr(),vsigmabb.ptr(),vsigmaab.ptr(),
                        v2rhoa2.ptr(),v2rhob2.ptr(),v2rhoab.ptr(),
                        v2rhoasigmaaa.ptr(),v2rhoasigmaab.ptr(),v2rhoasigmabb.ptr(),
                        v2rhobsigmabb.ptr(),v2rhobsigmaab.ptr(),v2rhobsigmaaa.ptr(),
                        v2sigmaaa2.ptr(),v2sigmaaaab.ptr(),v2sigmaaabb.ptr(),
                        v2sigmaab2.ptr(),v2sigmaabbb.ptr(),v2sigmabb2.ptr());
            result = qx + qc;
        } 
        else if(gga) {
            MADNESS_ASSERT(t.size() == 1);
//vama3
//vama3
//vama3
//vama3
//vama3
        } 
    }
    else {
    //    MADNESS_ASSERT(t.size() == 1);

        const double * restrict rhoa = rho[0].ptr();

        madness::Tensor<double> vrhoa(np);
        madness::Tensor<double> vsigmaaa(np);
        madness::Tensor<double> v2rhoa2(np);
        madness::Tensor<double> v2rhoasigmaaa(np);
        madness::Tensor<double> v2sigmaaa2(np);

        madness::Tensor<double> qx(3L, t[0].dims(), false);
        madness::Tensor<double> qc(3L, t[0].dims(), false);

        if (lda) {
            madness::Tensor<double> sigmaaa(np);
            rks_x_lda_(&deriv, &np, rhoa, sigmaaa.ptr(), qx.ptr(), 
                        vrhoa.ptr(), vsigmaaa.ptr(), v2rhoa2.ptr(), v2rhoasigmaaa.ptr(), v2sigmaaa2.ptr());
            rks_c_vwn5_(&deriv, &np, rhoa, sigmaaa.ptr(), qc.ptr(), 
                        vrhoa.ptr(), vsigmaaa.ptr(), v2rhoa2.ptr(), v2rhoasigmaaa.ptr(), v2sigmaaa2.ptr());
            result = qx + qc;
        } 
        else if(gga) {
//vama3            MADNESS_ASSERT(t.size() == 1);
            const double * restrict sig = sigma[0].ptr();

//            rks_x_pbe_(&deriv, &np, rhoa, sig, qx.ptr(), 
//                        vrhoa.ptr(), vsigmaaa.ptr(), v2rhoa2.ptr(), v2rhoasigmaaa.ptr(), v2sigmaaa2.ptr());
//            rks_c_pbe_(&deriv, &np, rhoa, sig, qc.ptr(), 
//                        vrhoa.ptr(), vsigmaaa.ptr(), v2rhoa2.ptr(), v2rhoasigmaaa.ptr(), v2sigmaaa2.ptr());
            rks_x_b88_(&deriv, &np, rhoa, sig, qx.ptr(), 
                        vrhoa.ptr(), vsigmaaa.ptr(), v2rhoa2.ptr(), v2rhoasigmaaa.ptr(), v2sigmaaa2.ptr());
            rks_c_lyp_(&deriv, &np, rhoa, sig, qc.ptr(), 
                        vrhoa.ptr(), vsigmaaa.ptr(), v2rhoa2.ptr(), v2rhoasigmaaa.ptr(), v2sigmaaa2.ptr());
            result = qx + qc;
        } 
    }
    return result;
}

madness::Tensor<double> XCfunctional::vxc(const std::vector< madness::Tensor<double> >& t, const int ispin, const int what) const
{
    //MADNESS_ASSERT(what == 0);

    std::vector<madness::Tensor<double> > rho;
    std::vector<madness::Tensor<double> > sigma;
    std::vector<madness::Tensor<double> > delrho;

    make_xc_args(t, rho, sigma, delrho);

    const int  np = t[0].size();

    madness::Tensor<double> result(3L, t[0].dims(), false);
    double * restrict res = result.ptr();


    int deriv = 1;
    
    int nvsig, nvrho;
    if (spin_polarized) {
        nvrho = 2;
        nvsig = 3;
    }
    else {
        nvrho = 1;
        nvsig = 1;
    }

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


    if (spin_polarized) {
        MADNESS_ASSERT(t.size() == 2);

        madness::Tensor<double> vrhoax(3L, t[0].dims(), false);
        madness::Tensor<double> vrhobx(3L, t[0].dims(), false);
        madness::Tensor<double> vrhoac(3L, t[0].dims(), false);
        madness::Tensor<double> vrhobc(3L, t[0].dims(), false);

        const double * restrict rhoa = rho[0].ptr();
        const double * restrict rhob = rho[1].ptr();

        if (lda) {

//           std::cout << " hola out\n" ;
            madness::Tensor<double> sigmaaa1(np);
            madness::Tensor<double> sigmabb1(np);
            madness::Tensor<double> sigmaab1(np);

            madness::Tensor<double> qx(3L, t[0].dims(), false);
            madness::Tensor<double> qc(3L, t[0].dims(), false);
    
            uks_x_lda_(&deriv, &np, rhoa, rhob,
                        sigmaaa1.ptr(),sigmabb1.ptr(),sigmaab1.ptr(),
                        qx.ptr(),vrhoax.ptr(),vrhobx.ptr(),vsigmaaa.ptr(),vsigmabb.ptr(),vsigmaab.ptr(),
                        v2rhoa2.ptr(),v2rhob2.ptr(),v2rhoab.ptr(),
                        v2rhoasigmaaa.ptr(),v2rhoasigmaab.ptr(),v2rhoasigmabb.ptr(),
                        v2rhobsigmabb.ptr(),v2rhobsigmaab.ptr(),v2rhobsigmaaa.ptr(),
                        v2sigmaaa2.ptr(),v2sigmaaaab.ptr(),v2sigmaaabb.ptr(),
                        v2sigmaab2.ptr(),v2sigmaabbb.ptr(),v2sigmabb2.ptr());
            uks_c_vwn5_(&deriv, &np, rhoa, rhob,
                        sigmaaa1.ptr(),sigmabb1.ptr(),sigmaab1.ptr(),
                        qc.ptr(),vrhoac.ptr(),vrhobc.ptr(),vsigmaaa.ptr(),vsigmabb.ptr(),vsigmaab.ptr(),
                        v2rhoa2.ptr(),v2rhob2.ptr(),v2rhoab.ptr(),
                        v2rhoasigmaaa.ptr(),v2rhoasigmaab.ptr(),v2rhoasigmabb.ptr(),
                        v2rhobsigmabb.ptr(),v2rhobsigmaab.ptr(),v2rhobsigmaaa.ptr(),
                        v2sigmaaa2.ptr(),v2sigmaaaab.ptr(),v2sigmaaabb.ptr(),
                        v2sigmaab2.ptr(),v2sigmaabbb.ptr(),v2sigmabb2.ptr());
    
            if (ispin == 0) {
               result = vrhoax + vrhoac ;
            }
            else {
               result = vrhobx + vrhobc ;
            }
        }
        else if(gga) {
            MADNESS_ASSERT(t.size() == 1);
//vama3
//vama3
//vama3
//vama3
        }
    }
    else {
//        MADNESS_ASSERT(t.size() == 1);

        madness::Tensor<double> v2rhoa2(np);
        madness::Tensor<double> v2rhoasigmaaa(np);
        madness::Tensor<double> v2sigmaaa2(np);

        madness::Tensor<double> qx(3L, t[0].dims(), false);
        madness::Tensor<double> qc(3L, t[0].dims(), false);

        madness::Tensor<double> vrhoax(3L, t[0].dims(), false);
        madness::Tensor<double> vrhoac(3L, t[0].dims(), false);

        const double * restrict rhoa = rho[0].ptr();

        if (lda) {
            madness::Tensor<double> sigmaaa(np);
            madness::Tensor<double> vsigmaaa(np);

            rks_x_lda_( &deriv, &np, rhoa, sigmaaa.ptr(), qx.ptr(), 
                        vrhoax.ptr(), vsigmaaa.ptr(), v2rhoa2.ptr(), v2rhoasigmaaa.ptr(), v2sigmaaa2.ptr());
    
            rks_c_vwn5_( &deriv, &np, rhoa, sigmaaa.ptr(), qc.ptr(), 
                        vrhoac.ptr(), vsigmaaa.ptr(), v2rhoa2.ptr(), v2rhoasigmaaa.ptr(), v2sigmaaa2.ptr());
    
            result = vrhoax + vrhoac;
        }
        else if(gga) {

//            madness::Tensor<double> vsigmaaa(nvsig*np);
            madness::Tensor<double> vsigmaaax(3L, t[0].dims(), false);
            madness::Tensor<double> vsigmaaac(3L, t[0].dims(), false);

//            double * restrict vrax = vrhoax.ptr();
//            double * restrict vrac = vrhoac.ptr();

            const double * restrict sig = sigma[0].ptr();

//            rks_x_pbe_(&deriv, &np, rhoa, sig, qx.ptr(), 
//                        vrhoax.ptr(), vsigmaaax.ptr(), v2rhoa2.ptr(), v2rhoasigmaaa.ptr(), v2sigmaaa2.ptr());
//            rks_c_pbe_(&deriv, &np, rhoa, sig, qc.ptr(), 
//                        vrhoac.ptr(), vsigmaaac.ptr(), v2rhoa2.ptr(), v2rhoasigmaaa.ptr(), v2sigmaaa2.ptr());
            rks_x_b88_(&deriv, &np, rhoa, sig, qx.ptr(), 
                        vrhoax.ptr(), vsigmaaax.ptr(), v2rhoa2.ptr(), v2rhoasigmaaa.ptr(), v2sigmaaa2.ptr());
            rks_c_lyp_(&deriv, &np, rhoa, sig, qc.ptr(), 
                        vrhoac.ptr(), vsigmaaac.ptr(), v2rhoa2.ptr(), v2rhoasigmaaa.ptr(), v2sigmaaa2.ptr());

            double * restrict vsax = vsigmaaax.ptr();
            double * restrict vsac = vsigmaaac.ptr();
            const double * restrict delx = delrho[0].ptr();
            const double * restrict dely = delrho[1].ptr();
            const double * restrict delz = delrho[2].ptr();


            if (what == 0) {
                  result = vrhoax + vrhoac;
            }
            else{
                  double * restrict res = result.ptr();
                  for (long j=0; j<np; j++) {
                      res[j] = (vsax[j]+vsac[j])*(delx[j]+dely[j]+delz[j])*2.0;
                      if (isnan_x(res[j])) throw "ouch";
                  }
            }
        } 
    }
    return result;
}


#endif
