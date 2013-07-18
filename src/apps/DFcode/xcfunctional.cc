#include <madness_config.h>

#ifndef MADNESS_HAS_LIBXC

//#ifndef HAS_LIBMADXC

#include <DFcode/xcfunctional.h>
#include <tensor/tensor.h>
#include <sstream>

#include <libMADxc.h>


#define  XC_LDA_X               1  /* Exchange                                                   */
#define  XC_LDA_C_RPA           3  /* Random Phase Approximation                                 */
#define  XC_LDA_C_VWN           7  /* Vosko, Wilk, & Nussair                                     */
#define  XC_LDA_C_PZ            9  /* Perdew & Zunger                                            */
#define  XC_LDA_C_PW           12  /* Perdew & Wang                                              */
#define  XC_GGA_X_PBE         101  /* Perdew, Burke & Ernzerhof exchange                         */
#define  XC_GGA_X_B88         106  /* Becke 88                                                   */
#define  XC_GGA_X_PW91        109  /* Perdew & Wang 91                                           */
#define  XC_GGA_X_FT97_B      115  /* Filatov & Thiel 97 (version B)                             */
#define  XC_GGA_C_PBE         130  /* Perdew, Burke & Ernzerhof correlation                      */
#define  XC_GGA_C_LYP         131  /* Lee, Yang & Parr                                           */
#define  XC_GGA_C_P86         132  /* Perdew 86                                                  */
#define  XC_GGA_C_PW91        134  /* Perdew & Wang 91                                           */
#define  XC_GGA_XC_HCTH_93    161  /* HCTH functional fitted to  93 molecules                    */
#define  XC_GGA_XC_HCTH_120   162  /* HCTH functional fitted to 120 molecules                    */
#define  XC_GGA_XC_HCTH_147   163  /* HCTH functional fitted to 147 molecules                    */
#define  XC_GGA_XC_HCTH_407   164  /* HCTH functional fitted to 147 molecules                    */
#define  XC_GGA_XC_EDF1       165  /* Empirical functionals from Adamson, Gill, and Pople        */
#define  XC_GGA_XC_B97        167  /* Becke 97                                                   */
#define  XC_GGA_XC_B97_1      168  /* Becke 97-1                                                 */
#define  XC_GGA_XC_B97_2      169  /* Becke 97-2                                                 */
#define  XC_HYB_GGA_XC_B3LYP  402  /* The (in)famous B3LYP                                       */
#define  LDA_FUNC               0
#define  GGA_FUNC               1  
#define  MGGA_FUNC              2  




struct xc_name_map {
    const std::string name;
    const int id;
    const int typefunc;
};


static xc_name_map map[] = {
    {"LDA_X", XC_LDA_X, LDA_FUNC},
    {"LDA_C_VWN", XC_LDA_C_VWN, LDA_FUNC},
    {"LDA_C_RPA", XC_LDA_C_RPA, LDA_FUNC},
    {"LDA_C_PZ", XC_LDA_C_PZ, LDA_FUNC},
    {"LDA_C_PW", XC_LDA_C_PW, LDA_FUNC},
    {"HYB_GGA_XC_B3LYP", XC_HYB_GGA_XC_B3LYP, GGA_FUNC},
    {"GGA_XC_HCTH_93", XC_GGA_XC_HCTH_93, GGA_FUNC},
    {"GGA_XC_HCTH_120", XC_GGA_XC_HCTH_120, GGA_FUNC},
    {"GGA_XC_HCTH_147", XC_GGA_XC_HCTH_147, GGA_FUNC},
    {"GGA_XC_HCTH_407", XC_GGA_XC_HCTH_407, GGA_FUNC},
    {"GGA_XC_EDF1", XC_GGA_XC_EDF1, GGA_FUNC},
    {"GGA_XC_B97_1", XC_GGA_XC_B97_1, GGA_FUNC},
    {"GGA_XC_B97_2", XC_GGA_XC_B97_2, GGA_FUNC},
    {"GGA_XC_B97", XC_GGA_XC_B97, GGA_FUNC},
    {"GGA_X_PW91", XC_GGA_X_PW91, GGA_FUNC},
    {"GGA_X_PBE", XC_GGA_X_PBE, GGA_FUNC},
    {"GGA_X_FT97_B", XC_GGA_X_FT97_B, GGA_FUNC},
    {"GGA_X_B88", XC_GGA_X_B88, GGA_FUNC},
    {"GGA_X_PW91", XC_GGA_X_PW91, GGA_FUNC},
    {"GGA_C_PBE", XC_GGA_C_PBE, GGA_FUNC},
    {"GGA_C_P86", XC_GGA_C_P86, GGA_FUNC},
    {"GGA_C_LYP", XC_GGA_C_LYP, GGA_FUNC},

    {"", -1,0} // MUST BE LAST
};

//vama3static int lookup_name(const std::string& name) {
//vama3    const xc_name_map* map = ::map;
//vama3    while (map->id > 0) {
//vama3        if (name == map->name) return map->id;
//vama3        map++;
//vama3    }
//vama3    return -1;
//vama3}
//vama3
//vama3static xc_func_type* make_func(int id, bool polarized) {
//vama3    xc_func_type* func = new xc_func_type;
//vama3    int POLARIZED = polarized ? XC_POLARIZED : XC_UNPOLARIZED;
//vama3    MADNESS_ASSERT(xc_func_init(func, id, POLARIZED) == 0);
//vama3    return func;
//vama3}
//vama3
//vama3static xc_func_type* lookup_func(const std::string& name, bool polarized) {
//vama3    int id = lookup_name(name);
//vama3    MADNESS_ASSERT(id > 0);
//vama3    return make_func(id, polarized);
//vama3}

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

//vama3void XCfunctional::initialize(const std::string& input_line, bool polarized)
//vama3{
//vama3    double factor;
//vama3    spin_polarized = polarized;
//vama3
//vama3    std::stringstream line(input_line);
//vama3    std::string name;
//vama3
//vama3    nderiv = 0;
//vama3    hf_coeff = 0.0;
//vama3    funcs.clear();
//vama3
//vama3    while (line >> name) {
//vama3        std::transform(name.begin(), name.end(), name.begin(), ::toupper);
//vama3        if (name == "LDA") {
//vama3            if (! (line >> factor)) factor = 1.0;
//vama3            funcs.push_back(std::make_pair(lookup_func("LDA_X",polarized),factor));
//vama3            funcs.push_back(std::make_pair(lookup_func("LDA_C_VWN",polarized),factor));
//vama3        }
//vama3        else if (name == "HF" || name == "HF_X") {
//vama3            if (! (line >> factor)) factor = 1.0;
//vama3            hf_coeff = factor;
//vama3        }
//vama3        else {
//vama3            if (! (line >> factor)) factor = 1.0;
//vama3            funcs.push_back(std::make_pair(lookup_func(name,polarized), factor));
//vama3        }
//vama3    }
//vama3
//vama3    for (unsigned int i=0; i<funcs.size(); i++) {
//vama3        if (funcs[i].first->info->family == XC_FAMILY_GGA) nderiv = std::max(nderiv,1);
//vama3        if (funcs[i].first->info->family == XC_FAMILY_HYB_GGA) nderiv = std::max(nderiv,1);
//vama3        if (funcs[i].first->info->family == XC_FAMILY_MGGA)nderiv = std::max(nderiv,2);
//vama3    }
//vama3}

XCfunctional::~XCfunctional() {}

//vama3bool XCfunctional::is_lda() const {
//vama3    return nderiv == 0;
//vama3}
//vama3
//vama3bool XCfunctional::is_gga() const {
//vama3    return nderiv == 1;
//vama3}
//vama3
//vama3bool XCfunctional::is_meta() const {
//vama3    return nderiv == 2;
//vama3}
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
//            MADNESS_ASSERT(t.size() == 2);
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
//           MADNESS_ASSERT(t.size() == 5);
            const double * restrict rhoa  = t[0].ptr();
            const double * restrict rhob  = t[1].ptr();

            const double * restrict sigaa = t[2].ptr();
            const double * restrict sigab = t[3].ptr();
            const double * restrict sigbb = t[4].ptr();

            const double * restrict drhoax = t[5].ptr();
            const double * restrict drhoay = t[6].ptr();
            const double * restrict drhoaz = t[7].ptr();

            const double * restrict drhobx = t[8].ptr();
            const double * restrict drhoby = t[9].ptr();
            const double * restrict drhobz = t[10].ptr();

            madness::Tensor<double> rhoa_t (np);
            madness::Tensor<double> rhob_t (np);
            madness::Tensor<double> sigmaa_t (np);
            madness::Tensor<double> sigmab_t (np);
            madness::Tensor<double> sigmbb_t (np);
            madness::Tensor<double> delrhoax_t (np);
            madness::Tensor<double> delrhoay_t (np);
            madness::Tensor<double> delrhoaz_t (np);

            madness::Tensor<double> delrhobx_t (np);
            madness::Tensor<double> delrhoby_t (np);
            madness::Tensor<double> delrhobz_t (np);

            double * restrict densa = rhoa_t.ptr();
            double * restrict densb = rhob_t.ptr();
            double * restrict sig_aa = sigmaa_t.ptr();
            double * restrict sig_ab = sigmab_t.ptr();
            double * restrict sig_bb = sigmbb_t.ptr();
            double * restrict ddensax = delrhoax_t.ptr();
            double * restrict ddensay = delrhoay_t.ptr();
            double * restrict ddensaz = delrhoaz_t.ptr();

            double * restrict ddensbx = delrhobx_t.ptr();
            double * restrict ddensby = delrhoby_t.ptr();
            double * restrict ddensbz = delrhobz_t.ptr();

            for (long i=0; i<np; i++) {
                double ra=rhoa[i], rb=rhob[i], saa=sigaa[i], sab=sigab[i], sbb=sigbb[i];

                double dax=drhoax[i], day=drhoay[i], daz=drhoaz[i];
                double dbx=drhobx[i], dby=drhoby[i], dbz=drhobz[i];

                munge5_der(ra, rb, saa, sab, sbb, dax, day, daz, dbx, dby, dbz);
                densa[i] = ra;
                densb[i] = rb;

                sig_aa[i] = saa;
                sig_ab[i] = sab;
                sig_bb[i] = sbb;

                ddensax[i] = dax;
                ddensay[i] = day;
                ddensaz[i] = daz;

                ddensbx[i] = dbx;
                ddensby[i] = dby;
                ddensbz[i] = dbz;
            }

            rho.push_back(rhoa_t);
            rho.push_back(rhob_t);
            sigma.push_back(sigmaa_t);
            sigma.push_back(sigmab_t);
            sigma.push_back(sigmbb_t);
            delrho.push_back(delrhoax_t);
            delrho.push_back(delrhoay_t);
            delrho.push_back(delrhoaz_t);
            delrho.push_back(delrhobx_t);
            delrho.push_back(delrhoby_t);
            delrho.push_back(delrhobz_t);
        }
        else {
            throw "not yet";
       }
    }
    else {
        if (is_lda()) {
//            MADNESS_ASSERT(t.size() == 1);

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

    madness::Tensor<double> result(3L, t[0].dims());
    //double * restrict res = result.ptr();

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
    //    MADNESS_ASSERT(t.size() == 2);

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
            //MADNESS_ASSERT(t.size() == 1);

            const double * restrict sigaa = sigma[0].ptr();
            const double * restrict sigab = sigma[1].ptr();
            const double * restrict sigbb = sigma[2].ptr();

            madness::Tensor<double> qx(3L, t[0].dims(), false);
            madness::Tensor<double> qc(3L, t[0].dims(), false);
    
//            uks_x_pbe_  (&deriv, &np, rhoa, rhob,
//                        sigaa, sigab, sigbb,
//                        qx.ptr(),vrhoa.ptr(),vrhob.ptr(),vsigmaaa.ptr(),vsigmabb.ptr(),vsigmaab.ptr(),
//                        v2rhoa2.ptr(),v2rhob2.ptr(),v2rhoab.ptr(),
//                        v2rhoasigmaaa.ptr(),v2rhoasigmaab.ptr(),v2rhoasigmabb.ptr(),
//                        v2rhobsigmabb.ptr(),v2rhobsigmaab.ptr(),v2rhobsigmaaa.ptr(),
//                        v2sigmaaa2.ptr(),v2sigmaaaab.ptr(),v2sigmaaabb.ptr(),
//                        v2sigmaab2.ptr(),v2sigmaabbb.ptr(),v2sigmabb2.ptr());
//            uks_c_pbe_  (&deriv, &np, rhoa, rhob,
//                        sigaa, sigab, sigbb,
//                        qc.ptr(),vrhoa.ptr(),vrhob.ptr(),vsigmaaa.ptr(),vsigmabb.ptr(),vsigmaab.ptr(),
//                        v2rhoa2.ptr(),v2rhob2.ptr(),v2rhoab.ptr(),
//                        v2rhoasigmaaa.ptr(),v2rhoasigmaab.ptr(),v2rhoasigmabb.ptr(),
//                        v2rhobsigmabb.ptr(),v2rhobsigmaab.ptr(),v2rhobsigmaaa.ptr(),
//                        v2sigmaaa2.ptr(),v2sigmaaaab.ptr(),v2sigmaaabb.ptr(),
//                        v2sigmaab2.ptr(),v2sigmaabbb.ptr(),v2sigmabb2.ptr());
            uks_x_b88_  (&deriv, &np, rhoa, rhob,
                        sigaa, sigab, sigbb,
                        qx.ptr(),vrhoa.ptr(),vrhob.ptr(),vsigmaaa.ptr(),vsigmabb.ptr(),vsigmaab.ptr(),
                        v2rhoa2.ptr(),v2rhob2.ptr(),v2rhoab.ptr(),
                        v2rhoasigmaaa.ptr(),v2rhoasigmaab.ptr(),v2rhoasigmabb.ptr(),
                        v2rhobsigmabb.ptr(),v2rhobsigmaab.ptr(),v2rhobsigmaaa.ptr(),
                        v2sigmaaa2.ptr(),v2sigmaaaab.ptr(),v2sigmaaabb.ptr(),
                        v2sigmaab2.ptr(),v2sigmaabbb.ptr(),v2sigmabb2.ptr());
            uks_c_lyp_  (&deriv, &np, rhoa, rhob,
                        sigaa, sigab, sigbb,
                        qc.ptr(),vrhoa.ptr(),vrhob.ptr(),vsigmaaa.ptr(),vsigmabb.ptr(),vsigmaab.ptr(),
                        v2rhoa2.ptr(),v2rhob2.ptr(),v2rhoab.ptr(),
                        v2rhoasigmaaa.ptr(),v2rhoasigmaab.ptr(),v2rhoasigmabb.ptr(),
                        v2rhobsigmabb.ptr(),v2rhobsigmaab.ptr(),v2rhobsigmaaa.ptr(),
                        v2sigmaaa2.ptr(),v2sigmaaaab.ptr(),v2sigmaaabb.ptr(),
                        v2sigmaab2.ptr(),v2sigmaabbb.ptr(),v2sigmabb2.ptr());
            result = qx + qc;
            //result = qx ;
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
            //result = qx ;
//            double * restrict work = qx.ptr();
//            double * restrict res = result.ptr();
//            for (long j=0; j<np; j++) {
//                res[j] += work[j]*rhoa[j];
//            }
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

    madness::Tensor<double> result(3L, t[0].dims());
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


        //   std::cout << " hola vcx\n" ;
    if (spin_polarized) {
       // MADNESS_ASSERT(t.size() == 2);

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
            //MADNESS_ASSERT(t.size() == 1);
            madness::Tensor<double> vsigmaaax(3L, t[0].dims(), false);
            madness::Tensor<double> vsigmaabx(3L, t[0].dims(), false);
            madness::Tensor<double> vsigmabbx(3L, t[0].dims(), false);

            madness::Tensor<double> vsigmaaac(3L, t[0].dims(), false);
            madness::Tensor<double> vsigmaabc(3L, t[0].dims(), false);
            madness::Tensor<double> vsigmabbc(3L, t[0].dims(), false);

            const double * restrict sigaa = sigma[0].ptr();
            const double * restrict sigab = sigma[1].ptr();
            const double * restrict sigbb = sigma[2].ptr();

            madness::Tensor<double> qx(3L, t[0].dims(), false);
            madness::Tensor<double> qc(3L, t[0].dims(), false);
    
//            uks_x_pbe_  (&deriv, &np, rhoa, rhob,
//                        sigaa, sigab, sigbb,
//                        qx.ptr(),vrhoa.ptr(),vrhob.ptr(),vsigmaaa.ptr(),vsigmabb.ptr(),vsigmaab.ptr(),
//                        v2rhoa2.ptr(),v2rhob2.ptr(),v2rhoab.ptr(),
//                        v2rhoasigmaaa.ptr(),v2rhoasigmaab.ptr(),v2rhoasigmabb.ptr(),
//                        v2rhobsigmabb.ptr(),v2rhobsigmaab.ptr(),v2rhobsigmaaa.ptr(),
//                        v2sigmaaa2.ptr(),v2sigmaaaab.ptr(),v2sigmaaabb.ptr(),
//                        v2sigmaab2.ptr(),v2sigmaabbb.ptr(),v2sigmabb2.ptr());
//            uks_c_pbe_  (&deriv, &np, rhoa, rhob,
//                        sigaa, sigab, sigbb,
//                        qc.ptr(),vrhoa.ptr(),vrhob.ptr(),vsigmaaa.ptr(),vsigmabb.ptr(),vsigmaab.ptr(),
//                        v2rhoa2.ptr(),v2rhob2.ptr(),v2rhoab.ptr(),
//                        v2rhoasigmaaa.ptr(),v2rhoasigmaab.ptr(),v2rhoasigmabb.ptr(),
//                        v2rhobsigmabb.ptr(),v2rhobsigmaab.ptr(),v2rhobsigmaaa.ptr(),
//                        v2sigmaaa2.ptr(),v2sigmaaaab.ptr(),v2sigmaaabb.ptr(),
//                        v2sigmaab2.ptr(),v2sigmaabbb.ptr(),v2sigmabb2.ptr());
            uks_x_b88_  (&deriv, &np, rhoa, rhob,
                        sigaa, sigab, sigbb,
                        qx.ptr(),vrhoax.ptr(),vrhobx.ptr(),vsigmaaax.ptr(),vsigmabbx.ptr(),vsigmaabx.ptr(),
                        v2rhoa2.ptr(),v2rhob2.ptr(),v2rhoab.ptr(),
                        v2rhoasigmaaa.ptr(),v2rhoasigmaab.ptr(),v2rhoasigmabb.ptr(),
                        v2rhobsigmabb.ptr(),v2rhobsigmaab.ptr(),v2rhobsigmaaa.ptr(),
                        v2sigmaaa2.ptr(),v2sigmaaaab.ptr(),v2sigmaaabb.ptr(),
                        v2sigmaab2.ptr(),v2sigmaabbb.ptr(),v2sigmabb2.ptr());
            uks_c_lyp_  (&deriv, &np, rhoa, rhob,
                        sigaa, sigab, sigbb,
                        qc.ptr(),vrhoac.ptr(),vrhobc.ptr(),vsigmaaac.ptr(),vsigmabbc.ptr(),vsigmaabc.ptr(),
                        v2rhoa2.ptr(),v2rhob2.ptr(),v2rhoab.ptr(),
                        v2rhoasigmaaa.ptr(),v2rhoasigmaab.ptr(),v2rhoasigmabb.ptr(),
                        v2rhobsigmabb.ptr(),v2rhobsigmaab.ptr(),v2rhobsigmaaa.ptr(),
                        v2sigmaaa2.ptr(),v2sigmaaaab.ptr(),v2sigmaaabb.ptr(),
                        v2sigmaab2.ptr(),v2sigmaabbb.ptr(),v2sigmabb2.ptr());

            double * restrict vsaax = vsigmaaax.ptr();
            double * restrict vsabx = vsigmaabx.ptr();
            double * restrict vsbbx = vsigmabbx.ptr();

            double * restrict vsaac = vsigmaaac.ptr();
            double * restrict vsabc = vsigmaabc.ptr();
            double * restrict vsbbc = vsigmabbc.ptr();

            const double * restrict delax = delrho[0].ptr();
            const double * restrict delay = delrho[1].ptr();
            const double * restrict delaz = delrho[2].ptr();
            const double * restrict delbx = delrho[0].ptr();
            const double * restrict delby = delrho[1].ptr();
            const double * restrict delbz = delrho[2].ptr();
            if (what == 0) {
                  if (ispin == 0){
                      result = vrhoax + vrhoac;
                      //result = vrhoax ;
                  }
                  else{
                  //     result = vrhobx ;
                  result = vrhobx +vrhobc;
                  }
            }
            else{
                  if (ispin == 0){
                      for (long j=0; j<np; j++) {
                          res[j] += (vsaax[j]+vsaac[j])*(delax[j]+delay[j]+delaz[j])*2.0;
                          res[j] += (vsabx[j]+vsabc[j])*(delbx[j]+delby[j]+delbz[j]);
                          //res[j] += (vsaax[j])*(delax[j]+delay[j]+delaz[j])*2.0;
                          //res[j] += (vsabx[j])*(delbx[j]+delby[j]+delbz[j]);
                          if (isnan_x(res[j])) throw "ouch";
                      }
                  }
                  else{
                      for (long j=0; j<np; j++) {
                          res[j] += (vsbbx[j]+vsbbc[j])*(delbx[j]+delby[j]+delbz[j])*2.0;
                          res[j] += (vsabx[j]+vsabc[j])*(delbx[j]+delby[j]+delbz[j]);
                          //res[j] += (vsbbx[j])*(delbx[j]+delby[j]+delbz[j])*2.0;
                          //res[j] += (vsabx[j])*(delbx[j]+delby[j]+delbz[j]);
                          if (isnan_x(res[j])) throw "ouch";
                      }
                  }
            }
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
                  //result = vrhoax ;
                  result = vrhoax + vrhoac;
            }
            else{
                  for (long j=0; j<np; j++) {
                      //res[j] += (vsax[j])*(delx[j]+dely[j]+delz[j])*2.0;
                      res[j] = (vsax[j]+vsac[j])*(delx[j]+dely[j]+delz[j])*2.0;
                      if (isnan_x(res[j])) throw "ouch";
                  }
            }
        } 
    }
    return result;
}


#endif
