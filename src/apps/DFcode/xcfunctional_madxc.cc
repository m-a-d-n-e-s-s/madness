#include <madness_config.h>


#include <DFcode/xcfunctional.h>
#include <tensor/tensor.h>
#include <sstream>

#include <libMADxc.h>


struct xc_name_map {
    const std::string name;
    const int id;
    const int typefunc;
};


static xc_name_map map[] = {
    {"LDA_X",             XC_LDA_X,             LDA_FUNC},
    {"LDA_C_RPA",         XC_LDA_C_RPA,         LDA_FUNC},
    {"LDA_C_VWN",         XC_LDA_C_VWN,         LDA_FUNC},
    {"LDA_C_PZ",          XC_LDA_C_PZ,          LDA_FUNC},
    {"LDA_C_PW",          XC_LDA_C_PW,          LDA_FUNC},
    {"GGA_X_PBE",         XC_GGA_X_PBE,         GGA_FUNC},
    {"GGA_X_B88",         XC_GGA_X_B88,         GGA_FUNC},
    {"GGA_X_PW91",        XC_GGA_X_PW91,        GGA_FUNC},
    {"GGA_X_FT97_B",      XC_GGA_X_FT97_B,      GGA_FUNC},
    {"GGA_C_PBE",         XC_GGA_C_PBE,         GGA_FUNC},
    {"GGA_C_LYP",         XC_GGA_C_LYP,         GGA_FUNC},
    {"GGA_C_P86",         XC_GGA_C_P86,         GGA_FUNC},
    {"GGA_C_PW91",        XC_GGA_C_PW91,        GGA_FUNC},
    {"GGA_XC_HCTH_93",    XC_GGA_XC_HCTH_93,    GGA_FUNC},
    {"GGA_XC_HCTH_120",   XC_GGA_XC_HCTH_120,   GGA_FUNC},
    {"GGA_XC_HCTH_147",   XC_GGA_XC_HCTH_147,   GGA_FUNC},
    {"GGA_XC_HCTH_407",   XC_GGA_XC_HCTH_407,   GGA_FUNC},
    {"GGA_XC_EDF1",       XC_GGA_XC_EDF1,       GGA_FUNC},
    {"GGA_XC_B97",        XC_GGA_XC_B97,        GGA_FUNC},
    {"GGA_XC_B97_1",      XC_GGA_XC_B97_1,      GGA_FUNC},
    {"GGA_XC_B97_2",      XC_GGA_XC_B97_2,      GGA_FUNC},
    {"GGA_XC_PBE",        XC_GGA_XC_PBE,        GGA_FUNC},
    {"HYB_GGA_XC_B3LYP",  XC_HYB_GGA_XC_B3LYP,  GGA_FUNC},

    {"", -1,0} // MUST BE LAST
};

static int lookup_name(const std::string& name) {
    const xc_name_map* map = ::map;
    while (map->id > 0) {
        if (name == map->name) {
           return map->id;
        }
        map++;
    }
    return -1;
}

static int lookup_type(const std::string& name) {
    const xc_name_map* map = ::map;
    while (map->id > 0) {
        if (name == map->name) {
           return map->typefunc;
        }
        map++;
    }
    return -1;
}

//vama4static xc_func_type* make_func(int id, bool polarized) {
static xc_func_type* make_func(int id, int type,bool polarized) {
    xc_func_type* func = new xc_func_type;
//vama4    int POLARIZED = polarized ? XC_POLARIZED : XC_UNPOLARIZED;
    //MADNESS_ASSERT(xc_func_init(func, id, POLARIZED) == 0);
    func->number = id;
    func->family = type;
    func->nspin = polarized;
    return func;
}

//vama4static xc_func_type* lookup_func(const std::string& name, bool polarized) {
static xc_func_type* lookup_func(const std::string& name, bool polarized) {
    int id = lookup_name(name);
    int type = lookup_type(name);
    MADNESS_ASSERT(id > 0);
//    MADNESS_ASSERT(type > 0);
    return make_func(id, type, polarized);
    //return type;
}

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

//vama4void XCfunctional::initialize(const std::string& input_line, bool polarized) 
//vama4{
//vama4    spin_polarized = polarized;
//vama4    
//vama4 
//vama4    std::stringstream s(input_line);
//vama4    std::string token;
//vama4    while (s >> token) {
//vama4        if (token == "lda") {
//vama4            lda = true;
//vama4        }
//vama4        if (token == "pbe") {
//vama4            gga = true;
//vama4        }
//vama4        else if (token == "hf") {
//vama4            hf_coeff = 1.0;
//vama4        }
//vama4    }
//vama4    hf_coeff=0.0;
//vama4    lda=false;
//vama4    gga=false;
//vama4
//vama4    gga=true;
//vama4}

void XCfunctional::initialize(const std::string& input_line, bool polarized)
{
    double factor;
    spin_polarized = polarized;

    std::stringstream line(input_line);
    std::string name;

    //i= 0;
    hf_coeff = 0.0;
    funcs.clear();
    nderiv = 0;

    std::cout << "hola\n";
    while (line >> name) {
        std::transform(name.begin(), name.end(), name.begin(), ::toupper);
        if (name == "LDA") {
            if (! (line >> factor)) factor = 1.0;
            funcs.push_back(std::make_pair(lookup_func("LDA_X",polarized),factor));
            funcs.push_back(std::make_pair(lookup_func("LDA_C_VWN",polarized),factor));
            std::cout << "factor = "<< factor<<"\n";
        }
        else if (name == "HF" || name == "HF_X") {
            if (! (line >> factor)) factor = 1.0;
            hf_coeff = factor;
            std::cout << "factor = "<< factor<<"\n";
        }
        else {
            if (! (line >> factor)) factor = 1.0;
            funcs.push_back(std::make_pair(lookup_func(name,polarized), factor));
            std::cout << "factor = "<< factor<<"\n";
        }
    }

    for (unsigned int i=0; i<funcs.size(); i++) {
//vama5            std::cout << "qqueue = "<< funcs[i].first->number<<"\n";
        if (funcs[i].first->family == GGA_FUNC) nderiv = std::max(nderiv,1);
//vama4        if (funcs[i].first->info->family == MGGA_FUNC)nderiv = std::max(nderiv,2);
            std::cout << "familgia = "<< funcs[i].first->family<<"\n";
    }
}

XCfunctional::~XCfunctional() {}

bool XCfunctional::is_lda() const {
    return nderiv == 0;
}

bool XCfunctional::is_gga() const {
    return nderiv == 1;
}

bool XCfunctional::is_meta() const {
    return nderiv == 2;
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

            for (unsigned int i=0; i<np; i++) {
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

            for (unsigned int i=0; i<np; i++) {
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

            for (unsigned int i=0; i<np; i++) {
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

            for (unsigned int i=0; i<np; i++) {
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

    //make_xc_args(t, rho, sigma, delrho);

    const int  np = t[0].size();

    madness::Tensor<double> result(3L, t[0].dims(),false);
    ///change
    double * restrict res = result.ptr();
    for (unsigned int j=0; j<np; j++) res[j] =  0.0 ;

//    madness::Tensor<double> tmp(3L, t[0].dims(), false);
//    double * restrict work = tmp.ptr();

    int deriv = 0;

    double * vsigmaaa;
    double * vsigmabb;
    double * vsigmaab;

    double * v2rhoa2;
    double * v2rhob2;
    double * v2rhoab;

    double * v2rhoasigmaaa;
    double * v2rhoasigmaab;
    double * v2rhoasigmabb;

    double * v2rhobsigmabb;
    double * v2rhobsigmaab;
    double * v2rhobsigmaaa;

    double * v2sigmaaa2;
    double * v2sigmaaaab;
    double * v2sigmaaabb;

    double * v2sigmaab2;
    double * v2sigmaabbb;
    double * v2sigmabb2;

    double * vrhoa;
    double * vrhob;

    make_xc_args(t, rho, sigma, delrho);

    madness::Tensor<double> qx(3L, t[0].dims(), false);
    //double * restrict qp = qx.ptr();
    //const double * restrict dens = rho.ptr();

    if (spin_polarized) {
    //    MADNESS_ASSERT(t.size() == 2);
        const double * restrict rhoa = rho[0].ptr();
        const double * restrict rhob = rho[1].ptr();

        for (unsigned int i=0; i<funcs.size(); i++) {
            switch(funcs[i].first->family) {
                case LDA_FUNC:
                {
                    double * sigmaaa;
                    double * sigmabb;
                    double * sigmaab;
        
                    xc_lda (funcs[i].first->number, spin_polarized,
                                &deriv, &np, rhoa, rhob,
                                sigmaaa, sigmabb, sigmaab,
                                qx.ptr(), vrhoa, vrhob, vsigmaaa, vsigmabb, vsigmaab,
                                v2rhoa2, v2rhob2, v2rhoab,
                                v2rhoasigmaaa, v2rhoasigmaab, v2rhoasigmabb,
                                v2rhobsigmabb, v2rhobsigmaab, v2rhobsigmaaa,
                                v2sigmaaa2, v2sigmaaaab, v2sigmaaabb,
                                v2sigmaab2, v2sigmaabbb, v2sigmabb2);
                    result +=  qx *funcs[i].second;
                 }
                    break;
                case GGA_FUNC:
                {
    //vama3            MADNESS_ASSERT(t.size() == 1);
                    const double * restrict sigmaaa = sigma[0].ptr();
                    const double * restrict sigmabb = sigma[1].ptr();
                    const double * restrict sigmaab = sigma[2].ptr();
    
                    xc_gga (funcs[i].first->number, spin_polarized,
                                &deriv, &np, rhoa, rhob,
                                sigmaaa, sigmabb, sigmaab,
                                qx.ptr(), vrhoa, vrhob, vsigmaaa, vsigmabb, vsigmaab,
                                v2rhoa2, v2rhob2, v2rhoab,
                                v2rhoasigmaaa, v2rhoasigmaab, v2rhoasigmabb,
                                v2rhobsigmabb, v2rhobsigmaab, v2rhobsigmaaa,
                                v2sigmaaa2, v2sigmaaaab, v2sigmaaabb,
                                v2sigmaab2, v2sigmaabbb, v2sigmabb2);
    
                    result +=  qx *funcs[i].second;
                 }
                    break;
                default:
                    throw "UGH!";
             }
        }
    }
    else {
    //    MADNESS_ASSERT(t.size() == 1);

        const double * restrict rhoa = rho[0].ptr();
        const double * restrict rhob;

        double * vrhoa;
        double * v2rhoa2;
        double * v2rhoasigmaaa;
        double * v2sigmaaa2;

    for (unsigned int i=0; i<funcs.size(); i++) {
        switch(funcs[i].first->family) {
            case LDA_FUNC:
                {
                     double * sigmaaa;
                     double * sigmabb;
                     double * sigmaab;
         
                     xc_lda (funcs[i].first->number, spin_polarized,
                                 &deriv, &np, rhoa, rhob,
                                 sigmaaa, sigmabb, sigmaab,
                                 qx.ptr(), vrhoa, vrhob, vsigmaaa, vsigmabb, vsigmaab,
                                 v2rhoa2, v2rhob2, v2rhoab,
                                 v2rhoasigmaaa, v2rhoasigmaab, v2rhoasigmabb,
                                 v2rhobsigmabb, v2rhobsigmaab, v2rhobsigmaaa,
                                 v2sigmaaa2, v2sigmaaaab, v2sigmaaabb,
                                 v2sigmaab2, v2sigmaabbb, v2sigmabb2);
                     result +=  qx *funcs[i].second;
                }
                break;
    
            case GGA_FUNC:
                {
     //vama3            MADNESS_ASSERT(t.size() == 1);
                     const double * restrict sigmaaa = sigma[0].ptr();
                     double * sigmabb;
                     double * sigmaab;
     
     
                     xc_gga (funcs[i].first->number, spin_polarized,
                                 &deriv, &np, rhoa, rhob,
                                 sigmaaa, sigmabb, sigmaab,
                                 qx.ptr(), vrhoa, vrhob, vsigmaaa, vsigmabb, vsigmaab,
                                 v2rhoa2, v2rhob2, v2rhoab,
                                 v2rhoasigmaaa, v2rhoasigmaab, v2rhoasigmabb,
                                 v2rhobsigmabb, v2rhobsigmaab, v2rhobsigmaaa,
                                 v2sigmaaa2, v2sigmaaaab, v2sigmaaabb,
                                 v2sigmaab2, v2sigmaabbb, v2sigmabb2);
     
                     result +=  qx *funcs[i].second;
                }
                break;
    
            default:
                throw "UGH!";
        }
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

    for (unsigned int j=0; j<np; j++) res[j] =  0.0 ;

    int deriv = 1;

    double * v2rhoa2;
    double * v2rhob2;
    double * v2rhoab;

    double * v2rhoasigmaaa;
    double * v2rhoasigmaab;
    double * v2rhoasigmabb;

    double * v2rhobsigmabb;
    double * v2rhobsigmaab;
    double * v2rhobsigmaaa;

    double * v2sigmaaa2;
    double * v2sigmaaaab;
    double * v2sigmaaabb;

    double * v2sigmaab2;
    double * v2sigmaabbb;
    double * v2sigmabb2;


        madness::Tensor<double> qx(3L, t[0].dims(), false);
    if (spin_polarized) {
       // MADNESS_ASSERT(t.size() == 2);

        madness::Tensor<double> vrhoa(3L, t[0].dims(), false);
        madness::Tensor<double> vrhob(3L, t[0].dims(), false);

        const double * restrict rhoa = rho[0].ptr();
        const double * restrict rhob = rho[1].ptr();

        for (unsigned int i=0; i<funcs.size(); i++) {
            switch(funcs[i].first->family) {
                case LDA_FUNC:
                     {
                     double * sigmaaa, * sigmabb, * sigmaab;
         
                     double * vsigmaaa, * vsigmabb, * vsigmaab;
         
                     xc_lda (funcs[i].first->number, spin_polarized,
                                 &deriv, &np, rhoa, rhob,
                                 sigmaaa, sigmabb, sigmaab,
                                 qx.ptr(), vrhoa.ptr(), vrhob.ptr(), vsigmaaa, vsigmabb, vsigmaab,
                                 v2rhoa2, v2rhob2, v2rhoab,
                                 v2rhoasigmaaa, v2rhoasigmaab, v2rhoasigmabb,
                                 v2rhobsigmabb, v2rhobsigmaab, v2rhobsigmaaa,
                                 v2sigmaaa2, v2sigmaaaab, v2sigmaaabb,
                                 v2sigmaab2, v2sigmaabbb, v2sigmabb2);

                         if (ispin == 0) {
                            result += vrhoa*funcs[i].second;
                         }
                         else {
                            result += vrhob*funcs[i].second;
                         }
                    }
                    break;
                case GGA_FUNC:
                {
    //vama3            MADNESS_ASSERT(t.size() == 1);
                    madness::Tensor<double> vsigmaaa(3L, t[0].dims(), false);
                    madness::Tensor<double> vsigmaab(3L, t[0].dims(), false);
                    madness::Tensor<double> vsigmabb(3L, t[0].dims(), false);
        
                    const double * restrict sigmaaa = sigma[0].ptr();
                    const double * restrict sigmabb = sigma[1].ptr();
                    const double * restrict sigmaab = sigma[2].ptr();
    
                    xc_gga (funcs[i].first->number, spin_polarized,
                                &deriv, &np, rhoa, rhob,
                                sigmaaa, sigmabb, sigmaab,
                                qx.ptr(), vrhoa.ptr(), vrhob.ptr(), 
                                vsigmaaa.ptr(), vsigmaab.ptr(), vsigmabb.ptr(),
                                v2rhoa2, v2rhob2, v2rhoab,
                                v2rhoasigmaaa, v2rhoasigmaab, v2rhoasigmabb,
                                v2rhobsigmabb, v2rhobsigmaab, v2rhobsigmaaa,
                                v2sigmaaa2, v2sigmaaaab, v2sigmaaabb,
                                v2sigmaab2, v2sigmaabbb, v2sigmabb2);
    
                 ///   result +=  qx *funcs[i].second;

        //vama5            const double * restrict delax = delrho[0].ptr();
        //vama5            const double * restrict delay = delrho[1].ptr();
        //vama5            const double * restrict delaz = delrho[2].ptr();
        //vama5            const double * restrict delbx = delrho[3].ptr();
        //vama5            const double * restrict delby = delrho[4].ptr();
        //vama5            const double * restrict delbz = delrho[5].ptr();

                    double * restrict vsaaa = vsigmaaa.ptr();
                    double * restrict vsaab = vsigmaab.ptr();
                    double * restrict vsabb = vsigmabb.ptr();

                    if (what == 0) {
                          if (ispin == 0){
                            result += vrhoa*funcs[i].second;
                         }
                         else {
                            result += vrhob*funcs[i].second;
                         }
                    }
                    else {
                          if (ispin == 0){
                                  for (unsigned  int j=0; j<np; j++) {
                                          res[j] += vsaaa[j]*2.0;
                                          res[j] += vsaab[j];
                                          res[j] *= funcs[i].second;
                                          if (isnan_x(res[j])) throw "ouch";
                                  }
                          }
                          else {
                                  for (unsigned int j=0; j<np; j++) {
                                          res[j] += vsabb[j]*2.0;
                                          res[j] += vsaab[j];
                                          res[j] *= funcs[i].second;
                                          if (isnan_x(res[j])) throw "ouch";
                                  }
                          }
                    }
                 }
                    break;
        
                default:
                    throw "UGH!";
             }
        }

    }
    else {
//        MADNESS_ASSERT(t.size() == 1);

        madness::Tensor<double> vrhoa(3L, t[0].dims(), false);

        const double * restrict rhoa = rho[0].ptr();
        const double * rhob;
        double * vrhob;

        for (unsigned int i=0; i<funcs.size(); i++) {
            switch(funcs[i].first->family) {
                case LDA_FUNC:
                {
                    double * sigmaaa, * sigmabb, * sigmaab;
                    double * vsigmaaa, * vsigmabb, * vsigmaab;
        
                    xc_lda (funcs[i].first->number, spin_polarized,
                                &deriv, &np, rhoa, rhob,
                                sigmaaa, sigmabb, sigmaab,
                                qx.ptr(), vrhoa.ptr(), vrhob, vsigmaaa, vsigmabb, vsigmaab,
                                v2rhoa2, v2rhob2, v2rhoab,
                                v2rhoasigmaaa, v2rhoasigmaab, v2rhoasigmabb,
                                v2rhobsigmabb, v2rhobsigmaab, v2rhobsigmaaa,
                                v2sigmaaa2, v2sigmaaaab, v2sigmaaabb,
                                v2sigmaab2, v2sigmaabbb, v2sigmabb2);
                    result += vrhoa*funcs[i].second;
                }
                    break;
                case GGA_FUNC:
                {
                    madness::Tensor<double> vsigmaaa(3L, t[0].dims(), false);
        
                    const double * restrict sigmaaa = sigma[0].ptr();
                    double * sigmabb;
                    double * sigmaab;
               
                    double * vsigmabb;
                    double * vsigmaab;

                    xc_gga (funcs[i].first->number, spin_polarized,
                                &deriv, &np, rhoa, rhob,
                                sigmaaa, sigmabb, sigmaab,
                                qx.ptr(), vrhoa.ptr(), vrhob, vsigmaaa.ptr(), vsigmabb, vsigmaab,
                                v2rhoa2, v2rhob2, v2rhoab,
                                v2rhoasigmaaa, v2rhoasigmaab, v2rhoasigmabb,
                                v2rhobsigmabb, v2rhobsigmaab, v2rhobsigmaaa,
                                v2sigmaaa2, v2sigmaaaab, v2sigmaaabb,
                                v2sigmaab2, v2sigmaabbb, v2sigmabb2);

                    if (what == 0) {
                          result += vrhoa*funcs[i].second*1.;
                    }
                    else {
        
                          double * restrict vsaaa = vsigmaaa.ptr();
        //                  const double * restrict delx = delrho[0].ptr();
        //                  const double * restrict dely = delrho[1].ptr();
        //                  const double * restrict delz = delrho[2].ptr();
                          for (unsigned int j=0; j<np; j++) {
        //                   double grad = delx[j] + dely[j] + delz[j] ;
        
        //vama6                   double grad = sqrt(delx[j]*delx[j] +
        //vama6                                      dely[j]*dely[j] +
        //vama6                                      delz[j]*delz[j] );
        
        //yyp                     res[j] = (vsax[j])*(grad);
                                  res[j] = vsaaa[j]*funcs[i].second;
                            //res[j] += (vsax[j])*(delx[j]+dely[j]+delz[j])*2.0;
                                  if (isnan_x(res[j])) throw "ouch";
                          }
                    }
                }
                    break;
                default:
                    throw "UGH!";
             }
        }
    }
    return result;
}


