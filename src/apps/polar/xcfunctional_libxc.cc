#include <madness/madness_config.h>
#include <chem/xcfunctional.h>
#include <madness/tensor/tensor.h>
#include <iostream>
#include <string>
#include <sstream>
#include <utility>
#include <madness/world/madness_excecption.h>
#include <madness/world/MADworld.h>
#include <xc.h>
#include <xc_funcs.h>

// This function is here because compiling with -ffast-math breaks the C
// function isnan. Also, there is a bug in some compilers where isnan is
// undefined when <cmath> is included.
namespace {
    inline int isnan_x(double x) {
        volatile double y = x;
        return x != y;
    }
}

namespace madness {
struct xc_name_map {
    const std::string name;
    const int id;
};

static xc_name_map map[] = {
    {"LDA_X",                   XC_LDA_X},
    {"LDA_C_WIGNER",            XC_LDA_C_WIGNER},
    {"LDA_C_RPA",               XC_LDA_C_RPA},
    {"LDA_C_HL",                XC_LDA_C_HL},
    {"LDA_C_GL",                XC_LDA_C_GL},
    {"LDA_C_XALPHA",            XC_LDA_C_XALPHA},
    {"LDA_C_VWN",               XC_LDA_C_VWN},
    {"LDA_C_VWN_RPA",           XC_LDA_C_VWN_RPA},
    {"LDA_C_PZ",                XC_LDA_C_PZ},
    {"LDA_C_PZ_MOD",            XC_LDA_C_PZ_MOD},
    {"LDA_C_OB_PZ",             XC_LDA_C_OB_PZ},
    {"LDA_C_PW",                XC_LDA_C_PW},
    {"LDA_C_PW_MOD",            XC_LDA_C_PW_MOD},
    {"LDA_C_OB_PW",             XC_LDA_C_OB_PW},
    {"LDA_C_2D_AMGB",           XC_LDA_C_2D_AMGB},
    {"LDA_C_2D_PRM",            XC_LDA_C_2D_PRM},
    {"LDA_C_vBH",               XC_LDA_C_vBH},
    {"LDA_C_1D_CSC",            XC_LDA_C_1D_CSC},
    {"LDA_X_2D",                XC_LDA_X_2D},
    {"LDA_XC_TETER93",          XC_LDA_XC_TETER93},
    {"LDA_X_1D",                XC_LDA_X_1D},
    {"LDA_C_ML1",               XC_LDA_C_ML1},
    {"LDA_C_ML2",               XC_LDA_C_ML2},
    {"LDA_C_GOMBAS",            XC_LDA_C_GOMBAS},
    {"LDA_K_TF",                XC_LDA_K_TF},
    {"LDA_K_LP",                XC_LDA_K_LP},
    {"GGA_X_PBE",               XC_GGA_X_PBE},
    {"GGA_X_PBE_R",             XC_GGA_X_PBE_R},
    {"GGA_X_B86",               XC_GGA_X_B86},
    {"GGA_X_HERMAN",            XC_GGA_X_HERMAN},
    {"GGA_X_B86_MGC",           XC_GGA_X_B86_MGC},
    {"GGA_X_B88",               XC_GGA_X_B88},
    {"GGA_X_G96",               XC_GGA_X_G96},
    {"GGA_X_PW86",              XC_GGA_X_PW86},
    {"GGA_X_PW91",              XC_GGA_X_PW91},
    {"GGA_X_OPTX",              XC_GGA_X_OPTX},
    {"GGA_X_DK87_R1",           XC_GGA_X_DK87_R1},
    {"GGA_X_DK87_R2",           XC_GGA_X_DK87_R2},
    {"GGA_X_LG93",              XC_GGA_X_LG93},
    {"GGA_X_FT97_A",            XC_GGA_X_FT97_A},
    {"GGA_X_FT97_B",            XC_GGA_X_FT97_B},
    {"GGA_X_PBE_SOL",           XC_GGA_X_PBE_SOL},
    {"GGA_X_RPBE",              XC_GGA_X_RPBE},
    {"GGA_X_WC",                XC_GGA_X_WC},
    {"GGA_X_AM05",              XC_GGA_X_AM05},
    {"GGA_X_PBEA",              XC_GGA_X_PBEA},
    {"GGA_X_MPBE",              XC_GGA_X_MPBE},
    {"GGA_X_XPBE",              XC_GGA_X_XPBE},
    {"GGA_X_2D_B86_MGC",        XC_GGA_X_2D_B86_MGC},
    {"GGA_X_BAYESIAN",          XC_GGA_X_BAYESIAN},
    {"GGA_X_PBE_JSJR",          XC_GGA_X_PBE_JSJR},
    {"GGA_X_2D_B88",            XC_GGA_X_2D_B88},
    {"GGA_X_2D_B86",            XC_GGA_X_2D_B86},
    {"GGA_X_2D_PBE",            XC_GGA_X_2D_PBE},
    {"GGA_C_PBE",               XC_GGA_C_PBE},
    {"GGA_C_LYP",               XC_GGA_C_LYP},
    {"GGA_C_P86",               XC_GGA_C_P86},
    {"GGA_C_PBE_SOL",           XC_GGA_C_PBE_SOL},
    {"GGA_C_PW91",              XC_GGA_C_PW91},
    {"GGA_C_AM05",              XC_GGA_C_AM05},
    {"GGA_C_XPBE",              XC_GGA_C_XPBE},
    {"GGA_C_LM",                XC_GGA_C_LM},
    {"GGA_C_PBE_JRGX",          XC_GGA_C_PBE_JRGX},
    {"GGA_X_OPTB88_VDW",        XC_GGA_X_OPTB88_VDW},
    {"GGA_X_PBEK1_VDW",         XC_GGA_X_PBEK1_VDW},
    {"GGA_X_OPTPBE_VDW",        XC_GGA_X_OPTPBE_VDW},
    {"GGA_X_RGE2",              XC_GGA_X_RGE2},
    {"GGA_C_RGE2",              XC_GGA_C_RGE2},
    {"GGA_X_RPW86",             XC_GGA_X_RPW86},
    {"GGA_X_KT1",               XC_GGA_X_KT1},
    {"GGA_XC_KT2",              XC_GGA_XC_KT2},
    {"GGA_C_WL",                XC_GGA_C_WL},
    {"GGA_C_WI",                XC_GGA_C_WI},
    {"GGA_X_MB88",              XC_GGA_X_MB88},
    {"GGA_X_LB",                XC_GGA_X_LB},
    {"GGA_XC_HCTH_93",          XC_GGA_XC_HCTH_93},
    {"GGA_XC_HCTH_120",         XC_GGA_XC_HCTH_120},
    {"GGA_XC_HCTH_147",         XC_GGA_XC_HCTH_147},
    {"GGA_XC_HCTH_407",         XC_GGA_XC_HCTH_407},
    {"GGA_XC_EDF1",             XC_GGA_XC_EDF1},
    {"GGA_XC_XLYP",             XC_GGA_XC_XLYP},
    {"GGA_XC_B97",              XC_GGA_XC_B97},
    {"GGA_XC_B97_1",            XC_GGA_XC_B97_1},
    {"GGA_XC_B97_2",            XC_GGA_XC_B97_2},
    {"GGA_XC_B97_D",            XC_GGA_XC_B97_D},
    {"GGA_XC_B97_K",            XC_GGA_XC_B97_K},
    {"GGA_XC_B97_3",            XC_GGA_XC_B97_3},
    {"GGA_XC_PBE1W",            XC_GGA_XC_PBE1W},
    {"GGA_XC_MPWLYP1W",         XC_GGA_XC_MPWLYP1W},
    {"GGA_XC_PBELYP1W",         XC_GGA_XC_PBELYP1W},
    {"GGA_XC_SB98_1a",          XC_GGA_XC_SB98_1a},
    {"GGA_XC_SB98_1b",          XC_GGA_XC_SB98_1b},
    {"GGA_XC_SB98_1c",          XC_GGA_XC_SB98_1c},
    {"GGA_XC_SB98_2a",          XC_GGA_XC_SB98_2a},
    {"GGA_XC_SB98_2b",          XC_GGA_XC_SB98_2b},
    {"GGA_XC_SB98_2c",          XC_GGA_XC_SB98_2c},
    {"GGA_X_LBM",               XC_GGA_X_LBM},
    {"GGA_X_OL2",               XC_GGA_X_OL2},
    {"GGA_X_APBE",              XC_GGA_X_APBE},
    {"GGA_K_APBE",              XC_GGA_K_APBE},
    {"GGA_C_APBE",              XC_GGA_C_APBE},
    {"GGA_K_TW1",               XC_GGA_K_TW1},
    {"GGA_K_TW2",               XC_GGA_K_TW2},
    {"GGA_K_TW3",               XC_GGA_K_TW3},
    {"GGA_K_TW4",               XC_GGA_K_TW4},
    {"GGA_K_VW",                XC_GGA_K_VW},
    {"GGA_K_GE2",               XC_GGA_K_GE2},
    {"GGA_K_GOLDEN",            XC_GGA_K_GOLDEN},
    {"GGA_K_YT65",              XC_GGA_K_YT65},
    {"GGA_K_BALTIN",            XC_GGA_K_BALTIN},
    {"GGA_K_LIEB",              XC_GGA_K_LIEB},
    {"GGA_K_ABSR1",             XC_GGA_K_ABSR1},
    {"GGA_K_ABSR2",             XC_GGA_K_ABSR2},
    {"GGA_K_GR",                XC_GGA_K_GR},
    {"GGA_K_LUDENA",            XC_GGA_K_LUDENA},
    {"GGA_K_GP85",              XC_GGA_K_GP85},
    {"GGA_K_PEARSON",           XC_GGA_K_PEARSON},
    {"GGA_K_OL1",               XC_GGA_K_OL1},
    {"GGA_K_OL2",               XC_GGA_K_OL2},
    {"GGA_K_FR_B88",            XC_GGA_K_FR_B88},
    {"GGA_K_FR_PW86",           XC_GGA_K_FR_PW86},
    {"GGA_K_DK",                XC_GGA_K_DK},
    {"GGA_K_PERDEW",            XC_GGA_K_PERDEW},
    {"GGA_K_VSK",               XC_GGA_K_VSK},
    {"GGA_K_VJKS",              XC_GGA_K_VJKS},
    {"GGA_K_ERNZERHOF",         XC_GGA_K_ERNZERHOF},
    {"GGA_K_LC94",              XC_GGA_K_LC94},
    {"GGA_K_LLP",               XC_GGA_K_LLP},
    {"GGA_K_THAKKAR",           XC_GGA_K_THAKKAR},
    {"HYB_GGA_XC_B3PW91",       XC_HYB_GGA_XC_B3PW91},
    {"HYB_GGA_XC_B3LYP",        XC_HYB_GGA_XC_B3LYP},
    {"HYB_GGA_XC_B3P86",        XC_HYB_GGA_XC_B3P86},
    {"HYB_GGA_XC_O3LYP",        XC_HYB_GGA_XC_O3LYP},
    {"HYB_GGA_XC_mPW1K",        XC_HYB_GGA_XC_mPW1K},
    {"HYB_GGA_XC_PBEH",         XC_HYB_GGA_XC_PBEH},
    {"HYB_GGA_XC_B97",          XC_HYB_GGA_XC_B97},
    {"HYB_GGA_XC_B97_1",        XC_HYB_GGA_XC_B97_1},
    {"HYB_GGA_XC_B97_2",        XC_HYB_GGA_XC_B97_2},
    {"HYB_GGA_XC_X3LYP",        XC_HYB_GGA_XC_X3LYP},
    {"HYB_GGA_XC_B1WC",         XC_HYB_GGA_XC_B1WC},
    {"HYB_GGA_XC_B97_K",        XC_HYB_GGA_XC_B97_K},
    {"HYB_GGA_XC_B97_3",        XC_HYB_GGA_XC_B97_3},
    //{"HYB_GGA_XC_mPW3PW",       XC_HYB_GGA_XC_mPW3PW},
    {"HYB_GGA_XC_B1LYP",        XC_HYB_GGA_XC_B1LYP},
    {"HYB_GGA_XC_B1PW91",       XC_HYB_GGA_XC_B1PW91},
    {"HYB_GGA_XC_mPW1PW",       XC_HYB_GGA_XC_mPW1PW},
    //{"HYB_GGA_XC_mPW3LYP",      XC_HYB_GGA_XC_mPW3LYP},
    {"HYB_GGA_XC_SB98_1a",      XC_HYB_GGA_XC_SB98_1a},
    {"HYB_GGA_XC_SB98_1b",      XC_HYB_GGA_XC_SB98_1b},
    {"HYB_GGA_XC_SB98_1c",      XC_HYB_GGA_XC_SB98_1c},
    {"HYB_GGA_XC_SB98_2a",      XC_HYB_GGA_XC_SB98_2a},
    {"HYB_GGA_XC_SB98_2b",      XC_HYB_GGA_XC_SB98_2b},
    {"HYB_GGA_XC_SB98_2c",      XC_HYB_GGA_XC_SB98_2c},
    {"MGGA_X_LTA",              XC_MGGA_X_LTA},
    {"MGGA_X_TPSS",             XC_MGGA_X_TPSS},
    //{"MGGA_X_M06L",             XC_MGGA_X_M06L},
    {"MGGA_X_GVT4",             XC_MGGA_X_GVT4},
    {"MGGA_X_TAU_HCTH",         XC_MGGA_X_TAU_HCTH},
    {"MGGA_X_BR89",             XC_MGGA_X_BR89},
    {"MGGA_X_BJ06",             XC_MGGA_X_BJ06},
    {"MGGA_X_TB09",             XC_MGGA_X_TB09},
    {"MGGA_X_RPP09",            XC_MGGA_X_RPP09},
    {"MGGA_X_2D_PRHG07",        XC_MGGA_X_2D_PRHG07},
    {"MGGA_X_2D_PRHG07_PRP10",  XC_MGGA_X_2D_PRHG07_PRP10},
    {"MGGA_C_TPSS",             XC_MGGA_C_TPSS},
    {"MGGA_C_VSXC",             XC_MGGA_C_VSXC},

    //{"GGA_X_mPW91",             XC_GGA_X_mPW91},
    //{"GGA_X_mPW91", XC_GGA_X_mPW91},
    //{"HYB_GGA_XC_mPW3PW", XC_HYB_GGA_XC_mPW3PW},
    //{"HYB_GGA_XC_mPW3LYP", XC_HYB_GGA_XC_mPW3LYP},
    //{"MGGA_X_M06L", XC_MGGA_X_M06L},

    {"", -1} // MUST BE LAST
};

static int lookup_name(const std::string& name) {
    const xc_name_map* map = madness::map;
    while (map->id > 0) {
        if (name == map->name) return map->id;
        map++;
    }
    return -1;
}


static xc_func_type* make_func(int id, bool polarized) {
    xc_func_type* func = new xc_func_type;
    int POLARIZED = polarized ? XC_POLARIZED : XC_UNPOLARIZED;
    MADNESS_ASSERT(xc_func_init(func, id, POLARIZED) == 0);
    return func;
}

static xc_func_type* lookup_func(const std::string& name, bool polarized) {
    int id = lookup_name(name);
    MADNESS_ASSERT(id > 0);
    return make_func(id, polarized);
}

//XCfunctional::XCfunctional() {}
//XCfunctional::XCfunctional() : hf_coeff(0.0) {std::printf("Construct XC Functional from LIBXC Library");}
XCfunctional::XCfunctional() : hf_coeff(0.0) {
    rhotol=1e-7; rhomin=0.0; sigtol=1e-10; sigmin=1e-10; // default values
}

void XCfunctional::initialize(const std::string& input_line, bool polarized, World& world)
{
    rhotol=1e-7; rhomin=0.0; sigtol=1e-10; sigmin=1e-10; // default values

    double factor;
    spin_polarized = polarized;


    std::stringstream line(input_line);
    std::string name;

    nderiv = 0;
    hf_coeff = 0.0;
    funcs.clear();

    if (world.rank() == 0) std::printf("Construct XC Functional from LIBXC Library");
    while (line >> name) {
        std::transform(name.begin(), name.end(), name.begin(), ::toupper);
        if (world.rank() == 0) std::cout <<"!NAME! "<< name << "pol " << polarized << std::endl;
        if (name == "LDA") {
            //if (! (line >> factor)) factor = 1.0;
            funcs.push_back(std::make_pair(lookup_func("LDA_X",polarized),1.0));
            funcs.push_back(std::make_pair(lookup_func("LDA_C_VWN",polarized),1.0));
        }
        else if (name == "RHOMIN") {
        	std::cout << "hello" << std::endl;
            line >> rhomin;
        }
        else if (name == "RHOTOL") {
            line >> rhotol;
        }
        else if (name == "SIGMIN") {
            line >> sigmin;
        }
        else if (name == "SIGTOL") {
            line >> sigtol;
        }
        else if (name == "HF" || name == "HF_X") {
            if (! (line >> factor)) factor = 1.0;
            hf_coeff = factor;
        }
        else {
            if (! (line >> factor)) factor = 1.0;
            funcs.push_back(std::make_pair(lookup_func(name,polarized), factor));
        }
    }

    for (unsigned int i=0; i<funcs.size(); i++) {
        if (funcs[i].first->info->family == XC_FAMILY_GGA) nderiv = std::max(nderiv,1);
        if (funcs[i].first->info->family == XC_FAMILY_HYB_GGA) nderiv = std::max(nderiv,1);
        if (funcs[i].first->info->family == XC_FAMILY_MGGA)nderiv = std::max(nderiv,2);
 //       if (funcs[i].first->info->family == XC_FAMILY_LDA) nderiv = std::max(nderiv,0);
        if (world.rank() == 0) std::cout << "factor " << i << "  " << funcs[i].second << std::endl;
    }
    if (world.rank() == 0) std::cout << "rhotol " << rhotol << " rhomin " << rhomin << " factor " <<factor << "hfcorf" <<hf_coeff <<  " input line was " << input_line << std::endl;
}

XCfunctional::~XCfunctional() {
    for (unsigned int i=0; i<funcs.size(); i++) {
        xc_func_end(funcs[i].first);
        delete funcs[i].first;
    }
    funcs.clear();
}

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
    return false; // not thought about this yet
}

bool XCfunctional::has_kxc() const
{
    return false;
}


/// Allocates rho (and if GGA also sigma) and copies data from t[] into rho and sigma.
//vama5void XCfunctional::make_libxc_args(const std::vector< madness::Tensor<double> >& t,
//vama5                                   madness::Tensor<double>& rho, madness::Tensor<double>& sigma, madness::Tensor<double>& delrho) const
void XCfunctional::make_libxc_args(const std::vector< madness::Tensor<double> >& t,
                                   madness::Tensor<double>& rho, madness::Tensor<double>& sigma) const
{
    const int np = t[0].size();
    if (spin_polarized) {
        if (is_lda()) {
            MADNESS_ASSERT(t.size() == 2);
            const double * restrict rhoa = t[0].ptr();
            const double * restrict rhob = t[1].ptr();
            rho  = madness::Tensor<double>(np*2L);
            double * restrict dens = rho.ptr();
            for (long i=0; i<np; i++) {
                dens[2*i  ] = munge(rhoa[i]);
                dens[2*i+1] = munge(rhob[i]);
            }
        }
        else if (is_gga()) {
            MADNESS_ASSERT(t.size() == 5);
            const double * restrict rhoa  = t[0].ptr();
            const double * restrict rhob  = t[1].ptr();

            const double * restrict sigaa = t[2].ptr();
            const double * restrict sigab = t[3].ptr();
            const double * restrict sigbb = t[4].ptr();


            rho   = madness::Tensor<double>(np*2L);
            sigma = madness::Tensor<double>(np*3L);

            double * restrict dens = rho.ptr();
            double * restrict sig  = sigma.ptr();
            for (long i=0; i<np; i++) {
                double ra=rhoa[i], rb=rhob[i], saa=sigaa[i], sab=sigab[i], sbb=sigbb[i];

                munge5(ra, rb, saa, sab, sbb);
                dens[2*i  ] = ra;
                dens[2*i+1] = rb;

                sig[3*i  ] = saa;
                sig[3*i+1] = sab;
                sig[3*i+2] = sbb;

            }
        }
        else {
            throw "not yet";
        }
    }
    else {
        if (is_lda()) {
            MADNESS_ASSERT(t.size() == 1);
            rho  = madness::Tensor<double>(np);
            const double * restrict rhoa = t[0].ptr();
            double * restrict dens = rho.ptr();
            for (long i=0; i<np; i++) {
                dens[i] = munge(2.0*rhoa[i]);
            }
        }
        else if (is_gga()) {
            MADNESS_ASSERT(t.size() == 2);
            const double * restrict rhoa = t[0].ptr();
            const double * restrict sigaa = t[1].ptr();
            rho  = madness::Tensor<double>(np);
            sigma  = madness::Tensor<double>(np);
            double * restrict dens = rho.ptr();
            double * restrict sig = sigma.ptr();
            for (long i=0; i<np; i++) {
                double ra=2.0*rhoa[i], saa=4.0*sigaa[i];
                //double ra=rhoa[i], saa=sigaa[i];
                munge2(ra, saa);
                dens[i] = ra;
                sig[i] = saa;
            }
        }
        else {
            throw "not yet";
        }
    }
}


madness::Tensor<double> XCfunctional::exc(const std::vector< madness::Tensor<double> >& t, const int ispin) const
{
    madness::Tensor<double> rho, sigma;
    make_libxc_args(t, rho, sigma);

    const int np = t[0].size();
    const double * restrict dens = rho.ptr();
    const double * restrict sig = sigma.ptr();

    madness::Tensor<double> result(3L, t[0].dims());
    double * restrict res = result.ptr();
    for (long j=0; j<np; j++) res[j] = 0.0;

    for (unsigned int i=0; i<funcs.size(); i++) {
        madness::Tensor<double> zk(3L, t[0].dims(), false);
        double * restrict work = zk.ptr();

        switch(funcs[i].first->info->family) {
        case XC_FAMILY_LDA:
            xc_lda_exc(funcs[i].first, np, dens, work);
            break;
        case XC_FAMILY_GGA:
            xc_gga_exc(funcs[i].first, np, dens, sig, work);
            break;
        case XC_FAMILY_HYB_GGA:
            xc_gga_exc(funcs[i].first, np, dens, sig, work);
            break;
        default:
            throw "HOW DID WE GET HERE?";
        }
        if (spin_polarized) {
            for (long j=0; j<np; j++) {
                res[j] +=  work[j]*(dens[2*j+1] + dens[2*j])*funcs[i].second;
                //res[j] +=  work[j]*dens[2*j+ispin]*funcs[i].second;
                //res[j] +=  work[j];
            }
        }
        else {
            for (long j=0; j<np; j++) {
                res[j] += work[j]*dens[j]*funcs[i].second;
            }
        }
    }

    return result;
}

madness::Tensor<double> XCfunctional::vxc(const std::vector< madness::Tensor<double> >& t, const int ispin, const int what) const
{
    madness::Tensor<double> rho, sigma;
    make_libxc_args(t, rho, sigma);

    const int np = t[0].size();

    int nvsig, nvrho;
    if (spin_polarized) {
    	nvrho = 2;
    	nvsig = 3;
    }
    else {
    	nvrho = 1;
    	nvsig = 1;
    }

    madness::Tensor<double> result(3L, t[0].dims());
    double * restrict res = result.ptr();
    const double * restrict dens = rho.ptr();
	for (long j=0; j<np; j++) 
    			res[j] = 0.0;

    for (unsigned int i=0; i<funcs.size(); i++) {
       	switch(funcs[i].first->info->family) {
       	     case XC_FAMILY_LDA:
                {
                madness::Tensor<double> vrho(nvrho*np);
    		double * restrict vr = vrho.ptr();
    		xc_lda_vxc(funcs[i].first, np, dens, vr);
             	if (what < 2) {
                     for (long j=0; j<np; j++) {
    			     	res[j] += vr[nvrho*j+ispin]*funcs[i].second;
    			     	if (isnan_x(res[j])) throw "ouch";
    			     }
    			}
    		}

		break;

             case XC_FAMILY_GGA:
    		{
    	        madness::Tensor<double> vrho(nvrho*np), vsig(nvsig*np);
                double * restrict vr = vrho.ptr();
                double * restrict vs = vsig.ptr();
    	        const double * restrict sig = sigma.ptr();
	    	xc_gga_vxc(funcs[i].first, np, dens, sig, vr, vs);
	    	if (spin_polarized) {
        		if (what == 0) {
        			for (long j=0; j<np; j++) {
    					res[j] += vr[nvrho*j+ispin] * funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else if (what == 1) {
    			// Vaa
    				for (long j=0; j<np; j++) {
    					res[j] += vs[nvsig*j + 2*ispin] * funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else if (what == 2) {
    			// Vab
    				for (long j=0; j<np; j++) {
    					res[j] += vs[nvsig*j + 1] * funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else {
    				throw "ouch";
    			}
    		}
    		else {
    			if (what == 0) {
    				for (long j=0; j<np; j++) {
    					res[j] += vr[j]*funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else if (what == 1) {
    			// Vaa
    				for (long j=0; j<np; j++) {
    					res[j] += vs[j]*funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else {
    				throw "ouch";
    			}
    		}
   		}
    		break;
     	case XC_FAMILY_HYB_GGA:
    	{
    				throw "ouch";
    		madness::Tensor<double> vrho(nvrho*np), vsig(nvsig*np);
    		double * restrict vr = vrho.ptr();
    		double * restrict vs = vsig.ptr();
    		const double * restrict sig = sigma.ptr();
    		xc_gga_vxc(funcs[i].first, np, dens, sig, vr, vs);
	    	if (spin_polarized) {
        		if (what == 0) {
        			for (long j=0; j<np; j++) {
    					res[j] += vr[nvrho*j+ispin]*funcs[i].second*.0;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else if (what == 1) {
    			// Vaa
    				for (long j=0; j<np; j++) {
    					res[j] += vs[nvsig*j + 2*ispin]*funcs[i].second*0.0;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else if (what == 2) {
    			// Vab
    				for (long j=0; j<np; j++) {
    					res[j] += vs[nvsig*j + 1]*funcs[i].second*.0;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else {
    				throw "ouch";
    			}
    		}
    		else {
    			if (what == 0) {
    				for (long j=0; j<np; j++) {
    					res[j] += vr[j]*funcs[i].second*.0;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else if (what == 1) {
    				throw "ouch";
    			// Vaa
    				for (long j=0; j<np; j++) {
    					res[j] += vs[j]*funcs[i].second*.0;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else {
    				throw "ouch";
    			}
    		}
    	}
    	break;
    	default:
    	throw "UGH!";
    	}
    }
    return result;
}

madness::Tensor<double> XCfunctional::fxc(const std::vector< madness::Tensor<double> >& t, const int ispin, const int what) const
{
    madness::Tensor<double> rho, sigma;
    make_libxc_args(t, rho, sigma);

    const int np = t[0].size();

    int nvsig, nvrho;
    if (spin_polarized) {
    	nvrho = 2;
    	nvsig = 3;
    }
    else {
    	nvrho = 1;
    	nvsig = 1;
    }

    madness::Tensor<double> result(3L, t[0].dims());
    double * restrict res = result.ptr();
    const double * restrict dens = rho.ptr();
    for (long j=0; j<np; j++) 
    			res[j] = 0.0;

    if(what ==4){
       for (long j=0; j<np; j++) {
  	     	res[j] = dens[j]*.5;
  	     	if (isnan_x(res[j])) throw "ouch";
	     }
    }
    else {
    for (unsigned int i=0; i<funcs.size(); i++) {
       	switch(funcs[i].first->info->family) {
       	     case XC_FAMILY_LDA:
                {
                madness::Tensor<double> vrho(nvrho*np);
    		double * restrict vr = vrho.ptr();
    		xc_lda_fxc(funcs[i].first, np, dens, vr);
             	if (what < 2) {
                     for (long j=0; j<np; j++) {
    			     	res[j] += vr[nvrho*j+ispin]*funcs[i].second;
    			     	//res[j] += vr[nvrho*j+ispin]*funcs[i].second*dens[j];
    			     	//res[j] = vr[j];
    			     	//res[j] += vr[nvrho*j+ispin]*funcs[i].second;
    			     	if (isnan_x(res[j])) throw "ouch";
    			     }
    			}
    		}

		break;

             case XC_FAMILY_GGA:
    		{
    				throw "ouch";
   		}
    		break;
         	case XC_FAMILY_HYB_GGA:
        	{
       				throw "ouch";
        	}
    	break;
    	default:
    	throw "UGH!";
    	}
    }
    }
    return result;
}
}
