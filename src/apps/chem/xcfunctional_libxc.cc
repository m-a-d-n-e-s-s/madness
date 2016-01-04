#include <madness/madness_config.h>
#include <chem/xcfunctional.h>
#include <madness/tensor/tensor.h>
#include <iostream>
#include <string>
#include <sstream>
#include <utility>
#include <madness/world/madness_exception.h>
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

static std::string lookup_id(const int id) {
    const xc_name_map* map = madness::map;
    while (map->id > 0) {
        if (id == map->id) return map->name;
        map++;
    }
    return "Functional not found";
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
    rhotol=1e-7; rhomin=0.0; sigtol=1e-10; sigmin=1e-10;
    ggatol=1.e-4;
    munge_ratio=10.0;
    nderiv=0;
    spin_polarized=false;
}

void XCfunctional::initialize(const std::string& input_line, bool polarized,
        World& world, const bool verbose) {
    rhotol=1e-7; rhomin=0.0; sigtol=1e-10; sigmin=1e-10; // default values
    ggatol=1.e-4;

    bool printit=verbose and (world.rank()==0);
    double factor;      // weight factor for the various functionals
    spin_polarized = polarized;


    std::stringstream line(input_line);
    std::string name;

    nderiv = 0;
    hf_coeff = 0.0;
    funcs.clear();

    if (printit) print("\nConstruct XC Functional from LIBXC Library");
    while (line >> name) {
        std::transform(name.begin(), name.end(), name.begin(), ::toupper);
        if (name == "LDA") {
            // Slater exchange and VWN-5 correlation
            funcs.push_back(std::make_pair(lookup_func("LDA_X",polarized),1.0));
            funcs.push_back(std::make_pair(lookup_func("LDA_C_VWN",polarized),1.0));
        } else if ((name == "BP86") or (name=="BP")) {
            // Becke exchange, VWN-3 correlation, Perdew correction
            funcs.push_back(std::make_pair(lookup_func("GGA_X_B88",polarized),1.0));
            funcs.push_back(std::make_pair(lookup_func("GGA_C_P86",polarized),1.0));
        } else if (name == "PBE") {
            funcs.push_back(std::make_pair(lookup_func("GGA_X_PBE",polarized),1.0));
            funcs.push_back(std::make_pair(lookup_func("GGA_C_PBE",polarized),1.0));
        } else if (name == "PBE0") {
            funcs.push_back(std::make_pair(lookup_func("GGA_X_PBE",polarized),0.75));
            funcs.push_back(std::make_pair(lookup_func("GGA_C_PBE",polarized),1.0));
            hf_coeff=0.25;
        } else if (name == "B3LYP") {
            // VWN-3 correlation
            funcs.push_back(std::make_pair(lookup_func("HYB_GGA_XC_B3LYP",polarized),1.0));
            hf_coeff=0.2;
        } else if (name == "RHOMIN") {
            line >> rhomin;
        } else if (name == "RHOTOL") {
            line >> rhotol;
        } else if (name == "SIGMIN") {
            line >> sigmin;
        } else if (name == "SIGTOL") {
            line >> sigtol;
        } else if (name == "GGATOL") {
            line >> ggatol;
        } else if (name == "MUNGERATIO") {
            line >> munge_ratio;
        } else if (name == "HF" || name == "HF_X") {
            if (! (line >> factor)) factor = 1.0;
            hf_coeff = factor;
        } else {
            if (! (line >> factor)) factor = 1.0;
            funcs.push_back(std::make_pair(lookup_func(name,polarized), factor));
        }
    }

    for (unsigned int i=0; i<funcs.size(); i++) {
        if (funcs[i].first->info->family == XC_FAMILY_GGA) nderiv = std::max(nderiv,1);
        if (funcs[i].first->info->family == XC_FAMILY_HYB_GGA) nderiv = std::max(nderiv,1);
        if (funcs[i].first->info->family == XC_FAMILY_MGGA)nderiv = std::max(nderiv,2);
 //       if (funcs[i].first->info->family == XC_FAMILY_LDA) nderiv = std::max(nderiv,0);
    }
    if (printit) {
        print("\ninput line was:",input_line);
        for (std::size_t i=0; i<funcs.size(); ++i) {
            int id=funcs[i].first->info->number;
//            print(lookup_id(id),"with weight:",funcs[i].second);
            printf(" %4.3f %s \n",funcs[i].second,lookup_id(id).c_str());
        }
        if (hf_coeff>0.0) printf(" %4.3f %s \n",hf_coeff,"HF exchange");
        print("\nscreening parameters");
        print(" rhotol, rhomin",rhotol,rhomin);
        print(" sigtol, sigmin",sigtol,sigmin);
        print("         ggatol",ggatol);
        print("    munge_ratio",munge_ratio);
        if (printit) print("polarized ",polarized,"\n");

    }
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
//    return (is_lda() || is_gga() || is_meta());
    return (funcs.size()>0);
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


/// Items in the vector argument \c t are interpreted as follows
///  - Spin un-polarized
///    - \c t[0] = \f$ \rho_{\alpha} \f$
///    - \c t[1] = \f$ \sigma_{\alpha\alpha} = \nabla \rho_{\alpha}.\nabla \rho_{\alpha} \f$ (GGA only)
///    - \c t[2] = \f$ \tilde \rho \f$ (the (perturbed) density for xc_kernel_apply)
///  - Spin polarized
///    - \c t[0] = \f$ \rho_{\alpha} \f$
///    - \c t[1] = \f$ \rho_{\beta} \f$
///    - \c t[2] = \f$ \sigma_{\alpha\alpha} = \nabla \rho_{\alpha}.\nabla \rho_{\alpha} \f$ (GGA only)
///    - \c t[3] = \f$ \sigma_{\alpha\beta}  = \nabla \rho_{\alpha}.\nabla \rho_{\beta} \f$ (GGA only)
///    - \c t[4] = \f$ \sigma_{\beta\beta}   = \nabla \rho_{\beta}.\nabla \rho_{\beta} \f$ (GGA only)
///    - \c t[5] = \f$ \tilde \rho \f$ (the (perturbed) density for xc_kernel_apply)
///
/// output
///  - if spin-unpolarized:
///    - rho = 2 t[0] ( = 2 rho_alpha)
///    - sigma = 4 sigma_aa = \nabla \rho . \nabla \rho
///  - if spin-polarized
///    - rho[2*npt] = rho[a,b]
///    - sigma[3*npt] = sigma[aa,ab,bb]
///    leading dimension is [a,b] and [aa,ab,bb], respectively (why??!)
void XCfunctional::make_libxc_args_old(const std::vector< madness::Tensor<double> >& t,
           madness::Tensor<double>& rho, madness::Tensor<double>& sigma,
           const munging_type& munging) const {
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

            // might happen if there are no beta electrons
            madness::Tensor<double> dummy;
            if ((rhob==NULL) or (sigab==NULL) or (sigbb==NULL)) {
                dummy=madness::Tensor<double>(np);
            }
            if (rhob==NULL) rhob=dummy.ptr();
            if (sigab==NULL) sigab=dummy.ptr();
            if (sigbb==NULL) sigbb=dummy.ptr();

            rho   = madness::Tensor<double>(np*2L);
            sigma = madness::Tensor<double>(np*3L);

            double * restrict dens = rho.ptr();
            double * restrict sig  = sigma.ptr();
            for (long i=0; i<np; i++) {
                double ra=rhoa[i], rb=rhob[i], saa=sigaa[i], sab=sigab[i], sbb=sigbb[i];

                munge5(ra, rb, saa, sab, sbb, munging);
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
                munge2(ra, saa, munging);
                dens[i] = ra;
                sig[i] = saa;
            }
        }
        else {
            throw "not yet";
        }
    }
}


void XCfunctional::make_libxc_args(const std::vector< madness::Tensor<double> >& xc_args,
           madness::Tensor<double>& rho, madness::Tensor<double>& sigma,
           const munging_type& munging) const {
    const int np = xc_args[0].size();


    if (not spin_polarized) {
        if (is_lda()) {
            rho  = madness::Tensor<double>(np);
            const double * restrict rhoa = xc_args[enum_rhoa].ptr();
            double * restrict dens = rho.ptr();
            for (long i=0; i<np; i++) {
                dens[i] = munge(2.0*rhoa[i]);
            }
        }
        else if (is_gga()) {
            const double * restrict rhoa = xc_args[enum_rhoa].ptr();
            const double * restrict sigaa = xc_args[enum_saa].ptr();
            rho  = madness::Tensor<double>(np);
            sigma  = madness::Tensor<double>(np);
            double * restrict dens = rho.ptr();
            double * restrict sig = sigma.ptr();
            for (long i=0; i<np; i++) {
                double ra=2.0*rhoa[i], saa=4.0*sigaa[i];
                munge2(ra, saa, munging);
                dens[i] = ra;
                sig[i] = saa;
            }
        }
        else {
            MADNESS_EXCEPTION("only LDA and GGA available in xcfunctional",1);
        }

    } else if (spin_polarized) {
        if (is_lda()) {
            const double * restrict rhoa = xc_args[enum_rhoa].ptr();
            const double * restrict rhob = xc_args[enum_rhob].ptr();
            rho  = madness::Tensor<double>(np*2L);
            double * restrict dens = rho.ptr();
            for (long i=0; i<np; i++) {
                dens[2*i  ] = munge(rhoa[i]);
                dens[2*i+1] = munge(rhob[i]);
            }
        }
        else if (is_gga()) {
            const double * restrict rhoa  = xc_args[enum_rhoa].ptr();
            const double * restrict rhob  = xc_args[enum_rhob].ptr();

            const double * restrict sigaa = xc_args[enum_saa].ptr();
            const double * restrict sigab = xc_args[enum_sab].ptr();
            const double * restrict sigbb = xc_args[enum_sbb].ptr();
            const double * restrict sigtot = xc_args[enum_sigtot].ptr();

            // might happen if there are no beta electrons
            madness::Tensor<double> dummy;
            if ((rhob==NULL) or (sigab==NULL) or (sigbb==NULL)) {
                dummy=madness::Tensor<double>(np);
            }
            if (rhob==NULL) rhob=dummy.ptr();
            if (sigab==NULL) sigab=dummy.ptr();
            if (sigbb==NULL) sigbb=dummy.ptr();

            rho   = madness::Tensor<double>(np*2L);
            sigma = madness::Tensor<double>(np*3L);

            double * restrict dens = rho.ptr();
            double * restrict sig  = sigma.ptr();
            for (long i=0; i<np; i++) {
                double ra=rhoa[i], rb=rhob[i], r=ra+rb;
                double saa=sigaa[i], sab=sigab[i], sbb=sigbb[i], stot=sigtot[i];

                munge7(ra, rb, r, saa, sab, sbb, stot, munging);
                dens[2*i  ] = ra;
                dens[2*i+1] = rb;

                sig[3*i  ] = saa;
                sig[3*i+1] = sab;
                sig[3*i+2] = sbb;

            }
        }
        else {
            MADNESS_EXCEPTION("only LDA and GGA available in xcfunctional",1);
        }
    }
}


madness::Tensor<double> XCfunctional::exc(const std::vector< madness::Tensor<double> >& t) const {
    madness::Tensor<double> rho, sigma;
    make_libxc_args(t, rho, sigma, xc_potential);

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
//                std::cout << "exc: " << j << " " << res[j] << " " << work[j] << " "
//                << dens[j] << " " << sig[j] << " " << funcs[i].first << " " << funcs[i].second << std::endl;
            }
        }
    }
    //std::cout << "result " << result << std::endl;
    return result;
}


madness::Tensor<double> XCfunctional::vxc(const std::vector< madness::Tensor<double> >& t,
        const int ispin, const xc_contribution xc_contrib) const {
    madness::Tensor<double> rho, sigma;
    make_libxc_args(t, rho, sigma, xc_potential);

    // number of grid points
    const int np = t[0].size();

    // number of intermediates depends on the spin
    int nvsig=1, nvrho=1;
    if (spin_polarized) {
        nvrho = 2;
        nvsig = 3;
    }

    madness::Tensor<double> result(3L, t[0].dims());
    double * restrict res = result.ptr();
    const double * restrict dens = rho.ptr();
    for (long j=0; j<np; j++) res[j] = 0.0;

    for (unsigned int i=0; i<funcs.size(); i++) {
        switch(funcs[i].first->info->family) {
        case XC_FAMILY_LDA:
        {
            madness::Tensor<double> vrho(nvrho*np);
            double * restrict vr = vrho.ptr();
            xc_lda_vxc(funcs[i].first, np, dens, vr);

            for (long j=0; j<np; j++) res[j] += vr[nvrho*j+ispin]*funcs[i].second;
        }

        break;

        case XC_FAMILY_HYB_GGA:
        case XC_FAMILY_GGA:
        {
            madness::Tensor<double> vrho(nvrho*np), vsig(nvsig*np);
            double * restrict vr = vrho.ptr();
            double * restrict vs = vsig.ptr();
            const double * restrict sig = sigma.ptr();
            // in: funcs[i].first
            // in: np      number of points
            // in: dens    the density [a,b]
            // in: sig     contracted density gradients \nabla \rho . \nabla \rho [aa,ab,bb]
            // out: vr     \del e/\del \rho_alpha [a,b]
            // out: vs     \del e/\del sigma_alpha [aa,ab,bb]
            xc_gga_vxc(funcs[i].first, np, dens, sig, vr, vs);

            if (spin_polarized) {
                if (xc_contrib == potential_rho) {
                    for (long j=0; j<np; j++) {                 // Vrhoa
                        res[j] += vr[nvrho*j+ispin] * funcs[i].second;
                    }
                }
                else if (xc_contrib == potential_same_spin) {   // Vsigaa/Vsigbb
                    for (long j=0; j<np; j++) {
                        res[j] += vs[nvsig*j + 2*ispin] * funcs[i].second;
                    }
                }
                else if (xc_contrib == potential_mixed_spin) {  // Vsigab
                    for (long j=0; j<np; j++) {
                        res[j] += vs[nvsig*j + 1] * funcs[i].second;
                    }
                }
                else {
                    throw "ouch";
                }
            }
            else {
                if (xc_contrib == potential_rho) {
                    for (long j=0; j<np; j++) {                 // Vrhoa
                        res[j] += vr[j]*funcs[i].second;
                    }
                }
                else if (xc_contrib == potential_same_spin) {   // Vsigaa
                    for (long j=0; j<np; j++) {
                        res[j] += vs[j]*funcs[i].second;
                    }
                }
                else {
                    throw "ouch";
                }
            }
        }
        break;
        default:
            MADNESS_EXCEPTION("unknown XC_FAMILY xcfunctional::vxc",1);
        }
    }
    for (long j=0; j<np; j++) {
        if (isnan_x(res[j])) MADNESS_EXCEPTION("NaN in xcfunctional::vxc",1);
    }

    return result;
}


/// compute the derivative of the XC potential (2nd derivative of the XC energy)

/// @param[in]  t   vector of Tensors holding rho and sigma
/// @param[in]  ispin   the current spin (0=alpha, 1=beta)
/// @param[in]  xc_contrib    which term to compute
Tensor<double> XCfunctional::fxc_apply(const std::vector<Tensor<double> >& t,
        const int ispin, const xc_contribution xc_contrib) const {

    MADNESS_ASSERT(!spin_polarized);    // for now
    MADNESS_ASSERT(ispin==0);           // for now

    // copy quantities from t to rho and sigma
    Tensor<double> rho,sigma;   // rho=2rho_alpha, sigma=4sigma_alpha
    make_libxc_args(t, rho, sigma, xc_potential);

    // number of grid points
    const int np = t[0].size();

    // spin dimensions of the tensors
    const int nspin=(spin_polarized ? 2 : 1);   // rhf: 1; uhf: 2
    const int nspin2=nspin*(nspin+1)/2;         // rhf: 1; uhf: 3
    const int nspin3=nspin2*(nspin2+1)/2;       // rhf: 1; uhf: 6

    // intermediate tensors: partial derivatives of f_xc wrt rho/sigma
    Tensor<double> v2rho2(nspin2*np);       // lda, gga
    Tensor<double> v2rhosigma(nspin3*np);   // gga
    Tensor<double> v2sigma2(nspin3*np);     // gga
    Tensor<double> vrho(nspin*np);          // gga
    Tensor<double> vsigma(nspin2*np);       // gga

    // result tensor
    Tensor<double> result(3L, t[0].dims());

    for (unsigned int i=0; i<funcs.size(); i++) {
        switch(funcs[i].first->info->family) {
        case XC_FAMILY_LDA: {
            double * restrict vr = v2rho2.ptr();
            const double * restrict dens = rho.ptr();
            xc_lda_fxc(funcs[i].first, np, dens, vr);
        }
        break;

        case XC_FAMILY_HYB_GGA:
        case XC_FAMILY_GGA:
        {
            if ((xc_contrib == XCfunctional::kernel_second_semilocal) or
                    (xc_contrib== XCfunctional::kernel_second_local)) {   // partial second derivatives
                double * restrict vrr = v2rho2.ptr();
                double * restrict vrs = v2rhosigma.ptr();
                double * restrict vss = v2sigma2.ptr();
                const double * restrict sig = sigma.ptr();
                const double * restrict dens = rho.ptr();

                // in: funcs[i].first
                // in: np      number of points
                // in: dens    the density [a,b], or 2*\rho_alpha
                // in: sig     contracted density gradients \nabla \rho . \nabla \rho [aa,ab,bb]
                // out: vrr     \del^2 e/\del \rho^2_alpha [a,b]
                // out: vrs     \del^2 e/\del \sigma_alpha\rho [aa,ab,bb]
                // out: vss     \del^2 e/\del \sigma^2_alpha [aa,ab,bb]
                xc_gga_fxc(funcs[i].first, np, dens, sig, vrr, vrs, vss);

            } else if (xc_contrib == XCfunctional::kernel_first_semilocal) {   // partial first derivatives
                double * restrict vr = vrho.ptr();
                double * restrict vs = vsigma.ptr();
                const double * restrict sig = sigma.ptr();
                const double * restrict dens = rho.ptr();

                // in: funcs[i].first
                // in: np      number of points
                // in: dens    the density [a,b]
                // in: sig     contracted density gradients \nabla \rho . \nabla \rho [aa,ab,bb]
                // out: vr     \del e/\del \rho_alpha [a,b]
                // out: vs     \del e/\del sigma_alpha [aa,ab,bb]
                xc_gga_vxc(funcs[i].first, np, dens, sig, vr, vs);

            }
        }
        break;
        default:
            MADNESS_EXCEPTION("unknown XC_FAMILY xcfunctional::fxc",1);
        }

        Tensor<double> result1;

        // multiply the kernel with the various densities
        if (xc_contrib== XCfunctional::kernel_second_local) {  // local terms, second derivative
            Tensor<double> dens_pt=copy(t[enum_rho_pt]);
            Tensor<double> sigma_pt=2.0*copy(t[enum_sigma_pta]);   // factor 2 for closed shell
            munger m(rhotol,rhomin);
            dens_pt.unaryop(m);
            sigma_pt.unaryop(m);

            result1=v2rho2.emul(dens_pt);
            if (is_gga()) result1+= 2.0*v2rhosigma.emul(sigma_pt);

        } else if (xc_contrib== XCfunctional::kernel_second_semilocal) {   // semilocal terms, second derivative
//            const Tensor<double>& dens_pt=t[enum_rho_pt];
//            const Tensor<double>& sigma_pt=2.0*t[enum_sigma_pta];       // factor 2 for closed shell
            Tensor<double> dens_pt=copy(t[enum_rho_pt]);
            Tensor<double> sigma_pt=2.0*copy(t[enum_sigma_pta]);   // factor 2 for closed shell
            munger m(rhotol,rhomin);
            dens_pt.unaryop(m);
            sigma_pt.unaryop(m);


            result1=2.0*v2rhosigma.emul(dens_pt) + 4.0*v2sigma2.emul(sigma_pt);

        } else if (xc_contrib== XCfunctional::kernel_first_semilocal) {   // semilocal terms, first derivative
            result1=2.0*vsigma;
        }

        // accumulate into result tensor with proper weighting
        result+=result1*funcs[i].second;

    }

    // check for NaNs
    double * restrict res = result.ptr();
    for (long j=0; j<np; j++) if (isnan_x(res[j]))
        MADNESS_EXCEPTION("NaN in xcfunctional::fxc_apply",1);
    return result;
}


madness::Tensor<double> XCfunctional::fxc_old(const std::vector< madness::Tensor<double> >& t,
        const int ispin, const xc_contribution what) const {
/* Some useful formulas:
   sigma_st = grad rho_s . grad rho_t
   zk = energy density per unit particle
   vrho_s = d zk / d rho_s
   vsigma_st = d n*zk / d sigma_st
   v2rho2_st = d^2 n*zk / d rho_s d rho_t
   v2rhosigma_svx = d^2 n*zk / d rho_s d sigma_tv
   v2sigma2_stvx = d^2 n*zk / d sigma_st d sigma_vx
if nspin == 2
   rho(2) = (u, d)
   sigma(3) = (uu, du, dd)
   vrho(2) = (u, d)
   vsigma(3) = (uu, du, dd)
   v2rho2(3) = (uu, du, dd)
   v2rhosigma(6) = (u_uu, u_ud, u_dd, d_uu, d_ud, d_dd)
   v2sigma2(6) = (uu_uu, uu_ud, uu_dd, ud_ud, ud_dd, dd_dd)
*/
    madness::Tensor<double> rho, sigma;
    make_libxc_args(t, rho, sigma, xc_potential);

    const int np = t[0].size();

    int nv2rho2, nv2rhosigma, nv2sigma2;
    if (spin_polarized) {
    	nv2rho2 = 3;
    	nv2rhosigma = 6;
    	nv2sigma2 = 6;
    }
    else {
    	nv2rho2 = 1;
    	nv2rhosigma = 1;
    	nv2sigma2 = 1;
    }

    madness::Tensor<double> result(3L, t[0].dims());
    double * restrict res = result.ptr();
    const double * restrict dens = rho.ptr();
    for (long j=0; j<np; j++) 
    			res[j] = 0.0;

    if(what ==99){
// for debugging 
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
                madness::Tensor<double> v2rho2(nv2rho2*np);
    		double * restrict vr = v2rho2.ptr();
    		xc_lda_fxc(funcs[i].first, np, dens, vr);
             	if (what < 2) {
                   //fxc_lda
                     for (long j=0; j<np; j++) {
    			     	res[j] += vr[nv2rho2*j+ispin]*funcs[i].second;
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
    	        madness::Tensor<double> v2rho2(nv2rho2*np);
    	        madness::Tensor<double> v2rhosigma(nv2rhosigma*np);
    	        madness::Tensor<double> v2sigma2(nv2sigma2*np);
                double * restrict v2r2 = v2rho2.ptr();
                double * restrict v2rs = v2rhosigma.ptr();
                double * restrict v2s2 = v2sigma2.ptr();
    	        const double * restrict sig = sigma.ptr();
 //void xc_gga_fxc(xc_func_type *p, int np, double *rho, double *sigma, double *v2rho2, double *v2rhosigma, double *v2sigma2)
	      	xc_gga_fxc(funcs[i].first, np, dens, sig, v2r2, v2rs, v2s2);
	    	if (spin_polarized) {
        		if (what == 0) {
                        // V2rho2_aa 1 or V2rho2_bb
        			for (long j=0; j<np; j++) {
    					res[j] += v2r2[nv2rho2*j    ] * funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else if (what == 1) {
    			// V2rho2_ab 2
    				for (long j=0; j<np; j++) {
    					res[j] += v2r2[nv2rho2*j + 1] * funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else if (what == 3) {
    			// V2rho2_bb 3*
    				for (long j=0; j<np; j++) {
    					res[j] += v2r2[nv2rho2*j + 2] * funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else if (what == 4) {
    			// V2rhosigma_a_aa 1
    				for (long j=0; j<np; j++) {
    					res[j] += v2rs[nv2rhosigma*j    ] * funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else if (what == 5) {
    			// V2rhosigma_a_ab 2
    				for (long j=0; j<np; j++) {
    					res[j] += v2rs[nv2rhosigma*j + 1] * funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else if (what == 6) {
    			// V2rhosigma_a_bb 3
    				for (long j=0; j<np; j++) {
    					res[j] += v2rs[nv2rhosigma*j + 2] * funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else if (what == 7) {
    			// V2rhosigma_b_aa 4
    				for (long j=0; j<np; j++) {
    					res[j] += v2rs[nv2rhosigma*j + 3] * funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else if (what == 8) {
    			// V2rhosigma_b_ab 5
    				for (long j=0; j<np; j++) {
    					res[j] += v2rs[nv2rhosigma*j + 4] * funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else if (what == 9) {
    			// V2rhosigma_b_bb 6 **
    				for (long j=0; j<np; j++) {
    					res[j] += v2rs[nv2rhosigma*j + 5] * funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else if (what == 10) {
    			// V2sigma2_aa_aa 1
    				for (long j=0; j<np; j++) {
    					res[j] += v2s2[nv2sigma2*j    ] * funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else if (what == 11) {
    			// V2sigma2_aa_ab 2
    				for (long j=0; j<np; j++) {
    					res[j] += v2s2[nv2sigma2*j + 1] * funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else if (what == 12) {
    			// V2sigma2_aa_bb 3
    				for (long j=0; j<np; j++) {
    					res[j] += v2s2[nv2sigma2*j + 2] * funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else if (what == 13) {
    			// V2sigma2_ab_ab 4
    				for (long j=0; j<np; j++) {
    					res[j] += v2s2[nv2sigma2*j + 3] * funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else if (what == 14) {
    			// V2sigma2_ab_bb 5
    				for (long j=0; j<np; j++) {
    					res[j] += v2s2[nv2sigma2*j + 4] * funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else if (what == 15) {
    			// V2sigma2_bb_bb 6
    				for (long j=0; j<np; j++) {
    					res[j] += v2s2[nv2sigma2*j + 5] * funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch";
    				}
    			}
    			else {
    				throw "ouch";
    			}
    		}
    		else {
        		if (what == 0) {
                        // V2rho2_aa 1
        			for (long j=0; j<np; j++) {
    					res[j] += v2r2[nv2rho2*j] * funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch nan V2rho2_aa";
    				}
    			}
    			else if (what == 1) {
    			// V2rhosigma_a_aa 1
    				for (long j=0; j<np; j++) {
    					res[j] += v2rs[nv2rhosigma*j] * funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch nan V2rhosigma_a_aa";
    				}
    			}
    			else if (what == 2) {
    			// V2sigma2_aa_aa 1
    				for (long j=0; j<np; j++) {
    					res[j] += v2s2[nv2sigma2*j] * funcs[i].second;
    					if (isnan_x(res[j])) throw "ouch nan V2sigma2_aa_aa";
    				}
    			}
    		}
   		}
    		break;
         	case XC_FAMILY_HYB_GGA:
        	{
       				throw "ouch XC_FAMILY_HYB_GGA fxc disabled" ;
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
