#include <madness/madness_config.h>
#include<madness/chem/xcfunctional.h>
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

static int lookup_name(const std::string& name) {
    // Call libxc routine
    //return XC(functional_get_number(name.c_str()));
    return xc_functional_get_number(name.c_str());
}

static std::string lookup_id(const int id) {
    // Call libxc routine, needs memory handling
    //char *namep(XC(functional_get_name(id)));
    //std::string name = (namep==NULL) ? "Functional not found" : std::string(namep);
    //free(namep);
    //return name;
    return std::string(xc_functional_get_name(id));
}

static xc_func_type* make_func(int id, bool polarized) {
    xc_func_type* func = new xc_func_type;
    int POLARIZED = polarized ? XC_POLARIZED : XC_UNPOLARIZED;
    //MADNESS_ASSERT(xc_func_init(func, id, POLARIZED) == 0); // SHOULD BE CHECK
    if (xc_func_init(func, id, POLARIZED) != 0) throw "bad stuff!!!!!!!!!!";
    return func;
}

static xc_func_type* lookup_func(const std::string& name, bool polarized) {
    int id = lookup_name(name);
    //MADNESS_ASSERT(id > 0); // SHOULD BE CHECK
    if(id <= 0) throw "bad stuff xxx";
    return make_func(id, polarized);
}

//XCfunctional::XCfunctional() {}
//XCfunctional::XCfunctional() : hf_coeff(0.0) {std::printf("Construct XC Functional from LIBXC Library");}
XCfunctional::XCfunctional() : hf_coeff(0.0) {
    rhotol=1e-7; rhomin=0.0;
    ggatol=1.e-4;
    nderiv=0;
    spin_polarized=false;
}

void XCfunctional::initialize(const std::string& input_line, bool polarized,
        World& world, const bool verbose) {
    rhotol=1e-7; rhomin=0.0; // default values
    ggatol=1.e-4;
    tautol=1.e-12;

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
        } else if (name == "TPSS") {
            // meta-gga: bring these up before r2scan/scan, since TPSS is built on
            // z = tau_W/tau and does not carry the oscillation in grad(de/dtau) at
            // alpha = 1 that makes the SCAN family grid-hungry
            funcs.push_back(std::make_pair(lookup_func("MGGA_X_TPSS",polarized),1.0));
            funcs.push_back(std::make_pair(lookup_func("MGGA_C_TPSS",polarized),1.0));
        } else if (name == "R2SCAN") {
            funcs.push_back(std::make_pair(lookup_func("MGGA_X_R2SCAN",polarized),1.0));
            funcs.push_back(std::make_pair(lookup_func("MGGA_C_R2SCAN",polarized),1.0));
        } else if (name == "SCAN") {
            funcs.push_back(std::make_pair(lookup_func("MGGA_X_SCAN",polarized),1.0));
            funcs.push_back(std::make_pair(lookup_func("MGGA_C_SCAN",polarized),1.0));
        } else if (name == "B3LYP") {
            // VWN-3 correlation; the 0.2 exact exchange comes from the libxc query below
            funcs.push_back(std::make_pair(lookup_func("HYB_GGA_XC_B3LYP",polarized),1.0));
        } else if (name == "RHOMIN") {
            line >> rhomin;
        } else if (name == "RHOTOL") {
            line >> rhotol;
        } else if (name == "GGATOL") {
            line >> ggatol;
        } else if (name == "TAUTOL") {
            line >> tautol;
        } else if (name == "HF" || name == "HF_X") {
            if (! (line >> factor)) factor = 1.0;
            hf_coeff = factor;
        } else {
            if (! (line >> factor)) factor = 1.0;
            funcs.push_back(std::make_pair(lookup_func(name,polarized), factor));
        }
    }

    // nderiv is a rung index (0 lda, 1 gga, 2 meta-gga), not a count of density
    // derivatives -- a meta-gga needs FIRST derivatives of the orbitals (through
    // tau), not second derivatives of the density. Query it through
    // is_lda()/is_gga()/is_meta() and the needs_sigma()/needs_tau() capabilities.
    for (unsigned int i=0; i<funcs.size(); i++) {
        const int family = funcs[i].first->info->family;
        if (family == XC_FAMILY_GGA) nderiv = std::max(nderiv,1);
        if (family == XC_FAMILY_HYB_GGA) nderiv = std::max(nderiv,1);
        if (family == XC_FAMILY_MGGA) nderiv = std::max(nderiv,2);
        if (family == XC_FAMILY_HYB_MGGA) nderiv = std::max(nderiv,2);
 //       if (family == XC_FAMILY_LDA) nderiv = std::max(nderiv,0);
    }

    // Laplacian-level functionals would need nabla^2 rho as a separate input and
    // contribute +nabla^2(de/dnabla^2 rho) to the potential -- four net derivatives
    // of the density. Not supported; lapl is passed to libxc as NULL, so a
    // functional that needs it must be refused rather than silently given zeros.
    for (unsigned int i=0; i<funcs.size(); i++) {
        if (funcs[i].first->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) {
            MADNESS_EXCEPTION("Laplacian-level functionals are not supported: "
                              "no nabla^2 rho intermediate is available",1);
        }
    }

    // The exact-exchange admixture of a hybrid is a property of the functional,
    // not of the input line, so query it from libxc instead of hardcoding it.
    // libxc has corrected several of these coefficients over the years (CAP0 had
    // 75% instead of 25%, LC_BLYP used omega=0.3 instead of 0.33, ...), and a
    // hardcoded value silently diverges from the functional it claims to be.
    // Without this, a hybrid requested by its libxc name -- as opposed to one of
    // the aliases above -- ran with no exact exchange at all.
    for (unsigned int i=0; i<funcs.size(); i++) {
        const int family = funcs[i].first->info->family;
        const bool is_hybrid = (family == XC_FAMILY_HYB_LDA)
                            or (family == XC_FAMILY_HYB_GGA)
                            or (family == XC_FAMILY_HYB_MGGA);
        if (not is_hybrid) continue;

        // range-separated hybrids need an attenuated exchange operator, which the
        // exchange operator does not provide; treating them as global hybrids
        // would silently give the wrong answer.
        const int flags = funcs[i].first->info->flags;
        if ((flags & XC_FLAGS_HYB_CAM) or (flags & XC_FLAGS_HYB_CAMY)) {
            MADNESS_EXCEPTION("range-separated hybrids are not supported: no "
                              "attenuated exchange operator is available",1);
        }
        hf_coeff += funcs[i].second * xc_hyb_exx_coef(funcs[i].first);
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
        print("         ggatol",ggatol);
        if (needs_tau()) print("         tautol",tautol);
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

bool XCfunctional::needs_sigma() const {
    return nderiv >= 1;
}

bool XCfunctional::needs_tau() const {
    return nderiv >= 2;
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


void XCfunctional::make_libxc_args(const std::vector< madness::Tensor<double> >& xc_args,
           madness::Tensor<double>& rho, madness::Tensor<double>& sigma,
           madness::Tensor<double>& tau,
           madness::Tensor<double>& rho_pt, madness::Tensor<double>& sigma_pt,
           std::vector<madness::Tensor<double> >& drho,
           std::vector<madness::Tensor<double> >& drho_pt,
           bool need_response) const {

    // number of grid points in this box
    const int np = xc_args[0].size();


    if (not spin_polarized) {
        if (is_lda()) {
            rho  = madness::Tensor<double>(np);
            const double * MADNESS_RESTRICT rhoa = xc_args[enum_rhoa].ptr();
            double * MADNESS_RESTRICT dens = rho.ptr();
            for (long i=0; i<np; i++) {
                dens[i] = munge(2.0*rhoa[i]);   // full dens is twice alpha dens
            }

            // add perturbed density in response calculations
            // note rho_pt does not depend on the spin
            if (need_response) {
                rho_pt  = madness::Tensor<double>(np);
                const double * MADNESS_RESTRICT rho_pt1 = xc_args[enum_rho_pt].ptr();
                double * MADNESS_RESTRICT dens_pt = rho_pt.ptr();
                for (long i=0; i<np; i++) {
                    dens_pt[i] = binary_munge(rho_pt1[i],rhoa[i],rhotol); // no factor 2
                }
            }
        }
        else if (needs_sigma()) {
            // rho is the density
            // the reduced density gradient sigma is given by
            // sigma = rho * rho * chi
            const double * MADNESS_RESTRICT rhoa = xc_args[enum_rhoa].ptr();
            const double * MADNESS_RESTRICT chiaa = xc_args[enum_chi_aa].ptr();
            const double * MADNESS_RESTRICT zetaa_x = xc_args[enum_zetaa_x].ptr();
            const double * MADNESS_RESTRICT zetaa_y = xc_args[enum_zetaa_y].ptr();
            const double * MADNESS_RESTRICT zetaa_z = xc_args[enum_zetaa_z].ptr();

            // output
            rho  = madness::Tensor<double>(np);
            drho[0]  = madness::Tensor<double>(np);
            drho[1]  = madness::Tensor<double>(np);
            drho[2]  = madness::Tensor<double>(np);
            sigma  = madness::Tensor<double>(np);

            double * MADNESS_RESTRICT dens = rho.ptr();
            double * MADNESS_RESTRICT sig = sigma.ptr();
            double * MADNESS_RESTRICT ddensx = drho[0].ptr();
            double * MADNESS_RESTRICT ddensy = drho[1].ptr();
            double * MADNESS_RESTRICT ddensz = drho[2].ptr();

            for (long i=0; i<np; i++) {
                dens[i]=munge(2.0*rhoa[i]);     // full dens is twice alpha dens
                sig[i] = std::max(1.e-14,dens[i]*dens[i] * chiaa[i]);   // 2 factors 2 included in dens
                ddensx[i]=dens[i]*zetaa_x[i];
                ddensy[i]=dens[i]*zetaa_y[i];
                ddensz[i]=dens[i]*zetaa_z[i];
            }

            if (needs_tau()) {
                const double * MADNESS_RESTRICT taua = xc_args[enum_taua].ptr();
                madness::Tensor<double> taudummy;
                if (taua==NULL) {           // caller did not provide tau
                    taudummy=madness::Tensor<double>(np);
                    taua=taudummy.ptr();
                }
                tau = madness::Tensor<double>(np);
                double * MADNESS_RESTRICT t = tau.ptr();
                for (long i=0; i<np; i++) {
                    double ti = std::max(tautol,2.0*taua[i]);   // full tau is twice alpha tau
                    // tau >= tau_W = sigma/(8 rho) is exact, so the Fermi hole
                    // curvature stays positive and the iso-orbital indicators stay in
                    // range. Bounding tau from below rather than clamping sigma down
                    // (which is what libxc's XC_FLAGS_ENFORCE_FHC does) leaves the
                    // density gradient untouched.
                    if (dens[i] > rhotol) ti = std::max(ti,sig[i]/(8.0*dens[i]));
                    t[i] = ti;
                }
            }

            // add perturbed density and density gradients in response calculations
            // note rho_pt does not depend on the spin
            if (need_response) {

                // input
                const double * MADNESS_RESTRICT rho_pt1 = xc_args[enum_rho_pt].ptr();
                const double * MADNESS_RESTRICT sig_pt1 = xc_args[enum_sigma_pta_div_rho].ptr();
                const double * MADNESS_RESTRICT drho_ptx1 = xc_args[enum_ddens_ptx].ptr();
                const double * MADNESS_RESTRICT drho_pty1 = xc_args[enum_ddens_pty].ptr();
                const double * MADNESS_RESTRICT drho_ptz1 = xc_args[enum_ddens_ptz].ptr();

                // output
                rho_pt  = madness::Tensor<double>(np);
                sigma_pt  = madness::Tensor<double>(np);
                drho_pt[0]  = madness::Tensor<double>(np);
                drho_pt[1]  = madness::Tensor<double>(np);
                drho_pt[2]  = madness::Tensor<double>(np);

                double * MADNESS_RESTRICT ddens_ptx = drho_pt[0].ptr();
                double * MADNESS_RESTRICT ddens_pty = drho_pt[1].ptr();
                double * MADNESS_RESTRICT ddens_ptz = drho_pt[2].ptr();
                double * MADNESS_RESTRICT dens_pt = rho_pt.ptr();
                double * MADNESS_RESTRICT sig_pt = sigma_pt.ptr();

                for (long i=0; i<np; i++) {
                    dens_pt[i] = binary_munge(rho_pt1[i],rhoa[i],rhotol);  // no factor 2
                    sig_pt[i] = dens[i]*sig_pt1[i];
                    ddens_ptx[i] = binary_munge(drho_ptx1[i],rhoa[i],rhotol);  // no factor 2
                    ddens_pty[i] = binary_munge(drho_pty1[i],rhoa[i],rhotol);  // no factor 2
                    ddens_ptz[i] = binary_munge(drho_ptz1[i],rhoa[i],rhotol);  // no factor 2
                    // dens is munged and includes factor of 2 for full density
                }
            }

        }
        else {
            MADNESS_EXCEPTION("unknown functional rung in xcfunctional",1);
        }

    } else if (spin_polarized) {
        if (is_lda()) {
            const double * MADNESS_RESTRICT rhoa = xc_args[enum_rhoa].ptr();
            const double * MADNESS_RESTRICT rhob = xc_args[enum_rhob].ptr();
            rho  = madness::Tensor<double>(np*2L);
            double * MADNESS_RESTRICT dens = rho.ptr();

            // might happen if there are no beta electrons
            madness::Tensor<double> dummy;
            if (rhob==NULL) {
                dummy=madness::Tensor<double>(np);
                rhob=dummy.ptr();
            }

            for (long i=0; i<np; i++) {
                dens[2*i  ] = munge(rhoa[i]);
                dens[2*i+1] = munge(rhob[i]);
            }
            if (need_response) {
                MADNESS_EXCEPTION("no spin polarized DFT response in xcfunctional",1);
            }
        }
        else if (needs_sigma()) {
            // input
            const double * MADNESS_RESTRICT rhoa  = xc_args[enum_rhoa].ptr();
            const double * MADNESS_RESTRICT rhob  = xc_args[enum_rhob].ptr();

            const double * MADNESS_RESTRICT chiaa = xc_args[enum_chi_aa].ptr();
            const double * MADNESS_RESTRICT chiab = xc_args[enum_chi_ab].ptr();
            const double * MADNESS_RESTRICT chibb = xc_args[enum_chi_bb].ptr();

            const double * MADNESS_RESTRICT zetaa_x = xc_args[enum_zetaa_x].ptr();
            const double * MADNESS_RESTRICT zetaa_y = xc_args[enum_zetaa_y].ptr();
            const double * MADNESS_RESTRICT zetaa_z = xc_args[enum_zetaa_z].ptr();

            const double * MADNESS_RESTRICT zetab_x = xc_args[enum_zetab_x].ptr();
            const double * MADNESS_RESTRICT zetab_y = xc_args[enum_zetab_y].ptr();
            const double * MADNESS_RESTRICT zetab_z = xc_args[enum_zetab_z].ptr();

            // might happen if there are no beta electrons: prep_xc_args leaves the
            // beta intermediates unassigned, so their data pointers are NULL.
            // Substitute a zero tensor -- rho_beta is zero as well, so all beta
            // contributions vanish regardless of what zeta_beta is set to.
            madness::Tensor<double> dummy;
            if ((rhob==NULL) or (chiab==NULL) or (chibb==NULL)
                    or (zetab_x==NULL) or (zetab_y==NULL) or (zetab_z==NULL)) {
                dummy=madness::Tensor<double>(np);
            }
            if (rhob==NULL) rhob=dummy.ptr();
            if (chiab==NULL) chiab=dummy.ptr();
            if (chibb==NULL) chibb=dummy.ptr();
            if (zetab_x==NULL) zetab_x=dummy.ptr();
            if (zetab_y==NULL) zetab_y=dummy.ptr();
            if (zetab_z==NULL) zetab_z=dummy.ptr();

            rho   = madness::Tensor<double>(np*2L);
            drho[0]  = madness::Tensor<double>(np*2L);
            drho[1]  = madness::Tensor<double>(np*2L);
            drho[2]  = madness::Tensor<double>(np*2L);
            sigma = madness::Tensor<double>(np*3L);

            double * MADNESS_RESTRICT dens = rho.ptr();
            double * MADNESS_RESTRICT sig  = sigma.ptr();
            double * MADNESS_RESTRICT ddensx  = drho[0].ptr();
            double * MADNESS_RESTRICT ddensy  = drho[1].ptr();
            double * MADNESS_RESTRICT ddensz  = drho[2].ptr();


            for (long i=0; i<np; i++) {

                double ra=munge(rhoa[i]);
                double rb=munge(rhob[i]);

                dens[2*i  ] = ra;
                dens[2*i+1] = rb;
                sig[3*i  ]  = std::max(1.e-14,ra * ra * chiaa[i]);  // aa
                // sigma_ab = grad(rho_a).grad(rho_b) is bilinear and may legitimately
                // be negative; only its magnitude is bounded (by Cauchy-Schwarz on the
                // zeta vectors, which holds automatically here). Do not clamp the sign.
                sig[3*i+1]  = ra * rb * chiab[i];                   // ab
                sig[3*i+2]  = std::max(1.e-14,rb * rb * chibb[i]);  // bb

                ddensx[2*i  ]=ra * zetaa_x[i];
                ddensx[2*i+1]=rb * zetab_x[i];
                ddensy[2*i  ]=ra * zetaa_y[i];
                ddensy[2*i+1]=rb * zetab_y[i];
                ddensz[2*i  ]=ra * zetaa_z[i];
                ddensz[2*i+1]=rb * zetab_z[i];


            }

            if (needs_tau()) {
                const double * MADNESS_RESTRICT taua = xc_args[enum_taua].ptr();
                const double * MADNESS_RESTRICT taub = xc_args[enum_taub].ptr();
                madness::Tensor<double> taudummy;
                if ((taua==NULL) or (taub==NULL)) {
                    taudummy=madness::Tensor<double>(np);
                }
                if (taua==NULL) taua=taudummy.ptr();
                if (taub==NULL) taub=taudummy.ptr();

                tau = madness::Tensor<double>(np*2L);
                double * MADNESS_RESTRICT t = tau.ptr();
                for (long i=0; i<np; i++) {
                    double ta = std::max(tautol,taua[i]);
                    double tb = std::max(tautol,taub[i]);
                    // tau_s >= sigma_ss/(8 rho_s) in each spin channel, see above
                    if (dens[2*i  ] > rhotol) ta = std::max(ta,sig[3*i  ]/(8.0*dens[2*i  ]));
                    if (dens[2*i+1] > rhotol) tb = std::max(tb,sig[3*i+2]/(8.0*dens[2*i+1]));
                    t[2*i  ] = ta;
                    t[2*i+1] = tb;
                }
            }

            if (need_response) {
                MADNESS_EXCEPTION("no spin polarized DFT response in xcfunctional",1);
            }
        }
        else {
            MADNESS_EXCEPTION("unknown functional rung in xcfunctional",1);
        }
    }
}


madness::Tensor<double> XCfunctional::exc(const std::vector< madness::Tensor<double> >& t) const {
    madness::Tensor<double> rho, sigma, tau, rho_pt, sigma_pt;
    std::vector<Tensor<double> > ddens(3), ddens_pt(3);
    make_libxc_args(t, rho, sigma, tau, rho_pt, sigma_pt, ddens, ddens_pt, false);

    const int np = t[0].size();
    const double * MADNESS_RESTRICT dens = rho.ptr();
    const double * MADNESS_RESTRICT sig = sigma.ptr();
    const double * MADNESS_RESTRICT ktau = tau.ptr();
    // zeroed laplacian input for meta-ggas (see the MGGA case below)
    madness::Tensor<double> lapl_zero;
    if (needs_tau()) lapl_zero=madness::Tensor<double>(rho.size());

    madness::Tensor<double> result(3L, t[0].dims());
    double * MADNESS_RESTRICT res = result.ptr();
    for (long j=0; j<np; j++) res[j] = 0.0;

    for (unsigned int i=0; i<funcs.size(); i++) {
        madness::Tensor<double> zk(3L, t[0].dims(), false);
        double * MADNESS_RESTRICT work = zk.ptr();

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
        case XC_FAMILY_MGGA:
        case XC_FAMILY_HYB_MGGA:
            // Laplacian-level functionals are refused in initialize(), so the
            // laplacian input is irrelevant here -- but pass a real zeroed buffer
            // rather than NULL, because libxc's generated kernels do not guard
            // every laplacian access on the pointer
            xc_mgga_exc(funcs[i].first, np, dens, sig, lapl_zero.ptr(), ktau, work);
            break;
        default:
            MADNESS_EXCEPTION("unknown XC_FAMILY in xcfunctional::exc",1);
        }
        if (spin_polarized) {
            for (long j=0; j<np; j++) {
                res[j] +=  work[j]*(dens[2*j+1] + dens[2*j])*funcs[i].second;
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


std::vector<madness::Tensor<double> > XCfunctional::vxc(
        const std::vector< madness::Tensor<double> >& t, const int ispin) const {
    madness::Tensor<double> rho, sigma, tau, dummy;
    std::vector<Tensor<double> > drho(3), ddens_pt(3);
    make_libxc_args(t, rho, sigma, tau, dummy, dummy, drho, ddens_pt, false);

    // number of grid points
    const int np = t[0].size();

    // number of intermediates depends on the spin

//    dens[2*i  ] = ra;
//    dens[2*i+1] = rb;
//    sig[3*i  ]  = std::max(1.e-14,ra * ra * chiaa[i]);  // aa
//    sig[3*i+1]  = std::max(1.e-14,ra * rb * chiab[i]);  // ab
//    sig[3*i+2]  = std::max(1.e-14,rb * rb * chibb[i]);  // bb

    int nvsig=1, nvrho=1;
    if (spin_polarized) {
        nvrho = 2;
        nvsig = 3;
    }

    // must agree with xc_potential::get_result_size(), which is what the caller
    // allocated -- multi_to_multi_op_values only cross-checks the two under
    // MADNESS_ASSERT, which is compiled out in a release build
    const std::size_t result_size=xc_potential(*this,ispin).get_result_size();
    MADNESS_CHECK(result_size>0);

    Tensor<double> r(3L, t[0].dims());
    r=0.0;
    std::vector<Tensor<double> > result(result_size);
    for (Tensor<double>& rr : result) rr=copy(r);

    const double * MADNESS_RESTRICT dens = rho.ptr();   // nspin * np
    const double * MADNESS_RESTRICT ddensx = drho[0].ptr();  // nspin * np
    const double * MADNESS_RESTRICT ddensy = drho[1].ptr();  // nspin * np
    const double * MADNESS_RESTRICT ddensz = drho[2].ptr();  // nspin * np

    for (unsigned int i=0; i<funcs.size(); i++) {
        switch(funcs[i].first->info->family) {
        case XC_FAMILY_LDA:
        {
            madness::Tensor<double> vrho(nvrho*np);
            double * MADNESS_RESTRICT vr = vrho.ptr();
            xc_lda_vxc(funcs[i].first, np, dens, vr);
            double * MADNESS_RESTRICT r0 = result[0].ptr();

            for (long j=0; j<np; j++) r0[j] += vr[nvrho*j+ispin]*funcs[i].second;
        }

        break;

        case XC_FAMILY_HYB_MGGA:
        case XC_FAMILY_MGGA:
        case XC_FAMILY_HYB_GGA:
        case XC_FAMILY_GGA:
        {
            const int family=funcs[i].first->info->family;
            const bool meta=(family==XC_FAMILY_MGGA) or (family==XC_FAMILY_HYB_MGGA);

            madness::Tensor<double> vrho(nvrho*np), vsig(nvsig*np), vtau;
            double * MADNESS_RESTRICT vr = vrho.ptr();
            double * MADNESS_RESTRICT vs = vsig.ptr();
            double * MADNESS_RESTRICT vt = NULL;
            if (meta) {
                vtau=madness::Tensor<double>(nvrho*np);
                vt=vtau.ptr();
            }
            const double * MADNESS_RESTRICT sig = sigma.ptr();
            // in: funcs[i].first
            // in: np      number of points
            // in: dens    the density [a,b]
            // in: sig     contracted density gradients \nabla \rho . \nabla \rho [aa,ab,bb]
            // in: ktau    the kinetic energy density [a,b] (meta-gga only)
            // out: vr     \del e/\del \rho_alpha [a,b]
            // out: vs     \del e/\del sigma_alpha [aa,ab,bb]
            // out: vt     \del e/\del \tau_alpha [a,b] (meta-gga only)
            if (meta) {
                // Laplacian-level functionals are refused in initialize(), so both
                // the laplacian input and the de/dlaplacian output are irrelevant --
                // but they get real buffers rather than NULL, because libxc's
                // generated kernels do not guard every laplacian access on the pointer
                madness::Tensor<double> lapl_zero(nvrho*np), vlapl(nvrho*np);
                xc_mgga_vxc(funcs[i].first, np, dens, sig, lapl_zero.ptr(), tau.ptr(),
                            vr, vs, vlapl.ptr(), vt);
            } else {
                xc_gga_vxc(funcs[i].first, np, dens, sig, vr, vs);
            }

            if (spin_polarized) {
                double * MADNESS_RESTRICT r0 = result[0].ptr();
                double * MADNESS_RESTRICT r1 = result[1].ptr();
                double * MADNESS_RESTRICT r2 = result[2].ptr();
                double * MADNESS_RESTRICT r3 = result[3].ptr();
                double * MADNESS_RESTRICT r4 = result[4].ptr();
                double * MADNESS_RESTRICT r5 = result[5].ptr();
                double * MADNESS_RESTRICT r6 = result[6].ptr();

                for (long j=0; j<np; j++) {
                    // Vrhoa
                    r0[j] += vr[nvrho*j+ispin] * funcs[i].second;

                    // Vsigaa/Vsigbb * rho
                    r1[j] += 2.0 * vs[nvsig*j + 2*ispin] * funcs[i].second       // aa or bb in steps of 3
                            *ddensx[nvrho*j + ispin];                         // a or b in steps of 2
                    r2[j] += 2.0 * vs[nvsig*j + 2*ispin] * funcs[i].second       // aa or bb in steps of 3
                            *ddensy[nvrho*j + ispin];                         // a or b in steps of 2
                    r3[j] += 2.0 * vs[nvsig*j + 2*ispin] * funcs[i].second       // aa or bb in steps of 3
                            *ddensz[nvrho*j + ispin];                         // a or b in steps of 2

                    // Vsigab * rho_other_spin
                    r4[j] += vs[nvsig*j + 1] * funcs[i].second             // ab in steps of 3
                            *ddensx[nvrho*j + (1-ispin)];                     // b or a in steps of 2
                    r5[j] += vs[nvsig*j + 1] * funcs[i].second             // ab in steps of 3
                            *ddensy[nvrho*j + (1-ispin)];                     // b or a in steps of 2
                    r6[j] += vs[nvsig*j + 1] * funcs[i].second             // ab in steps of 3
                            *ddensz[nvrho*j + (1-ispin)];                     // b or a in steps of 2

                }
            }
            else {
                double * MADNESS_RESTRICT r0 = result[0].ptr();
                double * MADNESS_RESTRICT r1 = result[1].ptr();
                double * MADNESS_RESTRICT r2 = result[2].ptr();
                double * MADNESS_RESTRICT r3 = result[3].ptr();

                for (long j=0; j<np; j++) {
                    // Vrhoa
                    r0[j] += vr[j]*funcs[i].second;

                    // Vsigaa
                    r1[j] += 2.0 * vs[j]*funcs[i].second*ddensx[j];    // total density
                    r2[j] += 2.0 * vs[j]*funcs[i].second*ddensy[j];    // total density
                    r3[j] += 2.0 * vs[j]*funcs[i].second*ddensz[j];    // total density
                }
            }

            // de/dtau of the spin this operator acts on; the caller turns it into
            // the non-multiplicative operator -1/2 nabla.(de/dtau nabla psi)
            if (meta) {
                double * MADNESS_RESTRICT rt = result[spin_polarized ? 7 : 4].ptr();
                // Unlike the semilocal flux terms, which carry a factor rho*zeta and
                // are damped in the tail on their own, de/dtau multiplies grad(psi).
                // That does not vanish nearly as fast, and de/dtau itself diverges
                // where the density is negligible -- it reaches O(100) on the atomic
                // initial guess. Screen it on the density, as the response kernel does.
                for (long j=0; j<np; j++)
                    rt[j] += binary_munge(vt[nvrho*j+ispin]*funcs[i].second,
                                          dens[nvrho*j+ispin],ggatol);
            }
        }
        break;
        default:
            MADNESS_EXCEPTION("unknown XC_FAMILY xcfunctional::vxc",1);
        }
    }

    // check for NaNs
    for (Tensor<double>& rr : result) {
        double * MADNESS_RESTRICT res = rr.ptr();
        for (long j=0; j<np; j++) {
            if (isnan_x(res[j])) MADNESS_EXCEPTION("NaN in xcfunctional::vxc",1);
        }
    }

    return result;
}


std::vector<madness::Tensor<double> > XCfunctional::fxc_apply(
        const std::vector<Tensor<double> >& t, const int ispin) const {

    //MADNESS_CHECK(!spin_polarized);    // for now
    //MADNESS_CHECK(ispin==0);           // for now
    if (spin_polarized || ispin!=0) throw "bad stuff yyy";

    // copy quantities from t to rho and sigma
    Tensor<double> rho,sigma, rho_pt, sigma_pt;   // rho=2rho_alpha, sigma=4sigma_alpha
    std::vector<Tensor<double> > drho(3), drho_pt(3);
    madness::Tensor<double> tau;
    make_libxc_args(t, rho, sigma, tau, rho_pt, sigma_pt, drho, drho_pt, true);

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
    Tensor<double> r(3L, t[0].dims());
    int result_size= this->is_gga() ? 4 : 1;
    std::vector<Tensor<double> > result(result_size);
    for (Tensor<double>& rr : result) rr=copy(r);

    for (unsigned int i=0; i<funcs.size(); i++) {
        switch(funcs[i].first->info->family) {
        case XC_FAMILY_LDA: {
            double * MADNESS_RESTRICT vr = v2rho2.ptr();
            const double * MADNESS_RESTRICT dens = rho.ptr();
            xc_lda_fxc(funcs[i].first, np, dens, vr);

            // only local terms
            result[0]+=funcs[i].second*v2rho2.emul(rho_pt);

        }
        break;

        case XC_FAMILY_HYB_GGA:
        case XC_FAMILY_GGA:
        {
            const double * MADNESS_RESTRICT sig = sigma.ptr();
            const double * MADNESS_RESTRICT dens = rho.ptr();
            const double * MADNESS_RESTRICT sig_pt = sigma_pt.ptr();
            const double * MADNESS_RESTRICT dens_pt = rho_pt.ptr();
            const double * MADNESS_RESTRICT ddens_ptx = drho_pt[0].ptr();
            const double * MADNESS_RESTRICT ddens_pty = drho_pt[1].ptr();
            const double * MADNESS_RESTRICT ddens_ptz = drho_pt[2].ptr();
            const double * MADNESS_RESTRICT ddensx = drho[0].ptr();
            const double * MADNESS_RESTRICT ddensy = drho[1].ptr();
            const double * MADNESS_RESTRICT ddensz = drho[2].ptr();

            double * MADNESS_RESTRICT vr = vrho.ptr();
            double * MADNESS_RESTRICT vs = vsigma.ptr();
            double * MADNESS_RESTRICT vrr = v2rho2.ptr();
            double * MADNESS_RESTRICT vrs = v2rhosigma.ptr();
            double * MADNESS_RESTRICT vss = v2sigma2.ptr();

            double * MADNESS_RESTRICT r0 = result[0].ptr();
            double * MADNESS_RESTRICT r1 = result[1].ptr();
            double * MADNESS_RESTRICT r2 = result[2].ptr();
            double * MADNESS_RESTRICT r3 = result[3].ptr();

            // in: funcs[i].first
            // in: np      number of points
            // in: dens    the density [a,b], or 2*\rho_alpha
            // in: sig     contracted density gradients \nabla \rho . \nabla \rho [aa,ab,bb]
            // out: vrr     \del^2 e/\del \rho^2_alpha [a,b]
            // out: vrs     \del^2 e/\del \sigma_alpha\rho [aa,ab,bb]
            // out: vss     \del^2 e/\del \sigma^2_alpha [aa,ab,bb]
            xc_gga_fxc(funcs[i].first, np, dens, sig, vrr, vrs, vss);

            // in: funcs[i].first
            // in: np      number of points
            // in: dens    the density [a,b]
            // in: sig     contracted density gradients \nabla \rho . \nabla \rho [aa,ab,bb]
            // out: vr     \del e/\del \rho_alpha [a,b]
            // out: vs     \del e/\del sigma_alpha [aa,ab,bb]
            xc_gga_vxc(funcs[i].first, np, dens, sig, vr, vs);


            const double w=funcs[i].second;
            for (long j=0; j<np; j++) {

                // local terms
                r0[j]+=w*(vrr[j]*dens_pt[j] + 2.0*vrs[j] * sig_pt[j]);

                // semilocal terms -- x,y,z
                r1[j]+= w*binary_munge(
                          2.0*vrs[j] * dens_pt[j] * ddensx[j]
                          + 4.0 * vss[j] * sig_pt[j] * ddensx[j]
                          + 2.0 * vs[j]*ddens_ptx[j],
                        dens[j],ggatol);

                r2[j]+=w*binary_munge(
                        2.0*vrs[j] * dens_pt[j] * ddensy[j]
                        + 4.0 * vss[j] * sig_pt[j] * ddensy[j]
                        + 2.0 * vs[j]*ddens_pty[j],
                        dens[j],ggatol);

                r3[j]+=w*binary_munge(
                        2.0*vrs[j] * dens_pt[j] * ddensz[j]
                        + 4.0 * vss[j] * sig_pt[j] * ddensz[j]
                        + 2.0 * vs[j]*ddens_ptz[j],
                        dens[j],ggatol);

            }
        }
        break;
        default:
            MADNESS_EXCEPTION("unknown XC_FAMILY xcfunctional::fxc",1);
        }

        // the functional weight is already folded into each contribution above
        for (std::size_t j=0; j<result.size(); ++j) {

            // check for NaNs
            double * MADNESS_RESTRICT res = result[j].ptr();
            for (long jj=0; jj<np; jj++) if (isnan_x(res[jj]))
                MADNESS_EXCEPTION("NaN in xcfunctional::fxc_apply",1);

        }
    }

    return result;
}

}
