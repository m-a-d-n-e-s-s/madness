#include <madness/madness_config.h>
#include<madness/chem/xcfunctional.h>
#include <madness/tensor/tensor.h>
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
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
    if (xc_func_init(func, id, POLARIZED) != 0) {
        delete func;
        // runtime_error, not MADNESS_EXCEPTION: MadnessException stores the
        // char* without copying, so a message built on the fly would dangle
        throw std::runtime_error("libxc could not initialize functional id "
                                 + std::to_string(id));
    }
    return func;
}

static xc_func_type* lookup_func(const std::string& name, bool polarized) {
    int id = lookup_name(name);
    // the offending name has to be in the message: it comes straight from the
    // user's `xc` input line, which is free-form (`xc GGA_X_PBE 0.75 ...`)
    if (id <= 0) throw std::runtime_error(
        "unknown libxc functional in the xc input line: >>" + name + "<<");
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
        // CAM/CAMY are not the only markers: LC/LCY flag a long-range-corrected
        // hybrid, whose exact-exchange coefficient would otherwise be handed to
        // the ordinary (unattenuated) Exchange operator, quietly evaluating a
        // different functional than the one that was asked for.
        if ((flags & XC_FLAGS_HYB_CAM) or (flags & XC_FLAGS_HYB_CAMY)
                or (flags & XC_FLAGS_HYB_LC) or (flags & XC_FLAGS_HYB_LCY)) {
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


/// contract two zeta vectors at grid point i: chi_st = zeta_s . zeta_t

/// chi is the reduced density gradient in the logarithmic representation,
/// sigma_st = rho_s rho_t chi_st. It used to be carried as its own multiwavelet
/// function (dot(grad(ln rho_s), grad(ln rho_t)) projected and truncated), but a
/// separately represented product is not pointwise consistent with the zeta
/// components it came from: near the nuclear cusp the projection error in
/// |grad ln rho|^2 is O(1), so chi_aa -- a sum of squares -- came out negative,
/// sigma_aa/sigma_bb were clamped to their positivity floor while sigma_ab kept
/// its (now unbounded) value, and the total sigma_aa + 2 sigma_ab + sigma_bb went
/// negative. Measured against libxc alone: with the diagonals on the floor and
/// sigma_ab swept, exc stays finite and independent of sigma_ab while vsigma of
/// GGA_C_PBE and MGGA_C_TPSS turns NaN exactly when the total crosses zero. The
/// exchange kernels never fail, so this corrupts the potential and not the energy,
/// and it hits plain gga as hard as meta-gga. libxc's own two-sided clamp
/// |sigma_ab| <= (sigma_aa + sigma_bb)/2 keeps the total non-negative while the
/// diagonals are physical, which is why a floored diagonal is a precondition.
///
/// Contracting the zeta components here instead makes the sigma matrix the exact
/// Gram matrix of the gradients at every grid point, so non-negativity of the
/// diagonals, Cauchy-Schwarz, and non-negativity of the total all hold by
/// construction rather than up to projection error. It is also cheaper: three
/// multiwavelet products per spin pair disappear from prep_xc_args.
///
/// The pointers are deliberately not MADNESS_RESTRICT: with no beta electrons the
/// three beta components all point at the same zero-filled scratch tensor.
static inline double chi_of(const double* sx, const double* sy, const double* sz,
                            const double* tx, const double* ty, const double* tz,
                            const long i) {
    return sx[i]*tx[i] + sy[i]*ty[i] + sz[i]*tz[i];
}

/// contract a zeta vector with itself: chi_ss = |zeta_s|^2 >= 0
static inline double chi_of(const double* sx, const double* sy, const double* sz,
                            const long i) {
    return sx[i]*sx[i] + sy[i]*sy[i] + sz[i]*sz[i];
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
                ddensx[i]=dens[i]*zetaa_x[i];
                ddensy[i]=dens[i]*zetaa_y[i];
                ddensz[i]=dens[i]*zetaa_z[i];
                // sigma = |grad rho|^2 contracted from the very gradient handed to
                // libxc, not read from a separately represented chi -- see chi_of()
                sig[i] = std::max(1.e-14,dens[i]*dens[i]*chi_of(zetaa_x,zetaa_y,zetaa_z,i));
            }

            if (needs_tau()) {
                const double * MADNESS_RESTRICT taua = xc_args[enum_taua].ptr();
                // Substituting zeros here would evaluate the functional at the tau
                // floor and return a plausible-looking but wrong energy, so the
                // missing precondition has to be an error, not a default.
                if (taua==NULL) MADNESS_EXCEPTION("meta-gga without a kinetic energy "
                        "density: XCOperator::set_tau() must be called before the "
                        "functional is evaluated",1);
                tau = madness::Tensor<double>(np);
                double * MADNESS_RESTRICT t = tau.ptr();
                for (long i=0; i<np; i++) {
                    double ti = std::max(tautol,2.0*taua[i]);   // full tau is twice alpha tau
                    // tau >= tau_W is exact, so the Fermi hole curvature stays positive
                    // and the iso-orbital indicators stay in range. Bounding tau from
                    // below rather than clamping sigma down (which is what libxc's
                    // XC_FLAGS_ENFORCE_FHC does) leaves the density gradient untouched.
                    //
                    // tau_W = |grad rho|^2/(8 rho) = sigma/(8 rho), but sigma is rho^2
                    // chi, so tau_W = rho chi/8 -- a product. Forming it as sigma/(8 rho)
                    // reintroduces the division by the density that the chi
                    // representation exists to avoid, and picks up sigma's positivity
                    // floor, which then grows like 1/rho instead of vanishing. The
                    // product needs no guard against rho -> 0 either, so the bound
                    // applies everywhere rather than being skipped where rho is munged.
                    ti = std::max(ti,dens[i]*chi_of(zetaa_x,zetaa_y,zetaa_z,i)/8.0);
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
            const bool no_beta_density = (rhob==NULL);
            if ((rhob==NULL)
                    or (zetab_x==NULL) or (zetab_y==NULL) or (zetab_z==NULL)) {
                dummy=madness::Tensor<double>(np);
            }
            if (rhob==NULL) rhob=dummy.ptr();
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

                ddensx[2*i  ]=ra * zetaa_x[i];
                ddensx[2*i+1]=rb * zetab_x[i];
                ddensy[2*i  ]=ra * zetaa_y[i];
                ddensy[2*i+1]=rb * zetab_y[i];
                ddensz[2*i  ]=ra * zetaa_z[i];
                ddensz[2*i+1]=rb * zetab_z[i];

                // Contract the sigma matrix from the same zeta vectors, so that it is
                // the exact Gram matrix of grad(rho_a), grad(rho_b) at this point:
                // both diagonals are non-negative, |sigma_ab| obeys Cauchy-Schwarz,
                // and sigma_aa + 2 sigma_ab + sigma_bb = |grad rho|^2 >= 0. libxc
                // relies on all three -- the total in particular: the correlation
                // kernels return NaN for vsigma as soon as it goes negative.
                double saa = ra * ra * chi_of(zetaa_x,zetaa_y,zetaa_z,i);
                double sbb = rb * rb * chi_of(zetab_x,zetab_y,zetab_z,i);
                // sigma_ab = grad(rho_a).grad(rho_b) is bilinear and may legitimately
                // be negative; only its magnitude is bounded. Do not clamp the sign.
                double sab = ra * rb * chi_of(zetaa_x,zetaa_y,zetaa_z,
                                              zetab_x,zetab_y,zetab_z,i);
                // the positivity floor is for libxc's benefit; raising a diagonal
                // could break Cauchy-Schwarz on its own, so re-impose the bound after
                saa = std::max(1.e-14,saa);
                sbb = std::max(1.e-14,sbb);
                const double cs = std::sqrt(saa*sbb);
                sab = std::min(cs,std::max(-cs,sab));

                sig[3*i  ]  = saa;   // aa
                sig[3*i+1]  = sab;   // ab
                sig[3*i+2]  = sbb;   // bb
            }

            if (needs_tau()) {
                const double * MADNESS_RESTRICT taua = xc_args[enum_taua].ptr();
                const double * MADNESS_RESTRICT taub = xc_args[enum_taub].ptr();
                if (taua==NULL) MADNESS_EXCEPTION("meta-gga without a kinetic energy "
                        "density: XCOperator::set_tau() must be called before the "
                        "functional is evaluated",1);
                // a missing beta tau is legitimate only when there are no beta
                // electrons -- the same condition under which the beta density is
                // absent above, and then every beta contribution vanishes anyway
                madness::Tensor<double> taudummy;
                if (taub==NULL) {
                    MADNESS_CHECK_THROW(no_beta_density, "meta-gga: beta kinetic "
                            "energy density missing although beta electrons are "
                            "present -- pass bmo to XCOperator::set_tau()");
                    taudummy=madness::Tensor<double>(np);
                    taub=taudummy.ptr();
                }

                tau = madness::Tensor<double>(np*2L);
                double * MADNESS_RESTRICT t = tau.ptr();
                for (long i=0; i<np; i++) {
                    double ta = std::max(tautol,taua[i]);
                    double tb = std::max(tautol,taub[i]);
                    // tau_W,s = rho_s chi_ss/8 in each spin channel, as a product rather
                    // than sigma_ss/(8 rho_s) -- see the spin-restricted branch above
                    ta = std::max(ta,dens[2*i  ]*chi_of(zetaa_x,zetaa_y,zetaa_z,i)/8.0);
                    tb = std::max(tb,dens[2*i+1]*chi_of(zetab_x,zetab_y,zetab_z,i)/8.0);
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
        case XC_FAMILY_HYB_LDA:
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
        case XC_FAMILY_HYB_LDA:
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
    MADNESS_CHECK_THROW(not spin_polarized and ispin==0,
                        "XCfunctional::fxc_apply: only spin-restricted, ispin==0 is implemented");

    // copy quantities from t to rho and sigma
    Tensor<double> rho,sigma, rho_pt, sigma_pt;   // rho=2rho_alpha, sigma=4sigma_alpha
    std::vector<Tensor<double> > drho(3), drho_pt(3);
    // The response path prepares rho_pt and sigma_pt only -- there is no perturbed
    // kinetic energy density, and tau's variation is orbital-dependent rather than
    // a functional of the perturbed density. Reject meta-ggas here instead of
    // letting them reach the unknown-family default below, so the message names
    // the actual limitation.
    if (needs_tau()) MADNESS_EXCEPTION("meta-gga response is not implemented: the "
            "xc kernel has no perturbed kinetic energy density",1);
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
        case XC_FAMILY_HYB_LDA:
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
