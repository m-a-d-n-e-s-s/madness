/*
 * CCPotentials.cc
 *
 *  Created on: 4 Jan 2017
 *      Author: kottmanj
 */

#include "CCPotentials.h"


namespace madness {

// some functors used by tests
static double functor_x(const coord_3d& r) { return r[0]; }

static double functor_y(const coord_3d& r) { return r[1]; }

static double functor_r2(const coord_3d& r) { return (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]); }

CCPotentials::CCPotentials(World& world_,  std::shared_ptr<Nemo> nemo, const CCParameters& param)
        : world(world_),
          nemo_(nemo),
          parameters(param),
          //mo_ket_(make_mo_ket(nemo)),
          //mo_bra_(make_mo_bra(nemo)),
          //orbital_energies_(init_orbital_energies(nemo))
//          g12(std::shared_ptr<CCConvolutionOperator(world, OT_G12, param)), f12(world, OT_F12, param),
          corrfac(world, param.gamma(), 1.e-7, nemo->get_calc()->molecule),
          get_potentials(world, param),
          output(world) {
    g12=std::shared_ptr<CCConvolutionOperator<double,3>>(new CCConvolutionOperator<double,3>(world,OpType::OT_G12,param));
    f12=std::shared_ptr<CCConvolutionOperator<double,3>>(new CCConvolutionOperator<double,3>(world,OpType::OT_F12,param));
    output.debug = parameters.debug();
//    reset_nemo(nemo);
//    g12.update_elements(mo_bra_, mo_ket_);
//    g12.sanity();
//    f12.update_elements(mo_bra_, mo_ket_);
//    f12.sanity();
}

madness::CC_vecfunction
CCPotentials::make_mo_bra(const Nemo& nemo) const {
    vector_real_function_3d tmp = mul(world, nemo.ncf->square(), nemo.get_calc()->amo);
    set_thresh(world, tmp, parameters.thresh_3D());
    truncate(world, tmp);
    reconstruct(world, tmp);
    CC_vecfunction mo_bra(tmp, HOLE);
    return mo_bra;
}

madness::CC_vecfunction
CCPotentials::make_mo_ket(const Nemo& nemo) const {
    vector_real_function_3d tmp = nemo.get_calc()->amo;
    set_thresh(world, tmp, parameters.thresh_3D());
    truncate(world, tmp);
    reconstruct(world, tmp);
    CC_vecfunction mo_ket(tmp, HOLE);
    return mo_ket;
}

std::vector<double>
CCPotentials::init_orbital_energies(const Nemo& nemo) const {
    std::vector<double> eps;
    if (world.rank() == 0) std::cout << "SCF Orbital Energies are:\n";

    for (size_t i = 0; i < mo_ket_.size(); i++) {
        eps.push_back(nemo.get_calc()->aeps(i));
        if (world.rank() == 0) std::cout << nemo.get_calc()->aeps(i);
    }
    if (world.rank() == 0) std::cout << "\n" << std::endl;

    return eps;
}

madness::CCPair
CCPotentials::make_pair_gs(const real_function_6d& u, const CC_vecfunction& tau, const size_t i, const size_t j) const {
    CCTimer time(world, "make pair u" + std::to_string(int(i)) + std::to_string(int(j)));
    MADNESS_ASSERT(tau.type == PARTICLE || tau.type == HOLE);
    // for  MP2: tau is empty or Hole states, the function will give back mo_ket_
    // for freeze!=0 the function will give back (mo0,mo1,...,t_freeze,t_freeze+1,...)
    const CC_vecfunction t = make_t_intermediate(tau);
    // functions for the projector
    CC_vecfunction pt;
    if (!parameters.QtAnsatz()) pt = mo_ket_;
    else {
        pt = make_full_t_intermediate(tau);
    }
    std::vector<CCPairFunction<double,6>> functions;
    CCPairFunction<double,6> u_part(u);
    functions.push_back(u_part);
    if (parameters.decompose_Q()) {
        CCPairFunction<double,6> f_part(f12, t(i), t(j));
        functions.push_back(f_part);
        CCPairFunction<double,6> Ot1 = apply_Ot(f_part, pt, 1);
        CCPairFunction<double,6> Ot2 = apply_Ot(f_part, pt, 2);
        CCPairFunction<double,6> PQ = apply_Qt(Ot1, pt, 2, 0.5);
        CCPairFunction<double,6> QP = apply_Qt(Ot2, pt, 1, 0.5);
        // assign signs
        PQ.invert_sign();
        QP.invert_sign();
        functions.push_back(PQ);
        functions.push_back(QP);
    } else {
        // TODO: turn this into separated form, needed (only?) in the energy computation
        real_function_6d ftt = make_f_xy(t(i), t(j));
        real_function_6d Qftt = apply_Q12t(ftt, pt);
        Qftt.truncate();
        CCPairFunction<double,6> residual(Qftt);
        functions.push_back(residual);
    }
    CalcType ctype = CT_CC2;
    if (t.type == HOLE) ctype = CT_MP2;

    CCPair pair(i, j, GROUND_STATE, ctype, functions);
    if (ctype == CT_CC2 && parameters.QtAnsatz()) MADNESS_ASSERT(pt.type == MIXED);

    if (ctype == CT_CC2 && !parameters.QtAnsatz()) MADNESS_ASSERT(pt.type == HOLE);

    if (parameters.decompose_Q()) MADNESS_ASSERT(functions.size() == 4);
//    else
//        MADNESS_ASSERT(functions.size() == 2);

    const double bsh_eps = get_epsilon(i, j);
    pair.bsh_eps = bsh_eps;
    time.info();
    return pair;
}

///// compute the matrix element <ij | g12 Q12 f12 | phi^0>
///// @return 	the energy <ij | g Q f | kl>
//double CCPotentials::compute_gQf(const int i, const int j, ElectronPair& pair) const {
//
//    // for clarity of notation
//    const int k = pair.i;
//    const int l = pair.j;
//
//    // the ket space
//    const real_function_3d& ket_i = hf->nemo(i);
//    const real_function_3d& ket_j = hf->nemo(j);
//
//    // the bra space
//    const real_function_3d& bra_k = hf->R2orbital(k);
//    const real_function_3d& bra_l = hf->R2orbital(l);
//
//    // compute <ij| fg |kl>: do it in 3D as (ik| fg |jl)
//    // the operator fg can be rewritten as 1/r12 - f/r12
//    // i.e. as poisson kernel and a bsh kernel. Note the
//    // the bsh kernel includes a factor of 1/(4 pi)
//    const real_function_3d ik = ket_i * bra_k;
//    const real_function_3d jl = ket_j * bra_l;
//
//    // make all the operators that we need
//    const double fourpi = 4.0 * constants::pi;
//    real_convolution_3d fg = BSHOperator<3>(world, corrfac.gamma(), lo,
//                                            bsh_eps / fourpi);
//    real_convolution_3d gg = CoulombOperator(world, lo, bsh_eps);
//    real_convolution_3d slaterf12 = SlaterF12Operator(world, corrfac.gamma(),
//                                                      lo, bsh_eps / fourpi);
//
//    //  < ij | fg | kl >
//    const real_function_3d ik_fg = (gg)(ik) - fourpi * fg(ik);
//    const double a = inner(ik_fg, jl) / (2.0 * corrfac.gamma());
//    if (world.rank() == 0)
//        printf("<%d%d | f/r              | %d%d>  %12.8f\n", i, j, k, l, a);
//
//    // compute <ij| f (O1 + O2) g | ij>
//
//    // compute bra space xi(ik,j)^dagger, i.e. the hermitian conjugate of xi
//    // the index k is implicit in the vector of functions
//    // naming: xi _ orbitals _ operator _ hc
//    std::vector<real_function_3d> xi_ij_g_ket = make_xi(ket_i, ket_j, *poisson,
//                                                        false);       // xi_{i,m*},j
//    std::vector<real_function_3d> xi_ji_g_ket = make_xi(ket_j, ket_i, *poisson,
//                                                        false);       // xi_{j,m*},i
//
//    std::vector<real_function_3d> xi_ij_f_bra = make_xi(bra_k, bra_l, slaterf12,
//                                                        true);       // xi_{i*,m},j*
//    std::vector<real_function_3d> xi_ji_f_bra = make_xi(bra_l, bra_k, slaterf12,
//                                                        true);       // xi_{j*,m},i*
//
//    // in the following do NOT use antisymmetrized pair functions:
//    // |ij> -> 0.5 |ij - ji>
//
//    // < ij | f12 O1 g12 | kl >
//    //   = \sum_m <i(1) j(2) | f12 | m(1) >< m(3) | g23 | k(3) l(2)>
//    //   = \sum_m < chi^f_i*,m(2) j*(2) | chi^g_k,m*(2) l(2) >
//    //   = \sum_m < xi^f_im,j | xi^g_km,l >
//    const double o1a = inner(world, xi_ij_f_bra, xi_ij_g_ket).sum();
//    if (world.rank() == 0)
//        printf("<%d%d | f12 O1 g12       | %d%d>  %12.8f\n", i, j, k, l, o1a);
//
//    // < ij | f12 O2 g12 | kl >
//    //    = \sum_m <i(1) j(2) | f12 | m(2) >< m(3) | g13 | k(1) l(3)>
//    //    = \sum_m <chi^f_j*,m(1) i*(1) | chi^g_l,m*(1) k(1) >
//    //    = \sum_m < xi^f_jm,i | xi^g_lm,k >
//    const double o2a = inner(world, xi_ji_f_bra, xi_ji_g_ket).sum();
//    if (world.rank() == 0)
//        printf("<%d%d | f12 O2 g12       | %d%d>  %12.8f\n", i, j, k, l, o2a);
//
//    // compute <ij| f O1 O2 g | kl>  // why do I need to swap ij in g_ijkl??
//    const Tensor<double> f_ijmn = matrix_inner(world, xi_ij_f_bra, hf->nemos());
//    const Tensor<double> g_ijmn = matrix_inner(world, hf->R2orbitals(),
//                                               xi_ji_g_ket);
//    const double o12 = f_ijmn.trace(g_ijmn);
//    if (world.rank() == 0)
//        printf("<%d%d | f12 O12 g12      | %d%d>  %12.8f\n", i, j, k, l, o12);
//
//    const double e = a - o1a - o2a + o12;
//    if (world.rank() == 0)
//        printf("<%d%d | g Q12 f          | %d%d>  %12.8f\n", i, j, k, l, e);
//
//    return e;
//}

madness::CCPair
CCPotentials::make_pair_ex(const real_function_6d& u, const CC_vecfunction& tau, const CC_vecfunction& x,
                           const size_t i, const size_t j, const CalcType ctype) const {
    MADNESS_ASSERT(tau.type == PARTICLE || tau.type == HOLE);
    MADNESS_ASSERT(x.type == RESPONSE);
    MADNESS_ASSERT(ctype == CT_CISPD || ctype == CT_LRCC2 || ctype == CT_ADC2);
    MADNESS_ASSERT(!(i < parameters.freeze()));
    MADNESS_ASSERT(!(j < parameters.freeze()));
    // for  CIS(D): tau is empty or Hole states, the function will give back mo_ket_
    // for freeze!=0 the function will give back (mo0,mo1,...,t_freeze,t_freeze+1,...)
    const CC_vecfunction t = make_t_intermediate(tau).copy();
    // functions for the projector
    CC_vecfunction pt;
    if (!parameters.QtAnsatz()) pt = mo_ket_.copy();
    else {
        pt = make_full_t_intermediate(tau).copy();
    }
    MADNESS_ASSERT(pt.size() == mo_ket_.size());
    std::vector<CCPairFunction<double,6>> functions;
    CCPairFunction<double,6> u_part(u);
    functions.push_back(u_part);
    if (parameters.decompose_Q()) {
        CCPairFunction<double,6> f_xt(f12, x(i), t(j));
        functions.push_back(f_xt);
        CCPairFunction<double,6> f_tx(f12, t(i), x(j));
        functions.push_back(f_tx);
        {
            CCPairFunction<double,6> Ot1_xt = apply_Ot(f_xt, pt, 1);     // O1t(f|xt>)
            CCPairFunction<double,6> OtQt_xt = apply_Qt(Ot1_xt, pt, 2, 0.5);     // O1t(1-0.5*O2t)f|xt>
            functions.push_back(OtQt_xt.invert_sign());     // - "
        }
        {
            CCPairFunction<double,6> Ot2_xt = apply_Ot(f_xt, pt, 2);     // O2t(f|xt>)
            CCPairFunction<double,6> QtOt_xt = apply_Qt(Ot2_xt, pt, 1, 0.5);     // (1-0.5*O1t)O2t(f|xt>)
            functions.push_back(QtOt_xt.invert_sign());     // - "
        }
        {
            CCPairFunction<double,6> Ot1_tx = apply_Ot(f_tx, pt, 1);     // O1t(f|tx>)
            CCPairFunction<double,6> OtQt_tx = apply_Qt(Ot1_tx, pt, 2, 0.5);     // O1t(1-0.5*O2t)f|tx>
            functions.push_back(OtQt_tx.invert_sign());     // - "
        }
        {
            CCPairFunction<double,6> Ot2_tx = apply_Ot(f_tx, pt, 2);     // O2t(f|tx>)
            CCPairFunction<double,6> QtOt_tx = apply_Qt(Ot2_tx, pt, 1, 0.5);     // (1-0.5*O1t)O2t(f|tx>)
            functions.push_back(QtOt_tx.invert_sign());     // - "
        }
        if (parameters.QtAnsatz()) {
            CCPairFunction<double,6> ftt(f12, t(i), t(j));     // f|tt>
            CCPairFunction<double,6> O1x_tt = apply_Ot(ftt, x, 1);     // O1x(f|tt>)
            CCPairFunction<double,6> OxQt_tt = apply_Qt(O1x_tt, pt, 2);     // O1xQt(f|tt>)
            functions.push_back(OxQt_tt.invert_sign());     // - "
            CCPairFunction<double,6> O2x_tt = apply_Ot(ftt, x, 2);     // O2x(f|tt>)
            CCPairFunction<double,6> QtOx_tt = apply_Qt(O2x_tt, pt, 1);     // Q1tO2x(f|tt>)
            functions.push_back(QtOx_tt.invert_sign());     // - "
        }
    } else {
        // use 6D residuals: Qt*f*|xitj+tixj> - (QtOx+OxQt)|titj> or without QtAnsatz: Qf|xitj+tixj>
        if (parameters.QtAnsatz()) {
            real_function_6d ftx = make_f_xy(t(i), x(j));
            real_function_6d fxt = make_f_xy(x(i), t(j));
            real_function_6d tmp = ftx + fxt;
            real_function_6d function_part = apply_Q12t(tmp, pt);
            real_function_6d ftt = make_f_xy(t(i), t(j));
            Projector<double, 3> Ox(x.get_vecfunction());
            Projector<double, 3> Ot(pt.get_vecfunction());
            real_function_6d O1x = Ox(ftt, 1);
            real_function_6d O2x = Ox(ftt, 2);
            real_function_6d OtOx = Ot(O2x, 1);
            real_function_6d OxOt = Ot(O1x, 2);
            real_function_6d projector_part = (O1x - OxOt + O2x - OtOx);
            real_function_6d residual = function_part - projector_part;
            CCPairFunction<double,6> res(residual);
            functions.push_back(res);
        } else {
            real_function_6d ftx = make_f_xy(t(i), x(j));
            real_function_6d fxt = make_f_xy(x(i), t(j));
            real_function_6d tmp = ftx + fxt;
            real_function_6d residual = apply_Q12t(tmp, mo_ket_);
            residual.truncate();
            CCPairFunction<double,6> res(residual);
            functions.push_back(res);
        }
    }
    CCPair pair(i, j, EXCITED_STATE, ctype, functions);
    if (parameters.decompose_Q()) {
        if (parameters.QtAnsatz()) MADNESS_ASSERT(functions.size() == 9);
        else
            MADNESS_ASSERT(functions.size() == 7);
    } else
        MADNESS_ASSERT(functions.size() == 2);

    MADNESS_ASSERT(x.omega != 0.0);
    const double bsh_eps = get_epsilon(i, j) + x.omega;
    pair.bsh_eps = bsh_eps;
    return pair;
}

double
CCPotentials::compute_pair_correlation_energy(const CCPair& u, const CC_vecfunction& singles) const {
    CCTimer timer(world, "Compute Correlation Energy");
    MADNESS_ASSERT(u.type == GROUND_STATE);
    if (singles.functions.empty()) MADNESS_ASSERT(u.ctype == CT_MP2);

    const bool print_details=(world.rank()==0 and parameters.debug());
    if (parameters.debug()) output("Compute pair-correlation energy of pair " + u.name());
    double result = 0.0;
    const CCFunction<double,3>& mobi = mo_bra_(u.i);
    const CCFunction<double,3>& mobj = mo_bra_(u.j);
    const bool symmetric = (u.i == u.j);


    for (size_t mm = 0; mm < u.functions.size(); mm++) {
        double tmp = 0.0;
        const double part1 = make_xy_op_u(mobi, mobj, *g12, u.functions[mm]);
        if (symmetric) tmp = part1;
        else     //if(world.rank()==0) std::cout << std::fixed << std::setprecision(10) << part1 << "\n";
        {
            const double part2 = make_xy_op_u(mobj, mobi, *g12, u.functions[mm]);
            tmp = 2.0 * (2.0 * part1 - part2);     // non symmetric pairs -> offdiagonal -> count twice
        }
        result += tmp;
        if (print_details)
            std::cout << std::setfill(' ') << std::setw(15) << "from " + u.functions[mm].name() + "="
                      << std::setfill(' ') << std::fixed << std::setprecision(10) << tmp << "\n";
    }
    if (u.ctype == CT_CC2 && !singles.get_vecfunction().empty()) {
        MADNESS_ASSERT(singles.type == PARTICLE);
        const double omega_s = 2.0 * mobi.inner((*g12)(mobj, singles(u.j)) * singles(u.i).function) -
                               mobi.inner((*g12)(mobj, singles(u.i)) * singles(u.j).function);
        if (print_details)
            std::cout << std::setw(15) << "from singles=" << std::setfill(' ') << std::fixed << std::setprecision(10)
                      << omega_s << "\n\n";

        result += omega_s;
    }
    // if (world.rank() == 0) std::cout << "------------\n" << std::fixed << std::setprecision(10) << result << "\n\n";

    timer.info(parameters.debug());
    return result;
}

double
CCPotentials::compute_cc2_correlation_energy(const CC_vecfunction& singles, const Pairs<CCPair>& doubles) const {
    MADNESS_ASSERT(singles.type == PARTICLE);
    CCTimer time(world, "Computing CC2 Correlation Energy");
    output.section("Computing CC2 Correlation Energy");
    double result = 0.0;
    for (const auto& tmp : doubles.allpairs) {
        const size_t i = tmp.second.i;
        const size_t j = tmp.second.j;
        const double omega = compute_pair_correlation_energy(tmp.second, singles);
        result += omega;
        if (world.rank() == 0)
            std::cout << std::fixed << "omega  " << i << j << " =" << std::setprecision(10) << omega << "\n";
    }
    if (world.rank() == 0) std::cout << std::fixed << "sum      " << " =" << std::setprecision(10) << result << "\n";

    time.info();
    return result;
}

double
CCPotentials::compute_kinetic_energy(const vector_real_function_3d& xbra, const vector_real_function_3d& xket) const {
    Kinetic<double, 3> T(world);
    double kinetic = 0.0;
    for (size_t k = 0; k < xket.size(); k++)
        kinetic += T(xbra[k], xket[k]);
    return kinetic;
}

double
CCPotentials::compute_cis_expectation_value(const CC_vecfunction& x, const vector_real_function_3d& V,
                                            const bool print) const {
    const vector_real_function_3d xbra = make_bra(x);
    const vector_real_function_3d xket = x.get_vecfunction();
    const double kinetic = compute_kinetic_energy(xbra, xket);
    const double norm = sqrt(inner(world, xbra, xket).sum());
    double eps = 0.0;
    for (size_t k = 0; k < xket.size(); k++)
        eps -= get_orbital_energies()[k + parameters.freeze()] * xbra[k].inner(xket[k]);
    double potential = inner(world, xbra, V).sum();
    const double result = 1.0 / (norm * norm) * (potential + kinetic + eps);
    if (world.rank() == 0 && print) {
        std::cout << "CCS Expectation Value:\n--------\n";
        std::cout << "Kinetic-Energy  =" << std::fixed << std::setprecision(8) << kinetic << "\n";
        std::cout << "Potential-Energy=" << std::fixed << std::setprecision(8) << potential << "\n";
        std::cout << "ei*<xi|xi>      =" << std::fixed << std::setprecision(8) << eps << "\n";
        std::cout << "||x||           =" << std::fixed << std::setprecision(8) << norm << "\n";
        std::cout << "Expectationvalue=" << std::fixed << std::setprecision(8) << result << "\n--------\n";
    }
    return result;
}

double
CCPotentials::compute_excited_pair_energy(const CCPair& d, const CC_vecfunction& x) const {
    const CC_vecfunction xbra(make_bra(x), RESPONSE, parameters.freeze());
    const CCFunction<double,3>& xbi = xbra(d.i);
    const CCFunction<double,3>& mobj = mo_bra_(d.j);
    double result = 0.0;
    double s2b = 2.0 * make_xy_op_u(xbi, mobj, *g12, d.functions) - make_xy_op_u(mobj, xbi, *g12, d.functions);
    double s2c = 0.0;
    for (const auto& ktmp : x.functions) {
        const size_t k = ktmp.first;
        const real_function_3d j_igk = (*g12)(mo_bra_(d.i), mo_ket_(k)) * mo_bra_(d.j).function;
        s2c -= 2.0 * make_xy_u(xbra(k), j_igk, d.functions) - make_xy_u(j_igk, xbra(k), d.functions);
    }
    result = s2b + s2c;
    if (world.rank() == 0) {
        std::cout << std::fixed << std::setprecision(10) << "\nExcited Pair Energy: " << "S2b=" << s2b << ", S2c="
                  << s2c << ", Both=" << result << "\n\n";
    }
    return result;
}

double
CCPotentials::compute_cispd_energy(const CC_vecfunction& x, const Pairs<CCPair> mp2, const Pairs<CCPair> cispd) const {
    MADNESS_ASSERT(x.type == RESPONSE);
    output.section("Compute CIS(D) Energy");
    const CC_vecfunction xbra(make_bra(x), RESPONSE, parameters.freeze());
    double s2b = 0.0;
    double s2c = 0.0;
    double s4a = 0.0;
    double s4b = 0.0;
    double s4c = 0.0;
    for (const auto& itmp : x.functions) {
        const size_t i = itmp.first;
        for (const auto& ktmp : x.functions) {
            // s2b part: <xi|s2b_i>
            const size_t k = ktmp.first;
            s2b += 2.0 * make_xy_op_u(xbra(i), mo_bra_(k), *g12, get_pair_function(cispd, i, k)) -
                   make_xy_op_u(mo_bra_(k), xbra(i), *g12, get_pair_function(cispd, i, k));
            const real_function_3d kgi = (*g12)(mo_bra_(k), mo_ket_(i));
            const real_function_3d kgxi = (*g12)(mo_bra_(k), x(i));
            const real_function_3d kgxk = (*g12)(mo_bra_(k), x(k));
            for (const auto& ltmp : x.functions) {
                // s2c part: <xi|s2c_i>
                const size_t l = ltmp.first;
                const real_function_3d k_lgxk = mo_bra_(k).function * (*g12)(mo_bra_(l), x(k));
                const real_function_3d l_kgxk = mo_bra_(l).function * kgxk;
                const real_function_3d l_kgi = mo_bra_(l).function * kgi;
                const real_function_3d l_kgxi = mo_bra_(l).function * kgxi;
                s2c += 2.0 * make_xy_u(xbra(i), l_kgi, get_pair_function(cispd, k, l)) -
                       make_xy_u(l_kgi, xbra(i), get_pair_function(cispd, k, l));
                const double xil = xbra(i).function.inner(x(l).function);
                s4a += xil * (2.0 * make_xy_op_u(mo_bra_(l), mo_bra_(k), *g12, get_pair_function(mp2, i, k)) -
                              make_xy_op_u(mo_bra_(k), mo_bra_(l), *g12, get_pair_function(mp2, i, k)));
                s4b += 2.0 * make_xy_u(xbra(i), l_kgxi, get_pair_function(mp2, k, l)) -
                       make_xy_u(l_kgxi, xbra(i), get_pair_function(mp2, k, l));
                s4c += 4.0 * make_xy_u(xbra(i), l_kgxk, get_pair_function(mp2, i, l)) -
                       2.0 * make_xy_u(l_kgxk, xbra(i), get_pair_function(mp2, i, l)) -
                       2.0 * make_xy_u(xbra(i), k_lgxk, get_pair_function(mp2, i, l))
                       + make_xy_u(k_lgxk, xbra(i), get_pair_function(mp2, i, l));
            }
        }
    }
    const double result = s2b - s2c - s4a - s4b + s4c;
    if (world.rank() == 0) {
        std::cout << std::fixed << std::setprecision(10) << "CIS(D) Correction\n" << "s2b =" << s2b << "\n" << "s2c ="
                  << -s2c << "\n" << "s4a =" << -s4a << "\n" << "s4b =" << -s4b << "\n" << "s4c ="
                  << s4c << "\n" << "-----------------\n" << "sum =" << result << "\n";
    }
    return result;
}

double
CCPotentials::compute_cc2_excitation_energy(const CC_vecfunction& stau, const CC_vecfunction& sx,
                                            const Pairs<CCPair> dtau, const Pairs<CCPair> dx) const {
    vector_real_function_3d tmp = mul(world, nemo_->ncf->square(), sx.get_vecfunction());
    truncate(world, tmp);
    CC_vecfunction xbra(tmp, RESPONSE, parameters.freeze());
    const double xbrax = inner(world, xbra.get_vecfunction(), sx.get_vecfunction()).sum();
    double result = potential_energy_ex(xbra, stau, dtau, sx, dx, POT_s3a_);
    result += potential_energy_ex(xbra, stau, dtau, sx, dx, POT_s3b_);
    result += potential_energy_ex(xbra, stau, dtau, sx, dx, POT_s3c_);
    result += potential_energy_ex(xbra, stau, dtau, sx, dx, POT_s5b_);
    result += potential_energy_ex(xbra, stau, dtau, sx, dx, POT_s5c_);
    result += potential_energy_ex(xbra, stau, dtau, sx, dx, POT_s6_);
    result += potential_energy_ex(xbra, stau, dtau, sx, dx, POT_s2b_);
    result += potential_energy_ex(xbra, stau, dtau, sx, dx, POT_s2c_);
    result += potential_energy_ex(xbra, stau, dtau, sx, dx, POT_s4a_);
    result += potential_energy_ex(xbra, stau, dtau, sx, dx, POT_s4b_);
    result += potential_energy_ex(xbra, stau, dtau, sx, dx, POT_s4c_);
    return 1.0 / xbrax * result;
}

madness::real_function_6d
CCPotentials::fock_residue_6d(const CCPair& u) const {
    if (u.function().norm2() == 0.0) {
        output("Pair-Function is zero so this is the first iteration ... skipping Fock residue");
        return real_factory_6d(world);
    }
    // make the special points for the cusps in the nuclear potential
    std::vector<Vector<double, 3> > sp3d = nemo_->get_calc()->molecule.get_all_coords_vec();
    std::vector<Vector<double, 6> > sp6d;
    for (size_t i = 0; i < sp3d.size(); i++) {
        Vector<double, 6> tmp;
        for (size_t j = 0; j < 3; j++) {
            tmp[j] = sp3d[i][j];
            tmp[j + 3] = sp3d[i][j];
        }
        sp6d.push_back(tmp);
    }
    // make the coulomb and local Un part with the composite factory
    real_function_3d hartree_potential = real_factory_3d(world);
    for (const auto& tmp : mo_ket_.functions)
        hartree_potential += (*g12)(mo_bra_(tmp.first), mo_ket_(tmp.first));
    real_function_3d local_part = (2.0 * hartree_potential + nemo_->ncf->U2());
    if (parameters.debug()) local_part.print_size("vlocal");

    if (parameters.debug()) u.function().print_size(u.name());

    // Contruct the BSH operator in order to screen
    double bsh_eps = u.bsh_eps;
    real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2.0 * bsh_eps), parameters.lo(), parameters.thresh_bsh_6D());
    op_mod.modified() = true;
    // Make the CompositeFactory
    real_function_6d vphi = CompositeFactory<double, 6, 3>(world).ket(copy(u.function())).V_for_particle1(
            copy(local_part)).V_for_particle2(copy(local_part)).special_points(sp6d);
    // Screening procedure
    vphi.fill_nuclear_cuspy_tree(op_mod, 0);
    if (parameters.debug()) vphi.print_size("vlocal|u>");

    vphi.truncate().reduce_rank();
    if (parameters.debug()) vphi.print_size("vlocal|u>");

    real_function_6d Un1 = real_factory_6d(world);
    real_function_6d Un2 = real_factory_6d(world);
    for (int axis = 0; axis < 3; ++axis) {
        real_derivative_6d D = free_space_derivative<double, 6>(world, axis);
        const real_function_6d Du = D(u.function()).truncate();
        const real_function_3d U1_axis = nemo_->ncf->U1(axis);
        double tight_thresh = parameters.thresh_6D();
        real_function_6d x = CompositeFactory<double, 6, 3>(world).ket(copy(Du)).V_for_particle1(copy(U1_axis)).thresh(
                tight_thresh).special_points(sp6d);
        x.fill_nuclear_cuspy_tree(op_mod, 1);
        if (x.norm2() < tight_thresh) x.print_size("Un_axis_" + stringify(axis));

        if (x.norm2() < tight_thresh) output.warning("||Un|u>|| is below the threshold");

        Un1 += x;
        Un1.truncate().reduce_rank();
        if (parameters.debug()) Un1.print_size("Un1");
    }
    if (u.i == u.j) {
        output(u.name() + " is a diagonal pair: Exploting permutation symmetry");
        Un2 = swap_particles(Un1);
    } else {
        for (int axis = 3; axis < 6; ++axis) {
            real_derivative_6d D = free_space_derivative<double, 6>(world, axis);
            const real_function_6d Du = D(u.function()).truncate();
            const real_function_3d U1_axis = nemo_->ncf->U1(axis % 3);
            double tight_thresh = parameters.thresh_6D();
            real_function_6d x = CompositeFactory<double, 6, 3>(world).ket(copy(Du)).V_for_particle2(
                    copy(U1_axis)).thresh(tight_thresh).special_points(sp6d);
            x.fill_nuclear_cuspy_tree(op_mod, 2);
            if (x.norm2() < tight_thresh) x.print_size("Un_axis_" + stringify(axis));

            if (x.norm2() < tight_thresh) output.warning("||Un|u>|| is below the threshold");

            Un2 += x;
            Un2.truncate().reduce_rank();
            if (parameters.debug()) Un2.print_size("Un2");
        }
    }
    vphi += (Un1 + Un2);
    vphi.truncate().reduce_rank();
    if (parameters.debug()) vphi.print_size("(Un + J1 + J2)|u>");

    // Exchange Part
    vphi = (vphi - K(u.function(), u.i == u.j)).truncate().reduce_rank();
    if (parameters.debug()) vphi.print_size("Fock-Residue");

    return vphi;
}


madness::real_function_6d
CCPotentials::fock_residue_6d_macrotask(World& world, const CCPair& u, const CCParameters& parameters,
                                        const std::vector< madness::Vector<double,3> >& all_coords_vec,
                                        const std::vector<real_function_3d>& mo_ket,
                                        const std::vector<real_function_3d>& mo_bra,
                                        const std::vector<real_function_3d>& U1,
                                        const real_function_3d& U2) {
    if (u.function().norm2() == 0.0) {
        print("Pair-Function is zero so this is the first iteration ... skipping Fock residue\n");
        return real_factory_6d(world);
    }
    // make the special points for the cusps in the nuclear potential
    std::vector<Vector<double, 3> > sp3d = all_coords_vec;
    std::vector<Vector<double, 6> > sp6d;
    for (size_t i = 0; i < sp3d.size(); i++) {
        Vector<double, 6> tmp;
        for (size_t j = 0; j < 3; j++) {
            tmp[j] = sp3d[i][j];
            tmp[j + 3] = sp3d[i][j];
        }
        sp6d.push_back(tmp);
    }
    // make the coulomb and local Un part with the composite factory
    real_convolution_3d g12 = CoulombOperator(world, parameters.lo(), parameters.thresh_poisson());
    real_function_3d hartree_potential = g12(dot(world, mo_bra, mo_ket));
    //for (int i = 0; i < mo_ket.size(); i++)
    //    hartree_potential += g12(mo_bra[i], mo_ket[i]);
    real_function_3d local_part = (2.0 * hartree_potential + U2);
    if (parameters.debug()) local_part.print_size("vlocal");

    if (parameters.debug()) u.function().print_size(u.name());

    // Contruct the BSH operator in order to screen
    double bsh_eps = u.bsh_eps;
    real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2.0 * bsh_eps), parameters.lo(), parameters.thresh_bsh_6D());
    op_mod.modified() = true;
    // Make the CompositeFactory
    real_function_6d vphi = CompositeFactory<double, 6, 3>(world).ket(copy(u.function())).V_for_particle1(
            copy(local_part)).V_for_particle2(copy(local_part)).special_points(sp6d);
    // Screening procedure
    vphi.fill_nuclear_cuspy_tree(op_mod, 0);
    if (parameters.debug()) vphi.print_size("vlocal|u>");

    vphi.truncate().reduce_rank();
    if (parameters.debug()) vphi.print_size("vlocal|u>");

    real_function_6d Un1 = real_factory_6d(world);
    real_function_6d Un2 = real_factory_6d(world);
    for (int axis = 0; axis < 3; ++axis) {
        real_derivative_6d D = free_space_derivative<double, 6>(world, axis);
        const real_function_6d Du = D(u.function()).truncate();
        const real_function_3d U1_axis = U1[axis];
        double tight_thresh = parameters.thresh_6D();
        real_function_6d x = CompositeFactory<double, 6, 3>(world).ket(copy(Du)).V_for_particle1(copy(U1_axis)).thresh(
                tight_thresh).special_points(sp6d);
        x.fill_nuclear_cuspy_tree(op_mod, 1);

        if (parameters.debug()) x.print_size("Un_axis_" + stringify(axis));
        Un1 += x;
    }
    if (parameters.debug()) Un1.print_size("Un1");

    if (u.i == u.j) {
        print(u.name() + " is a diagonal pair: Exploting permutation symmetry\n");
        Un2 = madness::swap_particles<double>(Un1);
    } else {
        for (int axis = 3; axis < 6; ++axis) {
            real_derivative_6d D = free_space_derivative<double, 6>(world, axis);
            const real_function_6d Du = D(u.function()).truncate();
            const real_function_3d U1_axis = U1[axis % 3];
            double tight_thresh = parameters.thresh_6D();
            real_function_6d x = CompositeFactory<double, 6, 3>(world).ket(copy(Du)).V_for_particle2(
                    copy(U1_axis)).thresh(tight_thresh).special_points(sp6d);
             x.fill_nuclear_cuspy_tree(op_mod, 2);
            if (parameters.debug()) x.print_size("Un_axis_" + stringify(axis));
            Un2 += x;
        }
        if (parameters.debug()) Un2.print_size("Un2");
    }
    vphi += (Un1 + Un2);
    if (parameters.debug()) vphi.print_size("before truncation (Un + J1 + J2)|u>");
    vphi.truncate().reduce_rank();
    if (parameters.debug()) vphi.print_size("(Un + J1 + J2)|u>");

    // Exchange Part
    vphi = (vphi - K_macrotask(world, mo_ket, mo_bra, u.function(), u.i == u.j, parameters)).truncate().reduce_rank();
    if (parameters.debug()) vphi.print_size("Fock-Residue");

    return vphi;
}

madness::real_function_6d
CCPotentials::make_constant_part_mp2_macrotask(World& world, const CCPair& pair,
                                               const std::vector<real_function_3d>& mo_ket,
                                               const std::vector<real_function_3d>& mo_bra,
                                               const CCParameters& parameters, const real_function_3d& Rsquare,
                                               const std::vector<real_function_3d>& U1,
                                               const std::vector<std::string> argument) {
    MADNESS_ASSERT(pair.ctype == CT_MP2);
    MADNESS_ASSERT(pair.type == GROUND_STATE);
    const std::string i_name = "phi" + stringify(pair.i);
    const std::string j_name = "phi" + stringify(pair.j);
    const double epsilon = pair.bsh_eps;
    const FuncType i_type = HOLE;
    const FuncType j_type = HOLE;

    // load constant part if available
    //if (parameters.no_compute_mp2_constantpart()) {
    //    pair.constant_part=real_factory_6d(world);
    //    load(pair.constant_part,pair.name()+"_const");
    //}
    //if (pair.constant_part.is_initialized()) return false;

    // make screening Operator
    real_convolution_6d Gscreen = BSHOperator<6>(world, sqrt(-2.0 * epsilon),
                                                 parameters.lo(), parameters.thresh_bsh_6D());
    Gscreen.modified() = true;

    print("\nCalculating Constant Part of MP2 pair " + i_name + j_name);
    CCTimer time(world, "Calculating Constant Part of MP2");
    MADNESS_ASSERT(i_type == HOLE);
    MADNESS_ASSERT(j_type == HOLE);
    real_function_6d V = apply_Vreg_macrotask(world, mo_ket, mo_bra, parameters, Rsquare,
                                              U1, pair.i, pair.j, i_type, j_type, argument, &Gscreen);
    if (parameters.debug()) V.print_size("Vreg");

    //MADNESS_ASSERT(t.type == HOLE || t.type == MIXED);
    MADNESS_ASSERT(mo_ket.size() == mo_bra.size());
    StrongOrthogonalityProjector<double, 3> Q(world);
    Q.set_spaces(mo_bra, mo_ket, mo_bra, mo_ket);

//    V = Q(V);
//
//    V.print_size("QVreg");
    real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0 * epsilon), parameters.lo(),
                                           parameters.thresh_bsh_6D());
    G.destructive() = true;
//    real_function_6d GV = -2.0 * G(V);

    // save memory:
    // split application of the BSH operator into high-rank, local part U|ij>, and
    // low-rank, delocalized part (-O1 -O2 +O1O2) U|ij> by splitting the SO operator

    // delocalized part
    auto [left,right]=Q.get_vectors_for_outer_product(V);
    real_function_6d GV1=-2.0*G(left,right);
    GV1.truncate();

    // local part
    real_function_6d GV = -2.0 * G(V);      // note V is destroyed here
    GV.truncate();

    GV+=GV1;
    GV.truncate();
    if (parameters.debug()) GV.print_size("GVreg");
    //MADNESS_ASSERT(t.type == HOLE || t.type == MIXED);
    MADNESS_ASSERT(mo_ket.size() == mo_bra.size());
    GV = Q(GV);

    GV.print_size("GVreg");
    time.info();

    return GV;
}

real_function_6d
CCPotentials::update_pair_mp2_macrotask(World& world, const CCPair& pair, const CCParameters& parameters,
                                             const std::vector< madness::Vector<double,3> >& all_coords_vec,
                                             const std::vector<real_function_3d>& mo_ket,
                                             const std::vector<real_function_3d>& mo_bra,
                                             const std::vector<real_function_3d>& U1,
                                             const real_function_3d& U2, const real_function_6d& mp2_coupling) {

    if (world.rank()==0) print(assign_name(pair.ctype) + "-Microiteration\n");
    CCTimer timer_mp2(world, "MP2-Microiteration of pair " + pair.name());

    if (parameters.debug()) mp2_coupling.print_size("coupling in macrotask");

    double bsh_eps = pair.bsh_eps;
    real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0 * bsh_eps), parameters.lo(), parameters.thresh_bsh_6D());
    G.destructive() = true;

    CCTimer timer_mp2_potential(world, "MP2-Potential of pair " + pair.name());
    real_function_6d mp2_potential = -2.0 * CCPotentials::fock_residue_6d_macrotask(world, pair, parameters,
                                                                                all_coords_vec, mo_ket, mo_bra, U1, U2);
    // add coupling, note sign and factor
    //real_function_6d coupling = mp2_coupling(pair.i, pair.j);
    mp2_potential += 2.0 * mp2_coupling;

    if (parameters.debug()) mp2_potential.print_size(assign_name(pair.ctype) + " Potential");
    mp2_potential.truncate().reduce_rank();
    timer_mp2_potential.info(true, mp2_potential.norm2());

    CCTimer timer_G(world, "Apply Greens Operator on MP2-Potential of pair " + pair.name());
    const real_function_6d GVmp2 = G(mp2_potential);
    timer_G.info(true, GVmp2.norm2());

    //CCTimer timer_addup(world, "Add constant parts and update pair " + pair.name());
    real_function_6d unew = GVmp2 + pair.constant_part;
    if (parameters.debug()) unew.print_size("unew");

    StrongOrthogonalityProjector<double, 3> Q(world);
    Q.set_spaces(mo_bra, mo_ket, mo_bra, mo_ket);
    unew = Q(unew);

    if (parameters.debug())unew.print_size("truncated-unew");
    timer_mp2.info();

    real_function_6d residue = (pair.function() - unew);
    // if (parameters.debug()) residue.print_size("bsh residual");
    residue.truncate(FunctionDefaults<6>::get_thresh()*0.1);
    if (parameters.debug()) residue.print_size("bsh residual, truncated");

    return residue;
}

madness::real_function_6d
CCPotentials::make_constant_part_cc2_gs(const CCPair& u, const CC_vecfunction& tau,
                                        const real_convolution_6d *Gscreen) const {
    output.section("Calculating CC2 Constant Part of Pair " + u.name() + ": Q-Ansatz");
    CCTimer time(world, "Constant Term");
    MADNESS_ASSERT(!parameters.QtAnsatz());
    MADNESS_ASSERT(u.ctype == CT_CC2);
    MADNESS_ASSERT(u.type == GROUND_STATE);
    MADNESS_ASSERT(tau.type == PARTICLE);
    MADNESS_ASSERT(tau.size() != 0);
    // convenience
    const size_t i = u.i;
    const size_t j = u.j;
    const bool symmetric = (i == j);
    // make ti intermediates
    const CCFunction<double,3> ti = make_t_intermediate(tau(i));
    const CCFunction<double,3> tj = make_t_intermediate(tau(j));
    real_function_6d GV;
    real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0 * get_epsilon(ti.i, tj.i)), parameters.lo(),
                                           parameters.thresh_bsh_6D());
    G.destructive() = true;
    CCTimer time_GV(world, "G(Regularization Potential)");
    {
        real_function_6d V = apply_Vreg(ti, tj, Gscreen);
        print_size(V, "Vreg", parameters.debug());
        V = apply_Q12t(V, mo_ket_);
        print_size(V, "QVreg");
        GV = -2.0 * apply_G(V, G);
    }
    //print_size(GV,"GVreg",parameters.debug);
    //GV = apply_Qt(GV,mo_ket_);
    //print_size(GV,"QGVreg");
    time_GV.stop();
    // make Coulomb coupling Potential
    // which is
    // (-OtauQ-QOtau+OtauOtau)g|titj>
    //output.section("\nCalculating Coulomb Coupling Potential of CC2\n");
    CCTimer time_Vcc(world, "G(Coulomb Coupling Potential)");
    real_function_6d GVcc = real_factory_6d(world);
    // make the g12|titj> function as op_decomposed function (not constructed in 6D)
    CCPairFunction<double,6> gtt(g12, ti, tj);
    // make Otau(1)(g12|titj>)
    CCPairFunction<double,6> Otau1_gtt = apply_Ot(gtt, tau, 1);
    // make Otau1Q2 part and the Otau1Otau2. Otau1Otau2 part IS NOT used in the symmetry exploit
    CCPairFunction<double,6> OtauQ = apply_Qt(Otau1_gtt, mo_ket_, 2);
    CCPairFunction<double,6> OtauOtau = apply_Ot(Otau1_gtt, tau, 2);
    // apply the Greens Operator
    const real_function_6d GVcc_1 = -2.0 * apply_G(OtauQ, G);
    const real_function_6d GVcc_3 = -2.0 * apply_G(OtauOtau, G);
    // pair symmetry exploit
    real_function_6d GVcc_2;
    if (symmetric) GVcc_2 = swap_particles(GVcc_1);
    else {
        CCPairFunction<double,6> Otau2_gtt = apply_Ot(gtt, tau, 2);
        CCPairFunction<double,6> QOtau = apply_Qt(Otau2_gtt, mo_ket_, 1);
        GVcc_2 = -2.0 * apply_G(QOtau, G);
    }
    if (parameters.debug()) print_size(GVcc, "GVcc", parameters.debug());

    if (parameters.debug()) GVcc = apply_Q12t(GVcc, mo_ket_);

    if (parameters.debug()) print_size(GVcc, "QGVcc");

    time_Vcc.stop();
    GVcc = GVcc_3 - GVcc_1 - GVcc_2;
    real_function_6d result = GV + GVcc;     // sign is absorbed into GVcc
    if (parameters.debug()) result.print_size("constant-part");

    result = apply_Q12t(result, mo_ket_);
    if (parameters.debug()) result.print_size("constant-part");

    output.section("Constant Term Calculation of Pair " + u.name() + " ended");
    time_Vcc.info(true, GVcc.norm2());
    time_GV.info(true, GV.norm2());
    time.info(true, result.norm2());
    return result;
}

madness::real_function_6d
CCPotentials::make_constant_part_cc2_Qt_gs(const CCPair& u, const CC_vecfunction& tau,
                                           const real_convolution_6d *Gscreen) const {
    output.section("Calculating Constant Part of CC2: Qt-Ansatz");
    CCTimer time(world, "Calculating Constant Part of CC2: Qt-Ansatz");
    MADNESS_ASSERT(parameters.QtAnsatz());
    MADNESS_ASSERT(u.ctype == CT_CC2);
    MADNESS_ASSERT(u.type == GROUND_STATE);
    MADNESS_ASSERT(tau.type == PARTICLE);
    // convenience
    const size_t i = u.i;
    const size_t j = u.j;
    const bool symmetric = (i == j);
    // make ti intermediates
    const CCFunction<double,3> ti = make_t_intermediate(tau(i));
    const CCFunction<double,3> tj = make_t_intermediate(tau(j));
    const CC_vecfunction t = make_full_t_intermediate(tau);
    MADNESS_ASSERT(t.size() == mo_ket_.size());
    MADNESS_ASSERT(t.type == MIXED);
    real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0 * get_epsilon(ti.i, tj.i)), parameters.lo(),
                                           parameters.thresh_bsh_6D());
    G.destructive() = true;
    // calculate [F,Qt] commutator which is [F1,Q1t]Q2t + Q1t [F2,Q2t]
    // and [F1,Q1t] = - [F1,O1t] = - (F-e_k) |tk><k| = - (F-e_k) |tauk><k| = |Vk><k|
    // commutator is applied to f12|titj>
    output.section("Make [F,Qt] commutator");
    CCTimer time_comm(world, "commutator");
    const vector_real_function_3d Vtmp = get_potentials(tau, POT_singles_);
    const CC_vecfunction V(Vtmp, UNDEFINED, parameters.freeze());
    const CCPairFunction<double,6> ftt(f12, ti, tj);
    const CCPairFunction<double,6> O1ftt = apply_Ot(ftt, V, 1);
    const CCPairFunction<double,6> O1Q2ftt = apply_Qt(O1ftt, t, 2);
    const real_function_6d part1 = -2.0 * apply_G(O1Q2ftt, G);
    real_function_6d part2;
    if (symmetric) part2 = swap_particles(part1);
    else {
        const CCPairFunction<double,6> O2ftt = apply_Ot(ftt, V, 2);
        const CCPairFunction<double,6> Q1O2ftt = apply_Qt(O2ftt, t, 1);
        part2 = -2.0 * apply_G(Q1O2ftt, G);
    }
    const real_function_6d commutator = part1 + part2;
    time_comm.info();
    CCTimer time_GV(world, "GV");
    real_function_6d GV;
    {
        real_function_6d V = apply_Vreg(ti, tj, Gscreen);
        print_size(V, "Vreg", parameters.debug());
        V = apply_Q12t(V, t);
        print_size(V, "QVreg");
        GV = -2.0 * apply_G(V, G);
    }
    print_size(GV, "GVreg", parameters.debug());
    GV = apply_Q12t(GV, mo_ket_);
    print_size(GV, "QtGVreg");
    time_GV.info();
    real_function_6d result = GV + commutator;
    print_size(GV, "GVreg");
    print_size(commutator, "[F,Qt]");
    result = apply_Q12t(result, mo_ket_);
    print_size(result, "constant-part");
    time.info();
    return result;
}

madness::real_function_6d
CCPotentials::make_constant_part_cispd(const CCPair& u, const CC_vecfunction& x,
                                       const real_convolution_6d *Gscreen) const {
    output.section("Make Constant Part of " + assign_name(u.ctype) + " for pair " + u.name() + ": Q-Ansatz");
    CCTimer time(world, "Constant Part");
    MADNESS_ASSERT(u.ctype == CT_CISPD || u.ctype == CT_ADC2);
    const size_t i = u.i;
    const size_t j = u.j;
    const CCFunction<double,3>& xi = x(i);
    const CCFunction<double,3>& xj = x(j);
    const CCFunction<double,3>& moi = mo_ket_(i);
    const CCFunction<double,3>& moj = mo_ket_(j);
    const bool symmetric = (i == j);
    MADNESS_ASSERT(x.type == RESPONSE);
    MADNESS_ASSERT(u.bsh_eps == get_epsilon(i, j) + x.omega);
    //MADNESS_ASSERT(not parameters.QtAnsatz());
    if (parameters.QtAnsatz()) output.warning("Demanded Constant Part with Q Ansatz, but parameter QtAnsatz is true");

    // timers
    CCTimer time_fr(world, "Functional Response");
    CCTimer time_cr(world, "Coulomb Response");
    real_function_6d GV;
    real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0 * u.bsh_eps), parameters.lo(), parameters.thresh_bsh_6D());
    G.destructive() = true;
    {
        time_fr.start();
        real_function_6d V;
        const real_function_6d Vxm = apply_Vreg(xi, moj, Gscreen);
        if (symmetric) V = Vxm;
        else V = Vxm + apply_Vreg(moi, xj, Gscreen);

        print_size(V, "Vreg", parameters.debug());
        V = apply_Q12t(V, mo_ket_);
        print_size(V, "QVreg");
        const real_function_6d GVtmp = -2.0 * apply_G(V, G);
        if (symmetric) GV = GVtmp + swap_particles(GVtmp);
        else GV = GVtmp;

        time_fr.stop();
    }     // now the response of the CC2 coupling potential
    // - (O1xQ2 + Q1O2x)g12|ij>
    real_function_6d GVcc;
    {
        time_cr.start();
        const CCPairFunction<double,6> gij(g12, moi, moj);
        const CCPairFunction<double,6> O1x_gij = apply_Ot(gij, x, 1);
        const CCPairFunction<double,6> OQ_part = apply_Qt(O1x_gij, mo_ket_, 2);
        const real_function_6d GOQ = -2.0 * apply_G(OQ_part, G);
        if (symmetric) GVcc = GOQ + swap_particles(GOQ);
        else {
            const CCPairFunction<double,6> O2x_gij = apply_Ot(gij, x, 2);
            const CCPairFunction<double,6> QO_part = apply_Qt(O2x_gij, mo_ket_, 1);
            const real_function_6d GQO = -2.0 * apply_G(QO_part, G);
            GVcc = GOQ + GQO;
        }
        time_cr.stop();
    }
    real_function_6d result = GV - GVcc;
    if (parameters.debug()) print_size(result, "constant-part", parameters.debug());

    result = apply_Q12t(result, mo_ket_);
    if (parameters.debug()) result.print_size("QConstant-Part");

    time_fr.info(true, GV.norm2());
    time_cr.info(true, GVcc.norm2());
    time.info(true, result.norm2());
    return result;
}

madness::real_function_6d
CCPotentials::make_constant_part_cispd_Qt(const CCPair& u, const CC_vecfunction& x,
                                          const real_convolution_6d *Gscreen) const {
    output.section("Make Constant Part of " + assign_name(u.ctype) + " for pair " + u.name() + ": Qt-Ansatz");
    CCTimer time(world, "Constant Part");
    MADNESS_ASSERT(u.ctype == CT_CISPD || u.ctype == CT_ADC2);
    const size_t i = u.i;
    const size_t j = u.j;
    const CCFunction<double,3>& xi = x(i);
    const CCFunction<double,3>& xj = x(j);
    const CCFunction<double,3>& moi = mo_ket_(i);
    const CCFunction<double,3>& moj = mo_ket_(j);
    const bool symmetric = (i == j);
    MADNESS_ASSERT(x.type == RESPONSE);
    MADNESS_ASSERT(u.bsh_eps == get_epsilon(i, j) + x.omega);
    //MADNESS_ASSERT(parameters.QtAnsatz());
    if (!parameters.QtAnsatz()) output.warning("Demanded Constant Part with Qt Ansatz, but parameter QtAnsatz is false");

    CCTimer time_FR(world, "Functional Response");
    CCTimer time_PR(world, "Projector Response");
    CCTimer time_CFR(world, "Commutator Functional Response");
    CCTimer time_CPR(world, "Commutator Projector Response");
    real_function_6d FR;
    real_function_6d PR;
    real_function_6d CFR;
    real_function_6d CPR;
    real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0 * u.bsh_eps), parameters.lo(), parameters.thresh_bsh_6D());
    G.destructive() = true;
    // Response of the Projector from Q12tf12|titj> --> (-O1x - O2x + O1xO2 + O1O2x)Vreg|ij> = -(O1xQ+QO2x)Vreg|ij>
    {
        time_PR.start();
        // here is an inconsistency: The Vreg potential will apply (F-ei) to the two hole states but we have here (F-eij-omega)
        // in the future this part here is supposed to be entirely 3D and not use the 6D apply_Vreg function, so right now this is a workaround
        // however, we have to add the missing -omega*f12|ij>
        real_function_6d Vtmp = apply_Vreg(moi, moj, Gscreen);
        real_function_6d f12omegaij = make_f_xy(moi, moj);
        f12omegaij.scale(x.omega);
        f12omegaij.print_size("omega*f12|ij>");
        Vtmp = Vtmp - f12omegaij;
        Vtmp.truncate().reduce_rank();
        if (parameters.debug()) Vtmp.print_size("Vreg-omega*f12|ij>");

        const CCPairFunction<double,6> V(Vtmp);
        const CCPairFunction<double,6> O1x_V = apply_Ot(V, x, 1);
        const CCPairFunction<double,6> O1xQ2_V = apply_Qt(O1x_V, mo_ket_, 2);
        const real_function_6d GOQ = -2.0 * apply_G(O1xQ2_V, G);
        if (symmetric) PR = GOQ + swap_particles(GOQ);
        else {
            const CCPairFunction<double,6> O2x_V = apply_Ot(V, x, 2);
            const CCPairFunction<double,6> Q1O2x_V = apply_Qt(O2x_V, mo_ket_, 1);
            const real_function_6d GQO = -2.0 * apply_G(Q1O2x_V, G);
            PR = GOQ + GQO;
        }
        time_PR.stop();
    }        //      print_size(GQR,"projector-response",parameters.debug());
    //      GQR = apply_Qt(GQR,mo_ket_);
    //      print_size(GQR,"projector-response",parameters.debug());
    //      time_QR.info(true,GQR.norm2());
    // function response
    {
        time_FR.start();
        real_function_6d V;
        const real_function_6d Vxm = apply_Vreg(xi, moj, Gscreen);
        if (symmetric) V = Vxm;
        else V = Vxm + apply_Vreg(moi, xj, Gscreen);

        print_size(V, "Vreg", parameters.debug());
        V = apply_Q12t(V, mo_ket_);
        print_size(V, "QVreg");
        const real_function_6d GVtmp = -2.0 * apply_G(V, G);
        if (symmetric) FR = GVtmp + swap_particles(GVtmp);
        else FR = GVtmp;

        time_FR.stop();
    }
    //      print_size(GV,"GVreg",parameters.debug());
    //      GV=apply_Qt(GV,mo_ket_);
    //      print_size(GV,"GVreg",parameters.debug());
    //      time_GV.info(true,GV.norm2());
    {
        time_CFR.start();
        CFR = real_factory_6d(world);
        time_CFR.stop();
    }     // make Commutator Projector response Response: [F,d/dtau(Qt)] part of d/dtau{([F,Qt])f12|xitj + tixj>}
    // {-O1x[F,Q2t] - Q1t[F,O2x] - [F,O1x]Q2t - [F,Q1t]O2x , used d/dtau(Qt) = -Ox
    //  O1x[F,O2t] - Q1t[F,O2x] - [F,O1x]Q2t + [F,O1t]O2x ,  used [F,Qt] = -[F,Ot]
    //  -O1x*O2Vt + Q1t*(O2Vx - omega*O2x) + (O1Vx-omega*O1x)Q2t - O1Vt*O2x , used [F,Ot] = -OVt and [F,Ox] = -(OVx - omega*Ox)
    // part1 = 01xO2Vt*f12|titj>
    // part2 = (O1Vx-omega*O1x)Q2t*f12|titj>
    // and then the same for 1 and 2 switched
    // }f12|titj>
    // Vt is zero for CIS(D) since there are no ground state singles
    CCTimer time_cpr(world, "Commutator-Projector Response");
    {
        time_CPR.start();
        const vector_real_function_3d Vxtmp = sub(world, get_potentials(x, POT_singles_),
                                                  x.omega * x.get_vecfunction());
        const CC_vecfunction Vx(Vxtmp, UNDEFINED, parameters.freeze());
        CCPairFunction<double,6> ftt(f12, moi, moj);
        real_function_6d tmp1;
        real_function_6d tmp2;
        {
            CCPairFunction<double,6> OVx = apply_Ot(ftt, Vx, 1);
            CCPairFunction<double,6> OVxQt = apply_Qt(OVx, mo_ket_, 2);
            real_function_6d part2 = -2.0 * apply_G(OVxQt, G);
            tmp1 = part2;
        }
        if (symmetric) tmp2 = swap_particles(tmp1);
        else {
            CCPairFunction<double,6> OVx = apply_Ot(ftt, Vx, 2);
            CCPairFunction<double,6> OVxQt = apply_Qt(OVx, mo_ket_, 1);
            real_function_6d part2 = -2.0 * apply_G(OVxQt, G);
            tmp2 = part2;
        }
        CPR = tmp1 + tmp2;
        time_CPR.stop();
    }
    real_function_6d result = FR - PR + CFR + CPR;
    print_size(result, "constant-part", parameters.debug());
    result = apply_Q12t(result, mo_ket_);
    print_size(result, "constant-part", parameters.debug());
    output.section("Constant Term for Pair " + u.name() + " ended");
    time_FR.info(true, FR.norm2());
    time_PR.info(true, PR.norm2());
    time_CFR.info(true, CFR.norm2());
    time_CPR.info(true, CPR.norm2());
    result.print_size("Constant Term");
    time.info(true, result.norm2());
    return result;
}

madness::real_function_6d
CCPotentials::make_constant_part_cc2_ex(const CCPair& u, const CC_vecfunction& tau, const CC_vecfunction& x,
                                        const real_convolution_6d *Gscreen) {
    output.section("Make Constant Part of " + assign_name(u.ctype) + " for pair " + u.name() + ": Q-Ansatz");
    CCTimer time(world, "Constant Part");
    MADNESS_ASSERT(tau.type == PARTICLE);
    MADNESS_ASSERT(x.type == RESPONSE);
    MADNESS_ASSERT(u.type == EXCITED_STATE);
    MADNESS_ASSERT(u.ctype == CT_LRCC2);
    if (parameters.QtAnsatz()) output.warning("Demanded Constant Part with Q Ansatz, but parameter QtAnsatz is true");

    // convenience
    const size_t i = u.i;
    const size_t j = u.j;
    const CC_vecfunction t = make_full_t_intermediate(tau);
    const CCFunction<double,3>& ti = t(i);
    const CCFunction<double,3>& tj = t(j);
    const CCFunction<double,3>& xi = x(i);
    const CCFunction<double,3>& xj = x(j);
    // use pair symmetry for diagonal pairs
    const bool symmetric = (i == j);
    // the Greens operator
    real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0 * u.bsh_eps), parameters.lo(), parameters.thresh_bsh_6D());
    G.destructive() = true;
    // data and output
    CCTimer time_fr(world, "Functional Response");
    CCTimer time_cr(world, "Coulomb Response");
    // functional response part: G(Q*Vreg*|(xitj + tixj)>)
    real_function_6d functional_response;
    {
        time_fr.start();
        const real_function_6d Vxt = (apply_Vreg(xi, ti, Gscreen)).truncate().reduce_rank();
        if (symmetric) {
            real_function_6d V = apply_Q12t(Vxt, mo_ket_);
            const real_function_6d tmp = -2.0 * G(V);
            functional_response = tmp + swap_particles(tmp);
        } else {
            const real_function_6d Vtx = apply_Vreg(ti, xj, Gscreen);
            real_function_6d V = (Vtx + Vxt).truncate().reduce_rank();
            V = apply_Q12t(V, mo_ket_);
            functional_response = -2.0 * G(V);
        }
        time_fr.stop();
    }      // Coulomb Response part, theoverall minus sign is applied in the end
    // Functional response:
    // -(QOtau + OtauQ - OtauOtau)g12|xitj+tixj> (part1: QOtau, part2, OtauQ, part3 OtauOtau)
    // Projector response:
    // -(QtOx+OxQt)g12|titj> (QtOx is part 4, OxQt part is for second part)
    //
    // the overall minus sign is added in the end
    real_function_6d coulomb_response;
    {
        time_cr.start();
        real_function_6d tmp1;
        real_function_6d tmp2;
        // make the xt parts of the functional and the QtOx part of the projector response
        {
            CCPairFunction<double,6> gxt(g12, xi, tj);
            // make QOtau*g*|xt>
            CCPairFunction<double,6> O2tmp = apply_Ot(gxt, tau, 2);
            CCPairFunction<double,6> QO = apply_Qt(O2tmp, mo_ket_, 1);
            const real_function_6d part1 = -2.0 * apply_G(QO, G);
            // make OtauQ*g*|xt>
            CCPairFunction<double,6> O1tmp = apply_Ot(gxt, tau, 1);
            CCPairFunction<double,6> OQ = apply_Qt(O1tmp, mo_ket_, 2);
            const real_function_6d part2 = -2.0 * apply_G(OQ, G);
            // OtauOtau*g*|xt>
            CCPairFunction<double,6> OO = apply_Ot(O1tmp, tau, 2);
            const real_function_6d part3 = -2.0 * apply_G(OO, G);
            // QtOx*g|titj>
            CCPairFunction<double,6> gtt(g12, ti, tj);
            CCPairFunction<double,6> O2x = apply_Ot(gtt, x, 2);
            CCPairFunction<double,6> QtOx = apply_Qt(O2x, t, 1);
            const real_function_6d part4 = -2.0 * apply_G(QtOx, G);
            tmp1 = part1 + part2 - part3 + part4;     // overall minus sign applied in the end
        }
        if (symmetric) tmp2 = swap_particles(tmp1);
        else {
            CCPairFunction<double,6> gtx(g12, ti, xj);
            // make QOtau*g*|tx>
            CCPairFunction<double,6> O2tmp = apply_Ot(gtx, tau, 2);
            CCPairFunction<double,6> QO = apply_Qt(O2tmp, mo_ket_, 1);
            const real_function_6d part1 = -2.0 * apply_G(QO, G);
            // make OtauQ*g*|tx>
            CCPairFunction<double,6> O1tmp = apply_Ot(gtx, tau, 1);
            CCPairFunction<double,6> OQ = apply_Qt(O1tmp, mo_ket_, 2);
            const real_function_6d part2 = -2.0 * apply_G(OQ, G);
            // OtauOtau*g*|tx>
            CCPairFunction<double,6> OO = apply_Ot(O1tmp, tau, 2);
            const real_function_6d part3 = -2.0 * apply_G(OO, G);
            // OxQt*g|titj>
            CCPairFunction<double,6> gtt(g12, ti, tj);
            CCPairFunction<double,6> O1x = apply_Ot(gtt, x, 1);
            CCPairFunction<double,6> OxQt = apply_Qt(O1x, t, 2);
            const real_function_6d part4 = -2.0 * apply_G(OxQt, G);
            tmp1 = part1 + part2 - part3 + part4;     // overall minus sign applied in the end
        }
        coulomb_response = tmp1 + tmp2;
        time_cr.stop();
    }
    real_function_6d result = functional_response - coulomb_response;
    result = apply_Q12t(result, mo_ket_);
    output.section("Constant Term for Pair " + u.name() + " ended");
    if (parameters.debug()) functional_response.print_size("Functional Response");

    if (parameters.debug()) coulomb_response.print_size("Coulomb Response");

    if (parameters.debug()) result.print_size("Constant Term");

    time_fr.info(true, functional_response.norm2());
    time_cr.info(true, coulomb_response.norm2());
    time.stop();
    time.info(true, result.norm2());
    return result;
}

madness::real_function_6d
CCPotentials::make_constant_part_cc2_Qt_ex(const CCPair& u, const CC_vecfunction& tau, const CC_vecfunction& x,
                                           const real_convolution_6d *Gscreen) {
    output.section("Make Constant Part of " + assign_name(u.ctype) + " for pair " + u.name() + ": Qt-Ansatz");
    MADNESS_ASSERT(tau.type == PARTICLE);
    MADNESS_ASSERT(x.type == RESPONSE);
    MADNESS_ASSERT(u.type == EXCITED_STATE);
    MADNESS_ASSERT(u.ctype == CT_LRCC2);
    if (!parameters.QtAnsatz()) output.warning("Demanded Constant Part with Qt Ansatz, but parameter QtAnsatz is false");

    // convenience
    const size_t i = u.i;
    const size_t j = u.j;
    const CC_vecfunction t = make_full_t_intermediate(tau);
    MADNESS_ASSERT(t.type == MIXED);
    MADNESS_ASSERT(t.size() == mo_ket_.size());
    const CCFunction<double,3>& ti = t(i);
    const CCFunction<double,3>& tj = t(j);
    const CCFunction<double,3>& xi = x(i);
    const CCFunction<double,3>& xj = x(j);
    // use pair symmetry for diagonal pairs
    const bool symmetric = (i == j);
    // the Greens operator
    real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0 * u.bsh_eps), parameters.lo(), parameters.thresh_bsh_6D());
    G.destructive() = true;
    // data and output
    CCTimer time_fr(world, "Functional Response");
    CCTimer time_pr(world, "Projector Response");
    CCTimer time_cr(world, "Commutator Response");
    CCTimer time_cpr(world, "Commutator-Projector Response");
    // Make functional response part: G(QtVreg|xitj + tixj>)
    real_function_6d functional_response;
    {
        time_fr.start();
        const real_function_6d Vxt = (apply_Vreg(xi, tj, Gscreen)).truncate().reduce_rank();
        if (symmetric) {
            real_function_6d V = apply_Q12t(Vxt, t);
            const real_function_6d tmp = -2.0 * G(V);
            functional_response = tmp + swap_particles(tmp);
        } else {
            const real_function_6d Vtx = apply_Vreg(ti, xj, Gscreen);
            real_function_6d V = (Vtx + Vxt).truncate().reduce_rank();
            V = apply_Q12t(V, t);
            functional_response = -2.0 * G(V);
        }
        time_fr.stop();
    }     // make Projector Response: -G(OxQt+QtOx)Vreg|titj>
    real_function_6d projector_response;
    {
        time_pr.start();
        // here is an inconsistency: The Vreg potential will apply (F12-eij) to the |titj> state but we have here (F12-eij-omega)
        // in the future this part here is supposed to be entirely 3D and not use the 6D apply_Vreg function, so right now this is a workaround
        // however, we have to add the missing -omega|titj>
        real_function_6d Vtt_tmp = apply_Vreg(ti, tj, Gscreen);
        real_function_6d titj = make_f_xy(ti, tj);
        Vtt_tmp = Vtt_tmp - x.omega * titj;
        CCPairFunction<double,6> Vtt(Vtt_tmp);
        real_function_6d tmp1;
        real_function_6d tmp2;
        {
            CCPairFunction<double,6> Ox = apply_Ot(Vtt, x, 1);
            CCPairFunction<double,6> OxQt = apply_Qt(Ox, t, 2);
            tmp1 = -2.0 * apply_G(OxQt, G);
        }
        if (symmetric) tmp2 = swap_particles(tmp1);
        else {
            CCPairFunction<double,6> Ox = apply_Ot(Vtt, x, 2);
            CCPairFunction<double,6> QtOx = apply_Qt(Ox, t, 1);
            tmp2 = -2.0 * apply_G(QtOx, G);
        }
        projector_response = tmp1 + tmp2;
        time_pr.stop();
    }
    // make commutator response: [F12,Qt12]f12|xitj+tixj> = (O1VQ2t + Q1tO2V)f12|xitj+tixj>
    real_function_6d commutator_response;
    {
        time_cr.start();
        real_function_6d part1;     // the xt parts
        const vector_real_function_3d Vtmp = get_potentials(tau, POT_singles_);
        const CC_vecfunction V(Vtmp, UNDEFINED, parameters.freeze());
        {
            const CCPairFunction<double,6> fxt(f12, xi, tj);
            const CCPairFunction<double,6> O1V = apply_Ot(fxt, V, 1);
            const CCPairFunction<double,6> OQ = apply_Qt(O1V, t, 2);
            const CCPairFunction<double,6> O2V = apply_Ot(fxt, V, 2);
            const CCPairFunction<double,6> QO = apply_Qt(O2V, t, 1);
            const real_function_6d tmp1 = -2.0 * apply_G(OQ, G);
            const real_function_6d tmp2 = -2.0 * apply_G(QO, G);
            part1 = tmp1 + tmp2;
        }
        real_function_6d part2;     // the tx parts
        if (symmetric) part2 = swap_particles(part1);
        else {
            const CCPairFunction<double,6> ftx(f12, ti, xj);
            const CCPairFunction<double,6> O1V = apply_Ot(ftx, V, 1);
            const CCPairFunction<double,6> OQ = apply_Qt(O1V, t, 2);
            const CCPairFunction<double,6> O2V = apply_Ot(ftx, V, 2);
            const CCPairFunction<double,6> QO = apply_Qt(O2V, t, 1);
            const real_function_6d tmp1 = -2.0 * apply_G(OQ, G);
            const real_function_6d tmp2 = -2.0 * apply_G(QO, G);
            part2 = tmp1 + tmp2;
        }
        commutator_response = part1 + part2;
        time_cr.stop();
    }
    // make Commutator Projector response Response: [F,d/dtau(Qt)] part of d/dtau{([F,Qt])f12|xitj + tixj>}
    // {-O1x[F,Q2t] - Q1t[F,O2x] - [F,O1x]Q2t - [F,Q1t]O2x , used d/dtau(Qt) = -Ox
    //  O1x[F,O2t] - Q1t[F,O2x] - [F,O1x]Q2t + [F,O1t]O2x ,  used [F,Qt] = -[F,Ot]
    //  -O1x*O2Vt + Q1t*(O2Vx - omega*O2x) + (O1Vx-omega*O1x)Q2t - O1Vt*O2x , used [F,Ot] = -OVt and [F,Ox] = -(OVx - omega*Ox)
    // part1 = 01xO2Vt*f12|titj>
    // part2 = (O1Vx-omega*O1x)Q2t*f12|titj>
    // and then the same for 1 and 2 switched
    // }f12|titj>
    real_function_6d commutator_projector_response;
    {
        time_cpr.start();
        const vector_real_function_3d Vxtmp = sub(world, get_potentials(x, POT_singles_),
                                                  x.omega * x.get_vecfunction());
        const vector_real_function_3d Vttmp = get_potentials(tau, POT_singles_);
        const CC_vecfunction Vx(Vxtmp, UNDEFINED, parameters.freeze());
        const CC_vecfunction Vt(Vttmp, UNDEFINED, parameters.freeze());
        CCPairFunction<double,6> ftt(f12, ti, tj);
        real_function_6d tmp1;
        real_function_6d tmp2;
        {
            CCPairFunction<double,6> Ox = apply_Ot(ftt, x, 1);
            CCPairFunction<double,6> OxOVt = apply_Ot(Ox, Vt, 2);
            real_function_6d part1 = -2.0 * apply_G(OxOVt, G);
            CCPairFunction<double,6> OVx = apply_Ot(ftt, Vx, 1);
            CCPairFunction<double,6> OVxQt = apply_Qt(OVx, t, 2);
            real_function_6d part2 = -2.0 * apply_G(OVxQt, G);
            tmp1 = part2 - part1;
        }
        if (symmetric) tmp2 = swap_particles(tmp1);
        else {
            CCPairFunction<double,6> Ox = apply_Ot(ftt, x, 2);
            CCPairFunction<double,6> OxOVt = apply_Ot(Ox, Vt, 1);
            real_function_6d part1 = -2.0 * apply_G(OxOVt, G);
            CCPairFunction<double,6> OVx = apply_Ot(ftt, Vx, 2);
            CCPairFunction<double,6> OVxQt = apply_Qt(OVx, t, 1);
            real_function_6d part2 = -2.0 * apply_G(OVxQt, G);
            tmp2 = part2 - part1;
        }
        commutator_projector_response = tmp1 + tmp2;
        time_cpr.stop();
    }
    real_function_6d result =
            functional_response - projector_response + commutator_response + commutator_projector_response;
    result = apply_Q12t(result, mo_ket_);
    output.section("Constant Term for Pair " + u.name() + " ended");
    time_fr.info(true, functional_response.norm2());
    time_pr.info(true, projector_response.norm2());
    time_cr.info(true, commutator_response.norm2());
    time_cpr.info(true, commutator_projector_response.norm2());
    result.print_size("Constant Term");
    return result;
}

madness::real_function_6d
CCPotentials::apply_Vreg(const CCFunction<double,3>& ti, const CCFunction<double,3>& tj, const real_convolution_6d *Gscreen) const {
    output("Applying Vreg to |" + ti.name() + tj.name() + ">");
    CCTimer timer(world, "Vreg|" + ti.name() + tj.name() + ">");
    CCTimer time_f(world, "F-Part");
    const real_function_6d F_part = apply_reduced_F(ti, tj, Gscreen);
    time_f.stop();
    CCTimer time_u(world, "U-Part");
    const real_function_6d U_part = apply_transformed_Ue(ti, tj, Gscreen);
    time_u.stop();
    CCTimer time_k(world, "K-Part");
    const real_function_6d K_part = apply_exchange_commutator(ti, tj);     // maybe use screening later
    time_k.stop();
    const real_function_6d result = F_part + U_part - K_part;
    if (parameters.debug()) F_part.print_size("F-Part");

    if (parameters.debug()) U_part.print_size("U-Part");

    if (parameters.debug()) K_part.print_size("K-Part");

    time_f.info(true, F_part.norm2());
    time_u.info(true, U_part.norm2());
    time_k.info(true, K_part.norm2());
    if (parameters.debug()) result.print_size("Vreg|" + ti.name() + tj.name() + ">");

    timer.info(true, result.norm2());
    return result;
}

madness::real_function_6d
CCPotentials::apply_Vreg_macrotask(World& world, const std::vector<real_function_3d>& mo_ket,
                                   const std::vector<real_function_3d>& mo_bra,
                                   const CCParameters& parameters, const real_function_3d& Rsquare,
                                   const std::vector<real_function_3d>& U1, const size_t& i, const size_t& j,
                                   const FuncType& x_type, const FuncType& y_type, const std::vector<std::string> argument,
                                   const real_convolution_6d *Gscreen) {
    const real_function_3d& x_ket = mo_ket[i];
    const real_function_3d& y_ket = mo_ket[j];
    const real_function_3d& x_bra = mo_bra[i];
    const real_function_3d& y_bra = mo_bra[j];

    const std::string x_name = "phi" + stringify(i);
    const std::string y_name = "phi" + stringify(j);

    // print("Applying Vreg to |" + x_name + y_name + ">\n");
    CCTimer timer(world, "Vreg|" + x_name + y_name + ">");
    //CCTimer time_f(world, "F-Part");
    //const real_function_6d F_part = apply_reduced_F(ti, tj, Gscreen);
    //time_f.stop();
    real_function_6d result = real_factory_6d(world);
    bool do_Ue=std::find(argument.begin(), argument.end(), "Ue") != argument.end();
    bool do_KffK=std::find(argument.begin(), argument.end(), "KffK") != argument.end();
    bool do_f12phi=std::find(argument.begin(), argument.end(), "f12phi") != argument.end();
    MADNESS_CHECK(do_Ue or do_KffK or do_f12phi);
    MADNESS_CHECK(do_f12phi xor (do_Ue or do_KffK));
    if (do_Ue) {
        CCTimer time_u(world, "U-Part");
        const real_function_6d U_part = apply_transformed_Ue_macrotask(world, mo_ket, parameters, Rsquare,
                                                                       U1, i, j, x_type, y_type, Gscreen);
        if (parameters.debug()) U_part.print_size("U-Part");
        result+=U_part;
        time_u.stop();
        time_u.info(true, U_part.norm2());
    }
    if (do_KffK) {
        CCTimer time_k(world, "K-Part");
        const real_function_6d K_part = apply_exchange_commutator_macrotask(world, mo_ket, mo_bra, Rsquare,
                                                                            i, j, parameters, x_type, y_type, Gscreen);
        if (parameters.debug()) K_part.print_size("K-Part");
        result-=K_part;
        time_k.stop();
        time_k.info(true, K_part.norm2());
    }
    if (do_f12phi) {
        CCTimer time_f(world, "f12phi-Part");
        const real_function_6d f12phi = make_f_xy_macrotask(world, x_ket, y_ket, x_bra, y_bra, i, j, parameters, x_type, y_type, Gscreen);
        if (parameters.debug()) f12phi.print_size("f12phi-Part");
        result+=f12phi;
        time_f.stop();
        time_f.info(true, f12phi.norm2());
    }

    if (parameters.debug()) result.print_size("Vreg|" + x_name + y_name + ">");

    timer.info(true, result.norm2());
    return result;
}

madness::real_function_6d
CCPotentials::apply_reduced_F(const CCFunction<double,3>& ti, const CCFunction<double,3>& tj, const real_convolution_6d *Gscreen) const {
    //CC_Timer time(world,"(F-eij)|"+ti.name()+tj.name()+">");
    // get singles potential
    const bool symmetric = (ti.type == tj.type && ti.i == tj.i);
    const real_function_3d Vti = get_potentials(ti, POT_singles_);
    const real_function_3d Vtj = get_potentials(tj, POT_singles_);
    const real_function_6d Vt = make_f_xy(Vti, tj, Gscreen);
    real_function_6d tV;
    if (symmetric) tV = swap_particles(Vt);
    else tV = make_f_xy(ti, Vtj, Gscreen);

    const real_function_6d result = -1.0 * (Vt + tV);
    //result.print_size("(F-eij)|"+ti.name()+tj.name()+">");
    //time.info();
    return result;
}


madness::real_function_6d
CCPotentials::apply_transformed_Ue(const CCFunction<double,3>& x, const CCFunction<double,3>& y, const real_convolution_6d *Gscreen) const {
    if (parameters.debug()) output("\nComputing Ue|" + x.name() + y.name() + ">\n");

    const bool symmetric = (x.type == y.type && x.i == y.i);
    CCTimer time_Ue(world, "Ue|" + x.name() + y.name() + ">");
    double tight_thresh = parameters.thresh_6D();     // right now this is the std. thresh
    // check if screening operator is in modified NS Form
    if (Gscreen != NULL) {
        if (!Gscreen->modified()) error("Demanded Screening for Ue but given BSH Operator is not in modified NS form");
    }
    if (parameters.debug()) output("Applying transformed Ue to " + x.name() + y.name());

    if (parameters.debug() && symmetric) output("Exploiting Pair Symmetry");

    real_function_6d Uxy = real_factory_6d(world);
    Uxy.set_thresh(tight_thresh);
    // Apply the untransformed U Potential
    Uxy = corrfac.apply_U(x.function, y.function, *Gscreen, symmetric);
    Uxy.set_thresh(tight_thresh);
    // Apply the double commutator R^{-1}[[T,f,R]
    for (size_t axis = 0; axis < 3; axis++) {
        // Make the local parts of the Nuclear and electronic U potentials
        const real_function_3d Un_local = nemo_->ncf->U1(axis);
        const real_function_3d Un_local_x = (Un_local * x.function).truncate();
        real_function_3d Un_local_y;
        if (symmetric) Un_local_y = copy(Un_local_x);
        else Un_local_y = (Un_local * y.function).truncate();

        const real_function_6d Ue_local = corrfac.U1(axis);
        // Now add the Un_local_x part to the first particle of the Ue_local potential
        real_function_6d UeUnx = CompositeFactory<double, 6, 3>(world).g12(Ue_local).particle1(Un_local_x).particle2(
                copy(y.function)).thresh(tight_thresh);
        // Fill the Tree where it will be necessary
        UeUnx.fill_cuspy_tree(*Gscreen);
        // Set back the thresh
        UeUnx.set_thresh(FunctionDefaults<6>::get_thresh());
        print_size(UeUnx, "UeUnx", parameters.debug());
        // Now add the Un_local_y part to the second particle of the Ue_local potential
        real_function_6d UeUny;
        if (symmetric) UeUny = -1.0 * swap_particles(UeUnx);     // Ue_local is antisymmetric
        else {
            UeUny = CompositeFactory<double, 6, 3>(world).g12(Ue_local).particle1(copy(x.function)).particle2(
                    Un_local_y).thresh(tight_thresh);
            // Fill the Tree were it will be necessary
            UeUny.fill_cuspy_tree(*Gscreen);
            // Set back the thresh
            UeUny.set_thresh(FunctionDefaults<6>::get_thresh());
        }
        print_size(UeUny, "UeUny", parameters.debug());
        // Construct the double commutator part and add it to the Ue part
        real_function_6d diff = (UeUnx - UeUny).scale(-1.0);
        diff.truncate();
        Uxy = (Uxy + diff).truncate();
    }
    if (parameters.debug()) time_Ue.info();

    // sanity check: <xy|R2 [T,g12] |xy> = <xy |R2 U |xy> - <xy|R2 g12 | xy> = 0
    CCTimer time_sane(world, "Ue-Sanity-Check");
    real_function_6d tmp = CompositeFactory<double, 6, 3>(world).particle1(
            copy(x.function * nemo_->ncf->square())).particle2(copy(y.function * nemo_->ncf->square()));
    const double a = inner(Uxy, tmp);
    const real_function_3d xx = (x.function * x.function * nemo_->ncf->square());
    const real_function_3d yy = (y.function * y.function * nemo_->ncf->square());
    const real_function_3d gxx = (*g12)(xx);
    const double aa = inner(yy, gxx);
    const double error = std::fabs(a - aa);
    const double diff = a - aa;
    time_sane.info(parameters.debug(), error);
    if (world.rank() == 0) {
        std::cout << std::fixed << std::setprecision(10) << "<" << x.name() + y.name() << "|U_R|" << x.name() + y.name()
                  << "> =" << a << ", <" << x.name() + y.name() << "|g12|" << x.name() + y.name()
                  << "> =" << aa << ", diff=" << error << "\n";
        //printf("<xy| U_R |xy>  %12.8f\n",a);
        //printf("<xy|1/r12|xy>  %12.8f\n",aa);
        if (error > FunctionDefaults<6>::get_thresh() * 10.0) output.warning("Ue Potential plain wrong!");
        else if (error > FunctionDefaults<6>::get_thresh()) output.warning("Ue Potential wrong!!!!");
        else output("Ue seems to be sane, diff=" + std::to_string(diff));
    }
    return Uxy;
}


madness::real_function_6d
CCPotentials::apply_transformed_Ue_macrotask(World& world, const std::vector<real_function_3d>& mo_ket,
                                             const CCParameters& parameters, const real_function_3d& Rsquare,
                                             const std::vector<real_function_3d>& U1, const size_t& i, const size_t& j,
                                             const FuncType& x_type, const FuncType& y_type, const real_convolution_6d *Gscreen) {
    const std::string x_name = "phi" + stringify(i);
    const std::string y_name = "phi" + stringify(j);

    if (parameters.debug()) print("Computing Ue|" + x_name + y_name + ">");

    real_function_3d x_function=mo_ket[i];
    real_function_3d y_function=mo_ket[j];
    CorrelationFactor corrfac(world, parameters.gamma(), 1.e-7, parameters.lo());

    const bool symmetric = (x_type == y_type && i == j);
    CCTimer time_Ue(world, "Ue|" + x_name + y_name + ">");
    double tight_thresh = parameters.thresh_6D();     // right now this is the std. thresh
    // check if screening operator is in modified NS Form
    if (Gscreen != NULL) {
        if (!Gscreen->modified()) error("Demanded Screening for Ue but given BSH Operator is not in modified NS form");
    }
    if (parameters.debug()) print("Applying transformed Ue to \n" + x_name + y_name);

    if (parameters.debug() && symmetric) print("Exploiting Pair Symmetry\n");

    real_function_6d Uxy = real_factory_6d(world);
    Uxy.set_thresh(tight_thresh);
    // Apply the untransformed U Potential
    Uxy = corrfac.apply_U(x_function, y_function, *Gscreen, symmetric);
    Uxy.set_thresh(tight_thresh);
    // Apply the double commutator R^{-1}[[T,f,R]
    for (size_t axis = 0; axis < 3; axis++) {
        // Make the local parts of the Nuclear and electronic U potentials
        const real_function_3d Un_local = U1[axis];
        const real_function_3d Un_local_x = (Un_local * x_function).truncate();
        real_function_3d Un_local_y;
        if (symmetric) Un_local_y = copy(Un_local_x);
        else Un_local_y = (Un_local * y_function).truncate();

        const real_function_6d Ue_local = corrfac.U1(axis);
        // Now add the Un_local_x part to the first particle of the Ue_local potential
        real_function_6d UeUnx = CompositeFactory<double, 6, 3>(world).g12(Ue_local).particle1(Un_local_x).particle2(
                copy(y_function)).thresh(tight_thresh);
        // Fill the Tree where it will be necessary
        UeUnx.fill_cuspy_tree(*Gscreen);
        // Set back the thresh
        UeUnx.set_thresh(FunctionDefaults<6>::get_thresh());
//        print_size(UeUnx, "UeUnx", parameters.debug());
        // Now add the Un_local_y part to the second particle of the Ue_local potential
        real_function_6d UeUny;
        if (symmetric) UeUny = -1.0 * madness::swap_particles(UeUnx);     // Ue_local is antisymmetric
        else {
            UeUny = CompositeFactory<double, 6, 3>(world).g12(Ue_local).particle1(copy(x_function)).particle2(
                    Un_local_y).thresh(tight_thresh);
            // Fill the Tree were it will be necessary
            UeUny.fill_cuspy_tree(*Gscreen);
            // Set back the thresh
            UeUny.set_thresh(FunctionDefaults<6>::get_thresh());
        }
//        print_size(UeUny, "UeUny", parameters.debug());
        // Construct the double commutator part and add it to the Ue part
        real_function_6d diff = (UeUnx - UeUny).scale(-1.0);
        diff.truncate();
        Uxy = (Uxy + diff).truncate();
    }
    if (parameters.debug()) time_Ue.info();

    // sanity check: <xy|R2 [T,g12] |xy> = <xy |R2 U |xy> - <xy|R2 g12 | xy> = 0
    CCTimer time_sane(world, "Ue-Sanity-Check");
    real_function_6d tmp = CompositeFactory<double, 6, 3>(world).particle1(
            copy(x_function * Rsquare)).particle2(copy(y_function * Rsquare));
    const double a = inner(Uxy, tmp);
    const real_function_3d xx = (x_function * x_function * Rsquare);
    const real_function_3d yy = (y_function * y_function * Rsquare);
//    const real_function_3d gxx = g12(xx);
    real_convolution_3d poisson= CoulombOperator(world,parameters.lo(),parameters.thresh_3D());
    const real_function_3d gxx= poisson(xx);

    const double aa = inner(yy, gxx);
    const double error = std::fabs(a - aa);
    const double diff = a - aa;
    time_sane.info(parameters.debug(), error);
    if (world.rank() == 0) {
        std::cout << std::fixed << std::setprecision(10) << "<" << x_name + y_name << "|U_R|" << x_name + y_name
                  << "> =" << a << ", <" << x_name + y_name << "|g12|" << x_name + y_name
                  << "> =" << aa << ", diff=" << error << "\n";
        //printf("<xy| U_R |xy>  %12.8f\n",a);
        //printf("<xy|1/r12|xy>  %12.8f\n",aa);
        if (error > FunctionDefaults<6>::get_thresh() * 10.0) std::cout << ("Ue Potential plain wrong!\n");
        else if (error > FunctionDefaults<6>::get_thresh()) std::cout << ("Ue Potential wrong!!!!\n");
        else std::cout << ("Ue seems to be sane, diff=" + std::to_string(diff)) << std::endl;
    }
    return Uxy;
}

madness::real_function_6d
CCPotentials::apply_exchange_commutator_macrotask(World& world, const std::vector<real_function_3d>& mo_ket,
                                                  const std::vector<real_function_3d>& mo_bra, const real_function_3d& Rsquare,
                                                  const size_t& i, const size_t& j, const CCParameters& parameters,
                                                  const FuncType& x_type, const FuncType& y_type,
                                                  const real_convolution_6d *Gscreen) {
    real_function_3d x_ket = mo_ket[i];
    real_function_3d y_ket = mo_ket[j];
    real_function_3d x_bra = mo_bra[i];
    real_function_3d y_bra = mo_bra[j];
    const std::string x_name = "phi" + stringify(i);
    const std::string y_name = "phi" + stringify(j);

    //apply Kf
    if (parameters.debug()) print("\nComputing [K,f]|" + x_name + y_name + ">\n");

    CCTimer time(world, "[K,f]|" + x_name + y_name + ">");
    CCTimer part1_time(world, "Kf" + x_name + y_name + ">");

    bool symmetric_kf = false;
    if ((x_type == y_type) && (i == j)) symmetric_kf = true;

    // First make the 6D function f12|x,y>
    real_function_6d f12xy = make_f_xy_macrotask(world, x_ket, y_ket, x_bra, y_bra, i, j, parameters, x_type, y_type, Gscreen);
    f12xy.truncate().reduce_rank();
    // Apply the Exchange Operator
    real_function_6d Kfxy = K_macrotask(world, mo_ket, mo_bra, f12xy, symmetric_kf, parameters);

    if (parameters.debug()) part1_time.info();

    //apply fk
    CCTimer part2_time(world, "fK" + x_name + y_name + ">");

    const bool symmetric_fk = (x_type == y_type && i == j);
    const real_function_3d Kx = K_macrotask(world, mo_ket, mo_bra, x_ket, parameters);
    const FuncType Kx_type = UNDEFINED;
    const real_function_6d fKphi0b = make_f_xy_macrotask(world, Kx, y_ket, x_bra, y_bra, i, j, parameters, Kx_type, y_type, Gscreen);
    real_function_6d fKphi0a;
    if (symmetric_fk) fKphi0a = madness::swap_particles(fKphi0b);
    else {
        real_function_3d Ky = K_macrotask(world, mo_ket, mo_bra, y_ket, parameters);
        const FuncType Ky_type = UNDEFINED;
        fKphi0a = make_f_xy_macrotask(world, x_ket, Ky, x_bra, y_bra, i, j, parameters, x_type, Ky_type, Gscreen);
    }
    const real_function_6d fKxy = (fKphi0a + fKphi0b);

    if (parameters.debug()) part2_time.info();

    //final result
    Kfxy.print_size("Kf" + x_name + y_name);
    Kfxy.set_thresh(parameters.thresh_6D());
    Kfxy.truncate().reduce_rank();
    Kfxy.print_size("Kf after truncation" + x_name + y_name);
    fKxy.print_size("fK" + x_name + y_name);
    real_function_6d result = (Kfxy - fKxy);
    result.set_thresh(parameters.thresh_6D());
    result.print_size("[K,f]" + x_name + y_name);
    result.truncate().reduce_rank();
    result.print_size("[K,f]" + x_name + y_name);

    //sanity check
    CCTimer sanity(world, "[K,f] sanity check");
    // make the <xy| bra state which is <xy|R2
    const real_function_3d brax = (x_ket * Rsquare);
    const real_function_3d bray = (y_ket * Rsquare);
    real_function_3d xres = result.project_out(brax, 0);
    const double test = bray.inner(xres);
    const double diff = test;
    if (world.rank() == 0) {
        std::cout << std::fixed << std::setprecision(10)
                  << "<" << x_name << y_name << "[K,f]" << x_name << y_name << "> =" << test << "\n";
    }
    if (world.rank() == 0 && fabs(diff) > parameters.thresh_6D()) print("Exchange Commutator Plain Wrong");
    else print("Exchange Commutator seems to be sane, diff=" + std::to_string(diff));

    if (parameters.debug()) sanity.info(diff);

    if (parameters.debug()) print("\n");

    return result;
}


madness::real_function_6d
CCPotentials::apply_exchange_commutator(const CCFunction<double,3>& x, const CCFunction<double,3>& y,
                                        const real_convolution_6d *Gscreen) const {
    real_function_6d result = apply_exchange_commutator1(x, y, Gscreen);
    {
        CCTimer sanity(world, "[K,f] sanity check");
        // make the <xy| bra state which is <xy|R2
        const real_function_3d brax = (x.function * nemo_->ncf->square()).truncate();
        const real_function_3d bray = (y.function * nemo_->ncf->square()).truncate();
        real_function_3d xres = result.project_out(brax, 0);
        const double test = bray.inner(xres);
        const double diff = test;
        if (world.rank() == 0) {
            std::cout << std::fixed << std::setprecision(10)
                      <<     //	  << "<" << x.name() << y.name() << "|fK|" << x.name() << y.name() << "> =" << xyfKxy << ", "
                      //	  << "<" << x.name() << y.name() << "|Kf|" << x.name() << y.name() << "> =" << xyKfxy << ", diff=" << diff << "\n";
                      "<" << x.name() << y.name() << "[K,f]" << x.name() << y.name() << "> =" << test << "\n";
        }
        if (world.rank() == 0 && fabs(diff) > parameters.thresh_6D()) output.warning("Exchange Commutator Plain Wrong");
        else output("Exchange Commutator seems to be sane, diff=" + std::to_string(diff));

        if (parameters.debug()) sanity.info(diff);
    }
    if (parameters.debug()) output("\n");

    return result;
}

madness::real_function_6d
CCPotentials::apply_exchange_commutator1(const CCFunction<double,3>& x, const CCFunction<double,3>& y,
                                         const real_convolution_6d *Gscreen) const {
    if (parameters.debug()) output("\nComputing [K,f]|" + x.name() + y.name() + ">\n");

    CCTimer time(world, "[K,f]|" + x.name() + y.name() + ">");
    // make first part of commutator
    CCTimer part1_time(world, "Kf" + x.name() + y.name() + ">");
    real_function_6d Kfxy = apply_Kf(x, y);
    if (parameters.debug()) part1_time.info();

    // make the second part of the commutator
    CCTimer part2_time(world, "fK" + x.name() + y.name() + ">");
    real_function_6d fKxy = apply_fK(x, y, Gscreen);
    if (parameters.debug()) part2_time.info();

    Kfxy.print_size("Kf" + x.name() + y.name());
    fKxy.print_size("fK" + x.name() + y.name());
    real_function_6d result = (Kfxy - fKxy);
    result.set_thresh(parameters.thresh_6D());
    result.print_size("[K,f]" + x.name() + y.name());
    result.truncate().reduce_rank();
    result.print_size("[K,f]" + x.name() + y.name());
    return result;
}

double
CCPotentials::make_xy_gf_ab(const CCFunction<double,3>& x, const CCFunction<double,3>& y, const CCFunction<double,3>& a, const CCFunction<double,3>& b) const {
    const real_function_3d xa = (x.function * a.function).truncate();
    const real_function_3d x_gf_a = apply_gf(xa);
    const double result = y.function.inner(x_gf_a * b.function);
    return result;
}

madness::real_function_3d
CCPotentials::apply_gf(const real_function_3d& f) const {
    std::shared_ptr<real_convolution_3d> fBSH = std::shared_ptr<real_convolution_3d>(
            BSHOperatorPtr3D(world, parameters.gamma(), parameters.lo(), parameters.thresh_poisson()));
    double bsh_prefactor = 4.0 * constants::pi;
    double prefactor = 1.0 / (2.0 * parameters.gamma());
    return prefactor * ((*g12)(f) - bsh_prefactor * (*fBSH)(f)).truncate();
}

double
CCPotentials::make_xy_u(const CCFunction<double,3>& x, const CCFunction<double,3>& y, const std::vector<CCPairFunction<double,6>>& u) const {
    double result = 0.0;
    for (size_t mm = 0; mm < u.size(); mm++) {
        result += u[mm].make_xy_u(x, y);
    }
    return result;
}

double
CCPotentials::make_xy_op_u(const CCFunction<double,3>& x, const CCFunction<double,3>& y, const CCConvolutionOperator<double,3>& op,
                           const CCPairFunction<double,6>& u) const {
    double result = 0.0;
    if (u.component->is_pure()) {
        real_function_6d xy_op = CompositeFactory<double, 6, 3>(world).particle1(copy(x.function)).particle2(
                copy(y.function)).g12(op.get_kernel());
        result = inner(u.get_function(), xy_op);
    } else if (u.component->is_decomposed()) {
        if (u.component->has_operator()) {
            if (op.type() == OpType::OT_G12 and u.decomposed().get_operator_ptr()->type() == OpType::OT_F12)
                result = make_xy_gf_ab(x, y, u.decomposed().get_a()[0], u.decomposed().get_b()[0]);
            else if (op.type() == OpType::OT_F12 and u.decomposed().get_operator_ptr()->type() == OpType::OT_G12)
                result = make_xy_gf_ab(x, y, u.decomposed().get_a()[0], u.decomposed().get_b()[0]);
            else if (op.type() == OpType::OT_F12 and u.decomposed().get_operator_ptr()->type() == OpType::OT_F12)
                result = make_xy_ff_ab(x, y, u.decomposed().get_a()[0], u.decomposed().get_b()[0]);
            else MADNESS_EXCEPTION(("xy_" + op.name() + u.name() + " not implemented").c_str(), 1);
        } else {
            for (size_t i = 0; i < u.decomposed().get_a().size(); i++)
                result += (x.function * u.decomposed().get_a()[i]).inner(op(y, u.decomposed().get_b()[i]));
        }
    } else error("Unknown CCPairFunction type in make_xy_op_u");

    return result;
}

double
CCPotentials::make_xy_op_u(const CCFunction<double,3>& x, const CCFunction<double,3>& y, const CCConvolutionOperator<double,3>& op,
                           const std::vector<CCPairFunction<double,6>>& u) const {
    double result = 0.0;
    for (size_t mm = 0; mm < u.size(); mm++) {
        const double tmp = make_xy_op_u(x, y, op, u[mm]);
        result += tmp;
    }
    return result;
}

double
CCPotentials::make_xy_op_ab(const CCFunction<double,3>& x, const CCFunction<double,3>& y, const CCConvolutionOperator<double,3>& op,
                            const CCFunction<double,3>& a, const CCFunction<double,3>& b) const {
    double result = 0.0;
    if (x.type == HOLE) {
        real_function_3d xopa = op(x, a);
        result = y.function.inner(xopa * b.function);
    } else {
        real_function_3d yopb = op(y, b);
        result = x.function.inner(yopb * a.function);
    }
    return result;
}

std::vector<CCPairFunction<double,6>>
CCPotentials::get_pair_function(const Pairs<CCPair>& pairs, const size_t i, const size_t j) const {
    if (i > j) {
        return swap_particles(pairs(j, i).functions);
    } else {
        return pairs(i, j).functions;
    }
}

madness::real_function_3d
CCPotentials::apply_s2b_operation(const CCFunction<double,3>& bra, const CCPairFunction<double,6>& u, const size_t particle) const {
    real_function_3d result;
    MADNESS_ASSERT(particle == 1 || particle == 2);
    if (u.is_pure()) {
        result = u.dirac_convolution(bra, *g12, particle);
    } else if (u.is_decomposed_no_op()) {
        result = u.dirac_convolution(bra, *g12, particle);
    } else if (u.is_op_decomposed()) {
        // retunrns <x|g12f12|x(1)y(2)>_particle
        CCFunction<double,3> a;
        CCFunction<double,3> b;
        if (particle == 1) {
            a = u.get_a()[0];
            b = u.get_b()[0];
        } else {
            a = u.get_b()[0];
            b = u.get_a()[0];
        }
        const real_function_3d tmp = (bra.function * a.function).truncate();
        const real_function_3d tmp2 = apply_gf(tmp);
        real_function_3d tmp3 = tmp2 * b.function;
        tmp3.truncate();
        result = tmp3;
    } else MADNESS_EXCEPTION("apply_s2b_operation: unknown type", 1)

    ;
    return result;
}

double
CCPotentials::overlap(const CCPair& x) const {
    if (world.rank() == 0 && parameters.debug()) std::cout << "Norms of " << x.name() << "\n";

    const size_t size = x.functions.size();
    double result = 0.0;
    for (size_t i = 0; i < size; i++) {
        for (size_t j = i; j < size; j++) {
            double factor = 1.0;
            if (i != j) factor = 2.0;     // count off diagonal elements twice since we reach only the upper triangle

            double tmp = overlap(x.functions[i], x.functions[j]);
            result += factor * tmp;
            if (world.rank() == 0 && parameters.debug())
                std::cout << std::fixed << std::setprecision(4) << "<" << x.functions[i].name() << "|"
                          << x.functions[j].name() << "> =" << tmp << "\n";
        }
    }
    return result;
}

madness::vector_real_function_3d
CCPotentials::apply_projector(const CC_vecfunction& f, const CC_vecfunction& ket_) const {
    // construct the bra state from mos
    vector_real_function_3d ket = copy(world, ket_.get_vecfunction());
    vector_real_function_3d bra = copy(world, get_mo_bra(ket_));
    MADNESS_ASSERT(ket.size() == bra.size());
    // built projector
    Projector<double, 3> O(bra, ket);
    // apply projector
    vector_real_function_3d Of = O(f.get_vecfunction());
    return Of;
}

madness::real_function_6d
CCPotentials::apply_Q12t(const real_function_6d& f, const CC_vecfunction& t) const {
    MADNESS_ASSERT(t.type == HOLE || t.type == MIXED);
    MADNESS_ASSERT(t.size() == mo_bra_.size());
    StrongOrthogonalityProjector<double, 3> Q(world);
    Q.set_spaces(mo_bra_.get_vecfunction(), t.get_vecfunction(), mo_bra_.get_vecfunction(), t.get_vecfunction());
    return Q(f);
}

madness::CCPairFunction<double,6>
CCPotentials::apply_Qt(const CCPairFunction<double,6>& f, const CC_vecfunction& t, const size_t particle, const double c) const {
    MADNESS_ASSERT(particle == 1 || particle == 2);
    MADNESS_ASSERT(f.is_decomposed_no_op());     // pure type is not needed and op_deomposed type can not be because the result would be (1-Ot)f12|xy> = f12|xy> - \sum_a|a1a2> a subtraction with different types
    if (particle == 1) {
        CCPairFunction<double,6> result(apply_Qt(f.get_a(), t, c), f.get_b());
        return result;
    } else {
        CCPairFunction<double,6> result(f.get_a(), apply_Qt(f.get_b(), t, c));
        return result;
    }
}

madness::vector_real_function_3d
CCPotentials::apply_Qt(const CC_vecfunction& f, const CC_vecfunction& ket_, const double c) const {
    vector_real_function_3d ket = ket_.get_vecfunction();
    vector_real_function_3d bra = get_mo_bra(ket_);
    Projector<double, 3> O(bra, ket);
    vector_real_function_3d Of = O(f.get_vecfunction());
    vector_real_function_3d scaled_of = c * Of;
    const vector_real_function_3d result = sub(world, f.get_vecfunction(), scaled_of);
    return result;
}

madness::CCPairFunction<double,6>
CCPotentials::apply_Ot(const CCPairFunction<double,6>& f, const CC_vecfunction& t, const size_t particle) const {
    MADNESS_ASSERT(particle == 1 || particle == 2);
    // get the right bra
    CC_vecfunction mbra;
    if (t.size() == mo_bra_.size()) mbra = CC_vecfunction(copy(world, mo_bra_.get_vecfunction()), HOLE);
    else mbra = CC_vecfunction(copy(world, get_active_mo_bra()), HOLE, parameters.freeze());

    MADNESS_ASSERT(mbra.size() == t.size());
    if (f.is_pure()) {
        vector_real_function_3d projected;
        for (const auto& ktmp : t.functions) {
            const CCFunction<double,3>& bra = mbra(ktmp.first);
            const real_function_3d kf = f.project_out(bra, particle);
            projected.push_back(kf);
        }
        if (particle == 1) {
            return CCPairFunction<double,6>(copy(world, t.get_vecfunction()), projected);
        } else {
            return CCPairFunction<double,6>(projected, copy(world, t.get_vecfunction()));
        }

    } else if (f.is_decomposed_no_op()) {
        if (particle == 1) return CCPairFunction<double,6>(apply_projector(f.get_a(), t), f.get_b());
        else return CCPairFunction<double,6>(f.get_a(), apply_projector(f.get_b(), t));
    } else if (f.is_op_decomposed()) {
        if (particle == 1) {
            const vector_real_function_3d a = copy(world, t.get_vecfunction());
//            const vector_real_function_3d b = mul(world, f.get_b()[0], f.op->operator()(mbra, f.get_a()[0]));
            const vector_real_function_3d b = mul(world, f.get_b()[0], f.decomposed().get_operator_ptr()->operator()(mbra, f.get_a()[0]));
            return CCPairFunction<double,6>(a, b);
        } else {
            const vector_real_function_3d a = mul(world, f.get_a()[0], f.decomposed().get_operator_ptr()->operator()(mbra, f.get_b()[0]));
            const vector_real_function_3d b = copy(world, t.get_vecfunction());
            return CCPairFunction<double,6>(a, b);
        }
    } else MADNESS_EXCEPTION("Should not end up here", 1)

    ;
    return CCPairFunction<double,6>(vector_real_function_3d(), vector_real_function_3d());
}

madness::real_function_6d
CCPotentials::apply_G(const CCPairFunction<double,6>& u, const real_convolution_6d& G) const {
    CCTimer time(world, "Applying G on " + u.name());
    return apply(G,u).get_function();
//    real_function_6d result = real_function_6d(world);
//    if (u.is_pure()) {
//        result = G(u.pure().get_function());
//    } else if (u.is_decomposed_no_op()) {
//        MADNESS_ASSERT(u.get_a().size() == u.get_b().size());
//        if (u.get_a().size() == 0) output.warning("!!!!!!!in G(ab): a.size()==0 !!!!!!");
//
//        for (size_t k = 0; k < u.get_a().size(); k++) {
//            const real_function_6d tmp = G(u.get_a()[k], u.get_b()[k]);
//            result += tmp;
//        }
//    } else error("Apply_G to CCPairFunction<double,6> of type other than pure or decomposed");
//
//    time.info(true, result.norm2());
//    if (result.norm2() == 0.0) output.warning("Gab is Zero");
//
//    return result;
}

madness::vector_real_function_3d
CCPotentials::get_CC2_singles_potential_gs(const CC_vecfunction& singles, const Pairs<CCPair>& doubles) const {
    CCTimer time(world, "CC2 Singles potential");
    vector_real_function_3d fock_residue = potential_singles_gs(singles, doubles, POT_F3D_);
    // CC2 Singles potential: Q(S4c) + Qt(ccs+s2b+s2c)
    vector_real_function_3d Vccs = potential_singles_gs(singles, doubles, POT_ccs_);
    vector_real_function_3d Vs2b = potential_singles_gs(singles, doubles, POT_s2b_);
    vector_real_function_3d Vs2c = potential_singles_gs(singles, doubles, POT_s2c_);
    vector_real_function_3d Vs4b = potential_singles_gs(singles, doubles, POT_s4b_);
    vector_real_function_3d Vs4c = potential_singles_gs(singles, doubles, POT_s4c_);
    vector_real_function_3d Vs4a = apply_projector(Vs2b, singles);     // need to subtract
    vector_real_function_3d unprojected = add(world, Vccs, add(world, Vs2b, add(world, Vs2c, add(world, Vs4b,
                                                                                                 sub(world, Vs4c,
                                                                                                     Vs4a)))));
    vector_real_function_3d potential = apply_Qt(unprojected, mo_ket_);
    truncate(world, potential);
    get_potentials.insert(copy(world, potential), singles, POT_singles_);
    time.info(true, norm2(world, potential));
    const vector_real_function_3d result = add(world, potential, fock_residue);
    return result;
}

madness::vector_real_function_3d
CCPotentials::get_CCS_potential_ex(CC_vecfunction& x, const bool print) const {
    if (x.type != RESPONSE) error("get_CCS_response_potential: Wrong type of input singles");

    Pairs<CCPair> empty_doubles;
    CC_vecfunction empty_singles(PARTICLE);
    const vector_real_function_3d fock_residue = potential_singles_ex(empty_singles, empty_doubles, x, empty_doubles,
                                                                      POT_F3D_);
    vector_real_function_3d potential = potential_singles_ex(empty_singles, empty_doubles, x, empty_doubles, POT_cis_);
    // the fock residue does not get projected, but all the rest
    potential = apply_Qt(potential, mo_ket_);
    truncate(world, potential);
    get_potentials.insert(copy(world, potential), x, POT_singles_);
    vector_real_function_3d result = add(world, fock_residue, potential);
    truncate(world, result);
    const double omega = compute_cis_expectation_value(x, result, print);
    x.omega = omega;
    return result;
}

madness::vector_real_function_3d
CCPotentials::get_CC2_singles_potential_ex(const CC_vecfunction& gs_singles, const Pairs<CCPair>& gs_doubles,
                                           CC_vecfunction& ex_singles, const Pairs<CCPair>& response_doubles) const {
    MADNESS_ASSERT(gs_singles.type == PARTICLE);
    MADNESS_ASSERT(ex_singles.type == RESPONSE);
    const vector_real_function_3d fock_residue = potential_singles_ex(gs_singles, gs_doubles, ex_singles,
                                                                      response_doubles, POT_F3D_);
    vector_real_function_3d Vccs = potential_singles_ex(gs_singles, gs_doubles, ex_singles, response_doubles, POT_ccs_);
    vector_real_function_3d Vs2b = potential_singles_ex(gs_singles, gs_doubles, ex_singles, response_doubles, POT_s2b_);
    vector_real_function_3d Vs2c = potential_singles_ex(gs_singles, gs_doubles, ex_singles, response_doubles, POT_s2c_);
    vector_real_function_3d Vs4b = potential_singles_ex(gs_singles, gs_doubles, ex_singles, response_doubles, POT_s4b_);
    vector_real_function_3d Vs4c = potential_singles_ex(gs_singles, gs_doubles, ex_singles, response_doubles, POT_s4c_);
    // make low scaling s4a potential
    // -Otau(s2b_response) + -Ox(s2b_gs)
    // maybe store full s2b potential of gs
    // both need to be subtracted
    vector_real_function_3d s2b_gs = potential_singles_gs(gs_singles, gs_doubles, POT_s2b_);
    vector_real_function_3d Vs4a =
            -1.0 * add(world, apply_projector(s2b_gs, ex_singles), apply_projector(Vs2b, gs_singles));
    //add up
    vector_real_function_3d unprojected = add(world, Vccs, add(world, Vs2b, add(world, Vs2c, add(world, Vs4a,
                                                                                                 add(world, Vs4b,
                                                                                                     Vs4c)))));
    vector_real_function_3d potential = apply_Qt(unprojected, mo_ket_);
    if (parameters.debug()) {
        // debug
        vector_real_function_3d xbra = mul(world, nemo_->ncf->square(), ex_singles.get_vecfunction());
        const double ccs = inner(world, xbra, Vccs).sum();
        const double s2b = inner(world, xbra, Vs2b).sum();
        const double s2c = inner(world, xbra, Vs2c).sum();
        const double s4a = inner(world, xbra, Vs4a).sum();
        const double s4b = inner(world, xbra, Vs4b).sum();
        const double s4c = inner(world, xbra, Vs4c).sum();
        std::cout << std::fixed << std::setprecision(10) << "functional response energies:" << "\n<x|ccs>=" << ccs
                  << "\n<x|S2b>=" << s2b << "\n<x|S2c>=" << s2c << "\n<x|s4a>=" << s4a << "\n<x|s4b>="
                  << s4b << "\n<x|s4c>=" << s4c << "\n";
        // debug end
    }
    // storing potential
    get_potentials.insert(copy(world, potential), ex_singles, POT_singles_);
    vector_real_function_3d result = add(world, fock_residue, potential);
    truncate(world, result);
    const double omega = compute_cis_expectation_value(ex_singles, result);
    ex_singles.omega = omega;
    return result;
}

madness::vector_real_function_3d
CCPotentials::get_ADC2_singles_potential(const Pairs<CCPair>& gs_doubles, CC_vecfunction& ex_singles,
                                         const Pairs<CCPair>& response_doubles) const {
    MADNESS_ASSERT(ex_singles.type == RESPONSE);
    vector_real_function_3d zero = zero_functions<double, 3>(world, get_active_mo_ket().size());
    CC_vecfunction tau(zero, PARTICLE, parameters.freeze());
    const vector_real_function_3d result = get_CC2_singles_potential_ex(tau, gs_doubles, ex_singles, response_doubles);
    return result;
}

double
CCPotentials::potential_energy_gs(const CC_vecfunction& bra, const CC_vecfunction& singles,
                                  const Pairs<CCPair>& doubles, const PotentialType& name) const {
    // sanity check
    MADNESS_ASSERT(singles.type == PARTICLE);
    CCTimer timer(world, "potential energy of " + assign_name(name));
    double result = 0.0;
    if (name == POT_s3a_) {
        result = x_s3a(bra, singles);
    } else if (name == POT_s3b_) {
        result = x_s3b(bra, singles);
    } else if (name == POT_s3c_) {
        result = x_s3c(bra, singles);
    } else if (name == POT_s5b_) {
        result = x_s5b(bra, singles, singles);
    } else if (name == POT_s5c_) {
        result = x_s5c(bra, singles, singles);
    } else if (name == POT_s6_) {
        result = x_s6(bra, singles, singles, singles);
    } else if (name == POT_F3D_) {
        result = x_s3a(bra, singles) - compute_kinetic_energy(bra.get_vecfunction(), singles.get_vecfunction());
    } else if (name == POT_ccs_) {
        result = x_s3c(bra, singles) + x_s5b(bra, singles, singles) + x_s5c(bra, singles, singles) +
                 x_s6(bra, singles, singles, singles);
    } else if (name == POT_s2b_) {
        result = x_s2b(bra, doubles);
    } else if (name == POT_s2c_) {
        result = x_s2c(bra, doubles);
    } else if (name == POT_s4a_) {
        result = x_s4a(bra, singles, doubles);
    } else if (name == POT_s4b_) {
        result = x_s4b(bra, singles, doubles);
    } else if (name == POT_s4c_) {
        result = x_s4c(bra, singles, doubles);
    }

    const std::pair<double, double> time = timer.current_time();
    if (result == 0.0) output.warning("Result of <x" + assign_name(name) + "> is zero!");

    if (world.rank() == 0) {
        std::cout << std::fixed << std::setprecision(10) << "<x|" << assign_name(name) << ">=" << result << ", "
                  << time.first << " (wall), " << time.second << " (cpu)" << "\n";
    }
    return result;
}

madness::vector_real_function_3d
CCPotentials::potential_singles_gs(const CC_vecfunction& singles, const Pairs<CCPair>& doubles,
                                   const PotentialType& name) const {
    MADNESS_ASSERT(singles.type == PARTICLE);
    vector_real_function_3d result;
    CCTimer timer(world, "Singles-Potential:" + assign_name(name));
    if (name == POT_F3D_) {
        result = fock_residue_closed_shell(singles);
    } else if (name == POT_ccs_) {
        const CC_vecfunction t = make_t_intermediate(singles);
        result = apply_Qt(ccs_unprojected(t, singles),
                          t);     // this is not the full t projector, but the potential will be projeted afterwards and this will unclude th frozen mos
    } else if (name == POT_s2b_) {
        //	// calculate the s2b potential and afterwards the s4a potential from the s2b potential
        //	// because:  Qt(S2b) = S2b + S4a
        //	vecfuncT result_s2b= s2b(singles,doubles);
        //	// some information
        //	if(world.rank()==0) std::cout <<"||" << assign_name(pot_s2b_) << "||=" << norm2(world,result_s2b) << ", "
        //	    << get_size(world,result_s2b) << " (GB), " << timer.current_time().first << "s (wall), " << timer.current_time().second << "s (cpu)\n";
        //	// get s4a from s2b
        //	vecfuncT result_s4a= s4a_from_s2b(result_s2b,singles);
        //	// some information
        //	if(world.rank()==0) std::cout <<"||" << assign_name(pot_s4a_) << "||=" << norm2(world,result_s4a) << ", "
        //	    << get_size(world,result_s4a) << " (GB), " << timer.current_time().first << "s (wall), " << timer.current_time().second << "s (cpu)\n";
        //	result = add(world,result_s2b,result_s4a);
        // returns the s2b potential (unprojected)
        result = s2b(singles, doubles);
    } else if (name == POT_s2c_) {
        result = s2c(singles, doubles);
    } else if (name == POT_s4a_) {
        error("potential_singles: Demanded s4a potential -> this is calculated along with the s2b potential");
    } else if (name == POT_s4b_) {
        result = s4b(singles, doubles);
    } else if (name == POT_s4c_) {
        result = s4c(singles, doubles);
    } else MADNESS_EXCEPTION(("potential_singles: Unknown potential " + assign_name(name)).c_str(), 1)

    ;
    const double size = get_size(world, result);
    const double norm = norm2(world, result);
    const std::pair<double, double> time = timer.current_time();
    if (world.rank() == 0) {
        std::cout << "||" << assign_name(name) << "||=" << std::fixed << std::setprecision(6) << norm << std::scientific
                  << std::setprecision(1) << ", " << size << " (GB), " << time.first
                  << "s (wall), " << time.second << "s (cpu)\n";
    }
    truncate(world, result);
    return result;
}

double
CCPotentials::potential_energy_ex(const CC_vecfunction& bra, const CC_vecfunction& singles_gs,
                                  const Pairs<CCPair>& doubles_gs, const CC_vecfunction& singles_ex,
                                  const Pairs<CCPair>& doubles_ex,
                                  const PotentialType& name) const {
    // sanity check
    MADNESS_ASSERT(singles_gs.type == PARTICLE);
    MADNESS_ASSERT(singles_ex.type == RESPONSE);
    CCTimer timer(world, "potential energy of " + assign_name(name));
    double result = 0.0;
    if (name == POT_s3a_) {
        result = x_s3a(bra, singles_ex);
    } else if (name == POT_s3b_) {
        result = x_s3b(bra, singles_ex);
    } else if (name == POT_s3c_) {
        result = x_s3c(bra, singles_ex);
    } else if (name == POT_s5b_) {
        result = x_s5b(bra, singles_ex, singles_gs) + x_s5b(bra, singles_gs, singles_ex);
    } else if (name == POT_s5c_) {
        result = x_s5c(bra, singles_ex, singles_gs) + x_s5c(bra, singles_ex, singles_gs);
    } else if (name == POT_s6_) {
        result = x_s6(bra, singles_ex, singles_gs, singles_gs) + x_s6(bra, singles_gs, singles_ex, singles_gs) +
                 x_s6(bra, singles_gs, singles_gs, singles_ex);
    } else if (name == POT_F3D_) {
        result = x_s3a(bra, singles_ex) - compute_kinetic_energy(bra.get_vecfunction(), singles_ex.get_vecfunction());
    } else if (name == POT_ccs_) {
        result = x_s3c(bra, singles_ex) + x_s5b(bra, singles_ex, singles_gs) + x_s5c(bra, singles_ex, singles_gs) +
                 x_s6(bra, singles_ex, singles_gs, singles_gs) + x_s5b(bra, singles_gs, singles_ex)
                 + x_s5c(bra, singles_gs, singles_ex) + x_s6(bra, singles_gs, singles_ex, singles_gs) +
                 x_s6(bra, singles_gs, singles_gs, singles_ex);
    } else if (name == POT_s2b_) {
        result = x_s2b(bra, doubles_ex);
    } else if (name == POT_s2c_) {
        result = x_s2c(bra, doubles_ex);
    } else if (name == POT_s4a_) {
        result = x_s4a(bra, singles_gs, doubles_ex) + x_s4a(bra, singles_ex, doubles_gs);
    } else if (name == POT_s4b_) {
        result = x_s4b(bra, singles_gs, doubles_ex) + x_s4b(bra, singles_ex, doubles_gs);
    } else if (name == POT_s4c_) {
        result = x_s4c(bra, singles_gs, doubles_ex) + x_s4c(bra, singles_ex, doubles_gs);
    }

    const std::pair<double, double> time = timer.current_time();
    if (result == 0.0) output.warning("Result of <x" + assign_name(name) + "> is zero!");

    if (world.rank() == 0) {
        std::cout << std::fixed << std::setprecision(10) << "<x|" << assign_name(name) << ">=" << result << ", "
                  << time.first << " (wall), " << time.second << " (cpu)" << "\n";
    }
    return result;
}

madness::vector_real_function_3d
CCPotentials::potential_singles_ex(const CC_vecfunction& singles_gs, const Pairs<CCPair>& doubles_gs,
                                   const CC_vecfunction& singles_ex, const Pairs<CCPair>& doubles_ex,
                                   const PotentialType& name) const {
    //if(mo_ket_.size()>1) output.warning("Potential for ExSingles is not ready for more than one orbital");
    // sanity check
    MADNESS_ASSERT(singles_gs.type == PARTICLE);
    MADNESS_ASSERT(singles_ex.type == RESPONSE);
    vector_real_function_3d result;
    CCTimer timer(world, "timer-ex-potential");
    if (name == POT_F3D_) {
        result = fock_residue_closed_shell(singles_ex);
    } else if (name == POT_ccs_) {
        const CC_vecfunction t = make_t_intermediate(singles_gs);
        vector_real_function_3d part1 = apply_Qt(ccs_unprojected(t, singles_ex), t);
        vector_real_function_3d part2 = apply_Qt(ccs_unprojected(singles_ex, singles_gs), t);
        vector_real_function_3d part3 = apply_projector(ccs_unprojected(t, singles_gs), singles_ex);
        vector_real_function_3d tmp = add(world, part1, part2);
        result = sub(world, tmp, part3);
    } else if (name == POT_s2b_) {
        result = s2b(singles_ex, doubles_ex);
    } else if (name == POT_s2c_) {
        result = s2c(singles_ex, doubles_ex);
    } else if (name == POT_s4a_) {
        error("potential_singles: Demanded s4a potential -> this is calculated from the s2b potential");
    } else if (name == POT_s4b_) {
        vector_real_function_3d s4b_part1 = s4b(singles_gs, doubles_ex);
        vector_real_function_3d s4b_part2 = s4b(singles_ex, doubles_gs);
        result = add(world, s4b_part1, s4b_part2);
    } else if (name == POT_s4c_) {
        vector_real_function_3d s4c_part1 = s4c(singles_gs, doubles_ex);
        vector_real_function_3d s4c_part2 = s4c(singles_ex, doubles_gs);
        result = add(world, s4c_part1, s4c_part2);
    } else if (name == POT_cis_) {
        result = ccs_unprojected(CC_vecfunction(get_active_mo_ket(), HOLE, parameters.freeze()), singles_ex);
    } else MADNESS_EXCEPTION(("potential_singles: Unknown potential " + assign_name(name)).c_str(), 1)

    ;
    const double size = get_size(world, result);
    const double norm = norm2(world, result);
    const std::pair<double, double> time = timer.current_time();
    if (world.rank() == 0) {
        std::cout << "||" << assign_name(name) << "||=" << std::fixed << std::setprecision(6) << norm << ", "
                  << std::scientific << std::setprecision(1) << size << " (GB), " << time.first
                  << "s (wall), " << time.second << "s (cpu)\n";
    }
    if (result.empty()) MADNESS_EXCEPTION("Result is empty", 1);
    truncate(world, result);
    return result;
}

madness::vector_real_function_3d
CCPotentials::fock_residue_closed_shell(const CC_vecfunction& singles) const {
    //	vecfuncT tau = singles.get_vecfunction();
    CCTimer timer_J(world, "J");
    //	vecfuncT J = mul(world, intermediates_.get_hartree_potential(), tau);
    vector_real_function_3d J;
    for (const auto& tmpi : singles.functions) {
        const CCFunction<double,3>& taui = tmpi.second;
        real_function_3d hartree_potential = real_function_3d(world);
        for (const auto& tmpk : mo_ket_.functions)
            hartree_potential += (*g12)(mo_bra_(tmpk.first), tmpk.second);
        const real_function_3d Ji = hartree_potential * taui.function;
        J.push_back(Ji);
    }
    truncate(world, J);
    scale(world, J, 2.0);
    timer_J.info(true, norm2(world, J));
    CCTimer timer_K(world, "K");
    vector_real_function_3d vK;
    for (const auto& tmpi : singles.functions) {
        const CCFunction<double,3>& taui = tmpi.second;
        const real_function_3d Ki = K(taui);
        vK.push_back(Ki);
    }
    scale(world, vK, -1.0);
    timer_K.info(true, norm2(world, vK));
    // apply nuclear potential
    Nuclear<double, 3> Uop(world, nemo_.get());
    vector_real_function_3d Upot = Uop(singles.get_vecfunction());
    vector_real_function_3d KU = add(world, vK, Upot);
    return add(world, J, KU);
}

madness::real_function_6d
CCPotentials::K(const real_function_6d& u, const bool symmetric) const {
    real_function_6d result = real_factory_6d(world).compressed();
    // K(1) Part
    result += apply_K(u, 1);
    // K(2) Part
    if (symmetric) {
        result += swap_particles(result);
    } else result += apply_K(u, 2);

    result.print_size("K|u>");
    return (result.truncate(parameters.tight_thresh_6D()));
}

madness::real_function_6d
CCPotentials::K_macrotask(World& world, const std::vector<real_function_3d>& mo_ket,
                          const std::vector<real_function_3d>& mo_bra, const real_function_6d& u,
                          const bool symmetric, const CCParameters& parameters) {
    real_function_6d result = real_factory_6d(world).compressed();
    // K(1) Part
    result += apply_K_macrotask(world, mo_ket, mo_bra, u, 1, parameters);
    // K(2) Part
    if (symmetric) {
        result += madness::swap_particles(result);
    } else result += apply_K_macrotask(world, mo_ket, mo_bra, u, 2, parameters);

    if (parameters.debug()) result.print_size("K|u>");
    return (result.truncate(parameters.tight_thresh_6D()));
}

madness::real_function_3d
CCPotentials::K(const CCFunction<double,3>& f) const {
    real_function_3d result = real_factory_3d(world);
    for (const auto& k_iterator : mo_ket_.functions) {
        result += (*g12)(mo_bra_(k_iterator.first), f) * mo_ket_(k_iterator.first).function;
    }
    return result;
}

madness::real_function_3d
CCPotentials::K_macrotask(World& world, const std::vector<real_function_3d>& mo_ket,
                          const std::vector<real_function_3d>& mo_bra, const real_function_3d& f,
                          const CCParameters& parameters) {
    real_function_3d result = real_factory_3d(world);
    real_convolution_3d g12 = CoulombOperator(world, parameters.lo(), parameters.thresh_poisson());
    for (size_t k = 0; k < mo_ket.size(); k++) {
        result += ((g12)(mo_bra[k] * f).truncate()) * mo_ket[k];
    }
    return result;
}

madness::real_function_6d
CCPotentials::apply_K(const real_function_6d& u, const size_t& particle) const {
    MADNESS_ASSERT(particle == 1 || particle == 2);
    //poisson->particle()=particle;
    real_function_6d result = real_factory_6d(world).compressed();
    for (size_t k = 0; k < mo_ket_.size(); k++) {
        real_function_6d copyu = copy(u);
        real_function_6d X = (multiply(copyu, copy(mo_bra_(k).function), particle)).truncate();
        //      real_function_6d Y=(*poisson)(X);
        real_function_6d Y = (*g12)(X, particle);     // overwrite X to save space
        result += (multiply(copy(Y), copy(mo_ket_(k).function),
                            particle)).truncate();     // this will destroy X, but I d not intend to use it again so I choose here to save this copy
    }
    return result;
}

madness::real_function_6d
CCPotentials::apply_K_macrotask(World& world, const std::vector<real_function_3d>& mo_ket,
                                const std::vector<real_function_3d>& mo_bra,
                                const real_function_6d& u, const size_t& particle, const CCParameters& parameters) {
    MADNESS_ASSERT(particle == 1 || particle == 2);
    //poisson->particle()=particle;
    real_function_6d result = real_factory_6d(world).compressed();
    //const double lo = 1.e-6;
    //const double bsh_eps = 1.e-7;
    real_convolution_3d g12 = CoulombOperator(world, parameters.lo(), parameters.thresh_poisson());
    g12.particle() = particle;
    for (size_t k = 0; k < mo_ket.size(); k++) {
        real_function_6d copyu = copy(u);
        real_function_6d X = (multiply(copyu, copy(mo_bra[k]), particle)).truncate();
        //      real_function_6d Y=(*poisson)(X);
        real_function_6d Y = g12(X);     // overwrite X to save space
        result += (multiply(copy(Y), copy(mo_ket[k]),
                            particle)).truncate();     // this will destroy X, but I d not intend to use it again so I choose here to save this copy
    }
    return result.truncate(parameters.tight_thresh_3D()*3.0).reduce_rank(parameters.tight_thresh_6D()*3.0);
}

madness::real_function_6d
CCPotentials::apply_Kf(const CCFunction<double,3>& x, const CCFunction<double,3>& y) const {
    bool symmetric = false;
    if ((x.type == y.type) && (x.i == y.i)) symmetric = true;

    // First make the 6D function f12|x,y>
    real_function_6d f12xy = make_f_xy(x, y);
    f12xy.truncate().reduce_rank();
    // Apply the Exchange Operator
    real_function_6d result = K(f12xy, symmetric);
    return result;
}

madness::real_function_6d
CCPotentials::apply_fK(const CCFunction<double,3>& x, const CCFunction<double,3>& y, const real_convolution_6d *Gscreen) const {
    const bool symmetric = (x.type == y.type && x.i == y.i);
    const real_function_3d Kx = K(x);
    const real_function_6d fKphi0b = make_f_xy(CCFunction<double,3>(Kx, x.i, UNDEFINED), y, Gscreen);
    real_function_6d fKphi0a;
    if (symmetric) fKphi0a = swap_particles(fKphi0b);
    else {
        real_function_3d Ky = K(y);
        fKphi0a = make_f_xy(x, CCFunction<double,3>(Ky, y.i, UNDEFINED), Gscreen);
    }
    const real_function_6d fKphi0 = (fKphi0a + fKphi0b);
    return fKphi0;
}

madness::real_function_6d
CCPotentials::make_f_xy(const CCFunction<double,3>& x, const CCFunction<double,3>& y, const real_convolution_6d *Gscreen) const {
    std::string screen = "";
    if (Gscreen != NULL) screen = "screened ";

    CCTimer timer(world, "Making " + screen + "f|" + x.name() + "," + y.name() + "> ");
    real_function_6d fxy = CompositeFactory<double, 6, 3>(world).g12(corrfac.f()).particle1(copy(x.function)).particle2(
            copy(y.function));
    if (Gscreen == NULL) fxy.fill_tree().truncate().reduce_rank();
    else fxy.fill_cuspy_tree(*Gscreen).truncate().reduce_rank();

    timer.info(parameters.debug(), fxy.norm2());
    if (x.type != UNDEFINED && y.type != UNDEFINED) {
        CCTimer timer_db(world, "f|xy> sanity check");
        const double test1 = (mo_bra_(y.i).function).inner(fxy.project_out(mo_bra_(x.i).function, 0));
        const double test2 = (mo_bra_(y.i).function).inner((*f12)(mo_bra_(x.i), x) * y.function);
        const double sanity = test1 - test2;
        if (fabs(sanity) > FunctionDefaults<6>::get_thresh()) {
            if (world.rank() == 0)
                std::cout << std::fixed << std::setprecision(6) << "test1=" << test1 << "\ntest2=" << test2 << "\ndiff="
                          << sanity << "\n";

            output.warning("make f|xy> not accurate!");
        }
        timer_db.info(parameters.debug(), sanity);
    }
    return fxy;
}

madness::real_function_6d
CCPotentials::make_f_xy_macrotask(World& world, const real_function_3d& x_ket, const real_function_3d& y_ket,
                                  const real_function_3d& x_bra, const real_function_3d& y_bra,
                                  const size_t& i, const size_t& j, const CCParameters& parameters,
                                  const FuncType& x_type, const FuncType& y_type,
                                  const real_convolution_6d *Gscreen){
    CorrelationFactor corrfac(world, parameters.gamma(), 1.e-7, parameters.lo());
    const std::string x_name = "phi" + stringify(i);
    const std::string y_name = "phi" + stringify(j);

    std::string screen = "";
    if (Gscreen != NULL) screen = "screened ";

    CCTimer timer(world, "Making " + screen + "f|" + x_name + "," + y_name + "> ");
    real_function_6d fxy = CompositeFactory<double, 6, 3>(world).g12(corrfac.f()).
                                                    particle1(copy(x_ket)).particle2(copy(y_ket));
    if (Gscreen == NULL) fxy.fill_tree().truncate().reduce_rank();
    else fxy.fill_cuspy_tree(*Gscreen).truncate().reduce_rank();

    timer.info(parameters.debug(), fxy.norm2());
    if (x_type != UNDEFINED && y_type != UNDEFINED) {
        CCTimer timer_db(world, "f|xy> sanity check");
        //const double fourpi = 4.0 * constants::pi;
       // const double lo = 1.e-6;
        //const double bsh_eps = 1.e-7;
        real_convolution_3d slaterf12 = SlaterF12Operator(world, corrfac.gamma(),
                                                          parameters.lo(), parameters.thresh_poisson());
        const double test1 = y_bra.inner(fxy.project_out(x_bra, 0));
        const double test2 = y_bra.inner(slaterf12(x_bra * x_ket) * y_ket);
        const double sanity = test1 - test2;
        if (fabs(sanity) > FunctionDefaults<6>::get_thresh()) {
            if (world.rank() == 0)
                std::cout << std::fixed << std::setprecision(6) << "test1=" << test1
                          << "\ntest2=" << test2 << "\ndiff=" << sanity << "\n";
            std::cout << ("make f|xy> not accurate!\n");
        }
        timer_db.info(parameters.debug(), sanity);
    }
    return fxy;
}

madness::vector_real_function_3d
CCPotentials::ccs_unprojected(const CC_vecfunction& ti, const CC_vecfunction& tk) const {
    vector_real_function_3d result;
    for (const auto& itmp : ti.functions) {
        real_function_3d kgtk = real_factory_3d(world);
        for (const auto& ktmp : tk.functions)
            kgtk += (*g12)(mo_bra_(ktmp.first), ktmp.second);
        const real_function_3d kgtk_ti = kgtk * ti(itmp.first).function;
        real_function_3d kgti_tk = real_factory_3d(world);
        for (const auto& ktmp : tk.functions)
            kgti_tk += (*g12)(mo_bra_(ktmp.first), ti(itmp.first)) * tk(ktmp.first).function;
        const real_function_3d resulti = 2.0 * kgtk_ti - kgti_tk;
        result.push_back(resulti);
    }
    return result;
}

madness::real_function_3d
CCPotentials::make_density(const CC_vecfunction& x) const {
    real_function_3d result = real_factory_3d(world);
    for (const auto& ktmp : x.functions) {
        const size_t k = ktmp.first;
        result += 2.0 * mo_bra_(k).function * (x(k).function);
    }
    result.truncate();
    return result;
}

double
CCPotentials::x_s3a(const CC_vecfunction& x, const CC_vecfunction& t) const {
    MADNESS_ASSERT(x.size() == t.size());
    Nuclear<double, 3> Uop(world, nemo_.get());
    vector_real_function_3d Ut = Uop(t.get_vecfunction());
    const double nuc = inner(world, x.get_vecfunction(), Ut).sum();
    double pot = 0.0;
    for (const auto& itmp : x.functions) {
        const size_t i = itmp.first;
        for (const auto& ktmp : mo_ket_.functions) {
            // unfrozen summation !!!!!! important !!!!
            const size_t k = ktmp.first;
            const double gpart = make_xy_op_ab(x(i), mo_bra_(k), *g12, t(i), mo_ket_(k));
            const double xpart = make_xy_op_ab(x(i), mo_bra_(k), *g12, mo_ket_(k), t(i));
            pot += (2.0 * gpart - xpart);
        }
    }
    double kinetic = compute_kinetic_energy(x.get_vecfunction(), t.get_vecfunction());
    return kinetic + pot + nuc;
}

double
CCPotentials::x_s3b(const CC_vecfunction& x, const CC_vecfunction& t) const {
    MADNESS_ASSERT(x.size() == t.size());
    const size_t freeze = x.functions.cbegin()->first;
    Tensor<double> overlap = inner(world, x.get_vecfunction(), t.get_vecfunction());
    double result = 0.0;
    for (int i = 0; i < overlap.size(); i++) {
        result += get_orbital_energies()[i + freeze] * overlap(i);
    }
    return -1.0 * result;
}

double
CCPotentials::x_s3c(const CC_vecfunction& x, const CC_vecfunction& t) const {
    MADNESS_ASSERT(x.size() == t.size());
    double result = 0.0;
    for (const auto& itmp : x.functions) {
        const size_t i = itmp.first;
        for (const auto& ktmp : t.functions) {
            const size_t k = ktmp.first;
            result += (2.0 * make_xy_op_ab(x(i), mo_bra_(k), *g12, mo_ket_(i), t(k)) -
                       make_xy_op_ab(x(i), mo_bra_(k), *g12, t(k), mo_ket_(i)));
        }
    }
    return result;
}

double
CCPotentials::x_s5b(const CC_vecfunction& x, const CC_vecfunction& t1, const CC_vecfunction& t2) const {
    MADNESS_ASSERT(x.size() == t1.size());
    MADNESS_ASSERT(t1.size() == t2.size());
    double result = 0.0;
    for (const auto& itmp : x.functions) {
        const size_t i = itmp.first;
        for (const auto& ktmp : t1.functions) {
            const size_t k = ktmp.first;
            result += (2.0 * make_xy_op_ab(x(i), mo_bra_(k), *g12, t1(i), t2(k)) -
                       make_xy_op_ab(x(i), mo_bra_(k), *g12, t2(k), t1(i)));
        }
    }
    return result;
}

double
CCPotentials::x_s5c(const CC_vecfunction& x, const CC_vecfunction& t1, const CC_vecfunction& t2) const {
    MADNESS_ASSERT(x.size() == t1.size());
    MADNESS_ASSERT(t1.size() == t2.size());
    double result = 0.0;
    for (const auto& itmp : x.functions) {
        const size_t i = itmp.first;
        for (const auto& ktmp : t1.functions) {
            const size_t k = ktmp.first;
            for (const auto& ltmp : t2.functions) {
                const size_t l = ltmp.first;
                result += (2.0 * make_xy_op_ab(mo_bra_(l), mo_bra_(k), *g12, mo_ket_(i), t1(k)) -
                           make_xy_op_ab(mo_bra_(l), mo_bra_(k), *g12, t1(k), mo_ket_(i))) *
                          x(i).function.inner(t2(l).function);
            }
        }
    }
    return -1.0 * result;
}

double
CCPotentials::x_s6(const CC_vecfunction& x, const CC_vecfunction& t1, const CC_vecfunction& t2,
                   const CC_vecfunction& t3) const {
    MADNESS_ASSERT(x.size() == t1.size());
    MADNESS_ASSERT(t1.size() == t2.size());
    double result = 0.0;
    for (const auto& itmp : x.functions) {
        const size_t i = itmp.first;
        for (const auto& ktmp : t1.functions) {
            const size_t k = ktmp.first;
            for (const auto& ltmp : t2.functions) {
                const size_t l = ltmp.first;
                result += (2.0 * make_xy_op_ab(mo_bra_(l), mo_bra_(k), *g12, t3(i), t1(k)) -
                           make_xy_op_ab(mo_bra_(l), mo_bra_(k), *g12, t1(k), t3(i))) *
                          x(i).function.inner(t2(l).function);
            }
        }
    }
    return -1.0 * result;
}

double
CCPotentials::x_s2b(const CC_vecfunction& x, const Pairs<CCPair>& u) const {
    double result = 0.0;
    for (const auto& itmp : x.functions) {
        const size_t i = itmp.first;
        for (const auto& ktmp : x.functions) {
            const size_t k = ktmp.first;
            result += (2.0 * make_xy_op_u(x(i), mo_bra_(k), *g12, get_pair_function(u, i, k)) -
                       make_xy_op_u(mo_bra_(k), x(i), *g12, get_pair_function(u, i, k)));
        }
    }
    return result;
}

double
CCPotentials::x_s2c(const CC_vecfunction& x, const Pairs<CCPair>& u) const {
    double result = 0.0;
    for (const auto& itmp : x.functions) {
        const size_t i = itmp.first;
        for (const auto& ktmp : x.functions) {
            const size_t k = ktmp.first;
            const real_function_3d kgi = (*g12)(mo_bra_(k), mo_ket_(i));
            for (const auto& ltmp : x.functions) {
                const size_t l = ltmp.first;
                real_function_3d l_kgi = (mo_bra_(l).function * kgi).truncate();
                result += 2.0 * make_xy_u(x(i), l_kgi, get_pair_function(u, k, l)) -
                          make_xy_u(l_kgi, x(i), get_pair_function(u, k, l));
            }
        }
    }
    return -1.0 * result;
}

double
CCPotentials::x_s4a(const CC_vecfunction& x, const CC_vecfunction& t, const Pairs<CCPair>& u) const {
    double result = 0.0;
    MADNESS_ASSERT(x.size() == t.size());
    for (const auto& itmp : x.functions) {
        const size_t i = itmp.first;
        for (const auto& ktmp : x.functions) {
            const size_t k = ktmp.first;
            for (const auto& ltmp : x.functions) {
                const size_t l = ltmp.first;
                result += (2.0 * make_xy_op_u(mo_bra_(l), mo_bra_(k), *g12, get_pair_function(u, i, k)) -
                           make_xy_op_u(mo_bra_(k), mo_bra_(l), *g12, get_pair_function(u, i, k))) *
                          x(i).function.inner(t(l).function);
            }
        }
    }
    return -1.0 * result;
}

double
CCPotentials::x_s4b(const CC_vecfunction& x, const CC_vecfunction& t, const Pairs<CCPair>& u) const {
    double result = 0.0;
    MADNESS_ASSERT(x.size() == t.size());
    for (const auto& itmp : x.functions) {
        const size_t i = itmp.first;
        for (const auto& ktmp : x.functions) {
            const size_t k = ktmp.first;
            const real_function_3d kgti = (*g12)(mo_bra_(k), t(i));
            for (const auto& ltmp : x.functions) {
                const size_t l = ltmp.first;
                real_function_3d l_kgti = (mo_bra_(l).function * kgti).truncate();
                result += 2.0 * make_xy_u(x(i), l_kgti, get_pair_function(u, k, l)) -
                          make_xy_u(l_kgti, x(i), get_pair_function(u, k, l));
            }
        }
    }
    return -1.0 * result;
}

double
CCPotentials::x_s4c(const CC_vecfunction& x, const CC_vecfunction& t, const Pairs<CCPair>& u) const {
    double result = 0.0;
    MADNESS_ASSERT(x.size() == t.size());
    for (const auto& itmp : x.functions) {
        const size_t i = itmp.first;
        for (const auto& ktmp : x.functions) {
            const size_t k = ktmp.first;
            const real_function_3d kgtk = (*g12)(mo_bra_(k), t(k));
            for (const auto& ltmp : x.functions) {
                const size_t l = ltmp.first;
                const real_function_3d lgtk = (*g12)(mo_bra_(l), t(k));
                const real_function_3d k_lgtk = (mo_bra_(k).function * lgtk).truncate();
                const real_function_3d l_kgtk = (mo_bra_(l).function * kgtk).truncate();
                result += (4.0 * make_xy_u(x(i), l_kgtk, get_pair_function(u, i, l)) -
                           2.0 * make_xy_u(l_kgtk, x(i), get_pair_function(u, i, l)) -
                           2.0 * make_xy_u(x(i), k_lgtk, get_pair_function(u, i, l))
                           + 1.0 * make_xy_u(k_lgtk, x(i), get_pair_function(u, i, l)));
            }
        }
    }
    return result;
}

madness::vector_real_function_3d
CCPotentials::s2b(const CC_vecfunction& singles, const Pairs<CCPair>& doubles) const {
    vector_real_function_3d result;
    // see if we can skip the recalculation of the pure 6D part since this does not change during the singles iteration
    vector_real_function_3d result_u = get_potentials(singles, POT_s2b_);
    bool recalc_u_part = false;
    if (result_u.empty()) recalc_u_part = true;

    for (const auto& itmp : singles.functions) {
        const size_t i = itmp.first;
        real_function_3d resulti_u = real_factory_3d(world);
        real_function_3d resulti_r = real_factory_3d(world);
        for (const auto& ktmp : singles.functions) {
            const size_t k = ktmp.first;
            std::vector<CCPairFunction<double,6>> uik = get_pair_function(doubles, i, k);
            // check if the first function in the vector is really the pure 6D part
            MADNESS_ASSERT(uik[0].is_pure());
            if (recalc_u_part) {
                resulti_u += 2.0 * apply_s2b_operation(mo_bra_(k), uik[0],
                                                       2);     //2.0*uik[0].dirac_convolution(mo_bra_(k),g12,2);
                resulti_u -= apply_s2b_operation(mo_bra_(k), uik[0],
                                                 1);     //uik[0].dirac_convolution(mo_bra_(k),g12,1);
            } else {
                resulti_u = result_u[i - parameters.freeze()];
            }
            for (size_t mm = 1; mm < uik.size(); mm++) {
                resulti_r += 2.0 * apply_s2b_operation(mo_bra_(k), uik[mm],
                                                       2);     //2.0*uik[mm].dirac_convolution(mo_bra_(k),g12,2);
                resulti_r -= apply_s2b_operation(mo_bra_(k), uik[mm],
                                                 1);     //uik[mm].dirac_convolution(mo_bra_(k),g12,1);
            }
        }
        result.push_back(resulti_r + resulti_u);
        if (recalc_u_part) result_u.push_back(resulti_u);
    }
    if (recalc_u_part) get_potentials.insert(result_u, singles, POT_s2b_);

    return result;
}

madness::vector_real_function_3d
CCPotentials::s2c(const CC_vecfunction& singles, const Pairs<CCPair>& doubles) const {
    vector_real_function_3d result;
    // see if we can skip the recalculation of the pure 6D part since this does not change during the singles iteration
    vector_real_function_3d result_u = get_potentials(singles, POT_s2c_);
    bool recalc_u_part = false;
    if (result_u.empty()) recalc_u_part = true;

    for (const auto& itmp : singles.functions) {
        const size_t i = itmp.first;
        real_function_3d resulti_u = real_factory_3d(world);
        real_function_3d resulti_r = real_factory_3d(world);
        for (const auto& ktmp : singles.functions) {
            const size_t k = ktmp.first;
            const real_function_3d kgi = (*g12)(mo_bra_(k), mo_ket_(i));
            for (const auto& ltmp : singles.functions) {
                const size_t l = ltmp.first;
                const real_function_3d l_kgi = mo_bra_(l).function * kgi;
                std::vector<CCPairFunction<double,6>> ukl = get_pair_function(doubles, k, l);
                // check if the first function in the vector is really the pure 6D part
                MADNESS_ASSERT(ukl[0].is_pure());
                if (recalc_u_part) {
                    resulti_u += -2.0 * ukl[0].project_out(l_kgi, 2);
                    resulti_u += ukl[0].project_out(l_kgi, 1);
                } else {
                    resulti_u = result_u[i - parameters.freeze()];
                }
                for (size_t mm = 1; mm < ukl.size(); mm++) {
                    resulti_r += -2.0 * ukl[mm].project_out(l_kgi, 2);
                    resulti_r += ukl[mm].project_out(l_kgi, 1);
                }
            }
        }
        result.push_back(resulti_r + resulti_u);
        if (recalc_u_part) result_u.push_back(resulti_u);
    }
    if (recalc_u_part) get_potentials.insert(result_u, singles, POT_s2c_);

    return result;
}

madness::vector_real_function_3d
CCPotentials::s4a_from_s2b(const vector_real_function_3d& s2b, const CC_vecfunction& singles) const {
    MADNESS_ASSERT(singles.type == PARTICLE || singles.type == RESPONSE);
    if (s2b.empty()) output.warning("S2b-potential is empty --> S4a will be zero");

    vector_real_function_3d result;
    for (size_t i = 0; i < s2b.size(); i++) {
        real_function_3d resulti = real_factory_3d(world);
        const Tensor<double> ls2bi = inner(world, s2b[i], get_active_mo_bra());
        for (const auto& ltmp : singles.functions) {
            resulti -= ls2bi[ltmp.first - parameters.freeze()] * singles(ltmp.first).function;
        }
        result.push_back(resulti);
    }
    if (parameters.debug()) {
        vector_real_function_3d s4a_2 = apply_projector(s2b, singles);
        scale(world, s4a_2, -1.0);
        vector_real_function_3d diff = sub(world, result, s4a_2);
        if (world.rank() == 0)
            std::cout << std::fixed << std::setprecision(5) << "||S4a||=" << norm2(world, result) << ", ||S4a||="
                      << norm2(world, s4a_2) << ", ||diff||=" << norm2(world, diff) << "\n";
    }
    // make shure the unprojected S2b potential was used
    if (norm2(world, result) == 0.0) output.warning("S4a potential is zero!! Was the s2b potential Q-projected ?");

    return result;
}

madness::vector_real_function_3d
CCPotentials::s4b(const CC_vecfunction& singles, const Pairs<CCPair>& doubles) const {
    vector_real_function_3d result;
    const vector_real_function_3d active_mo_bra = get_active_mo_bra();
    for (const auto& itmp : singles.functions) {
        const size_t i = itmp.first;
        real_function_3d resulti = real_factory_3d(world);
        for (const auto& ktmp : singles.functions) {
            const size_t k = ktmp.first;
            const real_function_3d kgi = (*g12)(mo_bra_(k), singles(i));
            vector_real_function_3d l_kgi = mul_sparse(world, kgi, active_mo_bra, parameters.thresh_3D());
            truncate(world, l_kgi);
            for (const auto& ltmp : singles.functions) {
                const size_t l = ltmp.first;
                const std::vector<CCPairFunction<double,6>> ukl = get_pair_function(doubles, k, l);
                for (size_t mm = 0; mm < ukl.size(); mm++) {
                    resulti += -2.0 * ukl[mm].project_out(l_kgi[l - parameters.freeze()], 2);
                    resulti += ukl[mm].project_out(l_kgi[l - parameters.freeze()], 1);
                }
            }
        }
        result.push_back(resulti);
    }
    return result;
}

madness::vector_real_function_3d
CCPotentials::s4c(const CC_vecfunction& singles, const Pairs<CCPair>& doubles) const {
    vector_real_function_3d result;
    const vector_real_function_3d active_mo_bra = get_active_mo_bra();
    for (const auto& itmp : singles.functions) {
        const size_t i = itmp.first;
        real_function_3d resulti = real_factory_3d(world);
        real_function_3d part1 = real_factory_3d(world);
        real_function_3d part2 = real_factory_3d(world);
        real_function_3d part3 = real_factory_3d(world);
        real_function_3d part4 = real_factory_3d(world);
        real_function_3d kgtauk = real_factory_3d(world);
        for (const auto& ktmp : singles.functions) {
            const size_t k = ktmp.first;
            kgtauk += (*g12)(mo_bra_(k), singles(k));
        }
        vector_real_function_3d l_kgtauk = mul(world, kgtauk, active_mo_bra);
        truncate(world, l_kgtauk);
        for (const auto& ltmp : singles.functions) {
            const size_t l = ltmp.first;
            const std::vector<CCPairFunction<double,6>> uil = get_pair_function(doubles, i, l);
            for (size_t mm = 0; mm < uil.size(); mm++) {
                part1 += uil[mm].project_out(l_kgtauk[l - parameters.freeze()], 2);
                part2 += uil[mm].project_out(l_kgtauk[l - parameters.freeze()], 1);
            }
            for (const auto& ktmp : singles.functions) {
                const size_t k = ktmp.first;
                const real_function_3d k_lgtauk = (mo_bra_(k).function * (*g12)(mo_bra_(l), singles(k))).truncate();
                for (size_t mm = 0; mm < uil.size(); mm++) {
                    part3 += uil[mm].project_out(k_lgtauk, 2);
                    part4 += uil[mm].project_out(k_lgtauk, 1);
                }
            }
        }
        resulti = 4.0 * part1 - 2.0 * part2 - 2.0 * part3 + part4;
        result.push_back(resulti);
    }
    return result;
}

/// Plotting (convenience)
void CCPotentials::plot(const vector_real_function_3d& f, const std::string& msg) const {
    CCTimer plot_time(world, "plotting " + std::to_string(f.size()) + " functions: " + msg);
    for (size_t k = 0; k < f.size(); k++) plot(f[k], msg + "_" + std::to_string(k));
    plot_time.info();
}

/// Plotting (convenience)
void CCPotentials::plot(const real_function_3d& f, const std::string& msg, const bool doprint) const {
    CCTimer plot_time(world, "plotting " + msg);
    plot_plane(world, f, msg);
    plot_time.info(doprint);
}

/// makes the t intermediates
/// t_i = mo_ket_(i) + factor*tau(i)
/// if factor!=1 then we can not use intermediates and set the type to UNDEFINED
CC_vecfunction CCPotentials::make_t_intermediate(const CC_vecfunction& tau, const double factor) const {
    FuncType returntype = MIXED;
    if (factor != 1.0) returntype = UNDEFINED;

    if (tau.type == HOLE) {
        output("make_t_intermediate: returning hole states");
        return CC_vecfunction(get_active_mo_ket(), HOLE, parameters.freeze());
    }
    if (tau.size() == 0) {
        output("make_t_intermediate: empty tau-> returning hole states");
        return CC_vecfunction(get_active_mo_ket(), HOLE, parameters.freeze());
    }

    CC_vecfunction result(returntype);
    for (const auto& itmp:tau.functions) {
        const size_t i = itmp.first;
        CCFunction<double,3> t(mo_ket_(i).function + factor * tau(i).function, i, MIXED);
        result.insert(i, t);

    }
    return result;
}

/// makes the t intermediates
/// t_i = mo_ket_(i) + factor*tau(i)
/// if the core is frozen the core ti will just be mo_ket_
CC_vecfunction CCPotentials::make_full_t_intermediate(const CC_vecfunction& tau) const {
    FuncType returntype = MIXED;

    if (tau.type == HOLE) {
        output("make_t_intermediate: returning hole states");
        return mo_ket_;
    }
    if (tau.size() == 0) {
        output("make_t_intermediate: empty tau-> returning hole states");
        return mo_ket_;
    }

    CC_vecfunction result(returntype);
    for (size_t i = 0; i < mo_ket_.size(); i++) {
        if (int(i) < parameters.freeze()) {
            result.insert(i, mo_ket_(i));
        } else {
            CCFunction<double,3> t(mo_ket_(i).function + tau(i).function, i, MIXED);
            result.insert(i, t);
        }
    }
    return result;
}

/// makes the t intermediates
/// t_i = mo_ket_(i) + tau
/// i = tau.i
CCFunction<double,3> CCPotentials::make_t_intermediate(const CCFunction<double,3>& tau) const {
    MADNESS_ASSERT(tau.type == PARTICLE);
    const CCFunction<double,3> t(mo_ket_(tau.i).function + tau.function, tau.i, MIXED);
    return t;
}


/// forms the regularized functions from Q and Qt Ansatz for CIS(D) where tau=0 and t=mo so that Qt=Q
void CCPotentials::test_pair_consistency(const CCPairFunction<double,6>& u, const size_t i, const size_t j,
                                         const CC_vecfunction& x) const {
    if (parameters.QtAnsatz()) {
        // u(QAnsatz) = u(QtAnsatz) - OxQftt - QOxftt
        std::vector<CCPairFunction<double,6>> v1;
        v1.push_back(u);
        std::vector<CCPairFunction<double,6>> v2;
        v2.push_back(u);
        CCPairFunction<double,6> ftt(f12, mo_ket_(i), mo_ket_(j));
        CCPairFunction<double,6> O1xftt = apply_Ot(ftt, x, 1);
        CCPairFunction<double,6> OxQftt = apply_Qt(O1xftt, mo_ket_, 2);
        CCPairFunction<double,6> OxQ = OxQftt.invert_sign();
        v2.push_back(OxQ);
        CCPairFunction<double,6> O2xftt = apply_Ot(ftt, x, 2);
        CCPairFunction<double,6> QOxftt = apply_Qt(O2xftt, mo_ket_, 1);
        CCPairFunction<double,6> QOx = QOxftt.invert_sign();
        v2.push_back(QOx);

        CCPair p1(i, j, EXCITED_STATE, CT_CISPD, v1);
        CCPair p2(i, j, EXCITED_STATE, CT_CISPD, v2);
        double norm2_u1 = overlap(p1);
        double norm2_u2 = overlap(p2);

        if (world.rank() == 0)
            std::cout << std::fixed << std::setprecision(10) << "||u(QtAnsatz)||**2=" << norm2_u1
                      << "\n||u(QAnsatz)||**2=" << norm2_u2 << "\n";

    } else {
        // u(QtAnsatz) = u(QAnsatz) + OxQftt - QOxftt
        std::vector<CCPairFunction<double,6>> v1;
        v1.push_back(u);
        std::vector<CCPairFunction<double,6>> v2;
        v2.push_back(u);
        CCPairFunction<double,6> ftt(f12, mo_ket_(i), mo_ket_(j));
        CCPairFunction<double,6> O1xftt = apply_Ot(ftt, x, 1);
        CCPairFunction<double,6> OxQftt = apply_Qt(O1xftt, mo_ket_, 2);
        v2.push_back(OxQftt);
        CCPairFunction<double,6> O2xftt = apply_Ot(ftt, x, 2);
        CCPairFunction<double,6> QOxftt = apply_Qt(O2xftt, mo_ket_, 1);
        v2.push_back(QOxftt);

        CCPair p1(i, j, EXCITED_STATE, CT_CISPD, v1);
        CCPair p2(i, j, EXCITED_STATE, CT_CISPD, v2);
        double norm2_u1 = overlap(p1);
        double norm2_u2 = overlap(p2);

        if (world.rank() == 0)
            std::cout << std::fixed << std::setprecision(10) << "||u(QAnsatz)||**2=" << norm2_u1
                      << "\n||u(QtAnsatz)||**2=" << norm2_u2 << "\n";
    }

}

bool CCPotentials::test_compare_pairs(const CCPair& pair1, const CCPair& pair2) const {
    bool result = true;
    // test overlap
    double ovlp_1 = overlap(pair1);
    double ovlp_2 = overlap(pair2);
    double diff = ovlp_1 - ovlp_2;
    if (world.rank() == 0)
        std::cout << std::fixed << std::setprecision(10)
                  << "||" << pair1.name() << "||**2 =" << ovlp_1 << "\n"
                  << "||" << pair2.name() << "||**2 =" << ovlp_2 << "\n";
    if (fabs(diff) > parameters.thresh_6D()) {
        output.warning("Test Failed, diff=" + std::to_string(diff));
        result = false;
    } else output("Test Passed, diff=" + std::to_string(diff));

    // test energy integration <ij|Qf|ij>
    double energy_1 = make_xy_op_u(mo_bra_(pair1.i), mo_bra_(pair1.j), *g12, pair1.functions);
    double energy_2 = make_xy_op_u(mo_bra_(pair2.i), mo_bra_(pair2.j), *g12, pair2.functions);
    double diff_energy = energy_1 - energy_2;
    if (world.rank() == 0)
        std::cout << std::fixed << std::setprecision(10)
                  << "<ij|g|" << pair1.name() << "> =" << energy_1 << "\n"
                  << "<ij|g|" << pair2.name() << "> =" << energy_2 << "\n";
    if (fabs(diff_energy) > parameters.thresh_6D()) {
        output.warning("Test Failed, diff=" + std::to_string(diff_energy));
        result = false;
    } else output("Test Passed, diff=" + std::to_string(diff_energy));


    // make full 6D pair and compare
    real_function_6d pair1_6D = make_6D_pair(pair1);
    real_function_6d pair2_6D = make_6D_pair(pair2);
    real_function_6d pair_diff = pair1_6D - pair2_6D;
    pair1_6D.print_size(pair1.name() + "_6D");
    pair2_6D.print_size(pair2.name() + "_6D");
    pair_diff.print_size("diff");
    if (pair_diff.norm2() > parameters.thresh_6D()) {
        result = false;
        output.warning("Test Failed, difference of pairs is not zero");
    } else output("Test Passed");

    return result;
}

// make a single 6D functions of the pair
real_function_6d CCPotentials::make_6D_pair(const CCPair& pair) const {
    std::vector<CCPairFunction<double,6>> functions = pair.functions;
    real_function_6d result = real_factory_6d(world);
    for (const auto& f:functions) {
        if (f.is_pure()) result += f.pure().get_function();
        else if (f.is_decomposed_no_op()) {
            for (size_t i = 0; i < f.get_a().size(); i++) {
                real_function_6d ab = CompositeFactory<double, 6, 3>(world).particle1(copy(f.get_a()[i])).particle2(
                        copy(f.get_b()[i]));
                ab.fill_tree().truncate().reduce_rank();
                result += ab;
            }
        } else if (f.is_op_decomposed()) {
            MADNESS_ASSERT(f.get_operator().type() == OpType::OT_F12);
            real_function_6d fxy = make_f_xy(f.get_a()[0], f.get_b()[0]);
            result += fxy;
        } else MADNESS_EXCEPTION("Unknown type of CCPairFunction<double,6>", 1);
    }
    return result;
}

void CCPotentials::test_pairs() {

    {
        output.section("Testing GS Reg Tails");
        CCTimer time(world, "Testing GS Reg Tails");
        bool regrestest = true;
        std::vector<size_t> testorbs;
        // make test with homo and core orbitals
        real_function_6d u = real_factory_6d(world);
        if (mo_ket_.size() - parameters.freeze() == 1) testorbs.push_back(mo_ket_.size() - 1);
        else {
            testorbs.push_back(parameters.freeze());
            testorbs.push_back(mo_ket_.size() - 1);
        }
        const CC_vecfunction t(get_active_mo_ket(), HOLE, parameters.freeze());
        for (const auto i:testorbs) {
            for (const auto j:testorbs) {
                CCPair rr_3D = (make_pair_gs(u, mo_ket_, i, j));
                real_function_6d rr_tmp = make_f_xy(t(i), t(j));
                real_function_6d rr_6D0 = apply_Q12t(rr_tmp, mo_ket_);
                CCPairFunction<double,6> rr_6D1(rr_6D0);
                std::vector<CCPairFunction<double,6>> rr_6D2(1, rr_6D1);
                CCPair rr_6D(i, j, GROUND_STATE, CT_CC2, rr_6D2);
                regrestest = test_compare_pairs(rr_3D, rr_6D);
            }
        }
        if (regrestest) output("All Tests Passed");
        else output.warning("NOT all Tests Passed");
        time.info();
    }

    {
        output.section("Testing ES Reg Tails");
        CCTimer time(world, "Testing GS Reg Tails");
        real_function_3d xf = real_factory_3d(world).f(functor_r2);
        vector_real_function_3d xmo = mul(world, xf, get_active_mo_ket());
        xmo = apply_Qt(xmo, mo_ket_);
        const double norm = norm2(world, xmo);
        scale(world, xmo, 1 / norm);
        CC_vecfunction x(xmo, RESPONSE, parameters.freeze());
        const double omega = -0.9 * get_orbital_energies()[mo_ket_.size() - 1];
        x.omega = omega;
        bool regrestest = true;
        std::vector<size_t> testorbs;
        // make test with homo and core orbitals
        real_function_6d u = real_factory_6d(world);
        if (mo_ket_.size() - parameters.freeze() == 1) testorbs.push_back(mo_ket_.size() - 1);
        else {
            testorbs.push_back(parameters.freeze());
            testorbs.push_back(mo_ket_.size() - 1);
        }
        const CC_vecfunction t(get_active_mo_ket(), HOLE, parameters.freeze());
        for (const auto i:testorbs) {
            for (const auto j:testorbs) {
                CCPair rr_3D = (make_pair_ex(u, mo_ket_, x, i, j, CT_CISPD));
                real_function_6d rr_6D0;
                {
                    real_function_6d rr_ftx = make_f_xy(x(i), t(j));
                    real_function_6d rr_fxt = make_f_xy(t(i), x(j));
                    real_function_6d rr_f = rr_ftx + rr_fxt;
                    rr_6D0 = apply_Q12t(rr_f, mo_ket_);
                    if (parameters.QtAnsatz()) {
                        MADNESS_EXCEPTION("No Test for QtAnsatz right now", 1);
                    }
                }
                CCPairFunction<double,6> rr_6D1(rr_6D0);
                std::vector<CCPairFunction<double,6>> rr_6D2(1, rr_6D1);
                CCPair rr_6D(i, j, GROUND_STATE, CT_MP2, rr_6D2);
                regrestest = test_compare_pairs(rr_3D, rr_6D);
            }
        }
        if (regrestest) output("All Tests Passed");
        else output.warning("NOT all Tests Passed");
        time.info();
    }

}

void CCPotentials::test_singles_potential() const {

    output("Test LRCC2 Singles Potential with empty doubles and compare to CIS");
    {
        vector_real_function_3d emptyv = zero_functions<double, 3>(world, get_active_mo_ket().size());
        real_function_3d r = real_factory_3d(world).f(functor_y);
        const CC_vecfunction gs_singles(emptyv, PARTICLE, parameters.freeze());
        CC_vecfunction ex_singles(apply_Qt(mul(world, r, get_active_mo_ket()), mo_ket_), RESPONSE, parameters.freeze());
        Pairs<CCPair> gs_doubles;
        Pairs<CCPair> ex_doubles;

        vector_real_function_3d cis_potential = potential_singles_ex(gs_singles, gs_doubles, ex_singles, ex_doubles,
                                                                     POT_cis_);
        vector_real_function_3d ccs_potential = potential_singles_ex(gs_singles, gs_doubles, ex_singles, ex_doubles,
                                                                     POT_ccs_);
        vector_real_function_3d diff = sub(world, cis_potential, ccs_potential);
        const double d = norm2(world, diff);
        madness::print_size<double, 3>(world, diff, "difference in potentials");
        if (world.rank() == 0) std::cout << "||difference|| = " << d << "\n";

    }

    // Make functions
    real_function_3d x = real_factory_3d(world).f(functor_x);
    real_function_3d r = real_factory_3d(world).f(functor_y);
    const CC_vecfunction gs_singles(apply_Qt(mul(world, x, get_active_mo_ket()), mo_ket_), PARTICLE, parameters.freeze());
    CC_vecfunction ex_singles(apply_Qt(mul(world, r, get_active_mo_ket()), mo_ket_), RESPONSE, parameters.freeze());
    ex_singles.omega = -0.9 * get_orbital_energies().back();
    Pairs<CCPair> gs_doubles;
    Pairs<CCPair> ex_doubles;
    for (size_t i = parameters.freeze(); i < mo_ket_.size(); i++) {
        for (size_t j = i; j < mo_ket_.size(); j++) {
            real_function_6d tmp = real_factory_6d(world);// make_f_xy(gs_singles(i),gs_singles(j));
            real_function_6d tmpx = real_factory_6d(world);//make_f_xy(ex_singles(i),ex_singles(j));
            CCPair ptmp = make_pair_gs(tmp, gs_singles, i, j);
            CCPair xtmp = make_pair_ex(tmpx, gs_singles, ex_singles, i, j, CT_LRCC2);
            gs_doubles.insert(i, j, ptmp);
            ex_doubles.insert(i, j, xtmp);
        }
    }

    std::vector<PotentialType> pots = {POT_ccs_, POT_F3D_, POT_s2b_, POT_s2c_, POT_s4b_, POT_s4c_};

    output.section("Testing of CC2 Singles Ground State Potential");
    CCTimer time_gs(world, "CC2 Singles GS Test");

    vector_real_function_3d tmp = mul(world, nemo_->ncf->square(), ex_singles.get_vecfunction());
    truncate(world, tmp);
    const CC_vecfunction xbra(tmp, RESPONSE, parameters.freeze());

    for (const auto pot:pots) {
        const vector_real_function_3d potential = potential_singles_gs(gs_singles, gs_doubles, pot);
        const double xpot1 = inner(world, xbra.get_vecfunction(), potential).sum();
        const double xpot2 = potential_energy_gs(xbra, gs_singles, gs_doubles, pot);
        const double diff = xpot1 - xpot2;
        if (world.rank() == 0)
            std::cout << std::fixed << std::setprecision(10) <<
                      "Testing of Potential " << assign_name(pot) << ":" <<
                      std::setw(10) << "\n<x|" + assign_name(pot) + "> = " << xpot1 <<
                      std::setw(10) << "\n<x|" + assign_name(pot) + "> = " << xpot2 <<
                      std::setw(10) << "\ndiff = " << diff << "\n";
        if (fabs(diff) > parameters.thresh_6D()) output.warning("Test Failed");
        else output("Test Passed");
        if (pot == POT_s2b_) {
            const vector_real_function_3d pot_s4a = -1.0 * apply_projector(potential, gs_singles);
            const double xxpot1 = inner(world, xbra.get_vecfunction(), pot_s4a).sum();
            const double xxpot2 = potential_energy_gs(xbra, gs_singles, gs_doubles, POT_s4a_);
            const double xdiff = xxpot1 - xxpot2;
            if (world.rank() == 0)
                std::cout <<
                          std::setw(10) << "\n<x|" + assign_name(POT_s4a_) + "> = " << xxpot1 <<
                          std::setw(10) << "\n<x|" + assign_name(POT_s4a_) + "> = " << xxpot2 <<
                          std::setw(10) << "\ndiff = " << xdiff << "\n";
            if (fabs(xdiff) > parameters.thresh_6D()) output.warning("Test Failed");
            else output("Test Passed");
        }
    }

    time_gs.info();

    output.section("Testing of CC2 Singles Response Potential");
    CCTimer time_ex(world, "CC2 Singles Response Test");

    for (const auto pot:pots) {
        const vector_real_function_3d potential = potential_singles_ex(gs_singles, gs_doubles, ex_singles, ex_doubles,
                                                                       pot);
        const double xpot1 = inner(world, xbra.get_vecfunction(), potential).sum();
        const double xpot2 = potential_energy_ex(xbra, gs_singles, gs_doubles, ex_singles, ex_doubles, pot);
        const double diff = xpot1 - xpot2;
        if (world.rank() == 0)
            std::cout << std::fixed << std::setprecision(10) <<
                      "Testing of Response Potential " << assign_name(pot) << ":" <<
                      std::setw(10) << "\n<x|" + assign_name(pot) + "> = " << xpot1 <<
                      std::setw(10) << "\n<x|" + assign_name(pot) + "> = " << xpot2 <<
                      std::setw(10) << "\ndiff = " << diff << "\n";
        if (fabs(diff) > parameters.thresh_6D()) output.warning("Test Failed");
        else output("Test Passed");
        if (pot == POT_s2b_) {
            const vector_real_function_3d potential_gs = potential_singles_gs(gs_singles, gs_doubles, pot);
            const vector_real_function_3d pot_s4a = -1.0 * add(world, apply_projector(potential, gs_singles),
                                                               apply_projector(potential_gs, ex_singles));
            const double xxpot1 = inner(world, xbra.get_vecfunction(), pot_s4a).sum();
            const double xxpot2 = potential_energy_ex(xbra, gs_singles, gs_doubles, ex_singles, ex_doubles, POT_s4a_);
            const double xdiff = xxpot1 - xxpot2;
            if (world.rank() == 0)
                std::cout <<
                          std::setw(10) << "\n<x|" + assign_name(POT_s4a_) + "> = " << xxpot1 <<
                          std::setw(10) << "\n<x|" + assign_name(POT_s4a_) + "> = " << xxpot2 <<
                          std::setw(10) << "\ndiff = " << xdiff << "\n";
            if (fabs(xdiff) > parameters.thresh_6D()) output.warning("Test Failed");
            else output("Test Passed");
        }
    }

    time_ex.info();
}

void CCPotentials::test() {
    output.section("Testing enums");
    CalcType test2 = CT_MP2;
    FuncType test4 = HOLE;
    CCState test5 = GROUND_STATE;
    PotentialType test6 = POT_F3D_;
    assign_name(test2);
    assign_name(test4);
    assign_name(test5);
    assign_name(test6);

    test_singles_potential();
    output.section("Testing Scalar Multiplication");
    {
        CC_vecfunction test = mo_ket_ * 2.0;
        double norma = norm2(world, mo_ket_.get_vecfunction());
        double normb = norm2(world, test.get_vecfunction());
        if (2.0 * norma == normb) output("Test Passed");
        else output.warning("Test Failed: Norm1 = " + std::to_string(norma) + ", Norm2 = " + std::to_string(normb));
    }
    {
        CCFunction<double,3> mo = mo_ket_(0);
        CCFunction<double,3> mo1 = mo * 2.0;
        double norma = mo.function.norm2();
        double normb = mo1.function.norm2();
        if (2.0 * norma == normb) output("Test Passed");
        else output.warning("Test Failed: Norm1 = " + std::to_string(norma) + ", Norm2 = " + std::to_string(normb));
    }

    test_pairs();
    output.section("Testing destructiveness of project out function");
    {
        real_function_3d bra = mo_bra_(mo_bra_.size() - 1).function;
        real_function_6d ftt = make_f_xy(mo_ket_(mo_ket_.size() - 1), mo_ket_(mo_ket_.size() - 1));
        ftt.print_size("f|homo,homo>");
        bra.print_size("<homo|");
        real_function_3d test1 = ftt.project_out(bra, 0);
        test1.print_size("<homo|f12|homo,homo>_1");
        ftt.print_size("f|homo,homo>");
        bra.print_size("<homo|");
        real_function_3d test2 = ftt.project_out(bra, 1);
        test2.print_size("<homo|f12|homo,homo>_2");
        ftt.print_size("f|homo,homo>");
        bra.print_size("<homo|");
    }


    output.section("Testing Overlaps of CCPairFunction<double,6>");
    // f^2 = 1/(4y^2)(1 - 2*f2(y) + f2(2y)) , f2(2y) =f2(y)^2
    const double y = parameters.gamma();
    SeparatedConvolution<double, 3> f = SlaterF12Operator(world, y, parameters.lo(), parameters.thresh_bsh_3D());
    SeparatedConvolution<double, 3> f2 = SlaterOperator(world, y, parameters.lo(), parameters.thresh_bsh_3D());
    SeparatedConvolution<double, 3> ff = SlaterOperator(world, 2.0 * y, parameters.lo(), parameters.thresh_bsh_3D());

    real_function_3d a = mo_ket_(0).function;
    real_function_3d b = mo_ket_(0).function;
    real_function_6d fab_6d = CompositeFactory<double, 6, 3>(world).g12(corrfac.f()).particle1(copy(a)).particle2(
            copy(b));
    fab_6d.fill_tree().truncate().reduce_rank();
    {
        real_function_3d aR = mo_bra_(0).function;
        real_function_3d bR = mo_bra_(0).function;
        real_function_6d ab_6d = CompositeFactory<double, 6, 3>(world).particle1(copy(a)).particle2(copy(b));
        ab_6d.fill_tree().truncate().reduce_rank();

        const real_function_3d aa = (aR * a).truncate();
        const real_function_3d bb = (bR * b).truncate();
        const real_function_3d af2a = f2(aa);
        const real_function_3d affa = ff(aa);
        const real_function_3d afa = f(aa);

        bool passed_lo = true;
        bool passed_hi = true;
        const double lo = parameters.thresh_6D();
        const double hi = parameters.thresh_3D();
        {
            CCPairFunction<double,6> fab(f12, a, b);
            const double test1 = overlap(fab, fab);
            const double prefactor = 1.0 / (4 * y * y);
            const double test2 = prefactor * (aR.inner(a) * bR.inner(b) - 2.0 * bb.inner(af2a) + bb.inner(affa));
            const double diff = test1 - test2;
            if (fabs(diff) > lo) passed_lo = false;
            if (fabs(diff) > hi) passed_hi = false;


            if (world.rank() == 0)
                std::cout << "Overlap Test 1 : " << std::fixed << std::setprecision(10) << "result=" << test1
                          << ", test=" << test2 << ", diff=" << diff << "\n";
        }
        {
            CCPairFunction<double,6> ab(mo_ket_.get_vecfunction(), mo_ket_.get_vecfunction());
            const double test1 = overlap(ab, ab);
            const double test2 = double(mo_ket_.size()); // mos are normed
            const double diff = test1 - test2;
            if (world.rank() == 0)
                std::cout << "Overlap Test 2 : " << std::fixed << std::setprecision(10) << "result=" << test1
                          << ", test=" << test2 << ", diff=" << diff << "\n";
            if (fabs(diff) > lo) passed_lo = false;
            if (fabs(diff) > hi) passed_hi = false;
        }
        {
            CCPairFunction<double,6> ab(ab_6d);
            const double test1 = overlap(ab, ab);
            const double test2 = double(mo_ket_.size()); // mos are normed
            const double diff = test1 - test2;
            if (world.rank() == 0)
                std::cout << "Overlap Test 3 : " << std::fixed << std::setprecision(10) << "result=" << test1
                          << ", test=" << test2 << ", diff=" << diff << "\n";
            if (fabs(diff) > lo) passed_lo = false;
            if (fabs(diff) > hi) passed_hi = false;
        }
        {
            CCPairFunction<double,6> ab1(mo_ket_.get_vecfunction(), mo_ket_.get_vecfunction());
            CCPairFunction<double,6> ab2(ab_6d);
            const double test1 = overlap(ab1, ab2);
            const double test2 = double(mo_ket_.size()); // mos are normed
            const double diff = test1 - test2;
            if (world.rank() == 0)
                std::cout << "Overlap Test 4 : " << std::fixed << std::setprecision(10) << "result=" << test1
                          << ", test=" << test2 << ", diff=" << diff << "\n";
            if (fabs(diff) > lo) passed_lo = false;
            if (fabs(diff) > hi) passed_hi = false;
        }
        {
            // the next tests evaulate <ab|f|ab> in different ways
            CCPairFunction<double,6> fab(fab_6d);
            CCPairFunction<double,6> ab2(mo_ket_.get_vecfunction(), mo_ket_.get_vecfunction());
            const double test1 = overlap(fab, ab2);
            const double test2 = bb.inner(afa);
            const double diff = test1 - test2;
            if (world.rank() == 0)
                std::cout << "Overlap Test 5 : " << std::fixed << std::setprecision(10) << "result=" << test1
                          << ", test=" << test2 << ", diff=" << diff << "\n";
            if (fabs(diff) > lo) passed_lo = false;
            if (fabs(diff) > hi) passed_hi = false;
        }
        {
            CCPairFunction<double,6> fab(fab_6d);
            CCPairFunction<double,6> ab2(ab_6d);
            const double test1 = overlap(fab, ab2);
            const double test2 = bb.inner(afa);
            const double diff = test1 - test2;
            if (world.rank() == 0)
                std::cout << "Overlap Test 6 : " << std::fixed << std::setprecision(10) << "result=" << test1
                          << ", test=" << test2 << ", diff=" << diff << "\n";
            if (fabs(diff) > lo) passed_lo = false;
            if (fabs(diff) > hi) passed_hi = false;
        }
        {
            CCPairFunction<double,6> fab(f12, a, b);
            CCPairFunction<double,6> ab2(ab_6d);
            const double test1 = overlap(fab, ab2);
            const double test2 = bb.inner(afa);
            const double diff = test1 - test2;
            if (world.rank() == 0)
                std::cout << "Overlap Test 7 : " << std::fixed << std::setprecision(10) << "result=" << test1
                          << ", test=" << test2 << ", diff=" << diff << "\n";
            if (fabs(diff) > lo) passed_lo = false;
            if (fabs(diff) > hi) passed_hi = false;
        }
        {
            CCPairFunction<double,6> fab(f12, a, b);
            CCPairFunction<double,6> ab2(mo_ket_.get_vecfunction(), mo_ket_.get_vecfunction());
            const double test1 = overlap(fab, ab2);
            const double test2 = bb.inner(afa);
            const double diff = test1 - test2;
            if (world.rank() == 0)
                std::cout << "Overlap Test 8 : " << std::fixed << std::setprecision(10) << "result=" << test1
                          << ", test=" << test2 << ", diff=" << diff << "\n";
            if (fabs(diff) > lo) passed_lo = false;
            if (fabs(diff) > hi) passed_hi = false;
        }
        output.section("Overlap Testing Ended");
        if (not passed_lo) output.warning("There were Complications!");
        if (not passed_hi) output.warning("A bit shaky  ... Expected at the pure 6D level");


    }
}

} /* namespace madness */
