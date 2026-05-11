//
// exchange_commutator.cc — Phase A implementation.
//
// Every apply_KffK_* entry point is a thin wrapper around an existing
// CCPotentials static method, extended with uniform wall-time / memory /
// rank instrumentation packaged in a KffKResult.  The scalar split-α
// routine and the harmonic-basis diagnostic are extracted verbatim from
// benchmark_exchange_commutator.cc.  No existing behaviour is duplicated.
//

#include <madness/chem/exchange_commutator.h>

#include <madness/chem/CCPotentials.h>
#include <madness/chem/SCFOperators.h>   // madness::Exchange
#include <madness/mra/operator.h>        // SlaterF12OperatorPtr_ND / SlaterOperatorPtr_ND
#include <madness/world/timing_utilities.h>
#include <madness/world/timers.h>

namespace madness {

namespace {

/// 6D pair size in GByte: sum over decomposed factors or the pure 6D part.
double pair_size_gb(World& world, const std::vector<CCPairFunction<double,6>>& v) {
    double total = 0.0;
    for (const auto& cc : v) {
        if (!cc.is_assigned()) continue;
        if (cc.is_pure()) {
            total += get_size(cc.get_function());
        } else if (cc.is_decomposed()) {
            total += get_size(world, cc.get_a()) + get_size(world, cc.get_b());
        }
    }
    return total;
}

long pair_rank(const std::vector<CCPairFunction<double,6>>& v) {
    long R = 0;
    for (const auto& cc : v) {
        if (cc.is_decomposed() && cc.is_assigned()) R += cc.get_a().size();
    }
    return R;
}

/// Sum the three named pieces into the aggregate mem / rank fields.
void finalize_sizes(World& world, ExchangeCommutator::KffKResult& r) {
    r.mem_gb = pair_size_gb(world, r.Kf)
             + pair_size_gb(world, r.fK)
             + pair_size_gb(world, r.KffK);
    r.rank   = pair_rank(r.Kf) + pair_rank(r.fK) + pair_rank(r.KffK);
}

struct wall_timer {
    World& world;
    double start;
    explicit wall_timer(World& w) : world(w) {
        world.gop.fence();
        start = wall_time();
    }
    double elapsed() {
        world.gop.fence();
        return wall_time() - start;
    }
};

// Apply a partial-Coulomb K̂(phi) = Σ_k k_ket * op_partial(k_bra * phi).
// op_partial is a 3D SeparatedConvolution holding only a subset of the
// Coulomb GFit (e.g. just the medium Gaussians).  Mirrors K_macrotask but
// substitutes op_partial for the full Coulomb.
real_function_3d K_partial_3d(
        World& world,
        const std::vector<real_function_3d>& mo_ket,
        const std::vector<real_function_3d>& mo_bra,
        const real_function_3d& f,
        SeparatedConvolution<double, 3>& op_partial) {
    op_partial.particle() = 1;
    real_function_3d result = real_factory_3d(world);
    for (size_t k = 0; k < mo_ket.size(); ++k) {
        result += op_partial(mo_bra[k] * f).truncate() * mo_ket[k];
    }
    return result;
}

// Apply Kf|xy⟩ in 6D using a partial-Coulomb operator.  Body cloned from
// CCPotentials::apply_Kfxy with g12 → op_partial; debug test-lambdas
// dropped (they returned 0.0 unconditionally upstream).
real_function_6d apply_Kfxy_partial(
        World& world,
        const CCFunction<double, 3>& x,
        const CCFunction<double, 3>& y,
        const Info& info,
        const CCParameters& parameters,
        SeparatedConvolution<double, 3>& op_partial) {
    CorrelationFactor corrfac(world, parameters.gamma(), 1.e-7, parameters.lo());
    op_partial.destructive() = true;

    real_function_6d result = FunctionFactory<double, 6>(world);

    for (int particle : {1, 2}) {
        const auto& kbra = info.mo_bra;
        const auto k_arg = (particle == 1) ? kbra * x.function : kbra * y.function;

        const std::size_t batchsize = 3;
        for (std::size_t kbatch = 0; kbatch < info.mo_ket.size(); kbatch += batchsize) {
            for (std::size_t k = kbatch;
                 k < std::min(kbatch + batchsize, info.mo_ket.size()); ++k) {
                real_function_3d xx = (particle == 1) ? k_arg[k] : x.function;
                real_function_3d yy = (particle == 2) ? k_arg[k] : y.function;
                xx.truncate();
                yy.truncate();
                real_function_6d X = CompositeFactory<double, 6, 3>(world)
                                        .g12(corrfac.f())
                                        .particle1(copy(xx))
                                        .particle2(copy(yy));
                X.fill_cuspy_tree().truncate(parameters.tight_thresh_6D()).reduce_rank();

                op_partial.particle() = particle;
                real_function_6d Y = op_partial(X);
                auto tmp = (multiply(copy(Y), copy(info.mo_ket[k]), particle))
                                .truncate(parameters.tight_thresh_6D() * 3.0);
                result += tmp;
            }
            result.truncate(parameters.tight_thresh_6D())
                  .reduce_rank(parameters.tight_thresh_6D());
        }

        if (x.i == y.i) {
            result += madness::swap_particles(result);
            break;
        }
    }
    return result;
}

} // namespace

// ---------------------------------------------------------------------------
//  6D path — mirrors test_kcomm_6d in the benchmark file
// ---------------------------------------------------------------------------

ExchangeCommutator::KffKResult
ExchangeCommutator::apply_KffK_6d(
        World& world,
        const CCFunction<double,3>& phi_i,
        const CCFunction<double,3>& phi_j,
        const Info& info,
        const real_convolution_6d* Gscreen)
{
    KffKResult out;
    out.algo = "6d";
    wall_timer t(world);

    const CCParameters& parameters = info.parameters;

    // --- Kf|ij> via apply_Kfxy ---------------------------------------------
    real_function_6d Kfxy = CCPotentials::apply_Kfxy(world, phi_i, phi_j, info, parameters);
    Kfxy.truncate().reduce_rank();

    // --- fK|ij> via K_macrotask + make_f_xy_macrotask ----------------------
    const real_function_3d x_ket = phi_i.function;
    const real_function_3d y_ket = phi_j.function;
    const real_function_3d x_bra = (info.R_square * phi_i.function).truncate();
    const real_function_3d y_bra = (info.R_square * phi_j.function).truncate();

    const bool symmetric_fk = (phi_i == phi_j);
    const real_function_3d Kx = CCPotentials::K_macrotask(
            world, info.mo_ket, info.mo_bra, x_ket, parameters);
    const FuncType Kx_type = UNDEFINED;
    const real_function_6d fKphi0b = CCPotentials::make_f_xy_macrotask(
            world, Kx, y_ket, x_bra, y_bra, phi_i.i, phi_j.i,
            parameters, Kx_type, phi_j.type, Gscreen);
    real_function_6d fKphi0a;
    if (symmetric_fk) {
        fKphi0a = madness::swap_particles(fKphi0b);
    } else {
        const real_function_3d Ky = CCPotentials::K_macrotask(
                world, info.mo_ket, info.mo_bra, y_ket, parameters);
        const FuncType Ky_type = UNDEFINED;
        fKphi0a = CCPotentials::make_f_xy_macrotask(
                world, x_ket, Ky, x_bra, y_bra, phi_i.i, phi_j.i,
                parameters, phi_i.type, Ky_type, Gscreen);
    }
    real_function_6d fKxy = (fKphi0a + fKphi0b);
    fKxy.truncate().reduce_rank();

    out.Kf   = { CCPairFunction<double,6>(Kfxy) };
    out.fK   = { CCPairFunction<double,6>(fKxy) };
    out.KffK = { CCPairFunction<double,6>(Kfxy),
                 -1.0 * CCPairFunction<double,6>(fKxy) };
    finalize_sizes(world, out);
    out.t_wall = t.elapsed();
    return out;
}


// ---------------------------------------------------------------------------
//  Full 6D split-α assembly — loops over mo_ket, produces a KffKResult
//  that diagnose/print_report can consume.
// ---------------------------------------------------------------------------

ExchangeCommutator::KffKResult
ExchangeCommutator::apply_KffK_lowrank_split_alpha(
        World& world,
        const CCFunction<double, 3>& phi_i,
        const CCFunction<double, 3>& phi_j,
        const Info& info,
        const LowRankFunctionParameters& lrfparam,
        const SplitAlphaOptions& opt)
{
    constexpr std::size_t LDIM = 3;
    constexpr std::size_t NDIM = 2 * LDIM;
    KffKResult out;
    out.algo = "lrf-split-alpha";
    timer t(world);
    wall_timer t1(world);

    const double thresh = FunctionDefaults<LDIM>::get_thresh();
    const bool symmetric = (phi_i == phi_j);

    auto f12_cc = CCConvolutionOperatorPtr<double, LDIM>(
            world, OT_F12, info.parameters);
    auto f12ptr = f12_cc->get_op();

    // K̂ kernel (Coulomb 1/r) via OperatorInfo; small-α subset built below.
    GFit<double, LDIM> fit = GFit<double, LDIM>::CoulombFit(
            opt.lo, opt.hi, opt.eps_gfit, false);
    const Tensor<double> c_all = fit.coeffs();
    const Tensor<double> a_all = fit.exponents();
    const long M = c_all.size();
    t.tag("fit operator expansion");

    std::vector<double> cs, as;
    for (long mu = 0; mu < M; ++mu) {
        if (a_all[mu] <= opt.alpha_star) {
            MADNESS_ASSERT(c_all[mu] > 0.0);  // risk #3 in the plan
            cs.push_back(c_all[mu]);
            as.push_back(a_all[mu]);
        }
    }
    const long Mk = cs.size();
    Tensor<double> c_t(Mk), a_t(Mk);
    for (long mu = 0; mu < Mk; ++mu) { c_t[mu] = cs[mu]; a_t[mu] = as[mu]; }
    auto trunc_op = std::make_shared<SeparatedConvolution<double, LDIM>>(
            world, c_t, a_t, opt.lo, thresh);

    std::vector<Vector<double, LDIM>> origins=info.molecular_coordinates;

    // The two commutator pieces are each stored as a single CCPairFunction
    // stack, so they compose cleanly into Kf / fK / KffK.  The Kf piece is a
    // plain decomposed pair `g_ρ(1) · (phi_j · A_ρ)(2)` (the f12 is already
    // baked into A_ρ); the fK piece must be wrapped with f12 because its
    // factors `C_i(1)` and `phi_j(2)` carry no f12 themselves.  The minus
    // sign of the commutator is applied on the KffK stack only.

    // The nemo exchange kernel is K(r,r') = Σ_k k(r)·R²(r')·k(r')/|r-r'|,
    // i.e. the bra-side R²-weighted orbital lives on the *integration*
    // coordinate.  In the LRF picture (r,r') → (particle 1, particle 2),
    // particle 2 is the integration coord, so k_bra (= R²·k) must be the
    // particle-2 multiplier.  When R²=1 the two arguments coincide, which
    // hid this asymmetry until non-trivial NCFs were tested.
    LRFunctorF12<double, NDIM> functor(trunc_op, info.mo_ket, info.mo_bra);
    auto lrf_exchange_op = LowRankFunctionFactory<double, NDIM>(lrfparam, origins)
               .project(functor, FunctionDefaults<6>::get_thresh(), 0);
    t.tag("construct LRF of exchange kernel");

    auto& gvec = lrf_exchange_op.g;
    auto& hvec = lrf_exchange_op.h;
    double tol=FunctionDefaults<LDIM>::get_thresh();

    // piece 1 (K̂₁ f part):  Σ_ρ g_ρ(1) · (phi_j · A_ρ)(2)
    // A_ρ(r₂) = f12(h_ρ · phi_i)(r₂)
    auto j_A_pi = phi_j.function * apply(world, *f12ptr, hvec * phi_i.function);
    auto Kf=LowRankFunction<double,NDIM>(gvec,j_A_pi,tol);
    if (symmetric) {
        Kf+=LowRankFunction<double,NDIM>(j_A_pi,gvec,tol);
    } else {
        auto i_A_pj = phi_i.function * apply(world, *f12ptr, hvec * phi_j.function);
        Kf+=LowRankFunction<double,NDIM>(i_A_pj,gvec,tol);
    }
    t.tag("construct Kf |ij>");

    // piece 2 (f K̂₁ part):  Σ_ρ B_ρi · g_ρ(1) · f(1,2) · phi_j(2)
    //                       =  C_i(1) f(1,2) phi_j(2)
    // B_ρ     = ⟨h_ρ | phi_i⟩
    // NOTE: for internal consistency the fK term uses the truncated exchange kernel
    // instead of simply applying the exchange operator on the ket.
    auto B_pi = inner(world, phi_i.function, hvec);
    auto C_i=transform(world,gvec,B_pi);
    MADNESS_CHECK_THROW(C_i.size()==1,"invalid size for Ci");
    auto fK=LowRankFunction<double,6>(C_i,{phi_j.function},tol);
    if (symmetric) {
        fK+=LowRankFunction<double,NDIM>({phi_j.function},C_i,tol);
    } else {
        auto B_pj = inner(world, phi_j.function, hvec);
        auto C_j=transform(world,gvec,B_pj);
        fK+=LowRankFunction<double,NDIM>({phi_i.function},C_j,tol);
    }
    t.tag("construct fK |ij>");

    // Assemble named pair stacks.  Kf has f12 absorbed into the per-ρ factor
    // `j_A_pi = phi_j · f12(h_ρ · phi_i)`, so it is a plain decomposed pair.
    // fK's factors `C_i(1)` and `phi_j(2)` are bare orbitals; the pair must
    // be wrapped with f12_cc to represent C_i(1) · f(1,2) · phi_j(2).
    out.Kf = { CCPairFunction<double, 6>(Kf.get_g(), Kf.get_h()) };
    out.fK = { CCPairFunction<double, 6>(f12_cc,    fK.get_g(), fK.get_h()) };
    out.fK[0].convert_to_pure_no_op_inplace();
    t.tag("convert fK to pure");

    out.KffK.push_back(out.Kf[0]);
    out.KffK.push_back(-1.0 * out.fK[0]);

    finalize_sizes(world, out);
    out.rank   = gvec.size();   // LRF rank of the underlying decomposition
    out.t_wall = t1.elapsed();
    return out;
}

// ---------------------------------------------------------------------------
//  Three-range LRF assembly: diffuse + medium sub-LRFs, tight discarded.
//  Mirrors apply_KffK_lowrank_split_alpha but builds two LRFs (with separate
//  LowRankFunctionParameters) and concatenates their g/h vectors before
//  assembling Kf and fK.  See lrf_three_range_gfit.md.
// ---------------------------------------------------------------------------

ExchangeCommutator::KffKResult
ExchangeCommutator::apply_KffK_lowrank_three_range(
        World& world,
        const CCFunction<double, 3>& phi_i,
        const CCFunction<double, 3>& phi_j,
        const Info& info,
        const ThreeRangeOptions& opt)
{
    constexpr std::size_t LDIM = 3;
    constexpr std::size_t NDIM = 2 * LDIM;
    KffKResult out;
    out.algo = "lrf-three-range";
    timer t(world);
    wall_timer t1(world);

    const double thresh = FunctionDefaults<LDIM>::get_thresh();
    const bool symmetric = (phi_i == phi_j);

    auto f12_cc = CCConvolutionOperatorPtr<double, LDIM>(
            world, OT_F12, info.parameters);
    auto f12ptr = f12_cc->get_op();

    // GFit of 1/r and partition into three slabs by α.
    GFit<double, LDIM> fit = GFit<double, LDIM>::CoulombFit(
            opt.lo, opt.hi, opt.eps_gfit, false);
    const Tensor<double> c_all = fit.coeffs();
    const Tensor<double> a_all = fit.exponents();
    const long M = c_all.size();
    t.tag("fit operator expansion");

    std::vector<double> cs_diff, as_diff, cs_med, as_med;
    long n_tight = 0;
    for (long mu = 0; mu < M; ++mu) {
        const double a = a_all[mu];
        const double c = c_all[mu];
        if (a > opt.alpha_hi) {            // tight: discard
            ++n_tight;
            continue;
        }
        MADNESS_ASSERT(c > 0.0);
        if (a < opt.alpha_lo) {            // diffuse
            cs_diff.push_back(c);
            as_diff.push_back(a);
        } else {                           // medium
            cs_med.push_back(c);
            as_med.push_back(a);
        }
    }
    const long Md = cs_diff.size();
    const long Mm = cs_med.size();
    if (world.rank() == 0) {
        print("[three-range] GFit partition: M_total =", M,
              " diffuse =", Md, " medium =", Mm, " tight (discard) =", n_tight,
              " | alpha_lo =", opt.alpha_lo, " alpha_hi =", opt.alpha_hi);
    }
    MADNESS_CHECK_THROW(Md > 0 || Mm > 0,
                       "three-range: both diffuse and medium are empty");

    auto make_op = [&](const std::vector<double>& cs,
                       const std::vector<double>& as)
            -> std::shared_ptr<SeparatedConvolution<double, LDIM>> {
        const long Mk = cs.size();
        if (Mk == 0) return {};
        Tensor<double> c_t(Mk), a_t(Mk);
        for (long mu = 0; mu < Mk; ++mu) { c_t[mu] = cs[mu]; a_t[mu] = as[mu]; }
        return std::make_shared<SeparatedConvolution<double, LDIM>>(
                world, c_t, a_t, opt.lo, thresh);
    };
    auto op_diff = make_op(cs_diff, as_diff);
    auto op_med  = make_op(cs_med,  as_med);

    std::vector<Vector<double, LDIM>> origins = info.molecular_coordinates;

    // Build the two LRFs with their own parameter sets.  Each one represents
    // its sub-kernel · k(r₁) · k(r₂); the kernel splits linearly so the two
    // (g, h) sets simply concatenate into the joint decomposition used below.
    std::vector<Function<double, LDIM>> gvec, hvec;
    long rank_diff = 0, rank_med = 0;

    if (op_diff) {
        LRFunctorF12<double, NDIM> functor_diff(op_diff, info.mo_ket, info.mo_bra);
        auto lrf_diff = LowRankFunctionFactory<double, NDIM>(
                            opt.lrfparam_diffuse, origins)
                .project(functor_diff, FunctionDefaults<6>::get_thresh(), 0);
        rank_diff = lrf_diff.g.size();
        gvec.insert(gvec.end(), lrf_diff.g.begin(), lrf_diff.g.end());
        hvec.insert(hvec.end(), lrf_diff.h.begin(), lrf_diff.h.end());
        t.tag("construct LRF (diffuse)");
    }
    if (op_med && !opt.medium_use_6d) {
        LRFunctorF12<double, NDIM> functor_med(op_med, info.mo_ket, info.mo_bra);
        auto factory_med = LowRankFunctionFactory<double, NDIM>(
                                opt.lrfparam_medium, origins);
        auto lrf_med = opt.medium_use_taylor
            ? factory_med.project_from_operator(
                    functor_med,
                    FunctionDefaults<6>::get_thresh(),
                    opt.max_taylor_order)
            : factory_med.project(
                    functor_med,
                    FunctionDefaults<6>::get_thresh(), 0);
        rank_med = lrf_med.g.size();
        gvec.insert(gvec.end(), lrf_med.g.begin(), lrf_med.g.end());
        hvec.insert(hvec.end(), lrf_med.h.begin(), lrf_med.h.end());
        t.tag(opt.medium_use_taylor
              ? "construct LRF (medium, Taylor)"
              : "construct LRF (medium, random-Y)");
    }
    if (world.rank() == 0) {
        print("[three-range] LRF ranks: diffuse =", rank_diff,
              " medium =", rank_med, " total =", rank_diff + rank_med);
    }

    double tol = FunctionDefaults<LDIM>::get_thresh();

    // piece 1 (K̂₁ f part):  Σ_ρ g_ρ(1) · (phi_j · A_ρ)(2)
    // A_ρ(r₂) = f12(h_ρ · phi_i)(r₂)
    auto j_A_pi = phi_j.function * apply(world, *f12ptr, hvec * phi_i.function);
    auto Kf = LowRankFunction<double, NDIM>(gvec, j_A_pi, tol);
    if (symmetric) {
        Kf += LowRankFunction<double, NDIM>(j_A_pi, gvec, tol);
    } else {
        auto i_A_pj = phi_i.function * apply(world, *f12ptr, hvec * phi_j.function);
        Kf += LowRankFunction<double, NDIM>(i_A_pj, gvec, tol);
    }
    t.tag("construct Kf |ij>");

    // piece 2 (f K̂₁ part):  Σ_ρ B_ρi · g_ρ(1) · f(1,2) · phi_j(2)
    auto B_pi = inner(world, phi_i.function, hvec);
    auto C_i = transform(world, gvec, B_pi);
    MADNESS_CHECK_THROW(C_i.size() == 1, "invalid size for Ci");
    auto fK = LowRankFunction<double, 6>(C_i, {phi_j.function}, tol);
    if (symmetric) {
        fK += LowRankFunction<double, NDIM>({phi_j.function}, C_i, tol);
    } else {
        auto B_pj = inner(world, phi_j.function, hvec);
        auto C_j  = transform(world, gvec, B_pj);
        fK += LowRankFunction<double, NDIM>({phi_i.function}, C_j, tol);
    }
    t.tag("construct fK |ij>");

    out.Kf = { CCPairFunction<double, 6>(Kf.get_g(), Kf.get_h()) };
    out.fK = { CCPairFunction<double, 6>(f12_cc, fK.get_g(), fK.get_h()) };
    out.fK[0].convert_to_pure_no_op_inplace();
    t.tag("convert fK to pure");

    out.KffK.push_back(out.Kf[0]);
    out.KffK.push_back(-1.0 * out.fK[0]);

    // -----------------------------------------------------------------
    // 6D medium piece: only built when opt.medium_use_6d is set.  Uses
    // the partial-Coulomb operator op_med (medium Gaussians only) and
    // the standard 6D apply_Kfxy / K_macrotask machinery.  Result is a
    // pure 6D pair appended to Kf / fK / KffK.  Tight is still
    // discarded (handled above by the GFit partition).
    // -----------------------------------------------------------------
    if (op_med && opt.medium_use_6d) {
        // K̂_med f12 |ij⟩  ------------------------------------------------
        real_function_6d Kfxy_med = apply_Kfxy_partial(
                world, phi_i, phi_j, info, info.parameters, *op_med);
        Kfxy_med.truncate().reduce_rank();
        t.tag("construct Kf medium (6D)");

        // f12 K̂_med |ij⟩  -----------------------------------------------
        const real_function_3d x_ket = phi_i.function;
        const real_function_3d y_ket = phi_j.function;
        const real_function_3d x_bra = (info.R_square * phi_i.function).truncate();
        const real_function_3d y_bra = (info.R_square * phi_j.function).truncate();

        const real_function_3d Kx_med = K_partial_3d(
                world, info.mo_ket, info.mo_bra, x_ket, *op_med);
        const FuncType Kx_type = UNDEFINED;
        const bool symmetric_fk = (phi_i == phi_j);
        const real_function_6d fKphi0b = CCPotentials::make_f_xy_macrotask(
                world, Kx_med, y_ket, x_bra, y_bra, phi_i.i, phi_j.i,
                info.parameters, Kx_type, phi_j.type, /*Gscreen=*/nullptr);
        real_function_6d fKphi0a;
        if (symmetric_fk) {
            fKphi0a = madness::swap_particles(fKphi0b);
        } else {
            const real_function_3d Ky_med = K_partial_3d(
                    world, info.mo_ket, info.mo_bra, y_ket, *op_med);
            const FuncType Ky_type = UNDEFINED;
            fKphi0a = CCPotentials::make_f_xy_macrotask(
                    world, x_ket, Ky_med, x_bra, y_bra, phi_i.i, phi_j.i,
                    info.parameters, phi_i.type, Ky_type, /*Gscreen=*/nullptr);
        }
        real_function_6d fKxy_med = (fKphi0a + fKphi0b);
        fKxy_med.truncate().reduce_rank();
        t.tag("construct fK medium (6D)");

        out.Kf.push_back(CCPairFunction<double, 6>(Kfxy_med));
        out.fK.push_back(CCPairFunction<double, 6>(fKxy_med));
        out.KffK.push_back(CCPairFunction<double, 6>(Kfxy_med));
        out.KffK.push_back(-1.0 * CCPairFunction<double, 6>(fKxy_med));

        if (world.rank() == 0) {
            print("[three-range] medium=6D piece appended:",
                  " Kf size =",  get_size(Kfxy_med), "GB",
                  " fK size =",  get_size(fKxy_med), "GB");
        }
    }

    finalize_sizes(world, out);
    out.rank   = rank_diff + rank_med;
    out.t_wall = t1.elapsed();
    return out;
}

// ---------------------------------------------------------------------------
//  Diagnostics (extracted from benchmark's test_kcomm_accuracy)
// ---------------------------------------------------------------------------

ExchangeCommutator::Diagnostics
ExchangeCommutator::diagnose(
        World& world,
        const std::vector<Function<double,3>>& kvec,
        const std::vector<Function<double,3>>& R2kvec,
        const Function<double,3>& phi_i,
        const Function<double,3>& phi_j,
        const std::vector<CCPairFunction<double,6>>& Kf,
        const std::vector<CCPairFunction<double,6>>& fK,
        const std::vector<CCPairFunction<double,6>>& KffK,
        const LowRankFunctionParameters& obs_param,
        bool verbose,
        bool include_K2,
        const std::vector<Vector<double,3>>& centers_in) const
{
    Diagnostics d;
    wall_timer t(world);

    // Observer basis: prefer the AO basis stored on this instance.  When
    // empty, fall back to a canonical-orthonormalized harmonic basis built
    // from obs_param at centers_in (origin if empty).  Use *canonical*
    // orthonormalization with a linear-dependency cutoff: the iterative
    // Löwdin in plain orthonormalize() diverges to NaN when the raw
    // Cartesian-Gaussian set has near-zero overlap eigenvalues (which
    // happens routinely once multiple centers are stacked).
    std::vector<Function<double,3>> phi_a_owned;
    if (ao_basis.empty()) {
        const auto centers = centers_in.empty()
                ? std::vector<Vector<double,3>>({ Vector<double,3>(0.0) })
                : centers_in;
        phi_a_owned = LowRankFunctionFactory<double,6>::harmonic_basis(
                world, obs_param.tempered(), 2, centers);
        phi_a_owned = orthonormalize_canonical(phi_a_owned);
    }
    const auto& phi_a = ao_basis.empty() ? phi_a_owned : ao_basis;

    const bool symmetric_ij = (&phi_i == &phi_j);  // cheap fast path; full check below

    if (verbose) {
        const char* basis_kind = ao_basis.empty() ? "harmonic (fallback)" : "AO";
        print("[diagnose]", basis_kind, "observer basis: phi_a.size() =", phi_a.size(),
              " kvec.size() =", kvec.size(),
              " ||phi_i|| =", phi_i.norm2(),
              " ||phi_j|| =", phi_j.norm2(),
              " same-object i,j =", symmetric_ij);
        std::vector<double> norms_a(phi_a.size());
        for (std::size_t a = 0; a < phi_a.size(); ++a) norms_a[a] = phi_a[a].norm2();
        print("[diagnose] ||phi_a[a]||:", norms_a);
    }

    // K̂ in nemo formalism is non-self-adjoint: bra-side R² makes
    //   K̂  = set_bra_and_ket(R²k, k)   acts on a ket: K̂φ(r) = Σ_k k(r) ∫ R²k(r') φ(r') dr'
    //   K̂† = set_bra_and_ket(k, R²k)   acts on a bra: K̂†φ(r) = R²(r) Σ_k k(r) ∫ k(r') φ(r') dr'
    // The algorithm (apply_KffK_lowrank etc.) uses info.mo_bra = R²·k and
    // info.mo_ket = k internally, so the reference must respect this:
    //   * apply K̂  on phi_i / phi_j (ket side)
    //   * apply K̂† on phi_a         (bra side, when moved over from ⟨ab| via Hermiticity)
    madness::Exchange<double, 3> K(world, 1.e-6);
    K.set_bra_and_ket(R2kvec, kvec);
    madness::Exchange<double, 3> Kdagger(world, 1.e-6);
    Kdagger.set_bra_and_ket(kvec, R2kvec);

    auto f12ptr = std::shared_ptr<SeparatedConvolution<double,3>>(
            SlaterF12OperatorPtr_ND<3>(world, 1.0, 1.e-6,
                                       FunctionDefaults<3>::get_thresh()));
    auto& f12 = *f12ptr;

    // some helpful intermediates
    const auto f12_aj = f12(phi_a * phi_j);
    const auto f12_ai = f12(phi_a * phi_i);

    auto Kdagger_a = Kdagger(phi_a);
    auto K_i = K(phi_i);
    auto K_j = K(phi_j);

    if (verbose) {
        std::vector<double> nKda(Kdagger_a.size());
        for (std::size_t a = 0; a < Kdagger_a.size(); ++a) nKda[a] = Kdagger_a[a].norm2();
        std::vector<double> nfaj(f12_aj.size());
        for (std::size_t a = 0; a < f12_aj.size(); ++a) nfaj[a] = f12_aj[a].norm2();
        print("[diagnose] ||K_dagger(phi_a)||:", nKda);
        print("[diagnose] ||f12(phi_a*phi_j)||:", nfaj);
        print("[diagnose] ||K(phi_i)|| =", K(phi_i).norm2(),
              " ||K(phi_j)|| =", K(phi_j).norm2());
    }

    // Reference integrals (analytic ⟨ab | · | ij⟩) — computed per K̂_p so the
    // four pieces are individually inspectable.
    //   ref_Kf_K1[a,b] = ⟨ab | K̂₁ f | ij⟩  =  ⟨(K̂†a)·i | f₁₂ | b·j⟩
    //   ref_Kf_K2[a,b] = ⟨ab | K̂₂ f | ij⟩  =  ⟨a·i        | f₁₂ | (K̂†b)·j⟩
    //   ref_fK_K1[a,b] = ⟨ab | f K̂₁ | ij⟩  =  ⟨a·(Ki)     | f₁₂ | b·j⟩
    //   ref_fK_K2[a,b] = ⟨ab | f K̂₂ | ij⟩  =  ⟨a·i        | f₁₂ | b·(Kj)⟩
    // K̂ on the ket gets K(phi_{i/j}); moved over to the bra via Hermiticity
    // of the integral it becomes K̂†(phi_a).
    //
    // Symmetry note: when phi_i = phi_j the K̂₁ and K̂₂ matrices are
    // *transposes* of each other (ref_Kf_K2[a,b] = ref_Kf_K1[b,a]), not
    // identical — the harmonic-basis matrix is not (a,b)-symmetric because
    // K̂† carries the asymmetric R²-weighted bra.  Their Frobenius norms
    // coincide, but ‖K1 − K2‖ measures 2·‖antisymmetric part‖, not zero.
    // Equality K1 = K2 is only expected as a consistency check between
    // independent algorithm paths, not a property of either reference alone.
    d.ref_Kf_K1 = matrix_inner(world, Kdagger_a * phi_i, f12_aj);
    d.ref_Kf_K2 = matrix_inner(world, f12_ai, Kdagger_a * phi_j);
    d.ref_fK_K1 = matrix_inner(world, phi_a * K_i,      f12_aj);
    d.ref_fK_K2 = matrix_inner(world, f12_ai,           phi_a * K_j);

    d.ref_piece1 = copy(d.ref_Kf_K1);
    d.ref_piece2 = copy(d.ref_fK_K1);
    if (include_K2) {
        d.ref_piece1 += d.ref_Kf_K2;
        d.ref_piece2 += d.ref_fK_K2;
    }

    if (verbose) {
        // K̂₁ vs K̂₂ matrices are transposes when phi_i = phi_j (not equal);
        // their norms coincide but ‖K1 − K2‖ reports 2·‖antisymmetric part‖.
        const double sym_Kf = (d.ref_Kf_K1 - transpose(d.ref_Kf_K2)).normf();
        const double sym_fK = (d.ref_fK_K1 - transpose(d.ref_fK_K2)).normf();
        print("[diagnose] ||ref_Kf_K1|| =", d.ref_Kf_K1.normf(),
              " ||ref_Kf_K2|| =", d.ref_Kf_K2.normf(),
              " ||ref_Kf_K1 - ref_Kf_K2^T|| =", sym_Kf,
              "   (zero iff phi_i = phi_j)");
        print("[diagnose] ||ref_fK_K1|| =", d.ref_fK_K1.normf(),
              " ||ref_fK_K2|| =", d.ref_fK_K2.normf(),
              " ||ref_fK_K1 - ref_fK_K2^T|| =", sym_fK,
              "   (zero iff phi_i = phi_j)");
        print("[diagnose] ||ref_piece1|| =", d.ref_piece1.normf(),
              " ||ref_piece2|| =", d.ref_piece2.normf(),
              " ||ref_piece1 - ref_piece2|| =", (d.ref_piece1 - d.ref_piece2).normf());
    }

    // Describe one pair entry: kind, rank/size, factor norms.
    auto dump_pair = [&world](const std::string& tag,
                              const std::vector<CCPairFunction<double,6>>& v) {
        if (v.empty()) { print("[diagnose]", tag, ": empty"); return; }
        for (std::size_t k = 0; k < v.size(); ++k) {
            const auto& f = v[k];
            std::string kind = "unassigned";
            if (f.is_assigned()) {
                if (f.is_pure())                     kind = "pure-6d";
                else if (f.is_decomposed_no_op())    kind = "decomposed";
                else if (f.is_op_decomposed())       kind = "op-decomposed";
                else                                 kind = "mixed";
            }
            std::stringstream ss;
            ss << "[diagnose] " << tag << "[" << k << "] kind=" << kind;
            if (f.is_assigned() && (f.is_decomposed_no_op() || f.is_op_decomposed())) {
                auto a = f.get_a(); auto b = f.get_b();
                ss << " rank=" << a.size()
                   << " ||a||=" << get_size(world, a)
                   << " GB ||b||=" << get_size(world, b) << " GB";
            } else if (f.is_assigned() && f.is_pure()) {
                ss << " ||f||=" << f.get_function().norm2();
            }
            print(ss.str());
        }
    };

    if (verbose) {
        dump_pair("Kf",   Kf);
        dump_pair("fK",   fK);
        dump_pair("KffK", KffK);
    }

    // Project each pair entry onto the harmonic-basis bras via partial_inner
    // (integrate particle 1 against phi_a[a], then matrix_inner with phi_a
    // over particle 2).  Sum contributions from all entries in v.
    auto compute_error = [&phi_a, &world, verbose](
            const std::vector<CCPairFunction<double,6>>& v,
            const Tensor<double>& reference,
            const std::string& tag) {
        if (v.empty()) return -1.0;
        auto p1 = particle<3>::particle1();
        Tensor<double> result(phi_a.size(), phi_a.size());
        for (const auto& f : v) {
            if (!f.is_assigned()) continue;
            std::vector<Function<double,3>> tmp(phi_a.size());
            for (std::size_t a = 0; a < phi_a.size(); ++a) {
                tmp[a] = inner(f, phi_a[a], p1.get_tuple(), p1.get_tuple());
            }
            result += matrix_inner(world, tmp, phi_a);
        }
        if (verbose) {
            print("[diagnose]", tag, "||computed||=", result.normf(),
                  " ||reference||=", reference.normf(),
                  " ||diff||=", (reference - result).normf());
        }
        return (reference - result).normf();
    };

    if (!Kf.empty())   d.err_Kf   = compute_error(Kf,   d.ref_piece1,                     "Kf");
    if (!fK.empty())   d.err_fK   = compute_error(fK,   d.ref_piece2,                     "fK");
    if (!KffK.empty()) d.err_KffK = compute_error(KffK, d.ref_piece1 - d.ref_piece2,      "KffK");

    d.t_reference = t.elapsed();
    return d;
}

// ---------------------------------------------------------------------------
//  Uniform one-line reporting
// ---------------------------------------------------------------------------

void ExchangeCommutator::print_report(
        const KffKResult& result,
        const Diagnostics* diag)
{
    char buf[256];
    std::snprintf(buf, sizeof(buf),
                  "[Kcomm] algo=%-18s rank=%5ld  mem=%7.3f GB  wall=%8.2f s",
                  result.algo.c_str(), result.rank, result.mem_gb, result.t_wall);
    print(buf);
    if (diag) {
        std::snprintf(buf, sizeof(buf),
                      "[Kcomm]   err_Kf=%10.3e  err_fK=%10.3e  err_KffK=%10.3e  t_ref=%7.2f s",
                      diag->err_Kf, diag->err_fK, diag->err_KffK, diag->t_reference);
        print(buf);
    }
}

// ---------------------------------------------------------------------------
//  G·[K̂,f] Schwinger diagnostic
// ---------------------------------------------------------------------------

ExchangeCommutator::GKffKDiagnostics
ExchangeCommutator::diagnose_GKffK(
        World& world,
        const std::vector<CCPairFunction<double,6>>& GKf_cc,
        const std::vector<CCPairFunction<double,6>>& GfK_cc,
        const real_function_3d& phi_i,
        const real_function_3d& phi_j,
        const real_function_3d& Kphi_i,
        const real_function_3d& Kphi_j,
        const Info& info,
        double energy) const
{
    MADNESS_CHECK_THROW(energy < 0.0,    "diagnose_GKffK: energy must be negative");
    MADNESS_CHECK_THROW(!ao_basis.empty(),"diagnose_GKffK: ao_basis must be non-empty");

    const double thresh = FunctionDefaults<3>::get_thresh();
    const int nbasis = ao_basis.size();
    const double wall0 = wall_time();

    GKffKDiagnostics diag;
    diag.ref_GKf      = Tensor<double>(nbasis, nbasis);
    diag.ref_GfK      = Tensor<double>(nbasis, nbasis);
    diag.result_GKf   = Tensor<double>(nbasis, nbasis);
    diag.result_GfK   = Tensor<double>(nbasis, nbasis);

    // --- 6D reference: project G·Kf and G·fK onto ⟨ab| -------------------------
    // Uses the same partial_inner / matrix_inner machinery as diagnose()::compute_error,
    // which handles pure, decomposed, and op-decomposed CCPairFunctions uniformly.
    // Multiple entries in the vector are summed (mirrors the multi-piece KffKResult).
    auto project_cc = [&](const std::vector<CCPairFunction<double,6>>& v) {
        auto p1 = particle<3>::particle1();
        Tensor<double> result(nbasis, nbasis);
        for (const auto& f : v) {
            if (!f.is_assigned()) continue;

            // matrix elements for U
            for (int a = 0; a < ao_basis.size(); ++a) {
                for (int b = 0; b < ao_basis.size(); ++b) {
                    CCPairFunction<double,6> ab(ao_basis[a],ao_basis[b]);
                    result(a,b)+= inner(ab,f);
                }
            }
        }
        return result;
    };
    diag.ref_GKf   = project_cc(GKf_cc);
    diag.ref_GfK   = project_cc(GfK_cc);
    diag.ref_GKffK = diag.ref_GKf - diag.ref_GfK;

    // --- Schwinger quadrature ---------------------------------------------------
    // G = (T−E)⁻¹ ≈ Σₙ w_n^{6d} [g_{αₙ}*·] ⊗ [g_{αₙ}*·]
    //
    // 6D weight: w_n^{6d} = c_n^{bsh} · (αₙ/π)^{3/2}
    //   (6D heat-kernel Jacobian vs 3D BSH fit weight, see lrf_G_Ue_diagnostics.md)
    //
    // Kf piece — move K̂ to the bra via its adjoint K̂†:
    //   ⟨ã b̃ | K̂₁ f₁₂ | φᵢ φⱼ⟩ = inner((K̂†ã)·φᵢ, f₁₂(b̃·φⱼ))    ← particle 1
    //   ⟨ã b̃ | K̂₂ f₁₂ | φᵢ φⱼ⟩ = inner((K̂†b̃)·φⱼ, f₁₂(ã·φᵢ))    ← particle 2
    //
    // fK piece — K̂ acts on the ket orbital before f₁₂:
    //   ⟨ã b̃ | f₁₂ K̂₁ | φᵢ φⱼ⟩ = inner(ã·K̂φᵢ, f₁₂(b̃·φⱼ))        ← particle 1
    //   ⟨ã b̃ | f₁₂ K̂₂ | φᵢ φⱼ⟩ = inner(ã·φᵢ,   f₁₂(b̃·K̂φⱼ))      ← particle 2

    const double mu    = std::sqrt(-2.0 * energy);
    const double lo    = info.parameters.lo();
    const double hi    = FunctionDefaults<3>::get_cell_width().normf();
    const double gamma = info.parameters.gamma();

    auto fit   = GFit<double,3>::BSHFit(mu, lo, hi, thresh);
    auto c3d   = fit.coeffs();
    auto alpha = fit.exponents();
    const int nfit = c3d.dim(0);

    // K̂ (acts on ket) and K̂† (acts on bra), matching the nemo convention in diagnose():
    //   K̂  = set_bra_and_ket(R²·k, k):  K̂φ  = Σ_k k_ket  · Coulomb(k_bra·φ)
    //   K̂† = set_bra_and_ket(k, R²·k):  K̂†φ = Σ_k k_bra  · Coulomb(k_ket·φ)
    madness::Exchange<double,3> Kdagger(world, lo);
    Kdagger.set_bra_and_ket(info.mo_ket, info.mo_bra);

    // Slater f₁₂ operator — same normalization as in diagnose()
    auto f12ptr = std::shared_ptr<SeparatedConvolution<double,3>>(
            SlaterF12OperatorPtr_ND<3>(world, gamma, lo, thresh));
    auto& f12 = *f12ptr;

    for (int n = 0; n < nfit; ++n) {
        const double an  = alpha[n];
        const double w6d = c3d[n] * std::pow(an / constants::pi, 1.5);

        // Convolve each AO with exp(−αₙ r²) to get ã_n / b̃_n
        auto gauss = SeparatedConvolution<double,3>(world, OperatorInfo(an, lo, thresh, OT_GAUSS));
        std::vector<real_function_3d> conv(nbasis);
        for (int a = 0; a < nbasis; ++a)
            conv[a] = gauss(ao_basis[a]);

        // K̂†(ã_n) for all a — needed for the Kf piece
        auto Kdagger_conv = Kdagger(conv);

        // Precompute products with fixed orbitals (vectorised * broadcasts the scalar function)
        auto Kdagger_conv_phi_i = Kdagger_conv * phi_i;  // (K̂†ã)·φᵢ   [size nbasis]
        auto Kdagger_conv_phi_j = Kdagger_conv * phi_j;  // (K̂†ã)·φⱼ
        auto conv_Kphi_i        = conv * Kphi_i;         // ã·K̂φᵢ
        auto conv_phi_i         = conv * phi_i;          // ã·φᵢ

        // Precompute f₁₂ applications indexed by b (the particle-2 basis index)
        auto f12_conv_phi_j  = f12(conv * phi_j);   // f₁₂(b̃·φⱼ)   used by both pieces
        auto f12_conv_phi_i  = f12(conv * phi_i);   // f₁₂(b̃·φᵢ)   Kf K̂₂
        auto f12_conv_Kphi_j = f12(conv * Kphi_j);  // f₁₂(b̃·K̂φⱼ) fK K̂₂

        // Kf accumulation: two matrix_inner calls (K̂₁ and K̂₂)
        //   Kf K̂₁[a,b] = inner((K̂†ã_a)·φᵢ,   f₁₂(b̃_b·φⱼ))
        //   Kf K̂₂[a,b] = inner(f₁₂(ã_a·φᵢ),   (K̂†b̃_b)·φⱼ)
        //             = inner((K̂†b̃_b)·φⱼ, f₁₂(ã_a·φᵢ))  → rows/cols swapped → use transpose
        diag.result_GKf += w6d * matrix_inner(world, Kdagger_conv_phi_i, f12_conv_phi_j);
        diag.result_GKf += w6d * transpose(matrix_inner(world, Kdagger_conv_phi_j, f12_conv_phi_i));

        // fK accumulation: two matrix_inner calls (K̂₁ and K̂₂)
        //   fK K̂₁[a,b] = inner(ã_a·K̂φᵢ, f₁₂(b̃_b·φⱼ))
        //   fK K̂₂[a,b] = inner(ã_a·φᵢ,   f₁₂(b̃_b·K̂φⱼ))
        diag.result_GfK += w6d * matrix_inner(world, conv_Kphi_i, f12_conv_phi_j);
        diag.result_GfK += w6d * matrix_inner(world, conv_phi_i,  f12_conv_Kphi_j);
    }

    diag.result_GKffK = diag.result_GKf - diag.result_GfK;
    diag.error_GKf    = (diag.ref_GKf   - diag.result_GKf).normf();
    diag.error_GfK    = (diag.ref_GfK   - diag.result_GfK).normf();
    diag.error_GKffK  = (diag.ref_GKffK - diag.result_GKffK).normf();
    diag.time = wall_time() - wall0;
    return diag;
}

} // namespace madness
