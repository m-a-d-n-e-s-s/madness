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
//  LRF path (canonical) — wraps CCPotentials::apply_KffK_low_rank
// ---------------------------------------------------------------------------

ExchangeCommutator::KffKResult
ExchangeCommutator::apply_KffK_lowrank(
        World& world,
        const CCFunction<double,3>& phi_i,
        const CCFunction<double,3>& phi_j,
        const Info& info,
        const LowRankFunctionParameters& lrfparam,
        const real_convolution_6d* Gscreen)
{
    KffKResult out;
    out.algo = "lrf";
    wall_timer t(world);

    // CCPotentials::apply_KffK_low_rank returns { Kf, -fK }.
    auto result = CCPotentials::apply_KffK_low_rank(
            world, phi_i, phi_j, info, Gscreen, lrfparam);
    if (result.size() >= 1) out.Kf = { result[0] };
    if (result.size() >= 2) out.fK = { -1.0 * result[1] };   // undo the sign
    out.KffK = std::move(result);

    finalize_sizes(world, out);
    out.t_wall = t.elapsed();
    return out;
}

// ---------------------------------------------------------------------------
//  LRF path (per-k direct) — wraps CCPotentials::apply_KffK_low_rank_direct
// ---------------------------------------------------------------------------

ExchangeCommutator::KffKResult
ExchangeCommutator::apply_KffK_lowrank_direct(
        World& world,
        const CCFunction<double,3>& phi_i,
        const CCFunction<double,3>& phi_j,
        const Info& info,
        const LowRankFunctionParameters& lrfparam,
        const real_convolution_6d* Gscreen)
{
    KffKResult out;
    out.algo = "lrf-direct";
    wall_timer t(world);

    // apply_KffK_low_rank_direct builds the commutator as a single combined
    // CCPairFunction; Kf / fK are not produced separately.
    out.KffK = CCPotentials::apply_KffK_low_rank_direct(
            world, phi_i, phi_j, info, Gscreen, lrfparam);

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
    wall_timer t(world);

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

    Vector<double, LDIM> origin(0.0);
    std::vector<Vector<double, LDIM>> origins = { origin };

    // Accumulators for the two commutator pieces — each stored as a single
    // CCPairFunction stack, so they compose cleanly into Kf / fK / KffK.
    //
    //   Kf_g/Kf_h   : separable pair   g_ρ(1) · (phi_j · A_ρ)(2)
    //   fK_g/fK_h   : f12-wrapped pair (B_ρ · g_ρ)(1) · f(1,2) · phi_j(2)
    //
    // The minus sign of the commutator is applied on the KffK stack only.
    std::vector<Function<double, LDIM>> Kf_g, Kf_h;
    std::vector<Function<double, LDIM>> fK_g, fK_h;

    long rank_accum = 0;
    for (std::size_t k = 0; k < info.mo_ket.size(); ++k) {
        const auto& k_bra = info.mo_bra[k];
        const auto& k_ket = info.mo_ket[k];

        // LRF of k_bra(1) · K_smooth(1,1') · k_ket(1') — represents the
        // k-th occupied contribution to the exchange kernel.
        LRFunctorF12<double, NDIM> functor(trunc_op, k_bra, k_ket);
        auto lrf = LowRankFunctionFactory<double, NDIM>(lrfparam, origins)
                   .project(functor, opt.eps_lrf, 0);

        auto& gvec = lrf.g;
        auto& hvec = lrf.h;
        const long R = gvec.size();
        rank_accum += R;

        // A_ρ(r₂) = f12(h_ρ · phi_i)(r₂)
        auto A_vec = apply(world, *f12ptr, hvec * phi_i.function);
        // B_ρ     = ⟨h_ρ | phi_i⟩
        Tensor<double> B = inner(world, phi_i.function, hvec);

        // piece 1 (K̂₁ f part):  Σ_ρ g_ρ(1) · (phi_j · A_ρ)(2)
        std::vector<Function<double, LDIM>> h1 = phi_j.function * A_vec;
        for (long r = 0; r < R; ++r) {
            Kf_g.push_back(gvec[r]);
            Kf_h.push_back(h1[r]);
        }

        // piece 2 (f K̂₁ part):  Σ_ρ B_ρ · g_ρ(1) · f(1,2) · phi_j(2)
        // (phi_i does NOT appear here — it only enters via B_ρ = ⟨h_ρ | i⟩.)
        for (long r = 0; r < R; ++r) {
            Function<double, LDIM> scaled = B(r) * gvec[r];
            fK_g.push_back(scaled);
            fK_h.push_back(phi_j.function);
        }
    }

    // Assemble named pair stacks.  Each piece becomes a single CCPairFunction
    // with multi-column factors, and KffK bakes in the commutator's sign.
    if (!Kf_g.empty()) {
        out.Kf = { CCPairFunction<double, 6>(Kf_g, Kf_h) };
    }
    if (!fK_g.empty()) {
        out.fK = { CCPairFunction<double, 6>(f12_cc, fK_g, fK_h) };
    }
    if (!out.Kf.empty())  out.KffK.push_back(out.Kf[0]);
    if (!out.fK.empty())  out.KffK.push_back(-1.0 * out.fK[0]);

    // Symmetric case (phi_i == phi_j): particle-2 mirror via swap_particles.
    if (symmetric) {
        if (!out.Kf.empty()) {
            auto m = swap_particles(out.Kf);
            for (auto&& cc : m) out.Kf.push_back(std::move(cc));
        }
        if (!out.fK.empty()) {
            auto m = swap_particles(out.fK);
            for (auto&& cc : m) out.fK.push_back(std::move(cc));
        }
        if (!out.KffK.empty()) {
            auto m = swap_particles(out.KffK);
            for (auto&& cc : m) out.KffK.push_back(std::move(cc));
        }
    }

    finalize_sizes(world, out);
    out.rank   = rank_accum;   // LRF rank of the underlying decomposition
    out.t_wall = t.elapsed();
    return out;
}

// ---------------------------------------------------------------------------
//  Diagnostics (extracted from benchmark's test_kcomm_accuracy)
// ---------------------------------------------------------------------------

ExchangeCommutator::Diagnostics
ExchangeCommutator::diagnose(
        World& world,
        const std::vector<Function<double,3>>& kvec,
        const Function<double,3>& phi_i,
        const Function<double,3>& phi_j,
        const std::vector<CCPairFunction<double,6>>& Kf,
        const std::vector<CCPairFunction<double,6>>& fK,
        const std::vector<CCPairFunction<double,6>>& KffK,
        const LowRankFunctionParameters& obs_param)
{
    Diagnostics d;
    wall_timer t(world);

    // Harmonic-basis observer functions (solid harmonics on origin).
    auto centers = std::vector<Vector<double,3>>({ Vector<double,3>(0.0) });
    auto phi_a = LowRankFunctionFactory<double,6>::harmonic_basis(
            world, obs_param.tempered(), 2, centers);
    phi_a = orthonormalize(phi_a);

    madness::Exchange<double, 3> K(world, 1.e-6);
    K.set_bra_and_ket(kvec, kvec);
    auto f12ptr = std::shared_ptr<SeparatedConvolution<double,3>>(
            SlaterF12OperatorPtr_ND<3>(world, 1.0, 1.e-6,
                                       FunctionDefaults<3>::get_thresh()));
    auto& f12 = *f12ptr;

    // Reference integrals.
    // <a b | [K̂₁, f] | i j> = <K(a)·i | f(b·j)> + <f(a·i) | K(a)·j>
    d.ref_piece1  = matrix_inner(world, K(phi_a) * phi_i, f12(phi_a * phi_j));
    d.ref_piece1 += matrix_inner(world, f12(phi_a * phi_i), K(phi_a) * phi_j);
    d.ref_piece2  = matrix_inner(world, phi_a * K(phi_i), f12(phi_a * phi_j));
    d.ref_piece2 += matrix_inner(world, f12(phi_a * phi_i), phi_a * K(phi_j));

    // Sum over pair pieces when projecting onto the harmonic-basis bra.
    auto compute_error = [&phi_a, &world](
            const std::vector<CCPairFunction<double,6>>& v,
            const Tensor<double>& reference) {
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
        return (reference - result).normf();
    };

    if (!Kf.empty())   d.err_Kf   = compute_error(Kf,   d.ref_piece1);
    if (!fK.empty())   d.err_fK   = compute_error(fK,   d.ref_piece2);
    if (!KffK.empty()) d.err_KffK = compute_error(KffK, d.ref_piece1 - d.ref_piece2);

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

} // namespace madness
