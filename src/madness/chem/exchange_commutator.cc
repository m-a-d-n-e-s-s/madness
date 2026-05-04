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

} // namespace madness
