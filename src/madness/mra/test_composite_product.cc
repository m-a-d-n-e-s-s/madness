/*
  This file is part of MADNESS.

  Tests for the unified composite_product operator (composite_product_op_NS).

  Covers, at NDIM ∈ {2, 4, 6}, with LDIM = NDIM/2:
    function mode:
      A) one term, factor on particle 1                       (= multiply(f,g,1))
      B) one term, factor on particle 2                       (= multiply(f,g,2))
      C) one term, factor on each particle                    (= make_Vphi_ij_u: u*c*d)
      D) sum of three terms (v1·u, v2·u, eri·u)               (= make_Vphi)
      E) one term with extra on-demand factor
      H) ket as a sum of hartree pairs (no explicit NDIM ket)
    inner mode:
      F) single term contracted against on-demand bra         (= inner_ij_u_eri)
      G) two-term sum contracted against on-demand bra
*/

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/world/test_utilities.h>
#include <iomanip>

using namespace madness;

// ---------------------------------------------------------------------------
// Templated LDIM-Gaussian building blocks (function pointers usable by Factory)
// ---------------------------------------------------------------------------

template <std::size_t LDIM> double gauss_a (const Vector<double,LDIM>& r) { return std::exp(-1.0 * inner(r,r)); }
template <std::size_t LDIM> double gauss_b (const Vector<double,LDIM>& r) { return std::exp(-1.3 * inner(r,r)); }
template <std::size_t LDIM> double gauss_c (const Vector<double,LDIM>& r) { return std::exp(-0.7 * inner(r,r)); }
template <std::size_t LDIM> double gauss_d (const Vector<double,LDIM>& r) { return std::exp(-0.5 * inner(r,r)); }
template <std::size_t LDIM> double gauss_v1(const Vector<double,LDIM>& r) { return std::exp(-0.3 * inner(r,r)); }
template <std::size_t LDIM> double gauss_v2(const Vector<double,LDIM>& r) { return std::exp(-0.2 * inner(r,r)); }

// ---------------------------------------------------------------------------
// Templated NDIM analytic references (u = g_a(1) ⊗ g_b(2), eri = exp(-β r12²), β=0.5)
// ---------------------------------------------------------------------------

template <std::size_t NDIM>
double exact_u(const Vector<double,NDIM>& r) {
    constexpr std::size_t LDIM = NDIM/2;
    Vector<double,LDIM> r1, r2;
    for (std::size_t i=0; i<LDIM; ++i) { r1[i]=r[i]; r2[i]=r[i+LDIM]; }
    return std::exp(-1.0*inner(r1,r1) - 1.3*inner(r2,r2));
}
template <std::size_t NDIM>
double exact_uc(const Vector<double,NDIM>& r) {
    constexpr std::size_t LDIM = NDIM/2;
    Vector<double,LDIM> r1, r2;
    for (std::size_t i=0; i<LDIM; ++i) { r1[i]=r[i]; r2[i]=r[i+LDIM]; }
    return std::exp(-1.7*inner(r1,r1) - 1.3*inner(r2,r2));
}
template <std::size_t NDIM>
double exact_ud(const Vector<double,NDIM>& r) {
    constexpr std::size_t LDIM = NDIM/2;
    Vector<double,LDIM> r1, r2;
    for (std::size_t i=0; i<LDIM; ++i) { r1[i]=r[i]; r2[i]=r[i+LDIM]; }
    return std::exp(-1.0*inner(r1,r1) - 1.8*inner(r2,r2));
}
template <std::size_t NDIM>
double exact_ucd(const Vector<double,NDIM>& r) {
    constexpr std::size_t LDIM = NDIM/2;
    Vector<double,LDIM> r1, r2;
    for (std::size_t i=0; i<LDIM; ++i) { r1[i]=r[i]; r2[i]=r[i+LDIM]; }
    return std::exp(-1.7*inner(r1,r1) - 1.8*inner(r2,r2));
}
template <std::size_t NDIM>
double exact_Vphi(const Vector<double,NDIM>& r) {           // v1·u + v2·u + eri·u
    constexpr std::size_t LDIM = NDIM/2;
    Vector<double,LDIM> r1, r2;
    for (std::size_t i=0; i<LDIM; ++i) { r1[i]=r[i]; r2[i]=r[i+LDIM]; }
    const double ur  = std::exp(-1.0*inner(r1,r1) - 1.3*inner(r2,r2));
    const double v1  = std::exp(-0.3*inner(r1,r1));
    const double v2  = std::exp(-0.2*inner(r2,r2));
    const double eri = std::exp(-0.5 * inner(r1-r2, r1-r2));
    return ur * (v1 + v2 + eri);
}
template <std::size_t NDIM>
double exact_extra(const Vector<double,NDIM>& r) {          // u · c(1) · eri(1,2)
    constexpr std::size_t LDIM = NDIM/2;
    Vector<double,LDIM> r1, r2;
    for (std::size_t i=0; i<LDIM; ++i) { r1[i]=r[i]; r2[i]=r[i+LDIM]; }
    const double ur  = std::exp(-1.7*inner(r1,r1) - 1.3*inner(r2,r2));
    const double eri = std::exp(-0.5 * inner(r1-r2, r1-r2));
    return ur * eri;
}
template <std::size_t NDIM>
double exact_pair_sum(const Vector<double,NDIM>& r) {       // (a⊗b + c⊗d) · v1(1)
    constexpr std::size_t LDIM = NDIM/2;
    Vector<double,LDIM> r1, r2;
    for (std::size_t i=0; i<LDIM; ++i) { r1[i]=r[i]; r2[i]=r[i+LDIM]; }
    const double t_ab = std::exp(-(1.0+0.3)*inner(r1,r1) - 1.3*inner(r2,r2));
    const double t_cd = std::exp(-(0.7+0.3)*inner(r1,r1) - 0.5*inner(r2,r2));
    return t_ab + t_cd;
}
template <std::size_t NDIM>
double f12_lambda(const Vector<double,NDIM>& r) {           // exp(-β |r1-r2|²), β = 0.5
    constexpr std::size_t LDIM = NDIM/2;
    Vector<double,LDIM> r1, r2;
    for (std::size_t i=0; i<LDIM; ++i) { r1[i]=r[i]; r2[i]=r[i+LDIM]; }
    return std::exp(-0.5 * inner(r1-r2, r1-r2));
}

// ---------------------------------------------------------------------------
// One-NDIM test driver: builds inputs, runs the 8 cases, returns failure count
// ---------------------------------------------------------------------------

template <typename T, std::size_t NDIM>
int run_tests(World& world, const long k, const double thresh) {
    constexpr std::size_t LDIM = NDIM / 2;
    static_assert(NDIM == 2*LDIM, "NDIM must be even");

    FunctionDefaults<LDIM>::set_k(k);
    FunctionDefaults<NDIM>::set_k(k);
    FunctionDefaults<LDIM>::set_thresh(thresh*0.01);
    FunctionDefaults<NDIM>::set_thresh(thresh);
    FunctionDefaults<LDIM>::set_tensor_type(TT_FULL);
    FunctionDefaults<NDIM>::set_tensor_type(TT_2D);
    FunctionDefaults<LDIM>::set_cubic_cell(-8.0, 8.0);
    FunctionDefaults<NDIM>::set_cubic_cell(-8.0, 8.0);

    test_output t1("composite_product_op_NS: NDIM=" + std::to_string(NDIM));
    t1.set_cout_to_terminal();
    std::cout << std::setprecision(8) << std::scientific;

    // --- inputs ----------------------------------------------------------
    Function<T,LDIM> a  = FunctionFactory<T,LDIM>(world).f(gauss_a<LDIM>);
    Function<T,LDIM> b  = FunctionFactory<T,LDIM>(world).f(gauss_b<LDIM>);
    Function<T,LDIM> c  = FunctionFactory<T,LDIM>(world).f(gauss_c<LDIM>);
    Function<T,LDIM> d  = FunctionFactory<T,LDIM>(world).f(gauss_d<LDIM>);
    Function<T,LDIM> v1 = FunctionFactory<T,LDIM>(world).f(gauss_v1<LDIM>);
    Function<T,LDIM> v2 = FunctionFactory<T,LDIM>(world).f(gauss_v2<LDIM>);

    Function<T,NDIM> u   = hartree_product(a, b);
    Function<T,NDIM> eri = FunctionFactory<T,NDIM>(world).is_on_demand().f(f12_lambda<NDIM>);

    auto make_term = [](Function<T,NDIM> ket) {
        composite_term<T,NDIM> t;
        t.ket = ket;
        return t;
    };

    // Function-mode tolerance is loosened: Function::err() returns pointwise (L_inf) error,
    // which is typically a few × the L_2-bounded `thresh` and grows with NDIM.  Inner-mode
    // contractions integrate the error away and remain very tight.
    const double tol_fn    = thresh * 10.0;   // function-mode checkpoint tolerance
    const double tol_inner = thresh;          // inner-mode checkpoint tolerance

    // Per-case (new, old) timing + error table; printed at the end.
    struct Row { std::string name; double t_new, e_new, t_old, e_old; };
    std::vector<Row> rows;
    auto wt = []() { return wall_time(); };

    // Old reference implementations (multiply, CompositeFactory.fill_tree, inner_ij_u_eri)
    // rely on TT_2D — degenerate for NDIM=2 where LDIM=NDIM/2, so skip the side-by-side
    // there and only show the new column.
    constexpr bool compare_old = (NDIM > 2);

    // --- A: u * c(1)  vs  multiply(u, c, 1) -------------------------------
    {
        double t0 = wt();
        std::vector<composite_term<T,NDIM>> terms{make_term(u)}; terms[0].factor1 = c;
        auto r_new = composite_product(world, terms);
        const double t_new = wt() - t0;
        const double e_new = r_new.err(exact_uc<NDIM>);

        double t_old = -1.0, e_old = -1.0;
        if constexpr (compare_old) {
            t0 = wt();
            auto r_old = multiply(copy(u), c, 1);
            t_old = wt() - t0;
            e_old = r_old.err(exact_uc<NDIM>);
        }
        rows.push_back({"A: u * c(1)", t_new, e_new, t_old, e_old});
        t1.checkpoint(e_new, tol_fn, "A: u * c(1)");
    }
    // --- B: u * d(2)  vs  multiply(u, d, 2) -------------------------------
    {
        double t0 = wt();
        std::vector<composite_term<T,NDIM>> terms{make_term(u)}; terms[0].factor2 = d;
        auto r_new = composite_product(world, terms);
        const double t_new = wt() - t0;
        const double e_new = r_new.err(exact_ud<NDIM>);

        double t_old = -1.0, e_old = -1.0;
        if constexpr (compare_old) {
            t0 = wt();
            auto r_old = multiply(copy(u), d, 2);
            t_old = wt() - t0;
            e_old = r_old.err(exact_ud<NDIM>);
        }
        rows.push_back({"B: u * d(2)", t_new, e_new, t_old, e_old});
        t1.checkpoint(e_new, tol_fn, "B: u * d(2)");
    }
    // --- C: u * c(1) * d(2)  vs  multiply(multiply(u,c,1), d, 2) ----------
    {
        double t0 = wt();
        std::vector<composite_term<T,NDIM>> terms{make_term(u)};
        terms[0].factor1 = c;  terms[0].factor2 = d;
        auto r_new = composite_product(world, terms);
        const double t_new = wt() - t0;
        const double e_new = r_new.err(exact_ucd<NDIM>);

        double t_old = -1.0, e_old = -1.0;
        if constexpr (compare_old) {
            t0 = wt();
            auto r_old = multiply(multiply(copy(u), c, 1), d, 2);
            t_old = wt() - t0;
            e_old = r_old.err(exact_ucd<NDIM>);
        }
        rows.push_back({"C: u * c(1) * d(2)", t_new, e_new, t_old, e_old});
        t1.checkpoint(e_new, tol_fn, "C: u * c(1) * d(2)");
    }
    // --- D: v1·u + v2·u + eri·u  vs  CompositeFactory.fill_tree() ---------
    //       The old path lives in FunctionImpl::make_Vphi, invoked through fill_tree
    //       on an on-demand CompositeFunctorInterface.
    {
        double t0 = wt();
        std::vector<composite_term<T,NDIM>> terms;
        auto t_v1c = make_term(u); t_v1c.factor1 = v1; terms.push_back(t_v1c);
        auto t_v2c = make_term(u); t_v2c.factor2 = v2; terms.push_back(t_v2c);
        auto t_erc = make_term(u); t_erc.extra   = eri; terms.push_back(t_erc);
        auto r_new = composite_product(world, terms);
        const double t_new = wt() - t0;
        const double e_new = r_new.err(exact_Vphi<NDIM>);

        double t_old = -1.0, e_old = -1.0;
        if constexpr (compare_old) {
            t0 = wt();
            Function<T,NDIM> r_old = CompositeFactory<T,NDIM,LDIM>(world)
                .ket(copy(u))
                .V_for_particle1(copy(v1))
                .V_for_particle2(copy(v2))
                .g12(eri);
            r_old.fill_tree();
            t_old = wt() - t0;
            e_old = r_old.err(exact_Vphi<NDIM>);
        }
        rows.push_back({"D: v1·u + v2·u + eri·u", t_new, e_new, t_old, e_old});
        t1.checkpoint(e_new, tol_fn, "D: v1·u + v2·u + eri·u");
    }
    // --- E: u * c(1) * eri(1,2)  (no direct old equivalent) ---------------
    {
        double t0 = wt();
        std::vector<composite_term<T,NDIM>> terms{make_term(u)};
        terms[0].factor1 = c;  terms[0].extra = eri;
        auto r_new = composite_product(world, terms);
        const double t_new = wt() - t0;
        const double e_new = r_new.err(exact_extra<NDIM>);

        rows.push_back({"E: u * c(1) * eri(1,2)", t_new, e_new, -1.0, -1.0});
        t1.checkpoint(e_new, tol_fn, "E: u * c(1) * eri(1,2)");
    }
    // --- F: <eri | u·c·d>  vs  inner_ij_u_eri (existing single-term op) ---
    {
        const double beta = 0.5;
        auto conv = GaussOperator<LDIM>(world, beta);
        Function<T,LDIM> ac = a * c;
        Function<T,LDIM> bd = b * d;
        Function<T,LDIM> conv_bd = apply(conv, bd);
        const double ref = inner(ac, conv_bd);

        double t0 = wt();
        std::vector<composite_term<T,NDIM>> terms{make_term(u)};
        terms[0].factor1 = c; terms[0].factor2 = d;
        const double v_new = composite_inner(world, terms, eri);
        const double t_new = wt() - t0;
        const double e_new = std::abs(v_new - ref);

        double t_old = -1.0, e_old = -1.0;
        if constexpr (compare_old) {
            Function<T,NDIM> u_r = copy(u);   u_r.change_tree_state(redundant);
            Function<T,LDIM> c_r = copy(c);   c_r.change_tree_state(redundant);
            Function<T,LDIM> d_r = copy(d);   d_r.change_tree_state(redundant);

            using leaf_opT = Leaf_op<T,NDIM,SeparatedConvolution<double,NDIM>,Specialbox_op<T,NDIM>>;
            leaf_opT lop(u_r.get_impl().get());

            t0 = wt();
            const double v_old = u_r.get_impl()->template inner_ij_u_eri<leaf_opT,LDIM>(
                lop, u_r.get_impl().get(),
                c_r.get_impl().get(), d_r.get_impl().get(),
                eri.get_impl().get());
            t_old = wt() - t0;
            e_old = std::abs(v_old - ref);
        }
        rows.push_back({"F: <eri | u*c*d>", t_new, e_new, t_old, e_old});
        t1.checkpoint(e_new, tol_inner, "F: <eri | u*c*d>");
    }
    // --- G: <eri | u·c·d + u·c>  (no direct old equivalent; new capability) -
    {
        const double beta = 0.5;
        auto conv = GaussOperator<LDIM>(world, beta);
        Function<T,LDIM> ac      = a * c;
        Function<T,LDIM> bd      = b * d;
        Function<T,LDIM> conv_bd = apply(conv, bd);
        Function<T,LDIM> conv_b  = apply(conv, b);
        const double ref = inner(ac, conv_bd) + inner(ac, conv_b);

        double t0 = wt();
        std::vector<composite_term<T,NDIM>> terms;
        auto t_one = make_term(u); t_one.factor1 = c; t_one.factor2 = d; terms.push_back(t_one);
        auto t_two = make_term(u); t_two.factor1 = c;                    terms.push_back(t_two);
        const double v_new = composite_inner(world, terms, eri);
        const double t_new = wt() - t0;
        const double e_new = std::abs(v_new - ref);

        rows.push_back({"G: <eri | u*c*d + u*c>", t_new, e_new, -1.0, -1.0});
        t1.checkpoint(e_new, tol_inner, "G: <eri | u*c*d + u*c>");
    }
    // --- H: (a⊗b + c⊗d) * v1(1)  — pair-sum ket, no explicit NDIM ket -----
    {
        std::vector<composite_term<T,NDIM>> terms;
        composite_term<T,NDIM> term;
        term.p1.push_back(a); term.p2.push_back(b);
        term.p1.push_back(c); term.p2.push_back(d);
        term.factor1 = v1;
        terms.push_back(term);

        auto r = composite_product(world, terms);
        t1.checkpoint(r.err(exact_pair_sum<NDIM>), tol_fn, "H: (a⊗b + c⊗d) * v1(1)");
    }

    // ---- comparison table ----
    if (world.rank() == 0) {
        std::cout << "\n  NDIM=" << NDIM << "  new (composite_product/inner) vs old\n"
                  << "  " << std::string(80, '-') << "\n"
                  << "    case                          t_new[s]   err_new      t_old[s]   err_old\n"
                  << "  " << std::string(80, '-') << "\n";
        for (const auto& r : rows) {
            char buf[256];
            if (r.t_old < 0) {
                std::snprintf(buf, sizeof(buf),
                              "    %-28s  %7.3f   %.2e       —          —",
                              r.name.c_str(), r.t_new, r.e_new);
            } else {
                std::snprintf(buf, sizeof(buf),
                              "    %-28s  %7.3f   %.2e   %7.3f   %.2e",
                              r.name.c_str(), r.t_new, r.e_new, r.t_old, r.e_old);
            }
            std::cout << buf << "\n";
        }
        std::cout << "  " << std::string(80, '-') << "\n";
    }

    return t1.end();
}

// ---------------------------------------------------------------------------

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world, argc, argv);

    const long   k      = 6;
    const double thresh = 1.e-4;

    int rc = 0;
    rc |= run_tests<double, 2>(world, k, thresh);
    rc |= run_tests<double, 4>(world, k, thresh);
    rc |= run_tests<double, 6>(world, k, thresh);

    finalize();
    return rc;
}
