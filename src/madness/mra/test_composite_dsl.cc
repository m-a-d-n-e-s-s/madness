/*
  Tests for the composite_dsl DSL: operator-overload syntax over Function
  expressions, plus the builder API (refine_with / thresh / truncate_mode).

  Coverage at NDIM ∈ {2, 4, 6}, with LDIM = NDIM/2:

  --- Section 1: function-mode / inner-mode DSL ---
    A) u(p1,p2) * c(p1)                              → composite_product
    B) u(p1,p2) * d(p2)
    C) u(p1,p2) * c(p1) * d(p2)
    D) u(p1,p2) * (v(p1) + v(p2) + eri(p1,p2))       → make_Vphi-style sum

    S1) 2.0 * u(p1,p2) * c(p1)
    S2) -u(p1,p2)*c(p1) + 0.5*u(p1,p2)*d(p2)

    O1) apply(BSH, u(p1,p2)) * v(p1)                  → eager apply, then DSL

    I1) inner(eri(p1,p2), u(p1,p2)*c(p1)*d(p2))      → on-demand bra

  --- Section 2: builder API (LeafOpBase / EvalBuilder) ---
    R1) refine_with(default Leaf_op) matches plain .eval()
    R2) refine_with(ElectronCuspyBox_op) on Vphi
    R3) .thresh() reflected on result
    R4) .truncate_mode() reflected on result
    R5) chained refine_with + thresh + truncate_mode
    R6) chaining order is irrelevant
    R7) inner() with refine_with on RHS — value unchanged
    R8) composite_product(LeafOpBase) overload matches default

  Known limitation (NOT exercised here):
    Inner products with a non-on-demand bra (e.g. a hartree pair `i(p1)*j(p2)`)
    currently exhibit a destructor-ordering issue between successive GaussOperator
    LDIM operator constructions; the explicit-bra path itself (composite_inner_apply
    with `cbra_tracker` set) is implemented but not yet exercised end-to-end through
    the DSL.  When that issue is resolved, add:
       inner(i_(p1)*j_(p2), v(p1)*u(p1,p2))
       inner(i_(p1)*j_(p2), (eri(p1,p2) + v(p1)) * u(p1,p2))
*/

#include <madness/mra/mra.h>
#include <madness/mra/composite_dsl.h>
#include <madness/mra/leafop.h>
#include <madness/mra/operator.h>
#include <madness/world/test_utilities.h>
#include <iomanip>
#include <type_traits>

using namespace madness;

namespace {

// ---- LDIM-templated 1-particle leaves ---------------------------------------------
template <std::size_t LDIM> double gauss_a(const Vector<double,LDIM>& r) { return std::exp(-1.0 * inner(r,r)); }
template <std::size_t LDIM> double gauss_b(const Vector<double,LDIM>& r) { return std::exp(-1.3 * inner(r,r)); }
template <std::size_t LDIM> double gauss_c(const Vector<double,LDIM>& r) { return std::exp(-0.7 * inner(r,r)); }
template <std::size_t LDIM> double gauss_d(const Vector<double,LDIM>& r) { return std::exp(-0.5 * inner(r,r)); }
template <std::size_t LDIM> double gauss_v(const Vector<double,LDIM>& r) { return std::exp(-0.3 * inner(r,r)); }

// ---- NDIM-templated analytic references -------------------------------------------
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
double exact_Vphi(const Vector<double,NDIM>& r) {
    constexpr std::size_t LDIM = NDIM/2;
    Vector<double,LDIM> r1, r2;
    for (std::size_t i=0; i<LDIM; ++i) { r1[i]=r[i]; r2[i]=r[i+LDIM]; }
    const double ur  = std::exp(-1.0*inner(r1,r1) - 1.3*inner(r2,r2));
    const double v1  = std::exp(-0.3*inner(r1,r1));
    const double v2_ = std::exp(-0.3*inner(r2,r2));
    const double eri = std::exp(-0.5 * inner(r1-r2, r1-r2));
    return ur * (v1 + v2_ + eri);
}
template <std::size_t NDIM>
double f12_lambda(const Vector<double,NDIM>& r) {
    constexpr std::size_t LDIM = NDIM/2;
    Vector<double,LDIM> r1, r2;
    for (std::size_t i=0; i<LDIM; ++i) { r1[i]=r[i]; r2[i]=r[i+LDIM]; }
    return std::exp(-0.5 * inner(r1-r2, r1-r2));
}

// ---- compile-time ordering check (NDIM-parameterised) ------------------------------
// Term::refine_with(LeafOpBase&) MUST return EvalBuilder, never Function — so the
// builder phase is observable to the type system and `.eval().refine_with(...)` is
// not a valid path.
template <typename T, std::size_t N>
constexpr bool term_refine_returns_builder = std::is_same_v<
    decltype(std::declval<Term<T,N>>().refine_with(std::declval<LeafOpBase<T,N>&>())),
    EvalBuilder<T,N>>;
static_assert(term_refine_returns_builder<double, 2>, "Term::refine_with → EvalBuilder (NDIM=2)");
static_assert(term_refine_returns_builder<double, 4>, "Term::refine_with → EvalBuilder (NDIM=4)");
static_assert(term_refine_returns_builder<double, 6>, "Term::refine_with → EvalBuilder (NDIM=6)");

// EvalBuilder must not compose further: no operator* with Term.
template <typename, typename, typename = void>
struct can_multiply_builder_by_term : std::false_type {};
template <typename B, typename Tm>
struct can_multiply_builder_by_term<B, Tm, std::void_t<
    decltype(std::declval<B>() * std::declval<Tm>())>>
    : std::true_type {};

static_assert(!can_multiply_builder_by_term<EvalBuilder<double,2>, Term<double,2>>::value,
    "EvalBuilder must be terminal (NDIM=2)");
static_assert(!can_multiply_builder_by_term<EvalBuilder<double,4>, Term<double,4>>::value,
    "EvalBuilder must be terminal (NDIM=4)");
static_assert(!can_multiply_builder_by_term<EvalBuilder<double,6>, Term<double,6>>::value,
    "EvalBuilder must be terminal (NDIM=6)");

}  // namespace

// ---------------------------------------------------------------------------
// One-NDIM test driver — combines the DSL operator-overload tests (A,B,C,D,
// S1,S2,O1,I1) and the builder API tests (R1..R8).
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

    test_output t1("composite_dsl: DSL + builder, NDIM=" + std::to_string(NDIM));
    t1.set_cout_to_terminal();
    std::cout << std::setprecision(8) << std::scientific;

    Function<T,LDIM> a  = FunctionFactory<T,LDIM>(world).f(gauss_a<LDIM>);
    Function<T,LDIM> b  = FunctionFactory<T,LDIM>(world).f(gauss_b<LDIM>);
    Function<T,LDIM> c  = FunctionFactory<T,LDIM>(world).f(gauss_c<LDIM>);
    Function<T,LDIM> d  = FunctionFactory<T,LDIM>(world).f(gauss_d<LDIM>);
    Function<T,LDIM> v  = FunctionFactory<T,LDIM>(world).f(gauss_v<LDIM>);

    Function<T,NDIM> u   = hartree_product(a, b);
    Function<T,NDIM> eri = FunctionFactory<T,NDIM>(world).is_on_demand().f(f12_lambda<NDIM>);

    const double tol_fn    = thresh * 10.0;
    const double tol_inner = thresh;

    // ============================================================
    // Section 1: operator-overload DSL
    // ============================================================

    // A: u * c(p1)
    {
        auto r = (u(p1,p2) * c(p1)).eval();
        t1.checkpoint(r.err(exact_uc<NDIM>), tol_fn, "A: u(p1,p2) * c(p1)");
    }
    // B: u * d(p2)
    {
        auto r = (u(p1,p2) * d(p2)).eval();
        t1.checkpoint(r.err(exact_ud<NDIM>), tol_fn, "B: u(p1,p2) * d(p2)");
    }
    // C: u * c(p1) * d(p2)
    {
        auto r = (u(p1,p2) * c(p1) * d(p2)).eval();
        t1.checkpoint(r.err(exact_ucd<NDIM>), tol_fn, "C: u(p1,p2) * c(p1) * d(p2)");
    }
    // D: u · (v(p1) + v(p2) + eri(p1,p2))
    {
        auto r = (u(p1,p2) * (v(p1) + v(p2) + eri(p1,p2))).eval();
        t1.checkpoint(r.err(exact_Vphi<NDIM>), tol_fn, "D: u · (v(p1)+v(p2)+eri)");
    }

    // S1: scalar prefactor
    {
        auto r = (2.0 * u(p1,p2) * c(p1)).eval();
        auto exact_2uc = [](const Vector<double,NDIM>& x){ return 2.0 * exact_uc<NDIM>(x); };
        t1.checkpoint(r.err(exact_2uc), tol_fn, "S1: 2.0 * u * c(p1)");
    }
    // S2: unary minus + signed sum
    {
        auto r = (-u(p1,p2)*c(p1) + 0.5*u(p1,p2)*d(p2)).eval();
        auto exact_s2 = [](const Vector<double,NDIM>& x){
            return -exact_uc<NDIM>(x) + 0.5*exact_ud<NDIM>(x);
        };
        t1.checkpoint(r.err(exact_s2), tol_fn, "S2: -u*c(p1) + 0.5*u*d(p2)");
    }

    // O1: op-apply barrier
    {
        auto bsh = BSHOperator<NDIM>(world, /*mu*/ 1.0, /*lo*/ 1e-6, /*eps*/ thresh*0.1);
        auto r_dsl = (apply(bsh, u(p1,p2)) * v(p1)).eval();
        Function<T,NDIM> bsh_u = apply(bsh, u);
        Function<T,NDIM> r_manual = multiply(bsh_u, v, /*particle*/ 1);
        const double diff = (r_dsl - r_manual).norm2();
        t1.checkpoint(diff, tol_fn, "O1: apply(BSH,u)*v(p1) consistency");
    }

    // I1: inner against on-demand bra (eri).
    // <eri | u·c(p1)·d(p2)>  =  <a·c | conv_β | b·d>  (β = 0.5)
    {
        const double beta = 0.5;
        auto conv = GaussOperator<LDIM>(world, beta);
        Function<T,LDIM> ac = a * c;
        Function<T,LDIM> bd = b * d;
        Function<T,LDIM> conv_bd = apply(conv, bd);
        const double ref = inner(ac, conv_bd);

        const double val = inner(eri(p1,p2), u(p1,p2)*c(p1)*d(p2));
        t1.checkpoint(std::abs(val - ref), tol_inner, "I1: <eri | u·c(p1)·d(p2)>");
    }

    // ============================================================
    // Section 2: builder API (LeafOpBase / EvalBuilder)
    // ============================================================

    using SC      = SeparatedConvolution<double, NDIM>;
    using SBox    = Specialbox_op<T, NDIM>;
    using LeafT   = Leaf_op<T, NDIM, SC, SBox>;
    using CuspBox = ElectronCuspyBox_op<T, NDIM>;
    using CuspLeaf= Leaf_op<T, NDIM, SC, CuspBox>;

    // R1: refine_with(default Leaf_op) matches plain .eval()
    {
        LeafT lop(u.get_impl().get());
        auto plain  = (u(p1,p2) * c(p1)).eval();
        auto built  = (u(p1,p2) * c(p1)).refine_with(lop).eval();
        const double diff = (plain - built).norm2();
        t1.checkpoint(diff, tol_fn, "R1: refine_with(default Leaf_op) == .eval()");
    }

    // R2: refine_with(ElectronCuspyBox_op)-driven Leaf_op on Vphi
    {
        CuspBox sbox;
        CuspLeaf cusp_leaf(u.get_impl().get(), sbox);
        auto built  = (u(p1,p2) * (v(p1) + v(p2) + eri(p1,p2)))
                          .refine_with(cusp_leaf).eval();
        t1.checkpoint(built.err(exact_Vphi<NDIM>), tol_fn,
                      "R2: ElectronCuspyBox_op refinement on Vphi");
    }

    // R3: .thresh() reflected on result
    {
        const double tight = thresh * 0.1;
        auto r = (u(p1,p2) * c(p1)).thresh(tight).eval();
        const double got = r.get_impl()->get_thresh();
        const double diff = std::abs(got - tight);
        t1.checkpoint(diff, 1e-15,
                      "R3: thresh() reflected on result.get_impl()->get_thresh()");
    }

    // R4: .truncate_mode() reflected on result
    {
        const int default_mode = FunctionDefaults<NDIM>::get_truncate_mode();
        const int picked = (default_mode == 0) ? 1 : 0;
        auto r = (u(p1,p2) * c(p1)).truncate_mode(picked).eval();
        const int got = r.get_impl()->get_truncate_mode();
        t1.checkpoint(got, picked, "R4: truncate_mode() reflected on result");
    }

    // R5: chained refine_with + thresh + truncate_mode
    {
        LeafT lop(u.get_impl().get());
        const double tight = thresh * 0.1;
        auto r = (u(p1,p2) * c(p1))
                     .refine_with(lop)
                     .thresh(tight)
                     .truncate_mode(1)
                     .eval();
        const bool thresh_ok = std::abs(r.get_impl()->get_thresh() - tight) < 1e-15;
        const bool mode_ok   = (r.get_impl()->get_truncate_mode() == 1);
        const double err = r.err(exact_uc<NDIM>);
        const double composite_error = err + (thresh_ok ? 0.0 : 1.0) + (mode_ok ? 0.0 : 1.0);
        t1.checkpoint(composite_error, tol_fn,
                      "R5: refine_with + thresh + truncate_mode chained");
    }

    // R6: chaining order is irrelevant
    {
        LeafT lop(u.get_impl().get());
        const double tight = thresh * 0.1;
        auto a_order = (u(p1,p2) * c(p1)).refine_with(lop).thresh(tight).truncate_mode(1).eval();
        auto b_order = (u(p1,p2) * c(p1)).truncate_mode(1).thresh(tight).refine_with(lop).eval();
        const double diff = (a_order - b_order).norm2();
        t1.checkpoint(diff, tol_fn, "R6: builder order is irrelevant");
    }

    // R7: inner() with refine_with on RHS — value matches the un-refined inner.
    {
        LeafT lop(u.get_impl().get());
        const double plain = inner(eri(p1,p2), u(p1,p2)*c(p1)*d(p2));
        const double built = inner(eri(p1,p2), (u(p1,p2)*c(p1)*d(p2)).refine_with(lop));
        t1.checkpoint(std::abs(built - plain), tol_inner,
                      "R7: inner with refine_with on RHS matches plain inner");
    }

    // R8: composite_product(LeafOpBase) overload matches default
    {
        LeafT lop(u.get_impl().get());
        std::vector<composite_term<T,NDIM>> terms(1);
        terms[0].ket = u;
        terms[0].factor1 = c;
        auto plain  = composite_product<T,NDIM>(world, terms);

        std::vector<composite_term<T,NDIM>> terms2(1);
        terms2[0].ket = u;
        terms2[0].factor1 = c;
        auto polym  = composite_product<T,NDIM>(world, terms2, lop);
        const double diff = (plain - polym).norm2();
        t1.checkpoint(diff, tol_fn,
                      "R8: composite_product(LeafOpBase) overload matches default");
    }

    return t1.end();
}


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world, argc, argv);

    const long   k      = 6;
    const double thresh = 3.e-4;

    int rc = 0;
    rc |= run_tests<double, 2>(world, k, thresh);
    rc |= run_tests<double, 4>(world, k, thresh);
    rc |= run_tests<double, 6>(world, k, thresh);

    finalize();
    return rc;
}
