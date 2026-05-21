/*
  Tests for the composite_dsl DSL: operator-overload syntax over Function expressions.

  Coverage at NDIM=6 (LDIM=3):

    function mode:
      A) u(p1,p2) * c(p1)                              → composite_product
      B) u(p1,p2) * d(p2)
      C) u(p1,p2) * c(p1) * d(p2)
      D) u(p1,p2) * (v(p1) + v(p2) + eri(p1,p2))       → make_Vphi-style sum

    scalar / unary minus:
      S1) 2.0 * u(p1,p2) * c(p1)
      S2) -u(p1,p2)*c(p1) + 0.5*u(p1,p2)*d(p2)

    op-apply barrier:
      O1) apply(BSH, u(p1,p2)) * v(p1)                  → eager apply, then DSL

    inner mode (on-demand bra):
      I1) inner(eri(p1,p2), u(p1,p2)*c(p1)*d(p2))

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
#include <madness/mra/operator.h>
#include <madness/world/test_utilities.h>
#include <iomanip>

using namespace madness;

namespace {

constexpr std::size_t LDIM = 3;
constexpr std::size_t NDIM = 6;

double gauss_a (const Vector<double,LDIM>& r) { return std::exp(-1.0 * inner(r,r)); }
double gauss_b (const Vector<double,LDIM>& r) { return std::exp(-1.3 * inner(r,r)); }
double gauss_c (const Vector<double,LDIM>& r) { return std::exp(-0.7 * inner(r,r)); }
double gauss_d (const Vector<double,LDIM>& r) { return std::exp(-0.5 * inner(r,r)); }
double gauss_v (const Vector<double,LDIM>& r) { return std::exp(-0.3 * inner(r,r)); }

double exact_uc (const Vector<double,NDIM>& r) {
    Vector<double,LDIM> r1, r2;
    for (std::size_t i=0; i<LDIM; ++i) { r1[i]=r[i]; r2[i]=r[i+LDIM]; }
    return std::exp(-1.7*inner(r1,r1) - 1.3*inner(r2,r2));
}
double exact_ud (const Vector<double,NDIM>& r) {
    Vector<double,LDIM> r1, r2;
    for (std::size_t i=0; i<LDIM; ++i) { r1[i]=r[i]; r2[i]=r[i+LDIM]; }
    return std::exp(-1.0*inner(r1,r1) - 1.8*inner(r2,r2));
}
double exact_ucd (const Vector<double,NDIM>& r) {
    Vector<double,LDIM> r1, r2;
    for (std::size_t i=0; i<LDIM; ++i) { r1[i]=r[i]; r2[i]=r[i+LDIM]; }
    return std::exp(-1.7*inner(r1,r1) - 1.8*inner(r2,r2));
}
double exact_Vphi (const Vector<double,NDIM>& r) {
    Vector<double,LDIM> r1, r2;
    for (std::size_t i=0; i<LDIM; ++i) { r1[i]=r[i]; r2[i]=r[i+LDIM]; }
    const double ur  = std::exp(-1.0*inner(r1,r1) - 1.3*inner(r2,r2));
    const double v1  = std::exp(-0.3*inner(r1,r1));
    const double v2_ = std::exp(-0.3*inner(r2,r2));
    const double eri = std::exp(-0.5 * inner(r1-r2, r1-r2));
    return ur * (v1 + v2_ + eri);
}
double f12_lambda(const Vector<double,NDIM>& r) {
    Vector<double,LDIM> r1, r2;
    for (std::size_t i=0; i<LDIM; ++i) { r1[i]=r[i]; r2[i]=r[i+LDIM]; }
    return std::exp(-0.5 * inner(r1-r2, r1-r2));
}

}  // namespace

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world, argc, argv);

    const long k = 6;
    const double thresh = 3.e-4;

    FunctionDefaults<LDIM>::set_k(k);
    FunctionDefaults<NDIM>::set_k(k);
    FunctionDefaults<LDIM>::set_thresh(thresh*0.01);
    FunctionDefaults<NDIM>::set_thresh(thresh);
    FunctionDefaults<LDIM>::set_tensor_type(TT_FULL);
    FunctionDefaults<NDIM>::set_tensor_type(TT_2D);
    FunctionDefaults<LDIM>::set_cubic_cell(-8.0, 8.0);
    FunctionDefaults<NDIM>::set_cubic_cell(-8.0, 8.0);

    test_output t1("composite_dsl: operator-overload DSL");
    t1.set_cout_to_terminal();
    std::cout << std::setprecision(8) << std::scientific;

    using T = double;
    Function<T,LDIM> a  = FunctionFactory<T,LDIM>(world).f(gauss_a);
    Function<T,LDIM> b  = FunctionFactory<T,LDIM>(world).f(gauss_b);
    Function<T,LDIM> c  = FunctionFactory<T,LDIM>(world).f(gauss_c);
    Function<T,LDIM> d  = FunctionFactory<T,LDIM>(world).f(gauss_d);
    Function<T,LDIM> v  = FunctionFactory<T,LDIM>(world).f(gauss_v);

    Function<T,NDIM> u   = hartree_product(a, b);
    Function<T,NDIM> eri = FunctionFactory<T,NDIM>(world).is_on_demand().f(f12_lambda);

    const double tol_fn    = thresh * 10.0;
    const double tol_inner = thresh;

    // A: u * c(p1) ------------------------------------------------------
    {
        auto r = (u(p1,p2) * c(p1)).eval();
        t1.checkpoint(r.err(exact_uc), tol_fn, "A: u(p1,p2) * c(p1)");
    }
    // B: u * d(p2) ------------------------------------------------------
    {
        auto r = (u(p1,p2) * d(p2)).eval();
        t1.checkpoint(r.err(exact_ud), tol_fn, "B: u(p1,p2) * d(p2)");
    }
    // C: u * c(p1) * d(p2) ---------------------------------------------
    {
        auto r = (u(p1,p2) * c(p1) * d(p2)).eval();
        t1.checkpoint(r.err(exact_ucd), tol_fn, "C: u(p1,p2) * c(p1) * d(p2)");
    }
    // D: u · (v(p1) + v(p2) + eri(p1,p2)) ------------------------------
    {
        auto r = (u(p1,p2) * (v(p1) + v(p2) + eri(p1,p2))).eval();
        t1.checkpoint(r.err(exact_Vphi), tol_fn, "D: u · (v(p1)+v(p2)+eri)");
    }

    // S1: scalar prefactor ---------------------------------------------
    {
        auto r = (2.0 * u(p1,p2) * c(p1)).eval();
        auto exact_2uc = [](const Vector<double,NDIM>& x){ return 2.0 * exact_uc(x); };
        t1.checkpoint(r.err(exact_2uc), tol_fn, "S1: 2.0 * u * c(p1)");
    }
    // S2: unary minus + signed sum -------------------------------------
    {
        auto r = (-u(p1,p2)*c(p1) + 0.5*u(p1,p2)*d(p2)).eval();
        auto exact_s2 = [](const Vector<double,NDIM>& x){
            return -exact_uc(x) + 0.5*exact_ud(x);
        };
        t1.checkpoint(r.err(exact_s2), tol_fn, "S2: -u*c(p1) + 0.5*u*d(p2)");
    }

    // O1: op-apply barrier ---------------------------------------------
    {
        auto bsh = BSHOperator<NDIM>(world, /*mu*/ 1.0, /*lo*/ 1e-6, /*eps*/ thresh*0.1);
        auto r_dsl = (apply(bsh, u(p1,p2)) * v(p1)).eval();
        Function<T,NDIM> bsh_u = apply(bsh, u);
        Function<T,NDIM> r_manual = multiply(bsh_u, v, /*particle*/ 1);
        const double diff = (r_dsl - r_manual).norm2();
        t1.checkpoint(diff, tol_fn, "O1: apply(BSH,u)*v(p1) consistency");
    }

    // I1: inner against on-demand bra ----------------------------------
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

    return t1.end();
}
