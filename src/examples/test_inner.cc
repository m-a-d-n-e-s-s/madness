

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/mra/mra.h>

using namespace madness;

typedef std::shared_ptr< FunctionFunctorInterface<double,3> > functorT;

static const double R = 1.4;    // bond length
static const double L = 5.0*R; // box size
static const long k = 8;        // wavelet order
static const double thresh = 1e-6; // precision

static double alpha_func(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return (sin(x*x + y*y + z*z));
};

static double beta_func(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return (exp(- x*x - y*y - z*z));
};

class alpha_functor : public FunctionFunctorInterface<double,3> {
private:
    double coeff;
public:
    alpha_functor(double coeff=1.0) : coeff(coeff) {}

    virtual double operator()(const coord_3d& r) const {
        return (coeff * sin(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]));
    }
};

class beta_functor : public FunctionFunctorInterface<double,3> {
private:
    double coeff;
public:
    beta_functor(double coeff=1.0) : coeff(coeff) {}

    virtual double operator()(const coord_3d& r) const {
        return (coeff * exp(- r[0]*r[0] - r[1]*r[1] - r[2]*r[2]));
    }
};

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);

    startup(world,argc,argv);
    std::cout.precision(6);

    FunctionDefaults<3>::set_defaults(world);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_initial_level(5);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_cubic_cell(-L/2, L/2);

    real_function_3d alpha = real_factory_3d(world).f(alpha_func);
    real_function_3d beta = real_factory_3d(world).f(beta_func);
    real_functor_3d alpha_ffi = real_functor_3d(new alpha_functor());
    real_functor_3d beta_ffi = real_functor_3d(new beta_functor());

    if (world.rank() == 0) print("Testing compression and reconstruction.");
    double norm = beta.norm2();
    beta.compress();
    double normc = beta.norm2();
    beta.reconstruct();
    double normr = beta.norm2();

    if (world.rank() == 0) print("Calculating inner product of two Functions.");
    double aa = alpha.inner(alpha);
    double ab = alpha.inner(beta);
    double bb = beta.inner(beta);

    if (world.rank() == 0) print("Calculating inner product of Function with external function.");
    double aa_f = alpha.inner_ext(alpha_func);
    double ab_f = alpha.inner_ext(beta_func);
    double ba_f = beta.inner_ext(alpha_func);
    double bb_f = beta.inner_ext(beta_func);

    if (world.rank() == 0) print("Calculating inner product of Function with FunctionFunctorInterface.");
    double aa_ffi = alpha.inner_ext(alpha_ffi);
    double ab_ffi = alpha.inner_ext(beta_ffi);
    double ba_ffi = beta.inner_ext(alpha_ffi);
    double bb_ffi = beta.inner_ext(beta_ffi);

    if (world.rank() == 0) {
        print("Norms are ", norm, normc, normr);
        print("**********************************************************");
        printf("<a|a> (using inner() with Function) =                      %7.10f\n", aa);
        printf("<a|a> (using inner_ext() with external function) =         %7.10f\n", aa_f);
        printf("<a|a> (using inner_ext() with FunctionFunctor Interface) = %7.10f\n", aa_ffi);
        print("**********************************************************");
        printf("<b|b> (using inner() with Function) =                      %7.10f\n", bb);
        printf("<b|b> (using inner_ext() with external function) =         %7.10f\n", bb_f);
        printf("<b|b> (using inner_ext() with FunctionFunctor Interface) = %7.10f\n", bb_ffi);
        print("**********************************************************");
        printf("<a|b> (using inner() with Function) =                      %7.10f\n", ab);
        printf("<a|b> (using inner_ext() with external function) =         %7.10f\n", ab_f);
        printf("<a|b> (using inner_ext() with FunctionFunctor Interface) = %7.10f\n", ab_ffi);
        printf("<b|a> (using inner_ext() with external function) =         %7.10f\n", ba_f);
        printf("<b|a> (using inner_ext() with FunctionFunctor Interface) = %7.10f\n", ba_ffi);
        print("**********************************************************");
    }

    world.gop.fence();

    finalize();
    return 0;
}
