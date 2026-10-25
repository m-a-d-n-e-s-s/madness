#include <madness/mra/mra.h>
#include <madness/tensor/tensor.h>
#include <madness/world/print.h>
#include <cmath>
#include <iostream>
#include <functional>
#include <madness/mra/convolution1d.h>

#include <iomanip>
#include <string>
#include <vector>
#include <cassert>
#include <algorithm>

using namespace madness;

// Computes the k x k extension operators P^1 and P^{-1} that maintain the quadratic Taylor expansion of a function at the right and left boundaries of the box [0, 1].
// P^{+1} = P_right maps box [0, 1] coeffs -> box [1, 2] coeffs
// P^{-1} = P_left  maps box [0, 1] coeffs -> box [-1, 0] coeffs
void make_taylor_extension_operators(long k,
                                     Tensor<double>& P_right,
                                     Tensor<double>& P_left)
{
    // Allocate k x k MADNESS tensors
    P_right = Tensor<double>(k, k);
    P_left  = Tensor<double>(k, k);

    // Initialize all entries to zero (since rows i >= 3 must be strictly zero)
    P_right = 0.0;
    P_left  = 0.0;

    const double inv_sqrt3 = 1.0 / std::sqrt(3.0);
    const double inv_sqrt5 = 1.0 / std::sqrt(5.0);

    for (long j = 0; j < k; ++j) {
        double sj = std::sqrt(2.0 * j + 1.0);
        double j_val = static_cast<double>(j);
        double j_p1 = j_val + 1.0;
        double jj1 = j_val * j_p1;

        // Boundary value and derivative components at x = 1
        double d0 = sj;
        double d1 = jj1 * sj;
        double d2 = (j >= 2) ? 0.5 * (j_val - 1.0) * jj1 * (j_val + 2.0) * sj : 0.0;

        // Row 0: c^1_0 = v0 + (1/2)*v1 + (1/6)*v2
        double p0 = d0 + 0.5 * d1 + (1.0 / 6.0) * d2;

        // Row 1: c^1_1 = (1/(2*sqrt(3)))*v1 + (1/(4*sqrt(3)))*v2
        double p1 = (0.5 * inv_sqrt3) * d1 + (0.25 * inv_sqrt3) * d2;

        // Row 2: c^1_2 = (1/(12*sqrt(5)))*v2
        double p2 = (inv_sqrt5 / 12.0) * d2;

        // Store into P_right
        P_right(0, j) = p0;
        if (k > 1) P_right(1, j) = p1;
        if (k > 2) P_right(2, j) = p2;

        // By parity reflection: P_left(i, j) = (-1)^(i + j) * P_right(i, j)
        double parity_j = (j % 2 == 0) ? 1.0 : -1.0;

        P_left(0, j) =  parity_j * p0;
        if (k > 1) P_left(1, j) = -parity_j * p1;
        if (k > 2) P_left(2, j) =  parity_j * p2;
    }
}

#include <cmath>
#include <algorithm>
#include <vector>
#include <madness/tensor/tensor.h>

using madness::Tensor;

// Computes the k x k extension operators P^{+1} (P_right) and P^{-1} (P_left)
// that maintain the Taylor expansion of degree `max_degree` (order = max_degree + 1)
// at the right (x = 1) and left (x = 0) boundaries of the box [0, 1].
//
// Parameters:
//   k          - Number of Legendre scaling functions in the basis
//   max_degree - Polynomial degree of the Taylor approximation:
//                0 = constant, 1 = linear, 2 = quadratic, 3 = cubic, etc.
//                (Clamped internally to [0, k - 1])
//   P_right    - Output k x k tensor mapping box [0, 1] coeffs -> box [1, 2] coeffs
//   P_left     - Output k x k tensor mapping box [0, 1] coeffs -> box [-1, 0] coeffs
void make_taylor_extension_operatorsX(long k,
                                      long max_degree,
                                      Tensor<double>& P_right,
                                      Tensor<double>& P_left)
{
    // Allocate and zero-initialize k x k MADNESS tensors
    P_right = Tensor<double>(k, k);
    P_left  = Tensor<double>(k, k);
    P_right = 0.0;
    P_left  = 0.0;

    if (k <= 0 || max_degree < 0) return;

    // The active Taylor degree cannot exceed the maximum polynomial degree (k - 1)
    const long M = std::min(max_degree, k - 1);

    // Precompute the Taylor-to-Legendre projection coefficients C(i, m):
    // C(i, m) = int_0^1 (y^m / m!) * phi_i(y) dy
    //         = sqrt(2i + 1) * m! / ((m - i)! * (m + i + 1)!)  for 0 <= i <= m <= M
    std::vector<double> C((M + 1) * (M + 1), 0.0);
    for (long i = 0; i <= M; ++i) {
        const double si = std::sqrt(2.0 * i + 1.0);

        // Base case m = i: C(i, i) = sqrt(2i + 1) * i! / (2i + 1)!
        double c_im = si;
        for (long n = i + 1; n <= 2 * i + 1; ++n) {
            c_im /= static_cast<double>(n);
        }
        C[i * (M + 1) + i] = c_im;

        // Recurrence for m > i:
        // C(i, m) = C(i, m - 1) * m / ((m - i) * (m + i + 1))
        for (long m = i + 1; m <= M; ++m) {
            c_im *= static_cast<double>(m) /
                    (static_cast<double>(m - i) * static_cast<double>(m + i + 1));
            C[i * (M + 1) + m] = c_im;
        }
    }

    // Reusable buffer for boundary derivatives d_m = phi_j^{(m)}(1)
    std::vector<double> d(M + 1, 0.0);

    for (long j = 0; j < k; ++j) {
        const long max_m = std::min(M, j);

        // Compute derivatives d[m] = phi_j^{(m)}(1) for 0 <= m <= min(M, j)
        // d[0] = sqrt(2j + 1)
        // d[m] = d[m - 1] * (j + m) * (j - m + 1) / m
        d[0] = std::sqrt(2.0 * j + 1.0);
        for (long m = 1; m <= max_m; ++m) {
            double num = static_cast<double>(j + m) * static_cast<double>(j - m + 1);
            d[m] = d[m - 1] * (num / static_cast<double>(m));
        }

        const double parity_j = (j % 2 == 0) ? 1.0 : -1.0;

        // Assemble rows 0 <= i <= min(M, j)
        // (Note: for i > j, P_right(i, j) is strictly 0 because phi_j has degree j < i)
        for (long i = 0; i <= max_m; ++i) {
            double p_ij = 0.0;
            for (long m = i; m <= max_m; ++m) {
                p_ij += C[i * (M + 1) + m] * d[m];
            }

            P_right(i, j) = p_ij;

            // Parity reflection: P_left(i, j) = (-1)^(i + j) * P_right(i, j)
            const double parity_ij = (i % 2 == 0) ? parity_j : -parity_j;
            P_left(i, j) = parity_ij * p_ij;
        }
    }
}


// projects a function f(x) defined on [0, 1] onto the first k Legendre basis functions
Tensor<double> project(std::function<double(double)> f, long k) {
    Tensor<double> coeffs(k), phi(k), x(k), w(k);
    gauss_legendre(k, 0.0, 1.0, x.ptr(), w.ptr());
        
    for (long j = 0; j < k; ++j) {
        double wfj = w(j) * f(x(j));
        legendre_scaling_functions(x(j), k, phi.ptr());
        for (long i = 0; i < k; ++i) {
            coeffs(i) += wfj * phi(i);
        }
    }
    return coeffs;
}

// evaluates a function represented by Legendre coefficients at a point x in [0, 1]
double evaluate(const Tensor<double>& coeffs, double x) {
    long k = coeffs.size();
    Tensor<double> phi(k);
    legendre_scaling_functions(x, k, phi.ptr());
    return coeffs.trace(phi);
}

double F(double x) {
    return std::sin(x);
}

void test0() {
    print("\nTesting Taylor extension operators for sin(x) on [0,1]\n");
    long k = 10; // Polynomial order    

    Tensor<double> Pp, Pm;
    make_taylor_extension_operators(k, Pp, Pm);

    print("P+1:", Pp);
    print("P-1:", Pm);

    // Test with sin(x) on [0,1]
    auto c0 = project(F, k);

    // Test the projection
    double xtest = 0.5;
    double f_exact = F(xtest);
    double f_approx = evaluate(c0, xtest);
    print("f_exact at x =", xtest, ":", f_exact);
    print("f_approx at x =", xtest, ":", f_approx);

    // The expected exact coeffs in box [1, 2] can be computed by projecting the taylor series expansion of sin(x+1) on the right
    auto G = [](double x) { return std::sin(1.0) + std::cos(1.0)*x - 0.5*std::sin(1.0)*x*x; };
    auto c_1_exact = project(G, k);
    auto H = [](double x) { return x-1.0; }; // sin(x) is linear at the origin
    auto c_minus1_exact = project(H, k);

    // Compute coefficients in box [1, 2] via matrix-vector multiplication: c1 = Pp * c0 
    Tensor<double> c1 = inner(Pp, c0);
    Tensor<double> c_minus1 = inner(Pm, c0);

    print("                       c0 (coeffs in [0,1]):", c0);
    print("                c1 (coeffs in [1,2] via Pp):", c1);
    print("  c1_exact (coeffs in [1,2] via projection):", c_1_exact);
    print("             c-1 (coeffs in [-1,0] via P-1):", c_minus1);
    print("c-1_exact (coeffs in [-1,0] via projection):", c_minus1_exact);
}

/*
  Second derivative matrix --- no BC as if polyn extended across all space
    (\bm{D}_2)_{i,q} &= \begin{cases}
    4 \sqrt{2i+1}\sqrt{2q+1} \, \frac{(q-i)(i+q+1)}{2}, & q \ge i+2, \, q \equiv i \pmod 2,\\
    0, & \text{otherwise}.

 */
Tensor<double> make_D2(int k) {
    Tensor<double>D2(k,k);
    for (int i=0; i<k; i++) {
        for (int q=i+2; q<k; q+=2) {
            D2(i,q) = 4.0*std::sqrt(2*i+1)*std::sqrt(2*q+1)*(q-i)*(i+q+1)/2.0;
        }
    }
    return D2;
}

// Matrix representation of diagonal block of Gaussian convolution with no BC as if polyn extended across all space
/*

 */
Tensor<double> make_T00(int k, int n, double expnt) {
    Tensor<double> M(k,k);

    Tensor<double> term(k,k);
    for (int i=0; i<k; i++) term(i,i) = 1.0;

    double hn = std::pow(2.0, -n);

    Tensor<double> D2 = make_D2(k).scale(1.0 / (4.0*hn*hn*expnt));

    M += term;
    
    for (int r=1; r<=(k-1)/2; r++) {
        //for (int r=1; r<=1; r++) {
        term = inner(D2, term).scale(1.0 / r);
        M += term;
    }
    
    return M;
}

void testD2() {
    print("\nTesting second derivative matrix D2\n");

    // Second derivative of sin(x) is -sin(x)
    long k = 10; // Polynomial order    
    auto c0 = project(F, k);
    auto D2c0 = -inner(make_D2(k), c0);
    auto D4c0 = -inner(make_D2(k), D2c0);
    auto D6c0 = -inner(make_D2(k), D4c0);
    print(c0);
    print(D2c0);
    print(D4c0);
    print(D6c0);
    print("error |D2 c0 - c0|:", (D2c0 - c0).normf());
    
}

void test1() {
    print("\nTesting Taylor extension operators for sin(x) on [0,1] with convolution with a Gaussian kernel\n");
    // The extension seems to be working, now let's try the Gaussian convolution operator manually
    // int(g(x-y,a)*sin(y), y=-inf..inf) = sin(x) * exp(-1/4*a), where g is a normalized gaussian with exponent a
    // that correponds to width sigma = 1/sqrt(2*a) or a = 1/(2*sigma^2)
    // if sigma << 1 (box width), we test convolving with the extension rather than the original function

    //double expnt = 1.0 / (2 * 0.16667);

    long k = 8;

    print("D2"); print(make_D2(k));
    
    Tensor<double> Pp, Pm;
    make_taylor_extension_operators(k, Pp, Pm);
    auto c0 = project(F, k);
    
    print("T0 is the 00 block of the Gaussian convolution");
    print("R0 is the 00 block of the Gaussian convolution with quadratic extension to nearest neighbors");
    print("M0 is the 00 block of the Gaussian convolution with full extension over all space (analytic form from D2 expansion)");
    print("rate of convergence = log(err/errprev)/log(sigma/sigmaprev)\n");
    printf("scale  sigma        expnt         err     err/errp    rate      errM   errM/errMp  rateM     Tr(T0)/k  Tr(R0)/k   Tr(M0)/k   |R0-M0|\n");
    printf("----- --------   ------------   --------  --------  --------  --------  --------  --------   --------  --------   --------   -------\n");
    double sigprev = 2.0;
    double errprev = 1.0;
    double errMprev= 1.0;
    for (int iscale=0; iscale<=16; iscale++) {
        double sigma = std::pow(2.0, -0.5*iscale);
        double expnt = 1.0 / (2 * sigma *sigma);
        double coeff = std::sqrt(expnt / M_PI);
        GaussianConvolution1D<double> gc(k, coeff, expnt, 0, LatticeRange(0));
        Tensor<double> r0 = gc.rnlij(0, 0);
        Tensor<double> rm = gc.rnlij(0, -1);
        Tensor<double> rp = gc.rnlij(0, 1);
        // print("[ 0,1] <--  [1,2]");
        // print(rm);
        // print("[ 0,1] <--> [0,1]");
        // print(r0);
        // print("[-1,0]  --> [0,1]");
        // print(rp);

        Tensor<double> R0 = inner(rp, Pm) + r0 + inner(rm, Pp);
        Tensor<double> M0 = make_T00(k, 0, expnt);
        double errM0 = (R0 - M0).normf();
        // print("M0"); print(M0);
        // print("R0"); print(R0);

        double alpha = 0.0, Alpha = 0.0, alphaM = 0.0;
        for (long i = 0; i < k; ++i) {
            alpha += r0(i, i);
            Alpha += R0(i, i);
            alphaM += M0(i, i);
        }
        alpha /= k;
        Alpha /= k;
        alphaM /= k;

        // print("R0 = rp*P-1 + r0 + rm*Pp");
        // print(R0);

        Tensor<double> c0new = inner(R0, c0);
        Tensor<double> c0newM = inner(M0, c0);
        Tensor<double> c0new_exact = c0 * std::exp(-1.0 / (4.0 * expnt));
        double err = (c0new - c0new_exact).normf() / c0new_exact.normf();
        double errM= (c0newM- c0new_exact).normf() / c0new_exact.normf();
        double rate = std::log(err/errprev)/log(sigma/sigprev);
        double rateM = std::log(errM/errMprev)/log(sigma/sigprev);
        printf(" %2d   %.6f   %.6e   %.2e  %.6f  %.6f  %.2e  %.6f  %.6f   %.6f  %.6f   %.6f   %.2e\n", iscale, sigma, expnt, err, err/errprev, rate, errM, errM/errMprev, rateM, alpha, Alpha, alphaM, errM0);
        errprev = err;
        errMprev= errM;
        sigprev = sigma;
    }
}

void test2() {
    print("\nPlot the error to examine behavior at edges of the box");
    long k = 8;
    Tensor<double> Pp, Pm;
    make_taylor_extension_operators(k, Pp, Pm);
    double sigma = 0.0625;
    double expnt = 1.0 / (2 * sigma *sigma);
    double coeff = std::sqrt(expnt / M_PI);
    GaussianConvolution1D<double> gc(k, coeff, expnt, 0, LatticeRange(0));
    Tensor<double> r0 = gc.rnlij(0, 0);
    Tensor<double> rm = gc.rnlij(0, -1);
    Tensor<double> rp = gc.rnlij(0, 1);
    Tensor<double> R0 = inner(rp, Pm) + r0 + inner(rm, Pp);
    auto c0 = project(F, k);
    Tensor<double> c0new = inner(R0, c0);
    long npt = 101;
    double h = 1.0/(npt-1);
    Tensor<double> x(npt), f_exact(npt), f_approx(npt), err(npt);
    // open plot data file
    auto file = fopen("gaussian_convolution_error.dat", "w");
    fprintf(file, "# x    f_exact f_approx err\n");
    for (long i=0; i<npt; i++) {
        double x = i*h;
        f_exact(i) = F(x) * std::exp(-1.0 / (4.0 * expnt));
        f_approx(i) = evaluate(c0new, x);
        err(i) = std::abs(f_exact(i)-f_approx(i));
        fprintf(file, "%.6f %.6e %.6e %.6e\n", x, f_exact(i), f_approx(i), err(i));
    }
    fclose(file);
}

void test3() {
    print("\nExamine the norm of the NS operator at different scales\n");
    int k = 8;
    double sigma = 0.0625;
    double expnt = 1.0 / (2 * sigma *sigma);
    double coeff = std::sqrt(expnt / M_PI);
    GaussianConvolution1D<double> gc(k, coeff, expnt, 0, LatticeRange(0));
    Translation lx = 0;
    print("n l   sigma/h      Rnorm         Tnorm       Rnormf      Tnormf        NSnormf");
    for (Level n = 0; n<=8; n++) {
        const ConvolutionData1D<double>* p = gc.nonstandard(n, lx);
        double h = 1.0 / (1<<n);
        print(n, lx, sigma/h, p->Rnorm, p->Tnorm, p->Rnormf, p->Tnormf, p->NSnormf);
    }
}

void test4(World& world) {
    print("\nTesting 3D norm estimations for blocks of Gaussian convolution operator\n");
    int k = 8;
    double thresh = 1e-6;
    double L = 10.0;
    double lo = 1e-6;
    
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_autorefine(false);
    FunctionDefaults<3>::set_cubic_cell(-L, L);

    double sigma = 0.0625;
    Tensor<double> expnt(1), coeff(1);
    expnt[0] = 1.0 / (2 * sigma *sigma);
    coeff[0] = std::pow(expnt[0] / M_PI, 1.5);
    printf("sigma=%.6f, expnt=%.6f, coeff=%.6f\n\n", sigma, expnt[0], coeff[0]);

    SeparatedConvolution<double,3> G(world, coeff, expnt, lo, thresh);

    auto power_iter = [=](const SeparatedConvolution<double,3>& G, Level n, const Key<3>& d) {
        const int max_iter = 100;
        Key<3> source(n), dm(n,-1*d.translation());
        
        Tensor<double> c(2*k,2*k,2*k);
        c.fillrandom();
        c -= 0.5;
        double cnormprev = 1e300;
        double cnorm = c.normf();
        for (int iter=0; iter<max_iter; iter++) {
            c *= 1.0 / cnorm;
            c = G.apply(source, d,  c, 0.0);
            c = G.apply(source, dm, c, 0.0); // effectively the adjoint operator
            cnorm = c.normf();
            //print(iter,cnorm);
            if (iter >= 2 && std::abs(cnorm - cnormprev) <= 1e-2 * cnorm) {
                return std::sqrt(cnorm);
            }
            cnormprev = cnorm;
        }
        throw std::runtime_error("Power iteration did not converge");
    };

    Vector<Translation,3> t000{0,0,0}, t001{0,0,1}, t011{0,1,1}, t111{1,1,1}, t002{0,0,2};
    
    printf("                     d000                  d001                  d011                  d111                  d002\n");
    printf("              ------------------    ------------------    ------------------    ------------------    ------------------\n");
    printf("    Level      normf     opnorm      normf     opnorm      normf     opnorm      normf     opnorm      normf     opnorm\n");
    printf("    -----     ------------------    ------------------    ------------------    ------------------    ------------------\n");
    for (Level n=0; n<=10; n++) {
        Key<3> source(n);
        std::array<Key<3>,5> ds{Key<3>(n,t000), Key<3>(n,t001), Key<3>(n,t011), Key<3>(n,t111), Key<3>(n,t002)};
        printf("      %2d: ", n);
        for (const auto& d : ds) {
            printf("    %.2e  %.2e", G.norm(n,d,source), power_iter(G,n, d));
        }
        printf("\n");
    }

}

void test5(World& world) {
    print("\n Testing scale invariance of Gaussian convolution matrices --- r[n,l,t] = r[n+1,l,4t] \n");
    int k = 8;
    double thresh = 1e-6;
    double L = 10.0;
    double lo = 1e-6;
    
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_autorefine(false);
    FunctionDefaults<3>::set_cubic_cell(-L, L);

    double expnt = 16.0;
    double coeff = std::pow(expnt / M_PI, 1.5);

    GaussianConvolution1D<double> gc0(k, coeff, expnt, 0, LatticeRange(0));
    Tensor<double> r0 = gc0.rnlij(0, 0);
    Tensor<double> rm = gc0.rnlij(0, -1);
    Tensor<double> rp = gc0.rnlij(0, 1);
        
    // Compute 4x exponent 1 level down
    GaussianConvolution1D<double> gc1(k, 2.0*coeff, 4.0*expnt, 0, LatticeRange(0));
    Tensor<double> r0x = gc1.rnlij(1, 0);
    Tensor<double> rmx = gc1.rnlij(1, -1);
    Tensor<double> rpx = gc1.rnlij(1, 1);
    
    // Check scale invariance
    double rel_err0 = (r0x - r0).normf() / r0.normf();
    double rel_errm = (rmx - rm).normf() / rm.normf();
    double rel_errp = (rpx - rp).normf() / rp.normf();
    
    printf("rel_err0 = %.6e, rel_errm = %.6e, rel_errp = %.6e\n", rel_err0, rel_errm, rel_errp);
}

class NormalizedGaussian : public FunctionFunctorInterface<double,3> {
private:
    double expnt;
    double coeff;
    
public:
    NormalizedGaussian(double expnt)
        : expnt(expnt)
        , coeff(std::pow(expnt / M_PI, 1.5))
    {}
        
    double operator()(const coord_3d& r) const {
        double x = r[0]-0.1; // displace so that refine deeply near high-level dyadic divide
        return coeff * std::exp(-expnt * (x*x + r[1]*r[1] + r[2]*r[2]));
    }
};

void test6(World& world) {
    print("\n Testing convolution of Gaussians to expose noise from Gibbs oscillations\n");
    int k = 10;
    double thresh = 1e-8;
    double L = 10.0;
    double lo = 1e-6;
    
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_truncate_mode(0);
    FunctionDefaults<3>::set_autorefine(false);
    FunctionDefaults<3>::set_cubic_cell(-L, L);

    // The reference Gaussian --- want it fairly spiky so have deep adaptive refinement
    const double a = 16.0;
    const double sigma_a = 1.0 / std::sqrt(2.0*a);
    real_function_3d f = real_factory_3d(world).functor(NormalizedGaussian(a));
    f.truncate();
    f.reconstruct();
    Level depth = f.max_depth();
    double userh = (2*L) * std::pow(2.0, -depth);
    printf("a=%.6f sigma_a=%.6f norm=%.6f depth=%d userh=%.6f  input_err=%.2e\n\n", a, sigma_a, f.norm2(), depth, userh, f.err(NormalizedGaussian(a)));

    printf(" userh/sigma     expnt      |G*f-exact|   |newG*f-exact|   |newG2*f-exact|\n");
    printf("                                                    order:");
    for (int i=2; i<k; i+=2) {
        printf("    %2d         ", i);
    }
    printf("\n");
    printf(" -----------   ----------   -----------   --------------   ---------------\n");
    for (int iscale=0; iscale<=40; iscale++) {
        const double sigma = std::pow(2.0, -0.5*iscale);
        const double b = 1.0 / (2 * sigma *sigma);
        
        Tensor<double> expnt(1), coeff(1);
        expnt[0] = b;
        coeff[0] = std::pow(b / M_PI, 1.5);

        SeparatedConvolution<double,3> G(world, coeff, expnt, lo, thresh);
        real_function_3d ff = G(f);
        double err2 = ff.err(NormalizedGaussian(a*b/(a+b)));
        
        auto op = [&](const Key<3>& key, FunctionNode<double,3>& node) {
            if (node.has_coeff()) {
                Level n = key.level();
                Tensor<double> M = make_T00(k, n, b*L*L*4.0);
                node.coeff() = transform(node.coeff(),transpose(M));
            }
        };

        int extension_order=0;

        auto op2 = [&](const Key<3>& key, FunctionNode<double,3>& node) {
            if (node.has_coeff()) {
                Level n = key.level();
                Tensor<double> Pp, Pm;
                make_taylor_extension_operatorsX(k, extension_order, Pp, Pm);
                std::shared_ptr<Convolution1D<double>> op = G.get_ops()[0].getop(0); // relying on cubic cell in test code
                Tensor<double> r0 = op->rnlij(n, 0);
                Tensor<double> rm = op->rnlij(n,-1);
                Tensor<double> rp = op->rnlij(n, 1);
                Tensor<double> R0 = inner(rp, Pm) + r0 + inner(rm, Pp);
                node.coeff() = transform(node.coeff(),transpose(R0));
            }
        };

        f.reconstruct();
        ff = copy(f);
        ff.unaryop_node(op);
        double err2M = ff.err(NormalizedGaussian(a*b/(a+b)));
        printf("  %.4e   %.4e     %.2e       %.2e", userh/sigma, b, err2, err2M);
        for (extension_order=2; extension_order<k; extension_order+=2) {
            ff = copy(f);
            ff.unaryop_node(op2);
            err2 = ff.err(NormalizedGaussian(a*b/(a+b)));
            printf("       %.2e", err2);
        }
        printf("\n");
    }
}


// Example usage and verification
int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    startup(world, argc, argv);
    test0();
    testD2();
    test1();
    test2();
    test3();
    test4(world);
    test5(world);
    test6(world);

    return 0;
}
