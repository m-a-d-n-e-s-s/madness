#include "testpc.h"
#include "madness/tensor/clapack_fortran.h"
#include "madness/world/madness_exception.h"

double (*f)(double, double, double) = nullptr;
double (*exact)(double, double, double) = nullptr;

// Evaluate erf(a*r)/r with correct limit at origin
double erfaroverr(double a, double r) {
    if (a*r*r/3.0 < 1e-8) {
        return 2*a/std::sqrt(pi);
    }
    else {
        return std::erf(a*r)/r;
    }
}

// Solve M[m,n] c[n] = f[m] in least squares sense
// M is the fitting matrix corresponding to m observations and n parameters
// f is the vector of m observations
// return the vector c of n fitted parameters
std::vector<double> solve(size_t m, size_t n, const std::vector<double>& M, const std::vector<double>& f) {
    std::vector<double> c = f;
    std::vector<double> M1 = M;
    integer lwork = 2 * std::max(m,n);
    std::vector<double> work(lwork);
    integer mm = m, nn = n, nrhs = 1, lda = n, ldb = n, info = 0;
    char trans = 'T';
    dgels_(&trans,&nn,&mm,&nrhs,M1.data(),&lda,c.data(),&ldb,work.data(),&lwork,&info,1);
    MADNESS_ASSERT(info == 0);
    
    return std::vector<double>(c.data(), c.data()+n);
}

// Fit the m data points in f to a polynomial of order n in inverse powers of N
std::vector<double> fit(size_t m, size_t n, const std::vector<double> N, const std::vector<double>& f) {
    //print(1,n,N);
    //print(m,1,f);
    
    std::vector<double> M;
    for (size_t i=0; i<m; i++) {
        for (size_t j=0; j<n; j++) {
            M.push_back(pow(1.0/N[i], j));
        }
    }
    //print(m,n,M);
    return solve(m, n, M, f);
}

// A C++ procedure to tabulate f(x,y,z) at x = x(i), y = y(j), z = z(k)
std::vector<double> tabulate(double(*f)(double, double, double), std::vector<double> x, std::vector<double> y, std::vector<double> z) {
    const size_t nx = x.size();
    const size_t ny = y.size();
    const size_t nz = z.size();
    std::vector<double> F(nx * ny * nz);
    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t k = 0; k < nz; k++) {
                F[i*ny*nz + j*nz + k] = f(x[i], y[j], z[k]);
            }
        }
    }
    return F;
}

// Periodic sum of functions
double periodic_sum(double x, double y, double z, double(*f)(double, double, double), const int N = 50) {
    // N number of lattice points summed in each direction
    double sum = 0.0;
    for (int X=-N; X<=N; X++) {
        for (int Y=-N; Y<=N; Y++) {
            for (int Z=-N; Z<=N; Z++) {
                sum += f(x+X*L,y+Y*L,z+Z*L);
            }
        }
    }
    return sum;
}

// Periodic sum of functions returning vector of values of partial sums in order of increasing cube size
std::vector<double> periodic_sum_partial(double x, double y, double z, double(*f)(double, double, double), const int N) {
    std::vector<double> sums;
    double sum = f(x,y,z);
    sums.push_back(sum);
    int count = 1;
    // Dumb loop structure since attempt to be smart failed!
    for (int n = 1; n <= N; n++) {
        for (int X=-n; X<=n; X++) {
            for (int Y=-n; Y<=n; Y++) {
                for (int Z=-n; Z<=n; Z++) {
                    if (std::abs(X) == n || std::abs(Y) == n || std::abs(Z) == n) {
                        sum += f(x+X*L,y+Y*L,z+Z*L);
                        count++;
                    }
                }
            }
        }
        sums.push_back(sum);
    }
    return sums;
}

// Periodic sum of functions accelerated by fitting t to a polynomial in 1/N
double periodic_sum_accelerated(double x, double y, double z, double(*f)(double, double, double), const size_t N=9, const size_t p=7) {
    std::vector<double> sums = periodic_sum_partial(x,y,z,f,N);
    // Extract the last p terms of the series in v
    std::vector<double> v(&sums[N-p], &sums[N]);
    std::vector<double> n(p);
    for (size_t j=0; j<p; j++) {
        n[j] = N-p+j;
    }
    std::vector<double> c = fit(p, p, n, v);
    return c[0];
}

double distancesq(double x, double y, double z, double x0, double y0, double z0) {
    double dx = x - x0;
    double dy = y - y0;
    double dz = z - z0;
    return dx*dx + dy*dy + dz*dz;
}

// Gaussian spheropole test function
const double a = 100.0;
const double b = 200.0;
double f_spheropole(double x, double y, double z) {
    const double afac = std::pow(a/pi, 1.5);
    const double bfac = std::pow(b/pi, 1.5);
    const double rsq = distancesq(x,y,z,xshift,yshift,zshift);
    return afac*exp(-a*rsq) - bfac*exp(-b*rsq);
}

// No need for periodic summation since the potential is zero exterior to the charge density
double exact_spheropole(double x, double y, double z) {
    const double mean = pi*(1/b - 1/a)/(L*L*L); // the mean of the potential over the domain
    const double rsq = distancesq(x,y,z,xshift,yshift,zshift);
    const double r = std::sqrt(rsq);
    //return std::erf(std::sqrt(a)*r)/r - std::erf(std::sqrt(b)*r)/r - mean;
    return erfaroverr(std::sqrt(a),r) - erfaroverr(std::sqrt(b),r) - mean;
};

// Gaussian dipole test function
const double offset = 0.1;
double f_dipole(double x, double y, double z) {
    const double bfac = std::pow(b/pi, 1.5);
    const double r1sq = distancesq(x,y,z,xshift,yshift,zshift-offset);
    const double r2sq = distancesq(x,y,z,xshift,yshift,zshift+offset);
    return bfac*(std::exp(-b*r1sq) - std::exp(-b*r2sq));
}

// Potential due to single Gaussian dipole
double exact_dipoleX(double x, double y, double z) {
    const double r1sq = distancesq(x,y,z,xshift,yshift,zshift-offset);
    const double r2sq = distancesq(x,y,z,xshift,yshift,zshift+offset);
    const double r1 = std::sqrt(r1sq);
    const double r2 = std::sqrt(r2sq);
    //return std::erf(std::sqrt(b)*r1)/r1 - std::erf(std::sqrt(b)*r2)/r2;
    double sb = std::sqrt(b);
    return erfaroverr(sb,r1) - erfaroverr(sb,r2);
};

// Potential due to opposing electric field generated by FT to satisty the periodic boundary conditions and continuity
double opposing_field_potential(double x, double y, double z) {
    // center of charge is at (xshift,yshift,zshift)
    const double mu = 2*offset;
    return mu*4*pi/(3*L*L*L) * (z-zshift);
}

// Periodic sum of dipole potentials
double exact_dipole(double x, double y, double z) {
    // const double mean = 0.0; // the mean of the potential over the domain(zero due to dipole symmetry)
    return periodic_sum_accelerated(x, y, z, exact_dipoleX) + opposing_field_potential(x,y,z);
}

// Gaussian quadrupole test function in yz plane
double f_quadrupole(double x, double y, double z) {
    const double bfac = std::pow(b/pi, 1.5);
    const double r1sq = distancesq(x, y, z, xshift, yshift,        zshift-offset);
    const double r2sq = distancesq(x, y, z, xshift, yshift,        zshift+offset);
    const double r3sq = distancesq(x, y, z, xshift, yshift-offset, zshift);
    const double r4sq = distancesq(x, y, z, xshift, yshift+offset, zshift);
    return bfac*(std::exp(-b*r1sq) + std::exp(-b*r2sq) - std::exp(-b*r3sq) - std::exp(-b*r4sq)  );
}

// Potential due to single Gaussian quadrupole
double exact_quadrupoleX(double x, double y, double z) {
    const double r1sq = distancesq(x, y, z, xshift, yshift,        zshift-offset);
    const double r2sq = distancesq(x, y, z, xshift, yshift,        zshift+offset);
    const double r3sq = distancesq(x, y, z, xshift, yshift-offset, zshift);
    const double r4sq = distancesq(x, y, z, xshift, yshift+offset, zshift);
    const double r1 = std::sqrt(r1sq);
    const double r2 = std::sqrt(r2sq);
    const double r3 = std::sqrt(r3sq);
    const double r4 = std::sqrt(r4sq);
    //return std::erf(std::sqrt(b)*r1)/r1 + std::erf(std::sqrt(b)*r2)/r2  - std::erf(std::sqrt(b)*r3)/r3  - std::erf(std::sqrt(b)*r4)/r4;
    double sb = std::sqrt(b);
    return  erfaroverr(sb,r1) + erfaroverr(sb,r2) - erfaroverr(sb,r3) - erfaroverr(sb,r4);
};

// Periodic sum of quadrupole potentials
double exact_quadrupole(double x, double y, double z) {
    // const double mean = 0.0; // the mean of the potential over the domain(zero due to dipole symmetry)
    //return periodic_sum(x, y, z, exact_quadrupoleX);
    return periodic_sum_accelerated(x, y, z, exact_quadrupoleX);
}


// Cosine test function
const double wx = 1;
const double wy = 2;
const double wz = 3;
double f_cosine(double x, double y, double z) {
    return std::cos(wx*2*pi*(x-xshift)/L)*std::cos(wy*2*pi*(y-yshift)/L)*std::cos(wz*2*pi*(z-zshift)/L);
}

double exact_cosine(double x, double y, double z) {
    return f_cosine(x,y,z) / (pi*(wx*wx + wy*wy + wz*wz)/(L*L));
}

// Return a vector with n equally spaced values between a and b, optionally including the right endpoint
std::vector<double> linspace(double a, double b, size_t n, bool include_right_endpoint) {
    double h = (b - a) / (include_right_endpoint ? (n - 1) : n);
    std::vector<double> v(n);
    for (size_t i = 0; i < n; i++) {
        v[i] = a + i*h;
    }
    return v;
}

void set_test_case(int test_case) {
    if (test_case == 1) {
        f = f_spheropole;
        exact = exact_spheropole;
        std::cout << "Test case set to Gaussian spheropole" << std::endl;
    } else if (test_case == 2) {
        f = f_dipole;
        exact = exact_dipole;
        std::cout << "Test case set to Gaussian dipole" << std::endl;
    } else if (test_case == 3) {
        f = f_quadrupole;
        exact = exact_quadrupole;
        std::cout << "Test case set to Gaussian quadrupole" << std::endl;
    } else if (test_case == 4) {
        f = f_cosine;
        exact = exact_cosine;
        std::cout << "Test case set to cosine" << std::endl;
    } else {
        std::cerr << "Invalid test case number" << std::endl;
        throw "bad";
    }
}
