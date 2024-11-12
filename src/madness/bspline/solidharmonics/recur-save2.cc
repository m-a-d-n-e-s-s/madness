#include <cmath>
#include <iostream>
#include <fstream>
#include <tuple>
#include <utility>
#include <string>
#include <cstdio>

// Solid harmonics following JCP 104, p8003 (1996)

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

static const double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816;

// sum(2*l+1,l=0..L) = (L+1)^2
// 20!  is largest factorial representable in 64-bit size_t
// 170! is largest factorial representable in double without overflowing the exponent

class LebedevQuadrature {
    size_t order;
    size_t npts;
    std::vector<double> thetas;
    std::vector<double> phis;
    std::vector<double> sinthetas;
    std::vector<double> costhetas;
    std::vector<double> sinphis;
    std::vector<double> cosphis;
    std::vector<double> weights;

public:

    // From https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html
    // Adjust input values by
    // (1) convert to radian
    // (2) switch theta and phi (physics convention)
    // (3) quadrature range phi=0..2pi (instead of -pi..pi)
    // (4) incorporate factor of 4*pi in the weights
    LebedevQuadrature(int maxl) {
        for (auto n : {3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,35,41,47,53,59,65,71,77,83,89,95,101,113,119,125,131}) {
            if (n >= maxl) {
                order = n;
                char filename[16];
                sprintf(filename,"leb_%03d.txt",n);
                std::ifstream f(filename);
                
                f >> npts;
                //std::cout << filename << " " << npts << std::endl;
                
                thetas.resize(npts);
                phis.resize(npts);
                sinthetas.resize(npts);
                costhetas.resize(npts);
                sinphis.resize(npts);
                cosphis.resize(npts);
                weights.resize(npts);

                double wsum = 0.0;
                
                double deg_to_rad = pi/180.0;
                for (size_t i=0; i<npts; i++) {
                    f >> phis[i] >> thetas[i] >> weights[i];
                    phis[i] += 180.0;
                    phis[i] = phis[i]*deg_to_rad;
                    thetas[i] = thetas[i]*deg_to_rad;
                    weights[i] *= 4*pi;

                    sinthetas[i] = std::sin(thetas[i]);
                    costhetas[i] = std::cos(thetas[i]);
                    sinphis[i] = std::sin(phis[i]);
                    cosphis[i] = std::cos(phis[i]);
                    
                    //std::cout << i << " " << thetas[i] << " " << phis[i] << " " << weights[i] << std::endl;
                    wsum += weights[i];
                }
                //std::cout << "wsum " << (wsum-4*pi)/(4*pi) << std::endl;
                return;
            }
        }
        throw "maximum supported angular momentum is 131";
    }

    template <typename functionT>
    auto cartesian_integral(functionT f, double r=1.0) const {
        auto value = f(r,0.0,0.0); // just to deduce the return type using auto
        value = 0;
        for (size_t i=0; i<npts; i++) {
            double x = r*sinthetas[i]*cosphis[i];
            double y = r*sinthetas[i]*sinphis[i];
            double z = r*costhetas[i];
            double weight = weights[i];
            value += f(x,y,z) * weight;
        }
        value *= r*r;
        return value;
    }

    template <typename functionT>
    auto spherical_integral(functionT f, double r=1.0) const {
        auto value = f(0.0,0.0,0.0); // just for type deduction
        value = 0;
        for (size_t i=0; i<npts; i++) {
            double theta = thetas[i];
            double phi = phis[i];
            double weight = weights[i];
            value += weight * f(r,theta,phi);
        }
        return value * r * r;
    }

    // Returns cartesian points on unit sphere
    std::tuple<std::vector<double>,std::vector<double>,std::vector<double>>
    cartesian_pts() const {
        std::vector<double> x(npts);
        std::vector<double> y(npts);
        std::vector<double> z(npts);
        for (size_t i=0; i<npts; i++) {
            x[i] = sinthetas[i]*cosphis[i];
            y[i] = sinthetas[i]*sinphis[i];
            z[i] = costhetas[i];
        }
        return {x,y,z};
    }

    std::vector<double> wts() const {
        return weights;
    }
};    

double factorial(double a) {
    double value = 1.0;
    while (a>0) {value*=a--;}
    return value;
}

// // given l>=0 and -l<=m<=l compute index into packed array for l, m
// constexpr size_t lm_index(unsigned l, int m) {
//     return l*(l+1) + m;
// }

double Alm(unsigned l, int m) {
    double top = pi * factorial(l+m) * factorial(l-m);
    double bot = 2*l+1;
    return std::sqrt(top/bot);
}

class NlmContainer {
    unsigned LMAX;
    std::vector<double> N;

    // The following are layed out contiguously 
    // N+[0,0]
    // N-[1,1] N+[1,0] N+[1,1]
    // N-[2,2] N-[2,1] N+[2,0] N+[2,1] N+[2,2]
    // etc.
    // (N-[l,m] for m=l,l-1,..,1) and then (N+[l,m] for m=0..l)

    // Returns offset to data for l so that N(loff(l)+m) yields the desired value for -l<=m<=l
    size_t loff(unsigned l) const {
        return l*(l+1);
    }
public:
    NlmContainer(unsigned lmax)
        : LMAX(lmax)
        , N((lmax+1)*(lmax+1),0.0) // must initialize everything so that unused upper elements are valid for ops like addition
    {}

    unsigned lmax() const {return LMAX;}

    // Returns reference to N[l,m] with m>=0 being N+ and m<0 being N-
    double& value(unsigned l, int m) {
        return N[loff(l)+m];
    }

    // Returns value N[l,m] with m>=0 being N+ and m<0 being N-
    double value(unsigned l, int m) const {
        return N[loff(l)+m];
    }

    // Fill with scalar
    NlmContainer& operator=(double s) {
        for (size_t i=0; i<N.size(); i++) N[i] = s;
        return *this;
    }

    NlmContainer& operator+=(const NlmContainer& a) {
        if (LMAX != a.LMAX) throw "sizes must match";
        for (size_t i=0; i<N.size(); i++) N[i] += a.N[i];
        return *this;
    }

    NlmContainer operator+(const NlmContainer& a) const {
        if (LMAX != a.LMAX) throw "sizes must match";
        NlmContainer b(LMAX);
        for (size_t i=0; i<N.size(); i++) b.N[i] = N[i] + a.N[i];
        return b;
    }

    NlmContainer& operator-=(const NlmContainer& a) {
        if (LMAX != a.LMAX) throw "sizes must match";
        for (size_t i=0; i<N.size(); i++) N[i] -= a.N[i];
        return *this;
    }

    NlmContainer operator-(const NlmContainer& a) const {
        if (LMAX != a.LMAX) throw "sizes must match";
        NlmContainer b(LMAX);
        for (size_t i=0; i<N.size(); i++) b.N[i] = N[i] - a.N[i];
        return b;
    }

    NlmContainer& operator*=(const NlmContainer& a) {
        if (LMAX != a.LMAX) throw "sizes must match";
        for (size_t i=0; i<N.size(); i++) N[i] *= a.N[i];
        return *this;
    }

    NlmContainer operator*(const NlmContainer& a) const {
        if (LMAX != a.LMAX) throw "sizes must match";
        NlmContainer b(LMAX);
        for (size_t i=0; i<N.size(); i++) b.N[i] = N[i] * a.N[i];
        return b;
    }

    NlmContainer& operator*=(double a) {
        for (size_t i=0; i<N.size(); i++) N[i] *= a;
        return *this;
    }

    NlmContainer operator*(double a) const {
        NlmContainer b(LMAX);
        for (size_t i=0; i<N.size(); i++) b.N[i] = N[i]*a;
        return b;
    }

    static void test() {
        int lmax = 20;
        NlmContainer N(lmax);
        for (int l=0; l<lmax; l++) {
            for (int m=-l; m<=l; m++) {
                N.value(l,m) = 101*l + m;
            }
        }

        for (int l=0; l<lmax; l++) {
            for (int m=-l; m<=l; m++) {
                if (N.value(l,m) != 101*l + m) throw "NlmContainer test failed";
            }
        }
    }
};


std::ostream& operator<<(std::ostream& s, const NlmContainer& Nlm) {
    char buf[256];
    for (int l=0; l<=int(Nlm.lmax()); l++) {
        for (int m=-l; m<=l; m++) {
            sprintf(buf,"%12.4e",Nlm.value(l,m));
            s << buf;
        }
        std::cout << std::endl;
    }
    return s;
}


// Compute N+ and N- for given (l,|m|) which is not efficient if you want many values
std::pair<double,double> Nlm(unsigned l, unsigned m, double x, double y, double z) {
    double rsq = x*x + y*y + z*z;

    // First recur m up from 0 to m to make N+|m||m| and N-|m||m|
    double Np_mm=1.0, Nm_mm=0.0;
    for (unsigned i=1; i<=m; i++) { // i plays role of m
        double Np = (-1.0/(2*i))*(x*Np_mm - y*Nm_mm);
        double Nm = (-1.0/(2*i))*(y*Np_mm + x*Nm_mm);
        Np_mm = Np;
        Nm_mm = Nm;
    }
    //std::cout << "after mm " << Np_mm << " " << Nm_mm << std::endl;

    // Now recur l up from m to l to make N+l|m| and N-l|m|
    // Note that Nlm is treated as zero if l<m, so we start the recurrence
    // with l-1=m and hence the l-2 term starts as zero.
    double Np_l1m=Np_mm, Np_l2m=0.0; // l-1,m and l-2,m
    double Nm_l1m=Nm_mm, Nm_l2m=0.0;
    for (unsigned i=m+1; i<=l; i++) { // i plays role of l
        double fac = 1.0/((i+m)*(i-m));
        double Np = fac*((2*i-1)*z*Np_l1m - rsq*Np_l2m);
        double Nm = fac*((2*i-1)*z*Nm_l1m - rsq*Nm_l2m);
        Np_l2m = Np_l1m;
        Np_l1m = Np;
        Nm_l2m = Nm_l1m;
        Nm_l1m = Nm;
    }
    return {Np_l1m, Nm_l1m};
}


// Compute Nlm for all l and m up to lmax
NlmContainer Nlm(unsigned lmax, double x, double y, double z) {
    double rsq = x*x + y*y + z*z;

    NlmContainer N(lmax);

    // Initialize m=0 for which N- are missing
    N.value(0,0) = 1;
    double Np_l1m=1, Np_l2m=0.0; // l-1,m and l-2,m
    for (unsigned l=1; l<=lmax; l++) {
        double fac = 1.0/(l*l);
        double Np = fac*((2*l-1)*z*Np_l1m - rsq*Np_l2m);
        Np_l2m = Np_l1m;
        Np_l1m = Np;
        N.value(l,0) = Np;
    }
    
    // Recur m up to make N+|m||m| and N-|m||m|
    for (unsigned m=1; m<=lmax; m++) {
        double Np1 = N.value(m-1,m-1);
        double Nm1 = (m<=1) ? 0 : N.value(m-1,1-m);
        N.value(m, m) = (-1.0/(2*m))*(x*Np1 - y*Nm1);
        N.value(m,-m) = (-1.0/(2*m))*(y*Np1 + x*Nm1);
        
        // Recur l up from m to lmax to make N+l|m| and N-l|m|
        // Note that Nlm is treated as zero if l<m, so we start the recurrence
        // with l-1=m and hence the l-2 term starts as zero.
        double Np_l1m=N.value(m, m), Np_l2m=0.0; // l-1,m and l-2,m
        double Nm_l1m=N.value(m,-m), Nm_l2m=0.0;
        for (unsigned l=m+1; l<=lmax; l++) {
            double fac = 1.0/((l+m)*(l-m));
            double Np = fac*((2*l-1)*z*Np_l1m - rsq*Np_l2m);
            double Nm = fac*((2*l-1)*z*Nm_l1m - rsq*Nm_l2m);
            Np_l2m = Np_l1m;
            Np_l1m = Np;
            Nm_l2m = Nm_l1m;
            Nm_l1m = Nm;
            N.value(l, m) = Np;
            N.value(l,-m) = Nm;
        }
    }
    return N;
}

// Compute Nlm from their definition in terms of spherical harmonics.  Boost uses a different
// definition for spherical harmonics so some details differ from the original paper.
std::pair<double,double> Nlm_reference(unsigned l, unsigned m, double x, double y, double z) {
    double r = std::sqrt(x*x + y*y + z*z);
    double theta=0, phi=0;

    if (r > 1e-16) {
        theta = std::acos(z/r);
        phi = std::atan(y/x);
    }

    auto ylm = boost::math::spherical_harmonic(l, m, theta, phi);
    auto ylmm = boost::math::spherical_harmonic(l, -int(m), theta, phi);

    std::complex<double> I = {0.0,1.0};
    double mfac = std::pow(-1.0,m);
    std::complex<double> facp = I * mfac;
    std::complex<double> facm = 1 * -mfac;
            
    double NNp = (facp*pi/Alm(l,m)/(2.0*l+1)*std::pow(r,l)*(ylmm + mfac*ylm)).imag();
    double NNm = (facm*pi/Alm(l,m)/(2.0*l+1)*std::pow(r,l)*(ylmm - mfac*ylm)).imag();
    
    return {NNp, NNm};
}

double Nlm_normalization(unsigned l, int m) {
    double value = pi/(Alm(l,m)*(2*l+1));
    static const double sqrthalf = std::sqrt(0.5);
    if (m) return sqrthalf/value;
    else return 0.5/value;
}

NlmContainer Nlm_normalization(unsigned lmax) {
    NlmContainer Nnorm(lmax);
    for (unsigned l=0; l<=lmax; l++) {
        for (unsigned m=0; m<=l; m++) {
            double norm = Nlm_normalization(l, m);
            Nnorm.value(l,m) = norm;
            if (m > 0) Nnorm.value(l,-m) = norm;
        }
    }
    return Nnorm;
}

void test_Nlm() {
    {
        // Check the container internal storage logic
        NlmContainer::test();
    }
    
    {
        // Check the values computed in three different ways
        size_t lmax = 8*3;
        double x = 0.1;
        double y = 0.3;
        double z = -0.7;
        
        NlmContainer Nfull = Nlm(lmax, x, y, z);
        NlmContainer Nnorm = Nlm_normalization(lmax);
        NlmContainer N(lmax), Nref(lmax);
        for (size_t l=0; l<=lmax; l++) {
            for (size_t m=0; m<=l; m++) {
                {
                    auto [Np, Nm] = Nlm(l, m, x, y, z);
                    N.value(l,m) = Np;
                    if (m > 0) N.value(l,-m) = Nm;
                }
                {
                    auto [Np, Nm] = Nlm_reference(l, m, x, y, z);
                    Nref.value(l,m) = Np;
                    if (m > 0) Nref.value(l,-m) = Nm;
                }
            }
        }
        
        // Verify the normalized values
        N *= Nnorm;
        Nref *= Nnorm;
        Nfull *= Nnorm;
        
        NlmContainer err1 = Nfull-N;
        NlmContainer err2 = Nfull-Nref;
        double e1=0, e2=0;
        for (int l=0; l<int(lmax); l++) {
            for (int m=-l; m<=l; m++) {
                double d1 = err1.value(l,m), d2 = err2.value(l,m);
                e1 += d1*d1;
                e2 += d2*d2;
            }
        }
        e1 = std::sqrt(e1);
        e2 = std::sqrt(e2);
        if (e1 != 0.0) throw "e1 failed";
        if (e2 > 6e-16) throw "e2 failed";
        std::cout << "xxxx " << e1 << " " << e2 << std::endl;
        // std::cout << N << std::endl;
        // std::cout << Nref << std::endl;
        // std::cout << Nfull << std::endl;
        // std::cout << err1 << std::endl;
        // std::cout << err2 << std::endl;
    }
        
    // Verify orthonormality
    {
        size_t lmax = 8; // tried 21 and it was OK
        NlmContainer Nnorm = Nlm_normalization(lmax);
        LebedevQuadrature q(2*lmax);

        for (int l=0; l<=int(lmax); l++) {
            for (int m=-l; m<=l; m++) {
                auto f = [&](double x, double y, double z) {
                    NlmContainer N = Nlm(lmax, x, y, z);
                    double v = N.value(l,m);
                    v *= Nnorm.value(l,m);
                    N *= Nnorm;
                    return N*v;
                };
                auto N = q.cartesian_integral(f, 1.0);
                //std::cout << l << " " << m << std::endl << N << std::endl;
                for (int ll=0; ll<=int(lmax); ll++) {
                    for (int mm=-ll; mm<=ll; mm++) {
                        double v = N.value(ll,mm);
                        double err = std::abs(v);
                        if (ll==l && mm==m) err -= 1;
                        // <<<< clearly need to compute reference quantities in extended precision esp if taking sqrts
                        // 1) is this just lebedev fail? are we using the correct order?
                        //    - it's a lot worse for the 3j coeffs
                        // 2) or just the number of points in the sum along with cancellation?
                        // 3) or is the recursion unstable with the chosen normalization?
                        //    - Tests against boost ylm agree to close to machine precision for the normalized values
                        //    - If it is unstable, then we can precompute Nlm on the spherical and MRA quadrature grids in higher precision.
                        if (err > 2e-12) {
                            std::cout << "bad norm (" << l << "," << m << ")  (" << ll << "," << mm << ")  " << v << " " << err << std::endl;
                            std::cout << N << std::endl;
                            throw "bad norm";
                        }
                    }
                }
            }
        }
    }
}

// Production version will need extended precision and use of GaussLegendre instead of Lebedev
void threeJ() {
    int lmax1 = 8;
    int lmax2 = lmax1; // assumed equal below
    int lmax3 = lmax1 + lmax2; // Product order is sum of both
    
    NlmContainer Nnorm = Nlm_normalization(lmax3);
    LebedevQuadrature q(lmax3*2); // Quadrature is over products so double order

    std::cout << std::setprecision(6) << std::scientific;

    // Sparse storage for the 3j coeffs
    int lx4 = (lmax1+1)*(lmax1+1)*(lmax1+1)*(lmax1+1);
    std::vector<std::pair<int,int>> index(lx4); // first=index into j3+lms, second=count
    std::vector<std::tuple<int,int,double>> j3; // l3, m3, value
    j3.reserve(2*lx4);
    std::cout << "j3 size " << j3.size()  << std::endl;
    auto lmx = [=](int l1, int m1, int l2, int m2) { // offset into the index vector
        return (l1*(l1+1)+m1)*(lmax1+1)*(lmax1+1) + l2*(l2+1)+m2;
    };
        
    for (int l1=0; l1<=lmax1; l1++) {
        for (int m1=-l1; m1<=l1; m1++) {
            for (int l2=0; l2<=lmax2; l2++) {
                for (int m2=-l2; m2<=l2; m2++) {
                    
                    auto f = [&](double x, double y, double z) {
                        NlmContainer v3 = Nlm(lmax3, x, y, z);
                        double v1 = v3.value(l1,m1);
                        double v2 = v3.value(l2,m2);
                        v3 *= Nnorm; // To make f(x)-s*g(x) orthog to g(x) use s=<g|f>/<g|g>
                        v3 *= Nnorm;
                        return v3*(v2*v1);
                    };

                    auto C = q.cartesian_integral(f,1.0);
                    size_t off=j3.size();
                    for (int l3=0; l3<=lmax3; l3++) {
                        for (int m3=-l3; m3<=l3; m3++) {
                            // For screening small values need to normalize all functions
                            double test = Nnorm.value(l1,m1)*Nnorm.value(l2,m2)*Nnorm.value(l3,m3)*C.value(l3,m3);
                            if (std::abs(test) > 1e-8) { // <<<<<<<<<<<<<<<<<<<<<<<<< again need lots more precision ... why is it so bad?
                                //std::cout << std::setw(3) << l1 << " " << std::setw(3) << m1 << " " << std::setw(3) << l2 << " " << std::setw(3) << m2 << "   --> " << std::setw(3) << l3 << " " << std::setw(3) << m3 << " " << std::setw(15) << C.value(l3,m3) << " " << test << std::endl;
                                j3.push_back({l3,m3,C.value(l3,m3)});
                            }
                        }
                    }
                    //std::cout << off << " " << j3.size() << std::endl;
                    index[lmx(l1,m1,l2,m2)] = {off,int(j3.size()-off)};
                }
            }
        }
    }
    std::cout << "j3 size " << j3.size() << " " << lx4 << std::endl;

    // // Now test reconstructing Nlm from the coupling coefficients

    double x = 0.1;
    double y = 0.3;
    double z = -0.7;
    double r = std::sqrt(x*x + y*y + z*z);
    x /= r; // Put test points on the unit sphere
    y /= r;
    z /= r;
    r = 1;

    NlmContainer N = Nlm(lmax3, x, y, z);

    std::cout << N;

    double maxerr = 0.0, maxrelerr = 0.0;
    for (int l1=0; l1<=lmax1; l1++) {
        for (int m1=-l1; m1<=l1; m1++) {
            for (int l2=0; l2<=lmax2; l2++) {
                for (int m2=-l2; m2<=l2; m2++) {
                    auto [off,n] = index[lmx(l1,m1,l2,m2)];
                    double Ntest = 0.0;
                    for (int i=0; i<n; i++) {
                        auto [l3,m3,coeff] = j3[off+i];
                        Ntest += coeff*N.value(l3,m3);
                    }
                    double prod = N.value(l1,m1)*N.value(l2,m2);
                    double err =  prod - Ntest;
                    if (std::abs(err) > 5e-14) {
                        std::cout << std::setw(3) << l1 << " " << std::setw(3) << m1 << " " << std::setw(3) << l2 << " " << std::setw(3) << m2 << "  " << std::setw(15) << err << std::endl;
                    }
                    maxerr = std::max(err,maxerr);
                    maxrelerr = std::max(err/prod,maxrelerr);
                }
            }
        }
    }
    std::cout << "   maxerr in J3 reconstruction on unit sphere " << maxerr << std::endl;
    std::cout << "maxrelerr in J3 reconstruction on unit sphere " << maxrelerr << std::endl;
}


int main() {
    test_Nlm();
    
    threeJ();
    
    return 0;
}
