#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <tuple>
#include <utility>
#include <string>
#include <cstdio>
#include <gauleg.h>
#include <myqd.h>

// Solid harmonics following JCP 104, p8003 (1996)

template <typename T>
struct constants {
    static T pi;

    static void init() {
        pi = from_str<T>("3.141592653589793238462643383279502884197169399375105820974944592307816");
    }
};
template <> double constants<double>::pi = 0;
template <> dd_real constants<dd_real>::pi = 0;
template <> qd_real constants<qd_real>::pi = 0;


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

                double pi = constants<double>::pi;
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

template <typename T>
T factorial(size_t a) {
    T value = T(1);
    while (a>0) {value*=T(double(a--));}
    return value;
}

template <typename T>
T Alm(unsigned l, int m) {
    T top = constants<T>::pi * factorial<T>(l+m) * factorial<T>(l-m);
    T bot = 2.*l+1;
    return std::sqrt(top/bot);
}

template <typename T>
class NlmContainer {
    unsigned LMAX;
    std::vector<T> N;

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
        , N((lmax+1)*(lmax+1),T(0)) // must initialize everything so that unused upper elements are valid for ops like addition
    {}

    unsigned lmax() const {return LMAX;}

    // Returns reference to N[l,m] with m>=0 being N+ and m<0 being N-
    T& value(unsigned l, int m) {
        return N[loff(l)+m];
    }

    // Returns value N[l,m] with m>=0 being N+ and m<0 being N-
    T value(unsigned l, int m) const {
        return N[loff(l)+m];
    }

    // Fill with scalar
    NlmContainer& operator=(T s) {
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

    NlmContainer& operator*=(T a) {
        for (size_t i=0; i<N.size(); i++) N[i] *= a;
        return *this;
    }

    NlmContainer operator*(T a) const {
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

template <typename T>
std::ostream& operator<<(std::ostream& s, const NlmContainer<T>& Nlm) {
    char buf[256];
    for (int l=0; l<=int(Nlm.lmax()); l++) {
        for (int m=-l; m<=l; m++) {
            sprintf(buf,"%12.4e",to_double(Nlm.value(l,m)));
            s << buf;
        }
        std::cout << std::endl;
    }
    return s;
}


// Compute N+ and N- for given (l,|m|) which is not efficient if you want many values
template <typename T>
std::pair<T,T> Nlm(unsigned l, unsigned m, const T& x, const T& y, const T& z) {
    T rsq = x*x + y*y + z*z;

    // First recur m up from 0 to m to make N+|m||m| and N-|m||m|
    T Np_mm=1, Nm_mm=0;
    for (unsigned i=1; i<=m; i++) { // i plays role of m
        T fac = T(-1)/T(2.*i);
        T Np = fac*(x*Np_mm - y*Nm_mm);
        T Nm = fac*(y*Np_mm + x*Nm_mm);
        Np_mm = Np;
        Nm_mm = Nm;
    }
    //std::cout << "after mm " << Np_mm << " " << Nm_mm << std::endl;

    // Now recur l up from m to l to make N+l|m| and N-l|m|
    // Note that Nlm is treated as zero if l<m, so we start the recurrence
    // with l-1=m and hence the l-2 term starts as zero.
    T Np_l1m=Np_mm, Np_l2m=T(0); // l-1,m and l-2,m
    T Nm_l1m=Nm_mm, Nm_l2m=T(0);
    for (unsigned i=m+1; i<=l; i++) { // i plays role of l
        T fac = T(1.)/T(double((i+m)*(i-m)));
        T Np = fac*(T(2.*i-1)*z*Np_l1m - rsq*Np_l2m);
        T Nm = fac*(T(2.*i-1)*z*Nm_l1m - rsq*Nm_l2m);
        Np_l2m = Np_l1m;
        Np_l1m = Np;
        Nm_l2m = Nm_l1m;
        Nm_l1m = Nm;
    }
    return {Np_l1m, Nm_l1m};
}


// Compute Nlm for all l and m up to lmax
template <typename T>
NlmContainer<T> Nlm(unsigned lmax, T x, T y, T z) {
    T rsq = x*x + y*y + z*z;

    NlmContainer<T> N(lmax);

    // Initialize m=0 for which N- are missing
    N.value(0,0) = T(1);
    T Np_l1m=T(1), Np_l2m=T(0); // l-1,m and l-2,m
    for (unsigned l=1; l<=lmax; l++) {
        T fac = T(1)/T(int(l*l));
        T Np = fac*(T(2*l-1.0)*z*Np_l1m - rsq*Np_l2m);
        Np_l2m = Np_l1m;
        Np_l1m = Np;
        N.value(l,0) = Np;
    }
    
    // Recur m up to make N+|m||m| and N-|m||m|
    for (unsigned m=1; m<=lmax; m++) {
        T Np1 = N.value(m-1,m-1);
        T Nm1 = (m<=1) ? T(0) : N.value(m-1,1-m);
        T fac = T(-1)/T(2.*m);        
        N.value(m, m) = fac*(x*Np1 - y*Nm1);
        N.value(m,-m) = fac*(y*Np1 + x*Nm1);
        
        // Recur l up from m to lmax to make N+l|m| and N-l|m|
        // Note that Nlm is treated as zero if l<m, so we start the recurrence
        // with l-1=m and hence the l-2 term starts as zero.
        T Np_l1m=N.value(m, m), Np_l2m=T(0); // l-1,m and l-2,m
        T Nm_l1m=N.value(m,-m), Nm_l2m=T(0);
        for (unsigned l=m+1; l<=lmax; l++) {
            T fac = T(1)/T(1.0*(l+m)*(l-m));
            T Np = fac*(T(2.*l-1)*z*Np_l1m - rsq*Np_l2m);
            T Nm = fac*(T(2.*l-1)*z*Nm_l1m - rsq*Nm_l2m);
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

template <typename T>
T Nlm_normalization(unsigned l, int m) {
    T value = constants<T>::pi/(Alm<T>(l,m)*T(2.*l+1));
    static const T sqrthalf = std::sqrt(T(0.5));
    if (m) return sqrthalf/value;
    else return T(0.5)/value;
}

template <typename T>
NlmContainer<T> Nlm_normalization(unsigned lmax) {
    NlmContainer<T> Nnorm(lmax);
    for (unsigned l=0; l<=lmax; l++) {
        for (unsigned m=0; m<=l; m++) {
            T norm = Nlm_normalization<T>(l, m);
            Nnorm.value(l,m) = norm;
            if (m > 0) Nnorm.value(l,-m) = norm;
        }
    }
    return Nnorm;
}

template <typename T, typename R=qd_real>
void test_Nlm() {
    {
        // Check the container internal storage logic
        NlmContainer<T>::test();
    }
    
    {
        // Check the values computed in three different ways
        size_t lmax = 8;
        T x = 0.125;
        T y = 0.325;
        T z = -0.75;
        
        NlmContainer<T> Nfull = Nlm(lmax, x, y, z);
        NlmContainer<T> Nnorm = Nlm_normalization<T>(lmax);
        NlmContainer<T> N(lmax);

        NlmContainer<R> Nref = Nlm(lmax, R(x), R(y), R(z));
        NlmContainer<R> NnormR = Nlm_normalization<R>(lmax);
        
        for (size_t l=0; l<=lmax; l++) {
            for (size_t m=0; m<=l; m++) {
                {
                    auto [Np, Nm] = Nlm(l, m, x, y, z);
                    N.value(l,m) = Np;
                    if (m > 0) N.value(l,-m) = Nm;
                }
            }
        }
        
        // Verify the normalized values
        N *= Nnorm;
        Nfull *= Nnorm;
        Nref *= NnormR;
        
        NlmContainer<T> err1 = Nfull-N;
        NlmContainer<T> err2(lmax); ///  = Nfull-Nref; don't have type conversion on containers yet
        
        T e1=0, e2=0;
        T maxreferr=0, maxrefrelerr=0;
        for (int l=0; l<=int(lmax); l++) {
            for (int m=-l; m<=l; m++) {
                err2.value(l,m) = convert<T>(R(Nfull.value(l,m))-Nref.value(l,m));
                T d1 = err1.value(l,m), d2 = err2.value(l,m);
                maxreferr = std::max(std::abs(d2),maxreferr);
                maxrefrelerr = std::max(std::abs(d2)/std::abs(Nfull.value(l,m)),maxrefrelerr);
                e1 += d1*d1;
                e2 += d2*d2;
            }
        }
        e1 = std::sqrt(e1);
        e2 = std::sqrt(e2);
        if (e1 != 0) throw "e1 failed";
        if (e2 > 6e-16) throw "e2 failed";
        std::cout << "xxxx " << e1 << " " << e2 << std::endl;
        std::cout << "yyyy " << maxreferr << " " << maxrefrelerr << std::endl;
        std::cout << N << std::endl;
        std::cout << Nref << std::endl;
        std::cout << Nfull << std::endl;
        std::cout << err1 << std::endl;
        std::cout << err2 << std::endl;
    }
        
    // Verify orthonormality
    {
        size_t lmax = 8;
        NlmContainer<T> Nnorm = Nlm_normalization<T>(lmax);
        GaussLegendreSphere<T> q(2*lmax);

        T maxerr = 0;
        for (int l=0; l<=int(lmax); l++) {
            for (int m=-l; m<=l; m++) {
                auto f = [&](T x, T y, T z) {
                    NlmContainer<T> N = Nlm(lmax, x, y, z);
                    T v = N.value(l,m);
                    v *= Nnorm.value(l,m);
                    N *= Nnorm;
                    return N*v;
                };
                auto N = q.cartesian_integral(f, T(1.0));
                //std::cout << l << " " << m << std::endl << N << std::endl;
                for (int ll=0; ll<=int(lmax); ll++) {
                    for (int mm=-ll; mm<=ll; mm++) {
                        T v = N.value(ll,mm);
                        T err = std::abs(v);
                        if (ll==l && mm==m) err -= 1;
                        maxerr =std::max(maxerr,err);
                        if (err > T(2e-12)) { // <<<<<<<< change to numeric_limit for data type
                            std::cout << "bad norm (" << l << "," << m << ")  (" << ll << "," << mm << ")  " << v << " " << err << std::endl;
                            std::cout << N << std::endl;
                            throw "bad norm";
                        }
                    }
                }
            }
        }
        std::cout << "maxerr in orthonormality test " << maxerr << std::endl;
    }
}

template <typename T>
auto threeJ() {
    int lmax1 = 16;
    int lmax2 = lmax1; // assumed equal below
    int lmax3 = lmax1 + lmax2; // Product order is sum of both
    
    std::vector<GaussLegendreSphere<T>> qs; // Use lowest order rule possible
    for (int l=0; l<=2*lmax3; l++) qs.emplace_back(GaussLegendreSphere<T>(l));

    //std::vector<LebedevQuadrature> qs; // Use lowest order rule possible
    //for (int l=0; l<=2*lmax3; l++) qs.emplace_back(LebedevQuadrature(l));

    std::cout << std::setprecision(6) << std::scientific;

    // Sparse storage for the 3j coeffs
    int lx4 = (lmax1+1)*(lmax1+1)*(lmax1+1)*(lmax1+1);
    std::vector<std::pair<int,int>> index(lx4); // first=index into j3, second=count
    std::vector<std::tuple<int,int,T>> j3; // l3, m3, value
    j3.reserve((lmax1*lx4)/2); // seems to be a v close upper estimate of the size
    
    auto lmx = [=](int l1, int m1, int l2, int m2) { // offset into the index vector
        return (l1*(l1+1)+m1)*(lmax1+1)*(lmax1+1) + l2*(l2+1)+m2;
    };


    // N(l1,m1)*N(l2,m2) = sum(l3,m3) j(l1,m1,l2,m2,l3,m3) N(l3,m3)
    // int(N(l1,m1)*N(l2,m2)*N(l3,m3)) = j(l1,m1,l2,m2,l3,m3) int(N(l3,m3)**2) by orthogonality
    // --> j(l1,m1,l2,m2,l3,m3) = int(N(l1,m1)*N(l2,m2)*N(l3,m3)) * Nnorm(l3,m3)**2
    
    // But the function norms can be massive, so for screening small
    // coeffs imagine using and projecting with the normalized
    // functions
    //
    // N(l1,m1) N(l2,m2) Norm(l1,m1) Norm(l2,m2) = sum(l3,m3) j(l1,m1,l2,m2,l3,m3) N(l3,m3) Norm(l1,m1) Norm(l2,m2) Norm(l3,m3)
    // int(N(l1,m1)*N(l2,m2)*N(l3,m3)) Norm(l1,m1) Norm(l2,m2) Norm(l3,m3) = j(l1,m1,l2,m2,l3,m3) Norm(l1,m1) Norm(l2,m2) / Norm(l3,m3)
    // so the normalized 3j coeffs are
    // Jnormalized(l1,m1,l2,m2,l3,m3) = j(l1,m1,l2,m2,l3,m3) Norm(l1,m1) Norm(l2,m2) / Norm(l3,m3)
    // If this is small relative to unity, we can discard values.
        
    for (int l1=0; l1<=lmax1; l1++) {
        for (int m1=-l1; m1<=l1; m1++) {
            for (int l2=0; l2<=l1; l2++) { // <<<<<<<<<<<<<< restrict range to save space
                for (int m2=-l2; m2<=l2; m2++) {
                    int lmax4 = l1+l2;
                    NlmContainer<T> Nnorm = Nlm_normalization<T>(lmax4);
                    auto f = [&](T x, T y, T z) {
                        NlmContainer<T> v3 = Nlm(lmax4, x, y, z);
                        T v1 = v3.value(l1,m1);
                        T v2 = v3.value(l2,m2);
                        v3 *= Nnorm;
                        v3 *= Nnorm;
                        return v3*(v2*v1);
                    };

                    auto C = qs[2*lmax4].cartesian_integral(f,1.0); // use lowest order possible
                    
                    size_t off=j3.size();
                    for (int l3=std::abs(l1-l2); l3<=l1+l2; l3++) { // enforce known restrictions
                        for (int m3=-l3; m3<=l3; m3++) {
                            // For screening small values need to include normalization of l1 and l2
                            // But need to partially undo that for l3 due to projection onto unnormalized functions
                            T j3value = C.value(l3,m3);
                            T test = j3value*Nnorm.value(l1,m1)*Nnorm.value(l2,m2)/Nnorm.value(l3,m3);
                            if (std::abs(test) > 1e-8) {
                                //std::cout << std::setw(3) << l1 << " " << std::setw(3) << m1 << " " << std::setw(3) << l2 << " " << std::setw(3) << m2 << "   --> " << std::setw(3) << l3 << " " << std::setw(3) << m3 << " " << std::setw(15) << j3value << " " << test << std::endl;
                                j3.push_back({l3,m3,j3value});
                            }
                        }
                    }
                    index[lmx(l1,m1,l2,m2)] = {off,int(j3.size()-off)};
                    index[lmx(l2,m2,l1,m1)] = {off,int(j3.size()-off)};
                }
            }
        }
    }
    std::cout << "j3 size " << j3.size() << " " << (lmax1*lx4)/2 << std::endl;

    return std::make_tuple(lmax1,index,j3);
}

template <typename T>
void test_j3(int lmax1, const std::vector<std::pair<int,int>>& index, const std::vector<std::tuple<int,int,T>>& j3) {
    // Test reconstructing Nlm from the coupling coefficients
    std::mt19937 gen(5551212);
    std::uniform_int_distribution rand(-32768,32768);
    
    int lmax2 = lmax1;
    int lmax3 = lmax1+lmax2;
    auto lmx = [=](int l1, int m1, int l2, int m2) { // offset into the index vector
        return (l1*(l1+1)+m1)*(lmax1+1)*(lmax1+1) + l2*(l2+1)+m2;
    };
    NlmContainer<T> Nnorm = Nlm_normalization<T>(lmax3);

    T maxerr = 0.0, maxrelerr = 0, maxerrnorm = 0;
    for (int sample=0; sample<1000; sample++) {
        T x = rand(gen);
        T y = rand(gen);
        T z = rand(gen);
        T r = std::sqrt(x*x + y*y + z*z);
        x /= r; // Put test points on the unit sphere
        y /= r;
        z /= r;
        r = 1;

        NlmContainer<T> N = Nlm(lmax3, x, y, z);

        T thresh = 500*std::numeric_limits<T>::epsilon();
        
        for (int l1=0; l1<=lmax1; l1++) {
            for (int m1=-l1; m1<=l1; m1++) {
                for (int l2=0; l2<=lmax2; l2++) {
                    for (int m2=-l2; m2<=l2; m2++) {
                        auto [off,n] = index[lmx(l1,m1,l2,m2)];
                        T Ntest = 0.0;
                        for (int i=0; i<n; i++) {
                            auto [l3,m3,coeff] = j3[off+i];
                            Ntest += coeff*N.value(l3,m3);
                        }
                        T prod = N.value(l1,m1)*N.value(l2,m2);
                        T err =  prod - Ntest;
                        if (std::abs(err) > thresh) { //
                            std::cout << std::setw(3) << l1 << " " << std::setw(3) << m1 << " " << std::setw(3) << l2 << " " << std::setw(3) << m2 << "  " << std::setw(15) << err << std::endl;
                        }
                        maxerr = std::max(err,maxerr);
                        maxerrnorm = std::max(err*Nnorm.value(l1,m1)*Nnorm.value(l2,m2),maxerrnorm);
                        T prodnormed = prod*Nnorm.value(l1,m1)*Nnorm.value(l2,m2);
                        if (std::abs(prodnormed) > thresh)
                            maxrelerr = std::max(err/prod,maxrelerr);
                    }
                }
            }
        }
    }
    std::cout << "    maxerr in J3 reconstruction on unit sphere " << maxerr << std::endl;
    std::cout << " maxrelerr in J3 reconstruction on unit sphere " << maxrelerr << std::endl;
    std::cout << "maxnormerr in J3 reconstruction on unit sphere " << maxerrnorm << std::endl;

    /*
     lmax = 8
     Using GL
        Testing j3 with double
        j3 size 26113 26244
            maxerr in J3 reconstruction on unit sphere 3.053113e-16
         maxrelerr in J3 reconstruction on unit sphere 4.518745e-03
        maxnormerr in J3 reconstruction on unit sphere 1.123624e-14

        Testing j3 with double-double
        j3 size 26113 26244
            maxerr in J3 reconstruction on unit sphere 1.386670e-31
         maxrelerr in J3 reconstruction on unit sphere 2.778597e-12
        maxnormerr in J3 reconstruction on unit sphere 6.280683e-30

        Testing j3 with double-double coeffs but in double
            maxerr in J3 reconstruction on unit sphere 2.220446e-16
         maxrelerr in J3 reconstruction on unit sphere 4.738321e-04
        maxnormerr in J3 reconstruction on unit sphere 3.322992e-15

    Using Lebedev (which we only have in double)
        Testing j3 with double
        j3 size 26113 26244
            maxerr in J3 reconstruction on unit sphere 3.663736e-15
         maxrelerr in J3 reconstruction on unit sphere 3.880582e-02
        maxnormerr in J3 reconstruction on unit sphere 1.249279e-13

     lmax = 16
     Using GL
        Testing j3 with double
        j3 size 603487 668168
            maxerr in J3 reconstruction on unit sphere 3.053113e-16
         maxrelerr in J3 reconstruction on unit sphere 2.373996e+04
        maxnormerr in J3 reconstruction on unit sphere 1.380099e-08

        Testing j3 with double-double
        j3 size 603487 668168
            maxerr in J3 reconstruction on unit sphere 1.386670e-31
         maxrelerr in J3 reconstruction on unit sphere 1.417492e+20
        maxnormerr in J3 reconstruction on unit sphere 1.380100e-08

        Testing j3 with double-double coeffs but in double
            maxerr in J3 reconstruction on unit sphere 2.220446e-16
         maxrelerr in J3 reconstruction on unit sphere 2.373996e+04
        maxnormerr in J3 reconstruction on unit sphere 1.380100e-08

     Using Lebedev (which we only have in double)
            maxerr in J3 reconstruction on unit sphere 3.663736e-15
         maxrelerr in J3 reconstruction on unit sphere 2.373997e+04
        maxnormerr in J3 reconstruction on unit sphere 1.380101e-08
     */
}



int main() {
    constants<double>::init();
    constants<dd_real>::init();
    constants<qd_real>::init();
    test_Nlm<double>();
    test_Nlm<dd_real>();
    test_Nlm<qd_real>();

    {
        std::cout << "Testing j3 with double\n";
        auto [lmax1, index, j3] = threeJ<double>();
        test_j3(lmax1, index, j3);
    }

    {
        std::cout << "Testing j3 with double-double\n";
        auto [lmax1, index, j3] = threeJ<dd_real>();
        test_j3(lmax1, index, j3);

        std::cout << "Testing j3 with double-double coeffs but in double\n";
        std::vector<std::tuple<int,int,double>> j3d;
        j3d.reserve(j3.size());
        for (auto [l,m,v] : j3) {j3d.push_back({l,m,to_double(v)});}
        test_j3(lmax1, index, j3d);
    }        

    //std::ifstream f("gauleg.220bit");
 
    return 0;
}
