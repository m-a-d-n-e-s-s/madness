#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <limits>
#include <cmath>
#include <myqd.h>

template <typename T>
T myexpm1(const T& x) {
    if (std::abs(x)<T(0.5)) {
        T s = x;
        T t = x;
        // 9 for single
        // 15 for double
        // 25 for dd
        // 42 for qd
        for (size_t i=2; i<100; i++) {
            t *= x/i;
            s += t;
            if (std::abs(t) < std::numeric_limits<T>::epsilon()) break;
        }
        return s;
    }
    return exp(x)-T(1);
}

namespace std {
    dd_real expm1(const dd_real& x) {
        return ::myexpm1<dd_real>(x);
    }

    qd_real expm1(const qd_real& x) {
        return ::myexpm1<qd_real>(x);
    }
}

template <typename T>
struct factorial_cache {
    static constexpr size_t N = 171;   // more than 171 will overflow exponent
    static std::array<T,N+1> f;      // +1 for 0
    factorial_cache() {
        f[0] = 1;
        for (int i=1; i<=int(N); i++) f[i] = f[i-1]*i; // int since T might be dd_real/qd_real
    }
};

template <typename T> 
struct double_factorial_cache {
    static constexpr size_t N = 300; // more than 300 will overflow exponent
    static std::array<T,N+1> f; // +1 for 0
    double_factorial_cache() {
        f[0] = 1;
        f[1] = 1;
        for (int i=2; i<=int(N); i++) f[i] = f[i-2]*i; // int since T might be dd_real/qd_real
    }
};

template <typename T> std::array<T,172> factorial_cache<T>::f = {};
template <typename T> std::array<T,301> double_factorial_cache<T>::f = {};

// If you might overflow size_t use floating point number for T
template <typename T>
T factorial(size_t i) {
    return factorial_cache<T>::f[i];
}

// If you might overflow size_t use floating point number for T
template <typename T>
T double_factorial(size_t i) {
    return double_factorial_cache<T>::f[i];
}

// Just for debugging precision issues
template <typename T>
T mypow(const T& x, const int k) {
    T r = 1;
    T term = x;
    int kk = k;
    while (true) {
        if (kk & 0x1) r *= term;
        kk >>= 1;
        if (kk == 0) break;
        term *= term;
    }
    return r;
}

template <typename T>
T JSphericalSeries(const size_t l, const T& r) {
    auto c = [](size_t l, size_t i) {
        //return std::ldexp( (factorial<T>((i+l)/2)/factorial<T>((i-l)/2))/factorial<T>(i+l+1), int(l)); // correct
        size_t k = (i-l)>>1;
        return std::ldexp(T(1)/(factorial<T>(k)*double_factorial<T>(i+l+1)), -int(k)); // also correct and maybe less prone to underflow
    };

    // dp r=0.3 i<(l+12)
    // dp r=0.5 similar 14
    // dp r=1  i<(l+24)

    // dd r=0.5 i<(l+24)
    // dd r=1  i<(l+30)

    // qd r=0.5<(l+42)
    // qd r=1 i<(l+50)

    // for r <= 1
    // size_t maxi=7;
    // if constexpr (std::is_same<T,double>::value) maxi=15;
    // else if constexpr (std::is_same<T,dd_real>::value) maxi=27;
    // else if constexpr (std::is_same<T,qd_real>::value) maxi=47;

    // for r <= 0.5
    size_t maxi=7;
    if constexpr (std::is_same<T,double>::value) maxi=13;
    else if constexpr (std::is_same<T,dd_real>::value) maxi=23;
    else if constexpr (std::is_same<T,qd_real>::value) maxi=40;

    T s = 0;
    T ri = 1;
    T r2 = r*r;
    for (size_t i=l; i<(l+maxi); i+=2) {
        //T term = c(l,i)*std::pow(r,int(i-l));
        T term = c(l,i)*ri;
        s += term;
        //std::cout << i << " " << to_str(c(l,i)) << " " << ri << " " << to_str(term) << " " << to_str(s) << std::endl;
        ri *= r2;
        if (ri < std::numeric_limits<T>::epsilon()) break;
    }
    s *= std::pow(r,int(l));
    
    //return s*std::pow(r,int(l)); // Avoid using pow(float,float) --- pow(float,int) more accurate.
    return s;
}


// Computes exp(-r)*j(l,r) using Miller's algorithm and downward recursion.
// In Maple:
// j := (l, x) -> Re(I^(-l)*sqrt(-1/2*I*Pi/x)*BesselJ(l + 1/2, x*I))
template <typename T>
T Jsphericalscaled(size_t l, T r)
{
    if (r <= T(0.5)) {
        T ee = std::expm1(-r)+T(1);
        T jj = JSphericalSeries(l,r);
        //std::cout << "ee=" << to_str(ee) << " jj=" << to_str(jj) << " " << to_str(ee*jj) << std::endl;
        return (std::expm1(-r)+T(1))*JSphericalSeries(l,r);
    }
    
    T jj0 = -std::expm1(-T(2)*r)/(T(2)*r);
    if (l == 0) return jj0;

    // else if (l == 1) {
    //     if (r < T(0.02)) {
    //         return (T(1)/T(3) + (-T(1)/T(3) + (T(1)/T(5) + (-T(4)/T(45) + (T(2)/T(63) + (-T(1)/T(105) + r/T(405))*r)*r)*r)*r)*r)*r;
    //     }
    //     else {
    //         T s = std::expm1(-2*r);
    //         return (r*(2+s) + s)/(2*r*r);
    //     }
    // }

    // T j0=std::numeric_limits<T>::min(), j1=0;
    // size_t L = std::ceil(std::max(r + T(20), (T) (l + 10)));
    // while (l < L) {
    //     T jm1 = j0*T(2*L+1)/r + j1;
    //     j1 = j0;
    //     j0 = jm1;
    //     L -= 1;
    // }

    // Errors now nearly all at small r --- need Taylor series there
    // Start recurring down from L with j(L+1)=0
    size_t L = std::ceil(std::max(to_double(r) + 20, double(l + 12))); // for double
    if constexpr (std::is_same_v<T,dd_real>) L = std::ceil(std::max(to_double(r) + 50, double(l + 14)));
    else if constexpr (std::is_same_v<T,qd_real>) L = std::ceil(std::max(to_double(r) + 70, double(l + 18)));

    // j0 implicitly = 1 in this loop since we are computing the ratio j1 = j(l+1)/j(l)
    T j1=0;
    T rr = T(1)/r;
    while (l < L) {
        j1 = T(1)/(T(int(2*L+1))*rr + j1);
        //j1 = T(1)/(T(int(2*L+1))/r + j1);
        //j1 = r/(T(int(2*L+1)) + j1*r);
        L -= 1;
    }
    
    T jl=1, j0=1;
    while (L>0) {
        T jm1 = j0*T(int(2*L+1))/r + j1;
        j1 = j0;
        j0 = jm1;
        L -= 1;
    }
    return(jl*jj0/j0);
}

template <typename T>
class Matrix {
    size_t n, m;
    std::vector<T> data;

public:
    Matrix() : n(0), m(0) {}
    Matrix(size_t n, size_t m) : n(n), m(m), data(n*m) {}
    T& operator()(size_t i, size_t j) { return data[i*m+j]; }
    const T& operator()(size_t i, size_t j) const { return data[i*m+j]; }
    size_t rows() const { return n; }
    size_t cols() const { return m; }
};

template <typename T> struct bessel_traits {};
template <> struct bessel_traits<float>  {static constexpr size_t maxL=4;};
template <> struct bessel_traits<double> {static constexpr size_t maxL=32;};
template <> struct bessel_traits<dd_real> {static constexpr size_t maxL=32;}; // was 50
template <> struct bessel_traits<qd_real> {static constexpr size_t maxL=32;}; // was 50

template <typename T>
std::tuple<std::vector<T>, Matrix<T>, Matrix<T>>
load_bessel_test_data() {
    size_t maxL, nR;
    std::ifstream f("bessel.txt");
    if (!f) throw std::runtime_error("cannot open bessel.txt");
    f >> maxL; 
    f >> nR; 
    maxL = std::min(maxL, bessel_traits<T>::maxL);

    std::string line;
    std::getline(f, line); std::getline(f, line);

    Matrix<T> h(maxL+1, nR), j(maxL+1, nR);
    std::vector<T> r(nR);

    std::string rs, js, hs;
    
    // for each l and for each r read l, r, j, h
    for (size_t l=0; l<=maxL; l++) {
        for (size_t i=0; i<nR; i++) {
            size_t ll;
            f >> ll >> rs >> js >> hs;
            r[i] = from_str<T>(rs);
            j(l,i) = from_str<T>(js);
            h(l,i) = from_str<T>(hs);
            if (l != ll) throw std::runtime_error("l mismatch");
        }
    }
    return std::make_tuple(r, j, h);
}

template <typename T>
void test_bessel() {
    std::vector<T> r;
    Matrix<T> j, h;
    std::tie(r, j, h) = load_bessel_test_data<T>();
    const size_t nR=r.size(), maxL=j.rows()-1; 
    
    for (size_t l=0; l<=maxL; l+=1) {
        for (size_t i=0; i<nR; i+=1) {
            //for (size_t l=40; l<=40; l+=5) {
            //for (size_t i=300; i<=300; i+=20) {
            //std::cout << "r " << i << " " << r[i] << std::endl;
            T j0 = Jsphericalscaled(l, r[i]);
            T err = std::abs(j0-j(l,i));
            T relerr = err/j(l,i);
            // double gives 20eps
            // qd gives 10eps
            // dd gives 50eps nearly everywhere
            T fudge = 20;
            if constexpr (std::is_same_v<T, double>) fudge = 20;
            else if constexpr (std::is_same_v<T, qd_real>) fudge = 10;
            else if constexpr (std::is_same_v<T, dd_real>) fudge = 50;

            if (err > 2e-323) { // really tiny values will be denormed
                if (relerr > fudge*std::numeric_limits<T>::epsilon()) {
                    std::cout << "l=" << l << " r=" << r[i] << " j=" << to_str(j(l,i)) << " j0=" << to_str(j0) << " err=" << err << " relerr=" << relerr << " " << std::endl;
                }
            }
        }
    }
}

template <typename T>
void test() {
    std::string s = "2.1565466832389396343019344302735465581977651653240437376157577064561689e-191"; // smaller than -191 loses precision in the conversion
    //std::string s = "6.5198687841476257935730899787031683422353037643557853848016038462712578e-237";
    T x = from_str<T>(s);
    std::string t = x.to_string(); //to_str(x);
    std::cout << s << std::endl;
    std::cout << t << std::endl;

    char buf[100];
    int expn;
    x.to_digits(buf,expn);
    std::cout << buf << " " << expn << std::endl;
}

int main() {
    unsigned int oldcw; // For x86
    fpu_fix_start(&oldcw);

    {
        // static initialize of dd/qd does not work ... need to do it
        // in main code after FPU fix, and thread-safe initialization
        // of dd/qd package
        factorial_cache<double> f;
        factorial_cache<dd_real> fd;
        factorial_cache<qd_real> fq;
        double_factorial_cache<double> df;
        double_factorial_cache<dd_real> dfd;
        double_factorial_cache<qd_real> dfq;
    }

    std::cout << "7! " << factorial<double>(7) << std::endl;
    std::cout << "7! " << factorial<dd_real>(7) << std::endl;
    std::cout << "7! " << factorial<qd_real>(7) << std::endl;

    std::cout << "7!! " << double_factorial<double>(7) << std::endl;
    std::cout << "7!! " << double_factorial<dd_real>(7) << std::endl;
    std::cout << "7!! " << double_factorial<qd_real>(7) << std::endl;

    test<dd_real>();
    test<qd_real>();
    
    // double r = 0.5;
    // for (size_t l=0; l<10; l++)
    //     std::cout << l << " " << Jsphericalscaled(l, r) << std::endl;

    std::cout << "double " << std::numeric_limits<double>::epsilon() << std::endl;
    test_bessel<double>();
    std::cout << "dd_real " << std::numeric_limits<dd_real>::epsilon() << std::endl;
    test_bessel<dd_real>();
    std::cout << "qd_real " << std::numeric_limits<qd_real>::epsilon() << std::endl;
    test_bessel<qd_real>();

    // std::cout << JSphericalSeries(2,0.3) << std::endl;
    // std::cout << Jsphericalscaled(2,0.3)*std::exp(0.3) << std::endl;
    // test_bessel<double>();

    //JSphericalSeries<qd_real>();
    return 0;
}
