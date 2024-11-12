//#pragma GCC optimize "-fno-associative-math"

#include <typeinfo>
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <tuple>
#include <array>
#include <limits>
#include <cmath>
#include <myqd.h>

#include <algorithm>
#include <numeric>

template <typename T>
class Matrix {
    size_t n, m;
    std::vector<T> data;

public:
    Matrix() : n(0), m(0) {}
    Matrix(size_t n, size_t m) : n(n), m(m), data(n*m) {}
    T& operator()(size_t i, size_t j) { return data.at(i*m+j); }
    const T& operator()(size_t i, size_t j) const { return data.at(i*m+j); }
    size_t rows() const { return n; }
    size_t cols() const { return m; }
};


size_t ggi = 9999999;
std::vector<dd_real> ggr;
Matrix<dd_real> ggj, ggh;

// print a vector
template <typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<T>& v) {
    s << "[ ";
    for (const auto& x : v) s << x << " ";
    s << "]";
    return s;
}
// print a pair
template <typename T, typename U>
std::ostream& operator<<(std::ostream& s, const std::pair<T,U>& p) {
    s << "(" << p.first << ", " << p.second << ")";
    return s;
}


template<typename T>
T log10_estimate(T x) {
    int e;
    T m = std::frexp(x, &e);
    return (e-1)*T(0.30103); // log10(2) = 0.30103
}

double log10_estimate(double x) {
    int e =  (int((*((long*)(&x))) >> 52) & 0x7FF) -1023;    
    return e*0.30103; // log10(2) = 0.30103
}

float log10_estimate(float x) {
    std::cout << "float\n";
    int e =  (int((*((int*)(&x))) >> 23) & 0xFF) - 127;    
    return e*0.30103; // log10(2) = 0.30103
}

template <typename T>
T myexpm1(const T& x) {
    if (std::abs(x)<T(0.5)) {
        T s = x;
        T t = x;
        // 9 for single
        // 15 for double
        // 25 for dd
        // 42 for qd
        for (size_t i=2; i<42; i++) {
            t *= x/i;
            s += t;
            //if (std::abs(t) < std::numeric_limits<T>::epsilon()) break;
        }
        return s;
    }
    return std::exp(x)-T(1);
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
    static inline std::array<T,N+1> f = {};      // +1 for 0
    factorial_cache() {
        f[0] = 1;
        for (int i=1; i<=int(N); i++) f[i] = f[i-1]*i; // int since T might be dd_real/qd_real
    }
};

template <typename T> 
struct double_factorial_cache {
    static constexpr size_t N = 300; // more than 300 will overflow exponent
    static inline std::array<T,N+1> f = {}; // +1 for 0
    double_factorial_cache() {
        f[0] = 1;
        f[1] = 1;
        for (int i=2; i<=int(N); i++) f[i] = f[i-2]*i; // int since T might be dd_real/qd_real
    }
};

// not needed since c++ 17 inline variables
//template <typename T> std::array<T,172> factorial_cache<T>::f = {};
//template <typename T> std::array<T,301> double_factorial_cache<T>::f = {};

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

// Barycentric interpolation class using N uniformly spaced-intervals with m uniformly spaced points within each interval (including endpoints)
template <typename T>
class BarycentricX {
    size_t N; // Number of larger intervals
    size_t m; // Number of points in each interval
    T a; // Left endpoint of the interval
    T b; // Right endpoint of the interval
    T H; // Larger interval size
    T h; // Smaller interval size

    std::vector<T> y; // Vector of function values
    std::vector<T> w; // Vector of m Barycentric weights

    static size_t factorial(size_t n) {
        size_t f = 1;
        for (size_t i=1; i<=n; i++) f *= i;
        return f;
    }

    static std::vector<T> weights(const size_t m) {
        std::vector<T> w(m);
        size_t n = m-1;
        size_t nfac = factorial(n);
        T sign = 1;
        for (size_t j=0; j<=n; j++) {
            w[j] = sign*nfac/(factorial(j)*factorial(n-j));
            sign = -sign;
        }
        return w;
    }

    public:

    template <typename funcT>
    BarycentricX(size_t N, size_t m, T a, T b, funcT f) 
        : N(N)
        , m(m)
        , a(a)
        , b(b)
        , H((b-a)/N)
        , h(H/(m-1))
        , y(N*(m-1)+1)
        , w(weights(m))
    {
        assert(m>=2);
        assert(N>=1);
        for (size_t i=0; i<y.size(); i++) {
            T x = a + i*h;
            y[i] = f(x);
        }
        // std::cout << "H = " << H << std::endl;
        // std::cout << "h = " << h << std::endl;
        // std::cout << "y = " << y << std::endl;
        // std::cout << f(a) << " " << f(b) << std::endl;
        // std::cout << "w = " << w << std::endl;
    }

    T operator()(T x) {
        //std::cout << "input x = " << x << " " << N << " " << H << std::endl;
        assert(x>=a && x<b);
        x -= a;
        const size_t i = x/H;
        x -= i*H;
        //std::cout << "i = " << i << " " << x << std::endl;

        T num = 0;
        T den = 0;
        for (size_t j=0; j<m; j++) {
            T diff = x-j*h;
            if (diff == 0) return y[i*(m-1) + j];
            //std::cout << "diff = " << diff << " j " << j << std::endl;
            T wj = w[j]/diff;
            num += wj*y[i*(m-1)+j];
            den += wj;
        }
        return num/den;
    }

    std::vector<T> operator()(const std::vector<T>& x) {
        std::vector<T> r(x.size());
        for (size_t i=0; i<x.size(); i++) r[i] = operator()(x[i]);
    }    

    template <typename funcT>
    std::pair<T,T> maxerr(funcT f, size_t oversample=5) {
        T testh = h/oversample;
        T err = 0;
        T relerr = 0;
        T x = a;
        while (x < b) {
            T fx = f(x);
            T abserr = std::abs(fx - (*this)(x));
            err = std::max(err, abserr);
            if (fx != 0) relerr = std::max(relerr, abserr/std::abs(fx));
            x += testh;
        }
        return std::make_pair(err, relerr);
    }    
};

// Barycentric interpolation class using N uniformly spaced-intervals (so N+1 points) and order m interpolation
template <typename T>
class Barycentric {
    size_t N; // Number of larger intervals
    size_t m; // Number of points for interpolation
    T a; // Left endpoint of the interval
    T b; // Right endpoint of the interval
    T h; // Interval size

    std::vector<T> y; // Vector of function values
    std::vector<T> w; // Vector of m Barycentric weights

    static size_t factorial(size_t n) {
        size_t f = 1;
        for (size_t i=1; i<=n; i++) f *= i;
        return f;
    }

    static std::vector<T> weights(const size_t m) {
        std::vector<T> w(m);
        size_t n = m-1;
        size_t nfac = factorial(n);
        T sign = 1;
        for (size_t j=0; j<=n; j++) {
            w[j] = sign*nfac/(factorial(j)*factorial(n-j));
            sign = -sign;
        }
        return w;
    }

    public:

    Barycentric() {}

    template <typename funcT>
    Barycentric(size_t N, size_t m, T a, T b, funcT f) 
        : N(N)
        , m(m)
        , a(a)
        , b(b)
        , h((b-a)/N)
        , y(N+1)
        , w(weights(m))
    {
        assert(m>=2);
        assert(N>=m);
        for (size_t i=0; i<y.size(); i++) {
            T x = a + i*h;
            y[i] = f(x);
        }
        // std::cout << "h = " << h << std::endl;
        // std::cout << "y = " << y << std::endl;
        // std::cout << f(a) << " " << f(b) << std::endl;
        // std::cout << "w = " << w << std::endl;
    }

    T operator()(T x) {
        //std::cout << "input x = " << x << " " << N << " " << h << std::endl;
        assert(x>=a && x<b);
        x -= a;
        size_t i = std::min(size_t(x/h), N-(m+1)/2);
        //std::cout << "initial i = " << i << " " << x << std::endl;
        i -= std::min((m-1)/2, i);
        x -= i*h;
        //std::cout << "shifted i = " << i << " " << x << std::endl;

        T num = 0;
        T den = 0;
        for (size_t j=0; j<m; j++) {
            T diff = x-j*h;
            if (diff == 0) return y[i + j];
            T wj = w[j]/diff;
            //std::cout << "diff = " << diff << " j " << j << std::endl;
            num += wj*y[i+j];
            den += wj;
        }
        return num/den;
    }

    template <typename funcT>
    std::pair<T,T> maxerr(funcT f, size_t oversample=5) {
        T testh = h/oversample;
        T err = 0;
        T relerr = 0;
        T x = a;
        while (x < b) {
            T fx = f(x);
            T abserr = std::abs(fx - (*this)(x));
            err = std::max(err, abserr);
            if (fx != 0) relerr = std::max(relerr, abserr/std::abs(fx));
            if (relerr > 1e-14) {
                std::cout << "bad err " << x << " " << fx << " " << (*this)(x) << " " << abserr << " " << relerr << std::endl;
                break;
            }
            x += testh;
        }
        return std::make_pair(err, relerr);
    }    
};

Barycentric<double> b33;

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


// template <typename T>
// T JSphericalRatioSeries(const size_t l, const T& r) {
//     T jl = 0;  // j_l
//     T jl1 = 0; // j_{l-1}

//     T rk = 1;
//     T r2 = r*r;
//     for (size_t i=l; i<(l+maxi); i+=2) {

//         s += term;
//         //std::cout << i << " " << to_str(c(l,i)) << " " << ri << " " << to_str(term) << " " << to_str(s) << std::endl;
//         ri *= r2;
//         if (ri < std::numeric_limits<T>::epsilon()) break;
//     }

//     return (r*jl)/((2*l+1)*jl1);
// }

template <typename T>
T JsphericalscaledX(size_t l, T r);

template <typename T>
T Jsphericalscaled(size_t l, T r) {
    if (r == T(0.0)) return T(0.0);
    const T one = T(1.0);
    const T half = T(0.5);
    T rr = one/r;
    if (r < T(std::max(10.0,5.0*l))) {
        //size_t n = l+40; // dp ... these with r < T(std::max(10.0,5.0*l))
        //size_t n = l+70; // dd
        //size_t n = l+120; // qd

        size_t n;
        if constexpr (std::is_same<T,double>::value) n = l+40;
        else if constexpr (std::is_same<T,dd_real>::value) n = l+70;
        else if constexpr (std::is_same<T,qd_real>::value) n = l+120;
        T R;
        T j0;

        // // Make initial guess for downward recursion ratio R_n
        // if (r > std::max(2*l,size_t(10))) { // Recur up and use downward recursion to correct accumulated rounding errors
        //     T expm2r = std::exp(-T(2)*r);
        //     j0 = (one-expm2r)*half*rr;
        //     if (l == 0) return j0;
            
        //     R = (r-one+(r+one)*expm2r)/(r*(one-expm2r));
        //     size_t m = 1;
        //     while (m < n) {
        //         std::cout << "      m " << m << " " << to_str(R) << std::endl;
        //         R = T(1)/R - T(2*m+1.0)*rr;
        //         m += 1;
        //     }
        //     T Rtest = JsphericalscaledX(n,r)/JsphericalscaledX(n-1,r);
        //     std::cout << "Guess " << l << " " << to_str(r) << " " << to_str(R) << " " << to_str(Rtest) << " " << to_str((R-Rtest)/Rtest) << std::endl;
        // }
        // else {
        //     j0 = -myexpm1(-T(2)*r)*half*rr; // for small r use expm1 to avoid cancellation
        //     if (l == 0) return j0;
        //     T t = T(2*n+1.0)*rr;
        //     R = (std::sqrt(t*t + T(4)) - t)*half;
        // }

        j0 = -myexpm1(-T(2)*r)*half*rr; // for small r use expm1 to avoid cancellation
        if (l == 0) return j0;

        R = (r*(2*r*r + (2*r + 1 + n)*n)) / (2*r*r*r + (1 + (2 + 4*r)*r + (4*r + 3 + 2*n)*n)*n);
        // {
        //     T Rtest = JsphericalscaledX(n,r)/JsphericalscaledX(n-1,r);
        //     std::cout << "Guess " << l << " " << to_str(r) << " " << to_str(R) << " " << to_str(Rtest) << " " << to_str((R-Rtest)/Rtest) << std::endl;
        // }            
        
        // Recur downward
        while (--n>=l) {
            R = T(1.0)/(R + T(2*n+1.0)*rr); // Compute R_n from R_{n+1}
        }

        T s = R;
        for (int n=l-1; n>=1; n--) {
            R = T(1.0)/(R + T(2*n+1.0)*rr); // Compute R_n from R_{n+1}
            s *= R;
            //std::cout << "R " << n << " " << to_str(R) << " " << to_str(s) << std::endl;
        }
        //std::cout << "final s " << to_str(s) << std::endl;

        //std::cout << to_str(j0) << std::endl;
        //T j0test = from_str<T>("0.4637228398768963676128862664265477510043731236013416774891131525604899694555197938928784651025466370");
        //std::cout << "j0 err " << to_str((j0test - j0)/j0) << std::endl;
        return s*j0;
    }
    else {
        const T expm2r = std::exp(-T(2)*r);
        T j0 = (one-expm2r)*half*rr;
        if (l == 0) return j0;
        T j1 = ((one+expm2r)*half - j0)*rr;
        if (l == 1) return j1;
        for (size_t n=1; n<l; n++) {
            T j2 = j0 - T(2*n+1.0)*j1*rr;
            j0 = j1;
            j1 = j2;
        }
        return j1;
    }
}

template <typename T>
T Rsmallr(size_t l, T r) {
    size_t tl1 = 2*l+1; 
    size_t tl3 = 2*l+3;
    size_t tl3sq = tl3*tl3;
    size_t tl5 = 2*l+5;
    size_t tl7 = 2*l+7;
    size_t tl9 = 2*l+9;
    size_t tl1sq = tl1*tl1;
    size_t tl1cu = tl1sq*tl1;
    size_t tl1fo = tl1sq*tl1sq;
    size_t tl1fi = tl1sq*tl1cu;
    T rsq = r*r;

    //     (1/(2*l + 1) + (-1/((2*l + 1)^2*(2*l + 3)) + (2/((2*l + 5)*(2*l + 3)*(2*l + 1)^3) + ((-10*l - 17)/((2*l + 1)^4*(2*l + 7)*(2*l + 3)^2*(2*l + 5)) + (28*l + 62)*r^2/((2*l + 1)^5*(2*l + 7)*(2*l + 9)*(2*l + 3)^2*(2*l + 5)))*r^2)*r^2)*r^2)*r

    return (T(1)/tl1 + (-T(1)/(tl1sq*tl3) + (T(2)/(tl5*tl3*tl1cu) + (T(-10*l - 17)/(tl1fo*tl7*tl3sq*tl5) + T(28*l + 62)*rsq/(tl1fi*tl7*tl9*tl3sq*tl5))*rsq)*rsq)*rsq)*r;
}

template <typename T>
T RsmallrX(size_t l, T r) {
    T rr = T(1)/r;
    T r1 = r/T(2*l+1.0); T r1sq = r1*r1; T r1cub = r1sq*r1; T r1four = r1sq*r1sq; T r1five = r1four*r1;
    T r3 = r/T(2*l+3.0); T r3sq = r3*r3;
    T r5 = r/T(2*l+5.0);
    T r7 = r/T(2*l+7.0);
    T r9 = r/T(2*l+9.0);
    T r11 = r/T(2*l+11.0);
    //(-10*l - 17)*r^7/((2*l + 1)^4*(2*l + 7)*(2*l + 3)^2*(2*l + 5)) + (28*l + 62)*r^9/((2*l + 1)^5*(2*l + 7)*(2*l + 9)*(2*l + 3)^2*(2*l + 5))
    return ((r1 - r1sq*r3) + 2*r5*r3*r1cub) - T(10*l+17.0)*r1four*r7*r3sq*r5*rr + T(28*l+62.0)*r1five*r7*r9*r3sq*r5*rr;
    //- T(336*l*l*l + 2392*l*l + 5564*l + 4146.0)*r11*r5*r5*r3sq*r3*r7*r1four*r1sq*r9*rr*rr*rr;
}

template <typename T>
std::vector<T> JsphericalscaledVecGood(size_t maxl, T r) {
    assert(maxl == 32);
    if (r == T(0.0)) return std::vector<T>(maxl+1,T(0.0));
    const T eps = std::numeric_limits<T>::epsilon();
    const T min = std::numeric_limits<T>::min();

    const T one = T(1.0);
    const T two = T(2.0);
    const T half = T(0.5);
    
    T rr = one/r;
    T R;

    bool doprint = (ggi != 9999999) && (std::abs(r - 0.00191116) < 0.00001);
    //bool doprint = false;
    if (doprint) {
        std::cout << "printing " << r << " " << (ggj(33,ggi)>min) << " " << min << std::endl;
    }

    std::vector<T> j(maxl+1);

    T rcut = 20.0; // So that exp(-2*rcut) is much less than epsilon including gradual underflow
    if constexpr      (std::is_same<T,float>::value) rcut = 10.0;
    else if constexpr (std::is_same<T,double>::value) rcut = 20.0;
    else if constexpr (std::is_same<T,dd_real>::value) rcut = 40.0;
    else if constexpr (std::is_same<T,qd_real>::value) rcut = 80.0;

    auto RELERR = [&](const T& a, const T&b) {return (std::abs(a-b)/std::abs(b))/eps;};

    if (r < std::max(rcut,T(12.5*maxl))) {
        if (r < 0.5) { // was 2.5
            R = Rsmallr(48,r);
        }
        else {
            R = b33(r); // Very accurate guess for R_38, slight error at small r but going from 38 to 33 will fix it // was 38 now 48
        }
        size_t n = 48;
        while (--n>maxl) {
            R = T(1.0)/(R + T(2*n+1.0)*rr); // Compute R_n from R_{n+1}
        }

        if (doprint) {
            T Rtest = ggj(33,ggi)/ggj(32,ggi);
            T relerr = RELERR(R,Rtest);
            if (relerr > 0) {
                std::cout << 33 << " Rtest " << Rtest << " " << R << " " << ggi << " " << r << " " << ggr[ggi] << " " << ggj(33,ggi) << " " << ggj(32,ggi) << " " << relerr << std::endl;
            }
        }

        T Rsave = R;
        
        // for (int l=maxl; l>0; l--) {
        //     j[l] = R = T(1.0)/(R + T(2*l+1.0)*rr); // Compute R_n from R_{n+1}
        // }
        // std::cout << "OLD " << r << " " <<  maxl << std::endl;
        // std::cout << "OLD " <<  j << std::endl;
        // j[0] = -myexpm1(-two*r)*half*rr; // for small r use expm1 to avoid cancellation
        // for (size_t l=1; l<=maxl; l++) {
        //     j[l] *= j[l-1];
        // }
        // std::cout << "OLD " <<  j << std::endl;

        // j = std::vector<T>(maxl+1,T(-1.0));

        T j1 = min * std::max(one, Rsave*r/T(2*maxl+1.0)); //was Rsave ... now trying to avoid underflow in s below
        T j0 = j1/Rsave; // was 1
        j[maxl] = j0;

        if (doprint) {
            T R = Rsave;
            T Rtest = ggj(maxl+1,ggi)/ggj(maxl,ggi);
            T relerr = RELERR(R,Rtest);
            if (relerr > 0) {
                std::cout << maxl+1 << " Rtest " << Rtest << " " << R << " " << ggi << " " << r << " " << ggr[ggi] << " " << ggj(maxl+1,ggi) << " " << ggj(maxl,ggi) << " " << relerr << std::endl;
            }
        }
        
        for (int l=maxl; l>0; l--) {
            j[l-1] = j1 + (T(2*l+1.0)*rr)*j0;
            
            //std::cout << "NEW " << l << " " <<  j1 << " " << j0 << " " << j[l-1] << std::endl;
            j1 = j0;
            j0 = j[l-1];

            if (doprint) {
                T R = j1/j0;
                T Rtest = ggj(l,ggi)/ggj(l-1,ggi);
                T relerr = RELERR(R,Rtest);
                if (relerr > 0) {
                    std::cout << l << " Rtest " << Rtest << " " << R << " " << ggi << " " << r << " " << ggr[ggi] << " " << ggj(l,ggi) << " " << ggj(l-1,ggi) << " " << relerr << std::endl;
                }
            }

            
        }
        //std::cout << "NEW " <<  j << std::endl;

        if (doprint) {
            T R = j[32]/j[0];
            T Rtest = ggj(32,ggi)/ggj(0,ggi);
            T relerr = RELERR(R,Rtest);
            if (relerr > 0) {
                std::cout << 99 << " Rtestall " << Rtest << " " << R << " " << ggi << " " << r << " " << ggr[ggi] << " " << ggj(0,ggi) << " " << ggj(32,ggi) << " " << relerr << std::endl;
            }
        }
        
        
        j0 = -myexpm1(-two*r)*half*rr; // for small r use expm1 to avoid cancellation

        // if (ggi != 9999999) {
        //     T j0test = ggj(0,ggi);
        //     T relerr = RELERR(j0,j0test);
        //     if (relerr > 0) {
        //         std::cout << 0 << " Rtest00 " << ggi << " " << r << " " << ggr[ggi] << " " << j0 << " " << j0test << " " << relerr <<  std::endl;
        //     }
        // }

        T s = j0/j[0];
        if (doprint) std::cout << "s " << s << std::endl;
        
        if (s < min) {
            std::cout << "warning " << s << " " << min << std::endl;
        }
        j[0] = j0;

        double reps = to_double(dd_real(r) - dd_real(1.0)/dd_real(rr));
        
        for (size_t l=1; l<=maxl; l++) {
            j[l] *= s;

            j[l] += (j[l-1]-(one+T(l+1)*rr)*j[l])*reps; // correction due to inexact computing of 1/r
            
            if (doprint) {
                T R = j[l]/j[l-1];
                T Rtest = ggj(l,ggi)/ggj(l-1,ggi);
                T relerr = RELERR(R,Rtest);
                if (relerr > 0) {
                    std::cout << l << " RtestZ " << Rtest << " " << R << " " << ggi << " " << r << " " << ggr[ggi] << " " << ggj(l,ggi) << " " << ggj(l-1,ggi) << " " << relerr << std::endl;
                }
            }
        }
        //std::cout << s << std::endl;
        //std::cout << "NEW " <<  j << std::endl;
        //std::exit(0);
    }
    else {
        // const T expm2r = std::exp(-two*r);
        // j[0] = (one-expm2r)*half*rr;
        // j[1] = ((one+expm2r)*half - j[0])*rr;
        j[0] = half*rr;
        j[1] = (half - j[0])*rr;
        for (size_t l=1; l<maxl; l++) {
            j[l+1] = j[l-1] - T(2*l+1.0)*j[l]*rr;
        }
    }
    
    return j;
}

template <typename T>
std::vector<T> JsphericalscaledVecXX(size_t maxl, T r) {
    assert(maxl == 32);
    if (r == T(0.0)) return std::vector<T>(maxl+1,T(0.0));
    const T eps = std::numeric_limits<T>::epsilon();
    const T min = std::numeric_limits<T>::min();

    const T one = T(1.0);
    const T two = T(2.0);
    const T half = T(0.5);
    
    const T rr = one/r;
    T R;

    //bool doprint = (ggi != 9999999) && (std::abs(r - 0.00191116) < 0.00001);
    const bool doprint = false;
    if (doprint) {
        std::cout << "printing " << r << " " << (ggj(33,ggi)>min) << " " << min << std::endl;
    }

    std::vector<T> j(maxl+1);

    T rcut = 20.0; // So that exp(-2*rcut) is much less than epsilon including gradual underflow
    if constexpr      (std::is_same<T,float>::value) rcut = 10.0;
    else if constexpr (std::is_same<T,double>::value) rcut = 20.0;
    else if constexpr (std::is_same<T,dd_real>::value) rcut = 40.0;
    else if constexpr (std::is_same<T,qd_real>::value) rcut = 80.0;

    auto RELERR = [&](const T& a, const T&b) {return (std::abs(a-b)/std::abs(b))/eps;};

    if (r < std::max(rcut,T(12.5*maxl))) {
        if (r < 0.5) { // was 2.5
            R = Rsmallr(48,r);
        }
        else {
            R = b33(r); // Very accurate guess for R_38, slight error at small r but going from 38 to 33 will fix it // was 38 now 48
        }
        size_t n = 48;
        while (--n>maxl) {
            R = T(1.0)/(R + T(2*n+1.0)*rr); // Compute R_n from R_{n+1}
        }

        if (doprint) {
            T Rtest = ggj(33,ggi)/ggj(32,ggi);
            T relerr = RELERR(R,Rtest);
            if (relerr > 0) {
                std::cout << 33 << " Rtest " << Rtest << " " << R << " " << ggi << " " << r << " " << ggr[ggi] << " " << ggj(33,ggi) << " " << ggj(32,ggi) << " " << relerr << std::endl;
            }
        }

        T Rsave = R;
        
        T j1 = min * std::max(one, Rsave*r/T(2*maxl+1.0)); //was Rsave ... now trying to avoid underflow in s below
        T j0 = j1/Rsave; // was 1
        j[maxl] = j0;

        if (doprint) {
            T R = Rsave;
            T Rtest = ggj(maxl+1,ggi)/ggj(maxl,ggi);
            T relerr = RELERR(R,Rtest);
            if (relerr > 0) {
                std::cout << maxl+1 << " Rtest " << Rtest << " " << R << " " << ggi << " " << r << " " << ggr[ggi] << " " << ggj(maxl+1,ggi) << " " << ggj(maxl,ggi) << " " << relerr << std::endl;
            }
        }
        
        for (int l=maxl; l>0; l--) {
            j[l-1] = j1 + (T(2*l+1.0)*rr)*j0;
            
            //std::cout << "NEW " << l << " " <<  j1 << " " << j0 << " " << j[l-1] << std::endl;
            j1 = j0;
            j0 = j[l-1];

            if (doprint) {
                T R = j1/j0;
                T Rtest = ggj(l,ggi)/ggj(l-1,ggi);
                T relerr = RELERR(R,Rtest);
                if (relerr > 0) {
                    std::cout << l << " Rtest " << Rtest << " " << R << " " << ggi << " " << r << " " << ggr[ggi] << " " << ggj(l,ggi) << " " << ggj(l-1,ggi) << " " << relerr << std::endl;
                }
            }

            
        }
        //std::cout << "NEW " <<  j << std::endl;

        if (doprint) {
            T R = j[32]/j[0];
            T Rtest = ggj(32,ggi)/ggj(0,ggi);
            T relerr = RELERR(R,Rtest);
            if (relerr > 0) {
                std::cout << 99 << " Rtestall " << Rtest << " " << R << " " << ggi << " " << r << " " << ggr[ggi] << " " << ggj(0,ggi) << " " << ggj(32,ggi) << " " << relerr << std::endl;
            }
        }
        
        
        j0 = -myexpm1(-two*r)*half*rr; // for small r use expm1 to avoid cancellation

        // if (ggi != 9999999) {
        //     T j0test = ggj(0,ggi);
        //     T relerr = RELERR(j0,j0test);
        //     if (relerr > 0) {
        //         std::cout << 0 << " Rtest00 " << ggi << " " << r << " " << ggr[ggi] << " " << j0 << " " << j0test << " " << relerr <<  std::endl;
        //     }
        // }

        T s = j0/j[0];
        if (doprint) std::cout << "s " << s << std::endl;
        
        if (s < min) {
            std::cout << "warning " << s << " " << min << std::endl;
        }
        j[0] = j0;

        double reps = to_double(dd_real(r) - dd_real(1.0)/dd_real(rr));
        
        for (size_t l=1; l<=maxl; l++) {
            j[l] *= s;

            j[l] += (j[l-1]-(one+T(l+1)*rr)*j[l])*reps; // correction due to inexact computing of 1/r
            
            if (doprint) {
                T R = j[l]/j[l-1];
                T Rtest = ggj(l,ggi)/ggj(l-1,ggi);
                T relerr = RELERR(R,Rtest);
                if (relerr > 0) {
                    std::cout << l << " RtestZ " << Rtest << " " << R << " " << ggi << " " << r << " " << ggr[ggi] << " " << ggj(l,ggi) << " " << ggj(l-1,ggi) << " " << relerr << std::endl;
                }
            }
        }
        //std::cout << s << std::endl;
        //std::cout << "NEW " <<  j << std::endl;
        //std::exit(0);
    }
    else {
        // const T expm2r = std::exp(-two*r);
        // j[0] = (one-expm2r)*half*rr;
        // j[1] = ((one+expm2r)*half - j[0])*rr;
        j[0] = half*rr;
        j[1] = (half - j[0])*rr;
        for (size_t l=1; l<maxl; l++) {
            j[l+1] = j[l-1] - T(2*l+1.0)*j[l]*rr;
        }
    }
    
    return j;
}

template <typename T>
std::vector<T> JsphericalscaledVecDumbBAD(size_t maxl, T r) {
    const T zero = 0;
    std::vector<T> j(maxl+1,zero);
    if (r == zero) return j;
    
    const T eps = std::numeric_limits<T>::epsilon();
    const T min = std::numeric_limits<T>::min();
    const T one = 1;
    const T two = 2;
    const T half = T(0.5);
    const T rr = one/r;
    T R;

    bool doprint = false; //std::abs(r-T(9.536743e-07)) < 1e-11;
    if (doprint) std::cout << "r " << r << " eps " << eps << std::endl;

    T rcut = 20.0; // So that exp(-2*rcut) is much less than epsilon including gradual underflow
    if constexpr      (std::is_same<T,float>::value)   rcut = 10.0;
    else if constexpr (std::is_same<T,double>::value)  rcut = 20.0;
    else if constexpr (std::is_same<T,dd_real>::value) rcut = 40.0;
    else if constexpr (std::is_same<T,qd_real>::value) rcut = 80.0;

    if (r < std::max(rcut,T(12.5*maxl))) {
        size_t n;
        if constexpr      (std::is_same<T,float>::value)   n = maxl + 20;
        else if constexpr (std::is_same<T,double>::value)  n = maxl + 80;
        else if constexpr (std::is_same<T,dd_real>::value) n = maxl + 130;
        else if constexpr (std::is_same<T,qd_real>::value) n = maxl + 200;


        {
            double s = (2*n+1)/to_double(r);
            R = T((std::sqrt(s*s + 4) - s)*0.5);
        }
        
        // while (--n>maxl) {
        //     R = one/(R + (2*n+one)*rr); // Compute R_n from R_{n+1}
        // }

        // Reciprocal is slow so write back again in terms of function values instead of ratios
        {
            T j1 = min;
            T j0 = j1/R;
            while (--n>maxl) {
                T j = j1 + ((2*n+one)*rr)*j0;
                j0 = j1;
                j1 = j;
                if (j1 > one) {
                    j0 = std::ldexp(j0,-30);
                    j1 = std::ldexp(j1,-30);
                }
            }
            R = j1 / j0;
            std::cout << "TTT " << maxl << " " << r << " " << R << std::endl;
        }

        {
            double s = (2*n+1)/to_double(r);
            R = T((std::sqrt(s*s + 4) - s)*0.5);
        }
        
        while (--n>maxl) {
            R = one/(R + (2*n+one)*rr); // Compute R_n from R_{n+1}
        }
        std::cout << "YYY " << maxl << " " << r << " " << R << std::endl;
        
        
        // This helps avoid underflow but loses precision when computing j1/R
        //T j1 = min * std::max(one, R*r/T(2*maxl+one)); //was R ... now trying to avoid gradual underflow in s below
        //T j0 = j1/R; // was 1

        T j1 = R;
        T j0 = one;
        if constexpr (std::is_same<T,dd_real>::value || std::is_same<T,qd_real>::value) {
            double small = std::exp2(-128.0);
            j1 = mul_pwr2(R,small); // exact power of 2 so easy and exact to compute and apply
            j0 = small;
            // //double small = std::exp2(-128.0);
            // j1 = std::ldexp(j1,-128); //mul_pwr2(R,small); // exact power of 2 so easy and exact to compute and apply
            // j0 = std::exp2(-128); //small;
        }
        
        j[maxl] = j0;

        for (int l=maxl; l>0; l--) {
            j[l-1] = j1 + ((2*l+one)*rr)*j0;
            j1 = j0;
            j0 = j[l-1];
        }

        j0 = -std::expm1(-two*r)*half*rr; // for small r use expm1 to avoid cancellation

        T s = j0/j[0];
        
        if (s < min) {
            std::cout << "warning: s has underflowed " << s << " " << min << std::endl;
        }
        j[0] = j0;

        // Newton correction seems to help a lot for double, a little for float, but not at all for dd_real.
        T reps = 0;
        if constexpr      (std::is_same<T,float>::value)   reps = float(double(r) - double(1.0)/double(rr));
        else if constexpr (std::is_same<T,double>::value)  reps = to_double(dd_real(r) - dd_real(1.0)/dd_real(rr));
        else if constexpr (std::is_same<T,dd_real>::value) reps = to_dd_real(qd_real(r) - qd_real(1.0)/qd_real(rr));

        for (size_t l=1; l<=maxl; l++) {
            j[l] *= s;
            if constexpr (!std::is_same<T,qd_real>::value) {
                j[l] += (j[l-1]-(one+(l+one)*rr)*j[l])*reps; // Newton correction due to inexact computing of 1/r
            }

            // auto RELERR = [&](const T& a, const T&b) {return (std::abs(a-b)/std::abs(b))/eps;};

            // if (doprint) {
            //     T R = j[l]/j[l-1];
            //     T Rtest = ggj(l,ggi)/ggj(l-1,ggi);
            //     T Rrelerr = RELERR(R,Rtest);
            //     T jrelerr = RELERR(j[l],ggj(l,ggi));
            //     if (Rrelerr > 0) {
            //         std::cout << l << " RtestZ " << Rtest << " " << R << " " << Rrelerr << " " << j[l] << " " << ggj(l,ggi) << " " << jrelerr << std::endl;
            //     }
            // }
            
        }
    }
    else {
        j[0] = half*rr;
        j[1] = (half - j[0])*rr;
        for (size_t l=1; l<maxl; l++) {
            j[l+1] = j[l-1] - (2*l+one)*j[l]*rr;
        }
    }
    
    return j;
}

// Validated for l up to 32 and r from 1e-8 to 1e8
template <typename T>
std::vector<T> JsphericalscaledVecDumb(size_t maxl, T r) {
    const T zero = 0;
    std::vector<T> j(maxl+1,zero);
    if (r == zero) return j;
    
    const T eps = std::numeric_limits<T>::epsilon();
    const T min = std::numeric_limits<T>::min();
    const T one = 1;
    const T two = 2;
    const T half = T(0.5);
    const T rr = one/r;
    T R;

    bool doprint = false; //std::abs(r-T(9.536743e-07)) < 1e-11;
    if (doprint) std::cout << "r " << r << " eps " << eps << std::endl;

    T rcut = 20.0; // So that exp(-2*rcut) is much less than epsilon including gradual underflow
    if constexpr      (std::is_same<T,float>::value)   rcut = 10.0;
    else if constexpr (std::is_same<T,double>::value)  rcut = 20.0;
    else if constexpr (std::is_same<T,dd_real>::value) rcut = 40.0;
    else if constexpr (std::is_same<T,qd_real>::value) rcut = 80.0;

    if (r < std::max(rcut,T(12.5*maxl))) {
        size_t n;
        if constexpr      (std::is_same<T,float>::value)   n = maxl + 22;
        else if constexpr (std::is_same<T,double>::value)  n = maxl + 80;
        else if constexpr (std::is_same<T,dd_real>::value) n = maxl + 130;
        else if constexpr (std::is_same<T,qd_real>::value) n = maxl + 204;

        // Downward recursion to get R_n for n > maxl using ratios which is numerically stable with no overflow but slow due to reciprocals
        // {
        //     T s = (2*n+one)/r;
        //     R = (std::sqrt(s*s + T(4)) - s)*half; // note fix below for small r cancellation
        //     while (--n>maxl) {
        //         R = one/(R + (2*n+one)*rr); // Compute R_n from R_{n+1}
        //     }
        // }

        {
            // Since reciprocal is slow so rewrite again in terms of function values instead of ratios but now
            // have to manually handle possibility of overflow.  Tile the loop so that don't need to have
            // if test in the inner lopp ... this will be important in the vectorized version of this code.

            //T s = (2*n+1)*rr;
            //T RR = (std::sqrt(s*s + T(4)) - s + eps)*half; // note the eps to avoid RR=0 if s*s > 4/eps but this is still wrong
            T RR = (r*r + (n + 1)*r)/(r*r + (2*n + 1)*r + (2*n + 1)*(n + 1)); // initial guess that is slightly less accurate but no cancellation

            //T RR = (2*r*r*r + (3*n +3)*r*r + (2*n + 3)*(n + 1)*r)/(2*r*r*r + (5*n + 3)*r*r + (6*n + 3)*(n + 1)*r + (2*n + 3)*(2*n + 1)*(n + 1)); // no better
            
            T j1 = min;
            T j0 = min/RR;

            int step = 30;
            if constexpr (std::is_same<T,float>::value) step = 3;
            for (int i=n-1; i>int(maxl); i-=step) { // use int to avoid size_t wraparound
                size_t topl = std::max(i-step, int(maxl));
                for (size_t m=i; m>topl; m--) {
                    T j = j1 + ((2*m+1)*rr)*j0;
                    j1 = j0;
                    j0 = j;
                }
                if (j1 > one) {
                    j0 *= min;
                    j1 *= min; 
                }
            }
            RR = j1 / j0;
            //std::cout << "TTT " << maxl << " " << r << " " << R << " " << RR << " " << (R-RR)/eps << std::endl;
            R = RR;
        }
        
        // This helps avoid underflow but loses precision when computing j1/R
        //T j1 = min * std::max(one, R*r/T(2*maxl+one)); //was R ... now trying to avoid gradual underflow in s below
        //T j0 = j1/R; // was 1

        T j1 = R;
        T j0 = one;
        if constexpr (std::is_same<T,dd_real>::value || std::is_same<T,qd_real>::value) {
            j1 = std::ldexp(R,-128); 
            j0 = std::exp2(-128);
        }
        
        j[maxl] = j0;

        for (int l=maxl; l>0; l--) {
            j[l-1] = j1 + ((2*l+1)*rr)*j0;
            j1 = j0;
            j0 = j[l-1];
        }

        j0 = -std::expm1(-two*r)*half*rr; // for small r use expm1 to avoid cancellation

        T s = j0/j[0];
        
        if (s < min) {
            std::cout << "warning: s has underflowed " << s << " " << min << std::endl;
        }
        j[0] = j0;

        // Newton correction seems to help a lot for double, a little for float, but not at all for dd_real.
        T reps = 0;
        if constexpr      (std::is_same<T,float>::value)   reps = float(double(r) - double(1.0)/double(rr));
        else if constexpr (std::is_same<T,double>::value)  reps = to_double(dd_real(r) - dd_real(1.0)/dd_real(rr));
        else if constexpr (std::is_same<T,dd_real>::value) reps = to_dd_real(qd_real(r) - qd_real(1.0)/qd_real(rr));

        for (size_t l=1; l<=maxl; l++) {
            j[l] *= s;
            if constexpr (!std::is_same<T,qd_real>::value) {
                j[l] += (j[l-1]-(one+(l+one)*rr)*j[l])*reps; // Newton correction due to inexact computing of 1/r
            }

            // auto RELERR = [&](const T& a, const T&b) {return (std::abs(a-b)/std::abs(b))/eps;};

            // if (doprint) {
            //     T R = j[l]/j[l-1];
            //     T Rtest = ggj(l,ggi)/ggj(l-1,ggi);
            //     T Rrelerr = RELERR(R,Rtest);
            //     T jrelerr = RELERR(j[l],ggj(l,ggi));
            //     if (Rrelerr > 0) {
            //         std::cout << l << " RtestZ " << Rtest << " " << R << " " << Rrelerr << " " << j[l] << " " << ggj(l,ggi) << " " << jrelerr << std::endl;
            //     }
            // }
            
        }
    }
    else {
        j[0] = half*rr;
        j[1] = (half - j[0])*rr;
        for (size_t l=1; l<maxl; l++) {
            j[l+1] = j[l-1] - (2*l+one)*j[l]*rr;
        }
    }
    
    return j;
}

template <typename T>
std::vector<T> JsphericalscaledVecDumbGoodReciprocal(size_t maxl, T r) {
    const T zero = 0;
    std::vector<T> j(maxl+1,zero);
    if (r == zero) return j;
    
    const T eps = std::numeric_limits<T>::epsilon();
    const T min = std::numeric_limits<T>::min();
    const T one = 1;
    const T two = 2;
    const T half = T(0.5);
    const T rr = one/r;
    T R;

    bool doprint = false; //std::abs(r-T(9.536743e-07)) < 1e-11;
    if (doprint) std::cout << "r " << r << " eps " << eps << std::endl;

    T rcut = 20.0; // So that exp(-2*rcut) is much less than epsilon including gradual underflow
    if constexpr      (std::is_same<T,float>::value)   rcut = 10.0;
    else if constexpr (std::is_same<T,double>::value)  rcut = 20.0;
    else if constexpr (std::is_same<T,dd_real>::value) rcut = 40.0;
    else if constexpr (std::is_same<T,qd_real>::value) rcut = 80.0;

    if (r < std::max(rcut,T(12.5*maxl))) {
        size_t n;
        if constexpr      (std::is_same<T,float>::value)   n = maxl + 20;
        else if constexpr (std::is_same<T,double>::value)  n = maxl + 80;
        else if constexpr (std::is_same<T,dd_real>::value) n = maxl + 130;
        else if constexpr (std::is_same<T,qd_real>::value) n = maxl + 200;

        {
            T s = (2*n+one)/r;
            R = (std::sqrt(s*s + T(4)) - s)*half;
        }
        
        while (--n>maxl) {
            R = one/(R + (2*n+one)*rr); // Compute R_n from R_{n+1}
        }

        // This helps avoid underflow but loses precision when computing j1/R
        //T j1 = min * std::max(one, R*r/T(2*maxl+one)); //was R ... now trying to avoid gradual underflow in s below
        //T j0 = j1/R; // was 1

        T j1 = R;
        T j0 = one;
        if constexpr (std::is_same<T,dd_real>::value || std::is_same<T,qd_real>::value) {
            double small = std::exp2(-128.0);
            j1 = mul_pwr2(R,small); // exact power of 2 so easy and exact to compute and apply
            j0 = small;
        }
        
        j[maxl] = j0;

        for (int l=maxl; l>0; l--) {
            j[l-1] = j1 + ((2*l+one)*rr)*j0;
            j1 = j0;
            j0 = j[l-1];
        }

        j0 = -std::expm1(-two*r)*half*rr; // for small r use expm1 to avoid cancellation

        T s = j0/j[0];
        
        if (s < min) {
            std::cout << "warning: s has underflowed " << s << " " << min << std::endl;
        }
        j[0] = j0;

        // Newton correction seems to help a lot for double, a little for float, but not at all for dd_real.
        T reps = 0;
        if constexpr      (std::is_same<T,float>::value)   reps = float(double(r) - double(1.0)/double(rr));
        else if constexpr (std::is_same<T,double>::value)  reps = to_double(dd_real(r) - dd_real(1.0)/dd_real(rr));
        else if constexpr (std::is_same<T,dd_real>::value) reps = to_dd_real(qd_real(r) - qd_real(1.0)/qd_real(rr));

        for (size_t l=1; l<=maxl; l++) {
            j[l] *= s;
            if constexpr (!std::is_same<T,qd_real>::value) {
                j[l] += (j[l-1]-(one+(l+one)*rr)*j[l])*reps; // Newton correction due to inexact computing of 1/r
            }

            // auto RELERR = [&](const T& a, const T&b) {return (std::abs(a-b)/std::abs(b))/eps;};

            // if (doprint) {
            //     T R = j[l]/j[l-1];
            //     T Rtest = ggj(l,ggi)/ggj(l-1,ggi);
            //     T Rrelerr = RELERR(R,Rtest);
            //     T jrelerr = RELERR(j[l],ggj(l,ggi));
            //     if (Rrelerr > 0) {
            //         std::cout << l << " RtestZ " << Rtest << " " << R << " " << Rrelerr << " " << j[l] << " " << ggj(l,ggi) << " " << jrelerr << std::endl;
            //     }
            // }
            
        }
    }
    else {
        j[0] = half*rr;
        j[1] = (half - j[0])*rr;
        for (size_t l=1; l<maxl; l++) {
            j[l+1] = j[l-1] - (2*l+one)*j[l]*rr;
        }
    }
    
    return j;
}


template <typename T>
T JsphericalscaledY(size_t l, T r) {
    return JsphericalscaledVecDumb(32, r)[l];
}


// Computes exp(-r)*j(l,r) using Miller's algorithm and downward recursion.
// In Maple:
// j := (l, x) -> Re(I^(-l)*sqrt(-1/2*I*Pi/x)*BesselJ(l + 1/2, x*I))
template <typename T>
T JsphericalscaledX(size_t l, T r)
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
T Rlarger(size_t l, T r) {
    T one = T(1);
    T expm2r = std::exp(-T(2)*r);
    T rr = one/r;
    T R = (r-one+(r+one)*expm2r)/(r*(one-expm2r));
    size_t n = 1;
    while (n < l) {
        R = T(1)/R - T(2*n+1.0)*rr;
        n += 1;
    }
    return R;
}

template <typename T>
T Rguess(size_t l, T r) {
    //if (r < 0.5*T(double(l))) {
    double rcut = 10.0;
    if (r < T(std::min(rcut, 0.5*l))) {
        return Rsmallr(l+1,r);
    }
    else if (r > T(rcut)) {
        return Rlarger(l+1,r);
    }
    else {
        T s = T(2*l+1.0)/r;
        return (std::sqrt(s*s + T(4)) - s)*T(0.5);
    }
}

// Computes R_l using Miller's algorithm and downward recursion starting from l=n
template <typename T>
T Rdown(size_t n, size_t l, T r) {
    // R = \frac{-\frac{2n+3}{r} + \sqrt{\left(\frac{2n+3}{r}\right)^2 + 4}}{2}
    T rr = T(1)/r;

    T R = Rguess(n,r);
    while (n >= l) {
        R = T(1.0)/(R + T(2*n+1.0)*rr); // Compute R_n from R_{n+1}
        n--;
    }
    return R;
}

// For given l and r, compute the required starting value for downward recursion to obtain full precision in R_l(r) = j_l(r)/j_{l-1}(r)
// initializing with the guess described in the notes.
template <typename T>
size_t miller_start(size_t l, T r)
{
    T eps=std::numeric_limits<T>::epsilon();
    T R;

    // Find upper bound for n stepping forward by nstep
    size_t nstep = 10;
    size_t n = l+nstep;
    while (1) {
        R = Rdown(n,l,r);
        T R1 = Rdown(n+nstep,l,r);
        T err = std::abs(R-R1);
        std::cout << "n=" << n << " R=" << to_str(R) << " R1=" << to_str(R1) << " err=" << to_str(err) << std::endl;
        if (err < eps*R) break;
        n += nstep;
    }
    T Rbest = R;

    // Tighten the bound stepping back by 1
    while (n>l) {
        T R1 = Rdown(n-1,l,r);
        T err = std::abs(R1-Rbest);
        std::cout << "n=" << n << " R=" << to_str(R) << " R1=" << to_str(R1) << " err=" << to_str(err/R) << std::endl;
        if (err > eps*R) break;
        R = R1;
        n = n-1;
    }

    std::cout << "n=" << n << " R=" << to_str(R) << " R-Rguess" << to_str((R-Rguess(l,r))/R) << std::endl;
    
    return n;
}


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

    //Matrix<T> h(maxL+1, nR), j(maxL+1, nR);
    Matrix<T> h(maxL+2, nR), j(maxL+2, nR);
    std::vector<T> r(nR);

    std::string rs, js, hs;
    
    // for each l and for each r read l, r, j, h
    for (size_t l=0; l<=maxL+1; l++) {
        for (size_t i=0; i<nR; i++) {
            size_t ll;
            f >> ll >> rs >> js >> hs;
            r[i] = from_str<T>(rs);
            j(l,i) = from_str<T>(js);
            if (j(l,i) == 0) {
                std::cout << "j(" << l << "," << i << ") is zero" << std::endl;
            }
            h(l,i) = from_str<T>(hs);
            if (l != ll) throw std::runtime_error("l mismatch");
        }
    }

    // Perform index sort on r
    std::vector<size_t> idx(nR);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](size_t i, size_t j) {return r[i] < r[j];});

    // Reorder r, j, h
    std::vector<T> r1(nR);
    Matrix<T> j1(maxL+2, nR), h1(maxL+2, nR);
    for (size_t i=0; i<nR; i++) {
        r1[i] = r[idx[i]];
        for (size_t l=0; l<=maxL+1; l++) {
            j1(l,i) = j(l,idx[i]);
            h1(l,i) = h(l,idx[i]);
        }
    }

    return std::make_tuple(r1, j1, h1);
    //return std::make_tuple(r, j, h);
}

template <typename T>
void test_bessel() {
    std::vector<T> r;
    Matrix<T> j, h;
    std::tie(r, j, h) = load_bessel_test_data<T>();
    const size_t nR=r.size(), maxL=j.rows()-1;
    const T eps = std::numeric_limits<T>::epsilon();
    
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
                if (relerr > fudge*eps) {
                    std::cout << "l=" << l << " r=" << r[i] << " j=" << to_str(j(l,i)) << " j0=" << to_str(j0) << " err=" << err << " relerr=" << relerr/eps << " " << std::endl;
                }
            }
        }
    }
}

template <typename T>
void test_bessel2() {
    std::vector<T> r;
    Matrix<T> j, h;
    std::tie(r, j, h) = load_bessel_test_data<T>();
    const size_t nR=r.size(), maxL=j.rows()-2; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! was -1
    const T eps = std::numeric_limits<T>::epsilon();

    T maxabserr = 0;
    T avgabserr = 0;
    T avgsgnerr = 0;
    T count = 0;
    T worstr = 0;
        
    //std::cout << "r " << r;
    
    for (size_t i=0; i<nR; i+=1) {
        bool doprint = false; //std::abs(r[i]-T(9.536743e-07)) < 1e-11;
        ggi = i;
        //std::cout << "r " << i << " " << r[i] << std::endl;
        std::vector<T> js = JsphericalscaledVecDumb(maxL, r[i]);
        for (size_t l=0; l<=maxL; l+=1) {
            //for (size_t l=40; l<=40; l+=5) {
            //for (size_t i=300; i<=300; i+=20) {

            //double reps = to_double(dd_real(r[i]) - dd_real(1.0)/dd_real(1.0/r[i]));
            
            T j0 = js[l];
            T err = (j0-j(l,i));
            T relerr = (err/j(l,i))/eps;
            //T estrelerr = ((l/r[i])*reps)/eps;
            //T estrelerr2 = 0.0;
            //if (l > 0) {
            //estrelerr2 = (j(l-1,i)-(1+(l+1)/r[i])*j(l,i))*reps/(j(l,i)*eps);
            //            }

            if (j(l,i)>std::numeric_limits<T>::min()) { // really tiny values will be denormed so ignore them
                avgsgnerr += relerr;
                err = abs(err);
                avgabserr += std::abs(relerr);
                if (std::abs(relerr) > maxabserr) {
                    maxabserr = std::abs(relerr);
                    worstr = r[i];
                }
                maxabserr = std::max(maxabserr,std::abs(relerr));
                count += 1;
                // double gave 26eps with worst errors at short and intermediate distances ... this was due to test values of r not being representable exactly in doubles
                // 26 goes to 14 from using representable r, then 14 to 6 from gradient correction for err in computing 1/r, then to 4 from not doing 1/R before starting downward recursion
                
                // qd gives 2eps
                // dd gives better than 20eps nearly everywhere with worst about 46
                T fudge = 1;
                if constexpr (std::is_same_v<T, float>) fudge = 2;
                if constexpr (std::is_same_v<T, double>) fudge = 5;
                else if constexpr (std::is_same_v<T, dd_real>) fudge = 47;
                else if constexpr (std::is_same_v<T, qd_real>) fudge = 2;
                
                if (doprint || relerr > fudge) {
                    std::cout << "l=" << l << " r=" << r[i] << " j=" << to_str(j(l,i)) << " j0=" << to_str(j0) << " err=" << err << " relerr=" << relerr << std::endl;
                }
            }
        }
    }
    std::cout << "test bessel2: " << typeid(T).name() << ": " <<  avgsgnerr/count << " " << avgabserr/count << " " << maxabserr << " eps" << " " << worstr << std::endl;
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

    // for (int i=-330; i<=330; i+=1) {
    //     char buf[100];
    //     sprintf(buf, "1e%d", i);
    //     double x = atof(buf);
    //     double y = from_str<double>(buf);
    //     std::cout << i << " " << x << " " << y << std::endl;
    // }
    // return 0;
    
    unsigned int oldcw; // For x86
    fpu_fix_start(&oldcw);

    {
        // static initialize of dd/qd does not work ... need to do it
        // in main code after FPU fix, and thread-safe initialization
        // of dd/qd package
        factorial_cache<double> f;
        factorial_cache<dd_real> fd;
        factorial_cache<qd_real> fq;
        factorial_cache<size_t> fs;
        double_factorial_cache<double> df;
        double_factorial_cache<dd_real> dfd;
        double_factorial_cache<qd_real> dfq;
        double_factorial_cache<size_t> dfs;

        for (size_t i=0; i<21; i++) std::cout << i << " " << factorial<size_t>(i) << " " << double_factorial<size_t>(i) << std::endl;
    }

    // std::cout << "7! " << factorial<double>(7) << std::endl;
    // std::cout << "7! " << factorial<dd_real>(7) << std::endl;
    // std::cout << "7! " << factorial<qd_real>(7) << std::endl;

    // std::cout << "7!! " << double_factorial<double>(7) << std::endl;
    // std::cout << "7!! " << double_factorial<dd_real>(7) << std::endl;
    // std::cout << "7!! " << double_factorial<qd_real>(7) << std::endl;

    
    // float x = 0.001f;
    // for (int i=0; i<10; i++) {
    //     std::cout << x << " " << log10_estimate(x) << std::endl;
    //     x *= 10;
    // }

    // Small r is when r < n ... so Taylor series for R is rapidly convergent as is downward recursion
    // Large r is when r > 5*n ... so upwrard recursion for R is sufficiently accurate
    size_t L = 48;
    //auto f = [=](double r){return Jsphericalscaled(L, r)/Jsphericalscaled(L-1, r);};
    auto f33 = [=](double r){return r==0 ? 0 : to_double(Jsphericalscaled(L, dd_real(r))/Jsphericalscaled(L-1, dd_real(r)));};
    //BarycentricX<double> b(2000, 6, L, 12.5*L, f33);
    Barycentric<double> b(1250, 12, 0.5, 12.5*L, f33);
    b33 = b;

    // std::cout << b(20.018) << std::endl;
    // std::cout << f(20.018) << std::endl;
    // return 0;

    
    std::cout << "double " << std::numeric_limits<double>::epsilon() << " " << std::numeric_limits<double>::min() << std::endl;
    std::cout << b.maxerr(f33) << std::endl;


    std::tie(ggr, ggj, ggh) = load_bessel_test_data<dd_real>();
    
    test_bessel2<float>();
    test_bessel2<double>();
    test_bessel2<dd_real>();
    test_bessel2<qd_real>();

    return 0;
    
    // using T = qd_real;
    // //T x = from_str<T>("0.9");
    // //T tt = from_str<T>("0.00002382456605715974270203848124219941145285588565482387037934163044550358740304620994555378008103303005");
    // T x = from_str<T>("40.0");
    // T tt = from_str<T>("0.008555280151367187500000000000000000326580158489513959748694041681666003919215557451893157658276768298");
    // T told = Jsphericalscaled(5, x);
    // T tnew = JsphericalscaledX(5, x);
    // std::cout << to_str((told-tt)/tt) << std::endl;
    // std::cout << to_str((tnew-tt)/tt) << std::endl;
    // std::cout << to_str(tt) << std::endl;

    // return 0;
    
    // T step = 0.0625;
    // T r = 1.0;
    // while (r <= T(4.0)) {
    //     size_t n = miller_start(20, r);
    //     std::cout << "test " << to_str(r) << " " << n << std::endl;
    //     r += step;
    // }
    // return 0;

    

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
