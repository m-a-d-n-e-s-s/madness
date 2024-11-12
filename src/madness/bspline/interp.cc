#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

template <typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<T>& v) {
    s << "[ ";
    for (const auto& x : v) s << x << " ";
    s << "]";
    return s;
}

// Barycentric interpolation class using N uniformly spaced-intervals with m uniformly spaced points within each interval (including endpoints)
template <typename T>
class Barycentric {
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
    Barycentric(size_t N, size_t m, T a, T b, funcT f) 
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
        // std::cout << "input x = " << x << " " << N << " " << H << std::endl;
        assert(x>=a && x<b);
        x -= a;
        const size_t i = x/H;
        x -= i*H;
        // std::cout << "i = " << i << " " << x << std::endl;

        T num = 0;
        T den = 0;
        for (size_t j=0; j<m; j++) {
            T diff = x-j*h;
            if (diff == 0) return y[i*(m-1) + j];
            T wj = w[j]/diff;
            num += wj*y[i*(m-1)+j];
            den += wj;
        }
        return num/den;
    }

    template <typename funcT>
    T maxerr(funcT f, size_t oversample=5) {
        T testh = h/oversample;
        T err = 0;
        T x = a;
        while (x < b) {
            err = std::max(std::abs(f(x) - (*this)(x)), err);
            x += testh;
        }
        return err;
    }    
};

// Barycentric interpolation class using N uniformly spaced-intervals with m Chebyshev points in each interval
template <typename T>
class BarycentricChebyshev {
    size_t N; // Number of larger intervals
    size_t m; // Number of points in each interval
    T a; // Left endpoint of the interval
    T b; // Right endpoint of the interval
    T H; // Larger interval size

    std::vector<T> y; // Vector of m*N function values
    std::vector<T> p; // Vector of m sampling points in each interval
    std::vector<T> w; // Vector of m Barycentric weights

    std::vector<T> points(const size_t m) {
        size_t n = m-1;
        std::vector<T> p(m);
        for (size_t j=0; j<=n; j++) {
            p[j] = 0.5*H*(1-std::cos(M_PI*(j+0.5)/m));
        }
        return p;
    }

    std::vector<T> weights(const size_t m) {
        std::vector<T> w(m);
        size_t n = m-1;
        T sign = 1;
        for (size_t j=0; j<=n; j++) {
            w[j] = sign*std::sin(M_PI*(j+0.5)/m);
            sign = -sign;
        }
        return w;
    }

public:
    template <typename funcT>
    BarycentricChebyshev(size_t N, size_t m, T a, T b, funcT f) 
        : N(N)
        , m(m)
        , a(a)
        , b(b)
        , H((b-a)/N)
        , y(N*m)
        , p(points(m))
        , w(weights(m))
    {
        assert(m>=2);
        assert(N>=1);
        for (size_t i=0; i<N; i++) {
            for (size_t j=0; j<m; j++) {
                T x = a + i*H + p[j];
                y[i*m+j] = f(x);
            }
        }
        // std::cout << "H = " << H << std::endl;
        // std::cout << "y = " << y << std::endl;
        // std::cout << "p = " << p << std::endl;
        // std::cout << "w = " << w << std::endl;
        // std::cout << f(a) << " " << f(b) << std::endl;
    }

    T operator()(T x) {
        assert(x>=a && x<b);
        x -= a;
        const size_t i = x/H;
        x -= i*H;
        // std::cout << "i = " << i << " " << x << std::endl;

        T num = 0;
        T den = 0;
        for (size_t j=0; j<m; j++) {
            // std::cout << j << "  p[j] = " << p[j] << std::endl;
            T diff = x-p[j];
            if (diff == 0) return y[i*m + j];
            T wj = w[j]/diff;
            num += wj*y[i*m+j];
            den += wj;
        }
        return num/den;
    }

    template <typename funcT>
    T maxerr(funcT f, size_t oversample=5) {
        T testh = H/(oversample*m);
        T err = 0;
        T x = a;
        while (x < b) {
            err = std::max(std::abs(f(x) - (*this)(x)), err);
            x += testh;
        }
        return err;
    }    
};

int main() {
    auto f = [](double x) { return std::sin(x); };
    BarycentricChebyshev<double> c(10, 4, 0, 1, f);
    std::cout << c(0.51) << " " << std::sin(0.51) << std::endl;
    std::cout << c.maxerr(f) << std::endl;
    
    for (size_t N=1; N<=20; N++) {
        for (size_t m=2; m<=10; m++) {
            Barycentric<double> b(N, m, 0, 1, f);
            BarycentricChebyshev<double> c(N, m, 0, 1, f);
            std::cout << "N = " << N << " m = " << m << " maxerr = " << b.maxerr(f) << " " << c.maxerr(f) << std::endl;
        }
    }
    return 0;
}
