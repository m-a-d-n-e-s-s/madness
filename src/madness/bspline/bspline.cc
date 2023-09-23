#include <limits>
#include <memory>
#include <complex>
#include <algorithm>

#include <madness.h>
#include "gauleg.h"
#include <madness/misc/gnuplot.h>

// For indexing matrices MADNESS uses the notation A(i,j).
// For indexing vectors MADNESS uses the notation A(i).
// MADNESS indexes from 0.
// To index an inclusive range from a to b MADNESS uses the notation Slice(a,b).

// In this file building new classes to implement B-spline basis
// functions and knots and test them as radial basis functions for
// quantum chemistry applications.

template<typename T> struct is_complex : std::false_type {};
template<typename T> struct is_complex<std::complex<T>> : std::true_type {};
template<typename T> struct scalar_type {typedef T type; };
template<typename T> struct scalar_type<std::complex<T>> {typedef T type;};

using namespace madness; // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

// Knots classes generate (don't need to store) the knots and provide the interval function via a inline-able templated interface 
template <typename T>
class KnotsGeneral {
public:
    static constexpr const char* name = "general";
    typedef T value_type;

private:    
    const Tensor<T> _knots; // unique knot vector in ascending order
    const T xlo, xhi; // interval
    const size_t nknots; // number of knots

public:
    KnotsGeneral(const Tensor<T>& knots)
        : _knots(knots)
        , xlo(_knots[0])
        , xhi(_knots[knots.size()-1])
        , nknots(_knots.size())
    {
        for (size_t i=1; i<nknots; i++) {
            if (_knots[i] < _knots[i-1]) throw "KnotsGeneral: knots not in ascending order or not unique";
        }
    }

    size_t size() const {return nknots;}

    // Return (shallow copy) the knots
    const Tensor<T> knots() const {
        return _knots;
    }

    // Locate the leftmost index of the unique knot interval
    // containing x in the vector of unique knots. If x exactly equals
    // a knot return the knot index (unless it is the rightmost knot).
    //
    // This general version is dumb and slow ... linear search!!!!!  Bisection would be a better choice.
    size_t interval(T x) const {
#if !defined(NDEBUG)
        if (x < xlo || x > xhi) throw "interval: x not in range";
#endif
        for (size_t i=0; i<nknots-2; i++) {
            if (x >= _knots[i] && x < _knots[i+1]) return i;
        }
        return nknots-1;
    }
};

template <typename T>
class KnotsUniform {
public:
    static constexpr const char* name = "uniform";
    typedef T value_type;

private:
    const size_t nknots; // number of knots
    const T xlo, xhi; // interval
    const T h; // knot spacing
    const T rheps; // reciprocal of knot spacing

public:
    KnotsUniform(size_t nknots, T xlo, T xhi)
        : nknots(nknots)
        , xlo(xlo)
        , xhi(xhi)
        , h((xhi - xlo)/(nknots-1))
        , rheps(1.0/(h*(1+std::numeric_limits<T>::epsilon()))) // to avoid endpoint issues
    {}

    Tensor<T> knots() const {
        Tensor<T> pts(nknots);
        for (size_t i=0; i<nknots-1; i++) pts[i] = i*h + xlo;
        pts[nknots-1] = xhi; // to avoid rounding errors
        return pts;
    }

    size_t interval(T x) const {
#if !defined(NDEBUG)
        if (x < xlo || x > xhi) throw "Xinterval: x not in t";
#endif
        return size_t((x - xlo)*rheps);
    }

    size_t size() const {return nknots;}
};
    
// Geometric or exponential (logarithmic depending on your perspective) knot spacing in [xlo,xhi+xlo] shifted to zero on the left
// so resulting grid is [0,xhi].
template <typename T>
class KnotsGeometricShifted {
public:
    static constexpr const char* name = "geometric-shifted";
    typedef T value_type;

private:
    const size_t nknots; // number of knots
    const T xlo, xhi; // interval
    const T h; // logarithmic knot spacing
    const T rxlo; // reciprocal of xlo
    const T rlogh; // reciprocal of log(h)
public:
    KnotsGeometricShifted(size_t nknots, T xlo, T xhi)
        : nknots(nknots)
        , xlo(xlo)
        , xhi(xhi)
        , h(std::pow((xhi+xlo)/xlo,1.0/(nknots-1)))
        , rxlo(1.0/xlo)
        , rlogh((1.0-4*std::numeric_limits<T>::epsilon())/std::log(h)) // to avoid endpoint issues
    {
        print("KnotsGeometricShifted: h = ", h, " xlo = ", xlo, " xhi = ", xhi, " nknots = ", nknots, " rlogh = ", rlogh, " rxlo = ", rxlo);
    }

    Tensor<T> knots() const {
        Tensor<T> pts(nknots);
        pts[0] = 0;
        for (size_t i=1; i<nknots-1; i++) pts[i] = xlo*(std::pow(h,i) - 1);
        pts[nknots-1] = xhi; // to avoid roundoff error
        return pts;
    }

    size_t interval(T x) const {
#if !defined(NDEBUG)
        if (x < 0 || x > xhi) throw "Xinterval: x not in t";
#endif
        return size_t(log(x*rxlo + T(1))*rlogh);
    }
   
    size_t size() const {return nknots;}
};
    
// Chebyshev points *modified to include the endpoints*
template <typename T>
class KnotsChebyshev {
public:
    static constexpr const char* name = "chebyshev-modified";
    typedef T value_type;

private:
    const size_t nknots; // number of knots
    const T xlo, xhi; // interval
    const T h, rh; // grid spacing before cosine transform and its inverse
    const double rLhalf; // 2.0/(xhi-xlo)

public:
    KnotsChebyshev(size_t nknots, T xlo, T xhi)
        : nknots(nknots)
        , xlo(xlo)
        , xhi(xhi)
        , h(constants::pi/(nknots - 1)) // WHAT ABOUT dd/qd
        , rh(1.0/this->h) 
        , rLhalf((1.0-std::numeric_limits<T>::epsilon())*2.0/(xhi-xlo)) // to avoid endpoint issues
    {}
    
    // Generate the knots
    Tensor<T> knots() const {
        Tensor<T> pts(nknots);
        for (size_t i=0; i<nknots; i++) pts[nknots-1-i] = std::cos(i*h);
        pts = (pts + T(1.0))*(T(0.5)*(xhi-xlo)) + xlo;
        pts[0] = xlo; // to avoid roundoff error
        pts[nknots-1] = xhi; 
        return pts;
    }
    
    size_t interval(T x) const {
#if !defined(NDEBUG)
        if (x < xlo || x > xhi) throw "Xinterval: x not in t";
#endif
        T y = (x - xlo)*rLhalf - T(1);
        size_t k = size_t((nknots-T(1)) - std::acos(y)*rh);
        return k;
    }

    size_t size() const {return nknots;}
};
    
// Oversample the knots by inserting geometric mean of neighboring knots.  For knots on a geometric/exponential grid this is
// what you want, and for a uniform grid we have sqrt(a*(a+h))= a+h/2+O(h^2) where a is a knot and h is the knot spacing.
// So for a dense uniform grid it will be close to the geometric mean for a>>h.
template <typename T>
Tensor<T> oversample_knots(const Tensor<T>& knots) {
    long nknots = knots.size();
    Tensor<T> newknots(2*nknots - 1);
    for (size_t i=0; i<nknots-1; i++) {
        newknots[2*i] = knots[i];
        if (knots[i] == 0 || knots[i+1] == 0 || knots[i]*knots[i+1] < 0) {
            newknots[2*i+1] = 0.5*(knots[i] + knots[i+1]);
        }
        else {
            double sgn = (knots[i] > 0) ? 1 : -1;
            newknots[2*i+1] = sgn*std::sqrt(knots[i]*knots[i+1]);
        }
    }
    newknots[2*nknots-2] = knots[nknots-1];
    return newknots;
}

// Manages the nearly minimum amount of data to define the b-spline basis and to efficiently compute values and related matrices/operators
template <typename T, typename knotsT>
class BsplineBasis : protected knotsT {
    static_assert(std::is_same<T,typename knotsT::value_type>::value, "Bsplinebasis: T and knotsT::value_type must be the same");
    static_assert(std::is_floating_point<T>::value, "Bsplinebasis: T must be floating point");
    
    // Pad/repeat the knots at the ends to complete the basis
    static Tensor<T> pad_knots(size_t order, const Tensor<T>& knots) {
        long nknots = knots.size();
        Tensor<T> padded(nknots + 2*order - 2);
        for (size_t i=0; i<order-1; i++) padded[i] = knots[0];
        for (size_t i=0; i<nknots; i++) padded[i+order-1] = knots[i];
        for (size_t i=0; i<order-1; i++) padded[i+nknots+order-1] = knots[nknots-1];
        return padded;
    }

public:
    typedef T value_type;
    typedef knotsT knots_type;
    const size_t order; // the order of the basis
    const size_t p;     // the degree of the basis = order-1
    const size_t nknots;// the number of unique/unpadded knots
    const size_t npknots;// the number of padded knots
    const size_t nbasis;// the number of basis functions
    const Tensor<T> knots; // the unpadded knots
    Tensor<T> t; // the padded knots

    /// Construct a B-spline basis with the given order and (unpadded) knots.
    BsplineBasis(size_t order, const knotsT& knots)
        : knotsT(knots)
        , order(order)
        , p(order - 1)
        , nknots(knots.size())
        , npknots(nknots + 2*order - 2)
        , nbasis(nknots + order - 2)
        , knots(knots.knots())
        , t(pad_knots(order,knots.knots()))
    {
        if (order < 2) throw "BsplineBasis: order must be >= 2";
        if (order > 20) throw "BsplineBasis: order must be <= 20 (sanity check)";
        if (nknots < 2*order - 1) throw "BsplineBasis: knots must have at least 2*order - 1 elements";
        for (size_t i=1; i<nknots; i++) {
            if (this->knots[i] < this->knots[i-1]) throw "BsplineBasis: knots not in ascending order or not unique";
        }
        if (t.size() != npknots) throw "BsplineBasis: internal error with padded knots";
        if (npknots-order != nbasis) throw "BsplineBasis: internal error with number of basis functions";
    }

    // Returns the total number of basis functions corresponding to the length of a coefficient vector
    size_t size() const { return nbasis; }

    // Returns a const ref to the knots object
    const knotsT& knots_object() const {return *(knotsT*)this;}

    // in next few routines can k=interval(x)+p is usually used with k-p except where k-r.
    
    // Stably evaluate a spline expansion at a point, i.e., compute sum(i) c[i] b[i](x)
    template <typename R>
    R deBoor(T x, const Tensor<R>& c) const {
        MADNESS_ASSERT(c.dims()[0] == nbasis);
        size_t k = knotsT::interval(x) + p;
        Tensor<R> d = copy(c(Slice(k-p, k)));
        for (size_t r=0; r<p; r++) {
            for (size_t j=p; j>r; j--) {
                T a = (x - t[j+k-p])/(t[j+k-r] - t[j+k-p]);
                d[j] = (T(1)-a)*d[j-1] + a*d[j];
            }
        }
        return d[p];
    }

    // Stably evaluate a spline expansion at many points, i.e., compute f[j] = sum(i) c[i] b[i](x[j])
    // This is the unvectorized version.
    template <typename R>
    Tensor<R> deBoor_unvectorized(const Tensor<T>& x, const Tensor<R>& c) const {
        MADNESS_ASSERT(c.dims()[0] == nbasis);
        const size_t nx = x.dims()[0];
        std::vector<size_t> k(nx);
        Tensor<R> d(p+1,nx);
        for (size_t i=0; i<nx; i++) {
            k[i] = knotsT::interval(x[i]) + p;
            d(_,i) =  c(Slice(k[i]-p, k[i]));
        }

        for (size_t r=0; r<p; r++) {
            for (size_t j=p; j>r; j--) {
                for (size_t i=0; i<nx; i++) {
                    T a = (x[i] - t[j+k[i]-p])/(t[j+k[i]-r] - t[j+k[i]-p]);
                    d(j,i) = (T(1)-a)*d(j-1,i) + a*d(j,i);
                }
            }
        }
        return copy(d(p,_));
    }

    // Stably evaluate a spline expansion at many points, i.e., compute f[j] = sum(i) c[i] b[i](x[j])
    template <typename R>
    Tensor<R> deBoor(const Tensor<T>& x, const Tensor<R>& c) const {
        MADNESS_ASSERT(c.dims()[0] == nbasis);
        const size_t nx = x.dims()[0];
        Tensor<R> d(p+1,nx), s(2*p,nx); // zeroing not necessary 
        for (size_t i=0; i<nx; i++) { // non-sequential writes could be optimized?
            size_t k = knotsT::interval(x[i])+p;
            d(_,i) =  c(Slice(k-p, k));
            s(_,i) =  t(Slice(k-p+1, k+p));
        }

        for (size_t r=0; r<p; r++) {
            for (size_t j=p; j>r; j--) {
                for (size_t i=0; i<nx; i++) { // Does this really vectorize or do we need an IVDEP pragma?
                    T a = (x[i] - s(j-1,i))/(s(p+j-r-1,i) - s(j-1,i));
                    d(j,i) = (T(1)-a)*d(j-1,i) + a*d(j,i);
                }
            }
        }
        return copy(d(p,_));
    }

    // Evaluate the basis functions at the point x
    Tensor<T> bspline(T x) const {
        size_t k = knotsT::interval(x)+p;
        Tensor<T> bm1(nbasis);
        Tensor<T> b(nbasis);
        bm1[k] = T(1.0);
        for (size_t q=1; q<=p; q++) {
            size_t ilo = k - q, ihi = k;
            if (ilo + q <= order) {
                ilo = order - q;
                size_t i = ilo - 1;
                b[i] = (t[i+q+1] - x)*bm1[i+1]/(t[i+q+1] - t[i+1]);
            }
            if (k == nbasis-1) {
                ihi = k - 1;
                size_t i = k;
                b[i] = (x - t[i])*bm1[i]/(t[i+q] - t[i]);
            }
            for (size_t i=ilo; i<=ihi; i++) {
                b[i] = (x - t[i])*bm1[i]/(t[i+q] - t[i]) + (t[i+q+1] - x)*bm1[i+1]/(t[i+q+1] - t[i+1]);
            }
            std::swap(b, bm1);
        }
        return bm1;
    }

    // Make the matrix of bsplines evaluated at sample points A[i,j] = b[j](Xsample[i]) 
    Tensor<T> tabulate_basis(const Tensor<T>& Xsample) const {
        size_t nsample = Xsample.dims()[0];
        Tensor<T> A(nsample, nbasis);
        for (size_t i=0; i<nsample; i++) {
            A(i,_) = bspline(Xsample[i]);
        } 
        return A;
    }

    // Make the matrix that does LSQ fit from values on sample points (which must be strictly within basis support) to basis coefficients (c).
    // nzeroL/R are the number of basis functions to be neglected on the left/right side of the interval.
    // nzero=0 (default) fits all basis functions.
    // nzero=1 sets the value on the boundary to zero.
    // nzero=2 sets the value and first derivative on the boundary to zero.
    // highest permitted value is order (degree+1)
    // This is the pseudo-inverse of A[i,j] = b[j](Xsample[i]) computed using the non-zero singular values.
    Tensor<T> make_lsq_matrix(const Tensor<T>& Xsample, size_t nzeroL=0, size_t nzeroR=0) const {
        size_t nsample = Xsample.dims()[0];
        if (nsample < nbasis) throw "make_lsq_matrix: no. of samples must be >= no. of basis functions";
        MADNESS_ASSERT(nzeroL <= order);
        MADNESS_ASSERT(nzeroR <= order);
        
        Tensor<T> A = tabulate_basis(Xsample);
        A = A(_,Slice(nzeroL,nbasis-nzeroR-1));
        size_t nbasnew = nbasis-nzeroL-nzeroR;
        Tensor<T> u, s, vt;
        svd(A, u, s, vt);
        for (size_t i=0; i<nbasnew; i++) {
            T fac = 1.0/s[i];
            for (size_t j=0; j<nbasnew; j++) {
                vt(i,j) *= fac;
            }
        }
        Tensor<T> M(nbasis,nsample);
        M(Slice(nzeroL,nbasis-nzeroR-1),_) = transpose(inner(u,vt));
        return M;
    }

    template <typename R>
    Tensor<R> fit(const Tensor<T>& Xsample, const Tensor<R>& Ysample, size_t nzeroL=0, size_t nzeroR=0) const {
        MADNESS_ASSERT(Xsample.dims()[0] == Ysample.dims()[0]);
        Tensor<T> M = make_lsq_matrix(Xsample,nzeroL,nzeroR);
        return inner(M,Ysample);
    }

    // Compute the "exact" derivative of a bspline expansion in terms of splines one degree lower.
    // Note to evaluate the result need to constuct a basis of degree p-1.
    Tensor<T> deriv_exact(const Tensor<T>& c) const {
        MADNESS_ASSERT(c.dims()[0] == nbasis);
        Tensor<T> d(nbasis - 1);
        for (size_t j=0; j<nbasis - 1; j++) {
            d[j] = p*(c[j + 1] - c[j])/(t[j + p + 1] - t[j + 1]);
        }
        return d;
    }

    // Make the (bidiagonal) matrix represenation of the operator that computes the exact derivative in terms of splines one degree lower
    Tensor<T> make_deriv_exact_matrix() const {
        Tensor<T> dMat(nbasis - 1, nbasis);
        for (size_t j=0; j<nbasis - 1; j++) {
            T s = p/(t[j + p + 1] - t[j + 1]);
            dMat(j, j + 1) = s;
            dMat(j, j) = -s;
        }
        return dMat;
    }

    // Make the matrix that applies the derivative after projecting into a one higher order basis.
    // Gives us 1 order higher accuracy and produces a result in the same order basis as the input.
    Tensor<T> make_deriv_matrix(const Tensor<T>& Xsample) const {
        size_t nsample = Xsample.dims()[0];
        BsplineBasis<T,knotsT> B1 = BsplineBasis(order + 1, knots_object());
        Tensor<T> A = tabulate_basis(Xsample);
        Tensor<T> M = B1.make_lsq_matrix(Xsample);
        Tensor<T> dMat = B1.make_deriv_exact_matrix();
        return inner(inner(dMat,M),A);
    }

    // Make the matrix that applies the exact derivative and projects back into the original order basis.
    Tensor<T> make_deriv_matrixX(const Tensor<T>& Xsample) const {
        size_t nsample = Xsample.dims()[0];
        BsplineBasis<T,knotsT> B1 = BsplineBasis(order - 1, knots_object());
        Tensor<T> A = B1.tabulate_basis(Xsample);
        Tensor<T> M = make_lsq_matrix(Xsample);
        Tensor<T> dMat = make_deriv_exact_matrix();
        return inner(inner(M,A),dMat);
    }

    // Make the matrix that applies the second derivative after projecting into a two higher order basis.
    // Gives us 2 order higher accuracy and produces a result in the same order basis as the input.
    Tensor<T> make_deriv2_matrix(const Tensor<T>& Xsample) const {
        size_t nsample = Xsample.dims()[0];
        BsplineBasis<T,knotsT> B1 = BsplineBasis(order + 1, knots_object());
        BsplineBasis<T,knotsT> B2 = BsplineBasis(order + 2, knots_object());
        Tensor<T> A = tabulate_basis(Xsample);
        Tensor<T> M2 = B2.make_lsq_matrix(Xsample);
        Tensor<T> dMat2 = B2.make_deriv_exact_matrix();
        Tensor<T> dMat1 = B1.make_deriv_exact_matrix();
        return inner(inner(inner(dMat1,dMat2),M2),A);
    }

    // Make the matrix that applies the exact second derivative and then projects into a two higher order basis.
    // Produces a result in the same order basis as the input.
    Tensor<T> make_deriv2_matrixX(const Tensor<T>& Xsample) const {
        size_t nsample = Xsample.dims()[0];
        BsplineBasis<T,knotsT> B1 = BsplineBasis(order - 1, knots_object());
        BsplineBasis<T,knotsT> B2 = BsplineBasis(order - 2, knots_object());
        Tensor<T> A = B2.tabulate_basis(Xsample);
        Tensor<T> M = make_lsq_matrix(Xsample);
        Tensor<T> dMat2 = B1.make_deriv_exact_matrix();
        Tensor<T> dMat1 = make_deriv_exact_matrix();
        return inner(M,inner(A,inner(dMat2,dMat1)));
    }

    // Make the quadrature points and weights for a rule that employs n GL points between each unique knot.  So we will have n*(nknot-1) points.
    // Returns a pair of tensors, the first is the quadrature points, the second is the quadrature weights.
    //
    // Rationale for choosing n when computing matrix elements starts
    // from GL order n exactly integrating a polynomial of degree
    // 2n-1 and any function is presumed already resolved into the
    // b-spline basis.
    // 1. For <bi|f|bj> choose n >= (3p+1)/2  (or (3p+3)/2 with r**2 weight)
    // 2. For <bi|bj> choose n >= (2p+1)/2  (or (2p+3)/2 with r**2 weight)
    // 3. For <bi|bj"> chose n >= (2p-1)/2  (or (2p+1)/2 with r**2 weight)
    std::pair<Tensor<T>,Tensor<T>> make_spline_quadrature(size_t n) const {
        GaussLegendre<T> g(n);
        Tensor<T> x(n), w(n);
        for (size_t i=0; i<n; i++) {
            x[i] = g.pts()[i];
            w[i] = g.wts()[i];
        }
        print("x", x);
        Tensor<T> X(n*(nknots-1)), W(n*(nknots-1));
        for (size_t knot=0; knot<nknots-1; knot++) {
            T a=knots[knot], b=knots[knot+1];
            for (size_t i=0; i<n; i++) {
                X[i + n*knot] = T(0.5)*(a+b) + T(0.5)*(b-a)*x[i];
                W[i + n*knot] = T(0.5)*(b-a)*w[i];
            }
        }
        return {X, W};
    }

    // Test the basic class functionality --- this is not yet a unit test
    static bool test() {
        size_t nknots = 50;
        size_t order = 12; 
        T xlo = 0;//1e-5; //-1; //-constants::pi;
        T xhi = 25;  //2*constants::pi;
        T xtest = 1; //0.5*(xhi+xlo);
        auto knots = knotsT(nknots, xlo, xhi); print(knotsT::name, knots.knots());
        auto xsam = oversample_knots(knots.knots()); print("xsam", xsam);
        long nbasis = nknots + order - 2;
        size_t nsam = xsam.dims()[0];
        
        BsplineBasis<T,knotsT> B(order, knots);
        
        for (size_t i=0; i<nsam; i++) {
            print("interval:", xsam[i], knots.interval(xsam[i]));
        }

        // for (long i=0; i<nbasis; i++) {
        //     auto b = B.bspline(xtest);
        //     Tensor<T> c(nbasis); c(i) = 1.0;
        //     print("deBoor", B.deBoor(xtest, c) - b[i]);
        // }
        
        //print(B.tabulate_basis(knots));
        //print(B.tabulate_basis(xsam));

        Tensor<T> f(nsam), df(nsam), d2f(nsam);
        for (size_t i=0; i<nsam; i++) {
            f[i] = std::exp(-2*xsam[i]);
            df[i] = -2*f[i];
            d2f[i] = 4*f[i];
        }

        //for (size_t i=0; i<nsam; i++) {
        //f[i] = std::sin(xsam[i]);
        //df[i] = std::cos(xsam[i]);
        //d2f[i] = -std::sin(xsam[i]);
        //}
        
        // {
        //     for (size_t nzeroL=0; nzeroL<=2; nzeroL++) {
        //         for (size_t nzeroR=0; nzeroR<=2; nzeroR++) {
        //             auto M = B.make_lsq_matrix(xsam,nzeroL,nzeroR);
        //             auto c1 = inner(M,f);
        //             print("sin c", nzeroL, nzeroR, c1);
        //         }
        //     }
        // }
        auto c = B.fit(xsam,f,0,2);
        //print(std::sin(xtest), B.deBoor(xtest, c));
        print(std::exp(-2*xtest), B.deBoor(xtest, c));
        print("|fit-f|", (B.deBoor(xsam, c) - f).normf());

        // Compute first derivative in the lower order basis
        BsplineBasis<T,knotsT> B1(order - 1, knots);
        auto d = B.deriv_exact(c);
        print(-2*std::exp(-2*xtest), B1.deBoor(xtest, d));
        //print(std::cos(xtest), B1.deBoor(xtest, d));

        auto dMat = B.make_deriv_exact_matrix();
        auto d2 = inner(dMat,c);
        auto df2 = B1.deBoor(xsam, d2);
        Tensor<T> derr1,derr2,derr3;
        Tensor<T> d2err1,d2err2,d2err3;
        derr1 = df2-df;
        print("Err in 'exact' first derivative", (df2-df).normf(), derr1[0]);
        //print(derr1);

        // Compute second derivative in the two lower order basis by applying deriv_exact twice
        BsplineBasis<T,knotsT> B2(order - 2, knots);
        {
            auto dMat1 = B1.make_deriv_exact_matrix();
            auto d2 = inner(dMat1,inner(dMat,c));
            auto df2 = B2.deBoor(xsam, d2);
            d2err1 = df2-d2f;
            print("Err in 'exact' second derivative", (df2-d2f).normf(), d2err1[0]);
        }
        
        // Compute first derivative in the same order basis
        auto dMat2 = B.make_deriv_matrix(xsam);
        auto d3 = inner(dMat2,c);
        auto df3 = B.deBoor(xsam, d3);
        derr2 = df3-df;
        //print(derr2);
        print("Err in fit+1 first derivative", derr2.normf(), derr2[0]);
        {
            auto dMat2 = B.make_deriv_matrixX(xsam);
            auto d3 = inner(dMat2,c);
            auto df3 = B.deBoor(xsam, d3);
            derr3 = df3-df;
            //print(derr3);
            print("Err in fit-1 first derivative", derr3.normf(), derr3[0]);
        }
        {
            char buf[256];
            snprintf(buf, sizeof(buf), "set style data l; set xrange [%.6f:%.6f]; set logscale x; set title '1st derivatives'", xlo, xhi);
            Gnuplot g(buf); // ; set xrange 
            g.plot(xsam, derr1, derr2, derr3);
        }
        // Compute second derivative in the same order basis
        {
            auto d2Mat = B.make_deriv2_matrix(xsam);
            d2 = inner(d2Mat,c);
            df2 = B.deBoor(xsam, d2);
            d2err2 = df2-d2f;
            print("Err in fit+2 second derivative", (df2-d2f).normf(), d2err2[0]);
        }
        {
            auto d2Mat = B.make_deriv2_matrixX(xsam);
            d2 = inner(d2Mat,c);
            df2 = B.deBoor(xsam, d2);
            d2err3 = df2-d2f;
            print("Err in fit-2 second derivative", (df2-d2f).normf(), d2err3[0]);
        }
        {
            char buf[256];
            snprintf(buf, sizeof(buf), "set style data l; set xrange [%.6f:%.6f]; set logscale x; set title '2nd derivatives'", xlo, xhi);
            Gnuplot g(buf); // ; set xrange 
            g.plot(xsam, d2err1, d2err2, d2err3);
        }

        return true;
    }
};

// Manages all of the data, including multiple matrices, that are
// needed to compute with the b-spline basis, including the basis
// itself.  It is presently assumed that the basis is the same for all
// functions and we are computing on a radial grid, so we include
// factors of r**2 and 4pi as needed.

// T is the same type as used in the BsplineBasis and Knots
template <typename T>
class BsplineData {
    static_assert(std::is_floating_point<T>::value, "BsplineData: T must be floating point");
public:
    typedef T value_type;
    //typedef KnotsChebyshev<scalar_type> knotsT;
    typedef KnotsUniform<T> knotsT;
    typedef BsplineBasis<T,knotsT> basisT;
    typedef Tensor<T> tensorT;

private:
    static const BsplineData<T>* data; // pointer to the singleton global data

    const basisT B; // the basis
    const tensorT rsam; // sample points
    const tensorT M[3][3]; // the LSQ matrix M(nbasis,nsamples), the pseudoinverse of A index by [nzeroL][nzeroR]
    const tensorT A; // the tabluated basis functions A(nsamples,nbasis)
    const std::pair<tensorT,tensorT> XW; // the quadrature points and weights for matrix elements
    const tensorT D; // first derivative operator
    const tensorT D2; // second derivative operator

    BsplineData(size_t order, size_t nknots, T rlo, T rhi)
        : B(order, knotsT(nknots, rlo, rhi))
        , rsam(oversample_knots(B.knots))
        , M {{B.make_lsq_matrix(rsam, 0, 0), B.make_lsq_matrix(rsam, 0, 1), B.make_lsq_matrix(rsam, 0, 2)},
             {B.make_lsq_matrix(rsam, 1, 0), B.make_lsq_matrix(rsam, 1, 1), B.make_lsq_matrix(rsam, 1, 2)},
             {B.make_lsq_matrix(rsam, 2, 0), B.make_lsq_matrix(rsam, 2, 1), B.make_lsq_matrix(rsam, 2, 2)}}
        , A(B.tabulate_basis(rsam))
        , XW(B.make_spline_quadrature((3*B.p+4)/2))
        , D(B.make_deriv_matrixX(rsam))
        , D2(B.make_deriv2_matrixX(rsam))
    { }
    
public:
    static const BsplineData<T>* ptr() {
        if (!data) throw "you forgot to call BsplineData::init()";
        return data;
    }

    static void init(size_t order, size_t nknots, T rlo, T rhi) {
        if (BsplineData<T>::data) throw "BsplineData already initialized";
        BsplineData<T>::data = new BsplineData<T>(order, nknots, rlo, rhi);
    }

    static void clean() {
        if (BsplineData<T>::data) delete BsplineData<T>::data;
        BsplineData<T>::data = nullptr;
    }

    static const size_t nbasis() { return basis().nbasis; }

    static const size_t order() { return basis().order; }

    static const size_t nknots() { return basis().nknots; }

    static const T rlo() { return knots().xlo; }

    static const T rhi() { return knots().xhi; }

    static const basisT& basis() { return ptr()->B; }

    static const knotsT& knots() { return basis().knots; }

    static const knotsT& knots_object() { return *(knotsT*)(&basis()); }

    static const tensorT& lsq_matrix(size_t nzeroL=0, size_t nzeroR=0) {
        MADNESS_ASSERT(nzeroL<3 and nzeroR<3);
        return ptr()->M[nzeroL][nzeroR];
    }

    static const tensorT& basis_matrix() { return ptr()->A; }

    static const tensorT& deriv_matrix() { return ptr()->D; }

    static const tensorT& deriv2_matrix() { return ptr()->D2; }

    static const tensorT& rsample() { return ptr()->rsam; }

    static const std::pair<tensorT,tensorT> quadrature() { return ptr()->XW; }

    // Tabulate the function at the sample points
    template <typename funcT>
    static auto tabulate(const funcT& func) {
        using resultT = decltype(func(T(0)));
        const auto& rsam = rsample();
        Tensor<resultT> f(rsam.size());
        for (size_t i=0; i<rsam.size(); ++i) f[i] = func(rsam[i]);
        return f;
    }

    template <typename funcT>
    static auto project(const funcT& func, size_t nzeroL=0, size_t nzeroR=0) {
        return inner(lsq_matrix(nzeroL,nzeroR), tabulate(func));
    }    
};

template <typename T> const BsplineData<T>* BsplineData<T>::data = nullptr;

template <typename T>
class BsplineFunction {
public:
    typedef T value_type;
    typedef typename TensorTypeData<T>::scalar_type scalar_type;

private:
    typedef typename TensorTypeData<T>::scalar_type scalarT;
    typedef Tensor<T> tensorT;
    typedef Tensor<scalarT> stensorT;
    typedef BsplineFunction<T> bfunctionT;
    typedef BsplineData<T> bdataT;

    template <typename U> friend class BsplineFunction;

    Tensor<T> c;

public:
    BsplineFunction(const bfunctionT& other) // Deep copy constructor
        : c(copy(other.c))
    {}

    BsplineFunction(bfunctionT&& other) // Move constructor
        : c(std::move(other.c))
    {}

    BsplineFunction(const tensorT& c = tensorT()) // Copies the coefficients, default is zero
        : c(c.size()>0 ? copy(c) : tensorT(bdataT::nbasis()))
    {}

    BsplineFunction(tensorT&& c) // Moves the coefficients
        : c(std::move(c))
    {}

    // Projects function (taking scalar_type as argument) onto the
    // basis.  If a functor or a lambda is provided instead of a bare
    // function pointer or an std::function then it could be inlined.
    template <typename funcT>
    BsplineFunction (const funcT& func, size_t nzeroR=0, size_t nzeroL=0)
        : c(bdataT::project(func, nzeroR, nzeroL))
    {}

    BsplineFunction& operator=(const bfunctionT& other) { // Deep assignment
        if (this != &other) {
            c = copy(other.c);
        }
        return *this;
    }
    
    auto operator+(const bfunctionT& other) const {
        return bfunctionT(c + other.c);
    }

    bfunctionT operator+(const T& s) const {
        return bfunctionT(c + s);
    }

    void operator+=(const bfunctionT& other) {
        c += other.c;
    }

    void operator+=(const T& s) {
        c += s;
    }

    bfunctionT operator-(const bfunctionT& other) const {
        return BsplineFunction(c - other.c);
    }

    bfunctionT operator-(const T& s) const {
        return BsplineFunction(c - s);
    }

    void operator-=(const bfunctionT& other) {
        c -= other.c;
    }

    void operator-=(const T& s) {
        c -= s;
    }

    bfunctionT operator*(const bfunctionT& other) const {
        const stensorT& A = bdataT::basis_matrix();
        const stensorT& M = bdataT::lsq_matrix();
        auto f = inner(A,c);
        auto g = inner(A,other.c);
        auto h = inner(M,f.emul(g));
        return bfunctionT(h);
    }

    bfunctionT operator*(const T& s) const {
        return BsplineFunction(c*s);
    }

    void operator*=(const T& s) {
        c *= s;
    }

    // Evaulate the function a point x
    T operator()(scalarT x) const {
        return bdataT::basis().deBoor(x, c);
    }

    // Evaulate the function at multiple points
    tensorT operator()(const stensorT& x) const {
        return bdataT::basis().deBoor(x, c);
    }

    // Differentiate the function (order=1 or 2)
    bfunctionT D(size_t order) const {
        if (order==1) {
            return bfunctionT(inner(bdataT::deriv_matrix(), c));
        } else if (order==2) {
            return bfunctionT(inner(bdataT::deriv2_matrix(), c));
        } else {
            throw std::runtime_error("BsplineFunction::D: order must be 1 or 2");
        }
        return bfunctionT(inner(bdataT::deriv_matrix(), c));
    }

    static void test() {
        size_t nknots = 61;
        size_t order = 8;
        T xlo = 0.0;
        T xhi = 50.0;
        bdataT::init(order, nknots, xlo, xhi);

        auto f = [](scalarT x){ return std::exp(-x); };
        auto fdat = bdataT::tabulate(f);
        
        bfunctionT g(f);
        print(g(0.5), std::exp(-0.5));
        // {
        //     Gnuplot G("set style data lines; set title 'error in g=exp(-x)'; set xlabel 'x'; set ylabel 'g(x)-exp(-x)'; set key left top");
        //     G.plot(bdataT::rsample(), g(bdataT::rsample())-fdat);
        // }

        auto h = g*g;
        print(h(0.5), std::exp(-2*0.5));
        // {
        //     Gnuplot G("set style data lines; set title 'error in g**2'; set xlabel 'x'; set ylabel 'g(x)**2-exp(-2x)'; set key left top");
        //     G.plot(bdataT::rsample(), h(bdataT::rsample())-fdat.emul(fdat));
        // }

        g += 1;

        // Construct the overlap, kinetic, and potential matrices and diagonalize to verify that they are correct.
        // This is brute force using none of the massive sparsity.
        stensorT S, PE, KE, H;
        {
            auto [X, W] = bdataT::quadrature();
            //print("X");
            //print(X);
            auto nbasis = bdataT::nbasis();
            print(nbasis);
            auto npt = X.size();
            print(npt);
            auto b = bdataT::basis().tabulate_basis(X); // basis functions evaulated at the GL points
            // {
            //     Gnuplot G("set style data lines; set title 'basis functions'; set xlabel 'x'; set ylabel 'b(x)'; set key right top; set xrange [0:2]");
            //     G.plot(X, b(_,0), b(_,1), b(_,2), b(_,3), b(_,4), b(_,5), b(_,6), b(_,7));
            // }
            auto wb = copy(b); // with weights and r**2 applied
            for (size_t i=0; i<npt; i++) {
                T wi = W[i];
                T ri = X[i];
                for (size_t j=0; j<nbasis; j++) {
                    wb(i,j) *= wi*ri;
                }
            }
            PE = -inner(b, wb, 0, 0);
            //print("potential");
            //print(PE);
            for (size_t i=0; i<npt; i++) {
                T ri = X[i];
                for (size_t j=0; j<nbasis; j++) {
                    wb(i,j) *= ri;
                }
            }
            S = inner(b, wb, 0, 0);
            //print("overlap");
            //print(S);

            stensorT L;
            {
                BsplineBasis<T,typename bdataT::knotsT> B0(order    , bdataT::knots_object());
                BsplineBasis<T,typename bdataT::knotsT> B1(order - 1, bdataT::knots_object());
                BsplineBasis<T,typename bdataT::knotsT> B2(order - 2, bdataT::knots_object());
                Tensor<T> dMat1 = B0.make_deriv_exact_matrix();
                Tensor<T> dMat2 = B1.make_deriv_exact_matrix();
                stensorT A = B2.tabulate_basis(X);
                L = inner(inner(A, dMat2), dMat1);
                //print("L");
                //print(L);
            }
            //Tensor<T> L = inner(b, bdataT::deriv2_matrix(), 1, 0);
            for (size_t i=0; i<npt; i++) {
                T ri = X[i];
                for (size_t j=0; j<nbasis; j++) {
                    L(i,j) *= ri;
                }
            }
            stensorT D1;
            {
                BsplineBasis<T,typename bdataT::knotsT> B0(order    , bdataT::knots_object());
                BsplineBasis<T,typename bdataT::knotsT> B1(order - 1, bdataT::knots_object());
                Tensor<T> dMat1 = B0.make_deriv_exact_matrix();
                stensorT A = B1.tabulate_basis(X);
                D1 = inner(A, dMat1);
            }
            L += T(2)*D1;
            //L += inner(b, T(2)*bdataT::deriv_matrix(), 1, 0);
            for (size_t i=0; i<npt; i++) {
                T wi = W[i];
                T ri = X[i];
                for (size_t j=0; j<nbasis; j++) {
                    L(i,j) *= wi*ri;
                }
            }

            printf("%.16e\n", b(0,0));
            KE = T(-0.5)*inner(b, L, 0, 0);
            //print("kinetic");
            //print(KE);
            //print(KE-transpose(KE));
        }
        {
            Tensor<T> V, e;
            syev(S, V, e);
            print("eigenvalues of S", e);
        }
        {
            auto nbasis = bdataT::nbasis();
            Tensor<std::complex<T>> V(nbasis, nbasis);
            Tensor<std::complex<T>> e(nbasis);
            ggev(KE, S, V, e);
            print("eigenvalues of KE", e);
        }
        {
            auto nbasis = bdataT::nbasis();
            Tensor<std::complex<T>> V(nbasis, nbasis);
            Tensor<std::complex<T>> e(nbasis);
            stensorT H = KE+PE;
            //print("H");
            //print(H);
            ggev(H, S, V, e);
            print("eigenvalues of H", e);
            {
                stensorT e(nbasis-1);
                stensorT HH(H(Slice(0,-2), Slice(0,-2)));
                stensorT SS(S(Slice(0,-2), Slice(0,-2)));
                stensorT V;
                // print("HH");
                // print(HH);
                // print("SS");
                // print(SS);
                sygv(HH, SS, 1, V, e);
                print("eigenvalues of H", e);
            }
        }

        

        // Solve H atom
        // 1) Project guess
        // 2) Project potential
        // 3) Compute energy = <psi|T+V|psi>/<psi|psi>
        // 4) V * psi
        // 5) -2*G V*psi
        // 6) Compute residual norm
    }
};


// // Make the matrix that applies the second derivative after projecting into a a two higher order basis.
// // Gives us 1 order higher accuracy and produces a result in the same order basis as the input.
// // make_deriv2_matrix := proc(t, order, Xsample) local Npknots, Nbasis, t1, t2, A, M, dMat1, dMat2; Npknots := numelems(t); Nbasis := Npknots - order; t1 := Vector(Npknots + 2); t1[2 .. Npknots + 1] := t; t1[1] := t[1]; t1[Npknots + 2] := t[Npknots]; t2 := Vector(Npknots + 4); t2[2 .. Npknots + 3] := t1; t2[1] := t[1]; t2[Npknots + 4] := t[Npknots]; A := tabulate_basis(t, order, Xsample); M := make_lsq_matrix(t2, order + 2, Xsample); dMat2 := make_deriv_exact_matrix(t2, order + 2); dMat1 := make_deriv_exact_matrix(t1, order + 1); return ((dMat1 . dMat2) . M) . A; end proc;
// Tensor<T> make_deriv2_matrix(const Tensor<T>& t, size_t order, const Tensor<T>& Xsample) {
//     size_t Npknots = t.dims()[0];
//     size_t Nbasis = Npknots - order;
//     Tensor<T> t1(Npknots + 2);
//     t1[0] = t[0];
//     t1[1] = t[0];
//     for (size_t i=0; i<Npknots; i++) {
//         t1[i + 2] = t[i];
//     }
//     t1[Npknots + 1] = t[Npknots - 1];
//     Tensor<T> t2(Npknots + 4);
//     t2[0] = t[0];
//     t2[1] = t[0];
//     for (size_t i=0; i<Npknots + 2; i++) {
//         t2[i + 2] = t1[i];
//     }
//     t2[Npknots + 3] = t[Npknots - 1];
//     Tensor<T> A = tabulate_basis(t, order, Xsample);
//     Tensor<T> M = make_lsq_matrix(t2, order + 2, Xsample);
//     Tensor<T> dMat2 = make_deriv_exact_matrix(t2, order + 2);
//     Tensor<T> dMat1 = make_deriv_exact_matrix(t1, order + 1);
//     return ((dMat1.dot(dMat2)).dot(M)).dot(A);
// }

// // Perform quadrature on a function over [a,b] given a rule over [-1,1]
// // quadrature := (f, a, b, X, W) -> local i; evalf(1/2*(b - a)*add(f(1/2*a + 1/2*b + 1/2*X[i]*(b - a))*W[i], i = 1 .. numelems(X)));
// T quadrature(const std::function<T(T)>& f, T a, T b, const Tensor<T>& X, const Tensor<T>& W) {
//     T sum = 0.0;
//     for (size_t i=0; i<X.dims()[0]; i++) {
//         sum += f(0.5*a + 0.5*b + 0.5*X[i]*(b - a))*W[i];
//     }
//     return 0.5*(b - a)*sum;
// }

// Given quadrature points and weights on [-1,1] adaptively compute int(f(x),x=a..b) to accuracy eps
// adaptive_quadrature := proc(f, a, b, X, W, eps) local i0, i1; i0 := quadrature(f, a, b, X, W); i1 := quadrature(f, a, 1/2*a + 1/2*b, X, W) + quadrature(f, 1/2*a + 1/2*b, b, X, W); print(a, b, i0, i1, abs(i0 - i1)); if abs(i0 - i1) <= eps then return i1; else return adaptive_quadrature(f, a, 1/2*a + 1/2*b, X, W, 0.5*eps) + adaptive_quadrature(f, 1/2*a + 1/2*b, b, X, W, 0.5*eps); end if; end proc;


// make_overlap_matrix := proc(order, Knots, pKnots) local X, W, A, B, i; X, W := make_spline_quadrature(order + 1, Knots); A := tabulate_basis(pKnots, order, X); B := copy(Transpose(A)); for i to numelems(X) do A[i] := A[i]*W[i]; end do; return B . A; end proc;
// Tensor<T> make_overlap_matrix(size_t order, const Tensor<T>& Knots, const Tensor<T>& pKnots) {
//     Tensor<T> X, W;
//     make_spline_quadrature(order + 1, Knots, X, W);
//     Tensor<T> A = tabulate_basis(pKnots, order, X);
//     Tensor<T> B = A.transpose();
//     for (size_t i=0; i<X.dims()[0]; i++) {
//         A[i] *= W[i];
//     }
//     return B.dot(A);
// }

int main() {
    //BsplineBasis<float>::test();
    //BsplineBasis<double,KnotsUniform<double>>::test();
    //BsplineBasis<double,KnotsChebyshev<double>>::test();
    //BsplineBasis<double,KnotsGeometricShifted<double>>::test();
    BsplineFunction<double>::test();
    return 0;
}
