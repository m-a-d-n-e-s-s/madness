#include <limits>
#include <memory>
#include <complex>
#include <algorithm>
#include <map>

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


// Define a gnuplot data block from a 1-d or 2-d tensor
template <typename T>
void gnuplot_tensor_db(Gnuplot& g, const std::string& name, const Tensor<T>& t) {
    if (t.ndim() > 2) throw "gnuplot: too many dimensions";
    g("$",false);
    g(name,false);
    g(" << EOD");
    char buf[256];
    const size_t n=t.dims()[0], m=t.dims()[1];
    for (size_t i=0; i<n; ++i) {
        for (size_t j=0; j<m; ++j) {
            snprintf(buf,sizeof(buf),"%16.8e ",double(t(i,j)));
            g(buf,j==(m-1));
        }
    }
    g("EOD");
}
        
// Plot multiple columns of a tensor
template <typename T>
void gnuplot_tensor(Gnuplot& g, const Tensor<T>& x, const Tensor<T>& t) {
    const size_t n=t.dims()[0], m=t.dims()[1];
    {
        Tensor<T> data(n,m+1);
        data(_,0) = x;
        data(_,Slice(1,-1)) = t;
        gnuplot_tensor_db(g, "data", data);
    }
    char buf[256];
    std::string cmd = "plot ";
    for (size_t j=0; j<m; j++) {
        snprintf(buf,sizeof(buf),"$data using 1:%zu%s",j+2,(j==m-1)? "" : ", ");
        cmd += buf;
    }
    g(cmd);
}

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

// Quartic to linear rational form
template <typename T>
class KnotsRational {
public:
    static constexpr const char* name = "rational form";
    typedef T value_type;

private:
    const size_t nknots; // number of knots
    const T a;   // scaling parameter
    const T xhi; // [0,xhi] interval

    T Q(T i) const {
        return i*i*(T(1) + a)/(a*i + T(1));
    }

    T Qinv(T x) const {
        return (a*x + std::sqrt(a*a*x*x + T(4)*a*x + T(4)*x))/(T(2)*(T(1) + a));
    }
            

public:
    KnotsRational(size_t nknots, T a, T xhi)
        : nknots(nknots)
        , a(a)
        , xhi(xhi)
    {}

    Tensor<T> knots() const {
        Tensor<T> pts(nknots);
        for (size_t i=0; i<nknots-1; i++) {
            T ii = T(i)/T(nknots-1);
            pts[i] = xhi*Q(Q(ii));
        }
        pts[nknots-1] = xhi; // to avoid rounding errors
        return pts;
    }

    size_t interval(T x) const {
#if !defined(NDEBUG)
        if (x < 0 || x > xhi) throw "Xinterval: x not in t";
#endif
        T xx = (x/xhi) * (T(1)-2*std::numeric_limits<T>::epsilon());
        T ii = Qinv(Qinv(xx));
        size_t j = size_t(ii*(nknots-1));

#if !defined(NDEBUG)
        if (j >= (nknots-1)) throw "interval: result interval out of range";
#endif
        return j;
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
        , rlogh((1.0-10*std::numeric_limits<T>::epsilon())/std::log(h)) // to avoid endpoint issues ...
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
    const T rLhalf; // 2.0/(xhi-xlo)

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
            T sgn = (knots[i] > 0) ? 1 : -1;
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

    // Make the matrix that applies the derivative after projecting
    // into a one higher order basis.  Gives us 1 order higher
    // accuracy and produces a result in the same order basis as the
    // input.  nzeroL/R are the number of basis functions to be
    // neglected on the left/right side of the interval when fitting
    // the *input* function in the higher-order basis.
    Tensor<T> make_deriv_matrix(const Tensor<T>& Xsample, size_t nzeroL, size_t nzeroR) const {
        size_t nsample = Xsample.dims()[0];
        BsplineBasis<T,knotsT> B1 = BsplineBasis(order + 1, knots_object());
        Tensor<T> A = tabulate_basis(Xsample);
        Tensor<T> M = B1.make_lsq_matrix(Xsample, nzeroL, nzeroR);
        Tensor<T> dMat = B1.make_deriv_exact_matrix();
        return inner(inner(dMat,M),A);
    }

    // Make the matrix that applies the exact derivative and projects back into the original order basis.

    // nzL/nzR reflect the input function, so the fitting is done in a one-lower order basis.
    Tensor<T> make_deriv_matrixX(const Tensor<T>& Xsample, size_t nzeroL, size_t nzeroR) const {
        size_t nsample = Xsample.dims()[0];
        BsplineBasis<T,knotsT> B1 = BsplineBasis(order - 1, knots_object());
        Tensor<T> A = B1.tabulate_basis(Xsample);
        size_t nzL = std::max(nzeroL, 1ul)-1;
        size_t nzR = std::max(nzeroR, 1ul)-1;
        Tensor<T> M = make_lsq_matrix(Xsample,nzL,nzR);
        Tensor<T> dMat = make_deriv_exact_matrix();
        return inner(inner(M,A),dMat);
    }

    // Make the matrix that applies the second derivative after
    // projecting into a two higher order basis.  Gives us 2 order
    // higher accuracy and produces a result in the same order basis
    // as the input.  nzeroL/R are the number of basis functions to be
    // neglected on the left/right side of the interval when fitting
    // the *input* function in the higher-order basis.
    Tensor<T> make_deriv2_matrix(const Tensor<T>& Xsample, size_t nzeroL, size_t nzeroR) const {
        size_t nsample = Xsample.dims()[0];
        BsplineBasis<T,knotsT> B1 = BsplineBasis(order + 1, knots_object());
        BsplineBasis<T,knotsT> B2 = BsplineBasis(order + 2, knots_object());
        Tensor<T> A = tabulate_basis(Xsample);
        Tensor<T> M2 = B2.make_lsq_matrix(Xsample,nzeroL,nzeroR);
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
        size_t nknots = 15;
        size_t order = 6; 
        T xlo = 0;//1e-5; //-1; //-constants::pi;
        T xhi = 25;  //2*constants::pi;
        T xtest = 1; //0.5*(xhi+xlo);
        //auto knots = knotsT(nknots, xlo, xhi); print(knotsT::name, knots.knots());
        auto knots = knotsT(nknots, 0.1, xhi); print(knotsT::name, knots.knots());
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
        {
            auto xxsam = oversample_knots(oversample_knots(oversample_knots(oversample_knots(xsam))));
            auto A = B.tabulate_basis(xxsam);
            Gnuplot G("set style data lines; set title 'B-spline basis'; set xlabel 'x'; set ylabel 'b[i](x)'; set key off","basis.gnuplot");
            gnuplot_tensor(G, xxsam, A);
        }

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
        auto dMat2 = B.make_deriv_matrix(xsam,0,0);
        auto d3 = inner(dMat2,c);
        auto df3 = B.deBoor(xsam, d3);
        derr2 = df3-df;
        //print(derr2);
        print("Err in fit+1 first derivative", derr2.normf(), derr2[0]);
        {
            auto dMat2 = B.make_deriv_matrixX(xsam,0,0);
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
            auto d2Mat = B.make_deriv2_matrix(xsam,0,0);
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
    typedef KnotsRational<T> knotsT;
    //typedef KnotsChebyshev<T> knotsT;
    //typedef KnotsGeometricShifted<T> knotsT;
    //typedef KnotsUniform<T> knotsT;
    typedef BsplineBasis<T,knotsT> basisT;
    typedef Tensor<T> tensorT;

    //private:
    static const BsplineData<T>* data; // pointer to the singleton global data

    const basisT B;       // b-spline basis
    const tensorT ones;   // ones[0..nbasis-1] = 1 to make constructor code easier to read
    const tensorT rsam;   // sample points
    const std::pair<tensorT,tensorT> XW; // the quadrature points and weights for matrix elements
    const tensorT Asam;   // bf at the sample points A(mu,j)=b[j](x[mu])
    const tensorT A;      // bf at the GL quadrature points A(mu,j)=b[j](X[mu])
    const tensorT Ar;     // bf*r at the GL quadrature points A(mu,j)=b[j](X[mu])*X[mu]
    const tensorT AW;     // weighted bf at the GL quadrature points AW(mu,j) = b[j](X[mu])*w[mu]
    const tensorT AWr;    // weighted*r bf at the GL quadrature points AWr(mu,j) = b[j](X[mu])*w[mu]*X[mu]
    const tensorT AWr2;   // weighted*r*r bf at the GL quadrature points AWr(mu,j) = b[j](X[mu])*w[mu]*X[mu]**2
    const tensorT btrace; // btrace[i]=4*PI*int(r**2 b[i](r) dr, r=0..rmax)
    const tensorT S;      // overlap matrix 4*PI*int(r**2 b[i](r) b[j](r) dr, r=0..rmax)
    mutable tensorT KE;   // kinetic energy matrix (+1/2)*4*PI*int(d/dr(r b[i](r)) d/dr(r b[j](r)) dr, r=0..rmax) assuming zero RHS BC
    mutable std::map<std::pair<size_t,size_t>,tensorT> M; // LSQ matrix (pseudoinverse of A) indexed by [nzeroL][nzeroR]
    mutable std::map<std::pair<size_t,size_t>,tensorT> D; // first derivative matrix indexed by [nzeroL][nzeroR]
    mutable std::map<std::pair<size_t,size_t>,tensorT> D2;// second derivative matrix indexed by [nzeroL][nzeroR]

    BsplineData(size_t order, size_t nknots, T rlo, T rhi)
        : B(order, knotsT(nknots, rlo, rhi))
        , ones(Tensor<T>(B.nbasis,1).fill(1))
        , rsam(oversample_knots(B.knots))
        , XW(B.make_spline_quadrature((2*B.p+4)/2))
        , Asam(B.tabulate_basis(rsam))
        , A(B.tabulate_basis(XW.first))
        , Ar(copy(A).emul(outer(XW.first,ones)))
        , AW(copy(A).emul(outer(XW.second,ones)))
        , AWr(copy(AW).emul(outer(XW.first,ones)))
        , AWr2(copy(AWr).emul(outer(XW.first,ones)))
        , btrace(4*constants::pi*inner(AWr,XW.first,0,0))
          //, S(inner(A,AWr2,0,0)*(4*constants::pi))
        , S(inner(Ar,AWr,0,0)*(4*constants::pi))
        , KE()
        , M()
        , D()
        , D2()
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

    static const Tensor<T>& knots() { return basis().knots; }

    static const knotsT& knots_object() { return *(const knotsT*)(&basis()); }

    static const tensorT& lsq_matrix(size_t nzeroL=0, size_t nzeroR=0) {
        MADNESS_ASSERT(nzeroL<=order() and nzeroR<=order());
        if (!ptr()->M.count(std::make_pair(nzeroL,nzeroR))) {
            ptr()->M[std::make_pair(nzeroL,nzeroR)] = basis().make_lsq_matrix(rsample(),nzeroL,nzeroR);
        }
        return ptr()->M[std::make_pair(nzeroL,nzeroR)];
    }

    static const tensorT& basis_at_sample_points() { return ptr()->Asam; }

    static const tensorT& basis_at_GL_points() { return ptr()->A; }

    static const tensorT& overlap_matrix() { return ptr()->S; }

    static const tensorT& ke_matrix() {
        if (ptr()->KE.size() == 0) {
            // // d/dr(r f) = f + r df/dr = (A + Ar D) c
            auto Dexact = ptr()->basis().make_deriv_exact_matrix();
            auto B1 = BsplineBasis<T,knotsT>(ptr()->order()-1,ptr()->knots_object());
            auto A1 = B1.tabulate_basis(ptr()->quadrature().first);
            auto F = ptr()->A + inner(A1,Dexact).emul(outer(ptr()->quadrature().first,ptr()->ones));
            auto FW= copy(F).emul(outer(ptr()->quadrature().second,ptr()->ones));
            ptr()->KE = (0.5*4*constants::pi)*inner(FW,F,0,0);
        }
        return ptr()->KE;
    }

    static const tensorT& deriv_matrix(size_t nzeroL=0, size_t nzeroR=0) {
        if (!ptr()->D.count(std::make_pair(nzeroL,nzeroR))) {
            ptr()->D[std::make_pair(nzeroL,nzeroR)] = basis().make_deriv_matrixX(rsample(),nzeroL,nzeroR); // was X
        }
        return ptr()->D[std::make_pair(nzeroL,nzeroR)];
    }

    static const tensorT& deriv2_matrix(size_t nzeroL=0, size_t nzeroR=0) {
        if (!ptr()->D2.count(std::make_pair(nzeroL,nzeroR))) {
            ptr()->D2[std::make_pair(nzeroL,nzeroR)] = basis().make_deriv2_matrix(rsample(),nzeroL,nzeroR);
        }
        return ptr()->D2[std::make_pair(nzeroL,nzeroR)];
    }

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

double exact_energy(double Z) {
    const double c = 137.035999679;
    double Za = Z/c;
    double s = Za*Za / (1.0 - Za*Za);
    return c*c/std::sqrt(1.0 + s) - c*c;
}

// Here T (real or complex) is the type of the function/coefficients and scalarT (real) is the type of the basis functions and knots
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
    typedef BsplineData<scalarT> bdataT;

    template <typename U> friend class BsplineFunction;

    size_t nzeroL=0, nzeroR=0; // number of zero values/derivatives on the left and right
    Tensor<T> c; // The b-spline coefficients

public:
    BsplineFunction(const bfunctionT& other) // Deep copy constructor
        : nzeroL(other.nzeroL)
        , nzeroR(other.nzeroR)
        , c(copy(other.c))
    {print("in copy fun constructor");}

    BsplineFunction(bfunctionT&& other) // Move constructor
        : nzeroL(other.nzeroL)
        , nzeroR(other.nzeroR)
        , c(std::move(other.c))
    {print("in move fun constructor");}

    BsplineFunction(const tensorT& c = tensorT(), size_t nzeroL=0, size_t nzeroR=0) // Copies the coefficients, default is zero
        : nzeroL(nzeroL)
        , nzeroR(nzeroR)
        , c(c.size()>0 ? copy(c) : tensorT(bdataT::nbasis()))
    {print("in copy ten constructor");}

    BsplineFunction(tensorT&& c, size_t nzeroL=0, size_t nzeroR=0) // Moves the coefficients
        : nzeroL(nzeroL)
        , nzeroR(nzeroR)
        , c(std::move(c))
    {print("in move ten constructor");}

    // Projects function (taking scalar_type as argument) onto the
    // basis.  If a functor or a lambda is provided instead of a bare
    // function pointer or an std::function then it could be inlined.
    template <typename funcT>
    BsplineFunction (const funcT& func, size_t nzeroR=0, size_t nzeroL=0)
        : nzeroL(nzeroL)
        , nzeroR(nzeroR)
        , c(bdataT::project(func, nzeroR, nzeroL)) 
    {}

    size_t nzL() const { return nzeroL; }

    size_t nzR() const { return nzeroR; }

    const tensorT& coeffs() const { return c; }

    const T& trace() const {
        
    }

    const scalarT norm() const { return std::sqrt(c.trace_conj(::madness::inner(bdataT::overlap_matrix(),c))); }

    template <typename U>
    TENSOR_RESULT_TYPE(T,U) inner(const BsplineFunction<U>& other) const {
        return c.trace_conj(::madness::inner(bdataT::overlap_matrix(),other.c));
    }

    BsplineFunction& operator=(const bfunctionT& other) { // Deep assignment
        if (this != &other) {
            nzeroL = other.nzeroL;
            nzeroR = other.nzeroR;
            c = copy(other.c);
        }
        return *this;
    }
    
    auto operator+(const bfunctionT& other) const {
        // x^n + x^m = O(x^min(n,m))
        size_t nzL = std::min(nzeroL, other.nzeroL);
        size_t nzR = std::min(nzeroR, other.nzeroR);
        return bfunctionT(c + other.c, nzL, nzR);
    }

    bfunctionT operator+(const T& s) const {
        return bfunctionT(c + s, 0, 0);
    }

    void operator+=(const bfunctionT& other) {
        nzeroL = std::min(this->nzeroL, other.nzeroL);
        nzeroR = std::min(this->nzeroR, other.nzeroR);
        c += other.c;
    }

    void operator+=(const T& s) {
        nzeroL=0;
        nzeroR=0;
        c += s;
    }

    bfunctionT operator-(const bfunctionT& other) const {
        size_t nzL = std::min(this->nzeroL, other.nzeroL);
        size_t nzR = std::min(this->nzeroR, other.nzeroR);
        return BsplineFunction(c - other.c, nzL, nzR);
    }

    bfunctionT operator-(const T& s) const {
        return BsplineFunction(c - s, 0, 0);
    }

    void operator-=(const bfunctionT& other) {
        nzeroL = std::min(this->nzeroL, other.nzeroL);
        nzeroR = std::min(this->nzeroR, other.nzeroR);
        c -= other.c;
    }

    void operator-=(const T& s) {
        nzeroL=0;
        nzeroR=0;
        c -= s;
    }

    bfunctionT operator*(const bfunctionT& other) const {
        // x^n * x^m = x^(n+m)
        size_t nzL = std::min(bdataT::order(), nzeroL + other.nzeroL);
        size_t nzR = std::min(bdataT::order(), nzeroR + other.nzeroR);
        const stensorT& A = bdataT::basis_at_sample_points();
        const stensorT& M = bdataT::lsq_matrix(nzL, nzR);
        auto f = ::madness::inner(A,c);
        auto g = ::madness::inner(A,other.c);
        auto h = ::madness::inner(M,f.emul(g));
        return bfunctionT(h, nzL, nzR);
    }

    bfunctionT operator*(const T& s) const {
        return BsplineFunction(c*s, nzeroL, nzeroR);
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

    // Differentiate the function (Dorder=1 or 2)
    bfunctionT D(size_t Dorder=1) const {
        if (Dorder==1) {
            size_t nzL = nzeroL > 1 ? nzeroL-1 : 0;
            size_t nzR = nzeroR > 1 ? nzeroR-1 : 0;
            return bfunctionT(::madness::inner(bdataT::deriv_matrix(nzeroL,nzeroR), c), nzL, nzR);
        } else if (Dorder==2) {
            size_t nzL = nzeroL > 2? nzeroL-2 : 0;
            size_t nzR = nzeroR > 2? nzeroR-2 : 0;
            return bfunctionT(::madness::inner(bdataT::deriv2_matrix(nzeroL,nzeroR), c), nzL, nzR);
        } else {
            throw std::runtime_error("BsplineFunction::D: order must be 1 or 2");
        }
    }

    // static void test_fit() {
    //     const T Z = 10.0;
    //     size_t nknots = 61;
    //     size_t order = 20;
    //     T xlo = 0.1; // a for KnotsRational //.01 for 103 //1e-7;
    //     T xhi = 26.0;
        
    //     bdataT::init(order, nknots, xlo, xhi);

    //     print("knots", bdataT::knots());

    //     auto A = bdataT::basis_at_GL_points(); 
    //     auto [X, W] = bdataT::quadrature();
    //     size_t ltest = 1;
    //     ?????????????????????????????
        

    // }

    

    static void test() {
        const T Z = 10.0;
        size_t nknots = 61;
        size_t order = 20;
        T xlo = 0.1; // a for KnotsRational //.01 for 103 //1e-7;
        T xhi = 26.0;
        
        bdataT::init(order, nknots, xlo, xhi);

        print("knots", bdataT::knots());
        auto f = [](scalarT x){ return std::exp(-x); };
        auto fdat = bdataT::tabulate(f);
        
        // bfunctionT g(f,0,2);
        // print(g(0.5), std::exp(-0.5));
        // {
        //     Gnuplot G("set style data lines; set title 'error in g=exp(-x)'; set xlabel 'x'; set ylabel 'g(x)-exp(-x)'; set key left top");
        //     G.plot(bdataT::rsample(), g(bdataT::rsample())-fdat);
        // }

        // auto h = g*g;
        // print(h(0.5), std::exp(-2*0.5));
        
        // {
        //     Gnuplot G("set style data lines; set title 'error in g**2'; set xlabel 'x'; set ylabel 'g(x)**2-exp(-2x)'; set key left top");
        //     G.plot(bdataT::rsample(), h(bdataT::rsample())-fdat.emul(fdat));
        // }

        // print("g(0.5)   = ", g(0.5), f(0.5));
        // g += 1;
        // print("g(0.5)+1 = ", g(0.5), f(0.5)+1);
        // print("Dg(0.5)", g.D()(0.5), -f(0.5));

        // Construct the overlap, kinetic, and potential matrices and diagonalize to verify that they are correct.
        stensorT S  = bdataT::overlap_matrix();
        stensorT KE = bdataT::ke_matrix();
        stensorT PE = (-Z*4*constants::pi)*::madness::inner(bdataT::ptr()->AWr,bdataT::ptr()->A,0,0);

        // Now balance the matrices to try to improve the condition numbers --- jacobi preconditioning is equivalent to normalizing the basis functions
        size_t nbasis = bdataT::nbasis();
        stensorT d(nbasis);
        ITERATOR(d, d[_i] = 1.0/std::sqrt(S(_i,_i)));
        //auto dd = outer(d,d);
        //S.emul(dd);
        //KE.emul(dd);
        //PE.emul(dd);
        {
            stensorT V, e;
            syev(S, V, e);
            print("eigenvalues of S", e);
        }
        {
            stensorT V, e;
            stensorT H = bdataT::ke_matrix();
            sygv(H, S, 1, V, e);
            print("eigenvalues of KE", e);
        }
        stensorT c0;
        T S00;
        {
            stensorT V, e;
            stensorT H = bdataT::ke_matrix()+PE;
            sygv(H, S, 1, V, e);
            print("eigenvalues of H", e);
            c0 = V(_,0);
            print("c0", c0);
            S00 = c0.trace(::madness::inner(S,c0));
            print("S00", S00);
            print("EEEE",c0.trace(::madness::inner(H,c0))/S00);

            {
                    bfunctionT psi(c0);
                    auto rsam = bdataT::ptr()->rsample();
                    Gnuplot g("set title \"non-rel\"; set style data lines; set xrange [0:1]");
                    g.plot(rsam, psi(rsam));
                }
            
        }

        // {
        //     auto f = [](T x) { return 1.0; };
        //     bfunctionT F(f);
        //     print("fitting a constant");
        //     print(F.coeffs());
        // }

        // Solve H atom
        // 1) Project guess
        // 2) Project potential
        // 3) Compute energy = <psi|T+V|psi>/<psi|psi>
        // 4) V * psi
        // 5) -2*G V*psi
        // 6) Compute residual norm

        bfunctionT r([](T x){return x;}, 1, 0);
        bfunctionT rFexact11([](T x){return x*std::exp(-x);}, 1, 1);
        bfunctionT rFexact00([](T x){return x*std::exp(-x);}, 0, 0);
        bfunctionT drFexact([](T x){return (1-x)*std::exp(-x);}, 0, 0);
        bfunctionT F(f, 0, 1);
        auto rF = r*F; print("rF", rF.nzL(), rF.nzR());
        auto drF = rF.D(); print("drF", drF.nzL(), drF.nzR());
        auto drFX = F + r*F.D(); print("drFX", drFX.nzL(), drFX.nzR());
        print("err", (drF-drFX).norm());

        // print(" rF analytic", 0.3*f(0.3), rF(0.3));
        // print("DrF analytic", f(0.3)-0.3*f(0.3), drF(0.3), drFX(0.3));

        // auto dftest(::madness::inner(bdataT::ptr()->A + ::madness::inner(bdataT::ptr()->Ar,bdataT::deriv_matrix(0,1)), F.coeffs()));
        // {
        //     Gnuplot gG("set style data line; set key left top; set xlabel 'x'; set ylabel 'y'; set logscale x");
        //     auto X = bdataT::quadrature().first;
        //     gG.plot(X, drF(X)- drFexact(X), drFX(X)- drFexact(X), dftest- drFexact(X));
        // }
        
        {
            // In here implement SF Dirac
            const T c = 137.035999679;

            // <i|p.Vp|j> = <pi|V|pj> = <di/dr | V | dj/dr> = Z int( 4 pi r^2 i' 1/r j' ) = Z 4 pi int(r i' j') = Z 4 pi sum((R A D)^T W (R A D)) where A is in order-1 basis

            //  d/dx V(r) d/dx f(r)
            //= d/dx V(r)  x/r f'(r)
            //=  x/r V' x/r f' +  V  1/r f'(r) +  V(r) -x^2/r^3 f'(r) +  V(r)  x^2/r^2 f''
            //=  x^2/r^2 V' f' +  V f' / r - V f' x^2/r^3 + V f'' x^2/r^2
            //sum(x) --> V' f' +  2 V f' / r + V f''
            // now compute spherical integral
            // 4 pi int(r^2 g (V' f' +  2 V f' / r + V f''), r=0..R)
            //=4 pi int(r^2 g V' f', r=0..R) + 8 pi int(r g V f', r=0..R) + 4 pi int(r g V f'', r=0..R)
            // putting V=-Z/r
            //=- Z 4 pi int(r g f', r=0..R) - Z 8 pi int(g f', r=0..R) - Z 4 pi int(g f'', r=0..R)
            //...


            // sum(x) int(g(r) d/dx V(r) d/dx f(r) dxdydz, -inf..inf)
            //=sum(x)-int((d/dx g(r)) V(r) (d/dx f(r)) dxdydz, -inf..inf) assuming zero bc at infinity
            
            // <i|px.Vpx|j> = <px i| V | px j> ... but this assumes zero bc on both sides

            // px f(r) = (x/r) f'(r)
            // sum(x y z) (px f(r)) (px g(r)) = sum(x y z) (x**2/r**2) f' g' = f' g' which is the same as above

            // int(b[i](r) r^2 d/dr  
            
            const auto& [X,W] = bdataT::quadrature();
            auto D = bdataT::basis().make_deriv_exact_matrix();
            BsplineBasis<T,typename bdataT::knotsT> B1(bdataT::order()-1, bdataT::knots_object());
            auto ones = bdataT::ptr()->ones;
            auto AD = ::madness::inner(B1.tabulate_basis(X), D);
            auto RAD = copy(AD).emul(outer(X,ones));
            auto WRAD = copy(RAD).emul(outer(W,ones));
            print("Z", Z);
            auto PVP = ::madness::inner(WRAD,AD,0,0) * (-Z*0.25*4*constants::pi/(c*c));
            //PVP.emul(dd);

            T PVP00 = c0.trace(::madness::inner(PVP,c0))/S00;
            print("PVP00",PVP00,Z*Z*Z*Z/(4*c*c));
            {
                size_t nbasis = bdataT::nbasis();
                stensorT SX(2*nbasis,2*nbasis), HX(2*nbasis,2*nbasis);
                Slice s0(0,nbasis-1), s1(nbasis,2*nbasis-1);
                SX(s0,s0) = S(s0,s0);
                SX(s1,s1) = KE(s0,s0)*(0.5/(c*c));
                HX(s0,s0) = PE(s0,s0);
                HX(s0,s1) = KE(s0,s0);
                HX(s1,s0) = KE(s0,s0);
                HX(s1,s1) = PVP(s0,s0)-KE(s0,s0);
                stensorT V, e;
                sygv(HX, SX, 1, V, e);
                print("eigenvalues of H sf", e);
                print("exact", exact_energy(Z), e(nbasis), exact_energy(Z)-e(nbasis));

                {
                    print("rsam", bdataT::ptr()->rsample());
                    print("starting to plot");
                    bfunctionT psi(copy(V(s0,nbasis)));
                    T sgn = psi(bdataT::ptr()->rsample()[0]) > 0 ? 1.0 : -1.0;
                    psi *= sgn;
                    bfunctionT phi(sgn*copy(V(s1,nbasis)));
                    bfunctionT psinr(c0);
                    sgn = psinr(bdataT::ptr()->rsample()[0]) > 0 ? 1.0 : -1.0;
                    psinr *= sgn;
                    print("plotting");
                    auto rsam = bdataT::ptr()->rsample();
                    Gnuplot g("set style data lines; set xrange [1e-7:1]; set logscale x","sfdirac.gnuplot");
                    g.plot(rsam, psi(rsam), phi(rsam), psinr(rsam));
                }
            }
        }
    }
};

// The spherical component of the GF is 2 (H(r-s) h(r) j(s) + H(s-r) h(s) j(r)) with
// * H the Heaviside function
// * h(r) = exp(-mu*r)/r
// * j(s) = sinh(mu*s)/(mu*s) = (exp(mu*s) - exp(-mu*s))/(2*mu*s)
// When applying to a function f(r) the first term gives the integral
// (1/mu r) int_0^r s (exp(-mu*(r-s)) - exp(-mu*(r+s))) f(s) ds
// and the second term gives the integral
// (1/mu r) int_r^inf s (exp(-mu*(s-r)) - exp(-mu*(s+r))) f(s) ds
// In both of the these the arguments to the exponentials are negative so nothing blows up.
// These are again computed using GL quadrature.
// int_0^r g(r) = \sum^{whole knot intervals} + \sum^{partial knot interval}
// whole knot intervals are easy because they reuse the current quadrature points
// \int_{interval} = sum_{quadrature points} w_i g(x_i) with w_i and x_ from XW as usual
// For a partial interval, we need additional quadrature points since the range is now
// [lo,r] or [r,hi].  Note you dont need both since [r,hi] = [lo,hi] - [lo,r], the former being the full interval.
// 


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
    std::cout.precision(10);
    //BsplineBasis<float>::test();
    //BsplineBasis<double,KnotsUniform<double>>::test();
    //BsplineBasis<double,KnotsChebyshev<double>>::test();
    //BsplineBasis<double,KnotsGeometricShifted<double>>::test();
    //BsplineBasis<double,KnotsRational<double>>::test();
    BsplineFunction<double>::test();
    //BsplineFunction<double>::test_fit();
    return 0;
}
