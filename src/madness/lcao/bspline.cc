#include <limits>
#include <memory>
#include <algorithm>
#include <madness.h>
#include <madness/misc/gnuplot.h>

// In this file we are translating the Maple code from the paper into C++ using the MADNESS library.
// The Maple code is in the comments, and the C++ code is below it.
// For indexing matrices MADNESS uses the notation A(i,j) just the same as Maple.
// Also, MADNESS uses the notation A(i) to refer to the ith element of a vector, just like Maple.
// Also, MADNESS indexes from 0, not 1, so we have to be very careful about that.
// Also, to index a range MADNESS uses the notation Slice(a,b) which is the same as Maple's a..b.
// Also, MADNESS uses the notation A(i,j) = x to set the value of A(i,j) to x
// Also, MADNESS uses the notation A(i) = x to set the value of A(i) to x

// In this file we are also building a new class to implement the B-spline basis functions and knots.
// The class is called Bspline and it is defined below.

template<typename T> struct is_complex : std::false_type {};
template<typename T> struct is_complex<std::complex<T>> : std::true_type {};
template<typename T> struct scalar_type {typedef T type; };
template<typename T> struct scalar_type<std::complex<T>> {typedef T type;};

using namespace madness; // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

// Locate the interval containing x in the padded knot vector t (increasing order) using linear search.
// Works for general knots on the interval including the endpoints.  Since it using the padded knot vector
// it returns an index in the inclusive range order-1..Npknots-order-1.
template <typename T>
size_t xinterval(T x, const Tensor<T>& t, size_t order) {
    size_t Npknots = t.size();
    if (x < t[0] || x > t[Npknots-1]) {
        print(x, t[0], t[Npknots-1]);
        throw "Xinterval: x not in t";
    }
    for (size_t i=order-1; i<Npknots-order-1; i++) {
        if (x >= t[i] && x < t[i+1]) return i;
    }
    return Npknots - order - 1;
}

// Locate the interval containing x in the padded *uniform* knot vector t (increasing order)
// Note that the knot interval is inclusive of both sides, whereas the general function above is exclusive of the right side.
template <typename T>
size_t xinterval_uniform(T x, const Tensor<T>& t, size_t order) {
    size_t Npknots = t.size();
    if (x < t[order-1] || x > t[Npknots-order]) throw "Xinterval: x not in t";
    T h = ((t[Npknots-order] - t[order-1])/(Npknots - 2* order + 1)) * (1 + std::numeric_limits<T>::epsilon());
    size_t i = size_t((x - t[0])/h ) + order - 1;
    return i;
}

// Uniform points
template <typename T>
Tensor<T> uniform_knots(size_t nknots, T xlo, T xhi) {
    Tensor<T> pts(nknots);
    T h = (xhi - xlo)/(nknots - 1);
    for (size_t i=0; i<nknots; i++) pts[i] = i*h + xlo;
    return pts;
}
    
// Chebyshev points modifyed to include the endpoints
template <typename T>
Tensor<T> cheby_knots(size_t nknots, T xlo, T xhi) {
    Tensor<T> pts(nknots);
    T pi = std::atan(T(1))*4;
    for (size_t i=0; i<nknots; i++) pts[nknots-1-i] = std::cos(i*pi/(nknots - 1));
    pts = (pts + T(1.0))*(T(0.5)*(xhi-xlo)) + xlo;
    return pts;
}

// Oversample the knots by inserting geometric mean of neighboring knots.  For knots on a geometric/exponential grid this is
// what you want, and for a uniform grid we have sqrt(a*(a+h))= a+h/2+O(h^2) where a is a knot and h is the knot spacing.
// So for a dense grid it will also be close to the geometric mean for a>>h.
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

// Needs extending to permit use of optimized xinterval() functions that use knowledge of the knot spacing
template <typename T>
class BsplineBasis {
public:
    const size_t order; // the order of the basis
    const size_t p;     // the degree of the basis = order-1
    const size_t nknots;// the number of unique/unpadded knots
    const size_t npknots;// the number of padded knots
    const size_t nbasis;// the number of basis functions
    const Tensor<T> knots; // the unpadded knots
    const Tensor<T> t; // the padded knots

    // Pad/repeat the knots at the ends to complete the basis
    static Tensor<T> pad_knots(size_t order, const Tensor<T>& knots) {
        long nknots = knots.size();
        Tensor<T> padded(nknots + 2*order - 2);
        for (size_t i=0; i<order-1; i++) padded[i] = knots[0];
        for (size_t i=0; i<nknots; i++) padded[i+order-1] = knots[i];
        for (size_t i=0; i<order-1; i++) padded[i+nknots+order-1] = knots[nknots-1];
        return padded;
    }
    
    BsplineBasis() : order(0) {}

    /// Construct a B-spline basis with the given order and (unpadded) knots.
    BsplineBasis(size_t order, const Tensor<T>& knots)
        : order(order)
        , p(order - 1)
        , nknots(knots.size())
        , npknots(nknots + 2*order - 2)
        , nbasis(nknots + order - 2)
        , knots(copy(knots))
        , t(pad_knots(order, knots))
    {
        if (order < 1) throw "BsplineBasis: order must be >= 1";
        if (order > 10) throw "BsplineBasis: order must be <= 10 (for now)";
        if (knots.size() < 2*order - 1) throw "BsplineBasis: knots must have at least 2*order - 1 elements";
        if (knots[0] > knots[knots.size()-1]) throw "BsplineBasis: knots must be in increasing order";
        if (t.size() != npknots) throw "BsplineBasis: internal error with padded knots";
        if (npknots-order != nbasis) throw "BsplineBasis: internal error with number of basis functions";
    }

    // Stably evaluate a spline expansion at a point, i.e., compute sum(i) c[i] b[i](x)
    template <typename R>
    R deBoor(T x, const Tensor<R>& c) const {
        MADNESS_ASSERT(c.dims()[0] == nbasis);
        size_t k = xinterval(x, t, order);
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
    template <typename R>
    Tensor<R> deBoor(const Tensor<T>& x, const Tensor<R>& c) const {
        MADNESS_ASSERT(c.dims()[0] == nbasis);
        const size_t nsam = x.dims()[0];
        std::vector<size_t> k(nsam);
        Tensor<R> d(p+1,nsam);
        for (size_t i=0; i<nsam; i++) {
            k[i] = xinterval(x[i], t, order);
            d(_,i) =  c(Slice(k[i]-p, k[i]));
        }

        for (size_t r=0; r<p; r++) {
            for (size_t j=p; j>r; j--) {
                for (size_t i=0; i<nsam; i++) {
                    T a = (x[i] - t[j+k[i]-p])/(t[j+k[i]-r] - t[j+k[i]-p]);
                    d(j,i) = (T(1)-a)*d(j-1,i) + a*d(j,i);
                }
            }
        }
        return copy(d(p,_));
    }

    // Evaluate the basis functions at the point x
    Tensor<T> bspline(T x) const {
        size_t k = xinterval(x, t, order);
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
    // highest permitted value is p+1=order.
    // This is the pseudo-inverse of A[i,j] = b[j](Xsample[i]) computed using the non-zero singular values.
    // make_lsq_matrix := proc(t, order, Xsample) local Npknots, Nbasis, A, i, u, s, vt; Npknots := numelems(t); Nbasis := Npknots - order; A := tabulate_basis(t, order, Xsample); u, s, vt := SingularValues(A, output = ['U', 'S', 'Vt'], thin); for i to Nbasis do s[i] := 1.0/s[i]; end do; return Transpose((u . (DiagonalMatrix(s[1 .. Nbasis]))) . vt); end proc;
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
        Tensor<T> A = make_lsq_matrix(Xsample,nzeroL,nzeroR);
        return inner(A,Ysample);
    }

    // // Compute the "exact" derivative of a bspline expansion in terms of splines one degree lower!!!
    // // Note to evaluate this need new set of padded knots with one fewer padding at each end.
    
    // // Note that it is just differentiating the piecewise polyn so boundary conditions are not enforced.
    // // The right way to enfoce b.c.s is to remove the problematic basis functions at the boundary.
    // // deriv_exact := proc(c, t, order) local p, d, j; p := order - 1; d := Vector(numelems(c) - 1); for j to numelems(c) - 1 do d[j] := p*(c[j + 1] - c[j])/(t[j + p + 1] - t[j + 1]); end do; return d; end proc;
    Tensor<T> deriv_exact(const Tensor<T>& c) const {
        MADNESS_ASSERT(c.dims()[0] == nbasis);
        Tensor<T> d(nbasis - 1);
        for (size_t j=0; j<nbasis - 1; j++) {
            d[j] = p*(c[j + 1] - c[j])/(t[j + p + 1] - t[j + 1]);
        }
        return d;
    }

    // Make the (bidiagonal) matrix represenation of the operator that computes the exact derivative
    // make_deriv_exact_matrix := proc(t, order) local p, Nbasis, dMat, j, s; p := order - 1; Nbasis := numelems(t) - order; dMat := Matrix(Nbasis - 1, Nbasis); for j to Nbasis - 1 do s := p/(t[j + p + 1] - t[j + 1]); dMat[j, j + 1] := s; dMat[j, j] := -s; end do; return dMat; end proc;
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
    // make_deriv_matrix := proc(t, order, Xsample) local Npknots, Nbasis, t1, A, M, dMat; Npknots := numelems(t); Nbasis := Npknots - order; t1 := Vector(Npknots + 2); t1[2 .. Npknots + 1] := t; t1[1] := t[1]; t1[Npknots + 2] := t[Npknots]; A := tabulate_basis(t, order, Xsample); M := make_lsq_matrix(t1, order + 1, Xsample); dMat := make_deriv_exact_matrix(t1, order + 1); return (dMat . M) . A; end proc;
    Tensor<T> make_deriv_matrix(const Tensor<T>& Xsample) const {
        size_t nsample = Xsample.dims()[0];
        BsplineBasis B1 = BsplineBasis(order + 1, this->knots);
        Tensor<T> A = tabulate_basis(Xsample);
        Tensor<T> M = B1.make_lsq_matrix(Xsample);
        Tensor<T> dMat = B1.make_deriv_exact_matrix();
        return inner(inner(dMat,M),A);
    }

    // same riff as make_deriv2_matrixX ... need more understanding
    Tensor<T> make_deriv_matrixX(const Tensor<T>& Xsample) const {
        size_t nsample = Xsample.dims()[0];
        BsplineBasis B1 = BsplineBasis(order - 1, this->knots);
        Tensor<T> A = B1.tabulate_basis(Xsample);
        Tensor<T> M = make_lsq_matrix(Xsample);
        Tensor<T> dMat = make_deriv_exact_matrix();
        return inner(inner(M,A),dMat);
    }

    // Make the matrix that applies the second derivative after projecting into a two higher order basis.
    // Gives us 2 order higher accuracy and produces a result in the same order basis as the input.
    Tensor<T> make_deriv2_matrix(const Tensor<T>& Xsample) const {
        size_t nsample = Xsample.dims()[0];
        BsplineBasis B1 = BsplineBasis(order + 1, this->knots);
        BsplineBasis B2 = BsplineBasis(order + 2, this->knots);
        Tensor<T> A = tabulate_basis(Xsample);
        Tensor<T> M2 = B2.make_lsq_matrix(Xsample);
        Tensor<T> dMat2 = B2.make_deriv_exact_matrix();
        Tensor<T> dMat1 = B1.make_deriv_exact_matrix();
        return inner(inner(inner(dMat1,dMat2),M2),A);
    }

    // Make the matrix that applies the second derivative and then projects into a two higher order basis.
    // Produces a result in the same order basis as the input.
    // Need to understand why this sometimes seems better than make_deriv2_matrix which intuitively seems better.
    Tensor<T> make_deriv2_matrixX(const Tensor<T>& Xsample) const {
        size_t nsample = Xsample.dims()[0];
        BsplineBasis B1 = BsplineBasis(order - 1, this->knots);
        BsplineBasis B2 = BsplineBasis(order - 2, this->knots);
        Tensor<T> A = B2.tabulate_basis(Xsample);
        Tensor<T> M = make_lsq_matrix(Xsample);
        Tensor<T> dMat2 = B1.make_deriv_exact_matrix();
        Tensor<T> dMat1 = make_deriv_exact_matrix();
        return inner(M,inner(A,inner(dMat2,dMat1)));
    }

    // Test the basic class functionality --- this is not yet a unit test
    static bool test() {
        size_t nknots = 32;
        size_t order = 8;
        T xlo = -constants::pi;
        T xhi = 2*constants::pi;
        T xtest = 1.77;
        auto knots = uniform_knots(nknots, xlo, xhi);     print("uniform", knots);
        auto cheby = cheby_knots(nknots, xlo, xhi);       print("cheby", cheby);

        //knots = cheby;
        
        auto pad = BsplineBasis<T>::pad_knots(order, knots); print("pad", pad);
        auto xsam = oversample_knots(knots);              print("xsam", xsam);
        long nbasis = nknots + order - 2;
        size_t nsam = xsam.dims()[0];
        
        BsplineBasis<T> B(order, knots);
        
        for (long i=0; i<nbasis; i++) {
            auto b = B.bspline(xtest);
            Tensor<T> c(nbasis); c(i) = 1.0;
            print("deBoor", B.deBoor(xtest, c) - b[i]);
        }
        
        for (size_t i=0; i<nsam; i++) {
            print(xsam[i], xinterval(xsam[i], pad, order), xinterval_uniform(xsam[i], pad, order));
        }

        //print(B.tabulate_basis(knots));
        //print(B.tabulate_basis(xsam));

        Tensor<T> f(nsam), df(nsam), d2f(nsam);
        // for (size_t i=0; i<nsam; i++) {
        //     f[i] = std::exp(-2*xsam[i]);
        //     df[i] = -2*f[i];
        //     d2f[i] = 4*f[i];
        // }
        for (size_t i=0; i<nsam; i++) {
            f[i] = std::sin(xsam[i]);
            df[i] = std::cos(xsam[i]);
            d2f[i] = -std::sin(xsam[i]);
        }
        {
            for (size_t nzeroL=0; nzeroL<=2; nzeroL++) {
                for (size_t nzeroR=0; nzeroR<=2; nzeroR++) {
                    auto M = B.make_lsq_matrix(xsam,nzeroL,nzeroR);
                    auto c1 = inner(M,f);
                    print("sin c", nzeroL, nzeroR, c1);
                }
            }
        }
        auto c = B.fit(xsam,f);
        print(std::sin(1.5), B.deBoor(1.5, c));
        //print(std::exp(-2*1.5), B.deBoor(1.5, c));
        print("|fit-f|", (B.deBoor(xsam, c) - f).normf());

        // Compute first derivative in the lower order basis
        BsplineBasis<T> B1(order - 1, knots);
        auto d = B.deriv_exact(c);
        print(std::cos(1.5), B1.deBoor(1.5, d));
        //print(-2*std::exp(-2*1.5), B1.deBoor(1.5, d));

        auto dMat = B.make_deriv_exact_matrix();
        auto d2 = inner(dMat,c);
        auto df2 = B1.deBoor(xsam, d2);
        Tensor<T> derr1,derr2,derr3;
        derr1 = df2-df;
        print("Err in 'exact' first derivative", (df2-df).normf());
        print(derr1);

        // Compute second derivative in the two lower order basis by applying deriv_exact twice
        BsplineBasis<T> B2(order - 2, knots);
        {
            auto dMat1 = B1.make_deriv_exact_matrix();
            auto d2 = inner(dMat1,inner(dMat,c));
            auto df2 = B2.deBoor(xsam, d2);
            print("Err in 'exact' second derivative", (df2-d2f).normf());

        }
        
        // Compute first derivative in the same order basis
        auto dMat2 = B.make_deriv_matrix(xsam);
        auto d3 = inner(dMat2,c);
        auto df3 = B.deBoor(xsam, d3);
        derr2 = df3-df;
        print(derr2);
        print("Err in fit+1 first derivative", derr2.normf());
        {
            auto dMat2 = B.make_deriv_matrixX(xsam);
            auto d3 = inner(dMat2,c);
            auto df3 = B.deBoor(xsam, d3);
            derr3 = df3-df;
            print(derr3);
            print("Err in fit-1 first derivative", derr3.normf());
        }
        {
            Gnuplot g("set style data l; set xrange [-2.5:5.5]");
            g.plot(xsam, derr1, derr2, derr3);
        }
        // Compute second derivative in the same order basis
        {
            auto d2Mat = B.make_deriv2_matrix(xsam);
            d2 = inner(d2Mat,c);
            df2 = B.deBoor(xsam, d2);
            print("Err in fit+2 second derivative", (df2-d2f).normf());
        }
        {
            auto d2Mat = B.make_deriv2_matrixX(xsam);
            d2 = inner(d2Mat,c);
            df2 = B.deBoor(xsam, d2);
            print("Err in fit-2 second derivative", (df2-d2f).normf());
        }

        return true;
    }
};

// To dos:
// Perhaps Need to eliminate shared_ptr since will be inefficient for large number of functions with threading
// need to modify tests that B=other.B to accomodate type promotion such as float->double

template <typename T>
class BsplineFunction {
public:
    typedef typename TensorTypeData<T>::scalar_type scalar_type;
    typedef BsplineBasis<scalar_type> basisT;
    typedef std::shared_ptr<const basisT> basis_ptrT;
    typedef Tensor<T> tensorT;
    typedef BsplineFunction<T> bfunctionT;

    template <typename U> friend class BsplineFunction;

private:
    basis_ptrT B;
    Tensor<T> c;

public:
    BsplineFunction() = default;
    ~BsplineFunction() = default;

    BsplineFunction(const BsplineFunction<T>& other) // Deep copy constructor
        : B(other.B)
        , c(copy(other.c))
    {}

    BsplineFunction(BsplineFunction<T>&& other) // Move constructor
        : B(std::move(other.B))
        , c(std::move(other.c))
    {}

    template <typename U>
    BsplineFunction(const BsplineFunction<U>& other) // Deep copy constructor
        : B(other.B)
        , c(other.c) // with type conversion don't need to take a copy
    {}

    BsplineFunction& operator=(const bfunctionT& other) { // Deep assignment
        if (this != &other) {
            B = other.B;
            c = copy(other.c);
        }
        return *this;
    }
    
    BsplineFunction(basis_ptrT& B, const tensorT& c = tensorT()) // Copies the coefficients
        : B(B)
        , c(copy(c))
    {}

    BsplineFunction(basis_ptrT& B, tensorT&& c) // Moves the coefficients
        : B(B)
        , c(std::move(c))
    {}

    template <typename U>
    auto operator+(const BsplineFunction<U>& other) const {
        MADNESS_ASSERT(B);
        MADNESS_ASSERT(B == other.B);
        return BsplineFunction<TENSOR_RESULT_TYPE(T,U)>  (B, c + other.c);
    }

    template <typename U>
    BsplineFunction<T> operator+(const T& s) const {
        MADNESS_ASSERT(B);
        return BsplineFunction(B, c + s);
    }


    void operator+=(const BsplineFunction<T>& other) const {
        MADNESS_ASSERT(B);
        MADNESS_ASSERT(B == other.B);
        c += other.c;
    }

    void operator+=(const T& s) const {
        MADNESS_ASSERT(B);
        c += s;
    }

    BsplineFunction<T> operator-(const BsplineFunction<T>& other) const {
        MADNESS_ASSERT(B);
        MADNESS_ASSERT(B == other.B);
        return BsplineFunction(B, c - other.c);
    }

    BsplineFunction<T> operator-(const T& s) const {
        MADNESS_ASSERT(B);
        return BsplineFunction(B, c - s);
    }

    void operator-=(const BsplineFunction<T>& other) const {
        MADNESS_ASSERT(B);
        MADNESS_ASSERT(B == other.B);
        c -= other.c;
    }

    void operator-=(const T& s) const {
        MADNESS_ASSERT(B);
        c -= s;
    }

    BsplineFunction<T> operator*(const T& s) const {
        MADNESS_ASSERT(B);
        return BsplineFunction(B, c*s);
    }

    void operator*=(const T& s) const {
        MADNESS_ASSERT(B);
        c *= s;
    }

    static void test() {
        size_t nknots = 21;
        size_t order = 5;
        T xlo = 0.0;
        T xhi = 13.0;
        auto knots = cheby_knots(nknots, xlo, xhi);
        std::shared_ptr<const basisT> B(new basisT(order, knots));
        BsplineFunction<T> f(B,knots);
        f += f;
        f += 1;

        // Solve H atom
        // 1) Project guess
        // 2) Project potential
        // 3) Compute energy
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
// Make the quadrature points and weights for a rule that employs n GL points between each knot (without padding).  So we will have n*(nknot-1) points.
// make_spline_quadrature := proc(n, Knots) local nknots, x, w, X, W, knot, a, b, i; nknots := numelems(Knots); x, w := make_gauss_legendre(n); X := Vector(n*(nknots - 1)); W := Vector(n*(nknots - 1)); for knot to nknots - 1 do a := Knots[knot]; b := Knots[knot + 1]; for i to n do X[i + n*(knot - 1)] := 1/2*a + 1/2*b + 1/2*x[i]*(b - a); W[i + n*(knot - 1)] := 1/2*(b - a)*w[i]; end do; end do; return X, W; end proc;
// make_overlap_matrix := proc(order, Knots, pKnots) local X, W, A, B, i; X, W := make_spline_quadrature(order + 1, Knots); A := tabulate_basis(pKnots, order, X); B := copy(Transpose(A)); for i to numelems(X) do A[i] := A[i]*W[i]; end do; return B . A; end proc;

int main() {
    //BsplineBasis<float>::test();
    BsplineBasis<double>::test();
    //BsplineFunction<double>::test();
    return 0;
}
