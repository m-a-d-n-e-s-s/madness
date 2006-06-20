#include <iostream>
using std::cout;
using std::endl;

#include <cmath>
using std::abs;
#include <algorithm>

#include <mra/sepop.h>
#include <tensor/tensor.h>
#include <mra/twoscale.h>
#include <mra/legendre.h>

//#include "tensor_lapack.h"

namespace madness {

/// Computes and caches matrix elements for convolution with 1-d Gaussian
    GaussianConvolution::GaussianConvolution(int k, double coeff, double expnt) {
        this->coeff = coeff;
        this->expnt = expnt;
        this->k = k;
        if (!autoc(k,&c)) throw "need to make a consistent exception framework!";
        if (!two_scale_hg(k,&hgT)) throw "need to make a consistent exception framework2!";
        hgT = transpose(hgT);
        npt = k + 11;
        quad_x = Tensor<double>(npt);
        quad_w = Tensor<double>(npt);
        gauss_legendre(npt,0.0,1.0,quad_x.ptr(),quad_w.ptr());
    };

/// Compute the projection of the operator onto the double order polynomials

/// The returned reference is to a cached tensor ... if you want to
/// modify it, take a copy first.
///
/// Return in \c v[p] \c p=0..2*k-1
/// \code
/// r(n,l,p) = 2^(-n) * int(K(2^(-n)*(z+l)) * phi(p,z), z=0..1)
/// \endcode
/// The kernel is coeff*exp(-expnt*z^2).  This is equivalent to
/// \code
/// r(n,l,p) = 2^(-n)*coeff * int( exp(-beta*z^2) * phi(p,z-l), z=l..l+1)
/// \endcode
/// where
/// \code
/// beta = alpha * 2^(-2*n)
/// \endcode
    Tensor<double> GaussianConvolution::rnlp(long n, long l) {
        const int twok = 2*k;
        Tensor<double> v(twok);       // Can optimize this away by passing in

        long lkeep = l;
        if (l<0) l = -l-1;

        /* Apply high-order Gauss Legendre onto subintervals
           
        coeff*int(exp(-beta(x+l)**2) * phi[p](x),x=0..1);

        The translations internally considered are all +ve, so
        signficant pieces will be on the left.  Finish after things
        become insignificant.  

        The resulting coefficients are accurate to about 1e-20.
        */

        // Rescale expnt & coeff onto level n so integration range
        // is [l,l+1]
        double scaledcoeff = coeff*pow(sqrt(0.5),double(n));
        double beta = expnt * pow(0.25,double(n));

        // Subdivide interval into nbox boxes of length h ... estimate
        // appropriate size from the exponent.  A Gaussian with exponent
        // beta falls in magnitude by a factor of 1/e at x=1/sqrt(beta),
        // and by a factor of e^-49 ~ 5e-22 at x=7/sqrt(beta).  So, if we
        // use a box of size 1/sqrt(beta) we will need at most 7 boxes.
        // Incorporate the coefficient into the screening since it may be
        // large.  We can represent exp(-x^2) over [l,l+1] with a polynomial of
        // order 21 over [l,l+1] to a relative precision of better than
        // machine precision for l=0,1,2,3 and for l>3 the absolute error
        // is less than 1e-23.  We want to compute matrix elements with
        // polynomials of order 2*k-1, so the total order is 2*k+20,
        // which can be integrated with a quadrature rule of npt=k+11.
        // npt is set in the constructor.

        double h = 1.0/sqrt(beta);  // 2.0*sqrt(0.5/beta);
        long nbox = long(1.0/h);
        if (nbox < 1) nbox = 1;       // If the exponent is
        h = 1.0/nbox;

        // Find argmax such that h*scaledcoeff*exp(-argmax)=1e-22 ... if
        // beta*xlo*xlo is already greater than argmax we can neglect this
        // and subsequent boxes
        double argmax = fabs(log(1e-22/fabs(scaledcoeff*h)));

        double* phix = new double[twok];
        for (long box=0; box<nbox; box++) {
            double xlo = box*h + l;
            if (beta*xlo*xlo > argmax) break;
            for (long i=0; i<npt; i++) {
                double xx = xlo + h*quad_x(i);
                double ee = scaledcoeff*exp(-beta*xx*xx)*quad_w(i)*h;
                legendre_scaling_functions(xx-l,twok,phix);
                for (long p=0; p<twok; p++) v(p) += ee*phix[p];
            }
        }
        delete[] phix;

        if (lkeep < 0) {
            /* phi[p](1-z) = (-1)^p phi[p](z) */
            for (long p=1; p<twok; p+=2) v(p) = -v(p);
        }

        return v;
    };

/// Computes the transition matrix elements for the convolution for n,l

/// Returns the tensor
/// \code
///   r(i,j) = int(K(x-y) phi[n0](x) phi[nl](y), x=0..1, y=0..1)
/// \endcode
/// This is computed from the matrix elements over the correlation
/// function which in turn are computed from the matrix elements
/// over the double order legendre polynomials.
    const Tensor<double>& GaussianConvolution::rnlij(long n, long lx) {
        Tensor<double>* p=rnlij_cache.getptr(n,lx);
        if (p) return *p;

        long twok = 2*k;
        Tensor<double>  R(2*twok);
        R(Slice(0,twok-1)) = rnlp(n,lx-1);
        R(Slice(twok,2*twok-1)) = rnlp(n,lx);
        R.scale(pow(0.5,0.5*n));
        R = inner(c,R);
        // Enforce symmetry because it seems important
        if (lx == 0)
            for (int i=0; i<k; i++)
                for (int j=0; j<i; j++)
                    R(i,j) = R(j,i) = ((i+j)&1) ? 0.0 : 0.5*(R(i,j)+R(j,i));

        rnlij_cache.set(n,lx,R);
        return *rnlij_cache.getptr(n,lx);
    };

/// Returns a reference to the nonstandard form of the operator
    const Tensor<double>& GaussianConvolution::nonstandard(long n, long lx) {
        Tensor<double> *p = ns_cache.getptr(n,lx);
        if (p) return *p;

        long lx2 = lx*2;
        Tensor<double> R(2*k,2*k);
        Slice s0(0,k-1), s1(k,2*k-1);

        R(s0,s0) = R(s1,s1) = rnlij(n+1,lx2);
        R(s1,s0) = rnlij(n+1,lx2+1);
        R(s0,s1) = rnlij(n+1,lx2-1);

        R = transform(R,hgT);
        // Enforce symmetry because it seems important
        if (lx == 0)
            for (int i=0; i<2*k; i++)
                for (int j=0; j<i; j++)
                    R(i,j) = R(j,i) = ((i+j)&1) ? 0.0 : 0.5*(R(i,j)+R(j,i));

        ns_cache.set(n,lx,transpose(R));
        return *(ns_cache.getptr(n,lx));
    };

/// Returns a reference to the T block of the nonstandard form
    Tensor<double>& GaussianConvolution::nonstandard_T(long n, long lx) {
        Tensor<double> *p = ns_T_cache.getptr(n,lx);
        if (p) return *p;
        Slice s0(0,k-1);
        ns_T_cache.set(n,lx,copy(nonstandard(n,lx)(s0,s0)));
        return *(ns_T_cache.getptr(n,lx));
    };

/// Returns true if the block is expected to be small
    bool GaussianConvolution::issmall(long n, long lx) {
        double beta = expnt*(pow(0.25,double(n)));
        long ll;
        if (lx > 0)
            ll = lx - 1;
        else if (lx < 0)
            ll = -1 - lx;
        else
            ll = 0;

        return (beta*ll*ll > 49.0);      // 49 -> 5e-22
    };

    SeparatedConvolution::SeparatedConvolution(long _k,
            const Tensor<double>& coeffs,
            const Tensor<double>& expnts)
            : k(_k)
            , coeff(coeffs)
            , expnt(expnts)
            , rank(coeffs.dim[0])
            , ops(rank)
    , signs(rank) {
        for (int i=0; i<rank; i++) {
            ops[i] = GaussianConvolution(k,pow(std::abs(coeff[i]),1.0/3.0),expnt[i]);
            signs[i] = (coeff[i]>=0) ? 1 : -1;
        }
    }

/// Returns the norm of the specified block of the NS operator of the mu'th Gaussian
    double SeparatedConvolution::munorm(long mu, long n, long x, long y, long z) {
        long twok = 2*k;
        Tensor<double> f(twok,twok,twok), ff(twok,twok,twok);

        double tol = 1e-20;

        f.fillrandom();
        f.scale(1.0/f.normf());
        double evalp = 1.0, eval, ratio=99.0;
        for (int iter=0; iter<100; iter++) {
            ff.fill(0.0);
            muopxv(mu,n,x,y,z,f,ff,tol);
            f.fill(0.0);
            muopxvt(mu,n,x,y,z,ff,f,tol);

            eval = f.normf();
            if (eval == 0.0) break;
            f.scale(1.0/eval);
            eval = sqrt(eval);
            ratio = eval/evalp;
            //printf("munorm: %d %10.2e %10.2e %10.2e \n", iter, eval, evalp, ratio);
            if (iter>0 && ratio<1.2 && ratio>=1.0) break; // 1.2 was 1.02
            if (iter>10 && eval<1e-30) break;
            evalp = eval;
            if (iter == 99) throw "munorm failed";
        }
        return eval*ratio;
    }

/// Returns the norm of the specified block of the NS operator
    double SeparatedConvolution::norm(long n, long x, long y, long z) {
        long twok = 2*k;
        Tensor<double> f(twok,twok,twok), ff(twok,twok,twok);

        double tol = 1e-20;

        f.fillrandom();
        f.scale(1.0/f.normf());
        double evalp = 1.0, eval, ratio=99.0;
        for (int iter=0; iter<100; iter++) {
            ff.fill(0.0);
            opxv(n,x,y,z,f,ff,tol);
            f.fill(0.0);
            opxvt(n,x,y,z,ff,f,tol);

            eval = f.normf();
            if (eval == 0.0) break;
            f.scale(1.0/eval);
            eval = sqrt(eval);
            ratio = eval/evalp;
            //printf("norm: %d %10.2e %10.2e %10.2e \n", iter, eval, evalp, ratio);
            if (iter>0 && ratio<1.2 && ratio>=1.0) break; // 1.2 was 1.02
            if (iter>10 && eval<1e-30) break;
            evalp = eval;
            if (iter == 99) throw "norm failed";
        }
        return eval*ratio;
    }

/// Apply the specified block to a vector accumulating into the result
    void SeparatedConvolution::opxv(long n, long x, long y, long z,
                                    const Tensor<double>& f, Tensor<double>& result,
                                    double tol) {
        for (long mu=0; mu<rank; mu++)
            muopxv(mu,n,x,y,z,f,result,tol);
    }

/// Apply transpose of specified block to a vector accumulating into the result
    void SeparatedConvolution::opxvt(long n, long x, long y, long z,
                                     const Tensor<double>& f, Tensor<double>& result,
                                     double tol) {
        for (long mu=0; mu<rank; mu++)
            muopxvt(mu,n,x,y,z,f,result,tol);
    }

/// Apply the specified block of the mu'th term to a vector accumulating into the result
    void SeparatedConvolution::muopxv(const long mu,
                                      const long n, const long x, const long y, const long z,
                                      const Tensor<double>& f, Tensor<double>& result,
                                      const double tol) {
        const Tensor<double>& X = ops[mu].nonstandard(n,x);
        const Tensor<double>& Y = ops[mu].nonstandard(n,y);
        const Tensor<double>& Z = ops[mu].nonstandard(n,z);

        //   cout << "muopxv " << " " << mu << " " << n << " " <<  x << " " << y << " " << z << " " << f.normf() << endl;
        //   cout << X;
        //   cout << endl;
        //   cout << Y;
        //   cout << endl;
        //   cout << Z;
        //   cout << endl;

        //   cout << result.normf() << endl;

        //   Tensor<double> t0 = ::inner(f,X,0,0);
        //   Tensor<double> t1 = ::inner(t0,Y,0,0);
        //   Tensor<double> t2 = ::inner(t1,Z,0,0);

        //   cout << t0.normf() << " " << t1.normf() << " " << t2.normf() << endl;

        //   Tensor<double> tmp0(2*k,2*k,2*k);
        //   for (int p=0; p<2*k; p++)
        //     for (int q=0; q<2*k; q++)
        //       for (int r=0; r<2*k; r++)
        //    for (int s=0; s<2*k; s++)
        //      tmp0(p,q,r) += f(p,q,s)*Z(s,r);

        //   Tensor<double> tmp1(2*k,2*k,2*k);
        //   for (int p=0; p<2*k; p++)
        //     for (int q=0; q<2*k; q++)
        //       for (int r=0; r<2*k; r++)
        //    for (int s=0; s<2*k; s++)
        //      tmp1(p,q,r) += tmp0(p,s,r)*Y(s,q);

        //   Tensor<double> tmp2(2*k,2*k,2*k);
        //   for (int p=0; p<2*k; p++)
        //     for (int q=0; q<2*k; q++)
        //       for (int r=0; r<2*k; r++)
        //    for (int s=0; s<2*k; s++)
        //      tmp2(p,q,r) += tmp1(s,q,r)*X(s,p);

        //   cout << tmp0.normf() << " " << tmp1.normf() << " " << tmp2.normf() << endl;

        //   cout << "matrix norms " << f.normf() << " " << X.normf() << " " << Y.normf() << " " << Z.normf() << endl;
        //   Tensor<double> U, sval, VT;
        //   svd(X,&U,&sval,&VT);
        //   cout << "sval\n" << sval << endl;
        //   cout << "first sing vec\n" << U(_,0) << endl << VT(0,_) << endl;


        //   cout << "\n" << f;
        //   cout << "\n" << tmp2;
        //   cout << (t2-tmp2).normf() << endl;

        //   cout << "input transformed into SVD vectors\n";
        //   cout << ::inner(::inner(::inner(f,U,0,0),U,0,0),U,0,0);

        //   result.gaxpy(1.0,tmp2,1.0);


        if (signs[mu] > 0)
            inner_result(madness::inner(madness::inner(f,X,0,0),Y,0,0),Z,0,0,result);
        else
            throw "DUDE!  The sign!";

        //   cout << "muopxv intermediate norm " << tmp2.normf() << endl;

        if (n > 0) {
            Slice s0(0,k-1);
            Tensor<double> r_T = Tensor<double>(k,k,k);
            const Tensor<double>& X_T = ops[mu].nonstandard_T(n,x);
            const Tensor<double>& Y_T = ops[mu].nonstandard_T(n,y);
            const Tensor<double>& Z_T = ops[mu].nonstandard_T(n,z);

            if (signs[mu] > 0)
                inner_result(madness::inner(madness::inner(copy(f(s0,s0,s0)),X_T,0,0),Y_T,0,0),Z_T,0,0,r_T);
            else
                throw "DUDE!  The sign!";

            result(s0,s0,s0) -= r_T;
        }
    }


/// Apply transpose of specified block of the mu'th term to a vector accumulating into the result
    void SeparatedConvolution::muopxvt(long mu, long n, long x, long y, long z,
                                       const Tensor<double>& f, Tensor<double>& result,
                                       double tol) {
        const Tensor<double> X = transpose(ops[mu].nonstandard(n,x));
        const Tensor<double> Y = transpose(ops[mu].nonstandard(n,y));
        const Tensor<double> Z = transpose(ops[mu].nonstandard(n,z));

        if (signs[mu] > 0)
            madness::inner_result(madness::inner(madness::inner(f,X,0,0),Y,0,0),Z,0,0,result);
        else
            throw "DUDE!  The sign!";

        if (n > 0) {
            Slice s0(0,k-1);

            Tensor<double> r_T = Tensor<double>(k,k,k);
            const Tensor<double>& X_T = transpose(ops[mu].nonstandard_T(n,x));
            const Tensor<double>& Y_T = transpose(ops[mu].nonstandard_T(n,y));
            const Tensor<double>& Z_T = transpose(ops[mu].nonstandard_T(n,z));

            if (signs[mu] > 0)
                madness::inner_result(madness::inner(madness::inner(f(s0,s0,s0),X_T,0,0),Y_T,0,0),Z_T,0,0,r_T);
            else
                throw "DUDE!  The sign!";

            result(s0,s0,s0) -= r_T;
        }
    }

}
