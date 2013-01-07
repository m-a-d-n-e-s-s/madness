#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cassert>
#include <cstdio>
#include <qd/qd_real.h>
#include <ls.h>

using namespace std;

// template <typename FLOAT>
// FLOAT myexp(const FLOAT& x) {
//     FLOAT r = exp(x);
//     if (isnan(r)) r = numeric_limits<FLOAT>::safe_max();
//     return r;
//     // static const FLOAT xmax = log();
//     // static const FLOAT xmin = -xmax;
//     // if (x <= xmin) return 0;
//     // else if (x >= xmax) return numeric_limits<FLOAT>::safe_max();
//     // else return exp(x);
// }

template <typename FLOAT>
class matrix {
    std::vector<FLOAT> v;
    int n, m;
public:
    matrix() : v(), n(0), m(0) {}

    matrix(int n, int m, const FLOAT& value=0) : v(n*m), n(n), m(m) 
    {
        for (int i=0; i<n*m; i++) v[i] = value;
    }

    FLOAT& operator()(int i, int j) {
        return v[i*m + j];
    }

    const FLOAT& operator()(int i, int j) const {
        return v[i*m + j];
    }

    int get_n() const {return n;}

    int get_m() const {return m;}
};

template <typename t>
std::ostream& operator<<(std::ostream& s, const matrix<t>& c) {
    int n = c.get_n();
    int m = c.get_m();

    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            s << c(i,j) << " ";
        }
        if (i != (n-1)) s << endl;
    }

    return s;
}

/// easy printing of vectors
template <typename t>
std::ostream& operator<<(std::ostream& s, const std::vector<t>& c) {
    s << "[";
    typename std::vector<t>::const_iterator it = c.begin();
    while (it != c.end()) {
        s << *it;
        ++it;
        if (it != c.end()) s << ", ";
    };
    s << "]";
    return s;
}

template <typename FLOAT>
FLOAT b(int j, FLOAT z) {
    if (j == 0) return 0;
    else if (j&1) return z;
    else return 1;
}

template <typename FLOAT>
FLOAT a(int j, FLOAT z) {
    static const FLOAT v = 1;
    if (j==1) return exp(-z);
    else if (j&1) return j/2;
    else return v+(j-1)/2;
}


template <typename FLOAT>
FLOAT cfrac(FLOAT z) {
    static const FLOAT eps = 2.0*numeric_limits<FLOAT>::epsilon();
    static const FLOAT tiny = eps*eps;
    static const FLOAT one = 1;
    FLOAT f = b(0,z);
    if (f == 0) f = tiny;

    FLOAT cjm1 = f;
    FLOAT djm1 = 0;
    
    for (int j=1; j<100000; j++) {
        FLOAT aj = a(j,z), bj = b(j,z);
        FLOAT dj = bj + aj * djm1;
        if (dj == 0) dj = tiny;
        dj = one/dj;

        FLOAT cj = bj + aj/cjm1;
        if (cj == 0) cj = tiny;

        FLOAT delta = cj * dj;
        f *= delta;

        //cout << "iter " << j << " " << f << endl;
        if (fabs(delta-one) < eps) return f;

        cjm1 = cj;
        djm1 = dj;
    }
    throw "too many iterations";
}

template <typename FLOAT>
FLOAT series(FLOAT z) {
    static const FLOAT eps = 2.0*numeric_limits<FLOAT>::epsilon();
    static const FLOAT euler = 0.57721566490153286060651209008240243104215933593992359880576723488;
    FLOAT sum = -euler - log(z) + z;
    FLOAT term = z;
    int n=1;
    while (1) {
        n++;
        term *= -z/n;
        FLOAT r = term/n;
        sum += r;
        //cout << n << " " << sum << endl;
        if (fabs(r) < sum*eps) return sum;
    }
}

template <typename FLOAT>
FLOAT e1(FLOAT z) {
    static const FLOAT two = 2;
    static const FLOAT zero = 0;
    if (z>two) return cfrac(z);
    else if (z>zero) return series(z);
    else throw "e1(z) must have z>0";
}



template <typename t>
t convert(const char* c) {
    return c;
}

template <>
double convert<double>(const char* c) {
    return atof(c);
}

template<>
float convert<float>(const char* c) {
    return atof(c);
}


template <typename FLOAT>
void jacobi(int n, matrix<FLOAT>& a, vector<FLOAT>& e, matrix<FLOAT>& v) {
    static const FLOAT zero = convert<FLOAT>("0.0");
    static const FLOAT half = convert<FLOAT>("0.5");
    static const FLOAT one = convert<FLOAT>("1.0");
    static const FLOAT two = convert<FLOAT>("2.0");
    static const FLOAT rsqrt2 = one/sqrt(two);

    e = vector<FLOAT>(n,zero);
    v = matrix<FLOAT>(n,n);
    
    const FLOAT tolmin = 2.0*numeric_limits<FLOAT>::epsilon();
    //cout << "tolmin " << tolmin << endl;
    FLOAT tol = convert<FLOAT>("1e-2");

    const FLOAT tolscale = convert<FLOAT>("0.1");
    const FLOAT quart = convert<FLOAT>("0.25");

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            v(i,j) = zero;
        }
        v(i,i) = one;
    }

    FLOAT maxd = zero ;
    for (int i=0; i<n; i++) {
        maxd = max(fabs(a(i,i)),maxd) ;
    }
    //cout << "maxd " << maxd << endl;

    int nrot;
    int nrotsum = 0;
    //cout << "iter, tol, maxdaij, nrot" << endl;

    for (int iter=0; iter<50; iter++) {
        FLOAT maxdaij = zero;
        nrot = 0;
        for (int i=0; i<n; i++)  { 
            for (int j=i+1; j<n; j++) {
                FLOAT aii = a(i,i);
                FLOAT ajj = a(j,j);
                FLOAT aij = a(i,j);
                FLOAT daij = fabs(aij);
                maxdaij = max(maxdaij,daij/maxd);
                //cout << "testing " << daij << " " << tol*maxd << endl;
                if (daij > tol*maxd) { // screen small elements 
                    FLOAT s = ajj - aii;
                    FLOAT ds = fabs(s);
                    if (daij > (tolmin*ds)) { // check for sufficient precision
                        FLOAT c;
                        nrot++;
                        if (tolmin*daij > ds) {
                            c = s = rsqrt2;
                        } 
                        else { 
                            FLOAT t = aij/s;
                            FLOAT u = quart/sqrt(quart+t*t);
                            c = sqrt(half+u);
                            s = two*t*u/c;
                        }
                        for (int k=0; k<=i; k++) {
                            FLOAT t = a(k,i);
                            FLOAT u = a(k,j);
                            a(k,i) = c*t - s*u;
                            a(k,j) = s*t + c*u;
                        }

                        for (int k=i+1; k<j; k++) {
                            FLOAT t = a(i,k);
                            FLOAT u = a(k,j) ;
                            a(i,k) = c*t - s*u ;
                            a(k,j) = s*t + c*u;
                        }
                        a(j,j) = s*aij + c*ajj;
                        a(i,i) = c*a(i,i) - s*(c*aij - s*ajj);
                          
                        for (int k=j; k<n; k++) {
                            FLOAT t = a(i,k);
                            FLOAT u = a(j,k);
                            a(i,k) = c*t - s*u ;
                            a(j,k) = s*t + c*u ;
                        }
                        for (int k=0; k<n; k++) {
                            FLOAT t = v(i,k) ;
                            FLOAT u = v(j,k) ;
                            v(i,k) = c*t - s*u ;
                            v(j,k) = s*t + c*u ;
                        }
                        a(j,i) = a(i,j) = zero ;
                        maxd = max(max(maxd,fabs(a(i,i))),fabs(a(j,j)));
                    }
                }
            }
        }

        //cout << iter << " " <<  tol << " " <<  maxdaij << " " <<  nrot << endl;

        nrotsum += nrot;
        if (nrot == 0 && tol <= tolmin) break;

        tol = min(min(tol,maxdaij*tolscale),maxdaij*maxdaij);
        tol = max(tol,tolmin);
    }
    if (nrot != 0) { 
        cout << "jacobi iteration did not converge in 50 passes" << endl;
    }

    // sort and copy evals

    for (int i=0; i<n; i++) {
        e[i] = a(i,i);
        for (int j=0; j<i; j++) {
            if (e[j] > e[i]) {
                 swap(e[j],e[i]);
                 for (int k=0; k<n; k++) {
                     swap(v(j,k),v(i,k));
                 }
            }
        }
    }

}

// solves ax=b for symmetric a discarding eigen values less than tol
template <typename FLOAT>
void LSQ(int n, const matrix<FLOAT>& a, const vector<FLOAT>& b, vector<FLOAT>& x, const FLOAT& tol) {
    static const FLOAT zero = convert<FLOAT>("0.0");
    matrix<FLOAT> acopy(a);
    matrix<FLOAT> v;
    vector<FLOAT> s;

    jacobi(n, acopy, s, v);

    x = vector<FLOAT>(n,zero);

    //for (int i=0; i<n; i++) cout << "s " << i << " " << s[i] << endl;

    for (int k=0; k<n; k++) {
        if (fabs(s[k]) >= tol) {
            FLOAT sum = zero;
            for (int j=0; j<n; j++) {
                sum += v(k,j) * b[j];
            }
            sum /= s[k];
            for (int i=0; i<n; i++) {
                x[i] += v(k,i) * sum;
            }
        }
        else {
            cout << "lsq discarding singular value " << s[k] << endl;
        }
    }
}

// template <typename FLOAT>
// void mxm(int n, const FLOAT* a, const FLOAT* b, FLOAT* c) {
//     for (int i=0; i<n; i++) {
//         for (int j=0; j<n; j++) {
//             FLOAT sum = 0.0;
//             for (int k=0; k<n; k++) {
//                 sum += a[ind(i,k)] * b[ind(k,j)];
//             }
//             c[ind(i,j)] = sum;
//         }
//     }
// }

//#define RECIPX
#ifdef RECIPX
// evaluate the target function
template <typename FLOAT>
FLOAT target(const FLOAT& x) {
    return 1/x;
}

// evaluate the target function and its first & second derivativse
template <typename FLOAT>
void target(const FLOAT& x, FLOAT& f, FLOAT& d1f, FLOAT& d2f) {
    f = 1/x;
    d1f = -f*f;
    d2f = -2*f*d1f;
}

// evaluate the weight function
template <typename FLOAT>
FLOAT weight(const FLOAT& x) {
    return x;
}

// evaluate the weight function and its first & second derivativse
template <typename FLOAT>
void weight(const FLOAT& x, FLOAT& f, FLOAT& d1f, FLOAT& d2f) {
    f = x;
    d1f = 1.0;
    d2f = 0.0;
}
#else
// evaluate the target function
template <typename FLOAT>
FLOAT target(const FLOAT& x) {
    return 1/sqrt(x);
}

// evaluate the target function and its first & second derivativse
template <typename FLOAT>
void target(const FLOAT& x, FLOAT& f, FLOAT& d1f, FLOAT& d2f) {
    f = 1/sqrt(x);
    d1f = -0.5*f/x;
    d2f = -1.5*d1f/x;
}

// evaluate the weight function
template <typename FLOAT>
FLOAT weight(const FLOAT& x) {
    return sqrt(x);
}

// evaluate the weight function and its first & second derivativse
template <typename FLOAT>
void weight(const FLOAT& x, FLOAT& f, FLOAT& d1f, FLOAT& d2f) {
    f = sqrt(x);
    d1f = 0.5*f/x;
    d2f = -0.5*d1f/x;
}

#endif


// evaluate the function fitted to exponentials
template <typename FLOAT>
FLOAT fit(const FLOAT& x, const vector<FLOAT>& p) {
    static const FLOAT zero= convert<FLOAT>("0.0");
    int n = p.size()/2;
    FLOAT sum = zero;
    for (int i=0; i<n; i++) {
        sum += exp(p[i] - x*p[i+n]);
    }
    return sum;
}

// evaluate the function and first/second derivatives fitted to exponentials
template <typename FLOAT>
void fit(const FLOAT& x, const vector<FLOAT>& p, FLOAT& f, FLOAT& g, FLOAT& h) {
    static const FLOAT zero= convert<FLOAT>("0.0");
    int n = p.size()/2;
    f = zero;
    g = zero;
    h = zero;
    for (int i=0; i<n; i++) {
        FLOAT term = exp(p[i] - x*p[i+n]);
        f += term;
        g -= p[i+n]*term;
        h += p[i+n]*p[i+n]*term;
    }
}


template <typename FLOAT>
void plot(int npt, const FLOAT& a, const FLOAT& b, const vector<FLOAT>& p) 
{
    FLOAT h = pow(b/a,FLOAT(1)/(npt-1));
    FLOAT x = a;
    for (int i=0; i<npt; i++) {
        cout << x << " " <<  (fit(x,p) - target(x))*weight(x) << endl;
        x *= h;
    }
}


// basis function mu is g[mu](x) = exp(p[mu] - x*p[mu+n])
// 
// non-zero first derivatives
// dg[mu]/dp[mu] = g[mu]
// dg[mu]/dp[mu+n] = -x*g[mu]
//
// non-zero second derivatives 
// d2g[mu]/dp[mu]dp[mu] = g[mu]
// d2g[mu]/dp[mu]dp[mu+n] = -x*g[mu]
// d2g[mu]/dp[mu+n]dp[mu+n] = x^2*g[mu]
//
// g(i,mu) = g[mu](x[i])
// eps[i] = sum(mu, g(i,mu)) - f(x[i])
//
// d0 = sum(i=0..npt-2, (eps[i]+eps[i+1])^2)
// d1_mu = dd0/dp[mu] = 2 sum(i, eps[i] * 
// 
// x[2*n+1]
// p[2*n]
// d1[2*n]
// d2[2*n,2*n]
template <typename FLOAT>
void makedata(const vector<FLOAT>& x, const vector<FLOAT>& f, const vector<FLOAT>& p, 
              FLOAT& d0, vector<FLOAT>& d1, matrix<FLOAT>& d2, bool d0only=false) 
{
    static const FLOAT zero = convert<FLOAT>("0.0");
    const int n = p.size()/2; // number of exponentials
    const int npt = x.size();    // number of values x
    //assert((unsigned) npt == x.size());

    
    // make basis functions and error values at grid points
    matrix<FLOAT> g(npt,n);
    vector<FLOAT> eps(npt);
    vector<FLOAT> w(npt);
    for (int i=0; i<npt; i++) {
        FLOAT sum = zero;
        for (int mu=0; mu<n; mu++) {
            g(i,mu) = exp(p[mu] - x[i]*p[mu+n]); // note constraint of positive coeffs
            sum += g(i,mu);
        }
        w[i] = weight(x[i]);
        eps[i] = (sum - f[i])*w[i];
    }

    d0 = zero;
    for (int i=0; i<npt-1; i++) 
        d0 += (eps[i]+eps[i+1])*(eps[i]+eps[i+1]);

    if (d0only) return;

    // make derivatives of eps
    matrix<FLOAT> deps1(npt, 2*n, zero);  // d/dp[mu] 
    matrix<FLOAT> deps2a(npt, n, zero); // d2/dp[mu]^2
    matrix<FLOAT> deps2b(npt, n, zero); // d2/dp[mu]dp[mu+n]
    matrix<FLOAT> deps2c(npt, n, zero); // d2/dp[mu+n]dp[mu+n]
    for (int i=0; i<npt; i++) {
        for (int mu=0; mu<n; mu++) {
            deps1(i,mu) = g(i,mu)*w[i];
            deps1(i,mu+n) = -x[i]*g(i,mu)*w[i];
            deps2a(i,mu) = g(i,mu)*w[i];
            deps2b(i,mu) = -x[i]*g(i,mu)*w[i];
            deps2c(i,mu) = x[i]*x[i]*g(i,mu)*w[i];
        }
    }

    // assemble derivatives of the target function now
    // just using the chain rule
    d1 = vector<FLOAT>(2*n,zero);
    d2 = matrix<FLOAT>(2*n,2*n,zero);

    for (int i=0; i<(npt-1); i++) {
        for (int mu=0; mu<2*n; mu++) {
            d1[mu] += 2 * (eps[i]+eps[i+1]) * (deps1(i,mu) + deps1(i+1,mu));
            for (int nu=0; nu<2*n; nu++) {
                d2(mu,nu) += 2 * (deps1(i,mu)+deps1(i+1,mu)) * (deps1(i,nu) + deps1(i+1,nu));
            }
        }
    }
    
    for (int i=0; i<(npt-1); i++) {
        for (int mu=0; mu<n; mu++) {
            d2(mu,mu) += 2 * (eps[i]+eps[i+1]) * (deps2a(i,mu) + deps2a(i+1,mu));
            
            d2(mu+n,mu) += 2 * (eps[i]+eps[i+1]) * (deps2b(i,mu) + deps2b(i+1,mu));
            d2(mu,mu+n) += 2 * (eps[i]+eps[i+1]) * (deps2b(i,mu) + deps2b(i+1,mu));
            
            d2(mu+n,mu+n) += 2 * (eps[i]+eps[i+1]) * (deps2c(i,mu) + deps2c(i+1,mu));
        }
    }
}

template <typename FLOAT> 
class Func {
    const vector<FLOAT>& x;
    const vector<FLOAT>& f;
    const vector<FLOAT>& p;
    const vector<FLOAT>& dp;
    mutable vector<FLOAT> pnew;

public:
    Func(const vector<FLOAT>& x, const vector<FLOAT>&f, const vector<FLOAT>& p, const vector<FLOAT>& dp) 
        : x(x), f(f), p(p), dp(dp), pnew(dp.size())
    {}

    FLOAT operator()(const FLOAT& s) const {
        for (unsigned int i=0; i<p.size(); i++) {
            pnew[i] = p[i] - s*dp[i];
        }
        FLOAT d0;
        vector<FLOAT> d1;
        matrix<FLOAT> d2;
        makedata(x, f, pnew, d0, d1, d2, true);

        if (isnan(d0)) d0 = numeric_limits<FLOAT>::safe_max();

        return d0;
    }
};

template <typename FLOAT> 
vector<FLOAT> opt(const vector<FLOAT>& x, const vector<FLOAT>& f, const vector<FLOAT>& w, const vector<FLOAT>& guess, int maxiter)
{
    static const FLOAT zero= convert<FLOAT>("0.0");
    static const FLOAT one = convert<FLOAT>("1.0");
    static const FLOAT half = convert<FLOAT>("0.5");
    vector<FLOAT> p = guess;
    const int n = p.size();

    const int nprint = 100;
    
    for (int iter=0; iter<maxiter; iter++) {
        const bool print = (iter%nprint) == 0;
        FLOAT d0;
        vector<FLOAT> d1;
        matrix<FLOAT> d2;

        makedata(x, f, p, d0, d1, d2, false);
        
        if (print) {
            cout << "\n\niteration " << iter << " " << d0 << endl << endl;
            //cout << " p " << p << endl;
            //cout << " d " << d1 << endl << endl;
        }


        bool alternating = true;
        FLOAT errprev = 0;
        FLOAT convtol = 0;
        for (unsigned int i=0; i<x.size(); i++) {
            FLOAT err = (fit(x[i],p)-f[i])*w[i];
            if (i>0 && err*errprev>0) alternating=false;
            errprev = err;
            convtol += err*err;
            if (print) cout << " err " << i << " " << x[i] << " " << err << endl;
        }
        convtol *= 1e-5 / n;
        if (print) cout << iter << " " << d0 << " " << alternating << endl;

        if (iter > 10 && d0 < convtol && alternating) {
            cout << "\nconverged " << iter << " " << d0 << " " << convtol << endl << endl;
            break;
        }
        
        matrix<FLOAT> v;
        vector<FLOAT> e;
        jacobi(n, d2, e, v);

        //if (print) {cout << "eigenvalues " << e << endl;}

        vector<FLOAT> dp(n,zero);
        FLOAT tol;

        //if ((iter%10) == 0) tol = e[n/3];
        //else if ((iter%3) == 1) tol = e[2*n/3];
        //else 
            tol = 1000.0*numeric_limits<FLOAT>::epsilon();
        
        for (int k=0; k<n; k++) {
            if (e[k] < -tol) {
                //cout << " making negative eigenvalue positive " << k << " " << e[k] << endl;
                e[k] = one/min(one,-e[k]);
            }
            else if (e[k] < tol) {
                //cout << " small eigenvalue " << k << " " << e[k] << endl;
                e[k] = one/tol;
            }
            else {
                e[k] = one/e[k];
            }

            FLOAT sum = zero;
            for (int j=0; j<n; j++) {
                sum += v(k,j) * d1[j];
            }
            //if (print) {cout << " dmode : " << k << " " << sum << endl;}
            sum *= e[k];
            for (int i=0; i<n; i++) {
                dp[i] += v(k,i) * sum;
            }
        }

        //cout << "\ndp " << dp << endl << endl;

        FLOAT pnorm=zero, dpnorm=zero;
        for (int i=0; i<n; i++) {
            pnorm += p[i]*p[i];
            dpnorm += dp[i]*dp[i];
        }
        pnorm = sqrt(pnorm);
        dpnorm = sqrt(dpnorm);

        vector<FLOAT> pnew(n);
        FLOAT s = min(one, 0.05*pnorm/dpnorm);
        FLOAT d0new;
        
        // exponents must be positive
        for (int i=n/2; i<n; i++) {
            if (p[i] - s*dp[i] < 0) s = min(s,p[i]/(5*dp[i]));
        }
        
        Func<FLOAT> func(x, f, p, dp);

        do  {
            
            d0new = func(s);
            //cout << "        s " << s << "  d0 old " << d0 << "  d0 new " << d0new << endl;

            if (d0new < d0 && !isnan(d0new)) break;
            else s *= half;
        } while (s > 1e-20);

        FLOAT ax=0, bx=s, cx, fa=d0, fb=d0new, fc, snew;
        mnbrak(&ax, &bx, &cx, &fa, &fb, &fc, func);
        d0new = brent(ax, bx, cx, func, FLOAT(1e-6), &s);

        for (int i=0; i<n; i++) {
            pnew[i] = p[i] - s*dp[i];
        }


        p = pnew;
        if (print) cout << "      s " << s << "  d0 old " << d0 << "  d0 new " << d0new << endl;
    }
    return p;
}


template <typename FLOAT>
vector<FLOAT> updatex(const vector<FLOAT>& x, const vector<FLOAT>& p) 
{
    static const FLOAT one = 1;

    const FLOAT a = x[0];
    const FLOAT b = x[x.size()-1];
    const int npt = 30*x.size();
    const FLOAT h = pow(b/a, one/npt);
    cout << "update x: h = " << h << endl;

    // bracket min/max then use newton
    //  f(x-h) f(x) f(x+h)
    //        >    >       decreasing to right
    //        <    >       maximum
    //        >    <       minimum
    //        <    <       increasing to right
    
    vector<FLOAT> xnew;
    xnew.push_back(a);
    FLOAT f2m = fit(a,p)   - target(a);
    FLOAT f1m = fit(a*h,p) - target(a*h);
    FLOAT xx = a*h*h;
    while (xx < b) {
        FLOAT f0m = (fit(xx,p) - target(xx))*weight(xx);

        if ( (f2m<f1m && f1m>f0m) || (f2m>f1m && f1m<f0m) ) {
            FLOAT f, df1, df2, t, dt1, dt2, w, dw1, dw2, g, dg1, dg2;
            FLOAT yy = xx;
            // Three iterations of Newton to converge
            fit(yy, p, f, df1, df2);
            target(yy, t, dt1, dt2);
            weight(yy, w, dw1, dw2);
            g = (f - t)*w;
            dg1 = (df1 - dt1)*w + (f-t)*dw1;
            dg2 = (df2 - dt2)*w + (f-t)*dw2 + 2*(df1 - dt1)*dw1;

            //cout << "UPDATE A: " << yy << " " <<g << " " << dg1 << " " << dg2 << " " << dg1/dg2 << endl;

            yy -= dg1/dg2;

            fit(yy, p, f, df1, df2);
            target(yy, t, dt1, dt2);
            weight(yy, w, dw1, dw2);
            g = (f - t)*w;
            dg1 = (df1 - dt1)*w + (f-t)*dw1;
            dg2 = (df2 - dt2)*w + (f-t)*dw2 + 2*(df1 - dt1)*dw1;

            //cout << "UPDATE B: " << yy << " " <<g << " " << dg1 << " " << dg2 << " " << dg1/dg2 << endl;

            yy -= dg1/dg2;

            fit(yy, p, f, df1, df2);
            target(yy, t, dt1, dt2);
            weight(yy, w, dw1, dw2);
            g = (f - t)*w;
            dg1 = (df1 - dt1)*w + (f-t)*dw1;
            dg2 = (df2 - dt2)*w + (f-t)*dw2 + 2*(df1 - dt1)*dw1;

            //cout << "UPDATE C: " <<  yy << " " << g << " " << dg1 << " " << dg2 << " " << dg1/dg2 << endl;

            yy -= dg1/dg2;


            xnew.push_back(yy);
        }

        xx *= h;
        f2m = f1m;
        f1m = f0m;
    }
    xnew.push_back(b);

    if (xnew.size() != p.size()+1) {
        cout << x << endl;
        cout << xnew << endl;
        throw "darn it";
    }

    return xnew;
}

template <typename FLOAT>
void test() {
    const FLOAT zero= convert<FLOAT>("0.0");
    int n = 5;
    matrix<FLOAT> a(n,n), acopy(n,n);
    vector<FLOAT> x(n), xcopy(n), b(n);
    
    cout.precision(numeric_limits<FLOAT>::digits10);

    for (int i=0; i<n; i++) {
        xcopy[i] = sqrt(FLOAT(i+2));
        for (int j=0; j<=i; j++) {
            a(i,j) = a(j,i) = FLOAT(1 + i*j)/FLOAT(i*i + j*j + i*j + 1);
        }
    }

    for (int i=0; i<n; i++) {
        FLOAT sum = 0.0;
        for (int j=0; j<n; j++) {
            acopy(i,j) = a(i,j);
            sum += a(i,j)*xcopy[j];
        }
        b[i] = sum;
    }
    
    matrix<FLOAT> u;
    vector<FLOAT> e;
    jacobi(n, a, e, u);
    
    FLOAT err = zero;
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            FLOAT sum = zero;
            for (int k=0; k<n; k++) {
                sum += u(k,i) * e[k] * u(k,j);
            }
            err += fabs(sum - acopy(i,j));
        }
    }
    
    cout << "    error in recontructed matrix " << err << endl;
    
    // for (int i=0; i<n; i++) {
    //     cout << "eigenvalue " << i << " " << e[i] << endl;
    //     for (int j=0; j<n; j++) {
    //         cout << "    " << u[i*n+j] << endl;
    //     }
    // }

    const FLOAT tol = 100.0*numeric_limits<FLOAT>::epsilon();
    LSQ(n, acopy, b, x, tol);
    
    err = zero;
    for (int i=0; i<n; i++) err += fabs(x[i] - xcopy[i]);
    cout << "    error in linear solution " << err << endl;

    // 
    {
        static const FLOAT one = 1;
        static const FLOAT pi = 4*atan(FLOAT(1));
        cout << pi << endl;
        const FLOAT a = 1e-14;
        const FLOAT b = 1e0;
        const int nfunc = 20;
        const int nx = 2*nfunc + 1;
        vector<FLOAT> x(nx), p(2*nfunc), f(nx), w(nx);
        if (nx < 4) throw "nx must be >= 4";
#ifdef RECIPX
        FLOAT rx = pow(b/a,one/(nx - 1 - 2));
        FLOAT xxx;
        x[0] = a;
        x[1] = a*pow(rx,one/3);
        x[2] = xxx = a*rx;
        for (int i=2; i<nx-2; i++) {
            x[i] = xxx;
            xxx *= rx;
        }
        x[nx-1] = b;
        x[nx-2] = b/pow(rx,one/3);
        x[nx-3] = b/rx;
#else
        FLOAT rx = pow(b/a,one/(nx - 1));
        FLOAT xxx;
        xxx = x[0] = a;
        for (int i=1; i<nx-1; i++) {
            xxx *= rx;
            x[i] = xxx;
        }
        x[nx-1] = b;
#endif

        // cheby points
        // for (int i=1; i<=nx; i++) {
        //     x[i-1] = exp(0.5*(log(a)+log(b)) + 0.5*(log(a)-log(b))*cos((2*i-1)*pi/(2*nx)));
        // }
        //x[0] = a;
        //x[nx-1] = b;

        for (int i=0; i<nx; i++) {
            f[i] = target(x[i]);
            w[i] = weight(x[i]);
        }

        FLOAT tmax = 2.3/(x[1]-x[0]);
        FLOAT tmin = 0.11/(b-a);
        FLOAT dt = pow(tmax/tmin,one/(nfunc-1));
        
        for (int mu=0; mu<nfunc; mu++) {
#ifdef RECIPX
            p[mu] = log(tmin*dt);
#else
            p[mu] = log(tmin*dt)*0.5;
#endif            
            p[mu+nfunc] = tmin;
            tmin *= dt;
        }

        FLOAT d0;
        vector<FLOAT> d1;
        matrix<FLOAT> d2;
        makedata(x, f, p, d0, d1, d2, false);

        cout << "x" << endl;
        cout << x << endl << endl;
        cout << "f" << endl;
        cout << f << endl << endl;
        cout << "p" << endl;
        cout << p << endl << endl;
        // cout << "d0" << endl;
        // cout << d0 << endl << endl;
        // cout << "d1" << endl;
        // cout << d1 << endl << endl;
        // cout << "d2" << endl;
        // cout << d2 << endl << endl;
        
        // FLOAT h(0.01);//sqrt(numeric_limits<FLOAT>::epsilon()));
        // FLOAT plus, minus, pp, mp, pm, mm;

        // cout << "numerical d1\n";
        // for (int mu=0; mu<2*nfunc; mu++) {
        //     p[mu] += h;
        //     makedata(x, f, p, plus, d1, d2, true);
        //     p[mu] -= h+h;
        //     makedata(x, f, p, minus, d1, d2, true);
        //     p[mu] += h;
        //     //cout << plus << " " << minus << endl;
        //     cout << (plus-minus)/(2*h) << " ";
        // }
        // cout << endl << endl;

        // cout << "numerical d2\n";
        // for (int mu=0; mu<2*nfunc; mu++) {
        //     for (int nu=0; nu<2*nfunc; nu++) {
        //         if (mu == nu) {
        //             p[mu] += h;
        //             makedata(x, f, p, plus, d1, d2, true);
        //             p[mu] -= h+h;
        //             makedata(x, f, p, minus, d1, d2, true);
        //             p[mu] += h;
        //             cout << (plus+minus-2*d0)/(h*h) << " ";
        //         }
        //         else {
        //             p[mu] += h; p[nu] += h;
        //             makedata(x, f, p, pp, d1, d2, true);
        //             p[mu] -= h+h;
        //             makedata(x, f, p, mp, d1, d2, true);
        //             p[mu] += h+h;
        //             p[nu] -= h+h;
        //             makedata(x, f, p, pm, d1, d2, true);
        //             p[mu] -= h+h;
        //             makedata(x, f, p, mm, d1, d2, true);
        //             p[nu] += h;
        //             p[mu] += h;
        //             cout << (pp + mm - mp - pm)/(4*h*h) << " ";
        //         }
        //     }
        //     cout << endl;
        // }
        // cout << endl;

        for (int iter=0; iter<10; iter++) {
            p = opt(x, f, w, p, 2000000);
            cout << "Converged parameters" << endl;
            for (int mu=0; mu<nfunc; mu++) {
                cout << mu << " " << p[mu] << " " << p[mu+nfunc] << endl;
            }
            cout << endl;
            plot(1001, a, b, p);
            cout << endl;
            x = updatex(x, p);
            cout << "Updated x " << x << endl;
            for (int i=0; i<nx; i++) {
                f[i] = target(x[i]);
                w[i] = weight(x[i]);
            }
        }
    }
}

int main() {
    //std::cout << "FLOAT " << endl;
    //test<float>();
    //std::cout << "DOUBLE " << endl;
    //test<double>();
    std::cout << "DoubleDOUBLE " << endl;
    test<dd_real>();
    //std::cout << "QuadDOUBLE " << endl;
    //test<qd_real>();
    return 0;
}
