
//*****************************************************************
// Iterative template routine -- BiCGSTAB
//
// BiCGSTAB solves the unsymmetric linear system Ax = b 
// using the Preconditioned BiConjugate Gradient Stabilized method
//
// BiCGSTAB follows the algorithm described on p. 27 of the 
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//  
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//  
//*****************************************************************
// def BiCGSTAB( V, x, b, M, max_iter, tol):

#include <cstdio>
#include <cmath>

template<typename T> class KryFun { 
public:
	virtual void operator() (const T&, T&) =0;
	inline virtual ~KryFun() {};
};

template<typename T, typename Q=double> // T can be an array of Functions
class BiCGStabBase {
public:
	virtual void gaxpy(T&, const Q, const T&, const Q) =0;
	virtual double norm(const T&) =0;
	virtual Q innerp(const T&, const T&) =0;
	virtual void makesurecopy(T&) =0;
	inline virtual ~BiCGStabBase() {};
	int operator() (KryFun<T>& A, T& x, const T& b, int& max_iter, double& tol) {

                //World world(MPI::COMM_WORLD) ;

		// Real resid;
		// Vector rho_1(1), rho_2(1), alpha(1), beta(1), omega(1);
		// Vector p, p, s, shat, t, v;

		double normb = norm(b);

		printf(" in bicgstab, norm b %e\n", normb);//, b, b[0].norm2(),b[1].norm2(),b[2].norm2()

		//  r = b - A*x
		T r(x);
		A(x, r);
		
		//r = b-r; //this line is temporarily replaced by the next line;
		gaxpy(r,-1.0,b,1.0);
		//for (unsigned int i=0; i<r.size(); ++i) r[i]=b[i]-r[i];
		
		//vtrunc_pair(r)
		T rtilde(r);
		makesurecopy(rtilde);

		if (normb == 0.0) normb = 1.0;

		double resid = norm(r)/normb;
		double minr = resid;
		//print minr, normb;
		if (resid <= tol) {
			tol = resid;
			max_iter = 0;
			return 0;
		}

		Q rho_1, rho_2, beta, alpha, omega;
		T v(r), t(r), p;
		makesurecopy(v);
		makesurecopy(t);
		int i;
		for (i=0; i<max_iter; ++i) {
			rho_1 = innerp(r,rtilde);
			if (rho_1 == 0.0) {
				tol = norm(r) / normb;
				max_iter=i;
				return 2;
			}

			if (i == 0) 
				{p = r; makesurecopy(p);} 
			else {
				beta = (rho_1/rho_2) * (alpha/omega);
				//p = r + (p + v.scale(-omega)).scale(beta)
				// p = p.truncate()
				gaxpy(p, 1.0, v, -omega);
			}
			gaxpy(p, beta, r, 1.0);
			//vtrunc_pair(p)

			//    p = M.solve(p);
			// v = A * p
			//p = p
			A(p, v);
			//vtrunc_pair(v);

			alpha = rho_1 / innerp(v,rtilde);
			//s = r + v.scale(-alpha)
			//s = s.truncate()
			gaxpy(r,1.0, v, -alpha);

			resid = norm(r)/normb;
			if (resid < tol) {
				gaxpy(x,1.0, p, alpha);
				//vtrunc_pair(x)

				tol = resid;
				max_iter = i;
				return 0;
			}
			//    shat = M.solve(s)
			// t = A * shat
			//shat = s  using r to save memory
			A(r, t);
			//vtrunc_pair(t)

			omega = innerp(t,r)/innerp(t,t);
			//x = x+ p.scale(alpha) + r.scale(omega)
			//x = x.truncate(.)
			gaxpy(x, 1.0, p, alpha);
			gaxpy(x, 1.0, r, omega);
			//vtrunc_pair(x)

			//r = s + t.scale(-omega)
			//r = r.truncate(.)
			gaxpy(r, 1.0, t, -omega);
			//vtrunc_pair(r)

			rho_2 = rho_1;
			resid = norm(r)/normb;
			if (resid < tol) {
				tol = resid;
				max_iter = i;
				return -1;
			}
			if (omega == 0.0) {
				tol = norm(r) / normb;
				max_iter = i;
				return 3;
			}
			if (resid<minr) {
				minr = resid;
			}
			//printf( "iteration %d residual :%e\n",i, resid);
		}
		tol = resid;
		max_iter = i;
                //world.gop.fence() ;
		return 1;
	};
};

inline double inner(double a, double b) {
	return a*b;
}

template<typename T, typename Q>
class BiCG0: public BiCGStabBase<T,Q> {
public:
	Q innerp(const T& va, const T& vb) {
		return inner(va,vb);
	};
	double norm(const T& vp) {
		return vp.norm2();
	};

	void gaxpy(T& vx, const Q a, const T& vy, const Q b) {
		if (a==1.0) vx += b *vy;
		else {
			vx *= a;
			vx += b *vy;
		}
		
	};
	inline void makesurecopy(T& vx) {
		vx=vx*1.0;
	}
};

template<typename T, typename Q>
class BiCG1: public BiCGStabBase<T,Q> {
public:
        int n;
        void setn(int _n) {n = _n;}
        Q innerp(const T& va, const T& vb) {
                Q p=inner(va[0],vb[0]);
                for (int i=1; i<n; ++i) p += inner(va[i],vb[i]);
                return p;
        };
        double norm(const T& vp) {
                double norm=vp[0].norm2();
                norm *=norm;
                for (int i=1; i<n; ++i) {
                	double dum = vp[i].norm2();
                	norm += dum*dum;
                }
                return sqrt(norm);
        };

        void gaxpy(T& vx, const Q a, const T& vy, const Q b) {
                if (a==1.0) for (int i=0; i<n; ++i) vx[i] += b *vy[i];
                else for (int i=0; i<n; ++i) {
                        vx[i] *= a;
                        vx[i] += b *vy[i];
                }

        };
        inline void makesurecopy(T& vx) {
                for (int i=0; i<n; ++i) vx[i]=vx[i]*1.0;
        }
};

template<typename T, typename Q>
class BiCG2: public BiCGStabBase<T,Q> {
public:
        int n;
        void setn(int _n) {n = _n;}
        Q innerp(const T& va, const T& vb) {
                Q p=inner(va[0],vb[0]);
                for (int i=1; i<n; ++i) p += inner(va[i],vb[i]);
                return p;
        };
        double norm(const T& vp) {
                double norm=vp[0];
                norm *=norm;
                for (int i=1; i<n; ++i) {
                        double dum = vp[i];
                        norm += dum*dum;
                }
                return sqrt(norm);
        };

        void gaxpy(T& vx, const Q a, const T& vy, const Q b) {
                if (a==1.0) for (int i=0; i<n; ++i) vx[i] += b *vy[i];
                else for (int i=0; i<n; ++i) {
                        vx[i] *= a;
                        vx[i] += b *vy[i];
                }

        };
        inline void makesurecopy(T& vx) {
                for (int i=0; i<n; ++i) vx[i]=vx[i]*1.0;
        }
};


