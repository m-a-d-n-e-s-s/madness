/*
 This file is part of MADNESS.

 Copyright (C) 2007,2010 Oak Ridge National Laboratory

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

 For more information please contact:

 Robert J. Harrison
 Oak Ridge National Laboratory
 One Bethel Valley Road
 P.O. Box 2008, MS-6367

 email: harrisonrj@ornl.gov
 tel:   865-241-3937
 fax:   865-572-0680


 $Id: key.h 2907 2012-06-14 10:15:05Z 3ru6ruWu $
 */

#ifndef MADNESS_MRA_GFIT_H__INCLUDED
#define MADNESS_MRA_GFIT_H__INCLUDED

/// \file gfit.h
/// \brief fit isotropic functions to a set of Gaussians with controlled precision

//#include <iostream>

#include <cmath>
#include "../constants.h"
#include "../tensor/basetensor.h"
#include "../tensor/slice.h"
#include "../tensor/tensor.h"
#include "../tensor/tensor_lapack.h"
#include "../world/madness_exception.h"
#include "../world/print.h"
#include <madness/mra/operatorinfo.h>


namespace madness {

template<typename T, std::size_t NDIM>
class GFit {

public:

	/// default ctor does nothing
	GFit() = default;

    GFit(OperatorInfo info) {
        double mu = info.mu;
        double lo = info.lo;
        double hi = info.hi;
        MADNESS_CHECK_THROW(hi>0,"hi must be positive in gfit: U need to set it manually in operator.h");
        double eps = info.thresh;
        bool debug=info.debug;
        OpType type = info.type;


        if (type==OT_G12) {*this=CoulombFit(lo,hi,eps,debug);
        } else if (type==OT_SLATER) {*this=SlaterFit(mu,lo,hi,eps,debug);
        } else if (type==OT_GAUSS)  {*this=GaussFit(mu,lo,hi,eps,debug);
        } else if (type==OT_F12)    {*this=F12Fit(mu,lo,hi,eps,debug);
        } else if (type==OT_FG12)   {*this=FGFit(mu,lo,hi,eps,debug);
        } else if (type==OT_F212)   {*this=F12sqFit(mu,lo,hi,eps,debug);
        } else if (type==OT_F2G12)  {*this=F2GFit(mu,lo,hi,eps,debug);
        } else if (type==OT_BSH)    {*this=BSHFit(mu,lo,hi,eps,debug);
        } else {
            print("Operator type not implemented: ",type);
            MADNESS_EXCEPTION("Operator type not implemented: ",1);
        }

    }

    /// copy constructor
    GFit(const GFit& other) = default;

    /// assignment operator
    GFit& operator=(const GFit& other) {
        coeffs_ = other.coeffs_;
        exponents_ = other.exponents_;
        return *this;
    }

	/// return a fit for the Coulomb function

    ///  f(r) = 1/r
	static GFit CoulombFit(double lo, double hi, double eps, bool prnt=false) {
		GFit fit=BSHFit(0.0,lo,hi,eps/(4.0*constants::pi),prnt);
		fit.coeffs_.scale(4.0*constants::pi);   // BSHFit scales by 1/(4 pi), undo here
		return fit;
	}

	/// return a fit for the bound-state Helmholtz function

	/// the BSH function is defined by
	///  f(r) = 1/ (4 pi) exp(-\mu r)/r
	/// @param[in]	mu	the exponent of the BSH
	/// @param[in]	lo	the smallest length scale that needs to be precisely represented
	/// @param[in]	hi	the largest length scale that needs to be precisely represented
	/// @param[in]	eps	the precision threshold
	/// @parma[in]	prnt	print level
	static GFit BSHFit(double mu, double lo, double hi, double eps, bool prnt=false) {
		GFit fit;
        bool fix_interval=false;
		if (NDIM==3) bsh_fit(mu,lo,hi,eps,fit.coeffs_,fit.exponents_,prnt,fix_interval);
		else bsh_fit_ndim(NDIM,mu,lo,hi,eps,fit.coeffs_,fit.exponents_,prnt);

        if (prnt) {
            print("bsh fit");
            auto exact = [&mu](const double r) -> double { return 1.0/(4.0 * constants::pi) * exp(-mu * r)/r; };
            fit.print_accuracy(exact, lo, hi);
        }
		return fit;
	}

	/// return a fit for the Slater function

	/// the Slater function is defined by
	///  f(r) = exp(-\gamma r)
	/// @param[in]	gamma	the exponent of the Slater function
	/// @param[in]	lo	the smallest length scale that needs to be precisely represented
	/// @param[in]	hi	the largest length scale that needs to be precisely represented
	/// @param[in]	eps	the precision threshold
	/// @parma[in]	prnt	print level
	static GFit SlaterFit(double gamma, double lo, double hi, double eps, bool prnt=false) {
		GFit fit;
		slater_fit(gamma,lo,hi,eps,fit.coeffs_,fit.exponents_,false);
        if (prnt) {
            print("Slater fit");
            auto exact = [&gamma](const double r) -> double { return exp(-gamma * r); };
            fit.print_accuracy(exact, lo, hi);
        }
		return fit;
	}

    /// return a (trivial) fit for a single Gauss function

    /// the Gauss function is defined by
    ///  f(r) = exp(-\gamma r^2)
    /// @param[in]	gamma	the exponent of the Gauss function
    /// @param[in]	lo	the smallest length scale that needs to be precisely represented
    /// @param[in]	hi	the largest length scale that needs to be precisely represented
    /// @param[in]	eps	the precision threshold
    /// @parma[in]	prnt	print level
    static GFit GaussFit(double gamma, double lo, double hi, double eps, bool prnt=false) {
        GFit fit;
        fit.coeffs_=Tensor<T>(1);
        fit.exponents_=Tensor<double>(1);
        fit.coeffs_=1.0;
        fit.exponents_=gamma;
        if (prnt) {
            print("Gauss fit");
            auto exact = [&gamma](const double r) -> double { return exp(-gamma * r*r); };
            fit.print_accuracy(exact, lo, hi);
        }
        return fit;
    }

    /// return a fit for the F12 correlation factor

    /// the Slater function is defined by
    ///  f(r) = 1/(2 gamma) * (1 - exp(-\gamma r))
    /// @param[in]	gamma	the exponent of the Slater function
    /// @param[in]	lo	the smallest length scale that needs to be precisely represented
    /// @param[in]	hi	the largest length scale that needs to be precisely represented
    /// @param[in]	eps	the precision threshold
    /// @parma[in]	prnt	print level
    static GFit F12Fit(double gamma, double lo, double hi, double eps, bool prnt=false) {
        GFit fit;
        f12_fit(gamma,lo*0.1,hi,eps*0.01,fit.coeffs_,fit.exponents_,false);
        fit.coeffs_*=(0.5/gamma);
        if (prnt) {
            print("f12 fit");
            auto exact=[&gamma](const double r) -> double {return 0.5/gamma*(1.0-exp(-gamma*r));};
            fit.print_accuracy(exact,lo,hi);
        }
        return fit;
    }

    /// return a fit for the F12^2 correlation factor

    /// the Slater function square is defined by
    ///  f(r) = [ 1/(2 gamma) * (1 - exp(-\gamma r)) ] ^2
    /// @param[in]	gamma	the exponent of the Slater function
    /// @param[in]	lo	the smallest length scale that needs to be precisely represented
    /// @param[in]	hi	the largest length scale that needs to be precisely represented
    /// @param[in]	eps	the precision threshold
    /// @parma[in]	prnt	print level
    static GFit F12sqFit(double gamma, double lo, double hi, double eps, bool prnt=false) {
        GFit fit;
        f12sq_fit(gamma,lo*0.1,hi,eps*0.01,fit.coeffs_,fit.exponents_,false);
        fit.coeffs_*=(0.25/(gamma*gamma));
        if (prnt) {
            print("f12sq fit");
            auto exact=[&gamma](const double r) -> double {return std::pow(0.5/gamma*(1.0-exp(-gamma*r)),2.0);};
            fit.print_accuracy(exact,lo,hi);
        }
        return fit;
    }


    /// return a fit for the FG function

    /// fg = 1/(2 mu) * (1 - exp(-gamma r12))  / r12
    ///    = 1/(2 mu) *( 1/r12 - exp(-gamma r12)/r12)
    ///    = 1/(2 mu) * (coulomb - bsh)
    /// @param[in]	gamma	the exponent of the Slater function
    /// @param[in]	lo	the smallest length scale that needs to be precisely represented
    /// @param[in]	hi	the largest length scale that needs to be precisely represented
    /// @param[in]	eps	the precision threshold
    /// @parma[in]	prnt	print level
    static GFit FGFit(double gamma, double lo, double hi, double eps, bool prnt=false) {
        GFit bshfit,coulombfit;
        eps*=0.1;
        lo*=0.1;
//        bool restrict_interval=false;
        bool fix_interval=true;
        bsh_fit(gamma,lo,hi,eps,bshfit.coeffs_,bshfit.exponents_,false,fix_interval);
        bsh_fit(0.0,lo,hi,eps,coulombfit.coeffs_,coulombfit.exponents_,false,fix_interval);
        // check the exponents are identical
        auto diffexponents=(coulombfit.exponents() - bshfit.exponents());
        MADNESS_CHECK(diffexponents.normf()/coulombfit.exponents().size()<1.e-12);
        auto diffcoefficients=(coulombfit.coeffs() - bshfit.coeffs());
        GFit fgfit;
        fgfit.exponents_=bshfit.exponents_;
        fgfit.coeffs_=4.0*constants::pi*0.5/gamma*diffcoefficients;
        GFit<double,3>::prune_small_coefficients(eps,lo,hi,fgfit.coeffs_,fgfit.exponents_);

        if (prnt) {
            print("fg fit");
            auto exact=[&gamma](const double r) -> double {return 0.5/gamma*(1.0-exp(-gamma*r))/r;};
            fgfit.print_accuracy(exact,lo,hi);
        }
        return fgfit;
    }

    /// return a fit for the F2G function

    /// f2g = (1/(2 mu) * (1 - exp(-gamma r12)))^2  / r12
    ///    = 1/(4 mu^2) * [ 1/r12 - 2 exp(-gamma r12)/r12) + exp(-2 gamma r12)/r12 ]
    /// @param[in]	gamma	the exponent of the Slater function
    /// @param[in]	lo	the smallest length scale that needs to be precisely represented
    /// @param[in]	hi	the largest length scale that needs to be precisely represented
    /// @param[in]	eps	the precision threshold
    /// @parma[in]	prnt	print level
    static GFit F2GFit(double gamma, double lo, double hi, double eps, bool prnt=false) {
        GFit bshfit,coulombfit,bsh2fit;
        eps*=0.1;
        lo*=0.1;
//        bool restrict_interval=false;
        bool fix_interval=true;
        bsh_fit(gamma,lo,hi,eps,bshfit.coeffs_,bshfit.exponents_,false,fix_interval);
        bsh_fit(2.0*gamma,lo,hi,eps,bsh2fit.coeffs_,bsh2fit.exponents_,false,fix_interval);
        bsh_fit(0.0,lo,hi,eps,coulombfit.coeffs_,coulombfit.exponents_,false,fix_interval);

        // check the exponents are identical
        auto diffexponents=(coulombfit.exponents() - bshfit.exponents());
        MADNESS_CHECK(diffexponents.normf()/coulombfit.exponents().size()<1.e-12);
        auto diffexponents1=(coulombfit.exponents() - bsh2fit.exponents());
        MADNESS_CHECK(diffexponents1.normf()/coulombfit.exponents().size()<1.e-12);

        auto coefficients=(coulombfit.coeffs() - 2.0* bshfit.coeffs() + bsh2fit.coeffs());
        GFit f2gfit;
        f2gfit.exponents_=bshfit.exponents_;
        // additional factor 4 pi due to implementation of bsh_fit
        double fourpi=4.0*constants::pi;
        double fourmu2=4.0*gamma*gamma;
        f2gfit.coeffs_=fourpi/fourmu2*coefficients;
        GFit<double,3>::prune_small_coefficients(eps,lo,hi,f2gfit.coeffs_,f2gfit.exponents_);

        if (prnt) {
            print("fg fit");
            auto exact=[&gamma](const double r) -> double {
                return 0.25/(gamma*gamma)*(1.0-2.0*exp(-gamma*r)+exp(-2.0*gamma*r))/r;
            };
            f2gfit.print_accuracy(exact,lo,hi);
        }
        return f2gfit;
    }


	/// return a fit for a general isotropic function

	/// note that the error is controlled over a uniform grid, the boundaries
	/// will be poorly represented in general. Following Beylkin 2005
	static GFit GeneralFit() {
		MADNESS_EXCEPTION("General GFit still to be implemented",1);
		return GFit();
	}

	/// return the coefficients of the fit
	Tensor<T> coeffs() const {return coeffs_;}

	/// return the exponents of the fit
	Tensor<T> exponents() const {return exponents_;}

    void static prune_small_coefficients(const double eps, const double lo, const double hi,
                                  Tensor<double>& coeff, Tensor<double>& expnt) {
        double mid = lo + (hi-lo)*0.5;
        long npt=coeff.size();
        long i;
        for (i=npt-1; i>0; --i) {
            double cnew = coeff[i]*exp(-(expnt[i]-expnt[i-1])*mid*mid);
            double errlo = coeff[i]*exp(-expnt[i]*lo*lo) -
                           cnew*exp(-expnt[i-1]*lo*lo);
            double errhi = coeff[i]*exp(-expnt[i]*hi*hi) -
                           cnew*exp(-expnt[i-1]*hi*hi);
            if (std::max(std::abs(errlo),std::abs(errhi)) > 0.03*eps) break;
            npt--;
            coeff[i-1] = coeff[i-1] + cnew;
        }
        coeff = coeff(Slice(0,npt-1));
        expnt = expnt(Slice(0,npt-1));
    }

    void truncate_periodic_expansion(Tensor<double>& c, Tensor<double>& e,
			double L, bool discardG0) const {
		double tcut = 0.25/L/L;

		if (discardG0) {
			// Relies on expnts being in decreasing order
			for (int i=0; i<e.dim(0); ++i) {
				if (e(i) < tcut) {
					c = c(Slice(0,i));
					e = e(Slice(0,i));
					break;
				}
			}
		} else {
//			// Relies on expnts being in decreasing order
//			int icut = -1;
//			for (int i=0; i<e.dim(0); ++i) {
//				if (e(i) < tcut) {
//					icut = i;
//					break;
//				}
//			}
//			if (icut > 0) {
//				for (int i=icut+1; i<e.dim(0); ++i) {
//					c(icut) += c(i);
//				}
//				c = c(Slice(0,icut));
//				e = e(Slice(0,icut));
//			}
		}
	}

private:

	/// ctor taking an isotropic function

	/// the function will be represented with a uniform error on a uniform grid
	/// @param[in]	f	a 1d-function that implements T operator()
	template<typename funcT>
	GFit(const funcT& f) {

	}

    /// print coefficients and exponents, and values and errors

    /// @param[in]  op  the exact function, e.g. 1/r, exp(-mu r), etc
    /// @param[in]  lo  lower bound for the range r
    /// @param[in]  hi  higher bound for the range r
    template<typename opT>
    void print_accuracy(opT op, const double lo, const double hi) const {

        std::cout << "weights and roots" << std::endl;
        for (int i=0; i<coeffs_.size(); ++i)
            std::cout << i << " " << coeffs_[i] << " " << exponents_[i] << std::endl;

        long npt = 300;
        std::cout << "       x         value   abserr   relerr" << std::endl;
        std::cout << "  ------------  ------- -------- -------- " << std::endl;
        double step = exp(log(hi/lo)/(npt+1));
        for (int i=0; i<=npt; ++i) {
            double r = lo*(pow(step,i+0.5));
//            double exact = exp(-mu*r)/r/4.0/constants::pi;
            double exact = op(r);
            double test = 0.0;
            for (int j=0; j<coeffs_.dim(0); ++j)
                test += coeffs_[j]*exp(-r*r*exponents_[j]);
            double err = 0.0;
            if (exact) err = (exact-test)/exact;
            printf("  %.6e %8.1e %8.1e %8.1e\n",r, exact, exact-test, err);
        }
    }

	/// the coefficients of the expansion f(x) = \sum_m coeffs[m] exp(-exponents[m] * x^2)
	Tensor<T> coeffs_;

	/// the exponents of the expansion f(x) = \sum_m coeffs[m] exp(-exponents[m] * x^2)
	Tensor<T> exponents_;

	/// fit the function exp(-mu r)/r

	/// formulas taken from
	/// G. Beylkin and L. Monzon, On approximation of functions by exponential sums,
	/// Appl Comput Harmon A, vol. 19, no. 1, pp. 17-48, Jul. 2005.
	/// and
	/// R. J. Harrison, G. I. Fann, T. Yanai, and G. Beylkin,
	/// Multiresolution Quantum Chemistry in Multiwavelet Bases,
	/// Lecture Notes in Computer Science, vol. 2660, p. 707, 2003.
	static void bsh_fit(double mu, double lo, double hi, double eps,
			Tensor<double>& pcoeff, Tensor<double>& pexpnt, bool prnt, bool fix_interval) {

        if (mu < 0.0) throw "cannot handle negative mu in bsh_fit";
//        bool restrict_interval=(mu>0) and use_mu_for_restricting_interval;


		if ((mu > 0) and (not fix_interval)) {
//        if (restrict_interval) {
			// Restrict hi according to the exponential decay
			double r = -log(4*constants::pi*0.01*eps);
			r = -log(r * 4*constants::pi*0.01*eps);
			if (hi > r) hi = r;
		}

		double TT;
		double slo, shi;

		if (eps >= 1e-2) TT = 5;
		else if (eps >= 1e-4) TT = 10;
		else if (eps >= 1e-6) TT = 14;
		else if (eps >= 1e-8) TT = 18;
		else if (eps >= 1e-10) TT = 22;
		else if (eps >= 1e-12) TT = 26;
		else TT = 30;

        if ((mu > 0) and (not fix_interval)) {
//		if (restrict_interval) {
			slo = -0.5*log(4.0*TT/(mu*mu));
		}
		else {
			slo = log(eps/hi) - 1.0;
		}
		shi = 0.5*log(TT/(lo*lo));
                if (shi <= slo) throw "bsh_fit: logic error in slo,shi";

		// Resolution required for quadrature over s
		double h = 1.0/(0.2-.50*log10(eps)); // was 0.5 was 0.47

		// Truncate the number of binary digits in h's mantissa
		// so that rounding does not occur when performing
		// manipulations to determine the quadrature points and
		// to limit the number of distinct values in case of
		// multiple precisions being used at the same time.
		h = floor(64.0*h)/64.0;


		// Round shi/lo up/down to an integral multiple of quadrature points
		shi = ceil(shi/h)*h;
		slo = floor(slo/h)*h;

		long npt = long((shi-slo)/h+0.5);

		//if (prnt)
                //std::cout << "mu " << mu << " slo " << slo << " shi " << shi << " npt " << npt << " h " << h << " " << eps << std::endl;

		Tensor<double> coeff(npt), expnt(npt);

		for (int i=0; i<npt; ++i) {
			double s = slo + h*(npt-i);	// i+1
			coeff[i] = h*2.0/sqrt(constants::pi)*exp(-mu*mu*exp(-2.0*s)/4.0)*exp(s);
			coeff[i] = coeff[i]/(4.0*constants::pi);
			expnt[i] = exp(2.0*s);
		}

#if ONE_TERM
		npt=1;
		double s=1.0;
		coeff[0]=1.0;
		expnt[0] = exp(2.0*s);
		coeff=coeff(Slice(0,0));
		expnt=expnt(Slice(0,0));
		print("only one term in gfit",s,coeff[0],expnt[0]);


#endif

		// Prune large exponents from the fit ... never necessary due to construction

		// Prune small exponents from Coulomb fit.  Evaluate a gaussian at
		// the range midpoint, and replace it there with the next most
		// diffuse gaussian.  Then examine the resulting error at the two
		// end points ... if this error is less than the desired
		// precision, can discard the diffuse gaussian.

        if ((mu == 0.0) and (not fix_interval)) {
//		if (restrict_interval) {
            GFit<double,3>::prune_small_coefficients(eps,lo,hi,coeff,expnt);
		}

		// Modify the coeffs of the largest exponents to satisfy the moment conditions
		//
		// SETTING NMOM>1 TURNS OUT TO BE A BAD IDEA (AS CURRENTLY IMPLEMENTED)
		// [It is accurate and efficient for a one-shot application but it seems to
		//  introduce fine-scale noise that amplifies during iterative solution of
		//  the SCF and DFT equations ... the symptom is negative coeffs in the fit]
		//
		// SET NMOM=0 or 1 (1 recommended) unless you are doing a one-shot application
		//
		// Determine the effective range of the four largest exponents and compute
		// moments of the exact and remainder of the fit.  Then adjust the coeffs
		// to reproduce the exact moments in that volume.
		//
		// If nmom!=4 we assume that we will eliminate n=-1 which is stored first
		// in the moment list
		//
		// <r^i|gj> cj = <r^i|exact-remainder>
		const long nmom = 1;
		if (nmom > 0) {
			Tensor<double> q(4), qg(4);
			double range = sqrt(-log(1e-6)/expnt[nmom-1]);
			if (prnt) print("exponent(nmom-1)",expnt[nmom-1],"has range", range);

			bsh_spherical_moments(mu, range, q);
			Tensor<double> M(nmom,nmom);
			for (int i=nmom; i<npt; ++i) {
				Tensor<double> qt(4);
				gaussian_spherical_moments(expnt[i], range, qt);
				qg += qt*coeff[i];
			}
			if (nmom != 4) {
				q = q(Slice(1,nmom));
				qg = qg(Slice(1,nmom));
			}
			if (prnt) {
				print("moments", q);
				print("moments", qg);
			}
			q = q - qg;
			for (int j=0; j<nmom; ++j) {
				Tensor<double> qt(4);
				gaussian_spherical_moments(expnt[j], range, qt);
				if (nmom != 4) qt = qt(Slice(1,nmom));
				for (int i=0; i<nmom; ++i) {
					M(i,j) = qt[i];
				}
			}
			Tensor<double> ncoeff;
			gesv(M, q, ncoeff);
			if (prnt) {
				print("M\n",M);
				print("old coeffs", coeff(Slice(0,nmom-1)));
				print("new coeffs", ncoeff);
			}

			coeff(Slice(0,nmom-1)) = ncoeff;
		}

		pcoeff = coeff;
		pexpnt = expnt;
	}

	/// fit a Slater function using a sum of Gaussians

	/// formula inspired by the BSH fit, with the roles of r and mu exchanged
	/// see also Eq. (A3) in
	/// S. Ten-no, Initiation of explicitly correlated Slater-type geminal theory,
	/// Chem. Phys. Lett., vol. 398, no. 1, pp. 56-61, 2004.
	static void slater_fit(double gamma, double lo, double hi, double eps,
			Tensor<double>& pcoeff, Tensor<double>& pexpnt, bool prnt) {

        MADNESS_CHECK(gamma >0.0);
		// empirical number TT for the upper integration limit
		double TT;
		if (eps >= 1e-2) TT = 5;
		else if (eps >= 1e-4) TT = 10;
		else if (eps >= 1e-6) TT = 14;
		else if (eps >= 1e-8) TT = 18;
		else if (eps >= 1e-10) TT = 22;
		else if (eps >= 1e-12) TT = 26;
		else TT = 30;

		// integration limits for quadrature over s: slo and shi
        // slo and shi must not depend on gamma!!!
		double slo=0.5 * log(eps) - 1.0;
		double shi=log(TT/(lo*lo))*0.5;

		// Resolution required for quadrature over s
		double h = 1.0/(0.2-.5*log10(eps)); // was 0.5 was 0.47

		// Truncate the number of binary digits in h's mantissa
		// so that rounding does not occur when performing
		// manipulations to determine the quadrature points and
		// to limit the number of distinct values in case of
		// multiple precisions being used at the same time.
		h = floor(64.0*h)/64.0;

		// Round shi/lo up/down to an integral multiple of quadrature points
		shi = ceil(shi/h)*h;
		slo = floor(slo/h)*h;

		long npt = long((shi-slo)/h+0.5);

		Tensor<double> coeff(npt), expnt(npt);

		for (int i=0; i<npt; ++i) {
			const double s = slo + h*(npt-i);	// i+1
			coeff[i] = h*exp(-gamma*gamma*exp(2.0*s) + s);
			coeff[i]*=2.0*gamma/sqrt(constants::pi);
			expnt[i] = 0.25*exp(-2.0*s);
		}

		if (prnt) {
			std::cout << "weights and roots for a Slater function with gamma=" << gamma << std::endl;
			for (int i=0; i<npt; ++i)
				std::cout << i << " " << coeff[i] << " " << expnt[i] << std::endl;

			long npt = 300;
			//double hi = 1.0;
			//if (mu) hi = min(1.0,30.0/mu);
			std::cout << "       x         value   abserr   relerr" << std::endl;
			std::cout << "  ------------  ------- -------- -------- " << std::endl;
			double step = exp(log(hi/lo)/(npt+1));
			for (int i=0; i<=npt; ++i) {
				double r = lo*(pow(step,i+0.5));
				double exact = exp(-gamma*r);
				double test = 0.0;
				for (int j=0; j<coeff.dim(0); ++j)
					test += coeff[j]*exp(-r*r*expnt[j]);
				double err = 0.0;
				if (exact) err = (exact-test)/exact;
				printf("  %.6e %8.1e %8.1e %8.1e\n",r, exact, exact-test, err);
			}
		}
		pcoeff = coeff;
		pexpnt = expnt;
	}


    /// fit a correlation factor (1- exp(-mu r))

    /// use the Slater fit with an additional term: 1*exp(-0 r^2)
    static void f12_fit(double gamma, double lo, double hi, double eps,
                           Tensor<double>& pcoeff, Tensor<double>& pexpnt, bool prnt) {
        Tensor<double> coeff,expnt;
        slater_fit(gamma, lo, hi, eps, coeff, expnt, prnt);

        pcoeff=Tensor<double>(coeff.size()+1);
        pcoeff(Slice(1,-1,1))=-coeff(_);
        pexpnt=Tensor<double>(expnt.size()+1);
        pexpnt(Slice(1,-1,1))=expnt(_);

        pcoeff(0l)=1.0;
        pexpnt(0l)=1.e-10;
    }


    /// fit a correlation factor f12^2 = (1- exp(-mu r))^2 = 1 - 2 exp(-mu r) + exp(-2 mu r)

    /// no factor 1/(2 mu) or square of it included!
    /// use the Slater fit with an additional term: 1*exp(-0 r^2)
    static void f12sq_fit(double gamma, double lo, double hi, double eps,
                        Tensor<double>& pcoeff, Tensor<double>& pexpnt, bool prnt) {
        Tensor<double> coeff1,expnt1, coeff2, expnt2;
        slater_fit(gamma, lo, hi, eps, coeff1, expnt1, prnt);
        slater_fit(2.0*gamma, lo, hi, eps, coeff2, expnt2, prnt);

        // check exponents are the same
        MADNESS_CHECK((expnt1-expnt2).normf()/expnt1.size()<1.e-12);

        // add exponential terms
        auto coeff=coeff2-2.0*coeff1;
        pcoeff=Tensor<double>(coeff1.size()+1);
        pcoeff(Slice(1,-1,1))=coeff(_);
        pexpnt=Tensor<double>(expnt1.size()+1);
        pexpnt(Slice(1,-1,1))=expnt1(_);

        // add constant term
        pcoeff(0l)=1.0;
        pexpnt(0l)=1.e-10;
    }


    void static bsh_fit_ndim(int ndim, double mu, double lo, double hi, double eps,
			Tensor<double>& pcoeff, Tensor<double>& pexpnt, bool prnt) {

		if (mu > 0) {
			// Restrict hi according to the exponential decay
			double r = -log(4*constants::pi*0.01*eps);
			r = -log(r * 4*constants::pi*0.01*eps);
			if (hi > r) hi = r;
		}


		// Determine range of quadrature by estimating when
		// kernel drops to exp(-100)

		double slo, shi;
		if (mu > 0) {
			slo = -0.5*log(4.0*100.0/(mu*mu));
			slo = -0.5*log(4.0*(slo*ndim - 2.0*slo + 100.0)/(mu*mu));
		}
		else {
			slo = log(eps/hi) - 1.0;
		}
		shi = 0.5*log(100.0/(lo*lo));

		// Resolution required for quadrature over s
		double h = 1.0/(0.2-.50*log10(eps)); // was 0.5 was 0.47

		// Truncate the number of binary digits in h's mantissa
		// so that rounding does not occur when performing
		// manipulations to determine the quadrature points and
		// to limit the number of distinct values in case of
		// multiple precisions being used at the same time.
		h = floor(64.0*h)/64.0;


		// Round shi/lo up/down to an integral multiple of quadrature points
		shi = ceil(shi/h)*h;
		slo = floor(slo/h)*h;

		long npt = long((shi-slo)/h+0.5);

		if (prnt)
			std::cout << "bsh: mu " << mu << " lo " << lo << " hi " << hi
			<< " eps " << eps << " slo " << slo << " shi " << shi
			<< " npt " << npt << " h " << h << std::endl;


		// Compute expansion pruning small coeffs and large exponents
		Tensor<double> coeff(npt), expnt(npt);
		int nnpt=0;
		for (int i=0; i<npt; ++i) {
			double s = slo + h*(npt-i);	// i+1
			double c = exp(-0.25*mu*mu*exp(-2.0*s)+(ndim-2)*s)*0.5/pow(constants::pi,0.5*ndim);
			double p = exp(2.0*s);
			c = c*h;
			if (c*exp(-p*lo*lo) > eps) {
				coeff(nnpt) = c;
				expnt(nnpt) = p;
				++nnpt;
			}
		}
		npt = nnpt;
#if ONE_TERM
		npt=1;
		double s=1.0;
		coeff[0]=1.0;
		expnt[0] = exp(2.0*s);
		coeff=coeff(Slice(0,0));
		expnt=expnt(Slice(0,0));
		print("only one term in gfit",s,coeff[0],expnt[0]);

#endif


		// Prune small exponents from Coulomb fit.  Evaluate a gaussian at
		// the range midpoint, and replace it there with the next most
		// diffuse gaussian.  Then examine the resulting error at the two
		// end points ... if this error is less than the desired
		// precision, can discard the diffuse gaussian.

		if (mu == 0.0) {
			double mid = lo + (hi-lo)*0.5;
			long i;
			for (i=npt-1; i>0; --i) {
				double cnew = coeff[i]*exp(-(expnt[i]-expnt[i-1])*mid*mid);
				double errlo = coeff[i]*exp(-expnt[i]*lo*lo) -
						cnew*exp(-expnt[i-1]*lo*lo);
				double errhi = coeff[i]*exp(-expnt[i]*hi*hi) -
						cnew*exp(-expnt[i-1]*hi*hi);
				if (std::max(std::abs(errlo),std::abs(errhi)) > 0.03*eps) break;
				npt--;
				coeff[i-1] = coeff[i-1] + cnew;
			}
		}

		// Shrink array to correct size
		coeff = coeff(Slice(0,npt-1));
		expnt = expnt(Slice(0,npt-1));


		if (prnt) {
			for (int i=0; i<npt; ++i)
				std::cout << i << " " << coeff[i] << " " << expnt[i] << std::endl;

			long npt = 300;
			std::cout << "       x         value" << std::endl;
			std::cout << "  ------------  ---------------------" << std::endl;
			double step = exp(log(hi/lo)/(npt+1));
			for (int i=0; i<=npt; ++i) {
				double r = lo*(pow(step,i+0.5));
				double test = 0.0;
				for (int j=0; j<coeff.dim(0); ++j)
					test += coeff[j]*exp(-r*r*expnt[j]);
				printf("  %.6e %20.10e\n",r, test);
			}
		}

		pcoeff = coeff;
		pexpnt = expnt;
	}

	// Returns in q[0..4] int(r^2(n+1)*exp(-alpha*r^2),r=0..R) n=-1,0,1,2
	static void gaussian_spherical_moments(double alpha, double R, Tensor<double>& q) {
		q[0] = -(-0.1e1 + exp(-alpha * R*R)) / alpha / 0.2e1;
		q[1] = (-0.2e1 * R * pow(alpha, 0.3e1 / 0.2e1) + sqrt(constants::pi)
				* erf(R * sqrt(alpha)) * alpha * exp(alpha * R*R))
				* pow(alpha, -0.5e1 / 0.2e1) * exp(-alpha * R*R) / 0.4e1;
		q[2] = -(-0.1e1 + exp(-alpha * R*R) + exp(-alpha * R*R) * alpha * R*R)
				* pow(alpha, -0.2e1) / 0.2e1;
		q[3] = -(-0.3e1 * sqrt(constants::pi) * erf(R * sqrt(alpha)) * pow(alpha, 0.2e1)
				* exp(alpha * R*R) + 0.6e1 * R * pow(alpha, 0.5e1 / 0.2e1)
				+ 0.4e1 * pow(R, 0.3e1) * pow(alpha, 0.7e1 / 0.2e1))
				* pow(alpha, -0.9e1 / 0.2e1) * exp(-alpha * R*R) / 0.8e1;
	}

	// Returns in q[0..4] int(r^2(n+1)*exp(-mu*r)/(4*constants::pi*r),r=0..R) n=-1,0,1,2
	static void bsh_spherical_moments(double mu, double R, Tensor<double>& q) {
		if (mu == 0.0) {
			q[0] = R / constants::pi  / 0.4e1;
			q[1] = pow(R, 0.2e1) / constants::pi / 0.8e1;
			q[2] = pow(R, 0.3e1) / constants::pi / 0.12e2;
			q[3] = pow(R, 0.4e1) / constants::pi / 0.16e2;
		}
		else  {
			q[0] = (exp(mu * R) - 0.1e1) / mu * exp(-mu * R) / constants::pi / 0.4e1;
			q[1] = -(-exp(mu * R) + 0.1e1 + mu * R) * pow(mu, -0.2e1) / constants::pi
					* exp(-mu * R) / 0.4e1;
			q[2] = -(-0.2e1 * exp(mu * R) + 0.2e1 + 0.2e1 * mu * R + R*R *
					pow(mu, 0.2e1))*pow(mu, -0.3e1) / constants::pi * exp(-mu * R) / 0.4e1;
			q[3] = -(-0.6e1 * exp(mu * R) + 0.6e1 + 0.6e1 * mu * R + 0.3e1 * R*R
					* pow(mu, 0.2e1) + pow(R, 0.3e1) * pow(mu, 0.3e1))
					* pow(mu, -0.4e1) / constants::pi * exp(-mu * R) / 0.4e1;
		}
	}

};


}

#endif // MADNESS_MRA_GFIT_H__INCLUDED

