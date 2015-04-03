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

  $Id$
*/

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES


/*!
  \file examples/correlationfactor.h
  \brief class for regularizing singular potentials in the molecular
  Hamilton operator

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
  /trunk/src/apps/examples/correlationfactor.h>here</a>.



  \par Introduction

  The correlation factors are intended to represent the cusps (nuclear or
  electronic) in the molecular wave function. Their commutator over the
  kinetic energy operator give a potential operator that cancels the
  singular potential

  [T, f12]		= U + 1/r12
  R^-1[T, R]	= U_nuc	- sum_A Z_A/r1A

  The regularized potentials U and U_nuc contain two terms each, a local
  potential, (denoted U2), and a potential that is multiplied with the
  derivative operator \nabla (denoted U1)

  U 			= U1 . (\nabla_1 - \nabla_2) + U2
  U_nuc			= U1 . \nabla + U2

  with

  U2			= (\nabla^2 f12)
  U2			= R^{-1}(\nabla^2 R)

  \vec U1		= (\vec \nabla f12)
  \vec U1		= R^{-1}(\vec \nabla R)

  To construct a nuclear correlation factor write:

  std::shared_ptr<NuclearCorrelationFactor> nuclear_correlation
   = create_nuclear_correlation_factor(world,*calc);

  where calc is an SCF calculation which holds the molecule and the
  nuclear_corrfac parameter name.
*/


#ifndef MADNESS_CHEM_NUCLEARCORRELATIONFACTOR_H__INCLUDED
#define MADNESS_CHEM_NUCLEARCORRELATIONFACTOR_H__INCLUDED


#include <madness/mra/mra.h>
#include <madness/mra/lbdeux.h>
#include <chem/molecule.h>
#include <chem/potentialmanager.h>

namespace madness {
class SCF;

/// ABC for the nuclear correlation factors
class NuclearCorrelationFactor {
public:
	enum corrfactype {None, GaussSlater, LinearSlater, Polynomial,
		Slater, Two};
	typedef std::shared_ptr< FunctionFunctorInterface<double,3> > functorT;

	/// ctor

	/// @param[in]	world	the world
	/// @param[in]	mol molecule with the sites of the nuclei
	NuclearCorrelationFactor(World& world, const Molecule& mol)
		: world(world), vtol(FunctionDefaults<3>::get_thresh()*0.1)
		, molecule(mol) {}

	/// virtual destructor
	virtual ~NuclearCorrelationFactor() {};

	/// initialize the regularized potentials U1 and U2
	void initialize() {

		// construct the potential functions
		// keep tighter threshold for orthogonalization
		for (int axis=0; axis<3; ++axis) {
			functorT U1f=functorT(new U1_functor(this,axis));
			U1_function.push_back(real_factory_3d(world).thresh(vtol)
					.functor(U1f).truncate_on_project());
			U1_function.back().set_thresh(FunctionDefaults<3>::get_thresh());
		}

		// U2 is the term -S"/S - Z/r
		functorT U2f=functorT(new U2_functor(this));
		U2_function=real_factory_3d(world).thresh(vtol)
				.functor(U2f).truncate_on_project();
		U2_function.set_thresh(FunctionDefaults<3>::get_thresh());

		// U3 is the term SA'/SA . SB'/SB
		functorT U3f=functorT(new U3_functor(this));
		real_function_3d tmp=real_factory_3d(world).thresh(vtol)
				.functor(U3f).truncate_on_project();
		tmp.set_thresh(FunctionDefaults<3>::get_thresh());
		U2_function+=tmp;
		U2_function.truncate();
	}

	virtual corrfactype type() const = 0;

	/// apply the regularized potential U_nuc on a given function rhs
	virtual real_function_3d apply_U(const real_function_3d& rhs) const {

		// the purely local part
		real_function_3d result=(U2()*rhs).truncate();

		// the part with the derivative operators
        result.compress();
        for (int axis=0; axis<3; ++axis) {
            real_derivative_3d D = free_space_derivative<double,3>(world, axis);
            const real_function_3d Drhs=D(rhs).truncate();
            result+=U1(axis)*Drhs;
        }

        result.truncate();
        return result;
	}

	/// return the nuclear correlation factor
	virtual real_function_3d function() const {
		functorT Rf=functorT(new R_functor(this,1));
		real_function_3d r=real_factory_3d(world).thresh(vtol)
				.functor(Rf).truncate_on_project();
		return r;
	}

	/// return the square of the nuclear correlation factor
	virtual real_function_3d square() const {
		real_function_3d R2=real_factory_3d(world).thresh(vtol)
				.functor2(R_functor(this,2)).truncate_on_project();
		return R2;
	}

    /// return the square of the nuclear correlation factor multiplied with
	/// the nuclear potential for the specified atom

	/// @return R^2 * Z_A/r_{1A}
    virtual real_function_3d square_times_V(const Atom& atom) const {
        real_function_3d R2=real_factory_3d(world).thresh(vtol)
                .functor2(square_times_V_functor(this,atom)).truncate_on_project();
        return R2;
    }

    /// return the square of the nuclear correlation factor multiplied with
    /// the derivative of the nuclear potential for the specified atom

    /// @return R^2 * \frac{\partial Z_A/r_{1A}}{\partial X_A}
    virtual real_function_3d square_times_V_derivative(const int iatom, const int axis) const {
        real_function_3d R2=real_factory_3d(world).thresh(vtol)
                .functor2(square_times_V_derivative_functor(this,molecule,iatom,axis)).truncate_on_project();
        return R2;
    }

	/// return the inverse nuclear correlation factor
	virtual real_function_3d inverse() const {
		real_function_3d R_inverse=real_factory_3d(world).thresh(vtol)
				.functor2(R_functor(this,-1)).truncate_on_project();
		return R_inverse;
	}

	/// return the U1 term of the correlation function
	virtual real_function_3d U1(const int axis) const {
		return U1_function[axis];
	}

	/// return the U2 term of the correlation function
	virtual real_function_3d U2() const  {return U2_function;}

private:

	/// the world
	World& world;

	/// the threshold for initial projection
	double vtol;

	/// the molecule
	const Molecule& molecule;

	/// the three components of the U1 potential
	std::vector<real_function_3d> U1_function;

	/// the purely local U2 potential, having absorbed the nuclear pot V_nuc
	real_function_3d U2_function;

	/// the correlation factor S wrt a given atom

	/// @param[in]	r	the distance of the req'd coord to the nucleus
	/// @param[in]	Z	the nuclear charge
	/// @return		the nuclear correlation factor S_A(r_1A)
	virtual double S(const double& r, const double& Z) const = 0;

	/// the partial derivative of correlation factor S' wrt a given atom

	/// @param[in]	vr1A	the vector of the req'd coord to the nucleus
	/// @param[in]	Z	the nuclear charge
	/// @return		the gradient of the nuclear correlation factor S'_A(r_1A)
	virtual coord_3d Sp(const coord_3d& vr1A, const double& Z) const = 0;

	/// the regularized potential wrt a given atom

	/// this is:  -S"/S - Z/r
	/// @param[in]	r	the distance of the req'd coord to the nucleus
	/// @param[in]	Z	the nuclear charge
	/// @return 	the Laplacian of the nuclear correlation factor divided
	///				by the correlation factor minus the nuclear potential
	virtual double Spp_div_S(const double& r, const double& Z) const = 0;

public:

	class R_functor : public FunctionFunctorInterface<double,3> {
		const NuclearCorrelationFactor* ncf;
		int exponent;
	public:
		R_functor(const NuclearCorrelationFactor* ncf, const int e=1)
			: ncf(ncf), exponent(e) {}
		double operator()(const coord_3d& xyz) const {
			double result=1.0;
			for (int i=0; i<ncf->molecule.natom(); ++i) {
				const Atom& atom=ncf->molecule.get_atom(i);
				const coord_3d vr1A=xyz-atom.get_coords();
				const double r=vr1A.normf();
				result*=ncf->S(r,atom.q);
			}
			if (exponent==-1) return 1.0/result;
			else if (exponent==2) return result*result;
			else if (exponent==1) return result;
			else {
				return std::pow(result,double(exponent));
			}

		}
		std::vector<coord_3d> special_points() const {
			return ncf->molecule.get_all_coords_vec();
		}
	};

	/// functor for the local part of the U1 potential -- NOTE THE SIGN

	/// U1 = -S'/S
	class U1_functor : public FunctionFunctorInterface<double,3> {

		const NuclearCorrelationFactor* ncf;
		const int axis;

	public:
		U1_functor(const NuclearCorrelationFactor* ncf, const int axis)
			: ncf(ncf), axis(axis) {}

		double operator()(const coord_3d& xyz) const {
			double result=0.0;
			for (int i=0; i<ncf->molecule.natom(); ++i) {
				const Atom& atom=ncf->molecule.get_atom(i);
				const coord_3d vr1A=xyz-atom.get_coords();
				const double r=vr1A.normf();
				result+=(ncf->Sp(vr1A,atom.q)[axis]/ncf->S(r,atom.q));
			}
			return -1.0*result;
		}
		std::vector<coord_3d> special_points() const {
			return ncf->molecule.get_all_coords_vec();
		}
	};

	class U2_functor : public FunctionFunctorInterface<double,3> {
		const NuclearCorrelationFactor* ncf;
	public:
		U2_functor(const NuclearCorrelationFactor* ncf) : ncf(ncf) {}
		double operator()(const coord_3d& xyz) const {
			double result=0.0;
			for (int i=0; i<ncf->molecule.natom(); ++i) {
				const Atom& atom=ncf->molecule.get_atom(i);
				const coord_3d vr1A=xyz-atom.get_coords();
				const double r=vr1A.normf();
				result+=ncf->Spp_div_S(r,atom.q);
			}
			return result;
		}
		std::vector<coord_3d> special_points() const {
			return ncf->molecule.get_all_coords_vec();
		}
	};

	class U3_functor : public FunctionFunctorInterface<double,3> {
		const NuclearCorrelationFactor* ncf;
	public:
		U3_functor(const NuclearCorrelationFactor* ncf) : ncf(ncf) {}
		double operator()(const coord_3d& xyz) const {
			std::vector<coord_3d> all_terms(ncf->molecule.natom());
			for (int i=0; i<ncf->molecule.natom(); ++i) {
				const Atom& atom=ncf->molecule.get_atom(i);
				const coord_3d vr1A=xyz-atom.get_coords();
				const double r=vr1A.normf();
				all_terms[i]=ncf->Sp(vr1A,atom.q)*(1.0/ncf->S(r,atom.q));
			}

			double result=0.0;
			for (int i=0; i<ncf->molecule.natom(); ++i) {
				for (int j=0; j<i; ++j) {
					result+=all_terms[i][0]*all_terms[j][0]
					       +all_terms[i][1]*all_terms[j][1]
					       +all_terms[i][2]*all_terms[j][2];
				}
			}

			return -1.0*result;
		}
		std::vector<coord_3d> special_points() const {
			return ncf->molecule.get_all_coords_vec();
		}
	};

    class square_times_V_functor : public FunctionFunctorInterface<double,3> {
        const NuclearCorrelationFactor* ncf;
        const Atom& thisatom;
    public:
        square_times_V_functor(const NuclearCorrelationFactor* ncf,
                const Atom& atom1) : ncf(ncf), thisatom(atom1) {}
        double operator()(const coord_3d& xyz) const {
            double result=1.0;
            for (int i=0; i<ncf->molecule.natom(); ++i) {
                const Atom& atom=ncf->molecule.get_atom(i);
                const coord_3d vr1A=xyz-atom.get_coords();
                const double r=vr1A.normf();
                result*=ncf->S(r,atom.q);
            }
            const coord_3d vr1A=xyz-thisatom.get_coords();
            const double V=thisatom.atomic_number/(vr1A.normf()+1.e-6);
            return result*result*V;

        }
        std::vector<coord_3d> special_points() const {
            return ncf->molecule.get_all_coords_vec();
        }
    };


    class square_times_V_derivative_functor : public FunctionFunctorInterface<double,3> {
        const NuclearCorrelationFactor* ncf;
        const Molecule& molecule;
        const int iatom;
        const int axis;
    public:
        square_times_V_derivative_functor(const NuclearCorrelationFactor* ncf,
                const Molecule& molecule1, const int atom1, const int axis1)
            : ncf(ncf), molecule(molecule1), iatom(atom1), axis(axis1) {}
        double operator()(const coord_3d& xyz) const {
            double result=1.0;
            for (int i=0; i<ncf->molecule.natom(); ++i) {
                const Atom& atom=ncf->molecule.get_atom(i);
                const coord_3d vr1A=xyz-atom.get_coords();
                const double r=vr1A.normf();
                result*=ncf->S(r,atom.q);
            }
            const double Vprime=molecule.nuclear_attraction_potential_derivative(
                    iatom, axis, xyz[0], xyz[1], xyz[2]);
            return result*result*Vprime;

        }
        std::vector<coord_3d> special_points() const {
            return ncf->molecule.get_all_coords_vec();
        }
    };

};

/// A nuclear correlation factor class

/// The nuclear correlation factor is given by
/// \[f
/// 	R = \prod S_A	; S_A=exp(-Z_A r_{1A}) + ( 1 - exp(-r_{1A}^2) )
/// \]f
class GaussSlater : public NuclearCorrelationFactor {
public:
	/// ctor

	/// @param[in]	world	the world
	/// @param[in]	mol molecule with the sites of the nuclei
	GaussSlater(World& world, const Molecule& mol)
		: NuclearCorrelationFactor(world,mol) {

		if (world.rank()==0) {
			print("constructed nuclear correlation factor of the form");
			print("  R   = Prod_A S_A");
			print("  S_A = exp(-Z_A r_{1A}) + (1 - exp(-Z_A^2*r_{1A}^2))");
			print("which is of Gaussian-Slater type\n");
		}

		initialize();
	}

	corrfactype type() const {return NuclearCorrelationFactor::GaussSlater;}

private:

	/// the nuclear correlation factor
	double S(const double& r, const double& Z) const {
		const double rho=r*Z;
		return exp(-rho)+(1.0-exp(-(rho*rho)));
	}

	/// radial part first derivative of the nuclear correlation factor
	coord_3d Sp(const coord_3d& vr1A, const double& Z) const {

		const double r=sqrt(vr1A[0]*vr1A[0] +
				vr1A[1]*vr1A[1] + vr1A[2]*vr1A[2]);

		const double eA=exp(-Z*r);
		const double gA=exp(-Z*Z*r*r);
		coord_3d term=(2.0*gA*Z*Z*vr1A-Z*eA*n12(vr1A,1.e-8));
		return term;
	}

	/// second derivative of the nuclear correlation factor

	/// -1/2 S"/S - Z/r
	double Spp_div_S(const double& r, const double& Z) const {
		const double rho=Z*r;
    	if (rho<1.e-4) {
    		return Z*Z*(-3.5 - 4.0*rho + 6.0*rho*rho + 12.0*rho*rho*rho);
    	} else {
			const double e=exp(-rho);
			const double g=exp(-rho*rho);
			const double term1=-Z/r*(1.0-g);
			const double term2=-g*Z*Z*(3.0-2.0*Z*Z*r*r) - Z*Z/2.0*e;
			const double S_inv=exp(-rho)+(1.0-exp(-(rho*rho)));
			return (term1+term2)/S_inv;
    	}
	}

};



/// A nuclear correlation factor class

/// The nuclear correlation factor is given by
/// \[f
/// 	R = \prod S_A	; S_A= -Z_A r_{1A} exp(-Z_A r_{1A}) + 1
/// \]f
class LinearSlater : public NuclearCorrelationFactor {
public:
	/// ctor

	/// @param[in]	world	the world
	/// @param[in]	mol molecule with the sites of the nuclei
	LinearSlater(World& world, const Molecule& mol, const double a)
		: NuclearCorrelationFactor(world,mol), a_(1.0) {

		if (a!=0.0) a_=a;

		if (world.rank()==0) {
			print("constructed nuclear correlation factor of the form");
			print("  S_A = -Z_A r_{1A} exp(-Z_A r_{1A}) + 1");
			print("    a = ",a_);
			print("which is of linear Slater type\n");
		}
		initialize();
	}

	corrfactype type() const {return NuclearCorrelationFactor::LinearSlater;}

private:

	/// the length scale parameter a
	double a_;

	double a_param() const {return 1.0;}

	/// the nuclear correlation factor
	double S(const double& r, const double& Z) const {
		const double rho=r*Z;
		const double b=a_param();
		return (-rho)*exp(-b*rho)+1.0;
	}

	/// radial part first derivative of the nuclear correlation factor
	coord_3d Sp(const coord_3d& vr1A, const double& Z) const {

		const double b=a_param();
		const double r=sqrt(vr1A[0]*vr1A[0] +
				vr1A[1]*vr1A[1] + vr1A[2]*vr1A[2]);

		const double ebrz=exp(-b*r*Z);
		const coord_3d term=Z*ebrz*(b*Z*(vr1A) - n12(vr1A));
		return term;
	}

	/// second derivative of the nuclear correlation factor

	/// -1/2 S"/S - Z/r
	double Spp_div_S(const double& r, const double& Z) const {

		const double b=a_param();
		const double rho=Z*r;
    	if (rho<1.e-4) {
    		const double O0=1.0- 3.0* b;
    		const double O1=Z - 4.0*b*Z + 3.0*b*b*Z;
    		const double O2=Z*Z - 5.0*b*Z*Z + 6.5*b*b*Z*Z - 5.0/3.0*b*b*b*Z*Z;
    		return Z*Z*(O0 + O1*r + O2*r*r);

    	} else {
			const double ebrz=exp(-b*rho);
			const double num=Z* (ebrz - 1.0 + 0.5*ebrz*rho* (2.0 + b*(b*rho-4.0)));
			const double denom=r*(rho*ebrz-1.0);
			return -num/denom;
    	}
	}
};



/// A nuclear correlation factor class
class Slater : public NuclearCorrelationFactor {
public:
	/// ctor

	/// @param[in]	world	the world
	/// @param[in]	mol molecule with the sites of the nuclei
	Slater(World& world, const Molecule& mol, const double a)
		: NuclearCorrelationFactor(world,mol), a_(1.5) {

		if (a!=0.0) a_=a;

		if (world.rank()==0) {
			print("\nconstructed nuclear correlation factor of the form");
			print("  S_A = 1/(a-1) exp(-a Z_A r_{1A}) + 1");
			print("    a = ",a_);
			print("which is of Slater type\n");
		}
		initialize();
	}

	corrfactype type() const {return NuclearCorrelationFactor::Slater;}

private:

	/// the length scale parameter
	double a_;

	double a_param() const {return a_;}

    /// the nuclear correlation factor
    double S(const double& r, const double& Z) const {
    	const double a=a_param();
    	return 1.0+1.0/(a-1.0) * exp(-a*Z*r);
    }

    /// radial part first derivative of the nuclear correlation factor
    coord_3d Sp(const coord_3d& vr1A, const double& Z) const {
    	const double a=a_param();
		const double r=vr1A.normf();
    	return -(a*exp(-a*Z*r)*Z)/(a-1.0)*n12(vr1A);
    }

    /// second derivative of the nuclear correlation factor
    double Spp_div_S(const double& r, const double& Z) const {
    	const double a=a_param();

    	if (r*Z<1.e-4) {
    		const double O0=1.0-(1.5*a);
    		const double O1=(a-1.0)*(a-1.0)*Z;
    		const double O2=(1.0/12.0 * (a-1.0)*(12.0+a*(5*a-18.0)))*Z*Z;
    		return Z*Z*(O0 + O1*r + O2*r*r);

    	} else {
    		const double earz=exp(-a*r*Z);
    		const double num=Z*(-earz + a*earz - (a-1.0) - 0.5*a*a*r*Z*earz);
    		const double denom=(r*earz + (a-1.0) * r);
    		return num/denom;
    	}
    }

};

/// A nuclear correlation factor class

/// should reduce to quartic for N=4
/// @tparam	N	the exponent of the polynomial
template<std::size_t N>
class Polynomial : public NuclearCorrelationFactor {
public:
	/// ctor

	/// @param[in]	world	the world
	/// @param[in]	mol molecule with the sites of the nuclei
	Polynomial(World& world, const Molecule& mol, const double a)
		: NuclearCorrelationFactor(world,mol) {

		/// length scale parameter a, default chosen that linear terms in U2 vanish
		a_=(2. + (-2. + sqrt(-1. + N))*N)/(-2. + N);

		if (a!=0.0) a_=a;

		if (world.rank()==0) {
			print("constructed nuclear correlation factor of the form");
			print("  R   = Prod_A S_A");
			print("  S_A = 1 + a (r/b -1)^N  if  r<b, with  b= (N*a)/((1+a) Z)");
			print("      = 1                 else ");
			print("which is of polynomial type with exponent N = ",N);
		}
		initialize();
	}

	corrfactype type() const {return NuclearCorrelationFactor::Polynomial;}

private:

	/// length scale parameter a, default chosen that linear terms in U2 vanish
	double a_;

	double a_param() const {return a_;}

	/// the cutoff
	static double b_param(const double& a) {return N*a/(1.0+a);}

    /// the nuclear correlation factor
    double S(const double& r, const double& Z) const {

    	const double rho=r*Z;
    	const double a=Polynomial<N>::a_param();
    	const double b=Polynomial<N>::b_param(a);

    	if (rho<b) {
    		const double arg=-1.0 + rho/b;
    		return 1.0 + power<N>(-1.0) * a*power<N>(arg);
    	} else {
    		return 1.0;
    	}

    }

    /// radial part first derivative of the nuclear correlation factor
    coord_3d Sp(const coord_3d& vr1A, const double& Z) const {

		const double r=vr1A.normf();
    	const double rho=r*Z;
    	const double a=Polynomial<N>::a_param();
    	const double b=Polynomial<N>::b_param(a);

    	if (rho<b) {
    		return power<N>(-1.)*(1.+a)* Z* power<N-1>(-1.+rho/b)*n12(vr1A);
    	}
    	return coord_3d(0.0);
    }

    /// second derivative of the nuclear correlation factor

    /// -1/2 S"/S - Z/r
    double Spp_div_S(const double& r, const double& Z) const {

    	const double rho=r*Z;
    	const double a=Polynomial<N>::a_param();
    	const double b=Polynomial<N>::b_param(a);

    	if (rho<1.e-6) {
    		const double ap1=1.0+a;
    		const double c0=((3. *(1. + a) - (3. + a) * N))/(2.* a*N);
    		const double c1=((2.* ap1*ap1 - ap1* (3. + a)*N + N*N)*Z)/(a*a*N*N);
    		const double c2=((30.*ap1*ap1*ap1- ap1*ap1* (55 + 18* a)*N +
    				   30 *ap1 *N*N + (-5 + a* (8 + a)) *N*N*N)* Z*Z)/(12 *a*a*a*N*N*N);
    		return Z*Z*(c0 + c1*r + c2*r*r);

    	} else if (rho<b) {

    		const double num=Z* (2 + (power<N>(-1)* a* power<N>(-1 + rho/b)
    			    * (-2 *a*N*N + (1 + a) *N* (1 + a *(-3 + N) + N)* rho +
    			      2 *(1 + a)*(1+a)* rho*rho))/power<2>(a* N - (1 + a)*rho));

    		const double denom=2.* (r + power<N>(-1) *a* r* power<N>(-1 + rho/b));
    		return -num/denom;

    	} else {
    		return -Z*Z/rho;
    	}
    }

};


class PseudoNuclearCorrelationFactor : public NuclearCorrelationFactor {

public:
	/// ctor

	/// @param[in]	world	the world
	/// @param[in]	mol molecule with the sites of the nuclei
	PseudoNuclearCorrelationFactor(World& world, const Molecule& mol,
			const std::shared_ptr<PotentialManager> pot, const double fac)
		: NuclearCorrelationFactor(world,mol), potentialmanager(pot), fac(fac) {

		if (world.rank()==0) {
			print("constructed nuclear correlation factor of the form");
			print("    R   = ",fac);
			print("which means it's (nearly) a conventional calculation\n");
		}
		initialize();

		// add the missing -Z/r part to U2!
	}

	corrfactype type() const {return None;}

	/// return the U2 term of the correlation function

	/// overloading to avoid inconsistent state of U2, which needs the
	/// nuclear potential
	real_function_3d U2() const {

//		if (not U2_function.is_initialized()) {
			MADNESS_ASSERT(potentialmanager->vnuclear().is_initialized());
//		}
		return potentialmanager->vnuclear();
	}

	/// apply the regularized potential U_nuc on a given function rhs

	/// overload the base class method for efficiency
	real_function_3d apply_U(const real_function_3d& rhs) const {
        return (U2()*rhs).truncate();
	}


private:

	/// underlying potential (=molecule)
    std::shared_ptr<PotentialManager> potentialmanager;

    /// the factor of the correlation factor: R=fac;
	const double fac;

    /// the nuclear correlation factor
    double S(const double& r, const double& Z) const {
    	return fac;
    }

    /// radial part first derivative of the nuclear correlation factor
    coord_3d Sp(const coord_3d& vr1A, const double& Z) const {
    	return coord_3d(0.0);
    }

    /// second derivative of the nuclear correlation factor

    /// this is missing the -Z/r part!
    double Spp_div_S(const double& r, const double& Z) const {
    	return 0.0;
    }
};

std::shared_ptr<NuclearCorrelationFactor>
create_nuclear_correlation_factor(World& world, const SCF& calc);

/// a class holding the electronic correlation factor for R12 theory
class CorrelationFactor {

    World& world;
    double _gamma;      ///< the correlation factor exp(-gamma r12)
    double dcut;		///< the cutoff for the 1/r potential
    double lo;			///< smallest length scale to be resolved

public:

    /// ctor, use negative gamma for linear correlation factor r12
    CorrelationFactor(World& world) : world(world), _gamma(-1.0), dcut(1.e-10),
    	lo(1.e-10) {
    }

    /// ctor, use negative gamma for linear correlation factor r12
    CorrelationFactor(World& world, const double& gamma, const double dcut,
    		const Molecule& molecule) : world(world), _gamma(gamma), dcut(dcut) {
        lo = molecule.smallest_length_scale();
        if (world.rank()==0) {

        	if (gamma>0.0) print("constructed correlation factor with gamma=",gamma);
            else if (gamma==0.0) print("constructed linear correlation factor");
        }
    }

    /// copy ctor
    CorrelationFactor(const CorrelationFactor& other) : world(other.world) {
    	_gamma=other._gamma;
    	dcut=other.dcut;
    	lo=other.lo;
    }

    /// assignment; assume other's world is this world
    CorrelationFactor& operator=(const CorrelationFactor& other) {
    	_gamma=other._gamma;
    	dcut=other.dcut;
    	lo=other.lo;
    	return *this;
    }

    /// return the exponent of this correlation factor
    double gamma() const {return _gamma;}

    /// return the value of the correlation factor
    double operator()(const coord_6d& r) const {
        const double rr=r12(r);
        if (_gamma>0.0) return (1.0-exp(-_gamma*rr))/(2.0*_gamma);
        return 0.5*rr;
    }

    /// apply Kutzelnigg's regularized potential to an orbital product
    real_function_6d apply_U(const real_function_3d& phi_i, const real_function_3d& phi_j,
    		const double eps) const {
//        const double bsh_thresh=FunctionDefaults<6>::get_thresh*0.1;
        const double bsh_thresh=1.e-7;

        real_function_6d result=real_factory_6d(world);

        real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2*eps), lo,bsh_thresh);
        op_mod.modified()=true;

        for (int axis=0; axis<3; ++axis) {
            if (world.rank()==0) print("working on axis",axis);
            real_derivative_3d D = free_space_derivative<double,3>(world, axis);
            const real_function_3d Di=(D(phi_i)).truncate();
            const real_function_3d Dj=(D(phi_j)).truncate();

            const real_function_6d u=U1(axis);

            real_function_6d tmp1=CompositeFactory<double,6,3>(world)
                        .g12(u).particle1(copy(Di)).particle2(copy(phi_j));
            tmp1.fill_tree(op_mod).truncate();
            real_function_6d tmp2=CompositeFactory<double,6,3>(world)
                        .g12(u).particle1(copy(phi_i)).particle2(copy(Dj));
            tmp2.fill_tree(op_mod).truncate();
            if (world.rank()==0) print("done with fill_tree");

            result=result+(tmp1-tmp2).truncate();
            tmp1.clear();
            tmp2.clear();
            world.gop.fence();
            result.truncate().reduce_rank();

            if (world.rank()==0) printf("done with multiplication with U at ime %.1f\n",wall_time());
            result.print_size("result");
        }

//        load_balance(result,true);

        // include the purely local potential that (partially) cancels 1/r12
        if (_gamma>0.0) {
            real_function_6d fg3=real_factory_6d(world).functor2(fg_(_gamma,dcut)).is_on_demand();
            real_function_6d mul=CompositeFactory<double,6,3>(world)
                                .g12(fg3).particle1(copy(phi_i)).particle2(copy(phi_j));
            mul.fill_tree(op_mod).truncate();
            mul.print_size("mul");

            result=(result+mul).truncate().reduce_rank();
        }
        result.print_size("U * |ij>");
        return result;
    }

    /// return the U1 term of the correlation function
    real_function_6d U1(const int axis) const {
        const real_function_6d u1=real_factory_6d(world)
        		.functor2(U(_gamma,axis,dcut)).is_on_demand();
        return u1;
    }

    /// return the U1 term of the correlation function
    real_function_6d U2() const {
    	if (world.rank()==0) print("U2 for the electronic correlation factor");
    	if (world.rank()==0) print("is expensive -- do you really need it??");
    	MADNESS_EXCEPTION("U2() not implemented, since it might be expensive",1);
    	return real_factory_6d(world);
    }

    /// return the correlation factor as on-demand function
    real_function_6d f() const {
//        real_function_6d tmp=real_factory_6d(world).functor2(*this).is_on_demand();
    	double thresh=FunctionDefaults<3>::get_thresh();
        real_function_6d tmp=TwoElectronFactory(world)
        		.dcut(dcut).gamma(_gamma).f12().thresh(thresh);
        return tmp;
    }

    /// return f^2 as on-demand function
    real_function_6d f2() const {
        real_function_6d tmp=real_factory_6d(world).functor2(f2_(_gamma)).is_on_demand();
        return tmp;
    }

    /// return fg+sth as on-demand function
    real_function_6d fg() const {
        real_function_6d tmp=real_factory_6d(world).functor2(fg_(_gamma,dcut)).is_on_demand();
        return tmp;
    }

    /// return f/r as on-demand function
    real_function_6d f_over_r() const {
        real_function_6d tmp=real_factory_6d(world).functor2(f_over_r_(_gamma,dcut)).is_on_demand();
        return tmp;
    }

    /// return (\nabla f)^2 as on-demand functions
    real_function_6d nablaf2() const {
        real_function_6d tmp=real_factory_6d(world).functor2(nablaf2_(_gamma)).is_on_demand();
        return tmp;
    }

private:
    /// functor for the local potential (1-f12)/r12 + sth (doubly connected term of the commutator)
    struct fg_ {
        double gamma;
        double dcut;
        fg_(double gamma, double dcut) : gamma(gamma), dcut(dcut) {
        	MADNESS_ASSERT(gamma>0.0);
        }
        double operator()(const coord_6d& r) const {
            const double rr=r12(r);
            const double e=exp(-gamma*rr);
            return (1.0-e)*u(rr,dcut) + 0.5*gamma*e;
        }
    };

    /// functor for the local potential (1-f12)/r12
    struct f_over_r_ {
        double gamma;
        double dcut;
        f_over_r_(double gamma, double dcut) : gamma(gamma), dcut(dcut) {
        	MADNESS_ASSERT(gamma>0.0);
        }
        double operator()(const coord_6d& r) const {
            const double rr=r12(r);
            const double e=exp(-gamma*rr);
            return (1.0-e)*u(rr,dcut)/(2.0*gamma);
        }
    };

    /// functor for the local part of the regularized potential: f12/r12*(r1-r2)(D1-D2)
    struct U {
        double gamma;
        int axis;
        double dcut;
        U(double gamma, int axis, double dcut) : gamma(gamma), axis(axis),
        	dcut(dcut) {
            MADNESS_ASSERT(axis>=0 and axis<3);
        }
        double operator()(const coord_6d& r) const {
        	const double rr=r12(r);
        	const coord_3d vr12=vec(r[0]-r[3],r[1]-r[4],r[2]-r[5]);
        	const coord_3d N=n12(vr12);
        	if (gamma>0.0) return -0.5*exp(-gamma*rr)*N[axis];
        	MADNESS_EXCEPTION("no gamma in electronic corrfac::U1",1);
//        	const double rr=r12(r);
//            const double g12=u(rr,dcut);
//            double a=0.5;
//            if (gamma>0.0) a=0.5*exp(-gamma*rr);
//            return -a*x12(r,axis) * g12;
        }
    };

    /// functor for the local potential (1-f12)^2
    struct f2_ {
        double gamma;
        f2_(double gamma) : gamma(gamma) {MADNESS_ASSERT(gamma>0.0);}
        double operator()(const coord_6d& r) const {
            const double rr=r12(r);
            const double e=exp(-gamma*rr);
            const double f=(1.0-e)/(2.0*gamma);
            return f*f;
        }
    };

    /// functor for the local potential (\nabla f)^2
    struct nablaf2_ {
        double gamma;
        nablaf2_(double gamma) : gamma(gamma) {
        	MADNESS_ASSERT(gamma>0.0);
        	MADNESS_ASSERT(gamma==1.0);	
        }
        double operator()(const coord_6d& r) const {
            const double rr=r12(r);
            const double f=exp(-2.0*gamma*rr)/(4.0*gamma*gamma);
            return f;
        }
    };

    /// Smoothed 1/r potential (c is the smoothing distance)
    static double u(double r, double c) {
        r = r/c;
        double r2 = r*r, pot;
        if (r > 6.5){
            pot = 1.0/r;
        } else if (r > 1e-2) {
            pot = erf(r)/r + exp(-r2)*0.56418958354775630;
        } else{
            pot = 1.6925687506432689-r2*(0.94031597257959381-r2*(0.39493270848342941-0.12089776790309064*r2));
        }
        return pot/c;
    }

    static double r12(const coord_6d& r) {
        const double x12=r[0]-r[3];
        const double y12=r[1]-r[4];
        const double z12=r[2]-r[5];
        const double r12=sqrt(x12*x12 + y12*y12 + z12*z12);
        return r12;
    }
    static double x12(const coord_6d& r, const int axis) {
        return r[axis]-r[axis+3];
    }


};

/// a class holding the electronic correlation factor for R12 theory
class CorrelationFactor2 {

    World& world;
    double _gamma;      ///< the correlation factor exp(-gamma r12)
	typedef std::shared_ptr< FunctionFunctorInterface<double,6> > functorT;

public:

    double dcut;		///< the cutoff for the 1/r potential
    double lo;			///< smallest length scale to be resolved
    double vtol;		///< initial projection threshold


    /// ctor, use negative gamma for linear correlation factor r12
    CorrelationFactor2(World& world) : world(world), _gamma(0.5), dcut(1.e-10),
    	lo(1.e-10), vtol(FunctionDefaults<3>::get_thresh()*0.1) {
    	MADNESS_ASSERT(_gamma==0.5);
    }

    /// return the exponent of this correlation factor
    double gamma() const {return _gamma;}

    real_function_6d function() const {
    	functorT R=functorT(new R_functor(_gamma,1));
    	return real_factory_6d(world).functor(R).is_on_demand();
    }

    real_function_6d square() const {
    	functorT R2=functorT(new R_functor(_gamma,2));
    	return real_factory_6d(world).functor(R2).is_on_demand();
    }

    real_function_6d inverse() const {
    	functorT R=functorT(new R_functor(_gamma,-1));
    	return real_factory_6d(world).functor(R).is_on_demand();
    }

    /// return the U1 term of the correlation function
    real_function_6d U1(const int axis) const {
		functorT U1f=functorT(new U1_functor(_gamma,axis));
    	return real_factory_6d(world).functor(U1f).is_on_demand();
    }

    /// return the U2 term of the correlation function
    real_function_6d U2() const {
    	functorT U2f=functorT(new U2_functor(_gamma));
    	return real_factory_6d(world).functor(U2f).is_on_demand();
    }

    /// apply Kutzelnigg's regularized potential to an orbital product
    real_function_6d apply_U(const real_function_6d& psi, const double eps) const {
    	const double bsh_thresh=1.e-7;

    	real_function_6d result=real_factory_6d(world);

    	real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2*eps), lo,bsh_thresh);
    	op_mod.modified()=true;

    	for (int axis=0; axis<3; ++axis) {
    		if (world.rank()==0) print("working on axis",axis);
    		real_derivative_6d D1 = free_space_derivative<double,6>(world, axis);
    		real_derivative_6d D2 = free_space_derivative<double,6>(world, axis+3);
    		const real_function_6d Drhs1=D1(psi).truncate();
    		const real_function_6d Drhs2=D2(psi).truncate();

    		const real_function_6d u1=U1(axis);

    		real_function_6d tmp1=CompositeFactory<double,6,3>(world)
                        		 .g12(u1).ket(copy(Drhs1));
    		tmp1.fill_tree(op_mod).truncate();

    		real_function_6d tmp2=CompositeFactory<double,6,3>(world)
                        		 .g12(u1).ket(copy(Drhs2));
    		tmp2.fill_tree(op_mod).truncate();
    		if (world.rank()==0) print("done with fill_tree");

    		result=result+(tmp1-tmp2).truncate();
    		tmp1.clear();
    		tmp2.clear();
    		world.gop.fence();
    		result.truncate().reduce_rank();

    		if (world.rank()==0)
    			printf("done with multiplication with U at ime %.1f\n",wall_time());
    		result.print_size("result");
    	}

    	real_function_6d u2=U2();
    	real_function_6d r2=CompositeFactory<double,6,3>(world).ket(copy(psi))
     								.g12(u2);
    	r2.fill_tree(op_mod);
    	result=(result+r2).truncate();
    	return result;
    }


private:

    /// functor for the correlation factor R
    class R_functor : public FunctionFunctorInterface<double,6> {
    	double gamma;
    	int exponent;

    public:
    	R_functor(double gamma, int e=1) : gamma(gamma), exponent(e) {
    		MADNESS_ASSERT(gamma==0.5);
    	}

        // only valid for gamma=1
        double operator()(const coord_6d& r) const {
            const double rr=r12(r);
            double val=(1.0-0.5*exp(-gamma*rr));
            if (exponent==1) return val;
            else if (exponent==2) return val*val;
            else if (exponent==-1) return 1.0/val;
            else {
            	MADNESS_EXCEPTION("fancy exponent in correlationfactor2",1);
            }
        }
    };

    /// functor for the U2 local potential
    class U2_functor : public FunctionFunctorInterface<double,6> {
    	double gamma;

    public:
    	U2_functor(double gamma) : gamma(gamma) {
    		MADNESS_ASSERT(gamma==0.5);
    	}

        // only valid for gamma=1
        double operator()(const coord_6d& r) const {
        	const double rr=r12(r);
        	// Taylor expansion for small r
        	if (rr<1.e-4) {	// valid for gamma==0.5, otherwise singular
        		return 5./4.0 - rr + (35.0* rr*rr)/48.0 - (101.0*rr*rr*rr)/192.0;
        	}
        	const double egr=exp(-gamma*rr);
        	return -(-8.*egr + 8.0 + rr*egr)/(4.0 *rr*egr - 8 *rr);
        }
    };

    /// functor for the U1 = -\frac{\vec\nabla_1 f_{12}}{f_{12}}  potential

    /// the potential is given by
    /// U1 = -\frac{\vec\nabla_1 f_{12}}{f_{12}}
    ///    =  \frac{e^{-r12/2}{4-2e^{-r12/2}} \vec n12
    /// the derivative operators are not included
    class U1_functor : public FunctionFunctorInterface<double,6> {
        double gamma;
        int axis;

    public:
        U1_functor(double gamma, int axis) : gamma(gamma), axis(axis) {
        	MADNESS_ASSERT(gamma==0.5);
        	MADNESS_ASSERT(axis<3);
        }

        double operator()(const coord_6d& r) const {
        	const double rr=r12(r);
        	const coord_3d vr12=vec(r[0]-r[3],r[1]-r[4],r[2]-r[5]);
        	const coord_3d N=n12(vr12);
        	// Taylor expansion for small r
        	double val;
        	if (rr<1.e-4) {	// valid for gamma==0.5, otherwise singular
        		val = 0.5 - 0.5*rr + 0.125*(3.*rr*rr) - (13.* rr*rr*rr)/48.0;
        	} else {
            	const double egr=exp(-gamma*rr);
        		val=egr/(4.0-2.0*egr);
        	}
        	// NOTE the sign
        	return -val*N[axis];
        }
    };

    /// helper function
    static double r12(const coord_6d& r) {
        const double x12=r[0]-r[3];
        const double y12=r[1]-r[4];
        const double z12=r[2]-r[5];
        const double r12=sqrt(x12*x12 + y12*y12 + z12*z12);
        return r12;
    }

};


}

#endif /* NUCLEARCORRELATIONFACTOR_H_ */


