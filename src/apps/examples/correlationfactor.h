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

#define WORLD_INSTANTIATE_STATIC_TEMPLATES


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

*/


#ifndef NUCLEARCORRELATIONFACTOR_H_
#define NUCLEARCORRELATIONFACTOR_H_


#include <mra/mra.h>
#include <mra/lbdeux.h>
#include <moldft/molecule.h>
#include <moldft/potentialmanager.h>

using namespace madness;

namespace madness {


/// ABC for the nuclear correlation factors
class NuclearCorrelationFactor {
public:
	enum corrfactype {None, GaussSlater, Two};
	/// ctor

	/// @param[in]	world	the world
	/// @param[in]	molecule	molecule with the sites of the nuclei
	NuclearCorrelationFactor(World& world, const Calculation& calc)
		: world(world) {}

	virtual corrfactype type() const = 0;

	/// apply the regularized potential U_nuc on a given function rhs
	virtual real_function_3d apply_U(const real_function_3d& rhs) const = 0;

	/// return the nuclear correlation factor
	virtual real_function_3d function() const = 0;

	/// return the square of the nuclear correlation factor
	virtual real_function_3d square() const = 0;

	/// return the inverse nuclear correlation factor
	virtual real_function_3d inverse() const = 0;

	/// return the U1 term of the correlation function
	virtual real_function_3d U1(const int axis) const = 0;

	/// return the U2 term of the correlation function
	virtual real_function_3d U2() const = 0;

	/// the world
	World& world;

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
	/// @param[in]	molecule	molecule with the sites of the nuclei
	GaussSlater(World& world, const Calculation& calc)
		: NuclearCorrelationFactor(world,calc), molecule(calc.molecule)
		, vtol(FunctionDefaults<3>::get_thresh()*0.1) {

		if (world.rank()==0) {
			print("constructed nuclear correlation factor of the form");
			print("  R   = Prod_A S_A");
			print("  S_A = exp(-Z_A r_{1A}) + (1 - exp(-Z_A^2*r_{1A}^2))");
			print("which is of Gaussian-Slater type\n");
		}


		// construct the potential functions
		// keep tighter threshold for orthogonalization
		for (int axis=0; axis<3; ++axis) {
			U1_function.push_back(real_factory_3d(world).thresh(vtol)
					.functor2(U1_functor(molecule,axis)).truncate_on_project());
			U1_function.back().set_thresh(FunctionDefaults<3>::get_thresh());
		}

		U2_function=real_factory_3d(world).thresh(vtol)
				.functor2(U2_functor(molecule)).truncate_on_project();
		U2_function.set_thresh(FunctionDefaults<3>::get_thresh());
		real_function_3d tmp=real_factory_3d(world).thresh(vtol)
				.functor2(U3_functor(molecule)).truncate_on_project();
		tmp.set_thresh(FunctionDefaults<3>::get_thresh());
		U2_function+=tmp;
		U2_function.truncate();
	}

	corrfactype type() const {return NuclearCorrelationFactor::GaussSlater;}

	/// apply the regularized potential U_nuc on a given function rhs
	real_function_3d apply_U(const real_function_3d& rhs) const {

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

	real_function_3d function() const {
		real_function_3d r=real_factory_3d(world).thresh(vtol)
				.functor2(R(1,molecule)).truncate_on_project();
		return r;
	}

	real_function_3d square() const {
		real_function_3d R2=real_factory_3d(world).thresh(vtol)
				.functor2(R(2,molecule)).truncate_on_project();
		return R2;
	}

	real_function_3d inverse() const {
		real_function_3d R_inverse=real_factory_3d(world).thresh(vtol)
				.functor2(R(-1,molecule)).truncate_on_project();
		return R_inverse;
	}

	/// return the U1 term of the correlation function
	real_function_3d U1(const int axis) const {return U1_function[axis];}

	/// return the U2 term of the correlation function
	real_function_3d U2() const {return U2_function;}

private:

//	static double norm() {return sqrt(M_PI)*0.125;}
	static double norm() {return 1.0;}

	/// underlying molecule
	const Molecule& molecule;

	/// the components of the U1 potential
	std::vector<real_function_3d> U1_function;

	/// the purely local U2 potential, having absorbed the nuclear pot V_nuc
	real_function_3d U2_function;

	/// the threshold for initial projection
	double vtol;

	struct R {
		int exponent;	/// return R^{exponent}
		const Molecule& molecule;

		R(const int e, const Molecule& m) : exponent(e), molecule(m) {}
        double operator()(const coord_3d& r) const {
			const double x=r[0], y=r[1], z=r[2];
			double result=1.0;
			for (int i=0; i<molecule.natom(); ++i) {
				const Atom& atom=molecule.get_atom(i);
				const double& xA=atom.x;
				const double& yA=atom.y;
				const double& zA=atom.z;
				const double rr=sqrt((x-xA)*(x-xA)+(y-yA)*(y-yA)+(z-zA)*(z-zA));
				const double rho=atom.q*rr;
				const double rho2=rho*rho;

				result*=norm()*exp(-rho)+(1.0-exp(-rho2));
			}
			if (exponent==-1) return 1.0/result;
			else if (exponent==2) return result*result;
			else if (exponent==1) return result;
			else {
				return std::pow(result,double(exponent));
			}
		}
	};

	/// functor for the local part of the singly connected part [T,\rho]
    struct U1_functor {

    	const Molecule& molecule;
    	const int axis;

    	U1_functor(const Molecule& molecule, const int axis)
    		: molecule(molecule), axis(axis) {}

        double operator()(const coord_3d& r) const {
        	double result=0.0;
        	for (int i=0; i<molecule.natom(); ++i) {
    			const Atom& atom1=molecule.get_atom(i);
    			const double& ZA=atom1.q;
    			coord_3d vr1A=r-atom1.get_coords();
    			const double r1A=sqrt(vr1A[0]*vr1A[0] +
    					vr1A[1]*vr1A[1] + vr1A[2]*vr1A[2]);

				const double eA=norm()*exp(-ZA*r1A);
    			const double gA=exp(-ZA*ZA*r1A*r1A);
    			const double SA_inv=1.0/(1.0-gA+eA);
    			coord_3d termA=SA_inv*(2.0*gA*ZA*ZA*vr1A-ZA*eA*n12(vr1A,1.e-8));
    			result+=termA[axis];
        	}
			return -1.0*result;
        }
    };

    /// functor for the term -1/2 S"/S - V_nuc
    struct U2_functor {

    	const Molecule& molecule;
    	U2_functor(const Molecule& molecule) : molecule(molecule) {}

        double operator()(const coord_3d& r) const {
        	double result=0.0;
        	for (int i=0; i<molecule.natom(); ++i) {
    			const Atom& atom1=molecule.get_atom(i);
    			const double& Z=atom1.q;
    			const coord_3d vr1A=r-atom1.get_coords();
    			const double r1A=sqrt(vr1A[0]*vr1A[0] +
    					vr1A[1]*vr1A[1] + vr1A[2]*vr1A[2]);

    			const double e=norm()*exp(-Z*r1A);
    			const double g=exp(-Z*Z*r1A*r1A);

    			const double S_inv=1.0/(1.0-g+e);
    			const double term1=-Z/r1A*(1.0-g);
    			const double term2=-g*Z*Z*(3.0-2.0*Z*Z*r1A*r1A) - Z*Z/2.0*e;
    			result+=S_inv*(term1+term2);
    		}
        	return result;
        }
    };

    /// functor for the term -1/2 S'_A/S_A . S'_B/S_B
    struct U3_functor {

    	const Molecule& molecule;
    	U3_functor(const Molecule& molecule) : molecule(molecule) {}

        double operator()(const coord_3d& r) const {
        	double result=0.0;
        	for (int i=0; i<molecule.natom(); ++i) {
    			const Atom& atom1=molecule.get_atom(i);
    			const double& ZA=atom1.q;
    			const coord_3d vr1A=r-atom1.get_coords();
    			const double r1A=sqrt(vr1A[0]*vr1A[0] +
    					vr1A[1]*vr1A[1] + vr1A[2]*vr1A[2]);

    			const double eA=norm()*exp(-ZA*r1A);
    			const double gA=exp(-ZA*ZA*r1A*r1A);
    			const double SA_inv=1.0/(1.0-gA+eA);
    			const coord_3d termA=SA_inv*2.0*ZA*ZA*gA*vr1A - SA_inv*ZA*eA*n12(vr1A);

            	for (int j=0; j<i; ++j) {
        			const Atom& atom2=molecule.get_atom(j);
        			const double& ZB=atom2.q;
        			const coord_3d vr1B=r-atom2.get_coords();
        			const double r1B=sqrt(vr1B[0]*vr1B[0] +
        					vr1B[1]*vr1B[1] + vr1B[2]*vr1B[2]);

        			const double eB=norm()*exp(-ZB*r1B);
        			const double gB=exp(-ZB*ZB*r1B*r1B);
        			const double SB_inv=1.0/(1.0-gB+eB);

        			const coord_3d termB=SB_inv*(2.0*ZB*ZB*gB*vr1B-ZB*eB*n12(vr1B));
        			result+=(termA[0]*termB[0]+termA[1]*termB[1]+termA[2]*termB[2]);
        		}
    		}
        	return -1.0*result;
        }
    };

};


class PseudoNuclearCorrelationFactor : public NuclearCorrelationFactor {

	/// underlying potential (=molecule)
    std::shared_ptr<PotentialManager> potentialmanager;

public:
	/// ctor

	/// @param[in]	world	the world
	/// @param[in]	molecule	molecule with the sites of the nuclei
	PseudoNuclearCorrelationFactor(World& world, const Calculation& calc)
		: NuclearCorrelationFactor(world,calc),
		  potentialmanager(calc.potentialmanager) {

		if (world.rank()==0) {
			print("constructed nuclear correlation factor of the form 1");
			print("    R   = 1");
			print("which means it's a conventional calculation\n");
		}
	}

	corrfactype type() const {return None;}

	/// apply the regularized potential U_nuc on a given function rhs
	real_function_3d apply_U(const real_function_3d& rhs) const {
		const real_function_3d& Vnuc=potentialmanager->vnuclear();
        return (Vnuc*rhs).truncate();
	}

	real_function_3d function() const {
		real_function_3d R1=real_factory_3d(world)
				.functor2(R()).truncate_on_project();
		return R1;
	}

	/// the inverse of this is this
	real_function_3d inverse() const {
		real_function_3d R_inverse=real_factory_3d(world)
				.functor2(R()).truncate_on_project();
		return R_inverse;
	}

	/// the inverse of this is this
	real_function_3d square() const {
		real_function_3d R2=real_factory_3d(world)
				.functor2(R()).truncate_on_project();
		return R2;
	}

	/// return the U1 term of the correlation function

	/// the term is the local part of  R^{-1} (\vec\nabla R) \nabla,
	/// without the trailing grad operator
	/// this is an empty function [T,1]=0
	real_function_3d U1(const int axis) const {
		return real_factory_3d(world);
	}

	/// return the U2 term of the correlation function

	/// the term is R^{-1} (\nabla^2 R)
	/// this is the nuclear potential [T,1] + V_nuc = 0 + V_nuc
	real_function_3d U2() const {
		const real_function_3d Vnuc=potentialmanager->vnuclear();
		return Vnuc;
	}


private:

	struct R {
		double operator()(const coord_3d& r) const {return 1.0;}
	};

};



class TwoNuclearCorrelationFactor : public NuclearCorrelationFactor {

	/// underlying potential (=molecule)
    std::shared_ptr<PotentialManager> potentialmanager;

public:
	/// ctor

	/// @param[in]	world	the world
	/// @param[in]	molecule	molecule with the sites of the nuclei
	TwoNuclearCorrelationFactor(World& world, const Calculation& calc)
		: NuclearCorrelationFactor(world,calc),
		  potentialmanager(calc.potentialmanager) {

		if (world.rank()==0) {
			print("constructed nuclear correlation factor of the form 1");
			print("    R   = 2");
			print("which means it's nearly a conventional calculation\n");
		}
	}

	corrfactype type() const {return Two;}

	/// apply the regularized potential U_nuc on a given function rhs
	real_function_3d apply_U(const real_function_3d& rhs) const {
		const real_function_3d& Vnuc=potentialmanager->vnuclear();
        return (Vnuc*rhs).truncate();
	}

	real_function_3d function() const {
		real_function_3d R1=real_factory_3d(world)
				.functor2(R(2.0)).truncate_on_project();
		return R1;
	}

	/// the inverse of this is this
	real_function_3d inverse() const {
		real_function_3d R_inverse=real_factory_3d(world)
				.functor2(R(0.5)).truncate_on_project();
		return R_inverse;
	}

	/// the inverse of this is this
	real_function_3d square() const {
		real_function_3d R2=real_factory_3d(world)
				.functor2(R(4.0)).truncate_on_project();
		return R2;
	}

	/// return the U1 term of the correlation function

	/// the term is the local part of  R^{-1} (\vec\nabla R) \nabla,
	/// without the trailing grad operator
	/// this is an empty function [T,1]=0
	real_function_3d U1(const int axis) const {
		return real_factory_3d(world);
	}

	/// return the U2 term of the correlation function

	/// the term is R^{-1} (\nabla^2 R)
	/// this is the nuclear potential [T,1] + V_nuc = 0 + V_nuc
	real_function_3d U2() const {
		const real_function_3d Vnuc=potentialmanager->vnuclear();
		return Vnuc;
	}


private:

	struct R {
		double fac;
		R(const double f) : fac(f){}
		double operator()(const coord_3d& r) const {return fac;}
	};

};


/// a class holding the correlation factor for R12 theory
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
        real_function_6d tmp=real_factory_6d(world).functor2(*this).is_on_demand();
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
        	MADNESS_ASSERT(gamma=1.0);	// I don't think this is right
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


}

#endif /* NUCLEARCORRELATIONFACTOR_H_ */


