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
*/

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES

/*!
  \file apps/chem/correlationfactor.h
  \brief class for regularizing singular potentials in the molecular
  Hamilton operator

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


#ifndef MADNESS_CHEM_NUCLEARCORRELATIONFACTOR_H_
#define MADNESS_CHEM_NUCLEARCORRELATIONFACTOR_H_

#include <madness/mra/mra.h>
#include <madness/mra/lbdeux.h>
#include <chem/molecule.h>
#include <chem/potentialmanager.h>
#include <chem/atomutil.h>

using namespace madness;

namespace madness {

/// ABC for the nuclear correlation factors
class NuclearCorrelationFactor {
public:
	enum corrfactype {None, GradientalGaussSlater, GaussSlater, LinearSlater,
	    Polynomial, Slater, Two};
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
	    R_functor r(this,2);
		real_function_3d R2=real_factory_3d(world).thresh(vtol)
				.functor(r).truncate_on_project();
		return R2;
	}

    /// return the square of the nuclear correlation factor multiplied with
    /// the derivative of the nuclear potential for the specified atom

    /// @return R^2 * \frac{\partial Z_A/r_{1A}}{\partial X_A}
    virtual real_function_3d square_times_V_derivative(const int iatom, const int axis) const {
        square_times_V_derivative_functor func(this,molecule,iatom,axis);
        real_function_3d R2=real_factory_3d(world).thresh(vtol)
                .functor(func).truncate_on_project();
        return R2;
    }

	/// return the inverse nuclear correlation factor
	virtual real_function_3d inverse() const {
	    R_functor r(this,-1);
		real_function_3d R_inverse=real_factory_3d(world).thresh(vtol)
				.functor(r).truncate_on_project();
		return R_inverse;
	}

	/// return the U1 term of the correlation function
	virtual const real_function_3d U1(const int axis) const {
		return U1_function[axis];
	}

    /// return the U1 functions in a vector
	std::vector<real_function_3d> U1vec() const {
	    std::vector<real_function_3d> uvec(3);
	    uvec[0]=U1_function[0];
	    uvec[1]=U1_function[1];
	    uvec[2]=U1_function[2];
	    return uvec;
	}

	/// return the U2 term of the correlation function
	virtual const real_function_3d U2() const  {return U2_function;}

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

	/// the partial derivative of correlation factor S' wrt the cartesian coordinates

	/// @param[in]	vr1A	the vector of the req'd coord to the nucleus
	/// @param[in]	Z	the nuclear charge
	/// @return		the gradient of the nuclear correlation factor S'_A(r_1A)
	virtual coord_3d Sp(const coord_3d& vr1A, const double& Z) const = 0;

	/// the regularized potential wrt a given atom wrt the cartesian coordinate

	/// S" is the Cartesian Laplacian applied on the NCF. Note the difference
	/// to Srr_div_S, which is the second derivative wrt the distance rho.
	/// this is:  -S"/S - Z/r
	/// @param[in]	r	the distance of the req'd coord to the nucleus
	/// @param[in]	Z	the nuclear charge
	/// @return 	the Laplacian of the nuclear correlation factor divided
	///				by the correlation factor minus the nuclear potential
	virtual double Spp_div_S(const double& r, const double& Z) const = 0;

public:

	/// first derivative of the NCF with respect to the relative distance rho
	/// \f[
	///     \frac{\partial S(\rho)}{\partial \rho} \frac{1}{S(\rho)}
	/// \f]
	/// where the distance of the electron to the nucleus A is given by
	/// \f[
	///    \rho = |\vec r - \vec R_A |
	/// \f]
	virtual double Sr_div_S(const double& r, const double& Z) const = 0;

	/// second derivative of the NCF with respect to the relative distance rho
    /// \f[
    ///     \frac{\partial^2 S(\rho)}{\partial \rho^2} \frac{1}{S(\rho)}
    /// \f]
    /// where the distance of the electron to the nucleus A is given by
    /// \f[
    ///    \rho = |\vec r - \vec R_A |
    /// \f]
	virtual double Srr_div_S(const double& r, const double& Z) const = 0;

    /// third derivative of the NCF with respect to the relative distance rho
    /// \f[
    ///     \frac{\partial^3 S(\rho)}{\partial \rho^3} \frac{1}{S(\rho)}
    /// \f]
    /// where the distance of the electron to the nucleus A is given by
    /// \f[
    ///    \rho = |\vec r - \vec R_A |
    /// \f]
	virtual double Srrr_div_S(const double& r, const double& Z) const = 0;

    /// derivative of the U2 potential wrt nuclear coordinate X (spherical part)

    /// need to reimplement this for all derived classes due to the
    /// range for r -> 0, where the singular terms cancel. With
    /// \f[
    ///   \rho = \left| \vec r- \vec R_A \right|
    /// \f]
    /// returns the term in the parenthesis without the the derivative of rho
    /// \f[
    /// \frac{\partial U_2}{\partial X_A} = \frac{\partial \rho}{\partial X}
    ///           \left(-\frac{1}{2}\frac{S''' S - S'' S'}{S^2} + \frac{1}{\rho^2}\frac{S'}{S}
    ///           - \frac{1}{\rho} \frac{S''S - S'^2}{S^2} + \frac{Z_A}{\rho^2}\right)
    /// \f]
    virtual double U2X_spherical(const double& r, const double& Z, const double& rcut) const {
        if (world.rank()==0) {
            print("you can't compute the Hessian matrix");
            print("U2X_spherical is not implemented for the nuclear correlation factor");
        }
        MADNESS_EXCEPTION("do more implementation work",1);
    }

public:

	/// smoothed unit vector for the computation of the U1 potential

    /// note the identity for exchanging nuclear and electronic coordinates
    /// (there is a sign change, unlike for the smoothed potential)
    /// \f[
    ///     \vec n  = \frac{\partial \rho}{\partial x} = -\frac{\partial \rho}{\partial X}
    /// \f]
    /// \f[
    /// \vec n = \left\{\frac{x \mathrm{erf}\left(\frac{r}{s}\right)}{r},
    ///          \frac{y \mathrm{erf}\left(\frac{r}{s}\right)}{r},
    ///          \frac{z \mathrm{erf}\left(\frac{r}{s}\right)}{r}\right\}
    /// \f]
	coord_3d smoothed_unitvec(const coord_3d& xyz, double smoothing=0.0) const {
#if 0

        if (smoothing==0.0) smoothing=molecule.get_eprec();
        // TODO:need to test this
        // reduce the smoothing for the unitvector
        //if (not (this->type()==None or this->type()==Two)) smoothing=sqrt(smoothing);
        smoothing=sqrt(smoothing);
        const double r=xyz.normf();
        const double rs=r/smoothing;
        if (r<1.e-4) {
            const double sqrtpi=sqrt(constants::pi);
            double erfrs_div_r=2.0/(smoothing*sqrtpi)-2.0/3.0*rs*rs/(sqrtpi*smoothing);
            return erfrs_div_r*xyz;
        } else if (r<6.0) {
            return erf(rs)/r*xyz;
        } else {
            return 1.0/r*xyz;
        }



#else
        if (smoothing==0.0) smoothing=molecule.get_eprec();
        // TODO:need to test this
        // reduce the smoothing for the unitvector
        //if (not (this->type()==None or this->type()==Two)) smoothing=sqrt(smoothing);
        const double r=xyz.normf();
        const double cutoff=smoothing;
        if (r>cutoff) {
            return 1.0/r*xyz;
        } else {
            const double xi=r/cutoff;
            const double xi2=xi*xi;
            const double xi3=xi*xi*xi;
//            const double nu21=0.5+1./32.*(45.*xi - 50.*xi3 + 21.*xi*xi*xi*xi*xi);
            const double nu22=0.5 + 1./64.*(105* xi - 175 *xi3 + 147* xi2*xi3 - 45* xi3*xi3*xi);
//            const double nu40=0.5 + 1./128.*(225 *xi - 350 *xi3 + 189*xi2*xi3);
            const double kk=2.*nu22-1.0;
            return kk/r*xyz;
        }

#endif
	}

	/// derivative of smoothed unit vector wrt the *electronic* coordinate

	/// note the sign change for exchanging nuclear and electronic coordinates
	/// \f[
	///     \frac{\partial \vec n}{\partial x}  = -\frac{\partial \vec n}{\partial X}
	/// \f]
	/// the derivative wrt x is given by
	/// \f[
	/// \frac{\partial\vec n}{\partial x} =
	/// \left\{\frac{\left(r^2-x^2\right) \mathrm{erf}\left(\frac{r}{s}\right)}{r^3}
	///    +\frac{2 x^2 e^{-\frac{r^2}{s^2}}}{\sqrt{\pi } r^2 s},
	///  x y \left(\frac{2 e^{-\frac{r^2}{s^2}}}{\sqrt{\pi } r^2 s}
	///    -\frac{\mathrm{erf}\left(\frac{r}{s}\right)}{r^3}\right),
	///  x z \left(\frac{2 e^{-\frac{r^2}{s^2}}}{\sqrt{\pi } r^2 s}
	///    -\frac{\mathrm{erf}\left(\frac{r}{s}\right)}{r^3}\right)\right\}
	/// \f]
	coord_3d dsmoothed_unitvec(const coord_3d& xyz, const int axis,
            double smoothing=0.0) const {

	    const double r=xyz.normf();
        coord_3d result;
        if (smoothing==0.0) smoothing=molecule.get_eprec();

#if 1
        // TODO:need to test this
        // reduce the smoothing for the unitvector
        //if (not (this->type()==None or this->type()==Two)) smoothing=sqrt(smoothing);
        smoothing=sqrt(smoothing);

        const double rs=r/smoothing;
        const static double sqrtpi=sqrt(constants::pi);
        const double sqrtpis3=sqrtpi*smoothing*smoothing*smoothing;

        if (r<1.e-4) {
            // series expansion
            double p=-4.0/(3.0*sqrtpis3) + 4.0*rs*rs/(5.0*sqrtpis3);

            double erfrs_div_r=2.0/(smoothing*sqrtpi)-2.0/3.0*rs*rs/(sqrtpi*smoothing);
            result=xyz*xyz[axis]*p;
            result[axis]+=erfrs_div_r;

        } else if (r<6.0) {
            const double erfrs_div_r=erf(rs)/r;
            const double term1=2.0*exp(-rs*rs)/(sqrtpi*r*r*smoothing);
            result=xyz*xyz[axis]*(term1-erfrs_div_r/(r*r));
            result[axis]+=erfrs_div_r;
#else
        if (r<smoothing) {
            double r2=r*r;
            double s2=smoothing*smoothing;
            double s7=s2*s2*s2*smoothing;
            double x2=xyz[axis]*xyz[axis];

            double fac_offdiag=-(((135. *r2*r2 - 294.* r2 *s2
                    + 175.*s2*s2))/(16.* s7));
            double fac_diag=-((45.* r2*r2*r2 - 147.* r2*r2* s2
                    + 175.* r2*s2*s2 - 105.* s2*s2*s2 + 270.* r2*r2* x2
                    - 588.* r2* s2* x2 + 350.*s2* s2 *x2)/(32.* s7));

            result[0]=fac_offdiag*xyz[0]*xyz[axis];
            result[1]=fac_offdiag*xyz[1]*xyz[axis];
            result[2]=fac_offdiag*xyz[2]*xyz[axis];
            result[axis]=fac_diag;

#endif
        } else {
            result=xyz*(-xyz[axis]/(r*r*r));
            result[axis]+=1/r;
        }
        return result;
    }


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
				const double& Z=atom.q;
//				result-=(ncf->Sp(vr1A,Z)[axis]/ncf->S(r,Z));
				result-=ncf->Sr_div_S(r,Z)*ncf->smoothed_unitvec(vr1A)[axis];
			}
			return result;
		}
		std::vector<coord_3d> special_points() const {
			return ncf->molecule.get_all_coords_vec();
		}
	};

    /// U1 functor for a specific atom

	/// NOTE THE SIGN !!
	/// this is
	/// \f[
	///  -\frac{\partial \rho}{\partial X_A}\frac{\partial S}{\partial \rho}\frac{1}{S}
	/// \f]
    class U1_atomic_functor : public FunctionFunctorInterface<double,3> {

        const NuclearCorrelationFactor* ncf;
        const int iatom;
        const int axis;

    public:
        U1_atomic_functor(const NuclearCorrelationFactor* ncf, const int atom,
                const int axis) : ncf(ncf), iatom(atom), axis(axis) {}

        double operator()(const coord_3d& xyz) const {
            const Atom& atom=ncf->molecule.get_atom(iatom);
            const coord_3d vr1A=xyz-atom.get_coords();
            const double r=vr1A.normf();
            const double& Z=atom.q;
            return ncf->Sr_div_S(r,Z)*ncf->smoothed_unitvec(vr1A)[axis];
        }

        std::vector<coord_3d> special_points() const {
            std::vector< madness::Vector<double,3> > c(1);
            const Atom& atom=ncf->molecule.get_atom(iatom);
            c[0][0]=atom.x;
            c[0][1]=atom.y;
            c[0][2]=atom.z;
            return c;
        }
    };


    /// functor for a local U1 dot U1 potential

    /// the unit vector dotted with itself vanishes, so what's left is
    /// \f[
    ///  U1\dot U1 = \frac{\left(S^r\right)^2}{S^2}
    /// \f]
    /// with positive sign!
    class U1_dot_U1_functor : public FunctionFunctorInterface<double,3> {

        const NuclearCorrelationFactor* ncf;

    public:
        U1_dot_U1_functor(const NuclearCorrelationFactor* ncf) : ncf(ncf) {}

        double operator()(const coord_3d& xyz) const {
            double result=0.0;
            for (int i=0; i<ncf->molecule.natom(); ++i) {
                const Atom& atom=ncf->molecule.get_atom(i);
                const coord_3d vr1A=xyz-atom.get_coords();
                const double r=vr1A.normf();
                const double& Z=atom.q;
                const double tmp=ncf->Sr_div_S(r,Z);
                result+=tmp*tmp;
            }
            return result;
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
//				all_terms[i]=ncf->Sp(vr1A,atom.q)*(1.0/ncf->S(r,atom.q));
				all_terms[i]=ncf->Sr_div_S(r,atom.q)*ncf->smoothed_unitvec(vr1A);
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

    /// U2 functor for a specific atom
    class U2_atomic_functor : public FunctionFunctorInterface<double,3> {

        const NuclearCorrelationFactor* ncf;
        const int iatom;

    public:
        U2_atomic_functor(const NuclearCorrelationFactor* ncf, const int atom)
            : ncf(ncf), iatom(atom) {}

        double operator()(const coord_3d& xyz) const {
            const Atom& atom=ncf->molecule.get_atom(iatom);
            const coord_3d vr1A=xyz-atom.get_coords();
            const double r=vr1A.normf();
            return ncf->Spp_div_S(r,atom.q);
        }

        std::vector<coord_3d> special_points() const {
            std::vector< madness::Vector<double,3> > c(1);
            const Atom& atom=ncf->molecule.get_atom(iatom);
            c[0][0]=atom.x;
            c[0][1]=atom.y;
            c[0][2]=atom.z;
            return c;
        }
    };

    /// U3 functor for a specific atom
    class U3_atomic_functor : public FunctionFunctorInterface<double,3> {

        const NuclearCorrelationFactor* ncf;
        const int iatom;

    public:
        U3_atomic_functor(const NuclearCorrelationFactor* ncf, const int atom)
            : ncf(ncf), iatom(atom) {}

        double operator()(const coord_3d& xyz) const {
            const Atom& atomA=ncf->molecule.get_atom(iatom);
            const coord_3d vr1A=xyz-atomA.get_coords();
            const double rA=vr1A.normf();
            const coord_3d nA=ncf->smoothed_unitvec(vr1A);
            double Sr_div_SA=ncf->Sr_div_S(rA,atomA.q);

            double result=0.0;
            // sum over B
            for (int i=0; i<ncf->molecule.natom(); ++i) {
                if (i==iatom) continue; // restricted sum

                const Atom& atomB=ncf->molecule.get_atom(i);
                const coord_3d vr1B=xyz-atomB.get_coords();
                const double rB=vr1B.normf();
                const coord_3d nB=ncf->smoothed_unitvec(vr1B);
                double Sr_div_SB=ncf->Sr_div_S(rB,atomB.q);

                double dot=nA[0]*nB[0] + nA[1]*nB[1] + nA[2]*nB[2];
                result+=Sr_div_SB*Sr_div_SA*dot;
            }
            return -0.5*result;
        }

        std::vector<coord_3d> special_points() const {
            std::vector< madness::Vector<double,3> > c(1);
            const Atom& atom=ncf->molecule.get_atom(iatom);
            c[0][0]=atom.x;
            c[0][1]=atom.y;
            c[0][2]=atom.z;
            return c;
        }
    };

    class square_times_V_functor : public FunctionFunctorInterface<double,3> {
        const NuclearCorrelationFactor* ncf;
        const Molecule& molecule;
        const int iatom;
    public:
        square_times_V_functor(const NuclearCorrelationFactor* ncf,
                const Molecule& mol, const int iatom1)
            : ncf(ncf), molecule(mol), iatom(iatom1) {}
        double operator()(const coord_3d& xyz) const {
            double result=1.0;
            for (int i=0; i<ncf->molecule.natom(); ++i) {
                const Atom& atom=ncf->molecule.get_atom(i);
                const coord_3d vr1A=xyz-atom.get_coords();
                const double r=vr1A.normf();
                result*=ncf->S(r,atom.q);
            }
            const double V=-molecule.atomic_attraction_potential(
                                iatom, xyz[0], xyz[1], xyz[2]);
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

    /// compute the derivative of R wrt the displacement of atom A, coord axis
    class RX_functor : public FunctionFunctorInterface<double,3> {
        const NuclearCorrelationFactor* ncf;
        const Atom& thisatom;
        const int derivativeaxis;   /// direction of the derivative operator
        const int exponent;         /// 1 or 2 -> R^X or R^X R

    public:
        RX_functor(const NuclearCorrelationFactor* ncf, const Atom& atom1,
                const int daxis, const int exponent) : ncf(ncf), thisatom(atom1),
                derivativeaxis(daxis), exponent(exponent) {
            MADNESS_ASSERT((exponent==1) or (exponent==2) or (exponent==-1));
        }

        RX_functor(const NuclearCorrelationFactor* ncf, const int iatom,
                const int daxis, const int exponent) : ncf(ncf),
                thisatom(ncf->molecule.get_atom(iatom)),
                derivativeaxis(daxis), exponent(exponent) {
            MADNESS_ASSERT((exponent==1) or (exponent==2) or (exponent==-1));
        }

        double operator()(const coord_3d& xyz) const {

            // compute the R term
            double result=1.0;
            if ((exponent==1) or (exponent==2)) {
                for (int i=0; i<ncf->molecule.natom(); ++i) {
                    const Atom& atom=ncf->molecule.get_atom(i);
                    const coord_3d vr1A=xyz-atom.get_coords();
                    const double r=vr1A.normf();
                    result*=ncf->S(r,atom.q);
                }
                if (exponent==2) result=result*result;
            }

            // compute the derivative term
            {
                const coord_3d vr1A=xyz-thisatom.get_coords();
                const double r=vr1A.normf();
                const double& Z=thisatom.q;
                const double S1=-ncf->Sr_div_S(r,Z) // note the sign
                        *ncf->smoothed_unitvec(vr1A)[derivativeaxis];
                result*=S1;
            }
            return result;
        }

        std::vector<coord_3d> special_points() const {
            return ncf->molecule.get_all_coords_vec();
        }

    };


    /// compute the derivative of U1 wrt the displacement of atom A, coord axis
    class U1X_functor : public FunctionFunctorInterface<double,3> {
        const NuclearCorrelationFactor* ncf;
        const Atom& thisatom;
        const int U1axis;           /// U1x/U1y/U1z potential?
        const int derivativeaxis;   /// direction of the derivative operator
    public:
        U1X_functor(const NuclearCorrelationFactor* ncf, const Atom& atom1,
                const int U1axis, const int daxis) : ncf(ncf), thisatom(atom1),
                U1axis(U1axis), derivativeaxis(daxis) {
            double lo=1.0/thisatom.q;
            set_length_scale(lo);
        }

        U1X_functor(const NuclearCorrelationFactor* ncf, const int iatom,
                const int U1axis, const int daxis) : ncf(ncf),
                thisatom(ncf->molecule.get_atom(iatom)),
                U1axis(U1axis), derivativeaxis(daxis) {
            double lo=1.0/thisatom.q;
            set_length_scale(lo);
        }

        double operator()(const coord_3d& xyz) const {
            const coord_3d vr1A=xyz-thisatom.get_coords();
            const double r=vr1A.normf();
            const double& Z=thisatom.q;
            const double S1=ncf->Sr_div_S(r,Z);
            const double S2=ncf->Srr_div_S(r,Z);

            // note the sign change smoothed_unitvec due to the
            // change in the derivative variable x: electronic -> nuclear
            const double drhodx=-ncf->smoothed_unitvec(vr1A)[derivativeaxis];
            return drhodx*(S2-S1*S1)*ncf->smoothed_unitvec(vr1A)[U1axis]
                      -S1*(ncf->dsmoothed_unitvec(vr1A,derivativeaxis)[U1axis]);
        }

        std::vector<coord_3d> special_points() const {
            std::vector< madness::Vector<double,3> > c(1);
            c[0][0]=thisatom.x;
            c[0][1]=thisatom.y;
            c[0][2]=thisatom.z;
            return c;
        }

    };


    /// compute the derivative of U2 wrt the displacement of atom A
    class U2X_functor : public FunctionFunctorInterface<double,3> {
        const NuclearCorrelationFactor* ncf;
        const int iatom;
        const int axis;
    public:
        U2X_functor(const NuclearCorrelationFactor* ncf, const int& atom1,
                const int axis) : ncf(ncf), iatom(atom1), axis(axis) {
            const Atom& atom=ncf->molecule.get_atom(iatom);
            double lo=1.0/atom.q;
            set_length_scale(lo);
        }

        double operator()(const coord_3d& xyz) const {
            const Atom& atom=ncf->molecule.get_atom(iatom);
            const coord_3d vr1A=xyz-atom.get_coords();
            const double r=vr1A.normf();
            const double& Z=atom.q;
            const double rcut=ncf->molecule.get_rcut()[iatom];

            // note the sign change due to the change in the derivative
            // variable x: electronic -> nuclear in drho/dx
            const double drhodx=-ncf->smoothed_unitvec(vr1A)[axis];
            return drhodx*ncf->U2X_spherical(r,Z,rcut);
        }

        std::vector<coord_3d> special_points() const {
            std::vector< madness::Vector<double,3> > c(1);
            const Atom& atom=ncf->molecule.get_atom(iatom);
            c[0][0]=atom.x;
            c[0][1]=atom.y;
            c[0][2]=atom.z;
            return c;
        }
    };


    /// compute the derivative of U3 wrt the displacement of atom A, coord axis

    /// \f[
    /// U_3^{X_A} = -\sum_{B\neq A}\left(\frac{\vec S_A'}{S_A}\right)^X\cdot\left(\frac{\vec S_B'}{S_B}\right)
    /// \f]
    /// with
    /// \f[
    /// \left(\frac{\vec S_A'}{S_A}\right)^X =
    ///     \frac{\partial \rho}{\partial X}\left(\frac{S''_A}{S_A}
    ///          -\left(\frac{S'_A}{S_A}\right)^2\right)\vec n_{1A}
    ///     + \left(\frac{S'_A}{S_A}\right)\frac{\partial \vec n_{1A}}{\partial X}
    /// \f]
    class U3X_functor : public FunctionFunctorInterface<double,3> {
        const NuclearCorrelationFactor* ncf;
        const int iatom;
        const int axis;
    public:
        U3X_functor(const NuclearCorrelationFactor* ncf, const int iatom,
                const int axis) : ncf(ncf), iatom(iatom), axis(axis) {}

        double operator()(const coord_3d& xyz) const {
            const Atom& atomA=ncf->molecule.get_atom(iatom);
            const coord_3d vr1A=xyz-atomA.get_coords();
            const double r1A=vr1A.normf();
            const double& ZA=atomA.q;

            double S1A=ncf->Sr_div_S(r1A,ZA);
            double S2A=ncf->Srr_div_S(r1A,ZA);
            double termA=S2A-S1A*S1A;

            // unit vector \vec n_A = \vec r_{1A}/r_{1A}
            const coord_3d nA=ncf->smoothed_unitvec(vr1A);
            // derivative of the unit vector \frac{\partial \vec n_A}{\partial X}
            const coord_3d dnA=ncf->dsmoothed_unitvec(vr1A,axis)*(-1.0);
            // \frac{\partial \rho}{\partial X}
            const double drhodx=-nA[axis];

            double term=0.0;
            for (int jatom=0; jatom<ncf->molecule.natom(); ++jatom) {
                if (iatom==jatom) continue; // restricted sum B \neq A

                const Atom& atomB=ncf->molecule.get_atom(jatom);
                const coord_3d vr1B=xyz-atomB.get_coords();
                const double r1B=vr1B.normf();
                const double& ZB=atomB.q;

                double S1B=ncf->Sr_div_S(r1B,ZB);
                const coord_3d nB=ncf->smoothed_unitvec(vr1B);

                double dot=0.0;     // n_A.n_B
                double ddot=0.0;    // n'_A.n_B
                for (int i=0; i<3; ++i) {
                    ddot+=dnA[i]*nB[i];
                    dot+=nA[i]*nB[i];
                }
                term+=(+drhodx*termA*S1B*dot + S1A*S1B*ddot);

            }

            return term;
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
			print("with eprec ",mol.get_eprec());
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
		coord_3d term=(2.0*gA*Z*Z*vr1A-Z*eA*smoothed_unitvec(vr1A));
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

    double Sr_div_S(const double& r, const double& Z) const {
        const double Zr=r*Z;
        const double eA=exp(-Zr);
        const double gA=exp(-Zr*Zr);
        const double num=Z*(2.0*Zr*gA-eA);
        const double denom=1.0+eA-gA;
        return num/denom;
    }

    double Srr_div_S(const double& r, const double& Z) const {
        const double Zr=r*Z;
        const double eA=exp(-Zr);
        const double gA=exp(-Zr*Zr);
        const double num=Z*Z*(eA+gA*(2.0-4.0*Zr*Zr));
        const double denom=1.0+eA-gA;
        return num/denom;
    }

    double Srrr_div_S(const double& r, const double& Z) const {
        const double Zr=r*Z;
        const double eA=exp(-Zr);
        const double gA=exp(-Zr*Zr);
        const double num=Z*Z*Z*(-eA - 12.0*gA*Zr + 8.0*gA*Zr*Zr*Zr);
        const double denom=1.0+eA-gA;
        return num/denom;

    }

    /// derivative of the U2 potential wrt X (scalar part)

    /// with
    /// \f[
    ///   \rho = \left| \vec r- \vec R_A \right|
    /// \f]
    /// returns the term in the parenthesis without the the derivative of rho
    /// \f[
    /// \frac{\partial U_2}{\partial X_A} = \frac{\partial \rho}{\partial X}
    ///           \left(-\frac{1}{2}\frac{S''' S - S'' S'}{S^2} + \frac{1}{\rho^2}\frac{S'}{S}
    ///           - \frac{1}{\rho} \frac{S''S - S'^2}{S^2} + \frac{Z_A}{\rho^2}\right)
    /// \f]
    double U2X_spherical(const double& r, const double& Z, const double& rcut) const {

        double result=0.0;
        if (r*Z<1.e-4) {
            const double ZZ=Z*Z;
            const double ZZZ=ZZ*Z;
            const double Z4=ZZ*ZZ;
            const double r0=-4.0*ZZZ;
            const double r1=12.0*Z4;
            const double r2=36*Z4*Z;
            const double r3=-67.0/6.0*Z4*ZZ;
            result=(r0 + r*r1 + r*r*r2 + r*r*r*r3);

        } else {
            const double S1=Sr_div_S(r,Z);
            const double S2=Srr_div_S(r,Z);
            const double S3=Srrr_div_S(r,Z);
            const double term1=-0.5*(S3-S1*S2);
            const double term2=(S1+Z)/(r*r);
            const double term3=(S2-S1*S1)/r;
            result=term1+term2-term3;
        }
        return result;
    }


};

/// A nuclear correlation factor class

/// The nuclear correlation factor is given by
/// \[f
///     R = \prod S_A   ; S_A=exp(-Z_A r_{1A}) + ( 1 - exp(-r_{1A}^2) )
/// \]f
class GradientalGaussSlater : public NuclearCorrelationFactor {
public:
    /// ctor

    /// @param[in]  world   the world
    /// @param[in]  mol molecule with the sites of the nuclei
    GradientalGaussSlater(World& world, const Molecule& mol, const double a)
        : NuclearCorrelationFactor(world,mol), a(a) {

        if (world.rank()==0) {
            print("constructed nuclear correlation factor of the form");
            print("  R   = Prod_A S_A");
            print("  S_A = 1/sqrt{Z} exp(-Z_A r_{1A}) + (1 - exp(-a^2*Z_A^2*r_{1A}^2))");
            print("  a   = ",a);
            print("with eprec ",mol.get_eprec());
            print("which is of Gradiental Gaussian-Slater type\n");
        }

        initialize();
    }

    corrfactype type() const {return NuclearCorrelationFactor::GradientalGaussSlater;}

private:

    const double a;

    /// the nuclear correlation factor
    double S(const double& r, const double& Z) const {
        const double rho=r*Z;
        return 1/sqrt(Z) * exp(-rho)+(1.0-exp(-(a*a*rho*rho)));
    }

    /// radial part first derivative of the nuclear correlation factor
    coord_3d Sp(const coord_3d& vr1A, const double& Z) const {

        const double r=sqrt(vr1A[0]*vr1A[0] +
                vr1A[1]*vr1A[1] + vr1A[2]*vr1A[2]);

        const double rho=Z*r;
        const double sqrtz=sqrt(Z);
        const double term=-exp(-rho)*sqrtz + 2.0*a*a*exp(-a*a*rho*rho)*Z*rho;
        return term*smoothed_unitvec(vr1A);
    }

    /// second derivative of the nuclear correlation factor

    /// -1/2 S"/S - Z/r
    double Spp_div_S(const double& r, const double& Z) const {
        const double rho=Z*r;
        const double sqrtz=sqrt(Z);
        if (rho<1.e-4) {
            const double zfivehalf=Z*Z*sqrtz;
            const double a2=a*a;
            const double a4=a2*a2;
            return  -0.5*Z*Z
                    - 3. *a2 * zfivehalf
                    - 4.* a2 *rho* zfivehalf
                    - 2. *a2 * rho*rho*zfivehalf
                    + 5. *a4 *rho*rho*zfivehalf
                    + 3. *a4 *rho*rho*Z*Z*Z
                    -0.5 *a2 *rho*rho*rho*zfivehalf
                    +5.5 *a4 *rho*rho*rho*zfivehalf
                    +7.  *a4 *rho*rho*rho*Z*Z*Z;
        } else {
            const double e=exp(-rho);
            const double g=exp(-a*a*rho*rho);
            const double poly=(2.0-6.0*a*a*rho + 4.0*a*a*a*a*rho*rho*rho);
            const double num=Z*(-2.0 - e*r*sqrtz + g*poly);
            const double denom=2.0*r*(1.0-g+e/sqrtz);
            return num/denom;
        }
    }

    double Sr_div_S(const double& r, const double& Z) const {
        const double rZ=r*Z;
        const double e=exp(-rZ);
        const double g=exp(-a*a*rZ*rZ);
        const double sqrtz=sqrt(Z);
        const double num=-sqrtz*e + 2.0*a*a*g*Z*rZ;
        const double denom=1.0-g+e/sqrtz;
        return num/denom;
    }

    double Srr_div_S(const double& r, const double& Z) const {
        const double rZ=r*Z;
        const double e=exp(-rZ);
        const double g=exp(-a*a*rZ*rZ);
        const double sqrtz=sqrt(Z);
        const double num=e*Z*sqrtz + g*(2.0*a*a - 4.0*power<4>(a)*rZ*rZ)*Z*Z;
        const double denom=1.0-g+e/sqrtz;
        return num/denom;
    }

    double Srrr_div_S(const double& r, const double& Z) const {
        const double rZ=r*Z;
        const double e=exp(-rZ);
        const double g=exp(-a*a*rZ*rZ);
        const double sqrtz=sqrt(Z);
        const double num=e*power<3>(Z) + (12.0*power<4>(a)*g*rZ
                -8.0*power<6>(a)*g*power<3>(rZ))*sqrtz*power<3>(Z);
        const double denom=e+sqrtz-g*sqrtz;
        return -num/denom;
    }

    double U2X_spherical(const double& r, const double& Z, const double& rcut) const {

        double result=0.0;
        if (r*Z<1.e-4) {
            const double sqrtz=sqrt(Z);
            const double Z2=Z*Z;
            const double Z4=Z2*Z2;
            const double Z5=Z4*Z;
            const double Z6=Z5*Z;
            const double Z7=Z6*Z;
            const double a2=a*a;
            const double a4=a2*a2;

            const double r0=-4.* a2* sqrt(Z7);
            const double r1=2.* (-2.* a2* Z*sqrt(Z7)+ 5.* a4* Z*sqrt(Z7) + 3.* a4 *Z5) *r;
            const double r2=1.5 * (-a2* sqrtz*Z5 + 11.* a4* sqrtz*Z5 + 14.*a4* Z6)* r*r;
            const double r3=1./6.* (-a2* sqrtz*Z6 + 66.* a4*sqrtz*Z6 - 84.* a2*a4* sqrtz*Z6 +
                    180. *a4* Z7 - 156.*a2*a4* Z7 - 72.* a2*a4*sqrtz*Z7) *r*r*r;
            result=(r0 + r1 + r2 + r3);

        } else {
            const double S1=Sr_div_S(r,Z);
            const double S2=Srr_div_S(r,Z);
            const double S3=Srrr_div_S(r,Z);
            const double term1=-0.5*(S3-S1*S2);
            const double term2=(S1+Z)/(r*r);
            const double term3=(S2-S1*S1)/r;
            result=term1+term2-term3;
        }
        return result;
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
            print("with eprec ",mol.get_eprec());
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
		const coord_3d term=Z*ebrz*(b*Z*(vr1A) - smoothed_unitvec(vr1A));
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

	double Sr_div_S(const double& r, const double& Z) const {
	    const double& a=a_param();
	    const double earz=exp(-a*r*Z);
	    return Z*earz*(a*r*Z-1.0)/(1.0-r*Z*earz);
	}

    double Srr_div_S(const double& r, const double& Z) const {
        const double& a=a_param();
        const double earz=exp(-a*r*Z);
        return a*Z*Z*earz*(a*r*Z-2.0)/(-1.0+r*Z*earz);
    }

    double Srrr_div_S(const double& r, const double& Z) const {
        const double& a=a_param();
        const double earz=exp(-a*r*Z);
        return a*a*Z*Z*Z*earz*(a*r*Z-3.0)/(1.0-r*Z*earz);
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
            print("with eprec ",mol.get_eprec());
			print("which is of Slater type\n");
		}
		initialize();
	}

	corrfactype type() const {return NuclearCorrelationFactor::Slater;}

private:

	/// the length scale parameter
	double a_;

	double a_param() const {return a_;}

	/// first derivative of the correlation factor wrt (r-R_A)

	/// \f[
	///     Sr_div_S = \frac{1}{S(r)}\frac{\partial S(r)}{\partial r}
	/// \f]
	double Sr_div_S(const double& r, const double& Z) const {
	    const double& a=a_param();
	    return -a*Z/(1.0+(a-1.0)*exp(a*r*Z));
	}

    /// second derivative of the correlation factor wrt (r-R_A)

    /// \f[
    ///     result = \frac{1}{S(r)}\frac{\partial^2 S(r)}{\partial r^2}
    /// \f]
    double Srr_div_S(const double& r, const double& Z) const {
        const double& a=a_param();
        const double aZ=a*Z;
        return aZ*aZ/(1.0+(a-1.0)*exp(r*aZ));
    }

    /// third derivative of the correlation factor wrt (r-R_A)

    /// \f[
    ///    result = \frac{1}{S(r)}\frac{\partial^3 S(r)}{\partial r^3}
    /// \f]
    double Srrr_div_S(const double& r, const double& Z) const {
        const double& a=a_param();
        const double aZ=a*Z;
        return -aZ*aZ*aZ/(1.0+(a-1.0)*exp(r*aZ));
    }

    /// the nuclear correlation factor
    double S(const double& r, const double& Z) const {
    	const double a=a_param();
    	return 1.0+1.0/(a-1.0) * exp(-a*Z*r);
    }

    /// radial part first derivative of the nuclear correlation factor
    coord_3d Sp(const coord_3d& vr1A, const double& Z) const {
    	const double a=a_param();
		const double r=vr1A.normf();
    	return -(a*exp(-a*Z*r)*Z)/(a-1.0)*smoothed_unitvec(vr1A);
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


    /// derivative of the U2 potential wrt X (scalar part)

    /// with
    /// \f[
    ///   \rho = \left| \vec r- \vec R_A \right|
    /// \f]
    /// returns the term in the parenthesis without the the derivative of rho
    /// \f[
    /// \frac{\partial U_2}{\partial X_A} = \frac{\partial \rho}{\partial X}
    ///           \left(-\frac{1}{2}\frac{S''' S - S'' S'}{S^2} + \frac{1}{\rho^2}\frac{S'}{S}
    ///           - \frac{1}{\rho} \frac{S''S - S'^2}{S^2} + \frac{Z_A}{\rho^2}\right)
    /// \f]
    double U2X_spherical(const double& r, const double& Z, const double& rcut) const {
        const double a=a_param();

        double result=0.0;
        if (r*Z<1.e-4) {
            const double ZZ=Z*Z;
            const double ZZZ=ZZ*Z;
            const double a2=a*a;
            const double a4=a2*a2;
            const double r0=ZZZ*(1. - 2.* a + a2);
            const double r1=ZZ*ZZ/6.* (12.0 - 30.* a + 23. *a2 - 5.*a*a2);
            const double r2=1./8.*ZZ*ZZZ* (24. - 72.*a + 74.*a2 - 29.*a2*a + 3.*a4);
            const double r3=1./60.*ZZZ*ZZZ* (240. - 840.*a + 1080.*a2 - 610.*a2*a
                    + 137.*a2*a2 - 7.*a4*a);
            result=(r0 + r*r1 + r*r*r2 + r*r*r*r3);

        } else {
            const double S1=Sr_div_S(r,Z);
            const double S2=Srr_div_S(r,Z);
            const double S3=Srrr_div_S(r,Z);
            const double term1=-0.5*(S3-S1*S2);
            const double term2=(S1+Z)/(r*r);
            const double term3=(S2-S1*S1)/r;
            result=term1+term2-term3;
        }
        return result;
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
			print("with eprec ",mol.get_eprec());
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
    		return power<N>(-1.)*(1.+a)* Z* power<N-1>(-1.+rho/b)*smoothed_unitvec(vr1A);
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

    double Sr_div_S(const double& r, const double& Z) const {
        const double rho=r*Z;
        const double a=Polynomial<N>::a_param();
        const double b=Polynomial<N>::b_param(a);

        if (rho<b) {
            const double negn= power<N>(-1.0);
            const double num=(negn*(1 + a)*Z*power<N-1>(-1 + ((1 + a)*r*Z)/(a*N)));
            const double denom=(1 + negn*a*power<N>(-1 + ((1 + a)*r*Z)/(a*N)));
            return num/denom;
        } else {
            return 0.0;
        }

    }

    double Srr_div_S(const double& r, const double& Z) const {
        const double rho=r*Z;
        const double a=Polynomial<N>::a_param();
        const double b=Polynomial<N>::b_param(a);

        if (rho<b) {
            const double negn= power<N>(-1.0);
            return (negn*power<2>(1 + a)*(-1 + N)*power<2>(Z)*power<N-2>(-1 + ((1 + a)*r*Z)/(a*N)))/
                    (a*N*(1 + negn*a*power<N>(-1 + ((1 + a)*r*Z)/(a*N))));
        } else {
            return 0.0;
        }
    }

    double Srrr_div_S(const double& r, const double& Z) const {
        const double rho=r*Z;
        const double a=Polynomial<N>::a_param();
        const double b=Polynomial<N>::b_param(a);

        if (rho<b) {
            const double negn= power<N>(-1.0);
            return (negn*power<3>(1 + a)*(-2 + N)*(-1 + N)*power<3>(Z)*power<N-3>(-1 + ((1 + a)*r*Z)/(a*N)))/
                    (power<2>(a*N)*(1 + negn*a*power<N>(-1 + ((1 + a)*r*Z)/(a*N))));
        } else {
            return 0.0;
        }
    }

    double U2X_spherical(const double& r, const double& Z, const double& rcut) const {
        const double a=a_param();
        const double aopt=(2. + (-2. + sqrt(-1. + N))*N)/(-2. + N);
        if (fabs(a-aopt)>1.e-10) {
            MADNESS_EXCEPTION("U2X_spherical for polynomial ncf only with aopt",1);
        }

        double result=0.0;
        if (r*Z<1.e-4) {
            const double rn=sqrt(N-1);
            const double r0=0.0;
            const double r1=((2.*(-8. + 9.*rn) + N*(25. + 10.*rn + N))*r*power<4>(Z))/
                    (6.*power<2>(-2 + N)*rn);
            const double r2=((-4*(17 + 9*rn) + N*(92 + 80*rn +
                    N*(-29 - 33*rn + N*(4 + 7*rn + N))))*power<5>(Z))/
                            (8.*power<3>(-2 + N)*(-1 + N)*rn);
            result=(r0 + r*r1 + r*r*r2);

        } else {
            const double S1=Sr_div_S(r,Z);
            const double S2=Srr_div_S(r,Z);
            const double S3=Srrr_div_S(r,Z);
            const double term1=-0.5*(S3-S1*S2);
            const double term2=(S1+Z)/(r*r);
            const double term3=(S2-S1*S1)/r;
            result=term1+term2-term3;
        }
        return result;
    }

};

class PseudoNuclearCorrelationFactor : public NuclearCorrelationFactor {

public:
	/// ctor

	/// @param[in]	world	the world
	/// @param[in]	mol molecule with the sites of the nuclei
	PseudoNuclearCorrelationFactor(World& world, const Molecule& mol,
			const std::shared_ptr<PotentialManager> pot, const double fac)
		: NuclearCorrelationFactor(world,mol), potentialmanager(pot),
		  eprec(mol.get_eprec()), fac(fac) {

		if (world.rank()==0) {
			print("constructed nuclear correlation factor of the form");
			print("    R   = ",fac);
            print("with eprec ",mol.get_eprec());
			print("which means it's (nearly) a conventional calculation\n");
		}
		initialize();

		// add the missing -Z/r part to U2!
	}

	corrfactype type() const {return None;}

	/// return the U2 term of the correlation function

	/// overloading to avoid inconsistent state of U2, which needs the
	/// nuclear potential
	const real_function_3d U2() const {

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
    double eprec;

    /// the factor of the correlation factor: R=fac;
	const double fac;

    double Sr_div_S(const double& r, const double& Z) const {return 0.0;}

    double Srr_div_S(const double& r, const double& Z) const {return 0.0;}

    double Srrr_div_S(const double& r, const double& Z) const {return 0.0;}

    /// the nuclear correlation factor
    double S(const double& r, const double& Z) const {
    	return fac;
    }

    /// radial part first derivative of the nuclear correlation factor
    coord_3d Sp(const coord_3d& vr1A, const double& Z) const {
    	return coord_3d(0.0);
    }

    /// second derivative of the nuclear correlation factor
    double Spp_div_S(const double& r, const double& Z) const {
        double rcut= 1.0 / smoothing_parameter(Z, eprec);
        return - Z * smoothed_potential(r*rcut)*rcut;
    }

    double U2X_spherical(const double& r, const double& Z, const double& rcut) const {
        // factor -1 from the definition of the dsmoothed_potential as -1/r^2
        return -Z*dsmoothed_potential(r * rcut) * (rcut * rcut);
    }

};

class SCF;

std::shared_ptr<NuclearCorrelationFactor>
create_nuclear_correlation_factor(World& world, const SCF& calc);

}
#endif /* NUCLEARCORRELATIONFACTOR_H_ */
