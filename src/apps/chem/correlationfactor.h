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
	    Polynomial, Slater, Slater2, Slater3, Slater4, Slater5, poly4erfc, Two};
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
	void initialize(const double vtol1) {

	    // set threshold for projections
	    vtol=vtol1;

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

//		initialize();
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

//        initialize();
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
//		initialize();
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
		eprec_=mol.get_eprec();

		if (world.rank()==0) {
			print("\nconstructed nuclear correlation factor of the form");
			print("  S_A = 1/(a-1) exp(-a Z_A r_{1A}) + 1");
			print("    a = ",a_);
            print("with eprec ",eprec_);
			print("which is of Slater type\n");
		}
//		initialize();
	}

	corrfactype type() const {return NuclearCorrelationFactor::Slater;}

private:

	/// the length scale parameter
	double a_;
	double eprec_;

	double a_param() const {return a_;}
    double eprec_param() const {return eprec_;}

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
    	const double eprec=eprec_param();
    	return 1.0+1.0/(a-1.0) * exp(-a*Z*r);

//    	return 1.0 + 0.5/(a-1.0) *
//    	        (exp(-a*r*Z + 0.5*a*a*Z*Z*eprec) * erfc((-r+a*eprec*Z)/sqrt(2*eprec))
//    	         + exp(-a*r*Z + 0.5*a*Z*(4.0*r+a*eprec*Z)) * erfc((r+a*eprec*Z)/sqrt(2*eprec)));
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

/// This ncf make the first 3 moments of the inverse ncf vanish
/// int (1-1/S) r^(2+n) dr =0    for n=0,1,2
/// The first  moment of the orbitals is
/// int \phi * (1-1/S) r^(2) dr = -0.000767821
class Slater2 : public NuclearCorrelationFactor {
public:
    /// ctor

    /// @param[in]  world   the world
    /// @param[in]  mol molecule with the sites of the nuclei
    Slater2(World& world, const Molecule& mol, const double aa)
        : NuclearCorrelationFactor(world,mol), a(1.5) {

        if (aa!=0.0) a=aa;
        a=0.4;
        eprec_=mol.get_eprec();

        if (world.rank()==0) {
            print("\nconstructed nuclear correlation factor of the form");
            print("  S_A = ");
            print("    a = ",a);
            print("with eprec ",eprec_);
            print("which is of Slater2 type\n");
        }
        a0=0.00023607587702526534367;
        a1=-0.028860765047694765039;
        a2=1.3589075252693079170;
        a3=-3.8302828360986384173;
//      initialize();
    }

    corrfactype type() const {return NuclearCorrelationFactor::Slater2;}

private:

    /// the length scale parameter
    double a;
    double a0, a1, a2, a3;
    double eprec_;

    double a_param() const {return a;}
    double eprec_param() const {return eprec_;}

    /// first derivative of the correlation factor wrt (r-R_A)

    /// \f[
    ///     Sr_div_S = \frac{1}{S(r)}\frac{\partial S(r)}{\partial r}
    /// \f]
    double Sr_div_S(const double& r, const double& Z) const {
       const double x=r*Z;
       double result=0.0;
       if (x<0.5) {
           result=(-1.1979516162870975156908594435 + x*(1.0070586513226257661408223584 +
                   x*(19.236867748975166091477130269 + x*
                       (-16.639734049500545664179782348 +
                         x*(-38.15233183363957585088408789 +
                            x*(39.66333395908636073317182043 +
                               x*(16.167248815100997375252832353 +
                                  x*(-29.61684122644152074993968583 + 8.898173455627668083163411636*x))))))))/
              (1.1979516162870999616326674855 + x*(0.19089296496164172376540222621 +
                   x*(2.9733505394129816050169313175 + x*
                       (-6.636565088905259599724815599 +
                         x*(9.293774735658646227823596031 +
                            x*(-1.4011704699354161140569248118 +
                               x*(0.9153443813261299256693816943 + x*(-0.40663263966467702847337319933 + 1.*x))))))));
       } else if (x<1.0) {
           result=(-0.6996334787908438908384451809 + x*(1.6739090275472449488185073634 +
                   x*(10.746468402929903692838559804 + x*
                       (-38.13312390182369461972818434 +
                         x*(48.7794572230824205065048821 +
                            x*(-32.35770194320264406507748781 +
                               x*(12.114812980957341479412813182 +
                                  (-2.474602981419607904182347005 + 0.2192823207376552488598805961*x)*x)))))))/
              (0.7047879060549777406562327478 + x*(-1.0885606474007758521466208277 +
                   x*(2.3949527046868955817544996614 + x*
                       (-5.2326125795954354171803378863 +
                         x*(9.988928765335076162902037582 +
                            x*(-10.851226471799446913993735236 +
                               x*(8.228563890293023705695887457 + x*(-3.6740283803242544464629411187 + 1.*x))))))));

       } else if (x<2.0) {
           result=(-22.731091754791490965743643 + x*(198.13483052286515834419613578 +
                   x*(-690.0023018140463118174311835 + x*
                       (1350.3314091552579925748394386 +
                         x*(-1704.0063360514264876129401637 +
                            x*(1483.6118281692660776304793586 +
                               x*(-921.7259978895374365715474525 +
                                  x*(413.8060684165483072719561012 +
                                     x*(-133.67284649619509878961495647 +
                                        x*(30.37208097912233751366804345 +
                                           x*(-4.616164411082228495071778853 +
                                              (0.4220701699125528335924949335 - 0.017582171160843229589975954923*x)*x))))))
                            )))))/
              (11.95594483337282894496679944 + x*(-52.992289110192354102937515453 +
                   x*(115.77912914291025481108315676 + x*
                       (-153.71655196418822134407095039 +
                         x*(136.0779888549758265726955482 +
                            x*(-81.88552528320974307041649157 +
                               x*(33.023913851492798376215534639 + x*(-8.18791648359038111339716106 + 1.*x))))))));
       } else if (x<5) {
           result=(-15.110351601549743764398703643 + x*(39.19284505285046134869389895 +
                   x*(-44.49541217462206268897968033 + x*
                       (28.833235663865859358742627437 +
                         x*(-11.672397147090461991508191129 +
                            x*(3.0517969964993992381049960783 +
                               x*(-0.5149331851643674114165937083 +
                                  x*(0.05407842584625129663110476315 +
                                     (-0.003211389903251499354946305412 + 0.00008234207573839682757066626554*x)*x))))))))/
              (-4804.3735892287918171237312648 + x*(16372.94526140132941081050504 +
                   x*(-24784.660967375785019931130962 +
                      x*(21930.074605973660411015454008 +
                         x*(-12533.587957947537258083416154 +
                            x*(4813.9823109812998150012618667 +
                               x*(-1247.2221256795763422176221957 +
                                  x*(211.28585970267282866335366585 + x*(-21.383294141423722170975899564 + 1.*x)))))))));

       } else if (x<10.0) {
           result=(0.33049674500398019802879378424 + x*(-0.2634588149358958048067954871 +
                   x*(0.08765870931901400539885493498 +
                      x*(-0.015920587444049065607352144674 +
                         x*(0.0017139918426058863178470442334 +
                            x*(-0.0001097761225862026171293155748 +
                               (3.8824367182990408259477813013e-6 - 5.8593762991355372983885431664e-8*x)*x))))))/
              (-4656.1829962732859661807172159 + x*(3784.8832011740610827477157776 +
                   x*(-1316.5042476206107411914341875 +
                      x*(243.26488236896925862692841368 + x*(-23.673063458437062695405038615 + 1.*x)))));

       } else {
           result=0.0;
       }
       return result*Z;
    }

    /// second derivative of the correlation factor wrt (r-R_A)

    /// \f[
    ///     result = \frac{1}{S(r)}\frac{\partial^2 S(r)}{\partial r^2}
    /// \f]
    double Srr_div_S(const double& r, const double& Z) const {
        MADNESS_EXCEPTION("no Srr_div_S in Slater2 yet",0);
        const double x=r/Z;
        double result=0.0;

        if (x<0.5) {
            result=(-0.14665029161753243 + x*(-2.438331591325405 +
                    x*(9.244987909837434 + x*(-6.077373138457289 +
                          x*(-38.2823924333491 + x*(138.81099430976346 + x*(-176.59210562913472 + 76.71323219547432*x)))))))/
               (-0.073325145688678 + x*(0.055284783278928766 +
                    x*(-0.38171185598768664 + x*(0.98230994468608 +
                          x*(-1.8242209615661105 + x*(2.0483842387137976 + x*(-2.065626828932944 + x)))))));
        } else if (x<2.0) {
            result=(11.055991572871836 + x*(-49.56142555819504 +
                    x*(81.11532424072232 + x*(-62.5164069864349 +
                          x*(22.483889766111915 + x*(-1.7486962456989592 +
                                x*(-1.1965980618432754 + x*
                                    (0.3310280993656777 + (-0.021861646054351112 - 0.0005222859475288226*x)*x))))))))/
               (1.9146915680110965 + x*(-7.8191318497581 +
                    x*(11.609421279513054 + x*(1.9199063286845102 +
                          x*(-32.75191958774527 + x*(55.33307318968936 +
                                x*(-48.9961140002982 + x*(25.63844194943869 + x*(-7.550860919558675 + x)))))))));
        } else if (x<3.0) {
            result=(171.38945923492457 + x*(-571.7445731155906 +
                    x*(804.203889750064 + x*(-628.3706259898612 +
                          x*(299.66517050162764 + x*(-89.54399836703575 +
                                x*(16.406590261886663 + (-1.6888011183975509 + 0.0749336047302158*x)*x)))))))/
               (-20.810871372256862 + x*(30.25674456230311 +
                    x*(48.32722302741074 + x*(-177.96006488851197 +
                          x*(218.14223693942046 + x*(-145.51414612997576 +
                                x*(56.05346565184985 + x*(-11.74702483306771 + x))))))));
        } else if (x<10.0) {
            result=(-0.31469853422390565 + x*(0.17653080765857593 +
                    x*(-0.03498004026043187 + (0.0029624700390950953 - 0.00009150922551135921*x)*x)))/
               (-964.1111495364693 + x*(917.2628553541806 + x*(-286.2195080147528 + x*(27.031829359700403 + x))));
        } else {
            result=0.0;
        }
        return result*Z*Z;
    }

    /// third derivative of the correlation factor wrt (r-R_A)

    /// \f[
    ///    result = \frac{1}{S(r)}\frac{\partial^3 S(r)}{\partial r^3}
    /// \f]
    double Srrr_div_S(const double& r, const double& Z) const {
        MADNESS_EXCEPTION("no Srrr_div_S in Slater2 yet",0);
        return 0.0;
    }

    /// the nuclear correlation factor
    double S(const double& r, const double& Z) const {
        const double arZ2=a*a*r*r*Z*Z;

        return 1./(1. - a* (a3 *exp(-32.* arZ2) + a2 *exp(-16. *arZ2) +
            a1 *exp(-4. *arZ2) + a0* exp(-arZ2)) *r *Z);
    }

    /// radial part first derivative of the nuclear correlation factor
    coord_3d Sp(const coord_3d& vr1A, const double& Z) const {
        MADNESS_EXCEPTION("no Sp in Slater2 yet",0);
        return smoothed_unitvec(vr1A);;
    }

    /// second derivative of the nuclear correlation factor
    double Spp_div_S(const double& r, const double& Z) const {
        const double x=r*Z;
        double result=0.0;
        if (x<0.5) {
            result=(-0.010833147179093136088691 + x*(-0.2005540951884271210258 +
                    x*(0.215193550354380169011 + x*(-0.6013634873121505683164 +
                          x*(0.6290000143609286802339 + x*
                              (0.15504953775383744985661 +
                                x*(0.006723696817615123570055 + (2.521910278924732248555 - 3.62311004156467862173*x)*x))))))
                 )/(0.0054165735895465579921527 + x*(0.0034243050602957265391176 +
                    x*(0.047792190243762839156508 + x*(-0.044770230848656988094998 +
                          x*(0.20946109160068066500479 +
                             x*(-0.33752710414970930810664 +
                                x*(0.8614789109642606617794 +
                                   x*(-1.1989906523498382686473 +
                                      x*(2.1016558728558497512196 +
                                         x*(-2.4663973303357583706372 +
                                            x*(2.806269670725901720259 + x*(-1.8654566011104620086627 + 1.*x))))))))))));
        } else if (x<2.0) {
            result=(-0.3408474284999150241901 + x*(-23.20114558485359215082 +
                    x*(127.90118185490392782976 + x*(-303.3905892499787218405 +
                          x*(421.5681782070726183881 + x*
                              (-383.8054868132624549274 +
                                x*(237.0618799426735353961 +
                                   x*(-98.30314370526024000921 +
                                      x*(25.72795549114537061185 +
                                         (-3.653125635621062884781 + 0.16704755059765827776793*x)*x)))))))))/
               (0.35861819748491908821352 + x*(0.29418412211774117470687 +
                    x*(-4.9425982883540078322377 + x*(10.864497472121939642753 +
                          x*(-4.2180357101768293711788 +
                             x*(-19.994189523998234398462 +
                                x*(42.483390552511029623166 +
                                   x*(-41.51356008142083531496 +
                                      x*(23.185896412349280264989 + x*(-7.187719155267694998125 + 1.*x))))))))));

        } else if (x<5)  {
            result=(-304.51326582248542029078 + x*(1013.7720390512033153404 +
                    x*(-1573.7177496270198933538 + x*(1505.495281269806281295 +
                          x*(-993.9662953842932734088 + x*
                              (480.34380741644015361717 +
                                x*(-175.46038947594003912346 +
                                   x*(49.160823365197951025703 +
                                      x*(-10.589807840220615651219 +
                                         x*(1.739785866520397032575 +
                                            x*(-0.21412427551777738676739 +
                                               x*(0.019116778345018873388539 +
                                                  x*(-0.0011693079554266859866836 +
                                                   (0.000043838923967203957503451 - 7.598724762940409361711e-7*x)*x)))))))))
                          ))))/(123.34487670369728721052 +
                 x*(-308.35618558323958114146 + x*(331.76553347658813236246 +
                       x*(-193.83530435236648630731 + x*(65.14983453123805040807 + x*(-12.044800350064563513271 + 1.*x))))));
        } else if (x<8.0) {
            result=(-2.3215362136491996249812 + x*(0.7990997950054717134968 +
                    x*(-0.28194115146517833450165 + x*(0.06388962085634121400231 +
                          x*(-0.009634165597538066812885 +
                             x*(0.0009661572170495110638597 +
                                x*(-0.000062113469110723594190983 +
                                   (2.3225165365168763415415e-6 - 3.8481219231104544641674e-8*x)*x)))))))/
               (0.972787508998877579084 + 1.*x);

        } else {
            result=-1.0/x;
        }

        return result*(Z*Z);
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
        MADNESS_EXCEPTION("no U2X_spherical in Slater2 yet",0);
        return 0.0;
    }

};



/// A nuclear correlation factor class

/// This ncf make the first 3 moments of the inverse ncf vanish
/// int \phi (1-1/S) r^(2+n) dr =0    for n=0,1,2
/// The first  moment of the orbitals is
class Slater3 : public NuclearCorrelationFactor {
public:
    /// ctor

    /// @param[in]  world   the world
    /// @param[in]  mol molecule with the sites of the nuclei
    Slater3(World& world, const Molecule& mol, const double aa)
        : NuclearCorrelationFactor(world,mol), a(1.5) {

        if (aa!=0.0) a=aa;
        a=0.4;
        eprec_=mol.get_eprec();

        if (world.rank()==0) {
            print("\nconstructed nuclear correlation factor of the form");
            print("  S_A = ");
            print("    a = ",a);
            print("with eprec ",eprec_);
            print("which is of Slater3 type\n");
        }
        a0=0.0047309349808059;
        a1=-0.115440391059631;
        a2=2.14968834839076;
        a3=-4.53897889231193;
//      initialize();
    }

    corrfactype type() const {return NuclearCorrelationFactor::Slater3;}

private:

    /// the length scale parameter
    double a;
    double a0, a1, a2, a3;
    double eprec_;

    double a_param() const {return a;}
    double eprec_param() const {return eprec_;}

    /// first derivative of the correlation factor wrt (r-R_A)

    /// \f[
    ///     Sr_div_S = \frac{1}{S(r)}\frac{\partial S(r)}{\partial r}
    /// \f]
    double Sr_div_S(const double& r, const double& Z) const {
       const double x=r*Z;
       double result=0.0;
       if (x<0.5) {
           result=(-0.13191445937053195791469554343 + x*(0.09928973461204942027568722556 +
                   x*(2.393777742290063826791864138 + x*
                       (-1.8341554425704035905078672895 +
                         x*(-5.332331713623470812441216196 +
                            x*(4.5850609334456264295492114 +
                               x*(2.711038107058855340392915941 +
                                  (-3.129474489515047537713850926 + 0.5068847030910685977586590574*x)*x)))))))/
              (0.13191445937391481005586179159 + x*(0.032624722486463481660513562443 +
                   x*(0.32612860057952709073289334962 + x*(-0.80212861242716635437561315 + 1.*x))));
       } else if (x<1.0) {
           result=(-0.3705996720744641254410656311 + x*(5.517734311872378192036611794 +
                   x*(-22.04015943640008607477422254 + x*
                       (16.923367295344266634131250789 +
                         x*(50.16282905443221087893893698 +
                            x*(-113.530556016545852609488024 +
                               x*(94.9950115280579179277585287 +
                                  x*(-37.24416566747294510025117535 + 5.739617206449810692975370557*x))))))))/
              (0.20196846587516713712818777006 + x*(-2.6351780711962187847625743893 +
                   x*(5.2125426349301667587409648554 + x*(-4.672790823157926394295403362 + 1.*x))));
       } else if (x<2.0) {
           result=(-4.060071632890521020776028028 + x*(29.29476590630533436413225694 +
                   x*(-77.33930793878901483580655457 + x*
                       (104.45775847253396879841685332 +
                         x*(-82.20611815859813720156031973 +
                            x*(39.54098530966696880370047964 +
                               x*(-11.513673494828296910385630831 +
                                  (1.8717134809001977710593656281 - 0.1308015385306826163168046533*x)*x)))))))/
              (1.3841172326502928368790629431 + x*(-3.8518925603315649599260058031 +
                   x*(5.1282494596099883424656963236 + x*(-3.1658250781763381131448892233 + 1.*x))));
       } else if (x<5) {
           result=(-171.58411978465703200026501523 + x*(341.5024952676718062173219621 +
                   x*(-284.6945621913430319025515047 + x*
                       (132.03355137686042595828626718 +
                         x*(-37.59893937744486923882200406 +
                            x*(6.762552651436994192184286204 +
                               x*(-0.7520957529877642163873328149 +
                                  (0.04736882696195653016654533714 - 0.0012953919346778370977443633728*x)*x)))))))/
              (-601.4912658422814682506796936 + x*(736.3879983492235803604316809 +
                   x*(-276.9799477558714387956558959 + x*(39.636365833366267879413872907 + 1.*x))));
       } else if (x<10.0) {
           result=(42.764329344142927033566096 + x*(-42.72537108145272496336977868 +
                   x*(17.793613538853565104879031506 + x*
                       (-4.10978263512895636075182383 +
                         x*(0.5816928240486467790962401178 +
                            x*(-0.05198398203928500408518447067 +
                               x*(0.0028759205525729612966430590565 +
                                  (-0.00009029555056185369006841443973 + 1.2341550692875554624506226952e-6*x)*x)))))))/
              (-2732.7151415272664113338495931 + x*(1681.5601749000608060193850311 +
                   x*(-282.62117789397501758814092601 + x*(15.127551136919081784469195983 + 1.*x))));
       } else {
           result=0.0;
       }
       return result*Z;
    }

    /// second derivative of the correlation factor wrt (r-R_A)

    /// \f[
    ///     result = \frac{1}{S(r)}\frac{\partial^2 S(r)}{\partial r^2}
    /// \f]
    double Srr_div_S(const double& r, const double& Z) const {
        MADNESS_EXCEPTION("no Srr_div_S in Slater2 yet",0);
        return 0.0;
    }

    /// third derivative of the correlation factor wrt (r-R_A)

    /// \f[
    ///    result = \frac{1}{S(r)}\frac{\partial^3 S(r)}{\partial r^3}
    /// \f]
    double Srrr_div_S(const double& r, const double& Z) const {
        MADNESS_EXCEPTION("no Srrr_div_S in Slater2 yet",0);
        return 0.0;
    }

    /// the nuclear correlation factor
    double S(const double& r, const double& Z) const {
        const double arZ2=a*a*r*r*Z*Z;

        return 1./(1. - a* (a3 *exp(-32.* arZ2) + a2 *exp(-16. *arZ2) +
            a1 *exp(-4. *arZ2) + a0* exp(-arZ2)) *r *Z);
    }

    /// radial part first derivative of the nuclear correlation factor
    coord_3d Sp(const coord_3d& vr1A, const double& Z) const {
        MADNESS_EXCEPTION("no Sp in Slater2 yet",0);
        return smoothed_unitvec(vr1A);;
    }

    /// second derivative of the nuclear correlation factor
    double Spp_div_S(const double& r, const double& Z) const {
        const double x=r*Z;
        double result=0.0;
        if (x<0.5) {
            result=(-0.12979581847164067539457 + x*(-2.540729515966097981699 +
                    x*(5.911265742177835107542 + x*(-1.875754962201413857763 +
                          x*(-23.91977983402994951409 + x*
                              (72.46406499964457092885 +
                                x*(-38.1856595476431114647 +
                                   x*(-188.39560243054782895538 +
                                      x*(560.7247278088910493291 +
                                         x*(-841.676741910136440931 +
                                            x*(691.7929569477466015086 +
                                               x*(-255.2475889953042243532 + 19.506702205754410045516*x))))))))))))/
               (0.064897909235894819102954 + x*(-0.019247139649382412214985 +
                    x*(0.30211614007061175436546 + x*(-0.6254803832468069187977 +
                          x*(1.2292106471373584880133 + x*(-1.1140250700342440103933 + 1.*x))))));
        } else if (x<2.0) {
            result=(-3.862819825426305966186 + x*(4.072392850686869562679 +
                    x*(59.58680988054514443337 + x*(-249.7426872116496813701 +
                          x*(486.5252289133322656545 + x*
                              (-573.9285828106307853127 +
                                x*(439.0558201030174707455 +
                                   x*(-219.273684277393081375 +
                                      x*(68.21962334952048267107 +
                                         x*(-11.12938613204399395078 +
                                            x*(0.13749081164058245340099 +
                                               (0.2481027476019732232912 - 0.0280882777801874543113*x)*x)))))))))))/
               (0.8777752343485302002156 + x*(-4.3471757452009840140341 +
                    x*(10.044033634234093678765 + x*(-13.170871728324875051971 +
                          x*(10.382428392173839871857 + x*(-4.6439628383204890713794 + 1.*x))))));

        } else if (x<5)  {
            result=(-205.17724577332355375955 + x*(539.7414137381013365403 +
                    x*(-636.1387537773395917484 + x*(430.4407023412384833778 +
                          x*(-180.49036425321909670245 +
                             x*(46.49030738608819450429 +
                                x*(-6.310714685380701890625 +
                                   x*(-0.05380836405799401461819 +
                                      x*(0.19520171953956841432588 +
                                         x*(-0.03881369575855166104667 +
                                            x*(0.003933712218791920040035 +
                                               (-0.00021426398134729217952327 + 4.99661968480774366098e-6*x)*x)))))))))))/
               (110.02123018996372258071 + x*(-270.85913703373498571947 +
                    x*(296.72294216351111650754 + x*(-178.27146158738262808471 +
                          x*(61.813927520837713535784 + x*(-11.692231949742778953799 + 1.*x))))));
        } else if (x<8.0) {
            result=(-41.85822742007710025059 + x*(47.49617615166762770364 +
                    x*(-25.751900705558955115845 + x*(6.990304481212185846758 +
                          x*(-1.1494901780740027173797 +
                             x*(0.1088316039917333294884 +
                                x*(-0.00654885640922458565974 +
                                   (0.00022763091528252416026716 - 3.487509599595671998554e-6*x)*x)))))))/
               (129.66634558802494477998 + x*(-134.65891659982782689942 +
                    x*(58.281658190277663877198 + x*(-10.763041551599489410392 + 1.*x))));

        } else {
            result=-1.0/x;
        }

        return result*(Z*Z);
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
        MADNESS_EXCEPTION("no U2X_spherical in Slater2 yet",0);
        return 0.0;
    }

};


/// A nuclear correlation factor class

/// This ncf make the first 3 moments of the inverse ncf vanish
/// int \phi (1-1/S) r^(2+n) dr =0    for n=0,1,2
/// The first  moment of the orbitals is
class Slater4 : public NuclearCorrelationFactor {
public:
    /// ctor

    /// @param[in]  world   the world
    /// @param[in]  mol molecule with the sites of the nuclei
    Slater4(World& world, const Molecule& mol, const double aa)
        : NuclearCorrelationFactor(world,mol), a(1.5) {

        if (aa!=0.0) a=aa;
        a=0.4;
        eprec_=mol.get_eprec();

        if (world.rank()==0) {
            print("\nconstructed nuclear correlation factor of the form");
            print("  S_A = ");
            print("    a = ",a);
            print("with eprec ",eprec_);
            print("which is of Slater4 type\n");
        }
        a0=0.00026287307489446;
        a1=-0.0120339100572805;
        a2=0.53857889236675;
        a3=-2.19270769120950;
        a4=-0.333640065669944;
//      initialize();
    }

    corrfactype type() const {return NuclearCorrelationFactor::Slater4;}

private:

    /// the length scale parameter
    double a;
    double a0, a1, a2, a3, a4;
    double eprec_;

    double a_param() const {return a;}
    double eprec_param() const {return eprec_;}

    /// first derivative of the correlation factor wrt (r-R_A)

    /// \f[
    ///     Sr_div_S = \frac{1}{S(r)}\frac{\partial S(r)}{\partial r}
    /// \f]
    double Sr_div_S(const double& r, const double& Z) const {
       const double x=r*Z;
       double result=0.0;
       if (x<0.5) {
           result=(-0.07068386723702151110058487421 + x*(-0.6127011931833234239129721429 +
                   x*(2.020575792777072758649014659 + x*
                       (1.9160429910753408662112655038 +
                         x*(-7.069544588182586700580465219 +
                            x*(-2.303164116196953782990272502 +
                               x*(14.553912658711387361757479812 +
                                  x*(-6.737876485254911834619879927 +
                                     x*(-8.790325267798580744226022481 +
                                        (10.342280058209660548980531563 - 3.281171567778618174766354851*x)*x)))))))))/
              (0.0706838672405547511179433877 + x*(-0.041417264937098335326443021682 +
                   x*(0.6337351996272023096447983484 + x*(-1.0249691214916015200147300977 + 1.*x))));
       } else if (x<1.0) {
           result=(7.752102408126742027772948886 + x*(515.3884748779432070817451986 +
                   x*(-3785.114049053370363512790652 + x*
                       (11573.653269357007250174045683 +
                         x*(-20700.81469851908261357586363 +
                            x*(24662.02763089243672902995068 +
                               x*(-20591.67787156452744023588338 +
                                  x*(12038.03941901518026826380425 +
                                     x*(-4701.80450159526432748090831 +
                                        (1099.0257198632417098241110462 - 115.87264537565785256600298163*x)*x)))))))))/
              (-23.793040706832289366271368654 + x*(25.034817678232037129299932894 +
                   x*(-6.866774185977390001501825819 + x*(-25.770946827569526318370957436 + 1.*x))));
       } else if (x<2.0) {
           result=(-1.755667370136172118363611603 + x*(10.319739945306707519372112185 +
                   x*(-23.99472405708318817968431242 + x*
                       (29.73273047111509181230794292 +
                         x*(-21.94750315276875755994742732 +
                            x*(10.026282169753579379898945489 +
                               x*(-2.794363079609283252789355071 +
                                  (0.4371362006860917035272358603 - 0.02951582233408706750845856761*x)*x)))))))/
              (0.8398699517290579194589443717 + x*(-2.7087467919805029124917989773 +
                   x*(4.1104027044226414074289490124 + x*(-2.9448172918321876959000464921 + 1.*x))));
       } else if (x<5) {
           result=(-0.10686134510505682051544860432 + x*(0.12814827540963442970161100417 +
                   x*(-0.027969575915317509583936979614 +
                      x*(-0.026278772339311528081038517578 +
                         x*(0.01933646281571081322379656238 +
                            x*(-0.00574388408844287813741135756 +
                               x*(0.0009075527305671596790828358656 +
                                  (-0.00007521656823895024435859701107 + 2.5825398626394252357676740692e-6*x)*x)))))))/
              (16.845171087195944556267472535 + x*(-33.822654726517974665906059565 +
                   x*(24.687134051454911751612699966 + x*(-7.808576133882478605748578906 + 1.*x))));
       } else if (x<10.0) {
           result=(-3.5877773517354716940710975711 + x*(3.4618815162722506862171381004 +
                   x*(-1.4112748521090538774361688171 +
                      x*(0.32141571386613287138092633067 +
                         x*(-0.045056202203846345159551978286 +
                            x*(0.0039992024854323723325370619418 +
                               x*(-0.00022016853763322952894238106724 +
                                  (6.888142120772356058914269295e-6 - 9.390372181025431034413483053e-8*x)*x)))))))/
              (7058.512759633114080833185154 + x*(-4120.8623744273387137386480844 +
                   x*(819.4978703044841222323673651 + x*(-69.69870232272929667424579037 + 1.*x))));
       } else {
           result=0.0;
       }
       return result*Z;
    }

    /// second derivative of the correlation factor wrt (r-R_A)

    /// \f[
    ///     result = \frac{1}{S(r)}\frac{\partial^2 S(r)}{\partial r^2}
    /// \f]
    double Srr_div_S(const double& r, const double& Z) const {
        MADNESS_EXCEPTION("no Srr_div_S in Slater2 yet",0);
        return 0.0;
    }

    /// third derivative of the correlation factor wrt (r-R_A)

    /// \f[
    ///    result = \frac{1}{S(r)}\frac{\partial^3 S(r)}{\partial r^3}
    /// \f]
    double Srrr_div_S(const double& r, const double& Z) const {
        MADNESS_EXCEPTION("no Srrr_div_S in Slater2 yet",0);
        return 0.0;
    }

    /// the nuclear correlation factor
    double S(const double& r, const double& Z) const {
        const double arZ2=a*a*r*r*Z*Z;

        return 1./(1. + a4 *exp(-64.* arZ2) - a* (a3 *exp(-32.* arZ2) + a2 *exp(-16. *arZ2) +
            a1 *exp(-4. *arZ2) + a0* exp(-arZ2)) *r *Z);
    }

    /// radial part first derivative of the nuclear correlation factor
    coord_3d Sp(const coord_3d& vr1A, const double& Z) const {
        MADNESS_EXCEPTION("no Sp in Slater2 yet",0);
        return smoothed_unitvec(vr1A);;
    }

    /// second derivative of the nuclear correlation factor
    double Spp_div_S(const double& r, const double& Z) const {
        const double x=r*Z;
        double result=0.0;
        if (x<0.5) {
            result=(0.08160513692633673675281 + x*(-0.406803333555475293091 +
                    x*(-0.572329537607280547045 + x*(0.6530142198087404146841 +
                          x*(1.2158639720851370108142 + x*
                              (-2.95168393483182341841 +
                                x*(-5.063199434299873594318 +
                                   x*(8.605294288417414536709 +
                                      x*(112.2237152794543141797 +
                                         x*(-352.4166598414874451396 +
                                            x*(400.4483361289205222421 +
                                               x*(-192.60376274298263505615 + 29.81268625852861860996*x))))))))))))/
               (0.0060984866262568169264358 + x*(0.0031024241923564946491385 +
                    x*(0.1094934045603559664956 + x*(-0.07162597708051576694952 +
                          x*(0.52857164824376918943403 + x*(-0.8381870014256874117882 + 1.*x))))));
        } else if (x<2.0) {
            result=(6.386340744353768134488 + x*(-99.31249333197553419697 +
                    x*(608.3179916282667966875 + x*(-2102.082536915607920133 +
                          x*(4710.277376749569950416 + x*
                              (-7343.783874369994907127 +
                                x*(8288.885445872178344545 +
                                   x*(-6925.124582172778796735 +
                                      x*(4324.902889925594067138 +
                                         x*(-2016.926285562990127038 +
                                            x*(693.4364345143947297294 +
                                               x*(-170.82697499161096755703 +
                                                  x*(28.55620928991024482746 +
                                                   (-2.903663168855036778235 + 0.13566070907519143654634*x)*x)))))))))))))/
               (0.54425553116256354232282 + x*(-3.1735364595672018262563 +
                    x*(8.209095571201105198383 + x*(-11.812011741753344172501 +
                          x*(9.974729876936461875251 + x*(-4.6720798237814916493018 + 1.*x))))));
        } else if (x<5)  {
            result=(-193.32436689940410906234 + x*(644.6932004484338765208 +
                    x*(-993.7648057697846885171 + x*(935.0688833144045003041 +
                          x*(-602.83887432768007524181 +
                             x*(283.3764625172547364657 +
                                x*(-100.67379592800511194195 +
                                   x*(27.497521923985414346403 +
                                      x*(-5.7934483893197921063331 +
                                         x*(0.9336555722268030971358 +
                                            x*(-0.1130081615260589543487 +
                                               x*(0.009944398957014917931144 +
                                                  x*(-0.00060068727204431774196999 +
                                                   (0.000022276870637662989770293 - 3.8249607502195822497451e-7*x)*x))))))))
                             )))))/
               (104.07433887751706598507 + x*(-271.84733680877016654819 +
                    x*(303.31309716381772833855 + x*(-182.52939199786110784004 +
                          x*(62.8579306823696883177 + x*(-11.848264620838455961611 + 1.*x))))));
        } else if (x<8.0) {
            result=-0.790398790424240025553667 + x*(0.248144434031831201300197 +
                    x*(-0.0386835359768122294199345 + (0.00299476480219362694548187 - 0.0000921229213128911629073872*x)*x));
        } else {
            result=-1.0/x;
        }

        return result*(Z*Z);
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
        MADNESS_EXCEPTION("no U2X_spherical in Slater2 yet",0);
        return 0.0;
    }

};



/// A nuclear correlation factor class

/// This ncf make the first 3 moments of the inverse ncf vanish
/// int \phi (1-1/S) r^(2+n) dr =0    for n=0,1,2
/// The first  moment of the orbitals is
class Slater5 : public NuclearCorrelationFactor {
public:
    /// ctor

    /// @param[in]  world   the world
    /// @param[in]  mol molecule with the sites of the nuclei
    Slater5(World& world, const Molecule& mol, const double aa)
        : NuclearCorrelationFactor(world,mol), a(1.5) {

        if (aa!=0.0) a=aa;
        eprec_=mol.get_eprec();

        if (world.rank()==0) {
            print("\nconstructed nuclear correlation factor of the form");
            print("  S_A = ");
            print("    a = ",a);
            print("with eprec ",eprec_);
            print("which is of Slater5 type\n");
        }
        a0=(32. - 9.*constants::pi)/(-32. + 4.*a*sqrt(constants::pi) + 9.*constants::pi);
        a1=(-4.*sqrt(constants::pi))/(-32. + 4.*a*sqrt(constants::pi) + 9.*constants::pi);
        a2=(2*(-8 + 3*constants::pi))/(-32 + 4*a*sqrt(constants::pi) + 9*constants::pi);
        a3=1.0;
        a4=0.0;
//      initialize();
    }

    corrfactype type() const {return NuclearCorrelationFactor::Slater5;}

private:

    /// the length scale parameter
    double a;
    double a0, a1, a2, a3, a4;
    double eprec_;

    double a_param() const {return a;}
    double eprec_param() const {return eprec_;}

    /// first derivative of the correlation factor wrt (r-R_A)

    /// \f[
    ///     Sr_div_S = \frac{1}{S(r)}\frac{\partial S(r)}{\partial r}
    /// \f]
    double Sr_div_S(const double& r, const double& Z) const {
       const double x=r*Z;
       const double sqrtpi=sqrt(constants::pi);
       double result=0.0;
       if (x<0.5) {
           result=(-2.471915184518348613836660032 + x*(0.7126792244376174005800069018 +
                   x*(12.521422865736285071320672141 + x*
                       (-14.514166552473522395140007259 +
                         x*(1.2156233600876834434198345565 +
                            x*(5.436890184287457901996634366 +
                               x*(-5.244936646242548877872743618 +
                                  (3.857137256899563063365266899 - 1.276176888126795974157664777*x)*x)))))))/
              (2.4719151845176747442598629577 + x*(-4.1009807911712621965075680654 +
                   x*(5.9388028673718398824585590505 + x*(-1.7981240187405409463167012747 + 1.*x))));
       } else if (x<1.0) {
           result=-0.408813397534404706124750775015 + x*(-9.82016783251024515165438524123 +
                   x*(51.4058720484530402087095676174 + x*(-76.8170298920974922602185039337 +
                         x*(-349.220930066954103731954266095 +
                            x*(2623.44169569874383382908725858 +
                               x*(-8347.44116043409850099417026266 +
                                  x*(16204.8342316471972124560719328 +
                                     x*(-20941.9543382160545468632452438 +
                                        x*(18528.2636553890444498861878894 +
                                           x*(-11166.6980931656581014036702897 +
                                              x*(4406.53528141840581868467709781 +
                                                 x*(-1030.74393475495301599426836402 + 108.691787299456964602967032468*x))))))
                                  ))))));
       } else if (x<2.0) {
           result=21.4946711002191328836113309162 + x*(-208.358018764159297898424350088 +
                   x*(854.807421991910142973947647236 + x*(-2059.07478691176730558469672808 +
                         x*(3351.29153493949233714482880454 +
                            x*(-3961.59321365753068511534648115 +
                               x*(3522.48398054210471536373948905 +
                                  x*(-2382.20957950394602975930835527 +
                                     x*(1218.97253749082431962749090695 +
                                        x*(-463.36753670752266857306494521 +
                                           x*(126.567970640867902461603976223 +
                                              x*(-23.4464099021492838270721698632 +
                                                 (2.63504968352350040741514547447 - 0.135565200732828177192193110211*x)*x)))))
                                  ))))));
       } else if (x<5) {
           result=(-0.61293299824721639101617863337 + x*(1.7159447901393692898512054041 +
                   x*(-2.2009899438241157720415505532 + x*
                       (1.7102653095552973026058666205 +
                         x*(-0.8966040157666084046140551118 +
                            x*(0.33407072228375245774791419436 +
                               x*(-0.09070742076107965185901163327 +
                                  x*(0.018083060533399376129706976765 +
                                     x*(-0.0026267951569006281542156400416 +
                                        x*(0.00027114269213617506595420204211 +
                                           x*(-0.000018877151555489312161248681046 +
                                              (7.958637943954770417412713077e-7 - 1.5365833542091685343331202907e-8*x)*x))))
                                  )))))))/
              (95.41248708641002657127258502 + x*(-256.33973108648842134907149745 +
                   x*(292.27008700578003373324917243 + x*
                       (-180.58279785318853443503147244 +
                         x*(63.79195168219952007729827881 + x*(-12.23469202161081875321365099 + 1.*x))))));
       } else if (x<10.0) {
           result=0.0;
       } else {
           result=0.0;
       }
       return result*Z;
    }

    /// second derivative of the correlation factor wrt (r-R_A)

    /// \f[
    ///     result = \frac{1}{S(r)}\frac{\partial^2 S(r)}{\partial r^2}
    /// \f]
    double Srr_div_S(const double& r, const double& Z) const {
        MADNESS_EXCEPTION("no Srr_div_S in Slater2 yet",0);
        return 0.0;
    }

    /// third derivative of the correlation factor wrt (r-R_A)

    /// \f[
    ///    result = \frac{1}{S(r)}\frac{\partial^3 S(r)}{\partial r^3}
    /// \f]
    double Srrr_div_S(const double& r, const double& Z) const {
        MADNESS_EXCEPTION("no Srrr_div_S in Slater2 yet",0);
        return 0.0;
    }

    /// the nuclear correlation factor
    double S(const double& r, const double& Z) const {
        const double arZ2=a*a*r*r*Z*Z;

        return 1.0 + exp(- a3 *arZ2) * (a0 + a1*sqrt(arZ2) + a2 * arZ2);
    }

    /// radial part first derivative of the nuclear correlation factor
    coord_3d Sp(const coord_3d& vr1A, const double& Z) const {
        MADNESS_EXCEPTION("no Sp in Slater2 yet",0);
        return smoothed_unitvec(vr1A);;
    }

    /// second derivative of the nuclear correlation factor
    double Spp_div_S(const double& r, const double& Z) const {
        const double x=r*Z;
        double result=0.0;
        if (x<0.5) {
            result=(11.809850954220120025809 + x*(-120.50725051649522040935 +
                    x*(211.6581939748748187491 + x*(2.146276504985675190597 +
                          x*(-264.8489104887500673513 + x*
                              (166.42770285759379550123 +
                                x*(43.55699174422308217551 +
                                   x*(-73.65025086499743002264 +
                                      x*(42.41976785723404166436 + x*(-27.94579594715360849353 + 9.264627384753940992589*x)))
                                   )))))))/
               (7.589494262513949053356 + x*(-20.09262959056524715904 +
                    x*(28.197568894588136619277 + x*(-17.814651231428548630769 + x*(0.7864429300819081919902 + 1.*x)))));
        } else if (x<1.0) {
            result=(2.118023027542999184338 + x*(-20.89000925979305982979 +
                    x*(31.92869238376278088209 + x*(1.2196130271017273135482 +
                          x*(-39.62704743125086663774 + x*
                              (34.79908667325655900563 + x*
                                 (-11.598411270393951667691 + (0.9148919877212963001764 + 0.17782129587583002911975*x)*x)))))
                    ))/(1.3630612305245051344777 + x*(-3.1562115067046947606073 +
                    x*(4.8011339776073516275952 + x*(-3.0055450707600381389398 + 1.*x))));
        } else if (x<2.0)  {
            result=(-1.3709533078772464496852 + x*(19.698229852843597246291 +
                    x*(-158.44181734366018476062 + x*(466.6618679217673119632 +
                          x*(-716.7627291868396935349 + x*
                              (659.8135580917825557185 + x*
                                 (-386.4430511615393489175 +
                                   x*(145.9674185966482717967 +
                                      x*(-34.59746981559744482102 + (4.698555126069378501543 - 0.2798138855826339542861*x)*x)
                                      ))))))))/
               (1.3469488281289268313003 + x*(-2.8301220597504922749822 +
                    x*(4.626067576658674715774 + x*(-3.0369327366945672542512 + 1.*x))));
        } else if (x<5.0)  {
            result=(-586.77335121801773174473 + x*(2045.6767845999824679689 +
                    x*(-3152.7763681477718654566 + x*(2802.6952589837051473105 +
                          x*(-1572.066692259934262199 + x*
                              (570.27838541454454567527 +
                                x*(-130.94361802449937581313 +
                                   x*(17.497550669721172100695 +
                                      x*(-1.0689268436583516717463 +
                                         x*(0.00665591958083010175981 +
                                            x*(-0.00043671807356393128709827 +
                                               (0.000017456513416762374363243 - 3.2112916648191429882317e-7*x)*x)))))))))))/
               (26.509296446267914711487 + x*(497.66877063050784386019 +
                    x*(-1909.4561164127667423565 + x*(3026.9776762477333943983 +
                          x*(-2724.2259586417457526236 + x*
                              (1537.0981146424259014841 +
                                x*(-558.81078034516961755499 +
                                   x*(128.13891846973022697575 + x*(-16.985587891946482709663 + 1.*x)))))))));
        } else if (x<8.0) {
            result=-1.0/x;
        } else {
            result=-1.0/x;
        }

        return result*(Z*Z);
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
        MADNESS_EXCEPTION("no U2X_spherical in Slater2 yet",0);
        return 0.0;
    }

};

class poly4erfc : public NuclearCorrelationFactor {
public:
    /// ctor

    /// @param[in]  world   the world
    /// @param[in]  mol molecule with the sites of the nuclei
    poly4erfc(World& world, const Molecule& mol, const double aa)
        : NuclearCorrelationFactor(world,mol), a(1.5) {

        if (aa!=0.0) a=aa;
        eprec_=mol.get_eprec();

        if (world.rank()==0) {
            print("\nconstructed nuclear correlation factor of the form");
            print("  S_A = 1 + (a0 + a1 arZ + a2 (arZ)^2 + a3 (arZ)^3 + a4 (arZ)^4) erfc(arZ)");
            print("    a = ",a);
            print("with eprec ",eprec_);
            print("which is of poly4erfc type\n");
        }
        const double pi32=std::pow(constants::pi,1.5);
        const double sqrtpi=sqrt(constants::pi);
        const double Pi=constants::pi;

        if (a==0.5) {
            a0=0.5083397721116242769;
            a1=-2.4430795355664112811;
            a2=3.569312300653802680;
            a3=-1.9812471972342746507;
            a4=0.3641705622093696564;
        } else if (a==1.0) {
            a0=0.20265985404508529127;
            a1=-0.9739826967339938056;
            a2=1.4229779953809877198;
            a3=-0.7898639647077711196;
            a4=0.14518390461225107425;

        } else {
            print("invalid parameter a for poly4erfc: only 0.5 and 1.0 implemented");
            MADNESS_EXCEPTION("stupid you",1);
        }
        //      initialize();
    }

    corrfactype type() const {return NuclearCorrelationFactor::poly4erfc;}

private:

    /// the length scale parameter
    double a;
    double a0, a1, a2, a3, a4;
    double eprec_;

    double a_param() const {return a;}
    double eprec_param() const {return eprec_;}

    /// first derivative of the correlation factor wrt (r-R_A)

    /// \f[
    ///     Sr_div_S = \frac{1}{S(r)}\frac{\partial S(r)}{\partial r}
    /// \f]
    double Sr_div_S(const double& r, const double& Z) const {
       const double x=r*Z;

       double result=0.0;
       if (a==0.5) {
           if (x<1.0) {
                result=(-17.97663543396820624361586474 + x*(32.78290319470982346841067868 +
                        x*(-18.158574783628271659713638233 +
                           x*(2.472138984374094343735335913 +
                              x*(0.5516975358315341628276502285 +
                                 x*(-0.008573693952875097234391220137 +
                                    x*(-0.05596791202351071993748992739 +
                                       (0.002673799219133696315436690424 +
                                          0.0013386538660557369902632289083*x)*x)))))))/
                   (17.976635433967702922140233083 + x*
                      (-13.062204300852085089323266568 +
                        x*(16.397871971437618641239835027 + x*(-5.383337491559214163188757918 + 1.*x))));
            } else if (x<2.0) {
                result=(-16.53370050883888159389958126 + x*(30.04151304875517461538549269 +
                        x*(-16.692529697855029750948871179 +
                           x*(2.50341008323651011875249567 +
                              x*(0.3106921665634860719234742532 +
                                 x*(0.08721948207311506458903445571 +
                                    x*(-0.10041387133168708232852057948 +
                                       (0.02000987266876476192949541524 - 0.0012508983745483161308604975792*x)*
                                        x)))))))/
                   (16.532205048243702951212522516 + x*
                      (-11.89273747187279634240347945 +
                        x*(15.157537549656745468895369276 + x*(-4.9102960292797655978798640519 + 1.*x))));
            } else if (x<5.0) {
                result=(-2352.191894900273554810118278 + x*(5782.846962269399174183661793 +
                        x*(-5653.246084369776756298278851 +
                           x*(2948.18046377483570925427449 +
                              x*(-913.4583247839311453090142452 +
                                 x*(174.39391722588915386106331206 +
                                    x*(-20.22035127074332315930567933 +
                                       (1.3107321165711966663791114988 - 0.03655666729452579098523876463*x)*x))
                                 )))))/
                   (886.4859678528423041797649741 + x*(269.17746130370931387996124706 +
                        x*(-130.21383548057958115685397713 + x*(37.644499985765056273193347388 + 1.*x))));
            } else if (x<10) {
                result=(2.2759176275121988686860433041 + x*(-1.8014283464827425541637211503 +
                        x*(0.60570955276433317373251991152 +
                           x*(-0.11235368819003308411943926239 +
                              x*(0.01243211635244600976892077538 +
                                 x*(-0.0008211800260491381826149891865 +
                                    x*(0.000029973534470203049782417744015 +
                                       (-4.6423722605763162431293872646e-7 -
                                          1.3224615425412157194986675329e-10*x)*x)))))))/
                   (1039.800013929971888016478838 + x*(-702.5378531848183775210948787 +
                        x*(183.17476380259879599459789974 + x*(-21.74106003575315304073197254 + 1.*x))));
            } else {
                result=0.0;
            }
       } else if (a==1.0) {
           if (x<1.0) {
               result=(-1.6046958001953006847027538457 + x*(5.948945186367159977879486279 +
                       x*(-6.884321742840291285733040882 +
                          x*(2.296896506418919905405783368 +
                             x*(0.616939354810622973212914089 +
                                x*(-0.13679830198890803519235207564 +
                                   x*(-0.3356576872066501398893403439 +
                                      (0.14876798925798674488426727928 - 0.016049886728185297028755535226*x)*x
                                      )))))))/
                  (1.6046957999975126909719683196 + x*
                     (-0.8234878506688316215458988304 +
                       x*(3.4607641859010551501639903314 + x*(-1.8955210085531557309609670978 + 1.*x))));
           } else if (x<2.0) {
               result=(-7.143856421301985019580778813 + x*(33.35568129248075086686087865 +
                       x*(-60.0412343766569246898209957 +
                          x*(55.46407913315151830247939138 +
                             x*(-28.7874770749840240158264326 +
                                x*(8.379317837934083469035061852 +
                                   x*(-1.2103278392957399092107317741 +
                                      (0.04275186003977071121074860478 + 0.005247730112126140726731063205*x)*x
                                      )))))))/
                  (5.4509691924562038998044993827 + x*
                     (-2.2732931206867811068721699796 +
                       x*(4.1530634219989859344450742259 + x*(-2.0183662125874366044951391259 + 1.*x))));
           } else if (x<5.0) {
               result=(-0.4869290414611847276899694883 + x*(1.0513417728218375522016338562 +
                       x*(-0.9694851629317255156038942437 +
                          x*(0.5007460889402078102011673774 +
                             x*(-0.15897436623286639052954575417 +
                                x*(0.03185042477631079356985872518 +
                                   x*(-0.003940899912654543218183049969 +
                                      (0.0002757919409481032696079686825 -
                                         8.369138363906178282041501561e-6*x)*x)))))))/
                  (30.81121613134246115634780413 + x*(-49.989974505724725146933397436 +
                       x*(31.455643953462635691568128729 + x*(-8.992097794824270871044305786 + 1.*x))));
           } else {
               result=0.0;
           }
       }
       return result*Z;
    }

    /// second derivative of the correlation factor wrt (r-R_A)

    /// \f[
    ///     result = \frac{1}{S(r)}\frac{\partial^2 S(r)}{\partial r^2}
    /// \f]
    double Srr_div_S(const double& r, const double& Z) const {
        MADNESS_EXCEPTION("no Srr_div_S in Slater2 yet",0);
        return 0.0;
    }

    /// third derivative of the correlation factor wrt (r-R_A)

    /// \f[
    ///    result = \frac{1}{S(r)}\frac{\partial^3 S(r)}{\partial r^3}
    /// \f]
    double Srrr_div_S(const double& r, const double& Z) const {
        MADNESS_EXCEPTION("no Srrr_div_S in Slater2 yet",0);
        return 0.0;
    }

    /// the nuclear correlation factor
    double S(const double& r, const double& Z) const {
        const double arZ=a*r*Z;
        const double arZ2=a*a*r*r*Z*Z;

        return 1.0 + (a0 +a1* arZ+ a2 * arZ2 + a3*arZ*arZ2 + + a4*arZ2*arZ2) *erfc( a*r*Z);
    }

    /// radial part first derivative of the nuclear correlation factor
    coord_3d Sp(const coord_3d& vr1A, const double& Z) const {
        MADNESS_EXCEPTION("no Sp in Slater2 yet",0);
        return smoothed_unitvec(vr1A);;
    }

    /// second derivative of the nuclear correlation factor
    double Spp_div_S(const double& r, const double& Z) const {
        const double x=r*Z;
        double result=0.0;
        if (a==0.5) {
            if (x<1.0) {
                result=(-37.75186842343823465059 + x*(21.3476988467348903615 +
                        x*(0.014608707333424946750026 +
                           x*(-3.704945314726273312722 +
                              x*(0.08808104845382944292252 +
                                 x*(0.3362428305409206967066 +
                                    x*(-0.02288039625102549092766 +
                                       (-0.017240056622850307571001 + 0.002412740490618117536527*x)*x)))))))/
                   (17.595611361287293183447 + x*(-12.421065066789808273295 +
                        x*(15.957564552156320207938 + x*(-5.0239100389132464317907 + 1.*x))));
            } else if (x<2.0) {
                result=(-35.46082798344982189328 + x*(20.3565414350879797493 +
                        x*(-0.3046488975234456455871 +
                           x*(-3.298067641731523613442 +
                              x*(-0.10712268902574947466554 +
                                 x*(0.4829111116912435121877 +
                                    x*(-0.10089023589818206760275 +
                                       (0.0013899828749998182176948 + 0.0008305150868988335610678*x)*x)))))))/
                   (16.526499668833810286111 + x*(-11.79820182874024593006 +
                        x*(15.115407083582433322342 + x*(-4.8449266197426825174319 + 1.*x))));
            } else if (x<5.0) {
                result=(-414.5311559516023264104 + x*(615.8720414166440747539 +
                        x*(-422.5938440932793094888 + x*
                            (159.72497494584873155352 +
                              x*(-35.6790348104188081907 +
                                 x*(4.658777521872728594702 +
                                    x*(-0.328305094433759490678 +
                                       (0.009162754689309596905172 + 0.00005047926659010755662873*x)*x)))))))/
                   (112.7044543272820830484 + x*(-56.072518762714727894479 +
                        x*(28.188843903409059322224 + x*(-6.519554545057610040741 + 1.*x))));
            } else if (x<10.0) {
                result=(-146.68256559112012287314 + x*(188.52807059353309478385 +
                        x*(-82.07176992590032524431 + x*
                            (18.107697718347802322776 +
                              x*(-2.3221393933622638466979 +
                                 x*(0.18681223946803275939642 +
                                    x*(-0.009601011427143648072501 +
                                       (0.00028623295732553770894583 - 3.77364531155112782235e-6*x)*x)))))))/
                   (633.1319105785227043552 + x*(-574.29230199331494406798 +
                        x*(171.872612865808639376 + x*(-21.906495121260483674385 + 1.*x))));
            } else {
                result=-1.0/x;
            }
        } else if (a==1.0) {
            if (x<1.0) {
                result=(-8.85288955131414420808 + x*(11.359434597010419928515 +
                        x*(-2.982094072176256165405 + x*
                            (-4.924201512880445076103 +
                              x*(0.6773928043289287907009 +
                                 x*(2.061680233090141287213 +
                                    x*(-0.5528612000728412913713 +
                                       (-0.3024276764350212595121 + 0.10765958646570631264003*x)*x)))))))/
                   (1.6731803922436094206275 + x*(-0.8242052614481793487834 +
                        x*(3.5912910648514747558597 + x*(-1.9146644266104277222142 + 1.*x))));
            } else if (x<2.0) {
                result=(-8.538839434207355205544 + x*(-37.85611260277589742674 +
                        x*(155.08466711228234382211 + x*
                            (-241.1247881821171613212 +
                              x*(199.33293375918859716585 +
                                 x*(-96.80750216221383928113 +
                                    x*(27.71704080319043943288 +
                                       (-4.354251431065185902435 + 0.2906781051124817624589*x)*x)))))))/
                   (2.2233877901091612794108 + x*(3.2084406367698851844258 +
                        x*(1.021615116328133792993 + x*(-0.9819667717342250528772 + 1.*x))));
            } else if (x<5.0) {
                result=(-7.338671998425412491091 + x*(18.593867285988629508481 +
                        x*(-19.15344657249203844577 + x*
                            (10.072359942912659732126 +
                              x*(-2.99771628023376895151 +
                                 x*(0.5324416109232007541639 +
                                    x*(-0.05993976172902101376555 +
                                       (0.003887110757665478374745 - 0.00011079212872583089652945*x)*x)))))))/
                   (13.273604111590967002999 + x*(-28.012330064250556933961 +
                        x*(22.161347591665727407914 + x*(-7.617582259533338164664 + 1.*x))));
            } else if (x<10) {
                result=(132.78397650676876308801 + x*(-78.01898251977413923271 +
                        x*(15.292039515801252996504 + x*
                            (-1.000000154745351328125 +
                              x*(1.950300262619412492847e-8 +
                                 x*(-1.6337541591423364318692e-9 +
                                    x*(8.771557330523968791519e-11 +
                                       (-2.7388800346031121159372e-12 + 3.7894572289998844980568e-14*x)*x))))))
                     )/(-4.724326520908917985827e-6 + x*
                      (-132.78397108375581258096 + x*(78.01897976124964688489 +
                           x*(-15.292038699709498552356 + 1.*x))));
            } else {
                result=-1.0/x;
            }
        }
        return result*Z*Z;
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
        MADNESS_EXCEPTION("no U2X_spherical in Slater2 yet",0);
        return 0.0;
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
//		initialize();
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
//		initialize();

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
