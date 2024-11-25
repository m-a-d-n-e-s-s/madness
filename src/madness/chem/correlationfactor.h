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
#include<madness/chem/molecule.h>
#include<madness/chem/potentialmanager.h>
#include<madness/chem/atomutil.h>

namespace madness {

/// ABC for the nuclear correlation factors
class NuclearCorrelationFactor {
public:
	enum corrfactype {None, GradientalGaussSlater, GaussSlater, LinearSlater,
	    Polynomial, Slater, poly4erfc, Two, Adhoc};
	typedef std::shared_ptr< FunctionFunctorInterface<double,3> > functorT;

	/// ctor

	/// @param[in]	world	the world
	/// @param[in]	mol molecule with the sites of the nuclei
	NuclearCorrelationFactor(World& world, const Molecule& mol)
		: world(world), vtol(FunctionDefaults<3>::get_thresh()*0.1)
		, eprec(mol.get_eprec()), molecule(mol) {}

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

    /// smoothing of the potential/step function
    double eprec;

	/// the molecule
	const Molecule& molecule;

protected:
	/// the three components of the U1 potential
	std::vector<real_function_3d> U1_function;

	/// the purely local U2 potential, having absorbed the nuclear pot V_nuc
	real_function_3d U2_function;

private:
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
        if (smoothing==0.0) smoothing=eprec;
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
        if (smoothing==0.0) smoothing=eprec;

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
			for (size_t i=0; i<ncf->molecule.natom(); ++i) {
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
			for (size_t i=0; i<ncf->molecule.natom(); ++i) {
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
        const size_t iatom;
        const int axis;

    public:
        U1_atomic_functor(const NuclearCorrelationFactor* ncf, const size_t atom,
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

    /// \f[
    ///  U1\dot U1 = \frac{\left(S^r_A S^r_B\right)}{S_A S_B} n_A \cdot n_B
    /// \f]
    /// with positive sign!
    class U1_dot_U1_functor : public FunctionFunctorInterface<double,3> {

        const NuclearCorrelationFactor* ncf;

    public:
        U1_dot_U1_functor(const NuclearCorrelationFactor* ncf) : ncf(ncf) {}

        double operator()(const coord_3d& xyz) const {
			std::vector<double> Sr_div_S(ncf->molecule.natom());
			std::vector<coord_3d> unitvec(ncf->molecule.natom());
			for (size_t i=0; i<ncf->molecule.natom(); ++i) {
				const Atom& atom=ncf->molecule.get_atom(i);
				const coord_3d vr1A=xyz-atom.get_coords();
				const double r=vr1A.normf();
				Sr_div_S[i]=ncf->Sr_div_S(r,atom.q);
				unitvec[i]=ncf->smoothed_unitvec(vr1A);
			}

			double result=0.0;
			for (size_t i=0; i<ncf->molecule.natom(); ++i) {
				for (size_t j=0; j<ncf->molecule.natom(); ++j) {
					double tmp=Sr_div_S[i]*Sr_div_S[j];
					if (i!=j) tmp*=inner(unitvec[i],unitvec[j]);
					result+=tmp;
				}
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
			for (size_t i=0; i<ncf->molecule.natom(); ++i) {
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
			for (size_t i=0; i<ncf->molecule.natom(); ++i) {
				const Atom& atom=ncf->molecule.get_atom(i);
				const coord_3d vr1A=xyz-atom.get_coords();
				const double r=vr1A.normf();
//				all_terms[i]=ncf->Sp(vr1A,atom.q)*(1.0/ncf->S(r,atom.q));
				all_terms[i]=ncf->Sr_div_S(r,atom.q)*ncf->smoothed_unitvec(vr1A);
			}

			double result=0.0;
			for (size_t i=0; i<ncf->molecule.natom(); ++i) {
				for (size_t j=0; j<i; ++j) {
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
        const size_t iatom;

    public:
        U2_atomic_functor(const NuclearCorrelationFactor* ncf, const size_t atom)
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
        const size_t iatom;

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
            for (size_t i=0; i<ncf->molecule.natom(); ++i) {
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
        const size_t iatom;
    public:
        square_times_V_functor(const NuclearCorrelationFactor* ncf,
                const Molecule& mol, const size_t iatom1)
            : ncf(ncf), molecule(mol), iatom(iatom1) {}
        double operator()(const coord_3d& xyz) const {
            double result=1.0;
            for (size_t i=0; i<ncf->molecule.natom(); ++i) {
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
        const size_t iatom;
        const int axis;
    public:
        square_times_V_derivative_functor(const NuclearCorrelationFactor* ncf,
                const Molecule& molecule1, const size_t atom1, const int axis1)
            : ncf(ncf), molecule(molecule1), iatom(atom1), axis(axis1) {}
        double operator()(const coord_3d& xyz) const {
            double result=1.0;
            for (size_t i=0; i<ncf->molecule.natom(); ++i) {
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
                for (size_t i=0; i<ncf->molecule.natom(); ++i) {
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
        const size_t iatom;
        const int axis;
    public:
        U3X_functor(const NuclearCorrelationFactor* ncf, const size_t iatom,
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
            for (size_t jatom=0; jatom<ncf->molecule.natom(); ++jatom) {
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
    	//const double eprec=eprec_param();
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


class poly4erfc : public NuclearCorrelationFactor {
public:
    /// ctor

    /// @param[in]  world   the world
    /// @param[in]  mol molecule with the sites of the nuclei
    poly4erfc(World& world, const Molecule& mol, const double aa)
        : NuclearCorrelationFactor(world,mol), a(1.0) {

        if (aa!=0.0) a=aa;
        eprec_=mol.get_eprec();

        if (world.rank()==0) {
            print("\nconstructed nuclear correlation factor of the form");
            print("  S_A = 1 + (a0 + a1 arZ + a2 (arZ)^2 + a3 (arZ)^3 + a4 (arZ)^4) erfc(arZ)");
            print("    a = ",a);
            print("with eprec ",eprec_);
            print("which is of poly4erfc type\n");
        }
        //const double pi32=std::pow(constants::pi,1.5);
        //const double sqrtpi=sqrt(constants::pi);
        //const double Pi=constants::pi;

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


/// this ncf has no information about itself, only U2 and U1 assigned
class AdhocNuclearCorrelationFactor : public NuclearCorrelationFactor {

public:
	/// ctor

	/// @param[in]	world	the world
	/// @param[in]	mol molecule with the sites of the nuclei
	AdhocNuclearCorrelationFactor(World& world, const real_function_3d U2,
		const std::vector<real_function_3d>& U1)
		: NuclearCorrelationFactor(world,Molecule()) {

		U2_function=U2;
		U1_function=U1;

		if (world.rank()==0) {
			print("constructed ad hoc nuclear correlation factor");
		}
	}

	corrfactype type() const {return Adhoc;}

private:

    double Sr_div_S(const double& r, const double& Z) const {
    	MADNESS_EXCEPTION("no Sr_div_S() in AdhocNuclearCorrelationFactor",0);
	    return 0.0;
    }

    double Srr_div_S(const double& r, const double& Z) const {
    	MADNESS_EXCEPTION("no Srr_div_S() in AdhocNuclearCorrelationFactor",0);
	    return 0.0;
    }

    double Srrr_div_S(const double& r, const double& Z) const {
    	MADNESS_EXCEPTION("no Srrr_div_S() in AdhocNuclearCorrelationFactor",0);
	    return 0.0;
    }

    /// the nuclear correlation factor
    double S(const double& r, const double& Z) const {
    	MADNESS_EXCEPTION("no S() in AdhocNuclearCorrelationFactor",0);
    	return 0.0;
    }

    /// radial part first derivative of the nuclear correlation factor
    coord_3d Sp(const coord_3d& vr1A, const double& Z) const {
    	MADNESS_EXCEPTION("no Sp() in AdhocNuclearCorrelationFactor",0);
    	return coord_3d(0.0);
    }

    /// second derivative of the nuclear correlation factor
    double Spp_div_S(const double& r, const double& Z) const {
    	MADNESS_EXCEPTION("no Spp_div_S() in AdhocNuclearCorrelationFactor",0);
    	return 0.0;
    }
};


std::shared_ptr<NuclearCorrelationFactor>
create_nuclear_correlation_factor(World& world,
		const Molecule& molecule,
		const std::shared_ptr<PotentialManager> pm,
		const std::string inputline);

std::shared_ptr<NuclearCorrelationFactor>
create_nuclear_correlation_factor(World& world,
		const Molecule& molecule,
		const std::shared_ptr<PotentialManager> pm,
		const std::pair<std::string,double>& ncf);

}


#endif /* NUCLEARCORRELATIONFACTOR_H_ */
