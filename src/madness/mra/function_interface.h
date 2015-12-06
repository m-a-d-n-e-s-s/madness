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

  $Id: function_factory_and_interface.h 3422 2014-03-24 09:16:15Z 3ru6ruWu $
*/


#ifndef MADNESS_MRA_FUNCTION_INTERFACE_H__INCLUDED
#define MADNESS_MRA_FUNCTION_INTERFACE_H__INCLUDED

#include <madness/tensor/tensor.h>
#include <madness/tensor/gentensor.h>
#include <madness/mra/key.h>

// needed for the TwoElectronInterface
#include <madness/mra/gfit.h>
#include <madness/mra/convolution1d.h>
#include <madness/mra/function_common_data.h>
namespace madness {

	// forward declaration needed for CompositeFunctorInterface
	template<typename T, std::size_t NDIM>
	class FunctionImpl;

    template<typename T, std::size_t NDIM>
    Tensor<T> fcube(const Key<NDIM>&, T (*f)(const Vector<double,NDIM>&), const Tensor<double>&);



	/// Abstract base class interface required for functors used as input to Functions
	template<typename T, std::size_t NDIM>
	class FunctionFunctorInterface {
	public:

	    typedef GenTensor<T> coeffT;
	    typedef Key<NDIM> keyT;
	    typedef T value_type;

	    Level special_level_;

	    FunctionFunctorInterface() : special_level_(6) {}

	    /// adapt the special level to resolve the smallest length scale
	    void set_length_scale(double lo) {
	        double Lmax=FunctionDefaults<NDIM>::get_cell_width().max();
	        double lo_sim=lo/Lmax;  // lo in simulation coordinates;
	        special_level_=Level(-log2(lo_sim));
	    }

	    /// Can we screen this function based on the bounding box information?
	    virtual bool screened(const Vector<double,NDIM>& c1, const Vector<double,NDIM>& c2) const {
	        return false;
	    }

	    /// Does the interface support a vectorized operator()?
	    virtual bool supports_vectorized() const {return false;}

	    virtual void operator()(const Vector<double*,1>& xvals, T* fvals, int npts) const {
	        MADNESS_EXCEPTION("FunctionFunctorInterface: This function should not be called!", 0);
	    }

	    virtual void operator()(const Vector<double*,2>& xvals, T* fvals, int npts) const {
	        MADNESS_EXCEPTION("FunctionFunctorInterface: This function should not be called!", 0);
	    }

	    virtual void operator()(const Vector<double*,3>& xvals, T* fvals, int npts) const {
	        MADNESS_EXCEPTION("FunctionFunctorInterface: This function should not be called!", 0);
	    }

	    virtual void operator()(const Vector<double*,4>& xvals, T* fvals, int npts) const {
	        MADNESS_EXCEPTION("FunctionFunctorInterface: This function should not be called!", 0);
	    }

	    virtual void operator()(const Vector<double*,5>& xvals, T* fvals, int npts) const {
	        MADNESS_EXCEPTION("FunctionFunctorInterface: This function should not be called!", 0);
	    }

	    virtual void operator()(const Vector<double*,6>& xvals, T* fvals, int npts) const {
	        MADNESS_EXCEPTION("FunctionFunctorInterface: This function should not be called!", 0);
	    }

	    /// You should implement this to return \c f(x)
	    virtual T operator()(const Vector<double, NDIM>& x) const = 0;

	    /// Override this to return list of special points to be refined more deeply
	    virtual std::vector< Vector<double,NDIM> > special_points() const {
	        return std::vector< Vector<double,NDIM> >();
	    }

	    /// Override this change level refinement for special points (default is 6)
	    virtual Level special_level() {return special_level_;}

	    virtual ~FunctionFunctorInterface() {}

	    virtual coeffT coeff(const keyT&) const {
	        MADNESS_EXCEPTION("implement coeff for FunctionFunctorInterface",0);
	        return coeffT();
	    }

	    virtual coeffT values(const keyT& key, const Tensor<double>& tensor) const {
	        MADNESS_EXCEPTION("implement values for FunctionFunctorInterface",0);
	        return coeffT();
	    }

	    /// does this functor directly provide sum coefficients? or only function values?
	    virtual bool provides_coeff() const {
	        return false;
	    }

	};



	///forward declaration
	template <typename T, std::size_t NDIM>
	//    void FunctionImpl<T,NDIM>::fcube(const keyT& key, const FunctionFunctorInterface<T,NDIM>& f, const Tensor<double>& qx, tensorT& fval) const {
	void fcube(const Key<NDIM>& key, const FunctionFunctorInterface<T,NDIM>& f, const Tensor<double>& qx, Tensor<T>& fval);

	/// CompositeFunctorInterface implements a wrapper of holding several functions and functors

	/// Use this to "connect" several functions and/or functors and to return their coefficients
	/// e.g. connect f1 and f2 with an addition, you can request the coefficients of any node
	/// and they will be computed on the fly and returned. Mainly useful to connect a functor
	/// with a function, if the functor is too large to be represented in MRA (e.g. 1/r12)
	///
	/// as of now, the operation connecting the functions/functors is simply addition.
	/// need to implement expression templates, if I only knew what that was...
	template<typename T, std::size_t NDIM, std::size_t MDIM>
	class CompositeFunctorInterface : public FunctionFunctorInterface<T,NDIM> {

		typedef Vector<double, NDIM> coordT; ///< Type of vector holding coordinates
		typedef FunctionImpl<T,NDIM> implT;
		typedef FunctionImpl<T,MDIM> implL;
		typedef std::shared_ptr<implT> pimplT;
		typedef std::shared_ptr<implL> pimplL;

		World& world;

	public:
		/// various MRA functions of NDIM dimensionality
		std::shared_ptr<implT> impl_ket;	///< supposedly the pair function
		std::shared_ptr<implT> impl_eri;	///< supposedly 1/r12

		/// various MRA functions of MDIM dimensionality (e.g. 3, if NDIM==6)
		std::shared_ptr<implL> impl_m1;	///< supposedly 1/r1
		std::shared_ptr<implL> impl_m2;	///< supposedly 1/r2
		std::shared_ptr<implL> impl_p1;	///< supposedly orbital 1
		std::shared_ptr<implL> impl_p2;	///< supposedly orbital 2

	public:

		/// constructor takes its Factory
		CompositeFunctorInterface(World& world, pimplT ket, pimplT g12,
				pimplL v1, pimplL v2, pimplL p1, pimplL p2)
			: world(world), impl_ket(ket), impl_eri(g12)
			, impl_m1(v1), impl_m2(v2), impl_p1(p1), impl_p2(p2)
		{

			// some consistency checks
			// either a pair ket is provided, or two particles (tba)
			MADNESS_ASSERT(impl_ket or (impl_p1 and impl_p2));

			// prepare base functions that make this function
			if (impl_ket and (not impl_ket->is_on_demand())) impl_ket->make_redundant(false);
			if (impl_eri) {
				if (not impl_eri->is_on_demand()) impl_eri->make_redundant(false);
			}
			if (impl_m1 and (not impl_m1->is_on_demand())) impl_m1->make_redundant(false);
			if (impl_m2 and (not impl_m2->is_on_demand())) impl_m2->make_redundant(false);

			if (impl_p1 and (not impl_p1->is_on_demand())) impl_p1->make_redundant(false);
			if (impl_p2 and (not impl_p2->is_on_demand())) impl_p2->make_redundant(false);
			world.gop.fence();

		}

		/// return value at point x; fairly inefficient
		T operator()(const coordT& x) const {
			print("there is no operator()(coordT&) in CompositeFunctorInterface, for good reason");
			MADNESS_ASSERT(0);
			return T(0);
		};

		bool provides_coeff() const {
			return false;
		}

	};


	/// ElementaryInterface (formerly FunctorInterfaceWrapper) interfaces a c-function

	/// hard-code your favorite function and interface it with this; Does only
	/// provide function values, no MRA coefficients. Care must be taken if the
	/// function we refer to is a singular function, and a on-demand function
	/// at the same time, since direct computation of coefficients via mraimpl::project
	/// might suffer from inaccurate quadrature.
	template<typename T, std::size_t NDIM>
	class ElementaryInterface : public FunctionFunctorInterface<T,NDIM> {

	public:
		typedef Vector<double, NDIM> coordT; ///< Type of vector holding coordinates
        typedef GenTensor<T> coeffT;

		T (*f)(const coordT&);

		ElementaryInterface(T (*f)(const coordT&)) : f(f) {}

		T operator()(const coordT& x) const {return f(x);}

		coeffT values(const Key<NDIM>& key, const Tensor<double>& quad_x) const {
	        typedef Tensor<T> tensorT;
            tensorT fval=madness::fcube(key,f,quad_x);
            return coeffT(fval,FunctionDefaults<NDIM>::get_thresh(),TT_FULL);
		}
	};

    /// FunctorInterface interfaces a class or struct with an operator()()
    template<typename T, std::size_t NDIM, typename opT>
    class FunctorInterface : public FunctionFunctorInterface<T,NDIM> {

    public:
        typedef Vector<double, NDIM> coordT; ///< Type of vector holding coordinates
        typedef GenTensor<T> coeffT;

        opT op;

        FunctorInterface(const opT& op) : op(op) {}

        T operator()(const coordT& x) const {return op(x);}
    };

	/// FunctionInterface implements a wrapper around any class with the operator()()
	template<typename T, size_t NDIM, typename opT>
	class FunctionInterface : public FunctionFunctorInterface<T,NDIM> {

	    typedef GenTensor<T> coeffT;
        typedef Vector<double, NDIM> coordT; ///< Type of vector holding coordinates

        const opT op;

    public:
        FunctionInterface(const opT& op) : op(op) {}

        T operator()(const coordT& coord) const {return op(coord);}

        bool provides_coeff() const {return false;}

	};

	/// base class to compute the wavelet coefficients for an isotropic 2e-operator

	/// all classes that derive from this base class use the Gaussian fitting
	/// procedure that has been developed for the BSH operator. We simply
	/// reuse the wavelet coefficients that we get from there to avoid
	/// evaluating the functions themselves, since the quadrature of singular
	/// functions is imprecise and slow.
	template<typename T, std::size_t NDIM>
    class TwoElectronInterface : public FunctionFunctorInterface<T,NDIM> {
	public:

		typedef GenTensor<T> coeffT;

		/// constructor: cf the Coulomb kernel

		/// @param[in]	lo		the smallest length scale to be resolved
		/// @param[in]	eps		the accuracy threshold
		TwoElectronInterface(double lo, double eps,
				const BoundaryConditions<6>& bc=FunctionDefaults<6>::get_bc(),
				int kk=FunctionDefaults<6>::get_k())
				:rank(), k(kk), lo(lo), hi(1.0) {

			// Presently we must have periodic or non-periodic in all dimensions.
			for (std::size_t d=1; d<6; ++d) {MADNESS_ASSERT(bc(d,0)==bc(0,0));}

			const Tensor<double>& width = FunctionDefaults<6>::get_cell_width();
			hi = width.normf(); // Diagonal width of cell
			if (bc(0,0) == BC_PERIODIC) hi *= 100; // Extend range for periodic summation

		}

		bool provides_coeff() const {
			return true;
		}

		/// return the coefficients of the function in 6D (x1,y1,z1, x2,y2,z2)
		coeffT coeff(const Key<NDIM>& key) const {
			Tensor<double> c=make_coeff(key);
            return coeffT(map_coeff(c),FunctionDefaults<6>::get_thresh(),TT_FULL);
		}

		T operator()(const Vector<double, NDIM>& x) const {
			print("there is no operator()(coordT&) in TwoElectronInterface, for good reason");
			MADNESS_ASSERT(0);
			return T(0);
		}

	protected:

		/// make the coefficients from the 1d convolution
		Tensor<double> make_coeff(const Key<6>& key) const {
			const Level n=key.level();
			const Vector<Translation,6> l=key.translation();

			// get the displacements for all 3 dimensions: x12, y12, z12
			const Translation l0=(l[0]-l[3]);
			const Translation l1=(l[1]-l[4]);
			const Translation l2=(l[2]-l[5]);

			Tensor<double> scr1(rank,k*k), scr2(rank,k*k,k*k);

			// lump all the terms together
			for (long mu=0; mu<rank; mu++) {
				const Tensor<double> r0=(ops[mu].getop(0)->rnlij(n,l0)).reshape(k*k);
				const Tensor<double> r1=(ops[mu].getop(1)->rnlij(n,l1)).reshape(k*k);
				const Tensor<double> r2=(ops[mu].getop(2)->rnlij(n,l2)).reshape(k*k);

				// include weights in first vector
				scr1(mu,Slice(_))=r0*ops[mu].getfac();

				// merge second and third vector to scr(r,k1,k2)
				scr2(mu,Slice(_),Slice(_))=outer(r1,r2);
			}

			Tensor<double> c=inner(scr1,scr2,0,0);
			return c;
		}

		/// the dimensions are a bit confused (x1,x2, y1,y2, z1,z2) -> (x1,y1,z1, x2,y2,z2)
		Tensor<double> map_coeff(const Tensor<double>& c) const {
			std::vector<long> map(6);
			map[0]=0;	map[1]=3;	map[2]=1;
			map[3]=4;	map[4]=2;	map[5]=5;
			return copy(c.reshape(k,k,k,k,k,k).mapdim(map));
		}

		/// initialize the Gaussian fit; uses the virtual function fit() to fit
		void initialize(const double eps) {
			GFit<double,3> fit=this->fit(eps);
			Tensor<double> coeff=fit.coeffs();
			Tensor<double> expnt=fit.exponents();

			// set some parameters
			rank=coeff.dim(0);
			ops.resize(rank);
			const Tensor<double>& width = FunctionDefaults<6>::get_cell_width();

			// construct all the terms
			for (int mu=0; mu<rank; ++mu) {
				//                double c = std::pow(sqrt(expnt(mu)/pi),static_cast<int>(NDIM)); // Normalization coeff
				double c = std::pow(sqrt(expnt(mu)/constants::pi),3); // Normalization coeff

				// We cache the normalized operator so the factor is the value we must multiply
				// by to recover the coeff we want.
				ops[mu].setfac(coeff(mu)/c);

				// only 3 dimensions here!
				for (std::size_t d=0; d<3; ++d) {
					ops[mu].setop(d,GaussianConvolution1DCache<double>::get(k, expnt(mu)*width[d]*width[d], 0, false));
				}
			}
		}

		/// derived classes must implement this -- cf GFit.h
		virtual GFit<double,3> fit(const double eps) const = 0;

		/// storing the coefficients
		mutable std::vector< ConvolutionND<double,6> > ops;

		/// the number of terms in the Gaussian quadrature
		int rank;

		/// the wavelet order
		int k;

		/// the smallest length scale that needs to be represented
		double lo;

		/// the largest length scale that needs to be represented
		double hi;

	};

	/// a function like f(x)=1/x
	class ElectronRepulsionInterface : public TwoElectronInterface<double,6> {
	public:

		/// constructor: cf the Coulomb kernel

		/// @param[in]	lo		the smallest length scale to be resolved
		/// @param[in]	eps		the accuracy threshold
		ElectronRepulsionInterface(double lo,double eps,
				const BoundaryConditions<6>& bc=FunctionDefaults<6>::get_bc(),
				int kk=FunctionDefaults<6>::get_k())
		  : TwoElectronInterface<double,6>(lo,eps,bc,kk) {

			initialize(eps);
		}

	private:

		GFit<double,3> fit(const double eps) const {
			return GFit<double,3>::CoulombFit(lo,hi,eps,false);
		}
	};

	/// a function like f(x) = exp(-mu x)/x
	class BSHFunctionInterface : public TwoElectronInterface<double,6> {
	public:

		/// constructor: cf the Coulomb kernel

      /// @param[in]	mu		the exponent of the BSH/inverse Laplacian
		/// @param[in]	lo		the smallest length scale to be resolved
		/// @param[in]	eps		the accuracy threshold
		BSHFunctionInterface(double mu, double lo, double eps,
				const BoundaryConditions<6>& bc=FunctionDefaults<6>::get_bc(),
				int kk=FunctionDefaults<6>::get_k())
		  : TwoElectronInterface<double,6>(lo,eps,bc,kk), mu(mu) {

			initialize(eps);
		}

	private:

		double mu;

		GFit<double,3> fit(const double eps) const {
			return GFit<double,3>::BSHFit(mu,lo,hi,eps,false);
		}
	};

	/// a function like f(x)=exp(-mu x)
	class SlaterFunctionInterface : public TwoElectronInterface<double,6> {
	public:

		/// constructor: cf the Coulomb kernel

		/// @param[in]	mu		the exponent of the Slater function
		/// @param[in]	lo		the smallest length scale to be resolved
		/// @param[in]	eps		the accuracy threshold
		SlaterFunctionInterface(double mu, double lo, double eps,
				const BoundaryConditions<6>& bc=FunctionDefaults<6>::get_bc(),
				int kk=FunctionDefaults<6>::get_k())
		  : TwoElectronInterface<double,6>(lo,eps,bc,kk), mu(mu) {
			initialize(eps);
		}

	private:

		double mu;

		GFit<double,3> fit(const double eps) const {
			return GFit<double,3>::SlaterFit(mu,lo,hi,eps,false);
		}
	};

	/// a function like f(x) = (1 - exp(-mu x))/(2 gamma)
	class SlaterF12Interface : public TwoElectronInterface<double,6> {
	public:

		/// constructor: cf the Coulomb kernel

		/// @param[in]	mu		the exponent of the Slater function
		/// @param[in]	lo		the smallest length scale to be resolved
		/// @param[in]	eps		the accuracy threshold
		SlaterF12Interface(double mu, double lo, double eps,
				const BoundaryConditions<6>& bc=FunctionDefaults<6>::get_bc(),
				int kk=FunctionDefaults<6>::get_k())
		  : TwoElectronInterface<double,6>(lo,eps,bc,kk), mu(mu) {

			initialize(eps);
		}

		/// overload the function of the base class
		coeffT coeff(const Key<6>& key) const {

			Tensor<double> c=make_coeff(key);

			// subtract 1 from the (0,0,..,0) element of the tensor,
			// which is the 0th order polynomial coefficient
        	double one_coeff1=1.0*sqrt(FunctionDefaults<6>::get_cell_volume())
        			*pow(0.5,0.5*6*key.level());
            std::vector<long> v0(6,0L);
            c(v0)-=one_coeff1;

			c.scale(-0.5/mu);
            return coeffT(map_coeff(c),FunctionDefaults<6>::get_thresh(),TT_FULL);
		}

	private:

		double mu;

		GFit<double,3> fit(const double eps) const {
			return GFit<double,3>::SlaterFit(mu,lo,hi,eps,false);
		}
	};

// Not right
//	/// a function like f(x) = (1 - exp(-mu x))/x
//	class FGInterface : public TwoElectronInterface<double,6> {
//	public:
//
//		/// constructor: cf the Coulomb kernel
//
//		/// @param[in]	mu		the exponent of the Slater function
//		/// @param[in]	lo		the smallest length scale to be resolved
//		/// @param[in]	eps		the accuracy threshold
//		FGInterface(double mu, double lo, double eps,
//				const BoundaryConditions<6>& bc=FunctionDefaults<6>::get_bc(),
//				int kk=FunctionDefaults<6>::get_k())
//		  : TwoElectronInterface<double,6>(lo,eps,bc,kk), mu(mu) {
//
//			initialize(eps);
//		}
//
//	private:
//
//		double mu;
//
//		GFit<double,3> fit(const double eps) const {
//			return GFit<double,3>::SlaterFit(mu,lo,hi,eps,false);
//		}
//	};


#if 0

	/// ElectronRepulsionInterface implements the electron repulsion term 1/r12

	/// this is essentially just a wrapper around ElectronRepulsion
	template<typename T, std::size_t NDIM>
	class ElectronRepulsionInterface : public FunctionFunctorInterface<T,NDIM> {

		typedef GenTensor<T> coeffT;
		typedef Vector<double, NDIM> coordT; ///< Type of vector holding coordinates

		/// the class computing the coefficients
		ElectronRepulsion eri;

	public:

		/// constructor takes the same parameters as the Coulomb operator
		/// which it uses to compute the coefficients
		ElectronRepulsionInterface(World& world,double lo,double eps,
                const BoundaryConditions<NDIM>& bc=FunctionDefaults<NDIM>::get_bc(),
                int k=FunctionDefaults<NDIM>::get_k())
			: eri(ElectronRepulsion(eps,eps,bc,k)) {
		}


		/// return value at point x; fairly inefficient
		T operator()(const coordT& x) const {
			print("there is no operator()(coordT&) in ElectronRepulsionInterface, for good reason");
			MADNESS_ASSERT(0);
			return T(0);
		};


		/// return sum coefficients for imagined node at key
		coeffT coeff(const Key<NDIM>& key) const {
            return coeffT(this->eri.coeff(key),FunctionDefaults<NDIM>::get_thresh(),
                    TT_FULL);
		}

	};

	/// FGIntegralInterface implements the two-electron integral (1-exp(-gamma*r12))/r12

	/// this is essentially just a wrapper around ElectronRepulsion
	/// The integral expressed as:   1/r12 - exp(-gamma*r12)/r12
	/// which can be expressed with an eri and a bsh
	template<typename T, std::size_t NDIM>
	class FGIntegralInterface : public FunctionFunctorInterface<T,NDIM> {

		typedef GenTensor<T> coeffT;
		typedef Vector<double, NDIM> coordT; ///< Type of vector holding coordinates

		/// the class computing the coefficients
		ElectronRepulsion eri;
		BSHFunction bsh;

	public:

		/// constructor takes the same parameters as the Coulomb operator
		/// which it uses to compute the coefficients
		FGIntegralInterface(World& world, double lo, double eps, double gamma,
                const BoundaryConditions<NDIM>& bc=FunctionDefaults<NDIM>::get_bc(),
                int k=FunctionDefaults<NDIM>::get_k())
			: eri(ElectronRepulsion(eps,eps,0.0,bc,k))
			, bsh(BSHFunction(eps,eps,gamma,bc,k)) {
		}

		bool provides_coeff() const {
			return true;
		}

		/// return value at point x; fairly inefficient
		T operator()(const coordT& x) const {
			print("there is no operator()(coordT&) in FGIntegralInterface, for good reason");
			MADNESS_ASSERT(0);
			return T(0);
		};

		/// return sum coefficients for imagined node at key
		coeffT coeff(const Key<NDIM>& key) const {
	        typedef Tensor<T> tensorT;
			tensorT e_b=eri.coeff(key)-bsh.coeff(key);
            return coeffT(e_b,FunctionDefaults<NDIM>::get_thresh(),TT_FULL);
		}

	};

#endif

}

#endif // MADNESS_MRA_FUNCTION_INTERFACE_H__INCLUDED
