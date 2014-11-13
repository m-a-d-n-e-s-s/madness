/*
 * projector.h
 *
 *  Created on: Jan 24, 2014
 *      Author: fbischoff
 */

#ifndef MADNESS_CHEM_PROJECTOR_H__INCLUDED
#define MADNESS_CHEM_PROJECTOR_H__INCLUDED

#include <madness/mra/mra.h>

namespace madness {

    /// simple projector class for 1- and 2-particle projectors
    template<typename T, std::size_t NDIM>
    class Projector {

        int particle_;
        std::vector<Function<T,NDIM> > p_;

    public:

        Projector() : particle_(0), p_(std::vector<Function<T,NDIM> >()) {}

        /// simple constructor with only one orbital to project out
        Projector(const Function<T,NDIM>& p, const int particle=0)
        	: particle_(particle), p_(std::vector<Function<T,NDIM> >(1,p)) {
            MADNESS_ASSERT(particle_==0 or particle_==1);
        }

        /// constructor with a set of orbitals to project out
        Projector(const std::vector<Function<T,NDIM> >& p, const int particle=0)
        	: particle_(particle), p_(p) {
            MADNESS_ASSERT(particle_==0 or particle_==1);
        }

        int& particle() {return particle_;}
        const int& particle() const {return particle_;}

        /// get a const reference to the orbitals
        const std::vector<Function<T,NDIM> >& p() const {return p_;}

        /// project f on p: \f$ \left| result \right> =  \left| p \right>\left< p \right| f> \f$
        template<std::size_t FDIM>
        typename enable_if_c<NDIM==FDIM, Function<T,FDIM> >::type
        operator()(const Function<T,FDIM>& f) const {

        	World& world=f.world();

        	Function<T,NDIM> sum=FunctionFactory<T,NDIM>(f.world());

            compress(world, p_,false);	// don't fence
        	sum.compress(false);
            f.compress();				// fence

        	// the overlap of all orbitals with the rhs
        	Tensor<double> ovlp=inner(world,f,p_);

            for (std::size_t i=0; i<p_.size(); ++i) {
            	if (ovlp(i) != T(0.0)) sum.gaxpy(1.0,p_[i],ovlp(i),false);
            }
            world.gop.fence();
            sum.truncate();

            return sum;
        }

        /// project p out of f: \f$ \left| result(1,2) \right> = sum_p \left| p(1) \right> \left< p(1) \right| \left. f(1,2) \right> \f$
        template<std::size_t FDIM>
        typename enable_if_c<2*NDIM==FDIM, Function<T,FDIM> >::type
        operator()(const Function<T,FDIM>& f) const {
            real_function_6d sum=real_factory_6d(f.world());
            for (unsigned int i=0; i<p_.size(); ++i) {
                const real_function_3d pf2=f.projfirect_out(p_[i],particle_);
                real_function_6d tmp;
                MADNESS_EXCEPTION("Projector class: the hartree product is inaccurate -- don't use it",1);
                if (particle_==0) tmp=hartree_product(p_[i],pf2);
                else tmp=hartree_product(pf2,p_[i]);
                sum=(sum+tmp);
            }
            sum.truncate();
            return sum;
        }
    };

    /// a SO projector class

    /// The SO projector is defined as
    ///  Q12 = (1-O1)(1-O2),
    ///  O = \sum_i | i >< i |
    /// where O1 and O2 are projectors for electron 1 and 2 on the occupied space
    /// As a special case there might be a similarity transformed occupied
    /// space, resulting in different bras and kets
    template<typename T, std::size_t NDIM>
    class StrongOrthogonalityProjector {

    	typedef std::vector<Function<T,NDIM> > vecfuncT;

    public:

    	/// default ctor
    	StrongOrthogonalityProjector(World& world) : world(world) {}

    	/// set the same spaces for the projectors for particle 1 and 2
    	void set_spaces(const vecfuncT& p) {
    		ket1_=p;
    		bra1_=p;
    		ket2_=p;
    		bra2_=p;
    	}

    	/// set different spaces for the projectors for particle 1 and 2

    	/// the SO projector is
    	///  Q12 = (1 - O1) (1 - O2)
    	///  O1 = \sum_i | ket1_i >< bra1_i |
    	///  O2 = \sum_i | ket2_i >< bra2_i |
    	/// this case occurs for a similarity transformed Q12
    	void set_spaces(const vecfuncT& bra1, const vecfuncT& ket1,
        		const vecfuncT& bra2, const vecfuncT& ket2) {
    		ket1_=ket1;
    		bra1_=bra1;
    		ket2_=ket2;
    		bra2_=bra2;
    	}

    	/// return the orbital space for the ket of particle 1
    	vecfuncT ket1() const {return ket1_;}

    	/// return the orbital space for the bra of particle 1
    	vecfuncT bra1() const {return bra1_;}

    	/// return the orbital space for the ket of particle 2
    	vecfuncT ket2() const {return ket2_;}

    	/// return the orbital space for the bra of particle 2
    	vecfuncT bra2() const {return bra2_;}

        /// apply the strong orthogonality operator Q12 on a function f

    	/// notation of the equations follows
    	/// J. Chem. Phys., vol. 139, no. 11, p. 114106, 2013.
        Function<T,2*NDIM> operator()(const Function<T,2*NDIM>& f) const {

        	// simple and it works for higher accuracies, but might be
        	// imprecise for lower accuracies
//        	return (f-O1(f)-O2(f)+O1(O2(f))).truncate().reduce_rank();

        	const double thresh=FunctionDefaults<2*NDIM>::get_thresh();
        	const double tight_thresh=FunctionDefaults<2*NDIM>::get_thresh()*0.1;

        	// Eq. (A9): g_kl = < k(1) l(2) | f(1,2) >
        	// note no (kl) symmetry here!
        	Tensor<double> g_kl(bra1_.size(),bra2_.size());
        	for (size_t k=0; k<bra1_.size(); ++k) {
        		for (size_t l=0; l<bra2_.size(); ++l) {
        			Function<T,2*NDIM> kl=CompositeFactory<T,2*NDIM,NDIM>(world)
    	            		.particle1(copy(bra1_[k])).particle2(copy(bra2_[l]));
        			g_kl(k,l)=inner(f,kl);
        		}
        	}
//        	if (world.rank()==0) {print(g_kl);};

        	// Eq. (A12)
        	// project out the mainly first particle: O1 (1 - 1/2 O2)
        	Function<T,2*NDIM> r1=FunctionFactory<T,2*NDIM>(world).thresh(tight_thresh);
        	for (size_t k=0; k<bra1_.size(); ++k) {

        		// Eq. (A10)
        		Function<T,NDIM> h2=f.project_out(bra1_[k],0);

        		// Eq. (A12)
            	for (size_t l=0; l<ket2_.size(); ++l) {
            		h2-=0.5*g_kl(k,l)*ket2_[l];
            	}

            	// Eq. (A7), second term rhs
            	// the hartree product tends to be inaccurate; tighten threshold
            	FunctionDefaults<2*NDIM>::set_thresh(tight_thresh);
            	r1=(r1+hartree_product(ket1_[k],h2));
            	FunctionDefaults<2*NDIM>::set_thresh(thresh);
            	r1.set_thresh(thresh);
            	r1.print_size("r1"+stringify(k));
        	}

        	// project out the mainly second particle: O2 (1 - 1/2 O1)
        	Function<T,2*NDIM> r2=FunctionFactory<T,2*NDIM>(world).thresh(tight_thresh);
        	for (size_t l=0; l<ket2_.size(); ++l) {

        		// Eq. (A11)
        		Function<T,NDIM> h1=f.project_out(bra2_[l],1);

        		// Eq. (A13)
            	for (size_t k=0; k<ket1_.size(); ++k) {
            		h1-=0.5*g_kl(k,l)*ket1_[k];			// ordering g(k,l) is correct
            	}

            	// Eq. (A7), third term rhs
            	// the hartree product tends to be inaccurate; tighten threshold
            	FunctionDefaults<2*NDIM>::set_thresh(tight_thresh);
            	r2=(r2+hartree_product(h1,ket2_[l]));
            	r2.set_thresh(thresh);
            	FunctionDefaults<2*NDIM>::set_thresh(thresh);
            	r2.print_size("r2"+stringify(l));
        	}
        	FunctionDefaults<2*NDIM>::set_thresh(tight_thresh);
        	Function<T,2*NDIM> result=(f-r1-r2).truncate().reduce_rank();
        	FunctionDefaults<2*NDIM>::set_thresh(thresh);

//        	// for debugging purposes only: check orthogonality
//        	for (size_t k=0; k<hf->nocc(); ++k) {
//        		for (size_t l=0; l<hf->nocc(); ++l) {
//    	            real_function_6d kl=CompositeFactory<double,6,3>(world)
//    	            		.particle1(copy(O1_mos[k])).particle2(copy(O2_mos[l]));
//        			g_kl(k,l)=inner(result,kl);
//        		}
//        	}
//        	if (world.rank()==0) {print(g_kl);};

        	return result;
        }

    private:

        /// the world
        World& world;

        /// the spaces of the projector of particles 1 and 2
        std::vector<Function<T,NDIM> > ket1_, bra1_, ket2_, bra2_;

    };
}

#endif /* PROJECTOR_H_ */
