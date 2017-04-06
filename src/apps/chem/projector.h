/*
 * projector.h
 *
 *  Created on: Jan 24, 2014
 *      Author: fbischoff
 */

#ifndef MADNESS_CHEM_PROJECTOR_H__INCLUDED
#define MADNESS_CHEM_PROJECTOR_H__INCLUDED

#include <madness/mra/mra.h>
#include <type_traits>

namespace madness {

    /// simple projector class

    /// use this class to project a function or a set of functions on
    /// another space of function. The projector can handle different sets of
    /// functions for the bra and the ket space, e.g. in case of regularized
    /// orbitals: |f>  <->  <f|R^2
    template<typename T, std::size_t NDIM>
    class Projector {

        typedef Function<T,NDIM> funcT;
        typedef std::vector<funcT> vecfuncT;

        /// the space onto which the test functions will be projected: |ket>
        std::vector<Function<T,NDIM> > mo_ket_;
        /// the dual space onto which the test functions will be projected: <bra|
        std::vector<Function<T,NDIM> > mo_bra_;

    public:

        Projector() : mo_ket_(vecfuncT()), mo_bra_(vecfuncT()) {}

        /// simple constructor with only one orbital to project

        /// bra and ket spaces are symmetric
        Projector(const Function<T,NDIM>& mo) : mo_ket_(vecfuncT(1,mo)),
                mo_bra_(vecfuncT(1,mo)) {}

        /// simple constructor with only one orbital to project

        /// bra and ket spaces are not symmetric (e.g. |ket>^+ = <bra|R2 )
        Projector(const funcT& bra, const funcT& ket) : mo_ket_(vecfuncT(1,ket))
                , mo_bra_(vecfuncT(1,bra)) {}

        /// constructor with a set of orbitals to project out

        /// bra and ket spaces are symmetric
        Projector(const vecfuncT& p) : mo_ket_(p), mo_bra_(p) {}

        /// constructor with a set of orbitals to project out

        /// bra and ket spaces are not symmetric (e.g. |ket>^+ = <bra|R2 )
        Projector(const vecfuncT& bra, const vecfuncT& ket) : mo_ket_(ket),
                mo_bra_(bra) {}

        /// project f on p:

        /// \f[
        ///     | result > =  \sum_p | p > <p|f>
        /// \f]
        /// @param[in]  f   the function to be projected
        /// @return     the projection of f on the space of p
        funcT operator()(const funcT& f) const {
            return this->operator()(vecfuncT(1,f)).front();
        }

        /// project f on p:

        /// \f[
        ///     | result > =  \sum_p | p > <p|f>
        /// \f]
        /// @param[in]  f   the vector of functions to be projected
        /// @return     the projection of f on the space of p
        vecfuncT operator()(const vecfuncT& f) const {
            MADNESS_ASSERT(f.size()>0);
            World& world=f[0].world();
            Tensor<double> ovlp=matrix_inner(world,mo_bra_,f);
            vecfuncT result=transform(world,mo_ket_,ovlp,true);
            truncate(world,result);
            return result;
        }

        /// apply 3D Projector to one particle of a 6D function
        /// \f[
        /// |result> = \sum_p |p(particle)> <p(particle)|f(1,2)>_{particle}
        /// \f]
        /// @param[in] f the 6D function to be projected
        /// @param[in] the particle that is projected (1 or 2)
        /// @return the projected function
        real_function_6d operator()(const real_function_6d&f, const size_t particle)const{
          real_function_6d result = real_factory_6d(f.world());
          MADNESS_ASSERT(particle==1 or particle==2);
          for(size_t i=0;i<mo_ket_.size();i++){
            real_function_3d tmp1 = mo_ket_[i];
            real_function_3d tmp2 = f.project_out(mo_bra_[i],particle-1);
            real_function_6d tmp12;
            if(particle==1){
              tmp12 = CompositeFactory<double,6,3>(f.world()).particle1(copy(tmp1)).particle2(copy(tmp2));
              tmp12.fill_tree();
            }else{
              tmp12 = CompositeFactory<double,6,3>(f.world()).particle1(copy(tmp2)).particle2(copy(tmp1));
              tmp12.fill_tree();
            }
            result += tmp12;
          }
          return result;
        }

    };


    /// orthogonality projector

    /// projects out the space given in the constructor
    /// \f[
    ///   |result> = |f> - \sum_p |p><p|f>
    /// \f]
    template<typename T, std::size_t NDIM>
    class QProjector {
        typedef std::vector<Function<T,NDIM> > vecfuncT;

    public:

        /// constructor with symmetric bra and ket spaces
        QProjector(World& world, const vecfuncT& amo) : world(world), O(amo) {};

        /// constructor with asymmetric bra and ket spaces
        QProjector(World& world, const vecfuncT& bra, const vecfuncT& ket)
            : world(world), O(bra,ket) {};

        Function<T,NDIM> operator()(const Function<T,NDIM>& rhs) const {
            return (rhs-O(rhs)).truncate();
        }

        vecfuncT operator()(const vecfuncT& rhs) const {
            vecfuncT result=sub(world,rhs,O(rhs));
            truncate(world,result);
            return result;
        }

    private:
        World& world;
        Projector<T,NDIM> O;
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
