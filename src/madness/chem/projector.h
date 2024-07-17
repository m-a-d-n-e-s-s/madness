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


    class ProjectorBase {
    protected:
        /// a projector might work only on a subset of dimensions, e.g. P(1) | \psi(1,2) >
        int particle=-1;        // must only be 0 or 1!
    public:
        virtual ~ProjectorBase() {}
        virtual void set_particle(const int p)
        {
            MADNESS_CHECK_THROW(p==0 or p==1, "particle must be 0 or 1");
            particle=p;
        }
        virtual int get_particle() const {return particle;}
        virtual std::string type() const = 0;
    };

    template<typename T, std::size_t NDIM>
    class CCPairFunction;

    template<typename T, std::size_t NDIM>
    std::vector<CCPairFunction<T,NDIM>> apply(const ProjectorBase& P, const std::vector<CCPairFunction<T,NDIM>>& argument);

    /// simple projector class

    /// use this class to project a function or a set of functions on
    /// another space of function. The projector can handle different sets of
    /// functions for the bra and the ket space, e.g. in case of regularized
    /// orbitals: |f>  <->  <f|R^2
    template<typename T, std::size_t NDIM>
    class Projector : public ProjectorBase {

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
                , mo_bra_(vecfuncT(1,bra)) {
            MADNESS_CHECK_THROW(mo_bra_.size()==mo_ket_.size(), "bra and ket spaces must have the same size in projector");
        }

        /// constructor with a set of orbitals to project out

        /// bra and ket spaces are symmetric
        Projector(const vecfuncT& p) : mo_ket_(p), mo_bra_(p) {
            MADNESS_CHECK_THROW(mo_bra_.size()==mo_ket_.size(), "bra and ket spaces must have the same size in projector");
        }

        /// constructor with a set of orbitals to project out

        /// bra and ket spaces are not symmetric (e.g. |ket>^+ = <bra|R2 )
        Projector(const vecfuncT& bra, const vecfuncT& ket) : mo_ket_(ket),
                mo_bra_(bra) {
            MADNESS_CHECK_THROW(mo_bra_.size()==mo_ket_.size(), "bra and ket spaces must have the same size in projector");
        }

    	void set_spaces(const vecfuncT& p) {
            mo_bra_=p;
            mo_ket_=p;
        }

        void set_spaces(const vecfuncT& bra, const vecfuncT& ket) {
            mo_bra_=bra;
            mo_ket_=ket;
            MADNESS_CHECK_THROW(mo_bra_.size()==mo_ket_.size(), "bra and ket spaces must have the same size in projector");
        }

        virtual std::string type() const override {return "PProjector";}

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
            if (f.size()==0) return vecfuncT();
            World& world=f[0].world();
            Tensor<T> ovlp=matrix_inner(world,mo_bra_,f);
            vecfuncT result=transform(world,mo_ket_,ovlp,true);
            truncate(world,result);
            return result;
        }

        /// apply 3D Projector to one particle of a 6D function
        /// \f[
        /// |result> = \sum_p |p(particle)> <p(particle)|f(1,2)>_{particle}
        /// \f]
        /// @param[in] f the 6D function to be projected
        /// @param[in] the particle that is projected (0 or 1)
        /// @return the projected function
        template<std::size_t KDIM>
        typename std::enable_if<KDIM==2*NDIM, Function<T,KDIM> >::type
        operator()(const Function<T,KDIM>& f, int particle1=-1) const {
            if (particle1==-1) particle1=get_particle();
            MADNESS_CHECK_THROW(particle1 == 0 or particle1 == 1, "particle must be 0 or 1");
            auto [left,right]=get_vectors_for_outer_product(f);
            return hartree_product(left,right);
        }

        /// apply the projection parts of the operator on a function f

        /// The operator applied on f(1,2) is
        ///  O(1)f(1,2) = \sum_i |i(1) > <i(1) | f(1,2)>_1 = \sum_i |i(1) f_i(2)>
        /// return the lo-dim vectors i and f_i only, perform no outer product
        std::pair<std::vector<Function<T,NDIM>>,std::vector<Function<T,NDIM>>>
        get_vectors_for_outer_product(const Function<T,2*NDIM>& f) const {
            World& world=f.world();
            reconstruct(world, mo_bra_, false);
            f.reconstruct(false);
            reconstruct(world, mo_ket_, true);
            std::vector<Function<T,NDIM>> projected;
            for (const auto& i : mo_bra_) {
                projected.push_back(f.project_out(i,particle));
            }
            if (particle==0) return std::make_pair(mo_ket_,projected);
            else if (particle==1) return std::make_pair(projected,mo_ket_);
            else {
                MADNESS_EXCEPTION("confused particles in Projector::get_vector_for_outer_products",1);
            }
        }

        template<typename argT>
        typename std::enable_if<!std::is_same<argT,Function<T,2*NDIM> >::value, argT>::type
        operator()(const argT& argument) const {
            return madness::apply(*this,argument);
        }

        vecfuncT get_bra_vector() const {return mo_bra_;}

        vecfuncT get_ket_vector() const {return mo_ket_;}

    };


    /// orthogonality projector

    /// projects out the space given in the constructor
    /// \f[
    ///   |result> = |f> - \sum_p |p><p|f>
    /// \f]
    template<typename T, std::size_t NDIM>
    class QProjector : public ProjectorBase {
        typedef std::vector<Function<T,NDIM> > vecfuncT;

    public:

        /// default ctor
        QProjector() = default;

        /// constructor with symmetric bra and ket spaces
        [[deprecated]] QProjector(World& world, const vecfuncT& amo) : O(amo) {};

        /// constructor with asymmetric bra and ket spaces
        [[deprecated]] QProjector(World& world, const vecfuncT& bra, const vecfuncT& ket)
            : O(bra,ket) {};

        /// constructor with symmetric bra and ket spaces
        QProjector(const vecfuncT& amo) : O(amo) {};

        /// constructor with asymmetric bra and ket spaces
        QProjector(const vecfuncT& bra, const vecfuncT& ket)
            : O(bra,ket) {};

        /// copy ctor
        QProjector(const QProjector& other) = default;

        std::string type() const override {return "QProjector";}

        void set_spaces(const vecfuncT& p) {
            O.set_spaces(p);
        }

        void set_spaces(const vecfuncT& bra, const vecfuncT& ket) {
            O.set_spaces(bra,ket);
        }

        Function<T,NDIM> operator()(const Function<T,NDIM>& rhs) const {
            return (rhs-O(rhs)).truncate();
        }

        vecfuncT operator()(const vecfuncT& rhs) const {
        	if (rhs.size()==0) return vecfuncT();
            vecfuncT result=rhs-O(rhs);
            truncate(result[0].world(),result);
            return result;
        }

        Function<T,2*NDIM> operator()(const Function<T,2*NDIM>& f, const size_t particle=-1) const {
            return f-O(f,particle);
        }

        template<typename argT>
        argT operator()(const argT& argument) const {
            return madness::apply(*this,argument);
        }

        vecfuncT get_bra_vector() const {return O.get_bra_vector();}

        vecfuncT get_ket_vector() const {return O.get_ket_vector();}

        Projector<T,NDIM> get_P_projector() const {return O;}

        void set_particle(const int p) override {
            O.set_particle(p);
            particle=p;
        }

        int get_particle() const override {
            return O.get_particle();
        }

    private:
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
    class StrongOrthogonalityProjector : public ProjectorBase {

    	typedef std::vector<Function<T,NDIM> > vecfuncT;

    public:

    	/// default ctor
    	StrongOrthogonalityProjector(World& world) : world(world) {}

        std::string type() const override {return "SOProjector";}

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
            MADNESS_CHECK_THROW(ket1.size()==bra1.size(), "bra1 and ket1 spaces must have the same size in SOprojector");
            MADNESS_CHECK_THROW(ket2.size()==bra2.size(), "bra2 and ket2 spaces must have the same size in SOprojector");
    	}

    	/// return the orbital space for the ket of particle 1
    	vecfuncT ket1() const {return ket1_;}

    	/// return the orbital space for the bra of particle 1
    	vecfuncT bra1() const {return bra1_;}

    	/// return the orbital space for the ket of particle 2
    	vecfuncT ket2() const {return ket2_;}

    	/// return the orbital space for the bra of particle 2
    	vecfuncT bra2() const {return bra2_;}

        void set_particle(const int p) override {
            MADNESS_EXCEPTION("You cannot set a particle in the SO projector",1);
        }

        template<typename argT>
        argT operator()(const argT& argument) const {
            return madness::apply(*this,argument);
        }

        /// apply the projection parts of the strong orthogonality operator Q12 on a function f

        /// The SO operator is defined as 1-O1-O2+O1O2, where O1 and O2 are projectors
        /// return the term -O1-O2+O1O2 only, such that Q12 f = 1 + outer(result.first,result.second)
        std::pair<std::vector<Function<T,NDIM>>,std::vector<Function<T,NDIM>>>
        get_vectors_for_outer_product(const Function<T,2*NDIM>& f) const {
            // Eq. (A9): g_kl = < k(1) l(2) | f(1,2) >
            // note no (kl) symmetry here!
            reconstruct(world, bra1_, false);
            reconstruct(world, bra2_, true);
            Tensor<double> g_kl(bra1_.size(), bra2_.size());
            for (size_t k = 0; k < bra1_.size(); ++k) {
                for (size_t l = 0; l < bra2_.size(); ++l) {
                    Function<T, 2 * NDIM> kl = CompositeFactory<T, 2 * NDIM, NDIM>(world)
                            .particle1(bra1_[k]).particle2(bra2_[l]);
                    g_kl(k, l) = inner(f, kl);
                }
            }

            // Eq. (A12)
            // project out the mainly first particle: O1 (1 - 1/2 O2)
            std::vector<Function<T, NDIM>> h2(bra1_.size());
            std::vector<Function<T, NDIM>> h1(ket2_.size());
            reconstruct(world, bra1_, false);
            reconstruct(world, bra2_, true);
            for (size_t k = 0; k < bra1_.size(); ++k) {

                // project out the mainly first particle: O1 (1 - 1/2 O2): first term
                // Eq. (A10)
                h2[k] = f.project_out(bra1_[k], 0);

                // project out the mainly second particle: O2 (1 - 1/2 O1): first term
                // Eq. (A11)
                std::size_t l = k;
                h1[l] = f.project_out(bra2_[l], 1);

            }

            // Transforms a vector of functions according to new[i] = sum[j] old[j]*c[j,i]
            // project out the mainly first particle: O1 (1 - 1/2 O2): second term
            // Eq. (A12), (A13)
            h2 -= transform(world, ket2_, 0.5 * transpose(g_kl), false);   // ordering g(k,l) is correct
            h1 -= transform(world, ket1_, 0.5 * g_kl, true);               // ordering g(k,l) is correct
            // aka
            // 	for (size_t l=0; l<ket2_.size(); ++l) {
            // 		h2[k]-=0.5*g_kl(k,l)*ket2_[l];
            // 	}

            change_tree_state(h1, reconstructed, false);
            change_tree_state(h2, reconstructed, false);
            change_tree_state(ket1_, reconstructed, false);
            change_tree_state(ket2_, reconstructed, false);
            world.gop.fence();

    	    auto left=append(ket1_,h1);
    	    auto right=append(h2,ket2_);
    	    return std::make_pair(-1.0*left,right);
        }

        /// apply the strong orthogonality operator Q12 on a function f

    	/// notation of the equations follows
    	/// J. Chem. Phys., vol. 139, no. 11, p. 114106, 2013.
        Function<T,2*NDIM> operator()(const Function<T,2*NDIM>& f) const {

            // simple and it works for higher accuracies, but might be
            // imprecise for lower accuracies
//        	return (f-O1(f)-O2(f)+O1(O2(f))).truncate().reduce_rank();

    	    auto [left,right]=get_vectors_for_outer_product(f);

    	    // temporarily tighten the threshold
            double thresh=FunctionDefaults<2*NDIM>::get_thresh();
            double tight_thresh=thresh*0.1;
            FunctionDefaults<2*NDIM>::set_thresh(tight_thresh);

    	    auto tmp=hartree_product(left,right);
    	    tmp.truncate(thresh*0.3);

            FunctionDefaults<2*NDIM>::set_thresh(thresh);
            Function<T, 2 * NDIM> result = copy(f);
    	    result+=tmp;
            result.truncate(thresh).reduce_rank();
            return result;
        }

    private:

        /// the world
        World& world;

        /// the spaces of the projector of particles 1 and 2
        std::vector<Function<T,NDIM> > ket1_, bra1_, ket2_, bra2_;

    };


    /// an outer product of two projectors
    template<typename projT, typename projQ>
    class OuterProjector : public ProjectorBase {
        projT projector0;
        projQ projector1;
    public:

        OuterProjector() = default;
        OuterProjector(const projT& p0, const projQ& p1) : projector0(p0), projector1(p1) {
            static_assert(std::is_base_of<ProjectorBase,projT>::value, "projT must be a ProjectorBase");
            static_assert(std::is_base_of<ProjectorBase,projQ>::value, "projQ must be a ProjectorBase");
            projector0.set_particle(0);
            projector1.set_particle(1);
        }

        std::string type() const override {
            return "OuterProjector";
        }

        template<typename resultT>
        resultT operator()(const resultT& argument) const {

            if (projector0.type()=="PProjector") return projector1(projector0(argument));
            return projector0(projector1(argument));
        }
    };

//    template<typename projT, typename projQ>
//    OuterProjector<projT, projQ> outer(const projT& p0 , const projQ& p1) {
//        return OuterProjector<projT, projQ>(p0, p1);
//    }

    template<typename projT, typename projQ>
    typename std::enable_if<std::is_base_of<ProjectorBase,projT>::value, OuterProjector<projT,projQ>>::type
    outer(const projT& p0 , const projQ& p1) {
        return OuterProjector<projT, projQ>(p0, p1);
    }
}

#endif /* PROJECTOR_H_ */
