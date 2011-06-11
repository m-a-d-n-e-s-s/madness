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

#ifndef MADNESS_MRA_FUNCIMPL_H__INCLUDED
#define MADNESS_MRA_FUNCIMPL_H__INCLUDED

/// \file funcimpl.h
/// \brief Provides FunctionCommonData, FunctionImpl and FunctionFactory

#include <iostream>
#include <world/world.h>
#include <world/print.h>
#include <misc/misc.h>
#include <tensor/tensor.h>
#include <mra/key.h>
#include <mra/funcdefaults.h>
#include "mra/flonode.h"
#include "mra/function_factory_and_interface.h"
//#include "mra/sepreptensor.h"
#include "mra/gentensor.h"
#include <world/typestuff.h>

namespace madness {
    template <typename T, std::size_t NDIM>
    class DerivativeBase;

    template<typename T, std::size_t NDIM>
    class FunctionImpl;

    template<typename T, std::size_t NDIM>
    class FunctionNode;

    template<typename T, std::size_t NDIM>
    class Function;

    template<typename T, std::size_t NDIM>
    class FunctionFactory;

    template<typename T, std::size_t NDIM, std::size_t MDIM>
    class CompositeFunctorInterface;

    template<int D>
    class LoadBalImpl;

    template<int D>
    class LBTree;

    template<int D>
    class MyPmap;
}

namespace madness {



#if HAVE_GENTENSOR
#  define __Tensor GenTensor
#else
#  define __Tensor Tensor
#endif
    /// A simple process map soon to be supplanted by Rebecca's
    template<typename keyT>
    class SimpleMap : public WorldDCPmapInterface<keyT> {
    private:
        const int nproc;
        const ProcessID me;
        const int n;

    public:
        SimpleMap(World& world, int n = 4) :
                nproc(world.nproc()), me(world.rank()), n(n) {
        }

        ProcessID
        owner(const keyT& key) const {
            if (key.level() == 0) {
                return 0;
            }
            else if (key.level() <= n) {
                return hash(key) % nproc;
            }
            else {
                return hash(key.parent(key.level() - n)) % nproc;
            }
        }
    };

    /// FunctionCommonData holds all Function data common for given k

    /// Since Function assignment and copy constructors are shallow it
    /// greatly simplifies maintaining consistent state to have all
    /// (permanent) state encapsulated in a single class.  The state
    /// is shared between instances using a shared_ptr.  Also,
    /// separating shared from instance specific state accelerates the
    /// constructor, which is important for massive parallelism, and
    /// permitting inexpensive use of temporaries.  The default copy
    /// constructor and assignment operator are used but are probably
    /// never invoked.
    template<typename T, std::size_t NDIM>
    class FunctionCommonData {
    private:
        static const FunctionCommonData<T, NDIM>* data[MAXK];

        /// Private.  Initialize the twoscale coefficients
        void
        _init_twoscale();

        /// Private.  Do first use initialization via get.
        FunctionCommonData(int k) {
            this->k = k;
            npt = k;
            for (int i = 0; i < 4; ++i)
                s[i] = Slice(i * k, (i + 1) * k - 1);
            s0 = std::vector<Slice>(NDIM);
            sh = std::vector<Slice>(NDIM);
            vk = std::vector<long>(NDIM);
            vq = std::vector<long>(NDIM);
            v2k = std::vector<long>(NDIM);
            for (std::size_t i = 0; i < NDIM; ++i) {
                s0[i] = s[0];
                sh[i] = Slice(0, (k - 1) / 2);
                vk[i] = k;
                vq[i] = npt;
                v2k[i] = 2 * k;
            }
            key0 = Key<NDIM> (0, Vector<Translation, NDIM> (0));

            _init_twoscale();
            _init_quadrature(k, npt, quad_x, quad_w, quad_phi, quad_phiw,
                             quad_phit);
        }

    public:
        typedef Tensor<T> tensorT; ///< Type of tensor used to hold coeff

        int k; ///< order of the wavelet
        int npt; ///< no. of quadrature points
        Slice s[4]; ///< s[0]=Slice(0,k-1), s[1]=Slice(k,2*k-1), etc.
        std::vector<Slice> s0; ///< s[0] in each dimension to get scaling coeff
        std::vector<Slice> sh; ///< Slice(0,(k-1)/2) in each dimension for autorefine test
        std::vector<long> vk; ///< (k,...) used to initialize Tensors
        std::vector<long> v2k; ///< (2k,...) used to initialize Tensors
        std::vector<long> vq; ///< (npt,...) used to initialize Tensors

        Key<NDIM> key0; ///< Key for root node

        Tensor<double> quad_x; ///< quadrature points
        Tensor<double> quad_w; ///< quadrature weights
        Tensor<double> quad_phi; ///< quad_phi(i,j) = at x[i] value of phi[j]
        Tensor<double> quad_phit; ///< transpose of quad_phi
        Tensor<double> quad_phiw; ///< quad_phiw(i,j) = at x[i] value of w[i]*phi[j]

        Tensor<double> h0, h1, g0, g1; ///< The separate blocks of twoscale coefficients
        Tensor<double> hg, hgT; ///< The full twoscale coeff (2k,2k) and transpose
        Tensor<double> hgsonly; ///< hg[0:k,:]

        static const FunctionCommonData<T, NDIM>&
        get(int k) {
            MADNESS_ASSERT(k > 0 && k <= MAXK);
            if (!data[k-1]) data[k-1] = new FunctionCommonData<T,NDIM>(k);
            return *(data[k-1]);
        }

        /// Initialize the quadrature information

        /// Made public with all arguments thru interface for reuse in FunctionImpl::err_box
        static void
        _init_quadrature(int k, int npt, Tensor<double>& quad_x, Tensor<
                         double>& quad_w, Tensor<double>& quad_phi,
                         Tensor<double>& quad_phiw, Tensor<double>& quad_phit);
    };



    /// FunctionNode holds the coefficients, etc., at each node of the 2^NDIM-tree
    template<typename T, std::size_t NDIM>
    class FunctionNode {
    public:
    	typedef GenTensor<T> coeffT;
    	typedef Tensor<T> tensorT;
    private:
        // Should compile OK with these volatile but there should
        // be no need to set as volatile since the container internally
        // stores the entire entry as volatile

        coeffT _coeffs; ///< The coefficients, if any
        double _norm_tree; ///< After norm_tree will contain norm of coefficients summed up tree
        bool _has_children; ///< True if there are children
        mutable double _normf; ///< Frobenius norm of this

    public:
        typedef WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > dcT; ///< Type of container holding the nodes
        /// Default constructor makes node without coeff or children
        FunctionNode() :
                _coeffs(), _norm_tree(1e300), _has_children(false) {
        }

        /// Constructor from given coefficients with optional children

        /// Note that only a shallow copy of the coeff are taken so
        /// you should pass in a deep copy if you want the node to
        /// take ownership.
        explicit
        FunctionNode(const coeffT& coeff, bool has_children = false) :
                _coeffs(coeff), _norm_tree(1e300), _has_children(has_children) {
        }

        explicit
        FunctionNode(const coeffT& coeff, double norm_tree, bool has_children) :
            _coeffs(coeff), _norm_tree(norm_tree), _has_children(has_children) {
        }

        FunctionNode(const FunctionNode<T, NDIM>& other) {
            *this = other;
        }

        FunctionNode<T, NDIM>&
        operator=(const FunctionNode<T, NDIM>& other) {
            if (this != &other) {
                coeff() = copy(other.coeff());
                _norm_tree = other._norm_tree;
                _has_children = other._has_children;
            }
            return *this;
        }

        /// Copy with possible type conversion of coefficients, copying all other state

        /// Choose to not overload copy and type conversion operators
        /// so there are no automatic type conversions.
        template<typename Q>
        FunctionNode<Q, NDIM>
        convert() const {
            return FunctionNode<Q, NDIM> (copy(coeff()), _has_children);
        }

        /// Returns true if there are coefficients in this node
        bool
        has_coeff() const {
        	return _coeffs.has_data();
        }

        bool exists() const {return this->has_data();}

        /// Returns true if this node has children
        bool
        has_children() const {
            return _has_children;
        }

        /// Returns true if this does not have children
        bool
        is_leaf() const {
            return !_has_children;
        }

        /// Returns true if this node is invalid (no coeffs and no children)
        bool
        is_invalid() const {
            return !(has_coeff() || has_children());
        }

        /// Returns a non-const reference to the tensor containing the coeffs

        /// Returns an empty tensor if there are no coefficients.
        coeffT&
        coeff() {
            MADNESS_ASSERT(_coeffs.ndim() == -1 || (_coeffs.dim(0) <= 2
                                                    * MAXK && _coeffs.dim(0) >= 0));
            return const_cast<coeffT&>(_coeffs);
        }

        /// Returns a const reference to the tensor containing the coeffs

        /// Returns an empty tensor if there are no coefficeints.
        const coeffT&
        coeff() const {
            return const_cast<const coeffT&>(_coeffs);
        }


        /// Returns the number of coefficients in this node
        size_t size() const {
        	return _coeffs.size();
        }


        /// Returns an empty tensor if there are no coefficients.
        Tensor<T>
        full_tensor_copy() const {
#if HAVE_GENTENSOR
        	Tensor<T> result;
        	if (this->has_coeff()) {
				result=_coeffs.reconstruct_tensor();
        	} else {
        		result=Tensor<T>();
        	}
#else
        	Tensor<T> result=copy(coeff());
#endif

        	return result;
        }

        /// Returns an empty tensor if there are no coefficients.
        Tensor<T>&
        full_tensor_reference() {
#if HAVE_GENTENSOR
        	MADNESS_ASSERT(_coeffs.tensor_type()==TT_FULL);
        	return _coeffs.full_tensor();
#else
        	return coeff();
#endif
        }

        /// Returns an empty tensor if there are no coefficients.
        const Tensor<T>&
        full_tensor_reference() const {
#if HAVE_GENTENSOR
        	MADNESS_ASSERT(_coeffs.type()==TT_FULL);
        	return _coeffs.full_tensor();
#else
        	return coeff();
#endif
        }

    public:


        /// Transforms the coeffs to the SR
        FunctionNode<T, NDIM>&
        node_to_low_rank(const double& thresh=FunctionDefaults<NDIM>::get_thresh()) {
        	MADNESS_ASSERT(thresh>0.0);
#if HAVE_GENTENSOR
        	to_low_rank(_coeffs,thresh,FunctionDefaults<NDIM>::get_tensor_type());
#endif
        	return *this;
        }

        /// Transforms the coeffs to the SR
        FunctionNode<T, NDIM>&
        node_to_full_rank() {
#if HAVE_GENTENSOR
        	to_full_rank(_coeffs);
#endif
        	return *this;
        }

        /// reduces the rank of the coefficients (if applicable)
        void reduceRank(const double& eps) {
        	_coeffs.reduceRank(eps);
        }


        /// Transforms the coeffs to the SR
        FunctionNode<T, NDIM>&
        ftr2sr(const double& thresh) {
        	MADNESS_ASSERT(thresh>0.0);
#if HAVE_GENTENSOR
        	to_low_rank(_coeffs,thresh,FunctionDefaults<NDIM>::get_tensor_type());
#endif
        	return *this;
        }

        /// Sets \c has_children attribute to value of \c flag.
        Void
        set_has_children(bool flag) {
            _has_children = flag;
            return None;
        }

        /// Sets \c has_children attribute to true recurring up to ensure connected
        Void
        set_has_children_recursive(const typename FunctionNode<T,NDIM>::dcT& c,const Key<NDIM>& key) {
            //madness::print("   set_chi_recu: ", key, *this);
            PROFILE_MEMBER_FUNC(FunctionNode);
            if (!(has_children() || has_coeff() || key.level()==0)) {
                // If node already knows it has children or it has
                // coefficients then it must already be connected to
                // its parent.  If not, the node was probably just
                // created for this operation and must be connected to
                // its parent.
                Key<NDIM> parent = key.parent();
                const_cast<dcT&>(c).task(parent, &FunctionNode<T,NDIM>::set_has_children_recursive, c, parent, TaskAttributes::hipri());
                //madness::print("   set_chi_recu: forwarding",key,parent);
            }
            _has_children = true;
            return None;
        }

        /// Sets \c has_children attribute to value of \c !flag
        void set_is_leaf(bool flag) {
            _has_children = !flag;
        }

        /// Takes a \em shallow copy of the coeff --- same as \c this->coeff()=coeff
        void set_coeff(const coeffT& coeffs) {
            coeff() = coeffs;
            if ((_coeffs.dim(0) < 0) || (_coeffs.dim(0)>2*MAXK)) {
                print("set_coeff: may have a problem");
                print("set_coeff: coeff.dim[0] =", coeffs.dim(0), ", 2* MAXK =", 2*MAXK);
            }
            MADNESS_ASSERT(coeffs.dim(0)<=2*MAXK && coeffs.dim(0)>=0);
        }

        /// Clears the coefficients (has_coeff() will subsequently return false)
        void clear_coeff() {
#if HAVE_GENTENSOR
//        	const TensorType tt=_coeffs.tensor_type();
            coeff()=coeffT();
#else
            coeff() = Tensor<T>();
#endif
        }

        /// Scale the coefficients of this node
        template <typename Q>
        void scale(Q a) {
        	_coeffs.scale(a);
        }

        /// Sets the value of norm_tree
        Void set_norm_tree(double norm_tree) {
            _norm_tree = norm_tree;
            return None;
        }

        /// Gets the value of norm_tree
        double get_norm_tree() const {
            return _norm_tree;
        }


        /// General bi-linear operation --- this = this*alpha + other*beta

        /// This/other may not have coefficients.  Has_children will be
        /// true in the result if either this/other have children.
        template <typename Q, typename R>
        Void gaxpy_inplace(const T& alpha, const FunctionNode<Q,NDIM>& other, const R& beta) {
            PROFILE_MEMBER_FUNC(FuncNode);
            if (other.has_children())
                _has_children = true;
            if (has_coeff()) {
                if (other.has_coeff()) {
                    coeff().gaxpy(alpha,other.coeff(),beta);
                }
                else {
                    coeff().scale(alpha);
                }
            }
            else if (other.has_coeff()) {
                coeff() = other.coeff()*beta; //? Is this the correct type conversion?
            }
            return None;
        }

        /// Accumulate inplace and if necessary connect node to parent
        Void accumulate2(const tensorT& t, const typename FunctionNode<T,NDIM>::dcT& c,
        		const Key<NDIM>& key) {

            if (has_coeff()) {
            	MADNESS_ASSERT(coeff().tensor_type()==TT_FULL);
//            	if (coeff().type==TT_FULL) {
        		coeff() += coeffT(t,-1.0,TT_FULL);
//            	} else {
//            		tensorT cc=coeff().full_tensor_copy();;
//            		cc += t;
//            		coeff()=coeffT(cc,args);
//            	}
            }
            else {
                // No coeff and no children means the node is newly
                // created for this operation and therefore we must
                // tell its parent that it exists.
            	coeff() = coeffT(t,-1.0,TT_FULL);
//                coeff() = copy(t);
//                coeff() = coeffT(t,args);
                if ((!_has_children) && key.level()> 0) {
                    Key<NDIM> parent = key.parent();
                    const_cast<dcT&>(c).task(parent, &FunctionNode<T,NDIM>::set_has_children_recursive, c, parent);
                }
            }
            return None;
        }


        /// Accumulate inplace and if necessary connect node to parent
        Void accumulate(const coeffT& t, const typename FunctionNode<T,NDIM>::dcT& c,
        		const Key<NDIM>& key, const TensorArgs& args) {

            if (has_coeff()) {
                coeff().add_SVD(t,args.thresh);
//           		coeff() += t;
//           		t.accumulate_into(coeff(),args.thresh,1.0);
            } else {
                // No coeff and no children means the node is newly
                // created for this operation and therefore we must
                // tell its parent that it exists.
            	coeff() = copy(t);
                if ((!_has_children) && key.level()> 0) {
                    Key<NDIM> parent = key.parent();
                    const_cast<dcT&>(c).task(parent, &FunctionNode<T,NDIM>::set_has_children_recursive, c, parent);
                }
            }
            return None;
        }

        T trace_conj(const FunctionNode<T,NDIM>& rhs) const {
        	return this->_coeffs.trace_conj((rhs._coeffs));
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar & coeff() & _has_children & double(_norm_tree);
        }

    };

    template <typename T, std::size_t NDIM>
    std::ostream& operator<<(std::ostream& s, const FunctionNode<T,NDIM>& node) {
        s << "(has_coeff=" << node.has_coeff() << ", has_children=" << node.has_children() << ", norm=";
        double norm = node.has_coeff() ? node.coeff().normf() : 0.0;
        if (norm < 1e-12)
            norm = 0.0;
//        norm=node.get_norm_tree();
        s << norm << ")";
        return s;
    }


    /// FunctionImpl holds all Function state to facilitate shallow copy semantics

    /// Since Function assignment and copy constructors are shallow it
    /// greatly simplifies maintaining consistent state to have all
    /// (permanent) state encapsulated in a single class.  The state
    /// is shared between instances using a shared_ptr<FunctionImpl>.
    ///
    /// The FunctionImpl inherits all of the functionality of WorldContainer
    /// (to store the coefficients) and WorldObject<WorldContainer> (used
    /// for RMI and for its unqiue id).
    ///
    /// The class methods are public to avoid painful multiple friend template
    /// declarations for Function and FunctionImpl ... but this trust should not be
    /// abused ... NOTHING except FunctionImpl methods should mess with FunctionImplData.
    /// The LB stuff might have to be an exception.
    template <typename T, std::size_t NDIM>
    class FunctionImpl : public WorldObject< FunctionImpl<T,NDIM> > {
    private:
        typedef WorldObject< FunctionImpl<T,NDIM> > woT; ///< Base class world object type

    public:
        typedef FunctionImpl<T,NDIM> implT; ///< Type of this class (implementation)
        typedef std::shared_ptr< FunctionImpl<T,NDIM> > pimplT; ///< pointer to this class
        typedef Tensor<T> tensorT; ///< Type of tensor for anything but to hold coeffs
        typedef Vector<Translation,NDIM> tranT; ///< Type of array holding translation
        typedef Key<NDIM> keyT; ///< Type of key
        typedef FunctionNode<T,NDIM> nodeT; ///< Type of node
        typedef GenTensor<T> coeffT; ///< Type of tensor used to hold coeffs
        typedef WorldContainer<keyT,nodeT> dcT; ///< Type of container holding the coefficients
        typedef std::pair<const keyT,nodeT> datumT; ///< Type of entry in container
        typedef Vector<double,NDIM> coordT; ///< Type of vector holding coordinates

        //template <typename Q, int D> friend class Function;
        template <typename Q, std::size_t D> friend class FunctionImpl;

        friend class LoadBalImpl<NDIM>;
        friend class LBTree<NDIM>;

        World& world;

    private:
        int k; ///< Wavelet order
        double thresh; ///< Screening threshold
        int initial_level; ///< Initial level for refinement
        int max_refine_level; ///< Do not refine below this level
        int truncate_mode; ///< 0=default=(|d|<thresh), 1=(|d|<thresh/2^n), 1=(|d|<thresh/4^n);
        bool autorefine; ///< If true, autorefine where appropriate
        bool truncate_on_project; ///< If true projection inserts at level n-1 not n
        bool nonstandard; ///< If true, compress keeps scaling coeff
        TensorArgs targs; ///< type of tensor to be used in the FunctionNodes

        const FunctionCommonData<T,NDIM>& cdata;

        std::shared_ptr< FunctionFunctorInterface<T,NDIM> > functor;

        bool on_demand; ///< does this function have an additional functor?
        bool compressed; ///< Compression status
        bool redundant; ///< If true, function keeps sum coefficients on all levels

        dcT coeffs; ///< The coefficients

        // Disable the default copy constructor
        FunctionImpl(const FunctionImpl<T,NDIM>& p);

    public:

        /// Initialize function impl from data in factory
        FunctionImpl(const FunctionFactory<T,NDIM>& factory)
                : WorldObject<implT>(factory._world)
                , world(factory._world)
                , k(factory._k)
                , thresh(factory._thresh)
                , initial_level(factory._initial_level)
                , max_refine_level(factory._max_refine_level)
                , truncate_mode(factory._truncate_mode)
                , autorefine(factory._autorefine)
                , truncate_on_project(factory._truncate_on_project)
                , nonstandard(false)
                , targs(factory._thresh,FunctionDefaults<NDIM>::get_tensor_type())
                , cdata(FunctionCommonData<T,NDIM>::get(k))
                , functor(factory.get_functor())
                , on_demand(factory._is_on_demand)
                , compressed(false)
                , redundant(false)
                , coeffs(world,factory._pmap,false)
                //, bc(factory._bc)
            {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            // !!! Ensure that all local state is correctly formed
            // before invoking process_pending for the coeffs and
            // for this.  Otherwise, there is a race condition.
            MADNESS_ASSERT(k>0 && k<=MAXK);

            bool empty = (factory._empty or is_on_demand());
            bool do_refine = factory._refine;

            if (do_refine)
                initial_level = std::max(0,initial_level - 1);

            if (empty) { // Do not set any coefficients at all
            // additional functors are only evaluated on-demand
            } else if (functor) { // Project function and optionally refine
                insert_zero_down_to_initial_level(cdata.key0);
                typename dcT::const_iterator end = coeffs.end();
                for (typename dcT::const_iterator it=coeffs.begin(); it!=end; ++it) {
                    if (it->second.is_leaf())
                        woT::task(coeffs.owner(it->first), &implT::project_refine_op, it->first, do_refine,
                             functor->special_points());
                }
            }
            else { // Set as if a zero function
                initial_level = 1;
                insert_zero_down_to_initial_level(keyT(0));
            }

            coeffs.process_pending();
            this->process_pending();
            if (factory._fence && functor)
                world.gop.fence();

        }

        /// Copy constructor

        /// Allocates a \em new function in preparation for a deep copy
        ///
        /// By default takes pmap from other but can also specify a different pmap.
        /// Does \em not copy the coefficients ... creates an empty container.
        template <typename Q>
        FunctionImpl(const FunctionImpl<Q,NDIM>& other,
                     const std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >& pmap,
                     bool dozero)
                : WorldObject<implT>(other.world)
                , world(other.world)
                , k(other.k)
                , thresh(other.thresh)
                , initial_level(other.initial_level)
                , max_refine_level(other.max_refine_level)
                , truncate_mode(other.truncate_mode)
                , autorefine(other.autorefine)
                , truncate_on_project(other.truncate_on_project)
                , nonstandard(other.nonstandard)
                , targs(other.targs)
                , cdata(FunctionCommonData<T,NDIM>::get(k))
                , functor()
                , on_demand(false)	// since functor() is an default ctor
                , compressed(other.compressed)
                , redundant(other.redundant)
                , coeffs(world, pmap ? pmap : other.coeffs.get_pmap())
                //, bc(other.bc)
        {
            if (dozero) {
                initial_level = 1;
                insert_zero_down_to_initial_level(cdata.key0);
            }
            coeffs.process_pending();
            this->process_pending();
        }

        virtual ~FunctionImpl() { }

        const std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >& get_pmap() const {
            return coeffs.get_pmap();
        }

        /// Copy coeffs from other into self
        template <typename Q>
        void copy_coeffs(const FunctionImpl<Q,NDIM>& other, bool fence) {
            typename FunctionImpl<Q,NDIM>::dcT::const_iterator end = other.coeffs.end();
            for (typename FunctionImpl<Q,NDIM>::dcT::const_iterator it=other.coeffs.begin();
                    it!=end; ++it) {
                const keyT& key = it->first;
                const typename FunctionImpl<Q,NDIM>::nodeT& node = it->second;
                coeffs.replace(key,node. template convert<Q>());
            }
            if (fence)
                world.gop.fence();
        }


        template <typename Q, typename R>
        struct do_gaxpy_inplace {
            typedef Range<typename FunctionImpl<Q,NDIM>::dcT::const_iterator> rangeT;
            FunctionImpl<T,NDIM>* f;
            T alpha;
            R beta;
            do_gaxpy_inplace() {};
            do_gaxpy_inplace(FunctionImpl<T,NDIM>* f, T alpha, R beta) : f(f), alpha(alpha), beta(beta) {}
            bool operator()(typename rangeT::iterator& it) const {
                const keyT& key = it->first;
                const FunctionNode<Q,NDIM>& other_node = it->second;
                // Use send to get write accessor and automated construction if missing
                f->coeffs.send(key, &nodeT:: template gaxpy_inplace<Q,R>, alpha, other_node, beta);
                return true;
            }
            template <typename Archive>
            void serialize(Archive& ar) {}
        };

        /// Inplace general bilinear operation
        template <typename Q, typename R>
        void gaxpy_inplace(const T& alpha,const FunctionImpl<Q,NDIM>& other, const R& beta, bool fence) {
            MADNESS_ASSERT(get_pmap() == other.get_pmap());
            if (alpha != T(1.0)) scale_inplace(alpha,false);
            typedef Range<typename FunctionImpl<Q,NDIM>::dcT::const_iterator> rangeT;
            typedef do_gaxpy_inplace<Q,R> opT;
            world.taskq.for_each<rangeT,opT>(rangeT(other.coeffs.begin(), other.coeffs.end()), opT(this, T(1.0), beta));
            if (fence)
                world.gop.fence();
        }

        template <typename Archive>
        void load(Archive& ar) {
            // WE RELY ON K BEING STORED FIRST
            int kk;
            ar & kk;

            MADNESS_ASSERT(kk==k);

            // note that functor should not be (re)stored
            ar & thresh & initial_level & max_refine_level & truncate_mode
            & autorefine & truncate_on_project & nonstandard & compressed ; //& bc;

            ar & coeffs;
            world.gop.fence();
        }

        template <typename Archive>
        void store(Archive& ar) {
            // WE RELY ON K BEING STORED FIRST

            // note that functor should not be (re)stored
            ar & k & thresh & initial_level & max_refine_level & truncate_mode
            & autorefine & truncate_on_project & nonstandard & compressed ; //& bc;

            ar & coeffs;
            world.gop.fence();
        }

        /// Returns true if the function is compressed.
        bool is_compressed() const {
            return compressed;
        }

        /// Returns true if the function is redundant.
        bool is_redundant() const {
            return redundant;
        }

        bool is_nonstandard() const {return nonstandard;}

        void set_functor(const std::shared_ptr<FunctionFunctorInterface<T,NDIM> > functor1) {
        	this->on_demand=true;
        	functor=functor1;
        }

        std::shared_ptr<FunctionFunctorInterface<T,NDIM> > get_functor() {
        	MADNESS_ASSERT(this->functor);
        	return functor;
        }

        void unset_functor() {
        	this->on_demand=false;
        	functor.reset();
        }

        /// if this function is an on-demand function and its functor provides
        /// a muster, return that muster
        std::shared_ptr<FunctionImpl<T,NDIM> > get_muster() {

        	MADNESS_ASSERT(is_on_demand());
        	MADNESS_ASSERT(functor);
            CompositeFunctorInterface<T,NDIM,3>* func=
            		dynamic_cast<CompositeFunctorInterface<T,NDIM,3>* >(&(*functor));
            MADNESS_ASSERT(func);
            return func->get_muster();

        }

        bool is_on_demand() const {return on_demand;};

        TensorType get_tensor_type() const {return targs.tt;}
        TensorArgs get_tensor_args() const {return targs;}

        double get_thresh() const {return thresh;}

        void set_thresh(double value) {thresh = value;}

        bool get_autorefine() const {return autorefine;}

        void set_autorefine(bool value) {autorefine = value;}

        int get_k() const {return k;}

        const dcT& get_coeffs() const {return coeffs;}

        dcT& get_coeffs() {return coeffs;}

        /// Adds a constant to the function.  Local operation, optional fence

        /// In scaling function basis must add value to first polyn in
        /// each box with appropriate scaling for level.  In wavelet basis
        /// need only add at level zero.
        void add_scalar_inplace(T t, bool fence);

        /// Initialize nodes to zero function at initial_level of refinement.

        /// Works for either basis.  No communication.
        void insert_zero_down_to_initial_level(const keyT& key);

        /// Truncate according to the threshold with optional global fence

        /// If thresh<=0 the default value of this->thresh is used
        void truncate(double tol, bool fence) {
            // Cannot put tol into object since it would make a race condition
            if (tol <= 0.0)
                tol = thresh;
            if (world.rank() == coeffs.owner(cdata.key0))
                truncate_spawn(cdata.key0,tol);
            if (fence)
                world.gop.fence();
        }

        /// Returns true if after truncation this node has coefficients

        /// Assumed to be invoked on process owning key.  Possible non-blocking
        /// communication.
        Future<bool> truncate_spawn(const keyT& key, double tol);

        /// Actually do the truncate operation
        bool truncate_op(const keyT& key, double tol, const std::vector< Future<bool> >& v);

        /// Evaluate function at quadrature points in the specified box
        void fcube(const keyT& key, const FunctionFunctorInterface<T,NDIM>& f, const Tensor<double>& qx, tensorT& fval) const;
        void fcube(const keyT& key,  T (*f)(const coordT&), const Tensor<double>& qx, tensorT& fval) const;
        //~ template<typename FF>
        //~ void fcube(const keyT& key, const FF& f, const Tensor<double>& qx, tensorT& fval) const;

        const keyT& key0() const {
            return cdata.key0;
        }

        void print_tree(Level maxlevel = 10000) const;

        void do_print_tree(const keyT& key, Level maxlevel) const;


        /// Print a plane ("xy", "xz", or "yz") containing the point x to file

        /// works for all dimensions; we walk through the tree, and if a leaf node
        /// inside the sub-cell touches the plane we print it in pstricks format
        void print_plane(const std::string filename, const std::string plane, const coordT& x_user,
        		const Tensor<double>& cell_user) {

        	// translate verbose plane to something computer-readable
        	int dim0;
        	int dim1;
        	if (plane=="xy") {
        		dim0=0;
        		dim1=1;
        	} else if (plane=="xz") {
        		dim0=0;
        		dim1=2;
        	} else if (plane=="yz") {
        		dim0=1;
        		dim1=2;
        	} else {
        		print("unknown slice in WF::printSlice: ", plane);
        		MADNESS_ASSERT(0);
        	}

        	coordT x_sim;
        	user_to_sim<NDIM>(x_user,x_sim);
        	x_sim[2]+=1.e-10;

        	// prepare file
        	FILE * pFile;
        	pFile = fopen(filename.c_str(), "w");

        	fprintf(pFile,"\\psset{unit=10cm}\n");
        	fprintf(pFile,"\\begin{pspicture}(0,0)(1,1)\n");
        	fprintf(pFile,"\\pslinewidth=0.005pt\n");



        	typename dcT::iterator end = coeffs.end();
            for (typename dcT::iterator it=coeffs.begin(); it!=end; ++it) {
                const keyT& key = it->first;
                nodeT& node = it->second;

                // thisKeyContains ignores dim0 and dim1
                if (key.thisKeyContains(x_sim,dim0,dim1) and node.is_leaf() and (node.has_coeff())) {

                	// next: only key inside the given box
                	for (size_t i=0; i<NDIM; i++) {
                		if (x_user[i]<cell_user(i,0) or x_user[i]>cell_user(i,1)) continue;
                	}

                	Level n=key.level();
                	Vector<Translation,NDIM> l=key.translation();
                	// get the diametral edges of the node in the plotting plane
                	double scale=std::pow(0.5,double(n));
                	double xloleft = scale*l[dim0];
                	double yloleft = scale*l[dim1];
                	double xhiright = scale*(l[dim0]+1);
                	double yhiright = scale*(l[dim1]+1);

                	// color
                	const int rank=node.coeff().rank();
         	        const double maxrank=40.0;
           	        double hue=0.7-(0.7/maxrank)*(rank);
           	        double color= std::max(0.0,hue);

           	        // plot only if Box is within bounds
//                	if ((xloleft>=x0) and (yloleft>=y0)
//                			and (xhiright<=x1) and (yhiright<=y1)) {
           	        if (1) {

                		fprintf(pFile,"\\newhsbcolor{mycolor}{%8.4f 1.0 0.7}\n",color);
                		fprintf(pFile,"\\psframe["//linewidth=0.5pt,"
                				"fillstyle=solid,"
                				"fillcolor=mycolor]"
                				"(%12.8f,%12.8f)(%12.8f,%12.8f)\n",
                				xloleft,yloleft,xhiright,yhiright);
                	}
                }
            }
            fprintf(pFile,"\\end{pspicture}\n");
            fclose(pFile);
        }


        /// Compute by projection the scaling function coeffs in specified box
        tensorT project(const keyT& key) const;

        /// Returns the truncation threshold according to truncate_method
        double truncate_tol(double tol, const keyT& key) const {
            if (truncate_mode == 0) {
                return tol;
            }
            else if (truncate_mode == 1) {
                double L = FunctionDefaults<NDIM>::get_cell_min_width();
                return tol*std::min(1.0,pow(0.5,double(key.level()))*L);
            }
            else if (truncate_mode == 2) {
                double L = FunctionDefaults<NDIM>::get_cell_min_width();
                return tol*std::min(1.0,pow(0.25,double(key.level()))*L*L);
            }
            else {
                MADNESS_EXCEPTION("truncate_mode invalid",truncate_mode);
            }
        }

        /// Returns patch referring to coeffs of child in parent box
        std::vector<Slice> child_patch(const keyT& child) const {
            std::vector<Slice> s(NDIM);
            const Vector<Translation,NDIM>& l = child.translation();
            for (std::size_t i=0; i<NDIM; ++i)
                s[i] = cdata.s[l[i]&1]; // Lowest bit of translation
            return s;
        }

        /// Projection with optional refinement
        Void project_refine_op(const keyT& key, bool do_refine,
                const std::vector<Vector<double,NDIM> >& specialpts);

        /// Compute the Legendre scaling functions for multiplication

        /// Evaluate parent polyn at quadrature points of a child.  The prefactor of
        /// 2^n/2 is included.  The tensor must be preallocated as phi(k,npt).
        /// Refer to the implementation notes for more info.
        void phi_for_mul(Level np, Translation lp, Level nc, Translation lc, Tensor<double>& phi) const;

        /// Directly project parent coeffs to child coeffs

        /// Currently used by diff, but other uses can be anticipated
        const coeffT parent_to_child(const coeffT& s, const keyT& parent, const keyT& child) const;


    /// Returns the box at level n that contains the given point in simulation coordinates
        Key<NDIM> simpt2key(const coordT& pt, Level n) const {
            Vector<Translation,NDIM> l;
            double twon = std::pow(2.0, double(n));
            for (std::size_t i=0; i<NDIM; ++i) {
                l[i] = Translation(twon*pt[i]);
            }
            return Key<NDIM>(n,l);
        }

        /// Get the scaling function coeffs at level n starting from NS form
        // N=2^n, M=N/q, q must be power of 2
        // q=0 return coeffs [N,k] for direct sum
        // q>0 return coeffs [k,q,M] for fft sum
        tensorT coeffs_for_jun(Level n, long q=0) {
            MADNESS_ASSERT(compressed && nonstandard && NDIM<=3);
            tensorT r,r0;
            long N=1<<n;
            long M = (q ? N/q: N);
            if (q==0) {
                q = 1;
                long dim[2*NDIM];
                for (std::size_t d=0; d<NDIM; ++d) {
                    dim[d ] = N;
                    dim[d+NDIM] = cdata.k;
                }
                tensorT rr(2*NDIM,dim);
                r0=r=rr;
                //NNkk->MqMqkk, since fuse is not allowed. Now needs to move back to 2*NDIM, since tensor max dim is 6
                //for (int d=NDIM-1; d>=0; --d) r.splitdim_inplace_base(d,M,q);
            } else {
                long dim[2*NDIM];
                for (std::size_t d=0; d<NDIM; ++d) {
                    //dim[d+NDIM*2] = M;
                    dim[d+NDIM ] = N;
                    dim[d ] = cdata.k;
                }
                tensorT rr(2*NDIM,dim);
                r0=rr;
                /*vector<long> map(3*NDIM);
                 for (int d=0; d<NDIM; ++d) {
                 map[d]=d+2*NDIM;
                 map[NDIM+d]=2*d+1;
                 map[2*NDIM+d]=2*d;
                 }
                 r.mapdim_inplace_base(map);
                 //print(rr);
                 //for (int d=1; d<NDIM; ++d) rr.swapdim_inplace_base(2*NDIM+d,NDIM+d); //kkqqMM->kkqMqM
                 //print(rr);
                 //for (int d=0; d<NDIM; ++d) rr.swapdim_inplace_base(NDIM+2*d,NDIM+2*d-1); //kkqMqM->kkMqMq
                 //print(rr);
                 //for (int d=0; d<NDIM; ++d) rr.fusedim_inplace_base(NDIM+d); //->kkNN
                 //seems that this fuse is not allowed :(

                 //print(rr);
                 */
                r=rr.cycledim(NDIM,0,-1); //->NNkk or MqMqkk
            }
            print("faking done M q r(fake) r0(real)",M,q,"\n", std::vector<long> (r.dims(),r.dims()+6),std::vector<long> (r0.dims(),r0.dims()+6));
            ProcessID me = world.rank();
            Vector<long,NDIM> t(N);

            Vector<long,NDIM> powq, powN, powM;
            long NDIM1 = NDIM-1;
            powM[NDIM1]=powq[NDIM1]=powN[NDIM1]=1;
            for (int d=NDIM1-1; d>=0; --d) {
                powM[d] = powM[d+1]*M;
                powq[d] = powq[d+1]*q;
                powN[d] = powN[d+1]*N;
            }
            long powMNDIM = powM[0]*M;

            for (IndexIterator it(t); it; ++it) {
                keyT key(n, Vector<Translation,NDIM>(*it));
                if (coeffs.owner(key) == me) {
                    typename dcT::iterator it = coeffs.find(key).get();
                    coeffT qq;

                    if (it == coeffs.end()) {
                        // must get from above
                        typedef std::pair< keyT,coeffT > pairT;
                        Future<pairT> result;
                        sock_it_to_me(key, result.remote_ref(world));
                        const keyT& parent = result.get().first;
//                        const tensorT& t = result.get().second.full_tensor_copy();
                        const coeffT& t = result.get().second;

                        qq = (parent_to_child(t, parent, key));
                    } else {
                        qq = copy(it->second.coeff());
                    }
                    std::vector<Slice> s(NDIM*2);
                    long ll = 0;
                    for (std::size_t d=0; d<NDIM; ++d) {
                        Translation l = key.translation()[d];
                        long dum = long(float(l)/q);
                        ll += (l - dum*q)*powMNDIM*powq[d] + dum*powM[d];
                        //ll += (l % q)*powM[NDIM]*pow((double)q,NDIM-d-1) + (l/q)*pow((double)M,NDIM-d-1);

                        //print("translation",l);
                        //s[d       ] = Slice(l,l,0);
                        //s[d+NDIM  ] = Slice(l%q,l%q,0);
                        //s[d+NDIM] = Slice(0,k-1,1);
                    }
                    //long dum = ll;
                    for (std::size_t d=0; d<NDIM; ++d) {
                        Translation l = Translation(float(ll) / powN[d]);
                        //Translation l = ll / pow((double)N,NDIM-d-1);
                        s[d ] = Slice(l,l,0);
                        s[d+NDIM] = Slice(0,k-1,1);
                        ll = ll - l*powN[d];
                        //ll = ll % long(pow((double)N,NDIM-d-1));
                    }
                    //print(s, dum, key.translation());
                    coeffT qqq=qq(cdata.s0);
                    r(s) = qqq.full_tensor_copy();

                }
            }

            world.gop.fence();
            world.gop.sum(r0);
            //print(r,r0);

            return r0;
        }

        template <typename Q>
        __Tensor<Q> coeffs2values(const keyT& key, const __Tensor<Q>& coeff) const {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            double scale = pow(2.0,0.5*NDIM*key.level())/sqrt(FunctionDefaults<NDIM>::get_cell_volume());
            return transform(coeff,cdata.quad_phit).scale(scale);
        }

        template <typename Q>
//        Tensor<Q> values2coeffs(const keyT& key, const Tensor<Q>& values) const {
        coeffT values2coeffs(const keyT& key, const coeffT& values) const {
//        __Tensor<Q> values2coeffs(const keyT& key, const __Tensor<Q>& values) const {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            double scale = pow(0.5,0.5*NDIM*key.level())*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
            return transform(values,cdata.quad_phiw).scale(scale);
        }

        template <typename Q>
        Tensor<Q> values2coeffs(const keyT& key, const Tensor<Q>& values) const {
//        coeffT values2coeffs(const keyT& key, const coeffT& values) const {
//        __Tensor<Q> values2coeffs(const keyT& key, const __Tensor<Q>& values) const {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            double scale = pow(0.5,0.5*NDIM*key.level())*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
            return transform(values,cdata.quad_phiw).scale(scale);
        }


        /// Compute the function values for multiplication

        /// Given coefficients from a parent cell, compute the value of
        /// the functions at the quadrature points of a child
        template <typename Q>
        __Tensor<Q> fcube_for_mul(const keyT& child, const keyT& parent, const __Tensor<Q>& coeff) const {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            if (child.level() == parent.level()) {
                return coeffs2values(parent, coeff);
            }
            else if (child.level() < parent.level()) {
                MADNESS_EXCEPTION("FunctionImpl: fcube_for_mul: child-parent relationship bad?",0);
            }
            else {
                Tensor<double> phi[NDIM];
                for (size_t d=0; d<NDIM; d++) {
                    phi[d] = Tensor<double>(cdata.k,cdata.npt);
                    phi_for_mul(parent.level(),parent.translation()[d],
                                child.level(), child.translation()[d], phi[d]);
                }
                return general_transform(coeff,phi).scale(1.0/sqrt(FunctionDefaults<NDIM>::get_cell_volume()));;
            }
        }

        /// Compute the function values for multiplication

        /// Given coefficients from a parent cell, compute the value of
        /// the functions at the quadrature points of a child
        coeffT fcube_for_mul2(const keyT& child, const std::pair<keyT,coeffT>& arg) const {
            PROFILE_MEMBER_FUNC(FunctionImpl);

            const keyT& parent=arg.first;
            const coeffT& coeff=arg.second;

            if (child.level() == parent.level()) {
                return coeffs2values(parent, coeff);
            }
            else if (child.level() < parent.level()) {
                MADNESS_EXCEPTION("FunctionImpl: fcube_for_mul: child-parent relationship bad?",0);
            }
            else {
                Tensor<double> phi[NDIM];
                for (size_t d=0; d<NDIM; d++) {
                    phi[d] = Tensor<double>(cdata.k,cdata.npt);
                    phi_for_mul(parent.level(),parent.translation()[d],
                                child.level(), child.translation()[d], phi[d]);
                }
                return general_transform(coeff,phi).scale(1.0/sqrt(FunctionDefaults<NDIM>::get_cell_volume()));
            }
        }


        /// Compute the function values for multiplication

        /// same as fcube_for_mul, but first find the appropriate parent node
        ///
        /// Given coefficients from a parent cell, compute the value of
        /// the functions at the quadrature points of a child
        coeffT fcube_for_mul_too(const keyT& key) const {

            typedef std::pair<keyT,coeffT> argT;

			Future<argT> result;
			woT::task(coeffs.owner(key), &implT::sock_it_to_me, key, result.remote_ref(world));
			return woT::task(coeffs.owner(key), &implT::fcube_for_mul2, key, result);

        }


#if HAVE_GENTENSOR


        template <typename Q>
        Tensor<Q> coeffs2values(const keyT& key, const Tensor<Q>& coeff) const {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            double scale = pow(2.0,0.5*NDIM*key.level())/sqrt(FunctionDefaults<NDIM>::get_cell_volume());
            return transform(coeff,cdata.quad_phit).scale(scale);
        }

        /// Compute the function values for multiplication

        /// Given coefficients from a parent cell, compute the value of
        /// the functions at the quadrature points of a child
        template <typename Q>
        Tensor<Q> fcube_for_mul(const keyT& child, const keyT& parent, const Tensor<Q>& coeff) const {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            if (child.level() == parent.level()) {
                return coeffs2values(parent, coeff);
            }
            else if (child.level() < parent.level()) {
                MADNESS_EXCEPTION("FunctionImpl: fcube_for_mul: child-parent relationship bad?",0);
            }
            else {
                Tensor<double> phi[NDIM];
                for (std::size_t d=0; d<NDIM; ++d) {
                    phi[d] = Tensor<double>(cdata.k,cdata.npt);
                    phi_for_mul(parent.level(),parent.translation()[d],
                                child.level(), child.translation()[d], phi[d]);
                }
                return general_transform(coeff,phi).scale(1.0/sqrt(FunctionDefaults<NDIM>::get_cell_volume()));;
            }
        }
#endif

        /// Invoked as a task by mul with the actual coefficients
        template <typename L, typename R>
        Void do_mul(const keyT& key, const Tensor<L>& left, const std::pair< keyT, Tensor<R> >& arg) {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            const keyT& rkey = arg.first;
            const Tensor<R>& rcoeff = arg.second;
            //madness::print("do_mul: r", rkey, rcoeff.size());
            Tensor<R> rcube = fcube_for_mul(key, rkey, rcoeff);
            //madness::print("do_mul: l", key, left.size());
            Tensor<L> lcube = fcube_for_mul(key, key, left);

            Tensor<T> tcube(cdata.vk,false);
            TERNARY_OPTIMIZED_ITERATOR(T, tcube, L, lcube, R, rcube, *_p0 = *_p1 * *_p2;);
            double scale = pow(0.5,0.5*NDIM*key.level())*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
            tcube = transform(tcube,cdata.quad_phiw).scale(scale);
            coeffs.replace(key, nodeT(coeffT(tcube,targs),false));
            return None;
        }

        /// Invoked as a task by do_binary_op with the actual coefficients
        template <typename L, typename R, typename opT>
        Void do_binary_op(const keyT& key, const Tensor<L>& left,
                          const std::pair< keyT, Tensor<R> >& arg,
                          const opT& op) {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            const keyT& rkey = arg.first;
            const Tensor<R>& rcoeff = arg.second;
            Tensor<R> rcube = fcube_for_mul(key, rkey, rcoeff);
            Tensor<L> lcube = fcube_for_mul(key, key, left);

            Tensor<T> tcube(cdata.vk,false);
            op(key, tcube, lcube, rcube);
            double scale = pow(0.5,0.5*NDIM*key.level())*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
            tcube = transform(tcube,cdata.quad_phiw).scale(scale);
            coeffs.replace(key, nodeT(coeffT(tcube,targs),false));
            return None;
        }

        /// Invoked by result to perform result += alpha*left+beta*right in wavelet basis

        /// Does not assume that any of result, left, right have the same distribution.
        /// For most purposes result will start as an empty so actually are implementing
        /// out of place gaxpy.  If all functions have the same distribution there is
        /// no communication except for the optional fence.
        template <typename L, typename R>
        void gaxpy(T alpha, const FunctionImpl<L,NDIM>& left,
                   T beta, const FunctionImpl<R,NDIM>& right, bool fence) {
            // Loop over local nodes in both functions.  Add in left and subtract right.
            // Not that efficient in terms of memory bandwidth but ensures we do
            // not miss any nodes.
            typename FunctionImpl<L,NDIM>::dcT::const_iterator left_end = left.coeffs.end();
            for (typename FunctionImpl<L,NDIM>::dcT::const_iterator it=left.coeffs.begin();
                    it!=left_end;
                    ++it) {
                const keyT& key = it->first;
                const typename FunctionImpl<L,NDIM>::nodeT& other_node = it->second;
                coeffs.send(key, &nodeT:: template gaxpy_inplace<T,L>, 1.0, other_node, alpha);
            }
            typename FunctionImpl<R,NDIM>::dcT::const_iterator right_end = right.coeffs.end();
            for (typename FunctionImpl<R,NDIM>::dcT::const_iterator it=right.coeffs.begin();
                    it!=right_end;
                    ++it) {
                const keyT& key = it->first;
                const typename FunctionImpl<L,NDIM>::nodeT& other_node = it->second;
                coeffs.send(key, &nodeT:: template gaxpy_inplace<T,R>, 1.0, other_node, beta);
            }
            if (fence)
                world.gop.fence();
        }

        /// Unary operation applied inplace to the coefficients WITHOUT refinement, optional fence
        template <typename opT>
        void unary_op_coeff_inplace(const opT& op, bool fence) {
            typename dcT::iterator end = coeffs.end();
            for (typename dcT::iterator it=coeffs.begin(); it!=end; ++it) {
                const keyT& parent = it->first;
                nodeT& node = it->second;
                if (node.has_coeff()) {
//                    op(parent, node.coeff());
                	node.node_to_full_rank();
                    op(parent, node.full_tensor_reference());
                    node.node_to_low_rank();
//                	op(parent,node);
                }
            }
            if (fence)
                world.gop.fence();
        }

        /// Unary operation applied inplace to the coefficients WITHOUT refinement, optional fence
        template <typename opT>
        void unary_op_node_inplace(const opT& op, bool fence) {
            typename dcT::iterator end = coeffs.end();
            for (typename dcT::iterator it=coeffs.begin(); it!=end; ++it) {
                const keyT& parent = it->first;
                nodeT& node = it->second;
                op(parent, node);
            }
            if (fence)
                world.gop.fence();
        }

        /// Unary operation applied inplace to the coefficients WITHOUT refinement, optional fence
        template <typename opT>
        void flo_unary_op_node_inplace(const opT& op, bool fence) {

        	typedef Range<typename dcT::iterator> rangeT;
            typedef do_unary_op_value_inplace<opT> xopT;
            world.taskq.for_each<rangeT,opT>(rangeT(coeffs.begin(), coeffs.end()), op);
            if (fence)
                 world.gop.fence();
        }




        /// keep only the sum coefficients in each node
        struct do_keep_sum_coeffs {
             typedef Range<typename dcT::iterator> rangeT;
             implT* impl;

             /// constructor need impl for cdata
             do_keep_sum_coeffs(implT* impl) :impl(impl) {}

             bool operator()(typename rangeT::iterator& it) const {

            	 nodeT& node = it->second;
            	 coeffT s=copy(node.coeff()(impl->cdata.s0));
            	 node.coeff()=s;
                 return true;
             }
             template <typename Archive> void serialize(const Archive& ar) {}

        };


        /// reduce the rank of the nodes, optional fence
        struct do_reduce_rank {
            typedef Range<typename dcT::iterator> rangeT;

            // threshold for rank reduction / SVD truncation
            TensorArgs args;

            // constructor takes target precision
            do_reduce_rank() {}
            do_reduce_rank(const TensorArgs& targs) : args(targs) {}

            //
            bool operator()(typename rangeT::iterator& it) const {

                nodeT& node = it->second;
                node.reduceRank(args.thresh);
                return true;
            }
            template <typename Archive> void serialize(const Archive& ar) {}
        };



        /// change representation of nodes' coeffs to low rank, optional fence
        struct do_low_rank_inplace {
            typedef Range<typename dcT::iterator> rangeT;

            // threshold for rank reduction / SVD truncation
            double eps_;

            // constructor takes target precision
            do_low_rank_inplace() : eps_(-1.0) {}
            do_low_rank_inplace(const double& eps) : eps_(eps) {}

            //
            bool operator()(typename rangeT::iterator& it) const {

                nodeT& node = it->second;
				node.ftr2sr(eps_);
                return true;
            }
            template <typename Archive> void serialize(const Archive& ar) {}
        };

        template <typename opT>
        struct do_unary_op_value_inplace {
            typedef Range<typename dcT::iterator> rangeT;
            implT* impl;
            opT op;
            do_unary_op_value_inplace(implT* impl, const opT& op) : impl(impl), op(op) {}
            bool operator()(typename rangeT::iterator& it) const {
                const keyT& key = it->first;
                nodeT& node = it->second;
                if (node.has_coeff()) {
                	node.node_to_full_rank();
//                    coeffT& t= node.coeff();
                    tensorT& t= node.full_tensor_reference();
                    //double before = t.normf();
                    tensorT values = impl->fcube_for_mul(key, key, t);
                    op(key, values);
                    double scale = pow(0.5,0.5*NDIM*key.level())*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
                    t = transform(values,impl->cdata.quad_phiw).scale(scale);
                    node.node_to_low_rank();
                    //double after = t.normf();
                    //madness::print("XOP:", key, before, after);
                }
                return true;
            }
            template <typename Archive> void serialize(const Archive& ar) {}
        };

        template <typename Q, typename R>
        Void vtransform_doit(const std::shared_ptr< FunctionImpl<R,NDIM> >& right,
                             const Tensor<Q>& c,
                             const std::vector< std::shared_ptr< FunctionImpl<T,NDIM> > >& vleft,
                             double tol) {

            // To reduce crunch on vectors being transformed each task
            // does them in a random order
            std::vector<unsigned int> ind(vleft.size());
            for (unsigned int i=0; i<vleft.size(); ++i) {
                ind[i] = i;
            }
            for (unsigned int i=0; i<vleft.size(); ++i) {
                unsigned int j = RandomValue<int>()%vleft.size();
                std::swap(ind[i],ind[j]);
            }

            typename FunctionImpl<R,NDIM>::dcT::const_iterator end = right->coeffs.end();
            for (typename FunctionImpl<R,NDIM>::dcT::const_iterator it=right->coeffs.begin(); it != end; ++it) {
                if (it->second.has_coeff()) {
                    const Key<NDIM>& key = it->first;
                    const __Tensor<R>& r = it->second.coeff();
                    double norm = r.normf();
                    double keytol = truncate_tol(tol,key);

                    for (unsigned int j=0; j<vleft.size(); ++j) {
                        unsigned int i = ind[j]; // Random permutation
                        if (std::abs(norm*c(i)) > keytol) {
                            implT* left = vleft[i].get();
                            typename dcT::accessor acc;
                            bool newnode = left->coeffs.insert(acc,key);
                            if (newnode && key.level()>0) {
                                Key<NDIM> parent = key.parent();
                                left->coeffs.task(parent, &nodeT::set_has_children_recursive, left->coeffs, parent);
                            }
                            nodeT& node = acc->second;
                            if (!node.has_coeff())
                                node.set_coeff(coeffT(cdata.v2k,targs));
                            coeffT& t = node.coeff();
                            t.gaxpy(1.0, r, c(i));
                        }
                    }
                }
            }
            return None;
        }

        /// Transforms a vector of functions left[i] = sum[j] right[j]*c[j,i] using sparsity
        template <typename Q, typename R>
        void vtransform(const std::vector< std::shared_ptr< FunctionImpl<R,NDIM> > >& vright,
                        const Tensor<Q>& c,
                        const std::vector< std::shared_ptr< FunctionImpl<T,NDIM> > >& vleft,
                        double tol,
                        bool fence) {
            for (unsigned int j=0; j<vright.size(); ++j) {
                woT::task(world.rank(), &implT:: template vtransform_doit<Q,R>, vright[j], copy(c(j,_)), vleft, tol);
            }
            if (fence)
                world.gop.fence();
        }

        /// Unary operation applied inplace to the values with optional refinement and fence
        template <typename opT>
        void unary_op_value_inplace(const opT& op, bool fence) {
            typedef Range<typename dcT::iterator> rangeT;
            typedef do_unary_op_value_inplace<opT> xopT;
            world.taskq.for_each<rangeT,xopT>(rangeT(coeffs.begin(), coeffs.end()), xopT(this,op));
            if (fence)
                world.gop.fence();
        }

        // Multiplication assuming same distribution and recursive descent
        template <typename L, typename R>
        Void mulXXveca(const keyT& key,
                       const FunctionImpl<L,NDIM>* left, const Tensor<L>& lcin,
                       const std::vector<const FunctionImpl<R,NDIM>*> vrightin,
                       const std::vector< Tensor<R> >& vrcin,
                       const std::vector<FunctionImpl<T,NDIM>*> vresultin,
                       double tol) {
            typedef typename FunctionImpl<L,NDIM>::dcT::const_iterator literT;
            typedef typename FunctionImpl<R,NDIM>::dcT::const_iterator riterT;

            double lnorm = 1e99;
            Tensor<L> lc = lcin;
            if (lc.size() == 0) {
                literT it = left->coeffs.find(key).get();
                MADNESS_ASSERT(it != left->coeffs.end());
                lnorm = it->second.get_norm_tree();
                if (it->second.has_coeff())
                    lc = it->second.full_tensor_copy();
            }

            // Loop thru RHS functions seeing if anything can be multiplied
            std::vector<FunctionImpl<T,NDIM>*> vresult;
            std::vector<const FunctionImpl<R,NDIM>*> vright;
            std::vector< Tensor<R> > vrc;
            vresult.reserve(vrightin.size());
            vright.reserve(vrightin.size());
            vrc.reserve(vrightin.size());

            for (unsigned int i=0; i<vrightin.size(); ++i) {
                FunctionImpl<T,NDIM>* result = vresultin[i];
                const FunctionImpl<R,NDIM>* right = vrightin[i];
                Tensor<R> rc = vrcin[i];
                double rnorm;
                if (rc.size() == 0) {
                    riterT it = right->coeffs.find(key).get();
                    MADNESS_ASSERT(it != right->coeffs.end());
                    rnorm = it->second.get_norm_tree();
                    if (it->second.has_coeff())
                        rc = it->second.full_tensor_copy();
                }
                else {
                    rnorm = rc.normf();
                }

                if (rc.size() && lc.size()) { // Yipee!
                    result->task(world.rank(), &implT:: template do_mul<L,R>, key, lc, std::make_pair(key,rc));
                }
                else if (tol && lnorm*rnorm < truncate_tol(tol, key)) {
                    result->coeffs.replace(key, nodeT(coeffT(cdata.vk,targs),false)); // Zero leaf
                }
                else {
                    result->coeffs.replace(key, nodeT(coeffT(),true)); // Interior node
                    vresult.push_back(result);
                    vright.push_back(right);
                    vrc.push_back(rc);
                }
            }

            if (vresult.size()) {
                Tensor<L> lss;
                if (lc.size()) {
                    Tensor<L> ld(cdata.v2k);
                    ld(cdata.s0) = lc(___);
                    lss = left->unfilter(ld);
                }

                std::vector< Tensor<R> > vrss(vresult.size());
                for (unsigned int i=0; i<vresult.size(); ++i) {
                    if (vrc[i].size()) {
                        Tensor<R> rd(cdata.v2k);
                        rd(cdata.s0) = vrc[i](___);
                        vrss[i] = vright[i]->unfilter(rd);
                    }
                }

                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT& child = kit.key();
                    Tensor<L> ll;

                    std::vector<Slice> cp = child_patch(child);

                    if (lc.size())
                        ll = copy(lss(cp));

                    std::vector< Tensor<R> > vv(vresult.size());
                    for (unsigned int i=0; i<vresult.size(); ++i) {
                        if (vrc[i].size())
                            vv[i] = copy(vrss[i](cp));
                    }

                    woT::task(coeffs.owner(child), &implT:: template mulXXveca<L,R>, child, left, ll, vright, vv, vresult, tol);
                }
            }
            return None;
        }

        // Multiplication using recursive descent and assuming same distribution
        template <typename L, typename R>
        Void mulXXa(const keyT& key,
                    const FunctionImpl<L,NDIM>* left, const Tensor<L>& lcin,
                    const FunctionImpl<R,NDIM>* right,const Tensor<R>& rcin,
                    double tol) {
            typedef typename FunctionImpl<L,NDIM>::dcT::const_iterator literT;
            typedef typename FunctionImpl<R,NDIM>::dcT::const_iterator riterT;

            double lnorm=1e99, rnorm=1e99;

            Tensor<L> lc = lcin;
            if (lc.size() == 0) {
                literT it = left->coeffs.find(key).get();
                MADNESS_ASSERT(it != left->coeffs.end());
                lnorm = it->second.get_norm_tree();
                if (it->second.has_coeff())
                    lc = it->second.full_tensor_copy();
            }

            Tensor<R> rc = rcin;
            if (rc.size() == 0) {
                riterT it = right->coeffs.find(key).get();
                MADNESS_ASSERT(it != right->coeffs.end());
                rnorm = it->second.get_norm_tree();
                if (it->second.has_coeff())
                    rc = it->second.full_tensor_copy();
            }

            // both nodes are leaf nodes: multiply and return
            if (rc.size() && lc.size()) { // Yipee!
                do_mul<L,R>(key, lc, std::make_pair(key,rc));
                return None;
            }

            if (tol) {
                if (lc.size())
                    lnorm = lc.normf(); // Otherwise got from norm tree above
                if (rc.size())
                    rnorm = rc.normf();
                if (lnorm*rnorm < truncate_tol(tol, key)) {
                    coeffs.replace(key, nodeT(coeffT(cdata.vk,targs),false)); // Zero leaf node
                    return None;
                }
            }

            // Recur down
            coeffs.replace(key, nodeT(coeffT(),true)); // Interior node

            Tensor<L> lss;
            if (lc.size()) {
                Tensor<L> ld(cdata.v2k);
                ld(cdata.s0) = lc(___);
                lss = left->unfilter(ld);
            }

            Tensor<R> rss;
            if (rc.size()) {
                Tensor<R> rd(cdata.v2k);
                rd(cdata.s0) = rc(___);
                rss = right->unfilter(rd);
            }

            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                const keyT& child = kit.key();
                Tensor<L> ll;
                Tensor<R> rr;
                if (lc.size())
                    ll = copy(lss(child_patch(child)));
                if (rc.size())
                    rr = copy(rss(child_patch(child)));

                woT::task(coeffs.owner(child), &implT:: template mulXXa<L,R>, child, left, ll, right, rr, tol);
            }

            return None;
        }

        // Binary operation on values using recursive descent and assuming same distribution
        template <typename L, typename R, typename opT>
        Void binaryXXa(const keyT& key,
                       const FunctionImpl<L,NDIM>* left, const Tensor<L>& lcin,
                       const FunctionImpl<R,NDIM>* right,const Tensor<R>& rcin,
                       const opT& op) {
            typedef typename FunctionImpl<L,NDIM>::dcT::const_iterator literT;
            typedef typename FunctionImpl<R,NDIM>::dcT::const_iterator riterT;

            Tensor<L> lc = lcin;
            if (lc.size() == 0) {
                literT it = left->coeffs.find(key).get();
                MADNESS_ASSERT(it != left->coeffs.end());
                if (it->second.has_coeff())
                    lc = it->second.full_tensor_copy();
            }

            Tensor<R> rc = rcin;
            if (rc.size() == 0) {
                riterT it = right->coeffs.find(key).get();
                MADNESS_ASSERT(it != right->coeffs.end());
                if (it->second.has_coeff())
                    rc = it->second.full_tensor_copy();
            }

            if (rc.size() && lc.size()) { // Yipee!
                do_binary_op<L,R>(key, lc, std::make_pair(key,rc), op);
                return None;
            }

            // Recur down
            coeffs.replace(key, nodeT(coeffT(),true)); // Interior node

            Tensor<L> lss;
            if (lc.size()) {
                Tensor<L> ld(cdata.v2k);
                ld(cdata.s0) = lc(___);
                lss = left->unfilter(ld);
            }

            Tensor<R> rss;
            if (rc.size()) {
                Tensor<R> rd(cdata.v2k);
                rd(cdata.s0) = rc(___);
                rss = right->unfilter(rd);
            }

            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                const keyT& child = kit.key();
                Tensor<L> ll;
                Tensor<R> rr;
                if (lc.size())
                    ll = copy(lss(child_patch(child)));
                if (rc.size())
                    rr = copy(rss(child_patch(child)));

                woT::task(coeffs.owner(child), &implT:: template binaryXXa<L,R,opT>, child, left, ll, right, rr, op);
            }

            return None;
        }

        template <typename Q, typename opT>
        struct coeff_value_adaptor {
            typedef typename opT::resultT resultT;
            const FunctionImpl<Q,NDIM>* impl_func;
            opT op;

            coeff_value_adaptor() {};
            coeff_value_adaptor(const FunctionImpl<Q,NDIM>* impl_func,
                                const opT& op)
                    : impl_func(impl_func), op(op) {}

            Tensor<resultT> operator()(const Key<NDIM>& key, const Tensor<Q>& t) const {
                Tensor<Q> invalues = impl_func->coeffs2values(key, t);

                Tensor<resultT> outvalues = op(key, invalues);

                return impl_func->values2coeffs(key, outvalues);
            }

            template <typename Archive>
            void serialize(Archive& ar) {
                ar & impl_func & op;
            }
        };

        // Multiplication using recursive descent and assuming same distribution
        template <typename Q, typename opT>
        Void unaryXXa(const keyT& key,
                      const FunctionImpl<Q,NDIM>* func, const opT& op) {

//            const Tensor<Q>& fc = func->coeffs.find(key).get()->second.full_tensor_copy();
        	const Tensor<Q> fc = func->coeffs.find(key).get()->second.full_tensor_copy();

            if (fc.size() == 0) {
                // Recur down
                coeffs.replace(key, nodeT(coeffT(),true)); // Interior node
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT& child = kit.key();
                    woT::task(coeffs.owner(child), &implT:: template unaryXXa<Q,opT>, child, func, op);
                }
            }
            else {
            	tensorT t=op(key,fc);
                coeffs.replace(key, nodeT(coeffT(t,targs),false)); // Leaf node
            }

            return None;
        }

        template <typename L, typename R>
        void mulXX(const FunctionImpl<L,NDIM>* left, const FunctionImpl<R,NDIM>* right, double tol, bool fence) {
            if (world.rank() == coeffs.owner(cdata.key0))
                mulXXa(cdata.key0, left, Tensor<L>(), right, Tensor<R>(), tol);
            if (fence)
                world.gop.fence();

            //verify_tree();
        }

        template <typename L, typename R, typename opT>
        void binaryXX(const FunctionImpl<L,NDIM>* left, const FunctionImpl<R,NDIM>* right,
                      const opT& op, bool fence) {
            if (world.rank() == coeffs.owner(cdata.key0))
                binaryXXa(cdata.key0, left, Tensor<L>(), right, Tensor<R>(), op);
            if (fence)
                world.gop.fence();

            //verify_tree();
        }

        template <typename Q, typename opT>
        void unaryXX(const FunctionImpl<Q,NDIM>* func, const opT& op, bool fence) {
            if (world.rank() == coeffs.owner(cdata.key0))
                unaryXXa(cdata.key0, func, op);
            if (fence)
                world.gop.fence();

            //verify_tree();
        }

        template <typename Q, typename opT>
        void unaryXXvalues(const FunctionImpl<Q,NDIM>* func, const opT& op, bool fence) {
            if (world.rank() == coeffs.owner(cdata.key0))
                unaryXXa(cdata.key0, func, coeff_value_adaptor<Q,opT>(func,op));
            if (fence)
                world.gop.fence();

            //verify_tree();
        }

        template <typename L, typename R>
        void mulXXvec(const FunctionImpl<L,NDIM>* left,
                      const std::vector<const FunctionImpl<R,NDIM>*>& vright,
                      const std::vector<FunctionImpl<T,NDIM>*>& vresult,
                      double tol,
                      bool fence) {
            std::vector< Tensor<R> > vr(vright.size());
            if (world.rank() == coeffs.owner(cdata.key0))
                mulXXveca(cdata.key0, left, Tensor<L>(), vright, vr, vresult, tol);
            if (fence)
                world.gop.fence();
        }

        Future<double> get_norm_tree_recursive(const keyT& key) const;

        mutable long box_leaf[1000];
        mutable long box_interior[1000];

        // horrifically non-scalable
        Void put_in_box(ProcessID from, long nl, long ni) const {
            if (world.size()> 1000)
                throw "NO!";
            box_leaf[from] = nl;
            box_interior[from] = ni;
            return None;
        }

        /// Prints summary of data distribution
        void print_info() const {
            if (world.size() >= 1000)
                return;
            for (int i=0; i<world.size(); ++i)
                box_leaf[i] = box_interior[i] == 0;
            world.gop.fence();
            long nleaf=0, ninterior=0;
            typename dcT::const_iterator end = coeffs.end();
            for (typename dcT::const_iterator it=coeffs.begin(); it!=end; ++it) {
                const nodeT& node = it->second;
                if (node.is_leaf())
                    ++nleaf;
                else
                    ++ninterior;
            }
            this->send(0, &implT::put_in_box, world.rank(), nleaf, ninterior);
            world.gop.fence();
            if (world.rank() == 0) {
                for (int i=0; i<world.size(); ++i) {
                    printf("load: %5d %8ld %8ld\n", i, box_leaf[i], box_interior[i]);
                }
            }
            world.gop.fence();
        }

        /// Verify tree is properly constructed ... global synchronization involved

        /// If an inconsistency is detected, prints a message describing the error and
        /// then throws a madness exception.
        ///
        /// This is a reasonably quick and scalable operation that is
        /// useful for debugging and paranoia.
        void verify_tree() const;

        /// Walk up the tree returning pair(key,node) for first node with coefficients

        /// Three possibilities.
        ///
        /// 1) The coeffs are present and returned with the key of the containing node.
        ///
        /// 2) The coeffs are further up the tree ... the request is forwarded up.
        ///
        /// 3) The coeffs are futher down the tree ... an empty tensor is returned.
        ///
        /// !! This routine is crying out for an optimization to
        /// manage the number of messages being sent ... presently
        /// each parent is fetched 2^(n*d) times where n is the no. of
        /// levels between the level of evaluation and the parent.
        /// Alternatively, reimplement multiply as a downward tree
        /// walk and just pass the parent down.  Slightly less
        /// parallelism but much less communication.
        Void sock_it_to_me(const keyT& key,
                           const RemoteReference< FutureImpl< std::pair<keyT,coeffT> > >& ref) const;
        /// As above, except
        /// 3) The coeffs are constructed from the avg of nodes further down the tree
        Void sock_it_to_me_too(const keyT& key,
                               const RemoteReference< FutureImpl< std::pair<keyT,coeffT> > >& ref) const;

        Void plot_cube_kernel(archive::archive_ptr< Tensor<T> > ptr,
                              const keyT& key,
                              const coordT& plotlo, const coordT& plothi, const std::vector<long>& npt,
                              bool eval_refine) const;


        /// Evaluate a cube/slice of points ... plotlo and plothi are already in simulation coordinates

        /// No communications
        Tensor<T> eval_plot_cube(const coordT& plotlo,
                                 const coordT& plothi,
                                 const std::vector<long>& npt,
                                 const bool eval_refine = false) const;


        /// Evaluate function only if point is local returning (true,value); otherwise return (false,0.0)

        /// maxlevel is the maximum depth to search down to --- the max local depth can be
        /// computed with max_local_depth();
        std::pair<bool,T> eval_local_only(const Vector<double,NDIM>& xin, Level maxlevel) ;


        /// Evaluate the function at a point in \em simulation coordinates

        /// Only the invoking process will get the result via the
        /// remote reference to a future.  Active messages may be sent
        /// to other nodes.
        Void eval(const Vector<double,NDIM>& xin,
                  const keyT& keyin,
                  const typename Future<T>::remote_refT& ref);

        /// Get the depth of the tree at a point in \em simulation coordinates

        /// Only the invoking process will get the result via the
        /// remote reference to a future.  Active messages may be sent
        /// to other nodes.
        ///
        /// This function is a minimally-modified version of eval()
        Void evaldepthpt(const Vector<double,NDIM>& xin,
                  const keyT& keyin,
                  const typename Future<Level>::remote_refT& ref);

        /// Get the rank of leaf box of the tree at a point in \em simulation coordinates

        /// Only the invoking process will get the result via the
        /// remote reference to a future.  Active messages may be sent
        /// to other nodes.
        ///
        /// This function is a minimally-modified version of eval()
        Void evalR(const Vector<double,NDIM>& xin,
                  const keyT& keyin,
                  const typename Future<long>::remote_refT& ref);


        /// Computes norm of low/high-order polyn. coeffs for autorefinement test

        /// t is a k^d tensor.  In order to screen the autorefinement
        /// during multiplication compute the norms of
        /// ... lo ... the block of t for all polynomials of order < k/2
        /// ... hi ... the block of t for all polynomials of order >= k/2
        ///
        /// k=5   0,1,2,3,4     --> 0,1,2 ... 3,4
        /// k=6   0,1,2,3,4,5   --> 0,1,2 ... 3,4,5
        ///
        /// k=number of wavelets, so k=5 means max order is 4, so max exactly
        /// representable squarable polynomial is of order 2.
        void tnorm(const tensorT& t, double* lo, double* hi) const;

        // This invoked if node has not been autorefined
        Void do_square_inplace(const keyT& key);

        // This invoked if node has been autorefined
        Void do_square_inplace2(const keyT& parent, const keyT& child, const tensorT& parent_coeff);

        /// Always returns false (for when autorefine is not wanted)
        bool noautorefine(const keyT& key, const tensorT& t) const {
            return false;
        }

        /// Returns true if this block of coeffs needs autorefining
        bool autorefine_square_test(const keyT& key, const nodeT& t) const {
            double lo, hi;
            tnorm(t.full_tensor_copy(), &lo, &hi);
            double test = 2*lo*hi + hi*hi;
            //print("autoreftest",key,thresh,truncate_tol(thresh, key),lo,hi,test);
            return test> truncate_tol(thresh, key);
        }

        /// Pointwise squaring of function with optional global fence

        /// If not autorefining, local computation only if not fencing.
        /// If autorefining, may result in asynchronous communication.
        void square_inplace(bool fence);
        void abs_inplace(bool fence);
        void abs_square_inplace(bool fence);

        Void sum_down_spawn(const keyT& key, const coeffT& s) {
            typename dcT::accessor acc;
            coeffs.insert(acc,key);
            nodeT& node = acc->second;
            coeffT& c = node.coeff();

            //print(key,"received",s.normf(),c.normf(),node.has_children());

            if (s.size() > 0) {
                if (c.size() > 0)
                    c.gaxpy(1.0,s,1.0);
                else
                    c = s;
            }

            if (node.has_children()) {
                coeffT d;
                if (c.size() > 0) {
                    d = coeffT(cdata.v2k,targs);
                    d(cdata.s0) = c;
                    d = unfilter(d);
                    node.clear_coeff();
                }
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    coeffT ss;
                    const keyT& child = kit.key();
                    if (d.size() > 0) ss = copy(d(child_patch(child)));
                    //print(key,"sending",ss.normf(),"to",child);
                    woT::task(coeffs.owner(child), &implT::sum_down_spawn, child, ss);
                }
            }
            else {
                // Missing coeffs assumed to be zero
                if (c.size() <= 0) c = coeffT(cdata.vk,targs);
            }
            node.node_to_low_rank();
            return None;
        }

        /// After 1d push operator must sum coeffs down the tree to restore correct scaling function coefficients
        void sum_down(bool fence) {
            if (world.rank() == coeffs.owner(cdata.key0)) sum_down_spawn(cdata.key0, coeffT());

            if (fence) world.gop.fence();
        }

        /// given two functions of LDIM, perform the Hartree/Kronecker/outer product

        /// |Phi(1,2)> = |phi(1)> x |phi(2)>
        template<std::size_t LDIM>
        Void hartree_product(const FunctionImpl<T,LDIM>* p1, const FunctionImpl<T,LDIM>* p2, bool fence) {
        	MADNESS_ASSERT(p1->is_nonstandard());
        	MADNESS_ASSERT(p2->is_nonstandard());

			typedef std::pair<Key<LDIM>, FunctionNode<T,LDIM> > datumL;
			typedef FunctionImpl<T,LDIM> implL;

            if (world.rank() == coeffs.owner(cdata.key0)) {
            	Future<datumL> datum1=p1->task(p1->get_coeffs().owner(p1->cdata.key0), &FunctionImpl<T,LDIM>::find_datum,
            	        p1->cdata.key0,TaskAttributes::hipri());
            	Future<datumL> datum2=p2->task(p2->get_coeffs().owner(p2->cdata.key0), &FunctionImpl<T,LDIM>::find_datum,
            	        p2->cdata.key0,TaskAttributes::hipri());

            	ProcessID owner = coeffs.owner(cdata.key0);
            	woT::task(owner, &implT:: template hartree_product_spawn<LDIM>, p1, p2, cdata.key0,
            			datum1,datum2);
            }

        	this->compressed=false;
        	if (fence) world.gop.fence();
        	return None;
        }


        /// given two Nodes, perform the Hartree product and insert the result
        template<std::size_t LDIM>
        Void hartree_product_op(const FunctionImpl<T,LDIM>* p1, const FunctionImpl<T,LDIM>* p2,
        		const keyT& key,
        		const std::pair<Key<LDIM>, FunctionNode<T,LDIM> >& datum1,
        		const std::pair<Key<LDIM>, FunctionNode<T,LDIM> >& datum2) {

        	MADNESS_ASSERT(NDIM==6);
        	MADNESS_ASSERT(LDIM==3);

			// break key into particles (these are the child keys, with datum1/2 come the parent keys)
			const Vector<Translation, NDIM>& l=key.translation();
			const Vector<Translation, LDIM> l1=Vector<Translation,LDIM> (vec(l[0],l[1],l[2]));
			const Vector<Translation, LDIM> l2=Vector<Translation,LDIM> (vec(l[3],l[4],l[5]));
			const Key<LDIM> key1(key.level(),l1);
			const Key<LDIM> key2(key.level(),l2);

			// iterators point to nodes in nonstandard representation: get the sum coeffs
			const coeffT s1=datum1.second.coeff()(p1->cdata.s0);
			const coeffT s2=datum2.second.coeff()(p2->cdata.s0);

	        coeffT coeff1=p1->parent_to_child(s1,datum1.first,key1);
	        coeffT coeff2=p2->parent_to_child(s2,datum2.first,key2);

	        // new coeffs are simply the hartree/kronecker/outer product
	        tensorT tcube=outer(coeff1.full_tensor_copy(),coeff2.full_tensor_copy());

	        coeffs.replace(key,nodeT(coeffT(tcube,targs),false));	// leaf node w/o children

			return None;

        }


        /// Hartree product of two n-D functions to yield a 2n-D function

        /// @param[in]	p1	function of particle 1
        /// @param[in]	p2	function of particle 2
        /// @param[in]	it1	pointer to a valid node in p1
        /// @param[in]	it2	pointer to a valid node in p2
        /// @return		this
        template<std::size_t LDIM>
        Void hartree_product_spawn(const FunctionImpl<T,LDIM>* p1, const FunctionImpl<T,LDIM>* p2,
        		const Key<NDIM>& key,
        		const std::pair<Key<LDIM>, FunctionNode<T,LDIM> >& datum1,
        		const std::pair<Key<LDIM>, FunctionNode<T,LDIM> >& datum2) {

        	// for now
        	MADNESS_ASSERT(LDIM==3);
        	MADNESS_ASSERT(NDIM==6);

        	// for convenience
			typedef std::pair<Key<LDIM>, FunctionNode<T,LDIM> > datumL;
        	typedef FunctionNode<T,LDIM> nodeL;
        	typedef FunctionImpl<T,LDIM> implL;

        	const nodeL& node1=datum1.second;
        	const nodeL& node2=datum2.second;

            // if the final norm is small, perform the hartree product and return
        	const coeffT s1=node1.coeff()(p1->cdata.s0);
        	const coeffT s2=node2.coeff()(p2->cdata.s0);
            const double norm1=s1.normf();
            const double norm2=s2.normf();
            const double norm=norm1*norm2;	// computing the outer product
            if (norm < truncate_tol(thresh, key)) {

				ProcessID owner = coeffs.owner(key);
				woT::task(owner, &implT:: template hartree_product_op<LDIM>, p1, p2, key, datum1, datum2);
				return None;
            }

            // get the error of both functions and of the pair function
            const double error1=p1->compute_error(node1);
            const double error2=p2->compute_error(node2);
            const double error=sqrt(norm1*norm1*error2*error2 + error1*error1*norm2*norm2 + error1*error1*error2*error2);

            // if the expected error is small, perform the hartree product and return
            if (error < truncate_tol(thresh,key)) {
				// should be local anyways
				ProcessID owner = coeffs.owner(key);
				woT::task(owner, &implT:: template hartree_product_op<LDIM>, p1, p2, key, datum1, datum2);
				return None;
            }

            // norm and error are large: keep recurring
			for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {

				const keyT& child = kit.key();

				// break key into particles
				const Vector<Translation, NDIM> l=child.translation();
				const Vector<Translation, LDIM> l1=Vector<Translation,LDIM> (vec(l[0],l[1],l[2]));
				const Vector<Translation, LDIM> l2=Vector<Translation,LDIM> (vec(l[3],l[4],l[5]));
				const Key<LDIM> key1(child.level(),l1);
				const Key<LDIM> key2(child.level(),l2);

                // point to "outermost" leaf node
				Future<datumL> datum11, datum22;
				if (node1.has_children()) {
					datum11=p1->task(p1->coeffs.owner(key1), &FunctionImpl<T,LDIM>::find_datum, key1,
						TaskAttributes::hipri());
				} else {
					datum11.set(datum1);
				}

				if (node2.has_children()) {
					datum22=p2->task(p2->coeffs.owner(key2), &FunctionImpl<T,LDIM>::find_datum, key2,
						TaskAttributes::hipri());
				} else {
					datum22.set(datum2);
				}

			woT::task(world.rank(), &implT:: template hartree_product_spawn2<LDIM>, p1, p2, child,
                        datum11, datum22);
			}

            coeffs.replace(key, nodeT(coeffT(),true));  // empty internal node w/ children
			return None;

        }

        /// work-around for serializing the Future<datumL>
        template<std::size_t LDIM>
        Void hartree_product_spawn2(const FunctionImpl<T,LDIM>* p1, const FunctionImpl<T,LDIM>* p2,
                        const Key<NDIM>& key,
                        const std::pair<Key<LDIM>, FunctionNode<T,LDIM> >& datum1,
                        const std::pair<Key<LDIM>, FunctionNode<T,LDIM> >& datum2) {
                ProcessID owner = coeffs.owner(key);
                woT::task(owner, &implT:: template hartree_product_spawn<LDIM>, p1, p2, key,
                        datum1, datum2,TaskAttributes::hipri());
                return None;
        }


        template <typename opT, typename R>
        Void
        apply_1d_realspace_push_op(const archive::archive_ptr<const opT>& pop, int axis, const keyT& key, const Tensor<R>& c) {
            const opT* op = pop.ptr;
            const Level n = key.level();
            const double cnorm = c.normf();
            const double tol = truncate_tol(thresh, key)*0.1; // ??? why this value????

            Vector<Translation,NDIM> lnew(key.translation());
            const Translation lold = lnew[axis];
            const Translation maxs = Translation(1)<<n;

            int nsmall = 0; // Counts neglected blocks to terminate s loop
            for (Translation s=0; s<maxs; ++s) {
                int maxdir = s ? 1 : -1;
                for (int direction=-1; direction<=maxdir; direction+=2) {
                    lnew[axis] = lold + direction*s;
                    if (lnew[axis] >= 0 && lnew[axis] < maxs) { // NON-ZERO BOUNDARY CONDITIONS IGNORED HERE !!!!!!!!!!!!!!!!!!!!
                        const Tensor<typename opT::opT>& r = op->rnlij(n, s*direction, true);
                        double Rnorm = r.normf();

                        if (Rnorm == 0.0) {
                            return None; // Hard zero means finished!
                        }

                        if (s <= 1  ||  r.normf()*cnorm > tol) { // Always do kernel and neighbor
                            nsmall = 0;
                            tensorT result = transform_dir(c,r,axis);

                            if (result.normf() > tol*0.3) {
                                Key<NDIM> dest(n,lnew);
                                coeffs.task(dest, &nodeT::accumulate2, result, coeffs, dest, TaskAttributes::hipri());
                            }
                        }
                        else {
                            ++nsmall;
                        }
                    }
                    else {
                        ++nsmall;
                    }
                }
                if (nsmall >= 4) {
                    // If have two negligble blocks in
                    // succession in each direction interpret
                    // this as the operator being zero beyond
                    break;
                }
            }
            return None;
        }

        template <typename opT, typename R>
        void
        apply_1d_realspace_push(const opT& op, const FunctionImpl<R,NDIM>* f, int axis, bool fence) {
            MADNESS_ASSERT(!f->is_compressed());

            typedef typename FunctionImpl<R,NDIM>::dcT::const_iterator fiterT;
            typedef FunctionNode<R,NDIM> fnodeT;
            fiterT end = f->coeffs.end();
            ProcessID me = world.rank();
            for (fiterT it=f->coeffs.begin(); it!=end; ++it) {
                const fnodeT& node = it->second;
                if (node.has_coeff()) {
                    const keyT& key = it->first;
                    const Tensor<R>& c = node.full_tensor_copy();
                    woT::task(me, &implT:: template apply_1d_realspace_push_op<opT,R>,
                    		archive::archive_ptr<const opT>(&op), axis, key, c);
                }
            }
            if (fence) world.gop.fence();
        }

        Void do_diff1(const DerivativeBase<T,NDIM>* D,
                      const implT* f,
                      const keyT& key,
                      const std::pair<keyT,coeffT>& left,
                      const std::pair<keyT,coeffT>& center,
                      const std::pair<keyT,coeffT>& right) {
            return D->do_diff1(f,this,key,left,center,right);
        }


        // Called by result function to differentiate f
        void diff(const DerivativeBase<T,NDIM>* D, const implT* f, bool fence) {
            typedef std::pair<keyT,coeffT> argT;
            typename dcT::const_iterator end = f->coeffs.end();
            for (typename dcT::const_iterator it=f->coeffs.begin(); it!=end; ++it) {
                const keyT& key = it->first;
                const nodeT& node = it->second;
                if (node.has_coeff()) {
                    Future<argT> left  = D->find_neighbor(f, key,-1);
                    argT center(key,node.coeff());
                    Future<argT> right = D->find_neighbor(f, key, 1);
                    woT::task(world.rank(), &implT::do_diff1, D, f, key, left, center, right, TaskAttributes::hipri());
                }
                else {
                    coeffs.replace(key,nodeT(coeffT(),true)); // Empty internal node
                }
            }
            if (fence) world.gop.fence();
        }

        /// Returns key of general neighbor enforcing BC

        /// Out of volume keys are mapped to enforce the BC as follows.
        ///   * Periodic BC map back into the volume and return the correct key
        ///   * Zero BC - returns invalid() to indicate out of volume
        keyT neighbor(const keyT& key, const keyT& disp, const std::vector<bool>& is_periodic) const;

        /// find_me. Called by diff_bdry to get coefficients of boundary function
        Future< std::pair<keyT,coeffT> > find_me(const keyT& key) const;

        /// return the a Future to a std::pair<key, node>, which MUST exist
        std::pair<Key<NDIM>,FunctionNode<T,NDIM> > find_datum(keyT key) const {
        	MADNESS_ASSERT(coeffs.probe(key));
        	return std::pair<Key<NDIM>,FunctionNode<T,NDIM> >(key,coeffs.find(key).get()->second);
        }



        /// construct the FunctionNodes of this on-demand FunctionImpl

        /// descend the tree of this until we find a leaf node of muster; there
        /// we construct a function node, and possibly refine by one scale
        /// @param[in]  muster  take this function as a template to fill this' nodes
        /// @param[in]  do_refine   if tnorm tells us that the coefficients provided by
        ///                         the functor are inaccurate, refine by one (1) scale
        /// @return     nothing, but fills this' leaf nodes with sum coefficients
        void fill_on_demand_tree(const implT* muster, const bool do_refine, const bool fence=true) {

        	// some checks
        	MADNESS_ASSERT(functor);
        	MADNESS_ASSERT(muster);

        	// get rid of all what we have right now
        	coeffs.clear();

        	const keyT& key0=cdata.key0;
            woT::task(coeffs.owner(key0),&implT::fill_coeff_spawn,muster,key0,do_refine);

            if (fence) world.gop.fence();
        }

        /// recursive part of the fill_on_demand_tree function

        /// @param[in]  muster  take this function as a template to fill this' nodes
        /// @param[in]  key     the current key we are working on
        /// @param[in]  do_refine   if tnorm tells us that the coefficients by functor are
        ///                         inaccurate, refine by one (1) scale
        /// @return     nothing, but fills this' leaf nodes with sum coefficients
        Void fill_coeff_spawn(const implT* muster, const keyT key, const bool do_refine) {

            // key of muster exists and is local
            MADNESS_ASSERT(muster->coeffs.probe(key));
            const nodeT node=muster->coeffs.find(key).get()->second;

            // if this is a leaf node of muster fill this's node
            // otherwise descend down the tree of muster
            if (node.is_leaf()) {

                // note that fill_coeff can also spawn more tasks on the children
                // if we have to refine
                functor->fill_coeff(this,key,do_refine);
            } else {

                // insert an empty internal node
                this->coeffs.replace(key, nodeT(coeffT(),true));

                // descend down the tree
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT& child = kit.key();
                    woT::task(coeffs.owner(child),&implT::fill_coeff_spawn,muster,child,do_refine);
                }
            }

            return None;
        }


        /// fill the children of this node with coefficients provided by its functor
        Void fill_coeff_refine(const keyT& key) {

            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                const keyT& child = kit.key();
                functor->fill_coeff(this,child,false);    // no further refinement
            }
            return None;
        }



        /// assemble the provided coefficients and insert them into the tree at key

        /// check for the autorefine criterion, and if it is not met, refine by 1 scale
        void assemble_coeff(const bool do_refine, const keyT& key,
            const coeffT& val_ket, const coeffT& val_eri,
            const coeffT& val_pot1, const coeffT& val_pot2) {

            woT::task(coeffs.owner(key),  &implT:: template do_assemble_coeff<3>,
                    do_refine,key,val_ket,val_eri,val_pot1,val_pot2);

        }

        /// assemble the provided coefficients and insert them into the tree at key

        /// check for the autorefine criterion, and if it is not met, refine by 1 scale
        template<size_t LDIM>
        Void do_assemble_coeff(const bool do_refine, const keyT& key,
                const coeffT& val_ket, const coeffT& val_eri,
                const coeffT& val_pot1, const coeffT& val_pot2) {


            // sort of makes sense
            MADNESS_ASSERT(val_ket.has_data());

            tensorT ket=val_ket.full_tensor_copy();
            tensorT coeff_v(cdata.vk);      // 6D, all addends

            // include the one-electron potential
            if (val_pot1.has_data() and val_pot2.has_data()) {

                FunctionCommonData<T,LDIM> cdataL=FunctionCommonData<T,LDIM>::get(k);
                Tensor<T> identity=Tensor<double>(cdataL.vk,false);
                identity=1.0;

                // direct product: V(1) * E(2) + E(1) * V(2)
                coeff_v = outer(val_pot1.full_tensor_copy(),identity)
                          + outer(identity,val_pot2.full_tensor_copy());
            } else if (val_pot1.has_data()) {

                MADNESS_ASSERT(LDIM==NDIM);
                coeff_v = (val_pot1.full_tensor_copy());
            }


            // add 2-particle contribution
            if (val_eri.has_data()) coeff_v+=val_eri.full_tensor_copy();

            // check for need of refinement
            bool needs_refinement=false;
            if (do_refine and (val_eri.has_data() or val_pot1.has_data())) {

                double lo_ket, hi_ket, lo_pot, hi_pot;
                tnorm(ket, &lo_ket, &hi_ket);
                tnorm(coeff_v, &lo_pot, &hi_pot);
                double test = lo_ket*hi_pot + lo_pot*hi_ket + hi_ket*hi_pot;
                //print("autoreftest",key,thresh,truncate_tol(thresh, key),lo,hi,test);
                needs_refinement= (test> truncate_tol(thresh, key));

            }

            if (needs_refinement and do_refine) {
                // insert empty internal node and refine one scale
                coeffs.replace(key,nodeT(coeffT(),true));   // empty internal node
                this->fill_coeff_refine(key);

            } else {
                // multiply potential with ket
                if (val_eri.has_data() or val_pot1.has_data()) {

                    tensorT tcube(cdata.vk,false);

                    TERNARY_OPTIMIZED_ITERATOR(T, tcube, T, ket, T, coeff_v, *_p0 = *_p1 * *_p2;);
                    ket=tcube;
                }

                coeffs.replace(key,nodeT(coeffT(values2coeffs(key,ket),targs)));  // leaf node
            }
            return None;
        }

        /// Permute the dimensions according to map
        void mapdim(const implT& f, const std::vector<long>& map, bool fence);

        T eval_cube(Level n, coordT& x, const tensorT& c) const;

        /// Transform sum coefficients at level n to sums+differences at level n-1

        /// Given scaling function coefficients s[n][l][i] and s[n][l+1][i]
        /// return the scaling function and wavelet coefficients at the
        /// coarser level.  I.e., decompose Vn using Vn = Vn-1 + Wn-1.
        /// \code
        /// s_i = sum(j) h0_ij*s0_j + h1_ij*s1_j
        /// d_i = sum(j) g0_ij*s0_j + g1_ij*s1_j
        //  \endcode
        /// Returns a new tensor and has no side effects.  Works for any
        /// number of dimensions.
        ///
        /// No communication involved.
        tensorT filter(const tensorT& s) const {
            tensorT r(cdata.v2k,false);
            tensorT w(cdata.v2k,false);
            return fast_transform(s,cdata.hgT,r,w);
            //return transform(s,cdata.hgT);
        }

#if HAVE_GENTENSOR
        coeffT filter(const coeffT& s) const {
            coeffT result=transform(s,cdata.hgT);
        	return result;
        }
#endif

        ///  Transform sums+differences at level n to sum coefficients at level n+1

        ///  Given scaling function and wavelet coefficients (s and d)
        ///  returns the scaling function coefficients at the next finer
        ///  level.  I.e., reconstruct Vn using Vn = Vn-1 + Wn-1.
        ///  \code
        ///  s0 = sum(j) h0_ji*s_j + g0_ji*d_j
        ///  s1 = sum(j) h1_ji*s_j + g1_ji*d_j
        ///  \endcode
        ///  Returns a new tensor and has no side effects
        ///
        ///  If (sonly) ... then ss is only the scaling function coeff (and
        ///  assume the d are zero).  Works for any number of dimensions.
        ///
        /// No communication involved.
        tensorT unfilter(const tensorT& s) const {
            tensorT r(cdata.v2k,false);
            tensorT w(cdata.v2k,false);
            return fast_transform(s,cdata.hg,r,w);
            //return transform(s, cdata.hg);
        }

#if HAVE_GENTENSOR
        coeffT unfilter(const coeffT& s) const {
            return transform(s,cdata.hg);
        }
#endif
        /// Projects old function into new basis (only in reconstructed form)
        void project(const implT& old, bool fence) {
            long kmin = std::min(cdata.k,old.cdata.k);
            std::vector<Slice> s(NDIM,Slice(0,kmin-1));
            typename dcT::const_iterator end = old.coeffs.end();
            for (typename dcT::const_iterator it=old.coeffs.begin(); it!=end; ++it) {
                const keyT& key = it->first;
                const nodeT& node = it->second;
                if (node.has_coeff()) {
                    coeffT c(cdata.vk,targs);
                    c(s) += node.coeff()(s);
                    coeffs.replace(key,nodeT(c,false));
                }
                else {
                    coeffs.replace(key,nodeT(coeffT(),true));
                }
            }
            if (fence)
                world.gop.fence();
        }

        struct true_refine_test {
            bool operator()(const implT* f, const keyT& key, const nodeT& t) const {
                return true;
            }
            template <typename Archive> void serialize(Archive& ar) {}
        };

        template <typename opT>
        Void refine_op(const opT& op, const keyT& key) {
            // Must allow for someone already having autorefined the coeffs
            // and we get a write accessor just in case they are already executing
            typename dcT::accessor acc;
            MADNESS_ASSERT(coeffs.find(acc,key));
            nodeT& node = acc->second;
            if (node.has_coeff() && key.level() < max_refine_level && op(this, key, node)) {
                coeffT d(cdata.v2k,targs);
                d(cdata.s0) += copy(node.coeff());
                d = unfilter(d);
                node.clear_coeff();
                node.set_has_children(true);
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT& child = kit.key();
                    coeffT ss = copy(d(child_patch(child)));
                    ss.reduceRank(targs.thresh);
//                    coeffs.replace(child,nodeT(ss,-1.0,false).node_to_low_rank());
                    coeffs.replace(child,nodeT(ss,-1.0,false));
                    // Note value -1.0 for norm tree to indicate result of refinement
                }
            }
            return None;
        }

        template <typename opT>
        Void refine_spawn(const opT& op, const keyT& key) {
            nodeT& node = coeffs.find(key).get()->second;
            if (node.has_children()) {
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit)
                    woT::task(coeffs.owner(kit.key()), &implT:: template refine_spawn<opT>, op, kit.key(), TaskAttributes::hipri());
            }
            else {
                woT::task(coeffs.owner(key), &implT:: template refine_op<opT>, op, key);
            }
            return None;
        }

        // Refine in real space according to local user-defined criterion
        template <typename opT>
        void refine(const opT& op, bool fence) {
            if (world.rank() == coeffs.owner(cdata.key0))
                woT::task(coeffs.owner(cdata.key0), &implT:: template refine_spawn<opT>, op, cdata.key0, TaskAttributes::hipri());
            if (fence)
                world.gop.fence();
        }

        bool exists_and_has_children(const keyT& key) const {
            return coeffs.probe(key) && coeffs.find(key).get()->second.has_children();
        }

        bool exists_and_is_leaf(const keyT& key) const {
            return coeffs.probe(key) && (not coeffs.find(key).get()->second.has_children());
        }


        Void broaden_op(const keyT& key, const std::vector< Future <bool> >& v) {
            for (unsigned int i=0; i<v.size(); ++i) {
                if (v[i]) {
                    refine_op(true_refine_test(), key);
                    break;
                }
            }
            return None;
        }

        // For each local node sets value of norm tree to 0.0
        void zero_norm_tree() {
            typename dcT::iterator end = coeffs.end();
            for (typename dcT::iterator it=coeffs.begin(); it!=end; ++it) {
                it->second.set_norm_tree(0.0);
            }
        }

        // Broaden tree
        void broaden(std::vector<bool> is_periodic, bool fence) {
            typename dcT::iterator end = coeffs.end();
            for (typename dcT::iterator it=coeffs.begin(); it!=end; ++it) {
            	const keyT& key = it->first;
            	typename dcT::accessor acc;
            	MADNESS_ASSERT(coeffs.find(acc,key));
            	nodeT& node = acc->second;
                if (node.has_coeff() &&
                    node.get_norm_tree() != -1.0 &&
                    node.coeff().normf() >= truncate_tol(thresh,key)) {

                    node.set_norm_tree(-1.0); // Indicates already broadened or result of broadening/refining

                    //int ndir = std::pow(3,NDIM);
                    int ndir = static_cast<int>(std::pow(static_cast<double>(3), static_cast<int>(NDIM)));
                    std::vector< Future <bool> > v = future_vector_factory<bool>(ndir);
                    keyT neigh;
                    int i=0;
                    for (HighDimIndexIterator it(NDIM,3); it; ++it) {
                        Vector<Translation,NDIM> l(*it);
                        for (std::size_t d=0; d<NDIM; ++d) {
                            const int odd = key.translation()[d] & 0x1L; // 1 if odd, 0 if even
                            l[d] -= 1; // (0,1,2) --> (-1,0,1)
                            if (l[d] == -1)
                                l[d] = -1-odd;
                            else if (l[d] ==  1)
                                l[d] = 2 - odd;
                        }
                        keyT neigh = neighbor(key, keyT(key.level(),l), is_periodic);

                        if (neigh.is_valid()) {
                            v[i++] = this->send(coeffs.owner(neigh), &implT::exists_and_has_children, neigh);
                        }
                        else {
                            v[i++].set(false);
                        }
                    }
                    woT::task(world.rank(), &implT::broaden_op, key, v);
                }
            }
            // Reset value of norm tree so that can repeat broadening
            if (fence) {
                world.gop.fence();
                zero_norm_tree();
                world.gop.fence();
            }
        }

        void reconstruct(bool fence) {
            // Must set true here so that successive calls without fence do the right thing
            nonstandard = compressed = redundant = false;
            if (world.rank() == coeffs.owner(cdata.key0))
                woT::task(world.rank(), &implT::reconstruct_op, cdata.key0,coeffT());
            if (fence)
                world.gop.fence();
        }

        // Invoked on node where key is local
//        Void reconstruct_op(const keyT& key, const tensorT& s);
        Void reconstruct_op(const keyT& key, const coeffT& s);

        /// compress the wave function

        /// after application there will be sum coefficients at the root level,
        /// and difference coefficients at all other levels; furthermore:
        /// @param[in] nonstandard	keep sum coeffs at all other levels, except leaves
        /// @param[in] keepleaves	keep sum coeffs (but no diff coeffs) at leaves
        void compress(bool nonstandard, bool keepleaves, bool fence) {
            // Must set true here so that successive calls without fence do the right thing
            this->compressed = true;
            this->nonstandard = nonstandard;
            this->redundant = false;

//            this->print_tree();
            if (world.rank() == coeffs.owner(cdata.key0)) {

           		compress_spawn(cdata.key0, nonstandard, keepleaves);
            }
            if (fence)
                world.gop.fence();
        }

        // Invoked on node where key is local
//        Future<tensorT> compress_spawn(const keyT& key, bool nonstandard, bool keepleaves);
        Future<coeffT > compress_spawn(const keyT& key, bool nonstandard, bool keepleaves);

        void ftr2sr() {
    		if (world.rank()==0) print("ftr2sr on ",this->tree_size()," boxes");
    		flo_unary_op_node_inplace(do_low_rank_inplace(thresh),true);
        }

        /// convert this to redundant, i.e. have sum coefficients on all levels
        void make_redundant(const bool fence) {

        	if (is_redundant()) return;

        	// compress, but keep sum coefficients
        	if (not is_nonstandard()) {
        		if (is_compressed()) reconstruct(true);
        		compress(true,true,true);
        	}

        	// get rid of difference coefficients
        	redundant=true;
        	compressed=false;
        	nonstandard=false;
    		flo_unary_op_node_inplace(do_keep_sum_coeffs(this),true);

    		if (fence) world.gop.fence();
        }

        /// compute for each FunctionNode the norm of the function inside that node
        void norm_tree(bool fence) {
            if (world.rank() == coeffs.owner(cdata.key0))
                norm_tree_spawn(cdata.key0);
            if (fence)
                world.gop.fence();
        }

        double norm_tree_op(const keyT& key, const std::vector< Future<double> >& v) {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            double sum = 0.0;
            int i=0;
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
                double value = v[i].get();
                sum += value*value;
            }
            sum = sqrt(sum);
            coeffs.task(key, &nodeT::set_norm_tree, sum);
            //if (key.level() == 0) std::cout << "NORM_TREE_TOP " << sum << "\n";
            return sum;
        }

        Future<double> norm_tree_spawn(const keyT& key) {
            nodeT& node = coeffs.find(key).get()->second;
            if (node.has_children()) {
                std::vector< Future<double> > v = future_vector_factory<double>(1<<NDIM);
                int i=0;
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
                    v[i] = woT::task(coeffs.owner(kit.key()), &implT::norm_tree_spawn, kit.key());
                }
                return woT::task(world.rank(),&implT::norm_tree_op, key, v);
            }
            else {
//                return Future<double>(node.coeff().normf());
                const double norm=node.coeff().normf();
                // invoked locally anyways
            	node.set_norm_tree(norm);
                return Future<double>(norm);
            }
        }

#if 1
        /// calculate the wavelet coefficients using the sum coefficients of all child nodes

        /// @param[in] key 	this's key
        /// @param[in] v 	sum coefficients of the child nodes
        /// @return 		the sum coefficients
        coeffT compress_op(const keyT& key, const std::vector< Future<coeffT > >& v, bool nonstandard) {
            PROFILE_MEMBER_FUNC(FunctionImpl);

            // Copy child scaling coeffs into contiguous block
            tensorT d(cdata.v2k);
//            coeffT d(cdata.v2k,targs);
            int i=0;
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
//                d(child_patch(kit.key())) += v[i].get();
                d(child_patch(kit.key())) += v[i].get().full_tensor_copy();
            }

            d = filter(d);

            typename dcT::accessor acc;
            MADNESS_ASSERT(coeffs.find(acc, key));

            if (acc->second.has_coeff()) {
            	print(" stuff in compress_op");
//                const coeffT& c = acc->second.coeff();
                const tensorT c = acc->second.coeff().full_tensor_copy();
                if (c.dim(0) == k) {
                    d(cdata.s0) += c;
                }
                else {
                    d += c;
                }
            }

            // this is deep copy
//            coeffT s= copy(d(cdata.s0));
            tensorT s=copy(d(cdata.s0));

            if (key.level()> 0 && !nonstandard)
                d(cdata.s0) = 0.0;

//            d.reduceRank(thresh);
            coeffT dd=coeffT(d,targs);
            acc->second.set_coeff(dd);
//            s.reduceRank(thresh);
            coeffT ss=coeffT(s,targs);

            return ss;
        }
#else

        /// calculate the wavelet coefficients using the sum coefficients of all child nodes

        /// @param[in] key  this's key
        /// @param[in] v    sum coefficients of the child nodes
        /// @return         the sum coefficients
        tensorT compress_op(const keyT& key, const std::vector< Future<tensorT> >& v, bool nonstandard) {
            PROFILE_MEMBER_FUNC(FunctionImpl);

            // Copy child scaling coeffs into contiguous block
//            coeffT d(cdata.v2k,targs);
            tensorT d(cdata.v2k);
            int i=0;
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
                d(child_patch(kit.key())) += v[i].get();
            }

            d = filter(d);

            typename dcT::accessor acc;
            MADNESS_ASSERT(coeffs.find(acc, key));

            if (acc->second.has_coeff()) {
                print(" stuff in compress_op");
                const tensorT& c = acc->second.coeff().full_tensor_copy();
                if (c.dim(0) == k) {
                    d(cdata.s0) += c;
                }
                else {
                    d += c;
                }
            }

            // this is deep copy
//            coeffT s= copy(d(cdata.s0));
            tensorT s=copy(d(cdata.s0));

            if (key.level()> 0 && !nonstandard)
                d(cdata.s0) = 0.0;

//            d.reduceRank(thresh);
            coeffT dd=coeffT(d,targs);
            acc->second.set_coeff(dd);
//            s.reduceRank(thresh);
            coeffT ss=coeffT(s,targs);

            return s;
        }
#endif

        /// Changes non-standard compressed form to standard compressed form
        void standard(bool fence) {
            typename dcT::iterator end = coeffs.end();
            for (typename dcT::iterator it=coeffs.begin(); it!=end; ++it) {
                const keyT& key = it->first;
                nodeT& node = it->second;
                if (key.level()> 0 && node.has_coeff()) {
                    if (node.has_children()) {
                        // Zero out scaling coeffs
                    	node.coeff()(cdata.s0)=0.0;
                    	node.reduceRank(thresh);
                    }
                    else {
                        // Deleting both scaling and wavelet coeffs
                        node.clear_coeff();
                    }
                }
            }
            nonstandard = false;
            if (fence)
                world.gop.fence();
        }

        struct do_op_args {
            keyT key, d, dest;
            double tol, fac, cnorm;
            do_op_args() {}
            do_op_args(const keyT& key, const keyT& d, const keyT& dest, double tol, double fac, double cnorm)
                    : key(key), d(d), dest(dest), tol(tol), fac(fac), cnorm(cnorm) {}
            template <class Archive>
            void serialize(Archive& ar) {
                ar & archive::wrap_opaque(this,1);
            }
        };

        template <typename opT, typename R>
        Void do_apply_kernel(const opT* op, const Tensor<R>& c, const do_op_args& args) {
            tensorT result = op->apply(args.key, args.d, c, args.tol/args.fac/args.cnorm);

            //print("APPLY", key, d, opnorm, cnorm, result.normf());

            // Screen here to reduce communication cost of negligible data
            // and also to ensure we don't needlessly widen the tree when
            // applying the operator
            if (result.normf()> 0.3*args.tol/args.fac) {
                // OPTIMIZATION NEEDED HERE ... CHANGING THIS TO TASK NOT SEND REMOVED
                // BUILTIN OPTIMIZATION TO SHORTCIRCUIT MSG IF DATA IS LOCAL
                coeffs.task(args.dest, &nodeT::accumulate2, result, coeffs, args.dest, TaskAttributes::hipri());
            }
            return None;
        }


        /// same as do_apply_kernel, but use low rank tensors
        template <typename opT, typename R>
        Void do_apply_kernel2(const opT* op, const Tensor<R>& c, const do_op_args& args,
                const TensorArgs& apply_targs) {

            tensorT result_full = op->apply(args.key, args.d, c, args.tol/args.fac/args.cnorm);
            coeffT result=coeffT(result_full,apply_targs);

            // Screen here to reduce communication cost of negligible data
            // and also to ensure we don't needlessly widen the tree when
            // applying the operator
            double norm;
            MADNESS_ASSERT(result.tensor_type()==TT_FULL or result.tensor_type()==TT_2D);
            if (result.tensor_type()==TT_2D) norm=result.config().svd_normf();
            if (result.tensor_type()==TT_FULL) norm=result.normf();
            if (norm > 0.3*args.tol/args.fac) {
                // OPTIMIZATION NEEDED HERE ... CHANGING THIS TO TASK NOT SEND REMOVED
                // BUILTIN OPTIMIZATION TO SHORTCIRCUIT MSG IF DATA IS LOCAL
                coeffs.task(args.dest, &nodeT::accumulate, result, coeffs, args.dest, apply_targs,
                        TaskAttributes::hipri());
            }
            return None;
        }

        /// apply an operator on the coeffs c (at node key)

        /// the result is accumulated inplace to this's tree at various FunctionNodes
        /// @param[in] op	the operator to act on the source function
        /// @param[in] f	the source function (not used???)
        /// @param[in] key	key of the source FunctionNode of f which is processed
        /// @param[in] c	coeffs of the FunctionNode of f which is processed
        template <typename opT, typename R>
        Void do_apply(const opT* op, const FunctionImpl<R,NDIM>* f, const keyT& key, const Tensor<R>& c) {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            // insert timer here
            double fac = 10.0; //3.0; // 10.0 seems good for qmprop ... 3.0 OK for others
            double cnorm = c.normf();
            //const long lmax = 1L << (key.level()-1);

            const std::vector<keyT>& disp = op->get_disp(key.level());

            static const std::vector<bool> is_periodic(NDIM,false); // Periodic sum is already done when making rnlp

            for (typename std::vector<keyT>::const_iterator it=disp.begin(); it != disp.end(); ++it) {
                const keyT& d = *it;

                keyT dest = neighbor(key, d, is_periodic);

                if (dest.is_valid()) {
                    double opnorm = op->norm(key.level(), d);
                    // working assumption here is that the operator is isotropic and
                    // montonically decreasing with distance
                    double tol = truncate_tol(thresh, key);

                    //print("APP", key, dest, cnorm, opnorm, (cnorm*opnorm> tol/fac));

                    if (cnorm*opnorm> tol/fac) {

                        // Most expensive part is the kernel ... do it in a separate task
                        if (d.distsq()==0) {
                            // This introduces finer grain parallelism
                            ProcessID where = world.rank();
                            do_op_args args(key, d, dest, tol, fac, cnorm);
                            woT::task(where, &implT:: template do_apply_kernel<opT,R>, op, c, args);
                        } else {
                            tensorT result = op->apply(key, d, c, tol/fac/cnorm);
                            if (result.normf()> 0.3*tol/fac) {
                                coeffs.task(dest, &nodeT::accumulate2, result, coeffs, dest, TaskAttributes::hipri());
                            }
                        }
                    } else if (d.distsq() >= 1)
                        break; // Assumes monotonic decay beyond nearest neighbor
                }
            }
            return None;
        }

        /// apply an operator on f to return this
        template <typename opT, typename R>
        void apply(opT& op, const FunctionImpl<R,NDIM>& f, const std::vector<bool>& is_periodic, bool fence) {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            typename dcT::const_iterator end = f.coeffs.end();
            for (typename dcT::const_iterator it=f.coeffs.begin(); it!=end; ++it) {
                // looping through all the coefficients in the source
                const keyT& key = it->first;
                const FunctionNode<R,NDIM>& node = it->second;
                if (node.has_coeff()) {
                    if (node.coeff().dim(0) != k || op.doleaves) {
                        ProcessID p = FunctionDefaults<NDIM>::get_apply_randomize() ? world.random_proc() : coeffs.owner(key);
                        woT::task(p, &implT:: template do_apply<opT,R>, &op, &f, key, node.full_tensor_copy());
                    }
                }
            }
            if (fence)
                world.gop.fence();
        }


        /// apply an operator on the coeffs c (at node key)

        /// the result is accumulated inplace to this's tree at various FunctionNodes
        /// @param[in] op     the operator to act on the source function
        /// @param[in] key    key of the source FunctionNode of f which is processed
        /// @param[in] coeff  coeffs of FunctionNode being processed
        template <typename opT, typename R>
        Void do_apply_source_driven(const opT* op, const keyT& key, const coeffT& coeff) {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            // insert timer here

            // fac is the number of contributing neighbors (approx)
            const double tol = truncate_tol(thresh, key);

            double fac = 10.0; //3.0; // 10.0 seems good for qmprop ... 3.0 OK for others
            if (NDIM==6) fac=729; //100.0;
            double cnorm = coeff.normf();

            double wall0=wall_time();
            bool verbose=false;
            long neighbors=0;
            long generated_terms=0;


            // for accumulation: keep slightly tighter TensorArgs
            TensorArgs apply_targs(targs);
            apply_targs.thresh=tol/fac;

            // for the kernel it may be more efficient to do the convolution in full rank
            const long break_even=100;       // what should this number be ?!
            tensorT coeff_full;
            if (coeff.rank()>break_even) coeff_full=coeff.full_tensor_copy();
            if (coeff.tensor_type()==TT_FULL) coeff_full=coeff.full_tensor();


            const std::vector<keyT>& disp = op->get_disp(key.level());

            static const std::vector<bool> is_periodic(NDIM,false); // Periodic sum is already done when making rnlp

            for (typename std::vector<keyT>::const_iterator it=disp.begin(); it != disp.end(); ++it) {
                const keyT& d = *it;

                keyT dest = neighbor(key, d, is_periodic);

                if (dest.is_valid()) {
                    double opnorm = op->norm(key.level(), d);
                    // working assumption here is that the operator is isotropic and
                    // montonically decreasing with distance

                    //print("APP", key, dest, cnorm, opnorm, (cnorm*opnorm> tol/fac));

                    if (cnorm*opnorm> tol/fac) {
//                        double wall00=wall_time();
                        neighbors++;

                        double cost_ratio=op->estimate_costs(key, d, coeff, tol/fac/cnorm, tol/fac);

                        if (cost_ratio>0.0) {
                            // Most expensive part is the kernel ... do it in a separate task
                            if (coeff_full.has_data() and (cost_ratio<1.0)) {


                                // This introduces finer grain parallelism
                                ProcessID here = world.rank();
//                                ProcessID there =  world.random_proc();
//            			        ProcessID where = FunctionDefaults<NDIM>::get_apply_randomize() ? there : here;
                                do_op_args args(key, d, dest, tol, fac, cnorm);
                                woT::task(here, &implT:: template do_apply_kernel2<opT,R>, op, coeff_full,
                                        args,apply_targs);

                            } else {

                                // apply2 returns result in SVD form
                                coeffT result = op->apply2(key, d, coeff, tol/fac/cnorm, tol/fac);
                                double result_norm=-1.0;
                                if (result.tensor_type()==TT_2D) result_norm=result.config().svd_normf();
                                if (result.tensor_type()==TT_FULL) result_norm=result.normf();
                                MADNESS_ASSERT(result_norm>-0.5);
                                if (result_norm> 0.3*tol/fac) {

                                    generated_terms+=result.rank();

                                    // accumulate also expects result in SVD form
                                    coeffs.task(dest, &nodeT::accumulate, result, coeffs, dest, apply_targs,
                                            TaskAttributes::hipri());
                                }
                            }
                        }
//                        double wall11=wall_time();
//                        print("finished dest",dest,wall11-wall00);
                    } else if (d.distsq() >= 1)
                        break; // Assumes monotonic decay beyond nearest neighbor
                }
            }
            double wall1=wall_time();
            if (verbose) {
                print("done with source node",key,wall1-wall0, cnorm, neighbors,generated_terms,coeff.rank(),
                        coeff_full.has_data());
            }
            return None;
        }



        /// similar to apply, but for low rank coeffs
         template <typename opT, typename R>
         void apply_source_driven(opT& op, const FunctionImpl<R,NDIM>& f,
                 const std::vector<bool>& is_periodic, bool fence) {
             PROFILE_MEMBER_FUNC(FunctionImpl);

//             Vector<Translation, NDIM> l;
//             l[0]=4; l[1]=4; l[2]=4; l[3]=4; l[4]=4; l[5]=4;      // k=5, rank=995
//             l[0]=1; l[1]=1; l[2]=1; l[3]=1; l[4]=1; l[5]=1;      // k=5, rank=507
//             l[0]=2; l[1]=2; l[2]=2; l[3]=2; l[4]=2; l[5]=2;        // k=5, rank=507
//             l[0]=2; l[1]=1; l[2]=2; l[3]=1; l[4]=2; l[5]=2;        // k=5, rank=58
//             l[0]=1; l[1]=2; l[2]=2; l[3]=1; l[4]=3; l[5]=2;        // k=5, rank=5
//             const keyT key1(3,l);


             // looping through all the coefficients of the source f
             typename dcT::const_iterator end = f.get_coeffs().end();
             for (typename dcT::const_iterator it=f.get_coeffs().begin(); it!=end; ++it) {

                 const keyT& key = it->first;
                 const coeffT& coeff = it->second.coeff();

//                 if (key==key1) {

                     if (coeff.has_data() and (coeff.rank()!=0)) {
//                         ProcessID p = f.get_coeffs().owner(key);
                         ProcessID p = FunctionDefaults<NDIM>::get_apply_randomize() ? world.random_proc() : coeffs.owner(key);

                         woT::task(p, &implT:: template do_apply_source_driven<opT,R>, &op, key, coeff);
//                     } else {
//                         print("done with empty source node",key);
                     }
//                 }
             }

             if (fence)
                 world.gop.fence();

             // reduce the rank of the final nodes
             flo_unary_op_node_inplace(do_reduce_rank(targs),true);

             if (fence)
                 world.gop.fence();
         }

        /// compute the error for a (compressed) node using the wavelet coeffs

        /// @param[in]	node			the node whose error to be computed
        /// @return 	error			the error (zero if leaf node)
        double compute_error(const nodeT& node) const {

        	if (node.is_leaf()) return 0.0;

        	MADNESS_ASSERT(node.has_coeff());
        	MADNESS_ASSERT(node.coeff().dim(0)==2*this->get_k());

    		coeffT d = copy(node.coeff());
       		d(cdata.s0)=0.0;
        	return d.normf();
        }



        /// Returns the square of the error norm in the box labelled by key

        /// Assumed to be invoked locally but it would be easy to eliminate
        /// this assumption
        template <typename opT>
        double err_box(const keyT& key, const nodeT& node, const opT& func,
                       int npt, const Tensor<double>& qx, const Tensor<double>& quad_phit,
                       const Tensor<double>& quad_phiw) const {

            std::vector<long> vq(NDIM);
            for (std::size_t i=0; i<NDIM; ++i)
                vq[i] = npt;
            tensorT fval(vq,false), work(vq,false), result(vq,false);

            // Compute the "exact" function in this volume at npt points
            // where npt is usually this->npt+1.
            fcube(key, func, qx, fval);

            // Transform into the scaling function basis of order npt
            double scale = pow(0.5,0.5*NDIM*key.level())*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
            fval = fast_transform(fval,quad_phiw,result,work).scale(scale);

            // Subtract to get the error ... the original coeffs are in the order k
            // basis but we just computed the coeffs in the order npt(=k+1) basis
            // so we can either use slices or an iterator macro.
            const tensorT& coeff = node.full_tensor_copy();
            ITERATOR(coeff,fval(IND)-=coeff(IND););
            // flo note: we do want to keep a full tensor here!

            // Compute the norm of what remains
            double err = fval.normf();
            return err*err;
        }

        template <typename opT>
        class do_err_box {
            const implT* impl;
            const opT* func;
            int npt;
            Tensor<double> qx;
            Tensor<double> quad_phit;
            Tensor<double> quad_phiw;
        public:
            do_err_box() {}

            do_err_box(const implT* impl, const opT* func, int npt, const Tensor<double>& qx,
                       const Tensor<double>& quad_phit, const Tensor<double>& quad_phiw)
                    : impl(impl), func(func), npt(npt), qx(qx), quad_phit(quad_phit), quad_phiw(quad_phiw) {}

            do_err_box(const do_err_box& e)
                    : impl(e.impl), func(e.func), npt(e.npt), qx(e.qx), quad_phit(e.quad_phit), quad_phiw(e.quad_phiw) {}

            double operator()(typename dcT::const_iterator& it) const {
                const keyT& key = it->first;
                const nodeT& node = it->second;
                if (node.has_coeff())
                    return impl->err_box(key, node, *func, npt, qx, quad_phit, quad_phiw);
                else
                    return 0.0;
            }

            double operator()(double a, double b) const {
                return a+b;
            }

            template <typename Archive>
            void serialize(const Archive& ar) {
                throw "not yet";
            }
        };

        /// Returns the sum of squares of errors from local info ... no comms
        template <typename opT>
        double errsq_local(const opT& func) const {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            // Make quadrature rule of higher order
            const int npt = cdata.npt + 1;
            Tensor<double> qx, qw, quad_phi, quad_phiw, quad_phit;
            FunctionCommonData<T,NDIM>::_init_quadrature(k+1, npt, qx, qw, quad_phi, quad_phiw, quad_phit);

            typedef Range<typename dcT::const_iterator> rangeT;
            rangeT range(coeffs.begin(), coeffs.end());
            return world.taskq.reduce< double,rangeT,do_err_box<opT> >(range,
                    do_err_box<opT>(this, &func, npt, qx, quad_phit, quad_phiw));
        }

        /// Returns \c int(f(x),x) in local volume
        T trace_local() const;

        struct do_norm2sq_local {
            double operator()(typename dcT::const_iterator& it) const {
                const nodeT& node = it->second;
                if (node.has_coeff()) {
                    double norm = node.coeff().normf();
                    return norm*norm;
                }
                else {
                    return 0.0;
                }
            }

            double operator()(double a, double b) const {
                return (a+b);
            }

            template <typename Archive> void serialize(const Archive& ar) {
                throw "NOT IMPLEMENTED";
            }
        };

        /// Returns the square of the local norm ... no comms
        double norm2sq_local() const {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            typedef Range<typename dcT::const_iterator> rangeT;
            return world.taskq.reduce<double,rangeT,do_norm2sq_local>(rangeT(coeffs.begin(),coeffs.end()),
                    do_norm2sq_local());
        }

        /// Returns the inner product ASSUMING same distribution
        /// note that if (g==f) both functions might be reconstructed
        template <typename R>
        TENSOR_RESULT_TYPE(T,R) inner_local(const FunctionImpl<R,NDIM>& g) const {
            TENSOR_RESULT_TYPE(T,R) sum = 0.0;
            typename dcT::const_iterator end = coeffs.end();
            for (typename dcT::const_iterator it=coeffs.begin(); it!=end; ++it) {
                const nodeT& fnode = it->second;
                if (fnode.has_coeff()) {
                    if (g.coeffs.probe(it->first)) {
                        const FunctionNode<R,NDIM>& gnode = g.coeffs.find(it->first).get()->second;
                        if (gnode.has_coeff()) {
                            if (gnode.coeff().dim(0) != fnode.coeff().dim(0)) {
                                madness::print("INNER", it->first, gnode.coeff().dim(0),fnode.coeff().dim(0));
                                MADNESS_EXCEPTION("functions have different k or compress/reconstruct error", 0);
                            }
                            sum += fnode.coeff().trace_conj(gnode.coeff());
                        }
                    }
                }
            }
            return sum;
        }


        /// Returns the inner product ASSUMING same distribution
        /// note that if (g==f) both functions might be reconstructed
        /// same as inner_local, but this is an on-demand function
        /// (and therefore not const)
        template <typename R>
        TENSOR_RESULT_TYPE(T,R) inner_local2(const FunctionImpl<R,NDIM>& f) {
            MADNESS_ASSERT(not f.is_compressed());
            MADNESS_ASSERT(this->is_on_demand());

        	// looks a little weird: first make sure all the gnodes are constructed,
        	// afterwards use them to compute the inner product;
    		this->fill_on_demand_tree(&f,false);

            TENSOR_RESULT_TYPE(T,R) sum = 0.0;
            typename dcT::const_iterator end = f.coeffs.end();
            for (typename dcT::const_iterator it=f.coeffs.begin(); it!=end; ++it) {
            	const keyT key=it->first;
                const nodeT& fnode = it->second;
                if (fnode.has_coeff()) {

//                   	coeffT gcoeff=this->functor->coeff(key);
                	coeffT gcoeff=coeffs.find(key).get()->second.coeff();

                    // compute inner product
                    sum += fnode.coeff().trace_conj(gcoeff);
                }
            }
            return sum;
        }

        /// Returns the maximum local depth of the tree ... no communications.
        std::size_t max_local_depth() const {
            std::size_t maxdepth = 0;
            typename dcT::const_iterator end = coeffs.end();
            for (typename dcT::const_iterator it=coeffs.begin(); it!=end; ++it) {
                std::size_t N = (std::size_t) it->first.level();
                if (N> maxdepth)
                    maxdepth = N;
            }
            return maxdepth;
        }


        /// Returns the maximum depth of the tree ... collective ... global sum/broadcast
        std::size_t max_depth() const {
            std::size_t maxdepth  = max_local_depth();
            world.gop.max(maxdepth);
            return maxdepth;
        }

        /// Returns the max number of nodes on a processor
        std::size_t max_nodes() const {
            std::size_t maxsize = 0;
            maxsize = coeffs.size();
            world.gop.max(maxsize);
            return maxsize;
        }

        /// Returns the min number of nodes on a processor
        std::size_t min_nodes() const {
            std::size_t minsize = 0;
            minsize = coeffs.size();
            world.gop.min(minsize);
            return minsize;
        }

        /// Returns the size of the tree structure of the function ... collective global sum
        std::size_t tree_size() const {
            std::size_t sum = 0;
            sum = coeffs.size();
            world.gop.sum(sum);
            return sum;
        }

        /// Returns the number of coefficients in the function ... collective global sum
        std::size_t size() const {
            std::size_t sum = 0;
#if 1
            typename dcT::const_iterator end = coeffs.end();
            for (typename dcT::const_iterator it=coeffs.begin(); it!=end; ++it) {
                const nodeT& node = it->second;
                if (node.has_coeff())
                    sum+=node.size();
            }
//            print("proc",world.rank(),sum);
#else
            typename dcT::const_iterator end = coeffs.end();
            for (typename dcT::const_iterator it=coeffs.begin(); it!=end; ++it) {
                const nodeT& node = it->second;
                if (node.has_coeff())
                    ++sum;
            }
            if (is_compressed())
                for (std::size_t i=0; i<NDIM; ++i)
                    sum *= 2*cdata.k;
            else
                for (std::size_t i=0; i<NDIM; ++i)
                    sum *= cdata.k;
#endif
            world.gop.sum(sum);

            return sum;
        }


        /// print the number of configurations per node
        void print_stats() const {
            if (this->targs.tt==TT_FULL) return;
            int dim=NDIM/2;
            int k0=k;
            if (is_compressed()) k0=2*k;
            Tensor<long> n(int(std::pow(k0,dim)+1));
            if (world.rank()==0) print("n.size(),k0,dim",n.size(),k0,dim);
            typename dcT::const_iterator end = coeffs.end();
            for (typename dcT::const_iterator it=coeffs.begin(); it!=end; ++it) {
                const nodeT& node = it->second;
                if (node.has_coeff()) {
                    if (node.coeff().rank()>long(n.size())) {
                        print("large rank",node.coeff().rank());
                    } else if (node.coeff().rank()<0) {
                        print("small rank",node.coeff().rank());
                    } else {
                        n[node.coeff().rank()]++;
                    }
                }
            }

            world.gop.sum(n.ptr(), n.size());

            if (world.rank()==0) {
                print("configurations     number of nodes");
                for (unsigned int i=0; i<n.size(); i++) {
                    long m=n[i];
                    if (world.rank()==0) print("           ",i,"    ",m);
                }
            }
        }

        /// In-place scale by a constant
        void scale_inplace(const T q, bool fence);

        /// Out-of-place scale by a constant
        template <typename Q, typename F>
        void scale_oop(const Q q, const FunctionImpl<F,NDIM>& f, bool fence) {
            typedef typename FunctionImpl<F,NDIM>::nodeT fnodeT;
            typedef typename FunctionImpl<F,NDIM>::dcT fdcT;
            typename fdcT::const_iterator end = f.coeffs.end();
            for (typename fdcT::const_iterator it=f.coeffs.begin(); it!=end; ++it) {
                const keyT& key = it->first;
                const fnodeT& node = it->second;

                if (node.has_coeff()) {
                    coeffs.replace(key,nodeT(coeffT(node.full_tensor_copy()*q,targs),node.has_children()));
//                    coeffs.replace(key,nodeT(node.coeff()*q,node.has_children()));
                }
                else {
//                    coeffs.replace(key,nodeT(tensorT(),node.has_children()).node_to_low_rank());
//                    coeffs.replace(key,nodeT(coeffT(FunctionDefaults<NDIM>::get_tensor_type()),node.has_children()));
                    coeffs.replace(key,nodeT(coeffT(),node.has_children()));
                }
            }
            if (fence)
                world.gop.fence();
        }

    private:
        /// Assignment is not allowed ... not even possible now that we have reference members
        //FunctionImpl<T>& operator=(const FunctionImpl<T>& other);
    };

    namespace archive {
        template <class Archive, class T, std::size_t NDIM>
        struct ArchiveLoadImpl<Archive,const FunctionImpl<T,NDIM>*> {
            static void load(const Archive& ar, const FunctionImpl<T,NDIM>*& ptr) {
                uniqueidT id;
                ar & id;
                World* world = World::world_from_id(id.get_world_id());
                MADNESS_ASSERT(world);
                ptr = static_cast< const FunctionImpl<T,NDIM>*>(world->ptr_from_id< WorldObject< FunctionImpl<T,NDIM> > >(id));
                if (!ptr)
                    MADNESS_EXCEPTION("FunctionImpl: remote operation attempting to use a locally uninitialized object",0);
            }
        };

        template <class Archive, class T, std::size_t NDIM>
        struct ArchiveStoreImpl<Archive,const FunctionImpl<T,NDIM>*> {
            static void store(const Archive& ar, const FunctionImpl<T,NDIM>*const& ptr) {
                ar & ptr->id();
            }
        };

        template <class Archive, class T, std::size_t NDIM>
        struct ArchiveLoadImpl<Archive, FunctionImpl<T,NDIM>*> {
            static void load(const Archive& ar, FunctionImpl<T,NDIM>*& ptr) {
                uniqueidT id;
                ar & id;
                World* world = World::world_from_id(id.get_world_id());
                MADNESS_ASSERT(world);
                ptr = static_cast< FunctionImpl<T,NDIM>*>(world->ptr_from_id< WorldObject< FunctionImpl<T,NDIM> > >(id));
                if (!ptr)
                    MADNESS_EXCEPTION("FunctionImpl: remote operation attempting to use a locally uninitialized object",0);
            }
        };

        template <class Archive, class T, std::size_t NDIM>
        struct ArchiveStoreImpl<Archive, FunctionImpl<T,NDIM>*> {
            static void store(const Archive& ar, FunctionImpl<T,NDIM>*const& ptr) {
                ar & ptr->id();
            }
        };

        template <class Archive, class T, std::size_t NDIM>
        struct ArchiveLoadImpl<Archive, std::shared_ptr<const FunctionImpl<T,NDIM> > > {
            static void load(const Archive& ar, std::shared_ptr<const FunctionImpl<T,NDIM> >& ptr) {
                const FunctionImpl<T,NDIM>* f = NULL;
                ArchiveLoadImpl<Archive, const FunctionImpl<T,NDIM>*>::load(ar, f);
                ptr.reset(f, & detail::no_delete<const FunctionImpl<T,NDIM> >);
            }
        };

        template <class Archive, class T, std::size_t NDIM>
        struct ArchiveStoreImpl<Archive, std::shared_ptr<const FunctionImpl<T,NDIM> > > {
            static void store(const Archive& ar, const std::shared_ptr<const FunctionImpl<T,NDIM> >& ptr) {
                ArchiveStoreImpl<Archive, const FunctionImpl<T,NDIM>*>::store(ar, ptr.get());
            }
        };

        template <class Archive, class T, std::size_t NDIM>
        struct ArchiveLoadImpl<Archive, std::shared_ptr<FunctionImpl<T,NDIM> > > {
            static void load(const Archive& ar, std::shared_ptr<FunctionImpl<T,NDIM> >& ptr) {
                FunctionImpl<T,NDIM>* f = NULL;
                ArchiveLoadImpl<Archive, FunctionImpl<T,NDIM>*>::load(ar, f);
                ptr.reset(f, & detail::no_delete<FunctionImpl<T,NDIM> >);
            }
        };

        template <class Archive, class T, std::size_t NDIM>
        struct ArchiveStoreImpl<Archive, std::shared_ptr<FunctionImpl<T,NDIM> > > {
            static void store(const Archive& ar, const std::shared_ptr<FunctionImpl<T,NDIM> >& ptr) {
                ArchiveStoreImpl<Archive, FunctionImpl<T,NDIM>*>::store(ar, ptr.get());
            }
        };
    }

}

#endif // MADNESS_MRA_FUNCIMPL_H__INCLUDED
