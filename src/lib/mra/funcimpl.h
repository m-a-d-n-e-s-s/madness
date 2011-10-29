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
#include <world/scopedptr.h>
#include <misc/misc.h>
#include <tensor/tensor.h>
#include <mra/key.h>
#include <mra/funcdefaults.h>
#include "mra/function_factory_and_interface.h"
#include "mra/gentensor.h"
#include <world/typestuff.h>
#include "mra/function_common_data.h"

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

    /// A process map that uses only part of the key for determining the owner

    /// @tparam     NDIM dimension of the key
    /// @tparam     LDIM use only dim 0..LDIM-1 of NDIM to compute the owner
    template<size_t NDIM, size_t LDIM>
    class PartialKeyMap : public WorldDCPmapInterface<Key<NDIM> > {
    private:
        const int nproc;
        const ProcessID me;
        const int n;

    public:
        PartialKeyMap(World& world, int n = 4) :
                nproc(world.nproc()), me(world.rank()), n(n) {
        }

        ProcessID
        owner(const Key<NDIM>& key) const {
            if (key.level() == 0) {
                return 0;
            }
            else if (key.level() <= n) {
                Key<LDIM> key1;
                Key<NDIM-LDIM> key2;
                key.break_apart(key1,key2);
//                return hash(key1) % nproc;
                return key1.hash() % nproc;
            }
            else {
                return key.parent(key.level() - n).hash() % nproc;
            }
        }
    };



#if 0
    // moved to its own file

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

        Tensor<double> h0, h1, g0, g1;      ///< The separate blocks of twoscale coefficients
        Tensor<double> h0T, h1T, g0T, g1T;  ///< The separate blocks of twoscale coefficients
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
#endif


    /// returns true if the function has a leaf node at key (works only locally)
    template<typename T, std::size_t NDIM>
    struct leaf_op {
        typedef FunctionImpl<T,NDIM> implT;
        const implT* f;

        leaf_op() {}
        leaf_op(const implT* f) : f(f) {}

        /// pre/post-determination is the same here
        bool operator()(const Key<NDIM>& key, const GenTensor<T>& coeff=GenTensor<T>()) const {
            MADNESS_ASSERT(f->get_coeffs().is_local(key));
            return (not f->get_coeffs().find(key).get()->second.has_children());
        }

        template <typename Archive> void serialize (Archive& ar) {
            ar & f;
        }
    };


    /// returns true if the node is well represented compared to its parent
    template<typename T, std::size_t NDIM>
    struct error_leaf_op {
        typedef FunctionImpl<T,NDIM> implT;
        typedef GenTensor<T> coeffT;
        const implT* f;

        error_leaf_op() {}
        error_leaf_op(const implT* f) : f(f) {}

        /// post-determination

        /// @param[in]  key the FunctionNode which we want to determine if it's a leaf node
        /// @param[in]  coeff   the coeffs of key
        /// @param[in]  parent  the coeffs of key's parent node
        /// @return is the FunctionNode of key a leaf node?
        bool operator()(const Key<NDIM>& key, const coeffT& coeff, const coeffT& parent) const {
            if (parent.has_no_data()) return false;
            coeffT upsampled=f->upsample(key,parent);
            upsampled.scale(-1.0);
            upsampled+=coeff;
            const double dnorm=upsampled.normf();
            const bool is_leaf=(dnorm<f->truncate_tol(f->get_thresh(),key.level()));
            return is_leaf;
        }

        template <typename Archive> void serialize (Archive& ar) {ar & f;}
    };



//    /// returns true if the result of a multiplication is a leaf node
//    template<typename T, std::size_t NDIM>
//    struct mul_leaf_op {
//        typedef FunctionImpl<T,NDIM> implT;
//
//        const implT* f;
//
//        mul_leaf_op() {}
//        mul_leaf_op(const implT* f) : f(f) {}
//
//        /// return true if f is a leaf and the result is well-represented
//        bool operator()(const Key<NDIM>& key, const GenTensor<T>& fcoeff, const GenTensor<T>& gcoeff) const {
//            double flo,fhi,glo,ghi;
//            bool is_leaf=true;
//            f->tnorm(fcoeff,&flo,&fhi);
//            f->tnorm(gcoeff,&glo,&ghi);
//            double total_hi=glo*fhi + ghi*flo + fhi*ghi;
//            if (total_hi>f->truncate_tol(f->get_thresh(),key)) is_leaf=false;
//            return is_leaf;
//        }
//        template <typename Archive> void serialize (Archive& ar) {
//            ar & f;
//        }
//    };



    /// returns true if the result of a hartree_product is a leaf node (compute norm & error)
    template<typename T, size_t NDIM>
    struct hartree_leaf_op {

        typedef FunctionImpl<T,NDIM> implT;
        const FunctionImpl<T,NDIM>* f;
        long k;

        hartree_leaf_op() {}
        hartree_leaf_op(const implT* f, const long& k) : f(f), k(k) {}

        /// no pre-determination
        bool operator()(const Key<NDIM>& key) const {return false;}

        /// no post-determination
        bool operator()(const Key<NDIM>& key, const GenTensor<T>& coeff) const {
            MADNESS_EXCEPTION("no post-determination in hartree_leaf_op",1);
            return true;
        }

        /// post-determination:
        /// return true if f is a leaf and the result is well-represented

        /// @param[in]  key the hi-dimensional key (breaks into keys for f and g)
        /// @param[in]  fcoeff coefficients of f of its appropriate key in NS form
        /// @param[in]  gcoeff coefficients of g of its appropriate key in NS form
        bool operator()(const Key<NDIM>& key, const GenTensor<T>& fcoeff, const GenTensor<T>& gcoeff) const {
//        bool operator()(const Key<NDIM>& key, const GenTensor<T>& coeff) const {

            // we expect coeffs that are a direct product of the two orbital coefficients
//            MADNESS_ASSERT(coeff.rank()==1);
//            const GenTensor<T> fcoeff=GenTensor<T>(coeff.config().ref_vector(0).reshape(k,k,k),0.0,TT_FULL);
//            const GenTensor<T> gcoeff=GenTensor<T>(coeff.config().ref_vector(1).reshape(k,k,k),0.0,TT_FULL);
            Slice s = Slice(0,k-1);
            std::vector<Slice> s0(NDIM/2,s);

            const double tol=f->get_thresh();
            const double thresh=f->truncate_tol(tol, key);
            // include the wavelets in the norm, makes it much more accurate
            const double fnorm=fcoeff.normf();
            const double gnorm=gcoeff.normf();

            // norm of the scaling function coefficients
            const GenTensor<T> s1=fcoeff(s0);
            const GenTensor<T> s2=gcoeff(s0);
            const double sfnorm=s1.normf();
            const double sgnorm=s2.normf();

            // if the final norm is small, perform the hartree product and return
            const double norm=fnorm*gnorm;  // computing the outer product
            if (norm < thresh) return true;

            // get the error of both functions and of the pair function
            const double ferror=sqrt(fnorm*fnorm-sfnorm*sfnorm);
            const double gerror=sqrt(gnorm*gnorm-sgnorm*sgnorm);

            // if the expected error is small, perform the hartree product and return
            const double error=fnorm*gerror + ferror*gnorm + ferror*gerror;
            if (error < thresh) return true;
            return false;
        }
        template <typename Archive> void serialize (Archive& ar) {
            ar & f & k;
        }
    };

    /// returns true if the result of the convolution operator op with some provided
    /// coefficients will be small
    template<typename T, size_t NDIM, typename opT>
    struct op_leaf_op {
        typedef FunctionImpl<T,NDIM> implT;

        const opT* op;    ///< the convolution operator
        const implT* f;   ///< the source or result function, needed for truncate_tol
        op_leaf_op() {}
        op_leaf_op(const opT* op, const implT* f) : op(op), f(f) {}

        /// pre-determination: we can't know if this will be a leaf node before we got the final coeffs
        bool operator()(const Key<NDIM>& key) const {return true;}

        /// post-determination: return true if operator and coefficient norms are small
        bool operator()(const Key<NDIM>& key, const GenTensor<T>& coeff) const {
            if (key.level()<2) return false;
            const double cnorm=coeff.normf();
            const double thresh=f->truncate_tol(f->get_thresh(),key);
            const std::vector<Key<NDIM> >& disp = op->get_disp(key.level());
            const Key<NDIM>& d = *disp.begin();         // use the zero-displacement for screening
            const double opnorm = op->norm(key.level(), d, key);
            const double norm=opnorm*cnorm;
            return norm<thresh;

        }

        template <typename Archive> void serialize (Archive& ar) {
            ar & op & f;
        }

    };


    /// returns true if the result of a hartree_product is a leaf node
    /// criteria are error, norm and its effect on a convolution operator
    template<typename T, size_t NDIM, size_t LDIM, typename opT>
    struct hartree_convolute_leaf_op {

        typedef FunctionImpl<T,NDIM> implT;
        typedef FunctionImpl<T,LDIM> implL;

        const FunctionImpl<T,NDIM>* f;
        const implL* g;     // for use of its cdata only
        const opT* op;

        hartree_convolute_leaf_op() {}
        hartree_convolute_leaf_op(const implT* f, const implL* g, const opT* op)
            : f(f), g(g), op(op) {}

        /// no pre-determination
        bool operator()(const Key<NDIM>& key) const {return true;}

        /// no post-determination
        bool operator()(const Key<NDIM>& key, const GenTensor<T>& coeff) const {
            MADNESS_EXCEPTION("no post-determination in hartree_convolute_leaf_op",1);
            return true;
        }

        /// post-determination:
        /// return true if f is a leaf and the result is well-represented

        /// @param[in]  key the hi-dimensional key (breaks into keys for f and g)
        /// @param[in]  fcoeff coefficients of f of its appropriate key in NS form
        /// @param[in]  gcoeff coefficients of g of its appropriate key in NS form
        bool operator()(const Key<NDIM>& key, const GenTensor<T>& fcoeff, const GenTensor<T>& gcoeff) const {
//        bool operator()(const Key<NDIM>& key, const GenTensor<T>& coeff) const {

            if (key.level()<2) return false;
            // we expect coeffs that are a direct product of the two orbital coefficients
//            MADNESS_ASSERT(coeff.rank()==1);
//            const long k=f->get_k();
//            const GenTensor<T> fcoeff=GenTensor<T>(coeff.config().ref_vector(0).reshape(2*k,2*k,2*k),0.0,TT_FULL);
//            const GenTensor<T> gcoeff=GenTensor<T>(coeff.config().ref_vector(1).reshape(2*k,2*k,2*k),0.0,TT_FULL);

            const double tol=f->get_thresh();
            const double thresh=f->truncate_tol(tol, key);
            // include the wavelets in the norm, makes it much more accurate
            const double fnorm=fcoeff.normf();
            const double gnorm=gcoeff.normf();

            // norm of the scaling function coefficients
            const GenTensor<T> s1=fcoeff(g->get_cdata().s0);
            const GenTensor<T> s2=gcoeff(g->get_cdata().s0);
            const double sfnorm=s1.normf();
            const double sgnorm=s2.normf();

            // if the final norm is small, perform the hartree product and return
            const double norm=fnorm*gnorm;  // computing the outer product
            if (norm < thresh) return true;

            // get the error of both functions and of the pair function
            const double ferror=sqrt(fnorm*fnorm-sfnorm*sfnorm);
            const double gerror=sqrt(gnorm*gnorm-sgnorm*sgnorm);

            // if the expected error is small, perform the hartree product and return
            const double error=fnorm*gerror + ferror*gnorm + ferror*gerror;
            if (error < thresh) return true;

            // now check if the norm of this and the norm of the operator are significant
            const std::vector<Key<NDIM> >& disp = op->get_disp(key.level());
            const Key<NDIM>& d = *disp.begin();         // use the zero-displacement for screening
            const double opnorm = op->norm(key.level(), d, key);
            const double final_norm=opnorm*sfnorm*sgnorm;
            if (final_norm < thresh) return true;

            return false;
        }
        template <typename Archive> void serialize (Archive& ar) {
            ar & f & op;
        }
    };

    template<typename T, size_t NDIM>
    struct noop {

        bool operator()(const Key<NDIM>& key, const GenTensor<T>& fcoeff, const GenTensor<T>& gcoeff) const {
            MADNESS_EXCEPTION("in noop::operator()",1);
            return true;
        }
        template <typename Archive> void serialize (Archive& ar) {}

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

        /// reduces the rank of the coefficients (if applicable)
        void reduceRank(const double& eps) {
        	_coeffs.reduceRank(eps);
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
        double accumulate(const coeffT& t, const typename FunctionNode<T,NDIM>::dcT& c,
        		const Key<NDIM>& key, const TensorArgs& args) {
            double cpu0=cpu_time();
            if (has_coeff()) {

#if 1
                // always do low rank
                coeff().add_SVD(t,args.thresh);


#else
                // make it tensor type depending on what we already have
                if (coeff().tensor_type()==TT_FULL) {
                    if (t.tensor_type()==TT_FULL)  coeff().full_tensor()+=t.full_tensor();
                    else coeff().full_tensor()+=t.full_tensor_copy();
                } else {
                    if (t.tensor_type()==TT_FULL) {
                        tensorT c=coeff().full_tensor_copy()+t.full_tensor();
                        coeff()=coeffT(c,TensorArgs(0.0,TT_FULL));
                    } else {
                        coeff().add_SVD(t,args.thresh);
                    }
                }
#endif

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
            double cpu1=cpu_time();
            return cpu1-cpu0;
        }

        T trace_conj(const FunctionNode<T,NDIM>& rhs) const {
        	return this->_coeffs.trace_conj((rhs._coeffs));
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar & coeff() & _has_children & _norm_tree;
        }

    };

    template <typename T, std::size_t NDIM>
    std::ostream& operator<<(std::ostream& s, const FunctionNode<T,NDIM>& node) {
        s << "(has_coeff=" << node.has_coeff() << ", has_children=" << node.has_children() << ", norm=";
        double norm = node.has_coeff() ? node.coeff().normf() : 0.0;
        if (norm < 1e-12)
            norm = 0.0;
//        norm=node.get_norm_tree();
        s << norm << ", rank="<< node.coeff().rank()<<")";
        return s;
    }

    /// shallow-copy, pared-down version of FunctionNode, for special purpose only
    template<typename T, std::size_t NDIM>
    struct ShallowNode {
        typedef GenTensor<T> coeffT;
        coeffT _coeffs;
        bool _has_children;
        ShallowNode() : _coeffs(), _has_children(false) {}
        ShallowNode(const FunctionNode<T,NDIM>& node)
               : _coeffs(node.coeff()), _has_children(node.has_children()) {}
        ShallowNode(const ShallowNode<T,NDIM>& node)
               : _coeffs(node.coeff()), _has_children(node._has_children) {}

        const coeffT& coeff() const {return _coeffs;}
        coeffT& coeff() {return _coeffs;}
        bool has_children() const {return _has_children;}
        bool is_leaf() const {return not _has_children;}
        template <typename Archive>
        void serialize(Archive& ar) {
            ar & coeff() & _has_children;
        }
    };

    template<typename T, size_t NDIM>
    struct impl_and_arg {
        typedef std::pair<Key<NDIM>, ShallowNode<T,NDIM> > datumT;
        typedef FunctionImpl<T,NDIM> implT;
        const implT* impl;
        datumT datum;

        impl_and_arg() : impl() {}
        impl_and_arg(const implT* impl, const datumT& datum)
            : impl(impl), datum(datum) {}
        impl_and_arg(const impl_and_arg<T,NDIM>& other) : impl(other.impl), datum(other.datum) {}

        Future<impl_and_arg> make_child(const Key<NDIM>& child) const {
            if (not impl) return Future<impl_and_arg>(impl_and_arg());
            if (datum.second.is_leaf()) return Future<impl_and_arg>(*this);
//            Future<datumT> datum1=impl->task(impl->get_coeffs().owner(child),
//                    &FunctionImpl<T,NDIM>::outermost_child,child,datum,TaskAttributes::hipri());
            Future<datumT> datum1=impl->task(impl->get_coeffs().owner(child), &implT::find_datum,child,
                    TaskAttributes::hipri());
            // wait for the nodes to arrive, and construct a new operator
            return impl->world.taskq.add(*const_cast<impl_and_arg*> (this), &impl_and_arg::make_child2,
                    impl,datum1);

        }
        impl_and_arg make_child2(const implT* impl1, const datumT& datum1) const {
            return impl_and_arg(impl1,datum1);
        }

        template <typename Archive> void serialize(const Archive& ar) {
            ar & impl & datum;
        }
    };



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
        Timer timer_accumulate;
        Timer timer_filter;
        Timer timer_compress_svd;
        Timer timer_target_driven;

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

        /// perform: this= alpha*f + beta*g, invoked by result

        /// f and g are reconstructed, so we can save on the compress operation,
        /// just walk down the joint tree, and add leaf coefficients
        /// @param[in]  alpha   prefactor for f
        /// @param[in]  f       first addend
        /// @param[in]  beta    prefactor for g
        /// @param[in]  g       second addend
        /// @return     nothing, but leaves this's tree reconstructed and as sum of f and g
        void gaxpy_oop_reconstructed(const double alpha, const implT& f,
                const double beta, const implT& g, const bool fence) {

            MADNESS_ASSERT(not f.is_compressed());
            MADNESS_ASSERT(not g.is_compressed());

            typedef std::pair<Key<NDIM>, ShallowNode<T,NDIM> > datumT;
            typedef FunctionImpl<T,NDIM> implL;

            ProcessID owner = coeffs.owner(cdata.key0);
            if (world.rank() == owner) {

                Future<datumT> dat1=f.task(owner, &implT::find_datum,f.cdata.key0,TaskAttributes::hipri());
                Future<datumT> dat2=g.task(owner, &implT::find_datum,g.cdata.key0,TaskAttributes::hipri());

                // have to wait for this, but should be local..
                datumT fdatum=dat1.get();
                datumT gdatum=dat2.get();

                add_op op(&f,&g,fdatum,gdatum,alpha,beta);
                woT::task(owner, &implT:: template recursive_op<add_op>, op, cdata.key0);

            }

            this->compressed=false;
            if (fence) world.gop.fence();
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

        std::shared_ptr<FunctionFunctorInterface<T,NDIM> > get_functor() const {
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

        bool& is_on_demand() {return on_demand;};
        const bool& is_on_demand() const {return on_demand;};

        TensorType get_tensor_type() const {return targs.tt;}
        TensorArgs get_tensor_args() const {return targs;}

        double get_thresh() const {return thresh;}

        void set_thresh(double value) {thresh = value;}

        bool get_autorefine() const {return autorefine;}

        void set_autorefine(bool value) {autorefine = value;}

        int get_k() const {return k;}

        const dcT& get_coeffs() const {return coeffs;}

        dcT& get_coeffs() {return coeffs;}

        const FunctionCommonData<T,NDIM>& get_cdata() const {return cdata;}

        Void accumulate_timer(const double time) const {
            timer_accumulate.accumulate(time);
            return None;
        }


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
            if (world.rank() == coeffs.owner(cdata.key0)) {
                if (is_compressed()) {
                    truncate_spawn(cdata.key0,tol);
                } else {
                    truncate_reconstructed_spawn(cdata.key0,tol);
                }
            }
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

        /// convert a number [0,limit] to a hue color code [blue,red],
        /// or, if log is set, a number [1.e-10,limit]
        struct do_convert_to_color {
            double limit;
            bool log;
            static const double lower=1.e-10;
            do_convert_to_color() {};
            do_convert_to_color(const double limit, const bool log) : limit(limit), log(log) {}
            double operator()(double val) const {
                double color=0.0;

                if (log) {
                    double val2=log10(val) - log10(lower);        // will yield >0.0
                    double upper=log10(limit) -log10(lower);;
                    val2=0.7-(0.7/upper)*val2;
                    color= std::max(0.0,val2);
                    color= std::min(0.7,color);
                } else {
                    double hue=0.7-(0.7/limit)*(val);
                    color= std::max(0.0,hue);
                }
                return color;
            }
        };


        /// Print a plane ("xy", "xz", or "yz") containing the point x to file

        /// works for all dimensions; we walk through the tree, and if a leaf node
        /// inside the sub-cell touches the plane we print it in pstricks format
        void print_plane(const std::string filename, const std::string plane, const coordT& x_user) {

            // get the local information
            Tensor<double> localinfo=print_plane_local(plane,x_user);

            // lump all the local information together, and gather on node0
            std::vector<Tensor<double> > localinfo_vec(1,localinfo);
            std::vector<Tensor<double> > printinfo=world.gop.concat0(localinfo_vec);
            world.gop.fence();

            // do the actual print
            if (world.rank()==0) do_print_plane(filename,printinfo);

        }


        /// collect the data for a plot of the MRA structure locally on each node
        Tensor<double> print_plane_local(const std::string plane, const coordT& x_user) {

            // translate verbose plane to something computer-readable
            int dim0, dim1;
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

            // dimensions are: (# boxes)(hue, x lo left, y lo left, x hi right, y hi right)
            Tensor<double> plotinfo(coeffs.size(),5);
            long counter=0;

            // loop over local boxes, if the fit, add the info to the output tensor
            typename dcT::const_iterator end = coeffs.end();
            for (typename dcT::const_iterator it=coeffs.begin(); it!=end; ++it) {
                const keyT& key = it->first;
                const nodeT& node = it->second;

                // thisKeyContains ignores dim0 and dim1
                if (key.thisKeyContains(x_sim,dim0,dim1) and node.is_leaf() and (node.has_coeff())) {

                    Level n=key.level();
                    Vector<Translation,NDIM> l=key.translation();
                    // get the diametral edges of the node in the plotting plane
                    double scale=std::pow(0.5,double(n));
                    double xloleft = scale*l[dim0];
                    double yloleft = scale*l[dim1];
                    double xhiright = scale*(l[dim0]+1);
                    double yhiright = scale*(l[dim1]+1);

                    // do rank or do error
                    double color=0.0;
                    if (1) {

                        const double maxrank=40;
                        do_convert_to_color hue(maxrank,false);
                        color=hue(node.coeff().rank());
                    } else {

                        // Make quadrature rule of higher order
                        const int npt = cdata.npt + 1;
                        Tensor<double> qx, qw, quad_phi, quad_phiw, quad_phit;
                        FunctionCommonData<T,NDIM>::_init_quadrature(k+1, npt, qx, qw, quad_phi, quad_phiw, quad_phit);
                        do_err_box< FunctionFunctorInterface<T,NDIM> > op(this, this->get_functor().get(), npt, qx, quad_phit, quad_phiw);

                        do_convert_to_color hue(1000.0,true);
                        double error=op(it);
                        error=sqrt(error);//*pow(2,key.level()*6);
                        color=hue(error);
                    }

                    plotinfo(counter,0)=color;
                    plotinfo(counter,1)=xloleft;
                    plotinfo(counter,2)=yloleft;
                    plotinfo(counter,3)=xhiright;
                    plotinfo(counter,4)=yhiright;
                    ++counter;
                }
            }

            // shrink the info
            if (counter==0) plotinfo=Tensor<double>();
            else plotinfo=plotinfo(Slice(0,counter-1),Slice(_));
            return plotinfo;
        }

        /// print the MRA structure
        Void do_print_plane(const std::string filename, std::vector<Tensor<double> > plotinfo) {

            // invoke only on master node
            MADNESS_ASSERT(world.rank()==0);

            // prepare file
            FILE * pFile;
            pFile = fopen(filename.c_str(), "w");
            Tensor<double> cell=FunctionDefaults<NDIM>::get_cell();


            fprintf(pFile,"\\psset{unit=10cm}\n");
            fprintf(pFile,"\\begin{pspicture}(0,0)(1,1)\n");
            fprintf(pFile,"\\pslinewidth=0.005pt\n");

            for (std::vector<Tensor<double> >::const_iterator it=plotinfo.begin(); it!=plotinfo.end(); ++it) {

                Tensor<double> localinfo=*it;
                if (localinfo.has_data()) {

                    for (long i=0; i<localinfo.dim(0); ++i) {

                        fprintf(pFile,"\\newhsbcolor{mycolor}{%8.4f 1.0 0.7}\n",localinfo(i,0));
                        fprintf(pFile,"\\psframe["//linewidth=0.5pt,"
                                "fillstyle=solid,"
                                "fillcolor=mycolor]"
                                "(%12.8f,%12.8f)(%12.8f,%12.8f)\n",
                                localinfo(i,1),localinfo(i,2),localinfo(i,3),localinfo(i,4));
                    }
                }
            }


            fprintf(pFile,"\\end{pspicture}\n");
            fclose(pFile);

            return None;
        }


        /// Compute by projection the scaling function coeffs in specified box
        tensorT project(const keyT& key) const;

        /// Returns the truncation threshold according to truncate_method

        /// here is our handwaving argument:
        /// this threshold will give each FunctionNode an error of less than tol. The
        /// total error can then be as high as sqrt(#nodes) * tol. Therefore in order
        /// to account for higher dimensions: divide tol by about the root of number
        /// of siblings (2^NDIM) that have a large error when we refine along a deep
        /// branch of the tree.
        double truncate_tol(double tol, const keyT& key) const {
            const static double fac=1.0/std::pow(2,NDIM*0.5);
            tol*=fac;

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
        GenTensor<Q> coeffs2values(const keyT& key, const GenTensor<Q>& coeff) const {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            double scale = pow(2.0,0.5*NDIM*key.level())/sqrt(FunctionDefaults<NDIM>::get_cell_volume());
            return transform(coeff,cdata.quad_phit).scale(scale);
        }

        template <typename Q>
        GenTensor<Q> values2coeffs(const keyT& key, const GenTensor<Q>& values) const {
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
        GenTensor<Q> fcube_for_mul(const keyT& child, const keyT& parent, const GenTensor<Q>& coeff) const {
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


        /// Compute the function values for multiplication for child 

	/// call as fcube_for_mul3(key,key); walk up the tree to find an appropriate node 
	/// @param[in] child	the key for which we want the values
	/// @param[in] parent	the key where we start looking for it
	/// @param[out] 	values at the quadrature points of child
        Future<coeffT> fcube_for_mul3(const keyT child, const keyT parent) {
            
	    MADNESS_ASSERT(parent.is_valid());
            typedef std::pair<keyT,coeffT> argT;
	    if (coeffs.probe(parent)) {
        	const argT result= std::pair<keyT,coeffT>(parent,coeffs.find(parent).get()->second.coeff());
		return woT::task(coeffs.owner(parent), &implT::fcube_for_mul2, child, result,TaskAttributes::hipri());
	    } else {
                const keyT grandparent = parent.parent();
	        return woT::task(coeffs.owner(grandparent),&implT::fcube_for_mul3,child,grandparent,TaskAttributes::hipri());
	    }

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
                    TensorArgs full(-1.0,TT_FULL);
                    change_tensor_type(node.coeff(),full);
                    op(parent, node.coeff().full_tensor());
                    change_tensor_type(node.coeff(),targs);
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



        /// Unary operation applied inplace to the coefficients WITHOUT refinement, optional fence
        template <typename opT>
        void flo_unary_op_node_inplace(const opT& op, bool fence) const {

            typedef Range<typename dcT::const_iterator> rangeT;
            typedef do_unary_op_value_inplace<opT> xopT;
            world.taskq.for_each<rangeT,opT>(rangeT(coeffs.begin(), coeffs.end()), op);
            if (fence)
                 world.gop.fence();
        }

        /// truncate tree at a certain level
        Void erase(const Level& max_level) {
            this->make_redundant(true);

            typename dcT::iterator end = coeffs.end();
            for (typename dcT::iterator it= coeffs.begin(); it!=end; ++it) {
                keyT key=it->first;
                nodeT& node=it->second;
                if (key.level()>max_level) coeffs.erase(key);
                if (key.level()==max_level) node.set_has_children(false);
            }
            this->undo_redundant(true);
            return None;
        };


        /// Returns some asymmetry measure ... no comms
        double check_symmetry_local() const {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            typedef Range<typename dcT::const_iterator> rangeT;
            return world.taskq.reduce<double,rangeT,do_check_symmetry_local>(rangeT(coeffs.begin(),coeffs.end()),
                    do_check_symmetry_local(*this));
        }

        /// given an NS tree resulting from a convolution, truncate leafs if appropriate
        struct do_truncate_NS_leafs {
            typedef Range<typename dcT::iterator> rangeT;
            const implT* f;     // for calling its member functions

            do_truncate_NS_leafs(const implT* f) : f(f) {}

            bool operator()(typename rangeT::iterator& it) const {

                const keyT& key = it->first;
                nodeT& node = it->second;

                if (node.is_leaf() and node.coeff().has_data()) {
                    coeffT d = copy(node.coeff());
                    d(f->cdata.s0)=0.0;
                    const double error=d.normf();
                    const double tol=f->truncate_tol(f->get_thresh(),key);
                    if (error<tol) node.coeff()=copy(node.coeff()(f->cdata.s0));
                }
                return true;
            }
            template <typename Archive> void serialize(const Archive& ar) {}

        };

        struct do_stuff {
            typedef Range<typename dcT::iterator> rangeT;
            FunctionCommonData<T,NDIM> cdata;
            TensorType tt;
            int k;

            // constructor takes target precision
            do_stuff()
                : cdata(FunctionCommonData<T,NDIM>::get(FunctionDefaults<NDIM>::get_k()))
                , tt(FunctionDefaults<NDIM>::get_tensor_type())
                , k(FunctionDefaults<NDIM>::get_k())
            {
            }

            bool operator()(typename rangeT::iterator& it) const {

                nodeT& node = it->second;
                const keyT& key = it->first;
                if (0) print("key",key);
                if (node.coeff().rank()==0) {
//                    print("empty key in do_stuffq",key);
                    node.coeff()=coeffT(cdata.v2k,tt);
                } else {
//                    if (node.has_children()) {
                    MADNESS_ASSERT(node.coeff().dim(0)==k);
                        TensorType tt=node.coeff().tensor_type();
                        coeffT d(cdata.v2k,tt);
                        d(cdata.s0)+=node.coeff();
                        node.coeff()=d;
//                    }
                }
                return true;
            }
            template <typename Archive> void serialize(const Archive& ar) {}
        };



        /// remove all coefficients of internal nodes
        /// presumably to switch from redundant to reconstructed state
        struct remove_internal_coeffs {
             typedef Range<typename dcT::iterator> rangeT;

             /// constructor need impl for cdata
             remove_internal_coeffs() {}

             bool operator()(typename rangeT::iterator& it) const {

                 nodeT& node = it->second;
                 if (node.has_children()) node.clear_coeff();
                 return true;
             }
             template <typename Archive> void serialize(const Archive& ar) {}

        };


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
            do_reduce_rank(const double& thresh) {
                args.thresh=thresh;
            }

            //
            bool operator()(typename rangeT::iterator& it) const {

                nodeT& node = it->second;
                node.reduceRank(args.thresh);
                return true;
            }
            template <typename Archive> void serialize(const Archive& ar) {}
        };



        /// check symmetry wrt particle exchange
        struct do_check_symmetry_local {
            typedef Range<typename dcT::const_iterator> rangeT;
            const implT* f;
            do_check_symmetry_local() {}
            do_check_symmetry_local(const implT& f) : f(&f) {}

            /// return the norm of the difference of this node and its "mirror" node
            double operator()(typename rangeT::iterator& it) const {

                const keyT& key = it->first;
                const nodeT& fnode = it->second;

                if (f->world.size()>1) return 0.0;

                // exchange particles
                std::vector<long> map(NDIM);
                map[0]=3; map[1]=4; map[2]=5;
                map[3]=0; map[4]=1; map[5]=2;

                // make mapped key
                Vector<Translation,NDIM> l;
                for (std::size_t i=0; i<NDIM; ++i) l[map[i]] = key.translation()[i];
                const keyT mapkey(key.level(),l);

                // hope it's local
                MADNESS_ASSERT(f->get_coeffs().probe(mapkey));
                const nodeT& mapnode=f->get_coeffs().find(mapkey).get()->second;

                tensorT c1=fnode.coeff().full_tensor_copy();
                tensorT c2=mapnode.coeff().full_tensor_copy();

                if (c2.size()) c2 = copy(c2.mapdim(map));
                double norm=(c1-=c2).normf();
                return norm*norm;
            }

            double operator()(double a, double b) const {
                return (a+b);
            }

            template <typename Archive> void serialize(const Archive& ar) {
                MADNESS_EXCEPTION("no serialization of do_check_symmetry yet",1);
            }


        };

        /// map this on f
        struct do_mapdim {
            typedef Range<typename dcT::iterator> rangeT;

            std::vector<long> map;
            implT* f;

            do_mapdim() {};
            do_mapdim(const std::vector<long> map, implT& f) : map(map), f(&f) {}

            bool operator()(typename rangeT::iterator& it) const {

                const keyT& key = it->first;
                const nodeT& node = it->second;

                Vector<Translation,NDIM> l;
                for (std::size_t i=0; i<NDIM; ++i) l[map[i]] = key.translation()[i];
                tensorT c = node.full_tensor_copy();
                if (c.size()) c = copy(c.mapdim(map));
                coeffT cc(c,f->get_tensor_args());
                f->get_coeffs().replace(keyT(key.level(),l), nodeT(cc,node.has_children()));

                return true;
            }
            template <typename Archive> void serialize(const Archive& ar) {
                MADNESS_EXCEPTION("no serialization of do_mapdim",1);
            }

        };

        /// "put" this on g
        struct do_average {
            typedef Range<typename dcT::const_iterator> rangeT;

            implT* g;

            do_average() {}
            do_average(implT& g) : g(&g) {}

            /// iterator it points to this
            bool operator()(typename rangeT::iterator& it) const {

                const keyT& key = it->first;
                const nodeT& fnode = it->second;

                // fast return if rhs has no coeff here
                if (fnode.has_coeff()) {

                    // check if there is a node already existing
                    typename dcT::accessor acc;
                    if (g->get_coeffs().find(acc,key)) {
                        nodeT& gnode=acc->second;
                        if (gnode.has_coeff()) gnode.coeff()+=fnode.coeff();
                    } else {
                        g->get_coeffs().replace(key,fnode);
                    }
                }

                return true;
            }
            template <typename Archive> void serialize(const Archive& ar) {}
        };

        /// change representation of nodes' coeffs to low rank, optional fence
        struct do_change_tensor_type {
            typedef Range<typename dcT::iterator> rangeT;

            // threshold for rank reduction / SVD truncation
            TensorArgs targs;

            // constructor takes target precision
            do_change_tensor_type() {}
            do_change_tensor_type(const TensorArgs& targs) : targs(targs) {}

            //
            bool operator()(typename rangeT::iterator& it) const {

                nodeT& node = it->second;
                change_tensor_type(node.coeff(),targs);
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
                    const TensorArgs full_args(-1.0,TT_FULL);
                    change_tensor_type(node.coeff(),full_args);
//                	node.node_to_full_rank();
//                    coeffT& t= node.coeff();
                    tensorT& t= node.coeff().full_tensor();
                    //double before = t.normf();
                    tensorT values = impl->fcube_for_mul(key, key, t);
                    op(key, values);
                    double scale = pow(0.5,0.5*NDIM*key.level())*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
                    t = transform(values,impl->cdata.quad_phiw).scale(scale);
                    change_tensor_type(node.coeff(),impl->get_tensor_args());
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
                    const GenTensor<R>& r = it->second.coeff();
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

        /// Refine multiple functions down to the same finest level

        /// @param[v] is the vector of functions we are refining.
        /// @param[key] is the current node.
        /// @param[c] is the vector of coefficients passed from above.
        Void refine_to_common_level(const std::vector<FunctionImpl<T,NDIM>*>& v,
                                    const std::vector<tensorT>& c,
                                    const keyT key) {
            if (key == cdata.key0 && coeffs.owner(key)!=world.rank()) return None;

            // First insert coefficients from above ... also get write accessors here
            ScopedArray<typename dcT::accessor> acc(new typename dcT::accessor[v.size()]);
            for (unsigned int i=0; i<c.size(); i++) {
                MADNESS_ASSERT(v[i]->coeffs.get_pmap() == coeffs.get_pmap());
                MADNESS_ASSERT(v[i]->coeffs.owner(key) == world.rank());
                bool exists = ! v[i]->coeffs.insert(acc[i],key);
                if (c[i].size()) {
                    MADNESS_ASSERT(!exists);
                    acc[i]->second = nodeT(coeffT(c[i],targs),false);
                }
                else {
                    MADNESS_ASSERT(exists);
                }
            }

            // If everyone has coefficients we are done
            bool done = true;
            for (unsigned int i=0; i<v.size(); i++) {
                done &= acc[i]->second.has_coeff();
            }

            if (!done) {
                // Those functions with coefficients need to be refined down
                std::vector<tensorT> d(v.size());
                for (unsigned int i=0; i<v.size(); i++) {
                    if (acc[i]->second.has_coeff()) {
                        tensorT s(cdata.v2k);
//                        s(cdata.s0) = acc[i]->second.coeff()(___);
                        s(cdata.s0) = acc[i]->second.coeff().full_tensor_copy();
                        acc[i]->second.clear_coeff();
                        d[i] = unfilter(s);
                        acc[i]->second.set_has_children(true);
                    }
                }

                // Loop thru children and pass down
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT& child = kit.key();
                    std::vector<Slice> cp = child_patch(child);
                    std::vector<tensorT> childc(v.size());
                    for (unsigned int i=0; i<v.size(); i++) {
                        if (d[i].size()) childc[i] = copy(d[i](cp));
                    }
                    woT::task(coeffs.owner(child), &implT::refine_to_common_level, v, childc, child);
                }
            }

            return None;
        }

        template <typename opT>
        Void multiop_values_doit(const keyT& key, const opT& op, const std::vector<implT*>& v) {
            std::vector<tensorT> c(v.size());
            for (unsigned int i=0; i<v.size(); i++) {
                c[i] = coeffs2values(key, v[i]->coeffs.find(key).get()->second.coeff().full_tensor_copy()); // !!!!! gack
            }
            tensorT r = op(key, c);
            coeffs.replace(key, nodeT(coeffT(values2coeffs(key, r),targs),false));
            return None;
        }

        // assumes all functions have been refined down to the same level
        template <typename opT>
        void multiop_values(const opT& op, const std::vector<implT*>& v) {
            typename dcT::iterator end = v[0]->coeffs.end();
            for (typename dcT::iterator it=v[0]->coeffs.begin(); it!=end; ++it) {
                const keyT& key = it->first;
                if (it->second.has_coeff())
                    world.taskq.add(*this, &implT:: template multiop_values_doit<opT>, key, op, v);
                else
                    coeffs.replace(key, nodeT(coeffT(),true));
            }
            world.gop.fence();
        }

        /// Transforms a vector of functions left[i] = sum[j] right[j]*c[j,i] using sparsity
        template <typename Q, typename R>
        void vtransform(const std::vector< std::shared_ptr< FunctionImpl<R,NDIM> > >& vright,
                        const Tensor<Q>& c,
                        const std::vector< std::shared_ptr< FunctionImpl<T,NDIM> > >& vleft,
                        double tol,
                        bool fence) {
            for (unsigned int j=0; j<vright.size(); ++j) {
                world.taskq.add(*this, &implT:: template vtransform_doit<Q,R>, vright[j], copy(c(j,_)), vleft, tol);
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
        void tnorm(const coeffT& fcoeff, double* lo, double* hi) const {
            if (fcoeff.has_no_data()) {
                *lo=0.0;
                *hi=0.0;
            } else {

                double norm=fcoeff.normf();
                coeffT tlo=fcoeff(cdata.sh);
                *lo=tlo.normf();
                *hi=norm-*lo;
            }
        }


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
            return None;
        }

        /// After 1d push operator must sum coeffs down the tree to restore correct scaling function coefficients
        void sum_down(bool fence) {
            if (world.rank() == coeffs.owner(cdata.key0)) sum_down_spawn(cdata.key0, coeffT());

            if (fence) world.gop.fence();
        }

        /// perform this multiplication: h(1,2) = f(1,2) * g(1)
        template<size_t LDIM>
        struct multiply_op {

            typedef FunctionImpl<T,LDIM> implL;
            typedef std::pair<Key<LDIM>, ShallowNode<T,LDIM> > datumL;
            typedef std::pair<Key<NDIM>, ShallowNode<T,NDIM> > datumT;

            implT* h;     ///< the result function h(1,2) = f(1,2) * g(1)
            const implT* f;     ///< the function f(1,2) that will be multiplied with g(1)
            const implL* g;     ///< the function g(1) or g(2) will be multiplied with f(1,2)
            datumT fdatum;      ///< pointing to the most-leaf node of f
            datumL gdatum;      ///< pointing to the most-leaf node of g
            int particle;       ///< if g is g(1) or g(2)

            multiply_op()
                : particle(1) {
            };

            multiply_op(implT* h, const implT* f, const implL* g, const datumT& fdatum, const datumL& gdatum, const int particle)
                : h(h)
                , f(f)
                , g(g)
                , fdatum(fdatum)
                , gdatum(gdatum)
                , particle(particle) {
                MADNESS_ASSERT(gdatum.second.coeff().has_data());
            };

            /// return true if this will be a leaf node

            /// use generalization of tnorm for a GenTensor
            bool screen(coeffT& fcoeff, const coeffT& gcoeff) const {
                double glo, ghi, flo, fhi;
                MADNESS_ASSERT(gcoeff.tensor_type()==TT_FULL);
                g->tnorm(gcoeff.full_tensor(), &glo, &ghi);

                // this assumes intimate knowledge of how a GenTensor is organized!
                MADNESS_ASSERT(fcoeff.tensor_type()==TT_2D);
                const long rank=fcoeff.rank();
                const long maxk=fcoeff.dim(0);
                tensorT vec=fcoeff.config().ref_vector(particle-1).reshape(rank,maxk,maxk,maxk);
                for (long i=0; i<rank; ++i) {
                    double lo,hi;
                    tensorT c=vec(Slice(i,i),_,_,_).reshape(maxk,maxk,maxk);
                    g->tnorm(c, &lo, &hi);        // note we use g instead of h, since g is 3D
                    flo+=lo*fcoeff.config().weights(i);
                    fhi+=hi*fcoeff.config().weights(i);
                }
                double total_hi=glo*fhi + ghi*flo + fhi*ghi;
                return (total_hi<(h->get_thresh()));

            }


            /// apply this on a FunctionNode of f and g of Key key

            /// @param[in]  key key for FunctionNode in f and g, (g: broken into particles)
            /// @return <this node is a leaf, coefficients of this node>
            std::pair<bool,coeffT> operator()(const Key<NDIM>& key) const {

//                bool is_leaf=(not fdatum.second.has_children());
//                if (not is_leaf) return std::pair<bool,coeffT> (is_leaf,coeffT());

                // break key into particles (these are the child keys, with f/gdatum come the parent keys)
                Key<LDIM> key1,key2;
                key.break_apart(key1,key2);
                const Key<LDIM> gkey= (particle==1) ? key1 : key2;

                // get coefficients of this FunctionNode, or of a parent (if applicable)
                // iterators point to nodes in redundant representation
                const coeffT& fcoeff=fdatum.second.coeff();
                const coeffT& gcoeff=gdatum.second.coeff();

                // get coefficients of the actual FunctionNode
                coeffT coeff1=f->parent_to_child(fcoeff,fdatum.first,key);
                coeff1.normalize();
                const coeffT coeff2=g->parent_to_child(gcoeff,gdatum.first,gkey);

                coeffT hcoeff;

                bool is_leaf=screen(coeff1,coeff2);
                if (key.level()<2) is_leaf=false;

                if (is_leaf) {

                    // convert coefficients to values
                    coeffT hvalues=f->coeffs2values(key,coeff1);
                    coeffT gvalues=g->coeffs2values(gkey,coeff2);

                    // multiply one of the two vectors of f with g
                    // note shallow copy of Tensor<T>
                    MADNESS_ASSERT(hvalues.tensor_type()==TT_2D);
                    MADNESS_ASSERT(gvalues.tensor_type()==TT_FULL);
                    const long rank=hvalues.rank();
                    const long maxk=h->get_k();
                    MADNESS_ASSERT(maxk==coeff1.dim(0));
                    tensorT vec=hvalues.config().ref_vector(particle-1).reshape(rank,maxk,maxk,maxk);
                    for (long i=0; i<rank; ++i) {
                        tensorT c=vec(Slice(i,i),_,_,_);
                        c.emul(gvalues.full_tensor());
                    }

                    // convert values back to coefficients
                    hcoeff=h->values2coeffs(key,hvalues);
                }

                return std::pair<bool,coeffT> (is_leaf,hcoeff);
            }

            ///
            Future<multiply_op> make_child_op(const keyT& child) const {

                // break key into particles
                Key<LDIM> key1, key2;
                child.break_apart(key1,key2);
                const Key<LDIM> gkey= (particle==1) ? key1 : key2;

                // point to "outermost" leaf node
                Future<datumT> datum11=f->outermost_child(child,fdatum);
                Future<datumL> datum22=g->outermost_child(gkey,gdatum);

                return h->world.taskq.add(*const_cast<multiply_op *> (this), &multiply_op<LDIM>::make_op,
                        h,f,g,datum11,datum22,particle);
            }

            /// taskq-compatible constructor
            multiply_op make_op(implT* h, const implT* f, const implL* g,
                    const datumT& fdatum, const datumL& gdatum, const int particle) {
                return multiply_op(h,f,g,fdatum,gdatum,particle);
            }

            /// serialization
            template <typename Archive> void serialize(const Archive& ar) {
                 ar & h & f & g & fdatum & gdatum & particle;
            }

        };


        /// add two functions f and g: result=alpha * f  +  beta * g
        struct add_op {

            typedef std::pair<Key<NDIM>, ShallowNode<T,NDIM> > datumT;

            const implT* f;     ///< first addend
            const implT* g;     ///< second addend
            datumT fdatum;      ///< pointing to the outermost node of f
            datumT gdatum;      ///< pointing to the outermost node of g
            double alpha, beta; ///< prefactor for f, g

            add_op() {};
            add_op(const implT* f, const implT* g, const datumT& fdatum, const datumT& gdatum,
                    const double alpha, const double beta)
                : f(f)
                , g(g)
                , fdatum(fdatum)
                , gdatum(gdatum)
                , alpha(alpha)
                , beta(beta) {
            }

            /// if we are at the bottom of the trees, return the sum of the coeffs
            std::pair<bool,coeffT> operator()(const keyT& key) const {

                bool is_leaf=(fdatum.second.is_leaf() and gdatum.second.is_leaf());
                if (not is_leaf) return std::pair<bool,coeffT> (is_leaf,coeffT());

                coeffT fcoeff=f->parent_to_child(fdatum.second.coeff(),fdatum.first,key);
                coeffT gcoeff=g->parent_to_child(gdatum.second.coeff(),gdatum.first,key);
                fcoeff.gaxpy(alpha,gcoeff,beta);
                return std::pair<bool,coeffT> (is_leaf,fcoeff);
            }

            Future<add_op> make_child_op(const keyT& child) const {

                // point to "outermost" leaf node
                Future<datumT> ffdatum=f->outermost_child(child,fdatum);
                Future<datumT> ggdatum=g->outermost_child(child,gdatum);

                return f->world.taskq.add(*const_cast<add_op *> (this), &add_op::make_op,
                        f,g,ffdatum,ggdatum,alpha,beta);
            }

            /// taskq-compatible constructor
            add_op make_op(const implT* f,const implT* g,const datumT& datum1,const datumT& datum2,
                    const double alpha, const double beta) {
                return add_op(f,g,datum1,datum2,alpha,beta);
            }


            template <typename Archive> void serialize(const Archive& ar) {
                ar & f & g & fdatum & gdatum & alpha & beta;
            }


        };


        /// multiply f (a pair function of NDIM) with an orbital g (LDIM=NDIM/2)

        /// as in (with h(1,2)=*this) : h(1,2) = g(1) * f(1,2)
        /// use tnorm as a measure to determine if f (=*this) must be refined
        /// @param[in]  f           the NDIM function f=f(1,2)
        /// @param[in]  g           the LDIM function g(1) (or g(2))
        /// @param[in]  particle    1 or 2, as in g(1) or g(2)
        /// @return     this        the NDIM function h(1,2)
        template<size_t LDIM>
        Void multiply(const implT* f, const FunctionImpl<T,LDIM>* g, const int particle) {

            typedef std::pair<Key<NDIM>, ShallowNode<T,NDIM> > datumT;
            typedef std::pair<Key<LDIM>, ShallowNode<T,LDIM> > datumL;
            typedef FunctionImpl<T,LDIM> implL;

            if (world.rank() == g->get_coeffs().owner(g->cdata.key0)) {
                Future<datumT> datf=f->task(f->get_coeffs().owner(f->cdata.key0), &implT::find_datum,
                        f->cdata.key0,TaskAttributes::hipri());
                Future<datumL> datg=g->task(g->get_coeffs().owner(g->cdata.key0), &implL::find_datum,
                        g->cdata.key0,TaskAttributes::hipri());

                // have to wait for this, but should be local..
                datumT fdatum=datf.get();
                datumL gdatum=datg.get();

                ProcessID owner = coeffs.owner(cdata.key0);

                typedef multiply_op<LDIM> op_type;
                op_type multiply_op(this,f,g,fdatum,gdatum,particle);

                woT::task(owner, &implT:: template recursive_op<op_type>, multiply_op, cdata.key0);

            }

            this->compressed=false;
            return None;

        }


        /// perform this multiplication: h=f*g, with g being on-demand

        /// the criteria for refinement are:
        /// 1. at least a leaf node of f
        /// 2. compare to coeffs of parent node
        template<typename leaf_opT>
        struct mul_op {

            typedef std::pair<Key<NDIM>, ShallowNode<T,NDIM> > datumT;

            implT* h;           ///< the result function h = f * g
            const implT* f;     ///< the function f(1,2) that will be multiplied with g
            const implT* g;     ///< the function g (on-demand)
            datumT fdatum;      ///< pointing to the most-leaf node of f
            leaf_opT leaf_op;
            coeffT hparent;     ///< coeffs of h of the parent node

            mul_op() {
            };

            mul_op(implT* h, const implT* f, const implT* g, const datumT& fdatum, const leaf_opT& leaf_op,
                    const coeffT& hparent)
                : h(h), f(f), g(g), fdatum(fdatum), leaf_op(leaf_op), hparent(hparent) {
            };

            /// apply this on a FunctionNode of f and g of Key key

            /// @param[in]  key key for FunctionNode in f and g
            /// @return <this node is a leaf, coefficients of this node>
            std::pair<bool,coeffT> operator()(const Key<NDIM>& key) const {

                // pre-screen if this is a leaf (supposedly leaf_node of f)
                bool is_leaf=leaf_op(fdatum.first);
                if (not is_leaf) return std::pair<bool,coeffT> (is_leaf,coeffT());

#if 0
                // loop over all children, get and get the error (akin to project_refine_op)
                tensorT r = tensorT(h->cdata.v2k);
                for (KeyChildIterator<NDIM> it(key); it; ++it) {
                    const keyT& child = it.key();
                    r(h->child_patch(child)) = this->coeff(child);
                }
                // Filter then test difference coeffs at level n
                tensorT d = h->filter(r);
                tensorT hcoeff = copy(d(h->cdata.s0));
                d(h->cdata.s0) = 0.0;
                double dnorm = d.normf();
#else
//                const tensorT h_upsampled=h->upsample(key,hparent).full_tensor_copy();
                const tensorT hcoeff=this->coeff(key);
//                const tensorT diff=h_upsampled-hcoeff;
//                const double dnorm=diff.normf();
#endif
//                is_leaf=(dnorm<h->truncate_tol(h->get_thresh(),key.level()));
//                if (is_leaf) return std::pair<bool,coeffT> (is_leaf,coeffT(hcoeff,h->get_tensor_args()));
//                return std::pair<bool,coeffT> (is_leaf,coeffT());
                return std::pair<bool,coeffT> (is_leaf,coeffT(hcoeff,h->get_tensor_args()));
            }

        private:
            /// return the coeffs of h=f*g, for key, using fdatum of the parent key
            tensorT coeff(const keyT& key) const {

                const coeffT fcoeff=f->parent_to_child(fdatum.second.coeff(),fdatum.first,key);
                const tensorT gcoeff=g->project(key);

                // convert coefficients to values
                tensorT fvalues=f->coeffs2values(key,fcoeff).full_tensor_copy();
                tensorT gvalues=g->coeffs2values(key,gcoeff);

                // multiply f and g
                if (fvalues.has_data() and gvalues.has_data()) {
                    tensorT hvalues(h->cdata.vk,false);
                    TERNARY_OPTIMIZED_ITERATOR(T, hvalues, T, fvalues, T, gvalues, *_p0 = *_p1 * *_p2;);
                    tensorT hcoeff=h->values2coeffs(key,hvalues);
                    return hcoeff;
                } else {
                    return tensorT(h->cdata.vk);
                }
            }

        public:

            /// make operator for child of this
            Future<mul_op> make_child_op(const keyT& child) const {

                // point to "outermost" leaf node
                Future<datumT> fdatum1=f->outermost_child(child,fdatum);
                const coeffT hparent1=coeffT(this->coeff(child.parent()),h->get_tensor_args());

                return h->world.taskq.add(*const_cast<mul_op *> (this), &mul_op<leaf_opT>::make_op,
                        h,f,g,fdatum1,leaf_op,hparent1);
            }

            /// taskq-compatible constructor
            mul_op make_op(implT* h, const implT* f, const implT* g,
                    const datumT& fdatum, const leaf_opT& leaf_op, const coeffT& hparent) {
                return mul_op(h,f,g,fdatum,leaf_op,hparent);
            }

            /// serialization
            template <typename Archive> void serialize(const Archive& ar) {
                 ar & h & f & g & fdatum & leaf_op & hparent;
            }

        };


        /// multiply f with an on-demand function, invoked by result

        /// @param[in]  leaf_op     use this as a (pre-) measure to determine if f needs refinement
        /// @param[in]  f           the NDIM function f=f(1,2)
        /// @param[in]  g           the on-demand function (provides coeffs via FunctorInterface)
        /// @return     this        f * g
        template<typename leaf_opT>
        Void multiply(const leaf_opT& leaf_op, const implT* f, const implT* g, const bool fence) {

            typedef std::pair<Key<NDIM>, ShallowNode<T,NDIM> > datumT;

            if (world.rank() == g->get_coeffs().owner(g->cdata.key0)) {
                Future<datumT> datf=f->task(f->get_coeffs().owner(f->cdata.key0), &implT::find_datum,
                        f->cdata.key0,TaskAttributes::hipri());

                // have to wait for this, but should be local..
                datumT fdatum=datf.get();

                typedef mul_op<leaf_opT> op_type;
                op_type mul_op(this,f,g,fdatum,leaf_op,fdatum.second.coeff());

                ProcessID owner = coeffs.owner(cdata.key0);
                woT::task(owner, &implT:: template recursive_op<op_type>, mul_op, cdata.key0);

            }
            if (fence) world.gop.fence();
            this->compressed=false;
            return None;

        }


        /// Hartree product of two LDIM functions to yield a NDIM = 2*LDIM function
        template<size_t LDIM, typename leaf_opT>
        struct hartree_op {

            typedef hartree_op<LDIM,leaf_opT> this_type;
            typedef FunctionImpl<T,LDIM> implL;
            typedef FunctionNode<T,LDIM> nodeL;
            typedef std::pair<Key<LDIM>, ShallowNode<T,LDIM> > datumL;
            typedef std::pair<Key<NDIM>, ShallowNode<T,NDIM> > datumT;

            static const double safety=1.0;

            const FunctionImpl<T,NDIM>* result;       ///< where to construct the pair function
            impl_and_arg<T,LDIM> p1;
            impl_and_arg<T,LDIM> p2;
            leaf_opT leaf_op;                   ///< determine if a given node will be a leaf node

            // ctor
            hartree_op() {}
            hartree_op(const implT* result, const implL* p1, const implL* p2,
                    const datumL& datum1, const datumL& datum2, const leaf_opT& leaf_op)
                 : result(result)
                 , p1(impl_and_arg<T,LDIM>(p1,datum1))
                 , p2(impl_and_arg<T,LDIM>(p2,datum2))
                 , leaf_op(leaf_op)
            {
                MADNESS_ASSERT(LDIM+LDIM==NDIM);
            }
            hartree_op(const implT* result, const impl_and_arg<T,LDIM>& p1, const impl_and_arg<T,LDIM>& p2,
                    const leaf_opT& leaf_op)
                : result(result)
                , p1(p1)
                , p2(p2)
                , leaf_op(leaf_op)
            {
                MADNESS_ASSERT(LDIM+LDIM==NDIM);
            }

            std::pair<bool,coeffT> operator()(const Key<NDIM>& key) const {

                const coeffT fcoeff=p1.datum.second.coeff();
                const coeffT gcoeff=p2.datum.second.coeff();
                bool is_leaf=leaf_op(key,fcoeff,gcoeff);
                if (not is_leaf) return std::pair<bool,coeffT> (is_leaf,coeffT());

                // break key into particles (these are the child keys, with datum1/2 come the parent keys)
                Key<LDIM> key1,key2;
                key.break_apart(key1,key2);

                // iterators point to nodes in nonstandard representation: get the sum coeffs
                const coeffT s1=p1.datum.second.coeff()(p1.impl->cdata.s0);
                const coeffT s2=p2.datum.second.coeff()(p2.impl->cdata.s0);

                coeffT coeff1=p1.impl->parent_to_child(s1,p1.datum.first,key1);
                coeffT coeff2=p2.impl->parent_to_child(s2,p2.datum.first,key2);

                // new coeffs are simply the hartree/kronecker/outer product --
                // construct the new 6D SVD tensors directly
                const unsigned int dim=6;
                const unsigned int maxk=result->get_k();
                MADNESS_ASSERT(maxk==coeff1.dim(0));

                // normalize
                Tensor<double> weights(1);
                weights=1.0;
                double norm1=coeff1.normf();
                double norm2=coeff2.normf();
                weights*=norm1*norm2;
                coeff1.scale(1.0/norm1);
                coeff2.scale(1.0/norm2);

                // need full_tensor_copy() here, b/c coeff1/2 are in NS form
                const SRConf<T> srconf(weights,coeff1.full_tensor_copy().reshape(1,maxk*maxk*maxk),
                        coeff2.full_tensor_copy().reshape(1,maxk*maxk*maxk),dim,maxk);
                const coeffT coeff(srconf);

                // no post-determination
//                is_leaf=leaf_op(key,coeff);
                return std::pair<bool,coeffT>(is_leaf,coeff);
            }

            Future<hartree_op> make_child_op(const keyT& child) const {

                // break key into particles
                Key<LDIM> key1, key2;
                child.break_apart(key1,key2);

                // point to "outermost" leaf node
                Future<impl_and_arg<T,LDIM> > p11=p1.make_child(key1);
                Future<impl_and_arg<T,LDIM> > p22=p2.make_child(key2);
                return result->world.taskq.add(*const_cast<hartree_op *> (this),
                        &hartree_op<LDIM,leaf_opT>::make_op,
                        result,p11,p22,leaf_op);
            }

            this_type make_op(const implT* result, const impl_and_arg<T,LDIM>& p11, const impl_and_arg<T,LDIM>& p22,
                    const leaf_opT& leaf_op) {
                return hartree_op(result,p11,p22,leaf_op);
            }

            template <typename Archive> void serialize(const Archive& ar) {
//                ar & result & p1 & p2 & datum1 & datum2 & leaf_op;
                ar & result & p1 & p2 & leaf_op;
            }
        };



        /// walk down the tree and perform an operation on each node

        /// might or might not insert coefficients into the tree, depending
        /// on the leaf_op
        /// @param[in]  op  the operator that creates the coefficients,
        ///                 and determines if they are leaf nodes
        /// @param[in]  key current FunctionNode we are working on
        /// @return nothing, but will insert (non-) empty coefficients in this' tree
        template<typename opT>
        Void recursive_op(const opT& op, const keyT& key) {

            // op returns <is_leaf, coeff>
            std::pair<bool,coeffT> datum=op(key);
            const bool is_leaf=datum.first;
            const coeffT& coeff=datum.second;

            // insert result into this' tree
            const bool has_children=(not is_leaf);
            coeffs.replace(key,nodeT(coeff,has_children));

            // descend if needed
            if (has_children) {
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT& child = kit.key();
                    Future<opT> child_op=op.make_child_op(child);
                    woT::task(world.rank(), &implT:: template forward_op<opT>, child_op, child);
                }
            }
            return None;
        }

        /// walk down the tree and perform an operation on each node
        template<typename opT>
        Void forward_op(const opT& op, const keyT& key) {
            woT::task(coeffs.owner(key), &implT:: template recursive_op<opT>, op, key);
            return None;
        }

        /// given two functions of LDIM, perform the Hartree/Kronecker/outer product

        /// |Phi(1,2)> = |phi(1)> x |phi(2)>
        template<std::size_t LDIM, typename leaf_opT>
        Void hartree_product(const FunctionImpl<T,LDIM>* p1, const FunctionImpl<T,LDIM>* p2,
                const leaf_opT& leaf_op, bool fence) {
            MADNESS_ASSERT(p1->is_nonstandard());
            MADNESS_ASSERT(p2->is_nonstandard());

            typedef std::pair<Key<LDIM>, ShallowNode<T,LDIM> > datumL;
            typedef FunctionImpl<T,LDIM> implL;

            if (world.rank() == p1->get_coeffs().owner(p1->cdata.key0)) {
                Future<datumL> dat1=p1->task(p1->get_coeffs().owner(p1->cdata.key0), &FunctionImpl<T,LDIM>::find_datum,
                        p1->cdata.key0,TaskAttributes::hipri());
                Future<datumL> dat2=p2->task(p2->get_coeffs().owner(p2->cdata.key0), &FunctionImpl<T,LDIM>::find_datum,
                        p2->cdata.key0,TaskAttributes::hipri());

                // have to wait for this, but should be local..
                datumL datum1=dat1.get();
                datumL datum2=dat2.get();

                typedef hartree_op<LDIM,leaf_opT> op_type;
                op_type hartree_op(this,p1,p2,datum1,datum2,leaf_op);

                ProcessID owner = coeffs.owner(cdata.key0);
                woT::task(owner, &implT:: template recursive_op<op_type>, hartree_op, cdata.key0);

            }

            this->compressed=false;
            if (fence) world.gop.fence();
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
                    world.taskq.add(*this, &implT::do_diff1, D, f, key, left, center, right, TaskAttributes::hipri());
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
        std::pair<Key<NDIM>,ShallowNode<T,NDIM> > find_datum(keyT key) const {
        	MADNESS_ASSERT(coeffs.probe(key));
        	ShallowNode<T,NDIM> snode(coeffs.find(key).get()->second);
        	return std::pair<Key<NDIM>,ShallowNode<T,NDIM> >(key,snode);
        }

        /// walk down the tree until a leaf node is hit
        Future<std::pair<keyT,ShallowNode<T,NDIM> > >
        outermost_child(const keyT& key, const std::pair<keyT,ShallowNode<T,NDIM> >& in) const {

            const ShallowNode<T,NDIM>& node=in.second;
            if (not node.has_children()) return Future<std::pair<keyT,ShallowNode<T,NDIM> > >(in);
            return woT::task(coeffs.owner(key), &implT::find_datum, key, TaskAttributes::hipri());
        }



        /// given a ket and the 1- and 2-electron potentials, construct the function V phi
        template<typename opT, size_t LDIM>
        struct Vphi_op {

            typedef Vphi_op<opT,LDIM> this_type;
            typedef FunctionImpl<T,LDIM> implL;
            typedef std::pair<Key<LDIM>, ShallowNode<T,LDIM> > datumL;
            typedef std::pair<Key<NDIM>, ShallowNode<T,NDIM> > datumT;
            typedef impl_and_arg<T,NDIM> iaT;
            typedef impl_and_arg<T,LDIM> iaL;

//            implT* result;             ///< where to construct the V phi
            impl_and_arg<T,NDIM> result;    ///< where to construct Vphi, and parent node
            opT leaf_op;               ///< deciding if a given FunctionNode will be a leaf node
            impl_and_arg<T,NDIM> iaket;
            impl_and_arg<T,LDIM> iap1;
            impl_and_arg<T,LDIM> iap2;
            impl_and_arg<T,LDIM> iav1;
            impl_and_arg<T,LDIM> iav2;
            const implT* eri;          ///< holding the 2-electron potential generator (on-demand)

            // ctor
            Vphi_op() {}

            Vphi_op(const iaT& result, const opT& leaf_op, const implT* ket, const implT* eri,
                    const implL* p1, const implL* p2, const implL* v1, const implL* v2,
                    const datumL& datum_p1, const datumL& datum_p2,
                    const datumL& datum_v1, const datumL& datum_v2,
                    const datumT& datum_ket)
                : result(result)
                , leaf_op(leaf_op)
                , iaket(ket,datum_ket)
                , iap1(p1,datum_p1)
                , iap2(p2,datum_p2)
                , iav1(v1,datum_v1)
                , iav2(v2,datum_v2)
                , eri(eri) {
            }

            Vphi_op(const iaT& result, const opT& leaf_op, const impl_and_arg<T,NDIM>& iaket,
                    const impl_and_arg<T,LDIM>& iap1, const impl_and_arg<T,LDIM>& iap2,
                    const impl_and_arg<T,LDIM>& iav1, const impl_and_arg<T,LDIM>& iav2,
                    const implT* eri)
                : result(result)
                , leaf_op(leaf_op)
                , iaket(iaket)
                , iap1(iap1)
                , iap2(iap2)
                , iav1(iav1)
                , iav2(iav2)
                , eri(eri) {
            }


            /// assemble the coefficients
            std::pair<bool,coeffT> operator()(const Key<NDIM>& key) const {

                // pre-determination: fast return if possible
                bool is_leaf=leaf_op(key);
                if (not is_leaf) return std::pair<bool,coeffT> (is_leaf,coeffT());


                // break key into particles (these are the child keys, with datum1/2 come the parent keys)
                Key<LDIM> key1, key2;
                key.break_apart(key1,key2);

                // values for ket: get it from its function, or construct it using orbitals
                tensorT val_ket;
                if (iaket.impl) {
                    val_ket=iaket.impl->fcube_for_mul(key,iaket.datum.first,
                            iaket.datum.second.coeff()).full_tensor_copy();
                } else {
                    MADNESS_ASSERT(iap1.impl and iap2.impl);
                    hartree_leaf_op<T,NDIM> hlop(result.impl,result.impl->get_k());
                    hartree_op<LDIM,hartree_leaf_op<T,NDIM> > op(result.impl,iap1,iap2,hlop);
                    const coeffT coeff_ket=op(key).second;
                    val_ket=result.impl->coeffs2values(key,coeff_ket).full_tensor_copy();
                }
                MADNESS_ASSERT(val_ket.has_data());

                // values for 1e-potentials
                const coeffT val_pot1= (iav1.impl)
                        ? iav1.impl->fcube_for_mul(key1,iav1.datum.first,iav1.datum.second.coeff())
                        : coeffT();
                const coeffT val_pot2= (iav2.impl)
                        ? iav2.impl->fcube_for_mul(key2,iav2.datum.first,iav2.datum.second.coeff())
                        : coeffT();

                // values for eri
                tensorT val_eri=tensorT();
                if (eri and eri->get_functor()->provides_coeff()) {
                    val_eri=eri->coeffs2values(key,eri->get_functor()->coeff(key)).full_tensor();
                } else if (eri) {
                    val_eri=tensorT(result.impl->cdata.vk);
                    result.impl->fcube(key, *(eri->get_functor()), eri->cdata.quad_x, val_eri);
                }

                // assemble all contributions
                tensorT coeff_v;
                if (eri or iav1.impl or iav2.impl) coeff_v=tensorT(result.impl->cdata.vk);      // 6D, all addends

                // include the one-electron potential
                if (val_pot1.has_data() and val_pot2.has_data()) {

                    FunctionCommonData<T,LDIM> cdataL=iav1.impl->cdata;
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
                if (val_eri.has_data()) coeff_v+=val_eri;

                // check for need of refinement
                bool needs_refinement=false;
                if (coeff_v.has_data()) {

                    double lo_ket, hi_ket, lo_pot, hi_pot;
                    result.impl->tnorm(val_ket, &lo_ket, &hi_ket);
                    result.impl->tnorm(coeff_v, &lo_pot, &hi_pot);
                    double test = lo_ket*hi_pot + lo_pot*hi_ket + hi_ket*hi_pot;
                    //print("autoreftest",key,thresh,truncate_tol(thresh, key),lo,hi,test);
                    needs_refinement= (test> result.impl->truncate_tol(result.impl->get_thresh(), key));

                }

                // multiply potential with ket
                if (coeff_v.has_data()) {
                    tensorT tcube(result.impl->cdata.vk,false);
                    TERNARY_OPTIMIZED_ITERATOR(T, tcube, T, val_ket, T, coeff_v, *_p0 = *_p1 * *_p2;);
                    val_ket=tcube;
                }

                const TensorArgs& targs=result.impl->get_tensor_args();
                const coeffT coeff_ket=coeffT(result.impl->values2coeffs(key,val_ket),targs);

                // post-determination
                is_leaf=leaf_op(key,coeff_ket);

                // compare to parent
                error_leaf_op<T,NDIM> elop(result.impl);
                is_leaf=elop(key,coeff_ket,result.datum.second.coeff());
                return std::pair<bool,coeffT> (is_leaf,coeff_ket);
            }

            /// given the current operator, construct the child operator, which means
            /// to provide pointers to the appropriate nodes from which the new
            /// result node will be constructed
            Future<Vphi_op> make_child_op(const keyT& child) const {

                // break key into particles
                Key<LDIM> key1, key2;
                child.break_apart(key1,key2);

                Future<iaT> result_=result.make_child(child.parent());
                Future<iaT> iaket_=iaket.make_child(child);
                Future<iaL> iap1_=iap1.make_child(key1);
                Future<iaL> iap2_=iap2.make_child(key2);
                Future<iaL> iav1_=iav1.make_child(key1);
                Future<iaL> iav2_=iav2.make_child(key2);

                // wait for the nodes to arrive, and construct a new operator
                return result.impl->world.taskq.add(*const_cast<Vphi_op *> (this), &this_type::make_op,
                        result_,leaf_op,iaket_,iap1_,iap2_,iav1_,iav2_,eri);
            }

            /// taskq-compatible child constructor
            Vphi_op make_op(const iaT& result, const opT& leaf_op, const impl_and_arg<T,NDIM>& iaket,
                    const impl_and_arg<T,LDIM>& iap1, const impl_and_arg<T,LDIM>& iap2,
                    const impl_and_arg<T,LDIM>& iav1, const impl_and_arg<T,LDIM>& iav2,
                    const implT* eri) const {
                return Vphi_op(result,leaf_op,iaket,iap1,iap2,iav1,iav2,eri);
            }

            /// serialize this (needed for use in recursive_op)
            template <typename Archive> void serialize(const Archive& ar) {
                ar & iaket & eri & result & leaf_op & iap1 & iap2 & iav1 & iav2;
            }
        };


        /// assemble the function V*phi using V and phi given from the functor

        /// this function must have been constructed using the CompositeFunctorInterface.
        /// The interface provides one- and two-electron potentials, and the ket, which are
        /// assembled to give V*phi. The MRA structure of the result is the same as the
        /// MRA structure of ket as given by the functor.
        /// @param[in]  leaf_op  operator to decide if a given node is a leaf node
        /// @param[in]  fence   global fence
        template<typename opT>
        void make_Vphi(const opT& leaf_op, const bool fence=true) {

            const size_t LDIM=3;
            typedef std::pair<Key<NDIM>, ShallowNode<T,NDIM> > datumT;
            typedef std::pair<Key<LDIM>, ShallowNode<T,LDIM> > datumL;

            CompositeFunctorInterface<T,NDIM,LDIM>* func=
                    dynamic_cast<CompositeFunctorInterface<T,NDIM,LDIM>* >(&(*functor));
            coeffs.clear();
            const keyT& key0=cdata.key0;


            const FunctionImpl<T,NDIM>* ket=func->impl_ket.get();
            const FunctionImpl<T,NDIM>* eri=func->impl_eri.get();
            const FunctionImpl<T,LDIM>* v1=func->impl_m1.get();
            const FunctionImpl<T,LDIM>* v2=func->impl_m2.get();
            const FunctionImpl<T,LDIM>* p1=func->impl_p1.get();
            const FunctionImpl<T,LDIM>* p2=func->impl_p2.get();

            if (world.rank() == coeffs.owner(key0)) {
                Future<datumL> dat1= (v1)
                        ? v1->task(v1->get_coeffs().owner(v1->cdata.key0),
                                &FunctionImpl<T,LDIM>::find_datum,v1->cdata.key0,TaskAttributes::hipri())
                        : Future<datumL>(datumL());
                Future<datumL> dat2= (v2)
                        ? v2->task(v2->get_coeffs().owner(v2->cdata.key0),
                                &FunctionImpl<T,LDIM>::find_datum,v2->cdata.key0,TaskAttributes::hipri())
                        : Future<datumL>(datumL());
                Future<datumL> datp1= (p1)
                        ? p1->task(p1->get_coeffs().owner(p1->cdata.key0),
                                &FunctionImpl<T,LDIM>::find_datum,p1->cdata.key0,TaskAttributes::hipri())
                        : Future<datumL>(datumL());
                Future<datumL> datp2= (p2)
                        ? p2->task(p2->get_coeffs().owner(p2->cdata.key0),
                                &FunctionImpl<T,LDIM>::find_datum,p2->cdata.key0,TaskAttributes::hipri())
                        : Future<datumL>(datumL());

                Future<datumT> dat_ket= (ket)
                        ? ket->task(ket->get_coeffs().owner(key0),
                                &implT::find_datum,key0,TaskAttributes::hipri())
                        : Future<datumT>(datumT());
                print("using walker for make_Vphi");
                if (eri) MADNESS_ASSERT(eri->is_on_demand());

                // have to wait for this..
                datumL datum1=dat1.get();
                datumL datum2=dat2.get();
                datumL datum_p1=datp1.get();
                datumL datum_p2=datp2.get();
                datumT datum_ket=dat_ket.get();

                // insert an empty internal node for comparison
                this->coeffs.replace(key0,nodeT(coeffT(),true));
                impl_and_arg<T,NDIM> iaresult(this,this->find_datum(key0));

                typedef Vphi_op<opT,LDIM> op_type;
                op_type op(iaresult,leaf_op,ket,eri,p1,p2,v1,v2,datum_p1,datum_p2,datum1,datum2,datum_ket);

                ProcessID owner = coeffs.owner(cdata.key0);
                woT::task(owner, &implT:: template recursive_op<op_type>, op, cdata.key0);

            }
            world.gop.fence();

            // remove internal coefficients
            this->redundant=true;
            this->undo_redundant(false);

            // set right state
            this->compressed=false;
            this->on_demand=false;
            this->redundant=false;
            this->nonstandard=false;
            this->unset_functor();
            if (fence) world.gop.fence();

        }


        /// Permute the dimensions of f according to map, result on this
        void mapdim(const implT& f, const std::vector<long>& map, bool fence) {

            PROFILE_MEMBER_FUNC(FunctionImpl);
            const_cast<implT*>(&f)->flo_unary_op_node_inplace(do_mapdim(map,*this),fence);

        }


        /// take the average of two functions, similar to: this=0.5*(this+rhs)

        /// works in either basis and also in nonstandard form
        void average(const implT& rhs) {

            rhs.flo_unary_op_node_inplace(do_average(*this),true);
            this->scale_inplace(0.5,true);
            flo_unary_op_node_inplace(do_reduce_rank(targs),true);
        }

        /// change the tensor type of the coefficients in the FunctionNode

        /// @param[in]  targs   target tensor arguments (threshold and full/low rank)
        void change_tensor_type1(const TensorArgs& targs, bool fence) {
            flo_unary_op_node_inplace(do_change_tensor_type(targs),fence);
        }

        /// reduce the rank of the coefficients tensors

        /// @param[in]  targs   target tensor arguments (threshold and full/low rank)
        void reduce_rank(const TensorArgs& targs, bool fence) {
            flo_unary_op_node_inplace(do_reduce_rank(targs),fence);
        }


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

        coeffT filter(const coeffT& s) const {
            coeffT result=transform(s,cdata.hgT);
        	return result;
        }

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

        coeffT unfilter(const coeffT& s) const {
            return transform(s,cdata.hg);
        }

        /// downsample the sum coefficients of level n+1 to sum coeffs on level n

        /// specialization of the filter method, will yield only the sum coefficients
        /// @param[in]  key key of level n
        /// @param[in]  v   vector of sum coefficients of level n+1
        /// @param[in]  args    TensorArguments for possible low rank approximations
        /// @return     sum coefficients on level n in full tensor format
        tensorT downsample(const keyT& key, const std::vector< Future<coeffT > >& v) const {

            tensorT result(cdata.vk);

            // the twoscale coefficients: for downsampling use h0/h1; see Alpert Eq (3.34a)
            const tensorT h[2] = {cdata.h0T, cdata.h1T};
            tensorT matrices[NDIM];

            // loop over all child nodes, transform and accumulate
            long i=0;
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {

                // get the appropriate twoscale coefficients for each dimension
                for (size_t ii=0; ii<NDIM; ++ii) matrices[ii]=h[kit.key().translation()[ii]%2];

                // transform and accumulate on the result
                result+=general_transform(v[i].get(),matrices).full_tensor_copy();

            }
            return result;
        }

        /// upsample the sum coefficients of level 1 to sum coeffs on level n+1

        /// specialization of the unfilter method, will transform only the sum coefficients
        /// @param[in]  key     key of level n+1
        /// @param[in]  coeff   sum coefficients of level n (does NOT belong to key!!)
        /// @param[in]  args    TensorArguments for possible low rank approximations
        /// @return     sum     coefficients on level n+1
        coeffT upsample(const keyT& key, const coeffT& coeff) const {

            // the twoscale coefficients: for upsampling use h0/h1; see Alpert Eq (3.35a/b)
            // note there are no difference coefficients; if you want that use unfilter
            const tensorT h[2] = {cdata.h0, cdata.h1};
            tensorT matrices[NDIM];

            // get the appropriate twoscale coefficients for each dimension
            for (size_t ii=0; ii<NDIM; ++ii) matrices[ii]=h[key.translation()[ii]%2];

            // transform and accumulate on the result
            const coeffT result=general_transform(coeff,matrices);
            return result;
        }


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

        /// sum all the contributions from all scales after applying an operator in mod-NS form
        void trickle_down(bool fence) {
            MADNESS_ASSERT(is_redundant());
            nonstandard = compressed = redundant = false;
            this->print_size("in trickle_down");
            if (world.rank() == coeffs.owner(cdata.key0))
                woT::task(world.rank(), &implT::trickle_down_op, cdata.key0,coeffT());
            if (fence)
                world.gop.fence();
        }

        /// sum all the contributions from all scales after applying an operator in mod-NS form

        /// cf reconstruct_op
        Void trickle_down_op(const keyT& key, const coeffT& s) {
            // Note that after application of an integral operator not all
            // siblings may be present so it is necessary to check existence
            // and if absent insert an empty leaf node.
            //
            // If summing the result of an integral operator (i.e., from
            // non-standard form) there will be significant scaling function
            // coefficients at all levels and possibly difference coefficients
            // in leaves, hence the tree may refine as a result.
            typename dcT::iterator it = coeffs.find(key).get();
            if (it == coeffs.end()) {
                coeffs.replace(key,nodeT(coeffT(),false));
                it = coeffs.find(key).get();
            }
            nodeT& node = it->second;

            // The integral operator will correctly connect interior nodes
            // to children but may leave interior nodes without coefficients
            // ... but they still need to sum down so just give them zeros
            if (node.coeff().has_no_data()) node.coeff()=coeffT(cdata.vk,targs);

//            if (node.has_children() || node.has_coeff()) { // Must allow for inconsistent state from transform, etc.
            if (node.has_children()) { // Must allow for inconsistent state from transform, etc.
                coeffT d = node.coeff();
                if (key.level() > 0) d += s; // -- note accumulate for NS summation
                node.clear_coeff();
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT& child = kit.key();
                    coeffT ss= upsample(child,d);
                    ss.reduceRank(thresh);
                    PROFILE_BLOCK(recon_send);
                    woT::task(coeffs.owner(child), &implT::trickle_down_op, child, ss);
                }
            }
            else {
                node.coeff()+=s;
                node.coeff().reduceRank(thresh);
            }
            return None;
        }

        void reconstruct(bool fence) {
            // Must set true here so that successive calls without fence do the right thing
            MADNESS_ASSERT(not is_redundant());
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
        /// @param[in] redundant    keep only sum coeffs at all levels, discard difference coeffs
        void compress(bool nonstandard, bool keepleaves, bool redundant, bool fence) {
            MADNESS_ASSERT(not is_redundant());
            // Must set true here so that successive calls without fence do the right thing
            this->compressed = true;
            this->nonstandard = nonstandard;
            this->redundant = redundant;

            // these two are exclusive
            MADNESS_ASSERT(not (redundant and nonstandard));
            // otherwise we loose information
            if (redundant) {MADNESS_ASSERT(keepleaves);}

//            this->print_tree();
            if (world.rank() == coeffs.owner(cdata.key0)) {

           		compress_spawn(cdata.key0, nonstandard, keepleaves, redundant);
            }
            if (fence)
                world.gop.fence();
        }

        // Invoked on node where key is local
//        Future<tensorT> compress_spawn(const keyT& key, bool nonstandard, bool keepleaves);
        Future<coeffT > compress_spawn(const keyT& key, bool nonstandard, bool keepleaves, bool redundant);

        /// convert this to redundant, i.e. have sum coefficients on all levels
        void make_redundant(const bool fence) {

            // fast return if possible
        	if (is_redundant()) return;

            // this is easy: just get rid of difference coefficients
            if (is_nonstandard()) {
                redundant=true;
                compressed=false;
                nonstandard=false;
                flo_unary_op_node_inplace(do_keep_sum_coeffs(this),fence);

            } else {
        		if (is_compressed()) reconstruct(true);
        		compress(false,true,true,fence);
        		compressed=false;
        	}

        }

        /// convert this from redundant to standard reconstructed form
        void undo_redundant(const bool fence) {

            if (!is_redundant()) return;
            redundant = compressed = nonstandard = false;
            flo_unary_op_node_inplace(remove_internal_coeffs(),fence);
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

        /// truncate using a tree in reconstructed form

        /// must be invoked where key is local
        Future<coeffT> truncate_reconstructed_spawn(const keyT& key, const double tol) {
            MADNESS_ASSERT(coeffs.probe(key));
            nodeT& node = coeffs.find(key).get()->second;

            // if this is a leaf node just return the sum coefficients
            if (not node.has_children()) return Future<coeffT>(node.coeff());

            // if this is an internal node, wait for all the children's sum coefficients
            // and use them to determine if the children can be removed
            std::vector<Future<coeffT> > v = future_vector_factory<coeffT>(1<<NDIM);
            int i=0;
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
                v[i] = woT::task(coeffs.owner(kit.key()), &implT::truncate_reconstructed_spawn, kit.key(),tol,TaskAttributes::hipri());
            }

            // will return (possibly empty) sum coefficients
            return woT::task(world.rank(),&implT::truncate_reconstructed_op,key,v,tol,TaskAttributes::hipri());

        }

        /// given the sum coefficients of all children, truncate or not

        /// @return     new sum coefficients (empty if internal, not empty, if new leaf); might delete its children
        coeffT truncate_reconstructed_op(const keyT& key, const std::vector< Future<coeffT > >& v, const double tol) {

            MADNESS_ASSERT(coeffs.probe(key));

            // the sum coefficients might be empty, which means they come from an internal node
            // and we must not truncate; so just return empty coeffs again
            for (size_t i=0; i<v.size(); ++i) if (v[i].get().has_no_data()) return coeffT();

            typename dcT::accessor acc;
            MADNESS_ASSERT(coeffs.find(acc, key));

            // the sum coefficients on this level, and their norm
            const tensorT s=downsample(key,v);
            const double snorm=s.normf();

            // get the norm of all child coefficients
            double dnorm=0.0;
            for (size_t i=0; i<v.size(); ++i) {
                const double d=v[i].get().normf();
                dnorm+=d*d;
            }

            // the error; equivalent to the norm of the wavelet coefficients
            const double error=sqrt(dnorm-snorm*snorm);

            nodeT& node = coeffs.find(key).get()->second;

            if (error < truncate_tol(tol,key)) {
                node.set_has_children(false);
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    coeffs.erase(kit.key());
                }
                // "replace" children with new sum coefficients
                coeffT ss=coeffT(s,targs);
                acc->second.set_coeff(ss);
                return ss;
            } else {
                return coeffT();
            }
        }

        /// calculate the wavelet coefficients using the sum coefficients of all child nodes

        /// @param[in] key 	this's key
        /// @param[in] v 	sum coefficients of the child nodes
        /// @param[in] nonstandard  keep the sum coefficients with the wavelet coefficients
        /// @param[in] redundant    keep only the sum coefficients, discard the wavelet coefficients
        /// @return 		the sum coefficients
        coeffT compress_op(const keyT& key, const std::vector< Future<coeffT > >& v, bool nonstandard, bool redundant) {
            PROFILE_MEMBER_FUNC(FunctionImpl);

            double cpu0=cpu_time();
            // Copy child scaling coeffs into contiguous block
            tensorT d(cdata.v2k);
//            coeffT d(cdata.v2k,targs);
            int i=0;
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
//                d(child_patch(kit.key())) += v[i].get();
                d(child_patch(kit.key())) += v[i].get().full_tensor_copy();
            }

            d = filter(d);
            double cpu1=cpu_time();
            timer_filter.accumulate(cpu1-cpu0);
            cpu0=cpu1;

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

            // tighter thresh for internal nodes
            TensorArgs targs2=targs;
            targs2.thresh*=0.1;

            // need the deep copy for contiguity
            coeffT ss=coeffT(copy(d(cdata.s0)),targs2);

            if (key.level()> 0 && !nonstandard)
                d(cdata.s0) = 0.0;

            // insert either sum or difference coefficients
            if (redundant) {
                acc->second.set_coeff(ss);
            } else {
                coeffT dd=coeffT(d,targs2);
                acc->second.set_coeff(dd);
            }
            cpu1=cpu_time();
            timer_compress_svd.accumulate(cpu1-cpu0);

            // return sum coefficients
            return ss;
        }

        /// similar to compress_op, but insert only the sum coefficients in the tree

        /// @param[in] key  this's key
        /// @param[in] v    sum coefficients of the child nodes
        /// @return         the sum coefficients
        coeffT make_redundant_op(const keyT& key, const std::vector< Future<coeffT > >& v) {

            // get the sum coefficients of this level given the sum coefficients of level n+1
            TensorArgs targs2=targs;
            targs2.thresh*=0.1;
            coeffT s(this->downsample(key,v),targs2);

            // insert sum coefficients into tree
            typename dcT::accessor acc;
            MADNESS_ASSERT(coeffs.find(acc, key));
            MADNESS_ASSERT(not (acc->second.has_coeff()));
            acc->second.set_coeff(s);

            return s;
        }

        /// Changes non-standard compressed form to standard compressed form
        void standard(bool fence) {

            flo_unary_op_node_inplace(do_standard(this),fence);
            nonstandard = false;
        }

        /// Changes non-standard compressed form to standard compressed form
        struct do_standard {
            typedef Range<typename dcT::iterator> rangeT;

            // threshold for rank reduction / SVD truncation
            implT* impl;

            // constructor takes target precision
            do_standard() {}
            do_standard(implT* impl) : impl(impl) {}

            //
            bool operator()(typename rangeT::iterator& it) const {

                const keyT& key = it->first;
                nodeT& node = it->second;
                if (key.level()> 0 && node.has_coeff()) {
                    if (node.has_children()) {
                        // Zero out scaling coeffs
                        node.coeff()(impl->cdata.s0)=0.0;
                        node.reduceRank(impl->targs.thresh);
                    } else {
                        // Deleting both scaling and wavelet coeffs
                        node.clear_coeff();
                    }
                }
                return true;
            }
            template <typename Archive> void serialize(const Archive& ar) {
                MADNESS_EXCEPTION("no serialization of do_standard",1);
            }
        };


        /// laziness
        template<size_t OPDIM>
        struct do_op_args {
            Key<OPDIM> key,d;
            keyT dest;
            double tol, fac, cnorm;
            do_op_args() {}
            do_op_args(const Key<OPDIM>& key, const Key<OPDIM>& d, const keyT& dest, double tol, double fac, double cnorm)
                    : key(key), d(d), dest(dest), tol(tol), fac(fac), cnorm(cnorm) {}
            template <class Archive>
            void serialize(Archive& ar) {
                ar & archive::wrap_opaque(this,1);
            }
        };

        /// for fine-grain parallelism: call the apply method of an operator in a separate task

        /// @param[in]  op      the operator working on our function
        /// @param[in]  c       full rank tensor holding the NS coefficients
        /// @param[in]  args    laziness holding norm of the coefficients, displacement, destination, ..
        /// @return     nothing, but accumulate the result tensor into the destination node
        template <typename opT, typename R, size_t OPDIM>
        Void do_apply_kernel(const opT* op, const Tensor<R>& c, const do_op_args<OPDIM>& args) {

            tensorT result = op->apply(args.key, args.d, c, args.tol/args.fac/args.cnorm);

            // Screen here to reduce communication cost of negligible data
            // and also to ensure we don't needlessly widen the tree when
            // applying the operator
            if (result.normf()> 0.3*args.tol/args.fac) {
                // OPTIMIZATION NEEDED HERE ... CHANGING THIS TO TASK NOT SEND REMOVED
                // BUILTIN OPTIMIZATION TO SHORTCIRCUIT MSG IF DATA IS LOCAL
                Future<double> time=coeffs.task(args.dest, &nodeT::accumulate2, result, coeffs, args.dest, TaskAttributes::hipri());
                woT::task(world.rank(),&implT::accumulate_timer,time,TaskAttributes::hipri());
            }

            return None;
        }


        /// same as do_apply_kernel, but use full rank tensors as input and low rank tensors as output

        /// @param[in]  op      the operator working on our function
        /// @param[in]  c       full rank tensor holding the NS coefficients
        /// @param[in]  args    laziness holding norm of the coefficients, displacement, destination, ..
        /// @param[in]  apply_targs TensorArgs with tightened threshold for accumulation
        /// @return     nothing, but accumulate the result tensor into the destination node
        template <typename opT, typename R, size_t OPDIM>
        Void do_apply_kernel2(const opT* op, const Tensor<R>& c, const do_op_args<OPDIM>& args,
                const TensorArgs& apply_targs) {

            const TensorArgs args_full(0.0,TT_FULL);

            tensorT result_full = op->apply(args.key, args.d, c, args.tol/args.fac/args.cnorm);
//            coeffT result=coeffT(result_full,args_full);
            coeffT result=coeffT(result_full,apply_targs);

            // Screen here to reduce communication cost of negligible data
            // and also to ensure we don't needlessly widen the tree when
            // applying the operator
            double norm;
            MADNESS_ASSERT(result.tensor_type()==TT_FULL or result.tensor_type()==TT_2D);
            if (result.tensor_type()==TT_2D) norm=result.config().svd_normf();
            if (result.tensor_type()==TT_FULL) norm=result.normf();
            if (norm > 0.1*args.tol/args.fac) {
                // OPTIMIZATION NEEDED HERE ... CHANGING THIS TO TASK NOT SEND REMOVED
                // BUILTIN OPTIMIZATION TO SHORTCIRCUIT MSG IF DATA IS LOCAL
                Future<double> time=coeffs.task(args.dest, &nodeT::accumulate, result, coeffs, args.dest, apply_targs,
                        TaskAttributes::hipri());
                woT::task(world.rank(),&implT::accumulate_timer,time,TaskAttributes::hipri());
            }
            return None;
        }


        /// same as do_apply_kernel2, but use low rank tensors as input and low rank tensors as output

        /// @param[in]  op      the operator working on our function
        /// @param[in]  c       full rank tensor holding the NS coefficients
        /// @param[in]  args    laziness holding norm of the coefficients, displacement, destination, ..
        /// @param[in]  apply_targs TensorArgs with tightened threshold for accumulation
        /// @return     nothing, but accumulate the result tensor into the destination node
        template <typename opT, typename R, size_t OPDIM>
        Void do_apply_kernel3(const opT* op, const GenTensor<R>& coeff, const do_op_args<OPDIM>& args,
                const TensorArgs& apply_targs) {

            coeffT result;
            if (2*OPDIM==NDIM) result= op->apply2_lowdim(args.key, args.d, coeff, args.tol/args.fac/args.cnorm, args.tol/args.fac);
            if (OPDIM==NDIM) result = op->apply2(args.key, args.d, coeff, args.tol/args.fac/args.cnorm, args.tol/args.fac);
            double result_norm=-1.0;
            if (result.tensor_type()==TT_2D) result_norm=result.config().svd_normf();
            if (result.tensor_type()==TT_FULL) result_norm=result.normf();
            MADNESS_ASSERT(result_norm>-0.5);

            if (result_norm> 0.3*args.tol/args.fac) {

                // accumulate also expects result in SVD form
                Future<double> time=coeffs.task(args.dest, &nodeT::accumulate, result, coeffs, args.dest, apply_targs,
                        TaskAttributes::hipri());
                woT::task(world.rank(),&implT::accumulate_timer,time,TaskAttributes::hipri());

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

            typedef typename opT::keyT opkeyT;
            static const size_t opdim=opT::opdim;

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
                    double opnorm = op->norm(key.level(), d, key);
                    // working assumption here is that the operator is isotropic and
                    // montonically decreasing with distance
                    double tol = truncate_tol(thresh, key);

                    //print("APP", key, dest, cnorm, opnorm, (cnorm*opnorm> tol/fac));

                    if (cnorm*opnorm> tol/fac) {

                        // Most expensive part is the kernel ... do it in a separate task
                        if (d.distsq()==0) {
                            // This introduces finer grain parallelism
                            ProcessID where = world.rank();
                            do_op_args<opdim> args(key, d, dest, tol, fac, cnorm);
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

            if (op.modified) flo_unary_op_node_inplace(do_stuff(),true);
            this->compressed=true;
            this->nonstandard=true;
            this->redundant=false;

        }


        /// apply an operator on the coeffs c (at node key)

        /// invoked by result; the result is accumulated inplace to this's tree at various FunctionNodes
        /// @param[in] op     the operator to act on the source function
        /// @param[in] key    key of the source FunctionNode of f which is processed (see "source")
        /// @param[in] coeff  coeffs of FunctionNode being processed
        template <typename opT, typename R>
        Void do_apply_source_driven(const opT* op, const keyT& key, const coeffT& coeff) {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            // insert timer here

            typedef typename opT::keyT opkeyT;
            static const size_t opdim=opT::opdim;
            Key<NDIM-opdim> nullkey(key.level());

            // source is that part of key that corresponds to those dimensions being processed
            opkeyT source;
            Key<NDIM-opdim> dummykey;
            if (op->particle()==1) key.break_apart(source,dummykey);
            if (op->particle()==2) key.break_apart(dummykey,source);

             // fac is the number of contributing neighbors (approx)
            const double tol = truncate_tol(thresh, key);

            double fac = 10.0; //3.0; // 10.0 seems good for qmprop ... 3.0 OK for others
            if (NDIM==6) fac=729; //100.0;
            if (op->modified()) fac*=10.0;
            double cnorm = coeff.normf();

            double wall0=wall_time();
            bool verbose=false;
            long neighbors=0;
            long generated_terms=0;


            // for accumulation: keep slightly tighter TensorArgs
            TensorArgs apply_targs(targs);
            apply_targs.thresh=tol/fac*0.01;
//            apply_targs.tt=TT_FULL;

            // for the kernel it may be more efficient to do the convolution in full rank
            tensorT coeff_full;


            const std::vector<opkeyT>& disp = op->get_disp(key.level());

            static const std::vector<bool> is_periodic(NDIM,false); // Periodic sum is already done when making rnlp
            double opnorm_old=100.0;
            double opnorm=100.0;

            for (typename std::vector<opkeyT>::const_iterator it=disp.begin(); it != disp.end(); ++it) {
                const opkeyT& d = *it;

                keyT disp1;
                if (op->particle()==1) disp1=it->merge_with(nullkey);
                if (op->particle()==2) disp1=nullkey.merge_with(*it);

                keyT dest = neighbor(key, disp1, is_periodic);

                if (dest.is_valid()) {
                    opnorm_old=opnorm;
                    opnorm = op->norm(key.level(), d, source);
//                    MADNESS_ASSERT(opnorm_old+1.e-10>=opnorm);

                    // working assumption here is that the operator is isotropic and
                    // montonically decreasing with distance

                    //print("APP", key, dest, cnorm, opnorm, (cnorm*opnorm> tol/fac));

                    if (cnorm*opnorm> tol/fac) {
//                    if (d.distsq()<2) {
//                        double wall00=wall_time();
                        neighbors++;

                        double cost_ratio=op->estimate_costs(source, d, coeff, tol/fac/cnorm, tol/fac);
//                        cost_ratio=1.5;     // force low rank
//                        cost_ratio=0.5;     // force full rank

                        if (cost_ratio>0.0) {

                            ProcessID here = world.rank();
//                            ProcessID there =  world.random_proc();
//                            ProcessID where = FunctionDefaults<NDIM>::get_apply_randomize() ? there : here;
                            do_op_args<opdim> args(source, d, dest, tol, fac, cnorm);

                            // Most expensive part is the kernel ... do it in a separate task -- full rank
                            if (cost_ratio<1.0) {

                                if (not coeff_full.has_data()) coeff_full=coeff.full_tensor_copy();

                                woT::task(here, &implT:: template do_apply_kernel2<opT,R,opdim>, op, coeff_full,
                                        args,apply_targs,TaskAttributes::hipri());

//                                world.taskq.add(*this,&implT:: template do_apply_kernel2<opT,R>, op, coeff_full,
//                                        args,apply_targs,TaskAttributes::hipri());

                            } else {
                                woT::task(here, &implT:: template do_apply_kernel3<opT,R,opdim> ,op,coeff,
                                        args,apply_targs,TaskAttributes::hipri());
                            }
                        }
                    } else if (d.distsq() >= 12) {
                        break; // Assumes monotonic decay beyond nearest neighbor
                    }
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

            bool do_target_driven=(NDIM==opT::opdim);
//            do_target_driven=false;
            if (do_target_driven) {

//                // make the MRA structure of the result tree
//                typename dcT::const_iterator end = f.get_coeffs().end();
//                for (typename dcT::const_iterator it=f.get_coeffs().begin(); it!=end; ++it) {
//                    const keyT& key = it->first;
//                    const coeffT& coeff = it->second.coeff();
//
//                    if (coeff.has_data() and (coeff.rank()!=0)) {
//                        const double cnorm=coeff.normf();
//                        ProcessID p = coeffs.owner(key);
//                        woT::task(p, &implT:: template do_apply_dry_run<opT>, &op, key, cnorm);
//                    }
//                }
//                world.gop.fence();
//
//                // consolidate the result tree
//                typename dcT::const_iterator end2 = get_coeffs().end();
//                for (typename dcT::const_iterator it=get_coeffs().begin(); it!=end2; ++it) {
//                    const keyT& key = it->first;
//                    coeffs.task(key, &FunctionNode<T,NDIM>::set_has_children_recursive, coeffs, key);
//                }
//                world.gop.fence();
//
//                // run the target_driven convolution ON ALL NODES
//                for (typename dcT::const_iterator it=get_coeffs().begin(); it!=end2; ++it) {
//                    const keyT& key = it->first;
//                    for (ProcessID p=0; p<world.size(); ++p) {
//                        woT::task(p, &implT:: template do_apply_target_driven<opT,R>, &op, &f, key);
//                    }
//                }


                const keyT key0=cdata.key0;
                const ProcessID p=coeffs.owner(key0);
                if (world.rank()==p) {
                    woT::task(p, &implT:: template do_apply_target_driven<opT,R>, &op, &f, key0);
                    std::vector<Future<bool> > has_contrib = future_vector_factory<bool>(world.size());
                    for (int i=0; i<world.size(); ++i) has_contrib[i].set(true);
                    woT::task(p, &implT:: template do_apply_target_driven_recursive<opT,R>,
//                            &op, &f, key0, has_contrib, TaskAttributes::hipri());
                            &op, &f, key0, TaskAttributes::hipri());
                }

            } else {


                // looping through all the coefficients of the source f
                typename dcT::const_iterator end = f.get_coeffs().end();
                for (typename dcT::const_iterator it=f.get_coeffs().begin(); it!=end; ++it) {

                    const keyT& key = it->first;
                    const coeffT& coeff = it->second.coeff();

                    if (coeff.has_data() and (coeff.rank()!=0)) {
                        ProcessID p = FunctionDefaults<NDIM>::get_apply_randomize() ? world.random_proc() : coeffs.owner(key);
                        woT::task(p, &implT:: template do_apply_source_driven<opT,R>, &op, key, coeff);
                    }
                }
            }
            world.gop.fence();

            // reduce the rank of the final nodes, leave full tensors unchanged
            flo_unary_op_node_inplace(do_reduce_rank(targs.thresh*0.01),true);
            flo_unary_op_node_inplace(do_reduce_rank(targs),true);

            // change TT_FULL to low rank
            flo_unary_op_node_inplace(do_change_tensor_type(targs),true);

            // truncate leaf nodes to avoid excessive tree refinement
            flo_unary_op_node_inplace(do_truncate_NS_leafs(this),true);

            if (op.modified()) {
                this->compressed=false;
                this->nonstandard=false;
                this->redundant=true;
            } else {
                this->compressed=true;
                this->nonstandard=true;
                this->redundant=false;
            }

            if (fence)
                world.gop.fence();
        }


        /// invoked by result; get the MRA structure

        /// @param[in] op     the operator to act on the source function
        /// @param[in] key    key of the source FunctionNode of f which is processed
        /// @param[in] cnorm  norm of the coefficients of key
        template <typename opT>
        typename enable_if_c<NDIM==opT::opdim, Void>::type
        do_apply_dry_run(const opT* op, const keyT& key, const double& cnorm) {

            const double tol = truncate_tol(thresh, key);
            // fac is the number of contributing neighbors (approx)
            double fac = 10.0; //3.0; // 10.0 seems good for qmprop ... 3.0 OK for others
            if (NDIM==6) fac=729; //100.0;

            // for accumulation: keep slightly tighter TensorArgs
            TensorArgs apply_targs(targs);
            apply_targs.thresh=tol/fac*0.01;

            const std::vector<keyT>& disp = op->get_disp(key.level());
            static const std::vector<bool> is_periodic(NDIM,false); // Periodic sum is already done when making rnlp

            for (typename std::vector<keyT>::const_iterator it=disp.begin(); it != disp.end(); ++it) {
                const keyT dest = neighbor(key, *it, is_periodic);
                if (dest.is_valid()) {
                    const double opnorm = op->norm(key.level(), *it, key);
                    if (cnorm*opnorm> tol/fac) {
                        coeffs.replace(dest,nodeT());
                    }
                }
            }
            return None;
        }

        template <typename opT>
        typename disable_if_c<NDIM==opT::opdim, Void>::type
        do_apply_dry_run(const opT* op, const keyT& key, const double& cnorm) {
            MADNESS_EXCEPTION("in dummy function do_apply_dry_run",1);
            return None;
        }


        /// submit a task to apply the operator on the children iff the parent had contributions

        /// invoked only on one node, but will spawn tasks on all nodes
        /// @param[in]  f   the source function
        /// @param[in]  key a FunctionNode that has already been processed; work on its children
        /// @param[in]  has_contrib for each process: has the node "key" had any contributions
        template <typename opT, typename R>
        Void do_apply_target_driven_recursive(opT* op, const FunctionImpl<R,NDIM>* f,
//                const keyT& key, const std::vector<Future<bool> >& has_contrib) const{
                const keyT& key) const{
            MADNESS_ASSERT(f->get_coeffs().probe(key));

            // fast return if there is no contribution from any node
            bool do_something=false;
//            for (size_t i=0; i<has_contrib.size(); ++i) {
//                if (has_contrib[i].get()) {
//                    do_something=true;
//                }
//            }

            do_something=f->get_coeffs().find(key).get()->second.has_children();
//            MADNESS_ASSERT(world.size()==1);
//            do_something=f->get_coeffs().probe(key);

            if (do_something) {
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT& child = kit.key();

                    ProcessID pp=f->get_coeffs().owner(child);
                    woT::task(pp, &implT:: template do_apply_target_driven_recursive<opT,R>,
//                            op, f, child, has_contrib2);
                            op, f, child, TaskAttributes::hipri());

//                    print("starting work on node ",child," at time",wall_time());
                    std::vector<Future<bool> > has_contrib2 = future_vector_factory<bool>(world.size());
                    for (ProcessID p=0; p<world.size(); ++p) {
                        has_contrib2[p]=woT::task(p, &implT:: template do_apply_target_driven<opT,R>,
                                op, f, child);
                    }
                }
            }
            return None;

        }

        /// apply an operator on all local boxes that contribute to target, then accumulate remotely

        /// invoked by result, for all possible FunctionNodes by all processes
        /// for now assume that the result function and the source function have the same tree structure
        /// @param[in]  op      operator acting on f
        /// @param[in]  f       the function op acts on, assuming SVD or full rank coefficients
        /// @param[in]  target  target FunctionNode (not necessarily local!)
        template <typename opT, typename R>
        typename enable_if_c<NDIM==opT::opdim, bool>::type
        do_apply_target_driven(opT* op, const FunctionImpl<R,NDIM>* f, const keyT& target) {

            double cpu0=cpu_time();

            // some prep stuff
            const double tol = truncate_tol(thresh, target);
            static const std::vector<bool> is_periodic(NDIM,false); // Periodic sum is already done when making rnlp
            // fac is the number of contributing neighbors (approx)
            double fac=(NDIM==6) ? 729.0 : 10.0; //100.0;
            if (op->modified()) fac*=10;

            // for accumulation: keep slightly tighter TensorArgs
            TensorArgs apply_targs(targs);
            apply_targs.thresh=tol/fac*0.01;

            // here we accumulate
            tensorT result_full(cdata.v2k);
            coeffT result_low;

            // loop over all local source FunctionNodes contributing to target
            const std::vector<keyT>& all_disp = op->get_disp(target.level());
            for (typename std::vector<keyT>::const_iterator it=all_disp.begin(); it != all_disp.end(); ++it) {

                // need to mirror displacement: source+disp=target
                const keyT mirror_disp(it->level(),it->translation()*(-1));
                const keyT disp(*it);
                const keyT source=neighbor(target, mirror_disp, is_periodic);

                if (source.is_valid() and f->get_coeffs().probe(source)) {
                    MADNESS_ASSERT(source.neighbor(disp)==target);

                    const nodeT& fnode=f->get_coeffs().find(source).get()->second;
                    const coeffT& coeff=fnode.coeff();  // this is local
                    tensorT coeff_full;
                    if (coeff.tensor_type()==TT_FULL) coeff_full=coeff.full_tensor();   // avoid deep copy below

//                    const double cnorm=fnode.normf();
                    const double cnorm=coeff.svd_normf();
//                    MADNESS_ASSERT(std::abs(cnorm-cnorm2)<1.e-14);

                    const double opnorm = op->norm(source.level(), disp, source);

                    if (cnorm*opnorm> tol/fac) {

                        double cost_ratio=op->estimate_costs(source, disp, coeff, tol/fac/cnorm, tol/fac);
//                            cost_ratio=1.5;     // force low rank
//                            cost_ratio=0.5;     // force full rank

                        if (cost_ratio>0.0) {

                            if (cost_ratio<1.0) {

                                if (not coeff_full.has_data()) coeff_full=coeff.full_tensor_copy();
                                result_full += op->apply(source,disp,coeff_full,tol/fac/cnorm);

                            } else {
                                const coeffT result=op->apply2(source,disp,coeff,tol/fac/cnorm,tol/fac);
                                result_low.add_SVD(result,apply_targs.thresh);
                            }
                        }
                    } else if (disp.distsq() >= 12) {
                        break; // Assumes monotonic decay beyond nearest neighbor
                    }
                }
            }

            // now that we have all the local contributions, send them out
            // accumulate also expects result in SVD form
            result_low.add_SVD(coeffT(result_full,apply_targs),apply_targs.thresh);
            Future<double> time=coeffs.task(target, &nodeT::accumulate, result_low, coeffs, target, apply_targs,
                    TaskAttributes::hipri());
            woT::task(world.rank(),&implT::accumulate_timer,time,TaskAttributes::hipri());
            double cpu1=cpu_time();
            timer_target_driven.accumulate(cpu1-cpu0);

            return (result_low.svd_normf()>tol/fac);
        }


        template <typename opT, typename R>
        typename disable_if_c<NDIM==opT::opdim, bool>::type
        do_apply_target_driven(opT* op, const FunctionImpl<R,NDIM>* f, const keyT key) {
            MADNESS_EXCEPTION("in dummy method do_apply_target_driven",1);
            return true;
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
            const tensorT coeff = node.coeff().full_tensor_copy();
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

        template<size_t LDIM>
        struct do_trace3functions_local {

            const implT* g;
            const implT* h;

            do_trace3functions_local() {}
            do_trace3functions_local(const implT* g, const implT* h) : g(g), h(h) {
                MADNESS_ASSERT(LDIM+LDIM==NDIM);
            }

            /// compute the trace f(1,2)g(1,3)h(2,3)

            /// @param[in]  it  iterator to a FunctionNode of f
            double operator()(typename dcT::const_iterator& it) const {

                Key<LDIM> fkey1,gkey1;
                Key<LDIM> fkey2,gkey2;
                it->first.break_apart(fkey1,fkey2);
                double result=0.0;

                // loop over g's keys
                typename dcT::const_iterator end = g->coeffs.end();
                for (typename dcT::const_iterator git=g->coeffs.begin(); git!=end; ++git) {

                    git->first.break_apart(gkey1,gkey2);

                    // find matching pairs for the trace f(1,2)g(1,3):
                    // first half of the key must be the same
                    if (fkey1==gkey1) {

                        // make the coeffs for h(2,3)
                        Key<NDIM> hkey=fkey2.merge_with(gkey2);
                        const coeffT hcoeff=h->get_functor()->coeff(hkey);
                        const coeffT& fcoeff=it->second.coeff();
                        const coeffT& gcoeff=git->second.coeff();
                        result+=inner3way(fcoeff,gcoeff,hcoeff);
                    }
                }
                return result;
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
            bool leaves_only=(this->is_redundant());
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
                            if (leaves_only) {
                                if (gnode.is_leaf() or fnode.is_leaf()) {
                                    sum += fnode.coeff().trace_conj(gnode.coeff());
                                }
                            } else {
                                sum += fnode.coeff().trace_conj(gnode.coeff());
                            }
                        }
                    }
                }
            }
            return sum;
        }

        /// Computes the trace of three functions of the form int dr1 dr2 dr3 f(1,2) g(1,3) h(2,3)

        /// requires all FunctionNodes of f and g with the same first half of their keys to be local!
        template <typename R>
        typename enable_if_c<NDIM==6, TENSOR_RESULT_TYPE(T,R)>::type
        trace3functions_local(const FunctionImpl<R,NDIM>& g,
                const FunctionImpl<T,NDIM>& h) const {

            typedef Range<typename dcT::const_iterator> rangeT;
            return world.taskq.reduce<double,rangeT,do_trace3functions_local<3> >(rangeT(coeffs.begin(),coeffs.end()),
                    do_trace3functions_local<3>(&g,&h));

        }

        /// project the low-dim function g on the hi-dim function f: this(x) = <f(x,y) | g(y)>

        /// invoked by result, a function of NDIM
        /// @param[in]  f   hi-dim function of LDIM+NDIM
        /// @param[in]  g   lo-dim function of LDIM
        /// @param[in]  dim over which dimensions to be integrated: 0..LDIM or LDIM..LDIM+NDIM-1
        /// @return this, with contributions on all scales
        template<size_t LDIM>
        void project_out(const FunctionImpl<T,LDIM+NDIM>* f, const FunctionImpl<T,LDIM>* g, const int dim) {

            typedef typename FunctionImpl<T,LDIM>::dcT::const_iterator giterator;
            typedef typename FunctionImpl<T,NDIM+LDIM>::dcT::const_iterator fiterator;

            // loop over all nodes of hi-dim f, compute the inner products with all
            // appropriate nodes of g, and accumulate in result
            fiterator end = f->get_coeffs().end();
            for (fiterator it=f->get_coeffs().begin(); it!=end; ++it) {
                const Key<LDIM+NDIM> key=it->first;
                const FunctionNode<T,LDIM+NDIM> fnode=it->second;
                const coeffT fcoeff=fnode.coeff();

                if (fnode.is_leaf() and fcoeff.has_data()) {

                    // break key into particle: over key1 will be summed, over key2 will be
                    // accumulated, or vice versa, depending on dim
                    Key<NDIM> key1;
                    Key<LDIM> key2;
                    key.break_apart(key1,key2);

                    if (dim==0) {
                        const Future<giterator> git=g->coeffs.find(key1);
                        woT::task(world.rank(),&implT::do_project_out<LDIM>,fcoeff,git,key2,dim);
                    } else if (dim==1) {
                        const Future<giterator> git=g->coeffs.find(key2);
                        woT::task(world.rank(),&implT::do_project_out<LDIM>,fcoeff,git,key1,dim);
                    } else {
                        MADNESS_EXCEPTION("confused dim in project_out",1);
                    }
                }
            }
            this->compressed=false;
            this->nonstandard=false;
            this->redundant=true;
        }

        /// compute the inner product of two nodes of only some dimensions and accumulate on result

        /// invoked by result
        /// @param[in]  fcoeff  coefficients of high dimension LDIM+NDIM
        /// @param[in]  git     iterator to node with coeffs of low dimension LDIM
        /// @param[in]  dest    destination node for the result
        /// @param[in]  dim     which dimensions should be contracted: 0..LDIM-1 or LDIM..NDIM+LDIM-1
        template<size_t LDIM>
        Void do_project_out(const coeffT& fcoeff, const typename FunctionImpl<T,LDIM>::dcT::const_iterator git,
                const Key<NDIM>& dest, const int dim) const {

            const coeffT gcoeff=git->second.coeff();

            // fast return if possible
            if (fcoeff.has_no_data() or gcoeff.has_no_data()) return None;

            // let's specialize for the time being on SVD tensors for f and full tensors of half dim for g
            MADNESS_ASSERT(gcoeff.tensor_type()==TT_FULL);
            MADNESS_ASSERT(fcoeff.tensor_type()==TT_2D);
            const tensorT gtensor=gcoeff.full_tensor();
            tensorT result(cdata.vk);

            const int otherdim=(dim+1)%2;

            // do the actual contraction
            for (int r=0; r<fcoeff.rank(); ++r) {
                const tensorT& contracted_tensor=fcoeff.config().ref_vector(dim);
                const tensorT& other_tensor=fcoeff.config().ref_vector(otherdim);
                const double ovlp= gtensor.trace_conj(contracted_tensor);
                const double fac=ovlp * fcoeff.config().weights(r);
                result+=fac*other_tensor;
            }

            // accumulate the result
            coeffs.task(dest, &nodeT::accumulate2, result, coeffs, dest, TaskAttributes::hipri());
            return None;
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

        /// print tree size and size
        void print_size(const std::string name) const {
            const size_t tsize=this->tree_size();
            const size_t size=this->size();
            const double wall=wall_time();
            if (this->world.rank()==0) {
                printf("%s at time %.1fs: treesize: %zu, size: %6.3f GByte\n",(name.c_str()),wall,
                        tsize,double(size)/(1024*1024*128));
            }
        }

        /// print the number of configurations per node
        void print_stats() const {
            if (this->targs.tt==TT_FULL) return;
            int dim=NDIM/2;
            int k0=k;
            if (is_compressed()) k0=2*k;
            Tensor<long> n(int(std::pow(k0,dim)+1));
            long n_full=0;
            long n_large=0;

            if (world.rank()==0) print("n.size(),k0,dim",n.size(),k0,dim);
            typename dcT::const_iterator end = coeffs.end();
            for (typename dcT::const_iterator it=coeffs.begin(); it!=end; ++it) {
                const nodeT& node = it->second;
                if (node.has_coeff()) {
                    if (node.coeff().rank()>long(n.size())) {
                        ++n_large;
                    } else if (node.coeff().rank()==-1) {
                        ++n_full;
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
                if (world.rank()==0) print("        full rank    ",n_full);
                for (unsigned int i=0; i<n.size(); i++) {
                    long m=n[i];
                    if (world.rank()==0) print("           ",i,"    ",m);
                }
                if (world.rank()==0) print("       large rank    ",n_large);
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
                    coeffs.replace(key,nodeT(node.coeff()*q,node.has_children()));
                }
                else {
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
                bool exists;
                ar & exists;
                if (exists) {
                    uniqueidT id;
                    ar & id;
                    World* world = World::world_from_id(id.get_world_id());
                    MADNESS_ASSERT(world);
                    ptr = static_cast< const FunctionImpl<T,NDIM>*>(world->ptr_from_id< WorldObject< FunctionImpl<T,NDIM> > >(id));
                    if (!ptr)
                        MADNESS_EXCEPTION("FunctionImpl: remote operation attempting to use a locally uninitialized object",0);
                } else {
                    ptr=NULL;
                }
            }
        };

        template <class Archive, class T, std::size_t NDIM>
        struct ArchiveStoreImpl<Archive,const FunctionImpl<T,NDIM>*> {
            static void store(const Archive& ar, const FunctionImpl<T,NDIM>*const& ptr) {
                bool exists=(ptr) ? true : false;
                ar & exists;
                if (exists) ar & ptr->id();
            }
        };

        template <class Archive, class T, std::size_t NDIM>
        struct ArchiveLoadImpl<Archive, FunctionImpl<T,NDIM>*> {
            static void load(const Archive& ar, FunctionImpl<T,NDIM>*& ptr) {
                bool exists;
                ar & exists;
                if (exists) {
                    uniqueidT id;
                    ar & id;
                    World* world = World::world_from_id(id.get_world_id());
                    MADNESS_ASSERT(world);
                    ptr = static_cast< FunctionImpl<T,NDIM>*>(world->ptr_from_id< WorldObject< FunctionImpl<T,NDIM> > >(id));
                    if (!ptr)
                        MADNESS_EXCEPTION("FunctionImpl: remote operation attempting to use a locally uninitialized object",0);
                } else {
                    ptr=NULL;
                }
            }
        };

        template <class Archive, class T, std::size_t NDIM>
        struct ArchiveStoreImpl<Archive, FunctionImpl<T,NDIM>*> {
            static void store(const Archive& ar, FunctionImpl<T,NDIM>*const& ptr) {
                bool exists=(ptr) ? true : false;
                ar & exists;
                if (exists) ar & ptr->id();
//                ar & ptr->id();
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
