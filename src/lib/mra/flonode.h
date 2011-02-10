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

#define HAVE_FLONODE 1

#ifndef MADNESS_MRA_FLONODE_H__INCLUDED
#define MADNESS_MRA_FLONODE_H__INCLUDED


/// \file flonode.h
/// \brief Provides a node/box for a MRA function, using LowRankTensor; to be replaced asap

#include <mra/sepreptensor.h>
#include <mra/funcdefaults.h>

namespace madness {


    /// FloNode holds the coefficients, etc., at each node of the 2^NDIM-tree
    template<typename T, int NDIM>
    class FloNode {
    private:
        // Should compile OK with these volatile but there should
        // be no need to set as volatile since the container internally
        // stores the entire entry as volatile

//        SepRepTensor<T>* _coeffs; ///< The coefficients, if any
        GenTensor<T> _coeffs; ///< The coefficients, if any
        double _norm_tree; ///< After norm_tree will contain norm of coefficients summed up tree
        bool _has_children; ///< True if there are children

    public:
        typedef WorldContainer<Key<NDIM> , FloNode<T, NDIM> > dcT; ///< Type of container holding the nodes

        /// Default constructor makes node without coeff or children
        FloNode() :
                _coeffs(0), _norm_tree(1e300), _has_children(false) {
        }

        /// Constructor from given coefficients with optional children

        /// Note that only a shallow copy of the coeff are taken so
        /// you should pass in a deep copy if you want the node to
        /// take ownership.

        /// construct a node with a full tensor as input;
        /// node will be in tensor product representation TPR
        explicit
        FloNode(const Tensor<T>& coeff, bool has_children = false) :
                _coeffs(0), _norm_tree(1e300), _has_children(has_children) {
        	_coeffs=new FullTensor<T>(coeff);
        }

        explicit
        FloNode(const Tensor<T>& coeff, double norm_tree, bool has_children) :
            _coeffs(0), _norm_tree(norm_tree), _has_children(has_children) {
        	_coeffs=new FullTensor<T>(coeff);
        }

        FloNode(const FloNode<T, NDIM>& other) {
            *this = other;
        }

        ~FloNode() {_coeffs.clear();};

        /// assignment will keep the data structure of the coeffs (TPR, SR)
        /// deep copy!
        FloNode<T, NDIM>&
        operator=(const FloNode<T, NDIM>& other) {
            if (this != &other) {
//                coeff() = copy(other.coeff());
            	_coeffs.clear();
            	// use virtual constructor
            	if (other.has_coeff()) _coeffs=copy(other.tensor());
                _norm_tree = other._norm_tree;
                _has_children = other._has_children;
            }
            return *this;
        }

        /// Copy with possible type conversion of coefficients, copying all other state

        /// Choose to not overload copy and type conversion operators
        /// so there are no automatic type conversions.
        template<typename Q>
        FloNode<Q, NDIM>
        convert() const {
        	MADNESS_EXCEPTION("no FloNode::convert",0);
//            return FloNode<Q, NDIM> (copy(coeff()), _has_children);
        }

        /// Returns true if there are coefficients in this node
        bool
        has_coeff() const {
        	return _coeffs.exists();
        };

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

        /// Returns the number of coefficients in this node
        size_t size() const {
        	return _coeffs.size();
        }

       /// Returns a const pointer to the SepRepTensor
        const GenTensor<T>
        tensor() const {
        	return _coeffs;
        }

        /// Returns an empty tensor if there are no coefficients.
        Tensor<T>
        full_tensor_copy() const {
        	Tensor<T> result;
        	if (this->has_coeff()) {
				if (_coeffs.type()==TT_2D) result=_coeffs.reconstruct_tensor();
				else if (_coeffs.type()==TT_3D) result=_coeffs.reconstruct_tensor();
				else if (_coeffs.type()==TT_FULL) result=_coeffs.full_tensor();
				else {
					MADNESS_EXCEPTION("unknown tensor type in full_tensor_copy",_coeffs.type());
				}
        	} else {
        		result=Tensor<T>();
        	}

        	return result;
        }

        /// Returns an empty tensor if there are no coefficients.
        Tensor<T>&
        full_tensor_reference() {
        	MADNESS_ASSERT(_coeffs.type()==TT_FULL);
        	return _coeffs.full_tensor();
        }

        /// Returns an empty tensor if there are no coefficients.
        const Tensor<T>&
        full_tensor_reference() const {
        	MADNESS_ASSERT(_coeffs->type()==TT_FULL);
        	return _coeffs->fullTensor();
        }

    private:
//        /// Returns an empty tensor if there are no coefficients.
//        Tensor<T>
//        coeff() {
////            MADNESS_ASSERT(_coeffs->ndim() == -1 || (_coeffs->dim(0) <= 2
////                                                    * MAXK && _coeffs->dim(0) >= 0));
////            return const_cast<Tensor<T>&>(_coeffs);
//            MADNESS_ASSERT(_coeffs==0 || (_coeffs->dim(0) <= 2
//                                                    * MAXK && _coeffs->dim(0) >= 0));
//            if (_coeffs==0) return Tensor<T> ();
//            return _coeffs->fullTensor();
//        }
//
//        /// Returns an empty tensor if there are no coefficients.
//        const Tensor<T>&
//        coeff() const {
////            return const_cast<const Tensor<T>&>(_coeffs);
//            MADNESS_ASSERT(_coeffs==0 || (_coeffs->dim(0) <= 2
//                                                    * MAXK && _coeffs->dim(0) >= 0));
//            if (_coeffs==0) return Tensor<T> ();
//            return _coeffs->fullTensor();
//        }

    public:


        /// Transforms the coeffs to the SR
        FloNode<T, NDIM>&
        node_to_low_rank(const double& thresh=FunctionDefaults<NDIM>::get_thresh()) {
        	MADNESS_ASSERT(thresh>0.0);
//        	static int counter=0;
        	to_low_rank(_coeffs,thresh,FunctionDefaults<NDIM>::get_tensor_type());
        	print("eps, new rank", thresh, this->_coeffs.rank(),this->_coeffs.size());
//        	print("ftr2sr calls",counter++,_coeffs->rank(),_coeffs->dim(0));
        	return *this;
        }

        /// Transforms the coeffs to the SR
        FloNode<T, NDIM>&
        node_to_full_rank() {
        	to_full_rank(_coeffs);
        	return *this;
        }



        /// Transforms the coeffs to the SR
        FloNode<T, NDIM>&
        ftr2sr(const double& thresh) {
        	MADNESS_ASSERT(thresh>0.0);
//        	static int counter=0;
        	to_low_rank(_coeffs,thresh,FunctionDefaults<NDIM>::get_tensor_type());
        	print("eps, new rank", thresh, this->_coeffs.rank(),this->_coeffs.size());
//        	print("ftr2sr calls",counter++,_coeffs->rank(),_coeffs->dim(0));
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
        set_has_children_recursive(const typename FloNode<T,NDIM>::dcT& c,const Key<NDIM>& key) {
            //madness::print("   set_chi_recu: ", key, *this);
            PROFILE_MEMBER_FUNC(FloNode);
            if (!(has_children() || has_coeff() || key.level()==0)) {
                // If node already knows it has children or it has
                // coefficients then it must already be connected to
                // its parent.  If not, the node was probably just
                // created for this operation and must be connected to
                // its parent.
                Key<NDIM> parent = key.parent();
                const_cast<dcT&>(c).task(parent, &FloNode<T,NDIM>::set_has_children_recursive, c, parent, TaskAttributes::hipri());
                //madness::print("   set_chi_recu: forwarding",key,parent);
            }
            _has_children = true;
            return None;
        }

        /// Sets \c has_children attribute to value of \c !flag
        void set_is_leaf(bool flag) {
            _has_children = !flag;
        }

        /// temporary template specializtion
        void set_coeff(const Tensor<T>& coeffs) {
        	_coeffs.clear();
        	_coeffs = new FullTensor<T>(coeffs);
        	to_low_rank(_coeffs,FunctionDefaults<NDIM>::get_thresh(),
        			FunctionDefaults<NDIM>::get_tensor_type());
            if ((_coeffs.dim(0) < 0) || (_coeffs.dim(0)>2*MAXK)) {
                print("set_coeff: may have a problem");
                print("set_coeff: coeff.dim[0] =", coeffs.dim(0), ", 2* MAXK =", 2*MAXK);
                print("set_coeff: coeff.dim[0] =", _coeffs.dim(0), ", 2* MAXK =", 2*MAXK);
                MADNESS_ASSERT(0);
            }
            MADNESS_ASSERT(coeffs.dim(0)<=2*MAXK && coeffs.dim(0)>=0);

        }

        /// Takes a \em shallow copy of the coeff --- same as \c this->coeff()=coeff
        void set_coeff(const GenTensor<T>& coeffs) {
        	_coeffs.clear();
        	_coeffs=coeffs;
//        	_coeffs = new FullTensor<T>(coeffs);
//        	to_low_rank(_coeffs,FunctionDefaults<NDIM>::get_thresh(),
//        			FunctionDefaults<NDIM>::get_tensor_type());
            if ((_coeffs.dim(0) < 0) || (_coeffs.dim(0)>2*MAXK)) {
                print("set_coeff: may have a problem");
                print("set_coeff: coeff.dim[0] =", coeffs.dim(0), ", 2* MAXK =", 2*MAXK);
                print("set_coeff: coeff.dim[0] =", _coeffs.dim(0), ", 2* MAXK =", 2*MAXK);
                MADNESS_ASSERT(0);
            }
            MADNESS_ASSERT(coeffs.dim(0)<=2*MAXK && coeffs.dim(0)>=0);
        }

        /// Clears the coefficients (has_coeff() will subsequently return false)
        void clear_coeff() {
        	_coeffs.clear();
        }

        /// Scale the coefficients of this node
        template <typename Q>
        void scale(Q a) {
        	_coeffs->scale(a);
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
        Void gaxpy_inplace(const T& alpha, const FloNode<Q,NDIM>& other, const R& beta) {
            PROFILE_MEMBER_FUNC(FuncNode);
            if (other.has_children())
                _has_children = true;
            if (has_coeff()) {
                if (other.has_coeff()) {
//                    coeff().gaxpy(alpha,other.coeff(),beta);
                	this->node_to_full_rank();
                    this->full_tensor_reference().gaxpy(alpha,other.full_tensor_copy(),beta);
                	this->node_to_low_rank();
                }
                else {
                    _coeffs.scale(alpha);
                }
            }
            else if (other.has_coeff()) {
            	this->node_to_full_rank();
            	this->full_tensor_reference() = other.full_tensor_copy()*beta; //? Is this the correct type conversion?
            	this->node_to_low_rank();
            }
            return None;
        }

        /// Accumulate inplace and if necessary connect node to parent
        Void accumulate(const Tensor<T>& t, const typename FloNode<T,NDIM>::dcT& c, const Key<NDIM>& key) {
            if (has_coeff()) {
            	this->node_to_full_rank();
                this->full_tensor_reference() += t;
                this->node_to_low_rank();
            }
            else {
                // No coeff and no children means the node is newly
                // created for this operation and therefore we must
                // tell its parent that it exists.
            	this->node_to_full_rank();
                this->full_tensor_reference() = copy(t);
                this->node_to_low_rank();

                if ((!_has_children) && key.level()> 0) {
                    Key<NDIM> parent = key.parent();
                    const_cast<dcT&>(c).task(parent, &FloNode<T,NDIM>::set_has_children_recursive, c, parent);
                }
            }
            return None;
        }

        T trace_conj(const FloNode<T,NDIM>& rhs) const {
        	return this->_coeffs.trace_conj((rhs._coeffs));
        }

        template <typename Archive>
        void serialize(Archive& ar) {

//            ar & coeff() & _has_children & double(_norm_tree);
        	Tensor<T> coeffs=this->full_tensor_copy();
            ar & coeffs & _has_children & double(_norm_tree);
        }
    };

    template <typename T, int NDIM>
    std::ostream& operator<<(std::ostream& s, const FloNode<T,NDIM>& node) {
        s << "(has_coeff=" << node.has_coeff() << ", has_children=" << node.has_children() << ", norm=";
        double norm = node.has_coeff() ? node.tensor().normf() : 0.0;
        if (norm < 1e-12)
            norm = 0.0;
        s << norm << ")";
        return s;
    }




}

#endif // MADNESS_MRA_FLONODE_H__INCLUDED
