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

#ifndef MADNESS_MRA_FUNCIMPL_H__INCLUDED
#define MADNESS_MRA_FUNCIMPL_H__INCLUDED

/// \file funcimpl.h
/// \brief Provides FunctionCommonData, FunctionImpl and FunctionFactory

#include <madness/world/MADworld.h>
#include <madness/world/print.h>
#include <madness/misc/misc.h>
#include <madness/tensor/tensor.h>
#include <madness/tensor/gentensor.h>

#include <madness/mra/function_common_data.h>
#include <madness/mra/indexit.h>
#include <madness/mra/key.h>
#include <madness/mra/funcdefaults.h>
#include <madness/mra/function_factory.h>
#include <madness/mra/displacements.h>

#include <madness/mra/leafop.h>

#include <array>
#include <iostream>
#include <type_traits>

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

}

namespace madness {


    /// A simple process map
    template<typename keyT>
    class SimplePmap : public WorldDCPmapInterface<keyT> {
    private:
        const int nproc;
        const ProcessID me;

    public:
        SimplePmap(World& world) : nproc(world.nproc()), me(world.rank())
        { }

        ProcessID owner(const keyT& key) const {
            if (key.level() == 0)
                return 0;
            else
                return key.hash() % nproc;
        }
    };

    /// A pmap that locates children on odd levels with their even level parents
    template <typename keyT>
    class LevelPmap : public WorldDCPmapInterface<keyT> {
    private:
        const int nproc;
    public:
        LevelPmap() : nproc(0) {};

        LevelPmap(World& world) : nproc(world.nproc()) {}

        /// Find the owner of a given key
        ProcessID owner(const keyT& key) const {
            Level n = key.level();
            if (n == 0) return 0;
            hashT hash;
            if (n <= 3 || (n&0x1)) hash = key.hash();
            else hash = key.parent().hash();
            return hash%nproc;
        }
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
        coeffT buffer; ///< The coefficients, if any
        double dnorm=-1.0;	///< norm of the d coefficients, also defined if there are no d coefficients
        double snorm=-1.0;	///< norm of the s coefficients

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

        explicit
        FunctionNode(const coeffT& coeff, double norm_tree, double snorm, double dnorm, bool has_children) :
	  _coeffs(coeff), _norm_tree(norm_tree), _has_children(has_children), dnorm(dnorm), snorm(snorm) {
        }

        FunctionNode(const FunctionNode<T, NDIM>& other) :
            _coeffs(other._coeffs), _norm_tree(other._norm_tree), _has_children(other._has_children),
            dnorm(other.dnorm), snorm(other.snorm) {
        }

        FunctionNode<T, NDIM>&
        operator=(const FunctionNode<T, NDIM>& other) {
            if (this != &other) {
                coeff() = copy(other.coeff());
                _norm_tree = other._norm_tree;
                _has_children = other._has_children;
                dnorm=other.dnorm;
                snorm=other.snorm;
                _norm_tree=other._norm_tree;
            }
            return *this;
        }

        /// Copy with possible type conversion of coefficients, copying all other state

        /// Choose to not overload copy and type conversion operators
        /// so there are no automatic type conversions.
        template<typename Q>
        FunctionNode<Q, NDIM>
        convert() const {
            return FunctionNode<Q, NDIM> (madness::convert<Q,T>(coeff()), _norm_tree, snorm, dnorm, _has_children);
        }

        /// Returns true if there are coefficients in this node
        bool
        has_coeff() const {
            return _coeffs.has_data();
        }


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

    public:

        /// reduces the rank of the coefficients (if applicable)
        void reduceRank(const double& eps) {
            _coeffs.reduce_rank(eps);
        }

        /// Sets \c has_children attribute to value of \c flag.
        void set_has_children(bool flag) {
            _has_children = flag;
        }

        /// Sets \c has_children attribute to true recurring up to ensure connected
        void set_has_children_recursive(const typename FunctionNode<T,NDIM>::dcT& c,const Key<NDIM>& key) {
            //madness::print("   set_chi_recu: ", key, *this);
            //PROFILE_MEMBER_FUNC(FunctionNode); // Too fine grain for routine profiling
            if (!(has_children() || has_coeff() || key.level()==0)) {
                // If node already knows it has children or it has
                // coefficients then it must already be connected to
                // its parent.  If not, the node was probably just
                // created for this operation and must be connected to
                // its parent.
                Key<NDIM> parent = key.parent();
                // Task on next line used to be TaskAttributes::hipri()) ... but deferring execution of this
                // makes sense since it is not urgent and lazy connection will likely mean that less forwarding
                // will happen since the upper level task will have already made the connection.
                const_cast<dcT&>(c).task(parent, &FunctionNode<T,NDIM>::set_has_children_recursive, c, parent);
                //const_cast<dcT&>(c).send(parent, &FunctionNode<T,NDIM>::set_has_children_recursive, c, parent);
                //madness::print("   set_chi_recu: forwarding",key,parent);
            }
            _has_children = true;
        }

        /// Sets \c has_children attribute to value of \c !flag
        void set_is_leaf(bool flag) {
            _has_children = !flag;
        }

        /// Takes a \em shallow copy of the coeff --- same as \c this->coeff()=coeff
        void set_coeff(const coeffT& coeffs) {
            coeff() = coeffs;
            if ((_coeffs.has_data()) and ((_coeffs.dim(0) < 0) || (_coeffs.dim(0)>2*MAXK))) {
                print("set_coeff: may have a problem");
                print("set_coeff: coeff.dim[0] =", coeffs.dim(0), ", 2* MAXK =", 2*MAXK);
            }
            MADNESS_ASSERT(coeffs.dim(0)<=2*MAXK && coeffs.dim(0)>=0);
        }

        /// Clears the coefficients (has_coeff() will subsequently return false)
        void clear_coeff() {
            coeff()=coeffT();
        }

        /// Scale the coefficients of this node
        template <typename Q>
        void scale(Q a) {
            _coeffs.scale(a);
        }

        /// Sets the value of norm_tree
        void set_norm_tree(double norm_tree) {
            _norm_tree = norm_tree;
        }

        /// Gets the value of norm_tree
        double get_norm_tree() const {
            return _norm_tree;
        }

        /// return the precomputed norm of the (virtual) d coefficients
        double get_dnorm() const {
        	return dnorm;
        }

        /// set the precomputed norm of the (virtual) s coefficients
        void set_snorm(const double sn) {
            snorm=sn;
        }

        /// set the precomputed norm of the (virtual) d coefficients
        void set_dnorm(const double dn) {
        	dnorm=dn;
        }

        /// get the precomputed norm of the (virtual) s coefficients
        double get_snorm() const {
            return snorm;
        }

        void recompute_snorm_and_dnorm(const FunctionCommonData<T,NDIM>& cdata) {
            snorm = 0.0;
            dnorm = 0.0;
            if (coeff().size() == 0) { ;
            } else if (coeff().dim(0) == cdata.vk[0]) {
                snorm = coeff().normf();

            } else if (coeff().is_full_tensor()) {
                Tensor<T> c = copy(coeff().get_tensor());
                snorm = c(cdata.s0).normf();
                c(cdata.s0) = 0.0;
                dnorm = c.normf();

            } else if (coeff().is_svd_tensor()) {
                coeffT c= coeff()(cdata.s0);
                snorm = c.normf();
                double norm = coeff().normf();
                dnorm = sqrt(norm * norm - snorm * snorm);

            } else {
                MADNESS_EXCEPTION("cannot use compute_dnorm", 1);
            }
        }


        /// General bi-linear operation --- this = this*alpha + other*beta

        /// This/other may not have coefficients.  Has_children will be
        /// true in the result if either this/other have children.
        template <typename Q, typename R>
        void gaxpy_inplace(const T& alpha, const FunctionNode<Q,NDIM>& other, const R& beta) {
            //PROFILE_MEMBER_FUNC(FuncNode);  // Too fine grain for routine profiling
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
        }

        /// Accumulate inplace and if necessary connect node to parent
        void accumulate2(const tensorT& t, const typename FunctionNode<T,NDIM>::dcT& c,
                           const Key<NDIM>& key) {
	  // double cpu0=cpu_time();
            if (has_coeff()) {
            	MADNESS_ASSERT(coeff().is_full_tensor());
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
		    if (c.is_local(parent))
			const_cast<dcT&>(c).send(parent, &FunctionNode<T,NDIM>::set_has_children_recursive, c, parent);
		    else
		      const_cast<dcT&>(c).task(parent, &FunctionNode<T,NDIM>::set_has_children_recursive, c, parent);
                }
            }
            //double cpu1=cpu_time();
        }


        /// Accumulate inplace and if necessary connect node to parent
        void accumulate(const coeffT& t, const typename FunctionNode<T,NDIM>::dcT& c,
                          const Key<NDIM>& key, const TensorArgs& args) {
            if (has_coeff()) {
                coeff().add_SVD(t,args.thresh);
                if (buffer.rank()<coeff().rank()) {
                    if (buffer.has_data()) {
                        buffer.add_SVD(coeff(),args.thresh);
                    } else {
                        buffer=copy(coeff());
                    }
                    coeff()=coeffT();
                }

            } else {
                // No coeff and no children means the node is newly
                // created for this operation and therefore we must
                // tell its parent that it exists.
            	coeff() = copy(t);
                if ((!_has_children) && key.level()> 0) {
                    Key<NDIM> parent = key.parent();
		    if (c.is_local(parent)) 
		      const_cast<dcT&>(c).send(parent, &FunctionNode<T,NDIM>::set_has_children_recursive, c, parent);
		    else
		      const_cast<dcT&>(c).task(parent, &FunctionNode<T,NDIM>::set_has_children_recursive, c, parent);
                }
            }
        }

        void consolidate_buffer(const TensorArgs& args) {
            if ((coeff().has_data()) and (buffer.has_data())) {
                coeff().add_SVD(buffer,args.thresh);
            } else if (buffer.has_data()) {
                coeff()=buffer;
            }
            buffer=coeffT();
        }

        T trace_conj(const FunctionNode<T,NDIM>& rhs) const {
            return this->_coeffs.trace_conj((rhs._coeffs));
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar & coeff() & _has_children & _norm_tree & dnorm & snorm;
        }

        /// like operator<<(ostream&, const FunctionNode<T,NDIM>&) but
        /// produces a sequence JSON-formatted key-value pairs
        /// @warning enclose the output in curly braces to make
        /// a valid JSON object
        void print_json(std::ostream& s) const {
            s << "\"has_coeff\":" << this->has_coeff()
              << ",\"has_children\":" << this->has_children() << ",\"norm\":";
            double norm = this->has_coeff() ? this->coeff().normf() : 0.0;
            if (norm < 1e-12)
                norm = 0.0;
            double nt = this->get_norm_tree();
            if (nt == 1e300)
                nt = 0.0;
            s << norm << ",\"norm_tree\":" << nt << ",\"snorm\":"
              << this->get_snorm() << ",\"dnorm\":" << this->get_dnorm()
              << ",\"rank\":" << this->coeff().rank();
            if (this->coeff().is_assigned())
                s << ",\"dim\":" << this->coeff().dim(0);
        }

    };

    template <typename T, std::size_t NDIM>
    std::ostream& operator<<(std::ostream& s, const FunctionNode<T,NDIM>& node) {
        s << "(has_coeff=" << node.has_coeff() << ", has_children=" << node.has_children() << ", norm=";
        double norm = node.has_coeff() ? node.coeff().normf() : 0.0;
        if (norm < 1e-12)
            norm = 0.0;
        double nt = node.get_norm_tree();
        if (nt == 1e300) nt = 0.0;
        s << norm << ", norm_tree, s/dnorm =" << nt << ", " << node.get_snorm() << " " << node.get_dnorm() << "), rank="<< node.coeff().rank()<<")";
        if (node.coeff().is_assigned()) s << " dim " << node.coeff().dim(0) << " ";
        return s;
    }


    /// returns true if the result of a hartree_product is a leaf node (compute norm & error)
    template<typename T, size_t NDIM>
    struct hartree_leaf_op {

        typedef FunctionImpl<T,NDIM> implT;
        const FunctionImpl<T,NDIM>* f;
        long k;
        bool do_error_leaf_op() const {return false;}

        hartree_leaf_op() = default;
        hartree_leaf_op(const implT* f, const long& k) : f(f), k(k) {}

        /// no pre-determination
        bool operator()(const Key<NDIM>& key) const {return false;}

        /// no post-determination
        bool operator()(const Key<NDIM>& key, const GenTensor<T>& coeff) const {
            MADNESS_EXCEPTION("no post-determination in hartree_leaf_op",1);
            return true;
        }

        /// post-determination: true if f is a leaf and the result is well-represented

        /// @param[in]  key the hi-dimensional key (breaks into keys for f and g)
        /// @param[in]  fcoeff coefficients of f of its appropriate key in NS form
        /// @param[in]  gcoeff coefficients of g of its appropriate key in NS form
        bool operator()(const Key<NDIM>& key, const Tensor<T>& fcoeff, const Tensor<T>& gcoeff) const {

            if (key.level()<2) return false;
            Slice s = Slice(0,k-1);
            std::vector<Slice> s0(NDIM/2,s);

            const double tol=f->get_thresh();
            const double thresh=f->truncate_tol(tol, key)*0.3;      // custom factor to "ensure" accuracy
            // include the wavelets in the norm, makes it much more accurate
            const double fnorm=fcoeff.normf();
            const double gnorm=gcoeff.normf();

            // if the final norm is small, perform the hartree product and return
            const double norm=fnorm*gnorm;  // computing the outer product
            if (norm < thresh) return true;

            // norm of the scaling function coefficients
            const double sfnorm=fcoeff(s0).normf();
            const double sgnorm=gcoeff(s0).normf();

            // get the error of both functions and of the pair function;
            // need the abs for numerics: sfnorm might be equal fnorm.
            const double ferror=sqrt(std::abs(fnorm*fnorm-sfnorm*sfnorm));
            const double gerror=sqrt(std::abs(gnorm*gnorm-sgnorm*sgnorm));

            // if the expected error is small, perform the hartree product and return
            const double error=fnorm*gerror + ferror*gnorm + ferror*gerror;
            //            const double error=sqrt(fnorm*fnorm*gnorm*gnorm - sfnorm*sfnorm*sgnorm*sgnorm);

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
        bool do_error_leaf_op() const {return true;}

        op_leaf_op() = default;
        op_leaf_op(const opT* op, const implT* f) : op(op), f(f) {}

        /// pre-determination: we can't know if this will be a leaf node before we got the final coeffs
        bool operator()(const Key<NDIM>& key) const {return true;}

        /// post-determination: return true if operator and coefficient norms are small
        bool operator()(const Key<NDIM>& key, const GenTensor<T>& coeff) const {
            if (key.level()<2) return false;
            const double cnorm=coeff.normf();
            return this->operator()(key,cnorm);
        }

        /// post-determination: return true if operator and coefficient norms are small
        bool operator()(const Key<NDIM>& key, const double& cnorm) const {
            if (key.level()<2) return false;

            typedef Key<opT::opdim> opkeyT;
            const opkeyT source=op->get_source_key(key);

            const double thresh=f->truncate_tol(f->get_thresh(),key);
            const std::vector<opkeyT>& disp = op->get_disp(key.level());
            const opkeyT& d = *disp.begin();         // use the zero-displacement for screening
            const double opnorm = op->norm(key.level(), d, source);
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
        bool do_error_leaf_op() const {return false;}

        hartree_convolute_leaf_op() = default;
        hartree_convolute_leaf_op(const implT* f, const implL* g, const opT* op)
            : f(f), g(g), op(op) {}

        /// no pre-determination
        bool operator()(const Key<NDIM>& key) const {return true;}

        /// no post-determination
        bool operator()(const Key<NDIM>& key, const GenTensor<T>& coeff) const {
            MADNESS_EXCEPTION("no post-determination in hartree_convolute_leaf_op",1);
            return true;
        }

        /// post-determination: true if f is a leaf and the result is well-represented

        /// @param[in]  key the hi-dimensional key (breaks into keys for f and g)
        /// @param[in]  fcoeff coefficients of f of its appropriate key in NS form
        /// @param[in]  gcoeff coefficients of g of its appropriate key in NS form
        bool operator()(const Key<NDIM>& key, const Tensor<T>& fcoeff, const Tensor<T>& gcoeff) const {
            //        bool operator()(const Key<NDIM>& key, const GenTensor<T>& coeff) const {

            if (key.level()<2) return false;

            const double tol=f->get_thresh();
            const double thresh=f->truncate_tol(tol, key);
            // include the wavelets in the norm, makes it much more accurate
            const double fnorm=fcoeff.normf();
            const double gnorm=gcoeff.normf();

            // norm of the scaling function coefficients
            const double sfnorm=fcoeff(g->get_cdata().s0).normf();
            const double sgnorm=gcoeff(g->get_cdata().s0).normf();

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
    	void operator()(const Key<NDIM>& key, const GenTensor<T>& coeff, const bool& is_leaf) const {}
        bool operator()(const Key<NDIM>& key, const GenTensor<T>& fcoeff, const GenTensor<T>& gcoeff) const {
            MADNESS_EXCEPTION("in noop::operator()",1);
            return true;
        }
        template <typename Archive> void serialize (Archive& ar) {}

    };

    /// insert/replaces the coefficients into the function
    template<typename T, std::size_t NDIM>
    struct insert_op {
    	typedef FunctionImpl<T,NDIM> implT;
    	typedef Key<NDIM> keyT;
    	typedef GenTensor<T> coeffT;
    	typedef FunctionNode<T,NDIM> nodeT;

    	implT* impl;
    	insert_op() : impl() {}
    	insert_op(implT* f) : impl(f) {}
    	insert_op(const insert_op& other) : impl(other.impl) {}
    	void operator()(const keyT& key, const coeffT& coeff, const bool& is_leaf) const {
    		MADNESS_ASSERT(impl->get_coeffs().is_local(key));
            impl->get_coeffs().replace(key,nodeT(coeff,not is_leaf));
    	}
        template <typename Archive> void serialize (Archive& ar) {
            ar & impl;
        }

    };

    /// inserts/accumulates coefficients into impl's tree

    /// NOTE: will use buffer and will need consolidation after operation ended !!  NOTE !!
    template<typename T, std::size_t NDIM>
    struct accumulate_op {
        typedef GenTensor<T> coeffT;
        typedef FunctionNode<T,NDIM> nodeT;

        FunctionImpl<T,NDIM>* impl;
        accumulate_op() = default;
        accumulate_op(FunctionImpl<T,NDIM>* f) : impl(f) {}
        accumulate_op(const accumulate_op& other) = default;
        void operator()(const Key<NDIM>& key, const coeffT& coeff, const bool& is_leaf) const {
            if (coeff.has_data())
                impl->get_coeffs().task(key, &nodeT::accumulate, coeff, impl->get_coeffs(), key, impl->get_tensor_args());
        }
        template <typename Archive> void serialize (Archive& ar) {
            ar & impl;
        }

    };


template<size_t NDIM>
    struct true_op {

    	template<typename T>
        bool operator()(const Key<NDIM>& key, const T& t) const {return true;}

    	template<typename T, typename R>
        bool operator()(const Key<NDIM>& key, const T& t, const R& r) const {return true;}
        template <typename Archive> void serialize (Archive& ar) {}

    };

    /// shallow-copy, pared-down version of FunctionNode, for special purpose only
    template<typename T, std::size_t NDIM>
    struct ShallowNode {
        typedef GenTensor<T> coeffT;
        coeffT _coeffs;
        bool _has_children;
        double dnorm=-1.0;
        ShallowNode() : _coeffs(), _has_children(false) {}
        ShallowNode(const FunctionNode<T,NDIM>& node)
            : _coeffs(node.coeff()), _has_children(node.has_children()),
			  dnorm(node.get_dnorm()) {}
        ShallowNode(const ShallowNode<T,NDIM>& node)
            : _coeffs(node.coeff()), _has_children(node._has_children),
			  dnorm(node.dnorm) {}

        const coeffT& coeff() const {return _coeffs;}
        coeffT& coeff() {return _coeffs;}
        bool has_children() const {return _has_children;}
        bool is_leaf() const {return not _has_children;}
        template <typename Archive>
        void serialize(Archive& ar) {
            ar & coeff() & _has_children & dnorm;
        }
    };


    /// a class to track where relevant (parent) coeffs are

    /// E.g. if a 6D function is composed of two 3D functions their coefficients must be tracked.
    /// We might need coeffs from a box that does not exist, and to avoid searching for
    /// parents we track which are their required respective boxes.
    ///  - CoeffTracker will refer either to a requested key, if it exists, or to its
    ///    outermost parent.
    ///  - Children must be made in sequential order to be able to track correctly.
    ///
    /// Usage: 	1. make the child of a given CoeffTracker.
    ///			   If the parent CoeffTracker refers to a leaf node (flag is_leaf)
    ///            the child will refer to the same node. Otherwise it will refer
    ///            to the child node.
    ///			2. retrieve its coefficients (possible communication/ returns a Future).
    ///            Member variable key always refers to an existing node,
    ///            so we can fetch it. Once we have the node we can determine
    ///            if it has children which allows us to make a child (see 1. )
    template<typename T, size_t NDIM>
    class CoeffTracker {

    	typedef FunctionImpl<T,NDIM> implT;
    	typedef Key<NDIM> keyT;
    	typedef GenTensor<T> coeffT;
        typedef std::pair<Key<NDIM>,ShallowNode<T,NDIM> > datumT;
        enum LeafStatus {no, yes, unknown};

        /// the funcimpl that has the coeffs
    	const implT* impl;
    	/// the current key, which must exists in impl
    	keyT key_;
    	/// flag if key is a leaf node
    	LeafStatus is_leaf_;
    	/// the coefficients belonging to key
    	coeffT coeff_;
    	/// norm of d coefficients corresponding to key
    	double dnorm_=-1.0;

    public:

    	/// default ctor
    	CoeffTracker() : impl(), key_(), is_leaf_(unknown), coeff_() {}

    	/// the initial ctor making the root key
    	CoeffTracker(const implT* impl) : impl(impl), key_(), is_leaf_(no) {
            if (impl) key_=impl->get_cdata().key0;
    	}

    	/// ctor with a pair<keyT,nodeT>
    	explicit CoeffTracker(const CoeffTracker& other, const datumT& datum)
    		: impl(other.impl), key_(other.key_), coeff_(datum.second.coeff()),
    			dnorm_(datum.second.dnorm) {
            if (datum.second.is_leaf()) is_leaf_=yes;
            else is_leaf_=no;
    	}

    	/// copy ctor
    	CoeffTracker(const CoeffTracker& other) : impl(other.impl), key_(other.key_),
    			is_leaf_(other.is_leaf_), coeff_(other.coeff_), dnorm_(other.dnorm_) {};

    	/// const reference to impl
    	const implT* get_impl() const {return impl;}

    	/// const reference to the coeffs
    	const coeffT& coeff() const {return coeff_;}

    	/// const reference to the key
    	const keyT& key() const {return key_;}

    	/// return the coefficients belonging to the passed-in key

    	/// if key equals tracked key just return the coeffs, otherwise
    	/// make the child coefficients.
    	/// @param[in]	key		return coeffs corresponding to this key
    	/// @return 	coefficients belonging to key
    	coeffT coeff(const keyT& key) const {
            MADNESS_ASSERT(impl);
            if (impl->is_compressed() or impl->is_nonstandard() or
            		impl->is_nonstandard_with_leaves())
                return impl->parent_to_child_NS(key,key_,coeff_);
            return impl->parent_to_child(coeff_,key_,key);
    	}

    	/// return the s and dnorm belonging to the passed-in key
    	double dnorm(const keyT& key) const {
    		if (key==key_) return dnorm_;
    		MADNESS_ASSERT(key.is_child_of(key_));
    		return 0.0;
    	}

    	/// const reference to is_leaf flag
    	const LeafStatus& is_leaf() const {return is_leaf_;}

    	/// make a child of this, ignoring the coeffs
    	CoeffTracker make_child(const keyT& child) const {

            // fast return
            if ((not impl) or impl->is_on_demand()) return CoeffTracker(*this);

            // can't make a child without knowing if this is a leaf -- activate first
            MADNESS_ASSERT((is_leaf_==yes) or (is_leaf_==no));

            CoeffTracker result;
            if (impl) {
                result.impl=impl;
                if (is_leaf_==yes) result.key_=key_;
                if (is_leaf_==no) {
                    result.key_=child;
                    // check if child is direct descendent of this, but root node is special case
                    if (child.level()>0) MADNESS_ASSERT(result.key().level()==key().level()+1);
                }
                result.is_leaf_=unknown;
            }
            return result;
    	}

    	/// find the coefficients

    	/// this involves communication to a remote node
    	/// @return	a Future<CoeffTracker> with the coefficients that key refers to
    	Future<CoeffTracker> activate() const {

            // fast return
            if (not impl) return Future<CoeffTracker>(CoeffTracker());
            if (impl->is_on_demand()) return Future<CoeffTracker>(CoeffTracker(impl));

            // this will return a <keyT,nodeT> from a remote node
            ProcessID p=impl->get_coeffs().owner(key());
            Future<datumT> datum1=impl->task(p, &implT::find_datum,key_,TaskAttributes::hipri());

            // construct a new CoeffTracker locally
            return impl->world.taskq.add(*const_cast<CoeffTracker*> (this),
                                         &CoeffTracker::forward_ctor,*this,datum1);
    	}

    private:
    	/// taskq-compatible forwarding to the ctor
    	CoeffTracker forward_ctor(const CoeffTracker& other, const datumT& datum) const {
            return CoeffTracker(other,datum);
    	}

    public:
    	/// serialization
        template <typename Archive> void serialize(const Archive& ar) {
            int il=int(is_leaf_);
            ar & impl & key_ & il & coeff_ & dnorm_;
            is_leaf_=LeafStatus(il);
        }
    };

    template<typename T, std::size_t NDIM>
    std::ostream&
    operator<<(std::ostream& s, const CoeffTracker<T,NDIM>& ct) {
        s << ct.key() << ct.is_leaf() << " " << ct.get_impl();
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
        typedef T typeT;
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

        World& world;

        /// getter
        int get_initial_level()const{return initial_level;}
        int get_special_level()const{return special_level;}
        const std::vector<Vector<double,NDIM> >& get_special_points()const{return special_points;}

    private:
        int k; ///< Wavelet order
        double thresh; ///< Screening threshold
        int initial_level; ///< Initial level for refinement
        int special_level; ///< Minimium level for refinement on special points
        std::vector<Vector<double,NDIM> > special_points; ///< special points for further refinement (needed for composite functions or multiplication)
        int max_refine_level; ///< Do not refine below this level
        int truncate_mode; ///< 0=default=(|d|<thresh), 1=(|d|<thresh/2^n), 2=(|d|<thresh/4^n);
        bool autorefine; ///< If true, autorefine where appropriate
        bool truncate_on_project; ///< If true projection inserts at level n-1 not n
        TensorArgs targs; ///< type of tensor to be used in the FunctionNodes

        const FunctionCommonData<T,NDIM>& cdata;

        std::shared_ptr< FunctionFunctorInterface<T,NDIM> > functor;
        TreeState tree_state;

        dcT coeffs; ///< The coefficients

        // Disable the default copy constructor
        FunctionImpl(const FunctionImpl<T,NDIM>& p);

    public:
        Timer timer_accumulate;
        Timer timer_change_tensor_type;
        Timer timer_lr_result;
        Timer timer_filter;
        Timer timer_compress_svd;
        Timer timer_target_driven;
        AtomicInt small;
        AtomicInt large;

        /// Initialize function impl from data in factory
        FunctionImpl(const FunctionFactory<T,NDIM>& factory)
            : WorldObject<implT>(factory._world)
            , world(factory._world)
            , k(factory._k)
            , thresh(factory._thresh)
            , initial_level(factory._initial_level)
	    , special_level(factory._special_level)
	    , special_points(factory._special_points)
            , max_refine_level(factory._max_refine_level)
            , truncate_mode(factory._truncate_mode)
            , autorefine(factory._autorefine)
            , truncate_on_project(factory._truncate_on_project)
//		  , nonstandard(false)
            , targs(factory._thresh,FunctionDefaults<NDIM>::get_tensor_type())
            , cdata(FunctionCommonData<T,NDIM>::get(k))
            , functor(factory.get_functor())
//		  , on_demand(factory._is_on_demand)
//		  , compressed(factory._compressed)
//		  , redundant(false)
		  , tree_state(factory._tree_state)
            , coeffs(world,factory._pmap,false)
            //, bc(factory._bc)
        {
            // PROFILE_MEMBER_FUNC(FunctionImpl); // No need to profile this
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
                // set the union of the special points of functor and the ones explicitly given to FunctionFactory
                std::vector<coordT> functor_special_points=functor->special_points();
                if (!functor_special_points.empty()) special_points.insert(special_points.end(), functor_special_points.begin(), functor_special_points.end());
                // near special points refine as deeply as requested by the factory AND the functor
                special_level = std::max(special_level, functor->special_level());

                typename dcT::const_iterator end = coeffs.end();
                for (typename dcT::const_iterator it=coeffs.begin(); it!=end; ++it) {
                    if (it->second.is_leaf())
                        woT::task(coeffs.owner(it->first), &implT::project_refine_op, it->first, do_refine,
                                  special_points);
                }
            }
            else { // Set as if a zero function
                initial_level = 1;
                insert_zero_down_to_initial_level(keyT(0));
            }

            coeffs.process_pending();
            this->process_pending();
            if (factory._fence && (functor || !empty)) world.gop.fence();
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
                , special_level(other.special_level)
                , special_points(other.special_points)
                , max_refine_level(other.max_refine_level)
                , truncate_mode(other.truncate_mode)
                , autorefine(other.autorefine)
                , truncate_on_project(other.truncate_on_project)
                , targs(other.targs)
                , cdata(FunctionCommonData<T,NDIM>::get(k))
                , functor()
                , tree_state(other.tree_state)
                , coeffs(world, pmap ? pmap : other.coeffs.get_pmap())
        {
            if (dozero) {
                initial_level = 1;
                insert_zero_down_to_initial_level(cdata.key0);
                //world.gop.fence(); <<<<<<<<<<<<<<<<<<<<<<   needs a fence argument
            }
            coeffs.process_pending();
            this->process_pending();
        }

        virtual ~FunctionImpl() { }

        const std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >& get_pmap() const;

        void replicate(bool fence=true) {
        	coeffs.replicate(fence);
        }

        void distribute(std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > > newmap) const {
        	auto currentmap=coeffs.get_pmap();
        	currentmap->redistribute(world,newmap);
        }


        /// Copy coeffs from other into self
        template <typename Q>
        void copy_coeffs(const FunctionImpl<Q,NDIM>& other, bool fence) {
            typename FunctionImpl<Q,NDIM>::dcT::const_iterator end = other.coeffs.end();
            for (typename FunctionImpl<Q,NDIM>::dcT::const_iterator it=other.coeffs.begin();
                 it!=end; ++it) {
                const keyT& key = it->first;
                const typename FunctionImpl<Q,NDIM>::nodeT& node = it->second;
                coeffs.replace(key,node. template convert<T>());
            }
            if (fence)
                world.gop.fence();
        }

        /// perform inplace gaxpy: this = alpha*this + beta*other
        /// @param[in]	alpha	prefactor for this
        /// @param[in]	beta	prefactor for other
        /// @param[in]	g       the other function, reconstructed
        /// @return     *this = alpha*this + beta*other, in either reconstructed or redundant_after_merge state
        template<typename Q, typename R>
        void gaxpy_inplace_reconstructed(const T& alpha, const FunctionImpl<Q,NDIM>& g, const R& beta, const bool fence) {
            // merge g's tree into this' tree
            gaxpy_inplace(alpha,g,beta,fence);
            tree_state=redundant_after_merge;
            // this->merge_trees(beta,g,alpha,fence);
            // tree is now redundant_after_merge
            // sum down the sum coeffs into the leafs if possible to keep the state most clean
            if (fence) sum_down(fence);
        }

        /// merge the trees of this and other, while multiplying them with the alpha or beta, resp

        /// first step in an inplace gaxpy operation for reconstructed functions; assuming the same
        /// distribution for this and other

        /// on output, *this = alpha* *this + beta * other
        /// @param[in]	alpha	prefactor for this
        /// @param[in]	beta	prefactor for other
        /// @param[in]	other	the other function, reconstructed
        template<typename Q, typename R>
        void merge_trees(const T alpha, const FunctionImpl<Q,NDIM>& other, const R beta, const bool fence=true) {
            MADNESS_ASSERT(get_pmap() == other.get_pmap());
            this->set_tree_state(redundant_after_merge);
            other.flo_unary_op_node_inplace(do_merge_trees<Q,R>(alpha,beta,*this),fence);
        }

        /// merge the trees of this and other, while multiplying them with the alpha or beta, resp

        /// result and rhs do not have to have the same distribution or live in the same world
        /// result+=alpha* this
        /// @param[in]	alpha	prefactor for this
        template<typename Q, typename R>
        void accumulate_trees(FunctionImpl<Q,NDIM>& result, const R alpha, const bool fence=true) const {
            flo_unary_op_node_inplace(do_accumulate_trees<Q,R>(result,alpha),fence);
        }

        /// perform: this= alpha*f + beta*g, invoked by result

        /// f and g are reconstructed, so we can save on the compress operation,
        /// walk down the joint tree, and add leaf coefficients; effectively refines
        /// to common finest level.

        /// nothing returned, but leaves this's tree reconstructed and as sum of f and g
        /// @param[in]  alpha   prefactor for f
        /// @param[in]  f       first addend
        /// @param[in]  beta    prefactor for g
        /// @param[in]  g       second addend
        void gaxpy_oop_reconstructed(const double alpha, const implT& f,
                                     const double beta, const implT& g, const bool fence);

        /// functor for the gaxpy_inplace method
        template <typename Q, typename R>
        struct do_gaxpy_inplace {
            typedef Range<typename FunctionImpl<Q,NDIM>::dcT::const_iterator> rangeT;
            FunctionImpl<T,NDIM>* f; ///< prefactor for current function impl
            T alpha; ///< the current function impl
            R beta; ///< prefactor for other function impl
            do_gaxpy_inplace() = default;
            do_gaxpy_inplace(FunctionImpl<T,NDIM>* f, T alpha, R beta) : f(f), alpha(alpha), beta(beta) {}
            bool operator()(typename rangeT::iterator& it) const {
                const keyT& key = it->first;
                const FunctionNode<Q,NDIM>& other_node = it->second;
                // Use send to get write accessor and automated construction if missing
                f->coeffs.send(key, &nodeT:: template gaxpy_inplace<Q,R>, alpha, other_node, beta);
                return true;
            }
            template <typename Archive>
            void serialize(Archive& ar) {
                ar & f & alpha & beta;
            }
        };

        /// Inplace general bilinear operation

        /// this's world can differ from other's world
        /// this = alpha * this + beta * other
        /// @param[in]  alpha   prefactor for the current function impl
        /// @param[in]  other   the other function impl
        /// @param[in]  beta    prefactor for other
        template <typename Q, typename R>
        void gaxpy_inplace(const T& alpha,const FunctionImpl<Q,NDIM>& other, const R& beta, bool fence) {
//            MADNESS_ASSERT(get_pmap() == other.get_pmap());
            if (alpha != T(1.0)) scale_inplace(alpha,false);
            typedef Range<typename FunctionImpl<Q,NDIM>::dcT::const_iterator> rangeT;
            typedef do_gaxpy_inplace<Q,R> opT;
            other.world.taskq. template for_each<rangeT,opT>(rangeT(other.coeffs.begin(), other.coeffs.end()), opT(this, T(1.0), beta));
            if (fence)
                other.world.gop.fence();
        }

        // loads a function impl from persistence
        // @param[in] ar   the archive where the function impl is stored
        template <typename Archive>
        void load(Archive& ar) {
            // WE RELY ON K BEING STORED FIRST
            int kk = 0;
            ar & kk;

            MADNESS_ASSERT(kk==k);

            // note that functor should not be (re)stored
            ar & thresh & initial_level & max_refine_level & truncate_mode
                & autorefine & truncate_on_project & tree_state;//nonstandard & compressed ; //& bc;

            ar & coeffs;
            world.gop.fence();
        }

        // saves a function impl to persistence
        // @param[in] ar   the archive where the function impl is to be stored
        template <typename Archive>
        void store(Archive& ar) {
            // WE RELY ON K BEING STORED FIRST

            // note that functor should not be (re)stored
            ar & k & thresh & initial_level & max_refine_level & truncate_mode
                & autorefine & truncate_on_project & tree_state;//nonstandard & compressed ; //& bc;

            ar & coeffs;
            world.gop.fence();
        }

        /// Returns true if the function is compressed.
        bool is_compressed() const;

        /// Returns true if the function is compressed.
        bool is_reconstructed() const;

        /// Returns true if the function is redundant.
        bool is_redundant() const;

        /// Returns true if the function is redundant_after_merge.
        bool is_redundant_after_merge() const;

        bool is_nonstandard() const;

        bool is_nonstandard_with_leaves() const;

        bool is_on_demand() const;

        bool has_leaves() const;

        void set_tree_state(const TreeState& state) {
        	tree_state=state;
        }

        TreeState get_tree_state() const {return tree_state;}

        void set_functor(const std::shared_ptr<FunctionFunctorInterface<T,NDIM> > functor1);

        std::shared_ptr<FunctionFunctorInterface<T,NDIM> > get_functor();

        std::shared_ptr<FunctionFunctorInterface<T,NDIM> > get_functor() const;

        void unset_functor();


        TensorType get_tensor_type() const;

        TensorArgs get_tensor_args() const;
        void set_tensor_args(const TensorArgs& t);

        double get_thresh() const;

        void set_thresh(double value);

        bool get_autorefine() const;

        void set_autorefine(bool value);

        int get_k() const;

        const dcT& get_coeffs() const;

        dcT& get_coeffs();

        const FunctionCommonData<T,NDIM>& get_cdata() const;

        void accumulate_timer(const double time) const; // !!!!!!!!!!!!  REDUNDANT !!!!!!!!!!!!!!!

        void print_timer() const;

        void reset_timer();

        /// Adds a constant to the function.  Local operation, optional fence

        /// In scaling function basis must add value to first polyn in
        /// each box with appropriate scaling for level.  In wavelet basis
        /// need only add at level zero.
        /// @param[in]  t   the scalar to be added
        void add_scalar_inplace(T t, bool fence);

        /// Initialize nodes to zero function at initial_level of refinement.

        /// Works for either basis.  No communication.
        void insert_zero_down_to_initial_level(const keyT& key);

        /// Truncate according to the threshold with optional global fence

        /// If thresh<=0 the default value of this->thresh is used
        /// @param[in]  tol   the truncation tolerance
        void truncate(double tol, bool fence);

        /// Returns true if after truncation this node has coefficients

        /// Assumed to be invoked on process owning key.  Possible non-blocking
        /// communication.
        /// @param[in]  key   the key of the current function node
        Future<bool> truncate_spawn(const keyT& key, double tol);

        /// Actually do the truncate operation
        /// @param[in] key the key to the current function node being evaluated for truncation
        /// @param[in] tol the tolerance for thresholding
        /// @param[in] v vector of Future<bool>'s that specify whether the current nodes children have coeffs
        bool truncate_op(const keyT& key, double tol, const std::vector< Future<bool> >& v);

        /// Evaluate function at quadrature points in the specified box

        /// @param[in] key the key indicating where the quadrature points are located
        /// @param[in] f the interface to the elementary function
        /// @param[in] qx quadrature points on a level=0 box
        /// @param[out] fval values
        void fcube(const keyT& key, const FunctionFunctorInterface<T,NDIM>& f, const Tensor<double>& qx, tensorT& fval) const;

        /// Evaluate function at quadrature points in the specified box

        /// @param[in] key the key indicating where the quadrature points are located
        /// @param[in] f the interface to the elementary function
        /// @param[in] qx quadrature points on a level=0 box
        /// @param[out] fval values
        void fcube(const keyT& key,  T (*f)(const coordT&), const Tensor<double>& qx, tensorT& fval) const;

        /// Returns cdata.key0
        const keyT& key0() const;

        /// Prints the coeffs tree of the current function impl
        /// @param[in] maxlevel the maximum level of the tree for printing
        /// @param[out] os the ostream to where the output is sent
        void print_tree(std::ostream& os = std::cout, Level maxlevel = 10000) const;

        /// Functor for the do_print_tree method
        void do_print_tree(const keyT& key, std::ostream& os, Level maxlevel) const;

        /// Prints the coeffs tree of the current function impl (using GraphViz)
        /// @param[in] maxlevel the maximum level of the tree for printing
        /// @param[out] os the ostream to where the output is sent
        void print_tree_graphviz(std::ostream& os = std::cout, Level maxlevel = 10000) const;

        /// Functor for the do_print_tree method (using GraphViz)
        void do_print_tree_graphviz(const keyT& key, std::ostream& os, Level maxlevel) const;

        /// Same as print_tree() but in JSON format
        /// @param[out] os the ostream to where the output is sent
        /// @param[in] maxlevel the maximum level of the tree for printing
        void print_tree_json(std::ostream& os = std::cout, Level maxlevel = 10000) const;

        /// Functor for the do_print_tree_json method
        void do_print_tree_json(const keyT& key, std::multimap<Level, std::tuple<tranT, std::string>>& data, Level maxlevel) const;

        /// convert a number [0,limit] to a hue color code [blue,red],
        /// or, if log is set, a number [1.e-10,limit]
        struct do_convert_to_color {
            double limit;
            bool log;
            static double lower() {return 1.e-10;};
            do_convert_to_color() = default;
            do_convert_to_color(const double limit, const bool log) : limit(limit), log(log) {}
            double operator()(double val) const {
                double color=0.0;

                if (log) {
                    double val2=log10(val) - log10(lower());        // will yield >0.0
                    double upper=log10(limit) -log10(lower());
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
        void print_plane(const std::string filename, const int xaxis, const int yaxis, const coordT& el2);

        /// collect the data for a plot of the MRA structure locally on each node

        /// @param[in]	xaxis	the x-axis in the plot (can be any axis of the MRA box)
        /// @param[in]	yaxis	the y-axis in the plot (can be any axis of the MRA box)
        /// @param[in]	el2     needs a description
        /// \todo Provide a description for el2
        Tensor<double> print_plane_local(const int xaxis, const int yaxis, const coordT& el2);

        /// Functor for the print_plane method
        /// @param[in] filename the filename for the output
        /// @param[in] plotinfo plotting parameters
        /// @param[in]	xaxis	the x-axis in the plot (can be any axis of the MRA box)
        /// @param[in]	yaxis	the y-axis in the plot (can be any axis of the MRA box)
        void do_print_plane(const std::string filename, std::vector<Tensor<double> > plotinfo,
                            const int xaxis, const int yaxis, const coordT el2);

        /// print the grid (the roots of the quadrature of each leaf box)
        /// of this function in user xyz coordinates
        /// @param[in] filename the filename for the output
        void print_grid(const std::string filename) const;

        /// return the keys of the local leaf boxes
        std::vector<keyT> local_leaf_keys() const;

        /// print the grid in xyz format

        /// the quadrature points and the key information will be written to file,
        /// @param[in]	filename	where the quadrature points will be written to
        /// @param[in]	keys		all leaf keys
        void do_print_grid(const std::string filename, const std::vector<keyT>& keys) const;

        /// read data from a grid

        /// @param[in]	keyfile		file with keys and grid points for each key
        /// @param[in]	gridfile 	file with grid points, w/o key, but with same ordering
        /// @param[in]	vnuc_functor	subtract the values of this functor if regularization is needed
        template<size_t FDIM>
        typename std::enable_if<NDIM==FDIM>::type
        read_grid(const std::string keyfile, const std::string gridfile,
                  std::shared_ptr< FunctionFunctorInterface<double,NDIM> > vnuc_functor) {

            std::ifstream kfile(keyfile.c_str());
            std::ifstream gfile(gridfile.c_str());
            std::string line;

            long ndata,ndata1;
            if (not (std::getline(kfile,line))) MADNESS_EXCEPTION("failed reading 1st line of key data",0);
            if (not (std::istringstream(line) >> ndata)) MADNESS_EXCEPTION("failed reading k",0);
            if (not (std::getline(gfile,line))) MADNESS_EXCEPTION("failed reading 1st line of grid data",0);
            if (not (std::istringstream(line) >> ndata1)) MADNESS_EXCEPTION("failed reading k",0);
            MADNESS_CHECK(ndata==ndata1);
            if (not (std::getline(kfile,line))) MADNESS_EXCEPTION("failed reading 2nd line of key data",0);
            if (not (std::getline(gfile,line))) MADNESS_EXCEPTION("failed reading 2nd line of grid data",0);

            // the quadrature points in simulation coordinates of the root node
            const Tensor<double> qx=cdata.quad_x;
            const size_t npt = qx.dim(0);

            // the number of coordinates (grid point tuples) per box ({x1},{x2},{x3},..,{xNDIM})
            long npoints=power<NDIM>(npt);
            // the number of boxes
            long nboxes=ndata/npoints;
            MADNESS_ASSERT(nboxes*npoints==ndata);
            print("reading ",nboxes,"boxes from file",gridfile,keyfile);

            // these will be the data
            Tensor<T> values(cdata.vk,false);

            int ii=0;
            std::string gline,kline;
            //            while (1) {
            while (std::getline(kfile,kline)) {

                double x,y,z,x1,y1,z1,val;

                // get the key
                long nn;
                Translation l1,l2,l3;
                // line looks like: # key:      n      l1   l2   l3
                kline.erase(0,7);
                std::stringstream(kline) >>  nn >> l1 >> l2 >> l3;
                //				kfile >> s >>  nn >> l1 >> l2 >> l3;
                const Vector<Translation,3> ll{ l1,l2,l3 };
                Key<3> key(nn,ll);

                // this is borrowed from fcube
                const Vector<Translation,3>& l = key.translation();
                const Level n = key.level();
                const double h = std::pow(0.5,double(n));
                coordT c; // will hold the point in user coordinates
                const Tensor<double>& cell_width = FunctionDefaults<NDIM>::get_cell_width();
                const Tensor<double>& cell = FunctionDefaults<NDIM>::get_cell();


                if (NDIM == 3) {
                    for (size_t i=0; i<npt; ++i) {
                        c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                        for (size_t j=0; j<npt; ++j) {
                            c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                            for (size_t k=0; k<npt; ++k) {
                                c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                                //								fprintf(pFile,"%18.12f %18.12f %18.12f\n",c[0],c[1],c[2]);
                                auto& success1 = std::getline(gfile,gline); MADNESS_CHECK(success1);
                                auto& success2 = std::getline(kfile,kline); MADNESS_CHECK(success2);
                                std::istringstream(gline) >> x >> y >> z >> val;
                                std::istringstream(kline) >> x1 >> y1 >> z1;
                                MADNESS_CHECK(std::fabs(x-c[0])<1.e-4);
                                MADNESS_CHECK(std::fabs(x1-c[0])<1.e-4);
                                MADNESS_CHECK(std::fabs(y-c[1])<1.e-4);
                                MADNESS_CHECK(std::fabs(y1-c[1])<1.e-4);
                                MADNESS_CHECK(std::fabs(z-c[2])<1.e-4);
                                MADNESS_CHECK(std::fabs(z1-c[2])<1.e-4);

                                // regularize if a functor is given
                                if (vnuc_functor) val-=(*vnuc_functor)(c);
                                values(i,j,k)=val;
                            }
                        }
                    }
                } else {
                    MADNESS_EXCEPTION("only NDIM=3 in print_grid",0);
                }

                // insert the new leaf node
                const bool has_children=false;
                coeffT coeff=coeffT(this->values2coeffs(key,values),targs);
                nodeT node(coeff,has_children);
                coeffs.replace(key,node);
                const_cast<dcT&>(coeffs).send(key.parent(), &FunctionNode<T,NDIM>::set_has_children_recursive, coeffs, key.parent());
                ii++;
            }

            kfile.close();
            gfile.close();
            MADNESS_CHECK(ii==nboxes);

        }


        /// read data from a grid

        /// @param[in]	gridfile		file with keys and grid points and values for each key
        /// @param[in]	vnuc_functor	subtract the values of this functor if regularization is needed
        template<size_t FDIM>
        typename std::enable_if<NDIM==FDIM>::type
        read_grid2(const std::string gridfile,
                   std::shared_ptr< FunctionFunctorInterface<double,NDIM> > vnuc_functor) {

            std::ifstream gfile(gridfile.c_str());
            std::string line;

            long ndata;
            if (not (std::getline(gfile,line))) MADNESS_EXCEPTION("failed reading 1st line of grid data",0);
            if (not (std::istringstream(line) >> ndata)) MADNESS_EXCEPTION("failed reading k",0);
            if (not (std::getline(gfile,line))) MADNESS_EXCEPTION("failed reading 2nd line of grid data",0);

            // the quadrature points in simulation coordinates of the root node
            const Tensor<double> qx=cdata.quad_x;
            const size_t npt = qx.dim(0);

            // the number of coordinates (grid point tuples) per box ({x1},{x2},{x3},..,{xNDIM})
            long npoints=power<NDIM>(npt);
            // the number of boxes
            long nboxes=ndata/npoints;
            MADNESS_CHECK(nboxes*npoints==ndata);
            print("reading ",nboxes,"boxes from file",gridfile);

            // these will be the data
            Tensor<T> values(cdata.vk,false);

            int ii=0;
            std::string gline;
            //            while (1) {
            while (std::getline(gfile,gline)) {

                double x1,y1,z1,val;

                // get the key
                long nn;
                Translation l1,l2,l3;
                // line looks like: # key:      n      l1   l2   l3
                gline.erase(0,7);
                std::stringstream(gline) >>  nn >> l1 >> l2 >> l3;
                const Vector<Translation,3> ll{ l1,l2,l3 };
                Key<3> key(nn,ll);

                // this is borrowed from fcube
                const Vector<Translation,3>& l = key.translation();
                const Level n = key.level();
                const double h = std::pow(0.5,double(n));
                coordT c; // will hold the point in user coordinates
                const Tensor<double>& cell_width = FunctionDefaults<NDIM>::get_cell_width();
                const Tensor<double>& cell = FunctionDefaults<NDIM>::get_cell();


                if (NDIM == 3) {
                    for (int i=0; i<npt; ++i) {
                        c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                        for (int j=0; j<npt; ++j) {
                            c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                            for (int k=0; k<npt; ++k) {
                                c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z

                                auto& success = std::getline(gfile,gline);
                                MADNESS_CHECK(success);
                                std::istringstream(gline) >> x1 >> y1 >> z1 >> val;
                                MADNESS_CHECK(std::fabs(x1-c[0])<1.e-4);
                                MADNESS_CHECK(std::fabs(y1-c[1])<1.e-4);
                                MADNESS_CHECK(std::fabs(z1-c[2])<1.e-4);

                                // regularize if a functor is given
                                if (vnuc_functor) val-=(*vnuc_functor)(c);
                                values(i,j,k)=val;
                            }
                        }
                    }
                } else {
                    MADNESS_EXCEPTION("only NDIM=3 in print_grid",0);
                }

                // insert the new leaf node
                const bool has_children=false;
                coeffT coeff=coeffT(this->values2coeffs(key,values),targs);
                nodeT node(coeff,has_children);
                coeffs.replace(key,node);
                const_cast<dcT&>(coeffs).send(key.parent(),
                                              &FunctionNode<T,NDIM>::set_has_children_recursive,
                                              coeffs, key.parent());
                ii++;
            }

            gfile.close();
            MADNESS_CHECK(ii==nboxes);

        }


        /// Compute by projection the scaling function coeffs in specified box
        /// @param[in] key the key to the current function node (box)
        tensorT project(const keyT& key) const;

        /// Returns the truncation threshold according to truncate_method

        /// here is our handwaving argument:
        /// this threshold will give each FunctionNode an error of less than tol. The
        /// total error can then be as high as sqrt(#nodes) * tol. Therefore in order
        /// to account for higher dimensions: divide tol by about the root of number
        /// of siblings (2^NDIM) that have a large error when we refine along a deep
        /// branch of the tree.
        double truncate_tol(double tol, const keyT& key) const;


        /// Returns patch referring to coeffs of child in parent box
        /// @param[in] child the key to the child function node (box)
        std::vector<Slice> child_patch(const keyT& child) const;

        /// Projection with optional refinement w/ special points
        /// @param[in] key the key to the current function node (box)
        /// @param[in] do_refine should we continue refinement?
        /// @param[in] specialpts vector of special points in the function where we need
        ///            to refine at a much finer level
        void project_refine_op(const keyT& key, bool do_refine,
                               const std::vector<Vector<double,NDIM> >& specialpts);

        /// Compute the Legendre scaling functions for multiplication

        /// Evaluate parent polyn at quadrature points of a child.  The prefactor of
        /// 2^n/2 is included.  The tensor must be preallocated as phi(k,npt).
        /// Refer to the implementation notes for more info.
        /// @todo Robert please verify this comment. I don't understand this method.
        /// @param[in] np level of the parent function node (box)
        /// @param[in] nc level of the child function node (box)
        /// @param[in] lp translation of the parent function node (box)
        /// @param[in] lc translation of the child function node (box)
        /// @param[out] phi tensor of the legendre scaling functions
        void phi_for_mul(Level np, Translation lp, Level nc, Translation lc, Tensor<double>& phi) const;

        /// Directly project parent coeffs to child coeffs

        /// Currently used by diff, but other uses can be anticipated

        /// @todo is this documentation correct?
        /// @param[in]	child	the key whose coeffs we are requesting
        /// @param[in]	parent	the (leaf) key of our function
        /// @param[in]	s	the (leaf) coeffs belonging to parent
        /// @return 	coeffs
        const coeffT parent_to_child(const coeffT& s, const keyT& parent, const keyT& child) const;

        /// Directly project parent NS coeffs to child NS coeffs

        /// return the NS coefficients if parent and child are the same,
        /// or construct sum coeffs from the parents and "add" zero wavelet coeffs
        /// @param[in]	child	the key whose coeffs we are requesting
        /// @param[in]	parent	the (leaf) key of our function
        /// @param[in]	coeff	the (leaf) coeffs belonging to parent
        /// @return 	coeffs in NS form
        coeffT parent_to_child_NS(const keyT& child, const keyT& parent,
                                  const coeffT& coeff) const;

        /// Return the values when given the coeffs in scaling function basis
        /// @param[in] key the key of the function node (box)
        /// @param[in] coeff the tensor of scaling function coefficients for function node (box)
        /// @return function values for function node (box)
        template <typename Q>
        GenTensor<Q> coeffs2values(const keyT& key, const GenTensor<Q>& coeff) const {
            // PROFILE_MEMBER_FUNC(FunctionImpl); // Too fine grain for routine profiling
            double scale = pow(2.0,0.5*NDIM*key.level())/sqrt(FunctionDefaults<NDIM>::get_cell_volume());
            return transform(coeff,cdata.quad_phit).scale(scale);
        }

        /// convert S or NS coeffs to values on a 2k grid of the children

        /// equivalent to unfiltering the NS coeffs and then converting all child S-coeffs
        /// to values in their respective boxes. If only S coeffs are provided d coeffs are
        /// assumed to be zero. Reverse operation to values2NScoeffs().
        /// @param[in]	key	the key of the current S or NS coeffs, level n
        /// @param[in]	coeff coeffs in S or NS form; if S then d coeffs are assumed zero
        /// @param[in]	s_only	sanity check to avoid unintended discard of d coeffs
        /// @return		function values on the quadrature points of the children of child (!)
        template <typename Q>
        GenTensor<Q> NScoeffs2values(const keyT& key, const GenTensor<Q>& coeff,
                                     const bool s_only) const {
            // PROFILE_MEMBER_FUNC(FunctionImpl); // Too fine grain for routine profiling

            // sanity checks
            MADNESS_ASSERT((coeff.dim(0)==this->get_k()) == s_only);
            MADNESS_ASSERT((coeff.dim(0)==this->get_k()) or (coeff.dim(0)==2*this->get_k()));

            // this is a block-diagonal matrix with the quadrature points on the diagonal
            Tensor<double> quad_phit_2k(2*cdata.k,2*cdata.npt);
            quad_phit_2k(cdata.s[0],cdata.s[0])=cdata.quad_phit;
            quad_phit_2k(cdata.s[1],cdata.s[1])=cdata.quad_phit;

            // the transformation matrix unfilters (cdata.hg) and transforms to values in one step
            const Tensor<double> transf = (s_only)
                ? inner(cdata.hg(Slice(0,k-1),_),quad_phit_2k)	// S coeffs
                : inner(cdata.hg,quad_phit_2k);					// NS coeffs

            // increment the level since the coeffs2values part happens on level n+1
            const double scale = pow(2.0,0.5*NDIM*(key.level()+1))/
                sqrt(FunctionDefaults<NDIM>::get_cell_volume());

           return transform(coeff,transf).scale(scale);
        }

        /// Compute the function values for multiplication

        /// Given S or NS coefficients from a parent cell, compute the value of
        /// the functions at the quadrature points of a child
        /// currently restricted to special cases
        /// @param[in]	child	key of the box in which we compute values
        /// @param[in]	parent	key of the parent box holding the coeffs
        ///	@param[in]	coeff	coeffs of the parent box
        /// @param[in]	s_only	sanity check to avoid unintended discard of d coeffs
        /// @return		function values on the quadrature points of the children of child (!)
        template <typename Q>
        GenTensor<Q> NS_fcube_for_mul(const keyT& child, const keyT& parent,
                                      const GenTensor<Q>& coeff, const bool s_only) const {
            // PROFILE_MEMBER_FUNC(FunctionImpl); // Too fine grain for routine profiling

            // sanity checks
            MADNESS_ASSERT((coeff.dim(0)==this->get_k()) == s_only);
            MADNESS_ASSERT((coeff.dim(0)==this->get_k()) or (coeff.dim(0)==2*this->get_k()));

            // fast return if possible
            //            if (child.level()==parent.level()) return NScoeffs2values(child,coeff,s_only);

            if (s_only) {

                Tensor<double> quad_phi[NDIM];
                // tmp tensor
                Tensor<double> phi1(cdata.k,cdata.npt);

                for (std::size_t d=0; d<NDIM; ++d) {

                    // input is S coeffs (dimension k), output is values on 2*npt grid points
                    quad_phi[d]=Tensor<double>(cdata.k,2*cdata.npt);

                    // for both children of "child" evaluate the Legendre polynomials
                    // first the left child on level n+1 and translations 2l
                    phi_for_mul(parent.level(),parent.translation()[d],
                                child.level()+1, 2*child.translation()[d], phi1);
                    quad_phi[d](_,Slice(0,k-1))=phi1;

                    // next the right child on level n+1 and translations 2l+1
                    phi_for_mul(parent.level(),parent.translation()[d],
                                child.level()+1, 2*child.translation()[d]+1, phi1);
                    quad_phi[d](_,Slice(k,2*k-1))=phi1;
                }

                const double scale = 1.0/sqrt(FunctionDefaults<NDIM>::get_cell_volume());
                return general_transform(coeff,quad_phi).scale(scale);
            }
            MADNESS_EXCEPTION("you should not be here in NS_fcube_for_mul",1);
            return GenTensor<Q>();
        }

        /// convert function values of the a child generation directly to NS coeffs

        /// equivalent to converting the function values to 2^NDIM S coeffs and then
        /// filtering them to NS coeffs. Reverse operation to NScoeffs2values().
        /// @param[in]	key		key of the parent of the generation
        /// @param[in]	values	tensor holding function values of the 2^NDIM children of key
        /// @return		NS coeffs belonging to key
        template <typename Q>
        GenTensor<Q> values2NScoeffs(const keyT& key, const GenTensor<Q>& values) const {
            //PROFILE_MEMBER_FUNC(FunctionImpl);  // Too fine grain for routine profiling

            // sanity checks
            MADNESS_ASSERT(values.dim(0)==2*this->get_k());

            // this is a block-diagonal matrix with the quadrature points on the diagonal
            Tensor<double> quad_phit_2k(2*cdata.npt,2*cdata.k);
            quad_phit_2k(cdata.s[0],cdata.s[0])=cdata.quad_phiw;
            quad_phit_2k(cdata.s[1],cdata.s[1])=cdata.quad_phiw;

            // the transformation matrix unfilters (cdata.hg) and transforms to values in one step
            const Tensor<double> transf=inner(quad_phit_2k,cdata.hgT);

            // increment the level since the values2coeffs part happens on level n+1
            const double scale = pow(0.5,0.5*NDIM*(key.level()+1))
                *sqrt(FunctionDefaults<NDIM>::get_cell_volume());

            return transform(values,transf).scale(scale);
        }

        /// Return the scaling function coeffs when given the function values at the quadrature points
        /// @param[in] key the key of the function node (box)
        /// @return function values for function node (box)
        template <typename Q>
        Tensor<Q> coeffs2values(const keyT& key, const Tensor<Q>& coeff) const {
            // PROFILE_MEMBER_FUNC(FunctionImpl); // Too fine grain for routine profiling
            double scale = pow(2.0,0.5*NDIM*key.level())/sqrt(FunctionDefaults<NDIM>::get_cell_volume());
            return transform(coeff,cdata.quad_phit).scale(scale);
        }

        template <typename Q>
        GenTensor<Q> values2coeffs(const keyT& key, const GenTensor<Q>& values) const {
            // PROFILE_MEMBER_FUNC(FunctionImpl); // Too fine grain for routine profiling
            double scale = pow(0.5,0.5*NDIM*key.level())*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
            return transform(values,cdata.quad_phiw).scale(scale);
        }

        template <typename Q>
        Tensor<Q> values2coeffs(const keyT& key, const Tensor<Q>& values) const {
            // PROFILE_MEMBER_FUNC(FunctionImpl); // Too fine grain for routine profiling
            double scale = pow(0.5,0.5*NDIM*key.level())*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
            return transform(values,cdata.quad_phiw).scale(scale);
        }

        /// Compute the function values for multiplication

        /// Given coefficients from a parent cell, compute the value of
        /// the functions at the quadrature points of a child
        /// @param[in] child the key for the child function node (box)
        /// @param[in] parent the key for the parent function node (box)
        /// @param[in] coeff the coefficients of scaling function basis of the parent box
        template <typename Q>
        Tensor<Q> fcube_for_mul(const keyT& child, const keyT& parent, const Tensor<Q>& coeff) const {
            // PROFILE_MEMBER_FUNC(FunctionImpl); // Too fine grain for routine profiling
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


        /// Compute the function values for multiplication

        /// Given coefficients from a parent cell, compute the value of
        /// the functions at the quadrature points of a child
        /// @param[in] child the key for the child function node (box)
        /// @param[in] parent the key for the parent function node (box)
        /// @param[in] coeff the coefficients of scaling function basis of the parent box
        template <typename Q>
        GenTensor<Q> fcube_for_mul(const keyT& child, const keyT& parent, const GenTensor<Q>& coeff) const {
            // PROFILE_MEMBER_FUNC(FunctionImpl); // Too fine grain for routine profiling
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


        /// Functor for the mul method
        template <typename L, typename R>
        void do_mul(const keyT& key, const Tensor<L>& left, const std::pair< keyT, Tensor<R> >& arg) {
            // PROFILE_MEMBER_FUNC(FunctionImpl); // Too fine grain for routine profiling
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
        }


        /// multiply the values of two coefficient tensors using a custom number of grid points

        /// note both coefficient tensors have to refer to the same key!
        /// @param[in]	c1	a tensor holding coefficients
        /// @param[in]	c2	another tensor holding coeffs
        /// @param[in]	npt	number of grid points (optional, default is cdata.npt)
        /// @return		coefficient tensor holding the product of the values of c1 and c2
        template<typename R>
        Tensor<TENSOR_RESULT_TYPE(T,R)> mul(const Tensor<T>& c1, const Tensor<R>& c2,
                                            const int npt, const keyT& key) const {
            typedef TENSOR_RESULT_TYPE(T,R) resultT;

            const FunctionCommonData<T,NDIM>& cdata2=FunctionCommonData<T,NDIM>::get(npt);

            // construct a tensor with the npt coeffs
            Tensor<T> c11(cdata2.vk), c22(cdata2.vk);
            c11(this->cdata.s0)=c1;
            c22(this->cdata.s0)=c2;

            // it's sufficient to scale once
            double scale = pow(2.0,0.5*NDIM*key.level())/sqrt(FunctionDefaults<NDIM>::get_cell_volume());
            Tensor<T> c1value=transform(c11,cdata2.quad_phit).scale(scale);
            Tensor<R> c2value=transform(c22,cdata2.quad_phit);
            Tensor<resultT> resultvalue(cdata2.vk,false);
            TERNARY_OPTIMIZED_ITERATOR(resultT, resultvalue, T, c1value, R, c2value, *_p0 = *_p1 * *_p2;);

            Tensor<resultT> result=transform(resultvalue,cdata2.quad_phiw);

            // return a copy of the slice to have the tensor contiguous
            return copy(result(this->cdata.s0));
        }


        /// Functor for the binary_op method
        template <typename L, typename R, typename opT>
	  void do_binary_op(const keyT& key, const Tensor<L>& left,
			    const std::pair< keyT, Tensor<R> >& arg,
			    const opT& op) {
            //PROFILE_MEMBER_FUNC(FunctionImpl); // Too fine grain for routine profiling
	  const keyT& rkey = arg.first;
	  const Tensor<R>& rcoeff = arg.second;
	  Tensor<R> rcube = fcube_for_mul(key, rkey, rcoeff);
	  Tensor<L> lcube = fcube_for_mul(key, key, left);

	  Tensor<T> tcube(cdata.vk,false);
	  op(key, tcube, lcube, rcube);
	  double scale = pow(0.5,0.5*NDIM*key.level())*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
	  tcube = transform(tcube,cdata.quad_phiw).scale(scale);
	  coeffs.replace(key, nodeT(coeffT(tcube,targs),false));
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
        /// @param[in] op the unary operator for the coefficients
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
        /// @param[in] op the unary operator for the coefficients
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

        /// Integrate over one particle of a two particle function and get a one particle function
        /// bsp \int g(1,2) \delta(2-1) d2 = f(1)
        /// The overall dimension of g should be even

		/// The operator
        template<std::size_t LDIM>
		void dirac_convolution_op(const keyT &key, const nodeT &node, FunctionImpl<T,LDIM>* f) const {
			// fast return if the node has children (not a leaf node)
			if(node.has_children()) return;

			const implT* g=this;

			// break the 6D key into two 3D keys (may also work for every even dimension)
			Key<LDIM> key1, key2;
			key.break_apart(key1,key2);

			// get the coefficients of the 6D function g
			const coeffT& g_coeff = node.coeff();

			// get the values of the 6D function g
			coeffT g_values = g->coeffs2values(key,g_coeff);

			// Determine rank and k
			const long rank=g_values.rank();
			const long maxk=f->get_k();
			MADNESS_ASSERT(maxk==g_coeff.dim(0));

			// get tensors for particle 1 and 2 (U and V in SVD)
			tensorT vec1=copy(g_values.get_svdtensor().ref_vector(0).reshape(rank,maxk,maxk,maxk));
			tensorT vec2=g_values.get_svdtensor().ref_vector(1).reshape(rank,maxk,maxk,maxk);
			tensorT result(maxk,maxk,maxk);  // should give zero tensor
			// Multiply the values of each U and V vector
			for (long i=0; i<rank; ++i) {
				tensorT c1=vec1(Slice(i,i),_,_,_); // shallow copy (!)
				tensorT c2=vec2(Slice(i,i),_,_,_);
				c1.emul(c2); // this changes vec1 because of shallow copy, but not the g function because of the deep copy made above
				double singular_value_i = g_values.get_svdtensor().weights(i);
				result += (singular_value_i*c1);
			}

			// accumulate coefficients (since only diagonal boxes are used the coefficients get just replaced, but accumulate is needed to create the right tree structure
			tensorT f_coeff = f->values2coeffs(key1,result);
			f->coeffs.task(key1, &FunctionNode<T,LDIM>::accumulate2, f_coeff, f->coeffs, key1, TaskAttributes::hipri());
//             coeffs.task(dest, &nodeT::accumulate2, result, coeffs, dest, TaskAttributes::hipri());


			return;
		}


        template<std::size_t LDIM>
        void do_dirac_convolution(FunctionImpl<T,LDIM>* f, bool fence) const {
            typename dcT::const_iterator end = this->coeffs.end();
            for (typename dcT::const_iterator it=this->coeffs.begin(); it!=end; ++it) {
                // looping through all the leaf(!) coefficients in the NDIM function ("this")
                const keyT& key = it->first;
                const FunctionNode<T,NDIM>& node = it->second;
                if (node.is_leaf()) {
                	// only process the diagonal boxes
        			Key<LDIM> key1, key2;
        			key.break_apart(key1,key2);
        			if(key1 == key2){
                        ProcessID p = coeffs.owner(key);
                        woT::task(p, &implT:: template dirac_convolution_op<LDIM>, key, node, f);
        			}
        		}
            }
            world.gop.fence(); // fence is necessary if trickle down is used afterwards
            // trickle down and undo redundand shouldnt change anything if only the diagonal elements are considered above -> check this
            f->trickle_down(true); // fence must be true otherwise undo_redundant will have trouble
//            f->undo_redundant(true);
            f->verify_tree();
            //if (fence) world.gop.fence(); // unnecessary, fence is activated in undo_redundant

        }


        /// Unary operation applied inplace to the coefficients WITHOUT refinement, optional fence
        /// @param[in] op the unary operator for the coefficients
        template <typename opT>
        void flo_unary_op_node_inplace(const opT& op, bool fence) {
            typedef Range<typename dcT::iterator> rangeT;
//            typedef do_unary_op_value_inplace<opT> xopT;
            world.taskq.for_each<rangeT,opT>(rangeT(coeffs.begin(), coeffs.end()), op);
            if (fence) world.gop.fence();
        }

        /// Unary operation applied inplace to the coefficients WITHOUT refinement, optional fence
        /// @param[in] op the unary operator for the coefficients
        template <typename opT>
        void flo_unary_op_node_inplace(const opT& op, bool fence) const {
            typedef Range<typename dcT::const_iterator> rangeT;
//            typedef do_unary_op_value_inplace<opT> xopT;
            world.taskq.for_each<rangeT,opT>(rangeT(coeffs.begin(), coeffs.end()), op);
            if (fence)
                world.gop.fence();
        }

        /// truncate tree at a certain level
        /// @param[in] max_level truncate tree below this level
        void erase(const Level& max_level);

        /// Returns some asymmetry measure ... no comms
        double check_symmetry_local() const;

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

        /// remove all coefficients of internal nodes
        struct remove_internal_coeffs {
            typedef Range<typename dcT::iterator> rangeT;

            /// constructor need impl for cdata
            remove_internal_coeffs() = default;

            bool operator()(typename rangeT::iterator& it) const {

                nodeT& node = it->second;
                if (node.has_children()) node.clear_coeff();
                return true;
            }
            template <typename Archive> void serialize(const Archive& ar) {}

        };

        /// remove all coefficients of leaf nodes
        struct remove_leaf_coeffs {
            typedef Range<typename dcT::iterator> rangeT;

            /// constructor need impl for cdata
            remove_leaf_coeffs() = default;

            bool operator()(typename rangeT::iterator& it) const {
                nodeT& node = it->second;
                if (not node.has_children()) node.clear_coeff();
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
            do_reduce_rank() = default;
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
            do_check_symmetry_local() : f(0) {}
            do_check_symmetry_local(const implT& f) : f(&f) {}

            /// return the norm of the difference of this node and its "mirror" node
            double operator()(typename rangeT::iterator& it) const {

	      // Temporary fix to GCC whining about out of range access for NDIM!=6
  	      if constexpr(NDIM==6) {
                const keyT& key = it->first;
                const nodeT& fnode = it->second;

                // skip internal nodes
                if (fnode.has_children()) return 0.0;

                if (f->world.size()>1) return 0.0;

                // exchange particles
                std::vector<long> map(NDIM);
                map[0]=3; map[1]=4; map[2]=5;
                map[3]=0; map[4]=1; map[5]=2;

                // make mapped key
                Vector<Translation,NDIM> l;
                for (std::size_t i=0; i<NDIM; ++i) l[map[i]] = key.translation()[i];
                const keyT mapkey(key.level(),l);

                double norm=0.0;


                // hope it's local
                if (f->get_coeffs().probe(mapkey)) {
                    MADNESS_ASSERT(f->get_coeffs().probe(mapkey));
                    const nodeT& mapnode=f->get_coeffs().find(mapkey).get()->second;

//                    bool have_c1=fnode.coeff().has_data() and fnode.coeff().config().has_data();
//                    bool have_c2=mapnode.coeff().has_data() and mapnode.coeff().config().has_data();
                    bool have_c1=fnode.coeff().has_data();
                    bool have_c2=mapnode.coeff().has_data();

                    if (have_c1 and have_c2) {
                        tensorT c1=fnode.coeff().full_tensor_copy();
                        tensorT c2=mapnode.coeff().full_tensor_copy();
                        c2 = copy(c2.mapdim(map));
                        norm=(c1-c2).normf();
                    } else if (have_c1) {
                        tensorT c1=fnode.coeff().full_tensor_copy();
                        norm=c1.normf();
                    } else if (have_c2) {
                        tensorT c2=mapnode.coeff().full_tensor_copy();
                        norm=c2.normf();
                    } else {
                        norm=0.0;
                    }
                } else {
                    norm=fnode.coeff().normf();
                }
                return norm*norm;
	      }
	      else {
		      MADNESS_EXCEPTION("ONLY FOR DIM 6!", 1);
	      }
            }

            double operator()(double a, double b) const {
                return (a+b);
            }

            template <typename Archive> void serialize(const Archive& ar) {
                MADNESS_EXCEPTION("no serialization of do_check_symmetry yet",1);
            }


        };

        /// merge the coefficent boxes of this into result's tree

        /// result+= alpha*this
        /// this and result don't have to have the same distribution or live in the same world
        /// no comm, and the tree should be in an consistent state by virtue
        template<typename Q, typename R>
        struct do_accumulate_trees{
            typedef Range<typename dcT::const_iterator> rangeT;
            FunctionImpl<Q,NDIM>* result=0;
            T alpha=T(1.0);
            do_accumulate_trees() = default;
            do_accumulate_trees(FunctionImpl<Q,NDIM>& result, const T alpha)
                : result(&result), alpha(alpha) {}

            /// return the norm of the difference of this node and its "mirror" node
            bool operator()(typename rangeT::iterator& it) const {

                const keyT& key = it->first;
                const nodeT& node = it->second;
                if (node.has_coeff()) result->get_coeffs().task(key, &nodeT::accumulate,
                    alpha*node.coeff(), result->get_coeffs(), key, result->targs);
                return true;
            }

            template <typename Archive> void serialize(const Archive& ar) {
                MADNESS_EXCEPTION("no serialization of do_accumulate_trees",1);
            }
        };


        /// merge the coefficient boxes of this into other's tree

        /// no comm, and the tree should be in an consistent state by virtue
        /// of FunctionNode::gaxpy_inplace
        template<typename Q, typename R>
        struct do_merge_trees {
            typedef Range<typename dcT::const_iterator> rangeT;
            FunctionImpl<Q,NDIM>* other;
            T alpha;
            R beta;
            do_merge_trees() : other(0) {}
            do_merge_trees(const T alpha, const R beta, FunctionImpl<Q,NDIM>& other)
                : other(&other), alpha(alpha), beta(beta) {}

            /// return the norm of the difference of this node and its "mirror" node
            bool operator()(typename rangeT::iterator& it) const {

                const keyT& key = it->first;
                const nodeT& fnode = it->second;

                // if other's node exists: add this' coeffs to it
                // otherwise insert this' node into other's tree
                typename dcT::accessor acc;
                if (other->get_coeffs().find(acc,key)) {
                    nodeT& gnode=acc->second;
                    gnode.gaxpy_inplace(beta,fnode,alpha);
                } else {
                    nodeT gnode=fnode;
                    gnode.scale(alpha);
                    other->get_coeffs().replace(key,gnode);
                }
                return true;
            }

            template <typename Archive> void serialize(const Archive& ar) {
                MADNESS_EXCEPTION("no serialization of do_merge_trees",1);
            }
        };


        /// map this on f
        struct do_mapdim {
            typedef Range<typename dcT::iterator> rangeT;

            std::vector<long> map;
            implT* f;

            do_mapdim() : f(0) {};
            do_mapdim(const std::vector<long> map, implT& f) : map(map), f(&f) {}

            bool operator()(typename rangeT::iterator& it) const {

                const keyT& key = it->first;
                const nodeT& node = it->second;

                Vector<Translation,NDIM> l;
                for (std::size_t i=0; i<NDIM; ++i) l[map[i]] = key.translation()[i];
                tensorT c = node.coeff().reconstruct_tensor();
                if (c.size()) c = copy(c.mapdim(map));
                coeffT cc(c,f->get_tensor_args());
                f->get_coeffs().replace(keyT(key.level(),l), nodeT(cc,node.has_children()));

                return true;
            }
            template <typename Archive> void serialize(const Archive& ar) {
                MADNESS_EXCEPTION("no serialization of do_mapdim",1);
            }

        };

        /// mirror dimensions of this, write result on f
        struct do_mirror {
            typedef Range<typename dcT::iterator> rangeT;

            std::vector<long> mirror;
            implT* f;

            do_mirror() : f(0) {};
            do_mirror(const std::vector<long> mirror, implT& f) : mirror(mirror), f(&f) {}

            bool operator()(typename rangeT::iterator& it) const {

                const keyT& key = it->first;
                const nodeT& node = it->second;

                // mirror translation index: l_new + l_old = l_max
                Vector<Translation,NDIM> l=key.translation();
                Translation lmax = (Translation(1)<<key.level()) - 1;
                for (std::size_t i=0; i<NDIM; ++i) {
                	if (mirror[i]==-1) l[i]= lmax - key.translation()[i];
                }

                // mirror coefficients: multiply all odd-k slices with -1
                tensorT c = node.coeff().full_tensor_copy();
            	if (c.size()) {
            		std::vector<Slice> s(___);

                	// loop over dimensions and over k
                	for (size_t i=0; i<NDIM; ++i) {
                		std::size_t kmax=c.dim(i);
                		if (mirror[i]==-1) {
                			for (size_t k=1; k<kmax; k+=2) {
                				s[i]=Slice(k,k,1);
                				c(s)*=(-1.0);
                			}
                			s[i]=_;
                		}
                	}
                }
                coeffT cc(c,f->get_tensor_args());
                f->get_coeffs().replace(keyT(key.level(),l), nodeT(cc,node.has_children()));

                return true;
            }
            template <typename Archive> void serialize(const Archive& ar) {
                MADNESS_EXCEPTION("no serialization of do_mirror",1);
            }

        };

        /// mirror dimensions of this, write result on f
        struct do_map_and_mirror {
            typedef Range<typename dcT::iterator> rangeT;

            std::vector<long> map,mirror;
            implT* f;

            do_map_and_mirror() = default;
            do_map_and_mirror(const std::vector<long> map, const std::vector<long> mirror, implT& f)
            		: map(map), mirror(mirror), f(&f) {}

            bool operator()(typename rangeT::iterator& it) const {

                const keyT& key = it->first;
                const nodeT& node = it->second;

                tensorT c = node.coeff().full_tensor_copy();
                Vector<Translation,NDIM> l=key.translation();

                // do the mapping first (if present)
                if (map.size()>0) {
                	Vector<Translation,NDIM> l1=l;
                	for (std::size_t i=0; i<NDIM; ++i) l1[map[i]] = l[i];
                	std::swap(l,l1);
                	if (c.size()) c = copy(c.mapdim(map));
                }

                if (mirror.size()>0) {
					// mirror translation index: l_new + l_old = l_max
                	Vector<Translation,NDIM> l1=l;
					Translation lmax = (Translation(1)<<key.level()) - 1;
					for (std::size_t i=0; i<NDIM; ++i) {
						if (mirror[i]==-1) l1[i]= lmax - l[i];
					}
                	std::swap(l,l1);

                	// mirror coefficients: multiply all odd-k slices with -1
					if (c.size()) {
						std::vector<Slice> s(___);

						// loop over dimensions and over k
						for (size_t i=0; i<NDIM; ++i) {
							std::size_t kmax=c.dim(i);
							if (mirror[i]==-1) {
								for (size_t k=1; k<kmax; k+=2) {
									s[i]=Slice(k,k,1);
									c(s)*=(-1.0);
								}
								s[i]=_;
							}
						}
					}
                }

                coeffT cc(c,f->get_tensor_args());
                f->get_coeffs().replace(keyT(key.level(),l), nodeT(cc,node.has_children()));
                return true;
            }
            template <typename Archive> void serialize(const Archive& ar) {
                MADNESS_EXCEPTION("no serialization of do_mirror",1);
            }

        };



        /// "put" this on g
        struct do_average {
            typedef Range<typename dcT::const_iterator> rangeT;

            implT* g;

            do_average() : g(0) {}
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
            implT* f;

            // constructor takes target precision
            do_change_tensor_type() = default;
            // do_change_tensor_type(const TensorArgs& targs) : targs(targs) {}
            do_change_tensor_type(const TensorArgs& targs, implT& g) : targs(targs), f(&g) {}

            //
            bool operator()(typename rangeT::iterator& it) const {

            	double cpu0=cpu_time();
                nodeT& node = it->second;
                change_tensor_type(node.coeff(),targs);
            	double cpu1=cpu_time();
                f->timer_change_tensor_type.accumulate(cpu1-cpu0);

                return true;

            }
            template <typename Archive> void serialize(const Archive& ar) {}
        };

        struct do_consolidate_buffer {
            typedef Range<typename dcT::iterator> rangeT;

            // threshold for rank reduction / SVD truncation
            TensorArgs targs;

            // constructor takes target precision
            do_consolidate_buffer() = default;
            do_consolidate_buffer(const TensorArgs& targs) : targs(targs) {}
            bool operator()(typename rangeT::iterator& it) const {
                it->second.consolidate_buffer(targs);
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
                    tensorT& t= node.coeff().full_tensor();
                    //double before = t.normf();
                    tensorT values = impl->fcube_for_mul(key, key, t);
                    op(key, values);
                    double scale = pow(0.5,0.5*NDIM*key.level())*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
                    t = transform(values,impl->cdata.quad_phiw).scale(scale);
                    node.coeff()=coeffT(t,impl->get_tensor_args());
                    //double after = t.normf();
                    //madness::print("XOP:", key, before, after);
                }
                return true;
            }
            template <typename Archive> void serialize(const Archive& ar) {}
        };

        template <typename Q, typename R>
        /// @todo I don't know what this does other than a trasform
        void vtransform_doit(const std::shared_ptr< FunctionImpl<R,NDIM> >& right,
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
				if (left->coeffs.is_local(parent))
				  left->coeffs.send(parent, &nodeT::set_has_children_recursive, left->coeffs, parent);
				else
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
        }

        /// Refine multiple functions down to the same finest level

        /// @param v the vector of functions we are refining.
        /// @param key the current node.
        /// @param c the vector of coefficients passed from above.
        void refine_to_common_level(const std::vector<FunctionImpl<T,NDIM>*>& v,
                                    const std::vector<tensorT>& c,
                                    const keyT key);

        /// Inplace operate on many functions (impl's) with an operator within a certain box
        /// @param[in] key the key of the current function node (box)
        /// @param[in] op the operator
        /// @param[in] v the vector of function impl's on which to be operated
        template <typename opT>
        void multiop_values_doit(const keyT& key, const opT& op, const std::vector<implT*>& v) {
            std::vector<tensorT> c(v.size());
            for (unsigned int i=0; i<v.size(); i++) {
                if (v[i]) {
                    coeffT cc = coeffs2values(key, v[i]->coeffs.find(key).get()->second.coeff());
                    c[i]=cc.full_tensor();
                }
            }
            tensorT r = op(key, c);
            coeffs.replace(key, nodeT(coeffT(values2coeffs(key, r),targs),false));
        }

        /// Inplace operate on many functions (impl's) with an operator within a certain box
        /// Assumes all functions have been refined down to the same level
        /// @param[in] op the operator
        /// @param[in] v the vector of function impl's on which to be operated
        template <typename opT>
        void multiop_values(const opT& op, const std::vector<implT*>& v) {
            // rough check on refinement level (ignore non-initialized functions
            for (std::size_t i=1; i<v.size(); ++i) {
                if (v[i] and v[i-1]) {
                    MADNESS_ASSERT(v[i]->coeffs.size()==v[i-1]->coeffs.size());
                }
            }
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

        /// Inplace operate on many functions (impl's) with an operator within a certain box

        /// @param[in] key the key of the current function node (box)
        /// @param[in] op the operator
        /// @param[in] vin the vector of function impl's on which to be operated
        /// @param[out] vout the resulting vector of function impl's
        template <typename opT>
        void multi_to_multi_op_values_doit(const keyT& key, const opT& op,
                const std::vector<implT*>& vin, std::vector<implT*>& vout) {
            std::vector<tensorT> c(vin.size());
            for (unsigned int i=0; i<vin.size(); i++) {
                if (vin[i]) {
                    coeffT cc = coeffs2values(key, vin[i]->coeffs.find(key).get()->second.coeff());
                    c[i]=cc.full_tensor();
                }
            }
            std::vector<tensorT> r = op(key, c);
            MADNESS_ASSERT(r.size()==vout.size());
            for (std::size_t i=0; i<vout.size(); ++i) {
                vout[i]->coeffs.replace(key, nodeT(coeffT(values2coeffs(key, r[i]),targs),false));
            }
        }

        /// Inplace operate on many functions (impl's) with an operator within a certain box

        /// Assumes all functions have been refined down to the same level
        /// @param[in] op the operator
        /// @param[in] vin the vector of function impl's on which to be operated
        /// @param[out] vout the resulting vector of function impl's
        template <typename opT>
        void multi_to_multi_op_values(const opT& op, const std::vector<implT*>& vin,
                std::vector<implT*>& vout, const bool fence=true) {
            // rough check on refinement level (ignore non-initialized functions
            for (std::size_t i=1; i<vin.size(); ++i) {
                if (vin[i] and vin[i-1]) {
                    MADNESS_ASSERT(vin[i]->coeffs.size()==vin[i-1]->coeffs.size());
                }
            }
            typename dcT::iterator end = vin[0]->coeffs.end();
            for (typename dcT::iterator it=vin[0]->coeffs.begin(); it!=end; ++it) {
                const keyT& key = it->first;
                if (it->second.has_coeff())
                    world.taskq.add(*this, &implT:: template multi_to_multi_op_values_doit<opT>,
                            key, op, vin, vout);
                else {
                    // fill result functions with empty box in this key
                    for (implT* it2 : vout) {
                        it2->coeffs.replace(key, nodeT(coeffT(),true));
                    }
                }
            }
            if (fence) world.gop.fence();
        }

        /// Transforms a vector of functions left[i] = sum[j] right[j]*c[j,i] using sparsity
        /// @param[in] vright vector of functions (impl's) on which to be transformed
        /// @param[in] c the tensor (matrix) transformer
        /// @param[in] vleft vector of of the *newly* transformed functions (impl's)
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
        /// @param[in] op the unary operator for the values
        template <typename opT>
        void unary_op_value_inplace(const opT& op, bool fence) {
            typedef Range<typename dcT::iterator> rangeT;
            typedef do_unary_op_value_inplace<opT> xopT;
            world.taskq.for_each<rangeT,xopT>(rangeT(coeffs.begin(), coeffs.end()), xopT(this,op));
            if (fence)
                world.gop.fence();
        }

        // Multiplication assuming same distribution and recursive descent
        /// Both left and right functions are in the scaling function basis
        /// @param[in] key the key to the current function node (box)
        /// @param[in] left the function impl associated with the left function
        /// @param[in] lcin the scaling function coefficients associated with the
        ///            current box in the left function
        /// @param[in] vrightin the vector of function impl's associated with
        ///            the vector of right functions
        /// @param[in] vrcin the vector scaling function coefficients associated with the
        ///            current box in the right functions
        /// @param[out] vresultin the vector of resulting functions (impl's)
        template <typename L, typename R>
        void mulXXveca(const keyT& key,
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
                    lc = it->second.coeff().full_tensor_copy();
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
                        rc = it->second.coeff().full_tensor_copy();
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
                else {  // Interior node
                    result->coeffs.replace(key, nodeT(coeffT(),true));
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
        }

        /// Multiplication using recursive descent and assuming same distribution
        /// Both left and right functions are in the scaling function basis
        /// @param[in] key the key to the current function node (box)
        /// @param[in] left the function impl associated with the left function
        /// @param[in] lcin the scaling function coefficients associated with the
        ///            current box in the left function
        /// @param[in] right the function impl associated with the right function
        /// @param[in] rcin the scaling function coefficients associated with the
        ///            current box in the right function
        template <typename L, typename R>
        void mulXXa(const keyT& key,
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
                    lc = it->second.coeff().reconstruct_tensor();
            }

            Tensor<R> rc = rcin;
            if (rc.size() == 0) {
                riterT it = right->coeffs.find(key).get();
                MADNESS_ASSERT(it != right->coeffs.end());
                rnorm = it->second.get_norm_tree();
                if (it->second.has_coeff())
                    rc = it->second.coeff().reconstruct_tensor();
            }

            // both nodes are leaf nodes: multiply and return
            if (rc.size() && lc.size()) { // Yipee!
                do_mul<L,R>(key, lc, std::make_pair(key,rc));
                return;
            }

            if (tol) {
                if (lc.size())
                    lnorm = lc.normf(); // Otherwise got from norm tree above
                if (rc.size())
                    rnorm = rc.normf();
                if (lnorm*rnorm < truncate_tol(tol, key)) {
                    coeffs.replace(key, nodeT(coeffT(cdata.vk,targs),false)); // Zero leaf node
                    return;
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
        }


        // Binary operation on values using recursive descent and assuming same distribution
        /// Both left and right functions are in the scaling function basis
        /// @param[in] key the key to the current function node (box)
        /// @param[in] left the function impl associated with the left function
        /// @param[in] lcin the scaling function coefficients associated with the
        ///            current box in the left function
        /// @param[in] right the function impl associated with the right function
        /// @param[in] rcin the scaling function coefficients associated with the
        ///            current box in the right function
        /// @param[in] op the binary operator
        template <typename L, typename R, typename opT>
        void binaryXXa(const keyT& key,
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
                    lc = it->second.coeff().reconstruct_tensor();
            }

            Tensor<R> rc = rcin;
            if (rc.size() == 0) {
                riterT it = right->coeffs.find(key).get();
                MADNESS_ASSERT(it != right->coeffs.end());
                if (it->second.has_coeff())
                    rc = it->second.coeff().reconstruct_tensor();
            }

            if (rc.size() && lc.size()) { // Yipee!
                do_binary_op<L,R>(key, lc, std::make_pair(key,rc), op);
                return;
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
        }

        template <typename Q, typename opT>
        struct coeff_value_adaptor {
            typedef typename opT::resultT resultT;
            const FunctionImpl<Q,NDIM>* impl_func;
            opT op;

            coeff_value_adaptor() = default;
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

        /// Out of place unary operation on function impl
        /// The skeleton algorithm should resemble something like
        ///
        /// *this = op(*func)
        ///
        /// @param[in] key the key of the current function node (box)
        /// @param[in] func the function impl on which to be operated
        /// @param[in] op the unary operator
        template <typename Q, typename opT>
        void unaryXXa(const keyT& key,
                      const FunctionImpl<Q,NDIM>* func, const opT& op) {

            //            const Tensor<Q>& fc = func->coeffs.find(key).get()->second.full_tensor_copy();
            const Tensor<Q> fc = func->coeffs.find(key).get()->second.coeff().reconstruct_tensor();

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
        }

        /// Multiplies two functions (impl's) together. Delegates to the mulXXa() method
        /// @param[in] left pointer to the left function impl
        /// @param[in] right pointer to the right function impl
        /// @param[in] tol numerical tolerance
        template <typename L, typename R>
        void mulXX(const FunctionImpl<L,NDIM>* left, const FunctionImpl<R,NDIM>* right, double tol, bool fence) {
            if (world.rank() == coeffs.owner(cdata.key0))
                mulXXa(cdata.key0, left, Tensor<L>(), right, Tensor<R>(), tol);
            if (fence)
                world.gop.fence();

            //verify_tree();
        }

        /// Performs binary operation on two functions (impl's). Delegates to the binaryXXa() method
        /// @param[in] left pointer to the left function impl
        /// @param[in] right pointer to the right function impl
        /// @param[in] op the binary operator
        template <typename L, typename R, typename opT>
        void binaryXX(const FunctionImpl<L,NDIM>* left, const FunctionImpl<R,NDIM>* right,
                      const opT& op, bool fence) {
            if (world.rank() == coeffs.owner(cdata.key0))
                binaryXXa(cdata.key0, left, Tensor<L>(), right, Tensor<R>(), op);
            if (fence)
                world.gop.fence();

            //verify_tree();
        }

        /// Performs unary operation on function impl. Delegates to the unaryXXa() method
        /// @param[in] func function impl of the operand
        /// @param[in] op the unary operator
        template <typename Q, typename opT>
        void unaryXX(const FunctionImpl<Q,NDIM>* func, const opT& op, bool fence) {
            if (world.rank() == coeffs.owner(cdata.key0))
                unaryXXa(cdata.key0, func, op);
            if (fence)
                world.gop.fence();

            //verify_tree();
        }

        /// Performs unary operation on function impl. Delegates to the unaryXXa() method
        /// @param[in] func function impl of the operand
        /// @param[in] op the unary operator
        template <typename Q, typename opT>
        void unaryXXvalues(const FunctionImpl<Q,NDIM>* func, const opT& op, bool fence) {
            if (world.rank() == coeffs.owner(cdata.key0))
                unaryXXa(cdata.key0, func, coeff_value_adaptor<Q,opT>(func,op));
            if (fence)
                world.gop.fence();

            //verify_tree();
        }

        /// Multiplies a function (impl) with a vector of functions (impl's). Delegates to the
        /// mulXXveca() method.
        /// @param[in] left pointer to the left function impl
        /// @param[in] vright vector of pointers to the right function impl's
        /// @param[in] tol numerical tolerance
        /// @param[out] vresult vector of pointers to the resulting function impl's
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
        void put_in_box(ProcessID from, long nl, long ni) const;

        /// Prints summary of data distribution
        void print_info() const;

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
        /// @todo Robert .... help!
        void sock_it_to_me(const keyT& key,
                           const RemoteReference< FutureImpl< std::pair<keyT,coeffT> > >& ref) const;
        /// As above, except
        /// 3) The coeffs are constructed from the avg of nodes further down the tree
        /// @todo Robert .... help!
        void sock_it_to_me_too(const keyT& key,
                               const RemoteReference< FutureImpl< std::pair<keyT,coeffT> > >& ref) const;

        /// @todo help!
        void plot_cube_kernel(archive::archive_ptr< Tensor<T> > ptr,
                              const keyT& key,
                              const coordT& plotlo, const coordT& plothi, const std::vector<long>& npt,
                              bool eval_refine) const;


        /// Evaluate a cube/slice of points ... plotlo and plothi are already in simulation coordinates
        /// No communications
        /// @param[in] plotlo the coordinate of the starting point
        /// @param[in] plothi the coordinate of the ending point
        /// @param[in] npt the number of points in each dimension
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
        void eval(const Vector<double,NDIM>& xin,
                  const keyT& keyin,
                  const typename Future<T>::remote_refT& ref);

        /// Get the depth of the tree at a point in \em simulation coordinates

        /// Only the invoking process will get the result via the
        /// remote reference to a future.  Active messages may be sent
        /// to other nodes.
        ///
        /// This function is a minimally-modified version of eval()
        void evaldepthpt(const Vector<double,NDIM>& xin,
                         const keyT& keyin,
                         const typename Future<Level>::remote_refT& ref);

        /// Get the rank of leaf box of the tree at a point in \em simulation coordinates

        /// Only the invoking process will get the result via the
        /// remote reference to a future.  Active messages may be sent
        /// to other nodes.
        ///
        /// This function is a minimally-modified version of eval()
        void evalR(const Vector<double,NDIM>& xin,
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
        void static tnorm(const tensorT& t, double* lo, double* hi);

        void static tnorm(const GenTensor<T>& t, double* lo, double* hi);

        void static tnorm(const SVDTensor<T>& t, double* lo, double* hi, const int particle);

        // This invoked if node has not been autorefined
        void do_square_inplace(const keyT& key);

        // This invoked if node has been autorefined
        void do_square_inplace2(const keyT& parent, const keyT& child, const tensorT& parent_coeff);

        /// Always returns false (for when autorefine is not wanted)
        bool noautorefine(const keyT& key, const tensorT& t) const;

        /// Returns true if this block of coeffs needs autorefining
        bool autorefine_square_test(const keyT& key, const nodeT& t) const;

        /// Pointwise squaring of function with optional global fence

        /// If not autorefining, local computation only if not fencing.
        /// If autorefining, may result in asynchronous communication.
        void square_inplace(bool fence);
        void abs_inplace(bool fence);
        void abs_square_inplace(bool fence);

        /// is this the same as trickle_down() ?
        void sum_down_spawn(const keyT& key, const coeffT& s);

        /// After 1d push operator must sum coeffs down the tree to restore correct scaling function coefficients
        void sum_down(bool fence);

        /// perform this multiplication: h(1,2) = f(1,2) * g(1)
        template<size_t LDIM>
        struct multiply_op {

            static bool randomize() {return false;}
            typedef CoeffTracker<T,NDIM> ctT;
            typedef CoeffTracker<T,LDIM> ctL;
            typedef multiply_op<LDIM> this_type;

            implT* h;     ///< the result function h(1,2) = f(1,2) * g(1)
            ctT f;
            ctL g;
            int particle;       ///< if g is g(1) or g(2)

            multiply_op() : h(), f(), g(), particle(1) {}

            multiply_op(implT* h1, const ctT& f1, const ctL& g1, const int particle1)
                : h(h1), f(f1), g(g1), particle(particle1) {};

            /// return true if this will be a leaf node

            /// use generalization of tnorm for a GenTensor
            bool screen(const coeffT& fcoeff, const coeffT& gcoeff, const keyT& key) const {
                MADNESS_ASSERT(gcoeff.is_full_tensor());
                MADNESS_ASSERT(fcoeff.is_svd_tensor());
                MADNESS_ASSERT(g.get_impl());
                MADNESS_ASSERT(h);

                double glo=0.0, ghi=0.0, flo=0.0, fhi=0.0;
                g.get_impl()->tnorm(gcoeff.get_tensor(), &glo, &ghi);
                g.get_impl()->tnorm(fcoeff.get_svdtensor(),&flo,&fhi,particle);

                double total_hi=glo*fhi + ghi*flo + fhi*ghi;
                return (total_hi<h->truncate_tol(h->get_thresh(),key));

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

                // get coefficients of the actual FunctionNode
                coeffT coeff1=f.get_impl()->parent_to_child(f.coeff(),f.key(),key);
                coeff1.normalize();
                const coeffT coeff2=g.get_impl()->parent_to_child(g.coeff(),g.key(),gkey);

                // multiplication is done in TT_2D
                coeffT coeff1_2D=coeff1.convert(TensorArgs(h->get_thresh(),TT_2D));
                coeff1_2D.normalize();

                bool is_leaf=screen(coeff1_2D,coeff2,key);
                if (key.level()<2) is_leaf=false;

                coeffT hcoeff;
                if (is_leaf) {

                    // convert coefficients to values
                    coeffT hvalues=f.get_impl()->coeffs2values(key,coeff1_2D);
                    coeffT gvalues=g.get_impl()->coeffs2values(gkey,coeff2);

                    // perform multiplication
                    coeffT result_val=h->multiply(hvalues,gvalues,particle-1);

                    hcoeff=h->values2coeffs(key,result_val);

                    // conversion on coeffs, not on values, because it implies truncation!
                    if (not hcoeff.is_of_tensortype(h->get_tensor_type()))
                        hcoeff=hcoeff.convert(h->get_tensor_args());
                }

                return std::pair<bool,coeffT> (is_leaf,hcoeff);
            }

            this_type make_child(const keyT& child) const {

                // break key into particles
                Key<LDIM> key1, key2;
                child.break_apart(key1,key2);
                const Key<LDIM> gkey= (particle==1) ? key1 : key2;

                return this_type(h,f.make_child(child),g.make_child(gkey),particle);
            }

            Future<this_type> activate() const {
            	Future<ctT> f1=f.activate();
            	Future<ctL> g1=g.activate();
                return h->world.taskq.add(detail::wrap_mem_fn(*const_cast<this_type *> (this),
                                          &this_type::forward_ctor),h,f1,g1,particle);
            }

            this_type forward_ctor(implT* h1, const ctT& f1, const ctL& g1, const int particle) {
            	return this_type(h1,f1,g1,particle);
            }

            template <typename Archive> void serialize(const Archive& ar) {
                ar & h & f & g & particle;
            }
        };


        /// add two functions f and g: result=alpha * f  +  beta * g
        struct add_op {

            typedef CoeffTracker<T,NDIM> ctT;
            typedef add_op this_type;

            bool randomize() const {return false;}

            /// tracking coeffs of first and second addend
            ctT f,g;
            /// prefactor for f, g
            double alpha, beta;

            add_op() = default;
            add_op(const ctT& f, const ctT& g, const double alpha, const double beta)
                : f(f), g(g), alpha(alpha), beta(beta){}

            /// if we are at the bottom of the trees, return the sum of the coeffs
            std::pair<bool,coeffT> operator()(const keyT& key) const {

                bool is_leaf=(f.is_leaf() and g.is_leaf());
                if (not is_leaf) return std::pair<bool,coeffT> (is_leaf,coeffT());

                coeffT fcoeff=f.get_impl()->parent_to_child(f.coeff(),f.key(),key);
                coeffT gcoeff=g.get_impl()->parent_to_child(g.coeff(),g.key(),key);
                coeffT hcoeff=copy(fcoeff);
                hcoeff.gaxpy(alpha,gcoeff,beta);
                hcoeff.reduce_rank(f.get_impl()->get_tensor_args().thresh);
                return std::pair<bool,coeffT> (is_leaf,hcoeff);
            }

            this_type make_child(const keyT& child) const {
                return this_type(f.make_child(child),g.make_child(child),alpha,beta);
            }

            /// retrieve the coefficients (parent coeffs might be remote)
            Future<this_type> activate() const {
            	Future<ctT> f1=f.activate();
            	Future<ctT> g1=g.activate();
                return f.get_impl()->world.taskq.add(detail::wrap_mem_fn(*const_cast<this_type *> (this),
                                                     &this_type::forward_ctor),f1,g1,alpha,beta);
            }

            /// taskq-compatible ctor
            this_type forward_ctor(const ctT& f1, const ctT& g1, const double alpha, const double beta) {
            	return this_type(f1,g1,alpha,beta);
            }

            template <typename Archive> void serialize(const Archive& ar) {
                ar & f & g & alpha & beta;
            }

        };

        /// multiply f (a pair function of NDIM) with an orbital g (LDIM=NDIM/2)

        /// as in (with h(1,2)=*this) : h(1,2) = g(1) * f(1,2)
        /// use tnorm as a measure to determine if f (=*this) must be refined
        /// @param[in]  f           the NDIM function f=f(1,2)
        /// @param[in]  g           the LDIM function g(1) (or g(2))
        /// @param[in]  particle    1 or 2, as in g(1) or g(2)
        template<size_t LDIM>
        void multiply(const implT* f, const FunctionImpl<T,LDIM>* g, const int particle) {

                CoeffTracker<T,NDIM> ff(f);
                CoeffTracker<T,LDIM> gg(g);

                typedef multiply_op<LDIM> coeff_opT;
                coeff_opT coeff_op(this,ff,gg,particle);

                typedef insert_op<T,NDIM> apply_opT;
                apply_opT apply_op(this);

            keyT key0=f->cdata.key0;
			if (world.rank() == coeffs.owner(key0)) {
                ProcessID p= coeffs.owner(key0);
                woT::task(p, &implT:: template forward_traverse<coeff_opT,apply_opT>, coeff_op, apply_op, key0);
            }

            set_tree_state(reconstructed);
        }

        /// Hartree product of two LDIM functions to yield a NDIM = 2*LDIM function
        template<size_t LDIM, typename leaf_opT>
        struct hartree_op {
            bool randomize() const {return false;}

            typedef hartree_op<LDIM,leaf_opT> this_type;
            typedef CoeffTracker<T,LDIM> ctL;

            implT* result; 	    ///< where to construct the pair function
            ctL p1, p2;			///< tracking coeffs of the two lo-dim functions
            leaf_opT leaf_op;   ///< determine if a given node will be a leaf node

            // ctor
            hartree_op() : result() {}
            hartree_op(implT* result, const ctL& p11, const ctL& p22, const leaf_opT& leaf_op)
                : result(result), p1(p11), p2(p22), leaf_op(leaf_op) {
                MADNESS_ASSERT(LDIM+LDIM==NDIM);
            }

            std::pair<bool,coeffT> operator()(const Key<NDIM>& key) const {

                // break key into particles (these are the child keys, with datum1/2 come the parent keys)
                Key<LDIM> key1,key2;
                key.break_apart(key1,key2);

                // this returns the appropriate NS coeffs for key1 and key2 resp.
            	const coeffT fcoeff=p1.coeff(key1);
                const coeffT gcoeff=p2.coeff(key2);
                bool is_leaf=leaf_op(key,fcoeff.full_tensor(),gcoeff.full_tensor());
                if (not is_leaf) return std::pair<bool,coeffT> (is_leaf,coeffT());

                // extract the sum coeffs from the NS coeffs
                const coeffT s1=fcoeff(p1.get_impl()->cdata.s0);
                const coeffT s2=gcoeff(p2.get_impl()->cdata.s0);

                // new coeffs are simply the hartree/kronecker/outer product --
                coeffT coeff=outer(s1,s2,result->get_tensor_args());
                // no post-determination
                //                is_leaf=leaf_op(key,coeff);
                return std::pair<bool,coeffT>(is_leaf,coeff);
            }

            this_type make_child(const keyT& child) const {

                // break key into particles
                Key<LDIM> key1, key2;
                child.break_apart(key1,key2);

                return this_type(result,p1.make_child(key1),p2.make_child(key2),leaf_op);
            }

            Future<this_type> activate() const {
            	Future<ctL> p11=p1.activate();
            	Future<ctL> p22=p2.activate();
                return result->world.taskq.add(detail::wrap_mem_fn(*const_cast<this_type *> (this),
                                               &this_type::forward_ctor),result,p11,p22,leaf_op);
            }

            this_type forward_ctor(implT* result1, const ctL& p11, const ctL& p22, const leaf_opT& leaf_op) {
            	return this_type(result1,p11,p22,leaf_op);
            }

            template <typename Archive> void serialize(const Archive& ar) {
                ar & result & p1 & p2 & leaf_op;
            }
        };

        /// traverse a non-existing tree

        /// part II: activate coeff_op, i.e. retrieve all the necessary remote boxes (communication)
        /// @param[in]	coeff_op	operator making the coefficients that needs activation
        /// @param[in]	apply_op	just passing thru
        /// @param[in]	key			the key we are working on
        template<typename coeff_opT, typename apply_opT>
        void forward_traverse(const coeff_opT& coeff_op, const apply_opT& apply_op, const keyT& key) const {
            MADNESS_ASSERT(coeffs.is_local(key));
            Future<coeff_opT> active_coeff=coeff_op.activate();
            woT::task(world.rank(), &implT:: template traverse_tree<coeff_opT,apply_opT>, active_coeff, apply_op, key);
        }


        /// traverse a non-existing tree

        /// part I: make the coefficients, process them and continue the recursion if necessary
        /// @param[in]	coeff_op	operator making the coefficients and determining them being leaves
        /// @param[in]	apply_op	operator processing the coefficients
        /// @param[in]	key			the key we are currently working on
        template<typename coeff_opT, typename apply_opT>
        void traverse_tree(const coeff_opT& coeff_op, const apply_opT& apply_op, const keyT& key) const {
            MADNESS_ASSERT(coeffs.is_local(key));

            typedef typename std::pair<bool,coeffT> argT;
            const argT arg=coeff_op(key);
            apply_op.operator()(key,arg.second,arg.first);

            const bool has_children=(not arg.first);
            if (has_children) {
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT& child=kit.key();
                    coeff_opT child_op=coeff_op.make_child(child);
                    // spawn activation where child is local
                    ProcessID p=coeffs.owner(child);

                    void (implT::*ft)(const coeff_opT&, const apply_opT&, const keyT&) const = &implT::forward_traverse<coeff_opT,apply_opT>;

                    woT::task(p, ft, child_op, apply_op, child);
                }
            }
        }


        /// given two functions of LDIM, perform the Hartree/Kronecker/outer product

        /// |Phi(1,2)> = |phi(1)> x |phi(2)>
        /// @param[in]	p1	FunctionImpl of particle 1
        /// @param[in]	p2	FunctionImpl of particle 2
        /// @param[in]	leaf_op	operator determining of a given box will be a leaf
        template<std::size_t LDIM, typename leaf_opT>
        void hartree_product(const std::vector<std::shared_ptr<FunctionImpl<T,LDIM>>> p1,
                             const std::vector<std::shared_ptr<FunctionImpl<T,LDIM>>> p2,
                             const leaf_opT& leaf_op, bool fence) {
            MADNESS_CHECK_THROW(p1.size()==p2.size(),"hartree_product: p1 and p2 must have the same size");
            for (auto& p : p1) MADNESS_CHECK(p->is_nonstandard() or p->is_nonstandard_with_leaves());
            for (auto& p : p2) MADNESS_CHECK(p->is_nonstandard() or p->is_nonstandard_with_leaves());

            const keyT key0=cdata.key0;

            for (std::size_t i=0; i<p1.size(); ++i) {
                if (world.rank() == this->get_coeffs().owner(key0)) {

                    // prepare the CoeffTracker
                    CoeffTracker<T,LDIM> iap1(p1[i].get());
                    CoeffTracker<T,LDIM> iap2(p2[i].get());

                    // the operator making the coefficients
                    typedef hartree_op<LDIM,leaf_opT> coeff_opT;
                    coeff_opT coeff_op(this,iap1,iap2,leaf_op);

                    // this operator simply inserts the coeffs into this' tree
//                typedef insert_op<T,NDIM> apply_opT;
                    typedef accumulate_op<T,NDIM> apply_opT;
                    apply_opT apply_op(this);

                    woT::task(world.rank(), &implT:: template forward_traverse<coeff_opT,apply_opT>,
                              coeff_op, apply_op, cdata.key0);

                }
            }

            set_tree_state(redundant_after_merge);
            if (fence) world.gop.fence();
        }


        template <typename opT, typename R>
        void
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
                            return; // Hard zero means finished!
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
                    const Tensor<R>& c = node.coeff().full_tensor_copy();
                    woT::task(me, &implT:: template apply_1d_realspace_push_op<opT,R>,
                              archive::archive_ptr<const opT>(&op), axis, key, c);
                }
            }
            if (fence) world.gop.fence();
        }

        void forward_do_diff1(const DerivativeBase<T,NDIM>* D,
                              const implT* f,
                              const keyT& key,
                              const std::pair<keyT,coeffT>& left,
                              const std::pair<keyT,coeffT>& center,
                              const std::pair<keyT,coeffT>& right);

        void do_diff1(const DerivativeBase<T,NDIM>* D,
                      const implT* f,
                      const keyT& key,
                      const std::pair<keyT,coeffT>& left,
                      const std::pair<keyT,coeffT>& center,
                      const std::pair<keyT,coeffT>& right);

        // Called by result function to differentiate f
        void diff(const DerivativeBase<T,NDIM>* D, const implT* f, bool fence);

        /// Returns key of general neighbor enforcing BC

        /// Out of volume keys are mapped to enforce the BC as follows.
        ///   * Periodic BC map back into the volume and return the correct key
        ///   * non-periodic BC - returns invalid() to indicate out of volume
        keyT neighbor(const keyT& key, const keyT& disp, const array_of_bools<NDIM>& is_periodic) const;

        /// Returns key of general neighbor that resides in-volume

        /// Out of volume keys are mapped to invalid()
        keyT neighbor_in_volume(const keyT& key, const keyT& disp) const;

        /// find_me. Called by diff_bdry to get coefficients of boundary function
        Future< std::pair<keyT,coeffT> > find_me(const keyT& key) const;

        /// return the a std::pair<key, node>, which MUST exist
        std::pair<Key<NDIM>,ShallowNode<T,NDIM> > find_datum(keyT key) const;

        /// multiply the ket with a one-electron potential rr(1,2)= f(1,2)*g(1)

        /// @param[in]	val_ket	function values of f(1,2)
        /// @param[in]	val_pot	function values of g(1)
        /// @param[in]	particle	if 0 then g(1), if 1 then g(2)
        /// @return		the resulting function values
        coeffT multiply(const coeffT& val_ket, const coeffT& val_pot, int particle) const;


        /// given several coefficient tensors, assemble a result tensor

        /// the result looks like: 	(v(1,2) + v(1) + v(2)) |ket(1,2)>
        /// or 						(v(1,2) + v(1) + v(2)) |p(1) p(2)>
        /// i.e. coefficients for the ket and coefficients for the two particles are
        /// mutually exclusive. All potential terms are optional, just pass in empty coeffs.
        /// @param[in]	key			the key of the FunctionNode to which these coeffs belong
        /// @param[in]	coeff_ket	coefficients of the ket
        /// @param[in]	vpotential1	function values of the potential for particle 1
        /// @param[in]	vpotential2	function values of the potential for particle 2
        /// @param[in]	veri		function values for the 2-particle potential
        coeffT assemble_coefficients(const keyT& key, const coeffT& coeff_ket,
                                     const coeffT& vpotential1, const coeffT& vpotential2,
                                     const tensorT& veri) const;



        template<std::size_t LDIM>
        struct pointwise_multiplier {
        	coeffT val_lhs, coeff_lhs;
            long oversampling=1;
        	double error=0.0;
        	double lo=0.0, hi=0.0, lo1=0.0, hi1=0.0, lo2=0.0, hi2=0.0;

    	    pointwise_multiplier() = default;
        	pointwise_multiplier(const Key<NDIM> key, const coeffT& clhs) : coeff_lhs(clhs) {
                const auto& fcf=FunctionCommonFunctionality<T,NDIM>(coeff_lhs.dim(0));
        		val_lhs=fcf.coeffs2values(key,coeff_lhs);
        		error=0.0;
        		implT::tnorm(coeff_lhs,&lo,&hi);
                if (coeff_lhs.is_svd_tensor()) {
                    FunctionImpl<T,LDIM>::tnorm(coeff_lhs.get_svdtensor(),&lo1,&hi1,1);
                    FunctionImpl<T,LDIM>::tnorm(coeff_lhs.get_svdtensor(),&lo2,&hi2,2);
                }
        	}

        	/// multiply values of rhs and lhs, result on rhs, rhs and lhs are of the same dimensions
        	tensorT operator()(const Key<NDIM> key, const tensorT& coeff_rhs) {

				MADNESS_ASSERT(coeff_rhs.dim(0)==coeff_lhs.dim(0));
                auto fcf=FunctionCommonFunctionality<T,NDIM>(coeff_lhs.dim(0));

				// the tnorm estimate is not tight enough to be efficient, better use oversampling
				bool use_tnorm=false;
        		if (use_tnorm) {
					double rlo, rhi;
					implT::tnorm(coeff_rhs,&rlo,&rhi);
					error = hi*rlo + rhi*lo + rhi*hi;
					tensorT val_rhs=fcf.coeffs2values(key, coeff_rhs);
					val_rhs.emul(val_lhs.full_tensor_copy());
					return fcf.values2coeffs(key,val_rhs);
        		} else {	// use quadrature of order k+1

    	            auto& cdata=FunctionCommonData<T,NDIM>::get(coeff_rhs.dim(0));		// npt=k+1
                    auto& cdata_npt=FunctionCommonData<T,NDIM>::get(coeff_rhs.dim(0)+oversampling);		// npt=k+1
                	FunctionCommonFunctionality<T,NDIM> fcf_hi_npt(cdata_npt);

    	            // coeffs2values for rhs: k -> npt=k+1
		            tensorT coeff1(cdata_npt.vk);
		            coeff1(cdata.s0)=coeff_rhs;		// s0 is smaller than vk!
		            tensorT val_rhs_k1=fcf_hi_npt.coeffs2values(key,coeff1);

		            // coeffs2values for lhs: k -> npt=k+1
		            tensorT coeff_lhs_k1(cdata_npt.vk);
		            coeff_lhs_k1(cdata.s0)=coeff_lhs.full_tensor_copy();
		            tensorT val_lhs_k1=fcf_hi_npt.coeffs2values(key,coeff_lhs_k1);

		            // multiply
		            val_lhs_k1.emul(val_rhs_k1);

                    // values2coeffs: npt = k+1-> k
		            tensorT result1=fcf_hi_npt.values2coeffs(key,val_lhs_k1);

                    // extract coeffs up to k
                    tensorT result=copy(result1(cdata.s0));
                    result1(cdata.s0)=0.0;
                    error=result1.normf();
                    return result;
        		}
        	}

        	/// multiply values of rhs and lhs, result on rhs, rhs and lhs are of differnet dimensions
        	coeffT operator()(const Key<NDIM> key, const tensorT& coeff_rhs, const int particle) {
                Key<LDIM> key1, key2;
                key.break_apart(key1,key2);
                const long k=coeff_rhs.dim(0);
                auto& cdata=FunctionCommonData<T,NDIM>::get(k);
                auto& cdata_lowdim=FunctionCommonData<T,LDIM>::get(k);
            	FunctionCommonFunctionality<T,LDIM> fcf_lo(cdata_lowdim);
            	FunctionCommonFunctionality<T,NDIM> fcf_hi(cdata);
            	FunctionCommonFunctionality<T,LDIM> fcf_lo_npt(k+oversampling);
            	FunctionCommonFunctionality<T,NDIM> fcf_hi_npt(k+oversampling);


            	// make hi-dim values from lo-dim coeff_rhs on npt grid points
                tensorT ones=tensorT(fcf_lo_npt.cdata.vk);
                ones=1.0;

                tensorT coeff_rhs_npt1(fcf_lo_npt.cdata.vk);
                coeff_rhs_npt1(fcf_lo.cdata.s0)=coeff_rhs;
                tensorT val_rhs_npt1=fcf_lo_npt.coeffs2values(key1,coeff_rhs_npt1);

                TensorArgs targs(-1.0,TT_2D);
                coeffT val_rhs;
                if (particle==1) val_rhs=outer(val_rhs_npt1,ones,targs);
                if (particle==2) val_rhs=outer(ones,val_rhs_npt1,targs);

            	// make values from hi-dim coeff_lhs on npt grid points
	            coeffT coeff_lhs_k1(fcf_hi_npt.cdata.vk,coeff_lhs.tensor_type());
	            coeff_lhs_k1(fcf_hi.cdata.s0)+=coeff_lhs;
	            coeffT val_lhs_npt=fcf_hi_npt.coeffs2values(key,coeff_lhs_k1);

	            // multiply
	            val_lhs_npt.emul(val_rhs);

                // values2coeffs: npt = k+1-> k
	            coeffT result1=fcf_hi_npt.values2coeffs(key,val_lhs_npt);

                // extract coeffs up to k
	            coeffT result=copy(result1(cdata.s0));
                result1(cdata.s0)=0.0;
                error=result1.normf();
                return result;
        	}

            template <typename Archive> void serialize(const Archive& ar) {
                ar & error & lo & lo1 & lo2 & hi & hi1& hi2 & val_lhs & coeff_lhs;
            }


        };

        /// given a ket and the 1- and 2-electron potentials, construct the function V phi

        /// small memory footstep version of Vphi_op: use the NS form to have information
        /// about parent and children to determine if a box is a leaf. This will require
        /// compression of the constituent functions, which will lead to more memory usage
        /// there, but will avoid oversampling of the result function.
        template<typename opT, size_t LDIM>
        struct Vphi_op_NS {

            bool randomize() const {return true;}

            typedef Vphi_op_NS<opT,LDIM> this_type;
            typedef CoeffTracker<T,NDIM> ctT;
            typedef CoeffTracker<T,LDIM> ctL;

            implT* result;  		///< where to construct Vphi, no need to track parents
            opT leaf_op;    	    ///< deciding if a given FunctionNode will be a leaf node
            ctT iaket;				///< the ket of a pair function (exclusive with p1, p2)
            ctL iap1, iap2;			///< the particles 1 and 2 (exclusive with ket)
            ctL iav1, iav2;			///< potentials for particles 1 and 2
            const implT* eri;		///< 2-particle potential, must be on-demand

            bool have_ket() const {return iaket.get_impl();}
            bool have_v1() const {return iav1.get_impl();}
            bool have_v2() const {return iav2.get_impl();}
            bool have_eri() const {return eri;}

            void accumulate_into_result(const Key<NDIM>& key, const coeffT& coeff) const {
                result->get_coeffs().task(key, &nodeT::accumulate, coeff, result->get_coeffs(), key, result->targs);
            }

            // ctor
            Vphi_op_NS() : result(), eri() {}
            Vphi_op_NS(implT* result, const opT& leaf_op, const ctT& iaket,
                       const ctL& iap1, const ctL& iap2, const ctL& iav1, const ctL& iav2,
                       const implT* eri)
                    : result(result), leaf_op(leaf_op), iaket(iaket), iap1(iap1), iap2(iap2)
                    , iav1(iav1), iav2(iav2), eri(eri) {

                // 2-particle potential must be on-demand
                if (eri) MADNESS_ASSERT(eri->is_on_demand());
            }

            /// make and insert the coefficients into result's tree
            std::pair<bool,coeffT> operator()(const Key<NDIM>& key) const {

                MADNESS_ASSERT(result->get_coeffs().is_local(key));
                if(leaf_op.do_pre_screening()){
                    // this means that we only construct the boxes which are leaf boxes from the other function in the leaf_op
                    if(leaf_op.pre_screening(key)){
                        // construct sum_coefficients, insert them and leave
                        auto [sum_coeff, error]=make_sum_coeffs(key);
                        accumulate_into_result(key,sum_coeff);
                        return std::pair<bool,coeffT> (true,coeffT());
                    }else{
                        return continue_recursion(std::vector<bool>(1<<NDIM,false),tensorT(),key);
                    }
                }

                // this means that the function has to be completely constructed and not mirrored by another function

                // if the initial level is not reached then this must not be a leaf box
                size_t il = result->get_initial_level();
                if(FunctionDefaults<NDIM>::get_refine()) il+=1;
                if(key.level()<int(il)){
                    return continue_recursion(std::vector<bool>(1<<NDIM,false),tensorT(),key);
                }
                // if further refinement is needed (because we are at a special box, special point)
                // and the special_level is not reached then this must not be a leaf box
                if(key.level()<result->get_special_level() and leaf_op.special_refinement_needed(key)){
                    return continue_recursion(std::vector<bool>(1<<NDIM,false),tensorT(),key);
                }

                auto [sum_coeff,error]=make_sum_coeffs(key);

                // coeffs are leaf (for whatever reason), insert into tree and stop recursion
                if(leaf_op.post_screening(key,sum_coeff)){
                    accumulate_into_result(key,sum_coeff);
                    return std::pair<bool,coeffT> (true,coeffT());
                }

                // coeffs are accurate, insert into tree and stop recursion
                if(error<result->truncate_tol(result->get_thresh(),key)){
                    accumulate_into_result(key,sum_coeff);
                    return std::pair<bool,coeffT> (true,coeffT());
                }

                // coeffs are inaccurate, continue recursion
                std::vector<bool> child_is_leaf(1<<NDIM,false);
                return continue_recursion(child_is_leaf,tensorT(),key);
            }


            /// loop over all children and either insert their sum coeffs or continue the recursion

            /// @param[in]	child_is_leaf	for each child: is it a leaf?
            /// @param[in]	coeffs	coefficient tensor with 2^N sum coeffs (=unfiltered NS coeffs)
            /// @param[in]	key		the key for the NS coeffs (=parent key of the children)
            /// @return		to avoid recursion outside this return: std::pair<is_leaf,coeff> = true,coeffT()
            std::pair<bool,coeffT> continue_recursion(const std::vector<bool> child_is_leaf,
                                                      const tensorT& coeffs, const keyT& key) const {
                std::size_t i=0;
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit, ++i) {
                    keyT child=kit.key();
                    bool is_leaf=child_is_leaf[i];

                    if (is_leaf) {
                        // insert the sum coeffs
                        insert_op<T,NDIM> iop(result);
                        iop(child,coeffT(copy(coeffs(result->child_patch(child))),result->get_tensor_args()),is_leaf);
                    } else {
                        this_type child_op=this->make_child(child);
                        noop<T,NDIM> no;
                        // spawn activation where child is local
                        ProcessID p=result->get_coeffs().owner(child);

                        void (implT::*ft)(const Vphi_op_NS<opT,LDIM>&, const noop<T,NDIM>&, const keyT&) const = &implT:: template forward_traverse< Vphi_op_NS<opT,LDIM>, noop<T,NDIM> >;
                        result->task(p, ft, child_op, no, child);
                    }
                }
                // return e sum coeffs; also return always is_leaf=true:
                // the recursion is continued within this struct, not outside in traverse_tree!
                return std::pair<bool,coeffT> (true,coeffT());
            }

            tensorT eri_coeffs(const keyT& key) const {
                MADNESS_ASSERT(have_eri());
                PROFILE_BLOCK(eri_coeffs);
                if (eri->get_functor()->provides_coeff()) {
                    return eri->get_functor()->coeff(key).full_tensor();
                } else {
                    tensorT val_eri(eri->cdata.vk);
                    eri->fcube(key,*(eri->get_functor()),eri->cdata.quad_x,val_eri);
                    return eri->values2coeffs(key,val_eri);
                }
            }

            /// the error is computed from the d coefficients of the constituent functions

            /// the result is h_n = P_n(f g), computed as h_n \approx Pn(f_n g_n)
            /// its error is therefore
            /// h_n = (f g)_n = ((Pn(f) + Qn(f)) (Pn(g) + Qn(g))
            ///               = Pn(fn gn) + Qn(fn gn) + Pn(f) Qn(g) + Qn(f) Pn(g) + Qn(f) Pn(g)
            /// the first term is what we compute, the second term is estimated by tnorm (in another function),
            /// the third to last terms are estimated in this function by e.g.: Qn(f)Pn(g) < ||Qn(f)|| ||Pn(g)||
            double compute_error_from_inaccurate_refinement(const keyT& key,
                                                            const tensorT& ceri) const {
                double error = 0.0;
                Key<LDIM> key1, key2;
                key.break_apart(key1,key2);

                PROFILE_BLOCK(compute_error);
                double dnorm_ket, snorm_ket;
                if (have_ket()) {
                    snorm_ket=iaket.coeff(key).normf();
                    dnorm_ket=iaket.dnorm(key);
                } else {
                    double s1=iap1.coeff(key1).normf();
                    double s2=iap2.coeff(key2).normf();
                    double d1=iap1.dnorm(key1);
                    double d2=iap2.dnorm(key2);
                    snorm_ket=s1*s2;
                    dnorm_ket=s1*d2 + s2*d1 + d1*d2;
                }

                if (have_v1()) {
                    double snorm=iav1.coeff(key1).normf();
                    double dnorm=iav1.dnorm(key1);
                    error+=snorm*dnorm_ket + dnorm*snorm_ket + dnorm*dnorm_ket;
                }
                if (have_v2()) {
                    double snorm=iav2.coeff(key2).normf();
                    double dnorm=iav2.dnorm(key2);
                    error+=snorm*dnorm_ket + dnorm*snorm_ket + dnorm*dnorm_ket;
                }
                if (have_eri()) {
                    tensorT s_coeffs=ceri(result->cdata.s0);
                    double snorm=s_coeffs.normf();
                    tensorT d=copy(ceri);
                    d(result->cdata.s0)=0.0;
                    double dnorm=d.normf();
                    error+=snorm*dnorm_ket + dnorm*snorm_ket + dnorm*dnorm_ket;
                }

                bool no_potential=not ((have_v1() or have_v2() or have_eri()));
                if (no_potential) {
                    error=dnorm_ket;
                }
                return error;
            }

            /// make the sum coeffs for key
            std::pair<coeffT,double> make_sum_coeffs(const keyT& key) const {
                PROFILE_BLOCK(make_sum_coeffs);
                // break key into particles
                Key<LDIM> key1, key2;
                key.break_apart(key1,key2);

        		// bool printme=(int(key.translation()[0])==int(std::pow(key.level(),2)/2)) and
        		// 		(int(key.translation()[1])==int(std::pow(key.level(),2)/2)) and
			// 			(int(key.translation()[2])==int(std::pow(key.level(),2)/2));

//        		printme=false;

                // get/make all coefficients
                const coeffT coeff_ket = (iaket.get_impl()) ? iaket.coeff(key)
                                                            : outer(iap1.coeff(key1),iap2.coeff(key2),result->get_tensor_args());
                const coeffT cpot1 = (have_v1()) ? iav1.coeff(key1) : coeffT();
                const coeffT cpot2 = (have_v2()) ? iav2.coeff(key2) : coeffT();
                const tensorT ceri = (have_eri()) ? eri_coeffs(key) : tensorT();

                // compute first part of the total error
                double refine_error=compute_error_from_inaccurate_refinement(key,ceri);
                double error=refine_error;

        		// prepare the multiplication
        		pointwise_multiplier<LDIM> pm(key,coeff_ket);

        		// perform the multiplication, compute tnorm part of the total error
        		coeffT cresult(result->cdata.vk,result->get_tensor_args());
        		if (have_v1()) {
        			cresult+=pm(key,cpot1.get_tensor(),1);
        			error+=pm.error;
            	}
        		if (have_v2()) {
        			cresult+=pm(key,cpot2.get_tensor(),2);
        			error+=pm.error;
        		}

                if (have_eri()) {
                    tensorT result1=cresult.full_tensor_copy();
                    result1+=pm(key,copy(ceri(result->cdata.s0)));
                    cresult=coeffT(result1,result->get_tensor_args());
                    error+=pm.error;
                } else {
                    cresult.reduce_rank(result->get_tensor_args().thresh);
                }
                if ((not have_v1()) and (not have_v2()) and (not have_eri())) {
                    cresult=coeff_ket;
                }

                return std::make_pair(cresult,error);
            }

            this_type make_child(const keyT& child) const {

                // break key into particles
                Key<LDIM> key1, key2;
                child.break_apart(key1,key2);

                return this_type(result,leaf_op,iaket.make_child(child),
                                 iap1.make_child(key1),iap2.make_child(key2),
                                 iav1.make_child(key1),iav2.make_child(key2),eri);
            }

            Future<this_type> activate() const {
                Future<ctT> iaket1=iaket.activate();
                Future<ctL> iap11=iap1.activate();
                Future<ctL> iap21=iap2.activate();
                Future<ctL> iav11=iav1.activate();
                Future<ctL> iav21=iav2.activate();
                return result->world.taskq.add(detail::wrap_mem_fn(*const_cast<this_type *> (this),
                                                                   &this_type::forward_ctor),result,leaf_op,
                                               iaket1,iap11,iap21,iav11,iav21,eri);
            }

            this_type forward_ctor(implT* result1, const opT& leaf_op, const ctT& iaket1,
                                   const ctL& iap11, const ctL& iap21, const ctL& iav11, const ctL& iav21,
                                   const implT* eri1) {
                return this_type(result1,leaf_op,iaket1,iap11,iap21,iav11,iav21,eri1);
            }

            /// serialize this (needed for use in recursive_op)
            template <typename Archive> void serialize(const Archive& ar) {
                ar & iaket & eri & result & leaf_op & iap1 & iap2 & iav1 & iav2;
            }
        };

        /// assemble the function V*phi using V and phi given from the functor

        /// this function must have been constructed using the CompositeFunctorInterface.
        /// The interface provides one- and two-electron potentials, and the ket, which are
        /// assembled to give V*phi.
        /// @param[in]  leaf_op  operator to decide if a given node is a leaf node
        /// @param[in]  fence   global fence
        template<typename opT>
        void make_Vphi(const opT& leaf_op, const bool fence=true) {

            constexpr size_t LDIM=NDIM/2;
            MADNESS_CHECK_THROW(NDIM==LDIM*2,"make_Vphi only works for even dimensions");


            // keep the functor available, but remove it from the result
            // result will return false upon is_on_demand(), which is necessary for the
            // CoeffTracker to track the parent coeffs correctly for error_leaf_op
            std::shared_ptr< FunctionFunctorInterface<T,NDIM> > func2(this->get_functor());
            this->unset_functor();

            CompositeFunctorInterface<T,NDIM,LDIM>* func=
                    dynamic_cast<CompositeFunctorInterface<T,NDIM,LDIM>* >(&(*func2));
            MADNESS_ASSERT(func);

            // make sure everything is in place if no fence is requested
            if (fence) func->make_redundant(true);          // no-op if already redundant
            MADNESS_CHECK_THROW(func->check_redundant(),"make_Vphi requires redundant functions");

            // loop over all functions in the functor (either ket or particles)
            for (auto& ket : func->impl_ket_vector) {
                FunctionImpl<T,NDIM>* eri=func->impl_eri.get();
                FunctionImpl<T,LDIM>* v1=func->impl_m1.get();
                FunctionImpl<T,LDIM>* v2=func->impl_m2.get();
                FunctionImpl<T,LDIM>* p1=nullptr;
                FunctionImpl<T,LDIM>* p2=nullptr;
                make_Vphi_only(leaf_op,ket.get(),v1,v2,p1,p2,eri,false);
            }

            for (std::size_t i=0; i<func->impl_p1_vector.size(); ++i) {
                FunctionImpl<T,NDIM>* ket=nullptr;
                FunctionImpl<T,NDIM>* eri=func->impl_eri.get();
                FunctionImpl<T,LDIM>* v1=func->impl_m1.get();
                FunctionImpl<T,LDIM>* v2=func->impl_m2.get();
                FunctionImpl<T,LDIM>* p1=func->impl_p1_vector[i].get();
                FunctionImpl<T,LDIM>* p2=func->impl_p2_vector[i].get();
                make_Vphi_only(leaf_op,ket,v1,v2,p1,p2,eri,false);
            }

            // some post-processing:
            // - FunctionNode::accumulate() uses buffer -> add the buffer contents to the actual coefficients
            // - the operation constructs sum coefficients on all scales -> sum down to get a well-defined tree-state
            if (fence) {
                world.gop.fence();
                flo_unary_op_node_inplace(do_consolidate_buffer(get_tensor_args()),true);
                sum_down(true);
                set_tree_state(reconstructed);
            }


        }

        /// assemble the function V*phi using V and phi given from the functor

        /// this function must have been constructed using the CompositeFunctorInterface.
        /// The interface provides one- and two-electron potentials, and the ket, which are
        /// assembled to give V*phi.
        /// @param[in]  leaf_op  operator to decide if a given node is a leaf node
        /// @param[in]  fence   global fence
        template<typename opT, std::size_t LDIM>
        void make_Vphi_only(const opT& leaf_op, FunctionImpl<T,NDIM>* ket,
                            FunctionImpl<T,LDIM>* v1, FunctionImpl<T,LDIM>* v2,
                            FunctionImpl<T,LDIM>* p1, FunctionImpl<T,LDIM>* p2,
                            FunctionImpl<T,NDIM>* eri,
                            const bool fence=true) {

            // prepare the CoeffTracker
            CoeffTracker<T,NDIM> iaket(ket);
            CoeffTracker<T,LDIM> iap1(p1);
            CoeffTracker<T,LDIM> iap2(p2);
            CoeffTracker<T,LDIM> iav1(v1);
            CoeffTracker<T,LDIM> iav2(v2);

            // the operator making the coefficients
            typedef Vphi_op_NS<opT,LDIM> coeff_opT;
            coeff_opT coeff_op(this,leaf_op,iaket,iap1,iap2,iav1,iav2,eri);

            // this operator simply inserts the coeffs into this' tree
            typedef noop<T,NDIM> apply_opT;
            apply_opT apply_op;

            if (world.rank() == coeffs.owner(cdata.key0)) {
                woT::task(world.rank(), &implT:: template forward_traverse<coeff_opT,apply_opT>,
                          coeff_op, apply_op, cdata.key0);
            }

            set_tree_state(redundant_after_merge);
            if (fence) world.gop.fence();

        }

        /// Permute the dimensions of f according to map, result on this
        void mapdim(const implT& f, const std::vector<long>& map, bool fence);

        /// mirror the dimensions of f according to map, result on this
        void mirror(const implT& f, const std::vector<long>& mirror, bool fence);

        /// map and mirror the translation index and the coefficients, result on this

        /// first map the dimensions, the mirror!
        /// this = mirror(map(f))
        void map_and_mirror(const implT& f, const std::vector<long>& map,
        		const std::vector<long>& mirror, bool fence);

        /// take the average of two functions, similar to: this=0.5*(this+rhs)

        /// works in either basis and also in nonstandard form
        void average(const implT& rhs);

        /// change the tensor type of the coefficients in the FunctionNode

        /// @param[in]  targs   target tensor arguments (threshold and full/low rank)
        void change_tensor_type1(const TensorArgs& targs, bool fence);

        /// reduce the rank of the coefficients tensors

        /// @param[in]  targs   target tensor arguments (threshold and full/low rank)
        void reduce_rank(const double thresh, bool fence);


        /// remove all nodes with level higher than n
        void chop_at_level(const int n, const bool fence=true);

            /// compute norm of s and d coefficients for all nodes
        void compute_snorm_and_dnorm(bool fence=true);

        /// compute the norm of the wavelet coefficients
        struct do_compute_snorm_and_dnorm {
            typedef Range<typename dcT::iterator> rangeT;

            const FunctionCommonData<T,NDIM>& cdata;
            do_compute_snorm_and_dnorm(const FunctionCommonData<T,NDIM>& cdata) :
                    cdata(cdata) {}

            bool operator()(typename rangeT::iterator& it) const {
                auto& node=it->second;
                node.recompute_snorm_and_dnorm(cdata);
                return true;
            }
        };


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
        tensorT filter(const tensorT& s) const;

        coeffT filter(const coeffT& s) const;

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
        tensorT unfilter(const tensorT& s) const;

        coeffT unfilter(const coeffT& s) const;

        /// downsample the sum coefficients of level n+1 to sum coeffs on level n

        /// specialization of the filter method, will yield only the sum coefficients
        /// @param[in]  key key of level n
        /// @param[in]  v   vector of sum coefficients of level n+1
        /// @return     sum coefficients on level n in full tensor format
        tensorT downsample(const keyT& key, const std::vector< Future<coeffT > >& v) const;

        /// upsample the sum coefficients of level 1 to sum coeffs on level n+1

        /// specialization of the unfilter method, will transform only the sum coefficients
        /// @param[in]  key     key of level n+1
        /// @param[in]  coeff   sum coefficients of level n (does NOT belong to key!!)
        /// @return     sum     coefficients on level n+1
        coeffT upsample(const keyT& key, const coeffT& coeff) const;

        /// Projects old function into new basis (only in reconstructed form)
        void project(const implT& old, bool fence);

        struct true_refine_test {
            bool operator()(const implT* f, const keyT& key, const nodeT& t) const {
                return true;
            }
            template <typename Archive> void serialize(Archive& ar) {}
        };

        template <typename opT>
        void refine_op(const opT& op, const keyT& key) {
            // Must allow for someone already having autorefined the coeffs
            // and we get a write accessor just in case they are already executing
            typename dcT::accessor acc;
            const auto found = coeffs.find(acc,key);
            MADNESS_CHECK(found);
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
                    ss.reduce_rank(targs.thresh);
                    //                    coeffs.replace(child,nodeT(ss,-1.0,false).node_to_low_rank());
                    coeffs.replace(child,nodeT(ss,-1.0,false));
                    // Note value -1.0 for norm tree to indicate result of refinement
                }
            }
        }

        template <typename opT>
        void refine_spawn(const opT& op, const keyT& key) {
            nodeT& node = coeffs.find(key).get()->second;
            if (node.has_children()) {
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit)
                    woT::task(coeffs.owner(kit.key()), &implT:: template refine_spawn<opT>, op, kit.key(), TaskAttributes::hipri());
            }
            else {
                woT::task(coeffs.owner(key), &implT:: template refine_op<opT>, op, key);
            }
        }

        // Refine in real space according to local user-defined criterion
        template <typename opT>
        void refine(const opT& op, bool fence) {
            if (world.rank() == coeffs.owner(cdata.key0))
                woT::task(coeffs.owner(cdata.key0), &implT:: template refine_spawn<opT>, op, cdata.key0, TaskAttributes::hipri());
            if (fence)
                world.gop.fence();
        }

        bool exists_and_has_children(const keyT& key) const;

        bool exists_and_is_leaf(const keyT& key) const;


        void broaden_op(const keyT& key, const std::vector< Future <bool> >& v);

        // For each local node sets value of norm tree, snorm and dnorm to 0.0
        void zero_norm_tree();

        // Broaden tree
        void broaden(const array_of_bools<NDIM>& is_periodic, bool fence);

        /// sum all the contributions from all scales after applying an operator in mod-NS form
        void trickle_down(bool fence);

        /// sum all the contributions from all scales after applying an operator in mod-NS form

        /// cf reconstruct_op
        void trickle_down_op(const keyT& key, const coeffT& s);

        /// reconstruct this tree -- respects fence
        void reconstruct(bool fence);

        void change_tree_state(const TreeState finalstate, bool fence=true);

        // Invoked on node where key is local
        //        void reconstruct_op(const keyT& key, const tensorT& s);
        void reconstruct_op(const keyT& key, const coeffT& s, const bool accumulate_NS=true);

        /// compress the wave function

        /// after application there will be sum coefficients at the root level,
        /// and difference coefficients at all other levels; furthermore:
        /// @param[in] nonstandard	keep sum coeffs at all other levels, except leaves
        /// @param[in] keepleaves	keep sum coeffs (but no diff coeffs) at leaves
        /// @param[in] redundant    keep only sum coeffs at all levels, discard difference coeffs
//        void compress(bool nonstandard, bool keepleaves, bool redundant, bool fence);
        void compress(const TreeState newstate, bool fence);

        /// Invoked on node where key is local
        Future<std::pair<coeffT,double> > compress_spawn(const keyT& key, bool nonstandard, bool keepleaves,
        		bool redundant1);

        private:
        /// convert this to redundant, i.e. have sum coefficients on all levels
        void make_redundant(const bool fence);
        public:

        /// convert this from redundant to standard reconstructed form
        void undo_redundant(const bool fence);

        void remove_internal_coefficients(const bool fence);
        void remove_leaf_coefficients(const bool fence);


        /// compute for each FunctionNode the norm of the function inside that node
        void norm_tree(bool fence);

        double norm_tree_op(const keyT& key, const std::vector< Future<double> >& v);

        Future<double> norm_tree_spawn(const keyT& key);

        /// truncate using a tree in reconstructed form

        /// must be invoked where key is local
        Future<coeffT> truncate_reconstructed_spawn(const keyT& key, const double tol);

        /// given the sum coefficients of all children, truncate or not

        /// @return     new sum coefficients (empty if internal, not empty, if new leaf); might delete its children
        coeffT truncate_reconstructed_op(const keyT& key, const std::vector< Future<coeffT > >& v, const double tol);

        /// calculate the wavelet coefficients using the sum coefficients of all child nodes

        /// also compute the norm tree for all nodes
        /// @param[in] key 	this's key
        /// @param[in] v 	sum coefficients of the child nodes
        /// @param[in] nonstandard  keep the sum coefficients with the wavelet coefficients
        /// @param[in] redundant    keep only the sum coefficients, discard the wavelet coefficients
        /// @return 		the sum coefficients
        std::pair<coeffT,double> compress_op(const keyT& key, const std::vector< Future<std::pair<coeffT,double>> >& v, bool nonstandard);


        /// similar to compress_op, but insert only the sum coefficients in the tree

        /// also compute the norm tree for all nodes
        /// @param[in] key  this's key
        /// @param[in] v    sum coefficients of the child nodes
        /// @return         the sum coefficients
        std::pair<coeffT,double> make_redundant_op(const keyT& key,const std::vector< Future<std::pair<coeffT,double> > >& v);

        /// Changes non-standard compressed form to standard compressed form
        void standard(bool fence);

        /// Changes non-standard compressed form to standard compressed form
        struct do_standard {
            typedef Range<typename dcT::iterator> rangeT;

            // threshold for rank reduction / SVD truncation
            implT* impl;

            // constructor takes target precision
            do_standard() = default;
            do_standard(implT* impl) : impl(impl) {}

            //
            bool operator()(typename rangeT::iterator& it) const {

                const keyT& key = it->first;
                nodeT& node = it->second;
                if (key.level()> 0 && node.has_coeff()) {
                    if (node.has_children()) {
                        // Zero out scaling coeffs
                    	MADNESS_ASSERT(node.coeff().dim(0)==2*impl->get_k());
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

            do_op_args() = default;
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
        template <typename opT, typename R, size_t OPDIM>
        void do_apply_kernel(const opT* op, const Tensor<R>& c, const do_op_args<OPDIM>& args) {

            tensorT result = op->apply(args.key, args.d, c, args.tol/args.fac/args.cnorm);

            // Screen here to reduce communication cost of negligible data
            // and also to ensure we don't needlessly widen the tree when
            // applying the operator
            if (result.normf()> 0.3*args.tol/args.fac) {
                coeffs.task(args.dest, &nodeT::accumulate2, result, coeffs, args.dest, TaskAttributes::hipri());
                //woT::task(world.rank(),&implT::accumulate_timer,time,TaskAttributes::hipri());
                // UGLY BUT ADDED THE OPTIMIZATION BACK IN HERE EXPLICITLY/
                if (args.dest == world.rank()) {
                    coeffs.send(args.dest, &nodeT::accumulate, result, coeffs, args.dest);
                }
                else {
                    coeffs.task(args.dest, &nodeT::accumulate, result, coeffs, args.dest, TaskAttributes::hipri());
                }
            }
        }

        /// same as do_apply_kernel, but use full rank tensors as input and low rank tensors as output

        /// @param[in]  op      the operator working on our function
        /// @param[in]  c       full rank tensor holding the NS coefficients
        /// @param[in]  args    laziness holding norm of the coefficients, displacement, destination, ..
        /// @param[in]  apply_targs TensorArgs with tightened threshold for accumulation
        /// @return     nothing, but accumulate the result tensor into the destination node
        template <typename opT, typename R, size_t OPDIM>
        double do_apply_kernel2(const opT* op, const Tensor<R>& c, const do_op_args<OPDIM>& args,
                                const TensorArgs& apply_targs) {

            tensorT result_full = op->apply(args.key, args.d, c, args.tol/args.fac/args.cnorm);
            const double norm=result_full.normf();

            // Screen here to reduce communication cost of negligible data
            // and also to ensure we don't needlessly widen the tree when
            // applying the operator
            // OPTIMIZATION NEEDED HERE ... CHANGING THIS TO TASK NOT SEND REMOVED
            // BUILTIN OPTIMIZATION TO SHORTCIRCUIT MSG IF DATA IS LOCAL
            if (norm > 0.3*args.tol/args.fac) {

                small++;
                //double cpu0=cpu_time();
                coeffT result=coeffT(result_full,apply_targs);
                MADNESS_ASSERT(result.is_full_tensor() or result.is_svd_tensor());
                //double cpu1=cpu_time();
                //timer_lr_result.accumulate(cpu1-cpu0);

                coeffs.task(args.dest, &nodeT::accumulate, result, coeffs, args.dest, apply_targs,
                                                TaskAttributes::hipri());

                //woT::task(world.rank(),&implT::accumulate_timer,time,TaskAttributes::hipri());
            }
            return norm;
        }



        /// same as do_apply_kernel2, but use low rank tensors as input and low rank tensors as output

        /// @param[in]  op      the operator working on our function
        /// @param[in]  coeff   full rank tensor holding the NS coefficients
        /// @param[in]  args    laziness holding norm of the coefficients, displacement, destination, ..
        /// @param[in]  apply_targs TensorArgs with tightened threshold for accumulation
        /// @return     nothing, but accumulate the result tensor into the destination node
        template <typename opT, typename R, size_t OPDIM>
        double do_apply_kernel3(const opT* op, const GenTensor<R>& coeff, const do_op_args<OPDIM>& args,
                                const TensorArgs& apply_targs) {

            coeffT result;
            if (2*OPDIM==NDIM) result= op->apply2_lowdim(args.key, args.d, coeff,
                    args.tol/args.fac/args.cnorm, args.tol/args.fac);
            if (OPDIM==NDIM) result = op->apply2(args.key, args.d, coeff,
                    args.tol/args.fac/args.cnorm, args.tol/args.fac);

            const double result_norm=result.svd_normf();

            if (result_norm> 0.3*args.tol/args.fac) {
                small++;

                double cpu0=cpu_time();
                if (not result.is_of_tensortype(targs.tt)) result=result.convert(targs);
                double cpu1=cpu_time();
                timer_lr_result.accumulate(cpu1-cpu0);

                // accumulate also expects result in SVD form
                coeffs.task(args.dest, &nodeT::accumulate, result, coeffs, args.dest, apply_targs,
                                                TaskAttributes::hipri());
//                woT::task(world.rank(),&implT::accumulate_timer,time,TaskAttributes::hipri());

            }
            return result_norm;

        }

        // volume of n-dimensional sphere of radius R
        double vol_nsphere(int n, double R) {
            return std::pow(madness::constants::pi,n*0.5)*std::pow(R,n)/std::tgamma(1+0.5*n);
        }
        

        /// apply an operator on the coeffs c (at node key)

        /// the result is accumulated inplace to this's tree at various FunctionNodes
        /// @param[in] op	the operator to act on the source function
        /// @param[in] key	key of the source FunctionNode of f which is processed
        /// @param[in] c	coeffs of the FunctionNode of f which is processed
        template <typename opT, typename R>
        void do_apply(const opT* op, const keyT& key, const Tensor<R>& c) {
            PROFILE_MEMBER_FUNC(FunctionImpl);

          // working assumption here WAS that the operator is
          // isotropic and monotonically decreasing with distance
          // ... however, now we are using derivative Gaussian
          // expansions (and also non-cubic boxes) isotropic is
          // violated. While not strictly monotonically decreasing,
          // the derivative gaussian is still such that once it
          // becomes negligible we are in the asymptotic region.

          typedef typename opT::keyT opkeyT;
          constexpr auto opdim = opT::opdim;
          const opkeyT source = op->get_source_key(key);

          // Tuning here is based on observation that with
          // sufficiently high-order wavelet relative to the
          // precision, that only nearest neighbor boxes contribute,
          // whereas for low-order wavelets more neighbors will
          // contribute.  Sufficiently high is picked as
          // k>=2-log10(eps) which is our empirical rule for
          // efficiency/accuracy and code instrumentation has
          // previously indicated that (in 3D) just unit
          // displacements are invoked.  The error decays as R^-(k+1),
          // and the number of boxes increases as R^d.
          //
          // Fac is the expected number of contributions to a given
          // box, so the error permitted per contribution will be
          // tol/fac

          // radius of shell (nearest neighbor is diameter of 3 boxes, so radius=1.5)
          double radius = 1.5 + 0.33 * std::max(0.0, 2 - std::log10(thresh) -
                                                         k); // 0.33 was 0.5
          //double radius = 2.5;
          double fac = vol_nsphere(NDIM, radius);
          // previously fac=10.0 selected empirically constrained by qmprop

          double cnorm = c.normf();

          // BC handling:
          // - if operator is lattice-summed then treat this as nonperiodic (i.e. tell neighbor() to stay in simulation cell)
          // - if operator is NOT lattice-summed then obey BC (i.e. tell neighbor() to go outside the simulation cell along periodic dimensions)
          // - BUT user can force operator to treat its arguments as non-periodic (`op.set_domain_periodicity({true,true,true})`) so ... which dimensions of this function are treated as periodic by op?
          const array_of_bools<NDIM> this_is_treated_by_op_as_periodic =
              (op->particle() == 1)
                  ? array_of_bools<NDIM>{false}.or_front(
                        op->domain_is_periodic())
                  : array_of_bools<NDIM>{false}.or_back(
                        op->domain_is_periodic());

          const auto default_distance_squared = [&](const auto &displacement)
              -> std::uint64_t {
            return displacement.distsq_bc(op->lattice_summed());
          };
          const auto default_skip_predicate = [&](const auto &displacement)
              -> bool {
            return false;
          };
          const auto for_each = [&](const auto &displacements,
                                    const auto &distance_squared,
                                    const auto &skip_predicate) -> std::optional<std::uint64_t> {

            // used to screen estimated and actual contributions
            //const double tol = truncate_tol(thresh, key);
            //const double tol = 0.1*truncate_tol(thresh, key);
            const double tol = truncate_tol(thresh, key);

            // assume isotropic decaying kernel, screen in shell-wise fashion by
            // monitoring the decay of magnitude of contribution norms with the
            // distance ... as soon as we find a shell of displacements at least
            // one of each in simulation domain (see neighbor()) and
            // all in-domain shells produce negligible contributions, stop.
            // a displacement is negligible if ||op|| * ||c|| > tol / fac
            // where fac takes into account
            int nvalid = 1; // Counts #valid at each distance
            int nused = 1;  // Counts #used at each distance
            std::optional<std::uint64_t> distsq;

            // displacements to the kernel range boundary are typically same magnitude (modulo variation estimate the norm of the resulting contributions and skip all if one is too small
            // this
            if constexpr (std::is_same_v<std::decay_t<decltype(displacements)>,BoxSurfaceDisplacementRange<opdim>>) {
              const auto &probing_displacement =
                  displacements.probing_displacement();
              const double opnorm =
                  op->norm(key.level(), probing_displacement, source);
              if (cnorm * opnorm <= tol / fac) {
                return {};
              }
            }

            const auto disp_end = displacements.end();
            for (auto disp_it = displacements.begin(); disp_it != disp_end;
                 ++disp_it) {
              const auto &displacement = *disp_it;
              if (skip_predicate(displacement)) continue;

              keyT d;
              Key<NDIM - opdim> nullkey(key.level());
              MADNESS_ASSERT(op->particle() == 1 || op->particle() == 2);
              if (op->particle() == 1)
                d = displacement.merge_with(nullkey);
              else
                d = nullkey.merge_with(displacement);

              // shell-wise screening, assumes displacements are grouped into shells sorted so that operator decays with shell index N.B. lattice-summed decaying kernel is periodic (i.e. does decay w.r.t. r), so loop over shells of displacements sorted by distances modulated by periodicity (Key::distsq_bc)
              const uint64_t dsq = distance_squared(displacement);
              if (!distsq ||
                  dsq != *distsq) { // Moved to next shell of neighbors
                if (nvalid > 0 && nused == 0 && dsq > 1) {
                  // Have at least done the input box and all first
                  // nearest neighbors, and none of the last set
                  // of neighbors made significant contributions.  Thus,
                  // assuming monotonic decrease, we are done.
                  break;
                }
                nused = 0;
                nvalid = 0;
                distsq = dsq;
              }

              keyT dest = neighbor(key, d, this_is_treated_by_op_as_periodic);
              if (dest.is_valid()) {
                nvalid++;
                const double opnorm = op->norm(key.level(), displacement, source);

                if (cnorm * opnorm > tol / fac) {
                  tensorT result =
                      op->apply(source, displacement, c, tol / fac / cnorm);
                  if (result.normf() > 0.3 * tol / fac) {
                    if (coeffs.is_local(dest))
                      coeffs.send(dest, &nodeT::accumulate2, result, coeffs,
                                  dest);
                    else
                      coeffs.task(dest, &nodeT::accumulate2, result, coeffs,
                                  dest);
                    nused++;
                  }
                }
              }
            }

            return distsq;
          };

          // process "standard" displacements, screening assumes monotonic decay of the kernel
          // list of displacements sorted in order of increasing distance
          // N.B. if op is lattice-summed use periodic displacements, else use
          // non-periodic even if op treats any modes of this as periodic
          const std::vector<opkeyT> &disp = op->get_disp(key.level());
          const auto max_distsq_reached = for_each(disp, default_distance_squared, default_skip_predicate);

          // for range-restricted kernels displacements to the boundary of the kernel range also need to be included
          // N.B. hard range restriction will result in slow decay of operator matrix elements for the displacements
          // to the range boundary, should use soft restriction or sacrifice precision
          if (op->range_restricted() && key.level() >= 1) {

            std::array<std::optional<std::int64_t>, opdim> box_radius;
            std::array<std::optional<std::int64_t>, opdim> surface_thickness;
            auto &range = op->get_range();
            for (int d = 0; d != opdim; ++d) {
              if (range[d]) {
                box_radius[d] = range[d].N();
                surface_thickness[d] = range[d].finite_soft() ? 1 : 0;
              }
            }

            typename BoxSurfaceDisplacementRange<opdim>::Filter filter;
            // skip surface displacements that take us outside of the domain and/or were included in regular displacements
            // N.B. for lattice-summed axes the "filter" also maps the displacement back into the simulation cell
            if (max_distsq_reached)
              filter = BoxSurfaceDisplacementFilter<opdim>(/* domain_is_infinite= */ op->domain_is_periodic(), /* domain_is_periodic= */ op->lattice_summed(), range, default_distance_squared, *max_distsq_reached);

            // this range iterates over the entire surface layer(s), and provides a probing displacement that can be used to screen out the entire box
            auto opkey = op->particle() == 1 ? key.template extract_front<opdim>() : key.template extract_front<opdim>();
            BoxSurfaceDisplacementRange<opdim>
                range_boundary_face_displacements(opkey, box_radius,
                                                  surface_thickness,
                                                  op->lattice_summed(),  // along lattice-summed axes treat the box as periodic, make displacements to one side of the box
                                                  filter);
            for_each(
                range_boundary_face_displacements,
                // surface displacements are not screened, all are included
                [](const auto &displacement) -> std::uint64_t { return 0; },
                default_skip_predicate);
          }
        }


        /// apply an operator on f to return this
        template <typename opT, typename R>
        void apply(opT& op, const FunctionImpl<R,NDIM>& f, bool fence) {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            MADNESS_ASSERT(!op.modified());
            typename dcT::const_iterator end = f.coeffs.end();
            for (typename dcT::const_iterator it=f.coeffs.begin(); it!=end; ++it) {
                // looping through all the coefficients in the source
                const keyT& key = it->first;
                const FunctionNode<R,NDIM>& node = it->second;
                if (node.has_coeff()) {
                    if (node.coeff().dim(0) != k /* i.e. not a leaf */ || op.doleaves) {
                        ProcessID p = FunctionDefaults<NDIM>::get_apply_randomize() ? world.random_proc() : coeffs.owner(key);
//                        woT::task(p, &implT:: template do_apply<opT,R>, &op, key, node.coeff()); //.full_tensor_copy() ????? why copy ????
                        woT::task(p, &implT:: template do_apply<opT,R>, &op, key, node.coeff().reconstruct_tensor());
                    }
                }
            }
            if (fence)
                world.gop.fence();

            set_tree_state(nonstandard_after_apply);
//            this->compressed=true;
//            this->nonstandard=true;
//            this->redundant=false;

        }



        /// apply an operator on the coeffs c (at node key)

        /// invoked by result; the result is accumulated inplace to this's tree at various FunctionNodes
        /// @param[in] op     the operator to act on the source function
        /// @param[in] key    key of the source FunctionNode of f which is processed (see "source")
        /// @param[in] coeff  coeffs of FunctionNode being processed
        /// @param[in] do_kernel	true: do the 0-disp only; false: do everything but the kernel
        /// @return	   max norm, and will modify or include new nodes in this' tree
        template <typename opT, typename R>
        double do_apply_directed_screening(const opT* op, const keyT& key, const coeffT& coeff,
                                           const bool& do_kernel) {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            // insert timer here
            typedef typename opT::keyT opkeyT;

            // screening: contains all displacement keys that had small result norms
            std::list<opkeyT> blacklist;

            constexpr auto opdim=opT::opdim;
            Key<NDIM-opdim> nullkey(key.level());

            // source is that part of key that corresponds to those dimensions being processed
            const opkeyT source=op->get_source_key(key);

            const double tol = truncate_tol(thresh, key);

            // fac is the root of the number of contributing neighbors (1st shell)
            double fac=std::pow(3,NDIM*0.5);
            double cnorm = coeff.normf();

            // for accumulation: keep slightly tighter TensorArgs
            TensorArgs apply_targs(targs);
            apply_targs.thresh=tol/fac*0.03;

            double maxnorm=0.0;

            // for the kernel it may be more efficient to do the convolution in full rank
            tensorT coeff_full;
            // for partial application (exchange operator) it's more efficient to
            // do SVD tensors instead of tensortrains, because addition in apply
            // can be done in full form for the specific particle
            coeffT coeff_SVD=coeff.convert(TensorArgs(-1.0,TT_2D));
#ifdef HAVE_GENTENSOR
            coeff_SVD.get_svdtensor().orthonormalize(tol*GenTensor<T>::fac_reduce());
#endif

            // list of displacements sorted in order of increasing distance
            // N.B. if op is lattice-summed gives periodic displacements, else uses
            // non-periodic even if op treats any modes of this as periodic
            const std::vector<opkeyT>& disp = Displacements<opdim>().get_disp(key.level(), op->lattice_summed());

            for (typename std::vector<opkeyT>::const_iterator it=disp.begin(); it != disp.end(); ++it) {
                const opkeyT& d = *it;

                const int shell=d.distsq_bc(op->lattice_summed());
                if (do_kernel and (shell>0)) break;
                if ((not do_kernel) and (shell==0)) continue;

                keyT disp1;
                if (op->particle()==1) disp1=it->merge_with(nullkey);
                else if (op->particle()==2) disp1=nullkey.merge_with(*it);
                else {
                    MADNESS_EXCEPTION("confused particle in operator??",1);
                }

                keyT dest = neighbor_in_volume(key, disp1);

                if (not dest.is_valid()) continue;

                // directed screening
                // working assumption here is that the operator is isotropic and
                // monotonically decreasing with distance
                bool screened=false;
                typename std::list<opkeyT>::const_iterator it2;
                for (it2=blacklist.begin(); it2!=blacklist.end(); it2++) {
                    if (d.is_farther_out_than(*it2)) {
                        screened=true;
                        break;
                    }
                }
                if (not screened) {

                    double opnorm = op->norm(key.level(), d, source);
                    double norm=0.0;

                    if (cnorm*opnorm> tol/fac) {

                        double cost_ratio=op->estimate_costs(source, d, coeff_SVD, tol/fac/cnorm, tol/fac);
                        //                        cost_ratio=1.5;     // force low rank
                        //                        cost_ratio=0.5;     // force full rank

                        if (cost_ratio>0.0) {

                            do_op_args<opdim> args(source, d, dest, tol, fac, cnorm);
                            norm=0.0;
                            if (cost_ratio<1.0) {
                                if (not coeff_full.has_data()) coeff_full=coeff.full_tensor_copy();
                                norm=do_apply_kernel2(op, coeff_full,args,apply_targs);
                            } else {
                                if (2*opdim==NDIM) {    // apply operator on one particle only
                                    norm=do_apply_kernel3(op,coeff_SVD,args,apply_targs);
                                } else {
                                    norm=do_apply_kernel3(op,coeff,args,apply_targs);
                                }
                            }
                            maxnorm=std::max(norm,maxnorm);
                        }

                    } else if (shell >= 12) {
                        break; // Assumes monotonic decay beyond nearest neighbor
                    }
                    if (norm<0.3*tol/fac) blacklist.push_back(d);
                }
            }
            return maxnorm;
        }


        /// similar to apply, but for low rank coeffs
        template <typename opT, typename R>
        void apply_source_driven(opT& op, const FunctionImpl<R,NDIM>& f, bool fence) {
            PROFILE_MEMBER_FUNC(FunctionImpl);

            MADNESS_ASSERT(not op.modified());
            // looping through all the coefficients of the source f
            typename dcT::const_iterator end = f.get_coeffs().end();
            for (typename dcT::const_iterator it=f.get_coeffs().begin(); it!=end; ++it) {

                const keyT& key = it->first;
                const coeffT& coeff = it->second.coeff();

                if (coeff.has_data() and (coeff.rank()!=0)) {
                    ProcessID p = FunctionDefaults<NDIM>::get_apply_randomize() ? world.random_proc() : coeffs.owner(key);
                    woT::task(p, &implT:: template do_apply_directed_screening<opT,R>, &op, key, coeff, true);
                    woT::task(p, &implT:: template do_apply_directed_screening<opT,R>, &op, key, coeff, false);
                }
            }
            if (fence) world.gop.fence();
            set_tree_state(TreeState::nonstandard_after_apply);
        }

        /// after apply we need to do some cleanup;

        /// forces fence
        double finalize_apply();

        /// after summing up we need to do some cleanup;

        /// forces fence
        void finalize_sum();

        /// traverse a non-existing tree, make its coeffs and apply an operator

        /// invoked by result
        /// here we use the fact that the hi-dim NS coefficients on all scales are exactly
        /// the outer product of the underlying low-dim functions (also in NS form),
        /// so we don't need to construct the full hi-dim tree and then turn it into NS form.
        /// @param[in]	apply_op the operator acting on the NS tree
        /// @param[in]	fimpl    the funcimpl of the function of particle 1
        /// @param[in]	gimpl    the funcimpl of the function of particle 2
        template<typename opT, std::size_t LDIM>
        void recursive_apply(opT& apply_op, const FunctionImpl<T,LDIM>* fimpl,
                             const FunctionImpl<T,LDIM>* gimpl, const bool fence) {

            //print("IN RECUR2");
            const keyT& key0=cdata.key0;

            if (world.rank() == coeffs.owner(key0)) {

                CoeffTracker<T,LDIM> ff(fimpl);
                CoeffTracker<T,LDIM> gg(gimpl);

                typedef recursive_apply_op<opT,LDIM> coeff_opT;
                coeff_opT coeff_op(this,ff,gg,&apply_op);

                typedef noop<T,NDIM> apply_opT;
                apply_opT apply_op;

                ProcessID p= coeffs.owner(key0);
                woT::task(p, &implT:: template forward_traverse<coeff_opT,apply_opT>, coeff_op, apply_op, key0);

            }
            if (fence) world.gop.fence();
            set_tree_state(TreeState::nonstandard_after_apply);
        }

        /// recursive part of recursive_apply
        template<typename opT, std::size_t LDIM>
        struct recursive_apply_op {
            bool randomize() const {return true;}

            typedef recursive_apply_op<opT,LDIM> this_type;

            implT* result;
            CoeffTracker<T,LDIM> iaf;
            CoeffTracker<T,LDIM> iag;
            opT* apply_op;

            // ctor
            recursive_apply_op() = default;
            recursive_apply_op(implT* result,
                               const CoeffTracker<T,LDIM>& iaf, const CoeffTracker<T,LDIM>& iag,
                               const opT* apply_op) : result(result), iaf(iaf), iag(iag), apply_op(apply_op)
            {
                MADNESS_ASSERT(LDIM+LDIM==NDIM);
            }
            recursive_apply_op(const recursive_apply_op& other) : result(other.result), iaf(other.iaf),
                                                                  iag(other.iag), apply_op(other.apply_op) {}


            /// make the NS-coefficients and send off the application of the operator

            /// @return		a Future<bool,coeffT>(is_leaf,coeffT())
            std::pair<bool,coeffT> operator()(const Key<NDIM>& key) const {

                //            	World& world=result->world;
                // break key into particles (these are the child keys, with datum1/2 come the parent keys)
                Key<LDIM> key1,key2;
                key.break_apart(key1,key2);

                // the lo-dim functions should be in full tensor form
                const tensorT fcoeff=iaf.coeff(key1).full_tensor();
                const tensorT gcoeff=iag.coeff(key2).full_tensor();

                // would this be a leaf node? If so, then its sum coeffs have already been
                // processed by the parent node's wavelet coeffs. Therefore we won't
                // process it any more.
                hartree_leaf_op<T,NDIM> leaf_op(result,result->get_k());
                bool is_leaf=leaf_op(key,fcoeff,gcoeff);

                if (not is_leaf) {
                    // new coeffs are simply the hartree/kronecker/outer product --
                    const std::vector<Slice>& s0=iaf.get_impl()->cdata.s0;
                    const coeffT coeff = (apply_op->modified())
                        ? outer(copy(fcoeff(s0)),copy(gcoeff(s0)),result->targs)
                        : outer(fcoeff,gcoeff,result->targs);

                    // now send off the application
                    tensorT coeff_full;
                    ProcessID p=result->world.rank();
                    double norm0=result->do_apply_directed_screening<opT,T>(apply_op, key, coeff, true);

                    result->task(p,&implT:: template do_apply_directed_screening<opT,T>,
                                 apply_op,key,coeff,false);

                    return finalize(norm0,key,coeff);

                } else {
                    return std::pair<bool,coeffT> (is_leaf,coeffT());
                }
            }

            /// sole purpose is to wait for the kernel norm, wrap it and send it back to caller
            std::pair<bool,coeffT> finalize(const double kernel_norm, const keyT& key,
                                            const coeffT& coeff) const {
            	const double thresh=result->get_thresh()*0.1;
            	bool is_leaf=(kernel_norm<result->truncate_tol(thresh,key));
            	if (key.level()<2) is_leaf=false;
            	return std::pair<bool,coeffT> (is_leaf,coeff);
            }


            this_type make_child(const keyT& child) const {

                // break key into particles
                Key<LDIM> key1, key2;
                child.break_apart(key1,key2);

                return this_type(result,iaf.make_child(key1),iag.make_child(key2),apply_op);
            }

            Future<this_type> activate() const {
            	Future<CoeffTracker<T,LDIM> > f1=iaf.activate();
            	Future<CoeffTracker<T,LDIM> > g1=iag.activate();
                return result->world.taskq.add(detail::wrap_mem_fn(*const_cast<this_type *> (this),
                                               &this_type::forward_ctor),result,f1,g1,apply_op);
            }

            this_type forward_ctor(implT* r, const CoeffTracker<T,LDIM>& f1, const CoeffTracker<T,LDIM>& g1,
                                   const opT* apply_op1) {
            	return this_type(r,f1,g1,apply_op1);
            }

            template <typename Archive> void serialize(const Archive& ar) {
                ar & result & iaf & iag & apply_op;
            }
        };

        /// traverse an existing tree and apply an operator

        /// invoked by result
        /// @param[in]	apply_op the operator acting on the NS tree
        /// @param[in]	fimpl    the funcimpl of the source function
        /// @param[in]	rimpl    a dummy function for recursive_op to insert data
        template<typename opT>
        void recursive_apply(opT& apply_op, const implT* fimpl, implT* rimpl, const bool fence) {

            print("IN RECUR1");

            const keyT& key0=cdata.key0;

            if (world.rank() == coeffs.owner(key0)) {

                typedef recursive_apply_op2<opT> coeff_opT;
                coeff_opT coeff_op(this,fimpl,&apply_op);

                typedef noop<T,NDIM> apply_opT;
                apply_opT apply_op;

                woT::task(world.rank(), &implT:: template forward_traverse<coeff_opT,apply_opT>,
                          coeff_op, apply_op, cdata.key0);

            }
            if (fence) world.gop.fence();
            set_tree_state(TreeState::nonstandard_after_apply);
        }

        /// recursive part of recursive_apply
        template<typename opT>
        struct recursive_apply_op2 {
            bool randomize() const {return true;}

            typedef recursive_apply_op2<opT> this_type;
            typedef CoeffTracker<T,NDIM> ctT;
            typedef std::pair<bool,coeffT> argT;

            mutable implT* result;
            ctT iaf;			/// need this for randomization
            const opT* apply_op;

            // ctor
            recursive_apply_op2() = default;
            recursive_apply_op2(implT* result, const ctT& iaf, const opT* apply_op)
            	: result(result), iaf(iaf), apply_op(apply_op) {}

            recursive_apply_op2(const recursive_apply_op2& other) : result(other.result),
                                                                    iaf(other.iaf), apply_op(other.apply_op) {}


            /// send off the application of the operator

            /// the first (core) neighbor (ie. the box itself) is processed
            /// immediately, all other ones are shoved into the taskq
            /// @return		a pair<bool,coeffT>(is_leaf,coeffT())
            argT operator()(const Key<NDIM>& key) const {

            	const coeffT& coeff=iaf.coeff();

                if (coeff.has_data()) {

                    // now send off the application for all neighbor boxes
                    ProcessID p=result->world.rank();
                    result->task(p,&implT:: template do_apply_directed_screening<opT,T>,
                                 apply_op, key, coeff, false);

                    // process the core box
                    double norm0=result->do_apply_directed_screening<opT,T>(apply_op,key,coeff,true);

                    if (iaf.is_leaf()) return argT(true,coeff);
                    return finalize(norm0,key,coeff,result);

                } else {
                    const bool is_leaf=true;
                    return argT(is_leaf,coeffT());
                }
            }

            /// sole purpose is to wait for the kernel norm, wrap it and send it back to caller
            argT finalize(const double kernel_norm, const keyT& key,
                          const coeffT& coeff, const implT* r) const {
            	const double thresh=r->get_thresh()*0.1;
            	bool is_leaf=(kernel_norm<r->truncate_tol(thresh,key));
            	if (key.level()<2) is_leaf=false;
            	return argT(is_leaf,coeff);
            }


            this_type make_child(const keyT& child) const {
                return this_type(result,iaf.make_child(child),apply_op);
            }

            /// retrieve the coefficients (parent coeffs might be remote)
            Future<this_type> activate() const {
            	Future<ctT> f1=iaf.activate();

//                Future<ctL> g1=g.activate();
//                return h->world.taskq.add(detail::wrap_mem_fn(*const_cast<this_type *> (this),
//                                          &this_type::forward_ctor),h,f1,g1,particle);

                return result->world.taskq.add(detail::wrap_mem_fn(*const_cast<this_type *> (this),
                                               &this_type::forward_ctor),result,f1,apply_op);
            }

            /// taskq-compatible ctor
            this_type forward_ctor(implT* result1, const ctT& iaf1, const opT* apply_op1) {
            	return this_type(result1,iaf1,apply_op1);
            }

            template <typename Archive> void serialize(const Archive& ar) {
                ar & result & iaf & apply_op;
            }
        };

        /// Returns the square of the error norm in the box labeled by key

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
            do_err_box() = default;

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
                MADNESS_EXCEPTION("not yet", 1);
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
                MADNESS_EXCEPTION("NOT IMPLEMENTED", 1);
            }
        };


        /// Returns the square of the local norm ... no comms
        double norm2sq_local() const;

        /// compute the inner product of this range with other
        template<typename R>
        struct do_inner_local {
            const FunctionImpl<R,NDIM>* other;
            bool leaves_only;
            typedef TENSOR_RESULT_TYPE(T,R) resultT;

            do_inner_local(const FunctionImpl<R,NDIM>* other, const bool leaves_only)
            	: other(other), leaves_only(leaves_only) {}
            resultT operator()(typename dcT::const_iterator& it) const {

            	TENSOR_RESULT_TYPE(T,R) sum=0.0;
            	const keyT& key=it->first;
                const nodeT& fnode = it->second;
                if (fnode.has_coeff()) {
                    if (other->coeffs.probe(it->first)) {
                        const FunctionNode<R,NDIM>& gnode = other->coeffs.find(key).get()->second;
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
                return sum;
            }

            resultT operator()(resultT a, resultT b) const {
                return (a+b);
            }

            template <typename Archive> void serialize(const Archive& ar) {
                MADNESS_EXCEPTION("NOT IMPLEMENTED", 1);
            }
        };

        /// Returns the inner product ASSUMING same distribution

        /// handles compressed and redundant form
        template <typename R>
        TENSOR_RESULT_TYPE(T,R) inner_local(const FunctionImpl<R,NDIM>& g) const {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            typedef Range<typename dcT::const_iterator> rangeT;
            typedef TENSOR_RESULT_TYPE(T,R) resultT;

            // make sure the states of the trees are consistent
            MADNESS_ASSERT(this->is_redundant()==g.is_redundant());
            bool leaves_only=(this->is_redundant());
            return world.taskq.reduce<resultT,rangeT,do_inner_local<R> >
                (rangeT(coeffs.begin(),coeffs.end()),do_inner_local<R>(&g, leaves_only));
        }


        /// compute the inner product of this range with other
        template<typename R>
        struct do_inner_local_on_demand {
            const FunctionImpl<T,NDIM>* bra;
            const FunctionImpl<R,NDIM>* ket;
            bool leaves_only=true;
            typedef TENSOR_RESULT_TYPE(T,R) resultT;

            do_inner_local_on_demand(const FunctionImpl<T,NDIM>* bra, const FunctionImpl<R,NDIM>* ket,
                                     const bool leaves_only=true)
                    : bra(bra), ket(ket), leaves_only(leaves_only) {}
            resultT operator()(typename dcT::const_iterator& it) const {

                constexpr std::size_t LDIM=std::max(NDIM/2,std::size_t(1));

                const keyT& key=it->first;
                const nodeT& fnode = it->second;
                if (not fnode.has_coeff()) return resultT(0.0); // probably internal nodes

                // assuming all boxes (esp the low-dim ones) are local, i.e. the functions are replicated
                auto find_valid_parent = [](auto& key, auto& impl, auto&& find_valid_parent) {
                    MADNESS_CHECK(impl->get_coeffs().owner(key)==impl->world.rank()); // make sure everything is local!
                    if (impl->get_coeffs().probe(key)) return key;
                    auto parentkey=key.parent();
                    return find_valid_parent(parentkey, impl, find_valid_parent);
                };

                // returns coefficients, empty if no functor present
                auto get_coeff = [&find_valid_parent](const auto& key, const auto& v_impl) {
                    if ((v_impl.size()>0) and v_impl.front().get()) {
                        auto impl=v_impl.front();

//                    bool have_impl=impl.get();
//                    if (have_impl) {
                        auto parentkey = find_valid_parent(key, impl, find_valid_parent);
                        MADNESS_CHECK(impl->get_coeffs().probe(parentkey));
                        typename decltype(impl->coeffs)::accessor acc;
                        impl->get_coeffs().find(acc,parentkey);
                        auto parentcoeff=acc->second.coeff();
                        auto coeff=impl->parent_to_child(parentcoeff, parentkey, key);
                        return coeff;
                    } else {
                        // get type of vector elements
                        typedef typename std::decay_t<decltype(v_impl)>::value_type::element_type::typeT S;
//                        typedef typename std::decay_t<decltype(v_impl)>::value_type S;
                        return GenTensor<S>();
//                        return GenTensor<typename std::decay_t<decltype(*impl)>::typeT>();
                    }
                };

                auto make_vector = [](auto& arg) {
                    return std::vector<std::decay_t<decltype(arg)>>(1,arg);
                };


                Key<LDIM> key1,key2;
                key.break_apart(key1,key2);

                auto func=dynamic_cast<CompositeFunctorInterface<R,NDIM,LDIM>* >(ket->functor.get());
                MADNESS_ASSERT(func);

                MADNESS_CHECK_THROW(func->impl_ket_vector.size()==0 or func->impl_ket_vector.size()==1,
                                    "only one ket function supported in inner_on_demand");
                MADNESS_CHECK_THROW(func->impl_p1_vector.size()==0 or func->impl_p1_vector.size()==1,
                                    "only one p1 function supported in inner_on_demand");
                MADNESS_CHECK_THROW(func->impl_p2_vector.size()==0 or func->impl_p2_vector.size()==1,
                                    "only one p2 function supported in inner_on_demand");
                auto coeff_bra=fnode.coeff();
                auto coeff_ket=get_coeff(key,func->impl_ket_vector);
                auto coeff_v1=get_coeff(key1,make_vector(func->impl_m1));
                auto coeff_v2=get_coeff(key2,make_vector(func->impl_m2));
                auto coeff_p1=get_coeff(key1,func->impl_p1_vector);
                auto coeff_p2=get_coeff(key2,func->impl_p2_vector);

                // construct |ket(1,2)> or |p(1)p(2)> or |p(1)p(2) ket(1,2)>
                double error=0.0;
                if (coeff_ket.has_data() and coeff_p1.has_data()) {
                    pointwise_multiplier<LDIM> pm(key,coeff_ket);
                    coeff_ket=pm(key,outer(coeff_p1,coeff_p2,TensorArgs(TT_FULL,-1.0)).full_tensor());
                    error+=pm.error;
                } else if (coeff_ket.has_data() or coeff_p1.has_data()) {
                    coeff_ket = (coeff_ket.has_data()) ? coeff_ket : outer(coeff_p1,coeff_p2);
                } else { // not ket and no p1p2
                    MADNESS_EXCEPTION("confused ket/p1p2 in do_inner_local_on_demand",1);
                }

                // construct (v(1) + v(2)) |ket(1,2)>
                coeffT v1v2ket;
                if (coeff_v1.has_data()) {
                    pointwise_multiplier<LDIM> pm(key,coeff_ket);
                    v1v2ket = pm(key,coeff_v1.full_tensor(), 1);
                    error+=pm.error;
                    v1v2ket+= pm(key,coeff_v2.full_tensor(), 2);
                    error+=pm.error;
                } else {
                    v1v2ket = coeff_ket;
                }

                resultT result;
                if (func->impl_eri) {         // project bra*ket onto eri, avoid multiplication with eri
                    MADNESS_CHECK(func->impl_eri->get_functor()->provides_coeff());
                    coeffT coeff_eri=func->impl_eri->get_functor()->coeff(key).full_tensor();
                    pointwise_multiplier<LDIM> pm(key,v1v2ket);
                    tensorT braket=pm(key,coeff_bra.full_tensor_copy().conj());
                    error+=pm.error;
                    if (error>1.e-3) print("error in key",key,error);
                    result=coeff_eri.full_tensor().trace(braket);

                } else {                            // no eri, project ket onto bra
                    result=coeff_bra.full_tensor_copy().trace_conj(v1v2ket.full_tensor_copy());
                }
                return result;
            }

            resultT operator()(resultT a, resultT b) const {
                return (a+b);
            }

            template <typename Archive> void serialize(const Archive& ar) {
                MADNESS_EXCEPTION("NOT IMPLEMENTED", 1);
            }
        };

        /// Returns the inner product of this with function g constructed on-the-fly

        /// the leaf boxes of this' MRA tree defines the inner product
        template <typename R>
        TENSOR_RESULT_TYPE(T,R) inner_local_on_demand(const FunctionImpl<R,NDIM>& gimpl) const {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            MADNESS_CHECK(this->is_reconstructed());

            typedef Range<typename dcT::const_iterator> rangeT;
            rangeT range(coeffs.begin(), coeffs.end());
            return world.taskq.reduce<T, rangeT, do_inner_local_on_demand<T>>(range,
                               do_inner_local_on_demand<R>(this, &gimpl));
        }

        /// compute the inner product of this range with other
        template<typename R>
        struct do_dot_local {
            const FunctionImpl<R,NDIM>* other;
            bool leaves_only;
            typedef TENSOR_RESULT_TYPE(T,R) resultT;

            do_dot_local(const FunctionImpl<R,NDIM>* other, const bool leaves_only)
            	: other(other), leaves_only(leaves_only) {}
            resultT operator()(typename dcT::const_iterator& it) const {

            	TENSOR_RESULT_TYPE(T,R) sum=0.0;
            	const keyT& key=it->first;
                const nodeT& fnode = it->second;
                if (fnode.has_coeff()) {
                    if (other->coeffs.probe(it->first)) {
                        const FunctionNode<R,NDIM>& gnode = other->coeffs.find(key).get()->second;
                        if (gnode.has_coeff()) {
                            if (gnode.coeff().dim(0) != fnode.coeff().dim(0)) {
                                madness::print("DOT", it->first, gnode.coeff().dim(0),fnode.coeff().dim(0));
                                MADNESS_EXCEPTION("functions have different k or compress/reconstruct error", 0);
                            }
                            if (leaves_only) {
                                if (gnode.is_leaf() or fnode.is_leaf()) {
                                    sum += fnode.coeff().full_tensor().trace(gnode.coeff().full_tensor());
                                }
                            } else {
                                sum += fnode.coeff().full_tensor().trace(gnode.coeff().full_tensor());
                            }
                        }
                    }
                }
                return sum;
            }

            resultT operator()(resultT a, resultT b) const {
                return (a+b);
            }

            template <typename Archive> void serialize(const Archive& ar) {
                MADNESS_EXCEPTION("NOT IMPLEMENTED", 1);
            }
        };

        /// Returns the dot product ASSUMING same distribution

        /// handles compressed and redundant form
        template <typename R>
        TENSOR_RESULT_TYPE(T,R) dot_local(const FunctionImpl<R,NDIM>& g) const {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            typedef Range<typename dcT::const_iterator> rangeT;
            typedef TENSOR_RESULT_TYPE(T,R) resultT;

            // make sure the states of the trees are consistent
            MADNESS_ASSERT(this->is_redundant()==g.is_redundant());
            bool leaves_only=(this->is_redundant());
            return world.taskq.reduce<resultT,rangeT,do_dot_local<R> >
                (rangeT(coeffs.begin(),coeffs.end()),do_dot_local<R>(&g, leaves_only));
        }

        /// Type of the entry in the map returned by make_key_vec_map
        typedef std::vector< std::pair<int,const coeffT*> > mapvecT;

        /// Type of the map returned by make_key_vec_map
        typedef ConcurrentHashMap< keyT, mapvecT > mapT;

        /// Adds keys to union of local keys with specified index
        void add_keys_to_map(mapT* map, int index) const {
            typename dcT::const_iterator end = coeffs.end();
            for (typename dcT::const_iterator it=coeffs.begin(); it!=end; ++it) {
                typename mapT::accessor acc;
                const keyT& key = it->first;
                const FunctionNode<T,NDIM>& node = it->second;
                if (node.has_coeff()) {
                    [[maybe_unused]] auto inserted = map->insert(acc,key);
                    acc->second.push_back(std::make_pair(index,&(node.coeff())));
                }
            }
        }

        /// Returns map of union of local keys to vector of indexes of functions containing that key

        /// Local concurrency and synchronization only; no communication
        static
        mapT
        make_key_vec_map(const std::vector<const FunctionImpl<T,NDIM>*>& v) {
            mapT map(100000);
            // This loop must be parallelized
            for (unsigned int i=0; i<v.size(); i++) {
                //v[i]->add_keys_to_map(&map,i);
                v[i]->world.taskq.add(*(v[i]), &FunctionImpl<T,NDIM>::add_keys_to_map, &map, int(i));
            }
            if (v.size()) v[0]->world.taskq.fence();
            return map;
        }

#if 0
// Original
        template <typename R>
        static void do_inner_localX(const typename mapT::iterator lstart,
                                    const typename mapT::iterator lend,
                                    typename FunctionImpl<R,NDIM>::mapT* rmap_ptr,
                                    const bool sym,
                                    Tensor< TENSOR_RESULT_TYPE(T,R) >* result_ptr,
                                    Mutex* mutex) {
            Tensor< TENSOR_RESULT_TYPE(T,R) >& result = *result_ptr;
            Tensor< TENSOR_RESULT_TYPE(T,R) > r(result.dim(0),result.dim(1));
            for (typename mapT::iterator lit=lstart; lit!=lend; ++lit) {
                const keyT& key = lit->first;
                typename FunctionImpl<R,NDIM>::mapT::iterator rit=rmap_ptr->find(key);
                if (rit != rmap_ptr->end()) {
                    const mapvecT& leftv = lit->second;
                    const typename FunctionImpl<R,NDIM>::mapvecT& rightv =rit->second;
                    const int nleft = leftv.size();
                    const int nright= rightv.size();

                    for (int iv=0; iv<nleft; iv++) {
                        const int i = leftv[iv].first;
                        const GenTensor<T>* iptr = leftv[iv].second;

                        for (int jv=0; jv<nright; jv++) {
                            const int j = rightv[jv].first;
                            const GenTensor<R>* jptr = rightv[jv].second;

                            if (!sym || (sym && i<=j))
                                r(i,j) += iptr->trace_conj(*jptr);
                        }
                    }
                }
            }
            mutex->lock();
            result += r;
            mutex->unlock();
        }
#else
       template <typename R>
       static void do_inner_localX(const typename mapT::iterator lstart,
                                   const typename mapT::iterator lend,
                                   typename FunctionImpl<R,NDIM>::mapT* rmap_ptr,
                                   const bool sym,
                                   Tensor< TENSOR_RESULT_TYPE(T,R) >* result_ptr,
                                   Mutex* mutex) {
           Tensor< TENSOR_RESULT_TYPE(T,R) >& result = *result_ptr;
           //Tensor< TENSOR_RESULT_TYPE(T,R) > r(result.dim(0),result.dim(1));
           for (typename mapT::iterator lit=lstart; lit!=lend; ++lit) {
               const keyT& key = lit->first;
               typename FunctionImpl<R,NDIM>::mapT::iterator rit=rmap_ptr->find(key);
               if (rit != rmap_ptr->end()) {
                   const mapvecT& leftv = lit->second;
                   const typename FunctionImpl<R,NDIM>::mapvecT& rightv =rit->second;
                   const size_t nleft = leftv.size();
                   const size_t nright= rightv.size();

                   unsigned int size = leftv[0].second->size();
                   Tensor<T> Left(nleft, size);
                   Tensor<R> Right(nright, size);
                   Tensor< TENSOR_RESULT_TYPE(T,R)> r(nleft, nright);
                   for(unsigned int iv = 0; iv < nleft; ++iv) Left(iv,_) = (*(leftv[iv].second)).full_tensor();
                   for(unsigned int jv = 0; jv < nright; ++jv) Right(jv,_) = (*(rightv[jv].second)).full_tensor();
                   // call mxmT from mxm.h in tensor
                   if(TensorTypeData<T>::iscomplex) Left = Left.conj();  // Should handle complex case and leave real case alone
                   mxmT(nleft, nright, size, r.ptr(), Left.ptr(), Right.ptr());
                   mutex->lock();
                   for(unsigned int iv = 0; iv < nleft; ++iv) {
                       const int i = leftv[iv].first;
                       for(unsigned int jv = 0; jv < nright; ++jv) {
                         const int j = rightv[jv].first;
                         if (!sym || (sym && i<=j)) result(i,j) += r(iv,jv);
                       }
                   }
                   mutex->unlock();
               }
           }
       }
#endif

#if 0
// Original
        template <typename R, typename = std::enable_if_t<std::is_floating_point_v<R>>>
        static void do_dot_localX(const typename mapT::iterator lstart,
                                  const typename mapT::iterator lend,
                                  typename FunctionImpl<R, NDIM>::mapT* rmap_ptr,
                                  const bool sym,
                                  Tensor<TENSOR_RESULT_TYPE(T, R)>* result_ptr,
                                  Mutex* mutex) {
            if (TensorTypeData<T>::iscomplex) MADNESS_EXCEPTION("no complex trace in LowRankTensor, sorry", 1);
            Tensor<TENSOR_RESULT_TYPE(T, R)>& result = *result_ptr;
            Tensor<TENSOR_RESULT_TYPE(T, R)> r(result.dim(0), result.dim(1));
            for (typename mapT::iterator lit = lstart; lit != lend; ++lit) {
                const keyT& key = lit->first;
                typename FunctionImpl<R, NDIM>::mapT::iterator rit = rmap_ptr->find(key);
                if (rit != rmap_ptr->end()) {
                    const mapvecT& leftv = lit->second;
                    const typename FunctionImpl<R, NDIM>::mapvecT& rightv = rit->second;
                    const int nleft = leftv.size();
                    const int nright = rightv.size();

                    for (int iv = 0; iv < nleft; iv++) {
                        const int i = leftv[iv].first;
                        const GenTensor<T>* iptr = leftv[iv].second;

                        for (int jv = 0; jv < nright; jv++) {
                            const int j = rightv[jv].first;
                            const GenTensor<R>* jptr = rightv[jv].second;

                            if (!sym || (sym && i <= j))
                                r(i, j) += iptr->trace_conj(*jptr);
                        }
                    }
                }
            }
            mutex->lock();
            result += r;
            mutex->unlock();
        }
#else
       template <typename R>
       static void do_dot_localX(const typename mapT::iterator lstart,
                                 const typename mapT::iterator lend,
                                 typename FunctionImpl<R, NDIM>::mapT* rmap_ptr,
                                 const bool sym,
                                 Tensor<TENSOR_RESULT_TYPE(T, R)>* result_ptr,
                                 Mutex* mutex) {
           Tensor<TENSOR_RESULT_TYPE(T, R)>& result = *result_ptr;
           // Tensor<TENSOR_RESULT_TYPE(T, R)> r(result.dim(0), result.dim(1));
           for (typename mapT::iterator lit = lstart; lit != lend; ++lit) {
               const keyT& key = lit->first;
               typename FunctionImpl<R, NDIM>::mapT::iterator rit = rmap_ptr->find(key);
               if (rit != rmap_ptr->end()) {
                   const mapvecT& leftv = lit->second;
                   const typename FunctionImpl<R, NDIM>::mapvecT& rightv = rit->second;
                   const size_t nleft = leftv.size();
                   const size_t nright= rightv.size();

                   unsigned int size = leftv[0].second->size();
                   Tensor<T> Left(nleft, size);
                   Tensor<R> Right(nright, size);
                   Tensor< TENSOR_RESULT_TYPE(T, R)> r(nleft, nright);
                   for(unsigned int iv = 0; iv < nleft; ++iv) Left(iv, _) = (*(leftv[iv].second)).full_tensor();
                   for(unsigned int jv = 0; jv < nright; ++jv) Right(jv, _) = (*(rightv[jv].second)).full_tensor();
                   // call mxmT from mxm.h in tensor
                   mxmT(nleft, nright, size, r.ptr(), Left.ptr(), Right.ptr());
                   mutex->lock();
                   for(unsigned int iv = 0; iv < nleft; ++iv) {
                       const int i = leftv[iv].first;
                       for(unsigned int jv = 0; jv < nright; ++jv) {
                         const int j = rightv[jv].first;
                         if (!sym || (sym && i <= j)) result(i, j) += r(iv, jv);
                       }
                   }
                   mutex->unlock();
               }
           }
       }
#endif

        static double conj(float x) {
            return x;
        }

        static std::complex<double> conj(const std::complex<double> x) {
            return std::conj(x);
        }

        template <typename R>
        static Tensor< TENSOR_RESULT_TYPE(T,R) >
        inner_local(const std::vector<const FunctionImpl<T,NDIM>*>& left,
                    const std::vector<const FunctionImpl<R,NDIM>*>& right,
                    bool sym) {

            // This is basically a sparse matrix^T * matrix product
            // Rij = sum(k) Aki * Bkj
            // where i and j index functions and k index the wavelet coeffs
            // eventually the goal is this structure (don't have jtile yet)
            //
            // do in parallel tiles of k (tensors of coeffs)
            //    do tiles of j
            //       do i
            //          do j in jtile
            //             do k in ktile
            //                Rij += Aki*Bkj

            mapT lmap = make_key_vec_map(left);
            typename FunctionImpl<R,NDIM>::mapT rmap;
            auto* rmap_ptr = (typename FunctionImpl<R,NDIM>::mapT*)(&lmap);
            if ((std::vector<const FunctionImpl<R,NDIM>*>*)(&left) != &right) {
                rmap = FunctionImpl<R,NDIM>::make_key_vec_map(right);
                rmap_ptr = &rmap;
            }

            size_t chunk = (lmap.size()-1)/(3*4*5)+1;

            Tensor< TENSOR_RESULT_TYPE(T,R) > r(left.size(), right.size());
            Mutex mutex;

            typename mapT::iterator lstart=lmap.begin();
            while (lstart != lmap.end()) {
                typename mapT::iterator lend = lstart;
                advance(lend,chunk);
                left[0]->world.taskq.add(&FunctionImpl<T,NDIM>::do_inner_localX<R>, lstart, lend, rmap_ptr, sym, &r, &mutex);
                lstart = lend;
            }
            left[0]->world.taskq.fence();

            if (sym) {
                for (long i=0; i<r.dim(0); i++) {
                    for (long j=0; j<i; j++) {
                        TENSOR_RESULT_TYPE(T,R) sum = r(i,j)+conj(r(j,i));
                        r(i,j) = sum;
                        r(j,i) = conj(sum);
                    }
                }
            }
            return r;
        }

        template <typename R>
        static Tensor<TENSOR_RESULT_TYPE(T, R)>
        dot_local(const std::vector<const FunctionImpl<T, NDIM>*>& left,
                  const std::vector<const FunctionImpl<R, NDIM>*>& right,
                  bool sym) {

            // This is basically a sparse matrix * matrix product
            // Rij = sum(k) Aik * Bkj
            // where i and j index functions and k index the wavelet coeffs
            // eventually the goal is this structure (don't have jtile yet)
            //
            // do in parallel tiles of k (tensors of coeffs)
            //    do tiles of j
            //       do i
            //          do j in jtile
            //             do k in ktile
            //                Rij += Aik*Bkj

            mapT lmap = make_key_vec_map(left);
            typename FunctionImpl<R, NDIM>::mapT rmap;
            auto* rmap_ptr = (typename FunctionImpl<R, NDIM>::mapT*)(&lmap);
            if ((std::vector<const FunctionImpl<R, NDIM>*>*)(&left) != &right) {
                rmap = FunctionImpl<R, NDIM>::make_key_vec_map(right);
                rmap_ptr = &rmap;
            }

            size_t chunk = (lmap.size() - 1) / (3 * 4 * 5) + 1;

            Tensor<TENSOR_RESULT_TYPE(T, R)> r(left.size(), right.size());
            Mutex mutex;

            typename mapT::iterator lstart=lmap.begin();
            while (lstart != lmap.end()) {
                typename mapT::iterator lend = lstart;
                advance(lend, chunk);
                left[0]->world.taskq.add(&FunctionImpl<T, NDIM>::do_dot_localX<R>, lstart, lend, rmap_ptr, sym, &r, &mutex);
                lstart = lend;
            }
            left[0]->world.taskq.fence();

            // sym is for hermiticity
            if (sym) {
                for (long i = 0; i < r.dim(0); i++) {
                    for (long j = 0; j < i; j++) {
                        TENSOR_RESULT_TYPE(T, R) sum = r(i, j) + conj(r(j, i));
                        r(i, j) = sum;
                        r(j, i) = conj(sum);
                    }
                }
            }
            return r;
        }

        template <typename R>
        void print_type_in_compilation_error(R&&)
        {
            static_assert(!std::is_same<R, int>::value &&
                          std::is_same<R, int>::value,
                          "Compilation failed because you wanted to know the type; see below:");
        }

        /// invoked by result

        /// contract 2 functions f(x,z) = \int g(x,y) * h(y,z) dy
        /// @tparam CDIM: the dimension of the contraction variable (y)
        /// @tparam NDIM: the dimension of the result (x,z)
        /// @tparam LDIM: the dimension of g(x,y)
        /// @tparam KDIM: the dimension of h(y,z)
        template<typename Q, std::size_t LDIM, typename R, std::size_t KDIM,
                std::size_t CDIM = (KDIM + LDIM - NDIM) / 2>
        void partial_inner(const FunctionImpl<Q, LDIM>& g, const FunctionImpl<R, KDIM>& h,
                           const std::array<int, CDIM> v1, const std::array<int, CDIM> v2) {

            typedef std::multimap<Key<NDIM>, std::list<Key<CDIM>>> contractionmapT;
            //double wall_get_lists=0.0;
            //double wall_recur=0.0;
            //double wall_contract=0.0;
            std::size_t nmax=FunctionDefaults<CDIM>::get_max_refine_level();
            const double thresh=FunctionDefaults<NDIM>::get_thresh();

            // auto print_map = [](const auto& map) {
                // for (const auto& kv : map) print(kv.first,"--",kv.second);
            // };
            // logical constness, not bitwise constness
            FunctionImpl<Q,LDIM>& g_nc=const_cast<FunctionImpl<Q,LDIM>&>(g);
            FunctionImpl<R,KDIM>& h_nc=const_cast<FunctionImpl<R,KDIM>&>(h);

            std::list<contractionmapT> all_contraction_maps;
            for (std::size_t n=0; n<nmax; ++n) {

                // list of nodes with d coefficients (and their parents)
	      //double wall0 = wall_time();
                auto [g_ijlist, g_jlist] = g.get_contraction_node_lists(n, v1);
                auto [h_ijlist, h_jlist] = h.get_contraction_node_lists(n, v2);
                if ((g_ijlist.size() == 0) and (h_ijlist.size() == 0)) break;
                //double wall1 = wall_time();
                //wall_get_lists += (wall1 - wall0);
                //wall0 = wall1;
//                print("g_jlist");
//                for (const auto& kv : g_jlist) print(kv.first,kv.second);
//                print("h_jlist");
//                for (const auto& kv : h_jlist) print(kv.first,kv.second);

                // next lines will insert s nodes into g and h -> possible race condition!
                bool this_first = true;  // are the remaining indices of g before those of g: f(x,z) = g(x,y) h(y,z)
                // CDIM, NDIM, KDIM
                contractionmapT contraction_map = g_nc.recur_down_for_contraction_map(
                        g_nc.key0(), g_nc.get_coeffs().find(g_nc.key0()).get()->second, v1, v2,
                        h_ijlist, h_jlist, this_first, thresh);

                this_first = false;
                // CDIM, NDIM, LDIM
                auto hnode0=h_nc.get_coeffs().find(h_nc.key0()).get()->second;
                contractionmapT contraction_map1 = h_nc.recur_down_for_contraction_map(
                        h_nc.key0(), hnode0, v2, v1,
                        g_ijlist, g_jlist, this_first, thresh);

                // will contain duplicate entries
                contraction_map.merge(contraction_map1);
                // turn multimap into a map of list
                auto it = contraction_map.begin();
                while (it != contraction_map.end()) {
                    auto it_end = contraction_map.upper_bound(it->first);
                    auto it2 = it;
                    it2++;
                    while (it2 != it_end) {
                        it->second.splice(it->second.end(), it2->second);
                        it2 = contraction_map.erase(it2);
                    }
                    it = it_end;
                }
//                print("thresh               ",thresh);
//                print("contraction list size",contraction_map.size());

                // remove all double entries
                for (auto& elem: contraction_map) {
                    elem.second.sort();
                    elem.second.unique();
                }
                //wall1 = wall_time();
                //wall_recur += (wall1 - wall0);
//                if (n==2) {
//                    print("contraction map for n=", n);
//                    print_map(contraction_map);
//                }
                all_contraction_maps.push_back(contraction_map);

                long mapsize=contraction_map.size();
                if (mapsize==0) break;
            }


            // finally do the contraction
            for (const auto& contraction_map : all_contraction_maps) {
                for (const auto& key_list : contraction_map) {
                    const Key<NDIM>& key=key_list.first;
                    const std::list<Key<CDIM>>& list=key_list.second;
                    woT::task(coeffs.owner(key), &implT:: template partial_inner_contract<Q,LDIM,R,KDIM>,
                              &g,&h,v1,v2,key,list);
                }
            }
        }

        /// for contraction two functions f(x,z) = \int g(x,y) h(y,z) dy

        /// find all nodes with d coefficients and return a list of complete keys and of
        /// keys holding only the y dimension, also the maximum norm of all d for the j dimension
        /// @param[in]  n   the scale
        /// @param[in]  v   array holding the indices of the integration variable
        /// @return     ijlist: list of all nodes with d coeffs; jlist: j-part of ij list only
        template<std::size_t CDIM>
        std::tuple<std::set<Key<NDIM>>, std::map<Key<CDIM>,double>>
        get_contraction_node_lists(const std::size_t n, const std::array<int, CDIM>& v) const {

            const auto& cdata=get_cdata();
            auto has_d_coeffs = [&cdata](const coeffT& coeff) {
                if (coeff.has_no_data()) return false;
                return (coeff.dim(0)==2*cdata.k);
            };

            // keys to be contracted in g
            std::set<Key<NDIM>> ij_list;           // full key
            std::map<Key<CDIM>,double> j_list;      // only that dimension that will be contracted

            for (auto it=get_coeffs().begin(); it!=get_coeffs().end(); ++it) {
                const Key<NDIM>& key=it->first;
                const FunctionNode<T,NDIM>& node=it->second;
                if ((key.level()==n) and (has_d_coeffs(node.coeff()))) {
                    ij_list.insert(key);
                    Vector<Translation,CDIM> j_trans;
                    for (std::size_t i=0; i<CDIM; ++i) j_trans[i]=key.translation()[v[i]];
                    Key<CDIM> jkey(n,j_trans);
                    const double max_d_norm=j_list[jkey];
                    j_list.insert_or_assign(jkey,std::max(max_d_norm,node.get_dnorm()));
                    Key<CDIM> parent_jkey=jkey.parent();
                    while (j_list.count(parent_jkey)==0) {
                        j_list.insert({parent_jkey,1.0});
                        parent_jkey=parent_jkey.parent();
                    }
                }
            }
            return std::make_tuple(ij_list,j_list);
        }

        /// make a map of all nodes that will contribute to a partial inner product

        /// given the list of d coefficient-holding nodes of the other function:
        /// recur down h if snorm * dnorm > tol and key n−jx ∈ other−ij-list. Make s
        /// coefficients if necessary. Make list of nodes n − ijk as map(n-ik, list(j)).
        ///
        /// !! WILL ADD NEW S NODES TO THIS TREE THAT MUST BE REMOVED TO AVOID INCONSISTENT TREE STRUCTURE !!
        ///
        /// @param[in]  key     for recursion
        /// @param[in]  node    corresponds to key
        /// @param[in]  v_this  this' dimension that are contracted
        /// @param[in]  v_other other's dimension that are contracted
        /// @param[in]  ij_other_list  list of nodes of the other function that will be contracted (and their parents)
        /// @param[in]  j_other_list  list of column nodes of the other function that will be contracted (and their parents)
        /// @param[in]  max_d_norm  max d coeff norm of the nodes in j_list
        /// @param[in]  this_first  are the remaining coeffs of this functions first or last in the result function
        /// @param[in]  thresh  threshold for including nodes in the contraction: snorm*dnorm > thresh
        /// @tparam     CDIM    dimension to be contracted
        /// @tparam     ODIM    dimensions of the other function
        /// @tparam     FDIM    dimensions of the final function
        template<std::size_t CDIM, std::size_t ODIM, std::size_t FDIM=NDIM+ODIM-2*CDIM>
        std::multimap<Key<FDIM>, std::list<Key<CDIM>>> recur_down_for_contraction_map(
                const keyT& key, const nodeT& node,
                const std::array<int,CDIM>& v_this,
                const std::array<int,CDIM>& v_other,
                const std::set<Key<ODIM>>& ij_other_list,
                const std::map<Key<CDIM>,double>& j_other_list,
                bool this_first, const double thresh) {

            std::multimap<Key<FDIM>, std::list<Key<CDIM>>> contraction_map;

            // fast return if the other function has no d coeffs
            if (j_other_list.empty()) return contraction_map;

            // continue recursion if this node may be contracted with the j column
            // extract relevant node translations from this node
            const auto j_this_key=key.extract_key(v_this);

//            print("\nkey, j_this_key", key, j_this_key);
            const double max_d_norm=j_other_list.find(j_this_key)->second;
            const bool sd_norm_product_large = node.get_snorm() * max_d_norm > truncate_tol(thresh,key);
//            print("sd_product_norm",node.get_snorm() * max_d_norm, thresh);

            // end recursion if we have reached the final scale n
            // with which nodes from other will this node be contracted?
            bool final_scale=key.level()==ij_other_list.begin()->level();
            if (final_scale and sd_norm_product_large) {
                for (auto& other_key : ij_other_list) {
                    const auto j_other_key=other_key.extract_key(v_other);
                    if (j_this_key != j_other_key) continue;
                    auto i_key=key.extract_complement_key(v_this);
                    auto k_key=other_key.extract_complement_key(v_other);
//                    print("key, ij_other_key",key,other_key);
//                    print("i, k, j key",i_key, k_key, j_this_key);
                    Key<FDIM> ik_key=(this_first) ? i_key.merge_with(k_key) : k_key.merge_with(i_key);
//                    print("ik_key",ik_key);
//                    MADNESS_CHECK(contraction_map.count(ik_key)==0);
                    contraction_map.insert(std::make_pair(ik_key,std::list<Key<CDIM>>{j_this_key}));
                }
                return contraction_map;
            }

            bool continue_recursion = (j_other_list.count(j_this_key)==1);
            if (not continue_recursion) return contraction_map;


            // continue recursion if norms are large
            continue_recursion = (node.has_children() or sd_norm_product_large);

            if (continue_recursion) {
                // in case we need to compute children's coefficients: unfilter only once
                bool compute_child_s_coeffs=true;
                coeffT d = node.coeff();
//                print("continuing recursion from key",key);

                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    keyT child=kit.key();
                    typename dcT::accessor acc;

                    // make child's s coeffs if it doesn't exist or if is has no s coeffs
                    bool childnode_exists=get_coeffs().find(acc,child);
                    bool need_s_coeffs= childnode_exists ? (acc->second.get_snorm()<=0.0) : true;

                    coeffT child_s_coeffs;
                    if (need_s_coeffs and compute_child_s_coeffs) {
                        if (d.dim(0)==cdata.vk[0]) {        // s coeffs only in this node
                            coeffT d1(cdata.v2k,get_tensor_args());
                            d1(cdata.s0)+=d;
                            d=d1;
                        }
                        d = unfilter(d);
                        child_s_coeffs=copy(d(child_patch(child)));
                        child_s_coeffs.reduce_rank(thresh);
                        compute_child_s_coeffs=false;
                    }

                    if (not childnode_exists) {
                        get_coeffs().replace(child,nodeT(child_s_coeffs,false));
                        get_coeffs().find(acc,child);
                    } else if (childnode_exists and need_s_coeffs) {
                        acc->second.coeff()=child_s_coeffs;
                    }
                    bool exists= get_coeffs().find(acc,child);
                    MADNESS_CHECK(exists);
                    nodeT& childnode = acc->second;
                    if (need_s_coeffs) childnode.recompute_snorm_and_dnorm(get_cdata());
//                    print("recurring down to",child);
                    contraction_map.merge(recur_down_for_contraction_map(child,childnode, v_this, v_other,
                                                                         ij_other_list, j_other_list, this_first, thresh));
//                    print("contraction_map.size()",contraction_map.size());
                }

            }

            return contraction_map;
        }


        /// tensor contraction part of partial_inner

        /// @param[in]  g   rhs of the inner product
        /// @param[in]  h   lhs of the inner product
        /// @param[in]  v1  dimensions of g to be contracted
        /// @param[in]  v2  dimensions of h to be contracted
        /// @param[in]  key key of result's (this) FunctionNode
        /// @param[in]  j_key_list list of contraction index-j keys contributing to this' node
        template<typename Q, std::size_t LDIM, typename R, std::size_t KDIM,
                std::size_t CDIM = (KDIM + LDIM - NDIM) / 2>
        void partial_inner_contract(const FunctionImpl<Q, LDIM>* g, const FunctionImpl<R, KDIM>* h,
                               const std::array<int, CDIM> v1, const std::array<int, CDIM> v2,
                               const Key<NDIM>& key, const std::list<Key<CDIM>>& j_key_list) {

            Key<LDIM - CDIM> i_key;
            Key<KDIM - CDIM> k_key;
            key.break_apart(i_key, k_key);

            coeffT result_coeff(get_cdata().v2k, get_tensor_type());
            for (const auto& j_key: j_key_list) {

                auto v_complement = [](const auto& v, const auto& vc) {
                    constexpr std::size_t VDIM = std::tuple_size<std::decay_t<decltype(v)>>::value;
                    constexpr std::size_t VCDIM = std::tuple_size<std::decay_t<decltype(vc)>>::value;
                    std::array<int, VCDIM> result;
                    for (std::size_t i = 0; i < VCDIM; i++) result[i] = (v.back() + i + 1) % (VDIM + VCDIM);
                    return result;
                };
                auto make_ij_key = [&v_complement](const auto i_key, const auto j_key, const auto& v) {
                    constexpr std::size_t IDIM = std::decay_t<decltype(i_key)>::static_size;
                    constexpr std::size_t JDIM = std::decay_t<decltype(j_key)>::static_size;
                    static_assert(JDIM == std::tuple_size<std::decay_t<decltype(v)>>::value);

                    Vector<Translation, IDIM + JDIM> l;
                    for (std::size_t i = 0; i < v.size(); ++i) l[v[i]] = j_key.translation()[i];
                    std::array<int, IDIM> vc1;
                    auto vc = v_complement(v, vc1);
                    for (std::size_t i = 0; i < vc.size(); ++i) l[vc[i]] = i_key.translation()[i];

                    return Key<IDIM + JDIM>(i_key.level(), l);
                };

                Key<LDIM> ij_key = make_ij_key(i_key, j_key, v1);
                Key<KDIM> jk_key = make_ij_key(k_key, j_key, v2);

                MADNESS_CHECK(g->get_coeffs().probe(ij_key));
                MADNESS_CHECK(h->get_coeffs().probe(jk_key));
                const coeffT& gcoeff = g->get_coeffs().find(ij_key).get()->second.coeff();
                const coeffT& hcoeff = h->get_coeffs().find(jk_key).get()->second.coeff();
                coeffT gcoeff1, hcoeff1;
                if (gcoeff.dim(0) == g->get_cdata().k) {
                    gcoeff1 = coeffT(g->get_cdata().v2k, g->get_tensor_args());
                    gcoeff1(g->get_cdata().s0) += gcoeff;
                } else {
                    gcoeff1 = gcoeff;
                }
                if (hcoeff.dim(0) == g->get_cdata().k) {
                    hcoeff1 = coeffT(h->get_cdata().v2k, h->get_tensor_args());
                    hcoeff1(h->get_cdata().s0) += hcoeff;
                } else {
                    hcoeff1 = hcoeff;
                }

                // offset: 0 for full tensor, 1 for svd representation with rand being the first dimension (r,d1,d2,d3) -> (r,d1*d2*d3)
                auto fuse = [](Tensor<T> tensor, const std::array<int, CDIM>& v, int offset) {
                    for (std::size_t i = 0; i < CDIM - 1; ++i) {
                        MADNESS_CHECK((v[i] + 1) == v[i + 1]); // make sure v is contiguous and ascending
                        tensor = tensor.fusedim(v[0]+offset);
                    }
                    return tensor;
                };

                // use case: partial_projection of 2-electron functions in svd representation  f(1) = \int g(2) h(1,2) d2
                // c_i = \sum_j a_j b_ij = \sum_jr a_j b_rj b'_rj
                //                       = \sum_jr ( a_j b_rj) b'_rj )
                auto contract2 = [](const auto& svdcoeff, const auto& tensor, const int particle) {
#if HAVE_GENTENSOR
                    const int spectator_particle=(particle+1)%2;
                    Tensor<Q> gtensor = svdcoeff.get_svdtensor().make_vector_with_weights(particle);
                    gtensor=gtensor.reshape(svdcoeff.rank(),gtensor.size()/svdcoeff.rank());
                    MADNESS_CHECK(gtensor.ndim()==2);
                    Tensor<Q> gtensor_other = svdcoeff.get_svdtensor().ref_vector(spectator_particle);
                    Tensor<T> tmp1=inner(gtensor,tensor.flat(),1,0);          // tmp1(r) = sum_j a'_(r,j) b(j)
                    MADNESS_CHECK(tmp1.ndim()==1);
                    Tensor<T> tmp2=inner(gtensor_other,tmp1,0,0);          // tmp2(i) = sum_r a_(r,i) tmp1(r)
                    return tmp2;
#else
                    MADNESS_EXCEPTION("no partial_inner using svd without GenTensor",1);
                    return Tensor<T>();
#endif
                };

                if (gcoeff.is_full_tensor() and hcoeff.is_full_tensor() and result_coeff.is_full_tensor()) {
                    // merge multiple contraction dimensions into one
                    int offset = 0;
                    Tensor<Q> gtensor = fuse(gcoeff1.full_tensor(), v1, offset);
                    Tensor<R> htensor = fuse(hcoeff1.full_tensor(), v2, offset);
                    result_coeff.full_tensor() += inner(gtensor, htensor, v1[0], v2[0]);
                    if (key.level() > 0) {
                        gtensor = copy(gcoeff1.full_tensor()(g->get_cdata().s0));
                        htensor = copy(hcoeff1.full_tensor()(h->get_cdata().s0));
                        gtensor = fuse(gtensor, v1, offset);
                        htensor = fuse(htensor, v2, offset);
                        result_coeff.full_tensor()(get_cdata().s0) -= inner(gtensor, htensor, v1[0], v2[0]);
                    }
                }


                // use case: 2-electron functions in svd representation  f(1,3) = \int g(1,2) h(2,3) d2
                // c_ik = \sum_j a_ij b_jk = \sum_jrr' a_ri a'_rj b_r'j b_r'k
                //                         = \sum_jrr' ( a_ri (a'_rj b_r'j) )  b_r'k
                //                         = \sum_jrr' c_r'i  b_r'k
                else if (gcoeff.is_svd_tensor() and hcoeff.is_svd_tensor() and result_coeff.is_svd_tensor()) {
                    MADNESS_CHECK(v1[0]==0 or v1[CDIM-1]==LDIM-1);
                    MADNESS_CHECK(v2[0]==0 or v2[CDIM-1]==KDIM-1);
                    int gparticle=  v1[0]==0 ? 0 : 1;       // which particle to integrate over
                    int hparticle=  v2[0]==0 ? 0 : 1;       // which particle to integrate over
                    // merge multiple contraction dimensions into one
                    Tensor<Q> gtensor = gcoeff1.get_svdtensor().flat_vector_with_weights(gparticle);
                    Tensor<Q> gtensor_other = gcoeff1.get_svdtensor().flat_vector((gparticle+1)%2);
                    Tensor<R> htensor = hcoeff1.get_svdtensor().flat_vector_with_weights(hparticle);
                    Tensor<R> htensor_other = hcoeff1.get_svdtensor().flat_vector((hparticle+1)%2);
                    Tensor<T> tmp1=inner(gtensor,htensor,1,1);      // tmp1(r,r') = sum_j b(r,j) a(r',j)
                    Tensor<T> tmp2=inner(tmp1,gtensor_other,0,0);   // tmp2(r',i) = sum_r tmp1(r,r') a(r,i)
                    Tensor< typename Tensor<T>::scalar_type > w(tmp2.dim(0));
                    MADNESS_CHECK(tmp2.dim(0)==htensor_other.dim(0));
                    w=1.0;
                    coeffT result_tmp(get_cdata().v2k, get_tensor_type());
                    result_tmp.get_svdtensor().set_vectors_and_weights(w,tmp2,htensor_other);
                    if (key.level() > 0) {
                        GenTensor<Q> gcoeff2 = copy(gcoeff1(g->get_cdata().s0));
                        GenTensor<R> hcoeff2 = copy(hcoeff1(h->get_cdata().s0));
                        Tensor<Q> gtensor = gcoeff2.get_svdtensor().flat_vector_with_weights(gparticle);
                        Tensor<Q> gtensor_other = gcoeff2.get_svdtensor().flat_vector((gparticle+1)%2);
                        Tensor<R> htensor = hcoeff2.get_svdtensor().flat_vector_with_weights(hparticle);
                        Tensor<R> htensor_other = hcoeff2.get_svdtensor().flat_vector((hparticle+1)%2);
                        Tensor<T> tmp1=inner(gtensor,htensor,1,1);      // tmp1(r,r') = sum_j b(r,j) a(r',j)
                        Tensor<T> tmp2=inner(tmp1,gtensor_other,0,0);   // tmp2(r',i) = sum_r tmp1(r,r') a(r,i)
                        Tensor< typename Tensor<T>::scalar_type > w(tmp2.dim(0));
                        MADNESS_CHECK(tmp2.dim(0)==htensor_other.dim(0));
                        w=1.0;
                        coeffT result_coeff1(get_cdata().vk, get_tensor_type());
                        result_coeff1.get_svdtensor().set_vectors_and_weights(w,tmp2,htensor_other);
                        result_tmp(get_cdata().s0)-=result_coeff1;
                    }
                    result_coeff+=result_tmp;
                }

                // use case: partial_projection of 2-electron functions in svd representation  f(1) = \int g(2) h(1,2) d2
                // c_i = \sum_j a_j b_ij = \sum_jr a_j b_rj b'_rj
                //                       = \sum_jr ( a_j b_rj) b'_rj )
                else if (gcoeff.is_full_tensor() and hcoeff.is_svd_tensor() and result_coeff.is_full_tensor()) {
                    MADNESS_CHECK(v1[0]==0 and v1[CDIM-1]==LDIM-1);
                    MADNESS_CHECK(v2[0]==0 or v2[CDIM-1]==KDIM-1);
                    MADNESS_CHECK(LDIM==CDIM);
                    int hparticle=  v2[0]==0 ? 0 : 1;       // which particle to integrate over

                    Tensor<T> r=contract2(hcoeff1,gcoeff1.full_tensor(),hparticle);
                    if (key.level()>0) r(get_cdata().s0)-=contract2(copy(hcoeff1(h->get_cdata().s0)),copy(gcoeff.full_tensor()(g->get_cdata().s0)),hparticle);
                    result_coeff.full_tensor()+=r;
                }
                // use case: partial_projection of 2-electron functions in svd representation  f(1) = \int g(1,2) h(2) d2
                // c_i = \sum_j a_ij b_j = \sum_jr a_ri a'_rj b_j
                //                       = \sum_jr ( a_ri (a'_rj b_j) )
                else if (gcoeff.is_svd_tensor() and hcoeff.is_full_tensor() and result_coeff.is_full_tensor()) {
                    MADNESS_CHECK(v1[0]==0 or v1[CDIM-1]==LDIM-1);
                    MADNESS_CHECK(v2[0]==0 and v2[CDIM-1]==KDIM-1);
                    MADNESS_CHECK(KDIM==CDIM);
                    int gparticle=  v1[0]==0 ? 0 : 1;       // which particle to integrate over

                    Tensor<T> r=contract2(gcoeff1,hcoeff1.full_tensor(),gparticle);
                    if (key.level()>0) r(get_cdata().s0)-=contract2(copy(gcoeff1(g->get_cdata().s0)),copy(hcoeff.full_tensor()(h->get_cdata().s0)),gparticle);
                    result_coeff.full_tensor()+=r;

                } else {
                    MADNESS_EXCEPTION("unknown case in partial_inner_contract",1);
                }
            }

            MADNESS_CHECK(result_coeff.is_assigned());
            result_coeff.reduce_rank(get_thresh());

            if (coeffs.is_local(key))
                coeffs.send(key, &nodeT::accumulate, result_coeff, coeffs, key, get_tensor_args());
            else
                coeffs.task(key, &nodeT::accumulate, result_coeff, coeffs, key, get_tensor_args(), TaskAttributes::hipri());
        }

        /// Return the inner product with an external function on a specified function node.

        /// @param[in] key Key of the function node to compute the inner product on. (the domain of integration)
        /// @param[in] c Tensor of coefficients for the function at the function node given by key
        /// @param[in] f Reference to FunctionFunctorInterface. This is the externally provided function
        /// @return Returns the inner product over the domain of a single function node, no guarantee of accuracy.
        T inner_ext_node(keyT key, tensorT c, const std::shared_ptr< FunctionFunctorInterface<T,NDIM> > f) const {
            tensorT fvals = tensorT(this->cdata.vk);
            // Compute the value of the external function at the quadrature points.
            fcube(key, *(f), cdata.quad_x, fvals);
            // Convert quadrature point values to scaling coefficients.
            tensorT fc = tensorT(values2coeffs(key, fvals));
            // Return the inner product of the two functions' scaling coefficients.
            return c.trace_conj(fc);
        }

        /// Call inner_ext_node recursively until convergence.
        /// @param[in] key Key of the function node on which to compute inner product (the domain of integration)
        /// @param[in] c coeffs for the function at the node given by key
        /// @param[in] f Reference to FunctionFunctorInterface. This is the externally provided function
        /// @param[in] leaf_refine boolean switch to turn on/off refinement past leaf nodes
        /// @param[in] old_inner the inner product on the parent function node
        /// @return Returns the inner product over the domain of a single function, checks for convergence.
        T inner_ext_recursive(keyT key, tensorT c, const std::shared_ptr< FunctionFunctorInterface<T,NDIM> > f, const bool leaf_refine, T old_inner=T(0)) const {
            int i = 0;
            tensorT c_child, inner_child;
            T new_inner, result = 0.0;

            c_child = tensorT(cdata.v2k); // tensor of child coeffs
            inner_child = Tensor<double>(pow(2, NDIM)); // child inner products

            // If old_inner is default value, assume this is the first call
            // and compute inner product on this node.
            if (old_inner == T(0)) {
                old_inner = inner_ext_node(key, c, f);
            }

            if (coeffs.find(key).get()->second.has_children()) {
                // Since the key has children and we know the func is redundant,
                // Iterate over all children of this compute node, computing
                // the inner product on each child node. new_inner will store
                // the sum of these, yielding a more accurate inner product.
                for (KeyChildIterator<NDIM> it(key); it; ++it, ++i) {
                    const keyT& child = it.key();
                    tensorT cc = coeffs.find(child).get()->second.coeff().full_tensor_copy();
                    inner_child(i) = inner_ext_node(child, cc, f);
                }
                new_inner = inner_child.sum();
            } else if (leaf_refine) {
                // We need the scaling coefficients of the numerical function
                // at each of the children nodes. We can't use project because
                // there is no guarantee that the numerical function will have
                // a functor.  Instead, since we know we are at or below the
                // leaf nodes, the wavelet coefficients are zero (to within the
                // truncate tolerance). Thus, we can use unfilter() to
                // get the scaling coefficients at the next level.
                tensorT d = tensorT(cdata.v2k);
                d = T(0);
                d(cdata.s0) = copy(c);
                c_child = unfilter(d);

                // Iterate over all children of this compute node, computing
                // the inner product on each child node. new_inner will store
                // the sum of these, yielding a more accurate inner product.
                for (KeyChildIterator<NDIM> it(key); it; ++it, ++i) {
                    const keyT& child = it.key();
                    tensorT cc = tensorT(c_child(child_patch(child)));
                    inner_child(i) = inner_ext_node(child, cc, f);
                }
                new_inner = inner_child.sum();
            } else {
                // If we get to here, we are at the leaf nodes and the user has
                // specified that they do not want refinement past leaf nodes.
                new_inner = old_inner;
            }

            // Check for convergence. If converged...yay, we're done. If not,
            // call inner_ext_node_recursive on each child node and accumulate
            // the inner product in result.
            // if (std::abs(new_inner - old_inner) <= truncate_tol(thresh, key)) {
            if (std::abs(new_inner - old_inner) <= thresh) {
                result = new_inner;
            } else {
                i = 0;
                for (KeyChildIterator<NDIM> it(key); it; ++it, ++i) {
                    const keyT& child = it.key();
                    tensorT cc = tensorT(c_child(child_patch(child)));
                    result += inner_ext_recursive(child, cc, f, leaf_refine, inner_child(i));
                }
            }

            return result;
        }

        struct do_inner_ext_local_ffi {
            const std::shared_ptr< FunctionFunctorInterface<T, NDIM> > fref;
            const implT * impl;
            const bool leaf_refine;
            const bool do_leaves;   ///< start with leaf nodes instead of initial_level

            do_inner_ext_local_ffi(const std::shared_ptr< FunctionFunctorInterface<T,NDIM> > f,
                    const implT * impl, const bool leaf_refine, const bool do_leaves)
                    : fref(f), impl(impl), leaf_refine(leaf_refine), do_leaves(do_leaves) {};

            T operator()(typename dcT::const_iterator& it) const {
                if (do_leaves and it->second.is_leaf()) {
                    tensorT cc = it->second.coeff().full_tensor();
                    return impl->inner_adaptive_recursive(it->first, cc, fref, leaf_refine, T(0));
                } else if ((not do_leaves) and (it->first.level() == impl->initial_level)) {
                    tensorT cc = it->second.coeff().full_tensor();
                    return impl->inner_ext_recursive(it->first, cc, fref, leaf_refine, T(0));
                } else {
                    return 0.0;
                }
            }

            T operator()(T a, T b) const {
                return (a + b);
            }

            template <typename Archive> void serialize(const Archive& ar) {
                MADNESS_EXCEPTION("NOT IMPLEMENTED", 1);
            }
        };

        /// Return the local part of inner product with external function ... no communication.
        /// @param[in] f Reference to FunctionFunctorInterface. This is the externally provided function
        /// @param[in] leaf_refine boolean switch to turn on/off refinement past leaf nodes
        /// @return Returns local part of the inner product, i.e. over the domain of all function nodes on this compute node.
        T inner_ext_local(const std::shared_ptr< FunctionFunctorInterface<T,NDIM> > f, const bool leaf_refine) const {
            typedef Range<typename dcT::const_iterator> rangeT;

            return world.taskq.reduce<T, rangeT, do_inner_ext_local_ffi>(rangeT(coeffs.begin(),coeffs.end()),
                    do_inner_ext_local_ffi(f, this, leaf_refine, false));
        }

        /// Return the local part of inner product with external function ... no communication.
        /// @param[in] f Reference to FunctionFunctorInterface. This is the externally provided function
        /// @param[in] leaf_refine boolean switch to turn on/off refinement past leaf nodes
        /// @return Returns local part of the inner product, i.e. over the domain of all function nodes on this compute node.
        T inner_adaptive_local(const std::shared_ptr< FunctionFunctorInterface<T,NDIM> > f, const bool leaf_refine) const {
            typedef Range<typename dcT::const_iterator> rangeT;

            return world.taskq.reduce<T, rangeT, do_inner_ext_local_ffi>(rangeT(coeffs.begin(),coeffs.end()),
                    do_inner_ext_local_ffi(f, this, leaf_refine, true));
        }

        /// Call inner_ext_node recursively until convergence.
        /// @param[in] key Key of the function node on which to compute inner product (the domain of integration)
        /// @param[in] c coeffs for the function at the node given by key
        /// @param[in] f Reference to FunctionFunctorInterface. This is the externally provided function
        /// @param[in] leaf_refine boolean switch to turn on/off refinement past leaf nodes
        /// @param[in] old_inner the inner product on the parent function node
        /// @return Returns the inner product over the domain of a single function, checks for convergence.
        T inner_adaptive_recursive(keyT key, const tensorT& c,
                const std::shared_ptr< FunctionFunctorInterface<T,NDIM> > f,
                const bool leaf_refine, T old_inner=T(0)) const {

            // the inner product in the current node
            old_inner = inner_ext_node(key, c, f);
            T result=0.0;

            // the inner product in the child nodes

            // compute the sum coefficients of the MRA function
            tensorT d = tensorT(cdata.v2k);
            d = T(0);
            d(cdata.s0) = copy(c);
            tensorT c_child = unfilter(d);

            // compute the inner product in the child nodes
            T new_inner=0.0; // child inner products
            for (KeyChildIterator<NDIM> it(key); it; ++it) {
                const keyT& child = it.key();
                tensorT cc = tensorT(c_child(child_patch(child)));
                new_inner+= inner_ext_node(child, cc, f);
            }

            // continue recursion if needed
            const double tol=truncate_tol(thresh,key);
            if (leaf_refine and (std::abs(new_inner - old_inner) > tol)) {
                for (KeyChildIterator<NDIM> it(key); it; ++it) {
                    const keyT& child = it.key();
                    tensorT cc = tensorT(c_child(child_patch(child)));
                    result += inner_adaptive_recursive(child, cc, f, leaf_refine, T(0));
                }
            } else {
                result = new_inner;
            }
            return result;

        }


        /// Return the gaxpy product with an external function on a specified
        /// function node.
        /// @param[in] key Key of the function node on which to compute gaxpy
        /// @param[in] lc Tensor of coefficients for the function at the
        ///            function node given by key
        /// @param[in] f Pointer to function of type T that takes coordT
        ///            arguments. This is the externally provided function and
        ///            the right argument of gaxpy.
        /// @param[in] alpha prefactor of c Tensor for gaxpy
        /// @param[in] beta prefactor of fcoeffs for gaxpy
        /// @return Returns coefficient tensor of the gaxpy product at specified
        ///         key, no guarantee of accuracy.
        template <typename L>
        tensorT gaxpy_ext_node(keyT key, Tensor<L> lc, T (*f)(const coordT&), T alpha, T beta) const {
            // Compute the value of external function at the quadrature points.
            tensorT fvals = madness::fcube(key, f, cdata.quad_x);
            // Convert quadrature point values to scaling coefficients.
            tensorT fcoeffs = values2coeffs(key, fvals);
            // Return the inner product of the two functions' scaling coeffs.
            tensorT c2 = copy(lc);
            c2.gaxpy(alpha, fcoeffs, beta);
            return c2;
        }

        /// Return out of place gaxpy using recursive descent.
        /// @param[in] key Key of the function node on which to compute gaxpy
        /// @param[in] left FunctionImpl, left argument of gaxpy
        /// @param[in] lcin coefficients of left at this node
        /// @param[in] c coefficients of gaxpy product at this node
        /// @param[in] f pointer to function of type T that takes coordT
        ///            arguments. This is the externally provided function and
        ///            the right argument of gaxpy.
        /// @param[in] alpha prefactor of left argument for gaxpy
        /// @param[in] beta prefactor of right argument for gaxpy
        /// @param[in] tol convergence tolerance...when the norm of the gaxpy's
        ///            difference coefficients is less than tol, we are done.
        template <typename L>
        void gaxpy_ext_recursive(const keyT& key, const FunctionImpl<L,NDIM>* left,
                                 Tensor<L> lcin, tensorT c, T (*f)(const coordT&),
                                 T alpha, T beta, double tol, bool below_leaf) {
            typedef typename FunctionImpl<L,NDIM>::dcT::const_iterator literT;

            // If we haven't yet reached the leaf level, check whether the
            // current key is a leaf node of left. If so, set below_leaf to true
            // and continue. If not, make this a parent, recur down, return.
            if (not below_leaf) {
                bool left_leaf = left->coeffs.find(key).get()->second.is_leaf();
                if (left_leaf) {
                    below_leaf = true;
                } else {
                    this->coeffs.replace(key, nodeT(coeffT(), true));
                    for (KeyChildIterator<NDIM> it(key); it; ++it) {
                        const keyT& child = it.key();
                        woT::task(left->coeffs.owner(child), &implT:: template gaxpy_ext_recursive<L>,
                                  child, left, Tensor<L>(), tensorT(), f, alpha, beta, tol, below_leaf);
                    }
                    return;
                }
            }

            // Compute left's coefficients if not provided
            Tensor<L> lc = lcin;
            if (lc.size() == 0) {
                literT it = left->coeffs.find(key).get();
                MADNESS_ASSERT(it != left->coeffs.end());
                if (it->second.has_coeff())
                    lc = it->second.coeff().reconstruct_tensor();
            }

            // Compute this node's coefficients if not provided in function call
            if (c.size() == 0) {
                c = gaxpy_ext_node(key, lc, f, alpha, beta);
            }

            // We need the scaling coefficients of the numerical function at
            // each of the children nodes. We can't use project because there
            // is no guarantee that the numerical function will have a functor.
            // Instead, since we know we are at or below the leaf nodes, the
            // wavelet coefficients are zero (to within the truncate tolerance).
            // Thus, we can use unfilter() to get the scaling coefficients at
            // the next level.
            Tensor<L> lc_child = Tensor<L>(cdata.v2k); // left's child coeffs
            Tensor<L> ld = Tensor<L>(cdata.v2k);
            ld = L(0);
            ld(cdata.s0) = copy(lc);
            lc_child = unfilter(ld);

            // Iterate over children of this node,
            // storing the gaxpy coeffs in c_child
            tensorT c_child = tensorT(cdata.v2k); // tensor of child coeffs
            for (KeyChildIterator<NDIM> it(key); it; ++it) {
                const keyT& child = it.key();
                tensorT lcoeff = tensorT(lc_child(child_patch(child)));
                c_child(child_patch(child)) = gaxpy_ext_node(child, lcoeff, f, alpha, beta);
            }

            // Compute the difference coefficients to test for convergence.
            tensorT d = tensorT(cdata.v2k);
            d = filter(c_child);
            // Filter returns both s and d coefficients, so set scaling
            // coefficient part of d to 0 so that we take only the
            // norm of the difference coefficients.
            d(cdata.s0) = T(0);
            double dnorm = d.normf();

            // Small d.normf means we've reached a good level of resolution
            // Store the coefficients and return.
            if (dnorm <= truncate_tol(tol,key)) {
                this->coeffs.replace(key, nodeT(coeffT(c,targs), false));
            } else {
                // Otherwise, make this a parent node and recur down
                this->coeffs.replace(key, nodeT(coeffT(), true)); // Interior node

                for (KeyChildIterator<NDIM> it(key); it; ++it) {
                    const keyT& child = it.key();
                    tensorT child_coeff = tensorT(c_child(child_patch(child)));
                    tensorT left_coeff = tensorT(lc_child(child_patch(child)));
                    woT::task(left->coeffs.owner(child), &implT:: template gaxpy_ext_recursive<L>,
                              child, left, left_coeff, child_coeff, f, alpha, beta, tol, below_leaf);
                }
            }
        }

        template <typename L>
        void gaxpy_ext(const FunctionImpl<L,NDIM>* left, T (*f)(const coordT&), T alpha, T beta, double tol, bool fence) {
            if (world.rank() == coeffs.owner(cdata.key0))
                gaxpy_ext_recursive<L> (cdata.key0, left, Tensor<L>(), tensorT(), f, alpha, beta, tol, false);
            if (fence)
                world.gop.fence();
        }

        /// project the low-dim function g on the hi-dim function f: result(x) = <this(x,y) | g(y)>

        /// invoked by the hi-dim function, a function of NDIM+LDIM

        /// Upon return, result matches this, with contributions on all scales
        /// @param[in]  result  lo-dim function of NDIM-LDIM \todo Should this be param[out]?
        /// @param[in]  gimpl  	lo-dim function of LDIM
        /// @param[in]  dim over which dimensions to be integrated: 0..LDIM or LDIM..LDIM+NDIM-1
        template<size_t LDIM>
        void project_out(FunctionImpl<T,NDIM-LDIM>* result, const FunctionImpl<T,LDIM>* gimpl,
                         const int dim, const bool fence) {

            const keyT& key0=cdata.key0;

            if (world.rank() == coeffs.owner(key0)) {

                // coeff_op will accumulate the result
                typedef project_out_op<LDIM> coeff_opT;
                coeff_opT coeff_op(this,result,CoeffTracker<T,LDIM>(gimpl),dim);

                // don't do anything on this -- coeff_op will accumulate into result
                typedef noop<T,NDIM> apply_opT;
                apply_opT apply_op;

                woT::task(world.rank(), &implT:: template forward_traverse<coeff_opT,apply_opT>,
                          coeff_op, apply_op, cdata.key0);

            }
            if (fence) world.gop.fence();

        }


        /// project the low-dim function g on the hi-dim function f: result(x) = <f(x,y) | g(y)>
        template<size_t LDIM>
        struct project_out_op {
            bool randomize() const {return false;}

            typedef project_out_op<LDIM> this_type;
            typedef CoeffTracker<T,LDIM> ctL;
            typedef FunctionImpl<T,NDIM-LDIM> implL1;
            typedef std::pair<bool,coeffT> argT;

            const implT* fimpl;		///< the hi dim function f
            mutable implL1* result;	///< the low dim result function
            ctL iag;				///< the low dim function g
            int dim;				///< 0: project 0..LDIM-1, 1: project LDIM..NDIM-1

            // ctor
            project_out_op() = default;
            project_out_op(const implT* fimpl, implL1* result, const ctL& iag, const int dim)
                : fimpl(fimpl), result(result), iag(iag), dim(dim) {}
            project_out_op(const project_out_op& other)
                : fimpl(other.fimpl), result(other.result), iag(other.iag), dim(other.dim) {}


            /// do the actual contraction
            Future<argT> operator()(const Key<NDIM>& key) const {

            	Key<LDIM> key1,key2,dest;
            	key.break_apart(key1,key2);

            	// make the right coefficients
                coeffT gcoeff;
                if (dim==0) {
                    gcoeff=iag.get_impl()->parent_to_child(iag.coeff(),iag.key(),key1);
                    dest=key2;
                }
                if (dim==1) {
                    gcoeff=iag.get_impl()->parent_to_child(iag.coeff(),iag.key(),key2);
                    dest=key1;
                }

                MADNESS_ASSERT(fimpl->get_coeffs().probe(key));		// must be local!
                const nodeT& fnode=fimpl->get_coeffs().find(key).get()->second;
                const coeffT& fcoeff=fnode.coeff();

                // fast return if possible
                if (fcoeff.has_no_data() or gcoeff.has_no_data())
                    return Future<argT> (argT(fnode.is_leaf(),coeffT()));;

                MADNESS_CHECK(gcoeff.is_full_tensor());
                tensorT final(result->cdata.vk);
                const int k=fcoeff.dim(0);
                const int k_ldim=std::pow(k,LDIM);
                std::vector<long> shape(LDIM, k);

                if (fcoeff.is_full_tensor()) {
                    //  result_i = \sum_j g_j f_ji
                    const tensorT gtensor = gcoeff.full_tensor().reshape(k_ldim);
                    const tensorT ftensor = fcoeff.full_tensor().reshape(k_ldim,k_ldim);
                    final=inner(gtensor,ftensor,0,dim).reshape(shape);

                } else if (fcoeff.is_svd_tensor()) {
                    if (fcoeff.rank()>0) {

                        //  result_i = \sum_jr g_j a_rj w_r b_ri
                        const int otherdim = (dim + 1) % 2;
                        const tensorT gtensor = gcoeff.full_tensor().flat();
                        const tensorT atensor = fcoeff.get_svdtensor().flat_vector(dim);  // a_rj
                        const tensorT btensor = fcoeff.get_svdtensor().flat_vector(otherdim);
                        const tensorT gatensor = inner(gtensor, atensor, 0, 1);              // ga_r
                        tensorT weights = copy(fcoeff.get_svdtensor().weights_);
                        weights.emul(gatensor); // ga_r * w_r
                        // sum over all ranks of b, include new weights:
                        // result_i = \sum_r ga_r * w_r * b_ri
                        for (int r = 0; r < fcoeff.rank(); ++r) final += weights(r) * btensor(r, _);
                        final = final.reshape(shape);
                    }

                } else {
                    MADNESS_EXCEPTION("unsupported tensor type in project_out_op",1);
                }

                // accumulate the result
                result->coeffs.task(dest, &FunctionNode<T,LDIM>::accumulate2, final, result->coeffs, dest, TaskAttributes::hipri());

                return Future<argT> (argT(fnode.is_leaf(),coeffT()));
            }

            this_type make_child(const keyT& child) const {
            	Key<LDIM> key1,key2;
            	child.break_apart(key1,key2);
            	const Key<LDIM> gkey = (dim==0) ? key1 : key2;

                return this_type(fimpl,result,iag.make_child(gkey),dim);
            }

            /// retrieve the coefficients (parent coeffs might be remote)
            Future<this_type> activate() const {
            	Future<ctL> g1=iag.activate();
                return result->world.taskq.add(detail::wrap_mem_fn(*const_cast<this_type *> (this),
                                               &this_type::forward_ctor),fimpl,result,g1,dim);
            }

            /// taskq-compatible ctor
            this_type forward_ctor(const implT* fimpl1, implL1* result1, const ctL& iag1, const int dim1) {
            	return this_type(fimpl1,result1,iag1,dim1);
            }

            template <typename Archive> void serialize(const Archive& ar) {
                ar & result & iag & fimpl & dim;
            }

        };


        /// project the low-dim function g on the hi-dim function f: this(x) = <f(x,y) | g(y)>

        /// invoked by result, a function of NDIM

        /// @param[in]  f   hi-dim function of LDIM+NDIM
        /// @param[in]  g   lo-dim function of LDIM
        /// @param[in]  dim over which dimensions to be integrated: 0..LDIM or LDIM..LDIM+NDIM-1
        template<size_t LDIM>
        void project_out2(const FunctionImpl<T,LDIM+NDIM>* f, const FunctionImpl<T,LDIM>* g, const int dim) {

            typedef std::pair< keyT,coeffT > pairT;
            typedef typename FunctionImpl<T,NDIM+LDIM>::dcT::const_iterator fiterator;

            // loop over all nodes of hi-dim f, compute the inner products with all
            // appropriate nodes of g, and accumulate in result
            fiterator end = f->get_coeffs().end();
            for (fiterator it=f->get_coeffs().begin(); it!=end; ++it) {
                const Key<LDIM+NDIM> key=it->first;
                const FunctionNode<T,LDIM+NDIM> fnode=it->second;
                const coeffT& fcoeff=fnode.coeff();

                if (fnode.is_leaf() and fcoeff.has_data()) {

                    // break key into particle: over key1 will be summed, over key2 will be
                    // accumulated, or vice versa, depending on dim
                    if (dim==0) {
                        Key<NDIM> key1;
                        Key<LDIM> key2;
                        key.break_apart(key1,key2);

                        Future<pairT> result;
                        //                        sock_it_to_me(key1, result.remote_ref(world));
                        g->task(coeffs.owner(key1), &implT::sock_it_to_me, key1, result.remote_ref(world), TaskAttributes::hipri());
                        woT::task(world.rank(),&implT:: template do_project_out<LDIM>,fcoeff,result,key1,key2,dim);

                    } else if (dim==1) {
                        Key<LDIM> key1;
                        Key<NDIM> key2;
                        key.break_apart(key1,key2);

                        Future<pairT> result;
                        //                        sock_it_to_me(key2, result.remote_ref(world));
                        g->task(coeffs.owner(key2), &implT::sock_it_to_me, key2, result.remote_ref(world), TaskAttributes::hipri());
                        woT::task(world.rank(),&implT:: template do_project_out<LDIM>,fcoeff,result,key2,key1,dim);

                    } else {
                        MADNESS_EXCEPTION("confused dim in project_out",1);
                    }
                }
            }
            set_tree_state(redundant);
//            this->compressed=false;
//            this->nonstandard=false;
//            this->redundant=true;
        }


        /// compute the inner product of two nodes of only some dimensions and accumulate on result

        /// invoked by result
        /// @param[in]  fcoeff  coefficients of high dimension LDIM+NDIM
        /// @param[in]  gpair   key and coeffs of low dimension LDIM (possibly a parent node)
        /// @param[in]  gkey    key of actual low dim node (possibly the same as gpair.first, iff gnode exists)
        /// @param[in]  dest    destination node for the result
        /// @param[in]  dim     which dimensions should be contracted: 0..LDIM-1 or LDIM..NDIM+LDIM-1
        template<size_t LDIM>
        void do_project_out(const coeffT& fcoeff, const std::pair<keyT,coeffT> gpair, const keyT& gkey,
                            const Key<NDIM>& dest, const int dim) const {

            const coeffT gcoeff=parent_to_child(gpair.second,gpair.first,gkey);

            // fast return if possible
            if (fcoeff.has_no_data() or gcoeff.has_no_data()) return;

            // let's specialize for the time being on SVD tensors for f and full tensors of half dim for g
            MADNESS_ASSERT(gcoeff.tensor_type()==TT_FULL);
            MADNESS_ASSERT(fcoeff.tensor_type()==TT_2D);
            const tensorT gtensor=gcoeff.full_tensor();
            tensorT result(cdata.vk);

            const int otherdim=(dim+1)%2;
            const int k=fcoeff.dim(0);
            std::vector<Slice> s(fcoeff.config().dim_per_vector()+1,_);

            // do the actual contraction
            for (int r=0; r<fcoeff.rank(); ++r) {
                s[0]=Slice(r,r);
                const tensorT contracted_tensor=fcoeff.config().ref_vector(dim)(s).reshape(k,k,k);
                const tensorT other_tensor=fcoeff.config().ref_vector(otherdim)(s).reshape(k,k,k);
                const double ovlp= gtensor.trace_conj(contracted_tensor);
                const double fac=ovlp * fcoeff.config().weights(r);
                result+=fac*other_tensor;
            }

            // accumulate the result
            coeffs.task(dest, &nodeT::accumulate2, result, coeffs, dest, TaskAttributes::hipri());
        }




        /// Returns the maximum local depth of the tree ... no communications.
        std::size_t max_local_depth() const;


        /// Returns the maximum depth of the tree ... collective ... global sum/broadcast
        std::size_t max_depth() const;

        /// Returns the max number of nodes on a processor
        std::size_t max_nodes() const;

        /// Returns the min number of nodes on a processor
        std::size_t min_nodes() const;

        /// Returns the size of the tree structure of the function ... collective global sum
        std::size_t tree_size() const;

        /// Returns the number of coefficients in the function for each rank
        std::size_t size_local() const;

        /// Returns the number of coefficients in the function ... collective global sum
        std::size_t size() const;

        /// Returns the number of coefficients in the function ... collective global sum
        std::size_t nCoeff() const;

        /// Returns the number of coefficients in the function ... collective global sum
        std::size_t real_size() const;

        /// print tree size and size
        void print_size(const std::string name) const;

        /// print the number of configurations per node
        void print_stats() const;

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

        /// Hash a pointer to \c FunctionImpl

        /// \param[in] impl pointer to a FunctionImpl
        /// \return The hash.
        inline friend hashT hash_value(const FunctionImpl<T,NDIM>* pimpl) {
            hashT seed = hash_value(pimpl->id().get_world_id());
            detail::combine_hash(seed, hash_value(pimpl->id().get_obj_id()));
            return seed;
        }

        /// Hash a shared_ptr to \c FunctionImpl

        /// \param[in] impl pointer to a FunctionImpl
        /// \return The hash.
        inline friend hashT hash_value(const std::shared_ptr<FunctionImpl<T,NDIM>> impl) {
            return hash_value(impl.get());
        }
    };

    namespace archive {
        template <class Archive, class T, std::size_t NDIM>
        struct ArchiveLoadImpl<Archive,const FunctionImpl<T,NDIM>*> {
            static void load(const Archive& ar, const FunctionImpl<T,NDIM>*& ptr) {
                bool exists=false;
                ar & exists;
                if (exists) {
                    uniqueidT id;
                    ar & id;
                    World* world = World::world_from_id(id.get_world_id());
                    MADNESS_ASSERT(world);
                    auto ptr_opt = world->ptr_from_id< WorldObject< FunctionImpl<T,NDIM> > >(id);
                    if (!ptr_opt)
                      MADNESS_EXCEPTION("FunctionImpl: remote operation attempting to use a locally uninitialized object",0);
                    ptr = static_cast< const FunctionImpl<T,NDIM>*>(*ptr_opt);
                    if (!ptr)
                        MADNESS_EXCEPTION("FunctionImpl: remote operation attempting to use an unregistered object",0);
                } else {
                    ptr=nullptr;
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
                bool exists=false;
                ar & exists;
                if (exists) {
                    uniqueidT id;
                    ar & id;
                    World* world = World::world_from_id(id.get_world_id());
                    MADNESS_ASSERT(world);
                    auto ptr_opt = world->ptr_from_id< WorldObject< FunctionImpl<T,NDIM> > >(id);
                    if (!ptr_opt)
                      MADNESS_EXCEPTION("FunctionImpl: remote operation attempting to use a locally uninitialized object",0);
                    ptr = static_cast< FunctionImpl<T,NDIM>*>(*ptr_opt);
                    if (!ptr)
                      MADNESS_EXCEPTION("FunctionImpl: remote operation attempting to use an unregistered object",0);
                } else {
                    ptr=nullptr;
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
                const FunctionImpl<T,NDIM>* f = nullptr;
                ArchiveLoadImpl<Archive, const FunctionImpl<T,NDIM>*>::load(ar, f);
                ptr.reset(f, [] (const FunctionImpl<T,NDIM> *p_) -> void {});
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
                FunctionImpl<T,NDIM>* f = nullptr;
                ArchiveLoadImpl<Archive, FunctionImpl<T,NDIM>*>::load(ar, f);
                ptr.reset(f, [] (FunctionImpl<T,NDIM> *p_) -> void {});
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
