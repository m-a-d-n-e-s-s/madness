/*
  This file is part of MADNESS.

  Copyright (C) <2007> <Oak Ridge National Laboratory>

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


#ifndef MAD_FUNC_DATA
#define MAD_FUNC_DATA

/// \file funcimpl.h
/// \brief Provides FunctionCommonData, FunctionImpl and FunctionFactory

#include <world/world.h>
#include <misc/misc.h>
#include <tensor/mtrand.h>
#include <tensor/tensor.h>
#include <mra/key.h>
#include <mra/funcdefaults.h>

namespace madness {
    template <typename T, int NDIM> class FunctionImpl;
    template <typename T, int NDIM> class Function;
    template <int D> class LoadBalImpl;
//    template <int D> class LBTreeImpl;
    template <int D> class LBTree;
    template <int D> class MyPmap;
}


namespace madness {


    /// A simple process map soon to be supplanted by Rebecca's
    template <typename keyT>
    class SimpleMap : public WorldDCPmapInterface< keyT > {
    private:
        const int nproc;
        const ProcessID me;
        const int n;

    public:
        SimpleMap(World& world, int n = 4) : nproc(world.nproc()), me(world.rank()), n(n) {}

        ProcessID owner(const keyT& key) const {
            if (key.level() == 0) {
                return 0;
            }
            else if (key.level() <= n) {
                return hash(key)%nproc;
            }
            else {
                return hash(key.parent(key.level()-n))%nproc;
            }
        }
    };


    /// FunctionCommonData holds all Function data common for given k

    /// Since Function assignment and copy constructors are shallow it
    /// greatly simplifies maintaining consistent state to have all
    /// (permanent) state encapsulated in a single class.  The state
    /// is shared between instances using a SharedPtr.  Also,
    /// separating shared from instance specific state accelerates the
    /// constructor, which is important for massive parallelism, and
    /// permitting inexpensive use of temporaries.  The default copy
    /// constructor and assignment operator are used but are probably
    /// never invoked.
    template <typename T, int NDIM>
    class FunctionCommonData {
    private:
        static FunctionCommonData<T,NDIM> data[MAXK+1]; /// Declared in mra.cc, initialized on first use

        /// Private.  Make the level-0 blocks of the periodic central difference derivative operator
        void _make_dc_periodic();

        /// Private.  Initialize the twoscale coefficients
        void _init_twoscale();

        /// Private.  Do first use initialization
        void _initialize(int k) {
            this->k = k;
            npt = k;
            for (int i=0; i<4; i++)
                s[i] = Slice(i*k,(i+1)*k-1);
            s0 = std::vector<Slice>(NDIM);
            sh = std::vector<Slice>(NDIM);
            vk = std::vector<long>(NDIM);
            vq = std::vector<long>(NDIM);
            v2k = std::vector<long>(NDIM);
            for (int i=0; i<NDIM; i++) {
                s0[i] = s[0];
                sh[i] = Slice(0,(k-1)/2);
                vk[i] = k;
                vq[i] = npt;
                v2k[i] = 2*k;
            }
            work1 = tensorT(vk,false);
            work2 = tensorT(v2k,false);
            workq = tensorT(vq,false);
            zero_coeff = tensorT(vk);
            key0 = Key<NDIM>(0,Vector<Translation,NDIM>(0));

            _init_twoscale();
            _init_quadrature(k, npt, quad_x, quad_w, quad_phi, quad_phiw, quad_phit);
            _make_dc_periodic();
            initialized = true;
        }

        FunctionCommonData()
            : initialized(false)
        {}

        bool initialized;
    public:
        typedef Tensor<T> tensorT; ///< Type of tensor used to hold coeff

        int k;                  ///< order of the wavelet
        int npt;                ///< no. of quadrature points
        Slice s[4];             ///< s[0]=Slice(0,k-1), s[1]=Slice(k,2*k-1), etc.
        std::vector<Slice> s0;  ///< s[0] in each dimension to get scaling coeff
        std::vector<Slice> sh;  ///< Slice(0,(k-1)/2) in each dimension for autorefine test
        std::vector<long> vk;   ///< (k,...) used to initialize Tensors
        std::vector<long> v2k;  ///< (2k,...) used to initialize Tensors
        std::vector<long> vq;   ///< (npt,...) used to initialize Tensors

        mutable Tensor<T> work1;///< work space of size (k,...)
        mutable Tensor<T> work2;///< work space of size (2k,...)
        mutable Tensor<T> workq;///< work space of size (npt,...)
        tensorT zero_coeff;     ///< Zero (k,...) tensor for internal convenience of diff

        Key<NDIM> key0;         ///< Key for root node

        Tensor<double> quad_x;  ///< quadrature points
        Tensor<double> quad_w;  ///< quadrature weights
        Tensor<double> quad_phi; ///< quad_phi(i,j) = at x[i] value of phi[j]
        Tensor<double> quad_phit; ///< transpose of quad_phi
        Tensor<double> quad_phiw; ///< quad_phiw(i,j) = at x[i] value of w[i]*phi[j]

        Tensor<double> h0, h1, g0, g1; ///< The separate blocks of twoscale coefficients
        Tensor<double> hg, hgT; ///< The full twoscale coeff (2k,2k) and transpose
        Tensor<double> hgsonly; ///< hg[0:k,:]

        Tensor<double> rm, r0, rp;        ///< Blocks of the derivative operator
        Tensor<double> rm_left, rm_right, rp_left, rp_right; ///< Rank-1 forms rm & rp

        static const FunctionCommonData<T,NDIM>& get(int k) {
            MADNESS_ASSERT(k>0 && k<=MAXK);
            if (!data[k].initialized) data[k]._initialize(k);
            return data[k];
        }

        /// Initialize the quadrature information

        /// Made public with all arguments thru interface for reuse in FunctionImpl::err_box
        static void _init_quadrature(int k, int npt,
                                     Tensor<double>& quad_x, Tensor<double>& quad_w,
                                     Tensor<double>& quad_phi, Tensor<double>& quad_phiw,
                                     Tensor<double>& quad_phit);
    };

    /// Interface required for functors used as input to Functions
    template <typename T, int NDIM>
    class FunctionFunctorInterface {
    public:
        virtual T operator()(const Vector<double,NDIM>& x) const = 0;
	virtual ~FunctionFunctorInterface() {}
    };

    /// FunctionFactory implements the named-parameter idiom for Function

    /// C++ does not provide named arguments (as does, e.g., Python).
    /// This class provides something very close.  Create functions as follows
    /// \code
    /// double myfunc(const double x[]);
    /// Function<double,3> f = FunctionFactory<double,3>(world).f(myfunc).k(11).thresh(1e-9).debug()
    /// \endcode
    /// where the methods of function factory, which specify the non-default
    /// arguments eventually passed to the \c Function constructor, can be
    /// used in any order.
    ///
    /// Need to add a general functor for initial projection with a standard interface.
    template <typename T, int NDIM>
    class FunctionFactory {
        friend class FunctionImpl<T,NDIM>;
        typedef Vector<double,NDIM> coordT;            ///< Type of vector holding coordinates
    protected:
        World& _world;
        int _k;
        double _thresh;
        int _initial_level;
        int _max_refine_level;
        int _truncate_mode;
        bool _refine;
        bool _empty;
        bool _autorefine;
        bool _truncate_on_project;
        bool _fence;
        Tensor<int> _bc;
        SharedPtr< WorldDCPmapInterface< Key<NDIM> > > _pmap;
        SharedPtr< FunctionFunctorInterface<T,NDIM> > _functor;

        struct FunctorInterfaceWrapper : public FunctionFunctorInterface<T,NDIM> {
            T (*f)(const coordT&);

            FunctorInterfaceWrapper(T (*f)(const coordT&)) : f(f) {};

            T operator()(const coordT& x) const {return f(x);};
        };

    public:
        FunctionFactory(World& world)
            : _world(world)
            , _k(FunctionDefaults<NDIM>::get_k())
            , _thresh(FunctionDefaults<NDIM>::get_thresh())
            , _initial_level(FunctionDefaults<NDIM>::get_initial_level())
            , _max_refine_level(FunctionDefaults<NDIM>::get_max_refine_level())
            , _truncate_mode(FunctionDefaults<NDIM>::get_truncate_mode())
            , _refine(FunctionDefaults<NDIM>::get_refine())
            , _empty(false)
            , _autorefine(FunctionDefaults<NDIM>::get_autorefine())
            , _truncate_on_project(FunctionDefaults<NDIM>::get_truncate_on_project())
            , _fence(true)
            , _bc(FunctionDefaults<NDIM>::get_bc())
            , _pmap(FunctionDefaults<NDIM>::get_pmap())
            , _functor(0)
        {}
        FunctionFactory& functor(const SharedPtr< FunctionFunctorInterface<T,NDIM> >& functor) {
            _functor = functor;
            return *this;
        }
        FunctionFactory& f(T (*f)(const coordT&)) {
            functor(SharedPtr< FunctionFunctorInterface<T,NDIM> >(new FunctorInterfaceWrapper(f)));
            return *this;
        }
        FunctionFactory& k(int k) {
            _k = k;
            return *this;
        }
        FunctionFactory& thresh(double thresh) {
            _thresh = thresh;
            return *this;
        }
        FunctionFactory& initial_level(int initial_level) {
            _initial_level = initial_level;
            return *this;
        }
        FunctionFactory& max_refine_level(int max_refine_level) {
            _max_refine_level = max_refine_level;
            return *this;
        }
        FunctionFactory& truncate_mode(int truncate_mode) {
            _truncate_mode = truncate_mode;
            return *this;
        }
        FunctionFactory& refine(bool refine = true) {
            _refine = refine;
            return *this;
        }
        FunctionFactory& norefine(bool norefine = true) {
            _refine = !norefine;
            return *this;
        }
        FunctionFactory& bc(const Tensor<int>& bc) {
            _bc = copy(bc);
            return *this;
        }
        FunctionFactory& empty() {
            _empty = true;
            return *this;
        }
        FunctionFactory& autorefine() {
            _autorefine = true;
            return *this;
        }
        FunctionFactory& noautorefine() {
            _autorefine = false;
            return *this;
        }
        FunctionFactory& truncate_on_project() {
            _truncate_on_project = true;
            return *this;
        }
        FunctionFactory& notruncate_on_project() {
            _truncate_on_project = false;
            return *this;
        }
        FunctionFactory& fence(bool fence=true) {
            _fence = fence;
            return *this;
        }
        FunctionFactory& nofence() {
            _fence = false;
            return *this;
        }
        FunctionFactory& pmap(const SharedPtr< WorldDCPmapInterface< Key<NDIM> > >& pmap) {
            _pmap = pmap;
            return *this;
        }
    };

    /// FunctionNode holds the coefficients, etc., at each node of the 2^NDIM-tree
    template <typename T, int NDIM>
    class FunctionNode {
    private:
        Tensor<T> _coeffs;  ///< The coefficients, if any
        double _norm_tree;  ///< After norm_tree will contain norm of coefficients summed up tree
        bool _has_children; ///< True if there are children

    public:
        typedef WorldContainer< Key<NDIM>, FunctionNode<T,NDIM> > dcT;  ///< Type of container holding the nodes
        /// Default constructor makes node without coeff or children
        FunctionNode()
            : _coeffs()
            , _norm_tree(1e300)
            , _has_children(false)
        {}

        /// Constructor from given coefficients with optional children

        /// Note that only a shallow copy of the coeff are taken so
        /// you should pass in a deep copy if you want the node to
        /// take ownership.
        explicit FunctionNode(const Tensor<T>& coeff, bool has_children=false)
            : _coeffs(coeff)
            , _norm_tree(1e300)
            , _has_children(has_children)
        {}

        /// Copy with possible type conversion of coefficients, copying all other state

        /// Choose to not overload copy and type conversion operators
        /// so there are no automatic type conversions.
        template <typename Q>
        FunctionNode<Q,NDIM> convert() const {
            return FunctionNode<Q,NDIM>(copy(_coeffs),_has_children);
        }

        /// Returns true if there are coefficients in this node
        bool has_coeff() const {return (_coeffs.size>0);}

        /// Returns true if this node has children
        bool has_children() const {return _has_children;}

        /// Returns true if this does not have children
        bool is_leaf() const {return !_has_children;}

        /// Returns true if this node is invalid (no coeffs and no children)
        bool is_invalid() const {return !(has_coeff() || has_children());}

        /// Returns a non-const reference to the tensor containing the coeffs

        /// Returns an empty tensor if there are no coefficients.
        Tensor<T>& coeff() {
            MADNESS_ASSERT(_coeffs.dim[0]<=2*MAXK && _coeffs.dim[0]>=0);
            return _coeffs;}

        /// Returns a const reference to the tensor containing the coeffs

        /// Returns an empty tensor if there are no coefficeints.
        const Tensor<T>& coeff() const {return _coeffs;}

        /// Sets \c has_children attribute to value of \c flag.
        Void set_has_children(bool flag) {_has_children = flag; return None;}

        /// Sets \c has_children attribute to true recurring up to ensure connected
        Void set_has_children_recursive(const typename FunctionNode<T,NDIM>::dcT& c, const Key<NDIM>& key) {
            PROFILE_MEMBER_FUNC(FunctionNode);
            if (!(_has_children || has_coeff() || key.level()==0)) {
                // If node already knows it has children or it has
                // coefficients then it must already be connected to
                // its parent.  If not, the node was probably just
                // created for this operation and must be connected to
                // its parent.
                Key<NDIM> parent = key.parent();
                const_cast<dcT&>(c).send(parent, &FunctionNode<T,NDIM>::set_has_children_recursive, c, parent);
            }
            _has_children = true;
            return None;
        }

        /// Sets \c has_children attribute to value of \c !flag
        void set_is_leaf(bool flag) {_has_children = !flag;}

        /// Takes a \em shallow copy of the coeff --- same as \c this->coeff()=coeff
        void set_coeff(const Tensor<T>& coeff) {
            _coeffs = coeff;
	    if ((coeff.dim[0] < 0) || (coeff.dim[0]>2*MAXK)) {
	      print("set_coeff: may have a problem");
	      print("set_coeff: coeff.dim[0] =", coeff.dim[0], ", 2* MAXK =", 2*MAXK);
	    }
            MADNESS_ASSERT(coeff.dim[0]<=2*MAXK && coeff.dim[0]>=0);
        }

        /// Clears the coefficients (has_coeff() will subsequently return false)
        void clear_coeff() {_coeffs = Tensor<T>();}

        /// Sets the value of norm_tree
        Void set_norm_tree(double norm_tree) {_norm_tree = norm_tree; return None;}

        /// Gets the value of norm_tree
        double get_norm_tree() const {return _norm_tree;}

        /// General bi-linear operation --- this = this*alpha + other*beta

        /// This/other may not have coefficients.  Has_children will be
        /// true in the result if either this/other have children.
        template <typename Q, typename R>
        Void gaxpy_inplace(const T& alpha, const FunctionNode<Q,NDIM>& other, const R& beta) {
            PROFILE_MEMBER_FUNC(FuncNode);
            if (other.has_children()) _has_children = true;
            if (has_coeff()) {
                if (other.has_coeff()) {
                    coeff().gaxpy(alpha,other.coeff(),beta);
                }
                else {
                    coeff().scale(alpha);
                }
            }
            else if (other.has_coeff()) {
                _coeffs = other.coeff()*beta; //? Is this the correct type conversion?
            }
            return None;
        }


        template <typename Archive>
        void serialize(Archive& ar) {
            ar & _coeffs & _has_children & _norm_tree;
        }
    };

    template <typename T, int NDIM>
    std::ostream& operator<<(std::ostream& s, const FunctionNode<T,NDIM>& node) {
        s << "(" << node.has_coeff() << ", " << node.has_children() << ", ";
        double norm = node.has_coeff() ? node.coeff().normf() : 0.0;
        if (norm < 1e-12) norm = 0.0;
        s << norm << ")";
        return s;
    }

    /// ApplyTime is a class for finding out the time taken in the apply
    /// function.

    template <int NDIM>
    class ApplyTime {
      typedef Key<NDIM> keyT;
      typedef WorldContainer<keyT, double> dcT;
      typedef std::pair<const keyT, double> datumT;

    private:
      World& world;
      dcT hash_table;
      double decay_val;

    public:
      ApplyTime(World& world) 
	: world(world)
	, hash_table(dcT(world))
	, decay_val(0.9)
	{}

      void set(datumT data) {
	hash_table.insert(data);
      }

      void clear() {
	hash_table.clear();
      }

      double get(keyT& key) {
	typename dcT::iterator it = hash_table.find(key);
	if (it == hash_table.end()) {
	  return 0.0;
	}
	else {
	  return it->second;
	}
      }

      double get(const keyT& key) {
	typename dcT::iterator it = hash_table.find(key);
	if (it == hash_table.end()) {
	  return 0.0;
	}
	else {
	  return it->second;
	}
      }

      /* datumT get(keyT& key) { */
/* 	double result = this->get(key); */
/* 	return datumT(key, result); */
/*       } */

      void update(datumT data) {
	typename dcT::iterator it = hash_table.find(data.first).get();
	if (it == hash_table.end()) {
	  hash_table.insert(data);
	}
	else {
	  double s = it->second, y = data.second;
	  data.second = s + (y-s)*decay_val;
	  hash_table.insert(data);
	}
      }

      void update(keyT key, double d) {
	update(datumT(key, d));
      }

      void print() {
	for (typename dcT::iterator it = hash_table.begin(); it != hash_table.end(); ++it) {
	  madness::print(it->first, "  ", it->second);
	}
      }

    };



    /// FunctionImpl holds all Function state to facilitate shallow copy semantics

    /// Since Function assignment and copy constructors are shallow it
    /// greatly simplifies maintaining consistent state to have all
    /// (permanent) state encapsulated in a single class.  The state
    /// is shared between instances using a SharedPtr<FunctionImpl>.
    ///
    /// The FunctionImpl inherits all of the functionality of WorldContainer
    /// (to store the coefficients) and WorldObject<WorldContainer> (used
    /// for RMI and for its unqiue id).
    template <typename T, int NDIM>
    class FunctionImpl : public WorldObject< FunctionImpl<T,NDIM> > {
    public:
	//friend class Function<T,NDIM>;
        template <typename Q, int D> friend class Function;
        template <typename Q, int D> friend class FunctionImpl;
	friend class LoadBalImpl<NDIM>;
	friend class LBTree<NDIM>;

        typedef FunctionImpl<T,NDIM> implT;       ///< Type of this class (implementation)
        typedef Tensor<T> tensorT;                     ///< Type of tensor used to hold coeffs
        typedef Vector<Translation,NDIM> tranT;         ///< Type of array holding translation
        typedef Key<NDIM> keyT;                        ///< Type of key
        typedef FunctionNode<T,NDIM> nodeT;            ///< Type of node
        typedef WorldContainer<keyT,nodeT> dcT;  ///< Type of container holding the coefficients
        typedef std::pair<const keyT,nodeT> datumT;    ///< Type of entry in container
        typedef Vector<double,NDIM> coordT;            ///< Type of vector holding coordinates

        World& world;
    private:
        int k;                  ///< Wavelet order
        double thresh;          ///< Screening threshold
        int initial_level;      ///< Initial level for refinement
        int max_refine_level;   ///< Do not refine below this level
        int truncate_mode;      ///< 0=default=(|d|<thresh), 1=(|d|<thresh/2^n), 1=(|d|<thresh/4^n);
        bool autorefine;        ///< If true, autorefine where appropriate
        bool truncate_on_project; ///< If true projection inserts at level n-1 not n
        bool nonstandard;       ///< If true, compress keeps scaling coeff

        const FunctionCommonData<T,NDIM>& cdata;

        SharedPtr< FunctionFunctorInterface<T,NDIM> > functor;

        bool compressed;        ///< Compression status

        dcT coeffs;             ///< The coefficients

        const Tensor<int> bc;     ///< Type of boundary condition -- currently only zero or periodic

	SharedPtr<ApplyTime<NDIM> > apply_time;

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
            , cdata(FunctionCommonData<T,NDIM>::get(k))
            , functor(factory._functor)
            , compressed(false)
            , coeffs(world,factory._pmap,false)
            , bc(factory._bc)
	    , apply_time(SharedPtr<ApplyTime<NDIM> > ())
        {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            // !!! Ensure that all local state is correctly formed
            // before invoking process_pending for the coeffs and
            // for this.  Otherwise, there is a race condition.
            MADNESS_ASSERT(k>0 && k<=MAXK);

            bool empty = factory._empty;
            bool do_refine = factory._refine;

            if (do_refine) initial_level = std::max(0,initial_level - 1);

            if (empty) {        // Do not set any coefficients at all
            }
            else if (functor) { // Project function and optionally refine
                insert_zero_down_to_initial_level(cdata.key0);
                for(typename dcT::iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
                    if (it->second.is_leaf()) task(coeffs.owner(it->first),
                                                   &implT::project_refine_op,
                                                   it->first, do_refine);
                }
            }
            else {  // Set as if a zero function
                initial_level = 1;
                insert_zero_down_to_initial_level(keyT(0));
            }

            coeffs.process_pending();
            this->process_pending();

            if (factory._fence && functor) world.gop.fence();
        }

        /// Copy constructor

        /// Allocates a \em new function in preparation for a deep copy
        ///
        /// By default takes pmap from other but can also specify a different pmap.
        /// Does \em not copy the coefficients ... creates an empty container.
        template <typename Q>
        FunctionImpl(const FunctionImpl<Q,NDIM>& other,
                     const SharedPtr< WorldDCPmapInterface< Key<NDIM> > >& pmap,
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
            , cdata(FunctionCommonData<T,NDIM>::get(k))
            , functor()
            , compressed(other.compressed)
            , coeffs(world, pmap ? pmap : other.coeffs.get_pmap())
            , bc(other.bc)
	  , apply_time(other.apply_time)
        {
            if (dozero) {
                initial_level = 1;
                insert_zero_down_to_initial_level(cdata.key0);
            }
            coeffs.process_pending();
            this->process_pending();
        }

	const SharedPtr< WorldDCPmapInterface< Key<NDIM> > >& get_pmap() const {
	    return coeffs.get_pmap();
	}

        /// Copy coeffs from other into self
        template <typename Q>
	void copy_coeffs(const FunctionImpl<Q,NDIM>& other, bool fence) {
	    for(typename FunctionImpl<Q,NDIM>::dcT::const_iterator it=other.coeffs.begin();
                it!=other.coeffs.end(); ++it) {
                const keyT& key = it->first;
                const typename FunctionImpl<Q,NDIM>::nodeT& node = it->second;
		coeffs.insert(key,node. template convert<Q>());
	    }
	    if (fence) world.gop.fence();
	}


        /// Inplace general bilinear operation
        template <typename Q, typename R>
        void gaxpy_inplace(const T& alpha,const FunctionImpl<Q,NDIM>& other, const R& beta, bool fence) {
            // Loop over coefficients in other that are local and then send an AM to
            // coeffs in self ... this is so can efficiently add functions with
            // different distributions.  Use an AM rather than a task to reduce memory
            // footprint on the remote end.
	    for(typename FunctionImpl<Q,NDIM>::dcT::const_iterator it=other.coeffs.begin();
                it!=other.coeffs.end();
                ++it) {
                const keyT& key = it->first;
                const typename FunctionImpl<Q,NDIM>::nodeT& other_node = it->second;
                coeffs.send(key, &nodeT:: template gaxpy_inplace<Q,R>, alpha, other_node, beta);
	    }
	    if (fence) world.gop.fence();
	}

        template <typename Archive>
        void load(Archive& ar) {
            int kk;
            ar & kk;

            MADNESS_ASSERT(kk==k);

            // note that functor should not be (re)stored
            ar & thresh & initial_level & max_refine_level & truncate_mode
                & autorefine & truncate_on_project & nonstandard & compressed & bc;

            ar & coeffs;
        }


        template <typename Archive>
        void store(Archive& ar) {
            // note that functor should not be (re)stored
            ar & k & thresh & initial_level & max_refine_level & truncate_mode
                & autorefine & truncate_on_project & nonstandard & compressed & bc;

            ar & coeffs;
        }



        /// Returns true if the function is compressed.
        bool is_compressed() const {return compressed;}

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
            if (tol <= 0.0) tol = thresh;
            if (world.rank() == coeffs.owner(cdata.key0)) truncate_spawn(cdata.key0,tol);
            if (fence) world.gop.fence();
        }


        /// Returns true if after truncation this node has coefficients

        /// Assumed to be invoked on process owning key.  Possible non-blocking
        /// communication.
        Future<bool> truncate_spawn(const keyT& key, double tol);

        /// Actually do the truncate operation
        bool truncate_op(const keyT& key, double tol, const std::vector< Future<bool> >& v);

        /// Evaluate function at quadrature points in the specified box
        void fcube(const keyT& key, const FunctionFunctorInterface<T,NDIM>& f, const Tensor<double>& qx, tensorT& fval) const;

        const keyT& key0() const {
            return cdata.key0;
        }

        void print_tree(Level maxlevel = 10000) const;

        void do_print_tree(const keyT& key, Level maxlevel) const;

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

        /// !! Returns a reference to \em STATIC data !!
        /// Can be optimized away by precomputing.
        const std::vector<Slice>& child_patch(const keyT& child) const {
            static std::vector<Slice> s(NDIM);
            const Vector<Translation,NDIM>& l = child.translation();
            for (int i=0; i<NDIM; i++) s[i] = cdata.s[l[i]&1]; // Lowest bit of translation
            return s;
        }


        /// Projection with optional refinement
        Void project_refine_op(const keyT& key, bool do_refine);


        /// Compute the Legendre scaling functions for multiplication

        /// Evaluate parent polyn at quadrature points of a child.  The prefactor of
        /// 2^n/2 is included.  The tensor must be preallocated as phi(k,npt).
        /// Refer to the implementation notes for more info.
        void phi_for_mul(Level np, Translation lp, Level nc, Translation lc, Tensor<double>& phi) const;

        /// Directly project parent coeffs to child coeffs

        /// Currently used by diff, but other uses can be anticipated
        const tensorT parent_to_child(const tensorT& s, const keyT& parent, const keyT& child) const;


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
            	for (int d=0; d<NDIM; d++) {
                	dim[d     ] = N;
                	dim[d+NDIM] = cdata.k;
		}
		tensorT rr(2*NDIM,dim);
		r0=r=rr;
		//NNkk->MqMqkk, since fuse is not allowed. Now needs to move back to 2*NDIM, since tensor 
max dim is 6
		//for (int d=NDIM-1; d>=0; --d) r.splitdim_inplace_base(d,M,q);
            } else {
	  	long dim[2*NDIM];
                for (int d=0; d<NDIM; d++) {
                        //dim[d+NDIM*2] = M;
                        dim[d+NDIM  ] = N;
			dim[d       ] = cdata.k;
                }
		tensorT rr(2*NDIM,dim);
		r0=r=rr;
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
		r.cycledim_inplace_base(NDIM,0,-1); //->NNkk or MqMqkk
	    }
	    //print("faking done M q r(fake) r0(real)",M,q,"\n", r,r0);
            ProcessID me = world.rank();
            Vector<long,NDIM> t(N);
            for (IndexIterator it(t); it; ++it) {
                keyT key(n, Vector<Translation,NDIM>(*it));
                if (coeffs.owner(key) == me) {
                    typename dcT::iterator it = coeffs.find(key).get();
                    tensorT qq;

                    if (it == coeffs.end()) {
                        // must get from above
                        typedef std::pair< keyT,Tensor<T> > pairT;
                        Future<pairT> result;
                        sock_it_to_me(key,  result.remote_ref(world));
                        const keyT& parent = result.get().first;
                        const tensorT& t = result.get().second;

                        qq = parent_to_child(t, parent, key);
                    }
                    else {
                        qq = it->second.coeff();
                    }
                    std::vector<Slice> s(NDIM*2);
		    long ll = 0;
                    for (int d=0; d<NDIM; d++) {
                        Translation l = key.translation()[d];
			ll += (l % q)*pow((double)M,NDIM)*pow((double)q,NDIM-d-1) + 
(l/q)*pow((double)M,NDIM-d-1);
			//print(d,l,(l % q)*pow(M,NDIM)*pow(q,NDIM-d-1) + (l/q)*pow(M,NDIM-d-1));

			//print("translation",l);
			//s[d       ] = Slice(l,l,0);
                        //s[d+NDIM  ] = Slice(l%q,l%q,0);
                        //s[d+NDIM] = Slice(0,k-1,1);
                    }
		    //long dum = ll;
		    for (int d=0; d<NDIM; d++) {
		    	Translation l = ll / pow((double)N,NDIM-d-1);
			s[d     ] = Slice(l,l,0);
			s[d+NDIM] = Slice(0,k-1,1);
			ll = ll % long(pow((double)N,NDIM-d-1));
		    }
		    //print(s, dum, key.translation());
                    r(s) = qq(cdata.s0);
                }
            }

            world.gop.fence();
            world.gop.sum(r0);
	    //print(r,r0);

            return r0;
        }


        
        /// Compute the function values for multiplication

        /// Given coefficients from a parent cell, compute the value of
        /// the functions at the quadrature points of a child
        template <typename Q>
        Tensor<Q> fcube_for_mul(const keyT& child, const keyT& parent, const Tensor<Q>& coeff) const {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            if (child.level() == parent.level()) {
                double scale = pow(2.0,0.5*NDIM*parent.level())/sqrt(FunctionDefaults<NDIM>::get_cell_volume());
                return transform(coeff,cdata.quad_phit).scale(scale);
            }
            else if (child.level() < parent.level()) {
                MADNESS_EXCEPTION("FunctionImpl: fcube_for_mul: child-parent relationship bad?",0);
            }
            else {
                Tensor<double> phi[NDIM];
                for (int d=0; d<NDIM; d++) {
                    phi[d] = Tensor<double>(cdata.k,cdata.npt);
                    phi_for_mul(parent.level(),parent.translation()[d],
                                child.level(), child.translation()[d], phi[d]);
                }
                return general_transform(coeff,phi).scale(1.0/sqrt(FunctionDefaults<NDIM>::get_cell_volume()));;
            }
        }

        /// Invoked as a task by mul with the actual coefficients
        template <typename L, typename R>
        Void do_mul(const keyT& key, const Tensor<L>& left, const std::pair< keyT, Tensor<R> >& arg) {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            const keyT& rkey = arg.first;
            const Tensor<R>& rcoeff = arg.second;
            //madness::print("do_mul: r", rkey, rcoeff.size);
            Tensor<R> rcube = fcube_for_mul(key, rkey, rcoeff);
            //madness::print("do_mul: l", key, left.size);
            Tensor<L> lcube = fcube_for_mul(key,  key, left);

            Tensor<T> tcube(cdata.vk,false);
            TERNARY_OPTIMIZED_ITERATOR(T, tcube, L, lcube, R, rcube, *_p0 = *_p1 * *_p2;);
            double scale = pow(0.5,0.5*NDIM*key.level())*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
            tcube = transform(tcube,cdata.quad_phiw).scale(scale);
            coeffs.insert(key, nodeT(tcube,false));
            return None;
        }


        /// Invoked by result to perform result += alpha*left+beta*right in wavelet basis

        /// Does not assume that any of result, left, right have the same distribution.
        /// For most purposes result will start as an empty so actually are implementing
        /// out of place gaxpy.  If all functions have the same distribution there is
        /// no communication except for the optional fence.
        template <typename L, typename R>
        void gaxpy(T alpha, const FunctionImpl<L,NDIM>& left,
                   T beta,  const FunctionImpl<R,NDIM>& right, bool fence) {
            // Loop over local nodes in both functions.  Add in left and subtract right.
            // Not that efficient in terms of memory bandwidth but ensures we do
            // not miss any nodes.
	    for(typename FunctionImpl<L,NDIM>::dcT::const_iterator it=left.coeffs.begin();
                it!=left.coeffs.end();
                ++it) {
                const keyT& key = it->first;
                const typename FunctionImpl<L,NDIM>::nodeT& other_node = it->second;
                coeffs.send(key, &nodeT:: template gaxpy_inplace<T,L>, 1.0, other_node, alpha);
	    }
	    for(typename FunctionImpl<R,NDIM>::dcT::const_iterator it=right.coeffs.begin();
                it!=right.coeffs.end();
                ++it) {
                const keyT& key = it->first;
                const typename FunctionImpl<L,NDIM>::nodeT& other_node = it->second;
                coeffs.send(key, &nodeT:: template gaxpy_inplace<T,R>, 1.0, other_node, beta);
	    }
            if (fence) world.gop.fence();
        }


        template <typename testT>
        void conditional_refine_doit(const testT& test, const keyT& key) {
          nodeT& node = coeffs[key];
          if (node.has_coeff() && test(key, node.coeff())) {
            tensorT s(cdata.v2k);
            s(cdata.s0) = node.coeff();
            s = unfilter(s);
            node.clear_coeff();
            node.set_has_children(true);
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
               const keyT& child = kit.key();
               task(coeffs.owner(child),
                    &implT:: template conditional_refine_insert_doit<testT>,
                    test, child, copy(s(child_patch(child))));
            }
          }
        }

        template <typename testT>
        void conditional_refine(const testT& test, bool fence) {
          MADNESS_ASSERT(!compressed);
          for(typename dcT::iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
            const keyT& key = it->first;
            conditional_refine_doit(test, key);
          }
          if (fence) world.gop.fence();
        }

        template <typename opT>
        Void unary_op_coeff_inplace_child(const opT& op, const keyT& parent, const keyT& child, const tensorT& t) {
            coeffs.insert(child, nodeT(parent_to_child(t, parent, child), false));
            nodeT& node = coeffs[child];
            op(child, node.coeff());
            return None;
        }

        /// Unary operation applied inplace to the coefficients with optional refinement and fence
        template <typename opT>
        void unary_op_coeff_inplace(bool (implT::*refineop)(const keyT&, const tensorT&) const,
                                    const opT& op,
                                    bool fence)
        {
            for(typename dcT::iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
                const keyT& parent = it->first;
                nodeT& node = it->second;
                if (node.has_coeff()) {
                    tensorT t= node.coeff();
                    if ((this->*refineop)(parent, t)) {
                        MADNESS_ASSERT(!compressed);
                        node.clear_coeff();
                        node.set_has_children(true);
                        for (KeyChildIterator<NDIM> kit(parent); kit; ++kit) {
                            const keyT& child = kit.key();
                            // Crucial that this is a task so that nodes inserted as a result of refinement
                            // are not operated on  twice due to a race condition if in parallel
                            task(coeffs.owner(child), &implT:: template unary_op_coeff_inplace_child<opT>, op, parent, child, t);
                        }
                    }
                    else {
                        op(parent, t);
                    }
                }
            }
            if (fence) world.gop.fence();
        }


        template <typename opT>
        Void unary_op_value_inplace_child(const opT& op, const keyT& parent, const keyT& child, const tensorT& t) {
            tensorT values = fcube_for_mul(child, parent, t);

            op(child, values);

            double scale = pow(0.5,0.5*NDIM*child.level())*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
            tensorT r = transform(values,cdata.quad_phiw).scale(scale);

            if (parent == child)
                coeffs[child].set_coeff(r);
            else
                coeffs.insert(child,nodeT(r,false));

            return None;
        }

        /// Unary operation applied inplace to the values with optional refinement and fence
        template <typename opT>
        void unary_op_value_inplace(bool (implT::*refineop)(const keyT&, const tensorT&) const,
                                    const opT& op,
                                    bool fence)
        {
            for(typename dcT::iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
                const keyT& parent = it->first;
                nodeT& node = it->second;
                if (node.has_coeff()) {
                    tensorT t= node.coeff();
                    if ((this->*refineop)(parent, t)) {
                        MADNESS_ASSERT(!compressed);
                        node.clear_coeff();
                        node.set_has_children(true);
                        for (KeyChildIterator<NDIM> kit(parent); kit; ++kit) {
                            const keyT& child = kit.key();
                            // Crucial that this is a task so that nodes inserted as a result of refinement
                            // are not operated on  twice due to a race condition if in parallel
                            task(coeffs.owner(child), &implT:: template unary_op_value_inplace_child<opT>, op, parent, child, t);
                        }
                    }
                    else {
                        unary_op_value_inplace_child(op, parent, parent, t);
                    }
                }
            }
            if (fence) world.gop.fence();
        }


        /// Invoked by result to compute the pointwise product result=left*right

        /// This version requires all three functions have the same distribution.
        /// Should be straightforward to do an efficient version that does not
        /// require this but I have not thought about that yet.
        ///
        /// Possible non-blocking communication and optional fence.
        template <typename L, typename R>
        void mul(const FunctionImpl<L,NDIM>& left, const FunctionImpl<R,NDIM>& right, bool fence) {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            typedef std::pair< keyT,Tensor<R> > rpairT;
            typedef std::pair< keyT,Tensor<L> > lpairT;
            MADNESS_ASSERT(coeffs.get_pmap() == left.coeffs.get_pmap() && \
                           coeffs.get_pmap() == right.coeffs.get_pmap());
            // The three possibilities for the relative position of
            // the left and right coefficients in the tree are:
            //
            // 1.  left==right
            // 2.  left>right
            // 3.  left<right
            //
            // First loop thru local coeff in left.  Handle right at the same level or above.
	    for(typename FunctionImpl<L,NDIM>::dcT::const_iterator it=left.coeffs.begin();
                it != left.coeffs.end();
                ++it) {
                const keyT& key = it->first;
                const FunctionNode<L,NDIM>& left_node = it->second;

                if (left_node.has_coeff()) {
                    if (right.coeffs.probe(key)) {
                        const FunctionNode<R,NDIM>& right_node = right.coeffs.find(key).get()->second;
                        if (right_node.has_coeff()) {
                            task(world.rank(), &implT:: template do_mul<L,R>, key, left_node.coeff(),
                                 rpairT(key,right_node.coeff()));  // Case 1.
                        }
                    }
                    else { // If right node does not exist then it must be further up the tree
                        const keyT parent = key.parent();
                        Future<rpairT> arg;
                        right.send(coeffs.owner(parent), &FunctionImpl<R,NDIM>::sock_it_to_me,
                                   parent, arg.remote_ref(world));
                        task(world.rank(), &implT:: template do_mul<L,R>, key, left_node.coeff(), arg); // Case 2.
                    }
                }
                else if (!coeffs.probe(key)) {
                    // Interior node
                    coeffs.insert(key,nodeT(tensorT(),true));
                }

            }

            // Now loop thru local coeff in right and do case 3.
	    for(typename FunctionImpl<R,NDIM>::dcT::const_iterator it=right.coeffs.begin();
                it != right.coeffs.end();
                ++it) {
                const keyT& key = it->first;
                const FunctionNode<R,NDIM>& right_node = it->second;
                if (right_node.has_coeff()) {
                    if (!left.coeffs.probe(key)) {
                        Future<lpairT> arg;
                        const keyT& parent = key.parent();
                        left.send(coeffs.owner(parent), &FunctionImpl<L,NDIM>::sock_it_to_me,
                                  parent, arg.remote_ref(world));
                        task(world.rank(), &implT:: template do_mul<R,L>, key, right_node.coeff(), arg); // Case 3.
                    }
                }
                else if (!coeffs.probe(key)) {
                    // Interior node
                    coeffs.insert(key,nodeT(tensorT(),true));
                }

            }
            if (fence) world.gop.fence();
        }

        template <typename L, typename R>
        Void do_mul_sparse2(const keyT& key,
                            const std::pair< keyT,Tensor<L> >& larg,
                            const std::pair< keyT,Tensor<R> >& rarg,
                            const FunctionImpl<R,NDIM>* right) {
            PROFILE_MEMBER_FUNC(FunctionImpl);

            if (rarg.second.size > 0) {
                if (larg.first == key) {
                    //madness::print("L*R",key,larg.first,larg.second.size,rarg.first,rarg.second.size);
                    do_mul(key, larg.second, rarg);
                }
                else {
                    //madness::print("R*L",key,larg.first,larg.second.size,rarg.first,rarg.second.size);
                    do_mul(key, rarg.second, larg);
                }
            }
            else {
                coeffs.insert(key, nodeT(tensorT(), true));  // Insert interior node
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    typedef std::pair< keyT,Tensor<R> > rpairT;
                    Future<rpairT> rarg;
                    right->send(coeffs.owner(kit.key()), &FunctionImpl<R,NDIM>::sock_it_to_me,
                                kit.key(), rarg.remote_ref(world));


                    task(world.rank(), &implT:: template do_mul_sparse2<L,R>,
                         kit.key(),larg, rarg, right);
                }
            }
            return None;
        }

        template <typename L, typename R>
        Void do_mul_sparse(const Tensor<L>& left_coeff, const FunctionImpl<R,NDIM>* right, double tol,
                           const keyT& key, double right_norm) {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            if (left_coeff.normf()*right_norm > truncate_tol(tol,key)) {
                typedef std::pair< keyT,Tensor<R> > rpairT;
                typedef std::pair< keyT,Tensor<L> > lpairT;
                Future<rpairT> rarg;
                right->send(coeffs.owner(key), &FunctionImpl<R,NDIM>::sock_it_to_me,
                            key, rarg.remote_ref(world));
                task(world.rank(), &implT:: template do_mul_sparse2<L,R>,
                     key ,lpairT(key,left_coeff), rarg, right);
            }
            else {
                coeffs.insert(key, nodeT(tensorT(cdata.vk), false));  // Result is zero
            }
            return None;
        }

        template <typename L, typename R>
        void mul_sparse(const FunctionImpl<L,NDIM>& left, const FunctionImpl<R,NDIM>& right, double tol, bool fence) {
            // I think that this should distribution agnostic
            //MADNESS_ASSERT(coeffs.get_pmap() == left.coeffs.get_pmap() && coeffs.get_pmap() == right.coeffs.get_pmap());

            // Loop thru leaf nodes in left
	    for(typename FunctionImpl<L,NDIM>::dcT::const_iterator it=left.coeffs.begin();
                it != left.coeffs.end();
                ++it) {
                const keyT& key = it->first;
                const FunctionNode<L,NDIM>& left_node = it->second;

                if (left_node.is_leaf()) {
                    Future<double> rarg = right.send(right.coeffs.owner(key), &implT::get_norm_tree_recursive, key);
                    task(world.rank(), &implT:: template do_mul_sparse<L,R>, left_node.coeff(), &right, tol, key,
                         rarg);
                }
                else {
                    coeffs.insert(key, nodeT(tensorT(), true));  // Insert interior node
                }
            }
            if (fence) world.gop.fence();
        }


        Future<double> get_norm_tree_recursive(const keyT& key) const;



        mutable long box_leaf[10000];
        mutable long box_interior[10000];

        // horrifically non-scalable
        Void put_in_box(ProcessID from, long nl, long ni) const {
            box_leaf[from] = nl;
            box_interior[from] = ni;
            return None;
        }


        /// Prints summary of data distribution
        void print_info() const {
            MADNESS_ASSERT(world.size() < 10000);
            for (int i=0; i<world.size(); i++) box_leaf[i] = box_interior[i] == 0;
            world.gop.fence();
            long nleaf=0, ninterior=0;
            for(typename dcT::const_iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
                const nodeT& node = it->second;
                if (node.is_leaf()) nleaf++;
                else ninterior++;
            }
            send(0, &implT::put_in_box, world.rank(), nleaf, ninterior);
            world.gop.fence();
            if (world.rank() == 0) {
                for (int i=0; i<world.size(); i++) {
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
                           const RemoteReference< FutureImpl< std::pair<keyT,tensorT> > >& ref) const;

        /// Evaluate a cube/slice of points ... plotlo and plothi are already in simulation coordinates

        /// No communications
        Tensor<T> eval_plot_cube(const coordT& plotlo,
                                 const coordT& plothi,
                                 const std::vector<long>& npt) const;

        /// Evaluate the function at a point in \em simulation coordinates

        /// Only the invoking process will get the result via the
        /// remote reference to a future.  Active messages may be sent
        /// to other nodes.
        Void eval(const Vector<double,NDIM>& xin,
                  const keyT& keyin,
                  const typename Future<T>::remote_refT& ref);


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
        bool noautorefine(const keyT& key, const tensorT& t) const {return false;}

        /// Returns true if this block of coeffs needs autorefining
        bool autorefine_square_test(const keyT& key, const tensorT& t) const {
            double lo, hi;
            tnorm(t, &lo, &hi);
            double test = 2*lo*hi + hi*hi;
            //print("autoreftest",key,thresh,truncate_tol(thresh, key),lo,hi,test);
            return test > truncate_tol(thresh, key);
        }


        /// Pointwise squaring of function with optional global fence

        /// If not autorefining, local computation only if not fencing.
        /// If autorefining, may result in asynchronous communication.
        void square_inplace(bool fence);

        /// Differentiation of function with optional global fence

        /// Result of differentiating f is placed into this which will
        /// have the same process map, etc., as f
        void diff(const implT& f, int axis, bool fence);


        /// Returns key of neighbor enforcing BC

        /// Out of volume keys are mapped to enforce the BC as follows.
        ///   * Periodic BC map back into the volume and return the correct key
        ///   * Zero BC - returns invalid() to indicate out of volume
        keyT neighbor(const keyT& key, int axis, int step) const;


        /// Returns key of general neighbor enforcing BC

        /// Out of volume keys are mapped to enforce the BC as follows.
        ///   * Periodic BC map back into the volume and return the correct key
        ///   * Zero BC - returns invalid() to indicate out of volume
        keyT neighbor(const keyT& key, const keyT& disp) const;


        /// Called by diff to find key and coeffs of neighbor enforcing BC

        /// Should work for any (small) step but only tested for step=+/-1
        ///
        /// do_diff1 handles the adpative refinement.  If it needs to refine, it calls
        /// forward_do_diff1 to pass the task locally or remotely with high priority.
        /// Actual differentiation is performend by do_diff2.
        Future< std::pair<keyT,tensorT> > find_neighbor(const keyT& key, int axis, int step) const;

        /// Used by diff1 to forward calls to diff1 elsewhere

        /// We cannot send a not ready future to another process
        /// so we send an active message to schedule the remote task.
        ///
        /// Local tasks that might produce additional communication
        /// are scheduled with high priority and plain-old compute
        /// tasks get normal priority.  This is an attempt to get all
        /// of the communication and adaptive refinement happening
        /// in parallel to productive computation.
        Void forward_do_diff1(const implT* f, int axis, const keyT& key,
                      const std::pair<keyT,tensorT>& left,
                      const std::pair<keyT,tensorT>& center,
                              const std::pair<keyT,tensorT>& right);


        Void do_diff1(const implT* f, int axis, const keyT& key,
                      const std::pair<keyT,tensorT>& left,
                      const std::pair<keyT,tensorT>& center,
                      const std::pair<keyT,tensorT>& right);

        Void do_diff2(const implT* f, int axis, const keyT& key,
                      const std::pair<keyT,tensorT>& left,
                      const std::pair<keyT,tensorT>& center,
                      const std::pair<keyT,tensorT>& right);

        /// Permute the dimensions according to map
        void mapdim(const implT& f, const std::vector<long>& map, bool fence);


        T eval_cube(Level n, coordT x, const tensorT c) const;


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
            return fast_transform(s,cdata.hgT,r,cdata.work2);
            //return transform(s,cdata.hgT);
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
            return fast_transform(s,cdata.hg,r,cdata.work2);
            //return transform(s, cdata.hg);
        }

        /// Projects old function into new basis (only in reconstructed form)
        void project(const implT& old, bool fence) {
            vector<Slice> s(NDIM,Slice(0,old.cdata.k-1));
            for(typename dcT::const_iterator it=old.coeffs.begin(); it!=old.coeffs.end(); ++it) {
                const keyT& key = it->first;
                const nodeT& node = it->second;
                if (node.has_coeff()) {
                    tensorT c(cdata.vk);
                    c(s) = node.coeff();
                    coeffs.insert(key,nodeT(c,false));
                }
                else {
                    coeffs.insert(key,nodeT(tensorT(),true));
                }
            }
            if (fence) world.gop.fence();
        }

        // This needed extending to accomodate a user-defined criterion
        void refine(bool fence);

//         // Ensures the function is gradually refined before application of an operator
//         void ensure_gradual(bool fence);


        // Widens the support of the tree in preparation for integral operator
        void widen(bool fence, int ndiff);

        void reconstruct(bool fence) {
            // Must set true here so that successive calls without fence do the right thing
            nonstandard = compressed = false;
            if (world.rank() == coeffs.owner(cdata.key0)) reconstruct_op(cdata.key0,tensorT());
            if (fence) world.gop.fence();
        }

        // Invoked on node where key is local
        Void reconstruct_op(const keyT& key, const tensorT& s);

        void compress(bool nonstandard, bool keepleaves, bool fence) {
            // Must set true here so that successive calls without fence do the right thing
            this->compressed = true;
            this->nonstandard = nonstandard;
            if (world.rank() == coeffs.owner(cdata.key0)) compress_spawn(cdata.key0, nonstandard, keepleaves);
            if (fence) world.gop.fence();
        }


        // Invoked on node where key is local
        Future<tensorT> compress_spawn(const keyT& key, bool nonstandard, bool keepleaves);

        void norm_tree(bool fence) {
           if (world.rank() == coeffs.owner(cdata.key0)) norm_tree_spawn(cdata.key0);
           if (fence) world.gop.fence();
        }

        double norm_tree_op(const keyT& key, const vector< Future<double> >& v) {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            double sum = 0.0;
            int i=0;
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
                double value = v[i].get();
                sum += value*value;
            }
            sum = sqrt(sum);
            coeffs.send(key, &nodeT::set_norm_tree, sum);
            //if (key.level() == 0) std::cout << "NORM_TREE_TOP " << sum << "\n";
            return sum;
        }

        Future<double> norm_tree_spawn(const keyT& key) {
            nodeT& node = coeffs.find(key).get()->second;
            if (node.has_children()) {
                std::vector< Future<double> > v = future_vector_factory<double>(1<<NDIM);
                int i=0;
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
                    v[i] = send(coeffs.owner(kit.key()), &implT::norm_tree_spawn, kit.key());
                }
                return task(world.rank(),&implT::norm_tree_op, key, v);
            }
            else {
                return Future<double>(node.coeff().normf());
            }
        }

        tensorT compress_op(const keyT& key, const vector< Future<tensorT> >& v, bool nonstandard) {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            // Copy child scaling coeffs into contiguous block
            tensorT d(cdata.v2k);
            int i=0;
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
                d(child_patch(kit.key())) = v[i].get();
            }
            d = filter(d);
            tensorT s = copy(d(cdata.s0));
            if (key.level() > 0 && !nonstandard) d(cdata.s0) = 0.0;
            coeffs.insert(key, nodeT(d,true));
            return s;
        }

        /// Changes non-standard compressed form to standard compressed form
        void standard(bool fence) {
            for(typename dcT::iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
                const keyT& key = it->first;
                nodeT& node = it->second;
                if (key.level() > 0 && node.has_coeff()) {
                    if (node.has_children()) {
                        node.coeff()(cdata.s0) = 0.0;
                    }
                    else {
                        node.clear_coeff();
                    }
                }
            }
            if (fence) world.gop.fence();
        }

        // This routine MUST be executed in an AM handler for atomicity
        template <typename opT, typename R>
        Void recur_down_with_apply(const opT* op,
                                   const FunctionImpl<R,NDIM>* cf,
                                   const keyT& key,
                                   const keyT& target,
                                   const Tensor<R>& r) {

            PROFILE_MEMBER_FUNC(FunctionImpl);
            // We send the coeffs down in this routine so we have effectively
            // atomic insert+apply to eliminate a race condition leading to
            // double application of the operator.

            FunctionImpl<R,NDIM>* f = const_cast<FunctionImpl<R,NDIM>*>(cf);

            if (!f->coeffs.probe(key)) {
                //madness::print("not there", key);
                f->coeffs.insert(key,FunctionNode<R,NDIM>());
            }
            FunctionNode<R,NDIM>& node = f->coeffs[key];

            if (r.size) {
                // If r has data then it will be the coefficients for this node.
                // If we don't already have coeffs courtesy of someone else then
                // insert them and apply the operator.
                if (!node.has_coeff()) {
                    MADNESS_ASSERT(r.iscontiguous());
                    node.set_coeff(r);
                    //madness::print("EXTENDED APPLY", key, node.coeff().normf());
                    task(world.rank(), &implT:: template do_apply<opT,R>, op, cf, key, node.coeff());
                    if (key.level() == target.level()) return None; // Mission accomplished!
                    if (!target.is_child_of(key)) return None; // This is a sibling of the correct path
                }
            }

            // If r does not have data or we are not yet at our target then we
            // must look at the node to see what to do

            // - If key==target
            //   The coeffs should already have been made (and the operator applied)
            //   while someone else was making another node.
            //
            // - Otherwise
            // a) Node does not exist ... forward up.  Accessing the node in the manner
            //    above would have made an invalid node ... so this is captured by b)
            // b) Node exists but is invalid ... forward up (this means that someone else
            //    is already trying to make this node ... better would be to attach
            //    a callback so that when the coeffs are set this task is initiated).
            // c) Node exists and has children ... forward down
            // d) Node exists and has no children  ... recur down

            Tensor<R> empty;

            if (node.has_coeff()) { // d ... recur down if appropriate
                if (key.level() < target.level() && target.is_child_of(key)) {
                    const Tensor<R>& r = node.coeff();
                    Tensor<R> s;
                    if (r.dim[0] == k) {
                        Tensor<R> d(cdata.v2k);
                        d(cdata.s0) = node.coeff()(cdata.s0);
                        s = unfilter(d);
                    }
                    else if (r.dim[0] == 2*k) {
                        s = unfilter(node.coeff());
                    }
                    else {
                        MADNESS_EXCEPTION("Uh?",r.dim[0]);
                    }

                    node.set_has_children(true);
                    for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                        const keyT& child = kit.key();
                        Tensor<R> ss = copy(s(child_patch(child)));
                        //madness::print("EXTENDED DOWN-2", key, "to", child, ss.normf());
                        send(coeffs.owner(child), &implT:: template recur_down_with_apply<opT,R>,
                             op, cf, child, target, ss);
                    }
                }
            }
            else { // a and b ... forward up
                keyT parent = key.parent();
                //madness::print("EXTENDED UP", key, "to", parent);
                send(coeffs.owner(parent), &implT:: template recur_down_with_apply<opT,R>, op, cf, parent, target, empty);
            }

            return None;
        }


        template <typename opT, typename R>
        Void do_apply_acc(const opT* op, const FunctionImpl<R,NDIM>* f, const keyT& key, const Tensor<T>& t) {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            if (!coeffs.probe(key)) coeffs.insert(key, nodeT());
            nodeT& node = coeffs[key];

            // Accumulate into the box
            if (node.has_coeff()) {
                node.coeff().gaxpy(1.0,t,1.0);
            }
            else {
                node.set_coeff(copy(t));
                // No existing coeff and no children means the node is newly created for
                // this operation and we must tell its parent that it exists.
                if ((!node.has_children()) && (key.level() > 0)) {
                    Key<NDIM> parent = key.parent();
                    coeffs.send(parent, &FunctionNode<T,NDIM>::set_has_children_recursive, coeffs, parent);
                }

                if (op->dowiden1) {
                    typename FunctionImpl<R,NDIM>::dcT::const_iterator it = f->coeffs.find(key);
                    if ((it==f->coeffs.end() || it->second.is_invalid()) &&
                        (t.normf() > truncate_tol(thresh,key))) {
                        // We just made the first contribution to box that did not
                        // exist in the source.  Make the source box with any
                        // missing parents and apply the operator to each of them.
                        recur_down_with_apply(op, f, key.parent(), key, Tensor<R>());
                    }
                }

            }
            return None;
        }

        template <typename opT, typename R>
        Void do_apply(const opT* op, const FunctionImpl<R,NDIM>* f, const keyT& key, const Tensor<R>& c) {
            PROFILE_MEMBER_FUNC(FunctionImpl);
	    // insert timer here
	    double start_time = 0;// cpu_time();
	    double end_time = 0, cum_time = 0;
            double fac = 10.0; // 10.0 seems good for qmprop ... 3.0 OK for others
            double cnorm = c.normf();
	    const long lmax = 1L << (key.level()-1);
	    start_time = cpu_time();
            const std::vector<keyT>& disp = op->get_disp(key.level());
            for (typename std::vector<keyT>::const_iterator it=disp.begin();  it != disp.end(); ++it) {
                const keyT& d = *it;

                keyT dest = neighbor(key, d);

                // For periodic directions restrict translations to be no more than
                // half of the unit cell to avoid double counting.
                bool doit = true;
                for (int i=0; i<NDIM; i++) {
                    if (bc(i,0) == 1) {
                        if (d.translation()[i] > lmax || d.translation()[i] <= -lmax) doit = false;
                        break;
                    }
                }
                if (!doit) break;


                if (dest.is_valid()) {
                    double opnorm = op->norm(key.level(), d);
                    // working assumption here is that the operator is isotropic and
                    // montonically decreasing with distance
                    double tol = truncate_tol(thresh, key);

                    if (cnorm*opnorm > tol/fac) {
                        tensorT result = op->apply(key, d, c, tol/fac/cnorm);
                        //print("APPLY", key, d, opnorm, cnorm, result.normf());

                        // Screen here to reduce communication cost of negligible data
                        // and also to ensure we don't needlessly widen the tree when
                        // applying the operator
                        if (result.normf() > 0.3*tol/fac) {
                            //coeffs.send(dest, &nodeT::accumulate, result, coeffs, dest);
                            //madness::print("apply rrrrr       ", key, dest, result.normf());
			  //			  start_time=cpu_time();
                            send(coeffs.owner(dest), &implT:: template do_apply_acc<opT,R>, op, f, dest, result);
			    //			    end_time = cpu_time();
			    //			    cum_time += (end_time-start_time);
//                             if (op->dowiden0 && d.distsq() == 0) {
//                                 // Be sure that all touching neighbors have also applied the operator
//                                 for (int axis=0; axis<NDIM; axis++) {
//                                     for (int step=-1; step<=1; step+=2) {
//                                         keyT neigh = neighbor(key, axis, step);
//                                         if (neigh.is_valid()) {
//                                             //madness::print("checking neighbor", key, "-->", neigh);
//                                             send(coeffs.owner(neigh), &implT:: template recur_down_with_apply<opT,R>,
//                                                  op, f, neigh, neigh, Tensor<R>());
//                                         }
//                                     }
//                                 }
//                             }
                        }
                        else if (d.distsq() == 0) {
                            // If there is not a diagonal contribution there
                            // won't be off-diagonal stuff ... REALLY?????
                            break;
                        }

                    }
                    else if (d.distsq() >= 1) { // Assumes monotonic decay beyond nearest neighbor
                        break;
                    }
                }

            }
	    // update Apply_Time
	    end_time = cpu_time();
	    //	    madness::print("time for key", key, ":", end_time-start_time);
	    if (apply_time) {
	      cum_time = end_time - start_time;
	      apply_time->update(key, cum_time);
	    }
            return None;
        }

        template <typename opT, typename R>
        void apply(opT& op, const FunctionImpl<R,NDIM>& f, bool fence) {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            for(typename dcT::const_iterator it=f.coeffs.begin(); it!=f.coeffs.end(); ++it) {
                const keyT& key = it->first;
                const FunctionNode<R,NDIM>& node = it->second;
                if (node.has_coeff()) {
                    if (node.coeff().dim[0] != k || op.doleaves) {
                        ProcessID p;
                        if (FunctionDefaults<NDIM>::get_apply_randomize()) {
                            p = world.random_proc();
                        }
                        else {
                            p = coeffs.owner(key);
                        }
                        task(p, &implT:: template do_apply<opT,R>, &op, &f, key, node.coeff());
                    }
                }
            }
            if (fence) world.gop.fence();
        }

      /// accessor functions for apply_time
      // no good place to put them, so here seems as good as any
//       double get_apply_time(const keyT& key) {
// 	return apply_time.get(key);
//       }

//       void print_apply_time() {
// 	apply_time.print();
//       }

      void set_apply_time_ptr(SharedPtr<ApplyTime<NDIM> > ptr) {
	apply_time = ptr;
      }

        /// Returns the square of the error norm in the box labelled by key

        /// Assumed to be invoked locally but it would be easy to eliminate
        /// this assumption
        template <typename opT>
        double err_box(const keyT& key, const nodeT& node, const opT& func,
                       int npt, const Tensor<double>& qx, const Tensor<double>& quad_phit,
                       const Tensor<double>& quad_phiw) const {

            std::vector<long> vq(NDIM);
            for (int i=0; i<NDIM; i++) vq[i] = npt;
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
            const tensorT& coeff = node.coeff();
            ITERATOR(coeff,fval(IND)-=coeff(IND););

            // Compute the norm of what remains
            double err = fval.normf();
            return err*err;
        }

        /// Returns the sum of squares of errors from local info ... no comms
        template <typename opT>
        double errsq_local(const opT& func) const {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            // Make quadrature rule of higher order
            const int npt = cdata.npt + 1;
            Tensor<double> qx, qw, quad_phi, quad_phiw, quad_phit;
            FunctionCommonData<T,NDIM>::_init_quadrature(k+1, npt, qx, qw, quad_phi, quad_phiw, quad_phit);

            double sum = 0.0;
            for(typename dcT::const_iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
                const keyT& key = it->first;
                const nodeT& node = it->second;
                if (node.has_coeff()) sum += err_box(key, node, func, npt, qx, quad_phit, quad_phiw);
            }
            return sum;
        }


        /// Returns \c int(f(x),x) in local volume
        T trace_local() const;


        /// Returns the square of the local norm ... no comms
        double norm2sq_local() const {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            double sum = 0.0;
            for(typename dcT::const_iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
                const nodeT& node = it->second;
                if (node.has_coeff()) {
                    double norm = node.coeff().normf();
                    sum += norm*norm;
                }
            }
            return sum;
        }

	/// Returns the inner product ASSUMING same distribution
	template <typename R>
	  TENSOR_RESULT_TYPE(T,R) inner_local(const FunctionImpl<R,NDIM>& g) const {

	  TENSOR_RESULT_TYPE(T,R) sum = 0.0;
	  for(typename dcT::const_iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
	    const nodeT& fnode = it->second;
	    if (fnode.has_coeff()) {
	      if (g.coeffs.probe(it->first)) {
		const FunctionNode<R,NDIM>& gnode = g.coeffs.find(it->first).get()->second;
		if (gnode.has_coeff()) {
                  if (gnode.coeff().dim[0] != fnode.coeff().dim[0]) {
			madness::print("INNER", it->first, gnode.coeff().dim[0],fnode.coeff().dim[0]);
                        throw "adios";
                  }
		  sum += fnode.coeff().trace_conj(gnode.coeff());
		}
	      }
            }
	  }
	  return sum;
	}

        /// Returns the maximum depth of the tree
	std::size_t max_depth() const {
	    std::size_t maxdepth = 0;
	    for (typename dcT::const_iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
		std::size_t N = (std::size_t) it->first.level();
		if (N > maxdepth) maxdepth = N;
	    }
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
            for(typename dcT::const_iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
                const nodeT& node = it->second;
                if (node.has_coeff()) sum++;
            }
            if (is_compressed()) for (int i=0; i<NDIM; i++) sum *= 2*cdata.k;
            else                 for (int i=0; i<NDIM; i++) sum *=   cdata.k;

            world.gop.sum(sum);

            return sum;
        }

        /// In-place scale by a constant
        void scale_inplace(const T q, bool fence);

        /// Out-of-place scale by a constant
        template <typename Q, typename F>
        void scale_oop(const Q q, const FunctionImpl<F,NDIM>& f, bool fence) {
            typedef typename FunctionImpl<F,NDIM>::nodeT fnodeT;
            typedef typename FunctionImpl<F,NDIM>::dcT fdcT;
            for(typename fdcT::const_iterator it=f.coeffs.begin(); it!=f.coeffs.end(); ++it) {
                const keyT& key = it->first;
                const fnodeT& node = it->second;
                if (node.has_coeff()) {
                    coeffs.insert(key,nodeT(node.coeff()*q,node.has_children()));
                }
                else {
                    coeffs.insert(key,nodeT(tensorT(),node.has_children()));
                }
            }
            if (fence) world.gop.fence();
        }


    private:
        /// Assignment is not allowed ... not even possibloe now that we have reference members
        //FunctionImpl<T>& operator=(const FunctionImpl<T>& other);

        Void ensure_exists(const keyT& key);
        Void recur_down_with_fill(const keyT& target, const keyT& key);
    };

    namespace archive {
        template <class Archive, class T, int NDIM>
        struct ArchiveLoadImpl<Archive,const FunctionImpl<T,NDIM>*> {
            static void load(const Archive& ar, const FunctionImpl<T,NDIM>*& ptr) {
                uniqueidT id;
                ar & id;
                World* world = World::world_from_id(id.get_world_id());
                MADNESS_ASSERT(world);
                ptr = static_cast< const FunctionImpl<T,NDIM>* >(world->ptr_from_id< WorldObject< FunctionImpl<T,NDIM> > >(id));
                if (!ptr) MADNESS_EXCEPTION("FunctionImpl: remote operation attempting to use a locally uninitialized object",0);
            }
        };

        template <class Archive, class T, int NDIM>
        struct ArchiveStoreImpl<Archive,const FunctionImpl<T,NDIM>*> {
            static void store(const Archive& ar, const FunctionImpl<T,NDIM>*const& ptr) {
                ar & ptr->id();
            }
        };
    }



}

#endif
