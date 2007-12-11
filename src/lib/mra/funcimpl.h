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
/// \brief Provides FunctionDefaults, FunctionCommonData, FunctionImpl and FunctionFactory


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


    /// The maximum wavelet order presently supported (up to 31 should work)
    static const int MAXK = 17;

    /// FunctionDefaults holds default paramaters as static class members

    /// Declared and initialized in mra.cc and/or funcimpl::initialize.  Note that it is currently
    /// not possible to combine functions with different simulation
    /// cells and this is is not even checked for.  Therefore, the
    /// only recommended approach is to set the simulation cell once
    /// and for all at the start of a calculation and to use it for
    /// all Function instances.  A little work could improve upon this.
    ///
    /// N.B.  Ultimately, we may need to make these defaults specific to each
    /// world as should be all global state.
    template <int NDIM>
    class FunctionDefaults {
    public:
        static int k;                  ///< Wavelet order
        static double thresh;          ///< Truncation threshold
        static int initial_level;      ///< Initial level for fine scale projection
        static int max_refine_level;   ///< Level at which to stop refinement
        static int truncate_mode;    ///< Truncation method
        static bool refine;            ///< Whether to refine new functions
        static bool autorefine;        ///< Whether to autorefine in multiplication, etc.
        static bool debug;             ///< Controls output of debug info
        static Tensor<int> bc;         ///< bc[NDIM][2] Boundary conditions -- zero(0) or periodic(1)
        static Tensor<double> cell ;   ///< cell[NDIM][2] Simulation cell, cell(0,0)=xlo, cell(0,1)=xhi, ...
        static SharedPtr< WorldDCPmapInterface< Key<NDIM> > > pmap; ///< Default mapping of keys to processes

        static void set_defaults (World& world) {
            k = 7;
            thresh = 1e-5;
            initial_level = 2;
            max_refine_level = 30;
            truncate_mode = 0;
            refine = true;
            autorefine = true;
            debug = false;
            bc = Tensor<int>(NDIM,2);
            cell = Tensor<double>(NDIM,2);
            cell(_,1) = 1.0;
            //pmap = SharedPtr< WorldDCPmapInterface< Key<NDIM> > >(new WorldDCDefaultPmap< Key<NDIM> >(world));
            pmap = SharedPtr< WorldDCPmapInterface< Key<NDIM> > >(new MyPmap<NDIM>(world));
            //pmap = SharedPtr< WorldDCPmapInterface< Key<NDIM> > >(new SimpleMap< Key<NDIM> >(world));
        }
    };

    /// Holds info about displacement to neighbor for application of operators
    template <int NDIM>
    struct Displacement {
        Vector<int,NDIM> d;
        int distsq;
        Displacement() {};
        Displacement(const Vector<int, NDIM>& d) : d(d), distsq(0) {
            for (int i=0; i<NDIM; i++) distsq += d[i]*d[i];
        }

        bool operator<(const Displacement<NDIM>& other) const {
            return distsq < other.distsq;
        }

        int operator[](int i) const {return d[i];}

        template <typename Archive>
        inline void serialize(Archive& ar) {
            ar & d & distsq;
        }
    };

    template <int NDIM>
    std::ostream& operator<<(std::ostream& s, const Displacement<NDIM>& disp) {
        s << disp.d;
        return s;
    } 


    /// FunctionCommonData holds all Function data common for given k

    /// Since Function assignment and copy constructors are shallow it
    /// greatly simplifies maintaining consistent state to have all
    /// (permanent) state encapsulated in a single class.  The state
    /// is shared between instances using a SharedPtr.  Also,
    /// separating shared from instance specific state accelerates the
    /// constructor, which is important for massive parallelism and
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
        
        /// Private.  Initialize the displacements
        void _make_disp();
        
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
            _make_disp();
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

        std::vector< Displacement<NDIM> > disp; ///< Displacements in order of increasing distance

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
        T (*_f)(const coordT&);
        void (*_vf)(long, const double*, T* RESTRICT);
        int _k;
        double _thresh;
        int _initial_level;
        int _max_refine_level;
        int _truncate_mode;
        bool _refine;
        bool _empty;
        bool _autorefine;
        bool _do_fence;
        SharedPtr< WorldDCPmapInterface< Key<NDIM> > > _pmap;
        SharedPtr< FunctionFunctorInterface<T,NDIM> > _functor;
    public:
        FunctionFactory(World& world) 
            : _world(world)
            , _f(0)
            , _vf(0)
            , _k(FunctionDefaults<NDIM>::k)
            , _thresh(FunctionDefaults<NDIM>::thresh)
            , _initial_level(FunctionDefaults<NDIM>::initial_level)
            , _max_refine_level(FunctionDefaults<NDIM>::max_refine_level)
            , _truncate_mode(FunctionDefaults<NDIM>::truncate_mode)
            , _refine(FunctionDefaults<NDIM>::refine)
            , _empty(false)
            , _autorefine(FunctionDefaults<NDIM>::autorefine)
            , _do_fence(true)
            , _pmap(FunctionDefaults<NDIM>::pmap)
            , _functor(0)
        {}
        inline FunctionFactory& f(T (*f)(const coordT&)) {
            _f = f;
            return *this;
        }
        inline FunctionFactory& vf(void (*vf)(long, const double*, T* RESTRICT)) {
            _vf = vf;
            return *this;
        }
        inline FunctionFactory& k(int k) {
            _k = k;
            return *this;
        }
        inline FunctionFactory& thresh(double thresh) {
            _thresh = thresh;
            return *this;
        }
        inline FunctionFactory& initial_level(int initial_level) {
            _initial_level = initial_level;
            return *this;
        }
        inline FunctionFactory& max_refine_level(int max_refine_level) {
            _max_refine_level = max_refine_level;
            return *this;
        }
        inline FunctionFactory& truncate_mode(int truncate_mode) {
            _truncate_mode = truncate_mode;
            return *this;
        }
        inline FunctionFactory& refine(bool refine = true) {
            _refine = refine;
            return *this;
        }
        inline FunctionFactory& norefine(bool norefine = true) {
            _refine = !norefine;
            return *this;
        }
        inline FunctionFactory& empty() {
            _empty = true;
            return *this;
        }
        inline FunctionFactory& autorefine() {
            _autorefine = true;
            return *this;
        }
        inline FunctionFactory& noautorefine() {
            _autorefine = false;
            return *this;
        }
        inline FunctionFactory& do_fence(bool fence) {
            _do_fence = fence;
            return *this;
        }
        inline FunctionFactory& fence() {
            _do_fence = true;
            return *this;
        }
        inline FunctionFactory& nofence() {
            _do_fence = false;
            return *this;
        }
        inline FunctionFactory& pmap(const SharedPtr< WorldDCPmapInterface< Key<NDIM> > >& pmap) {
            _pmap = pmap;
            return *this;
        }
        inline FunctionFactory& functor(const SharedPtr< FunctionFunctorInterface<T,NDIM> >& functor) {
            _functor = functor;
            return *this;
        }
    };

    /// FunctionNode holds the coefficients, etc., at each node of the 2^NDIM-tree
    template <typename T, int NDIM>
    class FunctionNode {
    private:
        Tensor<T> _coeffs;  ///< The coefficients, if any
        bool _has_children; ///< True if there are children
        
    public:
        typedef WorldContainer< Key<NDIM>, FunctionNode<T,NDIM> > dcT;  ///< Type of container holding the nodes
        /// Default constructor makes node without coeff or children
        FunctionNode() 
            : _coeffs()
            , _has_children(false)
        {}

        /// Constructor from given coefficients with optional children

        /// Note that only a shallow copy of the coeff are taken so
        /// you should pass in a deep copy if you want the node to
        /// take ownership.
        explicit FunctionNode(const Tensor<T>& coeff, bool has_children=false) 
            : _coeffs(coeff)
            , _has_children(has_children)
        {}

        /// Type conversion of coefficients, copying all other state

        /// Choose to not overload copy and type conversion operators
        /// so there are no automatic type conversions.
        template <typename Q>
        FunctionNode<Q,NDIM> convert() const {
            return FunctionNode<Q,NDIM>(_coeffs,_has_children);
        }

        /// Returns true if there are coefficients in this node
        bool has_coeff() const {return (_coeffs.size>0);}

        /// Returns true if this node has children
        bool has_children() const {return _has_children;}

        /// Returns true if this does not have children
        bool is_leaf() const {return !_has_children;}

        /// Returns a non-const reference to the tensor containing the coeffs

        /// Returns an empty tensor if there are no coefficients.
        Tensor<T>& coeff() {return _coeffs;}

        /// Returns a const reference to the tensor containing the coeffs

        /// Returns an empty tensor if there are no coefficeints.
        const Tensor<T>& coeff() const {return _coeffs;}

        /// Sets \c has_children attribute to value of \c flag.
        Void set_has_children(bool flag) {_has_children = flag; return None;}

        /// Sets \c has_children attribute to value of \c !flag
        void set_is_leaf(bool flag) {_has_children = !flag;}

        /// Takes a \em shallow copy of the coeff --- same as \c this->coeff()=coeff
        void set_coeff(const Tensor<T>& coeff) {_coeffs = coeff;}

        /// Clears the coefficients (has_coeff() will subsequently return false)
        void clear_coeff() {_coeffs = Tensor<T>();}

        /// General bi-linear operation --- this = this*alpha + other*beta

        /// This/other may not have coefficients.  Has_children will be
        /// true in the result if either this/other have children.
        template <typename Q, typename R> 
        Void gaxpy_inplace(const T& alpha, const FunctionNode<Q,NDIM>& other, const R& beta) {
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

        /// Accumulate inplace and also connect node to parent
        Void accumulate(const Tensor<T>& t, const typename FunctionNode<T,NDIM>::dcT& c, const Key<NDIM>& key) {
            if (has_coeff()) {
                _coeffs += t;
            }
            else {
                _coeffs = copy(t);
                if (key.level() > 0) const_cast<dcT&>(c).send(key.parent(), &FunctionNode<T,NDIM>::set_has_children, true);
            }
            return None;
        }

        template <typename Archive>
        inline void serialize(Archive& ar) {
            ar & _coeffs & _has_children;
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
	friend class Function<T,NDIM>;
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
        int truncate_mode;      ///< 0=default=(|d|<thresh), 1=(|d|<thresh/2^n);
        bool autorefine;        ///< If true, autorefine where appropriate
        bool nonstandard;       ///< If true, compress keeps scaling coeff

        const FunctionCommonData<T,NDIM>& cdata;

        T (*f)(const coordT&); ///< Scalar interface to function to compress
        void (*vf)(long, const double*, T* restrict); ///< Vector interface to function to compress
        SharedPtr< FunctionFunctorInterface<T,NDIM> > functor;

        bool compressed;        ///< Compression status

        dcT coeffs;             ///< The coefficients


        // ... currently not clear how to best handle bc/cell stuff on a per function basis
        const Tensor<double> cell;///< Simulation cell (range over function is defined).
        const Tensor<int> bc;     ///< Type of boundary condition -- currently only zero or periodic
        const Tensor<double> cell_width;///< Width of simulation cell in each dimension
        Tensor<double> rcell_width; ///< Reciprocal of width
        const double cell_volume;       ///< Volume of simulation cell

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
            , nonstandard(false)
            , cdata(FunctionCommonData<T,NDIM>::get(k))
            , f(factory._f)
            , vf(factory._vf)
            , functor(factory._functor)
            , compressed(false)
            , coeffs(world,factory._pmap,false)
            , cell(FunctionDefaults<NDIM>::cell) 
            , bc(FunctionDefaults<NDIM>::bc)
            , cell_width(cell(_,1)-cell(_,0))
            , rcell_width(copy(cell_width))
            , cell_volume(cell_width.product())
        {
            // !!! Ensure that all local state is correctly formed
            // before invoking process_pending for the coeffs and
            // for this.  Otherwise, there is a race condition.
            MADNESS_ASSERT(k>0 && k<=MAXK);

            // Ultimately, must set cell, bc, etc. from the factory not the defaults
            for (int i=0; i<NDIM; i++) rcell_width(i) = 1.0/rcell_width(i);

            bool empty = factory._empty;
            bool do_refine = factory._refine;


            if (do_refine) initial_level = std::max(0,initial_level - 1);

            if (empty) {        // Do not set any coefficients at all
            }
            else if (f || vf || functor) { // Project function and optionally refine
                insert_zero_down_to_initial_level(cdata.key0);
                for(typename dcT::iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
                    if (it->second.is_leaf()) task(coeffs.owner(it->first), 
                                                   &implT::project_refine_op, 
                                                   it->first, do_refine);
                }
            }
            else {  // Set as if a zero function
                initial_level = 0;
                insert_zero_down_to_initial_level(keyT(0));
            }

            coeffs.process_pending(); 
            this->process_pending();
            
            if (factory._do_fence) world.gop.fence();
        }

        /// Copy constructor

        /// Allocates a \em new function in preparation for a deep copy
        ///
        /// By default takes pmap from other but can also specify a different pmap.
        /// Does \em not copy the coefficients ... creates an empty container.
        template <typename Q>
        FunctionImpl(const FunctionImpl<Q,NDIM>& other, 
                     const SharedPtr< WorldDCPmapInterface< Key<NDIM> > >& pmap = 0)
            : WorldObject<implT>(other.world)
            , world(other.world)
            , k(other.k)
            , thresh(other.thresh)
            , initial_level(other.initial_level)
            , max_refine_level(other.max_refine_level)
            , truncate_mode(other.truncate_mode)
            , autorefine(other.autorefine)
            , nonstandard(other.nonstandard)
            , cdata(other.cdata)
            , f(other.f)
            , vf(other.vf)
            , functor(other.functor)
            , compressed(other.compressed)
            , coeffs(world, pmap ? pmap : other.coeffs.get_pmap())
            , cell(other.cell)
            , bc(other.bc)
            , cell_width(other.cell_width)
            , rcell_width(other.rcell_width)
            , cell_volume(other.cell_volume)
        {
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


        /// Returns true if the function is compressed.  
        bool is_compressed() const {return compressed;}

        /// Adds a constant to the function.  Local operation, optional fence

        /// In scaling function basis must add value to first polyn in 
        /// each box with appropriate scaling for level.  In wavelet basis
        /// need only add at level zero.
        void add_scalar_inplace(T t, bool fence) {
            std::vector<long> v0(NDIM,0L);
            if (is_compressed()) {
                if (world.rank() == coeffs.owner(cdata.key0)) {
                    typename dcT::iterator it = coeffs.find(cdata.key0).get();
                    MADNESS_ASSERT(it != coeffs.end());
                    nodeT& node = it->second;
                    MADNESS_ASSERT(node.has_coeff());
                    node.coeff()(v0) += t*sqrt(cell_volume);
                }
            }
            else {
                for(typename dcT::iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
                    Level n = it->first.level();
                    nodeT& node = it->second;
                    if (node.has_coeff()) {
                        node.coeff()(v0) += t*sqrt(cell_volume*pow(0.5,double(NDIM*n)));
                    }
                }
            }
            if (fence) world.gop.fence();
        }
        

        /// Initialize nodes to zero function at initial_level of refinement. 

        /// Works for either basis.  No communication.
        void insert_zero_down_to_initial_level(const keyT& key) {
            if (coeffs.is_local(key)) {
                if (compressed) {
                    coeffs.insert(key, nodeT(tensorT(cdata.v2k), key.level()<initial_level));
                }
                else if (key.level()<initial_level) {
                    coeffs.insert(key, nodeT(tensorT(), true));
                }
                else {
                    coeffs.insert(key, nodeT(tensorT(cdata.vk), false));
                }
            }
            if (key.level() < initial_level) {
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    insert_zero_down_to_initial_level(kit.key());
                }
            }

        }


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
        Future<bool> truncate_spawn(const keyT& key, double tol) {
            const nodeT& node = coeffs.find(key).get()->second;   // Ugh!
            if (node.has_children()) {
                std::vector< Future<bool> > v = future_vector_factory<bool>(1<<NDIM);
                int i=0;
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
                    v[i] = send(coeffs.owner(kit.key()), &implT::truncate_spawn, kit.key(), tol);
                }
                return task(world.rank(),&implT::truncate_op, key, tol, v);
            }
            else {
                MADNESS_ASSERT(!node.has_coeff());  // In compressed form leaves should not have coeffs
                return Future<bool>(false); 
            }
        }

        /// Actually do the truncate operation
        bool truncate_op(const keyT& key, double tol, const std::vector< Future<bool> >& v) {
            // If any child has coefficients, a parent cannot truncate
            for (int i=0; i<(1<<NDIM); i++) if (v[i].get()) return true;
            nodeT& node = coeffs.find(key).get()->second;
            if (key.level() > 0) {
                double dnorm = node.coeff().normf();
                if (dnorm < truncate_tol(tol,key)) {
                    node.clear_coeff();
                    node.set_has_children(false);
                    for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                        coeffs.erase(kit.key());
                    }
                }
            }
            return node.has_coeff();
        }

        /// Evaluate function at quadrature points in the specified box
        template <typename F>
        void fcube(const keyT& key, const F& f, const Tensor<double>& qx, tensorT& fval) const {
            const Vector<Translation,NDIM>& l = key.translation();
            const Level n = key.level();
            const double h = std::pow(0.5,double(n)); 
            coordT c; // will hold the point in user coordinates
            const int npt = qx.dim[0];
            
            if (NDIM == 1) {
                for (int i=0; i<npt; i++) {
                    c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                    fval(i) = f(c);
                }
            }
            else if (NDIM == 2) {
                for (int i=0; i<npt; i++) {
                    c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                    for (int j=0; j<npt; j++) {
                        c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                        fval(i,j) = f(c);
                    }
                }
            }
            else if (NDIM == 3) {
                for (int i=0; i<npt; i++) {
                    c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                    for (int j=0; j<npt; j++) {
                        c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                        for (int k=0; k<npt; k++) {
                            c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                            fval(i,j,k) = f(c);
                        }
                    }
                }
            }
            else if (NDIM == 4) {
                for (int i=0; i<npt; i++) {
                    c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                    for (int j=0; j<npt; j++) {
                        c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                        for (int k=0; k<npt; k++) {
                            c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                            for (int m=0; m<npt; m++) {
                                c[3] = cell(3,0) + h*cell_width[3]*(l[3] + qx(m)); // xx
                                fval(i,j,k,m) = f(c);
                            }
                        }
                    }
                }
            }
            else if (NDIM == 5) {
                for (int i=0; i<npt; i++) {
                    c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                    for (int j=0; j<npt; j++) {
                        c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                        for (int k=0; k<npt; k++) {
                            c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                            for (int m=0; m<npt; m++) {
                                c[3] = cell(3,0) + h*cell_width[3]*(l[3] + qx(m)); // xx
                                for (int n=0; n<npt; n++) {
                                    c[4] = cell(4,0) + h*cell_width[4]*(l[4] + qx(n)); // yy
                                    fval(i,j,k,m,n) = f(c);
                                }
                            }
                        }
                    }
                }
            }
            else if (NDIM == 6) {
                for (int i=0; i<npt; i++) {
                    c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                    for (int j=0; j<npt; j++) {
                        c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                        for (int k=0; k<npt; k++) {
                            c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                            for (int m=0; m<npt; m++) {
                                c[3] = cell(3,0) + h*cell_width[3]*(l[3] + qx(m)); // xx
                                for (int n=0; n<npt; n++) {
                                    c[4] = cell(4,0) + h*cell_width[4]*(l[4] + qx(n)); // yy
                                    for (int p=0; p<npt; p++) {
                                        c[5] = cell(5,0) + h*cell_width[5]*(l[5] + qx(p)); // zz
                                        fval(i,j,k,m,n,p) = f(c);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else {
                MADNESS_EXCEPTION("FunctionImpl: fcube: confused about NDIM?",NDIM);
            }
        }

        const keyT& key0() const {
            return cdata.key0;
        }

        void print_tree(Level maxlevel = 10000) const {
            if (world.rank() == 0) do_print_tree(cdata.key0, maxlevel);
            world.gop.fence();
            if (world.rank() == 0) std::cout.flush();
            world.gop.fence();
        }

        void do_print_tree(const keyT& key, Level maxlevel) const {
            typename dcT::iterator it = coeffs.find(key).get();
            if (it == const_cast<dcT*>(&coeffs)->end()) {
                MADNESS_EXCEPTION("FunctionImpl: do_print_tree: null node pointer",0);
            }
            const nodeT& node = it->second;
            for (int i=0; i<key.level(); i++) std::cout << "  ";
            std::cout << key << "  " << node << " --> " << coeffs.owner(key) << "\n";
            if (key.level() < maxlevel  &&  node.has_children()) {
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    do_print_tree(kit.key(),maxlevel);
                }
            }
        }


        /// Compute by projection the scaling function coeffs in specified box
        tensorT project(const keyT& key) const {
            tensorT fval(cdata.vq,false); // this will be the returned result
            tensorT& work = cdata.work1; // initially evaluate the function in here

            if (vf) {
                MADNESS_EXCEPTION("FunctionImpl: project: Vector function interface not yet implemented",0);
            }
            else if (f) {
                fcube(key,f,cdata.quad_x,work);
            }
            else if (functor) {
                fcube(key,*functor,cdata.quad_x,work);
            }
            else {
                MADNESS_EXCEPTION("FunctionImpl: project: confusion about function?",0);
            }

            work.scale(sqrt(cell_volume*pow(0.5,double(NDIM*key.level()))));
            //return transform(work,cdata.quad_phiw);
            return fast_transform(work,cdata.quad_phiw,fval,cdata.workq);
        }


        /// Returns the truncation threshold according to truncate_method
        inline double truncate_tol(double tol, const keyT& key) const {
            if (truncate_mode == 0) {
                return tol;
            }
            else if (truncate_mode == 1) {
                return tol*pow(0.5,double(key.level()));
            }
            else if (truncate_mode == 2) {
                return tol*pow(0.5,0.5*key.level()*NDIM);
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
        Void project_refine_op(const keyT& key, bool do_refine) {
            if (do_refine) {
                // Make in r child scaling function coeffs at level n+1
                tensorT r(cdata.v2k);
                for (KeyChildIterator<NDIM> it(key); it; ++it) {
                    const keyT& child = it.key();
                    r(child_patch(child)) = project(child);
                }

                // Insert empty node for parent
                coeffs.insert(key,nodeT(tensorT(),true));

                // Filter then test difference coeffs at level n
                tensorT d = filter(r);
                d(cdata.s0) = T(0);
                if ( (d.normf() < truncate_tol(thresh,key.level())) && 
                     (key.level() < max_refine_level) ) {
                    for (KeyChildIterator<NDIM> it(key); it; ++it) {
                        const keyT& child = it.key();
                        coeffs.insert(child,nodeT(copy(r(child_patch(child))),false));
                    }
                }
                else {
                    if (key.level()==max_refine_level) print("MAX REFINE LEVEL",key);
                    for (KeyChildIterator<NDIM> it(key); it; ++it) {
                        const keyT& child = it.key();
                        //ProcessID p = coeffs.owner(child);
                        ProcessID p = world.random_proc();
                        task(p, &implT::project_refine_op, child, do_refine); // ugh
                    }
                }
            }
            else {
                coeffs.insert(key,nodeT(project(key),false));
            }
            return None;
        }


        /// Compute the Legendre scaling functions for multiplication

        /// Evaluate parent polyn at quadrature points of a child.  The prefactor of 
        /// 2^n/2 is included.  The tensor must be preallocated as phi(k,npt).
        /// Refer to the implementation notes for more info.
        void phi_for_mul(Level np, Translation lp, Level nc, Translation lc, Tensor<double>& phi) const {
            double p[200];
            double scale = pow(2.0,double(np-nc));
            for (int mu=0; mu<cdata.npt; mu++) {
                double xmu = scale*(cdata.quad_x(mu)+lc) - lp;
                MADNESS_ASSERT(xmu>-1e-15 && xmu<(1+1e-15));
                legendre_scaling_functions(xmu,cdata.k,p);
                for (int i=0; i<k; i++) phi(i,mu) = p[i];
            }
            phi.scale(pow(2.0,0.5*np));
        }

        /// Directly project parent coeffs to child coeffs

        /// Currently used by diff, but other uses can be anticipated
        const tensorT parent_to_child(const tensorT& s, const keyT& parent, const keyT& child) const {
            // An invalid parent/child means that they are out of the box
            // and it is the responsibility of the caller to worry about that
            // ... most likely the coefficiencts (s) are zero to reflect
            // zero B.C. so returning s makes handling this easy.
            if (parent == child || parent.is_invalid() || child.is_invalid()) return s;
            
            tensorT result = fcube_for_mul<T>(child, parent, s);
            result.scale(sqrt(cell_volume*pow(0.5,double(NDIM*child.level()))));
            result = transform(result,cdata.quad_phiw);

//             if (parent.level()+1 == child.level()) {
//                 tensorT d(cdata.v2k);
//                 d(cdata.s0) = s;
//                 d = unfilter(d);
//                 //print("dumb result\n",d(child_patch(child)));
//                 //print("clev result\n",result);
//                 double err = (d(child_patch(child))-result).normf();
//                 if (err > 1e-12) {
//                     print("bad parent_to_child", parent, child, err, " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
//                 }
//             }

            return result;
        }


        // Testing.

        // This routine was implemented to help test fcube_for_mul.  Instead
        // of computing the function values in a child cell directly from
        // the parent coefficients it first recurs the parent coefficients down one
        // level and then uses the standard transformation to generate the
        // values.  All looks good!
        template <typename Q>
        Tensor<Q> fcube_for_mul2(const keyT& child, const keyT& parent, const Tensor<Q> coeff) const {
            MADNESS_ASSERT(child.level()-parent.level() == 1);
            Tensor<Q> w(cdata.v2k);
            w(cdata.s0) = coeff(___);
            w = unfilter(w);
            w = copy(w(child_patch(child)));

            double scale = pow(2.0,0.5*NDIM*child.level())/sqrt(cell_volume);
            return transform(w,cdata.quad_phit).scale(scale);
        }


        /// Compute the function values for multiplication

        /// Given coefficients from a parent cell, compute the value of
        /// the functions at the quadrature points of a child
        template <typename Q>
        Tensor<Q> fcube_for_mul(const keyT& child, const keyT& parent, const Tensor<Q>& coeff) const {
            if (child.level() == parent.level()) {
                double scale = pow(2.0,0.5*NDIM*parent.level())/sqrt(cell_volume);
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
                return general_transform(coeff,phi).scale(1.0/sqrt(cell_volume));;
            }
        }
        
        /// Invoked as a task by mul with the actual coefficients
        template <typename L, typename R>
        Void do_mul(const keyT& key, const Tensor<L>& left, const std::pair< keyT, Tensor<R> >& arg) {
            const keyT& rkey = arg.first;
            const Tensor<R>& rcoeff = arg.second;
            Tensor<R> rcube = fcube_for_mul(key, rkey, rcoeff);
            Tensor<L> lcube = fcube_for_mul(key,  key, left);

            Tensor<T> tcube(cdata.vk,false);
            TERNARY_OPTIMIZED_ITERATOR(T, tcube, L, lcube, R, rcube, *_p0 = *_p1 * *_p2;);
            double scale = pow(0.5,0.5*NDIM*key.level())*sqrt(cell_volume);
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



        /// Invoked by result to compute the pointwise product result=left*right

        /// This version requires all three functions have the same distribution.
        /// Should be straightforward to do an efficient version that does not
        /// require this but I have not thought about that yet.
        ///
        /// Possible non-blocking communication and optional fence.
        template <typename L, typename R>
        void mul(const FunctionImpl<L,NDIM>& left, const FunctionImpl<R,NDIM>& right, bool fence) {
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
                        const_cast<FunctionImpl<R,NDIM>*>(&right)->
                            send(coeffs.owner(parent), &FunctionImpl<R,NDIM>::sock_it_to_me, 
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
                        const_cast<FunctionImpl<L,NDIM>*>(&left)->
                            send(coeffs.owner(parent), &FunctionImpl<L,NDIM>::sock_it_to_me, 
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
            //print_tree(2);
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
        }
        


        /// Verify tree is properly constructed ... global synchronization involved

        /// If an inconsistency is detected, prints a message describing the error and 
        /// then throws a madness exception.
        ///
        /// This is a reasonably quick and scalable operation that is
        /// useful for debugging and paranoia.
        void verify_tree() const {
            world.gop.fence();  // Make sure nothing is going on

            // Verify consistency of compression status, existence and size of coefficients,
            // and has_children() flag.
            for(typename dcT::const_iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
                const keyT& key = it->first;
                const nodeT& node = it->second;
                bool bad;

                if (is_compressed()) {
                    if (node.has_children()) {
                        bad = node.coeff().dim[0] != 2*cdata.k;
                    }
                    else {
                        bad = node.coeff().size != 0;
                    }
                }
                else {
                    if (node.has_children()) {
                        bad = node.coeff().size != 0;
                    }
                    else {
                        bad = node.coeff().dim[0] != cdata.k;
                    }
                }

                if (bad) {
                    print(world.rank(), "FunctionImpl: verify: INCONSISTENT TREE NODE, key =", key, ", node =", node,
                          ", dim[0] =",node.coeff().dim[0],", compressed =",is_compressed());
                    std::cout.flush();
                    MADNESS_EXCEPTION("FunctionImpl: verify: INCONSISTENT TREE NODE", 0);
                }
            }
                            
            // Ensure that parents and children exist appropriately
            for(typename dcT::const_iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
                const keyT& key = it->first;
                const nodeT& node = it->second;

                // NOTE MESSED UP CONSTNESS HERE DUE TO INCOMPLETE
                // INTERFACE TO DC ... NEEDS FIXING IN THE DC CLASS.
                if (key.level() > 0) {
                    const keyT parent = key.parent();
                    typename dcT::iterator pit = coeffs.find(parent).get();
                    if (pit == const_cast<dcT*>(&coeffs)->end()) {
                        print(world.rank(), "FunctionImpl: verify: MISSING PARENT",key,parent);
                        std::cout.flush();
                        MADNESS_EXCEPTION("FunctionImpl: verify: MISSING PARENT", 0);
                    }
                    const nodeT& pnode = pit->second;
                    if (!pnode.has_children()) {
                        print(world.rank(), "FunctionImpl: verify: PARENT THINKS IT HAS NO CHILDREN",key,parent);
                        std::cout.flush();
                        MADNESS_EXCEPTION("FunctionImpl: verify: PARENT THINKS IT HAS NO CHILDREN", 0);
                    }
                }

                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    typename dcT::iterator cit = coeffs.find(kit.key()).get();
                    if (cit == const_cast<dcT*>(&coeffs)->end()) {
                        if (node.has_children()) {
                            print(world.rank(), "FunctionImpl: verify: MISSING CHILD",key,kit.key());
                            std::cout.flush();
                            MADNESS_EXCEPTION("FunctionImpl: verify: MISSING CHILD", 0);
                        }
                    }
                    else {
                        if (! node.has_children()) {
                            print(world.rank(), "FunctionImpl: verify: UNEXPECTED CHILD",key,kit.key());
                            std::cout.flush();
                            MADNESS_EXCEPTION("FunctionImpl: verify: UNEXPECTED CHILD", 0);
                        }
                    }
                }
            }

            world.gop.fence();
        }


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
                           const RemoteReference< FutureImpl< std::pair<keyT,tensorT> > >& ref) {
            if (coeffs.probe(key)) {
                const nodeT& node = coeffs.find(key).get()->second;
                Future< std::pair<keyT,tensorT> > result(ref);
                if (node.has_coeff()) {
                    //madness::print("sock found it with coeff",key);
                    result.set(std::pair<keyT,tensorT>(key,node.coeff()));
                }
                else {
                    //madness::print("sock found it without coeff",key);
                    result.set(std::pair<keyT,tensorT>(key,tensorT()));
                }
            }
            else {
                keyT parent = key.parent();
                //madness::print("sock forwarding to parent",key,parent);
                send(coeffs.owner(parent), &FunctionImpl<T,NDIM>::sock_it_to_me, parent, ref);
            }
            return None;
        }


        /// Evaluate the function at a point in \em simulation coordinates

        /// Only the invoking process will get the result via the
        /// remote reference to a future.  Active messages may be sent
        /// to other nodes.
        Void eval(const Vector<double,NDIM>& xin, 
                  const keyT& keyin, 
                  const typename Future<T>::remote_refT& ref) {

            // This is ugly.  We must figure out a clean way to use 
            // owner computes rule from the container.
            Vector<double,NDIM> x = xin;
            keyT key = keyin;
            Vector<Translation,NDIM> l = key.translation();
            ProcessID me = world.rank();
            while (1) {
                ProcessID owner = coeffs.owner(key);
                if (owner != me) {
                    send(owner, &implT::eval, x, key, ref);
                    return None;
                }
                else {
                    typename dcT::futureT fut = coeffs.find(key);
                    typename dcT::iterator it = fut.get();
                    nodeT& node = it->second;
                    if (node.has_coeff()) {
                        Future<T>(ref).set(eval_cube(key.level(), x, node.coeff()));
                        return None;
                    }
                    else {
                        for (int i=0; i<NDIM; i++) {
                            double xi = x[i]*2.0;
                            int li = int(xi);
                            if (li == 2) li = 1;
                            x[i] = xi - li;
                            l[i] = 2*l[i] + li;
                        }
                        key = keyT(key.level()+1,l);
                    }
                }
            }
            //MADNESS_EXCEPTION("should not be here",0);
        }


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
        void tnorm(const tensorT& t, double* lo, double* hi) const {
            // Chosen approach looks stupid but it is more accurate
            // than the simple approach of summing everything and
            // subtracting off the low-order stuff to get the high
            // order (assuming the high-order stuff is small relative
            // to the low-order)
            cdata.work1(___) = t(___);
            tensorT tlo = cdata.work1(cdata.sh);
            *lo = tlo.normf();
            tlo.fill(0.0);
            *hi = cdata.work1.normf();
        }


        // This invoked if node has not been autorefined
        Void do_square_inplace(const keyT& key) {
            Level n = key.level();
            nodeT& node = coeffs[key];

            double scale = pow(2.0,0.5*NDIM*n)/sqrt(cell_volume);
            tensorT tval = transform(node.coeff(),cdata.quad_phit).scale(scale);

            tval.emul(tval);
            
            scale = pow(0.5,0.5*NDIM*n)*sqrt(cell_volume);
            tval = transform(tval,cdata.quad_phiw).scale(scale);

            node.set_coeff(tval);

            return None;
        }


        // This invoked if node has been autorefined
        Void do_square_inplace2(const keyT& parent, const keyT& child, const tensorT& parent_coeff) {
            tensorT tval = fcube_for_mul (child, parent, parent_coeff);

            tval.emul(tval);
            
            Level n = child.level();
            double scale = pow(0.5,0.5*NDIM*n)*sqrt(cell_volume);
            tval = transform(tval,cdata.quad_phiw).scale(scale);

            coeffs.insert(child,nodeT(tval,false));

            return None;
        }


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
        void square_inplace(bool fence) {
            for(typename dcT::iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
                const keyT& key = it->first;
                nodeT& node = it->second;
                if (node.has_coeff()) {
                    if (autorefine && autorefine_square_test(key, node.coeff())) {
                        tensorT t = node.coeff();
                        node.clear_coeff();
                        node.set_has_children(true);
                        for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {         
                            const keyT& child = kit.key();
                            task(coeffs.owner(child), &implT::do_square_inplace2, key, child, t);
                        }
                    }
                    else {
                        task(world.rank(), &implT::do_square_inplace, key);
                    }
                }
            }
            if (fence) world.gop.fence();
        }

        /// Differentiation of function with optional global fence

        /// Result of differentiating f is placed into this which will
        /// have the same process map, etc., as f
        void diff(const implT& f, int axis, bool fence) {
            typedef std::pair<keyT,tensorT> argT;
            for(typename dcT::const_iterator it=f.coeffs.begin(); it!=f.coeffs.end(); ++it) {
                const keyT& key = it->first;
                const nodeT& node = it->second;
                if (node.has_coeff()) {
                    Future<argT> left   = f.find_neighbor(key,axis,-1);
                    argT center(key,node.coeff());
                    Future<argT> right  = f.find_neighbor(key,axis, 1);
                    task(world.rank(), &implT::do_diff1, &f, axis, key, left, center, right, task_attr_generator());
                }
                else {
                    // Internal empty node can be safely inserted
                    coeffs.insert(key,nodeT(tensorT(),true));
                }
            }
            if (fence) world.gop.fence();
        }
  

        /// Returns key of neighbor enforcing BC

        /// Out of volume keys are mapped to enforce the BC as follows.
        ///   * Periodic BC map back into the volume and return the correct key
        ///   * Zero BC - returns invalid() to indicate out of volume
        keyT neighbor(const keyT& key, int axis, int step) const {
            Vector<Translation,NDIM> l = key.translation();
            const Translation two2n = Translation(1)<<key.level();

            // Translation is an unsigned type so need to be careful if results
            // or intermediates are possibly negative.
            if (step < 0  &&  l[axis] < (unsigned)(-step)) {
                if (bc(axis,0) == 0) {
                    return keyT::invalid(); // Zero BC
                }
                else if (bc(axis,0) == 1) {
                    l[axis] = (l[axis] + two2n) + step; // Periodic BC
                }
                else {
                    MADNESS_EXCEPTION("neighbor: confused left BC?",bc(axis,0));
                }
            }
            else if (l[axis]+step >= two2n) {
                if (bc(axis,1) == 0) {
                    return keyT::invalid(); // Zero BC
                }
                else if (bc(axis,1) == 1) {
                    l[axis] = (l[axis] + step) - two2n; // Periodic BC
                }
                else {
                    MADNESS_EXCEPTION("neighbor: confused BC right?",bc(axis,1));
                }
            }
            else {
                l[axis] += step;
            }
      
            return keyT(key.level(),l);
        }


        /// Returns key of general neighbor enforcing BC

        /// Out of volume keys are mapped to enforce the BC as follows.
        ///   * Periodic BC map back into the volume and return the correct key
        ///   * Zero BC - returns invalid() to indicate out of volume
        keyT neighbor(const keyT& key, const Displacement<NDIM>& d) const {
            Translation two2n = 1ul << key.level();
            Vector<Translation,NDIM> l = key.translation();
            for (int axis=0; axis<NDIM; axis++) {
                int step = d[axis];
                if (step < 0  &&  l[axis] < (unsigned)(-step)) {
                    if (bc(axis,0) == 0) {
                        return keyT::invalid(); // Zero BC
                    }
                    else if (bc(axis,0) == 1) {
                        l[axis] = (l[axis] + two2n) + step; // Periodic BC
                    }
                    else {
                        MADNESS_EXCEPTION("neighbor: confused left BC?",bc(axis,0));
                    }
                }
                else if (l[axis]+step >= two2n) {
                    if (bc(axis,1) == 0) {
                        return keyT::invalid(); // Zero BC
                    }
                    else if (bc(axis,1) == 1) {
                        l[axis] = (l[axis] + step) - two2n; // Periodic BC
                    }
                    else {
                        MADNESS_EXCEPTION("neighbor: confused BC right?",bc(axis,1));
                    }
                }
                else {
                    l[axis] += step;
                }
            }
            return keyT(key.level(),l);
        }


        /// Called by diff to find key and coeffs of neighbor enforcing BC

        /// Should work for any (small) step but only tested for step=+/-1
        ///
        /// do_diff1 handles the adpative refinement.  If it needs to refine, it calls
        /// forward_do_diff1 to pass the task locally or remotely with high priority.
        /// Actual differentiation is performend by do_diff2.
        Future< std::pair<keyT,tensorT> > find_neighbor(const keyT& key, int axis, int step) const {
            typedef std::pair<keyT,tensorT> argT;
            keyT neigh = neighbor(key, axis, step);
            if (neigh.is_invalid()) {
                return Future<argT>(argT(neigh,tensorT(cdata.vk))); // Zero bc
            }
            else {
                Future<argT> result;
                send(coeffs.owner(neigh), &implT::sock_it_to_me, neigh, result.remote_ref(world));
                return result;
            }
        }


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
                      const std::pair<keyT,tensorT>& right) {
            ProcessID owner = coeffs.owner(key);
            if (owner == world.rank()) {
                if (left.second.size == 0) {
                    task(owner, &implT::do_diff1, f, axis, key, 
                         f->find_neighbor(key,axis,-1), center, right, 
                         task_attr_generator());
                }
                else if (right.second.size == 0) {
                    task(owner, &implT::do_diff1, f, axis, key, 
                         left, center, f->find_neighbor(key,axis,1), 
                         task_attr_generator());
                }
                else {
                    task(owner, &implT::do_diff2, f, axis, key, 
                         left, center, right);
                }
            }
            else {
                send(owner, &implT::forward_do_diff1, f, axis, key, 
                     left, center, right);
            }
            return None;
        }


        Void do_diff1(const implT* f, int axis, const keyT& key, 
                      const std::pair<keyT,tensorT>& left,
                      const std::pair<keyT,tensorT>& center,
                      const std::pair<keyT,tensorT>& right) {
            typedef std::pair<keyT,tensorT> argT;

            MADNESS_ASSERT(axis>=0 && axis<NDIM);

            if (left.second.size==0 || right.second.size==0) {
                // One of the neighbors is below us in the tree ... recur down
                coeffs.insert(key,nodeT(tensorT(),true));
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT& child = kit.key();
                    if ((child.translation()[axis]&1) == 0) { 
                        // leftmost child automatically has right sibling
                        forward_do_diff1(f, axis, child, left, center, center);
                    }
                    else { 
                        // rightmost child automatically has left sibling
                        forward_do_diff1(f, axis, child, center, center, right);
                    }
                }
            }
            else {
                forward_do_diff1(f, axis, key, left, center, right);
            }
            return None;
        }

        Void do_diff2(const implT* f, int axis, const keyT& key, 
                      const std::pair<keyT,tensorT>& left,
                      const std::pair<keyT,tensorT>& center,
                      const std::pair<keyT,tensorT>& right) {
            typedef std::pair<keyT,tensorT> argT;
            
            tensorT d = madness::inner(cdata.rp, 
                                       parent_to_child(left.second, left.first, neighbor(key,axis,-1)).swapdim(axis,0), 
                                       1, 0);
            inner_result(cdata.r0, 
                         parent_to_child(center.second, center.first, key).swapdim(axis,0),
                         1, 0, d);
            inner_result(cdata.rm, 
                         parent_to_child(right.second, right.first, neighbor(key,axis,1)).swapdim(axis,0), 
                         1, 0, d);
            if (axis) d = copy(d.swapdim(axis,0)); // make it contiguous
            d.scale(rcell_width[axis]*pow(2.0,(double) key.level()));
            coeffs.insert(key,nodeT(d,false));
            return None;
        }

        /// Permute the dimensions according to map
        void mapdim(const implT& f, const std::vector<long>& map, bool fence) {
            for(typename dcT::const_iterator it=f.coeffs.begin(); it!=f.coeffs.end(); ++it) {
                const keyT& key = it->first;
                const nodeT& node = it->second;

                Vector<Translation,NDIM> l;
                for (int i=0; i<NDIM; i++) l[map[i]] = key.translation()[i];
                tensorT c = node.coeff();
                if (c.size) c = copy(c.mapdim(map));

                coeffs.insert(keyT(key.level(),l), nodeT(c,node.has_children()));
            }
            if (fence) world.gop.fence();
        }


        T eval_cube(Level n, coordT x, const tensorT c) const {
            const int k = cdata.k;
            double px[NDIM][k];
            T sum = T(0.0);

            for (int i=0; i<NDIM; i++) legendre_scaling_functions(x[i],k,px[i]);

            if (NDIM == 1) {
                for (int p=0; p<k; p++) 
                    sum += c(p)*px[0][p];
            }
            else if (NDIM == 2) {
                for (int p=0; p<k; p++) 
                    for (int q=0; q<k; q++) 
                        sum += c(p,q)*px[0][p]*px[1][q];
            }
            else if (NDIM == 3) {
                for (int p=0; p<k; p++) 
                    for (int q=0; q<k; q++) 
                        for (int r=0; r<k; r++) 
                        sum += c(p,q,r)*px[0][p]*px[1][q]*px[2][r];
            }
            else if (NDIM == 4) {
                for (int p=0; p<k; p++) 
                    for (int q=0; q<k; q++) 
                        for (int r=0; r<k; r++) 
                            for (int s=0; s<k; s++) 
                                sum += c(p,q,r,s)*px[0][p]*px[1][q]*px[2][r]*px[3][s];
            }
            else if (NDIM == 5) {
                for (int p=0; p<k; p++) 
                    for (int q=0; q<k; q++) 
                        for (int r=0; r<k; r++) 
                            for (int s=0; s<k; s++) 
                                for (int t=0; t<k; t++) 
                                sum += c(p,q,r,s,t)*px[0][p]*px[1][q]*px[2][r]*px[3][s]*px[4][t];
            }
            else if (NDIM == 6) {
                for (int p=0; p<k; p++) 
                    for (int q=0; q<k; q++) 
                        for (int r=0; r<k; r++) 
                            for (int s=0; s<k; s++) 
                                for (int t=0; t<k; t++) 
                                    for (int u=0; u<k; u++) 
                                    sum += c(p,q,r,s,t,u)*px[0][p]*px[1][q]*px[2][r]*px[3][s]*px[4][t]*px[5][u];
            }
            else {
                MADNESS_EXCEPTION("FunctionImpl:eval_cube:NDIM?",NDIM);
            }
            return sum*pow(2.0,0.5*NDIM*n)/sqrt(cell_volume);
        }


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
        inline tensorT filter(const tensorT& s) const {
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
        inline tensorT unfilter(const tensorT& s) const {
            tensorT r(cdata.v2k,false);
            return fast_transform(s,cdata.hg,r,cdata.work2);
            //return transform(s, cdata.hg);
        }

        void reconstruct(bool fence) {
            if (world.rank() == coeffs.owner(cdata.key0)) reconstruct_op(cdata.key0,tensorT());
            if (fence) world.gop.fence();
            compressed = false;
        }

        // Invoked on node where key is local
        Void reconstruct_op(const keyT& key, const tensorT& s) {
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
                coeffs.insert(key,nodeT(tensorT(),false));
                it = coeffs.find(key).get();
            }
            nodeT& node = it->second;
            
            // The integral operator will correctly connect interior nodes
            // to children but may leave interior nodes without coefficients
            // ... but they still need to sum down so just give them zeros
            if (node.has_children() && !node.has_coeff()) node.set_coeff(tensorT(cdata.v2k));

            if (node.has_coeff()) {
                tensorT d = node.coeff();
                if (key.level() > 0) d(cdata.s0) += s; // -- note accumulate for NS summation
                d = unfilter(d);
                node.clear_coeff();
                node.set_has_children(true);
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT& child = kit.key();
                    tensorT ss = copy(d(child_patch(child)));
                    task(coeffs.owner(child), &implT::reconstruct_op, child, ss);
                }
            }
            else {
                node.set_coeff(s);
            }
            return None;
        }
        
        void compress(bool nonstandard, bool fence) {
           if (world.rank() == coeffs.owner(cdata.key0)) compress_spawn(cdata.key0, nonstandard);
           if (fence) world.gop.fence();
           compressed = true;
        }


        // Invoked on node where key is local
        Future<tensorT> compress_spawn(const keyT& key, bool nonstandard) {
            MADNESS_ASSERT(coeffs.probe(key));
            nodeT& node = coeffs.find(key).get()->second;
            if (node.has_children()) {
                std::vector< Future<tensorT> > v = future_vector_factory<tensorT>(1<<NDIM);
                int i=0;
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
                    v[i] = send(coeffs.owner(kit.key()), &implT::compress_spawn, kit.key(), nonstandard);
                }
                return task(world.rank(),&implT::compress_op, key, v, nonstandard);
            }
            else {
                Future<tensorT> result(node.coeff());
                if (!nonstandard) node.clear_coeff();
                return result;
            }
        }

        tensorT compress_op(const keyT& key, const vector< Future<tensorT> >& v, bool nonstandard) {
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


        template <typename opT, typename R>
        Void do_apply(const opT* op, const keyT& key, const Tensor<R>& c) {
            // This scaling of the norm accounts for the number of 
            // possible contributions coming into a region which is roughly
            // the number of neighbors times the number of levels.
            //
            // In 1D estimate nneigh=10, 2D 5^2, 3D and higher 3^NDIM
            //
            // Not yet computing actual tree depth of the tree.  Estimate as 10.
            int nmax = 10;
            int nneigh = std::max(10/NDIM, 3);
            double fac = std::pow(double(nneigh),-1.0*NDIM) / nmax;
            double cnorm = fac*c.normf();
            for (typename std::vector< Displacement<NDIM> >::const_iterator it=cdata.disp.begin(); 
                 it != cdata.disp.end(); 
                 ++it) {
                const Displacement<NDIM>& d = *it;

                keyT dest = neighbor(key, d);
                if (dest.is_valid()) {
                    double opnorm = op->norm(key, d);
                    
                    // working assumption here is that the operator is isotropic and
                    // montonically decreasing with distance ... this needs to be validated.
                    
                    if (cnorm*opnorm < truncate_tol(thresh, key)) break;
                    
                    tensorT result = op->apply(key, d, c);
                    
                    coeffs.send(dest, &nodeT::accumulate, result, coeffs, dest);
                }

            }
            return None;
        }

        template <typename opT, typename R>
        void apply(opT& op, const FunctionImpl<R,NDIM>& f, bool fence) {
            for(typename dcT::const_iterator it=f.coeffs.begin(); it!=f.coeffs.end(); ++it) {
                const keyT& key = it->first;
                const FunctionNode<R,NDIM>& node = it->second;
                // SCREEN HERE ?
                task(world.rank(), &implT:: template do_apply<opT,R>, &op, key, node.coeff());
            }
            if (fence) world.gop.fence();
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
            double scale = pow(0.5,0.5*NDIM*key.level())*sqrt(cell_volume);
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

        /// Returns the square of the local norm ... no comms
        double norm2sq_local() const {
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


        /// Convert user coords (cell[][]) to simulation coords ([0,1]^ndim)
        inline void user_to_sim(const coordT& xuser, coordT& xsim) const {
            for (int i=0; i<NDIM; i++) 
                xsim[i] = (xuser[i] - cell(i,0)) * rcell_width[i];
        }


        /// Convert simulation coords ([0,1]^ndim) to user coords (cell[][])
        inline void sim_to_user(const coordT& xsim, coordT& xuser) const {
            for (int i=0; i<NDIM; i++)
                xuser[i] = xsim[i]*cell_width[i] + cell(i,0);
        }

        /// In-place scale by a constant
        template <typename Q>
        void scale_inplace(const Q q, bool fence) {
            for(typename dcT::iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
                nodeT& node = it->second;
                if (node.has_coeff()) node.coeff().scale(q);
            }
            if (fence) world.gop.fence();
        }

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
    };

    namespace archive {
        template <class Archive, class T, int NDIM>
        struct ArchiveLoadImpl<Archive,const FunctionImpl<T,NDIM>*> {
            static inline void load(const Archive& ar, const FunctionImpl<T,NDIM>*& ptr) {
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
            static inline void store(const Archive& ar, const FunctionImpl<T,NDIM>*const& ptr) {
                ar & ptr->id();
            }
        };
    }
    


}

#endif
