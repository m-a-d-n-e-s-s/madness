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
    template <typename T, int D> class LoadBalImpl;
    template <int D> class LBTree;
    template <int D> class MyPmap;
}


namespace madness {


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
        static bool compress;          ///< Whether to compress new functions
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
            compress = false;
            refine = true;
            autorefine = true;
            debug = false;
            bc = Tensor<int>(NDIM,2);
            cell = Tensor<double>(NDIM,2);
            cell(_,1) = 1.0;
            //pmap = SharedPtr< WorldDCPmapInterface< Key<NDIM> > >(new WorldDCDefaultPmap< Key<NDIM> >(world));
            pmap = SharedPtr< WorldDCPmapInterface< Key<NDIM> > >(new MyPmap<NDIM>(world));
        };
    };

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
        
        /// Private.  Initialize the quadrature information
        void _init_quadrature();

        /// Private.  Do first use initialization
        void _initialize(int k) {
            print("initializing commondata for",k);
            this->k = k;
            npt = k;
            for (int i=0; i<4; i++) 
                s[i] = Slice(i*k,(i+1)*k-1);
            s0 = std::vector<Slice>(NDIM);
            vk = std::vector<long>(NDIM);
            vq = std::vector<long>(NDIM);
            v2k = std::vector<long>(NDIM);
            for (int i=0; i<NDIM; i++) {
                s0[i] = s[0];
                vk[i] = k;
                vq[i] = npt;
                v2k[i] = 2*k;
            }
            work1 = tensorT(vk);
            work2 = tensorT(v2k);
            workq = tensorT(vq);
            zero_coeff = tensorT(vk);
            key0 = Key<NDIM>(0,Vector<Translation,NDIM>(0));

            _init_twoscale();
            _init_quadrature();
            _make_dc_periodic();
            initialized = true;
        };

        FunctionCommonData() 
            : initialized(false) 
        {};
        
        bool initialized;
    public:
        typedef Tensor<T> tensorT; ///< Type of tensor used to hold coeff

        int k;                  ///< order of the wavelet
        int npt;                ///< no. of quadrature points
        Slice s[4];             ///< s[0]=Slice(0,k-1), s[1]=Slice(k,2*k-1), etc.
        std::vector<Slice> s0;  ///< s[0] in each dimension to get scaling coeff
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
            print("common data getting",k,data[k].initialized);
            if (!data[k].initialized) data[k]._initialize(k);
            return data[k];
        };
    };

    /// FunctionFactory implements the named-parameter idiom for Function

    /// C++ does not provide named arguments (as does, e.g., Python).
    /// This class provides something very close.  Create functions as follows
    /// \code
    /// double myfunc(const double x[]);
    /// Function<double,3> f = FunctionFactory<double,3>(world).f(myfunc).k(11).thresh(1e-9).debug().nocompress()
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
        bool _compress;
        bool _refine;
        bool _empty;
        bool _autorefine;
        SharedPtr< WorldDCPmapInterface< Key<NDIM> > > _pmap;
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
            , _compress(FunctionDefaults<NDIM>::compress)
            , _refine(FunctionDefaults<NDIM>::refine)
            , _empty(false)
            , _autorefine(FunctionDefaults<NDIM>::autorefine)
            , _pmap(FunctionDefaults<NDIM>::pmap)
        {};
        inline FunctionFactory& f(T (*f)(const coordT&)) {
            _f = f;
            return *this;
        };
        inline FunctionFactory& vf(void (*vf)(long, const double*, T* RESTRICT)) {
            _vf = vf;
            return *this;
        };
        inline FunctionFactory& k(int k) {
            _k = k;
            return *this;
        };
        inline FunctionFactory& thresh(double thresh) {
            _thresh = thresh;
            return *this;
        };
        inline FunctionFactory& initial_level(int initial_level) {
            _initial_level = initial_level;
            return *this;
        };
        inline FunctionFactory& max_refine_level(int max_refine_level) {
            _max_refine_level = max_refine_level;
            return *this;
        };
        inline FunctionFactory& truncate_mode(int truncate_mode) {
            _truncate_mode = truncate_mode;
            return *this;
        };
        inline FunctionFactory& compress(bool compress = true) {
            _compress = compress;
            return *this;
        };
        inline FunctionFactory& nocompress(bool nocompress = true) {
            _compress = !nocompress;
            return *this;
        };
        inline FunctionFactory& refine(bool refine = true) {
            _refine = refine;
            return *this;
        };
        inline FunctionFactory& norefine(bool norefine = true) {
            _refine = !norefine;
            return *this;
        };
        inline FunctionFactory& empty() {
            _empty = true;
            return *this;
        };
        inline FunctionFactory& autorefine() {
            _autorefine = true;
            return *this;
        };
        inline FunctionFactory& noautorefine() {
            _autorefine = false;
            return *this;
        };
        inline FunctionFactory& pmap(const SharedPtr< WorldDCPmapInterface< Key<NDIM> > >& pmap) {
            _pmap = pmap;
        };
    };

    /// FunctionNode holds the coefficients, etc., at each node of the 2^NDIM-tree
    template <typename T, int NDIM>
    class FunctionNode {
    private:
        Tensor<T> _coeffs;  ///< The coefficients, if any
        bool _has_children; ///< True if there are children
        
    public:
        /// Default constructor makes node without coeff or children
        FunctionNode() 
            : _coeffs()
            , _has_children(false)
        {};

        /// Constructor from given coefficients with optional children

        /// Note that only a shallow copy of the coeff are taken so
        /// you should pass in a deep copy if you want the node to
        /// take ownership.
        FunctionNode(const Tensor<T>& coeff, bool has_children=false) 
            : _coeffs(coeff)
            , _has_children(has_children)
        {};

        /// Returns true if there are coefficients in this node
        bool has_coeff() const {return (_coeffs.size>0);};

        /// Returns true if this node has children
        bool has_children() const {return _has_children;};

        /// Returns true if this does not have children
        bool is_leaf() const {return !_has_children;};

        /// Returns a non-const reference to the tensor containing the coeffs

        /// Returns an empty tensor if there are no coefficeints.
        Tensor<T>& coeff() {return _coeffs;};

        /// Returns a const reference to the tensor containing the coeffs

        /// Returns an empty tensor if there are no coefficeints.
        const Tensor<T>& coeff() const {return _coeffs;};

        /// Sets \c has_children attribute to value of \c flag.
        void set_has_children(bool flag) {_has_children = flag;};

        /// Sets \c has_children attribute to value of \c !flag
        void set_is_leaf(bool flag) {_has_children = !flag;};

        /// Takes a \em shallow copy of the coeff --- same as \c this->coeff()=coeff
        void set_coeff(const Tensor<T>& coeff) {_coeffs = coeff;}

        /// Clears the coefficients (has_coeff() will subsequently return false)
        void clear_coeff() {_coeffs = Tensor<T>();};

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
        s << norm;
        return s;
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
	friend class LoadBalImpl<T,NDIM>;
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
        double truncate_thr;    ///< Tolerance for truncation, defaults to thresh
        double autorefine_thr;  ///< Tolerance for autorefine, defaults to thresh        
        int initial_level;      ///< Initial level for refinement
        int max_refine_level;   ///< Do not refine below this level
        int truncate_mode;      ///< 0=default=(|d|<thresh), 1=(|d|<thresh/2^n);
        bool autorefine;        ///< If true, autorefine where appropriate
        bool dorefine;          ///< If true, refine when constructed
        bool nonstandard;       ///< If true, compress keeps scaling coeff

        const FunctionCommonData<T,NDIM>& cdata;

        T (*f)(const coordT&); ///< Scalar interface to function to compress
        void (*vf)(long, const double*, T* restrict); ///< Vector interface to function to compress

        bool compressed;        ///< Compression status
        long nterminated;       ///< No. of boxes where adaptive refinement was too deep

        dcT coeffs;             ///< The coefficients


        // ... currently not clear how to best handle bc/cell stuff on a per function basis
        const Tensor<double> cell;///< Simulation cell (range over function is defined).
        const Tensor<int> bc;     ///< Type of boundary condition -- currently only zero or periodic
        const Tensor<double> cell_width;///< Width of simulation cell in each dimension
        Tensor<double> rcell_width; ///< Reciprocal of width
        const double cell_volume;       ///< Volume of simulation cell

    public:

        /// Initialize function impl from data in factory
        FunctionImpl(const FunctionFactory<T,NDIM>& factory) 
            : WorldObject<implT>(factory._world)
            , world(factory._world)
            , k(factory._k)
            , thresh(factory._thresh)
            , truncate_thr(factory._thresh)
            , autorefine_thr(factory._thresh)
            , initial_level(factory._initial_level)
            , max_refine_level(factory._max_refine_level)
            , truncate_mode(factory._truncate_mode)
            , autorefine(factory._autorefine)
            , dorefine(factory._refine)
            , nonstandard(false)
            , cdata(FunctionCommonData<T,NDIM>::get(k))
            , f(factory._f)
            , vf(factory._vf)
            , compressed(false)
            , nterminated(0)
            , coeffs(world,factory._pmap,false)
            , cell(FunctionDefaults<NDIM>::cell) 
            , bc(FunctionDefaults<NDIM>::bc)
            , cell_width(cell(_,1)-cell(_,0))
            , rcell_width(copy(cell_width))
            , cell_volume(cell_width.product())
        {
            MADNESS_ASSERT(k>0 && k<MAXK);

            // Ultimately, must set cell, bc, etc. from the factory not the defaults
            for (int i=0; i<NDIM; i++) rcell_width(i) = 1.0/rcell_width(i);

            bool empty = factory._empty;
            bool do_compress = factory._compress;
            if (dorefine) initial_level = std::max(0,initial_level - 1);

            if (empty) {        // Do not set any coefficients at all
                compressed = do_compress;
                // No need to process pending for coeffs
            }
            else if (f || vf) { // Project function and optionally compress
                compressed = false;
                insert_zero_down_to_initial_level(cdata.key0);
                print("Tree after zero down");
                print_tree();
                for(typename dcT::iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
                    print("iterating",it->first,it->second.is_leaf());
                    if (it->second.is_leaf()) task(coeffs.owner(it->first), 
                                                   &implT::project_refine_op, 
                                                   it->first);
                };
                coeffs.process_pending(); 
                world.gop.fence(); // !!! Ultimately this must be optional
                //if (do_compress) compress();
            }
            else {  // Set as if a zero function
                initial_level = 0;
                compressed = do_compress;
                insert_zero_down_to_initial_level(keyT(0));
                // No need to process pending
            };
        };

        /// Copy constructor

        /// Allocates a \em new function in preparation for a deep copy
        ///
        /// By default takes pmap from other but can also specify a different pmap.
        /// Does \em not copy the coefficients ... creates an empty container.
        FunctionImpl(const FunctionImpl<T,NDIM>& other, 
                     const SharedPtr< WorldDCPmapInterface< Key<NDIM> > >& pmap = 0)
            : WorldObject<implT>(other.world)
            , world(other.world)
            , k(other.k)
            , thresh(other.thresh)
            , truncate_thr(other.truncate_thr)
            , autorefine_thr(other.autorefine_thr)
            , initial_level(other.initial_level)
            , max_refine_level(other.max_refine_level)
            , truncate_mode(other.truncate_mode)
            , autorefine(other.autorefine)
            , dorefine(other.dorefine)
            , nonstandard(other.nonstandard)
            , cdata(other.cdata)
            , f(other.f)
            , vf(other.vf)
            , compressed(other.compressed)
            , nterminated(other.nterminated)
            , coeffs(world, pmap ? pmap : other.coeffs.get_pmap())
            , cell(other.cell)
            , bc(other.bc)
            , cell_width(other.cell_width)
            , rcell_width(other.rcell_width)
            , cell_volume(other.cell_volume)
        {};

	const SharedPtr< WorldDCPmapInterface< Key<NDIM> > >& get_pmap() const {
	    return coeffs.get_pmap();
	};

	void copy_coeffs(const implT& other) {
	    for(typename dcT::const_iterator it=other.coeffs.begin(); it!=other.coeffs.end(); ++it) {
		coeffs.insert(*it);
	    };
	    world.gop.fence();
	};

        /// Returns true if the function is compressed.  
        bool is_compressed() const {return compressed;};

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

        };


        /// Evaluate function at quadrature points in the specified box
        tensorT fcube(const keyT& key) const {
            MADNESS_ASSERT(NDIM == 3);
            MADNESS_ASSERT(f);

            tensorT fval(cdata.vq);
            const Vector<Translation,NDIM>& l = key.translation();
            const Level n = key.level();
            const double h = std::pow(0.5,double(n)); 
            const int npt = cdata.npt;
            coordT c; // will hold the point in user coordinates
            
            for (int i=0; i<npt; i++) {
                c[0] = cell(i,0) + h*cell_width[0]*(l[0] + cdata.quad_x(i)); // x
                for (int j=0; j<npt; j++) {
                    c[1] = cell(j,0) + h*cell_width[1]*(l[1] + cdata.quad_x(j)); // y
                    for (int k=0; k<npt; k++) {
                        c[2] = cell(k,0) + h*cell_width[2]*(l[2] + cdata.quad_x(k)); // z
                        fval(i,j,k) = f(c);
                    }
                }
            }

            return fval;
        };

        const keyT& key0() const {
            return cdata.key0;
        };

        void print_tree() const {
            if (world.rank() == 0) do_print_tree(cdata.key0);
            world.gop.fence();
        };

        void do_print_tree(const keyT& key) const {
            const nodeT& node = coeffs.find(key).get()->second;
            for (int i=0; i<key.level(); i++) std::cout << "  ";
            std::cout << key << "  " << node << "\n";
            if (node.has_children()) {
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    do_print_tree(kit.key());
                };
            }
        };


        /// Compute by projection the scaling function coeffs in specified box
        tensorT project(const keyT& key) const {
            tensorT fval;

            if (vf) 
                MADNESS_ASSERT(0); //fval = vfcube(key);
            else
                fval = fcube(key);

            fval.scale(sqrt(cell_volume*pow(0.5,double(NDIM*key.level()))));
            return transform(fval,cdata.quad_phiw);
        };


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
            };
        };


        /// Returns patch refering to coeffs of child in parent box

        /// !! Returns a reference to \em STATIC data !!
        /// Can be optimized away by precomputing.
        const std::vector<Slice>& child_patch(const keyT& child) const {
            static std::vector<Slice> s(NDIM);
            const Vector<Translation,NDIM>& l = child.translation();
            for (int i=0; i<NDIM; i++) s[i] = cdata.s[l[i]&1]; // Lowest bit of translation
            return s;
        };


        /// Projection with optional refinement
        Void project_refine_op(const keyT& key) {
            //print("DOING PROJECT REFINE",key);

            if (dorefine) {
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
                        task(coeffs.owner(child), &implT::project_refine_op, child); // ugh
                    }
                }
            }
            else {
                coeffs.insert(key,nodeT(project(key),false));
            }
            return None;
        };


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
                    print("got node",it->first,it->second.is_leaf(),it->second.has_coeff());
                    //nodeT& node = coeffs.find(key).get()->second;
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
            MADNESS_EXCEPTION("should not be here",0);
        };


        T eval_cube(Level n, coordT x, const tensorT s) const {
            double px[64], py[64], pz[64];
            MADNESS_ASSERT(k<64);
            MADNESS_ASSERT(NDIM==3);
            legendre_scaling_functions(x[0],k,px);
            legendre_scaling_functions(x[1],k,py);
            legendre_scaling_functions(x[2],k,pz);
            T sum = T(0.0);
            for (int p=0; p<k; p++) {
                for (int q=0; q<k; q++) {
                    for (int r=0; r<k; r++) {
                        sum += s(p,q,r)*px[p]*py[q]*pz[r];
                    }
                }
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
           return transform(s, cdata.hgT);
        };


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
        inline tensorT unfilter(const tensorT& ss,
                                bool sonly = false) const {
            if (sonly)
                //return transform(ss,cdata.hgsonly);
                MADNESS_EXCEPTION("unfilter: sonly : not yet",0);
            else {
                return transform(ss, cdata.hg);
            }
        };

        void reconstruct(bool fence) {
            if (world.rank() == coeffs.owner(cdata.key0)) reconstruct_op(cdata.key0,tensorT());
            if (fence) world.gop.fence();
            compressed = false;
        };

        // Invoked on node where key is local
        Void reconstruct_op(const keyT& key, const tensorT& s) {
            nodeT& node = coeffs.find(key).get()->second;
            if (node.has_coeff()) {
                tensorT d = node.coeff();
                if (key.level() > 0) d(cdata.s0) = s;
                d = unfilter(d);
                node.clear_coeff();
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT& child = kit.key();
                    task(coeffs.owner(child), &implT::reconstruct_op, child, copy(d(child_patch(child))));
                }
            }
            else {
                node.set_coeff(s);
            }
            return None;
        };
        
        void compress(bool fence) {
           if (world.rank() == coeffs.owner(cdata.key0)) compress_spawn(cdata.key0);
           if (fence) world.gop.fence();
           compressed = true;
        };


        // Invoked on node where key is local
        Future<tensorT> compress_spawn(const keyT& key) {
            nodeT& node = coeffs.find(key).get()->second;
            if (node.has_children()) {
                std::vector< Future<tensorT> > v = future_vector_factory<tensorT>(1<<NDIM);
                int i=0;
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
                    v[i] = send(coeffs.owner(kit.key()), &implT::compress_spawn, kit.key());
                }
                return task(world.rank(),&implT::compress_op, key, v);
            }
            else {
                Future<tensorT> result(node.coeff());
                node.clear_coeff();
                return result;
            }
        };

        tensorT compress_op(const keyT& key, const vector< Future<tensorT> >& v) {
            // Copy child scaling coeffs into contiguous block
            tensorT d(cdata.v2k);
            int i=0;
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
                d(child_patch(kit.key())) = v[i].get();
            };
            d = filter(d);
            tensorT s = copy(d(cdata.s0));
            if (key.level() > 0) d(cdata.s0) = 0.0;
            coeffs.insert(key, nodeT(d,true));
            return s;
        };


        /// Convert user coords (cell[][]) to simulation coords ([0,1]^ndim)
        inline void user_to_sim(const coordT& xuser, coordT& xsim) const {
            for (int i=0; i<NDIM; i++)
                xsim[i] = (xuser[i] - cell(i,0)) * rcell_width[i];
        };


        /// Convert simulation coords ([0,1]^ndim) to user coords (cell[][])
        inline void sim_to_user(const coordT& xsim, coordT& xuser) const {
            for (int i=0; i<NDIM; i++)
                xuser[i] = xsim[i]*cell_width[i] + cell(i,0);
        };

    private:
        /// Assignment is not allowed ... not even possibloe now that we have reference members
        //FunctionImpl<T>& operator=(const FunctionImpl<T>& other);
    };


}

#endif
