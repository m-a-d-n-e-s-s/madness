#ifndef MAD_FUNC_DATA
#define MAD_FUNC_DATA

/// \file funcimpl.h
/// \brief Provides FunctionDefaults, FunctionCommonData, FunctionImpl and FunctionFactory


namespace madness {
    /// FunctionDefaults holds default paramaters as static class members

    /// Declared and initialized in mra.cc.  Note that it is currently
    /// not possible to combine functions with different simulation
    /// cells and this is is not even checked for.  Therefore, the
    /// only recommended approach is to set the simulation cell once
    /// and for all at the start of a calculation and to use it for
    /// all Function instances.  A little work could improve upon this.
    template <int NDIM>
    class FunctionDefaults {
    public:
        static int k;                  ///< Wavelet order
        static double thresh;          ///< Truncation threshold
        static int initial_level;      ///< Initial level for fine scale projection
        static int max_refine_level;   ///< Level at which to stop refinement
        static int truncate_method;    ///< Truncation method
        static bool compress;          ///< Whether to compress new functions
        static bool refine;            ///< Whether to refine new functions
        static bool autorefine;        ///< Whether to autorefine in multiplication, etc.
        static bool debug;             ///< Controls output of debug info
        static int bc[NDIM][2];        ///< Type of boundary condition -- currently only zero or periodic
        static double cell[NDIM][2];   ///< Simulation cell, cell[0][0]=xlo, cell[0][1]=xhi, ...
    };

    /// FunctionCommonData holds all Function data common for given k

    /// Since Function assignment and copy constructors are shallow it
    /// greatly simplifies maintaining consistent state to have all
    /// (permanent) state encapsulated in a single class.  The state
    /// is shared between instances using a SharedPtr.  Also,
    /// separating shared from instance specific state accelerates the
    /// constructor, which is important for massive parallelism and
    /// permitting inexpensive use of temporaries.  The default copy
    /// constructor and assignment operator are used.
    template <typename T, int NDIM>
    class FunctionCommonData {
    private:
        /// Private.  Make the level-0 blocks of the periodic central difference derivative operator
        void _make_dc_periodic();

        /// Private.  Initialize the twoscale coefficients
        void _init_twoscale();

        /// Private.  Initialize the quadrature information
        void _init_quadrature();

    public:
        int k;                  ///< Wavelet order
        int npt;                ///< no. of quadrature points

        Slice s[4];              ///< s[0]=Slice(0,k-1), s[1]=Slice(k,2*k-1), etc.
        std::vector<Slice> s0;  ///< s[0] in each dimension to get scaling coeffs
        std::vector<long> vk;   ///< (k,...) used to initialize Tensors
        std::vector<long> v2k;  ///< (2k,...) used to initialize Tensors
        std::vector<long> vq;   ///< (npt,...) used to initialize Tensors

        typedef Tensor<T> tensorT; ///< Type of tensor used to hold coeffs
        mutable Tensor<T> work1;        ///< work space of size (k,...)
        mutable Tensor<T> work2;        ///< work space of size (2k,...)
        mutable Tensor<T> workq;        ///< work space of size (npt,...)
        Tensor<T> zero_tensor;  ///< Zero (k,...) tensor for internal convenience of diff
        Key<NDIM> key0; ///< Key for root node

        Tensor<double> quad_x;  ///< quadrature points
        Tensor<double> quad_w;  ///< quadrature weights
        Tensor<double> quad_phi; ///< quad_phi(i,j) = at x[i] value of phi[j]
        Tensor<double> quad_phit; ///< transpose of quad_phi
        Tensor<double> quad_phiw; ///< quad_phiw(i,j) = at x[i] value of w[i]*phi[j]

        Tensor<double> h0, h1, g0, g1; ///< The separate blocks of twoscale coefficients
        Tensor<double> hg, hgT; ///< The full twoscale coeffs (2k,2k) and transpose
        Tensor<double> hgsonly; ///< hg[0:k,:]

        Tensor<double> rm, r0, rp;        ///< Blocks of the derivative operator
        Tensor<double> rm_left, rm_right, rp_left, rp_right; ///< Rank-1 forms rm & rp

        /// Default constructor necessary for array/vector construction
        FunctionCommonData() {};

        /// Constructor presently forces choice npt=k
        FunctionCommonData(int k) : k(k), npt(k) {
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
            zero_tensor = tensorT(vk);
            key0 = Key<NDIM>(0,Array<Translation,NDIM>(0));

            _init_twoscale();
            _init_quadrature();
            _make_dc_periodic();
        };

    };
    

    template <typename T, int NDIM, typename Pmap> class FunctionImpl;
    template <typename T, int NDIM, typename Pmap> class Function;


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
    template <typename T, int NDIM, typename Pmap=DCDefaultProcmap<Key<NDIM> > > class FunctionFactory {
        friend class FunctionImpl<T,NDIM,Pmap>;
    protected:
        World& _world;
        T (*_f)(const double[NDIM]);
        void (*_vf)(long, const double*, T* restrict);
        int _k;
        double _thresh;
        int _initial_level;
        int _max_refine_level;
        int _truncate_method;
        bool _compress;
        bool _refine;
        bool _empty;
        bool _autorefine;
    public:
        FunctionFactory(World& world) 
            : _world(world)
            , _f(0)
            , _vf(0)
            , _k(FunctionDefaults<NDIM>::k)
            , _thresh(FunctionDefaults<NDIM>::thresh)
            , _initial_level(FunctionDefaults<NDIM>::initial_level)
            , _max_refine_level(FunctionDefaults<NDIM>::max_refine_level)
            , _truncate_method(FunctionDefaults<NDIM>::truncate_method)
            , _compress(FunctionDefaults<NDIM>::compress)
            , _refine(FunctionDefaults<NDIM>::refine)
            , _empty(false)
            , _autorefine(FunctionDefaults<NDIM>::autorefine)
        {};
        inline FunctionFactory& f(T (*f)(const double[NDIM])) {
            _f = f;
            return *this;
        };
        inline FunctionFactory& vf(void (*vf)(long, const double*, T* restrict)) {
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
        inline FunctionFactory& truncate_method(int truncate_method) {
            _truncate_method = truncate_method;
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
    };

    /// FunctionNode holds the coefficients, etc., at each node of the 2^NDIM-tree
    template <typename T, int NDIM>
    struct FunctionNode {
        Tensor<T> coeffs;  ///< The coefficients, if any
        bool has_coeff;    ///< True if there are coefficients
        bool has_children; ///< True if there are children

        /// Default constructor makes node without coeffs or children
        FunctionNode() 
            : coeffs()
            , has_coeff(false)
            , has_children(false)
        {};

        /// Construct from given coefficients with optional children
        FunctionNode(const Tensor<T>& coeffs, bool has_children=false) 
            : coeffs(coeffs)
            , has_coeff(true)
            , has_children(has_children)
        {};

        template <typename Archive>
        inline void serialize(Archive& ar) {
            ar & coeffs & has_coeff & has_children;
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
    template <typename T, int NDIM, typename Pmap=DCDefaultProcmap<Key<NDIM> > >
    class FunctionImpl : public WorldContainer< Key<NDIM>, FunctionNode<T,NDIM>, Pmap > {
    private:
        static const int MAXK = 17;
        static FunctionCommonData<T,NDIM> commondata[MAXK + 1]; ///< Declared in mra.cc
        static bool initialized;	///< Declared and initialized to false in mra.cc

    public:
        typedef Tensor<T> tensorT;                     ///< Type of tensor used to hold coeffs
        typedef Array<Translation,NDIM> tranT;         ///< Type of array holding translation
        typedef Key<NDIM> keyT;                        ///< Type of key
        typedef FunctionNode<T,NDIM> nodeT;            ///< Type of node
        typedef WorldContainer<keyT,nodeT, Pmap> dcT;  ///< Type of container holding the coefficients

        World& world;
        int k;                  ///< Wavelet order
        double thresh;          ///< Screening threshold
        double truncate_thr;    ///< Tolerance for truncation, defaults to thresh
        double autorefine_thr;  ///< Tolerance for autorefine, defaults to thresh        
        int initial_level;      ///< Initial level for refinement
        int max_refine_level;   ///< Do not refine below this level
        int truncate_method;    ///< 0=default=(|d|<thresh), 1=(|d|<thresh/2^n);
        bool autorefine;        ///< If true, autorefine where appropriate
        bool dorefine;          ///< If true, refine when constructed
        bool nonstandard;       ///< If true, compress keeps scaling coeff

        const FunctionCommonData<T,NDIM>* cdata;

        T (*f)(const double[NDIM]); ///< Scalar interface to function to compress
        void (*vf)(long, const double*, T* restrict); ///< Vector interface to function to compress

        bool compressed;        ///< Compression status
        long nterminated;       ///< No. of boxes where adaptive refinement was too deep

        // ... currently not clear how to best handle bc/cell stuff on a per function basis
        double cell_volume;       ///< Volume of simulation cell
        double cell[NDIM][2];     ///< Simulation cell (range over function is defined).
        double cell_width[NDIM];  ///< Width of simulation cell in each dimension
        double rcell_width[NDIM]; ///< Reciprocal of width
        int bc[NDIM][2];          ///< Type of boundary condition -- currently only zero or periodic


        /// Initialize function impl from data in factory
        FunctionImpl(const FunctionFactory<T,NDIM,Pmap>& factory) 
            : dcT(factory._world)
            , world(factory._world)
            , k(factory._k)
            , thresh(factory._thresh)
            , truncate_thr(factory._thresh)
            , autorefine_thr(factory._thresh)
            , initial_level(factory._initial_level)
            , max_refine_level(factory._max_refine_level)
            , truncate_method(factory._truncate_method)
            , autorefine(factory._autorefine)
            , dorefine(factory._refine)
            , nonstandard(false)
            , cdata(commondata+factory._k)
            , f(factory._f)
            , vf(factory._vf)
            , compressed(false)
            , nterminated(0)
        {
            MADNESS_ASSERT(k>0 && k<MAXK);
            if (!initialized) {
                FunctionImpl<T,NDIM,Pmap>::initialize();
                cdata = commondata+k;
            }

            // Ultimately these need to be set from the factory not the defaults
            cell_volume = 1.0;
            for (int i=0; i<NDIM; i++) {
                bc[i][0] = FunctionDefaults<NDIM>::bc[i][0];
                bc[i][1] = FunctionDefaults<NDIM>::bc[i][1];
                cell[i][0] = FunctionDefaults<NDIM>::cell[i][0];
                cell[i][1] = FunctionDefaults<NDIM>::cell[i][1];
                cell_width[i] = cell[i][1] - cell[i][0];
                rcell_width[i] = 1.0/cell_width[i];
                cell_volume *= cell_width[i];
            }

//             bool empty = factory._empty;

//             if (f || vf) {
//                 if (world.rank() == 0) {
//                     world.taskq.add(coeffs.owner(cdata->key0),
//                                     *this, &FunctionImpl<T,NDIM>::project_refine, cdata->key0);
//                 }
//             } else if (empty) {             // Do not set any coefficients at all
//                 //impl->compressed = compress;
//             } else {                        // Set coefficients so have zero function
//                 if (world.rank() == 0) {
//                     coeffs.insert(cdata->key0, nodeT(tensorT(cdata->vk)));
//                 }
//             }
//             world.gop.fence();  // Ultimately, this must be optimized away
        };

        void project_refine(const keyT& key) {
            print("DOING IT!");
        };


        /// Copy constructor

        /// Allocates a \em new function index in preparation for a deep copy
        ///
        /// !!! DOES NOT COPY THE COEFFICIENTS !!!
        FunctionImpl(const FunctionImpl<T,NDIM,Pmap>& other) 
            : dcT(other.world)
            , world(other.world)
            , k(other.k)
            , thresh(other.thresh)
            , truncate_thr(other.truncate_thr)
            , autorefine_thr(other.autorefine_thr)
            , initial_level(other.initial_level)
            , max_refine_level(other.max_refine_level)
            , truncate_method(other.truncate_method)
            , autorefine(other.autorefine)
            , dorefine(other.dorefine)
            , nonstandard(other.nonstandard)
            , cdata(other.cdata)
            , f(other.f)
            , vf(other.vf)
            , compressed(other.compressed)
            , nterminated(other.nterminated)
        { 
            cell_volume = other.cell_volume;
            for (int i=0; i<NDIM; i++) {
                bc[i][0] = other.bc[i][0];
                bc[i][1] = other.bc[i][1];
                cell[i][0] = other.cell[i][0];
                cell[i][1] = other.cell[i][1];
                cell_width[i] = other.cell_width[i];
                rcell_width[i] = other.rcell_width[i];
            }
        };


        /// Convert user coords (cell[][]) to simulation coords ([0,1])
        inline void user_to_sim(double x[NDIM]) const {
            for (int i=0; i<NDIM; i++)
                x[i] = (x[i] - cell[i][0]) * rcell_width[i];
        };


        /// Convert simulation coords ([0,1]) to user coords (cell[][])
        inline void sim_to_user(double x[NDIM]) const {
            for (int i=0; i<NDIM; i++)
                x[i] = x[i]*cell_width[i] + cell[i][0];
        };

    private:
        /// Initialize static data
        static void initialize() {
            for (int k = 1; k <= MAXK; k++) commondata[k] = FunctionCommonData<T,NDIM>(k);
            initialized = true;
        };

        /// Assignment is not allowed ... not at all now we have reference member
        //FunctionImpl<T>& operator=(const FunctionImpl<T>& other);
    };
}

#endif
