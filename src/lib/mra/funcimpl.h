#ifndef MAD_FUNC_DATA
#define MAD_FUNC_DATA

/// \file funcimpl.h
/// \brief Provides FunctionDefaults, FunctionCommonData, FunctionImpl and FunctionFactory


namespace madness {
    /// FunctionDefaults holds default paramaters as static class members

    /// Defined and initialized in mra.cc.  Note that it is currently
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
            cell_volume = 1.0;
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

            _init_twoscale();
            _init_quadrature();
            _make_dc_periodic();
        };

    };
    

    template <typename T, int NDIM> class FunctionImpl;
    template <typename T, int NDIM> class Function;


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
    template <typename T, int NDIM> class FunctionFactory {
        friend class FunctionImpl<T,NDIM>;
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
        inline FunctionFactory& f(void (*f)(const double[NDIM], T* restrict)) {
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
    template <typename T, int NDIM>
    class FunctionImpl {
    private:
        static const int MAXK = 17;
        static FunctionCommonData<T,NDIM> commondata[MAXK + 1]; ///< Defined in mra.cc
        static bool initialized;	///< Defined and initialized to false in mra.cc

    public:
        typedef Tensor<T> tensorT;                     ///< Type of tensor used to hold coeffs
        typedef Key<NDIM> keyT;                        ///< Type of key
        typedef FunctionNode<T,NDIM> nodeT;            ///< Type of node
        typedef Array<Translation,NDIM> tranT;         ///< Type of array holding translation
        typedef DistributedContainer<keyT,nodeT> dcT;  ///< Type of container holding the coefficients

        World& world;
        int k;                  ///< Wavelet order
        double thresh;          ///< Screening threshold
        double truncate_thr;    ///< Tolerance for truncation, defaults to thresh
        double autorefine_thr;  ///< Tolerance for autorefine, defaults to thresh        
        int initial_level;      ///< Initial level for refinement
        int max_refine_level;   ///< Do not refine below this level
        int truncate_method;    ///< 0=default=(|d|<thresh), 1=(|d|<thresh/2^n);
        bool autorefine;        ///< If true, autorefine where appropriate
        bool refine;            ///< If true, refine when constructed
        bool nonstandard;        ///< If true, compress keeps scaling coeff

        const FunctionCommonData<T,NDIM>* cdata;

        T (*f)(const double[NDIM]); ///< Scalar interface to function to compress
        void (*vf)(long, const double*, T* restrict); ///< Vector interface to function to compress

        bool compressed;        ///< Compression status
        long nterminated;       ///< No. of boxes where adaptive refinement was too deep

        //typedef DistributedContainer< Key<NDIM>, FunctionNode<T,NDIM>, FunctionProcMap<NDIM> > dcT;
        dcT coeffs;


        // ... currently not clear how to best handle cell stuff on a per function basis
        double cell[NDIM][2];   ///< Simulation cell (range over function is defined).
        double cell_width[NDIM];///< Width of simulation cell in each dimension
        double rcell_width[NDIM];///< Reciprocal of width
        double cell_volume;     ///< Volume of simulation cell


        /// Initialize function impl from data in factory
        FunctionImpl(const FunctionFactory<T,NDIM>& factory) 
            : world(factory._world)
            , k(factory._k)
            , thresh(factory._thresh)
            , truncate_thr(factory._thresh)
            , autorefine_thr(factory._thresh)
            , initial_level(factory._initial_level)
            , max_refine_level(factory._max_refine_level)
            , truncate_method(factory._truncate_method)
            , autorefine(factory._autorefine)
            , refine(factory._refine)
            , nonstandard(false)
            , cdata(commondata+factory._k)
            , f(factory._f)
            , vf(factory._vf)
            , compressed(false)
            , nterminated(0)
            , coeffs(factory._world)
        {
            MADNESS_ASSERT(k>0 && k<MAXK);
            // Ultimately this needs to be set from the factory not the defaults
            cell_volume = 1.0;
            for (int i=0; i<NDIM; i++) {
                cell[i][0] = FunctionDefaults<NDIM>::cell[i][0];
                cell[i][1] = FunctionDefaults<NDIM>::cell[i][1];
                cell_width[i] = cell[i][1] - cell[i][0];
                rcell_width[i] = 1.0/cell_width[i];
                cell_volume *= cell_width[i];
            if (!initialized) {
                this->initialize();
                cdata = commondata+k;
            }

            bool compress = factory._compress;
            bool empty = factory._empty;

            if (f || vf) {
                //project_refine();
                //if (compress) this->compress();
            } else if (empty) {             // Do not set any coefficients at all
                //impl->compressed = compress;
            } else {                        // Set coefficients so have zero function
                compressed = compress;
                if (world.rank() == 0) {
                    if (compress) 
                        coeffs.insert(keyT(0,tranT(0)), nodeT(tensorT(cdata->v2k)));
                    else    
                        coeffs.insert(keyT(0,tranT(0)), nodeT(tensorT(cdata->vk)));
                }
            }
        };

        /// Copy constructor

        /// Allocates a \em new function index in preparation for a deep copy
        FunctionImpl(const FunctionImpl<T,NDIM>& other) 
            : world(other.world)
            , k(other.k)
            , thresh(other.thresh)
            , truncate_thr(other.truncate_thr)
            , autorefine_thr(other.autorefine_thr)
            , initial_level(other.initial_level)
            , max_refine_level(other.max_refine_level)
            , truncate_method(other.truncate_method)
            , autorefine(other.autorefine)
            , refine(other.refine)
            , nonstandard(other.nonstandard)
            , cdata(other.cdata)
            , f(other.f)
            , vf(other.vf)
            , compressed(other.compressed)
            , nterminated(other.nterminated)
            , coeffs(other.coeffs)
        { };


        /// Convert user coords (cell[][]) to simulation coords ([0,1])
        inline void user_to_sim(double x[NDIM]) const {
            for (int i=0; i<NDIM; i++)
                x[i] = (x[i] - cell[i][0]) * rcell_width[0];
        };


        /// Convert simulation coords ([0,1]) to user coords (cell[][])
        inline void sim_to_user(double x[NDIM]) const {
            for (int i=0; i<NDIM; i++)
                x[i] = x[i]*cell_width[i] + cell[i][0];
        };

    private:
        /// Initialize static data
        void initialize() {
            for (int k = 1; k <= MAXK; k++) commondata[k] = FunctionCommonData<T,NDIM>(k);
            initialized = true;
        };

        /// Assignment is not allowed ... not at all now we have reference member
        //FunctionImpl<T>& operator=(const FunctionImpl<T>& other);
    };
}

#endif
