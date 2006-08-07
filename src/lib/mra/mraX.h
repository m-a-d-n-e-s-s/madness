#ifndef MRA_H
#define MRA_H

/// \file mra.h
/// \brief Header for Function and friends

#include <vector>
#include <string>
#include <functional>

#include <madness_config.h>
#include <tensor/tensor.h>
#include <serialize/archive.h>
//#include <serialize/mpiar.h>
//#include <octtree/octtree.h>
#include <octtree/sendrecv.h>
#include <mra/sepop.h>
#include <misc/communicator.h>

namespace std {
    /// This to make norm work as desired for both complex and real
    static inline double norm(const double& d) {
        return d*d;
    };
}

namespace madness {

    /**

    The distributed data structure and class hierarchy underlying the
    Function class is pieced together as follows.

    A FunctionOctTree combines an OctTree with management of indices.
    Each independent communicator (or data distribution) must have its
    own FunctionOctTree.  Not yet implemented (but coming soon) is the
    ability to move Functions between different FunctionOctTrees using
    MPI inter-communicators.  In the meantime, a program only needs
    to create one FunctionOctTree and store it as the default
    inside the FunctionDefaults structure.  

    The OctTree class template provides the overall oct-tree
    structure, data distribution (layout) and connectivity between
    processes.  The OctTree is also associated with a communicator.

    Each node of the OctTree contains a FunctionNode which is little
    more than a vaguely type-safe array of pointers to BaseTensor.

    Each new function is allocated a unique index that is used to
    access its coefficients stored in each FunctionNode.  As functions
    are dynamically created and destroyed, the indices are managed by
    the FunctionOctTree that is uniquely associated with the OctTree.

    Dynamic load balancing according to a variety of criteria is an 
    essential component that is yet to be implemented.

    */


    class FunctionNode;         ///< Forward definition
    typedef OctTree<FunctionNode> OctTreeT; ///< Type of OctTree used to hold coeffs

    /// Used to hold data for all functions at each node of the OctTree

    /// Mildly unsatisfactory feature is we are messing with pointers
    /// to the tensor base class which was envisoned primarily to
    /// support the interface to Python.  Currently, this is statically
    /// sized, but dynamic resizing is straightforward.
    class FunctionNode {
    public:
        static const int size = 64; ///< Size of new vectors.

    private:
        std::vector<BaseTensor *> v; ///< Pointers to base tensor, NULL if not set
        std::vector<bool> a;    ///< Flags for active

        /// Private. Copy constructor not supported
        FunctionNode(const FunctionNode& f);

        /// Private. Assignment not supported
//        FunctionNode& operator=(const FunctionNode& f);

    public:
        /// Constructor initializes all data pointers to NULL
        FunctionNode() : v(size), a(size) {
            // The explicit initialization may not be necessary, but
            // better safe than sorry.
            for (int i = 0; i < size; i++) {
                v[i] = 0;
                a[i] = false;
            };
        };

	FunctionNode(std::vector<BaseTensor *> va, std::vector<bool> aa) :
	v(va), a(aa) {};

	FunctionNode(std::vector<bool> aa) : v(size), a(aa) {};

        // I know, I know, naughty for doing this, but yanno, it needs to be done
/*
        FunctionNode(const FunctionNode& f) {
            int i;
            int n = f.v.size();
            for (i = 0; i < n; i++) {
                BaseTensor* t = f.v[i];
                if ((t)&&(t->id == TensorTypeData<double>::id)) {
                    Tensor<double> *d = new Tensor<double>();
                    *d = *(const Tensor<double> *) t;
		    v[i] = d;
                } else if (t) {
                    throw "not yet";
                }
		a[i] = f.a[i];
            };
        };
*/

        FunctionNode& operator=(const FunctionNode& fn) {
            if (this == &fn) return *this;
            int i;
            int n = fn.v.size();
            for (i = 0; i < n; i++) { 
                BaseTensor* t = fn.v[i];
                if ((t) && (t->id == TensorTypeData<double>::id)) {
                    Tensor<double> *d = new Tensor<double>();
                    *d = *(const Tensor<double> *) t;
		    v[i] = d;
                } else if (t) {
                    throw "not yet";
                }
		a[i] = fn.a[i];
            };
            return *this;
        };

        /// Puts vector of tensors and vector of bools into supplied arguments

        template <typename T>
        inline void getTensorList(std::vector<Tensor<T> > *tv, std::vector<bool> *av) {
            int i = 0, n = size;
//std::cout << "beginning of getTensorList, n = " << n << std::endl;
            for (i = 0; i < n; i++) {
//std::cout << "loop of getTensorList, i = " << i << std::endl;
                tv->push_back(*(get<T>(i)));
                av->push_back(isactive(i));
            }
//std::cout << "end of getTensorList, n = " << n << std::endl;
        }


	inline int getTensorListSize() {
	    return v.size();
	}

        /// Sets v and a using supplied arguments
/*
        template <typename T>
        inline void setTensorList(int n, std::vector<Tensor<T> > tv, std::vector<bool> av) {
            int i = 0;
            for (i = 0; i < n; i++) {
                set(i, tv[i]);
                set_active_status(i, av[i]);
            }
        }
*/

	/// Sets a to be equal to input vector
	inline void setActiveList(std::vector<bool> av)
	{
	    int avlen = av.size();
	    int n;
	    if (avlen < size) 
		n = avlen;
	    else
		n = size;
	
//	    std::cout << "setActiveList: n = " << n << std::endl;
	    for (int i = 0; i < n; i++)
	    {
		a[i] = av[i];
	    }
	    for (int i = n; i < size; i++)
	    {
		a[i] = false;
	    }
//	    std::cout << "setActiveList: all done" << std::endl;
//	    std::cout << "setActiveList: v.size() = " << v.size() << std::endl;
	}

	/// Returns a
	inline std::vector<bool> getActiveList()
	{
	    return a;
	}

        /// Returns a pointer to data for entry ind.  NULL indicates no data.
        template <typename T>
        inline Tensor<T>* get(int ind) const {
                if (ind < 0 || ind >= (int)v.size())
                    throw "FunctionNode: get: invalid index";
                return (Tensor<T> *) v[ind];
            }

        /// Sets data (and active) by taking a \em shallow copy of the tensor value

        /// Returns pointer to data as if from get()
        template <typename T>
        inline Tensor<T>* set(int ind, const Tensor<T>& value) {
//std::cout << "set: ind = " << ind << ", v.size() = " << v.size() << std::endl;
            if (ind < 0 || ind >= (int)v.size())
                throw "FunctionNode: set: invalid index";
            v[ind] = new Tensor<T>(value);
            set_active(ind);
            return (Tensor<T> *) v[ind];
        }

        /// Unsets the data pointer (if any) and frees any associated data

        /// This does \em not change the active status of the function in
        /// this node.
        inline void unset(int ind) {
            if (ind < 0 || ind >= size)
                throw "FunctionNode: unset: invalid index";
            if (v[ind]) {
                delete v[ind];
                v[ind] = 0;
            }
        };

        /// Return true if function ind is active in this node
        inline bool isactive(int ind) const {
            return a[ind];
        };

        /// Mark function ind as active in this node
        inline void set_active(int ind) {
//std::cout << "set_active: about to set index " << ind << " active" << std::endl;
            a[ind] = true;
//std::cout << "set_active: just set index " << ind << " active" << std::endl;
        };

        /// Set active status of function ind in this node
        inline void set_active_status(int ind, bool status) {
            a[ind] = status;
        };

        /// Mark function ind as inactive in this node
        inline void set_inactive(int ind) {
            a[ind] = false;
        };

        /// Destructor frees all data
        ~FunctionNode() {
            for (int i = 0; i < size; i++) unset(i);
        };


        template <class Archive>
        inline void serialize(const Archive& ar) {
            ar & v & a;
        }
    };

    /// A FunctionOctTree marries an OctTree<FunctionNode> with an index manager

    /// The constructor takes ownership of the provided pointer to an already
    /// constructed OctTree<FunctionNode>.  I.e., don't delete the OctTree.
    class FunctionOctTree {
    private:
        int _nfree;
        std::vector<int> _free;
        std::vector<SharedPtr<OctTreeT> > _treeList;

        /// Verboten
        FunctionOctTree(const FunctionOctTree&);

        /// Verboten
//        FunctionOctTree& operator=(const FunctionOctTree&);

        /// Cleans up data of function being freed
        void free_data(OctTreeT* tree, int ind) {
            tree->data().unset(ind);
            FOREACH_CHILD(OctTreeT, tree, free_data(child, ind););
        };

    public:
        FunctionOctTree(OctTree<FunctionNode>* tree)
                : _nfree(FunctionNode::size)
                , _free(_nfree)
        , _treeList(std::vector< SharedPtr<OctTree<FunctionNode> > >(SharedPtr<OctTree<FunctionNode> > (tree))) {
            for (int i = 0; i < _nfree; i++) _free[i] = _nfree - i - 1;
        };

        FunctionOctTree(std::vector<SharedPtr<OctTree<FunctionNode> > > tree)
                : _nfree(FunctionNode::size)
                , _free(_nfree)
        , _treeList(tree) {
            for (int i = 0; i < _nfree; i++) _free[i] = _nfree - i - 1;
        };

        FunctionOctTree(std::vector<SharedPtr<OctTree<FunctionNode> > > tree, int nfun)
                : _nfree(FunctionNode::size)
                , _free(_nfree)
        , _treeList(tree) {
            for (int i = 0; i < _nfree; i++) _free[i] = _nfree - i - 1;
	    alloc(nfun);
        };

	/// change treeList without messing with anything else
	inline void setTreeList(std::vector<SharedPtr<OctTree<FunctionNode> > > treeList)
	{
	    _treeList.clear();
	    int tlen = treeList.size();
	    for (int i = 0; i < tlen; i++)
	    {
		_treeList.push_back(treeList[i]);
	    }
	}

	/// Depth-first traversal for diagnostic purposes
	inline void depthFirstTraverse()
	{
	    int n = _treeList.size();
	    for (int i = 0; i < n; i++)
	    {
		std::cout << "Subtree " << i << "of " << n << ":" << std::endl;
	    	_treeList[i]->depthFirstTraverseParents();
	    }
	};

	inline int nTrees() {return _treeList.size();};

	inline std::vector<SharedPtr<OctTree<FunctionNode> > > treeList()
	{
	    return _treeList;
	};

        /// Return a pointer to the contained octtree
        inline SharedPtr<OctTree<FunctionNode> > tree(int index) {
            return SharedPtr<OctTree<FunctionNode> > (_treeList[index]);
        };

        /// Return a const pointer to the contained octtree
        inline const SharedPtr<OctTree<FunctionNode> > tree(int index) const {
            return SharedPtr<OctTree<FunctionNode> > (_treeList[index]);
        };

        /// Allocate an index for a new function
        int alloc() {
            _nfree--;
            if (_nfree < 0)
                throw "FunctionIndexManager: too many functions: no resize yet";
//            std::cout << "Allocated Function: " << _free[_nfree] << " " << std::endl;
            return _free[_nfree];
        };

        /// Allocate n indices for functions
        int alloc(int n) {
	    int i = 0;
	    for (i = 0; i < n; i++)
	    {
		alloc();
	    }
            return _free[_nfree];
        };

	int getNalloc() {
	    return _free[_nfree]+1;
	};

        /// Free an index for a now deleted function
        void free(int ind) {
            if (ind < 0 || ind >= (int)_free.size())
                throw "FunctionIndexManager: freeing invalid index";
            for (int i = 0; i < _nfree; i++)
                if (_free[i] == ind)
		{
		    std::cout << "free: " << ind << " is already free, doofus!" << std::endl;
                    throw "FunctionIndexManager: freeing already free index";
		}

	    for (int i = 0; i < _treeList.size(); i++)
            	free_data(_treeList[i], ind);

//            std::cout << "Deallocated Function: " << ind << std::endl;
            _free[_nfree] = ind;
            _nfree++;
        };
    };

    template <typename T> class FunctionFactory; // Forward definition to define factory before Function

    /// FunctionDefaults holds default paramaters as static class members

    /// Defined and initialized in mra.cc.  Note that it is currently
    /// not possible to combine functions with different simulation
    /// cells and this is is not even checked for.  Therefore, the
    /// only recommended approach is to set the simulation cell once
    /// and for all at the start of a calculation and to use it for
    /// all Function instances.  A little work could improve upon this.
    class FunctionDefaults {
    public:
        static SharedPtr<FunctionOctTree> tree; ///< The default tree/distribution used for new functions
        static int k;                  ///< Wavelet order
        static double thresh;          ///< Truncation threshold
        static int initial_level;      ///< Initial level for fine scale projection
        static int max_refine_level;   ///< Level at which to stop refinement
        static int truncate_method;    ///< Truncation method
        static bool compress;          ///< Whether to compress new functions
        static bool refine;            ///< Whether to refine new functions
        static bool autorefine;        ///< Whether to autorefine in multiplication, etc.
        static bool debug;             ///< Controls output of debug info
        static double cell[3][2];      ///< Simulation cell, cell[0][0]=xlo, cell[0][1]=xhi, ...
    };

    /// FunctionCommonData holds all Function data common for given k

    /// Since Function assignment and copy constructors are shallow it
    /// greatly simplifies maintaining consistent state to have all
    /// (permanent) state encapsulated in a single class.  The state
    /// is shared between instances using a SharedPtr.   Also, separating
    /// shared from instance specific state greatly accelerates the constructor.
    /// which is important for massive parallelism.  The default copy constructor and
    /// assignment operator are used.
    template <typename T>
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
        double cell[3][2];      ///< Simulation cell (range over function is defined).
        double cell_width[3];   ///< Size of simulation cell in each dimension
        double cell_volume;     ///< Volume of simulation cell

        Slice s[2];        ///< s[0]=Slice(0,k-1), s[1]=Slice(k,2*k-1)
        std::vector<Slice> s0;  ///< s[0] in each dimension to get scaling coeffs

        typedef Tensor<T> TensorT; ///< Type of tensor used to hold coeffs
        Tensor<T> work1;        ///< work space of size (k,k,k)
        Tensor<T> work2;        ///< work space of size (2k,2k,2k)
        Tensor<T> workq;        ///< work space of size (npt,npt,npt)
        TensorT zero_tensor;    ///< Zero (k,k,k) tensor for internal convenience of diff

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
        FunctionCommonData() {}
        ;

        /// Constructor presently forces choice npt=k
        FunctionCommonData(int k) : k(k), npt(k) {
            s[0] = Slice(0, k - 1);
            s[1] = Slice(k, 2 * k - 1);
            s0 = std::vector<Slice>(3);
            s0[0] = s0[1] = s0[2] = s[0];
            cell_volume = 1.0;
            for (int i = 0; i < 3; i++) {
                cell[i][0] = FunctionDefaults::cell[i][0];
                cell[i][1] = FunctionDefaults::cell[i][1];
                cell_width[i] = cell[i][1] - cell[i][0];
                cell_volume *= cell_width[i];
            }
            work1 = TensorT(k, k, k);
            work2 = TensorT(2 * k, 2 * k, 2 * k);
            zero_tensor = TensorT(k, k, k);

            _init_twoscale();
            _init_quadrature();
            _make_dc_periodic();
        };

        /// Convert user coords (cell[][]) to simulation coords ([0,1])
        inline void user_to_sim(double& x, double& y, double& z) const {
            x = (x - cell[0][0]) / cell_width[0];
            y = (y - cell[1][0]) / cell_width[1];
            z = (z - cell[2][0]) / cell_width[2];
        };

        /// Convert simulation coords ([0,1]) to user coords (cell[][])
        inline void sim_to_user(double& x, double& y, double& z) const {
            x = x * cell_width[0] + cell[0][0];
            x = x * cell_width[1] + cell[1][0];
            x = x * cell_width[2] + cell[2][0];
        };
    };

    /// FunctionData holds all Function state to facilitate shallow copy semantics

    /// Since Function assignment and copy constructors are shallow it
    /// greatly simplifies maintaining consistent state to have all
    /// (permanent) state encapsulated in a single class.  The state
    /// is shared between instances using a SharedPtr.
    template <typename T>
    class FunctionData {
    private:
        static const int MAXK = 17;
        static FunctionCommonData<T> commondata[MAXK + 1]; 	///< Defined in mra.cc
        static bool initialized;	///< Defined and initialized to false in mra.cc

    public:
        typedef Tensor<T> TensorT; ///< Type of tensor used to hold coeffs

        int k;                  ///< Wavelet order
        double thresh;          ///< Screening threshold
        int initial_level;      ///< Initial level for refinement
        int max_refine_level;   ///< Do not refine below this level
        int truncate_method;    ///< 0=default=(|d|<thresh), 1=(|d|<thresh/2^n);
        bool debug;             ///< If true, verbose printing ... unused?
        bool autorefine;        ///< If true, autorefine where appropriate

        FunctionCommonData<T>* cdata;

        T (*f)(double, double, double); ///< Scalar interface to function to compress
        void (*vf)(long, const double*, const double*, const double*,
                   T* restrict); ///< Vector interface to function to compress

        SharedPtr<FunctionOctTree> trees; ///< The subtrees of coeffs
        int ind;		 ///< The unique index into the tree for this function

        bool compressed;        ///< Compression status
        long nterminated;               ///< No. of boxes where adaptive refinement was too deep


        /// Initialize function data from data in factory
        FunctionData(const FunctionFactory<T>& factory) {
            // might be cleaner to assign function factory data
            // via a struct so that don't need to manually modify here
            // every time a new element is added.
            if (!initialized) initialize();
            k = factory._k;
            cdata = commondata + k;

            trees = factory._tree;
            ind = trees->alloc();

            thresh = factory._thresh;
            initial_level = factory._initial_level;
            max_refine_level = factory._max_refine_level;
            truncate_method = factory._truncate_method;
            debug = factory._debug;
            autorefine = factory._autorefine;
            f = factory._f;
            vf = factory._vf;

            compressed = false;
            nterminated = 0;
        };

        /// Copy constructor

        /// Allocates a \em new function index in preparation for a deep copy
        FunctionData(const FunctionData<T>& other) {
            if (!initialized) initialize();
            k = other.k;
            cdata = commondata + k;

            trees = other.trees;
            ind = trees->alloc();

            thresh = other.thresh;
            initial_level = other.initial_level;
            max_refine_level = other.max_refine_level;
            truncate_method = other.truncate_method;
            debug = other.debug;
            f = other.f;
            vf = other.vf;

            compressed = other.compressed;
            nterminated = other.nterminated;
        };

        ~FunctionData() {
            	trees->free(ind);
        };

    private:
        /// Initialize static data
        void initialize() {
            for (int k = 1; k <= MAXK; k++) commondata[k] = FunctionCommonData<T>(k);
            initialized = true;
        };

        /// Assignment is not allowed
        FunctionData<T>& operator=(const FunctionData<T>& other);
    };


    /// Private.  Helps with unary negation operator.
    template <typename T>
    static inline Tensor<T> _negate_helper(const Tensor<T>& t) {
        return -t;
    };

    /// Private.  Helps with scaling by a constant.
    template <typename T>
    class _scale_helper {
    private:
        T s;
    public:
        _scale_helper(T value) : s(value) {}
        ;

        inline Tensor<T> operator()(const Tensor<T>& t) const {
            return t*s;
        };
    };

    /// Private.  Helps with adding/subtracting a constant.
    template <typename T>
    class _add_helper {
    private:
        T s;
    public:
        _add_helper(T value) : s(value) {}
        ;

        inline Tensor<T> operator()(const Tensor<T>& t) const {
            return t + s;
        };
    };


    /// Multiresolution 3d function of given type
    template <typename T>
    class Function {
        typedef Tensor<T> TensorT; ///< Type of tensor used to hold coeffs

        /// Private.  This funky constructor makes a partial deep copy.

        /// It takes a deep copy of data, but does NOT copy the
        /// coefficients.  Used to initialize the result of application
        /// of an operator to the original function.  No communication.
        explicit Function(const Function<T>* t)
                : data(new FunctionData<T>(*(t->data)))
                , k(data->k)
        , ind(data->ind) {}

    public:
        SharedPtr< FunctionData<T> > data; ///< Holds all function data
        int k;                  ///< For convenience replicates k from data
        int ind;                ///< For convenience replicates ind from data

        /// Constructor from FunctionFactory provides named parameter idiom

        /// No communication is involved.  Default constructor makes a
        /// zero compressed function.
        Function(const FunctionFactory<T>& factory = FunctionFactory<T>())
                : data(new FunctionData<T>(factory))
                , k(data->k)
        , ind(data->ind) {
//	    std::cout << "Function: FunctionFactory constructor: about to _init" << std::endl;
            _init(factory);
//	    std::cout << "Function: FunctionFactory constructor: done with _init" << std::endl;
        };


        /// Print out a summary of the tree with norms
        void pnorms() {
	    int tlen = data->trees->nTrees();
	    for (int i = 0; i < tlen; i++)
	    {
//		std::cout << "pnorms: tree number " << i <<  " of " << tlen << ":" << std::endl;
            	if (isactive(tree(i))) 
		    _pnorms(tree(i));
		else
		{
//		    std::cout << "pnorms: tree " << tree(i)->n() << " "  << tree(i)->x() << " "  
//			<< tree(i)->y() << " "  << tree(i)->z() << "is inactive" << std::endl;
		}
//		std::cout << "pnorms: end of tree number " << i << std::endl << std::endl;
	    }
//	    std::cout << "pnorms: end of function" << std::endl << std::endl;
        };

        /// Compress function (scaling function to wavelet)

        /// Communication streams up the tree.
        /// Returns self for chaining.
        Function<T>& compress() {
//	    std::cout << "compress: at beginning" << std::endl;
            if (!data->compressed) {
		int tlen = data->trees->nTrees();
		for (int i = 0; i < tlen; i++)
		{
//	    	    std::cout << "compress: compressing tree " << i << " of " << tlen << std::endl;
//		    std::cout << "   local root: n = " << tree(i)->n() << ", (" << tree(i)->x() << ","
//			<< tree(i)->y() << "," << tree(i)->z() << ")" << std::endl;
                    if (isactive(tree(i))) 
		    {
			_compress(tree(i));
		    }
		    else
		    {
//			std::cout << "compress: this tree is not active" << std::endl;
		    }
		}
                data->compressed = true;
            }
//	    std::cout << "compress: at end" << std::endl;
            return *this;
        };


        /// Reconstruct compressed function (wavelet to scaling function)

        /// Communication streams down the tree
        /// Returns self for chaining.
        Function<T>& reconstruct() {
            if (data->compressed) {
		int tlen = data->trees->nTrees();
		for (int i = 0; i < tlen; i++)
		{
                    if (isactive(tree(i))) 
			_reconstruct(tree(i));
		}
                data->compressed = false;
            }
            return *this;
        };

        /// Return true if the function is compressed
        bool iscompressed() const {
            return data->compressed;
        };


    private:
        /// Private.  Initialization required constructing function from stratch

        /// No communication involved.
        void _init(const FunctionFactory<T>& factory);


        /// Private.  Projects function in given box

        /// No communication involved.
        void _project(OctTreeT* s);


        /// Private.  Projects function at given level.  No communication.

        /// No communication involved.  It is assumed that this is
        /// being done to an initially zero Function by _init
        /// immediately after construction
        void _fine_scale_projection(OctTreeT* coeff, Level initial_level);

        /// Evaluate f(x,y,z) on the quadrature grid in box (n,lx,ly,lz)

        /// No communication involved.
        void _fcube(long n, long lx, long ly, long lz,
                    T (*f)(double, double, double), TensorT& fcube);


        /// Vector interface to evaluating function on the quadrature grid in box (n,lx,ly,lz)

        /// No communication involved.
        void _vfcube(long n, long lx, long ly, long lz,
                     void (*vf)(long l, const double*, const double*, const double*,
                                T* RESTRICT),
                     TensorT& fcube);


        /// Private.  Transform sum coefficients at level n to sums+differences at level n-1

        /// Given scaling function coefficients s[n][l][i] and s[n][l+1][i]
        /// return the scaling function and wavelet coefficients at the
        /// coarser level.  I.e., decompose Vn using Vn = Vn-1 + Wn-1.
        ///
        /// s_i = sum(j) h0_ij*s0_j + h1_ij*s1_j
        /// d_i = sum(j) g0_ij*s0_j + g1_ij*s1_j
        ///
        /// Returns a new tensor and has no side effects.  Works for any
        /// number of dimensions.
        ///
        /// No communication involved.
        inline TensorT filter(const TensorT& s) const {
            return transform3d(s, data->cdata->hgT);
            //return transform(s,hgT);
        };

        /// Private: Optimized filter (inplace, contiguous, no err checking)

        /// Transforms coefficients in s returning result also in s.
        /// Uses work2 from common data to eliminate temporary creation and
        /// to increase cache locality.
        ///
        /// No communication involved.
        inline void filter_inplace(TensorT& s) {
//std::cout << "filter_inplace: at beginning" << std::endl;
            transform3d_inplace(s, data->cdata->hgT, data->cdata->work2);
//std::cout << "filter_inplace: at end" << std::endl;
        };


        /// Private: Optimized unfilter (see info about filter_inplace)

        /// No communication involved.
        inline void unfilter_inplace(TensorT& s) {
            transform3d_inplace(s, data->cdata->hg, data->cdata->work2);
        };



        /// Private.  Transform sums+differences at level n to sum coefficients at level n+1

        ///  Given scaling function and wavelet coefficients (s and d)
        ///  return the scaling function coefficients at the next finer
        ///  level.  I.e., reconstruct Vn using Vn = Vn-1 + Wn-1.
        ///
        ///  s0 = sum(j) h0_ji*s_j + g0_ji*d_j
        ///  s1 = sum(j) h1_ji*s_j + g1_ji*d_j
        ///
        ///  Returns a new tensor and has no side effects
        ///
        ///  If (sonly) ... then ss is only the scaling function coeffs (and
        ///  assume the d are zero).  Works for any number of dimensions.
        ///
        /// No communication involved.
        inline TensorT unfilter(const TensorT& ss,
                                bool sonly = false) const {
            if (sonly)
                //return transform(ss,data->cdata->hgsonly);
                throw "unfilter: sonly : not yet";
            else
                return transform3d(ss, data->cdata->hg);
            //return transform(ss,hg);
        };


        /// Private. Returns the tree of coeffs

        /// No communication involved.
        inline SharedPtr<OctTreeT> tree(int index) {
            return data->trees->tree(index);
        };


        /// Private. Returns the const tree of coeffs

        /// No communication involved.
        inline const SharedPtr<OctTreeT> tree(int index) const {
            return data->trees->tree(index);
        };


        /// Private.  Retrieve pointer to coeffs of this Function from tree pointer

        /// No communication involved.
        inline TensorT* coeff(OctTreeT* tree) const {
            return tree->data(). template get<T>(ind);
        };


        /// Private.  Retrieve const pointer to coeffs of this Function from tree pointer

        /// No communication involved.
        inline const TensorT* coeff(const OctTreeT* tree) const {
            return tree->data(). template get<T>(ind);
        };


        /// Private.  Unset any coeffs of this Function from tree pointer

        /// No communication involved.
        inline void unset_coeff(OctTreeT* tree) {
            return tree->data().unset(ind);
        };


        /// Private.  Set coeffs of this Function from tree pointer

        /// No communication involved.
        inline TensorT* set_coeff(OctTreeT* tree, const TensorT& t) {
            return tree->data().set(ind, t);
        };


        /// Private.  Set Function active in tree node

        /// No communication involved.
        inline void set_active(OctTreeT* tree) const {
            tree->data().set_active(ind);
        };


        /// Private.  Set Function inactive in tree node

        /// No communication involved.
        inline void set_inactive(OctTreeT* tree) const {
            tree->data().set_inactive(ind);
        };


        /// Private.  Determine if Function is active in tree node

        /// No communication involved.
        inline bool isactive(const OctTreeT* tree) const {
            return tree->data().isactive(ind);
        };


        /// Private.  For semantic similarity to isactive.

        /// No communication involved.
        inline bool islocal(const OctTreeT* tree) const {
            return tree->islocal();
        };


        /// Private.  For semantic similarity to isactive.

        /// No communication involved.
        inline bool isremote(const OctTreeT* tree) const {
            return tree->isremote();
        };


        /// Private.  Retrieve a pointer to the communicator

        /// No communication involved.
        inline Communicator* comm() const {
            return (Communicator *) (tree(0)->comm());
        };


        /// Private.  Recursive function to refine from initial projection.

        /// Communication streams all the way down the tree.
        void _refine(OctTreeT* tree);


        /// Private.  Recursive function to compress tree.

        /// Communication streams up the tree.
        void _compress(OctTreeT* tree);


        /// Private.  Recursive function to reconstruct the tree.

        /// Communication streams down the tree.
        void _reconstruct(OctTreeT* tree);


        /// Private.  Recur down the tree printing out the norm of
        /// the coefficients.
        void _pnorms(OctTreeT *tree) const {
//	    std::cout << "_pnorms: tree " << tree->n() << " " << tree->x() << " " << tree->y() << " "
//			<< tree->z() << std::endl;
            const TensorT *t = coeff(tree);
            if (t) {
//		std::cout << "_pnorms: this tensor exists" << std::endl;
                for (long i=0; i<tree->n(); i++) std::printf("  ");
                std::printf("%4ld (%8ld, %8ld, %8ld) %9.1e\n",
                            tree->n(),tree->x(),tree->y(),tree->z(),t->normf());
		std::fflush(NULL);
            }
	    else
	    {
//		std::cout << "_pnorms: this tensor does not exist" << std::endl;
	    }
            FOREACH_CHILD(OctTreeT, tree, if (isactive(child)) _pnorms(child););
//	    std::cout << "_pnorms: very end of function" << std::endl;
        }

        /// Private.  Applies truncation method to give truncation threshold

        /// No communication involved.
        double truncate_tol(double tol, Level n) {
            if (data->truncate_method == 0) {
                return tol;
            } else if (data->truncate_method == 1) {
                return tol / two_to_power(n);
            } else {
                throw "truncation_tol: unknown truncation method?";
            }
        };

        /// Private.  Recursive function to support inplace scaling.
        void _scale(OctTreeT *tree, T s) {
            if (coeff(tree)) coeff(tree)->scale(s);
            FOREACH_CHILD(OctTreeT, tree,
                          if (isactive(child)) _scale(child, s););
        };

    };

    /// FunctionFactory implements the named-parameter idiom for Function

    /// C++ does not provide named arguments (as does, e.g., Python).
    /// This class provides something very close.  Create functions as follows
    /// \code
    /// double myfunc(double x,double y,double z);
    /// Function f = FunctionFactory(myfunc).k(11).thresh(1e-9).debug().nocompress()
    /// \endcode
    /// where the methods of function factory, which specify the non-default
    /// arguments eventually passed to the \c Function constructor, can be
    /// used in any order.
    template <typename T> class FunctionFactory {
    protected:
        T (*_f)(double, double, double);
        void (*_vf)(long, const double*, const double*, const double*, T* restrict);
        int _k;

        double _thresh;
        int _initial_level;
        int _max_refine_level;
        int _truncate_method;
        bool _compress;
        bool _refine;
        bool _debug;
        bool _empty;
        bool _autorefine;
        SharedPtr<FunctionOctTree> _tree;
        friend class Function<T>;
        friend class FunctionData<T>;
    public:
        FunctionFactory(T (*f)(double, double, double) = 0)
                : _f(f)
                , _vf(0)
                , _k(FunctionDefaults::k)
                , _thresh(FunctionDefaults::thresh)
                , _initial_level(FunctionDefaults::initial_level)
                , _max_refine_level(FunctionDefaults::max_refine_level)
                , _truncate_method(FunctionDefaults::truncate_method)
                , _compress(FunctionDefaults::compress)
                , _refine(FunctionDefaults::refine)
                , _debug(FunctionDefaults::debug)
                , _empty(false)
                , _autorefine(FunctionDefaults::autorefine)
        , _tree(FunctionDefaults::tree) {}
/*
        , _tree(FunctionDefaults::tree) 
		{
		    std::cout << "FunctionFactory: basic constructor" << std::endl;
		    std::cout << "Length of treeList = " << _tree.size() << std::endl;
		}
*/
        ;
        inline FunctionFactory& vf(void (*vf)(long, const double*, const double*, const double*, T* restrict)) {
            _vf = vf;
            return *this;
        };
        inline FunctionFactory& k(int k) {
            _k = k;
            return *this;
        };
        inline FunctionFactory& thresh(double thresh) {
            _thresh = thresh;
//	    std::cout << "thresh: inside, about to return" << std::endl;
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
//	    std::cout << "compress: inside, about to return" << std::endl;
            return *this;
        };
        inline FunctionFactory& nocompress(bool nocompress = true) {
            _compress = !nocompress;
            return *this;
        };
        inline FunctionFactory& debug() {
            _debug = true;
            return *this;
        };
        inline FunctionFactory& nodebug() {
            _debug = false;
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
        inline FunctionFactory& tree(const SharedPtr<FunctionOctTree>& tree) {
            _tree = tree;
        };
    };

    void balanceFunctionOctTree(SharedPtr<FunctionOctTree> fnList);

    void setRemoteActive(std::vector< SharedPtr <OctTree <FunctionNode > > > *treeList);

    class ActiveRootList
    {
	public:
	RootList r;
	std::vector<bool> activeList;

	ActiveRootList(): r(RootList()), activeList(std::vector<bool>()) {};

	ActiveRootList(RootList rl, std::vector<bool> al)
	{
	    r = rl;
	    activeList = al;
	}

	template <typename T>
	ActiveRootList(OctTree<T> *t, std::vector<bool> al):
	    r(RootList(t, t->rank(), t->rank())), activeList(al) {};
	
	template <typename T>
	ActiveRootList(OctTree<T> *t, ProcessID p, std::vector<bool> al):
	    r(RootList(t, p, t->rank())), activeList(al) {};

        template <class Archive>
        inline void serialize(const Archive& ar) {
            ar & r & activeList;
        }

	inline bool equals(RootList rl)
	{
	    return (r.equals(rl));
	}

	inline friend bool operator < (const ActiveRootList &a1, const ActiveRootList &a2)
	{
	    return (a1.r < a2.r);
	}

    };


    struct sless {
        bool operator() (const ActiveRootList p1, const ActiveRootList p2) {
            if (p1.r.current_owner < p2.r.current_owner)
                return true;
            else
                return false;
        }
    };

    struct rless {
	bool operator() (const RootList p1, const RootList p2) {
	    if (p1.current_owner < p2.current_owner)
		return true;
	    else
		return false;
	}
    };


    inline void sortBySend(std::vector<ActiveRootList> *list)
    {
	sort(list->begin(), list->end(), sless());
    }

    inline void sortBySend(std::vector<RootList> *list)
    {
	sort(list->begin(), list->end(), rless());
    }



    void findRemoteChildrenList(std::vector<RootList> *childList, SharedPtr<OctTree<FunctionNode > > tree);

    void makeRemoteParentList(std::vector<ActiveRootList> *parentList, SharedPtr<OctTree<FunctionNode> > 
		tree);
}

#endif
