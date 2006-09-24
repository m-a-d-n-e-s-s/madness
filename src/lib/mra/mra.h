#ifndef MRA_H
#define MRA_H

/// \file mra.h
/// \brief Header for Function and friends

#include <vector>
#include <string>
#include <functional>
#include <cmath>

#include <madness_config.h>
#include <tensor/tensor.h>
#include <serialize/archive.h>
#include <serialize/mpiar.h>
#include <octtree/octtree.h>
//#include <octtree/sendrecv.h>
#include <mra/sepop.h>
#include <misc/communicator.h>
#include <misc/madexcept.h>
#include <tasks/sav.h>
#include <tasks/tasks.h>

#include <serialize/textfsar.h>
using madness::archive::TextFstreamInputArchive;
using madness::archive::TextFstreamOutputArchive;

#include <serialize/binfsar.h>
using madness::archive::BinaryFstreamInputArchive;
using madness::archive::BinaryFstreamOutputArchive;

#include <serialize/vecar.h>
using madness::archive::VectorInputArchive;
using madness::archive::VectorOutputArchive;

namespace std {
    // This to make norm work as desired for both complex and real
    static inline double norm(const double d) {
        return d*d;
    }
}

namespace madness {

    Communicator& startup(int argc, char** argv);

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

#define FOREACH_ACTIVE_CHILD(childptrT, t, expr) \
    do { \
      for (int i=0; i<2; i++) \
        for (int j=0; j<2; j++) \
	  for (int k=0; k<2; k++) { \
            childptrT child = ((t)->child(i,j,k)); \
            if (child && isactive(child)) {expr} \
          } \
    } while(0)


    template <typename T> class Function;
    class FunctionNode;         ///< Forward definition
    typedef OctTree<FunctionNode> OctTreeT; ///< Type of OctTree used to hold coeffs
    typedef SharedPtr<OctTreeT> OctTreeTPtr; ///< SharedPtr to OctTree

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
        std::vector<bool> acflag; ///< Flags for autocleaning temporary values from tree

   public:
        /// Constructor initializes all data pointers to NULL
        FunctionNode() : v(size), a(size), acflag(size) {
            // The explicit initialization may not be necessary, but better safe than sorry.
            for (int i = 0; i < size; i++) {
                v[i] = 0;
                a[i] = false;
                acflag[i] = false;
            };
        };

        /// Copy constructor needed for data migration ???
        FunctionNode(const FunctionNode& f) {
            int i;
            int n = f.v.size();
            for (i = 0; i < n; i++) {
                BaseTensor* t = f.v[i];
                if (t->id == TensorTypeData<double>::id) {
                    Tensor<double> *d = new Tensor<double>();
                    *d = *(const Tensor<double> *) t;
                    v.push_back(d);
                    a.push_back(f.a[i]);
                } else {
                    MADNESS_EXCEPTION("not yet",0);
                }
            };
        };

        /// Asignment ... seems to be buggy?
        FunctionNode& operator=(const FunctionNode& fn) {
            if (this == &fn) return *this;
            int i;
            int n = fn.v.size();
            v.clear();
            a.clear();
            MADNESS_EXCEPTION("RJH ... THIS SEEMS WRONG.  IS IT ACTUALLY BEING USED?????",0);
            for (i = 0; i < n; i++) {
                BaseTensor* t = fn.v[i];
                if (t->id == TensorTypeData<double>::id) {
                    Tensor<double> *d = new Tensor<double>();
                    *d = *(const Tensor<double> *) t;
                    v.push_back(d);
                    a.push_back(fn.a[i]);
                } else {
                    MADNESS_EXCEPTION("not yet",1);
                }
            };
            return *this;
        };

        /// Puts vector of tensors and vector of bools into supplied arguments
        template <typename T>
        inline void getTensorList(std::vector<Tensor<T> > *tv, std::vector<bool> *av) {
            int i = 0, n = size;
            for (i = 0; i < n; i++) {
                tv->push_back(*(get<T>(i)));
                av->push_back(isactive(i));
            }
        }

        /// Sets v and a using supplied arguments
        template <typename T>
        inline void setTensorList(int n, std::vector<Tensor<T> > tv, std::vector<bool> av) {
            int i = 0;
            for (i = 0; i < n; i++) {
                set(i, tv[i]);
                set_active_status(i, av[i]);
            }
        }

        /// Returns a pointer to data for entry ind.  NULL indicates no data.
        template <typename T>
        inline Tensor<T>* get(int ind) const {
            if (ind < 0 || ind >= (int)v.size())
                MADNESS_EXCEPTION("FunctionNode: get: invalid index",ind);
            return (Tensor<T> *) v[ind];
        }

        /// Sets data (and active) by taking a \em shallow copy of the tensor value

        /// Returns pointer to data as if from get()
        template <typename T>
        inline Tensor<T>* set(int ind, const Tensor<T>& value) {
            if (ind < 0 || ind >= (int)v.size())
                MADNESS_EXCEPTION("FunctionNode: set: invalid index",ind);
            if (v[ind]) delete v[ind];
            v[ind] = new Tensor<T>(value);
            set_active(ind);
            return (Tensor<T> *) v[ind];
        }

        /// Unsets the data pointer (if any) and frees any associated data

        /// This does \em not change the active status of the function in
        /// this node.
        inline void unset(int ind) {
            if (ind < 0 || ind >= size)
                MADNESS_EXCEPTION("FunctionNode: unset: invalid index",ind);
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
            a[ind] = true;
        };

        /// Set active status of function ind in this node
        inline void set_active_status(int ind, bool status) {
            a[ind] = status;
        };

        /// Mark function ind as inactive in this node
        inline void set_inactive(int ind) {
            a[ind] = false;
        };
        
        /// Set auto clean flag to value
        inline void set_acflag(int ind, bool value) {
            acflag[ind] = value;
        };
        
        inline bool get_acflag(int ind) const {
            return acflag[ind];
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

    /// A FunctionOctTree marries an OctTree<FunctionNode> with an index 

    /// The constructor takes ownership of the provided pointer to an already
    /// constructed OctTree<FunctionNode>.  I.e., please don't delete the OctTree.
    class FunctionOctTree {
    private:
        int _nfree;
        std::vector<int> _free;
        SharedPtr<OctTreeT> _tree;

        /// Verboten
        FunctionOctTree(const FunctionOctTree&);

        /// Verboten
        FunctionOctTree& operator=(const FunctionOctTree&);

        /// Cleans up data of function being freed
        void free_data(OctTreeTPtr& tree, int ind) {
            tree->data().unset(ind);
            tree->data().set_inactive(ind);
            tree->data().set_acflag(ind,false);
            FOREACH_CHILD(OctTreeTPtr, tree, free_data(child, ind););
        };

    public:
        FunctionOctTree(OctTreeT* tree)
                : _nfree(FunctionNode::size)
                , _free(_nfree)
                , _tree(tree) {
            for (int i = 0; i < _nfree; i++) _free[i] = _nfree - i - 1;
        };

        /// Return a pointer to the contained octtree
        inline OctTreeTPtr& tree() {
            return _tree;
        };

        /// Return a const pointer to the contained octtree
        inline const OctTreeTPtr& tree() const {
            return _tree;
        };

        /// Allocate an index for a new function and clean out junk
        int alloc() {
            _nfree--;
            if (_nfree < 0) 
                MADNESS_EXCEPTION("FunctionIndexManager: too many functions: no resize yet",FunctionNode::size);
            int ind = _free[_nfree];
            free_data(_tree, ind);  // Not necessary but safer?
            return ind;
        };

        /// Free an index for a now deleted function
        void free(int ind) {
            if (ind < 0 || ind >= (int)_free.size())
                MADNESS_EXCEPTION("FunctionIndexManager: freeing invalid index",ind);
            for (int i = 0; i < _nfree; i++)
                if (_free[i] == ind)
                    MADNESS_EXCEPTION("FunctionIndexManager: freeing already free index",ind);

            free_data(_tree, ind);
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

        Slice s[4];              ///< s[0]=Slice(0,k-1), s[1]=Slice(k,2*k-1), etc.
        std::vector<Slice> s0;  ///< s[0] in each dimension to get scaling coeffs
        std::vector<long> vk;   ///< (k,k,k) used to initialize Tensors
        std::vector<long> v2k;  ///< (2k,2k,2k) used to initialize Tensors

        typedef Tensor<T> TensorT; ///< Type of tensor used to hold coeffs
        Tensor<T> work1;        ///< work space of size (k,k,k)
        Tensor<T> work2;        ///< work space of size (2k,2k,2k)
        Tensor<T> workq;        ///< work space of size (npt,npt,npt)
        Tensor<T> zero_tensor1;    ///< Zero (k,k,k) tensor for internal convenience of diff

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
            for (int i=0; i<4; i++) s[i] = Slice(i*k,(i+1)*k-1);
            s0 = std::vector<Slice>(3);
            s0[0] = s0[1] = s0[2] = s[0];
            cell_volume = 1.0;
            for (int i = 0; i < 3; i++) {
                cell[i][0] = FunctionDefaults::cell[i][0];
                cell[i][1] = FunctionDefaults::cell[i][1];
                cell_width[i] = cell[i][1] - cell[i][0];
                cell_volume *= cell_width[i];
            }
            vk = vector_factory((long)k,(long)k,(long)k);
            v2k = vector_factory(2L*k,2L*k,2L*k);
            work1 = TensorT(vk);
            work2 = TensorT(v2k);
            zero_tensor1 = TensorT(vk);

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
    
    template <typename T> class FunctionData;

    class FunctionDataPointersBase {
    protected:
        friend Communicator& startup(int argc, char** argv);
        static void* p[FunctionNode::size];  // Declared and initialized to zero in startup.cc
    };    

    
    template <typename T>
    class FunctionDataPointers : public FunctionDataPointersBase {
    public:
        static FunctionData<T>* get(int i) {
            if (i<0 || i>=FunctionNode::size)
                MADNESS_EXCEPTION("FunctionDataPointers: get: out of range",i);
            if (!p[i]) 
                MADNESS_EXCEPTION("FunctionDataPointers: get: null pointer",i);
            return (FunctionData<T>*) p[i];
        };
        
        static void set(int i, FunctionData<T>* ptr) {
            if (i<0 || i>=FunctionNode::size)
                MADNESS_EXCEPTION("FunctionDataPointers: set: out of range",i);
            if (p[i]) 
                MADNESS_EXCEPTION("FunctionDataPointers: set: already in use",i);
            p[i] = (void *) ptr;
        };
        
        static void clear(int i) {
            if (i<0 || i>=FunctionNode::size)
                MADNESS_EXCEPTION("FunctionDataPointers: clear: out of range",i);
            p[i] = 0;
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
        bool refine;            ///< If true, refine when constructed
        GlobalTree<3> gtree;    ///< Global root knowledge ... needs relocating

        FunctionCommonData<T>* cdata;

        T (*f)(double, double, double); ///< Scalar interface to function to compress
        void (*vf)(long, const double*, const double*, const double*,
                   T* restrict); ///< Vector interface to function to compress

        SharedPtr<FunctionOctTree> tree; ///< The tree of coeffs
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

            tree = factory._tree;
            ind = tree->alloc();
            FunctionDataPointers<T>::set(ind,this);

            thresh = factory._thresh;
            initial_level = factory._initial_level;
            max_refine_level = factory._max_refine_level;
            truncate_method = factory._truncate_method;
            debug = factory._debug;
            autorefine = factory._autorefine;
            refine = factory._refine;
            f = factory._f;
            vf = factory._vf;

            compressed = false;
            nterminated = 0;

            taskq.global_fence();
        };

        /// Copy constructor

        /// Allocates a \em new function index in preparation for a deep copy
        FunctionData(const FunctionData<T>& other) {
            if (!initialized) initialize();
            k = other.k;
            cdata = commondata + k;

            tree = other.tree;
            ind = tree->alloc();
            FunctionDataPointers<T>::set(ind,this);

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
            tree->free(ind);
            FunctionDataPointers<T>::clear(ind);
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
    }
    
    

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
    
    template <typename T> class TaskAwaitCoeff;
    template <typename T> class TaskDiff;
    template <typename T> class TaskSquare;
    template <typename T> class TaskMult;
    template <typename T> class TaskAutorefine;
    template <typename T, typename Derived> class TaskLeaf;
    template <typename T> class TaskRecurDownToMakeLocal;
    template <typename T> class TaskProjectRefine;
    void mratask_register();
    
    /// Multiresolution 3d function of given type
    template <typename T>
    class Function {
    private:
        typedef Tensor<T> TensorT; ///< Type of tensor used to hold coeffs
        typedef SAV<TensorT> ArgT; ///< Type of single assignment variable for tasks

        friend class TaskDiff<T>;
        friend class TaskLeaf< T, TaskDiff<T> >;
        friend class TaskAutorefine<T>;
        friend class TaskLeaf< T, TaskAutorefine<T> >;
        friend class TaskSquare<T>;
        friend class TaskMult<T>;
        friend class TaskLeaf< T, TaskSquare<T> >;
        friend class TaskLeaf< T, TaskMult<T> >;
        friend class TaskAwaitCoeff<T>;
        friend class TaskRecurDownToMakeLocal<T>;
        friend class TaskProjectRefine<T>;
        friend Communicator& startup(int argc, char** argv);
        friend void mratask_register();
        
        static void set_active_handler(Communicator& comm, ProcessID src, const AMArg& arg); 
        static void recur_down_handler(Communicator& comm, ProcessID src, VectorInputArchive& ar);
        SAV< Tensor<T> > _get_scaling_coeffs2(OctTreeTPtr& t, int axis, int inc);
        void _sock_it_to_me(OctTreeTPtr& tree, Level n, const Translation l[3], SAV< Tensor<T> >& arg);
        static void _sock_it_to_me_handler(Communicator& comm, ProcessID src, const AMArg& arg);
        void _sock_it_to_me_forward_request(Level n, const Translation l[3], SAV< Tensor<T> >& arg, ProcessID dest);
        void _diff(Function<T>& df, OctTreeTPtr& tree, int axis); 
        void _dodiff_kernel(Function<T>& df, OctTreeTPtr& tree, int axis,
                            const TensorT& t0, const TensorT& tm, const TensorT& tp);
        void _recur_coeff_down(OctTreeT* tree, bool keep);
        void _recur_coeff_down2(OctTreeT* tree, bool keep);
        const Tensor<T>* _get_scaling_coeffs(OctTreeTPtr& t, int axis, int inc);
        static void recur_down_to_make_handler(Communicator& comm, ProcessID src, const AMArg& arg);
        void recur_down_to_make_forward_request(Level n, const Translation l[3], ProcessID dest);
        void recur_down_to_make(OctTreeT* p, Level n, const Translation l[3]);
        void fill_in_local_tree(OctTreeT* p, Level n, const Translation l[3]);
        
        void project_refine();
        void _doreconstruct(OctTreeTPtr tree, const Tensor<T>& ss);
        static void reconstruction_handler(Communicator& comm, ProcessID src, VectorInputArchive& ar);
        bool _doproject_refine(OctTreeTPtr& tree);
        
        /// Private.  This funky constructor makes a partial deep copy.

        /// It takes a deep copy of data, but does NOT copy the
        /// coefficients.  Used to initialize the result of application
        /// of an operator to the original function.  No communication.
        explicit Function(const Function<T>* t)
                : data(new FunctionData<T>(*(t->data)))
                , k(data->k)
                , ind(data->ind) {};
                
        /// Private.  This funky constructor is used for remote method invocation        

        /// It generates a shallow copy of a function from the function index
        explicit Function(int ind)
                : data(FunctionDataPointers<T>::get(ind),false,false)
                , k(data->k)
                , ind(ind) {
            MADNESS_ASSERT(data->ind == ind);
        };

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
            _init(factory);
        };


        /// Copy constructor is \em shallow

        /// No communication is involved, works in either basis.
        Function(const Function<T>& f)
                : data(f.data)
                , k(f.k)
                , ind(f.ind) {}
        ;

        /// Deep copy.

        /// No communication involved, works in either basis
        Function<T> copy() const {
            Function<T> result(this); // Partial deep copy

            typedef Tensor<T> (*theop)(const Tensor<T>& t);
            if (isactive(tree())) _unaryop<theop>(result, result.tree(), madness::copy<T>);
            return result;
        };

        /// Assignment is \em shallow.

        /// No communication is involved, works in either basis.
        Function<T>& operator=(const Function<T>& f) {
            if (this != &f) {
                data = f.data;
                k = f.k;
                ind = f.ind;
            }
            return *this;
        };


        /// Print out a summary of the tree with norms
        
        /// No communication is involved, works in either basis.
        void pnorms() {
            if (isactive(tree())) _pnorms(tree());
            std::cout.flush();
        };


        /// Unary negation.  Produces new function.

        /// No communication involved.  Works in either basis.
        Function<T> operator-() const {
            Function<T> result = FunctionFactory<T>().k(k).compress(iscompressed()); //.empty();
            if (isactive(tree()))
                _unaryop(result, result.tree(), _negate_helper<T>);
            return result;
        };

        /// Function times a scalar.  Produces new function.
        
        /// No communication involved.  Works in either basis.
        Function<T> operator*(T s) const {
            Function<T> result = FunctionFactory<T>().k(k).compress(iscompressed()).empty();
            if (isactive(tree()))
                _unaryop(result, result.tree(), _scale_helper<T>(s));
            return result;
        };
        

        /// Scalar * function (see Function<T>::operator*(T t))
        friend inline Function<T> operator*(T t, const Function<T>& f) {
            return f*t;
        };
        

        /// Multiplication of two functions.
        Function<T> operator*(Function<T>& other) {
            return mult(other);
        };
        
        Function<T> mult(Function<T>& other);

        /// Binary addition.  Works in wavelet basis generating new Function

        /// No communication is involved \em unless one of the
        /// functions is not already compressed, in which case
        /// compression is forced.
        Function<T> operator+(const Function<T>& other) const {
            // Constness of this and other is only logical
            if (!iscompressed()) const_cast<Function<T>*>(this)->compress();
            if (!other.iscompressed()) const_cast<Function<T>*>(&other)->compress();

            Function<T> result = FunctionFactory<T>().k(k).compress(true).empty();

            if (isactive(tree()) ||
                    other.isactive(other.tree()))
                _add_sub(result, *this, other, result.tree(), false);

            return result;
        };

        /// Binary addition.  Works in wavelet basis generating new Function

        /// No communication is involved \em unless one of the
        /// functions is not already compressed, in which case
        /// compression is forced.
        Function<T> operator-(const Function<T>& other) const {
            // Constness of this and other is only logical
            if (!iscompressed()) const_cast<Function<T>*>(this)->compress();
            if (!other.iscompressed()) const_cast<Function<T>*>(&other)->compress();

            Function<T> result = FunctionFactory<T>().k(k).compress(true).empty();

            if (isactive(tree()) ||
                    other.isactive(other.tree()))
                _add_sub(result, *this, other, result.tree(), true);

            return result;
        };


        /// Returns the square of the norm of the \em local subtree

        /// No communication involved.  Works in either the scaling
        /// function or wavelet basis.
        double norm2sq_local() const {
            return _norm2sq_local(tree());
        };

        /// Returns the square of the norm of the entire tree

        /// Global communication is implied.   Works in either the scaling
        /// function or wavelet basis.
        double norm2sq() {
            return comm()->global_sum(norm2sq_local());
        };
        
        /// Differentiation
        Function<T> diff(int axis); 
        
        /// Autorefine (e.g., for more accurate application of integral operators)
        
        /// Communication, if any, flows down the tree
        /// Returns self for chaining.
        Function<T>& autorefine();

        
        void ptree() {_ptree(tree());};

        Function<T>& compress();
        void _compress(OctTreeTPtr& tree, ArgT& parent);
        void _compressop(OctTreeTPtr& tree, ArgT args[2][2][2], ArgT& parent);
        ArgT input_arg(const OctTreeTPtr& consumer, const OctTreeTPtr& producer);

        Function<T>& truncate(double tol=0.0) {return *this;};
        /// Reconstruct compressed function (wavelet to scaling function)

        /// Communication streams down the tree
        /// Returns self for chaining.
        Function<T>& reconstruct() {
            if (data->compressed) {
                data->compressed = false;
                if (tree()->n() == 0) _doreconstruct(tree(),Tensor<T>());
                taskq.global_fence();
            }
            return *this;
        };

        /// Inplace scale by a constant

        /// No communication.
        /// Returns self for chaining.
        Function<T>& scale(T s) {
            if (isactive(tree())) _scale(tree(), s);
            return *this;
        };

        /// Inplace square (pointwise multiplication by self)

        /// If not reconstructed, reconstruct() is called which communicates.
        /// If not autorefining, no communication is involved in
        /// the squaring.  If autorefinining, some downward communication
        /// may be involved depending on the distribution of tree nodes.
        ///
        /// Returns self for chaining.
        Function<T>& square();


        /// Inplace Generalized SAXPY to scale and add two functions.

        /// Works with both functions in either the scaling function or
        /// wavelet basis.  If only one function is already compressed, the
        /// other will also be compressed.
        /// \code
        ///   a = a*alpha + b*beta
        /// \endcode
        /// where \c a() and \c b() are Functions and \c alpha and \c beta are constants
        ///
        /// !!!!! If the operation is performed in the scaling
        /// function basis, you must eventually compress the result in
        /// order to correctly sum the scaling function coefficients
        /// at multiple levels in the tree.  This is useful if you
        /// have to produce and accumulate multiple intermediates in the
        /// scaling function basis, so you only have to run the
        /// compression on the final result rather than all
        /// intermediates.
        ///
        /// If no compression is required, no communication is involved.
        /// Returns self for chaining.
        /// The Functions \em must share the same tree (data distribution).
        /// Copying between distributions can only be done with the
        /// redistribute method (yet to be implemented when needed).
        Function<T>& gaxpy(double alpha, const Function<T>& b, double beta) {
            //check_trees(b);
            if (!(this->iscompressed() ^ b.iscompressed())) {
                this->compress();
                const_cast<Function<T>*>(&b)->compress(); // constness of b is logical
            }
            if (isactive(tree()) || b.isactive(b.tree()))
                _gaxpy(*this, alpha, b, beta, tree());
            return *this;
        };


        /// Return true if the function is compressed
        bool iscompressed() const {
            return data->compressed;
        };

//        /// This class stores informations on each branch.
//        class localTreeMember {
//        public:
//            Translation x, y, z;
//            Level n;
//            ProcessID rank;
//            bool remote, active, have_child;
//            template <class Archive>
//            inline void serialize(const Archive& ar) {
//                ar & x & y & z & n & rank & remote & active & have_child;
//            }
//        };
//
//	/// Saving Function members into the file. This function is prepared for sequential job.
//	template <class Archive>
//	void save_local(Archive& ar);
//
//	/// Saving Function members into the file. This member was prepared for save_local.
//	template <class Archive>
//	void _save_local(Archive& ar, const OctTreeT* tree); 
//
//	/// save member is the member to serialize coefficients.
//        ///  This member is already parallelized. 
//	void save(const char* f, const long partLevel, Communicator& comm);
//
//	/// saveMain member controlls before and after treatments of saveManager. 
//	template <class Archive>
//	void saveMain(const char* f, Archive& ar, const long partLevel, Communicator& comm);
//
//	/// saveManager Function member serialize coefficients on local (rank 0 ) branch. 
//	template <class Archive>
//	void saveManager(const char* f, Archive& ar, const OctTreeT* tree,
//              const long partLevel, Communicator& commFunc);
////	void saveManager(const char* f, Archive& ar, const OctTreeT* tree, const long partLevel, Communicator& commFunc);
//
//	/// This member serializes coefficients on client branch.
//	template <class Archive>
//	void shadowManager_save(const char* f, Archive& ar, const OctTreeT* tree,
//              Level n, Translation x, Translation y, Translation z, 
//              ProcessID remoteRank, const long partLevel, Communicator& commFunc);
////	void shadowManager_save(const char* f, Archive& ar, const OctTreeT* tree, Level n, Translation x, Translation y, Translation z, 
////        ProcessID remoteRank, const long partLevel, Communicator& commFunc);
//
//        /// This member serialize client's coeffcients.
//	template <class Archive>
//        void _shadowManager_save(const char* f, Archive& ar, const OctTreeT* tree,
//              std::vector<localTreeMember>& subtreeList, ProcessID remoteRank,
//              const long partLevel, Communicator& commFunc, const int nRemoteBranch,
//              int& iter); 
//
//	/// This member cotrols distributes coefficients to client branch.
//	template <class Archive>
//	void shadowManager_load(const char* f, const Archive& ar, OctTreeT* tree, 
//              Level n, Translation x, Translation y, Translation z, 
//              ProcessID remoteRank, Communicator& commFunc, const long partLevel, 
//	      bool active_flag, bool have_child);
//
//	/// This member is worker to send/recieve msg from rank 0.
//	void saveLoadWorker(OctTreeT* tree, 
//			Communicator& commFunc, bool save);
//
//	/// This member is member for DiskDir class.
//	void saveLoadWorker4DD(Communicator& comm, bool save){
//             saveLoadWorker(tree(), comm, save);
//        } 
//
//	/// This member Returns localsubtreelist of client computer to Master.
//	void localTreeList(std::vector<localTreeMember> &subtreeList, 
//			const OctTreeT* tree);
//
//	/// Send Coefficients to Rank 0. 
//	void sendRecvDataWorker_save(OctTreeT* tree, 
//			Communicator& commFunc);
//
//	/// Receive Coefficients from Rank 0. 
//	void sendRecvDataWorker_load(const OctTreeT* treeclient,  
//			Communicator& commFunc);
//
//	/// Loading Function members from the file. This member is not parallelized.
//	template <class Archive>
//	void load_local(const Archive& ar); 
//
//	/// Loading Function members from the file.
//        ///  This member is prepared for called from load_local. 
//	template <class Archive>
//	void _load_local(const Archive& ar, OctTreeT* tree);
//
//	/// Loading Function members from the file. This member is already parallelized.
//	//template <class Archive>
//	//void load(const Archive& iar, Communicator& comm);
//	void load(const char* f, Communicator& comm);
//
//	/// Load Managing Function member. 
//	template <class Archive>
//	void loadManager(const char* f, const Archive& ar, OctTreeT* tree,
//              Communicator& commFunc, bool active_flag, const long partLevel,
//              bool have_child);
////	void loadManager(const char* f, const Archive& ar, OctTreeT* tree, Communicator& commFunc, bool active_flag, const long partLevel, bool have_child);
//
//	/// Load Managing Function member for DiskDir class. 
//	template <class Archive>
//	void loadManager4DD(const Archive& ar, Communicator& commFunc, bool active_flag){
//             loadManager(ar, tree(), commFunc, active_flag, true);
//        }
//
//	/// Making Archive class's file name for 2-layer Serialization in sequential calculations.
//	void produceNewFilename(const char* f, const long partLevel,
//               const OctTreeTPtr& tree, char ftest[256]);
//
//	/// Making Archive class's file name for 2-layer Serialization(save, parallel).
//	void produceNewFilename2(const char* f, const long partLevel,
//               localTreeMember& subtreeList, char ftest[256]);
//
//	/// Making Archive class's file name for 2-layer Serialization(load, parallel).
//	void produceNewFilename3(const char* f, const long partLevel, Level n,
//               Translation x, Translation y, Translation z, char ftest[256]);
//
//        /// Save Client brach's data on master.
//	template <class Archive>
//        void dataSaveInShaManSave(const char* f, Archive& ar, const OctTreeT* tree,
//               localTreeMember& subtreeList, ProcessID remoteRank, const long partLevel,
//               Communicator& commFunc); 
////        void dataSaveInShaManSave(const char* f, Archive& ar, const OctTreeT* tree, localTreeMember& subtreeList, ProcessID remoteRank, const long partLevel, Communicator& commFunc); 

    	/// Inplace truncation of small wavelet coefficients
        
//        /// Works in the wavelet basis and compression is performed if the
//        /// function is not already compressed.  Communication streams up the
//        /// tree.  The default truncation threshold is that used to construct
//        /// the function.  If the threshold is zero, nothing is done.
//        /// Returns self for chaining.
//    	Function<T>& truncate(double tol = -1.0) {
//            if (tol < 0.0) tol = this->data->thresh;
//    	    if (tol == 0.0) return *this;
//            compress();
//            if (isactive(tree())) _truncate(tol, tree());
//            return *this;
//    	};
        
        T inner_local(const Function<T>& other) const {
            const_cast<Function<T>*>(this)->compress();
            const_cast<Function<T>*>(&other)->compress();
            return _inner_local(*this, other, tree());
        };
    
        /// Inner product between two Function classess. 
            
        /// Works in the wavelet basis.  Global communication is implied.
        T inner(const Function<T>& other) const {
            return comm()->global_sum(inner_local(other));
        };
    
    	/// This member is the member called from inner member.
    	T _inner_local(const Function<T>& a, const Function<T>& b, const OctTreeTPtr& tree) const ;

        /// Local evaluation of the function at a point in user coordinates

        /// If the coordinate is local to this process, returns true
        /// and sets value to the function value.  Otherwise, returns
        /// false and sets value to zero.
        ///
        /// If the function is not reconstructed an \em exception is
        /// thrown.  Also, node 0 throws an exception if the coordinate is
        /// outside the box.
        ///
        /// No communication is involved.
        bool eval_local(double x, double y, double z, T* value) const;


        /// Evaluate the function at a point in user coordinates

        /// !!! This is a \em collective operation.  In addition,
        /// if the function is not in the scaling function basis
        /// it is first reconstructed.
        ///
        /// On a sequential computer, this is a perfectly sensible
        /// function to use.  On a massively parallel computer, this
        /// is a shockingly inefficient way to evaluate anything more
        /// than a few values.  Instead, consider eval_local which
        /// does no communication, or eval_list which evaluates many
        /// values and only communicates/synchronizes once.
        T operator()(double x, double y, double z) const {
            T value;
            eval_list(1L, &x, &y, &z, &value);
            return value;
        };

        /// Evaluate function at a list of points, returning results to everyone

        /// Involves global communication and possible reconstruction.
        void eval_list(long n, const double* x, const double* y, const double*z,
                       T *value) const {
            if (iscompressed()) const_cast< Function<T>* >(this)->reconstruct();          
            for (long i = 0; i < n; i++) {
                eval_local(x[i], y[i], z[i], value + i);
            }
            comm()->global_sum(value, n);
        };

    private:
        /// Private.  Initialization required constructing function from stratch

        /// No communication involved.
        void _init(const FunctionFactory<T>& factory);


        /// Private.  Applies truncation method to give truncation threshold

        /// No communication involved.
        /// Truncate method m scales tol by 2^(-n*m/2)
        inline double truncate_tol(double tol, Level n) const {
            if (data->truncate_method == 0) {
                return tol;
            } else {
            	return tol*std::pow(0.5,0.5*n*data->truncate_method);
            }
        };


        /// Private.  Projects function in given box

        /// No communication involved.
        void _project(OctTreeTPtr& s);


        /// Private.  Projects function at given level.  No communication.

        /// No communication involved.  It is assumed that this is
        /// being done to an initially zero Function by _init
        /// immediately after construction
        void _fine_scale_projection(OctTreeTPtr& coeff, Level initial_level);


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
            transform3d_inplace(s, data->cdata->hgT, data->cdata->work2);
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
                MADNESS_EXCEPTION("unfilter: sonly : not yet",0);
            else
                return transform3d(ss, data->cdata->hg);
            //return transform(ss,hg);
        };


        /// Private. Returns the tree of coeffs

        /// No communication involved.
        inline OctTreeTPtr& tree() {
            MADNESS_ASSERT(data->tree->tree());
            return data->tree->tree();
        };


        /// Private. Returns the const tree of coeffs

        /// No communication involved.
        inline const OctTreeTPtr& tree() const {
            MADNESS_ASSERT(data->tree->tree());
            return data->tree->tree();
        };


        
        /// Private.  Retrieve pointer to coeffs of this Function from tree pointer

        /// No communication involved.
        inline TensorT* coeff(const OctTreeTPtr& tree) {
            MADNESS_ASSERT(tree);
            return tree->data(). template get<T>(ind);
        };

        /// Private.  Retrieve pointer to coeffs of this Function from tree pointer

        /// No communication involved.
        inline TensorT* coeff(const OctTreeT* tree) {
            MADNESS_ASSERT(tree);
            return tree->data(). template get<T>(ind);
        };


        /// Private.  Retrieve const pointer to coeffs of this Function from tree pointer

        /// No communication involved.
        inline const TensorT* coeff(const OctTreeTPtr& tree) const {
            MADNESS_ASSERT(tree);
            return tree->data(). template get<T>(ind);
        };


        /// Private.  Retrieve const pointer to coeffs of this Function from tree pointer

        /// No communication involved.
        inline const TensorT* coeff(const OctTreeT* tree) const {
            MADNESS_ASSERT(tree);
            return tree->data(). template get<T>(ind);
        };
        
        
        /// Private.  Unset any coeffs of this Function from tree pointer

        /// No communication involved.
        inline void unset_coeff(OctTreeTPtr& tree) {
            MADNESS_ASSERT(tree);
            tree->data().unset(ind);
        };


        /// Private.  Unset any coeffs of this Function from tree pointer

        /// No communication involved.
        inline void unset_coeff(OctTreeT* tree) {
            MADNESS_ASSERT(tree);
            tree->data().unset(ind);
        };


        /// Private.  Set coeffs of this Function from tree pointer

        /// No communication involved.
        inline TensorT* set_coeff(OctTreeTPtr& tree, const TensorT& t) {
            MADNESS_ASSERT(tree);
            return tree->data().set(ind, t);
        };

        /// Private.  Set coeffs of this Function from tree pointer

        /// No communication involved.
        inline TensorT* set_coeff(OctTreeT* tree, const TensorT& t) {
            MADNESS_ASSERT(tree);
            return tree->data().set(ind, t);
        };

        /// Private.  Set Function active in tree node

        /// No communication involved.
        inline void set_active(OctTreeTPtr& tree) const { // const ????
            MADNESS_ASSERT(tree);
            tree->data().set_active(ind);
        };

        /// Private.  Set Function active in tree node

        /// No communication involved.
        inline void set_active(OctTreeT* tree) const { // const ????
            MADNESS_ASSERT(tree);
            tree->data().set_active(ind);
        };


        /// Private.  Set Function inactive in tree node

        /// No communication involved.
        inline void set_inactive(OctTreeTPtr& tree) const { // const ????
            MADNESS_ASSERT(tree);
            tree->data().set_inactive(ind);
        };

        /// Private.  Set Function inactive in tree node

        /// No communication involved.
        inline void set_inactive(OctTreeT* tree) const { // const ????
            MADNESS_ASSERT(tree);
            tree->data().set_inactive(ind);
        };

    public:
        /// Private.  Determine if Function is active in tree node

        /// No communication involved.
        inline bool isactive(const OctTreeTPtr& tree) const {
            MADNESS_ASSERT(tree);
            return tree->data().isactive(ind);
        };

        /// No communication involved.
        inline bool isactive(const OctTreeT* tree) const {
            MADNESS_ASSERT(tree);
            return tree->data().isactive(ind);
        };

    private:
        /// Private.  For semantic similarity to isactive.

        /// No communication involved.
        inline bool islocal(const OctTreeTPtr& tree) const {
            MADNESS_ASSERT(tree);
            return tree->islocal();
        };

        /// No communication involved.
        inline bool islocal(const OctTreeT* tree) const {
            MADNESS_ASSERT(tree);
            return tree->islocal();
        };


        /// Private.  For semantic similarity to isactive.

        /// No communication involved.
        inline bool isremote(const OctTreeTPtr& tree) const {
            MADNESS_ASSERT(tree);
            return tree->isremote();
        };
        
        /// Private.  For semantic similarity to isactive.

        /// No communication involved.
        inline bool isremote(const OctTreeT* tree) const {
            MADNESS_ASSERT(tree);
            return tree->isremote();
        };
        
        /// Private.  For semantic similarity to isactive.

        /// No communication involved.
        inline void set_acflag(OctTreeTPtr& tree, bool value) {
            MADNESS_ASSERT(tree);
            tree->data().set_acflag(ind,value);
        };
        
        /// Private.  For semantic similarity to isactive.

        /// No communication involved.
        inline bool get_acflag(const OctTreeTPtr& tree) const {
            MADNESS_ASSERT(tree);
            return tree->data().get_acflag(ind);
        };
        

        /// Private.  Retrieve a pointer to the communicator

        /// No communication involved.
        inline Communicator* comm() const {
            return (Communicator *) (tree()->comm());
        };


        /// Private.  Recursive implementation of local 2-normsq

        /// No communication involved.
        double _norm2sq_local(const OctTreeTPtr& tree) const;

        /// Private.  Recursive function to refine from initial projection.

//        /// Communication streams all the way down the tree.
//        void _refine(OctTreeTPtr& tree);


        /// Private.  Recursive function to compress tree.

//        /// Communication streams up the tree.
//        void _compress(OctTreeTPtr& tree);


        /// Private.  Recursive function to reconstruct the tree.

        /// Communication streams down the tree.
        void _reconstruct(OctTreeTPtr& tree);


        /// Private.  Recursive function to support local unary operations
        /// including deep copy, scaling by constant, addition of constant,
        /// negation, etc.  Works in either basis.  No communication.
        /// The operation should be a function or functor that takes a
        /// reference to a const Tensor<T> as an argument and returns a
        /// new Tensor<T> as a result.
        template <class Op>
        void _unaryop(Function<T>& result, OctTreeTPtr& tree, const Op op) const {
            MADNESS_ASSERT(tree);
            result.set_active(tree);
            result.set_acflag(tree,false);
            const TensorT *t = coeff(tree);
            if (t) {
                Tensor<T> q = op(*t);
                result.set_coeff(tree, q);
            }

            FOREACH_CHILD(OctTreeTPtr, tree,
                          if (isactive(child))
                          _unaryop(result, child, op););
        }

//        /// Private:  Recursive kernel for truncation
//        void _truncate(double tol, OctTreeTPtr& tree);


        /// Private.  Recur down the tree printing out the norm of
        /// the coefficients.
        void _pnorms(const OctTreeTPtr& tree) const {
            MADNESS_ASSERT(tree);
            const TensorT *t = coeff(tree);
            for (long i=0; i<tree->n(); i++) std::printf("  ");
            if (t) {
                std::printf("%4d (%8lu, %8lu, %8lu) = %9.1e\n",
                            tree->n(),tree->x(),tree->y(),tree->z(),t->normf());
            } else {
                std::printf("%4d (%8lu, %8lu, %8lu) None\n",
                            tree->n(),tree->x(),tree->y(),tree->z());
            }
            FOREACH_CHILD(OctTreeTPtr, tree, if (isactive(child)) _pnorms(child););
        };

        /// Private.  Recursive function to support inplace scaling.
        void _scale(OctTreeTPtr& tree, T s) {
            MADNESS_ASSERT(tree);
            if (coeff(tree)) coeff(tree)->scale(s);
            FOREACH_CHILD(OctTreeTPtr, tree,
                          if (isactive(child)) _scale(child, s););
        };

        /// Private:  Compute norms of polyn with order <k/2 & >=k/2.

        /// No communication is involved, scaling fn basis only
        void _tnorms(const Tensor<T>& t, double& lo, double& hi) const;
        
        /// Returns list of unique processes with active remote children of a node in the tree
        
        /// Return value is the no. of unique remote processes with 
        /// active children of the node tree, and a list of their MPI ranks.
        /// No communication is involved.
        int unique_active_child_procs(const OctTreeTPtr& tree, ProcessID ranks[8]) const;
        
        
        void _ptree(OctTreeTPtr& tree) {
            MADNESS_ASSERT(tree);
            tree->print_coords();
            FOREACH_ACTIVE_CHILD(OctTreeTPtr, tree, _ptree(child););
        };
        
        /// Returns true if this block of scaling coefficients needs autorefining
        
        /// Criteria is that higher order polyn must be below threshold.
        /// Each additional level of refinement will reduce them by at least 0.5^k/2.
        inline bool autorefine_test(const OctTreeTPtr& tree) const {
            double lo, hi;
            MADNESS_ASSERT(tree);
            MADNESS_ASSERT(coeff(tree));
            _tnorms(*coeff(tree), lo, hi);
            return hi > this->truncate_tol(data->thresh,tree->n());
        };
        
        /// Returns true if this block of scaling coefficients needs autorefining
        
        /// Criteria is that the error in the square due to higher order polyn must be below threshold.
        /// Each additional level of refinement will reduce them by about 0.5^k.
        inline bool autorefine_square_test(const OctTreeTPtr& tree) const {
            if (!data->autorefine) return false;
            double lo, hi;
            MADNESS_ASSERT(tree);
            MADNESS_ASSERT(coeff(tree));
            _tnorms(*coeff(tree), lo, hi);
            return sqrt(hi*hi+2.0*hi*lo) > this->truncate_tol(data->thresh,tree->n());
        };
        
        
        /// Returns true if this block of scaling coefficients needs autorefining
        
        /// Criteria is that the error in the square due to higher order polyn must be below threshold.
        /// Each additional level of refinement will reduce them by about 0.5^k.
        inline bool autorefine_mult_test(const OctTreeTPtr& tree, const Tensor<T>& a, const Tensor<T>& b) const {
            if (!data->autorefine) return false;
            double alo, ahi, blo, bhi;
            _tnorms(a, alo, ahi);
            _tnorms(b, blo, bhi);
            return sqrt(ahi*bhi+ahi*blo+blo*ahi) > this->truncate_tol(data->thresh,tree->n());
        };


        /// Private.  Squares a block of coefficients
        void do_square(OctTreeTPtr& tree);
        
        /// Private.  Multiply two blocks of coeffs and set the result
        void do_mult(OctTreeTPtr& tree, const Tensor<T>& a, const Tensor<T>& b);
        
        /// Private.  Transforms in place coeffs to function values on quadrature grid
        void do_transform_to_function_values(OctTreeTPtr& tree);

        /// Private.  Transforms in place from function values on quadrature grid to coeffs
        void do_transform_from_function_values(OctTreeTPtr& tree);

        
        /// Private.  Recursive function to support gaxpy
        void _gaxpy(Function<T>& afun, double alpha,
                    const Function<T>& bfun, double beta, OctTreeTPtr& tree);
                    

        /// Private.  Recursive function to support addition and subtraction
        void _add_sub(Function<T>& c, const Function<T>& a, const Function<T>& b,
                      OctTreeTPtr& tree, bool subtract) const;


//        /// Private.  Check that this function and another share the same OctTree
//        void check_trees(const Function<T>& other) const {
//            if (tree() != other.tree()) 
//                MADNESS_EXCEPTION("Function: check_trees: trees are different",0);
//        };


        /// Private.  Evaluates function in scaling function basis within a cube
        T _eval_cube(Level n,
                     double xx, double yy, double zz,
                     const Tensor<T>& s) const;
        

        
        void _auto_clean(OctTreeTPtr& tree) {
            if (isactive(tree)) {
                if (get_acflag(tree)) {
                    unset_coeff(tree);
                    set_inactive(tree);
                    set_acflag(tree,false);
                };
                FOREACH_CHILD(OctTreeTPtr, tree, _auto_clean(child););
            }
        };
        
        void build_global_tree() {
            data->gtree = madness::build_global_tree(tree());
        };
        
        void clear_global_tree() {
            data->gtree = GlobalTree<3>();
        };
        
        ProcessID find_owner(Level n, const Translation l[3]) const {
            return data->gtree.find_owner(n,l);
        };
                           
    };

    /// Deep copy function ... invokes f.copy()
    template <typename T>
    inline Function<T> copy(const Function<T>& f) {
        return f.copy();
    }


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
}

#endif
