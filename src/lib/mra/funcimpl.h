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
#include <map>
#include <tr1/tuple>
#include <world/world.h>
#include <world/print.h>
#include <world/scopedptr.h>
#include <misc/misc.h>
#include <tensor/tensor.h>
#include <mra/key.h>
#include <mra/funcdefaults.h>
#include <world/GPU_streams.h>
//#include <pthread.h>
//#include <world/cuda_streams.h>
static double tttt, ssss;
#define STARTt_TIMER  tttt=wall_time(); ssss=cpu_time()
#define ENDt_TIMER(msg) tttt=wall_time()-tttt; ssss=cpu_time()-ssss;  printf("timer: %20.20s %8.10fs %8.10fs\n", msg, ssss, tttt)

static pthread_barrier_t barr;

#include<math.h>
extern "C" void streams_synchronize(void **,unsigned int);
namespace madness {

#define IS_ODD(n) ((n)&0x1)
#define IS_UNALIGNED(p) (((unsigned long)(p))&0x7)

    template <typename T, std::size_t NDIM>
    class DerivativeBase;

    template<typename T, std::size_t NDIM>
    class FunctionImpl;

    template <typename T, std::size_t NDIM>
    class TensorNode;

    template<typename T, std::size_t NDIM>
    class Function;
    /*
    template<typename T, std::size_t NDIM>
    class FunctionNode;
    */
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

    /*
    class GlobalBarrier {
        public:
        static pthread_barrier_t* barr[1];

        static pthread_barrier_t* barrier(){
            if (!barr[0]) pthread_barrier_init(barr[0], NULL, NUMPAR_MEMCPYS);
            return (barr[0]);
        }
    };
    */

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

        Tensor<double> h0, h1, g0, g1; ///< The separate blocks of twoscale coefficients
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

    /// Interface required for functors used as input to Functions
    template<typename T, std::size_t NDIM>
    class FunctionFunctorInterface {
    public:
        /// You should implement this to return \c f(x)
        virtual T operator()(const Vector<double, NDIM>& x) const = 0;

        /// Override this to return list of special points to be refined more deeply
        virtual std::vector< Vector<double,NDIM> > special_points() const {
            return std::vector< Vector<double,NDIM> >();
        }

        /// Override this change level refinement for special points (default is 6)
        virtual Level special_level() {return 6;}

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
    template<typename T, std::size_t NDIM>
    class FunctionFactory {
        friend class FunctionImpl<T, NDIM> ;
        typedef Vector<double, NDIM> coordT; ///< Type of vector holding coordinates
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
        //Tensor<int> _bc;
        std::shared_ptr<WorldDCPmapInterface<Key<NDIM> > > _pmap;
        std::shared_ptr<FunctionFunctorInterface<T, NDIM> > _functor;

    public:
        struct FunctorInterfaceWrapper : public FunctionFunctorInterface<T,NDIM> {
            T (*f)(const coordT&);

            FunctorInterfaceWrapper(T (*f)(const coordT&)) : f(f) {}

            T operator()(const coordT& x) const {return f(x);}
        };

        FunctionFactory(World& world) :
                _world(world),
                _k(FunctionDefaults<NDIM>::get_k()),
                _thresh(FunctionDefaults<NDIM>::get_thresh()),
                _initial_level(
                    FunctionDefaults<NDIM>::get_initial_level()),
                _max_refine_level(
                    FunctionDefaults<NDIM>::get_max_refine_level()),
                _truncate_mode(
                    FunctionDefaults<NDIM>::get_truncate_mode()),
                _refine(FunctionDefaults<NDIM>::get_refine()),
                _empty(false),
                _autorefine(FunctionDefaults<NDIM>::get_autorefine()),
                _truncate_on_project(
                    FunctionDefaults<NDIM>::get_truncate_on_project()),
                _fence(true), // _bc(FunctionDefaults<NDIM>::get_bc()),
                _pmap(FunctionDefaults<NDIM>::get_pmap()), _functor() {
        }
        FunctionFactory&
        functor(
            const std::shared_ptr<FunctionFunctorInterface<T, NDIM> >& f) {
            _functor = f;
            return *this;
        }
        FunctionFactory&
        f(T
          (*f)(const coordT&)) {
            functor(std::shared_ptr<FunctionFunctorInterface<T, NDIM> > (
                        new FunctorInterfaceWrapper(f)));
            return *this;
        }
        FunctionFactory&
        k(int k) {
            _k = k;
            return *this;
        }
        FunctionFactory&
        thresh(double thresh) {
            _thresh = thresh;
            return *this;
        }
        FunctionFactory&
        initial_level(int initial_level) {
            _initial_level = initial_level;
            return *this;
        }
        FunctionFactory&
        max_refine_level(int max_refine_level) {
            _max_refine_level = max_refine_level;
            return *this;
        }
        FunctionFactory&
        truncate_mode(int truncate_mode) {
            _truncate_mode = truncate_mode;
            return *this;
        }
        FunctionFactory&
        refine(bool refine = true) {
            _refine = refine;
            return *this;
        }
        FunctionFactory&
        norefine(bool norefine = true) {
            _refine = !norefine;
            return *this;
        }

        FunctionFactory&
        empty() {
            _empty = true;
            return *this;
        }
        FunctionFactory&
        autorefine() {
            _autorefine = true;
            return *this;
        }
        FunctionFactory&
        noautorefine() {
            _autorefine = false;
            return *this;
        }
        FunctionFactory&
        truncate_on_project() {
            _truncate_on_project = true;
            return *this;
        }
        FunctionFactory&
        notruncate_on_project() {
            _truncate_on_project = false;
            return *this;
        }
        FunctionFactory&
        fence(bool fence = true) {
            _fence = fence;
            return *this;
        }
        FunctionFactory&
        nofence() {
            _fence = false;
            return *this;
        }
        FunctionFactory&
        pmap(const std::shared_ptr<WorldDCPmapInterface<Key<NDIM> > >& pmap) {
            _pmap = pmap;
            return *this;
        }
    };

    //used for multi-threaded memcpy in compressop_allCompute
    template<typename T, std::size_t NDIM>
    struct MemcpyCompressData {
        T* dptr;
        const double* pc;
        T* tstart;
        double* dstart;
        long tbuf_stepbytes;
        long dbuf_stepbytes;
    };

    /// FunctionNode holds the coefficients, etc., at each node of the 2^NDIM-tree
    template<typename T, std::size_t NDIM>
    class FunctionNode {
    private:
        // Should compile OK with these volatile but there should
        // be no need to set as volatile since the container internally
        // stores the entire entry as volatile

        Tensor<T> _coeffs; ///< The coefficients, if any
        double _norm_tree; ///< After norm_tree will contain norm of coefficients summed up tree
        bool _has_children; ///< True if there are children

    public:
        typedef WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > dcT; ///< Type of container holding the nodes
        typedef Key<NDIM> keyT;
        typedef Tensor<T> tensorT;
        typedef WorldContainer<keyT, TensorNode<T,NDIM>> containerT;
        /// Default constructor makes node without coeff or children
        FunctionNode() :
                _coeffs(), _norm_tree(1e300), _has_children(false) {
        }

        /// Constructor from given coefficients with optional children

        /// Note that only a shallow copy of the coeff are taken so
        /// you should pass in a deep copy if you want the node to
        /// take ownership.
        explicit
        FunctionNode(const Tensor<T>& coeff, bool has_children = false) :
                _coeffs(coeff), _norm_tree(1e300), _has_children(has_children) {
        }

        explicit
        FunctionNode(const Tensor<T>& coeff, double norm_tree, bool has_children) :
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
            return (_coeffs.size() > 0);
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
        Tensor<T>&
        coeff() {
            MADNESS_ASSERT(_coeffs.ndim() == -1 || (_coeffs.dim(0) <= 2
                                                    * MAXK && _coeffs.dim(0) >= 0));
            return const_cast<Tensor<T>&>(_coeffs);
        }

        /// Returns a const reference to the tensor containing the coeffs

        /// Returns an empty tensor if there are no coefficeints.
        const Tensor<T>&
        coeff() const {
            return const_cast<const Tensor<T>&>(_coeffs);
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
        void set_coeff(const Tensor<T>& coeffs) {
            coeff() = coeffs;
            if ((_coeffs.dim(0) < 0) || (_coeffs.dim(0)>2*MAXK)) {
                print("set_coeff: may have a problem");
                print("set_coeff: coeff.dim[0] =", coeffs.dim(0), ", 2* MAXK =", 2*MAXK);
            }
            MADNESS_ASSERT(coeffs.dim(0)<=2*MAXK && coeffs.dim(0)>=0);
        }

        /// Clears the coefficients (has_coeff() will subsequently return false)
        void clear_coeff() {
            coeff() = Tensor<T>();
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

        tensorT unfilter(const tensorT& s, const FunctionCommonData<T,NDIM>& cdata) const {
            tensorT r(cdata.v2k,false);
            tensorT w(cdata.v2k,false);
            return fast_transform(s,cdata.hg,r,w);
            //return transform(s, cdata.hg);
        }

        tensorT filter(const tensorT& s, const FunctionCommonData<T,NDIM>& cdata) const {
            tensorT r(cdata.v2k,false);
            tensorT w(cdata.v2k,false);
            return fast_transform(s,cdata.hgT,r,w);
            //return transform(s,cdata.hgT);
        }

        std::vector<Slice> child_patch(const keyT& child, FunctionCommonData<T,NDIM> cdata) const {
            std::vector<Slice> s(NDIM);
            const Vector<Translation,NDIM>& l = child.translation();
            for (std::size_t i=0; i<NDIM; ++i)
                s[i] = cdata.s[l[i]&1]; // Lowest bit of translation
            return s;
        }

        Void reconstruct_dc_task(const keyT& key, dcT dc, const tensorT& s, int k);

        Void top_down(const keyT& key, dcT dc, containerT ct, bool nonstandard, bool keepleaves, const keyT& parent, int k){
          keyT& parent_t = const_cast<keyT&>(parent);
          ct.update(key, &TensorNode<T,NDIM>::compress_dc, dc, parent_t, nonstandard, k);

          if (this->has_children()) {
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                dc.update(kit.key(), &FunctionNode<T,NDIM>::top_down, ct, nonstandard, keepleaves, key, k);
            }
          }
          else {
            tensorT result(this->coeff());
            if (!keepleaves) this->clear_coeff();
            //keyT& key_t = const_cast<keyT&>(key);
            ct.update(parent_t, &TensorNode<T,NDIM>::update_map, dc, key, result, nonstandard, k);
          }

          return None;
        }

          std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT> compressop_preprocess(//const keyT& key, dcT dc,
            const std::tr1::tuple< keyT, dcT, std::map<keyT, tensorT>, containerT, bool, int, keyT >& in
//const std::map<keyT, tensorT>& tensor_keys, containerT ct, const bool& nonstandard, const int& k, const keyT& parent
){
            // Copy child scaling coeffs into contiguous block
            keyT key = std::tr1::get<0>(in);
            dcT dc = std::tr1::get<1>(in);
            std::map<keyT, tensorT> tensor_keys = std::tr1::get<2>(in);
            containerT ct = std::tr1::get<3>(in); 
            bool nonstandard = std::tr1::get<4>(in); 
            int k = std::tr1::get<5>(in);
            keyT parent = std::tr1::get<6>(in); 
                        
            FunctionCommonData<T,NDIM> cdata = FunctionCommonData<T,NDIM>::get(k);
            std::map<keyT, tensorT>& tk = const_cast<std::map<keyT, tensorT>&>(tensor_keys);
            tensorT d(cdata.v2k,false);
            int i=0;
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
                d(child_patch(kit.key(),cdata)) = tk[kit.key()];
            }

            int k1 = k;
            bool nonstandard1 = nonstandard;
            keyT parent1 = parent;

           std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT> t1(d,k1,key,ct,nonstandard1,parent1,dc);
           return t1;
}

          std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT,long,long,Tensor<double>*,T*,T*,tensorT*,tensorT*,const double*> compressop_preprocessGPU(//const keyT& key, dcT dc,
            const std::tr1::tuple< keyT, dcT, std::map<keyT, tensorT>, containerT, bool, int, keyT >& in
//const std::map<keyT, tensorT>& tensor_keys, containerT ct, const bool& nonstandard, const int& k, const keyT& parent
){
            // Copy child scaling coeffs into contiguous block
            keyT key = std::tr1::get<0>(in);
            dcT dc = std::tr1::get<1>(in);
            std::map<keyT, tensorT> tensor_keys = std::tr1::get<2>(in);
            containerT ct = std::tr1::get<3>(in); 
            bool nonstandard = std::tr1::get<4>(in); 
            int k = std::tr1::get<5>(in);
            keyT parent = std::tr1::get<6>(in); 
                        
            FunctionCommonData<T,NDIM> cdata = FunctionCommonData<T,NDIM>::get(k);
            std::map<keyT, tensorT>& tk = const_cast<std::map<keyT, tensorT>&>(tensor_keys);
            tensorT d(cdata.v2k,false);
            int i=0;
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
                d(child_patch(kit.key(),cdata)) = tk[kit.key()];
            }

            int k1 = k;
            bool nonstandard1 = nonstandard;
            keyT parent1 = parent;

           //std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT> t1(d,k1,key,ct,nonstandard1,parent1,dc);
           
           //FunctionCommonData<T, NDIM> cdata = FunctionCommonData<T, NDIM>::get(k);
           //tensorT r = new Tensor(cdata.v2k,false);
           //tensorT w = new Tensor(cdata.v2k,false);

           tensorT t = d;
           Tensor<double> * c = new Tensor<double>(cdata.hgT);
           tensorT * result = new tensorT(cdata.v2k,false);
           tensorT * workspace = new tensorT(cdata.v2k,false);

           typedef T resultT;       
        
           const double *pc=c->ptr();
           resultT *t0=workspace->ptr(), *t1=result->ptr();
           if (t.ndim()&1) {
             t0 = result->ptr();
             t1 = workspace->ptr();
           }

           long dimj = c->dim(1);
           long dimi = 1;
           for (int n=1; n<t.ndim(); ++n) dimi *= dimj;
           //long nij = dimi*dimj;

           std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT,long,long,Tensor<double>*,resultT*,resultT*,tensorT*,tensorT*,const double*> tu1(d,k1,key,ct,nonstandard1,parent1,dc,dimi,dimj,c,t0,t1,result,workspace,pc);

           return tu1;
        //return fast_transform(d,cdata.hgT,r,w);

           //return t1;
        }

        Void compressop_backToCPU(std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT,long,long,Tensor<double>*,T*,T*,tensorT*,tensorT*,const double*> in){
        tensorT d = std::tr1::get<0>(in);
        int k = std::tr1::get<1>(in);
        FunctionCommonData<T, NDIM> cdata = FunctionCommonData<T, NDIM>::get(k);
        keyT key = std::tr1::get<2>(in);
        containerT ct = std::tr1::get<3>(in); 
        bool nonstandard = std::tr1::get<4>(in); 
        keyT parent = std::tr1::get<5>(in); 
        dcT dc = std::tr1::get<6>(in); 
        long dimi = std::tr1::get<7>(in);
        long dimj = std::tr1::get<8>(in);
        Tensor<double>* c1 = std::tr1::get<9>(in);
        T* t0 = std::tr1::get<10>(in);
        T* t1 = std::tr1::get<11>(in);
        tensorT* result = std::tr1::get<12>(in);
        tensorT* workspace = std::tr1::get<13>(in);
        const double* pc = std::tr1::get<14>(in);
        
        delete c1;
        delete result;
        delete workspace;

        d = filter(d,cdata);

        if (this->has_coeff()) {
            const tensorT& c = this->coeff();
            if (c.dim(0) == k) {
                d(cdata.s0) += c;
            }
            else {
                d += c;
            }
        }

        tensorT s = copy(d(cdata.s0));

        if (key.level()> 0 && !nonstandard)
            d(cdata.s0) = 0.0;

        this->set_coeff(d);

        if (key.level() > 0){
            ct.update(parent, &TensorNode<T, NDIM>::update_map, dc, key, s, nonstandard, k);
        }

        return None;
   }
         
        

        std::tr1::tuple<tensorT*,int,keyT,containerT,bool,keyT,dcT,tensorT*,Tensor<double>*> compressop_fasttransform(std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT,long,long,Tensor<double>*,T*,T*,tensorT*,tensorT*,const double*> in){
        tensorT d = std::tr1::get<0>(in);
        int k = std::tr1::get<1>(in);
        FunctionCommonData<T, NDIM> cdata = FunctionCommonData<T, NDIM>::get(k);
        keyT key = std::tr1::get<2>(in);
        containerT ct = std::tr1::get<3>(in); 
        bool nonstandard = std::tr1::get<4>(in); 
        keyT parent = std::tr1::get<5>(in); 
        dcT dc = std::tr1::get<6>(in); 
        long dimi = std::tr1::get<7>(in);
        long dimj = std::tr1::get<8>(in);
        Tensor<double>* c = std::tr1::get<9>(in);
        T* t0 = std::tr1::get<10>(in);
        T* t1 = std::tr1::get<11>(in);
        tensorT* result = std::tr1::get<12>(in);
        tensorT* workspace = std::tr1::get<13>(in);
        const double* pc = std::tr1::get<14>(in);

        tensorT t = d;
        long nij = dimi*dimj;
        if (IS_ODD(dimi) || IS_ODD(dimj) ||
                IS_UNALIGNED(pc) || IS_UNALIGNED(t0) || IS_UNALIGNED(t1)) {
            for (long i=0; i<nij; ++i) t0[i] = 0.0;
            mTxm(dimi, dimj, dimj, t0, t.ptr(), pc);

           //     cu_mTxmq(dimi, dimj, dimj, t0, t.ptr(), pc,GPU_streams[0],1,t.size());
           for (int n=1; n<t.ndim(); ++n) {
                for (long i=0; i<nij; ++i) t1[i] = 0.0;
                mTxm(dimi, dimj, dimj, t1, t0, pc);

           //     cu_mTxmq(dimi, dimj, dimj, t1, t0, pc,GPU_streams[0],t.ndim(),t.size());
                std::swap(t0,t1);
           }
        }
        else {
            mTxmq(dimi, dimj, dimj, t0, t.ptr(), pc);
           // print("CUDA KERNEL (dim = ",dimi,",",dimj,")\n");
           // cu_mTxmq(dimi, dimj, dimj, t0, t.ptr(), pc,GPU_streams[0],1,t.size());
            for (int n=1; n<t.ndim(); ++n) {
                mTxmq(dimi, dimj, dimj, t1, t0, pc);
             //   cu_mTxmq(dimi, dimj, dimj, t1, t0, pc,GPU_streams[0],t.ndim(),t.size());
                std::swap(t0,t1);
            }
        }

          std::tr1::tuple<tensorT*,int,keyT,containerT,bool,keyT,dcT,tensorT*,Tensor<double>*> t11(result,k,key,ct,nonstandard,parent,dc,workspace,c) ;
          return t11 ;
}

        std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT> compressop_compute(std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT> in){
          //std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT> in = in1.get();
          tensorT d = std::tr1::get<0>(in);
          int k = std::tr1::get<1>(in);
          FunctionCommonData<T, NDIM> cdata = FunctionCommonData<T, NDIM>::get(k);
          keyT key = std::tr1::get<2>(in);
          containerT ct = std::tr1::get<3>(in); 
          bool nonstandard = std::tr1::get<4>(in); 
          keyT parent = std::tr1::get<5>(in); 
          dcT dc = std::tr1::get<6>(in); 
          d = filter(d,cdata);
          std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT> t11(d,k,key,ct,nonstandard,parent,dc) ;
          //Future< std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT> > t1(t11) ;
          return t11;

          //dc.task(key, &FunctionNode<T,NDIM>::compressop_postprocess, t11);

          //return None;
        }

        Void compressopComputePostprocess(std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT> in){
          //std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT> in = in1.get();
          tensorT d = std::tr1::get<0>(in);
          int k = std::tr1::get<1>(in);
          FunctionCommonData<T, NDIM> cdata = FunctionCommonData<T, NDIM>::get(k);
          keyT key = std::tr1::get<2>(in);
          containerT ct = std::tr1::get<3>(in); 
          bool nonstandard = std::tr1::get<4>(in); 
          keyT parent = std::tr1::get<5>(in); 
          dcT dc = std::tr1::get<6>(in); 
          d = filter(d,cdata);

          if (this->has_coeff()) {
                const tensorT& c = this->coeff();
                if (c.dim(0) == k) {
                    d(cdata.s0) += c;
                }
                else {
                    d += c;
                }
            }

            tensorT s = copy(d(cdata.s0));

            if (key.level()> 0 && !nonstandard)
                d(cdata.s0) = 0.0;

            this->set_coeff(d);

            if (key.level() > 0){
              ct.update(parent, &TensorNode<T, NDIM>::update_map, dc, key, s, nonstandard, k);
            }

            return None;
        }

    /*   template <typename W>
       W* GPUtransfer_buffer(W* buf, unsigned int offset){
         return buf;
       }

       template <typename W>
       void GPUdelete_buffer(W* buf){}
*/
        std::vector< std::tr1::tuple<tensorT*,int,keyT,containerT,bool,keyT,dcT,tensorT*,Tensor<double>*> > compressop_allComputeGPU(std::vector< std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT,long,long,Tensor<double>*,T*,T*,tensorT*,tensorT*,const double*> > inArgs, std::vector< FunctionNode<T,NDIM>* > inObj){
            std::vector< std::tr1::tuple<tensorT*,int,keyT,containerT,bool,keyT,dcT,tensorT*,Tensor<double>*> > outArg(inArgs.size(),inObj.at(0)->compressop_fasttransform(inArgs.at(0)));
            print("inArgs.size() = ",inArgs.size());
STARTt_TIMER;
            long op_dimi = std::tr1::get<7>(inArgs.at(0));
            long op_dimj = std::tr1::get<8>(inArgs.at(0));
            long op_nij = op_dimi*op_dimj;
            long tbuf_size = inArgs.size()*op_dimi*op_dimj;
            long dbuf_size = inArgs.size()*op_dimj*op_dimj;
            //long tbuf_sizebytes = inArgs.size()*op_dimi*op_dimj*sizeof(T);
            //long dbuf_sizebytes = inArgs.size()*op_dimj*op_dimj*sizeof(double);
            long tbuf_step = op_dimi*op_dimj;
            long dbuf_step = op_dimj*op_dimj;
            long tbuf_stepbytes = op_dimi*op_dimj*sizeof(T);
            long dbuf_stepbytes = op_dimj*op_dimj*sizeof(double);
            /*
            unsigned int t0_off = 0;
            unsigned int tptr_off = 0;
            unsigned int pc_off = 0;
            unsigned int t1_off = 0;
            for (unsigned int i = 0; i < inArgs.size(); i++){
		long dimi = std::tr1::get<7>(inArgs.at(i));
		long dimj = std::tr1::get<8>(inArgs.at(i));
                print("dimi = ",dimi," dimj = ",dimj);
                //memcpy(t0_buf+t0_off*sizeof(T),t0,dimi*dimj*sizeof(T));
                t0_off += dimi*dimj;
                //memcpy(tptr_buf+tptr_off*sizeof(T),tptr,dimi*dimj*sizeof(T));
                tptr_off += dimi*dimj;
                //memcpy(pc_buf+pc_off*sizeof(T),pc,dimj*dimj*sizeof(T));
                pc_off += dimj*dimj;
                //memcpy(t1_buf+t1_off*sizeof(T),t1,dimi*dimj*sizeof(T));
                t1_off += dimi*dimj;
            }
            //These memory allocations are expensive (hundreds of microsec)
            T* t0_buf = new T[t0_off];
            T* tptr_buf = new T[tptr_off];
            double* pc_buf = new double[pc_off];
            T* t1_buf = new T[t1_off];
            */
            T* t0_buf = new T[/*inArgs.size()*dimi*dimj*/tbuf_size];
            T* tptr_buf = new T[/*inArgs.size()*dimi*dimj*/tbuf_size];
            double* pc_buf = new double[/*inArgs.size()*dimj*dimj*/dbuf_size];
            T* t1_buf = new T[/*inArgs.size()*dimi*dimj*/tbuf_size];
ENDt_TIMER("initial_access");


            std::vector<T*> dptrs(inArgs.size());
            std::vector<const double*> pcs(inArgs.size());
            for (unsigned int i = 0; i < inArgs.size(); i++){
                tensorT d = std::tr1::get<0>(inArgs.at(i));
		dptrs[i] = d.ptr();
                pcs[i] = std::tr1::get<14>(inArgs.at(i));
            }

           cpu_set_t cpuset[8];

           /* Set affinity mask to include CPUs 0 to 7 */

            for (unsigned int j = 0; j < 8; j++) {
               CPU_ZERO(&cpuset[j]);
               CPU_SET(j, &cpuset[j]);
           }

            void * status;
            pthread_t memcpy_threads[NUMPAR_MEMCPYS];   
            MemcpyCompressData<T,NDIM> memcpy_data[NUMPAR_MEMCPYS];            

            for (unsigned int i = 0; i < inArgs.size(); i+=NUMPAR_MEMCPYS){
                unsigned int j;  
                for (j = i; j < inArgs.size() && j < i+NUMPAR_MEMCPYS; j++){
                    memcpy_data[j%NUMPAR_MEMCPYS].dptr = dptrs[j];
                    memcpy_data[j%NUMPAR_MEMCPYS].pc = pcs[j];
                    memcpy_data[j%NUMPAR_MEMCPYS].tstart = tptr_buf + j * tbuf_step;
                    memcpy_data[j%NUMPAR_MEMCPYS].dstart = pc_buf + j * dbuf_step;
                    memcpy_data[j%NUMPAR_MEMCPYS].tbuf_stepbytes = tbuf_stepbytes;
                    memcpy_data[j%NUMPAR_MEMCPYS].dbuf_stepbytes = dbuf_stepbytes;
                }
                STARTt_TIMER;
                
                for (j = i; j < i+NUMPAR_MEMCPYS; j++){
                    int ret = pthread_create( &memcpy_threads[j%NUMPAR_MEMCPYS], NULL, &FunctionImpl<T,NDIM>::blank, &memcpy_data[j%NUMPAR_MEMCPYS]);
                    pthread_setschedprio( memcpy_threads[j%NUMPAR_MEMCPYS], 1000 );
                    //long s = pthread_setaffinity_np( memcpy_threads[j%NUMPAR_MEMCPYS], sizeof(cpu_set_t), &cpuset[j%8]);
                    //if (s != 0)
                    //    perror("pthread_setaffinity_np");
                    //memcpy(tptr_buf+/*tptr_off*/i*tbuf_step,dptrs[i],/*op_dimi*op_dimj*sizeof(T)*/tbuf_stepbytes);
                    //memcpy(pc_buf+/*pc_off*/i*dbuf_step,pcs[i],/*op_dimj*op_dimj*sizeof(double)*/dbuf_stepbytes);
                }
                for (j = i; j < i+NUMPAR_MEMCPYS; j++){
                    pthread_join( memcpy_threads[j%NUMPAR_MEMCPYS], &status );
                }
                ENDt_TIMER("thread create and join");
                STARTt_TIMER;
                  
                for (j = i; j < inArgs.size() && j < i+NUMPAR_MEMCPYS; j++){
                    int ret = pthread_create( &memcpy_threads[j%NUMPAR_MEMCPYS], NULL, &FunctionImpl<T,NDIM>::memcpyCompress,&memcpy_data[j%NUMPAR_MEMCPYS]);
                    pthread_setschedprio( memcpy_threads[j%NUMPAR_MEMCPYS], 1000 );
                    //long s = pthread_setaffinity_np( memcpy_threads[j%NUMPAR_MEMCPYS], sizeof(cpu_set_t), &cpuset[j%8]);
                    //if (s != 0)
                    //    perror("pthread_setaffinity_np");
                    //memcpy(tptr_buf+/*tptr_off*/i*tbuf_step,dptrs[i],/*op_dimi*op_dimj*sizeof(T)*/tbuf_stepbytes);
                    //memcpy(pc_buf+/*pc_off*/i*dbuf_step,pcs[i],/*op_dimj*op_dimj*sizeof(double)*/dbuf_stepbytes);
                }
                for (unsigned int ii = j; ii < i+NUMPAR_MEMCPYS; ii++){
                    int ret = pthread_create( &memcpy_threads[ii%NUMPAR_MEMCPYS], NULL, &FunctionImpl<T,NDIM>::blank,&memcpy_data[ii%NUMPAR_MEMCPYS]);
                    //long s = pthread_setaffinity_np( memcpy_threads[j%NUMPAR_MEMCPYS], sizeof(cpu_set_t), &cpuset[ii%8]);
                    //if (s != 0)
                    //    perror("pthread_setaffinity_np");
                    pthread_setschedprio( memcpy_threads[ii%NUMPAR_MEMCPYS], 1000 );

                }
                for (j = i; j < i+NUMPAR_MEMCPYS; j++){
                    pthread_join( memcpy_threads[j%NUMPAR_MEMCPYS], &status );
                }
                ENDt_TIMER("memcpys");
            }

            //t0_off = 0;
            //tptr_off = 0;
            //pc_off = 0;
            //t1_off = 0;
            
            //for (unsigned int i = 0; i < inArgs.size(); i++){
                //STARTt_TIMER;
		//tensorT d = std::tr1::get<0>(inArgs.at(i));
		//int k = std::tr1::get<1>(inArgs.at(i));
		//FunctionCommonData<T, NDIM> cdata = FunctionCommonData<T, NDIM>::get(k);
		//keyT key = std::tr1::get<2>(inArgs.at(i));
		//containerT ct = std::tr1::get<3>(inArgs.at(i)); 
		//bool nonstandard = std::tr1::get<4>(inArgs.at(i)); 
		//keyT parent = std::tr1::get<5>(inArgs.at(i)); 
		//dcT dc = std::tr1::get<6>(inArgs.at(i)); 
		//long dimi = std::tr1::get<7>(inArgs.at(i));
		//long dimj = std::tr1::get<8>(inArgs.at(i));
		//Tensor<double>* c = std::tr1::get<9>(inArgs.at(i));
		//T* t0 = std::tr1::get<10>(inArgs.at(i));
		//T* t1 = std::tr1::get<11>(inArgs.at(i));
		//tensorT* result = std::tr1::get<12>(inArgs.at(i));
		//tensorT* workspace = std::tr1::get<13>(inArgs.at(i));
		//const double* pc = std::tr1::get<14>(inArgs.at(i));
                //ENDt_TIMER("loop_access");

                //memcpy(t0_buf+t0_off,t0,dimi*dimj*sizeof(T));
                //print("memcpy: ",dimi*dimj*sizeof(T)," bytes, ",dimj*dimj*sizeof(T)," bytes");
                //STARTt_TIMER;
                //memcpy(tptr_buf+/*tptr_off*/i*tbuf_step,d.ptr(),/*op_dimi*op_dimj*sizeof(T)*/tbuf_stepbytes);
                //ENDt_TIMER("memcpy 1");
                //STARTt_TIMER;
                //memcpy(pc_buf+/*pc_off*/i*dbuf_step,pc,/*op_dimj*op_dimj*sizeof(double)*/dbuf_stepbytes);
                //ENDt_TIMER("memcpy 2");
                //t0_off += dimi*dimj;
                //tptr_off += dimi*dimj;
                //pc_off += dimj*dimj;
                //memcpy(t1_buf+t1_off,t1,dimi*dimj*sizeof(T));
                //t1_off += dimi*dimj;
            //}

STARTt_TIMER;
            T * start_t0 = GPUtransfer_buffer(t0_buf,/*t0_off*/tbuf_size,false);
            T * start_tptr = GPUtransfer_buffer(tptr_buf,/*tptr_off*/tbuf_size,true);
            double * start_pc = GPUtransfer_buffer(pc_buf,/*pc_off*/dbuf_size,true);               
            T * start_t1 = GPUtransfer_buffer(t1_buf,/*t1_off*/tbuf_size,false);
ENDt_TIMER("transfer"); 

            //t0_off = 0;
            //tptr_off = 0;
            //pc_off = 0;
            //t1_off = 0;

STARTt_TIMER;
            for (unsigned int i = 0; i < inArgs.size(); i++){
		tensorT d = std::tr1::get<0>(inArgs.at(i));
		tensorT t = d;
		//long dimi = std::tr1::get<7>(inArgs.at(i));
		//long dimj = std::tr1::get<8>(inArgs.at(i));
		//long nij = dimi*dimj;
                /*
		if (IS_ODD(dimi) || IS_ODD(dimj) ||
			IS_UNALIGNED(pc) || IS_UNALIGNED(t0) || IS_UNALIGNED(t1)) {
                //    print("NON-CUDA ",dimi," ",dimj," ",IS_UNALIGNED(pc)," ",IS_UNALIGNED(t0)," ",IS_UNALIGNED(t1));
		    for (long j=0; j<nij; ++j) t0[j] = 0.0;
		  //  mTxm(dimi, dimj, dimj, t0, t.ptr(), pc);
                mTxmq(dimi, dimj, dimj, start_t0 + t0_off, start_tptr + tptr_off, start_pc + pc_off);

                //cu_mTxmq(dimi, dimj, dimj, t0, t.ptr(), pc,GPU_streams[i%NUM_STREAMS],1,t.size());
                //cu_mTxmq(dimi, dimj, dimj, start_t0 + t0_off, start_tptr + tptr_off, start_pc + pc_off,GPU_streams[i%NUM_STREAMS],1,0);
		    for (int n=1; n<t.ndim(); ++n) {
			for (long j=0; j<nij; ++j) t1[j] = 0.0;
		//	mTxm(dimi, dimj, dimj, t1, t0, pc);
			
                //cu_mTxmq(dimi, dimj, dimj, t1, t0, pc,GPU_streams[i%NUM_STREAMS],t.ndim(),t.size());
                //cu_mTxmq(dimi, dimj, dimj, start_t1 + t1_off, start_t0 + t0_off, start_pc + pc_off,GPU_streams[i%NUM_STREAMS],t.ndim(),0);
                mTxmq(dimi, dimj, dimj, start_t1 + t1_off, start_t0 + t0_off, start_pc + pc_off);
			std::swap(t0,t1);
		    }
		}
		else {
                */
		   //mTxmq(dimi, dimj, dimj, t0, t.ptr(), pc);
		    //print("CUDA KERNEL (dim = ",dimi,",",dimj,")\n");
		    //cu_mTxmq(dimi, dimj, dimj, t0, t.ptr(), pc,GPU_streams[i%NUM_STREAMS],1,t.size());
                    //cu_mTxmq(dimi, dimj, dimj, start_t0 + t0_off, start_tptr + tptr_off, start_pc + pc_off,GPU_streams[i%NUM_STREAMS],1,0);
                    cu_mTxmqq(op_dimi, op_dimj, op_dimj, start_t0 + /*t0_off*/i*tbuf_step, start_tptr + /*tptr_off*/i*tbuf_step, start_pc + /*pc_off*/i*dbuf_step,GPU_streams[i%NUM_STREAMS],0,0,cublas_handle);
                    /*
                    mTxmq(dimi, dimj, dimj, start_t0 + t0_off, start_tptr + tptr_off, start_pc + pc_off);
                    */
                    T* ptr_t0 = start_t0 + /*t0_off*/i*tbuf_step;
                    T* ptr_t1 = start_t1 + /*t1_off*/i*tbuf_step;
		    //for (int n=1; n<t.ndim(); ++n) {
			//mTxmq(dimi, dimj, dimj, t1, t0, pc);
			//cu_mTxmq(dimi, dimj, dimj, t1, t0, pc,GPU_streams[i%NUM_STREAMS],t.ndim(),t.size());
                        //cu_mTxmq(dimi, dimj, dimj, start_t1 + t1_off, start_t0 + t0_off, start_pc + pc_off,GPU_streams[i%NUM_STREAMS],t.ndim(),0);
                        cu_mTxmqq(op_dimi, op_dimj, op_dimj, /*start_t1 + t1_off*/ptr_t1, /*start_t0 + t0_off*/ptr_t0, start_pc + /*pc_off*/i*dbuf_step,GPU_streams[i%NUM_STREAMS],t.ndim(),1,cublas_handle);
                        //mTxmq(dimi, dimj, dimj, /*start_t1 + t1_off*/ptr_t1, /*start_t0 + t0_off*/ptr_t0, start_pc + pc_off);
			//std::swap(t0,t1);
			//std::swap(start_t0 + t0_off, start_t1 + t1_off);
                        /*
                        T* temp = ptr_t0;
                        ptr_t0 = ptr_t1;
                        ptr_t1 = temp;
                        */
                        //T* temp_buf = new T[dimi*dimj*sizeof(T)];
                        //memcpy(temp_buf, start_t0 + t0_off, dimi*dimj*sizeof(T));
                        //memcpy(start_t0 + t0_off, start_t1 + t1_off, dimi*dimj*sizeof(T));
                        //memcpy(start_t1 + t1_off, temp_buf, dimi*dimj*sizeof(T));
                        //delete[] temp_buf;
		    //}
		//}
                
                //t0_off += dimi*dimj;
                //tptr_off += dimi*dimj;
                //pc_off += dimj*dimj;
                //t1_off += dimi*dimj;
                //std::tr1::tuple<tensorT*,int,keyT,containerT,bool,keyT,dcT,tensorT*,Tensor<double>*> temp = inObj.at(i)->compressop_fasttransform(inArgs.at(i))/*compressop_compute(inArgs.at(i))*/;
                //outArg[i] = temp;
	   }

            //synchronize streams
            streams_synchronize(GPU_streams,NUM_STREAMS);
ENDt_TIMER("sTREAMS");

STARTt_TIMER;
            //if (d.ndim()&1)
              //CPUtransfer_buffer(t0_buf,start_t0,t0_off); 
            //else
              CPUtransfer_buffer(t1_buf,start_t1,/*t1_off*/tbuf_size); 
ENDt_TIMER("transfer2");
            
STARTt_TIMER;
            delete[] t0_buf;           
            GPUdelete_buffer(start_t0);
            delete[] tptr_buf;           
            GPUdelete_buffer(start_tptr);
            GPUdelete_buffer(start_pc);
            delete[] pc_buf;           
            GPUdelete_buffer(start_t1);           
            //delete[] t1_buf;           
ENDt_TIMER("cleanup"); 

            //t1_off = 0;
            //t0_off = 0;
            for (unsigned int i = 0; i < inArgs.size(); i++){
STARTt_TIMER;
		tensorT d = std::tr1::get<0>(inArgs.at(i));
		int k = std::tr1::get<1>(inArgs.at(i));
		FunctionCommonData<T, NDIM> cdata = FunctionCommonData<T, NDIM>::get(k);
		keyT key = std::tr1::get<2>(inArgs.at(i));
		containerT ct = std::tr1::get<3>(inArgs.at(i)); 
		bool nonstandard = std::tr1::get<4>(inArgs.at(i)); 
		keyT parent = std::tr1::get<5>(inArgs.at(i)); 
		dcT dc = std::tr1::get<6>(inArgs.at(i)); 
		long dimi = std::tr1::get<7>(inArgs.at(i));
		long dimj = std::tr1::get<8>(inArgs.at(i));
		Tensor<double>* c = std::tr1::get<9>(inArgs.at(i));
		T* t0 = std::tr1::get<10>(inArgs.at(i));
		T* t1 = std::tr1::get<11>(inArgs.at(i));
		tensorT* result = std::tr1::get<12>(inArgs.at(i));
		tensorT* workspace = std::tr1::get<13>(inArgs.at(i));
		//const double* pc = std::tr1::get<14>(inArgs.at(i));
ENDt_TIMER("final loop");

                //if (d.ndim()&1)
                  //memcpy(/*t0*/result->ptr(),t0_buf + t0_off,dimi*dimj*sizeof(T));
                //else
STARTt_TIMER;
                  memcpy(/*t1*/result->ptr(),t1_buf + /*t1_off*/i*tbuf_step,/*dimi*dimj*sizeof(T)*/tbuf_stepbytes);
ENDt_TIMER("memcpy3");
                
                std::tr1::tuple<tensorT*,int,keyT,containerT,bool,keyT,dcT,tensorT*,Tensor<double>*> t11(result,k,key,ct,nonstandard,parent,dc,workspace,c) ;
                outArg[i] = t11;
                //t1_off += dimi*dimj;
                //t0_off += dimi*dimj;
            }                 
            //delete[] t0_buf;           
            delete[] t1_buf;           
            return outArg;
        }

        std::vector< std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT> > compressop_allCompute(std::vector< std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT> > inArgs, std::vector< FunctionNode<T,NDIM>* > inObj){
            std::vector< std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT> > outArg(inArgs.size(),inObj.at(0)->compressop_compute(inArgs.at(0)));
            //print("inArgs.size() = ",inArgs.size());
            for (unsigned int i = 0; i < inArgs.size(); i++){
                std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT> temp = inObj.at(i)->compressop_compute(inArgs.at(i));
                outArg[i] = temp;
            }

            return outArg;
        }       
 
        Void compressop_postprocess(std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT> in){
          //  std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT> in = in1.get();
          tensorT d = std::tr1::get<0>(in);
          int k = std::tr1::get<1>(in);
          FunctionCommonData<T, NDIM> cdata = FunctionCommonData<T, NDIM>::get(k);
          keyT key = std::tr1::get<2>(in); 
          containerT ct = std::tr1::get<3>(in); 
          bool nonstandard = std::tr1::get<4>(in);
          keyT parent = std::tr1::get<5>(in);
          dcT dc = std::tr1::get<6>(in);
 
          if (this->has_coeff()) {
                const tensorT& c = this->coeff();
                if (c.dim(0) == k) {
                    d(cdata.s0) += c;
                }
                else {
                    d += c;
                }
            }

            tensorT s = copy(d(cdata.s0));

            if (key.level()> 0 && !nonstandard)
                d(cdata.s0) = 0.0;

            this->set_coeff(d);

            if (key.level() > 0){
              ct.update(parent, &TensorNode<T, NDIM>::update_map, dc, key, s, nonstandard, k);
            }

            return None;
        }
        
        Void compressop_postprocessGPU(std::tr1::tuple<tensorT*,int,keyT,containerT,bool,keyT,dcT,tensorT*,Tensor<double>*> in){
          //  std::tr1::tuple<tensorT,int,keyT,containerT,bool,keyT,dcT> in = in1.get();
          tensorT* result = std::tr1::get<0>(in);
          tensorT d(*result); delete result;
          int k = std::tr1::get<1>(in);
          FunctionCommonData<T, NDIM> cdata = FunctionCommonData<T, NDIM>::get(k);
          keyT key = std::tr1::get<2>(in); 
          containerT ct = std::tr1::get<3>(in); 
          bool nonstandard = std::tr1::get<4>(in);
          keyT parent = std::tr1::get<5>(in);
          dcT dc = std::tr1::get<6>(in);
          tensorT* workspace = std::tr1::get<7>(in); delete workspace;
          Tensor<double> * c1 = std::tr1::get<8>(in); delete c1; 
 
          if (this->has_coeff()) {
                const tensorT& c = this->coeff();
                if (c.dim(0) == k) {
                    d(cdata.s0) += c;
                }
                else {
                    d += c;
                }
            }

            tensorT s = copy(d(cdata.s0));

            if (key.level()> 0 && !nonstandard)
                d(cdata.s0) = 0.0;

            this->set_coeff(d);

            if (key.level() > 0){
              ct.update(parent, &TensorNode<T, NDIM>::update_map, dc, key, s, nonstandard, k);
            }

            return None;
        }

        Void compress_op(const keyT& key, dcT dc, const std::map<keyT, tensorT>& tensor_keys, containerT ct , const bool& nonstandard, const int& k, const keyT& parent) {

            if (HAVE_GPU & BACKTO_CPU)
                dc.local_updateGPU(key, &FunctionNode<T, NDIM>::compressop_preprocessGPU, &FunctionNode<T, NDIM>::compressop_allComputeGPU, &FunctionNode<T, NDIM>::compressop_postprocessGPU, &FunctionNode<T, NDIM>::compressop_backToCPU, std::tr1::tuple<keyT, dcT, std::map<keyT, tensorT>, containerT, bool, int, keyT>(key, dc, tensor_keys, ct, nonstandard, k, parent));
            else if (HAVE_GPU)
                dc.local_updateGPU(key, &FunctionNode<T, NDIM>::compressop_preprocessGPU, &FunctionNode<T, NDIM>::compressop_allComputeGPU, &FunctionNode<T, NDIM>::compressop_postprocessGPU, std::tr1::tuple<keyT, dcT, std::map<keyT, tensorT>, containerT, bool, int, keyT>(key, dc, tensor_keys, ct, nonstandard, k, parent));
            else if (SIM_GPU)
                dc.local_updateGPU(key, &FunctionNode<T, NDIM>::compressop_preprocess, &FunctionNode<T, NDIM>::compressop_allCompute, &FunctionNode<T, NDIM>::compressop_postprocess, std::tr1::tuple<keyT, dcT, std::map<keyT, tensorT>, containerT, bool, int, keyT>(key, dc, tensor_keys, ct, nonstandard, k, parent));
            else if (JUST_AGG)
                dc.local_updateJustAgg(key, &FunctionNode<T, NDIM>::compressop_preprocess, &FunctionNode<T, NDIM>::compressopComputePostprocess, std::tr1::tuple<keyT, dcT, std::map<keyT, tensorT>, containerT, bool, int, keyT>(key, dc, tensor_keys, ct, nonstandard, k, parent));
            else if (THREE_SPLIT)
                dc.local_update(key, &FunctionNode<T, NDIM>::compressop_preprocess, &FunctionNode<T, NDIM>::compressop_compute, &FunctionNode<T, NDIM>::compressop_postprocess, std::tr1::tuple<keyT, dcT, std::map<keyT, tensorT>, containerT, bool, int, keyT>(key, dc, tensor_keys, ct, nonstandard, k, parent));
            /*
            * RESTORE THIS
            std::tr1::tuple<tensorT,int,keyT, containerT,bool,keyT,dcT> t1 = compressop_preprocess(key, dc, tensor_keys, ct , nonstandard, k, parent);
            std::tr1::tuple<tensorT,int,keyT, containerT,bool,keyT,dcT > t2 = compressop_compute(t1);
            */

            //Future< std::tr1::tuple<tensorT,int,keyT, containerT,bool,keyT,dcT> > t2 = dc.task(key, &FunctionNode<T,NDIM>::compressop_compute, t1);  
            /*
            * RESTORE THIS
            dc.update(key, &FunctionNode<T,NDIM>::compressop_postprocess, t2);
            */
            else{
		    // Copy child scaling coeffs into contiguous block
		    FunctionCommonData<T,NDIM> cdata = FunctionCommonData<T,NDIM>::get(k);
		    std::map<keyT, tensorT>& tk = const_cast<std::map<keyT, tensorT>&>(tensor_keys);
		    tensorT d(cdata.v2k,false);
		    int i=0;
		    for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
			d(child_patch(kit.key(),cdata)) = tk[kit.key()];
		    }
		    d = filter(d,cdata);

		    if (this->has_coeff()) {
			const tensorT& c = this->coeff();
			if (c.dim(0) == k) {
			    d(cdata.s0) += c;
			}
			else {
			    d += c;
			}
		    }

		    tensorT s = copy(d(cdata.s0));

		    if (key.level()> 0 && !nonstandard)
			d(cdata.s0) = 0.0;

		    this->set_coeff(d); 

		    if (key.level() > 0){
		      ct.update(parent, &TensorNode<T, NDIM>::update_map, dc, key, s, nonstandard, k);
		    }
            }
            
            return None;

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
        Void accumulate(const Tensor<T>& t, const typename FunctionNode<T,NDIM>::dcT& c, const Key<NDIM>& key) {
            if (has_coeff()) {
                coeff() += t;
            }
            else {
                // No coeff and no children means the node is newly
                // created for this operation and therefore we must
                // tell its parent that it exists.
                coeff() = copy(t);
                if ((!_has_children) && key.level()> 0) {
                    Key<NDIM> parent = key.parent();
                    const_cast<dcT&>(c).task(parent, &FunctionNode<T,NDIM>::set_has_children_recursive, c, parent);
                }
            }
            return None;
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar & coeff() & _has_children & double(_norm_tree);
        }
    };

    template <typename T, std::size_t NDIM>
    std::ostream& operator<<(std::ostream& s, const FunctionNode<T,NDIM>& node) {
        s << "(has_coeff=" << node.has_coeff() << ", has_children=" << node.has_children() << ", norm=";
        double norm = node.has_coeff() ? node.coeff().normf() : 0.0;
        if (norm < 1e-12)
            norm = 0.0;
        s << norm << ")";
        return s;
    }


    template <typename T, std::size_t NDIM>
    class TensorNode{
        public:
          typedef Key<NDIM> keyT;
          typedef Tensor<T> tensorT;
          typedef FunctionNode<T,NDIM> nodeT; ///< Type of node
          typedef WorldContainer<keyT, TensorNode> containerT;
          typedef WorldContainer<keyT, nodeT> dcT;
          typedef std::map<keyT, tensorT> mapT;
 
          std::map<keyT, tensorT> tensor_keys;
          std::size_t max;
          bool set;
          keyT parent;

          Spinlock ready;

          TensorNode() : set(false){
            max = (std::size_t)pow(2, NDIM);
          }

          TensorNode(const keyT& p) : parent(p), set(false){
            max = (std::size_t)pow(2, NDIM);
          }

          TensorNode(const TensorNode<T, NDIM>& other) {
            *this = other;
          }

          TensorNode<T, NDIM>&
          operator=(const TensorNode<T, NDIM>& other) {
            if (this != &other) {
              this->parent = other.parent;
              this->set = other.set;
              this->max = (std::size_t)pow(2, NDIM);
              typename mapT::const_iterator it;
              for (it = other.tensor_keys.begin(); it != other.tensor_keys.end(); it++){
                tensor_keys.insert(std::pair<keyT, tensorT>(it->first,it->second));
              }
            }
            return *this;
          }

          Void compress_dc(const keyT& key, containerT ct, dcT dc, const keyT& p, bool nonstandard, int k){
            ready.lock();
              this->parent = p;
              this->set = true;
              if (tensor_keys.size() == max){
                //dc.update(key, &nodeT::compress_op, tensor_keys, ct, nonstandard, k, this->parent);
                if (HAVE_GPU)
                  dc.local_updateGPU(key, &nodeT::compressop_preprocessGPU, &nodeT::compressop_allComputeGPU, &nodeT::compressop_postprocessGPU, std::tr1::tuple< keyT, dcT, std::map<keyT, tensorT>, containerT, bool, int, keyT >(key, dc, tensor_keys, ct, nonstandard, k, this->parent));
                else if (SIM_GPU)
                  dc.local_updateGPU(key, &nodeT::compressop_preprocess, &nodeT::compressop_allCompute, &nodeT::compressop_postprocess, std::tr1::tuple< keyT, dcT, std::map<keyT, tensorT>, containerT, bool, int, keyT >(key, dc, tensor_keys, ct, nonstandard, k, this->parent));
                else if (JUST_AGG)
                  dc.local_updateJustAgg(key, &nodeT::compressop_preprocess, &nodeT::compressopComputePostprocess, std::tr1::tuple< keyT, dcT, std::map<keyT, tensorT>, containerT, bool, int, keyT >(key, dc, tensor_keys, ct, nonstandard, k, this->parent));
                else if (THREE_SPLIT)
                  dc.local_update(key, &nodeT::compressop_preprocess, &nodeT::compressop_compute, &nodeT::compressop_postprocess, std::tr1::tuple< keyT, dcT, std::map<keyT, tensorT>, containerT, bool, int, keyT >(key, dc, tensor_keys, ct, nonstandard, k, this->parent));
                else
                  dc.update(key, &nodeT::compress_op, tensor_keys, ct, nonstandard, k, this->parent);
              }
            ready.unlock();

            return None;
          }

          Void update_map(const keyT& key, containerT ct, dcT dc, keyT childKey, tensorT t, bool nonstandard, int k){
            
            ready.lock();
            this->tensor_keys.insert(std::pair<keyT, tensorT>(childKey,t));
            if (tensor_keys.size() == max && this->set){
              dc.update(key, &nodeT::compress_op, tensor_keys, ct, nonstandard, k, this->parent);
            }
            ready.unlock();

            return None;
          }

          template <typename Archive>
          void serialize(Archive& ar) {
            ar & tensor_keys & parent & max;
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
        typedef Tensor<T> tensorT; ///< Type of tensor used to hold coeffs
        typedef Vector<Translation,NDIM> tranT; ///< Type of array holding translation
        typedef Key<NDIM> keyT; ///< Type of key
        typedef FunctionNode<T,NDIM> nodeT; ///< Type of node
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

        const FunctionCommonData<T,NDIM>& cdata;

        std::shared_ptr< FunctionFunctorInterface<T,NDIM> > functor;

        bool compressed; ///< Compression status

        dcT coeffs; ///< The coefficients

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
                //, bc(factory._bc)
            {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            // !!! Ensure that all local state is correctly formed
            // before invoking process_pending for the coeffs and
            // for this.  Otherwise, there is a race condition.
            MADNESS_ASSERT(k>0 && k<=MAXK);

            bool empty = factory._empty;
            bool do_refine = factory._refine;

            if (do_refine)
                initial_level = std::max(0,initial_level - 1);

            if (empty) { // Do not set any coefficients at all
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
   
            pthread_barrier_init(&barr, NULL, NUMPAR_MEMCPYS);         
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
                , cdata(FunctionCommonData<T,NDIM>::get(k))
                , functor()
                , compressed(other.compressed)
                , coeffs(world, pmap ? pmap : other.coeffs.get_pmap())
                //, bc(other.bc)
        {
            if (dozero) {
                initial_level = 1;
                insert_zero_down_to_initial_level(cdata.key0);
            }
            
            pthread_barrier_init(&barr, NULL, NUMPAR_MEMCPYS);         
            coeffs.process_pending();
            this->process_pending();
        }

        virtual ~FunctionImpl() { }


    static void* blank(void * data){
        MemcpyCompressData<T,NDIM> * memcpy_data = (MemcpyCompressData<T,NDIM> *)data;
        pthread_barrier_wait(&barr);

        return NULL;
    }
    
    static void* memcpyCompress(void * data){
        MemcpyCompressData<T,NDIM> * memcpy_data = (MemcpyCompressData<T,NDIM> *)data;
        pthread_barrier_wait(&barr);

        memcpy(memcpy_data->tstart,memcpy_data->dptr,memcpy_data->tbuf_stepbytes);
        memcpy(memcpy_data->dstart,memcpy_data->pc,memcpy_data->dbuf_stepbytes);

        return NULL;
    }
        Void reconstruct_dc_tasknoupdate(const keyT& key, dcT dc, const tensorT& s);

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

        bool is_nonstandard() const {return nonstandard;}

        double get_thresh() const {return thresh;}

        void set_thresh(double value) {thresh = value;}

        bool get_autorefine() const {return autorefine;}

        void set_autorefine(bool value) {autorefine = value;}

        int get_k() const {return k;}

        const dcT& get_coeffs() const {return coeffs;}

        dcT& get_coeffs() {return coeffs;}

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
            if (world.rank() == coeffs.owner(cdata.key0))
                truncate_spawn(cdata.key0,tol);
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
        const tensorT parent_to_child(const tensorT& s, const keyT& parent, const keyT& child) const;


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
                    tensorT qq;

                    if (it == coeffs.end()) {
                        // must get from above
                        typedef std::pair< keyT,Tensor<T> > pairT;
                        Future<pairT> result;
                        sock_it_to_me(key, result.remote_ref(world));
                        const keyT& parent = result.get().first;
                        const tensorT& t = result.get().second;

                        qq = parent_to_child(t, parent, key);
                    } else {
                        qq = it->second.coeff();
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
                    r(s) = qq(cdata.s0);

                }
            }

            world.gop.fence();
            world.gop.sum(r0);
            //print(r,r0);

            return r0;
        }

        template <typename Q>
        Tensor<Q> coeffs2values(const keyT& key, const Tensor<Q>& coeff) const {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            double scale = pow(2.0,0.5*NDIM*key.level())/sqrt(FunctionDefaults<NDIM>::get_cell_volume());
            return transform(coeff,cdata.quad_phit).scale(scale);
        }

        template <typename Q>
        Tensor<Q> values2coeffs(const keyT& key, const Tensor<Q>& values) const {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            double scale = pow(0.5,0.5*NDIM*key.level())*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
            return transform(values,cdata.quad_phiw).scale(scale);
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
            coeffs.replace(key, nodeT(tcube,false));
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
            coeffs.replace(key, nodeT(tcube,false));
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
                    op(parent, node.coeff());
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
                    tensorT& t= node.coeff();
                    //double before = t.normf();
                    tensorT values = impl->fcube_for_mul(key, key, t);
                    op(key, values);
                    double scale = pow(0.5,0.5*NDIM*key.level())*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
                    t = transform(values,impl->cdata.quad_phiw).scale(scale);
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
                    const Tensor<R>& r = it->second.coeff();
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
                                node.set_coeff(tensorT(cdata.v2k));
                            tensorT& t = node.coeff();
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
                    acc[i]->second = nodeT(c[i],false);
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
                        s(cdata.s0) = acc[i]->second.coeff()(___);
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
                c[i] = coeffs2values(key, v[i]->coeffs.find(key).get()->second.coeff()); // !!!!! gack
            }
            tensorT r = op(key, c);
            coeffs.replace(key, nodeT(values2coeffs(key, r),false));
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
                    coeffs.replace(key, nodeT(tensorT(),true));
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
                    lc = it->second.coeff();
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
                        rc = it->second.coeff();
                }
                else {
                    rnorm = rc.normf();
                }

                if (rc.size() && lc.size()) { // Yipee!
                    result->task(world.rank(), &implT:: template do_mul<L,R>, key, lc, std::make_pair(key,rc));
                }
                else if (tol && lnorm*rnorm < truncate_tol(tol, key)) {
                    result->coeffs.replace(key, nodeT(tensorT(cdata.vk),false)); // Zero leaf
                }
                else {
                    result->coeffs.replace(key, nodeT(tensorT(),true)); // Interior node
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
                    lc = it->second.coeff();
            }

            Tensor<R> rc = rcin;
            if (rc.size() == 0) {
                riterT it = right->coeffs.find(key).get();
                MADNESS_ASSERT(it != right->coeffs.end());
                rnorm = it->second.get_norm_tree();
                if (it->second.has_coeff())
                    rc = it->second.coeff();
            }

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
                    coeffs.replace(key, nodeT(tensorT(cdata.vk),false)); // Zero leaf node
                    return None;
                }
            }

            // Recur down
            coeffs.replace(key, nodeT(tensorT(),true)); // Interior node

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
                    lc = it->second.coeff();
            }

            Tensor<R> rc = rcin;
            if (rc.size() == 0) {
                riterT it = right->coeffs.find(key).get();
                MADNESS_ASSERT(it != right->coeffs.end());
                if (it->second.has_coeff())
                    rc = it->second.coeff();
            }

            if (rc.size() && lc.size()) { // Yipee!
                do_binary_op<L,R>(key, lc, std::make_pair(key,rc), op);
                return None;
            }

            // Recur down
            coeffs.replace(key, nodeT(tensorT(),true)); // Interior node

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

            const Tensor<Q>& fc = func->coeffs.find(key).get()->second.coeff();

            if (fc.size() == 0) {
                // Recur down
                coeffs.replace(key, nodeT(tensorT(),true)); // Interior node
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT& child = kit.key();
                    woT::task(coeffs.owner(child), &implT:: template unaryXXa<Q,opT>, child, func, op);
                }
            }
            else {
                coeffs.replace(key, nodeT(op(key, fc),false)); // Leaf node
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
                           const RemoteReference< FutureImpl< std::pair<keyT,tensorT> > >& ref) const;
        /// As above, except
        /// 3) The coeffs are constructed from the avg of nodes further down the tree
        Void sock_it_to_me_too(const keyT& key,
                               const RemoteReference< FutureImpl< std::pair<keyT,tensorT> > >& ref) const;

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
        bool noautorefine(const keyT& key, const tensorT& t) const {
            return false;
        }

        /// Returns true if this block of coeffs needs autorefining
        bool autorefine_square_test(const keyT& key, const tensorT& t) const {
            double lo, hi;
            tnorm(t, &lo, &hi);
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

        Void sum_down_spawn(const keyT& key, const tensorT& s) {
            typename dcT::accessor acc;
            coeffs.insert(acc,key);
            nodeT& node = acc->second;
            tensorT& c = node.coeff();

            //print(key,"received",s.normf(),c.normf(),node.has_children());

            if (s.size() > 0) {
                if (c.size() > 0)
                    c.gaxpy(1.0,s,1.0);
                else
                    c = s;
            }

            if (node.has_children()) {
                tensorT d;
                if (c.size() > 0) {
                    d = tensorT(cdata.v2k);
                    d(cdata.s0) = c;
                    d = unfilter(d);
                    node.clear_coeff();
                }
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    tensorT ss;
                    const keyT& child = kit.key();
                    if (d.size() > 0) ss = copy(d(child_patch(child)));
                    //print(key,"sending",ss.normf(),"to",child);
                    woT::task(coeffs.owner(child), &implT::sum_down_spawn, child, ss);
                }
            }
            else {
                // Missing coeffs assumed to be zero
                if (c.size() <= 0) c = tensorT(cdata.vk);
            }
            return None;
        }

        /// After 1d push operator must sum coeffs down the tree to restore correct scaling function coefficients
        void sum_down(bool fence) {
            if (world.rank() == coeffs.owner(cdata.key0)) sum_down_spawn(cdata.key0, tensorT());

            if (fence) world.gop.fence();
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
                                coeffs.task(dest, &nodeT::accumulate, result, coeffs, dest, TaskAttributes::hipri());
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
            for (fiterT it=f->coeffs.begin(); it!=end; ++it) {
                const fnodeT& node = it->second;
                if (node.has_coeff()) {
                    const keyT& key = it->first;
                    const Tensor<R>& c = node.coeff();
                    world.taskq.add(*this, &implT:: template apply_1d_realspace_push_op<opT,R>, archive::archive_ptr<const opT>(&op), axis, key, c);
                }
            }
            if (fence) world.gop.fence();
        }

        Void do_diff1(const DerivativeBase<T,NDIM>* D,
                      const implT* f,
                      const keyT& key,
                      const std::pair<keyT,tensorT>& left,
                      const std::pair<keyT,tensorT>& center,
                      const std::pair<keyT,tensorT>& right) {
            return D->do_diff1(f,this,key,left,center,right);
        }


        // Called by result function to differentiate f
        void diff(const DerivativeBase<T,NDIM>* D, const implT* f, bool fence) {
            typedef std::pair<keyT,tensorT> argT;
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
                    coeffs.replace(key,nodeT(tensorT(),true)); // Empty internal node
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
        Future< std::pair<keyT,tensorT> > find_me(const keyT& key) const;


        /// Permute the dimensions according to map
        void mapdim(const implT& f, const std::vector<long>& map, bool fence);

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

        /// Projects old function into new basis (only in reconstructed form)
        void project(const implT& old, bool fence) {
            long kmin = std::min(cdata.k,old.cdata.k);
            std::vector<Slice> s(NDIM,Slice(0,kmin-1));
            typename dcT::const_iterator end = old.coeffs.end();
            for (typename dcT::const_iterator it=old.coeffs.begin(); it!=end; ++it) {
                const keyT& key = it->first;
                const nodeT& node = it->second;
                if (node.has_coeff()) {
                    tensorT c(cdata.vk);
                    c(s) = node.coeff()(s);
                    coeffs.replace(key,nodeT(c,false));
                }
                else {
                    coeffs.replace(key,nodeT(tensorT(),true));
                }
            }
            if (fence)
                world.gop.fence();
        }

        struct true_refine_test {
            bool operator()(const implT* f, const keyT& key, const tensorT& t) const {
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
            if (node.has_coeff() && key.level() < max_refine_level && op(this, key, node.coeff())) {
                tensorT d(cdata.v2k);
                d(cdata.s0) = node.coeff();
                d = unfilter(d);
                node.clear_coeff();
                node.set_has_children(true);
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT& child = kit.key();
                    tensorT ss = copy(d(child_patch(child)));
                    coeffs.replace(child,nodeT(ss,-1.0,false)); // Note value -1.0 for norm tree to indicate result of refinement
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
                    world.taskq.add(*this, &implT::broaden_op, key, v);
                }
            }
            // Reset value of norm tree so that can repeat broadening
            if (fence) {
                world.gop.fence();
                zero_norm_tree();
                world.gop.fence();
            }
        }

        void reconstruct(bool fence) {
            // Must set true here so that successive calls without fence do the right thing
            //world.gop.fence();
            nonstandard = compressed = false;
            if (world.rank() == coeffs.owner(cdata.key0)){
                //world.taskq.add(*this, &implT::reconstruct_op, cdata.key0,tensorT());
                //world.taskq.add(*this, &implT::reconstruct_update, cdata.key0,tensorT());
                //print("rec_dc");
                //coeffs.update(cdata.key0, &FunctionNode<T,NDIM>::reconstruct_dc_task, tensorT(), cdata.k);
                world.taskq.add(*this,&FunctionImpl<T,NDIM>::reconstruct_dc_tasknoupdate, cdata.key0, coeffs, tensorT());
            }
                //reconstruct_update(cdata.key0,tensorT());
            if (fence)
                world.gop.fence();
        }

        // Invoked on node where key is local
        //Void reconstruct_op(const keyT& key, const tensorT& s);
    //template <typename T, std::size_t NDIM>
    inline Void /*FunctionImpl<T,NDIM>::*/reconstruct_update(const keyT& key, const tensorT& s){
        // Note that after application of an integral operator not all
        // siblings may be present so it is necessary to check existence
        // and if absent insert an empty leaf node.
        //
        // If summing the result of an integral operator (i.e., from
        // non-standard form) there will be significant scaling function
        // coefficients at all levels and possibly difference coefficients
        // in leaves, hence the tree may refine as a result.

       //print("CPS\n");

        //reconstruct_access_data(key, s);

        //return None;
      //}   

    //template <typename T, std::size_t NDIM>
    //inline Void /*FunctionImpl<T,NDIM>::*/reconstruct_access_data(const keyT& key, const tensorT& s){
        //woT::task(coeffs.owner(key), &implT::reconstruct_prepare_work, key, s);
        //return None;
    //}

    //template <typename T, std::size_t NDIM>
    //inline Void /*FunctionImpl<T,NDIM>::*/reconstruct_prepare_work(const keyT& key, const tensorT& s){

        typename dcT::iterator it = coeffs.find(key).get();
            nodeT node1;
            //bool inserted = false;

	    if (it == coeffs.end()) {
	        node1 = nodeT(tensorT(),false);
                //inserted = true;
            }
            else node1 = it->second;

            nodeT& node = node1;

	    // The integral operator will correctly connect interior nodes
	    // to children but may leave interior nodes without coefficients
	    // ... but they still need to sum down so just give them zeros
	    if (node.has_children() && !node.has_coeff()) {
		node.set_coeff(tensorT(cdata.v2k));
	    }

	    if (node.has_children() || node.has_coeff()) { // Must allow for inconsistent state from transform, etc.
		tensorT d = node.coeff();
		if (d.size() == 0) d = tensorT(cdata.v2k);
		if (key.level() > 0) d(cdata.s0) += s; // -- note accumulate for NS summation
		d = unfilter(d);
		node.clear_coeff();
		node.set_has_children(true);
		for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
		    const keyT& child = kit.key();
		    tensorT ss = copy(d(child_patch(child)));
                    world.taskq.add(*this, &implT::reconstruct_update, child, ss);
		}
	    }
	    else {
		if (key.level()) node.set_coeff(copy(s));
		else node.set_coeff(s);
	    }
	    //std::pair<nodeT, bool> out(node, true);
        //std::pair<nodeT, bool> out = reconstruct_do_work(key, it, s);
        //if (out.second){
          coeffs.replace(key, node);
        //}
        return None;
    }
        inline std::pair<nodeT, bool> reconstruct_do_work(const keyT& key, const typename dcT::iterator it, const tensorT& s){

            nodeT node;

	    if (it == coeffs.end()) {
	        node = nodeT(tensorT(),false);
            }
            else node = it->second;

	    // The integral operator will correctly connect interior nodes
	    // to children but may leave interior nodes without coefficients
	    // ... but they still need to sum down so just give them zeros
	    if (node.has_children() && !node.has_coeff()) {
		node.set_coeff(tensorT(cdata.v2k));
	    }

	    if (node.has_children() || node.has_coeff()) { // Must allow for inconsistent state from transform, etc.
		tensorT d = node.coeff();
		if (d.size() == 0) d = tensorT(cdata.v2k);
		if (key.level() > 0) d(cdata.s0) += s; // -- note accumulate for NS summation
		d = unfilter(d);
		node.clear_coeff();
		node.set_has_children(true);
		for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
		    const keyT& child = kit.key();
		    tensorT ss = copy(d(child_patch(child)));
                    world.taskq.add(*this, &implT::reconstruct_update, child, ss);
		}
	    }
	    else {
		if (key.level()) node.set_coeff(copy(s));
		else node.set_coeff(s);
	    }
	    return std::pair<nodeT, bool>(node, true);
	}
 
        //inline Void reconstruct_update(const keyT& key, const tensorT& s);
        //inline Void reconstruct_access_data(const keyT& key, const tensorT& s);
        //inline Void reconstruct_prepare_work(const keyT& key, const tensorT& s);
        Void reconstruct_op(const keyT& key, const tensorT& s);
	
        typedef WorldContainer<keyT, TensorNode<T, NDIM>> tensorContainerT;

        void compress(bool nonstandard, bool keepleaves, bool fence) {
            // Must set true here so that successive calls without fence do the right thing
            this->compressed = true;
            this->nonstandard = nonstandard;
            tensorContainerT tensorTree(world,coeffs.get_pmap());
            world.gop.fence();
            if (world.rank() == coeffs.owner(cdata.key0)){
                //compress_spawn(cdata.key0, nonstandard, keepleaves);
                print("com_dc");
                coeffs.update(cdata.key0, &nodeT::top_down, tensorTree, nonstandard, keepleaves, cdata.key0, cdata.k);
            }
            if (fence)
                world.gop.fence();
        }

        // Invoked on node where key is local
        Future<tensorT> compress_spawn(const keyT& key, bool nonstandard, bool keepleaves);

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
                return world.taskq.add(*this,&implT::norm_tree_op, key, v);
            }
            else {
                return Future<double>(node.coeff().normf());
            }
        }

        tensorT compress_op(const keyT& key, const std::vector< Future<tensorT> >& v, bool nonstandard) {
            PROFILE_MEMBER_FUNC(FunctionImpl);
            // Copy child scaling coeffs into contiguous block
            tensorT d(cdata.v2k,false);
            int i=0;
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
                d(child_patch(kit.key())) = v[i].get();
            }
            d = filter(d);

            typename dcT::accessor acc;
            MADNESS_ASSERT(coeffs.find(acc, key));

            if (acc->second.has_coeff()) {
                const tensorT& c = acc->second.coeff();
                if (c.dim(0) == k) {
                    d(cdata.s0) += c;
                }
                else {
                    d += c;
                }
            }

            tensorT s = copy(d(cdata.s0));

            if (key.level()> 0 && !nonstandard)
                d(cdata.s0) = 0.0;

            acc->second.set_coeff(d);

            return s;
        }

        /// Changes non-standard compressed form to standard compressed form
        void standard(bool fence) {
            typename dcT::iterator end = coeffs.end();
            for (typename dcT::iterator it=coeffs.begin(); it!=end; ++it) {
                const keyT& key = it->first;
                nodeT& node = it->second;
                if (key.level()> 0 && node.has_coeff()) {
                    if (node.has_children()) {
                        // Zero out scaling coeffs
                        node.coeff()(cdata.s0) = 0.0;
                    }
                    else {
                        // Deleting both scaling and wavelet coeffs
                        node.clear_coeff();
                    }
                }
            }
            nonstandard = false;
            if (fence)
                world.gop.fence();
        }

        struct do_op_args {
            keyT key, d, dest;
            double tol, fac, cnorm;
            do_op_args() {}
            do_op_args(const keyT& key, const keyT& d, const keyT& dest, double tol, double fac, double cnorm)
                    : key(key), d(d), dest(dest), tol(tol), fac(fac), cnorm(cnorm) {}
            template <class Archive>
            void serialize(Archive& ar) {
                ar & archive::wrap_opaque(this,1);
            }
        };

        template <typename opT, typename R>
        Void do_apply_kernel(const opT* op, const Tensor<R>& c, const do_op_args& args) {
            tensorT result = op->apply(args.key, args.d, c, args.tol/args.fac/args.cnorm);

            //print("APPLY", key, d, opnorm, cnorm, result.normf());

            // Screen here to reduce communication cost of negligible data
            // and also to ensure we don't needlessly widen the tree when
            // applying the operator
            if (result.normf()> 0.3*args.tol/args.fac) {
                // OPTIMIZATION NEEDED HERE ... CHANGING THIS TO TASK NOT SEND REMOVED
                // BUILTIN OPTIMIZATION TO SHORTCIRCUIT MSG IF DATA IS LOCAL
                coeffs.task(args.dest, &nodeT::accumulate, result, coeffs, args.dest, TaskAttributes::hipri());
            }
            return None;
        }

        template <typename opT, typename R>
        Void do_apply(const opT* op, const FunctionImpl<R,NDIM>* f, const keyT& key, const Tensor<R>& c) {
            PROFILE_MEMBER_FUNC(FunctionImpl);
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
                    double opnorm = op->norm(key.level(), d);
                    // working assumption here is that the operator is isotropic and
                    // montonically decreasing with distance
                    double tol = truncate_tol(thresh, key);

                    //print("APP", key, dest, cnorm, opnorm, (cnorm*opnorm> tol/fac));

                    if (cnorm*opnorm> tol/fac) {

                        // Most expensive part is the kernel ... do it in a separate task
                        if (d.distsq()==0) {
                            // This introduces finer grain parallelism
                            ProcessID where = world.rank();
                            do_op_args args(key, d, dest, tol, fac, cnorm);
                            woT::task(where, &implT:: template do_apply_kernel<opT,R>, op, c, args);
                        } else {
                            tensorT result = op->apply(key, d, c, tol/fac/cnorm);
                            if (result.normf()> 0.3*tol/fac) {
                                coeffs.task(dest, &nodeT::accumulate, result, coeffs, dest, TaskAttributes::hipri());
                            }
                        }
                    } else if (d.distsq() >= 1)
                        break; // Assumes monotonic decay beyond nearest neighbor
                }
            }
            return None;
        }

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
                        woT::task(p, &implT:: template do_apply<opT,R>, &op, &f, key, node.coeff());
                    }
                }
            }
            if (fence)
                world.gop.fence();
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
            const tensorT& coeff = node.coeff();
            ITERATOR(coeff,fval(IND)-=coeff(IND););

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
                return a+b;
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
        template <typename R>
        TENSOR_RESULT_TYPE(T,R) inner_local(const FunctionImpl<R,NDIM>& g) const {
            TENSOR_RESULT_TYPE(T,R) sum = 0.0;
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
                            sum += fnode.coeff().trace_conj(gnode.coeff());
                        }
                    }
                }
            }
            return sum;
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
            typename fdcT::const_iterator end = f.coeffs.end();
            for (typename fdcT::const_iterator it=f.coeffs.begin(); it!=end; ++it) {
                const keyT& key = it->first;
                const fnodeT& node = it->second;
                if (node.has_coeff()) {
                    coeffs.replace(key,nodeT(node.coeff()*q,node.has_children()));
                }
                else {
                    coeffs.replace(key,nodeT(tensorT(),node.has_children()));
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
                uniqueidT id;
                ar & id;
                World* world = World::world_from_id(id.get_world_id());
                MADNESS_ASSERT(world);
                ptr = static_cast< const FunctionImpl<T,NDIM>*>(world->ptr_from_id< WorldObject< FunctionImpl<T,NDIM> > >(id));
                if (!ptr)
                    MADNESS_EXCEPTION("FunctionImpl: remote operation attempting to use a locally uninitialized object",0);
            }
        };

        template <class Archive, class T, std::size_t NDIM>
        struct ArchiveStoreImpl<Archive,const FunctionImpl<T,NDIM>*> {
            static void store(const Archive& ar, const FunctionImpl<T,NDIM>*const& ptr) {
                ar & ptr->id();
            }
        };

        template <class Archive, class T, std::size_t NDIM>
        struct ArchiveLoadImpl<Archive, FunctionImpl<T,NDIM>*> {
            static void load(const Archive& ar, FunctionImpl<T,NDIM>*& ptr) {
                uniqueidT id;
                ar & id;
                World* world = World::world_from_id(id.get_world_id());
                MADNESS_ASSERT(world);
                ptr = static_cast< FunctionImpl<T,NDIM>*>(world->ptr_from_id< WorldObject< FunctionImpl<T,NDIM> > >(id));
                if (!ptr)
                    MADNESS_EXCEPTION("FunctionImpl: remote operation attempting to use a locally uninitialized object",0);
            }
        };

        template <class Archive, class T, std::size_t NDIM>
        struct ArchiveStoreImpl<Archive, FunctionImpl<T,NDIM>*> {
            static void store(const Archive& ar, FunctionImpl<T,NDIM>*const& ptr) {
                ar & ptr->id();
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
