#include <complex>
#include <octtree/octtree.h>
#include <mra/mra.h>
#include <mra/twoscale.h>
#include <mra/legendre.h>
#include <tasks/sav.h>
#include <tasks/tasks.h>
#include <misc/print.h>


/// \file mra/taskstuff.cc
/// \brief Implements Function stuff that uses tasks


namespace madness {

    /// A task tailored to the needs of compress
    template <typename T>
    class TaskCompress : public TaskInterface {
    public:
        typedef Tensor<T> TensorT;
        typedef SAV<TensorT> ArgT;
        typedef OctTree<FunctionNode> OctTreeT;
        typedef Function<T> FunctionT;
    private:
        FunctionT* f;
        OctTreeTPtr tree;
        ArgT outarg;
        ArgT inarg[2][2][2];
    public:
        TaskCompress(FunctionT* f, OctTreeTPtr tree, ArgT in[2][2][2], ArgT& out)
                : f(f)
                , tree(tree)
                , outarg(out) 
        {
            FORIJK(inarg[i][j][k] = in[i][j][k];);
        };

        bool probe() const {
            FOREACH_CHILD(OctTreeTPtr, tree,
                          if (f->isactive(child) && !inarg[i][j][k].probe())
                          return false;);
            return true;
        };

        void run() {
            f->_compressop(tree,inarg,outarg);
        };

        virtual ~TaskCompress() {};
    };
    
    /// Hash from function+node to a "unique" message tag 
    static inline int taghash(int ind, Level n, Translation x, Translation y, Translation z) {
        // Presently assume all 31 bits in tag can be used, but in reality MPI::TAG_UB is
        // the upper bound and some implementations seem to have a limit of 32k.
        // Also, the first 2048=2^11 tags are reserved.
        //
        // Must have all of ind included since operations on difference functions
        // may be interleaved.  However, operations for a given function are much 
        // more likely to be ordered and less likely to collide.  Less likely is 
        // not the same as never.  The only way to improve upon this is by registering
        // a unique tag for every remote child.  This is doable and we should switch to this.
        MADNESS_ASSERT(ind<2048);
        int n3 = n&7;
        int x2 = x&3;
        int y2 = y&3;
        int z2 = z&3;
        int tag = (ind<<9)|(n3<<6)|(x2<<4)|(y2<<2)|(z2);
        tag = tag << 11;
        return tag; 
    };

    /// Called by consumer to make a variable that will be set by producer
    template <typename T>
    SAV< Tensor<T> > Function<T>::input_arg(const OctTreeTPtr& consumer, const OctTreeTPtr& producer) {
        int tag = taghash(ind,producer->n(),producer->x(),producer->y(),producer->z());
        if (consumer->islocal() && producer->islocal())
            return SAV< Tensor<T> >();
        else if (consumer->islocal() && producer->isremote())
            return SAV< Tensor<T> >(producer->rank(), tag, true, this->data->cdata->vk);
        else if (consumer->isremote() && producer->islocal())
            return SAV< Tensor<T> >(consumer->rank(), tag, false, this->data->cdata->vk);
        else
            throw "input_arg: should not happen?";
    };

    /// Called by consumer to make a variable that will be set by producer
    template <typename T>
    SAV<bool> Function<T>::input_argb(const OctTreeTPtr& consumer, const OctTreeTPtr& producer) {
        int tag = taghash(ind,producer->n(),producer->x(),producer->y(),producer->z());
        if (consumer->islocal() && producer->islocal())
            return SAV<bool>();
        else if (consumer->islocal() && producer->isremote())
            return SAV<bool>(producer->rank(), tag, true);
        else if (consumer->isremote() && producer->islocal())
            return SAV<bool>(consumer->rank(), tag, false);
        else
            throw "input_argb: should not happen?";
    };


    /// Compress function (scaling function to wavelet)

    /// Communication streams up the tree.
    /// Returns self for chaining.
    template <typename T>
    Function<T>& Function<T>::compress(bool nonstandard) {
        if (data->compressed && data->nonstandard!=nonstandard) {
            madness::print("!!!! reconstructing before NS compress");
            reconstruct();
        }
        if (!data->compressed) {
            data->nonstandard = nonstandard;
            if (isactive(tree())) {
                SAV<TensorT> dummy;
                _compress(tree(),dummy);
            }
            taskq.global_fence();
            data->compressed = true;
        }
        return *this;
    };

    template <typename T>
    void Function<T>::_compress(OctTreeTPtr& tree, SAV<TensorT>& parent) {
        SAV<TensorT> args[2][2][2];
        FOREACH_ACTIVE_CHILD(OctTreeTPtr, tree,
                             args[i][j][k] = this->input_arg(tree,child);
                             if (child->islocal()) _compress(child,args[i][j][k]););

        if (tree->islocal()) {
            taskq.add_local(new TaskCompress<T>(this,tree,args,parent));
        }
    };

    template <typename T>
    void Function<T>::_compressop(OctTreeTPtr& tree, SAV<TensorT> args[2][2][2], SAV<TensorT>& parent) {
        Slice* s = data->cdata->s;      
        TensorT t = TensorT(2*k,2*k,2*k);
        int nchild=0;
        FOREACH_ACTIVE_CHILD(OctTreeTPtr, tree,
                             nchild++; 
                             t(s[i],s[j],s[k])=args[i][j][k].get(););
        if (nchild) {
            filter_inplace(t);
            if (coeff(tree)) t(s[0],s[0],s[0]) += *coeff(tree);
            parent.set(madness::copy(t(data->cdata->s0)));
            if (tree->n() > 0) t(data->cdata->s0)=0.0;
            set_coeff(tree,t);
        }
        else {
            parent.set(*coeff(tree));
            if (!data->nonstandard) unset_coeff(tree);
        }           
    };
    
    template <typename T>
    Function<T>& Function<T>::nsclean(bool leave_compressed) {
        MADNESS_ASSERT(iscompressed() && data->nonstandard);
        _nsclean(tree(),leave_compressed);
        data->compressed = leave_compressed;      
        return *this;      
    };
    
    template <typename T>
    void Function<T>::_nsclean(OctTreeTPtr& tree, bool leave_compressed) {
        if (isactive(tree) && islocal(tree)) {
            MADNESS_ASSERT(coeff(tree));
            Tensor<T>& t = *coeff(tree);
            if ( (leave_compressed && t.dim[0] == k) || 
                (!leave_compressed && t.dim[0] != k) ) unset_coeff(tree);
            FOREACH_ACTIVE_CHILD(OctTreeTPtr, tree, _nsclean(child,leave_compressed););
        }
    };
    
    /// A task tailored to the needs of truncate which is similar to compress
    template <typename T>
    class TaskTruncate : public TaskInterface {
    public:
        typedef SAV<bool> ArgT;
        typedef OctTree<FunctionNode> OctTreeT;
        typedef Function<T> FunctionT;
    private:
        FunctionT* f;
        OctTreeTPtr tree;
        ArgT outarg;
        ArgT inarg[2][2][2];
    public:
        TaskTruncate(FunctionT* f, OctTreeTPtr tree, ArgT in[2][2][2], ArgT& out)
                : f(f)
                , tree(tree)
                , outarg(out) 
        {
            FORIJK(inarg[i][j][k] = in[i][j][k];);
        };

        bool probe() const {
            FOREACH_CHILD(OctTreeTPtr, tree,
                          if (f->isactive(child) && !inarg[i][j][k].probe())
                          return false;);
            return true;
        };

        void run() {
            f->_truncateop(tree,inarg,outarg);
        };

        virtual ~TaskTruncate() {};
    };
    
    template <typename T>
    void Function<T>::_truncate(OctTreeTPtr& tree, SAV<bool>& parent) {
        SAV<bool> args[2][2][2];
        FOREACH_ACTIVE_CHILD(OctTreeTPtr, tree,
                             args[i][j][k] = this->input_argb(tree,child);
                             if (child->islocal()) _truncate(child,args[i][j][k]););

        if (tree->islocal()) {
            taskq.add_local(new TaskTruncate<T>(this,tree,args,parent));
        }
    };
    
    /// Inplace truncation

    /// Communication streams up the tree.  Returns self for chaining.
    template <typename T>
    Function<T>& Function<T>::truncate(double tol) {
        compress();
        if (tol <= 0.0) tol = this->data->thresh;
        this->data->truncate_thr = tol;
        if (isactive(tree())) {
            SAV<bool> dummy;
            _truncate(tree(),dummy);
        }
        taskq.global_fence();
        return *this;
    };
    
    
    template <typename T>
    void Function<T>::set_inactive_handler(Communicator& comm, ProcessID src, const AMArg& arg) {
        int ind = arg.arg0;
        Level n = arg.arg1;
        Translation x=arg.arg2, y=arg.arg3, z=arg.arg4;
        madness::print("set_inactive_handler",ind,n,x,y,z);
        Function<T> f = Function<T>(ind);
        OctTreeTPtr t = f.tree();
        MADNESS_ASSERT(t);
        if (t->n()==n && t->x()==x && t->y()==y && t->z()==z) {
            f.do_set_children_inactive(t);
        }
        else {
            MADNESS_EXCEPTION("set_inactive_handler: confused",0);
        }
    };
    
    
    template <typename T>   
    void Function<T>::do_set_children_inactive(OctTreeTPtr& tree) {
        FOREACH_ACTIVE_CHILD(OctTreeTPtr, tree, if (child->islocal()) set_inactive(child););
        // Gotta let those pesky away-from-home kids know also
        ProcessID ranks[8];
        int np = unique_active_child_procs(tree,ranks);
        if (np) {
            AMArg arg(ind, tree->n(), tree->x(), tree->y(), tree->z());
            for (int i=0; i<np; i++) {
                comm()->am_send(ranks[i], set_inactive_handler, arg);
            }
        }
    }
    
    template <typename T>
    void Function<T>::_truncateop(OctTreeTPtr& tree, SAV<bool> args[2][2][2], SAV<bool>& parent) {
        MADNESS_ASSERT(tree->islocal());
        int nchild=0; // Counts no. of active children with coefficients
        FOREACH_ACTIVE_CHILD(OctTreeTPtr, tree, if (args[i][j][k].get()) nchild++;);
        if (coeff(tree)  &&  nchild == 0  &&  tree->n()>0) {
            const Tensor<T>& t = *coeff(tree);
            double dnorm = t.normf();
            if (dnorm < truncate_tol(data->truncate_thr,tree->n())) {
                unset_coeff(tree);  // Note, we leave it active.
                do_set_children_inactive(tree);
            }
        }
        parent.set(coeff(tree));
    };   
    
    
    
    template <typename T>
    void Function<T>::_dodiff_kernel(Function<T>& df, OctTreeTPtr& tree, int axis,
                    const TensorT& t0, const TensorT& tm, const TensorT& tp) {
        // This is the slow but correct version.
        // There is a faster version that can be copied in from mad++-tree 
        // once this pile of s*** is working for sure.   
        const TensorT& rp_right = this->data->cdata->rp_right;
        const TensorT& rp_left = this->data->cdata->rp_left;
        const TensorT& rm_right = this->data->cdata->rm_right;
        const TensorT& rm_left = this->data->cdata->rm_left;
        const TensorT& r0 = this->data->cdata->r0;             
        
        const TensorT q0 = t0.swapdim(axis,0);
        const TensorT qp = tp.swapdim(axis,0);
        const TensorT qm = tm.swapdim(axis,0);
            
        Tensor<T> d = ::inner(r0, q0) + outer(rp_left,::inner(rp_right,qm)) + outer(rm_left,::inner(rm_right,qp));
        if (axis) d = ::copy(d.swapdim(axis,0));
        d.scale((double) two_to_power(tree->n()));
        df.set_coeff(tree,d);
    };
            

    template class Function<double>;
    template class Function< std::complex<double> >;
    template class FunctionData<double>;
    template class FunctionData< std::complex<double> >;
    template class FunctionCommonData<double>;
    template class TaskCompress<double>;
    template class TaskCompress< std::complex<double> >;    

}
