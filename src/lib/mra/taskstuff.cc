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

    static inline int taghash(Level n, Translation x, Translation y, Translation z) {
        int n4  = (n&15)<<27;
        int x9 = (x&511)<<18;
        int y9 = (y&511)<<9;
        int z9 = (z&511);
        int tag = n4|x9|y9|z9;
        if (tag < 2048) tag += 2048; // To avoid collision with registered tags
        return tag; 
    };

    /// Called by consumer to make a variable that will be set by producer
    template <typename T>
    typename Function<T>::ArgT Function<T>::input_arg(const OctTreeTPtr& consumer,
            const OctTreeTPtr& producer) {
        int tag = taghash(producer->n(),producer->x(),producer->y(),producer->z());
        if (consumer->islocal() && producer->islocal())
            return ArgT();
        else if (consumer->islocal() && producer->isremote())
            return ArgT(producer->rank(), tag, true, this->data->cdata->vk);
        else if (consumer->isremote() && producer->islocal())
            return ArgT(consumer->rank(), tag, false, this->data->cdata->vk);
        else
            throw "input_arg: should not happen?";
    };


    /// Compress function (scaling function to wavelet)

    /// Communication streams up the tree.
    /// Returns self for chaining.
    template <typename T>
    Function<T>& Function<T>::compress() {
        build_global_tree();
        if (!data->compressed) {
            if (isactive(tree())) {
                ArgT dummy;
                _compress(tree(),dummy);
            }
            taskq.global_fence();
            data->compressed = true;
        }
        return *this;
    };

    template <typename T>
    void Function<T>::_compress(OctTreeTPtr& tree, Function<T>::ArgT& parent) {
        ArgT args[2][2][2];
        FOREACH_ACTIVE_CHILD(OctTreeTPtr, tree,
                             args[i][j][k] = this->input_arg(tree,child);
                             if (child->islocal()) _compress(child,args[i][j][k]););

        if (tree->islocal()) {
            taskq.add_local(new TaskCompress<T>(this,tree,args,parent));
        }
    };

    template <typename T>
    void Function<T>::_compressop(OctTreeTPtr& tree, Function<T>::ArgT args[2][2][2], Function<T>::ArgT& parent) {
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
            unset_coeff(tree);
        }           
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
