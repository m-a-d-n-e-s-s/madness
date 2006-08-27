#include <iostream>
using std::cout;
using std::endl;

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

    /// A task tailored to the needs of compress2
    template <typename T>
    class TaskCompress : public TaskInterface {
    public:
        typedef Tensor<T> TensorT;
        typedef SAV<TensorT> ArgT;
        typedef OctTree<FunctionNode> OctTreeT;
        typedef Function<T> FunctionT;
    private:
        FunctionT* that;
        OctTreeTPtr tree;
        ArgT outarg;
        ArgT inarg[2][2][2];
    public:
        TaskCompress(FunctionT* that, OctTreeTPtr tree, ArgT in[2][2][2], ArgT& out)
                : that(that)
                , tree(tree)
                , outarg(out) 
        {
            FORIJK(inarg[i][j][k] = in[i][j][k];);
        };

        bool probe() const {
            FOREACH_CHILD(OctTreeTPtr, tree,
                          if (that->isactive(child) && !inarg[i][j][k].probe())
                          return false;);
            return true;
        };

        void run() {
            that->_compress2op(tree,inarg,outarg);
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
        if (!data->compressed) {
            if (isactive(tree())) {
                ArgT dummy;
                _compress2(tree(),dummy);
                taskq.local_fence();
            }
            data->compressed = true;
        }
        return *this;
    };

    template <typename T>
    void Function<T>::_compress2(OctTreeTPtr& tree, Function<T>::ArgT& parent) {
        ArgT args[2][2][2];
        FOREACH_ACTIVE_CHILD(OctTreeTPtr, tree,
                             args[i][j][k] = this->input_arg(tree,child);
                             if (child->islocal()) _compress2(child,args[i][j][k]););

        if (tree->islocal()) {
            //print(comm()->rank(),"adding task",tree->n(),tree->x(),tree->y(),tree->z());
            taskq.add_local(new TaskCompress<T>(this,tree,args,parent));
        }
    };

    template <typename T>
    void Function<T>::_compress2op(OctTreeTPtr& tree, Function<T>::ArgT args[2][2][2], Function<T>::ArgT& parent) {
        //print(comm()->rank(),"executing task",tree->n(),tree->x(),tree->y(),tree->z());
        Slice *s = data->cdata->s;
        if (!coeff(tree)) set_coeff(tree,TensorT(2*k,2*k,2*k));
        TensorT& t = *coeff(tree);
        FOREACH_ACTIVE_CHILD(OctTreeTPtr, tree,t(s[i],s[j],s[k]) += args[i][j][k].get(););
        filter_inplace(t);
        parent.set(madness::copy(t(data->cdata->s0)));
        if (tree->n() > 0) t(data->cdata->s0)=0.0;
    };
    
    template <typename T>
    void Function<T>::_dodiff_kernel(Function<T>& df, OctTreeTPtr& tree, int axis,
                    const TensorT& t0, const TensorT& tm, const TensorT& tp) {
        // Swap the dimensions so that we have effective axis = 0
        const TensorT q0 = t0.swapdim(axis,0);
        const TensorT qp = tp.swapdim(axis,0);
        const TensorT qm = tm.swapdim(axis,0);
        // Copy into contiguous array with dim0 arranged [-1,0,+1]
        TensorT q(4*k,2*k,2*k);
        q(Slice(0,k-1),_,_) = qm(Slice(k,2*k-1),_,_);
        q(Slice(k,3*k-1),_,_) = q0;
        q(Slice(3*k,4*k-1),_,_) = qp(Slice(0,k-1),_,_);
        // Loop over each block in the output, collect corresponding input blocks,
        // and apply the operator.  This is the slow but correct version. There is
        // a faster version that will be copied in once the parallel code is working.
        TensorT d(2*k,2*k,2*k); // The result
        Slice* s = data->cdata->s;
        const TensorT& rp_right = this->data->cdata->rp_right;
        const TensorT& rp_left = this->data->cdata->rp_left;
        const TensorT& rm_right = this->data->cdata->rm_right;
        const TensorT& rm_left = this->data->cdata->rm_left;
        const TensorT& r0 = this->data->cdata->r0;
        FORIJK(TensorT sm = q(s[i  ],s[j],s[k]);
               TensorT s0 = q(s[i+1],s[j],s[k]);
               TensorT sp = q(s[i+2],s[j],s[k]);
               d(s[i],s[j],s[k]) = ::inner(r0, s0) +
                    outer(rp_left,::inner(rp_right,sm)) +
                    outer(rm_left,::inner(rm_right,sp)););
            
        if (axis) d = ::copy(d.swapdim(axis,0));
        d.scale((double) two_to_power(tree->n()+1));
        df.set_coeff(tree,d);
    };
            

    template <typename T>
    void Function<T>::_recur_coeff_down(OctTreeT* tree, bool keep) {
        // 1,1,1 since this is the last entry that would be made
        // ensuring that all others are there also.
        if (tree->child(1,1,1) && isactive(tree->child(1,1,1))) {
            // This node has already been recurred down ... skip it.
            return;
        }
        
        if (tree->n() >= data->max_refine_level) 
            MADNESS_EXCEPTION("recur_coeff_down: exceeding max refine level",tree->n());
        
        // sonly on unfilter is a 3x optimization in 3d ... we need it!
        const Tensor<T>& c = *coeff(tree);
        long k2 = k*2;
        const std::vector<Slice>& s0 = data->cdata->s0;
        const Slice* s = data->cdata->s;
        FORIJK(OctTreeTPtr child = tree->child(i,j,k);
               if (!child) child = tree->insert_local_child(i,j,k);
               set_active(child);
               // If we are keeping the parent, the child will eventually be autocleaned
               set_acflag(child,keep);
               Tensor<T>*t = set_coeff(child, Tensor<T>(k2, k2, k2));
               (*t)(s0) = c(s[i],s[j],s[k]);
               unfilter_inplace(*t);  // sonly!
               //madness::print(child->n(),child->x(),child->y(),child->z(),"recurring");
               );
        if (!keep) unset_coeff(tree);
    }
    
    template <typename T>
    const Tensor<T>* Function<T>::_get_scaling_coeffs(OctTreeTPtr& t, int axis, int inc) {
        //madness::print("getting",t->n(),t->x(),t->y(),t->z());

        long xyz[3] = {t->x(), t->y(), t->z()};
        xyz[axis] += inc;
        long x=xyz[0], y=xyz[1], z=xyz[2];
        
        // Enforce boundary conditions (here use zero for external values)
        if (xyz[axis] < 0 || ((unsigned long) xyz[axis]) >= two_to_power(t->n())) return &data->cdata->zero_tensor;
        
        // Find node or its parent in the tree.
        OctTreeT* p = t->find(t->n(),x, y, z);
        // If not active, find the lowest active parent
        while (p && !isactive(p)) p = p->parent();

        if (! p) MADNESS_EXCEPTION("get_scaling_coeffs: failed to find parent?",0);
        if (! coeff(p)) return 0; //  The parent does not have data ... it is below.

        // Recur data down from parent to desired level.  Note, need to introduce
        // additional checking and mutex when we start running multithreaded. 
        while(p->n() < t->n()) {
            // MUTEX 
            // check again that someone else has not already made it
            _recur_coeff_down(p,true);
            // END MUTEX
            long nn = t->n() - p->n() - 1;
            long xx = (x>>nn)&1;
            long yy = (y>>nn)&1;
            long zz = (z>>nn)&1;
            p = p->child(xx,yy,zz).get();
        }

        return coeff(p);
    }
   
    
    template <typename T>
    void Function<T>::_dodiff(Function<T>& df, OctTreeTPtr& tree, int axis) {
        //madness::print("dodiff",tree->n(),tree->x(),tree->y(),tree->z());

        // This routine is applied at each leaf node ... i.e., where we have data in tree.
        df.set_active(tree);
        const TensorT* t0 = coeff(tree);
        const TensorT* tm = _get_scaling_coeffs(tree,axis,-1);
        const TensorT* tp = _get_scaling_coeffs(tree,axis, 1);
        if (tm && tp) {
            _dodiff_kernel(df, tree, axis, *t0, *tm, *tp);
        }
        else {
            _recur_coeff_down(tree.get(),true);
            FOREACH_CHILD(OctTreeTPtr, tree, _dodiff(df, child, axis);); 
        }          
    };
    

    template class Function<double>;
    template class Function< std::complex<double> >;
    template class FunctionData<double>;
    template class FunctionData< std::complex<double> >;
    template class FunctionCommonData<double>;
    template class TaskCompress<double>;
    template class TaskCompress< std::complex<double> >;    

}
