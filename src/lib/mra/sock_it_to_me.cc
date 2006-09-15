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


/// \file mra/sock_it_to_me.cc
/// \brief Implements parallel getting and recursion of scaling function coefficients


namespace madness {
    
    template <typename T>
    void Function<T>::set_active_handler(Communicator& comm, ProcessID src, const AMArg& arg) {
        int ind = arg.arg0;
        Level n = arg.arg1;
        Translation x=arg.arg2, y=arg.arg3, z=arg.arg4;
        Function<T> f = Function<T>(FunctionDataPointers<T>::get(ind));
        OctTreeTPtr t = f.tree();
        MADNESS_ASSERT(t);
        if (t->n()==n && t->x()==x && t->y()==y && t->z()==z) {
            f.set_active(t);
        }
        else {
            MADNESS_EXCEPTION("set_active_handler: confused",0);
        }
    };
    
    template <typename T>
    void Function<T>::recur_down_handler(Communicator& comm, ProcessID src, VectorInputArchive& ar) {
        Tensor<T> c;
        int ind;
        Level n;
        Translation x, y, z;
        bool keep;
        ar & ind & n & x & y & z & keep & c;
        Function<T> f = Function<T>(FunctionDataPointers<T>::get(ind));
        OctTreeTPtr t = f.tree();
        MADNESS_ASSERT(t);
        MADNESS_ASSERT(t->n()==n-1 && t->x()==(x>>1) && t->y()==(x>>1) && t->z()==(z>>1));
        int xx=x&1, yy=y&1, zz=z&1;
        OctTreeTPtr child = t->child(xx,yy,zz);
        if (!child) child = t->insert_local_child(xx,yy,zz);
        if (!f.isactive(t)) {
            madness::print("WARNING: recur_down_handler: parent was not active but should have been");
            f.set_active(t);
        }
        f.set_active(child);
        f.set_acflag(child,keep);
        long k2 = 2*f.k;
        Tensor<T>& s = *f.set_coeff(child, Tensor<T>(k2, k2, k2));
        s(f.data->cdata->s0) = c;
        f.unfilter_inplace(s); //!!!! sonly
    };
        
    
    template <typename T>
    void Function<T>::_recur_coeff_down2(OctTreeT* tree, bool keep) {
        MADNESS_ASSERT(tree);
        if (tree->child(1,1,1) && isactive(tree->child(1,1,1))) {
            // This node has already been recurred down ... skip it.
            // 1,1,1 since this is the last entry that would be made
            // ensuring that all others are there also.
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
               set_acflag(child,keep);  // If keeping parent, child will eventually be autocleaned
               if (islocal(child)) {
                    Tensor<T>*t = set_coeff(child, Tensor<T>(k2, k2, k2));
                    (*t)(s0) = c(s[i],s[j],s[k]);
                    unfilter_inplace(*t);  // sonly!
                    // Sigh.  If the child has any remote children of its
                    // own then they need to be notified of the change in
                    // their parent's status.  This only needs to be done
                    // if after this recur_down we are to perform an
                    // upward traversal of the tree (e.g., compress) without
                    // autocleaning which only makes sense if keep=false.
                    if (!keep) {
                        ProcessID ranks[8];
                        int np = unique_active_child_procs(child,ranks);
                        AMArg arg(ind, child->n(), child->x(), child->y(), child->z());
                        for (int i=0; i<np; i++) {
                            comm()->am_send(ranks[i], set_active_handler, arg);
                        }
                    }
               }
               else {
                    std::vector<unsigned char> v;
                    VectorOutputArchive ar(v);
                    ar & ind & child->n() & child->x() & child->y() & child->z() & keep & madness::copy(c(s[i],s[j],s[k]));
                    taskq.add_generic(child->rank(),recur_down_handler,v);
               }
               );
        if (!keep) unset_coeff(tree);
    }
    
    template <typename T>
    void Function<T>::_sock_it_to_me_handler(Communicator& comm, ProcessID src, const AMArg& arg) {
        int ind = arg.arg0;
        Level n = arg.arg1;
        Translation l[3] = {arg.arg2, arg.arg3, arg.arg4};
        ProcessID requestor = arg.arg5;
        int tag = arg.arg6;
        
        Function<T> f = Function<T>(FunctionDataPointers<T>::get(ind));
        SAV< Tensor<T> > result(requestor, tag, false);
        f._sock_it_to_me(f.tree(), n, l, result);
    }
    
    template <typename T>
    void Function<T>::_sock_it_to_me_forward_request(Level n, const Translation l[3], SAV< Tensor<T> >& arg, ProcessID dest) {
        MADNESS_EXCEPTION("Uh?",0);
        if (arg.islocal()) {
            arg = SAV< Tensor<T> >(MPI::ANY_SOURCE, comm()->unique_tag(), true, data->cdata->v2k);
        }
        ProcessID requestor;
        int tag;
        arg.msginfo(requestor, tag);
        AMArg amarg(ind, n, l[0], l[1], l[2], requestor, tag);
        comm()->am_send(dest, _sock_it_to_me_handler, amarg);     
    }    
                 
                 
    // simpler sock it to me logic might be to just forward the recur down request (to make n,l)
    // but not to forward the SAV beyond the actual owner.  In this way we can just have a task
    // that waits for the coeffs to appear (non-null pointer) and then assigns the SAV.
    // The latency will be a little greater but if it is inserted as an HP task with the
    // assignment done by the actual probe we will be in business.  The local SAV can 
    // also just wrap the pointer (assuming the tree exists).         
                 
                                                        
    template <typename T>
    void Function<T>::_sock_it_to_me(OctTreeTPtr& tree, Level n, const Translation l[3], SAV< Tensor<T> >& arg) {
        Translation x=l[0], y=l[1], z=l[2];
        // Most requests will resolve locally ... try to handle those efficiently
       
        // Find node or its parent in the tree.
        OctTreeT* p = tree->find(tree->n(), x, y, z);
        
        if (p) { // Found the node or a parent thereof
            if (p->n() == n  &&  isactive(p)) { // The node itself and is active
                if (coeff(p)) { // Has coeff ... success!
                    arg.set(*coeff(p));
                    return;
                }
                else if (islocal(p)) { // No coeff but active, so coeff are below ... failure
                    arg.set(Tensor<T>());
                    return;
                }
                else { // Active but not local ... forward request to owner
                    madness::print("forward");
                    _sock_it_to_me_forward_request(n,l,arg,p->rank());
                    return;
                }
            }
            
            // Look up the local tree for the closest parent with coeff
            while (!coeff(p) && p->parent()) p = p->parent();
            
            if (coeff(p)) { // Local parent has some coeff ... recur down
                Translation x=l[0], y=l[1], z=l[2];
                while(p->n() < n) {
                    _recur_coeff_down2(p,true);
                    long nn = n - p->n() - 1;
                    long xx = (x>>nn)&1;
                    long yy = (y>>nn)&1;
                    long zz = (z>>nn)&1;
                    p = p->child(xx,yy,zz).get();
                    if (isremote(p)) {
                        _sock_it_to_me_forward_request(n,l,arg,p->rank());
                        return;
                    }
                }
                arg.set(*coeff(p));
                return;
            }
            else if (isremote(p)) {
                _sock_it_to_me_forward_request(n,l,arg,p->rank());
                return;
            }
            else {
                MADNESS_EXCEPTION("sock_it_to_me: logic error",1);
            } 
        }
        
        MADNESS_EXCEPTION("not yet",0);
        // Don't have enough info, try the global tree
        ProcessID owner = find_owner(n,l);
        if (owner == comm()->rank()) {
            MADNESS_EXCEPTION("sock_it_to_me: logic error",2);
        }
        else {
            _sock_it_to_me_forward_request(n,l,arg,owner);
        }
    }
    
    template <typename T>
    SAV< Tensor<T> > Function<T>::_get_scaling_coeffs2(OctTreeTPtr& t, int axis, int inc) {
        MADNESS_ASSERT(t);
        Translation xyz[3] = {t->x(), t->y(), t->z()};
        xyz[axis] += inc;
        
        // !!!!!!!!!!! relying here on 0u-1 goes to big +ve no. ?????
        // Enforce boundary conditions (here use zero for external values)
        if (xyz[axis]<0 || xyz[axis]>=two_to_power(t->n())) 
            return SAV< Tensor<T> >(data->cdata->zero_tensor);
    
        SAV< Tensor<T> > result;
        _sock_it_to_me(t, t->n(), xyz, result);
        return result;
    };
    
    template void Function<double>::_recur_coeff_down2(OctTreeT* tree, bool keep);
    template void Function<double_complex>::_recur_coeff_down2(OctTreeT* tree, bool keep);
    template SAV< Tensor<double> > Function<double>::_get_scaling_coeffs2(OctTreeTPtr& t, int axis, int inc);
    template SAV< Tensor<double_complex> > Function<double_complex>::_get_scaling_coeffs2(OctTreeTPtr& t, int axis, int inc);
    
}

   
    
