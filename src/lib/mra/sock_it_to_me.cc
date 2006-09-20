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
#include <cmath>


/// \file mra/sock_it_to_me.cc
/// \brief Implements parallel getting and recursion of scaling function coefficients


namespace madness {
    
    template <typename T>
    void Function<T>::set_active_handler(Communicator& comm, ProcessID src, const AMArg& arg) {
        int ind = arg.arg0;
        Level n = arg.arg1;
        Translation x=arg.arg2, y=arg.arg3, z=arg.arg4;
        madness::print("set_active_handler",ind,n,x,y,z);
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
        madness::print("recur_down_handler",ind,n,x,y,z,keep);
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
    void Function<T>::recur_down_to_make_handler(Communicator& comm, ProcessID src, const AMArg& arg) {
        int ind = arg.arg0;
        Level n = arg.arg1;
        Translation l[3] = {arg.arg2, arg.arg3, arg.arg4};
        madness::print("recur_down_to_make_handler",ind, n, l[0], l[1], l[2],src);
        
        Function<T> f = Function<T>(FunctionDataPointers<T>::get(ind));
        OctTreeT* p = f.tree()->find(n,l[0],l[1],l[2]);
        f.recur_down_to_make(p, n, l);
    }
    
    template <typename T>
    void Function<T>::recur_down_to_make_forward_request(Level n, const Translation l[3], ProcessID dest) { 
        AMArg amarg(ind, n, l[0], l[1], l[2]);
        madness::print("rdtm forwarding",n, l[0], l[1], l[2]);
        comm()->am_send(dest, recur_down_to_make_handler, amarg);     
    }    

    
    /// Recur down to make specified location
    template <typename T>
    void Function<T>::recur_down_to_make(OctTreeT* p, Level n, const Translation l[3]) {
        MADNESS_ASSERT(p);
        Translation x=l[0], y=l[1], z=l[2];
        while (p->n() < n) {
            MADNESS_ASSERT(p->islocal());
            _recur_coeff_down2(p,true);
            long nn = n - p->n() - 1;
            long xx = (x>>nn)&1;
            long yy = (y>>nn)&1;
            long zz = (z>>nn)&1;    
            p = p->child(xx,yy,zz);
            if (p->isremote()) {
                // This needs to be combined with the forwarding
                // implicit in recur_down
                recur_down_to_make_forward_request(n,l,p->rank());
                return;
            }
        }
        madness::print("rdtm just made",n,l[0],l[1],l[2],p->x(),p->y(),p->z(),coeff(p)->normf());
    }     
    
    
    template <typename T>
    void Function<T>::_sock_it_to_me_handler(Communicator& comm, ProcessID src, const AMArg& arg) {
        int ind = arg.arg0;
        Level n = arg.arg1;
        Translation l[3] = {arg.arg2, arg.arg3, arg.arg4};
        int tag = arg.arg5;
        madness::print("in sock handler",n, l[0], l[1], l[2],src,tag);
        
        Function<T> f = Function<T>(FunctionDataPointers<T>::get(ind));
        SAV< Tensor<T> > result(src, tag, false, f.data->cdata->v2k);
        f._sock_it_to_me(f.tree(), n, l, result);
        madness::print("after f.sock",n, l[0], l[1], l[2], result.probe());
        //MADNESS_ASSERT(result.probe());
    }
    
    template <typename T>
    void Function<T>::_sock_it_to_me_forward_request(Level n, const Translation l[3], SAV< Tensor<T> >& arg, ProcessID dest) {
        // Reassign result to be set remotely via a message from process dest 
        int tag = comm()->unique_tag();
        madness::print("sock forwarding",n, l[0], l[1], l[2],tag);
        MADNESS_ASSERT(arg.islocal());
        arg = SAV< Tensor<T> >(dest, tag, true, data->cdata->v2k);
        AMArg amarg(ind, n, l[0], l[1], l[2], tag);
        comm()->am_send(dest, _sock_it_to_me_handler, amarg);     
    }    
                 
                 
                 
    /// Fill in missing local nodes down to specified location
    template <typename T>
    void Function<T>::fill_in_local_tree(OctTreeT* p, Level n, const Translation l[3]) {
        MADNESS_ASSERT(p);
        Translation x=l[0], y=l[1], z=l[2];
        while (p->n() < n) {
            MADNESS_ASSERT(p->islocal() || p->parent()==0);
            // All children must exist
            FORIJK(if (!p->child(i,j,k)) p->insert_local_child(i,j,k););
            long nn = n - p->n() - 1;
            long xx = (x>>nn)&1;
            long yy = (y>>nn)&1;
            long zz = (z>>nn)&1;    
            p = p->child(xx,yy,zz);
        }
    }     
    
    template <typename T>
    class TaskAwaitCoeff : public TaskInterface {
    private:
        Function<T>* f;
        OctTreeT* p;
        mutable SAV< Tensor<T> > result;
    public:
        TaskAwaitCoeff(Function<T>* f, OctTreeT* p, SAV< Tensor<T> >& result) 
        : f(f), p(p), result(result) 
        {}
        
        bool probe() const {
            if (!f->coeff(p)) return false;
            if (!result.probe()) result.set(*f->coeff(p));
            return true;
        }
        void run() {};   
    };   
    
                                                       
               
    template <typename T>
    class TaskRecurDownToMakeLocal : public TaskInterface {
    private:
	Function<T>* f;
	OctTreeT* p;
	Level n;
	SAV< Tensor<T> > result;
	Translation l[3];
	
    public:
	TaskRecurDownToMakeLocal(Function<T>* f, OctTreeT* p, Level n, const Translation l[3], SAV< Tensor<T> >& result) 
	    : f(f),p(p),n(n),result(result) {
	    for (int i=0; i<3; i++) this->l[i] = l[i];
	};
	void run() {
	    madness::print("TaskRecurDownToMakeLocal",f->ind, p->n(),p->x()<p->y(),p->z(),n,l[0],l[1],l[2]);
	    //madness::print("TaskRecurDownToMakeLocal",(void*)f,(void*)p,n,l[0],l[1],l[2]);
	    f->recur_down_to_make(p,n,l);
	    result.set(*f->coeff(p->find(n,l[0],l[1],l[2])));
	};
	bool probe() const {return true;};
    };
    
    
    template <typename T>
    void Function<T>::_sock_it_to_me(OctTreeTPtr& tree, Level n, const Translation l[3], SAV< Tensor<T> >& result) {
        Translation x=l[0], y=l[1], z=l[2];
        
        OctTreeT* p = tree->find(n, x, y, z);  // Look for node/parent in local tree.
        if (p && p->n()==n  &&  isactive(p)) { // The node itself and is active
            if (coeff(p)) { // Has coeff ... success!
                result.set(*coeff(p));
                return;
            }
            else if (islocal(p)) { // No coeff but active, so coeff are below ... failure
                //madness::print("SOCK FAILED",n, x, y, z);
                result.set(Tensor<T>());
                return;
            }
            else {
                _sock_it_to_me_forward_request(n,l,result,p->rank());
                return;
            }
        }


        ProcessID owner = find_owner(n,l);  // Don't have enough info, look in global tree
        if (owner == comm()->rank()) {
            if (p) {
                // Look up the local tree for the closest parent with coeff
                while (!isactive(p) && p->parent()) p = p->parent();
                fill_in_local_tree(p, n, l);
                if (coeff(p)) {
		    if (result.islocal()) {
			taskq.add_local(new TaskRecurDownToMakeLocal<T>(this,p,n,l,result));
		    }
		    else { // Immediate gratification for remote requests
			recur_down_to_make(p,n,l); 
			result.set(*coeff(p->find(n,x,y,z)));
		    }
                    return;
                }
                else if (isremote(p)) {
                    taskq.add_local(new TaskAwaitCoeff<T>(this,p->find(n,x,y,z),result));       
                    recur_down_to_make_forward_request(n,l,p->rank());
                    return;
                }
                else {
                    MADNESS_EXCEPTION("sock_it_to_me: logic error",1);
                } 
            }
            else {
                madness::print("Gtree thinks I own",n,l[0],l[1],l[2]);
                pnorms();            
                MADNESS_EXCEPTION("sock_it_to_me: logic error",2);
            }
        }
        else {
            //madness::print("Forwarding request after gtree",n,l[0],l[1],l[2]);
            _sock_it_to_me_forward_request(n,l,result,owner);
        }
    }
    
    template <typename T>
    SAV< Tensor<T> > Function<T>::_get_scaling_coeffs2(OctTreeTPtr& t, int axis, int inc) {

//        const Tensor<T>* ss = _get_scaling_coeffs(t, axis, inc);
//        if (ss) {
//            return SAV< Tensor<T> >(*ss);
//        }
//        else if (!ss) {
//            return SAV< Tensor<T> >(Tensor<T>());
//        }

        MADNESS_ASSERT(t);
        Translation xyz[3] = {t->x(), t->y(), t->z()};
        xyz[axis] += inc;
        
        // Enforce boundary conditions (here use zero for external values)
        // !!!!!!!!!!! relying here on 0u-1 goes to big +ve no. ?????
        //if (xyz[axis]<0 || xyz[axis]>=two_to_power(t->n())) {
        if (xyz[axis]>=two_to_power(t->n())) {
            //madness::print("BOUNDARY",t->n(), xyz[0], xyz[1], xyz[2]);
            return SAV< Tensor<T> >(data->cdata->zero_tensor);
        }
    
        SAV< Tensor<T> > result;
        _sock_it_to_me(t, t->n(), xyz, result);
        
//        const Tensor<T>* s = _get_scaling_coeffs(t, axis, inc);
//        
//        if (s) {
//            MADNESS_ASSERT(result.get().size);
//            if (fabs(s->normf()-result.get().normf()) > 1e-12) {
//                madness::print("OLD vs NEW",t->n(), xyz[0], xyz[1], xyz[2],s->normf(),result.get().normf());
//            }
//        }
//        else {
//            MADNESS_ASSERT(result.get().size == 0);
//        }
        return result;
    };
    
    template void Function<double>::_recur_coeff_down2(OctTreeT* tree, bool keep);
    template void Function<double_complex>::_recur_coeff_down2(OctTreeT* tree, bool keep);
    template SAV< Tensor<double> > Function<double>::_get_scaling_coeffs2(OctTreeTPtr& t, int axis, int inc);
    template SAV< Tensor<double_complex> > Function<double_complex>::_get_scaling_coeffs2(OctTreeTPtr& t, int axis, int inc);
    
}

   
    
