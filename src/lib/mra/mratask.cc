#include <complex>
#include <octtree/octtree.h>
#include <mra/mra.h>
#include <mra/twoscale.h>
#include <mra/legendre.h>
#include <tasks/sav.h>
#include <tasks/tasks.h>
#include <misc/print.h>


/// \file mra/mratask.cc
/// \brief Implements more Function stuff that uses tasks


namespace madness {
    
    void mratask_register();

    /// Used for forwardable Function tasks applied to scaling function coeffs
    
    /// Derived types need to implement probe(), apply(), and predicate()
    /// Note the derived type is used to parameterize the base class
    /// so we can have access to its constructor and AM handler.
    template <typename T, class Derived>
    class TaskLeaf : public TaskInterface {
        friend void mratask_register();
    protected:
        OctTreeTPtr tree;
        Function<T> in1;
        Function<T> in2;
        Function<T> out1;
        bool keep;
        int axis;
        Level n;
        Translation l[3];
        
        struct Payload {
            Translation l1[3];
            Translation l2[3];
            short in1;
            short in2;
            short out1;
            char axis;
            char n1;
            char n2;
        };
        
        static void forward_task_handler(Communicator& comm, ProcessID src, const AMArg& amarg) {
            Payload arg;
            amarg.copyout(&arg,sizeof(arg));
            madness::print("FORWARD TASK HANDLER",arg.in1,arg.in2,arg.out1,arg.axis,arg.n1);
            // This function is static! The local variables shadow the class data
            Function<T> in1 = Function<T>(arg.in1);
            Function<T> in2 = Function<T>(arg.in2);
            Function<T> out1 = Function<T>(arg.out1);
            int axis = arg.axis;
            Level n = arg.n1;
            Translation l[3]= {arg.l1[0], arg.l1[1], arg.l1[2]};
            Translation x=arg.l2[0], y=arg.l2[1], z=arg.l2[2];
            OctTreeTPtr tree = OctTreeT::find_down(in1.tree(),n,x,y,z);
            if (!tree) MADNESS_EXCEPTION("TaskLeaf: SOL looking for tree?",0);
            taskq.add_local_hp(new Derived(tree, in1, in2, out1, axis, n, l));            
        };
        
        void forward(const OctTreeTPtr& t) const {
            madness::print("FORWARDING!!!!!!!!!!");
            Payload arg;
            arg.in1 = in1.ind;
            arg.in2 = in2.ind;
            arg.out1 = out1.ind;
            arg.n1 = n;
            arg.l1[0]=l[0]; arg.l1[1]=l[1]; arg.l1[2]=l[2];
            arg.n2 = tree->n(); 
            arg.l2[0]= tree->x(); arg.l2[1]= tree->y(); arg.l2[2]= tree->z();        
            taskq.add_am(t->rank(), Derived::forward_task_handler, AMArg(&arg,sizeof(arg)));
        };
        
        
    public:
        TaskLeaf(OctTreeTPtr& tree,
                 Function<T>& in1, Function<T>& in2, Function<T>& out1,
                 int axis, Level n, const Translation *l) 
        : tree(tree), in1(in1), in2(in2), out1(out1), keep(true), axis(axis), n(n)         
        {
            if (l) for (int i=0; i<3; i++) this->l[i] = l[i];
        };
        
        virtual ~TaskLeaf() {};
        
        virtual void apply() = 0;
                 
        virtual bool predicate() const = 0;
        

        /// Applies operator where predicate is true otherwise recurs down
        
        /// Recursion generates more tasks which may be forwarded
        void run() {
            MADNESS_ASSERT(tree);
            out1.set_active(tree);
            if (predicate()) {
                apply();
            }
            else {
                if (in1.coeff(tree)) in1._recur_coeff_down2(tree,keep);
                if (in1.ind!=in2.ind  &&  in2.coeff(tree)) in2._recur_coeff_down2(tree,keep);
                FOREACH_CHILD(OctTreeTPtr, tree, 
                              if (tree->islocal()) { 
                                  taskq.add_local(new Derived(child, in1, in2, out1, axis, n, l));
                              }
                              else {
                                  out1.set_active(child);
                                  forward(child);
                              });
            };           
        };
        
        /// If in1 or in2 has coeff spawns a task, otherwise recurs down active, local children 
        static void generate_tasks(OctTreeTPtr& tree,
                                   Function<T>& in1, Function<T>& in2, Function<T>& out,
                                    int axis=0, Level n=0, const Translation *l=0) {
            comm_default->am_poll();
            MADNESS_ASSERT(tree);
            if (in1.coeff(tree) || in2.coeff(tree)) {
                taskq.add_local(new Derived(tree, in1, in2, out, axis, n, l));
            }
            else if (in1.isactive(tree) || in2.isactive(tree)) {
                out.set_active(tree);
                FOREACH_CHILD(OctTreeTPtr, tree, 
                              if (in1.isactive(child) || in2.isactive(child)) {
                                  if (child->islocal()) 
                                      generate_tasks(child,in1,in2,out,axis,n,l);
                                  else
                                      out.set_active(child);
                              });
            }    
        };      
    };
    
    template <typename T>
    class TaskDiff : public TaskLeaf< T, TaskDiff<T> > {
        friend void mratask_register();
    private:
        typedef TaskLeaf< T, TaskDiff<T> > baseT;
        SAV< Tensor<T> > left;
        SAV< Tensor<T> > right;

    public:        
        TaskDiff(OctTreeTPtr& tree,
                 Function<T>& in1, Function<T>& in2, Function<T>& out1,
                 int axis, Level n, const Translation *l)
        : baseT(tree, in1, in2, out1, axis, n, l) {
            MADNESS_ASSERT(tree);
            left  = in1._get_scaling_coeffs2(tree, axis, -1);
            right = in1._get_scaling_coeffs2(tree, axis,  1);        
        };
        
        bool probe() const {
            return left.probe() && right.probe();
        };
        
        bool predicate() const {
            return left.get().size>0 && right.get().size>0;
        };
        
        void apply() {
            MADNESS_ASSERT(baseT::tree);
            MADNESS_ASSERT(baseT::in1.coeff(baseT::tree));
            const Tensor<T>& t0 = *baseT::in1.coeff(baseT::tree);
            const Tensor<T>& tm = left.get(); 
            const Tensor<T>& tp = right.get();
            baseT::in1._dodiff_kernel(baseT::out1, baseT::tree, baseT::axis, t0, tm, tp);
        };
    };
    
    /// Differentiation high level interface
    template <typename T>
    Function<T> Function<T>::diff(int axis) {
        reconstruct();
        build_global_tree();
        Function<T> df = FunctionFactory<T>().k(k).compress(false).empty();
        TaskDiff<T>::generate_tasks(tree(), *this, *this, df, axis);
        taskq.global_fence();
        _auto_clean(tree()); 
        clear_global_tree();
        return df;
    };
    
    
    
    
    /// Task for fine_scale_projection and refinement
    template <typename T>
    class TaskProjectRefine : public TaskInterface {
        friend void mratask_register();
    private:
        OctTreeTPtr tree;
        Function<T> f;
        
        static void forward_task_handler(Communicator& comm, ProcessID src, const AMArg& arg) {
            int ind = arg.arg0;
            Level n = arg.arg1;
            Translation x=arg.arg2, y=arg.arg3, z=arg.arg4;
            // This function is static! The local variables shadow the class data
            Function<T> f = Function<T>(ind);
            OctTreeTPtr tree = f.tree();
            f.set_active(tree);
            tree = tree->find_down(tree,n,x,y,z);
            if (!tree) MADNESS_EXCEPTION("TaskProjectRefine: SOL looking for tree?",0);
            taskq.add_local(new TaskProjectRefine<T>(tree, f));  
        };
        
        void forward(const OctTreeTPtr& t) const {
            taskq.add_am(t->rank(), TaskProjectRefine<T>::forward_task_handler, 
                        AMArg(f.ind,t->n(),t->x(),t->y(),t->z()));
        };
        
        
    public:
        TaskProjectRefine(OctTreeTPtr& tree, Function<T>& f) 
        : tree(tree), f(f) {};
        
        bool probe() const {return true;};
        
        void run() {
            f.set_active(tree);
            if (f._doproject_refine(tree)) {
                FORIJK(OctTreeTPtr child = tree->child(i,j,k);
                       if (!child) child = tree->insert_local_child(i,j,k);
                       f.set_active(child);
                       if (child->islocal()) 
                           taskq.add_local(new TaskProjectRefine<T>(child, f));
                       else
                           forward(child););
            }
        };
        
        /// Generates tasks at level n of the tree 
        static void generate_tasks(OctTreeTPtr& tree, Function<T>& f, Level n) {
            madness::comm_default->am_poll();
            MADNESS_ASSERT(tree);
            if (tree->n() < n) {
                f.set_active(tree);
                FORIJK(OctTreeTPtr child = tree->child(i,j,k);
                       if (!child && tree->islocal()) 
                           child = tree->insert_local_child(i,j,k);
                       if (child) {
                           f.set_active(child);
                           if (f.islocal(child)) generate_tasks(child,f,n);
                       });
            }
            else if (tree->n() == n) {
               f.set_active(tree);
               if (tree->islocal()) taskq.add_local(new TaskProjectRefine<T>(tree, f));     
            }
            else {
                f.set_inactive(tree);
            }                
        };                     
    };
    
    /// Projection/refining high level interface
    template <typename T>
    void Function<T>::project_refine() {
        TaskProjectRefine<T>::generate_tasks(tree(), *this, data->initial_level);
        taskq.global_fence();
    };
    
    template <typename T>
    void Function<T>::_doreconstruct(OctTreeTPtr tree, const Tensor<T>& ss) {
        if (coeff(tree)) {
            Tensor<T>& d = *coeff(tree);
            Slice* s = data->cdata->s;
            if (tree->n()>0) d(s[0],s[0],s[0]) = ss;
            unfilter_inplace(d);
            FOREACH_CHILD(OctTreeTPtr, tree,
                          MADNESS_ASSERT(child);
                          MADNESS_ASSERT(isactive(child));
                          if (child->islocal()) {
                              _doreconstruct(child,madness::copy(d(s[i],s[j],s[k])));
                          }
                          else {
                              std::vector<unsigned char> v;
                              VectorOutputArchive ar(v);
                              ar & ind & child->n() & child->x() & child->y() & child->z() & madness::copy(d(s[i],s[j],s[k]));
                              taskq.add_generic(child->rank(), reconstruction_handler, v);
                          });
            unset_coeff(tree);
        }
        else if (tree->islocal()) {
            set_coeff(tree,ss);
        }
    };

    /// Generic task interface for reconstruction
    template <typename T>
    void Function<T>::reconstruction_handler(Communicator& comm, ProcessID src, VectorInputArchive& ar) {
        int ind;
        Level n;
        Translation x, y, z;
        Tensor<T> s;
        ar & ind & n & x & y & z & s;
        madness::print("RECONSTRUCTION HANDLER",ind,n,x,y,z,s.size);
        Function<T> f(ind);
        OctTreeTPtr tree = f.tree();
        f._doreconstruct(tree->find_down(tree,n,x,y,z),s);
    }
    
    
    template <typename T>
    class TaskAutorefine : public TaskLeaf< T, TaskAutorefine<T> > {
        friend void mratask_register();
    private:
        typedef TaskLeaf< T, TaskAutorefine<T> > baseT;

    public:        
        TaskAutorefine(OctTreeTPtr& tree,
                       Function<T>& in1, Function<T>& in2, Function<T>& out1,
                       int axis, Level n, const Translation *l)
        : baseT(tree, in1, in2, out1, axis, n, l) {
            MADNESS_ASSERT(tree);
            baseT::keep = false;
        };
        
        bool probe() const {return true;};
        
        bool predicate() const {return !baseT::in1.autorefine_test(baseT::tree);};
        
        void apply() {};
    };
    
    template <typename T>
    Function<T>& Function<T>::autorefine(double tol) {
        reconstruct();
        if (tol <= 0.0) tol = data->thresh;
        data->autorefine_thr = tol;
        TaskAutorefine<T>::generate_tasks(this->tree(), *this, *this, *this);
        taskq.global_fence();
        return *this;
    }
    

    template <typename T>
    class TaskSquare : public TaskLeaf< T, TaskSquare<T> > {
        friend void mratask_register();
    private:
        typedef TaskLeaf< T, TaskSquare<T> > baseT;

    public:        
        TaskSquare(OctTreeTPtr& tree,
                   Function<T>& in1, Function<T>& in2, Function<T>& out1,
                   int axis, Level n, const Translation *l)
        : baseT(tree, in1, in2, out1, axis, n, l) {
            MADNESS_ASSERT(tree);
            baseT::keep = false;
        };
        
        bool probe() const {return true;};
        
        bool predicate() const {return !baseT::in1.autorefine_square_test(baseT::tree);};
        
        void apply() {baseT::in1.do_square(baseT::tree);};
    };

    
    template <typename T>
    Function<T>& Function<T>::square(double tol) {
        reconstruct();
        if (tol <= 0.0) tol = data->thresh;
        data->autorefine_thr = tol;
        
        TaskSquare<T>::generate_tasks(this->tree(), *this, *this, *this);
        taskq.global_fence();
        return *this;
    }

    template <typename T>
    class TaskMult : public TaskLeaf< T, TaskMult<T> > {
        friend void mratask_register();
    private:
        typedef TaskLeaf< T, TaskMult<T> > baseT;
        SAV< Tensor<T> > a;
        SAV< Tensor<T> > b;
    public:        
        TaskMult(OctTreeTPtr& tree,
                 Function<T>& in1, Function<T>& in2, Function<T>& out1,
                 int axis, Level n, const Translation *l)
        : baseT(tree, in1, in2, out1, axis, n, l) {
            MADNESS_ASSERT(tree);

            Level nn = tree->n();
            Translation ll[3] = {tree->x(),tree->y(),tree->z()};

            if (in1.coeff(tree)) a.set(*in1.coeff(tree));
            else in1._sock_it_to_me(tree, nn, ll, a);

            if (in2.coeff(tree)) b.set(*in2.coeff(tree));
            else in2._sock_it_to_me(tree, nn, ll, b);
        };
        
        bool probe() const {return a.probe() && b.probe();};
        
        bool predicate() const {
            return a.get().size>0 && b.get().size>0 &&
                   !baseT::in1.autorefine_mult_test(baseT::tree,a.get(),b.get());
        };
        
        void apply() {baseT::out1.do_mult(baseT::tree, a.get(), b.get());};
    };
    
    template <typename T>
    Function<T> Function<T>::mult(Function<T>& other, double tol) {
        reconstruct();
        other.reconstruct();
        Function<T> result = FunctionFactory<T>().k(k).compress(false).empty();
        if (tol <= 0.) tol = result.data->thresh;
        result.data->autorefine_thr = tol;
        TaskMult<T>::generate_tasks(this->tree(), *this, other, result);
        taskq.global_fence();
        _auto_clean(tree());
        other._auto_clean(tree());
        return result;
    }


    
    void mratask_register() {
#define REGAM(f) print("          AM handler",#f,madness::comm_default->am_register(f))
#define REGGE(f) print("Generic task handler",#f,taskq.register_generic_op(f))
        
        REGAM(TaskProjectRefine<double>::forward_task_handler);
        REGAM(TaskProjectRefine<double_complex>::forward_task_handler); 
        
        REGGE(Function<double>::reconstruction_handler);
        REGGE(Function<double_complex>::reconstruction_handler);
    }  
  
    template Function<double> Function<double>::diff(int axis);
    template Function<double_complex> Function<double_complex>::diff(int axis);
    template void Function<double>::project_refine();
    template void Function<double_complex>::project_refine();
    template Function<double>& Function<double>::autorefine(double);
    template Function<double_complex>& Function<double_complex>::autorefine(double);
    template Function<double>& Function<double>::square(double);
    template Function<double_complex>& Function<double_complex>::square(double);
    template Function<double> Function<double>::mult(Function<double>& other, double);
    template Function<double_complex> Function<double_complex>::mult(Function<double_complex>& other, double);

}
