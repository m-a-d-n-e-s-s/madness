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


/// \file mra/mratask.cc
/// \brief Implements more Function stuff that uses tasks


namespace madness {

    /// Used for forwardable Function tasks applied to leaves
    
    /// Derived types need to implement probe(), apply(), and predicate()
    /// Note the derived type is used to parameterize the base class
    /// so we can have access to its constructor and AM handler.
    template <typename T, class Derived>
    class TaskLeaf : public TaskInterface {
    protected:
        OctTreeTPtr tree;
        Function<T> in1;
        Function<T> in2;
        Function<T> out1;
        int axis;
        Level n;
        Translation l[3];
        
        struct Payload {
            int in1;
            int in2;
            int out1;
            int axis;
            Level n1;
            Level n2;
            Translation l1[3];
            Translation l2[3];
        };
        
        static void forward_task_handler(Communicator& comm, ProcessID src, const AMArg& amarg) {
            const Payload& arg = (const Payload&) amarg;
            // This function is static! The local variables shadow the class data
            Function<T> in1 = Function<T>(FunctionDataPointers<T>::get(arg.in1));
            Function<T> in2 = Function<T>(FunctionDataPointers<T>::get(arg.in2));
            Function<T> out1 = Function<T>(FunctionDataPointers<T>::get(arg.out1));
            int axis = arg.axis;
            Level n = arg.n1;
            Translation l[3]= {arg.l1[0], arg.l1[1], arg.l1[2]};
            Translation x=arg.l2[0], y=arg.l2[1], z=arg.l2[2];
            OctTreeTPtr tree = OctTreeT::find_down(in1.tree(),n,x,y,z);
            if (!tree) MADNESS_EXCEPTION("TaskLeaf: SOL looking for tree?",0);
            taskq.add_local_hp(new Derived(tree, in1, in2, out1, axis, n, l));            
        };
        
        void forward(const OctTreeTPtr& t) const {
            Payload arg;
            arg.in1 = in1.ind;
            arg.in2 = in2.ind;
            arg.out1 = out1.ind;
            arg.n1 = n;
            arg.l1[0]=l[0]; arg.l1[1]=l[1]; arg.l1[2]=l[2];
            arg.n2 = tree->n(); 
            arg.l2[0]= tree->x(); arg.l2[1]= tree->y(); arg.l2[2]= tree->z();
        
            madness::comm_default->am_send(t->rank(), Derived::forward_task_handler, (AMArg&) arg);
        };
        
        
    public:
        TaskLeaf(OctTreeTPtr& tree,
                 Function<T>& in1, Function<T>& in2, Function<T>& out1,
                 int axis, Level n, const Translation *l) 
        : tree(tree), in1(in1), in2(in2), out1(out1), axis(axis), n(n)         
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
                //madness::print("APPLYING",tree->n(),tree->x(),tree->y(),tree->z(),(void*) in1.coeff(tree));
                apply();
                //madness::print("APPLYING DONE");
            }
            else {
                //madness::print("RECURRING",tree->n(),tree->x(),tree->y(),tree->z(),(void*) in1.coeff(tree));
                if (in1.coeff(tree)) in1._recur_coeff_down2(tree.get(),true);
                if (in1.ind!=in2.ind  &&  in2.coeff(tree)) in2._recur_coeff_down2(tree.get(),true);
                FOREACH_CHILD(OctTreeTPtr, tree, 
                              if (tree->islocal()) 
                                  taskq.add_local(new Derived(child, in1, in2, out1, axis, n, l));
                              else 
                                  forward(child););
                //madness::print("RECURRING DONE",tree->n(),tree->x(),tree->y(),tree->z(),(void*) in1.coeff(tree));
            };           
        };
        
        /// If in.coeff(tree) spawns task, otherwise recurs down active, local children 
        static void generate_tasks1(OctTreeTPtr& tree,
                                    Function<T>& in, Function<T>& out,
                                    int axis=0, Level n=0, const Translation *l=0) {
            MADNESS_ASSERT(tree);
            if (in.coeff(tree)) {
                taskq.add_local(new Derived(tree, in, in, out, axis, n, l));
            }
            else if (in.isactive(tree)) {
                out.set_active(tree);
                FOREACH_CHILD(OctTreeTPtr, tree, 
                              if (in.isactive(child) && in.islocal(child))
                                  generate_tasks1(child,in,out,axis,n,l););
            }    
        };                      
    };
    
    template <typename T>
    class TaskDiff : public TaskLeaf< T, TaskDiff<T> > {
    private:
        typedef TaskLeaf< T, TaskDiff<T> > baseT;
        SAV< Tensor<T> > left;
        SAV< Tensor<T> > right;
        //const Tensor<T>* left;
        //const Tensor<T>* right;

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
    Function<T> Function<T>::diff2(int axis) {
        reconstruct();
        Function<T> df = FunctionFactory<T>().k(k).compress(false).empty();
        madness::print("STARTING TASK GEN");
        TaskDiff<T>::generate_tasks1(tree(), *this, df, axis);
        madness::print("FINISHED TASK GEN");
        taskq.global_fence();
        madness::print("FINISHED TASK FENCE");
        _auto_clean(tree()); // Could do this on the way back up.
        return df;
    };
    
    
    
    
    
    template Function<double> Function<double>::diff2(int axis);
    template Function<double_complex> Function<double_complex>::diff2(int axis);
   
}
