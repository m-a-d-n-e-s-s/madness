#ifndef TASK_LEAF_H
#define TASK_LEAF_H

namespace madness {


    struct true_predicate {
        bool operator()(const FunctionImpl<T,NDIM>::keyT& key) const {return true;};
    };
        

    template <typename T, 
              int NDIM, 
              typename memfunT, 
              typename select_predicateT,
              typename apply_predicateT=predicate_true>
    class TaskGenerator :  WorldObject< TaskGenerator<T,NDIM,memfunT,select_predicateT,apply_predicateT> > {
    protected:
        typedef FunctionImpl<T,NDIM> implT;
        typedef implT::keyT keyT;
        implT* fa;
        World& world; 
        memfunT memfun;
        select_predicateT select_predicate;
        apply_predicateT apply_predicate;
        
        void recur_down(const keyT& key) {
            MADNESS_EXCEPTION("No recur down yet",0);
        };
        
        // Conditionally applies operator on LOCAL data or recurs
        // down and forwards task to owner of child.
        void conditional_apply(const keyT& key) {
            if (predicate(key)) {
                task(world.rank(),apply,fa->find(key));
            }
            else {
                foreach_child(key,recur_down);
            }
        };
        
        
    public:
        TaskGenerator(FunctionImpl<T,NDIM>* fa, opT& op, 
                      const select_predicateT& select_predicate,
                      const apply_predicateT& apply_predicate)
            : WorldObject< TaskGenerator<T,NDIM> >(fa->world)
              , fa(fa)
              , world(fa->world)
              , memfun(memfun)
              , select_predicate(select_predicate) 
              , apply_predicate(apply_predicate)
        {
            // Pre-generate list of leaves since recur down will
            // modify container and invalidate iterators.  To avoid same
            // issue and/or double application of operations, remote recur
            // operations will send messages to recur_down_handler which
            // will be handled after the local list has been generated.
            std::vector<keyT> leaves = fa->leaf_nodes();
            std::for_each(leaves.begin(),leaves.end(),conditional_apply);
            process_pending();
        };
        
        void apply(const iterator& it) {
            (fa->*memfun)(it->first,it->second);
        };
    };

    template <typename T, 
              int NDIM, 
              typename memfunT, 
              typename select_predicateT,
              typename apply_predicateT=predicate_true>
    TaskGenerator<T,NDIM,memfunT,select_predicateT,apply_predicateT> 
    task_generator(FunctionImpl<T,NDIM>* fa, memfunT memfun,
                   const select_predicateT& select_predicate,
                   const apply_predicateT& apply_predicate) {
        return TaskGenerator<T,NDIM,memfunT,select_predicateT,apply_predicateT> (fa,memfun,select_predicate,apply_predicate);
    };
}




        










#endif
