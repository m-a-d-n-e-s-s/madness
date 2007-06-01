#ifndef TASK_LEAF_H
#define TASK_LEAF_H

namespace madness {

    template <typename T, int NDIM, typename Pmap> class FunctionImpl;

    template <typename implT,
              typename keyT,
              typename memfunT, 
              typename select_predicateT,
              typename apply_predicateT>
    class TaskGenerator :  public WorldObject< TaskGenerator<implT,keyT,memfunT,select_predicateT,apply_predicateT> > {
    protected:
        implT* fa;        
        World& world; 
        memfunT memfun;
        apply_predicateT apply_predicate;
        
        void recur_down(const keyT& key) {
            MADNESS_EXCEPTION("No recur down yet",0);
        };
        
        // Applied to LOCAL node.  If apply predicate is true,
        // generates task to apply operator.  Otherwise generates task
        // to recur down.
        void conditional_apply(const datumT& d) {
            if (apply_predicate(d)) {
                task(world.rank(),d)
            }
            else {
                foreach_child(d.firstkey,recur_down);
            }
        };
        
        
    public:
        TaskGenerator(implT* fa, 
                      memfunT& memfun,
                      const select_predicateT& select_predicate,
                      const apply_predicateT& apply_predicate)
            : WorldObject< TaskGenerator<implT,keyT,memfunT,select_predicateT,apply_predicateT> >(fa->world)
            , fa(fa)
            , world(fa->world)
            , memfun(memfun)
            , apply_predicate(apply_predicate)
        {
            // Pre-generate list of leaves since recur down will
            // modify container and invalidate iterators.  To avoid same
            // issue and/or double application of operations, remote recur
            // operations will send messages to recur_down_handler which
            // will be handled after the local list has been generated.
            std::vector<keyT> keys = fa->keys(select_predicate);
            for (typename std::vector<keyT>::iterator it=keys.begin(); it != keys.end(); ++it) {
                this->conditional_apply(*it);
            };
            this->process_pending();
        };
        
        void apply(const datumT& d) {
            (fa->*memfun)(d.first,d.second);
        };
    };

    template <typename T, 
              int NDIM, 
              typename pmapT,
              typename memfunT, 
              typename select_predicateT,
              typename apply_predicateT>
    TaskGenerator<FunctionImpl<T,NDIM,pmapT>,Key<NDIM>,memfunT,select_predicateT,apply_predicateT> 
    task_generator(FunctionImpl<T,NDIM,pmapT>* fa, memfunT memfun,
                   const select_predicateT& select_predicate,
                   const apply_predicateT& apply_predicate) {
        return TaskGenerator<FunctionImpl<T,NDIM,pmapT>,Key<NDIM>,memfunT,select_predicateT,apply_predicateT> (fa,memfun,select_predicate,apply_predicate);
    };
}




        










#endif
