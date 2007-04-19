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
        
        // Conditionally applies operator on LOCAL data or recurs
        // down and forwards task to owner of child.
        void conditional_apply(const keyT& key) {
            if (apply_predicate(key)) {
                task(world.rank(),apply,fa->find(key));
            }
            else {
                foreach_child(key,recur_down);
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
            // was for_each but as usual std::bind1st cannot cut it ... waiting
            // for boost and/or tr1
            for (typename std::vector<keyT>::iterator it=keys.begin(); it != keys.end(); ++it) {
                this->conditional_apply(*it);
            };
            this->process_pending();
        };
        
        void apply(const typename implT::iterator& it) {
            (fa->*memfun)(it->first,it->second);
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
