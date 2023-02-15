#include <madness/madness_config.h>

#ifdef HAVE_PARSEC
#include "parsec.h"
#include <iostream>
#include "thread.h"
#include <parsec/parsec_internal.h>

// Here we initialize with the right child class
namespace madness {
    parsec_hook_return_t complete_madness_task_execution (parsec_execution_stream_t *es,
                                                          parsec_task_t *task) {
        (void)es;
        task->taskpool->tdm.module->taskpool_addto_nb_tasks(task->taskpool, -1);
        return PARSEC_HOOK_RETURN_DONE;
    }

    parsec_hook_return_t release_madness_task (parsec_execution_stream_t *es, 
                                               parsec_task_t *task) {
        PoolTaskInterface *c = ((PoolTaskInterface **)task->locals)[0];
        assert(c);
        (void)es;
        delete(c);
        return PARSEC_HOOK_RETURN_DONE;
    }

    parsec_hook_return_t run_madness_task(parsec_execution_stream_t *eu,
                                          parsec_task_t *task) {
        PoolTaskInterface *c = ((PoolTaskInterface **)task->locals)[0];
        assert(c);
        c->run(TaskThreadEnv(1, 0, 0));
        return PARSEC_HOOK_RETURN_DONE;
    }

    static uint64_t RUN_hash(const parsec_taskpool_t *tp, const parsec_assignment_t *a) {
        (void)a; (void)tp;
        return 0;
    }

    static parsec_hook_return_t empty_hook (parsec_execution_stream_t *es,
                                            parsec_task_t *task) {
        (void)es; (void)task;
        return PARSEC_HOOK_RETURN_DONE;
    }

    static const __parsec_chore_t __RUN_chores[] = {
            {.type = PARSEC_DEV_CPU, .evaluate = NULL, .hook = run_madness_task},
            {.type = PARSEC_DEV_NONE, .evaluate = NULL, .hook = (parsec_hook_t *) NULL},
    };

    static int madness_parsec_update_runtime_nb_tasks(parsec_taskpool_t *tp, int32_t nb_tasks) {
         return tp->tdm.module->taskpool_addto_nb_tasks(tp, nb_tasks);
    }

    const parsec_task_class_t madness_parsec_tc = {
            .name = (char*)"RUN",
            .flags = PARSEC_HAS_IN_IN_DEPENDENCIES | PARSEC_USE_DEPS_MASK,
            .task_class_id = 0,
            .nb_flows = 0,
            .nb_parameters = 0,
            .nb_locals = 0,
            .dependencies_goal = 0,
            .params = {NULL,NULL,NULL, NULL},
            .locals = {NULL,NULL,NULL, NULL},
            .in = {NULL, NULL, NULL, NULL},
            .out = {NULL, NULL, NULL, NULL},
            .priority = NULL,
            .properties = NULL,
            .initial_data = NULL,
            .final_data = NULL,
            .data_affinity = NULL,
            .key_generator = (parsec_functionkey_fn_t *) RUN_hash,
            .key_functions = NULL,
            .make_key = NULL,
#if defined(PARSEC_SIM)
            .sim_cost_fct = (parsec_sim_cost_fct_t*) NULL,
#endif
            .get_datatype = (parsec_datatype_lookup_t *) NULL,
            .prepare_input = empty_hook,
            .incarnations = __RUN_chores,
            .prepare_output = (parsec_hook_t *) NULL,

            .find_deps = (parsec_find_dependency_fn_t*)NULL,
            .update_deps = NULL,

            .iterate_successors = (parsec_traverse_function_t *) NULL,
            .iterate_predecessors = (parsec_traverse_function_t *) NULL,
            .release_deps = (parsec_release_deps_t *) NULL,
            .complete_execution = complete_madness_task_execution,
            .new_task = (parsec_new_task_function_t*) NULL,
            .release_task = release_madness_task, //object delete,
            .fini = (parsec_hook_t *) NULL,
    };

    const parsec_task_class_t* madness_parsec_tc_array[]= {&(madness::madness_parsec_tc), NULL};
    parsec_taskpool_t ParsecRuntime::tp = {
            .super = { 0x0, },
            .taskpool_id = 0,
            .taskpool_name = (char*)"MADNESS taskpool",
            .nb_tasks = 0,
            .taskpool_type = 0,
            .devices_index_mask = PARSEC_DEVICES_ALL,
            .nb_task_classes = 1,
            .priority = 0,
            .nb_pending_actions = 0,
            .context = NULL,
            .tdm = {0x0, },
            .startup_hook = NULL,
            .task_classes_array = madness_parsec_tc_array,
#if defined(PARSEC_PROF_TRACE)
            .profiling_array = NULL,
#endif
            .on_enqueue = NULL,
            .on_enqueue_data = NULL,
            .on_complete = NULL,
            .on_complete_data = NULL,
            .update_nb_runtime_task = madness_parsec_update_runtime_nb_tasks,
            .destructor = NULL,
            .dependencies_array = NULL,
            .repo_array = NULL
    };
    parsec_context_t *ParsecRuntime::ctx = nullptr;
#ifdef PARSEC_PROF_TRACE
    int ParsecRuntime::taskpool_profiling_array[2];
#endif
    ParsecRuntime::ParsecRuntime(int nb_threads) {
        assert(ctx == nullptr);
        /* Scheduler init*/
        int argc = 1;
        char *argv[2] = { (char *)"madness-app", nullptr };
        char **pargv = (char **)argv;
        ctx = parsec_init(nb_threads, &argc, &pargv);
        MPI_Comm parsec_comm  = MPI_COMM_SELF;
        parsec_remote_dep_set_ctx(ctx, (intptr_t)parsec_comm);
#ifdef PARSEC_PROF_TRACE
        ParsecRuntime::tp.profiling_array = ParsecRuntime::taskpool_profiling_array;
        parsec_profiling_add_dictionary_keyword("MADNESS TASK", "fill:CC2828", 0, "",
                                                (int *)&ParsecRuntime::tp.profiling_array[0],
                                                (int *)&ParsecRuntime::tp.profiling_array[1]);
#endif
        if( 0 != parsec_context_add_taskpool(ctx, &tp) ) {
            std::cerr << "ERROR: parsec_context_add_taskpool failed!!" << std::endl;
        }
        parsec_taskpool_update_runtime_nbtask(&tp, 1);
        if( 0 != parsec_context_start(ctx) ) {
            std::cerr << "ERROR: context_context_start failed!!" << std::endl;
        }
    }

    ParsecRuntime::~ParsecRuntime() {
        parsec_fini(&ctx);
        ctx = nullptr;
    }

    parsec_context_t *ParsecRuntime::context() { return ctx; }

    parsec_execution_stream_t *ParsecRuntime::execution_stream() {
        parsec_execution_stream_t *my_es;
        my_es = parsec_my_execution_stream();
        if(nullptr == my_es) {
            my_es = context()->virtual_processes[0]->execution_streams[0];
        }
        assert(nullptr != my_es);
        return my_es;
    }

    parsec_task_t *ParsecRuntime::task(bool is_high_priority, void *ptr) {
        parsec_execution_stream_t *my_es = execution_stream();
        parsec_task_t* parsec_task = static_cast<parsec_task_t*>( parsec_thread_mempool_allocate( my_es->context_mempool ) );
        PARSEC_LIST_ITEM_SINGLETON(parsec_task);
        parsec_task->taskpool   = &tp;
        parsec_task->task_class = &madness_parsec_tc;
        parsec_task->chore_id   = 0;
        parsec_task->status     = PARSEC_TASK_STATUS_NONE;
        parsec_task->priority   = is_high_priority ? 1000 : 0; // 1 & 0 would work as good
        ((void **)parsec_task->locals)[0] = ptr;
        return parsec_task;
    }

    void ParsecRuntime::delete_parsec_task(parsec_task_t *task) {
        parsec_thread_mempool_free( task->mempool_owner, task );
    }

    void ParsecRuntime::schedule(PoolTaskInterface *task) {
        parsec_task_t* parsec_task = task->parsec_task;
        PARSEC_LIST_ITEM_SINGLETON(parsec_task);
        tp.tdm.module->taskpool_addto_nb_tasks(&tp, 1);
        __parsec_schedule(execution_stream(), parsec_task, 0);
    }

    void ParsecRuntime::wait() {
        parsec_taskpool_update_runtime_nbtask(&tp, -1);
        parsec_context_wait(ctx);
    }

    int ParsecRuntime::test() {
        int rc = parsec_taskpool_test(&tp);
        assert(rc >= 0);
        return rc;
    }

    extern "C"{

#include <stdio.h>

    }
}
#endif
