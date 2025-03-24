#include <madness/madness_config.h>

#ifdef HAVE_PARSEC
#include "parsec.h"
#include <iostream>
#include "thread.h"
#include <parsec/parsec_internal.h>

#include <cinttypes>

// Here we initialize with the right child class
namespace madness {
#ifndef MADNESS_ASSERTIONS_DISABLE
    static thread_local bool within_madness_task = false;
#endif

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
#ifndef MADNESS_ASSERTIONS_DISABLE
        within_madness_task = true;
#endif
        c->run(TaskThreadEnv(1, 0, 0));
#ifndef MADNESS_ASSERTIONS_DISABLE
        within_madness_task = false;
#endif
        return PARSEC_HOOK_RETURN_DONE;
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

    static parsec_key_t madness_parsec_make_key(const parsec_taskpool_t *tp, const parsec_assignment_t *assignments) {
        return (parsec_key_t)( reinterpret_cast<intptr_t> ( *(const void**)assignments ) );
    }

    static int madness_parsec_key_equal(parsec_key_t a, parsec_key_t b, void *user_data) {
        return a == b;
    }

    static char *madness_parsec_key_print(char *buffer, size_t buffer_size, parsec_key_t k, void *user_data) {
        snprintf(buffer, buffer_size, "MADNESS_TASK(%" PRIxPTR ")", k);
        return buffer;
    }

    static uint64_t madness_parsec_key_hash(parsec_key_t k, void *user_data) {
        return (uint64_t)k;
    }

    static parsec_key_fn_t madness_parsec_key_fns = {
        .key_equal = madness_parsec_key_equal,
        .key_print = madness_parsec_key_print,
        .key_hash  = madness_parsec_key_hash
    };

    const parsec_task_class_t madness_parsec_tc = {
            .name = (char*)"RUN",
            .flags = PARSEC_HAS_IN_IN_DEPENDENCIES | PARSEC_USE_DEPS_MASK,
            .task_class_id = 0,
            .nb_flows = 0,
            .nb_parameters = 0,
            .nb_locals = 0,
            .task_class_type = 0,
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
            .key_functions = &madness_parsec_key_fns,
            .make_key = madness_parsec_make_key,
            .task_snprintf = NULL,
#if defined(PARSEC_PROF_TRACE)
            .profile_info = NULL,
#endif
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

    static int parsec_madness_taskpool_enter_wait(parsec_taskpool_t* tp, void*_)
    {
        assert(tp != NULL);
        assert(NULL != tp->tdm.module);
        (void)_;
        tp->tdm.module->taskpool_ready(tp);
        return PARSEC_SUCCESS;
    }

    static int parsec_madness_taskpool_leave_wait(parsec_taskpool_t* tp, void*_b)
    {
        assert(tp != NULL);
        assert(NULL != tp->tdm.module);

        bool restart = *(bool*)_b;

        if(!restart) {
            /* This is the final taskpool_wait: we don't want to restart the termination
             * detection and re-attach to the context. */
            return PARSEC_SUCCESS;
        }

        /* Reset termination detector, so we can start adding tasks again */
        tp->tdm.module->monitor_taskpool(tp, parsec_taskpool_termination_detected);
        tp->tdm.module->taskpool_set_nb_tasks(tp, 0);
        tp->tdm.module->taskpool_set_runtime_actions(tp, 0); 

        /* We are re-attached to the context */
        parsec_atomic_fetch_inc_int32(&tp->context->active_taskpools);
        return PARSEC_SUCCESS;
    }

    const parsec_task_class_t* madness_parsec_tc_array[]= {&(madness::madness_parsec_tc), NULL};
    parsec_taskpool_t *ParsecRuntime::tp = nullptr;
    parsec_context_t *ParsecRuntime::ctx = nullptr;
    std::optional<bool> ParsecRuntime::made_new_ctx{};
    parsec_execution_stream_t *ParsecRuntime::madness_comm_thread_es = nullptr;
    bool ParsecRuntime::parsec_restart_taskpool = true;
#ifdef PARSEC_PROF_TRACE
    int ParsecRuntime::taskpool_profiling_array[2];
#endif
    ParsecRuntime::ParsecRuntime(int nb_threads) {
        if (ctx == nullptr) {
          /* Scheduler init*/
          int argc = 1;
          char *argv[2] = {(char *)"madness-app", nullptr};
          char **pargv = (char **)argv;
          ctx = parsec_init(nb_threads, &argc, &pargv);
          MPI_Comm parsec_comm = MPI_COMM_SELF;
          parsec_remote_dep_set_ctx(ctx, (intptr_t)parsec_comm);
          made_new_ctx = true;
        }
        else {
          made_new_ctx = false;
        }
        parsec_restart_taskpool = true;

        tp = PARSEC_OBJ_NEW(parsec_taskpool_t);
        tp->taskpool_name = strdup("MADNESS taskpool");
        tp->devices_index_mask = PARSEC_DEVICES_ALL;
        tp->nb_task_classes = 1;
        tp->task_classes_array = madness_parsec_tc_array;
        tp->update_nb_runtime_task = madness_parsec_update_runtime_nb_tasks;
        tp->on_enter_wait = parsec_madness_taskpool_enter_wait;
        tp->on_leave_wait = parsec_madness_taskpool_leave_wait;
        tp->on_leave_wait_data = &parsec_restart_taskpool;

        parsec_termdet_open_module(tp, (char*)"local");
        tp->tdm.module->monitor_taskpool(tp, parsec_taskpool_termination_detected);
        tp->tdm.module->taskpool_set_nb_tasks(tp, 0);
        tp->tdm.module->taskpool_set_runtime_actions(tp, 0); 
#ifdef PARSEC_PROF_TRACE
        tp->profiling_array = ParsecRuntime::taskpool_profiling_array;
        parsec_profiling_add_dictionary_keyword("MADNESS TASK", "fill:CC2828", 0, "",
                                                (int *)&tp->profiling_array[0],
                                                (int *)&tp->profiling_array[1]);
#endif
        if( 0 != parsec_context_add_taskpool(ctx, tp) ) {
            std::cerr << "ERROR!!! parsec_context_add_taskpool failed!!" << std::endl;
            MADNESS_EXCEPTION("ParsecRuntime::ParsecRuntime(): MADNESS taskpool could not be added to the PaRSEC context", 0);
        }
        if( 0 != parsec_context_start(ctx) ) {
            std::cerr << "ERROR!!! context_context_start failed!!" << std::endl;
            MADNESS_EXCEPTION("ParsecRuntime::ParsecRuntime(): PaRSEC context could not be started", 0);
        }
    }

    void ParsecRuntime::initialize_with_existing_context(parsec_context_t* ctx) {
        ParsecRuntime::ctx = ctx;
    }


    ParsecRuntime::~ParsecRuntime() {
        parsec_restart_taskpool = false;
        parsec_taskpool_wait(tp);
        parsec_taskpool_free(tp);
        assert(made_new_ctx.has_value());
        if (*made_new_ctx) {
            parsec_context_wait(ctx);
            parsec_fini(&ctx);
            ctx = nullptr;
        }
        if (nullptr != madness_comm_thread_es) {
              /* madness_comm_thread_es is just a copy of ES[0]. Resources (including es->profiling_es) are
             * actually freed during parsec_fini. Just need to free memory allocated to store it. */
              free(madness_comm_thread_es);
              madness_comm_thread_es = nullptr;
        }
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
        parsec_task->taskpool   = tp;
        parsec_task->task_class = &madness_parsec_tc;
        parsec_task->chore_mask = 1<<0;
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
        parsec_task_t* vp_parsec_task[1] = {parsec_task};
        PARSEC_LIST_ITEM_SINGLETON(parsec_task);
        tp->tdm.module->taskpool_addto_nb_tasks(tp, 1);

        parsec_execution_stream_t *es = parsec_my_execution_stream();
        if(nullptr == es) {
#ifndef MADNESS_ASSERTIONS_DISABLE
            if(nullptr != madness_comm_thread_es) {
                std::cerr << "!!! ERROR: the PaRSEC task backend expected at most one MADNESS communication thread.\n";
                MADNESS_EXCEPTION("ParsecRuntime::schedule(): more than one MADNESS communication thread", 0);
            }
#endif
            /* We are the MADNESS comm thread, which doesn't have a parsec execution stream.
             * Let's assign one... It's mostly a copy of 0, except for the profiling if this is enabled */
            madness_comm_thread_es = (parsec_execution_stream_t *)calloc(1, sizeof(parsec_execution_stream_t));
            memcpy(madness_comm_thread_es, context()->virtual_processes[0]->execution_streams[0], sizeof(parsec_execution_stream_t));
            /* This thread doesn't get a scheduler object: schedule happens on computation threads */
            madness_comm_thread_es->scheduler_object = nullptr;
#if defined(PARSEC_PROF_TRACE)
            madness_comm_thread_es->es_profile = parsec_profiling_stream_init(2*1024*1024, "MADNESS Comm thread");
            parsec_profiling_set_default_thread(madness_comm_thread_es->es_profile);
#endif /* PARSEC_PROF_TRACE */
            parsec_set_my_execution_stream(madness_comm_thread_es);
            es = madness_comm_thread_es;
        }

        __parsec_schedule_vp(es, vp_parsec_task, is_main_thread() ? 1 : 0);
    }

    void ParsecRuntime::wait() {
#ifndef MADNESS_ASSERTIONS_DISABLE
        if(within_madness_task) {
            std::cerr << "!!! ERROR: the PaRSEC task backend does not support recursive calls to wait.\n"
                      << "!!! ERROR: this is usually a symptom that a task is forcing a future that is not ready.\n";
            MADNESS_EXCEPTION("ParsecRuntime::wait(): recursive call to wait/test detected", 0);
        }
#endif
        parsec_taskpool_wait(tp);
    }

    int ParsecRuntime::test() {
#ifndef MADNESS_ASSERTIONS_DISABLE
        if(within_madness_task) {
            std::cerr << "!!! ERROR: the PaRSEC task backend does not support recursive calls to wait.\n"
                      << "!!! ERROR: this is usually a symptom that a task is forcing a future that is not ready.\n";
            MADNESS_EXCEPTION("ParsecRuntime::wait(): recursive call to wait/test detected", 0);
        }
#endif
        int rc = parsec_taskpool_test(tp);
        assert(rc >= 0);
        return rc;
    }

    extern "C"{

#include <stdio.h>

    }
}
#endif
