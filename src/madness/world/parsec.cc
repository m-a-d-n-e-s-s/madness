#include <madness/madness_config.h>

#ifdef HAVE_PARSEC
#include "parsec.h"
#include <iostream>
#include "thread.h"
#include <parsec/parsec_internal.h>

// Here we initialize with the right child class
namespace madness {
    parsec_hook_return_t complete_madness_task_execution (parsec_execution_stream_t *es,
                                                          parsec_task_t *task)
    {
        (void)es;
        task->taskpool->tdm.module->taskpool_addto_nb_tasks(task->taskpool, -1);
        return PARSEC_HOOK_RETURN_DONE;
    }

    parsec_hook_return_t release_madness_task (parsec_execution_stream_t *es, 
                                               parsec_task_t *task)
    {
        PoolTaskInterface *c = ((PoolTaskInterface **)task->locals)[0];
        (void)es;
        delete(c);
        return PARSEC_HOOK_RETURN_DONE;
    }

    const parsec_task_class_t* PoolTaskInterface::parsec_tc = &madness_parsec_tc;

    parsec_hook_return_t pointer_call(parsec_execution_stream_t *eu, 
                                      parsec_task_t *task)
    {
        PoolTaskInterface *c = ((PoolTaskInterface **)task->locals)[0];
        c->run(TaskThreadEnv(1, 0, 0));
        return PARSEC_HOOK_RETURN_DONE;
    }


  extern "C"{

#include <stdio.h>

        static uint64_t RUN_hash(const parsec_taskpool_t *tp, const parsec_assignment_t *a)
        {
            (void)a; (void)tp;
            return 0;
        }

        static parsec_hook_return_t empty_hook (parsec_execution_stream_t *es, 
                                                parsec_task_t *task)
        {
            (void)es; (void)task;
            return PARSEC_HOOK_RETURN_DONE;
        }

        static const __parsec_chore_t __RUN_chores[] = {
            {.type = PARSEC_DEV_CPU,
             .evaluate = NULL,
             .hook = pointer_call},
            {.type = PARSEC_DEV_NONE,
             .evaluate = NULL,
             .hook = (parsec_hook_t *) NULL},
        };


        const parsec_task_class_t madness_parsec_tc = {
            .name = "RUN",
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

        parsec_taskpool_t madness_parsec_tp = {
            .super = { 0x0, },
            .taskpool_id = 0,
            .taskpool_name = "MADNESS taskpool",
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
            .update_nb_runtime_task = NULL,
            .destructor = NULL,
            .dependencies_array = NULL,
            .repo_array = NULL
        };
    }
}
#endif
