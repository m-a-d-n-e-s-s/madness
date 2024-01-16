#ifndef MADNESS_PARSEC_INCLUED
#define MADNESS_PARSEC_INCLUED

#include <madness/madness_config.h>

#ifdef HAVE_PARSEC

#include <parsec.h>
#include <parsec/parsec_config.h>
#include <parsec/parsec_internal.h>
#include <parsec/mca/device/device.h>
#include <parsec/execution_stream.h>
#include <parsec/scheduling.h>

#include <iostream>
#include <optional>

namespace madness{

    class PoolTaskInterface;

    class ParsecRuntime {
    private:
        static parsec_context_t *ctx;
        static std::optional<bool> made_new_ctx;
        static parsec_taskpool_t *tp;
        static bool parsec_restart_taskpool;
        static parsec_execution_stream_t *madness_comm_thread_es;
#ifdef PARSEC_PROF_TRACE
        static int               taskpool_profiling_array[2];
#endif /* PARSEC_PROF_TRACE */

    public:
        ParsecRuntime(int nb_threads);
        ~ParsecRuntime();

        static parsec_context_t* context();
        static void initialize_with_existing_context(parsec_context_t* ctx);
        static parsec_execution_stream_t *execution_stream();
        static void schedule(PoolTaskInterface* task);
        static int test();
        static void wait();
        static parsec_task_t *task(bool is_high_priority, void *ptr);
        static void delete_parsec_task(parsec_task_t *t);
    };
}

#endif // HAVE_PARSEC

#endif // MADNESS_PARSEC_INCLUED
