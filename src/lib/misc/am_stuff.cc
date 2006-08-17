

#include <misc/communicator.h>

namespace madness {
    Communicator* comm_default;

    void am_ndiff_handler(Communicator& comm, ProcessID proc, const AMArg& arg) {
        comm_default->am_ndiff_single_threaded(true,bind_mem_fun(comm_default,&Communicator::am_wait),noop);
    };
}


