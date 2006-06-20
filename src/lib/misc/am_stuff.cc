

#include <misc/communicator.h>

namespace madness {
    static long _am_nchild_registered=0;

    Communicator* comm_default;

    void am_barrier_handler(Communicator& comm, ProcessID proc, const AMArg& arg) {
        _am_nchild_registered++;
    };

    long am_barrier_nchild_registered() {
        return _am_nchild_registered;
    };

    void am_barrier_zero_nchild_registered() {
        _am_nchild_registered = 0;
    };
}


