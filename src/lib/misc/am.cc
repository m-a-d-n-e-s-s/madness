#include <iostream>
using std::cout;
using std::endl;

#include <misc/communicator.h>
using madness::Communicator;
using madness::AMArg;

#include <misc/print.h>

#include <misc/misc.h>
using madness::redirectio;

static unsigned long next = 1;

int myrand(void) {
    next = next * 1103515245 + 12345;
    return((unsigned)(next/65536) % 32768);
}

void mysrand(unsigned seed) {
    next = seed;
}

namespace madness {
    void task_add_am(am_handlerT op, ProcessID src, const AMArg& arg) {};
}

ProcessID random_proc(const Communicator& comm) {
    while (1) {
        ProcessID p = myrand()%comm.size();
        if (p != comm.rank()) return p;
    }
}

void handler(Communicator& comm, ProcessID from, const AMArg& arg) {
    comm.Send(arg.arg0+1, from, 33);
}

void hello(Communicator& comm, ProcessID from, const AMArg& arg) {
    madness::print(comm.rank(),"got hello from",from);
    cout.flush();
}

int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    Communicator comm;
    redirectio(comm);
    ProcessID me = comm.rank();
    long nproc = comm.nproc();
    //comm.set_debug(true);

    if (nproc == 1) throw "Gimme someone to talk to!";

    mysrand(me);
    
    //madness::xterm_debug(comm,0,0);

    try {
        comm.am_register(&handler);
        comm.am_register(&hello);
        comm.print_handlers();
        madness::print("ENTERING FENCE 0");
        comm.am_global_fence_spmd();
        madness::print("EXITING FENCE 0");

        for (int i=0; i<20; i++) {
            long reply = -1;
            ProcessID p = random_proc(comm);
            comm.am_send_recv(p,handler,AMArg(me+1000),&reply,sizeof(reply),p,33);
            if (reply != me+1001) throw "Ooops ...";
            comm.am_poll();
        }
        madness::print("ENTERING FENCE 1");
        cout.flush();
        comm.am_global_fence_spmd();
        madness::print("EXITING FENCE 1");
        cout.flush();
        comm.am_broadcast(hello,AMArg());  // Everyone says hello to everyone else
        madness::print("ENTERING FENCE 2");
        cout.flush();
        comm.am_global_fence_spmd();
        madness::print("EXITING FENCE 2");
        cout.flush();
    } catch (MPI::Exception e) {
        cout << " caught an MPI exception " << endl;
    } catch (const char* s) {
        cout << " caught a string exception " << s << endl;
    }

    if (me == 0) cout << "AM seems to be working" << endl;

    comm.close();
    MPI::Finalize();
    return 0;
}
