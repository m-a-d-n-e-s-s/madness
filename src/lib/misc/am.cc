#include <iostream>
using std::cout;
using std::endl;

#include <misc/communicator.h>
using madness::Communicator;
using madness::AMArg;

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

ProcessID random_proc(const Communicator& comm) {
    while (1) {
        ProcessID p = myrand()%comm.size();
        if (p != comm.rank()) return p;
    }
}

void handler(Communicator& comm, ProcessID from, const AMArg& arg) {
    comm.Send(arg.arg0+1, from, 33);
}

int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    Communicator comm;
    //redirectio(comm);
    ProcessID me = comm.rank();
    long nproc = comm.nproc();

    if (nproc == 1) throw "Gimme someone to talk to!";

    mysrand(me);

    try {
        int handle = comm.am_register(&handler);

        for (int i=0; i<20; i++) {
            long reply = -1;
            ProcessID p = random_proc(comm);
            comm.am_send_recv(AMArg(handle,me+1000),p,&reply,sizeof(reply),p,33);
            //MPI::Request req = comm.Irecv(reply, p, 33);
            //comm.am_send(AMArg(handle,me+1000),p);
            //comm.am_wait(req);
            if (reply != me+1001) throw "Ooops ...";
            comm.am_poll();
        }

        comm.am_barrier();
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
