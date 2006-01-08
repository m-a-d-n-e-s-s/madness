#include <iostream>
using std::cout;
using std::endl;

/// \file misc/test.cc

#include "print.h"
#include "communicator.h"
using namespace madness;

int main(int argc, char* argv[]) {
    PBEGIN_(argc, argv);
    Communicator comm;

    long hello=0;
    if (comm.rank() == 0)  hello = 1;
    comm.bcast(hello,0);

    print(comm.rank(),": hello=", hello);

    // Everyone sends a message to everyone else
    char msg[] = "Greetings";
    for (ProcessID from=0; from<comm.nproc(); from++) {
        for (ProcessID to=0; to<comm.nproc(); to++) {
            if (to != from) {
                char buf[256];
                buf[0] = 0;
                if (comm.rank() == from) {
                    comm.send(msg, sizeof(msg), to, 1);
                }
                else if (comm.rank() == to) {
                    comm.recv(buf, sizeof(buf), from, 1);
                    print(to,"received",buf,"from",from);
                }
            }
        }
    }

    // Accumulate a sum around a ring
    double sum[] = {0.0,1.0,0.0};
    if (comm.rank() == 0) {
        comm.send(sum, 3, 1, 99);
        comm.recv(sum, 3, comm.nproc()-1, 99);
        print("the final sum is",sum[0],sum[1],sum[2]);
    }
    else {
        comm.recv(sum, 3, comm.rank()-1, 99);
        sum[1]++;
        comm.send(sum, 3, (comm.rank()+1)%comm.nproc(), 99);
    }

    PEND_();

    return 0;
}
