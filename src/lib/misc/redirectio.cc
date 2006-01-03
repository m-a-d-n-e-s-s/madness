#include <iostream>
using std::cout;
using std::endl;

#include <fstream>

#include "communicator.h"

static std::ofstream fout;
namespace madness {
    void redirectio(const Communicator& comm) {
        if (comm.rank() != 0) {
            char filename[256];
            std::sprintf(filename,"log.%5.5ld",comm.rank());
            fout.open(filename);
            cout.rdbuf(fout.rdbuf());
            std::cerr.rdbuf(fout.rdbuf());
        }
    }
}




