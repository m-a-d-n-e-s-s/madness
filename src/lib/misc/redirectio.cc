#include <iostream>
using std::cout;
using std::endl;

#include <fstream>

#include <cstdio>

#include <misc/communicator.h>

static std::ofstream fout;
namespace madness {
    void redirectio(const Communicator& comm) {
        if (comm.rank() != 0) {
            char filename[256];
            std::sprintf(filename,"log.%5.5d",comm.rank());
            
            // Was using this but it seems as if it does not
            // redirect C stdio interface (at least on Cygwin)
            //fout.open(filename);
            //cout.rdbuf(fout.rdbuf());
            //std::cerr.rdbuf(fout.rdbuf());
            
            freopen(filename, "w", stdout);
        }
    }
}




