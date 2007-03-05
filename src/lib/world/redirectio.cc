#include <world/world.h>
#include <iostream>
using std::cout;
using std::endl;
#include <fstream>
#include <cstdio>

static std::ofstream fout;
namespace madness {
    void redirectio(World& world) {
        if (world.mpi.rank() != 0) {
            char filename[256];
            std::sprintf(filename,"log.%5.5d",world.mpi.rank());
            
            // Was using this but it seems as if it does not
            // redirect C stdio interface (at least on Cygwin)
            //fout.open(filename);
            //cout.rdbuf(fout.rdbuf());
            //std::cerr.rdbuf(fout.rdbuf());
            
            freopen(filename, "w", stdout);
            freopen(filename, "w", stderr);
        }
    }
}




