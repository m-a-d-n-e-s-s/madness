/// \file mra/startup.cc

#include <mra/mra.h>

namespace madness {
    void startup(World& world, int argc, char** argv) {
        redirectio(world);
        
        if (world.rank() == 0){
            print("The processor frequency is",cpu_frequency());
            print("there are",world.mpi.nproc(),"processes and I am process",world.mpi.rank());
        }

        for (int arg=1; arg<argc; arg++) {
            if (std::strcmp(argv[arg],"-dx")==0) 
                xterm_debug("world", 0);
            else if (std::strcmp(argv[arg],"-dn") ==0 && 
                     std::atoi(argv[arg+1])==world.rank()) 
                xterm_debug("world",0);
            else if (std::strcmp(argv[arg],"-dam")==0) 
                world.am.set_debug(true);
            else if (std::strcmp(argv[arg],"-dmpi")==0) 
                world.mpi.set_debug(true);
        }

        world.gop.fence();

        load_coeffs(world);
        load_quadrature(world);
        MADNESS_ASSERT(gauss_legendre_test());
        MADNESS_ASSERT(test_two_scale_coefficients());
        
        world.gop.fence();
    }
}
