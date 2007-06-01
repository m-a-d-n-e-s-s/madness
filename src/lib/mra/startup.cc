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

        std::cout << std::boolalpha;  // Pretty printing of booleans

#ifdef FUNCTION_INSTANTIATE_1
        FunctionDefaults<1>::set_defaults();
#endif
#ifdef FUNCTION_INSTANTIATE_2
        FunctionDefaults<2>::set_defaults();
#endif
#ifdef FUNCTION_INSTANTIATE_3
        FunctionDefaults<3>::set_defaults();
#endif
#ifdef FUNCTION_INSTANTIATE_4
        FunctionDefaults<4>::set_defaults();
#endif
#ifdef FUNCTION_INSTANTIATE_5
        FunctionDefaults<5>::set_defaults();
#endif
#ifdef FUNCTION_INSTANTIATE_6
        FunctionDefaults<6>::set_defaults();
#endif


        print("loading coeffs, etc.");

        load_coeffs(world);
        load_quadrature(world);

        print("testing coeffs, etc.");
        MADNESS_ASSERT(gauss_legendre_test());
        MADNESS_ASSERT(test_two_scale_coefficients());

        print("done with startup");

        world.gop.fence();
    }
}
