/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680


  $Id$
*/


/// \file mra/startup.cc

#include <madness/mra/mra.h>
#include <madness/tensor/tensor.h>
#include <madness/world/timers.h>
//#include <madness/mra/mraimpl.h> !!!!!!!!!!!!!!!!  NOOOOOOOOOOOOOOOOOOOOOOOOOO !!!!!!!!!!!!!!!!!!!!!!!
#include <iomanip>
#include <cstdlib>

namespace madness {


  // Quick check on the performance of both BLAS and fastest/slowest processor
  static void time_transform(World& world, int& mflopslo, int& mflopshi) {
    const long N=20, NLOOP=100;
    Tensor<double> s(N,N,N), r(N,N,N), w(N,N,N), c(N,N);
    s=1.0; r=1.0; c=1.0;
    for (long loop=0; loop<10; loop++) fast_transform(s, c, r, w);
    double used = cpu_time();
    for (long loop=0; loop<NLOOP; loop++) fast_transform(s, c, r, w);
    used = cpu_time()-used;
    mflopslo = mflopshi = int(1e-6*2.0*3.0*N*N*N*N*double(NLOOP)/used);
    world.gop.max(mflopshi);
    world.gop.min(mflopslo);
  }



    void startup(World& world, int argc, char** argv, bool doprint, bool make_stdcout_nice_to_reals) {
        const char* data_dir = MRA_DATA_DIR;

        // Process command line arguments
        for (int arg=1; arg<argc; ++arg) {
            if (strcmp(argv[arg],"-dx")==0)
                xterm_debug(argv[0], 0);
            else if (strcmp(argv[arg],"-lf")==0)
                redirectio(world);
            else if (strcmp(argv[arg],"-dn") ==0 &&
                     std::atoi(argv[arg+1])==world.rank())
                xterm_debug("world",0);
//             else if (strcmp(argv[arg],"-dam")==0)
//                 world.am.set_debug(true);
//             else if (strcmp(argv[arg],"-dmpi")==0)
//                 world.mpi.set_debug(true);
            else if (strcmp(argv[arg],"-rio")==0)
                redirectio(world);
        }

        // Process environment variables
        if (getenv("MRA_DATA_DIR")) data_dir = getenv("MRA_DATA_DIR");

        // Need to add an RC file ...

        world.gop.fence();

	init_tensor_lapack();

        // if allowed: configure std::cout to be useful for printing reals
        if (make_stdcout_nice_to_reals) {
          std::cout << std::boolalpha; // Pretty printing of booleans
          std::cout << std::scientific;
          std::cout << std::showpoint;
          // std::cout << std::showpos;
          std::cout << std::setprecision(6);
        }

#ifdef FUNCTION_INSTANTIATE_1
        FunctionDefaults<1>::set_defaults(world);
        Displacements<1> d1;
#endif
#ifdef FUNCTION_INSTANTIATE_2
        FunctionDefaults<2>::set_defaults(world);
        Displacements<2> d2;
#endif
#ifdef FUNCTION_INSTANTIATE_3
        FunctionDefaults<3>::set_defaults(world);
        Displacements<3> d3;
#endif
#ifdef FUNCTION_INSTANTIATE_4
        FunctionDefaults<4>::set_defaults(world);
        Displacements<4> d4;
#endif
#ifdef FUNCTION_INSTANTIATE_5
        FunctionDefaults<5>::set_defaults(world);
        Displacements<5> d5;
#endif
#ifdef FUNCTION_INSTANTIATE_6
        FunctionDefaults<6>::set_defaults(world);
        Displacements<6> d6;
#endif

        //if (world.rank() == 0) print("loading coeffs, etc.");

        load_coeffs(world, data_dir);

        //if (world.rank() == 0) print("loading quadrature, etc.");

        load_quadrature(world, data_dir);

        // This to init static data while single threaded
        initialize_legendre_stuff();

        //if (world.rank() == 0) print("testing coeffs, etc.");
        MADNESS_CHECK(gauss_legendre_test());
        MADNESS_CHECK(test_two_scale_coefficients());

	int mflopslo, mflopshi;
	time_transform(world, mflopslo, mflopshi);

        // print the configuration options
        if (doprint && world.rank() == 0) {
            print("");
            print("--------------------------------------------");
            print("   MADNESS",MADNESS_PACKAGE_VERSION, "multiresolution suite");
            print("--------------------------------------------");
            print("");
            print("   number of processors ...", world.size());
            print("    processor frequency ...", cpu_frequency());
            print("            host system ...", HOST_SYSTEM);
            print("          configured by ...", MADNESS_CONFIGURATION_USER);
            print("          configured on ...", MADNESS_CONFIGURATION_HOST);
            print("          configured at ...", MADNESS_CONFIGURATION_DATE);
            print("                    CXX ...", MADNESS_CONFIGURATION_CXX);
            print("               CXXFLAGS ...", MADNESS_CONFIGURATION_CXXFLAGS);
#ifdef OPTERON_TUNE
            print("             tuning for ...", "opteron");
#elif defined(CORE_DUO_TUNE)
            print("             tuning for ...", "core duo");
#else
            print("             tuning for ...", "default");
#endif
#ifdef BOUNDS_CHECKING
            print(" tensor bounds checking ...", "enabled");
#endif
#ifdef TENSOR_INSTANCE_COUNT
            print("  tensor instance count ...", "enabled");
#endif
#ifdef STUBOUTMPI
	    print("                    MPI ...", "stubbed out");
#elif MADNESS_MPI_THREAD_LEVEL == MPI_THREAD_SERIALIZED
	    print("                    MPI ...", "serialized");
#else
	    print("                    MPI ...", "multiple");
#endif
#if HAVE_INTEL_TBB
            print(" multi-threaded runtime ...", "Intel TBB");
#elif defined(HAVE_PARSEC)
            print(" multi-threaded runtime ...", "PaRSEC");
#else
            print(" multi-threaded runtime ...", "MADNESS ThreadPool");
#endif
#if HAVE_INTEL_MKL
            print("                   BLAS ...", "Intel MKL",  mflopslo, mflopshi, "MFLOP/s");
#elif defined(HAVE_FAST_BLAS)
            print("                   BLAS ...", "Fast BLAS not MKL",  mflopslo, mflopshi, "MFLOP/s");
#else
            print("                   BLAS ...", "Slow reference",  mflopslo, mflopshi, "MFLOP/s");
#endif
#ifdef HAVE_MTXMQ
            print("                   BLAS ...", "Have MTXMQ");
#endif
            print("               compiled ...",__TIME__," on ",__DATE__);

            //         print(" ");
            //         IndexIterator::test();
        }


        world.gop.fence();
    }


    std::string get_mra_data_dir() {
        const char* data_dir = MRA_DATA_DIR; // as defined at configure time
        if (getenv("MRA_DATA_DIR")) data_dir = getenv("MRA_DATA_DIR"); // override by environment variable
        return std::string(data_dir);
    };

}
