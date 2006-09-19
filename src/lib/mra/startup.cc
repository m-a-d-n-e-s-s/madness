/// \file mra/startup.cc

#include <iostream>
using std::cout;
using std::endl;

#include <cstring>
using std::strcmp;

#include <mra/mra.h>
#include <misc/misc.h>
#include <misc/communicator.h>
#include <mra/twoscale.h>
#include <mra/legendre.h>
#include <tensor/tensor.h>
#include <misc/madexcept.h>
#include <tensor/mtrand.h>
#include <tasks/tasks.h>

namespace madness {
    
    void* FunctionDataPointersBase::p[FunctionNode::size];

    
    Communicator& startup(int argc, char** argv) {
        // The following should be used to setup all calculations
        // 1) Initialize parallel environment
        // 2) Setup communication information
        // 3) Redirect standard output+err for parallel processes
        // 4) Load coeffs and quadrature information from file
        // 5) Setup default OctTreeLayout
        // 6) Sanity check
        // 7) Top level catching of exceptions
        MPI::Init(argc, argv);
        static Communicator comm;
        madness::comm_default = &comm;
        redirectio(comm);
        comm.print();
        load_coeffs(comm);
        load_quadrature(comm);
        
        FunctionDefaults::tree = SharedPtr<FunctionOctTree>(new FunctionOctTree(OctTree<FunctionNode>::create_default(comm,2)));
        if (!gauss_legendre_test()) comm.Abort();
        if (!test_two_scale_coefficients()) comm.Abort();
    
        for (int i=1; i<argc; i++) {
            if (strcmp(argv[i],"-d") == 0) xterm_debug(comm,0,0);
            if (strcmp(argv[i],"-t") == 0) comm.set_debug(true);
        }
        
        for (int i=0; i<FunctionNode::size; i++) FunctionDataPointersBase::p[i] = 0;
        
        comm.am_register(Function<double>::set_active_handler);
        comm.am_register(Function< std::complex<double> >::set_active_handler);
        comm.am_register(Function<double>::_sock_it_to_me_handler);
        comm.am_register(Function< std::complex<double> >::_sock_it_to_me_handler);
        comm.am_register(Function<double>::recur_down_to_make_handler);
        comm.am_register(Function< std::complex<double> >::recur_down_to_make_handler);
        
        taskq.register_generic_op(Function<double>::recur_down_handler);
        taskq.register_generic_op(Function< std::complex<double> >::recur_down_handler);
    
        return comm;
    }
}
