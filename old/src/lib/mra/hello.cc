#include <iostream>
using std::cout;
using std::endl;

#include "mra.h"
#include "misc.h"
#include "communicator.h"
#include "twoscale.h"
#include "legendre.h"
#include "tensor.h"

using namespace madness;

int main(int argc, char* argv[]) {
    // The following should be used to setup all calculations
    // 1) Initialize parallel environment
    // 2) Setup communication information
    // 3) Redirect standard output+err for parallel processes
    // 4) Load coeffs and quadrature information from file
    PBEGIN_(argc, argv);
    Communicator comm;
    redirectio(comm);
    comm.print();
    load_coeffs(comm);
    load_quadrature(comm);

    // To ensure reliable cleanup catch all C++ exceptions here
    try {
        // Calculation here
        print("hello from", comm.rank());
    } catch (char const* msg) {
        std::cerr << "Exception (string): " << msg << std::endl;
        comm.abort();
    } catch (std::exception& e) {
        std::cerr << "Exception (std): " << e.what() << std::endl;
        comm.abort();
    } catch (TensorException& e) {
        std::cerr << e << std::endl;
        comm.abort();
    }
#ifdef USE_MPI
    catch (MPI::Exception& e) {
        std::cerr << "Exception (mpi): code=" << e.Get_error_code()
        << ", class=" << e.Get_error_class()
        << ", string=" << e.Get_error_string() << std::endl;
        comm.abort();
    }
#endif
    catch (...) {
        std::cerr << "Exception (general)" << std::endl;
        comm.abort();
    }

    // The follwing should be used for succesful termination
    PEND_();
    return 0;
}

