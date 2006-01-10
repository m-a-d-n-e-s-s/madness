#include <iostream>
using std::cout;
using std::endl;

/// \file octtree/test.cc

#include <fstream>

#include <octtree/octtree.h>
using namespace madness;
    
int main (int argc, char **argv) {
    
    MADMPIInit(argc, argv);
    Communicator comm;
    std::ofstream fout;

    // Redirect cout and cerr for nodes other than process 0
    if (comm.rank()) {
        char filename[256];
        std::sprintf(filename,"log.%5.5ld",comm.rank());
        fout.open(filename);
        cout.rdbuf(fout.rdbuf());
        std::cerr.rdbuf(fout.rdbuf());
    }
    
    comm.print();
    OctTree<double>* layout=0;
    try {
        layout = OctTree<double>::create_default(comm, 2);
        layout->print();
    }
    catch (char const* msg) {
        std::cerr << "Exception (string): " << msg << std::endl;
        comm.abort();
    }
    catch (std::exception& e) {
        std::cerr << "Exception (std): " << e.what() << std::endl;
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

    // doesn't work at the moment
#if 0
    OctTree<double> p;
    cout << " p layout " << (void *) (p.layout());
    OctTree<double>::set_default_layout(layout);
    OctTree<double> q;
    cout << " q layout " << (void *) (q.layout()) << endl;
#endif
 
    for (int i=0; i<10; i++) {
        cout << "walking down" << endl;
        layout->walk_down();
        cout << "walking up" << endl;
        layout->walk_up();
    }
    
    MADMPIFinalize();
 
    return 0;
}
    
