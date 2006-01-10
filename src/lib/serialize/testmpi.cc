#include <iostream>
using std::cout;
using std::endl;

/// \file serialize/testmpi.cc
/// \brief Tests serialization over mpi

#define ARCHIVE_REGISTER_TYPE_INSTANTIATE_HERE

#include <serialize/mpiar.h>
using madness::archive::MPIInputArchive;
using madness::archive::MPIOutputArchive;
using madness::archive::wrap;

using namespace std;

#include <misc/misc.h>
using madness::Communicator;
using madness::redirectio;


int main(int argc, char** argv) {
    MADMPIInit(argc, argv);
    Communicator comm;
    redirectio(comm);
    comm.print();

    int nproc = comm.nproc();
    ProcessID rank = comm.rank();

    if (nproc < 2) return 1;

    MPIOutputArchive right(comm,(rank+1)%nproc);
    MPIInputArchive left(comm,(rank+nproc-1)%nproc);

    cout << "ranks " << rank << " " << ((rank+1)%nproc) << " " << ((rank+nproc-1)%nproc) << endl;
    
    int sum = 1;

    if (rank == 0) {
        right & sum;
        left & sum;
        cout << "final sum " << sum << " " << nproc << endl;
    }
    else {
        left & sum;
        sum++;
        right & sum;
    }
    

    MADMPIFinalize();
    return 0;
}
