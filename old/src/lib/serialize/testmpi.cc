#include <iostream>
using std::cout;
using std::endl;

#include <vector>
using std::vector;

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

template <typename T>
ostream& operator<<(ostream& s, const std::vector<T>& v) {
    for (unsigned int i=0; i<v.size(); i++) s << v[i] << " ";
    return s;
}

int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    Communicator comm;
    redirectio(comm);
    comm.print();

    int me = comm.rank();
    int nproc = comm.nproc();
    ProcessID rank = comm.rank();

    //comm.set_debug(true);

    if (nproc < 2) return 1;

    MPIOutputArchive right(comm,(rank+1)%nproc);
    MPIInputArchive left(comm,(rank+nproc-1)%nproc);

    cout << "ranks " << rank << " " << ((rank+1)%nproc) << " " << ((rank+nproc-1)%nproc) << endl;

    // Send an integer around a ring, accumulating onto it
    int sum = 1;
    if (rank == 0) {
        MPIOutputArchive(comm,(rank+1)%nproc) & sum;
        left & sum;
        cout << "final sum " << sum << " " << nproc << endl;
    } else {
        left & sum;
        sum++;
        right & sum;
        right.flush();
    }

    // Send a vector around a ring, appending to it
    vector<int> v(1);
    if (rank == 0) {
        v[0] = 0;
        right & v;
        right.flush();
        left & v;
        cout << "final vector " << v << endl;
    } else {
        left & v;
        v.push_back(me);
        right & v;
        right.flush();
    }

    comm.close();
    MPI::Finalize();
    return 0;
}
