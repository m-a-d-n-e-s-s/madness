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
    }
    else {
        left & sum;
        sum++;
        right & sum;
    }

    // Send a vector around a ring, appending to it
    vector<int> w;
    vector< vector<int> > v(10, w);
//    vector< vector<int> > v(1);
//int npasses = 16776000;
int npasses = 16776000*2;
//int npasses = 16;

  for (int i = 0; i < 10; i++)
  {
//    v.clear();
    if (rank == 0)
        v[i][0] = 0;
    if (rank == 0) {
//        v[0] = 0;
//	v.push_back(rank);
	v[i].insert(v[i].end(), npasses, rank);
        right & v[i];
        left & v[i];
//        cout << "final vector " << v << endl;
	cout << "pass number " << i << endl;
    }
    else {
        left & v[i];
//        v.push_back(me);
	v[i].insert(v[i].end(), npasses, me);
        right & v[i];
    }
  }

    comm.close();
    MPI::Finalize();
    return 0;
}
