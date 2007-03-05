#include <mpi.h>
#include <iostream>

using namespace std;

int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    int np = MPI::COMM_WORLD.Get_size();
    if (np != 2) throw "2 only";
    
    int me = MPI::COMM_WORLD.Get_rank();
    int other = me? 0 : 1;

    int a=0, b=-1;
    MPI::Request rsend = MPI::COMM_WORLD.Isend(&a, sizeof(a), MPI::BYTE, other, 1);
    MPI::Request rrecv = MPI::COMM_WORLD.Irecv(&b, sizeof(b), MPI::BYTE, other, 1);

    MPI::Status status;

    while (!rsend.Get_status(status));
    while (!rrecv.Get_status(status));
    rsend.Test(status);
    rrecv.Test(status);
    
    //while (!rsend.Test(status)) ;
    //while (!rrecv.Test(status)) ;

    cout << me << " got " << b << endl;

    MPI::Finalize();
    return 0;
}
