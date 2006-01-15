#include <mpi.h>
#include "madmpiwrap.h"

/*** 
     \file madmpiwrap.c
     \brief Wrappers for mingw that does not have C++ MPI interface
***/

int madmpi_any_source() {
  return MPI_ANY_SOURCE;
}

int madmpi_any_tag() {
  return MPI_ANY_TAG;
}

int madmpi_comm_world() {
    return MPI_COMM_WORLD;
}
int madmpi_status_count(const int* data) {
    return ((MPI_Status*) data)->count; /*MPI_COUNT;*/
}
int madmpi_status_source(const int* data) {
    return ((MPI_Status*) data)->MPI_SOURCE;
}
int madmpi_status_tag(const int* data) {
    return ((MPI_Status*) data)->MPI_TAG;
}
int madmpi_status_err(const int* data) {
    return ((MPI_Status*) data)->MPI_ERROR;
}
void madmpi_init(int* argc, char*** argv) {
    MPI_Init(argc, argv);
}
void madmpi_finalize() {
    MPI_Finalize();
}
int madmpi_comm_size(int comm) {
    int n;
    MPI_Comm_size(comm, &n);
    return n;
}
int madmpi_comm_rank(int comm) {
    int n;
    MPI_Comm_rank(comm, &n);
    return n;
}
void madmpi_comm_bcast(void* buf, int lenbuf, int root, int comm) {
    printf("buf=%p lenbuf=%d root=%d comm=%d\n", buf, lenbuf, root, comm);
    MPI_Bcast(buf, lenbuf, MPI_BYTE, root, comm);
};
void madmpi_comm_send(const void* buf, int lenbuf, int dest, int tag, int comm) {
    MPI_Send((void *) buf, lenbuf, MPI_BYTE, dest, tag, comm);
}
void madmpi_comm_recv(void* buf, int lenbuf, int dest, int tag, int comm, void* status) {
    MPI_Recv(buf, lenbuf, MPI_BYTE, dest, tag, comm, (MPI_Status*) status);
}

void madmpi_iprobe(int src, int tag, int comm, int* flag, void* status) {
    MPI_Probe(src, tag, comm, (MPI_Status*) status);
}

