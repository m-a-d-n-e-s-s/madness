#ifndef MAD_MPIWRAP_H
#define MAD_MPIWRAP_H

#if defined(__cplusplus)
extern "C" {
#endif

    int madmpi_any_source();
    int madmpi_any_tag();
    int madmpi_comm_world();
    int madmpi_status_source(const int* data);
    int madmpi_status_count(const int* data);
    int madmpi_status_tag(const int* data);
    int madmpi_status_err(const int* data);
    void madmpi_init(int* argc, char***argv);
    void madmpi_finalize();
    int madmpi_comm_size(int comm);
    int madmpi_comm_rank(int comm);
    void madmpi_comm_bcast(void* buf, int lenbuf, int root, int comm);
    void madmpi_comm_send(const void* buf, int lenbuf, int dest, int tag, int comm);
    void madmpi_comm_recv(void* buf, int lenbuf, int dest, int tag, int comm, void* status);
    void madmpi_iprobe(int src, int tag, int comm, int* flag, void* status);

#if defined(__cplusplus)
}
#endif

#endif
