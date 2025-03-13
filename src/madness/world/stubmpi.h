#ifndef MADNESS_STUBMPI_H
#define MADNESS_STUBMPI_H

#include <madness/madness_config.h>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <madness/world/madness_exception.h>
#include <madness/world/timers.h>

typedef int MPI_Group;
typedef int MPI_Request;
typedef struct MPI_Status {
    int count;
    int cancelled;
    int MPI_SOURCE;
    int MPI_TAG;
    int MPI_ERROR;

} MPI_Status;
#define MPI_STATUS_IGNORE ((MPI_Status *)1)
#define MPI_STATUSES_IGNORE ((MPI_Status *)1)

typedef int MPI_Comm;
#define MPI_COMM_WORLD (0x44000000)
#define MPI_UNDEFINED      (-32766)

typedef int MPI_Errhandler;

typedef int MPI_Info;

typedef std::ptrdiff_t MPI_Aint;

/* MPI's error classes */
/* these constants are consistent with MPICH2 mpi.h */
#define MPI_SUCCESS          0      /* Successful return code */
#define MPI_ERR_COMM         5      /* Invalid communicator */
#define MPI_ERR_ARG         12      /* Invalid argument */
#define MPI_ERR_IN_STATUS 999999
#define MPI_ERRORS_RETURN 888888
#define MPI_MAX_ERROR_STRING 1024

/* Results of the compare operations. */
#define MPI_IDENT     0
#define MPI_CONGRUENT 1
#define MPI_SIMILAR   2
#define MPI_UNEQUAL   3

/* MPI null objects */
#define MPI_COMM_NULL       ((MPI_Comm)0x04000000)
#define MPI_OP_NULL         ((MPI_Op)0x18000000)
#define MPI_GROUP_NULL      ((MPI_Group)0x08000000)
#define MPI_DATATYPE_NULL   ((MPI_Datatype)0x0c000000)
#define MPI_REQUEST_NULL    ((MPI_Request)0x2c000000)
#define MPI_ERRHANDLER_NULL ((MPI_Errhandler)0x14000000)

/* MPI thread support levels */
/* these constants are consistent with MPICH2 mpi.h */
#define MPI_THREAD_SINGLE     0
#define MPI_THREAD_FUNNELED   1
#define MPI_THREAD_SERIALIZED 2
#define MPI_THREAD_MULTIPLE   3
#define MPI_COMM_TYPE_SHARED  4

/* these constants are consistent with MPICH2 mpi.h */
#define MPI_IN_PLACE   ((void *) -1)
#define MPI_PROC_NULL  -1
#define MPI_ANY_SOURCE -2
#define MPI_ANY_TAG    -1

/* MPI data types */
/* these constants are consistent with MPICH2 mpi.h */
typedef int MPI_Datatype;
#define MPI_CHAR               ((MPI_Datatype)0x4c000101)
#define MPI_SIGNED_CHAR        ((MPI_Datatype)0x4c000118)
#define MPI_UNSIGNED_CHAR      ((MPI_Datatype)0x4c000102)
#define MPI_BYTE               ((MPI_Datatype)0x4c00010d)
#define MPI_WCHAR              ((MPI_Datatype)0x4c00040e)
#define MPI_SHORT              ((MPI_Datatype)0x4c000203)
#define MPI_UNSIGNED_SHORT     ((MPI_Datatype)0x4c000204)
#define MPI_INT                ((MPI_Datatype)0x4c000405)
#define MPI_UNSIGNED           ((MPI_Datatype)0x4c000406)
#define MPI_LONG               ((MPI_Datatype)0x4c000807)
#define MPI_UNSIGNED_LONG      ((MPI_Datatype)0x4c000808)
#define MPI_FLOAT              ((MPI_Datatype)0x4c00040a)
#define MPI_DOUBLE             ((MPI_Datatype)0x4c00080b)
#define MPI_LONG_DOUBLE        ((MPI_Datatype)0x4c00100c)
#define MPI_LONG_LONG_INT      ((MPI_Datatype)0x4c000809)
#define MPI_UNSIGNED_LONG_LONG ((MPI_Datatype)0x4c000819)
#define MPI_LONG_LONG          ((MPI_Datatype)0x4c000809)

inline int MPI_Type_get_extent(MPI_Datatype datatype, MPI_Aint *lb,
                               MPI_Aint *extent) {
  switch(datatype) {
  case MPI_CHAR:
    *extent = sizeof(char); break;
  case MPI_SIGNED_CHAR:
    *extent = sizeof(signed char); break;
  case MPI_UNSIGNED_CHAR:
    *extent = sizeof(unsigned char); break;
  case MPI_BYTE:
    *extent = 1; break;
  case MPI_WCHAR:
    *extent = sizeof(wchar_t); break;
  case MPI_SHORT:
    *extent = sizeof(short); break;
  case MPI_UNSIGNED_SHORT:
    *extent = sizeof(unsigned short); break;
  case MPI_INT:
    *extent = sizeof(int); break;
  case MPI_UNSIGNED:
    *extent = sizeof(unsigned); break;
  case MPI_LONG:
    *extent = sizeof(long); break;
  case MPI_UNSIGNED_LONG:
    *extent = sizeof(unsigned long); break;
  case MPI_FLOAT:
    *extent = sizeof(float); break;
  case MPI_DOUBLE:
    *extent = sizeof(double); break;
  case MPI_LONG_DOUBLE:
    *extent = sizeof(long double); break;
  case MPI_LONG_LONG_INT: // same as MPI_LONG_LONG
    *extent = sizeof(long long int); break;
  case MPI_UNSIGNED_LONG_LONG:
    *extent = sizeof(unsigned long long); break;
  default:
    *extent = MPI_UNDEFINED;
  }
  *lb = 0;
  return MPI_SUCCESS;
}

/* MPI Reduction operation */
/* these constants are consistent with MPICH2 mpi.h */
typedef int MPI_Op;
#define MPI_MAX     ((MPI_Op)0x58000001)
#define MPI_MIN     ((MPI_Op)0x58000002)
#define MPI_SUM     ((MPI_Op)0x58000003)
#define MPI_PROD    ((MPI_Op)0x58000004)
#define MPI_LAND    ((MPI_Op)0x58000005)
#define MPI_BAND    ((MPI_Op)0x58000006)
#define MPI_LOR     ((MPI_Op)0x58000007)
#define MPI_BOR     ((MPI_Op)0x58000008)
#define MPI_LXOR    ((MPI_Op)0x58000009)
#define MPI_BXOR    ((MPI_Op)0x5800000a)
#define MPI_MINLOC  ((MPI_Op)0x5800000b)
#define MPI_MAXLOC  ((MPI_Op)0x5800000c)
#define MPI_REPLACE ((MPI_Op)0x5800000d)

/* function type given to MPI_Op_create */
typedef void (MPI_User_function) ( void * a,
               void * b, int * len, MPI_Datatype * );

inline int MPI_Group_translate_ranks(MPI_Group, int, const int [],
                            MPI_Group, int ranks2[]) {
    ranks2[0] = 0;
    return MPI_SUCCESS;
}

/* TODO The NO-OP implementation of may not be adequate. */
inline int MPI_Group_incl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup) {
    return MPI_SUCCESS;
}

/* TODO The NO-OP implementation may not be adequate. */
inline int MPI_Group_free(MPI_Group *group) {
    return MPI_SUCCESS;
}

// Initialization and finalize functions
inline int MPI_Init(int *, char ***) { return MPI_SUCCESS; }
inline int MPI_Init_thread(int *, char ***, int, int *provided) { *provided = MADNESS_MPI_THREAD_LEVEL; return MPI_SUCCESS; }
inline int MPI_Initialized(int* flag) { *flag = 1; return MPI_SUCCESS; }
inline int MPI_Finalize() { return MPI_SUCCESS; }
inline int MPI_Finalized(int* flag) { *flag = 0; return MPI_SUCCESS; }
inline int MPI_Query_thread(int *provided) { *provided = MADNESS_MPI_THREAD_LEVEL; return MPI_SUCCESS; }

// Buffer functions (do nothing since no messages may be sent)
inline int MPI_Buffer_attach(void*, int) { return MPI_SUCCESS; }
inline int MPI_Buffer_detach(void* buffer, int* size) { return MPI_SUCCESS; }

inline int MPI_Test(MPI_Request *, int *flag, MPI_Status *) {
    *flag = 1;
    return MPI_SUCCESS;
}

inline int MPI_Testany(int, MPI_Request[], int* index, int *flag, MPI_Status*) {
    *index = MPI_UNDEFINED;
    *flag = 1;
    return MPI_SUCCESS;
}

inline int MPI_Testsome(int, MPI_Request*, int *outcount, int*, MPI_Status*) {
    *outcount = MPI_UNDEFINED;
    return MPI_SUCCESS;
}

inline int MPI_Get_count(MPI_Status *, MPI_Datatype, int *count) {
    *count = 0;
    return MPI_SUCCESS;
}

// Communicator rank and size
inline int MPI_Comm_rank(MPI_Comm, int* rank) { *rank = 0; return MPI_SUCCESS; }
inline unsigned int MPI_Comm_size(MPI_Comm, int* size) { *size = 1; return MPI_SUCCESS; }

// There is only one node so sending messages is not allowed. Always return MPI_ERR_COMM
inline int MPI_Isend(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *) { return MPI_ERR_COMM; }
inline int MPI_Issend(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *) { return MPI_ERR_COMM; }
inline int MPI_Send(void*, int, MPI_Datatype, int, int, MPI_Comm) { return MPI_ERR_COMM; }
inline int MPI_Ssend(void*, int, MPI_Datatype, int, int, MPI_Comm) { return MPI_ERR_COMM; }
inline int MPI_Bsend(void*, int, MPI_Datatype, int, int, MPI_Comm) { return MPI_ERR_COMM; }
inline int MPI_Irecv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request*) { return MPI_ERR_COMM; }
inline int MPI_Recv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status*) { return MPI_ERR_COMM; }

// Gather = copy
inline int MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                void *recvbuf, const int recvcounts[], const int displs[], MPI_Datatype recvtype,
                int root, MPI_Comm) {
  MPI_Aint recvtype_extent;
  MPI_Aint recvtype_lb;
  MPI_Type_get_extent(recvtype, &recvtype_lb, &recvtype_extent);
  MADNESS_ASSERT(recvtype_lb == 0);
  MPI_Aint sendtype_extent;
  MPI_Aint sendtype_lb;
  MPI_Type_get_extent(sendtype, &sendtype_lb, &sendtype_extent);
  MADNESS_ASSERT(sendtype_lb == 0);
  MADNESS_ASSERT(sendcount * sendtype_extent <= recvcounts[0] * recvtype_extent);
  std::memcpy(recvbuf, sendbuf, sendcount * sendtype_extent);
  return MPI_SUCCESS;
}
inline int MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
               void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
  const int recvcounts[1] = {recvcount};
  const int displs[1] = {0};
  return MPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm);
}

// Bcast does nothing but return MPI_SUCCESS
inline int MPI_Bcast(void*, int, MPI_Datatype, int, MPI_Comm) { return MPI_SUCCESS; }

// Reduce does memcpy and returns MPI_SUCCESS
inline int MPI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype, MPI_Op, int, MPI_Comm) {
    if(sendbuf != MPI_IN_PLACE) std::memcpy(recvbuf, sendbuf, count);
    return MPI_SUCCESS;
}
inline int MPI_Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype, MPI_Op, MPI_Comm) {
    if(sendbuf != MPI_IN_PLACE) std::memcpy(recvbuf, sendbuf, count);
    return MPI_SUCCESS;
}

inline int MPI_Comm_get_attr(MPI_Comm, int, void*, int*) { return MPI_ERR_COMM; }

inline int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm) {
  //*newcomm = MPI_COMM_NULL;
  *newcomm = MPI_COMM_WORLD; // So that can successfully split a 1 process communicator
  return MPI_SUCCESS;
}

inline int MPI_Comm_split_type(MPI_Comm comm, int split_type, int key, MPI_Info info, MPI_Comm *newcomm) {
  *newcomm = (split_type == MPI_UNDEFINED) ? MPI_COMM_NULL : comm;
  return MPI_SUCCESS;
}

inline int MPI_Abort(MPI_Comm, int code) { exit(code); return MPI_SUCCESS; }

inline int MPI_Barrier(MPI_Comm) { return MPI_SUCCESS; }

inline int MPI_Comm_create(MPI_Comm,  MPI_Group, MPI_Comm *newcomm) {
    //*newcomm = MPI_COMM_NULL;
    *newcomm = MPI_COMM_WORLD; // So that can successfully split a 1 process communicator
    return MPI_SUCCESS;
}

inline int MPI_Comm_group(MPI_Comm, MPI_Group* group) {
    *group = MPI_GROUP_NULL;
    return MPI_SUCCESS;
}

inline int MPI_Comm_free(MPI_Comm * comm) {
    *comm = MPI_COMM_NULL;
    return MPI_SUCCESS;
}

inline int MPI_Comm_compare(MPI_Comm comm1, MPI_Comm comm2, int *result) {
    if (comm1 == comm2) {
        *result = MPI_IDENT;
    } else {
        *result = MPI_UNEQUAL;
    }
    return MPI_SUCCESS;
}


inline int MPI_Error_string(int errorcode, char *string, int *resultlen) {
    switch(errorcode) {
        case MPI_SUCCESS:
            *resultlen = 8;
            std::strncpy(string, "Success", *resultlen);
            break;
        case MPI_ERR_COMM:
            *resultlen = 21;
            std::strncpy(string, "Invalid communicator", *resultlen);
            break;
        case MPI_ERR_ARG:
            *resultlen = 17;
            std::strncpy(string, "Invalid argument", *resultlen);
            break;
        default:
            return MPI_ERR_ARG;
            break;
    }

    return MPI_SUCCESS;
}

inline int MPI_Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler errhandler) {return MPI_SUCCESS;}

inline double MPI_Wtime() { return madness::wall_time(); }

inline int MPI_Op_create(MPI_User_function *user_fn, int commute, MPI_Op *op) {
  return MPI_SUCCESS;
}

inline int MPI_Op_free(MPI_Op *op) {
  return MPI_SUCCESS;
}

inline int MPI_Info_create (MPI_Info *info) { return MPI_SUCCESS; }
inline int MPI_Info_free (MPI_Info *info) { return MPI_SUCCESS; }

#endif
