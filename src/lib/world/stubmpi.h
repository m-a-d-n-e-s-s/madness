#ifndef MADNESS_STUBMPI_H
#define MADNESS_STUBMPI_H

#include <cstddef>
#include <cstdlib>
#include <cstring>

#define MPI_COMM_WORLD (0x44000000)
#define MPI_UNDEFINED      (-32766)

namespace MPI {

    /* these constants are consistent with MPICH2 mpi.h */
    const int THREAD_SINGLE     = 0;
    const int THREAD_FUNNELED   = 1;
    const int THREAD_SERIALIZED = 2;
    const int THREAD_MULTIPLE   = 3;

    /* these constants are consistent with MPICH2 mpi.h */
    const int IN_PLACE   = -1;
    const int PROC_NULL  = -1;
    const int ANY_SOURCE = -2;
    const int ANY_TAG    = -1;

    /* these constants are consistent with MPICH2 mpi.h */
    enum Datatype {
        CHAR               = 0x4c000101,
        SIGNED_CHAR        = 0x4c000118,
        UNSIGNED_CHAR      = 0x4c000102,
        BYTE               = 0x4c00010d,
        WCHAR              = 0x4c00040e,
        SHORT              = 0x4c000203,
        UNSIGNED_SHORT     = 0x4c000204,
        INT                = 0x4c000405,
        UNSIGNED           = 0x4c000406,
        LONG               = 0x4c000807,
        UNSIGNED_LONG      = 0x4c000808,
        FLOAT              = 0x4c00040a,
        DOUBLE             = 0x4c00080b,
        LONG_DOUBLE        = 0x4c00100c,
        LONG_LONG_INT      = 0x4c000809,
        UNSIGNED_LONG_LONG = 0x4c000819,
        LONG_LONG          = 0x4c000809,
    };

    /* these constants are consistent with MPICH2 mpi.h */
    enum Op {
        MAX     = 0x58000001,
        MIN     = 0x58000002,
        SUM     = 0x58000003,
        PROD    = 0x58000004,
        LAND    = 0x58000005,
        BAND    = 0x58000006,
        LOR     = 0x58000007,
        BOR     = 0x58000008,
        LXOR    = 0x58000009,
        BXOR    = 0x5800000a,
        MINLOC  = 0x5800000b,
        MAXLOC  = 0x5800000c,
        REPLACE = 0x5800000d,
    };

    struct Group {
        static void Translate_ranks(const Group& v1, int v2, const int* v3, const Group& v4, int* v5) {
          *v5=0;
        }
    };

    class Exception {
        
    };

    struct Status {
        int Get_source() {
            throw "not implemented";
            return -1;
        }
        int Get_count(const Datatype &v2) {
            throw "not implemented";
            return -1;
        }
    };

    void Finalize() {
        return;
    }

    bool Is_finalized() {
         return true;
    }

    int Init_thread(int &argc, char **&argv, int required) {
            return THREAD_SERIALIZED; /* none of the functions defined in this file have side-effects */
    }

#ifdef MADNESS_USE_BSEND_ACKS
    void Attach_buffer(void* buffer, int size) {
        return;
    }

    int Detach_buffer(void*& buffer) {
        return 0;
    }
#endif // MADNESS_USE_BSEND_ACKS

    struct Request {
        bool Test() {
            return false;
        }

        static bool Testany(int n, Request* request, int& ind) {
            return false;
        }

        static int Testsome(int n, Request* request, int* ind, MPI::Status* status) {
            return false;
        }

        static int Testsome(int n, Request* request, int* ind) {
            return false;
        }
    };

    struct Intracomm {
        int Get_rank() const {
            return 0;
        }

        int Get_size() const {
            return 1;
        }

        Request Isend(const void* buf, size_t count, const MPI::Datatype& datatype, int dest, int tag) const {
            throw "not implemented";
        }

        Request Irecv(void* buf, size_t count, const MPI::Datatype& datatype, int src, int tag) const {
            throw "not implemented";
        }

        void Send(const void* buf, size_t count, const MPI::Datatype& datatype, int dest, int tag) const {
            throw "not implemented";
        }

#ifdef MADNESS_USE_BSEND_ACKS
        void Bsend(const void* buf, size_t count, const MPI::Datatype& datatype, int dest, int tag) const {
            throw "not implemented";
        }
#endif // MADNESS_USE_BSEND_ACKS

        void Recv(void* buf, int count, const MPI::Datatype& datatype, int source, int tag, MPI::Status& status) const {
            throw "not implemented";
        }

        void Recv(void* buf, int count, const MPI::Datatype& datatype, int source, int tag) const {
            throw "not implemented";
        }

        void Bcast(void* buf, size_t count, const MPI::Datatype& datatype, int root) const {
            return;
        }

        void Reduce(void* sendbuf, void* recvbuf, int count, const MPI::Datatype& datatype, const MPI::Op& op, int root) const {
            /* TODO: support other build-in datatypes */
            if (datatype != BYTE) throw "die scum!";
            /* TODO use memmove or check for IN_PLACE */
            std::memcpy(recvbuf, sendbuf, count);
        }

        void Allreduce(void* sendbuf, void* recvbuf, int count, const MPI::Datatype& datatype, const MPI::Op& op) const {
            /* TODO: support other build-in datatypes */
            if (datatype != BYTE) throw "die scum!";
            /* TODO use memmove or check for IN_PLACE */
            std::memcpy(recvbuf, sendbuf, count);
        }

        void Get_attr(int key, void* value) const {
            throw "not implemented";
        }

        void Abort(int code=1) const {
            exit(code);
        }

        void Barrier() const {
            return;
        }

        Intracomm Create(const MPI::Group& group) const {
            return Intracomm();
        }

        Group Get_group() const {
            return Group();
        }


    };

   static Intracomm COMM_WORLD;
}

#endif
