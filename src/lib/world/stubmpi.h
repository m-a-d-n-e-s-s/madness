#ifndef MADNESS_STUBMPI_H
#define MADNESS_STUBMPI_H

#include <cstddef>
#include <cstdlib>
#include <cstring>

#define MPI_COMM_WORLD (0x44000000)
#define MPI_UNDEFINED      (-32766)

static void MPI_Abort(int comm, int status) {
    std::exit(1);
};

namespace MPI {

    const int THREAD_SERIALIZED = 12345;
    const int ANY_SOURCE = -2;

    enum Datatype {BYTE};
    struct Op {};
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

    void Finalize();

    int Init_thread(int &argc, char **&argv, int required);

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
            if (datatype != BYTE) throw "die scum!";
            std::memcpy(recvbuf, sendbuf, count);
        }

        void Allreduce(void* sendbuf, void* recvbuf, int count, const MPI::Datatype& datatype, const MPI::Op& op) const {
            if (datatype != BYTE) throw "die scum!";
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
