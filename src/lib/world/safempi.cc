/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  $Id$
*/
#include <world/safempi.h>
//#include <world/worldthread.h>
#include <world/worldexc.h>

namespace SafeMPI {

#ifdef SERIALIZE_MPI
    madness::SCALABLE_MUTEX_TYPE charon;
#endif

    bool Request::Test() {
        SAFE_MPI_GLOBAL_MUTEX;
        return MPI::Request::Test();
    }

    bool Request::Test_got_lock_already() {
        return MPI::Request::Test();
    }

    bool Request::Testany(int n, SafeMPI::Request* request, int& ind) {
        SAFE_MPI_GLOBAL_MUTEX;
        return MPI::Request::Testany(n, static_cast<MPI::Request*>(request), ind);
    }

    int Request::Testsome(int n, SafeMPI::Request* request, int* ind, MPI::Status* status) {
        SAFE_MPI_GLOBAL_MUTEX;
        return MPI::Request::Testsome(n, static_cast<MPI::Request*>(request), ind, status);
    }

    int Request::Testsome(int n, SafeMPI::Request* request, int* ind) {
        SAFE_MPI_GLOBAL_MUTEX;
        return MPI::Request::Testsome(n, static_cast<MPI::Request*>(request), ind);
    }

    int Request::Testsome(int n, MPI::Request* request, int* ind, MPI::Status* status) {
        SAFE_MPI_GLOBAL_MUTEX;
        return MPI::Request::Testsome(n, request, ind, status);
    }


    Intracomm::Intracomm(MPI::Intracomm& comm) : comm(comm) {
        SAFE_MPI_GLOBAL_MUTEX;
        me = comm.Get_rank();
        numproc = comm.Get_size();
    }

    ::SafeMPI::Request Intracomm::Isend(const void* buf, size_t count, const MPI::Datatype& datatype, int dest, int tag) const {
        SAFE_MPI_GLOBAL_MUTEX;
        return comm.Isend(buf,count,datatype,dest,tag);
    }

    ::SafeMPI::Request Intracomm::Irecv(void* buf, size_t count, const MPI::Datatype& datatype, int src, int tag) const {
        SAFE_MPI_GLOBAL_MUTEX;
        return comm.Irecv(buf, count, datatype, src, tag);
    }

    void Intracomm::Send(const void* buf, size_t count, const MPI::Datatype& datatype, int dest, int tag) const {
        SAFE_MPI_GLOBAL_MUTEX;
        comm.Send(buf,count,datatype,dest,tag);
    }

    void Intracomm::Recv(void* buf, int count, const MPI::Datatype& datatype, int source, int tag, MPI::Status& status) const {
        SAFE_MPI_GLOBAL_MUTEX;
        comm.Recv(buf,count,datatype,source,tag,status);
    }

    void Intracomm::Recv(void* buf, int count, const MPI::Datatype& datatype, int source, int tag) const {
        SAFE_MPI_GLOBAL_MUTEX;
        comm.Recv(buf,count,datatype,source,tag);
    }

    void Intracomm::Bcast(void* buf, size_t count, const MPI::Datatype& datatype, int root) const {
        SAFE_MPI_GLOBAL_MUTEX;
        return comm.Bcast(buf, count, datatype, root);
    }

    void Intracomm::Reduce(void* sendbuf, void* recvbuf, int count, const MPI::Datatype& datatype, const MPI::Op& op, int root) const {
        SAFE_MPI_GLOBAL_MUTEX;
        comm.Reduce(sendbuf, recvbuf, count, datatype, op, root);
    }

    void Intracomm::Allreduce(void* sendbuf, void* recvbuf, int count, const MPI::Datatype& datatype, const MPI::Op& op) const {
        SAFE_MPI_GLOBAL_MUTEX;
        comm.Allreduce(sendbuf, recvbuf, count, datatype, op);
    }
    void Intracomm::Get_attr(int key, void* value) const {
        SAFE_MPI_GLOBAL_MUTEX;
        comm.Get_attr(key, value);
    }

    void Intracomm::Abort(int code) const {
        comm.Abort(code);
    }

    void Intracomm::Barrier() const {
        SAFE_MPI_GLOBAL_MUTEX;
        comm.Barrier();
    }

    /// Returns a unique tag for temporary use (1023<tag<=4095)

    /// These tags are intended for one time use to avoid tag
    /// collisions with other messages around the same time period.
    /// It simply increments/wraps a counter and returns the next
    /// legal value.
    ///
    /// So that send and receiver agree on the tag all processes
    /// need to call this routine in the same sequence.
    int Intracomm::unique_tag() {
        SAFE_MPI_GLOBAL_MUTEX;
        static volatile int tag = 1024;
        int result = tag++;
        if (tag >= 4095) tag = 1024;
        return result;
    }

    /// Returns a unique tag reserved for long-term use (0<tag<1000)

    /// Get a tag from this routine for long-term/repeated use.
    ///
    /// Tags in [1000,1023] are statically assigned.
    int Intracomm::unique_reserved_tag() {
        SAFE_MPI_GLOBAL_MUTEX;
        static volatile int tag = 1;
        int result = tag++;
        if (result >= 1000) MADNESS_EXCEPTION( "too many reserved tags in use" , result );
        return result;
    }

    int Intracomm::rank() const { return Get_rank(); }

    int Intracomm::nproc() const { return Get_size(); }

    int Intracomm::size() const { return Get_size(); }

    /// Construct info about a binary tree with given root

    /// Constructs a binary tree spanning the communicator with
    /// process root as the root of the tree.  Returns the logical
    /// parent and children in the tree of the calling process.  If
    /// there is no parent/child the value -1 will be set.
    void Intracomm::binary_tree_info(int root, int& parent, int& child0, int& child1) {
        int np = nproc();
        int me = (rank()+np-root)%np;   // Renumber processes so root has me=0
        parent = (((me-1)>>1)+root)%np;        // Parent in binary tree
        child0 = (me<<1)+1+root;        // Left child
        child1 = (me<<1)+2+root;        // Right child
        if (child0 >= np && child0<(np+root)) child0 -= np;
        if (child1 >= np && child1<(np+root)) child1 -= np;

        if (me == 0) parent = -1;
        if (child0 >= np) child0 = -1;
        if (child1 >= np) child1 = -1;
    }

} //namespace SafeMPI
