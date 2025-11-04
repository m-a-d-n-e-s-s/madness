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
#ifndef MADNESS_WORLD_SAFEMPI_H__INCLUDED
#define MADNESS_WORLD_SAFEMPI_H__INCLUDED

/// \file safempi.h
/// \brief Serializes calls to MPI in case it does not support THREAD_MULTIPLE

#include <madness/madness_config.h>

#ifdef STUBOUTMPI
#include <madness/world/stubmpi.h>
#else

//#ifdef SEEK_SET
//#undef SEEK_SET
//#endif
//#ifdef SEEK_CUR
//#undef SEEK_CUR
//#endif
//#ifdef SEEK_END
//#undef SEEK_END
//#endif

#ifdef MADNESS_MPI_HEADER
# include MADNESS_MPI_HEADER
#else
# include <mpi.h>
#endif

#endif


#if MADNESS_MPI_THREAD_LEVEL == MPI_THREAD_SERIALIZED
#  define MADNESS_SERIALIZES_MPI
#endif



#include <madness/world/worldmutex.h>
#include <madness/world/type_traits.h>
#include <iostream>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <sstream>

#define MADNESS_MPI_TEST(condition) \
    {                                       \
        int mpi_error_code = condition;     \
        if(mpi_error_code != MPI_SUCCESS) { \
            char s[MPI_MAX_ERROR_STRING];                             \
            s[0] = '\0';                                                \
            int len = 0;                                                \
            MPI_Error_string(mpi_error_code, s, &len);       \
            std::cout<< "MPI ERROR in " << __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << " code " << mpi_error_code << " err string " << s << "\n"; \
            throw ::SafeMPI::Exception(mpi_error_code);                 \
        }                                                               \
    }

namespace SafeMPI {

    extern madness::SCALABLE_MUTEX_TYPE charon;      // Inside safempi.cc
#ifdef MADNESS_SERIALIZES_MPI
#define SAFE_MPI_GLOBAL_MUTEX madness::ScopedMutex<madness::SCALABLE_MUTEX_TYPE> obolus(SafeMPI::charon);
#else
#define SAFE_MPI_GLOBAL_MUTEX
#endif

    /// tags in [1,999] ... allocated once by unique_reserved_tag
    ///
    /// tags in [1000,1023] ... statically assigned here
    ///
    /// tags in [1024,4095] ... allocated round-robin by unique_tag
    ///
    /// tags in [4096,8191] ... reserved for huge msg exchange by RMI
    ///
    /// tags in [8192,MPI::TAG_UB] ... not used/managed by madness

    static const int RMI_TAG = 1023;
    static const int MPIAR_TAG = 1001;
    static const int DEFAULT_SEND_RECV_TAG = 1000;

    // Forward declarations
    class Intracomm;
    extern Intracomm COMM_WORLD;
    inline int Finalize();

    /// Check MPI initialization status

    /// \return \c true if MPI has been initialized, \c false otherwise.
    inline bool Is_initialized() {
        int initialized = 0;
        MPI_Initialized(&initialized);
        return (initialized != 0);
    }

    /// Check MPI finalization status

    /// \return \c true if MPI has been finalized, \c false otherwise.
    inline bool Is_finalized() {
        int flag = 0;
        MPI_Finalized(&flag);
        return flag != 0;
    }

    namespace detail {
        /// Initialize SafeMPI::COMM_WORLD
        inline void init_comm_world();

        inline void print_mpi_error(const int rc, const char* function,
                const int line, const char* file)
        {
            int len = 0;
            char error_string[MPI_MAX_ERROR_STRING];
            MPI_Error_string(rc, error_string, &len);
            std::cerr << "!!! MPI ERROR (" << rc << ") in " << function <<
                    " at " << file << "(" << line << "): " << error_string << "\n";
        }

    }  // namespace detail

    /// SafeMPI exception object

    /// This exception is thrown whenever an MPI error occurs.
    class Exception : public std::exception {
    private:
        char mpi_error_string_[MPI_MAX_ERROR_STRING];
        std::string mpi_statuses_error_string_;
    public:

        Exception(const int mpi_error) throw() {
            int len = 0;
            if(MPI_Error_string(mpi_error, mpi_error_string_, &len) != MPI_SUCCESS)
                std::strncpy(mpi_error_string_, "UNKNOWN MPI ERROR!", MPI_MAX_ERROR_STRING);
        }

        Exception(const int mpi_error, const int nstatuses,
                  const int* indices,
                  MPI_Status* const statuses) noexcept {
            try {
                if (mpi_error == MPI_ERR_IN_STATUS) {
                    std::ostringstream oss;
                    for(auto s=0; s!=nstatuses; ++s) {
                        int len = 0;
                        auto status_error = statuses[s].MPI_ERROR;
                        if (status_error != MPI_SUCCESS) {
                            oss << "request " << indices[s] << ":";
                            if (MPI_Error_string(status_error, mpi_error_string_, &len) != MPI_SUCCESS)
                                oss << " unknown error!" << std::endl;
                            else
                                oss << mpi_error_string_ << std::endl;
                        }
                    }
                    mpi_statuses_error_string_ = oss.str();
                }
            }
            catch (...) {}

            int len = 0;
            if(MPI_Error_string(mpi_error, mpi_error_string_, &len) != MPI_SUCCESS)
                std::strncpy(mpi_error_string_, "UNKNOWN MPI ERROR!", MPI_MAX_ERROR_STRING);
        }

        Exception(const Exception& other) throw() {
            std::strncpy(mpi_error_string_, other.mpi_error_string_, MPI_MAX_ERROR_STRING);
            try {
                mpi_statuses_error_string_ = other.mpi_statuses_error_string_;
            } catch(...) { mpi_statuses_error_string_.clear(); }
        }

        Exception& operator=(const Exception& other) {
            std::strncpy(mpi_error_string_, other.mpi_error_string_, MPI_MAX_ERROR_STRING);
            try {
                mpi_statuses_error_string_ = other.mpi_statuses_error_string_;
            } catch(...) { mpi_statuses_error_string_.clear(); }
            return *this;
        }

        virtual ~Exception() throw() { }

        virtual const char* what() const throw() { return mpi_error_string_; }
        bool can_elaborate() const noexcept {
            return !mpi_statuses_error_string_.empty();
        }
        const char* elaborate() const noexcept {
            return mpi_statuses_error_string_.c_str();
        }

        friend std::ostream& operator<<(std::ostream& os, const Exception& e) {
            os << e.what();
            if (e.can_elaborate()) {
              os << e.elaborate();
            }
            return os;
        }
    }; // class Exception


    class Status {
    private:
        MPI_Status status_;

    public:
        // Constructors
        Status(void) : status_() { }
        Status(const Status &other) : status_(other.status_) { }
        Status(MPI_Status other) : status_(other) { }

        // Assignment operators
        Status& operator=(const Status &other) {
            status_ = other.status_;
            return *this;
        }

        Status& operator=(const MPI_Status other) {
            status_ = other;
            return *this;
        }

        // C/C++ cast and assignment
        operator MPI_Status*() { return &status_; }

        operator MPI_Status() const { return status_; }

//        bool Is_cancelled() const {
//            int flag = 0;
//            MADNESS_MPI_TEST(MPI_Test_cancelled(const_cast<MPI_Status*>(&status_), &flag));
//            return flag != 0;
//        }
//
//        int Get_elements(const MPI_Datatype datatype) const {
//            int elements = 0;
//            MADNESS_MPI_TEST(MPI_Get_elements(const_cast<MPI_Status*>(&status_), datatype, &elements));
//            return elements;
//        }

        int Get_count(const MPI_Datatype datatype) const {
            int count = 0;
            MADNESS_MPI_TEST(MPI_Get_count(const_cast<MPI_Status*>(&status_), datatype, &count));
            return count;
        }

//        void Set_cancelled(bool flag) {
//            MADNESS_MPI_TEST(MPI_Status_set_cancelled(&status_, flag));
//        }
//
//        void Set_elements( const MPI_Datatype &v2, int v3 ) {
//            MADNESS_MPI_TEST(MPI_Status_set_elements(&status_, v2, v3 ));
//        }

        int Get_source() const { return status_.MPI_SOURCE; }

        int Get_tag() const { return status_.MPI_TAG; }

        int Get_error() const { return status_.MPI_ERROR; }

        void Set_source(int source) { status_.MPI_SOURCE = source; }

        void Set_tag(int tag) { status_.MPI_TAG = tag; }

        void Set_error(int error) { status_.MPI_ERROR = error; }
    }; // class Status

    class Request {
        // Note: This class was previously derived from MPI::Request, but this
        // was changed with the removal of the MPI C++ bindings. Now this class
        // only implements the minimum functionality required by MADNESS. Feel
        // free to add more functionality as needed.

    private:
        MPI_Request request_;

    public:

        // Constructors
        Request() : request_(MPI_REQUEST_NULL) { }
        Request(MPI_Request other) : request_(other) { }
        Request(const Request& other) : request_(other.request_) { }

        // Assignment operators
        Request& operator=(const Request &other) {
            request_ = other.request_;
            return *this;
        }

        Request& operator=(const MPI_Request& other) {
            request_ = other;
            return *this;
        }

        // logical
        bool operator==(const Request &other) { return (request_ == other.request_); }
        bool operator!=(const Request &other) { return (request_ != other.request_); }

        // C/C++ cast and assignment
        operator MPI_Request*() { return &request_; }
        operator MPI_Request() const { return request_; }

        static bool Testany(int count, Request* requests, int& index, Status& status) {
            MADNESS_ASSERT(requests != nullptr);
            int flag;
            std::unique_ptr<MPI_Request[]> mpi_requests(new MPI_Request[count]);

            // Copy requests to an array that can be used by MPI
            for(int i = 0; i < count; ++i)
                mpi_requests[i] = requests[i].request_;
            {
                SAFE_MPI_GLOBAL_MUTEX;
                MADNESS_MPI_TEST(MPI_Testany(count, mpi_requests.get(), &index, &flag, status));
            }
            // Copy results from MPI back to the original array
            for(int i = 0; i < count; ++i)
                requests[i].request_ = mpi_requests[i];
            return flag != 0;
        }

        static bool Testany(int count, Request* requests, int& index) {
            MADNESS_ASSERT(requests != nullptr);
            int flag;
            std::unique_ptr<MPI_Request[]> mpi_requests(new MPI_Request[count]);

            // Copy requests to an array that can be used by MPI
            for(int i = 0; i < count; ++i)
                mpi_requests[i] = requests[i].request_;
            {
                SAFE_MPI_GLOBAL_MUTEX;
                MADNESS_MPI_TEST(MPI_Testany(count, mpi_requests.get(), &index, &flag, MPI_STATUS_IGNORE));
            }
            // Copy results from MPI back to the original array
            for(int i = 0; i < count; ++i)
                requests[i] = mpi_requests[i];
            return flag != 0;
        }

        static int Testsome(int incount, Request* requests, int* indices, Status* statuses) {
            MADNESS_ASSERT(requests != nullptr);
            MADNESS_ASSERT(indices != nullptr);
            MADNESS_ASSERT(statuses != nullptr);

            int outcount = 0;
            std::unique_ptr<MPI_Request[]> mpi_requests(new MPI_Request[incount]);
            std::unique_ptr<MPI_Status[]> mpi_statuses(new MPI_Status[incount]);
            for(int i = 0; i < incount; ++i)
                mpi_requests[i] = requests[i].request_;
            {
                SAFE_MPI_GLOBAL_MUTEX;
                {  // print out the status vars for the failed requests
                  auto mpi_error_code =
                      MPI_Testsome(incount, mpi_requests.get(), &outcount,
                                   indices, mpi_statuses.get());
                  if (mpi_error_code != MPI_SUCCESS) {
                    throw ::SafeMPI::Exception(mpi_error_code, outcount, indices, mpi_statuses.get());
                  }
                }
            }
            for(int i = 0; i < incount; ++i) {
                requests[i] = mpi_requests[i];
                statuses[i] = mpi_statuses[i];
            }
            return outcount;
        }

        static int Testsome(int incount, Request* requests, int* indices) {
            int outcount = 0;
            std::unique_ptr<MPI_Request[]> mpi_requests(new MPI_Request[incount]);
            for(int i = 0; i < incount; ++i)
                mpi_requests[i] = requests[i].request_;
            {
                SAFE_MPI_GLOBAL_MUTEX;
                MADNESS_MPI_TEST( MPI_Testsome( incount, mpi_requests.get(), &outcount, indices, MPI_STATUSES_IGNORE));
            }
            for(int i = 0; i < incount; ++i)
                requests[i] = mpi_requests[i];
            return outcount;
        }


        bool Test_got_lock_already(MPI_Status& status) {
            int flag;
            MADNESS_MPI_TEST(MPI_Test(&request_, &flag, &status));
            return flag != 0;
        }

        bool Test(MPI_Status& status) {
            SAFE_MPI_GLOBAL_MUTEX;
            return Test_got_lock_already(status);
        }

        bool Test_got_lock_already() {
            int flag;
            MADNESS_MPI_TEST(MPI_Test(&request_, &flag, MPI_STATUS_IGNORE));
            return flag != 0;
        }

        bool Test() {
            SAFE_MPI_GLOBAL_MUTEX;
            return Test_got_lock_already();
        }
    }; // class Request

    ///  Wrapper around MPI_Group. Has a shallow copy constructor. Usually deep copy is not needed, but can be created
    ///  via Group::Incl().
    class Group {
    public:
        Group Incl(int n, const int* ranks) const {
            // MPI <3 interface lacks explicit const sanitation
            Group result(std::shared_ptr<Impl>(new Impl(*pimpl, n, const_cast<int*>(ranks))));
            return result;
        }

        void Translate_ranks(int nproc, const int* ranks1, const Group& grp2, int* ranks2) const {
            // MPI <3 interface lacks explicit const sanitation
            MADNESS_ASSERT(pimpl);
            MADNESS_MPI_TEST(MPI_Group_translate_ranks(pimpl->group, nproc,
                    const_cast<int*>(ranks1), grp2.pimpl->group, ranks2));
        }

        MPI_Group group() const {
            MADNESS_ASSERT(pimpl);
            return pimpl->group;
        }

        Group(const Group& other) : pimpl(other.pimpl) { }

    private:

        struct Impl {
            MPI_Group group;

            Impl(MPI_Comm comm) {
                MADNESS_MPI_TEST(MPI_Comm_group(comm, &group));
            }

            Impl(const Impl& other, int n, const int* ranks) {
                // MPI <3 interface lacks explicit const sanitation
                MADNESS_MPI_TEST(MPI_Group_incl(other.group, n,
                        const_cast<int*>(ranks), & group));
            }

            ~Impl() {
                if(Is_initialized()) {
                    const int mpi_error_code = MPI_Group_free(&group);
                    if(mpi_error_code != MPI_SUCCESS)
                        ::SafeMPI::detail::print_mpi_error(mpi_error_code,
                                "SafeMPI::Group::Impl::~Impl()", __LINE__, __FILE__);
                }
            }

        }; // struct Impl

        friend class Intracomm;

        Group() : pimpl() { }

        // only Intracomm will use this
        Group(MPI_Comm comm) : pimpl(new Impl(comm)) { }

        // only myself will use this
        Group(const std::shared_ptr<Impl>& p) : pimpl(p) { }

        std::shared_ptr<Impl> pimpl;
    }; // class Group

    ///  Wrapper around MPI_Comm. Has a shallow copy constructor; use Create(Get_group()) for deep copy
    class Intracomm {

        static bool Comm_compare(const MPI_Comm& comm1, const MPI_Comm& comm2) {
            int compare_result;
            const int result = MPI_Comm_compare(comm1, comm2, &compare_result);
            return ((result == MPI_SUCCESS) && (compare_result == MPI_IDENT));
        }

        struct Impl {
            MPI_Comm comm;
            int me;
            int numproc;
            bool owner;

            int utag; // Only used by main thread so no volatile or mutex needed
            int urtag;// Ditto

            Impl(const MPI_Comm& c, int m, int n, bool o) :
                comm(c), me(m), numproc(n), owner(o), utag(1024), urtag(1)
            { MADNESS_ASSERT(comm != MPI_COMM_NULL); }

            ~Impl() {
                if(owner && Is_initialized() && !Is_finalized() && !Comm_compare(comm, MPI_COMM_WORLD) && comm != MPI_COMM_NULL) {
                    MPI_Comm_free(&comm);
                }
            }

            /// Returns a unique tag for temporary use (1023<tag<=4095)

            /// These tags are intended for one time use to avoid tag
            /// collisions with other messages around the same time period.
            /// It simply increments/wraps a counter and returns the next
            /// legal value.
            ///
            /// So that send and receiver agree on the tag all processes
            /// need to call this routine in the same sequence.
            int unique_tag() {
                // RJH removed mutex since ordering requirement across processes means
                // there can never be any thread contention.
                // Cannot use MPI mutex for anything else!
                // It will preprocess to nothing for MPI_THREAD_MULTIPLE!
                //madness::ScopedMutex<madness::SCALABLE_MUTEX_TYPE> obolus(SafeMPI::charon);
                int result = utag++;
                if (utag >= 4095) utag = 1024;
                return result;
            }

            /// \return the period of repeat of unique tags produces by unique_tag()
            static int unique_tag_period() {
              const auto min_tag_value = 1024;
              const auto max_tag_value = 4094;
              return max_tag_value - min_tag_value + 1;
            }

            /// Returns a unique tag reserved for long-term use (0<tag<1000)

            /// Get a tag from this routine for long-term/repeated use.
            ///
            /// Tags in [1000,1023] are statically assigned.
            int unique_reserved_tag() {
                // RJH removed mutex since ordering requirement across processes means
                // Cannot use MPI mutex for anything else!
                // It will preprocess to nothing for MPI_THREAD_MULTIPLE!
                // madness::ScopedMutex<madness::SCALABLE_MUTEX_TYPE> obolus(SafeMPI::charon);
                int result = urtag++;
                if (result >= 1000) MADNESS_EXCEPTION( "too many reserved tags in use" , result );
                return result;
            }

        };
        std::shared_ptr<Impl> pimpl;

        friend void SafeMPI::detail::init_comm_world();
        friend int Finalize();

        // For internal use only. Do not try to call this constructor. It is
        // only used to construct Intarcomm in Create().
        Intracomm(const std::shared_ptr<Impl>& i) : pimpl(i) { }

        // Not allowed
        Intracomm& operator=(const Intracomm& other);

        // makes an uninitialized ptr
        Intracomm() : pimpl(nullptr) {}

    public:
        struct WorldInitObject;

        // For internal use only. Do not try to call this constructor. It is
        // only used to construct COMM_WORLD.
        Intracomm(const WorldInitObject&);

        explicit Intracomm(const MPI_Comm& comm, bool take_ownership_of_comm = true) :
            pimpl()
        {
            MADNESS_ASSERT(Is_initialized());
            int rank = -1, size = -1;
            MADNESS_MPI_TEST(MPI_Comm_rank(comm, &rank));
            MADNESS_MPI_TEST(MPI_Comm_size(comm, &size));
            take_ownership_of_comm =
                    take_ownership_of_comm && (! Comm_compare(comm, MPI_COMM_WORLD));
            pimpl.reset(new Impl(comm, rank, size, take_ownership_of_comm));
        }

        Intracomm(const Intracomm& other) : pimpl(other.pimpl) { }

        ~Intracomm() { }

        /**
         * This collective operation creates a new \c Intracomm from an
         * \c Intracomm::Group object. Must be called by all processes that
         * belong to this communicator, but not all must use the same \c group .
         * Thus this \c Intracomm can be partitioned into several \c Intracomm
         * objects with one call.
         *
         * @param group Intracomm::Group describing the Intracomm object to be
         *   created (\c Intracomm::Get_group() and \c Intracomm::Group::Incl() )
         * @return a new Intracomm object
         */
        Intracomm Create(Group group) const {
            MADNESS_ASSERT(pimpl);
            SAFE_MPI_GLOBAL_MUTEX;
            MPI_Comm group_comm;
            MADNESS_MPI_TEST(MPI_Comm_create(pimpl->comm, group.group(), &group_comm));
            int me; MADNESS_MPI_TEST(MPI_Comm_rank(group_comm, &me));
            int nproc; MADNESS_MPI_TEST(MPI_Comm_size(group_comm, &nproc));
            return Intracomm(std::shared_ptr<Impl>(new Impl(group_comm, me, nproc, true)));
        }

        static const int UNDEFINED_COLOR = MPI_UNDEFINED;
        /**
         * This collective operation creates a new \c Intracomm using
         * the MPI_Comm_split. Must be called by all processes that
         * belong to this communicator. Each caller must provide Color of the new Intracomm
         * and Key (this controls the rank within the new Intracomm;
         * ties are broken by the rank in this Intracomm).
         *
         * @param Color Specifies the new Intracomm that the calling process is to be assigned to.
         *              The value of color must be non-negative. If Color=UNDEFINED_COLOR then
         *              an uninitialized Intracomm object will be produced.
         * @param Key The relative rank of the calling process in the group of the new Intracomm.
         *            If omitted, each communicator's ranks will be determined by
         *            the rank in the host communicator.
         * @return a new Intracomm object
         */
        Intracomm Split(int Color, int Key = 0) const {
            MADNESS_ASSERT(pimpl);
            SAFE_MPI_GLOBAL_MUTEX;
            MPI_Comm group_comm;
            MADNESS_MPI_TEST(MPI_Comm_split(pimpl->comm, Color, Key, &group_comm));
            if (group_comm != MPI_COMM_NULL) {
              int me; MADNESS_MPI_TEST(MPI_Comm_rank(group_comm, &me));
              int nproc; MADNESS_MPI_TEST(MPI_Comm_size(group_comm, &nproc));
              return Intracomm(std::shared_ptr<Impl>(new Impl(group_comm, me, nproc, true)));
            }
            else
              return Intracomm();
        }

      static const int UNDEFINED_SPLIT_TYPE = MPI_UNDEFINED;
      static const int SHARED_SPLIT_TYPE = MPI_COMM_TYPE_SHARED;
      /**
       * This collective operation creates a new \c Intracomm using
       * the MPI_Comm_split_type. Must be called by all processes that
       * belong to this communicator. Each caller must provide the split type
       * and Key (this controls the rank within the new Intracomm;
       * ties are broken by the rank in this Intracomm).
       *
       * @param Type  Controls how this Intracomm will be split.
       *              The value can only be UNDEFINED_SPLIT_TYPE or SHARED_SPLIT_TYPE.
       *              If Type=UNDEFINED_SPLIT_TYPE then
       *              an uninitialized Intracomm object will be produced.
       * @param Key The relative rank of the calling process in the group of the new Intracomm.
       *            If omitted, each communicator's ranks will be determined by
       *            the rank in the host communicator.
       * @return a new Intracomm object
       */
      Intracomm Split_type(int Type, int Key = 0) const {
        MADNESS_ASSERT(pimpl);
        SAFE_MPI_GLOBAL_MUTEX;
        MPI_Comm group_comm;
        MPI_Info info;
        MPI_Info_create(&info);
        MADNESS_MPI_TEST(MPI_Comm_split_type(pimpl->comm, Type, Key, info, &group_comm));
        MPI_Info_free(&info);
        if (group_comm != MPI_COMM_NULL) {
          int me; MADNESS_MPI_TEST(MPI_Comm_rank(group_comm, &me));
          int nproc; MADNESS_MPI_TEST(MPI_Comm_size(group_comm, &nproc));
          return Intracomm(std::shared_ptr<Impl>(new Impl(group_comm, me, nproc, true)));
        }
        else
          return Intracomm();
      }

      /**
         * Clones this Intracomm object
         *
         * @return a (deep) copy of this Intracomm object
         */
        Intracomm Clone() const {
          return Create(this->Get_group());
        }

        bool operator==(const Intracomm& other) const {
            return (pimpl == other.pimpl) || ((pimpl && other.pimpl) &&
                    Comm_compare(pimpl->comm, other.pimpl->comm));
        }

        /**
         * This local operation returns the Intracomm::Group object corresponding to this intracommunicator
         * @return the Intracomm::Group object corresponding to this intracommunicator
         */
        Group Get_group() const {
            MADNESS_ASSERT(pimpl);
            SAFE_MPI_GLOBAL_MUTEX;
            Group group(pimpl->comm);
            return group;
        }

        MPI_Comm& Get_mpi_comm() const {
            MADNESS_ASSERT(pimpl);
            return pimpl->comm;
        }

        int Get_rank() const {
            MADNESS_ASSERT(pimpl);
            return pimpl->me;
        }

        int Get_size() const {
            MADNESS_ASSERT(pimpl);
            return pimpl->numproc;
        }

        Request Isend(const void* buf, const int count, const MPI_Datatype datatype, const int dest, const int tag) const {
            MADNESS_ASSERT(pimpl);
            SAFE_MPI_GLOBAL_MUTEX;
            Request request;
            MADNESS_MPI_TEST(MPI_Isend(const_cast<void*>(buf), count, datatype, dest,tag, pimpl->comm, request));
            return request;
        }

        Request Issend(const void* buf, const int count, const MPI_Datatype datatype, const int dest, const int tag) const {
            MADNESS_ASSERT(pimpl);
            SAFE_MPI_GLOBAL_MUTEX;
            Request request;
            MADNESS_MPI_TEST(MPI_Issend(const_cast<void*>(buf), count, datatype, dest,tag, pimpl->comm, request));
            return request;
        }

        Request Irecv(void* buf, const int count, const MPI_Datatype datatype, const int src, const int tag) const {
            MADNESS_ASSERT(pimpl);
            SAFE_MPI_GLOBAL_MUTEX;
            Request request;
            MADNESS_MPI_TEST(MPI_Irecv(buf, count, datatype, src, tag, pimpl->comm, request));
            return request;
        }

        void Send(const void* buf, const int count, const MPI_Datatype datatype, int dest, int tag) const {
            MADNESS_ASSERT(pimpl);
            SAFE_MPI_GLOBAL_MUTEX;
            MADNESS_MPI_TEST(MPI_Ssend(const_cast<void*>(buf), count, datatype, dest, tag, pimpl->comm));
        }

#ifdef MADNESS_USE_BSEND_ACKS
        void Bsend(const void* buf, size_t count, const MPI_Datatype datatype, int dest, int tag) const {
            MADNESS_ASSERT(pimpl);
            SAFE_MPI_GLOBAL_MUTEX;
            if (count>10 || datatype!=MPI_BYTE) MADNESS_EXCEPTION("Bsend: this protocol is only for 1-byte acks", count );
            MADNESS_MPI_TEST(MPI_Bsend(const_cast<void*>(buf), count, datatype, dest, tag, pimpl->comm));
        }
#endif // MADNESS_USE_BSEND_ACKS

        void Recv(void* buf, const int count, const MPI_Datatype datatype, const int source, const int tag, MPI_Status& status) const {
            MADNESS_ASSERT(pimpl);
            SAFE_MPI_GLOBAL_MUTEX;
            MADNESS_MPI_TEST(MPI_Recv(buf, count, datatype, source, tag, pimpl->comm, &status));
        }

        void Recv(void* buf, const int count, const MPI_Datatype datatype, const int source, const int tag) const {
            MADNESS_ASSERT(pimpl);
            SAFE_MPI_GLOBAL_MUTEX;
            MADNESS_MPI_TEST(MPI_Recv(buf, count, datatype, source, tag, pimpl->comm, MPI_STATUS_IGNORE));
        }

        void Bcast(void* buf, size_t count, const MPI_Datatype datatype, const int root) const {
            MADNESS_ASSERT(pimpl);
            SAFE_MPI_GLOBAL_MUTEX;
            MADNESS_MPI_TEST(MPI_Bcast(buf, count, datatype, root, pimpl->comm));
        }

        void Reduce(const void* sendbuf, void* recvbuf, const int count, const MPI_Datatype datatype, const MPI_Op op, const int root) const {
            MADNESS_ASSERT(pimpl);
            SAFE_MPI_GLOBAL_MUTEX;
            MADNESS_MPI_TEST(MPI_Reduce(const_cast<void*>(sendbuf), recvbuf, count, datatype, op, root, pimpl->comm));
        }

        void Allreduce(const void* sendbuf, void* recvbuf, const int count, const MPI_Datatype datatype, const MPI_Op op) const {
            MADNESS_ASSERT(pimpl);
            SAFE_MPI_GLOBAL_MUTEX;
            MADNESS_MPI_TEST(MPI_Allreduce(const_cast<void*>(sendbuf), recvbuf, count, datatype, op, pimpl->comm));
        }
        bool Get_attr(int key, void* value) const {
            MADNESS_ASSERT(pimpl);
            int flag = 0;
            SAFE_MPI_GLOBAL_MUTEX;
            MADNESS_MPI_TEST(MPI_Comm_get_attr(pimpl->comm, key, value, &flag));
            return flag != 0;
        }

        void Abort(int code=1) const {
            // if SIGABRT has a handler, call std::abort to allow it be caught,
            // else call MPI_Abort which does not seem to call abort at all,
            // instead sends SIGTERM followed by SIGKILL
            // if have a custom signal handler for SIGABRT (i.e. we are running under a
            // debugger) then call abort()
            struct sigaction sa;
            auto rc = sigaction(SIGABRT, NULL, &sa);
            if (rc == 0 && sa.sa_handler != SIG_DFL) {
              std::abort();
            } else {
              MPI_Abort(pimpl->comm, code);
            }
        }

        void Barrier() const {
            MADNESS_ASSERT(pimpl);
            SAFE_MPI_GLOBAL_MUTEX;
            MADNESS_MPI_TEST(MPI_Barrier(pimpl->comm));
        }

        /// Returns a unique tag for temporary use (1023<tag<4095)

        /// These tags are intended for one time use to avoid tag
        /// collisions with other messages around the same time period.
        /// It simply increments/wraps a counter and returns the next
        /// legal value.
        ///
        /// So that send and receiver agree on the tag all processes
        /// need to call this routine in the same sequence.
        int unique_tag() {
            MADNESS_ASSERT(pimpl);
            return pimpl->unique_tag();
        }

        /// \return the period of repeat of unique tags produces by unique_tag()
        static int unique_tag_period() {
          return Impl::unique_tag_period();
        }

        /// Returns a unique tag reserved for long-term use (0<tag<1000)

        /// Get a tag from this routine for long-term/repeated use.
        ///
        /// Tags in [1000,1023] are statically assigned.
        int unique_reserved_tag() {
            MADNESS_ASSERT(pimpl);
            return pimpl->unique_reserved_tag();
        }

        /// Construct info about a binary tree with given root

        /// Constructs a binary tree spanning the communicator with
        /// process root as the root of the tree.  Returns the logical
        /// parent and children in the tree of the calling process.  If
        /// there is no parent/child the value -1 will be set.
        void binary_tree_info(int root, int& parent, int& child0, int& child1);

    }; // class Intracomm

    namespace detail {

        /// Initialize SafeMPI::COMM_WORLD
        inline void init_comm_world() {
            MADNESS_MPI_TEST(MPI_Comm_rank(COMM_WORLD.pimpl->comm, & COMM_WORLD.pimpl->me));
            MADNESS_MPI_TEST(MPI_Comm_size(COMM_WORLD.pimpl->comm, & COMM_WORLD.pimpl->numproc));
            MADNESS_MPI_TEST(MPI_Comm_set_errhandler(COMM_WORLD.pimpl->comm, MPI_ERRORS_RETURN));
        }

    }  // namespace detail


    /// Analogous to MPI_Init_thread

    /// \param argc the number of arguments in argv
    /// \param argv the vector of command-line arguments
    /// \param requested the desired thread level
    /// \return provided thread level
    inline int Init_thread(int & argc, char **& argv, int requested) {
        int provided = 0;
        MADNESS_MPI_TEST(MPI_Init_thread(&argc, &argv, requested, &provided));
        detail::init_comm_world();
        return provided;
    }

    /// Analogous to MPI_Init_thread

    /// \param requested the desired thread level
    /// \return provided thread level
    inline int Init_thread(int requested) {
        int argc = 0;
        char** argv = nullptr;
        return SafeMPI::Init_thread(argc, argv, requested);
    }

    /// Analogous to MPI_Init

    /// \param argc The number of arguments in argv
    /// \param argv The vector of command-line arguments
    inline void Init(int &argc, char **&argv) {
        MADNESS_MPI_TEST(MPI_Init(&argc, &argv));
        SafeMPI::detail::init_comm_world();
    }

    /// Analogous to MPI_Init
    inline void Init() {
        int argc = 0;
        char** argv = nullptr;
        SafeMPI::Init(argc,argv);
    }

    /// Analogous to MPI_Finalize

    /// This returns status rather than throw an
    /// exception upon failure because this is a "destructor", and throwing from
    /// destructors is evil.
    /// \return 0 if successful, nonzero otherwise (see MPI_Finalize() for the
    /// return codes).
    inline int Finalize() {
        SafeMPI::COMM_WORLD.pimpl.reset();
        const int result = MPI_Finalize();
        return result;
    }

    /// Analogous to MPI_Query_thread

    /// \return the MPI thread level provided by SafeMPI::Init_thread()
    inline int Query_thread() {
        int provided;
        MADNESS_MPI_TEST(MPI_Query_thread(&provided));
        return provided;
    }

    /// Wall time

    /// \return The current wall time
    inline double Wtime() { return MPI_Wtime(); }

    /// Set buffer for \c Bsend .

    /// \param buffer The buffer to be used by Bsend
    /// \param size The size of the buffer in Bytes
    inline void Attach_buffer(void* buffer, int size) {
        MADNESS_MPI_TEST(MPI_Buffer_attach(buffer, size));
    }

    /// Unset the \c Bsend buffer.

    /// \param[out] buffer The buffer that was used by Bsend
    inline int Detach_buffer(void *&buffer) {
        int size = 0;
        MPI_Buffer_detach(&buffer, &size);
        return size;
    }

    /// Analogous to MPI_Op_create
    inline MPI_Op Op_create(MPI_User_function *user_fn, int commute) {
      MPI_Op result;
      MADNESS_MPI_TEST(MPI_Op_create(user_fn, commute, &result));
      return result;
    }

    /// Analogous to MPI_Op_free
    inline void Op_free(MPI_Op op) {
      MADNESS_MPI_TEST(MPI_Op_free(&op));
    }

} // namespace SafeMPI

#endif // MADNESS_WORLD_SAFEMPI_H__INCLUDED
