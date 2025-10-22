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
*/

#ifndef MADNESS_WORLD_WORLDDC_H__INCLUDED
#define MADNESS_WORLD_WORLDDC_H__INCLUDED

/*!
  \file worlddc.h
  \brief Implements WorldContainer
  \addtogroup worlddc
  @{

*/

#include <functional>
#include <set>
#include <unordered_set>


#include <madness/world/parallel_archive.h>
#include <madness/world/worldhashmap.h>
#include <madness/world/mpi_archive.h>
#include <madness/world/world_object.h>
#include <madness/world/ranks_and_hosts.h>

namespace madness
{

    template <typename keyT, typename valueT, typename hashfunT>
    class WorldContainer;

    template <typename keyT, typename valueT, typename hashfunT>
    class WorldContainerImpl;

    template <typename keyT, typename valueT, typename hashfunT>
    void swap(WorldContainer<keyT, valueT, hashfunT> &, WorldContainer<keyT, valueT, hashfunT> &);

    template <typename keyT>
    class WorldDCPmapInterface;

    template <typename keyT>
    class WorldDCRedistributeInterface
    {
    public:
        virtual std::size_t size() const = 0;
        virtual void redistribute_phase1(const std::shared_ptr<WorldDCPmapInterface<keyT>> &newmap) = 0;
        virtual void redistribute_phase2() = 0;
        virtual void redistribute_phase3() = 0;
        virtual ~WorldDCRedistributeInterface() {};
    };

    /// some introspection of how data is distributed
    enum DistributionType {
        Distributed,            ///< no replication of the container, the container is distributed over the world
        RankReplicated,         ///< replicate the container over all world ranks
        NodeReplicated          ///< replicate the container over all hosts (compute nodes), once per node,
    ///< even if there are several ranks per node
    };

    template<typename T=long>
    std::ostream& operator<<(std::ostream& os, const DistributionType type) {
        if (type==DistributionType::Distributed) os << "Distributed";
        if (type==DistributionType::RankReplicated) os << "RankReplicated";
        if (type==DistributionType::NodeReplicated) os << "NodeReplicated";
        return os;
    }

    template<typename T=long>
    std::string to_string(const DistributionType type) {
        std::stringstream ss; ss << type;
        return ss.str();
    }

    /// Convert string to DistributionType with case/separator normalization
    inline DistributionType distribution_type_from_string(std::string s) {
        // Normalize: lowercase and replace spaces/dashes with underscores
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
        std::replace(s.begin(), s.end(), ' ', '_');
        std::replace(s.begin(), s.end(), '-', '_');

        if (s == "distributed") return Distributed;
        if (s == "rank_replicated" || s == "rankreplicated" || s == "rank") return RankReplicated;
        if (s == "node_replicated" || s == "nodereplicated" || s == "node" || 
            s == "host_replicated" || s == "hostreplicated" || s == "host") return NodeReplicated;

        throw std::invalid_argument("Unknown DistributionType: " + s);
    }

    /// Wrapper for implicit conversion from string to DistributionType
    struct DistributionTypeFromString {
        DistributionType value;
        DistributionTypeFromString(const std::string& s) : value(distribution_type_from_string(s)) {}
        DistributionTypeFromString(const char* s) : value(distribution_type_from_string(std::string(s))) {}
        operator DistributionType() const noexcept { return value; }
    };




    /// Interface to be provided by any process map

    /// NOTE: if the map is not distributed, but replicated, you must override the distribution_type() method.
    /// \ingroup worlddc
    template <typename keyT>
    class WorldDCPmapInterface
    {
    public:
        typedef WorldDCRedistributeInterface<keyT> *ptrT;


    private:
        std::set<ptrT> ptrs;

    public:
        /// Maps key to processor

        /// @param[in] key Key for container
        /// @return Processor that logically owns the key
        virtual ProcessID owner(const keyT &key) const = 0;

        virtual ~WorldDCPmapInterface() {}

        virtual void print() const {}

        /// by default the map is distributed
        virtual DistributionType distribution_type() const
        {
            return Distributed;
        }

        /// Registers object for receipt of redistribute callbacks

        /// @param[in] ptr Pointer to class derived from WorldDCRedistributedInterface
        void register_callback(ptrT ptr)
        {
            ptrs.insert(ptr);
        }

        /// Deregisters object for receipt of redistribute callbacks

        /// @param[in] ptr Pointer to class derived from WorldDCRedistributedInterface
        void deregister_callback(ptrT ptr)
        {
            ptrs.erase(ptr);
        }

        /// Invoking this switches all registered objects from this process map to the new one

        /// After invoking this routine all objects will be registered with the
        /// new map and no objects will be registered in the current map.
        /// @param[in] world The associated world
        /// @param[in] newpmap The new process map
        void redistribute(World &world, const std::shared_ptr<WorldDCPmapInterface<keyT>> &newpmap)
        {
            print_data_sizes(world, "before redistributing");
            world.gop.fence();
            for (typename std::set<ptrT>::iterator iter = ptrs.begin();
                 iter != ptrs.end();
                 ++iter)
            {
                (*iter)->redistribute_phase1(newpmap);
            }
            world.gop.fence();
            for (typename std::set<ptrT>::iterator iter = ptrs.begin();
                 iter != ptrs.end();
                 ++iter)
            {
                (*iter)->redistribute_phase2();
                newpmap->register_callback(*iter);
            }
            world.gop.fence();
            for (typename std::set<ptrT>::iterator iter = ptrs.begin();
                 iter != ptrs.end();
                 ++iter)
            {
                (*iter)->redistribute_phase3();
            }
            world.gop.fence();
            ptrs.clear();
            newpmap->print_data_sizes(world, "after redistributing");
        }

        /// Counts global number of entries in all containers associated with this process map

        /// Collective operation with global fence
        std::size_t global_size(World &world) const
        {
            world.gop.fence();
            std::size_t sum = local_size();
            world.gop.sum(sum);
            world.gop.fence();
            return sum;
        }

        /// Counts local number of entries in all containers associated with this process map
        std::size_t local_size() const
        {
            std::size_t sum = 0;
            for (typename std::set<ptrT>::iterator iter = ptrs.begin(); iter != ptrs.end(); ++iter)
            {
                sum += (*iter)->size();
            }
            return sum;
        }

        /// Prints size info to std::cout

        /// Collective operation with global fence
        void print_data_sizes(World &world, const std::string msg = "") const
        {
            world.gop.fence();
            std::size_t total = global_size(world);
            std::vector<std::size_t> sizes(world.size());
            sizes[world.rank()] = local_size();
            world.gop.sum(&sizes[0], world.size());
            if (world.rank() == 0)
            {
                madness::print("data distribution info", msg);
                madness::print("   total: ", total);
                std::cout << "   procs: ";
                for (int i = 0; i < world.size(); i++)
                    std::cout << sizes[i] << " ";
                std::cout << std::endl;
            }
            world.gop.fence();
        }
    };

    /// Default process map is "random" using madness::hash(key)

    /// \ingroup worlddc
    template <typename keyT, typename hashfunT = Hash<keyT>>
    class WorldDCDefaultPmap : public WorldDCPmapInterface<keyT>
    {
    private:
        const int nproc;
        hashfunT hashfun;

    public:
        WorldDCDefaultPmap(World &world, const hashfunT &hf = hashfunT()) : nproc(world.mpi.nproc()),
                                                                            hashfun(hf)
        {
        }

        ProcessID owner(const keyT &key) const
        {
            if (nproc == 1)
                return 0;
            return hashfun(key) % nproc;
        }
    };

    /// Local process map will always return the current process as owner

    /// \ingroup worlddc
    template <typename keyT, typename hashfunT = Hash<keyT>>
    class WorldDCLocalPmap : public WorldDCPmapInterface<keyT>
    {
    private:
        ProcessID me;

    public:
        WorldDCLocalPmap(World &world) : me(world.rank()) {}
        ProcessID owner(const keyT &key) const override
        {
            return me;
        }

        DistributionType distribution_type() const override
        {
            return RankReplicated;
        }
    };

    /// node-replicated map will return the lowest rank on the node as owner
    ///
    /// \ingroup worlddc
    template <typename keyT, typename hashfunT = Hash<keyT>>
    class WorldDCNodeReplicatedPmap : public WorldDCPmapInterface<keyT> {
        ProcessID myowner;

    public:
        /// ctor makes a map of all ranks to their owners (lowest rank on the host)
        /// calls a fence
        /// @param[in] world the associated world
        explicit WorldDCNodeReplicatedPmap(World& world, const std::map<std::string,std::vector<long>> ranks_per_host) {
            myowner=lowest_rank_on_host_of_rank(ranks_per_host,world.rank());
        }

        /// owner is the lowest rank on the node, same for all keys
        ProcessID owner(const keyT &key) const override {
            return myowner;
        }

        DistributionType distribution_type() const override
        {
            return NodeReplicated;
        }
    };

    /// check distribution type of WorldContainer -- global communication
    template <typename dcT>
    DistributionType validate_distribution_type(const dcT& dc)
    {
        // assume pmap distribution type is the correct result
        World& world=dc.get_world();
        auto result= dc.get_pmap()->distribution_type();

        auto local_size=dc.size();
        auto global_size=local_size;
        world.gop.sum(global_size);

        auto number_of_duplicates = [](const std::vector<hashT>& v) {
            std::unordered_map<hashT, size_t> counts;
            size_t duplicates = 0;
            for (const auto& elem : v) {
                if (++counts[elem] > 1) duplicates++;
            }
            return duplicates;
        };

        // collect all hash vales and determine the number of duplicates
        std::vector<hashT> all_hashes(local_size);
        const auto& hashfun=dc.get_hash();
        int i=0;
        for (auto it=dc.begin();it!=dc.end();++it,++i) all_hashes[i]=hashfun(it->first);
        all_hashes=world.gop.concat0(all_hashes);
        std::size_t ndup=number_of_duplicates(all_hashes);
        world.gop.broadcast(ndup,0);
        // print("rank, local, global, duplicates", world.rank(),local_size,global_size,ndup);

        // consistency checks
        if (result==Distributed) {
            // all keys should exist only on one process
            MADNESS_CHECK_THROW(ndup==0,"WorldDC inconsistent -- distributed has duplicates");
        }
        else if (result==RankReplicated) {
            // all keys should exist on all processes
            std::size_t nrank=world.size();
            MADNESS_CHECK_THROW(global_size==nrank*local_size,"WorldDC inconsistent");
            MADNESS_CHECK_THROW(ndup==local_size*(nrank-1),"WorldDC inconsistent - duplicates");
        }
        else if (result==NodeReplicated) {
            // all keys should exist on all nodes, not on all ranks
            auto ranks_per_host1=ranks_per_host(world);
            std::vector<int> primary_ranks=primary_ranks_per_host(world,ranks_per_host1);
            world.gop.broadcast_serializable(primary_ranks,0);
            world.gop.fence();
            print("primary_rank_per_host",primary_ranks);

            std::size_t nnodes=primary_ranks.size();
            bool is_primary=(std::find(primary_ranks.begin(),primary_ranks.end(),world.rank())!=primary_ranks.end());

            if (is_primary) {
                MADNESS_CHECK_THROW(global_size==(nnodes*local_size),"WorldDC inconsistent - global size");
                MADNESS_CHECK_THROW(ndup==local_size*(nnodes-1),"WorldDC inconsistent - duplicates");
            } else {
                MADNESS_CHECK_THROW(local_size==0,"WorldDC inconsistent -- secondary");
            }
        }
        return result;
    }


    /// Iterator for distributed container wraps the local iterator

    /// \ingroup worlddc
    template <class internal_iteratorT>
    class WorldContainerIterator
    {
    public:
        typedef typename std::iterator_traits<internal_iteratorT>::iterator_category iterator_category;
        typedef typename std::iterator_traits<internal_iteratorT>::value_type value_type;
        typedef typename std::iterator_traits<internal_iteratorT>::difference_type difference_type;
        typedef typename std::iterator_traits<internal_iteratorT>::pointer pointer;
        typedef typename std::iterator_traits<internal_iteratorT>::reference reference;

    private:
        internal_iteratorT it; ///< Iterator from local container
        // TODO: Convert this to a scoped pointer.
        mutable value_type *value; ///< holds the remote values

    public:
        /// Default constructor makes a local uninitialized value
        explicit WorldContainerIterator()
            : it(), value(nullptr) {}

        /// Initializes from a local iterator
        explicit WorldContainerIterator(const internal_iteratorT &it)
            : it(it), value(nullptr) {}

        /// Initializes to cache a remote value
        explicit WorldContainerIterator(const value_type &v)
            : it(), value(nullptr)
        {
            value = new value_type(v);
        }

        WorldContainerIterator(const WorldContainerIterator &other)
            : it(), value(nullptr)
        {
            copy(other);
        }

        template <class iteratorT>
        WorldContainerIterator(const WorldContainerIterator<iteratorT> &other)
            : it(), value(nullptr)
        {
            copy(other);
        }

        ~WorldContainerIterator()
        {
            delete value;
        }

        /// Assignment
        WorldContainerIterator &operator=(const WorldContainerIterator &other)
        {
            copy(other);
            return *this;
        }

        /// Determines if two iterators are identical
        bool operator==(const WorldContainerIterator &other) const
        {
            return (((!is_cached()) && (!other.is_cached())) && it == other.it) ||
                   ((is_cached() && other.is_cached()) && value->first == other.value->first);
        }

        /// Determines if two iterators are different
        bool operator!=(const WorldContainerIterator &other) const
        {
            return !(*this == other);
        }

        /// Pre-increment of an iterator (i.e., ++it) --- \em local iterators only

        /// Trying to increment a remote iterator will throw
        WorldContainerIterator &operator++()
        {
            MADNESS_ASSERT(!is_cached());
            ++it;
            return *this;
        }

        WorldContainerIterator operator++(int)
        {
            MADNESS_ASSERT(!is_cached());
            WorldContainerIterator<internal_iteratorT> result(*this);
            ++it;
            return result;
        }

        /// Iterators dereference to std::pair<const keyT,valueT>
        pointer operator->() const
        {
            return (is_cached() ? value : it.operator->());
        }

        /// Iterators dereference to std::pair<const keyT,valueT>
        reference operator*() const
        {
            return (is_cached() ? *value : *it);
        }

        /// Private: (or should be) Returns iterator of internal container
        const internal_iteratorT &get_internal_iterator() const
        {
            return it;
        }

        /// Returns true if this is non-local or cached value
        bool is_cached() const
        {
            return value != nullptr;
        }

        template <typename Archive>
        void serialize(const Archive &)
        {
            MADNESS_EXCEPTION("Serializing DC iterator ... why?", false);
        }

    private:
        template <class iteratorT>
        friend class WorldContainerIterator;

        template <class iteratorT>
        void copy(const WorldContainerIterator<iteratorT> &other)
        {
            if (static_cast<const void *>(this) != static_cast<const void *>(&other))
            {
                delete value;
                if (other.is_cached())
                {
                    value = new value_type(*other.value);
                    it = internal_iteratorT();
                }
                else
                {
                    it = other.it;
                    value = nullptr;
                }
            }
        }
    };

    /// Internal implementation of distributed container to facilitate shallow copy

    /// \ingroup worlddc
    template <typename keyT, typename valueT, typename hashfunT>
    class WorldContainerImpl
        : public WorldObject<WorldContainerImpl<keyT, valueT, hashfunT>>,
          public WorldDCRedistributeInterface<keyT>
#ifndef MADNESS_DISABLE_SHARED_FROM_THIS
        ,
          public std::enable_shared_from_this<WorldContainerImpl<keyT, valueT, hashfunT>>
#endif // MADNESS_DISABLE_SHARED_FROM_THIS
    {
    public:
        typedef typename std::pair<const keyT, valueT> pairT;
        typedef const pairT const_pairT;
        typedef WorldContainerImpl<keyT, valueT, hashfunT> implT;

        typedef ConcurrentHashMap<keyT, valueT, hashfunT> internal_containerT;

        // typedef WorldObject< WorldContainerImpl<keyT, valueT, hashfunT> > worldobjT;

        typedef typename internal_containerT::iterator internal_iteratorT;
        typedef typename internal_containerT::const_iterator internal_const_iteratorT;
        typedef typename internal_containerT::accessor accessor;
        typedef typename internal_containerT::const_accessor const_accessor;
        typedef WorldContainerIterator<internal_iteratorT> iteratorT;
        typedef WorldContainerIterator<internal_iteratorT> iterator;
        typedef WorldContainerIterator<internal_const_iteratorT> const_iteratorT;
        typedef WorldContainerIterator<internal_const_iteratorT> const_iterator;

        friend class WorldContainer<keyT, valueT, hashfunT>;

        //         template <typename containerT, typename datumT>
        //         inline
        //         static
        //         typename containerT::iterator replace(containerT& c, const datumT& d) {
        //             std::pair<typename containerT::iterator,bool> p = c.insert(d);
        //             if (!p.second) p.first->second = d.second;   // Who's on first?
        //             return p.first;
        //         }

    private:
        WorldContainerImpl(); // Inhibit default constructor

        std::shared_ptr<WorldDCPmapInterface<keyT>> pmap; ///< Function/class to map from keys to owning process
        const ProcessID me;                               ///< My MPI rank
        internal_containerT local;                        ///< Locally owned data
        std::vector<keyT> *move_list;                     ///< Tempoary used to record data that needs redistributing

        /// Handles find request
        void find_handler(ProcessID requestor, const keyT &key, const RemoteReference<FutureImpl<iterator>> &ref)
        {
            internal_iteratorT r = local.find(key);
            if (r == local.end())
            {
                // print("find_handler: failure:", key);
                this->send(requestor, &implT::find_failure_handler, ref);
            }
            else
            {
                // print("find_handler: success:", key, r->first, r->second);
                this->send(requestor, &implT::find_success_handler, ref, *r);
            }
        }

        /// Handles successful find response
        void find_success_handler(const RemoteReference<FutureImpl<iterator>> &ref, const pairT &datum)
        {
            FutureImpl<iterator> *f = ref.get();
            f->set(iterator(datum));
            // print("find_success_handler: success:", datum.first, datum.second, f->get()->first, f->get()->second);
            //  Todo: Look at this again.
            //            ref.reset(); // Matching inc() in find() where ref was made
        }

        /// Handles unsuccessful find response
        void find_failure_handler(const RemoteReference<FutureImpl<iterator>> &ref)
        {
            FutureImpl<iterator> *f = ref.get();
            f->set(end());
            // print("find_failure_handler");
            //  Todo: Look at this again.
            //            ref.reset(); // Matching inc() in find() where ref was made
        }

    public:
        WorldContainerImpl(World &world,
                           const std::shared_ptr<WorldDCPmapInterface<keyT>> &pm,
                           const hashfunT &hf)
            : WorldObject<WorldContainerImpl<keyT, valueT, hashfunT>>(world), pmap(pm), me(world.mpi.rank()), local(5011, hf)
        {
            pmap->register_callback(this);
        }

        virtual ~WorldContainerImpl()
        {
            pmap->deregister_callback(this);
        }

        const std::shared_ptr<WorldDCPmapInterface<keyT>> &get_pmap() const
        {
            return pmap;
        }

        std::shared_ptr<WorldDCPmapInterface<keyT>> &get_pmap()
        {
            return pmap;
        }

        void reset_pmap_to_local()
        {
            pmap->deregister_callback(this);
            pmap.reset(new WorldDCLocalPmap<keyT>(this->get_world()));
            pmap->register_callback(this);
        }

        /// replicates this WorldContainer on all ranks (ProcessIDs) and generates a
        /// ProcessMap where all ranks are local
        void replicate(bool fence) {
            World &world = this->get_world();
            pmap->deregister_callback(this);
            pmap.reset(new WorldDCLocalPmap<keyT>(world));
            pmap->register_callback(this);

            do_replicate(world);
            if (fence) world.gop.fence();
        }

        /// replicates this WorldContainer on all nodes (hosts) and generates a
        /// ProcessMap where all nodes are node-local (not rank-local)
        /// will always fence
        void replicate_on_hosts(bool fence) {
            MADNESS_CHECK(fence);

            /// print in rank-order
//            auto oprint = [&](World& world, auto &&... args) {
//                world.gop.fence();
//                for (int r=0; r<world.size(); ++r) {
//                    if (r==world.rank()) {
//                        std::cout << "rank " << world.rank() << ": ";
//                        print(std::forward<decltype(args)>(args)...);
//                    }
//                    world.gop.fence();
//                }
//            };

            World &world = this->get_world();

            // find primary ranks per host (lowest rank on each host)
            auto ranks_per_host1=ranks_per_host(world);
            std::vector<int> primary_ranks=primary_ranks_per_host(world,ranks_per_host1);
            world.gop.broadcast_serializable(primary_ranks,0);
            world.gop.fence();

//            auto sizes =[&](std::string msg) {
//                world.gop.fence();
//                auto local_size=size();
//                auto global_size=local_size;
//                world.gop.sum(global_size);
//                oprint(world,"rank, local, global",msg, world.rank(),local_size,global_size);
//                world.gop.fence();
//            };

            // change pmap to replicated
            pmap->deregister_callback(this);
            pmap.reset(new WorldDCLocalPmap<keyT>(world));
            pmap->register_callback(this);

            // shortcut: replace pmap and return
            if (world.size()==1) {
                // change pmap to node replicated
                pmap->deregister_callback(this);
                pmap.reset(new WorldDCNodeReplicatedPmap<keyT>(world,ranks_per_host1));
                pmap->register_callback(this);
                return;
            }
            world.gop.fence();

            // get a list of all other ranks that are not primary
            std::vector<int> secondary_ranks;
            for (int r=0; r<world.size(); ++r) {
                if (std::find(primary_ranks.begin(),primary_ranks.end(),r)==primary_ranks.end())
                    secondary_ranks.push_back(r);
            }

            // phase 1: for all ranks send data to the lowest rank on host

            // step 1-2: send data to lowest rank on host (which will become the owner)
            long myowner = lowest_rank_on_host_of_rank(ranks_per_host1, world.rank());
            // oprint(world,"my owner, size:", myowner,size());
            if (world.rank() != myowner) {
                // send data to myowner
                for (auto it = begin(); it != end(); ++it) {
                    keyT key = it->first;
                    valueT value = it->second;
                    this->send(myowner,&implT::insert,pairT(key,value));
                    // insert(pairT(key,value));        // this won't work with LocalPmap
                }
                // remove all local data after sending
                // clear();
            }
            // need a fence here to make sure send is finished
            world.gop.fence();
            // sizes("after step 1, before clear");
            // world.gop.fence();
            if (world.rank()!=myowner) clear();
            world.gop.fence();
            // sizes("after step 1");

            // change pmap to replicated
            pmap->deregister_callback(this);
            pmap.reset(new WorldDCNodeReplicatedPmap<keyT>(world,ranks_per_host1));
            pmap->register_callback(this);

            // check if this rank is in the primary list
            bool i_am_in_primary_list=world.rank()==myowner;
            if (i_am_in_primary_list) {

                // step 2-1: create a world with only the primary ranks and replicate there (see test_world.cc)

                SafeMPI::Group primary_group = world.mpi.comm().Get_group().Incl(primary_ranks.size(), &primary_ranks[0]);
                SafeMPI::Intracomm comm_primary = world.mpi.comm().Create(primary_group);
                // step 2-2: replicate in the primary world
                {
                    World world_primary(comm_primary);
                    // auto ranks_per_host1=ranks_per_host(world_primary);
                    // if (world_primary.rank()==0) {
                    //     print("host/rank map in primary world:");
                    //     for (auto& p : ranks_per_host1) print(p.first, p.second);
                    // }

                    world_primary.gop.fence();              // this fence seems necessary, why??
                    do_replicate(world_primary);
                    world_primary.gop.fence();
                }
            } else {
                // need this to avoid deadlock in MPI_Comm_create (why??)
                SafeMPI::Group secondary_group = world.mpi.comm().Get_group().Incl(secondary_ranks.size(), &secondary_ranks[0]);
                SafeMPI::Intracomm comm_secondary = world.mpi.comm().Create(secondary_group);
            }

            // phase 3: done
            if (fence) world.gop.fence();
            // validate_distribution_type(*this);
        }

        void do_replicate(World& world) {
            for (ProcessID rank = 0; rank < world.size(); rank++)
            {
                if (rank == world.rank())
                {
                    std::size_t sz = size();
                    world.gop.broadcast_serializable(sz, rank);

                    for (auto it = begin(); it != end(); ++it)
                    {
                        keyT key = it->first;
                        valueT value = it->second;
                        world.gop.broadcast_serializable(key, rank);
                        world.gop.broadcast_serializable(value, rank);
                    }
                }
                else
                {
                    size_t sz;
                    world.gop.broadcast_serializable(sz, rank);
                    for (size_t i = 0; i < sz; i++)
                    {
                        keyT key;
                        valueT value;
                        world.gop.broadcast_serializable(key, rank);
                        world.gop.broadcast_serializable(value, rank);
                        insert(pairT(key, value));
                    }
                }
            }
        }

        const hashfunT &get_hash() const { return local.get_hash(); }

        bool is_local(const keyT &key) const
        {
            return owner(key) == me;
        }

        ProcessID owner(const keyT &key) const
        {
            return pmap->owner(key);
        }

        bool probe(const keyT &key) const
        {
            ProcessID dest = owner(key);
            if (dest == me)
                return local.find(key) != local.end();
            else
                return false;
        }

        std::size_t size() const
        {
            return local.size();
        }

        void insert(const pairT &datum)
        {
            ProcessID dest = owner(datum.first);
            if (dest == me)
            {
                // Was using iterator ... try accessor ?????
                accessor acc;
                // N.B. key might already exist if want to simply replace
                [[maybe_unused]] auto inserted = local.insert(acc, datum.first);
                acc->second = datum.second;
            }
            else
            {
                // Must be send (not task) for sequential consistency (and relies on single-threaded remote server)
                this->send(dest, &implT::insert, datum);
            }
        }

        bool insert_acc(accessor &acc, const keyT &key)
        {
            MADNESS_ASSERT(owner(key) == me);
            return local.insert(acc, key);
        }

        bool insert_const_acc(const_accessor &acc, const keyT &key)
        {
            MADNESS_ASSERT(owner(key) == me);
            return local.insert(acc, key);
        }

        void clear()
        {
            local.clear();
        }

        void erase(const keyT &key)
        {
            ProcessID dest = owner(key);
            if (dest == me)
            {
                [[maybe_unused]] auto erased = local.try_erase(key);
                MADNESS_ASSERT(erased);
            }
            else
            {
                void (implT::*eraser)(const keyT &) = &implT::erase;
                this->send(dest, eraser, key);
            }
        }

        template <typename InIter>
        void erase(InIter it)
        {
            MADNESS_ASSERT(!it.is_cached());
            MADNESS_ASSERT(it != end());
            erase(it->first);
        }

        template <typename InIter>
        void erase(InIter first, InIter last)
        {
            InIter it = first;
            do
            {
                first++;
                erase(it->first);
                it = first;
            } while (first != last);
        }

        iterator begin()
        {
            return iterator(local.begin());
        }

        const_iterator begin() const
        {
            return const_iterator(local.begin());
        }

        iterator end()
        {
            return iterator(local.end());
        }

        const_iterator end() const
        {
            return const_iterator(local.end());
        }

        Future<const_iterator> find(const keyT &key) const
        {
            // Ugliness here to avoid replicating find() and
            // associated handlers for const.  Assumption is that
            // const and non-const iterators are identical except for
            // const attribute ... at some point probably need to do
            // the right thing.
            Future<iterator> r = const_cast<implT *>(this)->find(key);
            return *(Future<const_iterator> *)(&r);
        }

        Future<iterator> find(const keyT &key)
        {
            ProcessID dest = owner(key);
            if (dest == me)
            {
                return Future<iterator>(iterator(local.find(key)));
            }
            else
            {
                Future<iterator> result;
                this->send(dest, &implT::find_handler, me, key, result.remote_ref(this->get_world()));
                return result;
            }
        }

        bool find(accessor &acc, const keyT &key)
        {
            if (owner(key) != me)
                return false;
            return local.find(acc, key);
        }

        bool find(const_accessor &acc, const keyT &key) const
        {
            if (owner(key) != me)
                return false;
            return local.find(acc, key);
        }

        // Used to forward call to item member function
        template <typename memfunT>
        MEMFUN_RETURNT(memfunT)
        itemfun(const keyT &key, memfunT memfun)
        {
            accessor acc;
            // N.B. key may already exist, this is just to ensure lock is held by acc
            [[maybe_unused]] auto inserted = local.insert(acc, key);
            return (acc->second.*memfun)();
        }

        // Used to forward call to item member function
        template <typename memfunT, typename arg1T>
        MEMFUN_RETURNT(memfunT)
        itemfun(const keyT &key, memfunT memfun, const arg1T &arg1)
        {
            accessor acc;
            // N.B. key may already exist, this is just to ensure lock is held by acc
            [[maybe_unused]] auto inserted = local.insert(acc, key);
            return (acc->second.*memfun)(arg1);
        }

        // Used to forward call to item member function
        template <typename memfunT, typename arg1T, typename arg2T>
        MEMFUN_RETURNT(memfunT)
        itemfun(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2)
        {
            accessor acc;
            // N.B. key may already exist, this is just to ensure lock is held by acc
            [[maybe_unused]] auto inserted = local.insert(acc, key);
            return (acc->second.*memfun)(arg1, arg2);
        }

        // Used to forward call to item member function
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        MEMFUN_RETURNT(memfunT)
        itemfun(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3)
        {
            accessor acc;
            // N.B. key may already exist, this is just to ensure lock is held by acc
            [[maybe_unused]] auto inserted = local.insert(acc, key);
            return (acc->second.*memfun)(arg1, arg2, arg3);
        }

        // Used to forward call to item member function
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        MEMFUN_RETURNT(memfunT)
        itemfun(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3, const arg4T &arg4)
        {
            accessor acc;
            // N.B. key may already exist, this is just to ensure lock is held by acc
            [[maybe_unused]] auto inserted = local.insert(acc, key);
            return (acc->second.*memfun)(arg1, arg2, arg3, arg4);
        }

        // Used to forward call to item member function
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        MEMFUN_RETURNT(memfunT)
        itemfun(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3, const arg4T &arg4, const arg5T &arg5)
        {
            accessor acc;
            // N.B. key may already exist, this is just to ensure lock is held by acc
            [[maybe_unused]] auto inserted = local.insert(acc, key);
            return (acc->second.*memfun)(arg1, arg2, arg3, arg4, arg5);
        }

        // Used to forward call to item member function
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T>
        MEMFUN_RETURNT(memfunT)
        itemfun(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3, const arg4T &arg4, const arg5T &arg5, const arg6T &arg6)
        {
            accessor acc;
            // N.B. key may already exist, this is just to ensure lock is held by acc
            [[maybe_unused]] auto inserted = local.insert(acc, key);
            return (acc->second.*memfun)(arg1, arg2, arg3, arg4, arg5, arg6);
        }

        // Used to forward call to item member function
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T>
        MEMFUN_RETURNT(memfunT)
        itemfun(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3,
                const arg4T &arg4, const arg5T &arg5, const arg6T &arg6, const arg7T &arg7)
        {
            accessor acc;
            // N.B. key may already exist, this is just to ensure lock is held by acc
            [[maybe_unused]] auto inserted = local.insert(acc, key);
            return (acc->second.*memfun)(arg1, arg2, arg3, arg4, arg5, arg6, arg7);
        }

        // First phase of redistributions changes pmap and makes list of stuff to move
        void redistribute_phase1(const std::shared_ptr<WorldDCPmapInterface<keyT>> &newpmap)
        {
            pmap = newpmap;
            move_list = new std::vector<keyT>();
            for (typename internal_containerT::iterator iter = local.begin(); iter != local.end(); ++iter)
            {
                if (owner(iter->first) != me)
                    move_list->push_back(iter->first);
            }
        }

        struct P2Op
        {
            implT *impl;
            typedef Range<typename std::vector<keyT>::const_iterator> rangeT;
            P2Op(implT *impl) : impl(impl) {}
            P2Op(const P2Op &p) : impl(p.impl) {}
            bool operator()(typename rangeT::iterator &iterator) const
            {
                typename internal_containerT::iterator iter = impl->local.find(*iterator);
                MADNESS_ASSERT(iter != impl->local.end());

                // impl->insert(*iter);
                impl->task(impl->owner(*iterator), &implT::insert, *iter);

                impl->local.erase(iter); // delete local copy of the data
                return true;
            }
        };

        // Second phase moves data
        void redistribute_phase2()
        {
            this->get_world().taskq.for_each(typename P2Op::rangeT(move_list->begin(), move_list->end()), P2Op(this));
            // std::vector<keyT>& mvlist = *move_list;
            // for (unsigned int i=0; i<move_list->size(); ++i) {
            //     typename internal_containerT::iterator iter = local.find(mvlist[i]);
            //     MADNESS_ASSERT(iter != local.end());
            //     insert(*iter);
            //     local.erase(iter);
            // }
            // delete move_list;
        }

        // Third phase cleans up
        void redistribute_phase3()
        {
            delete move_list;
        }
    };

    /// Makes a distributed container with specified attributes

    /// \ingroup worlddc
    ///
    /// There is no communication or syncronization associated with
    /// making a new container, but every process must invoke the
    /// constructor for each container in the same order.  This is so
    /// that we can assign each container a unique ID without any
    /// communication.  Since remotely invoked operations may start
    /// happening before local construction, messages on not yet
    /// constructed containers are buffered pending construction.
    ///
    /// Similarly, when a container is destroyed, the actual
    /// destruction is deferred until a synchronization point
    /// (world.gop.fence()) in order to eliminate the need to fence
    /// before destroying every container.
    ///
    /// The distribution of data between processes is controlled by
    /// the process map (Pmap) class.  The default is uniform
    /// hashing based upon a strong (Bob Jenkins, lookup3) bytewise
    /// hash of the key.
    ///
    /// All operations, including constructors and destructors, are
    /// non-blocking and return immediately.  If communication occurs
    /// it is asynchronous, otherwise operations are local.
    template <typename keyT, typename valueT, typename hashfunT = Hash<keyT>>
    class WorldContainer : public archive::ParallelSerializableObject
    {
    public:
        // access keyT and valueT types for serialization
        typedef keyT key_type;
        typedef WorldContainer<keyT, valueT, hashfunT> containerT;
        typedef WorldContainerImpl<keyT, valueT, hashfunT> implT;
        typedef typename implT::pairT pairT;
        typedef typename implT::iterator iterator;
        typedef typename implT::const_iterator const_iterator;
        typedef typename implT::accessor accessor;
        typedef typename implT::const_accessor const_accessor;
        typedef Future<iterator> futureT;
        typedef Future<const_iterator> const_futureT;

    private:
        std::shared_ptr<implT> p;

        inline void check_initialized() const
        {
            MADNESS_ASSERT(p);
        }

    public:
        /// Makes an uninitialized container (no communication)

        /// The container is useless until assigned to from a fully
        /// constructed container.  There is no need to worry about
        /// default constructors being executed in order.
        WorldContainer()
            : p()
        {
        }

        /// Makes an initialized, empty container with default data distribution (no communication)

        /// A unique ID is associated with every distributed container
        /// within a world.  In order to avoid synchronization when
        /// making a container, we have to assume that all processes
        /// execute this constructor in the same order (does not apply
        /// to the non-initializing, default constructor).
        WorldContainer(World &world, bool do_pending = true, const hashfunT &hf = hashfunT())
            : p(new implT(world,
                          std::shared_ptr<WorldDCPmapInterface<keyT>>(new WorldDCDefaultPmap<keyT, hashfunT>(world, hf)),
                          hf))
        {
            if (do_pending)
                p->process_pending();
        }

        /// Makes an initialized, empty container (no communication)

        /// A unique ID is associated with every distributed container
        /// within a world.  In order to avoid synchronization when
        /// making a container, we have to assume that all processes
        /// execute this constructor in the same order (does not apply
        /// to the non-initializing, default constructor).
        WorldContainer(World &world,
                       const std::shared_ptr<WorldDCPmapInterface<keyT>> &pmap,
                       bool do_pending = true,
                       const hashfunT &hf = hashfunT())
            : p(new implT(world, pmap, hf))
        {
            if (do_pending)
                p->process_pending();
        }

        /// Copy constructor is shallow (no communication)

        /// The copy refers to exactly the same container as other
        /// which must be initialized.
        WorldContainer(const WorldContainer &other)
            : p(other.p)
        {
            check_initialized();
        }

        /// Assignment is shallow (no communication)

        /// The copy refers to exactly the same container as other
        /// which must be initialized.
        containerT &operator=(const containerT &other)
        {
            if (this != &other)
            {
                other.check_initialized();
                p = other.p;
            }
            return *this;
        }

        /// return the way data is distributed
        DistributionType get_distribution_type() const {
            if (!p) MADNESS_EXCEPTION("Uninitialized container", false);
            return p->get_pmap()->distribution_type();
        }

        bool is_distributed() const {
            return get_distribution_type()==Distributed;
        }

        // Canonical naming with "rank_replication" and "node_replication"
        bool rank_replication() const {
            return get_distribution_type()==RankReplicated;
        }

        bool node_replication() const {
            return get_distribution_type()==NodeReplicated;
        }

        // Deprecated aliases for backward compatibility
        [[deprecated("Use rank_replication() instead")]]
        bool is_replicated() const {
            return rank_replication();
        }

        [[deprecated("Use node_replication() instead")]]
        bool is_host_replicated() const {
            return node_replication();
        }

        /// Returns the world associated with this container
        World &get_world() const
        {
            check_initialized();
            return p->get_world();
        }

        std::shared_ptr<WorldDCPmapInterface<keyT>> &get_impl()
        {
            check_initialized();
            return p;
        }

        // Canonical naming: replicate_on_ranks and replicate_on_nodes
        /// replicates this WorldContainer on all ProcessIDs (ranks)
        void replicate_on_ranks(bool fence = true)
        {
            p->replicate(fence);
        }

        /// replicates this WorldContainer on all nodes (one PID per node/host)
        void replicate_on_nodes(bool fence = true)
        {
            p->replicate_on_hosts(fence);
        }

        // Deprecated aliases for backward compatibility
        [[deprecated("Use replicate_on_ranks() instead")]]
        void replicate(bool fence = true)
        {
            replicate_on_ranks(fence);
        }

        [[deprecated("Use replicate_on_nodes() instead")]]
        void replicate_on_hosts(bool fence = true)
        {
            replicate_on_nodes(fence);
        }

        /// Inserts/replaces key+value pair (non-blocking communication if key not local)
        void replace(const pairT &datum)
        {
            check_initialized();
            p->insert(datum);
        }

        /// Inserts/replaces key+value pair (non-blocking communication if key not local)
        void replace(const keyT &key, const valueT &value)
        {
            replace(pairT(key, value));
        }

        /// Write access to LOCAL value by key. Returns true if found, false otherwise (always false for remote).
        bool find(accessor &acc, const keyT &key)
        {
            check_initialized();
            return p->find(acc, key);
        }

        /// Read access to LOCAL value by key. Returns true if found, false otherwise (always false for remote).
        bool find(const_accessor &acc, const keyT &key) const
        {
            check_initialized();
            return p->find(acc, key);
        }

        /// Write access to LOCAL value by key. Returns true if inserted, false if already exists (throws if remote)
        bool insert(accessor &acc, const keyT &key)
        {
            check_initialized();
            return p->insert_acc(acc, key);
        }

        /// Read access to LOCAL value by key. Returns true if inserted, false if already exists (throws if remote)
        bool insert(const_accessor &acc, const keyT &key)
        {
            check_initialized();
            return p->insert_acc(acc, key);
        }

        /// Inserts pairs (non-blocking communication if key(s) not local)
        template <typename input_iterator>
        void replace(input_iterator &start, input_iterator &end)
        {
            check_initialized();
            using std::placeholders::_1;
            std::for_each(start, end, std::bind(this, std::mem_fn(&containerT::insert), _1));
        }

        /// Returns true if local data is immediately available (no communication)
        bool probe(const keyT &key) const
        {
            check_initialized();
            return p->probe(key);
        }

        /// Returns processor that logically owns key (no communication)

        /// Local remapping may have changed its physical location, but all
        /// operations should forward correctly.
        inline ProcessID owner(const keyT &key) const
        {
            check_initialized();
            return p->owner(key);
        }

        /// Returns true if the key maps to the local processor (no communication)
        bool is_local(const keyT &key) const
        {
            check_initialized();
            return p->is_local(key);
        }

        /// Returns a future iterator (non-blocking communication if key not local)

        /// Like an std::map an iterator "points" to an std::pair<const keyT,valueT>.
        ///
        /// Refer to Future for info on how to avoid blocking.
        Future<iterator> find(const keyT &key)
        { //
            check_initialized();
            return p->find(key);
        }

        /// Returns a future iterator (non-blocking communication if key not local)

        /// Like an std::map an iterator "points" to an std::pair<const keyT,valueT>.
        ///
        /// Refer to Future for info on how to avoid blocking.
        Future<const_iterator> find(const keyT &key) const
        {
            check_initialized();
            return const_cast<const implT *>(p.get())->find(key);
        }

        /// Returns an iterator to the beginning of the \em local data (no communication)
        iterator begin()
        {
            check_initialized();
            return p->begin();
        }

        /// Returns an iterator to the beginning of the \em local data (no communication)
        const_iterator begin() const
        {
            check_initialized();
            return const_cast<const implT *>(p.get())->begin();
        }

        /// Returns an iterator past the end of the \em local data (no communication)
        iterator end()
        {
            check_initialized();
            return p->end();
        }

        /// Returns an iterator past the end of the \em local data (no communication)
        const_iterator end() const
        {
            check_initialized();
            return const_cast<const implT *>(p.get())->end();
        }

        /// Erases entry from container (non-blocking comm if remote)

        /// Missing keys are quietly ignored.
        ///
        /// Note that erasing an entry may invalidate iterators on the
        /// remote end.  This is just the same as what happens when
        /// using STL iterators on an STL container in a sequential
        /// algorithm.
        void erase(const keyT &key)
        {
            check_initialized();
            p->erase(key);
        }

        /// Erases entry corresponding to \em local iterator (no communication)
        void erase(const iterator &it)
        {
            check_initialized();
            p->erase(it);
        }

        /// Erases range defined by \em local iterators (no communication)
        void erase(const iterator &start, const iterator &finish)
        {
            check_initialized();
            p->erase(start, finish);
        }

        /// Clears all \em local data (no communication)

        /// Invalidates all iterators
        void clear()
        {
            check_initialized();
            p->clear();
        }

        /// Returns the number of \em local entries (no communication)
        std::size_t size() const
        {
            check_initialized();
            return p->size();
        }

        /// Returns shared pointer to the process mapping
        inline const std::shared_ptr<WorldDCPmapInterface<keyT>> &get_pmap() const
        {
            check_initialized();
            return p->get_pmap();
        }

        /// Returns shared pointer to the process mapping
        inline void reset_pmap_to_local()
        {
            p->reset_pmap_to_local();
        }

        /// Returns a reference to the hashing functor
        const hashfunT &get_hash() const
        {
            check_initialized();
            return p->get_hash();
        }

        /// Process pending messages

        /// If the constructor was given \c do_pending=false then you
        /// \em must invoke this routine in order to process both
        /// prior and future messages.
        inline void process_pending()
        {
            check_initialized();
            p->process_pending();
        }

        /// Sends message "resultT memfun()" to item (non-blocking comm if remote)

        /// If item does not exist it is made with the default constructor.
        ///
        /// Future arguments must be ready for remote messages.
        ///
        /// Returns a future result (Future<void> may be ignored).
        ///
        /// The method executes with a write lock on the item.
        template <typename memfunT>
        Future<MEMFUN_RETURNT(memfunT)>
        send(const keyT &key, memfunT memfun)
        {
            check_initialized();
            MEMFUN_RETURNT(memfunT)
            (implT::*itemfun)(const keyT &, memfunT) = &implT::template itemfun<memfunT>;
            return p->send(owner(key), itemfun, key, memfun);
        }

        /// Sends message "resultT memfun(arg1T)" to item (non-blocking comm if remote)

        /// If item does not exist it is made with the default constructor.
        ///
        /// Future arguments must be ready for remote messages.
        ///
        /// Returns a future result (Future<void> may be ignored).
        ///
        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        send(const keyT &key, const memfunT &memfun, const arg1T &arg1)
        {
            check_initialized();
            // To work around bug in g++ 4.3.* use static cast as alternative mechanism to force type deduction
            MEMFUN_RETURNT(memfunT)
            (implT::*itemfun)(const keyT &, memfunT, const arg1T &) = &implT::template itemfun<memfunT, arg1T>;
            return p->send(owner(key), itemfun, key, memfun, arg1);
            /*return p->send(owner(key),
                           static_cast<MEMFUN_RETURNT(memfunT)(implT::*)(const keyT&, memfunT, const arg1T&)>(&implT:: template itemfun<memfunT,arg1T>),
                           key, memfun, arg1);*/
        }

        /// Sends message "resultT memfun(arg1T,arg2T)" to item (non-blocking comm if remote)

        /// If item does not exist it is made with the default constructor.
        ///
        /// Future arguments must be ready for both local and remote messages.
        ///
        /// Returns a future result (Future<void> may be ignored).
        ///
        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        send(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2)
        {
            check_initialized();
            // To work around bug in g++ 4.3.* use static cast as alternative mechanism to force type deduction
            MEMFUN_RETURNT(memfunT)
            (implT::*itemfun)(const keyT &, memfunT, const arg1T &, const arg2T &) = &implT::template itemfun<memfunT, arg1T, arg2T>;
            return p->send(owner(key), itemfun, key, memfun, arg1, arg2);
            /*return p->send(owner(key),
                           static_cast<MEMFUN_RETURNT(memfunT)(implT::*)(const keyT&, memfunT, const arg1T&, const arg2T&)>(&implT:: template itemfun<memfunT,arg1T,arg2T>), key, memfun, arg1, arg2);*/
        }

        /// Sends message "resultT memfun(arg1T,arg2T,arg3T)" to item (non-blocking comm if remote)

        /// If item does not exist it is made with the default constructor.
        ///
        /// Future arguments must be ready for both local and remote messages.
        ///
        /// Returns a future result (Future<void> may be ignored).
        ///
        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        send(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3)
        {
            check_initialized();
            MEMFUN_RETURNT(memfunT)
            (implT::*itemfun)(const keyT &, memfunT, const arg1T &, const arg2T &, const arg3T &) = &implT::template itemfun<memfunT, arg1T, arg2T, arg3T>;
            return p->send(owner(key), itemfun, key, memfun, arg1, arg2, arg3);
        }

        /// Sends message "resultT memfun(arg1T,arg2T,arg3T,arg4T)" to item (non-blocking comm if remote)

        /// If item does not exist it is made with the default constructor.
        ///
        /// Future arguments must be ready for both local and remote messages.
        ///
        /// Returns a future result (Future<void> may be ignored).
        ///
        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        send(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3, const arg4T &arg4)
        {
            check_initialized();
            MEMFUN_RETURNT(memfunT)
            (implT::*itemfun)(const keyT &, memfunT, const arg1T &, const arg2T &, const arg3T &, const arg4T &) = &implT::template itemfun<memfunT, arg1T, arg2T, arg3T, arg4T>;
            return p->send(owner(key), itemfun, key, memfun, arg1, arg2, arg3, arg4);
        }

        /// Sends message "resultT memfun(arg1T,arg2T,arg3T,arg4T,arg5T)" to item (non-blocking comm if remote)

        /// If item does not exist it is made with the default constructor.
        ///
        /// Future arguments must be ready for both local and remote messages.
        ///
        /// Returns a future result (Future<void> may be ignored).
        ///
        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        send(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3, const arg4T &arg4, const arg5T &arg5)
        {
            check_initialized();
            MEMFUN_RETURNT(memfunT)
            (implT::*itemfun)(const keyT &, memfunT, const arg1T &, const arg2T &, const arg3T &, const arg4T &, const arg5T &) = &implT::template itemfun<memfunT, arg1T, arg2T, arg3T, arg4T, arg5T>;
            return p->send(owner(key), itemfun, key, memfun, arg1, arg2, arg3, arg4, arg5);
        }

        /// Sends message "resultT memfun(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T)" to item (non-blocking comm if remote)

        /// If item does not exist it is made with the default constructor.
        ///
        /// Future arguments must be ready for both local and remote messages.
        ///
        /// Returns a future result (Future<void> may be ignored).
        ///
        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        send(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3, const arg4T &arg4, const arg5T &arg5, const arg6T &arg6)
        {
            check_initialized();
            MEMFUN_RETURNT(memfunT)
            (implT::*itemfun)(const keyT &, memfunT, const arg1T &, const arg2T &, const arg3T &, const arg4T &, const arg5T &, const arg6T &) = &implT::template itemfun<memfunT, arg1T, arg2T, arg3T, arg4T, arg5T, arg6T>;
            return p->send(owner(key), itemfun, key, memfun, arg1, arg2, arg3, arg4, arg5, arg6);
        }

        /// Sends message "resultT memfun(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T)" to item (non-blocking comm if remote)

        /// If item does not exist it is made with the default constructor.
        ///
        /// Future arguments must be ready for both local and remote messages.
        ///
        /// Returns a future result (Future<void> may be ignored).
        ///
        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        send(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3, const arg4T &arg4,
             const arg5T &arg5, const arg6T &arg6, const arg7T &arg7)
        {
            check_initialized();
            MEMFUN_RETURNT(memfunT)
            (implT::*itemfun)(const keyT &, memfunT, const arg1T &, const arg2T &, const arg3T &, const arg4T &, const arg5T &, const arg6T &, const arg7T &) = &implT::template itemfun<memfunT, arg1T, arg2T, arg3T, arg4T, arg5T, arg6T, arg7T>;
            return p->send(owner(key), itemfun, key, memfun, arg1, arg2, arg3, arg4, arg5, arg6, arg7);
        }

        /// Sends message "resultT memfun() const" to item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        send(const keyT &key, memfunT memfun) const
        {
            return const_cast<containerT *>(this)->send(key, memfun);
        }

        /// Sends message "resultT memfun(arg1T) const" to item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        send(const keyT &key, memfunT memfun, const arg1T &arg1) const
        {
            return const_cast<containerT *>(this)->send(key, memfun, arg1);
        }

        /// Sends message "resultT memfun(arg1T,arg2T) const" to item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        send(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2) const
        {
            return const_cast<containerT *>(this)->send(key, memfun, arg1, arg2);
        }

        /// Sends message "resultT memfun(arg1T,arg2T,arg3T) const" to item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        send(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3) const
        {
            return const_cast<containerT *>(this)->send(key, memfun, arg1, arg2, arg3);
        }

        /// Sends message "resultT memfun(arg1T,arg2T,arg3T,arg4T) const" to item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        send(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3, const arg4T &arg4) const
        {
            return const_cast<containerT *>(this)->send(key, memfun, arg1, arg2, arg3, arg4);
        }

        /// Sends message "resultT memfun(arg1T,arg2T,arg3T,arg4T,arg5T) const" to item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        send(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3, const arg4T &arg4, const arg5T &arg5) const
        {
            return const_cast<containerT *>(this)->send(key, memfun, arg1, arg2, arg3, arg4, arg5);
        }

        /// Sends message "resultT memfun(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T) const" to item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        send(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3,
             const arg4T &arg4, const arg5T &arg5, const arg6T &arg6) const
        {
            return const_cast<containerT *>(this)->send(key, memfun, arg1, arg2, arg3, arg4, arg5, arg6);
        }

        /// Sends message "resultT memfun(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T) const" to item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        send(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3,
             const arg4T &arg4, const arg5T &arg5, const arg6T &arg6, const arg7T &arg7) const
        {
            return const_cast<containerT *>(this)->send(key, memfun, arg1, arg2, arg3, arg4, arg5, arg6, arg7);
        }

        /// Adds task "resultT memfun()" in process owning item (non-blocking comm if remote)

        /// If item does not exist it is made with the default constructor.
        ///
        /// Future arguments for local tasks can generate dependencies, but for remote
        /// tasks all futures must be ready.
        ///
        /// Returns a future result (Future<void> may be ignored).
        ///
        /// The method executes with a write lock on the item.
        template <typename memfunT>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        task(const keyT &key, memfunT memfun, const TaskAttributes &attr = TaskAttributes())
        {
            check_initialized();
            MEMFUN_RETURNT(memfunT)
            (implT::*itemfun)(const keyT &, memfunT) = &implT::template itemfun<memfunT>;
            return p->task(owner(key), itemfun, key, memfun, attr);
        }

        /// Adds task "resultT memfun(arg1T)" in process owning item (non-blocking comm if remote)

        /// If item does not exist it is made with the default constructor.
        ///
        /// Future arguments for local tasks can generate dependencies, but for remote
        /// tasks all futures must be ready.
        ///
        /// Returns a future result (Future<void> may be ignored).
        ///
        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        task(const keyT &key, memfunT memfun, const arg1T &arg1, const TaskAttributes &attr = TaskAttributes())
        {
            check_initialized();
            typedef REMFUTURE(arg1T) a1T;
            MEMFUN_RETURNT(memfunT)
            (implT::*itemfun)(const keyT &, memfunT, const a1T &) = &implT::template itemfun<memfunT, a1T>;
            return p->task(owner(key), itemfun, key, memfun, arg1, attr);
        }

        /// Adds task "resultT memfun(arg1T,arg2T)" in process owning item (non-blocking comm if remote)

        /// If item does not exist it is made with the default constructor.
        ///
        /// Future arguments for local tasks can generate dependencies, but for remote
        /// tasks all futures must be ready.
        ///
        /// Returns a future result (Future<void> may be ignored).
        ///
        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        task(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const TaskAttributes &attr = TaskAttributes())
        {
            check_initialized();
            typedef REMFUTURE(arg1T) a1T;
            typedef REMFUTURE(arg2T) a2T;
            MEMFUN_RETURNT(memfunT)
            (implT::*itemfun)(const keyT &, memfunT, const a1T &, const a2T &) = &implT::template itemfun<memfunT, a1T, a2T>;
            return p->task(owner(key), itemfun, key, memfun, arg1, arg2, attr);
        }

        /// Adds task "resultT memfun(arg1T,arg2T,arg3T)" in process owning item (non-blocking comm if remote)

        /// If item does not exist it is made with the default constructor.
        ///
        /// Future arguments for local tasks can generate dependencies, but for remote
        /// tasks all futures must be ready.
        ///
        /// Returns a future result (Future<void> may be ignored).
        ///
        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        task(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3, const TaskAttributes &attr = TaskAttributes())
        {
            check_initialized();
            typedef REMFUTURE(arg1T) a1T;
            typedef REMFUTURE(arg2T) a2T;
            typedef REMFUTURE(arg3T) a3T;
            MEMFUN_RETURNT(memfunT)
            (implT::*itemfun)(const keyT &, memfunT, const a1T &, const a2T &, const a3T &) = &implT::template itemfun<memfunT, a1T, a2T, a3T>;
            return p->task(owner(key), itemfun, key, memfun, arg1, arg2, arg3, attr);
        }

        /// Adds task "resultT memfun(arg1T,arg2T,arg3T,arg4T)" in process owning item (non-blocking comm if remote)

        /// If item does not exist it is made with the default constructor.
        ///
        /// Future arguments for local tasks can generate dependencies, but for remote
        /// tasks all futures must be ready.
        ///
        /// Returns a future result (Future<void> may be ignored).
        ///
        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        task(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3, const arg4T &arg4, const TaskAttributes &attr = TaskAttributes())
        {
            check_initialized();
            typedef REMFUTURE(arg1T) a1T;
            typedef REMFUTURE(arg2T) a2T;
            typedef REMFUTURE(arg3T) a3T;
            typedef REMFUTURE(arg4T) a4T;
            MEMFUN_RETURNT(memfunT)
            (implT::*itemfun)(const keyT &, memfunT, const a1T &, const a2T &, const a3T &, const a4T &) = &implT::template itemfun<memfunT, a1T, a2T, a3T, a4T>;
            return p->task(owner(key), itemfun, key, memfun, arg1, arg2, arg3, arg4, attr);
        }

        /// Adds task "resultT memfun(arg1T,arg2T,arg3T,arg4T,arg5T)" in process owning item (non-blocking comm if remote)

        /// If item does not exist it is made with the default constructor.
        ///
        /// Future arguments for local tasks can generate dependencies, but for remote
        /// tasks all futures must be ready.
        ///
        /// Returns a future result (Future<void> may be ignored).
        ///
        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        task(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3, const arg4T &arg4, const arg5T &arg5, const TaskAttributes &attr = TaskAttributes())
        {
            check_initialized();
            typedef REMFUTURE(arg1T) a1T;
            typedef REMFUTURE(arg2T) a2T;
            typedef REMFUTURE(arg3T) a3T;
            typedef REMFUTURE(arg4T) a4T;
            typedef REMFUTURE(arg5T) a5T;
            MEMFUN_RETURNT(memfunT)
            (implT::*itemfun)(const keyT &, memfunT, const a1T &, const a2T &, const a3T &, const a4T &, const a5T &) = &implT::template itemfun<memfunT, a1T, a2T, a3T, a4T, a5T>;
            return p->task(owner(key), itemfun, key, memfun, arg1, arg2, arg3, arg4, arg5, attr);
        }

        /// Adds task "resultT memfun(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T)" in process owning item (non-blocking comm if remote)

        /// If item does not exist it is made with the default constructor.
        ///
        /// Future arguments for local tasks can generate dependencies, but for remote
        /// tasks all futures must be ready.
        ///
        /// Returns a future result (Future<void> may be ignored).
        ///
        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        task(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3, const arg4T &arg4, const arg5T &arg5, const arg6T &arg6, const TaskAttributes &attr = TaskAttributes())
        {
            check_initialized();
            typedef REMFUTURE(arg1T) a1T;
            typedef REMFUTURE(arg2T) a2T;
            typedef REMFUTURE(arg3T) a3T;
            typedef REMFUTURE(arg4T) a4T;
            typedef REMFUTURE(arg5T) a5T;
            typedef REMFUTURE(arg6T) a6T;
            MEMFUN_RETURNT(memfunT)
            (implT::*itemfun)(const keyT &, memfunT, const a1T &, const a2T &, const a3T &, const a4T &, const a5T &, const a6T &) = &implT::template itemfun<memfunT, a1T, a2T, a3T, a4T, a5T, a6T>;
            return p->task(owner(key), itemfun, key, memfun, arg1, arg2, arg3, arg4, arg5, arg6, attr);
        }

        /// Adds task "resultT memfun(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T)" in process owning item (non-blocking comm if remote)

        /// If item does not exist it is made with the default constructor.
        ///
        /// Future arguments for local tasks can generate dependencies, but for remote
        /// tasks all futures must be ready.
        ///
        /// Returns a future result (Future<void> may be ignored).
        ///
        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        task(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3, const arg4T &arg4, const arg5T &arg5, const arg6T &arg6, const arg7T &arg7, const TaskAttributes &attr = TaskAttributes())
        {
            check_initialized();
            typedef REMFUTURE(arg1T) a1T;
            typedef REMFUTURE(arg2T) a2T;
            typedef REMFUTURE(arg3T) a3T;
            typedef REMFUTURE(arg4T) a4T;
            typedef REMFUTURE(arg5T) a5T;
            typedef REMFUTURE(arg6T) a6T;
            typedef REMFUTURE(arg7T) a7T;
            MEMFUN_RETURNT(memfunT)
            (implT::*itemfun)(const keyT &, memfunT, const a1T &, const a2T &, const a3T &, const a4T &, const a5T &, const a6T &, const a7T &) = &implT::template itemfun<memfunT, a1T, a2T, a3T, a4T, a5T, a6T, a7T>;
            return p->task(owner(key), itemfun, key, memfun, arg1, arg2, arg3, arg4, arg5, arg6, arg7, attr);
        }

        /// Adds task "resultT memfun() const" in process owning item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        task(const keyT &key, memfunT memfun, const TaskAttributes &attr = TaskAttributes()) const
        {
            return const_cast<containerT *>(this)->task(key, memfun, attr);
        }

        /// Adds task "resultT memfun(arg1T) const" in process owning item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        task(const keyT &key, memfunT memfun, const arg1T &arg1, const TaskAttributes &attr = TaskAttributes()) const
        {
            return const_cast<containerT *>(this)->task(key, memfun, arg1, attr);
        }

        /// Adds task "resultT memfun(arg1T,arg2T) const" in process owning item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        task(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const TaskAttributes &attr = TaskAttributes()) const
        {
            return const_cast<containerT *>(this)->task(key, memfun, arg1, arg2, attr);
        }

        /// Adds task "resultT memfun(arg1T,arg2T,arg3T) const" in process owning item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        task(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3, const TaskAttributes &attr = TaskAttributes()) const
        {
            return const_cast<containerT *>(this)->task(key, memfun, arg1, arg2, arg3, attr);
        }

        /// Adds task "resultT memfun(arg1T,arg2T,arg3T, arg4T) const" in process owning item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        task(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3, const arg4T &arg4, const TaskAttributes &attr = TaskAttributes()) const
        {
            return const_cast<containerT *>(this)->task(key, memfun, arg1, arg2, arg3, arg4, attr);
        }

        /// Adds task "resultT memfun(arg1T,arg2T,arg3T,arg4T,arg5T) const" in process owning item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        task(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3, const arg4T &arg4, const arg5T &arg5, const TaskAttributes &attr = TaskAttributes()) const
        {
            return const_cast<containerT *>(this)->task(key, memfun, arg1, arg2, arg3, arg4, arg5, attr);
        }

        /// Adds task "resultT memfun(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T) const" in process owning item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        task(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3, const arg4T &arg4, const arg5T &arg5, const arg6T &arg6, const TaskAttributes &attr = TaskAttributes()) const
        {
            return const_cast<containerT *>(this)->task(key, memfun, arg1, arg2, arg3, arg4, arg5, arg6, attr);
        }

        /// Adds task "resultT memfun(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T) const" in process owning item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        task(const keyT &key, memfunT memfun, const arg1T &arg1, const arg2T &arg2, const arg3T &arg3, const arg4T &arg4, const arg5T &arg5, const arg6T &arg6, const arg7T &arg7, const TaskAttributes &attr = TaskAttributes()) const
        {
            return const_cast<containerT *>(this)->task(key, memfun, arg1, arg2, arg3, arg4, arg5, arg6, arg7, attr);
        }

        /// (de)Serialize --- *Local* data only to/from anything *except* Buffer*Archive and Parallel*Archive

        /// Advisable for *you* to fence before and after this to ensure consistency
        template <typename Archive>
        void serialize(const Archive &ar)
        {
            //
            // !! If you change the format of this stream make sure that
            // !! the parallel in/out archive below is compatible
            //
            const long magic = 5881828; // Sitar Indian restaurant in Knoxville
            unsigned long count = 0;
            check_initialized();

            if (Archive::is_output_archive)
            {
                ar & magic;
                for (iterator it = begin(); it != end(); ++it)
                    count++;
                ar & count;
                for (iterator it = begin(); it != end(); ++it)
                    ar &*it;
            }
            else
            {
                long cookie = 0l;
                ar & cookie;
                MADNESS_ASSERT(cookie == magic);
                ar & count;
                while (count--)
                {
                    pairT datum;
                    ar & datum;
                    replace(datum);
                }
            }
        }

        /// (de)Serialize --- !! ONLY for purpose of interprocess communication

        /// This just writes/reads the unique id to/from the Buffer*Archive.
        void serialize(const archive::BufferOutputArchive &ar)
        {
            check_initialized();
            ar &static_cast<WorldObject<implT> *>(p.get());
        }

        /// (de)Serialize --- !! ONLY for purpose of interprocess communication

        /// This just writes/reads the unique id to/from the Buffer*Archive.
        void serialize(const archive::BufferInputArchive &ar)
        {
            WorldObject<implT> *ptr = nullptr;
            ar & ptr;
            MADNESS_ASSERT(ptr);

#ifdef MADNESS_DISABLE_SHARED_FROM_THIS
            p.reset(static_cast<implT *>(ptr), [](implT *p_) -> void{});
#else
            p = static_cast<implT *>(ptr)->shared_from_this();
#endif // MADNESS_DISABLE_SHARED_FROM_THIS
        }

        /// Returns the associated unique id ... must be initialized
        const uniqueidT &id() const
        {
            check_initialized();
            return p->id();
        }

        /// Destructor passes ownership of implementation to world for deferred cleanup
        virtual ~WorldContainer()
        {
            detail::deferred_cleanup(p->get_world(), p);
        }

        friend void swap<>(WorldContainer &, WorldContainer &);
    };

    /// Swaps the content of two WorldContainer objects. It should be called on all nodes.

    /// \ingroup worlddc
    template <typename keyT, typename valueT, typename hashfunT>
    void swap(WorldContainer<keyT, valueT, hashfunT> &dc0, WorldContainer<keyT, valueT, hashfunT> &dc1)
    {
        std::swap(dc0.p, dc1.p);
    }


    namespace archive
    {

        /// Write container to parallel archive

        /// specialization for parallel serialization of a WorldContainer:
        /// all threads on each process serialize some values into a buffer, which gets concatenated
        /// and finally serialized to localarchive (aka VectorOutputArchive).
        template <class keyT, class valueT>
        struct ArchiveStoreImpl<ParallelOutputArchive<VectorOutputArchive>, WorldContainer<keyT, valueT>>
        {
            static void store(const ParallelOutputArchive<VectorOutputArchive> &ar, const WorldContainer<keyT, valueT> &t)
            {
                using localarchiveT = VectorOutputArchive;
                const long magic = -5881828; // Sitar Indian restaurant in Knoxville (negative to indicate parallel!)
                typedef WorldContainer<keyT, valueT> dcT;
                using const_iterator = typename dcT::const_iterator;
                int count = t.size(); // Must be INT for MPI and NOT const since we'll do a global sum eventually

                // Strategy:
                // 1. Serialize local data to a buffer in parallel over threads
                //    a) Compute the size of the buffer needed by each task
                //    b) Sum sizes and allocate the buffer of exact sizes needed for all threads
                //    c) Serialize the data into the buffer in parallel over threads
                // 2. Gather all buffers to process 0

                World *world = ar.get_world();
                world->gop.fence(); // Global fence here

                class op_inspector : public TaskInterface
                {
                    const_iterator start, end;
                    size_t &size;

                public:
                    op_inspector(const_iterator start, const_iterator end, size_t &size)
                        : start(start), end(end), size(size) {}
                    void run(World &world)
                    {
                        BufferOutputArchive bo;
                        for (const_iterator it = start; it != end; ++it)
                            bo &*it;
                        size = bo.size();
                    }
                };

                class op_executor : public TaskInterface
                {
                    const_iterator start, end;
                    unsigned char *buf;
                    const size_t size;

                public:
                    op_executor(const_iterator start, const_iterator end, unsigned char *buf, size_t size)
                        : start(start), end(end), buf(buf), size(size) {}
                    void run(World &world)
                    {
                        BufferOutputArchive bo(buf, size);
                        for (const_iterator it = start; it != end; ++it)
                        {
                            bo &*it;
                        }
                        MADNESS_CHECK(size == bo.size());
                    }
                };

                // No need for LOCAL fence here since only master thread is busy
                double wall0 = wall_time();
                const size_t ntasks = std::min(size_t(count), std::max(size_t(1), ThreadPool::size()));
                size_t local_size = 0;
                double wall1 = wall0;
                unsigned char* buf = 0;
                if (ntasks > 0)
                {
                    const size_t max_items_per_task = (std::max(1, count) - 1) / ntasks + 1;
                    // Compute the size of the buffer needed by each task
                    std::vector<const_iterator> starts(ntasks), ends(ntasks);
                    std::vector<size_t> local_sizes(ntasks);
                    const_iterator start = t.begin();
                    size_t nleft = count;
                    for (size_t taskid = 0; taskid < ntasks; taskid++)
                    {
                        const_iterator end = start;
                        if (taskid == (ntasks - 1))
                        {
                            end = t.end();
                        }
                        else
                        {
                            size_t nitems = std::min(max_items_per_task, nleft);
                            std::advance(end, max_items_per_task);
                            nleft -= nitems;
                        }
                        starts[taskid] = start;
                        ends[taskid] = end;
                        world->taskq.add(new op_inspector(start, end, local_sizes[taskid])); // Be sure to pass iterators by value!!
                        start = end;
                    }
                    world->taskq.fence(); // just need LOCAL fence
                    wall1 = wall_time();
                    // if (world->rank() == 0)
                        // printf("time in op_inspector: %8.4fs\n", wall1 - wall0);
                    wall0 = wall1;

                    // total size over all threads
                    for (size_t taskid = 0; taskid < ntasks; taskid++)
                    {
                        local_size += local_sizes[taskid];
                        // print("taskid",taskid,"size",local_sizes[taskid]);
                    }

                    // Allocate the buffer for all threads
                    buf = new unsigned char[local_size];

                    // Now execute the serialization
                    size_t offset = 0;
                    for (size_t taskid = 0; taskid < ntasks; taskid++)
                    {
                        world->taskq.add(new op_executor(starts[taskid], ends[taskid], buf + offset, local_sizes[taskid]));
                        offset += local_sizes[taskid];
                    }
                    world->taskq.fence(); // just need LOCAL fence

                    wall1 = wall_time();
                    // if (world->rank() == 0)
                        // printf("time in op_executor: %8.4fs\n", wall1 - wall0);
                    wall0 = wall1;
                }
                // VERify that the serialization worked!!
                // {
                //     BufferInputArchive bi(buf, local_size);
                //     for (int item=0; item<count; item++) {
                //         std::pair<keyT, valueT> datum;
                //         bi & datum;
                //         print("deserializing",datum.first);
                //     }
                // }

                // Gather all buffers to process 0
                // first gather all of the sizes and counts to a vector in process 0
                const int size = local_size;
                std::vector<int> sizes(world->size());
                MPI_Gather(&size, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, world->mpi.comm().Get_mpi_comm());
                world->gop.sum(count); // just need total number of elements

                // print("time 3",wall_time());
                //  build the cumulative sum of sizes
                std::vector<int> offsets(world->size());
                offsets[0] = 0;
                for (int i = 1; i < world->size(); ++i)
                    offsets[i] = offsets[i - 1] + sizes[i - 1];
                size_t total_size = offsets.back() + sizes.back();
                // if (world->rank() == 0)
                    // print("total_size", total_size);

                // print("time 4",wall_time());
                // gather the vector of data v from each process to process 0
                unsigned char *all_data = 0;
                if (world->rank() == 0)
                {
                    all_data = new unsigned char[total_size];
                }
                MPI_Gatherv(buf, local_size, MPI_BYTE, all_data, sizes.data(), offsets.data(), MPI_BYTE, 0, world->mpi.comm().Get_mpi_comm());

                wall1 = wall_time();
                // if (world->rank() == 0)
                    // printf("time in gather+gatherv: %8.4fs\n", wall1 - wall0);
                wall0 = wall1;

                delete[] buf;

                // print("time 5",wall_time());
                if (world->rank() == 0)
                {
                    auto &localar = ar.local_archive();
                    localar & magic & 1; // 1 client
                    // localar & t;
                    ArchivePrePostImpl<localarchiveT, dcT>::preamble_store(localar);
                    localar & -magic &(unsigned long)(count);
                    localar.store(all_data, total_size);
                    ArchivePrePostImpl<localarchiveT, dcT>::postamble_store(localar);
                    wall1 = wall_time();
                    // if (world->rank() == 0)
                        // printf("time in final copy on node 0: %8.4fs\n", wall1 - wall0);

                    delete[] all_data;
                }
                world->gop.fence();
                // print("time 6",wall_time());
            }
        };

        /// Write container to parallel archive with optional fence

        /// \ingroup worlddc
        /// Each node (process) is served by a designated IO node.
        /// The IO node has a binary local file archive to which is
        /// first written a cookie and the number of servers.  The IO
        /// node then loops thru all of its clients and in turn tells
        /// each to write its data over an MPI stream, which is copied
        /// directly to the output file.  The stream contents are then
        /// cookie, no. of clients, foreach client (usual sequential archive).
        ///
        /// If ar.dofence() is true (default) fence is invoked before and
        /// after the IO. The fence is optional but it is of course
        /// necessary to be sure that all updates have completed
        /// before doing IO, and that all IO has completed before
        /// subsequent modifications. Also, there is always at least
        /// some synchronization between a client and its IO server.
        template <class keyT, class valueT, class localarchiveT>
        struct ArchiveStoreImpl<ParallelOutputArchive<localarchiveT>, WorldContainer<keyT, valueT>>
        {
            static void store(const ParallelOutputArchive<localarchiveT> &ar, const WorldContainer<keyT, valueT> &t)
            {
                const long magic = -5881828; // Sitar Indian restaurant in Knoxville (negative to indicate parallel!)
                typedef WorldContainer<keyT, valueT> dcT;
                // typedef typename dcT::const_iterator iterator; // unused?
                typedef typename dcT::pairT pairT;
                World *world = ar.get_world();
                Tag tag = world->mpi.unique_tag();
                ProcessID me = world->rank();
                if (ar.dofence())
                    world->gop.fence();
                if (ar.is_io_node())
                {
                    auto &localar = ar.local_archive();
                    localar & magic & ar.num_io_clients();
                    for (ProcessID p = 0; p < world->size(); ++p)
                    {
                        if (p == me)
                        {
                            localar & t;
                        }
                        else if (ar.io_node(p) == me)
                        {
                            world->mpi.Send(int(1), p, tag); // Tell client to start sending
                            archive::MPIInputArchive source(*world, p);
                            long cookie = 0l;
                            unsigned long count = 0ul;

                            ArchivePrePostImpl<localarchiveT, dcT>::preamble_store(localar);

                            source & cookie & count;
                            localar & cookie & count;
                            while (count--)
                            {
                                pairT datum;
                                source & datum;
                                localar & datum;
                            }

                            ArchivePrePostImpl<localarchiveT, dcT>::postamble_store(localar);
                        }
                    }
                }
                else
                {
                    ProcessID p = ar.my_io_node();
                    int flag;
                    world->mpi.Recv(flag, p, tag);
                    MPIOutputArchive dest(*world, p);
                    dest & t;
                    dest.flush();
                }
                if (ar.dofence())
                    world->gop.fence();
            }
        };

        template <class keyT, class valueT, class localarchiveT>
        struct ArchiveLoadImpl<ParallelInputArchive<localarchiveT>, WorldContainer<keyT, valueT>>
        {
            /// Read container from parallel archive

            /// \ingroup worlddc
            /// See store method above for format of file content.
            /// !!! We presently ASSUME that the number of writers and readers are
            /// the same.  This is frustrating but not a show stopper since you
            /// can always run a separate job to copy to a different number.
            ///
            /// The IO node simply reads all data and inserts entries.
            static void load(const ParallelInputArchive<localarchiveT> &ar, WorldContainer<keyT, valueT> &t)
            {
                const long magic = -5881828; // Sitar Indian restaurant in Knoxville (negative to indicate parallel!)
                // typedef WorldContainer<keyT,valueT> dcT; // unused
                // typedef typename dcT::iterator iterator; // unused
                // typedef typename dcT::pairT pairT; // unused
                World *world = ar.get_world();
                if (ar.dofence())
                    world->gop.fence();
                if (ar.is_io_node())
                {
                    long cookie = 0l;
                    int nclient = 0;
                    auto &localar = ar.local_archive();
                    localar & cookie & nclient;
                    MADNESS_CHECK(cookie == magic);
                    while (nclient--)
                    {
                        localar & t;
                    }
                }
                if (ar.dofence())
                    world->gop.fence();
            }
        };
    }

}

///@}

#endif // MADNESS_WORLD_WORLDDC_H__INCLUDED
