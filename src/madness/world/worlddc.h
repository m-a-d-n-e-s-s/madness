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

#include <madness/world/parallel_archive.h>
#include <madness/world/worldhashmap.h>
#include <madness/world/mpi_archive.h>
#include <madness/world/world_object.h>

namespace madness {

    template <typename keyT, typename valueT, typename hashfunT>
    class WorldContainer;

    template <typename keyT, typename valueT, typename hashfunT>
    class WorldContainerImpl;

    template <typename keyT, typename valueT, typename hashfunT>
    void swap(WorldContainer<keyT, valueT, hashfunT>&, WorldContainer<keyT, valueT, hashfunT>&);

    template <typename keyT>
    class WorldDCPmapInterface;

    template <typename keyT>
    class WorldDCRedistributeInterface {
    public:
        virtual std::size_t size() const = 0;
        virtual void redistribute_phase1(const std::shared_ptr< WorldDCPmapInterface<keyT> >& newmap) = 0;
        virtual void redistribute_phase2() = 0;
        virtual void redistribute_phase3() = 0;
	virtual ~WorldDCRedistributeInterface() {};
    };


    /// Interface to be provided by any process map

    /// \ingroup worlddc
    template <typename keyT>
    class WorldDCPmapInterface {
    public:
        typedef WorldDCRedistributeInterface<keyT>* ptrT;
    private:
        std::set<ptrT> ptrs;
    public:
        /// Maps key to processor

        /// @param[in] key Key for container
        /// @return Processor that logically owns the key
        virtual ProcessID owner(const keyT& key) const = 0;

        virtual ~WorldDCPmapInterface() {}

        virtual void print() const {}

        /// Registers object for receipt of redistribute callbacks

        /// @param[in] ptr Pointer to class derived from WorldDCRedistributedInterface
        void register_callback(ptrT ptr) {
            ptrs.insert(ptr);
        }

        /// Deregisters object for receipt of redistribute callbacks

        /// @param[in] ptr Pointer to class derived from WorldDCRedistributedInterface
        void deregister_callback(ptrT ptr) {
            ptrs.erase(ptr);
        }

        /// Invoking this switches all registered objects from this process map to the new one

        /// After invoking this routine all objects will be registered with the
        /// new map and no objects will be registered in the current map.
        /// @param[in] world The associated world
        /// @param[in] newpmap The new process map
        void redistribute(World& world, const std::shared_ptr< WorldDCPmapInterface<keyT> >& newpmap) {
            print_data_sizes(world, "before redistributing");
            world.gop.fence();
            for (typename std::set<ptrT>::iterator iter = ptrs.begin();
                 iter != ptrs.end();
                 ++iter) {
                (*iter)->redistribute_phase1(newpmap);
            }
            world.gop.fence();
            for (typename std::set<ptrT>::iterator iter = ptrs.begin();
                 iter != ptrs.end();
                 ++iter) {
                (*iter)->redistribute_phase2();
                newpmap->register_callback(*iter);
            }
            world.gop.fence();
            for (typename std::set<ptrT>::iterator iter = ptrs.begin();
                 iter != ptrs.end();
                 ++iter) {
	         (*iter)->redistribute_phase3();
            }
            world.gop.fence();
            ptrs.clear();
            newpmap->print_data_sizes(world, "after redistributing");
        }

        /// Counts global number of entries in all containers associated with this process map

        /// Collective operation with global fence
        std::size_t global_size(World& world) const {
            world.gop.fence();
            std::size_t sum = local_size();
            world.gop.sum(sum);
            world.gop.fence();
            return sum;
        }

        /// Counts local number of entries in all containers associated with this process map
        std::size_t local_size() const {
            std::size_t sum = 0;
            for (typename std::set<ptrT>::iterator iter = ptrs.begin(); iter != ptrs.end(); ++iter) {
                sum += (*iter)->size();
            }
            return sum;
        }

        /// Prints size info to std::cout

        /// Collective operation with global fence
        void print_data_sizes(World& world, const std::string msg="") const {
            world.gop.fence();
            std::size_t total = global_size(world);
            std::vector<std::size_t> sizes(world.size());
            sizes[world.rank()] = local_size();
            world.gop.sum(&sizes[0],world.size());
            if (world.rank() == 0) {
                madness::print("data distribution info", msg);
                madness::print("   total: ", total);
                std::cout << "   procs: ";
                for (int i=0; i<world.size(); i++) std::cout << sizes[i] << " ";
                std::cout << std::endl;
            }
            world.gop.fence();
        }
    };

    /// Default process map is "random" using madness::hash(key)

    /// \ingroup worlddc
    template <typename keyT, typename hashfunT = Hash<keyT> >
    class WorldDCDefaultPmap : public WorldDCPmapInterface<keyT> {
    private:
        const int nproc;
        hashfunT hashfun;
    public:
        WorldDCDefaultPmap(World& world, const hashfunT& hf = hashfunT()) :
            nproc(world.mpi.nproc()),
            hashfun(hf)
        { }

        ProcessID owner(const keyT& key) const {
            if (nproc == 1) return 0;
            return hashfun(key)%nproc;
        }
    };

    /// Local process map will always return the current process as owner

    /// \ingroup worlddc
    template <typename keyT, typename hashfunT = Hash<keyT> >
    class WorldDCLocalPmap : public WorldDCPmapInterface<keyT> {
    private:
    	ProcessID me;
    public:
    	WorldDCLocalPmap(World& world) : me(world.rank())  { }
    	ProcessID owner(const keyT& key) const {
    		return me;
    	}
    };

    /// Iterator for distributed container wraps the local iterator

    /// \ingroup worlddc
    template <class internal_iteratorT>
    class WorldContainerIterator {
    public:
      typedef typename std::iterator_traits<internal_iteratorT>::iterator_category iterator_category;
      typedef typename std::iterator_traits<internal_iteratorT>::value_type value_type;
      typedef typename std::iterator_traits<internal_iteratorT>::difference_type difference_type;
      typedef typename std::iterator_traits<internal_iteratorT>::pointer pointer;
      typedef typename std::iterator_traits<internal_iteratorT>::reference reference;

    private:
        internal_iteratorT  it;       ///< Iterator from local container
        // TODO: Convert this to a scoped pointer.
        mutable value_type* value;    ///< holds the remote values

    public:
        /// Default constructor makes a local uninitialized value
        explicit WorldContainerIterator()
                : it(), value(nullptr) {}

        /// Initializes from a local iterator
        explicit WorldContainerIterator(const internal_iteratorT& it)
                : it(it), value(nullptr) {}

        /// Initializes to cache a remote value
        explicit WorldContainerIterator(const value_type& v)
                : it(), value(nullptr)
        {
            value = new value_type(v);
        }

        WorldContainerIterator(const WorldContainerIterator& other)
                : it(), value(nullptr)
        {
            copy(other);
        }

        template <class iteratorT>
        WorldContainerIterator(const WorldContainerIterator<iteratorT>& other)
                : it(), value(nullptr)
        {
            copy(other);
        }

        ~WorldContainerIterator() {
            delete value;
        }

        /// Assignment
        WorldContainerIterator& operator=(const WorldContainerIterator& other) {
            copy(other);
            return *this;
        }

        /// Determines if two iterators are identical
        bool operator==(const WorldContainerIterator& other) const {
            return (((!is_cached()) && (!other.is_cached())) && it == other.it) ||
                ((is_cached() && other.is_cached()) && value->first == other.value->first);
        }


        /// Determines if two iterators are different
        bool operator!=(const WorldContainerIterator& other) const {
            return !(*this == other);
        }


        /// Pre-increment of an iterator (i.e., ++it) --- \em local iterators only

        /// Trying to increment a remote iterator will throw
        WorldContainerIterator& operator++() {
            MADNESS_ASSERT( !is_cached() );
            ++it;
            return *this;
        }

        WorldContainerIterator operator++(int) {
            MADNESS_ASSERT( !is_cached() );
            WorldContainerIterator<internal_iteratorT> result(*this);
            ++it;
            return result;
        }

        /// Iterators dereference to std::pair<const keyT,valueT>
        pointer operator->() const {
            return (is_cached() ? value : it.operator->() );
        }

        /// Iterators dereference to std::pair<const keyT,valueT>
        reference operator*() const {
            return (is_cached() ? *value : *it );
        }

        /// Private: (or should be) Returns iterator of internal container
        const internal_iteratorT& get_internal_iterator() const {
            return it;
        }

        /// Returns true if this is non-local or cached value
        bool is_cached() const {
            return value != nullptr;
        }

        template <typename Archive>
        void serialize(const Archive&) {
            MADNESS_EXCEPTION("Serializing DC iterator ... why?", false);
        }

    private:
        template <class iteratorT>
        friend class WorldContainerIterator;

        template <class iteratorT>
        void copy(const WorldContainerIterator<iteratorT>& other) {
            if (static_cast<const void*>(this) != static_cast<const void*>(&other)) {
                delete value;
                if(other.is_cached()) {
                    value = new value_type(* other.value);
                    it = internal_iteratorT();
                } else {
                    it = other.it;
                    value = nullptr;
                }
            }
        }
    };

    /// Internal implementation of distributed container to facilitate shallow copy

    /// \ingroup worlddc
    template <typename keyT, typename valueT, typename hashfunT >
    class WorldContainerImpl
        : public WorldObject< WorldContainerImpl<keyT, valueT, hashfunT> >
        , public WorldDCRedistributeInterface<keyT>
#ifndef MADNESS_DISABLE_SHARED_FROM_THIS
        , public std::enable_shared_from_this<WorldContainerImpl<keyT, valueT, hashfunT> >
#endif // MADNESS_DISABLE_SHARED_FROM_THIS
    {
    public:
        typedef typename std::pair<const keyT,valueT> pairT;
        typedef const pairT const_pairT;
        typedef WorldContainerImpl<keyT,valueT,hashfunT> implT;

        typedef ConcurrentHashMap< keyT,valueT,hashfunT > internal_containerT;

	//typedef WorldObject< WorldContainerImpl<keyT, valueT, hashfunT> > worldobjT;

        typedef typename internal_containerT::iterator internal_iteratorT;
        typedef typename internal_containerT::const_iterator internal_const_iteratorT;
        typedef typename internal_containerT::accessor accessor;
        typedef typename internal_containerT::const_accessor const_accessor;
        typedef WorldContainerIterator<internal_iteratorT> iteratorT;
        typedef WorldContainerIterator<internal_iteratorT> iterator;
        typedef WorldContainerIterator<internal_const_iteratorT> const_iteratorT;
        typedef WorldContainerIterator<internal_const_iteratorT> const_iterator;

        friend class WorldContainer<keyT,valueT,hashfunT>;

//         template <typename containerT, typename datumT>
//         inline
//         static
//         typename containerT::iterator replace(containerT& c, const datumT& d) {
//             std::pair<typename containerT::iterator,bool> p = c.insert(d);
//             if (!p.second) p.first->second = d.second;   // Who's on first?
//             return p.first;
//         }

    private:

        WorldContainerImpl();   // Inhibit default constructor

        std::shared_ptr< WorldDCPmapInterface<keyT> > pmap;///< Function/class to map from keys to owning process
        const ProcessID me;                      ///< My MPI rank
        internal_containerT local;               ///< Locally owned data
        std::vector<keyT>* move_list;            ///< Tempoary used to record data that needs redistributing

        /// Handles find request
        void find_handler(ProcessID requestor, const keyT& key, const RemoteReference< FutureImpl<iterator> >& ref) {
            internal_iteratorT r = local.find(key);
            if (r == local.end()) {
                //print("find_handler: failure:", key);
                this->send(requestor, &implT::find_failure_handler, ref);
            }
            else {
                //print("find_handler: success:", key, r->first, r->second);
                this->send(requestor, &implT::find_success_handler, ref, *r);
            }
        }

        /// Handles successful find response
        void find_success_handler(const RemoteReference< FutureImpl<iterator> >& ref, const pairT& datum) {
            FutureImpl<iterator>* f = ref.get();
            f->set(iterator(datum));
            //print("find_success_handler: success:", datum.first, datum.second, f->get()->first, f->get()->second);
            // Todo: Look at this again.
//            ref.reset(); // Matching inc() in find() where ref was made
        }

        /// Handles unsuccessful find response
        void find_failure_handler(const RemoteReference< FutureImpl<iterator> >& ref) {
            FutureImpl<iterator>* f = ref.get();
            f->set(end());
            //print("find_failure_handler");
            // Todo: Look at this again.
//            ref.reset(); // Matching inc() in find() where ref was made
        }

    public:

        WorldContainerImpl(World& world,
                           const std::shared_ptr< WorldDCPmapInterface<keyT> >& pm,
                           const hashfunT& hf)
                : WorldObject< WorldContainerImpl<keyT, valueT, hashfunT> >(world)
                , pmap(pm)
                , me(world.mpi.rank())
                , local(5011, hf) {
            pmap->register_callback(this);
        }

        virtual ~WorldContainerImpl() {
            pmap->deregister_callback(this);
        }

        const std::shared_ptr< WorldDCPmapInterface<keyT> >& get_pmap() const {
            return pmap;
        }

        std::shared_ptr< WorldDCPmapInterface<keyT> >& get_pmap() {
            return pmap;
        }

        void reset_pmap_to_local() {
            pmap->deregister_callback(this);
            pmap.reset(new WorldDCLocalPmap<keyT>(this->get_world()));
            pmap->register_callback(this);
        }

        /// replicates this WorldContainer on all ProcessIDs and generates a
        /// ProcessMap where all nodes are local
        void replicate(bool fence) {

        	World& world=this->get_world();
        	pmap->deregister_callback(this);
        	pmap.reset(new WorldDCLocalPmap<keyT>(world));
        	pmap->register_callback(this);

        	for (ProcessID rank=0; rank<world.size(); rank++) {
        		if (rank == world.rank()) {
        			std::size_t sz = size();
        			world.gop.broadcast_serializable(sz, rank);

        			for (auto it=begin(); it!=end(); ++it) {
        				keyT key = it->first;
        				valueT value = it->second;
        				world.gop.broadcast_serializable(key, rank);
        				world.gop.broadcast_serializable(value, rank);
        			}
        		}
        		else {
        			size_t sz;
        			world.gop.broadcast_serializable(sz, rank);
        			for (size_t i=0; i<sz; i++) {
        				keyT key;
        				valueT value;
        				world.gop.broadcast_serializable(key, rank);
        				world.gop.broadcast_serializable(value, rank);
        				insert(pairT(key,value));
        			}
        		}
        	}
        	if (fence) world.gop.fence();
        }

        hashfunT& get_hash() const { return local.get_hash(); }

        bool is_local(const keyT& key) const {
            return owner(key) == me;
        }

        ProcessID owner(const keyT& key) const {
            return pmap->owner(key);
        }

        bool probe(const keyT& key) const {
            ProcessID dest = owner(key);
            if (dest == me)
                return local.find(key) != local.end();
            else
                return false;
        }

        std::size_t size() const {
            return local.size();
        }

        void insert(const pairT& datum) {
            ProcessID dest = owner(datum.first);
            if (dest == me) {
                // Was using iterator ... try accessor ?????
                accessor acc;
                // N.B. key might already exist if want to simply replace
                [[maybe_unused]] auto inserted = local.insert(acc,datum.first);
                acc->second = datum.second;
            }
            else {
  	        // Must be send (not task) for sequential consistency (and relies on single-threaded remote server)
                this->send(dest, &implT::insert, datum);
            }
        }

        bool insert_acc(accessor& acc, const keyT& key) {
            MADNESS_ASSERT(owner(key) == me);
            return local.insert(acc,key);
        }

        bool insert_const_acc(const_accessor& acc, const keyT& key) {
            MADNESS_ASSERT(owner(key) == me);
            return local.insert(acc,key);
        }

        void clear() {
            local.clear();
        }


        void erase(const keyT& key) {
            ProcessID dest = owner(key);
            if (dest == me) {
                [[maybe_unused]] auto erased = local.try_erase(key);
                MADNESS_ASSERT(erased);
            }
            else {
                void(implT::*eraser)(const keyT&) = &implT::erase;
                this->send(dest, eraser, key);
            }
        }

        template <typename InIter>
        void erase(InIter it) {
            MADNESS_ASSERT(!it.is_cached());
            MADNESS_ASSERT(it != end());
            erase(it->first);
        }

        template <typename InIter>
        void erase(InIter first, InIter last) {
            InIter it = first;
            do {
                first++;
                erase(it->first);
                it = first;
            } while(first != last);
        }

        iterator begin() {
            return iterator(local.begin());
        }

        const_iterator begin() const {
            return const_iterator(local.begin());
        }

        iterator end() {
            return iterator(local.end());
        }

        const_iterator end() const {
            return const_iterator(local.end());
        }

        Future<const_iterator> find(const keyT& key) const {
            // Ugliness here to avoid replicating find() and
            // associated handlers for const.  Assumption is that
            // const and non-const iterators are identical except for
            // const attribute ... at some point probably need to do
            // the right thing.
            Future<iterator> r = const_cast<implT*>(this)->find(key);
            return *(Future<const_iterator>*)(&r);
        }


        Future<iterator> find(const keyT& key) {
            ProcessID dest = owner(key);
            if (dest == me) {
                return Future<iterator>(iterator(local.find(key)));
            } else {
                Future<iterator> result;
                this->send(dest, &implT::find_handler, me, key, result.remote_ref(this->get_world()));
                return result;
            }
        }

        bool find(accessor& acc, const keyT& key) {
            if (owner(key) != me) return false;
            return local.find(acc,key);
        }


        bool find(const_accessor& acc, const keyT& key) const {
            if (owner(key) != me) return false;
            return local.find(acc,key);
        }


        // Used to forward call to item member function
        template <typename memfunT>
        MEMFUN_RETURNT(memfunT)
        itemfun(const keyT& key, memfunT memfun) {
            accessor acc;
            // N.B. key may already exist, this is just to ensure lock is held by acc
            [[maybe_unused]] auto inserted = local.insert(acc, key);
            return (acc->second.*memfun)();
        }

        // Used to forward call to item member function
        template <typename memfunT, typename arg1T>
        MEMFUN_RETURNT(memfunT)
        itemfun(const keyT& key, memfunT memfun, const arg1T& arg1) {
            accessor acc;
            // N.B. key may already exist, this is just to ensure lock is held by acc
            [[maybe_unused]] auto inserted = local.insert(acc, key);
            return (acc->second.*memfun)(arg1);
        }

        // Used to forward call to item member function
        template <typename memfunT, typename arg1T, typename arg2T>
        MEMFUN_RETURNT(memfunT)
        itemfun(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2) {
            accessor acc;
            // N.B. key may already exist, this is just to ensure lock is held by acc
            [[maybe_unused]] auto inserted = local.insert(acc, key);
            return (acc->second.*memfun)(arg1,arg2);
        }

        // Used to forward call to item member function
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        MEMFUN_RETURNT(memfunT)
        itemfun(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3) {
            accessor acc;
            // N.B. key may already exist, this is just to ensure lock is held by acc
            [[maybe_unused]] auto inserted = local.insert(acc, key);
            return (acc->second.*memfun)(arg1,arg2,arg3);
        }

        // Used to forward call to item member function
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        MEMFUN_RETURNT(memfunT)
        itemfun(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4) {
            accessor acc;
            // N.B. key may already exist, this is just to ensure lock is held by acc
            [[maybe_unused]] auto inserted = local.insert(acc, key);
            return (acc->second.*memfun)(arg1,arg2,arg3,arg4);
        }

        // Used to forward call to item member function
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        MEMFUN_RETURNT(memfunT)
        itemfun(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5) {
            accessor acc;
            // N.B. key may already exist, this is just to ensure lock is held by acc
            [[maybe_unused]] auto inserted = local.insert(acc, key);
            return (acc->second.*memfun)(arg1,arg2,arg3,arg4,arg5);
        }

        // Used to forward call to item member function
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T>
        MEMFUN_RETURNT(memfunT)
        itemfun(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5, const arg6T& arg6) {
            accessor acc;
            // N.B. key may already exist, this is just to ensure lock is held by acc
            [[maybe_unused]] auto inserted = local.insert(acc, key);
            return (acc->second.*memfun)(arg1,arg2,arg3,arg4,arg5,arg6);
        }

        // Used to forward call to item member function
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T>
        MEMFUN_RETURNT(memfunT)
        itemfun(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3,
				const arg4T& arg4, const arg5T& arg5, const arg6T& arg6, const arg7T& arg7) {
            accessor acc;
            // N.B. key may already exist, this is just to ensure lock is held by acc
            [[maybe_unused]] auto inserted = local.insert(acc, key);
            return (acc->second.*memfun)(arg1,arg2,arg3,arg4,arg5,arg6,arg7);
        }

        // First phase of redistributions changes pmap and makes list of stuff to move
        void redistribute_phase1(const std::shared_ptr< WorldDCPmapInterface<keyT> >& newpmap) {
            pmap = newpmap;
            move_list = new std::vector<keyT>();
            for (typename internal_containerT::iterator iter=local.begin(); iter!=local.end(); ++iter) {
                if (owner(iter->first) != me) move_list->push_back(iter->first);
            }
        }

	struct P2Op {
	  implT * impl;
	  typedef Range<typename std::vector<keyT>::const_iterator> rangeT;
	  P2Op(implT* impl) : impl(impl) {}
    	  P2Op(const P2Op& p) : impl(p.impl) {}
	  bool operator()(typename rangeT::iterator& iterator) const {
	    typename internal_containerT::iterator iter = impl->local.find(*iterator);
	    MADNESS_ASSERT(iter != impl->local.end());

	    //impl->insert(*iter);
	    impl->task(impl->owner(*iterator), &implT::insert, *iter);

	    impl->local.erase(iter); // delete local copy of the data
	    return true;
	  }
	};

        // Second phase moves data
        void redistribute_phase2() {
	  this->get_world().taskq.for_each(typename P2Op::rangeT(move_list->begin(), move_list->end()), P2Op(this));
	    //std::vector<keyT>& mvlist = *move_list;
            //for (unsigned int i=0; i<move_list->size(); ++i) {
            //    typename internal_containerT::iterator iter = local.find(mvlist[i]);
            //    MADNESS_ASSERT(iter != local.end());
            //    insert(*iter);
            //    local.erase(iter);
            //}
            //delete move_list;
        }

        // Third phase cleans up
        void redistribute_phase3() {
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
    template <typename keyT, typename valueT, typename hashfunT = Hash<keyT> >
    class WorldContainer : public archive::ParallelSerializableObject {
    public:
        typedef WorldContainer<keyT,valueT,hashfunT> containerT;
        typedef WorldContainerImpl<keyT,valueT,hashfunT> implT;
        typedef typename implT::pairT pairT;
        typedef typename implT::iterator iterator;
        typedef typename implT::const_iterator const_iterator;
        typedef typename implT::accessor accessor;
        typedef typename implT::const_accessor const_accessor;
        typedef Future<iterator> futureT;
        typedef Future<const_iterator> const_futureT;

    private:
        std::shared_ptr<implT> p;

        inline void check_initialized() const {
            MADNESS_ASSERT(p);
        }
    public:

        /// Makes an uninitialized container (no communication)

        /// The container is useless until assigned to from a fully
        /// constructed container.  There is no need to worry about
        /// default constructors being executed in order.
        WorldContainer()
                : p()
        {}


        /// Makes an initialized, empty container with default data distribution (no communication)

        /// A unique ID is associated with every distributed container
        /// within a world.  In order to avoid synchronization when
        /// making a container, we have to assume that all processes
        /// execute this constructor in the same order (does not apply
        /// to the non-initializing, default constructor).
        WorldContainer(World& world, bool do_pending=true, const hashfunT& hf = hashfunT())
            : p(new implT(world,
                          std::shared_ptr< WorldDCPmapInterface<keyT> >(new WorldDCDefaultPmap<keyT, hashfunT>(world, hf)),
                          hf))
        {
            if(do_pending)
                p->process_pending();
        }

        /// Makes an initialized, empty container (no communication)

        /// A unique ID is associated with every distributed container
        /// within a world.  In order to avoid synchronization when
        /// making a container, we have to assume that all processes
        /// execute this constructor in the same order (does not apply
        /// to the non-initializing, default constructor).
        WorldContainer(World& world,
                       const std::shared_ptr< WorldDCPmapInterface<keyT> >& pmap,
                       bool do_pending=true,
                       const hashfunT& hf = hashfunT())
            : p(new implT(world, pmap, hf))
        {
            if(do_pending)
                p->process_pending();
        }


        /// Copy constructor is shallow (no communication)

        /// The copy refers to exactly the same container as other
        /// which must be initialized.
        WorldContainer(const WorldContainer& other)
            : p(other.p)
        {
            check_initialized();
        }

        /// Assignment is shallow (no communication)

        /// The copy refers to exactly the same container as other
        /// which must be initialized.
        containerT& operator=(const containerT& other) {
            if (this != &other) {
                other.check_initialized();
                p = other.p;
            }
            return *this;
        }

        /// Returns the world associated with this container
        World& get_world() const {
            check_initialized();
            return p->get_world();
        }

        std::shared_ptr< WorldDCPmapInterface<keyT> >& get_impl() {
            check_initialized();
            return p;
        }

        /// replicates this WorldContainer on all ProcessIDs
        void replicate(bool fence=true) {
        	p->replicate(fence);
        }

        /// Inserts/replaces key+value pair (non-blocking communication if key not local)
        void replace(const pairT& datum) {
            check_initialized();
            p->insert(datum);
        }


        /// Inserts/replaces key+value pair (non-blocking communication if key not local)
        void replace(const keyT& key, const valueT& value) {
            replace(pairT(key,value));
        }


        /// Write access to LOCAL value by key. Returns true if found, false otherwise (always false for remote).
        bool find(accessor& acc, const keyT& key) {
            check_initialized();
            return p->find(acc,key);
        }


        /// Read access to LOCAL value by key. Returns true if found, false otherwise (always false for remote).
        bool find(const_accessor& acc, const keyT& key) const {
            check_initialized();
            return p->find(acc,key);
        }


        /// Write access to LOCAL value by key. Returns true if inserted, false if already exists (throws if remote)
        bool insert(accessor& acc, const keyT& key) {
            check_initialized();
            return p->insert_acc(acc,key);
        }


        /// Read access to LOCAL value by key. Returns true if inserted, false if already exists (throws if remote)
        bool insert(const_accessor& acc, const keyT& key) {
            check_initialized();
            return p->insert_acc(acc,key);
        }


        /// Inserts pairs (non-blocking communication if key(s) not local)
        template <typename input_iterator>
        void replace(input_iterator& start, input_iterator& end) {
            check_initialized();
            using std::placeholders::_1;
            std::for_each(start,end,std::bind(this,std::mem_fn(&containerT::insert),_1));
        }


        /// Returns true if local data is immediately available (no communication)
        bool probe(const keyT& key) const {
            check_initialized();
            return p->probe(key);
        }


        /// Returns processor that logically owns key (no communication)

        /// Local remapping may have changed its physical location, but all
        /// operations should forward correctly.
        inline ProcessID owner(const keyT& key) const {
            check_initialized();
            return p->owner(key);
        }


        /// Returns true if the key maps to the local processor (no communication)
        bool is_local(const keyT& key) const {
            check_initialized();
            return p->is_local(key);
        }


        /// Returns a future iterator (non-blocking communication if key not local)

        /// Like an std::map an iterator "points" to an std::pair<const keyT,valueT>.
        ///
        /// Refer to Future for info on how to avoid blocking.
        Future<iterator> find(const keyT& key) {          //
            check_initialized();
            return p->find(key);
        }


        /// Returns a future iterator (non-blocking communication if key not local)

        /// Like an std::map an iterator "points" to an std::pair<const keyT,valueT>.
        ///
        /// Refer to Future for info on how to avoid blocking.
        Future<const_iterator> find(const keyT& key) const {
            check_initialized();
            return const_cast<const implT*>(p.get())->find(key);
        }


        /// Returns an iterator to the beginning of the \em local data (no communication)
        iterator begin() {
            check_initialized();
            return p->begin();
        }


        /// Returns an iterator to the beginning of the \em local data (no communication)
        const_iterator begin() const {
            check_initialized();
            return const_cast<const implT*>(p.get())->begin();
        }

        /// Returns an iterator past the end of the \em local data (no communication)
        iterator end() {
            check_initialized();
            return p->end();
        }

        /// Returns an iterator past the end of the \em local data (no communication)
        const_iterator end() const {
            check_initialized();
            return const_cast<const implT*>(p.get())->end();
        }

        /// Erases entry from container (non-blocking comm if remote)

        /// Missing keys are quietly ignored.
        ///
        /// Note that erasing an entry may invalidate iterators on the
        /// remote end.  This is just the same as what happens when
        /// using STL iterators on an STL container in a sequential
        /// algorithm.
        void erase(const keyT& key) {
            check_initialized();
            p->erase(key);
        }

        /// Erases entry corresponding to \em local iterator (no communication)
        void erase(const iterator& it) {
            check_initialized();
            p->erase(it);
        }

        /// Erases range defined by \em local iterators (no communication)
        void erase(const iterator& start, const iterator& finish) {
            check_initialized();
            p->erase(start,finish);
        }


        /// Clears all \em local data (no communication)

        /// Invalidates all iterators
        void clear() {
            check_initialized();
            p->clear();
        }

        /// Returns the number of \em local entries (no communication)
        std::size_t size() const {
            check_initialized();
            return p->size();
        }

        /// Returns shared pointer to the process mapping
        inline const std::shared_ptr< WorldDCPmapInterface<keyT> >& get_pmap() const {
            check_initialized();
            return p->get_pmap();
        }

        /// Returns shared pointer to the process mapping
        inline void reset_pmap_to_local() {
            p->reset_pmap_to_local();
        }

        /// Returns a reference to the hashing functor
        hashfunT& get_hash() const {
            check_initialized();
            return p->get_hash();
        }

        /// Process pending messages

        /// If the constructor was given \c do_pending=false then you
        /// \em must invoke this routine in order to process both
        /// prior and future messages.
        inline void process_pending() {
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
        Future< MEMFUN_RETURNT(memfunT) >
        send(const keyT& key, memfunT memfun) {
            check_initialized();
            MEMFUN_RETURNT(memfunT)(implT::*itemfun)(const keyT&, memfunT) = &implT:: template itemfun<memfunT>;
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
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(const keyT& key, const memfunT& memfun, const arg1T& arg1) {
            check_initialized();
            // To work around bug in g++ 4.3.* use static cast as alternative mechanism to force type deduction
            MEMFUN_RETURNT(memfunT) (implT::*itemfun)(const keyT&, memfunT, const arg1T&) = &implT:: template itemfun<memfunT,arg1T>;
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
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2) {
            check_initialized();
            // To work around bug in g++ 4.3.* use static cast as alternative mechanism to force type deduction
            MEMFUN_RETURNT(memfunT) (implT::*itemfun)(const keyT&, memfunT, const arg1T&, const arg2T&) = &implT:: template itemfun<memfunT,arg1T,arg2T>;
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
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3) {
            check_initialized();
            MEMFUN_RETURNT(memfunT)(implT::*itemfun)(const keyT&, memfunT, const arg1T&, const arg2T&, const arg3T&) = &implT:: template itemfun<memfunT,arg1T,arg2T,arg3T>;
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
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4) {
            check_initialized();
            MEMFUN_RETURNT(memfunT)(implT::*itemfun)(const keyT&, memfunT, const arg1T&, const arg2T&, const arg3T&, const arg4T&) = &implT:: template itemfun<memfunT,arg1T,arg2T,arg3T,arg4T>;
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
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5) {
            check_initialized();
            MEMFUN_RETURNT(memfunT)(implT::*itemfun)(const keyT&, memfunT, const arg1T&, const arg2T&, const arg3T&, const arg4T&, const arg5T&) = &implT:: template itemfun<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T>;
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
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5, const arg6T& arg6) {
            check_initialized();
            MEMFUN_RETURNT(memfunT)(implT::*itemfun)(const keyT&, memfunT, const arg1T&, const arg2T&, const arg3T&, const arg4T&, const arg5T&, const arg6T&) = &implT:: template itemfun<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T,arg6T>;
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
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4,
		     const arg5T& arg5, const arg6T& arg6, const arg7T& arg7) {
            check_initialized();
            MEMFUN_RETURNT(memfunT)(implT::*itemfun)(const keyT&, memfunT, const arg1T&, const arg2T&, const arg3T&, const arg4T&, const arg5T&, const arg6T&, const arg7T&) = &implT:: template itemfun<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T>;
            return p->send(owner(key), itemfun, key, memfun, arg1, arg2, arg3, arg4, arg5, arg6, arg7);
        }


        /// Sends message "resultT memfun() const" to item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(const keyT& key, memfunT memfun) const {
            return const_cast<containerT*>(this)->send(key,memfun);
        }

        /// Sends message "resultT memfun(arg1T) const" to item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(const keyT& key, memfunT memfun, const arg1T& arg1) const {
            return const_cast<containerT*>(this)->send(key,memfun,arg1);
        }

        /// Sends message "resultT memfun(arg1T,arg2T) const" to item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2) const {
            return const_cast<containerT*>(this)->send(key,memfun,arg1,arg2);
        }


        /// Sends message "resultT memfun(arg1T,arg2T,arg3T) const" to item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3) const {
            return const_cast<containerT*>(this)->send(key,memfun,arg1,arg2,arg3);
        }

        /// Sends message "resultT memfun(arg1T,arg2T,arg3T,arg4T) const" to item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4) const {
            return const_cast<containerT*>(this)->send(key,memfun,arg1,arg2,arg3,arg4);
        }

        /// Sends message "resultT memfun(arg1T,arg2T,arg3T,arg4T,arg5T) const" to item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5) const {
            return const_cast<containerT*>(this)->send(key,memfun,arg1,arg2,arg3,arg4,arg5);
        }

        /// Sends message "resultT memfun(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T) const" to item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3,
			 const arg4T& arg4, const arg5T& arg5, const arg6T& arg6) const {
            return const_cast<containerT*>(this)->send(key,memfun,arg1,arg2,arg3,arg4,arg5,arg6);
        }

        /// Sends message "resultT memfun(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T) const" to item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3,
			 const arg4T& arg4, const arg5T& arg5, const arg6T& arg6, const arg7T& arg7) const {
            return const_cast<containerT*>(this)->send(key,memfun,arg1,arg2,arg3,arg4,arg5,arg6,arg7);
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
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(const keyT& key, memfunT memfun, const TaskAttributes& attr = TaskAttributes()) {
            check_initialized();
            MEMFUN_RETURNT(memfunT)(implT::*itemfun)(const keyT&, memfunT) = &implT:: template itemfun<memfunT>;
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
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(const keyT& key, memfunT memfun, const arg1T& arg1, const TaskAttributes& attr = TaskAttributes()) {
            check_initialized();
            typedef REMFUTURE(arg1T) a1T;
            MEMFUN_RETURNT(memfunT)(implT::*itemfun)(const keyT&, memfunT, const a1T&) = &implT:: template itemfun<memfunT,a1T>;
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
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const TaskAttributes& attr = TaskAttributes()) {
            check_initialized();
            typedef REMFUTURE(arg1T) a1T;
            typedef REMFUTURE(arg2T) a2T;
            MEMFUN_RETURNT(memfunT)(implT::*itemfun)(const keyT&, memfunT, const a1T&, const a2T&) = &implT:: template itemfun<memfunT,a1T,a2T>;
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
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const TaskAttributes& attr = TaskAttributes()) {
            check_initialized();
            typedef REMFUTURE(arg1T) a1T;
            typedef REMFUTURE(arg2T) a2T;
            typedef REMFUTURE(arg3T) a3T;
            MEMFUN_RETURNT(memfunT)(implT::*itemfun)(const keyT&, memfunT, const a1T&, const a2T&, const a3T&) = &implT:: template itemfun<memfunT,a1T,a2T,a3T>;
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
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const TaskAttributes& attr = TaskAttributes()) {
            check_initialized();
            typedef REMFUTURE(arg1T) a1T;
            typedef REMFUTURE(arg2T) a2T;
            typedef REMFUTURE(arg3T) a3T;
            typedef REMFUTURE(arg4T) a4T;
            MEMFUN_RETURNT(memfunT)(implT::*itemfun)(const keyT&, memfunT, const a1T&, const a2T&, const a3T&, const a4T&) = &implT:: template itemfun<memfunT,a1T,a2T,a3T,a4T>;
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
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5, const TaskAttributes& attr = TaskAttributes()) {
            check_initialized();
            typedef REMFUTURE(arg1T) a1T;
            typedef REMFUTURE(arg2T) a2T;
            typedef REMFUTURE(arg3T) a3T;
            typedef REMFUTURE(arg4T) a4T;
            typedef REMFUTURE(arg5T) a5T;
            MEMFUN_RETURNT(memfunT)(implT::*itemfun)(const keyT&, memfunT, const a1T&, const a2T&, const a3T&, const a4T&, const a5T&) = &implT:: template itemfun<memfunT,a1T,a2T,a3T,a4T,a5T>;
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
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5, const arg6T& arg6, const TaskAttributes& attr = TaskAttributes()) {
            check_initialized();
            typedef REMFUTURE(arg1T) a1T;
            typedef REMFUTURE(arg2T) a2T;
            typedef REMFUTURE(arg3T) a3T;
            typedef REMFUTURE(arg4T) a4T;
            typedef REMFUTURE(arg5T) a5T;
            typedef REMFUTURE(arg6T) a6T;
            MEMFUN_RETURNT(memfunT)(implT::*itemfun)(const keyT&, memfunT, const a1T&, const a2T&, const a3T&, const a4T&, const a5T&, const a6T&) = &implT:: template itemfun<memfunT,a1T,a2T,a3T,a4T,a5T,a6T>;
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
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5, const arg6T& arg6, const arg7T& arg7, const TaskAttributes& attr = TaskAttributes()) {
            check_initialized();
            typedef REMFUTURE(arg1T) a1T;
            typedef REMFUTURE(arg2T) a2T;
            typedef REMFUTURE(arg3T) a3T;
            typedef REMFUTURE(arg4T) a4T;
            typedef REMFUTURE(arg5T) a5T;
            typedef REMFUTURE(arg6T) a6T;
            typedef REMFUTURE(arg7T) a7T;
            MEMFUN_RETURNT(memfunT)(implT::*itemfun)(const keyT&, memfunT, const a1T&, const a2T&, const a3T&, const a4T&, const a5T&, const a6T&, const a7T&) = &implT:: template itemfun<memfunT,a1T,a2T,a3T,a4T,a5T,a6T,a7T>;
            return p->task(owner(key), itemfun, key, memfun, arg1, arg2, arg3, arg4, arg5, arg6, arg7, attr);
        }

        /// Adds task "resultT memfun() const" in process owning item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(const keyT& key, memfunT memfun,  const TaskAttributes& attr = TaskAttributes()) const {
            return const_cast<containerT*>(this)->task(key,memfun,attr);
        }

        /// Adds task "resultT memfun(arg1T) const" in process owning item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(const keyT& key, memfunT memfun, const arg1T& arg1,  const TaskAttributes& attr = TaskAttributes()) const {
            return const_cast<containerT*>(this)->task(key,memfun,arg1,attr);
        }

        /// Adds task "resultT memfun(arg1T,arg2T) const" in process owning item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2,  const TaskAttributes& attr = TaskAttributes()) const {
            return const_cast<containerT*>(this)->task(key,memfun,arg1,arg2,attr);
        }

        /// Adds task "resultT memfun(arg1T,arg2T,arg3T) const" in process owning item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const TaskAttributes& attr = TaskAttributes()) const {
            return const_cast<containerT*>(this)->task(key,memfun,arg1,arg2,arg3,attr);
        }

        /// Adds task "resultT memfun(arg1T,arg2T,arg3T, arg4T) const" in process owning item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const TaskAttributes& attr = TaskAttributes()) const {
            return const_cast<containerT*>(this)->task(key,memfun,arg1,arg2,arg3,arg4,attr);
        }

        /// Adds task "resultT memfun(arg1T,arg2T,arg3T,arg4T,arg5T) const" in process owning item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5, const TaskAttributes& attr = TaskAttributes()) const {
            return const_cast<containerT*>(this)->task(key,memfun,arg1,arg2,arg3,arg4,arg5,attr);
        }

        /// Adds task "resultT memfun(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T) const" in process owning item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5, const arg6T& arg6, const TaskAttributes& attr = TaskAttributes()) const {
            return const_cast<containerT*>(this)->task(key,memfun,arg1,arg2,arg3,arg4,arg5,arg6,attr);
        }

        /// Adds task "resultT memfun(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T) const" in process owning item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5, const arg6T& arg6, const arg7T& arg7, const TaskAttributes& attr = TaskAttributes()) const {
            return const_cast<containerT*>(this)->task(key,memfun,arg1,arg2,arg3,arg4,arg5,arg6,arg7,attr);
        }


        /// (de)Serialize --- *Local* data only to/from anything *except* Buffer*Archive and Parallel*Archive

        /// Advisable for *you* to fence before and after this to ensure consistency
        template <typename Archive>
        void serialize(const Archive& ar) {
            //
            // !! If you change the format of this stream make sure that
            // !! the parallel in/out archive below is compatible
            //
            const long magic = 5881828; // Sitar Indian restaurant in Knoxville
            unsigned long count = 0;
            check_initialized();

            if (Archive::is_output_archive) {
                ar & magic;
                for (iterator it=begin(); it!=end(); ++it) count++;
                ar & count;
                for (iterator it=begin(); it!=end(); ++it) ar & *it;
            }
            else {
                long cookie = 0l;
                ar & cookie;
                MADNESS_ASSERT(cookie == magic);
                ar & count;
                while (count--) {
                    pairT datum;
                    ar & datum;
                    replace(datum);
                }
            }
        }

        /// (de)Serialize --- !! ONLY for purpose of interprocess communication

        /// This just writes/reads the unique id to/from the Buffer*Archive.
        void serialize(const archive::BufferOutputArchive& ar) {
            check_initialized();
            ar & static_cast<WorldObject<implT>*>(p.get());
        }

        /// (de)Serialize --- !! ONLY for purpose of interprocess communication

        /// This just writes/reads the unique id to/from the Buffer*Archive.
        void serialize(const archive::BufferInputArchive& ar) {
            WorldObject<implT>* ptr = nullptr;
            ar & ptr;
            MADNESS_ASSERT(ptr);

#ifdef MADNESS_DISABLE_SHARED_FROM_THIS
            p.reset(static_cast<implT*>(ptr), [] (implT *p_) -> void {});
#else
            p = static_cast<implT*>(ptr)->shared_from_this();
#endif // MADNESS_DISABLE_SHARED_FROM_THIS
        }

        /// Returns the associated unique id ... must be initialized
        const uniqueidT& id() const {
            check_initialized();
            return p->id();
        }

        /// Destructor passes ownership of implementation to world for deferred cleanup
        virtual ~WorldContainer() {
            detail::deferred_cleanup(p->get_world(), p);
        }

        friend void swap<>(WorldContainer&, WorldContainer&);
    };

    /// Swaps the content of two WorldContainer objects. It should be called on all nodes.

    /// \ingroup worlddc
    template <typename keyT, typename valueT, typename hashfunT>
    void swap(WorldContainer<keyT, valueT, hashfunT>& dc0, WorldContainer<keyT, valueT, hashfunT>& dc1) {
      std::swap(dc0.p, dc1.p);
    }

    namespace archive {

        /// Write container to parallel archive

        /// specialization for parallel serialization of a WorldContainer:
        /// all threads on each process serialize some values into a buffer, which gets concatenated
        /// and finally serialized to localarchive (aka VectorOutputArchive).
        template <class keyT, class valueT>
        struct ArchiveStoreImpl< ParallelOutputArchive<VectorOutputArchive>, WorldContainer<keyT,valueT> > {
            static void store(const ParallelOutputArchive<VectorOutputArchive>& ar, const WorldContainer<keyT,valueT>& t) {
                using localarchiveT = VectorOutputArchive;
                const long magic = -5881828; // Sitar Indian restaurant in Knoxville (negative to indicate parallel!)
                typedef WorldContainer<keyT,valueT> dcT;
                using const_iterator = typename dcT::const_iterator;

                const size_t default_size = 100*1024*1024;

                World* world = ar.get_world();
                world->gop.fence();

                std::vector<unsigned char> v;
                v.reserve(default_size);
                size_t count = 0;

                class op : public TaskInterface {
                    const size_t ntasks;
                    const size_t taskid;
                    const dcT& t;
                    std::vector<unsigned char>& vtotal;
                    size_t& total_count;
                    Mutex& mutex;

                public:
                    op(size_t ntasks, size_t taskid, const dcT& t, std::vector<unsigned char>& vtotal, size_t& total_count, Mutex& mutex)
                        : ntasks(ntasks), taskid(taskid), t(t), vtotal(vtotal), total_count(total_count), mutex(mutex) {}
                    void run(World& world) {
                        std::vector<unsigned char> v;
                        v.reserve(std::max(size_t(1024*1024),vtotal.capacity()/ntasks));
                        VectorOutputArchive var(v);
                        const_iterator it=t.begin();
                        size_t count = 0;
                        size_t n = 0;
                        /// threads serialize round-robin over the container
                        while (it!=t.end()) {
                            if ((n%ntasks) == taskid) {
                                var & *it;
                                ++count;
                            }
                            ++it;
                            n++;
                        }

                        // concatenate the buffers from each thread
                        if (count) {
                            mutex.lock();
                            vtotal.insert(vtotal.end(), v.begin(), v.end());
                            total_count += count;
                            mutex.unlock();
                        }
                    }
                };

                Mutex mutex;
                size_t ntasks = std::max(size_t(1), ThreadPool::size());
                for (size_t taskid=0; taskid<ntasks; taskid++)
                    world->taskq.add(new op(ntasks, taskid, t, v, count, mutex));
                world->gop.fence();

                // Gather all buffers to process 0
                // first gather all of the sizes and counts to a vector in process 0
                int size = v.size();
                std::vector<int> sizes(world->size());
                MPI_Gather(&size, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, world->mpi.comm().Get_mpi_comm());
                world->gop.sum(count); // just need total number of elements

                // build the cumulative sum of sizes
                std::vector<int> offsets(world->size());
                offsets[0] = 0;
                for (int i=1; i<world->size(); ++i) offsets[i] = offsets[i-1] + sizes[i-1];
                int total_size = offsets.back() + sizes.back();

                // gather the vector of data v from each process to process 0
                unsigned char* all_data=0;
                if (world->rank() == 0) {
                    all_data = new unsigned char[total_size];
                }
                MPI_Gatherv(v.data(), v.size(), MPI_BYTE, all_data, sizes.data(), offsets.data(), MPI_BYTE, 0, world->mpi.comm().Get_mpi_comm());

                if (world->rank() == 0) {
                    auto& localar = ar.local_archive();
                    localar & magic & 1; // 1 client
                    // localar & t;
                    ArchivePrePostImpl<localarchiveT,dcT>::preamble_store(localar);
                    localar & -magic & count;
                    localar.store(all_data, total_size);
                    ArchivePrePostImpl<localarchiveT,dcT>::postamble_store(localar);

                    delete[] all_data;
                }
                world->gop.fence();
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
        struct ArchiveStoreImpl< ParallelOutputArchive<localarchiveT>, WorldContainer<keyT,valueT> > {
            static void store(const ParallelOutputArchive<localarchiveT>& ar, const WorldContainer<keyT,valueT>& t) {
                const long magic = -5881828; // Sitar Indian restaurant in Knoxville (negative to indicate parallel!)
                typedef WorldContainer<keyT,valueT> dcT;
                // typedef typename dcT::const_iterator iterator; // unused?
                typedef typename dcT::pairT pairT;
                World* world = ar.get_world();
                Tag tag = world->mpi.unique_tag();
                ProcessID me = world->rank();
                if (ar.dofence()) world->gop.fence();
                if (ar.is_io_node()) {
                    auto& localar = ar.local_archive();
                    localar & magic & ar.num_io_clients();
                    for (ProcessID p=0; p<world->size(); ++p) {
                        if (p == me) {
                            localar & t;
                        }
                        else if (ar.io_node(p) == me) {
                            world->mpi.Send(int(1),p,tag); // Tell client to start sending
                            archive::MPIInputArchive source(*world, p);
                            long cookie = 0l;
                            unsigned long count = 0ul;

                            ArchivePrePostImpl<localarchiveT,dcT>::preamble_store(localar);

                            source & cookie & count;
                            localar & cookie & count;
                            while (count--) {
                                pairT datum;
                                source & datum;
                                localar & datum;
                            }

                            ArchivePrePostImpl<localarchiveT,dcT>::postamble_store(localar);
                        }
                    }
                }
                else {
                    ProcessID p = ar.my_io_node();
                    int flag;
                    world->mpi.Recv(flag,p,tag);
                    MPIOutputArchive dest(*world, p);
                    dest & t;
                    dest.flush();
                }
                if (ar.dofence()) world->gop.fence();
            }
        };

        template <class keyT, class valueT, class localarchiveT>
        struct ArchiveLoadImpl< ParallelInputArchive<localarchiveT>, WorldContainer<keyT,valueT> > {
            /// Read container from parallel archive

            /// \ingroup worlddc
            /// See store method above for format of file content.
            /// !!! We presently ASSUME that the number of writers and readers are
            /// the same.  This is frustrating but not a show stopper since you
            /// can always run a separate job to copy to a different number.
            ///
            /// The IO node simply reads all data and inserts entries.
            static void load(const ParallelInputArchive<localarchiveT>& ar, WorldContainer<keyT,valueT>& t) {
                const long magic = -5881828; // Sitar Indian restaurant in Knoxville (negative to indicate parallel!)
                // typedef WorldContainer<keyT,valueT> dcT; // unused
                // typedef typename dcT::iterator iterator; // unused
                // typedef typename dcT::pairT pairT; // unused
                World* world = ar.get_world();
                if (ar.dofence()) world->gop.fence();
                if (ar.is_io_node()) {
                    long cookie = 0l;
                    int nclient = 0;
                    auto& localar = ar.local_archive();
                    localar & cookie & nclient;
                    MADNESS_CHECK(cookie == magic);
                    while (nclient--) {
                        localar & t;
                    }
                }
                if (ar.dofence()) world->gop.fence();
            }
        };
    }

}

///@}

#endif // MADNESS_WORLD_WORLDDC_H__INCLUDED
