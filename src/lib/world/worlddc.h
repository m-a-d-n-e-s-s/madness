/*
  This file is part of MADNESS.
  
  Copyright (C) <2007> <Oak Ridge National Laboratory>
  
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

  
#ifndef WORLD_DC
#define WORLD_DC

/// \file worlddc.h
/// \brief Implements WorldContainer

#include <world/parar.h>
#include <world/worldhashmap.h>

namespace madness {
        
    template <typename keyT, 
              typename valueT>
    class WorldContainer;
    
    template <typename keyT, 
              typename valueT>
    class WorldContainerImpl;

    /// Interface to be provided by any process map
    template <typename keyT>
    class WorldDCPmapInterface {
    public:
        virtual ProcessID owner(const keyT& key) const = 0;
        virtual ~WorldDCPmapInterface() {}
	virtual void print() const {}
    };

    /// Default process map is "random" using madness::hash(key)
    template <typename keyT>
    class WorldDCDefaultPmap : public WorldDCPmapInterface<keyT> {
    private:
        int nproc;
        
    public:
        WorldDCDefaultPmap(World& world) {
            nproc = world.mpi.nproc();
        }
        
        ProcessID owner(const keyT& key) const {
            if (nproc == 1) return 0;
            return hash(key)%nproc;
        }
    };
    
    
    /// Iterator for distributed container wraps the local iterator
    template <class internal_iteratorT, class pairT>
    class WorldContainerIterator {
    private:
        internal_iteratorT  it;       ///< Iterator from local container
        bool is_local;                ///< If true we are using the local container
        pairT value;                  ///< If (!is_local) holds the remote value

    public:
        /// Default constructor makes a local uninitialized value
        explicit WorldContainerIterator() 
            : it(), is_local(true), value()
        {}

        /// Initializes from a local iterator
        explicit WorldContainerIterator(const internal_iteratorT& it)
            : it(it), is_local(true), value() 
        {}

        /// Initializes to cache a remote value
        explicit WorldContainerIterator(const pairT& value)
            : it(), is_local(false), value(value)
        {}

        WorldContainerIterator(const WorldContainerIterator& other) 
        {
            *this = other;
        }

        /// Assignment 
        WorldContainerIterator& operator=(const WorldContainerIterator& other) {
            if (this != &other) {
                is_local = other.is_local;
                if (other.is_local) {
                    it = other.it;
                }
                else {
                    // Sigh ... there does not seem a kosher way to do this
                    *const_cast<REMCONST(typename pairT::first_type)*>(&value.first) = other.value.first;
                    *const_cast<REMCONST(typename pairT::second_type)*>(&value.second) = other.value.second;
                }
            }
            return *this;
        }
        

        /// Determines if two iterators are identical
        bool operator==(const WorldContainerIterator& other) const {
            if (is_local) {
                if (other.is_local) {
                    return it == other.it;
                }
                else {
                    return false;
                }
            }
            else if (other.is_local) {
                return false;
            }
            else {
                return value.first == other.value.first;
            }
        }
        

        /// Determines if two iterators are different
        bool operator!=(const WorldContainerIterator& other) const {
            return !(*this == other);
        }
        

        /// Pre-increment of an iterator (i.e., ++it) --- \em local iterators only
        
        /// Trying to increment a remote iterator will throw
        WorldContainerIterator& operator++() {
            MADNESS_ASSERT(is_local);
            ++it;
            return *this;
        }
        
        /// Iterators dereference to std::pair<const keyT,valueT>
        const pairT* operator->() const {
            if (is_local) return it.operator->();
            else return &value;
        }
        
        /// Iterators dereference to std::pair<const keyT,valueT>
        pairT* operator->() {
            if (is_local) return it.operator->();
            else return &value;
            
        }
        
        /// Iterators dereference to const std::pair<const keyT,valueT>
        const pairT& operator*() const {
            if (is_local) return *it;
            else return value; 
       }
        
        /// Iterators dereference to std::pair<const keyT,valueT>
        pairT& operator*() {
            if (is_local) return *it;
            else return value;
        }
        
        /// Private: (or should be) Returns iterator of internal container
        const internal_iteratorT& get_internal_iterator() const {
            return it;
        }
        
        /// Returns true if this is non-local or cached value
        bool is_cached() const {
            return !is_local;
        }
        
        virtual ~WorldContainerIterator() 
        {}
        
        template <typename Archive>
        void serialize(const Archive& ar) {
            throw "Serializing DC iterator ... why?";
        }
    };

    template <class internal_iteratorT, class pairT>
    std::ostream& operator<<(std::ostream& s, const WorldContainerIterator<internal_iteratorT, pairT>& it) {
        s << "WCIterator(" << *it << ")";
        return s;
    }

    
    /// Implementation of distributed container to enable PIMPL
    template <typename keyT, 
              typename valueT>
    class WorldContainerImpl 
        : public WorldObject< WorldContainerImpl<keyT, valueT> >
        , private NO_DEFAULTS 
    {
    public:
        typedef typename std::pair<const keyT,valueT> pairT;
        typedef const pairT const_pairT;
        typedef WorldContainerImpl<keyT,valueT> implT;
        
        template <typename T> 
        struct DCLocalHash {
            std::size_t operator()(const T& t) const {
                return hash(t);
            }
        };
        typedef ConcurrentHashMap< keyT,valueT,DCLocalHash<keyT> > internal_containerT;
        
        typedef typename internal_containerT::iterator internal_iteratorT;
        typedef typename internal_containerT::const_iterator internal_const_iteratorT;
        typedef typename internal_containerT::accessor accessor;
        typedef typename internal_containerT::const_accessor const_accessor;
        typedef WorldContainerIterator<internal_iteratorT,pairT> iteratorT;
        typedef WorldContainerIterator<internal_iteratorT,pairT> iterator;
        typedef WorldContainerIterator<internal_const_iteratorT, const_pairT> const_iteratorT;
        typedef WorldContainerIterator<internal_const_iteratorT, const_pairT> const_iterator;
        
        friend class WorldContainer<keyT,valueT>;
        
    private:
        
        WorldContainerImpl();   // Inhibit default constructor
        
        const SharedPtr< WorldDCPmapInterface<keyT> > pmap;       ///< Function/class to map from keys to owning process
        const ProcessID me;           ///< My MPI rank
        internal_containerT local;    ///< Locally owned data
        const iterator end_iterator;          ///< For fast return of end
        const const_iterator end_const_iterator; ///< For fast return of end
        
	template <typename containerT, typename datumT> 
	inline
	static
        typename containerT::iterator replace(containerT& c, const datumT& d) {
            std::pair<typename containerT::iterator,bool> p = c.insert(d);
            if (!p.second) p.first->second = d.second;   // Who's on first?
            return p.first;
	}

        /// Handles find request
        Void find_handler(ProcessID requestor, const keyT& key, const RemoteReference< FutureImpl<iterator> >& ref) {
            internal_iteratorT r = local.find(key);
            if (r == local.end()) {
                //print("find_handler: failure:", key);
                send(requestor, &implT::find_failure_handler, ref);
            }
            else {
                //print("find_handler: success:", key, r->first, r->second);
                send(requestor, &implT::find_success_handler, ref, *r);
            }
            return None;
        }
        
        /// Handles successful find response
        Void find_success_handler(const RemoteReference< FutureImpl<iterator> >& ref, const pairT& datum) {
            FutureImpl<iterator>* f = ref.get();
            f->set(iterator(datum));
            //print("find_success_handler: success:", datum.first, datum.second, f->get()->first, f->get()->second);
            ref.dec(); // Matching inc() in find() where ref was made
            return None;
        }
        
        /// Handles unsuccessful find response
        Void find_failure_handler(const RemoteReference< FutureImpl<iterator> >& ref) {
            FutureImpl<iterator>* f = ref.get();
            f->set(end());
            //print("find_failure_handler");
            ref.dec(); // Matching inc() in find() where ref was made
            return None;
        }
        
    public:
        
        WorldContainerImpl(World& world, 
                           const SharedPtr< WorldDCPmapInterface<keyT> >& pmap, 
                           bool do_pending)
            : WorldObject< WorldContainerImpl<keyT, valueT> >(world)
            , pmap(pmap)
            , me(world.mpi.rank())
            , local()
            , end_iterator(local.end())
            , end_const_iterator(const_cast<const internal_containerT&>(local).end())
        {
            if (do_pending) this->process_pending();
	}
        
        
        ~WorldContainerImpl() {
        }

	const SharedPtr< WorldDCPmapInterface<keyT> >& get_pmap() const {
	    return pmap;
	}

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

        Void insert(const pairT& datum) {
            ProcessID dest = owner(datum.first);
            if (dest == me) {
  	        replace(local,datum);
            }
            else {
                send(dest, &implT::insert, datum);
            }
            return None;
        }
        
        
        void clear() {
            local.clear();
        }
        
        
        Void erase(const keyT& key) {
            ProcessID dest = owner(key);
            if (dest == me) {
                local.erase(key);
            }
            else {
		Void (implT::*eraser)(const keyT&) = &implT::erase;
                send(dest, eraser, key);
            }                
            return None;
        }
        
        
        void erase(const iterator& it) {
            MADNESS_ASSERT(!it.is_cached());
            if (it == end()) {
                return;
            }
            else { 
                erase(it->first);
            }
        }
        
        void erase(const iterator& start, const iterator& finish) {
            iterator it = start;
            while (it!=finish && it!=end()) {
                iterator last=it;
                ++it;
                erase(last);
            }
        }
        
        iterator begin() {
            return iterator(local.begin());
        }
        
        const_iterator begin() const {
            return const_iterator(local.begin());
        }
        
        const iterator& end() {
            return end_iterator;
        }
        
        const const_iterator& end() const {
            return end_const_iterator;
        }
        
        
        Future<const_iterator> find(const keyT& key) const {
            // Ugliness here to avoid replicating find() and
            // associated handlers for const.  Assumption is that
            // const and non-const iterators are idential except for
            // const attribute ... at some point probably need to do
            // the right thing.
            Future<iterator> r = const_cast<implT*>(this)->find(key);
            return *(Future<const_iterator>*)(&r);
        }


        Future<iterator> find(const keyT& key) {
            ProcessID dest = owner(key);
            if (dest == me) {
                return Future<iterator>(iterator(local.find(key)));
            }
            else {
                Future<iterator> result; 
                send(dest, &implT::find_handler, me, key, result.remote_ref(this->world));
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
            local.insert(acc, key);
            return (acc->second.*memfun)();
        }
        
        // Used to forward call to item member function
        template <typename memfunT, typename arg1T>
        MEMFUN_RETURNT(memfunT)        
        itemfun(const keyT& key, memfunT memfun, const arg1T& arg1) {
            accessor acc;
            local.insert(acc, key);
            return (acc->second.*memfun)(arg1);
        }
        
        // Used to forward call to item member function
        template <typename memfunT, typename arg1T, typename arg2T>
        MEMFUN_RETURNT(memfunT)        
        itemfun(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2) {
            accessor acc;
            local.insert(acc, key);
            return (acc->second.*memfun)(arg1,arg2);
        }
        
        // Used to forward call to item member function
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        MEMFUN_RETURNT(memfunT)        
        itemfun(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3) {
            accessor acc;
            local.insert(acc, key);
            return (acc->second.*memfun)(arg1,arg2,arg3);
        }
        
        // Used to forward call to item member function
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        MEMFUN_RETURNT(memfunT)        
        itemfun(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg3T& arg4) {
            accessor acc;
            local.insert(acc, key);
            return (acc->second.*memfun)(arg1,arg2,arg3,arg4);
        }
        
        // Used to forward call to item member function
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        MEMFUN_RETURNT(memfunT)        
        itemfun(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg3T& arg4, const arg5T& arg5) {
            accessor acc;
            local.insert(acc, key);
            return (acc->second.*memfun)(arg1,arg2,arg3,arg4,arg5);
        }
    };
    
    
    /// Makes a distributed container with specified attributes
    
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
    template <typename keyT, 
              typename valueT>
    class WorldContainer : public ParallelSerializableObject {
    public:
        typedef WorldContainer<keyT,valueT> containerT;
        typedef WorldContainerImpl<keyT,valueT> implT;
        typedef typename implT::pairT pairT;
        typedef typename implT::iterator iterator;
        typedef typename implT::const_iterator const_iterator;
        typedef typename implT::accessor accessor;
        typedef typename implT::const_accessor const_accessor;
        typedef Future<iterator> futureT;
        typedef Future<const_iterator> const_futureT;
        
    private:
        SharedPtr<implT> p;
        
        inline void check_initialized() const {
            MADNESS_ASSERT(p);
        }
    public:
        
        /// Makes an uninitialized container (no communication)
        
        /// The container is useless until assigned to from a fully
        /// constructed container.  There is no need to worry about
        /// default constructors being executed in order.
        WorldContainer() 
            : p(0)
        {}
        
        
        /// Makes an initialized, empty container with default data distribution (no communication)
        
        /// A unique ID is associated with every distributed container
        /// within a world.  In order to avoid synchronization when
        /// making a container, we have to assume that all processes
        /// execute this constructor in the same order (does not apply
        /// to the non-initializing, default constructor).
        WorldContainer(World& world, bool do_pending=true) 
            : p(new implT(world, 
                          SharedPtr< WorldDCPmapInterface<keyT> >(new WorldDCDefaultPmap<keyT>(world)), 
                          do_pending))
        {}
        
        /// Makes an initialized, empty container (no communication)
        
        /// A unique ID is associated with every distributed container
        /// within a world.  In order to avoid synchronization when
        /// making a container, we have to assume that all processes
        /// execute this constructor in the same order (does not apply
        /// to the non-initializing, default constructor).
        WorldContainer(World& world, 
                       const SharedPtr< WorldDCPmapInterface<keyT> >& pmap, 
                       bool do_pending=true) 
            : p(new implT(world, pmap, do_pending))
        {}
        
        
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
            return p->world;
        }
        
        
        /// Inserts key+value pair (non-blocking communication if key not local)
        void replace(const pairT& datum) {
            check_initialized();
            p->insert(datum);
        }
        
        
        /// Inserts key+value pair (non-blocking communication if key not local)
        void replace(const keyT& key, const valueT& value) {
            replace(pairT(key,value));
        }

        
        /// Provides write access to LOCAL value by key ... always returns false for remote
        bool find(accessor& acc, const keyT& key) {
            check_initialized();
            return p->find(acc,key);
        }
        

        /// Provides read access to LOCAL value by key ... always returns false for remote
        bool find(const_accessor& acc, const keyT& key) const {
            check_initialized();
            return p->find(acc,key);
        }


        /// Provides write access to LOCAL value by key ... always returns false for remote
        bool insert(accessor& acc, const keyT& key) {
            check_initialized();
            return p->insert(acc,key);
        }


        /// Provides read access to LOCAL value by key ... always returns false for remote
        bool insert(const_accessor& acc, const keyT& key) {
            check_initialized();
            return p->insert(acc,key);
        }


        /// Inserts pairs (non-blocking communication if key(s) not local)
        template <typename input_iterator>
        void replace(input_iterator& start, input_iterator& end) {
            check_initialized();
            for_each(start,end,bind1st(mem_fun(&containerT::insert),this));
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
        const iterator& end() {
            check_initialized();
            return p->end();
        }
        
        /// Returns an iterator past the end of the \em local data (no communication)
        const const_iterator& end() const {
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
	inline const SharedPtr< WorldDCPmapInterface<keyT> >& get_pmap() const {
	    return p->get_pmap();
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
            MEMFUN_RETURNT(memfunT) (implT::*itemfun)(const keyT&, memfunT) = &implT:: template itemfun<memfunT>;
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
            //MEMFUN_RETURNT(memfunT) (implT::*itemfun)(const keyT&, memfunT, const arg1T&) = &implT:: template itemfun<memfunT,arg1T>;
            //return p->send(owner(key), itemfun, key, memfun, arg1);
            return p->send(owner(key), 
                           static_cast<MEMFUN_RETURNT(memfunT) (implT::*)(const keyT&, memfunT, const arg1T&)>(&implT:: template itemfun<memfunT,arg1T>), 
                           key, memfun, arg1);
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
            //MEMFUN_RETURNT(memfunT) (implT::*itemfun)(const keyT&, memfunT, const arg1T&, const arg2T&) = &implT:: template itemfun<memfunT,arg1T,arg2T>;
            //return p->send(owner(key), itemfun, key, memfun, arg1, arg2);
            return p->send(owner(key), 
                           static_cast<MEMFUN_RETURNT(memfunT) (implT::*)(const keyT&, memfunT, const arg1T&, const arg2T&)>(&implT:: template itemfun<memfunT,arg1T,arg2T>), 
                           key, memfun, arg1, arg2);
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
            MEMFUN_RETURNT(memfunT) (implT::*itemfun)(const keyT&, memfunT, const arg1T&, const arg2T&, const arg3T&) = &implT:: template itemfun<memfunT,arg1T,arg2T,arg3T>;
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
            MEMFUN_RETURNT(memfunT) (implT::*itemfun)(const keyT&, memfunT, const arg1T&, const arg2T&, const arg3T&, const arg4T&) = &implT:: template itemfun<memfunT,arg1T,arg2T,arg3T,arg4T>;
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
            MEMFUN_RETURNT(memfunT) (implT::*itemfun)(const keyT&, memfunT, const arg1T&, const arg2T&, const arg3T&, const arg4T&, const arg5T&) = &implT:: template itemfun<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T>;
            return p->send(owner(key), itemfun, key, memfun, arg1, arg2, arg3, arg4, arg5);
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
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > 
        send(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4) const {
            return const_cast<containerT*>(this)->send(key,memfun,arg1,arg2,arg3,arg4);
        }

        /// Sends message "resultT memfun(arg1T,arg2T,arg3T) const" to item (non-blocking comm if remote)

        /// The method executes with a write lock on the item.
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > 
        send(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5) const {
            return const_cast<containerT*>(this)->send(key,memfun,arg1,arg2,arg3,arg4,arg5);
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
            MEMFUN_RETURNT(memfunT) (implT::*itemfun)(const keyT&, memfunT) = &implT:: template itemfun<memfunT>;
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
            MEMFUN_RETURNT(memfunT) (implT::*itemfun)(const keyT&, memfunT, const a1T&) = &implT:: template itemfun<memfunT,a1T>;
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
            MEMFUN_RETURNT(memfunT) (implT::*itemfun)(const keyT&, memfunT, const a1T&, const a2T&) = &implT:: template itemfun<memfunT,a1T,a2T>;
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
            MEMFUN_RETURNT(memfunT) (implT::*itemfun)(const keyT&, memfunT, const a1T&, const a2T&, const a3T&) = &implT:: template itemfun<memfunT,a1T,a2T,a3T>;
            return p->task(owner(key), itemfun, key, memfun, arg1, arg2, arg3, attr);
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
                long cookie;
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
        void serialize(const BufferOutputArchive& ar) {
            check_initialized();
            ar & static_cast<WorldObject<implT>*>(p.get());
        }

        /// (de)Serialize --- !! ONLY for purpose of interprocess communication
        
        /// This just writes/reads the unique id to/from the Buffer*Archive.
        void serialize(const BufferInputArchive& ar) {
            WorldObject<implT>* ptr;
            ar & ptr;
            MADNESS_ASSERT(ptr);
            p = SharedPtr<implT>(static_cast<implT*>(ptr),false);  // use_count will be 0, which is good
        }

        /// Returns the associated unique id ... must be initialized
        const uniqueidT& id() const {
            check_initialized();
            return p->id();
        }
        
        /// Destructor passes ownership of implementation to world for deferred cleanup
        virtual ~WorldContainer() {
            if (p) p->world.deferred_cleanup(p);
        }
    };

    namespace archive {
        template <class keyT, class valueT>
        struct ArchiveStoreImpl< BufferOutputArchive, WorldContainer<keyT,valueT>* > {
            static void store(const BufferOutputArchive& ar, WorldContainer<keyT,valueT>* const & t) {
                ar & *t;
            }
        };

        template <class keyT, class valueT>
        struct ArchiveLoadImpl< BufferInputArchive, WorldContainer<keyT,valueT>* > {
            static void load(const BufferInputArchive& ar, WorldContainer<keyT,valueT>*& t) {
                ar & *t;
            }
        };


        /// Write container to parallel archive with optional fence

        /// Each node (process) is served by a designated IO node.
        /// The IO node has a binary local file archive to which is
        /// first written a cookie and the number of servers.  The IO
        /// node then loops thru all of its clients and in turn tells
        /// each to write its data over an MPI stream, which is copied
        /// directly to the output file.  The stream contents are then
        /// cookie, #clients, foreach client (usual sequential archive).
        ///
        /// If ar.dofence() is true (default) fence is invoked before and
        /// after the IO. The fence is optional but it is of course
        /// necessary to be sure that all updates have completed
        /// before doing IO, and that all IO has completed before
        /// subsequent modifications. Also, there is always at least
        /// some synchronization between a client and its IO server.
        template <class keyT, class valueT>
        struct ArchiveStoreImpl< ParallelOutputArchive, WorldContainer<keyT,valueT> > {
            static void store(const ParallelOutputArchive& ar, const WorldContainer<keyT,valueT>& t) {
                const long magic = -5881828; // Sitar Indian restaurant in Knoxville (negative to indicate parallel!)
                typedef WorldContainer<keyT,valueT> dcT;
                typedef typename dcT::const_iterator iterator;
                typedef typename dcT::pairT pairT;
                World* world = ar.get_world();
                Tag tag = world->mpi.unique_tag();
                ProcessID me = world->rank();
                if (ar.dofence()) world->gop.fence();
                if (ar.is_io_node()) {
                    BinaryFstreamOutputArchive& localar = ar.local_archive();
                    localar & magic & ar.num_io_clients();
                    for (ProcessID p=0; p<world->size(); p++) {
                        if (p == me) {
                            localar & t;
                        }
                        else if (ar.io_node(p) == me) {
                            world->mpi.Send(int(1),p,tag); // Tell client to start sending
                            MPIInputArchive source(*world, p);
                            long cookie;
                            unsigned long count;

                            ArchivePrePostImpl<BinaryFstreamOutputArchive,dcT>::preamble_store(localar);

                            source & cookie & count;
                            localar & cookie & count;
                            while (count--) {
                                pairT datum;
                                source & datum;
                                localar & datum;
                            }

                            ArchivePrePostImpl<BinaryFstreamOutputArchive,dcT>::postamble_store(localar);
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

        template <class keyT, class valueT>
        struct ArchiveLoadImpl< ParallelInputArchive, WorldContainer<keyT,valueT> > {
            /// Read container from parallel archive

            /// See store method above for format of file content.
            /// !!! We presently ASSUME that the number of writers and readers are
            /// the same.  This is frustrating but not a show stopper since you
            /// can always run a separate job to copy to a different number.
            /// 
            /// The IO node simply reads all data and inserts entries.
            static void load(const ParallelInputArchive& ar, WorldContainer<keyT,valueT>& t) {
                const long magic = -5881828; // Sitar Indian restaurant in Knoxville (negative to indicate parallel!)
                typedef WorldContainer<keyT,valueT> dcT;
                typedef typename dcT::iterator iterator;
                typedef typename dcT::pairT pairT;
                World* world = ar.get_world();
                if (ar.dofence()) world->gop.fence();
                if (ar.is_io_node()) {
                    long cookie;
                    int nclient;
                    BinaryFstreamInputArchive& localar = ar.local_archive();
                    localar & cookie & nclient;
                    MADNESS_ASSERT(cookie == magic);
                    while(nclient--) {
                        localar & t;
                    }
                }
                if (ar.dofence()) world->gop.fence();
            }
        };
    }

}

#endif
