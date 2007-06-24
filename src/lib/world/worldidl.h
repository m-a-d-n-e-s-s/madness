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

  $LastChangedDate$
  $Rev$
*/

  
#ifndef WORLD_IDL_H
#define WORLD_IDL_H

/// \file worldidl.h
/// \brief Provides global equivalent of SharedPtr with remote method invocation 

/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
/// DEPRECATED ... THIS STUFF IS ABOUT TO DIE
///
/// SharedPtr makes local pointers "safe" by reference counting
/// all uses and only deleting the associated object when there are no
/// more references.  This is the first third of what WorldSharedPtr
/// does.  Locally, it behaves pretty much like SharedPtr but
/// references are counted world-wide, not just locally.  You can pass
/// a WorldSharedPtr to another process and the local reference count
/// will be incremented and when the remote pointer is finally
/// destroyed, the local reference count will be decremented.
///
/// The second third of WorldSharedPtr's capability is a standard interface
/// for remote access to objects (prefetch(), flush_cache(), clear_cache(), probe(),
/// is_local()) so that a remote object may be accessed efficiently.
///
/// The final third is a standard set of macros to enable you to declare
/// new methods in the class and have their arguments and results 
/// forwarded to/from remote objects.  There are some limitations on
/// the arguments and results of methods, but essentially if it 
/// is serializable it should work.



// Nesting of macros here to force C preprocessor prescan to evaluate __LINE__
#define METHOD_HANDLER(method) UNIQUE_NAME(method##_handler)
#define UNIQUE_NAME(method) _UNIQUE_NAME(method,__LINE__)
#define _UNIQUE_NAME(method,line) __UNIQUE_NAME(method,line)
#define __UNIQUE_NAME(method,line) method##_##line

#define STRING(A) #A


// Empty template to force use of the macros
template <typename T> 
class WorldSharedPtr {};

template <typename T> 
struct WorldSharedPtrClassNames {};

// Forward declaration
template <typename T> class WorldSharedPtrImpl;


// Forward declaration
template <typename T> class WorldSharedPtrBase;


/// Holds the data for the WorldSharedPtr class

/// We want to use Pimpl for WorldSharedPtr in order to have a shallow
/// copy and to enable full reference counting of both local and
/// remote references.  However, this is tedious to do inside the
/// macros, so we put all of the data into this class and all the
/// implementation into WorldSharedPtr/Base which contains a shared
/// pointer to a WorldSharedPtrImpl.  Actually, even without the
/// macros we need to separate the data from the implementation so
/// that we can do reference counting on the data.
template <typename T>
class WorldSharedPtrImpl : NO_DEFAULTS {
    friend class WorldSharedPtr<T>;
    friend class WorldSharedPtrBase<T>;

public:
    typedef RemoteReference< WorldSharedPtrImpl<T> > refT;

private:
    World* world;         //< 42
    SharedPtr<T> local;   //< Holds pointer to any local data
    refT remote;          //< Reference to any remote data
    bool prefetch_pending;//< True if a prefetch is in progress

    WorldSharedPtrImpl(); // Inhibit default constructor

public:
    /// Takes ownership of a local pointer
    explicit WorldSharedPtrImpl(World& world, T* t) 
        : world(&world), local(t), remote(), prefetch_pending(false) {};

    /// Takes a reference to a local shared pointer
    explicit WorldSharedPtrImpl(World& world, const SharedPtr<T>& t) 
        : world(&world), local(t), remote(), prefetch_pending(false) {};

    /// Takes ownership of a reference to a remote object

    /// The remote reference must be initialized otherwise identifying
    /// the associated world will fail.
    explicit WorldSharedPtrImpl(const refT& ref) 
        : world(World::world_from_id(ref.world_id())), local(0), remote(ref), prefetch_pending(false) {
    };

    /// Returns true if data is local
    bool is_local() const {
        return !remote;
    };

    virtual ~WorldSharedPtrImpl() {
        if (remote) remote.dec();
    };
};

/// Implements standard functionality for WorldSharedPtr

/// Easier to put it here than in the macros which are then only used
/// to provide remote access to methods of class T.
template <class T>
class WorldSharedPtrBase {                                  
protected:
    static bool debug;
    SharedPtr< WorldSharedPtrImpl<T> > p;                                          

    static const char* classname() {
        return WorldSharedPtrClassNames<T>::classname();
    };

    /// Handles incoming prefetch requests
    static void prefetch_request_handler(World& world, ProcessID src, const AmArg& inarg) { 
        refT target, source;
        inarg.unstuff(target,source);
        if (debug) madness::print(classname(), "prefetch_han");
        implT* impl = target.get();        
        MADNESS_ASSERT(impl);
        if (impl->is_local()) {
            LongAmArg* out = new LongAmArg;
            size_t nbyte = out->stuff(source, *impl->local);
            world.am.send_long_managed(source.owner(), 
                                       baseT::prefetch_reply_handler, 
                                       out, 
                                       nbyte);
        }
        else {
            // Forward the request with reply going to original requestor
            MADNESS_ASSERT(impl->remote);
            if (debug) madness::print(classname(),"prefetch_for");
            world.am.send(impl->remote.owner(), 
                          baseT::prefetch_request_handler, 
                          AmArg(impl->remote,source));
            world.am.fence();
        }
    };                                                           
                                                                 
    /// Handles replies to prefetch requests
    static void prefetch_reply_handler(World& world, ProcessID src, void *buf, size_t nbyte) { 
        LongAmArg* arg = (LongAmArg*) buf;
        BufferInputArchive ar(arg->buf,nbyte-sizeof(arg->header));
        RemoteReference<implT> ref;
        ar & ref; 
        implT* impl = ref.get();
        if (debug) print(classname(),"prefetch_got",ref.get(),ref.owner(),"from",src);
        MADNESS_ASSERT(impl);
        ar & *impl->local;
        impl->prefetch_pending = false;
        ref.dec();                    // Matching inc when taking reference in prefetch()
    };                                                           


    /// Private: Takes a reference to an existing world shared ptr via internal pointer

    /// What is this being used for, if anything?
    explicit WorldSharedPtrBase< T > (SharedPtr< WorldSharedPtrImpl<T> >& p)
        : p(p)
    {};


    /// Private: Takes ownership of a reference to a remote object

    /// See RemoteReference for more info.
    explicit WorldSharedPtrBase< T > (const RemoteReference< WorldSharedPtrImpl<T> >& remote_ref) 
        : p(new WorldSharedPtrImpl<T>(remote_ref))
    {};                                       


    /// Returns a remote reference wrapping the local object

    /// Read the RemoteReference documentation to understand the
    /// reponsibilities of the owner of a remote reference.
    ///
    /// You can only get a remote reference to an initialized object.
    ///
    /// It is the LOCAL object that is pointed to by the reference.
    /// This is true even if the local object is referring to a remote
    /// object.  All of the object methods will forward correctly.
    /// However, if you want a reference to the original remote object
    /// (why?) you should call direct_remote_ref() which has to
    /// contact the original object in order to increment the
    /// reference count and get a new remote_ref.  Of course
    /// direct_remote_ref is not yet implemented.
    ///
    /// This method should probably be made public but this part of
    /// the design is still in progress
    RemoteReference< WorldSharedPtrImpl<T> >  new_remote_ref() {
        MADNESS_ASSERT(p);
        return refT(*p->world, p);
    };


public:                                                          
    typedef WorldSharedPtrImpl<T> implT;                       
    typedef RemoteReference<implT> refT;                          
    typedef WorldSharedPtrBase<T> baseT;                         
    typedef WorldSharedPtr<T> objT;                          
    typedef T userT;                                             

    /// Makes an unitialized local pointer
    WorldSharedPtrBase< T > ()                                   
        : p(0)                                                   
    {};                                                          
                                                                 

    /// Takes ownership of a local pointer just like a SharedPtr does

    /// I.e., it must have been allocated with new T(), and this
    /// class assumes the responsibility of eventually deleting it.
    explicit WorldSharedPtrBase< T > (World& world, T* t)        
        : p(new implT(world,t))                                  
    {};                                                          
                                                                 

    /// Takes a reference to a local shared pointer
    explicit WorldSharedPtrBase< T > (World& world, SharedPtr< T >& t)    
        : p(new implT(world,t))                                  
    {};                                                          
                                                                 

    /// Copy constructor is shallow

    /// Copying an unconstructed object throws an exception
    explicit WorldSharedPtrBase< T > (const WorldSharedPtrBase<T>& other) {
        MADNESS_ASSERT(other);
        p = other.p;
    };


    /// Assignment is shallow

    /// Assigning from an unconstructed object throws an exception
    WorldSharedPtrBase< T >& operator=(const WorldSharedPtrBase<T>& other) {
        if (this != &other) {
            MADNESS_ASSERT(other);
            p = other.p;
        }
        return *this;
    };

    /// Set debug flag to new value and return old value

    /// Debugging applies to all instances of this class and
    /// will cause method handlers/wrappers to print hopefully
    /// useful stuff.
    static bool set_debug(bool value) {
        bool status = debug;
        debug = value;
        return status;
    };
                                                                 

    /// True if the pointer is initialized (non-null)
    operator bool() const {
        return p;
    };


    /// Returns rank of owning process, or -1 if not initialized
    inline ProcessID owner() const {                             
        if (p) return p->remote.owner();
        else return -1;
    };                                                           
                                                                 

    /// Returns true if points to local data or if unitialized
    inline bool is_local() const {     
        return p && p->is_local();
    };                                                           
                                                                 

    /// Returns true if local data is available

    /// False means either data not locally available or not
    /// initialized.
    inline bool probe() const {                                  
        return (p && p->local && !p->prefetch_pending);
    };                                                           


    /// Prefetch --- no-op if have local data or prefetch in progress
    inline void prefetch() {                                     
        MADNESS_ASSERT(p);
        if (p->local) return;   // Data already here or prefetch already pending
        MADNESS_ASSERT(p->remote);
        p->prefetch_pending = true;
        p->local = SharedPtr<T>(new T());
        if (debug) madness::print(classname(),"prefetch    ");
        p->world->am.send(p->remote.owner(), 
                          baseT::prefetch_request_handler, 
                          AmArg(p->remote, new_remote_ref()));
        p->world->am.fence();
    };                                                           
                                                                 

    /// Returns pointer to local/cached data, waiting for communication

    /// If the pointer is not initialized, returns 0.  If local data
    /// is available, returns immediately.  Othewise, it issues a
    /// prefetch (if none has been issued already) and waits until the
    /// data is available.
    ///
    /// If probe() returns true, get() will immediately return with data.
    T* get() {
        if (!p) return 0;
        else if (p->local) return p->local.get();

        MADNESS_ASSERT(p->remote);
        if (!p->prefetch_pending) prefetch();
        p->world->await(bind_nullary_mem_fun(this,&WorldSharedPtrBase<T>::probe));
        while(!probe()) World::poll_all();
        return p->local;
    };


    /// Deletes any local cached copy of remote object
    void clear_cache() {                                        
        if (p && p->remote && p->local) p->local = SharedPtr<T>(0);
    };


    /// Flushes cached object to owner, keeping local copy

    /// Does not keep track of any modications ... it always sends.
    void flush_cache() {
        if (p && p->remote && p->local) {
            throw "flush not yet";
        }
    };


    /// Receives via MPI a pointer to a remote object from the specified process

    /// The pointer should have been sent by invoking the remote
    /// objects send_ptr method.
    static WorldSharedPtr<T> recv_ptr(World& world, ProcessID src, Tag tag=71) {  
        refT ref;
        MPI::Request req = world.mpi.Irecv(&ref, sizeof(ref), MPI::BYTE, src, tag); 
        world.await(req);                                     
        return WorldSharedPtr<T>(ref);
    };                                                           
                                                                 
    /// Send a reference to this object to another process

    /// The remote process should call WorldSharedPtr<T>::recv_ptr.
    void send_ptr(ProcessID dest, Tag tag=71) {
        MADNESS_ASSERT(p);
        refT ref = new_remote_ref();
        MPI::Request req = p->world->mpi.Isend(&ref, sizeof(ref), MPI::BYTE, dest, tag);  
        p->world->await(req);                                    
    };                                                           
                                                                 

    /// Convenience routine for exchanging pointers with another process

    /// If the destination and soruce processors are the same (which is
    /// the default) then this exchanges pointers.  But the source and
    /// destination can be different, which is useful, e.g., to pass
    /// pointers around a ring via
    /// \code
    ///    WorldSharedPtr< T > right_obj = obj.send_recv(left,right)
    /// \endcode
    WorldSharedPtr< T > send_recv_ptr(ProcessID dest, ProcessID src=-1, Tag tag=71) {     
        if (src == -1) src = dest;
        refT recv_ref, send_ref = new_remote_ref();
        MPI::Request recv_req = p->world->mpi.Irecv(&recv_ref, sizeof(recv_ref), MPI::BYTE, src, tag); 
        MPI::Request send_req = p->world->mpi.Isend(&send_ref, sizeof(send_ref), MPI::BYTE, dest, tag);  

        p->world->await(recv_req);                                    
        p->world->await(send_req);                                    

        return WorldSharedPtr<T>(recv_ref);
    };                                                           

    /// Broadcast pointer from root to all other processes

    /// This is INHERENTLY NOT SCALABLE since it is unavoidable that
    /// once every process has a reference that each will eventually
    /// send an active message to the owner to decrement the reference
    /// count (unless we also implement a collective global_destructor
    /// on a tree).
    static void broadcast_ptr(World& world, WorldSharedPtr<T>& obj, ProcessID root) {
        if (world.mpi.nproc() == 1) return;

        refT ref;
        if (world.mpi.rank() == root) {
            MADNESS_ASSERT(obj.p);
            MADNESS_ASSERT(obj.p->world->id() == world.id());
            // We need to get nproc-1 references
            for (int i=0; i<world.mpi.nproc()-1; i++) ref = new_remote_ref();
        }
        world.am.broadcast(&ref,sizeof(ref),root);
        if (world.mpi.rank() != root) {
            obj = WorldSharedPtr<T>(ref);
        }
    };

    virtual ~WorldSharedPtrBase() {};
};


#define START_REMOTE_CLASS_INTERFACE(T)                                \
                                                                       \
template <>                                                            \
struct WorldSharedPtrClassNames< T > {                                 \
    static const char* classname() {                                   \
        return "WorldSharedPtrBase<" #T ">";                           \
    };                                                                 \
};                                                                     \
                                                                       \
                                                                       \
template <>                                                            \
class WorldSharedPtr<T> : public WorldSharedPtrBase<T> {               \
public:                                                                \
                                                                       \
    WorldSharedPtr()                                                   \
        : WorldSharedPtrBase< T > () {};                               \
                                                                       \
    explicit WorldSharedPtr(World& world, T* t)                        \
        : WorldSharedPtrBase< T > (world, t) {};                       \
                                                                       \
    explicit WorldSharedPtr(World& world, SharedPtr< T >& t)           \
        :  WorldSharedPtrBase< T > (world, t) {};                      \
                                                                       \
    explicit WorldSharedPtr(const refT& remote_ref)                    \
        :  WorldSharedPtrBase< T > (remote_ref) {};                    \
                                                                       \


#define REMOTE_LONG_METHOD_V1(method, arg1T) \
    static void METHOD_HANDLER(method) (World& world, ProcessID src, void* buf, size_t nbyte) { \
        refT target;                                                   \
        arg1T arg1;                                                    \
        ((LongAmArg*) buf)->unstuff(nbyte, target, arg1);              \
        implT* p = target.get();                                       \
        MADNESS_ASSERT(p);                                             \
        if (p->is_local()) {                                           \
           p->local-> method (arg1);                                   \
        }                                                              \
        else {                                                         \
           MADNESS_ASSERT(p->remote);                                  \
           LongAmArg* out = new LongAmArg;                             \
           size_t nbyte = out->stuff(p->remote, arg1);                 \
           world.am.send_long_managed(p->remote.owner(),               \
                                      objT :: METHOD_HANDLER(method),  \
                                      out,                             \
                                      nbyte);                          \
        }                                                              \
    };                                                                 \
                                                                       \
    inline void method ( arg1T arg1) {                                 \
      MADNESS_ASSERT(p);                                               \
      if (p->is_local()) {                                             \
          p->local-> method (arg1);                                    \
      }                                                                \
      else {                                                           \
          MADNESS_ASSERT(p->remote);                                   \
          LongAmArg* out = new LongAmArg;                              \
          size_t nbyte = out->stuff(p->remote, arg1);                  \
          p->world->am.send_long_managed(p->remote.owner(),            \
                                         objT:: METHOD_HANDLER(method),\
                                         out,                          \
                                         nbyte);                       \
      }                                                                \
    };

#define REMOTE_METHOD_V1(method, arg1T) REMOTE_LONG_METHOD_V1(method, arg1T)

#define END_REMOTE_CLASS_INTERFACE(T) };

#ifdef WORLD_INSTANTIATE_STATIC_TEMPLATES
template <typename T> bool WorldSharedPtrBase<T>::debug = false;
#endif


#endif
