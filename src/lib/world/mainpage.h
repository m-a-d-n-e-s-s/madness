/// \file mainpage.h
/// \brief MADNESS Parallel Programming Environment

/** \mainpage

\section Overview

The MADNESS parallel programming environment combines several
successful elements from other models and aims to provide a
rich and scalable framework for massively parallel computing while
seamlessly integrating with legacy applications and libraries.
It includes
 - Distributed sparse containers with one-sided access to items,
   transparent remote method invocation, owner-computes task model,
   and optional user control over placement/distribution.
 - Futures to facilitate composition of latency tolerant
   algorithms and expression of depenencies between tasks.
 - Active messages, currently via polling on MPI but ports to GPC, GASNET,
   and vendor provided libraries are planned.
 - User-space threading or continuations to eliminate the need for explicit
   continuations or use of Futures when composing 
   latency tolerant algorithms (coming).
 - Kernel-space threading for use of multi-core processors (coming).
 - Globally accessible task queues in each process which 
   can be used individually or collectively to provide a single global
   task queue.
 - Work stealing for dynamic load balancing (coming v. soon).
 - Facile management of computations on processor sub-groups.
 - Integration with MPI
 - Optional integration with Global Arrays (J. Nieplocha, 
   http://www.emsl.pnl.gov/docs/global).

\section Introduction

All information about and all functionality of the parallel
environment is accessed through an instance of the class
World which is instantiated by wrapping an MPI communicator.
Multiple worlds may exist, overlap, or be dynamically created
and destroyed.

The World class has members
 - mpi - an instance of WorldMPIInterface,
 - am - an instance of WorldAMInterface,
 - taskq - an instance of WorldTaskQueue, and
 - ga - an instance of WorldGAInterface, and
 - others coming.

Distributed containers (currently associative arrays or hash tables)
may be constructed from a world instance.

The recommended approach to develop scalable and latency tolerant
parallel algorithms is to place shared data into one or more
distributed containers and to express computation as tasks with
dependencies managed via futures.  Placement of data and
scheduling/placement of computation can be delgated to the container
and task queue, unless there are spefic performance concerns in which
case the application can have full knowledge and control of these.

Items in a container may be accessed pretty much as if in a standard
STL container, but what is returned is a Future<item_type>, which is
a container for the result of a possibly unevaluated expression.  If the
requested item is local, the result is immediately available.
However, if the item is remote, it may take some time before the data
is made available locally.  You could immediately try to use the
future, which would work but with the downside of having to wait
for all of the communication to occur.  Much better is to keep on
computing with available data and only use the future when it is
ready. 

Aside:
  - To avoid a potentially unbounded nested invocation
    of tasks which could overflow the stack, new tasks
    are not presently started while blocking for communication.
    This will be relaxed in the near future which will reduce
    the negative impact of blocking for an unready future.
  - Once fibers or user-space threads are integrated, multiple
    tasks will always be scheduled and blocking will merely
    schedule the next fiber.

By far the best way to compute with futures is to give them as
arguments to a new task.  Once the futures are ready, the task will be
automatically scheduled as ready for execution.  Tasks that produce a
result also return it as a future, so this same mechanism may be used
to express dependencies between tasks.

Thus, a very natural expression of a parallel algorithm is as a
sequence of dependent tasks.  For example, in MADNESS many of the
algorithms working on distributed, multidimension trees start with
just a single task working on the root of the tree, with all other
processes waiting for something to do.  That one task starts
recursively (depth or breadth first) traversing the tree and
generating new tasks for each node.  These in turn generate more tasks
on their sub-trees.

The \c World.am member provides active message functionality, which is
the foundation on which everything else is built.  We do not recommend
that applications make routine or direct use of active messages.
Instead, try to compose applications using distributed containers, the
task queue(s), and messaging between objects in containers.

The \c World.mpi member is the preferred way to use MPI since it has a growing
amount of instrumentation and debugging capability, though MPI
routines may be called directly if necessary.  However, MPI is again a
low-level model and we do not encourage its direct use.  It is there
since it is the portable standard for communication and to facilitate
integration with legacy applications.

The \c World.ga member provides access to the capabilities of the
Global Array library (this is still being developed).

Discussion points to add
 -# Sequential memory consistency (read/write ordering)
 -# Sequential execution consistency (bounded buffer problem)
 -# Throttling task production (bounded buffer problem)
 -# Multiscale approach to task production
 -# Virtualization of resources
 -# Task stealing
 -# Controlling distribution in containers
 -# Caching in containers
 -# Computing with continuations (user space fibers)
 -# Why arguments to tasks and AM via DC or taskQ are passed
    by value or by const-ref (for remote operations this
    should be clear; for local operations it is to enable
    tasks to be stealable).  Is there a way to circumvent it?

\section Distributed Containers

The only currently provided containers are associative arrays or maps
that are almost directly equivalent to the STL map or the GNU
hash_map.  Indeed, the implementation can use either of these for the
local storage, though the GNU hash_map is to be preferred for
performance reasons and is the only one discussed here.

A map generalizes the concept of an array (which maps an integer index
in a dense range to a value) by mapping an arbitrary key to a value.
This is a very natural, general and efficient mechanism for storing
sparse data structures.  The distribution of items in the container
between processes is based upon a function which maps the key
to a process.  There is a default mapping which is essentially 
a pseudo-random uniform mapping, but the user can provide their own
(possibly data-dependent) operator to control the distribution.  

Although it will almost always be the case that all processes agree on
the mapping of a key to a process, this does not have to be the case
since the implementation supports forwarding
of remote requests.  \em NOT YET COMPLETED ... but it will be cool
when it is finished!

The keys and values associated with containers must be serializble
by the MADNESS archive mechanism.
Please refer to world/archive.h and documentation therein for
information about this.  In addition, the keys must support
 - testing for equality, either by overloading \c == or by
   specializing \c std::equal_to<key_type>, and
 - computing a hash value by invoking \c madness::hash(key),
   which can be done either by providing the member
   function with signature
\code
   hashT hash() const;
\endcode
   or by specializing \c madness::Hash<key_type>.

\c hashT is presently an unsigned 32-bit integer.  MADNESS provides
hash operations for all fundamental types, and variable and fixed
dimension arrays of the same.  Since having a good hash is important,
we are using Bob Jenkin's "lookup v3" hash from
http://www.burtleburtle.net/bob/c/lookup3.c.

Here is an example of a key that might be used in an octtree.
\code
   struct Key {
       typedef unsigned long ulong;
       ulong n, i, j, k;
       hashT hashval;

       Key() {};

       // Precompute the hash function for speed
       Key(ulong n, ulong i, ulong j, ulong k)
           : n(n), i(i), j(j), k(k), hashval(madness::hash(&this->n,4,0)) {};

       hashT hash() const {
           return hashval;
       };

       template <typename Archive>
       void serialize(const Archive& ar) {
           ar & n & i & j & k & hashval;
       }

       bool operator==(const Key& b) const {
           // Different keys will probably have a different hash
           return hashval==b.hashval && n==b.n && i==b.i && j==b.j && k==b.k;
       };
   };
\endcode

To be added
 - discussion of chaining hashes using initval optional argument
 - discussion of overriding the distribution across processes


\section Compilation and linking

\subsection C-preprocessor predefined macros

A list of C-preprocessor macros to be used for managing machine
dependencies that are either defined by the system or by
the MADNESS build process.

 - \c __GNUG__ --- The GNU processor and C++
 - \c X8632 --- A 32-bit x86 CPU.
 - \c X8664 --- A 64-bit x86 CPU.
 - \c UINT64_T --- The type for a 64-bit unsigned integer which is used to
    typedef uint64_t.  The macro is not defined
   if uint64_t is already a valid type.
 - \c _CRAY - Any Cray system (though currently we only support XT3/4).


\subsection Static data, etc., for templated classes

Several of the templated classes (currently just the
DistributedContainer, Future and RemoteReference classes) have static
data or helper functions associated with them.  These must be defined
in one and only one file.  To facilitate this definition, the
necessary templates have been wrapped in C-preprocessor conditional
block so that they are only enabled if \c
WORLD_INSTANTIATE_STATIC_TEMPLATES is defined.  In one of your files,
define this macro \em before including \c world.h, and then
instantiate the templates that you are using.


\section Gotchas

\subsection Futures and STL vectors (e.g., \c vectors<Future<int>> )

A common misconception is that STL containers initialize their
contents by \invoking the default constructor of each item in
the container since we are told that the items must be default
constructable.  But this is \em incorrect.  The items are initialized
by invoking the copy constructor for each element on a \em single
object made with the default constructor.   For futures this 
is a very bad problem.  For instance,
\code
   vector< Future<double> > v(3);
\endcode
is equivalent to the following with an array of three elements
\code
   Future<double> junk;
   Future<double> v[3] = {junk,junk,junk};
\endcode
Since the Future copy constructor is by necessity shallow, each
element of \c v ends up referring to the future implementation that
underlies \c junk.  When you assign to an element of \c v, you'll also
be assigning to junk.  But since futures are single assignment
variables, you can only do that once.  Hence, when you assign a
second element of \c v you'll get a runtime exception.

The fix (other than using arrays) is to initialize STL vectors and
other containers from the special element returned by
\c Future<T>::default_initializer() which if passed into the copy
constructor will cause it to behave just like the default contructor.
Thus, the following code is what you actually need to use an STL
vector of futures
\code
   vector< Future<double> > v(3,Future<double>::default_initializer());
\endcode
which sucks, so we provide the factory function
\code
   template <typename T>
   vector< Future<T> > future_vector_factory(std::size_t n);
\endcode
which enables you to write
\code
   vector< Future<double> > v = future_vector_factory<double>(3);
\endcode
which merely blows instead of sucking.

\section Development to-do list

 - Test multiple worlds
 - Cache semantics and testing for DC
 - Prefetch marker for DC cache
 - Verify reference counting for DC cache
 - Forwarding for DC
 - Cache clean for DC
 - Integration with user-space thread/fiber scheduling
 - Performance profiling with Tau
 - ASM clock for PPC and BGL
 - Test with pathscale compiler
 - Test with xlCC
 - What's happening with temporary DC's and deferred destruction?

*/
