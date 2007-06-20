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
   transparent remote method invocation, an owner-computes task model,
   and optional user control over placement/distribution.
 - Distributed objects that can be globally addressed.
 - Futures (results of unevaluated expressions) for composition of latency tolerant
   algorithms and expression of dependencies between tasks.
 - Globally accessible task queues in each process which 
   can be used individually or collectively to provide a single global
   task queue.
 - Work stealing for dynamic load balancing (coming v. soon).
 - Facile management of computations on processor sub-groups.
 - Integration with MPI
 - Optional integration with Global Arrays (J. Nieplocha, 
   http://www.emsl.pnl.gov/docs/global).
 - Active messages to items in a container, distributed objects, 
   and processes (currently via polling on MPI but ports to GPC,
   GASNET, and vendor provided libraries are planned).
 - User-space threading or continuations to eliminate the need for explicit
   continuations or use of Futures when composing latency tolerant algorithms (coming).
 - Kernel-space threading for use of multi-core processors (coming).

\section Motivations and attributions

There were several motivations for developing this environment.
 -# The rapid evolution of machines from hundreds (pre-2000), to
    millions (post-2008) of processors demonstrates the need to abandon
    process-centric models of computation and move to paradigms that
    virtualize or even hide the concept of a process.  
    The success of applications using the 
    Charm++ environment to scale raplidy to 30+K processes and the enormous effort
    required to scale most process-centric applications are the central examples.
 -# The arrival of multi-core processes and the associated needs of
    expressing much more concurrency and adopting techniques for
    latency hiding motivate the use of light weight work queues to
    capture much more concurrency and the use of futures for
    latency hiding.
 -# The complexity of composing irregular applications in partitioned, global-address space
    (PGAS) models using only MPI and/or one-sided memory access (GA, UPC, SHMEM, co-Array) 
    motivates the use of an object-centric active-message or remote method invocation (RMI) model 
    so that computation may be moved to the data with the same ease as 
    which data can be moved.  This greatly simplifies the task of maintaining
    and using distributed data structures.
 -# Interoperability with existing programming models to leverage existing
    functionality and to provide an evolutionary path forward.

The two main early influences for this work were Cilk (Kuszmaul,
http://supertech.csail.mit.edu/cilk) and Charm++ (Kale,
http://charm.cs.uiuc.edu).  Subsequently, ACE (Schmidt,
http://www.cs.wustl.edu/~schmidt/ACE.html), STAPL (Rauchwerger and
Amato, http://parasol.tamu.edu/groups/rwergergroup/research/stapl), and
the HPCS language projects (X10,
http://domino.research.ibm.com/comm/research_projects.nsf/pages/x10.index.html
; Chapel, http://chapel.cs.washington.edu ; Fortress, http://fortress.sunsource.net )
and the amazingly talented teams and individuals developing these.


\section Introduction

The entire parallel environment is encapsulated in an instance of the
class World which is instantiated by wrapping an MPI communicator.
Multiple worlds may exist, overlap, or be dynamically created and
destroyed.

The World class has members
 - mpi - an instance of WorldMPIInterface,
 - am - an instance of WorldAMInterface,
 - taskq - an instance of WorldTaskQueue, and
 - ga - an instance of WorldGAInterface, and
 - others coming.

Distributed containers (currently associative arrays or hash tables)
and distributed objects may be constructed from a world instance.

The recommended approaches to develop scalable and latency tolerant
parallel algorithms are either object- or task-centric decompositions
rather than the process-centric approach usually forced upon MPI
applications.  The object-centric approach uses distributed containers
(or distributed objects) to store application data.  Computation is
expressed by sending tasks or messages to objects, using the task
queue to automatically manage dependencies expressed via futures.
Placement of data and scheduling/placement of computation can be
delgated to the container and task queue, unless there are spefic
performance concerns in which case the application can have full
knowledge and control of these.

Items in a container may be accessed largely as if in a standard STL
container, but instead of returning an iterator, accessors instead
return a Future<iterator>. A future is a container for the result of a
possibly unevaluated expression.  In the case of an accessor, if the
requested item is local then the result is immediately
available. However, if the item is remote, it may take some time
before the data is made available locally.  You could immediately try
to use the future, which would work but with the downside of
internally waiting for all of the communication to occur.  Much better
is to keep on working and only use the future when it is ready.


Aside:
  - To avoid a potentially unbounded nested invocation
    of tasks which could overflow the stack and also be the source
    of live/deadlocks, new tasks
    are not presently started while blocking for communication.
    This will be relaxed in the near future which will reduce
    the negative impact of blocking for an unready future as long
    as there is work to perform in the task queue.
  - Once fibers or user-space threads are integrated, multiple
    tasks will always be scheduled and blocking will merely
    schedule the next fiber.

By far the best way to compute with futures is to pass them as
arguments to a new task.  Once the futures are ready, the task will be
automatically scheduled for execution.  Tasks that produce a
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

The \c World.am member provides inter-process active message functionality, which is
the foundation on which everything else is built.  We do not recommend
that applications make routine or direct use of inter-process active messages.
Instead, try to compose applications using messaging 
to/between items in distributed containers and the local
task queue(s).

The \c World.mpi member is the preferred way to use MPI since it has a growing
amount of instrumentation and debugging capability, though MPI
routines may be called directly if necessary.  However, MPI is again a
low-level model and we do not encourage its direct use.  It is there
since it is the portable standard for communication and to facilitate
integration with legacy applications.

The \c World.ga member provides access to the capabilities of the
Global Array library (this is still unfolding).

The execution model is sequentially consistent.  That is, 
from the perspective of a single thread of execution, operations 
on the same local/remote object behave as if executed sequentially
in the same order as programmed.   This means that performing
a read after a write/modify returns the modified value, as expected.
Such behavior applies only to the view of a single thread --- 
the execution of multiple threads and active messages from different
threads may be interleaved arbitrarily.

Creating, executing, and reaping a local, null task with 
no arguments or results presently takes about 350ns (Centos 4, 3GHz 
Core2, Pathscale 3.0 compiler, -Ofast).  The time
is dominated by \c new and and \c delete of the
task structure, and as such is unlikely to get any faster
except by the application caching and reusing the task structures.   
Creating and then executing a chain of
dependent tasks with the result of one task fed as the argument
of the next task (i.e., the input argument is an unevaluated future 
which is assigned by the next task) requires about 2000ns per
task, which we believe can be redcued
to about 1us (3 GHz Core2).  

Creating a remote task adds the
overhead of interprocess communication which is on the scale of 1-3us
(Cray XT).  Note that this is not the actual wall-time latency since
everything is presently performed using asynchronous messaging and
polling via MPI.  The wall-time latency, which is largely irrelevant
to the application if it has expressed enough parallelism, is mostly
determined by the polling interval which is dynamically adjusted
depending upon the amount of local work available to reduce the
overhead from polling.  We can improve the runtime software through better
agregation of messages and use of deeper message queues to reduce the
overhead of remote task creation to essentially that of a local task.

Thus, circa 1us defines the ganularity above which it is worth
considering encapsulating work (c.f., Hockney's n1/2).  However, this
is just considering the balance between overhead incurred v.s. useful
work performed.  The automatic scheduling of tasks dependent upon
future arguments confers many benefits, including
 - hiding the wall-time latency of remote data access,
 - removing from the programmer the burden of correct scheduling
   of dependent tasks, 
 - expressing all parallelism at all scales of the algorithm 
   for facile scaling to heavily multi-core architectures and
   massively parallel computers, and
 - virtualizing the system resources for maximum 
   future portability and scalability.

Available memory limits the number of tasks that can be generated
before any are consumed.  In addition to application specific data,
each task consumes circa 64 bytes on a 64-bit computer.  Thus, a few
hundred thousand outstanding tasks per processor are eminently
feasible even on the IBM BG/L.  Rather than making the application
entirely responsible for throttling it's own task production (which it
can), if the system exceeds more than a user-settable number of
outstanding tasks, it starts to run ready tasks before accepting new
tasks.  The success of this strategy presupposes that there are ready
tasks and that these tasks on average produce less than one new task
with unsatisfied dependencies per task run.  Ultimately, similar to
the Cilk execution model, safe algorithms (in the same sense as safe
MPI programs) must express tasks so that dependencies can be satisfied
without unreasonable expectation of buffering.

In a multiscale approach to parallelism, coarse gain tasks are
first enqueued, and these generate finer-grain tasks, which 
in turn generate finer and finer grain work.   [Expand this discussion
and include examples along with work stealing discussion]

Discussion points to add
 -# Why arguments to tasks and AM via DC or taskQ are passed
    by value or by const-ref (for remote operations this
    should be clear; for local operations it is to enable
    tasks to be stealable).  Is there a way to circumvent it? Pointers.
 -# Virtualization of other resources
 -# Task stealing
 -# Controlling distribution in containers
 -# Caching in containers
 -# Computing with continuations (user space fibers)

\section Distributed Containers (WorldContainer)

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
Please refer to world/archive/archive.h and documentation therein for
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


\section Distributed Objects (WorldObject)

Distributed objects (WorldObject) provide all of the communication
and other resources necessary to build new distributed capabilities.
The distributed container class (WorldContainer) actually inherits 
most of its functionality from the WorldObject.

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
WORLD_INSTANTIATE_STATIC_TEMPLATES is defined.  In one of your source
(not header) files,
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

\subsection Reference counting and arguments to STL agorithms

It is not specified how arguments are passed to STL algorithms so
they are free to pass by value (a concrete example of this "feature"
is provided by gcc4.2).  For reference counted types (e.g.,
SharedPtr), this causes the count to be incremented (if passed by
value) or not incremented (if passed by reference) and hence you
cannot reliably use an STL algorithm to look for stuff whose reference
count has a specific value (e.g., one).  

E.g., the following does not work reliably due to the standard
being incomplete (i.e., borked).

\code
  struct refcnt_is_one {
      bool operator()(const SharedPtr<DeferredCleanupInterface>& p) const {
        return p.use_count() == 1;
      };
  };

  std::remove_if(deferred.begin(),deferred.end(),refcnt_is_one());
\endcode
Instead, you need the explicit loop
\code
  for (std::list< SharedPtr<DeferredCleanupInterface> >::iterator it = deferred.begin(); 
      it != deferred.end();) {
      if (it->use_count() == 1) 
          it = deferred.erase(it);
      else 
          ++it;
  }
\endcode



\section Development to-do list

 - Test multiple worlds
 - Cache semantics and testing for DC
 - Prefetch marker for DC cache
 - Verify reference counting for DC cache
 - Forwarding for DC
 - Integration with user-space thread/fiber scheduling
 - Performance profiling with Tau
 - ASM clock for PPC and BGL
 - Test with xlCC

*/
