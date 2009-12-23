/// \file mainpage.h
/// \brief MADNESS 

///  \defgroup configuration MADNESS installation and configuration
///  \defgroup libraries MADNESS libraries

///  \defgroup parallel_runtime Parallel programming environment
///  \ingroup libraries

///  \defgroup mra Multiresolution analaysis
///  \ingroup libraries

///  \defgroup funcplot Function plotting routines
///  \ingroup mra

///  \defgroup tensor Tensors or multidimension arrays
///  \ingroup libraries

///  \defgroup linalg Linear algebra (interface to LAPACK)
///  \ingroup libraries

///  \defgroup solvers Iterative solvers for linear/non-linear equations and optimizers
///  \ingroup libraries

///  \defgroup misc Miscellany
///  \ingroup libraries

///  \defgroup applications MADNESS applications

///  \defgroup examples Examples 
///  \ingroup applications


/*!

\mainpage


\section overview Overview

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
 - Work stealing for dynamic load balancing (prototype being tested)
 - Facile management of computations on processor sub-groups.
 - Integration with MPI
 - Optional integration with Global Arrays (J. Nieplocha, 
   http://www.emsl.pnl.gov/docs/global).
 - Active messages to items in a container, distributed objects, 
   and processes.
 - Efficient use of multicore processors using pthreads.

\subsection used_by Software used by MADNESS

Our deep gratitude to these other projects whosse software
we are using within MADNESS

 - TinyXML - http://sourceforge.net/projects/tinyxml

 - CFFT - http://www.librow.com/articles/article-10/appendix-a-2

 - muParser - http://muparser.sourceforge.net

 - libxc - http://www.tddft.org/programs/octopus/wiki/index.php/Libxc

\section motivations Motivations and attributions for the parallel runtime

There were several motivations for developing this environment.
 -  The rapid evolution of machines from hundreds (pre-2000), to
    millions (post-2008) of processors demonstrates the need to abandon
    process-centric models of computation and move to paradigms that
    virtualize or even hide the concept of a process.  
    The success of applications using the 
    Charm++ environment to scale rapidly to 30+K processes and the enormous effort
    required to scale most process-centric applications are the central examples.
 -  The arrival of multi-core processes and the associated needs of
    expressing much more concurrency and adopting techniques for
    latency hiding motivate the use of light weight work queues to
    capture much more concurrency and the use of futures for
    latency hiding.
 -  The complexity of composing irregular applications in partitioned, global-address space
    (PGAS) models using only MPI and/or one-sided memory access (GA, UPC, SHMEM, co-Array) 
    motivates the use of an object-centric active-message or remote method invocation (RMI) model 
    so that computation may be moved to the data with the same ease as 
    which data can be moved.  This greatly simplifies the task of maintaining
    and using distributed data structures.
 -  Interoperability with existing programming models to leverage existing
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

\section intro_to_runtime Introduction to the parallel runtime

The entire parallel environment is encapsulated in an instance of the
class \c World which is instantiated by wrapping an MPI communicator.
Multiple worlds may exist, overlap, or be dynamically created and
destroyed.  Distributed containers (currently associative arrays or
hash tables) and distributed objects may be constructed from a world
instance.

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

The \c World.gop member provides global operations that are internally
non-blocking, enabling the invoking thread to continue working.

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
task, which we believe can be redcued to about 1us (3 GHz Core2).  

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

\section dist_cont Distributed Containers (WorldContainer)

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

\section dist_obj Distributed Objects (WorldObject)

Distributed objects (WorldObject) provide all of the communication
and other resources necessary to build new distributed capabilities.
The distributed container class (WorldContainer) actually inherits 
most of its functionality from the WorldObject.


\subsection static_data Static data, etc., for templated classes

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


*/
