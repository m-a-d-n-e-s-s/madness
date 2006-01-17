/// \file mainpage.h
/// \brief Mainpage for Doxygen generated documentation ... no code.

/** \mainpage Multiresolution Adaptive Numerical Scientific Simulation (MADNESS) Programmers Manual.

    \section intro Introduction

    \section install Installation

\subsection prereqs Prerequisites
<ul>
<li> ISO Standard-conforming C++ compiler
<li> GNU autoconf (version 2.59)
<li> GNU make
<li> Perl (version 5 or higher)
<li> MPI library
<li> BLAS and LAPACK libraries
<li> (optional) Doxygen to produce web-based documentation
<li> (optional) LaTeX to produce printable reference manual
</ul>

\subsection inst Configuration and Installation
\b NOTE: It is recommended that MADNESS is compiled in a directory separate from
the source directory. "In-place" compilation should work too, but such method
is tested less often.

MADNESS package can be installed by following these steps (starting from the top
source directory):
<ol>
<li> Create configure: <tt>aclocal -I lib/autoconf; autoconf</tt>
<li> If compiling in separate object directory, enter that directory, else skip this step
<li> Configure package: <tt>configure <opts></tt>, where <tt><opts></tt> are various
configure command-line options (see below)
<li> Build package: <tt>make</tt> (or <tt>gmake</tt> on some systems)
<li> Build and run tests: <tt>make testbuild; make testrun</tt>
<li> (optional) make documentation: <tt>cd doc; make</tt>
<li> Install binary files and dependencies: <tt>make install</tt>
<li> Install libraries and header files: <tt>make install_devel</tt>
<li> (optional) Return package to its initial state: <tt>make distclean</tt>
</ol>
Configure script can take many command-line options, which can be optained by <tt>configure --help</tt>.

\subsection configopts Common configure problems and solutions

<ul>

<li> Most installations will require that the following command-line arguments
given to configure: <tt>--with-cxx --with-cxx-optflags --with-cppflags --with-libdirs 
--with-libs --with-f77symbol</tt>.

<li> When you use configure options <tt>--with-blas</tt> and <tt>--with-lapack</tt>,
the given BLAS and LAPACK libraries are not tested for usability. It is thus recommended
to use <tt>--with-libs</tt> and <tt>--with-libdirs</tt> options to specify location of
BLAS and LAPACK libraries.

<li> By default MADNESS builds static libraries. If you would like to build shared libraries,
you must enable libtool as well: <tt>--enable-libtool --enable-shared</tt>.
Shared library build may also require that linker is given to configure explicitly
via <tt>--with-ld</tt> option. In any case, C++ compiler is always used to link executables.

</ul>

\subsection configtodo To-Do list for build system

<ul>
  <li> replace remaining preprocessor switches in the code by configure tests
  <li> automate extraction of type information (sizeof, etc.)
  <li> strengthen include guards
  <li> test cross-compilation
  <li> test various template instantiation models
  <li> test enable-shared on more architectures
  <li> add README, INSTALL, LICENSE, WARRANTY, COPYRIGHT, official-looking headers
  <li> 
</ul>

    \section doc Documentation

    \section test Testing

    \section style Style guide

Document near the top of every file the name and purpose of the file
using the Doxygen commands
\verbatim
/// \file filename
/// \brief Brief description of contents

/// optional detailed description
\endverbatim

Document using Doxygen the purpose of all procedures, even static or
private ones, providing both brief and detailed descriptions.  If this
documentation does not describe all input and output arguments,
separately document the arguments where they are declared.

Place documentation, if possible, in the file containing the
implementation, rather than the header file defining the interface.
This is to increase the probability that you will keep the two
consistent.

All MADNESS related global names should be put into the \c madness
namespace, or subspaces thereof.  Once we get a lot of
application-like functionality, rather than MADNESS-related stuff, we
may have to add additional namespaces.  

Note that namespace definition nest --- so do not include header files
within a namespace unless you explicitly want to nest all definitions.

Do not put "using namespace ..." or "using space::item" into header
files.  Instead, in the header file, correctly resolve every use of an
external name into its namespace.  Otherwise, any source file including
your header immediately pollutes its own namespace in unexpected ways.

In source (not header) files do \em not use entire name spaces, 
not even the \c madness namespace.  I.e., do not use stuff
like this
\code
using namespace std;
using namespace madness;
\endcode
since it ends up being too hard to manage:
 - What are you actually using from each namespace?

Instead, either correctly resolve every use, or, place immediately
below the relevant \c #include \c using directives for each item used.
E.g.,
\code
#include <iostream>
using std::cout;
using std::endl;

#include <mra/mra.h>
using madness::Function;
\endcode
Your code is forced to document exactly what it is using from 
where, which is "a good thing."
The only exception to this should be importing third-party or
legacy code.

    \section pyinter Python-C++ interface

    \section notes Notes --- to-do list and random design thoughts

Don't take all of the entries here too seriously ... some are just 
neurons firing randomly.

\subsection test Testing

1) run this mess on the cray sooner rather than later ... boost is alleged to compile

1a) template specializations (not instantiations) must be declared in header file

1b) -h conform

1c) -h instantiate=used should be put on files that others will link against
    but not for all source

1d) how big are the basic types

2) run this mess on the ibm ... boost won't compile with xlC ... must use gcc, so 
   perhaps is less of a big deal

3) coverage analysis ... seem to need a very recent gcc toolset to get this
   to work ok.  need to introduce makefile targets to do this

4) quantify performance.

\subsection oldmranotes [Relating to OLD hash based version] mra notes

1) functions need  consistent error handling

2) exceptions from libraries incl. tensor need to
to be caught locally, and failing that globally.

3) confusion between hash3d(constructor args) and
   hash3d operator() ... is this now resolved?

5) should we make hash interface consistent with that of the
   stl map ... so you can drop in replace them ...
   not completely possible esp. with plan for shmem

6) specify behaviour of get_ptr on failure ... does it
   change ptr?  It should not.

7) hash3d should probably bite it ... just use hash
   interface since it provides uniformity and hash3d
   does not really simplify and does not optimize at all.

   ... no ... hash3d needs to provide a deep copy
   constructor and assignement ... so needs to copy
   the tensors ... however should probably remove
   inessential methods to simplify the interface.

8) need to think thru function const members and 
   template parameters.

9) max_refine_level should be a constant due to dimensioning of coeffs

10) on startup should run self check of quadrature, etc., for
    paranoia's sake.

11) defaulting of parameters is done at compile time.  a small test of
    the current mechanism of specifying the default for a parameter as
    a structure element that can be changed at runtime seemed to work.
    But, there are probably some gotchas lurking.

12) transform3 optimization ... succesful but needs to be applied to
    the quadrature step to speed up projection as well as
    compress/reconstruct.  Also, need to solve problem with breaking
    constness of (un)filter and also thread safety.

\subsection tensornotes Tensor notes

5) document what operations are needed to be specialized for 
   new types to be added to tensor
   - random number generator
   - tensor::norm
   - type information

6) on assignments of scalars for fill, sometimes an explicit
   cast is necessary, e.g.,

   Tensor<float_complex> t(1,2,3)
   t = 0;               // fails to compile confusing null pointer and zero.

   t = float_complex(0) // OK

   long x = 0
   t = x                // OK

   !! This should be fixed now that I have eliminated the non-PC indexing but
   I have not checked this.

7) document for user integer values don't auto convert to complex

8) overlapping slices ... behavior undefined ?  OK.

9) recall why templated methods must be defined inline ... essentially
   because of instantiation problem for general case ???  Or could I
   just not figute out how to instantiate all of these methods in 
   tensor.cc and tesnoriter.cc?

11) need to be able to wrap a pointer to memory owned by someone else,
    and optionally free it in the usual manner.

    Wrapping pointers with ownership might enable to optimize
    away some copy constructors due to the return constructor
    optimization.

15) documentation ... coming slowly

18) complex conjugate, real, imag for complex types, initialize
    from separate real & imag tensors.  clearly, complex types
    need a lot more thought ... defer this until needed.

19) with few/no exceptions, class methods (except for the arithmetic
    operators) should operate inplace or return a view.  functions
    operating on tensors can return new tensors ... is this now done?

21) change traits templates to use boost traits classes ?  ...
    no clear benefit ?

22) increase consistency with stl containers including vector and
    iterators ... the foward iterator for a single tensor could
    be wrapped to be compatible.

23) trace should probably be renamed dot ... but should it appropriately
    conjugate complex arguments???

25) converters to/from multiarrays and numarrays ... crucial to
    reuse functionality and import data

27) boost has a slice class ... with stride as the second arg 
    so it is not defaultable and there is no mechanism to 
    handle dimensions of unkown size (i.e., no convention to 
    refer to the end).  boost also provides range.

28) The one past end for an iterator end is a pretty wide convention,
    and it is appealing to type Slice(0,k) rather than Slice(0,k-1)
    all the time.  So, should we change slice back to python convention?

    Problem is with specifying the end.  We choose -1 to mean the last
    element included is dim-1.  A python slice with -1 would actually
    give dim-2 as the last element included.  Python -1 as a scalar
    index gives dim-1.  In Python, the end of a dimension in a slice
    is specified by leaving the end field blank.  We don't have that
    luxury.  So, if we want Python interface compatibility, we will
    need to introduce a symbol and special value to indicate the end
    of a dimension of unknown size.

    Another issue with the Python convention is that if iterating
    foward, -1 means ndim-2, but if iterating backward it is probably
    intended to mean one past the beginnning end, i.e., 0.

    All in all, as long as the conventions are clearly documented,
    I'm inclined to stick with the current approach until more 
    experience is accumulated from external users.  It is not too
    hard to refactor the code to accomodate this change.

    Pehaps to reduce confusion we should rename Slice as Patch.

29) need to be more extensive in our use of STL containers.  all of
    the design and range checking etc. is a good thing ... though
    vector does not range check most indexing operations.

30) are we missing some tensor iterators unoptimized?

31) can we make tensor::id a constant .. it should be.

32) all integer arguments should be longs.  this avoids any problems
    on 64-bit machines in adressing large memories, and eliminates
    any confusion about what is expected.

33) clean up behavior of tensors allocated with default constructor
    ... using them should not cause crashes.  Document methods that
    unavoidably might do this (e.g, indexing) without introducing
    undesirable runtime overhead.

34) header file naming convention ... underscores or not betwen words?

35) how on earth does one easily and efficiently make and initialize
    an STL vector?  had to introduce vector_factory to make this easy.

36) extracted the mxm routines from tensor.cc and enabled calling of
    the Goto blas ... but it was slower ... disabled the blas and 
    it was still slower ... put the blas back into tensor.cc and
    it was still slower ... put the static qualifier on the mxm
    routine names and back to the original speed ... IPO, I guess.

\subsection mra tree code coverage
Missing coverage for
 - default constructor ?
 - reconstructing redundant from noncompressed tree
 - trace
 - inner
 - diff
 - operator+=
 - operator-=
 - _init default constructor for zero compressed function
 - _get_scaling_coeffs
 - _norm_tree


\subsection tensor code coverage

Major missing parts are
 - Tensor init
   - handling of failure in new memory allocation
   - p=0 assigment ... when is this triggered?
 - operator+(tensor)
 - operator-(tensor)
 - operator+(scalar)
 - operator-() unary negation
 - operator-(scalar)
 - operator+=(scalar)
 - operator(vector) ... general indexing including bounds checking
 - operator(slice)
 - operator(long,slice)
 - operator(slice,slice)
 - operator(slice,slice,long)
 - operator(long,long,slice), and variants
 - operator(slice,slice,slice,slice) and 5d and 6d
 - reshape(vector)
 - reshape(long,long,long) and 4d, 5d variants
 - splitdim(long,long,long) 
 - fusedim(long)
 - screen(double)
 - sumsq()
 - product()
 - min() - index generation ?
 - max() - index generation ?
 - absmin() - ditto
 - absmax() - ditto
 - unaryop
 - complex min, max
 - copy ?????
 - transform ???
 - abs()
 - arg()
 - real()
 - imag()
 - conj()

*/
