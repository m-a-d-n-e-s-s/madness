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

  
/*

Copyright (c) 2006, R.J. Harrison, UT/ORNL
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright 
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the names of Oak Ridge National Laboratory and the
      University of Tennessee, nor the names of its contributors may
      be used to endorse or promote products derived from this
      software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef FIBER_H
#define FIBER_H

/* A fiber is a non-preemptively scheduled execution context, each
   fiber having an independent instruction stream and stack, but
   otherwise sharing all of the other state of the containing thread.
   Essentially, a fiber is a user-space thread with FIFO scheduling.

   There are only two important routines for most uses  
   -  fiber_create() is used to create fibers, and 
   -  fiber_yield() is used to yield control to the next fiber in the 
      FIFO queue.  

   The remaining routines give the number of currently active fibers,
   provide a hook for a function to be invoked on every context
   switch, and set/get the stack size for each fiber.

   New fibers start execution by invoking a user-specified function
   with prototype

   void f(void *data)
  
   where data is a user-specified pointer (which may be null).  The
   fiber terminates by returning from f().

   Control can only be passed to another fiber by the fiber either
   terminating (return from f()) or by yielding control via a call to
   fiber_yield().  Execution will be resumed when the fiber is next
   scheduled by return from fiber_yield().  Note that the original
   or main fiber of execution is included in the scheduling (and 
   if it returns or exits the entire thread with all of the fibers
   will be terminated).

   The non-preemptive, user-space scheduling implies that there is NO
   CONCURRENCY or overlap in the execution of fibers in the same
   kernel-level thread.  Since execution of fibers is strictly
   interleaved, you do not need any mutex between fibers in the same
   thread.  There is also no danger of multiple fibers in the same
   thread entering non-reentrant library code (unless it calls back to
   your code).

   The code is not yet thread safe, mostly because if you've got
   threads you've probably already got fibers, and we don't have
   threads yet on the Cray XT3.  But it is not hard to do this.  Just
   need to put the already identified global/static structures into
   thread-private data.  */


/* Create a new fiber which will invoke f(data).  Returns true on
   success, false on failure */
int fiber_create(void (*f)(void *), void *data);


/* Yield to the next fiber in the FIFO queue.  Execution will resume
   by returning from fiber_yield() */
void fiber_yield();


/* Returns the number of currently active fibers created by
   fiber_create() */
int fiber_count();


/* Specifies a function and associated data to be invoked on every
   fiber context switch.  I.e., hook(data) is invoked every time a fiber
   terminates or yields.  The function has prototype

   void hook(void *data)

   and MUST NOT directly or indirectly call fiber_yield(), or 
   you'll have an infinite loop.

   Natural uses for the hook are to gather statistics or to implement
   a polling mechanism with the granulatiry of fiber self-scheduling.

   To unset the hook specify a null pointer for the function, i.e.,

   fiber_set_switch_hook(0,0) */
void fiber_set_switch_hook(void (*f)(void *), void *data);


/* Get the default stack size in bytes */
size_t fiber_get_stacksize();

/* Set the default stack size in bytes.  Default is 67kb with a
   minimum of 16k.  Alignment is forced to 128 bytes, as is the stack
   length.

   Consider using a prime no. of kb to reduce cache collisions.

   To speed fiber creation, 64 stacks are cached for reuse.

   Call this routine before making any fibers since we don't bother to
   free already cached stacks.  */
void fiber_set_stacksize(size_t s);

#endif
