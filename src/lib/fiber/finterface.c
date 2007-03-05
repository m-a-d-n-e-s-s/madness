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



#include <stdio.h>
#include <stdlib.h>
#include "fiber.h"

/*
 * interface to fortran
 */

#ifndef FORTRAN_INTEGER
#define FORTRAN_INTEGER int
#endif

typedef FORTRAN_INTEGER integer;


integer fiber_create__(void (*f) (void *), void *data)
{
	return(fiber_create( f, data ));
}

integer fiber_create_(void (*f) (void *), void *data)
{
	return(fiber_create(f,data));
}

void fiber_yield__()
{
	fiber_yield();
}

void fiber_yield_()
{
	fiber_yield();
}

integer fiber_count__()
{
	return( fiber_count() );
}

integer fiber_count_()
{
	return( fiber_count() );
}

void fiber_set_stacksize__(integer *s )
{
        size_t ss = *s; 
	fiber_set_stacksize(ss);
}

void fiber_set_stacksize_(integer *s )
{
        size_t ss = *s; 
	fiber_set_stacksize(ss);
}
