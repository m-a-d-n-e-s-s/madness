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

#ifdef _CRAY
#include <catamount/dclock.h>

double walltime() {
    return dclock();
}

double cputime() {
    return dclock();
}
#else

#include <time.h>
#include <sys/time.h>

double walltime() {
    struct timeval tv;
    gettimeofday(&tv,0);
    return tv.tv_sec + 1e-6*tv.tv_usec;
}

double cputime() {
    return clock()*1.0/CLOCKS_PER_SEC;
}
#endif

void set_vec(int n, int *a, int val) {
    int i;
    for (i=0; i<n; i++) a[i] = val;
}

void check_vec(int n, int *a, int val) {
    int i;
    for (i=0; i<n; i++) if (a[i] != val) {
        printf("bad vector value\n");
        exit(1);
    }
}

void hello(void *vmsg) {
    char *msg = (char *) vmsg;
    /* The set/check attempts to verify that stacks are not colliding */
    int a[4096];
    set_vec(4096, a, msg[0]);
    printf("Hello from %s\n", (char *) msg);
    fiber_yield();
    
    check_vec(4096, a, msg[0]);
    printf("I love you from %s\n", (char *) msg);
    fiber_yield();
    
    check_vec(4096, a, msg[0]);
    printf("I've changed my mind from %s\n", (char *) msg);
    fiber_yield();
    
    check_vec(4096, a, msg[0]);
    printf("Goodbye from %s\n", (char *) msg);
}

void noop(void *d) {}

void countdown(void *vn) {
    while (--(*((long *) vn))) 
        fiber_yield();
}

void count_switches(void *data) {
    (*(int *)data)++;
}

int main() {
    int count=0;
    fiber_set_switch_hook(count_switches,&count);
    fiber_create(hello,"Dave");
    fiber_create(hello,"Mary");
    fiber_create(hello,"Fred");
    fiber_create(hello,"Suzy");
    fiber_create(hello,"Nick");
    fiber_create(hello,"Emma");
    
    while (fiber_count())  {
        printf("\nMain thread yielding ...\n");
        fiber_yield();
    }
    printf("\nThere were %d fiber switches.\n",count);
    fiber_set_switch_hook(0,0);

#ifndef DEBUG
    {
        long i,n=100000;
	volatile long vn;
        double wall_used, cpu_used;

        /* Time create+yield+return+destroy sequence */
        printf("\n");
        while (1) {
            wall_used = walltime();
            cpu_used = cputime();
            for (i=0; i<n; i++) {
                fiber_create(noop,0);
                fiber_yield();
            }
            wall_used = walltime() - wall_used;
            cpu_used = cputime() - cpu_used;
            printf("n=%ld wall=%.1e cpu=%.1e\n",n,wall_used,cpu_used);
            if (wall_used > 1.0) break;
            else if (wall_used > 0.01) n = 1.04*n/wall_used;
            else  n = n*10;
        }
        printf("\n time per create+yield+return+destroy = %.1es (wall) %.1es (cpu)\n",
               wall_used/n, cpu_used/n);

        /* Time cost of yield */
	vn = n;
	fiber_create(countdown,(void *) &vn);
	wall_used = walltime();
	cpu_used = cputime();

	while(vn) 
            fiber_yield();

	wall_used = walltime() - wall_used;
	cpu_used = cputime() - cpu_used;
        printf("\n time per yield+countdown+yield = %.1es (wall) %.1es (cpu)\n",
               wall_used/n, cpu_used/n);


        /* Time creation of 100,000 independent fibers */
	wall_used = walltime();
	cpu_used = cputime();
        n = 100000;
        while (n--) fiber_create(noop,0);
	wall_used = walltime() - wall_used;
	cpu_used = cputime() - cpu_used;
        printf("\n time to create 100K independent fibers %.1es = %.1es/fiber\n",
               cpu_used, cpu_used*1e-5);

	wall_used = walltime();
	cpu_used = cputime();
        n = 100000;
        while (fiber_count()) fiber_yield();
	wall_used = walltime() - wall_used;
	cpu_used = cputime() - cpu_used;
        printf("\n time to yield+destroy 100K independent fibers %.1es = %.1es/fiber\n",
               cpu_used, cpu_used*1e-5);
    }
#endif

    return 0;
}
