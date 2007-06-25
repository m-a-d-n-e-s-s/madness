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

  
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <fiber.h>

/*
  Example program to demonstrate/test the use of fibers with MPI.
  
  A crude distributed-memory Cilk-like model is implemented with the
  spawn operation and simple futures.  Execution starts with a single
  fiber running in MPI process 0 in the user provided function
  
  int user_main(int argc, char** argv)
  
  When you return from user_main(), the program will end.
  
  You can spawn new activities using this simple call
  
  future* spawn(double (*f)(double), double arg)
  
  f = function to be invoked
  arg = input argument
  
  spawn() creates a new fiber to execute the function in another
  *random* MPI process and returns a pointer to a future in which the
  result will eventually appear.
  
  spawn() is non-blocking, so your current fiber happily continues
  executing, but when you try to access the result (via the future)
  using future_force() you will yield unless the result is 
  already available.
 
  double future_force(future* future)
  
  future_force() returns the result, yielding if it is not yet
  available (you can only call this once per future since this simple
  code also frees the underlying data structure).

  [
   I also added as byproduct of the implementation

   future* spawn_at(double (*f)(double), double arg, int rank)

   which will spawn an activity in a specific MPI process, and

   int future_probe(future *f)

   which returns true if the future is ready, false otherwise.
  ]
  
  ------  
  
  Simplifying assumptions that a better implementation would eliminate
  
  - input argument and the result are doubles whereas a better
  implementation would provide for serialization of general data
  structures
  
  - functions are assumed to have the same address in all MPI
  processes instead of providing a registration mechanism
  
  - random process assignment instead of dynamic load balancing.
  There is also no dependency analysis, migration of fibers between
  ready/not-ready queues, or scheduling of fibers as futures
  become ready.
  
  - unlimited system resources (e.g., outstanding MPI operations)
  which could be fixed by throttling the creation of new fibers 
  with a task queue.

  - dedicated processors - over subscribing the processors will
  cause exponential slowdown due to no backoff while spinning

  - immediately sends all messages rather than attempting aggregation
  or other optimizations

  - stack overflow is not detected
*/

/* A little bit of global state */
static int nproc;  /* The number of MPI processes */
static int me;     /* The MPI rank of this process */ 
static volatile int all_finished = 0; 
static int ntask_done = 0;


typedef struct {
    int rank;             /* Requesting process */
    int tag;              /* Unique tag for reply */
    double arg;           /* Input argument */
    double (*f)(double);  /* Function to invoke */
} request_buf;


typedef struct {
    int cookie;
    MPI_Request handle;    /* From MPI_Irecv */
    int status;            /* 2=ready(assigned), 1=irecv pending,  0=local-op pending */
    double value;          /* The result */
    request_buf request;   /* For convenience keep the original request */
} future;


static int unique_tag() {
    /* Returns a "unique" MPI tag ...  99 > tag < 32868 */
    static int tag = 0;
    return ((tag++)&0x7fff) + 100;
}

    
static void error(char *msg) {
    fprintf(stderr,"%d: fatal error: %s\n",me,msg);
    MPI_Abort(MPI_COMM_WORLD,1);
}


static void poll(void *);

/* Wait for MPI request to complete while polling ... do this where cannot yield */
static void await(MPI_Request *handle) {
    int flag;
    MPI_Status status;
    while (1) {
        MPI_Test(handle, &flag, &status);
        if (flag) 
            break;
        else 
            poll(0);
    }
}


static future *new_future() {
    future *f = malloc(sizeof(future));
    if (!f) error("new_future: failed allocating a future");
    f->cookie = 5551212;
    return f;
}


static void free_future(future *f) {
    if (f->cookie != 5551212) error("free_future: trying to free something that is not a future");
    f->cookie = 0;
    free(f);
}


static void wrapper(void *vr) {
    future *result = vr;
    result->value = result->request.f(result->request.arg);
    ntask_done++;
    if (result->request.rank == me) {      /* Assigning to local future */
        result->status = 2;
    }
    else { /* Remote future ... send value and free local copy */
        static MPI_Request handle;
	static int send_pending = 0;
	if (send_pending) await(&handle);
	if (MPI_Isend(&result->value, sizeof(double), MPI_BYTE, result->request.rank, 
                 result->request.tag, MPI_COMM_WORLD, &handle) != MPI_SUCCESS)
	    error("wrapper: failed posting isend");
	send_pending = 1;
        free_future(result);
    }
}


/* Returns random process != me (unless we only have one process) */
static int random_rank() {
    double ranmax = 32768.0;
    int rank = nproc*((random()&32767)/ranmax);
    return (rank+me)%nproc;
}


/* This routine services all pending requests */
static void poll(void *junk) {
    static int initialized = 0;
    static request_buf request;
    static MPI_Request handle;
    if (!initialized) {
      if (MPI_Irecv(&request, sizeof(request_buf), MPI_BYTE, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD, &handle)
	  != MPI_SUCCESS) error("poll: failed posting irecv(a)");
        initialized = 1;
    }
    
    while (1) {
        future *result;
        int flag;
        MPI_Status status;
        MPI_Test(&handle, &flag, &status);
        if (!flag) return;
        
        if (request.f == 0) {
            all_finished = 1;
            return;
        }
        else {
  	    result = new_future();
            result->request = request;
            result->status = 1;
            if (!fiber_create(wrapper,(void *)result)) error("poll: failed spawning new fiber");
	    if (MPI_Irecv(&request, sizeof(request_buf), MPI_BYTE, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD, &handle)
		!= MPI_SUCCESS) error("poll: failed posting irecv(a)");
        }
    }
}      


future* spawn_at(double (*f)(double), double arg, int rank) {
    future* result = new_future();
    result->request.rank = me;
    result->request.tag = unique_tag();
    result->request.arg = arg;
    result->request.f = f;

    if (nproc > 1) poll(0); /* Should ensure progress */
    
    if (rank == me) {
        result->status = 0;
        if (!fiber_create(wrapper,(void *)result)) error("spawn_at: failed spawning new fiber");
    }
    else {
        static MPI_Request handle;
	static int send_pending = 0;
	if (send_pending) await(&handle);
        result->status = 1;
        if (MPI_Irecv(&(result->value), sizeof(double), MPI_BYTE, rank, 
		       result->request.tag, MPI_COMM_WORLD, &result->handle) != MPI_SUCCESS)
	    error("spawn_at: failed posting irecv");
	if (MPI_Isend(&(result->request), sizeof(request_buf), MPI_BYTE, rank, 
                 99, MPI_COMM_WORLD, &handle) != MPI_SUCCESS) 
	    error("spawn_at: failed posting isend");
	 send_pending = 1;
    }
    return result;
}

future* spawn(double (*f)(double), double arg) {
    return spawn_at(f, arg, random_rank());
}


/* Returns true if the future is ready, false otherwise */
int future_probe(future *f) {
    if (f->status == 1) {
        int flag;
        MPI_Status status;
        MPI_Test(&(f->handle), &flag, &status);
        if (flag) f->status = 2;
    }
    return  (f->status == 2);
}


double future_force(future *f) {
    double result;
    while (!future_probe(f)) fiber_yield();
    result = f->value;
    free_future(f);
    return result;
}


static void end_the_world() {
    if (nproc > 1) {
        int rank;
        request_buf request;
        request.f = 0;
        for (rank=0; rank<nproc; rank++)
            MPI_Send(&request, sizeof(request), MPI_BYTE, rank, 99, MPI_COMM_WORLD);
    }
}

extern int user_main(int argc, char **argv);

int main(int argc, char **argv) {
    /* Create the world */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);

    /* Start the server */
    fiber_set_stacksize(127*1024);
    if (nproc > 1) {
      fiber_set_switch_hook(poll,0);  
      fiber_yield();
    }
    
    /* Execute the user code on process 0 */
    if (me == 0) {
        user_main(argc,argv);
        end_the_world();
    }
    else {
        while (!all_finished) fiber_yield();
    }
    
    printf("Process %d did %d tasks\n",me,ntask_done);
    MPI_Finalize();
    
    return 0;
}

/* **********************************************************
   Above is library code ... below is user code.  First look
   at the main program which starts in a single thread.
   ********************************************************** */

/*
  We will compute using adaptive numerical quadrature the integral
  
  int(f(x), x=0..1)
  
  The arguments to each task are n and l, where n is the level of
  refinement and l is the index of the box at that level
  (l = 0,1,...,2^n-1).  A task estimates the value of the integral
  
  int(f(x), x=l/2^n..(l+1)/2^n)
  
  using both the trapezoid rule and Simpson rules (yes, there are
  better ways of doing this).  If the difference is less than some
  threshold the Simpson result is returned.  Otherwise, the box is
  divided in two and a new task spawned to run one of them.  Via the
  futures, the sum is returned to the original single thread of
  execution.  About 8200 intervals are computed with 4100 tasks.  */

double f(double x) {
    return 100.0*exp(-100*fabs(x-0.5))*cos(12.0*(x-0.5));
}

double task(double arg) {
#define k 20
    int l = arg/256;
    int n = arg - l*256;
    double v[2*k+1];
    double diff, tk, t2k, s, h=(1.0/(1<<n));
    double xlo=h*l, tol=1e-8*h;
    int i;
    h = 0.5*h/k;
    
    for (i=0; i<2*k+1; i++)  v[i] = f(xlo + i*h);
    
    tk = t2k = (v[0] + v[2*k])*0.5;
    for (i=1; i<2*k; i++) t2k += v[i];
    for (i=2; i<2*k; i+=2) tk += v[i];
    t2k *= h;
    tk *= 2*h;
    s = (4.0*t2k - tk)/3.0;
    diff = fabs(s-t2k);
    if (fabs(s) > 1.0) tol *= fabs(s);
    
    /*printf("rank=%d n=%d l=%d trap(k)=%.12f trap(2k)=%.12f simp=%.12f diff=%.12f\n",
      me,n,l,tk,t2k,s,fabs(s-t2k)); */

    /* trapezoid err is O(h^2) and simpson's is O(h^4) */
    if (diff < tol) {
        return s;
    }
    else {
        future *left;
        double right;
        n = n+1;  l *= 2;
        left = spawn(task,n+l*256);    /*  <----- spawn */
        right = task(n+(l+1)*256);
        return right+future_force(left);  /*  <----- force */
    }
}

double greet(double d) {
  return -d;
}

void hello_world() {
  /* Spawn 10 tasks to show args/results are passed correctly */
  future *f[10];
  int i;
  for (i=0; i<10; i++) 
      f[i] = spawn(greet,i);
  for (i=0; i<10; i++) 
      printf("Reply from task %d is %.0f\n",i,future_force(f[i]));
}  

int user_main(int argc, char** argv) {
    future *sum;
    
    printf("Initial user process has started\n");

    hello_world();
    
    sum = spawn(task,0);           /*  <----- spawn first task with n=l=0*/
    
    printf("The numerical integral is %.17f\n",future_force(sum)); /* <----- force */
    printf("          Exact result is.%.17f\n",1.9716088328075709779);
    
    return 0;
}
