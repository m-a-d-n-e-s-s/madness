#include <math.h>
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
#include <string.h>
#include "fiber.h"

#ifdef DEBUG
#define DEBUGCODE(a) a
#else
#define DEBUGCODE(a)
#endif

/* Align stack on cache line */
#define STACK_ALIGNMENT 7

/* Default size for the stack for each fiber - 67kb */
#define STACK_DEFAULT_SIZE 68608
#define STACK_MIN_SIZE 16384
/* Maximum no. of unused contexts and stacks that will be cached */
#define MAX_FREE_CONTEXTS 64


/* To enable many threads with small stacks and for lack of any
   concrete performance evidence I removed the code that offsets the
   stack of each fiber to reduce cache conflicts.  However, I did make
   the default stack size 67kb (a prime number, rather than 64kb) so
   that up to 64 fibers started in a tight loop are offset in L2 cache
   (Intel Xeon) by 1k.  If you suspect cache collisions are giving a
   performance problem, let me know.
*/

#ifdef HAVE_UCONTEXT
#include <ucontext.h>
typedef ucontext_t context;
static void setstack(context *c, void *s, size_t stack_size) {
    c->uc_stack.ss_sp = s;
    c->uc_stack.ss_size = stack_size;
    c->uc_link = 0; /* Does not really belong here but is convenient */
}

#elif defined(HAVE_X86_ASM)

/* This and the 64-bit version are not a general purpose
   implementations of the ucontext interface.  They are tailored to
   suit our purpose only. */

typedef struct context_struct {
    unsigned long *sp;    /* The stack pointer */
} context;

#define getcontext(a) 

#define setstack(c, s, sz)  (c)->sp = (unsigned long *) (((char *) s) + sz)

#define makecontext(c, f, argc) \
    (c)->sp[-5]= (unsigned long) f; \
    (c)->sp[-4]= (unsigned long) ((c)->sp) ; \
    (c)->sp[-3] = 0; \
    (c)->sp[-2] = 0; \
    (c)->sp[-1] = 0; \
    (c)->sp -= 5

#define swapcontext(old,new) \
    __asm__ __volatile__ ( \
                 "pushl %%ebx;" \
                 "pushl %%esi;" \
                 "pushl %%edi;" \
                 "pushl %%ebp;" \
                 "leal L_pgcc_is_asleep, %%ebp;" \
                 "pushl %%ebp;" \
                 "movl  %%esp,(%0);" \
                 \
                 "movl %1,%%esp;" \
                 "popl %%eax;" \
                 "popl %%ebp;" \
                 "popl %%edi;" \
                 "popl %%esi;" \
                 "popl %%ebx;" \
                 \
                 "jmp *%%eax;" \
                 \
                 "L_pgcc_is_asleep:" \
                 : \
                 : "a"(&((old)->sp)),"d"((new)->sp) \
                 : "memory" \
                 )

#define setcontext(new) \
    __asm__ __volatile__ ( \
                 "movl %0,%%esp;" \
                 "popl %%eax;" \
                 "popl %%ebp;" \
                 "popl %%edi;" \
                 "popl %%esi;" \
                 "popl %%ebx;" \
                 "jmp *%%eax;" \
                 : \
                 : "a"((new)->sp) \
                 : "memory" \
                 )

#elif defined(HAVE_X8664_ASM)
typedef struct context_struct {
    unsigned long *sp;    /* The stack pointer */
} context;

#define getcontext(a) 

/* -8 since AMD64 (%rsp-8) should be a multiple of 16 */
#define setstack(c, s, sz)  (c)->sp = (unsigned long *) (((char *) s) + sz - 8)

#define makecontext(c, f, argc) \
    (c)->sp[-9]= (unsigned long) f; \
    (c)->sp[-8]= (unsigned long) ((c)->sp) ; \
    (c)->sp[-7] = 0; \
    (c)->sp[-6] = 0; \
    (c)->sp[-5] = 0; \
    (c)->sp[-4] = 0; \
    (c)->sp[-3] = 0; \
    (c)->sp[-2] = 0; \
    (c)->sp[-1] = 0; \
    (c)->sp -= 9

#define swapcontext(old,new) \
    __asm__ __volatile__ ( \
                 "pushq %%r12;" \
                 "pushq %%r13;" \
                 "pushq %%r14;" \
                 "pushq %%r15;" \
                 "pushq %%rbx;" \
                 "pushq %%rsi;" \
                 "pushq %%rdi;" \
                 "pushq %%rbp;" \
                 "leaq L_pgcc_is_asleep, %%rbp;" \
                 "pushq %%rbp;" \
                 "movq  %%rsp,(%%rax);" \
                 \
                 "movq %%rdx,%%rsp;" \
                 "popq %%rax;" \
                 "popq %%rbp;" \
                 "popq %%rdi;" \
                 "popq %%rsi;" \
                 "popq %%rbx;" \
                 "popq %%r15;" \
                 "popq %%r14;" \
                 "popq %%r13;" \
                 "popq %%r12;" \
                 \
                 "jmpq *%%rax;" \
                 \
                 "L_pgcc_is_asleep:" \
                 : \
                 : "a"(&((old)->sp)),"d"((new)->sp) \
                 : "memory" \
                 )

#define setcontext(new) \
    __asm__ __volatile__ ( \
                 "movq %%rax,%%rsp;" \
                 "popq %%rax;" \
                 "popq %%rbp;" \
                 "popq %%rdi;" \
                 "popq %%rsi;" \
                 "popq %%rbx;" \
                 "popq %%r15;" \
                 "popq %%r14;" \
                 "popq %%r13;" \
                 "popq %%r12;" \
                 "jmpq *%%rax;" \
                 : \
                 : "a"((new)->sp) \
                 : "memory" \
                 )
#elif defined(HAVE_PPC_ASM) 
/* 32 and 64-bit powerpc and, I think, power4 and 5 currently on AIX only */
extern void *fiber_make_context(void *, void(*)());
extern void fiber_swap_context(void **, void*); 
extern void fiber_set_context(void *);

typedef struct context_struct {
    unsigned long *sp;    /* The stack pointer */
} context;

#define getcontext(a) 
#define setstack(c, s, sz)  (c)->sp = (unsigned long *) (((char *) s) + sz)
#define makecontext(c, f, argc) (c)->sp = fiber_make_context((c)->sp, f)
#define swapcontext(old,new) fiber_swap_context((void *) &((old)->sp),(new)->sp)
#define setcontext(new) fiber_set_context((new)->sp)

#elif defined(HAVE_LONGJMP)

/* Incomplete ucontext interface on top of set/longjmp.  In particular, 
   no arguments and you CANNOT return from the function. */
#include <setjmp.h>
typedef struct ucontext_struct {jmp_buf env;} context;

static void setstack(context *c, void *s, size_t stack_size) {
    void **p = (void **) c->env;
    s = (char *)s + stack_size;
    
#ifdef CYGWIN
#define JB_SP 7
#define JB_PC 8
#endif

#ifdef JB_RSP
    p[JB_RSP] = s - 8; /* -8 to satisfy AMD64 alignment request */
#elif defined(JB_SP)
    p[JB_SP] = s;
#else
    /* Cannot figure out the jmp_buf structure */
    error this_is_very_unfortunate;
#endif
}
#define swapcontext(old,new) if (setjmp((old)->env) == 0) longjmp((new)->env, 1)
#define setcontext(new) longjmp((new)->env, 1)
#define getcontext(c) setjmp((c)->env)
#define makecontext(c, f, argc)  ((void **) (c)->env)[JB_PC] = (void *) (f)

/* The following procs are inlined via the above macros
static void swapcontext(context *old, context *new) {
    if (setjmp(old->env) == 0) longjmp(new->env, 1);
 }
static void getcontext(context *c) {setjmp(c->env);}
static void makecontext(context *c, void (*f)(), int argc) {
    void **p = (void **) c->env;
    p[JB_PC] = (void *) f;
 }
*/

#else
    /* Need a minimal implementation of the ucontext interface */
    error this_is_very_unfortunate;
#endif

typedef struct fiber_context_struct {
    struct fiber_context_struct *prev; 
    struct fiber_context_struct *next;
    void *stack;               /* Points to UNALIGNED stack for free()ing */
    void *actual_stack;        /* Points to ALIGNED stack for use */
    size_t actual_stack_size;  /* Actual space available in stack */
    void (*f)(void *);
    void *data;
    context c;
} fiber_context;


/* This stuff CANNOT be static once we map to multiple kernel threads.
   We will need some thread-private data to hold this which means
   participating in the thread creation or providing an initialization
   routine. */
static size_t default_stack_size = STACK_DEFAULT_SIZE;
static fiber_context *fibers = 0; /* Points to the current fiber */
static fiber_context *free_contexts[MAX_FREE_CONTEXTS]; /* Caches to unused structures */
static int nfree_contexts = 0; /* No. of free contexts */
static int fiber_num = 0; /* Counts no. of created fibers not yet completed */
static void (*fiber_switch_hook)(void *data) = 0;
static void *fiber_switch_hook_data = 0;

#ifdef DEBUG
static void fiber_error(const char *msg) {
    fprintf(stderr,"fiber: fatal internal error: %s\n",msg);
    exit(1);
}
#endif

void fiber_set_switch_hook(void (*f)(void *), void *data) {
    fiber_switch_hook = f;
    fiber_switch_hook_data = data;
}

static void fiber_initialize() {
    fibers = malloc(sizeof(fiber_context));
    DEBUGCODE(printf("fiber: initializing %p\n", (void *) fibers);)
    fibers->next = fibers->prev = fibers;
}

#ifdef DEBUG
static void fiber_print() {
    fiber_context *p;
    if (! fibers) fiber_initialize();
    printf("There are %d fibers in addition to the main fiber\n", fiber_num);
    p = fibers;
    while (1) {
        printf("   fiber %p <-- %p --> %p\n", (void *) p->prev, (void *) p, (void *) p->next);
        p = p->next;
        if (p == fibers) break;
    }
}
#endif

static void fiber_wrapper() {
    fiber_context *c;

    fibers->f(fibers->data);  /* Calls the user routine here */

    if (fiber_switch_hook) fiber_switch_hook(fiber_switch_hook_data);

#ifdef DEBUG
    printf("fiber: returned from %p and invoking %p\n", (void *) fibers, (void *) fibers->next);
    if (fibers == fibers->next) fiber_error("returning from last fiber - logic error?");
    if (fibers != fibers->next->prev) fiber_error("corrupt data structure in wrapper\n");
#endif

    fibers = fibers->next;

    /* Delete the just terminated fiber.  Cannot free the current so put in the cache  */
    c = fibers->prev;
    c->prev->next = c->next;
    c->next->prev = c->prev;
    if (nfree_contexts >= MAX_FREE_CONTEXTS) free(free_contexts[--nfree_contexts]->stack);
    free_contexts[nfree_contexts++] = c;

    fiber_num--;
    DEBUGCODE(if (fiber_num < 0) fiber_error("negative number of active fibers in wrapper");)

    /* Don't need to save state here since will never return */
    setcontext(&fibers->c);
}

/* Create a new fiber which will invoke f(data).
Returns true on success, false on failure */
int fiber_create(void (*f)(void *), void *data) {
    fiber_context *n;
    if (!fibers) fiber_initialize();
    
    /* Look for an unused context, or create a new one */
    if (nfree_contexts > 0) {
        nfree_contexts--;
        n = free_contexts[nfree_contexts];
        DEBUGCODE(printf("fiber: create: reusing %p\n", n);)
    }
    else {
        /* To only call malloc once, put the data structure before the fiber stack */
        int shift, size;
        char *stack = malloc(default_stack_size);
        if (!stack) return 0;

        /* Align the data structure on a cache line at the start of the stack */
        shift = ((size_t) stack) & ((1<<STACK_ALIGNMENT)-1);
        if (shift) shift = (1<<STACK_ALIGNMENT) - shift;
        n = (void *) (stack + shift);  
        n->stack = stack;  /* Save original pointer for freeing */
        stack = (char *) n;
        stack += sizeof(fiber_context);

        /* Align the actual stack on a cache line */
        shift = ((size_t) stack) & ((1<<STACK_ALIGNMENT)-1);
        if (shift) shift = (1<<STACK_ALIGNMENT) - shift;
        n->actual_stack = (void *) (stack + shift);

        /* Compute the actual stack size and make a multiple of a cache line so that 
           the start and end are aligned since on most machines the stack grows down */
        size = default_stack_size - ((char *) (n->actual_stack)-(char *) (n->stack));
        size = (size>>STACK_ALIGNMENT)<<STACK_ALIGNMENT;
        n->actual_stack_size = size;
        DEBUGCODE(printf("fiber: s=%p n=%p as=%p size=%ld\n",n->stack,n,n->actual_stack,n->actual_stack_size);)
    }

#ifdef DEBUG
    memset(n->actual_stack,127,n->actual_stack_size);
#endif

    DEBUGCODE(printf("fiber: creating %p\n", n);)

    n->f = f;
    n->data = data;
    n->next = fibers;
    n->prev = fibers->prev;
    n->prev->next = n;
    fibers->prev = n;

    /* We will actually call fiber_wrapper which calls the function.
       When the function returns into fiber_wrapper we clean up and
       then swap into the next context.  fiber_wrapper never
       actually returns.  */
    
    getcontext(&n->c);
    setstack(&n->c, n->actual_stack, n->actual_stack_size);
    makecontext(&n->c, fiber_wrapper, 0);
    
    DEBUGCODE(if (fiber_num < 0) fiber_error("negative number of active fibers in create");)
    fiber_num++;
    
    return 1;
}

/* Yield to the next fiber in the FIFO queue.  Execution will resume by returning
   from fiber_yield() */
void fiber_yield() {
    if (! fibers) fiber_initialize();
    if (fiber_switch_hook) fiber_switch_hook(fiber_switch_hook_data);
#ifdef DEBUG
    printf("fiber: yielding %p to %p\n", (void *) fibers, (void *) fibers->next);
    if (fibers != fibers->next->prev)  {
        fiber_print();
        fiber_error("corrupt data structure in yield\n");
    }
#endif
    if (fibers->next != fibers) {
        context *old=&(fibers->c);  /* Local variable for old fixes PGI -O2 optimization */
        fibers = fibers->next;
        swapcontext(old,&(fibers->c));
    }
}

int fiber_count() {return fiber_num;}

/* Get the default stack */
size_t fiber_get_stacksize() {return default_stack_size;}

/* Set the fiber stack size.  Default is 67k, and a minimum of 16k.
   Space for the fiber data structure and alignment are subtracted 
   from this value.
*/
void fiber_set_stacksize(size_t s) {
    if (s < STACK_MIN_SIZE) s = STACK_MIN_SIZE;
    s = (((s-1)>>STACK_ALIGNMENT) + 1)<<STACK_ALIGNMENT;
    default_stack_size = s;
    DEBUGCODE(printf("fiber: set_stack_size: %ul\n", s);)
}
