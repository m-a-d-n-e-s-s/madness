def _tester(i, m, trans=False):
    ret = """
    if (myid=={i}) {{
#ifdef HAVE_BGP
        HPM_Start("{m}");
#endif
        for (loop=0; loop<30; ++loop) {{
            double rate;
            start_timer();
            {m}(ni,nj,nk,c,a,b);"""
    if trans:
        ret += """
            {m}(ni,nj,nk,a,c,b);
            {m}(ni,nj,nk,c,a,b);"""
    ret +="""
            unsigned long long tt = stop_timer();
            rate = nflop/tt;
            if (rate > fastest) fastest = rate;
        }}

#ifdef HAVE_BGP
        HPM_Stop("{m}");
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        printf("{m}: %20s %3ld %3ld %3ld %8.2f\\n", s, ni, nj, nk, fastest);
#ifdef HAVE_BGP
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }}
"""
    return ret.format(i=i, m=m)


def _timer(mtxms, complex_c, complex_a, complex_b):
    ret = """void timer(const char* s, long ni, long nj, long nk, double {} *a, double {} *b, double {} *c) {{
    double fastest = 0.0;
    double nflop = 2.0{}{}*ni*nj*nk;
    long loop;
""".format(complex_a, complex_b, complex_c, complex_a and "*2.0", complex_b and "*2.0")

    for i, m in enumerate(mtxms):
        ret += _tester(i, m)

    ret += "}\n"
    return ret


def _transtimer(mtxms, complex_c, complex_a, complex_b):
    ret = """void trantimer(const char* s, long ni, long nj, long nk, double {} *a, double {} *b, double {} *c) {{
    double fastest = 0.0;
    double nflop = 3.0*2.0{}{}*ni*nj*nk;
    long loop;
""".format(complex_a, complex_b, complex_c, complex_a and "*2.0", complex_b and "*2.0")

    for i, m in enumerate(mtxms):
        ret += _tester(i, m, True)

    ret += "}\n"
    return ret


def _header(bg, complex_c, complex_a, complex_b):
    ret = ""
    if bg:
        ret += "#define HAVE_BGP 1\n"
    ret += """#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef HAVE_BGP
#include <builtins.h>
#include <complex>

#include "mpi.h"

extern "C" void HPM_Init(void);           // initialize the UPC unit
extern "C" void HPM_Start(char *label);   // start counting in a block marked by the label
extern "C" void HPM_Stop(char *label);    // stop counting in a block marked by the label
extern "C" void HPM_Print(void);          // print counters for all blocks
extern "C" void HPM_Print_Flops(void);

#else

#include <time.h>
#include <immintrin.h>
#include <complex.h>

#endif

int numprocs, myid;

double ran()
{
    static unsigned long seed = 76521;
    seed = seed *1812433253 + 12345;
    return ((double) (seed & 0x7fffffff)) * 4.6566128752458e-10;
}

void ran_fill(int n, double *a) {
    while (n--) *a++ = ran();
}

void mtxm(long dimi, long dimj, long dimk,
    """ + "double {} * c, const double {} * a, const double {} * b)".format(complex_c, complex_a, complex_b) + """ {
    int i, j, k;
    for (k=0; k<dimk; ++k) {
        for (j=0; j<dimj; ++j) {
            for (i=0; i<dimi; ++i) {
                c[i*dimj+j] += a[k*dimi+i]*b[k*dimj+j];
            }
        }
    }
}

unsigned long long _ts;

#if HAVE_BGP 

static __inline__ unsigned long long rdtsc(void)
{
  unsigned long long int result=0;
  unsigned long int upper, lower,tmp;
  __asm__ volatile(
                "0:            "
                "mftbu   %0    "
                "mftb    %1    "
                "mftbu   %2    "
                "cmpw    %2,%0 "
                "bne     0b    "
                : "=r"(upper),"=r"(lower),"=r"(tmp)
                );
  result = upper;
  result = result<<32;
  result = result|lower;

  return(result);
}

#else

static __inline__ unsigned long long rdtsc(void)
{
    unsigned long long int x;
    __asm__ volatile ("rdtsc; shlq $32, %%rdx; orq %%rdx, %0" : "=a" (x) : : "rdx");
    return x;
}

#endif

inline void start_timer() {
    _ts = rdtsc();
}

inline unsigned long long stop_timer() {
    return rdtsc() - _ts;
}
"""
    return ret

def _main(mtxms, complex_c, complex_a, complex_b):
    ret = """int main(int argc, char **argv) {
    const long nimax=30*30;
    const long njmax=100;
    const long nkmax=100;
    long ni, nj, nk, i, m;

    """ + "double {} *a;".format(complex_a) + """
    """ + "double {} *b;".format(complex_b) + """
    """ + "double {} *c;".format(complex_c) + """
    """ + "double {} *d;".format(complex_c) + """

#ifdef HAVE_BGP
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);

    HPM_Init();
#else
    numprocs = 1;
    myid = 0;
#endif

    posix_memalign((void **) &a, 16, nkmax*nimax*sizeof(double """ + complex_a + """));
    posix_memalign((void **) &b, 16, nkmax*njmax*sizeof(double """ + complex_b + """));
    posix_memalign((void **) &c, 16, nimax*njmax*sizeof(double """ + complex_c + """));
    posix_memalign((void **) &d, 16, nimax*njmax*sizeof(double """ + complex_c + """));

    ran_fill(nkmax*nimax""" + (complex_a and "*2") + """, (double*)a);
    ran_fill(nkmax*njmax""" + (complex_a and "*2") + """, (double*)b);

    for (ni=2; ni<60; ni+=2) {
        for (nj=2; nj<100; nj+=6) {
            for (nk=2; nk<100; nk+=6) {
                for (i=0; i<ni*nj; ++i) d[i] = c[i] = 0.0;
                mtxm(ni,nj,nk,c,a,b);
"""
    for i, f in enumerate(mtxms):
        ret += "                if (myid == {}) {}(ni,nj,nk,d,a,b);\n".format(i, f)
    ret += """
                for (i=0; i<ni*nj; ++i) {
                    double err = """ + (complex_c and "c" or "f") + """abs(d[i]-c[i]);
                    if (err > 1e-""" + (complex_c and "12" or "14") +""") {
"""
    for i, f in enumerate(mtxms):
        ret += """                        if (myid == {}) fprintf(stderr, "{} error %d: %ld %ld %ld %e\\n",myid,ni,nj,nk,err);\n""".format(i, f)
    ret += """
                    }
                }
            }
        }
    }
#ifdef HAVE_BGP
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    for (ni=2; ni<60; ni+=2) timer("(m*m)T*(m*m)", ni,ni,ni,a,b,c);
    for ( m=2; m<=30;  m+=2) timer("(m*m,m)T*(m*m)", m*m,m,m,a,b,c);
    """ + ((not (bool(complex_a) ^ bool(complex_b))) and """for ( m=2; m<=30;  m+=2) trantimer("tran(m,m,m)", m*m,m,m,a,b,c);""" or "") + """
    for ( m=2; m<=20;  m+=2) timer("(20*20,20)T*(20,m)", 20*20,m,20,a,b,c);

#ifdef HAVE_BGP
    HPM_Print();
    HPM_Print_Flops();

    MPI_Finalize();
#endif

    return 0;
}"""
    return ret

def tester_gen(f, mtxmspecs, mtxm_gen, complex_a=False, complex_b=False, bg=True):
    r = lambda x: x
    if bg:
        r = lambda x: x.replace("double complex", "__complex__ double").replace("cabs", "abs")
    complex_c = (complex_a or complex_b) and "complex" or ""
    complex_a = complex_a and "complex" or ""
    complex_b = complex_b and "complex" or ""
    print(r(_header(bg, complex_c, complex_a, complex_b)), file=f)
    for mtxm in mtxmspecs:
        mtxm_gen(f, *mtxm)
    mtxms = [x[-1] for x in mtxmspecs]
    print(r(_timer(mtxms, complex_c, complex_a, complex_b)), file=f)
    if not (bool(complex_a) ^ bool(complex_b)):
        print(r(_transtimer(mtxms, complex_c, complex_a, complex_b)), file=f)
    print(r(_main(mtxms, complex_c, complex_a, complex_b)), file=f)
