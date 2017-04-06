# TODO: MPI_Barrier(MPI_COMM_WORLD) around printf's
def _header(complex_c, complex_a, complex_b):
    ret = """
#include <stdio.h>
#include <hwi/include/bqc/A2_inlines.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <pthread.h>
#include <mpi.h>
#include <limits.h>
#include <complex>

#define NAPTIME 3
#define MAX_POSIX_THREADS 64

typedef struct _pt_data {
    int myid;
    int total_threads;
    uint64_t _tc;
    struct _pt_data **args;
    uint64_t tt;
} pt_data;

static pthread_t thread_pool[MAX_POSIX_THREADS];
pt_data *thread_args[MAX_POSIX_THREADS];

static int mpi_size, mpi_rank;
static int num_posix_threads;
pthread_barrier_t barr;

void ran_fill(int n, double *a) {
    while (n--) *a++ = ((double)rand())/RAND_MAX;
}

inline void start_timer(pt_data *pt) {
    pt->_tc = GetTimeBase();
}

inline uint64_t stop_timer(pt_data *pt) {
    return GetTimeBase() - pt->_tc;
}

void mtxm(long dimi, long dimj, long dimk, long extb,
""" + "double {} * c, const double {} * a, const double {} * b) {{".format(complex_c, complex_a, complex_b) + """
    int i, j, k;
    for (k=0; k<dimk; ++k) {
        for (j=0; j<dimj; ++j) {
            for (i=0; i<dimi; ++i) {
                c[i*dimj+j] += a[k*dimi+i]*b[k*extb+j];
            }
        }
    }
}
"""
    return ret

def _tester(i, m, trans=False):
    ret = """
    if (mpi_rank=={i}) {{
        for (loop=0; loop<30; ++loop) {{
            rc = pthread_barrier_wait(&barr);
            if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD) {{
                printf("Could not wait on barrier\\n");
                exit(-1);
            }}

            start_timer(pt);
            {m}(ni,nj,nk,nj,c,a,b);"""
    if trans:
        ret += """
            {m}(ni,nj,nk,nj,a,c,b);
            {m}(ni,nj,nk,nj,c,a,b);"""
    ret += """
            uint64_t tt = stop_timer(pt);
            pt->tt += tt;
        }}
        pt->tt /= 30;

        rc = pthread_barrier_wait(&barr);
        if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD) {{
            printf("Could not wait on barrier\\n");
            exit(-1);
        }}

        if (pt->myid == 0) {{
            uint64_t maxtt = 0; // longest time
            uint64_t mintt = ULONG_MAX; // shortest time
            uint64_t sumtt = 0;

            int i;
            for (i=0; i<pt->total_threads; i++) {{
                uint64_t itt;
                itt = pt->args[i]->tt;
                maxtt = maxtt > itt ? maxtt : itt;
                mintt = mintt < itt ? mintt : itt;
                sumtt += itt;
            }}

            double avgtt = (double)sumtt/(double)i;
            printf("{m} [%2d pthreads]: %20s %3ld %3ld %3ld %8.2f (avg: %8.2f min: %llu max: %llu)\\n", pt->total_threads, s, ni, nj, nk, nflop/maxtt, avgtt, mintt, maxtt);
        }}

        rc = pthread_barrier_wait(&barr);
        if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD) {{
            printf("Could not wait on barrier\\n");
            exit(-1);
        }}
    }}"""
    return ret.format(i=i, m=m)

def _timer(mtxms, complex_c, complex_a, complex_b):
    ret = """void timer(pt_data *pt, const char* s, long ni, long nj, long nk, double {} *a, double {} *b, double {} *c) {{
    double fastest = 0.0;
    double nflop = 2.0{}{}*ni*nj*nk*(pt->total_threads);
    long loop;
    int rc;
    pt->tt = 0;
""".format(complex_a, complex_b, complex_c, complex_a and "*2.0", complex_b and "*2.0")

    for i, m in enumerate(mtxms):
        ret += _tester(i, m)

    ret += "}\n"
    return ret

def _transtimer(mtxms, complex_c, complex_a, complex_b):
    ret = """void trantimer(pt_data *pt, const char* s, long ni, long nj, long nk, double {} *a, double {} *b, double {} *c) {{
    double fastest = 0.0;
    double nflop = 3.0*2.0{}{}*ni*nj*nk*(pt->total_threads);
    long loop;
    int rc;
    pt->tt = 0;
""".format(complex_a, complex_b, complex_c, complex_a and "*2.0", complex_b and "*2.0")

    for i, m in enumerate(mtxms):
        ret += _tester(i, m, True)

    ret += "}\n"
    return ret

def _main(mtxms, complex_c, complex_a, complex_b):
    ret = """void* entry(void* arg) {{
    pt_data *pt = (pt_data*)arg;

    const long nimax=30*30;
    const long njmax=100;
    const long nkmax=100;
    long ni, nj, nk, i, m;

    double {ca} *a;
    double {cb} *b;
    double {cc} *c;
    double {cc} *d;

    posix_memalign((void **) &a, 32, nkmax*nimax*sizeof(double {ca}));
    posix_memalign((void **) &b, 32, nkmax*njmax*sizeof(double {cb}));
    posix_memalign((void **) &c, 32, nimax*njmax*sizeof(double {cc}));
    posix_memalign((void **) &d, 32, nimax*njmax*sizeof(double {cc}));

    ran_fill(nkmax*nimax{cax}, (double*)a);
    ran_fill(nkmax*njmax{cbx}, (double*)b);

    if (pt->total_threads == 1) {{
        for (ni=2; ni<60; ni+=2) {{
            for (nj=2; nj<100; nj+=6) {{
                for (nk=2; nk<100; nk+=6) {{
                    for (i=0; i<ni*nj; ++i) d[i] = c[i] = 0.0;
                    mtxm(ni,nj,nk,nj,c,a,b);
"""
    for i, f in enumerate(mtxms):
        ret += "                    if (mpi_rank=={}) {}(ni,nj,nk,nj,d,a,b);\n".format(i, f)
    ret += """
                    for (i=0; i<ni*nj; ++i) {{
                        double err = abs(d[i]-c[i]);
                        if (err > 1e-{error}) {{
"""
    for i, f in enumerate(mtxms):
        ret += """                            if (mpi_rank=={}) fprintf(stderr, "{} error %d: %ld %ld %ld %e\\n",mpi_rank,ni,nj,nk,err);\n""".format(i, f)
    ret += """
                        }}
                    }}
                }}
            }}
        }}
    }}

    for (ni=2; ni<60; ni+=2) timer(pt, "(m*m)T*(m*m)", ni,ni,ni,a,b,c);
    for ( m=2; m<=30;  m+=2) timer(pt, "(m*m,m)T*(m*m)", m*m,m,m,a,b,c);
"""
    if (complex_a and complex_b) or (not complex_a and not complex_b):
        ret += """    for ( m=2; m<=30;  m+=2) trantimer(pt, "tran(m,m,m)", m*m,m,m,a,b,c);\n"""
    ret += """
    for ( m=2; m<=20;  m+=2) timer(pt, "(20*20,20)T*(20,m)", 20*20,m,20,a,b,c);

    free(a);
    free(b);
    free(c);
    free(d);

    return NULL;
}}

int main(int argc, char *argv[]) {{
    int i, rc;
    int provided;

    MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
    assert(provided == MPI_THREAD_MULTIPLE);

    MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

    srand(argc > 42 ? argc : 1);

    for (num_posix_threads = 1; num_posix_threads <= MAX_POSIX_THREADS; num_posix_threads*=2) {{
        MPI_Barrier(MPI_COMM_WORLD);

        for (i=0; i<num_posix_threads; i++) {{
            thread_args[i] = (pt_data*)malloc(sizeof(pt_data));
            thread_args[i]->myid = i;
            thread_args[i]->total_threads = num_posix_threads;
            thread_args[i]->args = thread_args;
        }}

        if(pthread_barrier_init(&barr, NULL, num_posix_threads)) {{
            printf("Could not create a barrier\\n");
            return -1;
        }}

        for (i=0 ; i<num_posix_threads ; i++) {{
            rc = pthread_create(&thread_pool[i], NULL, entry, thread_args[i]);
            assert(rc==0);
        }}

        sleep(NAPTIME);

        for (i=0 ; i<num_posix_threads ; i++) {{
            rc = pthread_join(thread_pool[i],NULL);
            assert(rc==0);
        }}

        MPI_Barrier(MPI_COMM_WORLD);

        for (i=0; i<num_posix_threads; i++) {{
            free(thread_args[i]);
        }}
    }}

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}}"""
    return ret.format(ca=complex_a, cb=complex_b, cc=complex_c, 
            cax=(complex_a and "*2"), cbx=(complex_b and "*2"), 
            error=(complex_c and "12" or "14"))

def tester_gen_bgq(f, mtxmspecs, mtxm_gen, complex_a, complex_b):
    r = lambda x: x.replace("double complex", "__complex__ double")
    complex_c = (complex_a or complex_b) and "complex" or ""
    complex_a = complex_a and "complex" or ""
    complex_b = complex_b and "complex" or ""
    print(r(_header(complex_c, complex_a, complex_b)), file=f)
    for mtxm in mtxmspecs:
        mtxm_gen(f, *mtxm)
    mtxms = [x[-1] for x in mtxmspecs]
    print(r(_timer(mtxms, complex_c, complex_a, complex_b)), file=f)
    if not (bool(complex_a) ^ bool(complex_b)):
        print(r(_transtimer(mtxms, complex_c, complex_a, complex_b)), file=f)
    print(r(_main(mtxms, complex_c, complex_a, complex_b)), file=f)
