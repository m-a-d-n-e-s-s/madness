#define WORLD_INSTANTIATE_STATIC_TEMPLATES

#include <madness/madness_config.h>
#include <madness/world/MADworld.h>
#include <madness/tensor/distributed_matrix.h>

using namespace madness;

double ij(int64_t i, int64_t j) {return ((i<<24) | j);}

void check(DistributedMatrix<double>& A) {
    A.fill(ij);

    // Verify local data and accessors
    int64_t ilo, ihi, jlo, jhi;
    A.local_colrange(ilo,ihi);
    A.local_rowrange(jlo,jhi);
    const Tensor<double>& t = A.data();
    for (int64_t i=ilo; i<=ihi; i++) {
        for (int64_t j=jlo; j<=jhi; j++) {
            MADNESS_ASSERT(A.get(i,j) == ij(i,j));
            MADNESS_ASSERT(t(i-ilo,j-jlo) == ij(i,j));
       }
    }

    // Verify ownership computation
    const ProcessID me = A.get_world().rank();
    const int64_t n=A.coldim(), m=A.rowdim();
    for (int64_t i=0; i<n; i++) {
        for (int64_t j=0; j<m; j++) {
            if (A.owner(i,j) == me) MADNESS_ASSERT(A.get(i,j) == ij(i,j));
        }
    }
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);

    const int64_t n=1031, m=977;

    {
        DistributedMatrix<double> A = row_distributed_matrix<double>(world, n, m, 13);
        check(A);
    }
    {
        DistributedMatrix<double> A = column_distributed_matrix<double>(world, n, m, 13);
        check(A);
    }

    world.gop.fence();
    finalize();
    return 0;
}
