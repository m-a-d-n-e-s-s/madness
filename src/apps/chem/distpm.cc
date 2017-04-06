#include <madness/madness_config.h>
#include <madness/mra/mra.h>
#include <madness/mra/lbdeux.h>
#include <madness/world/MADworld.h>
#include <madness/misc/ran.h>
#include <madness/tensor/tensor.h>
#include <madness/tensor/systolic.h>

#include <utility>
#include <vector>

//#ifdef HAVE_INTEL_TBB
//#define NTHREAD 1
//#else
#define NTHREAD ThreadPool::size()+1
//#endif

namespace madness {

typedef Tensor<double> tensorT;
typedef Function<double,3> functionT;
typedef std::vector<functionT> vecfuncT;

/// Simple (?) version of BLAS-1 DROT(N, DX, INCX, DY, INCY, DC, DS)
void drot(long n, double* restrict a, double* restrict b, double s, double c, long inc) {
    if (inc == 1) {
        for (long i=0; i<n; ++i) {
            double aa = a[i]*c - b[i]*s;
            double bb = b[i]*c + a[i]*s;
            a[i] = aa;
            b[i] = bb;
        }
    }
    else {
        for (long i=0; i<(n*inc); i+=inc) {
            double aa = a[i]*c - b[i]*s;
            double bb = b[i]*c + a[i]*s;
            a[i] = aa;
            b[i] = bb;
        }
    }
}

template <typename T, std::size_t NDIM>
void matrix_inner(DistributedMatrix<T>& A,
                  const std::vector< Function<T,NDIM> >& f,
                  const std::vector< Function<T,NDIM> >& g,
                  bool sym=false)
{
    const int64_t n = A.coldim();
    const int64_t m = A.rowdim();
    MADNESS_ASSERT(int64_t(f.size()) == n && int64_t(g.size()) == m);

    // Assume we can always create an ichunk*jchunk matrix locally
    const int ichunk = 1000;
    const int jchunk = 1000; // 1000*1000*8 = 8 MBytes
    for (int64_t ilo=0; ilo<n; ilo+=ichunk) {
        int64_t ihi = std::min(ilo + ichunk, n);
        std::vector< Function<T,NDIM> > ivec(f.begin()+ilo, f.begin()+ihi);
        for (int64_t jlo=0; jlo<m; jlo+=jchunk) {
            int64_t jhi = std::min(jlo + jchunk, m);
            std::vector< Function<T,NDIM> > jvec(g.begin()+jlo, g.begin()+jhi);

            Tensor<T> P = matrix_inner(A.get_world(),ivec,jvec);
            A.copy_from_replicated_patch(ilo, ihi-1, jlo, jhi-1, P);
        }
    }
}

// Computes sum(mu,nu on atom a) C[mu,i] S[mu,nu] C[nu,j]
static inline double PM_q(const tensorT & S, const double * restrict Ci, const double * restrict Cj, int lo, int nbf)
{
    double qij = 0.0;
    if (nbf == 1) { // H atom in STO-3G ... often lots of these!
        qij = Ci[lo]*S(0,0)*Cj[lo];
    }
    else {
        for(int mu = 0;mu < nbf;++mu){
            double Smuj = 0.0;
            for(int nu = 0;nu < nbf;++nu){
                Smuj += S(mu, nu) * Cj[nu + lo];
            }
            qij += Ci[mu + lo] * Smuj;
        }
    }

    return qij;
}

// Given a rotation matrix U, re-order rows to make the matrix more diagonal corresponding to least rotation
class SystolicFixOrbitalOrders : public SystolicMatrixAlgorithm<double> {
    AtomicInt nswitched;

public:
    SystolicFixOrbitalOrders(DistributedMatrix<double>& U, int tag=5556)
        : SystolicMatrixAlgorithm<double>(U, tag, NTHREAD)
    {}

    void start_iteration_hook(const TaskThreadEnv& env) {
        if (env.id() == 0) nswitched = 0;
    }

    void end_iteration_hook(const TaskThreadEnv& env) {
        if (env.id() == 0) {
            int nsw = nswitched;
            SystolicMatrixAlgorithm<double>::get_world().gop.sum(nsw);
            nswitched = nsw;
        }
    }

    bool converged(const TaskThreadEnv& env) const {
        return nswitched == 0;
    }

    void kernel(int i, int j, double * restrict Ui, double * restrict Uj) {
        const int m = get_rowdim();
        double sold = Ui[i]*Ui[i] + Uj[j]*Uj[j];
        double snew = Ui[j]*Ui[j] + Uj[i]*Uj[i];
        if (snew > sold) {
            nswitched++;
            //print("rotation", i, j, sold, snew);
            double tmp[m];
            memcpy(tmp, Ui, m*sizeof(double));
            memcpy(Ui, Uj, m*sizeof(double));
            memcpy(Uj, tmp, m*sizeof(double));
        }
        // While here fix the phases
        if (Ui[i] < 0) {
            //print("fixing sign", i);
            for (int i=0; i<m; i++) Ui[i]*=-1.0;
        }
        if (Uj[j] < 0) {
            //print("fixing sign", j);
            for (int j=0; j<m; j++) Uj[j]*=-1.0;
        }
    }
};

class SystolicPMOrbitalLocalize : public SystolicMatrixAlgorithm<double> {
    const std::vector<int>& set;
    const std::vector<int>& at_to_bf;
    const std::vector<int>& at_nbf;
    const std::vector<tensorT>& Svec;
    const double thresh;
    const double thetamax;
    double tol;
    const int natom;
    const int nao;
    const int nmo;
    int iter;
    AtomicInt ndone_iter;


    // Applies rotation between orbitals i and j for Pipek Mezy
    void localize_PM_ij(const int seti, const int setj,
                        double * restrict Ci, double * restrict Cj,
                        double * restrict Ui, double * restrict Uj)
    {
        if(seti == setj){
            std::vector<double> Qi(natom), Qj(natom);
            double ovij = 0.0;
            for(long a = 0; a < natom; ++a) {
                Qi[a] = PM_q(Svec[a], Ci, Ci, at_to_bf[a], at_nbf[a]);
                Qj[a] = PM_q(Svec[a], Cj, Cj, at_to_bf[a], at_nbf[a]);
                ovij += fabs(Qi[a] * Qj[a]);
            }

            if(ovij > tol * tol){
                double aij = 0.0;
                double bij = 0.0;
                for(long a = 0;a < natom;++a){
                    double qiia = Qi[a];
                    double qija = PM_q(Svec[a], Ci, Cj, at_to_bf[a], at_nbf[a]);
                    double qjja = Qj[a];
                    double d = qiia - qjja;
                    aij += qija * qija - 0.25 * d * d;
                    bij += qija * d;
                }

		// double theta = 0.25*bij/aij; // initial estimate for step restriction
		// if (theta > thetamax) {
		//   theta = thetamax;
		// }
		// else if (theta < -thetamax) {
		//   theta = -thetamax;
		// }
		// else {
		//   theta = 0.25*atan(bij/aij);
		// }

		double theta, fa=fabs(aij), fb=fabs(bij), r=fb/aij;
		// Full formula loses accuracy for b<<a. use taylor series instead
		if (fb < 1e-2*fa) {
                    theta = -0.25*r*(1.0 - r*r/3.0 + r*r*r*r/5.0);
                    //theta = -0.25*fabs(bij)/aij;
		}
		else {
		  theta = 0.25 * acos(-aij / sqrt(aij * aij + bij * bij));
		}
		
                if(bij > 0.0) theta = -theta;

                if(theta > thetamax)
                    theta = thetamax;
                else if(theta < -thetamax)
		    theta = -thetamax;

		if(fabs(theta) >= tol){
		    //print(theta, aij, bij);
		    ndone_iter++;
                    double c = cos(theta);
                    double s = sin(theta);
                    drot(nao, Ci, Cj, s, c, 1);
                    drot(nmo, Ui, Uj, s, c, 1);
                    // for(long a = 0;a < natom;++a){
                    //     Qi[a] = PM_q(Svec[a], Ci, Ci, at_to_bf[a], at_nbf[a]);
                    //     Qj[a] = PM_q(Svec[a], Cj, Cj, at_to_bf[a], at_nbf[a]);
                    // }
                }
            }
        }
    }


public:

    // A[i,...] = [ C[i,...],  U[i,...], Q[i,...] ]
    SystolicPMOrbitalLocalize(DistributedMatrix<double>& A,
                              const std::vector<int>& set,
                              const std::vector<int>& at_to_bf,
                              const std::vector<int>& at_nbf,
                              const std::vector<tensorT>& Svec,
                              double thresh,
                              double thetamax,
                              int natom,
                              int nao,
                              int nmo,
                              int tag=5555)
    : SystolicMatrixAlgorithm<double>(A, tag, NTHREAD),
          set(set),
          at_to_bf(at_to_bf),
          at_nbf(at_nbf),
          Svec(Svec),
          thresh(thresh),
          thetamax(thetamax),
          tol(0.1),
          natom(natom),
          nao(nao),
          nmo(nmo),
          iter(-1)
    {
        MADNESS_ASSERT(A.is_column_distributed());
        MADNESS_ASSERT(A.coldim() == nmo);
        MADNESS_ASSERT(A.rowdim() == nao + nmo);
    }

    void start_iteration_hook(const TaskThreadEnv& env) {
        if (env.id() == 0) {
            iter++;
            //if (iter > 0) tol = std::max(0.1 * std::min(maxtheta, tol), thresh);
            if (iter > 0) tol = std::max(0.333 * tol, thresh);
            ndone_iter = 0;
            //madness::print("start", SystolicMatrixAlgorithm::get_world().rank(),iter,tol);
        }
    }

    void end_iteration_hook(const TaskThreadEnv& env) {
        if(env.id() == 0) {
            int ndone = ndone_iter;
            SystolicMatrixAlgorithm<double>::get_world().gop.sum(ndone);
            ndone_iter = ndone;
            //madness::print("end", SystolicMatrixAlgorithm::get_world().rank(),iter,ndone);
        }
    }

    bool converged(const TaskThreadEnv& env) const {
        //if (env.id() == 0) madness::print("converged", SystolicMatrixAlgorithm::get_world().rank(),iter,int(ndone_iter), tol, thresh, (ndone_iter == 0 && tol == thresh));
        return (ndone_iter == 0 && tol == thresh);
    }

    void kernel(int i, int j, double * restrict rowi, double * restrict rowj) {

        double * restrict Ci = rowi;
        double * restrict Cj = rowj;
        double * restrict Ui = Ci + nao;
        double * restrict Uj = Cj + nao;

        localize_PM_ij(set[i], set[j],
                       Ci, Cj,
                       Ui, Uj);
    }
};


DistributedMatrix<double> distributed_localize_PM(World & world,
                                                  const vecfuncT & mo,
                                                  const vecfuncT & ao,
                                                  const std::vector<int> & set,
                                                  const std::vector<int> & at_to_bf,
                                                  const std::vector<int> & at_nbf,
                                                  const double thresh = 1e-9,
                                                  const double thetamax = 0.25,
                                                  const bool randomize = true,
                                                  const bool doprint = false)
{
    // Make Svec ... this can be much more efficient!
    tensorT S = matrix_inner(world, ao, ao, true);
    long nmo = mo.size();
    long nao = S.dim(0);
    long natom = at_to_bf.size();

    std::vector<tensorT> Svec(natom);
    for(long a = 0; a < natom; ++a){
        Slice as(at_to_bf[a], at_to_bf[a] + at_nbf[a] - 1);
        Svec[a] = copy(S(as, as));
    }
    S = tensorT();

    // Make initial matrices
    DistributedMatrix<double> dU = column_distributed_matrix<double>(world, nmo, nmo);
    dU.fill_identity();

    DistributedMatrix<double> dC = column_distributed_matrix<double>(world, nmo, nao);
    matrix_inner(dC, mo, ao);

    DistributedMatrix<double> dA = concatenate_rows(dC,dU);

    // Run the systolic algorithm
    world.taskq.add(new SystolicPMOrbitalLocalize(dA, set, at_to_bf, at_nbf, Svec, thresh, thetamax, natom, nao, nmo));
    world.taskq.fence();

    //print("DONE",world.rank());

    // Copy the data out
    Tensor<double> A(nmo, nao+nmo+natom);
    //    dA.copy_to_replicated(A);
    //U(___) = A(_,Slice(nao,nao+nmo-1));


    dA.extract_columns(nao,nmo+nao-1,dU);

    // Fix orbital orders in parallel
    world.taskq.add(new SystolicFixOrbitalOrders(dU));
    world.taskq.fence();

    return dU;

    //tensorT U(nmo, nmo);
    //dU.copy_to_replicated(U);
    //U = transpose(U);

    // if(world.rank() == 0){
    //     // Fix orbital orders
    //     bool switched = true;
    //     while (switched) {
    //         switched = false;
    //         for (int i=0; i<nmo; i++) {
    //     	for (int j=i+1; j<nmo; j++) {
    //                 if (set[i] == set[j]) {
    //                     double sold = U(i,i)*U(i,i) + U(j,j)*U(j,j);
    //                     double snew = U(i,j)*U(i,j) + U(j,i)*U(j,i);
    //                     if (snew > sold) {
    //                         tensorT tmp = copy(U(_,i));
    //                         U(_,i) = U(_,j);
    //                         U(_,j) = tmp;
    //                         switched = true;
    //                     }
    //                 }
    //     	}
    //         }
    //     }

    //     // Fix phases.
    //     for (long i=0; i<nmo; ++i) {
    //         if (U(i,i) < 0.0) U(_,i).scale(-1.0);
    //     }
    // }
    // world.gop.broadcast(U.ptr(), U.size(), 0);
    //return U;
}

}
