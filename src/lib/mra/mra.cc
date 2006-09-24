#include <iostream>
using std::cout;
using std::endl;

#include <cmath>
using std::abs;

#include <algorithm>
using std::max;
using std::sort;

#include <complex>
#include <cstring>
using std::strcpy;

#include "mra.h"
#include "twoscale.h"
#include "legendre.h"

#include <misc/madexcept.h>
#include <misc/shared_ptr.h>

#include <serialize/textfsar.h>
using madness::archive::TextFstreamInputArchive;
using madness::archive::TextFstreamOutputArchive;

#include <serialize/binfsar.h>
using madness::archive::BinaryFstreamInputArchive;
using madness::archive::BinaryFstreamOutputArchive;

#include <serialize/vecar.h>
using madness::archive::VectorInputArchive;
using madness::archive::VectorOutputArchive;

/// \file mra/mra.cc
/// \brief Implements non-task related stuff for Function and friends

namespace madness {

    // Definition and initialization of FunctionDefaults static members
    // It cannot be an instance of FunctionFactory since we want to
    // set the defaults independent of the data type.  The defaults
    // are only used in the FunctionFactory constructor except
    // for cell[][] which is used in FunctionData
    SharedPtr<FunctionOctTree> FunctionDefaults::tree(0);
    int FunctionDefaults::k = 7;
    double FunctionDefaults::thresh = 1e-5;
    int FunctionDefaults::initial_level = 2;
    int FunctionDefaults::max_refine_level = 30;
    int FunctionDefaults::truncate_method = 0;
    bool FunctionDefaults::compress = true;
    bool FunctionDefaults::refine = true;
    bool FunctionDefaults::autorefine = false;
    bool FunctionDefaults::debug = false;
    double FunctionDefaults::cell[3][2] = {{0.0,1.0}, {0.0,1.0}, {0.0,1.0}};
    
    template <typename T>
    FunctionCommonData<T> FunctionData<T>::commondata[FunctionData::MAXK+1];

    template <typename T>
    bool FunctionData<T>::initialized = false;

    template <typename T>
    void FunctionCommonData<T>::_make_dc_periodic() {
        // See ABGV for details
        r0 = Tensor<double>(k,k);
        rp = Tensor<double>(k,k);
        rm = Tensor<double>(k,k);

        double iphase = 1.0;
        for (int i=0; i<k; i++) {
            double jphase = 1.0;
            for (int j=0; j<k; j++) {
                double gammaij = sqrt(double((2*i+1)*(2*j+1)));
                double Kij;
                if (((i-j)>0) && (((i-j)%2)==1))
                    Kij = 2.0;
                else
                    Kij = 0.0;

                r0(i,j) = 0.5*(1.0 - iphase*jphase - 2.0*Kij)*gammaij;
                rm(i,j) = 0.5*jphase*gammaij;
                rp(i,j) =-0.5*iphase*gammaij;
                jphase = -jphase;
            }
            iphase = -iphase;
        }

        // Make the rank-1 forms of rm and rp
        rm_left = Tensor<double>(k);
        rm_right = Tensor<double>(k);
        rp_left = Tensor<double>(k);
        rp_right = Tensor<double>(k);

        iphase = 1.0;
        for (int i=0; i<k; i++) {
            double gamma = sqrt(0.5*(2*i+1));
            rm_left(i)  = rp_right(i) = gamma;
            rm_right(i) = rp_left(i)  = gamma*iphase;
            iphase *= -1.0;
        }
        rp_left.scale(-1.0);

//         Tensor<double> rm_test = outer(rm_left,rm_right);
//         Tensor<double> rp_test = outer(rp_left,rp_right);
    }

    template <typename T>
    void FunctionCommonData<T>::_init_twoscale() {
        if (! two_scale_hg(k, &hg)) throw "failed to get twoscale coefficients";
        hgT = transpose(hg);
        hgsonly = copy(hg(Slice(0,k-1),_));
    }

    template <typename T>
    void FunctionCommonData<T>::_init_quadrature() {
        quad_x = Tensor<double>(npt);
        quad_w = Tensor<double>(npt);
        quad_phi = Tensor<double>(npt,k);
        quad_phiw = Tensor<double>(npt,k);

        gauss_legendre(npt,0.0,1.0,quad_x.ptr(),quad_w.ptr());
        for (int i=0; i<npt; i++) {
            double phi[200];
            legendre_scaling_functions(quad_x(i),k,phi);
            for (int j=0; j<k; j++) {
                quad_phi(i,j) = phi[j];
                quad_phiw(i,j) = quad_w(i)*phi[j];
            }
        }
        quad_phit = transpose(quad_phi);
    }

    template <typename T>
    void Function<T>::_init(const FunctionFactory<T>& factory) {
        bool compress = factory._compress;
        bool empty = factory._empty;

        if (data->f || data->vf) {
            project_refine();
            if (compress) this->compress();
        } else if (empty) {             // Do not set any coefficients at all
            data->compressed = compress;
        } else {                        // Set coefficients so have zero function
            data->compressed = compress;
            if (tree()->n() == 0) {
                if (compress) 
                    set_coeff(tree(),TensorT(2*k,2*k,2*k));
                else    
                    set_coeff(tree(),TensorT(k,k,k));
            }
        }

        //tree()->print();
    }

    /// Private:  Projects onto wavelets at this level with optional refinement
    
    /// If refinement not necessary, sets scaling coeffs at level below
    /// and returns false.  Otherwise, returns true. 
    template <typename T>
    bool Function<T>::_doproject_refine(OctTreeTPtr& tree) {
        //madness::print("in doproject_refine");
        int npt = data->cdata->npt;
        Tensor<T> fval(npt,npt,npt);
        if (!data->refine) {
            Level n = tree->n();
            Translation x=tree->x(), y=tree->y(), z=tree->z();
            if (data->vf)
                  _vfcube(n, x, y, z, data->vf, fval);
               else
                  _fcube(n, x, y, z, data->f, fval);
            fval.scale(pow(0.125,0.5*n));
            set_coeff(tree,transform3d(fval,data->cdata->quad_phiw));
            return false;
        }
        else {             
            Tensor<T> r(2*k,2*k,2*k);
            Level n = tree->n()+1;
            Translation x=2*tree->x(), y=2*tree->y(), z=2*tree->z();
            const Slice* s = data->cdata->s;
            FORIJK(if (data->vf)
                      _vfcube(n, x+i, y+j, z+k, data->vf, fval);
                   else
                      _fcube(n, x+i, y+j, z+k, data->f, fval);
                   r(s[i],s[j],s[k]) = transform3d(fval,data->cdata->quad_phiw););
            r.scale(pow(0.125,0.5*n));
            
            filter_inplace(r);
            Tensor<T> ss = madness::copy(r(s[0],s[0],s[0]));
            r(s[0],s[0],s[0]) = 0.0;
            if (r.normf() > truncate_tol(data->thresh,tree->n())) return true;
            set_coeff(tree,ss);
            return false;
        }         
    }

    template <typename T>
    double Function<T>::_norm2sq_local(const OctTreeTPtr& tree) const {
        double sum = 0.0;
        if (coeff(tree)) {
            sum = coeff(tree)->normf();
            sum *= sum;
        }
        FOREACH_CHILD(const OctTreeTPtr, tree, 
                      if (isactive(child)) sum += _norm2sq_local(child););
  
        return sum;
    };


    template <typename T>
    void Function<T>::_add_sub(Function<T>& cfun,
                               const Function<T>& afun,
                               const Function<T>& bfun,
                               OctTreeTPtr& tree, bool subtract) const {
        // c = a + b,  (a == *this)
        const TensorT *bt=bfun.coeff(tree), *at=afun.coeff(tree);

        cfun.set_active(tree);
        if (bt) {
            if (at)
                if (subtract)
                    cfun.set_coeff(tree,*at-*bt);
                else
                    cfun.set_coeff(tree,*at+*bt);
            else
                if (subtract)
                    cfun.set_coeff(tree,-*bt);
                else
                    cfun.set_coeff(tree,madness::copy(*bt));

        } else if (at) {
            cfun.set_coeff(tree,madness::copy(*at));
        }

        FOREACH_CHILD(OctTreeTPtr, tree,
                      if (afun.isactive(child) || bfun.isactive(child))
                      _add_sub(cfun,afun,bfun,child,subtract););
    }



    template <typename T>
    void Function<T>::_gaxpy(Function<T>& afun, double alpha,
                             const Function<T>& bfun, double beta,
                             OctTreeTPtr& tree) {
        // a <= a*alpha + b*beta
        //
        // This routine will work if you are in either the scaling
        // function or wavelet basis.  However, in the scaling
        // function basis, you've got to remember to compress
        // afterwards in order to sum up the scaling function
        // coefficients at multiple levels.

        const TensorT *bt=bfun.coeff(tree);
        TensorT *at=afun.coeff(tree);

        afun.set_active(tree);
        if (bt) {
            if (at)
                at->gaxpy(alpha,*bt,beta);
            else
                afun.set_coeff(tree,(*bt)*beta);
        } else if (at) {
            if (alpha != 1.0) at->scale(alpha);
        }

        FOREACH_CHILD(OctTreeTPtr, tree,
                      if (afun.isactive(child) || bfun.isactive(child))
                      _gaxpy(afun,alpha,bfun,beta,child););
    }



    template <typename T>
    void Function<T>::_fcube(long n, long lx, long ly, long lz,
                             T (*f)(double,double,double), Tensor<T>& fcube) {
        double h = 1.0/two_to_power(n);
        double xlo = lx*h, ylo = ly*h, zlo = lz*h;

        int npt = data->cdata->npt;
        const Tensor<double>& quad_x = data->cdata->quad_x;

        for (int i=0; i<npt; i++) {
            double x = xlo + h*quad_x(i);
            for (int j=0; j<npt; j++) {
                double y = ylo + h*quad_x(j);
                for (int k=0; k<npt; k++) {
                    double z = zlo + h*quad_x(k);
                    fcube(i,j,k) = f(x,y,z);
                }
            }
        }
    }

    template <typename T>
    void Function<T>::_vfcube(long n, long lx, long ly, long lz,
                              void (*vf)(long l, const double*, const double*, const double*,
                                         T* RESTRICT),
                              Tensor<T>& fcube) {
        double h = 1.0/two_to_power(n);
        double xlo = lx*h, ylo = ly*h, zlo = lz*h;
        int npt = data->cdata->npt;
        const Tensor<double>& quad_x = data->cdata->quad_x;
        double x[64], y[64], z[64];
        MADNESS_ASSERT(npt<=64);

        int count = 0;
        for (int i=0; i<npt; i++) {
            double xx = xlo + h*quad_x(i);
            for (int j=0; j<npt; j++) {
                double yy = ylo + h*quad_x(j);
                for (int k=0; k<npt; k++) {
                    double zz = zlo + h*quad_x(k);
                    x[count] = xx;
                    y[count] = yy;
                    z[count] = zz;
                    count++;
                }
            }
        }

        vf(npt*npt*npt, x, y, z, fcube.ptr());
    }

    template <typename T>
    T Function<T>::_eval_cube(Level n,
                              double xx, double yy, double zz,
                              const Tensor<T>& s) const {

        double px[64], py[64], pz[64];
        MADNESS_ASSERT(k<64);
        legendre_scaling_functions(xx,k,px);
        legendre_scaling_functions(yy,k,py);
        legendre_scaling_functions(zz,k,pz);
        T sum = T(0.0);
        for (int p=0; p<k; p++) {
            for (int q=0; q<k; q++) {
                for (int r=0; r<k; r++) {
                    sum += s(p,q,r)*px[p]*py[q]*pz[r];
                }
            }
        }
        return sum*pow(2.0,1.5*n);
    }


    template <typename T>
    T Function<T>::_inner_local(const Function<T>& a, const Function<T>& b, const OctTreeTPtr& tree) const {
      T sum = 0.0;
      FOREACH_CHILD(OctTreeTPtr, tree, 
        if (a.isactive(child) && b.isactive(child)) sum += _inner_local(a, b, child);
      );
      if (a.coeff(tree) && b.coeff(tree)) sum += a.coeff(tree)->trace(*b.coeff(tree));

      return sum;
    }

    template <typename T>
    void Function<T>::_tnorms(const Tensor<T>& t, double& lo, double& hi) const {
        // t is a k*k*k tensor.  In order to screen the autorefinement
        // during squaring or multiplication, compute the norms of
        // ... lo ... the block of t for all polynomials of order < k/2
        // ... hi ... the block of t for all polynomials of order >= k/2

        // This routine exists to provide higher precision and in order to
        // optimize away the tensor slicing otherwise necessary to compute
        // the norm using
//        Slice s0(0,(k-1)/2);
//        double anorm = t.normf();
//        lo = t(s0,s0,s0).normf();
//        hi = sqrt(anorm*anorm - *lo**lo);

        // k=5   0,1,2,3,4     --> 0,1,2 ... 3,4
        // k=6   0,1,2,3,4,5   --> 0,1,2 ... 3,4,5

        // k=number of wavelets, so k=5 means max order is 4, so max exactly
        // representable squarable polynomial is of order 2 which is the third
        // polynomial.
        int k2 = (k - 1) / 2 + 1;
        double slo = 0.0, shi = 0.0;

        for (int p = 0; p < k2; p++)
            for (int q = 0; q < k2; q++)
                for (int r = 0; r < k2; r++)
                    slo += std::norm(t(p, q, r));

        for (int p = 0; p < k; p++)
            for (int q = k2; q < k; q++)
                for (int r = 0; r < k; r++)
                    shi += std::norm(t(p, q, r));

        for (int p = k2; p < k; p++)
            for (int q = 0; q < k2; q++)
                for (int r = 0; r < k; r++)
                    shi += std::norm(t(p, q, r));

        for (int p = 0; p < k2; p++)
            for (int q = 0; q < k2; q++)
                for (int r = k2; r < k; r++)
                    shi += std::norm(t(p, q, r));

        lo = sqrt(slo);
        hi = sqrt(shi);
    }
    
    template <typename T>
    int Function<T>::unique_active_child_procs(const OctTreeTPtr& tree, ProcessID ranks[8]) const {
        MADNESS_ASSERT(tree);
        int np = 0;
        FOREACH_REMOTE_CHILD(const OctTreeTPtr, tree,
            if (isactive(child)) {
                bool found = false;
                for (int p=0; p<np; p++) {
                    if (child->rank() == ranks[p]) {
                        found = true;
                        break;
                    }
                }
                if (!found) ranks[np++] = child->rank();
            }
        );
        return np;
    };
    
    template <typename T>
    void Function<T>::do_transform_to_function_values(OctTreeTPtr& tree) {
        Tensor<T>& t = *coeff(tree);
        double scale = std::pow(8.0, 0.5*tree->n());
        transform3d_inplace(t, data->cdata->quad_phit, data->cdata->work1);
        t.scale(scale);
    };
    

    template <typename T>
    void Function<T>::do_transform_from_function_values(OctTreeTPtr& tree) {
        transform3d_inplace(*coeff(tree), data->cdata->quad_phiw, data->cdata->work1);
    };

    
    template <typename T>
    void Function<T>::do_square(OctTreeTPtr& tree) {
        MADNESS_ASSERT(coeff(tree));
        Tensor<T>& t = *coeff(tree);
        double scale = std::pow(8.0, 0.5 * tree->n());
        transform3d_inplace(t, data->cdata->quad_phit, data->cdata->work1);
        t.emul(t).scale(scale);
        transform3d_inplace(t, data->cdata->quad_phiw, data->cdata->work1);
    };
    

    template <typename T>
    bool Function<T>::eval_local(double x, double y, double z, T* value) const {
        if (iscompressed())
            MADNESS_EXCEPTION("Function: eval_local: must not be compressed",ind);

        *value = T(0.0);
        OctTreeTPtr tree = this->tree();
        if (isactive(tree)) {
            // Determine if point might be local to this process
            data->cdata->user_to_sim(x, y, z);
            double twon = two_to_power(tree->n());
            x *= twon;
            y *= twon;
            z *= twon;
            x -= tree->x();
            y -= tree->y();
            z -= tree->z();

            if (x < 0.0 || y < 0.0 || z < 0.0 ||
                    x > 1.0 || y > 1.0 || z > 1.0) {
                if (tree->n() == 0) {
                    MADNESS_EXCEPTION("Function: eval_local: out of range point",ind);
                } else {
                    return false;
                }
            }

            // Recur down to find it ... inline for maximum efficiency
            while (tree && isactive(tree)) {
               const TensorT* t = coeff(tree);
                if (t) {
                    *value = _eval_cube(tree->n(), x, y, z, *t);
                     return true;
                }
                x *= 2.0;
                y *= 2.0;
                z *= 2.0;
                int lx = int(x);
                int ly = int(y);
                int lz = int(z);
                if (lx == 2) --lx;
                if (ly == 2) --ly;
                if (lz == 2) --lz;
                x -= lx;
                y -= ly;
                z -= lz; 
                tree = tree->child(lx, ly, lz);
            }
        }
        return false;
    };
    

    // Explicit instantiations for double and complex<double>

    template class Function<double>;
    template class Function< std::complex<double> >;
    template class FunctionData<double>;
    template class FunctionData< std::complex<double> >;
    template class FunctionCommonData<double>;
    template class FunctionCommonData<double_complex>;
    
//    template void Function<double>::_compress(SharedPtr< OctTree<madness::FunctionNode> >&);
//    template void Function<double_complex>::_compress(SharedPtr< OctTree<madness::FunctionNode> >&);
    
}
