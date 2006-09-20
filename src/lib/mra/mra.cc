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
/// \brief Implements Function and friends

/// Eventually, these notes provide a basic user's
/// programmer's manual.

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
        bool refine = factory._refine;
        bool empty = factory._empty;

        if (data->f || data->vf) {
            _fine_scale_projection(tree(), data->initial_level);
            if (refine) _refine(tree());
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

    template <typename T>
    void Function<T>::_fine_scale_projection(OctTreeTPtr& tree, Level initial_level) {
        // Recur down oct_tree to initial_level.  Refine locally
        // if necessary.  Project when get to desired level
        if (tree->n() < initial_level) {
            set_active(tree);
            FORIJK(OctTreeTPtr c = tree->child(i,j,k);
                   if (!c && islocal(tree)) c = tree->insert_local_child(i,j,k);
                   if (c) _fine_scale_projection(c,initial_level););
        } else if (tree->n()==initial_level) {
            set_active(tree);
            if (islocal(tree)) _project(tree);
        }
    }


    /// Private:  Projects function in given box
    template <typename T>
    void Function<T>::_project(OctTreeTPtr& tree) {
        // Evaluate all local child scaling coeffs of the current box
        int npt = data->cdata->npt;
        Tensor<T> r(2*k,2*k,2*k);
        Tensor<T> fval(npt,npt,npt);
        Level n = tree->n()+1;
        Translation x=2*tree->x(), y=2*tree->y(), z=2*tree->z();
        const Slice* s = data->cdata->s;
        for (int ix=0; ix<2; ix++) {
            for (int iy=0; iy<2; iy++) {
                for (int iz=0; iz<2; iz++) {
                    OctTreePtr child = tree->child(ix,iy,iz);
                    if (!child) child = tree->insert_local_child(ix,iy,iz);
                    if (data->vf) {
                        _vfcube(n, x+ix, y+iy, z+iz, data->vf, fval);
                    } else if (data->f) {
                        _fcube(n, x+ix, y+iy, z+iz, data->f, fval);
                    }
                    // Can optimize away the tensor construction in transform3d
                    r(s[ix],s[iy],s[iz]) = transform3d(fval,data->cdata->quad_phiw);
                }
            }
        }
        set_coeff(tree,r.scale(pow(0.125,0.5*n)));
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
    void Function<T>::_reconstruct(OctTreeTPtr& tree) {
        std::vector<Slice>& s0 = data->cdata->s0;
        Slice* s = data->cdata->s;
        TensorT& buf = data->cdata->work1;
        long k3 = k*k*k;

        if (isremote(tree)) {
            if (tree->islocalsubtreeparent()) {
                FOREACH_CHILD(OctTreeTPtr, tree,
                              if (isactive(child)) {
                              comm()->Recv(buf.ptr(),k3,tree->rank(),RECONSTRUCT_TAG);
              //                                  print("reconstruct: received",buf.normf(),"from",tree->rank());
                                  (*coeff(child))(s0) = buf;
                              }
                             );
            } else {
                throw "_reconstruct: should not be here?";
            }
        } else {
            int nchild = 0;
            TensorT* t = coeff(tree);
            unfilter_inplace(*t);
            FOREACH_CHILD(OctTreeTPtr, tree,
                          if (isactive(child)) {
                          nchild++;
                          if (islocal(child)) {
                                  (*coeff(child))(s0) = (*t)(s[i],s[j],s[k]);
                              } else {
                                  buf(s0) = (*t)(s[i],s[j],s[k]);
              //                                  print("reconstruct: sending",buf.normf(),"to",child->rank());
                                  comm()->Send(buf.ptr(),k3,child->rank(),RECONSTRUCT_TAG);
                              }
                          }
                         );
            if (nchild) unset_coeff(tree); // NEED TO KEEP IF WANT REDUNDANT TREE
        }

        FOREACH_CHILD(OctTreeTPtr, tree,
                      if (isactive(child) && islocal(child)) _reconstruct(child););
    }

    /// Private: Refine an existing block of coefficients
    template <typename T>
    void Function<T>::_refine(OctTreeTPtr& tree) {
        if (tree->islocalsubtreeparent() && isremote(tree)) {
            // This remote node is the parent of the local subtree.
            // Loop thru the local children and receive refinement
            // info from above.  The messages will come in order.
            FOREACH_CHILD(OctTreeTPtr, tree,
                          bool dorefine;
                          comm()->Recv(dorefine, tree->rank(), REFINE_TAG);
              //                          print("   refine: got message",dorefine,"for",child->n(),child->x(),child->y(),child->z());
                          if (dorefine) {
                          set_active(child);
                              set_active(tree);
                              _project(child);
                          }
                          _refine(child);
                         );
        } else {
            TensorT* t = coeff(tree);
            if (t) {
                // Local node with data
                if (isremote(tree)) throw "refine: remote node with data?";
                TensorT d = filter(*t);
                d(data->cdata->s0) = 0.0;
                bool dorefine = (d.normf() > truncate_tol(data->thresh,tree->n()));
//                if (dorefine)
//                    print("refine:",tree->n(),tree->x(),tree->y(),tree->z(),d.normf(),"dorefine =",dorefine);
                    if (dorefine) unset_coeff(tree);
                // First, send messages to remote children in order to get them working ASAP
                FOREACH_REMOTE_CHILD(OctTreeTPtr, tree,
                                     comm()->Send(dorefine, child->rank(), REFINE_TAG);
                                     set_active(child);
                     //                                     print("sent",dorefine,"message to",child->n(),child->x(),child->y(),child->z());
                                    );
                // Next, handle local stuff
                FORIJK(OctTreeTPtr child = tree->child(i,j,k);
                       // If child does not exist, OK to refine locally
                       if (!child && dorefine) child = tree->insert_local_child(i,j,k);
                       if ( child && islocal(child)) {
                           if (dorefine) {
                                   _project(child);
                                   set_active(child);
                               }
                               _refine(child);
                           }
                      );
            } else {
                // This node does not have data and is not a remote subtree parent.
                // It must be either an empty local node or a remote leaf.
                // Again, first send msgs to remote nodes, then treat local guys.
                FOREACH_REMOTE_CHILD(OctTreeTPtr, tree,
                                     comm()->Send(false, child->rank(), REFINE_TAG);
                     //                                     print("    sent false to",child->n(),child->x(),child->y(),child->z());
                                    );
                FOREACH_LOCAL_CHILD(OctTreeTPtr, tree, _refine(child););
            }
        }
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
#ifdef IBMXLC
        if (npt > 35) throw("hard dimension failure");
#endif
        double* x = new double[npt*npt*npt];
        double* y = new double[npt*npt*npt];
        double* z = new double[npt*npt*npt];

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

        delete[] x;
        delete[] y;
        delete[] z;
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
    void Function<T>::_tnorms(const Tensor<T>& t, double *lo, double *hi) {
        // t is a k*k*k tensor.  In order to screen the autorefinement
        // during squaring or multiplication, compute the norms of
        // ... lo ... the block of t for all polynomials of order < k/2
        // ... hi ... the block of t for all polynomials of order >= k/2

        // This routine exists to provide higher precision and in order to
        // optimize away the tensor slicing otherwise necessary to compute
        // the norm using
        Slice s0(0,(k-1)/2);
        double anorm = t.normf();
        *lo = t(s0,s0,s0).normf();
        *hi = sqrt(anorm*anorm - *lo**lo);
/*
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

            *lo = sqrt(slo);
            *hi = sqrt(shi);*/
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
    void Function<T>::_autorefine(OctTreeTPtr& tree) {
        MADNESS_ASSERT(tree);
        bool msg[2], refine=false;
        const Slice* s = data->cdata->s;
        const Slice& s0 = s[0];
        long k2 = k*2;
        Tensor<T>& work1 = data->cdata->work1;
        
        if (islocal(tree)) {
            ProcessID ranks[8];
            int np = tree->unique_child_procs(ranks);
            if (coeff(tree)) { // Test for refinement
                double tol = truncate_tol(data->thresh,tree->n()+1);
                const Tensor<T>& c = *coeff(tree);
                double lo,hi;
                FORIJK(_tnorms(c(s[i],s[j],s[k]),&lo,&hi);
                       refine = hi > tol;
                       if (refine) goto done;);
                done:;
            }

            // Tell remote clones what is going on
            msg[0] = isactive(tree); 
            msg[1] = refine;
            
            for (int i=0; i<np; i++) comm()->Send(msg,2,ranks[i],AUTOREF_TAG1);
            
            if (refine) { // refine, sending data as necessary;
                const Tensor<T>& c = *coeff(tree);
                FORIJK(OctTreeTPtr child = tree->child(i,j,k);
                       if (!child) child = tree->insert_local_child(i,j,k);
                       set_active(child);
                       if (child->islocal()) {
                          Tensor<T>*t = set_coeff(child, Tensor<T>(k2, k2, k2));
                          (*t)(s0, s0, s0) = c(s[i],s[j],s[k]);
                          unfilter_inplace(*t);  // sonly!
                       }
                       else {
                          work1(s0,s0,s0) = c(s[i],s[j],s[k]); // contig. copy
                          comm()->Send(work1.ptr(), work1.size, child->rank(), AUTOREF_TAG2);
                       }
                    );
                unset_coeff(tree);               
            }
        } else if (isremote(tree) && tree->islocalsubtreeparent()) {
            comm()->Recv(msg,2,tree->rank(),AUTOREF_TAG1);
            bool active=msg[0], refine=msg[1];
            if (!active && isactive(tree)) {
                MADNESS_EXCEPTION("Remote clone thinks it is inactive but I think I am active",0);
            } else if (active) {
                set_active(tree);
                if (refine) {
                    FORIJK(OctTreeTPtr child = tree->child(i,j,k);
                          if (!child) child = tree->insert_local_child(i,j,k);
                          set_active(child);
                          comm()->Recv(work1.ptr(), work1.size, tree->rank(), AUTOREF_TAG2);
                          Tensor<T>*c = set_coeff(child, Tensor<T>(k2, k2, k2));
                          (*c)(s0, s0, s0) = work1;
                          unfilter_inplace(*c); // sonly needed!
                    );
                }
            }
        }            
        FOREACH_CHILD(OctTreeTPtr, tree, if (islocal(child)) _autorefine(child););    
    };

    template <typename T>
    void Function<T>::_square(OctTreeTPtr& tree) {
         MADNESS_ASSERT(tree);
         FOREACH_ACTIVE_CHILD(OctTreeTPtr, tree, _square(child););

         if (coeff(tree)) {
            Tensor<T>& t = *coeff(tree);
            Tensor<T> r(k, k, k);
            const Slice* s = data->cdata->s;
            double scale = std::pow(8.0, 0.5 * (tree->n()+1));
            FORIJK(r(___) = t(s[i], s[j], s[k]);
                   transform3d_inplace(r, data->cdata->quad_phit,
                                       data->cdata->work1);
                   r.emul(r).scale(scale);
                   transform3d_inplace(r, data->cdata->quad_phiw,
                                       data->cdata->work1);
                   t(s[i], s[j], s[k]) = r;
                  );
         }
    };
    


    // Explicit instantiations for double and complex<double>

    template class Function<double>;
    template class Function< std::complex<double> >;
    template class FunctionData<double>;
    template class FunctionData< std::complex<double> >;
    template class FunctionCommonData<double>;
    template class FunctionCommonData<double_complex>;

}
