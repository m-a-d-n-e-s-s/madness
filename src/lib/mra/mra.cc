#include <iostream>
using std::cout;
using std::endl;

#include <cmath>
using std::abs;

#include <algorithm>
using std::max;
using std::sort;

#include <complex>

#include "mra.h"
#include "twoscale.h"
#include "legendre.h"

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
    SharedPtr<FunctionOctTree> FunctionDefaults::tree = 0;
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
        } else {                        // Set coefficients so have zero function
            if (tree()->n() == 0) set_coeff(tree(),TensorT(2*k,2*k,2*k));
        }
        data->compressed = compress; // Just to be sure we are consistent.

        //tree()->print();
    }

    template <typename T>
    void Function<T>::_fine_scale_projection(OctTreeT* tree, Level initial_level) {
        // Recur down oct_tree to initial_level.  Refine locally
        // if necessary.  Project when get to desired level
        if (tree->n() < initial_level) {
            set_active(tree);
            FORIJK(OctTreeT* c = tree->child(i,j,k);
                   if (!c && islocal(tree)) c = tree->insert_local_child(i,j,k);
                   if (c) _fine_scale_projection(c,initial_level););
        } else if (tree->n()==initial_level) {
            set_active(tree);
            if (islocal(tree)) _project(tree);
        }
    }


    /// Private:  Projects function in given box
    template <typename T>
    void Function<T>::_project(OctTreeT* tree) {
        // We are at level n in the of wavelet coefficients, which
        // corresponds to level n+1 in the scaling function
        // coefficients.  Evaluate all 8 blocks of the scaling
        // function coefficients which are stored together in the tree.
        int npt = data->cdata->npt;
        Tensor<T> r(2*k,2*k,2*k);
        Tensor<T> fval(npt,npt,npt);
        Level n = tree->n();
        Translation x=2*tree->x(), y=2*tree->y(), z=2*tree->z();
        const Slice* s = data->cdata->s;
        for (int ix=0; ix<2; ix++) {
            for (int iy=0; iy<2; iy++) {
                for (int iz=0; iz<2; iz++) {
                    if (data->vf) {
                        _vfcube(n+1, x+ix, y+iy, z+iz, data->vf, fval);
                    } else if (data->f) {
                        _fcube(n+1, x+ix, y+iy, z+iz, data->f, fval);
                    }
                    // Can optimize away the tensor construction in transform3d
                    r(s[ix],s[iy],s[iz]) = transform3d(fval,data->cdata->quad_phiw);
                }
            }
        }
        set_coeff(tree,r.scale(pow(0.125,0.5*n)));
    }

    template <typename T>
    double Function<T>::_norm2sq_local(const OctTreeT* tree) const {
        if (! isactive(tree)) return 0.0;
        double sum = 0.0;
        const TensorT* t = coeff(tree);
        if (t) {
            sum = t->normf();
            sum *= sum;
        }
        FOREACH_CHILD(const OctTreeT, tree,
                      if (isactive(child)) sum += _norm2sq_local(child););
        //sum += _norm2sq_local(child););
        return sum;
    };

    template <typename T>
    double Function<T>::_norm2sq(const OctTreeT* tree) const {
        double sum=0.0;
        FOREACH_CHILD(OctTreeT, tree,
                      if (isactive(child)) {
                      double tmpsum = _norm2sq(child);
                          sum += tmpsum;
                      }
                     );
        if ( isremote(tree) ) {
            if( (tree->parent()) && (tree->parent()->islocal()) ) {
                comm()->Recv(&sum, 1, tree->rank(), 4);
            }
	}
        else {
          const TensorT* t = coeff(tree);
          if (t) {
            double tsum = t->normf();
            sum += tsum*tsum;
          }
	  if ((tree->parent()) && (tree->parent()->isremote())) {
            comm()->Send(&sum, 1, tree->parent()->rank(), 4);
	  }
        }
        return sum;
    };

    template <typename T>
    void Function<T>::_compress(OctTreeT* tree) {
        FOREACH_CHILD(OctTreeT, tree,
                      if (isactive(child)) _compress(child););

        if (isremote(tree)) {
            long k3 = k*k*k;
            if (tree->islocalsubtreeparent()) {
                // Send data from active children to remote parent
                std::vector<Slice>& s0 = data->cdata->s0;
                TensorT& buf = data->cdata->work1;
                FOREACH_CHILD(OctTreeT, tree,
                              if (isactive(child)) {
                              const TensorT* c = coeff(child);
                                  if (!c) throw "compress: yo! show me the data(2)!";
                                  buf(s0) = (*c)(s0);
              //                                  print("compress: sending",buf.normf(),"to",tree->rank());
                                  comm()->Send(buf.ptr(), k3, tree->rank(), 2);
              //                                  print("child norm before",c->normf());
                                  (*c)(s0) = T(0.0);
              //                                  print("child norm after",c->normf());
                              }
                             );
            } else {
                // Receive data from a remote child and store as
                // ghost ... parent will eventually delete it.
                TensorT* t = set_coeff(tree,TensorT(k,k,k));
                comm()->Recv(t->ptr(), k3, tree->rank(), 2);
//                print("compress: received",t->normf(),"from",tree->rank());
            }
        } else {
            // Node is local.  May or may not have data.
            Slice *s = data->cdata->s;
            TensorT* t = coeff(tree);
            if (!t) t = set_coeff(tree,TensorT(2*k,2*k,2*k));
            FOREACH_CHILD(OctTreeT, tree,
                          if (isactive(child)) {
                          TensorT* c = coeff(child);
                              if (!c) throw "compress: yo! show me the data!";
                              (*t)(s[i],s[j],s[k]) += (*c)(s[0],s[0],s[0]);
                              if (isremote(child)) {
                                  unset_coeff(child);
                              } else {
                                  (*c)(s[0],s[0],s[0]) = T(0.0);
                              }
                          }
                         );
            //set_coeff(tree,filter(*t));
            filter_inplace(*t);
        }
    }

    template <typename T>
    void Function<T>::_add_sub(Function<T>& cfun,
                               const Function<T>& afun,
                               const Function<T>& bfun,
                               OctTreeT* tree, bool subtract) const {
        // c = a + b,  (a == *this)
        TensorT *bt=bfun.coeff(tree), *at=afun.coeff(tree);

        cfun.set_active(tree);
        if (bt) {
            if (at)
                if (subtract)
                    cfun.set_coeff(tree,*at-*bt);
                else
                    cfun.set_coeff(tree,*at+*bt);
            else
                if (subtract)
                    cfun.set_coeff(tree,madness::copy(-*bt));
                else
                    cfun.set_coeff(tree,madness::copy(*bt));

        } else if (at) {
            cfun.set_coeff(tree,madness::copy(*at));
        }

        FOREACH_CHILD(OctTreeT, tree,
                      if (afun.isactive(child) || bfun.isactive(child))
                      _add_sub(cfun,afun,bfun,child,subtract););
    }



    template <typename T>
    void Function<T>::_gaxpy(Function<T>& afun, double alpha,
                             const Function<T>& bfun, double beta,
                             OctTreeT* tree) {
        // a <= a*alpha + b*beta
        //
        // This routine will work if you are in either the scaling
        // function or wavelet basis.  However, in the scaling
        // function basis, you've got to remember to compress
        // afterwards in order to sum up the scaling function
        // coefficients at multiple levels.

        TensorT *bt=bfun.coeff(tree), *at=afun.coeff(tree);

        afun.set_active(tree);
        if (bt) {
            if (at)
                at->gaxpy(alpha,*bt,beta);
            else
                afun.set_coeff(tree,(*bt)*beta);
        } else if (at) {
            if (alpha != 1.0) at->scale(alpha);
        }

        FOREACH_CHILD(OctTreeT, tree,
                      if (afun.isactive(child) || bfun.isactive(child))
                      _gaxpy(afun,alpha,bfun,beta,child););
    }


    template <typename T>
    void Function<T>::_reconstruct(OctTreeT* tree) {
        std::vector<Slice>& s0 = data->cdata->s0;
        Slice* s = data->cdata->s;
        TensorT& buf = data->cdata->work1;
        long k3 = k*k*k;

        if (isremote(tree)) {
            if (tree->islocalsubtreeparent()) {
                FOREACH_CHILD(OctTreeT, tree,
                              if (isactive(child)) {
                              comm()->Recv(buf.ptr(),k3,tree->rank(),3);
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
            FOREACH_CHILD(OctTreeT, tree,
                          if (isactive(child)) {
                          nchild++;
                          if (islocal(child)) {
                                  (*coeff(child))(s0) = (*t)(s[i],s[j],s[k]);
                              } else {
                                  buf(s0) = (*t)(s[i],s[j],s[k]);
              //                                  print("reconstruct: sending",buf.normf(),"to",child->rank());
                                  comm()->Send(buf.ptr(),k3,child->rank(),3);
                              }
                          }
                         );
            if (nchild) unset_coeff(tree); // NEED TO KEEP IF WANT REDUNDANT TREE
        }

        FOREACH_CHILD(OctTreeT, tree,
                      if (isactive(child) && islocal(child)) _reconstruct(child););
    }

    /// Private: Refine an existing block of coefficients
    template <typename T>
    void Function<T>::_refine(OctTreeT* tree) {
        if (tree->islocalsubtreeparent() && isremote(tree)) {
            // This remote node is the parent of the local subtree.
            // Loop thru the local children and receive refinement
            // info from above.  The messages will come in order.
            FOREACH_CHILD(OctTreeT, tree,
                          bool dorefine;
                          comm()->Recv(dorefine, tree->rank(), 1);
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
                if (dorefine)
//                    print("refine:",tree->n(),tree->x(),tree->y(),tree->z(),d.normf(),"dorefine =",dorefine);
                    if (dorefine) unset_coeff(tree);
                // First, send messages to remote children in order to get them working ASAP
                FOREACH_REMOTE_CHILD(OctTreeT, tree,
                                     comm()->Send(dorefine, child->rank(), 1);
                                     set_active(child);
                     //                                     print("sent",dorefine,"message to",child->n(),child->x(),child->y(),child->z());
                                    );
                // Next, handle local stuff
                FORIJK(OctTreeT* child = tree->child(i,j,k);
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
                FOREACH_REMOTE_CHILD(OctTreeT, tree,
                                     comm()->Send(false, child->rank(), 1);
                     //                                     print("    sent false to",child->n(),child->x(),child->y(),child->z());
                                    );
                FOREACH_LOCAL_CHILD(OctTreeT, tree, _refine(child););
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
                              int lx, int ly, int lz,
                              const Tensor<T>& s) const {

        double px[64], py[64], pz[64];
        legendre_scaling_functions(xx,k,px);
        legendre_scaling_functions(yy,k,py);
        legendre_scaling_functions(zz,k,pz);

        lx *= k;
        ly *= k;
        lz *= k; // Offset into correct sub cube

        T sum = T(0.0);
        for (int p=0; p<k; p++) {
            for (int q=0; q<k; q++) {
                for (int r=0; r<k; r++) {
                    sum += s(p+lx,q+ly,r+lz)*px[p]*py[q]*pz[r];
                }
            }
        }
        return sum*pow(2.0,1.5*n);
    }

    template <typename T>
    void Function<T>::_truncate(double tol, OctTreeT *tree){
      FOREACH_CHILD(OctTreeT, tree, 
        if (isactive(child)) _truncate(tol, child);
      );
      const TensorT *t = coeff(tree);
      int nachildren = count_active_children(tree);
      bool normf_check = false;
      if(t) {
        double tnormf = t->normf();
        if(tnormf > tol) {
          normf_check = true;
        }
      }
      if ( !normf_check && nachildren == 0) {
        //unset_coeff(tree);
        if(!isremote(tree)) set_inactive(tree);
        if ((tree->parent()) && (tree->parent()->isremote())) {
          bool active;
          active = false;
          comm()->Send(&active, 1, tree->parent()->rank(), 5);
        }
      }
      else {
        set_active(tree);
        if ((tree->parent()) && (tree->parent()->isremote())) {
          bool active;
          active = true;
          comm()->Send(&active, 1, tree->parent()->rank(), 5);
        }
      }
    };

    template <typename T>
    int Function<T>::count_active_children(OctTreeT *tree) {
      int n = 0;
      bool active;
      FOREACH_LOCAL_CHILD(OctTreeT, tree, 
        if(isactive(child)) { n++; }
      );
      FOREACH_REMOTE_CHILD(OctTreeT, tree,
        comm()->Recv(&active, 1, child->rank(), 5);
        if (active) { 
          n++;
	}
	else {
          set_inactive(child);
	}
      );
      return n;
    }

    template <typename T>
    T Function<T>::_inner(const Function<T>& a, const Function<T>& b, OctTreeT* tree) const {
      T sum = 0.0;
      FOREACH_CHILD(OctTreeT, tree, 
        if (a.isactive(child) && b.isactive(child)) { 
          sum += _inner(a, b, child);
        }
      );
      if (isremote(tree)) {
        if( (tree->parent()) && (tree->parent()->islocal()) ) {
          comm()->Recv(&sum, 1, tree->rank(), 7);
        }
      }
      else {
        TensorT *c1 = a.coeff(tree);
        TensorT *c2 = b.coeff(tree);
        sum += c1->trace(*c2);
        if ( (tree->parent()) && (tree->parent()->isremote())) {
          comm()->Send(&sum, 1, tree->parent()->rank(), 7);
        }
      }
      return sum;
    }

    template <typename T>
    template <class Archive>
    void Function<T>::save_local(Archive& ar) {
      ar & FunctionDefaults::k;
      ar & FunctionDefaults::thresh;
      ar & FunctionDefaults::initial_level;
      ar & FunctionDefaults::max_refine_level;
      ar & FunctionDefaults::truncate_method;
      ar & FunctionDefaults::autorefine;
      _save_local(ar, tree());
    };

    template <typename T>
    template <class Archive>
    void Function<T>::_save_local(Archive& ar, const OctTreeT *tree) {
      ar & isactive(tree);
      ar & tree->n() & tree->x() & tree->y() & tree->z();
      if(isactive(tree)) {
        const TensorT *t = coeff(tree);
        ar & (t != 0);
        if(t) ar & *t;
      }
      FORIJK(
        ar & (tree->child(i,j,k)!=0);
        if (tree->child(i,j,k)) {
          _save_local(ar, tree->child(i,j,k));
        }
      );
      //}
    }

    template <typename T>
    template <class Archive>
    void Function<T>::load_local(const Archive& ar) {
      ar & FunctionDefaults::k;
      ar & FunctionDefaults::thresh;
      ar & FunctionDefaults::initial_level;
      ar & FunctionDefaults::max_refine_level;
      ar & FunctionDefaults::truncate_method;
      ar & FunctionDefaults::autorefine;
      bool active_flag;
      ar & active_flag;
      _load_local(ar, tree());
    }

    template <typename T>
    template <class Archive>
    void Function<T>::_load_local(const Archive& ar, OctTreeT *tree) {
      set_active(tree);
      Level n_local;
      Translation x_local, y_local, z_local;
      ar & n_local & x_local & y_local & z_local;
      bool have_coeffs;
      ar & have_coeffs;
      if(have_coeffs) {
        TensorT t;
        ar & t;
        set_coeff(tree, t);
      }
      FORIJK(
        bool have_child;
        ar & have_child;
        if (have_child) {
          OctTreeT* child = tree->child(i, j, k);
          if (!child) child = tree->insert_local_child(i,j,k);
          bool active_flag;
          ar & active_flag;
          if (active_flag) _load_local(ar, child);
        }
      );
    }

    template <typename T>
    void Function<T>::save(const char* f, Communicator& comm) {
      if (comm.rank() == 0) {
        TextFstreamOutputArchive oar(f);
        saveMain(oar, comm);
        oar.close();
      }
      else {
        saveLoadWorker(tree(), comm, true);
      }
    }

    template <typename T>
    template <class Archive>
    void Function<T>::saveMain(Archive& ar, Communicator& comm) {
      ar & FunctionDefaults::k;
      ar & FunctionDefaults::thresh;
      ar & FunctionDefaults::initial_level;
      ar & FunctionDefaults::max_refine_level;
      ar & FunctionDefaults::truncate_method;
      ar & FunctionDefaults::autorefine;
      saveManager(ar, tree(), comm);
      if( comm.size() != 1) {
        for ( int i = 1; i < comm.size(); i++) {
          archive::MPIOutputArchive arout(comm, i);
          arout & -1;
        }
      }
    };

    template <typename T>
    template <class Archive>
    void Function<T>::saveManager(Archive& ar, const OctTreeT *tree, Communicator& commFunc) {
      if (isremote(tree)) {
      shadowManager_save(ar, tree, tree->n(), tree->x(), tree->y(), 
               tree->z(), tree->rank(), commFunc);
      }
      else {
        ar & isactive(tree);
        ar & tree->n() & tree->x() & tree->y() & tree->z();
        cout << tree->n() << " " <<  tree->x() << " " <<  tree->y() << " " <<  tree->z() << endl;
        if(isactive(tree)) {
          const TensorT *t = coeff(tree);
          ar & (t != 0);
          if(t) ar & *t;
        }
        FORIJK( 
          OctTreeT *child = tree->child(i,j,k);
          if(child) {
            if(child->islocal()) {
              ar & (child!=0);
            }
            saveManager(ar, child, commFunc);
          }
          else {
            ar & (child!=0);
          }
        );
      }
    }

    template <typename T>
    template <class Archive>
    void Function<T>::shadowManager_save(Archive& ar, const OctTreeT *tree, Level n, Translation x, 
              Translation y, Translation z, 
              ProcessID remoteRank, Communicator& commFunc) 
    {
      int nRemoteBranch;
      madness::archive::MPIOutputArchive arout(commFunc, remoteRank);
      // Sending start branch data.
      arout & 1 & n & x & y & z;
      madness::archive::MPIInputArchive arin(commFunc, remoteRank);
      std::vector<localTreeMember> subtreeList; 
      // Getting local subtree list.
      arin & nRemoteBranch;
      subtreeList.resize(nRemoteBranch);
      for (int i = 0; i < nRemoteBranch; i++) {
        arin & subtreeList[i];
      }
      for (int i = 0; i < nRemoteBranch; i++) {
        ar & subtreeList[i].have_child;
        //cout << " local havechild " << i << "= " << subtreeList[i].have_child << endl;
        if(subtreeList[i].have_child) {
          if(subtreeList[i].remote) {
            if(subtreeList[i].rank==0) {
            // Going to other local tree on master.
              OctTreeT *treemaster = tree->find(subtreeList[i].n, 
               subtreeList[i].x, subtreeList[i].y, subtreeList[i].z);
               saveManager(ar, treemaster, commFunc);
            }
            else {
            // Going to other local tree on client.
              shadowManager_save(ar, tree, subtreeList[i].n, subtreeList[i].x, 
                subtreeList[i].y, subtreeList[i].z, 
                subtreeList[i].rank, commFunc);
            }
          }
          else {
            // Storing active flag.
            ar & subtreeList[i].active;
            ar & subtreeList[i].n & subtreeList[i].x & subtreeList[i].y & subtreeList[i].z;
            if(subtreeList[i].active) {
              // Inquiring Coefficients.
              bool coeffsPointer;
              commFunc.Recv(&coeffsPointer, 1, subtreeList[i].rank, 14);
              ar & (coeffsPointer != 0);
              // Storing Coefficients.
              if(coeffsPointer) {
                TensorT Coefficients(2*k, 2*k, 2*k);
                arin & Coefficients;
                ar & Coefficients;
              }
            }
          }
        }
      }
    }

    template <typename T>
    template <class Archive>
    void Function<T>::shadowManager_load(const Archive& ar, OctTreeT* tree, 
              Level n, 
              Translation x, Translation y, Translation z, 
              ProcessID remoteRank, Communicator& commFunc,
	      bool active_flag, bool have_child) 
    {
      //int nRemoteBranch;
      madness::archive::MPIOutputArchive arout(commFunc, remoteRank);
      // Send start branch data.
      arout & n & x & y & z;
      madness::archive::MPIInputArchive arin(commFunc, remoteRank);
      Level n_local;
      Translation x_local, y_local, z_local;
      ProcessID local_remoteRank;
      //bool remoteFlag;
      arout & have_child & active_flag;
      ar & n_local & x_local & y_local & z_local;
      if(active_flag) { 
        bool inquireCoeffs;
        ar & inquireCoeffs;
        arout & inquireCoeffs;
        if(inquireCoeffs) {
          Tensor<T> Coefficients(2*k, 2*k, 2*k);
          ar & Coefficients;
          arout & Coefficients;
        }
      }
      FORIJK( 
        //bool have_child;
        ar & have_child;
        arout & have_child;
        if(have_child) {
          ar & active_flag;
          arin & local_remoteRank & n_local & x_local & y_local & z_local;
          if(local_remoteRank == 0) {
            OctTreeT *treemaster = tree->find(n_local, x_local, y_local, z_local);
            loadManager(ar, treemaster, commFunc, active_flag, have_child);
          }
          else {
            if(remoteRank != local_remoteRank) {
              madness::archive::MPIOutputArchive arout2(commFunc, local_remoteRank);
              arout2 & 1;
            }
            shadowManager_load(ar, tree, n_local, 
              x_local, y_local, 
              z_local, local_remoteRank, commFunc, active_flag, have_child);
          }
        }
      );
    }

    template <typename T>
    void Function<T>::saveLoadWorker(OctTreeT *tree, 
		Communicator& commFunc, bool save) 
    {
      int msg;
      madness::archive::MPIInputArchive arrecv(commFunc, 0);
      madness::archive::MPIOutputArchive arsend(commFunc, 0);
      while (1)
      {
        arrecv & msg;
        if (msg == -1) {
          break;
        }
        else if (msg == 1) {
          if(save) {
            Level n_c;
            Translation x_c, y_c, z_c;
            arrecv & n_c & x_c & y_c & z_c;
            OctTreeT *treeclient = tree->find(n_c, x_c, y_c, z_c);
            std::vector<localTreeMember> subtreeList; 
            localTreeList(subtreeList, treeclient);
            int nRemoteBranch = subtreeList.size();
            arsend & nRemoteBranch;
            for (int i = 0; i < nRemoteBranch; i++) {
              arsend & subtreeList[i];
            }
            sendRecvDataWorker_save(treeclient, commFunc);
          }
          else {
            sendRecvDataWorker_load(tree, commFunc);
          }
        }
      }
    }

    template <typename T>
    void Function<T>::localTreeList(std::vector<localTreeMember> &subtreeList, 
		OctTreeT *tree)
    {
      localTreeMember branchList;
      if(tree) {
        branchList.have_child = true;
      }
      else {
        branchList.have_child = false;
      }
      if(branchList.have_child) {
        branchList.x      = tree->x();
        branchList.y      = tree->y();
        branchList.z      = tree->z();
        branchList.n      = tree->n();
        if(tree->rank() == -1) {
          branchList.rank   = comm()->rank();
        }
        else {
          branchList.rank   = tree->rank();
        }
        branchList.remote = isremote(tree);
        branchList.active = isactive(tree);
      }
      subtreeList.push_back(branchList);
      FORIJK( 
        OctTreeT *child = tree->child(i,j,k);
        if(child) {
          localTreeList(subtreeList, child);
        }
        else {
          localTreeMember branchList;
          branchList.have_child = false;
          branchList.x      = i;
          branchList.y      = j;
          branchList.z      = k;
          branchList.rank   = 0;
          branchList.remote   = false;
          branchList.active   = false;
          subtreeList.push_back(branchList);
        }
      );
    }

    template <typename T>
    void Function<T>::sendRecvDataWorker_save(OctTreeT *tree, 
		Communicator& commFunc)
    {
      madness::archive::MPIInputArchive arrecv(commFunc, 0);
      madness::archive::MPIOutputArchive arsend(commFunc, 0);
      if(isactive(tree)) {
        if ( tree->isremote() && !(tree->parent()) ) {
        }
        else {
          TensorT *t = coeff(tree);
          bool coeffsPointer = false;
          if(t) coeffsPointer = true;
          commFunc.Send(&coeffsPointer, 1, 0, 14);
          if(t) {
            arsend & *t;
          }
        }
      }
      FOREACH_LOCAL_CHILD(OctTreeT, tree, 
        sendRecvDataWorker_save(child, commFunc);
      );
    }

    template <typename T>
    void Function<T>::sendRecvDataWorker_load(OctTreeT *treeclient,  
		Communicator& commFunc)
    {
      madness::archive::MPIInputArchive arrecv(commFunc, 0);
      madness::archive::MPIOutputArchive arsend(commFunc, 0);
      Level n_c;
      Translation x_c, y_c, z_c;
      arrecv & n_c & x_c & y_c & z_c;
      OctTreeT *tree = treeclient->find(n_c, x_c, y_c, z_c);
      bool have_child, active_flag;
      arrecv & have_child & active_flag;
      if(active_flag) {
        set_active(tree);
        bool inquireCoeffs;
        arrecv & inquireCoeffs;
        if(inquireCoeffs) {
          Tensor<T> Coefficients(2*k, 2*k, 2*k);
          arrecv & Coefficients;
          set_coeff(tree, Coefficients);
        }
      }
      FORIJK( 
        OctTreeT* child = tree->child(i, j, k);
        arrecv & have_child;
        if(!child) {
          if(have_child) {
	   child = tree->insert_local_child(i,j,k);
          }
        }
        if(child) {
          if(child->rank() == -1) {
            arsend & comm()->rank() & child->n() & child->x() & child->y() & child->z();
          }
          else {
            arsend & child->rank() & child->n() & child->x() & child->y() & child->z();
          }
          if(islocal(child)) {
            sendRecvDataWorker_load(child, commFunc);
          }
        }
      );
    }

    template <typename T>
    void Function<T>::load(const char* f, Communicator& comm) {
      TextFstreamInputArchive iar(comm.rank() ? 0 : f);
      if (comm.rank() == 0) {
        iar & FunctionDefaults::k;
        iar & FunctionDefaults::thresh;
        iar & FunctionDefaults::initial_level;
        iar & FunctionDefaults::max_refine_level;
        iar & FunctionDefaults::truncate_method;
        iar & FunctionDefaults::autorefine;
      }
      comm.Bcast(FunctionDefaults::k, 0);
      comm.Bcast(FunctionDefaults::thresh, 0);
      comm.Bcast(FunctionDefaults::initial_level, 0);
      comm.Bcast(FunctionDefaults::max_refine_level, 0);
      comm.Bcast(FunctionDefaults::truncate_method, 0);
      comm.Bcast(FunctionDefaults::autorefine, 0);
      if (comm.rank() == 0) {
        bool active_flag;
        iar & active_flag;
        if(active_flag) loadManager(iar, tree(), comm, active_flag, true);
        for ( int i = 1; i < comm.size(); i++) {
          archive::MPIOutputArchive arout(comm, i);
          arout & -1;
        }
        iar.close();
      }
      else {
        saveLoadWorker(tree(), comm, false);
      }
    }

    template <typename T>
    template <class Archive>
    void Function<T>::loadManager(const Archive& ar, OctTreeT *tree, Communicator& commFunc, bool active_flag, bool have_child) {
      if (isremote(tree)) {
        madness::archive::MPIOutputArchive arout(commFunc, tree->rank());
        arout & 1;
        if(active_flag){
	 set_active(tree);
	}
	else {
	 set_inactive(tree);
	}
        shadowManager_load(ar, tree, tree->n(), tree->x(), tree->y(), 
            tree->z(), tree->rank(), commFunc, active_flag, have_child);
      }
      else {
        if(active_flag){
	 set_active(tree);
	}
	else {
	 set_inactive(tree);
	}
        Level n_local;
        Translation x_local, y_local, z_local;
        ar & n_local & x_local & y_local & z_local;
        cout << n_local << " " << x_local << " " << y_local << " " << z_local << endl;
	if(active_flag) {
          bool inquireCoeffs;
          ar & inquireCoeffs;
          if(inquireCoeffs) {
            TensorT t;
            ar & t;
            set_coeff(tree, t);
          }
	}
        FORIJK( 
          //bool have_child;
          ar & have_child;
          if(have_child) {
            OctTreeT* child = tree->child(i, j, k);
            if(!child) child = tree->insert_local_child(i,j,k);
            //bool active_flag;
            ar & active_flag;
            loadManager(ar, child, commFunc, active_flag, have_child);
          }
        );
      }
    }

/*
    template <typename T>
    class MulNodeInfo {
    public:
        const Tensor<T>* t;
        Level n;
        Translation lx, ly, lz;
        bool status;

        MulNodeInfo()
            : t(0),n(0),lx(0),ly(0),lz(0),status(false) {}

        void set(const Tensor<T>* t, Level n,
                 Translation lx, Translation ly, Translation lz) {
            this->t = t;
            this->n = n;
            this->lx = lx;
            this->ly = ly;
            this->lz = lz;
            status = true;
        template <typename T>
        class MulNodeInfo {
        public:
            const Tensor<T>* t;
            Level n;
            Translation lx, ly, lz;
            bool status;
     
            MulNodeInfo()
                : t(0),n(0),lx(0),ly(0),lz(0),status(false) {}
     
            void set(const Tensor<T>* t, Level n,
                     Translation lx, Translation ly, Translation lz) {
                this->t = t;
                this->n = n;
                this->lx = lx;
                this->ly = ly;
                this->lz = lz;
                status = true;
            };
     
            operator bool() {return status;};
     
            void send(Process rank);
     
            void recv(Process rank);
>>>>>>> .r80
        };
     
        /// At the quadrature points for box (n,l) evaluate the scaling functions (nn,ll)
        template <typename T>
        Tensor<double> Function<T>::mul_eval_sf(Level n, Translation l, Level nn, Translation ll) {
            int npt = data->cdata->npt;
            Tensor<double> phi(npt,k); // phi(i,j) = at  value of phi[j](x[i])
     
            Tensor<double>& quad_x = data->cdata->quad_x;
     
            double scale = 1.0/two_to_power(n);
            double fac = sqrt(double(two_to_power(nn)));
            if (n < nn) throw "mul_eval_sf: n < nn ?";
            double twonn = two_to_power(n-nn);
            double xlo = l - twonn*ll;
            for (int i=0; i<npt; i++) {
                double p[200];
                double x = scale*(xlo + quad_x(i));
                legendre_scaling_functions(x,k,p);
                for (int j=0; j<k; j++) {
                    phi(i,j) = phi[j]*fac
                }
            }
        }
     
        template <typename T>
        T Function<T>::_mul(Function<T>& afun, Function<T>& bfun,
                            Function<T>& cfun, OctTreeT* tree, 
                            MulNodeInfo<T>& ainfo, MulNodeInfo<T>& binfo) {
            cfun.set_active(tree);
     
            if (tree->isremote()) {
                if (tree->islocalsubtreeparent()) {
                    // If a function is locally labelled active, then the
                    // coefficients are somewhere below this parent node
                    // (perhaps remote).  If a function is inactive, then
                    // the coefficients must be above and we need to
                    // receive them.
                    if (!afun.isactive()) ainfo.recv(tree->rank());
                    if (!bfun.isactive()) binfo.recv(tree->rank());
                }
                else {
                    if (ainfo) ainfo.send(tree->rank());
                    if (binfo) binfo.send(tree->rank());
                    return;
                }
            }
     
            if (at) ainfo.set(at, tree->n(), tree->lx(), tree->ly(), tree->lz());
            if (bt) binfo.set(at, tree->n(), tree->lx(), tree->ly(), tree->lz());
            if (ainfo && binfo) {
                TensorT r(2*k,2*k,2*k);
                // All of the tensor constructors can be optimized away
                // with just a little work.
                FORIJKL(Level n = tree->n()+1;
                        Translation lx = tree->x()*2+i;
                        Translation ly = tree->y()*2+j;
                        Translation lz = tree->z()*2+k;
                        Tensor<double> aphix = mul_eval_sf(n,lx,a.n,a.lx);
                        Tensor<double> aphiy = mul_eval_sf(n,ly,a.n,a.ly);
                        Tensor<double> aphiz = mul_eval_sf(n,lz,a.n,a.lz);
                        Tensor<double> bphix = mul_eval_sf(n,lx,b.n,b.lx);
                        Tensor<double> bphiy = mul_eval_sf(n,ly,b.n,b.ly);
                        Tensor<double> bphiz = mul_eval_sf(n,lz,b.n,b.lz);
                        
                        TensorT aval = transform3d_3c(*ainfo.t,aphix,aphiy,aphiz);
                        TensorT bval = transform3d_3c(*binfo.t,bphix,bphiy,bphiz);
                        r(s[i],s[j],s[k]) = transform3d_inplace(aval.emul(bval),
                                                                data->cdata->quad_phiw,
                                                                data->cdata->work);
                        );
            }
            else {
                FOREACH_CHILD(OctTreeT, tree,
                              if (afun.isactive(child) || bfun.isactive(child)) {
                                  _mul(afun,bfun,cfun,child,ainfo,binfo);
                              });
            }
        }
    */


    // Explicit instantiations for double and complex<double>

    template class Function<double>;
    template class Function< std::complex<double> >;
    template class FunctionData<double>;
    template class FunctionData< std::complex<double> >;
    template class FunctionCommonData<double>;
    template void Function<double>::save_local<TextFstreamOutputArchive>(TextFstreamOutputArchive& ar);
    template void Function<double>::_save_local<TextFstreamOutputArchive>(TextFstreamOutputArchive& ar, const OctTreeT *tree); 
    template void Function<double>::load_local<TextFstreamInputArchive>(const TextFstreamInputArchive& ar); 
    template void Function<double>::_load_local<TextFstreamInputArchive>(const TextFstreamInputArchive& ar, OctTreeT *tree); 
    template void Function<double>::saveMain<TextFstreamOutputArchive>(TextFstreamOutputArchive& ar, Communicator& comm); 
    template void Function<double>::saveManager<TextFstreamOutputArchive>(TextFstreamOutputArchive& ar, const OctTreeT *tree,
              Communicator& commFunc); 
    template void Function<double>::shadowManager_save<TextFstreamOutputArchive>(TextFstreamOutputArchive& ar, const OctTreeT *tree, Level n,
              Translation x, Translation y, Translation z,
              ProcessID remoteRank, Communicator& commFunc); 
    template void Function<double>::shadowManager_load<TextFstreamInputArchive>(const TextFstreamInputArchive& ar, OctTreeT* tree, 
              Level n, Translation x, Translation y, Translation z, 
              ProcessID remoteRank, Communicator& commFunc,
	      bool active_flag, bool have_child); 
    template void Function<double>::loadManager<TextFstreamInputArchive>(const TextFstreamInputArchive& ar, OctTreeT *tree,
              Communicator& commFunc, bool active_flag, bool have_child); 
}

