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

  
#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>

extern "C" double round(double x);


/// \file mra.cc
/// \file Declaration and initialization of static data, some implementation, some instantiation

namespace madness {

    /// For sorting displacements into ascending order by distance
    template <int NDIM>
    static bool cmp_keys(const Key<NDIM>& a, const Key<NDIM>& b) {
        return a.distsq() < b.distsq();
    }


    // Definition and initialization of FunctionDefaults static members
    // It cannot be an instance of FunctionFactory since we want to
    // set the defaults independent of the data type.  


    template <typename T, int NDIM>
    void FunctionCommonData<T,NDIM>::_make_disp() {
        Vector<Translation,NDIM> d;
        int bmax;
        if (NDIM == 1) bmax = 7; // !! Make sure that SimpleCache is consistent!!!!
        else if (NDIM == 2) bmax = 7;
        else if (NDIM == 3) bmax = 3;
        else bmax = 3;

        int num = 1;
        for (int i=0; i<NDIM; i++) num *= (2*bmax + 1);
        disp = std::vector< Key<NDIM> >(num);

        num = 0;
        if (NDIM == 1) {
            for (d[0]=-bmax; d[0]<=bmax; d[0]++)
                disp[num++] = Key<NDIM>(0,d);
        }
        else if (NDIM == 2) {
            for (d[0]=-bmax; d[0]<=bmax; d[0]++)
                for (d[1]=-bmax; d[1]<=bmax; d[1]++)
                    disp[num++] = Key<NDIM>(0,d);
        }
        else if (NDIM == 3) {
            for (d[0]=-bmax; d[0]<=bmax; d[0]++)
                for (d[1]=-bmax; d[1]<=bmax; d[1]++)
                    for (d[2]=-bmax; d[2]<=bmax; d[2]++)
                        disp[num++] = Key<NDIM>(0,d);
        }
        else if (NDIM == 4) {
            for (d[0]=-bmax; d[0]<=bmax; d[0]++)
                for (d[1]=-bmax; d[1]<=bmax; d[1]++)
                    for (d[2]=-bmax; d[2]<=bmax; d[2]++)
                        for (d[3]=-bmax; d[3]<=bmax; d[3]++)
                            disp[num++] = Key<NDIM>(0,d);
        }
        else if (NDIM == 5) {
            for (d[0]=-bmax; d[0]<=bmax; d[0]++)
                for (d[1]=-bmax; d[1]<=bmax; d[1]++)
                    for (d[2]=-bmax; d[2]<=bmax; d[2]++)
                        for (d[3]=-bmax; d[3]<=bmax; d[3]++)
                            for (d[4]=-bmax; d[4]<=bmax; d[4]++)
                                disp[num++] = Key<NDIM>(0,d);
        }
        else if (NDIM == 6) {
            for (d[0]=-bmax; d[0]<=bmax; d[0]++)
                for (d[1]=-bmax; d[1]<=bmax; d[1]++)
                    for (d[2]=-bmax; d[2]<=bmax; d[2]++)
                        for (d[3]=-bmax; d[3]<=bmax; d[3]++)
                            for (d[4]=-bmax; d[4]<=bmax; d[4]++)
                                for (d[5]=-bmax; d[5]<=bmax; d[5]++)
                                    disp[num++] = Key<NDIM>(0,d);
        }
        else {
            MADNESS_EXCEPTION("_make_disp: hard dimension loop",NDIM);
        }

        std::sort(disp.begin(), disp.end(), cmp_keys<NDIM>);
    }


    template <typename T, int NDIM>
    void FunctionCommonData<T,NDIM>::_make_dc_periodic() {
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

    template <typename T, int NDIM>
    void FunctionCommonData<T,NDIM>::_init_twoscale() {
        if (! two_scale_hg(k, &hg)) throw "failed to get twoscale coefficients";
        hgT = transpose(hg);
        hgsonly = copy(hg(Slice(0,k-1),_));
    }

    template <typename T, int NDIM>
    void FunctionCommonData<T,NDIM>::_init_quadrature
    (int k, int npt, Tensor<double>& quad_x, Tensor<double>& quad_w, 
     Tensor<double>& quad_phi, Tensor<double>& quad_phiw, Tensor<double>& quad_phit) {
        quad_x = Tensor<double>(npt);
        quad_w = Tensor<double>(npt);
        quad_phi = Tensor<double>(npt,k);
        quad_phiw = Tensor<double>(npt,k);

        gauss_legendre(npt,0.0,1.0,quad_x.ptr(),quad_w.ptr());
        for (int mu=0; mu<npt; mu++) {
            double phi[200];
            legendre_scaling_functions(quad_x(mu),k,phi);
            for (int j=0; j<k; j++) {
                quad_phi(mu,j) = phi[j];
                quad_phiw(mu,j) = quad_w(mu)*phi[j];
            }
        }
        quad_phit = transpose(quad_phi);
    }


    template <typename T, int NDIM>
    void FunctionImpl<T,NDIM>::verify_tree() const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        world.gop.fence();  // Make sure nothing is going on
        
        // Verify consistency of compression status, existence and size of coefficients,
        // and has_children() flag.
        for(typename dcT::const_iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
            const keyT& key = it->first;
            const nodeT& node = it->second;
            bool bad;
            
            if (is_compressed()) {
                if (node.has_children()) {
                    bad = node.coeff().dim[0] != 2*cdata.k;
                }
                else {
                    bad = node.coeff().size != 0;
                }
            }
            else {
                if (node.has_children()) {
                    bad = node.coeff().size != 0;
                }
                else {
                    bad = node.coeff().dim[0] != cdata.k;
                }
            }
            
            if (bad) {
                print(world.rank(), "FunctionImpl: verify: INCONSISTENT TREE NODE, key =", key, ", node =", node,
                      ", dim[0] =",node.coeff().dim[0],", compressed =",is_compressed());
                std::cout.flush();
                MADNESS_EXCEPTION("FunctionImpl: verify: INCONSISTENT TREE NODE", 0);
            }
        }
        
        // Ensure that parents and children exist appropriately
        for(typename dcT::const_iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
            const keyT& key = it->first;
            const nodeT& node = it->second;
            
            if (key.level() > 0) {
                const keyT parent = key.parent();
                typename dcT::const_iterator pit = coeffs.find(parent).get();
                if (pit == coeffs.end()) {
                    print(world.rank(), "FunctionImpl: verify: MISSING PARENT",key,parent);
                    std::cout.flush();
                    MADNESS_EXCEPTION("FunctionImpl: verify: MISSING PARENT", 0);
                }
                const nodeT& pnode = pit->second;
                if (!pnode.has_children()) {
                    print(world.rank(), "FunctionImpl: verify: PARENT THINKS IT HAS NO CHILDREN",key,parent);
                    std::cout.flush();
                    MADNESS_EXCEPTION("FunctionImpl: verify: PARENT THINKS IT HAS NO CHILDREN", 0);
                }
            }
            
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                typename dcT::const_iterator cit = coeffs.find(kit.key()).get();
                if (cit == coeffs.end()) {
                    if (node.has_children()) {
                        print(world.rank(), "FunctionImpl: verify: MISSING CHILD",key,kit.key());
                        std::cout.flush();
                        MADNESS_EXCEPTION("FunctionImpl: verify: MISSING CHILD", 0);
                    }
                }
                else {
                    if (! node.has_children()) {
                        print(world.rank(), "FunctionImpl: verify: UNEXPECTED CHILD",key,kit.key());
                        std::cout.flush();
                        MADNESS_EXCEPTION("FunctionImpl: verify: UNEXPECTED CHILD", 0);
                    }
                }
            }
        }
        
        world.gop.fence();
    }

    template <typename T, int NDIM>
    T FunctionImpl<T,NDIM>::eval_cube(Level n, coordT x, const tensorT c) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        const int k = cdata.k;
        double px[NDIM][k];
        T sum = T(0.0);
        
        for (int i=0; i<NDIM; i++) legendre_scaling_functions(x[i],k,px[i]);
        
        if (NDIM == 1) {
            for (int p=0; p<k; p++) 
                sum += c(p)*px[0][p];
        }
        else if (NDIM == 2) {
            for (int p=0; p<k; p++) 
                for (int q=0; q<k; q++) 
                    sum += c(p,q)*px[0][p]*px[1][q];
        }
        else if (NDIM == 3) {
            for (int p=0; p<k; p++) 
                for (int q=0; q<k; q++) 
                    for (int r=0; r<k; r++) 
                        sum += c(p,q,r)*px[0][p]*px[1][q]*px[2][r];
        }
        else if (NDIM == 4) {
            for (int p=0; p<k; p++) 
                for (int q=0; q<k; q++) 
                    for (int r=0; r<k; r++) 
                        for (int s=0; s<k; s++) 
                            sum += c(p,q,r,s)*px[0][p]*px[1][q]*px[2][r]*px[3][s];
        }
        else if (NDIM == 5) {
            for (int p=0; p<k; p++) 
                for (int q=0; q<k; q++) 
                    for (int r=0; r<k; r++) 
                        for (int s=0; s<k; s++) 
                            for (int t=0; t<k; t++) 
                                sum += c(p,q,r,s,t)*px[0][p]*px[1][q]*px[2][r]*px[3][s]*px[4][t];
        }
        else if (NDIM == 6) {
            for (int p=0; p<k; p++) 
                for (int q=0; q<k; q++) 
                    for (int r=0; r<k; r++) 
                        for (int s=0; s<k; s++) 
                            for (int t=0; t<k; t++) 
                                for (int u=0; u<k; u++) 
                                    sum += c(p,q,r,s,t,u)*px[0][p]*px[1][q]*px[2][r]*px[3][s]*px[4][t]*px[5][u];
        }
        else {
            MADNESS_EXCEPTION("FunctionImpl:eval_cube:NDIM?",NDIM);
        }
        return sum*pow(2.0,0.5*NDIM*n)/sqrt(FunctionDefaults<NDIM>::get_cell_volume());
    }

    template <typename T, int NDIM>
    Void FunctionImpl<T,NDIM>::reconstruct_op(const keyT& key, const tensorT& s) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        // Note that after application of an integral operator not all
        // siblings may be present so it is necessary to check existence
        // and if absent insert an empty leaf node.
        //
        // If summing the result of an integral operator (i.e., from
        // non-standard form) there will be significant scaling function
        // coefficients at all levels and possibly difference coefficients
        // in leaves, hence the tree may refine as a result.
        typename dcT::iterator it = coeffs.find(key).get();
        if (it == coeffs.end()) {
            coeffs.insert(key,nodeT(tensorT(),false));
            it = coeffs.find(key).get();
        }
        nodeT& node = it->second;
        
        // The integral operator will correctly connect interior nodes
        // to children but may leave interior nodes without coefficients
        // ... but they still need to sum down so just give them zeros
        if (node.has_children() && !node.has_coeff()) {
            node.set_coeff(tensorT(cdata.v2k));
        }
        
        if (node.has_coeff()) {
            tensorT d = node.coeff();
            if (key.level() > 0) d(cdata.s0) += s; // -- note accumulate for NS summation
            d = unfilter(d);
            node.clear_coeff();
            node.set_has_children(true); 
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                const keyT& child = kit.key();
                tensorT ss = copy(d(child_patch(child)));
                task(coeffs.owner(child), &implT::reconstruct_op, child, ss);
            }
        }
        else {
            if (key.level()) node.set_coeff(copy(s));
            else node.set_coeff(s);
        }
        return None;
    }


    template <typename T, int NDIM>
    void FunctionImpl<T,NDIM>::fcube(const keyT& key, const FunctionFunctorInterface<T,NDIM>& f, const Tensor<double>& qx, tensorT& fval) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        const Vector<Translation,NDIM>& l = key.translation();
        const Level n = key.level();
        const double h = std::pow(0.5,double(n)); 
        coordT c; // will hold the point in user coordinates
        const int npt = qx.dim[0];
        
        const Tensor<double>& cell_width = FunctionDefaults<NDIM>::get_cell_width();
        const Tensor<double>& cell = FunctionDefaults<NDIM>::get_cell();

        if (NDIM == 1) {
            for (int i=0; i<npt; i++) {
                c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                fval(i) = f(c);
            }
        }
        else if (NDIM == 2) {
            for (int i=0; i<npt; i++) {
                c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                for (int j=0; j<npt; j++) {
                    c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                    fval(i,j) = f(c);
                }
            }
        }
        else if (NDIM == 3) {
            for (int i=0; i<npt; i++) {
                c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                for (int j=0; j<npt; j++) {
                    c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                    for (int k=0; k<npt; k++) {
                        c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                        fval(i,j,k) = f(c);
                    }
                }
            }
        }
        else if (NDIM == 4) {
            for (int i=0; i<npt; i++) {
                c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                for (int j=0; j<npt; j++) {
                    c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                    for (int k=0; k<npt; k++) {
                        c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                        for (int m=0; m<npt; m++) {
                            c[3] = cell(3,0) + h*cell_width[3]*(l[3] + qx(m)); // xx
                            fval(i,j,k,m) = f(c);
                        }
                    }
                }
            }
        }
        else if (NDIM == 5) {
            for (int i=0; i<npt; i++) {
                c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                for (int j=0; j<npt; j++) {
                    c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                    for (int k=0; k<npt; k++) {
                        c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                        for (int m=0; m<npt; m++) {
                            c[3] = cell(3,0) + h*cell_width[3]*(l[3] + qx(m)); // xx
                            for (int n=0; n<npt; n++) {
                                c[4] = cell(4,0) + h*cell_width[4]*(l[4] + qx(n)); // yy
                                fval(i,j,k,m,n) = f(c);
                            }
                        }
                    }
                }
            }
        }
        else if (NDIM == 6) {
            for (int i=0; i<npt; i++) {
                c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                for (int j=0; j<npt; j++) {
                    c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                    for (int k=0; k<npt; k++) {
                        c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                        for (int m=0; m<npt; m++) {
                            c[3] = cell(3,0) + h*cell_width[3]*(l[3] + qx(m)); // xx
                            for (int n=0; n<npt; n++) {
                                c[4] = cell(4,0) + h*cell_width[4]*(l[4] + qx(n)); // yy
                                for (int p=0; p<npt; p++) {
                                    c[5] = cell(5,0) + h*cell_width[5]*(l[5] + qx(p)); // zz
                                    fval(i,j,k,m,n,p) = f(c);
                                }
                            }
                        }
                    }
                }
            }
        }
        else {
            MADNESS_EXCEPTION("FunctionImpl: fcube: confused about NDIM?",NDIM);
        }
    }

    template <typename T, int NDIM>
    Void FunctionImpl<T,NDIM>::project_refine_op(const keyT& key, bool do_refine) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        if (do_refine) {
            // Make in r child scaling function coeffs at level n+1
            tensorT r(cdata.v2k);
            for (KeyChildIterator<NDIM> it(key); it; ++it) {
                const keyT& child = it.key();
                r(child_patch(child)) = project(child);
            }
            // Filter then test difference coeffs at level n
            tensorT d = filter(r);
            tensorT s0;
            if (truncate_on_project) s0 = copy(d(cdata.s0));
            d(cdata.s0) = T(0);
            if ( (d.normf() < truncate_tol(thresh,key.level())) && 
                 (key.level() < max_refine_level) ) {
                if (truncate_on_project) {
                    coeffs.insert(key,nodeT(s0,false));
                }
                else {
                    coeffs.insert(key,nodeT(tensorT(),true)); // Insert empty node for parent
                    for (KeyChildIterator<NDIM> it(key); it; ++it) {
                        const keyT& child = it.key();
                        coeffs.insert(child,nodeT(copy(r(child_patch(child))),false));
                    }
                }
            }
            else {
                coeffs.insert(key,nodeT(tensorT(),true)); // Insert empty node for parent
                if (key.level()==max_refine_level) print("MAX REFINE LEVEL",key);
                for (KeyChildIterator<NDIM> it(key); it; ++it) {
                    const keyT& child = it.key();
                    ProcessID p;
                    if (FunctionDefaults<NDIM>::get_project_randomize()) {
                        p = world.random_proc();
                    }
                    else {
                        p = coeffs.owner(child);
                    }
                    task(p, &implT::project_refine_op, child, do_refine); // ugh
                }
            }
        }
        else {
            coeffs.insert(key,nodeT(project(key),false));
        }
        return None;
    }

    template <typename T, int NDIM>
    void FunctionImpl<T,NDIM>::add_scalar_inplace(T t, bool fence) {
        std::vector<long> v0(NDIM,0L);
        if (is_compressed()) {
            if (world.rank() == coeffs.owner(cdata.key0)) {
                typename dcT::iterator it = coeffs.find(cdata.key0).get();
                MADNESS_ASSERT(it != coeffs.end());
                nodeT& node = it->second;
                MADNESS_ASSERT(node.has_coeff());
                node.coeff()(v0) += t*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
            }
        }
        else {
            for(typename dcT::iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
                Level n = it->first.level();
                nodeT& node = it->second;
                if (node.has_coeff()) {
                    node.coeff()(v0) += t*sqrt(FunctionDefaults<NDIM>::get_cell_volume()*pow(0.5,double(NDIM*n)));
                }
            }
        }
        if (fence) world.gop.fence();
    }

    template <typename T, int NDIM>
    void FunctionImpl<T,NDIM>::insert_zero_down_to_initial_level(const keyT& key) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        if (compressed) initial_level = std::max(initial_level,1); // Otherwise zero function is confused
        if (coeffs.is_local(key)) {
            if (compressed) {
                if (key.level() == initial_level) {
                    coeffs.insert(key, nodeT(tensorT(), false));
                }
                else {
                    coeffs.insert(key, nodeT(tensorT(cdata.v2k), true));
                }
            }
            else {
                if (key.level()<initial_level) {
                    coeffs.insert(key, nodeT(tensorT(), true));
                }
                else {
                    coeffs.insert(key, nodeT(tensorT(cdata.vk), false));
                }
            }
        }
        if (key.level() < initial_level) {
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                insert_zero_down_to_initial_level(kit.key());
            }
        }
        
    }
    
    
    template <typename T, int NDIM>
    Future<bool> FunctionImpl<T,NDIM>::truncate_spawn(const keyT& key, double tol) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        const nodeT& node = coeffs.find(key).get()->second;   // Ugh!
        if (node.has_children()) {
            std::vector< Future<bool> > v = future_vector_factory<bool>(1<<NDIM);
            int i=0;
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
                v[i] = send(coeffs.owner(kit.key()), &implT::truncate_spawn, kit.key(), tol);
            }
            return task(world.rank(),&implT::truncate_op, key, tol, v);
        }
        else {
            MADNESS_ASSERT(!node.has_coeff());  // In compressed form leaves should not have coeffs
            return Future<bool>(false); 
        }
    }
    

    template <typename T, int NDIM>
    bool FunctionImpl<T,NDIM>::truncate_op(const keyT& key, double tol, const std::vector< Future<bool> >& v) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        // If any child has coefficients, a parent cannot truncate
        for (int i=0; i<(1<<NDIM); i++) if (v[i].get()) return true;
        nodeT& node = coeffs.find(key).get()->second;
        if (key.level() > 1) { // >1 rather >0 otherwise reconstruct might get confused
            double dnorm = node.coeff().normf();
            if (dnorm < truncate_tol(tol,key)) {
                node.clear_coeff();
                node.set_has_children(false);
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    coeffs.erase(kit.key());
                }
            }
        }
        return node.has_coeff();
    }


    
    template <typename T, int NDIM>
    void FunctionImpl<T,NDIM>::print_tree(Level maxlevel) const {
        if (world.rank() == 0) do_print_tree(cdata.key0, maxlevel);
        world.gop.fence();
        if (world.rank() == 0) std::cout.flush();
        world.gop.fence();
    }
    
    
    template <typename T, int NDIM>
    void FunctionImpl<T,NDIM>::do_print_tree(const keyT& key, Level maxlevel) const {
        typename dcT::const_iterator it = coeffs.find(key).get();
        if (it == coeffs.end()) {
            MADNESS_EXCEPTION("FunctionImpl: do_print_tree: null node pointer",0);
        }
        const nodeT& node = it->second;
        for (int i=0; i<key.level(); i++) std::cout << "  ";
        std::cout << key << "  " << node << " --> " << coeffs.owner(key) << "\n";
        if (key.level() < maxlevel  &&  node.has_children()) {
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                do_print_tree(kit.key(),maxlevel);
            }
        }
    }
    
    template <typename T, int NDIM>
    Tensor<T> FunctionImpl<T,NDIM>::project(const keyT& key) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        tensorT fval(cdata.vq,false); // this will be the returned result
        tensorT& work = cdata.work1; // initially evaluate the function in here
        
        if (functor) {
            fcube(key,*functor,cdata.quad_x,work);
        }
        else {
            MADNESS_EXCEPTION("FunctionImpl: project: confusion about function?",0);
        }
        
        work.scale(sqrt(FunctionDefaults<NDIM>::get_cell_volume()*pow(0.5,double(NDIM*key.level()))));
        //return transform(work,cdata.quad_phiw);
        return fast_transform(work,cdata.quad_phiw,fval,cdata.workq);
    }

    template <typename T, int NDIM>
    Future<double> FunctionImpl<T,NDIM>::get_norm_tree_recursive(const keyT& key) const {
        if (coeffs.probe(key)) {
            return Future<double>(coeffs[key].get_norm_tree());
        }
        MADNESS_ASSERT(key.level());
        keyT parent = key.parent();
        return send(coeffs.owner(parent), &implT::get_norm_tree_recursive, parent);
    }


    template <typename T, int NDIM>
    Void FunctionImpl<T,NDIM>::sock_it_to_me(const keyT& key, 
                                             const RemoteReference< FutureImpl< std::pair<keyT,tensorT> > >& ref) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        if (coeffs.probe(key)) {
            const nodeT& node = coeffs.find(key).get()->second;
            Future< std::pair<keyT,tensorT> > result(ref);
            if (node.has_coeff()) {
                //madness::print("sock found it with coeff",key);
                result.set(std::pair<keyT,tensorT>(key,node.coeff()));
            }
            else {
                //madness::print("sock found it without coeff",key);
                result.set(std::pair<keyT,tensorT>(key,tensorT()));
            }
        }
        else {
            keyT parent = key.parent();
            //madness::print("sock forwarding to parent",key,parent);
            send(coeffs.owner(parent), &FunctionImpl<T,NDIM>::sock_it_to_me, parent, ref);
        }
        return None;
    }

    template <typename T, int NDIM>
    Void FunctionImpl<T,NDIM>::eval(const Vector<double,NDIM>& xin, 
                                    const keyT& keyin, 
                                    const typename Future<T>::remote_refT& ref) {
        
        PROFILE_MEMBER_FUNC(FunctionImpl);
        // This is ugly.  We must figure out a clean way to use 
        // owner computes rule from the container.
        Vector<double,NDIM> x = xin;
        keyT key = keyin;
        Vector<Translation,NDIM> l = key.translation();
        ProcessID me = world.rank();
        while (1) {
            ProcessID owner = coeffs.owner(key);
            if (owner != me) {
                send(owner, &implT::eval, x, key, ref);
                return None;
            }
            else {
                typename dcT::futureT fut = coeffs.find(key);
                typename dcT::iterator it = fut.get();
                nodeT& node = it->second;
                if (node.has_coeff()) {
                    Future<T>(ref).set(eval_cube(key.level(), x, node.coeff()));
                    return None;
                }
                else {
                    for (int i=0; i<NDIM; i++) {
                        double xi = x[i]*2.0;
                        int li = int(xi);
                        if (li == 2) li = 1;
                        x[i] = xi - li;
                        l[i] = 2*l[i] + li;
                    }
                    key = keyT(key.level()+1,l);
                }
            }
        }
        //MADNESS_EXCEPTION("should not be here",0);
    }
    
    template <typename T, int NDIM>
    void FunctionImpl<T,NDIM>::tnorm(const tensorT& t, double* lo, double* hi) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        // Chosen approach looks stupid but it is more accurate
        // than the simple approach of summing everything and
        // subtracting off the low-order stuff to get the high
        // order (assuming the high-order stuff is small relative
        // to the low-order)
        cdata.work1(___) = t(___);
        tensorT tlo = cdata.work1(cdata.sh);
        *lo = tlo.normf();
        tlo.fill(0.0);
        *hi = cdata.work1.normf();
    }
    
    namespace detail {
        template <typename A, typename B>
        struct noop {
            void operator()(const A& a, const B& b) const {};

            template <typename Archive> void serialize(Archive& ar) {}
        };

        template <typename T, int NDIM>
        struct scaleinplace {
            T q;
            scaleinplace() {}
            scaleinplace(T q) : q(q) {}
            void operator()(const Key<NDIM>& key, Tensor<T>& t) const {t.scale(q);}
            template <typename Archive> void serialize(Archive& ar) {ar & q;}
        };            

        template <typename T, int NDIM>
        struct squareinplace {
            void operator()(const Key<NDIM>& key, Tensor<T>& t) const {t.emul(t);}
            template <typename Archive> void serialize(Archive& ar) {}
        };
    }

    template <typename T, int NDIM>
    void FunctionImpl<T,NDIM>::refine(bool fence) 
    {
        //if (autorefine) {
        unary_op_coeff_inplace(&implT::autorefine_square_test, detail::noop<keyT,tensorT>(), fence);
        //}
        //else {
        //    unary_op_coeff_inplace(&implT::noautorefine, detail::noop<keyT,tensorT>(), fence);
        //}
    }

    template <typename T, int NDIM>
    void FunctionImpl<T,NDIM>::scale_inplace(const T q, bool fence) 
    {
        unary_op_coeff_inplace(&implT::noautorefine, detail::scaleinplace<T,NDIM>(q), fence);
    }
    
    template <typename T, int NDIM>
    void FunctionImpl<T,NDIM>::square_inplace(bool fence) {
        unary_op_value_inplace(&implT::autorefine_square_test, detail::squareinplace<T,NDIM>(), fence);
    }

    template <typename T, int NDIM>
    void FunctionImpl<T,NDIM>::phi_for_mul(Level np, Translation lp, Level nc, Translation lc, Tensor<double>& phi) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        double p[200];
        double scale = pow(2.0,double(np-nc));
        for (int mu=0; mu<cdata.npt; mu++) {
            double xmu = scale*(cdata.quad_x(mu)+lc) - lp;
            MADNESS_ASSERT(xmu>-1e-15 && xmu<(1+1e-15));
            legendre_scaling_functions(xmu,cdata.k,p);
            for (int i=0; i<k; i++) phi(i,mu) = p[i];
        }
        phi.scale(pow(2.0,0.5*np));
    }

    template <typename T, int NDIM>
    const Tensor<T> FunctionImpl<T,NDIM>::parent_to_child(const tensorT& s, const keyT& parent, const keyT& child) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        // An invalid parent/child means that they are out of the box
        // and it is the responsibility of the caller to worry about that
        // ... most likely the coefficiencts (s) are zero to reflect
        // zero B.C. so returning s makes handling this easy.
        if (parent == child || parent.is_invalid() || child.is_invalid()) return s;
        
        tensorT result = fcube_for_mul<T>(child, parent, s);
        result.scale(sqrt(FunctionDefaults<NDIM>::get_cell_volume()*pow(0.5,double(NDIM*child.level()))));
        result = transform(result,cdata.quad_phiw);
        
        return result;
    }


    template <typename T, int NDIM>
    T FunctionImpl<T,NDIM>::trace_local() const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        std::vector<long> v0(NDIM,0);
        T sum = 0.0;
        if (compressed) {
            if (world.rank() == coeffs.owner(cdata.key0)) {
                typename dcT::const_iterator it = coeffs.find(cdata.key0).get();
                if (it != coeffs.end()) {
                    const nodeT& node = it->second;
                    if (node.has_coeff()) sum = node.coeff()(v0);
                }
            }
        }
        else {
            for(typename dcT::const_iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
                const keyT& key = it->first;
                const nodeT& node = it->second;
                if (node.has_coeff()) sum += node.coeff()(v0)*pow(0.5,NDIM*key.level()*0.5);
            }
        }
        return sum*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
    }


    template <typename T, int NDIM>
    Void FunctionImpl<T,NDIM>::recur_down_with_fill(const keyT& target, const keyT& key) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        typename dcT::iterator it = coeffs.find(key).get();

        // a) Node does not exist ... forward request up with the correct target
        // b) Node exists and has coefficients ... start downward recursion
        // c) Node exists does not have children and does not have coefficients ... forward up
        // d) Node exists has children and does not have coefficients ... foward down

        if (it == coeffs.end() || (!it->second.has_coeff() && !it->second.has_children())) { // a and c
            keyT parent = key.parent();
            madness::print("rdwf: forwarding up ac", target, parent);
            send(coeffs.owner(parent), &implT::recur_down_with_fill, target, parent);
        }
        else if (it->second.has_coeff()) { // b
            // Have coefficients ... recur down to make all children telling
            // only the appropriate branch to also recur further
            tensorT s = tensorT(cdata.v2k);
            s(cdata.s0) = it->second.coeff();
            s = unfilter(s);
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                const keyT& child = kit.key();
                tensorT ss = copy(s(child_patch(child)));
                madness::print("rdwf: inserting", child);
                coeffs.insert(child,nodeT(ss,false));
                if (child.level() != target.level() && target.is_child_of(child)) {
                    send(coeffs.owner(child), &implT::recur_down_with_fill, target, child);
                }                    
            }
            it->second.set_has_children(true);
            it->second.clear_coeff();
        }
        else if (key.level() != (target.level()-1)) { // d
            //it->second.set_has_children(true);
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                const keyT& child = kit.key();
                if (target.is_child_of(child)) {
                    send(coeffs.owner(child), &implT::recur_down_with_fill, target, child);
                    break;
                }                    
            }
        }
    
        return None;
    }

    template <typename T, int NDIM>
    Void FunctionImpl<T,NDIM>::ensure_exists(const keyT& key) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        if (!coeffs.probe(key)) {
            keyT parent = key.parent();
            // If the node does not exist this implies that it will
            // be a leaf ... make it here so that we only send one
            // request up the tree to make it.
            coeffs.insert(key,nodeT(tensorT(),false));
            //madness::print("ensure_exists: sending recur up from", key, "to", parent);
            send(coeffs.owner(parent), &implT::recur_down_with_fill, key, parent);            
        }
        return None;
    }


    template <typename T, int NDIM>
    void FunctionImpl<T,NDIM>::widen(bool fence, int ndiff) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        double tol = std::min(1e3*thresh, sqrt(thresh));
        for(typename dcT::iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
            const keyT& key = it->first;
            const nodeT& node = it->second;
            if (node.is_leaf() && node.coeff().normf()>tol) {
                for (int axis=0; axis<NDIM; axis++) {
                    for (int step=-1; step<=1; step+=2) {
                        keyT neigh = neighbor(key, axis, step);
                        if (neigh.is_valid()) {
                            if (ndiff > 0) neigh = neigh.parent(ndiff);
                            send(coeffs.owner(neigh), &implT::ensure_exists, neigh);
                        }
                    }
                }
            }
        }
        if (fence) world.gop.fence();
    }


    template <typename T, int NDIM>
    void FunctionImpl<T,NDIM>::diff(const implT& f, int axis, bool fence) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        typedef std::pair<keyT,tensorT> argT;
        for(typename dcT::const_iterator it=f.coeffs.begin(); it!=f.coeffs.end(); ++it) {
            const keyT& key = it->first;
            const nodeT& node = it->second;
            if (node.has_coeff()) {
                Future<argT> left   = f.find_neighbor(key,axis,-1);
                argT center(key,node.coeff());
                Future<argT> right  = f.find_neighbor(key,axis, 1);
                task(world.rank(), &implT::do_diff1, &f, axis, key, left, center, right, task_attr_generator());
            }
            else {
                // Internal empty node can be safely inserted
                coeffs.insert(key,nodeT(tensorT(),true));
            }
        }
        if (fence) world.gop.fence();
    }
    
    static bool enforce_bc(int bc_left, int bc_right, Level n, Translation& l) {
        Translation two2n = 1ul << n;
        if (l < 0) {
            if (bc_left == 0) {
                return false; // Zero BC
            }
            else if (bc_left == 1) {
                l += two2n; // Periodic BC
            }
            else {
                MADNESS_EXCEPTION("enforce_bc: confused left BC?",bc_left);
            }
        }
        else if (l >= two2n) {
            if (bc_right == 0) {
                return false; // Zero BC
            }
            else if (bc_right == 1) {
                l -= two2n; // Periodic BC
            }
            else {
                MADNESS_EXCEPTION("enforce_bc: confused BC right?",bc_left);
            }
        }
        return true;
    }


    template <typename T, int NDIM>
    Key<NDIM> FunctionImpl<T,NDIM>::neighbor(const keyT& key, int axis, int step) const {
        Vector<Translation,NDIM> l = key.translation();
        
        l[axis] += step;

        if (!enforce_bc(bc(axis,0), bc(axis,1), key.level(), l[axis])) {
            return keyT::invalid();
        }
        else {
            return keyT(key.level(),l);
        }
    }
    
    template <typename T, int NDIM>
    Key<NDIM> FunctionImpl<T,NDIM>::neighbor(const keyT& key, const Key<NDIM>& disp) const {
        Vector<Translation,NDIM> l = key.translation();

        for (int axis=0; axis<NDIM; axis++) {
            l[axis] += disp.translation()[axis];

            if (!enforce_bc(bc(axis,0), bc(axis,1), key.level(), l[axis])) {
                return keyT::invalid();
            }
        }
        return keyT(key.level(),l);
    }
    
    template <typename T, int NDIM>
    Future< std::pair< Key<NDIM>,Tensor<T> > > 
    FunctionImpl<T,NDIM>::find_neighbor(const Key<NDIM>& key, int axis, int step) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        typedef std::pair< Key<NDIM>,Tensor<T> > argT;
        keyT neigh = neighbor(key, axis, step);
        if (neigh.is_invalid()) {
            return Future<argT>(argT(neigh,tensorT(cdata.vk))); // Zero bc
        }
        else {
            Future<argT> result;
            send(coeffs.owner(neigh), &implT::sock_it_to_me, neigh, result.remote_ref(world));
            return result;
        }
    }
    
    template <typename T, int NDIM>
    Void FunctionImpl<T,NDIM>::forward_do_diff1(const implT* f, int axis, const keyT& key, 
                                                const std::pair<keyT,tensorT>& left,
                                                const std::pair<keyT,tensorT>& center,
                                                const std::pair<keyT,tensorT>& right) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        ProcessID owner = coeffs.owner(key);
        if (owner == world.rank()) {
            if (left.second.size == 0) {
                task(owner, &implT::do_diff1, f, axis, key, 
                     f->find_neighbor(key,axis,-1), center, right, 
                     task_attr_generator());
            }
            else if (right.second.size == 0) {
                task(owner, &implT::do_diff1, f, axis, key, 
                     left, center, f->find_neighbor(key,axis,1), 
                     task_attr_generator());
            }
            else {
                task(owner, &implT::do_diff2, f, axis, key, 
                     left, center, right);
            }
        }
        else {
            send(owner, &implT::forward_do_diff1, f, axis, key, 
                 left, center, right);
        }
        return None;
    }
    
    template <typename T, int NDIM>
    Void FunctionImpl<T,NDIM>::do_diff1(const implT* f, int axis, const keyT& key, 
                                        const std::pair<keyT,tensorT>& left,
                                        const std::pair<keyT,tensorT>& center,
                                        const std::pair<keyT,tensorT>& right) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        typedef std::pair<keyT,tensorT> argT;
        
        MADNESS_ASSERT(axis>=0 && axis<NDIM);
        
        if (left.second.size==0 || right.second.size==0) {
            // One of the neighbors is below us in the tree ... recur down
            coeffs.insert(key,nodeT(tensorT(),true));
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                const keyT& child = kit.key();
                if ((child.translation()[axis]&1) == 0) { 
                    // leftmost child automatically has right sibling
                    forward_do_diff1(f, axis, child, left, center, center);
                }
                else { 
                    // rightmost child automatically has left sibling
                    forward_do_diff1(f, axis, child, center, center, right);
                }
            }
        }
        else {
            forward_do_diff1(f, axis, key, left, center, right);
        }
        return None;
    }
    
    template <typename T, int NDIM>
    Void FunctionImpl<T,NDIM>::do_diff2(const implT* f, int axis, const keyT& key, 
                                        const std::pair<keyT,tensorT>& left,
                                        const std::pair<keyT,tensorT>& center,
                                        const std::pair<keyT,tensorT>& right) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        typedef std::pair<keyT,tensorT> argT;
        
        tensorT d = madness::inner(cdata.rp, 
                                   parent_to_child(left.second, left.first, neighbor(key,axis,-1)).swapdim(axis,0), 
                                   1, 0);
        inner_result(cdata.r0, 
                     parent_to_child(center.second, center.first, key).swapdim(axis,0),
                     1, 0, d);
        inner_result(cdata.rm, 
                     parent_to_child(right.second, right.first, neighbor(key,axis,1)).swapdim(axis,0), 
                     1, 0, d);
        if (axis) d = copy(d.swapdim(axis,0)); // make it contiguous
        d.scale(FunctionDefaults<NDIM>::get_rcell_width()[axis]*pow(2.0,(double) key.level()));
        coeffs.insert(key,nodeT(d,false));
        return None;
    }
    
    template <typename T, int NDIM>
    void FunctionImpl<T,NDIM>::mapdim(const implT& f, const std::vector<long>& map, bool fence) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        for(typename dcT::const_iterator it=f.coeffs.begin(); it!=f.coeffs.end(); ++it) {
            const keyT& key = it->first;
            const nodeT& node = it->second;
            
            Vector<Translation,NDIM> l;
            for (int i=0; i<NDIM; i++) l[map[i]] = key.translation()[i];
            tensorT c = node.coeff();
            if (c.size) c = copy(c.mapdim(map));
            
            coeffs.insert(keyT(key.level(),l), nodeT(c,node.has_children()));
        }
        if (fence) world.gop.fence();
    }
    
    
    template <typename T, int NDIM>
    Future< Tensor<T> > FunctionImpl<T,NDIM>::compress_spawn(const Key<NDIM>& key, bool nonstandard, bool keepleaves) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        MADNESS_ASSERT(coeffs.probe(key));
        nodeT& node = coeffs.find(key).get()->second;
        if (node.has_children()) {
            std::vector< Future<tensorT> > v = future_vector_factory<tensorT>(1<<NDIM);
            int i=0;
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
                v[i] = send(coeffs.owner(kit.key()), &implT::compress_spawn, kit.key(), nonstandard, keepleaves);
            }
            return task(world.rank(),&implT::compress_op, key, v, nonstandard);
        }
        else {
            Future<tensorT> result(node.coeff());
            if (!keepleaves) node.clear_coeff();
            return result;
        }
    }

    template <typename T, int NDIM>
    Tensor<T> FunctionImpl<T,NDIM>::eval_plot_cube(const coordT& plotlo,
                                                   const coordT& plothi,
                                                   const std::vector<long>& npt) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        Tensor<T> r(NDIM, &npt[0]);
        //         r(___) = 99.0;
        MADNESS_ASSERT(!compressed);

        coordT h; // Increment between points in each dimension
        for (int i=0; i<NDIM; i++) {
            if (npt[i] > 1) {
                h[i] = (plothi[i]-plotlo[i])/(npt[i]-1);
            }
            else {
                MADNESS_ASSERT(plotlo[i] == plothi[i]);
                h[i] = 0.0;
            }           
        }
        //print("plot info", plotlo, plothi, npt, h);

        // Loop thru local boxes
        for(typename dcT::const_iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
            const keyT& key = it->first;
            const nodeT& node = it->second;
            if (node.has_coeff()) {
                //print("Looking at", key);
                // Determine the points, if any, of the plot grid that
                // are contained within this box
                coordT boxlo, boxhi;
                Vector<int,NDIM> boxnpt;
                double fac = pow(0.5,double(key.level()));
                int npttotal = 1;
                for (int d=0; d<NDIM; d++) {
                    // Coords of box
                    boxlo[d] = fac*key.translation()[d];
                    boxhi[d] = boxlo[d]+fac;
                    
                    if (boxlo[d] > plothi[d] || boxhi[d] < plotlo[d]) {
                        // Discard boxes out of the plot range
                        npttotal = boxnpt[d] = 0;
                        //print("OO range?");
                        break;
                    }
                    else if (npt[d] == 1) {
                        // This dimension is only a single point
                        boxlo[d] = boxhi[d] = plotlo[d];
                        boxnpt[d] = 1;
                    }
                    else {
                        // Restrict to plot range
//                         boxlo[d] = std::max(boxlo[d],plotlo[d]);
//                         boxhi[d] = std::min(boxhi[d],plothi[d]);
                        
                        // Round lo up to next plot point; round hi down
                        double xlo = long((boxlo[d]-plotlo[d])/h[d])*h[d] + plotlo[d];
                        if (xlo < boxlo[d]) xlo += h[d];
                        boxlo[d] =  xlo;
                        double xhi = long((boxhi[d]-plotlo[d])/h[d])*h[d] + plotlo[d];
                        if (xhi > boxhi[d]) xhi -= h[d];
                        // MADNESS_ASSERT(xhi >= xlo);  // nope
                        boxhi[d] = xhi;
                        boxnpt[d] = long(round((boxhi[d] - boxlo[d])/h[d])) + 1;
                    }
                    npttotal *= boxnpt[d];
                }
                //print("    box", boxlo, boxhi, boxnpt, npttotal);
                if (npttotal > 0) {
                    const tensorT& coeff = node.coeff();
                    const Level n = key.level();
                    const Vector<Translation,NDIM>& l = key.translation();
                    const double twon = pow(2.0,double(n));
                    long ind[NDIM];
                    coordT x; 
                    for (IndexIterator it(boxnpt); it; ++it) {
                        for (int d=0; d<NDIM; d++) {
                            double xd = boxlo[d] + it[d]*h[d]; // Sim. coords of point
                            x[d] = twon*xd - l[d]; // Offset within box
                            MADNESS_ASSERT(x[d]>=0.0 && x[d] <=1.0);  // sanity
                            if (npt[d] > 1) {
                                ind[d] = long(round((xd-plotlo[d])/h[d])); // Index of plot point
                            }
                            else {
                                ind[d] = 0;
                            }
                            MADNESS_ASSERT(ind[d]>=0 && ind[d]<npt[d]); // sanity
                        }
                        r(ind) = eval_cube(n, x, coeff);
                        //print("computing", n, x, ind, r(ind));
                    }                        
                }
            }
        }

        //        ITERATOR(r, if (r(IND) == 99.0) {print("BAD", IND); error("bad",0);});

        return r;
    }

    static void dxprintvalue(FILE* f, const double t) {
        fprintf(f,"%.6e\n",t);
    }

    static void dxprintvalue(FILE* f, const double_complex& t) {
        fprintf(f,"%.6e %.6e\n", t.real(), t.imag());
    }

    template <typename T, int NDIM>
    void plotdx(const Function<T,NDIM>& function,
                const char* filename, 
                const Tensor<double>& cell, 
                const std::vector<long>& npt,
                bool binary) {
        PROFILE_FUNC;
        MADNESS_ASSERT(NDIM<=6);
        const char* element[6] = {"lines","quads","cubes","cubes4D","cubes5D","cubes6D"};

        function.verify();
        World& world = const_cast< Function<T,NDIM>& >(function).world();
        FILE *f=0;
        if (world.rank() == 0) {
            f = fopen(filename, "w");
            if (!f) MADNESS_EXCEPTION("plotdx: failed to open the plot file", 0);

            fprintf(f,"object 1 class gridpositions counts ");
            for (int d=0; d<NDIM; d++) fprintf(f," %ld",npt[d]);
            fprintf(f,"\n");

            fprintf(f,"origin ");
            for (int d=0; d<NDIM; d++) fprintf(f, " %.6e", cell(d,0));
            fprintf(f,"\n");

            for (int d=0; d<NDIM; d++) {
                fprintf(f,"delta ");
                for (int c=0; c<d; c++) fprintf(f, " 0");
                double h = 0.0;
                if (npt[d]>1) h = (cell(d,1)-cell(d,0))/(npt[d]-1);
                fprintf(f," %.6e", h);
                for (int c=d+1; c<NDIM; c++) fprintf(f, " 0");
                fprintf(f,"\n");
            }
            fprintf(f,"\n");

            fprintf(f,"object 2 class gridconnections counts ");
            for (int d=0; d<NDIM; d++) fprintf(f," %ld",npt[d]);
            fprintf(f,"\n");
            fprintf(f, "attribute \"element type\" string \"%s\"\n", element[NDIM-1]);
            fprintf(f, "attribute \"ref\" string \"positions\"\n");
            fprintf(f,"\n");

            int npoint = 1;
            for (int d=0; d<NDIM; d++) npoint *= npt[d];
            const char* iscomplex = "";
            if (TensorTypeData<T>::iscomplex) iscomplex = "category complex";
            const char* isbinary = "";
            if (binary) isbinary = "binary";
            fprintf(f,"object 3 class array type double %s rank 0 items %d %s data follows\n",
                    iscomplex, npoint, isbinary);
        }

        world.gop.fence(); 
        Tensor<T> r = function.eval_cube(cell, npt);
        
        if (world.rank() == 0) {
            if (binary) {
                // This assumes that the values are double precision
                fflush(f);
                fwrite((void *) r.ptr(), sizeof(T), r.size, f);
                fflush(f);
            }
            else {
                for (IndexIterator it(npt); it; ++it) {
                    //fprintf(f,"%.6e\n",r(*it));
                    dxprintvalue(f,r(*it));
                }
            }
            fprintf(f,"\n");

            fprintf(f,"object \"%s\" class field\n",filename);
            fprintf(f,"component \"positions\" value 1\n");
            fprintf(f,"component \"connections\" value 2\n");
            fprintf(f,"component \"data\" value 3\n");
            fprintf(f,"\nend\n");
            fclose(f);
        }
        world.gop.fence(); 
    }

    template <int NDIM>
    void FunctionDefaults<NDIM>::set_defaults (World& world) {
            k = 7;
            thresh = 1e-5;
            initial_level = 2;
            max_refine_level = 30;
            truncate_mode = 0;
            refine = true;
            autorefine = true;
            debug = false;
            truncate_on_project = false;
            apply_randomize = false;
            project_randomize = false;
            bc = Tensor<int>(NDIM,2);
            cell = Tensor<double>(NDIM,2);
            cell(_,1) = 1.0;
            recompute_cell_info();

            //pmap = SharedPtr< WorldDCPmapInterface< Key<NDIM> > >(new WorldDCDefaultPmap< Key<NDIM> >(world));
            pmap = SharedPtr< WorldDCPmapInterface< Key<NDIM> > >(new MyPmap<NDIM>(world));
            //pmap = SharedPtr< WorldDCPmapInterface< Key<NDIM> > >(new SimpleMap< Key<NDIM> >(world));
        }


    // 
    // Below here we instantiate templates defined in this file
    //

    
    template <typename T, int NDIM>
    FunctionCommonData<T,NDIM> FunctionCommonData<T,NDIM>::data[MAXK+1];
    
    template <int NDIM> int FunctionDefaults<NDIM>::k;               
    template <int NDIM> double FunctionDefaults<NDIM>::thresh;       
    template <int NDIM> int FunctionDefaults<NDIM>::initial_level;   
    template <int NDIM> int FunctionDefaults<NDIM>::max_refine_level;
    template <int NDIM> int FunctionDefaults<NDIM>::truncate_mode; 
    template <int NDIM> bool FunctionDefaults<NDIM>::refine;         
    template <int NDIM> bool FunctionDefaults<NDIM>::autorefine;     
    template <int NDIM> bool FunctionDefaults<NDIM>::debug;          
    template <int NDIM> bool FunctionDefaults<NDIM>::truncate_on_project;          
    template <int NDIM> bool FunctionDefaults<NDIM>::apply_randomize;
    template <int NDIM> bool FunctionDefaults<NDIM>::project_randomize;
    template <int NDIM> Tensor<int> FunctionDefaults<NDIM>::bc;      
    template <int NDIM> Tensor<double> FunctionDefaults<NDIM>::cell;
    template <int NDIM> Tensor<double> FunctionDefaults<NDIM>::cell_width;
    template <int NDIM> Tensor<double> FunctionDefaults<NDIM>::rcell_width;
    template <int NDIM> double FunctionDefaults<NDIM>::cell_volume;
    template <int NDIM> double FunctionDefaults<NDIM>::cell_min_width;
    template <int NDIM> SharedPtr< WorldDCPmapInterface< Key<NDIM> > > FunctionDefaults<NDIM>::pmap;

#ifdef FUNCTION_INSTANTIATE_1
    template class FunctionDefaults<1>;
    template class Function<double, 1>;
    template class Function<std::complex<double>, 1>;
    template class FunctionImpl<double, 1>;
    template class FunctionImpl<std::complex<double>, 1>;
    template class FunctionCommonData<double, 1>;
    template class FunctionCommonData<double_complex, 1>;

    template void plotdx<double,1>(const Function<double,1>&, const char*, const Tensor<double>&, 
                                   const std::vector<long>&, bool binary);
    template void plotdx<double_complex,1>(const Function<double_complex,1>&, const char*, const Tensor<double>&, 
                                           const std::vector<long>&, bool binary);
#endif

#ifdef FUNCTION_INSTANTIATE_2
    template class FunctionDefaults<2>;
    template class Function<double, 2>;
    template class Function<std::complex<double>, 2>;
    template class FunctionImpl<double, 2>;
    template class FunctionImpl<std::complex<double>, 2>;
    template class FunctionCommonData<double, 2>;
    template class FunctionCommonData<double_complex, 2>;

    template void plotdx<double,2>(const Function<double,2>&, const char*, const Tensor<double>&, 
                                   const std::vector<long>&, bool binary);
    template void plotdx<double_complex,2>(const Function<double_complex,2>&, const char*, const Tensor<double>&, 
                                   const std::vector<long>&, bool binary);
#endif

#ifdef FUNCTION_INSTANTIATE_3
    template class FunctionDefaults<3>;
    template class Function<double, 3>;
    template class Function<std::complex<double>, 3>;
    template class FunctionImpl<double, 3>;
    template class FunctionImpl<std::complex<double>, 3>;
    template class FunctionCommonData<double, 3>;
    template class FunctionCommonData<double_complex, 3>;

    template void plotdx<double,3>(const Function<double,3>&, const char*, const Tensor<double>&, 
                                   const std::vector<long>&, bool binary);
    template void plotdx<double_complex,3>(const Function<double_complex,3>&, const char*, const Tensor<double>&, 
                                   const std::vector<long>&, bool binary);
#endif

#ifdef FUNCTION_INSTANTIATE_4
    template class FunctionDefaults<4>;
    template class Function<double, 4>;
    template class Function<std::complex<double>, 4>;
    template class FunctionImpl<double, 4>;
    template class FunctionImpl<std::complex<double>, 4>;
    template class FunctionCommonData<double, 4>;
    template class FunctionCommonData<double_complex, 4>;

    template void plotdx<double,4>(const Function<double,4>&, const char*, const Tensor<double>&, 
                                   const std::vector<long>&, bool binary);
    template void plotdx<double_complex,4>(const Function<double_complex,4>&, const char*, const Tensor<double>&, 
                                   const std::vector<long>&, bool binary);
#endif

#ifdef FUNCTION_INSTANTIATE_5
    template class FunctionDefaults<5>;
    template class Function<double, 5>;
    template class Function<std::complex<double>, 5>;
    template class FunctionImpl<double, 5>;
    template class FunctionImpl<std::complex<double>, 5>;
    template class FunctionCommonData<double, 5>;
    template class FunctionCommonData<double_complex, 5>;

    template void plotdx<double,5>(const Function<double,5>&, const char*, const Tensor<double>&, 
                                   const std::vector<long>&, bool binary);
    template void plotdx<double_complex,5>(const Function<double_complex,5>&, const char*, const Tensor<double>&, 
                                   const std::vector<long>&, bool binary);
#endif

#ifdef FUNCTION_INSTANTIATE_6
    template class FunctionDefaults<6>;
    template class Function<double, 6>;
    template class Function<std::complex<double>, 6>;
    template class FunctionImpl<double, 6>;
    template class FunctionImpl<std::complex<double>, 6>;
    template class FunctionCommonData<double, 6>;
    template class FunctionCommonData<double_complex, 6>;

    template void plotdx<double,6>(const Function<double,6>&, const char*, const Tensor<double>&, 
                                   const std::vector<long>&, bool binary);
    template void plotdx<double_complex,6>(const Function<double_complex,6>&, const char*, const Tensor<double>&, 
                                   const std::vector<long>&, bool binary);
#endif
}
