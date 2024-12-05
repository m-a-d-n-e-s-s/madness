/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

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
*/

#ifndef MADNESS_MRA_MRAIMPL_H__INCLUDED
#define MADNESS_MRA_MRAIMPL_H__INCLUDED

#ifndef MPRAIMPLX
#error "mraimpl.h should ONLY be included in one of the mraX.cc files (x=1..6)"
#endif

#include <memory>
#include <math.h>
#include <cmath>
#include <madness/world/world_object.h>
#include <madness/world/worlddc.h>
#include <madness/world/worldhashmap.h>
#include <madness/mra/function_common_data.h>

#include <madness/mra/funcimpl.h>
#include <madness/mra/displacements.h>

namespace std {
    template <typename T>
    bool isnan(const std::complex<T>& v) {
        return ::std::isnan(v.real()) || ::std::isnan(v.imag());
    }
}

/// \file mra/mraimpl.h
/// \brief Declaration and initialization of static data, some implementation, some instantiation

namespace madness {

    // Definition and initialization of FunctionDefaults static members
    // It cannot be an instance of FunctionFactory since we want to
    // set the defaults independent of the data type.

    template <typename T, std::size_t NDIM>
    void FunctionCommonData<T,NDIM>::_init_twoscale() {
        if (! two_scale_hg(k, &hg)) throw "failed to get twoscale coefficients";
        hgT = copy(transpose(hg));

        Slice sk(0,k-1), sk2(k,-1);
        hgsonly = copy(hg(Slice(0,k-1),_));

        h0 = copy(hg(sk,sk));
        h1 = copy(hg(sk,sk2));
        g0 = copy(hg(sk2,sk));
        g1 = copy(hg(sk2,sk2));

        h0T = copy(transpose(hg(sk,sk)));
        h1T = copy(transpose(hg(sk,sk2)));
        g0T = copy(transpose(hg(sk2,sk)));
        g1T = copy(transpose(hg(sk2,sk2)));

    }

    template <typename T, std::size_t NDIM>
    void FunctionCommonData<T,NDIM>::_init_quadrature
    (int k, int npt, Tensor<double>& quad_x, Tensor<double>& quad_w,
     Tensor<double>& quad_phi, Tensor<double>& quad_phiw, Tensor<double>& quad_phit) {
        quad_x = Tensor<double>(npt); // point
        quad_w = Tensor<double>(npt); // wheight
        quad_phi = Tensor<double>(npt,k);
        quad_phiw = Tensor<double>(npt,k);

        gauss_legendre(npt,0.0,1.0,quad_x.ptr(),quad_w.ptr());
        for (int mu=0; mu<npt; ++mu) {
            double phi[200];
            legendre_scaling_functions(quad_x(mu),k,phi);
            for (int j=0; j<k; ++j) {
                quad_phi(mu,j) = phi[j];
                quad_phiw(mu,j) = quad_w(mu)*phi[j];
            }
        }
        quad_phit = transpose(quad_phi);
    }


    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::verify_tree() const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        world.gop.fence();  // Make sure nothing is going on

        // Verify consistency of compression status, existence and size of coefficients,
        // and has_children() flag.
        for (typename dcT::const_iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
            const keyT& key = it->first;
            const nodeT& node = it->second;
            bool bad;

            if (is_compressed()) {
                if (node.has_children()) {
                    bad = (node.coeff().has_data()) and (node.coeff().dim(0) != 2*cdata.k);
                }
                else {
                    //                    bad = node.coeff().size() != 0;
                    bad = node.coeff().has_data();
                }
            }
            else {
                if (node.has_children()) {
                    //                    bad = node.coeff().size() != 0;
                    bad = node.coeff().has_data();
                }
                else {
                    bad = (node.coeff().has_data()) and ( node.coeff().dim(0) != cdata.k);
                }
            }

            if (bad) {
                print(world.rank(), "FunctionImpl: verify: INCONSISTENT TREE NODE, key =", key, ", node =", node,
                      ", dim[0] =",node.coeff().dim(0),", compressed =",is_compressed());
                std::cout.flush();
                MADNESS_EXCEPTION("FunctionImpl: verify: INCONSISTENT TREE NODE", 0);
            }
        }

        // Ensure that parents and children exist appropriately
        for (typename dcT::const_iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
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


    template <typename T, std::size_t NDIM>
    const std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >& FunctionImpl<T,NDIM>::get_pmap() const {
        return coeffs.get_pmap();
    }


    /// perform: this= alpha*f + beta*g, invoked by result

    /// f and g are reconstructed, so we can save on the compress operation,
    /// walk down the joint tree, and add leaf coefficients; effectively refines
    /// to common finest level.
    /// @param[in]  alpha   prefactor for f
    /// @param[in]  f       first addend
    /// @param[in]  beta    prefactor for g
    /// @param[in]  g       second addend
    /// @return     nothing, but leaves this's tree reconstructed and as sum of f and g
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::gaxpy_oop_reconstructed(const double alpha, const implT& f,
                                                       const double beta, const implT& g, const bool fence) {

        MADNESS_ASSERT(not f.is_compressed());
        MADNESS_ASSERT(not g.is_compressed());

        ProcessID owner = coeffs.owner(cdata.key0);
        if (world.rank() == owner) {

            CoeffTracker<T,NDIM> ff(&f);
            CoeffTracker<T,NDIM> gg(&g);

            typedef add_op coeff_opT;
            coeff_opT coeff_op(ff,gg,alpha,beta);
            typedef insert_op<T,NDIM> apply_opT;
            apply_opT apply_op(this);

            woT::task(world.rank(), &implT:: template forward_traverse<coeff_opT,apply_opT>,
                      coeff_op, apply_op, cdata.key0);

        }
        set_tree_state(reconstructed);
        if (fence) world.gop.fence();
    }

    /// Returns true if the function is compressed.
    template <typename T, std::size_t NDIM>
    bool FunctionImpl<T,NDIM>::is_compressed() const {
        return (tree_state==compressed);
    }

    /// Returns true if the function is reconstructed.
    template <typename T, std::size_t NDIM>
    bool FunctionImpl<T,NDIM>::is_reconstructed() const {
        return (tree_state==reconstructed);
    }

    /// Returns true if the function is redundant.
    template <typename T, std::size_t NDIM>
    bool FunctionImpl<T,NDIM>::is_redundant() const {
        return (tree_state==redundant);
    }

    template <typename T, std::size_t NDIM>
    bool FunctionImpl<T,NDIM>::is_nonstandard() const {
    	return (tree_state==nonstandard);
    }

    template <typename T, std::size_t NDIM>
    bool FunctionImpl<T,NDIM>::is_nonstandard_with_leaves() const {
    	return (tree_state==nonstandard_with_leaves);
    }

    template <typename T, std::size_t NDIM>
    bool FunctionImpl<T,NDIM>::is_on_demand() const {
    	return tree_state==on_demand;
    }

    template <typename T, std::size_t NDIM>
    bool FunctionImpl<T,NDIM>::has_leaves() const {
    	return (tree_state==nonstandard_with_leaves);
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::set_functor(const std::shared_ptr<FunctionFunctorInterface<T,NDIM> > functor1) {
        set_tree_state(on_demand);
//    	this->on_demand=true;
        functor=functor1;
    }

    template <typename T, std::size_t NDIM>
    std::shared_ptr<FunctionFunctorInterface<T,NDIM> > FunctionImpl<T,NDIM>::get_functor() {
        MADNESS_ASSERT(this->functor);
        return functor;
    }

    template <typename T, std::size_t NDIM>
    std::shared_ptr<FunctionFunctorInterface<T,NDIM> > FunctionImpl<T,NDIM>::get_functor() const {
        MADNESS_ASSERT(this->functor);
        return functor;
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::unset_functor() {
//        this->on_demand=false;
        set_tree_state(unknown);
        functor.reset();
    }

    template <typename T, std::size_t NDIM>
    TensorType FunctionImpl<T,NDIM>::get_tensor_type() const {return targs.tt;}

    template <typename T, std::size_t NDIM>
    TensorArgs FunctionImpl<T,NDIM>::get_tensor_args() const {return targs;}

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::set_tensor_args(const TensorArgs& t) {targs=t;}

    template <typename T, std::size_t NDIM>
    double FunctionImpl<T,NDIM>::get_thresh() const {return thresh;}

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::set_thresh(double value) {thresh = value;}

    template <typename T, std::size_t NDIM>
    bool FunctionImpl<T,NDIM>::get_autorefine() const {return autorefine;}

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::set_autorefine(bool value) {autorefine = value;}

    template <typename T, std::size_t NDIM>
    int FunctionImpl<T,NDIM>::get_k() const {return k;}

    template <typename T, std::size_t NDIM>
    const typename FunctionImpl<T,NDIM>::dcT& FunctionImpl<T,NDIM>::get_coeffs() const {return coeffs;}

    template <typename T, std::size_t NDIM>
    typename FunctionImpl<T,NDIM>::dcT& FunctionImpl<T,NDIM>::get_coeffs() {return coeffs;}

    template <typename T, std::size_t NDIM>
    const FunctionCommonData<T,NDIM>& FunctionImpl<T,NDIM>::get_cdata() const {return cdata;}

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::accumulate_timer(const double time) const {
        timer_accumulate.accumulate(time);
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::print_timer() const {
        if (world.rank()==0) {
            timer_accumulate.print("accumulate");
            timer_target_driven.print("target_driven");
            timer_lr_result.print("result2low_rank");
        }
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::reset_timer() {
        if (world.rank()==0) {
            timer_accumulate.reset();
            timer_target_driven.reset();
            timer_lr_result.reset();
        }
    }

    /// Truncate according to the threshold with optional global fence

    /// If thresh<=0 the default value of this->thresh is used
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::truncate(double tol, bool fence) {
        // Cannot put tol into object since it would make a race condition
        if (tol <= 0.0)
            tol = thresh;
        if (world.rank() == coeffs.owner(cdata.key0)) {
            if (is_compressed()) {
                truncate_spawn(cdata.key0,tol);
            } else {
                truncate_reconstructed_spawn(cdata.key0,tol);
            }
        }
        if (fence)
            world.gop.fence();
    }

    template <typename T, std::size_t NDIM>
    const typename FunctionImpl<T,NDIM>::keyT& FunctionImpl<T,NDIM>::key0() const {
        return cdata.key0;
    }

    /// Print a plane ("xy", "xz", or "yz") containing the point x to file

    /// works for all dimensions; we walk through the tree, and if a leaf node
    /// inside the sub-cell touches the plane we print it in pstricks format
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::print_plane(const std::string filename, const int xaxis, const int yaxis, const coordT& el2) {

        // get the local information
        Tensor<double> localinfo=print_plane_local(xaxis,yaxis,el2);

        // lump all the local information together, and gather on node0
        std::vector<Tensor<double> > localinfo_vec(1,localinfo);
        std::vector<Tensor<double> > printinfo=world.gop.concat0(localinfo_vec);
        world.gop.fence();

        // do the actual print
        if (world.rank()==0) do_print_plane(filename,printinfo,xaxis,yaxis,el2);
    }

    /// collect the data for a plot of the MRA structure locally on each node

    /// @param[in]	xaxis	the x-axis in the plot (can be any axis of the MRA box)
    /// @param[in]	yaxis	the y-axis in the plot (can be any axis of the MRA box)
    /// @param[in]	el2
    template <typename T, std::size_t NDIM>
    Tensor<double> FunctionImpl<T,NDIM>::print_plane_local(const int xaxis, const int yaxis, const coordT& el2) {
        coordT x_sim;
        user_to_sim<NDIM>(el2,x_sim);
        x_sim[0]+=1.e-10;

        // dimensions are: (# boxes)(hue, x lo left, y lo left, x hi right, y hi right)
        Tensor<double> plotinfo(coeffs.size(),5);
        long counter=0;

        // loop over local boxes, if the fit, add the info to the output tensor
        typename dcT::const_iterator end = coeffs.end();
        for (typename dcT::const_iterator it=coeffs.begin(); it!=end; ++it) {
            const keyT& key = it->first;
            const nodeT& node = it->second;

            // thisKeyContains ignores dim0 and dim1
            if (key.thisKeyContains(x_sim,xaxis,yaxis) and node.is_leaf() and (node.has_coeff())) {

                Level n=key.level();
                Vector<Translation,NDIM> l=key.translation();
                // get the diametral edges of the node in the plotting plane
                double scale=std::pow(0.5,double(n));
                double xloleft = scale*l[xaxis];
                double yloleft = scale*l[yaxis];
                double xhiright = scale*(l[xaxis]+1);
                double yhiright = scale*(l[yaxis]+1);

                // convert back to user coordinates
                Vector<double,4> user;
                user[0]=xloleft*FunctionDefaults<NDIM>::get_cell_width()[xaxis] + FunctionDefaults<NDIM>::get_cell()(xaxis,0);
                user[2]=xhiright*FunctionDefaults<NDIM>::get_cell_width()[xaxis] + FunctionDefaults<NDIM>::get_cell()(xaxis,0);
                user[1]=yloleft*FunctionDefaults<NDIM>::get_cell_width()[yaxis] + FunctionDefaults<NDIM>::get_cell()(yaxis,0);
                user[3]=yhiright*FunctionDefaults<NDIM>::get_cell_width()[yaxis] + FunctionDefaults<NDIM>::get_cell()(yaxis,0);


                //                    if ((xloleft<-5.0) or (yloleft<-5.0) or (xhiright>5.0) or (yhiright>5.0)) continue;
                if ((user[0]<-5.0) or (user[1]<-5.0) or (user[2]>5.0) or (user[3]>5.0)) continue;

                // do rank or do error
                double color=0.0;
                if (1) {

                    const double maxrank=40;
                    do_convert_to_color hue(maxrank,false);
                    color=hue(node.coeff().rank());
                } else {

                    // Make quadrature rule of higher order
                    const int npt = cdata.npt + 1;
                    Tensor<double> qx, qw, quad_phi, quad_phiw, quad_phit;
                    FunctionCommonData<T,NDIM>::_init_quadrature(k+1, npt, qx, qw, quad_phi, quad_phiw, quad_phit);
                    do_err_box< FunctionFunctorInterface<T,NDIM> > op(this, this->get_functor().get(), npt, qx, quad_phit, quad_phiw);

                    do_convert_to_color hue(1000.0,true);
                    double error=op(it);
                    error=sqrt(error);//*pow(2,key.level()*6);
                    color=hue(error);
                }

                plotinfo(counter,0)=color;
                plotinfo(counter,1)=user[0];
                plotinfo(counter,2)=user[1];
                plotinfo(counter,3)=user[2];
                plotinfo(counter,4)=user[3];
                ++counter;
            }
        }

        // shrink the info
        if (counter==0) plotinfo=Tensor<double>();
        else plotinfo=plotinfo(Slice(0,counter-1),Slice(_));
        return plotinfo;
    }

    /// print the MRA structure
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::do_print_plane(const std::string filename, std::vector<Tensor<double> > plotinfo,
                                              const int xaxis, const int yaxis, const coordT el2) {

        // invoke only on master node
        MADNESS_ASSERT(world.rank()==0);

        // prepare file
        FILE * pFile;
        pFile = fopen(filename.c_str(), "w");
        Tensor<double> cell=FunctionDefaults<NDIM>::get_cell();


        fprintf(pFile,"\\psset{unit=1cm}\n");
        fprintf(pFile,"\\begin{pspicture}(%4.2f,%4.2f)(%4.2f,%4.2f)\n",
                //            		cell(xaxis,0),cell(xaxis,1),cell(yaxis,0),cell(yaxis,1));
                -5.0,-5.0,5.0,5.0);
        fprintf(pFile,"\\pslinewidth=0.1pt\n");

        for (std::vector<Tensor<double> >::const_iterator it=plotinfo.begin(); it!=plotinfo.end(); ++it) {

            Tensor<double> localinfo=*it;
            if (localinfo.has_data()) {

                for (long i=0; i<localinfo.dim(0); ++i) {

                    fprintf(pFile,"\\newhsbcolor{mycolor}{%8.4f 1.0 0.7}\n",localinfo(i,0));
                    fprintf(pFile,"\\psframe["//linewidth=0.5pt,"
                            "fillstyle=solid,"
                            "fillcolor=mycolor]"
                            "(%12.8f,%12.8f)(%12.8f,%12.8f)\n",
                            localinfo(i,1),localinfo(i,2),localinfo(i,3),localinfo(i,4));
                }
            }
        }


        fprintf(pFile,"\\end{pspicture}\n");
        fclose(pFile);
    }

    /// print the grid (the roots of the quadrature of each leaf box)
    /// of this function in user xyz coordinates
    template <typename T, std::size_t NDIM>
    void  FunctionImpl<T,NDIM>::print_grid(const std::string filename) const {

        // get the local information
        std::vector<keyT> local_keys=local_leaf_keys();

        // lump all the local information together, and gather on node0
        std::vector<keyT> all_keys=world.gop.concat0(local_keys);
        world.gop.fence();

        // do the actual print
        if (world.rank()==0) do_print_grid(filename,all_keys);

    }

    /// return the keys of the local leaf boxes
    template <typename T, std::size_t NDIM>
    std::vector<typename FunctionImpl<T,NDIM>::keyT>  FunctionImpl<T,NDIM>::local_leaf_keys() const {

        // coeffs.size is maximum number of keys (includes internal keys)
        std::vector<keyT> keys(coeffs.size());

        // loop over local boxes, if they are leaf boxes add their quadrature roots
        // to the output tensor
        int i=0;
        typename dcT::const_iterator end = coeffs.end();
        for (typename dcT::const_iterator it=coeffs.begin(); it!=end; ++it) {
            const keyT& key = it->first;
            const nodeT& node = it->second;
            if (node.is_leaf()) keys[i++]=key;
        }

        // shrink the vector to number of leaf keys
        keys.resize(i);
        return keys;
    }

    /// print the grid in xyz format

    /// the quadrature points and the key information will be written to file,
    /// @param[in]	filename	where the quadrature points will be written to
    /// @param[in]	keys		all leaf keys
    template <typename T, std::size_t NDIM>
    void  FunctionImpl<T,NDIM>::do_print_grid(const std::string filename, const std::vector<keyT>& keys) const {
        // invoke only on master node
        MADNESS_ASSERT(world.rank()==0);

        // the quadrature points in simulation coordinates of the root node
        const Tensor<double> qx=cdata.quad_x;
        const size_t npt = qx.dim(0);

        // the number of coordinates (grid point tuples) per box ({x1},{x2},{x3},..,{xNDIM})
        long npoints=power<NDIM>(npt);
        // the number of boxes
        long nboxes=keys.size();

        // prepare file
        FILE * pFile;
        pFile = fopen(filename.c_str(), "w");

        fprintf(pFile,"%ld\n",npoints*nboxes);
        fprintf(pFile,"%ld points per box and %ld boxes \n",npoints,nboxes);

        // loop over all leaf boxes
        typename std::vector<keyT>::const_iterator key_it=keys.begin();
        for (key_it=keys.begin(); key_it!=keys.end(); ++key_it) {

            const keyT& key=*key_it;
            fprintf(pFile,"# key: %8d",key.level());
            for (size_t d=0; d<NDIM; d++) fprintf(pFile,"%8d",int(key.translation()[d]));
            fprintf(pFile,"\n");

            // this is borrowed from fcube
            const Vector<Translation,NDIM>& l = key.translation();
            const Level n = key.level();
            const double h = std::pow(0.5,double(n));
            coordT c; // will hold the point in user coordinates

            const Tensor<double>& cell_width = FunctionDefaults<NDIM>::get_cell_width();
            const Tensor<double>& cell = FunctionDefaults<NDIM>::get_cell();

            if (NDIM == 3) {
                for (size_t i=0; i<npt; ++i) {
                    c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                    for (size_t j=0; j<npt; ++j) {
                        c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                        for (size_t k=0; k<npt; ++k) {
                            c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                            // grid weights
                            //					            double scale = pow(0.5,0.5*NDIM*key.level())*
                            //					            		sqrt(FunctionDefaults<NDIM>::get_cell_volume());
                            //					            double w=cdata.quad_phiw[i]*cdata.quad_phiw[j]*cdata.quad_phiw[k];

                            fprintf(pFile,"%18.12f %18.12f %18.12f\n",c[0],c[1],c[2]);
                            //								fprintf(pFile,"%18.12e %18.12e %18.12e %18.12e\n",c[0],c[1],c[2],w*scale);
                        }
                    }
                }
            } else {
                MADNESS_EXCEPTION("only NDIM=3 in print_grid",0);
            }
        }
        fclose(pFile);
    }


    /// Returns the truncation threshold according to truncate_method
    template <typename T, std::size_t NDIM>
    double FunctionImpl<T,NDIM>::truncate_tol(double tol, const keyT& key) const {

        // RJH ... introduced max level here to avoid runaway
        // refinement due to truncation threshold going down to
        // intrinsic numerical error
        const int MAXLEVEL1 = 20; // 0.5**20 ~= 1e-6
        const int MAXLEVEL2 = 10; // 0.25**10 ~= 1e-6

        if (truncate_mode == 0) {
            return tol;
        }
        else if (truncate_mode == 1) {
            double L = FunctionDefaults<NDIM>::get_cell_min_width();
            return tol*std::min(1.0,pow(0.5,double(std::min(key.level(),MAXLEVEL1)))*L);
        }
        else if (truncate_mode == 2) {
            double L = FunctionDefaults<NDIM>::get_cell_min_width();
            return tol*std::min(1.0,pow(0.25,double(std::min(key.level(),MAXLEVEL2)))*L*L);
        }
        else if (truncate_mode == 3) {
            // similar to truncate mode 1, but with an additional factor to
            // account for an increased number of boxes in higher dimensions

            // here is our handwaving argument: this threshold will give each
            // FunctionNode an error of less than tol. The total error can
            // then be as high as sqrt(#nodes) * tol. Therefore in order to
            // account for higher dimensions: divide tol by about the root of
            // number of siblings (2^NDIM) that have a large error when we
            // refine along a deep branch of the tree. FAB
            //
            // Nope ... it can easily be as high as #nodes * tol.  The real
            // fix for this is an end-to-end error analysis of the larger
            // application and if desired to include this factor into the
            // threshold selected by the application. RJH
            const static double fac=1.0/std::pow(2,NDIM*0.5);
            tol*=fac;

            double L = FunctionDefaults<NDIM>::get_cell_min_width();
            return tol*std::min(1.0,pow(0.5,double(std::min(key.level(),MAXLEVEL1)))*L);

        } else {
            MADNESS_EXCEPTION("truncate_mode invalid",truncate_mode);
        }
    }

    /// Returns patch referring to coeffs of child in parent box
    template <typename T, std::size_t NDIM>
    std::vector<Slice> FunctionImpl<T,NDIM>::child_patch(const keyT& child) const {
        std::vector<Slice> s(NDIM);
        const Vector<Translation,NDIM>& l = child.translation();
        for (std::size_t i=0; i<NDIM; ++i)
            s[i] = cdata.s[l[i]&1]; // Lowest bit of translation
        return s;
    }

    /// Directly project parent NS coeffs to child NS coeffs

    template <typename T, std::size_t NDIM>
    typename FunctionImpl<T,NDIM>::coeffT FunctionImpl<T,NDIM>::parent_to_child_NS(
            const keyT& child, const keyT& parent, const coeffT& coeff) const {

        const implT* f=this;
        //        	MADNESS_ASSERT(coeff.tensor_type()==TT_FULL);
        coeffT result(f->cdata.v2k,coeff.tensor_type());

        // if the node for child is existent in f, and it is an internal node, we
        // automatically have the NS form; if it is a leaf node, we only have the
        // sum coeffs, so we take zero difference coeffs
        if (child==parent) {
            if (coeff.dim(0)==2*f->get_k()) result=coeff;		// internal node
            else if (coeff.dim(0)==f->get_k()) {			// leaf node
                result(f->cdata.s0)+=coeff;
            } else {
                MADNESS_EXCEPTION("confused k in parent_to_child_NS",1);
            }
        } else if (child.level()>parent.level()) {

            // parent and coeff should refer to a leaf node with sum coeffs only
            // b/c tree should be compressed with leaves kept.
            MADNESS_ASSERT(coeff.dim(0)==f->get_k());
            const coeffT scoeff=f->parent_to_child(coeff,parent,child);
            result(f->cdata.s0)+=scoeff;
        } else {
            MADNESS_EXCEPTION("confused keys in parent_to_child_NS",1);
        }
        return result;
    }

    /// truncate tree at a certain level
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::erase(const Level& max_level) {
        this->make_redundant(true);

        typename dcT::iterator end = coeffs.end();
        for (typename dcT::iterator it= coeffs.begin(); it!=end; ++it) {
            keyT key=it->first;
            nodeT& node=it->second;
            if (key.level()>max_level) coeffs.erase(key);
            if (key.level()==max_level) node.set_has_children(false);
        }
        this->undo_redundant(true);
    }


    /// Returns some asymmetry measure ... no comms
    template <typename T, std::size_t NDIM>
    double FunctionImpl<T,NDIM>::check_symmetry_local() const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        typedef Range<typename dcT::const_iterator> rangeT;
        return world.taskq.reduce<double,rangeT,do_check_symmetry_local>(rangeT(coeffs.begin(),coeffs.end()),
                                                                         do_check_symmetry_local(*this));
    }


    /// Refine multiple functions down to the same finest level

    /// @param[v] is the vector of functions we are refining.
    /// @param[key] is the current node.
    /// @param[c] is the vector of coefficients passed from above.
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::refine_to_common_level(const std::vector<FunctionImpl<T,NDIM>*>& v,
                                                      const std::vector<tensorT>& c,
                                                      const keyT key) {
        if (key == cdata.key0 && coeffs.owner(key)!=world.rank()) return;

        // First insert coefficients from above ... also get write accessors here
        std::unique_ptr<typename dcT::accessor[]> acc(new typename dcT::accessor[v.size()]);
        for (unsigned int i=0; i<c.size(); i++) {
            MADNESS_ASSERT(v[i]->coeffs.get_pmap() == coeffs.get_pmap());
            MADNESS_ASSERT(v[i]->coeffs.owner(key) == world.rank());
            bool exists = ! v[i]->coeffs.insert(acc[i],key);
            if (c[i].size()) {
                MADNESS_CHECK(!exists);
                acc[i]->second = nodeT(coeffT(c[i],targs),false);
            }
            else {
                MADNESS_ASSERT(exists);
            }
        }

        // If everyone has coefficients we are done
        bool done = true;
        for (unsigned int i=0; i<v.size(); i++) {
            done &= acc[i]->second.has_coeff();
        }

        if (!done) {
            // Those functions with coefficients need to be refined down
            std::vector<tensorT> d(v.size());
            for (unsigned int i=0; i<v.size(); i++) {
                if (acc[i]->second.has_coeff()) {
                    tensorT s(cdata.v2k);
                    //                        s(cdata.s0) = acc[i]->second.coeff()(___);
                    s(cdata.s0) = acc[i]->second.coeff().full_tensor_copy();
                    acc[i]->second.clear_coeff();
                    d[i] = unfilter(s);
                    acc[i]->second.set_has_children(true);
                }
            }

            // Loop thru children and pass down
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                const keyT& child = kit.key();
                std::vector<Slice> cp = child_patch(child);
                std::vector<tensorT> childc(v.size());
                for (unsigned int i=0; i<v.size(); i++) {
                    if (d[i].size()) childc[i] = copy(d[i](cp));
                }
                woT::task(coeffs.owner(child), &implT::refine_to_common_level, v, childc, child);
            }
        }
    }

    // horrifically non-scalable
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::put_in_box(ProcessID from, long nl, long ni) const {
        if (world.size()> 1000)
            throw "NO!";
        box_leaf[from] = nl;
        box_interior[from] = ni;
    }

    /// Prints summary of data distribution
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::print_info() const {
        if (world.size() >= 1000)
            return;
        for (int i=0; i<world.size(); ++i)
            box_leaf[i] = box_interior[i] == 0;
        world.gop.fence();
        long nleaf=0, ninterior=0;
        typename dcT::const_iterator end = coeffs.end();
        for (typename dcT::const_iterator it=coeffs.begin(); it!=end; ++it) {
            const nodeT& node = it->second;
            if (node.is_leaf())
                ++nleaf;
            else
                ++ninterior;
        }
        this->send(0, &implT::put_in_box, world.rank(), nleaf, ninterior);
        world.gop.fence();
        if (world.rank() == 0) {
            for (int i=0; i<world.size(); ++i) {
                printf("load: %5d %8ld %8ld\n", i, box_leaf[i], box_interior[i]);
            }
        }
        world.gop.fence();
    }

    template <typename T, std::size_t NDIM>
    bool FunctionImpl<T,NDIM>::noautorefine(const keyT& key, const tensorT& t) const {
        return false;
    }

    /// Returns true if this block of coeffs needs autorefining
    template <typename T, std::size_t NDIM>
    bool FunctionImpl<T,NDIM>::autorefine_square_test(const keyT& key, const nodeT& t) const {
        double lo, hi;
        tnorm(t.coeff().full_tensor_copy(), &lo, &hi);
        double test = 2*lo*hi + hi*hi;
        //print("autoreftest",key,thresh,truncate_tol(thresh, key),lo,hi,test);
        return test> truncate_tol(thresh, key);
    }


    /// is this the same as trickle_down() ?
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::sum_down_spawn(const keyT& key, const coeffT& s) {
        typename dcT::accessor acc;
        coeffs.insert(acc,key);
        nodeT& node = acc->second;
        coeffT& c = node.coeff();

        //print(key,"received",s.normf(),c.normf(),node.has_children());

        if (s.size() > 0) {
            if (c.size() > 0)
                c.gaxpy(1.0,s,1.0);
            else
                c = s;
        }

        if (node.has_children()) {
            coeffT d;
            if (c.has_data()) {
                d = coeffT(cdata.v2k,targs);
                d(cdata.s0) += c;
                d = unfilter(d);
                node.clear_coeff();
            }
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                coeffT ss;
                const keyT& child = kit.key();
                if (d.size() > 0) ss = copy(d(child_patch(child)));
                //print(key,"sending",ss.normf(),"to",child);
                woT::task(coeffs.owner(child), &implT::sum_down_spawn, child, ss);
            }
        }
        else {
            // Missing coeffs assumed to be zero
            if (c.size() <= 0) c = coeffT(cdata.vk,targs);
        }
    }

    /// After 1d push operator must sum coeffs down the tree to restore correct scaling function coefficients
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::sum_down(bool fence) {
        if (world.rank() == coeffs.owner(cdata.key0)) sum_down_spawn(cdata.key0, coeffT());

        if (fence) world.gop.fence();
    }


    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::forward_do_diff1(const DerivativeBase<T,NDIM>* D,
                                                const implT* f,
                                                const keyT& key,
                                                const std::pair<keyT,coeffT>& left,
                                                const std::pair<keyT,coeffT>& center,
                                                const std::pair<keyT,coeffT>& right) {
        D->forward_do_diff1(f,this,key,left,center,right);
    }


    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::do_diff1(const DerivativeBase<T,NDIM>* D,
                                        const implT* f,
                                        const keyT& key,
                                        const std::pair<keyT,coeffT>& left,
                                        const std::pair<keyT,coeffT>& center,
                                        const std::pair<keyT,coeffT>& right) {
        D->do_diff1(f,this,key,left,center,right);
    }


    // Called by result function to differentiate f
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::diff(const DerivativeBase<T,NDIM>* D, const implT* f, bool fence) {
        typedef std::pair<keyT,coeffT> argT;
        for (const auto& [key, node]: f->coeffs) {
            if (node.has_coeff()) {
                Future<argT> left  = D->find_neighbor(f, key,-1);
                argT center(key,node.coeff());
                Future<argT> right = D->find_neighbor(f, key, 1);
                world.taskq.add(*this, &implT::do_diff1, D, f, key, left, center, right, TaskAttributes::hipri());
            }
            else {
                coeffs.replace(key,nodeT(coeffT(),true)); // Empty internal node
            }
        }
        if (fence) world.gop.fence();
    }


    /// return the a std::pair<key, node>, which MUST exist
    template <typename T, std::size_t NDIM>
    std::pair<Key<NDIM>,ShallowNode<T,NDIM> > FunctionImpl<T,NDIM>::find_datum(keyT key) const {
        MADNESS_ASSERT(coeffs.probe(key));
        ShallowNode<T,NDIM> snode(coeffs.find(key).get()->second);
        return std::pair<Key<NDIM>,ShallowNode<T,NDIM> >(key,snode);
    }

    /// multiply the ket with a one-electron potential rr(1,2)= f(1,2)*g(1)

    /// @param[in]	val_ket	function values of f(1,2)
    /// @param[in]	val_pot	function values of g(1)
    /// @param[in]	particle	if 0 then g(1), if 1 then g(2)
    /// @return		the resulting function values
    template <typename T, std::size_t NDIM>
    typename FunctionImpl<T,NDIM>::coeffT FunctionImpl<T,NDIM>::multiply(const coeffT& val_ket,
            const coeffT& val_pot, int particle) const {

        MADNESS_ASSERT(particle==0 or particle==1);
        MADNESS_ASSERT(val_pot.is_full_tensor());
        MADNESS_ASSERT(val_ket.is_svd_tensor());

        std::vector<long> vkhalf=std::vector<long>(NDIM/2,cdata.vk[0]);
        tensorT ones=tensorT(vkhalf);
        ones=1.0;

        TensorArgs targs(-1.0,val_ket.tensor_type());
        coeffT pot12;
        if (particle==0) pot12=outer(val_pot.full_tensor(),ones,targs);
        else if (particle==1) pot12=outer(ones,val_pot.full_tensor(),targs);

        coeffT result=copy(val_ket);
        result.emul(pot12);

        return result;
    }


    /// given several coefficient tensors, assemble a result tensor

    /// the result looks like: 	(v(1,2) + v(1) + v(2)) |ket(1,2)>
    /// or 						(v(1,2) + v(1) + v(2)) |p(1) p(2)>
    /// i.e. coefficients for the ket and coefficients for the two particles are
    /// mutually exclusive. All potential terms are optional, just pass in empty coeffs.
    /// @param[in]	key			the key of the FunctionNode to which these coeffs belong
    /// @param[in]	cket		coefficients of the ket
    /// @param[in]	vpotential1	function values of the potential for particle 1
    /// @param[in]	vpotential2	function values of the potential for particle 2
    /// @param[in]	veri		function values for the 2-particle potential
    template <typename T, std::size_t NDIM>
    typename FunctionImpl<T,NDIM>::coeffT FunctionImpl<T,NDIM>::assemble_coefficients(
            const keyT& key, const coeffT& coeff_ket, const coeffT& vpotential1,
            const coeffT& vpotential2, const tensorT& veri) const {

        // take a shortcut if we are already done
        bool ket_only=(not (vpotential1.has_data() or vpotential2.has_data() or veri.has_data()));
        if (ket_only) return coeff_ket;

        // switch to values instead of coefficients
        coeffT val_ket=coeffs2values(key,coeff_ket);

        // the result tensor
        coeffT val_result=coeffT(val_ket.ndim(),val_ket.dims(),this->get_tensor_args().tt);
        coeffT coeff_result;

        // potential for particles 1 and 2, must be done in TT_2D
        if (vpotential1.has_data() or vpotential2.has_data()) {
            val_ket=val_ket.convert(TensorArgs(-1.0,TT_2D));
        }
        if (vpotential1.has_data()) val_result+=multiply(val_ket,vpotential1,0);
        if (vpotential2.has_data()) val_result+=multiply(val_ket,vpotential2,1);

        // values for eri: this must be done in full rank...
        if (veri.has_data()) {
            tensorT val_ket2=val_ket.full_tensor_copy().emul(veri);
            if (val_result.has_data()) val_ket2+=val_result.full_tensor_copy();
            // values2coeffs expensive (30%), coeffT() (relatively) cheap (8%)
            coeff_result=coeffT(values2coeffs(key,val_ket2),this->get_tensor_args());

        } else {

            // convert back to original tensor type
            val_ket=val_ket.convert(get_tensor_args());
            MADNESS_ASSERT(val_result.has_data());
            coeff_result=values2coeffs(key,val_result);
            coeff_result.reduce_rank(this->get_tensor_args().thresh);
        }

        return coeff_result;

    }

    /// Permute the dimensions of f according to map, result on this
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::mapdim(const implT& f, const std::vector<long>& map, bool fence) {

        PROFILE_MEMBER_FUNC(FunctionImpl);
        const_cast<implT*>(&f)->flo_unary_op_node_inplace(do_mapdim(map,*this),fence);

    }

    /// mirror the dimensions of f according to mirror, result on this
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::mirror(const implT& f, const std::vector<long>& mirrormap, bool fence) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        const_cast<implT*>(&f)->flo_unary_op_node_inplace(do_mirror(mirrormap,*this),fence);
    }

    /// map and mirror the translation index and the coefficients, result on this

    /// first map the dimensions, the mirror!
    /// this = mirror(map(f))
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::map_and_mirror(const implT& f, const std::vector<long>& map,
    		const std::vector<long>& mirror, bool fence) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        const_cast<implT*>(&f)->flo_unary_op_node_inplace(do_map_and_mirror(map,mirror,*this),fence);
    }



    /// take the average of two functions, similar to: this=0.5*(this+rhs)

    /// works in either basis and also in nonstandard form
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::average(const implT& rhs) {

        rhs.flo_unary_op_node_inplace(do_average(*this),true);
        this->scale_inplace(0.5,true);
        flo_unary_op_node_inplace(do_reduce_rank(targs),true);
    }

    /// change the tensor type of the coefficients in the FunctionNode

    /// @param[in]  targs   target tensor arguments (threshold and full/low rank)
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::change_tensor_type1(const TensorArgs& targs, bool fence) {
        flo_unary_op_node_inplace(do_change_tensor_type(targs,*this),fence);
    }

    /// reduce the rank of the coefficients tensors

    /// @param[in]  targs   target tensor arguments (threshold and full/low rank)
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::reduce_rank(const double thresh, bool fence) {
        flo_unary_op_node_inplace(do_reduce_rank(thresh),fence);
    }

    /// reduce the rank of the coefficients tensors

    /// @param[in]  targs   target tensor arguments (threshold and full/low rank)
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::chop_at_level(const int n, bool fence) {
        std::list<keyT> to_be_erased;
        for (auto it=coeffs.begin(); it!=coeffs.end(); ++it) {
            const keyT& key=it->first;
            nodeT& node=it->second;
            if (key.level()==n) node.set_is_leaf(true);
            if (key.level()>n) to_be_erased.push_back(key);
        }
        for (auto& key : to_be_erased) coeffs.erase(key);
    }


/// compute norm of s and d coefficients for all nodes

    /// @param[in]  targs   target tensor arguments (threshold and full/low rank)
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::compute_snorm_and_dnorm(bool fence) {
        //const auto& data=FunctionCommonData<T,NDIM>::get(get_k());
        flo_unary_op_node_inplace(
                do_compute_snorm_and_dnorm(cdata),fence);
    }


/// Transform sum coefficients at level n to sums+differences at level n-1

/// Given scaling function coefficients s[n][l][i] and s[n][l+1][i]
/// return the scaling function and wavelet coefficients at the
/// coarser level.  I.e., decompose Vn using Vn = Vn-1 + Wn-1.
    /// \code
    /// s_i = sum(j) h0_ij*s0_j + h1_ij*s1_j
    /// d_i = sum(j) g0_ij*s0_j + g1_ij*s1_j
    //  \endcode
    /// Returns a new tensor and has no side effects.  Works for any
    /// number of dimensions.
    ///
    /// No communication involved.
    template <typename T, std::size_t NDIM>
    typename FunctionImpl<T,NDIM>::tensorT FunctionImpl<T,NDIM>::filter(const tensorT& s) const {
        tensorT r(cdata.v2k,false);
        tensorT w(cdata.v2k,false);
        return fast_transform(s,cdata.hgT,r,w);
        //return transform(s,cdata.hgT);
    }

    template <typename T, std::size_t NDIM>
    typename FunctionImpl<T,NDIM>::coeffT FunctionImpl<T,NDIM>::filter(const coeffT& s) const {
        coeffT result=transform(s,cdata.hgT);
        return result;
    }

    ///  Transform sums+differences at level n to sum coefficients at level n+1

    ///  Given scaling function and wavelet coefficients (s and d)
    ///  returns the scaling function coefficients at the next finer
    ///  level.  I.e., reconstruct Vn using Vn = Vn-1 + Wn-1.
    ///  \code
    ///  s0 = sum(j) h0_ji*s_j + g0_ji*d_j
    ///  s1 = sum(j) h1_ji*s_j + g1_ji*d_j
    ///  \endcode
    ///  Returns a new tensor and has no side effects
    ///
    ///  If (sonly) ... then ss is only the scaling function coeff (and
    ///  assume the d are zero).  Works for any number of dimensions.
    ///
    /// No communication involved.
    template <typename T, std::size_t NDIM>
    typename FunctionImpl<T,NDIM>::tensorT FunctionImpl<T,NDIM>::unfilter(const tensorT& s) const {
        tensorT r(cdata.v2k,false);
        tensorT w(cdata.v2k,false);
        return fast_transform(s,cdata.hg,r,w);
        //return transform(s, cdata.hg);
    }

    template <typename T, std::size_t NDIM>
    typename FunctionImpl<T,NDIM>::coeffT FunctionImpl<T,NDIM>::unfilter(const coeffT& s) const {
        return transform(s,cdata.hg);
    }

    /// downsample the sum coefficients of level n+1 to sum coeffs on level n

    /// specialization of the filter method, will yield only the sum coefficients
    /// @param[in]  key key of level n
    /// @param[in]  v   vector of sum coefficients of level n+1
    /// @param[in]  args    TensorArguments for possible low rank approximations
    /// @return     sum coefficients on level n in full tensor format
    template <typename T, std::size_t NDIM>
    typename FunctionImpl<T,NDIM>::tensorT FunctionImpl<T,NDIM>::downsample(const keyT& key, const std::vector< Future<coeffT > >& v) const {

        tensorT result(cdata.vk);

        // the twoscale coefficients: for downsampling use h0/h1; see Alpert Eq (3.34a)
        const tensorT h[2] = {cdata.h0T, cdata.h1T};
        tensorT matrices[NDIM];

        // loop over all child nodes, transform and accumulate
        long i=0;
        for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {

            // get the appropriate twoscale coefficients for each dimension
            for (size_t ii=0; ii<NDIM; ++ii) matrices[ii]=h[kit.key().translation()[ii]%2];

            // transform and accumulate on the result
            result+=general_transform(v[i].get(),matrices).full_tensor_copy();

        }
        return result;
    }

    /// upsample the sum coefficients of level 1 to sum coeffs on level n+1

    /// specialization of the unfilter method, will transform only the sum coefficients
    /// @param[in]  key     key of level n+1
    /// @param[in]  coeff   sum coefficients of level n (does NOT belong to key!!)
    /// @param[in]  args    TensorArguments for possible low rank approximations
    /// @return     sum     coefficients on level n+1
    template <typename T, std::size_t NDIM>
    typename FunctionImpl<T,NDIM>::coeffT FunctionImpl<T,NDIM>::upsample(const keyT& key, const coeffT& coeff) const {

        // the twoscale coefficients: for upsampling use h0/h1; see Alpert Eq (3.35a/b)
        // note there are no difference coefficients; if you want that use unfilter
        const tensorT h[2] = {cdata.h0, cdata.h1};
        tensorT matrices[NDIM];

        // get the appropriate twoscale coefficients for each dimension
        for (size_t ii=0; ii<NDIM; ++ii) matrices[ii]=h[key.translation()[ii]%2];

        // transform and accumulate on the result
        const coeffT result=general_transform(coeff,matrices);
        return result;
    }


    /// Projects old function into new basis (only in reconstructed form)
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::project(const implT& old, bool fence) {
        long kmin = std::min(cdata.k,old.cdata.k);
        std::vector<Slice> s(NDIM,Slice(0,kmin-1));
        typename dcT::const_iterator end = old.coeffs.end();
        for (typename dcT::const_iterator it=old.coeffs.begin(); it!=end; ++it) {
            const keyT& key = it->first;
            const nodeT& node = it->second;
            if (node.has_coeff()) {
                coeffT c(cdata.vk,targs);
                c(s) += node.coeff()(s);
                coeffs.replace(key,nodeT(c,false));
            }
            else {
                coeffs.replace(key,nodeT(coeffT(),true));
            }
        }
        if (fence)
            world.gop.fence();
    }

    template <typename T, std::size_t NDIM>
    bool FunctionImpl<T,NDIM>::exists_and_has_children(const keyT& key) const {
        return coeffs.probe(key) && coeffs.find(key).get()->second.has_children();
    }

    template <typename T, std::size_t NDIM>
    bool FunctionImpl<T,NDIM>::exists_and_is_leaf(const keyT& key) const {
        return coeffs.probe(key) && (not coeffs.find(key).get()->second.has_children());
    }


    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::broaden_op(const keyT& key, const std::vector< Future <bool> >& v) {
        for (unsigned int i=0; i<v.size(); ++i) {
            if (v[i]) {
                refine_op(true_refine_test(), key);
                break;
            }
        }
    }

    // For each local node sets value of norm tree, snorm and dnorm to 0.0
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::zero_norm_tree() {
        typename dcT::iterator end = coeffs.end();
        for (typename dcT::iterator it=coeffs.begin(); it!=end; ++it) {
            it->second.set_norm_tree(0.0);
            it->second.set_snorm(0.0);
            it->second.set_dnorm(0.0);
        }
    }

    // Broaden tree
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::broaden(const std::array<bool, NDIM>& is_periodic, bool fence) {
        typename dcT::iterator end = coeffs.end();
        for (typename dcT::iterator it=coeffs.begin(); it!=end; ++it) {
            const keyT& key = it->first;
            typename dcT::accessor acc;
            const auto found = coeffs.find(acc,key);
            MADNESS_CHECK(found);
            nodeT& node = acc->second;
            if (node.has_coeff() &&
                node.get_norm_tree() != -1.0 &&
                node.coeff().normf() >= truncate_tol(thresh,key)) {

                node.set_norm_tree(-1.0); // Indicates already broadened or result of broadening/refining

                //int ndir = std::pow(3,NDIM);
                int ndir = static_cast<int>(std::pow(static_cast<double>(3), static_cast<int>(NDIM)));
                std::vector< Future <bool> > v = future_vector_factory<bool>(ndir);
                keyT neigh;
                int i=0;
                for (HighDimIndexIterator it(NDIM,3); it; ++it) {
                    Vector<Translation,NDIM> l(*it);
                    for (std::size_t d=0; d<NDIM; ++d) {
                        const int odd = key.translation()[d] & 0x1L; // 1 if odd, 0 if even
                        l[d] -= 1; // (0,1,2) --> (-1,0,1)
                        if (l[d] == -1)
                            l[d] = -1-odd;
                        else if (l[d] ==  1)
                            l[d] = 2 - odd;
                    }
                    keyT neigh = neighbor(key, keyT(key.level(),l), is_periodic);

                    if (neigh.is_valid()) {
                        v[i++] = this->task(coeffs.owner(neigh), &implT::exists_and_has_children, neigh);
                    }
                    else {
                        v[i++].set(false);
                    }
                }
                woT::task(world.rank(), &implT::broaden_op, key, v);
            }
        }
        // Reset value of norm tree so that can repeat broadening
        if (fence) {
            world.gop.fence();
            zero_norm_tree();
            world.gop.fence();
        }
    }

    /// sum all the contributions from all scales after applying an operator in mod-NS form
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::trickle_down(bool fence) {
        set_tree_state(reconstructed);
        if (world.rank() == coeffs.owner(cdata.key0))
            woT::task(world.rank(), &implT::trickle_down_op, cdata.key0,coeffT());
        if (fence) world.gop.fence();
    }

    /// sum all the contributions from all scales after applying an operator in mod-NS form

    /// cf reconstruct_op
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::trickle_down_op(const keyT& key, const coeffT& s) {
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
            coeffs.replace(key,nodeT(coeffT(),false));
            it = coeffs.find(key).get();
        }
        nodeT& node = it->second;

        // The integral operator will correctly connect interior nodes
        // to children but may leave interior nodes without coefficients
        // ... but they still need to sum down so just give them zeros
        if (node.coeff().has_no_data()) node.coeff()=coeffT(cdata.vk,targs);

        //            if (node.has_children() || node.has_coeff()) { // Must allow for inconsistent state from transform, etc.
        if (node.has_children()) { // Must allow for inconsistent state from transform, etc.
            coeffT d = node.coeff();
            if (key.level() > 0) d += s; // -- note accumulate for NS summation
            node.clear_coeff();
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                const keyT& child = kit.key();
                coeffT ss= upsample(child,d);
                ss.reduce_rank(thresh);
                PROFILE_BLOCK(recon_send);
                woT::task(coeffs.owner(child), &implT::trickle_down_op, child, ss);
            }
        }
        else {
            node.coeff()+=s;
            node.coeff().reduce_rank(thresh);
        }
    }

    /// change the tree state of this function, might or might not respect fence!
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::change_tree_state(const TreeState finalstate, bool fence) {

        TreeState current_state=get_tree_state();
        if (current_state==finalstate) return;

        // very special case
        if (get_tree_state()==nonstandard_after_apply) {
            MADNESS_CHECK(finalstate==reconstructed);
            reconstruct(fence);
            return;
        }
        MADNESS_CHECK_THROW(current_state!=TreeState::nonstandard_after_apply,"unknown tree state");
        bool must_fence=false;

        if (finalstate==reconstructed) {
            if (current_state==reconstructed) return;
            if (current_state==compressed) reconstruct(fence);
            if (current_state==nonstandard) reconstruct(fence);
            if (current_state==nonstandard_with_leaves) remove_internal_coefficients(fence);
            if (current_state==redundant) remove_internal_coefficients(fence);
            set_tree_state(reconstructed);
        } else if (finalstate==compressed) {
            if (current_state==reconstructed) compress(compressed,fence);
            if (current_state==compressed) return;
            if (current_state==nonstandard) standard(fence);
            if (current_state==nonstandard_with_leaves) standard(fence);
            if (current_state==redundant) {
                remove_internal_coefficients(true);
                must_fence=true;
                set_tree_state(reconstructed);
                compress(compressed,fence);
            }
            set_tree_state(compressed);
        } else if (finalstate==nonstandard) {
            if (current_state==reconstructed) compress(nonstandard,fence);
            if (current_state==compressed) {
                reconstruct(true);
                must_fence=true;
                compress(nonstandard,fence);
            }
            if (current_state==nonstandard) return;
            if (current_state==nonstandard_with_leaves) remove_leaf_coefficients(fence);
            if (current_state==redundant) {
                remove_internal_coefficients(true);
                must_fence=true;
                set_tree_state(reconstructed);
                compress(nonstandard,fence);
            }
            set_tree_state(nonstandard);
        } else if (finalstate==nonstandard_with_leaves) {
            if (current_state==reconstructed) compress(nonstandard_with_leaves,fence);
            if (current_state==compressed) {
                reconstruct(true);
                must_fence=true;
                compress(nonstandard_with_leaves,fence);
            }
            if (current_state==nonstandard) {
                standard(true);
                must_fence=true;
                reconstruct(true);
                compress(nonstandard_with_leaves,fence);
            }
            if (current_state==nonstandard_with_leaves) return;
            if (current_state==redundant) {
                remove_internal_coefficients(true);
                must_fence=true;
                set_tree_state(reconstructed);
                compress(nonstandard_with_leaves,fence);
            }
            set_tree_state(nonstandard_with_leaves);
        } else if (finalstate==redundant) {
            if (current_state==reconstructed) make_redundant(fence);
            if (current_state==compressed) {
                reconstruct(true);
                must_fence=true;
                make_redundant(fence);
            }
            if (current_state==nonstandard) {
                standard(true);
                must_fence=true;
                reconstruct(true);
                make_redundant(fence);
            }
            if (current_state==nonstandard_with_leaves) {
                remove_internal_coefficients(true);
                must_fence=true;
                set_tree_state(reconstructed);
                make_redundant(fence);
            }
            if (current_state==redundant) return;
            set_tree_state(redundant);
        } else {
            MADNESS_EXCEPTION("unknown/unsupported final tree state",1);
        }
        if (must_fence and world.rank()==0) {
            print("could not respect fence in change_tree_state");
        }
        if (fence && VERIFY_TREE) verify_tree(); // Must be after in case nonstandard
        return;
    }



    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::reconstruct(bool fence) {

        if (is_reconstructed()) return;

        if (is_redundant() or is_nonstandard_with_leaves()) {
            set_tree_state(reconstructed);
    		this->remove_internal_coefficients(fence);
    	} else if (is_compressed() or tree_state==nonstandard_after_apply) {
            // Must set true here so that successive calls without fence do the right thing
            set_tree_state(reconstructed);
            if (world.rank() == coeffs.owner(cdata.key0))
                woT::task(world.rank(), &implT::reconstruct_op, cdata.key0,coeffT(), true);
        } else if (is_nonstandard()) {
            // Must set true here so that successive calls without fence do the right thing
            set_tree_state(reconstructed);
            if (world.rank() == coeffs.owner(cdata.key0))
                woT::task(world.rank(), &implT::reconstruct_op, cdata.key0,coeffT(), false);
    	} else {
            MADNESS_EXCEPTION("cannot reconstruct this tree",1);
        }
        if (fence) world.gop.fence();

    }

    /// compress the wave function

    /// after application there will be sum coefficients at the root level,
    /// and difference coefficients at all other levels; furthermore:
    /// @param[in] nonstandard	keep sum coeffs at all other levels, except leaves
    /// @param[in] keepleaves	keep sum coeffs (but no diff coeffs) at leaves
    /// @param[in] redundant    keep only sum coeffs at all levels, discard difference coeffs
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::compress(const TreeState newstate, bool fence) {
        MADNESS_CHECK_THROW(is_reconstructed(),"impl::compress wants a reconstructe tree");
        // Must set true here so that successive calls without fence do the right thing
        set_tree_state(newstate);
        bool keepleaves1=(tree_state==nonstandard_with_leaves) or (tree_state==redundant);
        bool nonstandard1=(tree_state==nonstandard) or (tree_state==nonstandard_with_leaves);
        bool redundant1=(tree_state==redundant);

        if (world.rank() == coeffs.owner(cdata.key0)) {

            compress_spawn(cdata.key0, nonstandard1, keepleaves1, redundant1);
        }
        if (fence)
            world.gop.fence();
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::remove_internal_coefficients(const bool fence) {
        flo_unary_op_node_inplace(remove_internal_coeffs(),fence);
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::remove_leaf_coefficients(const bool fence) {
        flo_unary_op_node_inplace(remove_leaf_coeffs(),fence);
    }

    /// convert this to redundant, i.e. have sum coefficients on all levels
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::make_redundant(const bool fence) {

        // fast return if possible
        if (is_redundant()) return;
        MADNESS_CHECK_THROW(is_reconstructed(),"impl::make_redundant() wants a reconstructed tree");
        compress(redundant,fence);
    }

    /// convert this from redundant to standard reconstructed form
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::undo_redundant(const bool fence) {
        MADNESS_CHECK_THROW(is_redundant(),"impl::undo_redundant() wants a redundant tree");
        set_tree_state(reconstructed);
        flo_unary_op_node_inplace(remove_internal_coeffs(),fence);
    }


    /// compute for each FunctionNode the norm of the function inside that node
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::norm_tree(bool fence) {
        if (world.rank() == coeffs.owner(cdata.key0))
            norm_tree_spawn(cdata.key0);
        if (fence)
            world.gop.fence();
    }

    template <typename T, std::size_t NDIM>
    double FunctionImpl<T,NDIM>::norm_tree_op(const keyT& key, const std::vector< Future<double> >& v) {
        //PROFILE_MEMBER_FUNC(FunctionImpl);
        double sum = 0.0;
        int i=0;
        for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
            double value = v[i].get();
            sum += value*value;
        }
        sum = sqrt(sum);
        coeffs.task(key, &nodeT::set_norm_tree, sum); // why a task? because send is deprecated to keep comm thread free
        //if (key.level() == 0) std::cout << "NORM_TREE_TOP " << sum << "\n";
        return sum;
    }

    template <typename T, std::size_t NDIM>
    Future<double> FunctionImpl<T,NDIM>::norm_tree_spawn(const keyT& key) {
        nodeT& node = coeffs.find(key).get()->second;
        if (node.has_children()) {
            std::vector< Future<double> > v = future_vector_factory<double>(1<<NDIM);
            int i=0;
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
                v[i] = woT::task(coeffs.owner(kit.key()), &implT::norm_tree_spawn, kit.key());
            }
            return woT::task(world.rank(),&implT::norm_tree_op, key, v);
        }
        else {
            //                return Future<double>(node.coeff().normf());
            const double norm=node.coeff().normf();
            // invoked locally anyways
            node.set_norm_tree(norm);
            return Future<double>(norm);
        }
    }

    /// truncate using a tree in reconstructed form

    /// must be invoked where key is local
    template <typename T, std::size_t NDIM>
    Future<typename FunctionImpl<T,NDIM>::coeffT> FunctionImpl<T,NDIM>::truncate_reconstructed_spawn(const keyT& key, const double tol) {
        MADNESS_ASSERT(coeffs.probe(key));
        nodeT& node = coeffs.find(key).get()->second;

        // if this is a leaf node just return the sum coefficients
        if (not node.has_children()) return Future<coeffT>(node.coeff());

        // if this is an internal node, wait for all the children's sum coefficients
        // and use them to determine if the children can be removed
        std::vector<Future<coeffT> > v = future_vector_factory<coeffT>(1<<NDIM);
        int i=0;
        for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
            v[i] = woT::task(coeffs.owner(kit.key()), &implT::truncate_reconstructed_spawn, kit.key(),tol,TaskAttributes::hipri());
        }

        // will return (possibly empty) sum coefficients
        return woT::task(world.rank(),&implT::truncate_reconstructed_op,key,v,tol,TaskAttributes::hipri());

    }

    /// given the sum coefficients of all children, truncate or not

    /// @return     new sum coefficients (empty if internal, not empty, if new leaf); might delete its children
    template <typename T, std::size_t NDIM>
    typename FunctionImpl<T,NDIM>::coeffT FunctionImpl<T,NDIM>::truncate_reconstructed_op(const keyT& key, const std::vector< Future<coeffT > >& v, const double tol) {

        MADNESS_ASSERT(coeffs.probe(key));

        // the sum coefficients might be empty, which means they come from an internal node
        // and we must not truncate; so just return empty coeffs again
        for (size_t i=0; i<v.size(); ++i) if (v[i].get().has_no_data()) return coeffT();

        // do not truncate below level 1
        if (key.level()<2) return coeffT();

        // compute the wavelet coefficients from the child nodes
        typename dcT::accessor acc;
        const auto found = coeffs.find(acc, key);
        MADNESS_CHECK(found);
        int i=0;
        tensorT d(cdata.v2k);
        for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
            //                d(child_patch(kit.key())) += v[i].get();
            d(child_patch(kit.key())) += v[i].get().full_tensor_copy();
        }

        d = filter(d);
        tensorT s=copy(d(cdata.s0));
        d(cdata.s0) = 0.0;
        const double error=d.normf();

        nodeT& node = coeffs.find(key).get()->second;

        if (error < truncate_tol(tol,key)) {
            node.set_has_children(false);
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                coeffs.erase(kit.key());
            }
            // "replace" children with new sum coefficients
            coeffT ss=coeffT(s,targs);
            acc->second.set_coeff(ss);
            return ss;
        } else {
            return coeffT();
        }
    }

    /// calculate the wavelet coefficients using the sum coefficients of all child nodes

    /// @param[in] key 	this's key
    /// @param[in] v 	sum coefficients of the child nodes
    /// @param[in] nonstandard  keep the sum coefficients with the wavelet coefficients
    /// @param[in] redundant    keep only the sum coefficients, discard the wavelet coefficients
    /// @return 		the sum coefficients
    template <typename T, std::size_t NDIM>
    std::pair<typename FunctionImpl<T,NDIM>::coeffT,double> FunctionImpl<T,NDIM>::compress_op(const keyT& key,
    		const std::vector< Future<std::pair<coeffT,double> > >& v, bool nonstandard1) {
        //PROFILE_MEMBER_FUNC(FunctionImpl);

        double cpu0=cpu_time();
        // Copy child scaling coeffs into contiguous block
        tensorT d(cdata.v2k);
        //            coeffT d(cdata.v2k,targs);
        int i=0;
        double norm_tree2=0.0;
        for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
            //                d(child_patch(kit.key())) += v[i].get();
            d(child_patch(kit.key())) += v[i].get().first.full_tensor_copy();
            norm_tree2+=v[i].get().second*v[i].get().second;
        }

        d = filter(d);
        double cpu1=cpu_time();
        timer_filter.accumulate(cpu1-cpu0);
        cpu0=cpu1;

        typename dcT::accessor acc;
        const auto found = coeffs.find(acc, key);
        MADNESS_CHECK(found);
        MADNESS_CHECK_THROW(!acc->second.has_coeff(),"compress_op: existing coeffs where there should be none");

        // tighter thresh for internal nodes
        TensorArgs targs2=targs;
        targs2.thresh*=0.1;

        // need the deep copy for contiguity
        coeffT ss=coeffT(copy(d(cdata.s0)));
        double snorm=ss.normf();

        if (key.level()> 0 && !nonstandard1) d(cdata.s0) = 0.0;

        coeffT dd=coeffT(d,targs2);
        double dnorm=dd.normf();
        double norm_tree=sqrt(norm_tree2);

        acc->second.set_snorm(snorm);
        acc->second.set_dnorm(dnorm);
        acc->second.set_norm_tree(norm_tree);

        acc->second.set_coeff(dd);
        cpu1=cpu_time();
        timer_compress_svd.accumulate(cpu1-cpu0);

        // return sum coefficients
        return std::make_pair(ss,snorm);
    }

    /// similar to compress_op, but insert only the sum coefficients in the tree

    /// also sets snorm, dnorm and norm_tree for all nodes
    /// @param[in] key  this's key
    /// @param[in] v    sum coefficients of the child nodes
    /// @return         the sum coefficients
    template <typename T, std::size_t NDIM>
    std::pair<typename FunctionImpl<T,NDIM>::coeffT,double>
            FunctionImpl<T,NDIM>::make_redundant_op(const keyT& key, const std::vector< Future<std::pair<coeffT,double> > >& v) {

        tensorT d(cdata.v2k);
        int i=0;
        double norm_tree2=0.0;
        for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
            d(child_patch(kit.key())) += v[i].get().first.full_tensor_copy();
            norm_tree2+=v[i].get().second*v[i].get().second;
        }
        d = filter(d);
        double norm_tree=sqrt(norm_tree2);

        // tighter thresh for internal nodes
        TensorArgs targs2=targs;
        targs2.thresh*=0.1;

        // need the deep copy for contiguity
        coeffT s=coeffT(copy(d(cdata.s0)),targs2);
        d(cdata.s0)=0.0;
        double dnorm=d.normf();
        double snorm=s.normf();

        typename dcT::accessor acc;
        const auto found = coeffs.find(acc, key);
        MADNESS_CHECK(found);

        acc->second.set_coeff(s);
        acc->second.set_dnorm(dnorm);
        acc->second.set_snorm(snorm);
        acc->second.set_norm_tree(norm_tree);

        // return sum coefficients
        return std::make_pair(s,norm_tree);
    }

    /// Changes non-standard compressed form to standard compressed form
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::standard(bool fence) {

        if (is_compressed()) return;
        set_tree_state(compressed);
        flo_unary_op_node_inplace(do_standard(this),fence);
//        make_nonstandard = false;
    }


    /// after apply we need to do some cleanup;

    /// forces fence
    template <typename T, std::size_t NDIM>
    double FunctionImpl<T,NDIM>::finalize_apply() {
        bool print_timings=false;
        bool printme=(world.rank()==0 and print_timings);
        TensorArgs tight_args(targs);
        tight_args.thresh*=0.01;
        double begin=wall_time();
        double begin1=wall_time();
        flo_unary_op_node_inplace(do_consolidate_buffer(tight_args),true);
        double end1=wall_time();
        if (printme) printf("time in consolidate_buffer    %8.4f\n",end1-begin1);


        // reduce the rank of the final nodes, leave full tensors unchanged
        //            flo_unary_op_node_inplace(do_reduce_rank(tight_args.thresh),true);
        begin1=wall_time();
        flo_unary_op_node_inplace(do_reduce_rank(targs),true);
        end1=wall_time();
        if (printme) printf("time in do_reduce_rank        %8.4f\n",end1-begin1);

        // change TT_FULL to low rank
        begin1=wall_time();
        flo_unary_op_node_inplace(do_change_tensor_type(targs,*this),true);
        end1=wall_time();
        if (printme) printf("time in do_change_tensor_type %8.4f\n",end1-begin1);

        // truncate leaf nodes to avoid excessive tree refinement
        begin1=wall_time();
        flo_unary_op_node_inplace(do_truncate_NS_leafs(this),true);
        end1=wall_time();
        if (printme) printf("time in do_truncate_NS_leafs  %8.4f\n",end1-begin1);

        double end=wall_time();
        double elapsed=end-begin;
        set_tree_state(nonstandard_after_apply);
        world.gop.fence();
        return elapsed;
    }


    /// after summing up we need to do some cleanup;

    /// forces fence
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::finalize_sum() {
        world.gop.fence();
        flo_unary_op_node_inplace(do_consolidate_buffer(get_tensor_args()), true);
        sum_down(true);
        set_tree_state(reconstructed);
    }

    /// Returns the square of the local norm ... no comms
    template <typename T, std::size_t NDIM>
    double FunctionImpl<T,NDIM>::norm2sq_local() const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        typedef Range<typename dcT::const_iterator> rangeT;
        return world.taskq.reduce<double,rangeT,do_norm2sq_local>(rangeT(coeffs.begin(),coeffs.end()),
                                                                  do_norm2sq_local());
    }




    /// Returns the maximum local depth of the tree ... no communications.
    template <typename T, std::size_t NDIM>
    std::size_t FunctionImpl<T,NDIM>::max_local_depth() const {
        std::size_t maxdepth = 0;
        typename dcT::const_iterator end = coeffs.end();
        for (typename dcT::const_iterator it=coeffs.begin(); it!=end; ++it) {
            std::size_t N = (std::size_t) it->first.level();
            if (N> maxdepth)
                maxdepth = N;
        }
        return maxdepth;
    }


    /// Returns the maximum depth of the tree ... collective ... global sum/broadcast
    template <typename T, std::size_t NDIM>
    std::size_t FunctionImpl<T,NDIM>::max_depth() const {
        std::size_t maxdepth  = max_local_depth();
        world.gop.max(maxdepth);
        return maxdepth;
    }

    /// Returns the max number of nodes on a processor
    template <typename T, std::size_t NDIM>
    std::size_t FunctionImpl<T,NDIM>::max_nodes() const {
        std::size_t maxsize = 0;
        maxsize = coeffs.size();
        world.gop.max(maxsize);
        return maxsize;
    }

    /// Returns the min number of nodes on a processor
    template <typename T, std::size_t NDIM>
    std::size_t FunctionImpl<T,NDIM>::min_nodes() const {
        std::size_t minsize = 0;
        minsize = coeffs.size();
        world.gop.min(minsize);
        return minsize;
    }

    /// Returns the size of the tree structure of the function ... collective global sum
    template <typename T, std::size_t NDIM>
    std::size_t FunctionImpl<T,NDIM>::tree_size() const {
        std::size_t sum = 0;
        sum = coeffs.size();
        world.gop.sum(sum);
        return sum;
    }

    /// Returns the number of coefficients in the function ... collective global sum
    template <typename T, std::size_t NDIM>
    std::size_t FunctionImpl<T,NDIM>::size() const {
        std::size_t sum = 0;
#if 1
        typename dcT::const_iterator end = coeffs.end();
        for (typename dcT::const_iterator it=coeffs.begin(); it!=end; ++it) {
            const nodeT& node = it->second;
            if (node.has_coeff())
                sum+=node.size();
        }
        //            print("proc",world.rank(),sum);
#else
        typename dcT::const_iterator end = coeffs.end();
        for (typename dcT::const_iterator it=coeffs.begin(); it!=end; ++it) {
            const nodeT& node = it->second;
            if (node.has_coeff())
                ++sum;
        }
        if (is_compressed())
            for (std::size_t i=0; i<NDIM; ++i)
                sum *= 2*cdata.k;
        else
            for (std::size_t i=0; i<NDIM; ++i)
                sum *= cdata.k;
#endif
        world.gop.sum(sum);

        return sum;
    }

    /// Returns the number of coefficients in the function ... collective global sum
    template <typename T, std::size_t NDIM>
    std::size_t FunctionImpl<T,NDIM>::real_size() const {
        std::size_t sum = coeffs.size() * (sizeof(keyT) + sizeof(nodeT));
        typename dcT::const_iterator end = coeffs.end();
        for (typename dcT::const_iterator it=coeffs.begin(); it!=end; ++it) {
            const nodeT& node = it->second;
            if (node.has_coeff()) sum+=node.coeff().real_size();
        }
        world.gop.sum(sum);
        return sum;
    }

    /// Returns the number of coefficients in the function ... collective global sum
    template <typename T, std::size_t NDIM>
    std::size_t FunctionImpl<T,NDIM>::nCoeff() const {
        std::size_t sum = coeffs.size() * (sizeof(keyT) + sizeof(nodeT));
        typename dcT::const_iterator end = coeffs.end();
        for (typename dcT::const_iterator it=coeffs.begin(); it!=end; ++it) {
            const nodeT& node = it->second;
            if (node.has_coeff()) sum+=node.coeff().nCoeff();
        }
        world.gop.sum(sum);
        return sum;
    }


    /// print tree size and size
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::print_size(const std::string name) const {
        const size_t tsize=this->tree_size();
//        const size_t size=this->size();
        const size_t ncoeff=this->nCoeff();
        const double wall=wall_time();
        const double d=sizeof(T);
        const double fac=1024*1024*1024;

        double norm=0.0;
        {
            double local = norm2sq_local();
            this->world.gop.sum(local);
            this->world.gop.fence();
            norm=sqrt(local);
        }

        if (this->world.rank()==0) {

            constexpr std::size_t bufsize=128;
            char buf[bufsize];
            snprintf(buf, bufsize, "%40s at time %.1fs: norm/tree/#coeff/size: %7.5f %zu, %6.3f m, %6.3f GByte",
                   (name.c_str()), wall, norm, tsize,double(ncoeff)*1.e-6,double(ncoeff)/fac*d);
            print(std::string(buf));
        }
    }

    /// print the number of configurations per node
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::print_stats() const {
        if (this->targs.tt==TT_FULL) return;
        int dim=NDIM/2;
        int k0=k;
        if (is_compressed()) k0=2*k;
        Tensor<long> n(int(std::pow(double(k0),double(dim))+1));
        long n_full=0;
        long n_large=0;

        if (world.rank()==0) print("n.size(),k0,dim",n.size(),k0,dim);
        typename dcT::const_iterator end = coeffs.end();
        for (typename dcT::const_iterator it=coeffs.begin(); it!=end; ++it) {
            const nodeT& node = it->second;
            if (node.has_coeff()) {
                if (node.coeff().rank()>long(n.size())) {
                    ++n_large;
                } else if (node.coeff().rank()==-1) {
                    ++n_full;
                } else if (node.coeff().rank()<0) {
                    print("small rank",node.coeff().rank());
                } else {
                    n[node.coeff().rank()]++;
                }
            }
        }

        world.gop.sum(n.ptr(), n.size());

        if (world.rank()==0) {
            print("configurations     number of nodes");
            print("        full rank    ",n_full);
            for (unsigned int i=0; i<n.size(); i++) {
                print("           ",i,"    ",n[i]);
            }
            print("       large rank    ",n_large);

            // repeat for logarithmic scale: <3, <10, <30, <100, ..
            Tensor<long> nlog(6);
            nlog=0;
            for (unsigned int i=0; i<std::min(3l,n.size()); i++) nlog[0]+=n[i];
            for (unsigned int i=3; i<std::min(10l,n.size()); i++) nlog[1]+=n[i];
            for (unsigned int i=10; i<std::min(30l,n.size()); i++) nlog[2]+=n[i];
            for (unsigned int i=30; i<std::min(100l,n.size()); i++) nlog[3]+=n[i];
            for (unsigned int i=100; i<std::min(300l,n.size()); i++) nlog[4]+=n[i];
            for (unsigned int i=300; i<std::min(1000l,n.size()); i++) nlog[5]+=n[i];

            std::vector<std::string> slog={"3","10","30","100","300","1000"};
            for (unsigned int i=0; i<nlog.size(); i++) {
                print("          < ",slog[i],"    ",nlog[i]);
            }
            print("       large rank    ",n_large);

        }
    }

    template <typename T, std::size_t NDIM>
    T FunctionImpl<T,NDIM>::eval_cube(Level n, coordT& x, const tensorT& c) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        const int k = cdata.k;
        double px[NDIM][k];
        T sum = T(0.0);

        for (std::size_t i=0; i<NDIM; ++i) legendre_scaling_functions(x[i],k,px[i]);

        if (NDIM == 1) {
            for (int p=0; p<k; ++p)
                sum += c(p)*px[0][p];
        }
        else if (NDIM == 2) {
            for (int p=0; p<k; ++p)
                for (int q=0; q<k; ++q)
                    sum += c(p,q)*px[0][p]*px[1][q];
        }
        else if (NDIM == 3) {
            for (int p=0; p<k; ++p)
                for (int q=0; q<k; ++q)
                    for (int r=0; r<k; ++r)
                        sum += c(p,q,r)*px[0][p]*px[1][q]*px[2][r];
        }
        else if (NDIM == 4) {
            for (int p=0; p<k; ++p)
                for (int q=0; q<k; ++q)
                    for (int r=0; r<k; ++r)
                        for (int s=0; s<k; ++s)
                            sum += c(p,q,r,s)*px[0][p]*px[1][q]*px[2][r]*px[3][s];
        }
        else if (NDIM == 5) {
            for (int p=0; p<k; ++p)
                for (int q=0; q<k; ++q)
                    for (int r=0; r<k; ++r)
                        for (int s=0; s<k; ++s)
                            for (int t=0; t<k; ++t)
                                sum += c(p,q,r,s,t)*px[0][p]*px[1][q]*px[2][r]*px[3][s]*px[4][t];
        }
        else if (NDIM == 6) {
            for (int p=0; p<k; ++p)
                for (int q=0; q<k; ++q)
                    for (int r=0; r<k; ++r)
                        for (int s=0; s<k; ++s)
                            for (int t=0; t<k; ++t)
                                for (int u=0; u<k; ++u)
                                    sum += c(p,q,r,s,t,u)*px[0][p]*px[1][q]*px[2][r]*px[3][s]*px[4][t]*px[5][u];
        }
        else {
            MADNESS_EXCEPTION("FunctionImpl:eval_cube:NDIM?",NDIM);
        }
        return sum*pow(2.0,0.5*NDIM*n)/sqrt(FunctionDefaults<NDIM>::get_cell_volume());
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::reconstruct_op(const keyT& key, const coeffT& s, const bool accumulate_NS) {
        //PROFILE_MEMBER_FUNC(FunctionImpl);
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
            coeffs.replace(key,nodeT(coeffT(),false));
            it = coeffs.find(key).get();
        }
        nodeT& node = it->second;

        // The integral operator will correctly connect interior nodes
        // to children but may leave interior nodes without coefficients
        // ... but they still need to sum down so just give them zeros
        if (node.has_children() && !node.has_coeff()) {
            node.set_coeff(coeffT(cdata.v2k,targs));
        }

        if (node.has_children() || node.has_coeff()) { // Must allow for inconsistent state from transform, etc.
            coeffT d = node.coeff();
            if (!d.has_data()) d = coeffT(cdata.v2k,targs);
            if (accumulate_NS and (key.level() > 0)) d(cdata.s0) += s; // -- note accumulate for NS summation
            if (d.dim(0)==2*get_k()) {              // d might be pre-truncated if it's a leaf
                d = unfilter(d);
                node.clear_coeff();
                node.set_has_children(true);
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT& child = kit.key();
                    coeffT ss = copy(d(child_patch(child)));
                    ss.reduce_rank(thresh);
                    //PROFILE_BLOCK(recon_send); // Too fine grain for routine profiling
                    woT::task(coeffs.owner(child), &implT::reconstruct_op, child, ss, accumulate_NS);
                }
            } else {
                MADNESS_ASSERT(node.is_leaf());
                //                node.coeff()+=s;
                node.coeff().reduce_rank(targs.thresh);
            }
        }
        else {
            coeffT ss=s;
            if (s.has_no_data()) ss=coeffT(cdata.vk,targs);
            if (key.level()) node.set_coeff(copy(ss));
            else node.set_coeff(ss);
        }
    }

    template <typename T, std::size_t NDIM>
    Tensor<T> fcube(const Key<NDIM>& key, T (*f)(const Vector<double,NDIM>&), const Tensor<double>& qx) {
        //      fcube(key,typename FunctionFactory<T,NDIM>::FunctorInterfaceWrapper(f) , qx, fval);
        std::vector<long> npt(NDIM,qx.dim(0));
        Tensor<T> fval(npt);
        fcube(key,ElementaryInterface<T,NDIM>(f) , qx, fval);
        return fval;
    }

    template <typename T, std::size_t NDIM>
    Tensor<T> fcube(const Key<NDIM>& key, const FunctionFunctorInterface<T,NDIM>& f, const Tensor<double>& qx) {
        //      fcube(key,typename FunctionFactory<T,NDIM>::FunctorInterfaceWrapper(f) , qx, fval);
        std::vector<long> npt(NDIM,qx.dim(0));
        Tensor<T> fval(npt);
        fcube(key, f, qx, fval);
        return fval;
    }

    template <typename T, std::size_t NDIM>
    //    void FunctionImpl<T,NDIM>::fcube(const keyT& key, const FunctionFunctorInterface<T,NDIM>& f, const Tensor<double>& qx, tensorT& fval) const {
    void fcube(const Key<NDIM>& key, const FunctionFunctorInterface<T,NDIM>& f, const Tensor<double>& qx, Tensor<T>& fval) {
        //~ template <typename T, std::size_t NDIM> template< typename FF>
        //~ void FunctionImpl<T,NDIM>::fcube(const keyT& key, const FF& f, const Tensor<double>& qx, tensorT& fval) const {
        typedef Vector<double,NDIM> coordT;
        //PROFILE_MEMBER_FUNC(FunctionImpl);
        const Vector<Translation,NDIM>& l = key.translation();
        const Level n = key.level();
        const double h = std::pow(0.5,double(n));
        coordT c; // will hold the point in user coordinates
        const int npt = qx.dim(0);

        const Tensor<double>& cell_width = FunctionDefaults<NDIM>::get_cell_width();
        const Tensor<double>& cell = FunctionDefaults<NDIM>::get_cell();

        // Do pre-screening of the FunctionFunctorInterface, f, before calculating f(r) at quadrature points
        coordT c1, c2;
        for (std::size_t i = 0; i < NDIM; i++) {
          c1[i] = cell(i,0) + h*cell_width[i]*(l[i] + qx((long)0));
          c2[i] = cell(i,0) + h*cell_width[i]*(l[i] + qx(npt-1));
        }
        if (f.screened(c1, c2)) {
            fval(___) = 0.0;
            return;
        }

        Tensor<double> vqx;
        bool vectorized = f.supports_vectorized();
        if (vectorized) {
            T* fvptr = fval.ptr();
            if (NDIM == 1) {
                double* x1 = new double[npt];
                int idx = 0;
                for (int i=0; i<npt; ++i, ++idx) {
                    c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                    x1[idx] = c[0];
                }
                Vector<double*,1> xvals {x1};
                f(xvals, fvptr, npt);
                delete [] x1;
            }
            else if (NDIM == 2) {
                double* x1 = new double[npt*npt];
                double* x2 = new double[npt*npt];
                int idx = 0;
                for (int i=0; i<npt; ++i) {
                    c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                    for (int j=0; j<npt; ++j, ++idx) {
                        c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                        x1[idx] = c[0];
                        x2[idx] = c[1];
                    }
                }
                Vector<double*,2> xvals {x1, x2};
                f(xvals, fvptr, npt*npt);
                delete [] x1;
                delete [] x2;
            }
            else if (NDIM == 3) {
                double* x1 = new double[npt*npt*npt];
                double* x2 = new double[npt*npt*npt];
                double* x3 = new double[npt*npt*npt];
                int idx = 0;
                for (int i=0; i<npt; ++i) {
                    c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                    for (int j=0; j<npt; ++j) {
                        c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                        for (int k=0; k<npt; ++k, ++idx) {
                            c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                            x1[idx] = c[0];
                            x2[idx] = c[1];
                            x3[idx] = c[2];
                        }
                    }
                }
                Vector<double*,3> xvals {x1, x2, x3};
                f(xvals, fvptr, npt*npt*npt);
                delete [] x1;
                delete [] x2;
                delete [] x3;
            }
            else if (NDIM == 4) {
                double* x1 = new double[npt*npt*npt*npt];
                double* x2 = new double[npt*npt*npt*npt];
                double* x3 = new double[npt*npt*npt*npt];
                double* x4 = new double[npt*npt*npt*npt];
                int idx = 0;
                for (int i=0; i<npt; ++i) {
                    c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                    for (int j=0; j<npt; ++j) {
                        c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                        for (int k=0; k<npt; ++k) {
                            c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                            for (int m=0; m<npt; ++m, ++idx) {
                                c[3] = cell(3,0) + h*cell_width[3]*(l[3] + qx(m)); // xx
                                x1[idx] = c[0];
                                x2[idx] = c[1];
                                x3[idx] = c[2];
                                x4[idx] = c[3];
                            }
                        }
                    }
                }
                Vector<double*,4> xvals {x1, x2, x3, x4};
                f(xvals, fvptr, npt*npt*npt*npt);
                delete [] x1;
                delete [] x2;
                delete [] x3;
                delete [] x4;
            }
            else if (NDIM == 5) {
                double* x1 = new double[npt*npt*npt*npt*npt];
                double* x2 = new double[npt*npt*npt*npt*npt];
                double* x3 = new double[npt*npt*npt*npt*npt];
                double* x4 = new double[npt*npt*npt*npt*npt];
                double* x5 = new double[npt*npt*npt*npt*npt];
                int idx = 0;
                for (int i=0; i<npt; ++i) {
                    c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                    for (int j=0; j<npt; ++j) {
                        c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                        for (int k=0; k<npt; ++k) {
                            c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                            for (int m=0; m<npt; ++m) {
                                c[3] = cell(3,0) + h*cell_width[3]*(l[3] + qx(m)); // xx
                                for (int n=0; n<npt; ++n, ++idx) {
                                    c[4] = cell(4,0) + h*cell_width[4]*(l[4] + qx(n)); // yy
                                    x1[idx] = c[0];
                                    x2[idx] = c[1];
                                    x3[idx] = c[2];
                                    x4[idx] = c[3];
                                    x5[idx] = c[4];
                                }
                            }
                        }
                    }
                }
                Vector<double*,5> xvals {x1, x2, x3, x4, x5};
                f(xvals, fvptr, npt*npt*npt*npt*npt);
                delete [] x1;
                delete [] x2;
                delete [] x3;
                delete [] x4;
                delete [] x5;
            }
            else if (NDIM == 6) {
                double* x1 = new double[npt*npt*npt*npt*npt*npt];
                double* x2 = new double[npt*npt*npt*npt*npt*npt];
                double* x3 = new double[npt*npt*npt*npt*npt*npt];
                double* x4 = new double[npt*npt*npt*npt*npt*npt];
                double* x5 = new double[npt*npt*npt*npt*npt*npt];
                double* x6 = new double[npt*npt*npt*npt*npt*npt];
                int idx = 0;
                for (int i=0; i<npt; ++i) {
                    c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                    for (int j=0; j<npt; ++j) {
                        c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                        for (int k=0; k<npt; ++k) {
                            c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                            for (int m=0; m<npt; ++m) {
                                c[3] = cell(3,0) + h*cell_width[3]*(l[3] + qx(m)); // xx
                                for (int n=0; n<npt; ++n) {
                                    c[4] = cell(4,0) + h*cell_width[4]*(l[4] + qx(n)); // yy
                                    for (int p=0; p<npt; ++p, ++idx) {
                                        c[5] = cell(5,0) + h*cell_width[5]*(l[5] + qx(p)); // zz
                                        x1[idx] = c[0];
                                        x2[idx] = c[1];
                                        x3[idx] = c[2];
                                        x4[idx] = c[3];
                                        x5[idx] = c[4];
                                        x6[idx] = c[5];
                                    }
                                }
                            }
                        }
                    }
                }
                Vector<double*,6> xvals {x1, x2, x3, x4, x5, x6};
                f(xvals, fvptr, npt*npt*npt*npt*npt*npt);
                delete [] x1;
                delete [] x2;
                delete [] x3;
                delete [] x4;
                delete [] x5;
                delete [] x6;
            }
            else {
                MADNESS_EXCEPTION("FunctionImpl: fcube: confused about NDIM?",NDIM);
            }
        }
        else {
            if (NDIM == 1) {
                for (int i=0; i<npt; ++i) {
                    c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                    fval(i) = f(c);
                    MADNESS_ASSERT(!std::isnan(fval(i)));
                }
            }
            else if (NDIM == 2) {
                for (int i=0; i<npt; ++i) {
                    c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                    for (int j=0; j<npt; ++j) {
                        c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                        fval(i,j) = f(c);
                        MADNESS_ASSERT(!std::isnan(fval(i,j)));
                    }
                }
            }
            else if (NDIM == 3) {
                for (int i=0; i<npt; ++i) {
                    c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                    for (int j=0; j<npt; ++j) {
                        c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                        for (int k=0; k<npt; ++k) {
                            c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                            fval(i,j,k) = f(c);
                            MADNESS_ASSERT(!std::isnan(fval(i,j,k)));
                        }
                    }
                }
            }
            else if (NDIM == 4) {
                for (int i=0; i<npt; ++i) {
                    c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                    for (int j=0; j<npt; ++j) {
                        c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                        for (int k=0; k<npt; ++k) {
                            c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                            for (int m=0; m<npt; ++m) {
                                c[3] = cell(3,0) + h*cell_width[3]*(l[3] + qx(m)); // xx
                                fval(i,j,k,m) = f(c);
                                MADNESS_ASSERT(!std::isnan(fval(i,j,k,m)));
                            }
                        }
                    }
                }
            }
            else if (NDIM == 5) {
                for (int i=0; i<npt; ++i) {
                    c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                    for (int j=0; j<npt; ++j) {
                        c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                        for (int k=0; k<npt; ++k) {
                            c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                            for (int m=0; m<npt; ++m) {
                                c[3] = cell(3,0) + h*cell_width[3]*(l[3] + qx(m)); // xx
                                for (int n=0; n<npt; ++n) {
                                    c[4] = cell(4,0) + h*cell_width[4]*(l[4] + qx(n)); // yy
                                    fval(i,j,k,m,n) = f(c);
                                    MADNESS_ASSERT(!std::isnan(fval(i,j,k,m,n)));
                                }
                            }
                        }
                    }
                }
            }
            else if (NDIM == 6) {
                for (int i=0; i<npt; ++i) {
                    c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                    for (int j=0; j<npt; ++j) {
                        c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                        for (int k=0; k<npt; ++k) {
                            c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                            for (int m=0; m<npt; ++m) {
                                c[3] = cell(3,0) + h*cell_width[3]*(l[3] + qx(m)); // xx
                                for (int n=0; n<npt; ++n) {
                                    c[4] = cell(4,0) + h*cell_width[4]*(l[4] + qx(n)); // yy
                                    for (int p=0; p<npt; ++p) {
                                        c[5] = cell(5,0) + h*cell_width[5]*(l[5] + qx(p)); // zz
                                        fval(i,j,k,m,n,p) = f(c);
                                        MADNESS_ASSERT(!std::isnan(fval(i,j,k,m,n,p)));
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
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::fcube(const keyT& key, T (*f)(const coordT&), const Tensor<double>& qx, tensorT& fval) const {
        //      fcube(key,typename FunctionFactory<T,NDIM>::FunctorInterfaceWrapper(f) , qx, fval);
        madness::fcube(key,ElementaryInterface<T,NDIM>(f) , qx, fval);
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::fcube(const keyT& key, const FunctionFunctorInterface<T,NDIM>& f, const Tensor<double>& qx, tensorT& fval) const {
        madness::fcube(key,f,qx,fval);
    }


    /// project the functor into this functionimpl, and "return" a tree in reconstructed,
    /// rank-reduced form.

    /// @param[in]  key current FunctionNode
    /// @param[in]  do_refine
    /// @param[in]  specialpts  in case these are very spiky functions -- don't undersample
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::project_refine_op(const keyT& key,
                                                 bool do_refine,
                                                 const std::vector<Vector<double,NDIM> >& specialpts) {
        //PROFILE_MEMBER_FUNC(FunctionImpl);
        if (do_refine && key.level() < max_refine_level) {

            // Restrict special points to this box
            std::vector<Vector<double,NDIM> > newspecialpts;
            if (key.level() < functor->special_level() && specialpts.size() > 0) {
                BoundaryConditions<NDIM> bc = FunctionDefaults<NDIM>::get_bc();
                const auto bperiodic = bc.is_periodic();
                for (unsigned int i = 0; i < specialpts.size(); ++i) {
                    coordT simpt;
                    user_to_sim(specialpts[i], simpt);
                    Key<NDIM> specialkey = simpt2key(simpt, key.level());
                    if (specialkey.is_neighbor_of(key,bperiodic)) {
                        newspecialpts.push_back(specialpts[i]);
                    }
                }
            }

            // If refining compute scaling function coefficients and
            // norm of difference coefficients
            tensorT r, s0;
            double dnorm = 0.0;
            //////////////////////////if (newspecialpts.size() == 0)
            {
                // Make in r child scaling function coeffs at level n+1
                r = tensorT(cdata.v2k);
                for (KeyChildIterator<NDIM> it(key); it; ++it) {
                    const keyT& child = it.key();
                    r(child_patch(child)) = project(child);
                }
                // Filter then test difference coeffs at level n
                tensorT d = filter(r);
                if (truncate_on_project) s0 = copy(d(cdata.s0));
                d(cdata.s0) = T(0);
                dnorm = d.normf();
            }

            // If have special points always refine.  If don't have special points
            // refine if difference norm is big
            if (newspecialpts.size() > 0 || dnorm >=truncate_tol(thresh,key.level())) {
                coeffs.replace(key,nodeT(coeffT(),true)); // Insert empty node for parent
                for (KeyChildIterator<NDIM> it(key); it; ++it) {
                    const keyT& child = it.key();
                    ProcessID p;
                    if (FunctionDefaults<NDIM>::get_project_randomize()) {
                        p = world.random_proc();
                    }
                    else {
                        p = coeffs.owner(child);
                    }
                    //PROFILE_BLOCK(proj_refine_send); // Too fine grain for routine profiling
                    woT::task(p, &implT::project_refine_op, child, do_refine, newspecialpts);
                }
            }
            else {
                if (truncate_on_project) {
                    coeffT s(s0,thresh,FunctionDefaults<NDIM>::get_tensor_type());
                    coeffs.replace(key,nodeT(s,false));
                }
                else {
                    coeffs.replace(key,nodeT(coeffT(),true)); // Insert empty node for parent
                    for (KeyChildIterator<NDIM> it(key); it; ++it) {
                        const keyT& child = it.key();
                        coeffT s(r(child_patch(child)),thresh,FunctionDefaults<NDIM>::get_tensor_type());
                        coeffs.replace(child,nodeT(s,false));
                    }
                }
            }
        }
        else {
            coeffs.replace(key,nodeT(coeffT(project(key),targs),false));
        }
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::add_scalar_inplace(T t, bool fence) {
        std::vector<long> v0(NDIM,0L);
        std::vector<long> v1(NDIM,1L);
        std::vector<Slice> s(NDIM,Slice(0,0));
        const TensorArgs full_args(-1.0,TT_FULL);
        if (is_compressed()) {
            if (world.rank() == coeffs.owner(cdata.key0)) {
                typename dcT::iterator it = coeffs.find(cdata.key0).get();
                MADNESS_ASSERT(it != coeffs.end());
                nodeT& node = it->second;
                MADNESS_ASSERT(node.has_coeff());
                //                node.node_to_full_rank();
                //                node.full_tensor_reference()(v0) += t*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
                //                node.node_to_low_rank();
                change_tensor_type(node.coeff(),full_args);
                node.coeff().full_tensor()(v0) += t*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
                change_tensor_type(node.coeff(),targs);
            }
        }
        else {
            for (typename dcT::iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
                Level n = it->first.level();
                nodeT& node = it->second;
                if (node.has_coeff()) {
                    // this looks funny, but is necessary for GenTensor, since you can't access a
                    // single matrix element. Therefore make a (1^NDIM) tensor, convert to GenTensor, then
                    // add to the original one by adding a slice.
                    tensorT ttt(v1);
                    ttt=t*sqrt(FunctionDefaults<NDIM>::get_cell_volume()*pow(0.5,double(NDIM*n)));
                    coeffT tt(ttt,get_tensor_args());
                    node.coeff()(s) += tt;
                    // this was the original line:
                    // node.coeff().full_tensor()(v0) += t*sqrt(FunctionDefaults<NDIM>::get_cell_volume()*pow(0.5,double(NDIM*n)));

                }
            }
        }
        if (fence) world.gop.fence();
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::insert_zero_down_to_initial_level(const keyT& key) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        if (is_compressed()) initial_level = std::max(initial_level,1); // Otherwise zero function is confused
        if (coeffs.is_local(key)) {
            if (is_compressed()) {
                if (key.level() == initial_level) {
                    coeffs.replace(key, nodeT(coeffT(), false));
                }
                else {
                    coeffs.replace(key, nodeT(coeffT(cdata.v2k,targs), true));
                }
            }
            else {
                if (key.level()<initial_level) {
                    coeffs.replace(key, nodeT(coeffT(), true));
                }
                else {
                    coeffs.replace(key, nodeT(coeffT(cdata.vk,targs), false));
                }
            }
        }
        if (key.level() < initial_level) {
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                insert_zero_down_to_initial_level(kit.key());
            }
        }

    }


    template <typename T, std::size_t NDIM>
    Future<bool> FunctionImpl<T,NDIM>::truncate_spawn(const keyT& key, double tol) {
        //PROFILE_MEMBER_FUNC(FunctionImpl);
        typename dcT::iterator it = coeffs.find(key).get();
        if (it == coeffs.end()) {
            // In a standard tree all children would exist but some ops (transform)
            // can leave the tree in a messy state.  Just make the missing node as an
            // empty leaf.
            coeffs.replace(key,nodeT());
            it = coeffs.find(key).get();
        }
        nodeT& node = it->second;
        if (node.has_children()) {
            std::vector< Future<bool> > v = future_vector_factory<bool>(1<<NDIM);
            int i=0;
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
                v[i] = woT::task(coeffs.owner(kit.key()), &implT::truncate_spawn, kit.key(), tol, TaskAttributes::generator());
            }
            return woT::task(world.rank(),&implT::truncate_op, key, tol, v);
        }
        else {
            // In compressed form leaves should not have coeffs ... however the
            // transform op could leave the tree with leaves that do have coeffs
            // in which case we want something sensible to happen
            //MADNESS_ASSERT(!node.has_coeff());
            if (node.has_coeff() && key.level()>1) {
                double dnorm = node.coeff().normf();
                if (dnorm < truncate_tol(tol,key)) {
                    node.clear_coeff();
                }
            }
            return Future<bool>(node.has_coeff());
        }
    }


    template <typename T, std::size_t NDIM>
    bool FunctionImpl<T,NDIM>::truncate_op(const keyT& key, double tol, const std::vector< Future<bool> >& v) {
        //PROFILE_MEMBER_FUNC(FunctionImpl); // Too fine grain for routine profiling
        // If any child has coefficients, a parent cannot truncate
        for (int i=0; i<(1<<NDIM); ++i) if (v[i].get()) return true;
        nodeT& node = coeffs.find(key).get()->second;

        // Interior nodes should always have coeffs but transform might
        // leave empty interior nodes ... hence just force no coeffs to
        // be zero coeff unless it is a leaf.
        if (node.has_children() && !node.has_coeff()) node.set_coeff(coeffT(cdata.v2k,targs));

        if (key.level() > 1) { // >1 rather >0 otherwise reconstruct might get confused
            double dnorm = node.coeff().normf();
            if (dnorm < truncate_tol(tol,key)) {
                node.clear_coeff();
                if (node.has_children()) {
                    node.set_has_children(false);
                    for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                        coeffs.erase(kit.key());
                    }
                }
            }
        }
        return node.has_coeff();
    }


    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::print_tree(std::ostream& os, Level maxlevel) const {
        if (world.rank() == 0) do_print_tree(cdata.key0, os, maxlevel);
        world.gop.fence();
        if (world.rank() == 0) os.flush();
        world.gop.fence();
    }


    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::do_print_tree(const keyT& key, std::ostream& os, Level maxlevel) const {
        typename dcT::const_iterator it = coeffs.find(key).get();
        if (it == coeffs.end()) {
            //MADNESS_EXCEPTION("FunctionImpl: do_print_tree: null node pointer",0);
            for (int i=0; i<key.level(); ++i) os << "  ";
            os << key << "  missing --> " << coeffs.owner(key) << "\n";
        }
        else {
            const nodeT& node = it->second;
            for (int i=0; i<key.level(); ++i) os << "  ";
            os << key << "  " << node << " --> " << coeffs.owner(key) << "\n";
            if (key.level() < maxlevel  &&  node.has_children()) {
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    do_print_tree(kit.key(),os,maxlevel);
                }
            }
        }
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::print_tree_json(std::ostream& os, Level maxlevel) const {
        std::multimap<Level, std::tuple<tranT, std::string>> data;
        if (world.rank() == 0) do_print_tree_json(cdata.key0, data, maxlevel);
        world.gop.fence();
        if (world.rank() == 0) {
            for (Level level = 0; level != maxlevel; ++level) {
                if (data.count(level) == 0)
                    break;
                else {
                    if (level > 0)
                        os << ",";
                    os << "\"" << level << "\":{";
                    os << "\"level\": " << level << ",";
                    os << "\"nodes\":{";
                    auto range = data.equal_range(level);
                    for (auto it = range.first; it != range.second; ++it) {
                        os << "\"" << std::get<0>(it->second) << "\":"
                           << std::get<1>(it->second);
                        if (std::next(it) != range.second)
                            os << ",";
                    }
                    os << "}}";
                }
            }
            os.flush();
        }
        world.gop.fence();
    }


    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::do_print_tree_json(const keyT& key, std::multimap<Level, std::tuple<tranT, std::string>>& data, Level maxlevel) const {
        typename dcT::const_iterator it = coeffs.find(key).get();
        if (it == coeffs.end()) {
            MADNESS_EXCEPTION("FunctionImpl: do_print_tree_json: null node pointer",0);
        }
        else {
            const nodeT& node = it->second;
            std::ostringstream oss;
            oss << "{";
            node.print_json(oss);
            oss << ",\"owner\": " << coeffs.owner(key) << "}";
            auto node_json_str = oss.str();
            data.insert(std::make_pair(key.level(), std::make_tuple(key.translation(), node_json_str)));
            if (key.level() < maxlevel  &&  node.has_children()) {
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    do_print_tree_json(kit.key(),data, maxlevel);
                }
            }
        }
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::print_tree_graphviz(std::ostream& os, Level maxlevel) const {
        // aggregate data by level, thus collect data first, then dump
        if (world.rank() == 0) do_print_tree_graphviz(cdata.key0, os, maxlevel);
        world.gop.fence();
        if (world.rank() == 0) os.flush();
        world.gop.fence();
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::do_print_tree_graphviz(const keyT& key, std::ostream& os, Level maxlevel) const {

        struct uniqhash {
            static int64_t value(const keyT& key) {
                int64_t result = 0;
                for (int64_t j = 0; j <= key.level()-1; ++j) {
                    result += (1 << j*NDIM);
                }
                result += key.translation()[0];
                return result;
            }
        };

        typename dcT::const_iterator it = coeffs.find(key).get();
        if (it != coeffs.end()) {
            const nodeT& node = it->second;
            if (key.level() < maxlevel  &&  node.has_children()) {
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    os << uniqhash::value(key) << " -> " << uniqhash::value(kit.key()) << "\n";
                    do_print_tree_graphviz(kit.key(),os,maxlevel);
                }
            }
        }
    }

    template <typename T, std::size_t NDIM>
    Tensor<T> FunctionImpl<T,NDIM>::project(const keyT& key) const {
        //PROFILE_MEMBER_FUNC(FunctionImpl);

        if (not functor) MADNESS_EXCEPTION("FunctionImpl: project: confusion about function?",0);

        // if functor provides coeffs directly, awesome; otherwise use compute by yourself
        if (functor->provides_coeff()) return functor->coeff(key).full_tensor_copy();

        MADNESS_ASSERT(cdata.npt == cdata.k); // only necessary due to use of fast transform
        tensorT fval(cdata.vq,false); // this will be the returned result
        tensorT work(cdata.vk,false); // initially evaluate the function in here
        tensorT workq(cdata.vq,false); // initially evaluate the function in here

        // compute the values of the functor at the quadrature points and scale appropriately
        madness::fcube(key,*functor,cdata.quad_x,work);
        work.scale(sqrt(FunctionDefaults<NDIM>::get_cell_volume()*pow(0.5,double(NDIM*key.level()))));
        //return transform(work,cdata.quad_phiw);
        return fast_transform(work,cdata.quad_phiw,fval,workq);
    }

    template <typename T, std::size_t NDIM>
    Future<double> FunctionImpl<T,NDIM>::get_norm_tree_recursive(const keyT& key) const {
        if (coeffs.probe(key)) {
            return Future<double>(coeffs.find(key).get()->second.get_norm_tree());
        }
        MADNESS_ASSERT(key.level());
        keyT parent = key.parent();
        return woT::task(coeffs.owner(parent), &implT::get_norm_tree_recursive, parent, TaskAttributes::hipri());
    }


    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::sock_it_to_me(const keyT& key,
                                             const RemoteReference< FutureImpl< std::pair<keyT,coeffT> > >& ref) const {
        //PROFILE_MEMBER_FUNC(FunctionImpl);
        if (coeffs.probe(key)) {
            const nodeT& node = coeffs.find(key).get()->second;
            Future< std::pair<keyT,coeffT> > result(ref);
            if (node.has_coeff()) {
                //madness::print("sock found it with coeff",key);
                result.set(std::pair<keyT,coeffT>(key,node.coeff()));
            }
            else {
                //madness::print("sock found it without coeff",key);
                result.set(std::pair<keyT,coeffT>(key,coeffT()));
            }
        }
        else {
            keyT parent = key.parent();
            //madness::print("sock forwarding to parent",key,parent);
            //PROFILE_BLOCK(sitome_send); // Too fine grain for routine profiling
	    if (coeffs.is_local(parent)) 
	      woT::send(coeffs.owner(parent), &FunctionImpl<T,NDIM>::sock_it_to_me, parent, ref);
	    else
	      woT::task(coeffs.owner(parent), &FunctionImpl<T,NDIM>::sock_it_to_me, parent, ref, TaskAttributes::hipri());
        }
    }

    // like sock_it_to_me, but it replaces empty node with averaged coeffs from further down the tree
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::sock_it_to_me_too(const keyT& key,
                                                 const RemoteReference< FutureImpl< std::pair<keyT,coeffT> > >& ref) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        if (coeffs.probe(key)) {
            const nodeT& node = coeffs.find(key).get()->second;
            Future< std::pair<keyT,coeffT> > result(ref);
            if (node.has_coeff()) {
                result.set(std::pair<keyT,coeffT>(key,node.coeff()));
            }
            else {
                result.set(std::pair<keyT,coeffT>(key,nodeT(coeffT(project(key),targs),false).coeff()));
            }
        }
        else {
            keyT parent = key.parent();
            //PROFILE_BLOCK(sitome2_send); // Too fine grain for routine profiling
            woT::task(coeffs.owner(parent), &FunctionImpl<T,NDIM>::sock_it_to_me_too, parent, ref, TaskAttributes::hipri());
        }
    }


    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::eval(const Vector<double,NDIM>& xin,
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
                //PROFILE_BLOCK(eval_send); // Too fine grain for routine profiling
                woT::task(owner, &implT::eval, x, key, ref, TaskAttributes::hipri());
                return;
            }
            else {
                typename dcT::futureT fut = coeffs.find(key);
                typename dcT::iterator it = fut.get();
                nodeT& node = it->second;
                if (node.has_coeff()) {
                    Future<T>(ref).set(eval_cube(key.level(), x, node.coeff().full_tensor_copy()));
                    return;
                }
                else {
                    for (std::size_t i=0; i<NDIM; ++i) {
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


    template <typename T, std::size_t NDIM>
    std::pair<bool,T>
    FunctionImpl<T,NDIM>::eval_local_only(const Vector<double,NDIM>& xin, Level maxlevel) {
        Vector<double,NDIM> x = xin;
        keyT key(0);
        Vector<Translation,NDIM> l = key.translation();
        const ProcessID me = world.rank();
        while (key.level() <= maxlevel) {
            if (coeffs.owner(key) == me) {
                typename dcT::futureT fut = coeffs.find(key);
                typename dcT::iterator it = fut.get();
                if (it != coeffs.end()) {
                    nodeT& node = it->second;
                    if (node.has_coeff()) {
                        return std::pair<bool,T>(true,eval_cube(key.level(), x, node.coeff().full_tensor_copy()));
                    }
                }
            }
            for (std::size_t i=0; i<NDIM; ++i) {
                double xi = x[i]*2.0;
                int li = int(xi);
                if (li == 2) li = 1;
                x[i] = xi - li;
                l[i] = 2*l[i] + li;
            }
            key = keyT(key.level()+1,l);
        }
        return std::pair<bool,T>(false,0.0);
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::evaldepthpt(const Vector<double,NDIM>& xin,
                                           const keyT& keyin,
                                           const typename Future<Level>::remote_refT& ref) {

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
                //PROFILE_BLOCK(eval_send); // Too fine grain for routine profiling
                woT::task(owner, &implT::evaldepthpt, x, key, ref, TaskAttributes::hipri());
                return;
            }
            else {
                typename dcT::futureT fut = coeffs.find(key);
                typename dcT::iterator it = fut.get();
                nodeT& node = it->second;
                if (node.has_coeff()) {
                    Future<Level>(ref).set(key.level());
                    return;
                }
                else {
                    for (std::size_t i=0; i<NDIM; ++i) {
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

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::evalR(const Vector<double,NDIM>& xin,
                                     const keyT& keyin,
                                     const typename Future<long>::remote_refT& ref) {

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
                //PROFILE_BLOCK(eval_send); // Too fine grain for routine profiling
                woT::task(owner, &implT::evalR, x, key, ref, TaskAttributes::hipri());
                return;
            }
            else {
                typename dcT::futureT fut = coeffs.find(key);
                typename dcT::iterator it = fut.get();
                nodeT& node = it->second;
                if (node.has_coeff()) {
                    Future<long>(ref).set(node.coeff().rank());
                    return;
                }
                else {
                    for (std::size_t i=0; i<NDIM; ++i) {
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


    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::tnorm(const tensorT& t, double* lo, double* hi) {
        //PROFILE_MEMBER_FUNC(FunctionImpl); // Too fine grain for routine profiling
        auto& cdata=FunctionCommonData<T,NDIM>::get(t.dim(0));
        tensorT work = copy(t);
        tensorT tlo = work(cdata.sh);
        *lo = tlo.normf();
        tlo.fill(0.0);
        *hi = work.normf();
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::tnorm(const GenTensor<T>& t, double* lo, double* hi) {
        auto& cdata=FunctionCommonData<T,NDIM>::get(t.dim(0));
		coeffT shalf=t(cdata.sh);
		*lo=shalf.normf();
		coeffT sfull=copy(t);
		sfull(cdata.sh)-=shalf;
		*hi=sfull.normf();
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::tnorm(const SVDTensor<T>& t, double* lo, double* hi,
    		const int particle) {
    	*lo=0.0;
    	*hi=0.0;
        auto& cdata=FunctionCommonData<T,NDIM>::get(t.dim(0));
    	if (t.rank()==0) return;
    	const tensorT vec=t.flat_vector(particle-1);
    	for (long i=0; i<t.rank(); ++i) {
    		double lo1,hi1;
    		tensorT c=vec(Slice(i,i),_).reshape(cdata.vk);
    		tnorm(c, &lo1, &hi1);        // note we use g instead of h, since g is 3D
    		*lo+=lo1*t.weights(i);
    		*hi+=hi1*t.weights(i);
		}
    }


    namespace detail {
        template <typename A, typename B>
        struct noop {
            void operator()(const A& a, const B& b) const {};

            template <typename Archive> void serialize(Archive& ar) {}
        };

        template <typename T, std::size_t NDIM>
        struct scaleinplace {
            T q;
            scaleinplace() {}
	    // G++ 4.1.2 ICEs on BGP ... scaleinplace(T q) : q(q) {}
            scaleinplace(T q) {this->q = q;}
            void operator()(const Key<NDIM>& key, Tensor<T>& t) const {
                t.scale(q);
            }
            void operator()(const Key<NDIM>& key, FunctionNode<T,NDIM>& node) const {
                node.coeff().scale(q);
            }
            template <typename Archive> void serialize(Archive& ar) {
                ar & q;
            }
        };

        template <typename T, std::size_t NDIM>
        struct squareinplace {
            void operator()(const Key<NDIM>& key, Tensor<T>& t) const {
                t.emul(t);
            }
            template <typename Archive> void serialize(Archive& ar) {}
        };

        template <typename T, std::size_t NDIM>
        struct absinplace {
            void operator()(const Key<NDIM>& key, Tensor<T>& t) const {t=abs(t);}
            template <typename Archive> void serialize(Archive& ar) {}
        };

        template <typename T, std::size_t NDIM>
        struct abssquareinplace {
            void operator()(const Key<NDIM>& key, Tensor<T>& t) const {abs(t.emul(t));}
            template <typename Archive> void serialize(Archive& ar) {}
        };

    }

template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::scale_inplace(const T q, bool fence) {
        //        unary_op_coeff_inplace(detail::scaleinplace<T,NDIM>(q), fence);
        unary_op_node_inplace(detail::scaleinplace<T,NDIM>(q), fence);
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::square_inplace(bool fence) {
        //unary_op_value_inplace(&implT::autorefine_square_test, detail::squareinplace<T,NDIM>(), fence);
        unary_op_value_inplace(detail::squareinplace<T,NDIM>(), fence);
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::abs_inplace(bool fence) {
        unary_op_value_inplace(detail::absinplace<T,NDIM>(), fence);
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::abs_square_inplace(bool fence) {
        unary_op_value_inplace(detail::abssquareinplace<T,NDIM>(), fence);
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::phi_for_mul(Level np, Translation lp, Level nc, Translation lc, Tensor<double>& phi) const {
        //PROFILE_MEMBER_FUNC(FunctionImpl); // Too fine grain for routine profiling
        double p[200];
        double scale = pow(2.0,double(np-nc));
        for (int mu=0; mu<cdata.npt; ++mu) {
            double xmu = scale*(cdata.quad_x(mu)+lc) - lp;
            MADNESS_ASSERT(xmu>-1e-15 && xmu<(1+1e-15));
            legendre_scaling_functions(xmu,cdata.k,p);
            for (int i=0; i<k; ++i) phi(i,mu) = p[i];
        }
        phi.scale(pow(2.0,0.5*np));
    }

    template <typename T, std::size_t NDIM>

    const GenTensor<T> FunctionImpl<T,NDIM>::parent_to_child(const coeffT& s, const keyT& parent, const keyT& child) const {
        //PROFILE_MEMBER_FUNC(FunctionImpl); // Too fine grain for routine profiling
        // An invalid parent/child means that they are out of the box
        // and it is the responsibility of the caller to worry about that
        // ... most likely the coefficients (s) are zero to reflect
        // zero B.C. so returning s makes handling this easy.
        if (parent == child || parent.is_invalid() || child.is_invalid()) return s;

        coeffT result = fcube_for_mul<T>(child, parent, s);
        result.scale(sqrt(FunctionDefaults<NDIM>::get_cell_volume()*pow(0.5,double(NDIM*child.level()))));
        result = transform(result,cdata.quad_phiw);

        return result;
    }


    template <typename T, std::size_t NDIM>
    T FunctionImpl<T,NDIM>::trace_local() const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        std::vector<long> v0(NDIM,0);
        T sum = 0.0;
        if (is_compressed()) {
            if (world.rank() == coeffs.owner(cdata.key0)) {
                typename dcT::const_iterator it = coeffs.find(cdata.key0).get();
                if (it != coeffs.end()) {
                    const nodeT& node = it->second;
                    if (node.has_coeff()) sum = node.coeff().full_tensor_copy()(v0);
                }
            }
        }
        else {
            for (typename dcT::const_iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
                const keyT& key = it->first;
                const nodeT& node = it->second;
                if (node.has_coeff()) sum += node.coeff().full_tensor_copy()(v0)*pow(0.5,NDIM*key.level()*0.5);
            }
        }
        return sum*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
    }


    static inline bool enforce_bc(bool is_periodic, Level n, Translation& l) {
      const Translation two2n = 1ul << n;
      if (l < 0) {
        if (is_periodic) {
          do {
            l += two2n; // Periodic BC
          } while (l < 0);
        } else
          return false; // Zero BC
      } else if (l >= two2n) {
        if (is_periodic) {
          do {
            l -= two2n; // Periodic BC
          } while (l >= two2n);
        } else
          return false; // Zero BC
      }
      return true;
    }

    static inline bool enforce_in_volume(Level n, Translation& l) {
      Translation two2n = 1ul << n;
      return l >= 0 && l < two2n;
    }

    template <typename T, std::size_t NDIM>
    Key<NDIM> FunctionImpl<T,NDIM>::neighbor(const keyT& key, const Key<NDIM>& disp, const std::array<bool, NDIM>& is_periodic) const {
        Vector<Translation,NDIM> l = key.translation();

        for (std::size_t axis=0; axis<NDIM; ++axis) {
            l[axis] += disp.translation()[axis];

            //if (!enforce_bc(bc(axis,0), bc(axis,1), key.level(), l[axis])) {
            if (!enforce_bc(is_periodic[axis], key.level(), l[axis])) {
                return keyT::invalid();
            }
        }
        return keyT(key.level(),l);
    }

    template <typename T, std::size_t NDIM>
    Key<NDIM> FunctionImpl<T,NDIM>::neighbor_in_volume(const keyT& key, const Key<NDIM>& disp) const {
      Vector<Translation, NDIM> l = key.translation();

      for (std::size_t axis = 0; axis < NDIM; ++axis) {
        l[axis] += disp.translation()[axis];

        if (!enforce_in_volume(key.level(), l[axis])) {
          return keyT::invalid();
        }
      }
      return keyT(key.level(), l);
    }

    template <typename T, std::size_t NDIM>
    Future< std::pair< Key<NDIM>, GenTensor<T> > >
    FunctionImpl<T,NDIM>::find_me(const Key<NDIM>& key) const {
        //PROFILE_MEMBER_FUNC(FunctionImpl); // Too fine grain for routine profiling
        typedef std::pair< Key<NDIM>,coeffT > argT;
        Future<argT> result;
        //PROFILE_BLOCK(find_me_send); // Too fine grain for routine profiling
        woT::task(coeffs.owner(key), &implT::sock_it_to_me_too, key, result.remote_ref(world), TaskAttributes::hipri());
        return result;
    }


    /// will insert
    /// @return s coefficient and norm_tree for key
    template <typename T, std::size_t NDIM>
    Future< std::pair<GenTensor<T>,double> > FunctionImpl<T,NDIM>::compress_spawn(const Key<NDIM>& key,
				bool nonstandard1, bool keepleaves, bool redundant1) {
        if (!coeffs.probe(key)) print("missing node",key);
        MADNESS_ASSERT(coeffs.probe(key));

        // get fetches remote data (here actually local)
        nodeT& node = coeffs.find(key).get()->second;

        // internal node -> continue recursion
        if (node.has_children()) {
            std::vector< Future<std::pair<coeffT,double> > > v = future_vector_factory<std::pair<coeffT,double> >(1<<NDIM);
            int i=0;
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
                //PROFILE_BLOCK(compress_send); // Too fine grain for routine profiling
                // readily available
                v[i] = woT::task(coeffs.owner(kit.key()), &implT::compress_spawn, kit.key(),
                                 nonstandard1, keepleaves, redundant1, TaskAttributes::hipri());
            }
            if (redundant1) return woT::task(world.rank(),&implT::make_redundant_op, key, v);
            return woT::task(world.rank(),&implT::compress_op, key, v, nonstandard1);
        }

        // leaf node -> remove coefficients here and pass them back to parent for filtering
        // insert snorm, dnorm=0.0, normtree (=snorm)
        else {
            // special case: tree has only root node: keep sum coeffs and make zero diff coeffs
            if (key.level()==0) {
                coeffT result(node.coeff());
                coeffT sdcoeff(cdata.v2k,this->get_tensor_type());
                sdcoeff(cdata.s0)+=node.coeff();
                node.coeff()=sdcoeff;
                double snorm=node.coeff().normf();
                node.set_dnorm(0.0);
                node.set_snorm(snorm);
                node.set_norm_tree(snorm);
                return Future< std::pair<GenTensor<T>,double> >(std::make_pair(result,node.coeff().normf()));

            } else { // this is a leaf node
                Future<coeffT > result(node.coeff());
                if (not keepleaves) node.clear_coeff();

                auto snorm=(keepleaves) ? node.coeff().normf() : 0.0;
                node.set_norm_tree(snorm);
                node.set_snorm(snorm);
                node.set_dnorm(0.0);

                return Future< std::pair<GenTensor<T>,double> >(std::make_pair(result,snorm));
            }
        }
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::plot_cube_kernel(archive::archive_ptr< Tensor<T> > ptr,
                                                const keyT& key,
                                                const coordT& plotlo, const coordT& plothi, const std::vector<long>& npt,
                                                bool eval_refine) const {

        Tensor<T>& r = *ptr;

        coordT h; // Increment between points in each dimension
        for (std::size_t i=0; i<NDIM; ++i) {
            if (npt[i] > 1) {
                h[i] = (plothi[i]-plotlo[i])/(npt[i]-1);
            }
            else {
                MADNESS_ASSERT(plotlo[i] == plothi[i]);
                h[i] = 0.0;
            }
        }

        const Level n = key.level();
        const Vector<Translation,NDIM>& l = key.translation();
        const double twon = pow(2.0,double(n));
        const tensorT& coeff = coeffs.find(key).get()->second.coeff().full_tensor_copy(); // Ugh!
        //        const tensorT coeff = coeffs.find(key).get()->second.full_tensor_copy(); // Ugh!
        long ind[NDIM];
        coordT x;

        coordT boxlo, boxhi;
        Vector<int,NDIM> boxnpt;
        double fac = pow(0.5,double(key.level()));
        int npttotal = 1;
        for (std::size_t d=0; d<NDIM; ++d) {
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
                boxlo[d] = std::max(boxlo[d],plotlo[d]);
                boxhi[d] = std::min(boxhi[d],plothi[d]);

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
            for (IndexIterator it(boxnpt); it; ++it) {
                for (std::size_t d=0; d<NDIM; ++d) {
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
                if (eval_refine) {
                    r(ind) = n;
                }
                else {
                    T tmp = eval_cube(n, x, coeff);
                    r(ind) = tmp;
                    //print("    eval", ind, tmp, r(ind));
                }
            }
        }
    }

    /// Set plot_refine=true to get a plot of the refinment levels of
    /// the given function (defaulted to false in prototype).
    template <typename T, std::size_t NDIM>
    Tensor<T> FunctionImpl<T,NDIM>::eval_plot_cube(const coordT& plotlo,
                                                   const coordT& plothi,
                                                   const std::vector<long>& npt,
                                                   const bool eval_refine) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        Tensor<T> r(NDIM, &npt[0]);
        //r(___) = 99.0;
        MADNESS_ASSERT(is_reconstructed());

        for (typename dcT::const_iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
            const keyT& key = it->first;
            const nodeT& node = it->second;
            if (node.has_coeff()) {
                woT::task(world.rank(), &implT::plot_cube_kernel,
                          archive::archive_ptr< Tensor<T> >(&r), key, plotlo, plothi, npt, eval_refine);
            }
        }

        //        ITERATOR(r, if (r(IND) == 99.0) {print("BAD", IND); error("bad",0);});

        world.taskq.fence();
        world.gop.sum(r.ptr(), r.size());
        world.gop.fence();

        return r;
    }

    static inline void dxprintvalue(FILE* f, const double t) {
        fprintf(f,"%.6e\n",t);
    }

    static inline void dxprintvalue(FILE* f, const double_complex& t) {
        fprintf(f,"%.6e %.6e\n", t.real(), t.imag());
    }

    template <typename T, std::size_t NDIM>
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
            for (std::size_t d=0; d<NDIM; ++d) fprintf(f," %ld",npt[d]);
            fprintf(f,"\n");

            fprintf(f,"origin ");
            for (std::size_t d=0; d<NDIM; ++d) fprintf(f, " %.6e", cell(d,0));
            fprintf(f,"\n");

            for (std::size_t d=0; d<NDIM; ++d) {
                fprintf(f,"delta ");
                for (std::size_t c=0; c<d; ++c) fprintf(f, " 0");
                double h = 0.0;
                if (npt[d]>1) h = (cell(d,1)-cell(d,0))/(npt[d]-1);
                fprintf(f," %.6e", h);
                for (std::size_t c=d+1; c<NDIM; ++c) fprintf(f, " 0");
                fprintf(f,"\n");
            }
            fprintf(f,"\n");

            fprintf(f,"object 2 class gridconnections counts ");
            for (std::size_t d=0; d<NDIM; ++d) fprintf(f," %ld",npt[d]);
            fprintf(f,"\n");
            fprintf(f, "attribute \"element type\" string \"%s\"\n", element[NDIM-1]);
            fprintf(f, "attribute \"ref\" string \"positions\"\n");
            fprintf(f,"\n");

            int npoint = 1;
            for (std::size_t d=0; d<NDIM; ++d) npoint *= npt[d];
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
                fwrite((void *) r.ptr(), sizeof(T), r.size(), f);
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

    template <std::size_t NDIM>
    void FunctionDefaults<NDIM>::set_defaults(World& world) {
        k = 6;
        thresh = 1e-4;
        initial_level = 2;
        special_level = 3;
        max_refine_level = 30;
        truncate_mode = 0;
        refine = true;
        autorefine = true;
        debug = false;
        truncate_on_project = true;
        apply_randomize = false;
        project_randomize = false;
        bc = BoundaryConditions<NDIM>(BC_FREE);
        tt = TT_FULL;
        cell = make_default_cell();
        recompute_cell_info();
        set_default_pmap(world);
    }

    template <std::size_t NDIM>
    void FunctionDefaults<NDIM>::set_default_pmap(World& world) {
        //pmap = std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >(new WorldDCDefaultPmap< Key<NDIM> >(world));
        pmap = std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >(new madness::LevelPmap< Key<NDIM> >(world));
        //pmap = std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >(new SimplePmap< Key<NDIM> >(world));
    }


    template <std::size_t NDIM>
    void FunctionDefaults<NDIM>::print(){
    		std::cout << "Function Defaults:" << std::endl;
    		std::cout << "                      Dimension " <<  ": " <<  NDIM << std::endl;
    		std::cout << "                               k" <<  ": " << k << std::endl;
    		std::cout << "                          thresh" <<  ": " << thresh << std::endl;
    		std::cout << "                   initial_level" <<  ": " << initial_level << std::endl;
    		std::cout << "                   special_level" <<  ": " << special_level << std::endl;
    		std::cout << "                max_refine_level" <<  ": " << max_refine_level << std::endl;
    		std::cout << "                   truncate_mode" <<  ": " << truncate_mode << std::endl;
    		std::cout << "                          refine" <<  ": " << refine << std::endl;
    		std::cout << "                      autorefine" <<  ": " << autorefine << std::endl;
    		std::cout << "                           debug" <<  ": " << debug << std::endl;
    		std::cout << "             truncate_on_project" <<  ": " << truncate_on_project << std::endl;
    		std::cout << "                 apply_randomize" <<  ": " << apply_randomize << std::endl;
    		std::cout << "               project_randomize" <<  ": " << project_randomize << std::endl;
    		std::cout << "                              bc" <<  ": " << bc << std::endl;
    		std::cout << "                              tt" <<  ": " << tt << std::endl;
    		std::cout << "                            cell" <<  ": " << cell << std::endl;
    }

    template <typename T, std::size_t NDIM>
    const FunctionCommonData<T,NDIM>* FunctionCommonData<T,NDIM>::data[MAXK] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    // default values match those in FunctionDefaults::set_defaults(world)
    template <std::size_t NDIM> int FunctionDefaults<NDIM>::k = 6;
    template <std::size_t NDIM> double FunctionDefaults<NDIM>::thresh = 1e-4;
    template <std::size_t NDIM> int FunctionDefaults<NDIM>::initial_level = 2;
    template <std::size_t NDIM> int FunctionDefaults<NDIM>::special_level = 3;
    template <std::size_t NDIM> int FunctionDefaults<NDIM>::max_refine_level = 30;
    template <std::size_t NDIM> int FunctionDefaults<NDIM>::truncate_mode = 0;
    template <std::size_t NDIM> bool FunctionDefaults<NDIM>::refine = true;
    template <std::size_t NDIM> bool FunctionDefaults<NDIM>::autorefine = true;
    template <std::size_t NDIM> bool FunctionDefaults<NDIM>::debug = false;
    template <std::size_t NDIM> bool FunctionDefaults<NDIM>::truncate_on_project = true;
    template <std::size_t NDIM> bool FunctionDefaults<NDIM>::apply_randomize = false;
    template <std::size_t NDIM> bool FunctionDefaults<NDIM>::project_randomize = false;
    template <std::size_t NDIM> BoundaryConditions<NDIM> FunctionDefaults<NDIM>::bc = BoundaryConditions<NDIM>(BC_FREE);
    template <std::size_t NDIM> TensorType FunctionDefaults<NDIM>::tt = TT_FULL;
    template <std::size_t NDIM> Tensor<double> FunctionDefaults<NDIM>::cell = FunctionDefaults<NDIM>::make_default_cell();
    template <std::size_t NDIM> Tensor<double> FunctionDefaults<NDIM>::cell_width = FunctionDefaults<NDIM>::make_default_cell_width();
    template <std::size_t NDIM> Tensor<double> FunctionDefaults<NDIM>::rcell_width = FunctionDefaults<NDIM>::make_default_cell_width();
    template <std::size_t NDIM> double FunctionDefaults<NDIM>::cell_volume = 1.;
    template <std::size_t NDIM> double FunctionDefaults<NDIM>::cell_min_width = 1.;
    template <std::size_t NDIM> std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > > FunctionDefaults<NDIM>::pmap;

    template <std::size_t NDIM> std::vector< Key<NDIM> > Displacements<NDIM>::disp;
    template <std::size_t NDIM> std::vector< Key<NDIM> > Displacements<NDIM>::disp_periodic[64];
    template <std::size_t NDIM> array_of_bools<NDIM> Displacements<NDIM>::disp_periodic_axes{false};

}

#endif // MADNESS_MRA_MRAIMPL_H__INCLUDED
