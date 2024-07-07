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

#ifndef MADNESS_DERIVATIVE_H__INCLUDED
#define MADNESS_DERIVATIVE_H__INCLUDED

#include <iostream>
#include <iomanip>
#include <fstream>
#include <madness/world/MADworld.h>
#include <madness/world/worlddc.h>
#include <madness/world/print.h>
#include <madness/misc/misc.h>

#include <madness/tensor/tensor.h>
#include <madness/tensor/gentensor.h>

#include <madness/mra/key.h>
#include <madness/mra/funcdefaults.h>


/// \file mra/derivative.h
/// \brief Declaration and initialization of tree traversal functions and generic derivative
/// \ingroup mra

namespace madness {

	template<typename T, std::size_t NDIM>
	class FunctionNode;

    template<typename T, std::size_t NDIM>
    class Function;

}



namespace madness {


/// Tri-diagonal operator traversing tree primarily for derivative operator

    /// \ingroup mra
    template <typename T, std::size_t NDIM>
    class DerivativeBase : public WorldObject< DerivativeBase<T, NDIM> > {
        typedef WorldObject< DerivativeBase<T, NDIM> > woT;
    protected:
        World& world;
        const std::size_t axis      ;  ///< Axis along which the operation is performed
        const int k         ;  ///< Number of wavelets of the function
        const BoundaryConditions<NDIM> bc;
        const std::vector<long> vk; ///< (k,...) used to initialize Tensors

    public:
        friend class FunctionImpl<T, NDIM>;

        typedef Tensor<T>               tensorT  ;	///< regular tensors, like rm, etc
        typedef GenTensor<T>            coeffT   ;	///< holding the node's coeffs (possibly low rank)
        typedef Key<NDIM>               keyT     ;
        typedef std::pair<keyT,coeffT>  argT     ;
        typedef FunctionImpl<T,NDIM>    implT    ;
        typedef Function<T,NDIM>        functionT;
        typedef WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > dcT;
        typedef FunctionNode<T,NDIM> nodeT;


        DerivativeBase(World& world, std::size_t axis, int k, BoundaryConditions<NDIM> bc)
            : WorldObject< DerivativeBase<T, NDIM> >(world)
            , world(world)
            , axis(axis)
            , k(k)
            , bc(bc)
            , vk(NDIM,k)
        {
            // No!  Cannot process incoming messages until the *derived* class is constructed.
            // this->process_pending();
        }

        virtual ~DerivativeBase() { }

        void forward_do_diff1(const implT* f, implT* df, const keyT& key,
                              const argT& left,
                              const argT& center,
                              const argT& right)  const {

            const dcT& coeffs = f->get_coeffs();
            ProcessID owner = coeffs.owner(key);

            if (owner == world.rank()) {
                if (!left.second.has_data()) {
                    woT::task(owner, &madness::DerivativeBase<T,NDIM>::do_diff1,
                            f, df, key, find_neighbor(f, key,-1), center, right,
                            TaskAttributes::hipri());
                }
                else if (!right.second.has_data()) {
                    woT::task(owner, &madness::DerivativeBase<T,NDIM>::do_diff1,
                            f, df, key, left, center, find_neighbor(f, key,1),
                            TaskAttributes::hipri());
                }
                // Boundary node
                else if (left.first.is_invalid() || right.first.is_invalid()) {
                    woT::task(owner, &madness::DerivativeBase<T,NDIM>::do_diff2b,
                            f, df, key, left, center, right);
                }
                // Interior node
                else {
                    woT::task(owner, &madness::DerivativeBase<T,NDIM>::do_diff2i,
                            f, df, key, left, center, right);
                }
            }
            else {
                df->task(owner, &madness::FunctionImpl<T,NDIM>::forward_do_diff1,
                        this, f, key, left, center, right, TaskAttributes::hipri());
            }
        }

        void do_diff1(const implT* f, implT* df, const keyT& key,
                      const argT& left,
                      const argT& center,
                      const argT& right) const {
            MADNESS_ASSERT(axis<NDIM);

//            if (left.second.size()==0 || right.second.size()==0) {
            if ((!left.second.has_data()) || (!right.second.has_data())) {
                // One of the neighbors is below us in the tree ... recur down
                df->get_coeffs().replace(key,nodeT(coeffT(),true));
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT& child = kit.key();
                    if ((child.translation()[axis]&1) == 0) {
                        // leftmost child automatically has right sibling
                        forward_do_diff1(f, df, child, left, center, center);
                    }
                    else {
                        // rightmost child automatically has left sibling
                        forward_do_diff1(f, df, child, center, center, right);
                    }
                }
            }
            else {
                forward_do_diff1(f, df, key, left, center, right);
            }
        }

        virtual void do_diff2b(const implT* f, implT* df, const keyT& key,
                               const argT& left,
                               const argT& center,
                               const argT& right) const = 0;

        virtual void do_diff2i(const implT* f, implT* df, const keyT& key,
                               const argT& left,
                               const argT& center,
                               const argT& right) const = 0;


        /// Differentiate w.r.t. given coordinate (x=0, y=1, ...) with optional fence

        /// Returns a new function with the same distribution
        Function<T,NDIM>
        operator()(const functionT& f, bool fence=true) const {
            if (VERIFY_TREE) f.verify_tree();
            if (fence) f.change_tree_state(reconstructed);
            MADNESS_CHECK_THROW(f.is_reconstructed(),"diff: trying to diff a compressed function without fencing");

            functionT df;
            df.set_impl(f,false);

            df.get_impl()->diff(this, f.get_impl().get(), fence);
            return df;
        }


        static bool enforce_bc(int bc_left, int bc_right, Level n, Translation& l) {
            Translation two2n = 1ul << n;
            if (l < 0) {
                if (bc_left == BC_ZERO || bc_left == BC_FREE || bc_left == BC_DIRICHLET || bc_left == BC_ZERONEUMANN || bc_left == BC_NEUMANN) {
                    return false; // f=0 BC, or no BC, or nonzero f BC, or zero deriv BC, or nonzero deriv BC
                }
                else if (bc_left == BC_PERIODIC) {
                    l += two2n; // Periodic BC
                    MADNESS_ASSERT(bc_left == bc_right);   //check that both BCs are periodic
                }
                else {
                    MADNESS_EXCEPTION("enforce_bc: confused left BC?",bc_left);
                }
            }
            else if (l >= two2n) {
                if (bc_right == BC_ZERO || bc_right == BC_FREE || bc_right == BC_DIRICHLET || bc_right == BC_ZERONEUMANN || bc_right == BC_NEUMANN) {
                    return false; // f=0 BC, or no BC, or nonzero f BC, or zero deriv BC, or nonzero deriv BC
                }
                else if (bc_right == BC_PERIODIC) {
                    l -= two2n; // Periodic BC
                    MADNESS_ASSERT(bc_left == bc_right);   //check that both BCs are periodic
                }
                else {
                    MADNESS_EXCEPTION("enforce_bc: confused BC right?",bc_right);
                }
            }
            return true;
        }

        Key<NDIM> neighbor(const keyT& key, int step) const {
            Vector<Translation,NDIM> l = key.translation();
            l[axis] += step;
            if (!enforce_bc(bc(axis,0), bc(axis,1), key.level(), l[axis])) {
                return keyT::invalid();
            }
            else {
                return keyT(key.level(),l);
            }
        }

        Future<argT>
        find_neighbor(const implT* f, const Key<NDIM>& key, int step) const {
            keyT neigh = neighbor(key, step);
            if (neigh.is_invalid()) {
                return Future<argT>(argT(neigh,coeffT(vk,f->get_tensor_args()))); // Zero bc
            }
            else {
                Future<argT> result;
		if (f->get_coeffs().is_local(neigh))
		  f->send(f->get_coeffs().owner(neigh), &implT::sock_it_to_me, neigh, result.remote_ref(world));
		else
		  f->task(f->get_coeffs().owner(neigh), &implT::sock_it_to_me, neigh, result.remote_ref(world), TaskAttributes::hipri());
                return result;
            }
        }


        template <typename Archive> void serialize(const Archive& ar) const {
            throw "NOT IMPLEMENTED";
        }

    };  // End of the DerivativeBase class


    /// Implements derivatives operators with variety of boundary conditions on simulation domain
    template <typename T, std::size_t NDIM>
    class Derivative : public DerivativeBase<T, NDIM> {
    private:
        typedef DerivativeBase<T, NDIM> baseT;

    public:
        typedef Tensor<T>               tensorT  ;
        typedef GenTensor<T>            coeffT   ;	///< holding the node's coeffs (possibly low rank)
        typedef Key<NDIM>               keyT     ;
        typedef std::pair<keyT,coeffT>  argT     ;
        typedef FunctionImpl<T,NDIM>    implT    ;
        typedef Function<T,NDIM>        functionT;
        typedef WorldContainer< Key<NDIM> , FunctionNode<T, NDIM> > dcT;
        typedef FunctionNode<T,NDIM> nodeT;

    private:
        const functionT g1;  ///< Function describing the boundary condition on the right side
        const functionT g2;  ///< Function describing the boundary condition on the left side

        bool is_second;
        bool is_third;

        // Tensors for holding the modified coefficients
        Tensor<double> rm, r0, rp        ; ///< Blocks of the derivative operator
        Tensor<double> rmt, r0t, rpt     ; ///< Blocks of the derivative operator, transposed
        Tensor<double> left_rm, left_r0  ; ///< Blocks of the derivative for the left boundary
        Tensor<double> left_rmt, left_r0t  ; ///< Blocks of the derivative for the left boundary
        Tensor<double> right_r0, right_rp; ///< Blocks of the derivative for the right boundary
        Tensor<double> right_r0t, right_rpt; ///< Blocks of the derivative for the right boundary
        Tensor<double> bv_left, bv_right ; ///< Blocks of the derivative operator for the boundary contribution


        // Tensors for the bspline smoothed central difference operator
        Tensor<double> r0_bsp;
        Tensor<double> rm_bsp;
        Tensor<double> rp_bsp;
        Tensor<double> r0_bsp_t;
        Tensor<double> rm_bsp_t;
        Tensor<double> rp_bsp_t;

        void do_diff2b(const implT* f, implT* df, const keyT& key,
                       const argT& left,
                       const argT& center,
                       const argT& right) const {
            Vector<Translation,NDIM> l = key.translation();
            double lev   = (double) key.level();

            coeffT d;

            //left boundary
            if (l[this->axis] == 0) {

                coeffT tensor_right=df->parent_to_child(right.second, right.first, this->neighbor(key,1));
                coeffT tensor_center=df->parent_to_child(center.second, center.first, key);

                d= transform_dir(tensor_right,left_rmt,this->axis);
                d+=transform_dir(tensor_center,left_r0t,this->axis);
            }
            else {

                coeffT tensor_left=df->parent_to_child(left.second, left.first, this->neighbor(key,-1));
                coeffT tensor_center=df->parent_to_child(center.second, center.first, key);

                d= transform_dir(tensor_left,right_rpt,this->axis);
                d+=transform_dir(tensor_center,right_r0t,this->axis);
            }

            double fac = FunctionDefaults<NDIM>::get_rcell_width()[this->axis]*pow(2.0,lev);
            if (is_second) fac *= fac;
            else if (is_third) fac *= fac*fac;

            d.scale(fac);
            d.reduce_rank(df->get_thresh());
            df->get_coeffs().replace(key,nodeT(d,false));


            // This is the boundary contribution (formally in BoundaryDerivative)
            int bc_left  = this->bc(this->axis,0);
            int bc_right = this->bc(this->axis,1);

            Future<argT> found_argT;
            tensorT bf, bdry_t;
            //left boundary
            if (l[this->axis] == 0) {
                if (bc_left != BC_PERIODIC && bc_left != BC_FREE && bc_left != BC_ZERO && bc_left != BC_ZERONEUMANN) {
                    bf = copy(bv_left);
                    found_argT = g1.get_impl()->find_me(key);
                }
                else {
                    return;
                }
            }
            else { //right boundary
                if (bc_right != BC_PERIODIC && bc_right != BC_FREE && bc_right != BC_ZERO && bc_right != BC_ZERONEUMANN) {
                    bf = copy(bv_right);
                    found_argT = g2.get_impl()->find_me(key);
                }
                else {
                    return;
                }
            }
#ifdef HAVE_PARSEC
            std::cerr << "FATAL ERROR: PaRSEC does not support recursive task execution but Derivative::do_diff2b requires this. Use a different backend" << std::endl;
            abort();
#endif
            const auto& found_argT_value = found_argT.get();  // do not recursively execute tasks to avoid making PaRSEC sad
            tensorT gcoeffs = df->parent_to_child(found_argT_value.second, found_argT_value.first,key).full_tensor_copy();

            //if (this->bc.get_bc().dim(0) == 1) {
            if (NDIM == 1) {
                bdry_t = gcoeffs[0]*bf;
            }
            else {
                tensorT slice_aid(this->k);  //vector of zeros
                slice_aid[0] = 1;
                tensorT tmp = inner(slice_aid, gcoeffs, 0, this->axis);
                bdry_t = outer(bf,tmp);
                if (this->axis) bdry_t = copy(bdry_t.cycledim(this->axis,0,this->axis)); // make it contiguous
            }
            bdry_t.scale(FunctionDefaults<NDIM>::get_rcell_width()[this->axis]);

            if (l[this->axis]==0) {
                if (bc_left == BC_DIRICHLET)
                    bdry_t.scale( pow(2.0,lev));
				else if (bc_left ==BC_NEUMANN)
					bdry_t.scale(FunctionDefaults<NDIM>::get_cell_width()[this->axis]);
            }
            else {
                if (bc_right == BC_DIRICHLET)
                    bdry_t.scale( pow(2.0,lev));
				else if (bc_right ==BC_NEUMANN)
					bdry_t.scale(FunctionDefaults<NDIM>::get_cell_width()[this->axis]);
            }

            bdry_t += d.full_tensor_copy();;
            df->get_coeffs().replace(key,nodeT(coeffT(bdry_t,df->get_thresh(),df->get_tensor_type()),false));
        }

        void do_diff2i(const implT* f, implT*df, const keyT& key,
                       const argT& left,
                       const argT& center,
                       const argT& right) const
        {
//#if !HAVE_GENTENSOR
//            coeffT d = madness::inner(rp,
//                                       df->parent_to_child(left.second, left.first, baseT::neighbor(key,-1)).swapdim(this->axis,0),
//                                       1, 0);
//            inner_result(r0,
//                         df->parent_to_child(center.second, center.first, key).swapdim(this->axis,0),
//                         1, 0, d);
//            inner_result(rm,
//                         df->parent_to_child(right.second, right.first, baseT::neighbor(key,1)).swapdim(this->axis,0),
//                         1, 0, d);
//            // flo thinks this is wrong for higher dimensions -- need to cycledim
//            if (this->axis) d = copy(d.swapdim(this->axis,0)); // make it contiguous
//            d.scale(FunctionDefaults<NDIM>::get_rcell_width()[this->axis]*pow(2.0,(double) key.level()));
//            df->get_coeffs().replace(key,nodeT(d,false));
//
//#else
            coeffT tensor_left=df->parent_to_child(left.second, left.first, this->neighbor(key,-1));
            coeffT tensor_center=df->parent_to_child(center.second, center.first, key);
            coeffT tensor_right=df->parent_to_child(right.second, right.first, this->neighbor(key,1));

            coeffT d= transform_dir(tensor_left,rpt,this->axis);
            d+=transform_dir(tensor_center,r0t,this->axis);
            d+=transform_dir(tensor_right,rmt,this->axis);

            double fac = FunctionDefaults<NDIM>::get_rcell_width()[this->axis]*pow(2.0,(double) key.level());
            if (is_second) fac *= fac;
            else if (is_third) fac *= fac*fac;

            d.scale(fac);
            d.reduce_rank(df->get_thresh());
            df->get_coeffs().replace(key,nodeT(d,false));

//#endif

        }

        void initCoefficients()  {
            is_second = false;
            is_third = false;

            r0 = Tensor<double>(this->k,this->k);
            rp = Tensor<double>(this->k,this->k);
            rm = Tensor<double>(this->k,this->k);

            left_rm = Tensor<double>(this->k,this->k);
            left_r0 = Tensor<double>(this->k,this->k);

            right_r0 = Tensor<double>(this->k,this->k);
            right_rp = Tensor<double>(this->k,this->k);

            // These are the coefficients for the boundary contribution
            bv_left  = Tensor<double>(this->k);
            bv_right = Tensor<double>(this->k);

            int bc_left  = this->bc(this->axis,0);
            int bc_right = this->bc(this->axis,1);

            double kphase = -1.0;
            if (this->k%2 == 0) kphase = 1.0;
            double iphase = 1.0;
            for (int i=0; i<this->k; ++i) {
                double jphase = 1.0;
                for (int j=0; j<this->k; ++j) {
                    double gammaij = sqrt(double((2*i+1)*(2*j+1)));
                    double Kij;
                    if (((i-j)>0) && (((i-j)%2)==1))
                        Kij = 2.0;
                    else
                        Kij = 0.0;

                    r0(i,j) = 0.5*(1.0 - iphase*jphase - 2.0*Kij)*gammaij;
                    rm(i,j) = 0.5*jphase*gammaij;
                    rp(i,j) =-0.5*iphase*gammaij;

                    // Constraints on the derivative
                    if (bc_left == BC_ZERONEUMANN || bc_left == BC_NEUMANN) {
                        left_rm(i,j) = jphase*gammaij*0.5*(1.0 + iphase*kphase/this->k);

                        double phi_tmpj_left = 0;

                        for (int l=0; l<this->k; ++l) {
                            double gammalj = sqrt(double((2*l+1)*(2*j+1)));
                            double Klj;

                            if (((l-j)>0) && (((l-j)%2)==1))  Klj = 2.0;
                            else   Klj = 0.0;

                            phi_tmpj_left += sqrt(double(2*l+1))*Klj*gammalj;
                        }
                        phi_tmpj_left = -jphase*phi_tmpj_left;
                        left_r0(i,j) = (0.5*(1.0 + iphase*kphase/this->k) - Kij)*gammaij + iphase*sqrt(double(2*i+1))*phi_tmpj_left/pow(this->k,2.);
                    }
                    else if (bc_left == BC_ZERO || bc_left == BC_DIRICHLET || bc_left == BC_FREE) {
                        left_rm(i,j) = rm(i,j);

                        // B.C. with a function
                        if (bc_left == BC_ZERO || bc_left == BC_DIRICHLET)
                            left_r0(i,j) = (0.5 - Kij)*gammaij;

                        // No B.C.
                        else if (bc_left == BC_FREE)
                            left_r0(i,j) = (0.5 - iphase*jphase - Kij)*gammaij;
                    }

                    // Constraints on the derivative
                    if (bc_right == BC_ZERONEUMANN || bc_right == BC_NEUMANN) {
                        right_rp(i,j) = -0.5*(iphase + kphase / this->k)*gammaij;

                        double phi_tmpj_right = 0;
                        for (int l=0; l<this->k; ++l) {
                            double gammalj = sqrt(double((2*l+1)*(2*j+1)));
                            double Klj;
                            if (((l-j)>0) && (((l-j)%2)==1))  Klj = 2.0;
                            else   Klj = 0.0;
                            phi_tmpj_right += sqrt(double(2*l+1))*Klj*gammalj;
                        }
                        right_r0(i,j) = -(0.5*jphase*(iphase+ kphase/this->k) + Kij)*gammaij + sqrt(double(2*i+1))*phi_tmpj_right/pow(this->k,2.);
                    }
                    else if (bc_right == BC_ZERO || bc_right == BC_FREE || bc_right == BC_DIRICHLET) {
                        right_rp(i,j) = rp(i,j);

                        // Zero BC
                        if (bc_right == BC_ZERO || bc_right == BC_DIRICHLET)
                            right_r0(i,j) = -(0.5*iphase*jphase + Kij)*gammaij;

                        // No BC
                        else if (bc_right == BC_FREE)
                            right_r0(i,j) = (1.0 - 0.5*iphase*jphase - Kij)*gammaij;

                    }

                    jphase = -jphase;
                }
                iphase = -iphase;
            }

            // Coefficients for the boundary contributions
            iphase = 1.0;
            for (int i=0; i<this->k; ++i) {
                iphase = -iphase;

                if (bc_left == BC_DIRICHLET)
                    bv_left(i) = iphase*sqrt(double(2*i+1));            // vector for left dirichlet BC
                else if(bc_left == BC_NEUMANN)
                    bv_left(i) = -iphase*sqrt(double(2*i+1))/pow(this->k,2.);  // vector for left deriv BC
                else
                    bv_left(i) = 0.0;

                if (bc_right == BC_DIRICHLET)
                    bv_right(i) = sqrt(double(2*i+1));                  // vector for right dirichlet BC
                else if (bc_right == BC_NEUMANN)
                    bv_right(i) = sqrt(double(2*i+1))/pow(this->k,2.);         // vector for right deriv BC
                else
                    bv_right(i) = 0.0;
            }

            r0t = transpose(r0);
            rpt = transpose(rp);
            rmt = transpose(rm);

            right_r0t = transpose(right_r0);
            right_rpt = transpose(right_rp);

            left_rmt = transpose(left_rm);
            left_r0t = transpose(left_r0);

            //print(rm.normf(),r0.normf(),rp.normf(),left_rm.normf(),left_r0.normf(),right_r0.normf(),right_rp.normf(),bv_left.normf(),bv_right.normf());
        }

    public:
        typedef T opT;

        /// Constructs a derivative operator

        /// @param world The world
        /// @param axis The direction to differentiate
        /// @param bc Boundary conditions (default from FunctionDefaults)
        /// @param g1 Function providing left boundary value (default empty)
        /// @param g2 Function providing right boundary value (default empty)
        /// @param k Wavelet order (default from FunctionDefaults)
        Derivative(World& world,
                   std::size_t axis,
                   const BoundaryConditions<NDIM>& bc=FunctionDefaults<NDIM>::get_bc(),
                   const functionT g1=functionT(),
                   const functionT g2=functionT(),
                   int k=FunctionDefaults<NDIM>::get_k())
            :  DerivativeBase<T, NDIM>(world, axis, k, bc)
            , g1(g1)
            , g2(g2)
        {
            MADNESS_ASSERT(axis<NDIM);
            initCoefficients();
            g1.reconstruct();
            g2.reconstruct();

            this->process_pending();
        }

        virtual ~Derivative() { }

        void set_is_first() {is_second = false; is_third = false;}
        void set_is_second() {is_second = true; is_third=false;}
        void set_is_third() {is_second = false; is_third = true;}

        void set_bspline1() {
           int k = FunctionDefaults<NDIM>::get_k();
           if(k > 18) throw "Bspline derivatives are only available up to k=18";
           std::string filename = get_mra_data_dir() + "/b-spline-deriv1.txt";
           read_from_file(filename, 1);
        }

        void set_bspline2() {
           int k = FunctionDefaults<NDIM>::get_k();
           if(k > 18) throw "Bspline derivatives are only available up to k=18";
           std::string filename = get_mra_data_dir() + "/b-spline-deriv2.txt";
           read_from_file(filename, 2);
        }

        void set_bspline3() {
           int k = FunctionDefaults<NDIM>::get_k();
           if(k > 18) throw "Bspline derivatives are only available up to k=18";
           std::string filename = get_mra_data_dir() + "/b-spline-deriv3.txt";
           read_from_file(filename, 3);
        }

        void set_ble1() {
           int k = FunctionDefaults<NDIM>::get_k();
           if(k > 15) throw "BLE derivatives are only available up to k=15";
           std::string filename = get_mra_data_dir() + "/ble-first.txt";
           read_from_file(filename, 1);
        }

        void set_ble2() {
           int k = FunctionDefaults<NDIM>::get_k();
           if(k > 15) throw "BLE derivatives are only available up to k=15";
           std::string filename = get_mra_data_dir() + "/ble-second.txt";
           read_from_file(filename, 2);
        }

        void read_from_file(const std::string& filename, unsigned int order = 1) {

            Tensor<double> r0_bsp(this->k,this->k);
            Tensor<double> rp_bsp(this->k,this->k);
            Tensor<double> rm_bsp(this->k,this->k);

            std::ifstream f(filename);
            bool found=false;

            for (int m; f >> m; ) {
                if (m == this->k) {
                    for (int i=0; i<m; i++)
                        for (int j=0; j<m; j++)
                            MADNESS_CHECK(f >> rp_bsp(i,j));
                    for (int i=0; i<m; i++)
                        for (int j=0; j<m; j++)
                            MADNESS_CHECK(f >> r0_bsp(i,j));
                    for (int i=0; i<m; i++)
                        for (int j=0; j<m; j++)
                            MADNESS_CHECK(f >> rm_bsp(i,j));
                    found = true;
                    break;
                }
                else {
                    double junk;
                    for (int i=0; i<3*m*m; i++)
                        MADNESS_CHECK(f >> junk);
                }
            }
            MADNESS_CHECK(found);
            Tensor<double> r0_bsp_t = transpose(r0_bsp);
            Tensor<double> rp_bsp_t = transpose(rp_bsp);
            Tensor<double> rm_bsp_t = transpose(rm_bsp);

            r0=r0_bsp; r0t=r0_bsp_t; left_r0=r0_bsp; left_r0t=r0_bsp_t; right_r0=r0_bsp; right_r0t=r0_bsp_t;

            rp=rp_bsp; rpt=rp_bsp_t; right_rp=rp_bsp; right_rpt=rp_bsp_t;

            rm=rm_bsp; rmt=rm_bsp_t; left_rm=rm_bsp; left_rmt=rm_bsp_t;

            // Get scaling factor right for higher order derivatives
            if (order == 1) {
               set_is_first();
            }
            else if(order == 2) {
               set_is_second();
            }
            else if(order == 3) {
               set_is_third();
            }
        }
    };


    /// Convenience function returning derivative operator with free-space boundary conditions
    template <typename T, std::size_t NDIM>
    Derivative<T,NDIM>
    free_space_derivative(World& world, int axis, int k=FunctionDefaults<NDIM>::get_k()) {
        return Derivative<T, NDIM>(world, axis, BoundaryConditions<NDIM>(BC_FREE), Function<T,NDIM>(), Function<T,NDIM>(), k);
    }


    /// Conveinence function returning derivative operator with periodic boundary conditions
    template <typename T, std::size_t NDIM>
    Derivative<T,NDIM>
    periodic_derivative(World& world, int axis, int k=FunctionDefaults<NDIM>::get_k()) {
        return Derivative<T, NDIM>(world, axis, BoundaryConditions<NDIM>(BC_PERIODIC), Function<T,NDIM>(), Function<T,NDIM>(), k);
    }

    /// Applies derivative operator to function (for syntactic equivalence to integral operator apply)
    template <typename T, std::size_t NDIM>
    Function<T,NDIM>
    apply(const Derivative<T,NDIM>& D, const Function<T,NDIM>& f, bool fence=true) {
        return D(f,fence);
    }

    /// Convenience function returning vector of derivative operators implementing grad (\f$ \nabla \f$)

    /// This will only work for BC_ZERO, BC_PERIODIC, BC_FREE and
    /// BC_ZERONEUMANN since we are not passing in any boundary
    /// functions.
    template <typename T, std::size_t NDIM>
    std::vector< std::shared_ptr< Derivative<T,NDIM> > >
    gradient_operator(World& world,
                      const BoundaryConditions<NDIM>& bc = FunctionDefaults<NDIM>::get_bc(),
                      int k = FunctionDefaults<NDIM>::get_k()) {
        std::vector< std::shared_ptr< Derivative<T,NDIM> > > r(NDIM);
        for (std::size_t d=0; d<NDIM; ++d) {
            MADNESS_CHECK(bc(d,0)!=BC_DIRICHLET && bc(d,1)!=BC_DIRICHLET);
            MADNESS_CHECK(bc(d,0)!=BC_NEUMANN   && bc(d,1)!=BC_NEUMANN);
            r[d].reset(new Derivative<T,NDIM>(world,d,bc,Function<T,NDIM>(),Function<T,NDIM>(),k));
        }
        return r;
    }


    namespace archive {
        template <class Archive, class T, std::size_t NDIM>
        struct ArchiveLoadImpl<Archive,const DerivativeBase<T,NDIM>*> {
            static void load(const Archive& ar, const DerivativeBase<T,NDIM>*& ptr) {
                WorldObject< DerivativeBase<T,NDIM> >* p = nullptr;
                ar & p;
                ptr = static_cast< const DerivativeBase<T,NDIM>* >(p);
            }
        };

        template <class Archive, class T, std::size_t NDIM>
        struct ArchiveStoreImpl<Archive,const DerivativeBase<T,NDIM>*> {
            static void store(const Archive& ar, const DerivativeBase<T,NDIM>* const & ptr) {
                ar & ptr->id();
            }
        };
    }

}  // End of the madness namespace

#endif // MADNESS_MRA_DERIVATIVE_H_INCLUDED
