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


  $Id: mraimpl.h 1602 2009-12-27 19:53:06Z rjharrison $
*/

#ifndef MADNESS_DERIVATIVE_H__INCLUDED
#define MADNESS_DERIVATIVE_H__INCLUDED

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
//#include <mra/mra.h>
//#include <world/worldhashmap.h>

#include <iostream>
#include <world/world.h>
#include <world/print.h>
#include <misc/misc.h>
#include <tensor/tensor.h>
#include <mra/key.h>
#include <mra/funcdefaults.h>
#include <mra/funcimpl.h>
#include <mra/loadbal.h>

//extern "C" double round(double x);


/// \file mra/derivative.h
/// \brief Declaration and initialization of tree traversal functions and generic
///        derivative

namespace madness {

/**
 * BoundaryConds is a class for boundary conditions. It contains, Tensor<int> bc, which holds the 
 * Boundary Condition flags, and methods for checking that bc has the correct dimensions, and that
 * the integer codes in bc are supported and consistent.
 *
 * For periodic codes the is_valid_bc method checks that both boundaries have the periodic flag.
 *
 */


	template<int NDIM>
	class BoundaryConds { 

		private: 

			Tensor<int> bc; 

			bool is_valid_bc_code(const int code) {
				bool is_valid = true;
				if(!(code==0||code==1||code==2||code==3||code==4||code==5)) is_valid = false;
				return is_valid;
			}

			bool is_valid_bc(){
				bool is_valid = true;
            	if(!(bc.dim(0)==NDIM && bc.dim(1)==2 && bc.ndim()==2)) is_valid = false;
	        	else for(int i=0; i<NDIM ; ++i){
    	        	int tmp = bc(i,0);
        	    	if(!is_valid_bc_code(tmp)) is_valid = false;
	            	tmp = bc(i,1);
    	        	if(!is_valid_bc_code(tmp)) is_valid = false;
    	    	}
	        	return is_valid;
			}

			bool is_valid_bc(const Tensor<int> potential_bc){
				bool is_valid = true;
            	if(!(potential_bc.dim(0)==NDIM && potential_bc.dim(1)==2 && potential_bc.ndim()==2)) is_valid = false;
	        	else for(int i=0; i<NDIM ; ++i){
    	        	int tmp = potential_bc(i,0);
        	    	if(!(tmp==0||tmp==1||tmp==2||tmp==3||tmp==4||tmp==5)) is_valid = false;
	            	tmp = potential_bc(i,1);
    	        	if(!(tmp==0||tmp==1||tmp==2||tmp==3||tmp==4||tmp==5)) is_valid = false;
					if(tmp == 1){ 
						if (!(potential_bc(i,0) == 1)) is_valid = false; // both sides should be periodic
					}
    	    	}
	        	return is_valid;
			}

		public:

			// Constructor. Default boundary condition set to periodic
			
		    BoundaryConds(int code=1) {//=1?
				if (!is_valid_bc_code(code)){
		            MADNESS_EXCEPTION("BoundaryConds: invalid boundary condition",0);   /// no idea what the EXCEPTION code 0 represents
				} 
				else{
					long dim[2];
					dim[0]=NDIM;
					dim[1]=2;
			    	Tensor<int> tmp = Tensor<int>(2,dim);
					tmp = code;
					bc = tmp;
				}
		    };	

			// General constructor. 
			
		    BoundaryConds(Tensor<int> bc_tensor) {
				if(is_valid_bc(bc_tensor)){
					bc = copy(bc_tensor);   /// need to change this to a deep copy
					//bc = Tensor<int>(bc_tensor);   /// need to change this to a deep copy
					//bc = bc_tensor;   /// need to change this to a deep copy
				}
				else {
		            MADNESS_EXCEPTION("BoundaryConds: invalid boundary condition",0);   /// no idea what the EXCEPTION code 0 represents
				}
		    };	

			// get_bc for operators to extract the bc tensor
 		    const Tensor<int>& get_bc() const{
				return bc;
		    }


	        double operator()(int i, int j) const {
                     return bc(i,j) ;
	        }

			// assignment operator
			inline BoundaryConds<NDIM>&
	        operator=(const Tensor<int>& other) {
        		if(is_valid_bc(other)){
                	bc = copy(other);
			    }
			    else{
        	        MADNESS_EXCEPTION("operator= : invalid boundary condition",0);   /// no idea what the EXCEPTION code 0 represents
            	}
    	        return *this;
	        }

            template <typename Archive> void serialize(const Archive& ar) {};
		
	};






  template <typename T, int NDIM>
  class TreeTraversal : public WorldObject< TreeTraversal<T, NDIM> > {

    public:
      friend class FunctionImpl<T, NDIM> ;

      typedef Tensor<T>               tensorT   ;
      typedef Key<NDIM>               keyT      ;
      typedef std::pair<keyT,tensorT> argT      ;
      typedef FunctionImpl<T,NDIM>    implT     ;
      typedef Function<T,NDIM>        functionT ;
      typedef WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > dcT ;
      typedef FunctionNode<T,NDIM> nodeT ;

      World& world ; 

      BoundaryConds<NDIM> bc ;

      functionT    f ;  // Initial function
      functionT   df ;  // Derivative function

      int k          ;  // Number of wavelets of the function
      int axis       ;  // Axis along which the operation is performed

      std::vector<long> vk; ///< (k,...) used to initialize Tensors


      void impldiff(bool fence) {
          PROFILE_MEMBER_FUNC(TreeTraversal) ;
          typedef std::pair<keyT,tensorT> argT;
          dcT coeffs = f.get_impl()->get_coeffs() ;

          for (typename dcT::const_iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
              const keyT& key = it->first;
              const nodeT& node = it->second;
              if (node.has_coeff()) {
                  Future<argT> left = find_neighbor(key,-1);
                  argT center(key,node.coeff());
                  Future<argT> right  = find_neighbor(key, 1);
                  task(world.rank(), &madness::TreeTraversal<T, NDIM>::do_diff1, key, left, center, right, TaskAttributes::hipri());
              } 
              else {
                  // Internal empty node can be safely inserted
                  df.get_impl()->replace_coeff(key,nodeT(tensorT(),true)) ;
              }
          }
          if (fence) world.gop.fence();
      }

      Void forward_do_diff1(const keyT& key,
              const std::pair<keyT,tensorT>& left,
              const std::pair<keyT,tensorT>& center,
              const std::pair<keyT,tensorT>& right) {
          PROFILE_MEMBER_FUNC(TreeTraversal);

          dcT coeffs = f.get_impl()->get_coeffs() ;
          ProcessID owner = coeffs.owner(key);

          if (owner == world.rank()) {
              if (left.second.size() == 0) {
                  task(owner, &madness::TreeTraversal<T,NDIM>::do_diff1, key, find_neighbor(key,-1), center, right, TaskAttributes::hipri());
              }
              else if (right.second.size() == 0) {
                  task(owner, &madness::TreeTraversal<T,NDIM>::do_diff1, key, left, center, find_neighbor(key,1), TaskAttributes::hipri());
              }
              // Boundary node
              else if (left.first.is_invalid() || right.first.is_invalid()) { 
                  task(owner, &madness::TreeTraversal<T,NDIM>::do_diff2b, key, left, center, right);
              }
              // Interior node
              else { 
                  task(owner, &madness::TreeTraversal<T,NDIM>::do_diff2i, key, left, center, right);
              }
          }
          else {
              task(owner, &madness::TreeTraversal<T,NDIM>::forward_do_diff1, key, left, center, right, TaskAttributes::hipri());
          }
          return None;
      }

      Void do_diff1(const keyT& key,
                    const std::pair<keyT,tensorT>& left,
                    const std::pair<keyT,tensorT>& center,
                    const std::pair<keyT,tensorT>& right) {
          PROFILE_MEMBER_FUNC(TreeTraversal);
  
          MADNESS_ASSERT(axis>=0 && axis<NDIM);
  
          if (left.second.size()==0 || right.second.size()==0) {
              // One of the neighbors is below us in the tree ... recur down
              this->df.get_impl()->replace_coeff(key,nodeT(tensorT(),true));
              for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                  const keyT& child = kit.key();
                  if ((child.translation()[axis]&1) == 0) {
                      // leftmost child automatically has right sibling
                      forward_do_diff1(child, left, center, center);
                  }
                  else {
                      // rightmost child automatically has left sibling
                      forward_do_diff1(child, center, center, right);
                  }
              }
          }
          else {
            forward_do_diff1(key, left, center, right);
          }
          return None;
      }

      virtual Void do_diff2b(const keyT& key,
                             const std::pair<keyT,tensorT>& left,
                             const std::pair<keyT,tensorT>& center,
                             const std::pair<keyT,tensorT>& right) {return None;} ;

      virtual Void do_diff2i(const keyT& key,
                             const std::pair<keyT,tensorT>& left,
                             const std::pair<keyT,tensorT>& center,
                             const std::pair<keyT,tensorT>& right) {return None;} ;


      /// Differentiate w.r.t. given coordinate (x=0, y=1, ...) with optional fence

      /// Returns a new function with the same distribution
      Function<T,NDIM>
      operator()(functionT _f, bool fence=true) {
          PROFILE_FUNC;

          f = _f ;

          if (f.is_compressed()) {
              if (fence) {
                  f.reconstruct();
              }
              else {
                  MADNESS_EXCEPTION("diff: trying to diff a compressed function without fencing",0);
              }
          }
  
          if (VERIFY_TREE) f.verify_tree();

          df.set_impl(f,false) ;

          impldiff( fence) ;
          //impldiff(*f.get_impl(), *df.get_impl(), fence) ;
          //impldiff(&(f.get_impl().get()), &(df.get_impl().get()), fence) ;
          return df;
      }


      /// Take the derivative of a vector of functions
      /// THIS IS UNTESTED!!!

      /// Operates on a vector of functions
      std::vector< functionT> 
      operator()(std::vector<functionT> _vf, bool fence=true) {
          PROFILE_FUNC;

          std::vector<functionT> dvf(_vf.size() ) ; 

          for (unsigned int i=0; i<_vf.size(); i++)
          {
            dvf[i]= (*this)(_vf[i], false) ;
            if (((i+1) % VMRA_CHUNK_SIZE) == 0) world.gop.fence() ;
          }
          if (fence) world.gop.fence() ;
          return dvf ;
      }

      static bool enforce_bc(int bc_left, int bc_right, Level n, Translation& l) {
          Translation two2n = 1ul << n;
          if (l < 0) {
              if (bc_left == 0 || bc_left == 2 || bc_left ==3 || bc_left == 4 || bc_left == 5) {
                  return false; // f=0 BC, or no BC, or nonzero f BC, or zero deriv BC, or nonzero deriv BC
              }
              else if (bc_left == 1) {
                  l += two2n; // Periodic BC
                  MADNESS_ASSERT(bc_left == bc_right);   //check that both BCs are periodic
              }
              else {
                  MADNESS_EXCEPTION("enforce_bc: confused left BC?",bc_left);
              }
          }
          else if (l >= two2n) {
              if (bc_right == 0 || bc_right == 2 || bc_right == 3 || bc_right ==4 || bc_right == 5) {
                  return false; // f=0 BC, or no BC, or nonzero f BC, or zero deriv BC, or nonzero deriv BC
              }
              else if (bc_right == 1) {
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

      Future< std::pair< Key<NDIM>,Tensor<T> > >
      find_neighbor(const Key<NDIM>& key, int step) const {
          PROFILE_MEMBER_FUNC(TreeTraversal);
          keyT neigh = neighbor(key, step);
          if (neigh.is_invalid()) {
              return Future<argT>(argT(neigh,tensorT(vk))); // Zero bc
          }
          else {
              Future<argT> result;
              PROFILE_BLOCK(find_neigh_send);
              f.get_impl()->task(f.get_impl()->get_coeffs().owner(neigh), &implT::sock_it_to_me, neigh, result.remote_ref(world), TaskAttributes::hipri());
              return result;
          }
      }


      TreeTraversal(World& _world, int _k, int _axis, BoundaryConds<NDIM> _bc) :
         WorldObject< TreeTraversal<T, NDIM> >(_world) ,
         world(_world) ,
         bc(_bc) ,
         k(_k) ,
         axis(_axis) 
      {
            vk = std::vector<long>(NDIM);
 
            for (int i = 0; i < NDIM; i++) {
                vk[i] = k;
            }
      } 

       template <typename Archive> void serialize(const Archive& ar) {
           throw "NOT IMPLEMENTED";
       }

  };  // End of the TreeTraversal class


  template <typename T, int NDIM>
  class Derivative : public TreeTraversal<T, NDIM> {
    public:
      typedef Tensor<T>               tensorT   ;
      typedef Key<NDIM>               keyT      ;
      typedef std::pair<keyT,tensorT> argT      ;
      typedef FunctionImpl<T,NDIM>    implT     ;
      typedef Function<T,NDIM>        functionT ;
      typedef WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > dcT ;
      typedef FunctionNode<T,NDIM> nodeT ;

      functionT g1      ;  // Function describing the boundary condition on the right side
      functionT g2      ;  // Function describing the boundary condition on the left side

      // Tensors for holding the modified coefficients
      Tensor<double> rm, r0, rp         ; ///< Blocks of the derivative operator
      Tensor<double> left_rm, left_r0   ; ///< Blocks of the derivative for the left boundary
      Tensor<double> right_r0, right_rp ; ///< Blocks of the derivative for the right boundary
      Tensor<double> bv_left, bv_right  ; ///< Blocks of the derivative operator for the boundary contribution

      Void do_diff2b(const keyT& key,
                     const std::pair<keyT,tensorT>& left,
                     const std::pair<keyT,tensorT>& center,
                     const std::pair<keyT,tensorT>& right) {
          PROFILE_MEMBER_FUNC(Derivative);
          Vector<Translation,NDIM> l = key.translation();
          double lev   = (double) key.level() ;

          tensorT d ;
 
          //left boundary
          if (l[this->axis] == 0) {
             d = madness::inner(left_rm ,
                                this->df.get_impl()->parent_to_child(right.second, right.first, neighbor(key,1)).swapdim(this->axis,0),
                                1, 0);
             inner_result(left_r0,
                          this->df.get_impl()->parent_to_child(center.second, center.first, key).swapdim(this->axis,0),
                          1, 0, d);
          }
          else
          {
             d = madness::inner(right_rp,
                                this->df.get_impl()->parent_to_child(left.second, left.first, neighbor(key,-1)).swapdim(this->axis,0),
                                1, 0);
             inner_result(right_r0,
                          this->df.get_impl()->parent_to_child(center.second, center.first, key).swapdim(this->axis,0),
                          1, 0, d);
          }
          if (this->axis) d = copy(d.swapdim(this->axis,0)); // make it contiguous
          d.scale(FunctionDefaults<NDIM>::get_rcell_width()[this->axis]*pow(2.0,lev));
          this->df.get_impl()->replace_coeff(key,nodeT(d,false));


          // This is the boundary contribution (formally in BoundaryDerivative)
          //int k = this->f.get_impl()->get_k() ; 

          int bc_left  = this->bc(this->axis,0);
          int bc_right = this->bc(this->axis,1);
  
          Future<argT> found_argT ;
          tensorT bf, bdry_t ;
          //left boundary
          if (l[this->axis] == 0) 
          {
            if (bc_left != 1 && bc_left != 2)
            {
              bf = copy(bv_left) ;
              found_argT = g1.get_impl()->find_me(key) ;
            }
            else
              return None ;
          }
          //right boundary
          else 
          {
            if (bc_right != 1 && bc_right != 2)
            {
              bf = copy(bv_right) ;
              found_argT = g2.get_impl()->find_me(key) ;
            }
            else
              return None ;
          }

          tensorT gcoeffs = found_argT.get().second;

          if (this->bc.get_bc().dim(0) == 1)           
            bdry_t = gcoeffs[0]*bf;
          else
          {
            tensorT slice_aid(this->k);  //vector of zeros
            slice_aid[0] = 1;
            tensorT tmp = inner(slice_aid, gcoeffs, 0, this->axis);
            bdry_t = outer(bf,tmp);
            if (this->axis) bdry_t = copy(bdry_t.cycledim(this->axis,0,this->axis)); // make it contiguous
          }
          bdry_t.scale(FunctionDefaults<NDIM>::get_rcell_width()[this->axis]);

          if (l[this->axis]==0)
          {
             if (bc_left == 3)
               bdry_t.scale( pow(2.0,lev)) ;
          }
          else
          {
             if (bc_right == 3)
              bdry_t.scale( pow(2.0,lev)) ;
          }

          bdry_t = bdry_t + d ;

          this->df.get_impl()->replace_coeff(key,nodeT(bdry_t,false));
 
          return None;
      }

      Void do_diff2i(const keyT& key,
                  const std::pair<keyT,tensorT>& left,
                  const std::pair<keyT,tensorT>& center,
                  const std::pair<keyT,tensorT>& right) {
          PROFILE_MEMBER_FUNC(Derivative);
          tensorT d = madness::inner(rp,
                                     this->df.get_impl()->parent_to_child(left.second, left.first, neighbor(key,-1)).swapdim(this->axis,0),
                                     1, 0);
          inner_result(r0,
                       this->df.get_impl()->parent_to_child(center.second, center.first, key).swapdim(this->axis,0),
                       1, 0, d);
          inner_result(rm,
                       this->df.get_impl()->parent_to_child(right.second, right.first, neighbor(key,1)).swapdim(this->axis,0),
                       1, 0, d);
          if (this->axis) d = copy(d.swapdim(this->axis,0)); // make it contiguous
          d.scale(FunctionDefaults<NDIM>::get_rcell_width()[this->axis]*pow(2.0,(double) key.level()));
          this->df.get_impl()->replace_coeff(key,nodeT(d,false));
          return None;
      }


      //Derivative(World& _world, functionT _f, int _axis, BoundaryConds<NDIM>& _bc, functionT _g1, functionT _g2) :
      Derivative(World& _world, int _k, int _axis, const BoundaryConds<NDIM>& _bc, functionT _g1, functionT _g2) :
         TreeTraversal<T, NDIM>(_world, _k, _axis, _bc) ,
         g1(_g1) ,
         g2(_g2) 
      {
         r0 = Tensor<double>(this->k,this->k);
         rp = Tensor<double>(this->k,this->k);
         rm = Tensor<double>(this->k,this->k);

          left_rm = Tensor<double>(this->k,this->k);
          left_r0 = Tensor<double>(this->k,this->k);

         right_r0 = Tensor<double>(this->k,this->k);
         right_rp = Tensor<double>(this->k,this->k);

         int bc_left  = this->bc(this->axis,0);
         int bc_right = this->bc(this->axis,1);

         double kphase = -1.0;
         if (this->k%2 == 0) kphase = 1.0;
         double iphase = 1.0;
         for (int i=0; i<this->k; i++) {
             double jphase = 1.0;
             for (int j=0; j<this->k; j++) {
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
                 if (bc_left == 4 || bc_left == 5)
                 {
                   left_rm(i,j) = jphase*gammaij*0.5*(1.0 + iphase*kphase/this->k); 

                   double phi_tmpj_left = 0;

                   for (int l=0; l<this->k; l++) {
                     double gammalj = sqrt(double((2*l+1)*(2*j+1)));
                     double Klj;

                     if (((l-j)>0) && (((l-j)%2)==1))  Klj = 2.0;
                     else   Klj = 0.0;

                     phi_tmpj_left += sqrt(double(2*l+1))*Klj*gammalj;
                   }
                   phi_tmpj_left = -jphase*phi_tmpj_left;
                   left_r0(i,j) = (0.5*(1.0 + iphase*kphase/this->k) - Kij)*gammaij + iphase*sqrt(double(2*i+1))*phi_tmpj_left/pow(this->k,2.) ;
                 }
                 else if (bc_left == 0 || bc_left == 3 || bc_left == 2)
                 {
                   left_rm(i,j) = rm(i,j) ;
                 
                   // B.C. with a function
                   if (bc_left == 0 || bc_left == 3)
                     left_r0(i,j) = (0.5 - Kij)*gammaij;
         
                   // No B.C.
                   else if (bc_left == 2)
                     left_r0(i,j) = (0.5 - iphase*jphase - Kij)*gammaij;
                 }
                 
                 // Constraints on the derivative
                 if (bc_right == 4 || bc_right == 5)
                 {
                   right_rp(i,j) = -0.5*(iphase + kphase / this->k)*gammaij;

                   double phi_tmpj_right = 0;
                   for (int l=0; l<this->k; l++) {
                       double gammalj = sqrt(double((2*l+1)*(2*j+1)));
                       double Klj;
                       if (((l-j)>0) && (((l-j)%2)==1))  Klj = 2.0;
                       else   Klj = 0.0;
                       phi_tmpj_right += sqrt(double(2*l+1))*Klj*gammalj;
                    }
                    right_r0(i,j) = -(0.5*jphase*(iphase+ kphase/this->k) + Kij)*gammaij + sqrt(double(2*i+1))*phi_tmpj_right/pow(this->k,2.) ;
                 }
                 else if (bc_right == 0 || bc_right == 2 || bc_right == 3)
                 {
                   right_rp(i,j) = rp(i,j) ;
                 
                   // Zero BC
                   if (bc_right == 0 || bc_right == 3)
                     right_r0(i,j) = -(0.5*iphase*jphase + Kij)*gammaij;
                 
                   // No BC
                   else if (bc_right == 2)
                     right_r0(i,j) = (1.0 - 0.5*iphase*jphase - Kij)*gammaij; 
                 
                 }

                 jphase = -jphase;
             }
             iphase = -iphase;
         }


         // These are the coefficients for the boundary contribution
         bv_left  = Tensor<double>(this->k) ;
         bv_right = Tensor<double>(this->k) ;

         iphase = 1.0;
         for (int i=0; i<this->k; i++) {
             iphase = -iphase;

             if (bc_left == 3)
               bv_left(i) = iphase*sqrt(double(2*i+1));            // vector for left dirichlet BC
             else if(bc_left == 5)
               bv_left(i) = -iphase*sqrt(double(2*i+1))/pow(this->k,2);  // vector for left deriv BC
             else
               bv_left(i) = 0.0 ;

             if (bc_right == 3)
               bv_right(i) = sqrt(double(2*i+1));                  // vector for right dirichlet BC
             else if (bc_right == 5)
               bv_right(i) = sqrt(double(2*i+1))/pow(this->k,2);         // vector for right deriv BC 
             else
               bv_right(i) = 0.0 ;
         }

         return ;
      }
 
  } ;




    namespace archive {
        template <class Archive, class T, int NDIM>
        struct ArchiveLoadImpl<Archive,const TreeTraversal<T,NDIM>*> {
            static inline void load(const Archive& ar, const TreeTraversal<T,NDIM>*& ptr) {
                WorldObject< TreeTraversal<T,NDIM> >* p;
                ar & p;
                ptr = static_cast< const TreeTraversal<T,NDIM>* >(p);
            }
        };

    }

/*
  template <typename T, int NDIM>
  class FreeSpaceDerivative : public Deriative<T, NDIM> {

      FreeSpaceDerivative(World& _world, int _k, int _axis) : 
         TreeTraversal<T, NDIM>(_world, _k, _axis, 1) ,
         g1() ,
         g2() , 
      {
      }

  } ;
*/

}  // End of the madness namespace

#endif // MADNESS_MRA_DERIVATIVE_H_INCLUDED
