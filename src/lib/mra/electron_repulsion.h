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

  $Id$
*/


#ifndef MADNESS_MRA_ELECTRON_REPULSION_H__INCLUDED
#define MADNESS_MRA_ELECTRON_REPULSION_H__INCLUDED

#include <tensor/tensor.h>
#include <mra/convolution1d.h>
#include <mra/key.h>
#include <mra/funcdefaults.h>

/// Holds construction of the 6D MRA electron repulsion function 1/r12

/// The MRA function 1/r12 is too large to be kept in memory so we need to
/// construct single function nodes on the fly. However, for any internal
/// node the Gauss-Legendre quadrature is expensive and can be very inaccurate,
/// so we simply reuse the matrix elements of the kernel of the Coulomb
/// operator.

/// need to optimize, i.e. include the norm when computing the matrix elements

namespace madness {

	// how do I do this right??
	void flo_bsh_fit(double mu, double lo, double hi, double eps,
        Tensor<double> *pcoeff, Tensor<double> *pexpnt, bool prnt=false);


    template <std::size_t NDIM>
    class ElectronRepulsion {

    private:

        mutable std::vector< ConvolutionND<double,3> > ops;
        BoundaryConditions<6> bc;
        int rank;
        int k;
        double thresh;

    public:

        /// don't do nothing
        ElectronRepulsion() {}

    	/// constructor: cf the Coulomb kernel
    	ElectronRepulsion(World& world,double lo,double eps,
                const BoundaryConditions<6>& bc=FunctionDefaults<6>::get_bc(),
                int kk=FunctionDefaults<NDIM>::get_k())
    		: thresh(eps) {

    		// don't think anything else makes sense..
    		MADNESS_ASSERT(NDIM==6);

    		// Presently we must have periodic or non-periodic in all dimensions.
            for (std::size_t d=1; d<NDIM; ++d) {MADNESS_ASSERT(bc(d,0)==bc(0,0));}
            const Tensor<double>& width = FunctionDefaults<3>::get_cell_width();


            double hi = width.normf(); // Diagonal width of cell
            if (bc(0,0) == BC_PERIODIC) hi *= 100; // Extend range for periodic summation
            const double pi = constants::pi;

            // bsh_fit generates representation for 1/4Pir but we want 1/r
            // so have to scale eps by 1/4Pi
            Tensor<double> coeff, expnt;
            flo_bsh_fit(0.0, lo, hi, eps/(4.0*pi), &coeff, &expnt, false);
            coeff.scale(4.0*pi);

            // set some parameters
            rank=coeff.dim(0);
            ops.resize(rank);
            k=kk;

            // construct all the terms
            for (int mu=0; mu<rank; ++mu) {
//                double c = std::pow(sqrt(expnt(mu)/pi),static_cast<int>(NDIM)); // Normalization coeff
                double c = std::pow(sqrt(expnt(mu)/pi),3); // Normalization coeff

                // We cache the normalized operator so the factor is the value we must multiply
                // by to recover the coeff we want.
                ops[mu].setfac(coeff(mu)/c);

                // only 3 dimensions here!
                for (std::size_t d=0; d<3; ++d) {
                	ops[mu].setop(d,GaussianConvolution1DCache<double>::get(k, expnt(mu)*width[d]*width[d], 0, false));
                }
            }
    	}

    	/// return the coefficients of the function in 6D (x1,y1,z1, x2,y2,z2)
    	Tensor<double> coeff(const Key<NDIM>& key) const {

    		MADNESS_ASSERT(NDIM==6);
    		Tensor<double> c(k,k,k,k,k,k);

    		const Level n=key.level();
    		const Vector<Translation,NDIM> l=key.translation();

    		// get the displacements for all 3 dimensions: x12, y12, z12
    		const Translation l0=(l[0]-l[3]);
    		const Translation l1=(l[1]-l[4]);
    		const Translation l2=(l[2]-l[5]);

    		// the dimensions are a bit confused (x1,x2, y1,y2, z1,z2) -> (x1,y1,z1, x2,y2,z2)
    		std::vector<long> map(NDIM);
   			map[0]=0;
   			map[1]=3;
   			map[2]=1;
   			map[3]=4;
   			map[4]=2;
   			map[5]=5;

    		// accumulate the terms r^n_ll'_kk'
    		for (int mu=0; mu<rank; mu++) {

    			const Tensor<double>& r0=ops[mu].getop(0)->rnlij(n,l0);
    			const Tensor<double>& r1=ops[mu].getop(1)->rnlij(n,l1);
    			const Tensor<double>& r2=ops[mu].getop(2)->rnlij(n,l2);

    			const double norm0=r0.normf();
    			const double norm1=r1.normf();
    			const double norm2=r2.normf();
//    			if (norm0*norm1*norm2*ops[mu].getfac() > thresh/rank) {

					const Tensor<double> r01=outer(r0,r1*ops[mu].getfac());
					Tensor<double> r123=outer(r01,r2);
	//    			r123.scale(ops[mu].getfac());
					c+=r123.mapdim(map);
//    			}
    		}

    		return c;
    	}

    };
}

#endif // MADNESS_MRA_FLONODE_H__INCLUDED
