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

  
#include <iostream>

#include <algorithm>

#include "basetensor.h"

#include "tensorexcept.h"

/// \file basetensor.cc
/// \brief Implements BaseTensor

namespace madness {

#ifdef TENSOR_INSTANCE_COUNT
    MADATOMIC_INT BaseTensor::instance_count = MADATOMIC_INITIALIZE(0);
#endif

    void BaseTensor::set_dims_and_size(long nd, const long d[]) {
        ndim = nd;
        size = 1;
        if (ndim < 0) size=0;
        for (long i=ndim-1; i>=0; i--) {
            dim[i] = d[i];
            stride[i] = size;
            size *= d[i];
        }
        for (long i=std::max(ndim,0L); i<TENSOR_MAXDIM; i++) { // So can iterate over missing dimensions
            dim[i] = 1;
            stride[i] = 0;
        }
    }


    bool BaseTensor::conforms(const BaseTensor *t) const {
        /// Returns true if *this and *t are the same shape and size
        if (ndim != t->ndim) return false;
        for (long i=0; i<ndim; i++) {
            if (dim[i] != t->dim[i]) return false;
        }
        return true;
    }

    // bool BaseTensor::iscontiguous() const {
    //   /// Returns true if the tensor refers to contiguous memory locations.
    //   long size = 1;
    //   for (long i=ndim-1; i>=0; i--) {
    //     if (stride[i] != size) return false;
    //     size *= dim[i];
    //   }
    //   return true;
    // }

    void BaseTensor::reshape_inplace_base(const std::vector<long>& d) {
        /// Reshape the size and number of dimensions.

        /// Modifies the current tensor to have the number and size of
        /// dimensions as described in the \c vector \c d .  The total number
        /// of elements must be the same before and after, and the current
        /// tensor must be contiguous.

        TENSOR_ASSERT(iscontiguous(),
                      "cannot reshape non-contiguous tensor ... consider fuse/splitdim",
                      0,this);
        long newsize=1;
        int nd = d.size();
        for (long i=0; i<nd; i++) newsize *= d[i];
        TENSOR_ASSERT(size == newsize,"old and new sizes do not match",size,this);
        set_dims_and_size(nd,&(d[0]));
    }

    BaseTensor* BaseTensor::reshape_base(const std::vector<long>& d) const {
        /// Return a new view reshaping the size and number of dimensions.

        /// The new tensor will have the number and size of dimensions as
        /// described in the \c vector \c d .  The total number of elements
        /// must be the same before and after, and the current tensor must
        /// be contiguous.
        BaseTensor* result = new_shallow_copy_base();
        result->reshape_inplace_base(d);
        return result;
    }

    void BaseTensor::flat_inplace_base() {
        /// Reshape the current tensor to be the same size and 1-d. It must be contiguous.
        TENSOR_ASSERT(iscontiguous(),"not contiguous",0,this);
        long d[] = {size};
        set_dims_and_size(1,d);
    }

    BaseTensor* BaseTensor::flat_base() const {
        /// Return a 1-d view of the current tensor ... it must be contiguous.
        BaseTensor* result = new_shallow_copy_base();
        result->flat_inplace_base();
        return result;
    }

    void BaseTensor::splitdim_inplace_base(long i, long dimi0, long dimi1) {
        /// Split dimension i in two ... the product of the new dimensions must match the old.
        if (i < 0) i += ndim;
        TENSOR_ASSERT(i>=0 && i<ndim, "invalid dimension", i, this);
        TENSOR_ASSERT(dimi0*dimi1 == dim[i], "before & after sizes do not match",
                      dim[i], this);
        TENSOR_ASSERT(ndim+1 <= TENSOR_MAXDIM, "resulting tensor has too many dimensions",
                      ndim+1, this);
        for (long j=ndim-1; j>i; j--) {
            dim[j+1] = dim[j];
            stride[j+1] = stride[j];
        }
        dim[i+1] = dimi1;
        stride[i+1] = stride[i];
        dim[i] = dimi0;
        stride[i] *= dimi1;
        ndim++;
    }

    BaseTensor* BaseTensor::splitdim_base(long i, long dimi0, long dimi1) const {
        /// Split dimension i in two ... the product of the new dimensions must match the old.
        BaseTensor* result = new_shallow_copy_base();
        result->splitdim_inplace_base(i, dimi0, dimi1);
        return result;
    }

    void BaseTensor::fusedim_inplace_base(long i) { // fuse i,i+1 --> i
        /// Fuse the contiguous dimensions i and i+1.
        if (i < 0) i += ndim;
        TENSOR_ASSERT(i>=0 && i<(ndim-1) && ndim>1,"invalid dimension",i,this);
        TENSOR_ASSERT(stride[i] == dim[i+1]*stride[i+1],"dimensions are not contiguous",
                      i, this);
        dim[i] *= dim[i+1];
        stride[i] = stride[i+1];
        for (long j=i+1; j<=ndim-1; j++) {
            dim[j] = dim[j+1];
            stride[j] = stride[j+1];
        }
        ndim--;
        dim[ndim] = 1;		// So can iterate over missing dimensions
        stride[ndim] = 0;
    }

    BaseTensor* BaseTensor::fusedim_base(long i) const {
        /// Fuse the contiguous dimensions i and i+1.
        if (i < 0) i += ndim;
        BaseTensor* result = new_shallow_copy_base();
        result->fusedim_inplace_base(i);
        return result;
    }

    void BaseTensor::swapdim_inplace_base(long i, long j) {
        /// Swap the dimensions i and j.
        if (i < 0) i += ndim;
        if (j < 0) j += ndim;
        TENSOR_ASSERT(i>=0 && i<ndim,"invalid dimension i",i,this);
        TENSOR_ASSERT(j>=0 && j<ndim,"invalid dimension j",j,this);
        std::swap<long>(dim[i],dim[j]);
        std::swap<long>(stride[i],stride[j]);
    }

    BaseTensor* BaseTensor::swapdim_base(long idim, long jdim) const {
        /// Swap the dimensions i and j.
        BaseTensor* result = new_shallow_copy_base();
        result->swapdim_inplace_base(idim, jdim);
        return result;
    }

    void BaseTensor::cycledim_inplace_base(long nshift, long start, long end) {
        /// Cyclic shift by nshift places of the dimensions [start,....,end]
	/// right shift?

        /// start and end specify inclusive range of dimensions ...
        /// (PREVIOUS behavior was Pythonic and end was one larger.)
        long ndshift, dimtmp[TENSOR_MAXDIM], stridetmp[TENSOR_MAXDIM];
        if (start < 0) start += ndim; // Same convention for -ve as in Slice
        if (end < 0) end += ndim;
        TENSOR_ASSERT(start>=0 && start<ndim,"invalid start dimension",start,this);
        TENSOR_ASSERT(end>=0 && end>=start,"invalid end dimension",end,this);

        ndshift = end - start + 1;
        for (long i=start; i<=end; i++) {
            dimtmp[i] = dim[i];
            stridetmp[i] = stride[i];
        }
        for (long i=end; i>=start; i--) {
            long j = i + nshift;
            while (j > end) j -= ndshift;
            while (j < start) j += ndshift;
            dim[j] = dimtmp[i];
            stride[j] = stridetmp[i];
        }
    }

    BaseTensor* BaseTensor::cycledim_base(long nshift, long start, long end) const {
        /// Cyclic shift by nshift places of the dimensions [start,....,end]

        /// note defaults in prototype start=0 end=-1
        BaseTensor* result = new_shallow_copy_base();
        result->cycledim_inplace_base(nshift, start, end);
        return result;
    }

    void BaseTensor::mapdim_inplace_base(const std::vector<long>& map) {
        /// General permuation of the dimensions
        TENSOR_ASSERT(ndim == (int) map.size(),"map[] must include all dimensions",map.size(),this);
        long tmpd[TENSOR_MAXDIM], tmps[TENSOR_MAXDIM];
        for (long i=0; i<ndim; i++) {
            tmpd[map[i]] = dim[i];
            tmps[map[i]] = stride[i];
        }
        for (long i=0; i<ndim; i++) {
            dim[i] = tmpd[i];
            stride[i] = tmps[i];
        }
    }

    BaseTensor* BaseTensor::mapdim_base(const std::vector<long>& map) const {
        /// General permuation of the dimensions
        BaseTensor* result = new_shallow_copy_base();
        result->mapdim_inplace_base(map);
        return result;
    }

}
