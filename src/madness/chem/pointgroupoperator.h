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

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES


/*!
  \file chem/pointgroupoperator.h
  \brief implements point group operations

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
  /trunk/src/apps/chem/pointgroupoperator.h>here</a>.

*/

#ifndef SRC_APPS_CHEM_POINTGROUPOPERATOR_H_
#define SRC_APPS_CHEM_POINTGROUPOPERATOR_H_

#include <madness/tensor/vector_factory.h>

namespace madness {

template<typename T, std::size_t NDIM>
class Function;
/// This class implements the symmetry operations (not the point groups)
class pg_operator {

public:

	/// ctor with explicit mirror maps
    pg_operator(std::string symbol, std::string name, const std::vector<long> mm,
    		const std::vector<long> md)
		: symbol_(symbol), name_(name), mirrormap(mm), mapdim_(md) {}

    /// return the name of the symmetry operator
    std::string name() const {return name_;}

    std::string symbol() const {return symbol_;}

	/// apply the operator on an n-dimensional MRA function
	template<typename T, std::size_t NDIM>
	Function<T,NDIM> operator()(const Function<T,NDIM>& f, bool fence=true) const {
        Function<T,NDIM> result;

        if (name_=="identity") {
        	result=copy(f,fence);

        } else if (name_=="inversion"){

        	std::vector<long> mm(NDIM,-1);
        	result=mirror(f,mm,fence);

        } else {
        	result=map_and_mirror(f,mapdim_,mirrormap,fence);
        }

		return result;
	}

	/// apply the operator on an n-dimensional MRA function
	template<typename T, std::size_t NDIM>
	std::vector<Function<T,NDIM> > operator()(const std::vector<Function<T,NDIM> >& vf, bool fence=true) const {
        std::vector<Function<T,NDIM> > result(vf.size());

        // fast return
        if (vf.size()==0) return result;
        World& world=vf.begin()->world();

        for (size_t i=0; i<vf.size(); ++i) result[i]=this->operator()(vf[i],false);

        if (fence) world.gop.fence();
		return result;
	}


	friend
	std::ostream& operator<<(std::ostream& s, const pg_operator& pg_op) {
		s << "Symmetry operator " << pg_op.name() ;
		return s;
	}


private:

	std::string symbol_;
	std::string name_;
	std::vector<long> mirrormap;
    std::vector<long> mapdim_;

};


// n-dimensional operations
static inline pg_operator pg_identity() {
    return pg_operator("E","identity",std::vector<long>(), std::vector<long>());
}

static inline pg_operator pg_inversion() {
    return pg_operator("i","inversion",std::vector<long>(), std::vector<long>());
}

// 2D operations
static inline pg_operator pg_sigma_x() {
    return pg_operator("s","sigma_x",vector_factory<long>(1,-1),std::vector<long>());
}

static inline pg_operator pg_sigma_y() {
    return pg_operator("s","sigma_y",vector_factory<long>(-1,1),std::vector<long>());
}

static inline pg_operator pg_c2() {
    return pg_operator("C2","C_2",vector_factory<long>(-1,-1),std::vector<long>());
}

static inline pg_operator pg_c4() {
    return pg_operator("C4","C_4",vector_factory<long>(-1,1),vector_factory<long>(1,0));
}


// 3D operations
static inline pg_operator pg_sigma_xy() {
    return pg_operator("s","sigma_xy",vector_factory<long>(1,1,-1), std::vector<long>());
}

static inline pg_operator pg_sigma_xz() {
    return pg_operator("s","sigma_xz",vector_factory<long>(1,-1,1), std::vector<long>());
}

static inline pg_operator pg_sigma_yz() {
    return pg_operator("s","sigma_yz",vector_factory<long>(-1,1,1), std::vector<long>());
}

static inline pg_operator pg_c2x() {
	return pg_operator("C2","C_2(x)",vector_factory<long>(1,-1,-1), std::vector<long>());
}

static inline pg_operator pg_c2y() {
	return pg_operator("C2","C_2(y)",vector_factory<long>(-1,1,-1), std::vector<long>());
}

static inline pg_operator pg_c2z() {
	return pg_operator("C2","C_2(z)",vector_factory<long>(-1,-1,1), std::vector<long>());
}

static inline pg_operator pg_c4x() {
	return pg_operator("C4","C_4(x)",vector_factory<long>(1,-1,1), vector_factory<long>(0,2,1));
}

static inline pg_operator pg_c4y() {
	return pg_operator("C4","C_4(y)",vector_factory<long>(1,1,-1), vector_factory<long>(2,1,0));
}

static inline pg_operator pg_c4z() {
	return pg_operator("C4","C_4(z)",vector_factory<long>(-1,1,1), vector_factory<long>(1,0,2));
}


} // namespace madness

#endif /* SRC_APPS_CHEM_POINTGROUPOPERATOR_H_ */
