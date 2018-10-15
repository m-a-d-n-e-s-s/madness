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
  \file chem/pointgroupsymmetry.cc
  \brief implements point group operations

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
  /trunk/src/apps/chem/pointgroupsymmetry.cc>here</a>.

*/

#include <chem/pointgroupsymmetry.h>

namespace madness {



std::vector<std::string> projector_irrep::reduce(const std::vector<std::string> irreps) const {

	bool verbose=false;
	for (const std::string& r : irreps) {
		if (r=="null") return std::vector<std::string> (1,"null");
	}

	// make sure all reps in the input vector exist in the table
	for (auto r : irreps) {
		if (table_.irreps_.find(r)==table_.irreps_.end()) {
			print_character_table();
			print("Irrep", r);
			MADNESS_EXCEPTION("Could not find irrep in table",1);
		}
	}
	if (verbose) {
		print("input irreps in reduce");
		for (auto& i : irreps) print(i);
	}

	// the characters of the direct product of the input representation
	// is the product of their characters
	std::vector<int> reducible(table_.irreps_.find(irreps[0])->second);
	for (std::size_t irep=1; irep<irreps.size(); ++irep) {
		std::vector<int> rep_characters(table_.irreps_.find(irreps[irep])->second);
		for (std::size_t ichar=0; ichar<rep_characters.size(); ++ichar) {
			reducible[ichar]*=rep_characters[ichar];
		}
	}
	if (verbose) {
		print("characters of the reducible product representation");
		for (const int i : reducible) print(i);
	}

	// reduce the reducible representation
	std::vector<std::string> result;
	for (const std::string& irrep : get_all_irreps()) {
		int n_irrep=0;
		std::vector<int> ii=table_.irreps_.find(irrep)->second;
		for (int i=0; i<reducible.size(); ++i) {	// sum over all classes/operators
			n_irrep+=reducible[i]*ii[i];
		}
		MADNESS_ASSERT(n_irrep%table_.order_==0);
		n_irrep/=table_.order_;
		if (verbose) {
			print("found irrep",irrep, n_irrep," times");
		}
		if (n_irrep>1) MADNESS_EXCEPTION("cannot handle the same irrep multiple times in reduce",1);
		if (n_irrep==1) result.push_back(irrep);
	}
	return result;
}


charactertable projector_irrep::make_c1_table() const {
	charactertable c1;
	c1.schoenflies_="C1";
	c1.order_=1;
	c1.operators_.push_back(pg_identity());
	c1.irreps_["a"]=vector_factory<int>(1);
	c1.mullikan_=vector_factory<std::string>("a");
	return c1;
}

charactertable projector_irrep::make_cs_table() const {
	charactertable cs;
	cs.schoenflies_="Cs";
	cs.order_=2;
	cs.operators_.push_back(pg_identity());
	cs.operators_.push_back(pg_sigma_xy());
	cs.irreps_["a'"]=vector_factory<int>(1,1);
	cs.irreps_["a''"]=vector_factory<int>(1,-1);
	cs.mullikan_=vector_factory<std::string>("a'","a''");
	return cs;
}


charactertable projector_irrep::make_ci_table() const {
	charactertable cs;
	cs.schoenflies_="Ci";
	cs.order_=2;
	cs.operators_.push_back(pg_identity());
	cs.operators_.push_back(pg_inversion());
	cs.mullikan_=vector_factory<std::string>("ag","au");
	cs.irreps_["ag"]=vector_factory<int>(1,1);
	cs.irreps_["au"]=vector_factory<int>(1,-1);
	return cs;
}

charactertable projector_irrep::make_c2_table() const {
	charactertable cs;
	cs.schoenflies_="C2";
	cs.order_=2;
	cs.operators_.push_back(pg_identity());
	cs.operators_.push_back(pg_c2z());
	cs.mullikan_=vector_factory<std::string>("a","b");
	cs.irreps_["a"]=vector_factory<int>(1,1);
	cs.irreps_["b"]=vector_factory<int>(1,-1);
	return cs;
}

charactertable projector_irrep::make_c2v_table() const {
	charactertable cs;
	cs.schoenflies_="C2v";
	cs.order_=4;
	cs.operators_.push_back(pg_identity());
	cs.operators_.push_back(pg_c2z());
	cs.operators_.push_back(pg_sigma_xz());
	cs.operators_.push_back(pg_sigma_yz());
	cs.mullikan_=vector_factory<std::string>("a1","a2","b1","b2");
	cs.irreps_["a1"]=vector_factory<int>(1,1,1,1);
	cs.irreps_["a2"]=vector_factory<int>(1,1,-1,-1);
	cs.irreps_["b1"]=vector_factory<int>(1,-1,1,-1);
	cs.irreps_["b2"]=vector_factory<int>(1,-1,-1,1);
	return cs;
}


charactertable projector_irrep::make_c2h_table() const {
	charactertable cs;
	cs.schoenflies_="C2h";
	cs.order_=4;
	cs.operators_.push_back(pg_identity());
	cs.operators_.push_back(pg_c2z());
	cs.operators_.push_back(pg_inversion());
	cs.operators_.push_back(pg_sigma_xy());
	cs.mullikan_=vector_factory<std::string>("ag","bg","au","bu");
	cs.irreps_["ag"]=vector_factory<int>(1, 1, 1, 1);
	cs.irreps_["bg"]=vector_factory<int>(1,-1, 1,-1);
	cs.irreps_["au"]=vector_factory<int>(1, 1,-1,-1);
	cs.irreps_["bu"]=vector_factory<int>(1,-1,-1, 1);
	return cs;
}

charactertable projector_irrep::make_d2_table() const {
	charactertable cs;
	cs.schoenflies_="D2";
	cs.order_=4;
	cs.operators_.push_back(pg_identity());
	cs.operators_.push_back(pg_c2z());
	cs.operators_.push_back(pg_c2y());
	cs.operators_.push_back(pg_c2x());
	cs.mullikan_=vector_factory<std::string>("a1","b1","b2","b3");
	cs.irreps_["a1"]=vector_factory<int>(1, 1, 1, 1);
	cs.irreps_["b1"]=vector_factory<int>(1, 1,-1,-1);
	cs.irreps_["b2"]=vector_factory<int>(1,-1, 1,-1);
	cs.irreps_["b3"]=vector_factory<int>(1,-1,-1, 1);
	return cs;
}

charactertable projector_irrep::make_d2h_table() const {
	charactertable cs;
	cs.schoenflies_="D2h";
	cs.order_=8;
	cs.operators_.push_back(pg_identity());
	cs.operators_.push_back(pg_c2z());
	cs.operators_.push_back(pg_c2y());
	cs.operators_.push_back(pg_c2x());
	cs.operators_.push_back(pg_inversion());
	cs.operators_.push_back(pg_sigma_xy());
	cs.operators_.push_back(pg_sigma_xz());
	cs.operators_.push_back(pg_sigma_yz());
	cs.mullikan_=vector_factory<std::string>("ag","au","b1g","b1u","b2g","b2u","b3g","b3u");
	cs.irreps_["ag"] =vector_factory<int>(1, 1, 1, 1, 1, 1, 1, 1);
	cs.irreps_["b1g"]=vector_factory<int>(1, 1,-1,-1, 1 ,1,-1,-1);
	cs.irreps_["b2g"]=vector_factory<int>(1,-1, 1,-1, 1,-1, 1,-1);
	cs.irreps_["b3g"]=vector_factory<int>(1,-1,-1, 1, 1,-1,-1, 1);
	cs.irreps_["au"] =vector_factory<int>(1, 1, 1, 1,-1,-1,-1,-1);
	cs.irreps_["b1u"]=vector_factory<int>(1, 1,-1,-1,-1,-1, 1, 1);
	cs.irreps_["b2u"]=vector_factory<int>(1,-1, 1,-1,-1, 1,-1, 1);
	cs.irreps_["b3u"]=vector_factory<int>(1,-1,-1, 1,-1, 1, 1,-1);
	return cs;
}

} /* namespace madness */
