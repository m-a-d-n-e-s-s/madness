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

/// \file corepotential.cc
/// \brief Simple management of core potential and orbital information

#include <madness_config.h>
#include <constants.h>
#include <mra/mra.h>
#include <tinyxml/tinyxml.h>
#include <moldft/corepotential.h>
#include <cstdio>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <utility>
#include <map>
#include <set>
using std::string;
using std::vector;
using namespace madness;

typedef Vector<double,3> coordT;
typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef FunctionFactory<double,3> factoryT;
typedef vector<functionT> vecfuncT;

static const string dir = "coredata/";

static bool read_potential(TiXmlElement* elem, AtomCore& ac, double eprec) {
    TiXmlElement* p = elem->FirstChildElement("potential");
    if (!p) return false;

    std::istringstream iss(p->GetText());
    int l, n;
    double e, c;
    vector<int> vl, vn;
    vector<double> ve, vc;
    while (iss >> l) {
        iss >> n >> e >> c;
        if (l<0) continue;
        vl.push_back(l);
        vn.push_back(n);
        ve.push_back(e);
        vc.push_back(c);
    }
    ac.potential.l = vl;
    ac.potential.n = vn;
    ac.potential.A = vc;
    ac.potential.alpha = ve;
    ac.potential.eprec = eprec;
    int atn = ac.atomic_number;
    int ncore = ac.ncore;
    //ac.potential.rcut0 = 1.0/smoothing_parameter(ncore*2, eprec);
    ac.potential.rcut0 = 1.0/smoothing_parameter(ncore*2, 1.0);
    //ac.potential.rcut = 1.0/smoothing_parameter(atn-ncore*2, eprec);
    ac.potential.rcut = 1.0/smoothing_parameter(atn-ncore*2, 1.0);

    return true;
}

static bool read_orbital(TiXmlElement* e, AtomCore& ac) {
    TiXmlElement* p = e->FirstChildElement("core");
    if (!p) return false;

    std::istringstream i_num(p->Attribute("num"));
    i_num >> ac.ncore;

    vector<CoreOrbital> vc;

    for (TiXmlElement* node = p->FirstChildElement("orbital"); node; node = node->NextSiblingElement("orbital")) {
        int type;
        vector<double> coeff, expnt;
        double c, e;
        double bc;
        std::istringstream i_bc(node->Attribute("bc"));
        i_bc >> bc;
        std::istringstream i_type(node->Attribute("type"));
        i_type >> type;
        std::istringstream iss(node->GetText());
        while (iss >> e) {
            iss >> c;
            coeff.push_back(c);
            expnt.push_back(e);
        }
        CoreOrbital co(type, coeff, expnt, bc);
        vc.push_back(co);
    }
    ac.orbital = vc;
    if (vc.size() != ac.ncore) {
        MADNESS_EXCEPTION("CORE_INFO: read_orbital: inconsistent number of core.", -1);
    }

    return true;
}

static AtomCore read_atom(TiXmlElement* e, unsigned int atn, double eprec) {
    AtomCore ac;
    ac.atomic_number = atn;

    if (!read_orbital(e, ac)) {
        MADNESS_EXCEPTION("CORE_INFO: read_orbital failed.", -1);
    }
    if (!read_potential(e, ac, eprec)) {
        MADNESS_EXCEPTION("CORE_INFO: read_potential failed.", -1);
    }

    return ac;
}

void CorePotentialManager::read_file(string filename, std::set<unsigned int> atomset, double eprec) {
    core_type = filename;
    guess_filename = dir + filename + "_guess";

    TiXmlDocument doc(dir + core_type);
    if (!doc.LoadFile()) {
        MADNESS_EXCEPTION(("CORE_INFO: Failed to load core_info data file: " + dir + core_type).c_str(), -1);
        return;
    }
    TiXmlElement* core_info = doc.FirstChildElement();
    if (!core_info) {
        MADNESS_EXCEPTION("CORE_INFO: core_info data file is not valid.", -1);
        return;
    }

    for (TiXmlElement* node = core_info->FirstChildElement("atom"); node; node = node->NextSiblingElement("atom")) {
        unsigned int atn = symbol_to_atomic_number(node->Attribute("symbol"));
        if (atomset.find(atn) != atomset.end()) {
            AtomCore ac = read_atom(node, atn, eprec);
            if (ac.n_orbital() == 0) {
                MADNESS_EXCEPTION("CORE_INFO: read_atom Failed.", -1);
                return;
            }
            atom_core.insert(std::pair<unsigned int,AtomCore>(atn, ac));
        }
    }

    vector<unsigned int> atns;
    for (std::map<unsigned int, AtomCore>::iterator it = atom_core.begin(); it != atom_core.end(); it++) {
        atns.push_back(it->first);
    }
    madness::print("MCP parameters loaded for atomic numbers:", atns);
}

