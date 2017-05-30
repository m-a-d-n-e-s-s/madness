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

/// \file SCFProtocol.h
/// \brief solution protocol for SCF calculations


#ifndef MADNESS_CHEM_SCFPROTOCOL_H__INCLUDED
#define MADNESS_CHEM_SCFPROTOCOL_H__INCLUDED

namespace madness {

#include <chem/CalculationParameters.h>

/// struct for running a protocol of subsequently tightening precision
class SCFProtocol {
public:
    SCFProtocol(World& w, const CalculationParameters& param,
            const std::string name, const bool restart)
            : world(w), filename(name), restart(restart), start_prec(1.e-4),
            current_prec(start_prec), end_prec(param.econv),
            thresh(1.e-4), econv(1.e-4), dconv(1.e-3), user_dconv(1.e-20) {
        user_dconv=param.dconv;
    }

    World& world;

    std::string filename;   ///< name for the restart data file
    bool restart;           ///< do a restart from file (if possible)

    double start_prec;      ///< starting precision, typically 1.e-4
    double current_prec;    ///< current precision
    double end_prec;        ///< final precision

    double thresh;          ///< numerical precision of representing functions
    double econv;           ///< energy convergence of SCF calculations
    double dconv;           ///< density convergence of SCF calculations
    double user_dconv;      ///< density convergence provided by user

    void initialize() {

        // try to read restart data file
        if (restart) {
            std::ifstream f(filename.c_str());
            if(not f.fail()) {
                std::string s, tag;
                while (std::getline(f,s)) {
                    std::istringstream ss(s);
                    ss >> tag;
                    if (tag=="start_prec") ss >> start_prec;
                    if (tag=="end_prec") ss >> end_prec;
                    if (tag=="current_prec") ss >> current_prec;
                    if (tag=="econv") ss >> econv;
                    if (tag=="dconv") ss >> dconv;
                    if (tag=="thresh") ss >> thresh;
                    if (tag=="user_dconv") ss >> user_dconv;
                }
                if (world.rank()==0) print("tried restart but failed");
            }
        } else {
            current_prec=start_prec;
            infer_thresholds(current_prec);
        }
        if (world.rank()==0) {
            printf("\nstarting protocol at time %8.1fs \n",wall_time());
            print("precision steps ",start_prec," --> ",end_prec);
            print("protocol: thresh",thresh,"econv ",econv,"dconv",dconv);
        }
    }

    bool finished() const {
        return current_prec<0.0;
        return current_prec*0.9999<end_prec;   // account for noise
    }

    /// go to the next level
    SCFProtocol& operator++() {
        if (current_prec*0.9999>end_prec) {
            current_prec*=0.1;
            if (current_prec<end_prec) current_prec=end_prec;
            infer_thresholds(current_prec);
            if(world.rank()==0) print("protocol: thresh",thresh,"econv ",econv,"dconv",dconv);
        } else {
            current_prec=-1.0;
            if(world.rank()==0) printf("\nending protocol at time %8.1fs \n",wall_time());
        }

        // update restart data on file
        std::ofstream f(filename.c_str());
        f << "start_prec " << start_prec << std::endl;;
        f << "end_prec " << end_prec << std::endl;;
        f << "current_prec " << current_prec << std::endl;;
        f << "econv " << econv << std::endl;;
        f << "dconv " << dconv << std::endl;;
        f << "thresh " << thresh << std::endl;;
        f << "user_dconv " << user_dconv << std::endl;;

        return *this;
    }

    /// infer thresholds starting from a target precision
    void infer_thresholds(const double prec) {
        econv=prec;
        thresh=econv;
        dconv=std::min(1.e-3,sqrt(econv)*0.1);
//            dconv=std::min(1.e-3,econv*10.0);
        if (approx(current_prec,end_prec)) dconv=user_dconv;    // respect the user
    }

    /// compare two positive doubles to be equal
    bool approx(const double a, const double b) const {
        return (std::abs(a/b-1.0)<1.e-12);
    }
};



} // namespace madness


#endif /* SRC_APPS_CHEM_SCFPROTOCOL_H_ */
