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

#include<madness/chem/CalculationParameters.h>

namespace madness {

/// struct for running a protocol of subsequently tightening precision

/// The ladder is `CalculationParameters::protocol()` and nothing else -- the
/// same list moldft walks, so both engines refine through identical steps and
/// `k` follows from `SCF::set_protocol`'s thresh->k table. `econv` and `dconv`
/// are convergence criteria applied at each rung, relaxed to the rung's own
/// threshold the way `SCF::solve` already relaxes dconv: you cannot converge a
/// quantity tighter than the basis represents it.
class SCFProtocol {
public:
    SCFProtocol(World& w, const CalculationParameters& param)
            : world(w), converged(false), protocol(param.protocol()), index(0),
              user_econv(param.econv()), user_dconv(param.dconv()) {
        MADNESS_CHECK_THROW(not protocol.empty(), "empty protocol in SCFProtocol");
        start_prec = protocol.front();
        end_prec = protocol.back();
        current_prec = start_prec;
        infer_thresholds(current_prec);
    }

    World& world;

    bool converged;         ///< flag if protocol has converged

    double start_prec;      ///< starting precision, typically 1.e-4
    double current_prec;    ///< current precision
    double end_prec;        ///< final precision

    double thresh;          ///< numerical precision of representing functions
    double econv;           ///< energy convergence of SCF calculations
    double dconv;           ///< density convergence of SCF calculations

    /// number of rungs in the ladder
    std::size_t size() const {return protocol.size();}

    /// index of the rung this protocol will start at
    std::size_t start_index() const {return index;}

    /// start at a given rung rather than at the first

    /// used to skip rungs a restart has already converged through
    void set_start_index(const std::size_t i) {
        MADNESS_CHECK_THROW(i<protocol.size(), "start index beyond the protocol");
        index=i;
        start_prec=protocol[index];
        current_prec=start_prec;
        infer_thresholds(current_prec);
    }

    /// drop the last rung, i.e. stop one step short of the full precision

    /// used for cheap pre-iterations; never empties the ladder
    void drop_last_rung() {
        if (protocol.size()<2) return;
        protocol.pop_back();
        end_prec=protocol.back();
        if (index>=protocol.size()) set_start_index(protocol.size()-1);
    }

    void initialize() {

        // don't do anything if this protocol is already converged
        if (converged) return;

        current_prec=protocol[index];
        infer_thresholds(current_prec);

        if (world.rank()==0) {
            std::stringstream ss;
            ss <<"\nstarting protocol at time" << std::setw(8) << std::setprecision(2)
               << wall_time() << "s";
            print(ss.str());
            print("precision steps ",current_prec," --> ",end_prec,
                  " (rung",index,"of",protocol.size(),")");
            print("protocol: thresh",thresh,"econv ",econv,"dconv",dconv);
        }
    }

    bool finished() const {return converged;}

    /// go to the next rung of the ladder
    SCFProtocol& operator++() {
        if (index+1<protocol.size()) {
            ++index;
            current_prec=protocol[index];
            infer_thresholds(current_prec);
        } else {
            converged=true;
        }

        return *this;
    }

    /// true if the current rung is the last one
    bool on_last_rung() const {return index+1==protocol.size();}

    /// infer thresholds for a given rung of the ladder

    /// The rung sets the representation threshold. econv and dconv are the
    /// user's, but never tighter than the rung can support.
    ///
    /// dconv needs care. The BSH residual cannot fall much below the threshold
    /// the orbitals are represented at, and nemo tests `bsh_norm < dconv`
    /// strictly, so a dconv equal to the rung's threshold is unreachable and the
    /// rung burns maxiter without converging. Intermediate rungs therefore use
    /// the long-standing 0.1*sqrt(thresh) relaxation -- 1e-3 at thresh 1e-4,
    /// 1e-4 at thresh 1e-6 -- which is comfortably looser than the rung; they
    /// are only a stepping stone, so there is nothing to gain from converging
    /// them tightly. The last rung is the answer, so it honours the user's dconv
    /// (never demanding tighter than the representation supports).
    void infer_thresholds(const double prec) {
        thresh=prec;
        econv=std::max(prec,user_econv);
        if (on_last_rung()) dconv=std::max(prec,user_dconv);
        else dconv=std::max(user_dconv,std::min(1.e-3,sqrt(prec)*0.1));
    }

    /// compare two positive doubles to be equal
    bool approx(const double a, const double b) const {
        return (std::abs(a/b-1.0)<1.e-12);
    }

private:
    std::vector<double> protocol;   ///< the ladder, from CalculationParameters
    std::size_t index;              ///< current rung
    double user_econv;              ///< energy convergence provided by user
    double user_dconv;              ///< density convergence provided by user
};



} // namespace madness


#endif /* SRC_APPS_CHEM_SCFPROTOCOL_H_ */
