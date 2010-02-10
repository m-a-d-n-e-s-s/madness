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

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/sdf_domainmask.h>

using namespace madness;

/** \brief Li, Lowengrub, Ratz, Voight domain masking with a Gaussian for the
    surface function. */
class LLRVGaussianDomainMask : public DomainMaskInterface {
private:
    LLRVGaussianDomainMask() : epsilon(0.0) {} ///< Forbidden
        
protected:
    const double epsilon; ///< The width of the transition region
    
public:
    /** \brief Constructor for the domain mask

        \param[in] epsilon The effective width of the surface */
    LLRVGaussianDomainMask(double epsilon) 
        : epsilon(epsilon)
    {}

    /** \brief Value of characteristic function at normal distance d from
               the surface

        \param[in] d The signed distance.  Negative is ``inside,''
                     positive is ``outside.''
        \return The domain mask */
    double mask(double d) const {
        if (d > 8.0*epsilon) {
            return 0.0; // we're safely outside
        }
        else if (d < -8.0*epsilon) {
            return 1.0; // inside
        }
        else {
            return 0.5 * (1.0 - tanh(3.0 * d / epsilon));
        }
    }
    
    /** \brief Derivative of characteristic function with respect to the
               normal distance

        \param[in] d The signed distance
        \return The derivative */
    double dmask(double d) const {
        if (fabs(d) > 8.0*epsilon) {
            return 0.0; // we're safely outside or inside
        }
        else {
            double tanh3d = tanh(3.0*d/epsilon);
            return 1.5*(tanh3d*tanh3d - 1.0) / epsilon;
        }
    }
    
    /** \brief Value of surface function at distance d normal to surface

        \param[in] d The signed distance
        \return The value of the surface function */
    double surface(double d) const {
        return exp(-d*d*0.5/(epsilon*epsilon)) / (sqrt(2.0*constants::pi)
            * epsilon);
    }

    /** \brief Value of d(surface)/ddistance

        \param[in] d The signed distance
        \return The derivative of the surface function */
    double dsurface(double d) const {
        double phi = mask(d);
        double dphi = dmask(d);
        return 72.0*phi*(1.0-phi)*dphi*(1.0 - 2.0*phi)/epsilon;
    }
    
    virtual ~LLRVGaussianDomainMask() {}
};
