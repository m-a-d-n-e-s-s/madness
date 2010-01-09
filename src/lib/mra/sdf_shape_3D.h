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
/*!
  \file mra/sdf_shape.h
  \brief An abstract signed distance function for use when forming subdomains.
  \ingroup mrabcint

  This file defines SDF classes for common 3-D surfaces:
  - Plane
  - Sphere
  - Cone
  - Paraboloid
  - Box
  - Cube
  - Ellipsoid
  - Cylinder
  
  NOTE: The signed distance functions would be the shortest distance between
  a point and _any_ point on the surface.  This is hard to calculate in many
  cases, so we use contours here.  The surface layer may not be equally thick
  around all points on the surface.
*/  
  
#ifndef MADNESS_MRA_SDF_SHAPE_3D_H__INCLUDED
#define MADNESS_MRA_SDF_SHAPE_3D_H__INCLUDED

#include <mra/sdf_shape.h>

namespace madness {

    /// A plane surface (3 dimensions)
    class SDFPlane : public ShapeDFInterface<3> {
    protected:
        const coord_3d normal; ///< The normal vector pointing OUTSIDE the surface
        const coord_3d point; ///< A point in the plane

    public:
        
        /// SDF for a plane transecting the entire simulation volume

        /// @param normal The outward normal definining the plane
        /// @param point A point in the plane
        SDFPlane(const coord_3d& normal, const coord_3d& point)
            : normal(normal*(1.0/sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2])))
            , point(point) 
        {}

        /// Computes the normal distance

        /// @param pt Point at which to compute the distance from the surface
        /// @return The signed distance from the surface
        double sdf(const coord_3d& pt) const {
            return (pt[0]-point[0])*normal[0] + (pt[1]-point[1])*normal[1] + (pt[2]-point[2])*normal[2];
        }

        coord_3d grad_sdf(const coord_3d& pt) const {
            MADNESS_EXCEPTION("gradient method is not yet implemented for this shape",0);
        }
    };

    /// A spherical surface (3 dimensions)
    class SDFSphere : public ShapeDFInterface<3> {
    protected:
        const double radius; ///< Radius of sphere
        const coord_3d center; ///< Center of sphere

    public:
        /// SDF for a sphere partially or fully enclosed in the solution volume

        /// @param radius The radius of the sphere
        /// @param center The center of the sphere
        SDFSphere(const double radius, const coord_3d &center) 
            : radius(radius)
            , center(center) 
        {}

        double sdf(const coord_3d& pt) const {
            double temp, r;
            int i;
            
            r = 0.0;
            for(i = 0; i < 3; ++i) {
                temp = pt[i] - center[i];
                r += temp * temp;
            }
            
            return sqrt(r) - radius;
        }

        coord_3d grad_sdf(const coord_3d& pt) const {
            double x = pt[0] - center[0];
            double y = pt[1] - center[1];
            double z = pt[2] - center[2];
            double r = sqrt(x*x + y*y + z*z);
            coord_3d g;
            g[0] = x/r;
            g[1] = y/r;
            g[2] = z/r;
            return g;
        }
    };

    /// A cone (3 dimensions)

    /// The cone is defined by 
    /// \f[
    /// \sqrt{x^2 + y^2} - c * z = 0
    /// \f]
    /// where \f$ z \f$ is along the axis.
    class SDFCone : public ShapeDFInterface<3> {
    protected:
        const coord_3d apex; ///< The apex
        const double c; ///< The radius
        const coord_3d dir; ///< The direction of the axis, from the apex INSIDE

    public:
        /// Constructor for cone

        /// @param c Parameter \f$ c \f$ in the definition of the cone
        /// @param apex Apex of cone
        /// @param direc Oriented axis of the cone
        SDFCone(const double c, const coord_3d &apex, const coord_3d &direc) 
            : apex(apex)
            , c(c)
            , dir(direc*(1.0/sqrt(direc[0]*direc[0] + direc[1]*direc[1] + direc[2]*direc[2])))
        {}

        double sdf(const coord_3d& pt) const {
            coord_3d diff;
            unsigned int i;
            double dotp;
            
            for(i = 0; i < 3; ++i)
                diff[i] = pt[i] - apex[i];
            
            dotp = diff[0]*dir[0] + diff[1]*dir[1] + diff[2]*dir[2];
            
            for(i = 0; i < 3; ++i)
                diff[i] -= dotp * dir[i];
            
            return sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2])
                - c * dotp;
        }

        coord_3d grad_sdf(const coord_3d& pt) const {
            MADNESS_EXCEPTION("gradient method is not yet implemented for this shape",0);
        }
    };

    /// A paraboloid (3 dimensions)

    /// The surface is defined by 
    /// \f[
    ///  x^2 + y^2 - c * z == 0
    /// \f]
    /// where \f$ z \f$ is along the axis.
    class SDFParaboloid : public ShapeDFInterface<3> {
    protected:
        const coord_3d apex; ///< The apex
        const double c; ///< Curvature/radius of the surface
        const coord_3d dir;///< The direction of the axis, from the apex INSIDE

    public:
        SDFParaboloid(const double c, const coord_3d &apex, const coord_3d &direc) 
            : apex(apex)
            , c()
            , dir(direc*(1.0/sqrt(direc[0]*direc[0] + direc[1]*direc[1] + direc[2]*direc[2])))
        {}

        double sdf(const coord_3d& pt) const {
            coord_3d diff;
            unsigned int i;
            double dotp;

            for(i = 0; i < 3; ++i)
                diff[i] = pt[i] - apex[i];
            
            dotp = diff[0]*dir[0] + diff[1]*dir[1] + diff[2]*dir[2];
            
            for(i = 0; i < 3; ++i)
                diff[i] -= dotp * dir[i];
            
            return diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]
                - c * dotp;
        }

        coord_3d grad_sdf(const coord_3d& pt) const {
            MADNESS_EXCEPTION("gradient method is not yet implemented for this shape",0);
        }
    };

    /// A box (3 dimensions)

    /// LIMIT: the 3 primary axes must be x, y, and z
    class SDFBox : public ShapeDFInterface<3> {
    protected:
        const coord_3d lengths;  ///< Half the length of each side of the box
        const coord_3d center; ///< the center of the box

    public:
        SDFBox(const coord_3d& length, const coord_3d& center) 
            : lengths(length*0.5), center(center) 
        {}

        double sdf(const coord_3d& pt) const {
            double diff, max;
            int i;
            
            max = fabs(pt[0] - center[0]) - lengths[0];
            for(i = 1; i < 3; ++i) {
                diff = fabs(pt[i] - center[i]) - lengths[i];
                if(diff > max)
                    max = diff;
            }
            
            return max;
        }

        coord_3d grad_sdf(const coord_3d& pt) const {
            MADNESS_EXCEPTION("gradient method is not yet implemented for this shape",0);
        }
    };

    /// A cube (3 dimensions)

    /// LIMIT: the 3 primary axes must be x, y, and z
    class SDFCube : public SDFBox {
    public:
        SDFCube(const double length, const coord_3d& center) 
            : SDFBox(length,center)
        {}
    };

    /// An ellipsoid (3 dimensions)
    
    /// LIMIT: the 3 primary axes must be x, y, and z
    class SDFEllipsoid : public ShapeDFInterface<3> {
    protected:
        coord_3d radii; ///< the directional radii
        coord_3d center; ///< the center
        
    public:
        SDFEllipsoid(const coord_3d& radii, const coord_3d& center) 
            : radii(radii)
            , center(center) 
        {}
        
        double sdf(const coord_3d& pt) const {
            double quot, sum;
            int i;
            
            sum = 0.0;
            for(i = 0; i < 3; ++i) {
                quot = (pt[i] - center[i]) / radii[i];
                sum += quot * quot;
            }
            
            return sum - 1.0;
        }

        coord_3d grad_sdf(const coord_3d& pt) const {
            MADNESS_EXCEPTION("gradient method is not yet implemented for this shape",0);
        }
    };
    
    /// A cylinder (3 dimensions)
    class SDFCylinder : public ShapeDFInterface<3> {
    protected:
        double radius; ///< the radius of the cylinder
        double a; ///< half the length of the cylinder
        coord_3d center; ///< the central axial point of the cylinder (distance a from both ends)
        coord_3d axis; ///< the axial direction of the cylinder

    public:
        SDFCylinder(const double radius, const double length, const coord_3d& axpt, const coord_3d& axis) 
            : radius(radius)
            , a(length / 2.0)
            , center(axpt)
            , axis(axis*(1.0/sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]))) 
        {}
        
        double sdf(const coord_3d& pt) const {
            double dist;
            coord_3d rel, radial;
            int i;
            
            // axial distance
            dist = 0.0;
            for(i = 0; i < 3; ++i) {
                rel[i] = pt[i] - center[i];
                dist += rel[i] * axis[i];
            }
            
            // get the radial component
            for(i = 0; i < 3; ++i)
                radial[i] = rel[i] - dist * axis[i];
            
            return std::max(fabs(dist) - a, sqrt(radial[0]*radial[0] + radial[1]*radial[1]
                                                 + radial[2]*radial[2]) - radius);
        }

        coord_3d grad_sdf(const coord_3d& pt) const {
            MADNESS_EXCEPTION("gradient method is not yet implemented for this shape",0);
        }
    };

} // end of madness namespace

#endif // MADNESS_MRA_SDF_SHAPE_3D_H__INCLUDED
