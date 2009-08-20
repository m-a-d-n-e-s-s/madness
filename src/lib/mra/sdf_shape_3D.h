/// This file defines SDF classes for common 3-D surfaces:
///
/// Plane
/// Sphere
/// Cone
/// Paraboloid
/// Box
/// Cube
/// Ellipsoid
/// Cylinder
///
/// NOTE: The signed distance functions would be the shortest distance between
/// a point and _any_ point on the surface.  This is hard to calculate in many
/// cases, so we use contours here.  The surface layer may not be equally thick
/// around all points on the surface.

#ifndef __madness_sdf_shape_3d__
#define __madness_sdf_shape_3d__

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/sdf_shape.h>

namespace madness {

/// A plane surface (3 dimensions)
template <typename Q>
class SDF_Plane : public SurfaceLayer<Q, 3> {
	protected:
		// the normal vector pointing OUTSIDE the surface
		Vector<Q, 3> normal;
		// a point on the surface
		Vector<Q, 3> point;

	public:
		/// Needs the width and threshold, the outward surface normal, and a
		/// point on the plane
		SDF_Plane(const Q width, const Q thresh,
			const Vector<Q, 3> &nrm, const Vector<Q, 3> &apoint)
			: SurfaceLayer<Q, 3>(width, thresh), normal(nrm), point(apoint) {

			// normalize the normal vector
			Q norm = sqrt(normal[0] * normal[0] + normal[1] * normal[1] +
				normal[2] * normal[2]);
			if(norm < 1.0e-12)
				MADNESS_EXCEPTION("Plane normal must have non-zero magnitude", 0);
			normal *= 1.0 / norm;
		}

		Q sdf(const Vector<Q, 3> &pt) const {
			// get the normal distance with dotp
			return (pt[0] - point[0]) * normal[0] + (pt[1] - point[1]) * normal[1]
				+ (pt[2] - point[2]) * normal[2];
		}
};

/// A spherical surface (3 dimensions)
template <typename Q>
class SDF_Sphere : public SurfaceLayer<Q, 3> {
	protected:
		// the radius
		Q radius;
		// the center
		Vector<Q, 3> center;

	public:
		/// Needs the width and threshold, the radius, and the center
		SDF_Sphere(const Q width, const Q thresh, const Q rad,
			const Vector<Q, 3> &cen) : SurfaceLayer<Q, 3>(width, thresh),
			radius(rad), center(cen) {}

		Q sdf(const Vector<Q, 3> &pt) const {
			Q temp, r;
			int i;

			r = 0.0;
			for(i = 0; i < 3; ++i) {
				temp = pt[i] - center[i];
				r += temp * temp;
			}

			return sqrt(r) - radius;
		}
};

/// A cone (3 dimensions): Sqrt(x^2 + y^2) - c * z == 0
template <typename Q>
class SDF_Cone : public SurfaceLayer<Q, 3> {
	protected:
		// the apex
		Vector<Q, 3> apex;
		// the radius
		Q c;
		// the direction of the axis, from the apex INSIDE
		Vector<Q, 3> dir;

	public:
		/// Needs the width and threshold, the radius, the apex point,
		/// and the direction in which the cone opens
		SDF_Cone(const Q width, const Q thresh, const Q radius,
			const Vector<Q, 3> &apex_pt, const Vector<Q, 3> &direc) :
			SurfaceLayer<Q, 3>(width, thresh), apex(apex_pt), c(radius),
			dir(direc) {

			// normalize the direction vector
			Q norm = 1.0 / sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
			dir *= norm;
		}

		Q sdf(const Vector<Q, 3> &pt) const {
			Vector<Q, 3> diff;
			unsigned int i;
			Q dotp;

			for(i = 0; i < 3; ++i)
				diff[i] = pt[i] - apex[i];

			dotp = diff[0]*dir[0] + diff[1]*dir[1] + diff[2]*dir[2];

			for(i = 0; i < 3; ++i)
				diff[i] -= dotp * dir[i];

			return sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2])
				- c * dotp;
		}
};

/// A paraboloid (3 dimensions): x^2 + y^2 - c * z == 0
template <typename Q>
class SDF_Paraboloid : public SurfaceLayer<Q, 3> {
	protected:
		// the apex
		Vector<Q, 3> apex;
		// the curvature
		Q c;
		// the direction of the axis, from the apex INSIDE
		Vector<Q, 3> dir;

	public:
		/// Needs the width and threshold, the curvature, the apex point,
		/// and the direction in which the paraboloid opens
		SDF_Paraboloid(const Q width, const Q thresh, const Q curve,
			const Vector<Q, 3> &apex_pt, const Vector<Q, 3> &direc) :
			SurfaceLayer<Q, 3>(width, thresh), apex(apex_pt), c(curve),
			dir(direc) {

			// normalize the direction vector
			Q norm = 1.0 / sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
			dir *= norm;
		}

		Q sdf(const Vector<Q, 3> &pt) const {
			Vector<Q, 3> diff;
			unsigned int i;
			Q dotp;

			for(i = 0; i < 3; ++i)
				diff[i] = pt[i] - apex[i];

			dotp = diff[0]*dir[0] + diff[1]*dir[1] + diff[2]*dir[2];

			for(i = 0; i < 3; ++i)
				diff[i] -= dotp * dir[i];

			return diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]
				- c * dotp;
		}
};

/// A box (3 dimensions)
///
/// LIMIT: the 3 primary axes must be x, y, and z
template <typename Q>
class SDF_Box : public SurfaceLayer<Q, 3> {
	protected:
		// half the length of each side of the box
		Vector<Q, 3> lengths;
		// the center of the box
		Vector<Q, 3> center;

	public:
		/// Needs teh width and threshold, the lengths of the sides, and the
		/// center
		SDF_Box(const Q width, const Q thresh, const Vector<Q, 3> &length2,
			const Vector<Q, 3> &cen) : SurfaceLayer<Q, 3>(width, thresh),
			lengths(length2), center(cen) {

			// only want half the length
			lengths *= 0.5;
		}

		Q sdf(const Vector<Q, 3> &pt) const {
			Q diff, max;
			int i;

			max = fabs(pt[0] - center[0]) - lengths[0];
			for(i = 1; i < 3; ++i) {
				diff = fabs(pt[i] - center[i]) - lengths[i];
				if(diff > max)
					max = diff;
			}

			return max;
		}
};

/// A cube (3 dimensions)
///
/// LIMIT: the 3 primary axes must be x, y, and z
template <typename Q>
class SDF_Cube : public SurfaceLayer<Q, 3> {
	protected:
		// half the length of one side of the cube
		Q a;
		// the center of the cube
		Vector<Q, 3> center;

	public:
		/// Needs the width and threshold, the length of a side, and the center
		SDF_Cube(const Q width, const Q thresh, const Q length2,
			const Vector<Q, 3> &cen) : SurfaceLayer<Q, 3>(width, thresh),
			a(length2 / 2.0), center(cen) {}

		Q sdf(const Vector<Q, 3> &pt) const {
			Q diff, max;
			int i;

			max = 0.0;
			for(i = 0; i < 3; ++i) {
				diff = fabs(pt[i] - center[i]);
				if(diff > max)
					max = diff;
			}

			return max - a;
		}
};

/// An ellipsoid (3 dimensions)
///
/// LIMIT: the 3 primary axes must be x, y, and z
template <typename Q>
class SDF_Ellipsoid : public SurfaceLayer<Q, 3> {
	protected:
		// the directional radii
		Vector<Q, 3> radii;
		// the center
		Vector<Q, 3> center;

	public:
		/// Needs the width and threshold, the directional radii, and the center
		SDF_Ellipsoid(const Q width, const Q thresh, const Vector<Q, 3> &dirrad,
			const Vector<Q, 3> &cen) : SurfaceLayer<Q, 3>(width, thresh),
			radii(dirrad), center(cen) {}

		Q sdf(const Vector<Q, 3> &pt) const {
			Q quot, sum;
			int i;

			sum = 0.0;
			for(i = 0; i < 3; ++i) {
				quot = (pt[i] - center[i]) / radii[i];
				sum += quot * quot;
			}

			return sum - 1.0;
		}
};

/// A cylinder (3 dimensions)
template <typename Q>
class SDF_Cylinder : public SurfaceLayer<Q, 3> {
	protected:
		// the radius of the cylinder
		Q radius;
		// half the length of the cylinder
		Q a;
		// the central axial point of the cylinder (distance a from both ends)
		Vector<Q, 3> center;
		// the axial direction of the cylinder
		Vector<Q, 3> axis;

	public:
		/// Needs the width and threshold, the radius of the cylinder, the
		/// axial length, the axial point at the center, and the axial direction
		SDF_Cylinder(const Q width, const Q thresh, const Q rad, const Q length,
			const Vector<Q, 3> &axpt, const Vector<Q, 3> &axdir) :
			SurfaceLayer<Q, 3>(width, thresh), radius(rad), a(length / 2.0),
			center(axpt), axis(axdir) {

			// normalize the direction axis
			Q norm = sqrt(axis[0] * axis[0] + axis[1] * axis[1] +
				axis[2] * axis[2]);
			if(norm < 1.0e-12)
				MADNESS_EXCEPTION("Cylinder axial direction must have non-zero magnitude", 0);
			axis *= 1.0 / norm;
		}

		Q sdf(const Vector<Q, 3> &pt) const {
			Q dist;
			Vector<Q, 3> rel, radial;
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

			return max(fabs(dist) - a, sqrt(radial[0]*radial[0] + radial[1]*radial[1]
				+ radial[2]*radial[2]) - radius);
		}
};

} // end of madness namespace

#endif
