#ifndef MADNESS_MOLOPT_H
#define MADNESS_MOLOPT_H

#include <madness/chem/molecule.h>
#include <madness/tensor/solvers.h>
#include <madness/tensor/tensor.h>
#include <madness/tensor/tensor_lapack.h>

#include <algorithm>
#include <string>

namespace madness {

// Generalized Optimization using targetT to provide energy and gradient of a
// molecule targetT must have the following methods:
//  void energy_and_gradient(Molecule& molecule, double& energy, Tensor<double>&
//  gradient) double value(const Tensor<double>& x)  // value of target at x
//  Let's express that in code so we can get compile time errors if we don't
//  have the right methods
//
//
//

class MolOpt {
 private:
  const int maxiter;     //< Maximum number of iterations
  const double maxstep;  //< Maximum step in any one coordinate (currently
                         //Cartesian in a.u.)
  const double etol;     //< Convergence test for energy change
  const double gtol;     //< Convergence test for maximum gradient element
  const double xtol;     //< Convergence test for Cartesian step in a.u.
  const double energy_precision;    //< Assumed precision in energy
  const double gradient_precision;  //< Assumed precision in the gradient
  const int print_level;     //< print_level=0 is none; 1 is default; 2 is debug
  const std::string update;  //< update = "bfgs" (default) or "sr1"

  Tensor<double> hessian;

  /// Returns new search direction given gradient and projected/shifted hessian
  Tensor<double> new_search_direction(const Tensor<double>& g,
                                      const Tensor<double>& h) const {
    Tensor<double> dx, s;
    double tol = gradient_precision;  // threshold for small hessian eigenvalues
    double trust = std::min(maxstep * g.dim(0), 1.0);

    Tensor<double> v, e;
    syev(h, v, e);
    if (print_level > 1) print("hessian eigenvalues", e);

    // Transform gradient into spectral basis
    Tensor<double> gv = inner(g, v);
    if (print_level > 1) print("spectral gradient", gv);

    // Take step applying restrictions
    int nneg = 0, nsmall = 0, nrestrict = 0;
    for (int i = 0; i < e.dim(0); i++) {
      if (e[i] > 900.0) {
        // This must be a translation or rotation ... skip it
        if (print_level > 1) print("skipping redundant mode", i);
      } else if (e[i] <
                 -tol) {  // BGFS hessian should be positive ... SR1 may not be
        if (print_level > 0)
          printf("   forcing negative eigenvalue to be positive %d %.1e\n", i,
                 e[i]);
        nneg++;
        e[i] = tol;  // or -e[i] ??
      } else if (e[i] < tol) {
        if (print_level > 0)
          printf("   forcing small eigenvalue to be positive %d %.1e\n", i,
                 e[i]);
        nsmall++;
        e[i] = tol;
      }

      gv[i] = -gv[i] / e[i];  // Newton step

      if (std::abs(gv[i]) > trust) {  // Step restriction
        double gvnew = trust * std::abs(gv(i)) / gv[i];
        if (print_level > 0)
          printf("   restricting step in spectral direction %d %.1e --> %.1e\n",
                 i, gv[i], gvnew);
        nrestrict++;
        gv[i] = gvnew;
      }
    }
    if (print_level > 0 && (nneg || nsmall || nrestrict))
      printf("   nneg=%d nsmall=%d nrestrict=%d\n", nneg, nsmall, nrestrict);

    // Transform back from spectral basis
    gv = inner(v, gv);

    if (print_level > 1) print("cartesian dx before restriction", gv);

    // Now apply step restriction in real space
    bool printing = false;
    for (int i = 0; i < gv.dim(0); i++) {
      if (fabs(gv[i]) > maxstep) {
        gv[i] = maxstep * gv[i] / fabs(gv[i]);
        if (print_level > 0) {
          if (!printing) printf("  restricting step in Cartesian direction");
          printing = true;
          printf(" %d", i);
        }
      }
    }
    if (printing) printf("\n");

    return gv;
  }

  /// Makes the projector onto independent coordinates

  /// For Cartesians \code P*x removes the rotations and translations;
  /// eventually will add support for redundant internal coordinates
  Tensor<double> make_projector(const Molecule& molecule) {
    const int natom = molecule.natom();
    const Tensor<double> coords = molecule.get_all_coords();  // (natom,3)

    // Construct normalized vectors in V in the direction of the translations
    // and infinitesimal rotations
    Tensor<double> V(6, natom, 3);  // First 3 translations, second 3 rotations

    for (int k = 0; k < 3; k++)  // Translations already orthonormal
      V(k, _, k) = 1.0 / std::sqrt(static_cast<double>(natom));

    Tensor<double> centroid(3);
    for (int k = 0; k < 3; k++) centroid(k) = coords(_, k).sum() / natom;
    if (print_level > 1) print("centroid", centroid);

    for (int i = 0; i < natom; i++) {
      double x = coords(i, 0) - centroid[0];
      double y = coords(i, 1) - centroid[1];
      double z = coords(i, 2) - centroid[2];

      V(3, i, 0) = 0;  // Rotn about x axis
      V(3, i, 1) = z;
      V(3, i, 2) = -y;

      V(4, i, 0) = z;  // Rotn about y axis
      V(4, i, 1) = 0;
      V(4, i, 2) = -x;

      V(5, i, 0) = -y;
      V(5, i, 1) = x;
      V(5, i, 2) = 0;  // Rotn about z axis
    }

    V = V.reshape(6, 3 * natom);

    if (print_level > 1) print("V before orthonormal");
    if (print_level > 1) print(V);

    // Normalize rotations, orthonormalize rotns and translations,
    // noting may end up with a zero vector for linear molecules
    for (int i = 3; i < 6; i++) {
      V(i, _).scale(1.0 / V(i, _).normf());
      for (int j = 0; j < i; j++) {
        double s = V(i, _).trace(V(j, _));
        V(i, _) -= V(j, _) * s;
      }
      double vnorm = V(i, _).normf();
      if (vnorm > 1e-6) {
        V(i, _) *= 1.0 / vnorm;
      } else {
        V(i, _) = 0.0;
      }
    }

    // The projector is 1 - VT*V
    V = -inner(transpose(V), V);
    for (int i = 0; i < 3 * natom; i++) V(i, i) += 1.0;

    return V;
  }

  /// a1 is initial step (usually pick 1)
  /// energy0 is energy at zero step
  /// dxgrad is gradient projected onto search dir
  ///
  template <typename targetT>
  double line_search(Molecule& molecule, targetT& target,
                     const Tensor<double>& dx, double energy0, double dxgrad,
                     double a1 = 1.0) {
    double energy1;
    double hess, a2;
    const char* lsmode = "";

    Tensor<double> x = molecule.get_all_coords().flat();

    // Ensure we are walking downhill (BFGS should ensure that, but SR1 may not)
    if (dxgrad * a1 > 0.0) {
      if (print_level > 0) print("    line search gradient +ve ", a1, dxgrad);
      a1 = -a1;
    }

    // Compute energy at new point
    energy1 = target.value(x + a1 * dx);

    // Fit to a parabola using energy0, g0, energy1
    hess = 2.0 * (energy1 - energy0 - a1 * dxgrad) / (a1 * a1);
    a2 = -dxgrad / hess;  // Newton step

    if (std::abs(energy1 - energy0) <
        energy_precision) {  // Insufficient precision
      a2 = a1;
      lsmode = "fixed";
    } else if (hess > 0.0) {                           // Positive curvature
      if ((energy1 - energy0) <= -energy_precision) {  // a1 step went downhill
        lsmode = "downhill";
        if (std::abs(a2) >
            4.0 * std::abs(a1)) {  // Walking down hill but don't go too far
          lsmode = "restrict";
          a2 = 4.0 * a1;
        }
      } else {  // a1 step went uphill ... we have bracketed the minimum.
        lsmode = "bracket";
      }
    } else {                                         // Negative curvature
      if ((energy1 - energy0) < energy_precision) {  // keep walking down hill
        lsmode = "negative";
        a2 = 2e0 * a1;
      } else {
        lsmode = "punt";  // negative curvature but no apparent progress
        a2 = a1;
      }
    }

    if (std::abs(a2 - a1) <
        0.2 * std::abs(a1)) {  // Take full step to avoid reconverging SCF
      a2 = a1;
      lsmode = "fixed2";
    }

    // Predicted next energy
    double energy2 = energy0 + dxgrad * a2 + 0.5 * hess * a2 * a2;

    if (print_level > 0) {
      printf("\n   line search grad=%.2e hess=%.2e mode=%s newstep=%.3f\n",
             dxgrad, hess, lsmode, a2);
      printf("                      predicted %.12e\n\n", energy2);
    }

    return a2;
  }

 public:
  MolOpt(int maxiter = 20, double maxstep = 0.1, double etol = 1e-4,
         double gtol = 1e-3, double xtol = 1e-3, double energy_precision = 1e-5,
         double gradient_precision = 1e-4, int print_level = 1,
         std::string update = "BFGS")
      : maxiter(maxiter),
        maxstep(maxstep),
        etol(std::max(etol, energy_precision)),
        gtol(std::max(gtol, gradient_precision)),
        xtol(xtol),
        energy_precision(energy_precision),
        gradient_precision(gradient_precision),
        print_level(print_level),
        update(update)

  {
    if (print_level > 0) {
      std::cout << endl;
      print_justified("Molecular Geometry Optimization");
      std::cout << endl;
      print("       maximum iterations", maxiter);
      print("             maximum step", maxstep);
      print("       energy convergence", etol);
      print("     gradient convergence", gtol);
      print("    cartesian convergence", xtol);
      print("         energy precision", energy_precision);
      print("       gradient precision", gradient_precision);
      print("           hessian update", update);
    }
  }

  void set_hessian(const Tensor<double>& h) { hessian = h; }

  const Tensor<double>& get_hessian() const { return hessian; }

  void initialize_hessian(const Molecule& molecule) {
    const int N = 3 * molecule.natom();
    hessian = Tensor<double>(N, N);
    for (int i = 0; i < N; i++) hessian(i, i) = 0.5;
  }

  template <typename targetT>
  Molecule optimize(Molecule molecule,
                    targetT& target) {  ////!!!!!!! pass by value
    const int natom = molecule.natom();

    // Code structured so it will be straightforward to introduce redundant
    // internal coordinates

    if (hessian.size() == 0) initialize_hessian(molecule);

    double ep = 0.0;               // Previous energy
    Tensor<double> gp(3 * natom);  // Previous gradient
    Tensor<double> dx(3 * natom);  // Current search direction
    gp = 0.0;
    dx = 0.0;

    for (int iter = 0; iter < maxiter; iter++) {
      if (print_level > 0)
        print("\n\n Geometry optimization iteration", iter, "\n");
      if (print_level > 0) molecule.print();

      double e;
      Tensor<double> g;

      target.energy_and_gradient(molecule, e, g);

      double de = e - ep;
      double dxmax = dx.absmax();
      double gmax = g.absmax();

      bool dxconv = (iter > 0) && (dxmax < xtol);
      bool gconv = gmax < gtol;
      bool econv = (iter > 0) && (std::abs(de) < etol);
      bool converged = econv && dxconv && gconv;

      if (!converged && gmax < gradient_precision) {
        if (print_level > 0)
          print(
              "\nInsufficient precision in gradient to proceed further -- "
              "forcing convergence\n");
        converged = true;
      }

      if (print_level > 0) {
        const char* tf[] = {"F", "T"};
        print(" ");
        printf(
            "      energy        delta-e     max-dx     max-g     e  dx   g\n");
        printf(
            " ----------------  ---------  ---------  ---------  --- --- "
            "---\n");
        printf(" %15.6f   %9.2e  %9.2e  %9.2e   %s   %s   %s\n", e, de, dxmax,
               gmax, tf[econv], tf[dxconv], tf[gconv]);
        // print(e, de, econv, dxmax, dxconv, dxnorm, gmax, gconv, gnorm,
        // converged);
        print(" ");
      }

      if (converged) {
        if (print_level > 0) print("\n Geometry optimization converged!\n");
        if (print_level > 0) molecule.print();
        break;
      }

      // Construct projector
      Tensor<double> P = make_projector(molecule);

      // Project the gradient before updating Hessian
      g = inner(P, g);
      if (print_level > 1) print("gradient after projection", g);

      if (iter > 0) {
        if ((g - gp).absmax() < 2.0 * gradient_precision) {
          if (print_level > 0)
            print(
                "  skipping hessian update due to insufficient precision in "
                "gradient");
        } else if (update == "bfgs") {
          QuasiNewton::hessian_update_bfgs(dx, g - gp, hessian);
        } else if (update == "sr1") {
          QuasiNewton::hessian_update_sr1(dx, g - gp, hessian);
        } else {
          throw "unknown update";
        }
      }

      ep = e;
      gp = g;

      // Construct the projected and shifted hessian = PHP + shift*(1-P)
      const double shift = 1000.0;  // this value assumed in new_search_dir
      Tensor<double> PHPS = inner(P, inner(hessian, P)) - shift * P;
      if (print_level > 1) {
        print("projector");
        print(P);
      }
      for (int i = 0; i < 3 * natom; i++) PHPS(i, i) += shift;
      if (print_level > 1) {
        print("PHPS");
        print(PHPS);
      }

      // Construct new search direction by diagonalizing Hessian and taking
      // spectral step
      dx = new_search_direction(g, PHPS);
      if (print_level > 1) print("dx", dx);

      // Line search
      double alpha = line_search(molecule, target, dx, e, dx.trace(g), 1.0);
      if (print_level > 1) print("step", alpha);
      dx.scale(alpha);
      if (print_level > 1) print("scaled dx", dx);

      // Take the step
      Tensor<double> x = molecule.get_all_coords().flat();
      x += dx;
      molecule.set_all_coords(x.reshape(natom, 3));

      if (print_level > 1) print("new molecular coords");
    }
    // return the optimized molecule
    return molecule;
  }
};
}  // namespace madness
#endif  // MADNESS_MOLOPT_H
