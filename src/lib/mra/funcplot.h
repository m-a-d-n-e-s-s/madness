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
#ifndef MADNESS_MRA_FUNCPLOT_H__INCLUDED
#define MADNESS_MRA_FUNCPLOT_H__INCLUDED

#include <constants.h>

/*!

  \file mra/funcplot.h
  \brief Defines/implements plotting interface for functions
  \ingroup funcplot

  @{
 */

namespace madness {
    /// Writes an OpenDX format file with a cube/slice of points on a uniform grid

    /// Collective operation but only process 0 writes the file.  By convention OpenDX
    /// files end in ".dx" but this choice is up to the user.  The binary format is
    /// more compact and vastly faster to both write and load but is not as portable.
    ///
    /// Now follow some brief tips about how to look at files inside OpenDX.
    ///
    /// To view a 1D function \c file-selector-->import-->plot-->image.
    ///
    /// To view a 2D function as a colored plane \c file-selector-->import-->autocolor-->image.
    ///
    /// To view a 2D function as a 3D surface \c file-selector-->import-->rubbersheet-->image.
    ///
    /// To view a 3D function as an isosurface \c file-selector-->import-->isosurface-->image.
    ///
    /// To select the real/imaginary/absolute value of a complex number insert a compute
    /// element after the import.
    template <typename T, std::size_t NDIM>
    void plotdx(const Function<T,NDIM>& f,
                const char* filename,
                const Tensor<double>& cell = FunctionDefaults<NDIM>::get_cell(),
                const std::vector<long>& npt = std::vector<long>(NDIM,201L),
                bool binary=true);


    /// Writes the header information of a VTK file for plotting in an external
    /// post-processing package (such as Paraview)
    //
    /// @param world World communicator
    /// @param filename String containing the filename to export to
    /// @param plotlo Vector of double values indicating the minimum coordinate to plot to in each dimension
    /// @param plothi Vector of double values indicating the maximum coordinate to plot to in each dimension
    /// @param npt Vector of long integers indicating the number of points to plot in each dimension
    /// @param binary (optional) Boolean indicating whether to print in binary

    /// The VTK routines are also designed for SERIAL data, parallel coming...
    ///
    /// This header is templated by the dimension of the data.
    ///
    /// To plot with the plotvtk_* routines:
    ///    plotvtk_begin(...)
    ///    plotvtk_data(...)
    ///    plotvtk_data(...) ...
    ///    plotvtk_end(...)
    ///
    /// NOTE: Paraview expects the structured mesh points in a particular
    /// order, which is why the LowDimIndexIterator is used...
    template<std::size_t NDIM>
    void plotvtk_begin(World &world, const char *filename,
        const Vector<double, NDIM> &plotlo, const Vector<double, NDIM> &plothi,
        const Vector<long, NDIM> &npt, bool binary = false) {

        PROFILE_FUNC;
        MADNESS_ASSERT(NDIM>=1 && NDIM<=3); // how do we plot data in more than 3-D?

        Tensor<double> cell(NDIM, 2);
        std::size_t i;
        for(i = 0; i < NDIM; ++i) {
            cell(i, 0) = plotlo[i];
            cell(i, 1) = plothi[i];
        }

        FILE *f=0;
        if(world.rank() == 0) {
            f = fopen(filename, "w");
            if(!f)
                MADNESS_EXCEPTION("plotvtk: failed to open the plot file", 0);

            fprintf(f, "<VTKFile type=\"StructuredGrid\" version=\"0.1\"" \
                " byte_order=\"LittleEndian\" compressor=\"" \
                "vtkZLibDataCompressor\">\n");
            fprintf(f, "  <StructuredGrid WholeExtent=\"");
            for(i = 0; i < NDIM; ++i)
                fprintf(f, "0 %ld ", npt[i]-1);
            for(; i < 3; ++i)
                fprintf(f, "0 0 ");
            fprintf(f, "\">\n");
            fprintf(f, "    <Piece Extent=\"");
            for(i = 0; i < NDIM; ++i)
                fprintf(f, "0 %ld ", npt[i]-1);
            for(; i < 3; ++i)
                fprintf(f, "0 0 ");
            fprintf(f, "\">\n");
            fprintf(f, "      <Points>\n");
            fprintf(f, "        <DataArray NumberOfComponents=\"3\" " \
                "type=\"Float32\" format=\"ascii\">\n");

            Vector<double, NDIM> space;
            for(i = 0; i < NDIM; ++i) {
                if(npt[i] == 1)
                    space[i] = 0.0;
                else
                    space[i] = (cell(i, 1) - cell(i, 0)) / (npt[i] - 1);
            }

            // go through the grid
            for(LowDimIndexIterator it(npt); it; ++it) {
                for(i = 0; i < NDIM; ++i)
                    fprintf(f, "%f ", plotlo[i] + it[i]*space[i]);
                for(; i < 3; ++i)
                    fprintf(f, "0.0 ");
                fprintf(f, "\n");
            }

            fprintf(f, "        </DataArray>\n");
            fprintf(f, "      </Points>\n");
            fprintf(f, "      <PointData>\n");
            fclose(f);
        }
        world.gop.fence();
    }

    /// Generic VTK data writer. Specific type instances of this function are defined for
    /// both real and complex valued functions.
    //
    /// @param function Function (real or complex) that we wish to export the data of
    /// @param fieldname A string containing the name we wish to refer to this field as in the exported data
    /// @param world World communicator
    /// @param filename String containing the filename to export to
    /// @param plotlo Vector of double values indicating the minimum coordinate to plot to in each dimension
    /// @param plothi Vector of double values indicating the maximum coordinate to plot to in each dimension
    /// @param npt Vector of long integers indicating the number of points to plot in each dimension
    /// @param binary (optional) Boolean indicating whether to print in binary

    /// This templated function won't do anything except print a warning
    /// message.  Specialized versions of this function should be used.
    template<typename T, std::size_t NDIM>
    void plotvtk_data(const T &function, const char *fieldname, World &world,
        const char *filename, const Vector<double, NDIM> &plotlo,
        const Vector<double, NDIM> &plothi, const Vector<long, NDIM> &npt,
        bool binary = false) {

        MADNESS_EXCEPTION("plotvtk only supports madness::functions", 0);
    }

    /// VTK data writer for real-valued (not complex) madness::functions.

    /// Set plot_refine=true to get a plot of the refinement levels of
    /// the given function.
    template<typename T, std::size_t NDIM>
    void plotvtk_data(const Function<T, NDIM> &function, const char *fieldname,
        World &world, const char *filename, const Vector<double, NDIM> &plotlo,
        const Vector<double, NDIM> &plothi, const Vector<long, NDIM> &npt,
        bool binary = false, bool plot_refine = false) {

        PROFILE_FUNC;
        MADNESS_ASSERT(NDIM>=1 && NDIM<=3); // no plotting high-D functions, yet...

        Tensor<double> cell(NDIM, 2);
        std::size_t i;
        for(i = 0; i < NDIM; ++i) {
            cell(i, 0) = plotlo[i];
            cell(i, 1) = plothi[i];
        }
        std::vector<long> numpt(NDIM);
        for(i = 0; i < NDIM; ++i)
            numpt[i] = npt[i];

        world.gop.barrier();

        function.verify();
        FILE *f = 0;
        if(world.rank() == 0) {
            f = fopen(filename, "a");
            if(!f)
                MADNESS_EXCEPTION("plotvtk: failed to open the plot file", 0);

            fprintf(f, "        <DataArray Name=\"%s\" format=\"ascii\" " \
                "type=\"Float32\" NumberOfComponents=\"1\">\n", fieldname);
        }

        world.gop.fence();
        Tensor<T> tmpr = function.eval_cube(cell, numpt, plot_refine);
        world.gop.fence();

        if(world.rank() == 0) {
            for(LowDimIndexIterator it(numpt); it; ++it) {
                fprintf(f, "%.6e\n", tmpr(*it));
            }
            fprintf(f, "        </DataArray>\n");
            fclose(f);
        }
        world.gop.fence();
    }

    /// VTK data writer for complex-valued madness::functions.

    /// The complex-value is written as two reals (a vector from VTK's
    /// perspective.  The first (X) component is the real part and the second
    /// (Y) component is the imaginary part.
    /// Set plot_refine=true to get a plot of the refinement levels of
    /// the given function.
    template<typename T, std::size_t NDIM>
    void plotvtk_data(const Function<std::complex<T>, NDIM> &function,
        const char *fieldname, World &world, const char *filename,
        const Vector<double, NDIM> &plotlo, const Vector<double, NDIM> &plothi,
        const Vector<long, NDIM> &npt, bool binary = false,
        bool plot_refine = false) {

        // this is the same as plotvtk_data for real functions, except the
        // real and imaginary parts are printed on the same line (needed
        // to change NumberOfComponents in the XML tag)

        PROFILE_FUNC;
        MADNESS_ASSERT(NDIM>=1 && NDIM<=3); // no plotting high-D functions, yet...

        Tensor<double> cell(NDIM, 2);
        std::size_t i;
        for(i = 0; i < NDIM; ++i) {
            cell(i, 0) = plotlo[i];
            cell(i, 1) = plothi[i];
        }
        std::vector<long> numpt(NDIM);
        for(i = 0; i < NDIM; ++i)
            numpt[i] = npt[i];

        world.gop.barrier();

        function.verify();
        FILE *f = 0;
        if(world.rank() == 0) {
            f = fopen(filename, "a");
            if(!f)
                MADNESS_EXCEPTION("plotvtk: failed to open the plot file", 0);

            fprintf(f, "        <DataArray Name=\"%s\" format=\"ascii\" " \
                "type=\"Float32\" NumberOfComponents=\"2\">\n", fieldname);
        }

        world.gop.fence();
        Tensor<std::complex<T> > tmpr = function.eval_cube(cell, numpt,
                                                           plot_refine);
        world.gop.fence();

        if(world.rank() == 0) {
            for(LowDimIndexIterator it(numpt); it; ++it) {
                fprintf(f, "%.6e %.6e\n", real(tmpr(*it)), imag(tmpr(*it)));
            }
            fprintf(f, "        </DataArray>\n");
            fclose(f);
        }
        world.gop.fence();
    }

    /// Writes the footer information of a VTK file for plotting in an external
    /// post-processing package (such as Paraview)
    //
    /// @param world World communicator
    /// @param filename Name of VTK file
    /// @param binary (Optional) Boolean indicating whether to print in binary
    template<std::size_t NDIM>
    void plotvtk_end(World &world, const char *filename, bool binary = false) {
        PROFILE_FUNC;
        MADNESS_ASSERT(NDIM>=1 && NDIM<=3);

        FILE *f = 0;
        if(world.rank() == 0) {
            f = fopen(filename, "a");
            if(!f)
                MADNESS_EXCEPTION("plotvtk: failed to open the plot file", 0);

            fprintf(f, "      </PointData>\n");
            fprintf(f, "      <CellData>\n");
            fprintf(f, "      </CellData>\n");
            fprintf(f, "    </Piece>\n");
            fprintf(f, "  </StructuredGrid>\n");
            fprintf(f, "</VTKFile>\n");
            fclose(f);
        }
        world.gop.fence();
    }

    namespace detail {
        inline unsigned short htons_x(unsigned short a) {
            return (a>>8) | (a<<8);
        }
    }

    /// Writes a Povray DF3 format file with a cube of points on a uniform grid

    /// Collective operation but only process 0 writes the file.  By convention Povray
    /// files end in ".df3" but this choice is up to the user.  The dynamic range of
    /// function values is mapped onto [0,1] and values stored in 16-bit fixed precision.
    template <typename T>
    static void plotpovray(const Function<T,3>& function,
                           const char* filename,
                           const Tensor<double>& cell = FunctionDefaults<3>::get_cell(),
                           const std::vector<long>& npt = std::vector<long>(3,201L))
    {
        using detail::htons_x;

        MADNESS_ASSERT(npt.size() == 3);
        unsigned short dims[3] = {htons_x(npt[0]),htons_x(npt[1]),htons_x(npt[2])};

        World& world = const_cast< Function<T,3>& >(function).world();
        FILE *f=0;
        if (world.rank() == 0) {
            f = fopen(filename, "w");
            if (!f) MADNESS_EXCEPTION("plotdx: failed to open the plot file", 0);
            fwrite((void*) dims, sizeof(short), 3, f);
        }
        Tensor<T> r = function.eval_cube(cell, npt);
        if (world.rank() == 0) {
            double rmax = r.max();
            double rmin = r.min();
            double rrange = rmax + rmin;
            double rmean = rrange*0.5;
            double fac = 65535.0/rrange;

            printf("plot_povray: %s: min=%.2e(0.0) mean=%.2e(0.5) max=%.2e(1.0) range=%.2e\n",
                   filename,rmin,rmean,rmax,rrange);

            std::vector<unsigned short> d(npt[0]);
            for (unsigned int i2=0; i2<npt[2]; ++i2) {
                for (unsigned int i1=0; i1<npt[1]; ++i1) {
                    for (unsigned int i0=0; i0<npt[0]; ++i0) {
                        d[i0] = (unsigned short)(htons_x((unsigned short)(fac*(r(i0,i1,i2) - rmin))));
                        //printf("%d\n",htons_x(d[i0]));
                    }
                    fwrite((void*) &d[0], sizeof(short), npt[0], f);
                }
            }

            fclose(f);
        }
    }

    static inline void plot_line_print_value(FILE* f, double_complex v) {
        fprintf(f, "    %.14e %.14e   ", real(v), imag(v));
    }

    static inline void plot_line_print_value(FILE* f, double v) {
        fprintf(f, " %.14e", v);
    }

    /// Generates ASCII file tabulating f(r) at npoints along line r=lo,...,hi

    /// The ordinate is distance from lo
    template <typename T, std::size_t NDIM>
    void plot_line(const char* filename, int npt, const Vector<double,NDIM>& lo, const Vector<double,NDIM>& hi,
                   const Function<T,NDIM>& f) {
        typedef Vector<double,NDIM> coordT;
        coordT h = (hi - lo)*(1.0/(npt-1));

        double sum = 0.0;
        for (std::size_t i=0; i<NDIM; ++i) sum += h[i]*h[i];
        sum = sqrt(sum);

        World& world = f.world();
        f.reconstruct();
        if (world.rank() == 0) {
            FILE* file = fopen(filename,"w");
	    if(!file)
	      MADNESS_EXCEPTION("plot_line: failed to open the plot file", 0);
            for (int i=0; i<npt; ++i) {
                coordT r = lo + h*double(i);
                fprintf(file, "%.14e ", i*sum);
                plot_line_print_value(file, f.eval(r));
                fprintf(file,"\n");
            }
            fclose(file);
        }
        world.gop.fence();
    }

    /// Generates ASCII file tabulating f(r) and g(r) at npoints along line r=lo,...,hi

    /// The ordinate is distance from lo
    template <typename T, typename U, std::size_t NDIM>
    void plot_line(const char* filename, int npt, const Vector<double,NDIM>& lo, const Vector<double,NDIM>& hi,
                   const Function<T,NDIM>& f, const Function<U,NDIM>& g) {
        typedef Vector<double,NDIM> coordT;
        coordT h = (hi - lo)*(1.0/(npt-1));

        double sum = 0.0;
        for (std::size_t i=0; i<NDIM; ++i) sum += h[i]*h[i];
        sum = sqrt(sum);

        World& world = f.world();
        f.reconstruct();
        g.reconstruct();
        if (world.rank() == 0) {
            FILE* file = fopen(filename,"w");
	    if(!file)
	      MADNESS_EXCEPTION("plot_line: failed to open the plot file", 0);
            for (int i=0; i<npt; ++i) {
                coordT r = lo + h*double(i);
                fprintf(file, "%.14e ", i*sum);
                plot_line_print_value(file, f.eval(r));
                plot_line_print_value(file, g.eval(r));
                fprintf(file,"\n");
            }
            fclose(file);
        }
        world.gop.fence();
    }


    /// Generates ASCII file tabulating f(r), g(r), and a(r) at npoints along line r=lo,...,hi

    /// The ordinate is distance from lo
    template <typename T, typename U, typename V, std::size_t NDIM>
    void plot_line(const char* filename, int npt, const Vector<double,NDIM>& lo, const Vector<double,NDIM>& hi,
                   const Function<T,NDIM>& f, const Function<U,NDIM>& g, const Function<V,NDIM>& a) {
        typedef Vector<double,NDIM> coordT;
        coordT h = (hi - lo)*(1.0/(npt-1));

        double sum = 0.0;
        for (std::size_t i=0; i<NDIM; ++i) sum += h[i]*h[i];
        sum = sqrt(sum);

        World& world = f.world();
        f.reconstruct();
        g.reconstruct();
        a.reconstruct();
        if (world.rank() == 0) {
            FILE* file = fopen(filename,"w");
	    if(!file)
	      MADNESS_EXCEPTION("plot_line: failed to open the plot file", 0);
            for (int i=0; i<npt; ++i) {
                coordT r = lo + h*double(i);
                fprintf(file, "%.14e ", i*sum);
                plot_line_print_value(file, f.eval(r));
                plot_line_print_value(file, g.eval(r));
                plot_line_print_value(file, a.eval(r));
                fprintf(file,"\n");
            }
            fclose(file);
        }
        world.gop.fence();
    }

    /// Generates ASCII file tabulating f(r), g(r), a(r), b(r) at npoints along line r=lo,...,hi

    /// The ordinate is distance from lo
    template <typename T, typename U, typename V, typename W, std::size_t NDIM>
    void plot_line(const char* filename, int npt, const Vector<double,NDIM>& lo, const Vector<double,NDIM>& hi,
                   const Function<T,NDIM>& f, const Function<U,NDIM>& g, const Function<V,NDIM>& a, const Function<W,NDIM>& b) {
        typedef Vector<double,NDIM> coordT;
        coordT h = (hi - lo)*(1.0/(npt-1));

        double sum = 0.0;
        for (std::size_t i=0; i<NDIM; ++i) sum += h[i]*h[i];
        sum = sqrt(sum);

        World& world = f.world();
        f.reconstruct();
        g.reconstruct();
        a.reconstruct();
        b.reconstruct();
        if (world.rank() == 0) {
            FILE* file = fopen(filename,"w");
            for (int i=0; i<npt; ++i) {
                coordT r = lo + h*double(i);
                fprintf(file, "%.14e ", i*sum);
                plot_line_print_value(file, f.eval(r));
                plot_line_print_value(file, g.eval(r));
                plot_line_print_value(file, a.eval(r));
                plot_line_print_value(file, b.eval(r));
                fprintf(file,"\n");
            }
            fclose(file);
        }
        world.gop.fence();
    }


    template<typename T>
    static std::string stringify(T arg) {
    	std::ostringstream o;
    	if (!(o << arg))
    		MADNESS_EXCEPTION("stringify<T> failed",1);
    	return o.str();
    }


    typedef Vector<double,3> coord_3d;
    typedef Vector<double,6> coord_6d;

    // plot along this trajectory
    template<size_t NDIM>
    struct trajectory {

        typedef Vector<double,NDIM> coordT;

    	double lo;
    	double hi;
    	double radius;
    	long npt;
    	coordT start, end;
    	coord_3d el2;
    	coordT (*curve)(const coordT& lo, const coordT& hi, double radius, coord_3d el2, long npt, long ipt);

    	/// some tools for plotting MRA ranks of low order tensors

    	// return a hue number [0,0.7] according to the rank in relation to maxrank,
    	static double hueCode(const int rank) {
    			const double maxrank=40.0;
    			double hue=0.7-(0.7/maxrank)*(rank);
    			return std::max(0.0,hue);
    	}


        // print a dot of hue color at (x,y) in file f
        static void print_psdot(FILE *f, double x, double y, double color) {
        	fprintf(f,"\\newhsbcolor{mycolor}{%8.4f 1.0 0.7}\n",color);
            fprintf(f,"\\psdot[linecolor=mycolor](%12.8f,%12.8f)\n",x,y);
        }


        static coord_3d circle2(double lo, double hi, double radius, coord_3d el2, long npt, long ipt) {
        	double stepsize=constants::pi * 2.0 / npt;
        	double phi=ipt*stepsize;

        	// in the xz plane
        	coord_3d coord(0.0);
        	coord[0]=radius * sin(phi);
        	coord[1]=radius * cos(phi);
        	return coord;

        }

        static coord_6d circle_6d(const coord_6d& lo, const coord_6d& hi, double radius, coord_3d el2, long npt, long ipt) {
        	double stepsize=constants::pi * 2.0 / npt;

        	// start at phi=1.0
        	double phi=1.0+constants::pi+ipt*stepsize;

        	// in the xz plane
        	coord_6d coord(0.0);
        	coord[0]=radius * sin(phi);
        	coord[1]=radius * cos(phi);
        	coord[2]=0.0;
        	coord[3]=el2[0];
        	coord[4]=el2[1];
        	coord[5]=el2[2];

        	return coord;

        }


    //	typedef Vector<double,NDIM> (trajectory::circle_6d)(double lo, double hi, double radius, long npt, long ipt) const;

    	trajectory() {}
//    	// ctor for a straight line thru the origin
//    	trajectory(double lo, double hi, long npt) : lo(lo), hi(hi), npt(npt), curve(line) {
//    	}

    	// ctor for circle
    	trajectory(double radius, long npt) : radius(radius), npt(npt), curve(this->circle2) {
    	}

    	// ctor for circle with electron 2 fixed at coord_3d
    	trajectory(double radius, coord_3d el2, long npt) : radius(radius), npt(npt), el2(el2), curve(this->circle_6d) {
    	}


        static Vector<double, NDIM> line_internal(const coordT& lo, const coordT& hi, double radius, coord_3d el2, long npt, long ipt) {
            const coordT step=(hi-lo)*(1.0/npt);
            coordT coord=lo+step*ipt;
            return coord;
        }


        /// constructor for a line
//        static trajectory line(const Vector<double,NDIM>& lo, const Vector<double,NDIM>& hi, const long npt) {
        static trajectory line2(const coordT start, const coordT end, const long npt) {
            trajectory<NDIM> traj;
            traj.start=start;
            traj.end=end;
            traj.npt=npt;
            traj.curve=(trajectory::line_internal);
            return traj;
        }

        /// EZ ctor for a line a direction xyz={0,1,2,..,NDIM-1} thru the origin
        static trajectory line_xyz(const int xyz, const long npt) {
            double L=FunctionDefaults<NDIM>::get_cell_width()[0];
            coordT lo(0.0), hi(0.0);
            lo[xyz]=-L/2;
            hi[xyz]=L/2;
            return trajectory<NDIM>::line2(lo,hi,npt);
        }

        Vector<double,NDIM> operator()(int ipt) {
            return curve(start,end,radius,el2,npt,ipt);
        }

    };



    // plot along a line
    template<size_t NDIM>
    void plot_along(World& world, trajectory<NDIM> traj, const Function<double,NDIM>& function, std::string filename) {
    	 FILE *f=0;
    	 const int npt=traj.npt;

    	 const bool psdot=false;

    	 if(world.rank() == 0) {
    		 f = fopen(filename.c_str(), "w");
    		 if(!f) MADNESS_EXCEPTION("plot_along: failed to open the plot file", 0);

    		 if (psdot) {

                 fprintf(f,"\\psset{xunit=0.1cm}\n");
                 fprintf(f,"\\psset{yunit=10cm}\n");
                 fprintf(f,"\\begin{pspicture}(0,-0.3)(100,1.0)\n");
                 fprintf(f,"\\pslinewidth=0.05pt\n");
    		 }

    		 // walk along the line
    		 for (int ipt=0; ipt<npt; ipt++) {
    			 Vector<double,NDIM> coord=traj(ipt);
    			 if (psdot) {
    			     long rank=function.evalR(coord);
    			     trajectory<NDIM>::print_psdot(f,ipt,function(coord),trajectory<NDIM>::hueCode(rank));
    			 } else {
    			     fprintf(f,"%4i %12.6f\n",ipt, function(coord));
    			 }
    		 }


             if (psdot) fprintf(f,"\\end{pspicture}\n");

    		 fclose(f);
    	 }
    	 world.gop.fence();
    	 ;
    }


    // plot along a line
    template<size_t NDIM>
    void plot_along(World& world, trajectory<NDIM> traj, double (*ff)(const Vector<double,NDIM>&), std::string filename) {
    	 FILE *f=0;
    	 const int npt=traj.npt;

    	 const bool psdot=false;

    	 if(world.rank() == 0) {
    		 f = fopen(filename.c_str(), "w");
    		 if(!f) MADNESS_EXCEPTION("plotvtk: failed to open the plot file", 0);

    		 if (psdot) {
                 fprintf(f,"\\psset{xunit=0.05cm}\n");
                 fprintf(f,"\\psset{yunit=100cm}\n");
                 fprintf(f,"\\begin{pspicture}(0,0.25)(100,0.3)\n");
                 fprintf(f,"\\pslinewidth=0.005pt\n");
    		 }


    		 // walk along the line
    		 for (int ipt=0; ipt<npt; ipt++) {

    		     Vector<double,NDIM> coord=traj(ipt);
//    		     fprintf(f,"%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",coord[0],coord[1],coord[2],coord[3],coord[4],coord[5], ff(coord));
    			 if (psdot) {
                     // no hue code here
    			     //           long rank=ff.evalR(coord);
                     trajectory<NDIM>::print_psdot(f,ipt,ff(coord),trajectory<NDIM>::hueCode(0));
    			 } else {
    			     fprintf(f,"%4i %12.6f\n",ipt, ff(coord));
    			 }


    		 }
             if (psdot) fprintf(f,"\\end{pspicture}\n");

    		 fclose(f);
    	 }
    	 world.gop.fence();
    	 ;
    }



}

/* @} */
#endif // MADNESS_MRA_FUNCPLOT_H__INCLUDED
