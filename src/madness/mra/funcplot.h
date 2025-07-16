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
*/

#ifndef MADNESS_MRA_FUNCPLOT_H__INCLUDED
#define MADNESS_MRA_FUNCPLOT_H__INCLUDED

#include <madness/constants.h>
#include "QCCalculationParametersBase.h"
/*!

  \file mra/funcplot.h
  \brief Defines/implements plotting interface for functions
  \ingroup funcplot

  @{
 */

namespace madness {

    struct PlotParameters : public QCCalculationParametersBase {
        PlotParameters(World& world, const commandlineparser parser=commandlineparser(), const std::string tag="plot") : PlotParameters() {
            read_input_and_commandline_options(world,parser,tag);
        }

        PlotParameters() : QCCalculationParametersBase() {

            // initialize with: key, value, comment (optional), allowed values (optional)
            initialize<double>("zoom",2,"zoom into the simulation cell");
            initialize<long>("npoints",151,"number of plot points per dimension");
            initialize<std::vector<double>>("origin",{},"origin of the plot");
            initialize<std::vector<std::string>>("plane",{"x1","x2"},"plot plane: x1, x2, .., x6");
        }
        std::string get_tag() const override {
            return std::string("plot");
        }


        PlotParameters& set_zoom(const double z) {
            set_user_defined_value("zoom",z);
            return *this;
        }
        PlotParameters& set_npoints(const long n) {
            set_user_defined_value("npoints",n);
            return *this;
        }
        PlotParameters& set_plane(const std::vector<std::string> plane) {
            set_user_defined_value("plane",plane);
            return *this;
        }
        PlotParameters& set_origin(const std::vector<double> origin) {
            set_user_defined_value("origin",origin);
            return *this;
        }


        double zoom() const {return get<double>("zoom");}
        long npoints() const {return get<long>("npoints");}

        template<std::size_t NDIM>
        Vector<double,NDIM> origin() const {
            auto origin_vec=get<std::vector<double>>("origin");
            // fill in zeros if the default origin has fewer dimensions than the actual origin
            int missing=NDIM-origin_vec.size();
            for (auto i=0; i<missing; ++i) origin_vec.push_back(0.0);
            Vector<double,NDIM> o(origin_vec);
            return o;
        }

        std::vector<std::string> plane() const {return get<std::vector<std::string>>("plane");}


    };
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
    template <typename opT, std::size_t NDIM>
    void plot_line(World& world, const char* filename, int npt, const Vector<double,NDIM>& lo,
            const Vector<double,NDIM>& hi, const opT& op) {
        typedef Vector<double,NDIM> coordT;
        coordT h = (hi - lo)*(1.0/(npt-1));

        double sum = 0.0;
        for (std::size_t i=0; i<NDIM; ++i) sum += h[i]*h[i];
        sum = sqrt(sum);

        if (world.rank() == 0) {
            FILE* file = fopen(filename,"w");
        if(!file)
          MADNESS_EXCEPTION("plot_line: failed to open the plot file", 0);
            for (int i=0; i<npt; ++i) {
                coordT r = lo + h*double(i);
                fprintf(file, "%.14e ", i*sum);
                plot_line_print_value(file, op(r));
                fprintf(file,"\n");
            }
            fclose(file);
        }
        world.gop.fence();
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
    /// The ordinate is distance from lo
    template <typename T, std::size_t NDIM>
    void plot_line(const char* filename, int npt, const Vector<double,NDIM>& lo, const Vector<double,NDIM>& hi,
                   const std::vector<Function<T,NDIM>>& vf) {
        typedef Vector<double,NDIM> coordT;
        coordT h = (hi - lo)*(1.0/(npt-1));
        double sum = 0.0;
        for (std::size_t i=0; i<NDIM; ++i) sum += h[i]*h[i];
        sum = sqrt(sum);
        World& world = vf[0].world();// get world from first function
        // reconstruct each function in vf
        std::for_each(vf.begin(), vf.end(), [](const Function<T,NDIM>& f){f.reconstruct();});
        if (world.rank() == 0) {
            FILE* file = fopen(filename,"w");
            if(!file)
            MADNESS_EXCEPTION("plot_line: failed to open the plot file", 0);
            for (int i=0; i<npt; ++i) {
                coordT r = lo + h*double(i);
                fprintf(file, "%.14e ", i*sum);
                std::for_each(vf.begin(), vf.end(), [&](const Function<T,NDIM>& f){ plot_line_print_value(file, f.eval(r));});
                fprintf(file,"\n");
            }
            fclose(file);
        }
        world.gop.fence();
    }

    template<size_t NDIM>
    void plot_plane(World& world, const Function<double,NDIM>& function,
            const std::string name) {
        typedef std::vector<Function<double,NDIM> > vecfuncT;
        plot_plane(world,vecfuncT(1,function),name);
    }

    template<size_t NDIM>
    void plot_plane(World& world, const Function<double,NDIM>& function1,
            const Function<double,NDIM>& function2,
            const std::string name) {
        typedef std::vector<Function<double,NDIM> > vecfuncT;
        vecfuncT vf(2);
        vf[0]=function1;
        vf[1]=function2;
        plot_plane(world,vf,name);
    }

    template<size_t NDIM>
    void plot_plane(World& world, const Function<double,NDIM>& function1,
            const Function<double,NDIM>& function2,const Function<double,NDIM>& function3,
            const std::string name) {
        typedef std::vector<Function<double,NDIM> > vecfuncT;
        vecfuncT vf(3);
        vf[0]=function1;
        vf[1]=function2;
        vf[2]=function3;
        plot_plane(world,vf,name);
    }


    /// plot a 2-d slice of a given function and the according MRA structure
    /// FIXME: doesn't work for more than 1 rank

    /// the plotting parameters are taken from the input file "input" and its
    /// data group "plot", e.g. plotting the xy plane around (0,0,0.7):
    /// plot
    ///   plane x1 x2
    ///   zoom 2.0
    ///   npoints 100
    ///   origin 0.0 0.0 0.7
    /// end
    /// @param[in]	world	the world
    /// @param[in]	vfunction	the function to plot
    /// @param[in]	name		the output name
    template<size_t NDIM>
    void plot_plane(World& world, const std::vector<Function<double,NDIM> >& vfunction,
    		const std::string name, const std::string inputfile="input") {

		if (world.size()>1) return;
        // determine the ploting plane
    	std::string c1="x1", c2="x2";

    	// zoom factor
    	double zoom=1.0;

    	// output type: mathematica or gnuplot
    	std::string output_type="gnuplot";

    	// number of points in each direction
        int npoints=200;

        // the coordinates to be plotted
        Vector<double,NDIM> coord(0.0);
        Vector<double,NDIM> origin(0.0);

        try {
            std::ifstream f(inputfile);
            position_stream_to_word(f, "plot",'#',true,true);
            std::string s;
            while (f >> s) {
                if (s == "end") {
                    break;
                } else if (s == "plane") {
                    f >> c1 >> c2;
                } else if (s == "zoom") {
                    f >> zoom;
                } else if (s == "output") {
                    f >> output_type;
                } else if (s == "points") {
                    f >> npoints;
                } else if (s == "origin") {
                    for (std::size_t i=0; i<NDIM; ++i) f >> origin[i];
                }
            }
        } catch (...) {
            print("can't locate plot in file="+inputfile+" -- using default values");
        }
    	double scale=1.0/zoom;
    	coord=origin;

        // convert human to mad form
        size_t cc1=0, cc2=1;
        if (c1=="x1") cc1=0;
        if (c1=="x2") cc1=1;
        if (c1=="x3") cc1=2;
        if (c1=="x4") cc1=3;
        if (c1=="x5") cc1=4;
        if (c1=="x6") cc1=5;
        if (c2=="x1") cc2=0;
        if (c2=="x2") cc2=1;
        if (c2=="x3") cc2=2;
        if (c2=="x4") cc2=3;
        if (c2=="x5") cc2=4;
        if (c2=="x6") cc2=5;

        MADNESS_ASSERT(cc1<NDIM);
        MADNESS_ASSERT(cc2<NDIM);
        // output file name for the gnuplot data
        std::string filename="plane_"+c1+c2+"_"+name;
        // assume a cubic cell
        double lo=-FunctionDefaults<NDIM>::get_cell_width()[0]*0.5;
        lo=lo*scale;

        const double stepsize=FunctionDefaults<NDIM>::get_cell_width()[0]*scale/npoints;

        if(world.rank() == 0) {

        	// plot 3d plot
        	FILE *f =  0;
        	f=fopen(filename.c_str(), "w");
        	if(!f) MADNESS_EXCEPTION("plot_along: failed to open the plot file", 0);

        	for (int i0=0; i0<npoints; i0++) {
        		for (int i1=0; i1<npoints; i1++) {
        			// plot plane
        			coord[cc1]=lo+origin[cc1]+i0*stepsize;
        			coord[cc2]=lo+origin[cc2]+i1*stepsize;

        			// other electron
//        			fprintf(f,"%12.6f %12.6f %12.20f\n",coord[cc1],coord[cc2],
//        					function(coord));
                    fprintf(f,"%12.6f %12.6f",coord[cc1],coord[cc2]);
                    for (std::size_t ivec=0; ivec<vfunction.size(); ++ivec)
                        fprintf(f,"  %12.20f",vfunction[ivec](coord));
                    fprintf(f,"\n");

        		}
        		// additional blank line between blocks for gnuplot
        		if (output_type=="gnuplot") fprintf(f,"\n");
        	}
        	fclose(f);

        }

//        // plot mra structure
//    	filename="mra_structure_"+c1+c2+"_"+name;
//    	function.get_impl()->print_plane(filename.c_str(),cc1,cc2,coord);
    }


/// plot a 2-d slice of a given function and the according MRA structure

/// the plotting parameters are taken from the input file "input" and its
/// data group "plot", e.g. plotting the xy plane around (0,0,0.7):
/// plot
///   plane x1 x2
///   zoom 2.0
///   npoints 100
///   origin 0.0 0.0 0.7
/// end
/// @param[in]	world	the world
/// @param[in]	vfunction	the function to plot
/// @param[in]	name		the output name
template<size_t NDIM>
void plot_plane(World& world, const std::vector<Function<double,NDIM> >& vfunction,
                const std::string name, const PlotParameters param) {

    if (world.size()>1) return;

    auto plane=param.plane();
    std::string c1=plane[0];
    std::string c2=plane[1];
    auto npoints=param.npoints();
    auto origin=param.origin<NDIM>();
    auto coord=param.origin<NDIM>();
    double scale=1.0/param.zoom();
    std::string output_type="gnuplot";

    auto plane2dim = [](std::string c) {
        if (c=="x1") return 0;
        else if (c=="x2") return 1;
        else if (c=="x3") return 2;
        else if (c=="x4") return 3;
        else if (c=="x5") return 4;
        else if (c=="x6") return 5;
        else return -1;
    };

    // convert human to mad form
    std::size_t cc1=plane2dim(c1);
    std::size_t cc2=plane2dim(c2);

    MADNESS_ASSERT(cc1<NDIM);
    MADNESS_ASSERT(cc2<NDIM);

    // output file name for the gnuplot data
    std::string filename="plane_"+c1+c2+"_"+name;
    // assume a cubic cell
    double lo=-FunctionDefaults<NDIM>::get_cell_width()[0]*0.5;
    lo=lo*scale;

    const double stepsize=FunctionDefaults<NDIM>::get_cell_width()[0]*scale/npoints;

    if(world.rank() == 0) {

        // plot 3d plot
        FILE *f =  0;
        f=fopen(filename.c_str(), "w");
        if(!f) MADNESS_EXCEPTION("plot_along: failed to open the plot file", 0);

        for (int i0=0; i0<npoints; i0++) {
            for (int i1=0; i1<npoints; i1++) {
                // plot plane
                coord[cc1]=lo+origin[cc1]+i0*stepsize;
                coord[cc2]=lo+origin[cc2]+i1*stepsize;

                fprintf(f,"%12.6f %12.6f",coord[cc1],coord[cc2]);
                for (std::size_t ivec=0; ivec<vfunction.size(); ++ivec)
                    fprintf(f,"  %12.20f",vfunction[ivec](coord));
                fprintf(f,"\n");

            }
            // additional blank line between blocks for gnuplot
            if (output_type=="gnuplot") fprintf(f,"\n");
        }
        fclose(f);

    }
}


    template<size_t NDIM, typename opT>
    void plot_plane(World& world, const opT& op, const std::string name, const PlotParameters param) {

         if (world.size()>1) return;

        auto plane=param.plane();
        std::string c1=plane[0];
        std::string c2=plane[1];
        auto npoints=param.npoints();
        auto origin=param. template origin<NDIM>();
        auto coord=param. template origin<NDIM>();
        double scale=1.0/param.zoom();
        std::string output_type="gnuplot";

        auto plane2dim = [](std::string c) {
            if (c=="x1") return 0;
            else if (c=="x2") return 1;
            else if (c=="x3") return 2;
            else if (c=="x4") return 3;
            else if (c=="x5") return 4;
            else if (c=="x6") return 5;
            else return -1;
        };

        // convert human to mad form
        std::size_t cc1=plane2dim(c1);
        std::size_t cc2=plane2dim(c2);

        MADNESS_ASSERT(cc1<NDIM);
        MADNESS_ASSERT(cc2<NDIM);

         // output file name for the gnuplot data
         std::string filename="plane_"+c1+c2+"_"+name;
         // assume a cubic cell
         double lo=-FunctionDefaults<NDIM>::get_cell_width()[0]*0.5;
         lo=lo*scale;

         const double stepsize=FunctionDefaults<NDIM>::get_cell_width()[0]*scale/npoints;

         if(world.rank() == 0) {

             // plot 3d plot
             FILE *f =  0;
             f=fopen(filename.c_str(), "w");
             if(!f) MADNESS_EXCEPTION("plot_along: failed to open the plot file", 0);

             for (int i0=0; i0<npoints; i0++) {
                 for (int i1=0; i1<npoints; i1++) {
                     // plot plane
                     coord[cc1]=lo+origin[cc1]+i0*stepsize;
                     coord[cc2]=lo+origin[cc2]+i1*stepsize;

                     // other electron
 //                  fprintf(f,"%12.6f %12.6f %12.20f\n",coord[cc1],coord[cc2],
 //                          function(coord));
                     fprintf(f,"%12.6f %12.6f",coord[cc1],coord[cc2]);
                     fprintf(f,"  %12.20f\n",op(coord));

                 }
                 // additional blank line between blocks for gnuplot
                 if (output_type=="gnuplot") fprintf(f,"\n");
             }
             fclose(f);

         }
     }



    template<size_t NDIM>
    typename std::enable_if<NDIM==3,void>::type
    plot_cubefile(World& world, const Function<double,NDIM>& f, std::string filename,
            std::vector<std::string> molecular_info=std::vector<std::string>(), int npoints=100, double zoom=1.0,
            const Vector<double,NDIM> origin=Vector<double,NDIM>(0.0)) {

        if (world.size()>1) return;

        // dummy atom in the center
        if (molecular_info.size()==0)
        	molecular_info=std::vector<std::string>(1,"0 0 0.0 0.0 0.0\n");

        // the coordinates to be plotted
        // Vector<double,NDIM> origin(0.0);

        // number of points in each direction
        std::vector<int> npt(3,npoints);

        Tensor<double> cell=copy(FunctionDefaults<3>::get_cell());
        cell.scale(1.0/zoom);
        double xlen=cell(0,1)-cell(0,0);
        double ylen=cell(1,1)-cell(1,0);
        double zlen=cell(2,1)-cell(2,0);

        // plot file
        FILE *file =  0;
        file=fopen(filename.c_str(), "w");
        if(!file) MADNESS_EXCEPTION("plot_along: failed to open the plot file", 0);


        // print header
        fprintf(file,"cube file from MADNESS\n");
        fprintf(file,"comment line\n");

        // print the number of atoms if a calculation was provided
        fprintf(file,"%d %12.8f %12.8f %12.8f \n",int(molecular_info.size()),
                cell(0,0),cell(1,0),cell(2,0));

        // grid spacing for each dimension such that the cell edges are plotted
        const auto xdelta = xlen/(npt[0]-1);
        const auto ydelta = ylen/(npt[1]-1);
        const auto zdelta = zlen/(npt[2]-1);

        // print the cell constants
        fprintf(file,"%d %12.6f %12.6f %12.6f\n",npt[0],xdelta,0.0,0.0);
        fprintf(file,"%d %12.6f %12.6f %12.6f\n",npt[1],0.0,ydelta,0.0);
        fprintf(file,"%d %12.6f %12.6f %12.6f\n",npt[2],0.0,0.0,zdelta);

        // print the molecule
        for (const std::string& s : molecular_info) fprintf(file,"%s",s.c_str());


         Tensor<double> grid(npt[0], npt[1], npt[2]);
         long count_per_line = 0;
         for (int i = 0; i < npt[0]; ++i) {
             for (int j = 0; j < npt[1]; ++j) {
                 for (int k = 0; k < npt[2]; ++k) {
                     double x = cell(0, 0) + origin[0] + xdelta * i;
                     double y = cell(1, 0) + origin[1] + ydelta * j;
                     double z = cell(2, 0) + origin[2] + zdelta * k;
                     // the original format has up to 6 entries per line: https://paulbourke.net/dataformats/cube/
                     // stick with this, even though many codes can read an arbitrary number of entries per line
                     if (count_per_line < 6) {
                        fprintf(file, "%12.5e ", f(x, y, z));
                        count_per_line++;
                     } else {
                        fprintf(file, "%12.5e\n", f(x, y, z));
                        count_per_line = 0;
                     }
                 }
             }
         }
         fprintf(file, "\n");
         fclose(file);
     }

     template<typename T, size_t NDIM>
     void
     print_tree_jsonfile(World& world, const Function<T,NDIM>& f, std::string filename) {

         if (world.size() > 1)
             return;

         Tensor<double> cell = copy(FunctionDefaults<NDIM>::get_cell());
         std::ofstream os(filename.c_str());

         os << "{";
         os << "\"cell\":[";
         for (int xyz = 0; xyz != NDIM; ++xyz) {
             os << "[" << cell(xyz, 0) << "," << cell(xyz, 1) << "]";
             if (xyz != NDIM-1)
                 os << ",";
         }
         os << "],";

         os << "\"tree\":{";
         f.print_tree_json(os);
         os << "}}";
     }

    /// convenience to get plot_plane and plot_cubefile
    template<size_t NDIM>
    void plot(const std::vector<Function<double,NDIM> >& vf, const std::string& name, const std::vector<std::string>& header){
    	if(vf.empty()) return;
    	World& world=vf.front().world();
    	for(size_t i=0;i<vf.size();++i){
    		const std::string namei=name+"_"+std::to_string(i);
//    		vf[i].print_size("plot:"+namei);
    		plot_plane<NDIM>(world,vf[i],namei);
//    		plot_cubefile<NDIM>(world,vf[i],namei+".cube",header);
    	}
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
