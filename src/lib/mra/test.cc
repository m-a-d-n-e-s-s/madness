/*
  This file is part of MADNESS.
  
  Copyright (C) <2007> <Oak Ridge National Laboratory>
  
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

  
#include <mra/mra.h>
//#include <mra/loadbal.h>

using namespace madness;

double myfun(const double x[]) {
    double r2=0.0;
    for (int i=0; i < 3; i++)
	r2 += x[i]*x[i];
    return r2;
}

const double PI = 3.1415926535897932384;
const double myg_expnt = 120.0;

double myg(const Vector<double,3>& r) {
    /* A square-normalized gaussian */
    double fac = pow(2.0*myg_expnt/PI,0.75);
    double x = r[0]-0.7;
    double y = r[1]-0.5;
    double z = r[2]-0.3;
    return fac*exp(-myg_expnt*(x*x + y*y + z*z));
};

void vector_myg(long npt, const double *x, const double *y, 
                const double *z, double* RESTRICT f) {
    const double fac = pow(2.0*myg_expnt/PI,0.75);
    for (int i=0; i<npt; i++) {
        double xx = x[i] - 0.5;
        double yy = y[i] - 0.5;
        double zz = z[i] - 0.5;
        f[i] = fac*exp(-myg_expnt*(xx*xx + yy*yy + zz*zz));
    }
};

// double dmygdx(const double r[3]) {
//     /* Derivative of myg w.r.t. x */
//     return -2.0*myg_expnt*(r[0]-0.5)*myg(r);
// };

// double dmygdy(const double r[3]) {
//     /* Derivative of myg w.r.t. y */
//     return -2.0*myg_expnt*(r[1]-0.5)*myg(r);
// };

// double dmygdz(const double r[3]) {
//     /* Derivative of myg w.r.t. z */
//     return -2.0*myg_expnt*(r[2]-0.5)*myg(r);
// };
    

int main(int argc, char**argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);

    try {
        startup(world,argc,argv);

	double t0 = MPI::Wtime();
//        Function<double,3> f = FunctionFactory<double,3>(world).f(myg).k(5).thresh(1e-7).nocompress();
        Function<double,3> f = FunctionFactory<double,3>(world).f(myg).k(11).thresh(1e-12).nocompress();
	double t1 = MPI::Wtime();
	Function<double,3> g = copy(f);

	double t2 = MPI::Wtime();
	LoadBalImpl<double,3> lb(g);
	double t3 = MPI::Wtime();
	FunctionDefaults<3>::pmap = lb.loadBalance();
	double t4 = MPI::Wtime();
	f = copy(g);
	double t5 = MPI::Wtime();

	madness::print("Routine            |  Time");
	madness::print("-------------------+--------------");
	madness::print("Init f             | ", t1-t0);
	madness::print("Copy f to g        | ", t2-t1);
	madness::print("Create LoadBalImpl | ", t3-t2);
	madness::print("Load Balance       | ", t4-t3);
	madness::print("Copy g to f        | ", t5-t4);
	madness::print("-------------------+--------------");
	madness::print("Total Time         | ", t5-t0);
	

        //print("The tree after projection");
        //f.print_tree();
            

//         Vector<double,3> x = VectorFactory(0.5,0.5,0.5);
//         print("the result is",f.eval(x).get()-myg(x));

//         print("entering fence after eval");
//         world.gop.fence();

//         f.compress(false);
//         print("entering fence after compress");
//         world.gop.fence();
//         print("The tree after compress");
//         f.print_tree();
        
//         f.reconstruct(false);
//         print("entering fence after reconstruct");
//         world.gop.fence();
//         print("The tree after reconstruct");
//         f.print_tree();
        
	


//	xterm_debug("test", 0);

	//double t0, t1, t2, t3;
        
        //Function<double,3,MyProcmap<3> > f = FunctionFactory<double,3,MyProcmap<3> >(world).f(myfun).k(3).thresh(1e-2).nocompress();

// 	print("about to construct LoadBalImpl");
// 	t0 = wall_time();
// 	LoadBalImpl<double,3,MyProcmap<3> > lbi(f);
// 	t1 = wall_time();
// 	print("constructed LoadBalImpl, time =", t1-t0);
// 	t2 = wall_time();
// 	lbi.loadBalance();
// 	t3 = wall_time();
// 	print("load balanced, time =", t3-t2);
    } catch (const MPI::Exception& e) {
        print(e);
        error("caught an MPI exception");
    } catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    } catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    } catch (const char* s) {
        print(s);
        error("caught a string exception");
    } catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    } catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    } catch (...) {
        error("caught unhandled exception");
    }

    print("entering final fence");
    world.gop.fence();
    print("done with final fence");
    MPI::Finalize();

    return 0;
}
