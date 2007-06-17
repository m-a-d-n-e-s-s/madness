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
    double x = r[0]-0.5;
    double y = r[1]-0.5;
    double z = r[2]-0.5;
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

        Function<double,3> f = FunctionFactory<double,3>(world).f(myg).k(5).thresh(1e-7).nocompress();

        Vector<double,3> x = VectorFactory(0.5,0.5,0.5);
        print("the result is",f.eval(x).get()-myg(x));

        print("entering fence after eval");
        world.gop.fence();

        f.compress(false);
        print("entering fence after compress");
        world.gop.fence();
        



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
