#include <mra/mra.h>
#include <world/loadbal.h>

using namespace madness;

double myfun(const double x[]) {
    double r2=0.0;
    for (int i=0; i < 3; i++)
	r2 += x[i]*x[i];
    return r2;
}

int main(int argc, char**argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);

    try {
        startup(world,argc,argv);

/*
        Key<4> key(0,Array<Translation,4>(0));
        print("Initial key",key);
        for (KeyChildIterator<4> it(key); it; ++it) 
            print(it.key());
*/
	
        Function<double,3,MyProcmap<3> > f = FunctionFactory<double,3,MyProcmap<3> >(world).f(myfun).k(3).thresh(1e-2).nocompress();

    } catch (MPI::Exception e) {
        error("caught an MPI exception");
    } catch (madness::MadnessException e) {
        print(e);
        error("caught a MADNESS exception");
    } catch (const char* s) {
        print(s);
        error("caught a string exception");
    } catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    } catch (...) {
        error("caught unhandled exception");
    }

    print("entering final fence");
    world.gop.fence();
    print("done with final fence");
    MPI::Finalize();

    return 0;
}
