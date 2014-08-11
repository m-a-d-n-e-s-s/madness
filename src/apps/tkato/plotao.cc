//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/mra/mra.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <cmath>
using namespace madness;

typedef std::shared_ptr< WorldDCPmapInterface< Key<3> > > pmapT;
typedef Vector<double,3> coordT;
typedef std::shared_ptr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef std::vector<functionT> vecfuncT;
typedef Tensor<double> tensorT;
typedef FunctionFactory<double,3> factoryT;

class LevelPmap : public WorldDCPmapInterface< Key<3> > {
    private:
    const int nproc;
    public:
    LevelPmap() : nproc(0) {};

    LevelPmap(World& world) : nproc(world.nproc()) {}

    /// Find the owner of a given key
    ProcessID owner(const Key<3>& key) const {
        Level n = key.level();
        if (n == 0) return 0;
        hashT hash;
        if (n <= 3 || (n&0x1)) hash = key.hash();
        else hash = key.parent().hash();
        //hashT hash = key.hash();
        return hash%nproc;
    }
};

/// Given overlap matrix, return rotation with 3rd order error to orthonormalize the vectors
tensorT Q3(const tensorT& s) {
    tensorT Q = inner(s,s);
    Q.gaxpy(0.2,s,-2.0/3.0);
    for (int i=0; i<s.dim(0); ++i) Q(i,i) += 1.0;
    return Q.scale(15.0/8.0);
}

class AO {
public:
    const World& world;
    const char* filename;
    double vtol;

    bool spin_restricted;
    vecfuncT amo, bmo;
    std::vector<int> aset, bset;
    tensorT aocc, bocc;
    tensorT aeps, beps;
    unsigned int nmo;

    AO (World & world, const char* filename) : world(world), filename(filename)
    {
        double safety = 0.1;
        vtol = FunctionDefaults<3>::get_thresh() * safety;
        if (world.rank() == 0) {
            std::cout << "plotting orbitals contained in " << filename << std::endl;
        }
    }

    void load_mos (World & world)
    {
        const double thresh = FunctionDefaults<3>::get_thresh();
        const int k = FunctionDefaults<3>::get_k();
        unsigned int nmo;
        amo.clear(); bmo.clear();

        archive::ParallelInputArchive ar(world, filename);

        /*
           File format:

           bool spinrestricted --> if true only alpha orbitals are present

           unsigned int nmo_alpha;
           Tensor<double> aeps;
           Tensor<double> aocc;
           vector<int> aset;
           for i from 0 to nalpha-1:
           .   Function<double,3> amo[i]
           repeat for beta if !spinrestricted

        */

        ar & spin_restricted;

        ar & nmo;
        this->nmo = nmo;
        const double trantol = vtol / std::min(30.0, double(nmo));
        ar & aeps & aocc & aset;
        if (world.rank() == 0) {
            std::cout << " anmo :" << nmo << std::endl;
            std::cout << " aeps :" << aeps << std::endl;
            std::cout << " aocc :" << aocc << std::endl;
            std::cout << " aset :" << aset << std::endl;
        }
        amo.resize(nmo);
        for (unsigned int i=0; i<amo.size(); ++i) ar & amo[i];

        if (amo[0].k() != k) {
            reconstruct(world,amo);
            for(unsigned int i = 0;i < amo.size();++i) amo[i] = madness::project(amo[i], k, thresh, false);
            world.gop.fence();
        }
        normalize(world, amo);
        amo = transform(world, amo, Q3(matrix_inner(world, amo, amo)), trantol, true);
        truncate(world, amo);
        normalize(world, amo);

        if (!spin_restricted) {
            ar & nmo;
            if (nmo > 0) {
                this->nmo += nmo;
                ar & beps & bocc & bset;
                if (world.rank() == 0) {
                    std::cout << " bnmo :" << nmo << std::endl;
                    std::cout << " beps :" << beps << std::endl;
                    std::cout << " bocc :" << bocc << std::endl;
                    std::cout << " bset :" << bset << std::endl;
                }
                bmo.resize(nmo);
                for (unsigned int i=0; i<bmo.size(); ++i) ar & bmo[i];

                if (bmo[0].k() != k) {
                    reconstruct(world,bmo);
                    for(unsigned int i = 0;i < bmo.size();++i) bmo[i] = madness::project(bmo[i], k, thresh, false);
                    world.gop.fence();
                }

                normalize(world, bmo);
                bmo = transform(world, bmo, Q3(matrix_inner(world, bmo, bmo)), trantol, true);
                truncate(world, bmo);
                normalize(world, bmo);
            }
        }
        if (world.rank() == 0) {
            std::cout << "total nmo :" << this->nmo << std::endl;
        }
    }

    void plotline_z (const unsigned int npt, const double lo, const double hi, unsigned int index) const
    {
        coordT lov(0.0);
        lov[2] = lo;
        coordT hiv(0.0);
        hiv[2] = hi;

        static const char fnbase[] = "ao_%sz_%s_%02d.txt";
        char fn[64];
        if (amo.size() <= index && !spin_restricted) {
            index -= amo.size();
            sprintf(fn, fnbase, "", "beta", index);
            if (world.rank() == 0) { std::cout << "making " << fn << std::endl; }
            plot_line(fn, npt, lov, hiv, bmo[index]);
            sprintf(fn, fnbase, "-", "beta", index);
            if (world.rank() == 0) { std::cout << "making " << fn << std::endl; }
            plot_line(fn, npt, lov, hiv, (-1.0) * bmo[index]);
        }
        else {
            sprintf(fn, fnbase, "", "alpha", index);
            if (world.rank() == 0) { std::cout << "making " << fn << std::endl; }
            plot_line(fn, npt, lov, hiv, amo[index]);
            sprintf(fn, fnbase, "-", "alpha", index);
            if (world.rank() == 0) { std::cout << "making " << fn << std::endl; }
            plot_line(fn, npt, lov, hiv, (-1.0) * amo[index]);
        }
    }

    void plotradius (const unsigned int npt, const double lo, const double hi, unsigned int index) const
    {
        static const char fnbase[] = "ao_r_%s_%02d.txt";
        char fn[64];
        if (amo.size() <= index && !spin_restricted) {
            index -= amo.size();
            sprintf(fn, fnbase, "beta", index);
            if (world.rank() == 0) { std::cout << "making " << fn << "(norm:" << bmo[index].norm2() << ")" << std::endl; }
            plotradius_f(fn, npt, lo, hi, bmo[index]);
        }
        else {
            sprintf(fn, fnbase, "alpha", index);
            if (world.rank() == 0) { std::cout << "making " << fn << "(norm:" << amo[index].norm2() << ")" << std::endl; }
            plotradius_f(fn, npt, lo, hi, amo[index]);
        }
    }

private:
    void plotradius_printvalue (std::ofstream& f, double r, double v) const
    {
        f << r << " " << v << std::endl;
    }

    void plotradius_f (const char filename[], const unsigned int npt, double lo, double hi, const functionT & mo) const
    {
        if (lo > hi) std::swap(lo, hi);
        std::ofstream f(filename);
        f.precision(14);
        f << std::scientific;

        const double dr = (hi - lo) / (npt - 1);
        const double mid = 0.7 * (hi - lo);
        const unsigned int nint = static_cast<unsigned int>(sqrt(2.0 * constants::pi * constants::pi * 0.8 * mid * mid / (dr * dr)));
        const double dtheta = constants::pi / (nint - 1);
        const double dphi = 2.0 * constants::pi / nint;
        if (world.rank() == 0) {
            std::cout << "n of integral = " << nint << std::endl;
            std::cout << "(dr, dtheta, dphi) = (" << dr << ", " << dtheta << ", " << dphi << ")" << std::endl;
        }
        double sumoverall = 0.0;
        double sumsqoverall = 0.0;
        for (unsigned int k=0; k<npt; ++k) {
            double r = dr * k + lo;
            if (world.rank() == 0 && k%10 == 9) std::cout << ".";
            double sum = 0.0;
            double sumsq = 0.0;
            for (unsigned int i=0; i<nint; ++i) {
                double theta = dtheta * i;
                double z = r * cos(theta);
                double rsint = r * sin(theta);
                double rsqdrtp = r * rsint * dtheta * dphi;
                for (unsigned int j=0; j < nint; ++j) {
                    double phi = dphi * j;
                    double x = rsint * cos(phi);
                    double y = rsint * sin(phi);
                    coordT point;
                    point[0] = x;
                    point[1] = y;
                    point[2] = z;
                    double fvalue = mo.eval(point);
                    sum += fvalue * rsqdrtp;
                    sumsq += fvalue * fvalue * rsqdrtp;
                }
            }
            if (world.rank() == 0) plotradius_printvalue(f, r, sumsq);
            sumoverall += sum * dr;
            sumsqoverall += sumsq * dr;
        }
        if (world.rank() == 0) {
            std::cout << std::endl;
            std::cout << "\\iiint f(r,theta,phi) r^2 sin(theta) dr dtheta dphi = " << sumoverall << std::endl;
            std::cout << "\\iiint f(r,theta,phi)^2 r^2 sin(theta) dr dtheta dphi = " << sumsqoverall << std::endl;
        }
    }
};

void set_protocol(World & world, double thresh)
{
    int k;
    if(thresh >= 1e-2)
        k = 4;
    else if(thresh >= 1e-4)
        k = 6;
    else if(thresh >= 1e-6)
        k = 8;
    else if(thresh >= 1e-8)
        k = 10;
    else
        k = 12;

    FunctionDefaults<3>::set_cubic_cell(-50, 50);
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_initial_level(2);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_autorefine(false);
    FunctionDefaults<3>::set_apply_randomize(false);
    FunctionDefaults<3>::set_project_randomize(false);
    GaussianConvolution1DCache<double>::map.clear();
    if(world.rank() == 0){
        print("\nSolving with thresh", thresh, "    k", FunctionDefaults<3>::get_k(), "\n");
    }
}

int main (int argc, char* argv[]) {
    initialize(argc, argv);

    {
        World world(SafeMPI::COMM_WORLD);

        try {
            startup(world, argc, argv);
            FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap(world)));

            set_protocol(world, 1e-6);

            char filename[256];
            char* envfn = getenv("MAD_PLOTINPUT");
            if (envfn == NULL || strlen(envfn) < 1) {
                strcpy(filename, "restartdata");
            }
            else {
                strcpy(filename, envfn);
            }
            AO ao(world, filename);
            ao.load_mos(world);
            for (unsigned int i=0; i<ao.nmo; ++i) {
                if (world.rank() == 0) std::cout << std::endl;
                ao.plotline_z(1024, 0.0, 4.0, i);
                ao.plotradius(256, 0.0, 4.0, i);
            }
        }
        catch (const SafeMPI::Exception& e) {
            print(e);
            error("caught an MPI exception");
        }
        catch (const madness::MadnessException& e) {
            print(e);
            error("caught a MADNESS exception");
        }
        catch (const madness::TensorException& e) {
            print(e);
            error("caught a Tensor exception");
        }
        catch (char* s) {
            print(s);
            error("caught a string exception");
        }
        catch (const char* s) {
            print(s);
            error("caught a string exception");
        }
        catch (const std::string& s) {
            print(s);
            error("caught a string (class) exception");
        }
        catch (const std::exception& e) {
            print(e.what());
            error("caught an STL exception");
        }
        catch (...) {
            error("caught unhandled exception");
        }

        // Nearly all memory will be freed at this point
        world.gop.fence();
        world.gop.fence();
        print_stats(world);
    }

    finalize();

    return 0;
}

