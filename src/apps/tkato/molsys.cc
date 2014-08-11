#include <madness/mra/mra.h>
#include <tkato/scf.h>
#include <tkato/molsys.h>

void MolecularSystem::save (World & world, const char * filename)
{
    archive::ParallelOutputArchive ar(world, filename, nio);
    ar & spin_restricted;
    ar & (unsigned int)(amo.size());
    ar & aeps & aocc & aset;
    for (unsigned int i=0; i<amo.size(); ++i) ar & amo[i];
    if (!spin_restricted) {
        ar & (unsigned int)(bmo.size());
        ar & beps & bocc & bset;
        for (unsigned int i=0; i<bmo.size(); ++i) ar & bmo[i];
    }
}

void MolecularSystem::load (World & world, const char * filename, unsigned int nalpha, unsigned int nvalpha, unsigned int nbeta, unsigned int nvbeta, unsigned int n_core)
{
    double safety = 0.1;
    double vtol = FunctionDefaults<3>::get_thresh() * safety;
    const double trantol = vtol / std::min(30.0, static_cast<double>(nalpha));
    const double thresh = FunctionDefaults<3>::get_thresh();
    const int k = FunctionDefaults<3>::get_k();
    unsigned int nmo;
    bool spinrest;
    unsigned int nmo_alpha = nalpha + nvalpha;
    unsigned int nmo_beta = nbeta + nvbeta;
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

    // LOTS OF LOGIC MISSING HERE TO CHANGE OCCUPATION NO., SET,
    // EPS, SWAP, ... sigh

    ar & spinrest;

    ar & nmo;
    MADNESS_ASSERT(nmo >= nmo_alpha);
    ar & aeps & aocc & aset;
    amo.resize(nmo);
    for (unsigned int i=0; i<amo.size(); ++i) ar & amo[i];
    if (nmo > nmo_alpha) {
        aset = vector<int>(aset.begin()+n_core, aset.begin()+n_core+nmo_alpha);
        amo = vecfuncT(amo.begin()+n_core, amo.begin()+n_core+nmo_alpha);
        aeps = copy(aeps(Slice(n_core, n_core+nmo_alpha-1)));
        aocc = copy(aocc(Slice(n_core, n_core+nmo_alpha-1)));
    }

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

        if (spinrest) { // Only alpha spin orbitals were on disk
            MADNESS_ASSERT(nmo_alpha >= nmo_beta);
            bmo.resize(nmo_beta);
            bset.resize(nmo_beta);
            beps = copy(aeps(Slice(0,nmo_beta-1)));
            bocc = copy(aocc(Slice(0,nmo_beta-1)));
            for (unsigned int i=0; i<nmo_beta; ++i) bmo[i] = copy(amo[i]);
        }
        else {
            ar & nmo;
            ar & beps & bocc & bset;

            bmo.resize(nmo);
            for (unsigned int i=0; i<bmo.size(); ++i) ar & bmo[i];

            if (nmo > unsigned(nmo_beta)) {
                bset = vector<int>(bset.begin()+n_core, bset.begin()+n_core+nmo_beta);
                bmo = vecfuncT(bmo.begin()+n_core, bmo.begin()+n_core+nmo_beta);
                beps = copy(beps(Slice(n_core, n_core+nmo_beta-1)));
                bocc = copy(bocc(Slice(n_core, n_core+nmo_beta-1)));
            }

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
}

