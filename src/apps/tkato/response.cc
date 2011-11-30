#include <tkato/response.h>
using namespace madness;

static double ttt, sss;
static void START_TIMER(World& world) {
    world.gop.fence(); ttt=wall_time(); sss=cpu_time();
}

static void END_TIMER(World& world, const char* msg) {
    ttt=wall_time()-ttt; sss=cpu_time()-sss; if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
}

void CoupledPurturbation::guess_excite (int axis)
{
    // TODO: add selection to load stored functions
    START_TIMER(world);

    functionT dipole = factoryT(world).functor(functorT(new DipoleFunctor(axis)));
    rmo[axis] = vector<vecfuncT>(2*nXY);
    Vrmo[axis] = vector<vecfuncT>(2*nXY);
    for (int i=0; i<nXY; ++i) {
        rmo[axis][i*2] = mul(world, dipole, calc.amo);
        truncate(world, rmo[axis][i*2]);
        compress(world, rmo[axis][i*2]);
        if (!calc.param.spin_restricted && calc.param.nbeta) {
            rmo[axis][i*2+1] = mul(world, dipole, calc.bmo);
            truncate(world, rmo[axis][i*2+1]);
            compress(world, rmo[axis][i*2+1]);
        }
    }

    END_TIMER(world, "guess excite function");
}

void CoupledPurturbation::project_diag_space ()
{
    for (int i=0; i<nXY; ++i) {
        for (unsigned int j=0; j<rmo[axis][i*2].size(); ++j) {
            for (unsigned int k=0; k<calc.amo.size(); ++k) {
                double overlap = calc.amo[k].inner(rmo[axis][i*2][j]);
                rmo[axis][i*2][j].gaxpy(1.0, calc.amo[k], -overlap);
                // TODO: need to project out between alpha and beta ?
            }
        }
        if (!calc.param.spin_restricted && calc.param.nbeta) {
            for (unsigned int j=0; j<rmo[axis][i*2+1].size(); ++j) {
                for (unsigned int k=0; k<calc.bmo.size(); ++k) {
                    double overlap = calc.bmo[k].inner(rmo[axis][i*2+1][j]);
                    rmo[axis][i*2+1][j].gaxpy(1.0, calc.bmo[k], -overlap);
                    // TODO: need to project out between alpha and beta ?
                }
            }
        }
    }
}

double CoupledPurturbation::max_norm ()
{
    double max = 0.0;

    for (int i=0; i<nXY; ++i) {
        for (int j=0; j<rmo[axis][i*2].size(); ++j) {
            double norm = rmo[axis][i*2][j].norm2();
            if (max < norm) max = norm;
        }
        if (!calc.param.spin_restricted && calc.param.nbeta) {
            for (unsigned int j=0; j<rmo[axis][i*2+1].size(); ++j) {
                double norm = rmo[axis][i*2][j].norm2();
                if (max < norm) max = norm;
            }
        }
    }

    return max;
}

void CoupledPurturbation::make_vpsi (int axis)
{
    functionT rarho = factoryT(world);
    for (int i=0; i<nXY; ++i) {
        rarho += calc.make_density(world, calc.aocc, rmo[axis][i*2]);
    }
    rarho.truncate();
    functionT rbrho = rarho;
    if (!calc.param.spin_restricted && calc.param.nbeta) {
        rbrho = functionT(factoryT(world));
        for (int i=0; i<nXY; ++i) {
            rbrho += calc.make_density(world, calc.bocc, rmo[axis][i*2+1]);
        }
        rbrho.truncate();
    }
    functionT rrho = rarho + rbrho;
    rrho.truncate();
    functionT vcoul = apply(*calc.coulop, rrho);
    rrho.clear();
    functionT vlocal = vcoul + calc.vnuc;
    vlocal.truncate();
    functionT dummy;
    double exca = 0.0, excb = 0.0;
    for (int i=0; i<nXY; ++i) {
        Vrmo[axis][i*2] = calc.apply_potential(world, calc.aocc, rmo[axis][i*2], rarho, rbrho, dummy, dummy, vlocal, exca);
        if (calc.param.spin_restricted) {
            if (!calc.param.lda) excb = exca;
        }
        else if (calc.param.nbeta) {
            Vrmo = clac.apply_potential(world, calc.bocc, rmo[axis][i*2+1], rarho, rbrho, dummy, dummy, vlocal, excb);
        }
    }
}

vecfuncT CoupledPurturbation::calcBSH (vecfuncT & rmo, double eval)
{
    
}

// solve response linear equation
void CoupledPurturbation::solve (int axis, PurturbationOperator & purturb) {
    double scaling = 1.0;
    for (int i=0; i<calc.param.maxiter; i++) {
        if (world.rank() == 0) print("=== response iteration: ", i);

        // scaling factor
        double max = max_norm();
        if (max <= 0.9 || max >= 1.1) scaling *= max;

        project_diag_space();
        make_vpsi();
        double eval = freq / scaling;
    }
    START_TIMER(world);
    END_TIMER(world, "solve response");
}

void CoupledPurturbation::dipole_response ()
{
    vecfuncT dipoles(3);
    vector<PurturbationOperator> ops(3);
    for (int axis=0; axis<3; axis++) {
        dipoles[axis] = functionT(factoryT(world).functor(functorT(new DipoleFunctor(axis))));
        dipoles[axis].truncate();
        ops[axis] = PurturbationOperator(world, dipoles[axis], axis, freq);
        guess_excite(axis);
        solve(ops[axis]);
    }

    // print result
}

