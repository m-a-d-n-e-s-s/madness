#include "densities.h"
#include "loadbalance.h"
#include "states.h"

// Make vector of gaussians around nucleon coordinates
void make_MD(World& world, comp_vecfuncT& psi_n,
                           comp_vecfuncT& psi_p,
                           const double A,
                           const double Z,
                           const double L,
                           const double prec)
{
    double res = FunctionDefaults<3>::get_thresh()*prec;

    FILE *partn = fopen("MD/neutrons.tab", "r");
    if (world.rank() == 0) {std::cout << "- Begin: Read neutron data" << std::endl;}
    for (size_t i = 0; i < (A-Z) ; i++) {
        double x, y, z;
        fscanf(partn, "%lf \t%lf \t%lf\n", &x, &y, &z);
        x -= L;
        y -= L;
        z -= L;
        psi_n[i] = comp_factoryT(world).functor(comp_functorT(new MD(x, y, z, L)));
        psi_n[i].truncate(res);
    }
    fclose(partn);
    world.gop.fence();

    FILE *partp = fopen("MD/protons.tab", "r");
    if(world.rank() == 0) {std::cout << "- Begin: Read proton data" << std::endl;}
    for (size_t i = 0; i < Z ; i++) {
        double x, y, z;
        fscanf(partp, "%lf \t%lf \t%lf\n", &x, &y, &z);
        x -= L;
        y -= L;
        z -= L;
        psi_p[i] = comp_factoryT(world).functor(comp_functorT(new MD(x, y, z, L)));
        psi_p[i].truncate(res);
    }
    fclose(partp);
    world.gop.fence();
    if(world.rank() == 0) {std::cout << "- End: Read data" << std::endl;}
}



// Make vector of modified HO wavefunctions
void make_HO(World& world, comp_vecfuncT& u, const double A)
{
    double d = 0.5 * std::pow(A * 1.0,(1.0/3.0)) * 1.25;
    int ntot = u.size();
    int nstu = ntot/2;
    int nst, nps;

    for (int spin = 1; spin<3; spin++) {
        if( spin == 2) {nst = nstu; nps = ntot;}
        else{ nst = 0;  nps = nstu;}
        for (int ka = 0; ka<nps; ka++) {
            for (int k=0; k<nps; k++) {
                for (int j=0; j<nps; j++) {
                    for (int i=0; i<nps; i++) {
                        if (ka == (k + j + i)) {
                            if (nst<nps) {u[nst] = comp_factoryT(world).functor(comp_functorT(new HOm(i,j,k,d)));}
                            nst++;
                        }
                    }
                }
            }
        }
    }
    world.gop.fence();
}


// Make vector of plane waves wavefunctions
void make_Fermi(World& world, comp_vecfuncT& u,
                              const double A,
                              const double Z,
                              const double L,
                              const double prec)
{
    double res = FunctionDefaults<3>::get_thresh()*prec;
    int ntot   = u.size();
    int nstu   = ntot/2;
    int nst, nstold, nps, ii, jj, kk;

    for (int spin = 1; spin<3; spin++) {
        if (spin == 2) {nst = nstu; nps = ntot;}
        else {nst = 0; nps = nstu;}
        for (int ka=1; ka<nps; ka++) {
            if (nst < nps) {
                for (int k=nps-1; k>=0; k--) {
                for (int j=nps-1; j>=0; j--) {
                for (int i=nps-1; i>=0; i--) {

                if (ka == k * k + j * j + i * i) {
                    for (int pk = 0; pk < 2; pk++) {
                    for (int pj = 0; pj < 2; pj++) {
                    for (int pi = 0; pi < 2; pi++) {

                    if (nst<nps) {
                        ii = pow(-1.0, pi);
                        jj = pow(-1.0, pj);
                        kk = pow(-1.0, pk);
                        u[nst] = comp_factoryT(world).functor(comp_functorT(new Fermi(i, j, k, ii, jj, kk, L)));
                    }
                    nstold = nst;
                    nst++;

                    if ((i == 0 && j == 0 && k == 0)) {nst = nstold;}
                    else if (i == 0 && pi > 0) {nst = nstold;}
                    else if (j == 0 && pj > 0) {nst = nstold;}
                    else if (k == 0 && pk > 0) {nst = nstold;}
                    else if (nstold < nps) {
                        u[nstold].scale(1.0/u[nstold].norm2());
                        u[nstold].truncate(res);
                        if(world.rank() == 0) {print("State = ", nstold);}
                    }
                    }}}
                }
                }}}
            }
        }
    }

    if (u.size() > Z) {
        u[ntot-1] = comp_factoryT(world).functor(comp_functorT(new MD(-8.0,0.0,0.0, L)));
        u[ntot-2] = comp_factoryT(world).functor(comp_functorT(new MD(-6.0,0.0,0.0, L)));
        u[ntot-3] = comp_factoryT(world).functor(comp_functorT(new MD(-3.0,0.0,0.0, L)));
        u[ntot-4] = comp_factoryT(world).functor(comp_functorT(new MD(0.0,0.0,0.0, L)));
        u[ntot-5] = comp_factoryT(world).functor(comp_functorT(new MD(3.0,0.0,0.0, L)));
        u[ntot-6] = comp_factoryT(world).functor(comp_functorT(new MD(6.0,0.0,0.0, L)));
    }
    world.gop.fence();
}





// Normalizes two vectors of complex functions
void normalize_2v(World& world, comp_vecfuncT& psi_n, comp_vecfuncT& psi_p)
{
    std::vector<double> nnorm = norm2s(world, psi_n);
    std::vector<double> pnorm = norm2s(world, psi_p);

    world.gop.fence();
    for (unsigned int i = 0; i < psi_n.size(); i++) {
        psi_n[i].scale(double_complex(1.0/nnorm[i], 0.0), false);
    }
    for (unsigned int i = 0; i < psi_p.size(); i++) {
        psi_p[i].scale(double_complex(1.0/pnorm[i], 0.0), false);
    }
}


// Normalizes one vector of complex functions with spin up and down
void normalize_ud(World& world, comp_vecfuncT& psi_qu, comp_vecfuncT& psi_qd)
{
    std::vector<double> unorm = norm2s(world,psi_qu);
    std::vector<double> vnorm = norm2s(world,psi_qd);
    std::vector<double> normu(psi_qu.size());

    int nst = psi_qu.size();

    for (int i = 0; i < nst; i++) {
        normu[i] = std::sqrt(unorm[i] * unorm[i] + vnorm[i] * vnorm[i]);
    }
    world.gop.fence();

    for (int i = 0; i < nst; i++) {
        psi_qu[i].scale(double_complex(1.0/normu[i], 0.0), false);
    }
    world.gop.fence();

    for (int i = 0; i < nst; i++) {
        psi_qd[i].scale(double_complex(1.0/normu[i], 0.0), false);
    }
}


// Truncates two vectors of complex functions
void truncate2(World& world, comp_vecfuncT& psi_q1,
                             comp_vecfuncT& psi_q2,
                             const double prec)
{
    double res = FunctionDefaults<3>::get_thresh() * prec;
    truncate(world, psi_q1, res);
    truncate(world, psi_q2, res);
}


// Truncates one vector of complex functions
void truncate1(World& world, comp_vecfuncT& psi_q, const double prec)
{
    double res = FunctionDefaults<3>::get_thresh() * prec;
    truncate(world, psi_q, res);
}


// Splits one vector of complex functions into two for spin up and 
// spin down
void spin_split(World& world, const comp_vecfuncT& psi_q, 
                              comp_vecfuncT& psi_qu,
                              comp_vecfuncT& psi_qd)
{
    reconstruct(world, psi_q); 
    unsigned int nv  = psi_q.size();
    unsigned int nvu = nv/2.0;

    for (unsigned int i = 0; i < psi_q.size(); i++) {
        if (i < nvu) {psi_qu[i] = copy(psi_q[i]);}
        else         {psi_qd[i] = copy(psi_q[i]);}
    }
    world.gop.fence();
}
