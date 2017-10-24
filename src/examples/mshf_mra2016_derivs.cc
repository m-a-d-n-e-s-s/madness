// date 08/10/16, modified by IS
// // date 12/04/15, modified by IS
// // date 11/18/2015

/*
#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/vmra.h>
#include <mra/operator.h>
#include <mra/lbdeux.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
*/

// #define WORLD_INSTANTIATE_STATIC_TEMPLATES
#define NO_GENTENSOR
#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/mra/operator.h>
#include <madness/constants.h>
#include <madness/mra/lbdeux.h>
#include <iostream>
#include <iomanip>


/*
#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/mra/operator.h>
#include <madness/mra/lbdeux.h>
#include <madness/misc/ran.h>

#include <iostream>
#include <iostream>
#include <iomanip>
*/

// #include <mpreal.h>

const char* DATA_PATH =    "/gpfs/projects/rjh/mad-der/src/madness/mra/deriv-bsp.k=m+1.n=m+1";
//const char* DATA_PATH  = "/gpfs/projects/rjh/mad-der/src/madness/mra/ph-spline-deriv.txt";
//const char* DATA_PATH=   "/gpfs/projects/rjh/mad-der/src/madness/mra/prolates-joel";

using namespace madness;

// Definitions of functions, function-vectors, operators

typedef Vector<double,3> coordT;
typedef std::shared_ptr< FunctionFunctorInterface<double,3> > real_functorT;
typedef std::shared_ptr< FunctionFunctorInterface<double_complex,3> > comp_functorT;
typedef Function<double,3> real_functionT;
typedef Function<double_complex,3> comp_functionT;
typedef std::vector<real_functionT> real_vecfuncT;
typedef std::vector<comp_functionT> comp_vecfuncT;
typedef SeparatedConvolution<double,3> operatorT;
typedef std::shared_ptr<operatorT> poperatorT;
typedef Tensor<double> real_tensorT;
typedef Tensor<double_complex> comp_tensorT;
typedef FunctionFactory<double,3> real_factoryT;
typedef FunctionFactory<double_complex,3> comp_factoryT;
typedef std::shared_ptr< WorldDCPmapInterface< Key<3> > > pmapT;

double ttt, sss;
#define START_TIMER world.gop.fence(); ttt=wall_time(); sss=cpu_time()
#define END_TIMER(msg) ttt=wall_time()-ttt; sss=cpu_time()-sss; if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt)


// Input Parameters 
//-----------------

const int A          = 16;                                 // Mass number, n+p
const int Z          = 8;                                  // Charge number, p

const double Box     = 24.0;                               // Total box size
const double L       = Box/2.0;                             // 1/2 of box size (that is used in computation)

const int initial    = 2;                                   // Single-particle state initialization
                                                            // (1=MD input; 2=Harmonic oscillator; 3=Plane waves)
static const int HOshell = 30;                              // the HO basis shells can be up to 100

const int periodic   = 1;                                   // Periodic boundary conditions (1=yes, 0=no)
const int jellium    = 1;                                   // Jellium background (1=yes, 0=no)

const int spinorbit  = 1;                                   // Spin-orbit potential and density (1=yes, else no)
const int meff       = 1;                                   // Effective mass potential (1=yes, else no)

const int screening  = 0;                                   // Coulomb screening  (1=yes, 0=no)
const double screenl = 3.0;                                 // Screening length [fm]

const int avg_pot    = 0;                                   // Apply relaxation/averaging to Skyrme potential (1=yes)
const int avg_lap    = 1;                                   // Apply relaxation/averaging to laplacian of density (1=yes)
const int avg_wav    = 1;                                   // Apply relaxation/averaging to single particle states (1=yes)

const double prec    = 0.1;                                 // Additional factor for truncation. Truncation is done at thresh*prec	
const double d       = 0.5*std::pow(A*1.0,(1.0/3.0))*1.25;  // Radius of Gaussians for HO init.

const bool VTK_OUTPUT     = false;                          // no vtk
const bool PLOT_DENSITIES = true;                           // no no densities



// Directory for vts output files and checkpointing files
//------------------
//char direct[] = "/lustre/atlas/scratch/isagert/nph106/eos/GF/";
char direct[] = "./";

// Physical parameters
//------------------
const double e2        = 1.43989e0;             // square of elementary charge [MeV*fm]
const double_complex I = double_complex(0.0, 1.0);


// Skyrme Forces
//------------------------

// SLy4
/*
static const double t0    = -2488.91300, t1 = 486.81800, t2 = -546.39500, t3 = 13777.00000;
static const double x0    =  0.83400, x1 = -0.34400, x2 = -1.00000, x3 = 1.35400;
static const double alpha = 0.16667;
static const double t4    = 123.0000;
static const double bp4   = 61.50000;
const double kf = 20.73553;
*/

// SV-bas
static const double t0 = -1879.64001, t1 = 313.74933, t2 = 112.67627, t3 = 12527.38965, t4 = 124.63330;
static const double x0 = 0.25855, x1 = -0.38169, x2 = -2.82364, x3 = 0.12323, bp4 = 34.11170;
static const double alpha = 0.30;
static const double kfn = (20.72126e0 + 20.74982e0)/2.0;
const double kfp = (20.72126e0 + 20.74982e0)/2.0;

//-------------------------



static const double b0  = t0*(1.e0 + 0.5e0*x0);
static const double b1  = (t1 + 0.5e0*x1*t1 + t2 + 0.5e0*x2*t2)/4.e0;
static const double b2  = (3.e0*t1*(1.e0 + 0.5e0*x1) - t2*(1.e0 + 0.5e0*x2))/8.e0;
static const double b3  = t3*(1.e0 + 0.5e0*x3)/4.e0;
static const double b4  = 0.5e0*t4;
static const double bp0 = t0*(0.5e0 + x0);
static const double bp1 = (t1*(0.5e0 + x1) - t2*(0.5e0 + x2))/4.e0;
static const double bp2 = (3.e0*t1*(0.5e0 + x1) + t2*(0.5e0 + x2))/8.e0;
static const double bp3 = t3*(0.5e0 + x3)/4.e0;

static const double PI = 3.141592653589793238;



//////////////////////////////////////////////////////////////////////
//////////////////// DATA LOAD BALANCING /////////////////////////////
//////////////////// see dataloadbal.cc for details //////////////////
//////////////////////////////////////////////////////////////////////
// Load Balance
template <typename T, int NDIM>
struct LBCost 
{
    double leaf_value;
    double parent_value;
    LBCost(double leaf_value=1.0, double parent_value=0.0)
        : leaf_value(leaf_value), parent_value(parent_value){}
    double operator()(const Key<NDIM>& key, const FunctionNode<T,NDIM>& node) const 
    {
        if (key.level() <= 1){
        return 100.0*(leaf_value+parent_value);
        }   
        else if (node.is_leaf()) { return leaf_value; }
        else { return parent_value; }
    }
};


// Load balance of 3 real functions
void loadbalance_xxp(World& world, real_functionT& p_den, real_functionT& n_den, real_functionT& den3 )
{
    LoadBalanceDeux<3> lb(world);
    lb.add_tree(p_den, LBCost<double,3>(1.0,1.0), false);
    lb.add_tree(n_den, LBCost<double,3>(1.0,1.0), false);
    lb.add_tree(den3,  LBCost<double,3>(1.0,1.0), false);
    world.gop.fence();
    FunctionDefaults<3>::redistribute(world, lb.load_balance(6.0));
    world.gop.fence();
}


void loadbalance_v1(World& world, comp_vecfuncT& u1)
{
    LoadBalanceDeux<3> lb(world);
    for (unsigned int i=0; i<u1.size(); i++)
    {
        lb.add_tree(u1[i], LBCost<double_complex,3>(1.0,1.0), false);
    }
    world.gop.fence();  // FENCE
    FunctionDefaults<3>::redistribute(world, lb.load_balance(6.0));
    world.gop.fence();
}    



// Load-balance for two function-vectors
void loadbalance_v2(World& world, comp_vecfuncT& u1, comp_vecfuncT& u2)
{
    LoadBalanceDeux<3> lb(world);
    for (unsigned int i=0; i<u1.size(); i++)
    {
        lb.add_tree(u1[i], LBCost<double_complex,3>(1.0,1.0), false);
    }
    for (unsigned int i=0; i<u2.size(); i++)
    {
        lb.add_tree(u2[i], LBCost<double_complex,3>(1.0,1.0), false);
    }
    world.gop.fence();  // FENCE
    FunctionDefaults<3>::redistribute(world, lb.load_balance(6.0));
    world.gop.fence();
}



// Load-balance for three function-vectors
void loadbalance_v3(World& world, comp_vecfuncT& u1, comp_vecfuncT& u2,
comp_vecfuncT& u3)
{
    LoadBalanceDeux<3> lb(world);
    for (unsigned int i=0; i<u1.size(); i++)
    {
        lb.add_tree(u1[i], LBCost<double_complex,3>(1.0,1.0), false);
    }
    for (unsigned int i=0; i<u2.size(); i++)
    {
        lb.add_tree(u2[i], LBCost<double_complex,3>(1.0,1.0), false);
    }
    for (unsigned int i=0; i<u3.size(); i++)
    {
        lb.add_tree(u3[i], LBCost<double_complex,3>(1.0,1.0), false);
    }
    world.gop.fence();
    FunctionDefaults<3>::redistribute(world, lb.load_balance(6.0));
    world.gop.fence();
}


// Load-balance for four function-vectors
void loadbalance_v4(World& world, comp_vecfuncT& u1, comp_vecfuncT& u2,
comp_vecfuncT& u3, comp_vecfuncT& u4)
{
    LoadBalanceDeux<3> lb(world);
    for (unsigned int i=0; i<u1.size(); i++)
    {
        lb.add_tree(u1[i], LBCost<double_complex,3>(1.0,1.0), false);
    }
    for (unsigned int i=0; i<u2.size(); i++)
    {
        lb.add_tree(u2[i], LBCost<double_complex,3>(1.0,1.0), false);
    }
    for (unsigned int i=0; i<u3.size(); i++)
    {
        lb.add_tree(u3[i], LBCost<double_complex,3>(1.0,1.0), false);
    }
    for (unsigned int i=0; i<u4.size(); i++)
    {
        lb.add_tree(u4[i], LBCost<double_complex,3>(1.0,1.0), false);
    }
    world.gop.fence();
    FunctionDefaults<3>::redistribute(world, lb.load_balance(6.0));
    world.gop.fence();
}


void loadbalance_v6(World& world,
        comp_vecfuncT& u1,
        comp_vecfuncT& u2,
        comp_vecfuncT& u3,
        comp_vecfuncT& v1,
        comp_vecfuncT& v2,
        comp_vecfuncT& v3)
{
    LoadBalanceDeux<3> lb(world);
    for (unsigned int i=0; i<u1.size(); i++)
    {
        lb.add_tree(u1[i], LBCost<double_complex,3>(1.0,1.0), false);
    }
    for (unsigned int i=0; i<u2.size(); i++)
    {
        lb.add_tree(u2[i], LBCost<double_complex,3>(1.0,1.0), false);
    }
    for (unsigned int i=0; i<u3.size(); i++)
    {
        lb.add_tree(u3[i], LBCost<double_complex,3>(1.0,1.0), false);
    }
    for (unsigned int i=0; i<v1.size(); i++)
    {
        lb.add_tree(v1[i], LBCost<double_complex,3>(1.0,1.0), false);
    }
    for (unsigned int i=0; i<v2.size(); i++)
    {
        lb.add_tree(v2[i], LBCost<double_complex,3>(1.0,1.0), false);
    }
    for (unsigned int i=0; i<v3.size(); i++)
    {
        lb.add_tree(v3[i], LBCost<double_complex,3>(1.0,1.0), false);
    }
    world.gop.fence();
    FunctionDefaults<3>::redistribute(world, lb.load_balance(6.0));
    world.gop.fence();
}




////////////////////////////////////////////////////////////////////////
///////////////////////// INITIALIZATION ///////////////////////////////
////////////////////////////////////////////////////////////////////////
// Coordinates



double coord_x(const coordT& r) { const double x=r[0]; return x;}
double coord_y(const coordT& r) { const double y=r[1]; return y;}
double coord_z(const coordT& r) { const double z=r[2]; return z;}


double myf(const coord_3d& r) 
{
    //return exp(-1.0 * pow( 0.4*sqrt( r[0]*r[0] + r[1]*r[1] + r[2]*r[2] ) , 10.0 ) / (2.0*pow(3.0, 10.0))) ; 
    double ra = 1.25*pow(48.0,0.33);
    double x  = r[0] - 0.0*(1.0/3.0);
    double y  = r[1] - 0.0*(1.0/3.0);
    double z  = r[2] - 0.0*(1.0/3.0);
    return 65.0/(1.0 + exp( ( sqrt( x*x + y*y + z*z ) - ra) /0.55));
}



// Harmonic oscillator wavefunction
struct HO : FunctionFunctorInterface<double_complex,3>
{
    const int nx, ny, nz;
    HO(int nx, int ny, int nz) : nx(nx), ny(ny), nz(nz) {}
    double_complex operator()(const coordT& r) const 
    {
        double rsq = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        double psi = exp( (-1.0/(2.0*d*d))*( rsq*rsq ) );
        if (nx) for(int inx=0; inx<nx; inx++) {psi *= (1.0/d)*r[0];}
        if (ny) for(int iny=0; iny<ny; iny++) {psi *= (1.0/d)*r[1];}
        if (nz) for(int inz=0; inz<nz; inz++) {psi *= (1.0/d)*r[2];}
        return double_complex(psi,0.0);
     }
};


// Modified harmonic oscillator wavefunction
// Allows to have different gaussian widths in x,y,z direction 
struct HOm : FunctionFunctorInterface<double_complex,3> 
{
    const int nx, ny, nz;
    HOm(int nx, int ny, int nz) : nx(nx), ny(ny), nz(nz) {}
    double_complex operator()(const coordT& r) const {
    const double dx = d;  // width in x
    const double dy = dx; // width in y
    const double dz = dx; // width in z
    double rsq = std::sqrt( (r[0]*r[0]/(dx*dx)) + (r[1]*r[1]/(dy*dy)) + (r[2]*r[2]/(dz*dz)) );
    double psi = std::exp( (-1.0/2.0)*( rsq*rsq ) );

    if (nx) for(int inx=0; inx<nx; inx++) {psi *= r[0];}
    if (ny) for(int iny=0; iny<ny; iny++) {psi *= r[1];}
    if (nz) for(int inz=0; inz<nz; inz++) {psi *= r[2];}
    return double_complex(psi, 0.0);
    }
};


// Molecular dynamics input
// Uses nucleon coordinate from e.g. IUMD, transforms the 
// coordinates to a simulation space from -L,L , and 
// creates 27 mirror images due periodic boundary conditions.
// Folds Gaussian around nucleon corrdinate
struct MD : FunctionFunctorInterface<double_complex,3> 
{
    const double nx, ny, nz;
    MD(double nx, double ny, double nz) : nx(nx), ny(ny), nz(nz) {}
    double_complex operator()(const coordT& r) const 
    {
        const double dmd = 3.0; // width of Gaussian
        double psi = 0.0;

        for(int k=-1; k < 2; k++)
        {
            for(int j=-1; j < 2; j++) 
            {
                for(int i=-1; i < 2; i++)
                {          
                    const double xs = ( r[0] - ( nx + (2.0*L*i-0.0)) )/dmd;
                    const double ys = ( r[1] - ( ny + (2.0*L*j-0.0)) )/dmd;
                    const double zs = ( r[2] - ( nz + (2.0*L*k-0.0)) )/dmd;
                    psi += exp( (-1.0/2.0)*( xs*xs + ys*ys + zs*zs ) );
                }
            }
        }
        return double_complex(psi,0.0);
     }
};



// Makes fermi (plane wave) wavefunction
struct Fermi : FunctionFunctorInterface<double_complex,3>
{
  const int nx, ny, nz, px, py, pz;
  Fermi(int nx, int ny, int nz, int px, int py, int pz) : nx(nx), ny(ny), nz(nz), px(px), py(py), pz(pz) {}
  double_complex operator()(const coordT& r) const {
    double kx = px*nx*2.0*M_PI/(2.0*L);
    double_complex psi  = exp( I*kx*r[0] );
    double ky = py*ny*2.0*M_PI/(2.0*L);
    psi *= exp( I*ky*r[1] );
    double kz = pz*nz*2.0*M_PI/(2.0*L);
    psi *= exp( I*kz*r[2] );
    return psi;
  }
};


double hermite(const int n, const double x) 
{
    //  double HO[n+1], CHE[n+1];
    double ho[n+1], che[n+1]; // should set it somewhere
    ho[0]=1.0; ho[1]=2.0*x;
    che[0]=1.0; che[1]=1.0;
    // printf("hermit %d %f %f %f \n", n, x, ho[0], ho[1]);
    for(int i=2; i<(n+1); i++)
    {
        ho[i] = (2.0*x*ho[i-1]-2.0*(i-1.0)*ho[i-2]);
        che[i] = i*che[i-1];
        // printf("hermit i %d %f %f %f \n", i, ho[i], che[i]);
    }
    double r = ho[n]/sqrt(sqrt(PI)*pow(2.0,n)*che[n]);
    return r;
}

struct Hobasis : FunctionFunctorInterface<double_complex,3>
{
    const int nx, ny, nz;  // Ho eigenvalues

    Hobasis(int nx, int ny, int nz) 
        : nx(nx), ny(ny), nz(nz) {}

    double_complex operator()(const coordT& r) const
    {
        double bx = sqrt( 939.56*41.0/(197.3*197.3*pow(A, 1./3.)));
        double sbx = sqrt(bx), x = bx*r[0], y=bx*r[1], z=bx*r[2];
        double rp = sbx*hermite(nx,x)*exp(-0.5*x*x)
                    *sbx*hermite(ny,y)*exp(-0.5*y*y)*sbx*hermite(nz,z/1.4)*exp(-0.5*z*z/(1.4*1.4));
        return double_complex(rp,0.0);
    }
};



void make_ho2(World& world, comp_vecfuncT& u)
{
    double res = FunctionDefaults<3>::get_thresh();
    res *= prec; // ?

    int nst=0, nps;
    int ntot = (int) u.size();
    int nstu = ntot/2; // what happens if ntot = odd?

    for(int i=0; i <= HOshell; i++)
    {
        for(int k1=0; k1 <= i; k1++)
        {
            for(int k2=0; k2 <= i-k1; k2++)
            {
                int k3=i-k1-k2;
                if ((int) nst==u.size()) break;
                u[nst] =comp_functionT(comp_factoryT(world).functor(comp_functorT( new Hobasis(k1, k2, k3))).fence(false));  //
                if(world.rank() == 0)
                printf(" HO basis2 %d \n", nst);
                nst++;
            }
        }
    }
    world.gop.fence();  // FENCE
}



/*
// Make wavefunctions: Harmonic Oscillator states
void make_ho(World& world, comp_vecfuncT& u){
    double res = FunctionDefaults<3>::get_thresh()*prec;
    //    res *= prec;  
    
    int ntot = u.size();
    int nstu = ntot/2;
    int nst, nps;
    for (int spin = 1; spin<3; spin++)
      {
        if( spin == 2){ nst = nstu; nps = ntot;}
        else{ nst = 0;  nps = nstu;}        
        for (int ka = 0; ka<nps; ka++)
        {
            for(int k=0; k<nps; k++) {
            for(int j=0; j<nps; j++) {
            for(int i=0; i<nps; i++) {
            if( ka==(k+j+i) ) {
	      if( nst<nps) {
		//		u[nst] = comp_factoryT(world).functor(comp_functorT(new HOm(i,j,k))).nofence();
		//		u[nst] = comp_factoryT(world).functor(comp_functorT(new HOm(i,j,k)));
		u[nst] =comp_functionT(comp_factoryT(world).functor(comp_functorT( new HOm(i, j, k))).fence(false));  // s   
		if(world.rank() == 0)
		  print(" U:HO-0 basis", nst);
	      }
	      nst++;
	    } } } }
        }
    }
    world.gop.fence();  // FENCE
}
*/



//                                                              IS
void make_ho(World& world, comp_vecfuncT& u)
{
    double res = FunctionDefaults<3>::get_thresh();
    res *= prec;
    
    int ntot = u.size();
    int nstu = ntot/2;
    int nst, nps;
    for (int spin = 1; spin<3; spin++)
    {
        if( spin == 2){ nst = nstu; nps = ntot;}
        else{ nst = 0;  nps = nstu;}
        for (int ka = 0; ka<nps; ka++)
        {
            for(int k=0; k<nps; k++) {
            for(int j=0; j<nps; j++) {
            for(int i=0; i<nps; i++) {
            if( ka==(k+j+i) ) {
            if( nst<nps) {
                u[nst] = comp_factoryT(world).functor(comp_functorT(new HOm(i,j,k)));
        }
            nst++;
            } } } }
        }
    }

    world.gop.fence();  // FENCE
}
//                                                              IS





// Make wavefunctions: Fermi plane waes
void make_fermi(World& world, comp_vecfuncT& u)
{
    double res = FunctionDefaults<3>::get_thresh()*prec;
    int ntot   = u.size();
    int nstu   = ntot/2;
    int nst, nstold, nps, ii, jj, kk;

    for(int spin = 1; spin<3; spin++)
    {
        if( spin == 2){ nst = nstu; nps = ntot;}
        else{ nst = 0; nps = nstu;}
        for (int ka=1; ka<nps; ka++)
        {
            if( nst < nps) {
            for(int k=nps-1; k>=0; k--) {
            for(int j=nps-1; j>=0; j--) {
            for(int i=nps-1; i>=0; i--) {
            if(ka == k*k+j*j+i*i ) {
                for(int pk = 0; pk < 2; pk++) {
                for(int pj = 0; pj < 2; pj++) {
                for(int pi = 0; pi < 2; pi++) {
                if( nst<nps) {
		  ii = pow( -1.0, pi);
		  jj = pow( -1.0, pj);
		  kk = pow( -1.0, pk);
		  //		  u[nst] = comp_factoryT(world).functor( comp_functorT( new Fermi( i,j,k,ii,jj,kk ) ) ).thresh(.00001).nofence(); // should be nofence
		  //		  u[nst] = comp_factoryT(world).functor( comp_functorT( new Fermi( i,j,k,ii,jj,kk ) ) ).thresh(.00001); // should be nofence
		  u[nst] =comp_functionT(comp_factoryT(world).functor(comp_functorT( new Fermi(i, j, k, ii, jj, kk))).fence(false));  // s   

		  if(world.rank() == 0){
		    print(" State:", nst, "pi", pi, "pj", pj, "pk", pk,"i", i,"j", j);
		    print("k", k, "ka", ka,"spin", spin);
		  }

		}
		nstold = nst;
		nst++;

		if( (i == 0 && j == 0 && k == 0) ) { nst = nstold; }
		else if( i == 0 && pi > 0) { nst = nstold; }
		else if( j == 0 && pj > 0) { nst = nstold; }
		else if( k == 0 && pk > 0) { nst = nstold; }
		else if( nstold < nps ) 
		{
		    u[nstold].scale(1.0/u[nstold].norm2()); 
		    u[nstold].truncate(res);
		}
		//		if(world.rank() == 0){
		//print(" State:", nst, "pi", pi, "pj", pj, "pk", pk,"i", i,"j", j);
		//		  print("k", k, "ka", ka,"spin", spin);
		//}
		}
		}
		}
	    } } } } 
            }
        }
    }
    world.gop.fence();  //  FENCE

    //Put some disturbances into the system. Change some high-energy plane waves into 
    // Gaussians
    if( u.size() > Z)
    {
	// Sphere
	//u[ntot-1] = comp_factoryT(world).functor(comp_functorT(new MD(0.25*L,0.25*L,0.0)));	
	//u[ntot-2] = comp_factoryT(world).functor(comp_functorT(new MD(0.25*L,-0.25*L,0.0)));
	//u[ntot-3] = comp_factoryT(world).functor(comp_functorT(new MD(0.0,0.0,0.25*L)));
	//u[ntot-4] = comp_factoryT(world).functor(comp_functorT(new MD(-0.25*L,0.25*L,0.0)));
	//u[ntot-5] = comp_factoryT(world).functor(comp_functorT(new MD(-0.25*L,-0.25*L,0.0)));
	//u[ntot-6] = comp_factoryT(world).functor(comp_functorT(new MD(0.0,0.0,-0.25*L)));

	// Rod
      u[ntot-1] = comp_factoryT(world).functor(comp_functorT(new MD(-8.0,0.0,0.0))).nofence();
      u[ntot-2] = comp_factoryT(world).functor(comp_functorT(new MD(-6.0,0.0,0.0))).nofence();
      u[ntot-3] = comp_factoryT(world).functor(comp_functorT(new MD(-3.0,0.0,0.0))).nofence();
      u[ntot-4] = comp_factoryT(world).functor(comp_functorT(new MD(0.0,0.0,0.0))).nofence();
      u[ntot-5] = comp_factoryT(world).functor(comp_functorT(new MD(3.0,0.0,0.0))).nofence();
      u[ntot-6] = comp_factoryT(world).functor(comp_functorT(new MD(6.0,0.0,0.0))).nofence();
    }                            
    else
    {
        //u[0]    = comp_factoryT(world).functor(comp_functorT(new MD(0.75*L,-0.1*L,-0.05*L)));
        //u[nstu] = comp_factoryT(world).functor(comp_functorT(new MD(-0.65*L,-0.2*L,0.55*L)));    
    }                                                    
    world.gop.fence();  //  FENCE                                    
}





/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// NORMALIZE & TRUNCATE ///////////////////////////
/////////////////////////////////////////////////////////////////////////////////



// Normalizes two function-vectors, only used in the initialization
void normalize_2v(World& world, comp_vecfuncT& n, comp_vecfuncT& p)
{
    std::vector<double> nnorm = norm2s(world,n);
    std::vector<double> pnorm = norm2s(world,p); 
    world.gop.fence();  // FENCE // gif
    for (unsigned int i=0; i<n.size(); i++){ n[i].scale(double_complex(1.0/nnorm[i],0.0),false); }
    for (unsigned int i=0; i<p.size(); i++){ p[i].scale(double_complex(1.0/pnorm[i],0.0),false); }
    world.gop.fence();  // FENCE
}



// Normalizes one function-vector
void normalize_ud(World& world, comp_vecfuncT& u, comp_vecfuncT& v) 
{
    std::vector<double> unorm = norm2s(world,u);
    std::vector<double> vnorm = norm2s(world,v);
    std::vector<double> normu(u.size());
    for (unsigned int i=0; i<u.size(); i++) 
    {
        normu[i] = std::sqrt(unorm[i]*unorm[i] + vnorm[i]*vnorm[i]);
    }
    world.gop.fence();  // FENCE // gif
    for (unsigned int i=0; i<u.size(); i++)
    {
        u[i].scale(double_complex(1.0/normu[i],0.), false);
        v[i].scale(double_complex(1.0/normu[i],0.), false);
    }
    world.gop.fence();  // FENCE // gif
}



// Truncates two function-vectors
void truncate2(World& world, comp_vecfuncT& n, comp_vecfuncT& p)
{
    double res = FunctionDefaults<3>::get_thresh() * prec;
    truncate(world, n, res);
    truncate(world, p, res);
}


// Truncates one function-vector
void truncate1(World& world, comp_vecfuncT& u)
{
    double res = FunctionDefaults<3>::get_thresh() * prec;
    truncate(world, u, res);
}


// Splits one function-vector in two for spin up and spin down
void spin_split(World& world, const comp_vecfuncT& u, comp_vecfuncT& uu,
comp_vecfuncT& ud)
{
    reconstruct(world, u);  // Not in gf version
    int nv = u.size();
    int nvu = nv/2.0;
    for(unsigned int i=0; i < u.size(); i++)
    {
        if( i < nvu){ uu[i] = copy(u[i]); }
        else{ ud[i] = copy(u[i]); }
    }
    world.gop.fence();   // FENCE
}



// Spin split
void spin_split2(World& world, const comp_vecfuncT& u, comp_vecfuncT& uu,
comp_vecfuncT& ud)
{
    reconstruct(world, u);  // Not in gf version
    for(unsigned int i=0; i < u.size(); i++)
    {
        if( i==0 || i%2 == 0 ) { uu[i] = copy(u[i]); }
        else { ud[i] = copy(u[i]); }
    }
    world.gop.fence();   // FENCE
}



///////////////////////////////////////////////////////////////////////////
//////////////////////// POTENTIALS AND DENSITIES /////////////////////////
///////////////////////////////////////////////////////////////////////////


// Calculates abs of a real function
struct abs_func
{
    void operator()(const Key<3>& key, real_tensor ABS,
                    const real_tensor& rho,
                    const real_tensor& rho_u)
    const
    {
        ITERATOR(ABS, double d = rho(IND);
                      double p = rho(IND);
                      if( d < 0.0 ) d *= -1.0;
                      ABS(IND) = std::abs(d);
                );
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};




// Calculates abs of a real function
struct ln_func
{
    void operator()(const Key<3>& key, real_tensor LN,
                    const real_tensor& rho,
                    const real_tensor& rho_u)
    const
    {
        ITERATOR(LN,  double d = std::abs(rho(IND));
                      double p = rho(IND);
                      LN(IND) = std::log(d);
                );
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};



// Calculates abs of a real function
struct exp_func
{
    void operator()(const Key<3>& key, real_tensor EXPF,
                    const real_tensor& rho,
                    const real_tensor& rho_u)
    const
    {
        ITERATOR(EXPF, double d = rho(IND);
                       double p = rho(IND);
                       EXPF(IND) = std::exp(d);                       
                );
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};



// Calculates rho^alpha for Skyrme Energy (only used in energy calculations)
struct rho4
{
  void operator()(const Key<3>& key,
              real_tensor RHO,
              const real_tensor& rho,
              const real_tensor& rho_u)
    const
    {
    ITERATOR(RHO, double d = rho_u(IND);
                  double p = rho(IND);
                  if (p<1.e-10 || d<1.e-10) RHO(IND) = 0.0;
                  else RHO(IND) = std::pow(p,alpha);
         );
  }
  template <typename Archive>
  void serialize(Archive& ar) {}
};


// Calculates Fgamma factor for Coulomb exchange potential WITH screening
struct Fgamma
{
    void operator()(const Key<3>& key, real_tensor FGAMMA,
                    const real_tensor& rho,
                    const real_tensor& len)
    const
    {
        ITERATOR(FGAMMA, double length = len(IND);
                         double den    = rho(IND);
                         double gamma  = length/(std::pow(3.0*M_PI*M_PI*den,1.0/3.0));
                         double t1     = atan(2.0/gamma);
                         double ln1    = log(1.0 + 4.0/(gamma*gamma));
                         double g2     = gamma*gamma;
                         if (den<1.e-10) FGAMMA(IND) = 0.0;
                         else FGAMMA(IND) = 1.0 - (4.0/3.0)*gamma*t1 + 0.5*g2*ln1
                             + (1.0/6.0)*g2*(1 - 0.25*g2*ln1);
                );
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};



//Calculates rho^{alpha+1) for Skyrme potential
struct rho1
{
    void operator()(const Key<3>& key, 
                    real_tensor RHO31,
                    const real_tensor& rho,
                    const real_tensor& rho_u)
    const
    {
        ITERATOR(RHO31, double d = rho_u(IND);
                        double p = rho(IND);
                        if (p<=0.0) RHO31(IND) = 0.0;
                        else RHO31(IND) = std::pow(p,alpha+1.0);
                 );
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};



//Calculates rho_q*rho^{alpha) for Skyrme potential
struct rho2
{
    void operator()(const Key<3>& key, 
                    real_tensor RHO32,
                    const real_tensor& rho,
                    const real_tensor& rho_u)
    const
    {
        ITERATOR(RHO32, double d = rho_u(IND);
                        double p = rho(IND);
                        if (p<=0.0 || d<=0.0) RHO32(IND) = 0.0;
                        else RHO32(IND) = d*std::pow(p,alpha);
                );
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};



// Calculates rho_q*rho_q*rho^{alpha-1) for Skyrme potential
struct rho3
{
    void operator()(const Key<3>& key, 
                    real_tensor RHO33,
                    const real_tensor& rho,
                    const real_tensor& rho_u)
    const
    {
        ITERATOR(RHO33, double d = rho_u(IND);
                        double p = rho(IND);
                        if (p<=0.0 || d<=0.0) RHO33(IND) = 0.0;
                        else RHO33(IND) = d*d*std::pow(p,alpha-1.0);
                );
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};



//Calculates Coulomb exchange potential via Slater approximation
struct uex
{
    void operator()(const Key<3>& key, 
                    real_tensor C_ex,
                    const real_tensor& rho,
                    const real_tensor& rho_u)
    const
    {
        ITERATOR(C_ex, double d = rho_u(IND);
                       double p = rho(IND);
                       double coeff = 1.0/3.0;
                       double pre = (-1.0)*pow(3.0/M_PI,1.0/3.0);
                       if (p<=0.0 || d<=0.0) C_ex(IND) = 0.0;
                       else C_ex(IND) = pre*p*std::pow(p, coeff - 2.0);
                );
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};



//Calculates laplacian of rho from single particle states
real_functionT laplacian2(World& world, const comp_vecfuncT& dux, const comp_vecfuncT& ddx,
const comp_vecfuncT& duy, const comp_vecfuncT& ddy, const comp_vecfuncT& duz, const comp_vecfuncT& ddz, 
const comp_vecfuncT& uu, const comp_vecfuncT& ud)
{
    double ntol = FunctionDefaults<3>::get_thresh() * prec;     // IS
    comp_functionT lapl = comp_factoryT(world);
    complex_derivative_3d Dx(world, 0);

    // Bryan's edit
    // Change derivative type
    Dx.read_from_file(DATA_PATH);
 
    comp_vecfuncT duc = conj(world, dux);
    comp_vecfuncT ddc = conj(world, ddx);

    comp_vecfuncT dduc = apply(world, Dx, duc);
    truncate(world, dduc, ntol);                                // IS (change truncation with thresh to truncation with ntol)
    compress(world, dduc);
    comp_vecfuncT ddux = apply(world, Dx, dux);
    truncate(world, ddux, ntol);                                // IS
    compress(world, ddux);
    comp_vecfuncT dddx = apply(world, Dx, ddx);
    truncate(world, dddx, ntol);                                // IS
    compress(world, dddx);
    comp_vecfuncT dddc = apply(world, Dx, ddc);
    truncate(world, dddc, ntol);                                // IS
    compress(world, dddc);


    comp_vecfuncT term  = mul(world, dduc, uu);
    truncate(world, term, ntol);                                // IS
    compress(world, term);
    compress(world, duc);
    compress(world, ddc);

    gaxpy(world, 1.0, term, 1.0, mul(world, conj(world, uu), ddux));
    gaxpy(world, 1.0, term, 2.0, mul(world, duc, dux));
    gaxpy(world, 1.0, term, 1.0, mul(world, dddc, ud));
    gaxpy(world, 1.0, term, 1.0, mul(world, conj(world, ud), dddx));
    gaxpy(world, 1.0, term, 2.0, mul(world, ddc, ddx));
    truncate(world, term, ntol);

    world.gop.fence();
    dduc.clear();
    ddux.clear();
    dddx.clear();
    dddc.clear();
    world.gop.fence(); // clearing at this sync

    truncate(world, term, ntol);
    compress(world, term);

    complex_derivative_3d Dy(world, 1);

    // Bryan's edit
    // Change derivative type
    Dy.read_from_file(DATA_PATH);

    duc = conj(world, duy);
    ddc = conj(world, ddy);

    gaxpy(world, 1.0, term, 1.0, mul(world, apply(world, Dy, duc), uu));
    gaxpy(world, 1.0, term, 1.0, mul(world, conj(world, uu), apply(world, Dy, duy)));
    gaxpy(world, 1.0, term, 2.0, mul(world, duc, duy));
    gaxpy(world, 1.0, term, 1.0, mul(world, apply(world, Dy, ddc), ud));
    gaxpy(world, 1.0, term, 1.0, mul(world, conj(world, ud), apply(world, Dy, ddy)));
    gaxpy(world, 1.0, term, 2.0, mul(world, ddc, ddy));

    truncate(world, term, ntol);
    compress(world, term);

    complex_derivative_3d Dz(world, 2);

    // Bryan's edit
    // Change derivative type
    Dz.read_from_file(DATA_PATH);
  
    duc = conj(world, duz);
    ddc = conj(world, ddz);
    compress(world, duc);
    compress(world, ddc);

    gaxpy(world, 1.0, term, 1.0, mul(world, apply(world, Dz, duc), uu));
    gaxpy(world, 1.0, term, 1.0, mul(world, conj(world, uu), apply(world, Dz, duz)));
    gaxpy(world, 1.0, term, 2.0, mul(world, duc, duz));
    gaxpy(world, 1.0, term, 1.0, mul(world, apply(world, Dz, ddc), ud));
    gaxpy(world, 1.0, term, 1.0, mul(world, conj(world, ud), apply(world, Dz, ddz)));
    gaxpy(world, 1.0, term, 2.0, mul(world, ddc, ddz));

    world.gop.fence();
    duc.clear(); ddc.clear();
    truncate(world, term, ntol);
    compress(world, term);
    world.gop.fence();

    lapl = term[0];
    lapl.compress();
    for (unsigned int i=1; i<term.size(); i++) { lapl.gaxpy(double_complex(1.0, 0.), term[i], double_complex(1.0,0.), false); }

    world.gop.fence();
    return real(lapl);
}



//Calculates laplacian of rho from single particle states                                       // IS (this laplacian is a bit faster than laplacian2 )
real_functionT laplacian1(World& world, const comp_vecfuncT& dux, const comp_vecfuncT& ddx,
const comp_vecfuncT& duy, const comp_vecfuncT& ddy, const comp_vecfuncT& duz, const comp_vecfuncT& ddz,
const comp_vecfuncT& uu, const comp_vecfuncT& ud)
{
    comp_functionT laplx = comp_factoryT(world);
    comp_functionT laply = comp_factoryT(world);
    comp_functionT laplz = comp_factoryT(world);

    complex_derivative_3d Dx(world, 0);

    // Bryan's edit
    // Change derivative type
    Dx.read_from_file(DATA_PATH);

    comp_vecfuncT termx = mul(world, conj(world, dux), uu);
    compress(world, termx);
    world.gop.fence();

    gaxpy(world, 1.0, termx, 1.0, mul(world, conj(world, ddx), ud), false);
    gaxpy(world, 1.0, termx, 1.0, mul(world, conj(world, uu), dux), false);
    gaxpy(world, 1.0, termx, 1.0, mul(world, conj(world, ud), ddx), false);
    world.gop.fence();

    comp_vecfuncT Dtermx = apply(world, Dx, termx);

    complex_derivative_3d Dy(world, 1);

    // Bryan's edit
    // Change derivative type
    Dy.read_from_file(DATA_PATH);

    comp_vecfuncT termy = mul(world, conj(world, duy), uu);
    compress(world, termy);
    world.gop.fence();

    gaxpy(world, 1.0, termy, 1.0, mul(world, conj(world, ddy), ud), false);
    gaxpy(world, 1.0, termy, 1.0, mul(world, conj(world, uu), duy), false);
    gaxpy(world, 1.0, termy, 1.0, mul(world, conj(world, ud), ddy), false);
    world.gop.fence();

    comp_vecfuncT Dtermy = apply(world, Dy, termy);

    complex_derivative_3d Dz(world, 2);

    // Bryan's edit
    // Change derivative type
    Dz.read_from_file(DATA_PATH);

    comp_vecfuncT termz = mul(world, conj(world, duz), uu);
    compress(world, termz);
    world.gop.fence();

    gaxpy(world, 1.0, termz, 1.0, mul(world, conj(world, ddz), ud), false);
    gaxpy(world, 1.0, termz, 1.0, mul(world, conj(world, uu), duz), false);
    gaxpy(world, 1.0, termz, 1.0, mul(world, conj(world, ud), ddz), false);
    world.gop.fence();

    comp_vecfuncT Dtermz = apply(world, Dz, termz);

    world.gop.fence();
    termx.clear(); termy.clear(); termz.clear();
    world.gop.fence();

    compress(world, Dtermx);
    compress(world, Dtermy);
    compress(world, Dtermz);
    world.gop.fence();

    laplx = Dtermx[0];
    laply = Dtermy[0];
    laplz = Dtermz[0];

    for (unsigned int i=1; i < Dtermx.size(); i++)
    {
        laplx.gaxpy(1.0, Dtermx[i], 1.0, false);
        laply.gaxpy(1.0, Dtermy[i], 1.0, false);
        laplz.gaxpy(1.0, Dtermz[i], 1.0, false);
    }

    world.gop.fence();
    return real(laplx) + real(laply) + real(laplz);
}




//Calculates laplacian of rho via 2nd derivative of rho
real_functionT laplacian(World& world, const real_functionT& u, const double brad)
{
    double ntol = FunctionDefaults<3>::get_thresh() * prec;

    // Bryan's edits
    //real_functionT dux = real_derivative_3d(world,0)(u);
    //real_functionT duy = real_derivative_3d(world,1)(u);
    //real_functionT duz = real_derivative_3d(world,2)(u);

    //real_functionT d2ux = real_derivative_3d(world,0)(dux); 
    //real_functionT d2uy = real_derivative_3d(world,1)(duy); 
    //real_functionT d2uz = real_derivative_3d(world,2)(duz); 

    real_derivative_3d Dx(world, 0);
    Dx.read_from_file(DATA_PATH);

    real_derivative_3d Dy(world, 1);
    Dy.read_from_file(DATA_PATH);
 
    real_derivative_3d Dz(world, 2);
    Dz.read_from_file(DATA_PATH);

    real_functionT dux = Dx(u);
    real_functionT duy = Dy(u);
    real_functionT duz = Dz(u);

    real_functionT d2ux = Dx(dux);
    real_functionT d2uy = Dy(duy);
    real_functionT d2uz = Dz(duz); 

    real_functionT lap  = d2ux + d2uy + d2uz;
    lap.truncate(ntol);

    world.gop.fence();  // FENCE
    dux.clear(); duy.clear(); duz.clear();
    d2ux.clear(); d2uy.clear(); d2uz.clear();
    world.gop.fence();

    return lap;
}



//Calculates effective-mass potential 
void Umeff(World& world, comp_vecfuncT& Vu, comp_vecfuncT& Vd, const real_functionT& B,
const comp_vecfuncT& dux, const comp_vecfuncT& duy, const comp_vecfuncT& duz,
const comp_vecfuncT& ddx, const comp_vecfuncT& ddy, const comp_vecfuncT& ddz)
{
    double ntol = FunctionDefaults<3>::get_thresh() * prec;
    const double_complex one(1.0,0.0);
    const double_complex mone(-1.0,0.0);
    reconstruct(world, Vu);
    reconstruct(world, Vd);
    B.reconstruct();
    int nv = Vu.size();

    truncate2(world, Vu, Vd);
    compress(world, Vu);
    compress(world, Vd);

    // Bryan's edits
    complex_derivative_3d Dx(world, 0);
    Dx.read_from_file(DATA_PATH);

    complex_derivative_3d Dy(world, 1);
    Dy.read_from_file(DATA_PATH);

    complex_derivative_3d Dz(world, 2); 
    Dz.read_from_file(DATA_PATH);

    comp_vecfuncT Ueff1u = apply(world, Dx, mul(world, B, dux));
    comp_vecfuncT Ueff2u = apply(world, Dy, mul(world, B, duy));
    comp_vecfuncT Ueff3u = apply(world, Dz, mul(world, B, duz));

    truncate2(world, Ueff1u, Ueff2u);
    truncate(world, Ueff3u, ntol);

    compress(world, Ueff1u);
    compress(world, Ueff2u);
    compress(world, Ueff3u);
    world.gop.fence();

//    gaxpy(world, 1.0, Vu, -1.0, Ueff1u); // gif // need to make it double_complex
//    gaxpy(world, 1.0, Vu, -1.0, Ueff2u);
//    gaxpy(world, 1.0, Vu, -1.0, Ueff3u);

    gaxpy(world, one, Vu, mone, Ueff1u);    // IS (made double complex) 
    gaxpy(world, one, Vu, mone, Ueff2u);
    gaxpy(world, one, Vu, mone, Ueff3u);

    world.gop.fence();
    Ueff1u.clear(); Ueff2u.clear(); Ueff3u.clear();
    world.gop.fence();    

    comp_vecfuncT Ueff1d = apply(world, Dx, mul(world, B, ddx));
    comp_vecfuncT Ueff2d = apply(world, Dy, mul(world, B, ddy));
    comp_vecfuncT Ueff3d = apply(world, Dz, mul(world, B, ddz));

    truncate2(world, Ueff1d, Ueff2d);
    truncate(world, Ueff3d, ntol);

    compress(world, Ueff1d);
    compress(world, Ueff2d);
    compress(world, Ueff3d);
    world.gop.fence();

//    gaxpy(world, 1.0, Vd, -1.0, Ueff1d); // ditto
//    gaxpy(world, 1.0, Vd, -1.0, Ueff2d);
//    gaxpy(world, 1.0, Vd, -1.0, Ueff3d);

    gaxpy(world, one, Vd, mone, Ueff1d);    // IS (made double complex)
    gaxpy(world, one, Vd, mone, Ueff2d);
    gaxpy(world, one, Vd, mone, Ueff3d);

    world.gop.fence();
    Ueff1d.clear(); Ueff2d.clear(); Ueff3d.clear();
    world.gop.fence(); // FENCE
} 



//Calculates spin-orbit potential
void Uso(World& world, const comp_vecfuncT& u, const comp_vecfuncT& v, 
comp_vecfuncT& Vu, comp_vecfuncT& Vv, const real_functionT& W, const double brad,
const comp_vecfuncT& dux, const comp_vecfuncT& ddx, const comp_vecfuncT& duy, 
const comp_vecfuncT& ddy, const comp_vecfuncT& duz, const comp_vecfuncT& ddz)
{
    const double_complex one(1.0,0.0);
    const double lambda = -1.0;
    const double_complex lam(lambda,0.0);

    W.reconstruct();
    world.gop.fence();      // FENCE

    // Bryan's edits change D type here  
    //real_functionT Wx = real_derivative_3d(world,0)(W);
    //real_functionT Wy = real_derivative_3d(world,1)(W);
    //real_functionT Wz = real_derivative_3d(world,2)(W);

    real_derivative_3d derx(world, 0);
    derx.read_from_file(DATA_PATH);

    real_derivative_3d dery(world, 1);
    dery.read_from_file(DATA_PATH);
 
    real_derivative_3d derz(world, 2);
    derz.read_from_file(DATA_PATH);

    real_functionT Wx = derx(W);
    real_functionT Wy = dery(W);
    real_functionT Wz = derz(W);

    complex_derivative_3d Dx(world, 0); 
    complex_derivative_3d Dy(world, 1);
    complex_derivative_3d Dz(world, 2);

    // Bryan's edits
    Dx.read_from_file(DATA_PATH);
    Dy.read_from_file(DATA_PATH);
    Dz.read_from_file(DATA_PATH);

    compress(world, Vu);
    compress(world, Vv);

/*   
    gaxpy(world, one, Vu, -I*lam, mul(world, Wy, dux), false);
    gaxpy(world, one, Vv,  I*lam, mul(world, Wy, ddx), false);
    gaxpy(world, one, Vv,   -lam, mul(world, Wz, dux), false);
    gaxpy(world, one, Vu,    lam, mul(world, Wz, ddx), false);

    gaxpy(world, one, Vu,  I*lam, mul(world, Wx, duy), false);
    gaxpy(world, one, Vv, -I*lam, mul(world, Wx, ddy), false);
    gaxpy(world, one, Vv, -I*lam, mul(world, Wz, duy), false);
    gaxpy(world, one, Vu, -I*lam, mul(world, Wz, ddy), false);   

    gaxpy(world, one, Vu,   -lam, mul(world, Wx, ddz), false);
    gaxpy(world, one, Vv,    lam, mul(world, Wx, duz), false);
    gaxpy(world, one, Vu,  I*lam, mul(world, Wy, ddz), false);
    gaxpy(world, one, Vv,  I*lam, mul(world, Wy, duz), false); 
*/    

    gaxpy(world, one, Vu, -I*lam, mul(world, Wy, dux));
    gaxpy(world, one, Vv,  I*lam, mul(world, Wy, ddx));
    gaxpy(world, one, Vv,   -lam, mul(world, Wz, dux));
    gaxpy(world, one, Vu,    lam, mul(world, Wz, ddx));

    gaxpy(world, one, Vu,  I*lam, mul(world, Wx, duy));
    gaxpy(world, one, Vv, -I*lam, mul(world, Wx, ddy));
    gaxpy(world, one, Vv, -I*lam, mul(world, Wz, duy));
    gaxpy(world, one, Vu, -I*lam, mul(world, Wz, ddy));

    gaxpy(world, one, Vu,   -lam, mul(world, Wx, ddz));
    gaxpy(world, one, Vv,    lam, mul(world, Wx, duz));
    gaxpy(world, one, Vu,  I*lam, mul(world, Wy, ddz));
    gaxpy(world, one, Vv,  I*lam, mul(world, Wy, duz));

    world.gop.fence();              // FENCE
}




//Calculates matrix of kinetic energies
comp_tensorT Kmatrix(World& world, const comp_vecfuncT& u, const comp_vecfuncT& v,
const comp_vecfuncT& dux, const comp_vecfuncT& ddx, const comp_vecfuncT& duy,
const comp_vecfuncT& ddy, const comp_vecfuncT& duz, const comp_vecfuncT& ddz)
{
    double ntol = FunctionDefaults<3>::get_thresh() * prec;     // IS
    double kf;
    int ust = u.size();
    const double_complex one(1.0,0.0);
    comp_vecfuncT temp_u, temp_v;

    if(ust == Z) { kf = kfp; }
    else{ kf = kfn; }

    reconstruct(world, u);
    reconstruct(world, v);
            
    complex_derivative_3d Dx(world, 0);
    complex_derivative_3d Dy(world, 1);
    complex_derivative_3d Dz(world, 2);

    // Bryan's edits
    Dx.read_from_file(DATA_PATH);
    Dy.read_from_file(DATA_PATH);
    Dz.read_from_file(DATA_PATH);

    temp_u = apply(world, Dx, dux);                            
    truncate(world, temp_u, ntol);                              // IS
    gaxpy(world, one, temp_u, double_complex(1,0), apply(world, Dy, duy));
    truncate(world, temp_u, ntol);                              // IS
    gaxpy(world, one, temp_u, double_complex(1,0), apply(world, Dz, duz));
    truncate(world, temp_u, ntol);                              // IS

    temp_v = apply(world, Dx, ddx);
    truncate(world, temp_v, ntol);                              // IS  
    gaxpy(world, one, temp_v, double_complex(1,0), apply(world, Dy, ddy));
    truncate(world, temp_v, ntol);                              // IS  
    gaxpy(world, one, temp_v, double_complex(1,0), apply(world, Dz, ddz));
    truncate(world, temp_v, ntol);                              // IS  

    comp_tensorT r1 = matrix_inner(world, u, temp_u);
    comp_tensorT r2 = matrix_inner(world, v, temp_v);
    comp_tensorT r = r1+r2;

/*  
  comp_tensorT r1 = matrix_inner(world, u, apply(world, Dx, dux), true);
  comp_tensorT r2 = matrix_inner(world, v, apply(world, Dx, ddx), true);    
  comp_tensorT r3 = matrix_inner(world, u, apply(world, Dy, duy), true);
  comp_tensorT r4 = matrix_inner(world, v, apply(world, Dy, ddy), true);
  comp_tensorT r5 = matrix_inner(world, u, apply(world, Dz, duz), true);
  comp_tensorT r6 = matrix_inner(world, v, apply(world, Dz, ddz), true);
  comp_tensorT r = r1 + r2 + r3 + r4 + r5 + r6;
*/

    r = -kf*r;
    return r;
}



//Calculates spin-orbit density
//real_functionT so_density(World& world, const comp_vecfuncT& u, const comp_vecfuncT& v)
real_functionT so_density(World& world, const comp_vecfuncT& dux, const comp_vecfuncT& ddx,
const comp_vecfuncT& duy, const comp_vecfuncT& ddy, const comp_vecfuncT& duz, const comp_vecfuncT& ddz,
const comp_vecfuncT& u, const comp_vecfuncT& v)
{
    double ntol = FunctionDefaults<3>::get_thresh() * prec;     // IS
    comp_functionT djden = comp_factoryT(world);

    const double_complex one(1.0,0.0);

    complex_derivative_3d Dx(world, 0);
    complex_derivative_3d Dy(world, 1);
    complex_derivative_3d Dz(world, 2);

    // Bryan's edits
    Dx.read_from_file(DATA_PATH);
    Dy.read_from_file(DATA_PATH);
    Dz.read_from_file(DATA_PATH);

    comp_vecfuncT dJ = mul(world, conj(world, dux), duy);
    truncate(world, dJ, ntol);
    compress(world, dJ);

//    gaxpy(world, one, dJ,   I , mul(world, conj(world, dux), ddz), false);
//    gaxpy(world, one, dJ, -one, mul(world, conj(world, duy), dux), false);
//    gaxpy(world, one, dJ,  one, mul(world, conj(world, duy), ddz), false);
//    gaxpy(world, one, dJ,  -I , mul(world, conj(world, duz), ddx), false);
//    gaxpy(world, one, dJ, -one, mul(world, conj(world, duz), ddy), false);

    gaxpy(world, one, dJ,   I , mul(world, conj(world, dux), ddz));
    gaxpy(world, one, dJ, -one, mul(world, conj(world, duy), dux));
    gaxpy(world, one, dJ,  one, mul(world, conj(world, duy), ddz));
    gaxpy(world, one, dJ,  -I , mul(world, conj(world, duz), ddx));
    gaxpy(world, one, dJ, -one, mul(world, conj(world, duz), ddy));


//    gaxpy(world, one, dJ,  -I , mul(world, conj(world, ddx), duz), false);
//    gaxpy(world, one, dJ, -one, mul(world, conj(world, ddx), ddy), false);
//    gaxpy(world, one, dJ,  one, mul(world, conj(world, ddy), duz), false);
//    gaxpy(world, one, dJ,  one, mul(world, conj(world, ddy), ddx), false);
//    gaxpy(world, one, dJ,   I , mul(world, conj(world, ddz), dux), false);
//    gaxpy(world, one, dJ, -one, mul(world, conj(world, ddz), duy), false);

    gaxpy(world, one, dJ,  -I , mul(world, conj(world, ddx), duz));
    gaxpy(world, one, dJ, -one, mul(world, conj(world, ddx), ddy));
    gaxpy(world, one, dJ,  one, mul(world, conj(world, ddy), duz));
    gaxpy(world, one, dJ,  one, mul(world, conj(world, ddy), ddx));
    gaxpy(world, one, dJ,   I , mul(world, conj(world, ddz), dux));
    gaxpy(world, one, dJ, -one, mul(world, conj(world, ddz), duy));

//    world.gop.fence();

    truncate(world, dJ, ntol);                                // IS?
    compress(world, dJ);

    djden = dJ[0];      //new
    world.gop.fence();  // FENCE   
    for(unsigned int i=1; i < u.size(); i++)
    { 
        djden.gaxpy(double_complex(1.0,0.), dJ[i], double_complex(1.0,0.), false); 
    }
    world.gop.fence();  // FENCE   

    return real(-1.0*I*djden);
}  



//double abssq(comp_functionT u){
//  double stuff=u.real()*u.real()+u.imag()*u.image();
//  return stuff;
// }


//Calculates kinetic density tau
real_functionT kdensities(World& world, const comp_vecfuncT& dux, const comp_vecfuncT& ddx, 
const comp_vecfuncT& duy, const comp_vecfuncT& ddy, const comp_vecfuncT& duz, const comp_vecfuncT& ddz)
{
    comp_functionT kdensity = comp_factoryT(world);
    comp_vecfuncT kabs = add(world, mul(world, conj(world, dux), dux), mul(world, conj(world, ddx), ddx));

    compress(world, kabs);
    compress(world, duy);
    compress(world, ddy);
    compress(world, ddz);

    gaxpy(world, double_complex(1.0, 0.), kabs, double_complex(1.0, 0.), add(world, mul(world, conj(world, duy), duy), mul(world, conj(world, ddy), ddy)));
    gaxpy(world, double_complex(1.0, 0.), kabs, double_complex(1.0, 0.), add(world, mul(world, conj(world, duz), duz), mul(world, conj(world, ddz), ddz)));

    compress(world, kabs);
    kdensity.compress();
    world.gop.fence();

    for (unsigned int i=0; i<kabs.size(); i++) 
    {
        kdensity.gaxpy(double_complex(1.0, 0.), kabs[i], double_complex(1.0,0.), false);  //gif
    }
    world.gop.fence();  // FENCE

    kdensity.compress();
    real_functionT rkdensity = real(kdensity);                          // IS
    return binary_op(rkdensity, rkdensity, abs_func());                 // IS (makes sure that kinetic density is positive)
}



//Calculates number density rho
real_functionT densities(World& world, const comp_vecfuncT& uu, const comp_vecfuncT& ud )
{
    comp_functionT cdensity;
  
    compress(world, uu);
    compress(world, ud);

    comp_vecfuncT anorm= add(world, mul(world, conj(world, uu), uu), mul(world, conj(world, ud), ud));
    compress(world, anorm);
  
    cdensity = anorm[0];
    cdensity.compress();
    world.gop.fence();  // FENCE
    for(unsigned int i=1; i<uu.size(); i++)
    {
        cdensity.gaxpy( double_complex(1.0,0.), anorm[i], double_complex(1.0,0.), false); //?
    }

    world.gop.fence();  // FENCE    
    cdensity.compress();
    real_functionT density = real(cdensity);                            // IS
    return binary_op(density, density, abs_func());                     // IS (makes sure that number density is positive) 
}




//Orthonormalizes function-vector (nucleon wavefunctions) computing eigenvectors 
void orthonorm( World& world, comp_vecfuncT& uu, comp_vecfuncT& ud, comp_vecfuncT& Vuu,
        comp_vecfuncT& Vud, real_tensorT& E, double& Es, const double tol, const real_functionT& V,
        const comp_vecfuncT& dux, const comp_vecfuncT& ddx, const comp_vecfuncT& duy,
        const comp_vecfuncT& ddy, const comp_vecfuncT& duz, const comp_vecfuncT& ddz )
{
    int ust = uu.size();
    double kf;
    if( ust == Z){ kf = kfp;}
    else{ kf = kfn;}

    comp_tensorT C; 
    comp_tensorT H(ust,ust);
    comp_tensorT S(ust,ust);

    if(world.rank() == 0){print("           Kmatrix ... ");}
    comp_tensorT r = Kmatrix(world, uu, ud, dux, ddx, duy, ddy, duz, ddz);

    if(world.rank() == 0){print("           Hmatrix ... ");}
    S = matrix_inner(world, uu, uu, true) + matrix_inner(world, ud, ud, true);
    H = r + matrix_inner(world, uu, Vuu, true) + matrix_inner(world, ud, Vud, true);
    world.gop.fence();

    if(world.rank() == 0){print("           sygv ... ");}
    sygv(H, S, 1, C, E);  
    world.gop.fence();  // FENCE


    if(world.rank() == 0){print("           transform ... ");}

    if((spinorbit != 1) && (meff != 1))
    {
        E   = E(Slice(0,ust-1));
        world.gop.fence();
        uu  = transform(world, uu,  C(_,Slice(0,ust-1)), 1.e-3*tol, false);
        ud  = transform(world, ud,  C(_,Slice(0,ust-1)), 1.e-3*tol, false);
        world.gop.fence();
        Vuu = mul(world, V, uu);
        Vud = mul(world, V, ud);

//         uu  = transform(world, uu,  C(_,Slice(0,ust-1)), tol, false);
//         ud  = transform(world, ud,  C(_,Slice(0,ust-1)), tol, false);
//         Vuu = mul(world, V, uu, false);
//         Vud = mul(world, V, ud, false);
    }
    else
    {
        E   = E(Slice(0,ust-1));
        world.gop.fence();
        uu  = transform(world, uu,  C(_,Slice(0,ust-1)), 1.e-3*tol, false);
        ud  = transform(world, ud,  C(_,Slice(0,ust-1)), 1.e-3*tol, false);
        Vuu = transform(world, Vuu, C(_,Slice(0,ust-1)), 1.e-3*tol, false);
        Vud = transform(world, Vud, C(_,Slice(0,ust-1)), 1.e-3*tol, false);
        world.gop.fence();

//        uu  = transform(world, uu,  C(_,Slice(0,ust-1)), tol, false);
//        ud  = transform(world, ud,  C(_,Slice(0,ust-1)), tol, false);
//        Vuu = transform(world, Vuu, C(_,Slice(0,ust-1)), tol, false);
//        Vud = transform(world, Vud, C(_,Slice(0,ust-1)), tol, false);

    }
    world.gop.fence();  // FENCE

    double maxen = 0.0;
    double minen = 0.0;
    for(unsigned int i=0; i<uu.size(); i++) 
    { 
        maxen = std::max(maxen,E[i]);
        minen = std::min(minen,E[i]);
    }

//    world.gop.fence();
//    Es = std::max(0.5, 2.0*maxen);
    Es = std::max(15.0, 2.0*maxen);
    if( maxen > 0 ){if(world.rank() == 0){print("           Es = ", maxen);}}

    truncate2(world, uu, ud);
    truncate2(world, Vuu, Vud);
    normalize_ud(world, uu, ud);

    if(world.rank() == 0){print("           gaxpy ... ");}
    compress(world, Vuu);
    compress(world, Vud);
    compress(world, uu);
    compress(world, ud);

//    gaxpy(world, 1.0, Vuu, -1.0*Es, uu, false);
//    gaxpy(world, 1.0, Vud, -1.0*Es, ud, false);

    gaxpy(world, 1.0, Vuu, -1.0*Es, uu);
    gaxpy(world, 1.0, Vud, -1.0*Es, ud);

//    world.gop.fence();

    scale(world, Vuu, -1.0/kf);
    scale(world, Vud, -1.0/kf);
}




//Makes BSH Operators
std::vector<poperatorT> BSHoperators(World& world, const real_tensorT& evals, const double Es, 
const double tol)
{
    double ntol = FunctionDefaults<3>::get_thresh() * 0.01 ; //gf
    int n = evals.dim(0);
    double kf;
    if( n == Z){ kf = kfp; }
    else{ kf = kfn; }
    std::vector<poperatorT> ops(n);
    double lo = ntol * 0.1; // gc

    for (int i=0; i<n; i++)
    {
        double eps = evals(i);
        eps -= Es;
        if (eps > 0.0) { if(world.rank() == 0) { print( eps, Es );} eps = -0.05;}
        ops[i] = poperatorT(BSHOperatorPtr3D(world, sqrt(-(1.0/kf)*eps), lo, ntol)); //gf
    }   
    world.gop.fence();
    return ops;
}




//Iteration routine
void iterate(World& world,
        const real_functionT& U,
        comp_vecfuncT& uu,
        comp_vecfuncT& ud, 
        real_tensorT& eps,
        double& maxerr,
        const real_functionT& rho,
        real_functionT& rho_u, 
        const int iter,
        const double avg,
        const double tol,
        const double brad)
{
    double ntol = FunctionDefaults<3>::get_thresh() * prec;
    int nt = uu.size();
    double Es = 0.0;

    complex_derivative_3d Dx(world, 0);
    complex_derivative_3d Dy(world, 1);
    complex_derivative_3d Dz(world, 2);
    
    // Bryan's edits
    Dx.read_from_file(DATA_PATH);
    Dy.read_from_file(DATA_PATH);
    Dz.read_from_file(DATA_PATH);
    
    comp_vecfuncT dux = apply(world, Dx, uu);
    truncate(world, dux, ntol);                       
    comp_vecfuncT duy = apply(world, Dy, uu);
    truncate(world, duy, ntol);                       
    comp_vecfuncT duz = apply(world, Dz, uu);
    truncate(world, duz, ntol);                       

    loadbalance_v3(world, dux, duy, duz);

    comp_vecfuncT ddx = apply(world, Dx, ud);
    truncate(world, ddx, ntol);                       
    comp_vecfuncT ddy = apply(world, Dy, ud);
    truncate(world, ddy, ntol);                       
    comp_vecfuncT ddz = apply(world, Dz, ud);
    truncate(world, ddz, ntol);                       

    loadbalance_v3(world, ddx, ddy, ddz);
    
    if(world.rank() == 0){print("       orthonormalize ... ");}

    comp_vecfuncT Vuu = mul(world, U, uu);
    comp_vecfuncT Vud = mul(world, U, ud);


    if(spinorbit == 1)
    {	
        if(world.rank() == 0){print("           Uso ... ");}
        real_functionT W = b4*rho + bp4*rho_u;
        Uso(world, uu, ud, Vuu, Vud, W, brad, dux, ddx, duy, ddy, duz, ddz);
        truncate2(world, Vuu, Vud);
    }

    if(meff == 1)
    {
        if(world.rank() == 0){print("           Umeff ... ");}
        real_functionT B = b1*rho - bp1*rho_u;
        Umeff(world, Vuu, Vud, B, dux, duy, duz, ddx, ddy, ddz);
        truncate2(world, Vuu, Vud);
    }
    
    if(world.rank() == 0){print("           orthonorm ... ");}
    orthonorm(world, uu, ud, Vuu, Vud, eps, Es, tol, U, dux, ddx, duy, ddy, duz, ddz);
    loadbalance_v2(world, Vuu, Vud);
    truncate2(world, Vuu, Vud);

    world.gop.fence();  // FENCE
    dux.clear(); ddx.clear(); duy.clear(); ddy.clear(); duz.clear(); ddz.clear();
    world.gop.fence();

    if(world.rank() == 0){print("       apply BSH ... ");}
    std::vector<poperatorT> ops = BSHoperators(world, eps, Es, tol);

    comp_vecfuncT tmp_u = apply(world, ops, Vuu);
    comp_vecfuncT tmp_d = apply(world, ops, Vud);
    truncate2(world, tmp_u, tmp_d);

    world.gop.fence();  // FENCE // gif
    Vuu.clear(); ops.clear(); Vud.clear();
    world.gop.fence();  // FENCE // gif

    loadbalance_v2(world, tmp_u, tmp_d);

    if(world.rank() == 0){print("       update ... ");}
    normalize_ud(world, tmp_u, tmp_d);

    std::vector<double> rnormu = norm2s(world,sub(world, uu, tmp_u));
    std::vector<double> rnormd = norm2s(world,sub(world, ud, tmp_d));
    std::vector<double> rnorm(nt);

    compress(world, tmp_u);
    compress(world, tmp_d);
    compress(world, uu);
    compress(world, ud);

    if( avg_wav == 1)
    {
        gaxpy(world, 1.0-avg, uu, avg, tmp_u);
        gaxpy(world, 1.0-avg, ud, avg, tmp_d);
    }
    else
    {
        uu = copy(world, tmp_u);
        ud = copy(world, tmp_d);
    }

    world.gop.fence();  // FENCE // gif
    tmp_u.clear();
    tmp_d.clear();
    world.gop.fence();  // FENCE // gif

    for (int i=0; i<nt; i++)
    {
        rnorm[i] = std::sqrt(rnormu[i]*rnormu[i] + rnormd[i]*rnormd[i]);
    }
    world.gop.fence(); // FENCE

    int maxind = 0;
    for(unsigned int i = 0; i < uu.size(); i++)
    {
        maxerr = std::max(maxerr,rnorm[i]);
        if( maxerr == rnorm[i]){ maxind = i;}
    }
    world.gop.fence();  // FENCE

    normalize_ud(world, uu, ud);                    // IS
    rho_u = densities(world, uu, ud);               // IS

    truncate2(world, uu, ud);
    normalize_ud(world, uu, ud);

    world.gop.fence();  // FENCE
}




//Make a table of kinetic and potential energies for each state
void Energies(World& world, const real_functionT& potential,
        const comp_vecfuncT& uu,    const comp_vecfuncT& ud, const real_functionT& rho, 
        const real_functionT& rho_q, real_tensorT& Etot, const int iter, const comp_tensorT& r)
{
    reconstruct(world, uu);
    reconstruct(world, ud);
    int nv = uu.size();

    real_tensorT Ekin(nv);
    for( int i=0; i < nv; i++) { Ekin[i] += r(i,i).real(); }
    world.gop.fence();  // FENCE

    //(Only) In the 0th iteration - don't include spin-orbit potential 
    //and effective mass (quicker to calculate)
    if(iter == 0)
    {
        real_tensorT Epot(nv);
        world.gop.fence();  // FENCE
        for(int i=0; i < nv; i++)
        {
            double_complex epot = inner(uu[i], potential*uu[i])  
                                + inner(ud[i], potential*ud[i]); // gif ; need to vectorize
            Epot[i] = epot.real();
        }
        Etot = Ekin + Epot; 
    }

    if(world.rank() == 0)
    {
        std::ostringstream filename1;
        filename1 << "En_0" << std::setw(6) << std::setfill('0') << iter << ".txt";
        FILE * fid1 = fopen(filename1.str().c_str(), "a");
        fprintf(fid1,"\n");
        for(int i=0; i < nv; i++)
        {
            fprintf(fid1, "%u,\t%16.8f,\t%16.8f,\t%16.8f\n", i, Ekin[i], Etot[i] - Ekin[i], Etot[i]);
        }
        if(world.rank() == 0) { fprintf(fid1,"\n");}
        if(world.rank() == 0) { fclose(fid1); }
    }
    world.gop.fence();  // FENCE
}



// Calculate total binding energy
double Ebinding(World& world, const real_functionT& rho_p, const real_functionT& rho_n,
    const real_functionT& tau_p,  const real_functionT& tau_n, const real_functionT& lap_p, 
    const real_functionT& lap_n,  const real_functionT& dj_p,  const real_functionT& dj_n, 
    const real_functionT& U_c,    const real_functionT& U_e,   const comp_vecfuncT& uu,
    const comp_vecfuncT& ud,      const int iter)
{
    int ust = uu.size();
    double kf;
    if( ust == Z){ kf = kfp;}
    else{ kf = kfn;}

    real_functionT rho = rho_p + rho_n;
    real_functionT U01 =    0.5e0*b0*rho*rho;
    real_functionT U02 =  - 0.5e0*bp0*(rho_p*rho_p + rho_n*rho_n);
    double E0 = U01.trace() + U02.trace();

    real_functionT U11 =   b1*(rho_p + rho_n)*(tau_p + tau_n);
    real_functionT U12 = -1.e0*(bp1*(rho_p*tau_p) + bp1*(rho_n*tau_n));
    double E1 = U11.trace() + U12.trace();

    real_functionT U21 = - 0.5e0*b2*(rho_p + rho_n)*(lap_p + lap_n);
    real_functionT U22 = + 0.5e0*bp2*rho_p*lap_p + 0.5e0*bp2*rho_n*lap_n;
    double E2 = U21.trace() + U22.trace();

    real_functionT ho2 = binary_op(rho, rho, rho4());
    real_functionT U31 =  (1.e0/3.e0)*b3*ho2*rho*rho;
    real_functionT U32 = -(1.e0/3.e0)*bp3*ho2*(rho_p*rho_p + rho_n*rho_n);
    double E3 = U31.trace() + U32.trace();

    double E4 = 0.0;
    if(spinorbit == 1)
    {
        real_functionT U41 = -b4*rho*(dj_n + dj_p);
        real_functionT U42 = -bp4*(rho_p*dj_p + rho_n*dj_n);
        E4 = U41.trace() + U42.trace();
    }

    double Ec = 0.0;
    world.gop.fence();
    for(unsigned int i=0; i<uu.size(); i++)
    {
        double_complex ec = 0.5e0*inner(U_c*uu[i],uu[i]) + 0.5e0*inner(U_c*ud[i],ud[i]);
        Ec += ec.real();
    }
    world.gop.fence();


    double Ee = 0.0;
  /*
    real_functionT rho_e = real_factoryT(world);
    rho_e = 0.0*rho_e + (1.0*Z/(8.0*L*L*L));
    double_complex ee = 0.5*inner(U_c,rho_e);
    Ee = ee.real();
    Ee = 0.0;
  */

  //real_functionT U_ex = Uex(world, rho_p);
    real_functionT U_ex = binary_op(rho_p, rho_p, uex());
    U_ex = 0.75e0*e2*rho_p*rho_p*U_ex;
    double Eex = U_ex.trace();

    double Ekin  = kf*(tau_p.trace()) + kf*(tau_n.trace());
    double Ebind = E0 + E1 + E2 + E3 + E4 + Ec + Ee + Eex + Ekin;

    if( iter%10 == 0 || iter == 0)
    {
        if(world.rank() == 0)
        {
            std::ostringstream filename1;
            filename1 << "En_0" << std::setw(6) << std::setfill('0') << iter << ".txt";
            FILE * fid1 = fopen(filename1.str().c_str(), "w");
            fprintf(fid1,"\n");
            fprintf(fid1,"%s \t%u\n","Iteration:",iter);
            fprintf(fid1,"%s\n","Energies integrated from density functional: ");
            fprintf(fid1,"\n");
            fprintf(fid1,"%s \t%16.8f,\t%s \t%16.8f,\t%s \t%18.8f\n","E0 part:",E0,"E1 part:",E1,"E2 part:",E2);
            fprintf(fid1,"%s \t%16.8f,\t%s \t%16.8f\n","E3 part:",E3,"E4 part:",E4);
            fprintf(fid1,"%s \t%16.8f,\t%s \t%16.8f,\t%s \t%16.8f\n","Coulomb:",Ec + Eex,"Direct:",Ec,"Exchange:",Eex);
            fprintf(fid1,"%s \t%16.8f\n","Kinetic:",Ekin);
            fprintf(fid1,"%s \t%16.8f\n","Total  :",E0 + E1 + E2 + E3 + E4 + Ec + Ee + Eex + Ekin);
            fclose(fid1);
        }
    }
    world.gop.fence();
    return Ebind/A;
}



//Make output for plotting
void output(World& world,  const real_functionT& rho_p, const real_functionT& rho_n,
	    const real_functionT& tau, const real_functionT& lap_p, const real_functionT& lap_n,
	    const real_functionT& U,   const real_functionT& Uc,    const real_functionT& Uex,
	    const comp_vecfuncT& u,    const real_functionT& rho_1, const real_functionT& rho_2,
	    const real_functionT& rho_3, const real_functionT& djp, const real_functionT& djn, 
	    const int iter)
{
    if(world.rank() == 0){print(" "); print("Make output: ", iter);}
    double Lp = std::min(L , 24.0);

    real_functionT lap = lap_n + lap_p;
    real_functionT rho = rho_p + rho_n;

    //real_functionT test  = real_factory_3d(world).f(myf);
    //real_functionT dtestx = real_derivative_3d(world,0)(test);  
    //dtestx.truncate(1.0e-3);
    //dtestx.reconstruct();
    //real_functionT ddtestx = real_derivative_3d(world,0)(dtestx);


    // Paraview Output
    if (VTK_OUTPUT)
    {
        Vector<double, 3> plotlo, plothi;
        Vector<long, 3> npts;
        char filename[100];    
        sprintf(filename, "%s/paraview_%d.vts",direct,iter);
        world.gop.fence();  // FENCE

        for(int i = 0; i < 3; i++){plotlo[i] = -Lp; plothi[i] = Lp; npts[i] = 51;}
        world.gop.fence();  // FENCE    
        plotvtk_begin(world, filename, plotlo, plothi, npts);
        plotvtk_data(rho,   "rho", world, filename, plotlo, plothi, npts);
        plotvtk_data(rho_p, "rho_p", world, filename, plotlo, plothi, npts);
        plotvtk_data(rho_n, "rho_n", world, filename, plotlo, plothi, npts);
        plotvtk_data(tau,   "tau", world, filename, plotlo, plothi, npts);
        plotvtk_data(lap,   "lap", world, filename, plotlo, plothi, npts);
        plotvtk_data(U,     "U", world, filename, plotlo, plothi, npts);
        plotvtk_data(rho_1, "rho1", world, filename, plotlo, plothi, npts);
        plotvtk_data(rho_2, "rho2", world, filename, plotlo, plothi, npts);
        plotvtk_data(rho_3, "rho3", world, filename, plotlo, plothi, npts);
        plotvtk_end<3>(world, filename);
    }

    if ( PLOT_DENSITIES )
    {
        coord_3d lo,hi;
        // Ouput for profile along y
        lo[0] = 0.0; lo[1] = -Lp; lo[2] = 0.0;
        hi[0] = 0.0; hi[1] =  Lp; hi[2] = 0.0;
        char plotname[500];
        sprintf(plotname, "densities_y%d.txt", iter);
        plot_line(plotname, 501, lo, hi, rho, tau, lap);

        // Output for profile along x
        lo[0] = -Lp; lo[1] = 0.0; lo[2] = 0.0;
        hi[0] =  Lp; hi[1] = 0.0; hi[2] = 0.0;
        sprintf(plotname, "densities_x%d.txt", iter);
        plot_line(plotname, 5001, lo, hi, rho, tau, lap);

        // Output for profiles along z
        lo[0] = 0.0; lo[1] = 0.0; lo[2] = -Lp;
        hi[0] = 0.0; hi[1] = 0.0; hi[2] =  Lp;
        sprintf(plotname, "densities_z%d.txt", iter);
        plot_line(plotname, 501, lo, hi, rho, tau, lap);

        // Output of potential and states
        lo[0] = 0.0; lo[1] = -Lp; lo[2] = 0.0;
        hi[0] = 0.0; hi[1] =  Lp; hi[2] = 0.0;
        sprintf(plotname, "potential_%d.txt", iter);
        plot_line(plotname, 501, lo, hi, U, Uc, Uex);
    }

    world.gop.fence(); // FENCE
}





//Makes potential
void Potential(World& world,
           const comp_vecfuncT& pu,
           const comp_vecfuncT& pd, 
           const comp_vecfuncT& nu,
           const comp_vecfuncT& nd,
           real_tensorT& energy_p,
           real_tensorT& energy_n,
           real_functionT& potential_p,
           real_functionT& potential_n, 
           real_functionT& rho_p,
           real_functionT& rho_n,
           real_functionT& rho, 
           const int iter,
           double& BE,
           const double maxerr,
           const double tol,
           const double brad, 
           real_functionT& potential_p_old,
           real_functionT& potential_n_old, 
           real_functionT& lap_p_old,
           real_functionT& lap_n_old,
           const double avg)
{
    double ntol = FunctionDefaults<3>::get_thresh() * prec;  
    int np = pu.size();
    int nn = nu.size();


    // Make number density and kinetic densities
    if(world.rank() == 0){print("   Make densities");}
    if( BE == 0.0 )                                     // IS
    {                                                   // IS
//    START_TIMER;                                      // IS
        rho_p = densities(world, pu, pd);               // IS
        rho_n = densities(world, nu, nd);               // IS
//    END_TIMER("Densities");                           // IS
    }                                                   // IS

    double p_number = rho_p.trace();
    if(world.rank() == 0){print("   A =", p_number);}
    double n_number = rho_n.trace();
    if(world.rank() == 0){print("   A =", n_number);}
    world.gop.fence();  // FENCE


    if(world.rank() == 0){print("   Make derivatives 1");}
    complex_derivative_3d Dx(world, 0);
    
    // Bryan's edits
    // Change derivative type
    Dx.read_from_file(DATA_PATH);
 
    comp_vecfuncT dpux = apply(world, Dx, pu); 
    comp_vecfuncT dpdx = apply(world, Dx, pd);
//    truncate(world, dpux, ntol);                      // IS
//    truncate(world, dpdx, ntol);                      // IS

    if(world.rank() == 0){print("   Make derivatives 2");}
    complex_derivative_3d Dy(world, 1);

    // Bryan's edits
    // Change derivative type
    Dy.read_from_file(DATA_PATH);

    comp_vecfuncT dpuy = apply(world, Dy, pu);
    comp_vecfuncT dpdy = apply(world, Dy, pd);
//    truncate(world, dpuy, ntol);                      // IS
//    truncate(world, dpdy, ntol);                      // IS


    if(world.rank() == 0){print("   Make derivatives 3");}
    complex_derivative_3d Dz(world, 2);

    // Bryan's edits
    // Change derivative type
    Dz.read_from_file(DATA_PATH);

    comp_vecfuncT dpuz = apply(world, Dz, pu);
    comp_vecfuncT dpdz = apply(world, Dz, pd);
//    truncate(world, dpuz, ntol);                      // IS
//    truncate(world, dpdz, ntol);                      // IS


    if(world.rank() == 0){print("   Loadbalancing ");}
    //START_TIMER;  // IS
    loadbalance_v6(world, dpux, dpuy, dpuz, dpdx, dpdy, dpdz);
    //END_TIMER("Loadbalance6");    // IS

    if(world.rank() == 0){print("   Kdensities ");}
    //START_TIMER;  // IS
    real_functionT tau_p = kdensities(world, dpux, dpdx, dpuy, dpdy, dpuz, dpdz);
    //END_TIMER("Kdensities");  // IS
    if(world.rank() == 0){print("   Laplacian2 ");}
    //START_TIMER;  // IS
    real_functionT lap_p = laplacian1(world, dpux, dpdx, dpuy, dpdy, dpuz, dpdz, pu, pd);               // IS
    //real_functionT lap_p = laplacian(world, rho_p, brad);
    //real_functionT lap_p = laplacian2(world, dpux, dpdx, dpuy, dpdy, dpuz, dpdz, pu, pd);
    //END_TIMER("Laplacian2");  // IS

    real_functionT dJ_p  = 0.0*rho_p;

    if(spinorbit != 0)
    {    
        //START_TIMER;  // IS
        dJ_p  = so_density(world, dpux, dpdx, dpuy, dpdy, dpuz, dpdz, pu, pd);
        //END_TIMER("SO-p");    // IS
    }
    
    comp_tensorT rp(np,np);
    if( iter%100 == 0 || iter == 0)
    {rp = Kmatrix(world, pu, pd, dpux, dpdx, dpuy, dpdy, dpuz, dpdz);}

    world.gop.fence();
    dpux.clear(); dpdx.clear(); dpuy.clear(); dpdy.clear();
    dpuz.clear(); dpdz.clear();
    world.gop.fence();

    comp_vecfuncT dnux = apply(world, Dx, nu); 
//    truncate(world, dnux, ntol);                      // IS
    comp_vecfuncT dndx = apply(world, Dx, nd);
//    truncate(world, dndx, ntol);                      // IS

    comp_vecfuncT dnuy = apply(world, Dy, nu);
//    truncate(world, dnuy, ntol);                      // IS
    comp_vecfuncT dndy = apply(world, Dy, nd);
//    truncate(world, dndy, ntol);                      // IS

    comp_vecfuncT dnuz = apply(world, Dz, nu);
//    truncate(world, dnuz, ntol);                      // IS
    comp_vecfuncT dndz = apply(world, Dz, nd);
//    truncate(world, dndz, ntol);                      // IS

    //START_TIMER;
    loadbalance_v6(world, dnux, dnuy, dnuz, dndx, dndy, dndz);
    //END_TIMER("Loadbalance_v6");
    //START_TIMER;
    real_functionT tau_n = kdensities(world, dnux, dndx, dnuy, dndy, dnuz, dndz);
    real_functionT lap_n = laplacian1(world, dnux, dndx, dnuy, dndy, dnuz, dndz, nu, nd);               // IS
    //real_functionT lap_n = laplacian(world, rho_n, brad);
    //real_functionT lap_n = laplacian2(world, dnux, dndx, dnuy, dndy, dnuz, dndz, nu, nd);
    //END_TIMER("Tau_n and Lap_n");

    real_functionT dJ_n = real_factoryT(world);

    if(spinorbit != 0)
    {
        dJ_n  = so_density(world, dnux, dndx, dnuy, dndy, dnuz, dndz, nu, nd);
    }

    comp_tensorT rn(nn,nn);
    if( iter%100 == 0 || iter == 0)
    {rn = Kmatrix(world, nu, nd, dnux, dndx, dnuy, dndy, dnuz, dndz);}

    world.gop.fence();
    dnux.clear(); dndx.clear(); dnuy.clear(); dndy.clear();
    dnuz.clear(); dndz.clear();
    world.gop.fence();

    // Average laplacians if there is no Gaussian smoothing
    if( avg_lap == 1 && brad < 0.0 )
    {
        if(world.rank() == 0){print("   Avg. laplacian");}
        double avg_lap = 0.4;
        lap_p = avg_lap * lap_p + (1.0 - avg_lap) * lap_p_old;
        lap_n = avg_lap * lap_n + (1.0 - avg_lap) * lap_n_old;
    }

    // Make total densities
    rho = rho_p + rho_n;
    real_functionT tau = tau_p + tau_n;
    real_functionT lap = lap_p + lap_n;

    lap_p_old = copy(lap_p);
    lap_n_old = copy(lap_n);


    double nuc_number = rho.trace();
    if(world.rank() == 0){print("   A =", nuc_number);}
    world.gop.fence();  // FENCE


    // Gaussian blurring 
    const coordT origin(0.0);
    double sigma = brad;
    Tensor<double> exponent(1), coeff(1);
    exponent = 1.0/(2.0*sigma*sigma);
    coeff = pow(1.0/(sigma*sqrt(2.0*M_PI)),3.0);
    operatorT G(world, coeff, exponent);

    if( iter >= 0 && brad > 0.0)
    {
        if(world.rank() == 0) print("   Gaussian blurring");
        lap_p = G(lap_p); lap_n = G(lap_n); lap = G(lap);
        tau_p = G(tau_p); tau_n = G(tau_n); tau = G(tau);
        lap_p.truncate(ntol); lap_n.truncate(ntol); lap.truncate(ntol);
        tau_p.truncate(ntol); tau_n.truncate(ntol); tau.truncate(ntol);
        lap_p.reconstruct();  lap_n.reconstruct();  lap.reconstruct();
        tau_p.reconstruct();  tau_n.reconstruct();  tau.reconstruct();
    }
    world.gop.fence(); // FENCE

    real_functionT rho_1   = binary_op(rho, rho_p, rho1());
    real_functionT rho_2_n = binary_op(rho, rho_n, rho2());
    real_functionT rho_2_p = binary_op(rho, rho_p, rho2());
    real_functionT rho_3_n = binary_op(rho, rho_n, rho3());
    real_functionT rho_3_p = binary_op(rho, rho_p, rho3());
    world.gop.fence();  // FENCE


    //  real_convolution_3d cop = CoulombOperator(world, tol*0.1, tol);
    real_convolution_3d cop = CoulombOperator(world, ntol*0.001, ntol*0.01); //gc
    real_functionT U_c = e2*cop(rho_p);

    // Make Exchange potential, no screening
    real_functionT U_ex = binary_op(rho_p, rho_p, uex());
    U_ex = e2*rho_p*U_ex;

    real_functionT rho_e = real_factoryT(world);
    rho_e = 0.0*rho_e + (1.0*Z/(8.0*L*L*L));
    double mu = 1.0/screenl;
    real_convolution_3d bop = BSHOperator3D(world, mu, tol, ntol);

    // Make Exchange potential, yes screening
    if( screening == 1)
    {
        real_functionT l  = 0.0*rho_p + mu;
        real_functionT Fg = binary_op(rho_p, l, Fgamma());
        U_ex = e2*rho_p*Fg*binary_op(rho_p, rho_p, uex());
    }
    // Make Coulomb potential, jellium, no screening
    if( jellium == 1 && screening == 0)
    { 
        if(world.rank() == 0){print("   Jellium, no screening");}
        U_c = e2*cop(rho_p - rho_e);
    }
    // Make Coulomb potential, jellium, screening
    if( jellium == 1 && screening == 1)
    { 
        if(world.rank() == 0){print("   Jellium, yes screening");}
        U_c = e2*bop(4.0*M_PI*(rho_p - rho_e));
    }
    // Make Coulomb potential, no jellium, screening
    if( jellium == 0 && screening == 1)
    { 
        if(world.rank() == 0){print("   No jellium, yes screening");}
        U_c = e2*bop(4.0*M_PI*rho_p);
    }
    world.gop.fence(); // FENCE
 

    if(world.rank() == 0){print("   Make U_skyrme");}
    potential_n  = b0*rho - bp0*rho_n
                 + b3*((alpha+2.0)/3.0)*rho_1 - bp3*(2.0/3.0)*rho_2_n
                 - bp3*(alpha/3.0)*( rho_3_p + rho_3_n )
                 + b1*tau_p + (b1-bp1)*tau_n - b2*lap + bp2*lap_n
                 - b4*(dJ_p + dJ_n) - bp4*dJ_n;

    potential_p  = b0*rho - bp0*rho_p
                 + b3*((alpha+2.0)/3.0)*rho_1 - bp3*(2.0/3.0)*rho_2_p
                 - bp3*(alpha/3.0)*( rho_3_p + rho_3_n )
                 + (b1-bp1)*tau_p + b1*tau_n - b2*lap + bp2*lap_p
                 - b4*(dJ_p + dJ_n) - bp4*dJ_p
                 + U_c + U_ex;
    world.gop.fence(); // FENCE


    if( iter > 0 && avg_pot == 1)
    {
        if(world.rank() == 0){print("   Avg. potential");}
        potential_p = avg*potential_p + (1.0-avg)*potential_p_old;
        potential_n = avg*potential_n + (1.0-avg)*potential_n_old;
    }


    // Calculate all energy contributions (just for info - not used in iterations)
    if( iter%10 == 0 || BE == 0.0 )
    {
        if(world.rank() == 0){print("   Make Energies");}
        real_tensorT En(A-Z), Ep(Z);
        BE = Ebinding(world, rho_p, rho_n, tau_p, tau_n, lap_p, lap_n, dJ_p, dJ_n,
             U_c, U_ex, pu, pd, iter);
    }

    // Calculate single particle state energies (for info)
    if( iter%100 == 0 || iter == 0)
    {
        Energies(world, potential_n, nu, nd, rho, rho_n, energy_n, iter, rn);
        Energies(world, potential_p, pu, pd, rho, rho_p, energy_p, iter, rp);
    }

    // Make output
//    if( iter%10 == 0 || iter == 0 )
    {
        output(world, rho_p, rho_n, tau, lap_p, lap_n, potential_p, U_c, U_ex, 
            pu, rho_1, rho_2_p, rho_3_p, dJ_p, dJ_n, iter);
    }

    world.gop.fence();  // FENCE
    tau_p.clear(); tau_n.clear(); tau.clear(); lap_p.clear(); lap_n.clear(); lap.clear();
    rho_1.clear(); rho_2_p.clear(); rho_3_p.clear(); rho_2_n.clear(); rho_3_n.clear();
    dJ_p.clear(); dJ_n.clear(); U_c.clear(); U_ex.clear(); 
    world.gop.fence();  // FENCE
}







// Main routine
void doit(World& world)
{
    long k              = 8;            // wavelet order
    double thresh       = 1.e-4;        // precision
    double avg          = 0.6;          // averaging fraction between old and new wavefunctions         // IS
    double tol          = 1.e-5;        // tolerance (used in e.g. BSH and Coulomb)
    double brad         = 0.25;         // radius of the blurring gaussian
    double time_old     = 0.0;          // checkpointed wall_time
    const int IO_nodes  =   8;          // Number of nodes for checkpointing                            // IS


    // IS
    // Check, if there is a checkpointed file with updated values for above parameters
    FILE *tfile;
    char trecord_name1[100];
    sprintf(trecord_name1, "%s/t_record.00000",direct);
    tfile=fopen(trecord_name1,"r");
    if (tfile!=NULL)
    {
        char trecord_name2[100];
        sprintf(trecord_name2, "%s/t_record",direct);
        archive::ParallelInputArchive tin(world, trecord_name2, IO_nodes);
        tin & thresh;
        tin & k;
        tin & brad;
        tin & avg;
        tin & tol;
        tin & time_old;
        tin.close();
        fclose(tfile);
    }
    // IS



    //  FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap< Key<3> >(world)));

    std::cout.precision(8);
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_autorefine(false); // gif
    FunctionDefaults<3>::set_initial_level(5);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_cubic_cell(-L, L);

    if( periodic == 1)
    { 
        FunctionDefaults<3>::set_bc(BC_PERIODIC); 
        FunctionDefaults<3>::set_truncate_mode(1);
    }        

    FunctionDefaults<3>::set_apply_randomize(false);
    FunctionDefaults<3>::set_project_randomize(false);
    FunctionDefaults<3>::set_truncate_on_project(true);

    std::ofstream outfile;
       
    comp_vecfuncT n(A-Z), p(Z);        // Neutron and proton states
    comp_vecfuncT pu =  zero_functions<double_complex,3>(world,Z);
    comp_vecfuncT pd =  zero_functions<double_complex,3>(world,Z);
    comp_vecfuncT nu =  zero_functions<double_complex,3>(world,A-Z);
    comp_vecfuncT nd =  zero_functions<double_complex,3>(world,A-Z);

    real_tensorT energy_n(A-Z);         // Energy of states
    real_tensorT energy_p(Z);           // Energy of states

    real_tensorT nx(A-Z), ny(A-Z), nz(A-Z); 
    real_tensorT px(Z),   py(Z),   pz(Z);

    real_functionT pot_n, pot_p;
    real_functionT pot_n_old, pot_p_old;
    real_functionT lap_n_old, lap_p_old;
    real_functionT den_p, den_n, den;


    double maxerr       = -10.0;               // Maximum change is single particle states
    double maxerrold    =  20.0;
    double BE           =   0.0;               // Initialization of binding energy
    double BE_old       =   0.0;
    bool file_check     = false;               // Marker for checkpointed files




    // Initialization from functions or text file (only done when there is no checkpointed file)
    //--------------------

    // Check if a checkpoint file exists
    // only use one file to check for existance

    file_check=false;

    if (world.rank() == 0)
    {
        char grecord_name1[100];
        //gif ...this will beat up the servers ... 
        // need to have one node read and then broadcast the results
        sprintf(grecord_name1, "%s/g_record.00000",direct);
        FILE *file;
        file=fopen(grecord_name1,"r");
        if (file==NULL)
        {
            file_check = false;
        }
        else
        {
            file_check = true;
            fclose(file);
        }
    }  
    world.gop.broadcast(file_check);

    if (world.rank() == 0){ print(file_check); print(initial);} 

    // If there is no checkpoint file, initialize wavefunctions. Otherwise read from checkpoint file
    if( file_check == false)                // IS
    {
        if (world.rank() == 0){ print(""); print(""); print("Initialize single particle states");}

        // Read in from file or make initial wavefunctions via subroutine
        if ( initial == 1 )
        {
            FILE *partn = fopen("MD/neutrons.tab", "r");
            if(world.rank() == 0) {std::cout << "- Begin: Read neutron data" << std::endl;}
            for( size_t i = 0; i < (A-Z) ; i++)    
            {
                double x, y, z;
                fscanf(partn, "%lf \t%lf \t%lf\n", &x, &y, &z);
                x -= L;
                y -= L;
                z -= L;
                n[i] = comp_factoryT(world).functor(comp_functorT(new MD(x, y, z)));
                n[i].truncate(prec*thresh);
            }
            fclose(partn);
            world.gop.fence(); // FENCE

            FILE *partp = fopen("MD/protons.tab", "r");
            if(world.rank() == 0) {std::cout << "- Begin: Read proton data" << std::endl;}
            for( size_t i = 0; i < Z ; i++)
            {
                double x, y, z;
                fscanf(partp, "%lf \t%lf \t%lf\n", &x, &y, &z);
                x -= L;
                y -= L;
                z -= L;
                p[i] = comp_factoryT(world).functor(comp_functorT(new MD(x, y, z)));
                p[i].truncate(prec*thresh);
            }
            fclose(partp);
            world.gop.fence();  // FENCE
            if(world.rank() == 0) {std::cout << "- End: Read data" << std::endl;}
        } 
        else if (initial == 2) 
        {
            if (world.rank() == 0){ print("Harmonic Oscillator");}          // IS
            make_ho(world, p);                                              // IS
            make_ho(world, n);                                              // IS
            if (world.rank() == 0){ print("Done");}                         // IS

//	        START_TIMER;
//	        make_ho2(world, p);
//	        END_TIMER("HO for Proton");
//	        START_TIMER;
//	        loadbalance_v1(world, p);
//	        END_TIMER("Loadbalance p");
//	        START_TIMER; // sync in start_timer
//	        make_ho2(world, n);
//	        END_TIMER("HO for Neutron and Proton");
//	        world.gop.fence();  // FENCE
//	        START_TIMER;
//	        loadbalance_v2(world, p, n);
//	        END_TIMER("Loadbalance p and n");
	        if (world.rank() == 0) { print("Done"); }

        }    
        else if (initial == 3)
        {    
            if (world.rank() == 0){ print("Plane Waves");}
            make_fermi(world, p);
            truncate(world, p);
            make_fermi(world, n);
            truncate(world, n);
            normalize_2v(world, n, p);
            if (world.rank() == 0){ print("Done");}
        }
        else
        {
            if (world.rank() == 0){ print("No proper initalization!");}
        }

        truncate2(world, n, p);
        normalize_2v(world, n, p);    
      
        spin_split(world, p, pu, pd); // Split into spin-up and down
        spin_split(world, n, nu, nd); // Split into spin-up and down

        world.gop.fence();  // Irina ...don't get too eager to clear memory
        p.clear(); n.clear();
   }

  // ------------------------------
  // End of Initialization




    // Main iterations
    if(world.rank() == 0)
    {
        print("");print("");
        print("Main iteration starts");
        print("--------------");
    }


    int iter      = -1;
    int pindex    =  1;
    int terminate =  0;
    while( terminate == 0 )
    {    
        if( (sqrt(maxerr*maxerr) < 1.e-8) && (k >= 9) )
        {
            if(world.rank() == 0) { print(maxerr, k);}
            terminate = 1;
            world.gop.fence();
        }
        else
        {
            // Read-in of checkpointed files
            if( file_check == true)
            {    
                if(world.rank() == 0) { print(" "); print("Reading checkpoint");}

                char grecord_name2[100];
                sprintf(grecord_name2, "%s/g_record",direct);
                archive::ParallelInputArchive gin(world, grecord_name2, IO_nodes);
                gin & iter;
                gin & pindex;
                gin & maxerr;
                gin & maxerrold;
                gin & energy_n;
                gin & energy_p;
                gin & pot_n;
                gin & pot_p;
                gin & lap_n_old;
                gin & lap_p_old;
                gin.close();
                world.gop.fence();  // FENCE

                if(world.rank() == 0) { print("Reading neutrons");}
                char nrecord_name[100];
                sprintf(nrecord_name, "%s/n_record",direct);
                archive::ParallelInputArchive nind(world, nrecord_name, IO_nodes);
                for (unsigned int i=0; i< nu.size(); i++) { nind & nu[i]; nind & nd[i]; }
                nind.close();
                world.gop.fence();  // FENCE

                if(world.rank() == 0) { print("Reading protons");}
                char precord_name[100];
                sprintf(precord_name, "%s/p_record",direct);
                archive::ParallelInputArchive pind(world, precord_name, IO_nodes);
                for (unsigned int i=0; i< pu.size(); i++) { pind & pu[i]; pind & pd[i]; }
                pind.close();
                world.gop.fence();  // FENCE
                file_check = false;        
                if(world.rank() == 0) { print("Reading done"); print(" ");}

                //spin_split2(world, p, pu, pd);
                //spin_split2(world, n, nu, nd);
            }
            //END_TIMER("Reading Input");
        
            iter++;
            maxerr = -1.0;
            //double brad = -1.0;
            //double avg  = 0.40;
            double time = wall_time() + time_old;
            if(world.rank() == 0)
            {
                print("Iteration Nr.:   ", iter);
                print("wavelet number:  ", k);
                print("thresh:          ", thresh);
                print("Box size (2L):   ", 2.0*L);
                print("Averaging coeff: ", avg);
                if(brad > 0.0)
                {
                    print("Blurring radius: ", brad);
                }
                else { print("No blurring"); }
                print("wall time:       ", time);
                print(" ");
            }
        

            if( BE == 0.0)
            {
                den_p = densities(world, pu, pd);
                den_n = densities(world, nu, nd);
            }


            truncate2(world, pu, pd);
            truncate2(world, nu, nd);

        
            // Make Skyrme potential
            if(world.rank() == 0){print("Make potential");}
            if( iter > 0) 
            { 
                pot_p_old = copy(pot_p); 
                pot_n_old = copy(pot_n); 
            }

            //START_TIMER;
            Potential(world, pu, pd, nu, nd, energy_p, energy_n, pot_p, pot_n, 
                den_p, den_n, den, iter, BE, maxerr, tol, brad, pot_p_old, pot_n_old,
                lap_p_old, lap_n_old, avg);
            //END_TIMER("Potential");


            // Load-balance
            //START_TIMER;
            loadbalance_xxp(world, den_p, den_n, den);
            //END_TIMER("Loadbalance xxp");

            if(world.rank() == 0)
            {
                print(" ");
                print("Binding energy: ", BE); print(" ");
                print("Calculate new wave functions");
            }


            // Iterate neutrons
            if(world.rank() == 0){print("Neutrons: ");}
            //START_TIMER;
            iterate(world, pot_n, nu, nd, energy_n, maxerr, den, den_n, iter, avg, tol, brad);
            //END_TIMER("Iterate neutrons");

            // Iterate protons 
            if(world.rank() == 0){print("Protons: ");}
            //START_TIMER;
            iterate(world, pot_p, pu, pd, energy_p, maxerr, den, den_p, iter, avg, tol, brad);
            //END_TIMER("Iterate protons");

            world.gop.fence();

//            den_p.clear(); den_n.clear(); den.clear();                // IS
            pot_p_old.clear(); pot_n_old.clear();

            world.gop.fence();

            truncate2(world, pu, pd);
            truncate2(world, nu, nd);

            // Load-balance
            if(world.rank() == 0){print("Balance wf");}
            //START_TIMER;
            loadbalance_v4(world, pu, pd, nu, nd);
            //END_TIMER("loadbalance v4");

            // Print out maximum wavefunction change
            if(world.rank() == 0)
            {
                print(" "); print(" ");
                print("max_err(psi)=",maxerr);
            }

            // Output file to keep track of maxerr
            //START_TIMER;
            if( iter%1 == 0 || iter == 0)
            {    
                if(world.rank() == 0)
                {
                    outfile.open("log.txt", std::ios::app);
                    outfile << iter <<" "<< maxerr <<" "<< BE <<" "<< thresh <<" "<< time << std::endl;
                    outfile.close();
                }
            }

            //END_TIMER("Output outfile log.txt");
            world.gop.fence();  // FENCE gif

            BE_old = BE;


            // Project to higher wavelet number
            if( (maxerr <= 2.0 * A * thresh * prec && k<10) )
            {
                //START_TIMER;  // IS

                if( world.rank() == 0){ print("Project"); }

                thresh *= 1.e-1;
                k += 1;
                if ( k == 8 ) { brad = -1.0; }                                  // IS
                //brad = -1.00;
                //avg = 0.40;

                FunctionDefaults<3>::set_k(k);
                FunctionDefaults<3>::set_thresh(thresh);
                reconstruct(world, pu);
                reconstruct(world, pd);
                reconstruct(world, nu);
                reconstruct(world, nd);
                pot_n.reconstruct(); pot_p.reconstruct();
                //den_n.reconstruct(); den_p.reconstruct();                       // IS
                lap_n_old.reconstruct(); lap_p_old.reconstruct();

                for (unsigned int i=0; i<pu.size(); i++) 
                {
                    pu[i] = madness::project(pu[i], k, thresh);
                    pd[i] = madness::project(pd[i], k, thresh);
                }
                truncate2(world, pu, pd);
                normalize_ud(world, pu, pd);

                for (unsigned int i=0; i<nu.size(); i++)
                {
                    nu[i] = madness::project(nu[i], k, thresh);
                    nd[i] = madness::project(nd[i], k, thresh);
                }

                pot_n = madness::project(pot_n, k, thresh);
                pot_p = madness::project(pot_p, k, thresh);

                den_n = madness::project(den_n, k, thresh);                     // IS
                den_p = madness::project(den_p, k, thresh);                     // IS

                lap_n_old = madness::project(lap_n_old, k, thresh);
                lap_p_old = madness::project(lap_p_old, k, thresh);

                pindex++;
                world.gop.fence();  // FENCE gif, reorder 
                truncate2(world, nu, nd);
                normalize_ud(world, nu, nd);
                //END_TIMER("Reconstruct"); // IS
                den = den_p + den_n;                                            // IS
            }

            // Checkpointing
            //START_TIMER;  // IS

            if( iter%5 == 0 )
            {
                if( iter%10 == 0 && iter != 0)
                {
                    if( world.rank() == 0){ print(" "); print("Checkpoint"); print(" "); }
                    char trecord_name[100];
                    sprintf(trecord_name, "%s/tbck_record",direct);
                    archive::ParallelOutputArchive tout(world, trecord_name, IO_nodes);         // IS
                    tout & thresh;
                    tout & k;
                    tout & brad;
                    tout & avg;
                    tout & tol;
                    tout & time;
                    tout.close();
                    world.gop.fence();  // FENCE

                    char grecord_name[100];
                    sprintf(grecord_name, "%s/gbck_record",direct);
                    archive::ParallelOutputArchive gout(world, grecord_name, IO_nodes);         // IS
                    gout & iter;
                    gout & pindex;
                    gout & maxerr;
                    gout & maxerrold;
                    gout & energy_n;
                    gout & energy_p;
                    gout & pot_n;
                    gout & pot_p;
                    gout & lap_n_old;
                    gout & lap_p_old;
                    gout.close();
                    world.gop.fence();  // FENCE

                    char nrecord_name[100];
                    sprintf(nrecord_name, "%s/nbck_record",direct);
                    archive::ParallelOutputArchive noutc(world, nrecord_name, IO_nodes);        
                    for (unsigned int i=0; i<nu.size(); i++) { noutc & nu[i]; noutc & nd[i]; }
                    noutc.close();
                    world.gop.fence();  // FENCE

                    char precord_name[100];
                    sprintf(precord_name, "%s/pbck_record",direct);
                    archive::ParallelOutputArchive poutc(world, precord_name, IO_nodes);        
                    for (unsigned int i=0; i<pu.size(); i++) { poutc & pu[i]; poutc & pd[i]; }
                    poutc.close();
                    world.gop.fence();  // FENCE
                }
                else
                {
                    if( world.rank() == 0){ print(" "); print("Checkpoint"); print(" "); }
                    char trecord_name[100];
                    sprintf(trecord_name, "%s/t_record",direct);
                    archive::ParallelOutputArchive tout(world, trecord_name, IO_nodes);             // IS
                    tout & thresh;
                    tout & k;
                    tout & brad;
                    tout & avg;
                    tout & tol;
                    tout & time;
                    tout.close();
                    world.gop.fence();  // FENCE

                    char grecord_name[100];
                    sprintf(grecord_name, "%s/g_record",direct);
                    archive::ParallelOutputArchive gout(world, grecord_name, IO_nodes);         // IS
                    gout & iter;
                    gout & pindex;
                    gout & maxerr;
                    gout & maxerrold;
                    gout & energy_n;
                    gout & energy_p;
                    gout & pot_n;
                    gout & pot_p;
                    gout & lap_n_old;
                    gout & lap_p_old;
                    gout.close();
                    world.gop.fence();  // FENCE

                    char nrecord_name[100];
                    sprintf(nrecord_name, "%s/n_record",direct);
                    archive::ParallelOutputArchive noutc(world, nrecord_name, IO_nodes);
                    for (unsigned int i=0; i<nu.size(); i++) { noutc & nu[i]; noutc & nd[i];}

                    noutc.close();
                    world.gop.fence();  // FENCE

                    char precord_name[100];
                    sprintf(precord_name, "%s/p_record",direct);
                    archive::ParallelOutputArchive poutc(world, precord_name, IO_nodes);
                    for (unsigned int i=0; i<pu.size(); i++) { poutc & pu[i]; poutc & pd[i]; }
                    poutc.close();
                    world.gop.fence();  // FENCE
                }
            }
            //END_TIMER("Checkpointing");   // IS

            maxerrold = maxerr; 

            if(world.rank() == 0)
            {
                print("____________________________________________________");
                print("");
                print("");
            }
            world.gop.fence();
        }
    }
    world.gop.fence();
}




//Main Code
int main(int argc, char** argv)
{
  initialize(argc, argv);
  World world(SafeMPI::COMM_WORLD);
  try {
    startup(world,argc,argv);
    //print_meminfo(world.rank(), "startup");

    print("Deriv type:", DATA_PATH);

    FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap< Key<3> >(world)));
    std::cout.precision(8);
    
    doit(world);
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
  finalize();

  return 0;

}


