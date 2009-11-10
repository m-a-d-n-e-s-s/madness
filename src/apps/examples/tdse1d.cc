/// \file tdse1d.cc
/// \brief Example propagation of TDSE (translating atom) using various propagators


#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/qmprop.h>
#include <mra/operator.h>
#include <constants.h>
#include <string>

using namespace madness;

// typedefs to make life less verbose
typedef Vector<double,1> coordT;
typedef SharedPtr< FunctionFunctorInterface<double,1> > functorT;
typedef Function<double,1> functionT;
typedef FunctionFactory<double,1> factoryT;
typedef SeparatedConvolution<double,1> operatorT;

typedef SharedPtr< FunctionFunctorInterface<double_complex,1> > complex_functorT;
typedef Function<double_complex,1> complex_functionT;
typedef FunctionFactory<double_complex,1> complex_factoryT;
typedef SeparatedConvolution<double_complex,1> complex_operatorT;
typedef SharedPtr<complex_operatorT> pcomplex_operatorT;
 
// Simulation parameters
static const double L = 100.0; // Simulation in [-L,L]
static const double x0 = -L + 10.0; // Initial position of the atom
static const double energy_exact = -6.188788775728796797594788; // From Maple
static const long k = 16;        // wavelet order
static const double thresh = 1e-12; // precision
static const double velocity = 3.0;
//static const double eshift = energy_exact - 0.5*velocity*velocity; // Use this value to remove rotating phase
static const double eshift = 0.0;

static double current_time = 0.0; // Lazy but easier than making functors for everything

/////////////////////////////////// For quadrature rules///////////////////////////////
// global vars for the laziness
static const double_complex I = double_complex(0,1);
vector<double> B, tc;
pcomplex_operatorT G;
vector<pcomplex_operatorT> Gs, Gss;
const int maxiter = 20;
const double fix_iter_tol = 1e-5;
int np=3; // number of quadrature pts


// Position of atom at current time
double atom_position() {
    return x0 + velocity*current_time;
}

// Exact solution ... (H-E)psi is accurate to 2e-7 or better inside Maple
double_complex psi_exact(const coordT& r) {
    const double x = r[0] - atom_position();

    if (fabs(x) > 9.0) return 0.0;

    const double xsq = x*x;
    
    // Horner's form for stability ... yes, it is 70-order polyn ... don't panic ... all is OK
    const double psi = exp(-1.30*xsq)*(-1.02151632756275513018719494826+(.522210612113707231536059971069+(-.378478352719362210706386739834+(.183732263756009855463432582593+(-0.866826311335724362186706464428e-1+(0.364601910940641762284283418688e-1+(-0.144289291226899801775738667691e-1+(0.536464813679927807734520149659e-2+(-0.188945345474975346004237994967e-2+(0.628725522158030920219217207400e-3+(-0.195986657875763402765072721284e-3+(0.563993909330309881048156461300e-4+(-0.147273758530730072646826354339e-4+(0.343202525037692780236348165792e-5+(-7.03765391498970506694108123682e-7+(1.25577395245191089173671652503e-7+(-1.93270666918809750376161513191e-8+(2.54624395753990033043923369519e-9+(-2.84968491109847694778452937732e-10+(2.68398879018928855922002181223e-11+(-2.09811331054703124038096336404e-12+(1.32869596552058891546970129944e-13+(-6.47453843054578193348503832223e-15+(2.08146181239264250219157355910e-16+(-8.27692933283146365185698371000e-19+(-4.21076100567125673420604632290e-19+(3.34873683682880953223636900553e-20+(-1.62840449033428572820374528252e-21+(5.80781234060753332169593869292e-23+(-1.58754665363960354625758075480e-24+(3.34792415609741263184593450566e-26+(-5.37727687523701580992550859153e-28+(6.37706272777548158355242766766e-30+(-5.27154378837014394526580093435e-32+(2.71430533655911926094030668033e-34-6.55694230766452931026127498540e-37*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq);

    // Galilean translation factor
    const double arg = x*velocity - (energy_exact - velocity*velocity*0.5 - eshift)*current_time;
    const double_complex tranfac = exp(double_complex(0,arg));
    return psi*tranfac;
}

// Time-dependent potential for translating atom
double V(const coordT& r) {
    const double x = r[0] - atom_position();
    if (fabs(x) > 6.2) return 0.0;

    return -8.0*exp(-x*x) - eshift;
}

// (dV/dr)^2
double dVsq(const coordT& r) {
    const double x = r[0] - atom_position();
    if (fabs(x) > 4.5) return 0.0;

    double dv = 16.0*x*exp(-x*x);
    return dv*dv;
}
    

template<typename T, int NDIM>
struct unaryexp {
    void operator()(const Key<NDIM>& key, Tensor<T>& t) const {
        UNARY_OPTIMIZED_ITERATOR(T, t, *_p0 = exp(*_p0););
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};

// Uterm means include extra bit for Chin-Chen central potential
complex_functionT expV(World& world, const double t, bool Uterm=false) {
    functionT potn = factoryT(world).f(V);
    if (Uterm) {
        double delta = 1.5*t;
        functionT d = factoryT(world).f(dVsq).initial_level(10);
        potn.compress(); d.compress();
        potn.gaxpy(1.0, d, -delta*delta/48.0);
    }
    complex_functionT expV = double_complex(0.0,-t)*potn;
    expV.unaryop(unaryexp<double_complex,1>());
    expV.truncate();
    return expV;
}

// Evolve forward one time step using Chin-Chen
complex_functionT chinchen(World& world, const complex_operatorT& G, const complex_functionT& psi0, const double tstep) {
    complex_functionT psi;
    
    // G = exp(-itV(t)/6) G0(t/2) exp(-2itU(t/2)/3) G0(t/2) exp(-itV(0)/6)

    // U = V - t^2*(del V)^2 / 48
    
    psi = expV(world, tstep/6.0)*psi0;          psi.truncate();
    psi = apply(G, psi);                       psi.truncate();

    current_time += 0.5*tstep;
    
    psi = expV(world, 2.0*tstep/3.0, true)*psi; psi.truncate();

    current_time += 0.5*tstep;

    psi = apply(G,psi);                         psi.truncate();
    psi = expV(world, tstep/6.0)*psi;           psi.truncate();
    
    return psi;
}


// Evolve forward one time step using Trotter ... G = G0(tstep/2)
complex_functionT trotter(World& world, const complex_operatorT& G, const complex_functionT& psi0, const double tstep) {
    complex_functionT psi = apply(G, psi0);    psi.truncate();

    current_time += 0.5*tstep;

    psi = expV(world, tstep)*psi;              psi.truncate();

    current_time += 0.5*tstep;

    psi = apply(G,psi);                        psi.truncate();
    
    return psi;
}

void print_info(World& world, const complex_functionT& psi, int step) {
    functionT potn = factoryT(world).f(V).truncate_on_project();
    complex_functionT dpsi = diff(psi,0);
    double ke = inner(dpsi,dpsi).real() * 0.5;
    double pe = psi.inner(psi*potn).real();
    double norm = psi.norm2();

    complex_functionT psiX = complex_factoryT(world).f(psi_exact);
    double err = (psiX - psi).norm2();

    if ((step%40) == 0) {
        printf("\n");
        printf(" step    time      atom x        norm        kinetic    potential      energy      err norm   depth   size  \n");
        printf("------ -------- ------------ ------------ ------------ ------------ ------------ ------------ ----- ---------\n");
    }
    printf("%6d %8.2f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %4d %9ld\n",
           step, current_time, atom_position(), norm, ke, pe, ke+pe, err, int(psi.max_depth()), psi.size());
}


static void readin(int np) {
	B.resize(np);
	tc.resize(np);
	
	if (np==1) {
		tc[0] = 0.5;
		B[0] = 1;
	} else if (np==2) {
		tc[0] = .2113248654051871; 
		tc[1] = .7886751345948129; 
		B[0] = .5; 
		B[1] = .5;
	} else if (np==3) {
		tc[0] = .11270166537925831; 
		tc[1] = .5; 
		tc[2] = .88729833462074169; 
		B[0]  = .277777777777777777777777777778; 
		B[1]  = .444444444444444444444444444444; 
		B[2]  = .277777777777777777777777777778;
	} else if (np==4) {
		tc[0] = .06943184420297371; 
		tc[1] = .33000947820757187; 
		tc[2] = .66999052179242813; 
		tc[3] = .93056815579702629;
		B[0]  = .173927422568726928686531974611; 
		B[1]  = .326072577431273071313468025389;
		B[2]  = .326072577431273071313468025389; 
		B[3]  = .173927422568726928686531974611;
	} else {
		print("invalid np, throw now");
		throw;
	}
}

//j_th interpolating coefficients
double icoeff(const int np, const int j, const double t) {
	double dum = 1;
	for (int i=  0; i< j; ++i) dum *= (t-tc[i])/(tc[j]-tc[i]);
	for (int i=j+1; i<np; ++i) dum *= (t-tc[i])/(tc[j]-tc[i]);
	return dum;
}

//interpolate ps at t
template<typename T> T myp(const vector<T>& ps, const double t) {
	int np = ps.size();
	T p = ps[0]*icoeff(np, 0, t);
	for (int j=1; j<np; ++j) p += ps[j]*icoeff(np, j, t);
	return p;
}

// Evolve forward one time step using quadrature rules
complex_functionT q_r(World& world, const int np, const complex_functionT psi0, const double tstep) {
    //can be more vectorized.
	
	vector<complex_functionT> ps(np), ps1(np);
    for (int i=0; i<np; ++i) ps[i] = copy(psi0);
	
	
	double tdum = current_time;
	complex_functionT pdum;
	vector<complex_functionT> qs(np);
	vector<functionT> Vs(np);
	vector< vector<functionT> > Vss(np);
	for (int i=0; i<np; ++i) {
		current_time = tdum + tstep*tc[i];
		qs[i] = apply(*Gs[np - i - 1], psi0).truncate();
		Vs[i] = factoryT(world).f(V).truncate_on_project();
		Vs[i].truncate();
		Vss[i].resize(np);
		for (int k=0; k<np; ++k) {
			current_time = tdum + tstep*tc[i]*tc[k];
			Vss[i][k]=factoryT(world).f(V).truncate_on_project();
			Vss[i][k].truncate();
		}
	}
	
	// fix pt iterations for psi's on the quadrature pts
	printf("fix iters");
	double err;
	for (int j=0; j<maxiter; ++j) {
		err = 0;
		for (int i=0; i<np; ++i) {
			ps1[i] = copy(qs[i]);
			for (int k=0; k<np; ++k) ps1[i] -= apply(*Gss[i*np+k], Vss[i][k]*myp(ps,tc[i]*tc[k])).scale(tstep*tc[i]*I*B[k]);
			err += (ps1[i].truncate()-ps[i]).norm2();
		}
		err /= np;
		ps = copy(world, ps1);
		printf(" %6.0e",err);
		if (err <= fix_iter_tol) break;
	}
	printf("\n");
			
    // apply quadrature rule.
	pdum = apply(*G, psi0).truncate();
	
	for (int k=0; k<np; ++k) pdum -= apply(*Gs[k], Vs[k]*ps[k]).scale(I*B[k]*tstep);
	
	current_time = tdum + tstep;

    return pdum.truncate();
}


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);
    cout.precision(8);

    FunctionDefaults<1>::set_k(k);                 // Wavelet order
    FunctionDefaults<1>::set_thresh(thresh);       // Accuracy
    FunctionDefaults<1>::set_autorefine(false);
    FunctionDefaults<1>::set_cubic_cell(-L,L);
    FunctionDefaults<1>::set_initial_level(8);     // Initial projection level
    FunctionDefaults<1>::set_truncate_mode(0);

    // Estimate the bandwidth and largest practical time step using G0
    double ctarget = 20.0; // Estimated from FT of potential
    double c = 1.86*ctarget;
    double tcrit = 2*constants::pi/(c*c);

    print("Critical time step is", tcrit);

    complex_functionT psi = complex_factoryT(world).f(psi_exact);
    psi.truncate();

    string method("chinchen");

    if (method == "trotter") {
        double tstep = tcrit*3.0;
        int nstep = velocity==0 ? 100 : (L - 10 - x0)/velocity/tstep;
        print("No. of time steps is", nstep);
	complex_operatorT G0 = qm_free_particle_propagator<1>(world, k, c, 0.5*tstep);
        for (int step=0; step<nstep; step++) {
            print_info(world, psi,step);
            psi = trotter(world, G0, psi, tstep);
        } 
    }
    else if (method == "chinchen") {
        double tstep = tcrit*10.0;
        print("tstep", tstep);
        int nstep = velocity==0 ? 100 : (L - 10 - x0)/velocity/tstep;
        print("No. of time steps is", nstep);
	complex_operatorT G0 = qm_free_particle_propagator<1>(world, k, c, 0.5*tstep);
        for (int step=0; step<nstep; step++) {
            print_info(world, psi,step);
            psi = chinchen(world, G0, psi, tstep);
        } 
    }
    else {
        double tstep = tcrit*12;
	int nstep = velocity==0 ? 100 : (L - 10 - x0)/velocity/tstep;
        print("No. of time steps is", nstep);   
	readin(np);
	G = pcomplex_operatorT(qm_free_particle_propagatorPtr<1>(world, k, c, tstep));
        for (int i=0; i<np; ++i) Gs.push_back(pcomplex_operatorT(qm_free_particle_propagatorPtr<1>(world, k, c, (1-tc[i])*tstep)));
        for (int j=0; j<np; ++j) 
            for (int i=0; i<np; ++i) 
                Gss.push_back(pcomplex_operatorT(qm_free_particle_propagatorPtr<1>(world, k, c, (1-tc[i])*tstep*tc[j]))); //[j*np+i]

        for (int step=0; step<nstep; step++) {
            print_info(world, psi, step);
            psi = q_r(world, np, psi, tstep);
        }
    }
    
    finalize();
    return 0;
}

