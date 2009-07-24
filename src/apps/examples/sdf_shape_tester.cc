/// Demonstrates the use of signed distance function shapes in
/// <mra/sdf_shape_3D.h> for integrating shapes in MADNESS simulations.
/// After constructing the shape, it prints out the shape using the plotvtk
/// functions.

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/sdf_shape_3D.h>

using namespace madness;

typedef Vector<double,3> coordT3d;
typedef Function<double,3> functionT;
typedef FunctionFactory<double,3> factoryT;
typedef std::vector<functionT> funcVecT;
typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;

static const double L = 2.0;
static const double Lplot = L;

static double thresh = 1e-3;
static double thresh1 = thresh*0.1;

//*****************************************************************************
// Shape tester.
int main(int argc, char **argv) {

	MPI::Init(argc, argv);
	World world(MPI::COMM_WORLD);
	startup(world,argc,argv);

	// Function defaults
	int k = 5;
	FunctionDefaults<3>::set_k(k);
	FunctionDefaults<3>::set_cubic_cell(-L, L);
	FunctionDefaults<3>::set_thresh(thresh);

	Tensor<int> bc(3,2);
	bc(_,0) = 0;          // Dirichlet in all directions
	bc(_,1) = 0;
	FunctionDefaults<3>::set_bc(bc);

	// create the shape mask
	coordT3d pt, vec;
	double c;
	pt[0] = 0.0;
	pt[1] = 0.5;
	pt[2] = 0.0;
	vec[0] = 0.0;
	vec[1] = 0.0;
	vec[2] = 1.0;
	c = 0.5;
	functionT mask = factoryT(world).functor(functorT(
		//new SDF_Cube<double>(0.2, thresh, sqrt(2.0), pt)
		//new SDF_Cone<double>(0.2, thresh, c, pt, vec)
		//new SDF_Paraboloid<double>(0.2, thresh, c, pt, vec)
		//new SDF_Plane<double>(0.2, thresh, vec, pt)
		//new SDF_Sphere<double>(0.2, thresh, c, pt)
		//new SDF_Ellipsoid<double>(0.2, thresh, vec, pt)
		//new SDF_Box<double>(0.2, thresh, vec, pt)
		new SDF_Cylinder<double>(0.2, thresh, 0.75, 1.0, pt, vec)
		));

	// the following line permutes the "inside" and "outside"
	//mask.unaryop(&mask_complement<double, 3>);

	// print the shape to a file
	char filename[100];
	sprintf(filename, "shape.vts");
	Vector<double, 3> plotlo, plothi;
	Vector<long, 3> npts;
	for(int i = 0; i < 3; ++i) {
		plotlo[i] = -Lplot;
		plothi[i] = Lplot;
		npts[i] = 51;
	}
	plotvtk_begin(world, filename, plotlo, plothi, npts);
	plotvtk_data(mask, "mask", world, filename, plotlo, plothi, npts);
	plotvtk_end<3>(world, filename);

	MPI::Finalize();

	return 0;
}
