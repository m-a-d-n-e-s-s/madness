
//#define MAD_ARCHIVE_DEBUG_ENABLE
#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/world/MADworld.h>
#include <madness/world/worlddc.h>
#include <madness/world/parallel_dc_archive.h>
#include <madness/world/cloud.h>


/// helper type for storing data
struct test_cloud_pod {
	int i;
	double d;
	madness::WorldContainer<int,double> c;
	test_cloud_pod(madness::World& world) : i(), d(), c(world) {}
	template<typename Archive>
	void serialize(Archive& ar) {
		ar & i & d &c;
	}
};


/// procedure for testing:
/// 1. make a universe object with an initial value
/// 2. store it to the cloud
/// 3. load it into the subworlds,
/// 4. multiply with the color of the subworld ( [0,..,nworld-1])
/// 5. store it in the cloud,
/// 6. load it into the universe
/// 7. accumulate.
/// result is value * (0 + 1 + 2 + .. + nworld-1)

/// store a simple type (e.g. int, double, i.e. no distributed container)
template<typename T>
int simple_test(madness::Cloud& cloud, madness::World& universe, madness::World& subworld,
		int color,int nsubworld) {

	T initial(2);
	double factor = 0;
	for (int i=0; i<nsubworld; ++i) factor+=double(i);

	cloud.store(universe,initial,0);

	T i=cloud.load<T>(subworld,0);
	MADNESS_CHECK(i==initial);
	i*=color;
	subworld.gop.fence();
	cloud.store(subworld,i,color+1);
	universe.gop.fence();

	std::vector<T> j(nsubworld);;
	for (int icolor=0; icolor<nsubworld; ++icolor)
		cloud.load<T>(universe,j[icolor],icolor+1);	// loading is done by rank 0 of universe!
//	std::cout << "color, j in universe " << color << " " << j[color] << std::endl;
	T result=0;
	for (auto& i : j) result+=i;
	T err=std::abs(initial*factor-result);
//	if (universe.rank()==0) std::cout << "final result " << result << ", error " <<err<< std::endl;
	MADNESS_CHECK(std::abs(err)<1.e-10);
	return 0;
}



/// store a simple type (e.g. int, double, i.e. no distributed container)
template<typename T>
int simple_test_vector(madness::Cloud& cloud, madness::World& universe, madness::World& subworld,
		int color,int nsubworld) {

	std::vector<T> initial={T(2),T(3),T(5)};
	T initial_sum=T(2)+T(3)+T(5);
	double factor = 0;
	for (int i=0; i<nsubworld; ++i) factor+=double(i);
	cloud.store(universe,initial,0);

	std::vector<T> i;
	cloud.load<std::vector<T> >(subworld,i,0);
	MADNESS_CHECK(i==initial);
	for (auto& element : i) element*=color;
	subworld.gop.fence();
	cloud.store(subworld,i,color+1);
	universe.gop.fence();

	std::vector<std::vector<T> > j(nsubworld);;
	for (int icolor=0; icolor<nsubworld; ++icolor)
		cloud.load<std::vector<T> >(universe,j[icolor],icolor+1);	// loading is done by rank 0 of universe!
//	std::cout << "color, j in universe " << color << " " << j[color] << std::endl;
	T result=0;
	for (auto& i : j) {
		for (auto& element : i) result+=element;
	}
	T err=std::abs(initial_sum*factor-result);
//	if (universe.rank()==0) std::cout << "final result " << result << ", error " <<err<< std::endl;
	MADNESS_CHECK(std::abs(err)<1.e-10);
	return 0;
}


/// test distributed containers
template<typename T>
int worlddc_test(madness::Cloud& cloud, madness::World& universe, madness::World& subworld,
		int color,int nsubworld) {

	typedef madness::WorldContainer<int,T> dcT;

	double factor = 0;
	for (int i=0; i<nsubworld; ++i) factor+=double(i);

	dcT initial(universe);
	initial.replace(2,T(2));
	initial.replace(3,T(3));
	initial.replace(5,T(5));
	cloud.store(universe,initial,0);
	universe.gop.fence();

	{
		// work done in subworlds
		dcT sub_initial(subworld);
		cloud.load(subworld,sub_initial,0);
//		dcT sub_initial=cloud.load<dcT>(subworld,0);
//		MADNESS_CHECK(sub_initial.size()==3);		// don't check this: it returns the *local* size
		for (auto& i : sub_initial) i.second*=color;
		cloud.store(subworld,sub_initial,color+1);
		subworld.gop.fence();
	}
	subworld.gop.fence();		// <<<<<<--- this one is important: fence *after* destruction of sub_initial
	universe.gop.fence();

	// back to universe
	std::vector<dcT> vfinal(nsubworld);;
	for (int icolor=0; icolor<nsubworld; ++icolor) {
		vfinal[icolor]=dcT(universe);
		cloud.load<dcT>(universe,vfinal[icolor],icolor+1);	// loading is done by rank 0 of universe!
	}
//	std::cout << "color, j in universe " << color << " " << j[color] << std::endl;

	dcT result(universe);
	result.replace(2,T(0));
	result.replace(3,T(0));
	result.replace(5,T(0));
	universe.gop.fence();

	for (int i : {2,3,5}) {
        typename dcT::accessor acc;
        if (result.probe(i)) {	// check if local
			MADNESS_CHECK(result.find(acc,i));
			for (auto& f : vfinal) {
				acc->second+=f.find(i).get()->second;
			}
        }
	}
	universe.gop.fence();

	for (auto& element : result) {
		MADNESS_CHECK(std::abs(initial.find(element.first).get()->second*factor-element.second)<1.e-13);
	}

	universe.gop.fence();
	return 0;
}


/// test distributed containers, vector version
template<typename T>
int worlddc_test_vector(madness::Cloud& cloud, madness::World& universe, madness::World& subworld,
		int color,int nsubworld) {

	typedef madness::WorldContainer<int,T> dcT;
	const int nvec=3;

	auto init_dcTvector = [&nvec] (madness::World& w) {
		std::vector<dcT> initial(nvec);
		for (auto &  i : initial) i=dcT(w);
		return initial;
	};

	double factor = 0;
	for (int i=0; i<nsubworld; ++i) factor+=double(i);

	std::vector<dcT> initial=init_dcTvector(universe);
	for (auto &  i : initial) {
		i.replace(2,T(2));
		i.replace(3,T(3));
		i.replace(5,T(5));
	}
	cloud.store(universe,initial,0);
	universe.gop.fence();

	{
		std::vector<dcT> sub_initial;
		cloud.load(subworld,sub_initial,0);

		// work done in subworlds
		for (auto& i : sub_initial) {
			for (auto& element : i) {
				element.second*=color;
			}
		}
		cloud.store(subworld,sub_initial,color+1);
		subworld.gop.fence();
	}
	subworld.gop.fence();		// <<<<<<--- this one is important: fence *after* destruction of sub_initial
	universe.gop.fence();

	// back to universe: load each subworld vector into a vector element
	std::vector<std::vector<dcT> > vfinal(nsubworld);

	for (int icolor=0; icolor<nsubworld; ++icolor) {
		cloud.load<std::vector<dcT> >(universe,vfinal[icolor],icolor+1);	// loading is done by rank 0 of universe!
	}

	std::vector<dcT> result=init_dcTvector(universe);
	for (auto &  r : result) {
		r.replace(2,T(0));
		r.replace(3,T(0));
		r.replace(5,T(0));
	}
	universe.gop.fence();

	// loop over vector elements
	for (auto& i : vfinal) {		// loop over subworld results
		for (int ivec=0; ivec<nvec; ++ivec) {
			for (auto& element : i[ivec]) {
				result[ivec].find(element.first).get()->second+=element.second;
			}
		}
	}

	universe.gop.fence();
	for (int ivec=0; ivec<nvec; ++ivec) {
		for (auto& element : result[ivec]) {
			MADNESS_CHECK(std::abs(initial[ivec].find(element.first).get()->second*factor-element.second)<1.e-13);
		}
	}

	universe.gop.fence();
	return 0;
}



/// store a simple type (e.g. int, double, i.e. no distributed container)
int pod_test(madness::Cloud& cloud, madness::World& universe, madness::World& subworld,
		int color,int nsubworld) {

	test_cloud_pod initial(universe);
	initial.i=2;
	initial.d=3.0;
	double factor = 0;
	for (int i=0; i<nsubworld; ++i) factor+=double(i);
	cloud.store(universe,initial,0);

	{
		test_cloud_pod i=cloud.load<test_cloud_pod>(subworld,0);
		MADNESS_CHECK(i.i==initial.i);
		i.d*=color;
		i.i*=color;

		subworld.gop.fence();
		cloud.store(subworld,i,color+1);
	}
	subworld.gop.fence();
	universe.gop.fence();

	std::vector<test_cloud_pod> j(nsubworld,test_cloud_pod(universe));
	for (int icolor=0; icolor<nsubworld; ++icolor)
		cloud.load<test_cloud_pod>(universe,j[icolor],icolor+1);	// loading is done by rank 0 of universe!
//	std::cout << "color, j in universe " << color << " " << j[color] << std::endl;
	test_cloud_pod result(universe);
	for (auto& i : j) {
		result.i+=i.i;
		result.d+=i.d;
	}
	double err=std::abs(initial.i*factor-result.i) + std::abs(initial.d*factor-result.d);
//	if (universe.rank()==0) std::cout << "final result " << result << ", error " <<err<< std::endl;
	MADNESS_CHECK(std::abs(err)<1.e-10);
	return 0;
}


int test_cloud(madness::World& universe) {
	typedef madness::WorldContainer<int,double> dcT;

	madness::Cloud cloud(universe);
	cloud.set_debug(false);

	int success=0;
	for (int isub=0; isub<universe.size(); ++isub) {

		// Make subworlds
		int nsubworld=isub+1;
		int color = universe.rank() % nsubworld;
		SafeMPI::Intracomm comm = universe.mpi.comm().Split(color, universe.rank() / nsubworld);

		std::cout << "created " << nsubworld << " subworlds, color " << color << std::endl;

		madness::World subworld(comm);
		universe.gop.fence();
		subworld.gop.fence();

		success+=simple_test<int>(cloud,universe,subworld,color,nsubworld);
		success+=simple_test<double>(cloud,universe,subworld,color,nsubworld);
		success+=simple_test<std::complex<double> >(cloud,universe,subworld,color,nsubworld);
		if (universe.rank()==0) std::cout << "passed simple tests " << std::endl;

		success+=simple_test_vector<int>(cloud,universe,subworld,color,nsubworld);
		success+=simple_test_vector<double>(cloud,universe,subworld,color,nsubworld);
		success+=simple_test_vector<std::complex<double> >(cloud,universe,subworld,color,nsubworld);
		if (universe.rank()==0) std::cout << "passed simple vector tests " << std::endl;

		success+=worlddc_test<int>(cloud,universe,subworld,color,nsubworld);
		success+=worlddc_test<double>(cloud,universe,subworld,color,nsubworld);
		if (universe.rank()==0) std::cout << "passed dc tests " << std::endl;

		success+=worlddc_test_vector<int>(cloud,universe,subworld,color,nsubworld);
		success+=worlddc_test_vector<double>(cloud,universe,subworld,color,nsubworld);
		if (universe.rank()==0) std::cout << "passed dc vector tests " << std::endl;

		success+=pod_test(cloud,universe,subworld,color,nsubworld);
		if (universe.rank()==0) std::cout << "passed pod tests " << std::endl << std::endl;

	}
	universe.gop.fence();
	return success;
}


int main(int argc, char** argv) {
    madness::World& world = madness::initialize(argc, argv);

    int success=0;

	int nrepeat=3;
	for (int repeat=0; repeat<nrepeat; ++repeat) success+=test_cloud(world);
    
    madness::finalize();

    return 0;
}
