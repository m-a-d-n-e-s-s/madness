
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
//	std::cout << "color, i in universe " << color << " " << i << std::endl;
	cloud.store(subworld,i,color);
	std::vector<T> j(nsubworld);;
	for (int icolor=0; icolor<nsubworld; ++icolor)
		cloud.load<T>(universe,j[icolor],icolor);	// loading is done by rank 0 of universe!
//	std::cout << "color, j in universe " << color << " " << j[color] << std::endl;
	T result=0;
	for (auto& i : j) result+=i;
	MADNESS_CHECK(std::abs(initial*factor-result)<1.e-13);
	if (universe.rank()==0) std::cout << "final result " << result << std::endl;
	return 0;
}



/// procedure for testing:
/// 1. make a universe world object, fill it with data
/// 2. store it to the cloud
/// 3. load it into the subworlds,
/// 4. multiply with the color of the subworld ( [0,..,nworld-1])
/// 5. store it in the cloud,
/// 6. load it into the universe
/// 7. accumulate.
/// result is value * (0 + 1 + 2 + .. + nworld-1)
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
		{
		// work done in subworlds
		dcT sub_initial(subworld);
		cloud.load(subworld,sub_initial,0);
//		dcT sub_initial=cloud.load<dcT>(subworld,0);
//		MADNESS_CHECK(sub_initial.size()==3);		// don't check this: it returns the *local* size
		for (auto& i : sub_initial) i.second*=color;
		cloud.store(subworld,sub_initial,color);
		subworld.gop.fence();
		}
	}
	subworld.gop.fence();		// <<<<<<--- this one is important: fence *after* destruction of sub_initial
	universe.gop.fence();

	// back to universe
	std::vector<dcT> vfinal(nsubworld);;
	for (int icolor=0; icolor<nsubworld; ++icolor) {
		vfinal[icolor]=dcT(universe);
		cloud.load<dcT>(universe,vfinal[icolor],icolor);	// loading is done by rank 0 of universe!
	}
//	std::cout << "color, j in universe " << color << " " << j[color] << std::endl;

	dcT result(universe);
	result.replace(2,T(0));
	result.replace(3,T(0));
	result.replace(5,T(0));

	for (auto& i : vfinal) {
		for (auto& element : i) {
			result.find(element.first).get()->second+=element.second;
			std::cout << "element.second " << element.second << std::endl;
		}
	}
	if (universe.rank()==0) {
		for (auto& element : result) std::cout << element.first << " " << element.second << std::endl;
	}
	std::cout << result.size() << " result size " << std::endl;
	for (auto& element : result) MADNESS_CHECK(std::abs(initial.find(element.first).get()->second*factor-element.second)<1.e-13);

	universe.gop.fence();
	return 0;
}

void test_cloud(madness::World& universe) {
	typedef madness::WorldContainer<int,double> dcT;

	madness::Cloud cloud(universe);
	cloud.set_debug(false);

	int success=0;
	// Make subworlds
	for (int isub=0; isub<universe.size(); ++isub) {
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
		success+=worlddc_test<int>(cloud,universe,subworld,color,nsubworld);
//		success+=worlddc_test<double>(cloud,universe,subworld,color,nsubworld);
	}


	// WorldObjects
//	{
//		dcT w(universe);
//		w.replace(1,1.0);
//		cloud.store(universe,w,0);
//		dcT w1=cloud.load<dcT>(subworld,0);
//		dcT w2(subworld);
//		cloud.load(subworld,w2,0);
//	}

//        {
//        	test_cloud_pod p(universe), p2(subworld);
//    		p.c=dcT(universe);
//    		cloud.store(universe,p,0);
//
//    		test_cloud_pod p1=cloud.load<test_cloud_pod>(subworld,0);
//    		cloud.load(subworld,p2,0);
//
//    		cloud.store(universe,std::make_tuple(color,p.c),1);
//    		int i3; dcT c3(subworld);
//    		auto p3=std::tie(i3,c3);
//    		cloud.load(subworld,p3,1);
//        }
	{

//
//		// prepare universe functions
//		dcT a(universe);
//		for (int i=0; i<10; ++i) a.replace(i,0);
//		std::cout << "a worldid " << a.get_world().id() << std::endl;
//		for (int i=0; i<10; ++i) std::cout << "a " << a.find(i).get()->first << " " << a.find(i).get()->second  << std::endl;
//		std::cout << "a.size (local to universe.rank()= " << universe.rank() <<": " << a.size() << std::endl;
//		cloud.store(universe,a,1);
//
//		// do subworld work
//		dcT sub_a=cloud.load<dcT>(subworld,1);	// sub_a is now replicated
//		for (int i=0; i<10; ++i) sub_a.replace(i,double(i)*(color+1));
//		for (int i=0; i<10; ++i) std::cout << "sub_a " << sub_a.find(i).get()->first << " " << sub_a.find(i).get()->second  << std::endl;
//		subworld.gop.fence();
//		cloud.store(subworld,sub_a,color);		// sub_a is stored on different records
//
//
//		// collect data into universe
//		dcT b0=cloud.load<dcT>(universe,0);		// b0 is from subworld 1, b1 is from subworld 1000
//		std::cout << "b0.size (local to universe rank" << universe.rank() <<") " << b0.size() << std::endl;
//		dcT b1=cloud.load<dcT>(universe,1);
//		std::cout << "b1.size (local to universe rank" << universe.rank() <<") " << b1.size() << std::endl;
//
//		universe.gop.fence();
//		for (int i=0; i<10; ++i) std::cout << "b0 " << b0.find(i).get()->first << " " << b0.find(i).get()->second  << std::endl;
//		universe.gop.fence();
//		for (int i=0; i<10; ++i) std::cout << "b1 " << b1.find(i).get()->first << " " << b1.find(i).get()->second  << std::endl;
//
//		universe.gop.fence();

	}
	universe.gop.fence();

}

void test(madness::World& world)
{
    madness::WorldContainer<long,std::vector<unsigned char>> container(world);
    {
        madness::archive::ContainerRecordOutputArchive ar(world,container,7);
        madness::archive::ParallelOutputArchive<madness::archive::ContainerRecordOutputArchive> par(world, ar);
        par & 1;
    }
    
    {
        madness::archive::ContainerRecordInputArchive ar(world,container,7);
        madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(world, ar);
        int value;
        par & value;
        std::cout << value << std::endl;
    }
    world.gop.fence();
}


void test_multi_work(madness::World& universe, madness::World& subworld,
		madness::WorldContainer<long,std::vector<unsigned char>>& container) {

	long key = 1+(universe.rank()%2);
    int value;
    {
        madness::WorldContainer<int,double> data(subworld);
        madness::archive::ContainerRecordInputArchive ar(subworld,container,key);
        madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(subworld, ar);
        par & value & data;
        std::cout << subworld.id() << " : subworld read process " << subworld.rank() << " read " << value << std::endl;
        std::cout << "data.size() " << subworld.id() << " " << subworld.rank() << " " << data.size() << std::endl;
    }
    {
        madness::WorldContainer<int,double> subdata(subworld);
        subdata.replace(subworld.rank(), subworld.rank()*100.0);
        subworld.gop.fence();
        
        long newkey = key*1000;
        madness::archive::ContainerRecordOutputArchive ar(subworld,container,newkey);
        madness::archive::ParallelOutputArchive<madness::archive::ContainerRecordOutputArchive> par(subworld, ar);
        par & value*1000 & subdata;
    }
}

void test_multi(madness::World& universe) {

    // For testing passing of containers from universe to subworld
    madness::WorldContainer<int,double> data(universe);
    data.replace(universe.rank(), universe.rank()*100.0);
    universe.gop.fence();
    int sz=data.size();
    std::cout << "universe data size in rank " << universe.rank() << " "  << sz << std::endl;
    universe.gop.sum(sz);
    std::cout << "universe data size total   "   << sz << std::endl;
    for (auto& a : data) {
    	std::cout << " universe content " << a.first << " " << a.second << std::endl;
    }

    // From universe put data into container
    madness::WorldContainer<long,std::vector<unsigned char>> container(universe);
    {
        madness::archive::ContainerRecordOutputArchive ar(universe,container,1);
        madness::archive::ParallelOutputArchive<madness::archive::ContainerRecordOutputArchive> par(universe, ar);
        par & 1 & data;
    }
    {
        madness::archive::ContainerRecordOutputArchive ar(universe,container,2);
        madness::archive::ParallelOutputArchive<madness::archive::ContainerRecordOutputArchive> par(universe, ar);
        par & 2 & data;
    }
    universe.gop.fence(); // Rank 0 must finish writing before workers start

    // Make subworlds
    int color = universe.rank() % 2;
    SafeMPI::Intracomm comm = universe.mpi.comm().Split(color, universe.rank() / 2);
    madness::World subworld(comm);
    
    // In subworlds read data from universe container and put results back
    test_multi_work(universe, subworld, container);
    subworld.gop.fence();
    universe.gop.fence(); // Workers must finish before reading the results

    // In universe read back results from subworld computations
    {
        madness::WorldContainer<int,double> subdata(universe);
        madness::archive::ContainerRecordInputArchive ar(universe,container,1000);
        madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(universe, ar);
        int value;
        par & value & subdata;
        int sz=subdata.size();
        universe.gop.fence(); // Workers must finish before reading the results
        madness::print("first result", value, sz);
        for (auto& a : subdata) {
        	std::cout << " universe content after work " << a.first << " " << a.second << std::endl;
        }

    }
    {
        madness::WorldContainer<int,double> subdata(universe);
        madness::archive::ContainerRecordInputArchive ar(universe,container,2000);
        madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(universe, ar);
        int value;
        par & value & subdata;
        int sz=subdata.size();
        madness::print("second result", value,sz);
    }
    universe.gop.fence();
}

int main(int argc, char** argv) {
    madness::World& world = madness::initialize(argc, argv);

    madness::archive::xxxtest(world);
//    test(world);
//    test_multi(world);
    test_cloud(world);
    
    madness::finalize();

    return 0;
}
