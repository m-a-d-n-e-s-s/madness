/*
 * test_macrotask.cc
 *
 *  Created on: Dec 17, 2019
 *      Author: fbischoff
 */
//#define WORLD_INSTANTIATE_STATIC_TEMPLATES

#include <madness/mra/mra.h>
#include <iostream>
#include <madness/world/MADworld.h>
#include <madness/world/worlddc.h>
#include <random>
#include <madness/mra/funcimpl.h>

using namespace madness;
using namespace archive;

/**
 * 	Issues:
 * 	 - set_defaults<>(local_world) for loading (affects pmap)
 * 	 - serialization of task works only for int, double, .. but not for Function
 *   - turn data structure into tuple
 *   - prioritize tasks
 *
 */



double gaussian(const coord_4d& r) {
    double x=r[0], y=r[1], z=r[2], aa=r[3];
    return exp(-(x*x + y*y + z*z * aa*aa))*abs(sin(abs(2.0*x))) *cos(y);
}

/// for each process create a world using a communicator shared with other processes by round-robin
/// copy-paste from test_world.cc
std::shared_ptr<World> create_worlds(World& universe, const std::size_t nworld) {

	if (universe.size()<nworld) {
		print("trying to create ",nworld,"world with",universe.size(), "processes");
		MADNESS_EXCEPTION("increase number of processes",1);
	}

	if (universe.rank()==0) print("== multiple worlds created with Intracomm::Create()==",nworld);
    std::vector<std::vector<int> > process_list(nworld);
    std::shared_ptr<World> all_worlds;

	for (int i=0; i<universe.size(); ++i) process_list[i%nworld].push_back(i);
	if (universe.rank()<nworld) print("process_list",process_list[universe.rank()]);


	for (int i=0; i<process_list.size(); ++i) {
		const std::vector<int>& pl=process_list[i];
		bool found=(std::find(pl.begin(),pl.end(),universe.rank())!=pl.end());
		if (found) {
			print("assigning rank",universe.rank(),"to world group",pl);

			SafeMPI::Group group = universe.mpi.comm().Get_group().Incl(pl.size(), &pl[0]);
			SafeMPI::Intracomm comm_group = universe.mpi.comm().Create(group);

			all_worlds.reset(new World(comm_group));
		}
	}
	universe.gop.fence();
	return all_worlds;
}

/// base class
template<typename resultT, typename dataT>
class macro_task_base {
public:
	virtual resultT run(World& world, const dataT& data) = 0;
	virtual ~macro_task_base() {};

	double priority=0.0;

};

template<typename T, std::size_t NDIM>
Function<T,NDIM> localize(World& origin, World& destination, const Function<T,NDIM>& data, const long id) {

	origin.gop.fence();
	destination.gop.fence();
	std::string filename="smartie"+std::to_string(id);
	if (data.is_initialized() and data.world().id()==origin.id()) save(data,filename);

	destination.gop.fence();
	origin.gop.fence();

    auto pmap = std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >(new madness::LevelPmap< Key<NDIM> >(destination));
    FunctionDefaults<NDIM>::set_pmap(pmap);	// set default pmap to use only this world!
	Function<T,NDIM> result(destination);
	load(result,filename);

	destination.gop.fence();
	origin.gop.fence();

    archive::ParallelInputArchive ar2(destination, filename.c_str(), 1);
    ar2.remove();

	return result;
}

template<typename T, std::size_t NDIM>
struct data_type {
	typedef Function<T,NDIM> functionT;
	data_type() : i(), d(), f() {}
	data_type(const int& i, const double& d, Function<T,NDIM> f) : i(i), d(d), f(f),
			filename("dummy"+std::to_string(i)) {}
	data_type(const int& i, const double& d) : i(i), d(d), f(),
			filename("dummy"+std::to_string(i)) {}
//	data1(World& world) : i(), d(), f(Function<T,NDIM>) {}
	double d;
	int i;
	Function<T,NDIM> f;
	std::string filename="dummy";


    template <typename Archive>
    void serialize(const Archive& ar) {
    	bool fexist=f.is_initialized();
        ar & i & d & filename;
//		if (fexist) ar & f;

    }
	void save(World& world) const {
		world.gop.fence();
		archive::ParallelOutputArchive ar(world, filename.c_str() , 1);
//		print("saving to file",filename);
		ar & d & i & f;
	}

	void load(World& world) {
		world.gop.fence();
        auto pmap = std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >(new madness::LevelPmap< Key<NDIM> >(world));
        FunctionDefaults<NDIM>::set_pmap(pmap);	// set default pmap to use only this world!
//		print("loading from file",filename, world.id());
		archive::ParallelInputArchive ar(world, filename.c_str() , 1);
		ar & d & i & f;
	}

};



template<typename resultT, typename dataT>
class macro_task : public macro_task_base<resultT, dataT> {
public:
	typedef resultT result_type;
	typedef dataT data_type;

	macro_task() {}


	resultT run(World& world, const dataT& data) {
		const Function<double,4>& f=data.f;
		Function<double,4> g=real_factory_4d(world).f(gaussian);
		Function<double,4> f2=square(f)+g;
		return f2;
	}

	~macro_task() {}


    template <typename Archive>
    void serialize(const Archive& ar) {
        ar & this->priority;
    }
	void save(World& world, const dataT& data) const {
		data.save(world);
	}

	void load(World& world, dataT& data) {
		data.load(world);
	}

};



class MasterPmap : public WorldDCPmapInterface<long> {
public:
    MasterPmap() {}
    ProcessID owner(const long& key) const {return 0;}
};



template<typename taskT>
class macro_taskq : public WorldObject< macro_taskq<taskT> > {
    typedef macro_taskq<taskT> thistype;
    typedef WorldContainer<long,taskT> dcT;
    typedef typename taskT::result_type resultT;
    typedef typename taskT::data_type dataT;

    World& universe;
    std::shared_ptr<World> regional_ptr;
	WorldContainer<long,taskT> taskq;


public:

	World& get_regional() {return *regional_ptr;}

    /// create an empty taskq and initialize the regional world groups
	macro_taskq(World& universe, int nworld)
		  : universe(universe), WorldObject<thistype>(universe),
			taskq(universe,std::shared_ptr< WorldDCPmapInterface<long> > (new MasterPmap())) {

		regional_ptr=(create_worlds(universe,nworld));
	    World& regional=*(regional_ptr.get());
		this->process_pending();
		taskq.process_pending();
	}

	std::vector<resultT> run_all(taskT& task, const std::vector<dataT>& data) {
		std::vector<resultT> result(data.size());
		for (int i=0; i<data.size(); ++i) {
			add_task(std::make_tuple(i,task));
			data[i].save(universe);
		}

		World& regional=get_regional();

		bool working=true;
		while (working) {
			auto [key,task] =get_task_from_tasklist(regional);

			if (key>=0) {
				dataT regionaldata(data[key]);
				regionaldata.load(regional);

				double cpu0=cpu_time();
				result[key]=task.run(regional,regionaldata);
				double cpu1=cpu_time();
				printf("finished task %ld in %4.2fs\n",key,cpu1-cpu0);
			} else {
				working=false;
			}
		}
		universe.gop.fence();
		regional.gop.fence();

	    for (int i=0; i<result.size(); ++i) {
	    	result[i]=localize(regional,universe,result[i],i);
	    }

		return result;

	}

	void add_task(const std::tuple<long,taskT>& task) {
		int key=std::get<0>(task);
		if (universe.rank()==0) taskq.replace(key,std::get<1>(task));
	};

	std::pair<long,taskT> get_task_from_tasklist(World& regional) {

		// only root may pop from the task list
		std::pair<long,taskT> task;
		if (regional.rank()==0) task=pop();
		regional.gop.broadcast_serializable(task, 0);
		regional.gop.fence();
		return task;
	}

	std::size_t size() const {
		return taskq.size();
	}

	Future<std::pair<long,taskT> > pop() {
        return this->task(ProcessID(0), &macro_taskq<taskT>::pop_local, this->get_world().rank());
	}

	std::pair<long,taskT> pop_local(const int requested) {

		taskT result;
		long key;
		bool found=false;
		while (not found) {
			auto iter=taskq.begin();
			if (iter==taskq.end()) {
				print("taskq empty");
				return std::make_pair(-1,taskT());
			}

			key=iter->first;

			typename dcT::accessor acc;
			if (taskq.find(acc,key)) {
				result=acc->second;
				taskq.erase(key);
				break;
			} else {
				print("could not find key",key, "continue searching");
			}
		}

		return std::make_pair(key,result);
	}

};


int main(int argc, char** argv) {
//    madness::World& universe = madness::initialize(argc,argv);
    initialize(argc, argv);
    World universe(SafeMPI::COMM_WORLD);
    startup(universe,argc,argv);
    FunctionDefaults<4>::set_thresh(1.e-9);
    FunctionDefaults<4>::set_k(7);


    std::cout << "Hello from " << universe.rank() << std::endl;
    universe.gop.fence();
    int nworld=std::min(int(universe.size()),int(3));
    if (universe.rank()==0) print("creating nworld",nworld);


	/**
	 * 	vectorization model
	 *
	 *    	std::vector<data_type<double,3> > vinput(ntask);	// fill with input data
	 *    	macro_task<data_type<double,3> > task;				// implements run(World& world, const data_type& d);
	 *    	macro_taskq<taskT> taskq(universe,nworld);
	 *    	std::vector<Function<double,3> > result=taskq.run_all(task,vinput,fence=true);
	 *
	 */

    long ntask=20;

    // set up input data
    typedef data_type<double,4> dataT;
    std::vector<dataT> vdata;
    for (int i=0; i<ntask; ++i) {
    	Function<double,4> f(universe);
    	f.add_scalar(i);
    	vdata.push_back(dataT(i,i,f));
    }

    // set up taskq
    typedef macro_task<Function<double,4>, dataT> taskT;
    macro_taskq<taskT> taskq(universe,nworld);
    taskT task;
    std::vector<Function<double,4> > result=taskq.run_all(task,vdata);

    for (int i=0; i<result.size(); ++i) {
    	double val=result[i]({0,0,0,0});
    	if (universe.rank()==0) print("result",i,val);
    }


    madness::finalize();
    return 0;
}

template <> volatile std::list<detail::PendingMsg> WorldObject<macro_taskq<macro_task<Function<double,4>, data_type<double, 4ul> > > >::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<macro_taskq<macro_task<Function<double,4>, data_type<double, 4ul> > > >::pending_mutex(0);

template <> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<int, double, madness::Hash<int> > >::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<WorldContainerImpl<int, double, madness::Hash<int> > >::pending_mutex(0);

template <> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, macro_task<Function<double,4>, data_type<double, 4ul> >, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<WorldContainerImpl<long, macro_task<Function<double,4>, data_type<double, 4ul> >, madness::Hash<long> > >::pending_mutex(0);
