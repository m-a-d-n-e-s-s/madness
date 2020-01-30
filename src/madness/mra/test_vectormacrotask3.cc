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
#include <archive.h>

using namespace madness;
using namespace archive;

/**
 * 	Issues:
 * 	 - set_defaults<>(local_world) for loading (affects pmap)
 * 	 - serialization of task works only for int, double, .. but not for Function -> separate task from data
 * 	 - save/load of task data: must save data upon creation instead of consumption, because serialization of Function fails
 * 	 - macro_taskq as WorldObject??
 *   - turn data structure into tuple
 *   - prioritize tasks
 *   - submit tasks from within other tasks -> how to manage results?
 *
 */


/**
 *  - default map is OK as long as worlds (universe, subworlds) are disjoint
 *  - serialize Function pointer (cast to int64) using archive, serialize Function data using parallel archive
 *	- priority q on rank 0, rank 0 does/nt respond to requests for tasks, does bookkeeping.
 *
 *	- world.h load/store base pointer
 */


///
struct gaussian {
	double a;
	gaussian() : a() {};
	gaussian(double aa) : a(aa) {}
	double operator()(const coord_4d& r) const {
		double x=r[0], y=r[1], z=r[2], aa=r[3];
		return exp(-a*(x*x + y*y + z*z * aa*aa));//*abs(sin(abs(2.0*x))) *cos(y);
	}
};

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
        ar & i & d & fexist;

        if (fexist) {
			if (ar.is_output_archive) {
				ar & f.get_impl();
			}
			if (ar.is_input_archive) {
				FunctionImpl<T,NDIM>* impl;
    			ar & impl;
    			std::shared_ptr<FunctionImpl<T,NDIM> > fimpl(impl);
    			f.set_impl(fimpl);
    		}
    	}
    }

//    void localize(World& origin, World& destination) {
//    	f=::localize(origin,destination,f,i,"dummy");
//    }


	void save(World& world) const {
		world.gop.fence();
		ParallelOutputArchive ar(world, filename.c_str() , 1);
//		print("saving to file",filename);
		ar & d & i & f;
	}

	void load(World& world) {
		world.gop.fence();
        auto pmap = std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >(new madness::LevelPmap< Key<NDIM> >(world));
        FunctionDefaults<NDIM>::set_pmap(pmap);	// set default pmap to use only this world!
//		print("loading from file",filename, world.id());
		ParallelInputArchive ar(world, filename.c_str() , 1);
		ar & d & i & f;
	}

};



class MacroTaskBase;

typedef MacroTaskBase* (*inputfuntype)(const BufferInputArchive& ar);
typedef MacroTaskBase* (*outputfuntype)(BufferOutputArchive& ar);

/// base class
//template<typename resultT, typename dataT>
class MacroTaskBase {
public:
//	virtual resultT run(World& world, const dataT& data) = 0;
	MacroTaskBase() {}
	virtual ~MacroTaskBase() {};

	double priority=0.0;

	virtual void run(World& world) = 0;
    virtual inputfuntype get_allocate_and_deserialize_method() = 0;

    virtual void print_me(std::string s="") const {}
    virtual void store(const BufferOutputArchive& ar) = 0;

};


template<typename macrotaskT>
class MacroTaskIntermediate : public MacroTaskBase {
public:

	MacroTaskIntermediate() {}

	~MacroTaskIntermediate() {}

	void run(World& world) {
		dynamic_cast<macrotaskT*>(this)->run(world);
	}

    void store(const BufferOutputArchive& ar) {
    	dynamic_cast<macrotaskT*>(this)->serialize(ar);
    }

    static MacroTaskBase* allocate_and_deserialize(const BufferInputArchive& ar) {
    	macrotaskT* t = new macrotaskT;
        t->serialize(ar);
        return t;
    }

    inputfuntype get_allocate_and_deserialize_method() {
        return &allocate_and_deserialize;
    }

//
//	void save(World& world, const dataT& data) const {
//		data.save(world);
//	}
//
//	void load(World& world, dataT& data) {
//		data.load(world);
//	}

};

namespace madness {
namespace archive {

    template <typename, typename>
    struct ArchiveLoadImpl;
    template <typename, typename>
    struct ArchiveStoreImpl;

    /// Specialization of \c ArchiveLoadImpl for \c World pointers.

    /// Helps in archiving (reading) \c World objects.
    /// \tparam Archive The archive type.
    template <class Archive>
    struct ArchiveLoadImpl<Archive,MacroTaskBase* > {
        /// Loads a \c World from the specified archive.

        /// \param[in,out] ar The archive.
        static inline void load(const Archive& ar, MacroTaskBase*& mtb_ptr) {
        	mtb_ptr=NULL;
        	bool exist;
        	ar & exist;
        	if (exist) {
				inputfuntype voodoo;
				ar & voodoo;
				mtb_ptr=voodoo(ar);
				MADNESS_ASSERT(mtb_ptr);
        	}
        }
    }; // struct ArchiveLoadImpl<Archive,World*>

    /// Specialization of \c ArchiveStoreImpl for \c World pointers.

    /// Helps in archiving (writing) \c World objects.
    /// \tparam Archive The archive type.
    template <class Archive>
    struct ArchiveStoreImpl<Archive,MacroTaskBase*> {
        /// Writes a \c World to the specified archive.

        /// \param[in,out] ar The archive.
        static inline void store(const Archive& ar, MacroTaskBase* const & mtb_ptr) {
        	bool exist=(mtb_ptr!=NULL);
        	ar & exist;
        	if (exist) {
				auto voodoo=mtb_ptr->get_allocate_and_deserialize_method();
				ar & voodoo;
				mtb_ptr->store(ar);
        	}
        }
    }; // struct ArchiveStoreImpl<Archive,World*>
} // namespace archive
} // namespace madness





template<typename resultT, typename dataT>
class MacroTask : public MacroTaskIntermediate<MacroTask<resultT, dataT> > {

public:
	typedef resultT result_type;
	typedef dataT data_type;
	int stuff;
	dataT data;

	MacroTask() : stuff(1) {}
	MacroTask(const dataT& data) : stuff(2), data(data) {}

	void run(World& world) {
		const Function<double,4>& f=data.f;
		Function<double,4> g=real_factory_4d(world).functor(gaussian(data.d));
		Function<double,4> f2=square(f)+g;
		Derivative<double,4> D(world,1);
		Function<double,4> df2=(D(f2)).truncate();
		double trace=df2.trace();
	}

    template <typename Archive>
    void serialize(const Archive& ar) {

    	ar & stuff & data;
    }

    void print_me(std::string s="") const {
    	print("task",s,stuff, data.i,this);
    }


//
//	resultT run(World& world, const dataT& data) {
//		const Function<double,4>& f=data.f;
//		Function<double,4> g=real_factory_4d(world).functor(gaussian(data.d));
//		Function<double,4> f2=square(f)+g;
//		Derivative<double,4> D(world,1);
//		Function<double,4> df2=(D(f2)).truncate();
//		double trace=df2.trace();
//		return df2;
//	}
};


template<typename T, std::size_t NDIM>
Function<T,NDIM> localize(World& origin, World& destination, const Function<T,NDIM>& data,
		const long id, std::string filename="smartie") {

	origin.gop.fence();
	destination.gop.fence();
	filename+=std::to_string(id);
	if (data.is_initialized() and data.world().id()==origin.id()) save(data,filename);

	destination.gop.fence();
	origin.gop.fence();

    auto pmap = std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >(new madness::LevelPmap< Key<NDIM> >(destination));
    FunctionDefaults<NDIM>::set_pmap(pmap);	// set default pmap to use only this world!
	Function<T,NDIM> result(destination);
	load(result,filename);

	destination.gop.fence();
	origin.gop.fence();

    ParallelInputArchive ar2(destination, filename.c_str(), 1);
    ar2.remove();

	return result;
}

class MasterPmap : public WorldDCPmapInterface<long> {
public:
    MasterPmap() {}
    ProcessID owner(const long& key) const {return 0;}
};



template<typename taskT>
class macro_taskq : public WorldObject< macro_taskq<taskT> > {
    typedef macro_taskq<taskT> thistype;
    typedef typename taskT::result_type resultT;
    typedef typename taskT::data_type dataT;
    typedef MacroTaskBase* basetaskptr;

    World& universe;
    std::shared_ptr<World> subworld_ptr;
	std::priority_queue<std::shared_ptr<MacroTaskBase> > taskq;
	std::mutex taskq_mutex;


public:

	World& get_subworld() {return *subworld_ptr;}

    /// create an empty taskq and initialize the regional world groups
	macro_taskq(World& universe, int nworld)
		  : universe(universe), WorldObject<thistype>(universe),
			taskq() {

		subworld_ptr=(create_worlds(universe,nworld));
	    World& subworld=*(subworld_ptr.get());
		this->process_pending();
	}

	std::vector<resultT> run_all(std::vector<std::shared_ptr<MacroTaskBase> >& vtask) {

		std::vector<resultT> result(vtask.size());
		for (int i=0; i<vtask.size(); ++i) {
			add_task(vtask[i]);
		}
		myusleep(1*1.e6);
		print("taskq size",taskq.size());

		universe.gop.fence();
		World& subworld=get_subworld();
		print("\n\n\n");

		double total_cpu=0.0;
		double cp0=cpu_time();
		bool working=true;
		int counter=0;
		while (working) {
			std::shared_ptr<MacroTaskBase> task=get_task_from_tasklist(subworld);
			if (not task) break;
			print("taskq size",taskq.size());
			double cpu0=cpu_time();
			task->run(subworld);
			double cpu1=cpu_time();
			printf("finished task in %4.2fs\n",cpu1-cpu0);
			total_cpu+=(cpu1-cpu0);
		}
		universe.gop.fence();
		subworld.gop.fence();
		double cp1=cpu_time();
		double fence2fence=cp1-cp0;
		universe.gop.sum(total_cpu);
		if (universe.rank()==0) printf("wall time spent in task execution %2.2fs\n",fence2fence);
		if (universe.rank()==0) printf("CPU time spent in task execution  %2.2fs\n",total_cpu);

		double cpu0=cpu_time();
	    for (int i=0; i<result.size(); ++i) {
//	    	result[i]=localize(subworld,universe,result[i],i);
	    }
		double cpu1=cpu_time();
		double total_localize=cpu1-cpu0;
		universe.gop.sum(total_localize);

		if (universe.rank()==0) printf("time spent in localization %2.2fs\n",total_localize);

		return result;

	}

	void add_task(const std::shared_ptr<MacroTaskBase>& task) {
		ProcessID master=0;
		task->print_me("in add_task, universe.rank()="+std::to_string(universe.rank()));
		MacroTaskBase* taskptr=task.get();
		thistype::send(master,&thistype::add_task_local,taskptr);
	};

	void add_task_local(const basetaskptr& task) {
		MADNESS_ASSERT(universe.rank()==0);
		std::shared_ptr<MacroTaskBase> task1;
		task1.reset(task);
		task1->print_me("in add_task_local");
		taskq.push(task1);
	};

	std::shared_ptr<MacroTaskBase> get_task_from_tasklist(World& regional) {

		// only root may pop from the task list
		std::vector<unsigned char> buffer;
		if (regional.rank()==0) buffer=pop();
		regional.gop.broadcast_serializable(buffer, 0);
		regional.gop.fence();

		std::shared_ptr<MacroTaskBase> task;
		MacroTaskBase* task_ptr;
		BufferInputArchive ar(&buffer[0],buffer.size());
		ar & task_ptr;

		task.reset(task_ptr);
		return task;
	}

	std::size_t size() const {
		return taskq.size();
	}

	/// pass serialized task from universe.rank()==0 to world.rank()==0
	std::vector<unsigned char> pop() {
		return this->task(ProcessID(0), &macro_taskq<taskT>::pop_local);
	}

	/// pop highest-priority task and return it as serialized buffer
	std::vector<unsigned char> pop_local() {
		const std::lock_guard<std::mutex> lock(taskq_mutex);
		std::shared_ptr<MacroTaskBase> task(NULL);

		if (not taskq.empty()) {
			task=taskq.top();
			taskq.pop();
		}

		BufferOutputArchive ar_c;
		ar_c & task.get();
		long nbyte=ar_c.size();
		std::vector<unsigned char> buffer(nbyte);

		BufferOutputArchive ar2(&buffer[0],buffer.size());
		ar2 & task.get();

		return buffer;
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
    typedef MacroTask<real_function_4d, dataT> taskT;

    std::vector<std::shared_ptr<MacroTaskBase> > vtask;
    for (int i=0; i<ntask; ++i) {
    	Function<double,4> f(universe);
    	f.add_scalar(i);
    	vtask.push_back(std::shared_ptr<MacroTaskBase>(new taskT(dataT(i,i,f))));
    }

//    // set up taskq
    macro_taskq<taskT> taskq(universe,nworld);
    double cpu0=cpu_time();
    double wall0=wall_time();

//    std::vector<Function<double,4> > result=taskq.run_all(vtask);

    MacroTaskBase* task=vtask.begin()->get();
    task->print_me("initial task");
//    if (universe.rank()==1) task=NULL;
//    if (not task) print("nullified task on rank",universe.rank());

	BufferOutputArchive ar_c;
	ar_c & task;
	long nbyte=ar_c.size();
	char* buf=new char[nbyte];
	print("nbyte",nbyte);

    {

    	BufferOutputArchive ar(buf,nbyte);
    	ar & task;
        task->print_me("saving into buffer");

    }
	universe.gop.fence();
    {
    	task=NULL;
    	BufferInputArchive ar(buf,nbyte);
    	ar & task;
        task->print_me("loading from buffer");

    }

    print("\n\n\n");
//	universe.gop.broadcast_serializable(dynamic_cast<taskT&>(*task), 0);
	universe.gop.broadcast_serializable(task, 0);
	universe.gop.fence();
	if (task) task->print_me("broadcasted task, universe.rank()="+std::to_string(universe.rank()));




    double cpu1=cpu_time();
    double wall1=wall_time();

    if (universe.rank()==0) printf("done with macrotasks in %.2fs (cpu) %.2fs (wall)\n",cpu1-cpu0,wall1-wall0);
//    for (int i=0; i<result.size(); ++i) {
//    	double val=result[i]({0,0,0,0});
//    	if (universe.rank()==0) print("result",i,val);
//    }


    madness::finalize();
    return 0;
}

template <> volatile std::list<detail::PendingMsg> WorldObject<macro_taskq<MacroTask<Function<double,4>, data_type<double, 4ul> > > >::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<macro_taskq<MacroTask<Function<double,4>, data_type<double, 4ul> > > >::pending_mutex(0);

template <> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, MacroTask<Function<double,4>, data_type<double, 4ul> >, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<WorldContainerImpl<long, MacroTask<Function<double,4>, data_type<double, 4ul> >, madness::Hash<long> > >::pending_mutex(0);
