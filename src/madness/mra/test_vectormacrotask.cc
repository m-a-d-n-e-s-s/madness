/*
 * test_macrotask.cc
 *
 *  Created on: Dec 17, 2019
 *      Author: fbischoff
 */
//#define WORLD_INSTANTIATE_STATIC_TEMPLATES

#include <atomic>
static std::atomic<long> reading_time;	// in ms
static std::atomic<long> writing_time;	// in ms

#include <madness/mra/mra.h>
#include <iostream>
#include <madness/world/MADworld.h>
#include <madness/world/worlddc.h>
#include <random>
#include <madness/mra/funcimpl.h>
#include <archive.h>
#include <madness/world/cloud.h>

using namespace madness;
using namespace archive;


/**
 * @Robert:
 *   - copying Function to/from different worlds without fence in the source world
 *     currently this is done by writing/reading to disk
 *   - serializing tasks not used right now, but will be in the future (submit tasks from within a task)
 *
 * 	Issues:
 *   - default map is OK as long as worlds (universe, subworlds) are disjoint
 *   - turn data structure into tuple
 *   - serialize Function pointer (cast to int64) using archive, serialize Function data using parallel archive
 *   - priority q on rank 0, rank 0 does/nt respond to requests for tasks, does bookkeeping.
 *
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

	int nsubworld=nworld;
	int color = universe.rank() % nsubworld;
	SafeMPI::Intracomm comm = universe.mpi.comm().Split(color, universe.rank() / nsubworld);

	std::shared_ptr<World> all_worlds;
	all_worlds.reset(new World(comm));

	universe.gop.fence();
//	subworld.gop.fence();
	return all_worlds;

//	if (universe.size()<nworld) {
//		print("trying to create ",nworld,"world with",universe.size(), "processes");
//		MADNESS_EXCEPTION("increase number of processes",1);
//	}
//
//	if (universe.rank()==0) print("== multiple worlds created with Intracomm::Create()==",nworld);
//    std::vector<std::vector<int> > process_list(nworld);
//    std::shared_ptr<World> all_worlds;
//
//	for (int i=0; i<universe.size(); ++i) process_list[i%nworld].push_back(i);
//	if (universe.rank()<nworld) print("process_list",process_list[universe.rank()]);
//
//
//	for (int i=0; i<process_list.size(); ++i) {
//		const std::vector<int>& pl=process_list[i];
//		bool found=(std::find(pl.begin(),pl.end(),universe.rank())!=pl.end());
//		if (found) {
//			print("assigning rank",universe.rank(),"to world group",pl);
//
//			SafeMPI::Group group = universe.mpi.comm().Get_group().Incl(pl.size(), &pl[0]);
//			SafeMPI::Intracomm comm_group = universe.mpi.comm().Create(group);
//
//			all_worlds.reset(new World(comm_group));
//		}
//	}
//	universe.gop.fence();
//	return all_worlds;
}


/// base class
class MacroTaskBase {
public:
	MacroTaskBase() {}
	virtual ~MacroTaskBase() {};

	double priority=0.0;
	enum Status {Running, Waiting, Complete, Unknown} stat=Unknown;

	void set_complete() {stat=Complete;}
	void set_running() {stat=Running;}
	void set_waiting() {stat=Waiting;}

	bool is_complete() const {return stat==Complete;}
	bool is_running() const {return stat==Running;}
	bool is_waiting() const {return stat==Waiting;}

	virtual void run(World& world, Cloud& cloud) = 0;
    virtual void print_me(std::string s="") const {}

};


std::ostream& operator<<(std::ostream& os, const MacroTaskBase::Status s) {
	if (s==MacroTaskBase::Status::Running) os << "Running";
	if (s==MacroTaskBase::Status::Waiting) os << "Waiting";
	if (s==MacroTaskBase::Status::Complete) os << "Complete";
	if (s==MacroTaskBase::Status::Unknown) os << "Unknown";

	return os;
}

template<typename macrotaskT>
class MacroTaskIntermediate : public MacroTaskBase {
    
public:

	MacroTaskIntermediate() {}

	~MacroTaskIntermediate() {}

	void run(World& world, Cloud& cloud) {
		dynamic_cast<macrotaskT*>(this)->run(world,cloud);
		world.gop.fence();
	}

    macrotaskT& macrotask() {return dynamic_cast<macrotaskT&>(*this);}

};

static std::atomic<int> atomic_idx;



class MacroTask : public MacroTaskIntermediate<MacroTask> {

public:
	int inputrecord=0;
	int outputrecord=1000;
	int idx=0;
	double exponent=1.0;

	MacroTask(double e) : idx(atomic_idx++), exponent(e) {
	}

    template <typename Archive>
    void serialize(const Archive& ar) {
    	ar & inputrecord & outputrecord & idx;
    }

	void run(World& world, Cloud& cloud) {
		{
	//		print("doing something in world",world.id());
			if (world.rank()==0) print("task",idx,"is reading from record",inputrecord);
			Function<double,4> f=real_factory_4d(world);
			cloud.load(world,f,inputrecord);
			Function<double,4> g=real_factory_4d(world).functor(gaussian(exponent));
			Function<double,4> f2=square(f)+g;
			Derivative<double,4> D(world,1);
			Function<double,4> df2=(D(f2)).truncate();
			Function<double,4> df3=(D(df2)).truncate();
			double trace=df2.trace();
			world.gop.fence();
			if (world.rank()==0) print("task",idx,"is writing to record", outputrecord);
			cloud.store(world,f2,outputrecord);
		}
		world.gop.fence();
	}


    void print_me(std::string s="") const {
    	print("task",s, idx,this,this->stat);
    }

};


template<typename resultT, typename dataT>
class MacroTask1 : public MacroTaskIntermediate<MacroTask1<resultT, dataT> > {

public:
	int idx=1;
	int inputrecord=0;
	int outputrecord=1000;

	MacroTask1() {}

	void run(World& world, Cloud& cloud) {
//		print("doing something in world",world.id());
		Function<double,4> f=cloud.load<Function<double,4> > (world,inputrecord);
		Function<double,4> g=real_factory_4d(world).functor(gaussian(2.0));
		Function<double,4> f2=square(f)+g;
//		Derivative<double,4> D(world,1);
//		Function<double,4> df2=(D(f2)).truncate();
//		double trace=df2.trace();
		world.gop.fence();
		cloud.store(world,g,outputrecord);
	}

    template <typename Archive>
    void serialize(const Archive& ar) {
    	ar & inputrecord & outputrecord & idx;
    }

    void print_me(std::string s="") const {
    	print("task",s, idx,this,this->stat);
    }


};



/// TODO: replicate input data
class macro_taskq : public WorldObject< macro_taskq> {

    World& universe;
    std::shared_ptr<World> subworld_ptr;
	std::vector<std::shared_ptr<MacroTaskBase> > taskq;
	std::mutex taskq_mutex;

public:

	madness::Cloud cloud;
	World& get_subworld() {return *subworld_ptr;}

    /// create an empty taskq and initialize the regional world groups
	macro_taskq(World& universe, int nworld)
		  : universe(universe), WorldObject<macro_taskq>(universe), taskq(), cloud(universe) {

		subworld_ptr=create_worlds(universe,nworld);
		this->process_pending();
	}

	/// run all tasks, leave result in the tasks
	void run_all(std::vector<std::shared_ptr<MacroTaskBase> >& vtask) {

		for (auto t : vtask) if (universe.rank()==0) t->set_waiting();
		for (int i=0; i<vtask.size(); ++i) add_replicated_task(vtask[i]);
		print_taskq();

		universe.gop.fence();
		FunctionDefaults<4>::set_default_pmap(get_subworld());

		double cpu00=cpu_time();

		World& subworld=get_subworld();
		while (true){
			long element=get_scheduled_task_number(subworld);
			if (element<0) break;
			double cpu0=cpu_time();
			std::shared_ptr<MacroTaskBase> task=taskq[element];
                        
			task->run(subworld,cloud);

			double cpu1=cpu_time();
			set_complete(element);
			if (subworld.rank()==0) printf("completed task %3ld after %4.1fs\n",element,cpu1-cpu0);

		};
		universe.gop.fence();
		FunctionDefaults<4>::set_default_pmap(universe);
		double cpu11=cpu_time();
		printf("completed taskqueue after %4.1fs\n",cpu11-cpu00);


	}

//	/// run the task on the vector of input data, return vector of results
//	template<typename taskT>
//	std::vector<typename taskT::result_type> map(taskT& task1,
//			std::vector<typename taskT::data_type>& vdata) {
//
//		// create copies of the input task instance and fill with the data
//		std::vector<std::shared_ptr<MacroTaskBase> > vtask(vdata.size());
//		for (int i=0; i<vdata.size(); ++i) {
//			vtask[i]=std::shared_ptr<MacroTaskBase>(new taskT(task1));
//			dynamic_cast<taskT&>(*vtask[i].get()).set_data(vdata[i]);
//		}
//
//		// execute the task list
//		run_all(vtask);
//
//		// localize the result into universe
//		std::vector<typename taskT::result_type> vresult(vdata.size());
//		for (int i=0; i<vresult.size(); ++i) {
//			vtask[i]->load_result(universe,"result_of_task"+std::to_string(i));
//			vresult[i]=dynamic_cast<taskT&>(*(vtask[i].get())).get_result();
//		}
//		return vresult;
//	}

private:
	void add_replicated_task(const std::shared_ptr<MacroTaskBase>& task) {
		taskq.push_back(task);
	}

	void print_taskq() const {
		universe.gop.fence();
		if (universe.rank()==0) {
			print("taskq on universe rank",universe.rank());
			for (auto t : taskq) t->print_me();
		}
		universe.gop.fence();
	}

	/// scheduler is located on universe.rank==0
	long get_scheduled_task_number(World& subworld) {
		long number=0;
		if (subworld.rank()==0) number=this->task(ProcessID(0), &macro_taskq::get_scheduled_task_number_local);
		subworld.gop.broadcast_serializable(number, 0);
		subworld.gop.fence();
		return number;

	}

	long get_scheduled_task_number_local() {
		MADNESS_ASSERT(universe.rank()==0);
		std::lock_guard<std::mutex> lock(taskq_mutex);

		auto is_Waiting = [](std::shared_ptr<MacroTaskBase> mtb_ptr) {return mtb_ptr->is_waiting();};
		auto it=std::find_if(taskq.begin(),taskq.end(),is_Waiting);
		if (it!=taskq.end()) {
			it->get()->set_running();
			long element=it-taskq.begin();
			return element;
		}
//		print("could not find task to schedule");
		return -1;
	}

	/// scheduler is located on rank==0
	void set_complete(const long task_number) const {
		this->task(ProcessID(0), &macro_taskq::set_complete_local, task_number);
	}

	/// scheduler is located on rank==0
	void set_complete_local(const long task_number) const {
		MADNESS_ASSERT(universe.rank()==0);
		taskq[task_number]->set_complete();
	}

	std::size_t size() const {
		return taskq.size();
	}

//	void add_task(const std::shared_ptr<MacroTaskBase>& task) {
//		ProcessID master=0;
//		task->print_me("in add_task, universe.rank()="+std::to_string(universe.rank()));
//		MacroTaskBase* taskptr=task.get();
//		thistype::send(master,&thistype::add_task_local,taskptr);
//	};
//
//	void add_task_local(const basetaskptr& task) {
//		MADNESS_ASSERT(universe.rank()==0);
//		std::shared_ptr<MacroTaskBase> task1;
//		task1.reset(task);
//		task1->print_me("in add_task_local");
//		taskq.push(task1);
//	};
//
//	std::shared_ptr<MacroTaskBase> get_task_from_tasklist(World& regional) {
//
//		// only root may pop from the task list
//		std::vector<unsigned char> buffer;
//		if (regional.rank()==0) buffer=pop();
//		regional.gop.broadcast_serializable(buffer, 0);
//		regional.gop.fence();
//
//		std::shared_ptr<MacroTaskBase> task;
//		MacroTaskBase* task_ptr;
//		BufferInputArchive ar(&buffer[0],buffer.size());
//		ar & task_ptr;
//
//		task.reset(task_ptr);
//		return task;
//	}
//
//	/// pass serialized task from universe.rank()==0 to world.rank()==0
//	std::vector<unsigned char> pop() {
//		return this->task(ProcessID(0), &macro_taskq<taskT>::pop_local);
//	}
//
//	/// pop highest-priority task and return it as serialized buffer
//	std::vector<unsigned char> pop_local() {
//		const std::lock_guard<std::mutex> lock(taskq_mutex);
//		std::shared_ptr<MacroTaskBase> task(NULL);
//
//		if (not taskq.empty()) {
//			task=taskq.top();
//			taskq.pop();
//		}
//
//		BufferOutputArchive ar_c;
//		ar_c & task.get();
//		long nbyte=ar_c.size();
//		std::vector<unsigned char> buffer(nbyte);
//
//		BufferOutputArchive ar2(&buffer[0],buffer.size());
//		ar2 & task.get();
//
//		return buffer;
//	}

};


int main(int argc, char** argv) {
//    madness::World& universe = madness::initialize(argc,argv);
    madness::World& universe = madness::initialize(argc, argv);
    startup(universe,argc,argv);
    FunctionDefaults<4>::set_thresh(1.e-9);
    FunctionDefaults<4>::set_k(7);

    std::cout << "Hello from " << universe.rank() << std::endl;
    universe.gop.fence();
//    int nworld=std::min(int(universe.size()),int(3));
//    int nworld=universe.size();
    int nworld=1;
    if (universe.rank()==0) print("creating nworld",nworld);
    writing_time=0;
    reading_time=0;

    {
		macro_taskq taskq(universe,nworld);

		long ntask=35;

		Function<double,4> f=real_factory_4d(universe);
		f=f+1.0;
		taskq.cloud.store(universe,f,0);
		universe.gop.fence();

		std::vector<std::shared_ptr<MacroTaskBase> > vtask;
		for (int i=0; i<ntask; ++i) {
			MacroTask task(double(ntask-i)/2.0);
			task.inputrecord=0;
			task.outputrecord=i+1;
			vtask.push_back(std::shared_ptr<MacroTaskBase>(new MacroTask(task)));
		}
		universe.gop.fence();

		// try heterogeneous task list
	//    vtask.push_back(std::shared_ptr<MacroTaskBase>(new MacroTask1<real_function_4d,real_function_4d>(ff)));
		// set up taskq with a vector of tasks
		taskq.run_all(vtask);

		if (universe.rank()==0) {
			print("cloud storing time", writing_time*0.001);
			print("cloud reading time", reading_time*0.001);
		}
	//    MacroTask<real_function_4d, dataT> task;
	//    task.idx=3;
	//    std::vector<Function<double,4> > result=taskq.map(task,vdata);

	//    print_size(universe,result,"result after map");
    }
	universe.gop.fence();

    madness::finalize();
    return 0;
}

template <> volatile std::list<detail::PendingMsg> WorldObject<macro_taskq>::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<macro_taskq>::pending_mutex(0);

template <> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, MacroTask, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<WorldContainerImpl<long, MacroTask, madness::Hash<long> > >::pending_mutex(0);

template <> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending_mutex(0);
