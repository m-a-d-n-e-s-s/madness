#include <iostream>
#include <tasks/tasks.h>
#include <tasks/sav.h>
#include <misc/communicator.h>

using namespace std;
using namespace madness;

Communicator* comm_default;

template <typename T>
class SAVMPIImpl : public SAVImpl<T> {
private:
	int rank;
	int tag;
	bool input;
	mutable MPI::Request handle;

public:
    SAVMPIImpl() : SAVImpl<T>(), rank(-1), tag(-1), input(true) {};

    SAVMPIImpl(const T& t) : SAVImpl<T>(t), rank(-1), tag(-1), input(true) {};

	SAVMPIImpl(int rank, int tag, bool input) 
		: SAVImpl<T>(), rank(rank), tag(tag), input(input) {
		if (input) handle = comm_default->Irecv(*(this->get_ptr()), rank, tag);
	};
	
    void assign_action() const {
    	if (this->rank>=0 && !this->input) 
    		comm_default->Send(this->get(), rank, tag);
    };

    void probe_action() const {
    	if (this->rank>=0 && this->input && this->handle.Test()) 
    		this->set_assigned();
    };

    virtual ~SAVMPIImpl() {};
};

 
template <class T>
class SAVMPI : public SAVInterface<T> {
public:
	typedef madness::SharedPtr< SAVMPIImpl<T> > ptrT;

private:
    ptrT p;

public:
	/// Constructs a local unassigned variable which will be assigned locally
    SAVMPI() : p(new SAVMPIImpl<T>()) {};
    
    /// Constructs a local assigned variable
    SAVMPI(const T& t) : p(new SAVMPIImpl<T>(t)) {};
    
    /// Constructs an unassigned variable associated with message passing
    
    /// If input is true, it is an input variable that will be assigned via a message 
    /// from rank with given tag.  Otherwise, it is an output variable that when
    /// assigned to will send a message to rank with given tag.
    SAVMPI(int rank, int tag, bool input) : p(new SAVMPIImpl<T>(rank,tag,input)) {};
	
    inline void set(const T& value) {
        p->set(value);
    };

    inline const T& get() const {
        return p->get();
    };

    inline bool probe() const {
        return p->probe();
    };
    
    virtual ~SAVMPI() {};
};


void doit(int a) {
	cout << comm_default->rank() << " Doing " << a << endl;
}

void setit(SAVMPI<int>& a) {
	cout << comm_default->rank() << " Setting 99" << endl;
	a.set(99);
}

void chain(int a, SAVMPI<int>& b) {
	cout << comm_default->rank() << " Chaining " << a << endl;
	b.set(a-1);
}

template <typename T>
inline
SAVMPI<T> input() {return SAVMPI<T>();};

template <typename T>
inline
SAVMPI<T> input(int rank, int tag) {return SAVMPI<T>(rank, tag, true);};

template <typename T>
inline
SAVMPI<T> input(const T& t) {return SAVMPI<T>(t);};

template <typename T>
inline
SAVMPI<T> output(int rank, int tag) {return SAVMPI<T>(rank, tag, false);};


int main(int argc, char* argv[]) {
    MPI::Init(argc, argv);
    Communicator comm;
    comm_default = &comm;	
    TaskQueue q;
    
	typedef Task1in< void (*)(int), SAVMPI<int> > doitT;
	typedef Task1out< void (*)(SAVMPI<int>&), SAVMPI<int> > setitT;
	typedef Task1in1out< void (*)(int,SAVMPI<int>&), SAVMPI<int>, SAVMPI<int> > chainT;
    
    cout << "Hello from " << comm.rank() << endl;
    
    	
	// A bunch of tasks on node 0 all depending on the same 
	// assignment at node 1;
	
	if (comm.rank() == 0) {
		SAVMPI<int> a(1,1,true);
		q.add(new doitT(&doit, a));
		q.add(new doitT(&doit, a));
		q.add(new doitT(&doit, a));
		q.add(new doitT(&doit, a));
		q.add(new doitT(&doit, a));
		q.add(new doitT(&doit, a));
		q.add(new doitT(&doit, a));
		q.wait();
	}
	else if(comm.rank() == 1) {
		SAVMPI<int> a(0,1,false);
		a.set(33);
	}
	
	// A chain of tasks interleaved between 0 and 1
	if (comm.rank() == 0) {
		SAVMPI<int> a(1,1,true),b(1,2,false);
		SAVMPI<int> c(1,3,true),d(1,4,false);
		SAVMPI<int> e(1,5,true),f(1,6,false);
		q.add(new chainT(&chain, a, b));
		q.add(new chainT(&chain, c, d));
		q.add(new chainT(&chain, e, f));
	}
	else {
		SAVMPI<int> a(0,1,false),b(0,2,true);
		SAVMPI<int> c(0,3,false),d(0,4,true);
		SAVMPI<int> e(0,5,false),f(0,6,true);
		q.add(new chainT(&chain, b, c));
		q.add(new chainT(&chain, d, e));
		q.add(new doitT(&doit, f));
		a.set(10);
	}
	q.wait();
	    
    comm.close(); 
    MPI::Finalize();
	return 0;
}


