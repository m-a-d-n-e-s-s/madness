#include <iostream>
#include <typestuff.h>
#include <tasks/sav.h>
#include <list>
#include <vector>

using namespace std;
using namespace madness;

// more efficient than the current polling loop(s)
// would be a linking of input arguments
// and tasks so that assigning an input argument
// decrements a dependency count in the task,
// and a zero count triggers an automatic move
// to the ready queue.  not hard to do ... but
// since it does not change the public nterface 
// we can do it later.

class TaskInterface {
public:
	virtual bool probe() const = 0;
	virtual void run() = 0;
	virtual ~TaskInterface() {};
};

/// A task with one input argument.
template <typename OpT, typename T>
class Task1in : public TaskInterface {
private:
    OpT op;
    T arg;
public:
	Task1in(OpT op, const T& arg) : op(op), arg(arg) {};
	
	bool probe() const {return arg.probe();};
	
	void run() {op(arg.get());};
};

/// A task with one output argument.
template <typename OpT, typename T>
class Task1out : public TaskInterface {
private:
    OpT op;
    T arg;
public:
	Task1out(OpT op, T& arg) : op(op), arg(arg) {};
	
	bool probe() const {return true;};
	
	void run() {op(arg);};
};

/// A task with one input and one output argument.
template <typename OpT, typename inT, typename outT>
class Task1in1out : public TaskInterface {
private:
    OpT op;
    inT inarg;
    outT outarg;
public:
	Task1in1out(OpT op, inT& inarg, outT& outarg) 
	: op(op), inarg(inarg), outarg(outarg) {};
	
	bool probe() const {return inarg.probe();};
	
	void run() {op(inarg.get(),outarg);};
};

class TaskQueue {
private:
	// Is there a more efficient choice than list?
    list< TaskInterface* > ready;
    list< TaskInterface* > pending;
    
public:
	/// Add a new task ... the task queue takes ownership of the pointer.
	void add(TaskInterface* t) {
		if (t->probe()) ready.push_back(t);
		else pending.push_back(t);
	};
	
	/// Probe pending tasks and move the first ready one to ready queue.
	
	/// Returns true if a ready task was found, false otherwise.
	bool probe() {
		if (!pending.empty()) {
			for (list<TaskInterface *>::iterator p = pending.begin(); 
				 p != pending.end(); ++p) {
				TaskInterface *tp = *p;
			 	if (tp->probe()) {
			 		ready.push_back(tp);
			 		pending.erase(p);
			 		return true;
			 	}
			}
		}
		return false;
	};	 
	
	/// Runs the next ready task if there is one
	void run_next_ready_task() {
		if (!ready.empty()) {
			TaskInterface *p = ready.front();
			p->run();
			ready.pop_front();
		}
	};
	
	// Need a probe all to ensure progress of multistep stuff
	
	// Need hooks to permit use of MPI_Testany instead of busy wait
	
	/// Runs until all tasks have been completed
	void wait() {
		while (!(pending.empty() && ready.empty())) {
			probe();
			run_next_ready_task();
		}
	};
};

void doit(int a) {
	cout << "Doing " << a << endl;
}

void setit(SAV<int>& a) {
	cout << "Setting 99" << endl;
	a.set(99);
};

void chain(int a, SAV<int>& b) {
	cout << "Chaining " << a << endl;
	b.set(a-1);
}

int main() {
	TaskQueue q;
	typedef Task1in< void (*)(int), SAV<int> > doitT;
	typedef Task1out< void (*)(SAV<int>&), SAV<int> > setitT;
	typedef Task1in1out< void (*)(int,SAV<int>&), SAV<int>, SAV<int> > chainT;
	
	// A bunch of tasks all depending on the same assignment
	SAV<int> a;
	q.add(new doitT(&doit, a));
	q.add(new doitT(&doit, a));
	q.add(new doitT(&doit, a));
	q.add(new doitT(&doit, a));
	q.add(new doitT(&doit, a));
	q.add(new doitT(&doit, a));
	q.add(new doitT(&doit, a));
	q.add(new setitT(&setit, a));
	q.wait();
	
	// A chain of tasks triggered by the final assignment
	SAV<int> b,c,d,e,f,g;
	q.add(new chainT(&chain, b,c));
	q.add(new chainT(&chain, c,d));
	q.add(new chainT(&chain, d,e));
	q.add(new chainT(&chain, e,f));
	q.add(new chainT(&chain, f,g));
	q.add(new doitT(&doit, g));
	q.add(new setitT(&setit,b));
	q.wait();
	
	return 0;
}


