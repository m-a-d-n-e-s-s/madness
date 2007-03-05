#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <world.h>

using namespace madness;
using namespace std;

class Node {
private:
    int data;
    bool c[2][2];

    void nochildren() {
	for (int i=0; i<2; i++)
	    for (int j=0; j<2; j++)
		c[i][j] = false;
    };

public:
    Node(int d=0) : data(d) {
	nochildren();
    };

    bool has_child(int i, int j) const {
	return c[i][j];
    };

    void set_child(int i, int j) {
	c[i][j] = true;
    };
    
    void set_data(int d) {
	data = d;
    };

    int get_data() const {
	return data;
    };

    template <typename Archive>
    void serialize(const Archive& ar) {
	ar & data & c;
    }
};

struct Key {
    const unsigned int n;
    const unsigned int i;
    const unsigned int j;
    const hashT hashval;

    Key() 
        : n(0), i(0), j(0), hashval(0) {};

    Key(int n, int i, int j) 
        : n(n), i(i), j(j), hashval(madness::hash(&this->n,3,0)) {};

    template <typename Archive>
    void serialize(const Archive& ar) {
	ar & n & i & j & hashval;
    }

    hashT hash() const {return hashval;};

    bool operator==(const Key& b) const {
        return hashval==b.hashval && n==b.n && i==b.i && j==b.j;
    };
};


template <typename keyT>
class MyProcmap {
private:
    const ProcessID owner;
public:
    MyProcmap(ProcessID owner) : owner(owner) {};

    ProcessID operator()(const keyT& key) const {
	return owner;
    };
};

typedef DistributedContainer< Key,Node,MyProcmap<Key> > treeT;

void build_tree(treeT& tree, const Key& key) {
    int data = key.n*1000000+key.i*10000+key.j;
    Node parent = Node(data);
    if (key.n < 3) {
	for (int p=0; p<2; p++) {
	    for (int q=0; q<2; q++) {
		parent.set_child(p,q);
		build_tree(tree,Key(key.n+1,2*key.i+p,2*key.j+q));
	    }
	}
    }
    tree.insert(key,parent);
}


void print_tree(treeT& tree, const Key& key) {
    treeT::iterator it = tree.find(key);
    if (it!=tree.end()) {
	const Node& node = it->second;
	for (int i=0; i<(int)key.n; i++) cout << "   ";
	print(key.n,key.i,key.j,"owner",tree.owner(key));

	for (int p=0; p<2; p++) {
	    for (int q=0; q<2; q++) {
		if (node.has_child(p,q)) print_tree(tree,Key(key.n+1,2*key.i+p,2*key.j+q));
	    }
	}
    }
}


int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);
    ProcessID me = world.rank();


    try {
	treeT tree(world,MyProcmap<Key>(1));
	if (me == 0) build_tree(tree,Key(0,0,0));
	world.gop.fence();
	if (me == 1) print_tree(tree,Key(0,0,0));
    } catch (MPI::Exception e) {
        error("caught an MPI exception");
    } catch (madness::MadnessException e) {
        print(e);
        error("caught a MADNESS exception");
    } catch (const char* s) {
        print(s);
        error("caught a string exception");
    } catch (...) {
        error("caught unhandled exception");
    }

    world.gop.fence();
    MPI::Finalize();
    return 0;
}
