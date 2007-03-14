#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <world/world.h>

using namespace madness;
using namespace std;

template <typename Data, unsigned int D>
class Node {
private:
    Data data;
    std::vector<bool> c;

    void nochildren() {
        c.clear();
        c.assign(dim, false);
    };

public:
    static unsigned int dim;

    Node() {
	data = Data();
	nochildren();
    }
    Node(Data d) : data(d) {
	nochildren();
    };

    unsigned int compute_index(vector<unsigned int> v) {
	unsigned int vlen = v.size();
	unsigned int index = 0, twoD = 1;
	for (unsigned int i = vlen-1; i >= 0; i--) {
	    index+= (twoD*v[i]);
	    twoD*=2;
	}
	return index;
    };
	
    unsigned int compute_index(vector<int> v) {
	unsigned int vlen = v.size();
	unsigned int index = 0, twoD = 1;
	for (unsigned int i = vlen-1; i >= 0; i--) {
	    index+= (twoD*v[i]);
	    twoD*=2;
	}
	return index;
    };
	
    bool has_children() const {
	for (unsigned int i = 0; i < dim; i++)
	    if (c[i]) return true;
	return false;
    };

    bool has_child(unsigned int i) const {
	return c[i];
    };

    bool has_child(int i) const {
	return c[i];
    };

    bool has_child(vector<int> vi) const {
	return c[compute_index(vi)];
    };

    void set_child(int i, bool setto = true) {
	c[i] = setto;
    };
    
    void set_child(vector<int> vi, bool setto = true) {
	c[compute_index(vi)] = setto;
    };
    
    void set_data(Data d) {
	data = d;
    };

    Data get_data() const {
	return data;
    };

    template <typename Archive>
    void serialize(const Archive& ar) {
	ar & data & c;
    }
};


template <typename Data, unsigned int D>
unsigned int Node<Data,D>::dim = (unsigned int) pow(2.0, (int) D);


class NodeData {
public:
    int cost;
    int subcost;
    bool istaken;
    NodeData(int c = 1, int s = 1, bool i = false) : cost(c), subcost(s), istaken(i) {};
    template <typename Archive>
    void serialize(const Archive& ar) {
	ar & cost & subcost & istaken;
    };
    void print() {
	cout << "cost = " << cost << ", subcost = " << subcost << ", istaken = " << istaken << endl;
    }
};

template <unsigned int D>
class Key {
public:
//    const unsigned int n;
    unsigned int n;
    mutable vector<unsigned int> L;
//    const std::size_t hashval; // Saving this here is a big optimization
    hashT hashval; // Saving this here is a big optimization

    Key() {};

    Key(const Key& k) : n(k.n) {
	for (unsigned int i = 0; i < D; i++) {
	    L.push_back(k.L[i]);
	}
       hashval = k.hashval;
    };

    Key(int n, int i, int j) : n(n) {
      L.push_back(i); L.push_back(j);
       hashval = madness::hash(&L[0],D,madness::hash(n));
   }

    Key(int n, vector<int> v) : n(n) {
	for (int i = 0; i < v.size(); i++) {
	    L.push_back(v[i]);
	}
	for (int i = v.size(); i < D; i++) {
	    L.push_back(0);
	}
       hashval = madness::hash(&L[0],D,madness::hash(n));
    };

    Key(unsigned int n, vector<unsigned int> v) : n(n) {
	for (unsigned int i = 0; i < v.size(); i++) {
	    L.push_back(v[i]);
	}
	for (unsigned int i = v.size(); i < D; i++) {
	    L.push_back(0);
	}
       hashval = madness::hash(&L[0],D,madness::hash(n));
    };

    hashT hash() const {
	return hashval;
    };


    Key myChild(int k) const {
	vector<unsigned int> LL;
	for (unsigned int i = 0; i < D; i++) {
	    LL.push_back(2*L[i] + k%2);
	    k/=2;
	}
	return Key(n+1, LL);
    };

    Key myParent(int k=1) const {
	if (k == 0) return Key(*this);
	vector<unsigned int> LL;
	unsigned int twotok = (unsigned int) pow(2.0, k);
	for (unsigned int i = 0; i < D; i++) {
	    LL.push_back(L[i]/twotok);
	}
	return Key(n-k, LL);
    };


    bool isChildOf(const Key& key) const {
	if (*this == key) {
	    return true; // I am child of myself
	}
	else if (this->n <= key.n) {
	    return false; // I can't be child of something lower on the tree
	}
	else {
	    unsigned int dn = this->n - key.n;
	    Key mama = this->myParent(dn);
	    return (mama == key);
	}
    };

    bool isParentOf(const Key& key) const {
	return key.isChildOf(*this);
    };
	    

    int ordering(const Key& k1, const Key& k2) const {
	bool egalite = true;
	vector<int> dL(D, 0);
//	cout << "ordering: comparing ";
//	k1.print();
//	cout << " and ";
//	k2.print();
//	cout << endl;
	for (unsigned int i = 0; i < D; i++) {
	    dL[i] = k1.L[i] - k2.L[i];
	    if (k1.L[i]/2 != k2.L[i]/2) {
		egalite = false;
//		cout << "ordering: k1 and k2 do not have same parent" << endl;
	    }
	}
	if (!egalite) {
	    return (ordering(k1.myParent(), k2.myParent()));
	}
	else {
	    for (unsigned int i = 0; i < D; i++) {
		if (dL[i] > 0) {
//		    cout << "ordering: dL[" << i << "] > 0; returning -1" << endl;
		    return -1;
		}
		else if (dL[i] < 0) {
//		    cout << "ordering: dL[" << i << "] < 0; returning 1" << endl;
		    return 1;
		}
	    }
//	    cout << "ordering: no dL greater or less than zero; returning 0" << endl;
	    return 0;
	}
    };

    bool operator==(const Key& a) const {
	if (n != a.n) return false;
	for (unsigned int i = 0; i < D; i++) {
	    if (L[i] != a.L[i]) return false;
	}
	return true;
    };

    bool operator<(const Key& a) const {
	int ans;

	if (n == a.n) {
	    ans = ordering(*this, a);
	}
	else if (n > a.n) {
	    unsigned int dn = n - a.n;
	    Key newthis = this->myParent(dn);
	    if (newthis == a) {
		ans = 1;
	    }
	    else {
	    	ans = ordering(newthis, a);
	    }
	}
	else {
	    unsigned int dn = a.n - n;
	    Key newa = a.myParent(dn);
	    if (newa == *this) {
		ans = 0;
	    }
	    else {
	    	ans = ordering(*this, newa);
	    	if (ans < 0)
		    ans = 1;
	    	else
		    ans = 0;
	    }
	}
	if (ans > 0)
	    return true;
	else
	    return false;
    };

/*
    bool operator<(const Key& a) const {
        if (n < a.n)  return true;
        else if (n == a.n) {
	    for (unsigned int i = 0; i < D; i++) {
                if (L[i] < a.L[i]) return true;
                else if (L[i] > a.L[i]) return false;
            }
            return false;
        }
        else return false;
    };
*/

    void print() const {
	std::cout << "n = " << n << ", (";
	for (unsigned int i = 0; i < D-1; i++) {
	    std::cout << L[i] << ", ";
	}
	std::cout << L[D-1] << ")";
    }

    template <typename Archive>
    void serialize(const Archive& ar) {
	ar & n & L & hashval;
    }
};

typedef Key<2> KeyD;

// To enable gnu hash_map to work
namespace std {
    template <> struct equal_to<KeyD> {
        bool operator()(const KeyD& a, const KeyD& b) const {
	    bool retval = (a.n==b.n);
	    int D = a.L.size();
	    for (int i = 0; i < D; i++) {
		if (!retval) break;
		retval = (retval)&&(a.L[i]==b.L[i]);
	    }
            return retval;
        };
    };
}

// To enable gnu hash_map to work

namespace __gnu_cxx {
    template <> struct hash<KeyD> { 
        std::size_t operator()(const KeyD& key) const {
            return key.hashval;
        };
    };
}

struct TreeCoords {
    KeyD key;
    ProcessID owner;

    TreeCoords(const KeyD k, ProcessID o) : key(KeyD(k)), owner(o) {};
    TreeCoords(const TreeCoords& t) : key(KeyD(t.key)), owner(t.owner) {};
    TreeCoords() : key(KeyD()), owner(-1) {};
    void print() const {
	key.print();
	std::cout << "    owner = " << owner << std::endl;
    };

    bool operator< (const TreeCoords t) const {
	return (this->key < t.key);
    };
};


struct Tree {
    TreeCoords data;
    vector<Tree*> children;
    Tree* parent;

    Tree() {};
    Tree(TreeCoords d) : data(d), parent(0) {};
    Tree(TreeCoords d, Tree* p) : data(d), parent(p) {};

    void insertChild(TreeCoords d) {
	Tree* c = new Tree(d, this);
	children.push_back(c);
    };

    void print() {
	data.print();
	int csize = children.size();
	for (int j = 0; j < csize; j++) {
	    children[j]->print();
	}
    };

    bool isForeparentOf(KeyD key) const {
	return (this->data.key.isParentOf(key));
    };

    void findOwner(const KeyD key, ProcessID *ow) const {
	if (this->isForeparentOf(key)) {
	    *ow = this->data.owner;
	    int csize = children.size();
	    for (int j = 0; j < csize; j++) {
		children[j]->findOwner(key, ow);
	    }
	}
    };

    bool fill(vector<TreeCoords> v, int j) {
	bool success = false;
	if (this->isForeparentOf(v[j].key)) {
	    int csize = children.size();
	    for (int i = 0; i < csize; i++) {
		if (children[i]->isForeparentOf(v[j].key)) {
		    success = children[i]->fill(v, j);
		    
		}
	    }
	    if (!success) {
		this->insertChild(v[j]);
		success = true;
/*
cout << "successfully? " << success << " placed ";
v[j].key.print();
cout << " as child of ";
this->data.key.print();
cout << endl;
*/
	    }
	}
	return success;
    }
};



template <typename keyT>
class MyProcmap {
private:
    int whichmap;
    const ProcessID owner;
    Tree treeMap;


    ProcessID getOwner(const keyT& key) const {
	ProcessID owner;
	treeMap.findOwner(key, &owner);
	return owner;
    };


    void buildTreeMap(vector<TreeCoords> v) {
	sort(v.begin(), v.end());
	int vlen = v.size();
/*
print("sorted vector:");
for (int i = 0; i < vlen; i++) {
    v[i].print();
}
cout << endl;
*/
	if (vlen == 0) throw "empty map!!!";
	treeMap = Tree(v[vlen-1]);
	for (int j = vlen-2; j >= 0; j--) {
	    treeMap.fill(v, j);
	}
    };
	

public:
    MyProcmap() : whichmap(0), owner(0) {};
    MyProcmap(ProcessID owner) : whichmap(0), owner(owner) {};

    MyProcmap(vector<TreeCoords> v) : whichmap(1), owner(1) {
	buildTreeMap(v);
	treeMap.print();
    };

    ProcessID operator()(const keyT& key) const {
        //return madness::hash(key)%2;
	if (whichmap == 0)
	    return owner;
	else
	    return getOwner(key);
    };

};

typedef Node<NodeData,2> NodeD;
typedef DistributedContainer< KeyD,NodeD,MyProcmap<KeyD> > treeT;
//typedef DistributedContainer< KeyD,NodeD > treeT;

void build_tree(treeT& tree, const KeyD& key) {
//    cout << "at beginning of build_tree with key ";
//    key.print();
//    cout << endl;
    NodeData data(1,1,false);  
    NodeD parent(data);
    if (key.n < 5) {
	for (int p=0; p<2; p++) {
	    for (int q=0; q<2; q++) {
		parent.set_child(p+2*q);
		build_tree(tree,KeyD(key.n+1,2*key.L[0]+p,2*key.L[1]+q));
	    }
	}
    }
    else if ((key.n <= 9)&&(key.L[0] == key.L[1])) {
	for (int p=0; p<2; p++) {
	    for (int q=0; q<2; q++) {
		parent.set_child(p+2*q);
		build_tree(tree,KeyD(key.n+1,2*key.L[0]+p,2*key.L[1]+q));
	    }
	}
    }
    tree.insert(key,parent);
}


void print_tree(treeT& tree, const KeyD& key) {
    treeT::iterator it = tree.find(key);
    if (it!=tree.end()) {
	const NodeD& node = it->second;
    	// no longer need iterator
//    	it = tree.end();
	NodeData d = node.get_data();
	for (int i=0; i<(int)key.n; i++) cout << "   ";
	print(key.n,key.L[0],key.L[1],"owner",tree.owner(key),"cost",d.cost,"subcost", d.subcost);

	for (int p = 0; p < 4; p++) {
	    if (node.has_child(p)) { 
		KeyD mykey(key.n+1, 2*key.L[0]+p%2, 2*key.L[1]+p/2);
		print_tree(tree, mykey);
	    }
	}
    }
    else {
	cout << "print_tree: sorry, couldn't find key ";
	key.print();
	cout << endl;
    }
}

typedef int Cost;

Cost computeCost(treeT& tree, const KeyD& key) {
    Cost cost = 0;
    treeT::iterator it = tree.find(key);
    if (it == tree.end()) return cost;

    NodeD node = it->second;
    // no longer need iterator
//    it = tree.end();
    for (unsigned int i = 0; i < node.dim; i++) {
	if (node.has_child(i)) {
	    KeyD k = key.myChild(i);
	    cost += computeCost(tree,k);
	}
    }
    NodeData d = node.get_data();
    cost += d.cost;
    
    d.subcost = cost;
    node.set_data(d);
    tree.erase(key);
    tree.insert(key,node);
    return cost;
}

void meld(treeT& tree, const KeyD& key) {
    Cost cheapest = 0;
    treeT::iterator it = tree.find(key);
    if (it == tree.end()) return;
//    cout << "meld: at beginning, key = ";
//    key.print();
//    cout << endl;

    vector<unsigned int> mylist;

    NodeD node = it->second;
    // no longer need iterator
//    it = tree.end();
    for (unsigned int i = 0; i < node.dim; i++)
    {
	if (node.has_child(i)) {
	    KeyD k = key.myChild(i);
	    treeT::iterator itc = tree.find(k);
            if (itc == tree.end()) return;
            NodeD c = itc->second;
	    // no longer need iterator
//	    itc = tree.end();
            bool haskids = false;
            for (unsigned int j = 0; j < c.dim; j++) {
                if (c.has_child(j)) {
                    haskids = true;
                    break;
                }
            }
            if (!haskids) {
                Cost cost = c.get_data().cost;
                if ((cost < cheapest) || (cheapest == 0)) {
                    cheapest = cost;
                    mylist.clear();
                    mylist.push_back(i);
                }
                else if (cost == cheapest) {
                    mylist.push_back(i);
                }
            }
	}
    }

    if (cheapest == 0) {
	NodeData d = node.get_data();
	d.istaken = false;
	node.set_data(d);
	tree.erase(key);
	tree.insert(key,node);
	return;
    }

    NodeData d = node.get_data();

    for (unsigned int i = 0; i < mylist.size(); i++) {
        d.cost += cheapest;
        tree.erase(key.myChild(mylist[i]));
	node.set_child(mylist[i], false);
//cout << "meld: set child " << mylist[i] << " to be false" << endl;
    }
    d.istaken = false;
    node.set_data(d);
    tree.erase(key);
    tree.insert(key,node);
    treeT::iterator itd = tree.find(key);
    NodeD noded = itd->second;
    // no longer need iterator
//    itd = tree.end();
    NodeData dd = noded.get_data();
//    cout << "meld: at end, node has these values for children " << noded.has_child(0) << ",";
//    cout << noded.has_child(1) << "," << noded.has_child(2) << "," << noded.has_child(3) << endl;
//    cout << "meld: and cost = " << dd.cost << ", subcost = " << dd.subcost << endl;
}

void rollup(treeT tree, KeyD key) {
//    cout << "rollup: at beginning" << endl;
    treeT::iterator it = tree.find(key);
    if (it == tree.end()) return;

//    cout << "rollup: about to get node associated with key ";
//    key.print();
//    cout << endl;
    NodeD node = it->second;
    // no longer need iterator
//    it = tree.end();
    if (!node.has_children()) {
//	cout << "rollup: this node has no children; returning" << endl;
	return; // no rolling to be done here.
    }
//    cout << "rollup: this node has children" << endl;
    bool hasleafchild = false;
    for (unsigned int i = 0; i < node.dim; i++) {
	KeyD k = key.myChild(i);
	treeT::iterator itc = tree.find(k);
	if (itc != tree.end()) {
//	    cout << "rollup: found child ";
//	    k.print();
//	    cout << endl;
	    NodeD c = itc->second;
	    // no longer need iterator
//	    itc = tree.end();
	    if (c.has_children()) {
//		cout << "rollup: child ";
//		k.print();
//		cout << " has children" << endl;
		rollup(tree, k);
	    }
	    else {
//		cout << "rollup: child is leaf" << endl;
		hasleafchild = true;
	    }
	}
    }
    if (hasleafchild) {
//	cout << "rollup: about to meld" << endl;
	meld(tree,key);
    }
    it = tree.find(key);
    node = it->second;
    NodeData d = node.get_data();
    d.istaken = false;
//cout << "rollup: print tree (before setting data and erasing key)" << endl;
//print_tree(tree, key);
//cout << "rollup: end print tree" << endl;
    node.set_data(d);
    tree.erase(key);
    tree.insert(key,node);
//cout << "rollup: print tree" << endl;
//print_tree(tree, key);
//cout << "rollup: end print tree" << endl;
}

Cost fixCost(treeT tree, KeyD key) {
    treeT::iterator it = tree.find(key);
    if (it == tree.end()) return 0;

    NodeD node = it->second;
    // no longer need iterator
//    it = tree.end();
    NodeData d = node.get_data();
    d.subcost = d.cost;
    if (node.has_children())
    {
	for (unsigned int i = 0; i < node.dim; i++)
	{
	    d.subcost += fixCost(tree, key.myChild(i));
	}
    }
    node.set_data(d);
//cout << "about to insert key = ";
//key.print();
//cout << ", ";
//node.get_data().print();
    tree.erase(key);
    tree.insert(key,node);
/*
    treeT::iterator itt = tree.find(key);
    if (itt != tree.end()) {
    	NodeD noded = itt->second;
	// no longer need iterator
//	itt = tree.end();
    	NodeData dd = noded.get_data();
    	cout << "cost and subcost of ";
    	key.print();
    	cout << " = " << dd.cost << ", " << dd.subcost << endl;
    }
    else {
	print("uh oh, no node at all!");
    }
*/
    return d.subcost;
}

Cost computePartitionSize(Cost cost, unsigned int parts) {
    return (Cost) ceil(((double) cost)/((double) parts));
}


Cost makePartition(treeT tree, KeyD key, vector<KeyD>* klist, Cost partitionSize, 
	bool lastPartition, Cost usedUp = 0, bool *atleaf = false);

Cost depthFirstPartition(treeT tree, KeyD key, vector<TreeCoords>* klist, unsigned int npieces, 
	Cost totalcost = 0, Cost *maxcost = 0) {
//print("depthFirstPartition: at very beginning");
    if (totalcost == 0) {
	totalcost = computeCost(tree, key);
    }
//print("depthFirstPartition: totalcost =", totalcost);

    Cost costLeft = totalcost;
    int partsLeft = npieces;
    *maxcost = 0;
    Cost partitionSize = 0;

    for (int i = npieces-1; i >= 0; i--) {
	cout << endl << "Beginning partition number " << i << endl;
	vector<KeyD> tmplist;
	Cost tpart = computePartitionSize(costLeft, partsLeft);
	if (tpart > partitionSize) {
	    partitionSize = tpart;
	}
//print("depthFirstPartition: partitionSize =", partitionSize);
	Cost usedUp = 0;
	bool atleaf = false;
	usedUp = makePartition(tree, key, &tmplist, partitionSize, (i==0), usedUp, &atleaf);
	if (*maxcost < usedUp) *maxcost = usedUp;
	costLeft -= usedUp;
	partsLeft--;
	for (unsigned int j = 0; j < tmplist.size(); j++) {
	    klist->push_back(TreeCoords(KeyD(tmplist[j]), i)); 
	}
    }
    return totalcost;
}

void removeCost(treeT tree, KeyD key, Cost c) {
//cout << "removeCost: key ";
//key.print();
//cout << endl;
    if (((int) key.n) < 0) return;
    treeT::iterator it = tree.find(key);
//print("removeCost: found key");
    if (it == tree.end()) return;
    NodeD node = it->second;
    // no longer need iterator
//    it = tree.end();
    NodeData d = node.get_data();
//print("removeCost: got data");
    d.subcost -= c;
    if (key.n > 0) {
    	removeCost(tree, key.myParent(), c);
    }
//cout << "removeCost: before setting, data = ";
//d.print();
    node.set_data(d);
//cout << "removeCost: after setting, data = ";
//node.get_data().print();
    tree.erase(key);
    tree.insert(key,node);
//cout << "removeCost: after inserting, data = ";
//node.get_data().print();
/*
    treeT::iterator it2 = tree.find(key);
    NodeD node2 = it2->second;
    // no longer need iterator
//    it2 = tree.end();
    NodeData d2 = node2.get_data();
cout << "removeCost: key ";
key.print();
cout << ", retrieved data = ";
d2.print();
cout << endl;
*/
}


Cost makePartition(treeT tree, KeyD key, vector<KeyD>* klist, Cost partitionSize, bool lastPartition, Cost usedUp, bool *atleaf)
{
//    cout << "at beginning of makePartition: atleaf = ";
//    cout << *atleaf << endl;
    double fudgeFactor = 0.1;
    Cost maxAddl = (Cost) (fudgeFactor*partitionSize);

    treeT::iterator it = tree.find(key);
    if (it == tree.end()) { 
	return usedUp;
    }

    NodeD node = it->second;
    NodeData d = node.get_data();

    // no longer need iterator
//    it = tree.end();

//    cout << "data for key ";
//    key.print();
//    cout << ": cost = " << d.cost << ", subcost = " << d.subcost << endl;
//    cout << "partitionSize = " << partitionSize << ", lastPartition = " << lastPartition <<
//	", usedUp = " << usedUp << std::endl;

    if (d.istaken) {
//	cout << "this key is taken" << endl; 
	return usedUp;
    }

//    cout << "back to key ";
//    key.print();
//    cout << endl;

    // if either we're at the last partition, the partition is currently empty
    // and this is a single item, or there is still room in the partition and
    // adding this to it won't go above the fudge factor,
    // then add this piece to the partition.
    if ((lastPartition) || ((usedUp == 0) && (!node.has_children())) || 
	((usedUp < partitionSize) && (d.subcost+usedUp <= partitionSize+maxAddl))) {
	// add to partition
//	cout << "adding to partition ";
//	key.print();
//	cout << endl;
	klist->push_back(KeyD(key));
	d.istaken = true;
	usedUp += d.subcost;
	// REMOVE COST FROM FOREPARENTS (implement this)
	removeCost(tree, key.myParent(), d.subcost);
	node.set_data(d);
	tree.erase(key);
	tree.insert(key,node);
    }
    else if (usedUp < partitionSize) {
	// try this node's children (if any) 
	if (node.has_children()) {
	    for (unsigned int i = 0; i < node.dim; i++) {
	    	if (node.has_child(i)) {
	    	    KeyD k = key.myChild(i);
//		    key.print();
//	    	    std::cout << " recursively calling ";
//	    	    k.print();
//	    	    cout << endl;
	    	    usedUp = makePartition(tree, k, klist, partitionSize, lastPartition, usedUp, atleaf);
	    	    if ((*atleaf) || (usedUp >= partitionSize)) {
//			cout << "at leaf = " << *atleaf << ", usedup >= partitionSize? " << 
//				(usedUp >=partitionSize) << endl;
		        break;
		    }
		}
	    }
	}
	else {
//	    cout << "about to set atleaf = true" << endl;
	    *atleaf = true;
	}
    }
    return usedUp;
}

typedef double CompCost;

CompCost computeCompCost(Cost c, int n) {
    CompCost compcost;
    CompCost cfactor = 0.1, nfactor = 1.0;
    compcost = cfactor*c + nfactor*n;
    return compcost;
}

void findBestPartition(treeT tree, KeyD key, vector<TreeCoords>* klist, unsigned int npieces) {
    bool notdone = true;
    int count = 0;
    vector<vector<TreeCoords> > listoflist;
    vector<TreeCoords> emptylist;
    vector<Cost> costlist;

    listoflist.push_back(emptylist);
    costlist.push_back(0);
    Cost totalCost = 0;

//print("findBestPartition: about to fixCost");

    fixCost(tree, key);
//    print_tree(tree,KeyD(0,0,0));
//print("findBestPartition: about to depthFirstPartition");
    totalCost = depthFirstPartition(tree, key, &listoflist[count], npieces, totalCost, &costlist[count]);
//print("findBestPartition: after depthFirstPartition");
    int size = listoflist[count].size();
    std::cout << "Partitioned tree " << count << ":" << std::endl;
    for (int i = 0; i < size; i++)
	listoflist[count][i].print();
    std::cout << "Max cost for this tree = " << costlist[count] << std::endl;
    std::cout << std::endl;
    if (listoflist[count].size() < npieces)
	notdone = false;
    count++;

    while (notdone) {
	fixCost(tree, key); 
	rollup(tree, key);
//	print_tree(tree,KeyD(0,0,0));
	listoflist.push_back(emptylist);
	costlist.push_back(0);
	depthFirstPartition(tree, key, &listoflist[count], npieces, totalCost, &costlist[count]);
	int size = listoflist[count].size();
	std::cout << "Partitioned tree " << count << ":" << std::endl;
	for (int i = 0; i < size; i++)
	    listoflist[count][i].print();
	std::cout << "Max cost for this tree = " << costlist[count] << std::endl;
	std::cout << std::endl;
	
    	treeT::iterator it = tree.find(key);
    	if (it == tree.end()) return;
    	NodeD node = it->second;
    	// no longer need iterator
//    	it = tree.end();
	if (!(node.has_children()) || (listoflist[count].size() < npieces)) {
	    notdone = false;
	}
	if (listoflist[count].size() < npieces) {
	    listoflist.erase(listoflist.begin()+count);
	    break;
	}
	count++;
    }
    unsigned int shortestList = 0, SL_index, LB_index;
    Cost loadBalCost = 0;
    vector<unsigned int> len;
    for (int i = 0; i < count; i++) {
	len.push_back(listoflist[i].size());
	if ((len[i] < shortestList) || (shortestList == 0)) {
	    shortestList = len[i];
	    SL_index = i;
	}
	else if ((len[i] == shortestList) && (costlist[i] < costlist[SL_index])) {
	// all things being equal, prefer better balance
	    shortestList = len[i];
	    SL_index = i;
	}
	if ((costlist[i] < loadBalCost) || (loadBalCost == 0)) {
	    loadBalCost = costlist[i];
	    LB_index = i;
	}
	else if ((costlist[i] == loadBalCost) && (len[i] < listoflist[LB_index].size())) {
	// all things being equal, prefer fewer cuts
	    loadBalCost = costlist[i];
	    LB_index = i;
	}
    }

    cout << "The load balance with the fewest broken links has cost " << costlist[SL_index] << 
	", and " << shortestList-1 << " broken links" << std::endl;
    for (unsigned int i = 0; i < shortestList; i++) {
	listoflist[SL_index][i].print();
    }
    std::cout << std::endl;
    cout << "The load balance with the best balance has cost " << loadBalCost << ", and " <<
	listoflist[LB_index].size()-1 << " broken links" << std::endl;
    for (unsigned int i = 0; i < listoflist[LB_index].size(); i++) {
	listoflist[LB_index][i].print();
    }
    std::cout << std::endl;

    CompCost ccleast = 0;
    int cc_index;
    for (int i = 0; i < count; i++) {
	CompCost cctmp = computeCompCost(costlist[i], len[i]-1);
	if ((i==0) || (cctmp < ccleast)) {
	    ccleast = cctmp;
	    cc_index = i;
	}
    }
    cout << "The load balance with the best overall computational cost has cost " <<
	costlist[cc_index] << " and " << len[cc_index]-1 << " broken links" << std::endl;
    for (unsigned int i = 0; i < len[cc_index]; i++) {
	listoflist[cc_index][i].print();
    }
    for (unsigned int i = 0; i < len[cc_index]; i++) {
	klist->push_back(listoflist[cc_index][i]);
    }
}

int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);
    redirectio(world);
//    xterm_debug("world", 0);
    ProcessID me = world.rank();


    try {
	vector<TreeCoords> v;
	v.push_back(TreeCoords(KeyD(0,0,0),0));
	v.push_back(TreeCoords(KeyD(1,0,1),1));
	v.push_back(TreeCoords(KeyD(1,1,1),1));
/*
	v.push_back(TreeCoords(KeyD(1,0,1),2));
	v.push_back(TreeCoords(KeyD(0,0,0),6));
	v.push_back(TreeCoords(KeyD(1,1,0),4));
	v.push_back(TreeCoords(KeyD(1,0,0),1));
	v.push_back(TreeCoords(KeyD(2,2,2),5));
	v.push_back(TreeCoords(KeyD(2,1,0),0));
	v.push_back(TreeCoords(KeyD(2,2,0),3));
*/
	KeyD root(0,0,0);
	treeT tree(world,MyProcmap<KeyD>(v));
//	treeT tree(world,MyProcmap<KeyD>(1));
	print("Made tree");
	if (me == 0) { 
	    print("About to build tree");
	    build_tree(tree,root);
	}
	print("Built tree");
	world.gop.fence();
        print("Done fencing");
	print("Now we're going to print the tree");
	print("");
	if (me == 1) {
	    print("About to print tree");
//	    print_tree(tree,root);
	    print("Printed tree");
	}
	print("Done printing tree");
	print("");
	if (me == 1) {
//	    Cost cost = computeCost(tree,KeyD(0,0,0));
//	    print("cost of tree =", cost);
	    vector<TreeCoords> klist;
	    unsigned int npieces = world.nproc();
	    findBestPartition(tree, KeyD(0,0,0), &klist, npieces);
	}
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

    print("final fence");
    world.gop.fence();
    print("done final fence");
    MPI::Finalize();
    return 0;
}
