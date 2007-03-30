#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <world/loadbal.h>

//using namespace madness;
using namespace std;

namespace madness {

typedef int Cost;
typedef double CompCost;


/*
// To enable gnu hash_map to work
namespace std {
    template <> struct equal_to<typename madness::DClass<3>::KeyD > {
        bool operator()(typename madness::DClass<3>::KeyDConst& a, typename madness::DClass<3>::KeyDConst& b) const {
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
    template <> struct hash<typename madness::DClass<3>::KeyD > { 
        std::size_t operator()(typename madness::DClass<3>::KeyDConst& key) const {
            return key.hashval;
        };
    };
}

    template <> struct equal_to<typename madness::DClass<2>::KeyD > {
        bool operator()(typename madness::DClass<2>::KeyDConst& a, typename madness::DClass<2>::KeyDConst& b) const {
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
    template <> struct hash<madness::DClass<2>::KeyD > { 
        std::size_t operator()(madness::DClass<2>::KeyDConst& key) const {
            return key.hashval;
        };
    };
}
*/


template <unsigned int D>
void build_tree(typename DClass<D>::treeT& tree, typename DClass<D>::KeyDConst& key) {
    NodeData data(1,1,false);  
    typename DClass<D>::NodeD parent(data);
    int twotoD = (int) pow(2.0,(int)D);
    if (key.n < 2) {
	for (int i = 0; i < twotoD; i++) {
	    parent.set_child(i);
	    build_tree<D>(tree, key.myChild(i));
	}
    }
/*
    else if ((key.n <= 9)&&(key.L[0] == key.L[1])) {
	for (int i = 0; i < twotoD; i++) {
	    parent.set_child(i);
	    build_tree<D>(tree, key.myChild(i));
	}
    }
*/
/*
    if (key.n < 5) {
	for (int p=0; p<2; p++) {
	    for (int q=0; q<2; q++) {
		parent.set_child(p+2*q);
		build_tree<D>(tree,typename DClass<D>::KeyD(key.n+1,2*key.L[0]+p,2*key.L[1]+q));
	    }
	}
    }
    else if ((key.n <= 9)&&(key.L[0] == key.L[1])) {
	for (int p=0; p<2; p++) {
	    for (int q=0; q<2; q++) {
		parent.set_child(p+2*q);
		build_tree<D>(tree,typename DClass<D>::KeyD(key.n+1,2*key.L[0]+p,2*key.L[1]+q));
	    }
	}
    }
*/
    tree.insert(key,parent);
}


template <unsigned int D>
void print_tree(typename DClass<D>::treeT& tree, typename DClass<D>::KeyDConst& key) {
    typename DClass<D>::treeT::iterator it = tree.find(key);
    if (it!=tree.end()) {
	typename DClass<D>::NodeDConst& node = it->second;
	NodeData d = node.get_data();
	for (int i=0; i<(int)key.n; i++) cout << "   ";
	print(key.n,key.L[0],key.L[1],"owner",tree.owner(key),"cost",d.cost,"subcost", d.subcost);

	for (int p = 0; p < (int)node.dim; p++) {
	    typename DClass<D>::KeyD mykey = key.myChild(p);
	    print_tree<D>(tree, mykey);
	}
    }
}


template <unsigned int D>
Cost computeCost(typename DClass<D>::treeT& tree, typename DClass<D>::KeyDConst& key) {
    Cost cost = 0;
    typename DClass<D>::treeT::iterator it = tree.find(key);
    if (it == tree.end()) return cost;

    typename DClass<D>::NodeD node = it->second;
    for (unsigned int i = 0; i < node.dim; i++) {
	typename DClass<D>::KeyD k = key.myChild(i);
	cost += computeCost<D>(tree,k);
    }
    NodeData d = node.get_data();
    cost += d.cost;
    
    d.subcost = cost;
    node.set_data(d);
//    tree.erase(key);
    tree.insert(key,node);
    return cost;
}

template <unsigned int D>
void meld(typename DClass<D>::treeT& tree, typename DClass<D>::KeyDConst& key) {
    Cost cheapest = 0;
    typename DClass<D>::treeT::iterator it = tree.find(key);
    if (it == tree.end()) return;
//    cout << "meld: at beginning, key = ";
//    key.print();
//    cout << endl;

    vector<unsigned int> mylist;

    typename DClass<D>::NodeD node = it->second;
    for (unsigned int i = 0; i < node.dim; i++)
    {
	if (node.has_child(i)) {
	    typename DClass<D>::KeyD k = key.myChild(i);
	    typename DClass<D>::treeT::iterator itc = tree.find(k);
            if (itc == tree.end()) return;
            typename DClass<D>::NodeD c = itc->second;
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
//	tree.erase(key);
	tree.insert(key,node);
	return;
    }

    NodeData d = node.get_data();

    for (unsigned int i = 0; i < mylist.size(); i++) {
        d.cost += cheapest;
//        tree.erase(key.myChild(mylist[i]));
	node.set_child(mylist[i], false);
//cout << "meld: set child " << mylist[i] << " to be false" << endl;
    }
    d.istaken = false;
    node.set_data(d);
//    tree.erase(key);
    tree.insert(key,node);
}

template <unsigned int D>
void rollup(typename DClass<D>::treeT tree, typename DClass<D>::KeyD key) {
//    cout << "rollup: at beginning" << endl;
    typename DClass<D>::treeT::iterator it = tree.find(key);
    if (it == tree.end()) return;

//    cout << "rollup: about to get node associated with key ";
//    key.print();
//    cout << endl;
    typename DClass<D>::NodeD node = it->second;
    if (!node.has_children()) {
//	cout << "rollup: this node has no children; returning" << endl;
	return; // no rolling to be done here.
    }
//    cout << "rollup: this node has children" << endl;
    bool hasleafchild = false;
    for (unsigned int i = 0; i < node.dim; i++) {
	typename DClass<D>::KeyD k = key.myChild(i);
	typename DClass<D>::treeT::iterator itc = tree.find(k);
	if (itc != tree.end()) {
//	    cout << "rollup: found child ";
//	    k.print();
//	    cout << endl;
	    typename DClass<D>::NodeD c = itc->second;
	    if (c.has_children()) {
//		cout << "rollup: child ";
//		k.print();
//		cout << " has children" << endl;
		rollup<D>(tree, k);
	    }
	    else {
//		cout << "rollup: child is leaf" << endl;
		hasleafchild = true;
	    }
	}
    }
    if (hasleafchild) {
//	cout << "rollup: about to meld" << endl;
	meld<D>(tree,key);
    }
    it = tree.find(key);
    node = it->second;
    NodeData d = node.get_data();
    d.istaken = false;
//cout << "rollup: print tree (before setting data and erasing key)" << endl;
//print_tree(tree, key);
//cout << "rollup: end print tree" << endl;
    node.set_data(d);
//    tree.erase(key);
    tree.insert(key,node);
//cout << "rollup: print tree" << endl;
//print_tree(tree, key);
//cout << "rollup: end print tree" << endl;
}

template <unsigned int D>
Cost fixCost(typename DClass<D>::treeT tree, typename DClass<D>::KeyD key) {
//    cout << "fixCost: key = ";
//    key.print();
//    print(" is about to be looked for");
    typename DClass<D>::treeT::iterator it = tree.find(key);
//    cout << "fixCost: key = ";
//    key.print();
//    print(" was found (looked for),", (it == tree.end()));
    if (it == tree.end()) return 0;
//    print("fixCost: tree it was found (exists)");

    typename DClass<D>::NodeD node = it->second;
//    print("fixCost: got node");
    NodeData d = node.get_data();
//    print("fixCost: got data from node");
    d.subcost = d.cost;
//    print("fixCost: assigned node cost to subcost");
    if (node.has_children())
    {
//	print("fixCost: node has children");
	for (unsigned int i = 0; i < node.dim; i++)
	{
	    d.subcost += fixCost<D>(tree, key.myChild(i));
//	    typename DClass<D>::KeyD child = key.myChild(i);
//	    cout << "fixCost: about to call on child ";
//	    child.print();
//	    cout << endl;
//	    d.subcost += fixCost<D>(tree, child);
//	    print("fixCost: returned from recursive call of fixCost");
	}
    }
    node.set_data(d);
//cout << "about to insert key = ";
//key.print();
//cout << ", ";
//node.get_data().print();
//tree.erase(key);
    tree.insert(key,node);
//print("fixCost: inserted node");
    return d.subcost;
}

Cost computePartitionSize(Cost cost, unsigned int parts) {
    return (Cost) ceil(((double) cost)/((double) parts));
}


template <unsigned int D>
Cost depthFirstPartition(typename DClass<D>::treeT tree, typename DClass<D>::KeyD key, vector<typename DClass<D>::TreeCoords>* klist, 
	unsigned int npieces, Cost totalcost = 0, Cost *maxcost = 0) {
//print("depthFirstPartition: at very beginning");
    if (totalcost == 0) {
	totalcost = computeCost<D>(tree, key);
    }
//print("depthFirstPartition: totalcost =", totalcost);

    Cost costLeft = totalcost;
    int partsLeft = npieces;
    *maxcost = 0;
    Cost partitionSize = 0;

    for (int i = npieces-1; i >= 0; i--) {
	cout << endl << "Beginning partition number " << i << endl;
	vector<typename DClass<D>::KeyD> tmplist;
	Cost tpart = computePartitionSize(costLeft, partsLeft);
	if (tpart > partitionSize) {
	    partitionSize = tpart;
	}
//print("depthFirstPartition: partitionSize =", partitionSize);
	Cost usedUp = 0;
	bool atleaf = false;
	usedUp = makePartition<D>(tree, key, &tmplist, partitionSize, (i==0), usedUp, &atleaf);
	if (*maxcost < usedUp) *maxcost = usedUp;
	costLeft -= usedUp;
	partsLeft--;
	for (unsigned int j = 0; j < tmplist.size(); j++) {
	    klist->push_back(typename DClass<D>::TreeCoords(typename DClass<D>::KeyD(tmplist[j]), i)); 
	}
    }
    return totalcost;
}

template <unsigned int D>
void removeCost(typename DClass<D>::treeT tree, typename DClass<D>::KeyD key, Cost c) {
//cout << "removeCost: key ";
//key.print();
//cout << endl;
    if (((int) key.n) < 0) return;
    typename DClass<D>::treeT::iterator it = tree.find(key);
//print("removeCost: found key");
    if (it == tree.end()) return;
    typename DClass<D>::NodeD node = it->second;
    NodeData d = node.get_data();
//print("removeCost: got data");
    d.subcost -= c;
    if (key.n > 0) {
    	removeCost<D>(tree, key.myParent(), c);
    }
//cout << "removeCost: before setting, data = ";
//d.print();
    node.set_data(d);
//cout << "removeCost: after setting, data = ";
//node.get_data().print();
//    tree.erase(key);
    tree.insert(key,node);
//cout << "removeCost: after inserting, data = ";
//node.get_data().print();
}


template <unsigned int D>
Cost makePartition(typename DClass<D>::treeT tree, typename DClass<D>::KeyD key, vector<typename DClass<D>::KeyD>* klist, 
	Cost partitionSize, bool lastPartition, Cost usedUp, bool *atleaf) {
//    cout << "at beginning of makePartition: atleaf = ";
//    cout << *atleaf << endl;
    double fudgeFactor = 0.1;
    Cost maxAddl = (Cost) (fudgeFactor*partitionSize);

    typename DClass<D>::treeT::iterator it = tree.find(key);
    if (it == tree.end()) { 
	return usedUp;
    }

    typename DClass<D>::NodeD node = it->second;
    NodeData d = node.get_data();


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
	klist->push_back(typename DClass<D>::KeyD(key));
	d.istaken = true;
	usedUp += d.subcost;
	// REMOVE COST FROM FOREPARENTS (implement this)
	removeCost<D>(tree, key.myParent(), d.subcost);
	node.set_data(d);
//	tree.erase(key);
	tree.insert(key,node);
    }
    else if (usedUp < partitionSize) {
	// try this node's children (if any) 
	if (node.has_children()) {
	    for (unsigned int i = 0; i < node.dim; i++) {
	    	if (node.has_child(i)) {
	    	    typename DClass<D>::KeyD k = key.myChild(i);
//		    key.print();
//	    	    std::cout << " recursively calling ";
//	    	    k.print();
//	    	    cout << endl;
	    	    usedUp = makePartition<D>(tree, k, klist, partitionSize, lastPartition, usedUp, atleaf);
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


CompCost computeCompCost(Cost c, int n) {
    CompCost compcost;
    CompCost cfactor = 0.1, nfactor = 1.0;
    compcost = cfactor*c + nfactor*n;
    return compcost;
}

template <unsigned int D>
void findBestPartition(typename DClass<D>::treeT tree, typename DClass<D>::KeyD key, vector<typename DClass<D>::TreeCoords>* klist, unsigned int npieces) {
    bool notdone = true;
    int count = 0;
    vector<vector<typename DClass<D>::TreeCoords> > listoflist;
    vector<typename DClass<D>::TreeCoords> emptylist;
    vector<Cost> costlist;

    listoflist.push_back(emptylist);
    costlist.push_back(0);
    Cost totalCost = 0;

//print("findBestPartition: about to fixCost");

    fixCost<D>(tree, key);
    print_tree<D>(tree,typename DClass<D>::KeyD(0,0,0));
print("findBestPartition: about to depthFirstPartition");
    totalCost = depthFirstPartition<D>(tree, key, &listoflist[count], npieces, totalCost, &costlist[count]);
//print("findBestPartition: after depthFirstPartition");
    int size = listoflist[count].size();
    cout << "Partitioned tree " << count << ":" << endl;
    for (int i = 0; i < size; i++)
	listoflist[count][i].print();
    cout << "Max cost for this tree = " << costlist[count] << endl;
    cout << endl;
    if (listoflist[count].size() < npieces)
	notdone = false;
    count++;

    while (notdone) {
	fixCost<D>(tree, key); 
	rollup<D>(tree, key);
	print_tree<D>(tree,typename DClass<D>::KeyD(0,0,0));
	listoflist.push_back(emptylist);
	costlist.push_back(0);
	depthFirstPartition<D>(tree, key, &listoflist[count], npieces, totalCost, &costlist[count]);
	int size = listoflist[count].size();
	cout << "Partitioned tree " << count << ":" << endl;
	for (int i = 0; i < size; i++)
	    listoflist[count][i].print();
	cout << "Max cost for this tree = " << costlist[count] << endl;
	cout << endl;
	
    	typename DClass<D>::treeT::iterator it = tree.find(key);
    	if (it == tree.end()) return;
    	typename DClass<D>::NodeD node = it->second;
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
	", and " << shortestList-1 << " broken links" << endl;
    for (unsigned int i = 0; i < shortestList; i++) {
	listoflist[SL_index][i].print();
    }
    cout << endl;
    cout << "The load balance with the best balance has cost " << loadBalCost << ", and " <<
	listoflist[LB_index].size()-1 << " broken links" << endl;
    for (unsigned int i = 0; i < listoflist[LB_index].size(); i++) {
	listoflist[LB_index][i].print();
    }
    cout << endl;

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
	costlist[cc_index] << " and " << len[cc_index]-1 << " broken links" << endl;
    for (unsigned int i = 0; i < len[cc_index]; i++) {
	listoflist[cc_index][i].print();
    }
    for (unsigned int i = 0; i < len[cc_index]; i++) {
	klist->push_back(listoflist[cc_index][i]);
    }
}


template <unsigned int D>
void migrate_data(typename DClass<D>::treeT tfrom, typename DClass<D>::treeT tto, typename DClass<D>::KeyD key) {
    typename DClass<D>::treeT::iterator it = tfrom.find(key);
    if (it == tfrom.end()) return;

    typename DClass<D>::NodeD node = it->second;

    if (node.has_children()) {
	for (unsigned int i = 0; i < node.dim; i++) {
	    typename DClass<D>::KeyD child = key.myChild(i);
	    migrate_data<D>(tfrom, tto, child);
	}
    }
    tto.insert(key, node);
}


template <unsigned int D>
void migrate(typename DClass<D>::treeT tfrom, typename DClass<D>::treeT tto) {
    typename DClass<D>::KeyD root(0,0,0);
    migrate_data<D>(tfrom, tto, root);
}


// convert tree from templated form to tree to be used for load balancing
/*
template <typename Q, unsigned int N>
void convert_node(DistributedContainer<KeyD,Node<Q,N>,MyProcmap<KeyD> > orig, treeT skel, KeyD key) {
    typename DistributedContainer<KeyD,Node<Q,N>,MyProcmap<KeyD> >::iterator it = orig.find(key);

    if (it == orig.end()) return;

    Node<Q,N> node = it->second;

    if (node.has_children()) {
	for (unsigned int i = 0; i < node.dim; i++) {
	    KeyD child = key.myChild(i);
	    convert_node(orig, skel, child);
	}
    }
    
    DClass<D>::NodeD noded(NodeData());
    skel.insert(key, noded);
}

template <typename Q, unsigned int N>
void convert_tree(DistributedContainer<KeyD,Node<Q,N>,MyProcmap<KeyD> > orig, treeT skel) {
    KeyD root(0,0,0);
    convert_node<Q,N>(orig, skel, root);
}
*/

// Explicit instantiations for D=2

//template struct DClass<2>;

template void build_tree<2>(DClass<2>::treeT& tree, DClass<2>::KeyDConst& key);
template void print_tree<2>(DClass<2>::treeT& tree, DClass<2>::KeyDConst& key);
template Cost computeCost<2>(DClass<2>::treeT& tree, DClass<2>::KeyDConst& key);
template void meld<2>(DClass<2>::treeT& tree, DClass<2>::KeyDConst& key);
template void rollup<2>(DClass<2>::treeT tree, DClass<2>::KeyD key);
template Cost fixCost<2>(DClass<2>::treeT tree, DClass<2>::KeyD key);
template Cost makePartition<2>(DClass<2>::treeT tree, DClass<2>::KeyD key, 
	vector<DClass<2>::KeyD>* klist, Cost partitionSize, bool lastPartition, 
	Cost usedUp, bool *atleaf);
template Cost depthFirstPartition<2>(DClass<2>::treeT tree, DClass<2>::KeyD key,
        vector<DClass<2>::TreeCoords>* klist, unsigned int npieces, Cost totalcost, 
	Cost *maxcost);
template void removeCost<2>(DClass<2>::treeT tree, DClass<2>::KeyD key, Cost c);
template void findBestPartition<2>(DClass<2>::treeT tree, DClass<2>::KeyD key,
        vector<DClass<2>::TreeCoords>* klist, unsigned int npieces);
template void migrate_data<2>(DClass<2>::treeT tfrom, DClass<2>::treeT tto, 
	DClass<2>::KeyD key);
template void migrate<2>(DClass<2>::treeT tfrom, DClass<2>::treeT tto);
}
