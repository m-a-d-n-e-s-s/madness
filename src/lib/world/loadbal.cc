#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <world/loadbal.h>

//using namespace madness;
using namespace std;

namespace madness {

typedef int Cost;
typedef double CompCost;


template <int D>
void build_tree(typename DClass<D>::treeT& tree, typename DClass<D>::KeyDConst& key) {
    NodeData data(1,1,false);  
    typename DClass<D>::NodeD parent(data);
    int twotoD = (int) power<D>();
    if (key.level() < 3) {
	for (int i = 0; i < twotoD; i++) {
	    parent.set_child(i);
	}
	for (KeyChildIterator<D> kit(key); kit; ++kit) {
	    build_tree<D>(tree, kit.key());
	}
    }
//    else if ((key.level() <= 9)&&(key.translation()[0] == key.translation()[1])) {}
    else if ((key.level() <= 5)&&(key.translation()[0] == key.translation()[1])) {
	for (int i = 0; i < twotoD; i++) {
	    parent.set_child(i);
	}
	for (KeyChildIterator<D> kit(key); kit; ++kit) {
	    build_tree<D>(tree, kit.key());
	}
    }
    tree.insert(key,parent);
}


template <int D>
void print_tree(typename DClass<D>::treeT& tree, typename DClass<D>::KeyDConst& key) {
    typename DClass<D>::treeT::iterator it = tree.find(key);
    if (it!=tree.end()) {
	typename DClass<D>::NodeDConst& node = it->second;
	NodeData d = node.get_data();
	for (int i=0; i<(int)key.level(); i++) cout << "   ";
//	print(key,"owner",tree.owner(key),d);
	print(key,"owner",tree.owner(key),node);
	for (KeyChildIterator<D> kit(key); kit; ++kit) {
	    print_tree<D>(tree, kit.key());
	}
    }
}


template <int D>
Cost computeCost(typename DClass<D>::treeT& tree, typename DClass<D>::KeyDConst& key) {
    Cost cost = 0;
    typename DClass<D>::treeT::iterator it = tree.find(key);
    if (it == tree.end()) return cost;

    typename DClass<D>::NodeD node = it->second;
    for (KeyChildIterator<D> kit(key); kit; ++kit) {
	cost += computeCost<D>(tree,kit.key());
    }
    NodeData d = node.get_data();
    cost += d.cost;
    
    d.subcost = cost;
    node.set_data(d);
    tree.insert(key,node);
    return cost;
}


template <int D>
void meld(typename DClass<D>::treeT& tree, typename DClass<D>::KeyDConst& key) {
//    print("meld: at beginning, finding key", key);
    Cost cheapest = 0;
    typename DClass<D>::treeT::iterator it = tree.find(key);
    if (it == tree.end()) return;

    vector<unsigned int> mylist;

    typename DClass<D>::NodeD node = it->second;
    unsigned int i = 0;
//    print("meld: about to iterate over children of key", key);
    for (KeyChildIterator<D> kit(key); kit; ++kit) {
//    	print("    meld: iterating over child", i);
	if (node.has_child(i)) {
	    typename DClass<D>::treeT::iterator itc = tree.find(kit.key());
            if (itc == tree.end()) return;
            typename DClass<D>::NodeD c = itc->second;
            if (!c.has_children()) {
//		print("    meld: child",i,"has no children");
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
	i++;
    }

    if (cheapest == 0) {
//	print("meld: this node has no leaf children");
	NodeData d = node.get_data();
	d.istaken = false;
	node.set_data(d);
	tree.insert(key,node);
	return;
    }

    NodeData d = node.get_data();

    i = 0;
    int j = 0, mlsize = mylist.size();
    for (KeyChildIterator<D> kit(key); kit; ++kit) {
	if (mylist[j] == i) {
//	    print("meld: found a match, mylist[",j,"] =",i);
	    tree.erase(kit.key());
	    node.set_child(mylist[j], false);
	    d.cost += cheapest;
	    j++;
	}
	i++;
	if (j == mlsize) break;
    }
    d.istaken = false;
    node.set_data(d);
    tree.insert(key,node);
//    print("meld: inserted node back into tree; goodbye!");
}


template <int D>
void rollup(typename DClass<D>::treeT tree, typename DClass<D>::KeyD key) {
//    print("rollup: at beginning");
    typename DClass<D>::treeT::iterator it = tree.find(key);
    if (it == tree.end()) return;

//    print("rollup: about to get node associated with key",key);
    typename DClass<D>::NodeD node = it->second;
    if (!node.has_children()) {
//	print("rollup: this node has no children; returning");
	return; // no rolling to be done here.
    }
//    print("rollup: this node has children");
    bool hasleafchild = false;
    for (KeyChildIterator<D> kit(key); kit; ++kit) {
	typename DClass<D>::treeT::iterator itc = tree.find(kit.key());
	if (itc != tree.end()) {
//	    print("rollup: found child", kit.key());
	    typename DClass<D>::NodeD c = itc->second;
	    if (!c.has_children()) {
//		print("rollup: child is leaf");
		hasleafchild = true;
		break;
	    }
	    else {
//		print("rollup: child", kit.key(), "has children");
	    }
	}
    }
    if (hasleafchild) {
//	print("rollup: about to meld with key",key);
	meld<D>(tree,key);
    }
    for (KeyChildIterator<D> kit(key); kit; ++kit) {
	typename DClass<D>::treeT::iterator itc = tree.find(kit.key());
	if (itc != tree.end()) {
//	    print("rollup: found child", kit.key());
	    typename DClass<D>::NodeD c = itc->second;
	    if (c.has_children()) {
//		print("rollup: child", kit.key(), "has children");
		rollup<D>(tree, kit.key());
	    }
	}
    }
    it = tree.find(key);
    node = it->second;
    NodeData d = node.get_data();
    if (d.istaken) {
    	d.istaken = false;
    	node.set_data(d);
    	tree.insert(key,node);
    }
}

template <int D>
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
	for (KeyChildIterator<D> kit(key); kit; ++kit) {
	    d.subcost += fixCost<D>(tree, kit.key());
	}
    }
    node.set_data(d);
//cout << "about to insert key = ";
//key.print();
//cout << ", ";
//node.get_data().print();
    tree.insert(key,node);
//print("fixCost: inserted node");
    return d.subcost;
}

Cost computePartitionSize(Cost cost, unsigned int parts) {
    return (Cost) ceil(((double) cost)/((double) parts));
}


template <int D>
Cost depthFirstPartition(typename DClass<D>::treeT tree, typename DClass<D>::KeyD key, 
	vector<typename DClass<D>::TreeCoords>* klist, unsigned int npieces, 
	Cost totalcost = 0, Cost *maxcost = 0) {
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

template <int D>
void removeCost(typename DClass<D>::treeT tree, typename DClass<D>::KeyD key, Cost c) {
//cout << "removeCost: key ";
//key.print();
//cout << endl;
    if (((int) key.level()) < 0) return;
    typename DClass<D>::treeT::iterator it = tree.find(key);
//print("removeCost: found key");
    if (it == tree.end()) return;
    typename DClass<D>::NodeD node = it->second;
    NodeData d = node.get_data();
//print("removeCost: got data");
    d.subcost -= c;
    if (key.level() > 0) {
    	removeCost<D>(tree, key.parent(), c);
    }
//cout << "removeCost: before setting, data = ";
//d.print();
    node.set_data(d);
//cout << "removeCost: after setting, data = ";
//node.get_data().print();
    tree.insert(key,node);
//cout << "removeCost: after inserting, data = ";
//node.get_data().print();
}


template <int D>
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
	removeCost<D>(tree, key.parent(), d.subcost);
	node.set_data(d);
	tree.insert(key,node);
    }
    else if (usedUp < partitionSize) {
	// try this node's children (if any) 
	if (node.has_children()) {
	    int i = 0;
	    for (KeyChildIterator<D> kit(key); kit; ++kit) {
		if (node.has_child(i)) {
//		    print(key, "recursively calling", kit.key());
		    usedUp = makePartition<D>(tree, kit.key(), klist, partitionSize, lastPartition,
			usedUp, atleaf);
		    if ((*atleaf) || (usedUp >= partitionSize)) {
			break;
		    }
		}
		i++;
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

template <int D>
void findBestPartition(typename DClass<D>::treeT tree, typename DClass<D>::KeyD key, 
	vector<typename DClass<D>::TreeCoords>* klist, unsigned int npieces) {
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
//    print_tree<D>(tree,typename DClass<D>::KeyD(0));
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
//    	print_tree<D>(tree,typename DClass<D>::KeyD(0));
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


template <int D>
void migrate_data(typename DClass<D>::treeT tfrom, typename DClass<D>::treeT tto, typename DClass<D>::KeyD key) {
//    print("migrate_data: about to migrate from", key);
    typename DClass<D>::treeT::iterator it = tfrom.find(key);
    if (it == tfrom.end()) return;

//    print("migrate_data: about to get node for", key);
    typename DClass<D>::NodeD node = it->second;

    if (node.has_children()) {
//	print("migrate_data:", key, "has children");
        for (KeyChildIterator<D> kit(key); kit; ++kit) {
	    migrate_data<D>(tfrom, tto, kit.key());
//	    print("migrate_data: back from", kit.key());
	}
    }
//    print("migrate_data: about to insert", key);
    tto.insert(key, node);
//    print("migrate_data: just inserted", key);
}


template <int D>
void migrate(typename DClass<D>::treeT tfrom, typename DClass<D>::treeT tto) {
    typename DClass<D>::KeyD root(0);
//    print("about to migrate from root,", root);
    migrate_data<D>(tfrom, tto, root);
}


// Explicit instantiations for D=2 and D=3

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

template void build_tree<3>(DClass<3>::treeT& tree, DClass<3>::KeyDConst& key);
template void print_tree<3>(DClass<3>::treeT& tree, DClass<3>::KeyDConst& key);
template Cost computeCost<3>(DClass<3>::treeT& tree, DClass<3>::KeyDConst& key);
template void meld<3>(DClass<3>::treeT& tree, DClass<3>::KeyDConst& key);
template void rollup<3>(DClass<3>::treeT tree, DClass<3>::KeyD key);
template Cost fixCost<3>(DClass<3>::treeT tree, DClass<3>::KeyD key);
template Cost makePartition<3>(DClass<3>::treeT tree, DClass<3>::KeyD key, 
	vector<DClass<3>::KeyD>* klist, Cost partitionSize, bool lastPartition, 
	Cost usedUp, bool *atleaf);
template Cost depthFirstPartition<3>(DClass<3>::treeT tree, DClass<3>::KeyD key,
        vector<DClass<3>::TreeCoords>* klist, unsigned int npieces, Cost totalcost, 
	Cost *maxcost);
template void removeCost<3>(DClass<3>::treeT tree, DClass<3>::KeyD key, Cost c);
template void findBestPartition<3>(DClass<3>::treeT tree, DClass<3>::KeyD key,
        vector<DClass<3>::TreeCoords>* klist, unsigned int npieces);
template void migrate_data<3>(DClass<3>::treeT tfrom, DClass<3>::treeT tto, 
	DClass<3>::KeyD key);
template void migrate<3>(DClass<3>::treeT tfrom, DClass<3>::treeT tto);
}
