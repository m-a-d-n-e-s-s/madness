#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/loadbal.h>

using namespace std;

namespace madness {

typedef int Cost;
typedef double CompCost;

template <typename T, int D, typename Pmap>
vector<typename DClass<D>::TreeCoords> LoadBalImpl<T,D,Pmap>::findBestPartition() { 
    unsigned int npieces = this->world().nproc(); 	// check this!!!
    vector<typename DClass<D>::TreeCoords> klist;
    bool notdone = true;
    int count = 0;
    vector<vector<typename DClass<D>::TreeCoords> > listoflist;
    vector<typename DClass<D>::TreeCoords> emptylist;
    vector<Cost> costlist;

    listoflist.push_back(emptylist);
    costlist.push_back(0);
    Cost totalCost = 0;

//print("findBestPartition: about to fixCost");

    typename DClass<D>::KeyD root(0);
//    this->skeltree.fixCost<D>(root);
    this->skeltree->template fixCost(root);
print("findBestPartition: about to depthFirstPartition");
    this->skeltree->print(root);
    totalCost = this->skeltree->template depthFirstPartition(root, &listoflist[count], npieces, 
	totalCost, &costlist[count]);
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
//	this->skeltree.fixCost<D>(root); 
	this->skeltree->template fixCost(root); 
//	this->skeltree.rollup<D>(root);
	this->skeltree->template rollup(root);
	listoflist.push_back(emptylist);
	costlist.push_back(0);
//	this->skeltree.depthFirstPartition<D>(root, &listoflist[count], npieces, totalCost, &costlist[count]);
	this->skeltree->template depthFirstPartition(root, &listoflist[count], npieces, totalCost, &costlist[count]);
	int size = listoflist[count].size();
	cout << "Partitioned tree " << count << ":" << endl;
	for (int i = 0; i < size; i++)
	    listoflist[count][i].print();
	cout << "Max cost for this tree = " << costlist[count] << endl;
	cout << endl;
	
//    	typename DClass<D>::treeT::iterator it = this->skeltree.find(root);
    	typename DClass<D>::treeT::iterator it = this->skeltree->find(root);
//    	if (it == this->skeltree.end()) return klist;
    	if (it == this->skeltree->end()) return klist;
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
	klist.push_back(listoflist[cc_index][i]);
    }

    return klist;
}


template <int D, typename Pmap>
Cost LBTree<D,Pmap>::fixCost(typename DClass<D>::KeyDConst& key) {
//    cout << "fixCost: key = ";
//    key.print();
//    print(" is about to be looked for");
    typename DClass<D>::treeT::iterator it = this->find(key);
//    cout << "fixCost: key = ";
//    key.print();
//    print(" was found (looked for),", (it == this->end()));
    if (it == this->end()) return 0;
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
	    d.subcost += this->template fixCost(kit.key());
	}
    }
    node.set_data(d);
//cout << "about to insert key = ";
//key.print();
//cout << ", ";
//node.get_data().print();
    this->insert(key,node);
//print("fixCost: inserted node");
    return d.subcost;
}


template <int D, typename Pmap>
Cost LBTree<D,Pmap>::depthFirstPartition(typename DClass<D>::KeyDConst& key, 
	vector<typename DClass<D>::TreeCoords>* klist, unsigned int npieces, 
	Cost totalcost, Cost *maxcost) {
//print("depthFirstPartition: at very beginning");
    if (totalcost == 0) {
	totalcost = this->template computeCost(key);
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
	usedUp = this->template makePartition(key, &tmplist, partitionSize, (i==0), usedUp, &atleaf);
	if (*maxcost < usedUp) *maxcost = usedUp;
	costLeft -= usedUp;
	partsLeft--;
	for (unsigned int j = 0; j < tmplist.size(); j++) {
	    klist->push_back(typename DClass<D>::TreeCoords(typename DClass<D>::KeyD(tmplist[j]), i)); 
	}
    }
    return totalcost;
}

template <int D, typename Pmap>
void LBTree<D,Pmap>::rollup(typename DClass<D>::KeyDConst& key) {
//    print("rollup: at beginning");
    typename DClass<D>::treeT::iterator it = this->find(key);
    if (it == this->end()) return;

//    print("rollup: about to get node associated with key",key);
    typename DClass<D>::NodeD node = it->second;
    if (!node.has_children()) {
//	print("rollup: this node has no children; returning");
	return; // no rolling to be done here.
    }
//    print("rollup: this node has children");
    bool hasleafchild = false;
    for (KeyChildIterator<D> kit(key); kit; ++kit) {
	typename DClass<D>::treeT::iterator itc = this->find(kit.key());
	if (itc != this->end()) {
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
	this->template meld(key);
    }
    for (KeyChildIterator<D> kit(key); kit; ++kit) {
	typename DClass<D>::treeT::iterator itc = this->find(kit.key());
	if (itc != this->end()) {
//	    print("rollup: found child", kit.key());
	    typename DClass<D>::NodeD c = itc->second;
	    if (c.has_children()) {
//		print("rollup: child", kit.key(), "has children");
		this->template rollup(kit.key());
	    }
	}
    }
    it = this->find(key);
    node = it->second;
    NodeData d = node.get_data();
    if (d.istaken) {
    	d.istaken = false;
    	node.set_data(d);
    	this->insert(key,node);
    }
}

template <int D, typename Pmap>
void LBTree<D,Pmap>::meld(typename DClass<D>::KeyDConst& key) {
//    print("meld: at beginning, finding key", key);
    Cost cheapest = 0;
    typename DClass<D>::treeT::iterator it = this->find(key);
    if (it == this->end()) return;

    vector<unsigned int> mylist;

    typename DClass<D>::NodeD node = it->second;
    unsigned int i = 0;
//    print("meld: about to iterate over children of key", key);
    for (KeyChildIterator<D> kit(key); kit; ++kit) {
//    	print("    meld: iterating over child", i);
	if (node.has_child(i)) {
	    typename DClass<D>::treeT::iterator itc = this->find(kit.key());
            if (itc == this->end()) return;
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
	this->insert(key,node);
	return;
    }

    NodeData d = node.get_data();

    i = 0;
    int j = 0, mlsize = mylist.size();
    for (KeyChildIterator<D> kit(key); kit; ++kit) {
	if (mylist[j] == i) {
//	    print("meld: found a match, mylist[",j,"] =",i);
	    this->erase(kit.key());
	    node.set_child(mylist[j], false);
	    d.cost += cheapest;
	    j++;
	}
	i++;
	if (j == mlsize) break;
    }
    d.istaken = false;
    node.set_data(d);
    this->insert(key,node);
//    print("meld: inserted node back into tree; goodbye!");
}


template <int D, typename Pmap>
Cost LBTree<D,Pmap>::computeCost(typename DClass<D>::KeyDConst& key) {
    Cost cost = 0;
    typename DClass<D>::treeT::iterator it = this->find(key);
    if (it == this->end()) return cost;

    typename DClass<D>::NodeD node = it->second;
    for (KeyChildIterator<D> kit(key); kit; ++kit) {
	cost += this->template computeCost(kit.key());
    }
    NodeData d = node.get_data();
    cost += d.cost;
    
    d.subcost = cost;
    node.set_data(d);
    this->insert(key,node);
    return cost;
}


template <int D, typename Pmap>
Cost LBTree<D,Pmap>::makePartition(typename DClass<D>::KeyDConst& key, 
	vector<typename DClass<D>::KeyD>* klist, Cost partitionSize, bool lastPartition, 
	Cost usedUp, bool *atleaf) {
//    cout << "at beginning of makePartition: atleaf = ";
//    cout << *atleaf << endl;
    double fudgeFactor = 0.1;
    Cost maxAddl = (Cost) (fudgeFactor*partitionSize);

    typename DClass<D>::treeT::iterator it = this->find(key);
    if (it == this->end()) { 
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
	this->template removeCost(key.parent(), d.subcost);
	node.set_data(d);
	this->insert(key,node);
    }
    else if (usedUp < partitionSize) {
	// try this node's children (if any) 
	if (node.has_children()) {
	    int i = 0;
	    for (KeyChildIterator<D> kit(key); kit; ++kit) {
		if (node.has_child(i)) {
//		    print(key, "recursively calling", kit.key());
		    usedUp = this->template makePartition(kit.key(), klist, partitionSize, lastPartition,
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

template <int D, typename Pmap>
void LBTree<D,Pmap>::removeCost(typename DClass<D>::KeyDConst& key, Cost c) {
//cout << "removeCost: key ";
//key.print();
//cout << endl;
    if (((int) key.level()) < 0) return;
    typename DClass<D>::treeT::iterator it = this->find(key);
//print("removeCost: found key");
    if (it == this->end()) return;
    typename DClass<D>::NodeD node = it->second;
    NodeData d = node.get_data();
//print("removeCost: got data");
    d.subcost -= c;
    if (key.level() > 0) {
    	this->template removeCost(key.parent(), c);
    }
//cout << "removeCost: before setting, data = ";
//d.print();
    node.set_data(d);
//cout << "removeCost: after setting, data = ";
//node.get_data().print();
    this->insert(key,node);
//cout << "removeCost: after inserting, data = ";
//node.get_data().print();
}





Cost computePartitionSize(Cost cost, unsigned int parts) {
    return (Cost) ceil(((double) cost)/((double) parts));
}


CompCost computeCompCost(Cost c, int n) {
    CompCost compcost;
    CompCost cfactor = 0.1, nfactor = 1.0;
    compcost = cfactor*c + nfactor*n;
    return compcost;
}


template <typename T, int D, typename Pmap>
void migrate_data(SharedPtr<FunctionImpl<T,D,Pmap> > tfrom, SharedPtr<FunctionImpl<T,D,Pmap> > tto, 
	typename DClass<D>::KeyD key) {
    typename FunctionImpl<T,D,Pmap>::iterator it = tfrom->find(key);
    if (it == tfrom->end()) return;

    FunctionNode<T,D> node = it->second;

    if (node.has_children()) {
        for (KeyChildIterator<D> kit(key); kit; ++kit) {
	    migrate_data<T,D,Pmap>(tfrom, tto, kit.key());
	}
    }
    tto->insert(key, node);
}


template <typename T, int D, typename Pmap>
void migrate(SharedPtr<FunctionImpl<T,D,Pmap> > tfrom, SharedPtr<FunctionImpl<T,D,Pmap> > tto) {
    typename DClass<D>::KeyD root(0);
    migrate_data<T,D,Pmap>(tfrom, tto, root);
}

// Explicit instantiations for D=1:6
template void migrate<double,3,MyProcmap<3> >(SharedPtr<FunctionImpl<double,3,MyProcmap<3> > > tfrom, 
	SharedPtr<FunctionImpl<double,3,MyProcmap<3> > > tto);

template void migrate_data<double,3>(SharedPtr<FunctionImpl<double,3,MyProcmap<3> > > tfrom, 
	SharedPtr<FunctionImpl<double,3,MyProcmap<3> > > tto, DClass<3>::KeyD key);


// Who knows why these aren't cooperating, so commented out for now
/*
template void migrate_data<double,1>(SharedPtr<FunctionImpl<double,1,MyProcmap<1> > > tfrom, 
	SharedPtr<FunctionImpl<double,1,MyProcmap<1> > > tto, DClass<1>::KeyD key);
template void migrate_data<double,2>(SharedPtr<FunctionImpl<double,2,MyProcmap<2> > > tfrom, 
	SharedPtr<FunctionImpl<double,2,MyProcmap<2> > > tto, DClass<2>::KeyD key);
template void migrate_data<double,3>(SharedPtr<FunctionImpl<double,3,MyProcmap<3> > > tfrom, 
	SharedPtr<FunctionImpl<double,3,MyProcmap<3> > > tto, DClass<3>::KeyD key);
template void migrate_data<double,4>(SharedPtr<FunctionImpl<double,4,MyProcmap<4> > > tfrom, 
	SharedPtr<FunctionImpl<double,4,MyProcmap<4> > > tto, DClass<4>::KeyD key);
template void migrate_data<double,5>(SharedPtr<FunctionImpl<double,5,MyProcmap<5> > > tfrom, 
	SharedPtr<FunctionImpl<double,5,MyProcmap<5> > > tto, DClass<5>::KeyD key);
template void migrate_data<double,6>(SharedPtr<FunctionImpl<double,6,MyProcmap<6> > > tfrom, 
	SharedPtr<FunctionImpl<double,6,MyProcmap<6> > > tto, DClass<6>::KeyD key);

template void migrate_data<std::complex<double>,1>(SharedPtr<FunctionImpl<std::complex<double>,1,MyProcmap<1> > > tfrom, 
	SharedPtr<FunctionImpl<std::complex<double>,1,MyProcmap<1> > > tto, DClass<1>::KeyD key);
template void migrate_data<std::complex<double>,2>(SharedPtr<FunctionImpl<std::complex<double>,2,MyProcmap<2> > > tfrom, 
	SharedPtr<FunctionImpl<std::complex<double>,2,MyProcmap<2> > > tto, DClass<2>::KeyD key);
template void migrate_data<std::complex<double>,3>(SharedPtr<FunctionImpl<std::complex<double>,3,MyProcmap<3> > > tfrom, 
	SharedPtr<FunctionImpl<std::complex<double>,3,MyProcmap<3> > > tto, DClass<3>::KeyD key);
template void migrate_data<std::complex<double>,4>(SharedPtr<FunctionImpl<std::complex<double>,4,MyProcmap<4> > > tfrom, 
	SharedPtr<FunctionImpl<std::complex<double>,4,MyProcmap<4> > > tto, DClass<4>::KeyD key);
template void migrate_data<std::complex<double>,5>(SharedPtr<FunctionImpl<std::complex<double>,5,MyProcmap<5> > > tfrom, 
	SharedPtr<FunctionImpl<std::complex<double>,5,MyProcmap<5> > > tto, DClass<5>::KeyD key);
template void migrate_data<std::complex<double>,6>(SharedPtr<FunctionImpl<std::complex<double>,6,MyProcmap<6> > > tfrom, 
	SharedPtr<FunctionImpl<std::complex<double>,6,MyProcmap<6> > > tto, DClass<6>::KeyD key);

template void migrate<double,1,MyProcmap<1> >(SharedPtr<FunctionImpl<double,1,MyProcmap<1> > > tfrom, 
	SharedPtr<FunctionImpl<double,1,MyProcmap<1> > > tto);
template void migrate<double,2,MyProcmap<2> >(SharedPtr<FunctionImpl<double,2,MyProcmap<2> > > tfrom, 
	SharedPtr<FunctionImpl<double,2,MyProcmap<2> > > tto);
template void migrate<double,3,MyProcmap<3> >(SharedPtr<FunctionImpl<double,3,MyProcmap<3> > > tfrom, 
	SharedPtr<FunctionImpl<double,3,MyProcmap<3> > > tto);
template void migrate<double,4,MyProcmap<4> >(SharedPtr<FunctionImpl<double,4,MyProcmap<4> > > tfrom, 
	SharedPtr<FunctionImpl<double,4,MyProcmap<4> > > tto);
template void migrate<double,5,MyProcmap<5> >(SharedPtr<FunctionImpl<double,5,MyProcmap<5> > > tfrom, 
	SharedPtr<FunctionImpl<double,5,MyProcmap<5> > > tto);
template void migrate<double,6,MyProcmap<6> >(SharedPtr<FunctionImpl<double,6,MyProcmap<6> > > tfrom, 
	SharedPtr<FunctionImpl<double,6,MyProcmap<6> > > tto);

template void migrate<std::complex<double>,1,MyProcmap<1> >(SharedPtr<FunctionImpl<std::complex<double>,1,MyProcmap<1> > tfrom, 
	SharedPtr<FunctionImpl<std::complex<double>,1,MyProcmap<1> > > tto);
template void migrate<std::complex<double>,2,MyProcmap<2> >(SharedPtr<FunctionImpl<std::complex<double>,2,MyProcmap<2> > > tfrom, 
	SharedPtr<FunctionImpl<std::complex<double>,2,MyProcmap<2> > > tto);
template void migrate<std::complex<double>,3,MyProcmap<3> >(SharedPtr<FunctionImpl<std::complex<double>,3,MyProcmap<3> > > tfrom, 
	SharedPtr<FunctionImpl<std::complex<double>,3,MyProcmap<3> > > tto);
template void migrate<std::complex<double>,4,MyProcmap<4> >(SharedPtr<FunctionImpl<std::complex<double>,4,MyProcmap<4> > > tfrom, 
	SharedPtr<FunctionImpl<std::complex<double>,4,MyProcmap<4> > > tto);
template void migrate<std::complex<double>,5,MyProcmap<5> >(SharedPtr<FunctionImpl<std::complex<double>,5,MyProcmap<5> > > tfrom, 
	SharedPtr<FunctionImpl<std::complex<double>,5,MyProcmap<5> > > tto);
template void migrate<std::complex<double>,6,MyProcmap<6> >(SharedPtr<FunctionImpl<std::complex<double>,6,MyProcmap<6> > > tfrom, 
	SharedPtr<FunctionImpl<std::complex<double>,6,MyProcmap<6> > > tto);
*/

template class LoadBalImpl<double,1,MyProcmap<1> >;
template class LoadBalImpl<double,2,MyProcmap<2> >;
template class LoadBalImpl<double,3,MyProcmap<3> >;
template class LoadBalImpl<double,4,MyProcmap<4> >;
template class LoadBalImpl<double,5,MyProcmap<5> >;
template class LoadBalImpl<double,6,MyProcmap<6> >;

template class LoadBalImpl<std::complex<double>,1,MyProcmap<1> >;
template class LoadBalImpl<std::complex<double>,2,MyProcmap<2> >;
template class LoadBalImpl<std::complex<double>,3,MyProcmap<3> >;
template class LoadBalImpl<std::complex<double>,4,MyProcmap<4> >;
template class LoadBalImpl<std::complex<double>,5,MyProcmap<5> >;
template class LoadBalImpl<std::complex<double>,6,MyProcmap<6> >;

template class LBTree<1,MyProcmap<1> >;
template class LBTree<2,MyProcmap<2> >;
template class LBTree<3,MyProcmap<3> >;
template class LBTree<4,MyProcmap<4> >;
template class LBTree<5,MyProcmap<5> >;
template class LBTree<6,MyProcmap<6> >;
}
