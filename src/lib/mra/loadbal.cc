/*
  This file is part of MADNESS.
  
  Copyright (C) <2007> <Oak Ridge National Laboratory>
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
  
  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov 
  tel:   865-241-3937
  fax:   865-572-0680

  
  $Id$
*/

/// \file loadbal.cc
/// \brief Implements class methods associated with load balancing.
  
#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>


namespace madness {
    /// find_best_partition takes the result of find_partitions, determines
    /// which is the best partition, and broadcasts that to all processors
    template <typename T, int D>
    std::vector<typename DClass<D>::TreeCoords> LoadBalImpl<T,D>::find_best_partition() {
	std::vector<typename DClass<D>::TreeCoords> klist;
//	typename DClass<D>::KeyDConst root(0);
	bool manager = false;

	if (skeltree->world.mpi.nproc() == 1) {
	    klist.push_back(typename DClass<D>::TreeCoords(skeltree->root, 0));
	    return klist;
	}

	if (skeltree->world.mpi.rank() == skeltree->owner(skeltree->root)) manager = true;
	madness::print("find_best_partition: just starting out");
//	skeltree->fix_cost();
	madness::print("find_best_partition: fixed cost");
	std::vector< std::vector<typename DClass<D>::TreeCoords> > list_of_list;
	std::vector<Cost> costlist;
	skeltree->find_partitions(this->pi, list_of_list, costlist);
	madness::print("find_best_partition: before fence right after being out of find_partitions");
	skeltree->world.gop.fence();
	madness::print("find_best_partition: just finished find_partitions");


	if (manager) {
 	    madness::print("find_best_partition: I am the manager");
 	    unsigned int shortest_list = 0, sl_index, lb_index;
 	    Cost load_bal_cost = 0;
 	    int count = list_of_list.size();
 	    madness::print("find_best_partition: length of list_of_list =", count);
 	    std::vector<unsigned int> len;
 	    for (int i = 0; i < count; i++) {
 		len.push_back(list_of_list[i].size());
 		if ((len[i] < shortest_list) || (shortest_list == 0)) {
 		    shortest_list = len[i];
 		    sl_index = i;
 		} else if ((len[i] == shortest_list) && (costlist[i] < costlist[sl_index])) {
 		    // all things being equal, prefer better balance
 		    shortest_list = len[i];
 		    sl_index = i;
 		}
 		if ((costlist[i] < load_bal_cost) || (load_bal_cost == 0)) {
 		    load_bal_cost = costlist[i];
 		    lb_index = i;
 		} else if ((costlist[i] == load_bal_cost) && (len[i] < list_of_list[lb_index].size())) {
 		    // all things being equal, prefer fewer cuts
 		    load_bal_cost = costlist[i];
 		    lb_index = i;
 		}
 	    }
	    
 	    CompCost ccleast = 0;
 	    int cc_index = 0;
 	    for (int i = 0; i < count; i++) {
 		CompCost cctmp = compute_comp_cost(costlist[i], len[i]-1);
 		if ((i==0) || (cctmp < ccleast)) {
 		    ccleast = cctmp;
 		    cc_index = i;
 		}
 	    }
	    
 	    madness::print("The load balance with the fewest broken links has cost", 
 			   costlist[sl_index], "and", shortest_list-1, "broken links");
 	    for (unsigned int i = 0; i < shortest_list; i++) {
 		list_of_list[sl_index][i].print();
 	    }
 	    madness::print("");
 	    madness::print("The load balance with the best balance has cost", 
 			   load_bal_cost, "and", list_of_list[lb_index].size()-1, 
 			   "broken links");
 	    for (unsigned int i = 0; i < list_of_list[lb_index].size(); i++) {
 		list_of_list[lb_index][i].print();
 	    }
 	    madness::print("");
 	    madness::print("The load balance with the best overall computational cost has cost",
 			   costlist[cc_index], "and", len[cc_index]-1, "broken links");
 	    for (unsigned int i = 0; i < len[cc_index]; i++) {
 		list_of_list[cc_index][i].print();
 	    }
	    
 	    for (unsigned int i = 0; i < len[cc_index]; i++) {
 		klist.push_back(list_of_list[cc_index][i]);
 	    }
 	    unsigned int ksize;
 	    ksize = klist.size();
 	    skeltree->world.gop.template broadcast<unsigned int>(ksize);
 	    for (unsigned int i=0; i < ksize; i++) {
 		skeltree->world.gop.template broadcast<typename DClass<D>::TreeCoords>(klist[i]);
 	    }
	    madness::print("find_best_partition: number of broken links =",
		klist.size()-1);
	}
	else {
 		madness::print("find_best_partition: receiving broadcast");
 		std::vector<typename DClass<D>::TreeCoords> klist;
 		unsigned int ksize;
 		skeltree->world.gop.template broadcast<unsigned int>(ksize);
 		for (unsigned int i=0; i < ksize; i++) {
 		    skeltree->world.gop.template broadcast<typename DClass<D>::TreeCoords>(klist[i]);
 		}
	}
        return klist;	
    };



    /// Compute the cost of a given configuration: a weighted sum of the cost of the
    /// maximally-loaded process and the total number of broken links.
    /// In the future, the factors should be calibrated for a given machine, either 
    /// during configuration and setup or in real time at the beginning of the program
    /// or upon the creation of a LoadBalImpl.
    /// Arguments: Cost c -- maximum amount of cost assigned to a node
    ///            int n -- number of broken links
    /// Return: CompCost -- the cost associated with that partitioning
    template <typename T, int D>
    CompCost LoadBalImpl<T,D>::compute_comp_cost(Cost c, int n) {
        CompCost compcost;
	int k = f.k();
	CompCost k_to_D = pow((CompCost) k,D);
	CompCost twok_to_Dp1 = pow((CompCost) 2.0*k, D+1);
	compcost = c*(flop_time*D*twok_to_Dp1) + n*(comm_bandw*k_to_D + comm_latency);
        return compcost;
    }



    /// find_partitions performs the "melding" algorithm for load balancing: it recursively melds 
    /// and partitions the tree until it has found all possible configurations.
    template <int D> 
    void LBTree<D>::find_partitions(PartitionInfo<D>& lbi, std::vector< std::vector<typename DClass<D>::TreeCoords> >& list_of_list,
				    std::vector<Cost>& costlist) {
//	typename DClass<D>::KeyDConst root(0);
	bool manager = false, keep_going = true, first_time = true;
	madness::print("find_partitions: at beginning");
	this->world.gop.fence();

	if (this->world.mpi.rank() == this->impl.owner(root)) manager = true;
	
// 	madness::print("find_partitions: about to test_Void");
// 	if (manager) test_Void(root);
// 	this->world.gop.fence();
// 	madness::print("find_partitions: done with test_Void");

	while (keep_going) {
	    madness::print("find_partitions: the verdict is to keep_going");
	    meld_all(first_time);
	    madness::print("find_partitions: after meld_all");
	    if (manager) {
		madness::print("find_partitions: about to launch_make_partition");
		keep_going = this->launch_make_partition(lbi, list_of_list, costlist, first_time);
		madness::print("find_partitions: returned from launch_make_partition");
		madness::print("find_partitions: lbi =", lbi);
	    }
	    first_time = false;
	    this->world.gop.fence();
	    this->world.gop.template broadcast<bool>(keep_going);
	}
    }


    template <int D>
    Void LBTree<D>::test_Void(typename DClass<D>::KeyDConst& key) {
	madness::print("test_Void: at beginning, key =", key);
       	typename DClass<D>::treeT::iterator it = impl.find(key);
	if (it == impl.end()) return None;
	typename DClass<D>::NodeD node = it->second;
	if (node.has_children()) {
	    typename DClass<D>::KeyD child = first_child(key, node);
	    madness::print("test_Void: about to call for child", child, "owned by", owner(child));
	    send(impl.owner(child), &LBTree::test_Void, child);
	    madness::print("test_Void:", key, "right back out of there");
	} else {
	    typename DClass<D>::KeyD sib = next_sibling(key);
	    if (sib == typename DClass<D>::KeyD(0)) {
		madness::print("test_Void: back to root; all done");
		return None;
	    }
	    madness::print("test_Void: about to call for sib", sib, "owned by", owner(sib));
	    send(impl.owner(sib), &LBTree::test_Void, sib);
	    madness::print("test_Void:", key, "right back out of there");
	}
	madness::print("test_Void: key", key, "about to return None");
	return None;
    }



    template <int D>
    bool LBTree<D>::launch_make_partition(PartitionInfo<D>& lbi, 
					  std::vector< std::vector<typename DClass<D>::TreeCoords> >& list_of_list, 
					  std::vector<Cost>& costlist, bool first_time = false) {
	madness::print("launch_make_partition: at beginning");
	bool keep_going = true;
//	typename DClass<D>::KeyDConst root(0);
//	std::vector< std::vector<typename DClass<D>::TreeCoords> > list_of_list(1, std::vector<typename DClass<D>::TreeCoords>());
	unsigned int npieces;
	Cost used_up;
	// just manager process	    
	madness::print("launch_make_partition: it's good to be the manager");
	
	if (!first_time) {
//	    list_of_list.push_back(std::vector<typename DClass<D>::TreeCoords>());
//	    lbi.costlist.push_back(0);
	    lbi.step_num++;
	}
	
	npieces = world.nproc();
	lbi.partition_number = npieces-1;
	Cost tpart = compute_partition_size(lbi.skel_cost, npieces);
	used_up = 0;
	madness::print("launch_make_partition: lbi =", lbi);
	send(impl.owner(root), &LBTree<D>::make_partition, root, tpart, used_up, lbi, true);
	madness::print("launch_make_partition: back from send make_partition");
	unsigned int secret_signal;
	int tag = 3675384;
	int source;
	int bslen = (sizeof(int)+MPI::BSEND_OVERHEAD)*2;
	void *bshort = malloc(bslen);
	world.mpi.Attach_buffer(bshort, bslen);
	
	MPI::Request req = world.mpi.Irecv<int>(source, MPI::ANY_SOURCE, tag);
	tag++;

	madness::print("launch_make_partition: about to await request");
	world.await(req);
	madness::print("launch_make_partition: request finally granted, source =", source);
	world.mpi.Detach_buffer((void *&) bshort);
	free(bshort);

	bool recv_from_self = (source == world.mpi.rank());
	int bs2len = (sizeof(unsigned int)+MPI::BSEND_OVERHEAD)*2;

	void *bshort2 = malloc(bs2len);
	world.mpi.Attach_buffer(bshort2, bs2len);
	world.mpi.Recv<unsigned int>(secret_signal, source, tag);
	world.mpi.Detach_buffer((void *&) bshort2);
	free(bshort2);

	tag++;
	if (recv_from_self) {
	    int bsize = secret_signal*sizeof(typename DClass<D>::TreeCoords) + 3*sizeof(Cost) + 5*sizeof(unsigned int) + MPI::BSEND_OVERHEAD*3;
	    void *buffer = malloc(bsize);
	    world.mpi.Attach_buffer(buffer, bsize);
	    typename DClass<D>::TreeCoords tc;
	    std::vector<typename DClass<D>::TreeCoords> tmplist;
	    for (unsigned int i = 0; i < secret_signal; i++) {
		world.mpi.Recv<typename DClass<D>::TreeCoords>(tc, source, tag);
		madness::print("launch_make_partition: received part_list[",i,"]");
		tag++;
		tmplist.push_back(tc);
	    }
	    madness::print("launch_make_partition: tmplist =", tmplist);
	    Cost tmpcost;
	    world.mpi.Recv<Cost>(tmpcost, source, tag);
	    tag++;
	    madness::print("launch_make_partition: received maxcost");
	    costlist.push_back(tmpcost);
	    world.mpi.Recv<Cost>(lbi.cost_left, source, tag);
	    tag++;
	    madness::print("launch_make_partition: received cost_left");
	    world.mpi.Recv<Cost>(lbi.skel_cost, source, tag);
	    tag++;
	    madness::print("launch_make_partition: received skel_cost");
	    world.mpi.Recv<unsigned int>(lbi.partition_number, source, tag);
	    tag++;
	    madness::print("launch_make_partition: received partition_number");
	    world.mpi.Recv<unsigned int>(lbi.step_num, source, tag);
	    madness::print("launch_make_partition: received step_num");
	    list_of_list.push_back(tmplist);
	    madness::print("launch_make_partition: received lbi:", lbi);
	    free(buffer);
	} else {
	    typename DClass<D>::TreeCoords tc;
	    std::vector<typename DClass<D>::TreeCoords> tmplist;
	    for (unsigned int i = 0; i < secret_signal; i++) {
		world.mpi.Recv<typename DClass<D>::TreeCoords>(tc, source, tag);
		madness::print("launch_make_partition: received part_list[",i,"]");
		tag++;
		tmplist.push_back(tc);
	    }
	    madness::print("launch_make_partition: tmplist =", tmplist);
	    Cost tmpcost;
	    world.mpi.Recv<Cost>(tmpcost, source, tag);
	    tag++;
	    madness::print("launch_make_partition: received maxcost");
	    costlist.push_back(tmpcost);
	    world.mpi.Recv<Cost>(lbi.cost_left, source, tag);
	    tag++;
	    madness::print("launch_make_partition: received cost_left");
	    world.mpi.Recv<Cost>(lbi.skel_cost, source, tag);
	    tag++;
	    madness::print("launch_make_partition: received skel_cost");
	    world.mpi.Recv<unsigned int>(lbi.partition_number, source, tag);
	    tag++;
	    madness::print("launch_make_partition: received partition_number");
	    world.mpi.Recv<unsigned int>(lbi.step_num, source, tag);
	    madness::print("launch_make_partition: received step_num");
	    list_of_list.push_back(tmplist);
	    madness::print("launch_make_partition: received lbi:", lbi);
	}
/*
	if (recv_from_self) {
	    tag++;
	    world.mpi.Recv<void*>(buffer, source, tag);
	    world.mpi.Detach_buffer((void *&) buffer);
	    delete[] buffer;
	}
*/
	

	if (lbi.partition_number == 0) {
	    // make sure current partition is valid.  If not, quit.
	    // first, add root to partition
	    
	    madness::print("launch_make_partition: add root to partition and be done");
//	    lbi.list_of_list[lbi.step_num].push_back(typename DClass<D>::TreeCoords(typename DClass<D>::KeyD(0), lbi.partition_number));
	    madness::print("launch_make_partition: size of list_of_list[", lbi.step_num, "] =", list_of_list[lbi.step_num].size());
	    int count = lbi.step_num;
	    int lolcsize = list_of_list[count].size();
	    int npieces = this->world.nproc()-1;
	    madness::print("launch_make_partition: after resetting some stuff, lbi =", lbi);
	    
	    // Making sure we don't have an invalid partition, e.g. a case where one (or
	    // more!) processor(s) do(es)n't have any work.
	    int m = npieces;
	    bool invalid_partition = false;
	    for (int k = 0; k < lolcsize; k++) {
		int difff = m-list_of_list[count][k].owner;
		if (difff == 1) {
		    m--;
		} else if (difff > 1) {
		    invalid_partition = true;
		    break;
		}
	    }
	    if ((lolcsize < npieces) || (invalid_partition)) {
		list_of_list.erase(list_of_list.begin()+count);
		keep_going = false;
		count-=1;
	    }
	    count++;
	    madness::print("launch_make_partition: the verdict is that keep_going =", keep_going);
	
	}
	return keep_going;
    }


    template <int D>
    Void LBTree<D>::meld_all(bool first_time = false) {
	madness::print("meld_all: after launching on everybody, about to go to fix_cost");
	this->world.gop.barrier();
	madness::print("meld_all: about to fix_cost");
	if (!first_time) {
	    this->fix_cost();
	}
	madness::print("meld_all: done with fix_cost");
	this->world.gop.barrier();
	madness::print("meld_all: about to reset");
	this->reset(true);
	this->world.gop.barrier();
	madness::print("meld_all: about to rollup");
	this->rollup();
	this->world.gop.barrier();
	madness::print("meld_all: about to reset");
	this->reset(false);
	madness::print("meld_all: done with reset");
//	if (me == 0) this->world.gop.fence();
	this->world.gop.barrier();
	madness::print("meld_all: done with barrier; returning None");
	return None;
    }



    /// fix_cost resets the tree after the load balancing and melding have been performed, before the next
    /// round of load balancing.
    /// Arguments: none
    /// Return: Cost of entire tree
    /// Side effect: subcost (the cost of the subtree rooted at key) is reset
    /// Communication: none except communication by methods it calls


    template <int D>
    Cost LBTree<D>::fix_cost() {
	init_fix_cost();
	madness::print("fix_cost: done with init_fix_cost");
	this->world.gop.fence();
//	this->world.gop.barrier();
	madness::print("fix_cost: about to fix_cost_spawn");
	fix_cost_spawn();
	madness::print("fix_cost: done with fix_cost_spawn, now for fence");
	this->world.gop.fence();
//	this->world.gop.barrier();
	madness::print("AFTER FIXING COST");
	for (typename DClass<D>::treeT::iterator it = impl.begin(); it != impl.end(); ++it) {
	    madness::print(it->first, it->second);
	}
	madness::print("DONE FIXING IT");
	if (this->world.mpi.rank() == impl.owner(typename DClass<D>::KeyD(0))) {
	    typename DClass<D>::treeT::iterator it = impl.find(typename DClass<D>::KeyD(0));
	    return it->second.get_data().subcost;
	} else {
	    return 0;
	}
    }


    /// init_fix_cost resets and zeroes out nrecvd in each node
    /// Arguments: none
    /// Return: none
    /// Side effects: none
    /// Communication: none

    template <int D>
    void LBTree<D>::init_fix_cost() {
	for (typename DClass<D>::treeT::iterator it = impl.begin(); it != impl.end(); ++it) {
	    typename DClass<D>::KeyDConst& key = it->first;
	    typename DClass<D>::NodeD& node = it->second;

	    int dim = node.dim;
	    NodeData d = node.get_data();
	    d.subcost = d.cost;
	    node.nrecvd = dim - node.get_num_children();
	    node.set_data(d);
	    impl.insert(key,node);
	}
    }


    /// fix_cost_spawn launches sum up tree, beginning at leaf nodes
    /// Arguments: none
    /// Return: none
    /// Side effects: fixes subcost on each node
    /// Communication: sends cost of each leaf node to owner of its parent

    template <int D>
    void LBTree<D>::fix_cost_spawn() {
	for (typename DClass<D>::treeT::iterator it = impl.begin(); it != impl.end(); ++it) {
	    typename DClass<D>::KeyDConst& key = it->first;
	    typename DClass<D>::NodeD& node = it->second;
	    if (!node.has_children()) {
		typename DClass<D>::KeyD parent = key.parent();
		Cost c = node.get_data().cost;
//		madness::print("fix_cost_spawn: key", key, "is leaf child; sending", c,
//	       		       "to parent", parent, "at processor", impl.owner(parent));
		send(impl.owner(parent), &LBTree<D>::fix_cost_sum, parent, c);
	    }
	}
    }


    /// fix_cost_sum receives node cost from child, adds to this node's subcost,
    /// and, if it's added in all the costs from below, sends its cost to its parent
    /// Arguments: const Key<D> key: the key of the node in question
    ///            Cost c: the cost of the subtree rooted by the child
    /// Return: none
    /// Side effects: none
    /// Communication: may send its subtree cost to parent

    template <int D>
    Void LBTree<D>::fix_cost_sum(typename DClass<D>::KeyDConst& key, Cost c) {
	typename DClass<D>::treeT::iterator it = impl.find(key);
	typename DClass<D>::NodeD node = it->second;
	NodeData d = node.get_data();
	d.subcost += c;
//	madness::print("fix_cost_sum:", key, "received number", node.nrecvd+1, "cost", 
//		       c, " subtotal =", d.subcost);
	node.nrecvd++;
	node.set_data(d);
	impl.insert(key, node);
	if ((node.nrecvd == node.dim) && (key.level()!=0)) {
	    typename DClass<D>::KeyD parent = key.parent();
//	    madness::print("fix_cost_sum:", key, "sending cost", d.subcost, "to parent", parent);
	    task(impl.owner(parent), &LBTree<D>::fix_cost_sum, parent, d.subcost);
	}
	return None;
	
    }


    /// depth_first_partition finds partitions of trees by calling make_partition.  It figures
    /// out the size of the partitions to be made by make_partition at each iteration.
    /// Arguments: PartitionInfo<D> lbi -- contains partitioning info, including skel_cost total cost of
    ///                skeleton tree and other useful accounting-type variables
    /// Return: none
    /// Side effect: lbi updated by make_partition

/*
    template <int D> 
    void LBTree<D>::depth_first_partition(PartitionInfo<D> lbi) {
	typename DClass<D>::KeyDConst root(0);
	if (this->world.mpi.rank() == impl.owner(root)) {
	    Cost total = lbi.skel_cost;
	    unsigned int npieces = world.nproc();
	    Cost tpart = compute_partition_size(total, npieces);
	    Cost used_up = 0;
	    send(impl.owner(root), &LBTree<D>::make_partition, root, tpart, used_up, lbi);
//	    make_partition(root, tpart, used_up, lbi);
	}
    }
*/   

    /// rollup traverses the tree, calling meld upon Nodes that have leaf children
    /// Arguments: const Key<D> key -- node at which we begin
    /// Side effect: Nodes are changed by meld
    /// Communication: just finding the nodes that match a given key
    template <int D>
    void LBTree<D>::rollup() {
	for (typename DClass<D>::treeT::iterator it = impl.begin(); it != impl.end(); ++it) {
	    typename DClass<D>::KeyD key = it->first;
	    typename DClass<D>::NodeD node = it->second;
	    if (node.has_children()) {
                // First, check to see if it has any leaf children
	        bool has_leaf_child = false;
	        for (KeyChildIterator<D> kit(key); kit; ++kit) {
            	    typename DClass<D>::treeT::iterator itc = impl.find(kit.key());
		    if (itc != impl.end()) {
		        typename DClass<D>::NodeD c = itc->second;
		        NodeData d = c.get_data();
		        if ((!c.has_children()) && (d.is_taken)) {
			    has_leaf_child = true;
			    break;
		        }
		    }
	        }
	        if (has_leaf_child) {
                    // If there is at least one leaf child, then this node gets melded.
		    this->meld(it);
		    node = it->second;
	        }
	        NodeData d = node.get_data();
	        if (d.is_taken) {
                    // Setting to false, to signify that this node has been worked on.
		    d.is_taken = false;
		    node.set_data(d);
		    impl.insert(key,node);
	        }
	    }
	}
    }

    /// reset sets the is_taken variable within all local nodes to the value specified
    /// Arguments: bool taken -- value to set is_taken to
    /// Communication: none (local iterator)
    template <int D>
    void LBTree<D>::reset(bool taken) {
	for (typename DClass<D>::treeT::iterator it = impl.begin(); it != impl.end(); ++it) {
	    typename DClass<D>::KeyD key = it->first;
	    typename DClass<D>::NodeD node = it->second;
	    NodeData d = node.get_data();

	    d.is_taken = taken;
	    node.set_data(d);
	    impl.insert(key,node);
	}
    }

    /// meld fuses leaf child(ren) to parent and deletes the leaf child(ren) in question
    /// Arguments: const Key<D> key -- node at which we begin
    /// Side effect: parent nodes are updated, and leaf nodes are deleted
    /// Communication: find and insert 
    template <int D>
    void LBTree<D>::meld(typename DClass<D>::treeT::iterator it) {
	typename DClass<D>::KeyD key = it->first;
	typename DClass<D>::NodeD node = it->second;
	std::vector<unsigned int> mylist;
	unsigned int i = 0;
	Cost cheapest = 0;
	bool not_yet_found = true;

	for (KeyChildIterator<D> kit(key); kit; ++kit) {
	    if (node.has_child(i)) {
		typename DClass<D>::treeT::iterator itc = impl.find(kit.key());
		typename DClass<D>::NodeD c = itc->second;
		NodeData d = c.get_data();
		// if the child has no children and the is_taken flag is set to true, 
		// then this child is eligible to be melded into parent
		if ((!c.has_children()) && (d.is_taken)) {
		    Cost cost = d.cost;
		    if ((cost < cheapest) || (not_yet_found)) {
                        // if this child has the cheapest cost yet, then clear the 
			// list and addthis child to the list of children to be 
			// melded to the parent
			not_yet_found = false;
			cheapest = cost;
			mylist.clear();
			mylist.push_back(i);
		    } else if (cost == cheapest) {
                        // if this child's cost is equal to the cheapest cost found 
			// so far, then add this child to the list of children to be 
			// melded into parent
			mylist.push_back(i);
		    }
		}
	    }
	    i++;
	}
	if (not_yet_found) {
	    // this node has no leaf children
	    return;
	}
        // Now we do the actual melding
	NodeData d = node.get_data();
	i = 0;
	int j = 0, mlsize = mylist.size();
	for (KeyChildIterator<D> kit(key); kit; ++kit) {
	    if (mylist[j] == i) {
		impl.erase(kit.key());
		node.set_child(mylist[j], false);
		d.cost += cheapest;
		j++;
		if (j == mlsize) break;
	    }
	    i++;
	}
	node.set_data(d);
	impl.insert(key, node);
    }


    /// compute_cost resets the subtree cost value for a tree and its descendants.
    /// Arguments: const Key<D> key -- node at which we begin
    /// Return: Cost -- the Cost of the subtree headed by key
    /// Side effect: subcost is updated
    /// Communication: find and insert
/*
    template <int D>
    Cost LBTree<D>::compute_cost(typename DClass<D>::KeyDConst& key) {
        Cost cost = 0;
        typename DClass<D>::treeT::iterator it = impl.find(key);
        if (it == impl.end()) return cost;

        typename DClass<D>::NodeD node = it->second;
        for (KeyChildIterator<D> kit(key); kit; ++kit) {
            cost += this->compute_cost(kit.key());
        }
        NodeData d = node.get_data();
        cost += d.cost;

        d.subcost = cost;
        node.set_data(d);
        impl.insert(key,node);
        return cost;
    }
*/

    /// make_partition creates a partition.  It's called by depth_first_partition to actually do all the dirty 
    /// work for each partition. 
    /// Arguments: const Key<D> key -- node at which we begin
    ///            vector<Key<D> >* klist -- list of subtree root nodes obtained from partitioning
    ///            Cost partition_size -- the target size for the partition
    ///            bool last_partition -- is this the final partition
    ///            Cost used_up -- the cost used up so far in this partition
    ///            bool *at_leaf -- are we at a leaf node
    /// Return: Cost -- the Cost of what was used up in this partition 
    /// Side effect: klist and at_leaf are updated
    /// Communication: find and insert


    template <int D>
    Void LBTree<D>::make_partition(typename DClass<D>::KeyDConst& key, Cost partition_size, Cost used_up, PartitionInfo<D> lbi, bool downward = false) {

        // The fudge factor is the fraction by which you are willing to let the
        // partitions exceed the ideal partition size
        double fudge_factor = 0.1;
        Cost maxAddl = (Cost) (fudge_factor*partition_size);

	madness::print("make_partition: at beginning for key =", key);
	madness::print("make_partition: at beginning, lbi =", lbi);

	typename DClass<D>::treeT::iterator it = impl.find(key);
	if (it == impl.end()) return None;
	typename DClass<D>::NodeD node = it->second;

//	typename DClass<D>::KeyDConst root(0);
	typename DClass<D>::KeyDConst parent = key.parent();
	NodeData d = node.get_data();
	madness::print("make_partition: just found", key, node);


        // if either the partition is currently empty and this is a single item, or there is still 
	// room in the partition and adding this to it won't go above the fudge factor, then add this 
	// piece to the partition.
        if ((downward) && (((used_up == 0) && (!node.has_children())) ||
                ((used_up < partition_size) && (d.subcost+used_up <= partition_size+maxAddl)))) {
            // add to partition
	    madness::print("make_partition: adding", key, "to partition");
	    used_up += d.subcost;
//            send(impl.owner(key.parent()), &LBTree<D>::remove_cost, key.parent(), d.subcost);
	    lbi.part_list.push_back(TreeCoords<D>(key, lbi.partition_number));
//	    send(0, &LBTree<D>::add_to_partition, key, lbi);
	    madness::print("make_partition: size of part_list =", lbi.part_list.size());
	    if (key == root) {
		madness::print("make_partition: OMG root has no children and was added to partition!!!");
		totally_reset(lbi);
	    } else {
		send(impl.owner(parent), &LBTree::make_partition, parent, partition_size, used_up, lbi);
	    }
	    return None;
	}

	if (node.has_children()) {
	    madness::print("make_partition:", key, "too big for partition");
	    if (downward) {
		madness::print("make_partition: onward, downward");
		// ADJUST THIS TO ++ until an existing child is found or to end if no more children!!!!
		node.rpit = KeyChildIterator<D>(key);
		madness::print("make_partition: set rpit to", node.rpit.key(), node.rpit);
		impl.insert(key,node);
	    } else {
		// ADJUST THIS TO ++ until an existing child is found or to end if no more children!!!!
		++(node.rpit);
		madness::print("make_partition: incremented rpit to", node.rpit.key(), node.rpit);
		impl.insert(key,node);
	    }
	    if (node.rpit) {
		typename DClass<D>::KeyDConst& child = node.rpit.key();
		madness::print("make_partition: passing the torch to", child, "owned by", impl.owner(child));
		send(impl.owner(child), &LBTree::make_partition, child, partition_size, used_up, lbi, true);
		impl.insert(key,node);
		return None;
	    }
	}

	// All done with children, or I am a leaf
        if (((used_up == 0) && (!node.has_children())) ||
                ((used_up < partition_size) && (d.cost+used_up <= partition_size+maxAddl))) {
            // add to partition
	    madness::print("make_partition: adding", key, "to partition");
	    used_up += d.cost;
	    lbi.part_list.push_back(TreeCoords<D>(key, lbi.partition_number));
	    madness::print("make_partition: size of part_list =", lbi.part_list.size());
	    if (key == root) {
		madness::print("make_partition: added root at the end, now let's quit");
	    } else {
		madness::print("make_partition: now sending back to parent", parent);
		send(impl.owner(parent), &LBTree::make_partition, parent, partition_size, used_up, lbi);
	    }
	} else {
	    madness::print("make_partition: uh oh!  Need to reset partition!");
	    bool continue_as_normal = reset_partition(key, partition_size, used_up, lbi);
	    if (continue_as_normal) {
		send(impl.owner(key), &LBTree::make_partition, key, partition_size, used_up, lbi, downward);
	    } else {
		// totally reset
		madness::print("make_partition: time to totally reset");
		totally_reset(lbi);
	    }
	}
	return None;
    }

    template <int D>
    void LBTree<D>::totally_reset(PartitionInfo<D> lbi) {
	bool send_to_self = (world.mpi.rank() == impl.owner(root));
	int tag = 3675384;	
	unsigned int listsize = lbi.part_list.size();
	if (send_to_self) {
	    int bsize = listsize*sizeof(typename DClass<D>::TreeCoords) + 3*sizeof(Cost) + 5*sizeof(unsigned int) + MPI::BSEND_OVERHEAD*3;
	    void *buffer = malloc(bsize);
	    madness::print("totally_reset: about to attach buffer", buffer);
	    world.mpi.Attach_buffer(buffer, bsize);
	    madness::print("totally_reset: just finished attaching buffer");
	    world.mpi.Ibsend<int>(world.mpi.rank(), impl.owner(root), tag);
	    tag++;
	    world.mpi.Ibsend<unsigned int>(listsize, impl.owner(root), tag);
	    tag++;
	    madness::print("totally_reset: sent part_list_size =", listsize);
	    for (unsigned int i = 0; i < listsize; i++) {
		world.mpi.Ibsend<typename DClass<D>::TreeCoords>(lbi.part_list[i], impl.owner(root), tag);
		tag++;
	    }
	    madness::print("totally_reset: sent part_list =", lbi.part_list);
	    world.mpi.Ibsend<Cost>(lbi.maxcost, impl.owner(root), tag);
	    tag++;
	    madness::print("totally_reset: sent maxcost =", lbi.maxcost);
	    world.mpi.Ibsend<Cost>(lbi.cost_left, impl.owner(root), tag);
	    tag++;
	    madness::print("totally_reset: sent cost_left =", lbi.cost_left);
	    world.mpi.Ibsend<Cost>(lbi.skel_cost, impl.owner(root), tag);
	    tag++;
	    madness::print("totally_reset: sent skel_cost =", lbi.skel_cost);
	    world.mpi.Ibsend<unsigned int>(lbi.partition_number, impl.owner(root), tag);
	    tag++;
	    madness::print("totally_reset: sent partition_number =", lbi.partition_number);
	    world.mpi.Ibsend<unsigned int>(lbi.step_num, impl.owner(root), tag);
	    madness::print("totally_reset: sent step_num =", lbi.step_num);
	    world.mpi.Detach_buffer((void *&) buffer);
	    free(buffer);
	    madness::print("totally_reset: sent", impl.owner(root), "lbi: ", lbi);
	} else {
	    world.mpi.Send<int>(world.mpi.rank(), impl.owner(root), tag);
	    tag++;
	    world.mpi.Send<unsigned int>(listsize, impl.owner(root), tag);
	    tag++;
	    madness::print("totally_reset: sent part_list_size =", listsize);
	    for (unsigned int i = 0; i < listsize; i++) {
		world.mpi.Send<typename DClass<D>::TreeCoords>(lbi.part_list[i], impl.owner(root), tag);
		tag++;
	    }
	    madness::print("totally_reset: sent part_list =", lbi.part_list);
	    world.mpi.Send<Cost>(lbi.maxcost, impl.owner(root), tag);
	    tag++;
	    madness::print("totally_reset: sent maxcost =", lbi.maxcost);
	    world.mpi.Send<Cost>(lbi.cost_left, impl.owner(root), tag);
	    tag++;
	    madness::print("totally_reset: sent cost_left =", lbi.cost_left);
	    world.mpi.Send<Cost>(lbi.skel_cost, impl.owner(root), tag);
	    tag++;
	    madness::print("totally_reset: sent skel_cost =", lbi.skel_cost);
	    world.mpi.Send<unsigned int>(lbi.partition_number, impl.owner(root), tag);
	    tag++;
	    madness::print("totally_reset: sent partition_number =", lbi.partition_number);
	    world.mpi.Send<unsigned int>(lbi.step_num, impl.owner(root), tag);
	    madness::print("totally_reset: sent step_num =", lbi.step_num);
	    madness::print("totally_reset: sent", impl.owner(root), "lbi: ", lbi);
	}	    
    }

    template <int D>
    typename DClass<D>::KeyD LBTree<D>::first_child(typename DClass<D>::KeyDConst& key, const typename DClass<D>::NodeD& node) {
	for (int i = 0; i < 1<<D; i++) {
	    if (node.has_child(i)) {
		int j = 0;
		for (KeyChildIterator<D> kit(key); kit; ++kit, ++j) {
		    if (j == i) {
			return kit.key();
		    }
		}
	    }
	}
	return key;
    }


    template <int D>
    typename DClass<D>::KeyD LBTree<D>::next_sibling(typename DClass<D>::KeyDConst& key) {
//	madness::print("next_sibling: at beginning");
	if (key.level() == 0) return key;
	typename DClass<D>::KeyD parent = key.parent();
//	madness::print("next_sibling: computed parent", parent);
//	typename DClass<D>::NodeD node = send(impl.owner(parent), &LBTree<D>::gimme_node, parent);
//	madness::print("next_sibling: got back from gimme_node");
//	typename DClass<D>::treeT::iterator it = impl.find(parent);
//	madness::print("next_sibling: just did find of parent");
//	if (it == impl.end()) { 
//	    madness::print("next_sibling: could not find parent", parent, "returning", key);
//	    return key;
//	}
//	typename DClass<D>::NodeD node = it->second;
	int j = 0, moi = 1<<(D+1);
	for (KeyChildIterator<D> kit(parent); kit; ++kit, ++j) {
//	    madness::print("next_sibling: looking at", kit.key());
	    if ((kit.key()) == key) {
		moi = j;
	    }
//	    if ((moi < j) && (node.has_child(j)))
	    if (moi < j)
		return kit.key();
	}
	if (key.level() != 0) {
//	    madness::print("next_sibling: no existing siblings at this level, look for sibling of parent");
	    return next_sibling(parent);
	}
	else {
	    return key;
	}
    }


/*
    template <int D>
    Void LBTree<D>::add_to_partition(typename DClass<D>::KeyDConst& key, PartitionInfo<D> lbi) {
	madness::print("add_to_partition: about to push_back", key, "Proc:", lbi.partition_number);
	lbi.list_of_list[lbi.step_num].push_back(typename DClass<D>::TreeCoords(key, lbi.partition_number));
	madness::print("add_to_partition: partition number", lbi.step_num, "now has the following:");
	for (unsigned int i = 0; i < lbi.list_of_list[lbi.step_num].size(); i++) {
	    madness::print(lbi.list_of_list[lbi.step_num][i]);
	}
	return None;
    }
*/

    template <int D>
    bool LBTree<D>::reset_partition(typename DClass<D>::KeyDConst& key, Cost& partition_size, Cost& used_up, PartitionInfo<D>& lbi) {
	lbi.partition_number--;
	if (lbi.partition_number > 0) {
	    if (used_up > lbi.maxcost) {
		lbi.maxcost = used_up;
	    }
	    lbi.cost_left -= used_up;
	    used_up = 0;
	    Cost tpart = compute_partition_size(lbi.cost_left, lbi.partition_number+1);
	    if ((tpart > partition_size) || (tpart*lbi.facter < partition_size)) {
		partition_size = tpart;
	    }
	    return true;
	} else {
	    lbi.part_list.push_back(typename DClass<D>::TreeCoords(root, lbi.partition_number));
	    // totally reset!
	    return false;
	}
    }

/*
    template <int D>
    PartitionInfo<D> LBTree<D>::reset_partition(typename DClass<D>::KeyDConst& key, Cost partition_size, Cost used_up, PartitionInfo<D> lbi) {
//	typename DClass<D>::KeyDConst root(0);
	// reduce the partition number
	lbi.partition_number--;
	madness::print("reset_partition: decremented partition_number to", lbi.partition_number);
	if (lbi.partition_number == 0) {
	    madness::print("reset_partition: partition_number == 0, launch_make_partition");
//	    this->world.gop.fence();
//	    launch_make_partition(lbi);
	    return lbi;
	} else {
            // reset the partition size
	    lbi.cost_left -= used_up;
	    used_up = 0;
	    Cost tpart = compute_partition_size(lbi.cost_left, lbi.partition_number+1);
	    // Reconsider partition size at every step.  If too small, leaves too much work for P0.
	    // If too large, leaves NO work for final processors.  So, strike a balance.
	    if ((tpart > partition_size) || (tpart*lbi.facter < partition_size)) {
		partition_size = tpart;
	    }
	    madness::print("reset_partition: continue to make_partition with", key);        
//	    lbi = send(impl.owner(key), &LBTree<D>::make_partition, key, partition_size, used_up, lbi);
	    send(impl.owner(key), &LBTree<D>::make_partition, key, partition_size, used_up, lbi);
	} 
	madness::print("reset_partition: about to return None");
	return lbi;
    }
*/
      

    /// Remove evidence of a claimed subtree from its foreparents' subtree cost

    /// Arguments: const Key<D> key -- node at which we begin
    ///            Cost c -- cost of what is to be removed
    /// Side effect: Node's subcost is updated
    template <int D>
    Void LBTree<D>::remove_cost(typename DClass<D>::KeyDConst& key, Cost c) {
        if (((int) key.level()) < 0) return None;
        typename DClass<D>::treeT::iterator it = impl.find(key);
        if (it == impl.end()) return None;
        typename DClass<D>::NodeD node = it->second;
        NodeData d = node.get_data();
        d.subcost -= c;
        node.set_data(d);
        impl.insert(key,node);
        if (key.level() > 0) {
            send(impl.owner(key.parent()), &LBTree<D>::remove_cost, key.parent(), c);
        }

	return None;
    }


    /// Compute the partition size: a straight quotient of the cost by the number of
    /// remaining partitions
    Cost compute_partition_size(Cost cost, unsigned int parts) {
        return (Cost) ceil(((double) cost)/((double) parts));
    }



     // Explicit instantiations for D=1:6

    template class LoadBalImpl<double,1>;
    template class LoadBalImpl<double,2>;
    template class LoadBalImpl<double,3>;
    template class LoadBalImpl<double,4>;
    template class LoadBalImpl<double,5>;
    template class LoadBalImpl<double,6>;

    template class LoadBalImpl<std::complex<double>,1>;
    template class LoadBalImpl<std::complex<double>,2>;
    template class LoadBalImpl<std::complex<double>,3>;
    template class LoadBalImpl<std::complex<double>,4>;
    template class LoadBalImpl<std::complex<double>,5>;
    template class LoadBalImpl<std::complex<double>,6>;

    template class LBTree<1>;
    template class LBTree<2>;
    template class LBTree<3>;
    template class LBTree<4>;
    template class LBTree<5>;
    template class LBTree<6>;
}
