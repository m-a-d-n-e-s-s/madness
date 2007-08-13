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
	std::vector<std::vector<typename DClass<D>::TreeCoords> > list_of_list;
	std::vector<typename DClass<D>::TreeCoords> klist;
	SharedPtr<std::vector<Cost> > costlist = SharedPtr<std::vector<Cost> >(new std::vector<Cost>());

	list_of_list = skeltree->find_partitions(costlist);

	if (this->f.get_impl()->world.mpi.rank() == 0) {
	    unsigned int shortest_list = 0, sl_index, lb_index;
	    Cost load_bal_cost = 0;
	    int count = list_of_list.size();
	    std::vector<unsigned int> len;
	    for (int i = 0; i < count; i++) {
		len.push_back(list_of_list[i].size());
		if ((len[i] < shortest_list) || (shortest_list == 0)) {
		    shortest_list = len[i];
		    sl_index = i;
		} else if ((len[i] == shortest_list) && ((*costlist)[i] < (*costlist)[sl_index])) {
		    // all things being equal, prefer better balance
		    shortest_list = len[i];
		    sl_index = i;
		}
		if (((*costlist)[i] < load_bal_cost) || (load_bal_cost == 0)) {
		    load_bal_cost = (*costlist)[i];
		    lb_index = i;
		} else if (((*costlist)[i] == load_bal_cost) && (len[i] < list_of_list[lb_index].size())) {
		    // all things being equal, prefer fewer cuts
		    load_bal_cost = (*costlist)[i];
		    lb_index = i;
		}
	    }
/*
	    madness::print("The load balance with the fewest broken links has cost", costlist[sl_index], "and",
		shortest_list-1, "broken links");
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
*/
	    CompCost ccleast = 0;
	    int cc_index = 0;
	    for (int i = 0; i < count; i++) {
		CompCost cctmp = compute_comp_cost((*costlist)[i], len[i]-1);
		if ((i==0) || (cctmp < ccleast)) {
		    ccleast = cctmp;
		    cc_index = i;
		}
	    }
/*
	    madness::print("The load balance with the best overall computational cost has cost",
			   costlist[cc_index], "and", len[cc_index]-1, "broken links");
	    for (unsigned int i = 0; i < len[cc_index]; i++) {
		list_of_list[cc_index][i].print();
	    }
*/
	    for (unsigned int i = 0; i < len[cc_index]; i++) {
		klist.push_back(list_of_list[cc_index][i]);
	    }
	    unsigned int ksize = klist.size();
	    this->f.get_impl()->world.gop.template broadcast<unsigned int>(ksize);
	    for (unsigned int i=0; i < ksize; i++) {
		this->f.get_impl()->world.gop.template broadcast<typename DClass<D>::TreeCoords>(klist[i]);
	    }
	} else {
            unsigned int ksize;
	    // Receive the broadcast of the final partitioning
            this->f.get_impl()->world.gop.template broadcast<unsigned int>(ksize);
            for (unsigned int i = 0; i < ksize; i++) {
                typename DClass<D>::TreeCoords t;
                this->f.get_impl()->world.gop.template broadcast<typename DClass<D>::TreeCoords>(t);
                klist.push_back(t);
            }
	}
//            madness::print("done with broadcast");

//	double t7 = MPI::Wtime();

//	madness::print("find_best_partition: total time =", t7-t0);
//	madness::print("find_best_partition: time in loop =", t6-t1);

	if (this->f.get_impl()->world.mpi.rank() == 0) {
	    madness::print("find_best_partition: number of broken links =",
		klist.size()-1);
	}
        return klist;	
    };

    /// find_partitions performs the "melding" algorithm for load balancing: it recursively melds 
    /// and partitions the tree until it has found all possible configurations.
    template <int D>
    std::vector< std::vector<typename DClass<D>::TreeCoords> > LBTree<D>::find_partitions(SharedPtr<std::vector<Cost> > costlist) {
//	double t0 = MPI::Wtime();
	bool keep_going = true;
	std::vector< std::vector<typename DClass<D>::TreeCoords> > klists;
	if (this->world().nproc() == 1) {
	    klists.push_back(std::vector<TreeCoords<D> >(1,TreeCoords<D>(Key<D>(0), 0)));
	    return klists;
	}
        if (this->world().mpi.rank() != 0) {
	    while (keep_going) {
// Worker processes fence while manager queries them for fix_cost
	        this->world().gop.fence();
// Set all elements to true before rollup
		this->reset(true);
// Fence so that everybody starts rollup at the same time
	        this->world().gop.fence();
        	this->rollup();
// Fence so nobody resets until everyone's done with rollup
	        this->world().gop.fence();
		this->reset(false);
        	this->world().gop.fence();
// Fence while manager queries for depth_first_partition
        	this->world().gop.fence();
// Find out whether we're in for another round of load balancing
                this->world().gop.template broadcast<bool>(keep_going);
	    }
	    klists.push_back(std::vector<TreeCoords<D> >(1,TreeCoords<D>(Key<D>(0), 0)));
            return klists;
	}

	// The manager process coordinates the melding algorithm for load balancing, keeping a list of
	// lists of the configurations suggested by the melding algorithm, and selecting the best
	// configuration at the end.
        int npieces = this->world().nproc();
        int count = 0;
        std::vector<std::vector<typename DClass<D>::TreeCoords> > list_of_list;
        std::vector<typename DClass<D>::TreeCoords> emptylist;

        Cost totalCost = 0;
        typename DClass<D>::KeyD root(0);

//	double t1 = MPI::Wtime();
	while (keep_going) {
//	    double t2 = MPI::Wtime();
            this->fix_cost(root);
	    this->world().gop.fence();
//	    double t3 = MPI::Wtime();
	    this->reset(true);
	    this->world().gop.fence();
            this->rollup();
	    this->world().gop.fence();
	    this->reset(false);
            this->world().gop.fence();
//	    double t4 = MPI::Wtime();
            list_of_list.push_back(emptylist);
            costlist->push_back(0);
            totalCost = this->depth_first_partition(root, &list_of_list[count], npieces,
                    totalCost, &(*costlist)[count]);
//	    double t5 = MPI::Wtime();

//	    madness::print("find_best_partition: time for fix_cost number", count, "=", t3-t2);
//	    madness::print("find_best_partition: time for rollup number", count, "=", t4-t3);
//	    madness::print("find_best_partition: time for depth_first_partition number", count, "=", t5-t4);
	    int lolcsize = list_of_list[count].size();
/*
	    madness::print("find_best_partition: partition number", count, ":");
	    for (int k = 0; k < lolcsize; k++) {
		list_of_list[count][k].print();
	    }
	    madness::print("");
*/
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
                keep_going = false;
                list_of_list.erase(list_of_list.begin()+count);
                keep_going = false;
		count-=1;
            }
            count++;
            this->world().gop.fence();
            this->world().gop.template broadcast<bool>(keep_going);
	}
//	double t6 = MPI::Wtime();

	return list_of_list;
    }


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


    /// fix_cost resets the tree after the load balancing and melding have been performed, before the next
    /// round of load balancing.
    /// Argument: const Key<D> key -- the node to be reset
    /// Return: Cost of subtree rooted at key
    /// Side effect: subcost (the cost of the subtree rooted at key) is reset
    /// Communication: Just what's required for find and insert

    template <int D>
    Cost LBTree<D>::fix_cost(typename DClass<D>::KeyDConst& key) {
        typename DClass<D>::treeT::iterator it = this->find(key);
        if (it == this->end()) return 0;

        typename DClass<D>::NodeD node = it->second;
        NodeData d = node.get_data();
        d.subcost = d.cost;
// Recursively set the subcost of this node's children, and add it to
// this node's subcost.
        if (node.has_children()) {
            for (KeyChildIterator<D> kit(key); kit; ++kit) {
                d.subcost += this->fix_cost(kit.key());
            }
        }
        node.set_data(d);
        this->insert(key,node);
        return d.subcost;
    }


    /// depth_first_partition finds partitions of trees by calling make_partition.  It figures
    /// out the size of the partitions to be made by make_partition at each iteration.
    /// Arguments: const Key<D> key -- node at which we begin
    ///            vector<TreeCoords<D> >* klist -- list of subtree root nodes obtained from partitioning
    ///            unsigned int npieces -- number of pieces to partition the tree into
    ///            Cost totalcost -- cost of tree rooted at key
    ///            Cost *maxcost -- the maximum cost allowed in a partition
    /// Return: Cost -- the Cost of what has been used up in the partition so far
    /// Side effect: klist and maxcost are updated
    template <int D>
    Cost LBTree<D>::depth_first_partition(typename DClass<D>::KeyDConst& key,
            std::vector<typename DClass<D>::TreeCoords>* klist, unsigned int npieces,
            Cost totalcost, Cost *maxcost) {
        if (totalcost == 0) {
            totalcost = this->compute_cost(key);
        }

        Cost cost_left = totalcost;
        int parts_left = npieces;
        *maxcost = 0;
        Cost partition_size = 0;
	double facter = 1.1;

        for (int i = npieces-1; i >= 0; i--) {
//	    madness::print("");
//	    madness::print("Beginning partition number", i);
            std::vector<typename DClass<D>::KeyD> tmplist;
            Cost tpart = compute_partition_size(cost_left, parts_left);
	    // Reconsider partition size at every step.  If too small, leaves too much work for P0.
	    // If too large, leaves NO work for final processors.  So, strike a balance.
            if ((tpart > partition_size) || (tpart*facter < partition_size)) {
                partition_size = tpart;
            }
//            madness::print("depth_first_partition: partition_size =", partition_size);
            Cost used_up = 0;
            bool at_leaf = false;
            used_up = this->make_partition(key, &tmplist, partition_size, (i==0), used_up, &at_leaf);
            if (*maxcost < used_up) *maxcost = used_up;
            cost_left -= used_up;
            parts_left--;
            for (unsigned int j = 0; j < tmplist.size(); j++) {
                klist->push_back(typename DClass<D>::TreeCoords(typename DClass<D>::KeyD(tmplist[j]), i));
            }
        }
        return totalcost;
    }


    /// rollup traverses the tree, calling meld upon Nodes that have leaf children
    /// Arguments: const Key<D> key -- node at which we begin
    /// Side effect: Nodes are changed by meld
    /// Communication: just finding the nodes that match a given key
    template <int D>
    void LBTree<D>::rollup() {
	for (typename DClass<D>::treeT::iterator it = this->begin(); it != this->end(); ++it) {
	    typename DClass<D>::KeyD key = it->first;
	    typename DClass<D>::NodeD node = it->second;
	    if (node.has_children()) {
// First, check to see if it has any leaf children
	        bool has_leaf_child = false;
	        for (KeyChildIterator<D> kit(key); kit; ++kit) {
            	    typename DClass<D>::treeT::iterator itc = this->find(kit.key());
		    if (itc != this->end()) {
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
		    this->insert(key,node);
	        }
	    }
	}
    }

    /// reset sets the is_taken variable within all local nodes to the value specified
    /// Arguments: bool taken -- value to set is_taken to
    /// Communication: none (local iterator)
    template <int D>
    void LBTree<D>::reset(bool taken) {
	for (typename DClass<D>::treeT::iterator it = this->begin(); it != this->end(); ++it) {
	    typename DClass<D>::KeyD key = it->first;
	    typename DClass<D>::NodeD node = it->second;
	    NodeData d = node.get_data();

	    d.is_taken = taken;
	    node.set_data(d);
	    this->insert(key,node);
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
		typename DClass<D>::treeT::iterator itc = this->find(kit.key());
		typename DClass<D>::NodeD c = itc->second;
		NodeData d = c.get_data();
// if the child has no children and the is_taken flag is set to true, then
// this child is eligible to be melded into parent
		if ((!c.has_children()) && (d.is_taken)) {
		    Cost cost = d.cost;
		    if ((cost < cheapest) || (not_yet_found)) {
// if this child has the cheapest cost yet, then clear the list and add
// this child to the list of children to be melded to the parent
			not_yet_found = false;
			cheapest = cost;
			mylist.clear();
			mylist.push_back(i);
		    } else if (cost == cheapest) {
// if this child's cost is equal to the cheapest cost found so far, then
// add this child to the list of children to be melded into parent
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
		this->erase(kit.key());
		node.set_child(mylist[j], false);
		d.cost += cheapest;
		j++;
		if (j == mlsize) break;
	    }
	    i++;
	}
	node.set_data(d);
	this->insert(key, node);
    }

    /// compute_cost resets the subtree cost value for a tree and its descendants.
    /// Arguments: const Key<D> key -- node at which we begin
    /// Return: Cost -- the Cost of the subtree headed by key
    /// Side effect: subcost is updated
    /// Communication: find and insert
    template <int D>
    Cost LBTree<D>::compute_cost(typename DClass<D>::KeyDConst& key) {
        Cost cost = 0;
        typename DClass<D>::treeT::iterator it = this->find(key);
        if (it == this->end()) return cost;

        typename DClass<D>::NodeD node = it->second;
        for (KeyChildIterator<D> kit(key); kit; ++kit) {
            cost += this->compute_cost(kit.key());
        }
        NodeData d = node.get_data();
        cost += d.cost;

        d.subcost = cost;
        node.set_data(d);
        this->insert(key,node);
        return cost;
    }


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
    Cost LBTree<D>::make_partition(typename DClass<D>::KeyDConst& key,
	 			       std::vector<typename DClass<D>::KeyD>* klist, Cost partition_size, 
				       bool last_partition, Cost used_up, bool *at_leaf) {
// The fudge factor is the fraction by which you are willing to let the
// partitions exceed the ideal partition size
        double fudge_factor = 0.1;
        Cost maxAddl = (Cost) (fudge_factor*partition_size);
	std::vector<bool> my_children_status; 
	Cost my_current_cost = used_up;

        typename DClass<D>::treeT::iterator it = this->find(key);
        if (it == this->end()) {
            return used_up;
        }

        typename DClass<D>::NodeD node = it->second;
        NodeData d = node.get_data();

        if (d.is_taken) {
// this key is taken and cannot be used to fill any additional partitions
            return used_up;
        }

        // if either we're at the last partition, the partition is currently empty
        // and this is a single item, or there is still room in the partition and
        // adding this to it won't go above the fudge factor,
        // then add this piece to the partition.
        if ((last_partition) || ((used_up == 0) && (!node.has_children())) ||
                ((used_up < partition_size) && (d.subcost+used_up <= partition_size+maxAddl))) {
            // add to partition
            klist->push_back(typename DClass<D>::KeyD(key));
            d.is_taken = true;
            used_up += d.subcost;
            // REMOVE COST FROM FOREPARENTS 
            this->remove_cost(key.parent(), d.subcost);
            node.set_data(d);
            this->insert(key,node);
        } else if (used_up < partition_size) {
            // try this node's children (if any)
            if (node.has_children()) {
                int i = 0;
                for (KeyChildIterator<D> kit(key); kit; ++kit) {
                    if (node.has_child(i)) {
                        used_up = this->make_partition(kit.key(), klist, partition_size, last_partition,
                                                              used_up, at_leaf);
			if (used_up == 0) {
// if you've made no progress by going to the children (used_up == 0 means
// you've not put anything in the partition)
			    my_children_status.push_back(false);
			} else {
// You added a node's child to the partition when you called make_partition
// on that child
			    my_children_status.push_back(true);
			}
                        if ((*at_leaf) || (used_up >= partition_size)) {
// if we're at the leaf or we have used up the entire partition, we're outta here
                            break;
                        }
                    } else if (used_up == 0) {
// child kit.key() does not exist but the partition is currently empty
			my_children_status.push_back(false);
		    }
                    i++;
                }
		if (used_up == 0) {
		    int mcs_size = my_children_status.size();
		    bool no_avail_children = true;
		    for (int i = 0; i < mcs_size; i++) {
			if (my_children_status[i]) {
// if this key has an existing, not-previously-taken, child, then there is an
// available child, so don't add this node to partition
			    no_avail_children = false;
			    break;
			}
		    }
		    if (no_avail_children) {
// there are no available children, meaning that the children could exist but
// are all taken by other partitions already
			klist->push_back(typename DClass<D>::KeyD(key));
			d.is_taken = true;
			used_up += d.subcost;
			this->remove_cost(key.parent(), d.subcost);
			node.set_data(d);
			this->insert(key,node);
		    }
		} else if (used_up != my_current_cost) {
// Effectively, we're at a leaf, because we've used up all the children of
// this node but not the node itself
		    *at_leaf = true;
		}
            } else {
                *at_leaf = true;
            }
        }
        return used_up;
    }

    /// Remove evidence of a claimed subtree from its foreparents' subtree cost

    /// Arguments: const Key<D> key -- node at which we begin
    ///            Cost c -- cost of what is to be removed
    /// Side effect: Node's subcost is updated
    template <int D>
    void LBTree<D>::remove_cost(typename DClass<D>::KeyDConst& key, Cost c) {
        if (((int) key.level()) < 0) return;
        typename DClass<D>::treeT::iterator it = this->find(key);
        if (it == this->end()) return;
        typename DClass<D>::NodeD node = it->second;
        NodeData d = node.get_data();
        d.subcost -= c;
        if (key.level() > 0) {
            this->remove_cost(key.parent(), c);
        }
        node.set_data(d);
        this->insert(key,node);
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
