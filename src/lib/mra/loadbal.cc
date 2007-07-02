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

    typedef int Cost;
    typedef double CompCost;

    /// find_best_partition performs the "melding" algorithm for load balancing: it recursively melds 
    /// and partitions the tree until it has found all possible configurations.
    template <typename T, int D>
    std::vector<typename DClass<D>::TreeCoords> LoadBalImpl<T,D>::find_best_partition() {
//	madness::print("find_best_partition: at beginning");
	double t0 = MPI::Wtime();
	bool keep_going = true;
	std::vector<typename DClass<D>::TreeCoords> klist;
        if (this->f.get_impl()->world.mpi.rank() != 0) {
//	    madness::print("find_best_partition: just before while loop");
	    while (keep_going) {
//		madness::print("find_best_partition: just before first fence");
	        this->f.get_impl()->world.gop.fence();
//		madness::print("find_best_partition: about to rollup()");
		this->skeltree->reset(true);
	        this->f.get_impl()->world.gop.fence();
        	this->skeltree->rollup();
	        this->f.get_impl()->world.gop.fence();
		this->skeltree->reset(false);
//		madness::print("find_best_partition: finished with rollup()");
        	this->f.get_impl()->world.gop.fence();
//		madness::print("find_best_partition: finished with first fence in a row");
        	this->f.get_impl()->world.gop.fence();
//		madness::print("find_best_partition: finished with second fence in a row");
                this->f.get_impl()->world.gop.template broadcast<bool>(keep_going);
//		madness::print("find_best_partition: received keep_going =", keep_going);
	    }
            unsigned int ksize;
	    // Then, they receive the broadcast of the final partitioning
            this->f.get_impl()->world.gop.template broadcast<unsigned int>(ksize);
            for (unsigned int i = 0; i < ksize; i++) {
                typename DClass<D>::TreeCoords t;
                this->f.get_impl()->world.gop.template broadcast<typename DClass<D>::TreeCoords>(t);
                klist.push_back(t);
            }
//            madness::print("done with broadcast");
            return klist;
	}

	// The manager process coordinates the melding algorithm for load balancing, keeping a list of
	// lists of the configurations suggested by the melding algorithm, and selecting the best
	// configuration at the end.
        unsigned int npieces = this->f.get_impl()->world.nproc();
        int count = 0;
        std::vector<std::vector<typename DClass<D>::TreeCoords> > list_of_list;
        std::vector<typename DClass<D>::TreeCoords> emptylist;
        std::vector<Cost> costlist;
        Cost totalCost = 0;
        typename DClass<D>::KeyD root(0);

//	madness::print("find_best_partition: just before while loop");
	double t1 = MPI::Wtime();
	while (keep_going) {
	    double t2 = MPI::Wtime();
//	    madness::print("find_best_partition: just before fix_cost");
            this->skeltree->fix_cost(root);
//            madness::print("find_best_partition: just finished fix_cost");
	    this->f.get_impl()->world.gop.fence();
//            madness::print("find_best_partition: about to rollup");
	    double t3 = MPI::Wtime();
	    this->skeltree->reset(true);
	    this->f.get_impl()->world.gop.fence();
            this->skeltree->rollup();
	    this->f.get_impl()->world.gop.fence();
	    this->skeltree->reset(false);
//            madness::print("find_best_partition: just finished rollup");
            this->f.get_impl()->world.gop.fence();
	    double t4 = MPI::Wtime();
            list_of_list.push_back(emptylist);
//	    madness::print("find_best_partition: about to depth_first_partition this tree:");
//	    this->skeltree->print(root);
            costlist.push_back(0);
//            madness::print("find_best_partition: about to depth_first_partition");
            totalCost = this->skeltree->depth_first_partition(root, &list_of_list[count], npieces,
                    totalCost, &costlist[count]);
//            madness::print("find_best_partition: finished with depth_first_partition");
	    double t5 = MPI::Wtime();

	    madness::print("find_best_partition: time for fix_cost number", count, "=", t3-t2);
	    madness::print("find_best_partition: time for rollup number", count, "=", t4-t3);
	    madness::print("find_best_partition: time for depth_first_partition number", count, "=", t5-t4);
//	    int lolcsize = list_of_list[count].size();
//	    madness::print("find_best_partition: partition number", count, ":");
//	    for (int k = 0; k < lolcsize; k++) {
//		list_of_list[count][k].print();
//	    }
//	    madness::print("");
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
//	    madness::print("find_best_partition: before fence");
            this->f.get_impl()->world.gop.fence();
//	    madness::print("find_best_partition: after fence, about to broadcast keep_going =", keep_going);
            this->f.get_impl()->world.gop.template broadcast<bool>(keep_going);
	}
	double t6 = MPI::Wtime();
        unsigned int shortest_list = 0, sl_index, lb_index;
        Cost load_bal_cost = 0;
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

	madness::print("The load balance with the fewest broken links has cost", costlist[sl_index], "and",
		shortest_list-1, "broken links");
        for (unsigned int i = 0; i < shortest_list; i++) {
            list_of_list[sl_index][i].print();
        }
	madness::print("");
	madness::print("The load balance with the best balance has cost", load_bal_cost, "and", 
		list_of_list[lb_index].size()-1, "broken links");
        for (unsigned int i = 0; i < list_of_list[lb_index].size(); i++) {
            list_of_list[lb_index][i].print();
        }
	madness::print("");

        CompCost ccleast = 0;
        int cc_index;
        for (int i = 0; i < count; i++) {
            CompCost cctmp = compute_comp_cost(costlist[i], len[i]-1);
            if ((i==0) || (cctmp < ccleast)) {
                ccleast = cctmp;
                cc_index = i;
            }
        }
	madness::print("The load balance with the best overall computational cost has cost",
		costlist[cc_index], "and", len[cc_index]-1, "broken links");
        for (unsigned int i = 0; i < len[cc_index]; i++) {
            list_of_list[cc_index][i].print();
        }
        for (unsigned int i = 0; i < len[cc_index]; i++) {
            klist.push_back(list_of_list[cc_index][i]);
        }

//        madness::print("about to do broadcast");
        unsigned int ksize = klist.size();
        this->f.get_impl()->world.gop.template broadcast<unsigned int>(ksize);
        for (unsigned int i=0; i < ksize; i++) {
            this->f.get_impl()->world.gop.template broadcast<typename DClass<D>::TreeCoords>(klist[i]);
        }
//        madness::print("done with broadcast");

	double t7 = MPI::Wtime();

	madness::print("find_best_partition: total time =", t7-t0);
	madness::print("find_best_partition: time in loop =", t6-t1);

        return klist;
    }


    /// fix_cost resets the tree after the load balancing and melding have been performed, before the next
    /// round of load balancing.
    /// Argument: const Key<D> key -- the node to be reset
    /// Return: Cost of subtree rooted at key
    /// Side effect: subcost (the cost of the subtree rooted at key) is reset
    /// Communication: Just what's required for find and insert

    template <int D>
    Cost LBTree<D>::fix_cost(typename DClass<D>::KeyDConst& key) {
//         madness::print("fix_cost: key =", key, " is about to be looked for");
        typename DClass<D>::treeT::iterator it = this->find(key);
//         madness::print("fix_cost: key =", key, "exists:", !(it == this->end()));
        if (it == this->end()) return 0;

        typename DClass<D>::NodeD node = it->second;
//         madness::print("fix_cost: for key", key, "got node", node);
        NodeData d = node.get_data();
//         madness::print("fix_cost: got data from node");
        d.subcost = d.cost;
//         madness::print("fix_cost: assigned node cost to subcost");
        if (node.has_children()) {
//             madness::print("fix_cost: node has children");
            for (KeyChildIterator<D> kit(key); kit; ++kit) {
                d.subcost += this->fix_cost(kit.key());
            }
        }
        node.set_data(d);
//         madness::print("fix_cost: about to insert key =", key, ",", node.get_data());
        this->insert(key,node);
//         madness::print("fix_cost: inserted node");
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
//madness::print("depth_first_partition: at very beginning");
        if (totalcost == 0) {
            totalcost = this->compute_cost(key);
        }
//        madness::print("depth_first_partition: totalcost =", totalcost);

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
//	    madness::print("rollup: key", key, "is being rolled upon");
	    typename DClass<D>::NodeD node = it->second;
	    if (node.has_children()) {
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
//		    madness::print("rollup: key", key, "has leaf child");
//		    madness::print("rollup: key", key, "has cost", node.get_data().cost);
		    this->meld(it);
		    node = it->second;
//		    madness::print("rollup: after melding, key", key, "has cost", node.get_data().cost);
	        }
		else {
//		    madness::print("rollup: key", key, "doesn't have leaf child");
		}
	        NodeData d = node.get_data();
	        if (d.is_taken) {
		    d.is_taken = false;
		    node.set_data(d);
		    this->insert(key,node);
	        }
	    }
	}
    }

    template <int D>
    void LBTree<D>::reset(bool taken) {
	for (typename DClass<D>::treeT::iterator it = this->begin(); it != this->end(); ++it) {
	    typename DClass<D>::KeyD key = it->first;
	    typename DClass<D>::NodeD node = it->second;
	    NodeData d = node.get_data();
//	    madness::print("reset: key", key, "is being set to is_taken = false");
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

//	madness::print("meld: melding on key", key);

	for (KeyChildIterator<D> kit(key); kit; ++kit) {
	    if (node.has_child(i)) {
		typename DClass<D>::treeT::iterator itc = this->find(kit.key());
		typename DClass<D>::NodeD c = itc->second;
		NodeData d = c.get_data();
//		madness::print("meld: key", key, "has child", kit.key(), "with data", c);
		if ((!c.has_children()) && (d.is_taken)) {
		    Cost cost = d.cost;
//		    madness::print("meld: key", key, "has eligible leaf child", kit.key(), "with cost", cost);
		    if ((cost < cheapest) || (not_yet_found)) {
			not_yet_found = false;
			cheapest = cost;
			mylist.clear();
//			madness::print("meld: cleared list and then pushed child", i, "onto list");
			mylist.push_back(i);
		    } else if (cost == cheapest) {
			mylist.push_back(i);
//			madness::print("meld: pushed child", i, "onto list");
		    }
		}
	    }
	    i++;
	}
	if (not_yet_found) {
	    // this node has no leaf children
//	    madness::print("meld: key", key, "actually doesn't have any eligible leaf children!");
	    return;
	}
//	madness::print("meld: contents of mylist = ", mylist);
	NodeData d = node.get_data();
	i = 0;
	int j = 0, mlsize = mylist.size();
	for (KeyChildIterator<D> kit(key); kit; ++kit) {
	    if (mylist[j] == i) {
//		madness::print("meld: erasing", key, "'s child number", mylist[j], ":", kit.key());
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
//	madness::print("meld: key", key, "has cost", node.get_data().cost);
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
//    madness::print("at beginning of make_partition: at_leaf =", *at_leaf);
        double fudgeFactor = 0.1;
        Cost maxAddl = (Cost) (fudgeFactor*partition_size);
	std::vector<bool> my_children_status; 
	Cost my_current_cost = used_up;

        typename DClass<D>::treeT::iterator it = this->find(key);
        if (it == this->end()) {
            return used_up;
        }

        typename DClass<D>::NodeD node = it->second;
        NodeData d = node.get_data();

//        madness::print("make_partition: data for key", key, ":", d);
//        madness::print("make_partition: partition_size =", partition_size, ", last_partition =", last_partition, ", used_up =", used_up);

        if (d.is_taken) {
//            madness::print("make_partition: this key is taken");
            return used_up;
        }

//        madness::print("make_partition: back to key", key);

        // if either we're at the last partition, the partition is currently empty
        // and this is a single item, or there is still room in the partition and
        // adding this to it won't go above the fudge factor,
        // then add this piece to the partition.
        if ((last_partition) || ((used_up == 0) && (!node.has_children())) ||
                ((used_up < partition_size) && (d.subcost+used_up <= partition_size+maxAddl))) {
            // add to partition
//            madness::print("make_partition: adding to partition", key);
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
//                        madness::print("make_partition:", key, "recursively calling", kit.key());
                        used_up = this->make_partition(kit.key(), klist, partition_size, last_partition,
                                                              used_up, at_leaf);
			if (used_up == 0) {
//			    madness::print("make_partition: my child", kit.key(), "is taken, so it is not available for me to use");
			    my_children_status.push_back(false);
			} else {
//			    madness::print("make_partition: my child", kit.key(), "is not taken");
			    my_children_status.push_back(true);
			}
                        if ((*at_leaf) || (used_up >= partition_size)) {
                            break;
                        }
                    } else if (used_up == 0) {
//			madness::print("make_partition: my child", kit.key(), "does not exist but used_up == 0");
			my_children_status.push_back(false);
		    }
                    i++;
                }
		if (used_up == 0) {
		    int mcs_size = my_children_status.size();
		    bool no_avail_children = true;
		    for (int i = 0; i < mcs_size; i++) {
			if (my_children_status[i]) {
//			    madness::print("make_partition: child number", i, "is available, so don't do anything");
			    no_avail_children = false;
			    break;
			}
		    }
		    if (no_avail_children) {
//			madness::print("make_partition: parent with only taken children, adding", key);
			klist->push_back(typename DClass<D>::KeyD(key));
			d.is_taken = true;
			used_up += d.subcost;
			this->remove_cost(key.parent(), d.subcost);
			node.set_data(d);
			this->insert(key,node);
		    }
		} else if (used_up != my_current_cost) {
//		    madness::print("make_partition: setting at_leaf = true because we've used up all my children but not myself");
		    *at_leaf = true;
		}
            } else {
//                madness::print("make_partition: about to set at_leaf = true");
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
//        madness::print("remove_cost: key", key, "owner =", owner(key));
//        this->get_mypmap().print();
        if (((int) key.level()) < 0) return;
        typename DClass<D>::treeT::iterator it = this->find(key);
//        madness::print("remove_cost: found key");
        if (it == this->end()) return;
        typename DClass<D>::NodeD node = it->second;
        NodeData d = node.get_data();
//        madness::print("remove_cost: got data");
        d.subcost -= c;
        if (key.level() > 0) {
            this->remove_cost(key.parent(), c);
        }
//        madness::print("remove_cost: before setting, data =", d);
        node.set_data(d);
//        madness::print("remove_cost: after setting, data =", node.get_data());
        this->insert(key,node);
//        madness::print("remove_cost: after inserting, data = ", node.get_data());
    }


    /// Compute the partition size: a straight quotient of the cost by the number of
    /// remaining partitions
    Cost compute_partition_size(Cost cost, unsigned int parts) {
        return (Cost) ceil(((double) cost)/((double) parts));
    }


    /// Compute the cost of a given configuration: a weighted sum of the cost of the
    /// maximally-loaded process and the total number of broken links.
    /// In the future, the factors should be calibrated for a given machine, either 
    /// during configuration and setup or in real time at the beginning of the program
    /// or upon the creation of a LoadBalImpl.
    /// Arguments: Cost c -- maximum amount of cost assigned to a node
    ///            int n -- number of broken links
    /// Return: Cost -- the Cost of what has been used up in the partition so far
    CompCost compute_comp_cost(Cost c, int n) {
        CompCost compcost;
        CompCost cfactor = 0.01, nfactor = 1.0;
        compcost = cfactor*c + nfactor*n;
        return compcost;
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
