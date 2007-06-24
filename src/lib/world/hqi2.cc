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

  $LastChangedDate$
  $Rev$
*/

  
#define WORLD_INSTANTIATE_STATIC_TEMPLATES
//#ifndef LOADBAL_H
#include <world/loadbal.h>
//#endif

using namespace madness;
using namespace std;


int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);
    redirectio(world);
    //xterm_debug("hqi2", 0);
    ProcessID me = world.rank();


    try {
	Vector<unsigned long,2> k0, k1, k2;
	k0[0] = 0; k0[1] = 0;
	k1[0] = 0; k1[1] = 1;
	k2[0] = 1; k2[1] = 1;
	vector<DClass<2>::TreeCoords> v;
	v.push_back(DClass<2>::TreeCoords(DClass<2>::KeyD(0,k0),0));
	v.push_back(DClass<2>::TreeCoords(DClass<2>::KeyD(1,k1),1));
	v.push_back(DClass<2>::TreeCoords(DClass<2>::KeyD(1,k2),1));

	DClass<2>::KeyD root(0,k0);
	DClass<2>::treeT tree(world,DClass<2>::MyPmap(world,v));
//	DClass<2>::treeT tree(world,MyPmap<2>(v));
//	treeT tree(world,MyPmapD(1));
	print("Made tree");
	print("About to make tree1");
//	DClass<2>::treeT tree1(world,MyPmap<2>(v));
	DClass<2>::treeT tree1(world,DClass<2>::MyPmap(world,v));
	if (me == 0) { 
	    print("About to build tree");
	    build_tree<2>(tree,root);
	    print("About to convert to tree1");
	    migrate<2>(tree,tree1);
	    print("built tree1");
//	    print_tree<2>(tree1,root);
	    print("printed tree1");
	    print("");
	}
	print("Built tree");
	world.gop.fence();
        print("Done fencing");
	print("Now we're going to print the tree");
	print("");
	if (me == 1) {
	    print("About to print tree");
	    print_tree<2>(tree,root);
	    print("Printed tree");
	}
	print("Done printing tree");
	print("");
	vector<DClass<2>::TreeCoords> klist;
	if (me == 0) {
	    unsigned int npieces = world.nproc();
	    print("about to find best partition");
	    findBestPartition<2>(tree, root, &klist, npieces);
	    print("");
	}
	world.gop.fence();
	// Now, broadcast klist
	if (me == 0) {
	    unsigned int ksize = klist.size();
	    world.gop.broadcast<unsigned int>(ksize);
	    print("broadcasted ksize");
	    for (unsigned int i = 0; i < ksize; i++) {
		world.gop.broadcast<DClass<2>::TreeCoords>(klist[i]);
		klist[i].print();
	    }
	    print("done broadcasting klist");
	}
	else {
	    unsigned int ksize;
	    world.gop.broadcast<unsigned int>(ksize);
	    print("broadcasted ksize");
	    for (unsigned int i = 0; i < ksize; i++) {
		DClass<2>::TreeCoords t;
		world.gop.broadcast<DClass<2>::TreeCoords>(t);
		klist.push_back(t);
		t.print();
	    }
	    print("done broadcasting klist");
	}
	DClass<2>::treeT tree2(world,DClass<2>::MyPmap(world,klist));
	if (me == 0) {
	    migrate<2>(tree1, tree2);
	    print("copied tree1 to tree2");
	    print_tree<2>(tree2, root);
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
