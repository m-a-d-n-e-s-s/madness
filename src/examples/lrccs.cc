/*
 * lrccs.cc
 *
 *  Created on: 4 Jan 2017
 *      Author: kottmanj
 */

/*
/*
 * lrccs.cc
 *
 *  Created on: Aug 10, 2015
 *      Author: kottmanj
 */
/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

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

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES


/*!
  \file examples/tdhf.cc
  \brief compute the time-dependent HF equations (currently CIS approximation)

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
  /trunk/src/apps/examples/tdhf.cc>here</a>.

 */
#include <chem/TDHF.h>

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && defined(HAVE_UNISTD_H)
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
static inline int file_exists(const char * inpname)
{
    struct stat buffer;
    int rc = stat(inpname, &buffer);
    return (rc==0);
}
#endif


using namespace madness;

// test routine which tests h2o molecule
bool test_lrccs(World& world){

	std::vector<bool> results;
	std::vector<std::string> nemo_param={"none", "slater"};
	for(const auto nemo_type:nemo_param){
	bool passed=true;
	// print the input
	// primitive version since we hopefuly have a nice data format in the future
	if(world.rank()==0){
		const std::string filename ="lrccs_test_input_h2o_"+nemo_type;
		if(world.rank()==0) std::cout << "writing test input into file: " << filename << "\n";
		std::ofstream outfile(filename);
		  if(!outfile.is_open()) {
		    std::cerr << "Couldn't open " << filename  << std::endl;
		    return false;
		  }

		  outfile << "plot\nend"
				  << "\n\n"
				  << "\ndft\n canon \n xc hf\n econv 1.e-4\n dconv 1.e-3\n protocol 1.e-4\n L 30.0"
				  << "\n nuclear_corrfac " << nemo_type
				  <<"\nend"
		  	  	  << "\n\n"
		  	  	  << "\nresponse \n thresh 1.e-4\n excitations 4\n freeze 1 \nend"
				  << "\n\n"
				  << "\ngeometry"
				  << "\n o   0.00000000000000      0.00000000000000     -0.74803583254128"
				  << "\n h   1.43358660382183      0.00000000000000      0.37401791627063"
				  << "\n h  -1.43358660382183      0.00000000000000      0.37401791627063"
		  	  	  << "\nend\n\n";
		  outfile.close();



		  // ground state calculation
		  std::shared_ptr<SCF> calc(new SCF(world,filename.c_str()));
		  Nemo nemo(world,calc);
		  if (world.rank()==0) {
	    		calc->molecule.print();
	    		print("\n");
	    		calc->param.print(world);
		  }
		  double hf_energy = nemo.value();

		  // LR calculation
		  TDHF tdhf(world,nemo,filename);
		  std::vector<CC_vecfunction> roots=tdhf.solve_cis();

		  std::vector<double> expected_results = {0.3182483205,0.3794405844,0.3999497675,0.4088869291};
		  std::vector<double> results;
		  for(const auto& x:roots) results.push_back(x.omega);


		  if(world.rank()==0){
			  std::cout << "\n\n NEMO=" << nemo_type << "TEST FOR H2O ENDED:\n";

			  for(size_t i=0;i<results.size();++i){
				  const double err = expected_results[i]-results[i];
				  if(fabs(err)<1.e-3) std::cout << "found root " << i << "\n";
				  else{
					  std::cout << "cound not find root " << i << "\n";
					  passed=false;
				  }
			  }

			 if(passed) std::cout << "H2O TEST PASSED! :)\n";
			 else std::cout << "H2O TEST FAILED! :(\n";
			 results.push_back(passed);
		  }
	}
	}
bool result=true;
for(const auto x:results) if(x==false) result=false;
if(world.rank()==0){
	std::cout << "Test Results:\n";
	for(size_t i=0;i<results.size();++i) std::cout << "nemo=" << nemo_param[i] << " : " << results[i] << "\n";
}

return result;
}

int main(int argc, char** argv) {

	initialize(argc, argv);

	World world(SafeMPI::COMM_WORLD);

	// read out command line arguments
	bool analyze_only=false;		// default: compute the excitations

    for(int i = 1; i < argc; i++) {
        const std::string arg=argv[i];
        if (arg=="--analyze") analyze_only=true;
    }


	if (world.rank() == 0) {
		print("\n  MADNESS LINEAR-RESPONSE SUITE  \n");
		printf("starting at time %.1f\n", wall_time());
		print("\nmain() compiled at ",__TIME__," on ",__DATE__);

	}
	startup(world,argc,argv);
	std::cout.precision(6);
	FunctionDefaults<3>::set_truncate_mode(1);
	print("Truncate mode set to ",FunctionDefaults<3>::get_truncate_mode());

#ifdef GITREVISION
	const  char* gitrev =  GITREVISION;
	const std::string gitrevision(gitrev);
	if (world.rank()==0) {
		print("           git revision ...",gitrevision);
	}
#endif

    // Process 0 reads input information and broadcasts
    const char * inpname = "input";
    bool do_test=false;
    if(world.rank()==0) std::cout << "given arguments with function call:\n";
    for (int i=1; i<argc; i++) {
        if (argv[i][0] != '-') {
            inpname = argv[i];
            break;
        }
        if(world.rank()==0) std::cout << argv[i] << "\n";
        do_test=(std::string(argv[i])=="-test");
    }
    if (world.rank() == 0) print("input filename: ", inpname);
    if (!file_exists(inpname)) {
        throw "input file not found!";
    }

    if(world.rank()==0 and do_test) std::cout << "\n\n\n\n\n\n !!!!!!DEMADED TO PERFORM TESTS!!!!! \n\n\n\n\n\n\n";
    if(do_test){
    	test_lrccs(world);
    }else{
    	// Make reference
    	const std::string input = inpname;
    	//SCF calc(world,input.c_str());
    	std::shared_ptr<SCF> calc(new SCF(world,input.c_str()));
    	Nemo nemo(world,calc);
    	if (world.rank()==0) {
    		calc->molecule.print();
    		print("\n");
    		calc->param.print(world);
    	}
    	double hf_energy = nemo.value();
    	if(world.rank()==0) std::cout << "\n\n\n\n\n\n Reference Calclation Ended\n SCF Energy is: " << hf_energy
    			<<"\n current wall-time: " << wall_time()
    			<<"\n current cpu-time: " << cpu_time()<< "\n\n\n";


    	TDHF tdhf(world,nemo,input);

    	// solve the CIS equations
    	std::vector<CC_vecfunction> roots=tdhf.solve_cis();

    	// analyze the results
    	tdhf.analyze(roots);

    }



	typedef std::vector<real_function_3d> vecfuncT;



	if(world.rank() == 0) printf("\nfinished at time %.1fs\n\n", wall_time());
	world.gop.fence();
	finalize();

	return 0;
}// end main


