//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include "version.h"
#include <apps/chem/SCF.h>
#include <apps/chem/nemo.h>
#include "PNO.h"

using namespace madness;

// DEFINE PARAMETER TAGS
const std::string TAG_PNO = "pno";
const std::string TAG_F12 = "f12";
const std::string TAG_CP = "computeprotocol";

int main(int argc, char** argv) {
	initialize(argc, argv);
	World world(SafeMPI::COMM_WORLD);
	if (world.rank() == 0) printf("starting at time %.1f\n", wall_time());

	if(world.rank()==0) std::cout << "-------------MRA-PNO-MP2-F12---------------\n";
	if(world.rank()==0) std::cout << "Git Hash           :"<< GIT_COMMIT_HASH << "\n";
	if(world.rank()==0) std::cout << "Git Branch         :" << GIT_BRANCH << "\n";
	//if(world.rank()==0) std::cout << "MAD_ROOT_DIR       :" << MADNESS_ROOT_DIR << "\n";
	if(world.rank()==0) std::cout << "CXX_FLAGS=         :" << CXX_FLAGS << "\n";
	if(world.rank()==0) std::cout << "C_FLAGS            :" << C_FLAGS << "\n";
	if(world.rank()==0) std::cout << "RR_CHOLESKY        :" << RR_CHOLESKY << "\n";
	if(world.rank()==0) std::cout << "configured at      :" << CONFIGURE_DATE << "\n";
	if(world.rank()==0) std::cout << "--------------------------------------------\n";



	startup(world,argc,argv,true);
	print_meminfo(world.rank(), "startup");
	std::cout.precision(6);

	const std::string input = (argc > 1) ? argv[1] : "input";
	std::shared_ptr<SCF> calc(new SCF(world, input));
	Nemo nemo(world, calc, input);
	nemo.get_calc()->param.print();
	const double energy = nemo.value();
	if (world.rank() == 0) print("nemo energy: ", energy);
	if (world.rank() == 0) printf(" at time %.1f\n", wall_time());

	const vector_real_function_3d nemos = nemo.get_calc()->amo;
	const vector_real_function_3d R2nemos =
			mul(world, nemo.nuclear_correlation->square(), nemos);


	ComputeProtocol protocol;
	try{
		protocol=ComputeProtocol(world,input,TAG_CP);
	}catch(...){
		if(world.rank()==0) std::cout << "compute protocol not given\n...using defaults";
	}

	if(protocol.thresh().empty()){
		if(world.rank()==0) std::cout << "no compute protocol found --- You have to set the keywords for the PNO class yourself\n";
		PNOParameters parameters(world,input,TAG_PNO);
		F12Parameters paramf12(world, input, parameters, TAG_F12);
		PNO pno(world, nemo, parameters, paramf12);
		pno.solve();
	}else{
		if(world.rank()==0){
			std::cout << "Compute Protocol is:\n";
		}
		protocol.print(TAG_CP,"end");

		const auto amos=madness::copy(world,nemo.get_calc()->amo);
		print_size(world,nemo.get_calc()->amo,"mos original ");
		std::vector<PNOPairs> all_pairs;
		PNOParameters parameters; // default parameters
		try{
			parameters = PNOParameters(world,input, TAG_PNO);
		}catch(...){
			if(world.rank()==0) std::cout << "no pno input found: using defaults\n";
		}

		for(size_t i=0;i<protocol.thresh().size();++i){
			const double thresh=protocol.thresh()[i];
			if(world.rank()==0) std::cout << "Invoce Protocol for thresh="<< thresh << "\n";

			parameters.set_user_defined_value("thresh", thresh);
			parameters.set_user_defined_value("op_thresh", 0.1*thresh);
			if(protocol.op_thresh().size()>i) parameters.set_user_defined_value("op_thresh", protocol.op_thresh()[i]);
			if(protocol.econv_micro().size()>i) parameters.set_user_defined_value("econv_micro", protocol.econv_micro()[i]);
			if(protocol.econv_macro().size()>i) parameters.set_user_defined_value("econv_macro", protocol.econv_macro()[i]);
			if(protocol.maxiter_micro().size()>i) parameters.set_user_defined_value("maxiter_micro", protocol.maxiter_micro()[i]);
			if(protocol.maxiter_macro().size()>i) parameters.set_user_defined_value("maxiter_macro", protocol.maxiter_macro()[i]);
			if(protocol.exop().size()>i) parameters.set_user_defined_value("exop", protocol.exop()[i]);


			vector_real_function_3d amos_tr=madness::copy(world,amos);
			madness::truncate(world,amos_tr,thresh);
			nemo.get_calc()->amo=amos_tr;
			nemo.get_calc()->set_protocol<3>(world,thresh);
			print_size(world,nemo.get_calc()->amo,"mos truncated");
			nemo.get_calc()->make_nuclear_potential(world);

			if(world.rank()==0) std::cout << "\n\n\n" << "STARTING COMPUTE PROTOCOL WITH THRESH=" << thresh << "\n\n\n";

			F12Parameters paramf12(parameters); // default parameters
			try{
				paramf12 = F12Parameters(world,input, parameters, TAG_PNO);
			}catch(...){
				if(world.rank()==0) std::cout << "no f12 input found: using defaults\n";
			}

			PNO pno(world, nemo, parameters, paramf12);
			pno.solve(all_pairs);
			if(world.rank()==0) std::cout << "\n\n\n" << "ENDING COMPUTE PROTOCOL WITH THRESH=" << thresh << " at time " << wall_time() << "\n\n\n";
		}
	}



	world.gop.fence();
	if (world.rank() == 0) printf("finished at time %.1f\n", wall_time());
	world.gop.fence();
	print_stats(world);
	world.gop.fence();
	finalize();
	return 0;
}
