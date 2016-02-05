/**
 \file TDA_guess.h
 \brief Helper Function which gives back pre defined guess operators in string form
 */

namespace madness{

struct polynomial_guess{
public:
static std::vector<std::string> make_predefined_guess_strings(const std::string what){
	std::vector<std::string> exop_strings;
	if(what == "dipole"){
		exop_strings.resize(3);
		exop_strings[0] = "x 1.0";
		exop_strings[1] = "y 1.0";
		exop_strings[2] = "z 1.0";
	}else if(what == "quadrupole"){
		exop_strings.resize(9);
		exop_strings[0] = "x 1.0";
		exop_strings[1] = "y 1.0";
		exop_strings[2] = "z 1.0";
		exop_strings[3] = "x 2.0";
		exop_strings[4] = "y 2.0";
		exop_strings[5] = "z 2.0";
		exop_strings[6] = "x 1.0 y 1.0";
		exop_strings[7] = "x 1.0 z 1.0";
		exop_strings[8] = "y 1.0 z 1.0";
	}else if(what == "dipole+"){
		exop_strings.resize(4);
		exop_strings[0] = "x 1.0";
		exop_strings[1] = "y 1.0";
		exop_strings[2] = "z 1.0";
		exop_strings[3] = "x 2.0 , y 2.0 , z 2.0";
	}else if(what == "dipole+diffuse"){
		exop_strings.resize(4);
		exop_strings[0] = "x 3.0 , x 1.0 y 2.0 , x 1.0 z 2.0";
		exop_strings[1] = "x 2.0 y 1.0 , y 3.0 ,y 1.0 z 2.0";
		exop_strings[2] = "x 2.0 z 1.0 , y 2.0 z 1.0 , z 3.0";
		exop_strings[3] = "x 4.0 , 4 2.0 , 4 2.0, x 2.0 y 2.0, x 2.0 z 2.0, y 2.0 z 2.0";
	}else if(what == "dipole+diffuse_big"){
		exop_strings.resize(8);
		exop_strings[0] = "x 1.0";
		exop_strings[1] = "y 1.0";
		exop_strings[2] = "z 1.0";
		exop_strings[3] = "x 2.0 , y 2.0 , z 2.0";
		exop_strings[4] = "x 3.0 , x 1.0 y 2.0 , x 1.0 z 2.0";
		exop_strings[5] = "x 2.0 y 1.0 , y 3.0 ,y 1.0 z 2.0";
		exop_strings[6] = "x 2.0 z 1.0 , y 2.0 z 1.0 , z 3.0";
		exop_strings[7] = "x 4.0 , 4 2.0 , 4 2.0, x 2.0 y 2.0, x 2.0 z 2.0, y 2.0 z 2.0";
	}else if(what == "c2v"){
		exop_strings.resize(4);
		exop_strings[0] = "z 1.0 , z 3.0 , x 2.0 z 1.0 , y 2.0 z 1.0 , x 2.0 , y 2.0 , z 2.0 , x 4.0 , y 4.0 , z 4.0 , x 2.0 y 2.0 , x 2.0 z 2.0 , y 2.0 z 2.0";
		exop_strings[1] = "x 1.0 y 1.0 , x 3.0 y 1.0 , x 1.0 y 3.0 , x 1.0 y 1.0 z 1.0 , x 1.0 y 1.0 z 2.0";
		exop_strings[2] = "x 1.0 , x 1.0 z 1.0 , x 1.0 z 2.0 , x 3.0 , x 3.0 z 1.0 , x 1.0 z 3.0 , x 1.0 y 2.0 , x 1.0 y 2.0 z 1.0";
		exop_strings[3] = "y 1.0 , y 1.0 z 1.0 , y 1.0 z 2.0 , y 3.0 z 1.0 , y 1.0 z 3.0 , y 3.0, x 2.0 y 1.0 , x 2.0 y 1.0 z 1.0 ";
	}else if(what == "c2v_big"){
		exop_strings.resize(8);
		exop_strings[0] = "z 1.0 , z 3.0 , x 2.0 z 1.0 , y 2.0 z 1.0";
		exop_strings[1] = "x 2.0 , y 2.0 , z 2.0 , x 4.0 , y 4.0 , z 4.0 , x 2.0 y 2.0 , x 2.0 z 2.0 , y 2.0 z 2.0";
		exop_strings[2] = "x 1.0 y 1.0 , x 3.0 y 1.0 , x 1.0 y 3.0";
		exop_strings[3] = "x 1.0 y 1.0 z 1.0 , x 1.0 y 1.0 z 2.0";
		exop_strings[4] = "x 1.0 , x 1.0 z 1.0 , x 1.0 z 2.0 , x 3.0 , x 3.0 z 1.0 , x 1.0 z 3.0";
		exop_strings[5] = "x 1.0 y 2.0 , x 1.0 y 2.0 z 1.0";
		exop_strings[6] = "y 1.0 , y 1.0 z 1.0 , y 1.0 z 2.0 , y 3.0 z 1.0 , y 1.0 z 3.0 , y 3.0";
		exop_strings[7] = "x 2.0 y 1.0 , x 2.0 y 1.0 z 1.0";
	}else if(what == "big_fock"){exop_strings = make_auto_polynom_guess(6);
	}else if(what == "small_fock"){exop_strings = make_auto_polynom_guess(4);
	}else if(what == "big_fock_2"){exop_strings = make_auto_polynom_guess(2);
	}else if(what == "big_fock_3"){exop_strings = make_auto_polynom_guess(3);
	}else if(what == "big_fock_4"){exop_strings = make_auto_polynom_guess(4);
	}else if(what == "big_fock_5"){exop_strings = make_auto_polynom_guess(5);
	}else if(what == "big_fock_6"){exop_strings = make_auto_polynom_guess(6);
	}else if(what == "big_fock_7"){exop_strings = make_auto_polynom_guess(7);
	}else if(what == "big_fock_8"){exop_strings = make_auto_polynom_guess(8);
	}else if(what == "big_fock_9"){exop_strings = make_auto_polynom_guess(9);
	}else std::cout << "Keyword " << what << " is not known" << std::endl;
return exop_strings;
}

static std::vector<std::string> make_auto_polynom_guess(const size_t order){
	std::vector<std::string> exop_strings;
	for(size_t i=0; i<order+1; i++){
		for(size_t j=0; j<order+1 ; j++){
			for(size_t k=0;k<order+1 ; k++){
				if(i+j+k > order) ; // do nothing
				else if(i==0 and j==0 and k==0) ; // do nothing
				else{
					if(i==0 and j!=0 and k!=0) exop_strings.push_back(" y " + madness::stringify(j) + " z " + madness::stringify(k) );
					else if(j==0 and i!=0 and k!=0) exop_strings.push_back("x " + madness::stringify(i) + " z " + madness::stringify(k) );
					else if(k==0 and i!=0 and j!=0) exop_strings.push_back("x " + madness::stringify(i) + " y " + madness::stringify(j));
					else if(i==0 and j==0) exop_strings.push_back(" z " + madness::stringify(k) );
					else if(i==0 and k==0) exop_strings.push_back(" y " + madness::stringify(j));
					else if(j==0 and k==0) exop_strings.push_back("x " + madness::stringify(i));
					else exop_strings.push_back("x " + madness::stringify(i) + " y " + madness::stringify(j) + " z " + madness::stringify(k) );
				}
			}
		}
	}
	return exop_strings;
}

};
}
