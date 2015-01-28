/// Helper Function which gives back pre defined guess operators in string form

static std::vector<std::string> make_predefined_guess_strings(const std::string what){
	std::vector<std::string> exop_strings;
	if(what == "dipole"){
		exop_strings.resize(3);
		exop_strings[0] = "x 1.0";
		exop_strings[1] = "y 1.0";
		exop_strings[2] = "z 1.0";
	}else if(what == "quadrupole"){
		exop_strings.resize(1);
		exop_strings[0] = "x 1.0";
		exop_strings[0] = "y 1.0";
		exop_strings[0] = "z 1.0";
		exop_strings[0] = "x 2.0";
		exop_strings[0] = "y 2.0";
		exop_strings[0] = "z 2.0";
		exop_strings[0] = "x 1.0 y 1.0";
		exop_strings[0] = "x 1.0 z 1.0";
		exop_strings[0] = "y 1.0 z 1.0";
	}else if(what == "dipole+"){
		exop_strings.resize(4);
		exop_strings[0] = "x 1.0";
		exop_strings[1] = "y 1.0";
		exop_strings[2] = "z 1.0";
		exop_strings[3] = "x 2.0 , y 2.0 , z 2.0";
	}else if(what == "c2v"){
		exop_strings.resize(4);
		exop_strings[0] = "z 1.0 , z 3.0 , x 2.0 z 1.0 , y 2.0 z 1.0 , x 2.0 , y 2.0 , z 2.0 , x 4.0 , y 4.0 , z 4.0 , x 2.0 y 2.0 , x 2.0 z 2.0 , y 2.0 z 2.0";
		exop_strings[1] = "x 1.0 y 1.0 , x 3.0 y 1.0 , x 1.0 y 3.0 , x 1.0 y 1.0 z 1.0 , x 1.0 y 1.0 z 2.0";
		exop_strings[2] = "x 1.0 , x 1.0 z 1.0 , x 1.0 z 2.0 , x 3.0 , x 3.0 z 1.0 , x 1.0 z 3.0 , x 1.0 y 2.0 , x 1.0 y 2.0 z 1.0";
		exop_strings[3] = "y 1.0 , y 1.0 z 1.0 , y 1.0 z 2.0 , y 3.0 z 1.0 , y 1.0 z 3.0 , y 3.0, x 2.0 y 1.0 , x 2.0 y 1.0 z 1.0 ";
	}else if(what == "c2v_big"){
		exop_strings.resize(8);
		exop_strings[1] = "z 1.0 , z 3.0 , x 2.0 z 1.0 , y 2.0 z 1.0";
		exop_strings[2] = "x 2.0 , y 2.0 , z 2.0 , x 4.0 , y 4.0 , z 4.0 , x 2.0 y 2.0 , x 2.0 z 2.0 , y 2.0 z 2.0";
		exop_strings[3] = "x 1.0 y 1.0 , x 3.0 y 1.0 , x 1.0 y 3.0";
		exop_strings[4] = "x 1.0 y 1.0 z 1.0 , x 1.0 y 1.0 z 2.0";
		exop_strings[5] = "x 1.0 , x 1.0 z 1.0 , x 1.0 z 2.0 , x 3.0 , x 3.0 z 1.0 , x 1.0 z 3.0";
		exop_strings[6] = "x 1.0 y 2.0 , x 1.0 y 2.0 z 1.0";
		exop_strings[7] = "y 1.0 , y 1.0 z 1.0 , y 1.0 z 2.0 , y 3.0 z 1.0 , y 1.0 z 3.0 , y 3.0";
		exop_strings[8] = "x 2.0 y 1.0 , x 2.0 y 1.0 z 1.0";
	} else std::cout << "Keyword " << what << " is not known" << std::endl;
return exop_strings;
}
