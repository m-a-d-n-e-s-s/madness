/*
 * CCStructures.cc
 *
 *  Created on: Jan 4, 2017
 *      Author: kottmanj
 */

#include"CCStructures.h"
#include<madness/chem/CCPotentials.h>

namespace madness {

void
CCMessenger::output(const std::string& msg) const {
    if (scientific) std::cout << std::scientific;
    else std::cout << std::fixed;

    std::cout << std::setprecision(output_prec);
    if (world.rank() == 0) std::cout << msg << std::endl;
}

void
CCMessenger::section(const std::string& msg) const {
    if (world.rank() == 0) {
        std::cout << "\n" << std::setw(msg.size() + 10) << std::setfill('*') << "\n";
        std::cout << std::setfill(' ');
        output(msg);
        std::cout << std::setw(msg.size() + 10) << std::setfill('*') << "\n\n";
        std::cout << std::setfill(' ');
    }
}

void
CCMessenger::subsection(const std::string& msg) const {
    if (world.rank() == 0) {
        std::cout << "\n" << std::setw(msg.size() + 5) << std::setfill('-') << "\n";
        std::cout << std::setfill(' ');
        output(msg);
        std::cout << std::setw(msg.size() + 5) << std::setfill('-') << "\n";
        std::cout << std::setfill(' ');
    }
}

void
CCMessenger::warning(const std::string& msg) const {
    std::string tmp = "!!!!!WARNING:" + msg + "!!!!!!";
    output(tmp);
    warnings.push_back(msg);
}

void
CCTimer::info(const bool debug, const double norm) {
    if (debug == true) {
        update_time();
        std::string s_norm = "";
        if (norm != 12345.6789) s_norm = ", ||result||=" + std::to_string(norm);

        if (world.rank() == 0) {
            std::cout << std::setfill(' ') << std::scientific << std::setprecision(2)
                      << "Timer: " << time_wall << " (Wall), " << time_cpu << " (CPU)" << s_norm
                      << ", (" + operation + ")" << "\n";
        }
    }
}



void
CC_vecfunction::print_size(const std::string& msg) const {
    if (functions.size() == 0) {
        std::cout << "CC_vecfunction " << msg << " is empty\n";
    } else {
        std::string msg2;
        if (msg == "!?not assigned!?") msg2 = "";
        else msg2 = "_(" + msg + ")";

        for (auto x : functions) {
            x.second.function.print_size(x.second.name() + msg2);
        }
    }
}

void
CCPair::info() const {
    if (constant_part.world().rank() == 0) {
        std::cout << "\nInformation about electron pair: " << name() << "\n";
    }
    constant_part.print_size("ConstantPart");
    for (size_t k = 0; k < functions.size(); k++)
        functions[k].print_size();
    if (constant_part.world().rank() == 0) {
        std::cout << "\n";
    }
}

madness::vector_real_function_3d
CCIntermediatePotentials::get_potential(const PotentialType& ptype, const FuncType& ftype, const bool throw_if_empty) const {
    vector_real_function_3d result;
    if (parameters.debug()) print("Getting " , assign_name(ptype) , " for " , madness::name(ftype,0));
    if (ptype == POT_singles_ and (ftype == PARTICLE or ftype == MIXED)) result= current_singles_potential_gs_;
    else if (ptype == POT_singles_ and ftype == RESPONSE) result= current_singles_potential_ex_;
    else if (ptype == POT_s2b_ and ftype == PARTICLE) result= current_s2b_potential_gs_;
    else if (ptype == POT_s2b_ and ftype == RESPONSE) result= current_s2b_potential_ex_;
    else if (ptype == POT_s2c_ and ftype == PARTICLE) result= current_s2c_potential_gs_;
    else if (ptype == POT_s2c_ and ftype == RESPONSE) result= current_s2c_potential_ex_;
    else {
        MADNESS_EXCEPTION("unknown potential in CCIntermediatePotentials::get_potential", 1);
    }
    if (result.empty() and throw_if_empty) {
        std::string errmsg="CCIntermediatePotential was not computed/stored "+assign_name(ptype) + " " +assign_name(ftype);
        errmsg+="\n --> you might need to iterate the corresponding singles";
        print(errmsg);
        MADNESS_EXCEPTION(errmsg.c_str(),1);
    }

    if (not result.empty() and parameters.debug()) {
        World& world=result.front().world();
        if (parameters.debug()) print_size(world,result, "potential");
    }
    return result;
}


madness::vector_real_function_3d
CCIntermediatePotentials::operator()(const CC_vecfunction& f, const PotentialType& type,
    const bool throw_if_empty) const {
    return get_potential(type,f.type,throw_if_empty);
}

madness::real_function_3d
CCIntermediatePotentials::operator()(const CCFunction<double,3>& f, const PotentialType& type,
    const bool throw_if_empty) const {
    vector_real_function_3d result=get_potential(type,f.type,throw_if_empty);
    long iact=f.i-parameters.freeze();  // active index
    MADNESS_CHECK_THROW(size_t(iact)<result.size(),"potential not found for active occupied index iact");
    return result[iact];
}

void
CCIntermediatePotentials::insert(const vector_real_function_3d& potential, const CC_vecfunction& f,
                                 const PotentialType& type) {
    World& world=potential.front().world();
    if (world.rank()==0) output("Storing potential: " + assign_name(type) + " for " + f.name(0));
    if (parameters.debug()) {
        print_size(world, potential, "potential");
    }
    MADNESS_ASSERT(!potential.empty());
    if (type == POT_singles_ && (f.type == PARTICLE || f.type == MIXED)) current_singles_potential_gs_ = potential;
    else if (type == POT_singles_ && f.type == RESPONSE) current_singles_potential_ex_ = potential;
    else if (type == POT_s2b_ && f.type == PARTICLE) current_s2b_potential_gs_ = potential;
    else if (type == POT_s2b_ && f.type == RESPONSE) current_s2b_potential_ex_ = potential;
    else if (type == POT_s2c_ && f.type == PARTICLE) {
        current_s2c_potential_gs_ = potential;
    } else if (type == POT_s2c_ && f.type == RESPONSE) {
        current_s2c_potential_ex_ = potential;
    }
}

void CCParameters::set_derived_values() {
    if (not kain()) set_derived_value("kain_subspace",std::size_t(0));

    if (response()==true) set_derived_value("excitations",std::vector<std::size_t>({0}));

    // set all parameters that were not explicitly given
    set_derived_value("tight_thresh_6d",thresh_6D()*0.1);
    set_derived_value("thresh_3d",thresh_6D()*0.01);
    set_derived_value("tight_thresh_3d",thresh_3D()*0.1);
    set_derived_value("thresh_ue",tight_thresh_6D());
    set_derived_value("dconv_6d",3.0*thresh_6D());
    set_derived_value("dconv_3d",0.3*thresh_6D());
    set_derived_value("econv",0.1*dconv_6D());
    set_derived_value("econv_pairs",econv());


    set_derived_value("no_compute_gs",no_compute());
    set_derived_value("no_compute_mp2",no_compute() and no_compute_gs());
    set_derived_value("no_compute_cc2",no_compute() and no_compute_gs());
    set_derived_value("no_compute_cispd",no_compute() and no_compute_response());
    set_derived_value("no_compute_response",no_compute());
    set_derived_value("restart",no_compute() == true and restart() == false);

    if (thresh_3D() < 1.1e-1) set_derived_value("output_prec",std::size_t(3));
    if (thresh_3D() < 1.1e-2) set_derived_value("output_prec",std::size_t(4));
    if (thresh_3D() < 1.1e-3) set_derived_value("output_prec",std::size_t(5));
    if (thresh_3D() < 1.1e-4) set_derived_value("output_prec",std::size_t(6));
    if (thresh_3D() < 1.1e-5) set_derived_value("output_prec",std::size_t(7));
    if (thresh_3D() < 1.1e-6) set_derived_value("output_prec",std::size_t(8));
    std::cout.precision(output_prec());
}

void CCParameters::information(World& world) const {
    if (world.rank()==0) {
//        print("cc2","end");
        if (calc_type() != CT_LRCCS and calc_type() != CT_TDHF) {
            std::cout << "The Ansatz for the Pair functions |tau_ij> is: ";
            if (QtAnsatz()) std::cout << "(Qt)f12|titj> and response: (Qt)f12(|tixj> + |xitj>) - (OxQt + QtOx)f12|titj>";
            else std::cout << "Qf12|titj> and response: Qf12(|xitj> + |tixj>)" << std::endl;
        }
    }
}

void CCParameters::sanity_check(World& world) const {
    size_t warnings = 0;
    if (FunctionDefaults<3>::get_thresh() > 0.01 * FunctionDefaults<6>::get_thresh())
        warnings += warning(world, "3D Thresh is too low, should be 0.01*6D_thresh");
    if (FunctionDefaults<3>::get_thresh() > 0.1 * FunctionDefaults<6>::get_thresh())
        warnings += warning(world, "3D Thresh is way too low, should be 0.01*6D_thresh");
    if (FunctionDefaults<3>::get_cell_min_width() != FunctionDefaults<6>::get_cell_min_width())
        warnings += warning(world, "3D and 6D Cell sizes differ");
    if (FunctionDefaults<3>::get_k() != FunctionDefaults<6>::get_k())
        warnings += warning(world, "k-values of 3D and 6D differ ");
    if (FunctionDefaults<3>::get_truncate_mode() != 3) warnings += warning(world, "3D Truncate mode is not 3");
    if (FunctionDefaults<6>::get_truncate_mode() != 3) warnings += warning(world, "6D Truncate mode is not 3");
    if (dconv_3D() < FunctionDefaults<3>::get_thresh())
        warnings += warning(world, "Demanded higher convergence than threshold for 3D");
    if (dconv_6D() < FunctionDefaults<6>::get_thresh())
        warnings += warning(world, "Demanded higher convergence than threshold for 6D");
    if (thresh_3D() != FunctionDefaults<3>::get_thresh())
        warnings += warning(world, "3D thresh set unequal 3D thresh demanded");
    if (thresh_6D() != FunctionDefaults<6>::get_thresh())
        warnings += warning(world, "6D thresh set unequal 6D thresh demanded");
    if (econv() < FunctionDefaults<3>::get_thresh())
        warnings += warning(world, "Demanded higher energy convergence than threshold for 3D");
    if (econv() < FunctionDefaults<6>::get_thresh())
        warnings += warning(world, "Demanded higher energy convergence than threshold for 6D");
    if (econv() < 0.1 * FunctionDefaults<3>::get_thresh())
        warnings += warning(world,
                            "Demanded higher energy convergence than threshold for 3D (more than factor 10 difference)");
    if (econv() < 0.1 * FunctionDefaults<6>::get_thresh())
        warnings += warning(world,
                            "Demanded higher energy convergence than threshold for 6D (more than factor 10 difference)");
    // Check if the 6D thresholds are not too high
    if (thresh_6D() < 1.e-3) warnings += warning(world, "thresh_6D is smaller than 1.e-3");
    if (thresh_6D() < tight_thresh_6D()) warnings += warning(world, "tight_thresh_6D is larger than thresh_6D");
    if (thresh_6D() < tight_thresh_3D()) warnings += warning(world, "tight_thresh_3D is larger than thresh_3D");
    if (thresh_6D() < 1.e-3) warnings += warning(world, "thresh_6D is smaller than 1.e-3");
    if (thresh_Ue() < 1.e-4) warnings += warning(world, "thresh_Ue is smaller than 1.e-4");
    if (thresh_Ue() > 1.e-4) warnings += warning(world, "thresh_Ue is larger than 1.e-4");
    if (thresh_3D() > 0.01 * thresh_6D())
        warnings += warning(world, "Demanded 6D thresh is to precise compared with the 3D thresh");
    if (thresh_3D() > 0.1 * thresh_6D())
        warnings += warning(world, "Demanded 6D thresh is to precise compared with the 3D thresh");
    if (kain() and kain_subspace() == 0)
        warnings += warning(world, "Demanded Kain solver but the size of the iterative subspace is set to zero");
    if (warnings > 0) {
        if (world.rank() == 0) std::cout << warnings << "Warnings in parameters sanity check!\n\n";
    } else {
        if (world.rank() == 0) std::cout << "Sanity check for parameters passed\n\n" << std::endl;
    }
    if (restart() == false and no_compute() == true) {
        warnings += warning(world, "no_compute flag detected but no restart flag");
    }
}

template<typename T, std::size_t NDIM>
Function<T,NDIM>
CCConvolutionOperator<T,NDIM>::operator()(const CCFunction<T,NDIM>& bra, const CCFunction<T,NDIM>& ket, const bool use_im) const {
    Function<T,NDIM> result;
    if (not use_im) {
        if (world.rank() == 0)
            std::cout << "Recalculating <" << bra.name() << "|" << name() << "|" << ket.name()
                      << ">\n";
        result = ((*op)(bra.function * ket.function)).truncate();
    } else if (bra.type == HOLE and ket.type == HOLE and not imH.allpairs.empty()) result = imH(bra.i, ket.i);
    else if (bra.type == HOLE and ket.type == RESPONSE and not imR.allpairs.empty()) result = imR(bra.i, ket.i);
    else if (bra.type == HOLE and ket.type == PARTICLE and not imP.allpairs.empty()) result = imP(bra.i, ket.i);
    else if (bra.type == HOLE and ket.type == MIXED and (not imP.allpairs.empty() and not imH.allpairs.empty()))
        result = (imH(bra.i, ket.i) + imP(bra.i, ket.i));
    else {
        //if(world.rank()==0) std::cout <<"No Intermediate found for <" << bra.name()<<"|"<<assign_name(operator_type) <<"|"<<ket.name() <<"> ... recalculate \n";
        MADNESS_ASSERT(op);
        result = ((*op)(bra.function * ket.function)).truncate();
    }
    return result;
}

template<typename T, std::size_t NDIM>
Function<T,2*NDIM> CCConvolutionOperator<T,NDIM>::operator()(const Function<T,2*NDIM>& u, const size_t particle) const {
    MADNESS_CHECK(particle == 1 or particle == 2);
    MADNESS_CHECK(op);
    op->particle() = particle;
    return (*op)(u);
}

template<typename T, std::size_t NDIM>
Function<T,NDIM>
CCConvolutionOperator<T,NDIM>::operator()(const CCFunction<T,NDIM>& bra, const Function<T,2*NDIM>& u, const size_t particle) const {
    MADNESS_CHECK(particle == 1 or particle == 2);
    MADNESS_CHECK(op);
    const Function<T,2*NDIM> tmp = multiply(copy(u), copy(bra.function), particle);
    op->particle() = particle;
    const Function<T,2*NDIM> g_tmp = (*op)(tmp);
    const Function<T,NDIM> result = g_tmp.dirac_convolution();
    return result;
}

template<typename T, std::size_t NDIM>
void CCConvolutionOperator<T,NDIM>::update_elements(const CC_vecfunction& bra, const CC_vecfunction& ket) {
    if constexpr (NDIM==3) {
        const std::string operation_name = "<" + assign_name(bra.type) + "|" + name() + "|" + assign_name(ket.type) + ">";
        if (world.rank() == 0)
            std::cout << "updating operator elements: " << operation_name << " (" << bra.size() << "x" << ket.size() << ")"
                      << std::endl;
        if (bra.type != HOLE)
            error("Can not create intermediate of type " + operation_name + " , bra-element has to be of type HOLE");
        op.reset(init_op(type(), parameters));
        intermediateT<T,NDIM> xim;
        for (auto tmpk : bra.functions) {
            const CCFunction<T,NDIM>& k = tmpk.second;
            for (auto tmpl : ket.functions) {
                const CCFunction<T,NDIM>& l = tmpl.second;
                Function<T,NDIM> kl = (bra(k).function * l.function);
                Function<T,NDIM> result = ((*op)(kl)).truncate();
                result.reconstruct(); // for sparse multiplication
                xim.insert(k.i, l.i, result);
            }
        }
        if (ket.type == HOLE) imH = xim;
        else if (ket.type == PARTICLE) imP = xim;
        else if (ket.type == RESPONSE) imR = xim;
        else error("Can not create intermediate of type <" + assign_name(bra.type) + "|op|" + assign_name(ket.type) + ">");
    } else {
        std::string msg="update_elements not implemented for NDIM="+std::to_string(NDIM);
        MADNESS_EXCEPTION(msg.c_str(),1);
    }
}


template<typename T, std::size_t NDIM>
void CCConvolutionOperator<T,NDIM>::clear_intermediates(const FuncType& type) {
    if (world.rank() == 0)
        std::cout << "Deleting all <HOLE|" << name() << "|" << assign_name(type) << "> intermediates \n";
    switch (type) {
        case HOLE : {
            imH.allpairs.clear();
            break;
        }
        case PARTICLE: {
            imP.allpairs.clear();
            break;
        }
        case RESPONSE: {
            imR.allpairs.clear();
            break;
        }
        default:
            error("intermediates for " + assign_name(type) + " are not defined");
    }
}

template<typename T, std::size_t NDIM>
size_t CCConvolutionOperator<T,NDIM>::info() const {
    const size_t size_imH = size_of(imH);
    const size_t size_imP = size_of(imP);
    const size_t size_imR = size_of(imR);
    if (world.rank() == 0) {
        std::cout << "Size of " << name() << " intermediates:\n";
        std::cout << std::setw(5) << "(" << imH.allpairs.size() << ") x <H|" + name() + "H>=" << std::scientific
                  << std::setprecision(1) << size_imH << " (Gbyte)\n";
        std::cout << std::setw(5) << "(" << imP.allpairs.size() << ") x <H|" + name() + "P>=" << std::scientific
                  << std::setprecision(1) << size_imH << " (Gbyte)\n";
        std::cout << std::setw(5) << "(" << imR.allpairs.size() << ") x <H|" + name() + "R>=" << std::scientific
                  << std::setprecision(1) << size_imH << " (Gbyte)\n";
    }
    return size_imH + size_imP + size_imR;
}

template<typename T, std::size_t NDIM>
SeparatedConvolution<T, NDIM> *
CCConvolutionOperator<T,NDIM>::init_op(const OpType& type, const CCConvolutionOperator<T,NDIM>::Parameters& parameters) const {
    bool debug=false;
    bool printme=(world.rank()==0) and debug;
    if (printme) print("init_op: creating",type,"with thresh, lo, gamma",parameters.thresh_op,parameters.lo,parameters.gamma);
    return new SeparatedConvolution<T,NDIM>(world,OperatorInfo(parameters.gamma,parameters.lo,parameters.thresh_op,type));
}

/// Assigns strings to enums for formated output
std::string
assign_name(const CCState& input) {
    switch (input) {
        case GROUND_STATE:
            return "Ground State";
        case EXCITED_STATE:
            return "Excited State";
        default: {
            MADNESS_EXCEPTION("Unvalid enum assignement!", 1);
            return "undefined";
        }
    }
    MADNESS_EXCEPTION("assign_name:pairtype, should not end up here", 1);
    return "unknown pairtype";
}

/// Assigns enum to string
CalcType
assign_calctype(const std::string name) {
    if (name == "mp2") return CT_MP2;
    else if (name == "cc2") return CT_CC2;
    else if (name == "lrcc2" or name == "cc2_response") return CT_LRCC2;
    else if (name == "cispd") return CT_CISPD;
    else if (name == "cis" or name == "ccs" or name == "ccs_response" or name == "lrccs") return CT_LRCCS;
    else if (name == "experimental") return CT_TEST;
    else if (name == "adc2" or name == "adc(2)") return CT_ADC2;
    else if (name == "tdhf") return CT_TDHF;
    else {
        std::string msg = "CALCULATION OF TYPE: " + name + " IS NOT KNOWN!!!!";
        MADNESS_EXCEPTION(msg.c_str(), 1);
    }
}

/// Assigns strings to enums for formated output
std::string
assign_name(const CalcType& inp) {
    switch (inp) {
        case CT_CC2:
            return "CC2";
        case CT_MP2:
            return "MP2";
        case CT_LRCC2:
            return "LRCC2";
        case CT_CISPD:
            return "CISpD";
        case CT_LRCCS:
            return "LRCCS";
        case CT_ADC2:
            return "ADC2";
        case CT_TDHF:
            return "TDHF";
        case CT_TEST:
            return "experimental";
        default: {
            MADNESS_EXCEPTION("Unvalid enum assignement!", 1);
            return "undefined";
        }
    }
    return "unknown";
}

/// Assigns strings to enums for formated output
std::string
assign_name(const PotentialType& inp) {
    switch (inp) {
        case POT_F3D_:
            return "F3D";
        case POT_s3a_:
            return "s3a";
        case POT_s3b_:
            return "s3b";
        case POT_s3c_:
            return "s3c";
        case POT_s5a_:
            return "s5a";
        case POT_s5b_:
            return "s5b";
        case POT_s5c_:
            return "s5c";
        case POT_s6_:
            return "s6";
        case POT_s2b_:
            return "s2b";
        case POT_s2c_:
            return "s2c";
        case POT_s4a_:
            return "s4a";
        case POT_s4b_:
            return "s4b";
        case POT_s4c_:
            return "s4c";
        case POT_ccs_:
            return "ccs";
        case POT_cis_:
            return "cis-potential";
        case POT_singles_:
            return "singles potential";
        default: {
            MADNESS_EXCEPTION("Unvalid enum assignement!", 1);
            return "undefined";
        }
    }
    return "undefined";
}

/// Assigns strings to enums for formated output
std::string
assign_name(const FuncType& inp) {
    switch (inp) {
        case HOLE:
            return "Hole";
        case PARTICLE:
            return "Particle";
        case MIXED:
            return "Mixed";
        case RESPONSE:
            return "Response";
        case UNDEFINED:
            return "Undefined";
        default: {
            MADNESS_EXCEPTION("Unvalid enum assignement!", 1);
            return "undefined";
        }
    }
    return "???";
}

/// make a CCPair without the 6d function and some bookkeeping information
CCPair CCPairBuilder::make_bare_pair(const int i, const int j) const {
    MADNESS_ASSERT(i>=info.parameters.freeze() && size_t(i) < info.mo_bra.size());
    MADNESS_ASSERT(j>=info.parameters.freeze() && size_t(j) < info.mo_ket.size());

    CCPair pair(i, j, cc_state(ctype), ctype);
    pair.bsh_eps=CCPotentials::get_epsilon(i,j,info);
    if (cc_state(ctype)==EXCITED_STATE) {
        MADNESS_ASSERT(ex_singles.omega != 0.0);
        pair.bsh_eps += ex_singles.omega;
    }
    return pair;
}

/// make a CCPair with the given function
CCPair CCPairBuilder::make_pair(const int i, const int j, const std::vector<CCPairFunction<double,6>>& u) const {
    CCPair pair=make_bare_pair(i,j);
    // a lot of logic depends on the first function being the 6d function!
    if (u.size()>0) MADNESS_CHECK_THROW(u.front().is_pure(),"missing pure 6d function in CCPairBuilder::make_pair");
    pair.functions+=u;
    return pair;
}

/// make a CCPair with the 6d function only and some bookkeeping information
CCPair CCPairBuilder::make_bare_pair_from_file(const int i, const int j) const {
    CCPair pair=make_bare_pair(i,j);

    // load the 6d function u and the constant part from file
    std::string name = pair.name();
    real_function_6d utmp=load_function<double,6>(name, info.parameters.debug());
    real_function_6d const_part=load_function<double,6>(name + "_const", info.parameters.debug());

    // first term is the 6d function u, then follows Q12 f12 |ij>, which is added later
    if (utmp.is_initialized()) pair.functions+=CCPairFunction<double,6>(utmp);
    if (const_part.is_initialized()) pair.constant_part = const_part;

    return pair;
}


CCPair CCPairBuilder::complete_pair_with_low_rank_parts(const CCPair& pair) const {
    CCPair result=pair;

    // a lot of logic depends on the first function being the 6d function!
//    if (result.functions.size()==0) {
//        real_function_6d f;
//        result.functions.push_back(CCPairFunction<double,6>(f));
//    }
    MADNESS_CHECK_THROW(result.functions.size()==1,"missing pure 6d function in CCPairBuilder::complete_pair_with_low_rank_parts");
    MADNESS_CHECK_THROW(result.functions.front().is_pure(),"pure 6d function not pure in CCPairBuilder::complete_pair_with_low_rank_parts");
    long nact=info.get_active_mo_bra().size();
    if (result.ctype==CT_MP2) {
        // nothing to do
    } else if (result.ctype==CT_CC2) {
        MADNESS_CHECK_THROW(gs_singles.size()==size_t(nact),"missing gs_singles for completing the CC2 pair function");
    } else if (result.ctype==CT_LRCC2) {
        MADNESS_CHECK_THROW(gs_singles.size()==size_t(nact),"missing gs_singles for completing the LRCC2 pair function");
        MADNESS_CHECK_THROW(ex_singles.size()==size_t(nact),"missing ex_singles for completing the LRCC2 pair function");
    } else {
        print("unknown ctype in complete_pair_with_low_rank_parts",assign_name(result.ctype));
        MADNESS_EXCEPTION("unknown ctype",1);
    }

    timer t1(world);
    auto phi=info.mo_ket;
    auto phi_bra=info.mo_bra;
    StrongOrthogonalityProjector<double,3> Q12(world);
    auto f12=CCConvolutionOperatorPtr<double,3>(world,OT_F12,info.parameters);

    if (result.ctype==CT_MP2) {
        // ansatz is Q12 f12 |ij>
        Q12.set_spaces(phi_bra,phi,phi_bra,phi);
        CCPairFunction<double,6> fij(f12, phi[result.i], phi[result.j]);
        std::vector<CCPairFunction<double,6>> tmp=Q12(std::vector<CCPairFunction<double,6>>(1,fij));
        result.functions+=tmp;

    } else if (result.ctype==CT_CC2) {
        // ansatz is Qt12 f12 |t_i t_j>
        auto t=CCPotentials::make_full_t_intermediate(gs_singles,info).get_vecfunction();
        Q12.set_spaces(phi_bra,t,phi_bra,t);

        CCPairFunction<double,6> fij(f12, t[result.i], t[result.j]);
        std::vector<CCPairFunction<double,6>> tmp=Q12(std::vector<CCPairFunction<double,6>>(1,fij));
        result.functions+=tmp;
    } else if (result.ctype==CT_LRCC2) {
        // ansatz is Qt12 f12 ( |t_i x_j> + |x_i t_j> ) - dQt12 f12 |t_i t_j>
        result=CCPotentials::make_pair_lrcc2(world,result.ctype,result.function(),gs_singles,ex_singles,result.i,result.j,info, true);
    } else {
        MADNESS_EXCEPTION("unknown ctype",1);
    }
    t1.tag("make low-rank parts in make_pair_cc2 for pair("+stringify(result.i)+stringify(result.j)+")");
    return result;
}

std::vector<real_function_6d>
//MacroTaskMp2ConstantPart::operator() (const std::vector<CCPair>& pair, const std::vector<real_function_3d>& mo_ket,
//                                      const std::vector<real_function_3d>& mo_bra, const CCParameters& parameters,
//                                      const real_function_3d& Rsquare, const std::vector<real_function_3d>& U1,
//                                      const std::vector<std::string>& argument) const {
MacroTaskMp2ConstantPart::operator() (const std::vector<CCPair>& pair, const Info& info,
                                      const std::vector<std::string>& argument) const {
    World& world =info.mo_ket[0].world();
    resultT result = zero_functions_compressed<double, 6>(world, pair.size());
    for (size_t i = 0; i < pair.size(); i++) {
        result[i] = CCPotentials::make_constant_part_mp2_macrotask(world, pair[i], info.mo_ket, info.mo_bra,
                                    info.parameters, info.R_square, info.U1, argument);
    }
    return result;
}

std::vector<real_function_6d>
MacroTaskConstantPart::operator() (const std::vector<CCPair>& pair,
                                   const std::vector<Function<double,3>> & gs_singles,
                                   const std::vector<Function<double,3>> & ex_singles,
                                   const Info& info) const {

    World& world =info.mo_ket[0].world();
    CC_vecfunction singles(gs_singles, PARTICLE, info.parameters.freeze());
    CC_vecfunction exsingles(ex_singles, RESPONSE, info.parameters.freeze());


    resultT result = zero_functions_compressed<double, 6>(world, pair.size());
    for (size_t i = 0; i < pair.size(); i++) {
        result[i] = CCPotentials::make_constant_part_macrotask(world, pair[i], singles, exsingles, info);
    }
    return result;
}


std::vector<real_function_6d>
MacroTaskMp2UpdatePair::operator() (const std::vector<CCPair> &pair,
                                    const std::vector<real_function_6d> &mp2_coupling,
                                    const std::vector<madness::Vector<double, 3>> &all_coords_vec,
                                    const Info& info) const {
    World& world = info.mo_ket[0].world();
    resultT result = zero_functions_compressed<double, 6>(world, pair.size());
    print("in MacroTaskMp2UpdatePair::operator()", "pair.size()=", pair.size(), "batch=", this->batch);

    for (size_t i = 0; i < pair.size(); i++) {
        print("in loop of batch",this->batch);
        //(i, j) -> j*(j+1) + i
        result[i] = CCPotentials::update_pair_mp2_macrotask(world, pair[i], info.parameters, all_coords_vec, info.mo_ket,
                                                            info.mo_bra, info.U1, info.U2, mp2_coupling[i]);
    }
    return result;
}

std::vector<real_function_6d>
MacroTaskIteratePair::operator()(const std::vector<CCPair>& pair,
        const std::vector<real_function_6d>& local_coupling,
        const CC_vecfunction& gs_singles,
        const CC_vecfunction& ex_singles,
        const Info& info,
        const std::size_t& maxiter) const {
    World& world = info.mo_ket[0].world();
    resultT result = zero_functions_compressed<double, 6>(world, pair.size());

    for (size_t i = 0; i < pair.size(); i++) {
        result[i]=  CCPotentials::iterate_pair_macrotask(world, pair[i], gs_singles, ex_singles,
            local_coupling[i], info, maxiter).function();
    }
    return result;

}

/// convenience function


std::tuple<std::vector<real_function_3d>, std::vector<real_function_3d>>
MacroTaskSinglesPotentialEx::operator()(const std::vector<int>& result_index,
                                      const CC_vecfunction& singles_gs,
                                      const std::vector<CCPair>& doubles_gs,
                                      const CC_vecfunction& singles_ex,
                                      const std::vector<CCPair>& doubles_ex,
                                      const int& name,
                                      const Info& info) {
    World& world=singles_ex.get_vecfunction().front().world();

    auto triangular_map=PairVectorMap::triangular_map(info.parameters.freeze(),info.mo_ket.size());
    auto doubles_gs1=Pairs<CCPair>::vector2pairs(doubles_gs,triangular_map);
    auto doubles_ex1=Pairs<CCPair>::vector2pairs(doubles_ex,triangular_map);

//    // the doubles currently only contain the full 6d function -> complete it with the Q12 f12 |ti tj> part
//    for (auto& x : doubles_gs1.allpairs) {
//        auto& tau=x.second;
//        MADNESS_CHECK_THROW(tau.functions.size()==1,"doubles in MacroTaskSinglesPotentialsEx should only contain one function");
//        bool compute_Q12_F12=(PotentialType(name)==POT_s2b_ or PotentialType(name)==POT_s2c_);
//        x.second=CCPotentials::make_pair_cc2(world,tau.function(),singles_gs,tau.i,tau.j, info, compute_Q12_F12);
//    }
//    // the doubles currently only contain the full 6d function -> complete it with the Q12 f12 |ti tj> part
//    for (auto& x : doubles_ex1.allpairs) {
//        auto& tau=x.second;
//        MADNESS_CHECK_THROW(tau.functions.size()==1,"doubles in MacroTaskSinglesPotentialsEx should only contain one function");
//        x.second=tau=CCPotentials::make_pair_lrcc2(world,tau.ctype,tau.function(),singles_gs,singles_ex,tau.i,tau.j, info, true);
//    }

    resultT result=CCPotentials::potential_singles_ex(world,
                result_index,
                singles_gs,
                doubles_gs1,
                singles_ex,
                doubles_ex1,
                PotentialType(name),
                info);
    // if the second element of the tuple is empty, fill it with empty functions
    // to that "insert_batch" is not confused
    if (std::get<1>(result).empty()) std::get<1>(result)=zero_functions<double,3>(world,result_index.size());
    return result;
}

std::tuple<std::vector<real_function_3d>, std::vector<real_function_3d>>
MacroTaskSinglesPotentialGs::operator()(const std::vector<int>& result_index,
                                      const CC_vecfunction& singles_gs,
                                      const std::vector<CCPair>& doubles_gs,
                                      const int& name,
                                      const Info& info) {
    World& world=singles_gs.get_vecfunction().front().world();
    auto triangular_map=PairVectorMap::triangular_map(info.parameters.freeze(),info.mo_ket.size());
    auto doubles_gs1=Pairs<CCPair>::vector2pairs(doubles_gs,triangular_map);

//    // the doubles currently only contain the full 6d function -> complete it with the Q12 f12 |ti tj> part
//    for (auto& x : doubles_gs1.allpairs) {
//        auto& tau=x.second;
//        MADNESS_CHECK_THROW(tau.functions.size()==1,"doubles in MacroTaskSinglesPotentialsGS should only contain one function");
//        bool compute_Q12_F12=(PotentialType(name)==POT_s2b_ or PotentialType(name)==POT_s2c_);
//        tau=CCPotentials::make_pair_cc2(world,tau.function(),singles_gs,tau.i,tau.j,info, compute_Q12_F12);
//    }

    resultT result=CCPotentials::potential_singles_gs(world, result_index,
                singles_gs, doubles_gs1, PotentialType(name), info);

    // if the second element of the tuple is empty, fill it with empty functions
    // to that "insert_batch" is not confused
    auto& intermediate=std::get<1>(result);
    if (intermediate.empty()) intermediate=zero_functions<double,3>(world,result_index.size());
    // if (intermediate.empty()) intermediate.resize(result_index.size());
    MADNESS_CHECK_THROW(std::get<0>(result).size()==std::get<1>(result).size(),"result size mismatch 1 in MacroTaskSinglesPotentialGS");
    MADNESS_CHECK_THROW(std::get<0>(result).size()==result_index.size(),"result size mismatch 2 in MacroTaskSinglesPotentialGS");
    return result;

}

std::vector<ScalarResult<double>>
MacroTaskComputeCorrelationEnergy::operator()(const std::vector<CCPair>& pairs,
                                                  const CC_vecfunction& singles_gs,
                                                  const Info& info) const {
     World &world = pairs[0].function().world();
     auto result=scalar_result_vector<double>(world,pairs.size());
     CalcType ctype=pairs[0].ctype;
     for (size_t i=0; i<pairs.size(); ++i) {
         if (ctype==CT_MP2) {
             // when serialized the Qf12 |ij> part is not stored in the cloud, so recompute it here
             auto pair=CCPotentials::make_pair_mp2(world,pairs[i].function(),pairs[i].i,pairs[i].j,info, true);
             result[i]=CCPotentials::compute_pair_correlation_energy(world,pair,singles_gs,info);
        } else if (ctype==CT_CC2) {
             auto pair=CCPotentials::make_pair_cc2(world,pairs[i].function(),singles_gs,pairs[i].i,pairs[i].j,info, true);
             result[i]=CCPotentials::compute_pair_correlation_energy(world,pair,singles_gs,info);
        } else {
             MADNESS_EXCEPTION("MacroTaskComputeCorrelationEnergy: unknown ctype",1);
        }
     }
    return result;
}

template class CCConvolutionOperator<double,3>;
template class CCConvolutionOperator<double,2>;
template class CCConvolutionOperator<double,1>;

template class CCFunction<double,3>;
template class CCFunction<double,2>;
template class CCFunction<double,1>;

}// end namespace madness


