/*
 * CCOperators.cc
 *
 *  Created on: Jul 6, 2015
 *      Author: kottmanj
 */

#include "CCOperators.h"

#include <cmath>
#include "../../madness/constants.h"
#include "../../madness/mra/derivative.h"
#include "../../madness/mra/funcdefaults.h"
#include "../../madness/mra/funcimpl.h"
#include "../../madness/mra/funcplot.h"
#include "../../madness/mra/function_factory.h"
#include "../../madness/mra/functypedefs.h"
#include "../../madness/mra/mra.h"
#include "../../madness/mra/operator.h"
#include "../../madness/mra/vmra.h"
#include "../../madness/tensor/srconf.h"
#include "../../madness/tensor/tensor.h"
#include "../../madness/world/madness_exception.h"
#include "../../madness/world/parallel_archive.h"
#include "../../madness/world/print.h"
#include "../../madness/world/world.h"
#include "electronic_correlation_factor.h"
#include "TDA.h"

namespace madness {






  real_function_6d
  CC_Operators::make_cc2_coulomb_parts(const CC_function &ti,const CC_function &tj,const CC_vecfunction &singles,const double omega) const {
    output("Make Screened Coulomb Potential with " + ti.name() + tj.name()+ " and projector is mixed with functions of type " + assign_name(singles.type));

    // failsafe
    if(ti.type==PARTICLE) warning("Particle function has entered cc2_coulomb parts, should be mixed or response or hole");
    if(tj.type==PARTICLE) warning("Particle function has entered cc2_coulomb parts, should be mixed or response or hole");
    if(ti.type!=RESPONSE and tj.type!=RESPONSE and omega!=0.0) error("make_cc2_coulomb_parts omega is not zero but non of the functions is of response type");

    real_convolution_6d G=BSHOperator<6>(world,sqrt(-2.0 * get_epsilon(ti.i,tj.i)+omega),parameters.lo,parameters.thresh_bsh_6D);
    G.destructive()=true;
    // first do the O1 and O2 parts which are
    // Otau1(g|titj) = |tauk><k|(1)g|titj> = kgti(2)|tauktj>
    // same for Otau2 = kgtj(1)|titauk>
    real_function_6d G_O1tau_part=real_factory_6d(world);
    real_function_6d G_O2tau_part=real_factory_6d(world);
    for(const auto& ktmp : singles.functions){
      const size_t k=ktmp.first;
      const CC_function &tauk=ktmp.second;

      real_function_3d kgti_tj=g12(mo_bra_(k),ti) * tj.function;
      real_function_3d kgtj_ti=g12(mo_bra_(k),tj) * ti.function;

      real_function_3d l1=real_factory_3d(world);
      real_function_3d l2=real_factory_3d(world);
      for(const auto& ltmp : singles.functions){
	const CC_function& mo_bra_l=mo_bra_(ltmp.first);
	const CC_function& taul=ltmp.second;
	l1+=mo_bra_l.function.inner(kgtj_ti) * taul.function;
	l2+=mo_bra_l.function.inner(kgti_tj) * taul.function;
      }

      real_function_3d part1=-1.0 * kgtj_ti + 0.5 * l1;
      real_function_3d part2=-1.0 * kgti_tj + 0.5 * l2;
      part1=projector_Q(part1);
      part2=projector_Q(part2);
      G_O1tau_part+=-2.0 * G(copy(tauk.function),part2);
      G_O2tau_part+=-2.0 * G(part1,copy(tauk.function));
    }

    return G_O1tau_part + G_O2tau_part;
  }

  real_function_6d CC_Operators::make_cc2_residue_sepparated(const CC_function &taui,const CC_function &tauj,const double omega) const {

    calctype ctype=CC2_;
    const bool symmetric=(taui.i == tauj.i);
    if(make_norm(taui) < parameters.thresh_3D and make_norm(tauj) < parameters.thresh_3D){
      output("Singles are zero: Current Calculation is MP2");
      ctype=MP2_;
    }
    CC_function ti;
    CC_function tj;
    if(ctype == CC2_){
      ti=make_t_intermediate(taui);
      tj=make_t_intermediate(tauj);
    }else if(ctype == MP2_){
      ti=mo_ket_(taui.i);
      tj=mo_ket_(tauj.i);
    }

    const double epsij=get_epsilon(taui.i,tauj.i);
    const double epsi=get_orbital_energies()[taui.i];
    const double epsj=get_orbital_energies()[tauj.i];
    if((epsi + epsj != epsij)) warning("Error in epsilon values: (epsi+epsj-epsij)=" + stringify(epsi + epsj - epsij));
    // Greens operator to apply later:
    real_convolution_6d G=BSHOperator<6>(world,sqrt(-2.0 * epsij),parameters.lo,parameters.thresh_bsh_6D);
    G.destructive()=true;
    // Greens operator to screen
    real_convolution_6d Gscreen=BSHOperator<6>(world,sqrt(-2.0 * epsij),parameters.lo,parameters.thresh_bsh_6D);
    Gscreen.modified()=true;

    real_function_3d F_ti=real_factory_3d(world);
    real_function_3d F_tj=real_factory_3d(world);
    if(ctype == CC2_){
      F_ti=(apply_F(ti) - epsi * ti.function).truncate();
      if(symmetric) F_tj=copy(F_ti);
      else F_tj=(apply_F(tj) - epsj * tj.function).truncate();
    }

    output_section("CC2-Residue-Unprojected-Part");
    CC_Timer time_unprojected(world,"CC2-Residue:Unprojected-Part");
    real_function_6d unprojected_result;
    real_function_6d unprojected_potential;
    {
      real_function_6d fFeij_part=real_factory_6d(world);
      if(ctype == CC2_){
	fFeij_part=make_f_xy_screened(F_ti,tj,Gscreen) + make_f_xy_screened(ti,F_tj,Gscreen);
      }
      const real_function_6d Uepot_part=apply_transformed_Ue(ti,tj);
      const real_function_6d Kf_fK_part=apply_exchange_commutator(ti,tj);

      const real_function_6d V=(fFeij_part + Uepot_part - Kf_fK_part).truncate().reduce_rank();
      unprojected_potential=copy(V);     // necessary because G is detructive
      Kf_fK_part.print_size("[K,f]" + ti.name() + tj.name() + "   ");
      Uepot_part.print_size("Ue" + ti.name() + tj.name() + "      ");
      fFeij_part.print_size("f(F-eij)" + ti.name() + tj.name() + "");
      V.print_size("-2.0(F-eij+Ue-[K,f])" + ti.name() + tj.name());
      const real_function_6d tmp=G(-2.0 * V);
      unprojected_result=tmp;
      unprojected_result.print_size("G(-2.0(F-eij+Ue-[K,f]))" + ti.name() + tj.name());
    }
    time_unprojected.info();

    output_section("CC2-Residue-Projected-Part");
    CC_Timer time_projected(world,"CC2-Residue:Projected-Part");
    const double tight_thresh=parameters.tight_thresh_6D;
    real_function_6d projected_result=real_factory_6d(world);
    projected_result.set_thresh(tight_thresh);
    output("Tighten thresh to " + stringify(tight_thresh));
    FunctionDefaults<6>::set_thresh(tight_thresh);
    {
      // the f(F-eij+K) operator is of type A12 = f12(A1+A2)
      // (O1+O1-O12)(A12) = k(1)*[(<k|A|x>(2)*y(2) - 1/2 <kl|A|xy> l(2)] + []*l(2)
      // 					= |k> (x) (kAxy_1 - 1/2 im_k) + (kAxy_2 - 1/2 im_k)(x)|k>
      // im_k = \sum_l <kl|A|xy> |l>
      //
      vecfuncT kAxy_1;
      vecfuncT kAxy_2;
      vecfuncT im_k1;
      vecfuncT im_k2;
      for(const auto & ktmp : mo_bra_.functions){
	const CC_function & k=ktmp.second;
	const real_function_3d kAxy1=unprojected_potential.project_out(k.function,0);
	const real_function_3d kAxy2=unprojected_potential.project_out(k.function,1);
	real_function_3d imk1=real_factory_3d(world);
	real_function_3d imk2=real_factory_3d(world);
	for(const auto & ltmp : mo_bra_.functions){
	  const CC_function & l=ltmp.second;
	  imk1+=l.inner(kAxy1) * mo_ket_(l).function;
	  imk2+=l.inner(kAxy2) * mo_ket_(l).function;
	}
	kAxy_1.push_back(kAxy1);
	kAxy_2.push_back(kAxy2);
	imk1.truncate();
	imk2.truncate();
	im_k1.push_back(imk1);
	im_k2.push_back(imk2);
      }

      for(const auto & ktmp : mo_ket_.functions){
	const CC_function & k=ktmp.second;
	const real_function_3d k1=copy(k.function);      // necessary because G is detructive
	const real_function_3d k2=copy(k.function);      // necessary because G is detructive
	const real_function_3d tmp1=kAxy_1[k.i] - 0.5 * im_k1[k.i];
	const real_function_6d part1=G(-2.0 * k1,tmp1);
	const real_function_3d tmp2=kAxy_2[k.i] - 0.5 * im_k2[k.i];
	const real_function_6d part2=G(tmp2,-2.0 * k2);
	projected_result+=(part1 + part2).truncate(tight_thresh);
      }
      projected_result.print_size("-2.0G[(O1+O2-O12)(fF-feij+Ue-[K,f])|" + ti.name() + tj.name() + ">]");
    }
    time_projected.info();
    output("Lowering thresh back to " + stringify(parameters.thresh_6D));
    FunctionDefaults<6>::set_thresh(parameters.thresh_6D);
    real_function_6d cc2_residue=unprojected_result - projected_result;
    cc2_residue.print_size("cc2_residue");
    apply_Q12(cc2_residue,"cc2_residue");
    cc2_residue.print_size("Q12cc2_residue");
    return cc2_residue;
  }

  double
  CC_Operators::compute_mp2_pair_energy(CC_Pair &pair) const {
    real_function_3d uninitialized_dummy;
    return compute_cc2_pair_energy(pair,uninitialized_dummy,uninitialized_dummy);
  }

  // The Fock operator is partitioned into F = T + Vn + R
  // the fock residue R= 2J-K+Un for closed shell is computed here
  // J_i = \sum_k <k|r12|k> |tau_i>
  // K_i = \sum_k <k|r12|tau_i> |k>
  vecfuncT
  CC_Operators::fock_residue_closed_shell(const CC_vecfunction &singles) const {
    //	vecfuncT tau = singles.get_vecfunction();
    //CC_Timer timer_J(world,"J");
    //	vecfuncT J = mul(world, intermediates_.get_hartree_potential(), tau);
    vecfuncT J;
    for(const auto& tmpi : singles.functions){
      const CC_function& taui=tmpi.second;
      real_function_3d hartree_potential= real_function_3d(world);
      for(const auto& tmpk :mo_ket_.functions) hartree_potential+=g12(mo_bra_(tmpk.first),tmpk.second);
      const real_function_3d Ji = hartree_potential*taui.function;
      J.push_back(Ji);
    }
    truncate(world,J);
    scale(world,J,2.0);
    //timer_J.info();
    //CC_Timer timer_K(world,"K");
    vecfuncT vK;
    for(const auto& tmpi : singles.functions){
      const CC_function& taui=tmpi.second;
      const real_function_3d Ki=K(taui);
      vK.push_back(Ki);
    }
    scale(world,vK,-1.0);
    //timer_K.info();

    // apply nuclear potential
    Nuclear Uop(world,&nemo);
    vecfuncT Upot=Uop(singles.get_vecfunction());
    vecfuncT KU=add(world,vK,Upot);

    return add(world,J,KU);
  }

  vecfuncT
  CC_Operators::ccs_potential(const CC_vecfunction &tau) const {
    // first form the intermediate t-functions: ti = i + taui
    const CC_vecfunction tfunctions=make_t_intermediate(tau);
    vecfuncT result;

    // get the perturbed hartree_potential: kgtk = sum_k <k|g|\tau_k>
    real_function_3d kgtauk=real_function_3d(world);
    for(const auto& ktmp:tau.functions) kgtauk += g12(mo_bra_(ktmp.first),ktmp.second);

    for(const auto& ttmp : tfunctions.functions){
      const CC_function& ti=ttmp.second;
      real_function_3d resulti=real_factory_3d(world);

      const real_function_3d kgtauk_ti=kgtauk * ti.function;
      real_function_3d kgti_tauk=real_factory_3d(world);
      for(const auto &ktmp : tau.functions){
	const CC_function& tauk=ktmp.second;
	const real_function_3d kgti=g12(mo_bra_(tauk.i),ti);
	kgti_tauk+=kgti * tauk.function;
      }

      real_function_3d l_kgtauk_ti_taul=real_function_3d(world);
      real_function_3d l_kgti_tauk_taul=real_function_3d(world);
      for(const auto &ltmp : tau.functions){
	const CC_function& taul=ltmp.second;
	l_kgtauk_ti_taul+=mo_bra_(taul).inner(kgtauk_ti) * taul.function;
	l_kgti_tauk_taul+=mo_bra_(taul).inner(kgti_tauk) * taul.function;
      }

      resulti=2.0 * kgtauk_ti - kgti_tauk - 2.0 * l_kgtauk_ti_taul + l_kgti_tauk_taul;
      result.push_back(resulti);
    }
    return result;
  }

  /// the ccs potential without terms from the fock operator
  /// returns: \f$ (1-|\tau_k><k|)(2 <k|g|tau_k> |t_i> - <k|g|t_i> |\tau_k>)  \f$
  vecfuncT
  CC_Operators::ccs_potential_new(const CC_vecfunction &tau) const {
    const CC_vecfunction ti=make_t_intermediate(tau);
    vecfuncT unprojected = ccs_unprojected(ti,tau);
    CC_Timer projected_time(world,"ccs-projected");
    vecfuncT projected = P(unprojected,tau);
    scale(world,projected,-1.0);
    projected_time.info();
    return add(world,unprojected,projected);
  }

  // returns 2kgtk|ti> - kgti|tk>
  // the ccs potential: ti = ti and tk = tauk
  vecfuncT CC_Operators::ccs_unprojected(const CC_vecfunction &ti, const CC_vecfunction &tk)const{
    vecfuncT result;
    for(const auto& itmp:ti.functions){
      real_function_3d kgtk = real_factory_3d(world);
      for(const auto& ktmp:tk.functions) kgtk+=g12(mo_bra_(ktmp.first),ktmp.second);
      const real_function_3d kgtk_ti = kgtk*ti(itmp.first).function;

      real_function_3d kgti_tk = real_factory_3d(world);
      for(const auto& ktmp:tk.functions) kgti_tk+=g12(mo_bra_(ktmp.first),ti(itmp.first))*tk(ktmp.first).function;

      const real_function_3d resulti= 2.0*kgtk_ti - kgti_tk;
      result.push_back(resulti);
    }
    return result;
  }
  vecfuncT CC_Operators::ccs_response_potential(const CC_vecfunction &singles, const CC_vecfunction &response)const{
    const CC_vecfunction ti=make_t_intermediate(singles);
    const vecfuncT tmp1 = ccs_unprojected(ti,response);
    const vecfuncT tmp2 = ccs_unprojected(response,singles);
    const vecfuncT part1 = add(world,tmp1,tmp2);
    const vecfuncT tmp3 = ccs_unprojected(ti,singles);
    const vecfuncT pt_part1 = P(part1,singles);
    const vecfuncT part2 = P(tmp3,response);
    const vecfuncT ccs_response_tmp = sub(world,part1,pt_part1);
    const vecfuncT ccs_response = sub(world,ccs_response_tmp,part2);
    return ccs_response;
  }

  vecfuncT
  CC_Operators::S2b_u_part(const Pairs<CC_Pair> &doubles,const CC_vecfunction &singles) const {
    vecfuncT result;
    if(singles.type==PARTICLE) result=copy(world,current_s2b_u_part_gs);
    else if(singles.type==RESPONSE) result=copy(world,current_s2b_u_part_response);
    else error("singles of type " + assign_name(singles.type) +" in S2b_u_part");

    if(not result.empty()){
      output("S2b_u_part already calculated");
    }else{
      for(const auto& itmp : singles.functions){
	const size_t i=itmp.first;
	real_function_3d resulti=real_factory_3d(world);
	for(const auto& ktmp : singles.functions){
	  const size_t k=ktmp.first;
	  const real_function_6d uik=get_pair_function(doubles,i,k);
	  // S2b u-part
	  resulti += 2.0* g12(mo_bra_(k),uik,2);
	  // S2b u-part-exchange
	  resulti -= g12(mo_bra_(k),uik,1);
	}
	result.push_back(resulti);
      }
      if(singles.type==PARTICLE) current_s2b_u_part_gs=copy(world,result);
      else if(singles.type==RESPONSE) current_s2b_u_part_response = copy(world,result);

    }
    return result;
  }

  vecfuncT
  CC_Operators::S2c_u_part(const Pairs<CC_Pair> &doubles,const CC_vecfunction &singles) const {
    vecfuncT result;
    if(singles.type==PARTICLE) result=copy(world,current_s2c_u_part_gs);
    else if(singles.type==RESPONSE) result=copy(world,current_s2c_u_part_response);
    else error("singles of type " + assign_name(singles.type) +" in S2c_u_part");

    if(not result.empty()){
      output("S2c_u_part already calculated");
    }else{
      for(const auto& itmp : singles.functions){
	const size_t i=itmp.first;
	real_function_3d resulti=real_factory_3d(world);
	for(const auto& ktmp : singles.functions){
	  const size_t k=ktmp.first;
	  const real_function_3d kgi=g12(mo_bra_(k),mo_ket_(i));
	  for(const auto& ltmp : singles.functions){
	    const size_t l=ltmp.first;
	    const real_function_6d ukl=get_pair_function(doubles,k,l);
	    const real_function_3d l_kgi=mo_bra_(l).function * kgi;
	    resulti+=-2.0 * ukl.project_out(l_kgi,1);     // 1 means second particle
	    resulti+=ukl.project_out(l_kgi,0);
	  }
	}
	result.push_back(resulti);
      }
      if(singles.type==PARTICLE) current_s2c_u_part_gs=copy(world,result);
      else if(singles.type==RESPONSE) current_s2c_u_part_response = copy(world,result);
    }
    return result;
  }

  /// The Part of the CC2 singles potential which depends on singles and doubles (S4a, S4b, S4c)
  vecfuncT
  CC_Operators::S4a_u_part(const Pairs<CC_Pair> &doubles,const CC_vecfunction &singles) const {
    // S4a can be computed from the S2b potential
    // (-2<lk|g|uik> + <kl|g|uik>)|tau_l> =( <l( (-2)*<k|g|uik>_2) + <l| (<k|g|uik>_1) )|tau_l> = <l|s2b_u_part>*|tau_l> = - \sum_l <l|s2b_i> |l> important: minus sign and the fact that the s2b potential needs to be unprojected
    vecfuncT s4a;
    vecfuncT s2b;
    if(singles.type==PARTICLE) s2b=copy(world,current_s2b_u_part_gs);
    else if(singles.type==RESPONSE) s2b=copy(world,current_s2b_u_part_response);
    else warning("S4a_u_part: Singles are of type "+assign_name(singles.type));
    if(s2b.empty()){
    for(const auto& itmp : singles.functions){
      const size_t i=itmp.first;
      real_function_3d s4ai=real_factory_3d(world);
      real_function_3d s4ai_consistency=real_factory_3d(world);     // here the unprojected s2b result will be used to check consistency since this is not expensive this will be used everytime the s2b part was stored
      for(const auto& ltmp : singles.functions){
	const size_t l=ltmp.first;
	const CC_function& taul=ltmp.second;
	for(const auto& ktmp : singles.functions){
	  const size_t k=ktmp.first;
	  s4ai+=(-2.0 * make_ijgu(l,k,get_pair_function(doubles,i,k)) + make_ijgu(k,l,get_pair_function(doubles,i,k))) * taul.function;
	}
      }
      s4a.push_back(s4ai);
    }
    }else{
      output("S4a from stored S2b potential");
      const vecfuncT active_mo_bra = get_active_mo_bra();
      for(const auto& itmp : singles.functions){
	const Tensor<double> ls2b = inner(world,s2b[itmp.first-parameters.freeze],active_mo_bra);
	real_function_3d resulti = real_factory_3d(world);
	for(const auto& ltmp : singles.functions){
	  resulti -= ls2b(ltmp.first-parameters.freeze)*ltmp.second.function;
	}
	s4a.push_back(resulti);
      }
    }
    return s4a;
  }

  // result: -\sum_k( <l|kgtaui|ukl>_2 - <l|kgtaui|ukl>_1) | kgtaui = <k|g|taui>
  vecfuncT
  CC_Operators::S4b_u_part(const Pairs<CC_Pair> &doubles,const CC_vecfunction &singles) const {
    vecfuncT result;
    const vecfuncT active_mo_bra = get_active_mo_bra();
    for(const auto& itmp : singles.functions){
      const size_t i=itmp.first;
      real_function_3d resulti=real_factory_3d(world);
      for(const auto& ktmp : singles.functions){
	const size_t k=ktmp.first;
	const real_function_3d kgi = g12(mo_bra_(k),singles(i));
	vecfuncT l_kgi = mul_sparse(world,kgi,active_mo_bra,parameters.thresh_3D);
	truncate(world,l_kgi);
	for(const auto& ltmp : singles.functions){
	  const size_t l=ltmp.first;
	  const real_function_6d ukl=get_pair_function(doubles,k,l);
	  //const real_function_3d l_kgi=mo_bra_(l).function * kgi;
	  resulti+=-2.0 * ukl.project_out(l_kgi[l-parameters.freeze],1);     // 1 means second particle
	  resulti+=ukl.project_out(l_kgi[l-parameters.freeze],0);
	}
      }
      resulti.print_size("s4b_"+std::to_string(i));
      result.push_back(resulti);
    }
    return result;
  }

  vecfuncT
  CC_Operators::S4c_u_part(const Pairs<CC_Pair> &doubles,const CC_vecfunction &singles) const {
    vecfuncT result;
    const vecfuncT active_mo_bra = get_active_mo_bra();
    for(const auto& itmp : singles.functions){
      const size_t i=itmp.first;
      real_function_3d resulti=real_factory_3d(world);
      real_function_3d part1=real_factory_3d(world);
      real_function_3d part2=real_factory_3d(world);
      real_function_3d part3=real_factory_3d(world);
      real_function_3d part4=real_factory_3d(world);

      real_function_3d kgtauk= real_factory_3d(world);
      for(const auto& ktmp : singles.functions){
	const size_t k=ktmp.first;
	kgtauk += g12(mo_bra_(k),singles(k));
      }
      vecfuncT l_kgtauk = mul(world,kgtauk,active_mo_bra);
      truncate(world,l_kgtauk);

      for(const auto& ltmp : singles.functions){
	const size_t l=ltmp.first;
	const real_function_6d uil=get_pair_function(doubles,i,l);
	part1+=uil.project_out(l_kgtauk[l-parameters.freeze],1);
	part2+=uil.project_out(l_kgtauk[l-parameters.freeze],0);

	for(const auto& ktmp : singles.functions){
	  const size_t k=ktmp.first;
	  const real_function_3d k_lgtauk=(mo_bra_(k).function * g12(mo_bra_(l),singles(k))).truncate();
	  part3+=uil.project_out(k_lgtauk,1);
	  part4+=uil.project_out(k_lgtauk,0);
	}
      }
      resulti=4.0 * part1 - 2.0 * part2 - 2.0 * part3 + part4;
      result.push_back(resulti);
    }
    return result;
  }

  // the s2b potential with u=Qf|reg1,reg2>
  vecfuncT
  CC_Operators::S2b_reg_part(const CC_vecfunction &reg1, const CC_vecfunction &reg2) const {
    vecfuncT result;
    real_function_3d ktk=real_factory_3d(world);
    for(const auto& ktmp:reg2.functions) ktk+=mo_bra_(ktmp.first).function*ktmp.second.function;

    const real_function_3d kgftk=apply_gf(ktk);
    for(const auto& itmp : reg1.functions){			// convenience
      const CC_function& ti=itmp.second;				// convenience
      real_function_3d resulti=real_factory_3d(world);				// this will be the result
      real_function_3d Ipart=real_factory_3d(world);
      real_function_3d Ipartx=real_factory_3d(world);
      real_function_3d O1part=real_factory_3d(world);
      real_function_3d O1partx=real_factory_3d(world);
      real_function_3d O2part=real_factory_3d(world);
      real_function_3d O2partx=real_factory_3d(world);
      real_function_3d O12part=real_factory_3d(world);
      real_function_3d O12partx=real_factory_3d(world);
      Ipart+=2.0 * kgftk * ti.function;     // part1
      for(const auto& ktmp : reg2.functions){
	const size_t k=ktmp.first;
	const CC_function& tk=ktmp.second;
	const real_function_3d kti=mo_bra_(k).function * ti.function;
	const real_function_3d kgfti=apply_gf(kti);
	Ipartx+=-1.0 * kgfti * tk.function;     // part1x

	for(const auto& mtmp : mo_ket_.functions){
	  const size_t m=mtmp.first;
	  const CC_function& mom=mtmp.second;
	  const real_function_3d mftk= f12(mo_bra_(m),tk);
	  const real_function_3d mfti= f12(mo_bra_(m),ti);
	  const real_function_3d kgm= g12(mo_bra_(k),mo_ket_(m));
	  const real_function_3d mfti_tk=mfti * tk.function;
	  const real_function_3d mftk_ti=mftk * ti.function;
	  O2part-=(2.0 * kgm * mftk_ti);     //part3
	  O2partx-=(-1.0 * kgm * mfti_tk);
	  const real_function_3d k_mfti_tk=mo_bra_(k).function * mfti_tk;
	  const real_function_3d k_gmfti_tk=g12(k_mfti_tk);
	  const real_function_3d k_mftk_ti=mo_bra_(k).function * mftk_ti;
	  const real_function_3d k_gmftk_ti=g12(k_mftk_ti);
	  O1part-=(2.0 * k_gmfti_tk * mom.function);     //part2
	  O1partx-=(-1.0 * k_gmftk_ti * mom.function);
	  for(const auto& ntmp : mo_ket_.functions){
	    const CC_function& mon=ntmp.second;
	    const double nmftitk=mo_bra_(mon).inner(mftk_ti);
	    const double nmftkti=mo_bra_(mon).inner(mfti_tk);
	    O12part+=(2.0 * nmftitk * kgm * mon.function);
	    O12partx+=(-1.0 * nmftkti * kgm * mon.function);
	  }     // end n
	}     // end m
      }     // end k
      resulti=Ipart + Ipartx + O1part + O1partx + O2part + O2partx + O12part + O12partx;
      result.push_back(resulti);
    }     // end i
    return result;
  }

  vecfuncT
  CC_Operators::S2c_reg_part(const CC_vecfunction &reg1, const CC_vecfunction &reg2) const {
    vecfuncT result;
    for(const auto& itmp : reg1.functions){
      const size_t i=itmp.first;
      real_function_3d resulti=real_factory_3d(world);
      for(const auto& ktmp : reg1.functions){
	const CC_function tk=ktmp.second;
	const size_t k=ktmp.first;
	for(const auto& ltmp : reg2.functions){
	  const CC_function tl=ltmp.second;
	  const real_function_3d l_kgi_tmp=mo_bra_(tl).function * g12(mo_bra_(k),mo_ket_(i));
	  const CC_function l_kgi(l_kgi_tmp,99,UNDEFINED);
	  resulti-=(2.0 * convolute_x_Qf_yz(l_kgi,tk,tl) - convolute_x_Qf_yz(l_kgi,tl,tk));
	}
      }
      result.push_back(resulti);
    }
    return result;
  }

  vecfuncT
  CC_Operators::S4a_reg_part(const CC_vecfunction &reg1, const CC_vecfunction &reg2, const CC_vecfunction &singles) const {
    vecfuncT result;
    for(const auto& itmp : reg1.functions){
      const CC_function& ti=itmp.second;
      real_function_3d resulti=real_factory_3d(world);

      for(const auto& ktmp : reg2.functions){
	const CC_function& tk=ktmp.second;
	const size_t k=ktmp.first;

	for(const auto& ltmp : singles.functions){
	  const CC_function& taul=ltmp.second;
	  const size_t l=ltmp.first;

	  const double lkgQftitk=make_ijgQfxy(l,k,ti,tk);
	  const double klgQftitk=make_ijgQfxy(k,l,ti,tk);
	  const double factor =(2.0 * lkgQftitk - klgQftitk);
	  resulti-= factor* taul.function;
	}
      }
      result.push_back(resulti);
    }
    return result;
  }

  /// result: -\sum_{kl}( 2 <l|kgtaui|Qftktl> - <l|kgtaui|Qftltk>
  /// this is the same as S2c with taui instead of i
  vecfuncT
  CC_Operators::S4b_reg_part(const CC_vecfunction &reg1, const CC_vecfunction &reg2,const CC_vecfunction &singles) const {
    vecfuncT result;
    for(const auto& itmp : singles.functions){
      const CC_function taui=itmp.second;
      real_function_3d resulti=real_factory_3d(world);

      for(const auto& ktmp : reg1.functions){
	const CC_function tk=ktmp.second;
	for(const auto& ltmp : reg2.functions){
	  const CC_function tl=ltmp.second;
	  const real_function_3d l_kgi_tmp=msparse(mo_bra_(tl).function,g12(mo_bra_(tk.i),taui));
	  const CC_function l_kgi(l_kgi_tmp,99,UNDEFINED);
	  resulti-=(2.0 * convolute_x_Qf_yz(l_kgi,tk,tl) - convolute_x_Qf_yz(l_kgi,tl,tk));
	}
      }
      result.push_back(resulti);
    }
    return result;
  }

  /// result: 4<l|kgtauk|Qftitl> - 2<l|kgtauk|Qftlti> - 2<k|lgtauk|Qftitl> + <k|lgtauk|Qftlti>
  vecfuncT
  CC_Operators::S4c_reg_part(const CC_vecfunction &reg1, const CC_vecfunction &reg2,const CC_vecfunction &singles) const {
    vecfuncT result;
    for(const auto& itmp : reg1.functions){
      const CC_function& ti=itmp.second;
      real_function_3d resulti=real_factory_3d(world);

      real_function_3d kgtauk=real_factory_3d(world);
      for(const auto& ktmp : singles.functions){
	const CC_function& tauk=ktmp.second;
	kgtauk += g12(mo_bra_(tauk.i),tauk);
      }

      // first two parts
      real_function_3d part1=real_factory_3d(world);
      real_function_3d part2=real_factory_3d(world);
      for(const auto& ltmp : reg2.functions){
	const CC_function& tl=ltmp.second;
	const size_t l=ltmp.first;
	const real_function_3d l_kgtauk=msparse(mo_bra_(l).function,kgtauk);
	part1+=convolute_x_Qf_yz(CC_function(l_kgtauk,99,UNDEFINED),ti,tl);
	part2+=convolute_x_Qf_yz(CC_function(l_kgtauk,99,UNDEFINED),tl,ti);
      }

      // second two parts
      real_function_3d part3=real_factory_3d(world);
      real_function_3d part4=real_factory_3d(world);
      for(const auto& ktmp : singles.functions){
	const CC_function& tauk=ktmp.second;
	const size_t k=ktmp.first;

	for(const auto& ltmp : reg2.functions){
	  const CC_function& tl=ltmp.second;
	  const size_t l=ltmp.first;

	  const real_function_3d k_lgtauk=msparse(mo_bra_(k).function,g12(mo_bra_(l),tauk));
	  part3+=convolute_x_Qf_yz(CC_function(k_lgtauk,99,UNDEFINED),ti,tl);
	  part4+=convolute_x_Qf_yz(CC_function(k_lgtauk,99,UNDEFINED),tl,ti);
	}
      }
      resulti=4.0 * part1 - 2.0 * part2 - 2.0 * part3 + part4;
      result.push_back(resulti);
    }
    return result;
  }

  real_function_6d CC_Operators::apply_regularization_potential(const CC_function &a, const CC_function &b, const double omega)const{
    if((a.type != RESPONSE and b.type != RESPONSE) and omega !=0.0) error("Apply_regularization_potential: omega is not zero, but none of the functions has response type");
    const real_function_3d Fa=apply_F(a) - (get_orbital_energies()[a.i]+0.5*omega) * a.function;
    const real_function_3d Fb=apply_F(b) - (get_orbital_energies()[b.i]+0.5*omega) * b.function;

    // make the Fock operator part:  f(F-eij)|titj> = (F1+F2-ei-ej)|titj> = (F1-ei)|ti>|tj> + |ti>(F2-ei)|tj>
    const real_function_6d fFab=make_f_xy(Fa,b,guess_thresh(Fa,b)) + make_f_xy(a,Fb,guess_thresh(a,Fb));

    output("Applying Regularization Potential");

    // make the (U-[K,f])|titj> part
    // first the U Part
    const real_function_6d Uab=apply_transformed_Ue(a,b);
    // then the [K,f] part
    const real_function_6d KffKab=apply_exchange_commutator(a,b);

    real_function_6d V=(fFab + Uab - KffKab);
    V.print_size("Vreg"+a.name()+b.name());
    fFab.print_size("Vreg:f12(F-eij)|" + a.name() + b.name() + ">");
    Uab.print_size("          Vreg:U|" + a.name() + b.name() + ">");
    KffKab.print_size("      Vreg:[K,f]|" + a.name() + b.name() + ">");
    return V;
  }


  /// Make the CC2 Residue which is:  Q12f12(T-eij + 2J -K +Un )|titj> + Q12Ue|titj> - [K,f]|titj>  with |ti> = |\taui>+|i>
  /// @param[in] \tau_i which will create the |t_i> = |\tau_i>+|i> intermediate
  /// @param[in] \tau_j
  /// @param[in] u, the uij pair structure which holds the consant part of MP2
  /// @param[out] Q12f12(F-eij)|titj> + Q12Ue|titj> - [K,f]|titj>  with |ti> = |\taui>+|i>
  /// Right now Calculated in the decomposed form: |titj> = |i,j> + |\taui,\tauj> + |i,\tauj> + |\taui,j>
  /// The G_Q_Ue and G_Q_KffK part which act on |ij> are already calculated and stored as constant_term in u (same as for MP2 calculations) -> this should be the biggerst (faster than |titj> form)
  real_function_6d
  CC_Operators::make_regularization_residue(const CC_function &a,const CC_function &b, const calctype &ctype, const double omega) const {
    output("Calculating GV_reg|"+a.name()+b.name()+">");
    consistency_check(a,b,omega);
    const bool symmetric(a.i==b.i);

    real_function_6d Vreg;
    if(ctype==CC2_){
      Vreg = apply_regularization_potential(make_t_intermediate(a),make_t_intermediate(b),0.0);
    }else if(ctype==MP2_){
      Vreg = apply_regularization_potential(mo_ket_(a.i),mo_ket_(b.i),0.0);
    }else if(ctype==CISpD_){
      Vreg = apply_regularization_potential(a,mo_ket_(b.i),omega);
      if(symmetric){
	output("Exploiting Symmetry for Diagonal pairing: " +a.name()+mo_ket_(b.i).name() +" = P12" + mo_ket_(a.i).name()+b.name());
	Vreg = Vreg + swap_particles(Vreg);
      }
      else Vreg = Vreg = Vreg + apply_regularization_potential(mo_ket_(a.i),b,omega);
    }else if(ctype==CC2_response_){
      error("CC2_response residue not yet implemented");
    }else error("Error in make_cc2_residue: Unknown ctype:"+assign_name(ctype));

    Vreg.scale(-2.0);
    apply_Q12(Vreg,"Vreg");
    Vreg.truncate().reduce_rank();
    Vreg.print_size("-2.0*Q12Vreg");

    real_convolution_6d G=BSHOperator<6>(world,sqrt(-2.0 * get_epsilon(a.i,b.i)+omega),parameters.lo,parameters.thresh_bsh_6D);
    G.destructive()=true;
    real_function_6d GV=G(Vreg);
    apply_Q12(GV,"CC2-Residue:G(V)");
    return GV;
  }



  /// The 6D Fock residue on the cusp free pair function u_{ij}(1,2) is: (2J - Kn - Un)|u_{ij}>
  real_function_6d
  CC_Operators::fock_residue_6d(const CC_Pair &u) const {
    const double eps=get_epsilon(u.i,u.j);
    // make the coulomb and local Un part with the composite factory
    real_function_3d hartree_potential = real_factory_3d(world);
    for(const auto& tmp:mo_ket_.functions) hartree_potential += g12(mo_bra_(tmp.first),mo_ket_(tmp.first));
    real_function_3d local_part=(2.0 * hartree_potential + nemo.nuclear_correlation->U2());
    local_part.print_size("vlocal");
    u.function.print_size(u.name());

    // Contruct the BSH operator in order to screen

    real_convolution_6d op_mod=BSHOperator<6>(world,sqrt(-2.0 * eps),parameters.lo,parameters.thresh_bsh_6D);
    op_mod.modified()=true;
    // Make the CompositeFactory
    real_function_6d vphi=CompositeFactory<double, 6, 3>(world).ket(copy(u.function)).V_for_particle1(copy(local_part)).V_for_particle2(copy(local_part));
    // Screening procedure
    vphi.fill_tree(op_mod);

    vphi.print_size("vlocal|u>");

    // the part with the derivative operators: U1
    for(int axis=0; axis < 6; ++axis){
      real_derivative_6d D=free_space_derivative<double, 6>(world,axis);
      // Partial derivative of the pari function
      const real_function_6d Du=D(u.function).truncate();

      // % operator gives division rest (modulo operator)
      if(world.rank() == 0) print("axis, axis^%3, axis/3+1",axis,axis % 3,axis / 3 + 1);
      const real_function_3d U1_axis=nemo.nuclear_correlation->U1(axis % 3);

      double tight_thresh=parameters.tight_thresh_6D;
      if(tight_thresh > 1.e-4) warning("tight_thresh_6D is too low for Un potential");
      real_function_6d x;
      if(axis / 3 + 1 == 1){
	x=CompositeFactory<double, 6, 3>(world).ket(copy(Du)).V_for_particle1(copy(U1_axis)).thresh(tight_thresh);

      }else if(axis / 3 + 1 == 2){
	x=CompositeFactory<double, 6, 3>(world).ket(copy(Du)).V_for_particle2(copy(U1_axis)).thresh(tight_thresh);
      }
      x.fill_tree(op_mod);
      x.set_thresh(FunctionDefaults<6>::get_thresh());
      x.print_size("Un_axis_" + stringify(axis));
      vphi+=x;
      vphi.truncate().reduce_rank();
    }

    vphi.print_size("(Un + J1 + J2)|u>");

    // Exchange Part
    vphi=(vphi - K(u.function,u.i == u.j)).truncate().reduce_rank();
    vphi.print_size("(Un + J - K)|u>");
    return vphi;

  }

  /// Echange Operator on 3D function
  /// !!!!Prefactor (-1) is not included
  real_function_3d
  CC_Operators::K(const CC_function &x) const {
    return apply_K(x);
  }
  real_function_3d
  CC_Operators::K(const real_function_3d &x) const {
    const CC_function tmp(x,99,UNDEFINED);
    return apply_K(tmp);
  }

  /// Exchange Operator on Pair function: -(K(1)+K(2))u(1,2)
  /// if i==j in uij then the symmetry will be exploited
  /// !!!!Prefactor (-1) is not included here!!!!
  real_function_6d
  CC_Operators::K(const real_function_6d &u,const bool symmetric,const double thresh) const {
    real_function_6d result=real_factory_6d(world).compressed();
    // K(1) Part
    result+=apply_K(u,1,thresh);
    // K(2) Part
    if(symmetric) result+=swap_particles(result);
    else result+=apply_K(u,2,thresh);

    result.print_size("K(fxy)_untruncated");
    return (result.truncate(parameters.tight_thresh_6D));
  }

  /// Exchange Operator on Pair function: -(K(1)+K(2))u(1,2)
  /// K(1)u(1,2) = \sum_k <k(3)|g13|u(3,2)> |k(1)>
  /// 1. X(3,2) = bra_k(3)*u(3,2)
  /// 2. Y(1,2) = \int X(3,2) g13 d3
  /// 3. result = Y(1,2)*ket_k(1)
  /// !!!!Prefactor (-1) is not included here!!!!
  real_function_6d
  CC_Operators::apply_K(const real_function_6d &u,const size_t &particle,const double thresh) const {
    MADNESS_ASSERT(particle == 1 or particle == 2);
    //poisson->particle()=particle;
    real_function_6d result=real_factory_6d(world).compressed();
    for(size_t k=0; k < mo_ket_.size(); k++){
      real_function_6d copyu=copy(u);
      copyu.set_thresh(thresh);
      real_function_6d X=(multiply(copyu,copy(mo_bra_(k).function),particle)).truncate(thresh);
      X.set_thresh(thresh);
      //      real_function_6d Y=(*poisson)(X);
      real_function_6d Y=g12(X,particle);
      Y.set_thresh(thresh);
      result+=(multiply(copy(Y),copy(mo_ket_(k).function),particle)).truncate(thresh);
    }
    return result;
  }

  // the K operator runs over ALL orbitals (also the frozen ones)
  real_function_3d
  CC_Operators::apply_K(const CC_function &f) const {
    //    if(parameters.debug and world.rank() == 0) std::cout << "apply K on " << assign_name(f.type) << " function" << std::endl;
    //    if(parameters.debug and world.rank() == 0) std::cout << "K" << f.name() << "=";
    //    real_function_3d result=real_factory_3d(world);
    //    switch(f.type){
    //      case HOLE:
    //	for(auto k_iterator : mo_ket_.functions){
    //	  const CC_function& k=k_iterator.second;
    //	  const real_function_3d tmp=intermediates_.get_EX(k,f);
    //	  result+=(tmp * k.function);
    //	  if(parameters.debug and world.rank() == 0) std::cout << "+ <" << k.name() << "|g|" << f.name() << ">*" << k.name();
    //	}
    //	break;
    //      case PARTICLE:
    //	for(auto k_iterator : mo_ket_.functions){
    //	  const CC_function& k=k_iterator.second;
    //	  result+=(intermediates_.get_pEX(k,f) * k.function);
    //	  if(parameters.debug and world.rank() == 0) std::cout << "+ <" << k.name() << "|g|" << f.name() << ">*" << k.name();
    //	}
    //	break;
    //      case MIXED:
    //	for(auto k_iterator : mo_ket_.functions){
    //	  const CC_function& k=k_iterator.second;
    //	  result+=(intermediates_.get_EX(k,f) + intermediates_.get_pEX(k,f)) * k.function;
    //	  if(parameters.debug and world.rank() == 0) std::cout << "+ <" << k.name() << "|g|t" << f.i << ">*" << k.name();
    //	}
    //	break;
    //      default:
    //	for(auto k_iterator : mo_ket_.functions){
    //	  const CC_function& k=k_iterator.second;
    //	  real_function_3d tmp=((*poisson)(mo_bra_(k).function * f.function)).truncate();
    //	  result+=tmp * k.function;
    //	  if(parameters.debug and world.rank() == 0) std::cout << "+ poisson(mo_bra_" << k.i << "*" << f.name() << ")|mo_ket_" << k.i << ">" << std::endl;
    //	}
    //	break;
    //    }
    //    return result;
    real_function_3d result=real_factory_3d(world);
    for(auto k_iterator : mo_ket_.functions){
      result+=g12(mo_bra_(k_iterator.first),f)*mo_ket_(k_iterator.first).function;
    }
    return result;
  }

  /// Apply Ue on a tensor product of two 3d functions: Ue(1,2) |x(1)y(2)> (will be either |ij> or |\tau_i\tau_j> or mixed forms)
  /// The Transformed electronic regularization potential (Kutzelnigg) is R_{12}^{-1} U_e R_{12} with R_{12} = R_1*R_2
  /// It is represented as: R_{12}^{-1} U_e R_{12} = U_e + R^-1[Ue,R]
  /// where R^-1[Ue,R] = R^-1 [[T,f],R] (see: Regularizing the molecular potential in electronic structure calculations. II. Many-body
  /// methods, F.A.Bischoff)
  /// The double commutator can be evaluated as follows:  R^-1[[T,f],R] = -Ue_{local}(1,2)*(Un_{local}(1) - Un_{local}(2))
  /// @param[in] x the 3D function for particle 1
  /// @param[in] y the 3D function for particle 2
  /// @param[in] i the first index of the current pair function (needed to construct the BSH operator for screening)
  /// @param[in] j the second index of the current pair function
  /// @param[out]  R^-1U_eR|x,y> the transformed electronic smoothing potential applied on |x,y> :
  real_function_6d
  CC_Operators::apply_transformed_Ue(const CC_function &x,const CC_function &y) const {
    // make shure the thresh is high enough
    CC_Timer time_Ue(world,"Ue|" + x.name() + y.name() + ">");
    const size_t i=x.i;
    const size_t j=y.i;
    double tight_thresh=parameters.tight_thresh_6D;
    output("Applying transformed Ue with 6D thresh = " + stringify(tight_thresh));

    real_function_6d Uxy=real_factory_6d(world);
    Uxy.set_thresh(tight_thresh);
    // Apply the untransformed U Potential
    const double eps=get_epsilon(i,j);
    Uxy=corrfac.apply_U(x.function,y.function,eps);
    Uxy.set_thresh(tight_thresh);

    // Get the 6D BSH operator in modified-NS form for screening
    real_convolution_6d op_mod=BSHOperator<6>(world,sqrt(-2.0 * eps),parameters.lo,parameters.thresh_bsh_6D);
    op_mod.modified()=true;

    // Apply the double commutator R^{-1}[[T,f,R]
    for(size_t axis=0; axis < 3; axis++){
      // Make the local parts of the Nuclear and electronic U potentials
      const real_function_3d Un_local=nemo.nuclear_correlation->U1(axis);
      const real_function_3d Un_local_x=(Un_local * x.function).truncate();
      const real_function_3d Un_local_y=(Un_local * y.function).truncate();
      const real_function_6d Ue_local=corrfac.U1(axis);
      // Now add the Un_local_x part to the first particle of the Ue_local potential
      real_function_6d UeUnx=CompositeFactory<double, 6, 3>(world).g12(Ue_local).particle1(Un_local_x).particle2(copy(y.function)).thresh(tight_thresh);
      // Fill the Tree were it will be necessary
      UeUnx.fill_tree(op_mod);
      // Set back the thresh
      UeUnx.set_thresh(FunctionDefaults<6>::get_thresh());

      //UeUnx.print_size("UeUnx");

      // Now add the Un_local_y part to the second particle of the Ue_local potential
      real_function_6d UeUny=CompositeFactory<double, 6, 3>(world).g12(Ue_local).particle1(copy(x.function)).particle2(Un_local_y).thresh(tight_thresh);
      // Fill the Tree were it will be necessary
      UeUny.fill_tree(op_mod);
      // Set back the thresh
      UeUny.set_thresh(FunctionDefaults<6>::get_thresh());

      //UeUny.print_size("UeUny");

      // Construct the double commutator part and add it to the Ue part
      real_function_6d diff=(UeUnx - UeUny).scale(-1.0);
      diff.truncate();
      Uxy=(Uxy + diff).truncate();
    }
    time_Ue.info();

    // sanity check: <xy|R2 [T,g12] |xy> = <xy |R2 U |xy> - <xy|R2 g12 | xy> = 0
    CC_Timer time_sane(world,"Ue-Sanity-Check");
    real_function_6d tmp=CompositeFactory<double, 6, 3>(world).particle1(copy(x.function * nemo.nuclear_correlation->square())).particle2(copy(y.function * nemo.nuclear_correlation->square()));
    const double a=inner(Uxy,tmp);
    const real_function_3d xx=(x.function * x.function * nemo.nuclear_correlation->square());
    const real_function_3d yy=(y.function * y.function * nemo.nuclear_correlation->square());
    const real_function_3d gxx=g12(xx);
    const double aa=inner(yy,gxx);
    const double error=std::fabs(a - aa);
    time_sane.info();
    if(world.rank() == 0 and error > FunctionDefaults<6>::get_thresh()){
      printf("<xy| U_R |xy>  %12.8f\n",a);
      printf("<xy|1/r12|xy>  %12.8f\n",aa);
      warning("Ue Potential Inaccurate!");
      if(error > FunctionDefaults<6>::get_thresh() * 10.0) warning("Ue Potential wrong !!!!");
    }else output("Ue seems to be sane");
    return Uxy;
  }

  /// Apply the Exchange Commutator [K,f]|xy>
  real_function_6d
  CC_Operators::apply_exchange_commutator(const CC_function &x,const CC_function &y,const double thresh) const {
    output("Computing [K,f]|" + x.name() + y.name() + "> with thresh=" + stringify(thresh));
    CC_Timer time(world,"[K,f]|" + x.name() + y.name() + ">");
    // make first part of commutator
    CC_Timer part1_time(world,"Kf" + x.name() + y.name() + ">");
    real_function_6d Kfxy=apply_Kf(x,y,thresh);
    part1_time.info();
    // make the second part of the commutator
    CC_Timer part2_time(world,"fK" + x.name() + y.name() + ">");
    real_function_6d fKxy=apply_fK(x,y,thresh);
    part2_time.info();
    Kfxy.print_size("Kf" + x.name() + y.name());
    fKxy.print_size("fK" + x.name() + y.name());
    real_function_6d result=(Kfxy - fKxy);

    time.info();
    // sanity check
    // <psi|[A,B]|psi> = <psi|AB|psi> - <psi|BA|psi> = <Apsi|Bpsi> - <Bpsi|Apsi> = 0 (if A,B hermitian)
    {
      CC_Timer sanity(world,"[K,f] sanity check");
      // make the <xy| bra state which is <xy|R2
      const real_function_3d brax=x.function * nemo.nuclear_correlation->square();
      const real_function_3d bray=y.function * nemo.nuclear_correlation->square();
      const real_function_6d xy=make_xy(CC_function(brax,x.i,x.type),CC_function(bray,y.i,y.type));
      const double xyfKxy=xy.inner(fKxy);
      const double xyKfxy=xy.inner(Kfxy);
      const double diff=xyfKxy - xyKfxy;
      if(world.rank() == 0 and fabs(diff) > FunctionDefaults<6>::get_thresh()){
	std::cout << std::setprecision(parameters.output_prec);
	std::cout << "<" << x.name() << y.name() << "|fK|" << x.name() << y.name() << "> =" << xyfKxy << std::endl;
	std::cout << "<" << x.name() << y.name() << "|Kf|" << x.name() << y.name() << "> =" << xyKfxy << std::endl;
	std::cout << "difference = " << diff << std::endl;
	warning("Exchange Commutator Plain Wrong");
      }

      else if(fabs(diff) > FunctionDefaults<6>::get_thresh() * 0.1) warning("Exchange Commutator critical");
      else output("Exchange Commutator seems to be sane");

      sanity.info();

      if(fabs(diff) > parameters.thresh_6D){
	if(thresh == parameters.tight_thresh_6D) output("Will not raise thresh anymore, already at tight_thresh");
	else{
	  output("Recalculating [K,f]" + x.name() + y.name() + " with thresh=" + stringify(parameters.tight_thresh_6D));
	  result=apply_exchange_commutator(x,y,parameters.tight_thresh_6D);
	}
      }
    }
    return result;
  }

  /// Apply the Exchange operator on a tensor product multiplied with f12
  /// !!! Prefactor of (-1) is not inclued in K here !!!!
  real_function_6d
  CC_Operators::apply_Kf(const CC_function &x,const CC_function &y,const double thresh) const {
    bool symmetric=false;
    if((x.type == y.type) and (x.i == y.i)) symmetric=true;

    // First make the 6D function f12|x,y>
    real_function_6d f12xy=make_f_xy(x,y,thresh);
    // Apply the Exchange Operator
    real_function_6d result=K(f12xy,symmetric,thresh);
    return result;
  }

  /// Apply fK on a tensor product of two 3D functions
  /// fK|xy> = fK_1|xy> + fK_2|xy>
  /// @param[in] x, the first 3D function in |xy>, structure holds index i and type (HOLE, PARTICLE, MIXED, UNDEFINED)
  /// @param[in] y, the second 3D function in |xy>  structure holds index i and type (HOLE, PARTICLE, MIXED, UNDEFINED)
  real_function_6d
  CC_Operators::apply_fK(const CC_function &x,const CC_function &y,const double thresh) const {

    const real_function_3d Kx=K(x);
    const real_function_3d Ky=K(y);
    const real_function_6d fKphi0a=make_f_xy(x,CC_function(Ky,y.i,UNDEFINED),thresh);
    const real_function_6d fKphi0b=make_f_xy(CC_function(Kx,x.i,UNDEFINED),y,thresh);
    const real_function_6d fKphi0=(fKphi0a + fKphi0b);
    return fKphi0;

  }

  real_function_3d
  CC_Operators::apply_laplacian(const real_function_3d &x) const {
    // make gradient operator for new k and with new thresh
    size_t high_k=8;
    double high_thresh=1.e-6;
    std::vector < std::shared_ptr<Derivative<double, 3> > > gradop(3);
    for(std::size_t d=0; d < 3; ++d){
      gradop[d].reset(new Derivative<double, 3>(world,d,FunctionDefaults<3>::get_bc(),Function<double, 3>(),Function<double, 3>(),high_k));
    }

    // project the function to higher k grid
    real_function_3d f=project(x,high_k);
    f.set_thresh(high_thresh);
    f.refine();

    // apply laplacian
    real_function_3d empty=real_factory_3d(world);
    real_function_3d laplace_f=project(empty,high_k);
    laplace_f.set_thresh(high_thresh);
    for(size_t i=0; i < gradop.size(); i++){
      real_function_3d tmp=(*gradop[i])(f);
      real_function_3d tmp2=(*gradop[i])(tmp);
      laplace_f+=tmp2;
    }

    // project laplace_f back to the normal grid
    real_function_3d result=project(laplace_f,FunctionDefaults<3>::get_k());
    result.set_thresh(FunctionDefaults<3>::get_thresh());

    // debug and failsafe: make inverse of laplacian and apply
    real_convolution_3d G=BSHOperator<3>(world,0.0,parameters.lo,parameters.thresh_bsh_3D);
    real_function_3d Gresult=-1.0 * G(result);
    real_function_3d difference=x - Gresult;
    double diff=difference.norm2();
    plot_plane(world,difference,"Laplacian_error_iteration_" + stringify(performance_D.current_iteration));
    if(world.rank() == 0) std::cout << "Apply Laplace:\n" << "||x - G(Laplace(x))||=" << diff << std::endl;
    if(diff > FunctionDefaults<6>::get_thresh()) warning("Laplacian Error above 6D thresh");

    return result;
  }

  vecfuncT
  CC_Operators::apply_F(const CC_vecfunction &x) const {
    vecfuncT result;
    for(const auto itmp : x.functions){
      const CC_function& xi=itmp.second;
      result.push_back(apply_F(xi));
    }
    return result;
  }

  real_function_3d
  CC_Operators::apply_F(const CC_function &x) const {

    if(x.type == HOLE){
      return get_orbital_energies()[x.i] * x.function;
    }else if(x.type == PARTICLE and not current_singles_potential_gs.empty()){
      const real_function_3d singles_potential=current_singles_potential_gs[x.i - parameters.freeze];
      return (get_orbital_energies()[x.i] * x.function - singles_potential);
    }else if(x.type == MIXED and not current_singles_potential_gs.empty()){
      const real_function_3d singles_potential=current_singles_potential_gs[x.i - parameters.freeze];
      return (get_orbital_energies()[x.i] * x.function - singles_potential);     // for mixed: eps(i)*x.i = epsi*(moi + taui)
    }else{
      real_function_3d refined_x=copy(x.function).refine();
      // kinetic part
      CC_Timer T_time(world,"apply_T");
      std::vector < std::shared_ptr<real_derivative_3d> > gradop;
      gradop=gradient_operator<double, 3>(world);
      real_function_3d laplace_x=apply_laplacian(x.function);
      real_function_3d Tx=laplace_x.scale(-0.5).truncate();
      T_time.info();

      CC_Timer J_time(world,"apply_J");
      real_function_3d hartree_potential= real_function_3d(world);
      for(const auto& tmpk :mo_ket_.functions) hartree_potential+=g12(mo_bra_(tmpk.first),tmpk.second);
      const real_function_3d Jx = hartree_potential*x.function;
      J_time.info();

      CC_Timer K_time(world,"apply_K");
      const real_function_3d Kx=K(x);
      K_time.info();

      CC_Timer U_time(world,"apply_U");
      real_function_3d U2x=(nemo.nuclear_correlation->U2() * x.function).truncate();
      real_function_3d U1x=real_factory_3d(world);
      for(size_t axis=0; axis < 3; axis++){
	const real_function_3d U1_axis=nemo.nuclear_correlation->U1(axis);
	const real_function_3d dx=(*gradop[axis])(x.function);
	U1x+=(U1_axis * dx).truncate();
      }
      U_time.info();

      return (Tx + 2.0 * Jx - Kx + U2x + U1x).truncate();
    }
    error("apply_F: should not end up here");
    return real_factory_3d(world);
  }

  /// swap particles 1 and 2

  /// param[in]	f	a function of 2 particles f(1,2)
  /// return	the input function with particles swapped g(1,2) = f(2,1)
  real_function_6d
  CC_Operators::swap_particles(const real_function_6d& f) const {
   // CC_Timer timer_swap(world,"swap particles");
    // this could be done more efficiently for SVD, but it works decently
    std::vector<long> map(6);
    map[0]=3;
    map[1]=4;
    map[2]=5;     // 2 -> 1
    map[3]=0;
    map[4]=1;
    map[5]=2;     // 1 -> 2
   // timer_swap.info();
    return mapdim(f,map);
  }

  // Calculate the CC2 energy equation which is
  // \omega = \sum_{ij} 2<ij|g|\tau_{ij}> - <ij|g|\tau_{ji}> + 2 <ij|g|\tau_i\tau_j> - <ij|g|\tau_j\tau_i>
  // with \tau_{ij} = u_{ij} + Q12f12|ij> + Q12f12|\tau_i,j> + Q12f12|i,\tau_j> + Q12f12|\tau_i\tau_j>
  double
  CC_Operators::get_CC2_correlation_energy() const {
    MADNESS_EXCEPTION("get_cc2_correlation_energy not implemented yet",1);
    return 0.0;
  }


  // E = <x_i|S2b_i> + <x_i|S2c_i>
  //   = 2<x_ik|g|tauik> - <kx_i|g|tauik> + <x_i|S2c_i>
  double CC_Operators::compute_cispd_energy(const Pairs<CC_Pair> &u, const Pairs<CC_Pair> mp2_doubles, const CC_vecfunction x){
    remove_stored_singles_potentials();
    const CC_vecfunction active_mo = make_t_intermediate(x);
    const vecfuncT s2b_u = S2b_u_part(u,active_mo); // amo is only for size
    const vecfuncT s2b_r = add(world,S2b_reg_part(x,active_mo),S2b_reg_part(active_mo,x));
    const vecfuncT s2c_u = S2c_u_part(u,active_mo); // amo is only for size
    const vecfuncT s2c_r = add(world,S2c_reg_part(x,active_mo),S2c_reg_part(active_mo,x));
    const vecfuncT brax = mul(world,nemo.nuclear_correlation->square(),x.get_vecfunction());
    const Tensor<double> s2but=inner(world,brax,s2b_u);
    const Tensor<double> s2brt=inner(world,brax,s2b_r);
    const Tensor<double> s2cut=inner(world,brax,s2c_u);
    const Tensor<double> s2crt=inner(world,brax,s2c_r);
    const double ws2bu = s2but.sum();
    const double ws2br = s2brt.sum();
    const double ws2cu = s2cut.sum();
    const double ws2cr = s2crt.sum();
    const double tmp = ws2bu + ws2br + ws2cu + ws2cr;
    remove_stored_singles_potentials();
    const double constant_part = compute_cispd_energy_constant_part(mp2_doubles,x);
    remove_stored_singles_potentials();
    const double result = constant_part + tmp;
    if(world.rank()==0){
      std::cout << " CIS(D) Energy Calculation:\n";
      std::cout << " constant part: <xi|S4(MP2,CIS)> =" << constant_part << "\n";
      std::cout << " doubles  part: <xi|S2(CIS(D))>  =" << tmp << "\n";
      std::cout << " result = " << result << "\n";
    }
    return result;
  }

  // return <x|S4(u,x)>
  double CC_Operators::compute_cispd_energy_constant_part(const Pairs<CC_Pair> &u, const CC_vecfunction x)const{
    const CC_vecfunction active_mo = make_t_intermediate(x);
    const vecfuncT S4a = add(world,S4a_u_part(u,x),S4a_reg_part(active_mo,active_mo,x));
    const vecfuncT S4b = add(world,S4b_u_part(u,x),S4b_reg_part(active_mo,active_mo,x));
    const vecfuncT S4c = add(world,S4c_u_part(u,x),S4c_reg_part(active_mo,active_mo,x));
    const vecfuncT tmp = add(world,S4a,S4b);
    vecfuncT V = add(world,tmp,S4c);
    V=apply_Q(V,"cispd-energy-constant-part");
    const vecfuncT brax = mul(world,nemo.nuclear_correlation->square(),x.get_vecfunction());
    const Tensor<double> xV = inner(world,brax,V);
    return xV.sum();
  }

  double CC_Operators::compute_cc2_pair_energy(const CC_Pair &u,const CC_function &taui,const CC_function &tauj) const {
    calctype type;
    if(taui.function.is_initialized() and tauj.function.is_initialized()){
      type=CC2_;
      MADNESS_ASSERT(u.i == taui.i);
      MADNESS_ASSERT(u.j == tauj.i);
    }else{
      type=MP2_;
    }

    double omega=0.0;
    const size_t i=u.i;
    const size_t j=u.j;
    double u_part=0.0;
    double ij_part=0.0;
    double mixed_part=0.0;
    double titj_part=0.0;
    double singles_part=0.0;
    double tight_thresh=parameters.tight_thresh_6D;
    // Contribution from u itself, we will calculate <uij|g|ij> instead of <ij|g|uij> and then just make the inner product (see also mp2.cc)
    {
      real_function_6d coulomb=TwoElectronFactory(world).dcut(tight_thresh);
      real_function_6d g_ij=CompositeFactory<double, 6, 3>(world).particle1(copy(mo_bra_(i).function)).particle2(copy(mo_bra_(j).function)).g12(coulomb).thresh(tight_thresh);

      real_function_6d g_ji=CompositeFactory<double, 6, 3>(world).particle1(copy(mo_bra_(j).function)).particle2(copy(mo_bra_(i).function)).g12(coulomb).thresh(tight_thresh);
      const double uij_g_ij=inner(u.function,g_ij);
      const double uij_g_ji=inner(u.function,g_ji);     // =uji_g_ij
      u_part=2.0 * uij_g_ij - uij_g_ji;
    }
    // Contribution from the mixed f12(|\tau_i,j>+|i,\tau_j>) part
    if(type == CC2_){
      mixed_part+=2.0 * make_ijgQfxy(u.i,u.j,mo_ket_(i),tauj);
      mixed_part+=2.0 * make_ijgQfxy(u.i,u.j,taui,mo_ket_(j));
      mixed_part-=make_ijgQfxy(u.j,u.i,mo_ket_(i),tauj);
      mixed_part-=make_ijgQfxy(u.j,u.i,taui.function,mo_ket_(j));
    }
    // Contribution from the f12|ij> part, this should be calculated in the beginning
    {
      ij_part+=(2.0 * u.ij_gQf_ij - u.ji_gQf_ij);
    }
    // Contribution from the f12|\tau_i\tau_j> part
    if(type == CC2_){
      titj_part+=2.0 * make_ijgQfxy(u.i,u.j,taui,tauj);
      titj_part-=make_ijgQfxy(u.i,u.j,tauj,taui);
    }
    // Singles Contribution
    if(type == CC2_){
      // I should use intermediates later because the t1 integrals are also needed for the CC2 potential
      //omega += 2.0*intermediates_.get_integrals_t1()(u.i,u.j,u.i,u.j); //<ij|g|\taui\tauj>
      singles_part+=2.0 * make_integral(u.i,u.j,taui.function,tauj.function,g12_);
      //omega -= intermediates_.get_integrals_t1()(u.i,u.j,u.j,u.i);     //<ij|g|\tauj\taui>
      singles_part-=make_integral(u.i,u.j,tauj.function,taui.function,g12_);
    }

    omega=u_part + ij_part + mixed_part + titj_part + singles_part;
    if(world.rank() == 0){
      std::cout << "\n\nEnergy Contributions to the" + assign_name(type) + "-correlation energy of pair " << i << j << "\n";
      std::cout << std::fixed << std::setprecision(parameters.output_prec);
      std::cout << "from   |u" << i << j << "            |: " << u_part << "\n";
      std::cout << "from Qf|HH" << i << j << "           |: " << ij_part << "\n";
      if(type == CC2_) std::cout << "from Qf|HP" << i << j << "           |: " << mixed_part << "\n";
      if(type == CC2_) std::cout << "from Qf|PPu" << i << j << "          |: " << titj_part << "\n";
      if(type == CC2_) std::cout << "from   |tau" << i << ",tau" << j << "|: " << singles_part << "\n";
      std::cout << "all together = " << omega << "\n";
      std::cout << "\n\n";
    }
    return omega;
  }


  // make shure that the correct bra elements are inserted ( bra_x = x*R^2)
  double CC_Operators::make_ijgQfxy(const  CC_function &i, const  CC_function &j, const CC_function &x, const CC_function &y) const{
    // part1: No projector, <ij|gf|xy>
    const real_function_3d ix = i.function*x.function;
    const real_function_3d jy = j.function*y.function;
    const real_function_3d gfix = apply_gf(ix);
    const double part1 = jy.inner(gfix);

    // part2: O1-part, <ij|gO1f|xy> = <j|igm*mfx|y>
    // part3: O2-part, <ij|gO3f|xy> = <i|jgm*mfy|x>
    // make part 4 from part2 ims
    // part4: O12-part,<ij|gO12f|xy>= ijgmn*mnfyx = <j|igm|n>*<n|mfx|y>
    double part2 = 0.0;
    double part3 = 0.0;
    double part4 = 0.0;
    for(const auto& mtmp:mo_ket_.functions){
      const CC_function& mket =mo_ket_(mtmp.first);
      const CC_function& mbra =mo_bra_(mtmp.first);
      const real_function_3d igm = g12(mket,i); //g12(i,mket); this way intermediates can be used
      const real_function_3d mfx = f12(mbra,x);
      part2 += jy.inner(igm*mfx);
      part3 += ix.inner(g12(j,mket)*f12(mbra,y));
      for(const auto& ntmp:mo_ket_.functions){
	const CC_function &nket = mo_ket_(ntmp.first);
	const CC_function &nbra = mo_bra_(ntmp.first);
	const real_function_3d jn = j.function*nket.function;
	const real_function_3d ny = nbra.function*y.function;
	part4 += jn.inner(igm)*ny.inner(mfx);
      }
    }
//    if(world.rank()==0){
//      std::cout << "<" << i.name() << j.name() <<"|gQf|" << x.name() << y.name() << "\n";
//      std::cout << "part1=" << part1 << "\n";
//      std::cout << "part2=" << part2 << "\n";
//      std::cout << "part3=" << part3 << "\n";
//      std::cout << "part4=" << part4 << "\n";
//    }
    return part1-part2-part3+part4;
  }

  /// General Function to make the intergral <ij|gQf|xy>
  double
  CC_Operators::make_ijgQfxy(const size_t &i,const size_t &j,const CC_function &x,const CC_function &y) const {
    // convenience
    const real_function_3d brai=mo_bra_(i).function;
    const real_function_3d braj=mo_bra_(j).function;
    // part 1, no projector: <ij|gf|xy>
    const real_function_3d jy=msparse(braj,y.function);
    const real_function_3d ix=msparse(brai,x.function);
    const real_function_3d jgfy=apply_gf(jy);
    const double part1=ix.inner(jgfy);
    // part 2, projector on particle 1 <j|igm*mfx|y> = jy.inner(igm*mfx)
    double part2=0.0;
    for(const auto& mtmp : mo_ket_.functions){
      const CC_function& mom=mtmp.second;
      const size_t m=mtmp.first;
      const real_function_3d igm=g12(mo_bra_(i),mom);
      const real_function_3d mfx=f12(mo_bra_(m),x);
      part2-=jy.inner(msparse(igm.reconstruct(),mfx.reconstruct()));
    }
    // part3, projector on particle 2 <i|jgn*nfy|x>
    double part3=0.0;
    for(const auto& ntmp : mo_ket_.functions){
      const CC_function& mon=ntmp.second;
      const size_t n=ntmp.first;
      const real_function_3d jgn=g12(mo_bra_(j),mon);
      const real_function_3d nfy=f12(mo_bra_(n),y);
      part3-=ix.inner(msparse(jgn,nfy));
    }
    // part4, projector on both particles <ij|g|mn><mn|f|xy>
    double part4=0.0;
    for(const auto& mtmp : mo_ket_.functions){
      const CC_function& mom=mtmp.second;
      const size_t m=mtmp.first;
      const real_function_3d igm=g12(mo_bra_(i),mom);
      const real_function_3d mfx=f12(mo_bra_(m),x);
      for(const auto& ntmp : mo_ket_.functions){
	const CC_function& mon=ntmp.second;
	const size_t n=ntmp.first;
	const real_function_3d jn=msparse(braj,mon.function);
	const real_function_3d ny=msparse(mo_bra_(n).function,y.function);
	const double ijgmn=jn.inner(igm);
	const double mnfxy=ny.inner(mfx);
	part4+=ijgmn * mnfxy;
      }
    }
//    if(world.rank()==0){
//      std::cout << "<" << mo_ket_(i).name() << mo_ket_(j).name() <<"|gQf|" << x.name() << y.name() << "\n";
//      std::cout << "part1=" << part1 << "\n";
//      std::cout << "part2=" << part2 << "\n";
//      std::cout << "part3=" << part3 << "\n";
//      std::cout << "part4=" << part4 << "\n";
//    }
    double result = part1 + part2 + part3 + part4;
    //double test = make_ijgQfxy(mo_bra_(i),mo_bra_(j),x,y);
    //std::cout << " debug ijgQfxy:\nold =" << result <<"\nnew =" << test << "\ndiff=" << test-result << std::endl;
    //if(fabs(test-result)>1.e-5) warning("WARNING FOR IJGQFXY");
    return result;

  }

  double
  CC_Operators::make_ijgfxy(const size_t &i,const size_t &j,const real_function_3d &x,const real_function_3d &y) const {
    real_function_3d jy=mo_bra_(j).function * y;
    real_function_3d ix=mo_bra_(j).function * x;
    // I12 Part:
    double ijgfxy=(ix).inner(apply_gf(jy));
    return ijgfxy;
  }
  //
  //  /// General Function to make the two electron integral <ij|g|xy>
  //  /// For Debugging -> Expensive without intermediates
  //  double
  //  CC_Operators::make_ijgxy(const size_t &i,const size_t &j,const CC_function &x,const CC_function &y) const {
  //    real_function_3d igx=g12(mo_bra_(i),x)
  //    real_function_3d jy=(mo_bra_(j).function * y).truncate();
  //    return jy.inner(igx);
  //  }

  double
  CC_Operators::make_integral(const size_t &i,const size_t &j,const CC_function &x,const CC_function&y,const optype type) const {
    //    if(x.type == HOLE){
    //      real_function_3d igx_y=(intermediates_.get_EX(i,x.i) * y.function).truncate();
    //      return mo_bra_(j).function.inner(igx_y);
    //    }else if(x.type == PARTICLE){
    //      if(y.type == HOLE){
    //	real_function_3d jgy_x=(intermediates_.get_EX(j,y.i) * x.function);
    //	return mo_bra_(i).function.inner(jgy_x);
    //      }else if(y.type == PARTICLE){
    //	real_function_3d jgy_x=(intermediates_.get_pEX(j,y.i) * x.function);
    //	return mo_bra_(i).function.inner(jgy_x);
    //      }
    //    }else if(x.type == MIXED or y.type == MIXED){
    //      real_function_3d igx=((*poisson)(mo_bra_(i).function * x.function)).truncate();
    //      double result=mo_bra_(j).function.inner(igx * y.function);
    //      return result;
    //    }else if(x.type == UNDEFINED or y.type == UNDEFINED){
    //      real_function_3d igx=((*poisson)(mo_bra_(i).function * x.function)).truncate();
    //      double result=mo_bra_(j).function.inner(igx * y.function);
    //      return result;
    //    }else{
    //      error("ERROR in make_integrals ... should not end up here");
    //      return 0.0;
    //    }
    //    error("ERROR in make_integrals ... should not end up here");
    //    return 0.0;
    switch(type){
      case g12_ : {
	const real_function_3d igx= g12(mo_bra_(i),x);
	const real_function_3d jy=mo_bra_(j).function*y.function;
	return jy.inner(igx);
	break;
      }
      case f12_: {
	const real_function_3d igx= f12(mo_bra_(i),x);
	const real_function_3d jy=mo_bra_(j).function*y.function;
	return jy.inner(igx);
	break;
      }
    }
    MADNESS_EXCEPTION("make_integral: how did we get here?",1);
  }

  /// General Function to make two electron integrals with pair functions (needed for energy)
  double
  CC_Operators::make_ijgu(const size_t &i,const size_t &j,const CC_Pair &u) const {
    return make_ijgu(i,j,u.function);
  }
  double
  CC_Operators::make_ijgu(const size_t &i,const size_t &j,const real_function_6d &u) const {
    real_function_6d g=TwoElectronFactory(world).dcut(parameters.lo);
    real_function_6d ij_g=CompositeFactory<double, 6, 3>(world).particle1(copy(mo_bra_(i).function)).particle2(copy(mo_bra_(j).function)).g12(g);

    // compute < ij | g12 | u >
    const double ij_g_u=inner(u,ij_g);
    return ij_g_u;
  }

  double
  CC_Operators::make_ijgu(const CC_function &phi_i,const CC_function &phi_j,const CC_Pair &u) const {
    const real_function_3d brai=nemo.nuclear_correlation->square()*phi_i.function;
    const real_function_3d braj=nemo.nuclear_correlation->square()*phi_j.function;
    real_function_6d g=TwoElectronFactory(world).dcut(parameters.lo);
    real_function_6d ij_g=CompositeFactory<double, 6, 3>(world).particle1(copy(brai)).particle2(copy(braj)).g12(g);
    // compute < ij | g12 | u >
    const double ij_g_u=inner(u.function,ij_g);
    return ij_g_u;
  }

  /// General Function to make two electorn integrals with pair function and orbitals and the BSH Operator (needed for gf = g - bsh)
  double
  CC_Operators::make_ijGu(const size_t &i,const size_t &j,const CC_Pair &u) const {
    real_function_6d G=TwoElectronFactory(world).BSH().gamma(corrfac.gamma()).dcut(parameters.lo);
    double bsh_prefactor=4.0 * constants::pi;
    real_function_6d ij_G=CompositeFactory<double, 6, 3>(world).particle1(copy(mo_bra_(i).function)).particle2(copy(mo_bra_(j).function)).g12(G);

    // compute < ij | g12 | u >
    const double ij_G_u=inner(u.function,ij_G);
    return bsh_prefactor * ij_G_u;
  }

  real_function_3d
  CC_Operators::convolute_x_Qf_yz(const CC_function &x,const CC_function &y,const CC_function &z) const {
    // !!!!!!!!!!!!! bra element: make shure x is actually R2*x when given to the function
    real_function_3d xfz=f12(msparse(x.function,z.function));
    xfz.truncate();
    xfz.reconstruct();
    const real_function_3d xfz_y=msparse(xfz,y.function).truncate();
    const real_function_3d part1=msparse(xfz,y.function);

    real_function_3d part2=real_factory_3d(world);
    real_function_3d part3tmp=real_factory_3d(world);
    real_function_3d part4=real_factory_3d(world);
    for(const auto& mtmp : mo_ket_.functions){
      const CC_function& mom=mtmp.second;
      const double mxfyz=mo_bra_(mom).function.inner(xfz_y);
      part2-=mxfyz * mom.function;

      const double xm=x.function.inner(mom.function);

      const real_function_3d mfz=f12(mo_bra_(mom),z);
      const real_function_3d mfz_y=msparse(mfz,y.function);

      part3tmp-=xm * mfz;

      for(const auto& ntmp : mo_ket_.functions){
	const CC_function& mon=ntmp.second;
	const double nmfyz=mo_bra_(mon).function.inner(mfz_y);
	part4+=xm * nmfyz * mon.function;
      }

    }
    const real_function_3d part3=msparse(part3tmp,y.function);
    real_function_3d result=part1 + part2 + part3 + part4;
    result.truncate();
    return result;
  }

  /// apply the operator gf = 1/(2\gamma)*(Coulomb - 4\pi*BSH_\gamma)
  /// works only if f = (1-exp(-\gamma*r12))/(2\gamma)
  real_function_3d
  CC_Operators::apply_gf(const real_function_3d &f) const {
    double bsh_prefactor=4.0 * constants::pi;
    double prefactor=1.0 / (2.0 * corrfac.gamma());
    return prefactor * (g12(f) - bsh_prefactor * (*fBSH)(f)).truncate();
  }

  real_function_6d
  CC_Operators::make_xy(const CC_function &x,const CC_function &y) const {
    double thresh=guess_thresh(x,y);
    if(thresh < parameters.thresh_6D) thresh=parameters.tight_thresh_6D;
    else thresh=parameters.thresh_6D;
    CC_Timer timer(world,"Making |" + x.name() + "," + y.name() + "> with 6D thresh=" + stringify(thresh));
    real_function_6d xy=CompositeFactory<double, 6, 3>(world).particle1(copy(x.function)).particle2(copy(y.function)).thresh(thresh);
    xy.fill_tree().truncate().reduce_rank();
    timer.info();
    return xy;
  }

  real_function_6d
  CC_Operators::make_f_xy(const CC_function &x,const CC_function &y,const double thresh) const {
    CC_Timer timer(world,"Making f|" + x.name() + "," + y.name() + "> with 6D thresh=" + stringify(thresh));
    real_function_6d fxy=CompositeFactory<double, 6, 3>(world).g12(corrfac.f()).particle1(copy(x.function)).particle2(copy(y.function)).thresh(thresh);
    fxy.fill_tree().truncate().reduce_rank();
    timer.info();
    return fxy;
  }

  real_function_6d
  CC_Operators::make_f_xy_screened(const CC_function &x,const CC_function &y,const real_convolution_6d &screenG) const {
    double thresh=guess_thresh(x,y);
    if(thresh < parameters.thresh_6D) thresh=parameters.tight_thresh_6D;
    else thresh=parameters.thresh_6D;
    CC_Timer timer(world,"Making f|" + x.name() + "," + y.name() + "> with 6D thresh=" + stringify(thresh));
    real_function_6d fxy=CompositeFactory<double, 6, 3>(world).g12(corrfac.f()).particle1(copy(x.function)).particle2(copy(y.function)).thresh(thresh);
    fxy.fill_tree(screenG).truncate().reduce_rank();
    timer.info();
    return fxy;
  }

  /// Calculation is done in 4 steps over: Q12 = 1 - O1 - O2 + O12
  /// 1. <x|f12|z>*|y>
  /// 2. -\sum_m <x|m> <m|f12|z>*|y> = -(\sum_m <x|m> <m|f12|z>) *|y>
  /// 3. -\sum_n <nx|f12|zy> * |n>
  /// 4. +\sum_{mn} <x|n> <mn|f12|yz> * |m>
  /// Metric from nuclear cusp is not taken into account -> give the right bra elements to the function
  //real_function_3d CC_Operators::convolute_x_Qf_yz(const real_function_3d &x, const real_function_3d &y, const real_function_3d &z)const{
  //	// make intermediates
  //	vecfuncT moz = mul(world,z,mo_bra_);
  //	vecfuncT moy = mul(world,y,mo_bra_);
  //	// Do <x|f12|z>*|y>
  //	real_function_3d part1_tmp = (*f12op)(x*z);
  //	real_function_3d part1 = (part1_tmp*y).truncate();
  //	// Do -\sum_m <x|m> <m|f12|z>*|y>
  //	real_function_3d part2 = real_factory_3d(world);
  //	{
  //		Tensor<double> xm = inner(world,x,mo_ket_);
  //		vecfuncT f12mz = apply(world,*f12op,moz);
  //		real_function_3d tmp = real_factory_3d(world);
  //		for(size_t m=0;m<mo_bra_.size();m++) tmp += xm[m]*f12mz[m];
  //		part2 = tmp*y;
  //	}
  //	// Do -\sum_n <nx|f12|zy> * |n> |  <nx|f12|zy> = <n| xf12y |z>
  //	real_function_3d part3 = real_factory_3d(world);
  //	{
  //		real_function_3d xf12y = (*f12op)((x*y).truncate());
  //		for(size_t n=0;n<mo_bra_.size();n++){
  //			double nxfzy = xf12y.inner(mo_bra_[n]*z);
  //			part3 += nxfzy*mo_ket_[n];
  //		}
  //	}
  //	// Do +\sum_{mn} <x|n> <mn|f12|yz> * |m>
  //	real_function_3d part4 = real_factory_3d(world);
  //	{
  //		Tensor<double> xn = inner(world,x,mo_ket_);
  //		vecfuncT nf12z = apply(world,*f12op,moz);
  //		Tensor<double> mnfyz = inner(world,moy,nf12z);
  //		for(size_t m=0;m<mo_bra_.size();m++){
  //			for(size_t n=0;n<mo_ket_.size();n++){
  //				part4 += xn(n)*mnfyz(m,n)*mo_ket_[m];
  //			}
  //		}
  //	}
  //	real_function_3d result = part1 - part2 - part3 + part4;
  //	result.truncate();
  //
  //	if(parameters.debug){
  //		CC_Timer function_debug(world,"Debug-Time for <k|Qf|xy>");
  //		real_function_6d test_tmp = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(y)).particle2(copy(z));
  //		test_tmp.fill_tree().truncate().reduce_rank();
  //		real_function_6d test_1 = copy(test_tmp);
  //		real_function_6d test_4 = copy(Q12(test_tmp));
  //		real_function_3d test1 = test_1.project_out(x,1);
  //		real_function_3d test4 = test_4.project_out(x,1);
  //		double check1 = (test1 - part1).norm2();
  //		double check4 = (test4 - result).norm2();
  //		if(world.rank()==0) std::cout << std::setprecision(parameters.output_prec) << "<k|Qf|xy> debug, difference to check1 value is: " << check1 << std::endl;
  //		if(world.rank()==0) std::cout << std::setprecision(parameters.output_prec) << "<k|Qf|xy> debug, difference to check4 value is: " << check4 << std::endl;
  //		if(check1 > FunctionDefaults<6>::get_thresh()) warning("<k|Qf|xy> check1 failed");
  //		if(check4 > FunctionDefaults<6>::get_thresh()) warning("<k|Qf|xy> check4 failed");
  //		function_debug.info();
  //	}
  //
  //	return result;
  //}

  void
  CC_Operators::test_singles_potential() {

    output_section("Singles Potential Consistency Check with r*|i> singles and Q12f12|ij> doubles");
    // make test singles from mos: |taui> = r*|i>
    // make test doubles from mos: |uij>  = Q12f12|titj>
    real_function_3d r=real_factory_3d(world).f(f_r);
    vecfuncT singles_tmp;
    for(size_t i=parameters.freeze; i < mo_ket_.size(); i++){
      real_function_3d tmp=r * mo_ket_(i).function;
      tmp=projector_Q(tmp);
      double norm=tmp.norm2();
      tmp.scale(1.0 / norm);
      tmp.scale(0.5);
      tmp.print_size("TestSingle: r|" + stringify(i) + ">");
      singles_tmp.push_back(tmp);
    }
    CC_vecfunction singles(singles_tmp,PARTICLE,parameters.freeze);
    update_intermediates(singles);


    // test ccs
    output("\n\n Testing CCS Potential\n");
    vecfuncT ccs=ccs_potential(singles);
    vecfuncT new_ccs=ccs_potential_new(singles);
    for(size_t i=0; i < ccs.size(); i++){
      std::cout << "||ccs-potential    ||_" << i << "=" << (ccs[i]).norm2() << std::endl;
      std::cout << "||ccs-potential_new||_" << i << "=" << (new_ccs[i]).norm2() << std::endl;
      std::cout << "||difference       ||_" << i << "=" << (ccs[i] - new_ccs[i]).norm2() << std::endl;
    }
    ccs=apply_Q(ccs,"ccs");
    new_ccs=apply_Q(new_ccs,"new_ccs");
    std::cout << "Applied Q" << std::endl;
    for(size_t i=0; i < ccs.size(); i++){
      std::cout << "||ccs-potential    ||_" << i << "=" << (ccs[i]).norm2() << std::endl;
      std::cout << "||ccs-potential_new||_" << i << "=" << (new_ccs[i]).norm2() << std::endl;
      std::cout << "||difference       ||_" << i << "=" << (ccs[i] - new_ccs[i]).norm2() << std::endl;
    }
    vecfuncT trunc_ccs = copy(world,ccs);
    truncate(world,trunc_ccs);
    vecfuncT trunc_ccs_new = copy(world,new_ccs);
    truncate(world,trunc_ccs_new);
    for(size_t i=0; i < ccs.size(); i++){
      std::cout << "||ccs_truncation_error||_" << i << "=" << (ccs[i]-trunc_ccs[i]).norm2() << std::endl;
      std::cout << "||ccs-new_trunca_error||_" << i << "=" << (new_ccs[i]-trunc_ccs_new[i]).norm2() << std::endl;
    }


    Pairs<CC_Pair> doubles=make_reg_residues(singles);



    CC_data dummy;

    output("\n\n Checking u-parts and r-parts of singles potentials with doubles\n\n");
    const potentialtype_s u_parts_tmp[]={pot_S4a_u_, pot_S4b_u_, pot_S4c_u_, pot_S2b_u_, pot_S2c_u_};
    const potentialtype_s r_parts_tmp[]={pot_S4a_r_, pot_S4b_r_, pot_S4c_r_, pot_S2b_r_, pot_S2c_r_};
    std::vector < std::pair<std::string, double> > results;
    for(size_t pot=0; pot < 5; pot++){
      const potentialtype_s current_u=u_parts_tmp[pot];
      const potentialtype_s current_r=r_parts_tmp[pot];
      const std::string name=assign_name(current_u);
      double largest_error=0.0;
      output("\n\nConsistency Check of Singles Potential " + assign_name(current_u) + " with " + assign_name(current_r));
      const vecfuncT u=potential_singles(doubles,singles,current_u);
      const vecfuncT r=potential_singles(doubles,singles,current_r);
      const vecfuncT diff=sub(world,u,r);
      const double normdiff=norm2(world,u) - norm2(world,r);
      if(world.rank() == 0) std::cout << std::setw(20) << "||" + assign_name(current_u) + "||-||" + assign_name(current_r) + "||" << std::setfill(' ') << "=" << normdiff << std::endl;
      for(const auto d : diff){
	const double norm=d.norm2();
	if(norm > largest_error) largest_error=norm;
	if(world.rank() == 0) std::cout << std::setw(20) << "||" + assign_name(current_u) + "-" + assign_name(current_r) + "||" << std::setfill(' ') << "=" << norm << std::endl;
      }
      results.push_back(std::make_pair(name,largest_error));
      if(current_u == pot_S2b_u_){
	output("Making Integration Test for S2b potential:");
	// integrate the s2b potential against a function which not in the hole space = \sum_k 2<X,k|g|uik> - <k,X|g|uik>, with X=QX
	real_function_3d X=real_factory_3d(world);
	for(const auto&s : singles.functions){
	  X+=s.second.function;
	}
	X=projector_Q(X);
	X=X * nemo.nuclear_correlation->square();
	Tensor<double> xs2b=inner(world,X,u);
	Tensor<double> xs2b_reg=inner(world,X,r);
	std::vector<double> control_6D;
	for(auto& itmp : singles.functions){
	  const size_t i=itmp.first;
	  double controli_6D=0.0;
	  for(auto& ktmp : singles.functions){
	    const size_t k=ktmp.first;
	    real_function_6d g=TwoElectronFactory(world).dcut(parameters.lo);
	    real_function_6d Xk_g=CompositeFactory<double, 6, 3>(world).particle1(copy(X)).particle2(copy(mo_bra_(k).function)).g12(g);
	    real_function_6d g2=TwoElectronFactory(world).dcut(parameters.lo);
	    real_function_6d kX_g=CompositeFactory<double, 6, 3>(world).particle1(copy(mo_bra_(k).function)).particle2(copy(X)).g12(g2);
	    const double tmp_6D=2.0 * Xk_g.inner(get_pair_function(doubles,i,k)) - kX_g.inner(get_pair_function(doubles,i,k));
	    controli_6D+=tmp_6D;
	  }
	  control_6D.push_back(controli_6D);
	}
	for(size_t i=0; i < control_6D.size(); i++){
	  const double diff=xs2b(i) - control_6D[i];
	  const double diff2=xs2b_reg(i) - control_6D[i];
	  std::cout << "||diffu||_" << i << "=" << fabs(diff) << std::endl;
	  std::cout << "||diffr||_" << i << "=" << fabs(diff2) << std::endl;
	  if(fabs(diff) > FunctionDefaults<6>::get_thresh()) warning("Integration Test of S2b failed !!!!!");
	}
      }
    }

    bool problems=false;
    for(const auto res : results){
      std::string status="... passed!";
      if(res.second > FunctionDefaults<6>::get_thresh()){
	status="... failed!";
	problems=true;
      }
      if(world.rank() == 0) std::cout << std::setw(10) << res.first << status << " largest error was " << res.second << std::endl;
    }
    if(problems) warning("Inconsistencies in Singles Potential detected!!!!");
    else output("Singles Potentials seem to be consistent");
    output("\n\n Ending Testing Section\n\n");
  }

}
