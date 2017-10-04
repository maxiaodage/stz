/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "stzHEVPFlowRatePowerLawJ2.h"

template<>
InputParameters validParams<stzHEVPFlowRatePowerLawJ2>()
{
  InputParameters params = validParams<HEVPFlowRateUOBase>();
  params.addParam<Real>("reference_flow_rate",0.001,"Reference flow rate for rate dependent flow");
  params.addParam<Real>("flow_rate_exponent",10.0,"Power law exponent in flow rate equation");
  params.addParam<Real>("flow_rate_tol", 1e5, "Tolerance for flow rate");
  params.addParam<Real>("tol_plastic", 1e2, "Tolerance for Plasticity to kick in");

  // Adding Pressure
  params.addParam<Real>("grain_pressure",25e6,"Grain pressure");
  // Adding Chi
  params.addParam<Real>("effective_temp",0.12,"Effective Temperature");
  // Adding time scale tau
  params.addParam<Real>("t_tau",8e-7,"Time scale");
  // Density
  params.addParam<Real>("rho",1600,"Density");
  // Biot
  params.addParam<Real>("biot",0.0,"Biot coefficient");


   // Adding Coupling Chiv
  params.addCoupledVar("efftemp",0.000000001,"effective temperature") ;
  // Adding Coupling a Grain size
  params.addCoupledVar("grainsize",1e-4,"grain size") ;
  // Adding Coupling Pore pressure
  params.addCoupledVar("porepress",0.0,"Pore Pressure");

  //
   params.addClassDescription("User object to evaluate power law flow rate and flow direction based on J2");

   return params;
 }

 stzHEVPFlowRatePowerLawJ2::stzHEVPFlowRatePowerLawJ2(const InputParameters & parameters) :
     HEVPFlowRateUOBase(parameters),
     _ref_flow_rate(getParam<Real>("reference_flow_rate")),
     _flow_rate_exponent(getParam<Real>("flow_rate_exponent")),
     _flow_rate_tol(getParam<Real>("flow_rate_tol")),
     _tolplastic(getParam<Real>("tol_plastic")),
     // Adding the getting parameters Grain Pressure
      _grain_pressure(getParam<Real>("grain_pressure")),
      // Adding chi
      _chi(getParam<Real>("effective_temp")),
      //Adding time scale tau
      _t_tau(getParam<Real>("t_tau")),
      // Density
      _rho(getParam<Real>("rho")),
      // Biot
      _biot(getParam<Real>("biot")),
     // Adding coupled chiv
     _chiv(coupledValue("efftemp")),
     // Adding coupled Grainsize
    _agrain(coupledValue("grainsize")),
     // Adding coupled porepress
     _pf(coupledValue("porepress"))


 {
 }

 bool
 stzHEVPFlowRatePowerLawJ2::computeValue(unsigned int qp, Real & val) const
 {
   RankTwoTensor pk2_dev = computePK2Deviatoric(_pk2[qp], _ce[qp]);
   Real eqv_stress = computeEqvStress(pk2_dev, _ce[qp]);
   Real _pressure_eff=_grain_pressure-_biot*_pf[qp];
  //  Real _t_taunew = _agrain[qp]/std::sqrt(_grain_pressure/_rho);
   Real _t_taunew = _agrain[qp]/std::sqrt(_pressure_eff/_rho);
   Real _strengthnew;

     if (_t== _dt)
     _strengthnew=100e6;
     else
    _strengthnew=_strength[qp];

   // val = std::pow(eqv_stress/_strength[qp], _flow_rate_exponent) * _ref_flow_rate;
   //Adding Expression exp(sbar/P)
   // val = std::pow(eqv_stress/_strength[qp], _flow_rate_exponent) * _ref_flow_rate*std::exp(_grain_pressure/_grain_pressure);
   //  val = std::pow(eqv_stress/_strength[qp], _flow_rate_exponent) * _ref_flow_rate*std::exp(_grain_pressure/_grain_pressure)*(eqv_stress>_strength[qp]);
     if (eqv_stress>=_strengthnew)
  //     val = 1.0/_t_tau*std::exp(-1.0/_chiv[qp])*std::exp(eqv_stress/_grain_pressure/_chiv[qp])*(1.0-_strength[qp]/eqv_stress);
  //val = 1.0/_t_taunew*std::exp(-1.0/_chiv[qp])*std::exp(eqv_stress/_grain_pressure/_chiv[qp])*(1.0-_strength[qp]/eqv_stress);
  // Original
  val = 1.0/_t_taunew*std::exp(-1.0/_chiv[qp])*std::exp(eqv_stress/_pressure_eff/_chiv[qp])*(1.0-_strength[qp]/eqv_stress);
  // Modified
  //val = 1.0/_t_taunew*std::exp(-1.0/_chiv[qp])*(1.0-_strength[qp]/eqv_stress);


   else
    val = 0.0;
     //    std::cout<<"haha"<<val<<"\n";





     //std::cout <<"Yiedl strength\t"<< _strength[qp]<<"\n"<< "Von Mises\t"<< eqv_stress<<"\n"<<(eqv_stress<_strength[qp])<<"\n";

   // stz Dpl equation
   // val = 1/_t_tau*std::exp(-1/_chi)*std::exp(eqv_stress/_grain_pressure/_chi)*(1-_strength[qp]/(eqv_stress+1))*(eqv_stress>_strength[qp]);
   // Variable Chiv
   // val = 1/_t_tau*std::exp(-1/_chiv[qp])*std::exp(eqv_stress/_grain_pressure/_chiv[qp])*(1-_strength[qp]/(eqv_stress+1))*(eqv_stress>_strength[qp]);

   //val = 0;

    //
    //  std::cout <<"Val\t"<< val<<"\n"<< "Von Mises\t"<< eqv_stress<<"\n"<<(exp(-1/_chi))<<"\n"<<exp(eqv_stress/_grain_pressure/_chi)<<"\n";

   //val = 0;
   //std::cout <<"Yiedl strength\t"<< _strength[qp]<<"\n"<< "Von Mises\t"<< eqv_stress<<"\n"<<(eqv_stress<_strength[qp])<<"\n";
   //std::cout << _chi;
  if (val > _flow_rate_tol)
   {
// #ifdef DEBUG
//      mooseWarning("Flow rate greater than " << _flow_rate_tol << " " << val << " " << eqv_stress << " " << _strength[qp]);
// #endif
     return false;
   }
   return true;
 }

 bool
 stzHEVPFlowRatePowerLawJ2::computeDirection(unsigned int qp, RankTwoTensor & val) const
 {
   RankTwoTensor pk2_dev = computePK2Deviatoric(_pk2[qp], _ce[qp]);
   Real eqv_stress = computeEqvStress(pk2_dev, _ce[qp]);

   val.zero();
   if (eqv_stress > _strength[qp])
     {//   val = 1.5/eqv_stress * _ce[qp] * pk2_dev ;
        val = 1.5/eqv_stress*pk2_dev*_ce[qp];
      //
      //
     }

   return true;
 }

 bool
  stzHEVPFlowRatePowerLawJ2::computeDerivative(unsigned int qp, const std::string & coupled_var_name, Real & val) const
 {
   val = 0.0;

   if (_strength_prop_name == coupled_var_name)
   {
     RankTwoTensor pk2_dev = computePK2Deviatoric(_pk2[qp], _ce[qp]);
     Real eqv_stress = computeEqvStress(pk2_dev, _ce[qp]);
     Real _pressure_eff=_grain_pressure-_biot*_pf[qp];
    //  Real _t_taunew = _agrain[qp]/std::sqrt(_grain_pressure/_rho);
    Real _t_taunew = _agrain[qp]/std::sqrt(_pressure_eff/_rho);
    // Real _t_taunew = _agrain[qp]/std::sqrt(_grain_pressure/_rho);
     //  std::cout<<eqv_stress<<"\n";
     //  val = - _ref_flow_rate * _flow_rate_exponent * std::pow(eqv_stress/_strength[qp],_flow_rate_exponent)/_strength[qp];
     //stz
     if(eqv_stress>=_strength[qp])
        //   val=1.0/_t_tau*std::exp(-1.0/_chiv[qp])*std::exp(eqv_stress/_grain_pressure/_chiv[qp])*(-1.0/eqv_stress);
        //   val=1.0/_t_taunew*std::exp(-1.0/_chiv[qp])*std::exp(eqv_stress/_grain_pressure/_chiv[qp])*(-1.0/eqv_stress);
        //Original
          //val=1.0/_t_taunew*std::exp(-1.0/_chiv[qp])*std::exp(eqv_stress/_pressure_eff/_chiv[qp])*(-1.0/eqv_stress);
       // Modified
          val = 0.0;
       //  val=1.0/_t_tau*std::exp(-1.0/_chi)*std::exp(eqv_stress/_grain_pressure/_chi)*(-1.0/eqv_stress);

       //
       //  val=1/_t_tau*std::exp(-1/_chiv[qp])*std::exp(eqv_stress/_grain_pressure/_chiv[qp])*(-1/(eqv_stress+1))*(eqv_stress>_strength[qp]);
     // out put
  // std::cout <<"dirVal\t"<< val<<"\n"<<"*******************\n" ;
   }

   return true;
 }

 bool
 stzHEVPFlowRatePowerLawJ2::computeTensorDerivative(unsigned int qp, const std::string & coupled_var_name, RankTwoTensor & val) const
 {
   val.zero();

   if (_pk2_prop_name == coupled_var_name)
   {
     RankTwoTensor pk2_dev = computePK2Deviatoric(_pk2[qp], _ce[qp]);
     Real eqv_stress = computeEqvStress(pk2_dev, _ce[qp]);
     Real _pressure_eff=_grain_pressure-_biot*_pf[qp];
    //  Real _t_taunew = _agrain[qp]/std::sqrt(_grain_pressure/_rho);
    Real _t_taunew = _agrain[qp]/std::sqrt(_pressure_eff/_rho);
     //   Real dflowrate_dseqv = _ref_flow_rate * _flow_rate_exponent * std::pow(eqv_stress/_strength[qp],_flow_rate_exponent-1.0)/_strength[qp];
 //  stz
     //   Real dflowrate_dseqv = 1/_t_tau*std::exp(-1/_chiv[qp])*std::exp(eqv_stress/_grain_pressure/_chiv[qp])* \
     //	    (1/(_grain_pressure*_chiv[qp])*(1-_strength[qp]/eqv_stress)+_strength[qp]/(eqv_stress*eqv_stress))*(eqv_stress>_strength[qp]);
//
       Real dflowrate_dseqv;
      if (eqv_stress>_strength[qp])
	  //   dflowrate_dseqv = 1/_t_tau*std::exp(-1/_chiv[qp])*std::exp(eqv_stress/_grain_pressure/_chiv[qp])* \
	  (1.0/(_grain_pressure*_chiv[qp])*(1.0-_strength[qp]/eqv_stress)+_strength[qp]/(eqv_stress*eqv_stress));
  //   Original Stuff
  {
    dflowrate_dseqv = 1/_t_taunew*std::exp(-1/_chiv[qp])*std::exp(eqv_stress/_pressure_eff/_chiv[qp])* \
    (1.0/(_pressure_eff*_chiv[qp])*(1.0-_strength[qp]/eqv_stress)+_strength[qp]/(eqv_stress*eqv_stress));
  } // Modified
  //  {
      //dflowrate_dseqv = 1.0/_t_taunew*std::exp(-1.0/_chiv[qp])*(-1.0*_strength[qp])/(eqv_stress*eqv_stress);
  //  }	//  dflowrate_dseqv = 1/_t_tau*std::exp(-1/_chi)*std::exp(eqv_stress/_grain_pressure/_chi)* \
	      //  (1.0/(_grain_pressure*_chi)*(1.0-_strength[qp]/eqv_stress)+_strength[qp]/(eqv_stress*eqv_stress));

     RankTwoTensor tau = pk2_dev * _ce[qp];
     RankTwoTensor dseqv_dpk2dev;

     dseqv_dpk2dev.zero();
     if (eqv_stress > _strength[qp])
       dseqv_dpk2dev = 1.5/eqv_stress * tau * _ce[qp];

     RankTwoTensor ce_inv = _ce[qp].inverse();

     RankFourTensor dpk2dev_dpk2;
     for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
       for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
         for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
           for (unsigned int l = 0; l < LIBMESH_DIM; ++l)
           {
             dpk2dev_dpk2(i, j, k, l) = 0.0;
             if (i==k && j==l)
               dpk2dev_dpk2(i, j, k, l) = 1.0;
             dpk2dev_dpk2(i, j, k, l) -= ce_inv(i, j) * _ce[qp](k, l)/3.0;
           }
     val = dflowrate_dseqv * dpk2dev_dpk2.transposeMajor() * dseqv_dpk2dev;
   }
   return true;
 }

RankTwoTensor
 stzHEVPFlowRatePowerLawJ2::computePK2Deviatoric(const RankTwoTensor & pk2, const RankTwoTensor & ce) const
{
  return pk2 - (pk2.doubleContraction(ce) * ce.inverse())/3.0;
}

Real
stzHEVPFlowRatePowerLawJ2::computeEqvStress(const RankTwoTensor & pk2_dev, const RankTwoTensor & ce) const
{
  RankTwoTensor sdev = pk2_dev * ce;
  Real val = sdev.doubleContraction(sdev.transpose());
  return std::pow(1.5 * val, 0.5);
}
