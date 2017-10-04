/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "stzVOLFlowRate.h"

template<>
InputParameters validParams<stzVOLFlowRate>()
{
  InputParameters params = validParams<HEVPFlowRateUOBase>();
// volume dilatcny co
   // Adding Coupling Chiv
  //
   params.addClassDescription("User object to evaluate stzvolumetric deformation on J2");
   return params;
 }

 stzVOLFlowRate::stzVOLFlowRate(const InputParameters & parameters) :
     HEVPFlowRateUOBase(parameters),
     // Adding the getting parameters Grain Pressure
    // _volplc(getMaterialProperty<Real>("volmetricpe"))
    _volplc(getMaterialProperty<Real>("volmetricpe"))


 {
 }

 bool
 stzVOLFlowRate::computeValue(unsigned int qp, Real & val) const
 {
   RankTwoTensor pk2_dev = computePK2Deviatoric(_pk2[qp], _ce[qp]);
   Real eqv_stress = computeEqvStress(pk2_dev, _ce[qp]);
   Real _strengthnew;

     if (_t== _dt)
     _strengthnew=100e6;
     else
    _strengthnew=_strength[qp];
   // val = std::pow(eqv_stress/_strength[qp], _flow_rate_exponent) * _ref_flow_rate;
   //Adding Expression exp(sbar/P)
   // val = std::pow(eqv_stress/_strength[qp], _flow_rate_exponent) * _ref_flow_rate*std::exp(_grain_pressure/_grain_pressure);
   //  val = std::pow(eqv_stress/_strength[qp], _flow_rate_exponent) * _ref_flow_rate*std::exp(_grain_pressure/_grain_pressure)*(eqv_stress>_strength[qp]);
     if(eqv_stress>_strengthnew)
     {
       val =_volplc[qp];
     }
   else
    val = 0.0;

    return true;
 }

 bool
 stzVOLFlowRate::computeDirection(unsigned int qp, RankTwoTensor & val) const
 {
   RankTwoTensor pk2_dev = computePK2Deviatoric(_pk2[qp], _ce[qp]);
   Real eqv_stress = computeEqvStress(pk2_dev, _ce[qp]);
   RankTwoTensor pk2_hydro = computePK2Hydrostatic(_pk2[qp], _ce[qp]);


   val.zero();
   if (eqv_stress > _strength[qp])
     {   val.addIa(1.0);

      //
      //
     }

   return true;
 }

 bool
  stzVOLFlowRate::computeDerivative(unsigned int qp, const std::string & coupled_var_name, Real & val) const
 {
   val = 0.0;

   if (_strength_prop_name == coupled_var_name)
   {
     RankTwoTensor pk2_dev = computePK2Deviatoric(_pk2[qp], _ce[qp]);
     RankTwoTensor pk2_hydo = computePK2Hydrostatic(_pk2[qp], _ce[qp]);
     Real eqv_stress = computeEqvStress(pk2_dev, _ce[qp]);

    //  Real _t_taunew = _agrain[qp]/std::sqrt(_grain_pressure/_rho);
    // Real _t_taunew = _agrain[qp]/std::sqrt(_grain_pressure/_rho);
     //  std::cout<<eqv_stress<<"\n";
     //  val = - _ref_flow_rate * _flow_rate_exponent * std::pow(eqv_stress/_strength[qp],_flow_rate_exponent)/_strength[qp];
     //stz
     if(eqv_stress>=_strength[qp])
        //   val=1.0/_t_tau*std::exp(-1.0/_chiv[qp])*std::exp(eqv_stress/_grain_pressure/_chiv[qp])*(-1.0/eqv_stress);
        //   val=1.0/_t_taunew*std::exp(-1.0/_chiv[qp])*std::exp(eqv_stress/_grain_pressure/_chiv[qp])*(-1.0/eqv_stress);
          val=0.0;


       //  val=1.0/_t_tau*std::exp(-1.0/_chi)*std::exp(eqv_stress/_grain_pressure/_chi)*(-1.0/eqv_stress);

       //
       //  val=1/_t_tau*std::exp(-1/_chiv[qp])*std::exp(eqv_stress/_grain_pressure/_chiv[qp])*(-1/(eqv_stress+1))*(eqv_stress>_strength[qp]);
     // out put
  // std::cout <<"dirVal\t"<< val<<"\n"<<"*******************\n" ;
   }

   return true;
 }

 bool
 stzVOLFlowRate::computeTensorDerivative(unsigned int qp, const std::string & coupled_var_name, RankTwoTensor & val) const
 {
   val.zero();


   if (_pk2_prop_name == coupled_var_name)
   {
     RankTwoTensor pk2_dev = computePK2Deviatoric(_pk2[qp], _ce[qp]);
     Real eqv_stress = computeEqvStress(pk2_dev, _ce[qp]);
    //  Real _t_taunew = _agrain[qp]/std::sqrt(_grain_pressure/_rho);
     //   Real dflowrate_dseqv = _ref_flow_rate * _flow_rate_exponent * std::pow(eqv_stress/_strength[qp],_flow_rate_exponent-1.0)/_strength[qp];
 //  stz
     //   Real dflowrate_dseqv = 1/_t_tau*std::exp(-1/_chiv[qp])*std::exp(eqv_stress/_grain_pressure/_chiv[qp])* \
     //	    (1/(_grain_pressure*_chiv[qp])*(1-_strength[qp]/eqv_stress)+_strength[qp]/(eqv_stress*eqv_stress))*(eqv_stress>_strength[qp]);
//
       Real dflowrate_dseqv;
      if (eqv_stress>_strength[qp])
	  //   dflowrate_dseqv = 1/_t_tau*std::exp(-1/_chiv[qp])*std::exp(eqv_stress/_grain_pressure/_chiv[qp])* \
	  (1.0/(_grain_pressure*_chiv[qp])*(1.0-_strength[qp]/eqv_stress)+_strength[qp]/(eqv_stress*eqv_stress));
    dflowrate_dseqv = 0.0;
	//  dflowrate_dseqv = 1/_t_tau*std::exp(-1/_chi)*std::exp(eqv_stress/_grain_pressure/_chi)* \
	      //  (1.0/(_grain_pressure*_chi)*(1.0-_strength[qp]/eqv_stress)+_strength[qp]/(eqv_stress*eqv_stress));

     RankTwoTensor tau = pk2_dev * _ce[qp];
     RankTwoTensor dseqv_dpk2dev;

     dseqv_dpk2dev.zero();
     if (eqv_stress > 0.0)
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
 stzVOLFlowRate::computePK2Deviatoric(const RankTwoTensor & pk2, const RankTwoTensor & ce) const
{
  return pk2 - (pk2.doubleContraction(ce) * ce.inverse())/3.0;
}
RankTwoTensor
 stzVOLFlowRate::computePK2Hydrostatic(const RankTwoTensor & pk2, const RankTwoTensor & ce) const
{
  return (pk2.doubleContraction(ce) * ce.inverse())/3.0*ce;
}

Real
stzVOLFlowRate::computeEqvStress(const RankTwoTensor & pk2_dev, const RankTwoTensor & ce) const
{
  RankTwoTensor sdev = pk2_dev * ce;
  Real val = sdev.doubleContraction(sdev.transpose());
  return std::pow(1.5 * val, 0.5);
}
