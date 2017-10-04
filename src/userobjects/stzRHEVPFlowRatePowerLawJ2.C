/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "stzRHEVPFlowRatePowerLawJ2.h"

template<>
InputParameters validParams<stzRHEVPFlowRatePowerLawJ2>()
{
  InputParameters params = validParams<HEVPFlowRateUOBase>();
  params.addParam<Real>("reference_flow_rate",0.001,"Reference flow rate for rate dependent flow");
  params.addParam<Real>("flow_rate_exponent",10.0,"Power law exponent in flow rate equation");
  params.addParam<Real>("flow_rate_tol", 1e5, "Tolerance for flow rate");
  params.addParam<Real>("tol_plastic", 1e2, "Tolerance for Plasticity to kick in");

  // v
  params.addParam<Real>("v",1.0e8,"coefficient v");


   // Adding Coupling Chiv
  params.addCoupledVar("efftemp",0.000000001,"effective temperature") ;

  //
   params.addClassDescription("User object to evaluate power law flow rate and flow direction based on J2");

   return params;
 }

 stzRHEVPFlowRatePowerLawJ2::stzRHEVPFlowRatePowerLawJ2(const InputParameters & parameters) :
     HEVPFlowRateUOBase(parameters),
     _ref_flow_rate(getParam<Real>("reference_flow_rate")),
     _flow_rate_exponent(getParam<Real>("flow_rate_exponent")),
     _flow_rate_tol(getParam<Real>("flow_rate_tol")),
     _v(getParam<Real>("v")),
     _chiv(coupledValue("efftemp"))
 {
 }

 bool
 stzRHEVPFlowRatePowerLawJ2::computeValue(unsigned int qp, Real & val) const
 {
   RankTwoTensor pk2_dev = computePK2Deviatoric(_pk2[qp], _ce[qp]);
   Real eqv_stress = computeEqvStress(pk2_dev, _ce[qp]);
   Real qfunc;
   if (eqv_stress>1.0)
   {
     qfunc = 1.0/eqv_stress*(eqv_stress-1.0)*(eqv_stress-1.0);
   }
   else
   {
     qfunc=0.0;
   }

     val = _v*std::exp(-1.0/_chiv[qp])*qfunc;

  if (val > _flow_rate_tol)
   {
#ifdef DEBUG
     mooseWarning("Flow rate greater than " , _flow_rate_tol ," " , val , " " ,eqv_stress ," " , _strength[qp]);
#endif
     return false;
   }
   return true;
 }

 bool
 stzRHEVPFlowRatePowerLawJ2::computeDirection(unsigned int qp, RankTwoTensor & val) const
 {
   RankTwoTensor pk2_dev = computePK2Deviatoric(_pk2[qp], _ce[qp]);
   Real eqv_stress = computeEqvStress(pk2_dev, _ce[qp]);

   val.zero();
   if (eqv_stress >1.0)
     {//   val = 1.5/eqv_stress * _ce[qp] * pk2_dev ;
        val = 1.0/eqv_stress*pk2_dev*_ce[qp];
      //
      //
     }

   return true;
 }

 bool
  stzRHEVPFlowRatePowerLawJ2::computeDerivative(unsigned int qp, const std::string & coupled_var_name, Real & val) const
 {
   val = 0.0;
   RankTwoTensor pk2_dev = computePK2Deviatoric(_pk2[qp], _ce[qp]);
   Real eqv_stress = computeEqvStress(pk2_dev, _ce[qp]);
   if (_strength_prop_name == coupled_var_name)
   {


     if(eqv_stress>=1.0)
        {
          val=0.0;
        }


   }

   return true;
 }

 bool
 stzRHEVPFlowRatePowerLawJ2::computeTensorDerivative(unsigned int qp, const std::string & coupled_var_name, RankTwoTensor & val) const
 {
   val.zero();

   if (_pk2_prop_name == coupled_var_name)
   {
     RankTwoTensor pk2_dev = computePK2Deviatoric(_pk2[qp], _ce[qp]);
     Real eqv_stress = computeEqvStress(pk2_dev, _ce[qp]);
       Real dflowrate_dseqv;
      if (eqv_stress>_strength[qp])
	  //   dflowrate_dseqv = 1/_t_tau*std::exp(-1/_chiv[qp])*std::exp(eqv_stress/_grain_pressure/_chiv[qp])* \
	  (1.0/(_grain_pressure*_chiv[qp])*(1.0-_strength[qp]/eqv_stress)+_strength[qp]/(eqv_stress*eqv_stress));
    dflowrate_dseqv = _v*std::exp(-1/_chiv[qp])*(1.0-1.0/(eqv_stress*eqv_stress));


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
 stzRHEVPFlowRatePowerLawJ2::computePK2Deviatoric(const RankTwoTensor & pk2, const RankTwoTensor & ce) const
{
  return pk2 - (pk2.doubleContraction(ce) * ce.inverse())/3.0;
}

Real
stzRHEVPFlowRatePowerLawJ2::computeEqvStress(const RankTwoTensor & pk2_dev, const RankTwoTensor & ce) const
{
  RankTwoTensor sdev = pk2_dev * ce;
  Real val = sdev.doubleContraction(sdev.transpose());
  return std::pow(1.0 * val, 0.5);
}
