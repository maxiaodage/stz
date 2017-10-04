/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "testExampleMaterial.h"

template<>
InputParameters validParams<testExampleMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addParam<Real>("aa", 0.56, "Initial amplitude");
  params.addParam<Real>("q0", 1.0, "Finial amplitude");
  params.addParam<Real>("tau", 8e-7, "Initial time");
  params.addParam<Real>("chihatst", 0.08, "Initial time");
  params.addParam<Real>("qthreshold", 1.0e20, "Initial time");
  params.addParam<Real>("co", 0.0, "Volumetric dilatcny");
  params.addParam<Real>("alphaconst", 0.0, "Volumetric dilatcny");
  params.addParam<Real>("cohesion", 5e6, "Cohesion");
  params.addParam<Real>("beta", 0.0, "Control Beta For dilatancy");
  params.addCoupledVar("chicur",1e-4,"current compacitivty") ;




  return params;
}

testExampleMaterial::testExampleMaterial(const InputParameters & parameters) :
    Material(parameters),
    // Declare that this material is going to provide a Real
    // valued property named "diffusivity" that Kernels can use.
    _sourcet(declareProperty<Real>("sourcet")),
     _diffusivity(declareProperty<Real>("diffusivity")),
     //adding rate depedent chihat
     _chihat(declareProperty<Real>("chihat")),
     _volpl(declareProperty<Real>("volmetricpe")),
     _ratio(declareProperty<Real>("ratio")),
     _alpha_prime(declareProperty<Real>("alpha_prime")),
     _alpha_primeold(declarePropertyOld<Real>("alpha_prime")),
     _sonew(declareProperty<Real>("sonew")),
     _smises(declareProperty<Real>("vommisesnew")),
     _nondeq(declareProperty<Real>("nondimensionps")),
     _dpl(declareProperty<RankTwoTensor>("plrate")),
//###############
     _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _pk2_prop_name(_base_name + "pk2"),
    _pk2(getMaterialProperty<RankTwoTensor>(_pk2_prop_name)),
    _ce(getMaterialProperty<RankTwoTensor>(_base_name + "ce")),
    _fp(getMaterialProperty<RankTwoTensor>(_base_name + "fp")),
    _fp_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "fp")),
    _stensor(getMaterialProperty<RankTwoTensor>("stress")),
    _aa(getParam<Real>("aa")),
    _q0(getParam<Real>("q0")),
    _tau(getParam<Real>("tau")),
    _chihatst(getParam<Real>("chihatst")),
    _qthre(getParam<Real>("qthreshold")),
    _co(getParam<Real>("co")),
    _alphaconst(getParam<Real>("alphaconst")),
    _cohesion(getParam<Real>("cohesion")),
    _beta(getParam<Real>("beta")),
    _chidot(coupledDot("chicur"))

{}

void
testExampleMaterial::computeQpProperties()
{
  RankTwoTensor pk2_dev = computePK2Deviatoric(_pk2[_qp], _ce[_qp]);
  RankTwoTensor s_dev = computeSdev(pk2_dev, _ce[_qp]);
  RankTwoTensor s_devia = _stensor[_qp].deviatoric();
  Real          s_hydro =1.0/3.0*_stensor[_qp].trace();
  // Real s_mises = std::pow(3.0 / 2.0 * s_devia.doubleContraction(s_devia), 0.5);
  _smises[_qp]=std::pow(3.0 / 2.0 * s_devia.doubleContraction(s_devia), 0.5);
  RankTwoTensor fp_dot=computefpdot(_fp[_qp],_fp_old[_qp]);
  RankTwoTensor lp=computeLp(fp_dot,_fp[_qp]);
  RankTwoTensor dpl=0.5*(lp+lp.transpose());
  Real sdpl=s_dev.doubleContraction(dpl.transpose());
  Real dpldb= dpl.doubleContraction(dpl.transpose());
  Real eqv_ps= std::pow(1.5 * dpldb, 0.5);
  _sourcet[_qp] =sdpl;
  _diffusivity[_qp]=eqv_ps;
  _volpl[_qp]=_co*_chidot[_qp];
  _dpl[_qp] = dpl;
  if (_diffusivity[_qp]>1.0e-6)
  {
    _ratio[_qp]=std::abs(_volpl[_qp]/_diffusivity[_qp]);
    _alpha_prime[_qp]=_alphaconst+_beta*_ratio[_qp];
  }
  else
  {
    _ratio[_qp] = 0.0;
    _alpha_prime[_qp]=_alphaconst;
  }
  _sonew[_qp]=_cohesion-(_alpha_primeold[_qp])*s_hydro;
   if ((_sonew[_qp])<_smises[_qp])
   {
    //  _alpha_prime[_qp]=_alpha_primeold[_qp];
     _sonew[_qp]=_cohesion-(_alpha_prime[_qp])*s_hydro;
       //std::cout<<"alpha="<<_alpha_prime[_qp]<<"\n";
   }
   else
   {
     _alpha_prime[_qp] = _alpha_primeold[_qp];
   }

   // chi_hat depdent on chi hand
   // Non dimensional
   _nondeq[_qp]=eqv_ps*_tau;
  if (_nondeq[_qp]<_qthre)
  {
    _chihat[_qp]=_chihatst;
  }
  else
  {
    _chihat[_qp]=_aa/log10(_q0/(_nondeq[_qp]));
  }
    // std::cout <<"haha"<<sdpl/25e6<<"\n";
}
void
testExampleMaterial::initQpStatefulProperties()
{
  // init the diffusivity property (this will become
  // _diffusivity_old in the first call of computeProperties)
  _alpha_prime[_qp] = _alphaconst;
}

RankTwoTensor
testExampleMaterial::computeLp(const RankTwoTensor & fp_dot, const RankTwoTensor & fp) const
{
  return fp_dot*fp.inverse();
}

RankTwoTensor
testExampleMaterial::computefpdot(const RankTwoTensor & fp, const RankTwoTensor & fp_old) const
{
  return (fp-fp_old)/_dt;
}

RankTwoTensor
testExampleMaterial::computePK2Deviatoric(const RankTwoTensor & pk2, const RankTwoTensor & ce) const
{
  return pk2 - (pk2.doubleContraction(ce) * ce.inverse())/3.0;
}

RankTwoTensor
testExampleMaterial::computeSdev(const RankTwoTensor & pk2_dev, const RankTwoTensor & ce) const
{
  return pk2_dev * ce;
  // Real val = sdev.doubleContraction(sdev.transpose());
  // return std::pow(1.5 * val, 0.5);
}
