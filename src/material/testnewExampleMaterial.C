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

#include "testnewExampleMaterial.h"

template<>
InputParameters validParams<testnewExampleMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addParam<Real>("chihatst", 0.08, "Initial time");
  params.addCoupledVar("chicur",1e-4,"current compacitivty") ;




  return params;
}

testnewExampleMaterial::testnewExampleMaterial(const InputParameters & parameters) :
    Material(parameters),
    // Declare that this material is going to provide a Real
    // valued property named "diffusivity" that Kernels can use.
    _sourcet(declareProperty<Real>("sourcet")),
     _diffusivity(declareProperty<Real>("diffusivity")),
     //adding rate depedent chihat
     _chihat(declareProperty<Real>("chihat")),
     _dpl(declareProperty<RankTwoTensor>("plrate")),
//###############
     _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _pk2_prop_name(_base_name + "pk2"),
    _pk2(getMaterialProperty<RankTwoTensor>(_pk2_prop_name)),
    _ce(getMaterialProperty<RankTwoTensor>(_base_name + "ce")),
    _fp(getMaterialProperty<RankTwoTensor>(_base_name + "fp")),
    _fp_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "fp")),
    _stensor(getMaterialProperty<RankTwoTensor>("stress")),
    _chihatst(getParam<Real>("chihatst"))

{}

void
testnewExampleMaterial::computeQpProperties()
{
  RankTwoTensor pk2_dev = computePK2Deviatoric(_pk2[_qp], _ce[_qp]);
  RankTwoTensor s_dev = computeSdev(pk2_dev, _ce[_qp]);
  RankTwoTensor s_devia = _stensor[_qp].deviatoric();
  Real          s_hydro =1.0/3.0*_stensor[_qp].trace();
  RankTwoTensor fp_dot=computefpdot(_fp[_qp],_fp_old[_qp]);
  RankTwoTensor lp=computeLp(fp_dot,_fp[_qp]);
  RankTwoTensor dpl=0.5*(lp+lp.transpose());
  Real sdpl=s_dev.doubleContraction(dpl.transpose());
  Real dpldb= dpl.doubleContraction(dpl.transpose());
  Real eqv_ps= std::pow(1.5 * dpldb, 0.5);
  _sourcet[_qp] =sdpl;
  _diffusivity[_qp]=eqv_ps;
  _dpl[_qp] = dpl;
  _chihat[_qp]=_chihatst;
    // std::cout <<"haha"<<sdpl/25e6<<"\n";
}

RankTwoTensor
testnewExampleMaterial::computeLp(const RankTwoTensor & fp_dot, const RankTwoTensor & fp) const
{
  return fp_dot*fp.inverse();
}

RankTwoTensor
testnewExampleMaterial::computefpdot(const RankTwoTensor & fp, const RankTwoTensor & fp_old) const
{
  return (fp-fp_old)/_dt;
}

RankTwoTensor
testnewExampleMaterial::computePK2Deviatoric(const RankTwoTensor & pk2, const RankTwoTensor & ce) const
{
  return pk2 - (pk2.doubleContraction(ce) * ce.inverse())/3.0;
}

RankTwoTensor
testnewExampleMaterial::computeSdev(const RankTwoTensor & pk2_dev, const RankTwoTensor & ce) const
{
  return pk2_dev * ce;
  // Real val = sdev.doubleContraction(sdev.transpose());
  // return std::pow(1.5 * val, 0.5);
}
