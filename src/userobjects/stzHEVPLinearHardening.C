/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "stzHEVPLinearHardening.h"
#include "AuxKernel.h"
#include "RankTwoTensor.h"

template<>
InputParameters validParams<stzHEVPLinearHardening>()
{
  InputParameters params = validParams<HEVPStrengthUOBase>();
  params.addParam<Real>("yield_stress",0.0, "Yield strength");
  params.addParam<Real>("slope",0.0, "Linear hardening slope");
  params.addParam<Real>("friction_angle",0.0, "Frictonal angle of friction");
  params.addParam<Real>("ratio_turnon",0.0,"The ratio term trun on coefficent ");
  params.addClassDescription("User Object for linear hardening");
  return params;
}

stzHEVPLinearHardening::stzHEVPLinearHardening(const InputParameters & parameters) :
    HEVPStrengthUOBase(parameters),
    _sig0(getParam<Real>("yield_stress")),
    _slope(getParam<Real>("slope")),
    _tanfric(getParam<Real>("friction_angle")),
// Adding this tensor stress
    _tensor(getMaterialProperty<RankTwoTensor>("stress")),
    _volplc(getMaterialProperty<Real>("volmetricpe")),
    _diffc(getMaterialProperty<Real>("diffusivity")),
    _beta(getParam<Real>("ratio_turnon")),
    _alphanew(getMaterialProperty<Real>("alpha_prime")),
    _sohydro(getMaterialProperty<Real>("sonew"))



{
}

bool
stzHEVPLinearHardening::computeValue(unsigned int qp, Real & val) const
{
  //val = _sig0 +  _slope * _intvar[qp];
  Real invI1 = _tensor[qp].trace();
  Real hydro = 1.0/3.0*invI1;
  RankTwoTensor s_dev= _tensor[qp].deviatoric();
  Real vmises = std::pow(3.0 / 2.0 * s_dev.doubleContraction(s_dev), 0.5);
  Real bb;
  Real sonew;


// MaterialProperty<Real> ap_copy = _alpha_prime[qp].set();
// ap_copy = _alpha_primeold[qp];
// if (sonew < _alpha_primeold[qp]*hydro)
//   _alpha_prime[qp].set(_alpha_primeold[qp]);
val= _sohydro[qp];



//std::cout<<"rationew="<<rationew<<"**"<<"val="<<val<<"**"<<"vol="<<volnew<<"iffnew="<<diffnew<<"**"<<"trueratio="<<volnew/diffnew<<"\n";

//    std::cout<<"haha"<<val<<"\n";



  return true;
}

bool
stzHEVPLinearHardening::computeDerivative(unsigned int /*qp*/, const std::string & coupled_var_name, Real & val) const
{
  val = 0;

  if (_intvar_prop_name == coupled_var_name)
    val = _slope;

  return true;
}
