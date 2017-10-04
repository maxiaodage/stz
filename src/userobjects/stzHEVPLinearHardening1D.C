/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "stzHEVPLinearHardening1D.h"
#include "AuxKernel.h"
#include "RankTwoTensor.h"

template<>
InputParameters validParams<stzHEVPLinearHardening1D>()
{
  InputParameters params = validParams<HEVPStrengthUOBase>();
  params.addParam<Real>("yield_stress",0.0, "Yield strength");
  params.addParam<Real>("s1l",0.98,"s1l in the strength expression");
  params.addParam<Real>("s2l",-1.6,"s1l in the strength expression");
  params.addParam<Real>("p",80e6,"confining pressure");
  params.addCoupledVar("couple_chi",0.3,"effective temperature") ;




  params.addClassDescription("User Object for linear hardening");
  return params;
}

stzHEVPLinearHardening1D::stzHEVPLinearHardening1D(const InputParameters & parameters) :
    HEVPStrengthUOBase(parameters),
    _sig0(getParam<Real>("yield_stress")),
    _s1l(getParam<Real>("s1l")),
    _s2l(getParam<Real>("s2l")),
    _p(getParam<Real>("p")),
    _chiv(coupledValue("couple_chi"))




{
}

bool
stzHEVPLinearHardening1D::computeValue(unsigned int qp, Real & val) const
{
  //val = _sig0 +  _slope * _intvar[qp];



// MaterialProperty<Real> ap_copy = _alpha_prime[qp].set();
// ap_copy = _alpha_primeold[qp];
// if (sonew < _alpha_primeold[qp]*hydro)
//   _alpha_prime[qp].set(_alpha_primeold[qp]);
//val= _sohydro[qp];
val = _s1l*_p+_s2l*_p*_chiv[qp];



//std::cout<<"rationew="<<rationew<<"**"<<"val="<<val<<"**"<<"vol="<<volnew<<"iffnew="<<diffnew<<"**"<<"trueratio="<<volnew/diffnew<<"\n";

//    std::cout<<"haha"<<val<<"\n";



  return true;
}

bool
stzHEVPLinearHardening1D::computeDerivative(unsigned int /*qp*/, const std::string & coupled_var_name, Real & val) const
{
  val = 0;

  if (_intvar_prop_name == coupled_var_name)
  {
    val = 0.0;
  }
  //  val = _slope;

  return true;
}
