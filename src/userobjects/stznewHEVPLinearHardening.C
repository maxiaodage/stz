/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "stznewHEVPLinearHardening.h"
#include "AuxKernel.h"
#include "RankTwoTensor.h"
#include "RankTwoScalarTools.h"
#include "RankTwoScalarAux.h"


template<>
InputParameters validParams<stznewHEVPLinearHardening>()
{
  InputParameters params = validParams<HEVPStrengthUOBase>();
  params.addParam<Real>("friction_angle",30.0, "Frictonal angle of friction");
  params.addParam<Real>("cohesion",1.0e6, "cohesion for yield stress");

  params.addClassDescription("User Object for linear hardening");
  return params;
}

stznewHEVPLinearHardening::stznewHEVPLinearHardening(const InputParameters & parameters) :
    HEVPStrengthUOBase(parameters),
    _phi(getParam<Real>("friction_angle")),
    _cohesion(getParam<Real>("cohesion")),
// Adding this tensor stress
    _tensor(getMaterialProperty<RankTwoTensor>("stress"))



{
}

bool
stznewHEVPLinearHardening::computeValue(unsigned int qp, Real & val) const
{
  //val = _sig0 +  _slope * _intvar[qp];
  Real invI1 = _tensor[qp].trace();
  Real s_xx = _tensor[qp](0,0);
  Real s_yy = _tensor[qp](1,1);
  Real tau_xy= _tensor[qp](0,1);
  Real s_mean = 1.0/2.0*(s_xx+s_yy);
  Real PI = 3.1415926;
  Real A_c;
  Real B_c;
  A_c = (6.0*_cohesion*std::cos(_phi*PI/180.0))/(std::pow(3,0.5)*(3.0-std::sin(_phi*PI/180.0)));
  B_c = (2.0*std::sin(_phi*PI/180.0))/(std::pow(3,0.5)*(3.0-std::sin(_phi*PI/180.0)));

  // Real s_yield =(_cohesion*std::cos(_phi*PI/180.0))-s_mean*std::sin(_phi*PI/180.0);
  Real s_yield = A_c-B_c*invI1;
  //So = (c + sigma_mean x tan(phi)) [cos (phi)]^2 ;


val= s_yield;



//std::cout<<"rationew="<<rationew<<"**"<<"val="<<val<<"**"<<"vol="<<volnew<<"iffnew="<<diffnew<<"**"<<"trueratio="<<volnew/diffnew<<"\n";

//    std::cout<<"haha"<<val<<"\n";



  return true;
}

bool
stznewHEVPLinearHardening::computeDerivative(unsigned int /*qp*/, const std::string & coupled_var_name, Real & val) const
{
  val = 0;

  // if (_intvar_prop_name == coupled_var_name)
  //   val = _slope;

  return true;
}
