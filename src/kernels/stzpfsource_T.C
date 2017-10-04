/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "stzpfsource_T.h"
#include "Function.h"

template<>
InputParameters validParams<stzpfsource_T>()
{
  InputParameters params = validParams<BodyForce>();
  params.addClassDescription("Kernel that defines T dot term in Pore Pressure");
  params.addParam<Real>("Lambda",1.0,"Lambda parameter in front of T_dot in pore pressure");
  params.addCoupledVar("couple_T",100,"Coupling Temperature for using in the T_dot in pore pressure") ;



  return params;
}

stzpfsource_T::stzpfsource_T(const InputParameters & parameters) :
    BodyForce(parameters),
    _Lambda(getParam<Real>("Lambda")),
    _Tdot(coupledDot("couple_T"))


{
}

Real
stzpfsource_T::computeQpResidual()
{




  // return -_test[_i][_qp]*(1.0-_u[_qp]/_chihat)*bigc;
  return -_test[_i][_qp]*(_Lambda*_Tdot[_qp]);












}
