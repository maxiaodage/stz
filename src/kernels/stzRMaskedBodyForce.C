/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "stzRMaskedBodyForce.h"
#include "Function.h"

template<>
InputParameters validParams<stzRMaskedBodyForce>()
{
  InputParameters params = validParams<BodyForce>();
  params.addClassDescription("Kernel that defines a body force modified by a material mask");
  params.addParam<Real>("v",1.0e8,"Grain Pressure");
  params.addParam<Real>("co",0.01,"Grain Pressure");
  params.addParam<Real>("chihat",0.08,"Grain Pressure");

  return params;
}

stzRMaskedBodyForce::stzRMaskedBodyForce(const InputParameters & parameters) :
    BodyForce(parameters),
    _smises(getMaterialProperty<Real>("smises")),
    _v(getParam<Real>("v")),
    _co(getParam<Real>("co")),
    _chihat(getParam<Real>("chihat"))

{}

Real
stzRMaskedBodyForce::computeQpResidual()
{
  Real qfunc;
  if (_smises[_qp]>1.0)
  {
    qfunc = 1.0/_smises[_qp]*(_smises[_qp]-1.0)*(_smises[_qp]-1.0);
  }
  else
  {
    qfunc=0.0;
  }
   return -_test[_i][_qp]*(_v/_co)*std::exp(-1/_u[_qp])*(_chihat-_u[_qp])*_smises[_qp]*qfunc;
}
