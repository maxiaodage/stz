/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "stzagrainMaskedBodyForce.h"
#include "Function.h"

template<>
InputParameters validParams<stzagrainMaskedBodyForce>()
{
  InputParameters params = validParams<BodyForce>();
  params.addClassDescription("Kernel that defines a body force modified by a material mask");
  params.addParam<Real>("pressure",25e6,"Grain Pressure");
  params.addParam<Real>("ko",0.1,"Grain Pressur");
  params.addParam<Real>("co",0.01,"Grain Pressure");
  params.addParam<Real>("Pth",21.11e6,"Grain Pressure");
  params.addParam<Real>("chihat",0.12,"Steady state chi");
  params.addParam<Real>("gammag",1.0,"Surface Energy");
  params.addParam<Real>("Eyoungs",7.0e9,"Youngs Modulus");
  return params;
}

stzagrainMaskedBodyForce::stzagrainMaskedBodyForce(const InputParameters & parameters) :
    BodyForce(parameters),
    _sdpl(getMaterialProperty<Real>("sourcet")),
    _pressure(getParam<Real>("pressure")),
    _ko(getParam<Real>("ko")),
    _co(getParam<Real>("co")),
    _Pth(getParam<Real>("Pth")),
    _chihat(getParam<Real>("chihat")),
    _gammag(getParam<Real>("gammag")),
    _Eyoungs(getParam<Real>("Eyoungs"))
{
}

Real
stzagrainMaskedBodyForce::computeQpResidual()
{
  //  return BodyForce::computeQpResidual()*_mask[_qp];
  // return -_test[_i][_qp]*(1-_u[_qp]/_mask[_qp]);
  // Critile Pressure
  Real _Pthnew= std::sqrt(2.0*_Eyoungs*_gammag/(_u[_qp]*3.14));
//  Real _Pthnew= _Pth;
  //   std::cout<<"haha"<<_Pthnew<<"\n";

  Real bigc=_ko*std::exp(-_Pthnew/_pressure)*_sdpl[_qp]/_gammag;

   //Real bigc=1.0/_pressure;
   return _test[_i][_qp]*_u[_qp]*_u[_qp]*bigc;
}
