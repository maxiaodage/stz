/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "stzMaskedBodyForceCouple.h"
#include "Function.h"

template<>
InputParameters validParams<stzMaskedBodyForceCouple>()
{
  InputParameters params = validParams<BodyForce>();
  params.addClassDescription("Kernel that defines a body force modified by a material mask");
  params.addParam<Real>("chihat",0.12,"Steady State Chi value");
  params.addParam<Real>("pressure",25e6,"Grain Pressure");
  params.addParam<Real>("ko",0.1,"Grain Pressure");
  params.addParam<Real>("co",0.01,"Grain Pressure");
  params.addParam<Real>("Pth",21.11e6,"Grain Pressure");
//  params.addParam<Real>("chihat",0.12,"Grain Pressure");
  params.addParam<Real>("gammag",1.0,"Surface Energy");
  params.addParam<Real>("Eyoungs",7.0e9,"Youngs Modulus");
  // Biot
  params.addParam<Real>("biot",0.0,"Biot coefficient");
  params.addCoupledVar("grainsize",1e-4,"grain size") ;
  // Adding Coupling Pore pressure
  params.addCoupledVar("porepress",0.0,"Pore Pressure");
  return params;
}

stzMaskedBodyForceCouple::stzMaskedBodyForceCouple(const InputParameters & parameters) :
    BodyForce(parameters),
    _sdplcouple(getMaterialProperty<Real>("sdplcouple")),
  //  _chihat(getMaterialProperty<Real>("chihat")),
    _chihat(getParam<Real>("chihat")),
    _pressure(getParam<Real>("pressure")),
    _ko(getParam<Real>("ko")),
    _co(getParam<Real>("co")),
    _Pth(getParam<Real>("Pth")),
  //  _chihat(getParam<Real>("chihat")),
    _gammag(getParam<Real>("gammag")),
    _Eyoungs(getParam<Real>("Eyoungs")),
    // Biot
    _biot(getParam<Real>("biot")),
    // Adding coupled Grainsize
   _agrain(coupledValue("grainsize")),
   // Adding coupled porepress
   _pf(coupledValue("porepress"))
{
}

Real
stzMaskedBodyForceCouple::computeQpResidual()
{
  //  return BodyForce::computeQpResidual()*_mask[_qp];
  // return -_test[_i][_qp]*(1-_u[_qp]/_mask[_qp]);
  Real _Pthnew= std::sqrt(2.0*_Eyoungs*_gammag/(_agrain[_qp]*3.14));
//  std::cout<<"haha"<<_Pthnew<<"\n";
Real _pressure_eff=_pressure-_biot*_pf[_qp];

//  Real bigc=1.0/_pressure *((1-_ko*std::exp(-_Pth/_pressure))/_co);
//  Real bigc=1.0/_pressure_eff *((1.0-_ko*std::exp(-_Pthnew/_pressure_eff))/_co);

  Real bigc=1.0/(_pressure*_co);
  // return -_test[_i][_qp]*(1.0-_u[_qp]/_chihat[_qp])*_sdpl[_qp]*bigc;
  return -_test[_i][_qp]*(1.0-_u[_qp]/_chihat)*_sdplcouple[_qp]*bigc;

}
