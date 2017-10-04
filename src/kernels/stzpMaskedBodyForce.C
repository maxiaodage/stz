/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "stzpMaskedBodyForce.h"
#include "Function.h"

template<>
InputParameters validParams<stzpMaskedBodyForce>()
{
  InputParameters params = validParams<BodyForce>();
  params.addClassDescription("Kernel that defines a body force modified by a material mask");
  params.addParam<Real>("co",0.01,"dilatancy ratio");
  params.addCoupledVar("chicur",1e-4,"current compacitivty") ;
//  params.addCoupledVar("chiold",1e-4,"old compacitivty") ;

  return params;
}

stzpMaskedBodyForce::stzpMaskedBodyForce(const InputParameters & parameters) :
    BodyForce(parameters),
    _co(getParam<Real>("co")),
    // Adding coupled Grainsize
   _chidot(coupledDot("chicur"))
{
}

Real
stzpMaskedBodyForce::computeQpResidual()
{
  //Real _chidot=(_chicur[_qp]-_chiold[_qp])/_dt;
  //  return BodyForce::computeQpResidual()*_mask[_qp];
  // return -_test[_i][_qp]*(1-_u[_qp]/_mask[_qp]);
//  std::cout<<"haha"<<_Pthnew<<"\n";
   return _test[_i][_qp]*_co*_chidot[_qp];
}
