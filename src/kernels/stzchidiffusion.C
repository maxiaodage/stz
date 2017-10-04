/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "stzchidiffusion.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<stzchidiffusion>()
{
  InputParameters params = validParams<Diffusion>();
  params.addClassDescription("Compute thermal conductivity");
    params.addParam<Real>("D0",0.0,"D0 in the chi diffusion");
    params.addParam<Real>("a",350.0e-6,"a grain size");
    params.addCoupledVar("couple_rpl",0.0,"Couping the plastic strain rate");

  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

stzchidiffusion::stzchidiffusion(const InputParameters & parameters) :
    Diffusion(parameters),
    _D0(getParam<Real>("D0")),
    _a(getParam<Real>("a")),
    _couple_rpl(coupledValue("couple_rpl"))
{
}

Real
stzchidiffusion::computeQpResidual()
{
//  return _agrain[_qp]*_agrain[_qp]*_diffusivity[_qp]*Diffusion::computeQpResidual();
//  return _D0*_a*_a*_couple_rpl[_qp]*Diffusion::computeQpResidual();

Real diffusivity = 0.0 ;
if (_couple_rpl[_qp]>0.0)
{
  diffusivity = 1.0e-2;
}
  return diffusivity*Diffusion::computeQpResidual();

}

Real
stzchidiffusion::computeQpJacobian()
{
  Real diffusivity = 0.0 ;
  if (_couple_rpl[_qp]>0.0)
  {
    diffusivity = 1.0e-2;
  }
// Real jac = _agrain[_qp]*_agrain[_qp]*_diffusivity[_qp] * Diffusion::computeQpJacobian();
 //Real jac = _D0*_a*_a*_couple_rpl[_qp]*Diffusion::computeQpJacobian();
 Real jac = diffusivity*Diffusion::computeQpJacobian();


  return jac;
  // if (_diffusion_coefficient_dT)
  // jac += (*_diffusion_coefficient_dT)[_qp] * _phi[_j][_qp] * Diffusion::computeQpResidual();
  // return jac;
}
