/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "stzTdiffusion.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<stzTdiffusion>()
{
  InputParameters params = validParams<Diffusion>();
  params.addClassDescription("Compute thermal conductivity");
    params.addParam<Real>("alpha_th",0.0,"Tmperature diffusivity alpha_th in front of the Temperature diffusion term");
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

stzTdiffusion::stzTdiffusion(const InputParameters & parameters) :
    Diffusion(parameters),
    _alpha_th(getParam<Real>("alpha_th"))
{
}

Real
stzTdiffusion::computeQpResidual()
{
//  return _agrain[_qp]*_agrain[_qp]*_diffusivity[_qp]*Diffusion::computeQpResidual();
//  return _D0*_a*_a*_couple_rpl[_qp]*Diffusion::computeQpResidual();

  return _alpha_th*Diffusion::computeQpResidual();

}

Real
stzTdiffusion::computeQpJacobian()
{

// Real jac = _agrain[_qp]*_agrain[_qp]*_diffusivity[_qp] * Diffusion::computeQpJacobian();
 //Real jac = _D0*_a*_a*_couple_rpl[_qp]*Diffusion::computeQpJacobian();
 Real jac = _alpha_th*Diffusion::computeQpJacobian();


  return jac;
  // if (_diffusion_coefficient_dT)
  // jac += (*_diffusion_coefficient_dT)[_qp] * _phi[_j][_qp] * Diffusion::computeQpResidual();
  // return jac;
}
