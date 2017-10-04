/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "stzpfdiffusion.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<stzpfdiffusion>()
{
  InputParameters params = validParams<Diffusion>();
  params.addClassDescription("Compute thermal conductivity");
    params.addParam<Real>("alpha_hy",0.0,"Pf diffusivity alpha_hy in front of the Pf diffusion term");
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

stzpfdiffusion::stzpfdiffusion(const InputParameters & parameters) :
    Diffusion(parameters),
    _alpha_hy(getParam<Real>("alpha_hy"))
{
}

Real
stzpfdiffusion::computeQpResidual()
{
//  return _agrain[_qp]*_agrain[_qp]*_diffusivity[_qp]*Diffusion::computeQpResidual();
//  return _D0*_a*_a*_couple_rpl[_qp]*Diffusion::computeQpResidual();

  return _alpha_hy*Diffusion::computeQpResidual();

}

Real
stzpfdiffusion::computeQpJacobian()
{

// Real jac = _agrain[_qp]*_agrain[_qp]*_diffusivity[_qp] * Diffusion::computeQpJacobian();
 //Real jac = _D0*_a*_a*_couple_rpl[_qp]*Diffusion::computeQpJacobian();
 Real jac = _alpha_hy*Diffusion::computeQpJacobian();


  return jac;
  // if (_diffusion_coefficient_dT)
  // jac += (*_diffusion_coefficient_dT)[_qp] * _phi[_j][_qp] * Diffusion::computeQpResidual();
  // return jac;
}
