/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "stzHeatConduction.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<stzHeatConduction>()
{
  InputParameters params = validParams<Diffusion>();
  params.addClassDescription("Compute thermal conductivity");
    params.addParam<Real>("diffterm",0.0,"Diffusionterm");
    params.addCoupledVar("grainsize",1e-4,"grain size") ;

  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

stzHeatConduction::stzHeatConduction(const InputParameters & parameters) :
    Diffusion(parameters),
    _diffusivity(getMaterialProperty<Real>("diffusivity")),
    _diffterm(getParam<Real>("diffterm")),
    // Adding coupled Grainsize
   _agrain(coupledValue("grainsize"))
    //    _grainsize(getParam<Real>("grainsize"))
{
}

Real
stzHeatConduction::computeQpResidual()
{
//  return _agrain[_qp]*_agrain[_qp]*_diffusivity[_qp]*Diffusion::computeQpResidual();
  return _diffterm*Diffusion::computeQpResidual();
}

Real
stzHeatConduction::computeQpJacobian()
{
// Real jac = _agrain[_qp]*_agrain[_qp]*_diffusivity[_qp] * Diffusion::computeQpJacobian();
 Real jac =_diffterm * Diffusion::computeQpJacobian();
  return jac;
  // if (_diffusion_coefficient_dT)
  // jac += (*_diffusion_coefficient_dT)[_qp] * _phi[_j][_qp] * Diffusion::computeQpResidual();
  // return jac;
}
