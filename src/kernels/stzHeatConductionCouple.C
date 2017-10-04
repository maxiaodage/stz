/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "stzHeatConductionCouple.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<stzHeatConductionCouple>()
{
  InputParameters params = validParams<Diffusion>();
  params.addClassDescription("Compute thermal conductivity");
    params.addParam<Real>("diffterm",0.0,"Diffusion term on or not 1.0 = on 0.0 = off");
    params.addParam<Real>("grainsize",1.0e-4,"Grain Size");
//    params.addCoupledVar("grainsize",1e-4,"grain size") ;

  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

stzHeatConductionCouple::stzHeatConductionCouple(const InputParameters & parameters) :
    Diffusion(parameters),
    //_diffusivity(getMaterialProperty<Real>("diffusivity")),
    _equiprate(getMaterialProperty<Real>("equiprate")),
    _diffterm(getParam<Real>("diffterm")),
    _grainsize(getParam<Real>("grainsize"))
    // Adding coupled Grainsize
  // _agrain(coupledValue("grainsize"))
{
}

Real
stzHeatConductionCouple::computeQpResidual()
{
//  return _agrain[_qp]*_agrain[_qp]*_diffusivity[_qp]*Diffusion::computeQpResidual();
//  return _diffterm*Diffusion::computeQpResidual();
 return _diffterm*_grainsize*_grainsize*_equiprate[_qp]*Diffusion::computeQpResidual();
 //return _diffterm*Diffusion::computeQpResidual();

}

Real
stzHeatConductionCouple::computeQpJacobian()
{
// Real jac = _agrain[_qp]*_agrain[_qp]*_diffusivity[_qp] * Diffusion::computeQpJacobian();
 //Real jac =_diffterm * Diffusion::computeQpJacobian();

 Real jac = _diffterm*_grainsize*_grainsize*_equiprate[_qp] * Diffusion::computeQpJacobian();
// Real jac = _diffterm* Diffusion::computeQpJacobian();

  return jac;
  // if (_diffusion_coefficient_dT)
  // jac += (*_diffusion_coefficient_dT)[_qp] * _phi[_j][_qp] * Diffusion::computeQpResidual();
  // return jac;
}
