/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "RtestExampleMaterial.h"

template<>
InputParameters validParams<RtestExampleMaterial>()
{
  InputParameters params = validParams<Material>();
  return params;
}

RtestExampleMaterial::RtestExampleMaterial(const InputParameters & parameters) :
    Material(parameters),
    // Declare that this material is going to provide a Real
    // valued property named "diffusivity" that Kernels can use.
    _smises(declareProperty<Real>("smises")),
//###############
    _stensor(getMaterialProperty<RankTwoTensor>("stress"))

{}

void
RtestExampleMaterial::computeQpProperties()
{
  RankTwoTensor s_devia = _stensor[_qp].deviatoric();
  _smises[_qp]=std::pow(s_devia.doubleContraction(s_devia), 0.5);

}
