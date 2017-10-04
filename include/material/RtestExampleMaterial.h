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

#ifndef RTESTEXAMPLEMATERIAL_H
#define RTESTEXAMPLEMATERIAL_H

#include "Material.h"
#include "LinearInterpolation.h"
#include "RankTwoTensor.h"

//Forward Declarations
class RtestExampleMaterial;

template<>
InputParameters validParams<RtestExampleMaterial>();

/**
 * Example material class that defines a few properties.
 */
class RtestExampleMaterial : public Material
{
public:
  RtestExampleMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  
  MaterialProperty<Real> & _smises;
  const MaterialProperty<RankTwoTensor> & _stensor;
};

#endif //TESTEXAMPLEMATERIAL_H
