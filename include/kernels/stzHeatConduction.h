/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZHEATCONDUCTION_H
#define STZHEATCONDUCTION_H

#include "Diffusion.h"
#include "Material.h"

//Forward Declarations
class stzHeatConduction;

template<>
InputParameters validParams<stzHeatConduction>();

/**
 * Note: This class is named HeatConductionKernel instead of HeatConduction
 * to avoid a clash with the HeatConduction namespace.  It is registered
 * as HeatConduction, which means it can be used by that name in the input
 * file.
 */
class stzHeatConduction : public Diffusion
{
public:

  stzHeatConduction(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  const MaterialProperty<Real> & _diffusivity;
  Real _diffterm;

  //Real _grainsize;
  const VariableValue &_agrain;


  //private:
  // const unsigned _dim;
  // const MaterialProperty<Real> & _diffusion_coefficient;
  // const MaterialProperty<Real> * const _diffusion_coefficient_dT;
};

#endif //STZHEATCONDUCTION_H
