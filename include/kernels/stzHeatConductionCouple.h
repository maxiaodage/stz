/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZHEATCONDUCTIONCOUPLE_H
#define STZHEATCONDUCTIONCOUPLE_H

#include "Diffusion.h"
#include "Material.h"

//Forward Declarations
class stzHeatConductionCouple;

template<>
InputParameters validParams<stzHeatConductionCouple>();

/**
 * Note: This class is named HeatConductionKernel instead of HeatConduction
 * to avoid a clash with the HeatConduction namespace.  It is registered
 * as HeatConduction, which means it can be used by that name in the input
 * file.
 */
class stzHeatConductionCouple : public Diffusion
{
public:

  stzHeatConductionCouple(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  //const MaterialProperty<Real> & _diffusivity;
  const MaterialProperty<Real> & _equiprate;

  Real _diffterm;

  Real _grainsize;
  //const VariableValue &_agrain;


  //private:
  // const unsigned _dim;
  // const MaterialProperty<Real> & _diffusion_coefficient;
  // const MaterialProperty<Real> * const _diffusion_coefficient_dT;
};

#endif //STZHEATCONDUCTIONCOUPLE_H
