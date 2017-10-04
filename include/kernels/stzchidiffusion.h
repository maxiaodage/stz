/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZCHIDIFFUSION_H
#define STZCHIDIFFUSION_H

#include "Diffusion.h"
#include "Material.h"

//Forward Declarations
class stzchidiffusion;

template<>
InputParameters validParams<stzchidiffusion>();

/**
 * Note: This class is named HeatConductionKernel instead of HeatConduction
 * to avoid a clash with the HeatConduction namespace.  It is registered
 * as HeatConduction, which means it can be used by that name in the input
 * file.
 */
class stzchidiffusion : public Diffusion
{
public:

  stzchidiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  Real _D0;
  Real _a;
  const VariableValue &_couple_rpl;

};

#endif //STZCHIDIFFUSION_H
