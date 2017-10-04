/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZPFDIFFUSION_H
#define STZPFDIFFUSION_H

#include "Diffusion.h"
#include "Material.h"

//Forward Declarations
class stzpfdiffusion;

template<>
InputParameters validParams<stzpfdiffusion>();

/**
 * Note: This class is named HeatConductionKernel instead of HeatConduction
 * to avoid a clash with the HeatConduction namespace.  It is registered
 * as HeatConduction, which means it can be used by that name in the input
 * file.
 */
class stzpfdiffusion : public Diffusion
{
public:

  stzpfdiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  Real _alpha_hy;

};

#endif //STZCHIDIFFUSION_H
