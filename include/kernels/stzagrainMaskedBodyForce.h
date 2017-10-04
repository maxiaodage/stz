/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZAGRAINMASKEDBODYFORCE_H
#define STZAGRAINMASKEDBODYFORCE_H

#include "BodyForce.h"

//Forward Declarations
class stzagrainMaskedBodyForce;
class Function;

template<>
InputParameters validParams<stzagrainMaskedBodyForce>();

/**
 * This kernel creates a body force that is modified by a mask defined
 * as a material. Common uses of this would be to turn off or change the
 * body force in certain regions of the mesh.
 */

class stzagrainMaskedBodyForce : public BodyForce
{
public:

  stzagrainMaskedBodyForce(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  const MaterialProperty<Real> & _sdpl;
  // Adding the Eq[3] constants
  Real _pressure;
  Real _ko;
  Real _co;
  Real _Pth;
  Real _chihat;
  Real _gammag;
  Real _Eyoungs;

};

#endif
