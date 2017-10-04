/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZMASKEDBODYFORCE_H
#define STZMASKEDBODYFORCE_H

#include "BodyForce.h"

//Forward Declarations
class stzMaskedBodyForce;
class Function;

template<>
InputParameters validParams<stzMaskedBodyForce>();

/**
 * This kernel creates a body force that is modified by a mask defined
 * as a material. Common uses of this would be to turn off or change the
 * body force in certain regions of the mesh.
 */

class stzMaskedBodyForce : public BodyForce
{
public:

  stzMaskedBodyForce(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  const MaterialProperty<Real> & _sdpl;
  const MaterialProperty<Real> & _chihat;

  // Adding the Eq[3] constants
  Real _pressure;
  Real _ko;
  Real _co;
  Real _Pth;
  //Real _chihat;
  Real _gammag;
  Real _Eyoungs;
  Real _biot;
  const VariableValue &_agrain;
  const VariableValue &_pf;

};

#endif
