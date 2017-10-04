/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZRMASKEDBODYFORCE_H
#define STZRMASKEDBODYFORCE_H

#include "BodyForce.h"

//Forward Declarations
class stzRMaskedBodyForce;
class Function;

template<>
InputParameters validParams<stzRMaskedBodyForce>();

/**
 * This kernel creates a body force that is modified by a mask defined
 * as a material. Common uses of this would be to turn off or change the
 * body force in certain regions of the mesh.
 */

class stzRMaskedBodyForce : public BodyForce
{
public:

  stzRMaskedBodyForce(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  const MaterialProperty<Real> & _smises;
  Real _v;
  Real _co;
  Real _chihat;
};

#endif
