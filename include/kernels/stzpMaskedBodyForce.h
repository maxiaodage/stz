/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZPMASKEDBODYFORCE_H
#define STZPMASKEDBODYFORCE_H

#include "BodyForce.h"

//Forward Declarations
class stzpMaskedBodyForce;
class Function;

template<>
InputParameters validParams<stzpMaskedBodyForce>();

/**
 * This kernel creates a body force that is modified by a mask defined
 * as a material. Common uses of this would be to turn off or change the
 * body force in certain regions of the mesh.
 */

class stzpMaskedBodyForce : public BodyForce
{
public:

  stzpMaskedBodyForce(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  Real _co;
  const VariableValue &_chidot;
};

#endif
