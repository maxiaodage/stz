/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZHEATSOURCE_H
#define STZHEATSOURCE_H

#include "BodyForce.h"

//Forward Declarations
class stzheatsource;
class Function;

template<>
InputParameters validParams<stzheatsource>();

/**
 * This kernel creates a body force that is modified by a mask defined
 * as a material. Common uses of this would be to turn off or change the
 * body force in certain regions of the mesh.
 */

class stzheatsource: public BodyForce
{
public:

  stzheatsource(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  const MaterialProperty<Real> & _sdpl;
//  const MaterialProperty<Real> & _chihat;

  // Adding the Eq[3] constants
  Real _chihat;
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
