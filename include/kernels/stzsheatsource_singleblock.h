/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZSHEATSOURCE_SINGLEBLOCK_H
#define STZSHEATSOURCE_SINGLEBLOCK_H

#include "BodyForce.h"

//Forward Declarations
class stzsheatsource_singleblock;
class Function;

template<>
InputParameters validParams<stzsheatsource_singleblock>();

/**
 * This kernel creates a body force that is modified by a mask defined
 * as a material. Common uses of this would be to turn off or change the
 * body force in certain regions of the mesh.
 */

class stzsheatsource_singleblock: public BodyForce
{
public:

  stzsheatsource_singleblock(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();



  // Adding the Eq[3] constants
  Real _dbgswitch;
  Real _shearG;
  Real _imprate;
  Real _e0;
  Real _ez;
  Real _e1;
  Real _p;
  Real _a;
  Real _density;
  Real _rho;
  Real _s1l;
  Real _s2l;
  Real _R0;
  Real _L;
  const VariableValue &_chi;
  const VariableValue &_rpl;

};

#endif
