/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZPFSOURCES_T_H
#define STZPFSOURCES_T_H

#include "BodyForce.h"

//Forward Declarations
class stzpfsource_T;
class Function;

template<>
InputParameters validParams<stzpfsource_T>();

/**
 * This kernel creates a body force that is modified by a mask defined
 * as a material. Common uses of this would be to turn off or change the
 * body force in certain regions of the mesh.
 */

class stzpfsource_T : public BodyForce
{
public:

  stzpfsource_T(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();


  // Adding the Eq[3] constants
  Real _Lambda;
  const VariableValue &_Tdot;


};

#endif
