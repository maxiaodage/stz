/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZPFSOURCES_CHI_H
#define STZPFSOURCES_CHI_H

#include "BodyForce.h"

//Forward Declarations
class stzpfsource_chi;
class Function;

template<>
InputParameters validParams<stzpfsource_chi>();

/**
 * This kernel creates a body force that is modified by a mask defined
 * as a material. Common uses of this would be to turn off or change the
 * body force in certain regions of the mesh.
 */

class stzpfsource_chi : public BodyForce
{
public:

  stzpfsource_chi(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();


  // Adding the Eq[3] constants
  Real _alpha;
  Real _beta;
  const VariableValue &_chidot;



};

#endif
