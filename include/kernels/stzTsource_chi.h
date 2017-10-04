/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZTSOURCES_CHI_H
#define STZTSOURCES_CHI_H

#include "BodyForce.h"

//Forward Declarations
class stzTsource_chi;
class Function;

template<>
InputParameters validParams<stzTsource_chi>();

/**
 * This kernel creates a body force that is modified by a mask defined
 * as a material. Common uses of this would be to turn off or change the
 * body force in certain regions of the mesh.
 */

class stzTsource_chi : public BodyForce
{
public:

  stzTsource_chi(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();


  // Adding the Eq[3] constants
  Real _rhoc;
  Real _chihat_const;
  Real _r_pl_thre;
  Real _q0;
  Real _p;
  Real _a;
  Real _density;
  Real _alpha_biot;
  Real _chiswitch;
  const VariableValue &_s;
  const VariableValue &_couple_chi;
  const VariableValue &_couple_rpl;
  const VariableValue &_couple_pf;



};

#endif
