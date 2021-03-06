/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZCHISOURCE_H
#define STZCHISOURCE_H

#include "BodyForce.h"

//Forward Declarations
class stzchisource;
class Function;

template<>
InputParameters validParams<stzchisource>();

/**
 * This kernel creates a body force that is modified by a mask defined
 * as a material. Common uses of this would be to turn off or change the
 * body force in certain regions of the mesh.
 */

class stzchisource : public BodyForce
{
public:

  stzchisource(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();


  // Adding the Eq[3] constants
  Real _dbgswitch;
  Real _chihat_const;
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
  Real _eps0;
  Real _tfr;
  Real _r_pl_thre;
  Real _q0;
  const VariableValue &_s;
  const VariableValue &_couple_rpl;


};

#endif
