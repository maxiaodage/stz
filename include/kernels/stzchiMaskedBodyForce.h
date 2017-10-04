/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZCHIMASKEDBODYFORCE_H
#define STZCHIMASKEDBODYFORCE_H

#include "BodyForce.h"

//Forward Declarations
class stzchiMaskedBodyForce;
class Function;

template<>
InputParameters validParams<stzchiMaskedBodyForce>();

/**
 * This kernel creates a body force that is modified by a mask defined
 * as a material. Common uses of this would be to turn off or change the
 * body force in certain regions of the mesh.
 */

class stzchiMaskedBodyForce : public BodyForce
{
public:

  stzchiMaskedBodyForce(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  const MaterialProperty<Real> & _sdpl;
  const MaterialProperty<Real> & _diffusivity;
  // Adding the Eq[3] constants
  Real _pressure;
  Real _chihat;
  Real _chi_0;
  Real _tau;
  Real _epsilon_0;
  Real _epsilon_1;
  Real _miu_1;
  Real _miu_2;
  Real _tau_f;
  Real _ksi_0;

};

#endif
