/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef STZPLASTICRATEAUX_H
#define STZPLASTICRATEAUX_H

#include "AuxKernel.h"

// Forward Declarations
class stzplasticrateAux;

template <>
InputParameters validParams<stzplasticrateAux>();

/**
 * Coupled auxiliary value
 */
class stzplasticrateAux : public AuxKernel
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  stzplasticrateAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

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
  const VariableValue & _coupled_s;
  const VariableValue & _coupled_chi;


  Real _value;
};

#endif // STZPLASTICRATEAUX_H
