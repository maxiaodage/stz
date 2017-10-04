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

#ifndef STZCOUPLEDSCALARAUX_H
#define STZCOUPLEDSCALARAUX_H

#include "AuxKernel.h"

// Forward Declarations
class stzCoupledScalarAux;

template <>
InputParameters validParams<stzCoupledScalarAux>();

/**
 * Coupled auxiliary scalar value
 */
class stzCoupledScalarAux : public AuxKernel
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  stzCoupledScalarAux(const InputParameters & parameters);

protected:
  virtual Real computeValue();

  /// Coupled variable
  VariableValue & _coupled_val;

  /// The component of the scalar variable
  unsigned int _component;
};

#endif // COUPLEDSCALARAUX_H
