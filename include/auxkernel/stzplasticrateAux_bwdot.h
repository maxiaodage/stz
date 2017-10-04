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

#ifndef STZPLASTICRATEAUX_BWDOT_H
#define STZPLASTICRATEAUX_BWDOT_H

#include "AuxKernel.h"

// Forward Declarations
class stzplasticrateAux_bwdot;

template <>
InputParameters validParams<stzplasticrateAux_bwdot>();

/**
 * Coupled auxiliary value
 */
class stzplasticrateAux_bwdot : public AuxKernel
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  stzplasticrateAux_bwdot(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
  const VariableValue & _coupled_cur;
  const VariableValue & _coupled_old;


  Real _value;
};

#endif // STZPLASTICRATEAUX_H
