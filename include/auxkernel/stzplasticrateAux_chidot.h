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

#ifndef STZPLASTICRATEAUX_CHIDOT_H
#define STZPLASTICRATEAUX_CHIDOT_H

#include "AuxKernel.h"

// Forward Declarations
class stzplasticrateAux_chidot;

template <>
InputParameters validParams<stzplasticrateAux_chidot>();

/**
 * Coupled auxiliary value
 */
class stzplasticrateAux_chidot : public AuxKernel
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  stzplasticrateAux_chidot(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
  const VariableValue & _coupled_chidot;



};

#endif // STZPLASTICRATEAUX_H
