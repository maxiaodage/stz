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

#include "stzplasticrateAux_chidot.h"

template <>
InputParameters
validParams<stzplasticrateAux_chidot>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addCoupledVar("couple_chidot",0.0,"Coupling the pore Pressure pf");


  return params;
}

stzplasticrateAux_chidot::stzplasticrateAux_chidot(const InputParameters & parameters)
  : AuxKernel(parameters),
    _coupled_chidot(coupledDot("couple_chidot"))

{
}

/**
 * Auxiliary Kernels override computeValue() instead of computeQpResidual().  Aux Variables
 * are calculated either one per elemenet or one per node depending on whether we declare
 * them as "Elemental (Constant Monomial)" or "Nodal (First Lagrange)".  No changes to the
 * source are necessary to switch from one type or the other.
 */
Real
stzplasticrateAux_chidot::computeValue()
{

//  return _coupled_s[_qp] + _value;
  return _coupled_chidot[_qp];
}
