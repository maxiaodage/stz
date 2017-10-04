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

#include "stzplasticrateAux_bwdot.h"

template <>
InputParameters
validParams<stzplasticrateAux_bwdot>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addCoupledVar("coupled_cur",0.0, "Coupled variable current");
  params.addCoupledVar("coupled_old",0.0, "Coupled variable old");


  return params;
}

stzplasticrateAux_bwdot::stzplasticrateAux_bwdot(const InputParameters & parameters)
  : AuxKernel(parameters),
    _coupled_cur(coupledValue("coupled_cur")),
    _coupled_old(coupledValueOld("coupled_old"))

{
}

/**
 * Auxiliary Kernels override computeValue() instead of computeQpResidual().  Aux Variables
 * are calculated either one per elemenet or one per node depending on whether we declare
 * them as "Elemental (Constant Monomial)" or "Nodal (First Lagrange)".  No changes to the
 * source are necessary to switch from one type or the other.
 */
Real
stzplasticrateAux_bwdot::computeValue()
{
  Real var_dot = (_coupled_cur[_qp]-_coupled_old[_qp])/_dt;


//  return _coupled_s[_qp] + _value;
  return var_dot;
}
