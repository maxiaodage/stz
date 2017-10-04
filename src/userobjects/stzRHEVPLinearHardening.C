/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "stzRHEVPLinearHardening.h"
#include "AuxKernel.h"
#include "RankTwoTensor.h"

template<>
InputParameters validParams<stzRHEVPLinearHardening>()
{
  InputParameters params = validParams<HEVPStrengthUOBase>();
  params.addParam<Real>("yield_stress",0.0, "Yield strength");
  params.addParam<Real>("slope",0.0, "Linear hardening slope");
  params.addClassDescription("User Object for linear hardening");
  return params;
}

stzRHEVPLinearHardening::stzRHEVPLinearHardening(const InputParameters & parameters) :
    HEVPStrengthUOBase(parameters),
    _sig0(getParam<Real>("yield_stress")),
    _slope(getParam<Real>("slope"))

{
}

bool
stzRHEVPLinearHardening::computeValue(unsigned int qp, Real & val) const
{
val= _sig0;
  return true;
}

bool
stzRHEVPLinearHardening::computeDerivative(unsigned int /*qp*/, const std::string & coupled_var_name, Real & val) const
{
  val = 0;

  if (_intvar_prop_name == coupled_var_name)
    val = _slope;

  return true;
}
