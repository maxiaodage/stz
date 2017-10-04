/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "anaHEVPEqvPlasticStrain.h"

template<>
InputParameters validParams<anaHEVPEqvPlasticStrain>()
{
  InputParameters params = validParams<HEVPInternalVarUOBase>();
  params.addClassDescription("User Object to integrate equivalent plastic strain");
  return params;
}

anaHEVPEqvPlasticStrain::anaHEVPEqvPlasticStrain(const InputParameters & parameters) :
    HEVPInternalVarUOBase(parameters)
{
}

bool
anaHEVPEqvPlasticStrain::computeValue(unsigned int qp, Real dt, Real & val) const
{
  val = _this_old[qp] + _intvar_rate[qp] * dt;
  return true;
}

bool
anaHEVPEqvPlasticStrain::computeDerivative(unsigned int /*qp*/, Real dt, const std::string & coupled_var_name, Real & val) const
{
  val = 0;

  if (_intvar_rate_prop_name == coupled_var_name)
    val = dt;

  return true;
}
