/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "anaHEVPEqvPlasticStrainRate.h"

template<>
InputParameters validParams<anaHEVPEqvPlasticStrainRate>()
{
  InputParameters params = validParams<HEVPInternalVarRateUOBase>();
  params.addParam<Real>("h_scaling", 1.0, "Scaling parameter");
  params.addClassDescription("User Object computing equivalent plastic strain rate");
  return params;
}

anaHEVPEqvPlasticStrainRate::anaHEVPEqvPlasticStrainRate(const InputParameters & parameters) :
    HEVPInternalVarRateUOBase(parameters),
    _h(getParam<Real>("h_scaling"))
{
}

bool
anaHEVPEqvPlasticStrainRate::computeValue(unsigned int qp, Real & val) const
{
  val = _h * _flow_rate[qp];
  return true;
}

bool
anaHEVPEqvPlasticStrainRate::computeDerivative(unsigned int /*qp*/, const std::string & coupled_var_name, Real & val) const
{
  val = 0;

  if (_flow_rate_prop_name == coupled_var_name)
    val = _h;

  return true;
}
