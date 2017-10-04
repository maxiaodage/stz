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

#include "stzPostprocessorDT.h"

template <>
InputParameters
validParams<stzPostprocessorDT>()
{
  InputParameters params = validParams<TimeStepper>();
  params.addRequiredParam<PostprocessorName>("postprocessor_s","The name of the postprocessor that computes the dt");
  params.addRequiredParam<PostprocessorName>("postprocessor_s0","The name of the postprocessor that computes the dt");
  params.addRequiredParam<PostprocessorName>("postprocessor_rpl","The name of the postprocessor that computes the dt");
  params.addParam<Real>("dt", "Initial value of dt");
  params.addParam<Real>("scale", 1, "Multiple scale and supplied postprocessor value.");
  params.addParam<Real>("factor", 0, "Add a factor to the supplied postprocessor value.");
  return params;
}

stzPostprocessorDT::stzPostprocessorDT(const InputParameters & parameters)
  : TimeStepper(parameters),
    PostprocessorInterface(this),
    _pps_s(getPostprocessorValue("postprocessor_s")),
    _pps_s0(getPostprocessorValue("postprocessor_s0")),
    _pps_rpl(getPostprocessorValue("postprocessor_rpl")),
    _has_initial_dt(isParamValid("dt")),
    _initial_dt(_has_initial_dt ? getParam<Real>("dt") : 0.),
    _scale(getParam<Real>("scale")),
    _factor(getParam<Real>("factor"))
{
}

Real
stzPostprocessorDT::computeInitialDT()
{
  if (_has_initial_dt)
    return _initial_dt;
  else
    return computeDT();
}

Real
stzPostprocessorDT::computeDT()
{
  Real s = _pps_s;
  Real s0= _pps_s0;
  Real rpl = _pps_rpl;
  if ((s<s0)&&(rpl>0.0))
  {
    return getCurrentDT()*_scale;
  }
  else
  {
    return getCurrentDT();
  }
}
