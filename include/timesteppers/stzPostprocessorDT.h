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

#ifndef STZPOSTPROCESSORDT_H
#define STZPOSTPROCESSORDT_H

#include "TimeStepper.h"
#include "PostprocessorInterface.h"

class stzPostprocessorDT;

template <>
InputParameters validParams<stzPostprocessorDT>();

/**
 * Computes the value of dt based on a postprocessor value
 */
class stzPostprocessorDT : public TimeStepper, public PostprocessorInterface
{
public:
  stzPostprocessorDT(const InputParameters & parameters);

protected:
  virtual Real computeInitialDT() override;
  virtual Real computeDT() override;

  const PostprocessorValue & _pps_s;
  const PostprocessorValue & _pps_s0;
  const PostprocessorValue & _pps_rpl;

  bool _has_initial_dt;
  Real _initial_dt;

  /// Multiplier applied to the postprocessor value
  const Real & _scale;

  /// Factor added to the postprocessor value
  const Real & _factor;
};

#endif /* POSTPROCESSORDT_H */
