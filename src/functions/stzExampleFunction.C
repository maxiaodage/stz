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

#include "stzExampleFunction.h"

template<>
InputParameters validParams<stzExampleFunction>()
{
  InputParameters params = validParams<Function>();
  params.addParam<Real>("a0", 1.0, "Initial amplitude");
  params.addParam<Real>("a1", 1.0, "Finial amplitude");
  params.addParam<Real>("t0", 1.0, "Initial time");
  params.addParam<Real>("t1", 1.0, "Final time");

  return params;
}

stzExampleFunction::stzExampleFunction(const InputParameters & parameters) :
    Function(parameters),
    _a0(getParam<Real>("a0")),
    _a1(getParam<Real>("a1")),
    _t0(getParam<Real>("t0")),
    _t1(getParam<Real>("t1"))


{}

Real
stzExampleFunction::value(Real t, const Point & p)
{
  Real ksi= (t-_t0)/(_t1-_t0);

  if (t>=_t0&&t<=_t1){
  return _a0+(_a1-_a0)*pow(ksi,3)*(10-15*ksi+6*pow(ksi,2));
  }
  else{
    return _a1;
  }
}
