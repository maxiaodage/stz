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

#include "stzshearExampleFunction.h"

template<>
InputParameters validParams<stzshearExampleFunction>()
{
  InputParameters params = validParams<Function>();
  params.addParam<Real>("a0", 1.0, "Initial amplitude");
  params.addParam<Real>("a1", 1.0, "Finial amplitude");
  params.addParam<Real>("t0", 1.0, "Initial time");
  params.addParam<Real>("t1", 1.0, "Final time");
  params.addParam<Real>("t2", 1.0, "Initial time");
  params.addParam<Real>("t3", 1.0, "Final time");


  return params;
}

stzshearExampleFunction::stzshearExampleFunction(const InputParameters & parameters) :
    Function(parameters),
    _a0(getParam<Real>("a0")),
    _a1(getParam<Real>("a1")),
    _t0(getParam<Real>("t0")),
    _t1(getParam<Real>("t1")),
    _t2(getParam<Real>("t2")),
    _t3(getParam<Real>("t3"))




{}

Real
stzshearExampleFunction::value(Real t, const Point & p)
{
  // Real ksi= (t-_t0)/(_t1-_t0);
  // Real ksi2 = (-t-(_t2+_t1)-_t0)/(_t1-_t0);
  // if (t<=_t0)
  // {
  //     return 0.0;
  // }
  // else if (t>_t0&&t<=_t1)
  // {
  // return (_a0+(_a1-_a0)*pow(ksi,3)*(10-15*ksi+6*pow(ksi,2)));
  // }
  // else if (t>_t1&&t<=_t2)
  // {
  //   return _a1;
  // }
  // else if (t>_t2&&t<=_t3)
  // {
  //   return (_a0+(_a1-_a0)*pow(ksi2,3)*(10-15*ksi2+6*pow(ksi2,2)));
  // }
  // else
  // {
  //   return _a0;
  // }
  if (t<=_t0)
  {
    return _a0;
  }
  else if (t>_t0&&t<_t1)
  {
    return (_a1-_a0)/(_t1-_t0)*(t-_t1)+_a1;
  }
  else if (t>=_t1&&t<=_t2)
  {
    return  _a1;
  }
  else if (t>_t2&&t<_t3)
  {
    return (_a0-_a1)/(_t3-_t2)*(t-_t2)+_a1;
  }
  else
  {
    return _a0;
  }
}
