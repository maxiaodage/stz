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

#ifndef STZVELSTEPFUNCTION_H
#define STZVELSTEPFUNCTION_H

#include "Function.h"

class stzvelstepFunction;

template<>
InputParameters validParams<stzvelstepFunction>();

class stzvelstepFunction : public Function
{
public:
  stzvelstepFunction(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p);

protected:
  Real _a0;
  Real _a1;
  Real _av1;
  Real _av2;
  Real _t0;
  Real _t1;
  Real _tv1;
  Real _tv2;

};

#endif //STZVELSTEPFUNCTION_H
