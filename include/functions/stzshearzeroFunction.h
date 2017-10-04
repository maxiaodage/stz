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

#ifndef STZSHEARZEROFUNCTION_H
#define STZSHEARZEROFUNCTION_H

#include "Function.h"

class stzshearzeroFunction;

template<>
InputParameters validParams<stzshearzeroFunction>();

class stzshearzeroFunction : public Function
{
public:
  stzshearzeroFunction(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p);

protected:
  Real _a0;
  Real _a1;
  Real _t0;
  Real _t1;
};

#endif //STZSHEARZEROFUNCTION_H
