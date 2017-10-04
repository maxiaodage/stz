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

#ifndef STZEXAMPLEFUNCTION_H
#define STZEXAMPLEFUNCTION_H

#include "Function.h"

class stzExampleFunction;

template<>
InputParameters validParams<stzExampleFunction>();

class stzExampleFunction : public Function
{
public:
  stzExampleFunction(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p);

protected:
  Real _a0;
  Real _a1;
  Real _t0;
  Real _t1;
};

#endif //STZEXAMPLEFUNCTION_H
