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

#include "stzvelstepFunction.h"

template<>
InputParameters validParams<stzvelstepFunction>()
{
  InputParameters params = validParams<Function>();
  params.addParam<Real>("a0", 1.0, "Initial amplitude");
  params.addParam<Real>("a1", 1.0, "Finial amplitude");
  params.addParam<Real>("av1", 1.0, "First Stepping to Value");
  params.addParam<Real>("av2", 1.0, "Second Stepping to Value");

  params.addParam<Real>("t0", 1.0, "Initial time");
  params.addParam<Real>("t1", 1.0, "Final time");
  params.addParam<Real>("tv1", 1.0, "Time For 1st velocity stepping");
  params.addParam<Real>("tv2", 1.0, "Time For 2nd velocity stepping");



  return params;
}

stzvelstepFunction::stzvelstepFunction(const InputParameters & parameters) :
    Function(parameters),
    _a0(getParam<Real>("a0")),
    _a1(getParam<Real>("a1")),
    _av1(getParam<Real>("av1")),
    _av2(getParam<Real>("av2")),
    _t0(getParam<Real>("t0")),
    _t1(getParam<Real>("t1")),
    _tv1(getParam<Real>("tv1")),
    _tv2(getParam<Real>("tv2"))
{}

Real
stzvelstepFunction::value(Real t, const Point & p)
{
  Real ksi= (t-_t0)/(_t1-_t0);
   // Real A=0.4;
    //Real sigma_0=1e-5;
   // Real _t0 =1e-4;


    //return _a0+_a1*((_t1-_t0)/2.0)*std::log(std::cosh((t-_t0)/((_t1-_t0)/2.0)));
  //    return _a0*_a1*std::log(std::exp((t-_t0)/_a0)+1);

   if (t<=_t0)
   {
       return 0;

   }
    else if (t>_t0&&t<=_t1)
   {
    return _a0*t+2.5*(_a1-_a0)*pow(ksi,3)*(t-_t0)+3.0*(_a0-_a1)*(t-_t0)*pow(ksi,4)-(_a0-_a1)*(t-_t0)*pow(ksi,5);
   }
    else if (t>_t1&&t<=_tv1)
   {
    return _a0*t+2.5*(_a1-_a0)*(_t1-_t0)+3.0*(_a0-_a1)*(_t1-_t0)-(_a0-_a1)*(_t1-_t0)+_a1*(t-_t1);
   }
   else if (t>_tv1&&t<=_tv2)
   {
    return _a0*t+2.5*(_a1-_a0)*(_t1-_t0)+3.0*(_a0-_a1)*(_t1-_t0)-(_a0-_a1)*(_t1-_t0)+_a1*(_tv1-_t1)+_av1*(t-_tv1);
   }
   else
   {
     Real _dispbf =_a0*t+2.5*(_a1-_a0)*(_t1-_t0)+3.0*(_a0-_a1)*(_t1-_t0)-(_a0-_a1)*(_t1-_t0)+_a1*(_tv1-_t1)+_av1*(_tv2-_tv1);
     return _dispbf+_av2*(t-_tv2);
   }
  //  else
  //  {
  //   return _a1*t;
  //  }
}
