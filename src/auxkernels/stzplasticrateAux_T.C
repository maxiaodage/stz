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

#include "stzplasticrateAux_T.h"

template <>
InputParameters
validParams<stzplasticrateAux_T>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addParam<Real>("value", 0.0, "Scalar value used for our auxiliary calculation");
  params.addParam<Real>("e0",1.5,"epsiolon_0 in the plastic strain rate term");
  params.addParam<Real>("ez",0.5,"epsiolon_z in the plastic strain rate term");
  params.addParam<Real>("e1",0.3,"epsiolon_1 in the plastic strain rate term");
  params.addParam<Real>("p",80e6," Pressure ");
  params.addParam<Real>("a",350e-6,"Characteristic grain size");
  params.addParam<Real>("density",1600.0,"Density");
  params.addParam<Real>("rho",5e-4,"Vibration intensity");
  params.addParam<Real>("s1l",0.98,"Hardening parameters s1l first parameter");
  params.addParam<Real>("s2l",-1.6,"Hardening parameters s2l 2nd parameter");
  params.addParam<Real>("R0",1.0,"STZ flipping rate");
  params.addParam<Real>("alpha_biot",1.0,"The biot Coeffieicnt for Pore Pressure");
  params.addParam<Real>("s0expswitch",1.0,"s0 expression switcher 1.0= use exp experession, 0.0 use linear expression ");
  params.addParam<Real>("chihat_const",0.3,"Steady state chi value");
  params.addCoupledVar("coupled_chi",0.0, "Coupled variable chi");
  params.addCoupledVar("coupled_s", 0.0,"Coupled variable s");
  params.addCoupledVar("couple_pf",0.0,"Coupling the pore Pressure pf");


  return params;
}

stzplasticrateAux_T::stzplasticrateAux_T(const InputParameters & parameters)
  : AuxKernel(parameters),

    // We can couple in a value from one of our kernels with a call to coupledValueAux
    _e0(getParam<Real>("e0")),
    _ez(getParam<Real>("ez")),
    _e1(getParam<Real>("e1")),
    _p(getParam<Real>("p")),
    _a(getParam<Real>("a")),
    _density(getParam<Real>("density")),
    _rho(getParam<Real>("rho")),
    _s1l(getParam<Real>("s1l")),
    _s2l(getParam<Real>("s2l")),
    _R0(getParam<Real>("R0")),
    _alpha_biot(getParam<Real>("alpha_biot")),
    _s0expswtich(getParam<Real>("s0expswitch")),
    _chihat_const(getParam<Real>("chihat_const")),
    _coupled_s(coupledScalarValue("coupled_s")),
    _coupled_chi(coupledValue("coupled_chi")),
    _coupled_pf(coupledValue("couple_pf"))

{
}

/**
 * Auxiliary Kernels override computeValue() instead of computeQpResidual().  Aux Variables
 * are calculated either one per elemenet or one per node depending on whether we declare
 * them as "Elemental (Constant Monomial)" or "Nodal (First Lagrange)".  No changes to the
 * source are necessary to switch from one type or the other.
 */
Real
stzplasticrateAux_T::computeValue()
{

  Real _peff = _p- _alpha_biot*_coupled_pf[_qp];

  Real s = _coupled_s[0];

  Real s0;
  if (_s0expswtich)
  {
    s0= _s1l*_peff+_s2l*_peff*std::exp(1.0/_coupled_chi[_qp]);
  }
  else
  {
    s0= _s1l*_peff+_s2l*_peff*std::min(_coupled_chi[_qp],_chihat_const);
  }
//  Real s0= _s1l*_peff+_s2l*_peff*std::exp(1.0/_coupled_chi[_qp]);
//  Real s0= _s1l*_p+_s2l*_p*0.2;

  // if (1.0)
  // {
  //   s0 = _s1l*_p+_s2l*_p*0.3;
  // }
  Real T_current = std::tanh(_e0*s/(_peff*_ez*_coupled_chi[_qp]));
  Real tscale = _a*std::sqrt(_density/_peff);
  Real temp = 1.0+s/s0*T_current+_rho/(2.0*_R0);
  Real sqrtn = std::abs(temp*temp-4.0*s*T_current/s0);
  //Real  m = s0/(2.0*s)*temp-s0/(2.0*s)*std::sqrt(sqrtn);
  Real m = 0.0;
  if (s/s0*T_current<1.0)
  {
    m = T_current;
  }
  else
  {
    m = s0/s;
  }
  Real mechnoise = 2.0*s/s0*_R0*(T_current-m);

  Real plasticStrainRate = 4.0*_e0*std::exp(-1.0/_coupled_chi[_qp])*_R0*(T_current-m)/tscale;

//  return _coupled_s[_qp] + _value;
  return T_current;
}
