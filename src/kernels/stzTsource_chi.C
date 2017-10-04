/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "stzTsource_chi.h"
#include "Function.h"

template<>
InputParameters validParams<stzTsource_chi>()
{
  InputParameters params = validParams<BodyForce>();
  params.addClassDescription("Kernel that defines Chi source term in Temperature equations");
  params.addParam<Real>("rhoc",1.0,"Specific heat Capacity rho_c term in the Temeprateure source term");
  params.addParam<Real>("chihat_const",0.3,"Steady state chi value");
  params.addParam<Real>("r_pl_thre",10.0,"r* the threshold for the chi hat formulation");
  params.addParam<Real>("q0",1.0,"The q0 term in the chi hat formulation");
  params.addParam<Real>("p",80e6," Pressure ");
  params.addParam<Real>("a",350e-6,"Characteristic grain size");
  params.addParam<Real>("density",1600.0,"Density");
  params.addParam<Real>("alpha_biot",1.0,"The biot Coeffieicnt for Pore Pressure");
  params.addParam<Real>("chiswitch",1.0,"1.0 = Chi source term is on, 0.0 = off");
  params.addCoupledVar("couples",0.0,"Coupling s");
  params.addCoupledVar("couple_chi",0.0,"Couping chi");
  params.addCoupledVar("couple_rpl",0.0,"Couping the plastic strain rate");
  params.addCoupledVar("couple_pf",0.0,"Couping the pore pressure pf");
  return params;
}

stzTsource_chi::stzTsource_chi(const InputParameters & parameters) :
    BodyForce(parameters),
    _rhoc(getParam<Real>("rhoc")),
    _chihat_const(getParam<Real>("chihat_const")),
    _q0(getParam<Real>("q0")),
    _p(getParam<Real>("p")),
    _a(getParam<Real>("a")),
    _density(getParam<Real>("density")),
    _chiswitch(getParam<Real>("chiswitch")),
    _s(coupledScalarValue("couples")),
    _couple_chi(coupledScalarValue("couple_chi")),
    _couple_rpl(coupledValue("couple_rpl")),
    _couple_pf(coupledValue("couple_pf"))

{
}

Real
stzTsource_chi::computeQpResidual()
{
   Real _peff = _p- _alpha_biot*_couple_pf[_qp];

  Real tscale = _a*std::sqrt(_density/_peff);

Real chihat = 0.0;
  if (_couple_rpl[_qp]<_r_pl_thre)
  {
    chihat = _chihat_const;
  }
  else
  {
    Real A0  = _chihat_const*std::log(_q0/(_r_pl_thre*tscale));
    chihat = A0/(std::log(_q0/(_couple_rpl[_qp]*tscale)));
  }
Real term= _chiswitch*_s[_qp]*_couple_rpl[_qp]/(_rhoc)*(_couple_chi[_qp])/chihat;


  // return -_test[_i][_qp]*(1.0-_u[_qp]/_chihat)*bigc;
  return -_test[_i][_qp]*(term);












}
