/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "stzchisource_singleblock.h"
#include "Function.h"

template<>
InputParameters validParams<stzchisource_singleblock>()
{
  InputParameters params = validParams<BodyForce>();
  params.addClassDescription("Kernel that defines a body force modified by a material mask");
  params.addParam<Real>("dbgswitch",1.0,"dbgswitch for output 1.0");
  params.addParam<Real>("chihat_const",0.3,"Steady state chi value");
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
  params.addParam<Real>("eps0",4.0e-3,"Frictional Constant");
  params.addParam<Real>("tfr",1.3e-2,"Frictional Time Scale");
  params.addParam<Real>("r_pl_thre",10.0,"r* the threshold for the chi hat formulation");
  params.addParam<Real>("q0",1.0,"The q0 term in the chi hat formulation");
  params.addParam<Real>("alpha_biot",1.0,"The biot Coeffieicnt for Pore Pressure");
  params.addParam<Real>("s0expswitch",1.0,"s0 expression switcher 1.0= use exp experession, 0.0 use linear expression ");
  params.addCoupledVar("couples",0.0,"Coupling s");
  params.addCoupledVar("couple_rpl",0.0,"Couping the plastic strain rate");
  params.addCoupledVar("couple_pf",0.0,"Coupling the pore Pressure pf");



  return params;
}

stzchisource_singleblock::stzchisource_singleblock(const InputParameters & parameters) :
    BodyForce(parameters),
    _dbgswitch(getParam<Real>("dbgswitch")),
    _chihat_const(getParam<Real>("chihat_const")),
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
    _eps0(getParam<Real>("eps0")),
    _tfr(getParam<Real>("tfr")),
    _r_pl_thre(getParam<Real>("r_pl_thre")),
    _q0(getParam<Real>("q0")),
    _alpha_biot(getParam<Real>("alpha_biot")),
    _s0expswtich(getParam<Real>("s0expswitch")),
   // Adding coupled s
   _s(coupledValue("couples")),
   _couple_rpl(coupledValue("couple_rpl")),
   _couple_pf(coupledValue("couple_pf"))

{
}

Real
stzchisource_singleblock::computeQpResidual()
{

 Real _peff = _p- _alpha_biot*_couple_pf[_qp];
  Real s = _s[_qp];

Real s0;
  if (_s0expswtich)
  {
    s0= _s1l*_peff+_s2l*_peff*std::exp(1.0/_u[_qp]);
  }
  else
  {
    s0= _s1l*_peff+_s2l*_peff*std::min(_u[_qp],_chihat_const);
  }
// Real s0= _s1l*_p+_s2l*_p*std::exp(1.0/_u[_qp]);
 //Real s0= _s1l*_peff+_s2l*_peff*0.2;
  // if (1.)
  // {
  //   s0 = _s1l*_p+_s2l*_p*0.3;
  // }

  Real T_current = std::tanh(_e0*s/(_peff*_ez*_u[_qp]));
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
  Real eps = _eps0*std::tanh((_tfr*_couple_rpl[_qp])*(_tfr*_couple_rpl[_qp]));
  Real bigc=2.0*_e0*s0/(_peff*tscale*_e1)*std::exp(-1.0/_u[_qp])*mechnoise;

  Real chihat=_chihat_const;

  if (_couple_rpl[_qp]<_r_pl_thre)
  {
    chihat = _chihat_const;
  }
  else
  {
    Real A0  = _chihat_const*std::log(_q0/(_r_pl_thre*tscale));
    chihat = A0/(std::log(_q0/(_couple_rpl[_qp]*tscale)));
  }


//  eps = eps0*tanh((t_fric*plasticStrainRate).^2);
 //2.0*e0*s0./(p-alpha*pf).*exp(-1.0./chi)/(tscale*e1)

  if (_dbgswitch)
  {
    std::cout<<"******************"<<"\n"\
    <<"CHIsource"<<"\n"\
    <<"s="<<s<<"\n"\
    <<"s0="<<s0<<"\n"\
    <<"T_current="<<T_current<<"\n"\
    <<"tscale="<<tscale<<"\n"\
    <<"temp="<<temp<<"\n"\
    <<"sqrtn="<<sqrtn<<"\n"\
    <<"m="<<m<<"\n"\
    <<"R0="<<_R0<<"\n"\
    <<"mechnoise="<<mechnoise<<"\n"\
    <<"bigc"<<bigc<<"\n"\
    <<"return="<<-_test[_i][_qp]*(1.0-_u[_qp]/chihat)*bigc<<"\n"\
    <<"peff"<<_peff<<"\n";
  }

  // return -_test[_i][_qp]*(1.0-_u[_qp]/_chihat)*bigc;
  return -_test[_i][_qp]*(2.0*_e0*s0/(_peff*tscale*_e1)*std::exp(-1.0/_u[_qp])*(mechnoise*(1.0-_u[_qp]/chihat)-eps*_u[_qp]/chihat));










//  Real r_dotpl = _e0/_tau0*std::exp(_sv[_qp]/(_sc*_u[_qp]))*std::exp(-1.0/_u[_qp])*(1.0-_s0/_sv[_qp]);


}
