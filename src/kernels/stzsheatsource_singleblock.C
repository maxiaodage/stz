/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "stzsheatsource_singleblock.h"
#include "Function.h"

template<>
InputParameters validParams<stzsheatsource_singleblock>()
{
  InputParameters params = validParams<BodyForce>();
  params.addClassDescription("Kernel that defines a body force modified by a material mask");
  params.addParam<Real>("dbgswitch",1.0,"dbgswitch for output 1.0");

  params.addParam<Real>("shearG",109.6e6,"Shear modulus");
  params.addParam<Real>("imprate",100.0,"Imposing strain rate (V/H)");
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
  params.addParam<Real>("L",2.0,"Length of the domain");
  params.addCoupledVar("efftemp",0.10,"coupling variable chi") ;
  params.addCoupledVar("rpl",0.0,"coupling variable rpl") ;
  return params;
}

stzsheatsource_singleblock::stzsheatsource_singleblock(const InputParameters & parameters) :
    BodyForce(parameters),
    _dbgswitch(getParam<Real>("dbgswitch")),
    _shearG(getParam<Real>("shearG")),
    _imprate(getParam<Real>("imprate")),
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
    _L(getParam<Real>("L")),
    _chi(coupledValue("efftemp")),
    _rpl(coupledValue("rpl"))
{
}

Real
stzsheatsource_singleblock::computeQpResidual()
{
  // // return -_test[_i][_qp]*(1.0-_u[_qp]/_chihat)*std::abs(_sdpl[_qp])*bigc;
  // Real s = _u[_qp];
  // Real s0= _s1l*_p+_s2l*_p*_chi[_qp];
  // Real T_current = std::tanh(_e0*s/(_p*_ez*_chi[_qp]));
  // Real tscale = _a*std::sqrt(_density/_p);
  // Real temp = 1.0+s/s0*T_current+_rho/(2.0*_R0);
  // Real sqrtn = std::abs(temp*temp-4.0*s*T_current/s0);
  // Real  m = s0/(2.0*s)*temp-s0/(2.0*s)*std::sqrt(sqrtn);
  // Real mechnoise = 2.0*s/s0*_R0*(T_current-m);
  //
  // Real plasticStrainRate = 4.0*_e0*std::exp(-1.0/_chi[_qp])*_R0*(T_current-m)/tscale;



    if (_dbgswitch)
    {
      std::cout<<"******************"<<"\n"\
      <<"Ssource"<<"\n"\
      <<"return="<<-_test[_i][_qp]*(_shearG*(_imprate-_rpl[_qp]))<<"\n";
    }

    //std::cout<<"r_dotpl"<<r_dotpl<<"\n";
  return -_test[_i][_qp]*(_shearG*(_imprate-_rpl[_qp]));








  // return 0.0;
}
