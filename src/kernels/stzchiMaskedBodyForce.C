/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "stzchiMaskedBodyForce.h"
#include "Function.h"

template<>
InputParameters validParams<stzchiMaskedBodyForce>()
{
  InputParameters params = validParams<BodyForce>();
  params.addClassDescription("Kernel that defines a body force modified by a material mask");
  params.addParam<Real>("pressure",25e6,"Grain Pressure");
  params.addParam<Real>("chihat",0.3,"Steady State Compactivity");
  params.addParam<Real>("chi_0",0.18,"Initial Compactivity");
  params.addParam<Real>("tau",8e-7,"Time Scale");
  params.addParam<Real>("epsilon_0",1.5,"Grain Pressure");
  params.addParam<Real>("epsilon_1",0.3,"Grain Pressure");
  params.addParam<Real>("miu_1",10.0e6,"Inital Yield Strength");
  params.addParam<Real>("miu_2",2.5e6,"Steady State Yield Strength");
  params.addParam<Real>("tau_f",0.013,"Time Scale at Compactional effect of interpaticle friction becomes prominent");
  params.addParam<Real>("ksi_0",0.0,"Grain Pressure");
  return params;
}

stzchiMaskedBodyForce::stzchiMaskedBodyForce(const InputParameters & parameters) :
    BodyForce(parameters),
    _sdpl(getMaterialProperty<Real>("sourcet")),
    _diffusivity(getMaterialProperty<Real>("diffusivity")),
    _pressure(getParam<Real>("pressure")),
    _chihat(getParam<Real>("chihat")),
    _chi_0(getParam<Real>("chi_0")),
    _tau(getParam<Real>("tau")),
    _epsilon_0(getParam<Real>("epsilon_0")),
    _epsilon_1(getParam<Real>("epsilon_1")),
    _miu_1(getParam<Real>("miu_1")),
    _miu_2(getParam<Real>("miu_2")),
    _tau_f(getParam<Real>("tau_f")),
    _ksi_0(getParam<Real>("ksi_0"))
{
}

Real
stzchiMaskedBodyForce::computeQpResidual()
{

Real _slope = (_miu_1-_miu_2)/(_chi_0-_chihat);
Real _miu_0 ;
   if (_u[_qp]<=_chi_0)
      {
      _miu_0 = _miu_1;
      }
  else if (_u[_qp]>_chi_0 && _u[_qp]<_chihat)
     {
       _miu_0 = _slope*_u[_qp]+ (_miu_2-_slope*_chihat);
     }
  else if (_u[_qp]>=_chihat)
    {
          _miu_0 = _miu_2 ;
    }
  // STZ density
Real Lamda = std::exp(-1/_u[_qp]);
Real bigW = (2*_epsilon_0*_miu_0*Lamda)/(_tau*_epsilon_1);
Real bigR = (_tau)/(2*_epsilon_0*_miu_0*Lamda);
Real bigC = 1.0/_pressure*_sdpl[_qp];
Real bigM = _ksi_0*std::tanh((_tau_f*_diffusivity[_qp])*(_tau_f*_diffusivity[_qp]));
   //Real bigc=1.0/_pressure;
   return -_test[_i][_qp]*(bigW*(bigR*bigC*(1-_u[_qp]/_chihat))-bigM*_u[_qp]/_chihat);
}
