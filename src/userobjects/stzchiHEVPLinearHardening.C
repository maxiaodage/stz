/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "stzchiHEVPLinearHardening.h"
#include "AuxKernel.h"
#include "RankTwoTensor.h"

template<>
InputParameters validParams<stzchiHEVPLinearHardening>()
{
  InputParameters params = validParams<HEVPStrengthUOBase>();
  params.addParam<Real>("chihat",0.3,"Steady State Compactivity");
  params.addParam<Real>("chi_0",0.18,"Initial Compactivity");
  params.addParam<Real>("miu_1",10.0e6,"Inital Yield Strength");
  params.addParam<Real>("miu_2",2.5e6,"Steady State Yield Strength");
  params.addCoupledVar("efftemp",0.000000001,"effective temperature") ;
  params.addClassDescription("User Object for linear hardening");
  return params;
}

stzchiHEVPLinearHardening::stzchiHEVPLinearHardening(const InputParameters & parameters) :
    HEVPStrengthUOBase(parameters),
    _chihat(getParam<Real>("chihat")),
    _chi_0(getParam<Real>("chi_0")),
    _miu_1(getParam<Real>("miu_1")),
    _miu_2(getParam<Real>("miu_2")),
    _chiv(coupledValue("efftemp"))

{
}

bool
stzchiHEVPLinearHardening::computeValue(unsigned int qp, Real & val) const
{
  Real _slope = (_miu_1-_miu_2)/(_chi_0-_chihat);
  Real _miu_0 ;
     if (_chiv[qp]<=_chi_0)
        {
        _miu_0 = _miu_1;
        }
    else if (_chiv[qp]>_chi_0 && _chiv[qp]<_chihat)
       {
         _miu_0 = _slope*_chiv[qp]+ (_miu_2-_slope*_chihat);
       }
    else if (_chiv[qp]>=_chihat)
      {
            _miu_0 = _miu_2 ;
      }
    val = _miu_0 ;

//    std::cout<<"haha"<<val<<"\n";



  return true;
}

bool
stzchiHEVPLinearHardening::computeDerivative(unsigned int /*qp*/, const std::string & coupled_var_name, Real & val) const
{
  val = 0;

  // if (_intvar_prop_name == coupled_var_name)
  //   val = _slope;

  return true;
}
