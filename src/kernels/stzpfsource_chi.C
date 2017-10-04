/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "stzpfsource_chi.h"
#include "Function.h"

template<>
InputParameters validParams<stzpfsource_chi>()
{
  InputParameters params = validParams<BodyForce>();
  params.addClassDescription("Kernel that defines Chi dot term in Pore Pressure");
  params.addParam<Real>("alpha",1.0,"alpha parameter in front of chi_dot");
  params.addParam<Real>("beta",1.0e9,"beta parameter in front of chi_dot");
  params.addCoupledVar("couple_chi",1e-4,"Coupling chi for using in the chi_dot in pore pressure") ;



  return params;
}

stzpfsource_chi::stzpfsource_chi(const InputParameters & parameters) :
    BodyForce(parameters),
    _alpha(getParam<Real>("alpha")),
    _beta(getParam<Real>("beta")),
    _chidot(coupledDot("couple_chi"))


{
}

Real
stzpfsource_chi::computeQpResidual()
{


 Real chidot= _chidot[_qp];
Real V_dot  = (_alpha*chidot);


return _test[_i][_qp]*(_beta*V_dot);








}
