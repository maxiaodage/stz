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

#include "stzImplicitODEs.h"

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template <>
InputParameters
validParams<stzImplicitODEs>()
{
  InputParameters params = validParams<ODEKernel>();
  //params.addCoupledVar("y", "variable Y coupled into this kernel");
  params.addParam<Real>("shearG",109.6e6,"Shear modulus");
  params.addParam<Real>("imprate",100.0,"Imposing strain rate (V/H)");
  params.addParam<Real>("L",0.02,"Length of the domain");
  params.addParam<PostprocessorName>("postprocessor_rpl", 0.0, "Postprocessor to coupling the integral of the plastic strainrate ");
  return params;
}

stzImplicitODEs::stzImplicitODEs(const InputParameters & parameters)
  : // You must call the constructor of the base class first
    ODEKernel(parameters),
    // get the coupled variable number and values
    // _y_var(coupledScalar("y")),
    // _y(coupledScalarValue("y")),
    _shearG(getParam<Real>("shearG")),
    _imprate(getParam<Real>("imprate")),
    _L(getParam<Real>("L")),
    _intrpl(getPostprocessorValue("postprocessor_rpl"))

{
}

Real
stzImplicitODEs::computeQpResidual()
{
  // the term of the ODE without the time derivative term
//  return -3. * _u[_i];
// if (_intrpl>0.0)
// {
// std::cout<<"_intrpl="<<_intrpl<<"\n";
// int a;
// std::cin >> a;
// }
 return -(_shearG*(_imprate-_intrpl/_L));


}

Real
stzImplicitODEs::computeQpJacobian()
{
  // dF/dx
  return 0.0;
}

Real
stzImplicitODEs::computeQpOffDiagJacobian(unsigned int jva)
{
  // if (jvar == _y_var)
  //   return -2.; // dF/dy
  // else
    return 0.; // everything else
}
