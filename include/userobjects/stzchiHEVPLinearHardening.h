/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZCHIHEVPLINEARHARDENING_H
#define STZCHIHEVPLINEARHARDENING_H

#include "HEVPStrengthUOBase.h"
#include "AuxKernel.h"
#include "RankTwoTensor.h"
class stzchiHEVPLinearHardening;

template<>
InputParameters validParams<stzchiHEVPLinearHardening>();

/**
 * This user object classs
 * Computes linear hardening
 */
class stzchiHEVPLinearHardening : public HEVPStrengthUOBase
{
public:
  stzchiHEVPLinearHardening(const InputParameters & parameters);

  virtual bool computeValue(unsigned int, Real &) const;
  virtual bool computeDerivative(unsigned int, const std::string &, Real &) const;

//protected:
  Real _chihat;
  Real _chi_0;
  Real _miu_1;
  Real _miu_2;
  //adding tensor stress
  const VariableValue &_chiv;

};

#endif
