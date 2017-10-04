/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZNEWHEVPLINEARHARDENING_H
#define STZNEWHEVPLINEARHARDENING_H

#include "HEVPStrengthUOBase.h"
#include "AuxKernel.h"
#include "RankTwoTensor.h"
#include "RankTwoScalarTools.h"
#include "RankTwoScalarAux.h"
class stznewHEVPLinearHardening;

template<>
InputParameters validParams<stznewHEVPLinearHardening>();

/**
 * This user object classs
 * Computes linear hardening
 */
class stznewHEVPLinearHardening : public HEVPStrengthUOBase
{
public:
  stznewHEVPLinearHardening(const InputParameters & parameters);

  virtual bool computeValue(unsigned int, Real &) const;
  virtual bool computeDerivative(unsigned int, const std::string &, Real &) const;

//protected:
  Real _phi;
  Real _cohesion;

  //adding tensor stress
 const MaterialProperty<RankTwoTensor> & _tensor;



};

#endif
