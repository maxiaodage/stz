/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZHEVPLINEARHARDENING1D_H
#define STZHEVPLINEARHARDENING1D_H

#include "HEVPStrengthUOBase.h"
#include "AuxKernel.h"
#include "RankTwoTensor.h"
class stzHEVPLinearHardening1D;

template<>
InputParameters validParams<stzHEVPLinearHardening1D>();

/**
 * This user object classs
 * Computes linear hardening
 */
class stzHEVPLinearHardening1D : public HEVPStrengthUOBase
{
public:
  stzHEVPLinearHardening1D(const InputParameters & parameters);

  virtual bool computeValue(unsigned int, Real &) const;
  virtual bool computeDerivative(unsigned int, const std::string &, Real &) const;

//protected:
  Real _sig0;
  Real _s1l;
  Real _s2l;
  Real _p;
  const VariableValue &_chiv;





};

#endif
