/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ANAHEVPLINEARHARDENING_H
#define ANAHEVPLINEARHARDENING_H

#include "HEVPStrengthUOBase.h"

class anaHEVPLinearHardening;

template<>
InputParameters validParams<anaHEVPLinearHardening>();

/**
 * This user object classs
 * Computes linear hardening
 */
class anaHEVPLinearHardening : public HEVPStrengthUOBase
{
public:
  anaHEVPLinearHardening(const InputParameters & parameters);

  virtual bool computeValue(unsigned int, Real &) const;
  virtual bool computeDerivative(unsigned int, const std::string &, Real &) const;

protected:
  Real _sig0;
  Real _slope;
};

#endif
