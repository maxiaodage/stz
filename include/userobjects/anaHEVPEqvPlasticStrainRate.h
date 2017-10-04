/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ANAHEVPEQVPLASTICSTRAINRATE_H
#define ANAHEVPEQVPLASTICSTRAINRATE_H

#include "HEVPInternalVarRateUOBase.h"

class anaHEVPEqvPlasticStrainRate;

template<>
InputParameters validParams<anaHEVPEqvPlasticStrainRate>();

/**
 * This user object classs
 * Computes equivalent plastic strain rate
 */
class anaHEVPEqvPlasticStrainRate : public HEVPInternalVarRateUOBase
{
public:
  anaHEVPEqvPlasticStrainRate(const InputParameters & parameters);

  virtual bool computeValue(unsigned int, Real &) const;
  virtual bool computeDerivative(unsigned int, const std::string &, Real &) const;

protected:
  Real _h;
};

#endif
