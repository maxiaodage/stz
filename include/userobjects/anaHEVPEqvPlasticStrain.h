/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ANAHEVPEQVPLASTICSTRAIN_H
#define ANAHEVPEQVPLASTICSTRAIN_H

#include "HEVPInternalVarUOBase.h"

class anaHEVPEqvPlasticStrain;

template<>
InputParameters validParams<anaHEVPEqvPlasticStrain>();

/**
 * This user object classs
 * Computes equivalent plastic strain
 */
class anaHEVPEqvPlasticStrain : public HEVPInternalVarUOBase
{
public:
  anaHEVPEqvPlasticStrain(const InputParameters & parameters);

  virtual bool computeValue(unsigned int, Real, Real &) const;
  virtual bool computeDerivative(unsigned int, Real, const std::string &, Real &) const;
};

#endif
