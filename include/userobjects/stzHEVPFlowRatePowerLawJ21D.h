/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZHEVPFLOWRATEPOWERLAWJ21D_H
#define STZHEVPFLOWRATEPOWERLAWJ21D_H

#include "HEVPFlowRateUOBase.h"


class stzHEVPFlowRatePowerLawJ21D;

template<>
InputParameters validParams<stzHEVPFlowRatePowerLawJ21D>();

/**
 * This user object classs
 * Computes flow rate based on power law and
 * Direction based on J2
 */
class stzHEVPFlowRatePowerLawJ21D : public HEVPFlowRateUOBase
{
public:
  stzHEVPFlowRatePowerLawJ21D(const InputParameters & parameters);

  virtual bool computeValue(unsigned int, Real &) const;
  virtual bool computeDirection(unsigned int, RankTwoTensor &) const;
  virtual bool computeDerivative(unsigned int, const std::string &, Real &) const;
  virtual bool computeTensorDerivative(unsigned int, const std::string &, RankTwoTensor &) const;

protected:
  Real _ref_flow_rate;
  Real _flow_rate_exponent;
  Real _flow_rate_tol;
  Real _tolplastic;
// Adding the grain pressure
  Real _grain_pressure;
// Adding the chi
  Real _chi;
// Adding the tau
  Real _t_tau;
// Density
  Real _rho;
// Biot coefficient
  Real _biot;

//
  const VariableValue &_chiv;
  const VariableValue &_agrain;
  const VariableValue &_pf;


  RankTwoTensor computePK2Deviatoric(const RankTwoTensor &, const RankTwoTensor &) const;
  Real computeEqvStress(const RankTwoTensor &, const RankTwoTensor &) const;


};

#endif
