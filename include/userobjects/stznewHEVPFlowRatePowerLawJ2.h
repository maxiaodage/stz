/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZNEWHEVPFLOWRATEPOWERLAWJ2_H
#define STZNEWHEVPFLOWRATEPOWERLAWJ2_H

#include "HEVPFlowRateUOBase.h"


class stznewHEVPFlowRatePowerLawJ2;

template<>
InputParameters validParams<stznewHEVPFlowRatePowerLawJ2>();

/**
 * This user object classs
 * Computes flow rate based on power law and
 * Direction based on J2
 */
class stznewHEVPFlowRatePowerLawJ2 : public HEVPFlowRateUOBase
{
public:
  stznewHEVPFlowRatePowerLawJ2(const InputParameters & parameters);

  virtual bool computeValue(unsigned int, Real &) const;
  virtual bool computeDirection(unsigned int, RankTwoTensor &) const;
  virtual bool computeDerivative(unsigned int, const std::string &, Real &) const;
  virtual bool computeTensorDerivative(unsigned int, const std::string &, RankTwoTensor &) const;

protected:
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

// Friction angle
  Real _phiangle;

//
  const VariableValue &_chiv;
  const VariableValue &_agrain;
  const VariableValue &_pf;

  //adding tensor stress
 const MaterialProperty<RankTwoTensor> & _tensor;


  RankTwoTensor computePK2Deviatoric(const RankTwoTensor &, const RankTwoTensor &) const;
  Real computeEqvStress(const RankTwoTensor &, const RankTwoTensor &) const;


};

#endif
