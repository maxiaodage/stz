/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZVOLFLOWRATE_H
#define STZVOLFLOWRATE_H

#include "HEVPFlowRateUOBase.h"

class stzVOLFlowRate;

template<>
InputParameters validParams<stzVOLFlowRate>();

/**
 * This user object classs
 * Computes flow rate based on power law and
 * Direction based on J2
 */
class stzVOLFlowRate : public HEVPFlowRateUOBase
{
public:
  stzVOLFlowRate(const InputParameters & parameters);

  virtual bool computeValue(unsigned int, Real &) const;
  virtual bool computeDirection(unsigned int, RankTwoTensor &) const;
  virtual bool computeDerivative(unsigned int, const std::string &, Real &) const;
  virtual bool computeTensorDerivative(unsigned int, const std::string &, RankTwoTensor &) const;

protected:

//
const MaterialProperty<Real> & _volplc;




  RankTwoTensor computePK2Deviatoric(const RankTwoTensor &, const RankTwoTensor &) const;
  RankTwoTensor computePK2Hydrostatic(const RankTwoTensor &, const RankTwoTensor &) const;
  Real computeEqvStress(const RankTwoTensor &, const RankTwoTensor &) const;


};

#endif
