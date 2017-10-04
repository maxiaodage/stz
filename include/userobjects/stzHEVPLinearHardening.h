/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZHEVPLINEARHARDENING_H
#define STZHEVPLINEARHARDENING_H

#include "HEVPStrengthUOBase.h"
#include "AuxKernel.h"
#include "RankTwoTensor.h"
class stzHEVPLinearHardening;

template<>
InputParameters validParams<stzHEVPLinearHardening>();

/**
 * This user object classs
 * Computes linear hardening
 */
class stzHEVPLinearHardening : public HEVPStrengthUOBase
{
public:
  stzHEVPLinearHardening(const InputParameters & parameters);

  virtual bool computeValue(unsigned int, Real &) const;
  virtual bool computeDerivative(unsigned int, const std::string &, Real &) const;

//protected:
  Real _sig0;
  Real _slope;
  Real _tanfric;
  //adding tensor stress
 const MaterialProperty<RankTwoTensor> & _tensor;
 const MaterialProperty<Real> & _volplc;
 const MaterialProperty<Real> & _diffc;
 Real _beta;
 const MaterialProperty<Real> & _alphanew;
 const MaterialProperty<Real> & _sohydro;



};

#endif
