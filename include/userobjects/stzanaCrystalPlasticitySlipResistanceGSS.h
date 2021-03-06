/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZANACRYSTALPLASTICITYSLIPRESISTANCEGSS_H
#define STZANACRYSTALPLASTICITYSLIPRESISTANCEGSS_H

#include "CrystalPlasticitySlipResistance.h"

class stzanaCrystalPlasticitySlipResistanceGSS;

template<>
InputParameters validParams<stzanaCrystalPlasticitySlipResistanceGSS>();

/**
 * Phenomenological constitutive model slip resistance userobject class.
 */
class stzanaCrystalPlasticitySlipResistanceGSS : public CrystalPlasticitySlipResistance
{
 public:
  stzanaCrystalPlasticitySlipResistanceGSS(const InputParameters & parameters);

  virtual bool calcSlipResistance(unsigned int qp, std::vector<Real> & val) const;

 protected:
  const MaterialProperty<std::vector<Real> > & _mat_prop_state_var;
};

#endif // STZANACRYSTALPLASTICITYSLIPRESISTANCEGSS_H
