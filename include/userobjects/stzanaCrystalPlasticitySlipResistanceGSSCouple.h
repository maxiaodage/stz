/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZANACRYSTALPLASTICITYSLIPRESISTANCEGSSCOUPLE_H
#define STZANACRYSTALPLASTICITYSLIPRESISTANCEGSSCOUPLE_H

#include "CrystalPlasticitySlipResistance.h"

class stzanaCrystalPlasticitySlipResistanceGSSCouple;

template<>
InputParameters validParams<stzanaCrystalPlasticitySlipResistanceGSSCouple>();

/**
 * Phenomenological constitutive model slip resistance userobject class.
 */
class stzanaCrystalPlasticitySlipResistanceGSSCouple : public CrystalPlasticitySlipResistance
{
 public:
  stzanaCrystalPlasticitySlipResistanceGSSCouple(const InputParameters & parameters);

  virtual bool calcSlipResistance(unsigned int qp, std::vector<Real> & val) const;

 protected:
  const MaterialProperty<std::vector<Real> > & _mat_prop_state_var;
};

#endif // STZANACRYSTALPLASTICITYSLIPRESISTANCEGSSCOUPLE_H
