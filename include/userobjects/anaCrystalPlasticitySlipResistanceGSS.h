/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ANACRYSTALPLASTICITYSLIPRESISTANCEGSS_H
#define ANACRYSTALPLASTICITYSLIPRESISTANCEGSS_H

#include "CrystalPlasticitySlipResistance.h"

class anaCrystalPlasticitySlipResistanceGSS;

template<>
InputParameters validParams<anaCrystalPlasticitySlipResistanceGSS>();

/**
 * Phenomenological constitutive model slip resistance userobject class.
 */
class anaCrystalPlasticitySlipResistanceGSS : public CrystalPlasticitySlipResistance
{
 public:
  anaCrystalPlasticitySlipResistanceGSS(const InputParameters & parameters);

  virtual bool calcSlipResistance(unsigned int qp, std::vector<Real> & val) const;

 protected:
  const MaterialProperty<std::vector<Real> > & _mat_prop_state_var;
};

#endif // ANACRYSTALPLASTICITYSLIPRESISTANCEGSS_H
