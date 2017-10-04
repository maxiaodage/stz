/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZANACRYSTALPLASTICITYSTATEVARRATECOMPONENTGSSCOUPLE_H
#define STZANACRYSTALPLASTICITYSTATEVARRATECOMPONENTGSSCOUPLE_H

#include "CrystalPlasticityStateVarRateComponent.h"
#include "RankTwoTensor.h"


class stzanaCrystalPlasticityStateVarRateComponentGSSCouple;

template<>InputParameters validParams<stzanaCrystalPlasticityStateVarRateComponentGSSCouple>();

/**
 * Phenomenological constitutive model state variable evolution rate component userobject class.
 */
class stzanaCrystalPlasticityStateVarRateComponentGSSCouple : public CrystalPlasticityStateVarRateComponent
{
 public:
  stzanaCrystalPlasticityStateVarRateComponentGSSCouple(const InputParameters & parameters);

  virtual bool calcStateVariableEvolutionRateComponent(unsigned int qp, std::vector<Real> & val) const;

 protected:
  const MaterialProperty<std::vector<Real> > &  _mat_prop_slip_rate;
  const MaterialProperty<std::vector<Real> > & _mat_prop_state_var;
  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankTwoTensor> & _pk2;

  const MaterialProperty<RankTwoTensor> & _pstrainrate;


  /// The hardening parameters in this class are read from .i file. The user can override to read from file.
  FileName _slip_sys_hard_prop_file_name;

  std::vector<Real> _hprops;

  Real _p;
  Real _g0;
  Real _eta_cv;
  Real _chihat;
  Real _c0;
  Real _mu_fricrate;
  Real _cohesion_fricrate_p1;
  Real _cohesion_fricrate_p2;
  Real _cohesion_fricrate_p3;

  Real _pressure;





};

#endif // STZANACRYSTALPLASTICITYSTATEVARRATECOMPONENTGSSCOUPLE_H
