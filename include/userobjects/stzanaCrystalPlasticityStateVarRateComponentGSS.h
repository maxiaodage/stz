/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZANACRYSTALPLASTICITYSTATEVARRATECOMPONENTGSS_H
#define STZANACRYSTALPLASTICITYSTATEVARRATECOMPONENTGSS_H

#include "CrystalPlasticityStateVarRateComponent.h"
#include "RankTwoTensor.h"


class stzanaCrystalPlasticityStateVarRateComponentGSS;

template<>InputParameters validParams<stzanaCrystalPlasticityStateVarRateComponentGSS>();

/**
 * Phenomenological constitutive model state variable evolution rate component userobject class.
 */
class stzanaCrystalPlasticityStateVarRateComponentGSS : public CrystalPlasticityStateVarRateComponent
{
 public:
  stzanaCrystalPlasticityStateVarRateComponentGSS(const InputParameters & parameters);

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

#endif // STZANACRYSTALPLASTICITYSTATEVARRATECOMPONENTGSS_H
