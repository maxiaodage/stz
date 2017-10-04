/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZANACRYSTALPLASTICITYSLIPRATEGSSCOUPLE_H
#define STZANACRYSTALPLASTICITYSLIPRATEGSSCOUPLE_H

#include "CrystalPlasticitySlipRate.h"
#include "RankTwoTensor.h"

class stzanaCrystalPlasticitySlipRateGSSCouple;

template<>
InputParameters validParams<stzanaCrystalPlasticitySlipRateGSSCouple>();

/**
 * Phenomenological constitutive model slip rate userobject class.
 */
class stzanaCrystalPlasticitySlipRateGSSCouple : public CrystalPlasticitySlipRate
{
 public:
  stzanaCrystalPlasticitySlipRateGSSCouple(const InputParameters & parameters);

  virtual bool calcSlipRate(unsigned int qp, Real dt, std::vector<Real> & val) const override;
//  virtual bool
  virtual bool calcSlipRateDerivative(unsigned int qp, Real /*dt*/, std::vector<Real> & val) const override;
  //virtual void calcFlowDirection(unsigned int qp, std::vector<RankTwoTensor> & flow_direction) const;
//  virtual
  void calcFlowDirection(unsigned int qp, std::vector<RankTwoTensor> & flow_direction) const;
  void calcFlowDirectionnormal(unsigned int qp, std::vector<RankTwoTensor> & flow_direction_normal) const;
  void calcFlowDirections(unsigned int qp, std::vector<RankTwoTensor> & flow_direction_s) const;


 protected:
  virtual void readFileFlowRateParams();
  virtual void getFlowRateParams();

  const MaterialProperty<std::vector<Real> > & _mat_prop_state_var;

  const MaterialProperty<RankTwoTensor> & _pk2;
  // Adding cauchy
  const MaterialProperty<RankTwoTensor> & _stress;


  DenseVector<Real> _a0;
  DenseVector<Real> _xm;

  const MaterialProperty<std::vector<RankTwoTensor> > & _flow_direction;
  const MaterialProperty<std::vector<RankTwoTensor> > & _flow_direction_normal;
  const MaterialProperty<std::vector<RankTwoTensor> > & _flow_direction_s;



  Real _mu;
  Real _c_cv;
  Real _b;
  Real _q;
  Real _eta_cv;
  Real _g0;
  Real _p;
  Real _v0;
  Real _m0;
  Real _mode;
  Real _e0;
  Real _tau0;
  Real _pressure_exp;
  Real _cohesion_p1;
  Real _cohesion_p2;
  Real _cohesion_p3;
  Real _exptermcap;
  const VariableValue &_chiv;






};

#endif // STZANACRYSTALPLASTICITYSLIPRATEGSSCOUPLE_H
