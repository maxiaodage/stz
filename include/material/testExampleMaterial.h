/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef TESTEXAMPLEMATERIAL_H
#define TESTEXAMPLEMATERIAL_H

#include "Material.h"
#include "LinearInterpolation.h"
#include "RankTwoTensor.h"

//Forward Declarations
class testExampleMaterial;

template<>
InputParameters validParams<testExampleMaterial>();

/**
 * Example material class that defines a few properties.
 */
class testExampleMaterial : public Material
{
public:
  testExampleMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;


  //private:

  MaterialProperty<Real> & _sourcet;
  MaterialProperty<Real> & _diffusivity;
  MaterialProperty<Real> & _chihat;
  MaterialProperty<Real> & _volpl;
  MaterialProperty<Real> & _ratio;
  MaterialProperty<Real> & _alpha_prime;
  MaterialProperty<Real> & _alpha_primeold;
  MaterialProperty<Real> & _sonew;
  MaterialProperty<Real> & _smises;
  MaterialProperty<Real> & _nondeq;
  MaterialProperty<RankTwoTensor> & _dpl;








  std::string _strength_prop_name;
  std::string _base_name;
  // const MaterialProperty<Real> & _strength;
  std::string _pk2_prop_name;
  const MaterialProperty<RankTwoTensor> & _pk2;
  const MaterialProperty<RankTwoTensor> & _ce;
  const MaterialProperty<RankTwoTensor> & _fp;
  const MaterialProperty<RankTwoTensor> & _fp_old;
  const MaterialProperty<RankTwoTensor> & _stensor;


  //
  RankTwoTensor computePK2Deviatoric(const RankTwoTensor &, const RankTwoTensor &) const;
  RankTwoTensor computefpdot(const RankTwoTensor &, const RankTwoTensor &) const;
  RankTwoTensor computeLp(const RankTwoTensor &, const RankTwoTensor &) const;
  RankTwoTensor computeSdev(const RankTwoTensor &, const RankTwoTensor &) const;
  Real _aa;
  Real _q0;
  Real _tau;
  Real _chihatst;
  Real _qthre;
  Real _co;
  Real _alphaconst;
  Real _cohesion;
  Real _beta;
  const VariableValue &_chidot;


};

#endif //TESTEXAMPLEMATERIAL_H
