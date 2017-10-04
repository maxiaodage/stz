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

#ifndef TESTEXAMPLEMATERIAL1D_H
#define TESTEXAMPLEMATERIAL1D_H

#include "Material.h"
#include "LinearInterpolation.h"
#include "RankTwoTensor.h"

//Forward Declarations
class testExampleMaterial1D;

template<>
InputParameters validParams<testExampleMaterial1D>();

/**
 * Example material class that defines a few properties.
 */
class testExampleMaterial1D : public Material
{
public:
  testExampleMaterial1D(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;


  //private:

  MaterialProperty<Real> & _sourcet;
  MaterialProperty<Real> & _diffusivity;
  MaterialProperty<RankTwoTensor> & _dpl;








  std::string _strength_prop_name;
  std::string _base_name;
  // const MaterialProperty<Real> & _strength;
  std::string _pk2_prop_name;
  const MaterialProperty<RankTwoTensor> & _pk2;
  const MaterialProperty<RankTwoTensor> & _ce;
  const MaterialProperty<RankTwoTensor> & _fp;
  const MaterialProperty<RankTwoTensor> & _fp_old;


  //
  RankTwoTensor computePK2Deviatoric(const RankTwoTensor &, const RankTwoTensor &) const;
  RankTwoTensor computefpdot(const RankTwoTensor &, const RankTwoTensor &) const;
  RankTwoTensor computeLp(const RankTwoTensor &, const RankTwoTensor &) const;
  RankTwoTensor computeSdev(const RankTwoTensor &, const RankTwoTensor &) const;



};

#endif //TESTEXAMPLEMATERIAL_H
