/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZSTRESSDIVERGENCETENSORS_H
#define STZSTRESSDIVERGENCETENSORS_H

#include "ALEKernel.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

//Forward Declarations
class stzStressDivergenceTensors;
class RankTwoTensor;
class RankFourTensor;

template<>
InputParameters validParams<stzStressDivergenceTensors>();

/**
 * StressDivergenceTensors mostly copies from StressDivergence.  There are small changes to use
 * RankFourTensor and RankTwoTensors instead of SymmElasticityTensors and SymmTensors.  This is done
 * to allow for more mathematical transparancy.
 */
class stzStressDivergenceTensors : public ALEKernel
{
public:
  stzStressDivergenceTensors(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  virtual void computeJacobian();
  virtual void computeOffDiagJacobian(unsigned int jvar);

  virtual void computeFiniteDeformJacobian();

  std::string _base_name;
  bool _use_finite_deform_jacobian;

  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankFourTensor> & _Jacobian_mult;

  std::vector<RankFourTensor> _finite_deform_Jacobian_mult;
  const MaterialProperty<RankTwoTensor> * _deformation_gradient;
  const MaterialProperty<RankTwoTensor> * _deformation_gradient_old;
  const MaterialProperty<RankTwoTensor> * _rotation_increment;
  // MaterialProperty<RankTwoTensor> & _d_stress_dT;

  const unsigned int _component;

  /// Coupled displacement variables
  unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;

  const bool _temp_coupled;

  const unsigned int _temp_var;
};

#endif //STRESSDIVERGENCETENSORS_H
