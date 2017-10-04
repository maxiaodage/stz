/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ANACOMPUTEELASTICITYTENSORCP_H
#define ANACOMPUTEELASTICITYTENSORCP_H

#include "ComputeElasticityTensor.h"
#include "ElementPropertyReadFile.h"
#include "RankTwoTensor.h"
#include "RotationTensor.h"

/**
 * ComputeElasticityTensorCP defines an elasticity tensor material object for crystal plasticity.
 */
class anaComputeElasticityTensorCP : public ComputeElasticityTensor
{
public:
  anaComputeElasticityTensorCP(const InputParameters & parameters);

protected:
  virtual void computeQpElasticityTensor();

  virtual void assignEulerAngles();

  /**
   * Element property read user object
   * Presently used to read Euler angles -  see test
   */
  const ElementPropertyReadFile * _read_prop_user_object;

  MaterialProperty<RealVectorValue> & _Euler_angles_mat_prop;

  /// Crystal Rotation Matrix
  MaterialProperty<RankTwoTensor> & _crysrot;

  /// Rotation matrix
  RotationTensor _R;
};

#endif //ANACOMPUTEELASTICITYTENSORCP_H
