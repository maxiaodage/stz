/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "anaHEVPFlowRatePowerLawJ2.h"

template<>
InputParameters validParams<anaHEVPFlowRatePowerLawJ2>()
{
  InputParameters params = validParams<HEVPFlowRateUOBase>();
  params.addParam<Real>("reference_flow_rate",0.001,"Reference flow rate for rate dependent flow");
  params.addParam<Real>("flow_rate_exponent",10.0,"Power law exponent in flow rate equation");
  params.addParam<Real>("flow_rate_tol", 1e3, "Tolerance for flow rate");
  params.addClassDescription("User object to evaluate power law flow rate and flow direction based on J2");

  return params;
}

anaHEVPFlowRatePowerLawJ2::anaHEVPFlowRatePowerLawJ2(const InputParameters & parameters) :
    HEVPFlowRateUOBase(parameters),
    _ref_flow_rate(getParam<Real>("reference_flow_rate")),
    _flow_rate_exponent(getParam<Real>("flow_rate_exponent")),
    _flow_rate_tol(getParam<Real>("flow_rate_tol")),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>(_base_name + "deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "deformation_gradient")),
    _rotation_increment(getMaterialProperty<RankTwoTensor>(_base_name + "rotation_increment"))
{
}

bool
anaHEVPFlowRatePowerLawJ2::computeValue(unsigned int qp, Real & val) const
{
  RankTwoTensor pk2_dev = computePK2Deviatoric(_pk2[qp], _ce[qp]);
  Real eqv_stress = computeEqvStress(pk2_dev, _ce[qp]);
  val = std::pow(eqv_stress/_strength[qp], _flow_rate_exponent) * _ref_flow_rate;

  if (val > _flow_rate_tol)
  {
// #ifdef DEBUG
//     mooseWarning("Flow rate greater than " , _flow_rate_tol << " " << val << " " << eqv_stress << " " << _strength[qp]);
// #endif
    return false;
  }
  return true;
}

bool
anaHEVPFlowRatePowerLawJ2::computeDirection(unsigned int qp, RankTwoTensor & val) const
{
  RankTwoTensor pk2_dev = computePK2Deviatoric(_pk2[qp], _ce[qp]);
  Real eqv_stress = computeEqvStress(pk2_dev, _ce[qp]);
  RankTwoTensor fe = _deformation_gradient[qp];
  Real fe_det = fe.det();

  RankTwoTensor cauchy = 1.0/fe_det*fe*_pk2[qp]*fe.transpose();

  val.zero();
  if (eqv_stress > 0.0)
    val = 1.5/eqv_stress * _ce[qp] * pk2_dev * _ce[qp];

    std::vector<Real> eigval(3, 0.0);
    RankTwoTensor eigvec;
    pk2_dev.symmetricEigenvaluesEigenvectors(eigval, eigvec);

    // // If the elastic strain is beyond the cracking strain, save the eigen vectors as
    // // the rotation tensor.
    // principal_strain(0, 0) = eigval[0];
    // principal_strain(1, 0) = eigval[1];
    // principal_strain(2, 0) = eigval[2];

    RankTwoTensor R = eigvec;
    RankTwoTensor ePrime = R.transpose() * pk2_dev * R ; //elastic_strain tensor in principial coordinate


  return true;
}

bool
anaHEVPFlowRatePowerLawJ2::computeDerivative(unsigned int qp, const std::string & coupled_var_name, Real & val) const
{
  val = 0.0;

  if (_strength_prop_name == coupled_var_name)
  {
    RankTwoTensor pk2_dev = computePK2Deviatoric(_pk2[qp], _ce[qp]);
    Real eqv_stress = computeEqvStress(pk2_dev, _ce[qp]);
    val = - _ref_flow_rate * _flow_rate_exponent * std::pow(eqv_stress/_strength[qp],_flow_rate_exponent)/_strength[qp];
  }

  return true;
}

bool
anaHEVPFlowRatePowerLawJ2::computeTensorDerivative(unsigned int qp, const std::string & coupled_var_name, RankTwoTensor & val) const
{
  val.zero();

  if (_pk2_prop_name == coupled_var_name)
  {
    RankTwoTensor pk2_dev = computePK2Deviatoric(_pk2[qp], _ce[qp]);
    Real eqv_stress = computeEqvStress(pk2_dev, _ce[qp]);
    Real dflowrate_dseqv = _ref_flow_rate * _flow_rate_exponent * std::pow(eqv_stress/_strength[qp],_flow_rate_exponent-1.0)/_strength[qp];

    RankTwoTensor tau = pk2_dev * _ce[qp];
    RankTwoTensor dseqv_dpk2dev;

    dseqv_dpk2dev.zero();
    if (eqv_stress > 0.0)
      dseqv_dpk2dev = 1.5/eqv_stress * tau * _ce[qp];

    RankTwoTensor ce_inv = _ce[qp].inverse();

    RankFourTensor dpk2dev_dpk2;
    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
      for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
        for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
          for (unsigned int l = 0; l < LIBMESH_DIM; ++l)
          {
            dpk2dev_dpk2(i, j, k, l) = 0.0;
            if (i==k && j==l)
              dpk2dev_dpk2(i, j, k, l) = 1.0;
            dpk2dev_dpk2(i, j, k, l) -= ce_inv(i, j) * _ce[qp](k, l)/3.0;
          }
    val = dflowrate_dseqv * dpk2dev_dpk2.transposeMajor() * dseqv_dpk2dev;
  }
  return true;
}

RankTwoTensor
anaHEVPFlowRatePowerLawJ2::computePK2Deviatoric(const RankTwoTensor & pk2, const RankTwoTensor & ce) const
{
  return pk2 - (pk2.doubleContraction(ce) * ce.inverse())/3.0;
}

Real
anaHEVPFlowRatePowerLawJ2::computeEqvStress(const RankTwoTensor & pk2_dev, const RankTwoTensor & ce) const
{
  RankTwoTensor sdev = pk2_dev * ce;
  Real val = sdev.doubleContraction(sdev.transpose());
  return std::pow(1.5 * val, 0.5);
}
