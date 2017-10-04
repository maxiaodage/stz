/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STZANAFINITESTRAINUOBASEDCP_H
#define STZANAFINITESTRAINUOBASEDCP_H

#include "ComputeStressBase.h"

// #include "stzanaFiniteStrainUObasedCP.h"
// #include "petscblaslapack.h"

// #include "anaCrystalPlasticitySlipRateGSS.h"
// #include "anaCrystalPlasticitySlipResistanceGSS.h"
// #include "anaCrystalPlasticityStateVariable.h"
// #include "anaCrystalPlasticityStateVarRateComponentGSS.h"

#include "stzanaCrystalPlasticitySlipRateGSS.h"
#include "stzanaCrystalPlasticitySlipResistanceGSS.h"
#include "stzanaCrystalPlasticityStateVariable.h"
#include "stzanaCrystalPlasticityStateVarRateComponentGSS.h"


class stzanaFiniteStrainUObasedCP;

template<>
InputParameters validParams<stzanaFiniteStrainUObasedCP>();

/**
 * FiniteStrainUObasedCP uses the multiplicative decomposition of deformation gradient
 * and solves the PK2 stress residual equation at the intermediate configuration to evolve the material state.
 * The internal variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 *
 * Involves 4 different types of user objects that calculates:
 * State variables - update state variable (derive from CrystalPlasticityStateVariable)
 * State variable evolution compoment - individual component of state variable incremental rate (derive from CrystalPlasticityStateVariableEvolutionRateComponent)
 * Slip resistance - calcuate slip resistances (derive from CrystalPlasticitySlipResistances)
 * Slip rates - calcuate flow direction and slip rates (derive from CrystalPlasticitySlipRates)
 */
class stzanaFiniteStrainUObasedCP : public ComputeStressBase
{
public:
  stzanaFiniteStrainUObasedCP(const InputParameters & parameters);

protected:
  /**
   * updates the stress at a quadrature point.
   */
  virtual void computeQpStress();

  /**
   * initializes the stateful properties such as
   * stress, plastic deformation gradient, slip system resistances, etc.
   */
  virtual void initQpStatefulProperties();

  /**
   * calls the residual and jacobian functions used in the
   * stress update algorithm.
   */
  virtual void calcResidJacob();

  /**
   * updates the slip system resistances and state variables.
   * override to modify slip system resistance and state variable evolution.
   */
  virtual void updateSlipSystemResistanceAndStateVariable();

  /**
   * set variables for stress and internal variable solve.
   */
  virtual void preSolveQp();

  /**
   * solve stress and internal variables.
   */
  virtual void solveQp();

  /**
   * update stress and internal variable after solve.
   */
  virtual void postSolveQp();

  /**
   * set variables for internal variable solve.
   */
  virtual void preSolveStatevar();

  /**
   * solve internal variables.
   */
  virtual void solveStatevar();

  /**
   * update internal variable after solve.
   */
  virtual void postSolveStatevar();

  /**
   * set variables for stress solve.
   */
  virtual void preSolveStress();

  /**
   * solves for stress, updates plastic deformation gradient.
   */
  virtual void solveStress();

  /**
   * update stress and plastic deformation gradient after solve.
   */
  virtual void postSolveStress();

  /**
   * calculate stress residual.
   */
  virtual void calcResidual();

  /**
   * calculate jacobian.
   */
  virtual void calcJacobian();

  /**
   * updates the slip rates.
  */
  virtual void getSlipRates();

  /**
   * calculate the tangent moduli for preconditioner.
   * Default is the elastic stiffness matrix.
   * Exact jacobian is currently implemented.
   * tan_mod_type can be modified to exact in .i file to turn it on.
   */
  virtual void calcTangentModuli();

  /**
   * calculate the elastic tangent moduli for preconditioner.
   */
  virtual void elasticTangentModuli();

  /**
   * calculate the exact tangent moduli for preconditioner.
   */
  virtual void elastoPlasticTangentModuli();

  /**
   * performs the line search update
   */
  bool lineSearchUpdate(const Real rnorm_prev, const RankTwoTensor);

  /**
   * evaluates convergence of state variables.
   */
  virtual bool isStateVariablesConverged();

  /// User objects that define the slip rate
  std::vector<const stzanaCrystalPlasticitySlipRateGSS *> _uo_slip_rates;

  /// User objects that define the slip resistance
  std::vector<const stzanaCrystalPlasticitySlipResistanceGSS *> _uo_slip_resistances;

  /// User objects that define the state variable
  std::vector<const stzanaCrystalPlasticityStateVariable *> _uo_state_vars;

  /// User objects that define the state variable evolution rate component
  std::vector<const stzanaCrystalPlasticityStateVarRateComponentGSS *> _uo_state_var_evol_rate_comps;

  /// Slip rates material property
  std::vector<MaterialProperty<std::vector<Real> > *> _mat_prop_slip_rates;

  /// Slip resistance material property
  std::vector<MaterialProperty<std::vector<Real> > *> _mat_prop_slip_resistances;

  /// State variable material property
  std::vector<MaterialProperty<std::vector<Real> > *> _mat_prop_state_vars;

  /// Old state variable material property
  std::vector<MaterialProperty<std::vector<Real> > *> _mat_prop_state_vars_old;

  /// State variable evolution rate component material property
  std::vector<MaterialProperty<std::vector<Real> > *> _mat_prop_state_var_evol_rate_comps;

  /// Number of slip rate user objects
  unsigned int _num_uo_slip_rates;

  /// Number of slip resistance user objects
  unsigned int _num_uo_slip_resistances;

  /// Number of state variable user objects
  unsigned int _num_uo_state_vars;

  /// Number of state variable evolution rate component user objects
  unsigned int _num_uo_state_var_evol_rate_comps;

  /// Local state variable
  std::vector<std::vector<Real> > _state_vars_old;

  /// Local old state variable
  std::vector<std::vector<Real> > _state_vars_prev;

  /// Stress residual equation relative tolerance
  Real _rtol;
  /// Stress residual equation absolute tolerance
  Real _abs_tol;
  /// Internal variable update equation tolerance
  Real _stol;
  /// Residual tolerance when variable value is zero. Default 1e-12.
  Real _zero_tol;

  /// Residual tensor
  RankTwoTensor _resid;

  /// Jacobian tensor
  RankFourTensor _jac;

  /// Maximum number of iterations for stress update
  unsigned int _maxiter;
  /// Maximum number of iterations for internal variable update
  unsigned int _maxiterg;

  /// Type of tangent moduli calculation
  MooseEnum _tan_mod_type;

  /// Maximum number of substep iterations
  unsigned int _max_substep_iter;

  /// Flag to activate line serach
  bool _use_line_search;

  /// Minimum line search step size
  Real _min_lsrch_step;

  /// Line search bisection method tolerance
  Real _lsrch_tol;

  /// Line search bisection method maximum iteration number
  unsigned int _lsrch_max_iter;

  /// Line search method
  MooseEnum _lsrch_method;

  MaterialProperty<RankTwoTensor> & _fp;
  MaterialProperty<RankTwoTensor> & _fp_old;
  MaterialProperty<RankTwoTensor> & _pk2;
  MaterialProperty<RankTwoTensor> & _pk2_old;
  MaterialProperty<RankTwoTensor> & _lag_e;
  MaterialProperty<RankTwoTensor> & _eurl_e;

  // Output Material Property eigen direction
  MaterialProperty<RankTwoTensor> & _egen_dir;
  MaterialProperty<RankTwoTensor> & _pstrainrate;
  MaterialProperty<Real> & _equiprate;
  MaterialProperty<Real> & _tau1;
  MaterialProperty<Real> & _sigma1;
  MaterialProperty<Real> & _tau2;
  MaterialProperty<Real> & _sigma2;
  MaterialProperty<Real> & _tau3;
  MaterialProperty<Real> & _sigma3;
  MaterialProperty<Real> & _yields1;
  MaterialProperty<Real> & _yields2;
  MaterialProperty<Real> & _yields3;
  MaterialProperty<Real> & _eig1;
  MaterialProperty<Real> & _eig2;
  MaterialProperty<Real> & _eig3;





  MaterialProperty<RankTwoTensor> & _update_rot;
  MaterialProperty<RankTwoTensor> & _update_rot_old;

  const MaterialProperty<RankTwoTensor> & _deformation_gradient;
  const MaterialProperty<RankTwoTensor> & _deformation_gradient_old;

  /// Crystal rotation
  const MaterialProperty<RankTwoTensor> & _crysrot;

  RankTwoTensor _dfgrd_tmp;
  RankTwoTensor _fe, _fp_old_inv, _fp_inv;
  DenseVector<Real> _tau;
  // Adding Normal
  DenseVector<Real> _sigma;

  std::vector<MaterialProperty<std::vector<RankTwoTensor> > * > _flow_direction;

  std::vector<MaterialProperty<std::vector<RankTwoTensor> > * > _flow_direction_normal;

  std::vector<MaterialProperty<std::vector<RankTwoTensor> > * > _flow_direction_s;




  /// Flag to check whether convergence is achieved
  bool _err_tol;

  /// Used for substepping; Uniformly divides the increment in deformation gradient
  RankTwoTensor _delta_dfgrd, _dfgrd_tmp_old;
  /// Scales the substepping increment to obtain deformation gradient at a substep iteration
  Real _dfgrd_scale_factor;

  Real _mu_fric;
  Real _cohesion_fric_p1;
  Real _cohesion_fric_p2;
  Real _cohesion_fric_p3;



};

#endif //STZANAFINITESTRAINUOBASEDCP_H
