/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ANACRYSTALPLASTICITYSTATEVARIABLE_H
#define ANACRYSTALPLASTICITYSTATEVARIABLE_H

#include "CrystalPlasticityStateVariable.h"

class anaCrystalPlasticityStateVariable;

template<>
InputParameters validParams<anaCrystalPlasticityStateVariable>();

/**
 * Crystal plasticity state variable userobject class.
 */
class anaCrystalPlasticityStateVariable : public CrystalPlasticityStateVariable
{
public:
  anaCrystalPlasticityStateVariable(const InputParameters & parameters);

  virtual bool updateStateVariable(unsigned int qp, Real dt, std::vector<Real> & val) const;
  // Adding qp points
  virtual void initSlipSysProps(std::vector<Real> & val, Real xcoord, Real ycoord, Real zcoord) const;

protected:
  virtual void readInitialValueFromFile(std::vector<Real> & val) const;

  virtual void readInitialValueFromInline(std::vector<Real> & val, Real xcoord, Real ycoord, Real zcoord) const;

  unsigned int _num_mat_state_var_evol_rate_comps;

  std::vector<const MaterialProperty<std::vector<Real> > * > _mat_prop_state_var_evol_rate_comps;

  const MaterialProperty<std::vector<Real> > &  _mat_prop_state_var;
  const MaterialProperty<std::vector<Real> > &  _mat_prop_state_var_old;

  /// File should contain initial values of the state variable.
  FileName _state_variable_file_name;

  /// Read from options for initial values of internal variables
  MooseEnum _intvar_read_type;

  /** The _groups variable is used to group slip systems and assign the initial values to each group.
   *  The format is taken as [start end)
   *  i.e. _groups = '0 4 8 11', it means three groups 0-3, 4-7 and 8-11
   */
  std::vector<unsigned int> _groups;

  /** The _group_values are the initial values corresponding to each group.
   *  i.e. _groups = '0 4 8 11', and _group_values = '1.0 2.0 3.0'
   *  it means that initial values of slip system 0-3 is 1.0 , 4-7 is 2.0 and 8-11 is 3.0
   */
  std::vector<Real> _group_values;

  /// Numerical zero for internal variable
  Real _zero;



  /// Scale factor of individual component
  std::vector<Real> _scale_factor;

  Real _xlength;
  Real _ylength;
  Real _zlength;
  Real _inclusionvalue;
  Real _backgroundvalue;




};

#endif // ANACRYSTALPLASTICITYSTATEVARIABLE_H
