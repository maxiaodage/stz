 /****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "stzanaCrystalPlasticityStateVariableCouple.h"

template<>
InputParameters validParams<stzanaCrystalPlasticityStateVariableCouple>()
{
  InputParameters params = validParams<CrystalPlasticityStateVariable>();
  params.addParam<FileName>("state_variable_file_name", "", "Name of the file containing the initial values of slip system resistances");
  MooseEnum intvar_read_options("file_input inline_input ","inline_input");
  params.addParam<MooseEnum>("intvar_read_type", intvar_read_options, "Read from options for initial value of internal variables: Default from .i file");
  params.addParam<std::vector<unsigned int> >("groups", "To group the initial values on different slip systems 'format: [start end)', i.e.'0 4 8 11' groups 0-3, 4-7 and 8-11 ");
  params.addParam<std::vector<Real> >("group_values", "The initial values correspoinding to each group, i.e. '0.0 1.0 2.0' means 0-4 = 0.0, 4-8 = 1.0 and 8-12 = 2.0 ");
  params.addParam<std::vector<std::string> >("uo_state_var_evol_rate_comp_name", "Name of state variable evolution rate component property: Same as state variable evolution rate component user object specified in input file.");
  params.addParam<Real>("zero", 0.0 ,"Numerical zero for interval variable");
  params.addParam<Real>("xlength", 0.0 ,"Inclusion x length");
  params.addParam<Real>("ylength", 0.0 ,"Inclusion y length");
  params.addParam<Real>("zlength", 0.0 ,"Inclusion z length");
  params.addParam<Real>("inclusionvalue", 0.08,"Inclusion value of internal variable");
  params.addParam<Real>("backgroundvalue", 0.04,"Background value of internal variable");




  params.addParam<std::vector<Real> >("scale_factor", "Scale factor of individual component.");
  params.addClassDescription("Crystal plasticity state variable class.  Override the virtual functions in your class");
  return params;
}

stzanaCrystalPlasticityStateVariableCouple::stzanaCrystalPlasticityStateVariableCouple(const InputParameters & parameters) :
    CrystalPlasticityStateVariable(parameters),
    _num_mat_state_var_evol_rate_comps(parameters.get<std::vector<std::string> >("uo_state_var_evol_rate_comp_name").size()),
    _mat_prop_state_var(getMaterialProperty<std::vector<Real> >(_name)),
    _mat_prop_state_var_old(getMaterialPropertyOld< std::vector<Real> >(_name)),
    _state_variable_file_name(getParam<FileName>("state_variable_file_name")),
    _intvar_read_type(getParam<MooseEnum>("intvar_read_type")),
    _groups(getParam<std::vector<unsigned int> >("groups")),
    _group_values(getParam<std::vector<Real> >("group_values")),
    _zero(getParam<Real>("zero")),
    _scale_factor(getParam<std::vector<Real> >("scale_factor")),
    _xlength(getParam<Real>("xlength")),
    _ylength(getParam<Real>("ylength")),
    _zlength(getParam<Real>("zlength")),
    _inclusionvalue(getParam<Real>("inclusionvalue")),
    _backgroundvalue(getParam<Real>("backgroundvalue"))



{
  if (_scale_factor.size() != _num_mat_state_var_evol_rate_comps)
    mooseError("anaCrystalPlasticityStateVariable: Scale factor should be have the same size of evolution rate components.");

  _mat_prop_state_var_evol_rate_comps.resize(_num_mat_state_var_evol_rate_comps);

  for (unsigned int i = 0; i < _num_mat_state_var_evol_rate_comps; ++i)
    _mat_prop_state_var_evol_rate_comps[i] = &getMaterialProperty<std::vector<Real> >(parameters.get<std::vector<std::string> >("uo_state_var_evol_rate_comp_name")[i]);
}

void
stzanaCrystalPlasticityStateVariableCouple::initSlipSysProps(std::vector<Real> & val, Real xcoord, Real ycoord, Real zcoord) const
{

  switch (_intvar_read_type)
  {
    case 0:
      readInitialValueFromFile(val);
      break;
    case 1:
      readInitialValueFromInline(val,xcoord,ycoord,zcoord);


      // if ( xcoord< 1 && ycoord<1)
      // {
      //   std::cout<<"xcoord"<<xcoord<<"y"<<ycoord<<"z"<<zcoord<<"\n";
      //
      // }

      break;
    default:
      mooseError("anaCrystalPlasticityStateVariable: Read option for initial value of internal variables is not supported.");
  }

  // for (unsigned int i = 0; i < _variable_size; ++i)
  //   if (val[i] <= 0.0)
  //     mooseError("anaCrystalPlasticityStateVariable: Value of state variables " << i  << " non positive");
}

void
stzanaCrystalPlasticityStateVariableCouple::readInitialValueFromFile(std::vector<Real> & val) const
{
  MooseUtils::checkFileReadable(_state_variable_file_name);


  std::ifstream file;
  file.open(_state_variable_file_name.c_str());

  for (unsigned int i = 0; i < _variable_size; ++i)
    if (!(file >> val[i]))
      mooseError("Error anaCrystalPlasticityStateVariable: Premature end of state_variable file");

  file.close();
}

void
stzanaCrystalPlasticityStateVariableCouple::readInitialValueFromInline(std::vector<Real> & val,Real xcoord, Real ycoord, Real zcoord) const
{
//  std::cout<<"xcoord"<<xcoord<<"y"<<ycoord<<"z"<<zcoord<<"\n";





  for (unsigned int i = 0; i <_variable_size; ++i)
  {
        //if ((std::abs(xcoord)<_xlength)&& (std::abs(ycoord)<_ylength))
     if ((std::pow(std::abs(xcoord),2.0)+std::pow(std::abs(ycoord),2.0))<_xlength*_xlength)
        {
          val[i] = _inclusionvalue;
        }
        else
        {
          val[i] = _backgroundvalue;
        }
  }

}

bool
stzanaCrystalPlasticityStateVariableCouple::updateStateVariable(unsigned int qp, Real dt, std::vector<Real> & val) const
{
  // Anand
  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    val[i] = 0.0;
    for (unsigned int j = 0; j < _num_mat_state_var_evol_rate_comps; j++)
    {      // Anand
       val[i] += (*_mat_prop_state_var_evol_rate_comps[j])[qp][i] * dt ;

      //  val[i] += (*_mat_prop_state_var_evol_rate_comps[j])[qp][i] * dt * _scale_factor[j];
    }
  }

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    if (_mat_prop_state_var_old[qp][i] < _zero && val[i] < 0.0)
      val[i] = _mat_prop_state_var_old[qp][i];
    else
      val[i] = _mat_prop_state_var_old[qp][i] + val[i];

    if (val[i] < 0.0)
      return false;
  }
  return true;
}
