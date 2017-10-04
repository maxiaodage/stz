/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "stzanaCrystalPlasticityStateVarRateComponentGSSCouple.h"

template<>
InputParameters validParams<stzanaCrystalPlasticityStateVarRateComponentGSSCouple>()
{
  InputParameters params = validParams<CrystalPlasticityStateVarRateComponent>();
  params.addParam<std::string>("uo_slip_rate_name", "Name of slip rate property: Same as slip rate user object specified in input file.");
  params.addParam<std::string>("uo_state_var_name", "Name of state variable property: Same as state variable user object specified in input file.");
  params.addParam<FileName>("slip_sys_hard_prop_file_name", "", "Name of the file containing the values of hardness evolution parameters");
  params.addParam<std::vector<Real> >("hprops", "Hardening properties");
  params.addParam<Real>("p", 0.8, "p parameter");
  params.addParam<Real>("g0",0.04, "g0 parameter");
  params.addParam<Real>("eta_cv", 100, "eta_cv parameter");
  params.addParam<Real>("chihat", 0.12, "STZ chi Steady State value");
  params.addParam<Real>("c0", 0.01, "STZ c0 volumetric expansion coefficient");
  params.addParam<Real>("mu_fricrate", 0.0, "Tangent of Friction angle");
  params.addParam<Real>("cohesion_fricrate_p1", 1.0e6, "Cohesion rate plane 12");
  params.addParam<Real>("cohesion_fricrate_p2", 1.0e6, "Cohesion rate plane 34");
  params.addParam<Real>("cohesion_fricrate_p3", 1.0e6, "Cohesion rate plane 56");
  params.addParam<Real>("pressure", 10.0e6, "Confining Pressre");






  params.addClassDescription("Phenomenological constitutive model state variable evolution rate component base class.  Override the virtual functions in your class");
  return params;
}

stzanaCrystalPlasticityStateVarRateComponentGSSCouple::stzanaCrystalPlasticityStateVarRateComponentGSSCouple(const InputParameters & parameters) :
    CrystalPlasticityStateVarRateComponent(parameters),
    _mat_prop_slip_rate(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_slip_rate_name"))),
    _mat_prop_state_var(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_state_var_name"))),
    _stress(getMaterialPropertyByName<RankTwoTensor>("stress")),
    _pk2(getMaterialPropertyByName<RankTwoTensor>("pk2")),

    _pstrainrate(getMaterialPropertyByName<RankTwoTensor>("pstrainrate")),
    _slip_sys_hard_prop_file_name(getParam<FileName>("slip_sys_hard_prop_file_name")),
    _hprops(getParam<std::vector<Real> >("hprops")),
    _p(getParam<Real>("p")),
    _g0(getParam<Real>("g0")),
    _eta_cv(getParam<Real>("eta_cv")),
    _chihat(getParam<Real>("chihat")),
    _c0(getParam<Real>("c0")),
    _mu_fricrate(getParam<Real>("mu_fricrate")),
    _cohesion_fricrate_p1(getParam<Real>("cohesion_fricrate_p1")),
    _cohesion_fricrate_p2(getParam<Real>("cohesion_fricrate_p2")),
    _cohesion_fricrate_p3(getParam<Real>("cohesion_fricrate_p3")),
    _pressure(getParam<Real>("pressure"))







{
}

bool
stzanaCrystalPlasticityStateVarRateComponentGSSCouple::calcStateVariableEvolutionRateComponent(unsigned int qp, std::vector<Real> & val) const
{
  val.assign(_variable_size, 0.0);
    RankTwoTensor s_dev;
    s_dev = _pk2[qp].deviatoric();
    // Real sdpl;
    // sdpl= s_dev.doubleContraction(_pstrainrate[qp]);
    Real sigma_mean;
    //sigma_mean= 1.0/3.0*_stress[qp].trace();
    sigma_mean= 1.0/3.0*_pk2[qp].trace();

      // std::cout<<"Cauchy Stress= "<<"\n"\
      //  <<_stress[qp](0,0)<<"|"<<_stress[qp](0,1)<<"|"<<_stress[qp](0,2)<<"\n"\
      //  <<_stress[qp](1,0)<<"|"<<_stress[qp](1,1)<<"|"<<_stress[qp](1,2)<<"\n"\
      //  <<_stress[qp](2,0)<<"|"<<_stress[qp](1,1)<<"|"<<_stress[qp](2,2)<<"\n"\
      //     <<"Pk2 Stress= "<<"\n"\
      //     <<_pk2[qp](0,0)<<"|"<<_pk2[qp](0,1)<<"|"<<_pk2[qp](0,2)<<"\n"\
      //     <<_pk2[qp](1,0)<<"|"<<_pk2[qp](1,1)<<"|"<<_pk2[qp](1,2)<<"\n"\
      //     <<_pk2[qp](2,0)<<"|"<<_pk2[qp](1,1)<<"|"<<_pk2[qp](2,2)<<"\n"\
      //     <<"Dpl= "<<"\n"\
      //     <<_pstrainrate[qp](0,0)<<"|"<<_pstrainrate[qp](0,1)<<"|"<<_pstrainrate[qp](0,2)<<"\n"\
      //     <<_pstrainrate[qp](1,0)<<"|"<<_pstrainrate[qp](1,1)<<"|"<<_pstrainrate[qp](1,2)<<"\n"\
      //     <<_pstrainrate[qp](2,0)<<"|"<<_pstrainrate[qp](1,1)<<"|"<<_pstrainrate[qp](2,2)<<"\n";
// Anand
//   Real vscalar;
//   vscalar = 0.0 ;
//    for (unsigned int i = 0; i < _variable_size; ++i)
//    {
// //     vscalar += _mat_prop_slip_rate[qp][i];
//    vscalar = vscalar+_mat_prop_slip_rate[qp][i];
//      //std::cout<<"i"<<i<<"vscalr"<<vscalar<<"v_alpha"<<_mat_prop_slip_rate[qp][i]<<"\n";
//
//    }
DenseVector<Real> tau_rate(_variable_size);
DenseVector<Real> sigma_rate(_variable_size);
DenseVector<Real> yield_rate(_variable_size);

std::vector<Real> eigval_ratenew(3, 0.0);
std::vector<Real> eigval_rate(3, 0.0);
//
RankTwoTensor eigvec_rate;
//
Real theta ;
Real PI = 3.1415926;
//Real mu = 0.04 ;
Real phi = std::atan(_mu_fricrate);
theta =PI/4.0+phi/2.0;
_pk2[qp].symmetricEigenvaluesEigenvectors(eigval_ratenew, eigvec_rate);
//  _stress[qp].symmetricEigenvaluesEigenvectors(eigval_new, eigvec);
//
eigval_rate[0]=eigval_ratenew[2];
eigval_rate[1]=eigval_ratenew[1];
eigval_rate[2]=eigval_ratenew[0];
// // Using cuachy stress
//
// //_stress[qp].symmetricEigenvaluesEigenvectors(eigval, eigvec);
//
// 1 2
tau_rate(0) = 0.5*std::sin(2*theta)*(eigval_rate[0]-eigval_rate[2]);
sigma_rate(0) = -0.5* (eigval_rate[0]+eigval_rate[2])+0.5*std::cos(2*theta)*(eigval_rate[0]-eigval_rate[2]);
yield_rate(0) = _cohesion_fricrate_p1+_mu_fricrate*sigma_rate(0);

tau_rate(1) = tau_rate(0);
sigma_rate(1) = sigma_rate(0);
yield_rate(1) = _cohesion_fricrate_p1+_mu_fricrate*sigma_rate(1);
// 3 4
tau_rate(2) = 0.5*std::sin(2*theta)*(eigval_rate[0]-eigval_rate[1]);
sigma_rate(2) = -0.5* (eigval_rate[0]+eigval_rate[1])+0.5*std::cos(2*theta)*(eigval_rate[0]-eigval_rate[1]);
yield_rate(2) = _cohesion_fricrate_p2+_mu_fricrate*sigma_rate(2);

tau_rate(3) = tau_rate(2);
sigma_rate(3) = sigma_rate(2);
yield_rate(3) = _cohesion_fricrate_p2+_mu_fricrate*sigma_rate(3);
// 5 6
tau_rate(4) = 0.5*std::sin(2*theta)*(eigval_rate[1]-eigval_rate[2]);
sigma_rate(4) = -0.5* (eigval_rate[1]+eigval_rate[2])+0.5*std::cos(2*theta)*(eigval_rate[1]-eigval_rate[2]);
yield_rate(4) = _cohesion_fricrate_p3+_mu_fricrate*sigma_rate(4);

tau_rate(5) = tau_rate(4);
sigma_rate(5) = sigma_rate(4);
yield_rate(5) = _cohesion_fricrate_p3+_mu_fricrate*sigma_rate(5);

// Avoding sigma_rate become tensio

// for (unsigned int i = 0; i<_variable_size; ++i)
// {
//   if (sigma_rate(i)<=0.0)
//   {
//     sigma_rate(i)=0.0;
//   }
// }

//
Real sdpl = 0.0;
//sdpl=0.0;
   for (unsigned int i = 0; i < _variable_size; ++i)
   {
//     vscalar += _mat_prop_slip_rate[qp][i];
//if ((tau_rate(i)-std::abs(_cohesion_fricrate+_mu_fricrate*sigma_rate(i)))>=0.0)

    if (yield_rate(i)<=0.0)
    {
      yield_rate(i) = 0.0;
    }

//if ((tau_rate(i))>(cohesion_fricrate+_mu_fricrate*sigma_rate(i)))
if ((tau_rate(i))>=(yield_rate(i)))
{
  sdpl = sdpl+_mat_prop_slip_rate[qp][i]*tau_rate(i);
  //   std::cout<<"tau_rate(i)="<<i<<"\n"\
  //   <<tau_rate(0)<<"|"<<tau_rate(1)<<"|"<<"\n"\
  //   <<tau_rate(2)<<"|"<<tau_rate(3)<<"|"<<"\n"\
  //   <<tau_rate(4)<<"|"<<tau_rate(5)<<"|"<<"\n"\
  //   <<"yield_rate(i)="<<"\n"\
  //  <<yield_rate(0)<<"|"<<yield_rate(1)<<"|"<<"\n"\
  //  <<yield_rate(2)<<"|"<<yield_rate(3)<<"|"<<"\n"\
  //  <<yield_rate(4)<<"|"<<yield_rate(5)<<"|"<<"\n";
}
// else
// {
//   sdpl = sdpl+0.0;
// }


//  std::cout<<"tau(i)="<<tau_rate(i)<<"mu*sigma(i="<<_mu_fricrate*sigma_rate(i)+_cohesion_fricrate<<"\n";
// if (qp==0)
// {
//   std::cout<<"tau(1)="<<tau_rate(0)<<"mu*sigma(1)+cohesion="<<_mu_fricrate*sigma_rate(0)+_cohesion_fricrate<<"\n";
//   //<<"tau(2)="<<tau_rate(2)<<"mu*sigma(2)+cohesion="<<_mu_fricrate*sigma_rate(2)+_cohesion_fricrate<<"\n"\
//   <<"tau(3)="<<tau_rate(4)<<"mu*sigma(3)+cohesion="<<_mu_fricrate*sigma_rate(4)+_cohesion_fricrate<<"\n";
// }


     //std::cout<<"i"<<i<<"vscalr"<<vscalar<<"v_alpha"<<_mat_prop_slip_rate[qp][i]<<"\n";

   }

   for (unsigned int i = 0; i < _variable_size; ++i)
   {
    //  val[i] =_g0 *std::pow((1.0-_mat_prop_state_var[qp][i]/_eta_cv),_p) * vscalar;
    //
    //  if ((_mat_prop_state_var[qp][i]-_eta_cv)>=0.0)
    //  {
    //    val[i]=0.0;
    //  }

  //  val[i] = (std::abs(sdpl)/(_c0*std::abs(sigma_mean)))*(1.0-_mat_prop_state_var[qp][i]/_chihat);
    // Debug
   val[i] =(std::abs(sdpl))/(_c0*_pressure)*(1.0-_mat_prop_state_var[qp][i]/_chihat);
  // val[i] =(1.0)/(_c0)*(1.0-_mat_prop_state_var[qp][i]/_chihat);

    //val[i] = 0.0;

  // std::cout<<"i="<<i<<"sdpl"<<sdpl<<"sigma_mean="<<sigma_mean<<"chi_rate="<<val[i]<<"chi="<<_mat_prop_state_var[qp][i]<<"\n";

    if (_mat_prop_state_var[qp][i]>=_chihat)
    {
      val[i]=0.0;
    }




   }

   //


  return true;
}
