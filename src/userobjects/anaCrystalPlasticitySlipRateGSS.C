/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "anaCrystalPlasticitySlipRateGSS.h"

template<>
InputParameters validParams<anaCrystalPlasticitySlipRateGSS>()
{
  InputParameters params = validParams<CrystalPlasticitySlipRate>();
  params.addParam<std::string>("uo_state_var_name", "Name of state variable property: Same as state variable user object specified in input file.");
  params.addClassDescription("Phenomenological constitutive model slip rate class.  Override the virtual functions in your class");
  // Biot
  params.addParam<Real>("mu",0.04,"Tangent of Friction Angle");
  params.addParam<Real>("c_cv",660.0e6,"Cohesion C_cv parameters");
  params.addParam<Real>("b",300.0e6,"Cohesion b parameters");
  params.addParam<Real>("q",1.2,"Cohesion q parameters");
  params.addParam<Real>("eta_cv",100,"Cohesion eta_cv parameters");
  params.addParam<Real>("g0",0.04,"Eta g0 parameters");
  params.addParam<Real>("p",0.8,"Eta p parameters");
  params.addParam<Real>("v0",0.001,"reference plastic strain rate");
  params.addParam<Real>("m0",0.005,"rate depedent power");
  params.addParam<Real>("mode",1.0,"1.0 for tension , 0,0 for compresson");





  return params;
}

anaCrystalPlasticitySlipRateGSS::anaCrystalPlasticitySlipRateGSS(const InputParameters & parameters) :
    CrystalPlasticitySlipRate(parameters),
    _mat_prop_state_var(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_state_var_name"))),
    _pk2(getMaterialPropertyByName<RankTwoTensor>("pk2")),
    // Adding cauchy
    _stress(getMaterialPropertyByName<RankTwoTensor>("stress")),
    _a0(_variable_size),
    _xm(_variable_size),
    _flow_direction(getMaterialProperty<std::vector<RankTwoTensor> >(_name + "_flow_direction")),
    _flow_direction_normal(getMaterialProperty<std::vector<RankTwoTensor> >(_name + "_flow_direction_normal")),
    _flow_direction_s(getMaterialProperty<std::vector<RankTwoTensor> >(_name + "_flow_direction_s")),
    _mu(getParam<Real>("mu")),
    _c_cv(getParam<Real>("c_cv")),
    _b(getParam<Real>("b")),
    _q(getParam<Real>("q")),
    _eta_cv(getParam<Real>("eta_cv")),
    _g0(getParam<Real>("g0")),
    _p(getParam<Real>("p")),
    _v0(getParam<Real>("v0")),
    _m0(getParam<Real>("m0")),
    _mode(getParam<Real>("mode"))






{
  if (_slip_sys_flow_prop_file_name.length() != 0)
    readFileFlowRateParams();
  else
    getFlowRateParams();
}

void
anaCrystalPlasticitySlipRateGSS::readFileFlowRateParams()
{
  MooseUtils::checkFileReadable(_slip_sys_flow_prop_file_name);

  std::ifstream file;
  file.open(_slip_sys_flow_prop_file_name.c_str());

  std::vector<Real> vec;
  vec.resize(_num_slip_sys_flowrate_props);

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    for (unsigned int j = 0; j < _num_slip_sys_flowrate_props; ++j)
      if (!(file >> vec[j]))
        mooseError("Error CrystalPlasticitySlipRateGSS: Premature end of slip_sys_flow_rate_param file");

    _a0(i)=vec[0];
    _xm(i)=vec[1];
  }

  file.close();
}

void
anaCrystalPlasticitySlipRateGSS::getFlowRateParams()
{
  if (_flowprops.size() <= 0)
    mooseError("CrystalPlasticitySlipRateGSS: Error in reading flow rate  properties: Specify input in .i file or a slip_sys_flow_prop_file_name");

  _a0.resize(_variable_size);
  _xm.resize(_variable_size);
  for (unsigned int i = 0; i < _variable_size; ++i)
 {
   _a0(i) = _v0;
   _xm(i) = _m0;
 }
}

void
anaCrystalPlasticitySlipRateGSS::calcFlowDirection(unsigned int qp, std::vector<RankTwoTensor> & flow_direction) const
{
  DenseVector<Real> mo(LIBMESH_DIM*_variable_size),no(LIBMESH_DIM*_variable_size);
  //
  std::vector<Real> eigval(3, 0.0);
  RankTwoTensor eigvec;
  RankTwoTensor eigvec_new;

  eigvec.zero();

  // if (_mode)
  // {
  //   // tesion
  //    eigvec(0,2)=1.0;
  //    eigvec(1,1)=1.0;
  //    eigvec(2,0)=1.0;
  // }
  // else
  // {
  //  compresson
    // eigvec_new(0,0)=1.0;
    // eigvec_new(1,1)=1.0;
    // eigvec_new(2,2)=1.0;
  // // }

  // Original
  // // tesion
  // // eigvec(0,2)=1.0;
  // // eigvec(1,1)=1.0;
  // // eigvec(2,0)=1.0;
  // // compresson
  // eigvec(0,0)=1.0;
  // eigvec(1,1)=1.0;
  // eigvec(2,2)=1.0;

  Real theta ;
  Real PI = 3.1415926;
  Real phi = std::atan(_mu);
  theta =PI/4.0- phi/2.0;
  _pk2[qp].symmetricEigenvaluesEigenvectors(eigval, eigvec_new);

  //First Coloum
  eigvec(0,0)=eigvec_new(0,2);
  eigvec(1,0)=eigvec_new(1,2);
  eigvec(2,0)=eigvec_new(2,2);

// Secodnary colum
eigvec(0,1)=eigvec_new(0,1);
eigvec(1,1)=eigvec_new(1,1);
eigvec(2,1)=eigvec_new(2,1);
  //Third Coluum
  eigvec(0,2)=eigvec_new(0,0);
  eigvec(1,2)=eigvec_new(1,0);
  eigvec(2,2)=eigvec_new(2,0);


  //  if (qp==0)
  //  {
  //    std::cout<<"Eigeval="<<"\n"\
  //    <<eigval[0]<<"|"<<eigval[1]<<"|"<<eigval[2]<<"\n"\
  //    <<"Eigvect_using="<<"\n"\
  //    <<eigvec(0,0)<<"|"<<eigvec(0,1)<<"|"<<eigvec(0,2)<<"\n"\
  //    <<eigvec(1,0)<<"|"<<eigvec(1,1)<<"|"<<eigvec(1,2)<<"\n"\
  //    <<eigvec(2,0)<<"|"<<eigvec(2,1)<<"|"<<eigvec(2,2)<<"\n"\
  //    <<"Eigvect_right="<<"\n"\
  //    <<eigvec_new(0,0)<<"|"<<eigvec_new(0,1)<<"|"<<eigvec_new(0,2)<<"\n"\
  //    <<eigvec_new(1,0)<<"|"<<eigvec_new(1,1)<<"|"<<eigvec_new(1,2)<<"\n"\
  //    <<eigvec_new(2,0)<<"|"<<eigvec_new(2,1)<<"|"<<eigvec_new(2,2)<<"\n"\
  //    <<"***********************************************\n";
  //  }





  //
  RankTwoTensor I_den;
  I_den.addIa(1.0);
  mo= _mo;
  no= _no;

  std::vector<Real> s1(3, 0.0), s2(3, 0.0),s3(3, 0.0), s4(3, 0.0), s5(3, 0.0), s6(3, 0.0);
  std::vector<Real> s1_new(3, 0.0), s2_new(3, 0.0) ,s3_new(3, 0.0), s4_new(3, 0.0), s5_new(3, 0.0) , s6_new(3, 0.0);
  std::vector<Real> m1(3, 0.0), m2(3, 0.0),m3(3, 0.0), m4(3, 0.0), m5(3, 0.0), m6(3, 0.0);
  std::vector<Real> m1_new(3, 0.0), m2_new(3, 0.0) ,m3_new(3, 0.0), m4_new(3, 0.0), m5_new(3, 0.0) , m6_new(3, 0.0);

  //s2(3, 0.0),s3(3,0.0),s4(3,0.0),s5(3,0.0),s6(3,0.0);
  Real c = std::cos(theta);
  Real s = std::sin(theta);

  // S1
  s1[0]= c , s1[1]= 0.0 , s1[2]= s , m1[0]= s , m1[1]= 0.0 , m1[2]= -c;
  // S2
  s2[0]= c ,  s2[1]= 0.0, s2[2]= -s, m2[0]= s , m2[1]= 0.0, m2[2]= c;
  // S3
  s3[0]= c ,  s3[1]= s, s3[2]=  0.0, m3[0]= s ,  m3[1]= -c, m3[2]=  0.0 ;
  // S4
  s4[0]= c ,  s4[1]= -s, s4[2]= 0.0, m4[0]= s ,  m4[1]= c, m4[2]=  0.0;
  // S5
  s5[0]= 0.0 ,  s5[1]= c, s5[2]= s, m5[0]= 0.0,  m5[1]= s, m5[2]=  -c;
  // S6
  s6[0]= 0.0 ,  s6[1]= c, s6[2]= -s, m6[0]= 0.0 ,  m6[1]= s, m6[2]= c;


  for (unsigned int i = 0; i < 3; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      // s_new;
      s1_new[i]=s1_new[i]+eigvec(i, j) * s1[j];
      s2_new[i]=s2_new[i]+eigvec(i, j) * s2[j];
      s3_new[i]=s3_new[i]+eigvec(i, j) * s3[j];
      s4_new[i]=s4_new[i]+eigvec(i, j) * s4[j];
      s5_new[i]=s5_new[i]+eigvec(i, j) * s5[j];
      s6_new[i]=s6_new[i]+eigvec(i, j) * s6[j];

      // m_new;
      m1_new[i]=m1_new[i]+eigvec(i, j) * m1[j];
      m2_new[i]=m2_new[i]+eigvec(i, j) * m2[j];
      m3_new[i]=m3_new[i]+eigvec(i, j) * m3[j];
      m4_new[i]=m4_new[i]+eigvec(i, j) * m4[j];
      m5_new[i]=m5_new[i]+eigvec(i, j) * m5[j];
      m6_new[i]=m6_new[i]+eigvec(i, j) * m6[j];

    }
  }
  // Adding new flowdirection
  std::vector<Real> _beta(_variable_size, 0.0);

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    // Anand
 _beta[i]  = _g0*(std::pow((1.0-_mat_prop_state_var[qp][i]/_eta_cv),_p));
//  _beta[i] = 0.0 ;
  //std::cout<<"qp"<<qp<<"i"<<i<<"state"<<_mat_prop_state_var[qp][i]<<"beta"<<_beta<<"\n";
    //
   }
  //
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
         {
        //   flow_direction[i](j,k) = mo(i*LIBMESH_DIM+j) * no(i*LIBMESH_DIM+k)+_beta*no(i*LIBMESH_DIM+j) * no(i*LIBMESH_DIM+k);
           flow_direction[0](j,k) = s1_new[j] * m1_new[k]+_beta[0]*m1_new[j] * m1_new[k];
           flow_direction[1](j,k) = s2_new[j] * m2_new[k]+_beta[1]*m2_new[j] * m2_new[k];
           flow_direction[2](j,k) = s3_new[j] * m3_new[k]+_beta[2]*m3_new[j] * m3_new[k];
           flow_direction[3](j,k) = s4_new[j] * m4_new[k]+_beta[3]*m4_new[j] * m4_new[k];
           flow_direction[4](j,k) = s5_new[j] * m5_new[k]+_beta[4]*m5_new[j] * m5_new[k];
           flow_direction[5](j,k) = s6_new[j] * m6_new[k]+_beta[5]*m6_new[j] * m6_new[k];


         }
    }


  // std::cout<<"s1="<<"\n"\
  // <<s1[0]<<"|"<<s1[1]<<"|"<<s1[2]<<"\n"\
  // <<"m1="<<"\n"\
  // <<m1[0]<<"|"<<m1[1]<<"|"<<m1[2]<<"\n"\

  //  std::cout<<"s1_new="<<"\n"
  //  <<s1_new[0]<<"|"<<s1_new[1]<<"|"<<s1_new[2]<<"\n"\
  //  <<"m1_new="<<"\n"
  //  <<m1_new[0]<<"|"<<m1_new[1]<<"|"<<m1_new[2]<<"\n";







  // Adding flowdirection right

  //
// //  unsigned int k=0;
// //  std::cout<<"qp="<<qp<<"\n"\
//   //<<mo(k*LIBMESH_DIM+0)<<"|"<<mo(k*LIBMESH_DIM+1)<<"|"<<mo(k*LIBMESH_DIM+2)<<"\n"\
//   //<<no(k*LIBMESH_DIM+0)<<"|"<<no(k*LIBMESH_DIM+1)<<"|"<<no(k*LIBMESH_DIM+2)<<"\n"\
//   <<"******************************"<<"\n"\
  //
  // unsigned int k=0;
  // std::cout<<"flowdirection="<<"\n"\
  // <<flow_direction[k](0,0)<<"|"<<flow_direction[k](0,1)<<"|"<<flow_direction[k](0,2)<<"\n"\
  // <<flow_direction[k](1,0)<<"|"<<flow_direction[k](1,1)<<"|"<<flow_direction[k](1,2)<<"\n"\
  // <<flow_direction[k](2,0)<<"|"<<flow_direction[k](2,1)<<"|"<<flow_direction[k](2,2)<<"\n";

}

void
anaCrystalPlasticitySlipRateGSS::calcFlowDirectionnormal(unsigned int qp, std::vector<RankTwoTensor> & flow_direction_normal) const
{
  DenseVector<Real> no(LIBMESH_DIM*_variable_size);
  no = _no;
  // Update slip direction and normal with crystal orientation
  std::vector<Real> eigval(3, 0.0);
  RankTwoTensor eigvec;

  _pk2[qp].symmetricEigenvaluesEigenvectors(eigval, eigvec);

  RankTwoTensor R = eigvec;
  RankTwoTensor I_den;
  I_den.addIa(1.0);
  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      no(i*LIBMESH_DIM+j) = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        {
        //  no(i*LIBMESH_DIM+j) = no(i*LIBMESH_DIM+j) + _crysrot[qp](j,k) * _no(i*LIBMESH_DIM+k);
        no(i*LIBMESH_DIM+j) = no(i*LIBMESH_DIM+j) + I_den(j,k) * _no(i*LIBMESH_DIM+k);

          //std::cout<<"qp="<<qp<<"rot="<<_crysrot[qp](j,k)<<"j="<<j<<"k"<<k<<"\n";
        }
    }
  }

  // Calculate Schmid tensor and resolved shear stresses
  for (unsigned int i = 0; i < _variable_size; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        flow_direction_normal[i](j,k) = no(i*LIBMESH_DIM+j) * no(i*LIBMESH_DIM+k);
        //Anand
}
void
anaCrystalPlasticitySlipRateGSS::calcFlowDirections(unsigned int qp, std::vector<RankTwoTensor> & flow_direction_s) const
{
  DenseVector<Real> mo(LIBMESH_DIM*_variable_size),no(LIBMESH_DIM*_variable_size);
  // Update slip direction and normal with crystal orientation
  std::vector<Real> eigval(3, 0.0);
  RankTwoTensor eigvec;

  _pk2[qp].symmetricEigenvaluesEigenvectors(eigval, eigvec);

  RankTwoTensor R = eigvec;
  RankTwoTensor I_den;
  I_den.addIa(1.0);
  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      mo(i*LIBMESH_DIM+j) = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
      {
        //mo(i*LIBMESH_DIM+j) = mo(i*LIBMESH_DIM+j) + _crysrot[qp](j,k) * _mo(i*LIBMESH_DIM+k);
        mo(i*LIBMESH_DIM+j) = mo(i*LIBMESH_DIM+j) + I_den(j,k) * _mo(i*LIBMESH_DIM+k);

      }
    }

    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      no(i*LIBMESH_DIM+j) = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
      {
        no(i*LIBMESH_DIM+j) = no(i*LIBMESH_DIM+j) +I_den(j,k) * _no(i*LIBMESH_DIM+k);

      }
    }
  }

  // Calculate Schmid tensor and resolved shear stresses
  for (unsigned int i = 0; i < _variable_size; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        flow_direction_s[i](j,k) = mo(i*LIBMESH_DIM+j) * no(i*LIBMESH_DIM+k);
}
bool
anaCrystalPlasticitySlipRateGSS::calcSlipRate(unsigned int qp, Real dt, std::vector<Real> & val) const
{
  DenseVector<Real> tau(_variable_size);
  DenseVector<Real> sigma(_variable_size);
  // Adding the eqautions 82
  DenseVector<Real> tau_new(_variable_size);
  DenseVector<Real> sigma_new(_variable_size);

  std::vector<Real> eigval(3, 0.0);
  std::vector<Real> eigval_new(3, 0.0);

  RankTwoTensor eigvec;
  RankTwoTensor eigvec_new;


  Real theta ;
  Real PI = 3.1415926;
//   Real mu = 0.04 ;
  Real phi = std::atan(_mu);
  theta =PI/4.0-phi/2.0;
  _pk2[qp].symmetricEigenvaluesEigenvectors(eigval_new, eigvec);
  eigval[0]=eigval_new[2];
  eigval[1]=eigval_new[1];
  eigval[2]=eigval_new[0];
  // using Cauchy
  //_stress[qp].symmetricEigenvaluesEigenvectors(eigval, eigvec);

  // std::cout<<"Cauchy="<<"\n"\
  // <<eigval[0]<<"|"<<eigval[1]<<"|"<<eigval[2]<<"\n"\
  // <<"Pk2="<<"\n"\
  // <<eigval_new[0]<<"|"<<eigval_new[1]<<"|"<<eigval_new[2]<<"\n";


// 1 2
  tau_new(0) = 0.5*std::sin(2*theta)*(eigval[0]-eigval[2]);
  sigma_new(0) = -0.5* (eigval[0]+eigval[2])+0.5*std::cos(2*theta)*(eigval[0]-eigval[2]);
  tau_new(1) = tau_new(0);
  sigma_new(1) = sigma_new(0);
// 3 4
  tau_new(2) = 0.5*std::sin(2*theta)*(eigval[0]-eigval[1]);
  sigma_new(2) = -0.5* (eigval[0]+eigval[1])+0.5*std::cos(2*theta)*(eigval[0]-eigval[1]);
  tau_new(3) = tau_new(2);
  sigma_new(3) = sigma_new(2);
// 5 6
  tau_new(4) = 0.5*std::sin(2*theta)*(eigval[1]-eigval[2]);
  sigma_new(4) = -0.5* (eigval[1]+eigval[2])+0.5*std::cos(2*theta)*(eigval[1]-eigval[2]);
  tau_new(5) = tau_new(4);
  sigma_new(5) = sigma_new(4);
  //


// Eq 82 Output

// std::cout<<"tau_new="<<"\n"\
// <<tau_new(0)<<"|"<<tau_new(1)<<"|"<<"\n"\
// <<tau_new(2)<<"|"<<tau_new(3)<<"|"<<"\n"\
// <<tau_new(4)<<"|"<<tau_new(5)<<"|"<<"\n"\
// <<"sigma_new="<<"\n"\
// <<sigma_new(0)<<"|"<<sigma_new(1)<<"|"<<"\n"\
// <<sigma_new(2)<<"|"<<sigma_new(3)<<"|"<<"\n"\
// <<sigma_new(4)<<"|"<<sigma_new(5)<<"|"<<"\n";


//
  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    // Anand
    Real _beta;
    _beta  = _g0*(std::pow((1.0-_mat_prop_state_var[qp][i]/_eta_cv),_p));
    Real _cohesion;
    Real _extra;
    _extra = _b*std::pow((1.0-_mat_prop_state_var[qp][i]/_eta_cv),_q);
     _cohesion = _c_cv+_extra;
         if ((std::abs(tau_new(i))-_beta*std::abs(sigma_new(i)))>0.0)
    {
      // Real _cohesion;
      //  _cohesion = _c_cv+_b*std::pow((1.0-_mat_prop_state_var[qp][i]/_eta_cv),_q);
     val[i] = _a0(i) * std::pow(std::abs(tau_new(i)) / (std::abs(_mu*sigma_new(i)+_cohesion)), 1.0 / _xm(i)) ;
    }
    else
    {
      val[i]=0.0;
    }
//  std::cout <<"no_c_cv"<<_b*std::pow((1.0-_mat_prop_state_var[qp][i]/_eta_cv),_q)<<"cohesion"<<_cohesion<<"mat_prop_state="<<_mat_prop_state_var[qp][i]<<"\n";

  //std::cout<<"qp"<<qp<<"i"<<i<<"cohesion"<<_cohesion<<"b"<<_b<<"_c_cv"<<_c_cv<<"eta"<<_mat_prop_state_var[qp][i]<<"eta_cv"<<_eta_cv<<"\n";

    //val[i] = _a0(i) * std::pow(std::abs(tau(i) / (_mu*sigma(i)+_mat_prop_state_var[qp][i])), 1.0 / _xm(i)) * copysign(1.0, tau(i));
    if (std::abs(val[i] * dt) > _slip_incr_tol)
    {
#ifdef DEBUG
      mooseWarning("Maximum allowable slip increment exceeded ",  std::abs(val[i])*dt);
#endif
      return false;
    }
  }

  return true;
}

bool
anaCrystalPlasticitySlipRateGSS::calcSlipRateDerivative(unsigned int qp, Real /*dt*/, std::vector<Real> & val) const
{
  DenseVector<Real> tau(_variable_size);
  DenseVector<Real> sigma(_variable_size);

  // Adding the eqautions 82
  DenseVector<Real> tau_new(_variable_size);
  DenseVector<Real> sigma_new(_variable_size);

  std::vector<Real> eigval(3, 0.0);
  std::vector<Real> eigval_new(3, 0.0);

  RankTwoTensor eigvec;

  Real theta ;
  Real PI = 3.1415926;
   //Real mu = 0.04 ;
  Real phi = std::atan(_mu);
  theta =PI/4.0- phi/2.0;
  _pk2[qp].symmetricEigenvaluesEigenvectors(eigval_new, eigvec);
  eigval[0]=eigval_new[2];
  eigval[1]=eigval_new[1];
  eigval[2]=eigval_new[0];
// Using cuachy stress

//_stress[qp].symmetricEigenvaluesEigenvectors(eigval, eigvec);

// 1 2
  tau_new(0) = 0.5*std::sin(2*theta)*(eigval[0]-eigval[2]);
  sigma_new(0) = -0.5* (eigval[0]+eigval[2])+0.5*std::cos(2*theta)*(eigval[0]-eigval[2]);
  tau_new(1) = tau_new(0);
  sigma_new(1) = sigma_new(0);
// 3 4
  tau_new(2) = 0.5*std::sin(2*theta)*(eigval[0]-eigval[1]);
  sigma_new(2) = -0.5* (eigval[0]+eigval[1])+0.5*std::cos(2*theta)*(eigval[0]-eigval[1]);
  tau_new(3) = tau_new(2);
  sigma_new(3) = sigma_new(2);
// 5 6
  tau_new(4) = 0.5*std::sin(2*theta)*(eigval[1]-eigval[2]);
  sigma_new(4) = -0.5* (eigval[1]+eigval[2])+0.5*std::cos(2*theta)*(eigval[1]-eigval[2]);
  tau_new(5) = tau_new(4);
  sigma_new(5) = sigma_new(4);
  //




  for (unsigned int i = 0; i < _variable_size; ++i)
  {
  //   Real _cohesion;
  //  _cohesion = _c_cv+_b*std::pow((1.0-_mat_prop_state_var[qp][i]/_eta_cv),_q);

   Real _cohesion;
   Real _extra;
   _extra = _b*std::pow((1.0-_mat_prop_state_var[qp][i]/_eta_cv),_q);
    _cohesion = _c_cv+_extra;
  //   if (_extra<=1.0)
  //   {
  //     _cohesion = _c_cv;
  //   }
  //  _cohesion =_c_cv;
//   val[i] = _a0(i)/_xm(i) * std::pow(tau(i) / (_mu*sigma(i)+_cohesion), 1.0 /_xm(i)-1.0)/(_mu*sigma(i)+_cohesion);
// Anand
   val[i] = _a0(i)/_xm(i) * std::pow(std::abs(tau_new(i)) / (std::abs(_mu*sigma_new(i)+_cohesion)), 1.0 /_xm(i)-1.0)/std::abs(_mu*sigma_new(i)+_cohesion);

//std::cout<<_mu*sigma_new(i)+_cohesion<<"\n";
  }

  //  val[i] = _a0(i) / _xm(i) * std::pow(std::abs(tau(i) / _mat_prop_state_var[qp][i]), 1.0 / _xm(i) - 1.0) / _mat_prop_state_var[qp][i];

  return true;
}
