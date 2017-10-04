#include "StzApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

//kernels
#include "stzStressDivergenceTensors.h"
//
#include "stzHeatConduction.h"
#include "stzHeatConductionTimeDerivative.h"
//#include "stzHeatSource.h"
#include "stzMaskedBodyForce.h"
#include "stzRMaskedBodyForce.h"
#include "stzagrainMaskedBodyForce.h"
#include "stzchiMaskedBodyForce.h"
#include "stzpMaskedBodyForce.h"

// STZ COUPLE kernels
#include "stzHeatConductionCouple.h"
#include "stzHeatConductionTimeDerivativeCouple.h"
#include "stzMaskedBodyForceCouple.h"


// STZ 1D
#include "stzMaskedBodyForce1D.h"



// Contacts
//functions
#include "stzExampleFunction.h"
#include "stzshearExampleFunction.h"
#include "stzshearzeroFunction.h"
#include "stzvelstepFunction.h"



//material
#include "stzFiniteStrainHyperElasticViscoPlastic.h"
#include "anaFiniteStrainHyperElasticViscoPlastic.h"
#include "anaFiniteStrainUObasedCP.h"
#include "anaComputeElasticityTensorCP.h"
// STZ +anand Model
#include "stzanaFiniteStrainUObasedCP.h"
// STZ coupled
#include "stzanaFiniteStrainUObasedCPCouple.h"


//#include "stzExampleMaterial.h"
#include "testExampleMaterial.h"
#include "testnewExampleMaterial.h"

#include "RtestExampleMaterial.h"

// 1D STZ material
#include "testExampleMaterial1D.h"




//userobjects
#include "stzHEVPFlowRatePowerLawJ2.h"
#include "stznewHEVPFlowRatePowerLawJ2.h"

#include "stzcHEVPFlowRatePowerLawJ2.h"
#include "stzRHEVPFlowRatePowerLawJ2.h"
#include "stzVOLFlowRate.h"
#include "stzHEVPLinearHardening.h"


// Anand Model

#include "anaHEVPEqvPlasticStrain.h"
#include "anaHEVPEqvPlasticStrainRate.h"
#include "anaHEVPFlowRatePowerLawJ2.h"
#include "anaHEVPLinearHardening.h"
#include "anaCrystalPlasticitySlipRateGSS.h"
#include "anaCrystalPlasticitySlipResistanceGSS.h"
#include "anaCrystalPlasticityStateVariable.h"
#include "anaCrystalPlasticityStateVarRateComponentGSS.h"

// STZ+Anand model
#include "stzanaCrystalPlasticitySlipRateGSS.h"
#include "stzanaCrystalPlasticitySlipResistanceGSS.h"
#include "stzanaCrystalPlasticityStateVariable.h"
#include "stzanaCrystalPlasticityStateVarRateComponentGSS.h"

//
#include "stznewHEVPLinearHardening.h"

#include "stzRHEVPLinearHardening.h"
#include "stzchiHEVPLinearHardening.h"

//STZ coupling model
#include "stzanaCrystalPlasticitySlipRateGSSCouple.h"
#include "stzanaCrystalPlasticitySlipResistanceGSSCouple.h"
#include "stzanaCrystalPlasticityStateVariableCouple.h"
#include "stzanaCrystalPlasticityStateVarRateComponentGSSCouple.h"

// STZ 1D User objects

#include "stzHEVPFlowRatePowerLawJ21D.h"
#include "stzHEVPLinearHardening1D.h"


// STZ 1D its own kernel
#include "stzsdot.h"
#include "stzsheatsource.h"
#include "stzchisource.h"
#include "stzplasticrateAux.h"
#include "stzchidiffusion.h"
#include "stzImplicitODEs.h"
#include "stzplasticrateAux_scalars.h"
#include "stzchisource_scalars.h"
#include "stzplasticrateAux_s0.h"
#include "stzplasticrateAux_m.h"
#include "stzplasticrateAux_T.h"
#include "stzpfsource_chi.h"
#include "stzpfdiffusion.h"
#include "stzTsource_chi.h"
#include "stzTdiffusion.h"
#include "stzpfsource_T.h"
#include "stzplasticrateAux_chidot.h"
#include "stzplasticrateAux_stzdensity.h"
#include "stzplasticrateAux_bwdot.h"
#include "stzPostprocessorDT.h"
#include "stzCoupledScalarAux.h"
#include "stzsheatsource_singleblock.h"
#include "stzchisource_singleblock.h"






template<>
InputParameters validParams<StzApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  params.set<bool>("use_legacy_output_syntax") = false;

  return params;
}

StzApp::StzApp(InputParameters parameters) :
    MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  StzApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  StzApp::associateSyntax(_syntax, _action_factory);
}

StzApp::~StzApp()
{
}

// External entry point for dynamic application loading
extern "C" void StzApp__registerApps() { StzApp::registerApps(); }
void
StzApp::registerApps()
{
  registerApp(StzApp);
}

// External entry point for dynamic object registration
extern "C" void StzApp__registerObjects(Factory & factory) { StzApp::registerObjects(factory); }
void
StzApp::registerObjects(Factory & factory)
{
    //Kernels
    registerKernel(stzStressDivergenceTensors);
//
  registerKernel(stzHeatConduction);
  registerKernel(stzHeatConductionTimeDerivative);
//  registerKernel(stzHeatSource);
  registerKernel(stzMaskedBodyForce);
  registerKernel(stzRMaskedBodyForce);
  registerKernel(stzagrainMaskedBodyForce);
  registerKernel(stzchiMaskedBodyForce);
  registerKernel(stzpMaskedBodyForce);
// Kernel STZ COUPLEs
registerKernel(stzHeatConductionCouple);
registerKernel(stzHeatConductionTimeDerivativeCouple);
registerKernel(stzMaskedBodyForceCouple);

// STZ 1D
registerKernel(stzMaskedBodyForce1D);


  // Contact

    //Material
  registerMaterial(stzFiniteStrainHyperElasticViscoPlastic);
  registerMaterial(anaFiniteStrainHyperElasticViscoPlastic);
  registerMaterial(anaFiniteStrainUObasedCP);
  registerMaterial(anaComputeElasticityTensorCP);

  // STZ + Anand Model
  registerMaterial(stzanaFiniteStrainUObasedCP);
  // STZ Coupled Material
  registerMaterial(stzanaFiniteStrainUObasedCPCouple);



  //registerKernel(stzExampleMaterial);
  registerMaterial(testExampleMaterial);
  registerMaterial(testnewExampleMaterial);

  registerMaterial(RtestExampleMaterial);

  // 1D STZ material
  registerMaterial(testExampleMaterial1D);


  // functions
  registerFunction(stzExampleFunction);
  registerFunction(stzshearExampleFunction);
  registerFunction(stzshearzeroFunction);
  registerFunction(stzvelstepFunction);


    //Userobjects
  registerUserObject(stzHEVPFlowRatePowerLawJ2);
  registerUserObject(stznewHEVPFlowRatePowerLawJ2);

  registerUserObject(stzcHEVPFlowRatePowerLawJ2);
  registerUserObject(stzRHEVPFlowRatePowerLawJ2);
  registerUserObject(stzVOLFlowRate);

  registerUserObject(stzHEVPLinearHardening);

  // Anand Model
  registerUserObject(anaHEVPEqvPlasticStrain);
  registerUserObject(anaHEVPEqvPlasticStrainRate);
  registerUserObject(anaHEVPFlowRatePowerLawJ2);
  registerUserObject(anaHEVPLinearHardening);
  registerUserObject(anaCrystalPlasticitySlipRateGSS);
  registerUserObject(anaCrystalPlasticitySlipResistanceGSS);
  registerUserObject(anaCrystalPlasticityStateVariable);
  registerUserObject(anaCrystalPlasticityStateVarRateComponentGSS);

  // STZ+Anand Model
  registerUserObject(stzanaCrystalPlasticitySlipRateGSS);
  registerUserObject(stzanaCrystalPlasticitySlipResistanceGSS);
  registerUserObject(stzanaCrystalPlasticityStateVariable);
  registerUserObject(stzanaCrystalPlasticityStateVarRateComponentGSS);




//
  registerUserObject(stznewHEVPLinearHardening);

  registerUserObject(stzRHEVPLinearHardening);
  registerUserObject(stzchiHEVPLinearHardening);

  // STZ coupling formulation
  registerUserObject(stzanaCrystalPlasticitySlipRateGSSCouple);
  registerUserObject(stzanaCrystalPlasticitySlipResistanceGSSCouple);
  registerUserObject(stzanaCrystalPlasticityStateVariableCouple);
  registerUserObject(stzanaCrystalPlasticityStateVarRateComponentGSSCouple);

  // STZ 1D user objects
  registerUserObject(stzHEVPFlowRatePowerLawJ21D);
  registerUserObject(stzHEVPLinearHardening1D);


  // STZ 1D its own kernels
  registerKernel(stzsdot);
  registerKernel(stzsheatsource);
  registerKernel(stzchisource);
  registerAux(stzplasticrateAux);
  registerKernel(stzchidiffusion);
  registerScalarKernel(stzImplicitODEs);
  registerAux(stzplasticrateAux_scalars);
  registerKernel(stzchisource_scalars);
  registerAux(stzplasticrateAux_s0);
  registerAux(stzplasticrateAux_m);
  registerAux(stzplasticrateAux_T);
  registerKernel(stzpfsource_chi);
  registerKernel(stzpfdiffusion);
  registerKernel(stzTsource_chi);
  registerKernel(stzTdiffusion);
  registerKernel(stzpfsource_T);
  registerAux(stzplasticrateAux_chidot);
  registerAux(stzplasticrateAux_stzdensity);
  registerAux(stzplasticrateAux_bwdot);
  registerTimeStepper(stzPostprocessorDT);
  registerAux(stzCoupledScalarAux);
  registerKernel(stzsheatsource_singleblock);
  registerKernel(stzchisource_singleblock);






}

// External entry point for dynamic syntax association
extern "C" void StzApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { StzApp::associateSyntax(syntax, action_factory); }
void
StzApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
