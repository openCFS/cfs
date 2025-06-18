#include "BasePDE.hh"
#include "Domain/Domain.hh"
#include "Driver/BaseDriver.hh"
#include "Driver/SolveSteps/BaseSolveStep.hh"
#include "PDE/MechPDE.hh"


namespace CoupledField {


  Enum<BasePDE::AnalysisType> BasePDE::analysisType;

  Enum<BasePDE::StabilisationType> BasePDE::stabilisationType;

  // ***********************
  //   Default Constructor
  // ***********************
  BasePDE::BasePDE( PtrParamNode paramNode,
                    PtrParamNode infoNode,
                    shared_ptr<SimState> simState,
                    Domain * ptDomain) :
    sequenceStep_(0),
    myParam_(paramNode),
    myInfo_(infoNode),
    pdename_(),
    simState_(simState),
    domain_(ptDomain),
    mp_(ptDomain->GetMathParser())
  {
	  approxSourceWithDeltaFnc_ = false;
  }

  // **************
  //   Destructor
  // **************
  BasePDE::~BasePDE() {
  }

  bool BasePDE::IsComplex() 
  {
    return domain_->GetDriver()->IsComplex();
  }
  
  void BasePDE::SetEnums()
  {
    analysisType.SetName("BasePDE::AnalysisType");
    analysisType.Add(NO_ANALYSIS, "undefined");
    analysisType.Add(STATIC, "static");
    analysisType.Add(TRANSIENT, "transient");
    analysisType.Add(HARMONIC, "harmonic");
    analysisType.Add(MULTIHARMONIC, "multiharmonic");
    analysisType.Add(HARMONIC, "paramIdent", false); // the value is not unique
    analysisType.Add(EIGENFREQUENCY, "eigenFrequency");
    analysisType.Add(INVERSESOURCE, "inverseSource");
    analysisType.Add(MULTI_SEQUENCE, "multiSequence");
    analysisType.Add(BUCKLING, "buckling");
    analysisType.Add(EIGENVALUE, "eigenValue");
    analysisType.Add(HARMONIC25D, "harmonic25d");


    MechPDE::testStrain.SetName("MechPDE::TestStrain");
    MechPDE::testStrain.Add(MechPDE::X, "x");
    MechPDE::testStrain.Add(MechPDE::Y, "y");
    MechPDE::testStrain.Add(MechPDE::Z, "z");
    MechPDE::testStrain.Add(MechPDE::YZ, "yz");
    MechPDE::testStrain.Add(MechPDE::XZ, "xz");
    MechPDE::testStrain.Add(MechPDE::XY, "xy");

    stabilisationType.SetName("BasePDE::StabilisationType");
    stabilisationType.Add(NO_STABILISATION, "none");
    stabilisationType.Add(SUPG, "SUPG");
    stabilisationType.Add(ARTIFICIAL_DIFFUSION, "ArtificialDiffusion");

  }
  
  std::map<FEMatrixType,Integer>  BasePDE::GetMatrixDerivativeMap(){
    std::map<FEMatrixType,Integer> retMap;
    retMap[MASS] = 2;
    retMap[MASS_UPDATE] = 2;
    retMap[DAMPING] = 1;
    retMap[DAMPING_UPDATE] = 1;
    retMap[STIFFNESS] = 0;
    retMap[STIFFNESS_UPDATE] = 0;
    retMap[AUXILIARY] = 0;
    retMap[CONVECTION] = 0;
    retMap[GEOMETRIC_STIFFNESS] = 0;


    return retMap;
  }

}
