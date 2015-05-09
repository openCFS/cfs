#include "BasePDE.hh"
#include "Domain/Domain.hh"
#include "Driver/BaseDriver.hh"
#include "Driver/SolveSteps/BaseSolveStep.hh"


namespace CoupledField {


  Enum<BasePDE::AnalysisType> BasePDE::analysisType;

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
    analysisType.Add(HARMONIC, "paramIdent", false); // the value is not unique
    analysisType.Add(EIGENFREQUENCY, "eigenFrequency");
    analysisType.Add(INVERSESOURCE, "inverseSource");
    analysisType.Add(MULTI_SEQUENCE, "multiSequence");
  }
  
  std::map<FEMatrixType,Integer>  BasePDE::GetMatrixDerivativeMap(){
    std::map<FEMatrixType,Integer> retMap;
    retMap[MASS] = 2;
    retMap[DAMPING] = 1;
    retMap[STIFFNESS] = 0;
    retMap[AUXILIARY] = 0;
    retMap[CONVECTION] = 0;

    return retMap;
  }

}
