#include "Domain/Domain.hh"
#include "BasePDE.hh"
#include "Driver/SolveSteps/BaseSolveStep.hh"


namespace CoupledField {


  Enum<BasePDE::AnalysisType> BasePDE::analysisType;

  // ***********************
  //   Default Constructor
  // ***********************
  BasePDE::BasePDE( PtrParamNode paramNode ) :
    sequenceStep_(0),
    myParam_(paramNode),
    pdename_()
  {
  }

  // **************
  //   Destructor
  // **************
  BasePDE::~BasePDE() {
  }

  bool BasePDE::IsComplex() 
  {
    return domain->GetDriver()->IsComplex();
  }
  
  void BasePDE::SetEnums()
  {
    analysisType.SetName("BasePDE::AnalysisType");
    analysisType.Add(STATIC, "static");
    analysisType.Add(TRANSIENT, "transient");
    analysisType.Add(HARMONIC, "harmonic");
    analysisType.Add(HARMONIC, "paramIdent", false); // the value is not unique
    analysisType.Add(EIGENFREQUENCY, "eigenFrequency");
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
