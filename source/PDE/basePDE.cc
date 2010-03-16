// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "basePDE.hh"
#include "Driver/baseSolveStep.hh"


namespace CoupledField {


  Enum<BasePDE::AnalysisType> BasePDE::analysisType;

  // ***********************
  //   Default Constructor
  // ***********************
  BasePDE::BasePDE( PtrParamNode paramNode ) :
    converged_(false),
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

  bool BasePDE::IsComplex(AnalysisType type) 
  {
    switch(type) 
    {
    case STATIC:
      return false;

    case TRANSIENT:
      return false;
      
    case HARMONIC:
      return true;
      
    case EIGENFREQUENCY:
      return false;
      break;

    default:
      EXCEPTION("type not implemented")
    };
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
  
  
  
  
}
