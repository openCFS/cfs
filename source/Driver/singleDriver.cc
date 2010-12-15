// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "singleDriver.hh"

#include "PDE/basePDE.hh"
#include "Domain/domain.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField{


  SingleDriver::SingleDriver( UInt sequenceStep,
                              bool isPartOfSequence )
    : BaseDriver()
      
  {
  

    sequenceStep_ = sequenceStep;
    isPartOfSequence_ = isPartOfSequence;
    ptPDE_ = NULL;
  }
  
  SingleDriver::~SingleDriver()
  {
  
  }

  void SingleDriver::InitializePDEs() {
   
       // read in pde data 
    if( ! isPartOfSequence_ ) {
      
      // Initialize pdes 
      domain->CreatePDEs( 1 );
      ptPDE_ = domain->GetBasePDE();

      domain->InitPDEs( 1 );
      // Trigger reading of restart file
      ReadRestart();

     
      Info->StartProgress ("Starting to solve problem", false);
    }
  }

  void SingleDriver::SetPDE( BasePDE *pde) {
    ptPDE_ = pde;
  }
      
} // end of namespace
