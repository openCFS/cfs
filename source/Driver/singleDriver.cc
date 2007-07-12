// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "singleDriver.hh"

#include "PDE/basePDE.hh"
#include "Domain/domain.hh"
#include "DataInOut/WriteInfo.hh"

namespace CoupledField{


  SingleDriver::SingleDriver( UInt sequenceStep,
                              bool isPartOfSequence )
    : BaseDriver()
      
  {
  
    ENTER_FCN( "SingleDriver::SingleDriver" );

    sequenceStep_ = sequenceStep;
    isPartOfSequence_ = isPartOfSequence;
    ptPDE_ = NULL;

  }
  
  SingleDriver::~SingleDriver()
  {
    ENTER_FCN( "SingleDriver::~SingleDriver" );
  
  }

  void SingleDriver::SolveProblem()
  {
    BaseDriver::SolveProblem();
  }
   

  void SingleDriver::InitializePDEs() {
    ENTER_FCN( "void InitializePDEs()" );
   
       // read in pde data 
    if( ! isPartOfSequence_ ) {
      
      // Initialize pdes 
      domain->CreatePDEs( 1 );
      ptPDE_ = domain->GetBasePDE();

      // Trigger reading of restart file
      ReadRestart();

      domain->InitPDEs( 1 );
     
      Info->StartProgress ("Starting to solve problem", false);
    }
  }

  void SingleDriver::SetPDE( BasePDE *pde) {
    ENTER_FCN( "SingleDriver::SetPDE" );
    ptPDE_ = pde;
  }
      
} // end of namespace
