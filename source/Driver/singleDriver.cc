// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stddef.h>
#include <iostream>

#include "Domain/domain.hh"
#include "Driver/basedriver.hh"
#include "General/environment.hh"
#include "Utils/mathParser/mathParser.hh"
#include "singleDriver.hh"

namespace CoupledField{


  SingleDriver::SingleDriver( UInt sequenceStep,
                              bool isPartOfSequence )
    : BaseDriver()
      
  {
    sequenceStep_ = sequenceStep;
    isPartOfSequence_ = isPartOfSequence;
    ptPDE_ = NULL;
    
    // Set current value of time step and time step size in the mathParser
    mathParser_ = domain->GetMathParser(); 
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "t", 0.0 );
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "dt", 0.0 );    
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "step", 0 );        
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "f", 0.0 );

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

      std::cout << "++ Starting to solve problem" << std::endl;
    }
  }

  void SingleDriver::SetPDE( BasePDE *pde) {
    ptPDE_ = pde;
  }
      
} // end of namespace
