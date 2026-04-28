#include "SingleDriver.hh"

#include "PDE/BasePDE.hh"
#include "Domain/Domain.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/mathParser/mathParser.hh"

namespace CoupledField{


  SingleDriver::SingleDriver( UInt sequenceStep,
                              bool isPartOfSequence,
                              shared_ptr<SimState> state, Domain* domain,
                              PtrParamNode paramNode, PtrParamNode infoNode)
    : BaseDriver(state, domain, paramNode, infoNode)
      
  {
    sequenceStep_ = sequenceStep;
    isPartOfSequence_ = isPartOfSequence;
    ptPDE_ = NULL;
    writeAllSteps_ = false;
    approxSourceWithDeltaFnc_ = false;
    
    // Set current value of time step and time step size in the mathParser
    mathParser_ = domain_->GetMathParser();
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "t", 0.0 );
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "t0", 0.0 );
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "dt", 0.0 );
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "f", 0.0 );
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "step", 0 );

  }
  
  void SingleDriver::InitializePDEs() {
    // read in pde data
    if( ! isPartOfSequence_ ) {
      // Initialize pdes 
      domain_->CreatePDEs( 1, info_->GetParent() );
      ptPDE_ = domain_->GetBasePDE();

      ptPDE_->SetSourceApproxType( approxSourceWithDeltaFnc_);

      // for optimization we init the PDEs in Domain::PostInit(), because we need finer control over the steps.
      if(domain_->GetOptimization() == nullptr)
        domain_->InitPDEs( 1 );
    }
  }

  void SingleDriver::SetPDE( BasePDE *pde) {
    ptPDE_ = pde;
  }
      
} // end of namespace
