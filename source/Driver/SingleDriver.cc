#include "SingleDriver.hh"

#include "PDE/BasePDE.hh"
#include "Domain/Domain.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

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
  }
  
  SingleDriver::~SingleDriver()
  {
  
  }

  void SingleDriver::InitializePDEs() {
   
       // read in pde data 
    if( ! isPartOfSequence_ ) {
      
      // Initialize pdes 
      domain_->CreatePDEs( 1, info_->GetParent() );
      ptPDE_ = domain_->GetBasePDE();

      domain_->InitPDEs( 1 );

      std::cout << "++ Starting to solve problem" << std::endl;
    }
  }

  void SingleDriver::SetPDE( BasePDE *pde) {
    ptPDE_ = pde;
  }
      
} // end of namespace
