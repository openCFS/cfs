#include "SingleDriver.hh"

#include "PDE/BasePDE.hh"
#include "Domain/Domain.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField{


  SingleDriver::SingleDriver( UInt sequenceStep,
                              bool isPartOfSequence,
                              shared_ptr<SimState> state, Domain* domain  )
    : BaseDriver(state, domain)
      
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
      domain_->CreatePDEs( 1 );
      ptPDE_ = domain_->GetBasePDE();

      domain_->InitPDEs( 1 );
      // Trigger reading of restart file
      ReadRestart();

      std::cout << "++ Starting to solve problem" << std::endl;
    }
  }

  void SingleDriver::SetPDE( BasePDE *pde) {
    ptPDE_ = pde;
  }
      
} // end of namespace
