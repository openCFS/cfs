#include "singleDriver.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "PDE/basePDE.hh"
#include "Domain/domain.hh"

namespace CoupledField{


  SingleDriver::SingleDriver(Domain * adomain, 
                             Integer stepOffset, 
                             Double timeOffset, 
                             std::string driverTag,
                             Boolean isPartOfSequence)
    : BaseDriver(adomain)
      
  {
  
    ENTER_FCN( "SingleDriver::SingleDriver" );

    stepOffset_ = stepOffset;
    timeOffset_ = timeOffset;
    driverTag_ = driverTag;
    isPartOfSequence_ = isPartOfSequence;
  }
  

  SingleDriver::~SingleDriver()
  {
    ENTER_FCN( "SingleDriver::~SingleDriver" );
  
  }

  void SingleDriver::SetPDE( BasePDE *pde) {
    ENTER_FCN( "SingleDriver::SetPDE" );
    ptPDE_ = pde;
  }

  void SingleDriver::GetMyPDEs()
  {
    ENTER_FCN( "SingleDriver::GetMyPDEs()" );

    StdVector<std::string> pdeNames;
    params->GetPDEList( pdeNames );
  
    // Initialize pdes with general Tag 'anyTag'
    StdVector<std::string> tags;
    tags.Resize(pdeNames.GetSize());
    tags.Init("anyTag");
    ptdomain_->InitPDEs(pdeNames,1,tags);
  
    ptPDE_ = ptdomain_->GetBasePDE();
  }
      
} // end of namespace
