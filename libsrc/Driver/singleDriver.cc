#include "singleDriver.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"


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
  timeOffset_ = timeOffset_;
  driverTag_ = driverTag;
  isPartOfSequence_ = isPartOfSequence;
}
  

SingleDriver::~SingleDriver()
{
  ENTER_FCN( "SingleDriver::~SingleDriver" );
  
}

void SingleDriver::GetMyPDEs()
{
  ENTER_FCN( "SingleDriver::GetMyPDEs()" );

  StdVector<std::string> pdeNames;
#ifndef XMLPARAMS
      conf->getliststr("list_pdes",pdeNames);
#else
      params->GetPDEList( pdeNames );
#endif

      // Initialize pdes with general Tag 'anyTag'
      pdes_.Resize(pdeNames.GetSize());
      StdVector<std::string> tags;
      tags.Resize(pdes_.GetSize());
      tags.Init("anyTag");
      ptdomain_->InitPDEs(pdeNames,1,tags);
      
      // Get pointers to PDEs back
      for (Integer iPDE=0; iPDE<pdeNames.GetSize(); iPDE++)
	pdes_[iPDE] = ptdomain_->GetPDE(pdeNames[iPDE]);
    
}
      
} // end of namespace
