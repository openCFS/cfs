#include "singleDriver.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "PDE/basePDE.hh"
#include "Domain/domain.hh"

namespace CoupledField{


  SingleDriver::SingleDriver( std::string driverTag,
                              bool isPartOfSequence )
    : BaseDriver()
      
  {
  
    ENTER_FCN( "SingleDriver::SingleDriver" );

    driverTag_ = driverTag;
    isPartOfSequence_ = isPartOfSequence;
    ptPDE_ = NULL;
  }
  

  SingleDriver::~SingleDriver()
  {
    ENTER_FCN( "SingleDriver::~SingleDriver" );
  
  }

  void SingleDriver::InitializePDEs() {
    ENTER_FCN( "void InitializePDEs()" );
   
       // read in pde data 
    if( ! isPartOfSequence_ ) {
      StdVector<std::string> pdeNames;
      params->GetPDEList( pdeNames );
      
      // Initialize pdes with general Tag 'anyTag'
      StdVector<std::string> tags;
      tags.Resize(pdeNames.GetSize());
      tags.Init( driverTag_ );
      domain->CreatePDEs( pdeNames, 1, tags );
      ptPDE_ = domain->GetBasePDE();

      // Trigger reading of restart file
      ReadRestart();

      domain->InitPDEs(1, tags );
     
      Info->StartProgress ("Starting to solve problem", false);
    }
  }

  void SingleDriver::SetPDE( BasePDE *pde) {
    ENTER_FCN( "SingleDriver::SetPDE" );
    ptPDE_ = pde;
  }

//   void SingleDriver::GetMyPDEs()
//   {
//     ENTER_FCN( "SingleDriver::GetMyPDEs()" );

//     StdVector<std::string> pdeNames;
//     params->GetPDEList( pdeNames );
  
//     // Initialize pdes with general Tag 'anyTag'
//     StdVector<std::string> tags;
//     tags.Resize(pdeNames.GetSize());
//     tags.Init("anyTag");
//     domain->CreatePDEs( pdeNames, 1, tags );
//     ptPDE_ = domain->GetBasePDE();
//   }
      
} // end of namespace
