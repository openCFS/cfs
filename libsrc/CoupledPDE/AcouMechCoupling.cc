#include "AcouMechCoupling.hh"


#include "Driver/assemble.hh"
#include "DataInOut/MaterialData.hh"
#include "DataInOut/LoadMaterialData.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "PDE/SinglePDE.hh"


// integrator (bi-)linear forms
#include "Forms/acouMechInt.hh"


namespace CoupledField
{

  AcouMechCoupling::AcouMechCoupling( SinglePDE *pde1, SinglePDE *pde2 )
    : BasePairCoupling( pde1, pde2 ) {
    ENTER_FCN( "AcouMechCoupling::AcouMechCoupling" );

    couplingName_ = "acouMechDirect";
    materialClass_ = "fluid";

    // determine the needed matrices
    matrixTypes_.insert(DAMPING);
  }


  AcouMechCoupling::~AcouMechCoupling() {
    ENTER_FCN( "AcouMechCoupling::~AcouMechCoupling" )
      }

  void AcouMechCoupling::PostProcess() {
    ENTER_FCN( "AcouMechCoupling::PostProcess" );
  }
  

  void AcouMechCoupling::DefineIntegrators() {
    ENTER_FCN( "AcouMechCoupling::DefineIntegrators" );
    

    //number of mechanical dofs = geometric dimension
    UInt dofs = ptGrid_->GetDim();
    Boolean isAxi = params->HasValue( "type", "axi", "geometry" );
    

    // iterate over all subdomains
    for ( UInt actSD = 0; actSD < surfRegions_.GetSize(); actSD++ ) {
       
      

      // ==============  add mass ========================================
      SurfForm * massInt = new AcouMechInt(dofs,isAxi);

      // Set information about two neighbouring volume regions
      massInt->SetFirstVoluInfo(pde1_->GetName(), pde1_->getPDE_subdoms(), 
                                pde1_->getPDEMaterialData());
      massInt->SetSecondVoluInfo(pde2_->GetName(), pde2_->getPDE_subdoms(), 
                                 pde2_->getPDEMaterialData());

      assemble_->AddSurfIntegrator(massInt, surfRegions_[actSD], DAMPING, FALSE);

      //if pressure formulation: RHS integrator
      
    }

  }
  
  void AcouMechCoupling::ReadStoreResults() {
    ENTER_FCN( "AcouMechCoupling::ReadStoreResults" );
    
  }

  void AcouMechCoupling::ReadMaterialData() {
    ENTER_FCN( "AcouMechCoupling::ReadMaterialData" );
    
  }
  


} // end of namespace
