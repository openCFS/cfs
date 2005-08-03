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
    
    // Abfrage ob pressure Formulierung, wenn ja
    // FEMatrixType destMat EFFECTIVE oder wie bisher DAMPING
    // Vorfaktoren? Faktor=1 oder Faktor für pressure
    // mit MASSINT setfaktor setzten
    // iterate over all subdomains

    // check for acoustic formulation

    std::string str;
    SolutionType formulation;
    params->Get( "formulation", str, "pdeList", pde2_->GetName() );
    String2Enum( str, formulation );

    if ( formulation ==  ACOU_PRESSURE) {
      for ( UInt actSD = 0; actSD < surfRegions_.GetSize(); actSD++ ) {

       	// ==============  add first mass integrator (mechanic equation)=========
	SurfForm * massInt1 = new AcouMechInt(dofs,isAxi);
	massInt1->SetFormulation(formulation);

	// Set information about two neighbouring volume regions
	massInt1->SetFirstVoluInfo(pde1_->GetName(), pde1_->getPDE_subdoms(), 
				  pde1_->getPDEMaterialData());
	massInt1->SetSecondVoluInfo(pde2_->GetName(), pde2_->getPDE_subdoms(), 
				   pde2_->getPDEMaterialData());
	
	IntegratorDescriptor * massDescr1 = 
	  new IntegratorDescriptor(massInt1, DAMPING, FALSE);
	massDescr1->SetPDEIds(pde1_, pde2_);      
	massDescr1->SetCounterPart(false);
	
	assemble_->AddSurfIntegrator(massDescr1, surfRegions_[actSD]);

       	// ==============  add second mass integrator (acoustic equation)=========
	SurfForm * massInt2 = new AcouMechInt(dofs,isAxi);
	massInt2->SetFormulation(formulation);

	// Set information about two neighbouring volume regions
	massInt2->SetFirstVoluInfo(pde2_->GetName(), pde2_->getPDE_subdoms(), 
				  pde2_->getPDEMaterialData());
	massInt2->SetSecondVoluInfo(pde1_->GetName(), pde1_->getPDE_subdoms(), 
				   pde1_->getPDEMaterialData());
	
	IntegratorDescriptor * massDescr2 = 
	  new IntegratorDescriptor(massInt2, DAMPING, FALSE);
	massDescr2->SetPDEIds(pde2_, pde1_);      
	massDescr2->SetCounterPart(false);
	
	assemble_->AddSurfIntegrator(massDescr2, surfRegions_[actSD]);
      }
    }
    else {      
      for ( UInt actSD = 0; actSD < surfRegions_.GetSize(); actSD++ ) {

	// ==============  add mass ========================================
	SurfForm * massInt = new AcouMechInt(dofs,isAxi);
	
	// Set information about two neighbouring volume regions
	massInt->SetFirstVoluInfo(pde1_->GetName(), pde1_->getPDE_subdoms(), 
				  pde1_->getPDEMaterialData());
	massInt->SetSecondVoluInfo(pde2_->GetName(), pde2_->getPDE_subdoms(), 
				   pde2_->getPDEMaterialData());
	
	IntegratorDescriptor * massDescr = 
	  new IntegratorDescriptor(massInt, DAMPING, FALSE);
	massDescr->SetPDEIds(pde1_, pde2_);      
	
	assemble_->AddSurfIntegrator(massDescr, surfRegions_[actSD]);
      }
    }

  }
  
  void AcouMechCoupling::ReadStoreResults() {
    ENTER_FCN( "AcouMechCoupling::ReadStoreResults" );
    
  }

  void AcouMechCoupling::ReadMaterialData() {
    ENTER_FCN( "AcouMechCoupling::ReadMaterialData" );
    
  }
  


} // end of namespace
