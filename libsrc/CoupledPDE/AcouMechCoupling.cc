#include "AcouMechCoupling.hh"


#include "Driver/assemble.hh"
#include "Materials/baseMaterial.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/acousticPDE.hh" 

// integrator (bi-)linear forms
#include "Forms/acouMechInt.hh"


namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  AcouMechCoupling::AcouMechCoupling( SinglePDE *pde1, SinglePDE *pde2 )
    : BasePairCoupling( pde1, pde2 ) {

    ENTER_FCN( "AcouMechCoupling::AcouMechCoupling" );

    couplingName_ = "acouMechDirect";
    materialClass_ = FLUID;

    // Check, if matrix storage is SparseNonSym
    std::string mMatString;
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    keyVec  = "linearSystems", "system", "matrix", "storage";
    attrVec = "", "name", "";
    valVec  = "", "direct", "";
    params->Get( keyVec, attrVec, valVec, mMatString );
//     OLAS::MatrixStorageType mType;
//     OLAS::String2Enum( mMatString, mType );

    AcousticPDE* acouPDE  = dynamic_cast<AcousticPDE*> (pde2);
    SolutionType formulation = acouPDE->GetFormulation();

    if ( formulation ==  ACOU_PRESSURE && mMatString != "sparseNonSym" ) {
      Error("For MechAcou-Coupling with pressure formulation we need matrix storage: SPARSE_NONSYM",
	    __FILE__,__LINE__);
    }
  }


  // **************
  //   Destructor
  // **************
  AcouMechCoupling::~AcouMechCoupling() {
    ENTER_FCN( "AcouMechCoupling::~AcouMechCoupling" );
  }


  // ***************
  //   PostProcess
  // ***************
  void AcouMechCoupling::PostProcess() {
    ENTER_FCN( "AcouMechCoupling::PostProcess" );
  }


  // ***************
  //   Output of Results
  // ***************
  void AcouMechCoupling::WriteResultsInFile(const UInt kstep,
                                            const Double asteptime,
                                            UInt stepOffset,
                                            Double timeOffset){
    ENTER_FCN( "AcouMechCoupling::WriteResultsInFile" );
  }



  // *********************
  //   DefineIntegrators
  // *********************
  void AcouMechCoupling::DefineIntegrators() {

    ENTER_FCN( "AcouMechCoupling::DefineIntegrators" );
    
    //number of mechanical dofs = geometric dimension
    UInt dofs = ptGrid_->GetDim();
    
    bool isAxi = params->HasValue( "type", "axi", "geometry" );
    
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

       // create new entity list
      shared_ptr<SurfElemList> actSDList( new SurfElemList(ptGrid_ ) );
      actSDList->SetRegion( surfRegions_[actSD] );
      
       // =========  add first damp integrator (mechanic equation) ======
       SurfForm * massInt1 = new AcouMechInt(dofs,isAxi);
       massInt1->SetFormulation(formulation);
       
       // Set information about two neighbouring volume regions
       massInt1->SetFirstVoluInfo(pde1_->GetName(), 
                                  pde1_->getPDEMaterialData());
       massInt1->SetSecondVoluInfo(pde2_->GetName(), 
                                   pde2_->getPDEMaterialData());
       
       BiLinFormContext * massContext1 = 
         new BiLinFormContext(massInt1, STIFFNESS );
       massContext1->SetPtPdes(pde1_, pde2_);  
       massContext1->SetResults( results1_[0], results2_[0],
                                 actSDList, actSDList );    
       massContext1->SetCounterPart(false);
       
       // ========  add second mass integrator (acoustic equation)=======
       SurfForm * massInt2 = new AcouMechInt(dofs,isAxi);
       massInt2->SetFormulation(formulation);
       
       // Set information about two neighbouring volume regions
       massInt2->SetFirstVoluInfo(pde2_->GetName(), 
                                  pde2_->getPDEMaterialData());
       massInt2->SetSecondVoluInfo(pde1_->GetName(),
                                   pde1_->getPDEMaterialData());
       
       BiLinFormContext * massContext2 = 
         new BiLinFormContext(massInt2, MASS );
       massContext2->SetPtPdes(pde2_, pde1_);      
       massContext2->SetResults( results2_[0], results1_[0],
                                 actSDList, actSDList );    
       massContext2->SetCounterPart(false);

       assemble_->AddBiLinearForm( massContext1 );
       assemble_->AddBiLinearForm( massContext2 );
       
       
     }
    }
    else {      
      for ( UInt actSD = 0; actSD < surfRegions_.GetSize(); actSD++ ) {

        shared_ptr<SurfElemList> actSDList( new SurfElemList(ptGrid_ ) );
        actSDList->SetRegion( surfRegions_[actSD] );        

        // ========  add damping coupling ===================================
        SurfForm * dampInt = new AcouMechInt(dofs,isAxi);
        
        // Set information about two neighbouring volume regions
        dampInt->SetFirstVoluInfo(pde1_->GetName(), 
    			  pde1_->getPDEMaterialData());
        dampInt->SetSecondVoluInfo(pde2_->GetName(),
                                   pde2_->getPDEMaterialData());
        
        BiLinFormContext * dampContext = 
          new BiLinFormContext( dampInt, DAMPING );

        // We also need to set the transposed of the coupling
        // matrix to the lower diagonal side
        dampContext->SetCounterPart( true );

        dampContext->SetPtPdes(pde1_, pde2_);      
        dampContext->SetResults( results1_[0], results2_[0],
                                 actSDList, actSDList );
        assemble_->AddBiLinearForm( dampContext );
      }
   }
  }


  // ********************
  //   ReadStoreResults
  // ********************
  void AcouMechCoupling::ReadStoreResults() {
    ENTER_FCN( "AcouMechCoupling::ReadStoreResults" );
  }


  // ********************
  //   ReadMaterialData
  // ********************
  void AcouMechCoupling::ReadMaterialData() {
    ENTER_FCN( "AcouMechCoupling::ReadMaterialData" );
  }

}
