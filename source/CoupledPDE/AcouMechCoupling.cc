// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "AcouMechCoupling.hh"
#include "CoupledPDE/BasePairCoupling.hh"
#include "Domain/entityList.hh"
#include "Domain/grid.hh"
#include "Driver/assemble.hh"
#include "Driver/formsContext.hh"
// integrator (bi-)linear forms
#include "Forms/acouMechInt.hh"
#include "Forms/acouMechNcInt.hh"
#include "General/environment.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/acousticPDE.hh" 
#include "Utils/StdVector.hh"


namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  AcouMechCoupling::AcouMechCoupling( SinglePDE *pde1, SinglePDE *pde2,
                                      PtrParamNode paramNode  )
    : BasePairCoupling( pde1, pde2, paramNode,"acouMechDirect")
  {
    materialClass_ = FLUID;


    // NOTE: The information if the definition of the
    // bilinearforms is in conformance with the matrix layout should
    // be performed within the CFSOLASParams-class!

  //   // Check, if matrix storage is SparseNonSym
//     std::string mMatString;
//     StdVector<std::string> keyVec;
//     StdVector<std::string> attrVec;
//     StdVector<std::string> valVec;
//     keyVec  = "linearSystems", "system", "matrix", "storage";
//     attrVec = "", "name", "";
//     valVec  = "", "direct", "";
//     params->Get( keyVec, attrVec, valVec, mMatString );
// //     OLAS::MatrixStorageType mType;
// //     OLAS::String2Enum( mMatString, mType );

//     AcousticPDE* acouPDE  = dynamic_cast<AcousticPDE*> (pde2);
//     SolutionType formulation = acouPDE->GetFormulation();

//     if ( formulation ==  ACOU_PRESSURE && mMatString != "sparseNonSym" ) {
//       EXCEPTION( "For MechAcou-Coupling with pressure formulation we need matrix storage: SPARSE_NONSYM" );
//     }
  }


  // **************
  //   Destructor
  // **************
  AcouMechCoupling::~AcouMechCoupling() {
  }


  // ***************
  //   CalcResults
  // ***************
  void AcouMechCoupling::CalcResults() {
  }


  // ***************
  //   Output of Results
  // ***************
  void AcouMechCoupling::WriteResultsInFile(const UInt kstep,
                                            const Double asteptime,
                                            UInt stepOffset,
                                            Double timeOffset){
  }



  // *********************
  //   DefineIntegrators
  // *********************
  void AcouMechCoupling::DefineIntegrators() {

    
    //number of mechanical dofs = geometric dimension
    UInt dofs = ptGrid_->GetDim();
    
    // Abfrage ob pressure Formulierung, wenn ja
    // FEMatrixType destMat EFFECTIVE oder wie bisher DAMPING
    // Vorfaktoren? Faktor=1 oder Faktor f�r pressure
    // mit MASSINT setfaktor setzten
    // iterate over all subdomains
    
    // check for acoustic formulation
    AcousticPDE* acouPDE  = dynamic_cast<AcousticPDE*> (pde2_);
    SolutionType formulation = acouPDE->GetFormulation();
    
    
    if ( formulation ==  ACOU_PRESSURE) {
      for ( UInt actSD = 0; actSD < entityLists_.GetSize(); actSD++ ) {

        shared_ptr<EntityList> actSDList = entityLists_[actSD];
        
        // =========  add first damp integrator (mechanic equation) ======
       AcouMechInt * massInt1 = new AcouMechInt(dofs,isaxi_);
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
        AcouMechInt * massInt2 = new AcouMechInt(dofs,isaxi_);
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

      // Iterate over nonconforming interfaces
      for (UInt actNCI = 0; actNCI < ncIfaces_.GetSize(); ++actNCI) {
        // get surface elements of coupling region
        shared_ptr<SurfElemList> actNCList(new SurfElemList(ptGrid_));
        actNCList->SetRegion(ncIfaces_[actNCI]);

        //====== add mass integrator (acoustic eq.) ======
        AcouMechNcInt *ncMassInt = new AcouMechNcInt(dofs, isaxi_);
        ncMassInt->SetFormulation(ACOU_PRESSURE);

        // set material data
        ncMassInt->SetFirstVoluInfo(pde2_->GetName(),
                                    pde2_->getPDEMaterialData());
        ncMassInt->SetSecondVoluInfo(pde1_->GetName(),
                                     pde1_->getPDEMaterialData());

        // create context
        NcBiLinFormContext *massContext =
          new NcBiLinFormContext(ncMassInt, MASS);
        massContext->SetPtPdes(pde2_, pde1_);
        massContext->SetResults(results2_[0], results1_[0],
                                actNCList, actNCList);
        massContext->SetCounterPart(false);

        //====== add stiffness integrator (mechanic eq.) ======
        AcouMechNcInt *ncStiffInt = new AcouMechNcInt(dofs, isaxi_);
        ncStiffInt->SetFormulation(ACOU_PRESSURE);

        // set material data
        ncStiffInt->SetFirstVoluInfo(pde1_->GetName(),
                                     pde1_->getPDEMaterialData());
        ncStiffInt->SetSecondVoluInfo(pde2_->GetName(),
                                      pde2_->getPDEMaterialData());

        // create context
        NcBiLinFormContext *stiffContext =
          new NcBiLinFormContext(ncStiffInt, STIFFNESS);
        stiffContext->SetPtPdes(pde1_, pde2_);
        stiffContext->SetResults(results1_[0], results2_[0],
                                 actNCList, actNCList);
        stiffContext->SetCounterPart(false);

        //====== hand over to assemble object ======
        assemble_->AddBiLinearForm(massContext);
        assemble_->AddBiLinearForm(stiffContext);
      }
    }
    else {      
      for ( UInt actSD = 0; actSD < entityLists_.GetSize(); actSD++ ) {

        shared_ptr<EntityList> actSDList = entityLists_[actSD];

        // ========  add damping coupling ===================================
        AcouMechInt * dampInt = new AcouMechInt(dofs,isaxi_);
        
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

      // Iterate over nonconforming interfaces
      for (UInt actNCI = 0; actNCI < ncIfaces_.GetSize(); ++actNCI) {
        // get surface elements of coupling region
        shared_ptr<SurfElemList> actNCList(new SurfElemList(ptGrid_));
        actNCList->SetRegion(ncIfaces_[actNCI]);

        // define integrator
        AcouMechNcInt *ncInt = new AcouMechNcInt(dofs, isaxi_);
        ncInt->SetFormulation(ACOU_POTENTIAL);

        // set material data
        ncInt->SetFirstVoluInfo(pde1_->GetName(),
                                pde1_->getPDEMaterialData());
        ncInt->SetSecondVoluInfo(pde2_->GetName(),
                                 pde2_->getPDEMaterialData());

        // define damping context
        NcBiLinFormContext *dampContext =
          new NcBiLinFormContext(ncInt, DAMPING);

        dampContext->SetCounterPart(true); // coupling is symmetric
        dampContext->SetPtPdes(pde1_, pde2_);
        dampContext->SetResults(results1_[0], results2_[0],
                                actNCList, actNCList);

        // hand over to assemble object
        assemble_->AddBiLinearForm(dampContext);
      }
    }
  }


  // ********************
  //   ReadStoreResults
  // ********************
  void AcouMechCoupling::ReadStoreResults() {
  }


  // ********************
  //   ReadMaterialData
  // ********************
  void AcouMechCoupling::ReadMaterialData() {
  }

}
